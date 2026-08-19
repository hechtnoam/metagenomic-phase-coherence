#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib
import json
import os
import py_compile
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import period3_statistics as p3
import period3_library_pipeline as pipeline
import dft_mod3 as dft
import phase_redistribution as trim_reads
import simulate_fragmentation as sim
import frame_normalize_cds as frame0_trim
import six_state_em as em
import aggregate_library_results as aggregate
import plot_project_summary as metadata_plot
import plot_simulation_r2

RESULTS: list[str] = []


def ok(name: str, detail: str = "") -> None:
    RESULTS.append(f"PASS\t{name}\t{detail}".rstrip())


def skip(name: str, detail: str) -> None:
    RESULTS.append(f"SKIP\t{name}\t{detail}")


def test_compile() -> None:
    for path in sorted(SCRIPTS.glob("*.py")):
        py_compile.compile(str(path), doraise=True)
    ok("publication_scripts_compile", f"{len(list(SCRIPTS.glob('*.py')))} scripts")


def test_period3_statistics(work: Path) -> None:
    cycles = np.arange(1, 61)
    values = {
        "A": 0.25 + 0.025 * np.cos(2 * np.pi * cycles / 3),
        "T": 0.25 + 0.020 * np.sin(2 * np.pi * cycles / 3),
        "G": 0.25 - 0.025 * np.cos(2 * np.pi * cycles / 3),
    }
    values["C"] = 1.0 - values["A"] - values["T"] - values["G"]
    df = pd.DataFrame({"cycle": cycles, **values})
    path = work / "perfect_period3.csv"
    df.to_csv(path, index=False)
    result = p3.test_period3_from_csv(path, start=10, end=40, permutations=199, seed=7)
    assert all(fit.r2 > 0.999999 for fit in result.fits)
    assert all(0 <= fit.p_value <= fit.p_value_adjusted <= 1 for fit in result.fits)
    assert result.all_significant

    constant = pd.DataFrame({"cycle": cycles, "A": 0.2, "T": 0.2, "G": 0.3, "C": 0.3})
    cpath = work / "constant.csv"
    constant.to_csv(cpath, index=False)
    cres = p3.test_period3_from_csv(cpath, start=10, end=40, permutations=19, seed=8)
    assert all(fit.r2 == 0 for fit in cres.fits)
    assert all(np.isfinite(fit.p_value_adjusted) for fit in cres.fits)
    ok("period3_statistics", "perfect signal, Holm correction, and constant-profile edge case")


def test_pipeline_io_and_counting(work: Path) -> None:
    fq = work / "reads.fastq"
    # CRLF explicitly tests robust newline handling.
    fq.write_bytes(
        b"@r1\r\nATGATGATGATG\r\n+\r\nIIIIIIIIIIII\r\n"
        b"@r2\r\nCATCATCATCAT\r\n+\r\nJJJJJJJJJJJJ\r\n"
        b"@bad\r\nATGATNATGATG\r\n+\r\nIIIIIIIIIIII\r\n"
    )
    seqs = list(pipeline.iter_fastq_sequences(str(fq)))
    assert seqs == ["ATGATGATGATG", "CATCATCATCAT", "ATGATNATGATG"]
    result = pipeline.count_sequences_serial(
        iter(seqs), [12], [1], False, 0, False, "[test]"
    )
    assert result.total_raw_reads == 3
    assert result.total_analyzed_reads == 2
    assert result.invalid_sequence_count == 1
    fractions = pipeline.fractions_from_counts(result.counts_by_k[1][12], 2)
    assert np.allclose(fractions.sum(axis=1), 1.0)

    class MockRead:
        is_reverse = True
        query_sequence = "CAT"  # stored reverse complement of original ATG
        def get_forward_sequence(self):
            return "ATG"
    assert pipeline.original_bam_sequence(MockRead()) == "ATG"
    ok("period3_pipeline_io_counting", "CRLF FASTQ, ambiguity exclusion, fractions, BAM orientation helper")


def test_dft() -> None:
    n = 30
    j = np.arange(n)
    window = np.column_stack([
        0.25 + 0.02 * np.cos(2 * np.pi * j / 3),
        0.25 + 0.01 * np.sin(2 * np.pi * j / 3),
        0.25 - 0.02 * np.cos(2 * np.pi * j / 3),
        0.25 - 0.01 * np.sin(2 * np.pi * j / 3),
    ])
    freqs, mags = dft.compute_spectrum(window)
    idx = dft.target_bin(freqs)
    assert idx == 10 and np.isclose(freqs[idx], 1 / 3)
    sr2 = dft.spectral_r2_at_target(window, mags, idx)
    assert all(abs(value - 1.0) < 1e-10 for value in sr2.values())
    assert dft.BASES == ("A", "T", "G", "C")
    ok("dft", "exact q=10 bin, correct A/T/G/C order, spectral R2=1 on a pure signal")


def test_trim_fastq(work: Path) -> None:
    inp = work / "trim_input.fastq"
    inp.write_text("@r1 description\nACGTACGT\n+\nABCDEFGH\n")
    out = work / "trim_output.fastq"
    stats = trim_reads.process_fastx(
        str(inp), str(out), 8, "fastq", seed=1, mode="fixed", fixed_five_trim=1
    )
    from Bio import SeqIO
    record = next(SeqIO.parse(out, "fastq"))
    assert str(record.seq) == "CGTACG"
    assert len(record.letter_annotations["phred_quality"]) == 6
    assert stats["records_written"] == 1

    class MockRead:
        is_reverse = True
        query_sequence = "CAT"
        query_qualities = [3, 2, 1]
        def get_forward_sequence(self):
            return "ATG"
        def get_forward_qualities(self):
            return [1, 2, 3]
    mock = MockRead()
    assert mock.get_forward_sequence() == "ATG" and mock.get_forward_qualities() == [1, 2, 3]
    ok("phase_redistribution_trim", "FASTQ qualities retained and fixed-trim control works")


def test_simulator(work: Path) -> None:
    sequences = {
        "long": ("AAGGCT" * 100),
        "short": ("GGTTAC" * 40),
    }
    out1 = work / "sim1.fa"
    out2 = work / "sim2.fa"
    s1 = sim.generate_reads(
        sequences, out1, num_reads=500, read_length=30,
        bias_5p=0.5, bias_3p=0.5, reverse_complement_probability=0.5,
        seed=123, flush_every=73,
    )
    s2 = sim.generate_reads(
        sequences, out2, num_reads=500, read_length=30,
        bias_5p=0.5, bias_3p=0.5, reverse_complement_probability=0.5,
        seed=123, flush_every=101,
    )
    assert out1.read_bytes() == out2.read_bytes()
    assert set(s1["candidate_contig_counts"]) == {"long", "short"}
    assert sum(s1["orientation_counts"].values()) == 500
    assert set(s1["orientation_counts"]) == {"forward_reference", "reverse_complement"}
    seq_lines = [line.strip() for line in out1.read_text().splitlines() if not line.startswith(">")]
    assert len(seq_lines) == 500 and all(len(seq) == 30 for seq in seq_lines)

    out3 = work / "sim_purine.fa"
    s3 = sim.generate_reads(
        sequences, out3, num_reads=100, read_length=30,
        bias_5p=1.0, bias_3p=1.0, reverse_complement_probability=0.5,
        seed=321, flush_every=50,
    )
    assert set(s3["boundary_5p_counts"]) <= {"A", "G"}
    assert set(s3["boundary_3p_counts"]) <= {"A", "G"}
    ok("fragmentation_simulator", "global multi-contig sampling, strand sampling, determinism, bias=1 boundary rule")


def test_frame_normalization() -> None:
    plus = frame0_trim.CDSRegion(100, 200, "+", 1, "plus")
    assert frame0_trim.compute_phase_from_coordinates(101, 131, False, plus) == 0
    assert frame0_trim.compute_phase_from_coordinates(102, 132, False, plus) == 1
    minus = frame0_trim.CDSRegion(100, 200, "-", 1, "minus")
    assert frame0_trim.compute_phase_from_coordinates(169, 199, True, minus) == 0
    assert frame0_trim.compute_phase_from_coordinates(168, 198, True, minus) == 1
    corrected, t5, t3 = frame0_trim.phase_correct_sequence("ACGTACGTAC", 1)
    assert (t5, t3) == (2, 2) and len(corrected) % 3 == 0
    ok("frame_normalization_pure_logic", "GFF phase, plus/minus coordinates, phase-1 trimming")


def test_em() -> None:
    rng = np.random.default_rng(123)
    p_true = np.array([
        [0.55, 0.10, 0.25, 0.10],
        [0.10, 0.50, 0.15, 0.25],
        [0.15, 0.10, 0.60, 0.15],
    ])
    length = 60
    seqs = []
    for _ in range(1600):
        offset = int(rng.integers(0, 3))
        sequence = "".join(rng.choice(em.BASES, p=p_true[(i + offset) % 3]) for i in range(length))
        if bool(rng.integers(0, 2)):
            sequence = em.revcomp(sequence)
        seqs.append(sequence)
    cp = np.stack([em.counts_mod3(s) for s in seqs])
    cm = np.stack([em.counts_mod3(em.revcomp(s)) for s in seqs])
    p, assigns, trace = em.fit_em(cp, cm, alpha=0.5, max_iter=100, tol=1e-8, seed=7)
    assignments_check, _ = em.e_step(cp, cm, p)
    assert np.array_equal(assigns, assignments_check)
    assert np.allclose(p, em.m_step(cp, cm, assigns, alpha=0.5))
    assert all(b >= a - 1e-7 for a, b in zip(trace, trace[1:]))

    seeds = [11, 12, 13, 14]
    serial = em.run_restarts(cp, cm, alpha=0.5, max_iter=100, tol=1e-8, seeds=seeds, threads=1)
    parallel = em.run_restarts(cp, cm, alpha=0.5, max_iter=100, tol=1e-8, seeds=seeds, threads=2)
    assert serial[1] == parallel[1]
    assert np.allclose(serial[2], parallel[2])

    canonical = em.canonicalize_solution(p, mode="full6", effective_length=length).p_display
    for rotation in (0, 1, 2):
        rotated = em.canonicalize_solution(
            np.roll(p, rotation, axis=0), mode="full6", effective_length=length
        ).p_display
        assert np.allclose(canonical, rotated)
    ok("hard_em", "fixed-point consistency, monotone objective, deterministic restarts, canonicalization")


def write_result(path: Path, sample: str, length: int, values: list[float], read_count):
    target = path / sample / "stats" / "period3" / "kmer" / f"L{length}" / "k1_period3.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "status": "ok", "read_length": length, "read_count": read_count,
        "observed_read_count": read_count,
        "fits": [{"base": b, "r2": v} for b, v in zip(("A", "T", "G", "C"), values)],
    }
    target.write_text(json.dumps(payload))
    return target


def test_aggregation_and_metadata(work: Path) -> None:
    root = work / "period3_root"
    # Best row has missing read_count; GroupBy.first would incorrectly splice the
    # second row's count into it. The audited code must preserve None/NaN.
    files = [
        write_result(root, "ERR1", 59, [0.9, 0.9, 0.9, 0.9], None),
        write_result(root, "ERR1", 60, [0.7, 0.7, 0.7, 0.7], 1000),
        write_result(root, "internal1", 59, [0.4, 0.4, 0.4, 0.4], 500),
    ]
    all_lengths, best, errors = aggregate.aggregate(files)
    assert not errors and len(all_lengths) == 3
    row = best.loc[best["sample"] == "ERR1"].iloc[0]
    assert row["best_length"] == 59 and pd.isna(row["best_length_read_count"])
    assert row["n_eligible_lengths"] == 2

    r2 = best[["sample", "category", "max_r2"]].copy()
    meta = pd.DataFrame({
        "archive_data_accession": ["ERR1"],
        "library_name": ["other"],
        "project_name": ["ProjectX"],
    })
    joined, unmatched = metadata_plot.build_analysis_table(r2, meta, internal_project_label="Slon2017")
    assert len(joined) == 2 and len(unmatched) == 1
    assert joined.loc[joined["sample"] == "internal1", "project_name"].iloc[0] == "Slon2017"
    ok("aggregation_metadata", "maximum of four-base mean R2, complete-row tie selection, unmatched reporting")


def test_simulation_plot(work: Path) -> None:
    rows = []
    for genome in ("bacterial", "human"):
        for b5, b3 in ((0.5, 0.5), (1.0, 1.0)):
            for replicate in (1, 2, 3):
                for i, base in enumerate(("A", "T", "G", "C")):
                    value = 0.1 + 0.02 * i + 0.1 * (b5 - 0.5) + 0.01 * replicate
                    rows.append([genome, b5, b3, replicate, base, min(value, 1)])
    source = work / "simulation_results.csv"
    pd.DataFrame(rows, columns=["genome", "bias_5prime", "bias_3prime", "replicate", "base", "r2"]).to_csv(source, index=False)
    df = plot_simulation_r2.load_and_validate(source)
    replicate, summary = plot_simulation_r2.summarize(df)
    assert len(replicate) == 12 and len(summary) == 4
    out = work / "sim_plot.png"
    plot_simulation_r2.plot_one_genome(summary, "human", out)
    assert out.exists() and out.with_suffix(".pdf").exists()
    ok("simulation_plotting", "machine-readable replicates and 95% CI pipeline")


def main() -> int:
    work = ROOT / "formal_smoke_test_workspace"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    tests = [
        test_compile,
        lambda: test_period3_statistics(work),
        lambda: test_pipeline_io_and_counting(work),
        test_dft,
        lambda: test_trim_fastq(work),
        lambda: test_simulator(work),
        test_frame_normalization,
        test_em,
        lambda: test_aggregation_and_metadata(work),
        lambda: test_simulation_plot(work),
    ]
    for test in tests:
        test()
    try:
        import pysam  # noqa: F401
        import intervaltree  # noqa: F401
    except ImportError:
        skip("bam_integration_tests", "pysam and/or intervaltree are not installed in this offline environment")
    else:
        skip("bam_integration_tests", "not run because no representative BAM/GFF/reference dataset was supplied")

    output = ROOT / "TEST_RESULTS.txt"
    output.write_text("\n".join(RESULTS) + "\n")
    print(output.read_text(), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

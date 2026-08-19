#!/usr/bin/env python3
"""Frequency-domain analysis of read-cycle nucleotide fractions.

The publication DFT uses a 30-cycle one-based window (cycles 10-39 by
default), so the period-3 component lies exactly on Fourier bin q=10,
frequency 1/3.  Numerical cycle fractions and spectra are always exported.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from period3_library_pipeline import AnalysisConfig, detect_input_mode, iter_sequences_for_input

BASES: Tuple[str, ...] = ("A", "T", "G", "C")
BASE_TO_COL = {base: idx for idx, base in enumerate(BASES)}
COLORS = {
    "A": "#1f77b4",
    "T": "#ff7f0e",
    "G": "#2ca02c",
    "C": "#d62728",
}
SCRIPT_VERSION = "2.2.0-audited"


def trim_sequence(sequence: str, trim5: int, trim3: int) -> str:
    """Trim sequenced 5' and 3' ends; return an empty string if impossible."""

    if trim5 < 0 or trim3 < 0:
        raise ValueError("trim5 and trim3 must be non-negative.")
    if trim5 + trim3 >= len(sequence):
        return ""
    end = len(sequence) - trim3 if trim3 else len(sequence)
    return sequence[trim5:end]


def is_unambiguous_dna(sequence: str) -> bool:
    return bool(sequence) and all(base in BASE_TO_COL for base in sequence)


def get_positional_frequencies(
    sequences: Iterable[str],
    *,
    target_length: int,
    trim5: int = 0,
    trim3: int = 0,
    max_reads: int = 0,
    collapse_duplicates: bool = False,
) -> Tuple[Optional[np.ndarray], Dict[str, int]]:
    """Count A/T/G/C at each retained cycle for one exact raw read length."""

    if target_length <= 0:
        raise ValueError("target_length must be positive.")
    effective_length = target_length - trim5 - trim3
    if effective_length <= 0:
        raise ValueError("trim5 + trim3 must be smaller than target_length.")

    counts = np.zeros((effective_length, len(BASES)), dtype=np.uint64)
    seen: Optional[set[str]] = set() if collapse_duplicates else None
    stats = {
        "records_scanned": 0,
        "raw_length_matched": 0,
        "invalid_or_ambiguous": 0,
        "duplicates_collapsed": 0,
        "reads_analyzed": 0,
    }

    for raw_sequence in sequences:
        stats["records_scanned"] += 1
        sequence = raw_sequence.upper()
        if len(sequence) != target_length:
            continue
        stats["raw_length_matched"] += 1
        retained = trim_sequence(sequence, trim5, trim3)
        if not is_unambiguous_dna(retained):
            stats["invalid_or_ambiguous"] += 1
            continue
        if seen is not None:
            if retained in seen:
                stats["duplicates_collapsed"] += 1
                continue
            seen.add(retained)

        encoded = np.frombuffer(retained.encode("ascii"), dtype="S1")
        counts[:, 0] += encoded == b"A"
        counts[:, 1] += encoded == b"T"
        counts[:, 2] += encoded == b"G"
        counts[:, 3] += encoded == b"C"
        stats["reads_analyzed"] += 1
        if max_reads and stats["reads_analyzed"] >= max_reads:
            break

    n_reads = stats["reads_analyzed"]
    if n_reads == 0:
        return None, stats
    frequencies = counts.astype(np.float64) / float(n_reads)
    return frequencies, stats


def select_window(
    cycle_frequencies: np.ndarray, start_cycle: int, end_cycle: int
) -> np.ndarray:
    """Select an inclusive one-based cycle window."""

    if start_cycle < 1 or end_cycle < start_cycle:
        raise ValueError("Require 1 <= window_start <= window_end.")
    if end_cycle > cycle_frequencies.shape[0]:
        raise ValueError(
            f"Window ends at cycle {end_cycle}, but retained reads have only "
            f"{cycle_frequencies.shape[0]} cycles."
        )
    return cycle_frequencies[start_cycle - 1 : end_cycle, :]


def compute_spectrum(window_frequencies: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return rFFT frequencies and |F|/N magnitudes for all four bases."""

    if window_frequencies.ndim != 2 or window_frequencies.shape[1] != len(BASES):
        raise ValueError("window_frequencies must have shape (cycles, 4).")
    n_cycles = window_frequencies.shape[0]
    if n_cycles < 3:
        raise ValueError("At least three cycles are required.")
    frequencies = np.fft.rfftfreq(n_cycles, d=1.0)
    magnitudes = np.abs(np.fft.rfft(window_frequencies, axis=0)) / n_cycles
    return frequencies, magnitudes


def target_bin(frequencies: np.ndarray, target_frequency: float = 1 / 3) -> int:
    """Return an exact target bin; fail rather than silently using a nearby bin."""

    matches = np.flatnonzero(np.isclose(frequencies, target_frequency, atol=1e-12, rtol=0.0))
    if len(matches) != 1:
        raise ValueError(
            f"Frequency {target_frequency:g} is not an exact DFT bin for this window. "
            "Use a cycle count divisible by three (30 cycles in the manuscript)."
        )
    return int(matches[0])


def spectral_r2_at_target(
    window_frequencies: np.ndarray,
    magnitudes: np.ndarray,
    target_index: int,
) -> Dict[str, float]:
    """Fraction of variance carried by the target conjugate Fourier pair."""

    n_cycles = window_frequencies.shape[0]
    nyquist_index = n_cycles // 2 if n_cycles % 2 == 0 else None
    pair_factor = 1.0 if target_index == nyquist_index else 2.0
    result: Dict[str, float] = {}
    for column, base in enumerate(BASES):
        variance = float(np.var(window_frequencies[:, column]))
        if variance <= 0.0:
            result[base] = 0.0
        else:
            value = pair_factor * float(magnitudes[target_index, column] ** 2) / variance
            result[base] = float(np.clip(value, 0.0, 1.0))
    return result


def write_cycle_frequencies(path: Path, frequencies: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cycle", *BASES])
        for cycle, row in enumerate(frequencies, start=1):
            writer.writerow([cycle, *map(float, row)])


def write_spectrum(path: Path, frequencies: np.ndarray, magnitudes: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["frequency", *[f"magnitude_{base}" for base in BASES]])
        for frequency, row in zip(frequencies, magnitudes):
            writer.writerow([float(frequency), *map(float, row)])


def plot_dft_spectrum(
    frequencies: np.ndarray,
    magnitudes: np.ndarray,
    *,
    n_reads: int,
    target_length: int,
    output_path: Path,
    ymax: Optional[float] = None,
) -> None:
    """Plot non-zero rFFT magnitudes with a correctly labeled A/T/G/C mapping."""

    fig, ax = plt.subplots(figsize=(8, 6))
    positive = frequencies > 0
    for column, base in enumerate(BASES):
        ax.plot(
            frequencies[positive],
            magnitudes[positive, column],
            marker="o",
            linewidth=2,
            markersize=4,
            color=COLORS[base],
            label=base,
        )
    ax.axvline(1 / 3, color="black", linestyle="--", linewidth=2, label="Period-3 (f=1/3)")
    ax.set_xlabel("Frequency (cycles per read position)")
    ax.set_ylabel("Normalized magnitude |F(q)| / N")
    ax.set_title(f"DFT spectrum of nucleotide frequencies\nL={target_length}, n={n_reads:,}")
    ax.set_xlim(0, 0.5)
    if ymax is not None:
        if ymax <= 0:
            raise ValueError("ymax must be positive.")
        ax.set_ylim(0, ymax)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    fig.savefig(output_path.with_suffix(".pdf"))
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="DFT analysis of three-base periodicity.")
    parser.add_argument("--input", nargs="+", required=True, help="One BAM/FASTA/FASTQ or paired FASTQs")
    parser.add_argument("--input-type", choices=["auto", "bam", "fasta", "fastq"], default="auto")
    parser.add_argument("--length", type=int, required=True, help="Exact raw read length")
    parser.add_argument("--out-prefix", required=True, help="Output path prefix")
    parser.add_argument("--trim5", type=int, default=0)
    parser.add_argument("--trim3", type=int, default=0)
    parser.add_argument("--window-start", type=int, default=10)
    parser.add_argument("--window-end", type=int, default=39)
    parser.add_argument("--max-reads", type=int, default=0, help="Maximum analyzed reads; 0 means all")
    parser.add_argument("--collapse-exact-sequences", action="store_true")
    parser.add_argument("--exclude-duplicate-flag", action="store_true", help="BAM/CRAM only")
    parser.add_argument("--fq-which", choices=("r1", "r2", "both"), default="r1")
    parser.add_argument("--ymax", type=float, default=None, help="Optional fixed y-axis maximum for matched panels")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if len(args.input) > 2:
        raise ValueError("At most two input paths are supported; two paths are treated as paired FASTQ.")
    input_path = args.input[0] if len(args.input) == 1 else ",".join(args.input)
    if args.input_type == "auto":
        detected = detect_input_mode(input_path)
    else:
        detected = args.input_type.upper()
        if len(args.input) == 2 and detected != "FASTQ":
            raise ValueError("Two input paths are supported only for paired FASTQ input.")

    config = AnalysisConfig(
        input_path=input_path,
        collapse_duplicates=False,
        exclude_duplicate_flag=args.exclude_duplicate_flag,
        fq_which=args.fq_which,
        max_reads=0,
        verbose=True,
    )
    sequence_iterator = iter_sequences_for_input(config, detected)
    cycle_frequencies, stats = get_positional_frequencies(
        sequence_iterator,
        target_length=args.length,
        trim5=args.trim5,
        trim3=args.trim3,
        max_reads=args.max_reads,
        collapse_duplicates=args.collapse_exact_sequences,
    )
    if cycle_frequencies is None:
        print("No qualifying reads found.", file=sys.stderr)
        return 1

    window = select_window(cycle_frequencies, args.window_start, args.window_end)
    frequencies, magnitudes = compute_spectrum(window)
    index = target_bin(frequencies)
    spectral_r2 = spectral_r2_at_target(window, magnitudes, index)

    prefix = Path(args.out_prefix)
    stem = f"{prefix.name}_L{args.length}"
    base_dir = prefix.parent
    cycle_csv = base_dir / f"{stem}_cycle_frequencies.csv"
    spectrum_csv = base_dir / f"{stem}_dft_spectrum.csv"
    plot_png = base_dir / f"{stem}_dft_spectrum.png"
    metadata_json = base_dir / f"{stem}_dft_metadata.json"
    write_cycle_frequencies(cycle_csv, cycle_frequencies)
    write_spectrum(spectrum_csv, frequencies, magnitudes)
    plot_dft_spectrum(
        frequencies,
        magnitudes,
        n_reads=stats["reads_analyzed"],
        target_length=args.length,
        output_path=plot_png,
        ymax=args.ymax,
    )

    metadata = {
        "script": "dft_mod3.py",
        "script_version": SCRIPT_VERSION,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input": args.input,
        "input_type": detected,
        "raw_read_length": args.length,
        "trim5": args.trim5,
        "trim3": args.trim3,
        "effective_length": int(cycle_frequencies.shape[0]),
        "window_one_based_inclusive": [args.window_start, args.window_end],
        "n_window_cycles": int(window.shape[0]),
        "normalization": "abs(rfft)/N_cycles",
        "target_frequency": 1 / 3,
        "target_bin_index": index,
        "target_bin_frequency": float(frequencies[index]),
        "magnitude_at_one_third": {
            base: float(magnitudes[index, column]) for column, base in enumerate(BASES)
        },
        "spectral_r2_at_one_third": spectral_r2,
        "collapse_exact_sequences": bool(args.collapse_exact_sequences),
        "exclude_duplicate_flag": bool(args.exclude_duplicate_flag) if detected == "BAM" else None,
        "read_counts": stats,
        "outputs": {
            "cycle_frequencies_csv": str(cycle_csv),
            "spectrum_csv": str(spectrum_csv),
            "plot_png": str(plot_png),
            "plot_pdf": str(plot_png.with_suffix(".pdf")),
        },
    }
    metadata_json.parent.mkdir(parents=True, exist_ok=True)
    metadata_json.write_text(json.dumps(metadata, indent=2))

    print(f"Analyzed {stats['reads_analyzed']:,} reads.")
    for base in BASES:
        print(
            f"{base}: magnitude(f=1/3)={metadata['magnitude_at_one_third'][base]:.8g}; "
            f"spectral R^2={spectral_r2[base]:.6f}"
        )
    print(f"Wrote {plot_png}, {spectrum_csv}, and {metadata_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

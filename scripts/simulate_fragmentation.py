#!/usr/bin/env python3
"""Simulate fixed-length fragments under a purine-associated boundary model.

Candidate molecules are sampled globally across all FASTA contigs in proportion
to the number of valid start positions.  Each candidate is assigned a forward or
reverse-complement orientation, tested at the read-oriented 5' outside base and
3' inside terminal base, and accepted according to the documented purine rule.

The original implementation assigned one contig to each worker and generated a
fixed accepted quota per contig.  On multi-contig references this could omit most
chromosomes and, more subtly, could distort the accepted contig distribution.
This audited implementation samples every *candidate* from the global reference
distribution before rejection, so the conditional distribution of accepted
fragments follows the stated model exactly.  Generation is intentionally
single-process and deterministic; ``--threads`` is retained only for command-line
compatibility and is recorded in metadata.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Dict, Optional, Sequence, TextIO, Tuple

import numpy as np

SCRIPT_VERSION = "3.0.0-audited"
DNA = frozenset("ACGT")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def complement_base(base: str) -> str:
    return base.translate(COMPLEMENT)


def open_text(path: str, mode: str = "rt") -> TextIO:
    return gzip.open(path, mode) if path.lower().endswith(".gz") else open(path, mode)


def read_fasta(path: str) -> Dict[str, str]:
    """Read a FASTA while rejecting duplicate or empty identifiers."""

    sequences: Dict[str, str] = {}
    header: Optional[str] = None
    chunks: list[str] = []
    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(chunks).upper()
                    if not sequence:
                        raise ValueError(f"FASTA record '{header}' has an empty sequence.")
                    if header in sequences:
                        raise ValueError(f"Duplicate FASTA identifier: {header}")
                    sequences[header] = sequence
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"Empty FASTA identifier at line {line_number}.")
                header = fields[0]
                chunks = []
            else:
                if header is None:
                    raise ValueError("FASTA sequence encountered before the first header.")
                chunks.append(line)
        if header is not None:
            sequence = "".join(chunks).upper()
            if not sequence:
                raise ValueError(f"FASTA record '{header}' has an empty sequence.")
            if header in sequences:
                raise ValueError(f"Duplicate FASTA identifier: {header}")
            sequences[header] = sequence
    return sequences


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def pyrimidine_acceptance_probability(bias: float) -> float:
    """Return P(accept | C/T); P(accept | A/G) is one.

    ``bias`` parameterizes the purine:pyrimidine acceptance odds.  At 0.5,
    purines and pyrimidines are accepted equally.  At 1.0, only purines are
    accepted.  It is not, in general, the final fraction of accepted boundaries
    that are purines because that fraction also depends on reference composition.
    """

    if not (0.5 <= bias <= 1.0):
        raise ValueError("bias must lie in [0.5, 1.0].")
    return (1.0 - bias) / bias


def accept_break_after_base(base: str, bias: float, rng: np.random.Generator) -> bool:
    """Accept A/G with probability one and C/T with (1-bias)/bias."""

    if base not in DNA:
        return False
    if base in ("A", "G"):
        return True
    return bool(rng.random() < pyrimidine_acceptance_probability(bias))


def valid_start_count(sequence_length: int, read_length: int) -> int:
    """Number of candidates with an available read-oriented 5' outside base."""

    return max(0, sequence_length - read_length)


def build_candidate_distribution(
    sequences: Dict[str, str], read_length: int
) -> Tuple[list[str], np.ndarray, np.ndarray]:
    """Return contig names, valid-start counts, and cumulative candidate counts."""

    names: list[str] = []
    counts: list[int] = []
    for name, sequence in sequences.items():
        count = valid_start_count(len(sequence), read_length)
        if count > 0:
            names.append(name)
            counts.append(count)
    if not names:
        raise ValueError(
            "No contig is long enough for the requested read length and the 5' outside-base test."
        )
    count_array = np.asarray(counts, dtype=np.int64)
    cumulative = np.cumsum(count_array, dtype=np.int64)
    return names, count_array, cumulative


def _choose_contig(cumulative_counts: np.ndarray, rng: np.random.Generator) -> int:
    slot = int(rng.integers(0, int(cumulative_counts[-1])))
    return int(np.searchsorted(cumulative_counts, slot, side="right"))


def generate_reads(
    sequences: Dict[str, str],
    output_path: Path,
    *,
    num_reads: int,
    read_length: int,
    bias_5p: float,
    bias_3p: float,
    reverse_complement_probability: float,
    seed: int,
    flush_every: int,
) -> dict:
    """Generate accepted reads from a global candidate distribution."""

    names, valid_counts, cumulative_counts = build_candidate_distribution(sequences, read_length)
    rng = np.random.default_rng(seed)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    attempts = 0
    accepted = 0
    max_attempts = max(100_000, num_reads * 10_000)
    invalid_candidates = 0
    boundary_5 = Counter()
    boundary_3 = Counter()
    orientations = Counter()
    accepted_contigs = Counter()
    candidate_contigs = Counter()
    lines: list[str] = []

    with output_path.open("w") as output:
        while accepted < num_reads:
            attempts += 1
            if attempts > max_attempts:
                raise RuntimeError(
                    f"Accepted only {accepted}/{num_reads} reads after {max_attempts} attempts. "
                    "Check ambiguous reference content or reduce boundary bias."
                )

            contig_index = _choose_contig(cumulative_counts, rng)
            contig_name = names[contig_index]
            sequence = sequences[contig_name]
            max_start = int(valid_counts[contig_index])
            candidate_contigs[contig_name] += 1

            reverse = bool(rng.random() < reverse_complement_probability)
            if reverse:
                # start in [0, max_start-1]; sequence[end] is the base just
                # outside the read-oriented 5' end on the reverse strand.
                start = int(rng.integers(0, max_start))
            else:
                # start in [1, max_start]; sequence[start-1] is the base just
                # outside the read-oriented 5' end on the forward strand.
                start = int(rng.integers(1, max_start + 1))
            end = start + read_length
            fragment = sequence[start:end]

            if reverse:
                base_5p = complement_base(sequence[end])
                base_3p = complement_base(sequence[start])
                emitted = reverse_complement(fragment)
                orientation = "reverse_complement"
            else:
                base_5p = sequence[start - 1]
                base_3p = sequence[end - 1]
                emitted = fragment
                orientation = "forward_reference"

            if (
                base_5p not in DNA
                or base_3p not in DNA
                or any(base not in DNA for base in fragment)
            ):
                invalid_candidates += 1
                continue
            if not accept_break_after_base(base_5p, bias_5p, rng):
                continue
            if not accept_break_after_base(base_3p, bias_3p, rng):
                continue

            header = (
                f">read_{accepted}|chrom={contig_name}|start={start}|end={end}"
                f"|base_5p_outside_oriented={base_5p}"
                f"|base_3p_inside_oriented={base_3p}|orientation={orientation}"
            )
            lines.append(f"{header}\n{emitted}\n")
            boundary_5[base_5p] += 1
            boundary_3[base_3p] += 1
            orientations[orientation] += 1
            accepted_contigs[contig_name] += 1
            accepted += 1

            if len(lines) >= flush_every:
                output.writelines(lines)
                lines.clear()

        if lines:
            output.writelines(lines)

    return {
        "accepted_reads": accepted,
        "attempts": attempts,
        "acceptance_fraction": accepted / attempts if attempts else 0.0,
        "invalid_candidates": invalid_candidates,
        "boundary_5p_counts": {key: int(value) for key, value in sorted(boundary_5.items())},
        "boundary_3p_counts": {key: int(value) for key, value in sorted(boundary_3.items())},
        "orientation_counts": {key: int(value) for key, value in sorted(orientations.items())},
        "candidate_contig_counts": {
            key: int(value) for key, value in sorted(candidate_contigs.items())
        },
        "accepted_contig_counts": {
            key: int(value) for key, value in sorted(accepted_contigs.items())
        },
        "valid_start_counts": {
            name: int(count) for name, count in zip(names, valid_counts)
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Purine-associated fixed-length fragmentation simulator"
    )
    parser.add_argument("-i", "--input", required=True, help="Reference FASTA[.gz]")
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("-n", "--num-reads", type=int, required=True)
    parser.add_argument("-l", "--read-length", type=int, required=True)
    parser.add_argument("--bias-5prime", type=float, default=0.5)
    parser.add_argument("--bias-3prime", type=float, default=0.5)
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help=(
            "Retained for compatibility. The audited global rejection sampler is "
            "single-process so output is deterministic and contig-unbiased."
        ),
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--reverse-complement-probability",
        type=float,
        default=0.5,
        help="Candidate probability of a reverse-oriented molecule (default: 0.5)",
    )
    parser.add_argument(
        "--chunk-reads",
        type=int,
        default=100_000,
        help="Number of accepted FASTA records buffered before writing",
    )
    parser.add_argument("--output-name", default="simulated_reads.fasta")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.num_reads <= 0 or args.read_length <= 0:
        raise ValueError("num-reads and read-length must be positive.")
    if args.threads <= 0 or args.chunk_reads <= 0:
        raise ValueError("threads and chunk-reads must be positive.")
    if not (0.0 <= args.reverse_complement_probability <= 1.0):
        raise ValueError("reverse-complement-probability must lie in [0, 1].")
    pyrimidine_acceptance_probability(args.bias_5prime)
    pyrimidine_acceptance_probability(args.bias_3prime)

    sequences = read_fasta(args.input)
    if not sequences:
        raise ValueError("No sequences found in FASTA.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / args.output_name
    summary = generate_reads(
        sequences,
        output_path,
        num_reads=args.num_reads,
        read_length=args.read_length,
        bias_5p=args.bias_5prime,
        bias_3p=args.bias_3prime,
        reverse_complement_probability=args.reverse_complement_probability,
        seed=args.seed,
        flush_every=args.chunk_reads,
    )

    metadata = {
        "script": "simulate_fragmentation.py",
        "script_version": SCRIPT_VERSION,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_fasta": args.input,
        "output_fasta": str(output_path),
        "num_reads": args.num_reads,
        "read_length": args.read_length,
        "bias_5prime": args.bias_5prime,
        "bias_3prime": args.bias_3prime,
        "acceptance_rule": {
            "A_or_G": 1.0,
            "C_or_T_5prime": pyrimidine_acceptance_probability(args.bias_5prime),
            "C_or_T_3prime": pyrimidine_acceptance_probability(args.bias_3prime),
        },
        "bias_parameter_interpretation": (
            "purine-versus-pyrimidine relative acceptance parameter; not the final "
            "purine fraction unless candidate composition is balanced"
        ),
        "boundary_definition": {
            "5prime": (
                "base immediately outside the sequenced 5-prime end, expressed in "
                "read orientation"
            ),
            "3prime": (
                "terminal base inside the sequenced 3-prime end, expressed in read orientation"
            ),
        },
        "candidate_sampling": (
            "global across all contigs, proportional to valid candidate start positions, "
            "before boundary rejection"
        ),
        "reverse_orientation_candidate_probability": args.reverse_complement_probability,
        "seed": args.seed,
        "threads_requested": args.threads,
        "execution_mode": "single-process deterministic global rejection sampling",
        "chunk_reads": args.chunk_reads,
        **summary,
    }
    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    metadata_path.write_text(json.dumps(metadata, indent=2))
    print(f"Generated {summary['accepted_reads']:,} reads in {output_path}")
    print(f"Metadata written to {metadata_path}")
    if args.threads != 1:
        print(
            "NOTE: --threads is recorded but not used by the audited global sampler; "
            "this avoids thread-dependent contig bias."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

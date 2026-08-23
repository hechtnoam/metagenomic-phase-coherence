#!/usr/bin/env python3
"""Create phase-redistributed or fixed-trim controls.

For the random control, each unique exact-length input sequence is assigned
deterministically (given the seed) to a 5' trim of 0, 1, or 2 nucleotides,
with a complementary 3' trim so that the total trim is always two nucleotides.
All duplicate copies of the same original sequence therefore receive the same
trim offset. A fixed-trim mode is provided as the matched control.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Optional, Sequence, TextIO, Tuple

from Bio import SeqIO

try:
    import pysam  # type: ignore
except ImportError:  # FASTA/FASTQ operation remains available
    pysam = None  # type: ignore[assignment]

SCRIPT_VERSION = "2.3.0-audited-sequence-level-randomization"
TOTAL_TRIM = 2
COMPLEMENT_TRANS = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT_TRANS)[::-1]


def trim_bounds(length: int, five_trim: int) -> Tuple[int, int]:
    """Return Python slice bounds for a constant two-base total trim."""

    if five_trim not in (0, 1, 2):
        raise ValueError("five_trim must be 0, 1, or 2.")
    if length <= TOTAL_TRIM:
        raise ValueError(f"Sequence length must exceed {TOTAL_TRIM}.")
    three_trim = TOTAL_TRIM - five_trim
    end = length - three_trim if three_trim else length
    return five_trim, end


def trim_sequence(sequence: str, five_trim: int) -> str:
    start, end = trim_bounds(len(sequence), five_trim)
    return sequence[start:end]


def open_text(path: str, mode: str) -> TextIO:
    if path.lower().endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def choose_five_trim(sequence: str, seed: int, mode: str, fixed_five_trim: int) -> int:
    """Choose the trim offset.

    In random mode the assignment is a deterministic function of the original
    sequence and seed. Consequently, exact duplicate input sequences always
    receive the same trim offset, preventing randomization from splitting one
    duplicate class into several distinct post-trim sequences.
    """

    if mode == "fixed":
        return fixed_five_trim
    payload = f"{seed}\0{sequence}".encode("ascii")
    digest = hashlib.blake2b(payload, digest_size=8, person=b"phase-trim").digest()
    return int.from_bytes(digest, byteorder="big", signed=False) % 3


def process_fastx(
    infile: str,
    outfile: str,
    target_length: int,
    file_format: str,
    seed: int,
    mode: str = "random",
    fixed_five_trim: int = 0,
) -> dict:
    """Trim FASTA/FASTQ records; SeqRecord slicing preserves FASTQ qualities."""

    stats = Counter(records_seen=0, length_matched=0, ambiguous_excluded=0, records_written=0)
    offsets = Counter()
    with open_text(infile, "r") as in_handle, open_text(outfile, "w") as out_handle:
        for record in SeqIO.parse(in_handle, file_format):
            stats["records_seen"] += 1
            if len(record.seq) != target_length:
                continue
            stats["length_matched"] += 1
            sequence = str(record.seq).upper()
            if any(base not in "ACGT" for base in sequence):
                stats["ambiguous_excluded"] += 1
                continue
            five_trim = choose_five_trim(sequence, seed, mode, fixed_five_trim)
            start, end = trim_bounds(len(record.seq), five_trim)
            trimmed_record = record[start:end]
            # Keep the original identifier exactly; add the selected offset to the description.
            trimmed_record.id = record.id
            trimmed_record.name = record.name
            trimmed_record.description = f"{record.description} phase_trim5={five_trim} trim3={TOTAL_TRIM-five_trim}"
            SeqIO.write(trimmed_record, out_handle, file_format)
            stats["records_written"] += 1
            offsets[five_trim] += 1
    return {**stats, "offset_counts": {str(k): int(offsets[k]) for k in (0, 1, 2)}}


def _copy_selected_tags(source, destination) -> None:
    """Preserve non-alignment metadata useful for tracing reads."""

    for tag in ("RG", "BC", "QT", "RX", "QX", "CB", "CR", "CY", "UB", "UR", "UY"):
        if source.has_tag(tag):
            destination.set_tag(tag, source.get_tag(tag))


def process_bam(
    infile: str,
    outfile: str,
    target_length: int,
    mapped_only: bool = False,
    min_mapq: Optional[int] = None,
    seed: int = 42,
    mode: str = "random",
    fixed_five_trim: int = 0,
) -> dict:
    """Write valid unmapped BAM records after sequence/quality trimming."""

    if pysam is None:
        raise ImportError("BAM processing requires pysam.")
    stats = Counter(records_seen=0, length_matched=0, ambiguous_excluded=0, records_written=0)
    offsets = Counter()

    with pysam.AlignmentFile(infile, "rb") as bam_in, pysam.AlignmentFile(
        outfile, "wb", template=bam_in
    ) as bam_out:
        for read in bam_in.fetch(until_eof=True):
            stats["records_seen"] += 1
            if read.is_secondary or read.is_supplementary or read.is_qcfail:
                continue
            if mapped_only and read.is_unmapped:
                continue
            if min_mapq is not None and (read.is_unmapped or read.mapping_quality < min_mapq):
                continue
            if hasattr(read, "get_forward_sequence"):
                sequence = read.get_forward_sequence()
            else:
                sequence = read.query_sequence
                if sequence is not None and read.is_reverse:
                    sequence = reverse_complement(sequence)
            if sequence is None or len(sequence) != target_length:
                continue
            stats["length_matched"] += 1
            sequence = sequence.upper()
            if any(base not in "ACGT" for base in sequence):
                stats["ambiguous_excluded"] += 1
                continue
            five_trim = choose_five_trim(sequence, seed, mode, fixed_five_trim)
            start, end = trim_bounds(len(sequence), five_trim)
            new_sequence = sequence[start:end]
            if hasattr(read, "get_forward_qualities"):
                qualities = read.get_forward_qualities()
            else:
                qualities = read.query_qualities
                if qualities is not None and read.is_reverse:
                    qualities = qualities[::-1]
            new_qualities = qualities[start:end] if qualities is not None else None

            out_read = pysam.AlignedSegment(bam_out.header)
            out_read.query_name = read.query_name
            out_read.query_sequence = new_sequence
            if new_qualities is not None:
                out_read.query_qualities = new_qualities
            flag = 0x4  # read unmapped
            if read.is_paired:
                flag |= 0x1 | 0x8  # paired, mate unmapped
                if read.is_read1:
                    flag |= 0x40
                if read.is_read2:
                    flag |= 0x80
            out_read.flag = flag
            out_read.reference_id = -1
            out_read.reference_start = -1
            out_read.mapping_quality = 0
            out_read.cigar = None
            out_read.next_reference_id = -1
            out_read.next_reference_start = -1
            out_read.template_length = 0
            _copy_selected_tags(read, out_read)
            out_read.set_tag("ZT", f"trim5:{five_trim};trim3:{TOTAL_TRIM-five_trim}", value_type="Z")
            bam_out.write(out_read)
            stats["records_written"] += 1
            offsets[five_trim] += 1

    return {**stats, "offset_counts": {str(k): int(offsets[k]) for k in (0, 1, 2)}}


def detect_format(filename: str) -> str:
    lower = filename.lower()
    for suffix in (".bam",):
        if lower.endswith(suffix):
            return "bam"
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if lower.endswith(suffix):
            return "fastq"
    for suffix in (".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna"):
        if lower.endswith(suffix):
            return "fasta"
    raise ValueError(f"Unsupported input extension: {filename}")


def input_stem(path: str) -> str:
    """Remove a recognized sequence-file suffix without truncating internal dots."""

    name = Path(path).name
    for suffix in (
        ".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz", ".fna.gz",
        ".fastq", ".fq", ".fasta", ".fa", ".fna", ".bam",
    ):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return Path(path).stem


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Random phase redistribution or a matched fixed two-base trim."
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("-l", "--length", type=int, required=True, help="Exact input read length")
    parser.add_argument("--mode", choices=("random", "fixed"), default="random")
    parser.add_argument(
        "--fixed-five-trim",
        type=int,
        choices=(0, 1, 2),
        default=0,
        help="5' trim in fixed mode; 3' trim is 2 minus this value",
    )
    parser.add_argument("--mapped-only", action="store_true", help="BAM only")
    parser.add_argument("--min-mapq", type=int, default=None, help="BAM only")
    parser.add_argument("--seed", type=int, default=42)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.length <= TOTAL_TRIM:
        raise ValueError(f"--length must exceed {TOTAL_TRIM}.")
    fmt = detect_format(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    input_path = Path(args.input)
    label = "phase_randomized" if args.mode == "random" else f"fixed_trim5_{args.fixed_five_trim}"

    if fmt == "bam":
        outfile = outdir / f"{input_path.stem}.{label}.bam"
        stats = process_bam(
            args.input,
            str(outfile),
            args.length,
            mapped_only=args.mapped_only,
            min_mapq=args.min_mapq,
            seed=args.seed,
            mode=args.mode,
            fixed_five_trim=args.fixed_five_trim,
        )
    else:
        compressed = args.input.lower().endswith(".gz")
        suffix = ".fastq" if fmt == "fastq" else ".fasta"
        if compressed:
            suffix += ".gz"
        outfile = outdir / f"{input_stem(args.input)}.{label}{suffix}"
        stats = process_fastx(
            args.input,
            str(outfile),
            args.length,
            fmt,
            args.seed,
            mode=args.mode,
            fixed_five_trim=args.fixed_five_trim,
        )

    manifest = {
        "script": "phase_redistribution.py",
        "script_version": SCRIPT_VERSION,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input": args.input,
        "output": str(outfile),
        "input_format": fmt,
        "target_input_length": args.length,
        "output_length": args.length - TOTAL_TRIM,
        "mode": args.mode,
        "seed": args.seed,
        "fixed_five_trim": args.fixed_five_trim if args.mode == "fixed" else None,
        "total_trim": TOTAL_TRIM,
        "randomization_unit": "unique original sequence" if args.mode == "random" else "fixed offset for every sequence",
        "randomization_method": "BLAKE2b(seed + original sequence) modulo 3" if args.mode == "random" else None,
        "ambiguous_input_policy": "exclude any input read containing symbols outside A/C/G/T before offset sampling",
        "bam_filters": {
            "mapped_only": args.mapped_only,
            "min_mapq": args.min_mapq,
            "secondary_supplementary_qcfail_excluded": True,
            "sequence_orientation": "original sequenced orientation restored before trimming",
        },
        "counts": stats,
    }
    manifest_path = outdir / f"{outfile.name}.manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"Output written to: {outfile}")
    print(f"Manifest written to: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

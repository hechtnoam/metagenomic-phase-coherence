#!/usr/bin/env python3
"""Normalize sense CDS reads to codon phase 0 and coding orientation.

Only primary, simple (M/= /X-only), fully CDS-contained alignments are retained.
The GFF3 phase field is respected.  Reverse-strand reads are converted back to
their sequenced orientation, which is also the CDS coding orientation for sense
alignments, before trimming the biological 5' end.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import pysam  # type: ignore
except ImportError:
    pysam = None  # type: ignore[assignment]

try:
    from intervaltree import IntervalTree  # type: ignore
except ImportError:
    IntervalTree = None  # type: ignore[assignment]

SCRIPT_VERSION = "2.1.0-audited"
DNA = frozenset("ACGT")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class CDSRegion:
    start: int  # zero-based inclusive
    end: int  # zero-based exclusive
    strand: str
    phase: int
    identifier: str


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def parse_gff_attributes(text: str) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for field in text.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            result[key] = value
    return result


def load_cds_intervals(gff_path: str):
    if IntervalTree is None:
        raise ImportError("frame_normalize_cds.py requires intervaltree.")
    cds_by_contig = collections.defaultdict(IntervalTree)
    with open(gff_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            chrom, _source, _feature, start, end, _score, strand, phase_text, attrs = fields
            if strand not in ("+", "-"):
                continue
            if phase_text not in ("0", "1", "2", "."):
                raise ValueError(f"Invalid GFF phase '{phase_text}' at line {line_number}.")
            phase = 0 if phase_text == "." else int(phase_text)
            start0 = int(start) - 1
            end0 = int(end)
            attributes = parse_gff_attributes(attrs)
            identifier = attributes.get("ID") or attributes.get("Parent") or f"{chrom}:{start}-{end}:{strand}"
            region = CDSRegion(start0, end0, strand, phase, identifier)
            cds_by_contig[chrom].addi(start0, end0, region)
    return cds_by_contig


def read_is_sense(read, cds: CDSRegion) -> bool:
    return (not read.is_reverse and cds.strand == "+") or (read.is_reverse and cds.strand == "-")


def compute_phase_from_coordinates(
    reference_start: int,
    reference_end: int,
    is_reverse: bool,
    cds: CDSRegion,
) -> int:
    """Return the read biological-5' codon phase in {0,1,2}."""

    if (not is_reverse and cds.strand != "+") or (is_reverse and cds.strand != "-"):
        raise ValueError("Read and CDS are not in the same (sense) orientation.")
    if reference_start < cds.start or reference_end > cds.end:
        raise ValueError("Read is not fully contained in the CDS.")
    if cds.strand == "+":
        first_complete_codon_base = cds.start + cds.phase
        biological_5p = reference_start
        return int((biological_5p - first_complete_codon_base) % 3)
    first_complete_codon_base = cds.end - 1 - cds.phase
    biological_5p = reference_end - 1
    return int((first_complete_codon_base - biological_5p) % 3)


def phase_correct_sequence(sequence_in_coding_orientation: str, phase: int) -> Tuple[str, int, int]:
    """Trim biological 5' to phase 0, then trim 3' to complete codons."""

    if phase not in (0, 1, 2):
        raise ValueError("phase must be 0, 1, or 2.")
    trim_5p = (-phase) % 3
    if len(sequence_in_coding_orientation) <= trim_5p:
        return "", trim_5p, 0
    corrected = sequence_in_coding_orientation[trim_5p:]
    trim_3p = len(corrected) % 3
    if trim_3p:
        corrected = corrected[:-trim_3p]
    return corrected, trim_5p, trim_3p


def has_simple_cigar(read) -> bool:
    """Allow only aligned match/equal/mismatch operations; reject clips and indels."""

    return bool(read.cigartuples) and all(operation in (0, 7, 8) for operation, _length in read.cigartuples)


def choose_unique_cds(read, tree) -> Tuple[Optional[CDSRegion], str]:
    overlaps = tree.overlap(read.reference_start, read.reference_end)
    if not overlaps:
        return None, "no_cds_overlap"
    containing = [
        interval.data
        for interval in overlaps
        if read.reference_start >= interval.begin
        and read.reference_end <= interval.end
        and read_is_sense(read, interval.data)
    ]
    # De-duplicate identical annotations that may appear more than once.
    unique = {(c.start, c.end, c.strand, c.phase, c.identifier): c for c in containing}
    if not unique:
        return None, "not_fully_contained_or_antisense"
    if len(unique) > 1:
        return None, "ambiguous_multiple_cds"
    return next(iter(unique.values())), "ok"


def process_bam_to_fasta(
    bam_path: str,
    gff_path: str,
    output_fasta: str,
    *,
    min_mapq: int = 25,
    include_duplicate_flag: bool = False,
) -> dict:
    if pysam is None:
        raise ImportError("frame_normalize_cds.py requires pysam.")
    cds_trees = load_cds_intervals(gff_path)
    stats = collections.Counter()
    phase_counts = collections.Counter()
    trim5_counts = collections.Counter()
    output_length_counts = collections.Counter()

    with pysam.AlignmentFile(bam_path, "rb") as bam, open(output_fasta, "w") as output:
        for read in bam.fetch(until_eof=True):
            stats["records_seen"] += 1
            if read.is_unmapped:
                stats["unmapped"] += 1
                continue
            if read.is_secondary or read.is_supplementary or read.is_qcfail:
                stats["non_primary_or_qcfail"] += 1
                continue
            if read.is_duplicate and not include_duplicate_flag:
                stats["duplicate_flag_excluded"] += 1
                continue
            if read.mapping_quality < min_mapq:
                stats["low_mapq"] += 1
                continue
            if not has_simple_cigar(read):
                stats["complex_cigar_excluded"] += 1
                continue
            chrom = read.reference_name
            if chrom not in cds_trees:
                stats["contig_without_cds"] += 1
                continue
            cds, reason = choose_unique_cds(read, cds_trees[chrom])
            if cds is None:
                stats[reason] += 1
                continue

            stored_sequence = read.query_sequence
            if stored_sequence is None:
                stats["missing_sequence"] += 1
                continue
            # pysam query_sequence is reverse-complemented in reverse alignments;
            # get_forward_sequence restores the sequence as originally read (5'->3').
            if hasattr(read, "get_forward_sequence"):
                forward_sequence = read.get_forward_sequence()
                if forward_sequence is None:
                    stats["missing_sequence"] += 1
                    continue
                coding_sequence = forward_sequence.upper()
            else:
                coding_sequence = reverse_complement(stored_sequence) if read.is_reverse else stored_sequence.upper()
            if not coding_sequence or any(base not in DNA for base in coding_sequence):
                stats["ambiguous_sequence"] += 1
                continue

            phase = compute_phase_from_coordinates(
                read.reference_start, read.reference_end, read.is_reverse, cds
            )
            corrected, trim_5p, trim_3p = phase_correct_sequence(coding_sequence, phase)
            if len(corrected) < 3:
                stats["too_short_after_trim"] += 1
                continue

            output.write(
                f">{read.query_name}|phase={phase}|trim5={trim_5p}|trim3={trim_3p}"
                f"|cds={cds.identifier}|cds_strand={cds.strand}\n{corrected}\n"
            )
            stats["written"] += 1
            phase_counts[phase] += 1
            trim5_counts[trim_5p] += 1
            output_length_counts[len(corrected)] += 1

    return {
        "counts": dict(stats),
        "input_phase_counts": {str(k): int(phase_counts[k]) for k in (0, 1, 2)},
        "five_prime_trim_counts": {str(k): int(trim5_counts[k]) for k in (0, 1, 2)},
        "output_length_counts": {str(k): int(v) for k, v in sorted(output_length_counts.items())},
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Normalize sense CDS reads to codon phase 0.")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--output", required=True, help="Output FASTA")
    parser.add_argument("--min-mapq", type=int, default=25)
    parser.add_argument("--include-duplicate-flag", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.min_mapq < 0:
        raise ValueError("min-mapq must be non-negative.")
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary = process_bam_to_fasta(
        args.bam,
        args.gff,
        str(output_path),
        min_mapq=args.min_mapq,
        include_duplicate_flag=args.include_duplicate_flag,
    )
    metadata = {
        "script": "frame_normalize_cds.py",
        "script_version": SCRIPT_VERSION,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "bam": args.bam,
        "gff": args.gff,
        "output": str(output_path),
        "min_mapq": args.min_mapq,
        "include_duplicate_flag": args.include_duplicate_flag,
        "filters": {
            "primary_only": True,
            "qcfail_excluded": True,
            "simple_cigar_only": True,
            "fully_contained_in_unique_sense_CDS": True,
            "ambiguous_bases_excluded": True,
        },
        "gff_phase_used": True,
        "reverse_strand_handling": "converted to sequenced/coding orientation before 5-prime trimming",
        **summary,
    }
    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    metadata_path.write_text(json.dumps(metadata, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"Output written to {output_path}")
    print(f"Metadata written to {metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
period3_library_pipeline.py

Publication-ready k-mer read-position-composition pipeline for the article:
"Widespread phase-coherent three-base periodicity in metagenomic sequencing reads".

The pipeline intentionally preserves the plotting style and period-3 statistical
rationale of the original scripts while making the requested publication changes:

* k=1 is the default analysis; k=2 and k=3 are optional via --k.
* default length selection uses every observed read length with at least 40,000
  reads instead of a fixed list of read lengths.
* per-length read counts are written to the run manifest, read-length stats, and
  period-3 JSON outputs.
* FASTQ/FASTA/BAM streaming counts support bounded multiprocessing.
* reads containing symbols outside A/C/G/T are excluded and reported.
* exact-sequence collapsing is opt-in and recorded in the run manifest.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import os
import subprocess
import sys
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, Future, ProcessPoolExecutor, wait
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

try:
    import period3_statistics as period3
except Exception:  # pragma: no cover - handled at runtime with a clear message
    period3 = None  # type: ignore[assignment]

try:
    import pysam  # type: ignore
except Exception:  # pragma: no cover - pysam is optional
    pysam = None  # type: ignore[assignment]

ARTICLE_TITLE = "Widespread phase-coherent three-base periodicity in metagenomic sequencing reads"
SCRIPT_VERSION = "3.2.0-audited"

DEFAULT_BASE_OUTDIR = Path("period3_output")
DEFAULT_K_VALUES: Tuple[int, ...] = (1,)
DEFAULT_MIN_READ_COUNT = 40_000

P3_START = 10
P3_END = 40
P3_ALPHA = 0.05
P3_PERMUTATIONS = 10_000
P3_SEED = 12_345

PROGRESS_EVERY = 200_000
MAX_AUTO_WORKERS = 64
BATCH_TARGET_BASES = 6_000_000
MAX_INFLIGHT_FACTOR = 2
MIN_INFLIGHT_BATCHES = 4

BASES: Tuple[str, ...] = ("A", "T", "G", "C")
BASE_TO_INDEX = {"A": 0, "T": 1, "G": 2, "C": 3}
BASE_TO_INDEX_BYTES = {b"A": 0, b"T": 1, b"G": 2, b"C": 3}
COMPLEMENT_TRANS = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class AnalysisConfig:
    """Runtime configuration parsed from the command line."""

    input_path: str
    base_outdir: Path = DEFAULT_BASE_OUTDIR
    threads: int = 0
    collapse_duplicates: bool = False
    exclude_duplicate_flag: bool = False
    image_format: str = "png"
    plot_style: str = "dots"
    skip_existing: bool = False
    dry_run: bool = False
    verbose: bool = True
    fq_which: str = "r1"
    max_reads: int = 0
    workers: int = 0
    min_read_count: int = DEFAULT_MIN_READ_COUNT
    period3_enabled: bool = True
    p3_start: int = P3_START
    p3_end: int = P3_END
    p3_alpha: float = P3_ALPHA
    p3_permutations: int = P3_PERMUTATIONS
    p3_seed: int = P3_SEED


@dataclass(frozen=True)
class AnalysisPlan:
    """Concrete work plan after input type and read lengths have been resolved."""

    input_mode: str
    k_values: Tuple[int, ...]
    read_lengths: Tuple[int, ...]
    observed_length_counts: Mapping[int, int]
    total_scanned_reads: int
    explicit_lengths: bool
    invalid_sequence_count: int = 0


@dataclass(frozen=True)
class CountResult:
    """In-memory k-mer counts and read-count metadata from one analysis pass."""

    counts_by_k: Dict[int, Dict[int, np.ndarray]]
    analyzed_read_counts: Dict[int, int]
    total_raw_reads: int
    total_analyzed_reads: int
    effective_workers: int
    duplicate_collapse_applied: bool
    invalid_sequence_count: int = 0
    duplicate_sequence_count: int = 0


# ---------------------------------------------------------------------------
# Input recognition and sequence streaming
# ---------------------------------------------------------------------------


def is_fastq_path(path: str) -> bool:
    lower = path.lower()
    return lower.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq"))


def is_fasta_path(path: str) -> bool:
    lower = path.lower()
    return lower.endswith((".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna"))


def is_bam_path(path: str) -> bool:
    lower = path.lower()
    return lower.endswith((".bam", ".cram"))


def split_paired_fastq_argument(input_arg: str) -> Tuple[str, Optional[str]]:
    """Return (R1, R2) for an input argument; R2 is None for single-end data."""

    if "," not in input_arg:
        return input_arg, None
    r1, r2 = input_arg.split(",", 1)
    return r1.strip(), r2.strip()


# Backward-compatible helper alias.
split_fastq_arg = split_paired_fastq_argument


def detect_input_mode(input_arg: str) -> str:
    """Classify input as FASTQ, FASTA, or BAM."""

    r1, r2 = split_paired_fastq_argument(input_arg)
    if is_fastq_path(r1) or (r2 is not None and is_fastq_path(r2)):
        if r2 is not None and not (is_fastq_path(r1) and is_fastq_path(r2)):
            raise ValueError("Paired FASTQ input must be provided as R1.fastq[.gz],R2.fastq[.gz].")
        return "FASTQ"
    if is_fasta_path(r1) and r2 is None:
        return "FASTA"
    if is_bam_path(r1) and r2 is None:
        return "BAM"
    raise ValueError(
        "Unrecognized input type. Provide BAM/CRAM, FASTQ(.gz), FASTA(.gz), "
        "or paired FASTQ as R1,R2."
    )


def derive_sample_name(input_arg: str) -> str:
    """Create a stable sample name from FASTQ/FASTA/BAM input path(s)."""

    r1, r2 = split_paired_fastq_argument(input_arg)
    if is_bam_path(r1):
        return Path(r1).stem
    if is_fasta_path(r1):
        name = Path(r1).name
        for suffix in (".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna"):
            if name.lower().endswith(suffix):
                return name[: -len(suffix)]
        return Path(r1).stem

    def fastq_stem(path: str) -> str:
        name = Path(path).name
        for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if name.lower().endswith(suffix):
                name = name[: -len(suffix)]
                break
        for mate_suffix in ("_R1", "_R2", "_1", "_2", ".R1", ".R2", ".1", ".2"):
            if name.endswith(mate_suffix):
                name = name[: -len(mate_suffix)]
                break
        return name

    stem1 = fastq_stem(r1)
    if r2 is None:
        return stem1
    stem2 = fastq_stem(r2)
    shared_chars = 0
    for a, b in zip(stem1, stem2):
        if a != b:
            break
        shared_chars += 1
    return stem1[:shared_chars] if shared_chars >= 3 else stem1


sample_name_from_input = derive_sample_name


def open_text_stream(path: str) -> Iterator[str]:
    """Yield text lines from plain or gzipped sequence files."""

    if path.lower().endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
            yield from (line.rstrip("\r\n") for line in handle)
    else:
        with open(path, "rt", encoding="utf-8", errors="strict") as handle:
            yield from (line.rstrip("\r\n") for line in handle)


def iter_fastq_sequences(path: str) -> Iterator[str]:
    """Stream validated uppercase read sequences from one FASTQ file."""

    stream = open_text_stream(path)
    record_number = 0
    while True:
        try:
            header = next(stream)
        except StopIteration:
            break
        record_number += 1
        try:
            sequence = next(stream)
            plus = next(stream)
            quality = next(stream)
        except StopIteration as exc:
            raise ValueError(
                f"Truncated FASTQ record {record_number} in {path}."
            ) from exc
        if not header.startswith("@"):
            raise ValueError(
                f"FASTQ record {record_number} in {path} has no '@' header."
            )
        if not plus.startswith("+"):
            raise ValueError(
                f"FASTQ record {record_number} in {path} has no '+' separator."
            )
        sequence = sequence.strip().upper()
        if len(sequence) != len(quality):
            raise ValueError(
                f"FASTQ record {record_number} in {path} has sequence/quality "
                f"lengths {len(sequence)}/{len(quality)}."
            )
        if sequence:
            yield sequence


def iter_paired_fastq_sequences(r1: str, r2: str, which: str) -> Iterator[str]:
    """Stream R1, R2, or both mates from paired FASTQ files."""

    if which == "r1":
        yield from iter_fastq_sequences(r1)
        return
    if which == "r2":
        yield from iter_fastq_sequences(r2)
        return
    yield from iter_fastq_sequences(r1)
    yield from iter_fastq_sequences(r2)


def iter_fasta_sequences(path: str) -> Iterator[str]:
    """Stream validated uppercase sequences from one FASTA file."""

    current: List[str] = []
    saw_header = False
    for line_number, line in enumerate(open_text_stream(path), start=1):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            saw_header = True
            if current:
                sequence = "".join(current).upper()
                if sequence:
                    yield sequence
                current = []
        else:
            if not saw_header:
                raise ValueError(
                    f"FASTA sequence encountered before the first header at line {line_number} in {path}."
                )
            current.append(line)
    if current:
        sequence = "".join(current).upper()
        if sequence:
            yield sequence

def reverse_complement(sequence: str) -> str:
    """Return the reverse complement while preserving N/n symbols."""

    return sequence.translate(COMPLEMENT_TRANS)[::-1]


def original_bam_sequence(read) -> Optional[str]:
    """Return the original sequenced orientation of a BAM/CRAM record."""

    if hasattr(read, "get_forward_sequence"):
        sequence = read.get_forward_sequence()
    else:
        sequence = read.query_sequence
        if sequence is not None and read.is_reverse:
            sequence = reverse_complement(sequence)
    return sequence.upper() if sequence else None


def iter_bam_sequences_pysam(path: str, exclude_duplicate_flag: bool = False) -> Iterator[str]:
    """Stream query sequences from BAM/CRAM with pysam, including unmapped reads."""

    if pysam is None:
        raise RuntimeError("pysam is not available")
    alignment_file = pysam.AlignmentFile(path, "rb")
    try:
        for read in alignment_file.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary or read.is_qcfail:
                continue
            if exclude_duplicate_flag and read.is_duplicate:
                continue
            sequence = original_bam_sequence(read)
            if sequence:
                yield sequence
    finally:
        alignment_file.close()


def iter_bam_sequences_samtools(path: str, threads: int, exclude_duplicate_flag: bool = False) -> Iterator[str]:
    """Stream query sequences from BAM/CRAM using samtools view."""

    excluded_flags = 2816 + (1024 if exclude_duplicate_flag else 0)
    command = ["samtools", "view", "-F", str(excluded_flags)]
    if threads > 0:
        command.extend(["-@", str(threads)])
    command.append(path)

    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
    except FileNotFoundError as exc:
        raise RuntimeError("samtools is not available in PATH and pysam could not be used.") from exc

    assert process.stdout is not None
    assert process.stderr is not None
    try:
        for line in process.stdout:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            sequence = fields[9]
            if sequence and sequence != "*":
                sequence = sequence.upper()
                flag = int(fields[1])
                if flag & 0x10:
                    sequence = reverse_complement(sequence)
                yield sequence
    finally:
        stderr = process.stderr.read()
        process.stdout.close()
        process.stderr.close()
        return_code = process.wait()
        if return_code != 0:
            message = stderr.strip() or f"samtools view exited with code {return_code}"
            raise RuntimeError(message)


def iter_bam_sequences(path: str, threads: int, verbose: bool, exclude_duplicate_flag: bool = False) -> Iterator[str]:
    """Stream BAM/CRAM sequences via pysam when possible, otherwise samtools."""

    if pysam is not None:
        try:
            # Open once here so fallback only happens before any reads are yielded.
            alignment_file = pysam.AlignmentFile(path, "rb")
        except Exception as exc:
            if verbose:
                print(f"[BAM] pysam open failed ({exc}); falling back to samtools view.", file=sys.stderr)
        else:
            try:
                for read in alignment_file.fetch(until_eof=True):
                    if read.is_secondary or read.is_supplementary or read.is_qcfail:
                        continue
                    if exclude_duplicate_flag and read.is_duplicate:
                        continue
                    sequence = original_bam_sequence(read)
                    if sequence:
                        yield sequence
                return
            finally:
                alignment_file.close()

    yield from iter_bam_sequences_samtools(path, threads, exclude_duplicate_flag)


def iter_sequences_for_input(config: AnalysisConfig, input_mode: str) -> Iterator[str]:
    """Return a fresh sequence iterator for the configured input."""

    r1, r2 = split_paired_fastq_argument(config.input_path)
    if input_mode == "FASTQ":
        if r2 is None:
            yield from iter_fastq_sequences(r1)
        else:
            yield from iter_paired_fastq_sequences(r1, r2, config.fq_which)
    elif input_mode == "FASTA":
        yield from iter_fasta_sequences(r1)
    elif input_mode == "BAM":
        yield from iter_bam_sequences(
            r1, config.threads, config.verbose, config.exclude_duplicate_flag
        )
    else:
        raise ValueError(f"Unsupported input mode: {input_mode}")


def is_unambiguous_dna(sequence: str) -> bool:
    """Return True only for non-empty A/C/G/T sequences."""

    return bool(sequence) and all(base in BASE_TO_INDEX for base in sequence)


# ---------------------------------------------------------------------------
# Length selection and metadata
# ---------------------------------------------------------------------------


def format_count(value: int) -> str:
    return f"{value:,}"


pretty_count = format_count


def scan_read_lengths(config: AnalysisConfig, input_mode: str) -> Tuple[Dict[int, int], int, int]:
    """Count retained A/C/G/T-only read lengths in a streaming pre-pass."""

    counts: Counter[int] = Counter()
    total = 0
    invalid = 0
    tag = f"[SCAN:{input_mode}]"
    for sequence in iter_sequences_for_input(config, input_mode):
        total += 1
        if not is_unambiguous_dna(sequence):
            invalid += 1
        else:
            counts[len(sequence)] += 1
        if config.max_reads and total >= config.max_reads:
            break
        if config.verbose and total % PROGRESS_EVERY == 0:
            print(f"{tag} processed {format_count(total)} reads ...")

    if config.verbose:
        print(
            f"{tag} done: {format_count(total)} records; "
            f"{format_count(invalid)} with non-ACGT bases excluded; "
            f"{len(counts)} retained lengths."
        )
    return dict(counts), total, invalid


def select_lengths_by_read_count(length_counts: Mapping[int, int], min_read_count: int) -> Tuple[int, ...]:
    """Select every read length whose observed count passes the threshold."""

    return tuple(sorted(length for length, count in length_counts.items() if count >= min_read_count))


def describe_top_lengths(length_counts: Mapping[int, int], limit: int = 10) -> str:
    top = sorted(length_counts.items(), key=lambda item: (-item[1], item[0]))[:limit]
    return ", ".join(f"L{length}={format_count(count)}" for length, count in top) or "none"


def resolve_k_values(raw_values: Sequence[int]) -> Tuple[int, ...]:
    values = tuple(sorted(set(int(value) for value in raw_values)))
    invalid = [value for value in values if value not in (1, 2, 3)]
    if invalid:
        raise ValueError(f"Only k=1, k=2, and k=3 are supported; got {invalid}.")
    if not values:
        raise ValueError("At least one k-mer size must be requested.")
    return values


# ---------------------------------------------------------------------------
# Counting logic and multiprocessing
# ---------------------------------------------------------------------------


def allocate_count_tables(
    k_values: Sequence[int], read_lengths: Sequence[int]
) -> Tuple[Dict[int, Dict[int, np.ndarray]], Dict[int, int]]:
    """Allocate zero-filled count matrices for requested k values and lengths."""

    counts_by_k: Dict[int, Dict[int, np.ndarray]] = {}
    for k in k_values:
        counts_by_k[k] = {
            length: np.zeros((length - k + 1, 4**k), dtype=np.uint64)
            for length in read_lengths
            if length >= k
        }
    analyzed_read_counts = {length: 0 for length in read_lengths}
    return counts_by_k, analyzed_read_counts


def update_k1_counts(matrix: np.ndarray, sequence: str) -> None:
    """Update per-position 1-mer counts for one sequence."""

    bases = np.frombuffer(sequence.encode("ascii"), dtype="S1")
    matrix[:, 0] += bases == b"A"
    matrix[:, 1] += bases == b"T"
    matrix[:, 2] += bases == b"G"
    matrix[:, 3] += bases == b"C"


def update_kmer_counts(matrix: np.ndarray, sequence: str, k: int) -> None:
    """Update per-position k-mer counts for k=2 or k=3."""

    bases = np.frombuffer(sequence.encode("ascii"), dtype="S1")
    length = len(sequence)
    for start in range(length - k + 1):
        code = 0
        is_valid = True
        for offset in range(k):
            value = BASE_TO_INDEX_BYTES.get(bases[start + offset])
            if value is None:
                is_valid = False
                break
            code = (code << 2) | value
        if is_valid:
            matrix[start, code] += 1


def count_sequence_into_tables(
    sequence: str,
    counts_by_k: Dict[int, Dict[int, np.ndarray]],
    analyzed_read_counts: Dict[int, int],
) -> bool:
    """Count one sequence if its length is selected; return True when counted."""

    length = len(sequence)
    if length not in analyzed_read_counts:
        return False
    analyzed_read_counts[length] += 1
    if 1 in counts_by_k and length in counts_by_k[1]:
        update_k1_counts(counts_by_k[1][length], sequence)
    if 2 in counts_by_k and length in counts_by_k[2]:
        update_kmer_counts(counts_by_k[2][length], sequence, 2)
    if 3 in counts_by_k and length in counts_by_k[3]:
        update_kmer_counts(counts_by_k[3][length], sequence, 3)
    return True


def batch_size_for_lengths(read_lengths: Sequence[int]) -> int:
    """Choose a bounded batch size for multiprocessing payloads."""

    max_length = max(read_lengths, default=150)
    return max(64, min(4096, BATCH_TARGET_BASES // max(max_length, 1)))


def resolve_worker_count(requested_workers: int) -> int:
    """Resolve --workers; 0 means automatic CPU count capped for safety."""

    if requested_workers <= 0:
        return max(1, min(os.cpu_count() or 1, MAX_AUTO_WORKERS))
    return max(1, requested_workers)


def _count_sequence_batch(
    payload: Tuple[List[str], Tuple[int, ...], Tuple[int, ...]],
) -> Tuple[Dict[int, Dict[int, np.ndarray]], Dict[int, int], int, int]:
    """Worker entry point for multiprocessing sequence counting."""

    sequences, read_lengths, k_values = payload
    counts_by_k, analyzed_counts = allocate_count_tables(k_values, read_lengths)
    analyzed = 0
    for sequence in sequences:
        if count_sequence_into_tables(sequence, counts_by_k, analyzed_counts):
            analyzed += 1
    return counts_by_k, analyzed_counts, len(sequences), analyzed


def merge_batch_counts(
    total_counts_by_k: Dict[int, Dict[int, np.ndarray]],
    total_analyzed_counts: Dict[int, int],
    batch_result: Tuple[Dict[int, Dict[int, np.ndarray]], Dict[int, int], int, int],
) -> int:
    """Merge one worker result into process-level count matrices."""

    batch_counts_by_k, batch_analyzed_counts, _batch_size, batch_analyzed = batch_result
    for k, length_tables in batch_counts_by_k.items():
        for length, matrix in length_tables.items():
            total_counts_by_k[k][length] += matrix
    for length, count in batch_analyzed_counts.items():
        total_analyzed_counts[length] += int(count)
    return int(batch_analyzed)


def count_sequences_serial(
    sequence_iter: Iterator[str],
    read_lengths: Sequence[int],
    k_values: Sequence[int],
    collapse_duplicates: bool,
    max_reads: int,
    verbose: bool,
    log_tag: str,
) -> CountResult:
    """Count sequences in the current process, optionally collapsing duplicates."""

    counts_by_k, analyzed_counts = allocate_count_tables(k_values, read_lengths)
    selected_lengths = set(read_lengths)
    seen_by_length: Optional[Dict[int, set[str]]] = None
    if collapse_duplicates:
        seen_by_length = {length: set() for length in read_lengths}

    total_raw = 0
    total_analyzed = 0
    invalid_sequences = 0
    duplicate_sequences = 0
    for sequence in sequence_iter:
        total_raw += 1
        if not is_unambiguous_dna(sequence):
            invalid_sequences += 1
            if max_reads and total_raw >= max_reads:
                break
            continue
        length = len(sequence)
        if length in selected_lengths:
            if seen_by_length is not None:
                seen = seen_by_length[length]
                if sequence in seen:
                    duplicate_sequences += 1
                    if max_reads and total_raw >= max_reads:
                        break
                    if verbose and total_raw % PROGRESS_EVERY == 0:
                        print(f"{log_tag} processed {format_count(total_raw)} reads ...")
                    continue
                seen.add(sequence)
            if count_sequence_into_tables(sequence, counts_by_k, analyzed_counts):
                total_analyzed += 1

        if max_reads and total_raw >= max_reads:
            break
        if verbose and total_raw % PROGRESS_EVERY == 0:
            print(f"{log_tag} processed {format_count(total_raw)} reads ...")

    return CountResult(
        counts_by_k=counts_by_k,
        analyzed_read_counts={length: int(count) for length, count in analyzed_counts.items()},
        total_raw_reads=total_raw,
        total_analyzed_reads=total_analyzed,
        effective_workers=1,
        duplicate_collapse_applied=collapse_duplicates,
        invalid_sequence_count=invalid_sequences,
        duplicate_sequence_count=duplicate_sequences,
    )


def count_sequences_parallel(
    sequence_iter: Iterator[str],
    read_lengths: Sequence[int],
    k_values: Sequence[int],
    max_reads: int,
    verbose: bool,
    requested_workers: int,
    log_tag: str,
) -> CountResult:
    """Count selected sequences with a bounded ProcessPoolExecutor."""

    workers = resolve_worker_count(requested_workers)
    if workers <= 1:
        return count_sequences_serial(
            sequence_iter=sequence_iter,
            read_lengths=read_lengths,
            k_values=k_values,
            collapse_duplicates=False,
            max_reads=max_reads,
            verbose=verbose,
            log_tag=log_tag,
        )

    counts_by_k, analyzed_counts = allocate_count_tables(k_values, read_lengths)
    read_lengths_tuple = tuple(read_lengths)
    k_values_tuple = tuple(k_values)
    selected_lengths = set(read_lengths)
    batch_limit = batch_size_for_lengths(read_lengths)
    max_inflight = max(MIN_INFLIGHT_BATCHES, workers * MAX_INFLIGHT_FACTOR)

    futures: set[Future] = set()
    batch: List[str] = []
    total_raw = 0
    total_analyzed = 0
    invalid_sequences = 0

    def submit_batch(executor: ProcessPoolExecutor) -> None:
        nonlocal batch
        if not batch:
            return
        futures.add(executor.submit(_count_sequence_batch, (batch, read_lengths_tuple, k_values_tuple)))
        batch = []

    def drain_finished(block_until_one_finishes: bool) -> None:
        nonlocal total_analyzed
        if not futures:
            return
        if block_until_one_finishes:
            done, _pending = wait(futures, return_when=FIRST_COMPLETED)
        else:
            done = {future for future in futures if future.done()}
        for future in done:
            futures.remove(future)
            total_analyzed += merge_batch_counts(counts_by_k, analyzed_counts, future.result())

    with ProcessPoolExecutor(max_workers=workers) as executor:
        for sequence in sequence_iter:
            total_raw += 1
            if not is_unambiguous_dna(sequence):
                invalid_sequences += 1
            elif len(sequence) in selected_lengths:
                batch.append(sequence)
                if len(batch) >= batch_limit:
                    submit_batch(executor)
                    while len(futures) >= max_inflight:
                        drain_finished(block_until_one_finishes=True)
                    drain_finished(block_until_one_finishes=False)

            if max_reads and total_raw >= max_reads:
                break
            if verbose and total_raw % PROGRESS_EVERY == 0:
                print(f"{log_tag} processed {format_count(total_raw)} reads ...")

        submit_batch(executor)
        while futures:
            drain_finished(block_until_one_finishes=True)

    return CountResult(
        counts_by_k=counts_by_k,
        analyzed_read_counts={length: int(count) for length, count in analyzed_counts.items()},
        total_raw_reads=total_raw,
        total_analyzed_reads=total_analyzed,
        effective_workers=workers,
        duplicate_collapse_applied=False,
        invalid_sequence_count=invalid_sequences,
        duplicate_sequence_count=0,
    )


# ---------------------------------------------------------------------------
# Output paths, CSV export, plotting, and statistics
# ---------------------------------------------------------------------------


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def should_skip(*paths: Path, skip_existing: bool) -> bool:
    return skip_existing and all(path.exists() for path in paths)


def build_output_roots(config: AnalysisConfig, sample: str) -> Tuple[Path, Path, Path, Path, Path]:
    sample_root = config.base_outdir / sample
    csv_root = sample_root / "csv"
    plots_root = sample_root / "plots"
    stats_root = sample_root / "stats"
    period3_root = stats_root / "period3"
    if not config.dry_run:
        csv_root.mkdir(parents=True, exist_ok=True)
        plots_root.mkdir(parents=True, exist_ok=True)
        period3_root.mkdir(parents=True, exist_ok=True)
    return sample_root, csv_root, plots_root, stats_root, period3_root


def path_for_kmer_outputs(
    csv_root: Path,
    plots_root: Path,
    k: int,
    read_length: int,
    image_format: str,
) -> Tuple[Path, Path]:
    csv_path = csv_root / "kmer" / f"L{read_length}" / f"k{k}.csv"
    plot_path = plots_root / "kmer" / f"L{read_length}" / f"k{k}.{image_format}"
    return csv_path, plot_path


def kmer_tokens(k: int) -> List[str]:
    return ["".join(parts) for parts in itertools.product(BASES, repeat=k)]


# Backward-compatible helper alias.
_k_tokens = kmer_tokens


def kmer_palette(k: int, n: int):
    tab20 = list(plt.cm.tab20.colors)
    tab20b = list(plt.cm.tab20b.colors)
    tab20c = list(plt.cm.tab20c.colors)
    tab10 = list(plt.cm.tab10.colors)
    palette = tab20 + tab20b + tab20c + tab10
    if len(palette) < n:
        palette += [plt.cm.hsv(i / max(n, 1)) for i in range(n - len(palette))]
    return palette[:n]


_palette_for_k = kmer_palette


def legend_on_right(ax, title: str = "k-mer"):
    return ax.legend(
        title=title,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0,
        frameon=True,
    )


_legend_side = legend_on_right


def plot_k1(path: Path, fractions: np.ndarray, sample: str, read_length: int, n_reads: int, style: str = "dots") -> None:
    """Save the original-style 1-mer read-position-composition plot."""

    ensure_parent(path)
    x = np.arange(1, fractions.shape[0] + 1)
    fig, ax = plt.subplots(figsize=(10, 5))

    labels = ["A", "T", "G", "C"]

    if style == "dots":
        for column, label in enumerate(labels):
            ax.scatter(x, fractions[:, column], s=25, label=label)
    else:
        for column, label in enumerate(labels):
            ax.plot(x, fractions[:, column], label=label, linewidth=1.6)
    ax.set_xlabel("Position")
    ax.set_ylabel("Fraction")
    ax.set_ylim(0.07, 0.52) 

    ax.set_title(f"{sample} (n={format_count(n_reads)})")
    legend = legend_on_right(ax, title="Base")
    fig.tight_layout(rect=(0, 0, 0.82, 1))
    fig.savefig(path, bbox_extra_artists=(legend,), bbox_inches="tight")
    plt.close(fig)


def plot_kN(path: Path, fractions: np.ndarray, sample: str, k: int, read_length: int, n_reads: int, style: str = "dots") -> None:
    """Save the original-style 2-mer or 3-mer position-composition plot."""

    ensure_parent(path)
    x = np.arange(1, fractions.shape[0] + 1)
    n_series = fractions.shape[1]
    tokens = kmer_tokens(k)
    palette = kmer_palette(k, n_series)

    fig, ax = plt.subplots(figsize=(12, 6))
    handles: List[Line2D] = []
    for series_index in range(n_series):
        color = palette[series_index]
        if style == "dots":
            ax.scatter(x, fractions[:, series_index], s=8, color=color)
        else:
            ax.plot(x, fractions[:, series_index], linewidth=0.7, alpha=0.9, color=color)
        label = tokens[series_index] if series_index < len(tokens) else f"{k}-mer #{series_index + 1}"
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="None",
                label=label,
                markersize=5.5,
                markerfacecolor=color,
                markeredgecolor=color,
            )
        )

    ax.set_xlabel("Read position (start of k-mer)")
    ax.set_ylabel("Fraction")
    ax.set_title(f"{sample} — {k}-mer  —  (n={format_count(n_reads)})")
    ncol = 1 if k == 2 else 3
    legend = ax.legend(
        handles=handles,
        title=f"{k}-mer",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        frameon=True,
        ncol=ncol,
        prop={"size": 8 if k == 3 else 9},
    )
    fig.tight_layout(rect=(0, 0, 0.70 if k == 3 else 0.82, 1))
    fig.savefig(path, bbox_extra_artists=(legend,), bbox_inches="tight")
    plt.close(fig)


_plot_k1 = plot_k1
_plot_kN = plot_kN


def fractions_from_counts(counts: np.ndarray, denominator: int) -> np.ndarray:
    if int(denominator) <= 0:
        raise ValueError("Cannot compute position fractions with a non-positive read count.")
    fractions = counts.astype(float)
    fractions /= int(denominator)
    return fractions


_to_fraction_table = fractions_from_counts


def write_kmer_csv(csv_path: Path, fractions: np.ndarray, k: int) -> None:
    ensure_parent(csv_path)
    tokens = list(BASES) if k == 1 else kmer_tokens(k)
    with open(csv_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["position", *tokens])
        for position_index in range(fractions.shape[0]):
            writer.writerow([position_index + 1, *[fractions[position_index, j] for j in range(fractions.shape[1])]])


def period3_fit_to_json(fit) -> Dict[str, object]:
    return {
        "base": fit.base,
        "n_positions": fit.n_positions,
        "intercept": fit.intercept,
        "cos_coef": fit.cos_coef,
        "sin_coef": fit.sin_coef,
        "r2": fit.r2,
        "p_value": fit.p_value,
        "p_value_adjusted": getattr(fit, "p_value_adjusted", fit.p_value),
        "significant": fit.significant,
        "amplitude": fit.amplitude,
        "phase_degrees": fit.phase_degrees() if hasattr(fit, "phase_degrees") else fit.phase_deg(),
    }


def write_period3_statistics(
    csv_path: Path,
    read_length: int,
    analyzed_read_count: int,
    observed_read_count: int,
    stats_root: Path,
    config: AnalysisConfig,
) -> None:
    """Run period-3 stats for one k=1 CSV and include read-count metadata."""

    if not config.period3_enabled:
        return
    if period3 is None:
        if config.verbose:
            print(f"[P3:SKIP] period3_statistics.py is not available for {csv_path.name}")
        return

    output_root = stats_root / "kmer" / f"L{read_length}"
    json_out = output_root / "k1_period3.json"
    phase_plot = output_root / "k1_period3_phase_profile.png"
    strength_plot = output_root / "k1_period3_strength.png"

    if should_skip(json_out, phase_plot, strength_plot, skip_existing=config.skip_existing):
        if config.verbose:
            print(f"[P3:SKIP] {json_out} (exists)")
        return

    ensure_parent(json_out)
    if config.verbose:
        print(f"[P3] {csv_path.name} (L={read_length}, n={format_count(analyzed_read_count)}) ...")

    length_seed = int(config.p3_seed + read_length * 10_007)
    base_payload: Dict[str, object] = {
        "article": ARTICLE_TITLE,
        "csv": str(csv_path),
        "read_length": int(read_length),
        "read_count": int(analyzed_read_count),
        "observed_read_count": int(observed_read_count),
        "window": [config.p3_start, config.p3_end],
        "alpha": float(config.p3_alpha),
        "permutations": int(config.p3_permutations),
        "seed": length_seed,
        "period3_statistics_version": getattr(period3, "SCRIPT_VERSION", None),
        "period3_statistics_sha256": (
            hashlib.sha256(Path(period3.__file__).read_bytes()).hexdigest()
            if getattr(period3, "__file__", None)
            else None
        ),
    }

    try:
        result = period3.test_period3_from_csv(
            csv_path=str(csv_path),
            start=config.p3_start,
            end=config.p3_end,
            alpha=config.p3_alpha,
            permutations=config.p3_permutations,
            seed=length_seed,
        )
    except ValueError as exc:
        payload = {
            **base_payload,
            "status": "skipped",
            "reason": str(exc),
            "all_significant": False,
            "fits": [],
        }
        json_out.write_text(json.dumps(payload, indent=2))
        if config.verbose:
            print(f"[P3:SKIP] {csv_path.name}: {exc}")
        return

    payload = {
        **base_payload,
        "status": "ok",
        "window": list(result.window),
        "all_significant": bool(result.all_significant),
        "fits": [period3_fit_to_json(fit) for fit in result.fits],
    }
    json_out.write_text(json.dumps(payload, indent=2))

    df = period3.load_one_mer_csv(str(csv_path))
    pattern = period3.phase_means_by_base(df, config.p3_start, config.p3_end)
    period3.plot_phase_profile_single(
        pattern=pattern,
        window=(config.p3_start, config.p3_end),
        title=csv_path.name,
        save_to=phase_plot,
    )
    period3.plot_strength_bars(
        fits=result.fits,
        window=(config.p3_start, config.p3_end),
        title=csv_path.name,
        save_to=strength_plot,
    )
    if config.verbose:
        print(f"[P3:OK] {json_out.name}, {phase_plot.name}, {strength_plot.name}")


def write_read_length_stats(
    stats_root: Path,
    plan: AnalysisPlan,
    count_result: CountResult,
    config: AnalysisConfig,
) -> None:
    """Write per-length observed/analyzed read counts for the run."""

    stats_root.mkdir(parents=True, exist_ok=True)
    all_lengths = sorted(set(plan.observed_length_counts) | set(plan.read_lengths))
    rows = []
    for length in all_lengths:
        observed_count = int(plan.observed_length_counts.get(length, 0))
        analyzed_count = int(count_result.analyzed_read_counts.get(length, 0))
        rows.append(
            {
                "length": length,
                "observed_read_count": observed_count,
                "analyzed_read_count": analyzed_count,
                "processed": length in set(plan.read_lengths),
                "passed_min_read_count": observed_count >= config.min_read_count,
            }
        )

    csv_path = stats_root / "read_lengths.csv"
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "length",
                "observed_read_count",
                "analyzed_read_count",
                "processed",
                "passed_min_read_count",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    json_path = stats_root / "read_lengths.json"
    json_path.write_text(
        json.dumps(
            {
                "article": ARTICLE_TITLE,
                "input": config.input_path,
                "min_read_count": config.min_read_count,
                "total_scanned_reads": plan.total_scanned_reads,
                "total_raw_reads_counted": count_result.total_raw_reads,
                "total_analyzed_reads": count_result.total_analyzed_reads,
                "invalid_sequence_count": count_result.invalid_sequence_count,
                "duplicate_sequence_count": count_result.duplicate_sequence_count,
                "lengths": rows,
            },
            indent=2,
        )
    )


def write_run_manifest(
    sample_root: Path,
    sample: str,
    plan: AnalysisPlan,
    count_result: CountResult,
    config: AnalysisConfig,
) -> None:
    """Write reproducibility metadata for the full run."""

    selected_lengths = [
        {
            "length": int(length),
            "observed_read_count": int(plan.observed_length_counts.get(length, 0)),
            "analyzed_read_count": int(count_result.analyzed_read_counts.get(length, 0)),
        }
        for length in plan.read_lengths
    ]
    manifest = {
        "article": ARTICLE_TITLE,
        "script_version": SCRIPT_VERSION,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input": config.input_path,
        "input_mode": plan.input_mode,
        "sample": sample,
        "base_outdir": str(config.base_outdir.resolve()),
        "image_format": config.image_format,
        "plot_style": config.plot_style,
        "skip_existing": config.skip_existing,
        "collapse_duplicates": count_result.duplicate_collapse_applied,
        "exclude_duplicate_flag": config.exclude_duplicate_flag if plan.input_mode == "BAM" else None,
        "bam_sequence_orientation": (
            "original sequenced orientation restored with get_forward_sequence or FLAG 0x10"
            if plan.input_mode == "BAM"
            else None
        ),
        "fq_which": config.fq_which if plan.input_mode == "FASTQ" else None,
        "threads": config.threads if plan.input_mode == "BAM" else None,
        "count_workers": count_result.effective_workers,
        "max_reads": config.max_reads,
        "min_read_count": config.min_read_count,
        "explicit_lengths": plan.explicit_lengths,
        "k_values": list(plan.k_values),
        "read_lengths_processed": selected_lengths,
        "total_scanned_reads": plan.total_scanned_reads,
        "total_raw_reads_counted": count_result.total_raw_reads,
        "total_analyzed_reads": count_result.total_analyzed_reads,
        "invalid_sequence_count_scan": plan.invalid_sequence_count,
        "invalid_sequence_count_counting": count_result.invalid_sequence_count,
        "duplicate_sequence_count": count_result.duplicate_sequence_count,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python_version": sys.version,
        "period3": {
            "enabled": config.period3_enabled and period3 is not None,
            "start": config.p3_start,
            "end": config.p3_end,
            "alpha": config.p3_alpha,
            "permutations": config.p3_permutations,
            "base_seed": config.p3_seed,
            "per_length_seed_rule": "base_seed + read_length * 10007",
            "multiple_testing": "Holm adjustment across A/T/G/C",
            "statistics_module_version": getattr(period3, "SCRIPT_VERSION", None) if period3 is not None else None,
            "statistics_module_sha256": (
                hashlib.sha256(Path(period3.__file__).read_bytes()).hexdigest()
                if period3 is not None and getattr(period3, "__file__", None)
                else None
            ),
        },
        "outputs": {
            "read_length_stats_csv": "stats/read_lengths.csv",
            "read_length_stats_json": "stats/read_lengths.json",
        },
    }
    (sample_root / "run_manifest.json").write_text(json.dumps(manifest, indent=2))


# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------


def execute_analysis(config: AnalysisConfig, plan: AnalysisPlan) -> None:
    """Execute the concrete k-mer analysis plan."""

    sample = derive_sample_name(config.input_path)
    sample_root, csv_root, plots_root, stats_root, period3_root = build_output_roots(config, sample)

    if config.verbose:
        length_message = ", ".join(str(length) for length in plan.read_lengths)
        k_message = ", ".join(str(k) for k in plan.k_values)
        print(f"[{plan.input_mode}] sample: {sample}")
        print(f"[{plan.input_mode}] k values: {k_message}")
        print(f"[{plan.input_mode}] lengths: {length_message}")

    if config.dry_run:
        print("[DRY-RUN] Planned analysis only; no files were written.")
        return

    collapse_for_counting = config.collapse_duplicates

    use_parallel = not collapse_for_counting and resolve_worker_count(config.workers) > 1
    if use_parallel and config.verbose:
        print(
            f"[{plan.input_mode}] counting with {resolve_worker_count(config.workers)} worker processes "
            f"(~{batch_size_for_lengths(plan.read_lengths)} reads/batch; reader is single-threaded)."
        )
    elif collapse_for_counting and config.verbose:
        print(f"[{plan.input_mode}] duplicate-collapse enabled; counting in one process.")

    sequence_iter = iter_sequences_for_input(config, plan.input_mode)
    if use_parallel:
        count_result = count_sequences_parallel(
            sequence_iter=sequence_iter,
            read_lengths=plan.read_lengths,
            k_values=plan.k_values,
            max_reads=config.max_reads,
            verbose=config.verbose,
            requested_workers=config.workers,
            log_tag=f"[{plan.input_mode}]",
        )
    else:
        count_result = count_sequences_serial(
            sequence_iter=sequence_iter,
            read_lengths=plan.read_lengths,
            k_values=plan.k_values,
            collapse_duplicates=collapse_for_counting,
            max_reads=config.max_reads,
            verbose=config.verbose,
            log_tag=f"[{plan.input_mode}]",
        )

    if config.verbose:
        print(
            f"[{plan.input_mode}] done reading: {format_count(count_result.total_raw_reads)} reads; "
            f"{format_count(count_result.total_analyzed_reads)} reads used for selected lengths."
        )

    for k in plan.k_values:
        for read_length in plan.read_lengths:
            if read_length < k:
                if config.verbose:
                    print(f"[SKIP] k={k}, L={read_length}: read length is shorter than k.")
                continue
            if read_length not in count_result.counts_by_k.get(k, {}):
                continue
            analyzed_n = int(count_result.analyzed_read_counts.get(read_length, 0))
            observed_n = int(plan.observed_length_counts.get(read_length, 0))
            csv_path, plot_path = path_for_kmer_outputs(
                csv_root, plots_root, k, read_length, config.image_format
            )
            if should_skip(csv_path, plot_path, skip_existing=config.skip_existing):
                if config.verbose:
                    print(f"[SKIP] k={k}, L={read_length} (existing)")
                if k == 1 and csv_path.exists():
                    write_period3_statistics(
                        csv_path=csv_path,
                        read_length=read_length,
                        analyzed_read_count=analyzed_n,
                        observed_read_count=observed_n,
                        stats_root=period3_root,
                        config=config,
                    )
                continue

            fractions = fractions_from_counts(count_result.counts_by_k[k][read_length], analyzed_n)
            write_kmer_csv(csv_path, fractions, k)
            if k == 1:
                plot_k1(plot_path, fractions, sample, read_length, analyzed_n, style=config.plot_style)
                write_period3_statistics(
                    csv_path=csv_path,
                    read_length=read_length,
                    analyzed_read_count=analyzed_n,
                    observed_read_count=observed_n,
                    stats_root=period3_root,
                    config=config,
                )
            else:
                plot_kN(plot_path, fractions, sample, k, read_length, analyzed_n, style=config.plot_style)

            if config.verbose:
                print(
                    f"[WRITE] k={k}, L={read_length}, n={format_count(analyzed_n)}: "
                    f"{csv_path} | {plot_path}"
                )

    write_read_length_stats(stats_root, plan, count_result, config)
    write_run_manifest(sample_root, sample, plan, count_result, config)
    if config.verbose:
        print(f"[DONE] Outputs written under {sample_root}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", type=str, help="BAM/CRAM, FASTQ, FASTA, or paired FASTQ as R1,R2")
    parser.add_argument(
        "--base-outdir",
        type=Path,
        default=DEFAULT_BASE_OUTDIR,
        help=f"Base output directory (default: {DEFAULT_BASE_OUTDIR})",
    )
    parser.add_argument("--threads", type=int, default=10, help="samtools threads for BAM/CRAM fallback")
    parser.add_argument(
        "--exclude-duplicate-flag",
        action="store_true",
        help="For BAM/CRAM, exclude records carrying the SAM duplicate flag (off by default).",
    )
    parser.add_argument(
        "--collapse-exact-sequences",
        action="store_true",
        help="Collapse identical full read sequences within each exact-length class (off by default).",
    )
    parser.add_argument(
        "--no-collapse",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--format", choices=("png", "pdf"), default="png", help="Image format for plots")
    parser.add_argument("--plot-style", choices=("dots", "lines"), default="dots", help="Plot with dots or lines")
    parser.add_argument("--skip-existing", action="store_true", help="Reuse existing outputs without checking whether their parameters match")
    parser.add_argument("--force", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--dry-run", action="store_true", help="Resolve input mode and lengths without writing outputs")
    parser.add_argument("--quiet", action="store_true", help="Reduce log verbosity")
    parser.add_argument(
        "--fq-which",
        choices=("r1", "r2", "both"),
        default="r1",
        help="For paired FASTQ, choose which mate(s) to analyze",
    )
    parser.add_argument("--max-reads", type=int, default=0, help="Optional cap on streamed reads (0 = all)")
    parser.add_argument(
        "--workers",
        type=int,
        default=10,
        help="FASTQ/FASTA/BAM counting processes when duplicate-collapse is off: 0=auto, 1=serial",
    )
    parser.add_argument(
        "--min-read-count",
        type=int,
        default=DEFAULT_MIN_READ_COUNT,
        help="Default length filter: process observed lengths with at least this many reads",
    )
    parser.add_argument("--no-period3", action="store_true", help="Disable automatic period-3 statistics for k=1")
    parser.add_argument("--p3-start", type=int, default=P3_START, help="Start read position for period-3 window")
    parser.add_argument("--p3-end", type=int, default=P3_END, help="End read position for period-3 window")
    parser.add_argument("--p3-alpha", type=float, default=P3_ALPHA, help="Alpha for period-3 significance")
    parser.add_argument("--p3-perms", type=int, default=P3_PERMUTATIONS, help="Permutation count for p-values")
    parser.add_argument("--p3-seed", type=int, default=P3_SEED, help="Random seed for period-3 permutations")
    parser.add_argument(
        "--k",
        type=int,
        nargs="+",
        default=list(DEFAULT_K_VALUES),
        metavar="K",
        help="k-mer sizes to run. Default is k=1; add 2 and/or 3 explicitly, e.g. --k 1 2 3.",
    )
    parser.add_argument(
        "--lengths",
        type=int,
        nargs="+",
        default=None,
        metavar="L",
        help="Explicit read lengths. If omitted, all observed lengths with at least --min-read-count reads are processed.",
    )


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run k-mer position-composition and period-3 analyses for ancient/modern DNA libraries. "
            "drop-T analyses are intentionally not part of this publication pipeline."
        )
    )
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    all_parser = subparsers.add_parser(
        "all",
        help="Run the publication k-mer pipeline. Defaults: k=1 and all lengths with >=40,000 reads.",
    )
    add_common_arguments(all_parser)

    kmer_parser = subparsers.add_parser(
        "kmer",
        help="Alias of 'all' kept for users who prefer an explicit k-mer command.",
    )
    add_common_arguments(kmer_parser)

    return parser.parse_args(argv)


def build_config_and_plan(ns: argparse.Namespace) -> Tuple[AnalysisConfig, AnalysisPlan]:
    config = AnalysisConfig(
        input_path=ns.input,
        base_outdir=ns.base_outdir,
        threads=ns.threads,
        collapse_duplicates=bool(ns.collapse_exact_sequences) and not ns.no_collapse,
        exclude_duplicate_flag=bool(ns.exclude_duplicate_flag),
        image_format=ns.format,
        plot_style=ns.plot_style,
        skip_existing=bool(ns.skip_existing) and not ns.force,
        dry_run=ns.dry_run,
        verbose=not ns.quiet,
        fq_which=ns.fq_which,
        max_reads=ns.max_reads,
        workers=ns.workers,
        min_read_count=ns.min_read_count,
        period3_enabled=not ns.no_period3,
        p3_start=ns.p3_start,
        p3_end=ns.p3_end,
        p3_alpha=ns.p3_alpha,
        p3_permutations=ns.p3_perms,
        p3_seed=ns.p3_seed,
    )

    if config.max_reads < 0 or config.min_read_count < 1:
        raise ValueError("--max-reads must be non-negative and --min-read-count must be positive.")
    if config.p3_start < 1 or config.p3_end < config.p3_start:
        raise ValueError("Require 1 <= --p3-start <= --p3-end.")
    if not (0.0 < config.p3_alpha < 1.0) or config.p3_permutations < 1:
        raise ValueError("Require 0 < --p3-alpha < 1 and --p3-perms >= 1.")
    input_mode = detect_input_mode(config.input_path)
    observed_counts, total_scanned, invalid_scanned = scan_read_lengths(config, input_mode)
    k_values = resolve_k_values(ns.k)

    if ns.lengths is None:
        read_lengths = select_lengths_by_read_count(observed_counts, config.min_read_count)
        explicit_lengths = False
        if not read_lengths:
            raise ValueError(
                f"No read length has at least {format_count(config.min_read_count)} reads. "
                f"Top observed lengths: {describe_top_lengths(observed_counts)}. "
                "Use --min-read-count or --lengths to override."
            )
    else:
        requested_lengths = tuple(sorted(set(int(length) for length in ns.lengths)))
        explicit_lengths = True
        if not requested_lengths or any(length <= 0 for length in requested_lengths):
            raise ValueError("--lengths must contain one or more positive integers.")
        missing_lengths = [length for length in requested_lengths if observed_counts.get(length, 0) == 0]
        read_lengths = tuple(length for length in requested_lengths if observed_counts.get(length, 0) > 0)
        if missing_lengths and config.verbose:
            print(
                "[PLAN] requested lengths with zero retained A/C/G/T-only reads were skipped: "
                + ", ".join(map(str, missing_lengths))
            )
        if not read_lengths:
            raise ValueError(
                "None of the explicitly requested lengths had retained A/C/G/T-only reads. "
                f"Top observed lengths: {describe_top_lengths(observed_counts)}."
            )

    if config.verbose:
        if explicit_lengths:
            print(f"[PLAN] explicit lengths: {', '.join(map(str, read_lengths))}")
        else:
            print(
                f"[PLAN] selected {len(read_lengths)} lengths with >= "
                f"{format_count(config.min_read_count)} reads: {', '.join(map(str, read_lengths))}"
            )
        print(f"[PLAN] k values: {', '.join(map(str, k_values))}")

    plan = AnalysisPlan(
        input_mode=input_mode,
        k_values=k_values,
        read_lengths=read_lengths,
        observed_length_counts=observed_counts,
        total_scanned_reads=total_scanned,
        explicit_lengths=explicit_lengths,
        invalid_sequence_count=invalid_scanned,
    )
    return config, plan


# Backward-compatible function name for users of the old script internals.
def build_plan(ns: argparse.Namespace) -> Tuple[AnalysisConfig, AnalysisPlan]:
    return build_config_and_plan(ns)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    config, plan = build_config_and_plan(args)
    execute_analysis(config, plan)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

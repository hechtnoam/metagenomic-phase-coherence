#!/usr/bin/env python3
"""
shift_bias_mle_6state_v3_7_unified_canon_mp.py

Unified 6-state hard-EM model for BAM, FASTA, FASTQ, and gzipped FASTA/FASTQ.

What this script does
---------------------
- Reads sequence data from BAM or FASTA/FASTQ(.gz)
- Trims and optionally de-duplicates reads/sequences
- Fits the 6-state hard-EM model independently for each requested raw read length
- Supports multi-restart fitting with multiprocessing
- Reports per-length base-probability matrices, 6-state proportions, pooled 3-shift
  proportions, and hard-assignment gap statistics
- Optionally canonicalizes the displayed solution via:
    * none      : keep raw fitted labels
    * rotation  : resolve only cyclic 0/1/2 row-label symmetry
    * full6     : resolve cyclic symmetry plus strand-complement symmetry

Canonicalization summary
------------------------
The fitted 3x4 base-probability matrix is not fully identifiable from the model alone.
For a fixed effective read length L_eff (= raw length after trimming), there are two
sources of label symmetry:

1) Cyclic row symmetry:
   rotating the 3 position-class rows and relabeling offsets 0/1/2 gives an equivalent
   solution.

2) Strand-complement symmetry:
   because the model already contains both + and - latent orientations, there is a second
   equivalent branch obtained by reverse-complementing the displayed matrix in an
   L_eff mod 3 dependent way and relabeling the 6 states accordingly.

The 6-state shift plot follows the exact same relabeling applied to the displayed base-
probability matrix.

Dependencies
------------
- numpy, pandas, matplotlib, scipy
- pysam only if BAM input is used
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import multiprocessing as mp
import os
import time
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import chi2

SCRIPT_VERSION = "4.2.0-audited"

BASES = ["A", "C", "G", "T"]
B2I = {b: i for i, b in enumerate(BASES)}
STATE_NAMES = ["+0", "+1", "+2", "-0", "-1", "-2"]
STATE_TO_INDEX = {st: i for i, st in enumerate(STATE_NAMES)}
COMPL = str.maketrans("ACGTNacgtn", "TGCANtgcan")
COL_COMP_IDX = np.array([3, 2, 1, 0], dtype=np.int32)  # A<-T, C<-G, G<-C, T<-A
FASTA_SUFFIXES = (".fa", ".fasta", ".fna")
FASTQ_SUFFIXES = (".fq", ".fastq")

# Globals for multiprocessing workers
_G_CP = None
_G_CM = None
_G_ALPHA = None
_G_MAX_ITER = None
_G_TOL = None


@dataclass(frozen=True)
class CanonicalizationResult:
    p_display: np.ndarray
    state_map_old_to_new: np.ndarray
    shift_map_old_to_new: np.ndarray
    mode: str
    rotation: int
    rc_flip: bool
    effective_length: int
    effective_length_mod3: int


# -----------------------------------------------------------------------------
# Sequence I/O and preprocessing
# -----------------------------------------------------------------------------

def revcomp(seq: str) -> str:
    return seq.translate(COMPL)[::-1]


def open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="strict") if path.lower().endswith(".gz") else open(path, "rt", encoding="utf-8", errors="strict")


def strip_gz_suffix(path: str) -> str:
    return path[:-3] if path.lower().endswith(".gz") else path


def infer_kind_from_extension(path: str) -> Optional[str]:
    p = strip_gz_suffix(path.lower())
    if p.endswith(".bam"):
        return "bam"
    if p.endswith(FASTA_SUFFIXES):
        return "fasta"
    if p.endswith(FASTQ_SUFFIXES):
        return "fastq"
    return None


def sniff_text_sequence_kind(path: str) -> Optional[str]:
    try:
        with open_text(path) as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    return "fasta"
                if line.startswith("@"):
                    return "fastq"
                return None
    except Exception:
        return None
    return None


def detect_input_kind(paths: Sequence[str], forced_kind: str) -> str:
    if forced_kind != "auto":
        return forced_kind

    kinds = set()
    for path in paths:
        kind = infer_kind_from_extension(path)
        if kind is None:
            kind = sniff_text_sequence_kind(path)
        if kind is None:
            raise ValueError(f"Could not infer input format for: {path}")
        kinds.add(kind)

    if len(kinds) != 1:
        raise ValueError(f"Mixed input kinds are not supported in one run. Detected: {sorted(kinds)}")
    return next(iter(kinds))


def iter_fasta(paths: Sequence[str]) -> Iterator[str]:
    for fp in paths:
        with open_text(fp) as handle:
            header = None
            chunks: List[str] = []
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        yield "".join(chunks).upper()
                    header = line[1:].strip()
                    chunks = []
                else:
                    chunks.append(line)
            if header is not None:
                yield "".join(chunks).upper()


def iter_fastq(paths: Sequence[str]) -> Iterator[str]:
    for fp in paths:
        with open_text(fp) as handle:
            record_number = 0
            while True:
                header = handle.readline()
                if not header:
                    break
                record_number += 1
                sequence = handle.readline()
                plus = handle.readline()
                quality = handle.readline()
                if not quality:
                    raise ValueError(f"Truncated FASTQ record {record_number} in {fp}.")
                if not header.startswith("@") or not plus.startswith("+"):
                    raise ValueError(f"Malformed FASTQ record {record_number} in {fp}.")
                sequence = sequence.strip().upper()
                quality = quality.rstrip("\r\n")
                if len(sequence) != len(quality):
                    raise ValueError(
                        f"FASTQ record {record_number} in {fp} has sequence/quality "
                        f"lengths {len(sequence)}/{len(quality)}."
                    )
                if sequence:
                    yield sequence


def iter_bam(
    paths: Sequence[str],
    *,
    include_secondary: bool,
    include_supplementary: bool,
    include_qcfail: bool,
    include_duplicates_flag: bool,
) -> Iterator[str]:
    try:
        import pysam  # type: ignore
    except ImportError as exc:
        raise ImportError("BAM input requires pysam. Install it with `pip install pysam`.") from exc

    for fp in paths:
        with pysam.AlignmentFile(fp, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if (not include_secondary) and read.is_secondary:
                    continue
                if (not include_supplementary) and read.is_supplementary:
                    continue
                if (not include_qcfail) and read.is_qcfail:
                    continue
                if (not include_duplicates_flag) and read.is_duplicate:
                    continue
                if hasattr(read, "get_forward_sequence"):
                    seq = read.get_forward_sequence()
                else:
                    seq = read.query_sequence
                    if seq is not None and read.is_reverse:
                        seq = revcomp(seq)
                if seq is None:
                    continue
                yield seq.upper()


def iter_sequences(
    paths: Sequence[str],
    input_kind: str,
    *,
    include_secondary: bool,
    include_supplementary: bool,
    include_qcfail: bool,
    include_duplicates_flag: bool,
) -> Iterator[str]:
    if input_kind == "bam":
        yield from iter_bam(
            paths,
            include_secondary=include_secondary,
            include_supplementary=include_supplementary,
            include_qcfail=include_qcfail,
            include_duplicates_flag=include_duplicates_flag,
        )
    elif input_kind == "fasta":
        yield from iter_fasta(paths)
    elif input_kind == "fastq":
        yield from iter_fastq(paths)
    else:
        raise ValueError(f"Unsupported input kind: {input_kind}")


def trim_seq(seq: str, trim5: int, trim3: int) -> str:
    if trim5 + trim3 >= len(seq):
        return ""
    return seq[trim5 : len(seq) - trim3]


def dedup_key(seq: str, mode: str):
    if mode == "none":
        return None
    if mode == "sequence":
        return seq
    if mode == "seq_rc":
        rc = revcomp(seq)
        return seq if seq <= rc else rc
    if mode == "hash64":
        return hashlib.blake2b(seq.encode("ascii", "ignore"), digest_size=8).digest()
    if mode == "hash64_rc":
        rc = revcomp(seq)
        c = seq if seq <= rc else rc
        return hashlib.blake2b(c.encode("ascii", "ignore"), digest_size=8).digest()
    raise ValueError(f"Unknown dedup mode: {mode}")


def counts_mod3(seq: str) -> np.ndarray:
    counts = np.zeros((3, 4), dtype=np.uint32)
    for i, ch in enumerate(seq):
        j = B2I.get(ch)
        if j is None:
            continue
        counts[i % 3, j] += 1
    return counts


# -----------------------------------------------------------------------------
# EM model core
# -----------------------------------------------------------------------------

def init_p_from_global(cp: np.ndarray, alpha: float = 0.5, seed: int = 1) -> np.ndarray:
    """Initialize three rows near the global base composition.

    ``alpha`` is accepted for backward compatibility but is deliberately not
    added to already-normalized probabilities.  Pseudocounts belong in the
    M-step, where they regularize counts; adding 0.5 to probabilities here
    would almost erase the random row perturbations used to distinguish
    restarts.
    """

    _ = alpha
    rng = np.random.default_rng(seed)
    total_by_base = cp.sum(axis=(0, 1)).astype(np.float64)
    total = float(total_by_base.sum())
    if total > 0:
        base_freq = total_by_base / total
    else:
        base_freq = np.ones(4, dtype=np.float64) / 4.0
    p = np.tile(base_freq, (3, 1)).astype(np.float64)
    p += rng.normal(0.0, 0.01, size=p.shape)
    p = np.clip(p, 1e-8, None)
    p /= np.sum(p, axis=1, keepdims=True)
    return p


def score_all_states(cp: np.ndarray, cm: np.ndarray, p: np.ndarray) -> np.ndarray:
    logp = np.log(p)
    cols = []
    for o in (0, 1, 2):
        logp_o = np.roll(logp, -o, axis=0)
        cols.append(np.einsum("nrb,rb->n", cp, logp_o, optimize=True))
    for o in (0, 1, 2):
        logp_o = np.roll(logp, -o, axis=0)
        cols.append(np.einsum("nrb,rb->n", cm, logp_o, optimize=True))
    return np.stack(cols, axis=1)


def e_step(cp: np.ndarray, cm: np.ndarray, p: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    ll = score_all_states(cp, cm, p)
    assigns = np.argmax(ll, axis=1).astype(np.int32)
    best = np.max(ll, axis=1).astype(np.float64)
    return assigns, best


def m_step(cp: np.ndarray, cm: np.ndarray, assigns: np.ndarray, alpha: float = 0.5) -> np.ndarray:
    total = np.zeros((3, 4), dtype=np.float64)
    for state_idx, state_name in enumerate(STATE_NAMES):
        mask = assigns == state_idx
        if not np.any(mask):
            continue
        sign = state_name[0]
        offset = int(state_name[1])
        cnt = cp[mask] if sign == "+" else cm[mask]
        cnt_shifted = np.roll(cnt, +offset, axis=1)  # axis=1 is the 3-row position-class axis
        total += cnt_shifted.sum(axis=0)
    p = total + alpha
    p = p / np.sum(p, axis=1, keepdims=True)
    return p


def fit_em(
    cp: np.ndarray,
    cm: np.ndarray,
    *,
    alpha: float = 0.5,
    max_iter: int = 100,
    tol: float = 1e-6,
    seed: int = 1,
) -> Tuple[np.ndarray, np.ndarray, List[float]]:
    """Fit hard EM to a self-consistent assignment/parameter fixed point.

    Convergence is declared only when the deterministic hard assignments are
    unchanged after an M-step.  This ensures that the returned parameter matrix
    is the M-step estimate for the returned assignments *and* that those
    assignments are optimal under the returned matrix.  The likelihood
    tolerance is used to detect numerical decreases, not to stop at an
    assignment-inconsistent plateau.
    """

    if cp.shape != cm.shape or cp.ndim != 3 or cp.shape[1:] != (3, 4):
        raise ValueError("cp and cm must have matching shape (n_reads, 3, 4).")
    if cp.shape[0] == 0:
        raise ValueError("At least one read is required.")
    if alpha <= 0:
        raise ValueError("alpha must be positive so log probabilities remain finite.")
    if max_iter < 1 or tol < 0:
        raise ValueError("max_iter must be positive and tol non-negative.")

    p = init_p_from_global(cp, alpha=alpha, seed=seed)
    ll_trace: List[float] = []
    previous_assigns: Optional[np.ndarray] = None
    previous_ll = -np.inf

    for _iteration in range(max_iter):
        assigns, best = e_step(cp, cm, p)
        current_ll = float(np.sum(best))
        if np.isfinite(previous_ll) and current_ll + tol * (1.0 + abs(previous_ll)) < previous_ll:
            raise RuntimeError(
                "Hard-EM complete-data score decreased beyond numerical tolerance: "
                f"{previous_ll:.12g} -> {current_ll:.12g}."
            )
        ll_trace.append(current_ll)

        # p was produced by the M-step for previous_assigns. If the E-step
        # returns the same assignments, the pair is a true fixed point.
        if previous_assigns is not None and np.array_equal(assigns, previous_assigns):
            return p, assigns, ll_trace

        p = m_step(cp, cm, assigns, alpha=alpha)
        previous_assigns = assigns
        previous_ll = current_ll

    raise RuntimeError(
        "Hard EM did not reach assignment stability within "
        f"max_iter={max_iter}. Increase --max-iter or inspect tied states."
    )


def _init_worker(cp: np.ndarray, cm: np.ndarray, alpha: float, max_iter: int, tol: float) -> None:
    global _G_CP, _G_CM, _G_ALPHA, _G_MAX_ITER, _G_TOL
    _G_CP = cp
    _G_CM = cm
    _G_ALPHA = alpha
    _G_MAX_ITER = max_iter
    _G_TOL = tol


def _run_restart(seed: int):
    p, assigns, ll_trace = fit_em(_G_CP, _G_CM, alpha=_G_ALPHA, max_iter=_G_MAX_ITER, tol=_G_TOL, seed=seed)
    return ll_trace[-1], seed, p, assigns, ll_trace


def run_restarts(
    cp: np.ndarray,
    cm: np.ndarray,
    *,
    alpha: float,
    max_iter: int,
    tol: float,
    seeds: Sequence[int],
    threads: int,
) -> Tuple[float, int, np.ndarray, np.ndarray, List[float]]:
    best = None

    if threads <= 1 or len(seeds) == 1:
        for sd in seeds:
            p, assigns, ll_trace = fit_em(cp, cm, alpha=alpha, max_iter=max_iter, tol=tol, seed=sd)
            ll = ll_trace[-1]
            if (
                best is None
                or ll > best[0] + 1e-12
                or (abs(ll - best[0]) <= 1e-12 and sd < best[1])
            ):
                best = (ll, sd, p, assigns, ll_trace)
    else:
        start_methods = mp.get_all_start_methods()
        ctx = mp.get_context("fork") if "fork" in start_methods else mp.get_context()
        with ctx.Pool(
            processes=min(threads, len(seeds)),
            initializer=_init_worker,
            initargs=(cp, cm, alpha, max_iter, tol),
        ) as pool:
            for ll, sd, p, assigns, ll_trace in pool.imap_unordered(_run_restart, seeds, chunksize=1):
                if (
                    best is None
                    or ll > best[0] + 1e-12
                    or (abs(ll - best[0]) <= 1e-12 and sd < best[1])
                ):
                    best = (ll, sd, p, assigns, ll_trace)

    if best is None:
        raise RuntimeError("No EM restart completed successfully.")
    return best


# -----------------------------------------------------------------------------
# Canonicalization helpers
# -----------------------------------------------------------------------------

def lex_key_for_p(p: np.ndarray) -> Tuple[float, ...]:
    return tuple(np.round(p, 12).reshape(-1).tolist())


def rotation_state_map(rotation: int) -> np.ndarray:
    out = np.empty(6, dtype=np.int32)
    for old_idx, st in enumerate(STATE_NAMES):
        sign = st[0]
        offset = int(st[1])
        new_offset = (offset - rotation) % 3
        out[old_idx] = STATE_TO_INDEX[f"{sign}{new_offset}"]
    return out


def rotation_shift_map(rotation: int) -> np.ndarray:
    return np.array([(offset - rotation) % 3 for offset in range(3)], dtype=np.int32)


def rc_base_state_map() -> np.ndarray:
    out = np.empty(6, dtype=np.int32)
    for old_idx, st in enumerate(STATE_NAMES):
        sign = st[0]
        offset = int(st[1])
        new_sign = "-" if sign == "+" else "+"
        new_offset = (-offset) % 3
        out[old_idx] = STATE_TO_INDEX[f"{new_sign}{new_offset}"]
    return out


def rc_base_shift_map() -> np.ndarray:
    return np.array([(-offset) % 3 for offset in range(3)], dtype=np.int32)


def compose_maps(first_old_to_mid: np.ndarray, second_mid_to_new: np.ndarray) -> np.ndarray:
    return second_mid_to_new[first_old_to_mid]


def rc_branch_transform(p: np.ndarray, effective_length_mod3: int) -> np.ndarray:
    # For displayed row q, take source row (L_eff - 1 - q) mod 3 and complement columns.
    row_src = [((effective_length_mod3 - 1) - q) % 3 for q in range(3)]
    return p[row_src][:, COL_COMP_IDX]


def canonicalize_solution(
    p: np.ndarray,
    *,
    mode: str,
    effective_length: int,
) -> CanonicalizationResult:
    if mode not in {"none", "rotation", "full6"}:
        raise ValueError(f"Unknown canonicalization mode: {mode}")

    identity_state = np.arange(6, dtype=np.int32)
    identity_shift = np.arange(3, dtype=np.int32)
    eff_mod3 = effective_length % 3

    if mode == "none":
        return CanonicalizationResult(
            p_display=p.copy(),
            state_map_old_to_new=identity_state,
            shift_map_old_to_new=identity_shift,
            mode=mode,
            rotation=0,
            rc_flip=False,
            effective_length=effective_length,
            effective_length_mod3=eff_mod3,
        )

    candidates = []

    for rotation in (0, 1, 2):
        p_rot = np.roll(p, -rotation, axis=0)
        candidates.append(
            (
                lex_key_for_p(p_rot),
                p_rot,
                rotation_state_map(rotation),
                rotation_shift_map(rotation),
                rotation,
                False,
            )
        )

    if mode == "full6":
        p_rc_base = rc_branch_transform(p, eff_mod3)
        state_map_rc = rc_base_state_map()
        shift_map_rc = rc_base_shift_map()
        for rotation in (0, 1, 2):
            p_rc_rot = np.roll(p_rc_base, -rotation, axis=0)
            state_map = compose_maps(state_map_rc, rotation_state_map(rotation))
            shift_map = compose_maps(shift_map_rc, rotation_shift_map(rotation))
            candidates.append(
                (
                    lex_key_for_p(p_rc_rot),
                    p_rc_rot,
                    state_map,
                    shift_map,
                    rotation,
                    True,
                )
            )

    best_key, best_p, best_state_map, best_shift_map, best_rotation, best_rc_flip = min(
        candidates, key=lambda x: x[0]
    )
    _ = best_key
    return CanonicalizationResult(
        p_display=best_p,
        state_map_old_to_new=best_state_map,
        shift_map_old_to_new=best_shift_map,
        mode=mode,
        rotation=best_rotation,
        rc_flip=best_rc_flip,
        effective_length=effective_length,
        effective_length_mod3=eff_mod3,
    )


# -----------------------------------------------------------------------------
# Summaries and statistics
# -----------------------------------------------------------------------------

def summarize_assigns(assigns: np.ndarray) -> Tuple[Dict[str, int], Dict[int, int], Dict[str, int]]:
    bc = np.bincount(assigns, minlength=6)
    c6 = {STATE_NAMES[i]: int(bc[i]) for i in range(6)}
    c3 = {0: int(bc[0] + bc[3]), 1: int(bc[1] + bc[4]), 2: int(bc[2] + bc[5])}
    cstrand = {"+": int(bc[0] + bc[1] + bc[2]), "-": int(bc[3] + bc[4] + bc[5])}
    return c6, c3, cstrand


def chisq_uniform(counts: Sequence[int], k: int) -> Tuple[float, float]:
    obs = np.array(counts, dtype=np.float64)
    if obs.sum() == 0:
        return np.nan, np.nan
    exp = np.ones(k, dtype=np.float64) * (obs.sum() / k)
    stat = np.sum((obs - exp) ** 2 / exp)
    pval = 1.0 - chi2.cdf(stat, df=k - 1)
    return float(stat), float(pval)


def summary_stats(arr: np.ndarray, prefix: str) -> Dict[str, float]:
    arr = np.asarray(arr, dtype=np.float64)
    if arr.size == 0:
        return {
            f"{prefix}_mean": np.nan,
            f"{prefix}_median": np.nan,
            f"{prefix}_sd": np.nan,
            f"{prefix}_min": np.nan,
            f"{prefix}_q05": np.nan,
            f"{prefix}_q25": np.nan,
            f"{prefix}_q75": np.nan,
            f"{prefix}_q95": np.nan,
            f"{prefix}_max": np.nan,
        }
    q = np.quantile(arr, [0.05, 0.25, 0.75, 0.95])
    return {
        f"{prefix}_mean": float(np.mean(arr)),
        f"{prefix}_median": float(np.median(arr)),
        f"{prefix}_sd": float(np.std(arr, ddof=0)),
        f"{prefix}_min": float(np.min(arr)),
        f"{prefix}_q05": float(q[0]),
        f"{prefix}_q25": float(q[1]),
        f"{prefix}_q75": float(q[2]),
        f"{prefix}_q95": float(q[3]),
        f"{prefix}_max": float(np.max(arr)),
    }


def summarize_assignment_gaps(
    ll_mat: np.ndarray,
    assigns_raw: np.ndarray,
    n_info: np.ndarray,
    raw_length: int,
    *,
    state_map_old_to_new: np.ndarray,
) -> Tuple[Dict[str, object], List[Dict[str, object]], List[Dict[str, object]]]:
    order = np.argsort(ll_mat, axis=1)
    best_idx_raw = order[:, -1].astype(np.int32)
    second_idx_raw = order[:, -2].astype(np.int32)
    best_ll = ll_mat[np.arange(ll_mat.shape[0]), best_idx_raw]
    second_ll = ll_mat[np.arange(ll_mat.shape[0]), second_idx_raw]
    gap = best_ll - second_ll
    gap_per_base = gap / np.maximum(n_info.astype(np.float64), 1.0)

    assigns = state_map_old_to_new[assigns_raw]
    second_idx = state_map_old_to_new[second_idx_raw]

    summary = {
        "length": int(raw_length),
        "n_reads": int(ll_mat.shape[0]),
        "mean_informative_bases": float(np.mean(n_info)) if n_info.size else np.nan,
        "fraction_tie_best_vs_runnerup": float(np.mean(gap <= 1e-12)) if gap.size else np.nan,
        "fraction_gap_lt_0p1": float(np.mean(gap < 0.1)) if gap.size else np.nan,
        "fraction_gap_lt_1": float(np.mean(gap < 1.0)) if gap.size else np.nan,
    }
    summary.update(summary_stats(gap, "gap_ll"))
    summary.update(summary_stats(gap_per_base, "gap_ll_per_base"))

    state_rows: List[Dict[str, object]] = []
    pair_rows: List[Dict[str, object]] = []

    for state_idx, state_name in enumerate(STATE_NAMES):
        mask = assigns == state_idx
        n_state = int(np.sum(mask))
        row = {"length": int(raw_length), "assigned_state": state_name, "n_reads": n_state}
        row.update(summary_stats(gap[mask], "gap_ll"))
        row.update(summary_stats(gap_per_base[mask], "gap_ll_per_base"))
        state_rows.append(row)

        denom = max(n_state, 1)
        for runner_idx, runner_name in enumerate(STATE_NAMES):
            if runner_idx == state_idx:
                continue
            cnt = int(np.sum(second_idx[mask] == runner_idx))
            pair_rows.append(
                {
                    "length": int(raw_length),
                    "assigned_state": state_name,
                    "runnerup_state": runner_name,
                    "count": cnt,
                    "prop_within_assigned": float(cnt / denom) if n_state else np.nan,
                }
            )

    return summary, state_rows, pair_rows


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def plot_ll_trace(ll_trace: Sequence[float], outprefix: str) -> None:
    plt.figure()
    plt.plot(np.arange(1, len(ll_trace) + 1), ll_trace)
    plt.xlabel("Iteration")
    plt.ylabel("Total log-likelihood")
    plt.tight_layout()
    plt.savefig(outprefix + ".png", dpi=200)
    plt.savefig(outprefix + ".pdf")
    plt.close()


def plot_base_probs(p: np.ndarray, outprefix: str, title: Optional[str] = None) -> None:
    df = pd.DataFrame(p, columns=BASES)
    df["pos"] = ["1", "2", "3"]
    x = np.arange(3)
    width = 0.18
    plt.figure()
    for j, base in enumerate(BASES):
        plt.bar(x + (j - 1.5) * width, df[base].values, width=width, label=base)
    plt.xticks(x, df["pos"])
    plt.ylim(0.0, 1.0)
    plt.xlabel("Position class")
    plt.ylabel("Estimated base probability")
    if title:
        plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outprefix + ".png", dpi=200)
    plt.savefig(outprefix + ".pdf")
    plt.close()


def plot_state6_props(df_shift: pd.DataFrame, outdir: str, canonicalize_mode: str) -> None:
    x = df_shift["length"].values
    plt.figure()
    for st in STATE_NAMES:
        plt.plot(x, df_shift[f"prop_{st}"].values, marker="o", label=st)
    plt.axhline(1 / 6, linestyle="--")
    plt.xlabel("Read length")
    plt.ylabel("Proportion")
    plt.title(f"Inferred 6-state distribution by length ({canonicalize_mode})")
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "state6_props_by_length.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "state6_props_by_length.pdf"))
    plt.close()


def plot_shift3_props(df_shift: pd.DataFrame, outdir: str, canonicalize_mode: str) -> None:
    x = df_shift["length"].values
    plt.figure()
    for shift in (0, 1, 2):
        plt.plot(x, df_shift[f"prop_shift{shift}"].values, marker="o", label=f"shift {shift}")
    plt.axhline(1 / 3, linestyle="--")
    plt.xlabel("Read length")
    plt.ylabel("Proportion")
    plt.title(f"Pooled 3-shift distribution by length ({canonicalize_mode})")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "shift3_props_by_length.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "shift3_props_by_length.pdf"))
    plt.close()


def plot_orientation(df_shift: pd.DataFrame, outdir: str, canonicalize_mode: str) -> None:
    x = df_shift["length"].values
    plt.figure()
    plt.plot(x, df_shift["prop_plus"].values, marker="o", label="+ orientation")
    plt.plot(x, df_shift["prop_minus"].values, marker="o", label="- orientation")
    plt.axhline(0.5, linestyle="--")
    plt.xlabel("Read length")
    plt.ylabel("Proportion")
    plt.title(f"Inferred orientation usage ({canonicalize_mode})")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "strand_balance_by_length.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "strand_balance_by_length.pdf"))
    plt.close()


def plot_state6_heatmap(df_shift: pd.DataFrame, outdir: str, canonicalize_mode: str) -> None:
    lengths = df_shift["length"].astype(int).tolist()
    mat = np.vstack([df_shift[f"prop_{st}"].values for st in STATE_NAMES]).T
    plt.figure(figsize=(7, max(2, 0.35 * len(lengths))))
    plt.imshow(mat, aspect="auto")
    plt.yticks(np.arange(len(lengths)), lengths)
    plt.xticks(np.arange(len(STATE_NAMES)), STATE_NAMES)
    plt.colorbar(label="Proportion")
    plt.title(f"6-state assignments by length ({canonicalize_mode})")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "state6_heatmap_by_length.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "state6_heatmap_by_length.pdf"))
    plt.close()


def plot_gap_summary(df_gap: pd.DataFrame, outdir: str) -> None:
    x = df_gap["length"].values
    plt.figure()
    plt.plot(x, df_gap["gap_ll_mean"].values, marker="o", label="mean top1-top2 gap")
    plt.plot(x, df_gap["gap_ll_median"].values, marker="o", label="median top1-top2 gap")
    plt.xlabel("Read length")
    plt.ylabel("Log-likelihood gap")
    plt.title("Hard-EM assignment gap by length")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "assignment_gap_by_length.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "assignment_gap_by_length.pdf"))
    plt.close()


def write_latex_table(df_shift: pd.DataFrame, outdir: str, canonicalize_mode: str) -> None:
    lines = [
        r"\begin{table}[t]",
        r"\centering",
        r"\begin{tabular}{r r r r r r r r}",
        r"\hline",
        r"Length & $n$ & +0 & +1 & +2 & -0 & -1 & -2 \\",
        r"\hline",
    ]
    for _, row in df_shift.iterrows():
        lines.append(
            f"{int(row['length'])} & {int(row['n_sequences'])} & "
            f"{row['prop_+0']:.4f} & {row['prop_+1']:.4f} & {row['prop_+2']:.4f} & "
            f"{row['prop_-0']:.4f} & {row['prop_-1']:.4f} & {row['prop_-2']:.4f} \\\\")
    lines.extend(
        [
            r"\hline",
            r"\end{tabular}",
            rf"\caption{{Inferred 6-state proportions under the hard-EM model, displayed with canonicalization mode `{canonicalize_mode}`.}}",
            r"\label{tab:state6_mle}",
            r"\end{table}",
        ]
    )
    with open(os.path.join(outdir, "latex_tables.tex"), "w") as handle:
        handle.write("\n".join(lines) + "\n")


# -----------------------------------------------------------------------------
# Main driver
# -----------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Unified 6-state hard-EM model for BAM and FASTA/FASTQ inputs.")
    parser.add_argument("--input", nargs="+", required=True, help="Input file(s): BAM, FASTA, FASTQ, FASTA.gz, or FASTQ.gz")
    parser.add_argument("--input-type", choices=["auto", "bam", "fasta", "fastq"], default="auto")
    parser.add_argument("--out", required=True)
    parser.add_argument("--lengths", required=True, help="Comma-separated raw read lengths to analyze, e.g. 59,60,61")
    parser.add_argument("--trim5", type=int, default=0)
    parser.add_argument("--trim3", type=int, default=0)
    parser.add_argument(
        "--dedup",
        choices=["none", "sequence", "seq_rc", "hash64", "hash64_rc"],
        default="none",
        help="De-duplication mode applied after trimming",
    )
    parser.add_argument("--max-unique-per-length", type=int, default=None)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--max-iter", type=int, default=100)
    parser.add_argument("--tol", type=float, default=1e-6)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--restarts", type=int, default=10)
    parser.add_argument("--threads", type=int, default=10)
    parser.add_argument("--progress-every", type=int, default=2000000)
    parser.add_argument("--canonicalize", choices=["none", "rotation", "full6"], default="full6")

    # BAM-specific filters; ignored for FASTA/FASTQ.
    parser.add_argument("--include-secondary", action="store_true")
    parser.add_argument("--include-supplementary", action="store_true")
    parser.add_argument("--include-qcfail", action="store_true")
    parser.add_argument("--include-duplicates-flag", action="store_true")
    parser.add_argument("--allow-ambiguous", action="store_true", help="Retain reads containing non-ACGT symbols; those positions are ignored")
    return parser


def collect_counts_by_length(args: argparse.Namespace, input_kind: str, lengths: Sequence[int]):
    os.makedirs(args.out, exist_ok=True)

    cp_lists = {L: [] for L in lengths}
    cm_lists = {L: [] for L in lengths}
    seen = {L: set() for L in lengths} if args.dedup != "none" else None

    def all_full() -> bool:
        if args.max_unique_per_length is None or args.dedup == "none":
            return False
        return all(len(seen[L]) >= args.max_unique_per_length for L in lengths)

    iterator = iter_sequences(
        args.input,
        input_kind,
        include_secondary=args.include_secondary,
        include_supplementary=args.include_supplementary,
        include_qcfail=args.include_qcfail,
        include_duplicates_flag=args.include_duplicates_flag,
    )

    t0 = time.time()
    scanned = 0
    skipped_ambiguous = 0
    skipped_duplicate = 0
    for seq in iterator:
        scanned += 1
        if args.progress_every and scanned % args.progress_every == 0:
            dt = time.time() - t0
            got = ", ".join([f"L{L}:{len(cp_lists[L]):,}" for L in lengths])
            print(f"[progress] scanned={scanned:,} ({dt/3600:.2f}h) collected: {got}", flush=True)

        raw_length = len(seq)
        if raw_length not in cp_lists:
            continue

        trimmed = trim_seq(seq, args.trim5, args.trim3)
        if not trimmed:
            continue
        if not args.allow_ambiguous and any(base not in B2I for base in trimmed):
            skipped_ambiguous += 1
            continue

        if args.dedup != "none":
            if args.max_unique_per_length is not None and len(seen[raw_length]) >= args.max_unique_per_length:
                if all_full():
                    print("[info] reached max unique for all lengths; stopping scan.", flush=True)
                    break
                continue
            key = dedup_key(trimmed, args.dedup)
            if key in seen[raw_length]:
                skipped_duplicate += 1
                continue
            seen[raw_length].add(key)

        cp_count = counts_mod3(trimmed)
        if int(cp_count.sum()) == 0:
            skipped_ambiguous += 1
            continue
        cp_lists[raw_length].append(cp_count)
        cm_lists[raw_length].append(counts_mod3(revcomp(trimmed)))

        if all_full():
            print("[info] reached max unique for all lengths; stopping scan.", flush=True)
            break

    collection_stats = {
        "records_scanned": scanned,
        "skipped_ambiguous": skipped_ambiguous,
        "skipped_duplicate": skipped_duplicate,
        "collected_by_length": {str(L): len(cp_lists[L]) for L in lengths},
    }
    return cp_lists, cm_lists, collection_stats


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    lengths = sorted({int(x.strip()) for x in args.lengths.split(",") if x.strip()})
    if not lengths or any(length <= 0 for length in lengths):
        raise ValueError("--lengths must contain positive integers.")
    if any(length - args.trim5 - args.trim3 <= 0 for length in lengths):
        raise ValueError("trim5 + trim3 must be smaller than every requested raw length.")
    if args.max_unique_per_length is not None:
        if args.max_unique_per_length <= 0:
            raise ValueError("--max-unique-per-length must be positive.")
        if args.dedup == "none":
            raise ValueError("--max-unique-per-length requires a non-'none' --dedup mode.")

    input_kind = detect_input_kind(args.input, args.input_type)
    os.makedirs(args.out, exist_ok=True)

    if args.trim5 < 0 or args.trim3 < 0:
        raise ValueError("trim5 and trim3 must be non-negative.")
    if args.alpha <= 0 or args.max_iter < 1 or args.tol < 0 or args.restarts < 1 or args.threads < 1:
        raise ValueError("Require alpha>0, max_iter>=1, tol>=0, restarts>=1, and threads>=1.")
    cp_lists, cm_lists, collection_stats = collect_counts_by_length(args, input_kind, lengths)

    rows_probs: List[Dict[str, object]] = []
    rows_shift: List[Dict[str, object]] = []
    rows_gap: List[Dict[str, object]] = []
    rows_gap_state: List[Dict[str, object]] = []
    rows_runnerup: List[Dict[str, object]] = []
    rows_canon: List[Dict[str, object]] = []

    for raw_length in lengths:
        if not cp_lists[raw_length]:
            continue

        cp = np.stack(cp_lists[raw_length], axis=0)
        cm = np.stack(cm_lists[raw_length], axis=0)
        seeds = [args.seed + (raw_length * 10007) + r for r in range(max(1, args.restarts))]

        best_ll, seed_used, p_raw, assigns_raw, ll_trace = run_restarts(
            cp,
            cm,
            alpha=args.alpha,
            max_iter=args.max_iter,
            tol=args.tol,
            seeds=seeds,
            threads=args.threads,
        )

        effective_length = raw_length - args.trim5 - args.trim3
        canon = canonicalize_solution(
            p_raw,
            mode=args.canonicalize,
            effective_length=effective_length,
        )

        assigns_display = canon.state_map_old_to_new[assigns_raw]
        c6, c3, cstrand = summarize_assigns(assigns_display)

        ll_mat = score_all_states(cp, cm, p_raw)
        n_info = cp.sum(axis=(1, 2)).astype(np.int32)
        gap_summary, gap_state_rows, runnerup_rows = summarize_assignment_gaps(
            ll_mat,
            assigns_raw,
            n_info,
            raw_length,
            state_map_old_to_new=canon.state_map_old_to_new,
        )

        plot_ll_trace(ll_trace, os.path.join(args.out, f"ll_trace_length_{raw_length}"))
        title = f"Estimated base probabilities (L={raw_length}; canonicalize={args.canonicalize})"
        plot_base_probs(canon.p_display, os.path.join(args.out, f"base_probs_length_{raw_length}"), title=title)

        for pos_idx in range(3):
            for base in BASES:
                rows_probs.append(
                    {
                        "length": raw_length,
                        "effective_length": effective_length,
                        "position_class": pos_idx + 1,
                        "base": base,
                        "prob": float(canon.p_display[pos_idx, B2I[base]]),
                        "canonicalization_mode": args.canonicalize,
                        "canonical_rotation": canon.rotation,
                        "canonical_rc_flip": canon.rc_flip,
                    }
                )

        n = int(cp.shape[0])
        stat6, p6 = chisq_uniform([c6[s] for s in STATE_NAMES], 6)
        stat3, p3 = chisq_uniform([c3[s] for s in (0, 1, 2)], 3)
        stat_orient, p_orient = chisq_uniform([cstrand[s] for s in ("+", "-")], 2)

        row = {
            "length": raw_length,
            "effective_length": effective_length,
            "n_sequences": n,
            "n_unique_sequences": n if args.dedup != "none" else None,
            "best_LL": float(best_ll),
            "seed_used": int(seed_used),
            "n_em_e_steps": int(len(ll_trace)),
            "em_converged_by_assignment_stability": True,
            "canonicalization_mode": args.canonicalize,
            "canonical_rotation": canon.rotation,
            "canonical_rc_flip": bool(canon.rc_flip),
            "effective_length_mod3": int(canon.effective_length_mod3),
            "chisq_6states": stat6,
            "p_6states": p6,
            "chisq_3shifts": stat3,
            "p_3shifts": p3,
            "chisq_orientation": stat_orient,
            "p_orientation": p_orient,
            "prop_plus": cstrand["+"] / n if n else np.nan,
            "prop_minus": cstrand["-"] / n if n else np.nan,
        }
        row.update({k: v for k, v in gap_summary.items() if k not in {"length", "n_reads"}})
        for st in STATE_NAMES:
            row[f"count_{st}"] = int(c6[st])
            row[f"prop_{st}"] = c6[st] / n if n else np.nan
        for shift in (0, 1, 2):
            row[f"count_shift{shift}"] = int(c3[shift])
            row[f"prop_shift{shift}"] = c3[shift] / n if n else np.nan
        rows_shift.append(row)

        rows_canon.append(
            {
                "length": raw_length,
                "effective_length": effective_length,
                "canonicalization_mode": args.canonicalize,
                "canonical_rotation": canon.rotation,
                "canonical_rc_flip": bool(canon.rc_flip),
                "effective_length_mod3": int(canon.effective_length_mod3),
            }
        )
        rows_gap.append(gap_summary)
        rows_gap_state.extend(gap_state_rows)
        rows_runnerup.extend(runnerup_rows)

    df_probs = pd.DataFrame(rows_probs).sort_values(["length", "position_class", "base"]) if rows_probs else pd.DataFrame()
    df_shift = pd.DataFrame(rows_shift).sort_values("length") if rows_shift else pd.DataFrame()
    df_gap = pd.DataFrame(rows_gap).sort_values("length") if rows_gap else pd.DataFrame()
    df_gap_state = pd.DataFrame(rows_gap_state).sort_values(["length", "assigned_state"]) if rows_gap_state else pd.DataFrame()
    df_runnerup = pd.DataFrame(rows_runnerup).sort_values(["length", "assigned_state", "runnerup_state"]) if rows_runnerup else pd.DataFrame()
    df_canon = pd.DataFrame(rows_canon).sort_values("length") if rows_canon else pd.DataFrame()

    if not df_probs.empty:
        df_probs.to_csv(os.path.join(args.out, "base_probs_by_length.tsv"), sep="\t", index=False)
    if not df_shift.empty:
        df_shift.to_csv(os.path.join(args.out, "shift_distribution_by_length.tsv"), sep="\t", index=False)
        df_shift.to_csv(os.path.join(args.out, "state_distribution_by_length.tsv"), sep="\t", index=False)
    if not df_gap.empty:
        df_gap.to_csv(os.path.join(args.out, "assignment_gap_summary_by_length.tsv"), sep="\t", index=False)
    if not df_gap_state.empty:
        df_gap_state.to_csv(os.path.join(args.out, "assignment_gap_by_state_by_length.tsv"), sep="\t", index=False)
    if not df_runnerup.empty:
        df_runnerup.to_csv(os.path.join(args.out, "assignment_runnerup_by_length.tsv"), sep="\t", index=False)
    if not df_canon.empty:
        df_canon.to_csv(os.path.join(args.out, "canonicalization_by_length.tsv"), sep="\t", index=False)

    if not df_shift.empty:
        plot_state6_props(df_shift, args.out, args.canonicalize)
        plot_shift3_props(df_shift, args.out, args.canonicalize)
        plot_orientation(df_shift, args.out, args.canonicalize)
        plot_state6_heatmap(df_shift, args.out, args.canonicalize)
        write_latex_table(df_shift, args.out, args.canonicalize)
    if not df_gap.empty:
        plot_gap_summary(df_gap, args.out)

    metadata = {
        "script": "six_state_em.py",
        "script_version": SCRIPT_VERSION,
        "script_sha256": hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
        "input": args.input,
        "input_type": input_kind,
        "out": args.out,
        "lengths": lengths,
        "trim5": args.trim5,
        "trim3": args.trim3,
        "dedup": args.dedup,
        "max_unique_per_length": args.max_unique_per_length,
        "alpha": args.alpha,
        "max_iter": args.max_iter,
        "tol": args.tol,
        "seed": args.seed,
        "restarts": args.restarts,
        "threads": args.threads,
        "progress_every": args.progress_every,
        "canonicalize": args.canonicalize,
        "include_secondary": args.include_secondary,
        "include_supplementary": args.include_supplementary,
        "include_qcfail": args.include_qcfail,
        "include_duplicates_flag": args.include_duplicates_flag,
        "bam_sequence_orientation": (
            "original sequenced orientation restored with get_forward_sequence"
            if input_kind == "bam"
            else None
        ),
        "allow_ambiguous": args.allow_ambiguous,
        "collection_stats": collection_stats,
        "initialization": (
            "global base composition plus Normal(0, 0.01), clipped and row-normalized; "
            "pseudocount alpha is applied only in the M-step"
        ),
        "restart_selection": "largest final complete-data hard-assignment log score",
    }
    with open(os.path.join(args.out, "run_metadata.json"), "w") as handle:
        json.dump(metadata, handle, indent=2)

    print(f"[done] Wrote results to: {args.out}")


if __name__ == "__main__":
    main()

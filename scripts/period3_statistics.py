#!/usr/bin/env python3
"""
period3_statistics.py

Statistical utilities for the article:
"Widespread phase-coherent three-base periodicity in metagenomic sequencing reads".

The routines in this module preserve the original statistical rationale:
period-3 structure is measured in 1-mer cycle-fraction CSV files by fitting
sine/cosine components with a fixed period of three cycles and estimating
p-values by permutation. The plotting helpers are intentionally kept visually
consistent with the original script.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from math import atan2, pi, sqrt
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ARTICLE_TITLE = "Widespread phase-coherent three-base periodicity in metagenomic sequencing reads"
SCRIPT_VERSION = "3.0.0-audited"
BASES: Tuple[str, ...] = ("A", "T", "G", "C")

# Stable classic colors retained for compatibility with the plotting modules.
DEFAULT_1MER_COLORS = {
    "A": "#1f77b4",
    "T": "#ff7f0e",
    "G": "#2ca02c",
    "C": "#d62728",
}


@dataclass(frozen=True)
class Period3Fit:
    """Period-3 model fit for one nucleotide base."""

    base: str
    start: int
    end: int
    n_cycles: int
    intercept: float
    cos_coef: float
    sin_coef: float
    amplitude: float
    phase_rad: float
    r2: float
    p_value: float
    p_value_adjusted: float
    significant: bool

    def phase_degrees(self) -> float:
        """Return phase angle in degrees within [0, 360)."""

        return (self.phase_rad * 180.0 / pi) % 360.0

    # Backward-compatible method name used by the previous script.
    def phase_deg(self) -> float:
        return self.phase_degrees()


@dataclass(frozen=True)
class Period3TestResult:
    """Result of testing one 1-mer CSV for period-3 structure."""

    csv_path: str
    alpha: float
    fits: List[Period3Fit]
    window: Tuple[int, int]
    all_significant: bool


@dataclass(frozen=True)
class Period3ComparisonResult:
    """Shape similarity between two period-3 phase profiles."""

    csv1: str
    csv2: str
    window: Tuple[int, int]
    best_shift: int
    similarity_percent: float
    per_base_similarity_percent: Dict[str, float]
    pattern1: Dict[str, List[float]]
    pattern2: Dict[str, List[float]]


# Backward-compatible dataclass aliases.
TestResult = Period3TestResult
CompareResult = Period3ComparisonResult


def load_one_mer_csv(csv_path: str | Path) -> pd.DataFrame:
    """Load and validate a 1-mer read-position CSV.

    Publication-facing files use 'position'. The legacy 'cycle' header is
    still accepted for backward compatibility and normalized internally.
    """

    df = pd.read_csv(csv_path)
    coordinate = (
        "position" if "position" in df.columns
        else "cycle" if "cycle" in df.columns
        else None
    )
    expected_bases = set(BASES)
    if coordinate is None or not expected_bases.issubset(df.columns):
        raise ValueError(
            "CSV must contain a coordinate column named 'position' (preferred) "
            f"or legacy 'cycle', plus {sorted(expected_bases)}; got {list(df.columns)}"
        )
    df = df.loc[:, [coordinate, *BASES]].copy()
    if coordinate != "cycle":
        df = df.rename(columns={coordinate: "cycle"})
    df["cycle"] = pd.to_numeric(df["cycle"], errors="raise").astype(int)
    if df["cycle"].duplicated().any():
        duplicate_positions = sorted(df.loc[df["cycle"].duplicated(), "cycle"].unique())
        raise ValueError(f"Read-position column contains duplicates: {duplicate_positions[:10]}")
    for base in BASES:
        df[base] = pd.to_numeric(df[base], errors="raise").astype(float)
        if not np.isfinite(df[base]).all():
            raise ValueError(f"Column {base} contains non-finite values.")
        if ((df[base] < 0.0) | (df[base] > 1.0)).any():
            raise ValueError(f"Column {base} contains values outside [0, 1].")
    return df.sort_values("cycle").reset_index(drop=True)

def slice_cycle_window(df: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """Return rows whose cycle is in the inclusive [start, end] window."""

    if start > end:
        raise ValueError(f"Window start ({start}) must not exceed end ({end}).")
    subset = df[(df["cycle"] >= start) & (df["cycle"] <= end)].copy()
    if subset.empty:
        raise ValueError(
            f"No cycles in requested window [{start}, {end}]. "
            f"Available: {df['cycle'].min()}..{df['cycle'].max()}."
        )
    return subset.sort_values("cycle").reset_index(drop=True)


def build_period3_design(n_cycles: int, offset: int = 0) -> np.ndarray:
    """Build cos/sin design columns for a period of three cycles."""

    if n_cycles <= 0:
        raise ValueError("n_cycles must be positive.")
    t = np.arange(n_cycles, dtype=float) + offset
    omega = 2.0 * np.pi / 3.0
    return np.column_stack([np.cos(omega * t), np.sin(omega * t)])


def _resolve_cycle_indices(
    y: np.ndarray,
    abs_start_cycle: Optional[int],
    cycle_indices: Optional[np.ndarray],
) -> np.ndarray:
    if cycle_indices is None:
        if abs_start_cycle is None:
            raise ValueError("Provide abs_start_cycle or cycle_indices.")
        return np.arange(abs_start_cycle, abs_start_cycle + len(y), dtype=float)
    cycles = np.asarray(cycle_indices, dtype=float)
    if cycles.ndim != 1 or len(cycles) != len(y):
        raise ValueError("cycle_indices must be one-dimensional and match y.")
    if not np.isfinite(cycles).all():
        raise ValueError("cycle_indices contains non-finite values.")
    return cycles


def fit_period3_component(
    y: np.ndarray,
    abs_start_cycle: Optional[int] = None,
    *,
    cycle_indices: Optional[np.ndarray] = None,
) -> Tuple[float, float, float, float]:
    """Regress values on intercept + cos(2*pi*k/3) + sin(2*pi*k/3)."""

    values = np.asarray(y, dtype=float)
    if values.ndim != 1 or values.size < 3:
        raise ValueError("y must be one-dimensional with at least three values.")
    if not np.isfinite(values).all():
        raise ValueError("y contains non-finite values.")
    abs_idx = _resolve_cycle_indices(values, abs_start_cycle, cycle_indices)
    omega = 2.0 * np.pi / 3.0
    design = np.column_stack(
        [np.ones(values.size), np.cos(omega * abs_idx), np.sin(omega * abs_idx)]
    )
    beta = np.linalg.lstsq(design, values, rcond=None)[0]
    y_hat = design @ beta
    ss_total = float(((values - values.mean()) ** 2).sum())
    ss_error = float(((values - y_hat) ** 2).sum())
    r2 = 0.0 if ss_total == 0.0 else 1.0 - ss_error / ss_total
    return float(beta[0]), float(beta[1]), float(beta[2]), float(np.clip(r2, 0.0, 1.0))


def permutation_p_value_for_r2(
    y: np.ndarray,
    abs_start_cycle: Optional[int],
    observed_r2: float,
    n_permutations: int = 10_000,
    seed: Optional[int] = None,
    *,
    cycle_indices: Optional[np.ndarray] = None,
) -> float:
    """Estimate a permutation p-value using R^2 as the test statistic."""

    if n_permutations < 1:
        raise ValueError("n_permutations must be at least 1.")
    values = np.asarray(y, dtype=float)
    cycles = _resolve_cycle_indices(values, abs_start_cycle, cycle_indices)
    rng = np.random.default_rng(seed)
    greater_or_equal = 0
    for _ in range(n_permutations):
        permuted = rng.permutation(values)
        _, _, _, r2_perm = fit_period3_component(
            permuted, cycle_indices=cycles
        )
        if r2_perm >= observed_r2 - 1e-15:
            greater_or_equal += 1
    return float((greater_or_equal + 1.0) / (n_permutations + 1.0))


def holm_adjust(p_values: Sequence[float]) -> List[float]:
    """Holm-adjust a family of p-values while preserving input order."""

    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1 or values.size == 0:
        return []
    if ((values < 0.0) | (values > 1.0) | ~np.isfinite(values)).any():
        raise ValueError("p-values must be finite and lie in [0, 1].")
    order = np.argsort(values)
    adjusted_sorted = np.empty(values.size, dtype=float)
    running = 0.0
    m = values.size
    for rank, idx in enumerate(order):
        candidate = (m - rank) * values[idx]
        running = max(running, candidate)
        adjusted_sorted[rank] = min(1.0, running)
    adjusted = np.empty(values.size, dtype=float)
    for rank, idx in enumerate(order):
        adjusted[idx] = adjusted_sorted[rank]
    return adjusted.tolist()


def test_period3_from_csv(
    csv_path: str | Path,
    start: int = 10,
    end: int = 40,
    alpha: float = 0.05,
    permutations: int = 10_000,
    seed: Optional[int] = 12345,
) -> Period3TestResult:
    """Test one 1-mer CSV for a period-3 component in each nucleotide."""

    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must lie strictly between 0 and 1.")
    df = load_one_mer_csv(csv_path)
    subset = slice_cycle_window(df, start, end)
    n_cycles = len(subset)
    if n_cycles < 9:
        raise ValueError(
            f"Too few cycles in window [{start},{end}] (n={n_cycles}). "
            "Need at least 9 for a stable test."
        )

    cycles = subset["cycle"].to_numpy(dtype=float)
    raw_rows = []
    seed_sequence = np.random.SeedSequence(seed)
    child_seeds = seed_sequence.spawn(len(BASES))
    for base, child_seed in zip(BASES, child_seeds):
        y = subset[base].to_numpy(dtype=float)
        intercept, cos_coef, sin_coef, r2 = fit_period3_component(
            y, cycle_indices=cycles
        )
        amplitude = sqrt(cos_coef**2 + sin_coef**2)
        phase = atan2(sin_coef, cos_coef)
        p_value = permutation_p_value_for_r2(
            y,
            abs_start_cycle=None,
            cycle_indices=cycles,
            observed_r2=r2,
            n_permutations=permutations,
            seed=int(child_seed.generate_state(1, dtype=np.uint64)[0]),
        )
        raw_rows.append((base, intercept, cos_coef, sin_coef, amplitude, phase, r2, p_value))

    adjusted = holm_adjust([row[-1] for row in raw_rows])
    fits: List[Period3Fit] = []
    for row, p_adjusted in zip(raw_rows, adjusted):
        base, intercept, cos_coef, sin_coef, amplitude, phase, r2, p_value = row
        fits.append(
            Period3Fit(
                base=base,
                start=int(cycles[0]),
                end=int(cycles[-1]),
                n_cycles=n_cycles,
                intercept=float(intercept),
                cos_coef=float(cos_coef),
                sin_coef=float(sin_coef),
                amplitude=float(amplitude),
                phase_rad=float(phase),
                r2=float(r2),
                p_value=float(p_value),
                p_value_adjusted=float(p_adjusted),
                significant=bool(p_adjusted < alpha),
            )
        )

    return Period3TestResult(
        csv_path=str(csv_path),
        alpha=float(alpha),
        fits=fits,
        window=(start, end),
        all_significant=all(fit.significant for fit in fits),
    )


# Backward-compatible function name used by the previous script.
def test_threefold_from_csv(*args, **kwargs) -> Period3TestResult:
    return test_period3_from_csv(*args, **kwargs)


def phase_means_by_base(df: pd.DataFrame, start: int, end: int) -> Dict[str, List[float]]:
    """Compute mean nucleotide fraction for each read-position modulo-3 phase."""

    subset = slice_cycle_window(df, start, end)
    phases = (subset["cycle"] % 3).map({0: 0, 1: 1, 2: 2}).to_numpy()
    pattern: Dict[str, List[float]] = {}
    for base in BASES:
        values = subset[base].to_numpy(dtype=float)
        pattern[base] = [float(values[phases == phase].mean()) for phase in (0, 1, 2)]
    return pattern


def concatenate_phase_pattern(
    pattern: Dict[str, List[float]], order: Tuple[str, ...] = BASES
) -> np.ndarray:
    """Concatenate base-by-phase means into one 12-dimensional vector."""

    return np.concatenate([np.asarray(pattern[base], dtype=float) for base in order], axis=0)


def zscore_phase_pattern_by_base(pattern: Dict[str, List[float]]) -> np.ndarray:
    """Standardize each base's 3-phase vector to isolate phase shape."""

    parts = []
    for base in BASES:
        vector = np.asarray(pattern[base], dtype=float)
        vector = vector - vector.mean()
        norm = np.linalg.norm(vector)
        if norm > 0:
            vector = vector / norm
        else:
            vector = np.zeros_like(vector)
        parts.append(vector)
    return np.concatenate(parts, axis=0)


def shift_phase_pattern(vector: np.ndarray, shift: int) -> np.ndarray:
    """Apply a global 0/1/2 phase shift within each nucleotide block."""

    shifted = np.empty_like(vector)
    for i in range(4):
        block = vector[i * 3 : (i + 1) * 3]
        shifted[i * 3 : (i + 1) * 3] = np.roll(block, shift)
    return shifted


def compare_period3_patterns(
    csv1: str | Path,
    csv2: str | Path,
    start: int = 10,
    end: int = 40,
) -> Period3ComparisonResult:
    """Compare phase-shape similarity between two 1-mer CSVs."""

    pattern1 = phase_means_by_base(load_one_mer_csv(csv1), start, end)
    pattern2 = phase_means_by_base(load_one_mer_csv(csv2), start, end)
    z1 = zscore_phase_pattern_by_base(pattern1)
    z2 = zscore_phase_pattern_by_base(pattern2)

    best_shift = 0
    best_r = -1.0
    for shift in (0, 1, 2):
        z2_shifted = shift_phase_pattern(z2, shift)
        if np.linalg.norm(z1) == 0 or np.linalg.norm(z2_shifted) == 0:
            r = 0.0
        else:
            r = float(np.corrcoef(z1, z2_shifted)[0, 1])
            if not np.isfinite(r):
                r = 0.0
        if r > best_r:
            best_r = r
            best_shift = shift

    similarity_percent = 100.0 * float(max(0.0, best_r) ** 2)
    per_base: Dict[str, float] = {}
    for i, base in enumerate(BASES):
        v1 = z1[i * 3 : (i + 1) * 3]
        v2 = np.roll(z2[i * 3 : (i + 1) * 3], best_shift)
        if np.linalg.norm(v1) == 0 or np.linalg.norm(v2) == 0:
            r = 0.0
        else:
            r = float(np.corrcoef(v1, v2)[0, 1])
            if not np.isfinite(r):
                r = 0.0
        per_base[base] = 100.0 * float(max(0.0, r) ** 2)

    return Period3ComparisonResult(
        csv1=str(csv1),
        csv2=str(csv2),
        window=(start, end),
        best_shift=int(best_shift),
        similarity_percent=float(similarity_percent),
        per_base_similarity_percent=per_base,
        pattern1=pattern1,
        pattern2=pattern2,
    )


# Backward-compatible function name used by the previous script.
def compare_threefold_patterns(*args, **kwargs) -> Period3ComparisonResult:
    return compare_period3_patterns(*args, **kwargs)


def ensure_output_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def plot_phase_profile_single(
    pattern: Dict[str, List[float]],
    window: Tuple[int, int],
    title: str,
    save_to: Path,
) -> str:
    """Save grouped bars of phase means for A/T/G/C."""

    ensure_output_directory(save_to.parent)
    fig, ax = plt.subplots(figsize=(8.5, 4.5), dpi=300)
    phases = [0, 1, 2]
    x_positions = np.arange(len(BASES))
    width = 0.22
    for j, phase in enumerate(phases):
        values = [pattern[base][phase] for base in BASES]
        ax.bar(x_positions + (j - 1) * width, values, width=width, label=f"phase {phase}")
    ax.set_xticks(x_positions, labels=BASES)
    ax.set_ylabel("Mean fraction (positions mod 3)")
    ax.set_title(f"{title}\nphase profile {window[0]}–{window[1]}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_to, bbox_inches="tight")
    plt.close(fig)
    return str(save_to)


def plot_strength_bars(
    fits: List[Period3Fit],
    window: Tuple[int, int],
    title: str,
    save_to: Path,
) -> str:
    """Save bar chart of period-3 R^2 values for A/T/G/C."""

    ensure_output_directory(save_to.parent)
    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=300)
    x_positions = np.arange(len(BASES))
    values = [next(fit.r2 for fit in fits if fit.base == base) for base in BASES]
    ax.bar(x_positions, values)
    ax.set_xticks(x_positions, labels=BASES)
    ax.set_ylabel("R² (period‑3 component)")
    ax.set_ylim(0, max(0.05, max(values) * 1.15))
    ax.set_title(f"{title}\nstrength (R²) {window[0]}–{window[1]}")
    fig.tight_layout()
    fig.savefig(save_to, bbox_inches="tight")
    plt.close(fig)
    return str(save_to)


def plot_compare_phase_profiles(
    pattern1: Dict[str, List[float]],
    pattern2: Dict[str, List[float]],
    shift: int,
    window: Tuple[int, int],
    labels: Tuple[str, str],
    save_to: Path,
) -> str:
    """Save side-by-side phase-profile comparison after phase alignment."""

    ensure_output_directory(save_to.parent)
    fig, ax = plt.subplots(figsize=(10.0, 4.8), dpi=320)
    x_labels = []
    values1 = []
    values2 = []
    for base in BASES:
        for phase in (0, 1, 2):
            x_labels.append(f"{base}·{phase}")
            values1.append(pattern1[base][phase])
            values2.append(pattern2[base][(phase - shift) % 3])
    x_positions = np.arange(len(x_labels))
    width = 0.40
    ax.bar(x_positions - width / 2, values1, width=width, label=labels[0])
    ax.bar(x_positions + width / 2, values2, width=width, label=f"{labels[1]} (shift {shift})")
    ax.set_xticks(x_positions, labels=x_labels, rotation=45, ha="right")
    ax.set_ylabel("Mean fraction (positions mod 3)")
    ax.set_title(f"Phase profiles {window[0]}–{window[1]} (shift‑aligned)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_to, bbox_inches="tight")
    plt.close(fig)
    return str(save_to)


def plot_similarity_bars(
    per_base_similarity_percent: Dict[str, float],
    overall_similarity_percent: float,
    save_to: Path,
) -> str:
    """Save per-base similarity bars for a pairwise comparison."""

    ensure_output_directory(save_to.parent)
    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=300)
    x_positions = np.arange(len(BASES))
    values = [per_base_similarity_percent[base] for base in BASES]
    ax.bar(x_positions, values)
    ax.set_xticks(x_positions, labels=BASES)
    ax.set_ylabel("Similarity (% r²)")
    ax.set_ylim(0, 100)
    ax.set_title(f"Per-base similarity (overall: {overall_similarity_percent:.1f}%)")
    fig.tight_layout()
    fig.savefig(save_to, bbox_inches="tight")
    plt.close(fig)
    return str(save_to)


# Backward-compatible private names used by the previous pipeline.
_load_1mer_csv = load_one_mer_csv
_window_df = slice_cycle_window
_design_period3 = build_period3_design
_fit_period3 = fit_period3_component
_perm_p_value_r2 = permutation_p_value_for_r2
_phase_means_by_base = phase_means_by_base
_concat_pattern = concatenate_phase_pattern
_zscore_by_base = zscore_phase_pattern_by_base
_shift_pattern_three = shift_phase_pattern
_plot_phase_profile_single = plot_phase_profile_single
_plot_strength_bars = plot_strength_bars
_plot_compare_phase_profiles = plot_compare_phase_profiles
_plot_similarity_bars = plot_similarity_bars


def _fit_to_json(fit: Period3Fit) -> Dict[str, object]:
    data = asdict(fit)
    data["phase_degrees"] = fit.phase_degrees()
    data.pop("phase_rad", None)
    return data


def _cmd_test_one(ns: argparse.Namespace) -> None:
    result = test_period3_from_csv(
        csv_path=ns.csv,
        start=ns.start,
        end=ns.end,
        alpha=ns.alpha,
        permutations=ns.perms,
        seed=ns.seed,
    )
    output: Dict[str, object] = {
        "article": ARTICLE_TITLE,
        "script_version": SCRIPT_VERSION,
        "csv": result.csv_path,
        "read_length": ns.read_length,
        "read_count": ns.read_count,
        "window": result.window,
        "alpha": result.alpha,
        "all_significant": result.all_significant,
        "fits": [_fit_to_json(fit) for fit in result.fits],
    }

    if ns.plot:
        outdir = Path(ns.plot_dir or Path(ns.csv).parent)
        df = load_one_mer_csv(ns.csv)
        pattern = phase_means_by_base(df, ns.start, ns.end)
        output["plot_phase_profile"] = plot_phase_profile_single(
            pattern=pattern,
            window=(ns.start, ns.end),
            title=Path(ns.csv).name,
            save_to=outdir / f"{Path(ns.csv).stem}_period3_phase_profile.{ns.plot_format}",
        )
        output["plot_strength"] = plot_strength_bars(
            fits=result.fits,
            window=(ns.start, ns.end),
            title=Path(ns.csv).name,
            save_to=outdir / f"{Path(ns.csv).stem}_period3_strength.{ns.plot_format}",
        )

    print(json.dumps(output, indent=2))


def _cmd_compare_two(ns: argparse.Namespace) -> None:
    result = compare_period3_patterns(
        csv1=ns.csv1,
        csv2=ns.csv2,
        start=ns.start,
        end=ns.end,
    )
    output: Dict[str, object] = {
        "article": ARTICLE_TITLE,
        "script_version": SCRIPT_VERSION,
        "csv1": result.csv1,
        "csv2": result.csv2,
        "window": result.window,
        "best_shift": result.best_shift,
        "similarity_percent": result.similarity_percent,
        "per_base_similarity_percent": result.per_base_similarity_percent,
    }

    if ns.plot:
        outdir = Path(ns.plot_dir or Path(ns.csv1).parent)
        output["plot_compare_phase_profiles"] = plot_compare_phase_profiles(
            pattern1=result.pattern1,
            pattern2=result.pattern2,
            shift=result.best_shift,
            window=(ns.start, ns.end),
            labels=(Path(ns.csv1).name, Path(ns.csv2).name),
            save_to=outdir / f"compare_{Path(ns.csv1).stem}_vs_{Path(ns.csv2).stem}_phase_profiles.{ns.plot_format}",
        )
        output["plot_similarity_bars"] = plot_similarity_bars(
            per_base_similarity_percent=result.per_base_similarity_percent,
            overall_similarity_percent=result.similarity_percent,
            save_to=outdir / f"compare_{Path(ns.csv1).stem}_vs_{Path(ns.csv2).stem}_similarity.{ns.plot_format}",
        )

    print(json.dumps(output, indent=2))


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Test and compare period-3 structure in 1-mer read-position fraction CSVs."
    )
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    test_parser = subparsers.add_parser("test", help="Significance test of period-3 structure")
    test_sub = test_parser.add_subparsers(dest="which", required=True)
    one_parser = test_sub.add_parser("one", help="Test one 1-mer CSV")
    one_parser.add_argument("csv", type=str, help="Path to 1-mer CSV with columns cycle,A,T,G,C (cycle is the machine-readable read-position field)")
    one_parser.add_argument("--start", type=int, default=10)
    one_parser.add_argument("--end", type=int, default=40)
    one_parser.add_argument("--alpha", type=float, default=0.05)
    one_parser.add_argument("--perms", type=int, default=10_000, help="Number of permutations")
    one_parser.add_argument("--seed", type=int, default=12345, help="Random seed")
    one_parser.add_argument("--read-length", type=int, default=None, help="Optional read length metadata")
    one_parser.add_argument("--read-count", type=int, default=None, help="Optional read count metadata")
    one_parser.add_argument("--plot", action="store_true", help="Save diagnostic plots")
    one_parser.add_argument("--plot-dir", type=str, default=None, help="Optional plot output directory")
    one_parser.add_argument("--plot-format", choices=("png", "pdf"), default="png")
    one_parser.set_defaults(func=_cmd_test_one)

    compare_parser = subparsers.add_parser("compare", help="Compare two period-3 patterns")
    compare_sub = compare_parser.add_subparsers(dest="which", required=True)
    two_parser = compare_sub.add_parser("two", help="Compare two 1-mer CSV files")
    two_parser.add_argument("csv1", type=str)
    two_parser.add_argument("csv2", type=str)
    two_parser.add_argument("--start", type=int, default=10)
    two_parser.add_argument("--end", type=int, default=40)
    two_parser.add_argument("--plot", action="store_true", help="Save comparison plots")
    two_parser.add_argument("--plot-dir", type=str, default=None)
    two_parser.add_argument("--plot-format", choices=("png", "pdf"), default="png")
    two_parser.set_defaults(func=_cmd_compare_two)

    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

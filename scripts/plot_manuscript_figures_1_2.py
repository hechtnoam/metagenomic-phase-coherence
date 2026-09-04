#!/usr/bin/env python3
"""Generate manuscript Figures 1 and 2 from frozen nucleotide-frequency source data.

Figure 1 shows the A2424 L=59 read-position profile, with the primary
period-3 regression window (positions 10--40 inclusive) shaded in gray.

Figure 2 compares two read-length pairs from different modulo-3 classes:
59/80 and 61/79. Panel A shows nucleotide
fractions by read position; panel B shows mean nucleotide fractions after
positions 10--40 are grouped by position modulo three.

The plotting code intentionally reads the read count from the matching period-3
JSON rather than hard-coding it, preventing figure/caption count drift.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASES: Tuple[str, ...] = ("A", "T", "G", "C")
COLORS = {"A": "tab:blue", "T": "tab:orange", "G": "tab:green", "C": "tab:red"}
DEFAULT_ANALYSIS_START = 10
DEFAULT_ANALYSIS_END = 40

ROOT = Path(__file__).resolve().parents[1]
SOURCE1 = ROOT / "source_data" / "figure_1"
SOURCE2 = ROOT / "source_data" / "figure_2"


@dataclass(frozen=True)
class Profile:
    length: int
    read_count: int
    observed_read_count: int | None
    fractions: pd.DataFrame


def _default_csv(length: int, figure: int) -> Path:
    source = SOURCE1 if figure == 1 else SOURCE2
    return source / f"A2424_L{length}_position_frequencies.csv"


def _default_json(length: int, figure: int) -> Path:
    source = SOURCE1 if figure == 1 else SOURCE2
    return source / f"A2424_L{length}_period3.json"


def load_profile(csv_path: Path, json_path: Path, expected_length: int) -> Profile:
    """Load and cross-check one frozen read-position profile."""
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Missing source-data CSV for L={expected_length}: {csv_path}. "
            "Regenerate that exact-length A2424 profile from the raw run before plotting."
        )
    if not json_path.exists():
        raise FileNotFoundError(
            f"Missing period-3 metadata for L={expected_length}: {json_path}. "
            "Regenerate that exact-length A2424 profile from the raw run before plotting."
        )

    df = pd.read_csv(csv_path)
    required = {"position", *BASES}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{csv_path}: missing required columns {sorted(missing)}")
    df = df.sort_values("position").reset_index(drop=True)

    metadata = json.loads(json_path.read_text())
    actual_length = int(metadata.get("read_length", expected_length))
    if actual_length != expected_length:
        raise ValueError(
            f"{json_path}: read_length={actual_length}, expected {expected_length}"
        )
    read_count = metadata.get("read_count")
    if read_count is None:
        raise ValueError(f"{json_path}: missing read_count")
    observed = metadata.get("observed_read_count")

    if int(df["position"].min()) != 1 or int(df["position"].max()) != expected_length:
        raise ValueError(
            f"{csv_path}: read-position range is {int(df['position'].min())}--"
            f"{int(df['position'].max())}, expected 1--{expected_length}"
        )

    return Profile(
        length=expected_length,
        read_count=int(read_count),
        observed_read_count=int(observed) if observed is not None else None,
        fractions=df,
    )


def shade_analysis_window(ax, start: int, end: int) -> None:
    """Shade the same inclusive analysis window used in Figure 3."""
    ax.axvspan(start - 0.5, end + 0.5, color="0.5", alpha=0.07, zorder=0)


def draw_read_position_profile(
    ax,
    profile: Profile,
    *,
    start: int,
    end: int,
    show_legend: bool = False,
    title_with_count: bool = True,
) -> None:
    shade_analysis_window(ax, start, end)
    for base in BASES:
        ax.scatter(
            profile.fractions["position"],
            profile.fractions[base],
            s=18,
            color=COLORS[base],
            edgecolors="none",
            label=base,
            zorder=2,
        )
    if title_with_count:
        ax.set_title(f"L={profile.length} (n={profile.read_count:,})", fontsize=10)
    else:
        ax.set_title(f"L={profile.length}", fontsize=10)
    ax.set_xlabel("Read position")
    ax.set_ylabel("Fraction")
    ax.set_ylim(0.07, 0.52)
    ax.set_xlim(0, profile.length + 1)
    if show_legend:
        ax.legend(title="Base", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)


def phase_means(profile: Profile, start: int, end: int) -> Mapping[str, Sequence[float]]:
    """Mean base fractions for positions in each modulo-3 class."""
    df = profile.fractions
    subset = df[(df["position"] >= start) & (df["position"] <= end)].copy()
    if subset.empty:
        raise ValueError(f"No read positions in requested window {start}--{end}")
    modulo = subset["position"].to_numpy(dtype=int) % 3
    result: Dict[str, Sequence[float]] = {}
    for base in BASES:
        values = subset[base].to_numpy(dtype=float)
        result[base] = [float(values[modulo == phase].mean()) for phase in (0, 1, 2)]
    return result


def draw_phase_profile(ax, profile: Profile, start: int, end: int) -> None:
    pattern = phase_means(profile, start, end)
    x = np.arange(len(BASES))
    width = 0.22
    for j, phase in enumerate((0, 1, 2)):
        ax.bar(
            x + (j - 1) * width,
            [pattern[base][phase] for base in BASES],
            width=width,
            label=f"phase {phase}",
        )
    ax.set_xticks(x, labels=BASES)
    ax.set_ylabel("Mean fraction")
    ax.set_title(f"L={profile.length}", fontsize=10)
    ax.set_ylim(0, 0.33)


def save_figure(fig, output_prefix: Path) -> None:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_prefix.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Wrote {output_prefix.with_suffix('.png')}")
    print(f"Wrote {output_prefix.with_suffix('.pdf')}")


def make_figure1(
    profile: Profile,
    output_prefix: Path,
    start: int,
    end: int,
) -> None:
    fig, ax = plt.subplots(figsize=(10, 5.4))
    draw_read_position_profile(
        ax,
        profile,
        start=start,
        end=end,
        show_legend=True,
        title_with_count=True,
    )
    ax.set_title(f"A2424 (n={profile.read_count:,})", fontsize=13)
    fig.tight_layout()
    save_figure(fig, output_prefix)
    plt.close(fig)


def _panel_letter(fig, text: str, x: float, y: float) -> None:
    fig.text(x, y, text, fontsize=16, fontweight="bold", va="top", ha="left")


def make_figure2(
    profiles: Mapping[int, Profile],
    output_prefix: Path,
    start: int,
    end: int,
) -> None:
    lengths = (59, 80, 61, 79)
    missing = [length for length in lengths if length not in profiles]
    if missing:
        raise ValueError(f"Figure 2 requires all four lengths; missing {missing}")

    # Build panels A and B as two separate 2x2 blocks.  Keeping the blocks in
    # separate nested grids leaves a dedicated header band above each panel,
    # so the panel title and shared legend cannot overlap one another or the
    # x-axis labels of the preceding block.
    fig = plt.figure(figsize=(14, 16.0))
    outer = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=(1, 1),
        hspace=0.26,
        top=0.90,
        bottom=0.055,
        left=0.09,
        right=0.98,
    )
    gs_a = outer[0].subgridspec(nrows=2, ncols=2, hspace=0.55, wspace=0.28)
    gs_b = outer[1].subgridspec(nrows=2, ncols=2, hspace=0.55, wspace=0.28)

    axes_a = [fig.add_subplot(gs_a[i // 2, i % 2]) for i in range(4)]
    axes_b = [fig.add_subplot(gs_b[i // 2, i % 2]) for i in range(4)]

    for ax, length in zip(axes_a, lengths):
        draw_read_position_profile(
            ax,
            profiles[length],
            start=start,
            end=end,
            show_legend=False,
            title_with_count=True,
        )

    for ax, length in zip(axes_b, lengths):
        draw_phase_profile(ax, profiles[length], start, end)

    # Shared legends avoid repeating the same legend eight times.
    base_handles = [
        plt.Line2D(
            [],
            [],
            linestyle="none",
            marker="o",
            markersize=5,
            color=COLORS[base],
            label=base,
        )
        for base in BASES
    ]
    phase_handles, phase_labels = axes_b[0].get_legend_handles_labels()

    # Draw once so the axes positions are final, then position the two panel
    # headers relative to the actual top of each 2x2 block.  This makes the
    # layout robust to later changes in figure size or subplot spacing.
    fig.canvas.draw()
    top_a = max(ax.get_position().y1 for ax in axes_a)
    top_b = max(ax.get_position().y1 for ax in axes_b)

    header_fontsize = 13
    panel_letter_fontsize = 16

    title_a_y = top_a + 0.060
    legend_a_y = top_a + 0.043
    title_b_y = top_b + 0.060
    legend_b_y = top_b + 0.038

    fig.text(
        0.5,
        title_a_y,
        "Nucleotide fractions by read position",
        ha="center",
        va="bottom",
        fontsize=header_fontsize,
    )
    fig.text(
        0.5,
        title_b_y,
        "Mean base fractions by position modulo three",
        ha="center",
        va="bottom",
        fontsize=header_fontsize,
    )

    fig.legend(
        handles=base_handles,
        labels=BASES,
        title="Base",
        loc="upper center",
        bbox_to_anchor=(0.5, legend_a_y),
        ncol=4,
        frameon=False,
    )
    fig.legend(
        handles=phase_handles,
        labels=phase_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, legend_b_y),
        ncol=3,
        frameon=False,
    )

    fig.text(
        0.035,
        title_a_y + 0.008,
        "A",
        fontsize=panel_letter_fontsize,
        fontweight="bold",
        va="bottom",
        ha="left",
    )
    fig.text(
        0.035,
        title_b_y + 0.008,
        "B",
        fontsize=panel_letter_fontsize,
        fontweight="bold",
        va="bottom",
        ha="left",
    )

    # Label the modulo-3 equivalence class for each row, using the actual row
    # centers rather than hard-coded vertical coordinates.
    row_centers = [
        np.mean([axes_a[0].get_position().y0, axes_a[0].get_position().y1]),
        np.mean([axes_a[2].get_position().y0, axes_a[2].get_position().y1]),
        np.mean([axes_b[0].get_position().y0, axes_b[0].get_position().y1]),
        np.mean([axes_b[2].get_position().y0, axes_b[2].get_position().y1]),
    ]

    save_figure(fig, output_prefix)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", choices=("1", "2", "both"), default="both")
    parser.add_argument("--analysis-start", type=int, default=DEFAULT_ANALYSIS_START)
    parser.add_argument("--analysis-end", type=int, default=DEFAULT_ANALYSIS_END)
    parser.add_argument(
        "--figure1-output-prefix",
        type=Path,
        default=ROOT / "figures" / "fig_a2424_introduction",
    )
    parser.add_argument(
        "--figure2-output-prefix",
        type=Path,
        default=ROOT / "figures" / "fig_read-position_phase_profile_a2424",
    )
    for length in (59, 61, 79, 80):
        parser.add_argument(
            f"--l{length}-csv",
            type=Path,
            default=_default_csv(length, 1 if length == 59 else 2),
        )
        parser.add_argument(
            f"--l{length}-json",
            type=Path,
            default=_default_json(length, 1 if length == 59 else 2),
        )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.analysis_start < 1 or args.analysis_end < args.analysis_start:
        raise ValueError("Require 1 <= --analysis-start <= --analysis-end")

    profile59_fig1 = load_profile(args.l59_csv, args.l59_json, 59)
    if args.figure in ("1", "both"):
        make_figure1(
            profile59_fig1,
            args.figure1_output_prefix,
            args.analysis_start,
            args.analysis_end,
        )

    if args.figure in ("2", "both"):
        profiles = {
            59: load_profile(args.l59_csv, args.l59_json, 59),
            61: load_profile(args.l61_csv, args.l61_json, 61),
            79: load_profile(args.l79_csv, args.l79_json, 79),
            80: load_profile(args.l80_csv, args.l80_json, 80),
        }
        make_figure2(
            profiles,
            args.figure2_output_prefix,
            args.analysis_start,
            args.analysis_end,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

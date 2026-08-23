#!/usr/bin/env python3
"""Create publication phase-redistribution dot figure."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

BASES = ("A", "T", "G", "C")
COLORS = {"A": "tab:blue", "T": "tab:orange", "G": "tab:green", "C": "tab:red"}

DEFAULTS = {
    "original": Path("results/libraries/A2424/csv/kmer/L59/k1.csv"),
    "randomized": Path("results/phase_redistribution/analysis/seed20/A2424.phase_randomized/csv/kmer/L57/k1.csv"),
    "fixed02": Path("results/phase_redistribution/analysis/A2424.fixed_trim5_0/csv/kmer/L57/k1.csv"),
    "fixed11": Path("results/phase_redistribution/analysis/A2424.fixed_trim5_1/csv/kmer/L57/k1.csv"),
    "fixed20": Path("results/phase_redistribution/analysis/A2424.fixed_trim5_2/csv/kmer/L57/k1.csv"),
}

PANELS = [
    ("original", "A", "Original, L=59", 306862),
    ("randomized", "B", "Phase randomized, L=57", 306853),
    ("fixed02", "C", "Fixed trim 0/2, L=57", 306857),
    ("fixed11", "D", "Fixed trim 1/1, L=57", 306853),
    ("fixed20", "E", "Fixed trim 2/0, L=57", 306854),
]

def load(path):
    df = pd.read_csv(path)
    required = {"cycle", *BASES}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    return df

def draw(ax, df, letter, title, n, ymin, ymax, start, end):
    ax.axvspan(start - 0.5, end + 0.5, color="0.5", alpha=0.07, zorder=0)
    for base in BASES:
        ax.scatter(df["cycle"], df[base], s=15, color=COLORS[base],
                   edgecolors="none", zorder=2)
    ax.set_title(f"{title} (n={n:,})", fontsize=10, pad=7)
    ax.set_xlabel("Position", fontsize=9)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0, max(60, int(df["cycle"].max()) + 1))
    ax.tick_params(labelsize=8)
    ax.text(-0.10, 1.12, letter, transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="top", ha="left")

def main():
    ap = argparse.ArgumentParser()
    for key, path in DEFAULTS.items():
        ap.add_argument(f"--{key}", type=Path, default=path)
    ap.add_argument("--output-prefix", type=Path,
                    default=Path("figures/fig_phase_redistribution_final"))
    ap.add_argument("--y-min", type=float, default=0.07)
    ap.add_argument("--y-max", type=float, default=0.51)
    ap.add_argument("--analysis-start", type=int, default=10)
    ap.add_argument("--analysis-end", type=int, default=40)
    ap.add_argument("--header", default="Phase redistribution disrupts read-coordinate period-3 coherence")
    args = ap.parse_args()

    data = {k: load(getattr(args, k)) for k in DEFAULTS}

    fig = plt.figure(figsize=(11.5, 7.4))
    gs = fig.add_gridspec(2, 6, hspace=0.40, wspace=0.45)
    axes = [
        fig.add_subplot(gs[0, 0:3]),
        fig.add_subplot(gs[0, 3:6]),
        fig.add_subplot(gs[1, 0:2]),
        fig.add_subplot(gs[1, 2:4]),
        fig.add_subplot(gs[1, 4:6]),
    ]

    for ax, (key, letter, title, n) in zip(axes, PANELS):
        draw(ax, data[key], letter, title, n, args.y_min, args.y_max,
             args.analysis_start, args.analysis_end)

    # One y-axis label for each row, as requested.
    fig.text(0.022, 0.675, "Fraction", rotation=90, va="center", ha="center", fontsize=11)
    fig.text(0.022, 0.255, "Fraction", rotation=90, va="center", ha="center", fontsize=11)

    handles = [
        plt.Line2D([], [], linestyle="none", marker="o", markersize=5,
                   color=COLORS[b], label=b) for b in BASES
    ]

    # Header, then shared legend beneath it.
    fig.suptitle(args.header, fontsize=14, y=0.992)
    fig.legend(handles=handles, labels=BASES, loc="upper center",
               bbox_to_anchor=(0.5, 0.958), ncol=4,
               frameon=False, fontsize=9)

    fig.subplots_adjust(top=0.88, bottom=0.08, left=0.075, right=0.985)

    out = args.output_prefix
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out.with_suffix('.png')}")
    print(f"Wrote {out.with_suffix('.pdf')}")

if __name__ == "__main__":
    main()

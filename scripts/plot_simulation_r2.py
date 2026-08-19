#!/usr/bin/env python3
"""Plot period-3 simulation results from machine-readable source data.

Required input columns:
    genome, bias_5prime, bias_3prime, replicate, base, r2

Each replicate/condition must contain exactly one row for each of A, T, G, and C.
The script first averages R^2 across the four bases within each replicate, then
plots the mean across independent replicates.  Error bars are 95% t confidence
intervals across replicate-level means (zero-width for a single replicate).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Optional, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import t

SCRIPT_VERSION = "2.0.0-audited"
BASES = ("A", "T", "G", "C")
REQUIRED = {"genome", "bias_5prime", "bias_3prime", "replicate", "base", "r2"}


def load_and_validate(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = REQUIRED - set(df.columns)
    if missing:
        raise ValueError(f"Input CSV is missing columns: {sorted(missing)}")
    df = df.loc[:, ["genome", "bias_5prime", "bias_3prime", "replicate", "base", "r2"]].copy()
    df["base"] = df["base"].astype(str).str.upper()
    if (~df["base"].isin(BASES)).any():
        bad = sorted(df.loc[~df["base"].isin(BASES), "base"].unique())
        raise ValueError(f"Unsupported bases: {bad}")
    for column in ("bias_5prime", "bias_3prime", "r2"):
        df[column] = pd.to_numeric(df[column], errors="raise")
    if ((df["r2"] < 0) | (df["r2"] > 1) | ~np.isfinite(df["r2"])).any():
        raise ValueError("R^2 values must be finite and lie in [0, 1].")
    if ((df["bias_5prime"] < 0.5) | (df["bias_5prime"] > 1.0)).any():
        raise ValueError("bias_5prime must lie in [0.5, 1.0].")
    if ((df["bias_3prime"] < 0.5) | (df["bias_3prime"] > 1.0)).any():
        raise ValueError("bias_3prime must lie in [0.5, 1.0].")

    keys = ["genome", "bias_5prime", "bias_3prime", "replicate", "base"]
    if df.duplicated(keys).any():
        raise ValueError("Duplicate genome/bias/replicate/base rows were found.")
    counts = df.groupby(keys[:-1])["base"].agg(lambda values: tuple(sorted(values)))
    expected = tuple(sorted(BASES))
    bad = counts[counts != expected]
    if not bad.empty:
        raise ValueError(
            "Every replicate/condition must contain A, T, G, and C exactly once; "
            f"invalid groups: {list(bad.index[:10])}"
        )
    return df


def summarize(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    replicate = (
        df.groupby(["genome", "bias_5prime", "bias_3prime", "replicate"], as_index=False)
        .agg(mean_r2_across_bases=("r2", "mean"), sd_r2_across_bases=("r2", "std"))
    )

    rows = []
    for keys, group in replicate.groupby(["genome", "bias_5prime", "bias_3prime"], sort=True):
        values = group["mean_r2_across_bases"].to_numpy(dtype=float)
        n = len(values)
        mean = float(values.mean())
        sd = float(values.std(ddof=1)) if n > 1 else 0.0
        sem = sd / np.sqrt(n) if n > 1 else 0.0
        half_width = float(t.ppf(0.975, df=n - 1) * sem) if n > 1 else 0.0
        rows.append(
            {
                "genome": keys[0],
                "bias_5prime": keys[1],
                "bias_3prime": keys[2],
                "n_replicates": n,
                "mean_r2": mean,
                "sd_between_replicates": sd,
                "sem_between_replicates": sem,
                "ci95_lower": max(0.0, mean - half_width),
                "ci95_upper": min(1.0, mean + half_width),
            }
        )
    return replicate, pd.DataFrame(rows)


def plot_one_genome(summary: pd.DataFrame, genome: str, output: Path) -> None:
    part = summary.loc[summary["genome"] == genome].copy()
    if part.empty:
        raise ValueError(f"No rows for genome '{genome}'.")
    part = part.sort_values(["bias_5prime", "bias_3prime"]).reset_index(drop=True)
    labels = [f"{row.bias_5prime:.2f}/{row.bias_3prime:.2f}" for row in part.itertuples()]
    values = part["mean_r2"].to_numpy(dtype=float)
    lower = values - part["ci95_lower"].to_numpy(dtype=float)
    upper = part["ci95_upper"].to_numpy(dtype=float) - values

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(part))
    ax.errorbar(x, values, yerr=np.vstack([lower, upper]), marker="o", linewidth=2, capsize=4)
    ax.set_xticks(x, labels=labels)
    ax.set_xlabel("5' bias / 3' bias")
    ax.set_ylabel("Mean period-3 R² across bases")
    ax.set_ylim(0, 1)
    ax.set_title(f"{genome}: simulated fragmentation")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=400, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    df = load_and_validate(args.input_csv)
    replicate, summary = summarize(df)
    args.outdir.mkdir(parents=True, exist_ok=True)
    replicate_path = args.outdir / "simulation_replicate_means.csv"
    summary_path = args.outdir / "simulation_condition_summary.csv"
    replicate.to_csv(replicate_path, index=False)
    summary.to_csv(summary_path, index=False)
    plot_paths = []
    for genome in summary["genome"].drop_duplicates():
        safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in str(genome))
        plot_path = args.outdir / f"simulation_r2_{safe}.png"
        plot_one_genome(summary, str(genome), plot_path)
        plot_paths.extend([str(plot_path), str(plot_path.with_suffix(".pdf"))])
    manifest = {
        "script": "plot_simulation_r2.py",
        "script_version": SCRIPT_VERSION,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_csv": str(args.input_csv.resolve()),
        "summary_rule": (
            "average R2 across A/T/G/C within each replicate, then mean and 95% t CI "
            "across independent replicate means"
        ),
        "input_rows": len(df),
        "replicate_rows": len(replicate),
        "condition_rows": len(summary),
        "outputs": {
            "replicate_means_csv": str(replicate_path),
            "condition_summary_csv": str(summary_path),
            "plots": plot_paths,
        },
    }
    manifest_path = args.outdir / "simulation_plot_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"Wrote source summaries and {summary['genome'].nunique()} genome plot(s) to {args.outdir}")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

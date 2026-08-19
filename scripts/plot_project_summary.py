#!/usr/bin/env python3
"""Join period-3 library summaries to AncientMetagenomeDir metadata and plot groups."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

SCRIPT_VERSION = "2.0.0-audited"
DEFAULT_PARAMETERS = (
    "project_name",
    "publication_year",
    "strand_type",
    "library_polymerase",
    "library_treatment",
    "instrument_model",
)
UNKNOWN_LABEL = "Unknown"
MISSING_LABELS = {
    "unknown", "unspecified", "not specified", "not available", "na", "n/a", "null", ""
}


def normalize_category_value(value: object) -> str:
    if pd.isna(value):
        return UNKNOWN_LABEL
    text = str(value).strip()
    return UNKNOWN_LABEL if text.lower() in MISSING_LABELS else text


def clean_categories(series: pd.Series, min_count: int) -> pd.Series:
    values = series.map(normalize_category_value)
    counts = values.value_counts()
    keep = set(counts[counts >= min_count].index)
    return values.where(values.isin(keep), "Other")


def _assert_unique_nonmissing(df: pd.DataFrame, column: str, table_name: str) -> None:
    if column not in df:
        raise ValueError(f"{table_name} is missing required column '{column}'.")
    duplicated = df[column].dropna().duplicated(keep=False)
    if duplicated.any():
        values = sorted(df.loc[df[column].dropna().index[duplicated], column].astype(str).unique())
        preview = ", ".join(values[:10])
        raise ValueError(f"{table_name}.{column} is not unique (examples: {preview}).")


def build_analysis_table(
    r2: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    internal_project_label: str | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"sample", "category", "max_r2"}
    missing = required - set(r2.columns)
    if missing:
        raise ValueError(f"R2 table is missing columns: {sorted(missing)}")

    _assert_unique_nonmissing(r2, "sample", "R2 table")
    _assert_unique_nonmissing(metadata, "archive_data_accession", "metadata")
    if "library_name" in metadata:
        nonmissing = metadata.loc[metadata["library_name"].notna()].copy()
        _assert_unique_nonmissing(nonmissing, "library_name", "metadata")

    external_source = r2.loc[r2["category"] == "external_ancient"].copy()
    external = external_source.merge(
        metadata,
        left_on="sample",
        right_on="archive_data_accession",
        how="left",
        validate="one_to_one",
        indicator="_metadata_merge",
        suffixes=("", "_metadata"),
    )

    internal_source = r2.loc[r2["category"] == "internal_ancient"].copy()
    if internal_source.empty:
        internal = internal_source
    elif "library_name" not in metadata:
        internal = internal_source.copy()
        internal["_metadata_merge"] = "left_only"
    else:
        internal = internal_source.merge(
            metadata,
            left_on="sample",
            right_on="library_name",
            how="left",
            validate="one_to_one",
            indicator="_metadata_merge",
            suffixes=("", "_metadata"),
        )

    if internal_project_label is not None and not internal.empty:
        if "project_name" not in internal:
            internal["project_name"] = internal_project_label
        else:
            missing_project = internal["project_name"].isna()
            internal.loc[missing_project, "project_name"] = internal_project_label

    combined = pd.concat([external, internal], ignore_index=True, sort=False)
    unmatched = combined.loc[combined["_metadata_merge"] != "both", ["sample", "category", "_metadata_merge"]]
    return combined, unmatched.reset_index(drop=True)


def create_dotplot(
    df: pd.DataFrame,
    column: str,
    outdir: Path,
    *,
    min_count: int,
    seed: int,
) -> Path:
    tmp = df.loc[df["max_r2"].notna()].copy()
    tmp[column] = clean_categories(tmp[column], min_count)
    categories = list(tmp[column].dropna().unique())
    if not categories:
        raise ValueError(f"No categories available for '{column}'.")

    order = sorted(categories, key=lambda c: (tmp.loc[tmp[column] == c, "max_r2"].median(), str(c)))
    counts = tmp[column].value_counts()
    palette = sns.color_palette("tab20", n_colors=max(3, len(order)))
    color_map = dict(zip(order, palette[: len(order)]))
    rng = np.random.default_rng(seed)

    fig_height = max(6.0, 0.38 * len(order))
    fig, ax = plt.subplots(figsize=(11, fig_height))
    for i, category in enumerate(order):
        values = tmp.loc[tmp[column] == category, "max_r2"].to_numpy(dtype=float)
        if values.size == 0:
            continue
        ax.scatter(
            values,
            i + rng.uniform(-0.22, 0.22, size=values.size),
            s=18,
            alpha=0.75,
            color=color_map[category],
            edgecolors="white",
            linewidths=0.3,
            zorder=2,
        )
        ax.plot(
            float(np.median(values)), i, marker="|", markersize=14,
            markeredgewidth=2, color="black", zorder=3,
        )

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{category}  (n={counts[category]})" for category in order])
    ax.set_xlabel("Maximum mean R² across eligible read lengths", fontsize=13)
    ax.set_ylabel(column.replace("_", " ").title(), fontsize=13)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_title(
        f"Period-3 variance explained\nGrouped by {column.replace('_', ' ')}",
        fontsize=15,
    )
    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / f"dotplot_by_{column}.png"
    fig.savefig(outfile, dpi=400, bbox_inches="tight")
    plt.close(fig)
    return outfile


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r2-csv", type=Path, required=True)
    parser.add_argument("--metadata-tsv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--parameters", nargs="+", default=list(DEFAULT_PARAMETERS))
    parser.add_argument("--min-libraries-per-category", type=int, default=5)
    parser.add_argument("--jitter-seed", type=int, default=0)
    parser.add_argument(
        "--internal-project-label",
        default="Slon2017",
        help="Project label used only when an internal sample lacks a metadata match; pass an empty string to disable.",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    if args.min_libraries_per_category < 1:
        raise ValueError("--min-libraries-per-category must be positive.")

    sns.set_theme(style="whitegrid", context="talk")
    r2 = pd.read_csv(args.r2_csv)
    metadata = pd.read_csv(args.metadata_tsv, sep="\t", low_memory=False)
    label = args.internal_project_label.strip() or None
    analysis, unmatched = build_analysis_table(r2, metadata, internal_project_label=label)

    args.outdir.mkdir(parents=True, exist_ok=True)
    joined_path = args.outdir / "period3_metadata_joined_source_data.csv"
    unmatched_path = args.outdir / "unmatched_libraries.csv"
    analysis.to_csv(joined_path, index=False)
    unmatched.to_csv(unmatched_path, index=False)

    print(f"Libraries in R2 table: {len(r2):,}")
    print(f"Libraries retained in joined analysis: {len(analysis):,}")
    print(f"Libraries without direct metadata match: {len(unmatched):,}")

    plotted = []
    skipped = []
    for parameter in args.parameters:
        if parameter not in analysis:
            print(f"WARNING: '{parameter}' not found in joined table; skipping.")
            skipped.append(parameter)
            continue
        outfile = create_dotplot(
            analysis,
            parameter,
            args.outdir,
            min_count=args.min_libraries_per_category,
            seed=args.jitter_seed,
        )
        plotted.append(str(outfile))
        print(f"Saved: {outfile}")

    manifest = {
        "script": "plot_project_summary.py",
        "script_version": SCRIPT_VERSION,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "r2_csv": str(args.r2_csv.resolve()),
        "metadata_tsv": str(args.metadata_tsv.resolve()),
        "parameters_requested": list(args.parameters),
        "parameters_skipped": skipped,
        "min_libraries_per_category": args.min_libraries_per_category,
        "jitter_seed": args.jitter_seed,
        "internal_project_label": label,
        "r2_rows": len(r2),
        "joined_rows": len(analysis),
        "unmatched_rows": len(unmatched),
        "outputs": {
            "joined_source_data_csv": str(joined_path),
            "unmatched_csv": str(unmatched_path),
            "plots": plotted,
        },
    }
    manifest_path = args.outdir / "metadata_plot_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()

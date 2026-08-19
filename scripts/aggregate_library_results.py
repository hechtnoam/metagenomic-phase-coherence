#!/usr/bin/env python3
"""Aggregate per-length period-3 JSON files into library-level summaries.

For every library, the primary historical statistic is retained exactly as used in
this project: mean R^2 across A/T/G/C at each eligible read length, followed by
the maximum over eligible lengths.  The output also records the number and range
of eligible lengths so the multiplicity implicit in that maximum is visible.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

import pandas as pd

SCRIPT_VERSION = "2.0.0-audited"
BASES = ("A", "T", "G", "C")
ACCESSION_RE = re.compile(r"^(?:SRR|ERR|DRR|PRJ)", re.IGNORECASE)


@dataclass(frozen=True)
class LengthResult:
    sample: str
    category: str
    read_length: int
    read_count: int | None
    observed_read_count: int | None
    mean_r2: float
    min_base_r2: float
    max_base_r2: float
    all_four_bases_present: bool
    source_json: str


def determine_category(sample_name: str) -> str:
    """Classify accession-like sample names as external, otherwise internal."""

    return "external_ancient" if ACCESSION_RE.match(sample_name) else "internal_ancient"


def _sample_from_json_path(path: Path) -> str:
    """Resolve the sample directory in <sample>/stats/period3/kmer/L*/... layouts."""

    for parent in path.parents:
        if parent.name == "stats":
            return parent.parent.name
    raise ValueError(f"Could not locate a 'stats' ancestor for {path}")


def parse_result(path: Path) -> LengthResult | None:
    with path.open() as handle:
        data = json.load(handle)

    if data.get("status") == "skipped":
        return None

    read_length = data.get("read_length")
    if read_length is None:
        match = re.fullmatch(r"L(\d+)", path.parent.name)
        if not match:
            raise ValueError("Missing read_length and parent directory is not L<number>.")
        read_length = int(match.group(1))
    read_length = int(read_length)

    fits = data.get("fits")
    if not isinstance(fits, list) or not fits:
        raise ValueError("Missing non-empty 'fits' list.")

    r2_by_base: dict[str, float] = {}
    for fit in fits:
        base = str(fit.get("base", "")).upper()
        if base in BASES and fit.get("r2") is not None:
            r2_by_base[base] = float(fit["r2"])

    missing = [base for base in BASES if base not in r2_by_base]
    if missing:
        raise ValueError(f"Missing R^2 values for bases: {', '.join(missing)}")

    values = [r2_by_base[base] for base in BASES]
    if any((not pd.notna(value)) or value < 0.0 or value > 1.0 + 1e-12 for value in values):
        raise ValueError(f"Invalid R^2 values: {values}")

    sample = _sample_from_json_path(path)
    return LengthResult(
        sample=sample,
        category=determine_category(sample),
        read_length=read_length,
        read_count=int(data["read_count"]) if data.get("read_count") is not None else None,
        observed_read_count=(
            int(data["observed_read_count"])
            if data.get("observed_read_count") is not None
            else None
        ),
        mean_r2=float(sum(values) / len(values)),
        min_base_r2=float(min(values)),
        max_base_r2=float(max(values)),
        all_four_bases_present=True,
        source_json=str(path.resolve()),
    )


def aggregate(json_files: Iterable[Path]) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    rows: list[dict[str, object]] = []
    errors: list[str] = []
    for path in sorted(json_files):
        try:
            parsed = parse_result(path)
            if parsed is not None:
                rows.append(asdict(parsed))
        except Exception as exc:  # continue across libraries but preserve an auditable error log
            errors.append(f"{path}: {type(exc).__name__}: {exc}")

    all_lengths = pd.DataFrame(rows)
    if all_lengths.empty:
        return all_lengths, pd.DataFrame(), errors

    # Guard against duplicated reruns occupying the same sample/read-length cell.
    duplicated = all_lengths.duplicated(["sample", "read_length"], keep=False)
    if duplicated.any():
        duplicate_rows = all_lengths.loc[duplicated, ["sample", "read_length", "source_json"]]
        raise ValueError(
            "Duplicate sample/read-length results were found:\n"
            + duplicate_rows.to_string(index=False)
        )

    all_lengths = all_lengths.sort_values(["sample", "read_length"]).reset_index(drop=True)

    # Deterministic tie-breaking: highest mean R^2, then most analyzed reads,
    # then shortest read length. Missing read counts sort last.
    ranked = all_lengths.copy()
    ranked["_read_count_rank"] = ranked["read_count"].fillna(-1)
    ranked = ranked.sort_values(
        ["sample", "mean_r2", "_read_count_rank", "read_length"],
        ascending=[True, False, False, True],
    )
    # GroupBy.first selects the first non-null value independently in each
    # column and can splice values from different ranked rows. Keep the first
    # complete ranked row instead.
    best = ranked.groupby("sample", sort=True, group_keys=False).head(1).copy()
    best = best.drop(columns=["_read_count_rank"])

    per_sample = all_lengths.groupby("sample", sort=True).agg(
        n_eligible_lengths=("read_length", "size"),
        min_eligible_length=("read_length", "min"),
        max_eligible_length=("read_length", "max"),
        median_mean_r2=("mean_r2", "median"),
        mean_mean_r2_across_lengths=("mean_r2", "mean"),
    )

    best = best.merge(per_sample, on="sample", how="left", validate="one_to_one")
    best = best.rename(
        columns={
            "read_length": "best_length",
            "read_count": "best_length_read_count",
            "observed_read_count": "best_length_observed_read_count",
            "mean_r2": "max_r2",
            "min_base_r2": "best_length_min_base_r2",
            "max_base_r2": "best_length_max_base_r2",
            "source_json": "best_length_source_json",
        }
    )
    best = best[
        [
            "sample",
            "category",
            "best_length",
            "best_length_read_count",
            "best_length_observed_read_count",
            "max_r2",
            "best_length_min_base_r2",
            "best_length_max_base_r2",
            "n_eligible_lengths",
            "min_eligible_length",
            "max_eligible_length",
            "median_mean_r2",
            "mean_mean_r2_across_lengths",
            "best_length_source_json",
        ]
    ].sort_values("sample")
    return all_lengths, best.reset_index(drop=True), errors


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-outdir",
        type=Path,
        required=True,
        help="Root containing <sample>/stats/period3/kmer/L*/k1_period3.json.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Library-level CSV (default: <base-outdir>/population_figures/best_length_per_library.csv).",
    )
    parser.add_argument(
        "--all-lengths-output",
        type=Path,
        default=None,
        help="Per-length source-data CSV (default: beside --output).",
    )
    parser.add_argument(
        "--errors-output",
        type=Path,
        default=None,
        help="Text error log (default: beside --output).",
    )
    parser.add_argument(
        "--manifest-output",
        type=Path,
        default=None,
        help="Aggregation manifest JSON (default: beside --output).",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    base = args.base_outdir.resolve()
    output = args.output or (base / "population_figures" / "best_length_per_library.csv")
    all_output = args.all_lengths_output or output.with_name("period3_all_eligible_lengths.csv")
    errors_output = args.errors_output or output.with_name("period3_aggregation_errors.txt")
    manifest_output = args.manifest_output or output.with_name("period3_aggregation_manifest.json")

    files = list(base.rglob("stats/period3/kmer/L*/k1_period3.json"))
    if not files:
        raise FileNotFoundError(f"No k1_period3.json files found under {base}")

    all_lengths, best, errors = aggregate(files)
    if best.empty:
        raise RuntimeError("No valid period-3 results were available after parsing.")

    output.parent.mkdir(parents=True, exist_ok=True)
    all_output.parent.mkdir(parents=True, exist_ok=True)
    errors_output.parent.mkdir(parents=True, exist_ok=True)
    manifest_output.parent.mkdir(parents=True, exist_ok=True)
    all_lengths.to_csv(all_output, index=False)
    best.to_csv(output, index=False)
    errors_output.write_text("\n".join(errors) + ("\n" if errors else ""))
    manifest = {
        "script": "aggregate_library_results.py",
        "script_version": SCRIPT_VERSION,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "base_outdir": str(base),
        "source_pattern": "stats/period3/kmer/L*/k1_period3.json",
        "selection_statistic": (
            "mean R2 across A/T/G/C at each eligible exact read length, then maximum "
            "over lengths; ties by larger analyzed read count then shorter length"
        ),
        "source_json_files_found": len(files),
        "valid_length_rows": len(all_lengths),
        "libraries": int(best["sample"].nunique()),
        "parse_errors": len(errors),
        "outputs": {
            "library_summary_csv": str(output),
            "all_lengths_csv": str(all_output),
            "errors_txt": str(errors_output),
        },
    }
    manifest_output.write_text(json.dumps(manifest, indent=2))

    print(f"Parsed {len(all_lengths):,} length-level results from {best['sample'].nunique():,} libraries.")
    print(f"Library summary: {output}")
    print(f"Per-length source data: {all_output}")
    print(f"Parse errors: {len(errors):,} ({errors_output})")
    print(f"Manifest: {manifest_output}")


if __name__ == "__main__":
    main()

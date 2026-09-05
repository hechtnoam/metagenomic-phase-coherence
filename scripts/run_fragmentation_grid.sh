#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON:-python3}"
OUTROOT="${OUTROOT:-$ROOT/results/fragmentation_simulation}"

: "${HUMAN_REF:?Set HUMAN_REF to the human reference FASTA}"
: "${PSEUDO_REF:?Set PSEUDO_REF to the Pseudomonas reference FASTA}"

N="${NUM_READS:-100000}"
L="${READ_LENGTH:-59}"
PERMS="${PERMS:-10000}"
WORKERS="${P3_WORKERS:-10}"

SIM="$ROOT/scripts/simulate_fragmentation.py"
P3="$ROOT/scripts/period3_library_pipeline.py"
COLLECT="$ROOT/scripts/collect_simulation_results.py"
PLOT="$ROOT/scripts/plot_simulation_r2.py"

SEEDS=(1 2 3 4 5 6 7 8 9 10)
CONDS=(
  "random:0.50:0.50"
  "b055_050:0.55:0.50"
  "b070_050:0.70:0.50"
  "b080_050:0.80:0.50"
  "b090_050:0.90:0.50"
  "b100_050:1.00:0.50"
  "b100_100:1.00:1.00"
)

for path in "$SIM" "$P3" "$COLLECT" "$PLOT" "$HUMAN_REF" "$PSEUDO_REF"; do
  if [[ ! -e "$path" ]]; then
    echo "Missing required file: $path" >&2
    exit 1
  fi
done

run_genome() {
  local genome="$1"
  local ref="$2"
  local spec label b5 b3 seed run fasta metadata p3json

  for spec in "${CONDS[@]}"; do
    IFS=: read -r label b5 b3 <<< "$spec"
    for seed in "${SEEDS[@]}"; do
      run="$OUTROOT/$genome/$label/replicate_$seed"
      fasta="$run/simulated_reads.fasta"
      metadata="$fasta.metadata.json"
      p3json="$run/period3/simulated_reads/stats/period3/kmer/L$L/k1_period3.json"
      mkdir -p "$run"

      echo "=== $genome | $label | bias=$b5/$b3 | replicate=$seed ==="

      if [[ ! -s "$fasta" || ! -s "$metadata" ]]; then
        "$PYTHON" "$SIM" \
          --input "$ref" \
          --output-dir "$run" \
          --num-reads "$N" \
          --read-length "$L" \
          --bias-5prime "$b5" \
          --bias-3prime "$b3" \
          --reverse-complement-probability 0.5 \
          --seed "$seed"
      else
        echo "[skip simulation] existing FASTA and metadata"
      fi

      if [[ ! -s "$p3json" ]]; then
        "$PYTHON" "$P3" all "$fasta" \
          --base-outdir "$run/period3" \
          --no-collapse \
          --lengths "$L" \
          --workers "$WORKERS" \
          --p3-start 10 \
          --p3-end 40 \
          --p3-perms "$PERMS"
      else
        echo "[skip period3] existing statistics"
      fi
    done
  done
}

mkdir -p "$OUTROOT"
run_genome human "$HUMAN_REF"
run_genome pseudomonas "$PSEUDO_REF"

"$PYTHON" "$COLLECT" \
  --root "$OUTROOT" \
  --output "$OUTROOT/simulation_results.csv" \
  --read-length "$L"

"$PYTHON" "$PLOT" \
  --input-csv "$OUTROOT/simulation_results.csv" \
  --outdir "$OUTROOT/final_plots"

echo "DONE: $OUTROOT/final_plots"

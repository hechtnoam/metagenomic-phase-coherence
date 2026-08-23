#!/usr/bin/env bash
set -euo pipefail
ROOT="$(pwd)"; OUTROOT="${OUTROOT:-$ROOT/results/fragmentation_simulation}"
HUMAN_REF="/home/MainNas/ncbi_genomes/sedaDNA/data/Human/GCF_000001405.26_GRCh38_genomic.fna"
PSEUDO_REF="/home/MainNas/ncbi_genomes/sedaDNA/data/Pseudomonas/pseudomonas_reference.fna"
N="${NUM_READS:-100000}"; L="${READ_LENGTH:-59}"; PERMS="${PERMS:-10000}"; WORKERS="${P3_WORKERS:-10}"
SEEDS=(1 2 3 4 5 6 7 8 9 10)
CONDS=("random:0.50:0.50" "b055_050:0.55:0.50" "b070_050:0.70:0.50" "b080_050:0.80:0.50" "b090_050:0.90:0.50" "b100_050:1.00:0.50" "b100_100:1.00:1.00")
for x in scripts/simulate_fragmentation.py scripts/period3_library_pipeline.py scripts/collect_simulation_results.py scripts/plot_simulation_r2.py "$HUMAN_REF" "$PSEUDO_REF"; do [[ -e "$x" ]] || { echo "Missing: $x" >&2; exit 1; }; done
run_genome(){ local genome="$1" ref="$2"; for spec in "${CONDS[@]}"; do IFS=: read -r label b5 b3 <<< "$spec"; for seed in "${SEEDS[@]}"; do
 d="$OUTROOT/$genome/$label/replicate_$seed"; f="$d/simulated_reads.fasta"; m="$f.metadata.json"; p="$d/period3/simulated_reads/stats/period3/kmer/L${L}/k1_period3.json"; mkdir -p "$d"
 echo "=== $genome $label ($b5/$b3) replicate $seed ==="
 if [[ ! -s "$f" || ! -s "$m" ]]; then python3 scripts/simulate_fragmentation.py -i "$ref" -o "$d" -n "$N" -l "$L" --bias-5prime "$b5" --bias-3prime "$b3" --reverse-complement-probability 0.5 --seed "$seed"; else echo "[skip simulation]"; fi
 if [[ ! -s "$p" ]]; then python3 scripts/period3_library_pipeline.py all "$f" --base-outdir "$d/period3" --no-collapse --lengths "$L" --workers "$WORKERS" --p3-start 10 --p3-end 40 --p3-perms "$PERMS"; else echo "[skip period3]"; fi
 done; done; }
mkdir -p "$OUTROOT"; run_genome human "$HUMAN_REF"; run_genome pseudomonas "$PSEUDO_REF"
python3 scripts/collect_simulation_results.py --root "$OUTROOT" --output "$OUTROOT/simulation_results.csv" --read-length "$L"
python3 scripts/plot_simulation_r2.py --input-csv "$OUTROOT/simulation_results.csv" --outdir "$OUTROOT/final_plots"
echo "DONE: $OUTROOT/final_plots"
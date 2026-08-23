#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"

SIM="$ROOT/scripts/simulate_fragmentation.py"
P3="$ROOT/scripts/period3_library_pipeline.py"

HUMAN="/home/MainNas/ncbi_genomes/sedaDNA/data/Human/GCF_000001405.26_GRCh38_genomic.fna"
PSEUDO="/home/MainNas/ncbi_genomes/sedaDNA/data/Pseudomonas/pseudomonas_reference.fna"

# IMPORTANT: completely separate from the definitive 140-run experiment.
OUT="$ROOT/results/fragmentation_forensic_forward_only"

READS=100000
LENGTH=59
PERMS=10000

SEEDS=(1 2 3 4 5 6 7 8 9 10)

# condition:bias5:bias3
CONDITIONS=(
  "random:0.50:0.50"
  "b100_050:1.00:0.50"
  "b100_100:1.00:1.00"
)

run_genome () {
    genome="$1"
    ref="$2"

    for spec in "${CONDITIONS[@]}"; do

        IFS=: read -r condition b5 b3 <<< "$spec"

        for seed in "${SEEDS[@]}"; do

            RUN="$OUT/$genome/$condition/replicate_${seed}"
            FASTA="$RUN/simulated_reads.fasta"
            META="$FASTA.metadata.json"
            P3OUT="$RUN/period3"

            mkdir -p "$RUN"

            echo
            echo "=================================================="
            echo "$genome | $condition | seed=$seed"
            echo "bias=$b5/$b3 | FORWARD ONLY"
            echo "=================================================="

            if [[ ! -s "$FASTA" || ! -s "$META" ]]; then

                python3 "$SIM" \
                  -i "$ref" \
                  -o "$RUN" \
                  -n "$READS" \
                  -l "$LENGTH" \
                  --bias-5prime "$b5" \
                  --bias-3prime "$b3" \
                  --reverse-complement-probability 0 \
                  --seed "$seed"

            else
                echo "[SKIP simulation] already exists"
            fi

            P3JSON="$P3OUT/simulated_reads/stats/period3/kmer/L59/k1_period3.json"

            if [[ ! -s "$P3JSON" ]]; then

                python3 "$P3" all "$FASTA" \
                  --base-outdir "$P3OUT" \
                  --no-collapse \
                  --lengths 59 \
                  --workers 10 \
                  --p3-start 10 \
                  --p3-end 40 \
                  --p3-perms "$PERMS"

            else
                echo "[SKIP period3] already exists"
            fi

        done
    done
}

run_genome human "$HUMAN"
run_genome pseudomonas "$PSEUDO"

echo
echo "=================================================="
echo "FORWARD-ONLY FORENSIC RUN COMPLETE"
echo "Results:"
echo "$OUT"
echo "=================================================="
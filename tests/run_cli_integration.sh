#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CODE="$ROOT/scripts"
WORK="$ROOT/cli_integration_current"
rm -rf "$WORK"
mkdir -p "$WORK"

python - "$WORK/input.fa" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1])
seqs=[]
for i in range(180):
    phase=i%3
    motif=("ATG","TGA","GAT")[phase]
    seq=(motif*20)[:60]
    seqs.append(f">r{i}\n{seq}\n")
p.write_text("".join(seqs))
PY

python "$CODE/period3_library_pipeline.py" all "$WORK/input.fa" \
  --base-outdir "$WORK/period3" --lengths 60 --min-read-count 1 \
  --workers 1 --p3-perms 99 --p3-seed 101 --format png --quiet

CSV="$WORK/period3/input/csv/kmer/L60/k1.csv"
python "$CODE/period3_statistics.py" test one "$CSV" --start 10 --end 40 --perms 99 --seed 101 > "$WORK/period3_standalone.json"

python "$CODE/dft_mod3.py" --input "$WORK/input.fa" --length 60 \
  --out-prefix "$WORK/dft/test" --window-start 10 --window-end 39 --max-reads 150

python "$CODE/phase_redistribution.py" --input "$WORK/input.fa" --outdir "$WORK/trim" \
  --length 60 --mode random --seed 9
python "$CODE/phase_redistribution.py" --input "$WORK/input.fa" --outdir "$WORK/trim" \
  --length 60 --mode fixed --fixed-five-trim 1 --seed 9

python "$CODE/six_state_em.py" --input "$WORK/input.fa" --input-type fasta \
  --out "$WORK/em" --lengths 60 --restarts 2 --threads 1 --max-iter 100 \
  --canonicalize full6 --progress-every 0

python - "$WORK/ref.fa" <<'PY'
from pathlib import Path
import sys
Path(sys.argv[1]).write_text(">contig1\n" + "AAGGCT"*200 + "\n>contig2\n" + "GGTTAC"*120 + "\n")
PY
python "$CODE/simulate_fragmentation.py" --input "$WORK/ref.fa" --output-dir "$WORK/sim" \
  --num-reads 200 --read-length 30 --bias-5prime 0.7 --bias-3prime 0.5 \
  --seed 22 --chunk-reads 37

python "$CODE/aggregate_library_results.py" --base-outdir "$WORK/period3" \
  --output "$WORK/aggregate/best.csv"

cat > "$WORK/meta.tsv" <<'TSV'
archive_data_accession	library_name	project_name
NA	input	SyntheticProject
TSV
python "$CODE/plot_project_summary.py" \
  --r2-csv "$WORK/aggregate/best.csv" --metadata-tsv "$WORK/meta.tsv" \
  --outdir "$WORK/meta_plots" --parameters project_name --min-libraries-per-category 1

python - "$WORK/simulation_results.csv" <<'PY'
import sys, pandas as pd
rows = []
ORDER = [('random', .5, .5), ('b055_050', .55, .5), ('b070_050', .7, .5), 
         ('b080_050', .8, .5), ('b090_050', .9, .5), ('b100_050', 1., .5), ('b100_100', 1., 1.)]
for genome in ("pseudomonas", "human"):
    for cond, b5, b3 in ORDER:
        for rep in (1, 2):
            for base in ("A", "T", "G", "C"):
                rows.append((genome, cond, b5, b3, rep, base, 0.5))
df = pd.DataFrame(rows, columns=["genome", "condition", "bias_5prime", "bias_3prime", "replicate", "base", "r2"])
df.to_csv(sys.argv[1], index=False)
PY
python "$CODE/plot_simulation_r2.py" --input-csv "$WORK/simulation_results.csv" --outdir "$WORK/sim_plot"

test -s "$WORK/period3/input/run_manifest.json"
test -s "$WORK/dft/test_L60_dft_metadata.json"
test -s "$WORK/trim/input.phase_randomized.fasta.manifest.json"
test -s "$WORK/em/run_metadata.json"
test -s "$WORK/sim/simulated_reads.fasta.metadata.json"
test -s "$WORK/aggregate/period3_aggregation_manifest.json"
test -s "$WORK/meta_plots/metadata_plot_manifest.json"
test -s "$WORK/sim_plot/simulation_plot_manifest.json"

BAM_TEST="$ROOT/tests/data/test_A2424.bam"
if [[ -f "$BAM_TEST" ]]; then
    python "$CODE/six_state_em.py" --input "$BAM_TEST" --input-type bam \
      --out "$WORK/em_bam" --lengths 59 --restarts 1 --threads 1 --max-iter 10 \
      --canonicalize full6 --progress-every 0
    test -s "$WORK/em_bam/run_metadata.json"
fi

printf 'PASS\tcli_integration\tall runnable reviewed command-line entry points produced expected outputs and manifests\n' > "$ROOT/CLI_INTEGRATION_RESULTS.txt"
cat "$ROOT/CLI_INTEGRATION_RESULTS.txt"

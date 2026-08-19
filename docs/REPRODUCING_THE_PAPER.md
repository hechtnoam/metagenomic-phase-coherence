# Reproducibility rerun plan

## 1. Freeze inputs and environment

- Create a complete accession-level table with library identifier, project, source URL, file checksum, input type, paired-end handling, read-count threshold, and any preprocessing.
- Record the exact AncientMetagenomeDir release or Git commit.
- Record every reference FASTA and GFF accession/version and its SHA-256 checksum.
- Create the environment from `environment.yml`, record the resolved package lock, and archive the exact code revision.
- Run `python tests/run_smoke_tests.py` and `bash tests/run_cli_integration.sh` in the production environment. Add BAM/GFF tests before any data rerun.

## 2. Empirical cycle composition

Template:

```bash
python reviewed_code/period3_library_pipeline.py all INPUT \
  --base-outdir OUTPUT_ROOT \
  --min-read-count 40000 \
  --workers N \
  --p3-start 10 --p3-end 40 \
  --p3-perms 10000 --p3-seed 12345
```

Do not enable exact-sequence collapse unless it is scientifically justified, prespecified, and reported. Retain every run manifest, cycle CSV, statistics JSON, exclusion count, and log.

## 3. Phase redistribution and fixed-trim controls

For every prespecified seed, run the random control and the three fixed controls on the identical filtered source read set:

```bash
python reviewed_code/trim_reads.py --input SOURCE --outdir OUT --length L --mode random --seed SEED
python reviewed_code/trim_reads.py --input SOURCE --outdir OUT --length L --mode fixed --fixed-five-trim 0 --seed SEED
python reviewed_code/trim_reads.py --input SOURCE --outdir OUT --length L --mode fixed --fixed-five-trim 1 --seed SEED
python reviewed_code/trim_reads.py --input SOURCE --outdir OUT --length L --mode fixed --fixed-five-trim 2 --seed SEED
```

Run the same composition/regression/DFT analyses on every output. Verify source and output counts from manifests.

## 4. DFT

```bash
python reviewed_code/dft_mod3.py \
  --input INPUT --length L \
  --window-start 10 --window-end 39 \
  --out-prefix OUTPUT_PREFIX \
  --ymax SHARED_LIMIT
```

Use one common y-axis limit for matched panels. Archive the cycle-frequency CSV, full spectrum CSV, JSON metadata, and PDF/PNG. Verify that the JSON reports the expected exact `q=10`, `f=1/3` bin and that A/T/G/C labels match source columns.

## 5. Library-level aggregation

```bash
python reviewed_code/aggregate_period3_results.py \
  --base-outdir OUTPUT_ROOT \
  --output population/best_length_per_library.csv
```

Publish both the all-eligible-length table and selected table. For inference, either:

- apply the entire maximum-over-length rule to each permuted/simulated null library, or
- prespecify a read length/range and use a summary that is directly comparable among libraries.

Then regenerate metadata plots with explicit paths and preserve merge diagnostics/source data.

## 6. Hard-EM validation and rerun

```bash
python reviewed_code/EM_model.py \
  --input INPUT --input-type auto \
  --out OUTDIR --lengths L1 L2 L3 \
  --dedup none --canonicalize full6 \
  --alpha 0.5 --restarts 10 --max-iter 100 \
  --threads N
```

Before interpreting biological data, benchmark phase-randomized/null inputs, known-state synthetic mixtures, restart stability, held-out likelihood, and external-library replication. Report score gaps as score gaps, not posterior probabilities.

## 7. Coding-frame normalization

```bash
python reviewed_code/frame0_trim.py \
  --bam ALIGNMENTS.bam --gff ANNOTATION.gff3 \
  --output frame0.fasta --min-mapq 25
```

Confirm FASTA/GFF contig naming, sense definition, full containment, GFF phase, and every exclusion count. Analyze normalized and matched unnormalized CDS reads through the same composition pipeline. Report amplitudes, not only near-ceiling R-squared values.

## 8. Fragmentation simulation

```bash
python reviewed_code/simulate_aDNA_reads.py \
  --input REFERENCE.fa --output-dir OUT \
  --num-reads N --read-length L \
  --bias-5prime B5 --bias-3prime B3 \
  --reverse-complement-probability 0.5 \
  --seed SEED
```

Use multiple independent seeds for every genome/bias combination. Analyze every replicate through the same composition pipeline and combine results in a CSV with columns:

`genome,bias_5prime,bias_3prime,replicate,base,r2`

Then run:

```bash
python reviewed_code/plot_simulation_r2.py --input-csv simulation_results.csv --outdir figures
```

Publish the per-base/per-replicate table and all simulator metadata.

## 9. Final consistency check

- Recompute every number in the manuscript from machine-readable source data.
- Link each figure/table to its generating command and manifest.
- Confirm all sample counts, read lengths, windows, seeds, control definitions, reference versions, and error-bar definitions.
- Resolve all declarations and every red author-action box.
- Deposit code/data with a persistent identifier and license before submission.

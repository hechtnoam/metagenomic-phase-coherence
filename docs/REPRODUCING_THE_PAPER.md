# Reproducing the manuscript analyses

This guide describes the analysis workflows and numerical source data deposited
with manuscript release v1.0.0. Run commands from the repository root.

Large primary sequencing files and reference genomes are not redistributed here.
Public sequencing accessions for the cross-library analysis are listed in
`supplementary/Supplementary_Table_S1.tsv` and
`supplementary/Supplementary_Table_S1.xlsx`. The exact AncientMetagenomeDir
environmental-library metadata table used for annotation is archived under
`metadata/`.

## 1. Environment and release checks

```bash
mamba env create -f environment.yml
conda activate period3-metagenomic-reads

python3 tests/run_smoke_tests.py
bash tests/run_cli_integration.sh
```

The tests create temporary output directories that are excluded from version
control.

## 2. Read-position nucleotide profiles and period-3 statistics

The main pipeline computes nucleotide frequencies at each read position within
exact-length read groups and automatically fits the period-3 sine/cosine model
for k=1.

```bash
python3 scripts/period3_library_pipeline.py all INPUT \
  --base-outdir results/libraries \
  --min-read-count 40000 \
  --workers N \
  --p3-start 10 --p3-end 40 \
  --p3-perms 10000 --p3-seed SEED
```

Use `--collapse-exact-sequences` for analyses specified in the manuscript as
being performed on distinct exact sequences; omit it when every retained read is
to contribute. Non-ACGT reads are excluded before exact-length analysis.

For the A2424 read lengths used in Figures 1 and 2:

```bash
python3 scripts/period3_library_pipeline.py all A2424.bam \
  --base-outdir results/libraries \
  --collapse-exact-sequences \
  --lengths 59 61 79 80 \
  --k 1 \
  --workers 8 \
  --p3-start 10 --p3-end 40 \
  --p3-perms 10000 --p3-seed 602758
```

The deposited Figure 1 and Figure 2 source tables are under
`source_data/figure_1/` and `source_data/figure_2/`. Regenerate the final plots
directly from those frozen tables with:

```bash
python3 scripts/plot_manuscript_figures_1_2.py --figure both
```

## 3. Phase redistribution and fixed-trim controls

For an exact input read length L, the randomized control assigns each distinct
input sequence to one of the three retained phase offsets. The fixed controls
apply the corresponding 0/2, 1/1, and 2/0 trims to every read.

```bash
python3 scripts/phase_redistribution.py --input SOURCE --outdir OUT/randomized \
  --length L --mode random --seed SEED
python3 scripts/phase_redistribution.py --input SOURCE --outdir OUT/fixed02 \
  --length L --mode fixed --fixed-five-trim 0
python3 scripts/phase_redistribution.py --input SOURCE --outdir OUT/fixed11 \
  --length L --mode fixed --fixed-five-trim 1
python3 scripts/phase_redistribution.py --input SOURCE --outdir OUT/fixed20 \
  --length L --mode fixed --fixed-five-trim 2
```

Analyze each resulting L-2 read set with `scripts/period3_library_pipeline.py`
using the same period-3 settings. Manuscript source data are stored in
`source_data/figure_3/` and `source_data/phase_redistribution_replicates/`.
Figure 3 plotting is implemented in `scripts/plot_phase_redistribution_figure.py`.

## 4. Frequency-domain representation

The DFT analysis uses read positions 10-39 inclusive, a 30-position window in
which frequency 1/3 is exactly Fourier bin q=10.

```bash
python3 scripts/dft_mod3.py --input INPUT --length L \
  --window-start 10 --window-end 39 --out-prefix OUTPUT_PREFIX
```

Matched panels may additionally use a common `--ymax`. Numerical spectra and
metadata are written alongside the plot. Figure 4 source data are deposited in
`source_data/figure_4/`. The DFT is a frequency-domain representation of the
same period-3 harmonic as the sine/cosine regression, not an independent test.

## 5. Cross-library aggregation

After per-library period-3 analyses have been generated:

```bash
python3 scripts/aggregate_library_results.py \
  --base-outdir results/libraries \
  --output results/population/best_length_per_library.csv
```

For each library the script computes mean R2 across A/T/G/C for every eligible
exact read length and records the length with the largest mean R2. Figure 5 uses
this statistic descriptively; it is not a multiplicity-adjusted library-level
significance test. Deposited source data are under `source_data/figure_5/`.

## 6. Coding-frame normalization

```bash
python3 scripts/frame_normalize_cds.py \
  --bam ALIGNMENTS.bam --gff ANNOTATION.gff3 \
  --output results/coding_frame/frame_normalized.fasta \
  --min-mapq 25
```

The script also writes a matched unnormalized FASTA containing the same accepted
read IDs before phase-dependent trimming.

```bash
python3 scripts/extract_matched_normalized_length.py \
  --normalized results/coding_frame/frame_normalized.fasta \
  --unnormalized results/coding_frame/frame_normalized.matched_unnormalized.fasta \
  --length 75 --outdir results/coding_frame/matched_L75
```

Analyze both matched FASTAs through the same period-3 pipeline. Deposited
numerical results are under `source_data/coding_frame_analysis/`.

## 7. Six-state hard-EM model

```bash
python3 scripts/six_state_em.py \
  --input INPUT --input-type auto \
  --out results/six_state_em \
  --lengths 59,60,61 \
  --dedup sequence --canonicalize full6 \
  --alpha 0.5 --restarts 20 --max-iter 100 --threads N
```

Each restart is retained in the diagnostic outputs. The inferred hard states are
model assignments, not posterior probabilities or directly observed biological
reading frames. Synthetic known-state validation is implemented in
`scripts/validate_six_state_em.py`; Figure 7 plotting is implemented in
`scripts/plot_six_state_em_figure.py`. Source data are in `source_data/figure_7/`.

## 8. Fragmentation simulations

The manuscript grid uses ten independent seeds for each of seven boundary
conditions in both the human and Pseudomonas references. Reference FASTA paths
are supplied through environment variables; no machine-specific paths are stored
in the release.

```bash
HUMAN_REF=/path/to/GCF_000001405.26_GRCh38_genomic.fna \
PSEUDO_REF=/path/to/pseudomonas_reference.fna \
bash scripts/run_fragmentation_grid.sh
```

Optional environment variables are `OUTROOT`, `NUM_READS`, `READ_LENGTH`,
`PERMS`, `P3_WORKERS`, and `PYTHON`.

The grid script runs `scripts/simulate_fragmentation.py`, analyzes every
replicate with `scripts/period3_library_pipeline.py`, combines results with
`scripts/collect_simulation_results.py`, and generates summary plots with
`scripts/plot_simulation_r2.py`. Deposited per-base, per-replicate results are
in `source_data/figure_6/fragmentation_simulation_results.csv`.

## 9. Release-level consistency checks

- Run both repository test suites successfully.
- Confirm that `git diff --check` reports no whitespace errors.
- Confirm that source-data files referenced by plotting scripts are present.
- Confirm the final manuscript title and version in `README.md` and `CITATION.cff`.
- Verify that the Git tag, GitHub release, and Zenodo archive refer to the same
  release snapshot.

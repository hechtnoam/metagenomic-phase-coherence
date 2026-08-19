# Period-3 periodicity in metagenomic sequencing reads

Code supporting the manuscript **“Widespread phase-coherent three-base periodicity in metagenomic sequencing reads.”**

This repository contains the publication-facing analysis code for read-coordinate nucleotide profiles, period-3 regression and permutation tests, DFT analysis, phase-redistribution controls, the six-state hard-EM model, coding-frame normalization, fragmentation simulations, and cross-library summaries. Exploratory scripts that are not part of the manuscript workflow are intentionally excluded.

## Installation

The recommended environment is Conda/Mamba:

```bash
mamba env create -f environment.yml
conda activate period3-metagenomic-reads
```

## Core scripts

| Analysis | Script |
|---|---|
| Cycle nucleotide profiles and main period-3 pipeline | `scripts/period3_library_pipeline.py` |
| Sinusoidal period-3 fits and permutation tests | `scripts/period3_statistics.py` |
| DFT spectrum and f=1/3 summaries | `scripts/dft_mod3.py` |
| Phase-randomization and fixed-trim controls | `scripts/phase_redistribution.py` |
| Six-state hard-EM model | `scripts/six_state_em.py` |
| CDS/codon-frame normalization | `scripts/frame_normalize_cds.py` |
| Purine-associated fragmentation simulation | `scripts/simulate_fragmentation.py` |
| Cross-library aggregation | `scripts/aggregate_library_results.py` |
| Project/metadata summary plots | `scripts/plot_project_summary.py` |
| Simulation R2 summary plots | `scripts/plot_simulation_r2.py` |

See `docs/REPRODUCING_THE_PAPER.md` for the rerun plan and `docs/FIGURE_SCRIPT_MAP.md` for the manuscript-to-code crosswalk.

## Quick example

```bash
python scripts/period3_library_pipeline.py all reads.fastq.gz \
  --base-outdir results \
  --min-read-count 40000 \
  --p3-perms 10000 \
  --p3-seed 12345
```

The default period-3 testing window is cycles 10–40 inclusive. The DFT script uses a 30-cycle window (10–39 by default) so that f=1/3 lies exactly on a Fourier bin.

## Data

Large primary sequencing files and reference genomes are not stored in Git. Put accession numbers, reference versions, and checksums in `docs/DATASETS.md`. Numerical data underlying final manuscript figures should be placed under `source_data/` once the final reruns are complete.

## Tests

```bash
python tests/run_smoke_tests.py
bash tests/run_cli_integration.sh
```

BAM/CRAM/GFF workflows additionally require representative aligned inputs for full integration testing.

## Reproducibility note

The scripts here are the audited publication-facing versions. Historical exploratory and superseded scripts are intentionally not included. Final manuscript figures should be regenerated from frozen inputs with the commands recorded in `docs/REPRODUCING_THE_PAPER.md` before creating the submission release.

## Citation

Citation metadata is provided in `CITATION.cff`. Update the DOI and publication metadata after acceptance/publication.

## License

MIT License. See `LICENSE`.

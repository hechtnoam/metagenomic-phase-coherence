# Period-3 structure in complex DNA sequencing libraries

Code and source data supporting the manuscript
**“Ubiquitous 3-Base Periodicity in Complex DNA Mixtures and a Fast, Alignment-Free Algorithm to Uncover Its Source.”**

This repository contains the analysis and validation workflows used to
characterize read-coordinate three-base periodicity in ancient and modern
metagenomic sequencing libraries. Analyses include exact-length nucleotide
composition profiling, period-3 regression and permutation testing, discrete
Fourier analysis, phase-redistribution controls, coding-frame normalization,
purine-associated fragmentation simulations, and a six-state hard-EM model.

## Installation

The recommended environment is Conda/Mamba:

```bash
mamba env create -f environment.yml
conda activate period3-metagenomic-reads
```

Alternatively, Python dependencies are listed in `requirements.txt`.

## Core scripts

| Analysis | Script |
|---|---|
| Cycle nucleotide profiles and period-3 pipeline | `scripts/period3_library_pipeline.py` |
| Sinusoidal period-3 fits and permutation tests | `scripts/period3_statistics.py` |
| DFT spectrum and frequency-1/3 summaries | `scripts/dft_mod3.py` |
| Phase redistribution and fixed-trim controls | `scripts/phase_redistribution.py` |
| Six-state hard-EM model | `scripts/six_state_em.py` |
| Synthetic validation of the six-state model | `scripts/validate_six_state_em.py` |
| CDS/coding-frame normalization | `scripts/frame_normalize_cds.py` |
| Extraction of matched coding-frame datasets | `scripts/extract_matched_normalized_length.py` |
| Purine-associated fragmentation simulation | `scripts/simulate_fragmentation.py` |
| Fragmentation simulation grid | `scripts/run_fragmentation_grid.sh` |
| Cross-library aggregation | `scripts/aggregate_library_results.py` |
| Project/metadata summary plots | `scripts/plot_project_summary.py` |
| Fragmentation summary figure | `scripts/plot_simulation_r2.py` |
| Six-state EM figure | `scripts/plot_six_state_em_figure.py` |
| Supplementary Table S1 generation | `scripts/generate_supplementary_table_S1.py` |

Additional reproduction information is provided in
`docs/REPRODUCING_THE_PAPER.md` and `docs/FIGURE_SCRIPT_MAP.md`.

## Quick example

```bash
python3 scripts/period3_library_pipeline.py all reads.fastq.gz \
  --base-outdir results \
  --min-read-count 40000 \
  --p3-start 10 \
  --p3-end 40 \
  --p3-perms 10000
```

The primary period-3 regression window is cycles 10–40 inclusive. The DFT
analysis uses cycles 10–39, giving a 30-cycle window in which frequency 1/3
falls exactly on a Fourier bin.

## Data and source data

Large primary sequencing datasets and reference genomes are not distributed
with this repository. Public sequencing accessions and dataset provenance for
the ancient-DNA cross-library analysis are provided in
`supplementary/Supplementary_Table_S1.xlsx`.

The exact AncientMetagenomeDir environmental-library metadata table used for
the cross-library analysis is archived under `metadata/`.

Compact numerical data underlying the manuscript figures and principal
reported analyses are provided under `source_data/`. Large generated
intermediate files, including simulated FASTA files and analysis output
directories, are intentionally excluded from version control and can be
regenerated using the supplied scripts.

## Supplementary material

`Supplementary_Table_S1.xlsx` contains accession-level provenance and selected
period-3 summary values for the ancient-DNA cross-library analysis. A
machine-readable TSV version and matching report are provided in the same
directory.

## Tests

```bash
python3 tests/run_smoke_tests.py
bash tests/run_cli_integration.sh
```

Small representative test inputs are included under `tests/data/`.

## Reproducibility

The release associated with the submitted manuscript is `v1.0.0`. Analysis
parameters and deterministic random seeds are recorded by the relevant
workflows and/or their output metadata.

The exact software release associated with the manuscript is permanently
archived at Zenodo. The DOI for the archived release will be added after the
release is deposited.

## Citation

Citation metadata are provided in `CITATION.cff`.

## License

This software is distributed under the MIT License. See `LICENSE`.
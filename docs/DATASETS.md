# Datasets and reference inputs

This repository does not redistribute large primary sequencing files or reference
genomes. Public accessions, frozen metadata, and manuscript source tables are
provided so that the analyzed inputs can be identified and retrieved.

## Cross-library sequencing datasets

Accession-level provenance for the libraries included in the cross-library
analysis is provided in:

- `supplementary/Supplementary_Table_S1.tsv`
- `supplementary/Supplementary_Table_S1.xlsx`

These tables record the manuscript library identifier, project information,
archive accession, source metadata, and the analysis status used for the
cross-library comparison.

The exact AncientMetagenomeDir environmental-library metadata table used for
annotation is archived as:

- `metadata/ancientmetagenome-environmental_libraries.tsv`

See `metadata/README.md` for the role of this frozen metadata file.

## A2424 analyses

A2424 is the manuscript's principal worked example. The repository contains
machine-readable numerical source data for the exact-length read-position
profiles, phase-redistribution controls, DFT analysis, coding-frame comparison,
and six-state model under `source_data/`. Raw reads and alignment files are not
redistributed.

## Reference genomes and annotations

The fragmentation simulations use the human GRCh38 reference
(`GCF_000001405.26`) and the *Pseudomonas aeruginosa* PAO1 reference
(`AE004091.2`). Reference FASTA paths are supplied at runtime to
`scripts/run_fragmentation_grid.sh` through the `HUMAN_REF` and `PSEUDO_REF`
environment variables.

The coding-frame normalization analysis uses the concatenated microbial/human
reference and the corresponding CDS annotations described in the manuscript.
Those large reference and annotation files are not redistributed here.

## Numerical source data

Numerical values underlying the manuscript figures and matched analyses are
deposited under `source_data/`. See `source_data/README.md` and
`docs/FIGURE_SCRIPT_MAP.md` for the mapping between source tables and manuscript
figures.

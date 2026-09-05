# Test status

The release candidate was tested on 2026-09-05 in the
`period3-metagenomic-reads` environment.

- `python3 tests/run_smoke_tests.py` passes all reported checks, including
  compilation of 17 publication scripts, period-3 statistics, FASTA/FASTQ/BAM
  input handling, DFT, phase redistribution, fragmentation simulation,
  coding-frame normalization logic, hard-EM fitting, aggregation, simulation
  plotting, and BAM integration.
- `bash tests/run_cli_integration.sh` is the end-to-end command-line integration
  test for the released scripts. The release should be tagged only after this
  command completes successfully.
- Generated test workspaces are excluded from version control by `.gitignore`.

Release-level test results should be regenerated after any code change affecting
the analysis pipeline.

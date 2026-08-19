# Test status

The GitHub-ready package was checked on 2026-08-19.

- All 10 publication-facing Python scripts compile.
- `python tests/run_smoke_tests.py` passes all runnable pure-Python tests.
- `bash tests/run_cli_integration.sh` passes for the runnable FASTA/CSV workflows and produces the expected manifests/outputs.
- BAM/CRAM/GFF integration remains environment/data dependent and should be rerun with representative aligned inputs before the submission release.

The test scripts create temporary/generated output directories that are ignored by `.gitignore`.

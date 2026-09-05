# Manuscript analysis-to-code map

| Manuscript analysis | Primary script(s) |
|---|---|
| Read-cycle nucleotide-frequency profiles | `scripts/period3_library_pipeline.py` |
| Period-3 sinusoidal fit, amplitude, R2, permutation testing | `scripts/period3_statistics.py` |
| DFT frequency-domain analysis | `scripts/dft_mod3.py` |
| Random phase redistribution / matched trimming controls | `scripts/phase_redistribution.py` |
| Six-state latent phase/orientation model | `scripts/six_state_em.py` |
| Coding-frame normalization | `scripts/frame_normalize_cds.py` |
| Purine-associated fragmentation simulations | `scripts/simulate_fragmentation.py` |
| Simulation R2 summaries | `scripts/plot_simulation_r2.py` |
| Library-level best-length aggregation | `scripts/aggregate_library_results.py` |
| Project/metadata-level summary plots | `scripts/plot_project_summary.py` |

# Manuscript figure-to-code map

| Figure / analysis | Primary script(s) | Source data |
|---|---|---|
| Figure 1: introductory A2424 read-position profile | `scripts/plot_manuscript_figures_1_2.py` | `source_data/figure_1/` |
| Figure 2: modulo-3 read-length phase structure | `scripts/plot_manuscript_figures_1_2.py` | `source_data/figure_2/` |
| Figure 3: phase redistribution and fixed trimming | `scripts/phase_redistribution.py`, `scripts/plot_phase_redistribution_figure.py` | `source_data/figure_3/`, `source_data/phase_redistribution_replicates/` |
| Figure 4: DFT representation of period-3 coherence | `scripts/dft_mod3.py` | `source_data/figure_4/` |
| Figure 5: cross-library period-3 summary | `scripts/aggregate_library_results.py`, `scripts/plot_project_summary.py` | `source_data/figure_5/` |
| Figure 6: nucleotide-dependent fragmentation simulation | `scripts/simulate_fragmentation.py`, `scripts/run_fragmentation_grid.sh`, `scripts/collect_simulation_results.py`, `scripts/plot_simulation_r2.py` | `source_data/figure_6/` |
| Figure 7: six-state hard-EM model | `scripts/six_state_em.py`, `scripts/validate_six_state_em.py`, `scripts/plot_six_state_em_figure.py` | `source_data/figure_7/` |
| Coding-frame normalization | `scripts/frame_normalize_cds.py`, `scripts/extract_matched_normalized_length.py` | `source_data/coding_frame_analysis/` |

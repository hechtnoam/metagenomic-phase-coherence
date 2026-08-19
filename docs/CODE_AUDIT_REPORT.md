# Code audit report for the Genome Biology manuscript

## Executive assessment

All 20 uploaded Python files were inspected for syntax, data flow, coordinate conventions, statistical definitions, reproducibility, and consistency with the manuscript. Ten scripts form the manuscript's core analytical workflow and were replaced with reviewed versions. The other uploaded scripts are exploratory, display-only, superseded, or insufficiently specified for use as manuscript evidence.

The audit found several defects capable of changing numerical or biological conclusions. Most importantly, the original DFT script swapped the labels of T and C; the original coding-frame normalization ignored the GFF phase field and did not enforce a unique, fully contained, sense CDS alignment; the original simulator sampled contigs per worker rather than candidates globally; and the original hard-EM routine could return a parameter matrix and assignments from different EM states. Consequently, the historical DFT, EM, coding-frame-normalization, and fragmentation-simulation figures must not be used until they are regenerated from frozen inputs with the reviewed code.

The reviewed code passed compilation, pure-function tests, FASTA/FASTQ tests, synthetic statistical tests, and command-line integration tests. No code review can establish that software is literally free of every bug. BAM/CRAM/GFF integration could not be executed in this offline environment because `pysam` and `intervaltree` were unavailable and representative BAM/GFF/reference files were not supplied. The reviewed BAM/GFF workflows therefore still require one production-environment integration run before submission.

## Critical findings in the uploaded code

### 1. DFT nucleotide labels were wrong

`dft_mod3(7).py` imported `BASES` from `EM_model(3).py`. The EM file defines the column order as `A, C, G, T`, while the DFT frequency matrix was populated in the order `A, T, G, C`. The DFT plotting loop then used the imported order to select columns. This exchanges T and C in the spectrum and invalidates the supplied DFT figures.

**Disposition:** corrected in `reviewed_code/dft_mod3.py`. The revised script declares `BASES = ("A", "T", "G", "C")` locally, exports the full numeric spectrum and metadata, uses cycles 10--39 by default, and identifies the exact 30-cycle `q=10` bin at `f=1/3`. All old DFT panels are quarantined in `legacy_figures_not_for_submission/`.

### 2. Phase-redistribution trimming could corrupt FASTQ/BAM records

The original FASTQ code changed `record.seq` before replacing the quality annotation, which can produce a sequence/quality length inconsistency. The original BAM path truncated `query_sequence`, assigned a simple all-match CIGAR, and marked the record unmapped without comprehensively clearing alignment-associated fields. It also trimmed the BAM-stored sequence directly rather than explicitly restoring original sequenced orientation first.

**Disposition:** corrected in `reviewed_code/trim_reads.py`. Biopython record slicing preserves sequence-quality synchronization. BAM records are emitted as valid unmapped records from the original sequenced sequence and qualities. The script supports the randomized 0/2, 1/1, 2/0 mixture and all three matched fixed-trim controls, records the seed and realized offsets, and writes a manifest with hashes and counts.

### 3. Fragmentation simulation did not sample the reference globally

The original simulator assigned one reference contig to each worker and required that worker to produce a fixed number of accepted reads. With a multi-contig reference, contig representation therefore depended on thread count and worker assignment rather than the number of candidate start sites. The original script also emitted only forward-reference fragments.

**Disposition:** corrected in `reviewed_code/simulate_aDNA_reads.py`. Candidate positions are drawn globally in proportion to the number of valid start sites across all contigs, then subjected to the stated boundary acceptance rule. Forward-reference and reverse-complement orientations are sampled explicitly. Metadata records candidate and accepted contig counts, orientations, boundary bases, attempts, input hash, and seed.

### 4. Coding-frame normalization did not implement the stated method

The original `frame0_trim.py` parsed the GFF phase field but discarded it, selected an arbitrary overlap with `next(iter(overlaps))`, accepted partial or antisense overlaps, and performed strand-dependent slicing on `query_sequence` without a robust original-sequence-orientation rule. These choices can change the assigned codon phase and the retained sequence.

**Disposition:** corrected in `reviewed_code/frame0_trim.py`. The reviewed workflow requires a primary, QC-passing, simple M/=/X alignment with sufficient mapping quality, full containment in exactly one sense CDS, and uses the GFF3 phase on both strands. It restores the original sequenced/coding orientation, trims 0, 2, or 1 nucleotide from the biological 5' end for phases 0, 1, or 2, and trims the coding 3' end to retain complete codons. Exclusion counts and run metadata are written.

### 5. Hard-EM outputs could be internally inconsistent

In the original `fit_em`, assignments were obtained in the E-step, a new parameter matrix was then generated in the M-step, and the function could return that new matrix together with assignments computed under the previous matrix. The original defaults also enabled hash-based deduplication and disabled canonicalization, although neither behavior was described in the manuscript.

**Disposition:** corrected in `reviewed_code/EM_model.py`. A run is accepted only after a self-consistent fixed point is reached: the returned assignments are the E-step optimum under the returned matrix and the matrix is the M-step update from those assignments. The objective trace is checked, restart selection is deterministic, exact-sequence deduplication is off by default, ambiguity is excluded by default, and full six-state canonicalization is the default display convention. Assignment gaps are explicitly described as score gaps, not posterior probabilities.

### 6. The project-level statistic involves selection over read lengths

The original aggregator first averaged the four nucleotide-specific R-squared values at each read length and then chose the maximum mean across eligible lengths. The project dot plot therefore represents a maximum-over-length statistic, not a prespecified single-length statistic. Its null distribution changes with the number and range of eligible lengths.

**Disposition:** `reviewed_code/aggregate_period3_results.py` preserves the historical definition for transparency but exports every eligible length, the number and range of lengths searched, counts, tie information, errors, and a manifest. The manuscript now labels the historical plot provisional and requires calibration of the full selection procedure by permutation/simulation or replacement with a prespecified summary.

### 7. Simulation figure values were hard-coded

`plot_r2.py` contains a literal array of human-simulation R-squared values and plots the standard deviation across four bases as if it were an uncertainty estimate. No machine-readable per-read/per-replicate source was supplied, and no corresponding bacterial source-data script was supplied.

**Disposition:** replaced by `reviewed_code/plot_simulation_r2.py`. It requires machine-readable rows containing genome, both bias values, replicate, base, and R-squared. It first averages across the four bases within each independent replicate, then reports the mean and 95% t confidence interval across independent replicates. Old simulation panels are quarantined.

### 8. Default filtering was not fully transparent

The original library pipeline enabled exact-sequence collapsing and skip-existing behavior by default. Those defaults can change effective sample size and can silently reuse outputs from a different code or parameter state. Ambiguous symbols could also affect denominators in ways not clearly reported.

**Disposition:** corrected in `reviewed_code/period3_library_pipeline.py`. Exact-sequence collapse and skip-existing behavior are opt-in. Non-ACGT reads are excluded and counted before length selection. Secondary, supplementary, and QC-failed BAM records are excluded; the SAM duplicate flag is retained unless the user explicitly excludes it. Original sequenced orientation is restored for reverse-strand BAM alignments. Run manifests include parameters, counts, hashes, Python version, and script/module versions.

## Statistical alignment with the manuscript

- **Period-3 regression:** cycles 10--40 inclusive; intercept plus sine and cosine at frequency 1/3; nucleotide-specific R-squared and absolute amplitude are both exported.
- **Permutation testing:** shuffles values within the selected cycle window; deterministic child seeds; 10,000 permutations by default; Holm adjustment across A, T, G, and C within one read-length class. This does not correct the separate search over read lengths.
- **DFT:** cycles 10--39 inclusive; 30 positions; normalized magnitude `|F(q)|/30`; exact period-3 bin `q=10`; zero-frequency term omitted from the displayed positive-frequency spectrum. The DFT is a frequency-domain representation of the same harmonic, not an independent statistical test of the sinusoidal regression.
- **Phase redistribution:** random 5' offset in {0,1,2} plus complementary 3' trimming for total trim 2; all outputs have length L-2; fixed 0/2, 1/1, and 2/0 controls are required on the same filtered source read set.
- **Library summary:** historical statistic is the maximum, across eligible lengths, of the mean of four nucleotide R-squared values.
- **Hard EM:** deterministic six-state hard assignments, not posterior probabilities or validated biological reading-frame calls.
- **Frame normalization:** fully contained, unique, sense CDS alignments with GFF phase and complete-codon output.
- **Fragmentation simulation:** candidate-site rejection sampling with purine acceptance parameter b in [0.5,1]; A/G accepted with probability 1 and C/T with probability `(1-b)/b`. This parameter is an acceptance rule, not necessarily the final purine fraction.

## Script-by-script disposition

| Uploaded file | Relevance | Disposition |
|---|---|---|
| `period3_library_pipeline(5).py` | Core | Replaced by audited pipeline |
| `period3_statistics(2).py` | Core | Replaced; validation, amplitude, Holm adjustment, deterministic seeds added |
| `dft_mod3(7).py` | Core | Replaced; critical T/C labeling bug fixed |
| `trim_reads(2).py` | Core | Replaced; FASTQ/BAM integrity and matched controls fixed |
| `EM_model(3).py` | Core | Replaced; fixed-point consistency and defaults corrected |
| `frame0_trim.py` | Core | Replaced; GFF phase, containment, sense, orientation, and filters corrected |
| `simulate_aDNA_reads.py` | Core | Replaced; global multi-contig and strand-aware sampling added |
| `aggregate_period3_results(1).py` | Core | Replaced; full provenance and selection metadata added |
| `generate_ancient_metadata_histograms(1).py` | Core figure support | Replaced; portable inputs, validated joins, source data, deterministic jitter, manifest |
| `plot_r2.py` | Core figure support | Superseded by machine-readable replicate plotter |
| `zoom_k1.py` | Display only | Not part of manuscript inference; retained only as original upload |
| `periodicity_stats.py` | Exploratory | Not used in manuscript; hard-coded workflow and denominator/reporting concerns |
| `plot_weighted_period3(1).py` | Exploratory | Not used; coordinate/CIGAR and weighting assumptions require a separate validated analysis plan |
| `plot_weighted_noncoding_control.py` | Exploratory | Not used; one-coordinate CDS exclusion and overlapping-3-mer background do not establish a noncoding control |
| `plot_weighted_codon_period3.py` | Exploratory | Not used; overlapping 3-mers are not GFF-phase codons and exceptions are broadly suppressed |
| `plot_weighted_cds_length_filter.py` | Exploratory | Not used; similar CDS/coordinate limitations |
| `plot_minus1_groups(1).py` | Exploratory | Not used; boundary coordinates assume simple alignment geometry |
| `plot_hierarchical_G.py` | Exploratory | Not used; grouping analysis lacks a frozen inferential specification |
| `analyze_adna.py` | Exploratory damage plot | Not used in current manuscript; needs explicit CIGAR/soft-clip/indel safeguards and provenance |
| `analyze_adna_compare.py` | Exploratory damage plot | Not used in current manuscript; same safeguards required |
| `dft_mod3_pipeline_version(1).py` | None | Placeholder text rather than executable Python; excluded from workflow |

## Tests performed

See `TEST_RESULTS.txt` and `CLI_INTEGRATION_RESULTS.txt` for the exact results. The tests cover:

- compilation of all reviewed modules;
- pure and constant period-3 signals, amplitude/R-squared, permutation p-values, and Holm adjustment;
- CRLF FASTQ decoding, ambiguity exclusion, cycle fractions, and BAM-orientation helper logic;
- the exact DFT `q=10` bin, explicit A/T/G/C mapping, and spectral R-squared for a pure signal;
- FASTQ sequence/quality integrity and fixed-trim controls;
- deterministic multi-contig simulation, both orientations, and `b=1` purine-only boundaries;
- plus/minus CDS coordinate and GFF phase logic;
- hard-EM fixed-point consistency, monotone objective, restart determinism, and canonicalization;
- complete-row library aggregation, metadata joins, unmatched reporting, and source-data export;
- machine-readable simulation summaries and confidence intervals;
- end-to-end CLI execution for every reviewed workflow that does not require BAM/GFF dependencies.

## What still must happen before submission

1. Install the pinned environment including `pysam` and `intervaltree`, and run the BAM/GFF integration suite on representative data.
2. Freeze exact input accession lists, reference FASTA/GFF versions, checksums, software versions, and command lines.
3. Regenerate every empirical and mechanistic output from the reviewed code. Do not reuse the quarantined DFT, EM, frame-normalization, or simulation figures.
4. For phase redistribution, use the identical filtered source set, all three fixed-trim controls, multiple prespecified random seeds, and report amplitude, R-squared, adjusted p-values, and DFT outputs.
5. Calibrate the maximum-over-length library statistic under a null that repeats the entire length-selection procedure, or replace it with a prespecified statistic.
6. Validate the hard-EM model on null simulations, phase-randomized reads, known-state synthetic data, restart stability, held-out likelihood, and external libraries.
7. Run coding-frame normalization with documented reference/GFF accessions and report every exclusion category plus a matched unnormalized control.
8. Run independent simulation replicates for every reference/bias condition and publish machine-readable per-base/per-replicate data.
9. Publish code, source data, manifests, accession table, and an archived release with a persistent identifier and license.
10. Replace every red author-action box in the manuscript; do not merely hide the boxes.

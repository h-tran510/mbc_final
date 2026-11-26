# deployed_scripts_v1

This folder contains the runnable TCGA-BRCA methylation and survival
pipeline. It brings together data preparation and the promoter analyses,
renumbered from
`scripts/` and `scripts/diagnostics/`, dropping superseded drafts,
abandoned directions, and pure debugging scripts. At consolidation, nothing
changed the analysis logic apart from the one bug fix under "Fixes applied"
below. Later pull requests did change it (for example PAM50 in the 08 baseline
and the Freedman-Lane null in 11, PR #10, and the one-aliquot rule in 02, PR
#23); `git log` and `docs/research_journal.md` record each change.

Exploratory/unvalidated follow-up work (e.g. the full-array variance-filter
bypass screen) intentionally stays out of this folder and lives under
`scripts/diagnostics/` instead — see the note at the bottom of "What's
excluded and why."

## How to run

Start with the [root README](../README.md#start-here) for setup and run
commands. Run them from the repository root.

The sections below record how the deployed pipeline was assembled and which
experiments stayed in the development scripts.

## What's included and why

| # | File | Originally | Purpose |
|---|------|-----------|---------|
| 00 | `00_utils.R` | `scripts/00_utils.R` | Shared helper functions |
| 01 | `01_download_tcga.R` | `scripts/01_download_tcga.R` | Download RNA-seq counts |
| 02 | `02_load_and_prepare_data.R` | `scripts/02_load_and_prepare_data.R` | Merge with clinical data |
| 03 | `03_normalization_and_qc.R` | `scripts/03_normalization_and_qc.R` | DESeq2 normalization, VST |
| 04 | `04_survival_analysis.R` | `scripts/04_survival_analysis.R` | Univariate Cox screen (RNA-seq) |
| 05 | `05_multivariate_cox.R` | `scripts/05_multivariate_cox.R` | LASSO-Cox risk score (RNA-seq) |
| 06 | `06_download_methylation.R` | `scripts/06_download_methylation.R` | Download HM450 methylation |
| 07 | `07_methylation_prepare.R` | `scripts/07_methylation_prepare.R` | QC, M-value transform |
| 08 | `08_methylation_survival.R` | `scripts/08_methylation_survival_v3.R` | Production methylation risk score |
| 09 | `09_pc1_biological_characterization.R` | `scripts/diagnostics/13_...R` | What is PC1, biologically? (B0-B7) |
| 10 | `10_subtype_purity_stratified_screen.R` | `scripts/diagnostics/14_...R` | Per-promoter naive vs. stratified screen |
| 11 | `11_stratified_screen_validation.R` | `scripts/diagnostics/15_...R` | Permutation + Her2-exclusion validation |
| 12 | `12_cbx_family_stratified_screen.R` | `scripts/diagnostics/20_...R` | CBX-family signature test (variance-filter bypass) |
| 13 | `13_cbx5_figures.R` | new | Forest plot, KM-by-subtype, permutation histogram for CBX5 |
| 14 | `14_sign_bias_diagnostics.R` | new | Sign bias vs PC1, sign-fraction null, Freedman-Lane permutation |
| 15 | `15_external_validation.R` | new | Pre-registered replication of the stratified-screen survivors in GSE141441 |

## What's excluded and why

- **`scripts/08_methylation_survival_v1.R`, `_v2.R`** — superseded drafts of
  what is now script 08. Kept in `scripts/` for history, not needed to
  reproduce the final result.
- **`scripts/diagnostics/08_methylation_survival_alignment.R`,
  `_correct_id.R`, `_direction.R`, `_elastic_test_1.R`, `_id_test.R`,
  `_signal.R`** — earlier debugging sessions (barcode alignment, sign
  conventions, elastic-net experiments) that fed into the current script 07
  and 08 but aren't themselves part of the final pipeline.
- **`scripts/diagnostics/16_epifactors_candidate_panel.R`** — an explored
  but abandoned direction (a broad EpiFactors-derived candidate panel). The
  project pivoted to the CBX-family test (script 12) instead. Kept in
  `diagnostics/` as a record of what was tried.
- **`scripts/diagnostics/17_cbx_family_panel_test.R`** — an earlier version
  of the CBX-family test that ran on `meth_m_matrix.rds` directly, before
  the variance-filter coverage gap was diagnosed. Coverage was incomplete
  (5 of 8 genes missing) at that point, so its combined-score results never
  actually ran. Fully superseded by script 12, which rebuilds full coverage
  from the raw data. Kept in `diagnostics/` for provenance.
- **`scripts/diagnostics/18_check_cbx_coverage.R`,
  `19_check_cbx_qc_dropout.R`** — the two diagnostic scripts that traced the
  CBX coverage gap to `07_methylation_prepare.R`'s variance filter (not
  array design, not cross-reactivity). Their findings are summarized in
  the root README's "Sample and probe choices" section; the scripts themselves are investigative
  breadcrumbs, not part of the reproducible pipeline. Kept in
  `diagnostics/` for provenance.
- **`scripts/diagnostics/21_full_variance_filter_bypass_screen.R`** —
  extends script 12's per-gene variance-filter bypass to the whole array
  (not just CBX), to check whether other promoters were silently excluded
  from the original screen the way CBX1/CBX3/etc. were. This is a large,
  unvalidated, exploratory re-screen (permutation-tested only by the exploratory
  `scripts/diagnostics/22_*`, much
  bigger multiple-testing burden) — deliberately kept out of this folder
  until any candidate it finds is validated with new checks of the kind the
  stratified-screen hits got: 11's permutation null and Her2 exclusion, and
  15's pre-registered external replication (11 is written for 07's filtered
  matrix, 15 for two named genes).

## Fixes applied during consolidation

**`rowMeans` vs `colMeans` bug in combined-score construction** (script 12,
and retroactively in `scripts/diagnostics/17_...R`): building a per-patient
combined score from several genes' promoter values requires averaging
*across genes* for each patient (`rowMeans` on a patients-by-genes matrix).
The original code used `colMeans`, which instead averages *across patients*
for each gene, producing a length-3 (or length-5) vector instead of one
score per patient. This silently made every downstream `coxph()` call fail
inside a `try()` wrapper, returning `NA` rather than an error. This is the
same class of bug caught once already in script 09's B1 proliferation-score
section (there as `colMeans(scale(t(...)))`, same fix). All three sites now
have a `stopifnot(length(...) == ncol(...))` guard to catch a recurrence
immediately instead of silently.

## Known limitations

- **Purity gate in script 08**: purity enters 08's models only when more than
  90% of samples have a CPE estimate (`purity_available`); below that it is
  dropped from every 08 model. Scripts 10 to 14 instead drop patients who lack
  purity (complete cases).
- **The genome-wide variance filter in script 07, step 6(d)**: applied
  before any gene-level analysis, this keeps only the top 10% most variable
  probes array-wide. It is why some genes (confirmed for the full CBX
  family) have thin or zero promoter coverage in `meth_m_matrix.rds`
  despite having real probes on the array. This is expected, documented
  pipeline behavior — script 12 shows the pattern for bypassing it on a
  small, pre-specified panel when needed — but it means any future
  candidate-gene follow-up should check coverage the same way scripts 18/19
  did before assuming a gene is absent or uninformative.

## Result summary

The measured results, each with the run it came from, are in
`docs/research_journal.md`, the source of record; the root `README.md` carries
a short summary that is updated with it.

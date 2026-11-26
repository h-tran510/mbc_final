# Code review — `deployed_scripts_v1/`, scripts 08–13

This review records the original consolidation. For current run instructions,
see the [root README](../README.md#running).

> **STATUS: fixes applied.** Every finding below (C1–C2, H1–H4, M1–M9,
> L1–L3, S1) has been implemented in the scripts. The findings are kept
> verbatim as the rationale for each change. Two things to know before you
> rerun anything:
>
> 1. **One deliberate behaviour change.** `build_promoter_gene_mat()` now
>    drops promoter probes whose `UCSC_RefGene_Name` is empty or `NA`. These
>    previously aggregated into a spurious `""` pseudo-gene that was screened
>    like a real gene and counted in the multiple-testing burden. Expect
>    `nrow(gene_mat)` to fall by 1 and BH-adjusted q-values to shift very
>    slightly. **Verify before trusting new numbers:** rerun script 10 and
>    diff `dim(gene_mat)` and `gene_mat["CBX5", ]` / `gene_mat["CBX2", ]`
>    against your current saved results. The per-gene values themselves
>    should be bit-identical; only the row count and the FDR denominator
>    should move.
> 2. **Run directories are now timestamped** (`%Y-%m-%d_%H%M%S_tag` rather
>    than `%Y-%m-%d_tag`), so same-day reruns no longer overwrite each other.
>    The "find the latest run" helpers in scripts 13/22/23 sort
>    lexicographically and are unaffected. Existing run directories keep
>    their old names and still sort correctly relative to new ones.
>
> S2 (per-patient identifiers written under `02_analysis/`) is a
> `.gitignore` question, not a code change — confirm `02_analysis/` is
> ignored. `05_multivariate_cox.R` also carried the M7 `1 - pchisq(df = 1)`
> bug and was fixed for consistency even though it sits outside the reviewed
> range.

Reviewed: `08_methylation_survival.R`, `09_pc1_biological_characterization.R`,
`10_subtype_purity_stratified_screen.R`, `11_stratified_screen_validation.R`,
`12_cbx_family_stratified_screen.R`, `13_cbx5_figures.R`, plus `00_utils.R`
where shared helpers are implicated.

Severity key: **C** = correctness risk that can silently produce wrong
results; **H** = high (crashes, lost compute, large resource waste); **M** =
medium; **L** = low / stylistic.

---

## C1. Annotation cache is unkeyed and can silently misalign every probe

**Files/lines:** `08:140–148`, `09:78–86`, `10:60–66`, `11:67–73`, `13:50–56`
(all call the same cache), plus `00_utils.R:40–49` (`cache_or_compute`).

```r
ann_cache_path <- "01_processed_data/ann_matched_cache.rds"
ann_matched <- cache_or_compute(ann_cache_path, function() {
  ann_df <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann_probe_idx <- match(rownames(m_mat), rownames(ann_df))
  ann_df[ann_probe_idx, ]
})
is_promoter <- grepl("TSS200|TSS1500", ann_matched$UCSC_RefGene_Group)
m_mat_promoters <- m_mat[is_promoter, ]
```

**Why it's an issue.** `cache_or_compute()` returns the cached object purely
on file existence — it never checks that the cache was built from the *same*
`m_mat`. `ann_matched` is a positional, row-for-row companion to `m_mat`;
its correctness depends entirely on `rownames(m_mat)` being unchanged since
the cache was written. If `07_methylation_prepare.R` is ever rerun with a
different filter (a different variance quantile, a different missingness
threshold, batch correction toggled), `m_mat` changes but the stale cache is
served anyway. `is_promoter` is then a logical vector of the wrong length or
— worse, if lengths coincidentally match — the *right length but wrong
probes*, so every promoter label is silently attached to the wrong CpG. No
error is raised. Every downstream result in 08–13 would be garbage that
looks completely normal.

This is not hypothetical for this project specifically: script 21's
whole-array rebuild produces exactly such an alternative `m_mat`, and the
CBX coverage investigation involved reasoning about which probes are in
`m_mat` at all.

**Fix (minimal, add after every `cache_or_compute` call):**

```r
ann_matched <- cache_or_compute(ann_cache_path, function() { ... })

stopifnot(
  "Cached annotation is stale: row count does not match m_mat. Delete
   01_processed_data/ann_matched_cache.rds and rerun." =
    nrow(ann_matched) == nrow(m_mat),
  "Cached annotation is stale: probe IDs do not match m_mat. Delete
   01_processed_data/ann_matched_cache.rds and rerun." =
    isTRUE(all.equal(rownames(ann_matched), rownames(m_mat)))
)
```

**Fix (structural, preferred — key the cache on its input):**

```r
# in 00_utils.R
cache_or_compute_keyed <- function(path_prefix, key, compute_fn, force = FALSE) {
  key_hash <- substr(digest::digest(key), 1, 10)
  path <- paste0(path_prefix, "_", key_hash, ".rds")
  cache_or_compute(path, compute_fn, force = force)
}

# call site
ann_matched <- cache_or_compute_keyed(
  "01_processed_data/ann_matched_cache",
  key = rownames(m_mat),
  compute_fn = function() { ... }
)
```

**Trade-offs.** The `stopifnot` version is zero-dependency and catches the
bug loudly, but forces a manual cache delete when inputs legitimately
change. The keyed version adds a `digest` dependency and leaves multiple
cache files on disk (a few MB each), but makes stale cache structurally
impossible and lets several `m_mat` variants coexist — which this project
now actually needs. I'd take the keyed version.

---

## C2. `rowsum()` / `table()` alignment is an unchecked implicit invariant

**Files/lines:** `08:151–156`, `09:88–93`, `10:68–72`, `11:75–79`,
`12:230–232` and `12:291`, `13:59–62`.

```r
gene_sums   <- rowsum(m_mat_promoters, group = gene_labels_clean)
gene_counts <- table(gene_labels_clean)
gene_mat    <- gene_sums / as.vector(gene_counts)   # positional division
```

**Why it's an issue.** This divides a matrix by a bare numeric vector,
relying on `rowsum()` and `table()` producing groups in *identical* order.
Both sort by group value today, so it works — but the code never asserts it.
`as.vector()` deliberately strips the names that would otherwise make a
mismatch visible. Any change that introduces `NA` group labels, a factor
with unused levels, or a locale-dependent collation difference silently
divides gene A's summed M-values by gene B's probe count. Result: wrong
promoter means, no error, plausible-looking output.

`12:291` compounds this by re-indexing `gene_counts[recovered_genes]` — that
one *is* name-safe, which makes the inconsistency in the same file more
confusing.

**Fix:**

```r
gene_sums   <- rowsum(m_mat_promoters, group = gene_labels_clean)
gene_counts <- table(gene_labels_clean)

stopifnot(
  "rowsum/table group ordering diverged -- promoter means would be misaligned" =
    identical(rownames(gene_sums), names(gene_counts))
)
gene_mat <- gene_sums / as.numeric(gene_counts)
```

Better still, factor the whole block into `00_utils.R` (it is copy-pasted
verbatim across six scripts) so the guard exists in exactly one place:

```r
# 00_utils.R
build_promoter_gene_mat <- function(m_mat, ann_matched) {
  stopifnot(isTRUE(all.equal(rownames(ann_matched), rownames(m_mat))))
  is_promoter <- grepl("TSS200|TSS1500", ann_matched$UCSC_RefGene_Group)
  labels <- sapply(strsplit(ann_matched$UCSC_RefGene_Name[is_promoter], ";"),
                   function(x) x[1])
  keep <- !is.na(labels) & nzchar(labels)
  sums   <- rowsum(m_mat[is_promoter, ][keep, ], group = labels[keep])
  counts <- table(labels[keep])
  stopifnot(identical(rownames(sums), names(counts)))
  sums / as.numeric(counts)
}
```

**Trade-offs.** Centralising is the right call for a six-way duplication and
also fixes C1 in one place — but it changes every script's data-prep block
at once, so it needs one careful verification run (compare `dim(gene_mat)`
and a few `gene_mat["CBX5", ]` values before/after) before you trust
downstream numbers.

Note the extra `keep <- !is.na(labels) & nzchar(labels)` in that helper:
the current code does not filter empty or `NA` gene labels, so probes
annotated with an empty `UCSC_RefGene_Name` become a spurious `""` "gene"
row in `gene_mat` and get screened like any other. Worth checking whether
`"" %in% rownames(gene_mat)` in your current results.

---

## H1. Clusters are not released on error — leaked workers holding GBs

**Files/lines:** `08:238–240` / `08:290`; `11:174` / `11:197`;
`13:277–279` / `13:290`.

```r
cl <- makeCluster(n_cores)
registerDoParallel(cl)
results <- foreach(...) %dopar% { ... }   # if this throws, stopCluster never runs
stopCluster(cl)
```

**Why it's an issue.** Any error inside the `foreach` body — a singular Cox
fit, an out-of-memory worker, a user interrupt — aborts before
`stopCluster()`. The forked/socket workers survive as orphan R processes,
each still holding its copy of the exported globals (see H2: that is roughly
300 MB–1 GB *each* here). Repeated failed runs in one session will exhaust
RAM and the failure mode looks like an unrelated memory problem later.

**Fix:**

```r
cl <- makeCluster(n_cores)
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
registerDoParallel(cl)
...
```

At top level `on.exit()` is attached to the script's evaluation frame, so
for `source()`-driven use wrap the parallel section in a function or use an
explicit handler:

```r
run_permutations <- function(gene_mat_c, surv_c, meth_c, N_PERM, n_cores) {
  cl <- makeCluster(n_cores)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  registerDoParallel(cl)
  foreach(r = seq_len(N_PERM), .packages = "survival", .combine = rbind) %dopar% {
    ...
  }
}
perm_results <- run_permutations(gene_mat_c, surv_c, meth_c, N_PERM, n_cores)
```

**Trade-offs.** Wrapping in a function is the robust form and also fixes H2
for free (only the arguments get exported, not the whole global
environment). Cost: the permutation block stops being a flat, readable
top-to-bottom script, which is a real readability loss in an analysis
codebase where linear narrative matters.

---

## H2. `foreach` exports the entire global environment, including `m_mat`

**Files/lines:** `08:244–247`, `11:170–174`, `13:280–282`.

**Why it's an issue.** `%dopar%` with doParallel exports every symbol the
loop body appears to reference from the enclosing environment, and in
practice for socket clusters this pulls a large slice of the global env into
*each* worker. At this point in all three scripts, `m_mat` (~48,000 × 755
doubles ≈ 290 MB) and `m_mat_promoters` are still resident, alongside
`gene_mat` (~25 MB). With `detectCores() - 1` = 7 workers that is on the
order of 2 GB of pure duplication for data the loop never touches — the loop
only needs `gene_mat_c`, `surv_c`, `meth_c`, and two functions.

This is very likely a contributor to the 27-minute runtime of script 11's
100-permutation run: serialisation of hundreds of MB per worker at cluster
startup, plus memory pressure.

**Fix (both halves):**

```r
# 1. drop what is no longer needed before spawning workers
rm(m_mat, m_mat_promoters, gene_sums, gene_labels, gene_labels_clean)
gc()

# 2. export explicitly instead of implicitly
perm_results <- foreach(
  r         = seq_len(N_PERM),
  .packages = "survival",
  .export   = c("gene_mat_c", "surv_c", "meth_c",
                "run_full_screen", "stratified_screen_one"),
  .noexport = ls(envir = .GlobalEnv)[!ls(envir = .GlobalEnv) %in%
                c("gene_mat_c", "surv_c", "meth_c",
                  "run_full_screen", "stratified_screen_one")],
  .combine  = rbind
) %dopar% { ... }
```

**Trade-offs.** Explicit `.export` is brittle in the opposite direction — if
you later reference a new object in the loop body and forget to add it,
you get a `could not find function/object` failure inside a worker, whose
error message is much less clear than a normal one. The function-wrapping
approach from H1 avoids both problems and is what I'd recommend instead of
hand-maintaining `.export`.

---

## H3. `write_manifest()` can abort at the very end of a multi-hour run

**File/lines:** `00_utils.R:23–37`, called at `08:443`, `09:653`, `10:271`,
`11:292`, `12:325`.

```r
git_hash <- system("git rev-parse --short HEAD", intern = TRUE)
```

**Why it's an issue.** If the working directory is not a git repo, `git` is
not on `PATH`, or the repo has no commits yet, `system(..., intern = TRUE)`
raises an error (or returns `character(0)` with a warning, which then makes
`paste0` produce a zero-length vector and `writeLines` write a truncated
manifest). This runs *after* all analysis — in script 11 that is after a
27-minute permutation run, in script 08 after a 10×4-fold nested CV. Losing
the manifest is minor; aborting the script before `saveRDS`/`write.csv` of
the actual results is not.

Note the ordering in `08`: `write_manifest()` at line 443 executes *before*
`saveRDS(cv_fit, ...)` at 445 and all `write.csv` calls at 448–465. A git
failure there discards every output of the run.

**Fix:**

```r
write_manifest <- function(run, params) {
  git_hash <- tryCatch(
    system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE),
    error   = function(e) "NOT_A_GIT_REPO",
    warning = function(w) "NOT_A_GIT_REPO"
  )
  if (length(git_hash) == 0) git_hash <- "NOT_A_GIT_REPO"
  ...
}
```

and in `08`, move `write_manifest()` to after the result-writing block.

**Trade-offs.** Swallowing the error means a run can complete with
`NOT_A_GIT_REPO` provenance, which is exactly the provenance you wanted the
manifest for. Mitigation: emit a loud `warning()` in that branch so it shows
up in the console, rather than failing silently.

---

## H4. Manifest parsing in script 13 can silently produce a wrong figure

**File/lines:** `13:243–254`, consumed at `13:295–306`.

```r
obs_n_hits <- as.integer(sub(".*= ", "", grep("observed_hits_fdr10", manifest_lines, value = TRUE)))
perm_p_hitcount <- as.numeric(sub(".*= ", "", grep("permutation_p_hitcount", manifest_lines, value = TRUE)))
```

**Why it's an issue.** Three failure modes, none of which raise an error:

1. If the located run's manifest predates those keys (an older script
   version), `grep(...)` returns `character(0)`, `as.integer(character(0))`
   is `integer(0)`, and downstream `abline(v = integer(0))` draws nothing
   while `sprintf("Observed = %d", integer(0))` yields `character(0)` — the
   legend and the observed-value marker vanish from the PDF with no message.
   You get a publication figure that is missing its single most important
   element.
2. `grep()` is unanchored, so a key like `observed_hits_fdr10_secondary`
   would also match and silently produce a length-2 vector.
3. `sub(".*= ", ...)` is greedy on `.*`; it happens to work for
   `key = value` lines but breaks on any value containing `= `.

**Fix — parse defensively and assert:**

```r
read_manifest_value <- function(lines, key) {
  hit <- grep(paste0("^", key, " = "), lines, value = TRUE)
  if (length(hit) != 1L) {
    stop(sprintf("Expected exactly one '%s' line in manifest, found %d.",
                 key, length(hit)))
  }
  sub(paste0("^", key, " = "), "", hit)
}

obs_n_hits      <- as.integer(read_manifest_value(manifest_lines, "observed_hits_fdr10"))
perm_p_hitcount <- as.numeric(read_manifest_value(manifest_lines, "permutation_p_hitcount"))
stopifnot(is.finite(obs_n_hits), is.finite(perm_p_hitcount))
```

**Trade-offs.** This converts a silent-wrong-figure into a hard stop, which
is correct for a figure destined for a paper. The deeper fix is to stop
round-tripping numbers through a human-readable text manifest at all — have
script 11 `saveRDS(list(obs_n_hits = ..., perm_p = ...), file.path(paths$models, "perm_summary.rds"))`
and have 13 read that. Text manifests are for humans; RDS is for code.

---

## M1. `purity_available` all-or-nothing gate (already documented)

**File/lines:** `08:95–96`, used at `08:227–232, 267–269, 312–316, 336–342`.

```r
purity_available <- !any(is.na(meth_clinical$purity))
```

**Why it's an issue.** Six patients missing a CPE estimate disable purity
adjustment for all 755. Since purity is a documented confounder in this
project (script 09's B3/B7), silently dropping it changes what the
production model is actually adjusting for. The header comment
acknowledges this; flagging it here because it is the highest-impact
*statistical* defect in 08 and the fix is small.

**Fix:**

```r
purity_missing <- sum(is.na(meth_clinical$purity))
if (purity_missing > 0) {
  warning(sprintf(
    "%d of %d samples lack a purity estimate; those rows will be dropped from
     purity-adjusted models by complete-case deletion (coxph default).",
    purity_missing, nrow(meth_clinical)))
}
purity_available <- mean(!is.na(meth_clinical$purity)) > 0.90
```

and let `coxph()`/`cv.glmnet()` handle the remainder — noting that
`cv.glmnet` does *not* do complete-case deletion, so for the LASSO path
you must either subset explicitly or median-impute with a missingness
indicator column.

**Trade-offs.** Complete-case deletion changes the analysis cohort between
models (see M2). Median-imputation-plus-indicator keeps n fixed at the cost
of a mild bias toward the null for those six rows. With 6/755 either is
defensible; what is not defensible is silently dropping the covariate
entirely.

---

## M2. Cross-model comparisons are made on different complete-case subsets

**File/lines:** `09:303–312` (B2, subtype-complete rows) vs `09:364–367`
(B3, purity-complete rows) vs `09:498–507` (B5 c-index comparison).

**Why it's an issue.** `m_plain`/`m_strat` are fit on `meth_clinical[sub_ok, ]`;
`m_pur` is fit on the full frame with `coxph()` dropping purity-missing rows;
B5's three models drop rows on whichever of `pc1`/`cimp` is missing. The HRs
and C-indices logged from these are then read side by side in the results
table as if commensurable. They are not — they describe different patient
sets. C-index in particular is not comparable across differing n.

**Fix:** define one analysis cohort up front, as scripts 10–13 already do,
and fit every model in the section on it:

```r
b_cohort <- !is.na(meth_clinical$subtype) & !is.na(meth_clinical$purity) &
            !is.na(meth_clinical$age_years) & !is.na(meth_clinical$stage)
dat <- meth_clinical[b_cohort, ]
m_plain <- coxph(Surv(time, event) ~ scale(pc1) + age_years + stage, data = dat)
m_strat <- coxph(Surv(time, event) ~ scale(pc1) + age_years + stage + strata(subtype), data = dat)
m_pur   <- coxph(Surv(time, event) ~ scale(pc1) + age_years + stage + scale(purity), data = dat)
```

**Trade-offs.** Costs sample size in the sections that did not need every
covariate (B1's proliferation subset is only 755→~740 patients with RNA-seq,
so a fully-common cohort shrinks B2/B3 unnecessarily). A reasonable middle
ground: keep per-test cohorts but log `n` alongside every logged statistic
so the results table cannot be misread — `log_result()` already has a `note`
field for exactly this.

---

## M3. `tapply()` over a factor with unused levels yields `NA`, then `if (NA)` errors

**File/lines:** `10:137–138`, and `10:280` (`min(event_counts)`).

```r
event_counts <- tapply(meth_c$event, meth_c$subtype, sum)
if (any(event_counts < 10)) { ... }
```

**Why it's an issue.** `meth_c$subtype` is a factor built from the full
PAM50 table; complete-case filtering can leave a level with zero rows (e.g.
a "Normal-like" call that no complete-case patient carries). `tapply` then
returns `NA` for that cell, `any(NA < 10)` is `NA`, and `if (NA)` raises
`missing value where TRUE/FALSE needed` — the script dies after the two
multi-minute screens have already run but before any output is written.
`min(event_counts)` at line 280 has the same exposure and would write `NA`
into the manifest.

**Fix:**

```r
meth_c$subtype <- droplevels(meth_c$subtype)
event_counts <- tapply(meth_c$event, meth_c$subtype, sum)
stopifnot("Empty subtype stratum after complete-case filtering" =
            all(!is.na(event_counts)))
if (any(event_counts < 10)) { ... }
```

`droplevels()` here also matters for `strata()`: an empty factor level in a
stratification term is silently ignored by `coxph` but changes nothing, so
this is defensive rather than result-changing — but the crash is real.

---

## M4. `sapply()` collapses to a vector when `n_repeats == 1`

**File/lines:** `08:292–295` (`n_repeats` is user-configurable at `08:37`).

```r
all_risk_scores <- sapply(results, function(x) x$risk)   # matrix if n_repeats > 1
mean_cv_risk    <- rowMeans(all_risk_scores)             # errors if it's a vector
```

**Why it's an issue.** With `n_repeats <- 1` (a natural choice for a quick
smoke test — and `05_multivariate_cox.R:143` already sets exactly that),
`sapply` returns a plain numeric vector and `rowMeans` fails with
`'x' must be an array of at least two dimensions`, after the full screen has
run.

**Fix:**

```r
all_risk_scores <- vapply(results, function(x) x$risk,
                          numeric(nrow(meth_clinical)))
if (is.null(dim(all_risk_scores))) all_risk_scores <- matrix(all_risk_scores, ncol = 1)
mean_cv_risk <- rowMeans(all_risk_scores)
```

`vapply` additionally asserts each element's length, catching a short/long
`risk` vector from a malformed worker result.

---

## M5. `all(cf == 0)` propagates `NA` into an `if()`

**File/lines:** `08:280–283`.

```r
cf <- as.matrix(coef(fit_k, s = "lambda.1se"))[colnames(X_k)[gene_cols_k], 1]
if (all(cf == 0)) { ... }
```

**Why it's an issue.** If any name in `colnames(X_k)[gene_cols_k]` is absent
from the coefficient matrix rownames, the subset yields `NA` for that
element; `all(c(0, NA) == 0)` is `NA`; `if (NA)` errors inside a parallel
worker, where the traceback is truncated and hard to attribute.

**Fix:**

```r
coef_mat <- as.matrix(coef(fit_k, s = "lambda.1se"))
gene_names_k <- colnames(X_k)[gene_cols_k]
stopifnot(all(gene_names_k %in% rownames(coef_mat)))
cf <- coef_mat[gene_names_k, 1]
if (all(cf == 0, na.rm = TRUE)) { ... }
```

---

## M6. Fixed evaluation times can exceed follow-up

**File/lines:** `08:355–367` and `08:425` (`times = c(365, 730, 1825, 3652)`).

**Why it's an issue.** `timeROC` and `riskRegression::Score` are asked for
estimates at 10 years unconditionally. If the cohort's maximum follow-up (or
maximum *event* time) is below a requested horizon, the IPCW estimate is
undefined and you get `NaN`/`NA` AUC and Brier values — or an error from
`Score()` — which then flow into `08:369` and the AUC-over-time plot as
blank points. With PFI as the endpoint, censoring is heavy and long-horizon
estimates are the least stable.

**Fix:**

```r
eval_times <- c(365, 730, 1825, 3652)
max_t <- max(meth_clinical$time[meth_clinical$event == 1], na.rm = TRUE)
dropped <- eval_times[eval_times >= max_t]
if (length(dropped) > 0) {
  warning("Dropping evaluation times beyond last event time (",
          round(max_t), " days): ", paste(dropped, collapse = ", "))
}
eval_times <- eval_times[eval_times < max_t]
stopifnot(length(eval_times) > 0)
```

then pass `eval_times` to both `timeROC()` and `Score()` and use it for the
x-axis at `08:425` instead of the hard-coded literal (which is currently
duplicated — a second, independent bug waiting to happen if the times are
ever edited in one place only).

---

## M7. Log-rank df is hard-coded and p-values lose precision

**File/lines:** `08:378–379`, `13:198–199`.

```r
lr_test <- survdiff(surv_obj ~ risk_group)
lr_pval <- 1 - pchisq(lr_test$chisq, df = 1)
```

**Why it's an issue.** Two separate problems. `df = 1` is correct only for a
two-group comparison; it is right today but breaks silently (wrong p, not an
error) if the split ever becomes tertiles. And `1 - pchisq(q, df)` underflows
to exactly `0` for large `q` — reporting `p = 0` in a figure caption is
indefensible in a paper.

**Fix:**

```r
lr_pval <- pchisq(lr_test$chisq, df = length(lr_test$n) - 1, lower.tail = FALSE)
```

Same change at `13:199`, where the `tryCatch` already handles the
single-group case by returning `NULL`.

---

## M8. Global median split applied to per-subtype panels can empty an arm

**File/lines:** `13:179–182` (global split) consumed at `13:193–214`
(per-subtype panels).

**Why it's an issue.** The High/Low threshold is the median of the *whole*
cohort, but panels are drawn per subtype. Because promoter methylation
differs systematically by subtype (this is B2's finding), a small stratum can
end up with all patients on one side. `survfit(sub_surv ~ cbx5_group)` then
returns a single curve while `plot(..., col = c("steelblue","tomato"))`
supplies two colours — the surviving arm is drawn in `steelblue` regardless
of which arm it is, and the legend claims both are present with an `n = 0`
entry. A reader sees a mislabelled figure; nothing errors.

**Fix:**

```r
for (st in subtypes_present) {
  sub_df <- meth_c[meth_c$subtype == st, ]
  sub_df$cbx5_group <- droplevels(sub_df$cbx5_group)
  if (nlevels(sub_df$cbx5_group) < 2) {
    plot.new(); title(main = sprintf("%s: single methylation group (n=%d) -- not shown",
                                     st, nrow(sub_df)))
    next
  }
  grp_cols <- c("Low methylation" = "steelblue", "High methylation" = "tomato")
  cols_here <- grp_cols[levels(sub_df$cbx5_group)]
  ...
  plot(km_fit, col = cols_here, ...)
}
```

**Trade-offs.** A per-subtype median split would guarantee two populated
arms and better match the `strata(subtype)` model, but it makes the panels
non-comparable to each other (a different threshold per panel) and to the
Cox model, which uses the continuous value. Keeping the global threshold and
handling the degenerate panel explicitly, as above, is the more honest
option.

---

## M9. Same-day reruns silently overwrite a previous run's outputs

**File/lines:** `00_utils.R:6–18`, called at `08:42`, `09:46`, `10:44`,
`11:38`, `12:73`.

```r
full_path <- file.path(base, paste0(Sys.Date(), "_", tag))
dir.create(file.path(full_path, subdir), recursive = TRUE, showWarnings = FALSE)
```

**Why it's an issue.** The run directory is keyed on date + tag only, and
`showWarnings = FALSE` suppresses the "already exists" notice. Two runs on
the same day with different parameters (`N_PERM`, `n_screen`, a changed
filter) write into the same directory, and the second silently overwrites
the first's tables, models, and manifest. Given `write_manifest()` exists
precisely to make runs reproducible and attributable, this undercuts it.
Script 13 and script 22's "reuse the latest run" helpers then pick up a
directory whose manifest may describe a *different* run than the artefacts
sitting next to it.

**Fix:**

```r
make_run_dir <- function(tag, base, allow_existing = FALSE) {
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  full_path <- file.path(base, paste0(stamp, "_", tag))
  if (dir.exists(full_path) && !allow_existing) {
    stop("Run directory already exists: ", full_path)
  }
  ...
}
```

**Trade-offs.** Timestamped directories break the `order(basename(...),
decreasing = TRUE)[1]` "latest run" lookups in `13:235–241` only in the sense
that they keep working — `YYYY-MM-DD_HHMMSS` still sorts lexicographically
in chronological order, so those helpers are unaffected. The real cost is
directory sprawl (one dir per run rather than per day), which is a fair
price for not losing results.

---

## L1. Row-wise `apply()` where a vectorised routine exists

**File/line:** `08:164`.

```r
probe_var <- apply(gene_mat, 1, var)
```

`apply()` over 4,141 rows splits the matrix into 4,141 vectors. Replace with
`matrixStats::rowVars(gene_mat)` (roughly an order of magnitude faster, and
allocation-free). Marginal here given the dominant cost is the Cox screen —
listed for completeness, not urgency.

---

## L2. `plate` extraction assumes barcode arity

**File/line:** `08:108` (only reached when `run_batch_correction = TRUE`).

```r
plate <- sapply(strsplit(meth_clinical$sample_barcode, "-"), `[`, 6)
```

Any barcode with fewer than 6 hyphen-separated fields yields `NA` silently,
and if lengths vary `sapply` may return a list rather than a character
vector, which `ComBat` will reject with an opaque error. Add:

```r
plate <- vapply(strsplit(meth_clinical$sample_barcode, "-"),
                function(x) if (length(x) >= 6) x[6] else NA_character_,
                character(1))
stopifnot("Could not parse plate from all barcodes" = !anyNA(plate))
```

---

## L3. Graphics devices are not exception-safe

**File/lines:** `08:381–433`, `13:145–164`, `13:189–216`, `13:295–304`.

Every `pdf(...) ... dev.off()` pair assumes no error in between. If a plot
call fails (e.g. M8's degenerate panel, or `ggforest` choking on a model
term), the device stays open: subsequent plots are written into the wrong
file and the PDF on disk is truncated/corrupt. Wrap as:

```r
pdf(path, width = 7, height = 5)
on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)   # inside a function
```

or, at top level, `tryCatch({ ... }, finally = dev.off())`.

---

## Security assessment

I looked for the categories you listed and want to be straight about what is
and isn't applicable here: this is a local, single-user analysis pipeline
that reads fixed local paths and writes to a local output tree. There is no
network listener, no authentication, no user-supplied input reaching an
interpreter, and no SQL. So there is no injection or auth surface to report,
and I'd rather say that than manufacture findings.

Two things that are genuinely worth noting:

**S1 — `system()` call with a fixed command (`00_utils.R:25`).** Not an
injection risk (nothing external is interpolated into the string), but it is
an unguarded shell dependency in the failure path described in H3. The fix
there covers it.

**S2 — Patient identifiers written into the repo tree (`08:456–461`).**

```r
write.csv(data.frame(sample_barcode = meth_clinical$sample_barcode,
                      patient_id = meth_clinical$patient_id, ...),
          file.path(paths$tables, "methylation_cv_risk_scores.csv"))
```

TCGA barcodes are public, pseudonymous identifiers, so this is low risk by
itself — but the file lands under `02_analysis/`, inside the git working
tree, alongside per-patient risk scores and event status. Confirm
`02_analysis/` is in `.gitignore` (the repo has one; I did not verify its
contents), and keep it that way, so that per-patient derived data is never
committed. If this pipeline is ever pointed at non-public clinical data, the
same code path becomes a real disclosure risk with no code change required —
which is the argument for fixing the habit now rather than later.

---

## Summary table

| ID | Severity | File:line | One-line |
|----|----------|-----------|----------|
| C1 | Critical | 08:140, 09:78, 10:60, 11:67, 13:50 | Unkeyed annotation cache can silently misalign all probes |
| C2 | Critical | 08:151, 10:68, 11:75, 12:230, 13:59 | `rowsum`/`table` ordering assumed, never asserted |
| H1 | High | 08:238/290, 11:174/197, 13:277/290 | `stopCluster` skipped on error → leaked workers |
| H2 | High | 08:244, 11:170, 13:280 | Whole global env (incl. 290 MB `m_mat`) exported per worker |
| H3 | High | 00_utils:25 (08:443) | `git` failure aborts run before results are written |
| H4 | High | 13:253, 13:295 | Fragile manifest parsing → figure silently missing observed marker |
| M1 | Medium | 08:95 | `purity_available` all-or-nothing gate |
| M2 | Medium | 09:303/364/498 | Models compared across different complete-case subsets |
| M3 | Medium | 10:137, 10:280 | `tapply` `NA` for empty stratum → `if (NA)` crash |
| M4 | Medium | 08:292 | `sapply` collapses when `n_repeats == 1` |
| M5 | Medium | 08:280 | `all(cf == 0)` with `NA` → crash in worker |
| M6 | Medium | 08:355, 08:425 | Fixed eval times may exceed follow-up; literal duplicated |
| M7 | Medium | 08:379, 13:199 | Hard-coded `df = 1`; `1 - pchisq` underflow |
| M8 | Medium | 13:179 | Global median split can empty a per-subtype arm |
| M9 | Medium | 00_utils:6 (all) | Same-day reruns overwrite prior outputs |
| L1 | Low | 08:164 | `apply(_, 1, var)` → `matrixStats::rowVars` |
| L2 | Low | 08:108 | Barcode plate extraction unguarded |
| L3 | Low | 08:381+, 13:145+ | PDF devices not exception-safe |
| S1 | Info | 00_utils:25 | Fixed-string `system()` call, no injection surface |
| S2 | Low | 08:456 | Per-patient identifiers written inside repo tree |

## Suggested order of work

1. **C1 and C2 first.** They are the only findings that can produce wrong
   numbers with no visible symptom, and both are cheap to fix. Do them
   together via the shared `build_promoter_gene_mat()` helper, then rerun
   one script and diff `gene_mat` dimensions and a few known values
   (`CBX5`, `CBX2`) against your current results to confirm nothing moved.
2. **H3 and H4** — small, and they protect work you have already paid for.
3. **H1/H2** as one refactor (wrap each parallel block in a function).
4. The M-series opportunistically; M3 and M4 are two-line changes that turn
   late-stage crashes into clean failures.

#!/bin/bash
# =============================================================================
# test_mode/run_parallel.sh
#
# Dependency-ordered runner for the pipeline, in two modes.
#
# TEST mode (default): scripts 02-15 on the sample dataset, independent
# scripts in parallel, everything stopped after LIMIT seconds. The default,
# 180 s, is the reference (unloaded runs take ~100-140 s). When other R jobs
# are running, pass LIMIT explicitly (e.g. LIMIT=600); that run is then not a
# timing measurement.
#
#   bash test_mode/run_parallel.sh            # from the project root
#   LIMIT=600 bash test_mode/run_parallel.sh   # other R jobs running: no timing
#
#   Needs the sample dataset first: Rscript test_mode/make_sample_dataset.R
#   Works in a fresh _test_run/ copy of the project layout, so the real
#   00_raw_data/, 01_processed_data/ and 02_analysis/ are never written.
#   Scripts 01 and 06 (GDC download) are replaced by the saved sample files,
#   and 02 reads the saved clinical table, so the run makes no network call.
#   Order (-> = waits for, | = runs at the same time):
#     aliquots_plant -> 02 -> aliquots_verify (test_mode/check_aliquots.R: duplicate
#     RNA aliquots planted in the sample SE copy, 02's one-per-patient choice checked)
#     -> 03 | (07 -> annotation cache) -> 04 | 05 | 08 | 09 | 10 | 11 | 12 -> 13
#     14 runs alongside 08-12 (same inputs as 11; does not read 10's or 11's output)
#     15 runs alongside 08-12 (reads only 07's M matrix, the annotation cache and GSE141441)
#
# FULL mode: the real cohort, in place from the project root (reads and writes
# the real 00_raw_data/, 01_processed_data/, 02_analysis/), MBC_SUBSET_N unset,
# no time limit unless LIMIT is set, logs in run_logs/full_<script>.log.
# A LIMIT hit kills scripts mid-write, which can leave a truncated real
# 01_processed_data/*.rds or 02_analysis/ output: rerun from 07 (or FROM=02)
# after a TIMEOUT.
#
#   FULL=1 bash test_mode/run_parallel.sh             # 07 -> 14
#   FULL=1 FROM=02 bash test_mode/run_parallel.sh     # 02-05 first, then 07 -> 14
#   FULL=1 DRY=1 bash test_mode/run_parallel.sh       # print commands + order, run nothing
#
#   01 and 06 (downloads) are never run. 02-05 run only with FROM=02, as ONE
#   R session sourcing 02, 03, 04, 05 in order (outside test mode 03 and 04
#   use objects the previous script left in the session). 02 still re-fetches
#   the GDC clinical table when 00_raw_data/gdc_clinical.rds is over 7 days old.
#   Order: [02-05 session ->] 07 -> annotation cache
#            -> 09 | 10 | 12 | 15 (light, uncapped) | 08 | 14 (heavy)
#            -> 11 (once 08 or 14 finishes; the others may still run) -> 13
#   The heavy scripts (08, 11, 13, 14) each use detectCores()-1 workers, so no
#   more than two of them run at once.
#
# DRY=1 (either mode) prints each Rscript command line (arguments quoted) in run
# order -- jobs of one parallel batch in any order -- instead of running it, and
# touches no file; it ends with "RESULT: dry run" and exit 0. Test-mode commands
# run inside _test_run/ with MBC_SUBSET_N=1 -- never paste them at the root.
# Exit codes: 0 all passed, 1 failed or skipped, 2 over LIMIT.
# Rscript: $RSCRIPT if non-empty (a path is taken from the directory you run
# from, a bare name is searched on PATH), else Rscript on PATH, else
# /c/Program Files/R/R-4.6.0/bin/Rscript.exe; the run stops if none exists (a
# DRY run needs none). Reported runs used R 4.6.0.
# =============================================================================

set -u
FULL="${FULL:-0}"
DRY="${DRY:-0}"
case "$FULL$DRY" in 00|01|10|11) ;; *) echo "FULL and DRY must be 0 or 1, got FULL='$FULL' DRY='$DRY'"; exit 1 ;; esac
# Rscript (see header). Resolved before the cd below, so a relative RSCRIPT is
# taken from the caller's directory; made absolute because test mode runs the
# steps from _test_run/.
WINR="/c/Program Files/R/R-4.6.0/bin/Rscript.exe"
if [ "$DRY" = 0 ]; then
  if [ -n "${RSCRIPT:-}" ]; then
    R=$(type -P -- "$RSCRIPT") || { echo "RSCRIPT='$RSCRIPT' not found (a path must be an executable file, a bare name must be on PATH); fix it, or unset it for the default lookup"; exit 1; }
  else
    R=$(type -P -- Rscript || type -P -- "$WINR") || { echo "Rscript not found (tried: Rscript on PATH, $WINR). Install R (reported runs used 4.6.0) or set RSCRIPT=/path/to/Rscript."; exit 1; }
  fi
  case "$R" in /*|\\*|[A-Za-z]:[/\\]*) ;; *) R="${PWD%/}/$R" ;; esac
  echo "Rscript: $R"
fi
cd "$(dirname "$0")/.." || exit 1
ROOT=$(pwd)
[ "$FULL" = 0 ] && [ -n "${FROM:-}" ] && { echo "FROM applies only with FULL=1"; exit 1; }
FROM="${FROM:-07}"
DRYNOTE=""
[ "$DRY" = 1 ] && DRYNOTE=", DRY run (nothing executed)"
# Script 15's GEO input, not fetched by any script. Without it only 15 is
# skipped (the run then reports incomplete); a DRY run still prints its line.
GSE=00_raw_data/external/GSE141441/GSE141441_series_matrix.txt.gz
RUN15=1
[ -f "$GSE" ] || { echo "note: missing $GSE -- 15_external_validation would be skipped"; [ "$DRY" = 1 ] || RUN15=0; }

if [ "$FULL" = 1 ]; then
  LIMIT="${LIMIT:-}"
  [ -z "$LIMIT" ] || [[ "$LIMIT" =~ ^[1-9][0-9]*$ ]] || { echo "LIMIT must be a whole number of seconds > 0, got '$LIMIT'"; exit 1; }
  case "$FROM" in 02|07) ;; *) echo "FROM must be 02 or 07, got '$FROM'"; exit 1 ;; esac
  (( BASH_VERSINFO[0] > 5 || (BASH_VERSINFO[0] == 5 && BASH_VERSINFO[1] >= 1) )) ||
    { echo "FULL=1 needs bash 5.1+ (wait -n on given jobs keeps the heavy cap), got $BASH_VERSION"; exit 1; }
  echo "MODE: FULL cohort, in place in $ROOT, from $FROM, limit ${LIMIT:-none}$DRYNOTE"
  META=00_raw_data/TCGA_BRCA_Methylation450_SE.rds
  if [ ! -f "$META" ]; then
    [ "$DRY" = 1 ] || { echo "FULL=1 refused: $META is missing (download it with 06 first)"; exit 1; }
    echo "note: $META is missing -- a real FULL run would refuse"
  fi
  # 07 reads 02's merged_data.rds, which is not tracked (a fresh clone has none)
  if [ "$FROM" = 07 ] && [ ! -f 01_processed_data/merged_data.rds ]; then
    [ "$DRY" = 1 ] || { echo "FULL=1 refused: 01_processed_data/merged_data.rds is missing -- run with FROM=02 (02 builds it)"; exit 1; }
    echo "note: 01_processed_data/merged_data.rds is missing -- a real FROM=07 run would refuse"
  fi
  unset MBC_SUBSET_N
  L="run_logs/full_"
  if [ "$DRY" != 1 ]; then
    mkdir -p run_logs
    # one FULL run per checkout: two would overwrite each other's in-place outputs
    mkdir run_logs/full_LOCK 2> /dev/null ||
      { echo "FULL=1 refused: another FULL run holds run_logs/full_LOCK (remove it if that run is dead)"; exit 1; }
    # on Ctrl-C the background pipeline can outlive this shell: keep the lock then
    trap 'if [ -n "${PIPE:-}" ] && kill -0 "$PIPE" 2> /dev/null; then echo "pipeline still running -- run_logs/full_LOCK kept"; else rmdir run_logs/full_LOCK 2> /dev/null; fi' EXIT
    rm -f "${L}FAILED" "${L}SKIPPED" "${L}TIMEOUT"
  fi
else
  LIMIT="${LIMIT:-180}"
  [ -z "$LIMIT" ] || [[ "$LIMIT" =~ ^[1-9][0-9]*$ ]] || { echo "LIMIT must be a whole number of seconds > 0, got '$LIMIT'"; exit 1; }
  echo "MODE: TEST (sample dataset in _test_run/), limit ${LIMIT}s$DRYNOTE"
  [ "$DRY" = 1 ] && echo "note: a real test run executes these inside _test_run/ with MBC_SUBSET_N=1"
  SAMPLE="$ROOT/00_raw_data/sample"
  W="$ROOT/_test_run"
  L="run_logs/"

  for f in TCGA_BRCA_SummarizedExperiment.rds TCGA_BRCA_Normal_SummarizedExperiment.rds \
           TCGA_BRCA_Methylation450_SE.rds; do
    [ -f "$SAMPLE/$f" ] || { echo "Missing $SAMPLE/$f -- run: Rscript test_mode/make_sample_dataset.R"; exit 1; }
  done
  [ -f 00_raw_data/gdc_clinical.rds ] || { echo "Missing 00_raw_data/gdc_clinical.rds -- run: Rscript test_mode/make_sample_dataset.R"; exit 1; }

  # - - - fresh workspace (not in a DRY run, which leaves _test_run/ as it is) - - -
  if [ "$DRY" != 1 ]; then
    rm -rf "$W"
    mkdir -p "$W"/00_raw_data/TCGAsupp "$W"/01_processed_data \
             "$W"/02_analysis "$W"/03_figures "$W"/run_logs "$W"/$(dirname "$GSE")
    cp -r deployed_scripts_v1 "$W"/
    mkdir -p "$W"/test_mode && cp test_mode/check_aliquots.R "$W"/test_mode/
    cp 00_raw_data/TCGAsupp/TCGA-CDR-SupplementalTableS1.xlsx "$W"/00_raw_data/TCGAsupp/
    cp "$SAMPLE"/TCGA_BRCA_*.rds "$W"/00_raw_data/
    cp 00_raw_data/gdc_clinical.rds "$W"/00_raw_data/   # fresh copy: 02 reads it, no GDC call
    [ -f "$GSE" ] && cp "$GSE" "$W/$GSE"
    cp 01_processed_data/PAM50_classification.txt 01_processed_data/aran_purity_tcga.xlsx \
       01_processed_data/epifactors_*.csv "$W"/01_processed_data/
    cd "$W" || exit 1
  fi

  export MBC_SUBSET_N="${MBC_SUBSET_N:-1}"   # turns test mode on in every script
fi
START=$(date +%s)

step() {  # step <name> <Rscript args...> : one Rscript session, log in ${L}<name>.log
  local name=$1; shift
  if [ "$DRY" = 1 ]; then
    echo "$(date +%T) START $name: Rscript$(printf " '%s'" "$@") > ${L}$name.log"
    return 0
  fi
  local t=$(date +%s)
  "$R" "$@" > "${L}$name.log" 2>&1
  local rc=$?
  echo "$(date +%T) $name rc=$rc $(( $(date +%s)-t ))s"
  [ $rc -eq 0 ] || echo "$name" >> "${L}FAILED"
  return $rc
}

run() { step "$1" "deployed_scripts_v1/$1.R"; }

warm_cache() {  # build the probe annotation cache once, before 08-13 read it in parallel,
                # then check it is keyed on this run's m_mat (a stale cache fails the step)
  # TEST mode only: plant a well-formed stale cache (same row count, other order)
  # first, so the key-mismatch rebuild is exercised on every test run
  local seed=""
  [ "$FULL" = 1 ] || seed='saveRDS(list(key = rev(rownames(m)), value = data.frame(Name = rev(rownames(m)))), "01_processed_data/ann_matched_cache.rds"); '
  step annotation_cache -e 'library(minfi); library(IlluminaHumanMethylation450kanno.ilmn12.hg19); source("deployed_scripts_v1/00_utils.R"); m <- readRDS("01_processed_data/meth_m_matrix.rds"); '"$seed"'invisible(load_matched_annotation(m)); assert_annotation_cache(m)'
}

skip() {  # record scripts that could not start because an input step failed
  for s in "$@"; do echo "$(date +%T) $s SKIPPED (input step failed)"; echo "$s" >> "${L}SKIPPED"; done
}

METH_ALL="08_methylation_survival 09_pc1_biological_characterization 10_subtype_purity_stratified_screen
          11_stratified_screen_validation 12_cbx_family_stratified_screen 13_cbx5_figures
          14_sign_bias_diagnostics 15_external_validation"

test_pipeline() {
  # plant duplicate RNA aliquots in the sample SE, then check 02 keeps one per patient
  local rest="03_normalization_and_qc 04_survival_analysis 05_multivariate_cox 07_methylation_prepare $METH_ALL"
  step aliquots_plant test_mode/check_aliquots.R plant ||
    { skip 02_load_and_prepare_data aliquots_verify $rest; return 1; }
  run 02_load_and_prepare_data || { skip aliquots_verify $rest; return 1; }
  step aliquots_verify test_mode/check_aliquots.R verify || { skip $rest; return 1; }
  run 03_normalization_and_qc &
  local p03=$!
  { run 07_methylation_prepare && warm_cache; } &
  local p07=$!
  wait $p03; local r03=$?
  wait $p07; local r07=$?

  local pids=()
  if [ $r03 -eq 0 ]; then
    run 04_survival_analysis & pids+=($!)
    run 05_multivariate_cox & pids+=($!)
  else
    skip 04_survival_analysis 05_multivariate_cox
  fi
  if [ $r07 -ne 0 ]; then
    skip $METH_ALL
  else
    run 08_methylation_survival & pids+=($!)
    run 09_pc1_biological_characterization & pids+=($!)
    run 10_subtype_purity_stratified_screen & pids+=($!)
    run 12_cbx_family_stratified_screen & pids+=($!)
    { run 11_stratified_screen_validation || { skip 13_cbx5_figures; false; } && run 13_cbx5_figures; } & pids+=($!)
    run 14_sign_bias_diagnostics & pids+=($!)
    if [ "$RUN15" = 1 ]; then run 15_external_validation & pids+=($!); else skip 15_external_validation; fi
  fi
  wait "${pids[@]}"
}

full_pipeline() {
  if [ "$FROM" = 02 ]; then   # 02 -> 03 -> 04 -> 05 in one R session, like the README's full run
    step 02-05_rna_session -e 'for (s in c("02_load_and_prepare_data", "03_normalization_and_qc", "04_survival_analysis", "05_multivariate_cox")) { message("== ", s); source(file.path("deployed_scripts_v1", paste0(s, ".R")), print.eval = TRUE) }' \
      || { skip 07_methylation_prepare $METH_ALL; return 1; }
  fi
  { run 07_methylation_prepare && warm_cache; } || { skip $METH_ALL; return 1; }

  local light=() heavy=()
  run 09_pc1_biological_characterization & light+=($!)   # light, uncapped
  run 10_subtype_purity_stratified_screen & light+=($!)
  run 12_cbx_family_stratified_screen & light+=($!)
  if [ "$RUN15" = 1 ]; then run 15_external_validation & light+=($!); else skip 15_external_validation; fi
  run 08_methylation_survival & heavy+=($!)              # heavy 1 of 2
  run 14_sign_bias_diagnostics & heavy+=($!)             # heavy 2 of 2
  wait -n "${heavy[@]}"   # a heavy slot is free: 11 -> 13 can take it (bash 5.1+)
  [ "$DRY" = 1 ] && echo "         (waits here for the first of heavy 08 / 14; the rest keep running)"
  # 11 reads only 07's outputs and the cache, so it runs whatever 08-14 returned (a
  # failed 08 leaves 11/13 no stale input). 13 also reads 11's newest
  # *_stratified_screen_validation run, so it runs only after 11 passed.
  if run 11_stratified_screen_validation; then run 13_cbx5_figures; else skip 13_cbx5_figures; fi
  wait "${heavy[@]}" "${light[@]}"
}

if [ "$FULL" = 1 ]; then full_pipeline & else test_pipeline & fi
PIPE=$!
WATCH=""
if [ -n "$LIMIT" ] && [ "$DRY" != 1 ]; then
  # sleep gets no stdout/stderr: a killed watcher must not leave an orphaned
  # sleep holding a pipe open (`... | tail` would hang until LIMIT)
  ( sleep "$LIMIT" > /dev/null 2>&1 < /dev/null
    kill -0 $PIPE 2> /dev/null || exit 0   # finished right at LIMIT: not a timeout
    echo "TIME LIMIT ${LIMIT}s reached -- stopping"; touch "${L}TIMEOUT"
    taskkill //F //T //PID "$(cat /proc/$PIPE/winpid)" > /dev/null 2>&1   # this run's process tree only
    kill $PIPE 2> /dev/null ) &
  WATCH=$!
fi
wait $PIPE 2> /dev/null
[ -n "$WATCH" ] && kill $WATCH 2> /dev/null

[ "$DRY" = 1 ] && { echo "RESULT: dry run -- nothing executed"; exit 0; }
echo "total $(( $(date +%s)-START ))s"
if [ -f "${L}TIMEOUT" ]; then
  echo "RESULT: over the ${LIMIT}s limit"
  [ "$FULL" = 1 ] && echo "outputs written in place may be truncated -- rerun before using them"
  exit 2
fi
[ -f "${L}SKIPPED" ] && echo "skipped: $(tr '\n' ' ' < "${L}SKIPPED")"
if [ -f "${L}FAILED" ]; then echo "RESULT: failed: $(tr '\n' ' ' < "${L}FAILED")"; exit 1; fi
if [ -f "${L}SKIPPED" ]; then echo "RESULT: incomplete -- scripts skipped"; exit 1; fi
if [ "$FULL" = 1 ]; then
  echo "RESULT: all passed -- logs in run_logs/full_*.log, results in 02_analysis and 03_figures"
else
  echo "RESULT: all passed -- logs in _test_run/run_logs, results in _test_run/02_analysis and _test_run/03_figures"
fi

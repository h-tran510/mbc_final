# Contributing

Thank you for looking at this project. Questions, bug reports and suggestions
are welcome as GitHub issues; please include the script name, the command you
ran, and the `run_manifest.txt` or log of the failing run.

## Changes

1. **Branch.** Never commit to `main`. Create a feature branch from the latest
   `main` (`fix/<topic>`, `feat/<topic>`, `docs/<topic>`, `chore/<topic>`).
2. **Test mode before the PR.** Run the test mode from the repository root and
   paste its last lines (the `RESULT:` line and timings) into the PR:

   ```bash
   Rscript test_mode/make_sample_dataset.R   # once, if 00_raw_data/sample/ is missing
   bash test_mode/run_parallel.sh
   ```

   The sample builder needs the manually downloaded TCGA-CDR table, the
   tracked PAM50 and purity files and the GDC Data Transfer Tool (README,
   Inputs).

   Without the optional GSE141441 file the run ends `RESULT: incomplete` because
   script 15 is skipped; say so in the PR. A change that only touches
   documentation says so instead. A change to a
   number-producing step also states the full-run numbers before and after, with
   the run folder they came from.
3. **Pull request.** Open a PR into `main` with a table of what changed, why,
   before and after, and name any check you did not run.
4. **Review.** Every PR gets a code review before it is merged; every finding is
   fixed or answered, and the review repeats until it is clean. The maintainer
   merges.

## Data

Never commit data or run outputs: `00_raw_data/`, `02_analysis/`, `run_logs/`,
`_test_run/`, any `.rds`, `.RData` or `.rda`, R session files. `.gitignore`
ignores all of these; still read `git status --porcelain` before each commit. The third-party files that are tracked
under `01_processed_data/` are listed, with their licences, in the README's
Inputs table; a new input needs a row there (source URL, version, bytes, md5,
licence) in the same PR. Measure bytes and md5 on the staged blob after
`git add` (`git cat-file -s :<path>`, `git show :<path> | md5sum`), not on a
Windows working copy, where line-ending conversion changes both.

## Conduct

Be respectful and specific. Report results as measured, including results that
do not support the hypothesis.

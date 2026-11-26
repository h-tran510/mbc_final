# Docker image

`Dockerfile` (repository root) builds the pipeline's R environment: R 4.6.0 on
`rocker/r-ver:4.6.0` (pinned by digest), Bioconductor 3.23 via BiocManager, and
every package at the version pinned in `renv.lock`, CRAN packages from a
date-pinned Posit Package Manager snapshot (`CRAN_DATE`, default
2026-09-28). The image holds no code and no data; mount the repository at
`/work`. The system libraries (apt) are not pinned: they come from the Ubuntu
archive at build time.

## Build

From the repository root (`.dockerignore` sends only `renv.lock` to the
build):

    docker build -t mbc26 .

Bioconductor packages build from source, so the first build takes a while.
A Bioconductor version superseded by a later 3.23 patch is fetched from the
release's archive.

## Inputs the repository does not carry

`00_raw_data/` is git-ignored. Test mode needs, in the mounted repository:

- `00_raw_data/TCGAsupp/TCGA-CDR-SupplementalTableS1.xlsx`: TCGA Pan-Cancer
  Clinical Data Resource, Table S1 (Liu et al., Cell 2018). No script
  downloads it.
- `00_raw_data/sample/` and `00_raw_data/gdc_clinical.rds`: the test-mode
  sample dataset and the GDC clinical table, both written by
  `test_mode/make_sample_dataset.R` (next section).
- optional, for script 15 only:
  `00_raw_data/external/GSE141441/GSE141441_series_matrix.txt.gz` (GEO);
  without it script 15 is skipped and the runner ends "RESULT: incomplete"
  (exit 1) even when every other script passes.

## Run in test mode

Test mode runs scripts 02-15 on a small, event-enriched sample
(`test_mode/run_parallel.sh`). The `--user` flag keeps files written into the
mounted repository owned by you rather than root; `HOME=/tmp` gives R a
writable home for its caches; `--init` lets Ctrl-C and `docker stop` end the
run.

Build the sample dataset once (it downloads from the GDC with `gdc-client`,
which the image does not include; TCGAbiolinks uses one from the repository
root or fetches its own). A `gdc-client` already in the root of a checkout
shared with Windows is the Windows build and will not run in the container:
move it aside first, or build the sample on the host with the same command
outside Docker. Likewise, a Linux `gdc-client` the container fetched into the
root will not run on a Windows host: move it aside before a host-side GDC
download.

    docker run --rm --init --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/work mbc26 \
      Rscript test_mode/make_sample_dataset.R

Then run the pipeline in test mode (the image's default command). The
runner's 180 s default `LIMIT` was measured on a native Windows run; a
container, above all one on a Windows drive mounted through WSL, is slower,
so give it more:

    docker run --rm --init --user "$(id -u):$(id -g)" -e HOME=/tmp -e LIMIT=900 -v "$PWD":/work mbc26

Outputs and logs go to `_test_run/` inside the mounted repository. Mount a
main clone, not a `git worktree`: a worktree's `.git` points to a host path
the container cannot reach, so the run manifests would record no commit.

`test_mode/run_parallel.sh` must have LF line endings for bash in the
container; `.gitattributes` gives every new checkout LF. A Windows clone made
before that keeps its CRLF copy: refresh it once with
`rm test_mode/run_parallel.sh && git checkout -- test_mode/run_parallel.sh`.

At `LIMIT` the runner kills its R jobs with Windows `taskkill`,
which Linux lacks: the jobs stop here only because the container exits when
the runner does. Start the runner as the container's command (above), not with
`docker exec` into a running container.

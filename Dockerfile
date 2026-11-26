# R environment for the mbc_26 pipeline: R 4.6.0 + Bioconductor 3.23, every
# package pinned by renv.lock. The image holds the environment only; the code
# and data are mounted at run time (see docker/README.md).
# rocker/r-ver:4.6.0, pinned by digest (multi-arch index, 2026-06-24)
FROM rocker/r-ver:4.6.0@sha256:6f05a1a8b8c52328f181593923909b01cbfd14c9caea93bf75ddc65e806d8eac

# Build tools and system libraries for the source builds in renv.lock
# (curl, openssl, xml2/XML, RCurl, Rhtslib, png, jpeg, nloptr); git for the
# commit hash and dirty flag the run manifests record (autocrlf input: a
# Windows checkout's CRLF files do not read as modified).
RUN apt-get update && apt-get install -y --no-install-recommends \
      git cmake pkg-config libbz2-dev libcurl4-openssl-dev libjpeg-dev liblzma-dev \
      libpng-dev libssl-dev libxml2-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* \
    && git config --system --add safe.directory /work \
    && git config --system core.autocrlf input

# CRAN from a date-pinned Posit Package Manager snapshot (the date the lock was
# read; every locked CRAN version is current or archived there): Linux binaries,
# and the same renv and BiocManager on every rebuild. Override with
# --build-arg CRAN_DATE=YYYY-MM-DD.
ARG CRAN_DATE=2026-09-28
WORKDIR /opt/mbc26
COPY renv.lock ./
# BiocManager pinned to the lockfile's Bioconductor release (3.23, matching R
# 4.6), then the exact versions from renv.lock into the site library, so no
# project activation is needed at run time. No renv cache: one copy per image.
RUN export RENV_CONFIG_CACHE_ENABLED=FALSE MAKEFLAGS="-j$(nproc)" \
    && export CRAN_REPO="https://p3m.dev/cran/__linux__/$(. /etc/os-release; echo $VERSION_CODENAME)/${CRAN_DATE}" \
    && R -q -e 'repo <- c(CRAN = Sys.getenv("CRAN_REPO")); \
         install.packages(c("BiocManager", "renv"), repos = repo); \
         BiocManager::install(version = renv::lockfile_read("renv.lock")$Bioconductor$Version, ask = FALSE, update = FALSE); \
         renv::restore(lockfile = "renv.lock", library = .libPaths()[1], repos = repo, prompt = FALSE)' \
    && rm -rf /root/.cache/R

# test_mode/run_parallel.sh finds Rscript on PATH.
WORKDIR /work
CMD ["bash", "test_mode/run_parallel.sh"]

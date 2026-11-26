# =============================================================================
# test_mode/make_sample_dataset.R
#
# Builds the fixed sample dataset that test_mode/run_parallel.sh runs on.
# Run once, from the project root. Every GDC answer is saved (gdc_cached() in
# 00_utils.R: retried, reused while younger than 7 days), so a rerun after a
# network failure resumes, and files already in 00_raw_data/GDCdata are reused:
#
#   Rscript test_mode/make_sample_dataset.R
#
# Patients are drawn only from those with every input scripts 08-13 need:
# primary-tumor RNA-seq, HM450 methylation, a PAM50 subtype, a CPE purity
# estimate, a valid TCGA-CDR PFI endpoint, stage and age. The sample is
# EVENT-ENRICHED (N_EVENT PFI events + N_CENSORED censored, vs ~13% events in
# the full cohort) so the survival models have enough events to fit. It is for
# testing the pipeline end to end, not for estimating anything.
#
# Output: 00_raw_data/sample/ (git-ignored with the rest of 00_raw_data/)
# =============================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(readxl)

source("deployed_scripts_v1/00_utils.R")

N_EVENT    <- 80
N_CENSORED <- 220
N_NORMAL   <- 10
set.seed(42)

dl_dir  <- "00_raw_data/GDCdata"
out_dir <- "00_raw_data/sample"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

patient <- function(barcode) substr(barcode, 1, 12)

# - - - 1. Queries (same as scripts 01 and 06) - - -
q_file <- function(name) file.path(out_dir, paste0("query_", name, ".rds"))
tumor_q <- gdc_cached(q_file("tumor"), function() GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts", sample.type = "Primary Tumor"))
normal_q <- gdc_cached(q_file("normal"), function() GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts", sample.type = "Solid Tissue Normal"))
meth_q <- gdc_cached(q_file("meth"), function() GDCquery(project = "TCGA-BRCA", data.category = "DNA Methylation",
                   data.type = "Methylation Beta Value",
                   platform = "Illumina Human Methylation 450", sample.type = "Primary Tumor"))

# - - - 2. Eligible patients - - -
cdr <- read_excel("00_raw_data/TCGAsupp/TCGA-CDR-SupplementalTableS1.xlsx")
cdr <- cdr[cdr$type == "BRCA" & is.na(cdr$Redaction) & !is.na(cdr$PFI) & !is.na(cdr$PFI.time) & cdr$PFI.time > 0 &
             !is.na(suppressWarnings(as.numeric(cdr$age_at_initial_pathologic_diagnosis))) &
             grepl("^Stage (I|II|III|IV)[ABC]?$", cdr$ajcc_pathologic_tumor_stage), ]

pam50  <- read.delim("01_processed_data/PAM50_classification.txt", check.names = FALSE)
purity <- read_excel("01_processed_data/aran_purity_tcga.xlsx", sheet = "Supp Data 1", skip = 2)
purity <- purity[!is.na(suppressWarnings(as.numeric(purity$CPE))), ]

eligible <- Reduce(intersect, list(
  cdr$bcr_patient_barcode,
  patient(pam50[[1]]),
  patient(purity$`Sample ID`),
  patient(getResults(tumor_q)$cases),
  patient(getResults(meth_q)$cases)
))
cdr <- cdr[cdr$bcr_patient_barcode %in% eligible, ]
cat("Eligible patients:", nrow(cdr), "| PFI events:", sum(cdr$PFI == 1), "\n")

events   <- cdr$bcr_patient_barcode[cdr$PFI == 1]
censored <- cdr$bcr_patient_barcode[cdr$PFI == 0]
chosen <- c(sample(events, min(N_EVENT, length(events))),
            sample(censored, min(N_CENSORED, length(censored))))
cat("Sample:", length(chosen), "patients,", sum(chosen %in% events), "events\n")

# - - - 3. Subset each query: one file per chosen patient - - -
keep_patients <- function(q, patients) {
  res <- q$results[[1]]
  res <- res[patient(res$cases) %in% patients, ]
  res <- res[order(res$cases), ]
  q$results[[1]] <- res[!duplicated(patient(res$cases)), ]
  q
}
tumor_q  <- keep_patients(tumor_q, chosen)
meth_q   <- keep_patients(meth_q, chosen)
normal_q$results[[1]] <- head(normal_q$results[[1]], N_NORMAL)

# - - - 4. Download anything missing, prepare, save under the names 01/06 use - - -
for (x in list(list(tumor_q,  "TCGA_BRCA_SummarizedExperiment.rds"),
               list(normal_q, "TCGA_BRCA_Normal_SummarizedExperiment.rds"),
               list(meth_q,   "TCGA_BRCA_Methylation450_SE.rds"))) {
  se <- gdc_cached(file.path(out_dir, x[[2]]), function() {
    GDCdownload(x[[1]], directory = dl_dir, method = "client", files.per.chunk = 30)
    GDCprepare(x[[1]], directory = dl_dir)
  })
  cat(x[[2]], ":", dim(se), "\n")
}

# Clinical table: the same file script 02 reads, so neither run calls GDC for it
gdc_cached("00_raw_data/gdc_clinical.rds",
           function() GDCquery_clinic(project = "TCGA-BRCA", type = "clinical"))

writeLines(sort(chosen), file.path(out_dir, "patients.txt"))
cat("Sample dataset written to", out_dir, "\n")

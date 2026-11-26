# =============================================================================
# 16_epifactors_candidate_panel.R
#
# Purpose: build a pre-specified candidate gene panel from EpiFactors, the
# methylation-machinery analog of what Xing et al. did with SeneQuest for
# senescence. Restricting the promoter screen to this panel (instead of all
# 4,141 promoters) trades an outcome-blind, mechanistically justified feature
# space for the whole-genome one -- fewer tests, smaller multiple-testing
# burden, and a defensible answer to "why these genes."
#
# TWO PANELS ARE BUILT:
#   panel_broad  -- all ~801 EpiFactors proteins (scale-matched to SeneQuest)
#   panel_narrow -- only the "DNA or RNA modification" category (DNMT/TET/
#                   MBD-family machinery; mechanistically tight to a
#                   methylation array specifically)
#
# CAUTION: do not filter on the string "methylation" alone. EpiFactors'
# category taxonomy separates DNA methylation from histone methylation
# (a much larger, unrelated family -- EZH2, SETD2, the KMT/KDM genes). Only
# panel_narrow's explicit category filter avoids pulling those in by
# accident.
#
# This script builds the candidate LIST and cross-references it against your
# own data -- it does not run a screen. Feed the resulting gene list into
# stratified_screen_one() from script 14/15, restricted to these promoters,
# as the next step.
# =============================================================================


# - - - 0. Setup - - -

source("scripts/00_utils.R")


# - - - 1. Download and inspect (cached -- EpiFactors updates infrequently) - - -

epi_path <- "01_processed_data/epifactors_v2.1_main.csv"

if (!file.exists(epi_path)) {
  download.file(
    "https://epifactors.autosome.org/public_data/v2.1/EpiGenes_main.csv",
    destfile = epi_path, mode = "wb"
  )
}

epi <- read.csv(epi_path, stringsAsFactors = FALSE, check.names = FALSE)

cat("EpiFactors table:", nrow(epi), "rows,", ncol(epi), "columns\n")
cat("\nColumn names -- INSPECT THIS before trusting the filters below.",
    "Column names have changed across EpiFactors versions before, and the",
    "exact strings used inside the function/category column are not",
    "guaranteed to match what's assumed here:\n")
print(colnames(epi))


# - - - 2. Locate the gene-symbol and function/category columns - - -
# Pattern-matched rather than hard-coded, the same defensive approach used
# for the PAM50 file in script 13 -- print what was found so a mismatch is
# caught immediately instead of silently producing an empty or wrong panel.

symbol_col <- grep("HGNC|gene.?symbol|gene.?name", colnames(epi),
                    value = TRUE, ignore.case = TRUE)[1]
function_col <- grep("function|category|type|classification", colnames(epi),
                      value = TRUE, ignore.case = TRUE)[1]

cat("\nUsing symbol column:  ", symbol_col, "\n")
cat("Using function column:", function_col, "\n")
stopifnot(
  "Could not identify a gene symbol column -- inspect colnames(epi) manually" =
    !is.na(symbol_col),
  "Could not identify a function/category column -- inspect colnames(epi) manually" =
    !is.na(function_col)
)

cat("\nUnique values in the function/category column",
    "(find the exact DNA-modification string here, it may not match",
    "\"DNA or RNA modification\" verbatim):\n")
print(sort(unique(epi[[function_col]])))


# - - - 3. Build the two panels - - -

panel_broad <- unique(epi[[symbol_col]])
panel_broad <- panel_broad[panel_broad != "" & !is.na(panel_broad)]

# Adjust this pattern after inspecting the printed unique values above --
# it is written to match "DNA" AND ("modification" OR "methylation") while
# EXCLUDING anything that also mentions "histone", specifically to avoid the
# histone-methylation trap described at the top of this script.
is_dna_mod <- grepl("DNA", epi[[function_col]], ignore.case = TRUE) &
  grepl("modif|methyl", epi[[function_col]], ignore.case = TRUE) &
  !grepl("histone", epi[[function_col]], ignore.case = TRUE)

panel_narrow <- unique(epi[[symbol_col]][is_dna_mod])
panel_narrow <- panel_narrow[panel_narrow != "" & !is.na(panel_narrow)]

cat("\npanel_broad: ", length(panel_broad), "genes\n")
cat("panel_narrow:", length(panel_narrow), "genes\n")
print(panel_narrow)

# Sanity check: the DNMT/TET/UHRF1 core machinery should be in panel_narrow.
# If none of these show up, the category filter above needs adjusting to
# match this version's actual category strings.
core_expected <- c("DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2", "TET3", "UHRF1")
cat("\nCore DNA methylation machinery found in panel_narrow:\n")
print(intersect(core_expected, panel_narrow))
missing_core <- setdiff(core_expected, panel_narrow)
if (length(missing_core) > 0) {
  cat("MISSING from panel_narrow (adjust the is_dna_mod filter above):",
      paste(missing_core, collapse = ", "), "\n")
}


# - - - 4. Cross-reference against your own promoter panel - - -
# EpiFactors gene symbols will not perfectly match gene_mat's rownames --
# HGNC aliasing, deprecated symbols, and the multi-gene TSS collapsing in
# your promoter-aggregation step (07/08/13) all cause dropout. Report the
# overlap explicitly rather than assuming it is complete.

gene_mat_promoters_path <- "01_processed_data/ann_matched_cache.rds"
if (file.exists(gene_mat_promoters_path)) {
  # If gene_mat itself has been saved from a prior run, load it directly.
  # Otherwise this section just reports against the promoter GENE LABELS,
  # which is enough to gauge overlap without rebuilding the full matrix.
  cat("\nTo check overlap against your actual gene_mat, run after loading",
      "gene_mat from script 13/14/15 (rebuild via section 1c of those",
      "scripts), then:\n\n")
  cat("  overlap_broad  <- intersect(panel_broad,  rownames(gene_mat))\n")
  cat("  overlap_narrow <- intersect(panel_narrow, rownames(gene_mat))\n")
  cat("  cat(length(overlap_broad), 'of', length(panel_broad),",
      "'broad-panel genes have a promoter probe in your data\\n')\n")
  cat("  cat(length(overlap_narrow), 'of', length(panel_narrow),",
      "'narrow-panel genes have a promoter probe in your data\\n')\n")
}


# - - - 5. Save both panels - - -

write.csv(data.frame(gene = panel_broad),
          "01_processed_data/epifactors_panel_broad.csv", row.names = FALSE)
write.csv(data.frame(gene = panel_narrow),
          "01_processed_data/epifactors_panel_narrow.csv", row.names = FALSE)

cat("\nPanels saved to 01_processed_data/. Next step: restrict the",
    "screen in 14/15 to gene_mat[intersect(panel_narrow, rownames(gene_mat)), ]",
    "instead of the full 4,141-promoter matrix.\n")

# =============================================================================
# 19_check_cbx_qc_dropout.R
#
# Script 18 confirmed CBX1, CBX3, CBX4, CBX6, CBX8 all have promoter probes
# in the raw 450K manifest, correctly listed first in UCSC_RefGene_Name --
# so they should NOT have been dropped by the aggregation step's first-gene
# convention. That check queried the full array, though, not your actual
# working matrix. This checks whether those probes survived the upstream QC
# filtering that produced meth_m_matrix.rds -- if they were flagged there
# (detection p-value failure, SNP-affected probe, cross-reactive probe, etc.)
# they'd disappear before gene aggregation ever sees them, which would fully
# explain script 17's "missing" report without any bug in scripts 13/14/15/17.
# =============================================================================

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

m_mat <- readRDS("01_processed_data/meth_m_matrix.rds")
ann_df <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

check_genes <- c("CBX1", "CBX3", "CBX4", "CBX6", "CBX8", "CBX2", "CBX5", "CBX7")

cat("Probes in m_mat (post-QC) out of probes on the full array:\n\n")

for (g in check_genes) {
  hits <- grepl(paste0("(^|;)", g, "(;|$)"), ann_df$UCSC_RefGene_Name)
  sub <- ann_df[hits, c("UCSC_RefGene_Name", "UCSC_RefGene_Group")]
  is_promoter <- grepl("TSS200|TSS1500", sub$UCSC_RefGene_Group)
  promoter_probe_ids <- rownames(sub)[is_promoter]

  n_total <- length(promoter_probe_ids)
  n_in_mmat <- sum(promoter_probe_ids %in% rownames(m_mat))

  cat(sprintf("%-6s %d/%d promoter probes survived into m_mat%s\n",
              g, n_in_mmat, n_total,
              ifelse(n_in_mmat == 0, "  <-- QC filtering dropped ALL of them, real gap",
                     ifelse(n_in_mmat < n_total, "  <-- partial dropout", "  <-- fully retained"))))
}

cat("\nIf a gene shows 0/N survived, its promoter probes were removed during",
    "whatever QC step upstream produced meth_m_matrix.rds -- check the",
    "filtering section of 08_methylation_survival_v3.R (detection p-value,",
    "SNP/cross-reactive probe masking) for what removed them, and whether",
    "any of the standard exclusion lists (e.g. Chen et al. 2013",
    "cross-reactive probes) happen to catch this gene's probes specifically.\n")

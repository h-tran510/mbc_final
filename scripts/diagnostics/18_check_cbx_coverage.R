# =============================================================================
# 18_check_cbx_coverage.R
#
# Purpose: script 17 reported CBX1, CBX3, CBX4, CBX6, CBX8 as "missing" from
# gene_mat. There are two very different explanations for that, and they
# matter for how to read the CBX2/CBX5/CBX7 results:
#
#   (a) the 450K array genuinely has no TSS200/TSS1500 probe near these
#       genes' promoters -- a real coverage gap, nothing to fix.
#   (b) a probe DOES exist, but the aggregation step in scripts 07/08/13/17
#       only keeps the FIRST gene name when a probe is annotated to multiple
#       genes (gene_labels_clean <- sapply(strsplit(x, ";"), function(x) x[1])).
#       If a probe is annotated e.g. "SOMEGENE;CBX1" with CBX1 listed second,
#       it silently drops out of gene_mat even though real data exists.
#
# This just checks the raw 450K annotation directly, before any aggregation,
# so (a) and (b) can be told apart.
# =============================================================================

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

ann_df <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

check_genes <- c("CBX1", "CBX3", "CBX4", "CBX6", "CBX8", "CBX2", "CBX5", "CBX7")

for (g in check_genes) {
  hits <- grepl(paste0("(^|;)", g, "(;|$)"), ann_df$UCSC_RefGene_Name)
  n_probes <- sum(hits)
  if (n_probes == 0) {
    cat(sprintf("%-6s NOT in 450K manifest at all -- true coverage gap\n", g))
    next
  }
  sub <- ann_df[hits, c("UCSC_RefGene_Name", "UCSC_RefGene_Group")]
  is_promoter <- grepl("TSS200|TSS1500", sub$UCSC_RefGene_Group)
  n_promoter_probes <- sum(is_promoter)
  if (n_promoter_probes == 0) {
    cat(sprintf("%-6s %d probe(s) on the array, but none in TSS200/TSS1500 -- true promoter coverage gap\n",
                g, n_probes))
    next
  }
  first_gene_ok <- sapply(sub$UCSC_RefGene_Name[is_promoter], function(x) {
    strsplit(x, ";")[[1]][1] == g
  })
  cat(sprintf("%-6s %d promoter probe(s) exist | listed FIRST on %d/%d of them%s\n",
              g, n_promoter_probes, sum(first_gene_ok), length(first_gene_ok),
              ifelse(sum(first_gene_ok) == 0,
                     "  <-- DROPPED BY THE FIRST-GENE-ONLY CONVENTION, fixable",
                     "")))
}

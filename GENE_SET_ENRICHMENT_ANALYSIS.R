# loading packages
library(missMethyl)
library(data.table)
library(dplyr)

# identifying file names for the aging EWAS results to be analyzed
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/")
load("results_final_EWAS_HEALS_combined_mval_specifySV_log2_UrAsgmCr_9SV.RData")

# importing hallmark geneset database
hallmark <-
  readRDS(url(
    "http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds"
  ))

# setting parameters for missMethyl
sig.cpg <- (results_final %>% filter(adj.P.Val < 0.05))$Name
all.cpg <- results_final$Name
GSEA_GO_FDR <-
  data.frame(gometh(
    sig.cpg,
    all.cpg,
    collection = ("GO"),
    array.type = c("EPIC")
  ))
GSEA_KEGG_FDR <-
  data.frame(gometh(
    sig.cpg,
    all.cpg,
    collection = ("KEGG"),
    array.type = c("EPIC")
  ))
GSEA_HALLMARK_FDR <-
  gsameth(sig.cpg,
          all.cpg,
          collection = hallmark,
          array.type = c("EPIC"))
GSEA_HALLMARK_FDR <- data.frame(GSEA_HALLMARK_FDR)
GSEA_HALLMARK_FDR$Description <- rownames(GSEA_HALLMARK_FDR)

write.table(
  GSEA_GO_FDR %>% arrange(P.DE) %>% filter(P.DE < 0.05),
  file = paste0(
    "/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/GSEA/Supplementary_Table_GSEA_GO_FDR.tsv"
  ),
  row.names = F,
  sep = "\t"
)
write.table(
  GSEA_KEGG_FDR %>% arrange(P.DE) %>% filter(P.DE < 0.05),
  file = paste0(
    "/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/GSEA/Supplementary_Table_GSEA_KEGG_FDR.tsv"
  ),
  row.names = F,
  sep = "\t"
)
write.table(
  GSEA_HALLMARK_FDR %>% arrange(P.DE) %>% filter(P.DE < 0.05),
  file = paste0(
    "/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/GSEA/Supplementary_Table_GSEA_HALLMARK_FDR.tsv"
  ),
  row.names = F,
  sep = "\t"
)

print(paste(
  "Number of Significant GO terms (FDR):",
  nrow(GSEA_GO_FDR %>% filter(FDR < 0.05))
))
print(paste(
  "Number of Significant KEGG terms (FDR):",
  nrow(GSEA_KEGG_FDR %>% filter(FDR < 0.05))
))
print(paste(
  "Number of Significant Hallmark GS (FDR):",
  nrow(GSEA_HALLMARK_FDR %>% filter(FDR < 0.05))
))

print(paste(
  "Number of Significant GO terms (nominal):",
  nrow(GSEA_GO_FDR %>% filter(P.DE < 0.05))
))
print(paste(
  "Number of Significant KEGG terms (nominal):",
  nrow(GSEA_KEGG_FDR %>% filter(P.DE < 0.05))
))
print(paste(
  "Number of Significant Hallmark GS (nominal):",
  nrow(GSEA_HALLMARK_FDR %>% filter(P.DE < 0.05))
))

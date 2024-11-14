library(data.table)
library(tidyr)
library(dplyr)
library(tidyverse)
library(sva)
library(mgcv)
library(nlme)
library(genefilter)
library(BiocParallel)
library(ggplot2)
library(scales)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# reading in background chromatin segmentation file
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/download_data")
coreMarks_index=args[1]
chromatin_file_name <- paste0(coreMarks_index,"_15_coreMarks_dense.bed.gz")
chromatin_background <- fread(chromatin_file_name,header=F) %>% select(V1,V2,V3,V4)

# identifying bonferroni significant CpGs from the pooled analysis that were identified in both the log2_UrAsgmCr and log2_WArsenic EWAS
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output")
# log2_UrAsgmCr
exposure_var = "log2_UrAsgmCr"
load(paste0("results_final_EWAS_HEALS_combined_mval_specifySV_",exposure_var,"_9SV.RData"))
# log2_UrAsgmCr_significant_CpG_list<-(results_final%>%filter(P.Value< (0.05/nrow(results_final))))$Name
log2_UrAsgmCr_significant_CpG_list<-(results_final%>%filter(adj.P.Val < 0.05))$Name
results_final$chr_character <- paste0(rep("chr",nrow(results_final)), results_final$CHR)

# adding a column for chromatin segmentation
results_final$chromatin_segmentation <- NA
for (i in c(1:nrow(results_final))) {
  if (i %% 10000 == 0) {
    print(paste("CpG:", i))
  }
  tmp_result <- (chromatin_segmentation = (chromatin_background %>% filter(V1 == results_final$chr_character[i],V2 < results_final$MAPINFO[i],V3 > results_final$MAPINFO[i]))$V4)
  if (length(tmp_result) != 0) {
    results_final$chromatin_segmentation[i] <- tmp_result
  } else {
    results_final$chromatin_segmentation[i] <- NA
  }
}

# computing log odds ratios for significant vs. not significant CpGs for each chromatin segmentation annotation, as well as associated chisq.test p-values
significant_chromatin <- results_final %>% filter(Name %in% log2_UrAsgmCr_significant_CpG_list) %>% select(Name,chromatin_segmentation)
notsignificant_chromatin <- results_final %>% filter(!(Name %in% log2_UrAsgmCr_significant_CpG_list)) %>% select(Name,chromatin_segmentation)

segmentation_list_annotation <- unique(results_final$chromatin_segmentation)
segmentation_list_annotation <- segmentation_list_annotation[!is.na(segmentation_list_annotation)]
segmentation_list_logOR <- vector(length=length(segmentation_list_annotation))
segmentation_list_pval <- vector(length=length(segmentation_list_annotation))

for (i in 1:length(segmentation_list_annotation)) {
value_a <- nrow(significant_chromatin %>% filter(chromatin_segmentation==segmentation_list_annotation[i]))
value_c <- nrow(significant_chromatin %>% filter(chromatin_segmentation!=segmentation_list_annotation[i]))
value_b <- nrow(notsignificant_chromatin %>% filter(chromatin_segmentation==segmentation_list_annotation[i]))
value_d <- nrow(notsignificant_chromatin %>% filter(chromatin_segmentation!=segmentation_list_annotation[i]))

segmentation_list_logOR[i] <- log(value_a*value_d / (value_b * value_c))
segmentation_list_pval[i] <- chisq.test(matrix(ncol=2,nrow=2,c(value_a,value_c,value_b,value_d)))$p.value

print(matrix(ncol=2,nrow=2,c(value_a,value_c,value_b,value_d)))
}

chromatin_logOR_pval <- data.frame(
  chromatin_segmentation_type = segmentation_list_annotation,
  logOR = segmentation_list_logOR,
  pval = segmentation_list_pval
)
chromatin_logOR_pval$stars <- cut(chromatin_logOR_pval$pval, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 
values <- c(0.01, 0.5, 1, 5, 10, 20) %>% log
breaks <- c(0.01, 0.5, 1, 5, 10, 20) %>% log
labels <- c("0.01", "0.5", "1", "5", "10", "20")
limits <- c(0.01, 20) %>% log

chromatin_logOR_pval$type <- rep("Differentially Methylated",nrow(chromatin_logOR_pval))

setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output")
output_file_name <- paste0("chromatin_logOR_pval_",coreMarks_index,".RData")
save(chromatin_logOR_pval,file=output_file_name)
print(paste("Finished outputting", output_file_name))

# plotting the enrichment of each chromatin segmentation annotation
chromatin_logOR_pval <- chromatin_logOR_pval %>% filter(logOR!=-Inf)
g1 <- ggplot(chromatin_logOR_pval, aes(x = type, factor(chromatin_segmentation_type))) +
    geom_tile(aes(fill = logOR))+
    geom_text(aes(label=stars), color="black", size=2)+
    theme_bw()+coord_equal()+
    scale_fill_gradientn(colours=c("orangered2",
                                   "lightsalmon",
                                   "white", "lightblue", "lightblue3", "lightblue4"),
                         values=rescale(values), 
                         guide="colorbar",
                         labels = labels,
                         limits = limits, 
                         breaks=breaks)+xlab("")+ylab("Chromatin segmentation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+labs(fill = "logOR")
g1
ggsave(plot = g1, paste0("tissue_chmm_stratified_enrichment_",coreMarks_index,".png"), width = 10, height = 12, dpi = 320)














# initializing packages
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

myargs = commandArgs(trailingOnly=TRUE)
exposure_var<-myargs[1]

# importing methylation and covariate data
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data")
load("combined_cov.RData")
load("combined_beta.RData")
load("combined_mval.RData")

if (exposure_var %in% c("log2_UrAsgmCr", "log2_WArsenic")) {
    #############################################################
    # EWAS MVAL 
    fmla <- as.formula(paste0("~ combined_cov$",exposure_var," + factor(combined_cov$Sex) + combined_cov$Age + factor(combined_cov$bmi_cat) + factor(combined_cov$cigsmoke) + factor(combined_cov$Batch_Plate) + combined_cov$B + combined_cov$NK + combined_cov$CD4T + combined_cov$CD8T + combined_cov$Mono"))
    design <- model.matrix(fmla)
    library(limma)
    fit <- lmFit(combined_mval, design)
    fit <- eBayes(fit)
    results <- topTable(fit, n=dim(combined_mval)[1], sort.by="P",adjust="BH", coef=paste0("combined_cov$",exposure_var), confint=TRUE)
    save(results, file=paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/results_EWAS_HEALS_combined_mval_celltype_noSV_",exposure_var,".RData"))

    #manifest file obtain from this
    # wget https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip
    # unzip infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip 
    manifest<-read.csv("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/annotation/manifest/manifest_4.csv", header=T,sep=",")
    results$Name<-rownames(results)
    results_annotated <- merge(results, manifest, by="Name")
    results_final <- results_annotated[order(results_annotated$P.Value), ]
    output_file<-paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/results_final_EWAS_HEALS_combined_mval_celltype_noSV_",exposure_var,".RData")
    save(results_final,file=output_file)


    #############################################################
    # EWAS BETA 
    #running limma
    fmla <- as.formula(paste0("~ combined_cov$",exposure_var," + factor(combined_cov$Sex) + combined_cov$Age + factor(combined_cov$bmi_cat) + factor(combined_cov$cigsmoke) + factor(combined_cov$Batch_Plate) + combined_cov$B + combined_cov$NK + combined_cov$CD4T + combined_cov$CD8T + combined_cov$Mono"))
    design <- model.matrix(fmla)
    library(limma)
    fit <- lmFit(combined_beta, design)
    fit <- eBayes(fit)
    results <- topTable(fit, n=dim(combined_beta)[1], sort.by="P",adjust="BH", coef=paste0("combined_cov$",exposure_var), confint=TRUE)
    save(results, file=paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/results_EWAS_HEALS_combined_beta_celltype_noSV_",exposure_var,".RData"))

    #manifest file obtain from this
    # wget https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip
    # unzip infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip 
    manifest<-read.csv("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/annotation/manifest/manifest_4.csv", header=T,sep=",")
    results$Name<-rownames(results)
    results_annotated <- merge(results, manifest, by="Name")
    results_final <- results_annotated[order(results_annotated$P.Value), ]
    output_file<-paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/results_final_EWAS_HEALS_combined_beta_celltype_noSV_",exposure_var,".RData")
    save(results_final,file=output_file)

} else {
    print("###### ERROR ######")
}

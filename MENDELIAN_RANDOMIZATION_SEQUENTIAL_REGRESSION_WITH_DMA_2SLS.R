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

args = commandArgs(trailingOnly=TRUE)
pval_threshold_method = args[1]

# identifying bonferroni significant CpGs from the pooled analysis that were identified in both the log2_UrAsgmCr and log2_WArsenic EWAS
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output")
# log2_WArsenic
exposure_var = "log2_WArsenic"
load(paste0("results_final_EWAS_HEALS_combined_mval_specifySV_",exposure_var,"_9SV.RData"))
log2_WArsenic_significant_CpG_list<-(results_final%>%filter(adj.P.Val< 0.05))$Name

# log2_UrAsgmCr
exposure_var = "log2_UrAsgmCr"
load(paste0("results_final_EWAS_HEALS_combined_mval_specifySV_",exposure_var,"_9SV.RData"))

if (pval_threshold_method=="bonferroni") {
  log2_UrAsgmCr_significant_CpG_list<-(results_final%>%filter(P.Value< (0.05/nrow(results_final))))$Name
} else if (pval_threshold_method=="fdr") { 
  log2_UrAsgmCr_significant_CpG_list<-(results_final%>%filter(adj.P.Val< 0.05))$Name
} else {
  print("error buddy!")
}

significant_CpG_list <- intersect(log2_UrAsgmCr_significant_CpG_list,log2_WArsenic_significant_CpG_list)
# need the latest results to be for log2_UrAsgmCr prior to running the below command 
significant_CpG_table <- results_final %>% filter(Name %in% significant_CpG_list)
saveRDS(significant_CpG_list, file = paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/significant_CpG_list_",pval_threshold_method,".rds"))

# loading in covariate and mval matrix files 
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data")
load("combined_cov.RData")
load("combined_mval.RData")

# obtaining list of individuals 
methylation_individuals <- combined_cov$Sample_Name

# setting up a combined data table containing CpG mvals (not corrected for batch) and covariates for CpGs that we're interested in
MR_CpGs<-combined_mval[(rownames(combined_mval) %in% significant_CpG_list),] 
t_MR_CpGs<-data.frame(t(MR_CpGs))
t_MR_CpGs$Sample_Name<-rownames(t_MR_CpGs)
MR_covariates_mval<-inner_join(t_MR_CpGs,combined_cov,by=c("Sample_Name"))

# combining the table with covariates and mvals for each sample with the SNP dosage / GP-DMA% / binary high efficiency SNP alleles table
SNP_dosages<-fread("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/snp_data/HEALS_chr10.104616500_chr10.104747534_chr10.104635687_chr21.47572887_DMAperc_PRS_DMApredicted_10.24.22.csv",header=T,sep=",")
SNP_dosages<-na.omit(SNP_dosages) 
SNP_dosages$Sample_Name<-paste0("X",SNP_dosages$IID)
SNP_dosages <- SNP_dosages %>% filter(Sample_Name %in% methylation_individuals)
SNP_dosages$DMAperc_predicted_OLD <- SNP_dosages$DMAperc_predicted
SNP_dosages$DMAperc_predicted <- NULL
# calculating GP-DMA% (comeback 2 SNP or 4 SNP)
linear_model <- lm(DMA_perc~
                      X10.104616500_C +
                      #X10.104747534_G +
                      #X10.104635687_T +
                      X21.47572887_T, data = SNP_dosages)
SNP_dosages$DMAperc_predicted <- predict(linear_model, newdata = SNP_dosages)
SNP_dosages <- SNP_dosages %>% select(-DMAperc_predicted_OLD)

# creating a binary SNP score for whether someone only has high efficiency alleles
SNP_dosages$binary_score <- rep(0,nrow(SNP_dosages))
SNP_dosages$binary_score[((SNP_dosages$X10.104616500_C==0) & (SNP_dosages$X21.47572887_T==0))] <- 1

MR_covariates_mval_SNPdosage<-inner_join(MR_covariates_mval,SNP_dosages,by=c("Sample_Name"))

MR_covariates_mval_SNPdosage<-MR_covariates_mval_SNPdosage[which(!is.na(MR_covariates_mval_SNPdosage$DMAperc_predicted)),]
save(MR_covariates_mval_SNPdosage, file="/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/MR_covariates_mval_SNPdosage_significantCpGs.RData")

#######################################
# MR start
load("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/MR_covariates_mval_SNPdosage_significantCpGs.RData")

#######################################
#######################################
#######################################
# regressing mval on DMA% 
last_CpG_col<-tail(which(grepl("cg",names(MR_covariates_mval_SNPdosage))),1)

# covariates to adjust for:
# sex, age, methylation batch+plate, and cigarette smoking status

cpg_name_vec<-vector(length=last_CpG_col)
coef_vec<-vector(length=last_CpG_col)
se_vec<-vector(length=last_CpG_col)
pval_vec<-vector(length=last_CpG_col)

# performing regression of CpG mval on DMA% for all our CpGs
for (i in 1:last_CpG_col) {
  # fitting regression 
  fmla <- as.formula(paste0(
    "MR_covariates_mval_SNPdosage[,i] ~ MR_covariates_mval_SNPdosage$DMA_perc + as.factor(MR_covariates_mval_SNPdosage$Sex) + MR_covariates_mval_SNPdosage$Age + as.factor(MR_covariates_mval_SNPdosage$cigsmoke) + as.factor(MR_covariates_mval_SNPdosage$Batch_Plate) + MR_covariates_mval_SNPdosage$B + MR_covariates_mval_SNPdosage$NK + MR_covariates_mval_SNPdosage$CD4T + MR_covariates_mval_SNPdosage$CD8T + MR_covariates_mval_SNPdosage$Mono")
  )
  fit_OLS<-(
    lm(
      fmla
    )
  )
  
  # determining which coefficient in the output corresponds to our desired p-value of DMA% as well the index of the standard error of this coefficient
  pval_index_reg_output<-((ncol(summary(fit_OLS)$coefficients)-1)*(nrow(summary(fit_OLS)$coefficients)))+2
  se_index_reg_output<-((ncol(summary(fit_OLS)$coefficients)-3)*(nrow(summary(fit_OLS)$coefficients)))+2
  coef_index_reg_output<-((ncol(summary(fit_OLS)$coefficients)-4)*(nrow(summary(fit_OLS)$coefficients)))+2
  
  #print(summary(fit_OLS)$coefficients[coef_index_reg_output])
  #print(summary(fit_OLS)$coefficients[se_index_reg_output])
  #print(summary(fit_OLS)$coefficients[pval_index_reg_output])
  
  # storing regression results of DMA% in vectors
  cpg_name_vec[i]<-colnames(MR_covariates_mval_SNPdosage)[i]
  coef_vec[i]<-summary(fit_OLS)$coefficients[coef_index_reg_output]
  se_vec[i]<-summary(fit_OLS)$coefficients[se_index_reg_output]
  pval_vec[i]<-summary(fit_OLS)$coefficients[pval_index_reg_output]
  
  print(i)
}


# constructing a data table for our MR estimates
MR_dt<-data.table(cpg_name_vec,coef_vec,se_vec,pval_vec)
colnames(MR_dt)[1]<-"Name"
MR_dt$direction[MR_dt$coef_vec>1] <- 1
MR_dt$direction[MR_dt$coef_vec<1] <- -1
MR_dt$sample_size<-nrow(MR_covariates_mval_SNPdosage)
save(MR_dt,file=paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/MR/",pval_threshold_method,"/MR_DMA_MR_dt.RData"))

# plotting histogram of p-values
setwd(paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/MR/",pval_threshold_method))
png("MR_DMA_pval_hist.png")
hist(pval_vec, main="Histogram of p-values of the association between \nDMA% and each CpG (M-value)", xlab="p-value", breaks=30)
dev.off()
# outputting a summary of the distribution of p-values
print(summary(pval_vec))
print(summary(p.adjust(pval_vec,method="bonferroni")))
print(summary(p.adjust(pval_vec,method="fdr")))

# merging this MR estimate data table with effect size estimates from the EWAS
MR_dt_EWAS<-inner_join(MR_dt,significant_CpG_table,by=c("Name"))
MR_dt_EWAS$Consistency<-"Inconsistent"
MR_dt_EWAS$Consistency[(MR_dt_EWAS$coef_vec>0 & MR_dt_EWAS$logFC<0) | (MR_dt_EWAS$coef_vec<0 & MR_dt_EWAS$logFC>0)]<-"Consistent"
MR_dt_EWAS<-MR_dt_EWAS[order(MR_dt_EWAS$coef_vec),]
MR_dt_EWAS$index<-1:nrow(MR_dt_EWAS)

# performing binomial test for consistency
binom.test.summary <- binom.test(x=nrow(MR_dt_EWAS%>%filter(Consistency=="Consistent")), n=nrow(MR_dt_EWAS), p = 0.5,
                                 alternative = c("greater"),
                                 conf.level = 0.95)
plot_p_str <- paste0("Binomial p-value: ", formatC(binom.test.summary$p.value, format = "e", digits = 2))
plot_prop_str <- paste0("Directional Consistency of CpGs: ", paste(round(100*binom.test.summary$estimate, 2), "%", sep=""))

annotation_p <- annotation_custom(grid::textGrob(label = plot_p_str,x = unit(0.35, "npc"), y = unit(0.8, "npc"),gp = grid::gpar(cex = 1)))
annotation_prop <- annotation_custom(grid::textGrob(label = plot_prop_str,x = unit(0.35, "npc"), y = unit(0.9, "npc"),gp = grid::gpar(cex = 1)))

# making an error bar plot of these MR estimates
png("MR_DMA_error_bar_plot.png",units="in",height=4,width=5.75,res=1800)
p<- ggplot(MR_dt_EWAS, aes(x=index, y=coef_vec, group=Consistency, color=Consistency)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=coef_vec-se_vec, ymax=coef_vec+se_vec), width=.2,
                position=position_dodge(0.05)) + 
  #ggtitle("Error bar plot of the effect of DMA% \non CpG methylation (M-value)") + 
  xlab("CpG Site Ordered by Association Estimate") + 
  #ylab(expression(atop("Association of measured AME (DMA%) on", paste("DNAm at As-associated CpG sites")))) +
  ylab("Association of measured AME (DMA%)") +
  theme_classic() + theme(legend.position = c(0.75,0.2), panel.border = element_rect(colour = "black", fill=NA)) + 
  annotation_p + annotation_prop + geom_hline(yintercept = 0,linetype="dotted")
p
dev.off()

#######################################
#######################################
#######################################
# regressing mval on GP-DMA% 
last_CpG_col<-tail(which(grepl("cg",names(MR_covariates_mval_SNPdosage))),1)

# covariates to adjust for:
# sex, age, methylation batch+plate, and cigarette smoking status

cpg_name_vec<-vector(length=last_CpG_col)
coef_vec<-vector(length=last_CpG_col)
se_vec<-vector(length=last_CpG_col)
pval_vec<-vector(length=last_CpG_col)

# performing regression of CpG mval on GP-DMA% for all our CpGs
for (i in 1:last_CpG_col) {
  # fitting regression 
  fmla <- as.formula(paste0(
    "MR_covariates_mval_SNPdosage[,i] ~ MR_covariates_mval_SNPdosage$DMAperc_predicted + as.factor(MR_covariates_mval_SNPdosage$Sex) + MR_covariates_mval_SNPdosage$Age + as.factor(MR_covariates_mval_SNPdosage$cigsmoke) + as.factor(MR_covariates_mval_SNPdosage$Batch_Plate) + MR_covariates_mval_SNPdosage$B + MR_covariates_mval_SNPdosage$NK + MR_covariates_mval_SNPdosage$CD4T + MR_covariates_mval_SNPdosage$CD8T + MR_covariates_mval_SNPdosage$Mono")
  )
  fit_OLS<-(
    lm(
      fmla
    )
  )
  
  # determining which coefficient in the output corresponds to our desired p-value of GP-DMA% as well the index of the standard error of this coefficient
  pval_index_reg_output<-((ncol(summary(fit_OLS)$coefficients)-1)*(nrow(summary(fit_OLS)$coefficients)))+2
  se_index_reg_output<-((ncol(summary(fit_OLS)$coefficients)-3)*(nrow(summary(fit_OLS)$coefficients)))+2
  coef_index_reg_output<-((ncol(summary(fit_OLS)$coefficients)-4)*(nrow(summary(fit_OLS)$coefficients)))+2
  
  #print(summary(fit_OLS)$coefficients[coef_index_reg_output])
  #print(summary(fit_OLS)$coefficients[se_index_reg_output])
  #print(summary(fit_OLS)$coefficients[pval_index_reg_output])
  
  # storing regression results of GP-DMA% in vectors
  cpg_name_vec[i]<-colnames(MR_covariates_mval_SNPdosage)[i]
  coef_vec[i]<-summary(fit_OLS)$coefficients[coef_index_reg_output]
  se_vec[i]<-summary(fit_OLS)$coefficients[se_index_reg_output]
  pval_vec[i]<-summary(fit_OLS)$coefficients[pval_index_reg_output]
  
  print(i)
}


# constructing a data table for our MR estimates
MR_dt<-data.table(cpg_name_vec,coef_vec,se_vec,pval_vec)
colnames(MR_dt)[1]<-"Name"
MR_dt$direction[MR_dt$coef_vec>1] <- 1
MR_dt$direction[MR_dt$coef_vec<1] <- -1
MR_dt$sample_size<-nrow(MR_covariates_mval_SNPdosage)
save(MR_dt,file=paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/MR/",pval_threshold_method,"/MR_GPDMA_MR_dt.RData"))

# plotting histogram of p-values
setwd(paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/MR/",pval_threshold_method))
png("MR_GPDMA_pval_hist.png")
hist(pval_vec, main="Histogram of p-values of the association between \nGP-DMA% and each CpG (M-value)", xlab="p-value", breaks=30)
dev.off()
# outputting a summary of the distribution of p-values
print(summary(pval_vec))
print(summary(p.adjust(pval_vec,method="bonferroni")))
print(summary(p.adjust(pval_vec,method="fdr")))

# merging this MR estimate data table with effect size estimates from the EWAS
MR_dt_EWAS<-inner_join(MR_dt,significant_CpG_table,by=c("Name"))
MR_dt_EWAS$Consistency<-"Inconsistent"
MR_dt_EWAS$Consistency[(MR_dt_EWAS$coef_vec>0 & MR_dt_EWAS$logFC<0) | (MR_dt_EWAS$coef_vec<0 & MR_dt_EWAS$logFC>0)]<-"Consistent"
MR_dt_EWAS<-MR_dt_EWAS[order(MR_dt_EWAS$coef_vec),]
MR_dt_EWAS$index<-1:nrow(MR_dt_EWAS)

# performing binomial test for consistency
binom.test.summary <- binom.test(x=nrow(MR_dt_EWAS%>%filter(Consistency=="Consistent")), n=nrow(MR_dt_EWAS), p = 0.5,
                                 alternative = c("greater"),
                                 conf.level = 0.95)
plot_p_str <- paste0("Binomial p-value: ", formatC(binom.test.summary$p.value, format = "e", digits = 2))
plot_prop_str <- paste0("Directional Consistency of CpGs: ", paste(round(100*binom.test.summary$estimate, 2), "%", sep=""))

annotation_p <- annotation_custom(grid::textGrob(label = plot_p_str,x = unit(0.35, "npc"), y = unit(0.8, "npc"),gp = grid::gpar(cex = 1)))
annotation_prop <- annotation_custom(grid::textGrob(label = plot_prop_str,x = unit(0.35, "npc"), y = unit(0.9, "npc"),gp = grid::gpar(cex = 1)))

# making an error bar plot of these MR estimates
png("MR_GPDMA_error_bar_plot.png",units="in",height=4,width=5.75,res=1800)
p<- ggplot(MR_dt_EWAS, aes(x=index, y=coef_vec, group=Consistency, color=Consistency)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=coef_vec-se_vec, ymax=coef_vec+se_vec), width=.2,
                position=position_dodge(0.05)) + 
  #ggtitle("Error bar plot of the effect of GP-DMA% \non CpG methylation (M-value)") + 
  xlab("CpG Site Ordered by Association Estimate") + 
  ylab(expression("Weighted GP-DMA% ("*beta*")")) +
  theme_classic() + theme(legend.position = c(0.75,0.2), panel.border = element_rect(colour = "black", fill=NA)) +
  annotation_p + annotation_prop + geom_hline(yintercept = 0,linetype="dotted")
p
dev.off()

#######################################
#######################################
#######################################
# regressing mval on Binary SNP Score 
last_CpG_col<-tail(which(grepl("cg",names(MR_covariates_mval_SNPdosage))),1)

# covariates to adjust for:
# sex, age, methylation batch, and cigarette smoking status

cpg_name_vec<-vector(length=last_CpG_col)
coef_vec<-vector(length=last_CpG_col)
se_vec<-vector(length=last_CpG_col)
pval_vec<-vector(length=last_CpG_col)

# performing regression of CpG mval on Binary SNP Score for all our CpGs
for (i in 1:last_CpG_col) {
  # fitting regression 
  fmla <- as.formula(paste0(
    "MR_covariates_mval_SNPdosage[,i] ~ MR_covariates_mval_SNPdosage$binary_score + as.factor(MR_covariates_mval_SNPdosage$Sex) + MR_covariates_mval_SNPdosage$Age + as.factor(MR_covariates_mval_SNPdosage$cigsmoke) + as.factor(MR_covariates_mval_SNPdosage$Batch_Plate) + MR_covariates_mval_SNPdosage$B + MR_covariates_mval_SNPdosage$NK + MR_covariates_mval_SNPdosage$CD4T + MR_covariates_mval_SNPdosage$CD8T + MR_covariates_mval_SNPdosage$Mono")
  )
  fit_OLS<-(
    lm(
      fmla
    )
  )
  
  # determining which coefficient in the output corresponds to our desired p-value of Binary SNP Score as well the index of the standard error of this coefficient
  pval_index_reg_output<-((ncol(summary(fit_OLS)$coefficients)-1)*(nrow(summary(fit_OLS)$coefficients)))+2
  se_index_reg_output<-((ncol(summary(fit_OLS)$coefficients)-3)*(nrow(summary(fit_OLS)$coefficients)))+2
  coef_index_reg_output<-((ncol(summary(fit_OLS)$coefficients)-4)*(nrow(summary(fit_OLS)$coefficients)))+2
  
  #print(summary(fit_OLS)$coefficients[coef_index_reg_output])
  #print(summary(fit_OLS)$coefficients[se_index_reg_output])
  #print(summary(fit_OLS)$coefficients[pval_index_reg_output])
  
  # storing regression results of Binary SNP Score in vectors
  cpg_name_vec[i]<-colnames(MR_covariates_mval_SNPdosage)[i]
  coef_vec[i]<-summary(fit_OLS)$coefficients[coef_index_reg_output]
  se_vec[i]<-summary(fit_OLS)$coefficients[se_index_reg_output]
  pval_vec[i]<-summary(fit_OLS)$coefficients[pval_index_reg_output]
  
  print(i)
}


# constructing a data table for our MR estimates
MR_dt<-data.table(cpg_name_vec,coef_vec,se_vec,pval_vec)
colnames(MR_dt)[1]<-"Name"
MR_dt$direction[MR_dt$coef_vec>1] <- 1
MR_dt$direction[MR_dt$coef_vec<1] <- -1
MR_dt$sample_size<-nrow(MR_covariates_mval_SNPdosage)
save(MR_dt,file=paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/MR/",pval_threshold_method,"/MR_BINSCORE_MR_dt.RData"))

# plotting histogram of p-values
setwd(paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/MR/",pval_threshold_method))
png("BINSCORE_pval_hist.png")
hist(pval_vec, main="Histogram of p-values of the association between \nBinary SNP Score and each CpG (M-value)", xlab="p-value", breaks=30)
dev.off()
# outputting a summary of the distribution of p-values
print(summary(pval_vec))
print(summary(p.adjust(pval_vec,method="bonferroni")))
print(summary(p.adjust(pval_vec,method="fdr")))

# merging this MR estimate data table with effect size estimates from the EWAS
MR_dt_EWAS<-inner_join(MR_dt,significant_CpG_table,by=c("Name"))
MR_dt_EWAS$Consistency<-"Inconsistent"
MR_dt_EWAS$Consistency[(MR_dt_EWAS$coef_vec>0 & MR_dt_EWAS$logFC<0) | (MR_dt_EWAS$coef_vec<0 & MR_dt_EWAS$logFC>0)]<-"Consistent"
MR_dt_EWAS<-MR_dt_EWAS[order(MR_dt_EWAS$coef_vec),]
MR_dt_EWAS$index<-1:nrow(MR_dt_EWAS)

# performing binomial test for consistency
binom.test.summary <- binom.test(x=nrow(MR_dt_EWAS%>%filter(Consistency=="Consistent")), n=nrow(MR_dt_EWAS), p = 0.5,
                                 alternative = c("greater"),
                                 conf.level = 0.95)
plot_p_str <- paste0("Binomial p-value: ", formatC(binom.test.summary$p.value, format = "e", digits = 2))
plot_prop_str <- paste0("Directional Consistency of CpGs: ", paste(round(100*binom.test.summary$estimate, 2), "%", sep=""))

annotation_p <- annotation_custom(grid::textGrob(label = plot_p_str,x = unit(0.35, "npc"), y = unit(0.8, "npc"),gp = grid::gpar(cex = 1)))
annotation_prop <- annotation_custom(grid::textGrob(label = plot_prop_str,x = unit(0.35, "npc"), y = unit(0.9, "npc"),gp = grid::gpar(cex = 1)))

# making an error bar plot of these MR estimates
png("MR_BINSCORE_error_bar_plot.png",units="in",height=4,width=5.75,res=1800)
p<- ggplot(MR_dt_EWAS, aes(x=index, y=coef_vec, group=Consistency, color=Consistency)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=coef_vec-se_vec, ymax=coef_vec+se_vec), width=.2,
                position=position_dodge(0.05)) + 
  #ggtitle("Error bar plot of the effect of Binary SNP Score \non CpG methylation (M-value)") + 
  xlab("CpG Site Ordered by Association Estimate") + 
  #ylab(expression(atop("Association of genetically determined AME on", paste("DNAm at As-associated CpG sites")))) +
  ylab("Association of genetically-determined AME") +
  theme_classic() + theme(legend.position = c(0.75,0.2), panel.border = element_rect(colour = "black", fill=NA)) +
  annotation_p + annotation_prop + geom_hline(yintercept = 0,linetype="dotted")
p
dev.off()

# examining magnitude of effect sizes vs. consistency
magnitude_df <- MR_dt_EWAS %>% select(coef_vec,Consistency) %>% mutate(magnitude = abs(coef_vec))
t.test(magnitude ~ Consistency, data = magnitude_df)


#######################################
#######################################
#######################################
library(AER)

# regressing mval on Binary SNP Score 
last_CpG_col<-tail(which(grepl("cg",names(MR_covariates_mval_SNPdosage))),1)

# covariates to adjust for:
# sex, age, methylation batch, and cigarette smoking status

cpg_name_vec<-vector(length=last_CpG_col)
coef_vec<-vector(length=last_CpG_col)
se_vec<-vector(length=last_CpG_col)
pval_vec<-vector(length=last_CpG_col)

# performing regression of CpG mval on Binary SNP Score for all our CpGs
for (i in 1:last_CpG_col) {
  # fitting 2SLS for IV analysis
  current_cpg_col_name <- colnames(MR_covariates_mval_SNPdosage)[i]
  
  fmla <- as.formula(paste0(
    current_cpg_col_name," ~ DMA_perc+as.factor(Sex)+Age+as.factor(cigsmoke)+as.factor(Batch_Plate)+B+NK+CD4T+CD8T+Mono | X10.104616500_C+X21.47572887_T+as.factor(Sex)+Age+as.factor(cigsmoke)+as.factor(Batch_Plate)+B+NK+CD4T+CD8T+Mono"
  ))
  
  fit_ivreg <- summary(ivreg(formula=fmla, data = MR_covariates_mval_SNPdosage))$coefficients
  
  # storing regression results of Binary SNP Score in vectors
  cpg_name_vec[i]<-colnames(MR_covariates_mval_SNPdosage)[i]
  coef_vec[i]<-fit_ivreg["DMA_perc","Estimate"]
  se_vec[i]<-fit_ivreg["DMA_perc","Std. Error"]
  pval_vec[i]<-fit_ivreg["DMA_perc","Pr(>|t|)"]
  print(i)
}

# constructing a data table for our MR estimates
MR_dt<-data.table(cpg_name_vec,coef_vec,se_vec,pval_vec)
colnames(MR_dt)[1]<-"Name"
MR_dt$direction[MR_dt$coef_vec>1] <- 1
MR_dt$direction[MR_dt$coef_vec<1] <- -1
MR_dt$sample_size<-nrow(MR_covariates_mval_SNPdosage)
save(MR_dt,file=paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/MR/",pval_threshold_method,"/MR_2SLS_MR_dt.RData"))

# plotting histogram of p-values
setwd(paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/MR/",pval_threshold_method))
png("2SLS_pval_hist.png")
hist(pval_vec, main="Histogram of p-values of the association between \nBinary SNP Score and each CpG (M-value)", xlab="p-value", breaks=30)
dev.off()
# outputting a summary of the distribution of p-values
print(summary(pval_vec))
print(summary(p.adjust(pval_vec,method="bonferroni")))
print(summary(p.adjust(pval_vec,method="fdr")))

# merging this MR estimate data table with effect size estimates from the EWAS
MR_dt_EWAS<-inner_join(MR_dt,significant_CpG_table,by=c("Name"))
MR_dt_EWAS$Consistency<-"Inconsistent"
MR_dt_EWAS$Consistency[(MR_dt_EWAS$coef_vec>0 & MR_dt_EWAS$logFC<0) | (MR_dt_EWAS$coef_vec<0 & MR_dt_EWAS$logFC>0)]<-"Consistent"
MR_dt_EWAS<-MR_dt_EWAS[order(MR_dt_EWAS$coef_vec),]
MR_dt_EWAS$index<-1:nrow(MR_dt_EWAS)

# performing binomial test for consistency
binom.test.summary <- binom.test(x=nrow(MR_dt_EWAS%>%filter(Consistency=="Consistent")), n=nrow(MR_dt_EWAS), p = 0.5,
                                 alternative = c("greater"),
                                 conf.level = 0.95)
plot_p_str <- paste0("Binomial p-value: ", formatC(binom.test.summary$p.value, format = "e", digits = 2))
plot_prop_str <- paste0("Directional Consistency of CpGs: ", paste(round(100*binom.test.summary$estimate, 2), "%", sep=""))

annotation_p <- annotation_custom(grid::textGrob(label = plot_p_str,x = unit(0.35, "npc"), y = unit(0.8, "npc"),gp = grid::gpar(cex = 1)))
annotation_prop <- annotation_custom(grid::textGrob(label = plot_prop_str,x = unit(0.35, "npc"), y = unit(0.9, "npc"),gp = grid::gpar(cex = 1)))

# making an error bar plot of these MR estimates
png("MR_2SLS_error_bar_plot.png",units="in",height=4,width=5.75,res=1800)
p<- ggplot(MR_dt_EWAS, aes(x=index, y=coef_vec, group=Consistency, color=Consistency)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=coef_vec-se_vec, ymax=coef_vec+se_vec), width=.2,
                position=position_dodge(0.05)) + 
  xlab("CpG Site Ordered by Association Estimate") + 
  ylab("Association of genetically-predicted AME") +
  theme_classic() + theme(legend.position = c(0.75,0.2), panel.border = element_rect(colour = "black", fill=NA)) +
  annotation_p + annotation_prop + geom_hline(yintercept = 0,linetype="dotted")
p
dev.off()

# examining magnitude of effect sizes vs. consistency
magnitude_df <- MR_dt_EWAS %>% select(coef_vec,Consistency) %>% mutate(magnitude = abs(coef_vec))
t.test(magnitude ~ Consistency, data = magnitude_df)















































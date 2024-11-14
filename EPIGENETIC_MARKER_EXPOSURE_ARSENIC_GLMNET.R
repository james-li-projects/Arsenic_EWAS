#########################
# INITIALIZING PACKAGES #
#########################
library(data.table)
library(tidyr)
library(tidyverse)
library(sva)
library(mgcv)
library(nlme)
library(genefilter)
library(BiocParallel)
library(ggplot2)
library(glmnet)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

######################################
# INITIALIZING VALUES FOR PARAMETERS #
######################################
sbatch_index = args[1]
set.seed(sbatch_index)
exposure_var<-"log2_UrAsgmCr"
specifySV<-9

###############################
# IMPORTING CPG MANIFEST FILE #
###############################
# annotating chromosome and position
manifest<-read.csv("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/annotation/manifest/manifest_4.csv", header=T,sep=",") 
manifest <- manifest %>% select(Name,CHR,MAPINFO)

############################################
# IMPORTING METHYLATION AND COVARIATE DATA #
############################################
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data")
load("combined_cov.RData")
load("combined_mval.RData")
# get rid of missing BMI individuals - comeback
# #combined_cov <- combined_cov %>% filter(bmi_cat != "Unknown")
# #combined_mval <- combined_mval[,combined_cov$Sample_Name]
# recoding NA skin lesion cases as missing 
combined_cov$anyskinlesion[is.na(combined_cov$anyskinlesion)] <- "missing"

#############################################################################################################
# SUBSETTING OUR DATA INTO TRAINING AND TESTING SETS [TESTING SET HAS A 1:2 RATIO OF SL CASES AND CONTROLS] #
#############################################################################################################

##############
# COVARIATES #
##############
# separating SL cases from testing/training set
SL_cases_Sample_Name <- (combined_cov %>% filter(anyskinlesion == 1))$Sample_Name
validation_data_SL_cases_cov <- combined_cov %>% filter(Sample_Name %in% SL_cases_Sample_Name)
combined_cov_minusSL <- combined_cov %>% filter(!(Sample_Name %in% SL_cases_Sample_Name))

# identifying individuals with missing SL data and assigning their covariates data for the model training set
training_missingSL_cov <- combined_cov_minusSL %>% filter(anyskinlesion == "missing")
combined_cov_minusSL <- combined_cov_minusSL %>% filter(anyskinlesion != "missing")

# subsetting remaining covariate data into training set and remaining samples for testing/validation
picked = sample(seq_len(nrow(combined_cov_minusSL)),size = 802+186)
training_testing_data_cov <- combined_cov_minusSL[picked,]
training_testing_data_cov <- rbind(training_testing_data_cov, training_missingSL_cov)
validation_data_SL_controls_cov <- combined_cov_minusSL[-picked,]

# sorting the rows in each covariate dataset by Sample_Name
validation_data_SL_cases_cov <- validation_data_SL_cases_cov[order((validation_data_SL_cases_cov$Sample_Name)),]
training_testing_data_cov <- training_testing_data_cov[order((training_testing_data_cov$Sample_Name)),]
validation_data_SL_controls_cov <- validation_data_SL_controls_cov[order((validation_data_SL_controls_cov$Sample_Name)),]

##############
# DNAm MVALS #
##############
# subsetting mval matrices 
training_testing_data_mval <- combined_mval[,colnames(combined_mval) %in% training_testing_data_cov$Sample_Name]
validation_data_SL_cases_mval <- combined_mval[,colnames(combined_mval) %in% validation_data_SL_cases_cov$Sample_Name]
validation_data_SL_controls_mval <- combined_mval[,colnames(combined_mval) %in% validation_data_SL_controls_cov$Sample_Name]

# sorting columns by Sample_Name
training_testing_data_mval <- training_testing_data_mval[, order((colnames(training_testing_data_mval)))]
validation_data_SL_cases_mval <- validation_data_SL_cases_mval[, order((colnames(validation_data_SL_cases_mval)))]
validation_data_SL_controls_mval <- validation_data_SL_controls_mval[, order((colnames(validation_data_SL_controls_mval)))]

# checking sorting of Sample_Name in both cov and mval files
identical(colnames(training_testing_data_mval),training_testing_data_cov$Sample_Name)
identical(colnames(validation_data_SL_cases_mval),validation_data_SL_cases_cov$Sample_Name)
identical(colnames(validation_data_SL_controls_mval),validation_data_SL_controls_cov$Sample_Name)

#############################################################
# Performing EWAS with M-values on the training_testing set #
#############################################################
load(paste0("sv_EWAS_HEALS_combined_svaSV_",exposure_var,".RData"))
V <- sv[,1:specifySV]
V <- V[rownames(V) %in% training_testing_data_cov$Sample_Name,]
xnam <- paste(" V[,", 1:dim(V)[2], "]", sep="")
fmla <- as.formula(paste0("~ training_testing_data_cov$",exposure_var," + factor(training_testing_data_cov$Sex) + training_testing_data_cov$Age + factor(training_testing_data_cov$bmi_cat) + factor(training_testing_data_cov$cigsmoke) + factor(training_testing_data_cov$Batch_Plate) + ", paste(xnam, collapse= "+")))
design <- model.matrix(fmla)
library(limma)
fit <- lmFit(training_testing_data_mval, design)
fit <- eBayes(fit)
results <- topTable(fit, n=dim(training_testing_data_mval)[1], sort.by="P",adjust="BH", coef=paste0("training_testing_data_cov$",exposure_var), confint=TRUE)
results$Name <- rownames(results)
results_annotated <- inner_join(results, manifest, by=c("Name"))
results_final <- results_annotated %>% arrange(Name) 
colnames(results_final)[which(grepl("CHR",colnames(results_final)))] <- "#CHROM"
colnames(results_final)[which(grepl("MAPINFO",colnames(results_final)))] <- "POS"
results_final$`#CHROM` <- as.numeric(results_final$`#CHROM`)
results_final <- results_final %>% arrange(`#CHROM`,POS)

# filtering for CpGs that pass an FDR of 0.05
significant_results_final <- results_final %>% filter(adj.P.Val < 0.05)
training_testing_data_mval <- training_testing_data_mval[significant_results_final$Name,]

##################################################
# performing 5-fold cross-validation with glmnet #
##################################################
# setting up our glmnet inputs
t_training_testing_data_mval <- data.frame(t(training_testing_data_mval))
t_training_testing_data_mval$Sample_Name <- rownames(t_training_testing_data_mval)
joined_training_testing_df <- inner_join(training_testing_data_cov %>% select(Sample_Name,log2_UrAsgmCr),t_training_testing_data_mval,by=c("Sample_Name")) 
joined_training_testing_df_NO_SAMPLE_NAME <- joined_training_testing_df %>% select(-Sample_Name)
# previewing this joined table
joined_training_testing_df_NO_SAMPLE_NAME[1,1:20]
# running glmnet
Y=joined_training_testing_df_NO_SAMPLE_NAME[,1]
X=as.matrix(joined_training_testing_df_NO_SAMPLE_NAME[,-1])
cvfit<-cv.glmnet(x=X,y=Y,nfolds=5)
# obtaining coefficients of final model
cvfit$lambda.min
cvfit_coef <- data.frame(as.matrix(coef(cvfit, s = "lambda.min")))
cvfit_coef$Name <- rownames(cvfit_coef)
cvfit_coef <- cvfit_coef %>% filter(s1 != 0) %>% mutate(logFC=s1) %>% select(Name,logFC)
cvfit_coef <- cvfit_coef[-1,]

####################################################################################
# ASSESSING OUR EPIGENETIC MARKER OF URINARY ARSENIC LEVELS IN THE VALIDATION SET #
####################################################################################

# combining data for SL cases and controls 
validation_data_mval <- cbind(validation_data_SL_cases_mval,validation_data_SL_controls_mval)
validation_data_cov <- rbind(validation_data_SL_cases_cov,validation_data_SL_controls_cov)
# sorting validation data mval by CpG name according to the cvfit_coef file 
validation_data_mval <- validation_data_mval[rownames(cvfit_coef),]
# checking sorting of CpGs
identical(rownames(validation_data_mval),rownames(cvfit_coef))
# checking sorting of Sample_Name in both mval and cov files
identical(colnames(validation_data_mval),validation_data_cov$Sample_Name)


############################################
############################################
# EVALUATING PERFORMANCE IN VALIDATION SET #
############################################
############################################
num_significant_cpg_validation <- nrow(cvfit_coef)

# computing epigenetic scores
estimated_log2_UrAsgmCr_validation <- c()
if (num_significant_cpg_validation > 1) {
  for (current_sample in validation_data_cov$Sample_Name) {
    estimated_log2_UrAsgmCr_validation <- c(estimated_log2_UrAsgmCr_validation,
                                            (sum(as.numeric(validation_data_mval[,current_sample]) * cvfit_coef$logFC))
    )
  }
} else if (num_significant_cpg_validation == 1) {
  for (current_sample in validation_data_cov$Sample_Name) {
    estimated_log2_UrAsgmCr_validation <- c(estimated_log2_UrAsgmCr_validation,
                                            (sum(as.numeric(validation_data_mval[current_sample]) * cvfit_coef$logFC))
    )
  }
} else {
  print("ERROR!!!!")
}

# assessing correlations between our epigenetic marker and observed urinary arsenic levels
validation_data_cov$SCORE <- estimated_log2_UrAsgmCr_validation
output_cor <- cor(validation_data_cov$log2_UrAsgmCr,validation_data_cov$SCORE)
print(paste("Correlation (R2) Between Predicted and Measured Arsenic in Validation Set:", output_cor^2))
save(output_cor, file=paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/EPIGENETIC_MARKER_ARSENIC_EXPOSURE/EPIGENETIC_MARKER_VALIDATION_R2.RData"))

# outputting a summary of the regression with covariates
print(summary(lm(data=validation_data_cov, formula = log2_UrAsgmCr ~ SCORE + Age + factor(Sex) + factor(bmi_cat) + factor(cigsmoke) + factor(Batch_Plate))))

# plotting scatter plot of epigenetically predicted and measured urinary arsenic levels in the HEALS validation subcohort
HEALS_predicted_measured_As <- data.frame(`Measured Urinary Arsenic Levels` = validation_data_cov$log2_UrAsgmCr, `Epigenetically Predicted Urinary Arsenic Levels` = estimated_log2_UrAsgmCr_validation)
colnames(HEALS_predicted_measured_As) <- c("Measured Urinary Arsenic Levels", "Epigenetically Predicted Urinary Arsenic Levels")
summary(lm(HEALS_predicted_measured_As, formula = `Measured Urinary Arsenic Levels` ~ `Epigenetically Predicted Urinary Arsenic Levels`))
png("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/EPIGENETIC_MARKER_ARSENIC_EXPOSURE/VALIDATION_CORRELATION_PLOT_HEALS.png")
p <- ggplot(data=HEALS_predicted_measured_As, aes(x=`Epigenetically Predicted Urinary Arsenic Levels`,y=`Measured Urinary Arsenic Levels`)) + geom_point() + geom_smooth(method = "lm") + xlab("Epigenetically Predicted Urinary Arsenic Score") + theme_minimal()
print(p)
dev.off()

# outputting the beta coefficients for the epigenetic marker
cvfit_coef$Name <- rownames(cvfit_coef)
EPIGENETIC_MARKER_DF <- head(cvfit_coef,num_significant_cpg_validation) %>% select(Name, logFC)
save(EPIGENETIC_MARKER_DF, file=paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/EPIGENETIC_MARKER_ARSENIC_EXPOSURE/EPIGENETIC_MARKER_DF.RData"))

###################################################################################################
# ASSESSING THE PREDICTIVE POWER OF EPIGENETIC MARKER ON SKIN LESION STATUS IN THE VALIDATION SET #
###################################################################################################
# evaluating the predictive power of our epigenetic marker on skin lesion status
library(pROC)
roc_score=roc(as.factor(validation_data_cov$anyskinlesion), validation_data_cov$SCORE) #AUC score
output_auc <- auc(roc_score)
print(paste("AUC of Baseline Skin Lesion Status:", toString(sort(as.numeric(ci(output_auc))))))
save(output_auc, file=paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/EPIGENETIC_MARKER_ARSENIC_EXPOSURE/EPIGENETIC_MARKER_VALIDATION_SL_AUC.RData"))

# covariate-adjusted AUC
library(ROCnReg)
validation_data_cov_BMI <- validation_data_cov %>% filter(bmi_cat %in% c("Underweight","Normal","Overweight"))
output_AROC.sp <- AROC.sp(
  formula.h = SCORE ~ Age + factor(Sex) + factor(bmi_cat) + factor(cigsmoke),
  group = "anyskinlesion",
  data = validation_data_cov_BMI,
  tag.h = 0)
print(paste("AUC of Baseline Skin Lesion Status:", toString(sort(as.numeric(output_AROC.sp$AUC)))))

# plotting the ROC curve of SL status based on our epigenetic biomarker
library(pROC)
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/EPIGENETIC_MARKER_ARSENIC_EXPOSURE")
png(paste0("EPIGENETIC_MARKER_ROC_SL.png"))
plot(roc_score, main ="ROC curve of Skin Lesion Status and Epigenetic Marker \nof Arsenic Exposure")
dev.off()

# print out seed and string
output_string = paste(sbatch_index,nrow(EPIGENETIC_MARKER_DF),nrow(significant_results_final),output_cor^2,as.numeric(output_AROC.sp$AUC[2]),as.numeric(output_AROC.sp$AUC[1]),as.numeric(output_AROC.sp$AUC[3]),sep=",")
output_string_command = paste0("echo ",output_string," >> /gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/EPIGENETIC_MARKER_ARSENIC_EXPOSURE/output.log")
system(output_string_command)

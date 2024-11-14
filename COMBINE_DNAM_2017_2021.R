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

####################################
# importing 2021 data and covariates
# setting working directory and loading in normalized and preprocessed beta matrix
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/champ_idat")
file_name<-paste0("IMPUTED_BMIQ_2021.RData")
load(file_name)
# normalizing beta values into M-values
beta <- data.matrix(get(paste0("IMPUTED_BMIQ_2021")))
mval4 <- log2(beta/(1 - beta)) 
mval4 <- mval4[, order(colnames(mval4))]
mval_2021 <- mval4
# comment out below if there are already X before each sample name
colnames(mval_2021) <- paste0("X",colnames(mval_2021))

#reading in outcome and covariates data
cov <- read.table("sample_annotation.csv",header = TRUE,sep=",")
cov$Sample_Name2<-paste0("X",cov$Sample_Name)
cov$Sample_Name<-cov$Sample_Name2
cov$Sample_Name2<-NULL
cov$log2_UrAsgmCr<-log2(cov$UrAsgmCr)
cov$log2_WArsenic<-log2(cov$WArsenic)
names.use <- colnames(mval_2021)[(colnames(mval_2021) %in% cov$Sample_Name)]
mval_2021 <- mval_2021[, names.use]
cov<-cov[cov$Sample_Name %in% colnames(mval_2021),]
cov_2021<-cov
cov_2021<-cov_2021%>%select(Sample_Name,DMA_pct,Age,Sex,bmi_cat,cigsmoke,log2_UrAsgmCr,log2_WArsenic,anyskinlesion,Plate)
# formatting sample_plate to be batch specific
cov_2021$Sample_Plate<-paste0("2021_",cov_2021$Plate); cov_2021$Plate<-NULL;

####################################
# importing 2017 data and covariates
# setting working directory and loading in normalized and preprocessed beta matrix
load("/gpfs/data/phs/groups/Projects/GEMS/Kathryn/heals_beta_normalized_2-19-2019.RData")
dim(heals_beta_normalized)
heals_beta_normalized <- heals_beta_normalized[,order(colnames(heals_beta_normalized))]
heals_beta_normalized <- heals_beta_normalized[order(rownames(heals_beta_normalized)), ]
beta<-heals_beta_normalized
mval4 <- log2(beta/(1 - beta)) 
mval4 <- mval4[, order(colnames(mval4))]
mval_2017<-mval4

#reading in exposure and covariates data
cov <- read.table("/gpfs/data/phs/groups/Projects/GEMS/Kathryn/covariates/covariates 4-12-2017.csv",header = TRUE,sep=",")
cov$Sample_Name2<-paste0("X",cov$Sample_Name)
cov$Sample_Name<-cov$Sample_Name2
cov$Sample_Name2<-NULL
cov$log2_UrAsgmCr<-log2(cov$UrAsgmCr)
cov$log2_WArsenic<-log2(cov$WArsenic)
names.use <- colnames(mval_2017)[(colnames(mval_2017) %in% cov$Sample_Name)]
mval_2017 <- mval_2017[, names.use]
cov<-cov[cov$Sample_Name %in% colnames(mval_2017),]
cov_2017<-cov
# importing skin lesion data and joining skin lesion data with 2017 pheno data
dt_lesion<-fread("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/pheno_data/LinTong/HEALS_SkinLesion.txt",header=T)
names(dt_lesion)[1]<-"SubjectID"
tmp_cov_2017_sl<-inner_join(cov_2017,dt_lesion,by=c("SubjectID"))
# selecting columns of interest from pheno data
cov_2017<-tmp_cov_2017_sl%>%select(Sample_Name,DMA_pct,Age,Sex,BMICat,CigNFC,log2_UrAsgmCr,log2_WArsenic,anyskinlesion,Sample_Plate)
# reformatting BMI category labels
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
for (i in 1:length(cov_2017$BMICat)) {
    cov_2017$BMICat[i]<-firstup(cov_2017$BMICat[i])
}
cov_2017$CigNFC<-cov_2017$CigNFC-1
colnames(cov_2017)[5]<-"bmi_cat"
colnames(cov_2017)[6]<-"cigsmoke"
# formatting sample_plate to be batch specific
cov_2017$Sample_Plate<-paste0("2017_",cov_2017$Sample_Plate)

# combining covariate tables from both years
combined_cov<-rbind(cov_2021,cov_2017)
col_rename_index<-which(grepl("Sample_Plate",colnames(combined_cov)))
colnames(combined_cov)[col_rename_index]<-"Batch_Plate"
combined_cov$Batch_Plate<-as.factor(combined_cov$Batch_Plate)
combined_cov$Sex<-combined_cov$Sex-1

# combining mval matrices from both years
shared_CpG_list_2017_2021<-intersect(rownames(mval_2021),rownames(mval_2017))
mval_2021<-mval_2021[(rownames(mval_2021) %in% shared_CpG_list_2017_2021),]
mval_2017<-mval_2017[(rownames(mval_2017) %in% shared_CpG_list_2017_2021),]
mval_2021<-mval_2021[ order((row.names(mval_2021))), ]
mval_2017<-mval_2017[ order((row.names(mval_2017))), ]
# testing if all CpGs match between the two batches
identical(row.names(mval_2021),row.names(mval_2017))
combined_mval<-cbind(mval_2021,mval_2017)
combined_mval<-data.matrix(combined_mval)

# converting combined_mval matrix into combined_beta matrix
tmp_mval2beta<-2^(combined_mval)
combined_beta<-tmp_mval2beta / (1 + tmp_mval2beta)
combined_beta<-data.matrix(combined_beta)

# saving combined_mval, combined_beta, and combined_cov files
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data")
save(combined_mval, file = "/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/combined_mval.RData")
save(combined_beta, file = "/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/combined_beta.RData")
save(combined_cov, file = "/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/combined_cov.RData")

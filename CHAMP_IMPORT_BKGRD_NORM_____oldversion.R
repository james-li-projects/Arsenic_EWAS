library("ChAMP")
library("dplyr")
library(data.table)
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/champ_idat")
########################################
################ IMPORT ################
########################################
# loading using minfi import in ChAMP with imputation
myLoad <- champ.load(method="minfi",arraytype="EPIC",population="SAS",force=TRUE,detPcut=0.01,SampleCutoff=0.1,ProbeCutoff=0.1,autoimpute=TRUE)

########################################
################ ssNOOB ################
########################################
# performing ssNoob background correction
noob<-preprocessNoob(myLoad$rgSet, offset = 15, dyeCorr = TRUE, verbose = FALSE, dyeMethod=c("single"))
noob_beta<-data.frame(getBeta(noob))
noob_beta$cpg_id<-rownames(noob_beta)

# retrieving CpGs that passed ChAMP's initial filtering step
filtered_champ_beta<-data.frame(myLoad$beta)
filtered_champ_beta$cpg_id<-rownames(myLoad$beta)
filtered_champ_CpG_list<-filtered_champ_beta%>%select(cpg_id)

# retrieving ssNoob normalized beta values only for CpGs that passed ChAMP's initial filtering step
library(dplyr)
minfi_beta<-inner_join(filtered_champ_CpG_list,noob_beta,by=c("cpg_id"))
rownames(minfi_beta)<-minfi_beta$cpg_id
minfi_beta$cpg_id<-NULL
detectCores()
save(minfi_beta,file="minfi_beta.RData")

########################################
################# BMIQ #################
########################################
# performing BMIQ normalization
library("ChAMP")
load("minfi_beta.RData")
noob_bmiq<-champ.norm(beta=minfi_beta, method="BMIQ", arraytype="EPIC",cores=16)
save(noob_bmiq,file="noob_bmiq.RData")
load("noob_bmiq.RData")

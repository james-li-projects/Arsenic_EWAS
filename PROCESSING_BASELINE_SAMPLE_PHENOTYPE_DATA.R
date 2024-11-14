library(data.table)
library(dplyr)

setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/pheno_data/LinTong")
# reading in skin lesion and other phenotype data
dt_lesion<-fread("HEALS_SkinLesion.txt",header=T)
names(dt_lesion)[1]<-"SubjectID"
dt_pheno<-fread("HEALS_pheno.txt",header=T)
dt_MERGED_pheno<-inner_join(dt_lesion,dt_pheno,by=c("SubjectID"))

# checking if the missing values in UrineAs and UrineCreat are the same
hi<-which(is.na(dt_MERGED_pheno$UrineAs))
hi2<-which(is.na(dt_MERGED_pheno$UrineCreat))
identical(hi,hi2)

# computing UrAsgmCr (creatinine adjusted urinary arsenic) for values that have NA values
dt_MERGED_pheno$UrAsgmCr_computed=(dt_MERGED_pheno$UrineAs*100)/dt_MERGED_pheno$UrineCreat
dt_MERGED_pheno$UrAsgmCr[is.na(dt_MERGED_pheno$UrAsgmCr)]<-dt_MERGED_pheno$UrAsgmCr_computed[is.na(dt_MERGED_pheno$UrAsgmCr)]

# extracting data only for 2021 methylation sample
methyl2021_id_list<-fread("/gpfs/data/phs/groups/Projects/GEMS/EPIC_array/EPIC_data_2021/sampleID_methylation.txt",header=T)
names(methyl2021_id_list)<-"SubjectID"
tmp_dt_2021<-inner_join(methyl2021_id_list,dt_MERGED_pheno)
names(tmp_dt_2021)[1]<-"Sample_ID"
tmp_dt_2021$Sample_ID<-as.character(tmp_dt_2021$Sample_ID)
# creating bmi category variable
tmp_dt_2021$bmi_cat<-"Unknown"
tmp_dt_2021$bmi_cat[tmp_dt_2021$bmi<18.5]<-"Underweight"
tmp_dt_2021$bmi_cat[tmp_dt_2021$bmi<25 & tmp_dt_2021$bmi>=18.5]<-"Normal"
tmp_dt_2021$bmi_cat[tmp_dt_2021$bmi<30 & tmp_dt_2021$bmi>=25]<-"Overweight"
tmp_dt_2021$bmi_cat[tmp_dt_2021$bmi>=30]<-"Obese"

# selecting covariates of interest
tmp_dt_2021_select<-tmp_dt_2021%>%select(Sample_ID,anyskinlesion,Sex,Age,UrAsgmCr,InAS_pct,MMA_pct,DMA_pct,WArsenic,bmi_cat,cigsmoke)
# recoding missing skin lesion data to "missing"
tmp_dt_2021_select$anyskinlesion[is.na(tmp_dt_2021_select$anyskinlesion)]<-"missing"
# just in case, drop samples that are missing any covariates, but check to make sure there are no dropped individuals, unless there is a good reason to drop
tmp_dt_2021<-na.omit(tmp_dt_2021_select)

# reading in sample well and plate annotations and joining them to our previous covariate table
sample_well_plate<-fread("../Infinium_MethylationEPIC_HEALS_Batch_01thru08_BP20210505_Sentrix.csv",header=T)
names(sample_well_plate)[1]<-"Sample_ID"
# dropping the Pool_ID column
sample_well_plate$Pool_ID<-NULL
dt_2021<-inner_join(tmp_dt_2021,sample_well_plate,by = c("Sample_ID"))
# renaming Sample_ID column to Sample_Name
names(dt_2021)[1]<-"Sample_Name"

## identifying duplicated samples and renaming one of the duplicates to have the suffix _DUP
#dup_index<-which(duplicated(dt_2021$Sample_Name))
#dups<-(dt_2021[c(dup_index),])$Sample_Name
## below unique_dups is redundant to the dups vector
#unique_dups<-unique(dt_2021[dt_2021$Sample_Name %in% dups,]$Sample_Name)
#all_dup_index<-which(dt_2021$Sample_Name %in% unique_dups)
#for (i in 1:length(unique_dups)) {
#    dup_indices<-unique_dups[i]
#    unchanged_index_dup<-(head(which(dt_2021$Sample_Name %in% dup_indices),1))
#    change_index_dup<-(tail(which(dt_2021$Sample_Name %in% dup_indices),1))
#    dt_2021$Sample_Name[change_index_dup]<-paste0(dt_2021$Sample_Name[change_index_dup],"_DUP")
#}
##verifying the Sample_Name for duplicates were changed
#dt_2021$Sample_Name[all_dup_index]

# identifying duplicated samples and removing one of the duplicates in each duplicated pair
dup_index<-which(duplicated(dt_2021$Sample_Name))
dups<-(dt_2021[c(dup_index),])$Sample_Name
# below unique_dups is redundant to the dups vector
unique_dups<-unique(dt_2021[dt_2021$Sample_Name %in% dups,]$Sample_Name)
all_dup_index<-which(dt_2021$Sample_Name %in% unique_dups)
for (i in 1:length(unique_dups)) {
    dup_indices<-unique_dups[i]
    unchanged_index_dup<-(head(which(dt_2021$Sample_Name %in% dup_indices),1))
    change_index_dup<-(tail(which(dt_2021$Sample_Name %in% dup_indices),1))
    dt_2021$Sample_Name[change_index_dup]<-paste0(dt_2021$Sample_Name[change_index_dup],"_DUP")
}
# checking duplicates were renamed
dt_2021$Sample_Name[all_dup_index]
# removing duplicates
nrow(dt_2021)
dt_2021<-dt_2021[!grepl("_DUP", dt_2021$Sample_Name),]
nrow(dt_2021)

#writing out sample annotation file in multiple locations
write.csv(dt_2021,file="/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/idat/sample_annotation.csv",quote=F,row.names=F)
write.csv(dt_2021,file="/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/champ_idat/sample_annotation.csv",quote=F,row.names=F)
write.csv(dt_2021,file="/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/sample_annotation.csv",quote=F,row.names=F)

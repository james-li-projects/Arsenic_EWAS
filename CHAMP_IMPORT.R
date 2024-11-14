library(ChAMP)
library(dplyr)
library(data.table)

setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/champ_idat")

# loading using minfi import in ChAMP
myLoad <- champ.load(method="minfi",arraytype="EPIC",population="SAS",force=TRUE,detPcut=0.01,SampleCutoff=0.1,ProbeCutoff=0.1,autoimpute=FALSE)

# performing imputation using the KNN method (k=10)
myImpute <- champ.impute(method="KNN",k=10,SampleCutoff=0.1,ProbeCutoff=0.1)
IMPUTED_ONLY_2021 <- myImpute$beta
save(IMPUTED_ONLY_2021,file="IMPUTED_ONLY_2021.RData")

print("Only Imputation")
print(myImpute$beta[894,1])
print(myImpute$beta[1616,1])
print(myImpute$beta[1812,1])
print(myImpute$beta[2626,1])
print(myImpute$beta[3762,1])
print(myImpute$beta[13550,1])

# performing BMIQ normalization
IMPUTED_BMIQ_2021<-champ.norm(beta=myImpute$beta, method="BMIQ", arraytype="EPIC",cores=16)
save(IMPUTED_BMIQ_2021,file="IMPUTED_BMIQ_2021.RData")

print("Imputation + BMIQ")
print(IMPUTED_BMIQ_2021[894,1])
print(IMPUTED_BMIQ_2021[1616,1])
print(IMPUTED_BMIQ_2021[1812,1])
print(IMPUTED_BMIQ_2021[2626,1])
print(IMPUTED_BMIQ_2021[3762,1])
print(IMPUTED_BMIQ_2021[13550,1])

# examining effect of BMIQ on imputed CpGs
NA_rowcol <- data.frame((which(is.na(myLoad$beta),arr.ind=TRUE)))
save(NA_rowcol, file = "NA_rowcol.RData")

vec_NA_IMPUTED_ONLY <- vector(length = nrow(NA_rowcol))
vec_NA_IMPUTED_BMIQ <- vector(length = nrow(NA_rowcol))

for (i in 1:nrow(NA_rowcol)) {
    index_row <- NA_rowcol$row[i]
    index_col <- NA_rowcol$col[i]
    vec_NA_IMPUTED_ONLY[i] <- IMPUTED_ONLY_2021[index_row,index_col]
    vec_NA_IMPUTED_BMIQ[i] <- IMPUTED_BMIQ_2021[index_row,index_col]
}
identical(vec_NA_IMPUTED_ONLY,vec_NA_IMPUTED_BMIQ)

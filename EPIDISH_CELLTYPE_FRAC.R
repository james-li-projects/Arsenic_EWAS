library(EpiDISH)
library(dplyr)
library(data.table)

data(centDHSbloodDMC.m)

setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data")
load("combined_cov.RData")
load("combined_beta.RData")
load("combined_mval.RData")

######################################
##### Saving cell type fractions #####
######################################
out.l <- epidish(beta.m = combined_beta, ref.m = centDHSbloodDMC.m, method = 'RPC')
frac.m <- out.l$estF
CT_frac <- data.frame(frac.m)
CT_frac$Sample_Name <- rownames(CT_frac)
save(CT_frac,file="/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/CT_frac.RData")

##################################
# updating the combined_cov file #
##################################
# loading the above cell type fraction output file
load("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/CT_frac.RData")
# removing granulocytes to avoid collinearity issues
CT_frac <- CT_frac %>% select(-Neutro,-Eosino)
# normalizing cell type proportions to be z-scores
stand_z_func <- function(x) {
    return( (x-mean(x)) / sd(x) )
}
CT_frac$B <- stand_z_func(CT_frac$B)
CT_frac$NK <- stand_z_func(CT_frac$NK)
CT_frac$CD4T <- stand_z_func(CT_frac$CD4T)
CT_frac$CD8T <- stand_z_func(CT_frac$CD8T)
CT_frac$Mono <- stand_z_func(CT_frac$Mono)

temp_combined_cov <- inner_join(combined_cov,CT_frac,by=c("Sample_Name"))
combined_cov <- temp_combined_cov
save(combined_cov,file="/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/combined_cov.RData")

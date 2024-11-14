library(data.table)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(survival)
library(ranger)
library(survival)

# import SHS general database
load("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/download_data/SHS/IMPORTED_COV.RData")
colnames(IMPORTED_COV)[3] <- "idno"
IMPORTED_COV <- IMPORTED_COV %>% select(idno,as.cr.ln,estimated_log2_UrAsgmCr,s1ldlest,s1ldlbq,sex,center,s1smoke,s1etoh,s1bmi,s1htnrx2,s1dmwhb,s1sbp,s1acr,s1ldlbq,s1hdl,s1age,s1ckdepi.new,s1hba1c,s1whr,s1g0,s1edu,ab.kations,sum.ias.cr,as.cr,s1ckdepi.new,s1htnhx)
# changing acr variable to a factor
IMPORTED_COV$s1acr <- as.factor(IMPORTED_COV$s1acr)
# recoding the education variable
IMPORTED_COV$s1edu_cat[IMPORTED_COV$s1edu<=8] <- "< high school"
IMPORTED_COV$s1edu_cat[IMPORTED_COV$s1edu<12 & IMPORTED_COV$s1edu>=9] <- "Some high school"
IMPORTED_COV$s1edu_cat[IMPORTED_COV$s1edu==12] <- "High school diploma"
IMPORTED_COV$s1edu_cat[IMPORTED_COV$s1edu>12] <- "More than high school"

# SHS CVD follow-up database
shs_cvd <- fread("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/download_data/SHS/shs_ccvd2017.csv",sep=",",header=T)
names(shs_cvd) <- casefold(names(shs_cvd))
shs_cvd$idno <- as.character(shs_cvd$idno)
shs_cvd <- merge(shs_cvd, IMPORTED_COV, by='idno', all=FALSE)

# Closing date
shs_cvd$closingdate <- as.Date(as.character(shs_cvd$closingdate), format="%m/%d/%Y")

# Visit 1 date
shs_cvd$s1examdt <- as.Date(as.character(shs_cvd$s1examdt), format="%m/%d/%Y")
# Date of death
shs_cvd$dod <- as.Date(as.character(shs_cvd$dod), format="%m/%d/%Y")

# comeback - set closing date to 2009
shs_cvd$closingdate[shs_cvd$closingdate=="2017-12-31"] <- "2009-12-31"
shs_cvd$dod[shs_cvd$dod>as.Date("2009-12-31")] <- NA

# Date of CVD (and subtypes) incidence
shs_cvd$anycvddts1 <- as.Date(as.character(shs_cvd$anycvddts1), format="%m/%d/%Y")
shs_cvd$anydstkdts1 <- as.Date(as.character(shs_cvd$anydstkdts1), format="%m/%d/%Y")
shs_cvd$anychddts1 <- as.Date(as.character(shs_cvd$anychddts1), format="%m/%d/%Y")
shs_cvd$anychfdts1 <- as.Date(as.character(shs_cvd$anychfdts1), format="%m/%d/%Y")

# DEFINE SURVIVAL ANALYSIS VARIABLES - ENTRY, EXIT, TIME METRIC, FOLLOW-UP TIME

# COMPOSITE CVD INCIDENCE 
shs_cvd$peryr.exam.anycvd <- NULL
shs_cvd$peryr.exam.anycvd <- ifelse(shs_cvd$closingdate > shs_cvd$anycvddts1,(shs_cvd$anycvddts1 - shs_cvd$s1examdt)/365.25,
                                    (shs_cvd$closingdate - shs_cvd$s1examdt)/365.25)
shs_cvd$peryr.exam.anycvd[is.na(shs_cvd$anycvddts1)&!is.na(shs_cvd$dod)] <- (shs_cvd$dod[is.na(shs_cvd$anycvddts1)&!is.na(shs_cvd$dod)] -
                                                                               shs_cvd$s1examdt[is.na(shs_cvd$anycvddts1)&!is.na(shs_cvd$dod)])/365.25
shs_cvd$peryr.exam.anycvd[is.na(shs_cvd$anycvddts1)&is.na(shs_cvd$peryr.exam.anycvd)] <-
  (shs_cvd$closingdate[is.na(shs_cvd$anycvddts1)&is.na(shs_cvd$peryr.exam.anycvd)] - shs_cvd$s1examdt[is.na(shs_cvd$anycvddts1)&is.na(shs_cvd$peryr.exam.anycvd)])/365.25
shs_cvd$peryr.age.anycvd <- shs_cvd$peryr.exam.anycvd + shs_cvd$s1age
shs_cvd$anycvd <-  ifelse(shs_cvd$anycvds1==1&!is.na(shs_cvd$anycvds1),1,0) 

# STROKE INCIDENCE
shs_cvd$peryr.exam.anystk <- NULL
shs_cvd$peryr.exam.anystk <- ifelse(shs_cvd$closingdate > shs_cvd$anydstkdts1, (shs_cvd$anydstkdts1 - shs_cvd$s1examdt)/365.25,
                                    (shs_cvd$closingdate - shs_cvd$s1examdt)/365.25)
shs_cvd$peryr.exam.anystk[is.na(shs_cvd$anydstkdts1)&!is.na(shs_cvd$dod)] <- (shs_cvd$dod[is.na(shs_cvd$anydstkdts1)&!is.na(shs_cvd$dod)] -
                                                                                shs_cvd$s1examdt[is.na(shs_cvd$anydstkdts1)&!is.na(shs_cvd$dod)])/365.25
shs_cvd$peryr.exam.anystk[is.na(shs_cvd$anydstkdts1)&is.na(shs_cvd$peryr.exam.anystk)] <-
  (shs_cvd$closingdate[is.na(shs_cvd$anydstkdts1)&is.na(shs_cvd$peryr.exam.anystk)] - shs_cvd$s1examdt[is.na(shs_cvd$anydstkdts1)&is.na(shs_cvd$peryr.exam.anystk)])/365.25
shs_cvd$peryr.age.anystk <- shs_cvd$peryr.exam.anystk + shs_cvd$s1age
shs_cvd$anystk <-  ifelse(shs_cvd$anydstks1==1&!is.na(shs_cvd$anydstks1),1,0)

# CORONARY HEART DISEASE INCIDENCE
shs_cvd$peryr.exam.anychd <- NULL
shs_cvd$peryr.exam.anychd <- ifelse(shs_cvd$closingdate > shs_cvd$anychddts1, (shs_cvd$anychddts1 - shs_cvd$s1examdt)/365.25,
                                    (shs_cvd$closingdate - shs_cvd$s1examdt)/365.25)
shs_cvd$peryr.exam.anychd[is.na(shs_cvd$anychddts1)&!is.na(shs_cvd$dod)] <- (shs_cvd$dod[is.na(shs_cvd$anychddts1)&!is.na(shs_cvd$dod)] -
                                                                               shs_cvd$s1examdt[is.na(shs_cvd$anychddts1)&!is.na(shs_cvd$dod)])/365.25
shs_cvd$peryr.exam.anychd[is.na(shs_cvd$anychddts1)&is.na(shs_cvd$peryr.exam.anychd)] <-
  (shs_cvd$closingdate[is.na(shs_cvd$anychddts1)&is.na(shs_cvd$peryr.exam.anychd)] - shs_cvd$s1examdt[is.na(shs_cvd$anychddts1)&is.na(shs_cvd$peryr.exam.anychd)])/365.25
shs_cvd$peryr.age.anychd <- shs_cvd$peryr.exam.anychd + shs_cvd$s1age
shs_cvd$anychd <-  ifelse(shs_cvd$anychds1==1&!is.na(shs_cvd$anychds1),1,0)

# CONGESTIVE HEART FAILURE INCIDENCE
shs_cvd$peryr.exam.anychf <- NULL
shs_cvd$peryr.exam.anychf <- ifelse(shs_cvd$closingdate > shs_cvd$anychfdts1, (shs_cvd$anychfdts1 - shs_cvd$s1examdt)/365.25,
                                    (shs_cvd$closingdate - shs_cvd$s1examdt)/365.25)
shs_cvd$peryr.exam.anychf[is.na(shs_cvd$anychfdts1)&!is.na(shs_cvd$dod)] <- (shs_cvd$dod[is.na(shs_cvd$anychfdts1)&!is.na(shs_cvd$dod)] -
                                                                               shs_cvd$s1examdt[is.na(shs_cvd$anychfdts1)&!is.na(shs_cvd$dod)])/365.25
shs_cvd$peryr.exam.anychf[is.na(shs_cvd$anychfdts1)&is.na(shs_cvd$peryr.exam.anychf)] <-
  (shs_cvd$closingdate[is.na(shs_cvd$anychfdts1)&is.na(shs_cvd$peryr.exam.anychf)] - shs_cvd$s1examdt[is.na(shs_cvd$anychfdts1)&is.na(shs_cvd$peryr.exam.anychf)])/365.25
shs_cvd$peryr.age.anychf <- shs_cvd$peryr.exam.anychf + shs_cvd$s1age
shs_cvd$anychf <-  ifelse(shs_cvd$anychfs1==1&!is.na(shs_cvd$anychfs1),1,0)

# MORTALITY F/U TIME
shs_cvd$peryr.exam.anycvd.death <- NULL
shs_cvd$peryr.exam.anycvd.death <- ifelse((shs_cvd$closingdate > shs_cvd$dod), (shs_cvd$dod - shs_cvd$s1examdt)/365.25,
                                          (shs_cvd$closingdate - shs_cvd$s1examdt)/365.25)
shs_cvd$peryr.exam.anycvd.death[is.na(shs_cvd$dod)] <- (shs_cvd$closingdate[is.na(shs_cvd$dod)] - shs_cvd$s1examdt[is.na(shs_cvd$dod)])/365.25
shs_cvd$peryr.age.anycvd.death <- shs_cvd$peryr.exam.anycvd.death + shs_cvd$s1age

shs_cvd$anycvd.death <-  ifelse(shs_cvd$fatalcvd==1&!is.na(shs_cvd$fatalcvd),1,0)
shs_cvd$anychd.death <-  ifelse(shs_cvd$fatalchd==1&!is.na(shs_cvd$fatalchd),1,0)
shs_cvd$anychf.death <-  ifelse(shs_cvd$fatalchf==1&!is.na(shs_cvd$fatalchf),1,0)
shs_cvd$anystk.death <-  ifelse(shs_cvd$fatalstk==1&!is.na(shs_cvd$fatalstk),1,0)

# New LDL cholesterol variable and sex
shs_cvd$s1ldl[!is.na(shs_cvd$s1ldlest)] <- shs_cvd$s1ldlest[!is.na(shs_cvd$s1ldlest)]
shs_cvd$s1ldl[is.na(shs_cvd$s1ldlest)]  <- shs_cvd$s1ldlbq[is.na(shs_cvd$s1ldlest)]
summary(shs_cvd$s1ldl,na.rm=TRUE, digits=3)
shs_cvd$sex <- ifelse(shs_cvd$sex=="F",1,0)

# coding follow-up time for people who passed away
shs_cvd$day_FU[!is.na(shs_cvd$dod)] <- shs_cvd$dod[!is.na(shs_cvd$dod)] -
  shs_cvd$s1examdt[!is.na(shs_cvd$dod)]

# coding follow-up time for people who did not pass away
shs_cvd$day_FU[is.na(shs_cvd$dod)] <- shs_cvd$closingdate[is.na(shs_cvd$dod)] -
  shs_cvd$s1examdt[is.na(shs_cvd$dod)] 

# obtaining year_FU and dead variables
shs_cvd <- shs_cvd %>% filter(!is.na(day_FU))
shs_cvd$year_FU <- shs_cvd$day_FU/365.25 
shs_cvd$dead[!is.na(shs_cvd$dod)] <- 1
shs_cvd$dead[is.na(shs_cvd$dod)] <- 0
shs_cvd$dead <- as.numeric(shs_cvd$dead)
shs_cvd$year_FU <- shs_cvd$year_FU + shs_cvd$s1age 

# Exclude baseline CVD and missings in relevant covariates
no.cvd <- shs_cvd[shs_cvd$s1cvdfree==1,]	
shs.cox <- no.cvd[             					 
  !is.na(no.cvd$s1age)&
    !is.na(no.cvd$sex)&
    !is.na(no.cvd$center)&
    !is.na(no.cvd$s1smoke)&                          
    !is.na(no.cvd$s1bmi)&
    !is.na(no.cvd$s1ldl)&
    !is.na(no.cvd$s1htnrx2)&
    !is.na(no.cvd$s1dmwhb)&
    !is.na(no.cvd$s1sbp)&
    !is.na(no.cvd$s1acr)&
    !is.na(no.cvd$s1hdl)&
    !is.na(no.cvd$as.cr.ln)&
    !is.na(no.cvd$estimated_log2_UrAsgmCr),]

dim(shs.cox)

###########################################
# correlation between measured (sum.ias.cr) and predicted arsenic
summary(lm(data = shs.cox, formula = estimated_log2_UrAsgmCr ~ log2(sum.ias.cr) + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1whr + s1ckdepi.new))

# SLR correlation
summary(lm(data = shs.cox, formula = log2(sum.ias.cr) ~ estimated_log2_UrAsgmCr))

# plotting correlation
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/SHS/")
png("SHS_CORRELATION.png",units="in",width=5,height=5,res=1200)
ggplot(data = shs.cox, aes(x=estimated_log2_UrAsgmCr,y=log2(sum.ias.cr))) + geom_point() + geom_smooth(formula = y ~ x) + theme_classic() + xlab("Epigenetically Predicted Urinary Arsenic Scores") + ylab("Measured Urinary Arsenic Levels") + annotate(geom = "text", x = -0.70, y=1.6, label=expression(paste(R^2,": 0.01")),size=5) + annotate(geom = "text", x = -0.70, y=1.1, label=expression(paste("p-value: 6.03E-06")),size=5) + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) 
dev.off()

###########################################

shs.cox$year_FU <- shs.cox$year_FU - shs.cox$s1age
shs.cox$peryr.age.anycvd.death <- shs.cox$peryr.age.anycvd.death - shs.cox$s1age
shs.cox$peryr.age.anycvd <- shs.cox$peryr.age.anycvd - shs.cox$s1age
shs.cox$peryr.age.anychf <- shs.cox$peryr.age.anychf - shs.cox$s1age
shs.cox$peryr.age.anystk <- shs.cox$peryr.age.anystk - shs.cox$s1age

###################################
############# MODEL 1 #############
###################################

# initializing vectors for storing association estimates
OR <- c()
Lower <- c()
Upper <- c()
Model <- c()
Exposure <- c()
Pvalue <- c()

##################################################
# survival outcomes for sum.ias.cr

# overall mortality
fit<-coxph(formula = Surv(year_FU,dead) ~ log2(sum.ias.cr) + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new, data = shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #1")
Exposure <- c(Exposure, "Measured Urinary Arsenic Levels")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# CVD
fit<-coxph(Surv(peryr.age.anycvd, event=anycvd) ~ log2(sum.ias.cr) + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #1")
Exposure <- c(Exposure, "Measured Urinary Arsenic Levels")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# CVD death
fit<-coxph(Surv(peryr.age.anycvd.death, event=anycvd.death) ~ log2(sum.ias.cr) + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #1")
Exposure <- c(Exposure, "Measured Urinary Arsenic Levels")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# CHF
fit<-coxph(Surv(peryr.age.anychf, event=anychf) ~ log2(sum.ias.cr) + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #1")
Exposure <- c(Exposure, "Measured Urinary Arsenic Levels")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# STK
fit<-coxph(Surv(peryr.age.anystk, event=anystk) ~ log2(sum.ias.cr) + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #1")
Exposure <- c(Exposure, "Measured Urinary Arsenic Levels")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

##################################################
# our epigenetic marker
# overall mortality 
fit<-coxph(formula = Surv(year_FU,dead) ~ estimated_log2_UrAsgmCr + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new, data = shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #1")
Exposure <- c(Exposure, "Epigenetically Predicted Urinary Arsenic Scores")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# CVD
fit<-coxph(Surv(peryr.age.anycvd, event=anycvd) ~ estimated_log2_UrAsgmCr + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #1")
Exposure <- c(Exposure, "Epigenetically Predicted Urinary Arsenic Scores")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# CVD death
fit<-coxph(Surv(peryr.age.anycvd.death, event=anycvd.death) ~ estimated_log2_UrAsgmCr + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #1")
Exposure <- c(Exposure, "Epigenetically Predicted Urinary Arsenic Scores")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# CHF
fit<-coxph(Surv(peryr.age.anychf, event=anychf) ~ estimated_log2_UrAsgmCr + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #1")
Exposure <- c(Exposure, "Epigenetically Predicted Urinary Arsenic Scores")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# STK
fit<-coxph(Surv(peryr.age.anystk, event=anystk) ~ estimated_log2_UrAsgmCr + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #1")
Exposure <- c(Exposure, "Epigenetically Predicted Urinary Arsenic Scores")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

###################################
############# MODEL 2 #############
###################################

##################################################
# survival outcomes for sum.ias.cr

# overall mortality
fit<-coxph(formula = Surv(year_FU,dead) ~ log2(sum.ias.cr) + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new + s1ldl + s1hdl + s1sbp + s1dmwhb + s1htnhx + s1acr, data = shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #2")
Exposure <- c(Exposure, "Measured Urinary Arsenic Levels")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# CVD
fit<-coxph(Surv(peryr.age.anycvd, event=anycvd) ~ log2(sum.ias.cr) + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new + s1ldl + s1hdl + s1sbp + s1dmwhb + s1htnhx + s1acr, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #2")
Exposure <- c(Exposure, "Measured Urinary Arsenic Levels")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# CVD death
fit<-coxph(Surv(peryr.age.anycvd.death, event=anycvd.death) ~ log2(sum.ias.cr) + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new + s1ldl + s1hdl + s1sbp + s1dmwhb + s1htnhx + s1acr, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #2")
Exposure <- c(Exposure, "Measured Urinary Arsenic Levels")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# CHF
fit<-coxph(Surv(peryr.age.anychf, event=anychf) ~ log2(sum.ias.cr) + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new + s1ldl + s1hdl + s1sbp + s1dmwhb + s1htnhx + s1acr, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #2")
Exposure <- c(Exposure, "Measured Urinary Arsenic Levels")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# STK
fit<-coxph(Surv(peryr.age.anystk, event=anystk) ~ log2(sum.ias.cr) + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new + s1ldl + s1hdl + s1sbp + s1dmwhb + s1htnhx + s1acr, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #2")
Exposure <- c(Exposure, "Measured Urinary Arsenic Levels")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

##################################################
# our epigenetic marker
# overall mortality 
fit<-coxph(formula = Surv(year_FU,dead) ~ estimated_log2_UrAsgmCr + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new + s1ldl + s1hdl + s1sbp + s1dmwhb + s1htnhx + s1acr, data = shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #2")
Exposure <- c(Exposure, "Epigenetically Predicted Urinary Arsenic Scores")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# CVD
fit<-coxph(Surv(peryr.age.anycvd, event=anycvd) ~ estimated_log2_UrAsgmCr + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new + s1ldl + s1hdl + s1sbp + s1dmwhb + s1htnhx + s1acr, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #2")
Exposure <- c(Exposure, "Epigenetically Predicted Urinary Arsenic Scores")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# CVD death
fit<-coxph(Surv(peryr.age.anycvd.death, event=anycvd.death) ~ estimated_log2_UrAsgmCr + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new + s1ldl + s1hdl + s1sbp + s1dmwhb + s1htnhx + s1acr, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #2")
Exposure <- c(Exposure, "Epigenetically Predicted Urinary Arsenic Scores")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# CHF
fit<-coxph(Surv(peryr.age.anychf, event=anychf) ~ estimated_log2_UrAsgmCr + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new + s1ldl + s1hdl + s1sbp + s1dmwhb + s1htnhx + s1acr, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #2")
Exposure <- c(Exposure, "Epigenetically Predicted Urinary Arsenic Scores")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

# STK
fit<-coxph(Surv(peryr.age.anystk, event=anystk) ~ estimated_log2_UrAsgmCr + strata(center) + as.factor(sex) + as.factor(s1smoke) + s1bmi + s1ckdepi.new + s1ldl + s1hdl + s1sbp + s1dmwhb + s1htnhx + s1acr, data=shs.cox)
fit
exp(confint(fit))
OR <- c(OR, exp(fit$coefficients)[1])
Lower <- c(Lower, as.numeric(data.frame(exp(confint(fit)))[1,][1]))
Upper <- c(Upper, as.numeric(data.frame(exp(confint(fit)))[1,][2]))
Model <- c(Model,"Model #2")
Exposure <- c(Exposure, "Epigenetically Predicted Urinary Arsenic Scores")
Pvalue <- c(Pvalue,(summary(fit)$coefficients)[1,5])

######################################
# plotting KM curves
# KM for overall survival - measured As
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/SHS")
shs.cox_km <- shs.cox
shs.cox_km$quantile <- ntile(shs.cox_km$sum.ias.cr, 2)
shs.cox_km$`MeasuredUrinaryAs`[shs.cox_km$quantile==2] <- "Upper Half"
shs.cox_km$`MeasuredUrinaryAs`[shs.cox_km$quantile==1] <- "Lower Half"
shs.cox_km_fit <- survfit(Surv(year_FU,dead) ~ MeasuredUrinaryAs, data=shs.cox_km)
png("KM_OS_MEASURED_AS.png",res=1200,width=5,height=5,units="in")
autoplot(shs.cox_km_fit, conf.int = FALSE, censor.size = 1) + xlab("Follow-up Time (years)") + ylab("Survival") + theme_classic() + labs(colour = "MeasuredUrinaryAs") + theme(legend.position = c(.75,.75))
dev.off()

# KM for overall survival - predicted As 
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/SHS")
shs.cox_km <- shs.cox
shs.cox_km$quantile <- ntile(shs.cox_km$estimated_log2_UrAsgmCr, 2)
shs.cox_km$`EpigeneticAsScores`[shs.cox_km$quantile==2] <- "Upper Half"
shs.cox_km$`EpigeneticAsScores`[shs.cox_km$quantile==1] <- "Lower Half"
shs.cox_km_fit <- survfit(Surv(year_FU,dead) ~ EpigeneticAsScores, data=shs.cox_km)
png("KM_OS_EPIGENETIC_AS.png",res=1200,width=5,height=5,units="in")
autoplot(shs.cox_km_fit, conf.int = FALSE, censor.size = 1) + xlab("Follow-up Time (years)") + ylab("Survival") + theme_classic() + labs(colour = "EpigeneticAsScores") + theme(legend.position = c(.75,.75))
dev.off()

#################################
# Plot forest plots
#Outcome <- c(rep(c("05 Overall Mortality","03 CVD Incidence","04 CVD Mortality","02 Congestive Heart Failure","01 Stroke"),4))
Outcome <- c(rep(c("Overall Mortality","CVD Incidence","CVD Mortality","CHF Incidence","Stroke"),4))
ASSOC_DF <- data.frame(Outcome,OR,Lower,Upper,Model,Exposure,Pvalue)
# ordering the survival outcome variable for plotting
ASSOC_DF$Outcome <- factor(ASSOC_DF$Outcome, levels = c("Stroke","CHF Incidence","CVD Incidence","CVD Mortality","Overall Mortality"))
# ordering the exposure variable
ASSOC_DF$Exposure <- factor(ASSOC_DF$Exposure, levels = c("Epigenetically Predicted Urinary Arsenic Scores","Measured Urinary Arsenic Levels"))
ASSOC_DF <- ASSOC_DF %>% mutate(Significance=ifelse(Pvalue<1e-3,"**",ifelse(Pvalue<0.05,"*","")))
ASSOC_DF$OR <- round(ASSOC_DF$OR,digits=2)
ASSOC_DF$Lower <- round(ASSOC_DF$Lower,digits=2)
ASSOC_DF$Upper <- round(ASSOC_DF$Upper,digits=2)
# assembling confidence interval string
ASSOC_DF$CI <- paste(sprintf("%.2f",ASSOC_DF$Lower),sprintf("%.2f",ASSOC_DF$Upper),sep = "-")
ASSOC_DF$Association <- paste0(sprintf("%.2f",ASSOC_DF$OR), " (", ASSOC_DF$CI, ")",ASSOC_DF$Significance)
ASSOC_DF$HR <- ASSOC_DF$OR

#define colours for dots and bars
dotCOLS = c("#a6d8f0","#f9b282")
barCOLS = c("#008fd5","#de6b35")

##########################################
# Plotting results from both regression models, only including mortality and CVD-related outcomes
ASSOC_DF <- ASSOC_DF %>% filter(!grepl("CHF Incidence",Outcome)) %>% filter(!grepl("Stroke",Outcome))
p <- ggplot(ASSOC_DF, aes(x=Outcome, y=OR, ymin=Lower, ymax=Upper,col=Exposure,fill=Exposure,label=OR)) + 
  #specify position here
  geom_linerange(size=5,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  #specify position here too
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  scale_x_discrete(name="Survival outcomes") +
  scale_y_continuous(name="Hazard ratio", limits = c(1/2, 2.3)) +
  coord_flip() +
  theme_classic() + 
  theme(legend.position = "bottom", panel.border = element_rect(colour = "black", fill=NA),legend.title=element_blank(),
  axis.text=element_text(size=14),
  axis.title=element_text(size=14,face="bold")) +  
  geom_text(aes(label = Association), color = "black", position = position_dodge(width = 0.5),hjust=-0.35) +
  guides(colour = guide_legend(reverse = TRUE),fill = guide_legend(reverse = TRUE)) +
  facet_wrap(~Model,nrow=1) + 
  theme(strip.text.x = element_text(size = 14,color="black",face="bold"))

# printing out forest plot
png(paste0("FORESTPLOT_Model_Both.png"),units="in",width=8.25,height=4,res=1500)
print(p)
dev.off()





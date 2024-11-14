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

exposure_var<-"log2_sumInMetAsCr"

# importing methylation and covariate data
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data")
load("combined_cov_w_sumInMetAsCr.RData")
combined_cov <- combined_cov_w_sumInMetAsCr
load("combined_beta.RData")
load("combined_mval.RData")
# restricting beta and mval matrix to samples
combined_mval<-combined_mval[,combined_cov$Sample_Name]
combined_beta<-combined_beta[,combined_cov$Sample_Name]

# EWAS MVAL 
mod0 <- model.matrix(~ as.factor(combined_cov$Sex) + combined_cov$Age + as.factor(combined_cov$bmi_cat) + as.factor(combined_cov$cigsmoke) + as.factor(combined_cov$Batch_Plate) )
mod1 <- model.matrix(~ as.factor(combined_cov$Sex) + combined_cov$Age + as.factor(combined_cov$bmi_cat) + as.factor(combined_cov$cigsmoke) + as.factor(combined_cov$Batch_Plate) + combined_cov[,exposure_var] )

svobj_tl_irw = sva(combined_mval,mod1, mod0, method="irw")
    sv <- data.frame(svobj_tl_irw$sv)
    colnames(sv) <- paste("SV", colnames(sv), sep="")
    rownames(sv) <- colnames(combined_mval)
    save(sv,file=paste0("sv_EWAS_HEALS_combined_svaSV_",exposure_var,".RData"))
    V <- sv
    xnam <- paste(" V[,", 1:dim(V)[2], "]", sep="")
    fmla <- as.formula(paste0("~ combined_cov$",exposure_var," + factor(combined_cov$Sex) + combined_cov$Age + factor(combined_cov$bmi_cat) + factor(combined_cov$cigsmoke) + factor(combined_cov$Batch_Plate) + ", paste(xnam, collapse= "+")))
    design <- model.matrix(fmla)
    library(limma)
    fit <- lmFit(combined_mval, design)
    fit <- eBayes(fit)
    results <- topTable(fit, n=dim(combined_mval)[1], sort.by="P",adjust="BH", coef=paste0("combined_cov$",exposure_var), confint=TRUE)
    save(results, file=paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/results_EWAS_HEALS_combined_mval_svaSV_",exposure_var,".RData"))

    #manifest file obtain from this
    # wget https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip
    # unzip infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip 
    manifest<-read.csv("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/annotation/manifest/manifest_4.csv", header=T,sep=",")
    results$Name<-rownames(results)
    results_annotated <- merge(results, manifest, by="Name")
    results_final <- results_annotated[order(results_annotated$P.Value), ]
    output_file<-paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/results_final_EWAS_HEALS_combined_mval_svaSV_",exposure_var,".RData")
    save(results_final,file=output_file)

    
    setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data")
    #########################################################################################################
    # computing inflation factors for including 2 SVs to the max number of SVs that the sva package outputs #
    #########################################################################################################
    # initializing increment list from SV #2 to the max number of SVs, as well as a vector to store inflation factors "inflation_list"
    incre_list<-c(2:ncol(sv))
    inflation_list<-vector(length=length(incre_list))
    num_sig_cpg_list<-vector(length=length(incre_list))
    
    for (incre in 1:length(incre_list)) {
      
      # printing out the number of SVs included for each incremental increase in the number of SVs used 
      incre_value<-incre_list[incre]
      print(incre_value)
      
      # building a model matrix to input into limma for our analysis
      V <- sv[,1:incre_value]
      xnam <- paste(" V[,", 1:dim(V)[2], "]", sep="")
      fmla <- as.formula(paste("~ combined_cov$",exposure_var," + factor(combined_cov$Sex) + combined_cov$Age + factor(combined_cov$bmi_cat) + factor(combined_cov$cigsmoke) + factor(combined_cov$Batch_Plate) + ", paste(xnam, collapse= "+")))
      design <- model.matrix(fmla)
      
      # utilizing limma to run our EWAS for the number of SVs denoted by the "incre_value" variable
      library(limma)
      fit <- lmFit(combined_mval, design)
      fit <- eBayes(fit)
      results <- topTable(fit, n=dim(combined_mval)[1], sort.by="P",adjust="BH", coef=paste0("combined_cov$",exposure_var), confint=TRUE)
      
      # computing the inflation factor for "incre_value" number of SVs and storing it into the "inflation_list" vector that contains SV #2 until the last SV index that the sva package outputs
      pvalue <- results$P.Value
      chisq <- qchisq(1 - pvalue, 1)
      lambda <- median(chisq) / qchisq(0.5, 1)
      inflation_list[incre] <- lambda
      
      # print number of significant hits
      num_sig_cpg_list[incre] <- (sum(results$P.Value < 0.05/nrow(results)))
    }
    
    ############################################
    # computing the inflation factor for SV #1 #
    ############################################
    V<-sv
    xnam=" V[,1]"
    fmla <- as.formula(paste("~ combined_cov$",exposure_var," + factor(combined_cov$Sex) + combined_cov$Age + factor(combined_cov$bmi_cat) + factor(combined_cov$cigsmoke) + factor(combined_cov$Batch_Plate) + ", paste(xnam, collapse= "+")))
    design <- model.matrix(fmla)
    library(limma)
    fit <- lmFit(combined_mval, design)
    fit <- eBayes(fit)
    results <- topTable(fit, n=dim(combined_mval)[1], sort.by="P",adjust="BH", coef=paste0("combined_cov$",exposure_var), confint=TRUE)
    pvalue <- results$P.Value
    chisq <- qchisq(1 - pvalue, 1)
    lambda <- median(chisq) / qchisq(0.5, 1)
    inflation_1_SV <- lambda
    # print number of significant hits
    num_sig_cpg_list <- c((sum(results$P.Value < 0.05/nrow(results))) , num_sig_cpg_list)
    
    #######################################################
    # computing the inflation factor when including 0 SVs #
    #######################################################
    fmla <- as.formula(paste("~ combined_cov$",exposure_var," + factor(combined_cov$Sex) + combined_cov$Age + factor(combined_cov$bmi_cat) + factor(combined_cov$cigsmoke) + factor(combined_cov$Batch_Plate)"))
    design <- model.matrix(fmla)
    library(limma)
    fit <- lmFit(combined_mval, design)
    fit <- eBayes(fit)
    results <- topTable(fit, n=dim(combined_mval)[1], sort.by="P",adjust="BH", coef=paste0("combined_cov$",exposure_var), confint=TRUE)
    pvalue <- results$P.Value
    chisq <- qchisq(1 - pvalue, 1)
    lambda <- median(chisq) / qchisq(0.5, 1)
    inflation_0_SV <- lambda
    # print number of significant hits
    num_sig_cpg_list <- c((sum(results$P.Value < 0.05/nrow(results))) , num_sig_cpg_list)
    
    # combining inflation factors for 0 and 1 SVs with the inflation factor vector containing SV #2 until the last SV index the sva package outputs
    inflation_list_including_0and1<-c(inflation_0_SV,inflation_1_SV,inflation_list)
    incre_list_including_0and1<-c(0:ncol(sv))
    
    #####################################################
    # Plotting the SV index plot with inflation factors #
    #####################################################
    pdf(paste0("inflation_num_SV_",length(incre_list)+1,"_",exposure_var,".pdf"))
    plot(x=incre_list_including_0and1,y=inflation_list_including_0and1,main="Inflation factors for different # of SVs",xlab="Number of SVs",ylab="Inflation Factor")
    dev.off()
    
    ###################################################################################
    # Outputting a data table containing the inflation factors for each number of SVs #
    ###################################################################################
    inflation_dt<-data.table(incre_list_including_0and1,inflation_list_including_0and1)
    write.csv(inflation_dt,file=paste0("inflation_dt_",exposure_var,".csv"),quote=F,row.names=F)
    
    
    ##############################################################
    # Plotting the SV index plot with number of significant cpgs #
    ##############################################################
    pdf(paste0("significant_hits_num_SV_",length(incre_list)+1,"_",exposure_var,".pdf"))
    plot(x=incre_list_including_0and1,y=num_sig_cpg_list,main="Number of significant CpGs for different # of SVs",xlab="Number of SVs",ylab="# of significant CpGs")
    dev.off()
    
    ############################################################################################
    # Outputting a data table containing the number of significant cpgs for each number of SVs #
    ############################################################################################
    num_significant_cpgs_dt<-data.table(incre_list_including_0and1,num_sig_cpg_list)
    write.csv(num_significant_cpgs_dt,file=paste0("num_significant_cpgs_dt_",exposure_var,".csv"),quote=F,row.names=F)
    
    
    
    
    
    
    
    
    
    
    
    
    
# running EWAS with specified number of SVs
    #############################################################
    # EWAS MVAL 
    specifySV=9
    load(paste0("sv_EWAS_HEALS_combined_svaSV_",exposure_var,".RData"))
    V <- sv[,1:specifySV]
    xnam <- paste(" V[,", 1:dim(V)[2], "]", sep="")
    fmla <- as.formula(paste0("~ combined_cov$",exposure_var," + factor(combined_cov$Sex) + combined_cov$Age + factor(combined_cov$bmi_cat) + factor(combined_cov$cigsmoke) + factor(combined_cov$Batch_Plate) + ", paste(xnam, collapse= "+")))
    design <- model.matrix(fmla)
    library(limma)
    fit <- lmFit(combined_mval, design)
    fit <- eBayes(fit)
    results <- topTable(fit, n=dim(combined_mval)[1], sort.by="P",adjust="BH", coef=paste0("combined_cov$",exposure_var), confint=TRUE)
    save(results, file=paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/results_EWAS_HEALS_combined_mval_specifySV_",exposure_var,"_",specifySV,"SV",".RData"))
    
    #manifest file obtain from this
    # wget https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip
    # unzip infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip 
    manifest<-read.csv("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/annotation/manifest/manifest_4.csv", header=T,sep=",")
    results$Name<-rownames(results)
    results_annotated <- merge(results, manifest, by="Name")
    results_final <- results_annotated[order(results_annotated$P.Value), ]
    output_file<-paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/results_final_EWAS_HEALS_combined_mval_specifySV_",exposure_var,"_",specifySV,"SV",".RData")
    save(results_final,file=output_file)
    
    
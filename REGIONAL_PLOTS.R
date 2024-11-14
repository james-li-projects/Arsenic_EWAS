library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)     

# set working directory
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/regional_plots")

# identifying mean DNAm beta values for each CpG
load("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/combined_beta.RData")
cpg_Name <- rownames(combined_beta)
cpg_Avg <- rowMeans(combined_beta)
cpg_beta_df <- data.frame(cpg_Name,cpg_Avg)
colnames(cpg_beta_df) <- c("Name","Avg_Beta")

# import ewas results log2_UrAsgmCr
load("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/results_final_EWAS_HEALS_combined_mval_specifySV_log2_UrAsgmCr_9SV.RData")
results_final <- inner_join(results_final,cpg_beta_df,by=c("Name"))
results_final$CHR <- as.numeric(results_final$CHR)
results_final$MAPINFO <- as.numeric(results_final$MAPINFO)

results_final$Label <- ""
#results_final$Label[results_final$P.Value < 0.05/nrow(results_final)] <- results_final$Name[results_final$P.Value < 0.05/nrow(results_final)]
results_final$Label[results_final$adj.P.Val < 0.05] <- results_final$Name[results_final$adj.P.Val < 0.05]

# adding in the CpG invovled in AxCT interactions
#results_final$Label[results_final$Name=="cg18693395"] <- "cg18693395"
#results_final$Label[results_final$Name=="cg09353563"] <- "cg09353563"

TABLE_S2 <- fread("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/TABLE_S2.tsv") %>% dplyr::select(Name,`Nearest Gene`)
colnames(TABLE_S2)[grepl("Gene",colnames(TABLE_S2))] <- "GeneName"
results_final <- full_join(results_final,TABLE_S2,by=c("Name"))

# specifying CpGs of interest and their window sizes
regional_cpg_name_list <- c(
  "cg19534475",
  "cg10003262",
  "cg18369034",
  "cg04689048",
  "cg05889085",
  "cg05923226",
  "cg05889085",
  "cg14937504",
  "cg07948401"
)
regional_cpg_window_list <- c(
  75000,
  75000,
  25000,
  25000,
  100000,
  5000,
  75000,
  75000,
  75000
)
# creating a data.frame of these CpGs of interest
regional_cpg_df <- data.frame(regional_cpg_name_list,regional_cpg_window_list)
colnames(regional_cpg_df) <- c("Name","Window")

# generating regional plots for CpGs of interest
for(cpg_i in 1:nrow(regional_cpg_df)) {
  regional_cpg_name <- regional_cpg_df$Name[cpg_i]
  print(paste("Plotting CpG",cpg_i,"of",length(regional_cpg_name_list)))
  regional_CHR<-(results_final %>% filter(Name==regional_cpg_name))$CHR
  regional_MAPINFO<-(results_final %>% filter(Name==regional_cpg_name))$MAPINFO
  regional_GeneName<-(results_final %>% filter(Name==regional_cpg_name))$GeneName
  
  regional_results <- results_final %>% filter(CHR==regional_CHR,MAPINFO>(regional_MAPINFO-regional_cpg_df$Window[cpg_i]),MAPINFO<(regional_MAPINFO+regional_cpg_df$Window[cpg_i])) %>% mutate(bp_mb=MAPINFO/1e6) %>% mutate(Sign=ifelse(logFC>0,"Positive","Negative"))
  lower_xlim <- min(regional_results$bp_mb)
  upper_xlim <- max(regional_results$bp_mb)
  
  # plotting out p-value regional plot
  png(paste("PVAL_",regional_GeneName,".png"),units="in",height=2,width=11,res=1500)
  p<-ggplot(data=regional_results, aes(x=`bp_mb`, y=-log10(P.Value),label=Label)) +
    geom_text_repel(box.padding = 0.75, max.overlaps = Inf, min.segment.length=0) +
    geom_point(data=subset(regional_results, Sign=="Negative"), color="#F8766D", size=4) + 
    geom_point(data=subset(regional_results, Sign=="Positive"), color="#619CFF", size=4) +
    geom_point(shape = 1,size = 4,colour = "black") + 
    scale_x_continuous(limits = c(lower_xlim,upper_xlim), expand = c(0, 0)) +
    xlab(paste0("Chromosome ",regional_CHR," Position (Mb)")) + 
    ylab("-log10(p)") + 
    theme_classic() + 
    theme(legend.position = "none", panel.border = element_rect(colour = "black", fill=NA), text = element_text(size = 15)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  print(p)
  dev.off()
  print(paste0(regional_CHR,":",lower_xlim*1e6,"-",upper_xlim*1e6))
  
  # plotting out average beta regional plot
  png(paste("AVGBETA_",regional_GeneName,".png"),units="in",height=2,width=11,res=1500)
  p<-ggplot(data=regional_results, aes(x=`bp_mb`, y=Avg_Beta,label=Label)) +
    geom_text_repel(box.padding = 0.75, max.overlaps = Inf, min.segment.length=0) +
    geom_point(data=subset(regional_results, Sign=="Negative"), color="#F8766D", size=4) + 
    geom_point(data=subset(regional_results, Sign=="Positive"), color="#619CFF", size=4) +
    geom_point(shape = 1,size = 4,colour = "black") + 
    scale_x_continuous(limits = c(lower_xlim,upper_xlim), expand = c(0, 0)) + 
    scale_y_continuous(limits = c(0,1)) +
    xlab(paste0("Chromosome ",regional_CHR," Position (Mb)")) + 
    ylab("Mean DNAm (beta)") + 
    theme_classic() + 
    theme(legend.position = "none", panel.border = element_rect(colour = "black", fill=NA), text = element_text(size = 15)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  print(p)
  dev.off()
}



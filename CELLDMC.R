library(EpiDISH)
library(dplyr)
library(data.table)

data(centDHSbloodDMC.m)

setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data")
load("combined_cov.RData")
load("combined_beta.RData")
load("combined_mval.RData")
exposure_var <- "log2_UrAsgmCr"

# computing CT fractions using EpiDISH (same code as the reference-based EWAS approach)
out.l <- epidish(beta.m = combined_beta, ref.m = centDHSbloodDMC.m, method = 'RPC')
frac.m <- out.l$estF
CT_frac <- data.frame(frac.m)
CT_frac$Sample_Name <- rownames(CT_frac)

pheno.v <- combined_cov$log2_UrAsgmCr
mod1 <- model.matrix(~ as.factor(combined_cov$Sex) + combined_cov$Age + as.factor(combined_cov$bmi_cat) + as.factor(combined_cov$cigsmoke) + as.factor(combined_cov$Batch_Plate))

print("About to run CellDMC!")
celldmc.o <- CellDMC(beta.m = combined_beta, pheno.v = pheno.v, frac.m = frac.m, adjPMethod = "fdr", cov.mod = mod1)
save("celldmc.o",file="/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/celldmc/celldmc.o.RData")
length(which(rowVars((celldmc.o$dmct))>0))

# load("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/celldmc/celldmc.o.RData")
celldmc.o$dmct <- as.matrix(celldmc.o$dmct)

print("log2_UrAsgmCr")
for (i in 2:ncol(celldmc.o$dmct)) {
  print(colnames(celldmc.o$dmct)[i])
  print(table(celldmc.o$dmct[,i]))
}

# parsing DMCT results
celldmc_parsed_table<-data.table(Name=as.character(),dmct_CT=as.character(),Estimate=as.numeric(),p=as.numeric())
cellCT_list <- setdiff(names(celldmc.o$coe),c("Eosino"))
for (currentCT in cellCT_list) {
  print(currentCT)
  currentCT_df <- celldmc.o$coe[[currentCT]] %>% filter(adjP < 0.05) 
  currentCT_df <- currentCT_df %>% mutate(dmct_CT = rep(currentCT,nrow(currentCT_df))) %>% mutate(Name = rownames(currentCT_df)) %>% dplyr::select(Name,dmct_CT,Estimate,p)
  celldmc_parsed_table <- rbind(celldmc_parsed_table, currentCT_df)
}
DMCT_df <- celldmc_parsed_table

# importing manifest
manifest<-read.csv("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/annotation/manifest/manifest_4.csv", header=T,sep=",")

# joining DMCT_df to manifest to add annotations
DMCT_df_annotated <- inner_join(DMCT_df, manifest,by=c("Name"))
DMCT_df_annotated <- DMCT_df_annotated %>% dplyr::select(Name,dmct_CT,CHR,MAPINFO,Estimate,p,Relation_to_UCSC_CpG_Island,UCSC_RefGene_Name)
colnames(DMCT_df_annotated) <- c("Name","AxCT Cell Type","Chromosome", "Position","AxCT Effect Size Estimate","AxCT p-value","CpG Location", "Nearest Gene")

# further annotating the table
DMCT_df_annotated$`CpG Location`[DMCT_df_annotated$`CpG Location`==""]<-"Non-CpG Island"
DMCT_df_annotated$`CpG Location`[DMCT_df_annotated$`CpG Location`=="Island"]<-"Island"
DMCT_df_annotated$`CpG Location`[DMCT_df_annotated$`CpG Location`=="N_Shelf"]<-"Shelf"
DMCT_df_annotated$`CpG Location`[DMCT_df_annotated$`CpG Location`=="S_Shelf"]<-"Shelf"
DMCT_df_annotated$`CpG Location`[DMCT_df_annotated$`CpG Location`=="N_Shore"]<-"Shore"
DMCT_df_annotated$`CpG Location`[DMCT_df_annotated$`CpG Location`=="S_Shore"]<-"Shore"

# annotating nearby genes with bumphunter
library(bumphunter)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

gr=GRanges(seqnames=paste0(rep("chr",nrow(DMCT_df_annotated)),DMCT_df_annotated$Chromosome),
           ranges=IRanges(start=DMCT_df_annotated$Position,end=DMCT_df_annotated$Position)
)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
subject <- annotateTranscripts(txdb, annotationPackage = NULL, by = c("tx","gene"), codingOnly=FALSE, verbose = TRUE, requireAnnotation = FALSE, mappingInfo = NULL, simplifyGeneID = FALSE)
nearest_gene_annotations <- data.table(matchGenes(gr,subject, type = c("any"))) %>% dplyr::select(name,region,distance)
nearest_gene_annotations$region <- as.character(nearest_gene_annotations$region)
colnames(nearest_gene_annotations)[1] <- "geneName"


DMCT_df_annotated <- cbind(DMCT_df_annotated,nearest_gene_annotations)
DMCT_df_annotated$`Nearest Gene` <- nearest_gene_annotations$geneName

DMCT_df_annotated$`Gene Feature`<-rep("",nrow(DMCT_df_annotated))
DMCT_df_annotated$`Gene Feature`[ DMCT_df_annotated$`Gene Feature`==""] <-  DMCT_df_annotated$region[ DMCT_df_annotated$`Gene Feature`==""]
DMCT_df_annotated$`Gene Feature`[grepl("TSS", DMCT_df_annotated$`Gene Feature`)] <-  DMCT_df_annotated$region[grepl("TSS", DMCT_df_annotated$`Gene Feature`)]

DMCT_df_annotated <- DMCT_df_annotated %>% dplyr::mutate(`Gene Feature`=ifelse(`Gene Feature`=="inside","Body",`Gene Feature`)) %>% dplyr::mutate(`Gene Feature`=ifelse(`Gene Feature`=="promoter","Promoter",`Gene Feature`)) %>% dplyr::mutate(`Gene Feature`=ifelse(`Gene Feature`=="upstream","Upstream",`Gene Feature`)) %>% dplyr::mutate(`Gene Feature`=ifelse(`Gene Feature`=="downstream","Downstream",`Gene Feature`)) %>% dplyr::mutate(`Gene Feature`=ifelse(`Gene Feature`=="1stExon","Body",`Gene Feature`)) %>% dplyr::select(-geneName,-region,-distance) # %>% dplyr::select(`Nearest Gene`,`Gene Feature`,geneName,region,distance)

# outputting table
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/celldmc")
write.table(DMCT_df_annotated,file="Supplementary Table Cell DMC.tsv",sep="\t",quote=F,row.names=F)

###############################################
# examining if there are any overlapping CpGs that were discovered in our primary EWAS vs the cell type-specific analysis
load("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/results_final_EWAS_HEALS_combined_mval_specifySV_log2_UrAsgmCr_9SV.RData")
significant_As_CpG <- results_final %>% filter(P.Value < 0.05/nrow(results_final))
intersect(significant_As_CpG$Name,DMCT_df_annotated$Name)
###############################################
# examine if there are any overlapping genes that were discovered in our primary EWAS vs the cell type-specific analysis
significant_As_CpG <- fread("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/TABLE2.tsv",header=T)
intersect(significant_As_CpG$`Nearest Gene`,DMCT_df_annotated$`Nearest Gene`)
DMCT_df_annotated %>% filter(`Nearest Gene` %in% c("PTGDR","UACA"))

###############################################
# exporting all cell type interaction outputs #
###############################################
#celltype_list<-c("B","NK","CD4T","CD8T","Mono","Neutro","Eosino")
results<-celldmc.o$coe$B
results$Name<-rownames(results)
results_annotated <- merge(results, manifest, by="Name")
results_final <- results_annotated[order(results_annotated$p), ]
output_file<-paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/celldmc/supplementary_data_file_celldmc_B_",exposure_var,".tsv")
write.table(results_final,file=output_file,quote=F,sep="\t",row.names=F)

results<-celldmc.o$coe$NK
results$Name<-rownames(results)
results_annotated <- merge(results, manifest, by="Name")
results_final <- results_annotated[order(results_annotated$p), ]
output_file<-paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/celldmc/supplementary_data_file_celldmc_NK_",exposure_var,".tsv")
write.table(results_final,file=output_file,quote=F,sep="\t",row.names=F)


results<-celldmc.o$coe$CD4T
results$Name<-rownames(results)
results_annotated <- merge(results, manifest, by="Name")
results_final <- results_annotated[order(results_annotated$p), ]
output_file<-paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/celldmc/supplementary_data_file_celldmc_CD4T_",exposure_var,".tsv")
write.table(results_final,file=output_file,quote=F,sep="\t",row.names=F)


results<-celldmc.o$coe$CD8T
results$Name<-rownames(results)
results_annotated <- merge(results, manifest, by="Name")
results_final <- results_annotated[order(results_annotated$p), ]
output_file<-paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/celldmc/supplementary_data_file_celldmc_CD8T_",exposure_var,".tsv")
write.table(results_final,file=output_file,quote=F,sep="\t",row.names=F)


results<-celldmc.o$coe$Mono
results$Name<-rownames(results)
results_annotated <- merge(results, manifest, by="Name")
results_final <- results_annotated[order(results_annotated$p), ]
output_file<-paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/celldmc/supplementary_data_file_celldmc_Mono_",exposure_var,".tsv")
write.table(results_final,file=output_file,quote=F,sep="\t",row.names=F)


results<-celldmc.o$coe$Neutro
results$Name<-rownames(results)
results_annotated <- merge(results, manifest, by="Name")
results_final <- results_annotated[order(results_annotated$p), ]
output_file<-paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/celldmc/supplementary_data_file_celldmc_Neutro_",exposure_var,".tsv")
write.table(results_final,file=output_file,quote=F,sep="\t",row.names=F)


results<-celldmc.o$coe$Eosino
results$Name<-rownames(results)
results_annotated <- merge(results, manifest, by="Name")
results_final <- results_annotated[order(results_annotated$p), ]
output_file<-paste0("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output/celldmc/supplementary_data_file_celldmc_Eosino_",exposure_var,".tsv")
write.table(results_final,file=output_file,quote=F,sep="\t",row.names=F)

library(data.table)
library(dplyr)
library(tidyr)

# import ewas results log2_WArsenic
setwd("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/output")
load("results_final_EWAS_HEALS_combined_mval_specifySV_log2_WArsenic_9SV.RData")
results_final$CHR <- as.numeric(results_final$CHR)
results_final$MAPINFO <- as.numeric(results_final$MAPINFO)
results_final <- results_final %>% mutate(Sign = ifelse(logFC>0,"Hypermethylation","Hypomethylation"))
snpsOfInterest_1_fdr <- (results_final %>% filter(adj.P.Val < 0.05,Sign=="Hypermethylation"))$Name
snpsOfInterest_2_fdr <- (results_final %>% filter(adj.P.Val < 0.05,Sign=="Hypomethylation"))$Name
snpsOfInterest_1_bonferroni <- (results_final %>% filter(P.Value < 0.05/nrow(results_final),Sign=="Hypermethylation"))$Name
snpsOfInterest_2_bonferroni <- (results_final %>% filter(P.Value < 0.05/nrow(results_final),Sign=="Hypomethylation"))$Name

# obtaining the p-value equivalent threshold of an FDR of 0.05
fdr_p_thresh <- (results_final %>% filter(adj.P.Val < 0.05)  %>% arrange(P.Value) %>% select(P.Value) %>%
                   tail(1))$P.Value

# make GWAS using ggplot2 to annotate hyper and hypomethylated sites
gwasResults <- results_final %>% select(CHR,MAPINFO,Name,P.Value,adj.P.Val,UCSC_RefGene_Name)
colnames(gwasResults) <- c("CHR","BP","SNP","P","adj.P.Val","Gene")
gwasResults <- gwasResults %>% separate(Gene,sep = ";",into = c("FirstGene"),remove = F) %>% mutate(Gene=FirstGene) %>% select(-FirstGene)

# manually adding ABR gene
gwasResults$Gene[gwasResults$SNP=="cg10003262"] <- "ABR"
gwasResults$Gene[gwasResults$SNP=="cg01912040"] <- "ABR"
gwasResults$Gene[gwasResults$SNP=="cg05962511"] <- "SEMA4G"
gwasResults$Gene[gwasResults$SNP=="cg05428706"] <- "SEMA4G"
gwasResults$Gene[gwasResults$SNP=="cg05889085"] <- "UNC5D"
gwasResults$Gene[gwasResults$SNP=="cg24845774"] <- "UNC5D"
gwasResults$Gene[gwasResults$Gene=="PNMAL2"] <- "PNMA8B"
gwasResults$Gene[gwasResults$SNP=="cg14937504"] <- "AC009518.4"
FOCUS_GENE_ANNOTATE <- c("ATP1B3","ABR","B3GALT5","GBAP1","PNMA8B","CCDC105","AC009518.4")
snpsOfInterest_annotate <- (gwasResults %>% filter(Gene %in% FOCUS_GENE_ANNOTATE, adj.P.Val < 0.05) %>% group_by(Gene) %>% filter(row_number()==1))$SNP
snpsOfInterest_annotate_outline <- (gwasResults %>% filter(Gene %in% FOCUS_GENE_ANNOTATE, adj.P.Val < 0.05) %>% group_by(Gene))$SNP

# List of SNPs to highlight are in the snpsOfInterest object
# We will use ggrepel for the annotation
library(ggrepel)

# Prepare the dataset
don <- gwasResults %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight_1_fdr=ifelse(SNP %in% snpsOfInterest_1_fdr, "yes", "no")) %>%
  mutate( is_highlight_2_fdr=ifelse(SNP %in% snpsOfInterest_2_fdr, "yes", "no")) %>% 
  mutate( is_highlight_1_bonferroni=ifelse(SNP %in% snpsOfInterest_1_bonferroni, "yes", "no")) %>%
  mutate( is_highlight_2_bonferroni=ifelse(SNP %in% snpsOfInterest_2_bonferroni, "yes", "no")) %>% 
  mutate( is_annotate=ifelse(SNP %in% snpsOfInterest_annotate, "yes", "no")) %>%
  # outline 
  mutate( is_annotate_outline=ifelse(SNP %in% snpsOfInterest_annotate_outline, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
p<-ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=1, size=1.3) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,ceiling(-log10(min(don$P)))+1)) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight_1_fdr=="yes"), color="#afcdff", size=2) +
  geom_point(data=subset(don, is_highlight_2_fdr=="yes"), color="#fcbcb8", size=2) +
  geom_point(data=subset(don, is_highlight_1_bonferroni=="yes"), color="#619CFF", size=2) +
  geom_point(data=subset(don, is_highlight_2_bonferroni=="yes"), color="#F8766D", size=2) +
  geom_point(data=subset(don, is_highlight_1_fdr=="yes" & is_annotate_outline=="yes"), shape = 1,size = 2,colour = "black",stroke = 1.1) +
  geom_point(data=subset(don, is_highlight_2_fdr=="yes" & is_annotate_outline=="yes"), shape = 1,size = 2,colour = "black",stroke = 1.1) +
  
  # Add label using ggrepel to avoid overlapping
  geom_text_repel(data=subset(don, is_annotate=="yes"), aes(label=Gene), size=6, point.padding = 3.1, direction="x", fontface="italic") +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
    panel.border = element_rect(colour = "black", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  geom_hline(yintercept = -log10(0.05/nrow(results_final)),color="black") +
  geom_hline(yintercept = -log10(fdr_p_thresh),color="black",linetype = 'dashed') +
  xlab("Chromosome") +
  #coord_cartesian(clip = "off")
  
  # outputting manhattan plot for EWAS of log2_WArsenic
  png("GGPLOT2_MANHATTAN_log2_WArsenic.png",units="in",height=4,width=12,res=1500)
print(p)
dev.off()

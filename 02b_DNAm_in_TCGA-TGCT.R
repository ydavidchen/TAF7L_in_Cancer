##########################################################################################################################
# Visualization of CpG methylation in TCGA-TGCT
# Script author: David Chen
# Original date: 06/29/2017
# Revision for manuscript: 07/30/2017; 08/02/2017
# Notes:
##########################################################################################################################

rm(list=ls())

library(ggplot2)
library(matrixStats)
library(pheatmap)
library(RColorBrewer); gradient_cols <- brewer.pal(12, "Paired")
library(doParallel); registerDoParallel(detectCores() - 1)

load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/TAF7L_in_Testicular_Cancer/062917_TCGA-TGCT_betas.RData")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19);
data(IlluminaHumanMethylation450kanno.ilmn12.hg19);
annot.450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## CpG information:
cpg_info <- annot.450k[grepl("TAF7L",annot.450k$UCSC_RefGene_Name), ]
cpg_info <- data.frame(
  cgID = cpg_info$Name,
  chr  = cpg_info$chr,
  mapInfo = cpg_info$pos,
  strand =cpg_info$strand,
  Island_Name = cpg_info$Islands_Name,
  Relation_to_Island = cpg_info$Relation_to_Island,
  Context = cpg_info$UCSC_RefGene_Group
)

## Update target annotation for heat map:
targets$stage[targets$tumor_stage %in% c("stage i", "stage ia", "stage ib")] <- "Stage I/Ia/Ib";
targets$stage[targets$tumor_stage %in% c("stage ii", "stage iia", "stage iib", "stage iic" )] <- "Stage II/IIa-c";
targets$stage[targets$tumor_stage %in% c("stage iii", "stage iiia", "stage iiib", "stage iiic" )] <- "Stage III/IIIa-c";
targets$stage[targets$tumor_stage == "is"] <- "IS";
targets$stage[targets$tumor_stage == "not reported"] <- "missing";

## Update capitalization of race info:
targets$race[targets$race=="white"] <- "White";
targets$race[targets$race=="black or african american"] <- "Black or African American";
targets$race[targets$race=="asian"] <- "Asian"
targets$race[targets$race=="not reported"] <- "Not reported"

## Load mutation, expression, CN data:
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/TAF7L_in_Testicular_Cancer/062817_TAF7L_TCGA-TGCT.RData")
CN.TAF7L <- as.data.frame(CN.TAF7L);
CN.TAF7L$Sample_ID  <- gsub(".", "-", rownames(CN.TAF7L), fixed=T); 
CN.TAF7L$Sample_ID  <- substr(CN.TAF7L$Sample_ID, 1, 12)
targets <- merge(targets, CN.TAF7L, by="Sample_ID", all.x=TRUE) #1st merge

rsem.TAF7L <- as.data.frame(rsem.TAF7L);
rsem.TAF7L$Sample_ID <- gsub(".", "-", rownames(rsem.TAF7L), fixed=T);
rsem.TAF7L$Sample_ID <- substr(rsem.TAF7L$Sample_ID, 1, 12)
targets <- merge(targets, rsem.TAF7L, by="Sample_ID", all.x=TRUE) #2nd merge

colnames(targets)[colnames(targets) %in% c("TAF7L.x", "TAF7L.y")] <- c("TAF7L.CN", "TAF7L.mRNA")
cut <- median(targets$TAF7L.mRNA, na.rm=TRUE)
targets$TAF7L.expression <- ifelse(targets$TAF7L.mRNA > cut, "above median", "below median")

heat_annot <- data.frame(
  row.names = targets$Sample_ID,
  stage     = targets$stage,
  race      = targets$race,
  TAF7L.mutation = ifelse(targets$Sample_ID == "TCGA-S6-A8JY", "yes", "no"),
  TAF7L.CN = factor(targets$TAF7L.CN),
  TAF7L.expression = targets$TAF7L.expression
);
row_annot <- data.frame(
  row.names = annot.450k$Name,
  TSS  = ifelse(grepl("TSS", annot.450k$UCSC_RefGene_Group), "yes", "no" ),
  Context = annot.450k$Relation_to_Island
)
ann_colors <- list(
  TAF7L.mutation = c(yes="black", no="lightgray"),
  TSS = c(yes="black", no="lightgray"),
  TAF7L.CN = c(`-1`="skyblue",`0`="lightgray",`1`="salmon"),
  TAF7L.expression = c(`above median`="salmon", `below median`="skyblue"),
  race = c(Asian="pink", `Black or African American`="green", `Not reported`="lightgray", White="orange"),
  Context = c(Island="purple", N_Shore=gradient_cols[4], S_Shore=gradient_cols[8], N_Shelf=gradient_cols[3],S_Shelf=gradient_cols[7],OpenSea=gradient_cols[1])
)
sum(rownames(tgct.betas) %in% cpg_info$cgID) #20

## All 20 TAF7L CpGs:
res <- pheatmap(
  tgct.betas[rownames(tgct.betas) %in% cpg_info$cgID, ],
  show_colnames = FALSE, #samples
  annotation_col = heat_annot, 
  annotation_row = row_annot, 
  annotation_colors =  ann_colors,
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  clustering_method = "average",
  border_color = NA,
  color = colorRampPalette(c("blue", "yellow"))(1024),
  fontsize = 7.5, 
  cutree_cols = 2
)
## 13 TAF7L promoter CpGs:
pheatmap(
  tgct.betas[rownames(tgct.betas) %in% cpg_info$cgID[grepl("TSS",cpg_info$Context)], ],
  show_colnames = FALSE, #samples
  annotation_col = heat_annot, 
  annotation_row = row_annot, 
  annotation_colors =  ann_colors,
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  clustering_method = "average",
  border_color = NA,
  color = colorRampPalette(c("blue","yellow"))(1024),
  fontsize = 7.5
)

## Extract group information:
samp.clust.TAF7L <- cutree(res$tree_col, k=2)
samp.clust.df <- data.frame(
  Sample_ID = names(samp.clust.TAF7L),
  methylation.cluster = samp.clust.TAF7L
)
targets <- merge(targets, samp.clust.df, by="Sample_ID")
heat_annot <- data.frame(
  row.names = targets$Sample_ID,
  stage     = targets$stage,
  race      = targets$race,
  TAF7L.mutation = ifelse(targets$Sample_ID == "TCGA-S6-A8JY", "yes", "no"),
  TAF7L.CN = factor(targets$TAF7L.CN),
  TAF7L.expression = targets$TAF7L.expression,
  methylation.cluster = as.factor(targets$methylation.cluster)
);
pheatmap(
  tgct.betas[rownames(tgct.betas) %in% cpg_info$cgID, ],
  show_colnames = FALSE, #samples
  annotation_col = heat_annot, 
  annotation_row = row_annot, 
  annotation_colors =  ann_colors,
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  clustering_method = "average",
  border_color = NA,
  color = colorRampPalette(c("blue", "yellow"))(1024),
  fontsize = 7.5
)

## Fisher test: 
sum(! is.na(targets$TAF7L.mRNA))
my_samples <- targets[!is.na(targets$TAF7L.mRNA), ] #1 sample missing
my_contigency <- matrix(
  c(sum(my_samples$methylation.cluster==1 & my_samples$TAF7L.mRNA > cut),
    sum(my_samples$methylation.cluster==2 & my_samples$TAF7L.mRNA > cut),
    sum(my_samples$methylation.cluster==1 & my_samples$TAF7L.mRNA <= cut),
    sum(my_samples$methylation.cluster==2 & my_samples$TAF7L.mRNA <= cut) ),
  byrow = TRUE, ncol=2
)
my_contigency
fisher.test(my_contigency)
# Fisher's Exact Test for Count Data
# 
# data:  my_contigency
# p-value = 1.331e-11
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.02808799 0.17521256
# sample estimates:
# odds ratio 
# 0.07265449 


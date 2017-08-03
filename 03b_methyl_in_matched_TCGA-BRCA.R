##########################################################################################################################
# TAF7L methylation in TCGA-BRCA primary tumor vs. matched normal adjacent
# Script author: David Chen
# Date: 06/30/2017
# Revision for manuscript: 07/30/2017
# Notes:
# 1. Goal: to provide evidence that TAF7L is an oncogene
##########################################################################################################################

rm(list=ls())

library(ggplot2)
library(matrixStats)
library(pheatmap)
library(RColorBrewer); gradient_cols <- brewer.pal(12, "Paired")
library(reshape2)
library(doParallel); registerDoParallel(detectCores() - 1)

load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/BRCA1ness_TCGA_all_types/063017_TCGA-BRCA_matched.RData")
gdac_profiles <- read.csv("~/Dropbox (Christensen Lab)/Testicular_Cancer_Data_Sets/063017_TAF7L_CN_and_mRNA.csv")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19);
data(IlluminaHumanMethylation450kanno.ilmn12.hg19);
annot.450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19);
cpg_info <- annot.450k[grepl("TAF7L",annot.450k$UCSC_RefGene_Name), ] #subset to TAF7

## Update target annotation for heat map:
targets$stage[targets$pathol_stage %in% c("stage i", "stage ia", "stage ib","stage iia", "stage iib")] <- "Stage I-II";
targets$stage[targets$pathol_stage %in% c("stage iiia", "stage iiib", "stage iiic","stage iv")] <- "Stage III-IV";

## Update capitalization of race info:
targets$race[targets$race=="white"] <- "White";
targets$race[targets$race=="black or african american"] <- "Black or African American";
targets$race[targets$race=="asian"] <- "Asian"

## Update capitalization of tissue context:
targets$tissue.definition <- tolower(targets$tissue.definition)

## Drop unmatched
## and one with duplicate sample: "TCGA-A7-A13G-01B-04D-A22R-05"
ind_unmatched <- which(! duplicated(targets$bcr_patient_barcode) & ! duplicated(targets$bcr_patient_barcode, fromLast=T))
ind_dup <- which(targets$cases == "TCGA-A7-A13G-01B-04D-A22R-05")
targets <- targets[-c(ind_unmatched,ind_dup), ] #drop unmatched
betas_matched <- betas_matched[ , colnames(betas_matched) %in% targets$Sample_ID] #subset

## Update target file with CN & mRNA profiles:
targets <- merge(targets, gdac_profiles, all.x=TRUE, by="Sample_ID")

median(targets$TAF7L.mRNA, na.rm=T)
plot(sort(targets$TAF7L.mRNA))
heat_annot <- data.frame(
  row.names   = targets$Sample_ID,
  stage       = targets$stage,
  race        = targets$race,
  tissue = targets$tissue.definition,
  TAF7L.CN= factor(targets$TAF7L.copy.num),
  TAF7L.expression = ifelse(targets$TAF7L.mRNA > 0, "above 0", "below 0")
);
row_annot <- data.frame(
  row.names = annot.450k$Name,
  TSS  = ifelse(grepl("TSS", annot.450k$UCSC_RefGene_Group), "yes", "no" ),
  context = annot.450k$Relation_to_Island
)
ann_colors <- list(
  TSS = c(yes="black", no="lightgray"),
  TAF7L.CN= c(`-1`="skyblue",`0`="lightgray",`1`="salmon"),
  TAF7L.expression = c(`above 0`="salmon", `below 0`="skyblue"),
  race = c(Asian="pink", `Black or African American`="green", `Not reported`="lightgray", White="orange"),
  tissue = c(`primary solid tumor`="purple", `solid tissue normal`=gradient_cols[4]),
  context = c(Island="purple", N_Shore=gradient_cols[4], S_Shore=gradient_cols[8], N_Shelf=gradient_cols[3],S_Shelf=gradient_cols[7],OpenSea=gradient_cols[1])
)
pheatmap(
  betas_matched[rownames(betas_matched) %in% cpg_info$Name, ],
  show_colnames = FALSE, #samples
  annotation_col = heat_annot, 
  annotation_row = row_annot, 
  annotation_colors =  ann_colors,
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  clustering_method = "average",
  border_color = NA,
  # main = 'Methylation (beta-values) of 20 TAF7L CpGs across 90 pairs of matched normal and primary solid breast tumors',
  color = colorRampPalette(c("blue", "yellow"))(1024),
  fontsize = 7.5
)
pheatmap(
  betas_matched[rownames(betas_matched) %in% cpg_info$Name[grepl("TSS",cpg_info$UCSC_RefGene_Group)], ],
  show_colnames = FALSE, #samples
  annotation_col = heat_annot, 
  annotation_row = row_annot, 
  annotation_colors =  ann_colors,
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  clustering_method = "average",
  border_color = NA,
  # main = 'Methylation (beta-values) of 13 promoter TAF7L CpGs across 90 pairs of matched normal and primary solid breast tumors',
  color = colorRampPalette(c("blue", "yellow"))(1024),
  fontsize = 7.5
)

## Most variable CpGs:
# plot(sort(rowVars(betas_matched)), cex=0.2, pch=16, bty="l);
# abline(h=0.05, lty=2)
# sum(rowVars(betas_matched) >= 0.05)
sele <- order(rowVars(betas_matched), decreasing=TRUE)[1:6000];
pheatmap(
  betas_matched[sele, ],
  show_colnames = FALSE, #samples
  show_rownames = FALSE, #CpGs
  annotation_col = heat_annot, 
  annotation_row = row_annot, 
  annotation_colors =  ann_colors,
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  clustering_method = "average",
  border_color = NA,
  color = colorRampPalette(c("blue", "yellow"))(1024),
  main = 'Most variable 6000 CpGs across 90 pairs of matched normal and primary solid breast tumors'
)

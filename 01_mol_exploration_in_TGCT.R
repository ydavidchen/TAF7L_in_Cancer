##############################################################################
# Exploring CN & mRNA of TAF7L in TCGA-TGCT
# Script author: David Chen
# Original date: 06/28/2017
# Revision for manuscript: 07/30/2017
# Note:
##############################################################################

rm(list=ls())

library(matrixStats)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(WriteXLS)

setwd("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/TAF7L_in_Testicular_Cancer/")

## Somatic mutation:
mut.TAF7L <- read.table("cBioPortal_query/cBio_TAF7L_somatic.txt", row.names=1, header=T, stringsAsFactors=F, strip.white=T)
mut.TAF7L$COMMON <- NULL;
mut.TAF7L[mut.TAF7L=="NaN"] <- NA;
# colnames(mut.TAF7L) <- gsub(".","-",colnames(mut.TAF7L), fixed = TRUE)
colnames(mut.TAF7L) [! is.na(mut.TAF7L[1, ]) ]
## [1] "TCGA.S6.A8JY.01"

## mRNA-seq Z-scores from cBio:
rsem.TAF7L <- read.table("cBioPortal_query/cBio_TAF7L_mRNA_Zscore.txt", row.names=1, header=T, stringsAsFactors=F, strip.white=T)
rsem.TAF7L$COMMON <- NULL;

## CN levels from cBio:
CN.TAF7L <- read.table("cBioPortal_query/cBio_TAF7L_CN.txt", row.names=1, header=T, stringsAsFactors=F, strip.white=T)
CN.TAF7L[CN.TAF7L =="NaN"] <- NA;
CN.TAF7L$COMMON <- NULL;
sum(is.na(CN.TAF7L[1, ]))

## Master molecular data file:
CN.TAF7L <- t(CN.TAF7L);
rsem.TAF7L <- t(rsem.TAF7L)

TAF7L_master <- merge(CN.TAF7L, rsem.TAF7L, by="row.names");
TAF7L_master <- TAF7L_master[complete.cases(TAF7L_master), ];
colnames(TAF7L_master) <- c("Sample","CN.TAF7L", "mRNA.TAF7L");
TAF7L_master$Sample <- gsub('.', '-', TAF7L_master$Sample, fixed=TRUE)
TAF7L_master$mutation <- ifelse(TAF7L_master$Sample=="TCGA-S6-A8JY-01", "yes", "no")

## mRNA expression by copy number:
ggplot(TAF7L_master, aes(x=reorder(Sample, mRNA.TAF7L), y=mRNA.TAF7L, color=mutation)) +
  geom_point(size=1) +
  geom_hline(yintercept=-0.3489, linetype="dashed") +
  theme_classic() +
  ggtitle("Rank of TAF7L gene expression") +
  labs(x="Sample", y="TAF7L mRNA") +
  theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
  annotate("text", x=30, y=-0.2, label="median expression = -0.34")

table(TAF7L_master$CN.TAF7L)
TAF7L_master$CN.TAF7L <- as.factor(TAF7L_master$CN.TAF7L)
plt.TAF7L <- melt(TAF7L_master)
ggplot(plt.TAF7L, aes(x=variable, y=value, fill=CN.TAF7L)) +
  geom_boxplot(outlier.size=0.5) +
  theme_classic() +
  scale_fill_manual(values=c("skyblue","lightgray","salmon")) +
  labs(x="", y="TAF7L mRNA expression") +
  theme(legend.position="top", axis.text.x=element_text(size=10), axis.text.y=element_text(size=10))

TAF7L_master$CN[TAF7L_master$Sample=="TCGA-S6-A8JY-01"]

## Export as RData file: 
# save(list=ls(), file="~/Dropbox (Christensen Lab)/Christensen Lab - 2017/TAF7L_in_Testicular_Cancer/062817_TAF7L_TCGA-TGCT.RData", compress=TRUE)

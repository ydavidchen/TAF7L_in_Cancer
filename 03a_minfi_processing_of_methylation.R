####################################################################################################
# minfi processing of TCGA-TGCT primary tumor data set
# Script author: David Chen <youdinghuan.chen.gr@dartmouth.edu>
# Date: 06/28/2017
# Notes:
####################################################################################################

rm(list=ls())

library(minfi)
library(matrixStats)
library(doParallel); registerDoParallel(detectCores() - 1)

## Load data and sample sheet (saved as a CSV in a minfi-compatible format):
setwd("~/Dropbox (Christensen Lab)/Testicular_Cancer_Data_Sets/Unpacked_450kIdat_TCGT/");
targets <- read.metharray.sheet(getwd());

gdc_primary_tgct <- read.metharray.exp(targets=targets);
gdc_primary_tgct

names    <- pData(gdc_primary_tgct)$Sample_ID;

## Density plots:
# par(mar=c(5,10,4,2));
# densityBeanPlot(gdc_primary_tgct, sampNames=subjects, main="Raw Intensities")
par(mar=c(5,4,4,2));
densityPlot(gdc_primary_tgct, main="All TCGA-TGCT primary solid tumors: raw intensities", cex=0.75)

## Outlier plot: Convert to a MethylSet: 
Mset <- preprocessRaw(gdc_primary_tgct);
Mset <- minfiQC(Mset, fixOutliers=TRUE, verbose=TRUE)
plotQC(Mset$qc);
title("QC of TCGA-TGCT primary solid tumors (N=134)")

## Save initial workspace:
# save(list=c( "gdc_primary_tgct","targets"), file="~/Dropbox (Christensen Lab)/Christensen Lab - 2017/TAF7L_in_Testicular_Cancer/062917_TCGA-TGCT_initial_workspace.RData", compress=TRUE)

#---------------------------------------Quality Control & Normalization---------------------------------------
## Normalization (Funnorm), background correction (methylumi.noob)
normalized450k <- preprocessFunnorm(gdc_primary_tgct);
dim(normalized450k)

## Remove probes failed to meet detection P-value threshold of 0.05:
pvals <- detectionP(gdc_primary_tgct) #from the original dataset
failedP <- (pvals > 0.05) #Be careful with the inequality sign!
summary(failedP)

## Remove failed probes based on the proportion of samples in which they failed:
fraction <- 0.20;
failedProbes <- rownames(failedP)[rowMeans(failedP) > fraction]; ##list of probes
sum(rowMeans(failedP) > fraction)

sum(! rownames(normalized450k) %in% failedProbes)
normalized450kfilter <- normalized450k[! rownames(normalized450k) %in% failedProbes];
dim(normalized450kfilter)

## Remove non-CpGs, control SNP probes, and polymorphic SNP probes:
norm450k.NoSNPs <- dropMethylationLoci(normalized450kfilter); ##drop technical SNP probes & ch. probes
norm450k.NoSNPs <- dropLociWithSnps(norm450k.NoSNPs); ##drop polymorphic SNPs; MAF is set to 0 by default
dim(norm450k.NoSNPs)
# [1] 464757    134

## Visualize the final set:
# par(mar=c(5,10,4,2));
# densityBeanPlot(getBeta(norm450k.NoSNPs), sampNames=names, main="Normalized, all detection P < 0.05, SNP-removed")
par(mar=c(5,4,4,2));
densityPlot(getBeta(norm450k.NoSNPs), main="Normalized, all detection P < 0.05, SNP-removed") 

#---------------------------------------Save beta values & minfi worksheet as R workspace---------------------------------------
tgct.betas <- getBeta(norm450k.NoSNPs);
if( identical(colnames(tgct.betas), targets$Sample_Name) ){
  print("Sample order matches! Proceed to column name updating...");
  colnames(tgct.betas) <- targets$Sample_ID #TCGA barcode
} else {
  stop("Stop! matching is required first!")
}
stopifnot( identical(colnames(tgct.betas), targets$Sample_ID) )
# save(list=c("tgct.betas", "targets"), file="~/Dropbox (Christensen Lab)/Christensen Lab - 2017/TAF7L_in_Testicular_Cancer/062917_TCGA-TGCT_betas.RData", compress=TRUE)


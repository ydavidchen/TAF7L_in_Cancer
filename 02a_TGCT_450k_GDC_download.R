####################################################################################################
# Query & download TCGA TGCT data set
# Script author: David Chen <youdinghuan.chen.gr@dartmouth.edu>
# Date: 06/28/2017
# Notes:
####################################################################################################

library(gdata)
library(DT)
library(doParallel); registerDoParallel(detectCores() - 1)
library(matrixStats)
library(pheatmap)
library(SummarizedExperiment)
library(TCGAbiolinks)

gdc_tgct <- GDCquery_clinic("TCGA-TGCT")
sum(duplicated(gdc_tgct$submitter_id))

setwd("~/Dropbox (Christensen Lab)/Testicular_Cancer_Data_Sets/")
query.meth <- GDCquery(
  project               = "TCGA-TGCT",
  data.category         = "Raw microarray data",
  data.type             = "Raw intensities", 
  experimental.strategy = "Methylation array", 
  legacy                = TRUE,
  file.type             = ".idat",
  platform              = "Illumina Human Methylation 450"
)

# GDCdownload(query.meth)

## In command line (e.g. local or computer cluster)
## cd Raw_intensities/
## find /global/scratch/ydchen/Raw_intensities -name "*.idat" -exec cp -t /global/scratch/ydchen/Unpacked_450kIdat_TCGT/ {} +

## Meta-data associated with downloaded files:
## Important: it includes matched cases-directory
tcgt_IDAT_queried <- getResults(query.meth)

## Save as backup:
# save(
#   list = c("gdc_tgct", "query.meth", "tcgt_IDAT_queried"),
#   file = "~/Dropbox (Christensen Lab)/Testicular_Cancer_Data_Sets/062817_TCGAbiolinks_GDC_query_record.RData",
#   compress = TRUE
# )

load( "~/Dropbox (Christensen Lab)/Testicular_Cancer_Data_Sets/062817_TCGAbiolinks_GDC_query_record.RData")

## Match idat files to patients:
tcgt_IDAT_queried$submitter_id <- substr(tcgt_IDAT_queried$cases, 1, 12);

tcgt_IDAT_queried$vial <- substr(tcgt_IDAT_queried$cases, 14, 16);
table(tcgt_IDAT_queried$vial )

tcgt_IDAT_queried.unique <- tcgt_IDAT_queried[grepl("Grn.idat", tcgt_IDAT_queried$file_name), ];
tcgt_IDAT_queried.unique$file_name <- gsub("_Grn.idat", "", tcgt_IDAT_queried.unique$file_name) #Important drop idat extensions
tcgt_IDAT_queried.unique$updated_datetime <- tcgt_IDAT_queried.unique$created_datetime <- tcgt_IDAT_queried.unique$state <- NULL;

## Identify patients with multiple samples (either different sample types or replicates)
tcgt_IDAT_queried.unique <- tcgt_IDAT_queried.unique[tcgt_IDAT_queried.unique$vial %in% c("01A","01B"), ]
sum(duplicated(tcgt_IDAT_queried.unique$Sample_ID)) 

GDC_tgct_master <- merge(gdc_tgct, tcgt_IDAT_queried.unique, by="submitter_id")

#-------------------------------------------Assemble minfi sample sheets-------------------------------------------
setwd("~/Dropbox (Christensen Lab)/Testicular_Cancer_Data_Sets/Unpacked_450kIdat_TCGT/")
idats <- substr(list.files(), 1, 17)
stopifnot( all(tcgt_IDAT_queried.unique$file_name %in% idats ) ) #checkpoint

## Export sample sheet into directory with IDAT files:
sample_sheet <- data.frame(
  ## Required information for minfi pre-processing:
  Sample_Name       = GDC_tgct_master$file_name,
  Sample_Well       = NA,
  Sample_ID         = GDC_tgct_master$submitter_id,
  Pool_ID           = NA,
  Sentrix_ID        = gsub("_R.*", "", GDC_tgct_master$file_name),
  Sentrix_Position  = gsub(".*_" , "", GDC_tgct_master$file_name),
  
  ## Additional covariates (tumor grade is missing):
  vial              = GDC_tgct_master$vial, 
  tumor_stage       = GDC_tgct_master$tumor_stage,
  morphology        = GDC_tgct_master$morphology,
  race              = GDC_tgct_master$race,
  ethnicity         = GDC_tgct_master$ethnicity,
  age_at_diagnosis  = GDC_tgct_master$age_at_diagnosis,
  year_of_birth     = GDC_tgct_master$year_of_birth,
  vital_status      = GDC_tgct_master$vital_status,
  tissue.definition = GDC_tgct_master$tissue.definition,
  gender            = GDC_tgct_master$gender,
  disease           = GDC_tgct_master$disease
)
# write.csv(sample_sheet, file="TCGA_tgct_minfi_sample_sheet.csv", row.names=FALSE, quote=FALSE)

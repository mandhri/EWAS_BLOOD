
getwd()
setwd("/home/mandhri/Data preprocess/BLOOD dataset analysis")
getwd()
""
# Loading libraries

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(magrittr)
library(HelpersMG)
library(limma)

#Making a directory for the files 
WORKING_DIR="GSE55763"
baseDir<-WORKING_DIR
dir.create("GSE55763")

#### DMPs ####

# Creating function called "create_summary" to put the resulted dmps in a summary
create_summary <- function(toptable = NULL,
                           dataset_label = NULL,
                           directory = getwd())
{
  CPG <- rownames(toptable)
  ALLELE1 <- rep(1,nrow(toptable))
  ALLELE2 <- rep(2,nrow(toptable))
  TESTSTAT <- toptable$t
  PVALUE <- toptable$P.Value
  EFFECTSIZE <- toptable$logFC
  SE <- toptable$SE
  results = data.frame(CPG,
                       ALLELE1,
                       ALLELE2,
                       TESTSTAT,
                       PVALUE,
                       EFFECTSIZE,
                       SE)
  write.table(results,
              file=paste0(directory,"/",dataset_label,".tbl"),
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE,
              sep="\t")
}










getwd("/home/ubuntu/")
setwd("/home/ubuntu")
getwd()


# Loading libraries

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(magrittr)
library(HelpersMG)
library(limma)
library(minfi)
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

#load the data 
phenotypes<-read.table("GSE55763 Phenotypes.txt")
B <- data.table::fread("Data_preprocess")


B <- as.data.frame(B)
dim(B)
B[1:6]
rownames(B) <- B$V1
class(B$V1)

#Drop the column V1 in B.
#Use dplyr:: since R mistake the select function to MASS package for dplyr
B <- B %>% dplyr::select(-V1)

#Calculating M values by getting logistic distribution to base 2
library("minfi")
M <- logit2(B)
M[1:6]
pheno <- read.delim("Data_preprocess")
glimpse(pheno)

## For the dataset GSE55763, cell type composition is not provided. Therefore champ.refbase will be used.

## Correcting cell type composition.

celltypes <- champ.refbase(beta = B,arraytype = "450K")
head(celltypes[[2]])
head(celltypes[[1]])

#Cell type[[1]] contains all the Sentrix IDs where as celltype [[2]], contains all the cell type fractions.
celltypes_1 <- celltypes[[2]]
class(celltypes_1)
head(celltypes_1)
celltypes_1 <- as.data.frame(celltypes_1)
###
celltypes_1$Sample_Name <- rownames(celltypes_1)

#Adding sample _name to the pheno 
pheno <- left_join(pheno, celltypes_1, by = "Sample_Name")

pheno[1:3,2707]
row.names(pheno)
colnames(pheno)
row.names(celltypes_1)
colnames(celltypes_1)

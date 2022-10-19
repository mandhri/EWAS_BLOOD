
getwd("/home/ubuntu/")
setwd("/home/ubuntu/")
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
B <- data.table::fread("GSE55763 beta filtered.txt")
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
pheno <- read.delim("GSE55763 beta filtered.txt")
glimpse(pheno)

## For the dataset GSE55763, cell type composition is not provided. Therefore champ.refbase will be used.

## Correcting cell type composition.

celltypes <- champ.refbase(beta = B,arraytype = "450K")
celltypes <- celltypes[[2]]
celltypes <- as.data.frame(celltypes)
celltypes$Sample_Name <- rownames(celltypes)
pheno <- left_join(pheno, celltypes, by = "Sample_Name")

#write.table(pheno, 
#            "GSE100264 Phenotypes.txt",
#            col.names = TRUE,
#            sep = "\t")

#DMPs
design=model.matrix(~age +
                      hcv_dx +
                      CD4T +
                      CD8T +
                      Mono +
                      NK +
                      Gran,
                    pheno)

fit1_M <- lmFit(M,
                design)
fit2_M <- eBayes(fit1_M)

fit1_B <- lmFit(B,
                design)
fit2_B <- eBayes(fit1_B)

coef = "age"
results <- topTable(fit2_M,
                    coef=coef,
                    number = Inf,
                    p.value = 1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('GSE100264_pvalhist_DMPsCTC.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for age DMPs with CTC",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue",bins = 30)
dev.off()

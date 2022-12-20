
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

B <- data.table::fread("GSE51032 filtered and imputed.txt")
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

pheno<-read.delim("GSE51032 Phenotypes normal samples only.txt")

#No need to remove extra columns. This pheno table is not updated with entrophy values

#Checking the content in the pheno
glimpse(pheno)

## For the dataset GSE55763, cell type composition is not provided. Therefore champ.refbase will be used.
## Correcting cell type composition.

celltypes <- champ.refbase(beta = B,arraytype = "450K")
head(celltypes[[2]])
head(celltypes[[1]])

#Cell type[[1]] contains all the Sentrix IDs where as celltype [[2]], contains all the cell type fractions.
celltypes <- celltypes[[2]]
class(celltypes)
head(celltypes)
celltypes<- as.data.frame(celltypes)
###
celltypes$geo_accession <- rownames(celltypes)
pheno$geo_accession <- as.character(pheno$geo_accession)
pheno <- left_join(pheno, celltypes, by = "geo_accession")



##EXTRA#################
###################################################
#naming the 1st column as simple name
#colnames(pheno1)[1] <- "Sample_Name"

##changing the dataframe rows to columns 
#pheno_m<-t(pheno)
#pheno_m<-as.data.frame(pheno_m)
##Creating a variable 
#x<-sub('.', '', pheno1$Sample_Name)
#pheno_m<-cbind(x,pheno_m)
#rownames(pheno_m) <- NULL

##########################################

#Making the model where males will be represented as 1 and females as 0 in pheno$sex. Since in my analysis, i am
#correcting for cell composition. I will be using different cell types in the model.
design=model.matrix(~age+
                      sex,
                    pheno)

##Linear models for series of Array to differential methylated cpgs/genes 
##identifying differential methylated genes that are associated with phenotype of interest (age).
#In here, the model is getting fitted using empirical Bayes across the differential methylated genes.
#M values

fit1_M <- lmFit(M,
                design)
fit2_M <- eBayes(fit1_M)

names(fit1_M)
names(fit2_M)

#B values
fit1_B <- lmFit(B,
                design)
fit2_B <- eBayes(fit1_B)

#Observing top differentially methyalted cpgs associated with age for both B values and M values
coef = "age"
results <- topTable(fit2_M,
                    coef=coef,
                    number = Inf,
                    p.value = 1)

results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)


#Differential cpg site,
#the log-fold change associated with the comparison we are interested in (DNAM change per unit age)
#the average intensity of the probe set/gene across all chips,
#the (moderated) t-statistic for the hypothesis that the log-fold change is zero (or equivalently, that the fold change is one),
#the associated raw and adjusted p-values for the t-statistic,
##an estimated log-odds ratio for DE.

#Substitute the log FC values of M to the log FC vales of Beta values. In here, log FC values will be normally 
#distributed. 


results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

setwd("/home/mandhri/ewas")
get
#Save p values for distribution of age DMPS with cell type composition in a histogram 
tiff('GSE51032_pvalhist_DMPsCTC.tiff',
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


#Number of DMPs
fdr = 0.005
results_age=topTable(fit2_M,
                     coef = "age",
                     number = nrow(M),
                     adjust.method = "BH",
                     p.value = fdr)
#22108 DMPs 

directory = ("/home/mandhri/ewas/")
create_summary(toptable = results,
               dataset_label = "GSE51032",
               directory = directory)

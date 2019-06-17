## Mehrnoosh Oghbaie
## 06/12/2019
## Define class and run the project

rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

CRAN.packages <- c("readr","readxl","data.table","reshape2","dplyr","magrittr",
                   "igraph","sqldf","stringr","corrplot","ggplot2","R6","ggridges",
                   "gridExtra","ggrepel","rgl","venn", 'writexl')
bioconductor.packages <- c("biomaRt","limma","qvalue","msa","ape","seqinr","ggseqlogo")


############################################################################################
## Check the values in input.info file

source("functions/All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)
source("functions/Data_preparation.R")
source("functions/Imputation.R")
source("functions/Anova.R")
source("functions/Normalize.R")
source("functions/Mutation.R")


run_order <- rbind(cbind("MT302_FLAG","MT302_ORF2"),
                   cbind("PA1_ORF1","PA1_ORF2"), 
                   cbind("MT302_FLAG","PA1_ORF2"),
                   cbind("MT302_ORF2","PA1_ORF2"),
                   cbind("PA1_ORF1","MT302_ORF2"),
                   cbind("PA1_ORF1","MT302_FLAG"))

eLife_list <- read_excel("../Input_data/eLife_list.xlsx", 
                         col_names = FALSE)
colnames(eLife_list) <- c("geneName","uniprotID","Description","Condition")

## Data_preparation.R
PA <- Template$new()
PA$removeContaminant()
PA$logTransformation()
PA$separatedList(run_order)

## Imputation.R
PA$removeAllZeros()
PA$imputeAll()

## Anova.R
PA$anovaAnalysis()
PA$anovaAnalysisPreImpute()
PA$drawVenDiagramSignificant()
PA$drawVolcanoPlot() #or run the following commands at the end (it might crash)


PA$calculateAverageNormalized(run_order)
PA$drawCommonVenndiagram()
PA$drawHeatmap()
PA$drawHeatmapUnnormalized()
PA$drawMutation()

#save(PA,file = "../backup.RData")
save(PA,file = "../backup_v4.RData")

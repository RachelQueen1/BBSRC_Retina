library(Spaniel)
library(scran)
library(scater)
library(batchelor)
library(dplyr)
library(Seurat)
library(SCFunctionsV3)
library(cowplot)
library(RColorBrewer)
library(SCFunctionsV3)
library(harmony)


allData <- "/data/rachel/Linda_Lako/Spatial/spaceranger_outs_2021_102"
tissue <- "fetal_retina"

projDir <- paste0("/data/rachel/Linda_Lako/Spatial/", tissue)
outputDir <- file.path(projDir, "/rObjects/")

## import all files, cluster replicates and and save
sampleNames <- c("8PCW", "13PCW", "12PCW")
for (sampleName in sampleNames){
  source(file.path(projDir, "rscripts/1_load_data.R"))
}


## find markers
source(file.path(projDir,"rscripts/2_Find_Markers.R"))



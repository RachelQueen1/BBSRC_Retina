library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(dplyr)
library(patchwork)
library(rtracklayer)
library(cowplot)


set.seed(1234)

sampleInfo <- readRDS("rObjects/SampleInfo.rds")
sObj_predicted <- readRDS("rObjects/sObj_predictions.rds")


for (i in 1:length(sObj_predicted)){
  sample <- names(sObj_predicted)[i]
  stage <- sampleInfo$Stage[sampleInfo$ID == sample]
  fName <- paste0("csvFiles/", sample, "_", stage, "_markers.csv")
  markers <- FindAllMarkers(sObj_predicted[[sample]]) %>% 
    arrange(cluster, desc(avg_log2FC))
  write.csv(markers, fName)
}
  
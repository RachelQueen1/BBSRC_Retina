library(dplyr)
library(Seurat)
library("ggplot2")
library(monocle3)
library(patchwork)
library(SeuratWrappers)
library(harmony)
library(ggplotify)
library(rmarkdown)



subsets <- readRDS("rObjects/subsets_list.rds")

for (i in 1:9){
  outName <- names(subsets)[i]
  seuratObj_subset <- SetIdent(seuratObj_subset, value = "annotation")
  markers <- FindAllMarkers(seuratObj_subset, only.pos = TRUE, logfc.threshold = 0.5 )
  markers <- markers %>% arrange(cluster, desc(avg_log2FC))
  
  fName <- paste0("csvFiles/markers_annotation_", outName, ".csv")
  write.csv(markers, fName)
  fName <- gsub("csvFiles", "rObjects", fName) %>% gsub(".csv", ".rds", .)
}
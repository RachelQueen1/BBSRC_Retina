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


anno_all <- read.csv("annotation/individual_annotations")
sObj_retina_only <- list()
plot_selected <- list()



for (sample in names(sObj_predicted)){
  anno <- anno_all[anno_all$ID == sample, ]
  sObj <- sObj_predicted[[sample]]
  sObj$CellType <- anno$CellType[match(sObj@active.ident, anno$cluster)]
  plot_selected[[sample]]  <- 
}



for (sample in names(sObj_predicted)){
  anno <- anno_all[anno_all$ID == sample, ]
  sObj <- sObj_predicted[[sample]]
  
  sObj$CellType <- anno$CellType[match(sObj@active.ident, anno$cluster)]
  
  p1 <- DimPlot(sObj, group.by = "CellType") + ggtitle(paste0("Annoatated_", sample))
  
  sObj_retina_only[[sample]] <- sObj[,!grepl("remove", sObj$CellType)]
  
  
  p2 <- DimPlot(sObj_retina_only[[sample]], group.by = "CellType") + ggtitle(paste0("Retina_only_", sample))
  plot_selected[[sample]] <- p1 + p2
}

saveRDS(sObj_retina_only, "rObjects/sObj_retina_only.rds")
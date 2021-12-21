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
  cds_subset_cc <- sobjToCds(seuratObj_subset)
  saveRDS(cds_subset_cc, paste0("rObjects/cds_", outName, ".rds"))
  p1 <- DimPlot(seuratObj_subset, label = TRUE, group.by = "annotation") 
  p2 <- plot_cells(cds_subset_cc, 
                   label_principal_points = TRUE, 
                   color_cells_by = "annotation"
  ) 
  p1 + p2 + ggsave(paste0("StartNodePlots/", outName, ".png"), width = 15, height = 6)
}
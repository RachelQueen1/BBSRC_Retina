## ---Libraries ------
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(dplyr)
library(patchwork)
library(rtracklayer)
library(cowplot)
#library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)

##--- Functions --------
howManyCells <- function(sObj){dim(sObj)[2]}

sObj <- readRDS("rObjects/seurat_PT_RPC_T1_T2_T3.rds")
sObj$CellType <- sObj$CellType %>% as.character()
sObj$CellType[grepl("RPC", sObj$CellType)] <- "RPCs"
sObj <- SetIdent(sObj, value = "CellType")
DefaultAssay(sObj) <- "chromvar"

differential.activity <- FindAllMarkers(
  object = sObj,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

write.csv(differential.activity, "csvFiles/chromvar_enriched_motifs_RPC_T1_T2_T3.csv")
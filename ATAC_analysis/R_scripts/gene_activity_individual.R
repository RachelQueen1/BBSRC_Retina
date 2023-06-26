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

sObj_Filtered <- sObj_cluster <- readRDS("rObjects/sObj_cluster.rds")
sampleInfo <- readRDS("rObjects/SampleInfo.rds")

## add gene activities to atac seq data
for (i in 1:length(sObj_Filtered)){
  sObj <- sObj_Filtered[[i]]
  gene.activities <- GeneActivity(sObj)
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  sObj[['RNA']] <- CreateAssayObject(counts = gene.activities)
  sObj <- NormalizeData(
    object = sObj,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(sObj$nCount_RNA)
  )
  sObj_Filtered[[i]] <- sObj  
}

saveRDS(sObj_Filtered, "rObjects/sObj_cluster.rds")
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(dplyr)
library(patchwork)
library(rtracklayer)
library(cowplot)


sObj_Filtered <- readRDS("rObjects/sObj_cluster.rds")
sampleInfo <- readRDS("rObjects/SampleInfo.rds")
names(sObj_Filtered) <- sampleInfo$ID

## remove "14823" "15424" (no matching RNAseq samples)
sObj_Filtered <- sObj_Filtered[!names(sObj_Filtered) %in% c("14823", "15424")]
sampleInfo <- sampleInfo[!sampleInfo$ID %in% c("14823", "15424"), ]



retina <- readRDS("../developmental_retina/rObjects/sObj_Filtered.rds")
eye <- readRDS("../Eye/rObjects/sObj_Filtered.rds")

anno_df_retina <- read.csv("../HighResFigures/retina_annotations.csv", row.names = 1)
anno_df_eye <- read.csv("annoations._eye_atac.txt")
anno_df <-rbind(anno_df_retina, anno_df_eye)

getSampleNames <- function(sObj){
  sample <- levels(sObj$orig.ident)
}
names(eye) <- sapply(eye, getSampleNames)

retina_short <- retina[names(retina) %in% sampleInfo$ID]
eye_short <- eye[names(eye) %in% sampleInfo$ID]

setdiff(sampleInfo$ID, c(names(eye_short), names(retina_short)))
ref_list <- c(retina_short, eye_short)

## re-order ref list to match sObj_Filtered list
ref_list <- ref_list[names(sObj_Filtered)]

## add anotations
for (i in 1:length(ref_list)){
  sample <- names(ref_list)[i]
  stage <- sampleInfo$Stage[sampleInfo$ID == sample]
  anno <- anno_df[anno_df$Sample == sample, ]

  mm <- match(ref_list[[i]]@active.ident, anno$Cluster)
  ref_list[[i]]$annotation <- anno$Annotation[mm]
}

saveRDS(ref_list, "rObjects/ref_list.rds")



for (i in 1:length(sObj_Filtered)){
  id  <- names(sObj_Filtered)[i]
  sObj_rna <- ref_list[[id]]
  sObj_atac <- sObj_Filtered[[id]]
  DefaultAssay(sObj_atac) <- 'RNA'
  
  
  transfer.anchors <- FindTransferAnchors(
    reference = sObj_rna,
    query = sObj_atac,
    reduction = 'cca'
  )
  
  predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = sObj_rna$annotation,
    weight.reduction = sObj_atac[['lsi']],
    dims = 2:30
  )
  
  sObj_atac <- AddMetaData(object = sObj_atac, 
                           metadata = predicted.labels)
  sObj_Filtered[[id]] <- sObj_atac 
}

saveRDS(sObj_Filtered, "rObjects/sObj_predictions.rds")
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
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)

##--- Functions --------
howManyCells <- function(sObj){dim(sObj)[2]}

sObj_filtered <- readRDS("rObjects/shared_sObj_merge.rds")
# sObj_13PCW <- readRDS("../mulitome/rObjects/sObj_shared_15761_13PCW.rds")
# sObj_filtered <- c(sObj_filtered, sObj_13PCW)

# ## --integrate samples -----
sObj_merge <- merge(sObj_filtered[[1]], sObj_filtered[-1])

topFeat <- FindTopFeatures(object = sObj_merge[['peaks']][])



close <- ClosestFeature(sObj_merge, regions = rownames(sObj_merge))

# retinaDir <- "/data/rachel/Linda_Lako/Retina/"
# sObj_RNA <- readRDS(paste0(retinaDir, "developmental_retina_and_eye_integrated/rObjects/seuratObj_retina_eye_mt10_noHB.rds"))
# vf <- VariableFeatures(sObj_RNA)
# saveRDS(vf, "rObjects/vf.rds")
vf <- readRDS("rObjects/vf.rds")


featuresUse <- close[close$gene_name %in%  vf , "query_region"]


## Normalization and linear dimensional reduction
sObj_merge <- RunTFIDF(sObj_merge)
sObj_merge <- FindTopFeatures(sObj_merge, min.cutoff = 'q0')
sObj_merge <- RunSVD(sObj_merge, features = featuresUse)

elbowDat <- ElbowPlot(sObj_merge, reduction = "lsi", ndims = 50)
corVals <- DepthCor(sObj_merge, n = 50)
componentUse <- corVals$data$Component[abs(corVals$data$nCount_peaks) < 0.4 &
                                         elbowDat$data$stdev > 1.4]



## non linear dimension reduction and clustering (no batch correction)
sObj_merge <- RunUMAP(object = sObj_merge,
                      reduction = 'lsi',
                      dims = componentUse,
                      n.neighbors = 30,
                      min.dist = 0.25,
                      n.epochs = 5000
                      )

DimPlot(sObj_merge, group.by = "CellType", label = TRUE)

sObj_merge <- FindNeighbors(object = sObj_merge, reduction = 'lsi', dims = componentUse)
sObj_merge <- FindClusters(object = sObj_merge, verbose = FALSE, algorithm = 3)

saveRDS(sObj_merge, "rObjects/shared_sObj_merge_with13PCW.rds")

### DA peaks  ------------------------------------------------
sObj_merge <- readRDS("rObjects/shared_sObj_merge_with13PCW.rds")
DA_peaks <- FindAllMarkers(sObj_merge,
                           only.pos = FALSE,
                           min.pct = 0.05,
                           test.use = 'LR',
                           latent.vars = 'peak_region_fragments'
)

saveRDS(DA_peaks, "rObjects/shared_DA_peaks.rds" )


# ### Add motif information ---------------------------------
DefaultAssay(sObj_merge) <- "peaks"
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

sObj_merge <- AddMotifs(
  object = sObj_merge,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

sObj_merge <- RunChromVAR(
  object = sObj_merge,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
saveRDS(sObj_merge,
         "rObjects/shared_sObj_merge_with13PCW.rds")

### Find enriched motifs  ---------------------------------
sObj_merge <- readRDS("rObjects/shared_sObj_merge_with13PCW.rds")

DefaultAssay(sObj_merge) <- "peaks"
DA_peaks <- readRDS("rObjects/shared_DA_peaks.rds")

# get top differentially accessible peaks
DA_peak.sig <- DA_peaks[DA_peaks$p_val < 0.005, ]

noMotifs <- DA_peak.sig$cluster %>% table() %>% data.frame()
colnames(noMotifs) <- c("cluster", "number")
noMotifs <- noMotifs[noMotifs$number > 10,]
clusters <- noMotifs$cluster


DA_peak.sig <- DA_peak.sig[DA_peak.sig$cluster %in% clusters, ]
top.da.peak <- rownames(DA_peak.sig)


# find peaks open in in each cluster
open.peaks <- AccessiblePeaks(sObj_merge, idents = clusters)

# match the overall GC content in the peak set
meta.feature <- GetAssayData(sObj_merge, assay = "peaks", slot = "meta.features")

open <- meta.feature[open.peaks, ]
query <- meta.feature[top.da.peak, ]#
query <- query[!is.na(query$count), ]

## background
peaks.matched <- MatchRegionStats(
  meta.feature = open,
  query.feature = query,
  n = 10000
)


motif_list <- list()

for (cluster in clusters){
  # get top differentially accessible peaks for the cluster
  DA_peaks_cluster <- DA_peak.sig[DA_peak.sig$cluster == cluster,]
  top.da.peak <- rownames(DA_peaks_cluster)
  top.da.peak <- top.da.peak[top.da.peak %in% rownames(sObj_merge)]

  motif_list[[cluster]] <- FindMotifs(
    object = sObj_merge,
    features = top.da.peak,
    background = peaks.matched
  )

}
library(openxlsx)
write.xlsx(motif_list, "csvFiles/shared_enriched_motifs_with13PCW.xls")
saveRDS(motif_list, "rObjects/shared_enriched_motifs_with13PCW.rds")


motif_list_sig <- lapply(motif_list, function(ml){ml[ml$pvalue < 0.05,]})
write.xlsx(motif_list_sig, "csvFiles/shared_enriched_motifs_with13PCW_significant.xls")



### DE genes  ------------------------------------------------
sObj_merge <- readRDS("rObjects/shared_sObj_merge_with13PCW.rds")
DefaultAssay(sObj_merge) <- "RNA"

markers <- FindAllMarkers(sObj_merge,
                          only.pos = TRUE,
                          logfc.threshold = 0.05) %>%
  arrange(cluster, desc(avg_log2FC))
saveRDS(markers, "rObjects/shared_markers_with13PCW.rds")
write.csv(markers, "csvFiles/shared_markers_with13PCW.csv")

##--end------


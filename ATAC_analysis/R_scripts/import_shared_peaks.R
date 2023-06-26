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

# ##--- Original Data -------
sObj_Orig_Retina <- readRDS("rObjects/sObj_retina_only.rds")

# ##--- Setup -----------
set.seed(1234)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
# 
# ##--Sample Info -----------
dataDir <- "shared_peaks"

sampleInfo <- read.csv("csvFiles/sampleInfo", sep = ",", header = TRUE, stringsAsFactors = FALSE)
sampleInfo <- sampleInfo %>% filter(fName!='14611_15PCW_Retina_Shared')


# 
# ##--Load data ------------
# sObj_List <- list()
# 
for (i in 1:nrow(sampleInfo)){
  c_data <- Read10X_h5(filename = (file.path(dataDir, sampleInfo$fName[i], "/filtered_peak_bc_matrix.h5")))
  m_data <- read.csv(file.path(dataDir, sampleInfo$fName[i], "/singlecell.csv"),  header = TRUE,  row.names = 1)
  chrom_assay <- CreateChromatinAssay(counts = c_data,
                                      sep = c(":", "-"),
                                      genome = "hg38",
                                      fragments = (file.path(dataDir,
                                                             sampleInfo$fName[i],
                                                             "/fragments.tsv.gz")) ,
                                      min.cells = 10, min.features = 200)
  sObj_List[[i]] <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = m_data)
  Annotation(sObj_List[[i]]) <- annotations

  ## add sample metadata
  sObj_List[[i]]$tissue <- sampleInfo[i,"Tissue"]
  sObj_List[[i]]$stage <- sampleInfo[i,"Stage"]
  sObj_List[[i]]$ID <- sampleInfo[i,"ID"]
  sObj_List[[i]]$fName <- sampleInfo[i,"fName"]

  ## add QC metrics
  sObj_List[[i]] <- NucleosomeSignal(object = sObj_List[[i]])
  sObj_List[[i]] <- TSSEnrichment(object = sObj_List[[i]], fast = FALSE)
  sObj_List[[i]]$pct_reads_in_peaks <- sObj_List[[i]]$peak_region_fragments / sObj_List[[i]]$passed_filters * 100
  sObj_List[[i]]$blacklist_ratio <- sObj_List[[i]]$blacklist_region_fragments / sObj_List[[i]]$peak_region_fragments



  }




names(sObj_List) <- sampleInfo$ID


saveRDS(sObj_List, "rObjects/shared_sObj_List.rds")

sObj_List <- readRDS("rObjects/shared_sObj_List.rds")

## --- Filter --------
sObj_filtered <- list()
ids <- intersect(names(sObj_List), names(sObj_Orig_Retina))

for (id in ids){
print(id)
sObj_orig <- sObj_Orig_Retina[[id]]
sObj_shared <- sObj_List[[id]]

cellsUse <- intersect(colnames(sObj_orig), colnames(sObj_shared))
sObj_orig <- sObj_orig[,cellsUse]
sObj_shared <- sObj_shared[,cellsUse]
sObj_shared$CellType <- sObj_orig$CellType
sObj_filtered[[id]] <- sObj_shared
}

## update cell numbers
rownames(sampleInfo) <- sampleInfo$ID
sampleInfo$shared_BeforeFiltering <- sapply(sObj_List, howManyCells)
sampleInfo$shared_AfterFiltering <- NA
after <- sapply(sObj_filtered, howManyCells)
sampleInfo[names(after), "shared_AfterFiltering"] <- after

saveRDS(sampleInfo, "rObjects/shared_SampleInfo.rds")





## --- Normalise, Dimension reduction, cluster
for (i in 1:length(sObj_filtered)){
  sObj_filtered[[i]] <- RunTFIDF(sObj_filtered[[i]])
  sObj_filtered[[i]] <- FindTopFeatures(sObj_filtered[[i]], min.cutoff = 'q0')
  sObj_filtered[[i]] <- RunSVD(sObj_filtered[[i]])
  sObj_filtered[[i]] <- RunUMAP(object = sObj_filtered[[i]], reduction = 'lsi', dims = 2:30)
  sObj_filtered[[i]] <- FindNeighbors(object = sObj_filtered[[i]], reduction = 'lsi', dims = 2:30)
  sObj_filtered[[i]] <- FindClusters(object = sObj_filtered[[i]], verbose = FALSE, algorithm = 3)
  }

saveRDS(sObj_filtered, "rObjects/shared_sObj_filtered.rds")




## --- Add Gene Expression  --------
sObj_filtered <- readRDS("rObjects/shared_sObj_filtered.rds")


for (i in 1:length(sObj_filtered)){
  sObj <- sObj_filtered[[i]]
  gene.activities <- GeneActivity(sObj)
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  sObj[['RNA']] <- CreateAssayObject(counts = gene.activities)
  sObj <- NormalizeData(
    object = sObj,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(sObj$nCount_RNA)
  )
  sObj_filtered[[i]] <- sObj
}


## --- save data --------
saveRDS(sObj_filtered, "rObjects/shared_sObj_filtered.rds")

## --integrate samples -----
sObj_merge <- merge(sObj_filtered[[1]], sObj_filtered[-1])

topFeat <- FindTopFeatures(object = sObj_merge[['peaks']][])
close <- ClosestFeature(sObj_merge, regions = rownames(sObj_merge))

retinaDir <- "/data/rachel/Linda_Lako/Retina/"
sObj_RNA <- sObj_Retina <- readRDS(paste0(retinaDir, "developmental_retina_and_eye_integrated/rObjects/seuratObj_retina_eye_mt10_noHB.rds"))
vf <- VariableFeatures(sObj_RNA)


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

saveRDS(sObj_merge, "rObjects/shared_sObj_merge.rds")

### DA peaks  ------------------------------------------------
sObj_merge <- readRDS("rObjects/shared_sObj_merge.rds")
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
# saveRDS(sObj_merge,
#         "rObjects/shared_sObj_merge.rds")

### Find enriched motifs  ---------------------------------
sObj_merge <- readRDS("rObjects/shared_sObj_merge.rds")

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
write.xlsx(motif_list, "csvFiles/shared_enriched_motifs.xls")
saveRDS(motif_list, "rObjects/shared_enriched_motifs.rds")



### DE genes  ------------------------------------------------
sObj_merge <- readRDS("rObjects/shared_sObj_merge.rds")
DefaultAssay(sObj_merge) <- "RNA"

markers <- FindAllMarkers(sObj_merge,
                          only.pos = TRUE,
                          logfc.threshold = 0.05) %>%
  arrange(cluster, desc(avg_log2FC))
saveRDS(markers, "rObjects/shared_markers.rds")
write.csv(markers, "csvFiles/shared_markers.csv")

##--end------














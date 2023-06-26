library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(TFBSTools)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

motifs <- sapply(pfm, function(x){x@name})
motif_df <- data.frame(motif = names(motifs),
                       name = unname(motifs)
)

#sObj <- readRDS("rObjects/shared_sObj_merge_with13PCW.rds")
sObj <- readRDS("rObjects/shared_sObj_merge_with13PCW_footprints.rds")

#### add cell type
anno <- read.csv("annotation/all_samples")

sObj$CellType <- anno$cell_type[match(sObj$peaks_snn_res.0.8, anno$cluster)]

sObj$CellType <- gsub("GABAaergic amacrine cells ", "Gabaergic amacrine cells", sObj$CellType)
levs <- c("early RPCs","late RPCs","T1","T2","T3","RGCs", "HCs", "glycinergic amacrine cells ", "Gabaergic amacrine cells", "starburst amacrine cells", "cones", "rods", "BPs", "microglia", "MG") 
sObj$CellType <- factor(sObj$CellType, levels = levs)

sObj <- SetIdent(sObj,value="CellType")

### read motifs
motifs <- readRDS("rObjects/shared_enriched_motifs_with13PCW_cellTypes.rds")

## top 5
# top_motifs <- lapply(motifs, head, n = 10) %>% 
#   lapply(tail, n = 5)
# lapply(pull, name = "motif.name") %>% 
#   unlist() %>% unname()

## next 5
top_motifs <- lapply(motifs, head, n = 10) %>% 
  lapply(tail, n = 5) %>%
  lapply(pull, name = "motif.name") %>% 
  unlist() %>% unname() %>% unique()


sObj <- Footprint(
  object = sObj,
  motif.name = top_motifs,
  in.peaks = TRUE,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# sObj <- Footprint(
#   object = sObj,
#   motif.name = "BHLHA15(var.2)",
#   in.peaks = TRUE,
#   genome = BSgenome.Hsapiens.UCSC.hg38
# )


saveRDS(sObj, "rObjects/shared_sObj_merge_with13PCW_footprints.rds")

# ggsave(PlotFootprint(sObj, features = motifs[[2]]), "footprint.tiff")
# 
# PlotFootprint(sObj_with_motifswithout_correction_fp, features = motifs[[15]])

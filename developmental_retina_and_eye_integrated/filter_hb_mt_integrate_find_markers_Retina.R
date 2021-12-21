library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)

### Functions
howManyCells <- function(sObj){dim(sObj)[2]}
## downsample to 4000
downsample <- function(sObj, no_cells = 4000){
  if (howManyCells(sObj) >  no_cells){
    set.seed(1234)
    s <- sample(1:howManyCells(sObj), no_cells, replace = FALSE)
    sObj<- sObj[,s]}
  return(sObj)
}

## read in QCed data
seuratObj_list <- readRDS("../fetal/rObjects/sObj_Filtered.rds")


## remove, HB cells, change mitchondrial threshold to 10% and downsample to 4000
sObj_filtered <- list()
## remove hb and set mitochondrial threshold to 10%
dontuse <- c("HBP1", "HBEGF", "HBS1L")



for (i in 1:length(seuratObj_list)){
  
  sObj <- seuratObj_list[[i]]
  HB_sample <- rownames(seuratObj_list[[i]])[grepl("^HB", rownames(seuratObj_list[[i]]))]
  HB_sample <- setdiff(HB_sample, dontuse)
  
  if (length(HB_sample) == 1){
    numberHB <- sObj@assays$RNA@counts[HB_sample, ]
  } else{
    numberHB <- sObj@assays$RNA@counts[HB_sample, ] %>% colSums()
  }
  
  ## remove high mitochondrial
  filter <- sObj$percent.mt < 10 & 
    numberHB == 0
  
  sObj <- sObj[,filter]
  
  
  # downsample to 4000
  sObj_filtered[[i]] <- downsample(sObj)
}


## merge objects
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sObj_ds <- merge(sObj_filtered[[1]], sObj_filtered[-1]) %>% 
  CellCycleScoring(s.features = s.genes, 
                   g2m.features = g2m.genes, 
                   set.ident = FALSE)

sObj_ds$version <- 3
sObj_ds$version[grepl(14149, sObj_ds$sample)]


seuratObj <- sObj_ds %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE, 
            vars.to.regress = c("percent.mt", 
                                "nCount_RNA", 
                                "nFeature_RNA")) %>% 
  RunPCA(pc.genes = seuratObj_ds@var.genes, npcs = 20, verbose = FALSE) %>%
  RunHarmony(c("sample", "version"), plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:10) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = c(0.2,1,2.2))


#seuratObj$stage.2 <- seuratObj$stage
seuratObj$stage <- seuratObj$stage %>% gsub("PCW.*", "PCW", .) 
saveRDS(seuratObj, "rObjects/seuratObj_retina_mt10_noHB.rds")



## Find Markers

markers <- FindAllMarkers(seuratObj, logfc.threshold = 0.7)## DE

markers <- markers %>% dplyr::arrange(cluster, desc(avg_log2FC))

fName <- paste0("csvFiles/", "retina_noHB_MT10", ".csv")
write.csv(markers, fName)
fName <- gsub("csvFiles", "rObjects", fName) %>% gsub(".csv", ".rds", .)
saveRDS(markers, fName)



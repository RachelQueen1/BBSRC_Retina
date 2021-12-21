library(dplyr)
library(Seurat)
library("ggplot2")
library(monocle3)
library(patchwork)
library(SeuratWrappers)
library(harmony)
library(ggplotify)
library(rmarkdown)

sobjToCds <- function(sObj){
  harmonyLoadings <- sObj@reductions[["harmony"]]@feature.loadings.projected
  sObj <- DietSeurat(sObj, graphs = "umap")
  cds <- as.cell_data_set(sObj)
  cds <- cluster_cells(cds, reduction_method = "UMAP", k = 15)
  cds <- learn_graph(cds, use_partition = FALSE)
  ## Add gene names into CDS
  rowData(cds)$gene_name <- rownames(cds)
  rowData(cds)$gene_short_name <- rowData(cds)$gene_name
  ## estimate size factors
  cds <- estimate_size_factors(
    cds,
    round_exprs = TRUE,
    method = c("mean-geometric-mean-total", "mean-geometric-mean-log-total")
  )
  # add harmony reduction
  cds@preprocess_aux$gene_loadings <- harmonyLoadings
  
  
  return(cds)
}

seuratObj <- readRDS("../rObjects/seuratObj_retina_eye_annotated.rds")
cellAnnoWide <- readRDS("rObjects/cellAnnoWide.rds")


seuratObj$consensus_e_l <- cellAnnoWide$consensus %>%  gsub("removed_", "", .) 
seuratObj$consensus <- seuratObj$consensus_e_l %>% gsub("early ", "", .) %>% gsub("late ", "", .) 

DimPlot(seuratObj, group.by = "consensus", label = TRUE)
DimPlot(seuratObj, group.by = "consensus_e_l", label = TRUE)

subsets <- list()
subsets[["RPC_T2_HCs_Acs"]] <- c("RPCs", "T1", "T2",  "HCs", "ACs")
subsets[["RPC_T1_T2_RGC_HC_ACs"]] <- c("RPCs", "T1", "T2", "RGCs", "HCs", "ACs")
subsets[["T1_T3_PR_BP"]] <- c("T1", "T3", "rod precursors", "rods", "cone precursors", "cones" , "BPs")

saveRDS(subsets, "rObjects/subsets_list_consensus.rds")


for (i in 2:5){
  clustersUse <- subsets[[i]]
  outName <- names(subsets)[i]
  
  
  
  
  cc_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
  noCC <- !rownames(seuratObj) %in% cc_genes
  
  seuratObj_subset <- seuratObj[noCC,
                                seuratObj$consensus %in% clustersUse]
  
  if (outName == "Early_Late_RPCs_MG"){
    seuratObj_subset$consensus <- seuratObj_subset$consensus_e_l
  }
  
  seuratObj_subset <- seuratObj_subset %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE, 
              vars.to.regress = c("percent.mt", 
                                  "nCount_RNA", 
                                  "nFeature_RNA",
                                  "S.Score", 
                                  "G2M.Score")) %>% 
    RunPCA(pc.genes = seuratObj_ds@var.genes, npcs = 20, verbose = FALSE) %>%
    RunHarmony(c("sample", "version"), plot_convergence = TRUE) 
  
  
  seuratObj_subset <- seuratObj_subset %>%
    RunUMAP(reduction = "harmony", dims = 1:10) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
    FindClusters(resolution = c(1))
  
  saveRDS(seuratObj_subset, paste0("rObjects/seurat_subset_", outName, ".rds"))
  seuratObj_subset <- SetIdent(seuratObj_subset, value = "consensus")
  markers <- FindAllMarkers(seuratObj_subset, logfc.threshold = 0.7)## DE

  markers <- markers %>% dplyr::arrange(cluster, desc(avg_log2FC))

  fName <- paste0("csvFiles/", outName, ".csv")
  write.csv(markers, fName)
  fName <- gsub("csvFiles", "rObjects", fName) %>% gsub(".csv", ".rds", .)
  saveRDS(markers, fName)
  
  
  
  ### Convert to cds
  cds_subset_cc <- sobjToCds(seuratObj_subset)
  saveRDS(cds_subset_cc, paste0("rObjects/cds_", outName, ".rds"))
  p1 <- DimPlot(seuratObj_subset, label = TRUE, group.by = "consensus") 
  p2 <- plot_cells(cds_subset_cc, 
                   label_principal_points = TRUE, 
                   color_cells_by = "consensus"
  ) 
  p1 + p2 + ggsave(paste0("StartNode/", outName, ".png"), width = 15, height = 6)
  
}


### AFTER setting starting node
## Compile Consensus Pseudotime
for (i in 2:2){
  thetitle <-  names(subsets)[i]
  outName <- names(subsets)[i]
  clustersUse <- subsets[i]
  render('Pseudotime_Template.Rmd', output_file = paste0("consensus_clusters_", thetitle))
  
}







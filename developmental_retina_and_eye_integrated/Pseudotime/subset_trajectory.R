library(dplyr)
library(Seurat)
library("ggplot2")
library(monocle3)
library(patchwork)
library(SeuratWrappers)
library(harmony)
library(ggplotify)
library(rmarkdown)



seuratObj <- readRDS("../rObjects/seuratObj_retina_eye_annotated.rds")


subsets <- list()

annotation <- read.csv("../csvFiles/retina_eye_annotations.txt", sep = "-", header = FALSE)
RPCs <- annotation$V1[annotation$V2 == "RPCs"]
T1 <- annotation$V1[annotation$V2 == "T1"]
T2 <- annotation$V1[annotation$V2 == "T2"]
T3 <- annotation$V1[annotation$V2 == "T3"]
HCs <- annotation$V1[annotation$V2 == "HCs"]
ACs <- annotation$V1[annotation$V2 == "ACs"]
AC_and_HC_precursors <- annotation$V1[annotation$V2 == "AC_and_HC_precursors"]
RGCs <- annotation$V1[annotation$V2 == "RGCs"]
MG <- annotation$V1[annotation$V2 == "MG"]
BP <- annotation$V1[annotation$V2 == "BPs"]

PR <- annotation$V1[annotation$V2 %in% c("rods", "cones", "photoreceptor_precursors")]


subsets <- list()
subsets[["T3_PR"]] <- c(T3, PR)
subsets[["T3_PR_BP"]] <- c(T3, PR, BP)
subsets[["RPC_T1_RGC"]] <- c(RPCs, T1, T2, T3, RGCs)
subsets[["RPC_T2_HC+AC_HCs"]] <- c(RPCs, T1, T2, T3, HCs, AC_and_HC_precursors)
subsets[["RPC_T2_HC+ACS_Acs"]] <- c(RPCs, T1, T2, T3, ACs, AC_and_HC_precursors)
subsets[["EarlyRPC_lateRPC_MG"]] <- c(RPCs, MG)
subsets[["T1_T3_PR_BP"]] <- c(T1, T3, PR, BP)
subsets[["RPC_T2_HC+ACS_HCs_Acs"]] <- c(RPCs, T1, T2, T3, HCs, ACs, AC_and_HC_precursors)
subsets[["RPC_T1_T2_RGC_HC_ACs"]] <- c(RPCs, T1, T2, RGCs, HCs, ACs, AC_and_HC_precursors)

saveRDS(subsets, "rObjects/subsets_list.rds")

forIDs <- c("T1_T3_PR_BP", 
            "RPC_T2_HC+AC_HCs", 
            "RPC_T2_HC+ACS_Acs", 
            "RPC_T1_RGC", 
            "EarlyRPC_lateRPC_MG")





subsets <- subsets[forIDs]

for (i in 1:9){
  
  clustersUse <- subsets[[i]]
  outName <- names(subsets)[i]
  
  
  
  
  cc_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
  noCC <- !rownames(seuratObj) %in% cc_genes 
  
  seuratObj_subset <- seuratObj[noCC,
                                seuratObj@active.ident %in% clustersUse]
  
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
  
  markers <- FindAllMarkers(seuratObj_subset, logfc.threshold = 0.7)## DE

  markers <- markers %>% dplyr::arrange(cluster, desc(avg_log2FC))

  fName <- paste0("csvFiles/", outName, ".csv")
  write.csv(markers, fName)
  fName <- gsub("csvFiles", "rObjects", fName) %>% gsub(".csv", ".rds", .)
  saveRDS(markers, fName)
  
  saveRDS(seuratObj_subset, paste0("rObjects/seurat_subset_part1", outName, ".rds"))
  
  # ### Compile the results
  thetitle <-  names(subsets)[i]
  outName <- names(subsets)[i]
  clustersUse <- subsets[i]
  render('subsets.rmd', output_file = thetitle)
  
  }







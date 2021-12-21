library(dplyr)
library(tidyverse)

seuratObj <- readRDS("../rObjects/seuratObj_retina_eye_annotated.rds")
subsets <- readRDS("rObjects/subsets_list.rds")
forIDs <- c("T1_T3_PR_BP", 
            "RPC_T2_HC+AC_HCs", 
            "RPC_T2_HC+ACS_Acs", 
            "RPC_T1_RGC", 
            "EarlyRPC_lateRPC_MG")

subsets <- subsets[forIDs]


cellAnnotations <- data.frame(cell = NULL,
                              cellType = NULL,
                              subset = NULL)

for(outName in names(subsets)){
annoFile <- paste0("annotations/", outName, "_annotation")

if (file.exists(annoFile)){
  cc_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
  noCC <- !rownames(seuratObj) %in% cc_genes 
  clustersUse <- subsets[[outName]]
  seuratObj_subset <- readRDS(paste0("rObjects/seurat_subset_part1", outName, ".rds"))
  
  
  
    annotations <- read.csv(annoFile)
    annotations$Cluster <- factor(annotations$Cluster)
    seuratObj_subset$annotation <- annotations$Cell_Type[match(seuratObj_subset@active.ident, annotations$Cluster)]
    
    tmp <-data.frame(cell = colnames(seuratObj_subset),
      cellType = seuratObj_subset$annotation,
      subset = outName) %>% filter(!is.na(cellType)) %>% tibble::remove_rownames()
    
    cellAnnotations <- rbind(cellAnnotations, tmp)

}

}

saveRDS(cellAnnotations, "rObjects/cellAnnotations.rds")



# cellAnnoWide <- cellAnnotations %>% 
#   spread(subset, cellType) %>% 
#   tibble::column_to_rownames("cell") %>%
#   .[colnames(seuratObj), ]
#   
# cellAnnoWide[is.na(cellAnnoWide)] <- 0
# 
# 
# ## make the same order as seurat
# cellAnnoWide <- cellAnnoWide[colnames(seuratObj), ]
# seuratObj_cons <- AddMetaData(seuratObj, cellAnnoWide)
# # seuratObj_cons$test <-cellAnnoWide
# # 
# # 
# DimPlot(seuratObj_cons, group.by = "RPC_T2_HC.AC_HCs", label = TRUE)
# DimPlot(seuratObj_cons, group.by = "RPC_T2_HC.ACS_Acs")
# DimPlot(seuratObj_cons, group.by = "RPC_T1_RGC")
# 
# # cellAnnoWide$EarlyRPC_lateRPC_MG <- gsub("mitotic", "removed_mitotic", 
# #                                          cellAnnoWide$EarlyRPC_lateRPC_MG)
# 

cellAnnoWide <- cellAnnotations %>% 
  spread(subset, cellType) %>% 
  tibble::column_to_rownames("cell") %>%
  .[colnames(seuratObj), ]



## make the same order as seurat
cellAnnoWide <- cellAnnoWide[colnames(seuratObj), ]
cellAnnoWide[is.na(cellAnnoWide)] <- 0




## use HC names
cellAnnoWide$consensus <- cellAnnoWide$`RPC_T2_HC+AC_HCs`
seuratObj_cons$consensus <- cellAnnoWide$consensus
DimPlot(seuratObj_cons, group.by = "consensus", label = TRUE)

#Add in ACs
cellAnnoWide$consensus[cellAnnoWide$consensus == 0 & 
                         cellAnnoWide$`RPC_T2_HC+ACS_Acs` == "ACs"] <- "ACs"

seuratObj_cons$consensus <- cellAnnoWide$consensus
DimPlot(seuratObj_cons, group.by = "consensus", label = TRUE)


#Add in RGCs
cellAnnoWide$consensus[cellAnnoWide$consensus == 0 & 
                         cellAnnoWide$RPC_T1_RGC == "RGCs"] <- "RGCs"

seuratObj_cons$consensus <- cellAnnoWide$consensus
DimPlot(seuratObj_cons, group.by = "consensus", label = TRUE)

## Remove removed_RGCs PRs and T3"
cellAnnoWide$consensus[cellAnnoWide$consensus == "removed_PRs"] <- 0
cellAnnoWide$consensus[cellAnnoWide$consensus == "removed_late RPCs"] <- 0
cellAnnoWide$consensus[cellAnnoWide$consensus == "T3"] <- 0
cellAnnoWide$consensus[cellAnnoWide$consensus == "removed_RGCs"] <- 0



## add small RGC back in
cellAnnoWide$consensus[seuratObj@active.ident == 30] <- "RGCs"

seuratObj_cons$consensus <- cellAnnoWide$consensus
DimPlot(seuratObj_cons, group.by = "consensus", label = TRUE)


#Add in MG
cellAnnoWide$consensus[cellAnnoWide$consensus == 0 & 
                         cellAnnoWide$EarlyRPC_lateRPC_MG == "MG"] <- "MG"
seuratObj_cons$consensus <- cellAnnoWide$consensus
DimPlot(seuratObj_cons, group.by = "consensus", label = TRUE)

#Add in T3 PR and BP cells
cellAnnoWide$T1_T3_PR_BP[grepl("removed", cellAnnoWide$T1_T3_PR_BP)] <- 0
cellAnnoWide$consensus[cellAnnoWide$consensus == 0] <- cellAnnoWide$T1_T3_PR_BP[cellAnnoWide$consensus == 0]


## change HC+AC precursors to T2
cellAnnoWide$consensus[cellAnnoWide$consensus == "HC+AC precursors"] <- "T2"


seuratObj_cons$consensus <- cellAnnoWide$consensus
DimPlot(seuratObj_cons, group.by = "consensus", label = TRUE)

cellAnnoWide$consensus[cellAnnoWide$consensus == 0] <- NA

seuratObj_cons$consensus <- cellAnnoWide$consensus
DimPlot(seuratObj_cons, group.by = "consensus", label = TRUE)


saveRDS(cellAnnoWide, "rObjects/cellAnnoWide.rds")






# cellAnnoWide$numberMitotic <- rowSums(cellAnnoWide == "removed_mitotic", na.rm = TRUE)
# sum(cellAnnoWide$numberMitotic == 0)
# 
# cellAnnoWide$numberMitotic == 0
# 
# mitotic <- function(tableColumn){
#   tableColumn <- tableColumn[!is.na(tableColumn)]
#   number_mitotic <- sum(grepl("mitotic", tableColumn))
#   not_mitotic <- sum(!grepl("mitotic", tableColumn), na.rm = TRUE)
#   total_cells <- length(tableColumn)
# }
# 
# 
# 
# tableRow <- cellAnnoWide[1,]
# 
# mostFrequent(tableRow){
# summary <- tableRow %>% .[!is.na(.)] %>% table()
# names(summary)[max(summary)]
# 
# }
# 

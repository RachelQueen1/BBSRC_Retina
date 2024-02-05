library(cisTopic)
library(Signac)
library(Seurat)
library(tidyverse)
library(arrow)



#https://htmlpreview.github.io/?https://github.com/aertslab/cisTopic/blob/master/vignettes/WarpLDA_CompleteAnalysis.html

## read in processed signac data
sObj_ATAC <- readRDS("../rObjects/sObj_with_motifs_fp_rpc.rds")
sObj_ATAC$CellType <- sObj_ATAC$CellType %>% as.character()
sObj_ATAC$CellType[grepl("RPC", sObj_ATAC$CellType)] <- "RPCs"

# cellorder <- c("early RPCs","late RPCs","T1","T2","T3","RGCs", "HCs", "glycinergic amacrine cells ", "Gabaergic amacrine cells", "starburst amacrine cells", "cones", "rods", "BPs", "microglia", "MG")
cellorder <- c("RPCs","T1","T2","T3","RGCs", "HCs", "glycinergic amacrine cells ", "Gabaergic amacrine cells", "starburst amacrine cells", "cones", "rods", "BPs", "microglia", "MG")

## create cisTopic from matrix
counts_AD3 <- sObj_ATAC@assays$ATAC@counts
rownames(counts_AD3) <- rownames(counts_AD3) %>% sub("-", ":", .)
# counts_AD3_feather <- as.data.frame(counts_AD3)
# write_feather(counts_AD3_feather, sink=paste0("feathers/", 'counts_AD3'))


cisTopicObject <- createcisTopicObject(counts_AD3, 
                                       project.name='AD3')

## add metadata
cellData_md <- sObj_ATAC[[]]

cisTopicObject <- addCellMetadata(cisTopicObject, 
                                  cell.data = cellData_md)


## build models
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(16:24), 
                                   seed=123, nCores=5, 
                                   addModels=FALSE, 
                                   tmp = "tmp_models")

saveRDS(cisTopicObject, "cisTopicObject.rds")
# ## Selection of the best model

cisTopicObject <- readRDS("cisTopicObject.rds")

# par(mfrow=c(3,3))
# cisTopicObject <- selectModel(cisTopicObject, type='maximum')
# cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject, type='derivative')
# 
# ## identification of cell states
# cisTopicObject <- runUmap(cisTopicObject, target='cell')
# par(mfrow=c(1,2))
# plotFeatures(cisTopicObject, method='Umap',
#              target='cell',
#              topic_contr=NULL, colorBy=c('CellType'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)

# par(mfrow=c(2,5))
# plotFeatures(cisTopicObject, 
#              method='Umap', target='cell', 
#              topic_contr='Probability', colorBy=NULL, 
#              cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
# 
# cellTopicHeatmap(cisTopicObject, method='Z-score', colorBy=c('CellType'))


# # Export data to python
library(arrow)
path <- "feathers/"
cisTopic_obj <- readRDS("cisTopicObject.rds")
modelMat <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
modelMat <- as.data.frame(modelMat)
write_feather(modelMat, sink=paste0(path, 'model_to_pycisTopic/cell_topic.feather'))
modelMat <- modelMatSelection(cisTopicObject, 'region', 'Probability', all.regions=TRUE)
modelMat <- as.data.frame(modelMat)
write_feather(modelMat, sink=paste0(path, 'model_to_pycisTopic/topic_region.feather'))
## cell data
celldata <- cisTopicObject@cell.data
celldata <- as.data.frame(celldata)
write_feather(celldata, sink=paste0(path, 'model_to_pycisTopic/celldata.feather'))

# You can also save the count.matrix
ct <- cisTopicObject@count.matrix
regions <- rownames(ct)
ct <- as.data.frame(ct)
write_feather(ct, sink=paste0(path, 'model_to_pycisTopic/count_matrix.feather'))
write_feather(regions, sink=paste0(path, 'model_to_pycisTopic/regions.feather'))

write.csv(regions, "regions.txt", quote = FALSE, row.names = FALSE)

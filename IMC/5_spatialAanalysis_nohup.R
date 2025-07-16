# For spatial analysis
## library
library(SpatialExperiment)
library(imcRtools)
library(igraph)

library(ggplot2)
library(ggridges)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(ggraph)
library(viridis)
library(patchwork)

library(dplyr)
library(tidyr)

library(BiocParallel)
library(lisaClust)
library(scales)
library(circlize)

set.seed(619)

date_time <- "0714"

# =============================================================================
# LOAD AND VALIDATE DATA
# =============================================================================
codeSpace <- "/home/lenislin/Experiment/projects/CRLM_2025/IMC"

# Set Working Directory and Source Utilities
workDir <- "/mnt/public/lyx/IMC_HE_Merge/CRLM"
setwd(workDir)

# Define directories for saving results and Steinbock outputs
saveDir <- file.path(workDir, "results")
steinResultDir <- file.path(workDir, "steinbock")
figureDir <- file.path(workDir, "figures","3_Spatial")

if(!dir.exists(figureDir)){
  dir.create(figureDir,recursive = T)
}

# Load your SpatialExperiment object
spe <- readRDS(file.path(saveDir,paste0("subanno_spe_",date_time,".rds")))
img_id_ <- "sample_id"

# 1.Cellular neighborhood analysis

graph_types <- c("knn","delaunay")
k_s <- c(10,20,30,50)

for(graph_type in graph_types){
  if(graph_type == "knn"){
    for(k_ in k_s){
      spe <- buildSpatialGraph(spe, img_id = img_id_, type = "knn", k = k_, name = paste0(graph_type,"_",k_))
    }
  }
  else{
    spe <- buildSpatialGraph(spe, img_id = img_id_, type = graph_type, name =graph_type)
  }
}

# 2.Cluster CN and interaction test
k_clusters <- c(8, 10, 12, 15) # Choose K clusters
Pairnames <- colPairNames(spe)

for(pairname_ in Pairnames){
  
  ## aggregate neighbor
  spe <- aggregateNeighbors(
    spe,ncolPairName = pairname_,
    aggregate_by = "metadata",ncount_by = "sub_celltype"
    )
  
  ## Iterative on number of clusters
  for(k_cluster_ in k_clusters){

    cn <- kmeans(spe$aggregatedNeighbors, centers = k_cluster_)
    
    # 创建动态列名避免覆盖
    cluster_colname <- paste0(pairname_,"_cluster_", k_cluster_)
    colData(spe)[[cluster_colname]] <- as.factor(cn$cluster)
  }
  
  # Interaction analysis
  out <- testInteractions(spe,
                          group_by = img_id_,
                          label = "sub_celltype",
                          colPairName = pairname_,
                          BPPARAM = SerialParam(RNGseed = 619))
  
  saveRDS(out,file.path(saveDir,paste0("Interaction_analysis_out_of_",pairname_,".rds")))
}

saveRDS(spe,file.path(saveDir,paste0("spatial_spe_",date_time,".rds")))

# 3. Patch detection
## calculate the minimum distance between each cell and tumor patch
# detectPatch_List <- list()
# detectPatch_List[["EndoMT"]] <- c("EndoMT","PDL1+ Endomt")  #### EndoMT
# detectPatch_List[["Endothelial"]] <- c("Endothelial") #### Endothelial
# 
# for(patch_name_ in names(detectPatch_List)){
# 
#   ## Assign column
#   colData(spe)[patch_name_] <- FALSE
#   colData(spe)[spe$sub_celltype %in% detectPatch_List[[patch_name_]],patch_name_] <- TRUE
# 
#   ## Patch detection
#   col_name_ = paste0(patch_name_,"_patch")
#   spe <- patchDetection(spe,
#                         patch_cells = colData(spe)[patch_name_][,1],
#                         name = `col_name_`,
#                         img_id = img_id_,
#                         expand_by = 20,
#                         min_patch_size = 10,
#                         colPairName = "knn_20",
#                         BPPARAM = MulticoreParam())
# 
# }

# saveRDS(spe,file.path(saveDir,"spatial_spe_0606.rds"))



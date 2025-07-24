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

date_time <- "0724"

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
spe <- readRDS(file.path(saveDir,paste0("spatial_spe_0722.rds")))
img_id_ <- "sample_id"

## Rename Subtype
subtypes_ <- as.character(spe$sub_celltype)
subtypes_[subtypes_ == "EC_Ki67_CAIX"] <- "EC_CAIX"
spe$sub_celltype <- subtypes_

# 2.Cluster CN + Spatial Context + interaction test
k_clusters <- c(8, 10, 12, 15) # Choose K clusters
Pairnames <- colPairNames(spe)

for(pairname_ in Pairnames){
  
  ## aggregate neighbor
  spe <- aggregateNeighbors(
    spe, colPairName = pairname_,
    aggregate_by = "metadata",count_by = "sub_celltype", name = "aggregatedNeighborhood"
  )
  
  ## Iterative on number of clusters
  for(k_cluster_ in k_clusters){
    
    cn <- kmeans(spe$aggregatedNeighborhood, centers = k_cluster_)
    
    # 创建动态列名避免覆盖
    cluster_colname <- paste0("CN_", pairname_,"_cluster_", k_cluster_)
    colData(spe)[[cluster_colname]] <- as.factor(cn$cluster)
    
    # Spatial context analysis
    spe <- aggregateNeighbors(spe, 
                              colPairName = pairname_,
                              aggregate_by = "metadata",
                              count_by = cluster_colname,
                              name = "aggregatedNeighborhood")
    
    context_colname <- paste0("Context_", pairname_,"_cluster_", k_cluster_)
    spe <- detectSpatialContext(spe, 
                                entry = "aggregatedNeighborhood",
                                threshold = 0.90,
                                name = context_colname)
    
    cat("✓ Calculate CN and Spatial Context of", cluster_colname, "done","\n")
  }
  
  colData(spe) <- colData(spe)[,-match(c("aggregatedNeighborhood"),colnames(colData(spe)))]
  
  # Interaction analysis
  out <- testInteractions(spe,
                          group_by = img_id_,
                          label = "sub_celltype",
                          colPairName = pairname_,
                          BPPARAM = SerialParam(RNGseed = 619))
  
  cat("✓ Perform permutation test on", pairname_, "done","\n")
  
  saveRDS(out,file.path(saveDir,paste0("Interaction_analysis_out_of_",pairname_,".rds")))
}

saveRDS(spe,file.path(saveDir,paste0("spatial_spe_",date_time,".rds")))
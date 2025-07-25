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
spe <- readRDS(file.path(saveDir,paste0("subanno_spe_0714.rds")))
img_id_ <- "sample_id"

## Rename Subtype
if(T){
  subtypes_ <- as.character(spe$sub_celltype)
  subtypes_[subtypes_ == "EC_Ki67_CAIX"] <- "EC_CAIX"
  spe$sub_celltype <- subtypes_
}

# 1.Cellular neighborhood analysis
graph_types <- c("knn","delaunay")
k_s <- c(10,20,40)

for(graph_type in graph_types){
  if(graph_type == "knn"){
    for(k_ in k_s){
      spe <- buildSpatialGraph(spe, img_id = img_id_, type = "knn", k = k_, name = paste0(graph_type,"_",k_))
      
      cat("✓ Construct graph through:", paste0(graph_type,"_",k_), "done","\n")
    }
  }
  else{
    spe <- buildSpatialGraph(spe, img_id = img_id_, type = graph_type, name =graph_type)
    cat("✓ Construct graph through:", graph_type, "done","\n")
  }
}

saveRDS(spe,file.path(saveDir,paste0("spatial_spe_",date_time,".rds")))

# 2.Cluster CN + interaction test
k_clusters <- c(10, 12) # Choose K clusters
print(colPairNames(spe))

Pairnames_for_CN <- c("delaunay","knn_10","knn_20")
cn_types <- c()

for(pairname_ in Pairnames_for_CN){
  
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
    
    cn_types <- c(cn_types,cluster_colname) ## combine for Spatial Context analysis
    cat("✓ Calculate CN of", cluster_colname, "done","\n")
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

# 3. Spatial Context
Pairnames_for_SC <- c("knn_40")
cn_types <- c("CN_knn_20_cluster_10", "CN_knn_20_cluster_12")

for(pairname_ in Pairnames_for_SC){
  for(cn_type_ in cn_types){
    
    # Compute the fraction of cellular neighborhoods around each cell
    spe <- aggregateNeighbors(
      spe, colPairName = pairname_,
      aggregate_by = "metadata",
      count_by = cn_type_, 
      name = "aggregatedNeighborhood"
    )
    
    # Detect spatial contexts
    context_colname <- paste0("Context_graph_", pairname_,"_from_", cn_type_)
    spe <- detectSpatialContext(spe, 
                                entry = "aggregatedNeighborhood",
                                threshold = 0.90,
                                name = context_colname)
    
    cat("✓ Calculate Spatial Context of", pairname_, "under", cn_type_,"done","\n")
    colData(spe) <- colData(spe)[,-match(c("aggregatedNeighborhood"),colnames(colData(spe)))]
  }
}

saveRDS(spe,file.path(saveDir,paste0("spatial_spe_",date_time,".rds")))

# 4. Patch detection
## calculate the minimum distance between each cell and tumor patch
detectPatch_List <- list()
celltypes <- unique(spe$sub_celltype)

detectPatch_List[["Tumor"]] <- celltypes[startsWith(celltypes,prefix = "EC")]
detectPatch_List[["Metabolism_activate_Tumor"]] <- c("EC_GLUT1","EC_CAIX","EC_Vimentin")
detectPatch_List[["Quiescent_Tumor"]] <- "EC_EpCAM"

for(patch_name_ in names(detectPatch_List)){

  ## Assign column
  colData(spe)[patch_name_] <- FALSE
  colData(spe)[spe$sub_celltype %in% detectPatch_List[[patch_name_]],patch_name_] <- TRUE

  ## Patch detection
  col_name_ = paste0(patch_name_,"_patch")
  spe <- patchDetection(spe,
                        patch_cells = colData(spe)[patch_name_][,1],
                        name = `col_name_`,
                        img_id = img_id_,
                        expand_by = 20,
                        min_patch_size = 5,
                        colPairName = "knn_20",
                        BPPARAM = MulticoreParam())
  
  cat("✓ Perform ", patch_name_, "patch detection on", "knn_20", "done","\n")
}

saveRDS(spe,file.path(saveDir,paste0("spatial_spe_",date_time,".rds")))



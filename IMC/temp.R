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
date_time <- "0728"
spe <- readRDS(file.path(saveDir,paste0("spatial_spe_","0724",".rds")))
img_id_ <- "sample_id"

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
                        BPPARAM = MulticoreParam(workers = 4))
  
  cat("✓ Perform ", patch_name_, "patch detection on", "knn_20", "done","\n")
}

saveRDS(spe,file.path(saveDir,paste0("spatial_spe_",date_time,".rds")))



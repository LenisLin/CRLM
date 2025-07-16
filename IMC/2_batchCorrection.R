## clustering and Annotation

library(imcRtools)
library(cytomapper)
library(SpatialExperiment)
library(scater)
library(CATALYST)
library(dplyr)
library(batchelor)
library(bluster) 
library(BiocParallel)
library(scran)

## plot
library(RColorBrewer)
library(dittoSeq)
library(viridis)
library(patchwork)
library(cowplot)
library(scuttle)
library(ggplot2)
library(ggrepel)

set.seed(619)

# ------------------------------
# Set Working Directory and Source Utilities
# ------------------------------
workDir <- "/mnt/public/lyx/IMC_HE_Merge/CRLM"

# Define directories for saving results and Steinbock outputs
setwd(workDir)
source("./code/utils.R")
source("./code/visualize.R")

saveDir <- file.path(workDir, "results")
steinResultDir <- file.path(workDir, "steinbock")
figureDir <- file.path(workDir, "figures","1_Anno")

# Create the results directory if it doesn't exist
if(!dir.exists(figureDir)){
  dir.create(figureDir,recursive = T)
}

# Load data
spe <- readRDS(file.path(saveDir,"norm_spe.rds"))

## Add color for patient ID
if(T){
  metadata(spe)$color_vectors[["patient_id"]] <- metadata(spe)$color_vectors$color_multi[1:length(unique(spe$patient_id))]
}

# Perform batch correction 
set.seed(619)
out <- fastMNN(spe, batch = spe$Batch,
               auto.merge = TRUE, subset.row = rowData(spe)$use_channel, assay.type = "exprs") 
stopifnot(all.equal(colnames(spe), colnames(out))) # Check that order of cells is the same

# Store corrected embeddings in SPE object 
reducedDim(spe, "fastMNN") <- reducedDim(out, "corrected")
spe <- runUMAP(spe, dimred= "fastMNN", name = "UMAP_mnnCorrected") 

# Save objects
saveRDS(spe, file = file.path(saveDir,"spe_correct.rds"))

# visualize patient id 
if(T){
  # Check for Batch effect
  spe_subset <- spe[,sample.int(n=ncol(spe),size = 5e4, replace = F)]

  spe_subset <- runUMAP(spe_subset, subset_row = (rowData(spe_subset)$use_channel), exprs_values = "exprs", name = "UMAP")
  spe_subset <- runUMAP(spe_subset, dimred= "fastMNN", name = "UMAP_mnnCorrected") # Compute UMAP on corrected embeddings 
  
  p1 <- dittoDimPlot(spe_subset, var = "Batch", 
                     reduction.use = "UMAP", size = 0.2) + 
    #scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
    ggtitle("Patient ID on UMAP before correction")
  p2 <- dittoDimPlot(spe_subset, var = "Batch", 
                     reduction.use = "UMAP_mnnCorrected", size = 0.2) + 
    #scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
    ggtitle("Patient ID on UMAP after correction")
  
  pdf(file.path(figureDir,"Batch_correction_after.pdf"),height = 8,width = 12)
  plot_grid(p1, p2)
  dev.off()
  
  rm(spe_subset,p1,p2) # skip the batch correction
  gc()
}
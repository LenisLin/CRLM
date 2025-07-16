## Enhanced Spatial Sequencing Data Visualization
library(imcRtools)
library(cytomapper)
library(RColorBrewer)
library(dittoSeq)
library(viridis)
library(ggplot2)
library(patchwork)
library(scales)

# Setup directories
workDir <- "~/IMCDataAnalysis/horizontalProject/IMC_20240725_63"
saveDir <- file.path(workDir,"data/")  
figuresDir <- file.path(workDir,"figures/")
RestulsDir <- file.path(workDir, "results")
setwd(workDir)
source("./code/utils.R")

# Enhanced directory creation with better organization
dirs_to_create <- c(saveDir, figuresDir, 
                    file.path(figuresDir, "high_res"),
                    file.path(figuresDir, "overview"),
                    file.path(figuresDir, "individual_samples"))

for(dir in dirs_to_create){
  if(!dir.exists(dir)){
    dir.create(dir, recursive = TRUE)
  }
}

## Load data
spe <- readRDS(file.path(saveDir,"spatial_spe_0424.rds"))

## Enhanced sample selection strategy
sample_ids <- unique(spe$sample_id)
# Select samples with good representation of cell types
samples_toPlot_ids <- sample_ids[sample(length(sample_ids), min(7, length(sample_ids)))]

# Setup image directories
imgDir <- file.path(workDir,"steinbock","img")
maskDir <- file.path(workDir,"steinbock","masks")
imgTemp <- file.path(workDir,"steinbock","img_select")
maskTemp <- file.path(workDir,"steinbock","masks_select")

if(!file.exists(imgTemp)){
  dir.create(imgTemp, recursive = TRUE)
  dir.create(maskTemp, recursive = TRUE)
}

## Copy images and masks
for(sample_id_ in samples_toPlot_ids){
  ori_img_path <- file.path(imgDir, paste0(sample_id_,".tiff"))
  target_img_path <- file.path(imgTemp, paste0(sample_id_,".tiff"))
  
  ori_mask_path <- file.path(maskDir, paste0(sample_id_,".tiff"))
  target_mask_path <- file.path(maskTemp, paste0(sample_id_,".tiff"))
  
  file.copy(from = ori_img_path, to = target_img_path, overwrite = TRUE, copy.mode = TRUE)  
  file.copy(from = ori_mask_path, to = target_mask_path, overwrite = TRUE, copy.mode = TRUE)
}

## Load images and masks
speTemp <- spe[,spe$sample_id %in% samples_toPlot_ids]
images <- loadImages(imgTemp)
masks <- loadImages(maskTemp, as.is = TRUE)

## Setup metadata
channelNames(images) <- rownames(spe)
all.equal(names(images), names(masks))
mcols(images) <- mcols(masks) <- DataFrame(sample_id = names(images))

## Enhanced normalization with better contrast
images <- cytomapper::normalize(images, separateImages = TRUE)
images <- cytomapper::normalize(images, inputRange = c(0, 0.2))

# ======================================================
# 1. ENHANCED SEGMENTATION VISUALIZATION
# ======================================================

# Create a more sophisticated color palette for markers
marker_colors <- list(
  CD20 = c("black", "#E31A1C"),      # Bright red for B cells
  CD3 = c("black", "#1F78B4"),       # Blue for T cells  
  E_Cadherin = c("black", "#33A02C"), # Green for epithelial
  Pan_CK = c("black", "#FF7F00"),     # Orange for cytokeratin
  VEGFR2 = c("black", "#6A3D9A"),     # Purple for endothelial
  DNA1 = c("black", "blue")       # Light blue for nuclei
)

# High-resolution segmentation plot
pdf(file.path(figuresDir, "segmentation_visualization_enhanced.pdf"), 
    height = 16, width = 25, useDingbats = FALSE)
plotPixels(
  images, 
  mask = masks, 
  img_id = "sample_id", 
  missing_colour = "white", 
  colour_by = c("CD20", "CD3", "E_Cadherin", "Pan_CK", "VEGFR2", "DNA1"), 
  colour = marker_colors,
  image_title = list(position = "top", font = 2, cex = 1.2),
  legend = list(
    colour_by.title.cex = 1.1, 
    colour_by.labels.cex = 1.0
  ),
  thick = TRUE
)
dev.off()

# ======================================================
# 2. ENHANCED CELL TYPE VISUALIZATION
# ======================================================

# Create a sophisticated color palette for cell types
unique_celltypes <- unique(spe$sub_celltype)
n_celltypes <- length(unique_celltypes)

# Use a combination of qualitative palettes for better distinction
if(n_celltypes <= 12) {
  color.panel <- setNames(brewer.pal(min(n_celltypes, 12), "Set3"), unique_celltypes)
} else {
  # For more cell types, use a combination of palettes
  color.panel <- setNames(
    c(brewer.pal(12, "Set3"), 
      brewer.pal(min(n_celltypes-12, 8), "Dark2"),
      rainbow(max(0, n_celltypes-20)))[1:n_celltypes], 
    unique_celltypes
  )
}

# High-resolution cell type visualization
pdf(file.path(figuresDir, "celltype_visualization_enhanced.pdf"), 
    height = 16, width = 25, useDingbats = FALSE)
plotCells(
  mask = masks,
  object = speTemp, 
  cell_id = "ObjectNumber", 
  img_id = "sample_id",
  colour_by = "sub_celltype",
  colour = list(sub_celltype = color.panel),
  exprs_values = "exprs",
  image_title = list(position = "top", font = 2, cex = 1.2),
  legend = list(
    colour_by.title.cex = 1.1, 
    colour_by.labels.cex = 0.9
  ),
  missing_colour = "grey90"
)
dev.off()

# ======================================================
# 3. INDIVIDUAL SAMPLE DETAILED VIEWS
# ======================================================

# Create individual high-quality plots for each sample
for(sample_id in samples_toPlot_ids) {
  
  # Subset data for current sample
  current_spe <- speTemp[, speTemp$sample_id == sample_id]
  current_images <- images[names(images) == sample_id]
  current_masks <- masks[names(masks) == sample_id]
  
  # Individual marker visualization
  png(file.path(figuresDir, "individual_samples", paste0("markers_", sample_id, ".png")), 
      height = 2000, width = 3000, res = 300)
  
  plotPixels(
    current_images, 
    mask = current_masks, 
    img_id = "sample_id", 
    missing_colour = "white", 
    colour_by = c("CD20", "CD3", "E_Cadherin", "Pan_CK", "VEGFR2", "DNA1"), 
    colour = marker_colors,
    image_title = list(text = paste("Sample:", sample_id), position = "top", font = 2, cex = 1.4),
    legend = list(colour_by.title.cex = 1.2, colour_by.labels.cex = 1.0),
    thick = TRUE
  )
  dev.off()
  
  # Individual cell type visualization
  png(file.path(figuresDir, "individual_samples", paste0("celltypes_", sample_id, ".png")), 
      height = 2000, width = 3000, res = 300)
  
  plotCells(
    mask = current_masks,
    object = current_spe, 
    cell_id = "ObjectNumber", 
    img_id = "sample_id",
    colour_by = "sub_celltype",
    colour = list(sub_celltype = color.panel),
    exprs_values = "exprs",
    image_title = list(text = paste("Cell Types -", sample_id), position = "top", font = 2, cex = 1.4),
    legend = list(colour_by.title.cex = 1.2, colour_by.labels.cex = 1.0),
    missing_colour = "grey90"
  )
  dev.off()
}
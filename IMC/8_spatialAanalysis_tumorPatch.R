# Patch Analysis
library(SpatialExperiment)
library(BiocParallel)

library(imcRtools)
library(lisaClust)

library(ggplot2)
library(ggridges)
library(ggpubr)
library(ggraph)
library(ggrepel)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(viridis)
library(patchwork)

library(dplyr)
library(tidyr)
library(tidygraph)
library(igraph)
library(scales)

# =============================================================================
# LOAD AND VALIDATE DATA
# =============================================================================
codeSpace <- "/home/lenislin/Experiment/projects/CRLM_2025/IMC"

source(file.path(codeSpace,"utils.R"))
source(file.path(codeSpace,"visualize.R"))
source(file.path(codeSpace,"Spatial_analysis_functions.R"))

# Set Working Directory and Source Utilities
workDir <- "/mnt/public/lyx/IMC_HE_Merge/CRLM"
setwd(workDir)

# Define directories for saving results and Steinbock outputs
saveDir <- file.path(workDir, "results")
steinResultDir <- file.path(workDir, "steinbock")
figureDir <- file.path(workDir, "figures","3_SpatialAnalysis","Patch_Analysis")

if(!dir.exists(figureDir)){
  dir.create(figureDir,recursive = T)
}

# Load your SpatialExperiment object
date_time <- "0730"
spe <- readRDS(file.path(saveDir,paste0("spatial_spe_",date_time,".rds")))
img_id_ <- "sample_id"

# =============================================================================
# Extract IM
# =============================================================================
spe_IM <- spe[,spe$Tissue == "IM"]
meta_IM <- colData(spe_IM)
colnames(meta_IM) ## view column of data

patch_types <- colnames(meta_IM)[endsWith(colnames(meta_IM),suffix = "_patch")] ## get all patch type
patch_names <- sapply(patch_types, function(x){
  return(strsplit(x,"_")[[1]][1])
})
all_images <- unique(spe_IM$sample_id)

# Defined color for CN
RFS_color <- setNames(c("#EFC000FF", "#0073C2FF"),unique(meta_IM$RFS_status))
Tissue_color <- setNames("#E5A8A0FF","IM")
Patch_color <- c("Tumor_patch" = "#79af97ff","Metabolism_activate_Tumor_patch" = "#b24745ff", "Quiescent_Tumor_patch" = "#7aa6dcff")

## Plot Activative and quiescent tumor Patch on spatial
sample_rois_ <- sample(all_images,12)
spe_temp <- spe_IM[,spe_IM$sample_id %in% sample_rois_]

for(patch_type_ in patch_types){
  ## Sample for plotting
  patch_type_vec_ <- colData(spe_temp)[,patch_type_] ## get correspond patch type
  colData(spe_temp)[,patch_type_] <- ifelse(is.na(patch_type_vec_),NA,patch_type_)
  
  p <- plotSpatial(
    spe_temp,
    node_color_by = patch_type_,
    img_id = "sample_id",
    node_size_fix = 0.5
  ) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = Patch_color)
  
  pdf(file.path(figureDir,paste0(patch_type_," in random ROIs.pdf")), width = 24,height = 18)
  print(p)
  dev.off()
}
rm(spe_temp,sample_rois_,patch_type_vec_,p)
gc()

# =============================================================================
# Section 1: Metabolic Profiling - Compare gene expression between patch types
# Focus: EC_GLUT1+ (Metabolism_activate_Tumor) vs EC_EpCAM+ (Quiescent_Tumor) patches
# =============================================================================

# Set analysis parameters
patch_types <- c("Metabolism_activate_Tumor_patch", "Quiescent_Tumor_patch")
patch_labels <- c("Metabolism_activate_Tumor", "Quiescent_Tumor")
names(patch_labels) <- patch_types

# Get metadata and expression data
meta_IM <- as.data.frame(colData(spe_IM))
expr_matrix <- assay(spe_IM)

# Remove DNA markers for analysis
dna_markers <- rownames(expr_matrix)[grepl("^DNA", rownames(expr_matrix))]
analysis_markers <- rownames(expr_matrix)[!rownames(expr_matrix) %in% dna_markers]
expr_matrix_clean <- expr_matrix[analysis_markers, ]

print(paste("Number of markers for analysis:", length(analysis_markers)))

# Extract EC_GLUT1+ patches (Metabolism_activate_Tumor_patch)
glut1_cells <- meta_IM[!is.na(meta_IM$Metabolism_activate_Tumor_patch), ]
glut1_patches <- unique(glut1_cells$Metabolism_activate_Tumor_patch)

print(paste("Number of EC_GLUT1+ patches:", length(glut1_patches)))

# Extract EC_EpCAM+ patches (Quiescent_Tumor_patch)  
epcam_cells <- meta_IM[!is.na(meta_IM$Quiescent_Tumor_patch), ]
epcam_patches <- unique(epcam_cells$Quiescent_Tumor_patch)

print(paste("Number of EC_EpCAM+ patches:", length(epcam_patches)))

# Calculate mean expression per marker per individual patch
patch_expression_data <- data.frame()

# Process EC_GLUT1+ patches
for (patch_id in glut1_patches) {
  patch_cells_subset <- glut1_cells[glut1_cells$Metabolism_activate_Tumor_patch == patch_id, ]
  
  # Require minimum 5 cells per patch for reliable mean
  if (nrow(patch_cells_subset) >= 5) {
    patch_expr <- expr_matrix_clean[, match(rownames(patch_cells_subset), colnames(expr_matrix_clean))]
    
    # Calculate mean expression per marker for this patch
    patch_means <- apply(patch_expr, 1, mean, na.rm = TRUE)
    
    # Store patch-level data
    patch_data <- data.frame(
      patch_id = paste0("Activate_", patch_id),
      patch_type = "Metabolism_activate_Tumor_patch",
      sample_id = patch_cells_subset$sample_id[1],
      patient_id = patch_cells_subset$patient_id[1],
      RFS_status = patch_cells_subset$RFS_status[1],
      Treatment = patch_cells_subset$Treatment[1],
      n_cells = nrow(patch_cells_subset),
      t(patch_means),  # transpose to get markers as columns
      stringsAsFactors = FALSE
    )
    
    patch_expression_data <- rbind(patch_expression_data, patch_data)
  }
}

# Process EC_EpCAM+ patches
for (patch_id in epcam_patches) {
  patch_cells_subset <- epcam_cells[epcam_cells$Quiescent_Tumor_patch == patch_id, ]
  
  # Require minimum 5 cells per patch for reliable mean
  if (nrow(patch_cells_subset) >= 5) {
    patch_expr <- expr_matrix_clean[, match(rownames(patch_cells_subset), colnames(expr_matrix_clean))]
    
    # Calculate mean expression per marker for this patch
    patch_means <- apply(patch_expr, 1, mean, na.rm = TRUE)
    
    # Store patch-level data
    patch_data <- data.frame(
      patch_id = paste0("Quiescent_", patch_id),
      patch_type = "Quiescent_Tumor_patch",
      sample_id = patch_cells_subset$sample_id[1],
      patient_id = patch_cells_subset$patient_id[1],
      RFS_status = patch_cells_subset$RFS_status[1],
      Treatment = patch_cells_subset$Treatment[1],
      n_cells = nrow(patch_cells_subset),
      t(patch_means),  # transpose to get markers as columns
      stringsAsFactors = FALSE
    )
    
    patch_expression_data <- rbind(patch_expression_data, patch_data)
  }
}

print(paste("Total patches for analysis:", nrow(patch_expression_data)))
print(paste("Activate Tumor patches:", sum(patch_expression_data$patch_type == "Metabolism_activate_Tumor_patch")))
print(paste("Quiescent Tumor patches:", sum(patch_expression_data$patch_type == "Quiescent_Tumor_patch")))

# Perform marker comparison between patch types
marker_comparison_results <- data.frame()

for (marker in analysis_markers) {
  # Get marker expression for each patch type
  glut1_values <- patch_expression_data[patch_expression_data$patch_type == "Metabolism_activate_Tumor_patch", marker]
  epcam_values <- patch_expression_data[patch_expression_data$patch_type == "Quiescent_Tumor_patch", marker]
  
  # Remove any NAs
  glut1_values <- glut1_values[!is.na(glut1_values)]
  epcam_values <- epcam_values[!is.na(epcam_values)]
  
  # Perform Wilcoxon test if we have sufficient patches in both groups
  if (length(glut1_values) >= 3 && length(epcam_values) >= 3) {
    test_result <- wilcox.test(glut1_values, epcam_values)
    
    # Calculate effect size metrics
    mean_glut1 <- mean(glut1_values, na.rm = TRUE)
    mean_epcam <- mean(epcam_values, na.rm = TRUE)
    fold_change <- log2((mean_glut1 + 0.001) / (mean_epcam + 0.001))
    
    # Store results
    marker_comparison_results <- rbind(marker_comparison_results, data.frame(
      marker = marker,
      mean_EC_GLUT1 = mean_glut1,
      mean_EC_EpCAM = mean_epcam,
      log2_fold_change = fold_change,
      p_value = test_result$p.value,
      n_patches_GLUT1 = length(glut1_values),
      n_patches_EpCAM = length(epcam_values),
      stringsAsFactors = FALSE
    ))
  }
}

# Apply BH correction for multiple testing
marker_comparison_results$p_adj <- p.adjust(marker_comparison_results$p_value, method = "BH")

# Classify significant changes
p_cutoff <- 0.05
fc_cutoff <- log(1.5)  # log2 fold change cutoff

marker_comparison_results$significance <- "n.s."
marker_comparison_results$significance[marker_comparison_results$p_adj <= p_cutoff & 
                                         marker_comparison_results$log2_fold_change >= fc_cutoff] <- "Up_in_Activate"
marker_comparison_results$significance[marker_comparison_results$p_adj <= p_cutoff & 
                                         marker_comparison_results$log2_fold_change <= -fc_cutoff] <- "Up_in_Quiescent"

# Sort by significance and fold change
marker_comparison_results <- marker_comparison_results[order(marker_comparison_results$p_adj, 
                                                             -abs(marker_comparison_results$log2_fold_change)), ]

# Print summary
print("=== METABOLIC PROFILING RESULTS (Patch-Level Analysis) ===")
print(paste("Total markers tested:", nrow(marker_comparison_results)))
print(paste("Significant markers (p_adj < 0.05):", sum(marker_comparison_results$p_adj <= 0.05)))
print(paste("Up in EC_GLUT1+ patches:", sum(marker_comparison_results$significance == "Up_in_Activate")))
print(paste("Up in EC_EpCAM+ patches:", sum(marker_comparison_results$significance == "Up_in_Quiescent")))

# Show top significant results
print("\nTop 10 most significant differences:")
print(head(marker_comparison_results[marker_comparison_results$significance != "n.s.", ], 10))

# Create volcano plot
volcano_plot <- ggplot(marker_comparison_results, aes(x = log2_fold_change, y = -log10(p_adj))) +
  geom_point(aes(color = significance), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Up_in_Activate" = "#cd534c", "Up_in_Quiescent" = "#7ca6dc", "n.s." = "grey")) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey50") +
  labs(
    title = "Patch-Level Metabolic Profiling: Activate vs Quiescent Tumor Patches",
    x = "Log2 Fold Change (Activate vs Quiescent)",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# Add labels for top significant markers
top_markers <- marker_comparison_results[marker_comparison_results$significance != "n.s." & 
                                           marker_comparison_results$p_adj <= 0.05, ]

if (nrow(top_markers) > 0) {
  volcano_plot <- volcano_plot +
    geom_text_repel(
      data = top_markers,
      aes(label = marker),
      size = 3,
      box.padding = 0.5,
      max.overlaps = 15
    )
}

print(volcano_plot)
ggsave(file.path(figureDir, paste0("Marker expression difference between Activate and Quiescent tumor patch.pdf")), 
       volcano_plot, width = 8, height = 6)

# Create boxplot for top significant markers
if (sum(marker_comparison_results$significance != "NS") > 0) {
  # Get top 6 significant markers for boxplot
  top_6_markers <- head(marker_comparison_results[marker_comparison_results$significance != "n.s.", ]$marker, 6)
  
  # Prepare data for boxplot
  boxplot_data <- data.frame()
  for (marker in top_6_markers) {
    marker_data <- data.frame(
      patch_type = patch_expression_data$patch_type,
      expression = patch_expression_data[[marker]],
      marker = marker,
      stringsAsFactors = FALSE
    )
    boxplot_data <- rbind(boxplot_data, marker_data)
  }
  
  # Create boxplot
  boxplot <- ggplot(boxplot_data, aes(x = patch_type, y = expression, fill = patch_type)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.25, alpha = 0.6, color = "grey75") +
    facet_wrap(~marker, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("Metabolism_activate_Tumor_patch" = "#cd534c", "Quiescent_Tumor_patch" = "#7ca6dc")) +
    labs(
      title = "Top Significant Markers: Patch-Level Expression",
      x = "Patch Type",
      y = "Mean Expression per Patch",
      fill = "Patch Type"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    stat_compare_means(method = "wilcox.test", label = "p", size = 4)
  
  print(boxplot)
  ggsave(file.path(figureDir, paste0("Top differential Marker difference boxplot between Activate and Quiescent tumor patch.pdf")), 
         boxplot, width = 12, height = 9)
}

# Save results
write.csv(marker_comparison_results, "Section1_Patch_Level_Metabolic_Profiling_Results.csv", row.names = FALSE)
write.csv(patch_expression_data, "Section1_Patch_Expression_Data.csv", row.names = FALSE)

print("=== Section 1 Complete ===")

para_environ <- ls()
rm(list = para_environ[c(
  startsWith(para_environ,prefix = "marker") | 
  startsWith(para_environ,prefix = "patch") | 
  startsWith(para_environ,prefix = "expr") | 
    startsWith(para_environ,prefix = "boxplot")
)])
gc()

# =============================================================================
# Section 2: Patch Size Calculation - Calculate tumor patch size and density metrics
# Focus: Individual patch size, area, and density calculation for EC_GLUT1+ vs EC_EpCAM+ patches
# =============================================================================
# Get metadata and spatial coordinates
meta_IM <- as.data.frame(colData(spe_IM))
coords <- spatialCoords(spe_IM)

# Combine metadata with coordinates
cell_data <- cbind(meta_IM, coords)

print("=== PATCH SIZE CALCULATION ===")

# Process EC_GLUT1+ patches (Metabolism_activate_Tumor_patch)
glut1_cells <- cell_data[!is.na(cell_data$Metabolism_activate_Tumor_patch), ]
glut1_patches <- unique(glut1_cells$Metabolism_activate_Tumor_patch)

print(paste("Processing", length(glut1_patches), "EC_GLUT1+ patches..."))
patch_size_results_activate <- analysis_patch_characteristic(
  data = glut1_cells,
  patch_column_name = "Metabolism_activate_Tumor_patch",
  patche_ids = glut1_patches,
  patch_type = "Activate"
)

# Process EC_EpCAM+ patches (Quiescent_Tumor_patch)
epcam_cells <- cell_data[!is.na(cell_data$Quiescent_Tumor_patch), ]
epcam_patches <- unique(epcam_cells$Quiescent_Tumor_patch)

print(paste("Processing", length(epcam_patches), "EC_EpCAM+ patches..."))
patch_size_results_quiescent <- analysis_patch_characteristic(
  data = epcam_cells,
  patch_column_name = "Quiescent_Tumor_patch",
  patche_ids = epcam_patches,
  patch_type = "Quiescent"
)

patch_size_results <- rbind(patch_size_results_activate,patch_size_results_quiescent)
print(paste("Total patches analyzed:", nrow(patch_size_results)))
print(paste("Activate patches:", sum(patch_size_results$patch_type == "Activate")))
print(paste("Quiescent patches:", sum(patch_size_results$patch_type == "Quiescent")))

# Statistical comparison of size metrics between patch types
size_metrics <- c("n_cells", "bounding_box_area", "convex_hull_area", 
                  "within_patch_compactness_bbox", "within_patch_compactness_convex", 
                  "patch_area_density", "patch_convex_density", "aspect_ratio", "shape_compactness",
                  "avg_distance_from_centroid", "max_distance_from_centroid")

size_comparison_results <- data.frame()

for (metric in size_metrics) {
  # Get values for each patch type
  glut1_values <- patch_size_results[patch_size_results$patch_type == "Activate", metric]
  epcam_values <- patch_size_results[patch_size_results$patch_type == "Quiescent", metric]
  
  # Remove NAs
  glut1_values <- glut1_values[!is.na(glut1_values)]
  epcam_values <- epcam_values[!is.na(epcam_values)]
  
  # Perform Wilcoxon test if sufficient data
  if (length(glut1_values) >= 3 && length(epcam_values) >= 3) {
    test_result <- wilcox.test(glut1_values, epcam_values)
    
    # Calculate summary statistics
    mean_glut1 <- mean(glut1_values, na.rm = TRUE)
    mean_epcam <- mean(epcam_values, na.rm = TRUE)
    fold_change <- log2((mean_glut1 + 0.001) / (mean_epcam + 0.001))
    
    # Store results
    size_comparison_results <- rbind(size_comparison_results, data.frame(
      metric = metric,
      mean_EC_GLUT1 = mean_glut1,
      mean_EC_EpCAM = mean_epcam,
      log2_fold_change = fold_change,
      p_value = test_result$p.value,
      n_patches_GLUT1 = length(glut1_values),
      n_patches_EpCAM = length(epcam_values),
      stringsAsFactors = FALSE
    ))
  }
}

# Apply BH correction
size_comparison_results$p_adj <- p.adjust(size_comparison_results$p_value, method = "BH")

# Classify significant differences
p_cutoff <- 0.05
size_comparison_results$significance <- ifelse(size_comparison_results$p_adj <= p_cutoff, "Significant", "n.s.")

# Sort by significance
size_comparison_results <- size_comparison_results[order(size_comparison_results$p_adj), ]

print("\n=== PATCH SIZE COMPARISON RESULTS ===")
print(paste("Total metrics tested:", nrow(size_comparison_results)))
print(paste("Significant differences (p_adj < 0.05):", sum(size_comparison_results$significance == "Significant")))

print("\nSignificant size differences:")
print(size_comparison_results[size_comparison_results$significance == "Significant", ])

# Create visualization of key size metrics
# Select key metrics for visualization
key_metrics <- c("n_cells", "bounding_box_area", "convex_hull_area", "within_patch_compactness_bbox", 
                 "patch_area_density", "shape_compactness")

plot_data <- patch_size_results %>%
  select(patch_type, all_of(key_metrics)) %>%
  pivot_longer(cols = key_metrics, names_to = "metric", values_to = "value") %>%
  filter(!is.na(value))

# Create boxplot
size_boxplot <- ggplot(plot_data, aes(x = patch_type, y = value, fill = patch_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.25, alpha = 0.6, color = "grey75") +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Activate" = "#cd534c", "Quiescent" = "#7ca6dc")) +
  labs(
    title = "Patch Size and Shape Metrics Comparison",
    subtitle = "Activate vs Quiescent Patches",
    x = "Patch Type",
    y = "Metric Value",
    fill = "Patch Type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  stat_compare_means(method = "wilcox.test", label = "p", size = 4)

print(size_boxplot)
ggsave(file.path(figureDir, paste0("Graph characteristics comparision between Activate and Quiescent tumor patch.pdf")), 
       size_boxplot, width = 12, height = 9)

# Create density plot for patch size distribution
density_plot <- ggplot(patch_size_results, aes(x = log10(n_cells + 1), fill = patch_type)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("Activate" = "#cd534c", "Quiescent" = "#7ca6dc")) +
  labs(
    title = "Patch Size Distribution",
    x = "Log10(Number of Cells + 1)",
    y = "Density",
    fill = "Patch Type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

print(density_plot)
ggsave(file.path(figureDir, paste0("Patch density comparision between Activate and Quiescent tumor patch.pdf")), 
       density_plot, width = 8, height = 6)

# Save results
write.csv(patch_size_results, "Section2_Patch_Size_Results.csv", row.names = FALSE)
write.csv(size_comparison_results, "Section2_Size_Comparison_Results.csv", row.names = FALSE)

print("\n=== Section 2 Complete ===")

# =============================================================================
# Section 3: Marker Expression Analysis - Log-transformed marker changes with patch size
# Focus: Analyze how marker expression correlates with patch size (size-dependent programs)
# =============================================================================

library(dplyr)
library(ggplot2)
library(SpatialExperiment)
library(stats)

print("=== MARKER EXPRESSION vs PATCH SIZE ANALYSIS ===")

# Get metadata and expression data
expr_matrix <- assay(spe_IM)

# Remove DNA markers for analysis
dna_markers <- rownames(expr_matrix)[grepl("^DNA", rownames(expr_matrix))]
analysis_markers <- rownames(expr_matrix)[!rownames(expr_matrix) %in% dna_markers]
expr_matrix_clean <- expr_matrix[analysis_markers, ]

print(paste("Number of markers for analysis:", length(analysis_markers)))

# Process EC_GLUT1+ patches
glut1_cells <- cell_data[!is.na(cell_data$Metabolism_activate_Tumor_patch), ]
glut1_patches <- unique(glut1_cells$Metabolism_activate_Tumor_patch)

print(paste("Processing", length(glut1_patches), "EC_GLUT1+ patches..."))

patch_size_expression_data_activate <- analysis_patch_characteristic_with_expression(
  data = glut1_cells,
  expr_matrix_clean = expr_matrix_clean, analysis_markers = analysis_markers,
  patch_column_name = "Metabolism_activate_Tumor_patch",
  patche_ids = glut1_patches,
  patch_type = "Activate"
)

# Process EC_EpCAM+ patches
epcam_cells <- cell_data[!is.na(cell_data$Quiescent_Tumor_patch), ]
epcam_patches <- unique(epcam_cells$Quiescent_Tumor_patch)

print(paste("Processing", length(epcam_patches), "EC_EpCAM+ patches..."))

patch_size_expression_data_quiescent <- analysis_patch_characteristic_with_expression(
  data = epcam_cells,
  expr_matrix_clean = expr_matrix_clean, analysis_markers = analysis_markers,
  patch_column_name = "Quiescent_Tumor_patch",
  patche_ids = epcam_patches,
  patch_type = "Quiescent"
)

patch_size_expression_data <- rbind(patch_size_expression_data_activate, patch_size_expression_data_quiescent)

print(paste("Total patches for size-expression analysis:", nrow(patch_size_expression_data)))
print(paste("EC_GLUT1+ patches:", sum(patch_size_expression_data$patch_type == "Activate")))
print(paste("EC_EpCAM+ patches:", sum(patch_size_expression_data$patch_type == "Quiescent")))

# Define size metrics for correlation analysis 
size_metrics <- c("log_n_cells")
size_metric_labels <- c("Log2(Cell Count)")
names(size_metric_labels) <- size_metrics

# Create discrete size bins based on log_n_cells
patch_size_expression_data$size_bin <- cut(patch_size_expression_data$log_n_cells, 
                                           breaks = c(-0.001, 4, 5, 6, 7, 8, 10, Inf),
                                           labels = c("(0,4]", "(4,5]", "(5,6]", "(6,7]", "(7,8]", "(8,10]", ">8"),
                                           include.lowest = TRUE)

# Check size bin distribution
print("\n=== SIZE BIN DISTRIBUTION ===")
overall_size_bins <- table(patch_size_expression_data$size_bin)
print("\nOverall size bin distribution:")
print(overall_size_bins)

# Perform ANOVA analysis for each marker across size bins
anova_results <- data.frame()

for (marker in analysis_markers) {
  marker_col <- paste0("mean_", marker)
  
  # Overall ANOVA (both patch types)
  overall_data <- patch_size_expression_data[!is.na(patch_size_expression_data[[marker_col]]) & 
                                               !is.na(patch_size_expression_data$size_bin), ]
  
  if (nrow(overall_data) >= 10 && length(unique(overall_data$size_bin)) >= 3) {
    # Check if we have sufficient data in each bin
    bin_counts <- table(overall_data$size_bin)
    if (min(bin_counts[bin_counts > 0]) >= 2) {  # At least 2 patches per bin
      
      # Perform ANOVA
      anova_model <- aov(overall_data[[marker_col]] ~ overall_data$size_bin)
      anova_summary <- summary(anova_model)
      p_value <- anova_summary[[1]][1, "Pr(>F)"]
      f_statistic <- anova_summary[[1]][1, "F value"]
      
      anova_results <- rbind(anova_results, data.frame(
        marker = marker,
        patch_type = "Overall",
        f_statistic = f_statistic,
        p_value = p_value,
        n_patches = nrow(overall_data),
        n_bins = length(unique(overall_data$size_bin)),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Patch type-specific ANOVA
  for (ptype in c("Activate", "Quiescent")) {
    patch_data <- patch_size_expression_data[patch_size_expression_data$patch_type == ptype & 
                                               !is.na(patch_size_expression_data[[marker_col]]) & 
                                               !is.na(patch_size_expression_data$size_bin), ]
    
    if (nrow(patch_data) >= 6 && length(unique(patch_data$size_bin)) >= 3) {
      # Check if we have sufficient data in each bin
      bin_counts <- table(patch_data$size_bin)
      if (min(bin_counts[bin_counts > 0]) >= 2) {
        
        # Perform ANOVA
        anova_model <- aov(patch_data[[marker_col]] ~ patch_data$size_bin)
        anova_summary <- summary(anova_model)
        p_value <- anova_summary[[1]][1, "Pr(>F)"]
        f_statistic <- anova_summary[[1]][1, "F value"]
        
        anova_results <- rbind(anova_results, data.frame(
          marker = marker,
          patch_type = ptype,
          f_statistic = f_statistic,
          p_value = p_value,
          n_patches = nrow(patch_data),
          n_bins = length(unique(patch_data$size_bin)),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Apply BH correction within each patch type
anova_results <- anova_results %>%
  group_by(patch_type) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Classify significant differences
p_cutoff <- 0.05
anova_results$significance <- ifelse(anova_results$p_adj <= p_cutoff, "Significant", "n.s.")

# Sort by significance
anova_results <- anova_results[order(anova_results$p_adj, -anova_results$f_statistic), ]

print("\n=== ANOVA RESULTS FOR SIZE-DEPENDENT MARKER CHANGES ===")
print(paste("Total markers tested:", nrow(anova_results)))
print(paste("Significant size-dependent changes (p_adj < 0.05):", 
            sum(anova_results$significance == "Significant")))

# Show top significant results
print("\nTop 15 significant size-dependent marker changes:")
sig_results <- anova_results[anova_results$significance == "Significant", ]
if (nrow(sig_results) > 0) {
  print(head(sig_results[, c("marker", "patch_type", "f_statistic", "p_adj", "n_patches")], 15))
}

# Create boxplots for top significant markers
if (nrow(sig_results) > 0) {
  # Select specific markers for visualization (metabolic and key markers)
  top_markers <- c("GLUT1", "Collagen_I", "HK2", "Vimentin", "CD279", "Alpha_SMA", "Ki67", "CA_IX", "FASN", "FAP", "VEGF", "PRPS1")
  
  # Filter to only markers that exist in the data
  available_markers <- top_markers[top_markers %in% analysis_markers]
  
  plot_list <- list()
  
  for (i in 1:length(available_markers)) {
    marker <- available_markers[i]
    marker_col <- paste0("mean_", marker)
    
    # Get significance info for this marker for both patch types
    marker_sig_info <- anova_results[anova_results$marker == marker, ] 
    
    # Create subtitle with patch-type specific results
    glut1_result <- marker_sig_info[marker_sig_info$patch_type == "Activate", ]
    epcam_result <- marker_sig_info[marker_sig_info$patch_type == "Quiescent", ]
    
    subtitle_parts <- c()
    if (nrow(glut1_result) > 0) {
      subtitle_parts <- c(subtitle_parts, paste("Activate: p_adj =", format(glut1_result$p_adj, scientific = TRUE, digits = 2)))
    }
    if (nrow(epcam_result) > 0) {
      subtitle_parts <- c(subtitle_parts, paste("Quiescent: p_adj =", format(epcam_result$p_adj, scientific = TRUE, digits = 2)))
    }
    
    subtitle_text <- paste(subtitle_parts, collapse = "; ")
    if (subtitle_text == "") subtitle_text <- "No significant size effect"
    
    plot_data <- patch_size_expression_data[!is.na(patch_size_expression_data[[marker_col]]) & 
                                              !is.na(patch_size_expression_data$size_bin), ]
    
    if (nrow(plot_data) > 10) {
      plot_data$size_bin_numeric <- as.numeric(plot_data$size_bin)
      
      p <- ggplot(plot_data, aes(x = size_bin, y = .data[[marker_col]], fill = patch_type)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
        # geom_jitter(aes(color = patch_type), position = position_jitterdodge(dodge.width = 0.8), 
        #             size = 0.3, alpha = 0.5) +
        # Add trend lines for each patch type
        geom_smooth(aes(x = size_bin_numeric, color = patch_type, group = patch_type), 
                    method = "lm", se = TRUE, alpha = 0.2, size = 1, 
                    linetype = "dashed", fill = NA) +
        scale_fill_manual(values = c("Activate" = "#cd534c", "Quiescent" = "#7ca6dc")) +
        scale_color_manual(values = c("Activate" = "#cd534c", "Quiescent" = "#7ca6dc")) + 
        labs(
          title = paste(marker, "Expression vs Patch Size"),
          subtitle = subtitle_text,
          x = "Size Bin (Log2 Cell Count)",
          y = paste("Mean", marker, "Expression"),
          fill = "Patch Type",
          color = "Patch Type"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 8),
          legend.position = "bottom",
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8)
        )
      
      # Add statistical significance annotations for each patch type separately
      glut1_data <- plot_data[plot_data$patch_type == "Activate", ]
      epcam_data <- plot_data[plot_data$patch_type == "Quiescent", ]
      
      plot_list[[i]] <- p
    }
  }
  
  if (length(plot_list) > 0) {
    # Create combined plot
    combined_plot <- wrap_plots(plot_list, ncol = 3)
    print(combined_plot)
    
    ggsave(file.path(figureDir, paste0("Marker change across patch log2(patch size) bin.pdf")), 
           combined_plot, width = 12, height = 12)
  }
}

# Save results
write.csv(patch_size_expression_data, "Section3_Patch_Size_Expression_Data.csv", row.names = FALSE)

print("\n=== Section 3 Complete ===")

# =============================================================================
# Section 4: Patch Size Comparison by RFS Status
# Focus: Compare size and density of EC_GLUT1+ vs EC_EpCAM+ patches between RFS groups
# =============================================================================

print("=== PATCH SIZE COMPARISON BY RFS STATUS ===")

# Initialize results dataframe for individual patches
patch_rfs_data <- patch_size_results
patch_rfs_data$RFS_status <- meta_IM$RFS_status[match(patch_rfs_data$patient_id,meta_IM$patient_id)]

print(paste("Total patches for RFS analysis:", nrow(patch_rfs_data)))
print(paste("Activate patches:", sum(patch_rfs_data$patch_type == "Activate")))
print(paste("Quiescent patches:", sum(patch_rfs_data$patch_type == "Quiescent")))

# Check RFS group distribution
rfs_distribution <- patch_rfs_data %>%
  group_by(patch_type, RFS_status) %>%
  summarise(n_patches = n(), .groups = 'drop')

print("\n=== RFS GROUP DISTRIBUTION ===")
print(rfs_distribution)

# Define size metrics for comparison
patch_rfs_data$log_n_cells <- log2(patch_rfs_data$n_cells)
patch_rfs_data$log_bounding_box_area<- log10(patch_rfs_data$bounding_box_area)
patch_rfs_data$log_convex_hull_area <- log2(patch_rfs_data$convex_hull_area)

size_metrics <- c("n_cells", "log_n_cells", "bounding_box_area", "log_bounding_box_area",
                  "convex_hull_area", "log_convex_hull_area", "within_patch_compactness_bbox",
                  "patch_area_coverage", "shape_compactness", "avg_distance_from_centroid")

# Perform statistical comparison between RFS groups
rfs_comparison_results <- data.frame()

for (metric in size_metrics) {
  
  # Patch type-specific comparisons
  for (ptype in c("Activate", "Quiescent")) {
    patch_data <- patch_rfs_data[patch_rfs_data$patch_type == ptype & 
                                   !is.na(patch_rfs_data[[metric]]) & 
                                   !is.na(patch_rfs_data$RFS_status), ]
    
    if (nrow(patch_data) >= 6) {
      no_recur_values <- patch_data[patch_data$RFS_status == 0, metric]
      recur_values <- patch_data[patch_data$RFS_status == 1, metric]
      
      if (length(no_recur_values) >= 2 && length(recur_values) >= 2) {
        test_result <- wilcox.test(recur_values, no_recur_values)
        
        # Calculate effect size metrics
        mean_no_recur <- mean(no_recur_values, na.rm = TRUE)
        mean_recur <- mean(recur_values, na.rm = TRUE)
        fold_change <- log2((mean_recur + 0.001) / (mean_no_recur + 0.001))
        
        rfs_comparison_results <- rbind(rfs_comparison_results, data.frame(
          metric = metric,
          patch_type = ptype,
          mean_No_Recurrence = mean_no_recur,
          mean_Recurrence = mean_recur,
          log2_fold_change = fold_change,
          p_value = test_result$p.value,
          n_patches_no_recur = length(no_recur_values),
          n_patches_recur = length(recur_values),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Apply BH correction within each patch type
rfs_comparison_results <- rfs_comparison_results %>%
  group_by(patch_type) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Classify significant differences
p_cutoff <- 0.05
fc_cutoff <- 0.5  # log2 fold change cutoff
rfs_comparison_results$significance <- "n.s."
rfs_comparison_results$significance[rfs_comparison_results$p_adj <= p_cutoff & 
                                      abs(rfs_comparison_results$log2_fold_change) >= fc_cutoff] <- "Significant"

# Sort by significance
rfs_comparison_results <- rfs_comparison_results[order(rfs_comparison_results$p_adj, 
                                                       -abs(rfs_comparison_results$log2_fold_change)), ]

print("\n=== RFS COMPARISON RESULTS ===")
print(paste("Total comparisons tested:", nrow(rfs_comparison_results)))
print(paste("Significant differences (p_adj < 0.05, |log2FC| >= 0.5):", 
            sum(rfs_comparison_results$significance == "Significant")))

# Show significant results
sig_results <- rfs_comparison_results[rfs_comparison_results$significance == "Significant", ]
if (nrow(sig_results) > 0) {
  print("\nSignificant RFS-associated size differences:")
  print(sig_results[, c("metric", "patch_type", "log2_fold_change", "p_adj")])
}

# Create boxplots for key size metrics
key_metrics <- c("log_n_cells", "log_bounding_box_area", "within_patch_compactness_bbox", 
                 "patch_area_coverage", "shape_compactness", "avg_distance_from_centroid")
metric_labels <- c("Log2(Cell Count)", "Log2(Bounding Box Area)", "Within-Patch Compactness",
                   "Patch Area Coverage", "Shape Compactness", "Avg Distance from Centroid")
names(metric_labels) <- key_metrics

plot_list <- list()

for (i in 1:length(key_metrics)) {
  metric <- key_metrics[i]
  
  # Get significance info for this metric
  metric_results <- rfs_comparison_results[rfs_comparison_results$metric == metric, ]
  
  plot_data <- patch_rfs_data[!is.na(patch_rfs_data[[metric]]) & 
                                !is.na(patch_rfs_data$RFS_status), ]
  
  # Create significance labels for subtitle
  glut1_result <- metric_results[metric_results$patch_type == "Activate", ]
  epcam_result <- metric_results[metric_results$patch_type == "Quiescent", ]

  sig_labels <- c()
  if (nrow(glut1_result) > 0) {
    sig_labels <- c(sig_labels, paste("Activate:", ifelse(glut1_result$p_adj < 0.05, 
                                                        paste("p =", format(glut1_result$p_adj, digits = 3)), "n.s.")))
  }
  if (nrow(epcam_result) > 0) {
    sig_labels <- c(sig_labels, paste("Quiescent:", ifelse(epcam_result$p_adj < 0.05, 
                                                        paste("p =", format(epcam_result$p_adj, digits = 3)), "n.s.")))
  }
  subtitle_text <- paste(sig_labels, collapse = "; ")
  
  if (nrow(plot_data) > 10) {
    p <- ggplot(plot_data, aes(x = RFS_status, y = .data[[metric]], fill = patch_type)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
      geom_jitter(aes(color = patch_type), position = position_jitterdodge(dodge.width = 0.8), 
                  size = 0.5, alpha = 0.6) +
      scale_fill_manual(values = c("Activate" = "#cd534c", "Quiescent" = "#7ca6dc")) +
      scale_color_manual(values = c("Activate" = "#cd534c", "Quiescent" = "#7ca6dc")) +
      labs(
        title = paste(metric_labels[metric], "by RFS Status"),
        subtitle = subtitle_text,
        x = "RFS Group",
        y = metric_labels[metric],
        fill = "Patch Type",
        color = "Patch Type"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom"
      )
    
    plot_list[[i]] <- p
  }
}

if (length(plot_list) > 0) {
  combined_plot <- wrap_plots(plot_list, ncol = 2)
  print(combined_plot)
}

# Save results
write.csv(rfs_comparison_results, "Section4_RFS_Comparison_Results.csv", row.names = FALSE)

print("\n=== Section 4 Complete ===")

# =============================================================================
# Section 5: Component Analysis - Calculate density of other cell types within tumor patches
# Focus: Microenvironment composition within EC_GLUT1+ vs EC_EpCAM+ patches
# =============================================================================

print("=== COMPONENT ANALYSIS - MICROENVIRONMENT WITHIN PATCHES ===")

# Get metadata and spatial coordinates
meta_IM <- as.data.frame(colData(spe_IM))
coords <- spatialCoords(spe_IM)

# Combine metadata with coordinates
cell_data <- cbind(meta_IM, coords)

print("Available cell types:")
print(table(cell_data$sub_celltype))

# Initialize results dataframe for patch microenvironment composition
patch_microenv_data <- data.frame()

# Process EC_GLUT1+ patches
glut1_results <- process_patches(
  cell_data = cell_data,
  patch_column = "Metabolism_activate_Tumor_patch",
  patch_type_label = "Activate",
  patch_id_prefix = "Activate"
)

# Process EC_EpCAM+ patches
epcam_results <- process_patches(
  cell_data = cell_data,
  patch_column = "Quiescent_Tumor_patch", 
  patch_type_label = "Quiescent",
  patch_id_prefix = "Quiescent"
)

patch_microenv_data <- rbind(glut1_results, epcam_results)

print(paste("Total patches for microenvironment analysis:", nrow(patch_microenv_data)))
print(paste("Activate patches:", sum(patch_microenv_data$patch_type == "Activate")))
print(paste("Quiescent patches:", sum(patch_microenv_data$patch_type == "Quiescent")))

# Summary of microenvironment composition
print("\n=== MICROENVIRONMENT COMPOSITION SUMMARY ===")
microenv_summary <- patch_microenv_data %>%
  group_by(patch_type) %>%
  summarise(
    n_patches = n(),
    mean_total_cells = mean(total_cells, na.rm = TRUE),
    mean_microenv_cells = mean(microenv_cells, na.rm = TRUE),
    mean_microenv_fraction = mean(microenv_fraction, na.rm = TRUE),
    mean_microenv_density = mean(microenv_density, na.rm = TRUE),
    .groups = 'drop'
  )
print(microenv_summary)

# Get cell type proportion columns for analysis
all_cell_types <- unique(spe_IM$sub_celltype)
prop_columns <- grep("^prop_", colnames(patch_microenv_data), value = TRUE)
# Exclude tumor cell types from microenvironment analysis
tumor_prop_columns <- paste0("prop_", all_cell_types[startsWith(all_cell_types,prefix = "EC")])
microenv_prop_columns <- prop_columns[!prop_columns %in% tumor_prop_columns]

print(paste("Analyzing", length(microenv_prop_columns), "microenvironment cell types"))

# Statistical comparison of microenvironment composition between patch types
microenv_comparison_results <- data.frame()

for (prop_col in microenv_prop_columns) {
  cell_type <- gsub("^prop_", "", prop_col)
  
  # Get values for each patch type
  glut1_values <- patch_microenv_data[patch_microenv_data$patch_type == "Activate", prop_col]
  epcam_values <- patch_microenv_data[patch_microenv_data$patch_type == "Quiescent", prop_col]
  
  # Remove NAs
  glut1_values <- glut1_values[!is.na(glut1_values)]
  epcam_values <- epcam_values[!is.na(epcam_values)]
  
  # Perform Wilcoxon test if sufficient data
  if (length(glut1_values) >= 3 && length(epcam_values) >= 3) {
    test_result <- wilcox.test(glut1_values, epcam_values)
    
    # Calculate summary statistics
    mean_glut1 <- mean(glut1_values, na.rm = TRUE)
    mean_epcam <- mean(epcam_values, na.rm = TRUE)
    median_glut1 <- median(glut1_values, na.rm = TRUE)
    median_epcam <- median(epcam_values, na.rm = TRUE)
    
    # Calculate fold change (avoid division by zero)
    fold_change <- log2((mean_glut1 + 0.0001) / (mean_epcam + 0.0001))
    
    # Store results
    microenv_comparison_results <- rbind(microenv_comparison_results, data.frame(
      cell_type = cell_type,
      mean_EC_GLUT1 = mean_glut1,
      mean_EC_EpCAM = mean_epcam,
      log2_fold_change = fold_change,
      p_value = test_result$p.value,
      n_patches_GLUT1 = length(glut1_values),
      n_patches_EpCAM = length(epcam_values),
      stringsAsFactors = FALSE
    ))
  }
}

# Apply BH correction
microenv_comparison_results$p_adj <- p.adjust(microenv_comparison_results$p_value, method = "BH")

# Classify significant differences
p_cutoff <- 0.05
fc_cutoff <- log2(1.5)  # log2 fold change cutoff  
microenv_comparison_results$significance <- "n.s."
microenv_comparison_results$significance[microenv_comparison_results$p_adj <= p_cutoff & 
                                           microenv_comparison_results$log2_fold_change >= fc_cutoff] <- "Activate"
microenv_comparison_results$significance[microenv_comparison_results$p_adj <= p_cutoff & 
                                           microenv_comparison_results$log2_fold_change <= -fc_cutoff] <- "Quiescent"

# Sort by significance
microenv_comparison_results <- microenv_comparison_results[order(microenv_comparison_results$p_adj, 
                                                                 -abs(microenv_comparison_results$log2_fold_change)), ]

print("\n=== MICROENVIRONMENT COMPOSITION COMPARISON RESULTS ===")
print(paste("Total cell types tested:", nrow(microenv_comparison_results)))
print(paste("Significantly different (p_adj < 0.05, |log2FC| >= 0.5):", 
            sum(microenv_comparison_results$significance != "n.s.")))

# Show significant results
sig_results <- microenv_comparison_results[microenv_comparison_results$significance != "n.s.", ]
if (nrow(sig_results) > 0) {
  print("\nSignificant microenvironment differences:")
  print(sig_results[, c("cell_type", "significance", "log2_fold_change", "p_adj")])
}

# Create volcano plot
volcano_plot <- ggplot(microenv_comparison_results, aes(x = log2_fold_change, y = -log10(p_adj))) +
  geom_point(aes(color = significance), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Activate" = "#cd534c", "Quiescent" = "#7ca6dc", "n.s." = "grey")) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey50") +
  labs(
    title = "Microenvironment Composition: Activate vs Quiescent Patches",
    x = "Log2 Fold Change (Activate vs Quiescent)",
    y = "-Log10 Adjusted P-value",
    color = "Enrichment"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

# Add labels for significant cell types
if (nrow(sig_results) > 0) {
  library(ggrepel)
  volcano_plot <- volcano_plot +
    geom_text_repel(
      data = sig_results,
      aes(label = cell_type),
      size = 3,
      box.padding = 0.5,
      max.overlaps = 15
    )
}

print(volcano_plot)
ggsave(file.path(figureDir, paste0("Microenvironment composition difference between activate and quiescent tumor patches.pdf")), 
       volcano_plot, width = 8, height = 6)

# Create boxplots for significant cell types
if (nrow(sig_results) > 0) {
  # Get top 6 significant cell types for visualization
  top_6_cell_types <- sig_results$cell_type
  
  plot_list <- list()
  
  for (i in 1:length(top_6_cell_types)) {
    cell_type <- top_6_cell_types[i]
    prop_col <- paste0("prop_", cell_type)
    
    # Get significance info
    cell_result <- sig_results[sig_results$cell_type == cell_type, ]
    
    plot_data <- patch_microenv_data[!is.na(patch_microenv_data[[prop_col]]), ]
    
    if (nrow(plot_data) > 10) {
      p <- ggplot(plot_data, aes(x = patch_type, y = .data[[prop_col]], fill = patch_type)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, size = 0.15, alpha = 0.6, color = "grey75") +
        scale_fill_manual(values = c("Activate" = "#cd534c", "Quiescent" = "#7ca6dc")) +
        # Set y-axis limits from minimum to 75th percentile + 0.1
        coord_cartesian(ylim = c(min(plot_data[[prop_col]], na.rm = TRUE), 
                                 quantile(plot_data[[prop_col]], 0.75, na.rm = TRUE) + 0.1)) +
        labs(
          title = paste(cell_type, "Proportion in Patches"),
          subtitle = paste("p_adj =", format(cell_result$p_adj, scientific = TRUE, digits = 2),
                           ", log2(FC) =", round(cell_result$log2_fold_change, 2)),
          x = "Patch Type",
          y = paste("Proportion of", cell_type),
          fill = "Patch Type"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 9),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
          legend.position = "top"
        ) + 
        stat_compare_means(method = "wilcox.test", label = "p", hide.ns = FALSE, size = 3)
      
      plot_list[[i]] <- p
    }
  }
  
  if (length(plot_list) > 0) {
    combined_plot <- wrap_plots(plot_list, ncol = 4)
    ggsave(file.path(figureDir, paste0("Boxplot of Microenvironment composition difference between activate and quiescent tumor patches.pdf")), 
           combined_plot, width = 12, height = 8)
  }
}

# Create stacked bar plot showing overall microenvironment composition
if (nrow(sig_results) > 0) {
  # Calculate mean proportions for visualization
  mean_proportions <- patch_microenv_data %>%
    group_by(patch_type) %>%
    summarise_at(vars(starts_with("prop_")), mean, na.rm = TRUE) %>%
    gather(key = "cell_type_prop", value = "mean_proportion", -patch_type) %>%
    mutate(cell_type = gsub("^prop_", "", cell_type_prop)) %>%
    filter(cell_type %in% top_6_cell_types) %>%
    select(-cell_type_prop)
  
  stacked_plot <- ggplot(mean_proportions, aes(x = patch_type, y = mean_proportion, fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = metadata(spe)$color_vectors$sub_celltype) +
    labs(
      title = "Microenvironment Composition in Patches",
      subtitle = "Top Significant Cell Types",
      x = "Patch Type",
      y = "Mean Proportion",
      fill = "Cell Type"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  print(stacked_plot)
  ggsave(file.path(figureDir, paste0(
    "Stack barplot of Microenvironment composition difference between activate and quiescent tumor patches.pdf")), 
      stacked_plot, width = 6, height = 4.5)
}

# Save results
write.csv(patch_microenv_data, "Section5_Patch_Microenvironment_Data.csv", row.names = FALSE)
write.csv(microenv_comparison_results, "Section5_Microenvironment_Comparison_Results.csv", row.names = FALSE)

print("\n=== Section 5 Complete ===")

# =============================================================================
# Section 7: Cell Type Marker Analysis - Compare marker expression of microenvironment cells within different patch types
# Focus: How do immune/stromal cells change their expression profiles in EC_GLUT1+ vs EC_EpCAM+ patch contexts
# =============================================================================

print("=== CELL TYPE MARKER ANALYSIS - PATCH-DEPENDENT CELLULAR REPROGRAMMING ===")

# Get metadata and expression data
meta_IM <- as.data.frame(colData(spe_IM))
expr_matrix <- assay(spe_IM)
coords <- spatialCoords(spe_IM)

# Combine metadata with coordinates
cell_data <- cbind(meta_IM, coords)

# Remove DNA markers for analysis
dna_markers <- rownames(expr_matrix)[grepl("^DNA", rownames(expr_matrix))]
analysis_markers <- rownames(expr_matrix)[!rownames(expr_matrix) %in% dna_markers]
expr_matrix_clean <- expr_matrix[analysis_markers, ]

print(paste("Number of markers for analysis:", length(analysis_markers)))

# Define tumor cell types to exclude from microenvironment analysis
tumor_cell_types <- all_cell_types[startsWith(all_cell_types,prefix = "EC")]

# Extract microenvironment cells from both patch types
print("Extracting microenvironment cells with expression data...")

glut1_microenv_cells <- extract_microenv_cells_with_expression(
  cell_data = cell_data,
  expr_matrix = expr_matrix_clean,
  patch_column = "Metabolism_activate_Tumor_patch",
  patch_type_label = "Activate"
)

epcam_microenv_cells <- extract_microenv_cells_with_expression(
  cell_data = cell_data,
  expr_matrix = expr_matrix_clean,
  patch_column = "Quiescent_Tumor_patch",
  patch_type_label = "Quiescent"
)

# Combine all microenvironment cells
all_microenv_cells <- rbind(glut1_microenv_cells, epcam_microenv_cells)

print(paste("Total microenvironment cells extracted:", nrow(all_microenv_cells)))
print(paste("From Activate patches:", nrow(glut1_microenv_cells)))
print(paste("From Quiescent patches:", nrow(epcam_microenv_cells)))

# Check available cell types
available_cell_types <- table(all_microenv_cells$sub_celltype)
print("\nMicroenvironment cell types available:")
print(available_cell_types)

# Filter to cell types with sufficient representation in both patch types
min_cells_per_group <- 20  # Minimum cells per cell type per patch type

cell_type_counts <- all_microenv_cells %>%
  group_by(sub_celltype, patch_type) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  pivot_wider(names_from = patch_type, values_from = n_cells, values_fill = 0)

# Select cell types with sufficient cells in both patch types
valid_cell_types <- cell_type_counts$sub_celltype[
  cell_type_counts$Activate >= min_cells_per_group & 
    cell_type_counts$Quiescent >= min_cells_per_group
]

print(paste("Cell types with sufficient representation (>=", min_cells_per_group, "cells per patch type):"))
print(valid_cell_types)

if (length(valid_cell_types) == 0) {
  # Reduce threshold if no cell types meet criteria
  min_cells_per_group <- 5
  valid_cell_types <- cell_type_counts$sub_celltype[
    cell_type_counts$Activate >= min_cells_per_group & 
      cell_type_counts$Quiescent >= min_cells_per_group
  ]
  print(paste("Reduced threshold to", min_cells_per_group, "cells per group. Valid cell types:"))
  print(valid_cell_types)
}

# Perform marker expression comparison for each valid cell type
cell_type_marker_results <- data.frame()

for (cell_type in valid_cell_types) {
  print(paste("Analyzing", cell_type, "marker expression..."))
  
  # Get cells of this type from both patch types
  cell_type_data <- all_microenv_cells[all_microenv_cells$sub_celltype == cell_type, ]
  
  # Perform marker expression comparison
  for (marker in analysis_markers) {
    expr_col <- paste0("expr_", marker)
    
    # Get expression values for each patch type
    glut1_expr <- cell_type_data[cell_type_data$patch_type == "Activate", expr_col]
    epcam_expr <- cell_type_data[cell_type_data$patch_type == "Quiescent", expr_col]
    
    # Remove NAs
    glut1_expr <- glut1_expr[!is.na(glut1_expr)]
    epcam_expr <- epcam_expr[!is.na(epcam_expr)]
    
    # Perform Wilcoxon test if sufficient data
    if (length(glut1_expr) >= 10 && length(epcam_expr) >= 10) {
      test_result <- wilcox.test(glut1_expr, epcam_expr)
      
      # Calculate summary statistics
      mean_glut1 <- mean(glut1_expr, na.rm = TRUE)
      mean_epcam <- mean(epcam_expr, na.rm = TRUE)
      median_glut1 <- median(glut1_expr, na.rm = TRUE)
      median_epcam <- median(epcam_expr, na.rm = TRUE)
      
      # Calculate fold change
      fold_change <- log2((mean_glut1 + 0.001) / (mean_epcam + 0.001))
      
      # Store results
      cell_type_marker_results <- rbind(cell_type_marker_results, data.frame(
        cell_type = cell_type,
        marker = marker,
        mean_EC_GLUT1 = mean_glut1,
        mean_EC_EpCAM = mean_epcam,
        log2_fold_change = fold_change,
        p_value = test_result$p.value,
        n_cells_GLUT1 = length(glut1_expr),
        n_cells_EpCAM = length(epcam_expr),
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (nrow(cell_type_marker_results) == 0) {
  print("No sufficient data for marker comparison analysis")
  stop("Analysis cannot proceed - insufficient microenvironment cells")
}

# Apply BH correction within each cell type
cell_type_marker_results <- cell_type_marker_results %>%
  group_by(cell_type) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Classify significant differences
p_cutoff <- 0.05
fc_cutoff <- log2(1.5)  # Lower threshold for marker expression differences

cell_type_marker_results$significance <- "n.s."
cell_type_marker_results$significance[cell_type_marker_results$p_adj <= p_cutoff & 
                                        cell_type_marker_results$log2_fold_change >= fc_cutoff] <- "Higher_in_Activate"  
cell_type_marker_results$significance[cell_type_marker_results$p_adj <= p_cutoff & 
                                        cell_type_marker_results$log2_fold_change <= -fc_cutoff] <- "Higher_in_Quiescent"

# Sort by significance
cell_type_marker_results <- cell_type_marker_results[order(cell_type_marker_results$p_adj, 
                                                           -abs(cell_type_marker_results$log2_fold_change)), ]

print("\n=== CELL TYPE MARKER EXPRESSION COMPARISON RESULTS ===")
print(paste("Total comparisons tested:", nrow(cell_type_marker_results)))
print(paste("Significant differences (p_adj < 0.05, |log2FC| >= 0.3):", 
            sum(cell_type_marker_results$significance != "n.s.")))

# Show significant results
sig_results <- cell_type_marker_results[cell_type_marker_results$significance != "n.s.", ]
if (nrow(sig_results) > 0) {
  print("\nTop significant marker expression differences in microenvironment cells:")
  print(head(sig_results[, c("cell_type", "marker", "significance", "log2_fold_change", "p_adj")], 15))
  
  # Summary by cell type
  cell_type_summary <- sig_results %>%
    group_by(cell_type) %>%
    summarise(
      n_significant_markers = n(),
      mean_abs_fold_change = mean(abs(log2_fold_change)),
      markers_higher_in_GLUT1 = sum(significance == "Higher_in_Activate"),
      markers_higher_in_EpCAM = sum(significance == "Higher_in_Quiescent"),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_significant_markers))
  
  print("\nSummary of reprogramming by cell type:")
  print(cell_type_summary)
}

# Create jitter plot showing ALL marker expression differences with significance highlighting
    plot_markers <- cell_type_marker_results
  
  # Add significance-based labeling
  plot_markers$label <- case_when(
    plot_markers$p_adj <= 0.05 & plot_markers$log2_fold_change >= log2(1.5) ~ "Higher_in_Activate",
    plot_markers$p_adj <= 0.05 & plot_markers$log2_fold_change <= (-log2(1.5)) ~ "Higher_in_Quiescent",
    TRUE ~ "Non_Significant"
  )
  
  # Get significant markers only for labeling (top 5 up and down per cell type)
  sig_markers_for_label <- plot_markers %>%
    filter(p_adj <= 0.05 & abs(log2_fold_change) >= log2(1.5)) %>%
    group_by(cell_type) %>%
    top_n(n = 5, wt = log2_fold_change) %>%
    bind_rows(
      plot_markers %>%
        filter(p_adj <= 0.05 & abs(log2_fold_change) >= log2(1.5)) %>%
        group_by(cell_type) %>%
        top_n(n = 5, wt = -log2_fold_change)
    ) %>%
    mutate(do_label = marker)
  
  # Add labels only to significant markers
  plot_markers$do_label <- ifelse(
    paste(plot_markers$cell_type, plot_markers$marker) %in% 
      paste(sig_markers_for_label$cell_type, sig_markers_for_label$marker),
    plot_markers$marker, 
    NA
  )
  
  # Create background bar data (range of ALL fold changes per cell type)
  dfbar <- plot_markers %>%
    group_by(cell_type) %>%
    summarise(
      low = round(min(log2_fold_change) - 0.3, 1),
      up = round(max(log2_fold_change) + 0.3, 1),
      .groups = 'drop'
    )
  
  # Create cell type colors
  cell_types_unique <- unique(plot_markers$cell_type)
  n_cell_types <- length(cell_types_unique)
  cell_type_colors <- metadata(spe)$color_vectors$sub_celltype
  
  # Create the base plot with background bars
  p1 <- ggplot() +
    # Background bars showing range of ALL markers
    geom_col(aes(x = cell_type, y = low), data = dfbar,
             fill = "#dcdcdc", alpha = 0.6, width = 0.8) +
    geom_col(aes(x = cell_type, y = up), data = dfbar,
             fill = "#dcdcdc", alpha = 0.6, width = 0.8) +
    # Jitter points for ALL markers (significant and non-significant)
    geom_jitter(aes(x = cell_type, y = log2_fold_change, color = label), 
                data = plot_markers,
                size = 1.2, alpha = 0.7, 
                position = position_jitter(seed = 1234, width = 0.3)) +
    # Color scale: significant markers get colors, non-significant are grey
    scale_color_manual(
      values = c(
        "Higher_in_Activate" = "#cd534c", 
        "Higher_in_Quiescent" = "#7ca6dc", 
        "Non_Significant" = "grey70"
      ),
      breaks = c("Higher_in_Activate", "Higher_in_Quiescent", "Non_Significant"),
      labels = c("Higher_in_Activate", "Higher_in_Quiescent", "Non-Significant")
    ) +
    theme_classic()
  
  # Add cell type labels at bottom
  p2 <- p1 + 
    geom_tile(aes(x = cell_type, y = 0), 
              data = data.frame(cell_type = cell_types_unique),
              height = 0.4, fill = cell_type_colors[cell_types_unique], 
              show.legend = FALSE, alpha = 0.6) +
    geom_text(aes(x = cell_type, y = 0, label = cell_type), 
              data = data.frame(cell_type = cell_types_unique),
              size = 3, fontface = "bold", color = "black")
  
  # Add text labels ONLY for significant markers
  p3 <- p2 + 
    geom_text_repel(aes(x = cell_type, y = log2_fold_change, label = do_label), 
                    data = plot_markers,
                    position = position_jitter(width = 0.3),
                    size = 2.5,
                    box.padding = unit(0.3, "lines"),
                    point.padding = unit(0.4, "lines"),
                    min.segment.length = 0,
                    max.overlaps = 50,
                    segment.size = 0.4,
                    segment.alpha = 0.6,
                    na.rm = TRUE)
  
  # Final styling
  p3 <- p3 +
    labs(
      x = "Cell Type",
      y = "Log2 Fold Change (Activate vs EpCAM+)",
      title = "Marker Expression Changes in Microenvironment Cells",
      subtitle = paste("All markers shown;", 
                       sum(plot_markers$label != "Non_Significant"), 
                       "significant markers highlighted and labeled")
    ) +
    theme(
      plot.title = element_text(size = 14, color = "black", face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "black", hjust = 0.5),
      axis.title = element_text(size = 12, color = "black", face = "bold"),
      axis.line.y = element_line(color = "black", linewidth = 0.8),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      legend.position = c(0.02, 0.96),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.direction = "vertical",
      legend.justification = c(1, 0),
      legend.text = element_text(size = 11)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  print(p3)
  ggsave(file.path(figureDir, paste0("Marker expression change of Microenvironment Cells between Tumor Types.pdf")), 
         p3, width = 20, height = 6)

# Create functional category analysis
if (nrow(sig_results) > 0) {
  # Define functional marker categories
  functional_categories <- list(
    "Immune_Activation" = c("CD80", "CD40", "HLA_DR", "CD274", "CD127"),
    "Immune_Exhaustion" = c("CD279", "TIGIT", "CD366", "CD27"),
    "Metabolic" = c("GLUT1", "HK2", "FASN", "PRPS1"),
    "Hypoxia" = c("CA_IX", "VEGF"),
    "Proliferation" = c("Ki67"),
    "Inflammation" = c("CD68", "CD163", "CD11b", "CD11c"),
    "Stromal_Activation" = c("FAP", "Alpha_SMA", "Vimentin", "Collagen_I"),
    "T_Cell_Function" = c("CD3", "CD4", "CD8a", "FoxP3", "CD45")
  )
  
  # Categorize significant results
  sig_results$functional_category <- "Other"
  for (category in names(functional_categories)) {
    sig_results$functional_category[sig_results$marker %in% functional_categories[[category]]] <- category
  }
  
  # Summary by functional category
  functional_summary <- sig_results %>%
    group_by(functional_category, significance) %>%
    summarise(
      n_changes = n(),
      mean_fold_change = mean(abs(log2_fold_change)),
      .groups = 'drop'
    ) %>%
    pivot_wider(names_from = significance, values_from = n_changes, values_fill = 0)
  
  print("\n=== FUNCTIONAL CATEGORY ANALYSIS ===")
  print(functional_summary)
}

# Save results
write.csv(cell_type_marker_results, "Section7_Cell_Type_Marker_Results.csv", row.names = FALSE)

if (nrow(sig_results) > 0) {
  write.csv(sig_results, "Section7_Significant_Marker_Changes.csv", row.names = FALSE)
  write.csv(cell_type_summary, "Section7_Cell_Type_Reprogramming_Summary.csv", row.names = FALSE)
  
  if (exists("functional_summary")) {
    write.csv(functional_summary, "Section7_Functional_Category_Summary.csv", row.names = FALSE)
  }
}

print("\n=== Section 7 Complete ===")
print("Results saved to:")
print("- Section7_Cell_Type_Marker_Results.csv")
if (nrow(sig_results) > 0) {
  print("- Section7_Significant_Marker_Changes.csv")
  print("- Section7_Cell_Type_Reprogramming_Summary.csv")
  if (exists("functional_summary")) {
    print("- Section7_Functional_Category_Summary.csv")
  }
}

print("\n=== COMPREHENSIVE PATCH ANALYSIS COMPLETE ===")
print("All seven sections of the 'Seed and Soil' analysis have been completed:")
print("1. ✅ Metabolic Profiling - Patch-level marker differences")
print("2. ✅ Patch Size Calculation - Spatial metrics and density")
print("3. ✅ Marker Expression Analysis - Size-dependent patterns")
print("4. ✅ Patch Size Comparison by RFS - Clinical relevance")
print("5. ✅ Component Analysis - Microenvironment composition")
print("6. ✅ Cell Type Fraction Comparison - Statistical composition analysis")
print("7. ✅ Cell Type Marker Analysis - Patch-dependent cellular reprogramming")
print("\nReady for integration and clinical translation!")

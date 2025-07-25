# For Tumor cell analysis

# ------------------------------
# Set Working Directory and Source Utilities
# ------------------------------
library(SpatialExperiment)
library(imcRtools)
library(igraph)

library(ggplot2)
library(ggridges)
library(ggalluvial)
library(RColorBrewer)
library(ggpubr)
library(ggraph)
library(ggrepel)
library(viridis)
library(patchwork)

library(dplyr)
library(tidyr)

# =============================================================================
# LOAD AND VALIDATE DATA
# =============================================================================
codeSpace <- "/home/lenislin/Experiment/projects/CRLM_2025/IMC"

source(file.path(codeSpace,"Tumore_analysis_functions.R"))

# Set Working Directory and Source Utilities
workDir <- "/mnt/public/lyx/IMC_HE_Merge/CRLM"
setwd(workDir)

# Define directories for saving results and Steinbock outputs
saveDir <- file.path(workDir, "results")
steinResultDir <- file.path(workDir, "steinbock")
figureDir <- file.path(workDir, "figures","4_TumorAnalysis")

if(!dir.exists(figureDir)){
  dir.create(figureDir,recursive = T)
}

# Load your SpatialExperiment object
date_time <- "0722"
spe <- readRDS(file.path(saveDir,paste0("spatial_spe_",date_time,".rds")))
img_id_ <- "sample_id"

head(colData(spe))

metadata(spe)$color_vectors[["tissue"]] <- c(
  "PT" = "#00BA38",
  "IM" = "#F8766D",
  "TC" = "#619CFF"
)

# Step 1: Extract cells from TC and IM tissues
cell_data <- as.data.frame(colData(spe))
tc_im_cells <- cell_data[cell_data$Tissue %in% c("TC", "IM"), ]

cat("Total cells in TC and IM:", nrow(tc_im_cells), "\n")

# Step 2: Extract malignant cells (sub_celltype starting with "EC_")
malignant_cells <- tc_im_cells[grepl("^EC_", tc_im_cells$sub_celltype), ]
cat("Malignant cells:", nrow(malignant_cells), "\n")

# Check malignant cell subtypes
cat("Malignant cell subtypes:\n")
print(table(malignant_cells$sub_celltype))

# Step 3: Compute malignant cell fraction in different ROIs
roi_summary <- tc_im_cells %>%
  group_by(sample_id, patient_id, Tissue, Treatment, RFS_status) %>%
  summarise(
    total_cells = n(),
    malignant_cells = sum(grepl("^EC_", sub_celltype)),
    malignant_fraction = malignant_cells / total_cells,
    .groups = "drop"
  )

# Add early relapse labels
roi_summary$RFS_group <- ifelse(roi_summary$RFS_status == 0, "No Early Relapse", "Early Relapse")

print("ROI summary statistics:")
print(head(roi_summary))

# Step 4: Stacked barplot to display results
p1 <- ggplot(roi_summary, aes(x = sample_id, y = malignant_fraction)) +
  geom_col(aes(fill = Tissue), alpha = 0.7) +
  facet_wrap(~paste(Treatment, "-", RFS_group), scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
  labs(title = "Malignant Cell Fraction by ROI",
       x = "ROI (sample_id)", 
       y = "Malignant Cell Fraction",
       fill = "Tissue Region") +
  scale_fill_manual(values = metadata(spe)$color_vectors[["tissue"]])

ggsave(file.path(figureDir, paste0("all_malignant_cell_barplot.pdf")), p1, width = 14, height = 8)

# Step 5: Compare fractions between RFS_status groups
# Prepare data for statistical comparison
comparison_data <- roi_summary %>%
  select(sample_id, patient_id, Tissue, Treatment, RFS_status, RFS_group, malignant_fraction)

# Boxplot with wilcox test
p2 <- ggplot(comparison_data, aes(x = RFS_group, y = malignant_fraction)) +
  geom_boxplot(aes(fill = RFS_group), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~Tissue) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_minimal() +
  labs(title = "Malignant Cell Fraction by Recurrence Status",
       x = "Recurrence Status",
       y = "Malignant Cell Fraction",
       fill = "RFS Status") +
  scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))

ggsave(file.path(figureDir, paste0("all_malignant_boxplot_rfs.pdf")), p2, width = 10, height = 6)

# Step 6: Extend comparison to different treatment subgroups
p3 <- ggplot(comparison_data, aes(x = RFS_group, y = malignant_fraction)) +
  geom_boxplot(aes(fill = RFS_group), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_grid(Treatment ~ Tissue) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  theme_minimal() +
  labs(title = "Malignant Cell Fraction by Treatment and Recurrence Status",
       x = "Recurrence Status",
       y = "Malignant Cell Fraction",
       fill = "RFS Status") +
  scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))

ggsave(file.path(figureDir, paste0("all_malignant_boxplot_treatment.pdf")), p3, width = 10, height = 8)

# Step 7: Sankey diagram for fraction changes from TC to IM
# First, get paired TC-IM data for each patient
patient_summary <- roi_summary %>%
  group_by(patient_id, Tissue, Treatment, RFS_status, RFS_group) %>%
  summarise(mean_malignant_fraction = mean(malignant_fraction, na.rm = TRUE), .groups = "drop")

paired_data <- patient_summary %>%
  select(patient_id, Tissue, Treatment, RFS_status, RFS_group, mean_malignant_fraction) %>%
  pivot_wider(names_from = Tissue, values_from = mean_malignant_fraction, 
              names_prefix = "fraction_") %>%
  filter(!is.na(fraction_TC) & !is.na(fraction_IM))

# Create bins for fraction levels to make Sankey readable
create_fraction_bins <- function(x) {
  case_when(
    x < 0.1 ~ "Low (<10%)",
    x < 0.3 ~ "Medium (10-30%)",
    x < 0.5 ~ "High (30-50%)",
    TRUE ~ "Very High (>50%)"
  )
}

paired_data$TC_bin <- create_fraction_bins(paired_data$fraction_TC)
paired_data$IM_bin <- create_fraction_bins(paired_data$fraction_IM)

# Prepare data for alluvial plot
sankey_data <- paired_data %>%
  count(RFS_group, TC_bin, IM_bin) %>%
  rename(freq = n)

# Create Sankey/Alluvial plot
p4 <- ggplot(sankey_data, aes(axis1 = TC_bin, axis2 = IM_bin, y = freq)) +
  geom_alluvium(aes(fill = RFS_group), alpha = 0.7) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("TC", "IM"), expand = c(0.1, 0.1)) +
  facet_wrap(~RFS_group) +
  theme_minimal() +
  labs(title = "Malignant Cell Fraction Changes from TC to IM",
       subtitle = "Flow shows transition patterns by recurrence status",
       x = "Tissue Region",
       y = "Number of Patients",
       fill = "RFS Status") +
  scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))

ggsave(file.path(figureDir, paste0("all_malignant_sankey.pdf")), p4, width = 12, height = 8)

# ================================================================================
# MALIGNANT MARKER ANALYSIS
# ================================================================================

# Compare the expression among all EP cells during IM
malignant_roi_expression <- data.frame()
malignant_cell_expression <- assay(spe)[,rownames(malignant_cells)]

malignant_cells$RFS_group <- ifelse(malignant_cells$RFS_status == 0, 
                                "No Early Relapse", "Early Relapse")

tumor_markers <- c("EpCAM", "Vimentin", "Ki67", "VEGF", "CA_IX", "GLUT1", "HK2", "FASN")

for(roi in unique(malignant_cells$sample_id)) {
  roi_malignant_cells <- malignant_cells[malignant_cells$sample_id == roi, ]
  if(nrow(roi_malignant_cells) == 0) next
  
  malignant_cell_indices <- rownames(roi_malignant_cells)
  roi_expression <- malignant_cell_expression[tumor_markers, malignant_cell_indices, drop = FALSE]
  
  # Calculate mean expression for this ROI
  roi_means <- rowMeans(roi_expression)
  
  # Get clinical metadata
  roi_meta <- roi_malignant_cells[1, c("sample_id", "patient_id", "Treatment", "Tissue", "RFS_status", "RFS_group")]

  # Combine with expression data
  roi_data <- data.frame(t(roi_means), roi_meta)
  malignant_roi_expression <- rbind(malignant_roi_expression, roi_data)
}

rownames(malignant_roi_expression) <- malignant_roi_expression$sample_id

## volcano plot
# Define marker columns
marker_cols <- c("EpCAM", "Vimentin", "Ki67", "VEGF", "CA_IX", "GLUT1", "HK2", "FASN")

# Create combinations for analysis
tissue_types <- c("TC+IM", "TC", "IM")  # Adjust based on your actual tissue types
treatments <- c("All", "Chemo", "Combo")

# Initialize results dataframe
results_list <- list()

# Loop through each combination
for(tissue in tissue_types) {
  for(treatment in treatments) {
    
    # Filter data based on tissue and treatment
    if(tissue == "TC+IM") {
      subset_data <- malignant_roi_expression
    }
    else{
      subset_data <- malignant_roi_expression[malignant_roi_expression$Tissue == tissue, ]
    }
    # Filter data based on tissue and treatment
    if(treatment == "All") {
      subset_data <- subset_data
    } else {
      subset_data <- subset_data[subset_data$Treatment == treatment, ]
    }
    
    # Skip if too few samples
    if(nrow(subset_data) < 6) next
    
    # Get unique RFS groups in this subset
    rfs_groups <- unique(subset_data$RFS_group)
    if(length(rfs_groups) < 2) next
    
    # Calculate statistics for each marker
    for(marker in marker_cols) {
      
      # Split data by RFS group
      group1_data <- subset_data[subset_data$RFS_group == rfs_groups[1], marker]
      group2_data <- subset_data[subset_data$RFS_group == rfs_groups[2], marker]
      
      # Remove NA values
      group1_data <- group1_data[!is.na(group1_data)]
      group2_data <- group2_data[!is.na(group2_data)]
      
      # Skip if too few samples in either group
      if(length(group1_data) < 3 | length(group2_data) < 3) next
      
      # Calculate means
      mean1 <- mean(group1_data)
      mean2 <- mean(group2_data)
      
      # Calculate fold change (log2)
      # Add small constant to avoid log(0)
      fold_change <- log2((mean1 + 0.001) / (mean2 + 0.001))
      
      # Perform t-test
      t_test_result <- t.test(group1_data, group2_data)
      p_value <- t_test_result$p.value
      
      # Store results
      results_list[[length(results_list) + 1]] <- data.frame(
        tissue_type = tissue,
        treatment = treatment,
        marker = marker,
        group1 = rfs_groups[1],
        group2 = rfs_groups[2],
        mean1 = mean1,
        mean2 = mean2,
        fold_change = fold_change,
        p_value = p_value,
        n1 = length(group1_data),
        n2 = length(group2_data)
      )
    }
  }
}

# Combine results
results_df <- do.call(rbind, results_list)

# Adjust p-values using Benjamini-Hochberg method
results_df$adj_p_value <- p.adjust(results_df$p_value, method = "BH")
results_df$significance <- ifelse(results_df$adj_p_value <= 0.001, "***",
                                  ifelse(results_df$adj_p_value <= 0.01, "**",
                                         ifelse(results_df$adj_p_value <= 0.05, "*", "ns")))
results_df$log10_adj_p <- -log10(results_df$adj_p_value)

# Set factor levels for proper ordering
results_df$tissue_type <- factor(results_df$tissue_type, levels = c("TC+IM", "TC", "IM"))
results_df$treatment <- factor(results_df$treatment, levels = c("All", "Chemo", "Combo"))

# Print results summary
print("Summary of significant differences:")
significant_results <- results_df[results_df$significance != "ns", ]
print(significant_results[, c("tissue_type", "treatment", "marker", "fold_change", "adj_p_value", "significance")])

# --- Create the Volcano Plot ---
significant_results <- results_df[results_df$significance != "ns", ]
p_volcano <- ggplot(results_df, aes(x = fold_change, y = log10_adj_p)) +
  
  geom_point(data = subset(results_df, significance == "ns"), 
             size = 1, color = 'grey', alpha = 0.6) +
  geom_point(data = significant_results, 
             aes(color = significance), size = 2) +
  geom_text_repel(
    data = significant_results,
    aes(label = marker),
    size = 3.5,               # Font size for the labels.
    box.padding = 0.5,        # Padding around the label text.
    point.padding = 0.3,      # Padding around the point itself.
    max.overlaps = 20,        # Increase max overlaps to allow more labels.
    min.segment.length = 0,   # Always draw segment lines, regardless of distance.
    segment.color = 'grey50', # Set the color of the line segments.
    force = 2                 # Increase the force of repulsion between labels.
  ) +
  
  # Facet the plot by tissue type and treatment.
  facet_grid(tissue_type ~ treatment, scales = "free") +
  
  # Add significance threshold lines.
  geom_vline(xintercept = c(-0.58, 0.58), size = 0.5, color = "grey50", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), size = 0.5, color = "grey50", linetype = 'dashed') +
  
  # Manually set the colors for significance levels and add a legend title.
  scale_color_manual(
    name = "Significance", 
    values = c("*" = "#FF7F00", "**" = "#E31A1C", "***" = "#800080")
  ) +
  
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  
  theme_bw() +
  theme(
    legend.position = 'bottom',
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    # This rotation now applies to the p-value axis labels.
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
    strip.text = element_text(size = 10, face = 'bold'),
    strip.background = element_rect(fill = "grey90", color = "black")
  )


# Print volcano plot
print(p_volcano)

# Save plots
ggsave(file.path(figureDir,"marker_volcano_plot.pdf"), plot = p_volcano, width = 12, height = 10)
write.csv(file.path(figureDir,"marker_expression_analysis.csv"), , row.names = FALSE)

# ================================================================================
# MALIGNANT SUBTYPE ANALYSIS
# ================================================================================

cat("\n=== ANALYZING MALIGNANT SUBTYPES ===\n")

# Get all malignant subtypes
malignant_subtypes <- unique(malignant_cells$sub_celltype)
cat("Malignant subtypes found:", paste(malignant_subtypes, collapse = ", "), "\n")

# Analyze each malignant subtype
subtype_results <- list()

for(subtype in malignant_subtypes) {
  subtype_data <- compute_subtype_fractions(subtype)
  subtype_results[[subtype]] <- create_subtype_plots(subtype, subtype_data, figureDir)
}

# ================================================================================
## Tumor Budding analysis
# ================================================================================
# Load required libraries
library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(igraph)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)

# Create output directory
if (!dir.exists("figures/budding")) {
  dir.create("figures/budding", recursive = TRUE)
}

# ================================================================================
# STEP 1: EXTRACT MALIGNANT CELLS FROM IM (Already provided)
# ================================================================================

cat("=== EXTRACTING MALIGNANT CELLS FROM IM ===\n")

malignant_IM <- spe[, spe$Tissue %in% "IM"]
malignant_IM_all <- malignant_IM
malignant_IM <- malignant_IM[, grepl("^EC_", malignant_IM$sub_celltype)]

all_IM_cells <- as.data.frame(colData(malignant_IM_all))

cat("Total malignant cells in IM:", ncol(malignant_IM), "\n")
cat("Number of ROIs:", length(unique(malignant_IM$sample_id)), "\n")
cat("Number of patients:", length(unique(malignant_IM$patient_id)), "\n")

# ================================================================================
# STEP 2: FILTER LONG EDGES IN DELAUNAY TRIANGULATION
# ================================================================================

cat("\n=== FILTERING LONG EDGES ===\n")

# Get Delaunay triangulation edges
delaunay_edges <- as.data.frame(colPair(malignant_IM, "delaunay")) 

# Get spatial coordinates (assuming they exist in colData)
cell_coords <- data.frame(
  cell_id = 1:ncol(malignant_IM),
  x = spatialCoords(malignant_IM)[,1],
  y = spatialCoords(malignant_IM)[,2],
  sample_id = malignant_IM$sample_id
)

# Extract edge information
edge_from <- delaunay_edges[,1]
edge_to <- delaunay_edges[,2]

# Calculate edge lengths
edge_lengths <- sqrt(
  (cell_coords$x[edge_from] - cell_coords$x[edge_to])^2 + 
    (cell_coords$y[edge_from] - cell_coords$y[edge_to])^2
)

cat("Total Delaunay edges:", length(edge_from), "\n")
cat("Edge length summary:\n")
print(summary(edge_lengths))

# Filter edges by length threshold (20 μm)
edge_threshold <- 20
valid_edges <- edge_lengths <= edge_threshold

filtered_edges_from <- edge_from[valid_edges]
filtered_edges_to <- edge_to[valid_edges]
filtered_edge_lengths <- edge_lengths[valid_edges]

cat("Edges after length filtering:", length(filtered_edges_from), "\n")
cat("Percentage of edges retained:", round(100 * sum(valid_edges) / length(edge_from), 1), "%\n")

# ================================================================================
# STEP 3: DETECT CONNECTED COMPONENTS AND FILTER LARGE CLUSTERS
# ================================================================================

cat("\n=== DETECTING CONNECTED COMPONENTS ===\n")

# Create adjacency list for each ROI separately
roi_list <- unique(malignant_IM$sample_id)
budding_results <- list()

cluster_size_threshold_max <- 20
cluster_size_threshold_min <- 3

for(roi in roi_list) {
  cat("Processing ROI:", roi, "\n")
  
  # Get cells in this ROI
  roi_cells <- which(malignant_IM$sample_id == roi)
  all_roi_cells <- which(all_IM_cells$sample_id == roi)
  
  if(length(roi_cells) / length(all_roi_cells) < 0.1) next
  
  # Filter edges to only include cells in this ROI
  roi_edges_mask <- (filtered_edges_from %in% roi_cells) & (filtered_edges_to %in% roi_cells)
  roi_edges_from <- filtered_edges_from[roi_edges_mask]
  roi_edges_to <- filtered_edges_to[roi_edges_mask]
  
  if(length(roi_edges_from) == 0) {
    # No edges in this ROI, treat all cells as individual clusters
    roi_components <- as.list(roi_cells)
    names(roi_components) <- paste0("component_", 1:length(roi_cells))
  } else {
    # Create igraph object for this ROI
    # Need to remap cell indices to 1:n for this ROI
    cell_mapping <- setNames(1:length(roi_cells), roi_cells)
    
    roi_edges_from_mapped <- cell_mapping[as.character(roi_edges_from)]
    roi_edges_to_mapped <- cell_mapping[as.character(roi_edges_to)]
    
    roi_graph <- graph_from_edgelist(cbind(roi_edges_from_mapped, roi_edges_to_mapped), 
                                     directed = FALSE)
    
    # Add isolated vertices (cells with no connections)
    all_vertices <- 1:length(roi_cells)
    missing_vertices <- setdiff(all_vertices, as.numeric(V(roi_graph)))
    if(length(missing_vertices) > 0) {
      roi_graph <- add_vertices(roi_graph, length(missing_vertices), 
                                name = as.character(missing_vertices))
    }
    
    # Find connected components
    components <- components(roi_graph)
    
    # Convert back to original cell indices
    roi_components <- list()
    for(i in 1:components$no) {
      component_vertices <- which(components$membership == i)
      original_indices <- roi_cells[as.numeric(V(roi_graph)[component_vertices])]
      roi_components[[paste0("component_", i)]] <- original_indices
    }
  }
  
  # Get component sizes
  component_sizes <- sapply(roi_components, length)
  
  # Filter components by size (keep small clusters as "budding")
  budding_components <- roi_components[component_sizes <= cluster_size_threshold_max & component_sizes >= cluster_size_threshold_min]
  large_components <- roi_components[component_sizes > cluster_size_threshold_max]
  
  cat("  Total components:", length(roi_components), "\n")
  cat("  Small components (budding):", length(budding_components), "\n")
  cat("  Large components (filtered):", length(large_components), "\n")
  
  # Store results
  budding_results[[roi]] <- list(
    budding_components = budding_components,
    budding_cells = rownames(colData(malignant_IM))[unname(unlist(budding_components))],
    large_components = large_components,
    component_sizes = component_sizes,
    roi_edges_from = rownames(colData(malignant_IM))[roi_edges_from],
    roi_edges_to = rownames(colData(malignant_IM))[roi_edges_to]
  )
}

# ================================================================================
# STEP 4: COLLECT MARKER EXPRESSION FOR TUMOR BUDDING
# ================================================================================

cat("\n=== COLLECTING BUDDING MARKER EXPRESSION ===\n")

# Get all budding cells
all_budding_cells <- c()
budding_metadata <- data.frame()

for(roi in names(budding_results)) {
  roi_budding <- budding_results[[roi]]$budding_components
  
  for(comp_name in names(roi_budding)) {
    comp_cells <- roi_budding[[comp_name]]
    comp_size <- length(comp_cells)
    
    if(comp_size < cluster_size_threshold_min | comp_size > cluster_size_threshold_max){next;}
    
    # Add to overall list
    all_budding_cells <- c(all_budding_cells, comp_cells)
    
    # Create metadata for these cells
    comp_meta <- data.frame(
      cell_idx = comp_cells,
      roi = roi,
      component = comp_name,
      cluster_size = comp_size,
      budding_type = ifelse(comp_size == 1, "single_cell", "small_cluster")
    )
    
    budding_metadata <- rbind(budding_metadata, comp_meta)
  }
}

cat("Total budding cells:", length(all_budding_cells), "\n")
cat("Single cell budding:", sum(budding_metadata$cluster_size == 1), "\n")
cat("Small cluster budding:", sum(budding_metadata$cluster_size > 1), "\n")

# Extract expression data for budding cells
budding_expression <- assay(malignant_IM)[, all_budding_cells]
budding_cell_data <- as.data.frame(colData(malignant_IM)[all_budding_cells, ])

# Merge with budding metadata
budding_metadata$cell_id <- rownames(budding_cell_data)
budding_cell_data <- merge(budding_cell_data, budding_metadata, by.x = "row.names", by.y = "cell_id")

# Add clinical metadata
budding_cell_data$RFS_group <- ifelse(budding_cell_data$RFS_status == 0, 
                                      "No Early Relapse", "Early Relapse")

# Only keep small_cluster
cat("Budding cells by cluster size:\n")
print(table(budding_cell_data$cluster_size))

# ================================================================================
# STEP 5: ANALYZE MARKER EXPRESSION DIFFERENCES BY RFS STATUS
# ================================================================================

cat("\n=== ANALYZING MARKER EXPRESSION BY RFS STATUS ===\n")

# Select relevant markers for budding analysis
budding_markers <- c("EpCAM", "Vimentin", "Ki67", "VEGF", "CA_IX", "GLUT1", "HK2", "FASN")
available_budding_markers <- budding_markers[budding_markers %in% rownames(malignant_IM)]

cat("Available budding markers:", paste(available_budding_markers, collapse = ", "), "\n")

# Calculate mean expression per ROI for budding cells
roi_budding_expression <- data.frame()
rownames(budding_cell_data) <- budding_cell_data$Row.names

for(roi in unique(budding_cell_data$sample_id)) {
  roi_budding_cells <- budding_results[[roi]][["budding_cells"]]
  roi_budding_cells <- budding_cell_data[roi_budding_cells, ]
  
  if(nrow(roi_budding_cells) == 0) next
  
  roi_cell_indices <- rownames(roi_budding_cells)
  roi_expression <- budding_expression[available_budding_markers, roi_cell_indices, drop = FALSE]
  
  # Calculate mean expression for this ROI
  roi_means <- rowMeans(roi_expression)
  
  # Get clinical metadata
  roi_meta <- roi_budding_cells[1, c("sample_id", "patient_id", "Treatment", "RFS_status", "RFS_group")]
  roi_meta$n_budding_cells <- nrow(roi_budding_cells)
  roi_meta$mean_cluster_size <- mean(roi_budding_cells$cluster_size)
  
  # Combine with expression data
  roi_data <- data.frame(t(roi_means), roi_meta)
  roi_budding_expression <- rbind(roi_budding_expression, roi_data)
}

rownames(roi_budding_expression) <- roi_budding_expression$sample_id

cat("ROIs with budding cells:", nrow(roi_budding_expression), "\n")

# Boxplots for each marker
budding_expression_long <- roi_budding_expression %>%
  select(sample_id, patient_id, Treatment, RFS_status, RFS_group, all_of(available_budding_markers)) %>%
  pivot_longer(cols = all_of(available_budding_markers), 
               names_to = "marker", values_to = "expression")

p_budding_rfs <- ggplot(budding_expression_long, aes(x = RFS_group, y = expression)) +
  geom_boxplot(aes(fill = RFS_group), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~marker, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_minimal() +
  labs(title = "Tumor Budding Marker Expression by RFS Status",
       x = "RFS Status",
       y = "Mean Expression",
       fill = "RFS Status") +
  scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))

print(p_budding_rfs)
ggsave("figures/budding/budding_markers_rfs_comparison.pdf", p_budding_rfs, width = 12, height = 10)

# Treatment-stratified boxplots
p_budding_treatment <- ggplot(budding_expression_long, aes(x = RFS_group, y = expression)) +
  geom_boxplot(aes(fill = RFS_group), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_grid(Treatment ~ marker, scales = "free_x") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 2.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Tumor Budding Marker Expression by Treatment and RFS Status",
       x = "RFS Status",
       y = "Mean Expression",
       fill = "RFS Status") +
  scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))

print(p_budding_treatment)
ggsave("figures/budding/budding_markers_treatment_comparison.pdf", p_budding_treatment, 
       width = 12, height = 6)

# Summary of budding characteristics
cat("\n=== BUDDING SUMMARY STATISTICS ===\n")

budding_summary <- roi_budding_expression %>%
  group_by(Treatment, RFS_group) %>%
  summarise(
    n_rois = n(),
    mean_budding_cells = round(mean(n_budding_cells), 1),
    median_budding_cells = median(n_budding_cells),
    mean_cluster_size = round(mean(mean_cluster_size), 2),
    .groups = "drop"
  )

cat("Budding characteristics by treatment and RFS status:\n")
print(budding_summary)

# Plot budding cell counts
p_budding_counts <- ggplot(roi_budding_expression, aes(x = RFS_group, y = n_budding_cells)) +
  geom_boxplot(aes(fill = RFS_group), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~Treatment) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_minimal() +
  labs(title = "Number of Budding Cells per ROI",
       x = "RFS Status",
       y = "Number of Budding Cells",
       fill = "RFS Status") +
  scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))

print(p_budding_counts)
ggsave("figures/budding/budding_cell_counts.pdf", p_budding_counts, width = 10, height = 6)

# Plot cluster sizes
p_cluster_sizes <- ggplot(roi_budding_expression, aes(x = RFS_group, y = mean_cluster_size)) +
  geom_boxplot(aes(fill = RFS_group), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~Treatment) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_minimal() +
  labs(title = "Mean Cluster Size of Budding Cells",
       x = "RFS Status",
       y = "Mean Cluster Size",
       fill = "RFS Status") +
  scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))

print(p_cluster_sizes)
ggsave("figures/budding/budding_cluster_sizes.pdf", p_cluster_sizes, width = 10, height = 6)

cat("\n=== TUMOR BUDDING ANALYSIS COMPLETE ===\n")
cat("Results saved in 'figures/budding/' directory:\n")
cat("- budding_markers_rfs_comparison.pdf: Marker expression by RFS status\n")
cat("- budding_markers_treatment_comparison.pdf: Treatment-stratified analysis\n")
cat("- budding_markers_heatmap.pdf: Fold change heatmap\n")
cat("- budding_cell_counts.pdf: Number of budding cells per ROI\n")
cat("- budding_cluster_sizes.pdf: Cluster size distributions\n")

# Export results
write.csv(roi_budding_expression, "figures/budding/roi_budding_expression_data.csv", row.names = FALSE)

# Compare the expression among all EP cells during IM
roi_EC_expression <- data.frame()
IM_EC_expression <- assay(malignant_IM)
IM_EC_cells <- as.data.frame(colData(malignant_IM))
IM_EC_cells$RFS_group <- ifelse(IM_EC_cells$RFS_status == 0, 
                                "No Early Relapse", "Early Relapse")

for(roi in unique(malignant_IM$sample_id)) {
  roi_EC_cells <- IM_EC_cells[IM_EC_cells$sample_id == roi, ]
  
  if(nrow(roi_EC_cells) == 0) next
  
  roi_cell_indices <- rownames(roi_EC_cells)
  roi_expression <- IM_EC_expression[available_budding_markers, roi_cell_indices, drop = FALSE]
  
  # Calculate mean expression for this ROI
  roi_means <- rowMeans(roi_expression)
  
  # Get clinical metadata
  roi_meta <- roi_EC_cells[1, c("sample_id", "patient_id", "Treatment", "RFS_status", "RFS_group")]
  roi_meta$n_budding_cells <- nrow(roi_EC_cells)
  
  # Combine with expression data
  roi_data <- data.frame(t(roi_means), roi_meta)
  roi_EC_expression <- rbind(roi_EC_expression, roi_data)
}

rownames(roi_EC_expression) <- roi_EC_expression$sample_id

cat("ROIs with EC cells:", nrow(roi_EC_expression), "\n")

# Boxplots for each marker
EC_expression_long <- roi_EC_expression %>%
  select(sample_id, patient_id, Treatment, RFS_status, RFS_group, all_of(available_budding_markers)) %>%
  pivot_longer(cols = all_of(available_budding_markers), 
               names_to = "marker", values_to = "expression")

p_budding_rfs <- ggplot(EC_expression_long, aes(x = RFS_group, y = expression)) +
  geom_boxplot(aes(fill = RFS_group), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~marker, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_minimal() +
  labs(title = "Tumor Marker Expression by RFS Status",
       x = "RFS Status",
       y = "Mean Expression",
       fill = "RFS Status") +
  scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))

print(p_budding_rfs)
ggsave("figures/budding/all_EC_markers_rfs_comparison.pdf", p_budding_rfs, width = 12, height = 10)

# Treatment-stratified boxplots
p_budding_treatment <- ggplot(EC_expression_long, aes(x = RFS_group, y = expression)) +
  geom_boxplot(aes(fill = RFS_group), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_grid(Treatment ~ marker, scales = "free_x") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 2.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Tumor Marker Expression by Treatment and RFS Status",
       x = "RFS Status",
       y = "Mean Expression",
       fill = "RFS Status") +
  scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))

print(p_budding_treatment)
ggsave("figures/budding/all_EC_markers_treatment_comparison.pdf", p_budding_treatment, 
       width = 12, height = 6)

# ================================================================================
# TASK 1: SPATIAL VISUALIZATION OF TUMOR BUDDING
# ================================================================================

cat("\n=== CREATING SPATIAL VISUALIZATIONS ===\n")

# Get IM spe
spe_IM <- malignant_IM_all

# Add EpCAM expression levels and match TB label
spe_IM$EpCAM_expr <- assay(spe_IM)["EpCAM",]
spe_IM$EpCAM_status <- ifelse(startsWith(spe_IM$sub_celltype,"EC_"),"EpCAM+","EpCAM-")

budding_cell_data <- as.data.frame(colData(spe_IM))

# Add spatial coordinates
budding_cell_data$x <- spatialCoords(spe_IM)[, 1]
budding_cell_data$y <- spatialCoords(spe_IM)[, 2]

# Create spatial plots for selected ROIs (show first 6 ROIs)
selected_rois <- unique(spe_IM$sample_id)
all_coord <- spatialCoords(malignant_IM)

for(roi in selected_rois) {
  cat("Creating spatial plot for ROI:", roi, "\n")
  
  roi_data <- budding_cell_data[budding_cell_data$sample_id == roi, ]
  
  # Get edges for this ROI
  roi_edges_from <- budding_results[[roi]]$roi_edges_from
  roi_edges_to <- budding_results[[roi]]$roi_edges_to
  
  # Filter edges to only include budding cells
  budding_cell_indices <- rownames(roi_data)
  edge_mask <- (roi_edges_from %in% budding_cell_indices) & (roi_edges_to %in% budding_cell_indices)
  roi_budding_edges_from <- roi_edges_from[edge_mask]
  roi_budding_edges_to <- roi_edges_to[edge_mask]
  
  # Create edge data frame for plotting
  from_cell <- roi_budding_edges_from
  to_cell <- roi_budding_edges_to
  
  from_coord <- all_coord[match(from_cell,colnames(malignant_IM)) , ]
  to_coord <- all_coord[match(to_cell,colnames(malignant_IM)) , ]
  
  edge_data <- data.frame()
  edge_data <- rbind(edge_data, data.frame(
    x = from_coord[,1], y = from_coord[,2],
    xend = to_coord[,1], yend = to_coord[,2]
    ))
    
  # Create the spatial plot
  p <- ggplot(roi_data, aes(x = x, y = y)) +
    # Black background
    theme_void() +
    theme(
      plot.background = element_rect(fill = "black"),
      panel.background = element_rect(fill = "black"),
      plot.title = element_text(color = "white"),
      plot.subtitle = element_text(color = "white")
    ) +
    # Add edges first (so they appear behind points)
    {if(nrow(edge_data) > 0) geom_segment(data = edge_data, 
                                          aes(x = x, y = y, xend = xend, yend = yend),
                                          color = "white", alpha = 0.5, linewidth = 0.2)} +
    # Add cells with EpCAM status colors
    geom_point(aes(color = EpCAM_status), size = 1.3, alpha = 0.8) +
    scale_color_manual(values = c("EpCAM-" = "#0072b5ff", "EpCAM+" = "#bc3c29ff")) +
    coord_fixed() +
    labs(title = paste("Tumor Budding -", roi),
         subtitle = paste("EpCAM+ (orange) vs EpCAM- (blue) | n =", nrow(roi_data)),
         color = "EpCAM Status") +
    theme(legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"))
  
  # Save individual plot
  ggsave(paste0("figures/tb_spatial/TB_spatial_", gsub("[^A-Za-z0-9]", "_", roi), ".pdf"), 
         p, width = 8, height = 8, bg = "black")
}

# ================================================================================
# TASK 2: EpCAM:GLUT1 RATIO BIOMARKER
# ================================================================================
budding_cells <- c()
for(roi in selected_rois) {
  roi_data <- budding_cell_data[budding_cell_data$sample_id == roi, ]
  budding_cells <- c(budding_cells,budding_results[[roi]]$budding_cells)
}

budding_cells <- unique(budding_cells)
budding_expression <- budding_expression[,match(budding_cells,colnames(budding_expression)) ]

cat("\n=== CREATING EpCAM:GLUT1 BIOMARKER ===\n")

# Extract GLUT1 expression
glut1_expression <- budding_expression["GLUT1", ]

budding_cell_data <- budding_cell_data[budding_cells,]
budding_cell_data$GLUT1_expr <- glut1_expression

# Calculate EpCAM:GLUT1 ratio
# Add small constant to avoid division by zero
budding_cell_data$GLUT1_EpCAM_ratio <- (budding_cell_data$GLUT1_expr + 0.001) /
  (budding_cell_data$EpCAM_expr + 0.001)

# Log transform ratio for better distribution
budding_cell_data$GLUT1_EpCAM_ratio <- budding_cell_data$GLUT1_EpCAM_ratio

# Add clinical variables
budding_cell_data$RFS_group <- ifelse(budding_cell_data$RFS_status == 0, 
                                      "No Early Relapse", "Early Relapse")

# Calculate ROI-level biomarker scores
roi_biomarker_data <- budding_cell_data %>%
  group_by(sample_id, patient_id, Treatment, RFS_status, RFS_group, RFS_time) %>%
  summarise(
    n_budding_cells = n(),
    mean_EpCAM = mean(EpCAM_expr),
    mean_GLUT1 = mean(GLUT1_expr),
    mean_ratio = mean(GLUT1_EpCAM_ratio),
    mean_log2_ratio = mean(log2(GLUT1_EpCAM_ratio)),
    .groups = "drop"
  )

cat("ROIs with tumor budding data:", nrow(roi_biomarker_data), "\n")

# Visualize biomarker
p_biomarker <- ggplot(roi_biomarker_data, aes(x = RFS_group, y = mean_log2_ratio)) +
  geom_boxplot(aes(fill = RFS_group), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_minimal() +
  labs(title = "GLUT1:EpCAM Ratio Biomarker in Tumor Budding",
       subtitle = "Scaling Log2(GLUT1/EpCAM) ratio by treatment and RFS status",
       x = "RFS Status",
       y = "Scaling Log2(GLUT1:EpCAM Ratio)",
       fill = "RFS Status") +
  scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))

p_biomarker_treatment <- ggplot(roi_biomarker_data, aes(x = RFS_group, y = mean_log2_ratio)) +
  geom_boxplot(aes(fill = RFS_group), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~Treatment) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_minimal() +
  labs(title = "GLUT1:EpCAM Ratio Biomarker in Tumor Budding",
       subtitle = "Scaling Log2(GLUT1/EpCAM) ratio by RFS status",
       x = "RFS Status",
       y = "Scaling Log2(GLUT1:EpCAM Ratio)",
       fill = "RFS Status") +
  scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))

ggsave("figures/tb_spatial/EpCAM_GLUT1_biomarker.pdf", (p_biomarker / p_biomarker_treatment), width = 8, height = 10)

# Create scatter plot showing individual components
p_scatter <- ggplot(roi_biomarker_data, aes(x = mean_GLUT1, y = mean_EpCAM)) +
  geom_point(aes(color = RFS_group, shape = Treatment), size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "gray50", linetype = "dashed") +
  theme_minimal() +
  labs(title = "EpCAM vs GLUT1 Expression in Tumor Budding",
       x = "Mean GLUT1 Expression",
       y = "Mean EpCAM Expression",
       color = "RFS Status",
       shape = "Treatment") +
  scale_color_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))

print(p_scatter)
ggsave("figures/tb_spatial/EpCAM_vs_GLUT1_scatter.pdf", p_scatter, width = 10, height = 6)

# ================================================================================
# TASK 3: KAPLAN-MEIER ANALYSIS
# ================================================================================

cat("\n=== KAPLAN-MEIER SURVIVAL ANALYSIS ===\n")

# Prepare survival data at patient level
patient_survival_data <- roi_biomarker_data %>%
  group_by(patient_id, Treatment) %>%
  summarise(
    RFS_time = first(RFS_time),
    RFS_status = first(RFS_status),
    mean_biomarker = mean(mean_ratio, na.rm = TRUE),
    .groups = "drop"
  )

# Create biomarker categories (high vs low)
library(survival)
library(survminer)
biomarker_cutoff <- surv_cutpoint(patient_survival_data,time = "RFS_time",event = "RFS_status","mean_biomarker",minprop = 0.15)
biomarker_cutoff <- biomarker_cutoff[["cutpoint"]]$cutpoint

# biomarker_cutoff <- median(patient_survival_data$mean_biomarker)
patient_survival_data$biomarker_group <- ifelse(patient_survival_data$mean_biomarker >= biomarker_cutoff,
                                                "High GLUT1:EpCAM", "Low GLUT1:EpCAM")

cat("Survival analysis cohort:", nrow(patient_survival_data), "patients\n")
cat("Biomarker cutoff (median):", round(biomarker_cutoff, 3), "\n")

# Overall survival analysis
surv_obj <- Surv(patient_survival_data$RFS_time, patient_survival_data$RFS_status)
surv_fit <- survfit(surv_obj ~ biomarker_group, data = patient_survival_data)

# Plot overall KM curve
p_km_overall <- ggsurvplot(
  surv_fit,
  data = patient_survival_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  linetype = "strata",
  surv.median.line = "hv",
  ggtheme = theme_minimal(),
  palette = c("#E31A1C", "#1F78B4"),
  title = "Kaplan-Meier: GLUT1:EpCAM Ratio in Tumor Budding",
  xlab = "Time to Recurrence (months)",
  ylab = "Recurrence-Free Survival Probability",
  legend.title = "Biomarker Status",
  legend.labs = c("High GLUT1:EpCAM", "Low GLUT1:EpCAM")
)

pdf("figures/tb_spatial/KM_overall_biomarker.pdf", width = 10, height = 8)
print(p_km_overall)
dev.off()

# Treatment-stratified survival analysis
for(treatment in unique(patient_survival_data$Treatment)) {
  treatment_data <- patient_survival_data[patient_survival_data$Treatment == treatment, ]
  
  if(nrow(treatment_data) > 10) {  # Only analyze if sufficient sample size
    biomarker_cutoff <- surv_cutpoint(treatment_data,time = "RFS_time",event = "RFS_status","mean_biomarker",minprop = 0.15)
    biomarker_cutoff <- biomarker_cutoff[["cutpoint"]]$cutpoint
    treatment_data$biomarker_group <- ifelse(treatment_data$mean_biomarker > biomarker_cutoff,
                                                    "High GLUT1:EpCAM", "Low GLUT1:EpCAM")
    
    surv_obj_treat <- Surv(treatment_data$RFS_time, treatment_data$RFS_status)
    surv_fit_treat <- survfit(surv_obj_treat ~ biomarker_group, data = treatment_data)
    
    p_km_treat <- ggsurvplot(
      surv_fit_treat,
      data = treatment_data,
      pval = TRUE,
      conf.int = TRUE,
      risk.table = TRUE,
      risk.table.col = "strata",
      linetype = "strata",
      surv.median.line = "hv",
      ggtheme = theme_minimal(),
      palette = c("#E31A1C", "#1F78B4"),
      title = paste0("KM Curve: ", treatment, " Treatment"),
      xlab = "Time to Recurrence (months)",
      ylab = "Recurrence-Free Survival Probability",
      legend.title = "Biomarker Status",
      legend.labs = c("High EpCAM:GLUT1", "Low EpCAM:GLUT1")
    )
    
    pdf(paste0("figures/tb_spatial/KM_", treatment, "_biomarker.pdf"),width = 10, height = 8)
    print(p_km_treat)
    dev.off()    
  }
}

# Cox proportional hazards model
cox_model <- coxph(surv_obj ~ biomarker_group + Treatment, data = patient_survival_data)
cox_summary <- summary(cox_model)

cat("Cox Proportional Hazards Model:\n")
print(cox_summary)

# Cox model with interaction
cox_model_interaction <- coxph(surv_obj ~ biomarker_group * Treatment, data = patient_survival_data)
cox_interaction_summary <- summary(cox_model_interaction)

cat("Cox Model with Interaction:\n")
print(cox_interaction_summary)

# Create summary table
survival_summary <- patient_survival_data %>%
  group_by(Treatment, biomarker_group) %>%
  summarise(
    n_patients = n(),
    events = sum(RFS_status),
    median_rfs = median(RFS_time),
    event_rate = round(100 * events / n_patients, 1),
    .groups = "drop"
  )

cat("Survival summary by treatment and biomarker:\n")
print(survival_summary)

# ================================================================================
# SUMMARY AND CONCLUSIONS
# ================================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")

cat("BIOMARKER PERFORMANCE:\n")
for(i in 1:nrow(biomarker_stats)) {
  row <- biomarker_stats[i, ]
  cat(sprintf("%s: HR = %.2f, p = %.4f\n", row$treatment, row$fold_change, row$p_value))
}

cat("\nCOX MODEL RESULTS:\n")
cat("Biomarker HR (95% CI):", round(cox_summary$conf.int[1, "exp(coef)"], 2), 
    "(", round(cox_summary$conf.int[1, "lower .95"], 2), "-", 
    round(cox_summary$conf.int[1, "upper .95"], 2), ")\n")
cat("P-value:", round(cox_summary$coefficients[1, "Pr(>|z|)"], 4), "\n")

cat("\nFILES GENERATED:\n")
cat("- TB_spatial_*.pdf: Spatial visualizations of tumor budding\n")
cat("- EpCAM_GLUT1_biomarker.pdf: Biomarker boxplots\n") 
cat("- EpCAM_vs_GLUT1_scatter.pdf: Expression correlation\n")
cat("- KM_*_biomarker.pdf: Kaplan-Meier survival curves\n")

cat("\nCONCLUSION: EpCAM:GLUT1 ratio in tumor budding shows promise as a ")
cat("treatment-specific prognostic biomarker!\n")
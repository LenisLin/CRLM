# For Tumor cell analysis

# ------------------------------
# Set Working Directory and Source Utilities
# ------------------------------
library(SpatialExperiment)
library(imcRtools)
library(igraph)
library(Hmisc)

library(ggplot2)
library(ggridges)
library(ggalluvial)
library(ggcorrplot)
library(ggpubr)
library(ggraph)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(pheatmap)
library(corrplot)

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
date_time <- "0730"
spe <- readRDS(file.path(saveDir,paste0("spatial_spe_",date_time,".rds")))
img_id_ <- "sample_id"

# head(colData(spe))

# Defined color for CN
RFS_color <- setNames(c("#EFC000FF", "#0073C2FF"),unique(spe$RFS_status))
Tissue_color <- setNames(metadata(spe)$color_vectors$color_20[1:3],unique(spe$Tissue)[1:3])

# Step 1: Extract cells from TC and IM tissues
spe <- spe[,spe$Tissue %in% c("IM","TC")]
tc_im_cells <- as.data.frame(colData(spe))

cat("Total cells in TC and IM:", nrow(tc_im_cells), "\n")

# Step 2: Extract malignant cells (sub_celltype starting with "EC_")
malignant_cells <- tc_im_cells[grepl("^EC_", tc_im_cells$sub_celltype), ]
cat("Malignant cells:", nrow(malignant_cells), "\n")

# Check malignant cell subtypes
cat("Malignant cell subtypes:\n")
print(table(malignant_cells$sub_celltype))

# ================================================================================
# Basic comparison of malignant cell in all cells
# ================================================================================
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
  scale_fill_manual(values = Tissue_color)

ggsave(file.path(figureDir, paste0("all_malignant_cell_barplot.pdf")), p1, width = 14, height = 8)

# Step 5: Compare fractions between RFS_status groups
# Prepare data for statistical comparison
comparison_data <- roi_summary %>%
  select(sample_id, patient_id, Tissue, Treatment, RFS_status, RFS_group, malignant_fraction)
comparison_data$RFS_status <- as.factor(comparison_data$RFS_status)

# Boxplot with wilcox test
p2 <- ggplot(comparison_data, aes(x = RFS_group, y = malignant_fraction, fill = RFS_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~Tissue) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_minimal() +
  labs(title = "Malignant Cell Fraction by Recurrence Status",
       x = "Recurrence Status",
       y = "Malignant Cell Fraction",
       fill = "RFS Status") +
  scale_fill_manual(values = RFS_color)

ggsave(file.path(figureDir, paste0("all_malignant_boxplot_rfs.pdf")), p2, width = 10, height = 6)

rm(p1,p2,comparison_data)
gc()

# Step 7: Paired Boxplot for fraction changes from TC -> IM
# First, get paired TC-IM data for each patient
paired_data <- roi_summary %>%
  group_by(patient_id, Tissue, Treatment, RFS_status, RFS_group) %>%
  summarise(mean_malignant_fraction = mean(malignant_fraction, na.rm = TRUE), .groups = "drop") %>%
  select(patient_id, Tissue, Treatment, RFS_status, RFS_group, mean_malignant_fraction) %>%
  filter(!is.na(mean_malignant_fraction))

# Prepare data for the three-panel plot
plot_data <- paired_data %>%
  mutate(
    # Create a panel variable for faceting
    Panel = case_when(
      Treatment == "Chemo" ~ "Chemo Only",
      Treatment == "Combo" ~ "Combo Only",
      TRUE ~ "Error"
    )
  ) %>%
  # Add "All Samples" version
  bind_rows(
    paired_data %>% mutate(Panel = "All Samples")
  ) %>%
  # Ensure proper factor ordering
  mutate(
    Panel = factor(Panel, levels = c("All Samples", "Chemo Only", "Combo Only")),
    Tissue = factor(Tissue, levels = c("IM", "TC")),
    RFS_group = factor(RFS_group)
  )

# Create Paired-wise boxplot
plot_data$RFS_status <- as.factor(plot_data$RFS_status)
p <- ggplot(plot_data, aes(x = Tissue, y = mean_malignant_fraction, fill = RFS_status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(aes(color = RFS_status) ,size = 1, alpha = 0.8,
             position = position_jitter(width = 0.1, seed = 123)) +
  
  facet_wrap(~Panel, scales = "free_x") +
    labs(
    title = "Paired Analysis: Malignant Fraction in IM vs TC",
    subtitle = "Connected lines show paired samples from the same patient",
    x = "Tissue Type",
    y = "Mean Malignant Fraction",
    fill = "RFS Group",
    color = "RFS Group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(values = RFS_color) +
  scale_color_manual(values = RFS_color) +
  stat_compare_means(method = "wilcox.test", label = "p")
  
ggsave(file.path(figureDir, paste0("malignant fraction comparison between treatments.pdf")), p, width = 10, height = 6)

# ================================================================================
## Tumor Budding Defining
# ================================================================================
# STEP 1: EXTRACT MALIGNANT CELLS FROM IM (Already provided)
cat("=== EXTRACTING MALIGNANT CELLS FROM IM ===\n")

spe_IM <- spe[, spe$Tissue %in% "IM"]
spe_malignant_IM <- spe_IM[, grepl("^EC_", spe_IM$sub_celltype)]
all_IM_cells <- as.data.frame(colData(spe_IM))

cat("Total malignant cells in IM:", ncol(spe_malignant_IM), "\n")
cat("Number of ROIs:", length(unique(spe_malignant_IM$sample_id)), "\n")
cat("Number of patients:", length(unique(spe_malignant_IM$patient_id)), "\n")

# STEP 2: FILTER LONG EDGES IN DELAUNAY TRIANGULATION
cat("\n=== FILTERING LONG EDGES ===\n")

# Get Delaunay triangulation edges
delaunay_edges <- as.data.frame(colPair(spe_malignant_IM, "delaunay")) 

# Get spatial coordinates (assuming they exist in colData)
cell_coords <- data.frame(
  cell_id = 1:ncol(spe_malignant_IM),
  x = spatialCoords(spe_malignant_IM)[,1],
  y = spatialCoords(spe_malignant_IM)[,2],
  sample_id = spe_malignant_IM$sample_id
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

# STEP 3: DETECT CONNECTED COMPONENTS AND FILTER LARGE CLUSTERS
cat("\n=== DETECTING CONNECTED COMPONENTS ===\n")

# Create adjacency list for each ROI separately
roi_list <- unique(spe_malignant_IM$sample_id)
budding_results <- list()

cluster_size_threshold_max <- 20
cluster_size_threshold_min <- 3

# Extended analysis loop
for(roi in roi_list) {
  cat("Processing ROI:", roi, "\n")
  
  # Get cells in this ROI
  roi_cells <- which(spe_malignant_IM$sample_id == roi)
  all_roi_cells <- which(spe_IM$sample_id == roi)
  
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
  large_components <- roi_components[component_sizes > cluster_size_threshold_max * 2]
  
  cat("  Total components:", length(roi_components), "\n")
  cat("  Small components (budding):", length(budding_components), "\n")
  cat("  Large components (core tumor):", length(large_components), "\n")
  
  # Distance calculation section ===
  budding_distances <- list()
  
  # Calculate distances from each small cluster to the top-2 largest clusters
  if(length(budding_components) > 0 && length(large_components) > 0) {
      # Get sizes of all large clusters and sort to find top-2 largest
      large_sizes <- sapply(large_components, length)
      sorted_indices <- order(large_sizes, decreasing = TRUE)
      
      # Determine how many large clusters to use (max 2)
      n_top_clusters <- min(2, length(large_components))
      top_large_indices <- sorted_indices[1:n_top_clusters]
      
      cat("  Large cluster sizes:", large_sizes[sorted_indices], "\n")
      cat("  Using top", n_top_clusters, "largest clusters for distance calculation\n")
      
      # Store distance information for each small cluster
      budding_distances_to_top2 <- list()
      
      for(i in seq_along(budding_components)) {
        cluster_name <- names(budding_components)[i]
        small_cluster <- budding_components[[i]]
        small_cluster_coords <- spatialCoords(spe_malignant_IM)[small_cluster, , drop = FALSE]
        
        # Calculate distances to each of the top large clusters
        distances_to_top_clusters <- list()
        distance_values <- numeric(n_top_clusters)
        
        for(j in 1:n_top_clusters) {
          large_idx <- top_large_indices[j]
          large_cluster_name <- names(large_components)[large_idx]
          large_cluster <- large_components[[large_idx]]
          large_cluster_coords <- spatialCoords(spe_malignant_IM)[large_cluster, , drop = FALSE]
          
          # Calculate distances using your existing function
          distances <- calculate_cluster_distances(small_cluster_coords, large_cluster_coords)
          distances_to_top_clusters[[large_cluster_name]] <- distances
          
          # Store the mean distance for averaging
          distance_values[j] <- distances$min
        }
        
        # Calculate distance
        mean_distance_to_top2 <- min(distance_values) # mean

        # Store comprehensive distance information
        budding_distances_to_top2[[cluster_name]] <- list(
          mean_distance_to_top2 = mean_distance_to_top2,
          top_cluster_sizes = large_sizes[top_large_indices]
        )
        
        cat("    ", cluster_name, "- Size:", length(small_cluster), 
            "- Mean dist to top-2:", round(mean_distance_to_top2, 2), "\n")
      }
  }
  
  ## Remove Bile duct
  bile_duct_compoments <- list()
  for(tb_name_ in names(budding_components)){
    if(budding_distances_to_top2[[tb_name_]]$mean_distance_to_top2>100){
      bile_duct_compoments[[tb_name_]] <- budding_components[[tb_name_]]
      budding_components[[tb_name_]] <- NULL
    }
  }
  
  # Store results (extended with distance information)
  budding_results[[roi]] <- list(
    budding_components = budding_components,
    budding_cells = rownames(colData(spe_malignant_IM))[unname(unlist(budding_components))],
    bile_duct_compoments = bile_duct_compoments,
    bile_duct_cells = rownames(colData(spe_malignant_IM))[unname(unlist(bile_duct_compoments))],
    large_components = large_components,
    large_components_cells = rownames(colData(spe_malignant_IM))[unname(unlist(large_components))],
    component_sizes = component_sizes,
    roi_edges_from = rownames(colData(spe_malignant_IM))[roi_edges_from],
    roi_edges_to = rownames(colData(spe_malignant_IM))[roi_edges_to]
  )
}

# Collect Results
TB_indexs <- c()
BD_indexs <- c()
TC_indexs <- c()

for(i in 1:length(budding_results)){
  TB_indexs <- c(TB_indexs,budding_results[[i]][["budding_cells"]])
  BD_indexs <- c(BD_indexs,budding_results[[i]][["bile_duct_cells"]])
  TC_indexs <- c(TC_indexs,budding_results[[i]][["large_components_cells"]])
}

cat("Length of Tumor Budding cells", length(TB_indexs), "\n")
cat("Length of Bile Duct cells:", length(BD_indexs), "\n")
cat("Length of Core Tumor cells", length(TC_indexs), "\n")

spe_malignant_IM$cell_position <- "Single"
spe_malignant_IM$cell_position[match(TB_indexs,colnames(spe_malignant_IM))] <- "TumorBudding"
spe_malignant_IM$cell_position[match(BD_indexs,colnames(spe_malignant_IM))] <- "BileDuct"
spe_malignant_IM$cell_position[match(TC_indexs,colnames(spe_malignant_IM))] <- "TumorCore"

saveRDS(spe_malignant_IM,file = file.path(figureDir,"IM_malignant_spe.rds"))
saveRDS(budding_results,file.path(figureDir,"TumorBudding_results.rds"))

# ================================================================================
# Compare the number of Tumor budding number and size
# ================================================================================
# Get all ROI names
roi_names <- names(budding_results)

# Initialize result data frame
results <- data.frame(
  ROI_ID = character(),
  budding_number = numeric(),
  total_size = numeric(),
  mean_size = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each ROI
for (roi in roi_names) {
  # Get budding components for this ROI
  roi_components <- budding_results[[roi]][["budding_components"]]
  
  # Calculate metrics for this ROI
  roi_metrics <- calculate_roi_budding_metrics(roi_components)
  
  # Add ROI ID and combine with results
  roi_metrics$ROI_ID <- roi
  results <- rbind(results, roi_metrics[c("ROI_ID", "budding_number", "total_size", "mean_size")])
}

results$patient_id <- spe$patient_id[match(results$ROI_ID,spe$sample_id)]
results$Treatment <- spe$Treatment[match(results$ROI_ID,spe$sample_id)]
results$RFS_status <- spe$RFS_status[match(results$ROI_ID,spe$sample_id)]

plotdf <- results %>% 
  group_by(patient_id, Treatment, RFS_status) %>%
  summarise(across(c(1:(ncol(results) - 3)), mean, na.rm = TRUE)) %>%
  select(-ROI_ID)

plotdf <- pivot_longer(plotdf,cols = 4:ncol(plotdf), values_to = "value", names_to = "Budding_features")

head(plotdf)

p <- ggplot(plotdf, aes(x = factor(RFS_status), y = value, fill = factor(RFS_status))) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
  geom_point(position = position_jitter(width = 0.2, seed = 123), 
             alpha = 0.6, size = 1) +
  facet_grid(Budding_features ~ Treatment, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     size = 3,
                     vjust = 0.5) +
  scale_fill_manual(values = c("0" = "#0073c2b2", "1" = "#efc000b2"),
                    labels = c("0" = "No Early Relapse", "1" = "Early Relapse")) +
  scale_x_discrete(labels = c("0" = "No\nRelapse", "1" = "Early\nRelapse")) +
  labs(
    title = paste0("Tumor Budding Size analysis"),
    x = "RFS Status",
    y = "Value",
    fill = "RFS Status"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    strip.text = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(figureDir, paste0("Boxplot_of_Tumor_Budding_charateristic_across_treatment.pdf")), p, width = 6, height = 8)
rm(p, plotdf, results,roi_components,roi_metrics)
gc()

# ================================================================================
# SPATIAL VISUALIZATION OF TUMOR BUDDING
# ================================================================================
cat("\n=== CREATING SPATIAL VISUALIZATIONS ===\n")

# Add EpCAM expression levels and match TB label
spe_IM$EpCAM_status <- "EpCAM-"
spe_IM$EpCAM_status[match(c(TB_indexs,TC_indexs),colnames(spe_IM))] <- "EpCAM+"

# Add spatial coordinates
budding_cell_data <- colData(spe_IM)

budding_cell_data$x <- spatialCoords(spe_IM)[, 1]
budding_cell_data$y <- spatialCoords(spe_IM)[, 2]

# Create spatial plots for selected ROIs (show first 6 ROIs)
selected_rois <- unique(spe_IM$sample_id)
all_coord <- spatialCoords(spe_malignant_IM)

figureDir_temp <- file.path(figureDir,"TB_spatialplot")
if(!file.exists(figureDir_temp)){
  dir.create(figureDir_temp,recursive = T)
}

for(roi in selected_rois) {
  cat("Creating spatial plot for ROI:", roi, "\n")
  
  roi_data <- budding_cell_data[budding_cell_data$sample_id == roi, ]
  
  # Get edges for this ROI
  roi_edges_from <- budding_results[[roi]]$roi_edges_from
  roi_edges_to <- budding_results[[roi]]$roi_edges_to
  
  # Filter edges to only include budding cells
  all_cell_indices <- rownames(roi_data)
  edge_mask <- (roi_edges_from %in% all_cell_indices) & (roi_edges_to %in% all_cell_indices)
  roi_edges_from <- roi_edges_from[edge_mask]
  roi_edges_to <- roi_edges_to[edge_mask]
  
  # Create edge data frame for plotting
  from_cell <- roi_edges_from
  to_cell <- roi_edges_to
  
  from_coord <- all_coord[match(from_cell,colnames(spe_malignant_IM)) , ]
  to_coord <- all_coord[match(to_cell,colnames(spe_malignant_IM)) , ]
  
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
  ggsave(file.path(figureDir_temp,paste0("TB_spatial_", gsub("[^A-Za-z0-9]", "_", roi), ".pdf")) , 
         p, width = 8, height = 8, bg = "black")
}

rm(budding_cell_data,budding_components,budding_distances,budding_distances_to_top2)
rm(bile_duct_compoments,cell_coords,components,delaunay_edges,distances,distances_to_top_clusters,large_components,large_cluster_coords,large_cluster)
rm(roi_summary,roi_data,roi_graph,tc_im_cells,edge_data,edge_from,edge_to,i,j,edge_lengths,edge_mask,distance_values,filtered_edge_lengths,filtered_edges_from,filtered_edges_to)
rm(to_coord,p,small_cluster_coords,from_coord,all_coord)
gc()

# ================================================================================
# Malignant subtype analysis - Sankey diagram for fraction changes from TC -> IM -> TB
# ================================================================================
## Collect Data
malignant_cell_meta <- as.data.frame(colData(spe))
TB_cells <- colnames(spe_malignant_IM)[spe_malignant_IM$cell_position == "TumorBudding"]
BD_cells <- colnames(spe_malignant_IM)[spe_malignant_IM$cell_position == "BileDuct"]

malignant_cell_meta$Tissue[match(TB_cells,rownames(malignant_cell_meta))] <- "TB"
malignant_cell_meta$Tissue[match(BD_cells,rownames(malignant_cell_meta))] <- "BD"

malignant_cell_meta <- malignant_cell_meta[grepl("^EC_", malignant_cell_meta$sub_celltype),]

roi_fractions <- malignant_cell_meta %>%
  select(sample_id, patient_id, Tissue, Treatment, RFS_status, sub_celltype) %>%
  group_by(sample_id, Tissue) %>%
  mutate(
    total_cells_in_tissue = n()  # Total cells in this sample_id + tissue combination
  ) %>%
  group_by(sample_id, Tissue, sub_celltype, .add = TRUE) %>%
  summarise(
    # Preserve metadata (take first value since they should be identical)
    patient_id = first(patient_id),
    Treatment = first(Treatment),
    RFS_status = first(RFS_status),
    # Calculate counts and fractions
    subcelltype_count = n(),                               # Count of this sub_celltype
    total_cells_in_tissue = first(total_cells_in_tissue),  # Total cells in this sample + tissue
    fraction = subcelltype_count / total_cells_in_tissue,  # Fraction within this tissue
    .groups = "drop"
  ) %>%
  arrange(sample_id, Tissue, sub_celltype)

roi_fractions <- roi_fractions[roi_fractions$Tissue != "BD",] ## remove bile duct

# Create bins for fraction levels to make Sankey readable
create_fraction_bins <- function(x) {
  case_when(
    x < 0.1 ~ "Low (<10%)",
    x < 0.3 ~ "Medium (10-30%)",
    x < 0.5 ~ "High (30-50%)",
    TRUE ~ "Very High (>50%)"
  )
}

for(target_celltype in unique(roi_fractions$sub_celltype)){
  # Prepare data for three-tissue Sankey plot
  sankey_prep_data <- roi_fractions %>%
    # 1. Select the chosen celltype
    filter(sub_celltype == target_celltype) %>%
    # 2. Select the columns we need
    select(patient_id, sample_id, Tissue, fraction, Treatment, RFS_status) %>%
    # 3. Calculate the mean fraction of specific tissue among same patient 
    # (one patient may have multiple samples for the same tissue)
    group_by(patient_id, Tissue, Treatment, RFS_status) %>%
    summarise(mean_fraction = mean(fraction, na.rm = TRUE), .groups = "drop") %>%
    # 4. Reshape to wide format to get one row per patient with tissue columns
    pivot_wider(
      names_from = Tissue,
      values_from = mean_fraction,
      names_prefix = "fraction_"
    ) %>%
    # Keep only patients with all three tissue types
    filter(!is.na(fraction_TC) & !is.na(fraction_IM) & !is.na(fraction_TB)) %>%
    # Calculate overall mean fraction across all tissues for each patient
    mutate(mean_fraction_across_tissues = rowMeans(select(., starts_with("fraction_")), na.rm = TRUE)) %>%
    # Create fraction bins AFTER calculating mean fractions
    mutate(
      fraction_bin_TC = create_fraction_bins(fraction_TC),
      fraction_bin_IM = create_fraction_bins(fraction_IM),
      fraction_bin_TB = create_fraction_bins(fraction_TB),
      mean_fraction_bin = create_fraction_bins(mean_fraction_across_tissues)
    ) %>%
    # Add RFS grouping
    mutate(RFS_group = case_when(
      RFS_status == 0 ~ "No Early Relapse",
      RFS_status == 1 ~ "Early Relapse",
      TRUE ~ "Unknown"
    ))
  
  # Display the structure of the prepared data
  cat("Columns in sankey_prep_data:\n")
  print(colnames(sankey_prep_data))
  
  # Create the three-tissue Sankey plot
  sankey_data <- sankey_prep_data %>%
    count(RFS_group, fraction_bin_TC, fraction_bin_IM, fraction_bin_TB, Treatment) %>%
    rename(freq = n) %>%
    filter(freq > 0)  # Remove empty combinations
  
  p_sankey_rfs <- ggplot(sankey_data, 
                         aes(axis1 = fraction_bin_TC, 
                             axis2 = fraction_bin_IM, 
                             axis3 = fraction_bin_TB, 
                             y = freq)) +
    geom_alluvium(aes(fill = RFS_group), alpha = 0.7, width = 0.3) +
    geom_stratum(width = 0.3, alpha = 0.8, color = "white", size = 0.5) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
              size = 2.5, angle = 0, color = "black", fontface = "bold") +
    scale_x_discrete(limits = c("TC", "IM", "TB"), 
                     expand = c(0.15, 0.15)) +
    facet_wrap(Treatment~RFS_group, scales = "free_y") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      panel.grid = element_blank()
    ) +
    labs(title = paste0("Tissue Progression Flow: ", target_celltype, " Fraction Changes"),
         subtitle = "TC → IM → TB transition patterns by recurrence status",
         x = "Tissue Region",
         y = "Number of Patients",
         fill = "RFS Status") +
    scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))
  
    ggsave(file.path(figureDir, paste0(target_celltype,"_change_acrosse_tissue_sankey.pdf")), p_sankey_rfs, width = 12, height = 6)
}

# Stream flow of cell types change
theme_data_rfs <- roi_fractions %>%
  mutate(RFS_group = ifelse(RFS_status == 1, "Early Relapse", "No Early Relapse")) %>%
  group_by(Tissue, sub_celltype, RFS_group) %>%
  summarise(mean_fraction = mean(fraction, na.rm = TRUE), .groups = "drop") %>%
  # Normalize fractions to sum to 1.0 within each tissue and RFS group
  group_by(Tissue, RFS_group) %>%
  mutate(
    total_fraction = sum(mean_fraction),
    normalized_fraction = mean_fraction / total_fraction
  ) %>%
  ungroup() %>%
  mutate(Tissue = factor(Tissue, levels = c("TC", "IM", "TB"))) %>%
  mutate(tissue_pos = as.numeric(Tissue)) %>%
  arrange(RFS_group, Tissue, sub_celltype) %>%
  group_by(RFS_group, Tissue) %>%
  mutate(
    cumulative_fraction = cumsum(normalized_fraction),
    cumulative_start = lag(cumulative_fraction, default = 0)
  ) %>%
  ungroup()

ribbon_data_rfs <- theme_data_rfs %>%
  select(tissue_pos, sub_celltype, cumulative_start, cumulative_fraction, RFS_group) %>%
  rename(ymin = cumulative_start, ymax = cumulative_fraction)


p_rfs <- ggplot(ribbon_data_rfs, aes(x = tissue_pos)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = sub_celltype), 
              alpha = 0.8, color = "white", size = 0.3) +
  facet_wrap(~RFS_group) +
  scale_x_continuous(
    breaks = 1:3,
    labels = c("TC", "IM", "TB"),
    expand = c(0.05, 0.05)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = c(0, 0)
  ) +
  labs(
    title = "Epithelial Subtype Relative Composition by RFS Group",
    subtitle = "ThemeRiver comparison showing normalized proportions between groups",
    x = "Tissue Region",
    y = "Cumulative Proportion",
    fill = "Epithelial Subtype"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  scale_fill_viridis_d(option = "plasma")

ggsave(file.path(figureDir, paste0("All_malignant_change_acrosse_tissue.pdf")), p_rfs, width = 12, height = 8)

rm(p_rfs,ribbon_data_rfs,theme_data_rfs,sankey_data,p_sankey)
gc()

# ================================================================================
# Malignant subtype analysis - Boxplot for fraction changes
# ================================================================================

# Prepare data and perform tests for each treatment group
all_data <- roi_fractions %>%
  mutate(RFS_group = ifelse(RFS_status == 1, "Early Relapse", "No Early Relapse"))

chemo_data <- all_data %>% filter(Treatment == "Chemo")
combo_data <- all_data %>% filter(Treatment == "Combo")

# Perform statistical tests
all_tests <- perform_wilcox_tests(all_data, "All Samples")
chemo_tests <- perform_wilcox_tests(chemo_data, "Chemo Only")
combo_tests <- perform_wilcox_tests(combo_data, "Combo Only")

# Combine test results
combined_tests <- bind_rows(all_tests, chemo_tests, combo_tests)

# Display test results
cat("Statistical Test Results (BH-adjusted):\n")
print(combined_tests %>% 
        select(group, Tissue, sub_celltype, n0, n1, p_value, p_adjusted, significance) %>%
        arrange(group, Tissue, sub_celltype))

# Create the three comparison plots
p1 <- create_comparison_plot(all_data, combined_tests, "All Samples", "All Samples")
ggsave(file.path(figureDir, paste0("Boxplot_of_all_malignant_subtypes_acrosee_tissues.pdf")), p1, width = 12, height = 10)

p2 <- create_comparison_plot(chemo_data, combined_tests, "Chemotherapy Only", "Chemo Only")
ggsave(file.path(figureDir, paste0("Boxplot_of_Chemo_malignant_subtypes_acrosee_tissues.pdf")), p2, width = 12, height = 10)

p3 <- create_comparison_plot(combo_data, combined_tests, "Combination Therapy Only", "Combo Only")
ggsave(file.path(figureDir, paste0("Boxplot_of_Combo_malignant_subtypes_acrosee_tissues.pdf")), p3, width = 12, height = 10)

rm(p1,p2,p3,all_data,chemo_data,combo_data,all_tests,chemo_tests,combo_tests)
gc()

# ================================================================================
# MALIGNANT MARKER ANALYSIS
# ================================================================================
head(malignant_cell_meta)
table(malignant_cell_meta$Tissue)
table(malignant_cell_meta$sub_celltype)

malignant_cells <- malignant_cell_meta[malignant_cell_meta$Tissue != "BD",]

# Compare the expression among all EP cells during IM
malignant_roi_expression <- data.frame()
malignant_cell_expression <- assay(spe)[,rownames(malignant_cells)]

malignant_cells$RFS_group <- ifelse(malignant_cells$RFS_status == 0, 
                                "No Early Relapse", "Early Relapse")

tumor_markers <- c("EpCAM", "Vimentin", "Ki67", "VEGF", "CA_IX", "GLUT1", "HK2", "FASN","CD279")
tissues <- unique(malignant_cells$Tissue)

for(tissue_ in tissues){
  roi_malignant_cells <- malignant_cells[malignant_cells$Tissue == tissue_,]
  
for(roi in unique(roi_malignant_cells$sample_id)) {
  roi_tissue_malignant_cells <- roi_malignant_cells[roi_malignant_cells$sample_id == roi, ]
    if(nrow(roi_tissue_malignant_cells) == 0) next
    
    malignant_cell_indices <- rownames(roi_tissue_malignant_cells)
    roi_expression <- malignant_cell_expression[tumor_markers, malignant_cell_indices, drop = FALSE]
    
    # Calculate mean expression for this ROI
    roi_means <- rowMeans(roi_expression)
    
    # Get clinical metadata
    roi_meta <- roi_tissue_malignant_cells[1, c("sample_id", "patient_id", "Treatment", "Tissue", "RFS_status", "RFS_group")]
    
    # Combine with expression data
    roi_data <- data.frame(t(roi_means), roi_meta)
    malignant_roi_expression <- rbind(malignant_roi_expression, roi_data)
  }
}

rownames(malignant_roi_expression) <- paste0(rownames(malignant_roi_expression),"_",malignant_roi_expression$Tissue)

## volcano plot
# Create combinations for analysis
tissue_types <- c("TC+IM+TB", "TC", "IM","TB")  # Adjust based on your actual tissue types
treatments <- c("All", "Chemo", "Combo")

# Initialize results dataframe
results_list <- list()

# Loop through each combination
for(tissue in tissue_types) {
  for(treatment in treatments) {
    
    # Filter data based on tissue and treatment
    if(tissue == "TC+IM+TB") {
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
    for(marker in tumor_markers) {
      
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
      
      # Perform wilcox test
      t_test_result <- wilcox.test(group1_data, group2_data)
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
results_df$tissue_type <- factor(results_df$tissue_type, levels = c("TC+IM+TB", "TC", "IM", "TB"))
results_df$treatment <- factor(results_df$treatment, levels = c("All", "Chemo", "Combo"))

# Print results summary
print("Summary of significant differences:")
significant_results <- results_df[results_df$significance != "ns", ]
print(significant_results[, c("tissue_type", "treatment", "marker", "fold_change", "adj_p_value", "significance")])

# --- Create the Volcano Plot ---
log2_fc_threshold <- log2(1.5)
significant_results <- results_df[results_df$significance != "ns", ]
significant_results <- significant_results[significant_results$fold_change >= log2_fc_threshold | significant_results$fold_change <= (-log2_fc_threshold),]
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
  facet_grid(treatment ~ tissue_type, scales = "free") +
  
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
  ) +
  labs(subtitle = "Right: Early Recurrence Left: Later Recurrence",
       plot.subtitle = element_text(hjust = 0.5, size = 11))

# Print volcano plot
print(p_volcano)

# Save plots
ggsave(file.path(figureDir,"marker_volcano_plot.pdf"), plot = p_volcano, width = 10, height = 7.5)
write.csv(results_df, file.path(figureDir,"marker_expression_analysis.csv"), row.names = FALSE)

# ================================================================================
# The correlation between marker expression
# ================================================================================
print("=== MARKER EXPRESSION CORRELATION ANALYSIS ACROSS TISSUES ===")

# Define marker columns (expression data)
marker_columns <- c("EpCAM", "Vimentin", "Ki67", "VEGF", "CA_IX", "GLUT1", "HK2", "FASN", "CD279")

print("Marker columns for analysis:")
print(marker_columns)

# Check tissue distribution
tissue_summary <- table(malignant_roi_expression$Tissue)
print("\nTissue distribution:")
print(tissue_summary)

# Calculate correlations for each tissue type
tissue_correlations <- list()
tissue_types <- unique(malignant_roi_expression$Tissue)

for (tissue in tissue_types) {
  print(paste("Calculating correlations for", tissue, "tissue..."))
  
  tissue_data <- malignant_roi_expression[malignant_roi_expression$Tissue == tissue, ]
  
  cor_result <- calculate_correlation_with_pval(tissue_data, method = "spearman")
  
  if (!is.null(cor_result)) {
    tissue_correlations[[tissue]] <- cor_result
    print(paste("  - Samples used:", cor_result$n_samples))
  } else {
    print(paste("  - Insufficient data for", tissue))
  }
}

# Create correlation heatmaps with significance for each tissue
plot_list <- list()

for (tissue in names(tissue_correlations)) {
  cor_data <- tissue_correlations[[tissue]]
  
  # Prepare data for ggcorrplot
  cor_matrix <- cor_data$correlations
  p_matrix <- cor_data$p_values
  
  # Get significance stars
  sig_stars <- get_significance_stars(p_matrix)
  
  # Create correlation plot with ggcorrplot
  p <- ggcorrplot(cor_matrix, 
                  method = "circle",
                  hc.order = FALSE,
                  # p.mat = p_matrix,
                  type = "lower",
                  lab = TRUE,
                  lab_size = 3,
                  show.diag = FALSE,
                  colors = c("#7ca6dc", "white", "#cd534c"),
                  title = paste("Spearman Correlations -", tissue, "Tissue ( n=",cor_data$n_samples,")"),
                  ggtheme = theme_bw()) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      legend.position = "bottom"
    )
  
  plot_list[[tissue]] <- p
}

# Create a combined plot
if (length(plot_list) > 1) {
  combined_plot <- wrap_plots(plot_list, ncol =3)
  print(combined_plot)
  
  ggsave(file.path(figureDir,"marker expression corretation within different tumor tissue.pdf"), 
         plot = combined_plot, width = 12, height = 6)
}

# Create detailed correlation table with significance
correlation_summary <- data.frame()

for (tissue in names(tissue_correlations)) {
  cor_data <- tissue_correlations[[tissue]]
  cor_matrix <- cor_data$correlations
  p_matrix <- cor_data$p_values
  
  # Convert matrices to long format for detailed analysis
  for (i in 1:(nrow(cor_matrix)-1)) {
    for (j in (i+1):ncol(cor_matrix)) {
      marker1 <- rownames(cor_matrix)[i]
      marker2 <- colnames(cor_matrix)[j]
      correlation <- cor_matrix[i, j]
      p_value <- p_matrix[i, j]
      
      correlation_summary <- rbind(correlation_summary, data.frame(
        Tissue = tissue,
        Marker1 = marker1,
        Marker2 = marker2,
        Correlation = round(correlation, 3),
        P_value = p_value,
        Significance = case_when(
          p_value <= 0.001 ~ "***",
          p_value <= 0.01 ~ "**", 
          p_value <= 0.05 ~ "*",
          TRUE ~ "NS"
        ),
        Abs_Correlation = abs(correlation),
        n_samples = cor_data$n_samples,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Apply multiple testing correction within each tissue
correlation_summary <- correlation_summary %>%
  group_by(Tissue) %>%
  mutate(P_adj = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    Significance_raw = case_when(
      P_value <= 0.001 ~ "***",
      P_value <= 0.01 ~ "**",
      P_value <= 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    Significance_adj = case_when(
      P_adj <= 0.001 ~ "***", 
      P_adj <= 0.01 ~ "**",
      P_adj <= 0.05 ~ "*",
      TRUE ~ "NS"
    )
  )

# Sort by absolute correlation strength
correlation_summary <- correlation_summary %>%
  arrange(Tissue, desc(Abs_Correlation))

print("\n=== TOP CORRELATIONS BY TISSUE ===")
for (tissue in tissue_types) {
  tissue_cors <- correlation_summary %>%
    filter(Tissue == tissue, Significance_adj != "NS") %>%
    head(10)
  
  if (nrow(tissue_cors) > 0) {
    print(paste("\nTop significant correlations in", tissue, "tissue:"))
    print(tissue_cors[, c("Marker1", "Marker2", "Correlation", "P_adj", "Significance_adj")])
  } else {
    print(paste("\nNo significant correlations found in", tissue, "tissue"))
  }
}

# GLUT1-Focused Chord Plot Visualization
# Create chord plots showing GLUT1 correlations with other markers across tissues
library(circlize)
library(grid)
library(gridExtra)

print("=== CREATING GLUT1-FOCUSED CHORD PLOTS ===")

# Define marker columns (expression data)
marker_columns <- c("EpCAM", "Vimentin", "Ki67", "VEGF", "CA_IX", "GLUT1", "HK2", "FASN", "CD279")

# Extract GLUT1 correlations for each tissue
tissue_types <- unique(malignant_roi_expression$Tissue)
glut1_correlations <- list()

for (tissue in tissue_types) {
  print(paste("Extracting GLUT1 correlations for", tissue, "tissue..."))
  
  tissue_data <- malignant_roi_expression[malignant_roi_expression$Tissue == tissue, ]
  
  cor_result <- extract_glut1_correlations(tissue_data, tissue)
  
  if (!is.null(cor_result)) {
    glut1_correlations[[tissue]] <- cor_result
    print(paste("  - Samples used:", cor_result$n_samples))
  } else {
    print(paste("  - Insufficient data for", tissue))
  }
}

# Create individual chord plots and save plot data
print("\nCreating chord plots for each tissue...")

# Set up for three-column layout
pdf(file.path(figureDir,"GLUT1_Chord_Plots_Three_Tissues.pdf"), width = 12, height = 6)
par(mfrow = c(1, 3))

plot_data_list <- list()

for (tissue in names(glut1_correlations)) {
  print(paste("Creating chord plot for", tissue, "..."))
  
  # Create chord plot
  plot_data <- create_glut1_chord_plot(glut1_correlations[[tissue]], tissue)
  plot_data_list[[tissue]] <- plot_data
  
  # Clear circos for next plot
  circos.clear()
}

dev.off()

# Save results
write.csv(correlation_summary, "Marker_Correlation_Analysis_by_Tissue.csv", row.names = FALSE)

print("\n=== CORRELATION ANALYSIS COMPLETE ===")

# ================================================================================
# Hypothesis
# GLUT1-high -> chemotherapy -> early recurrence
# FASN-high -> combination therapy -> recurrence
# ================================================================================
malignant_patient_expression <- malignant_roi_expression[,-match(c("sample_id","Tissue"),colnames(malignant_roi_expression))]
malignant_patient_expression <- malignant_patient_expression %>% 
  group_by(patient_id, Treatment, RFS_status, RFS_group) %>%
  summarise(across(c(1:(ncol(malignant_patient_expression) - 4)), mean, na.rm = TRUE))
malignant_patient_expression <- as.data.frame(malignant_patient_expression)
malignant_patient_expression$RFS_time <- spe$RFS_time[match(malignant_patient_expression$patient_id,spe$patient_id)]

# Check if RFS_time exists
if(!"RFS_time" %in% colnames(malignant_patient_expression)) {
  cat("ERROR: RFS_time column is missing from the data.\n")
  cat("Please add the time-to-event variable to perform survival analysis.\n")
  stop("Missing RFS_time variable")
}

# List of genes to analyze
gene_list <- c("EpCAM", "Vimentin", "Ki67", "VEGF", "CA_IX", "GLUT1", "HK2", "FASN")
figureDir_temp <- file.path(figureDir,"KM_malignant_single_gene")
if(!file.exists(figureDir_temp)){
  dir.create(figureDir_temp,recursive = T)
}

# Main analysis loop
for(gene in gene_list) {
  cat("ANALYZING GENE:", gene, "\n")
  cat(rep("=", 60), "\n")
  
  # Calculate thresholds
  mean_cutoff <- mean(malignant_patient_expression[[gene]], na.rm = TRUE)
  median_cutoff <- median(malignant_patient_expression[[gene]], na.rm = TRUE)
  optimal_cutoff <- find_optimal_cutoff(malignant_patient_expression, gene)
  
  cat("Thresholds - Mean:", round(mean_cutoff, 3), 
      "| Median:", round(median_cutoff, 3), 
      "| Optimal:", round(optimal_cutoff, 3), "\n")
  
  # Prepare treatment groups
  all_data <- malignant_patient_expression
  chemo_data <- malignant_patient_expression %>% filter(Treatment == "Chemo")
  combo_data <- malignant_patient_expression %>% filter(Treatment == "Combo")
  
  # Define thresholds
  thresholds <- list(
    "Mean" = mean_cutoff,
    "Median" = median_cutoff,
    "Optimal" = optimal_cutoff
  )
  
  # Loop through each threshold
  for(threshold_name in names(thresholds)) {
    
    cutoff_value <- thresholds[[threshold_name]]
    cat("\n--- ", threshold_name, " threshold (", round(cutoff_value, 3), ") ---\n")
    
    # Create groups
    current_data <- all_data %>%
      filter(!is.na(.data[[gene]])) %>%
      mutate(
        expression_group = ifelse(.data[[gene]] >= cutoff_value, "High", "Low"),
        expression_group = factor(expression_group, levels = c("Low", "High"))
      )

      # Create survival object
      current_data$RFS_time <- as.numeric(current_data$RFS_time)
      current_data$RFS_status <- as.numeric(current_data$RFS_status)
      
      surv_obj <- Surv(time = current_data$RFS_time, event = current_data$RFS_status)
      
      # Fit survival curve
      fit <- survfit(surv_obj ~ expression_group, data = current_data)
      
      # Log-rank test
      logrank_test <- survdiff(surv_obj ~ expression_group, data = current_data)
      p_value <- pchisq(logrank_test$chisq, df = 1, lower.tail = FALSE)
      
      # Sample sizes
      n_low <- sum(current_data$expression_group == "Low")
      n_high <- sum(current_data$expression_group == "High")
      
      # Create plot
      p_all <- ggsurvplot(
        fit,
        data = current_data,
        pval = TRUE,
        pval.method = TRUE,
        conf.int = FALSE,
        risk.table = TRUE,
        risk.table.height = 0.3,
        ggtheme = theme_bw(),
        palette = c("#3B4992FF", "#EE0000FF"),
        legend.title = paste0(gene, " Expression"),
        legend.labs = c(paste0("Low (n=", n_low, ")"), paste0("High (n=", n_high, ")")),
        title = paste0(gene, " - ", threshold_name, " threshold (", round(cutoff_value, 3), ") - All Patients"),
        xlab = "Time (months)",
        ylab = "Relapse-free survival probability"
      )
      cat("All patients: p =", format.pval(p_value), "| Low: n =", n_low, "| High: n =", n_high, "\n")
    
    
    # === CHEMO ONLY ===
    cat("Processing Chemo Only...\n")
    
    # Create groups
    current_data <- chemo_data %>%
      filter(!is.na(.data[[gene]])) %>%
      mutate(
        expression_group = ifelse(.data[[gene]] >= cutoff_value, "High", "Low"),
        expression_group = factor(expression_group, levels = c("Low", "High"))
      )

      # Create survival object
      surv_obj <- Surv(time = current_data$RFS_time, event = current_data$RFS_status)
      
      # Fit survival curve
      fit <- survfit(surv_obj ~ expression_group, data = current_data)
      
      # Log-rank test
      logrank_test <- survdiff(surv_obj ~ expression_group, data = current_data)
      p_value <- pchisq(logrank_test$chisq, df = 1, lower.tail = FALSE)
      
      # Sample sizes
      n_low <- sum(current_data$expression_group == "Low")
      n_high <- sum(current_data$expression_group == "High")
      
      # Create plot
      p_chemo <- ggsurvplot(
        fit,
        data = current_data,
        pval = TRUE,
        pval.method = TRUE,
        conf.int = FALSE,
        risk.table = TRUE,
        risk.table.height = 0.3,
        ggtheme = theme_bw(),
        palette = c("#3B4992FF", "#EE0000FF"),
        legend.title = paste0(gene, " Expression"),
        legend.labs = c(paste0("Low (n=", n_low, ")"), paste0("High (n=", n_high, ")")),
        title = paste0(gene, " - ", threshold_name, " threshold (", round(cutoff_value, 3), ") - Chemotherapy Only"),
        xlab = "Time (months)",
        ylab = "Relapse-free survival probability"
      )
      cat("Chemo only: p =", format.pval(p_value), "| Low: n =", n_low, "| High: n =", n_high, "\n")
    
    # === COMBO ONLY ===
    cat("Processing Combo Only...\n")
    
    # Create groups
    current_data <- combo_data %>%
      filter(!is.na(.data[[gene]])) %>%
      mutate(
        expression_group = ifelse(.data[[gene]] >= cutoff_value, "High", "Low"),
        expression_group = factor(expression_group, levels = c("Low", "High"))
      )
    

      # Create survival object
      surv_obj <- Surv(time = current_data$RFS_time, event = current_data$RFS_status)
      
      # Fit survival curve
      fit <- survfit(surv_obj ~ expression_group, data = current_data)
      
      # Log-rank test
      logrank_test <- survdiff(surv_obj ~ expression_group, data = current_data)
      p_value <- pchisq(logrank_test$chisq, df = 1, lower.tail = FALSE)
      
      # Sample sizes
      n_low <- sum(current_data$expression_group == "Low")
      n_high <- sum(current_data$expression_group == "High")
      
      # Create plot
      p_combo <- ggsurvplot(
        fit,
        data = current_data,
        pval = TRUE,
        pval.method = TRUE,
        conf.int = FALSE,
        risk.table = TRUE,
        risk.table.height = 0.3,
        ggtheme = theme_bw(),
        palette = c("#3B4992FF", "#EE0000FF"),
        legend.title = paste0(gene, " Expression"),
        legend.labs = c(paste0("Low (n=", n_low, ")"), paste0("High (n=", n_high, ")")),
        title = paste0(gene, " - ", threshold_name, " threshold (", round(cutoff_value, 3), ") - Combination Therapy Only"),
        xlab = "Time (months)",
        ylab = "Relapse-free survival probability"
      )
      cat("Combo only: p =", format.pval(p_value), "| Low: n =", n_low, "| High: n =", n_high, "\n")
      
      combined_plot <- arrange_ggsurvplots(list(p_all, p_chemo, p_combo), ncol = 3, nrow = 1)
      
      pdf(file.path(figureDir_temp,paste0(gene,"_km_with_cutoff_",threshold_name,".pdf")),width = 16,height = 5)
      print(combined_plot)
      dev.off()
  }
}

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE\n")
cat("Total genes analyzed:", length(gene_list), "\n")
cat(rep("=", 60), "\n")

# ================================================================================
# Hypothesis: High EC_EpCAM WITH COMPARISON SURVIVAL WITH COMBINATION THERAPY
# ================================================================================
malignant_patient_expression_chemo <- malignant_patient_expression[malignant_patient_expression$Treatment == "Chemo",]
malignant_patient_expression_combo <- malignant_patient_expression[malignant_patient_expression$Treatment == "Combo",]

# Calculate thresholds
gene<-"EpCAM"
thresholds_chemo <- find_optimal_cutoff(malignant_patient_expression_chemo, gene)
thresholds_combn <- find_optimal_cutoff(malignant_patient_expression_combo, gene)


  # Create groups
  current_data_chemo <- malignant_patient_expression_chemo %>%
    filter(!is.na(.data[[gene]])) %>%
    mutate(
      expression_group = ifelse(.data[[gene]] >= thresholds_chemo, "High", "Low")
    )
  
  current_data_combn <- malignant_patient_expression_combo %>%
    filter(!is.na(.data[[gene]])) %>%
    mutate(
      expression_group = ifelse(.data[[gene]] >= thresholds_combn, "High", "Low")
    )
  
  # Combine EC_EpCAM_high with combination therapy
  current_data <- rbind(current_data_chemo,current_data_combn)
  current_data <- current_data[current_data$expression_group == "High",]
  current_data$expression_group <- paste0(current_data$Treatment,"_",current_data$expression_group)
  
  # Create survival object
  current_data$RFS_time <- as.numeric(current_data$RFS_time)
  current_data$RFS_status <- as.numeric(current_data$RFS_status)
  
  surv_obj <- Surv(time = current_data$RFS_time, event = current_data$RFS_status)
  
  # Fit survival curve
  fit <- survfit(surv_obj ~ expression_group, data = current_data)
  
  # Log-rank test
  logrank_test <- survdiff(surv_obj ~ expression_group, data = current_data)
  p_value <- pchisq(logrank_test$chisq, df = 1, lower.tail = FALSE)
  
  # Sample sizes
  n_low <- sum(current_data$expression_group == "Combo_High")
  n_high <- sum(current_data$expression_group == "Chemo_High")
  
  # Create plot
  p_all <- ggsurvplot(
    fit,
    data = current_data,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.3,
    ggtheme = theme_bw(),
    palette = c("#8491b4ff", "#f39b7fff"),
    legend.title = paste0(gene, " Expression"),
    legend.labs = c(paste0("EC_EpCAM high with combination therapy (n=", n_low, ")"), paste0("EC_EpCAM high with chemotherapy (n=", n_high, ")")),
    title = paste0(gene, " - ", threshold_name, " threshold (", round(cutoff_value, 3), ")"),
    xlab = "Time (months)",
    ylab = "Relapse-free survival probability"
  )
  pdf(file.path(figureDir_temp,paste0(gene,"_km_comparsion_during_treatment_.pdf")),width = 8,height = 6)
  print(p_all)
  dev.off()
  
# ================================================================================
# Double marker: GLUT1 + FASN
# ================================================================================
  # Step 1: Combine the two datasets
  combined_data <- rbind(malignant_patient_expression_combo, malignant_patient_expression_chemo)
  
  # Step 1: Primary stratification by GLUT1 (chemo sufficiency indicator)
  glut1_cutoff <- median(combined_data$GLUT1)
  combined_data$GLUT1_group <- ifelse(combined_data$GLUT1 >= glut1_cutoff, "High_GLUT1", "Low_GLUT1")
  
  print(paste("GLUT1 median cutoff:", round(glut1_cutoff, 3)))
  print("Primary stratification by GLUT1:")
  print(table(combined_data$GLUT1_group, combined_data$Treatment))
  
  # Step 2: Secondary stratification by FASN (for high GLUT1 patients only)
  fasn_cutoff <- median(combined_data$FASN)
  combined_data$FASN_group <- ifelse(combined_data$FASN >= fasn_cutoff, "High_FASN", "Low_FASN")
  
  print(paste("FASN median cutoff:", round(fasn_cutoff, 3)))
  
  # Step 3: Create hierarchical treatment recommendation groups
  combined_data$treatment_recommendation <- case_when(
    combined_data$GLUT1_group == "Low_GLUT1" ~ "Chemotherapy_Sufficient",
    combined_data$GLUT1_group == "High_GLUT1" & combined_data$FASN_group == "High_FASN" ~ "Need_Target_Therapy",
    combined_data$GLUT1_group == "High_GLUT1" & combined_data$FASN_group == "Low_FASN" ~ "Alternative_Combination_Needed",
    TRUE ~ "Unclassified"
  )
  
  # Create actual treatment-recommendation combinations
  combined_data$treatment_recommendation_actual <- paste(combined_data$Treatment, combined_data$treatment_recommendation, sep = "_")
  
  # Display the hierarchical strategy results
  hierarchical_summary <- combined_data %>%
    group_by(treatment_recommendation, Treatment) %>%
    summarise(
      n_patients = n(),
      n_events = sum(RFS_status),
      event_rate = round(100 * n_events / n_patients, 1),
      median_rfs = round(median(RFS_time), 1),
      .groups = 'drop'
    )
  
  print("\n=== HIERARCHICAL STRATEGY RESULTS ===")
  print(hierarchical_summary)
  
  # Calculate treatment appropriateness
  treatment_appropriateness <- combined_data %>%
    mutate(
      appropriate_treatment = case_when(
        treatment_recommendation == "Chemotherapy_Sufficient" & Treatment == "Chemo" ~ "Appropriate",
        treatment_recommendation == "Need_Target_Therapy" & Treatment == "Combo" ~ "Appropriate",
        treatment_recommendation == "Alternative_Combination_Needed" ~ "Needs_Different_Approach",
        TRUE ~ "Inappropriate"
      )
    ) %>%
    group_by(appropriate_treatment) %>%
    summarise(
      n_patients = n(),
      n_events = sum(RFS_status),
      event_rate = round(100 * n_events / n_patients, 1),
      median_rfs = round(median(RFS_time), 1),
      .groups = 'drop'
    )
  
  print("\n=== TREATMENT APPROPRIATENESS ANALYSIS ===")
  print(treatment_appropriateness)
  
  # Detailed analysis by recommendation group
  print("\n=== DETAILED ANALYSIS BY RECOMMENDATION GROUP ===")
  
  # Group 1: Chemotherapy Sufficient (Low GLUT1)
  chemo_sufficient <- combined_data[combined_data$treatment_recommendation == "Chemotherapy_Sufficient", ]
  print(paste("Group 1 - Chemotherapy Sufficient (Low GLUT1):", nrow(chemo_sufficient), "patients"))
  
  if (nrow(chemo_sufficient) > 5) {
    chemo_sufficient_summary <- chemo_sufficient %>%
      group_by(Treatment) %>%
      summarise(
        n = n(),
        events = sum(RFS_status),
        event_rate = round(100 * events / n, 1),
        median_rfs = round(median(RFS_time), 1),
        .groups = 'drop'
      )
    print("Treatment outcomes in Chemotherapy Sufficient group:")
    print(chemo_sufficient_summary)
  }
  
  # Group 2: Need Target Therapy (High GLUT1 + High FASN)
  need_target <- combined_data[combined_data$treatment_recommendation == "Need_Target_Therapy", ]
  print(paste("Group 2 - Need Target Therapy (High GLUT1 + High FASN):", nrow(need_target), "patients"))
  
  if (nrow(need_target) > 5) {
    need_target_summary <- need_target %>%
      group_by(Treatment) %>%
      summarise(
        n = n(),
        events = sum(RFS_status),
        event_rate = round(100 * events / n, 1),
        median_rfs = round(median(RFS_time), 1),
        .groups = 'drop'
      )
    print("Treatment outcomes in Need Target Therapy group:")
    print(need_target_summary)
  }
  
  # Group 3: Alternative Combination Needed (High GLUT1 + Low FASN)
  alt_combination <- combined_data[combined_data$treatment_recommendation == "Alternative_Combination_Needed", ]
  print(paste("Group 3 - Alternative Combination Needed (High GLUT1 + Low FASN):", nrow(alt_combination), "patients"))
  
  if (nrow(alt_combination) > 5) {
    alt_combination_summary <- alt_combination %>%
      group_by(Treatment) %>%
      summarise(
        n = n(),
        events = sum(RFS_status),
        event_rate = round(100 * events / n, 1),
        median_rfs = round(median(RFS_time), 1),
        .groups = 'drop'
      )
    print("Treatment outcomes in Alternative Combination Needed group:")
    print(alt_combination_summary)
  }
  
  # Statistical testing of the hierarchical strategy
  print("\n=== STATISTICAL VALIDATION ===")
  
  surv_obj <- Surv(time = combined_data$RFS_time, event = combined_data$RFS_status)
  
  # Test 1: Overall difference between recommendation groups
  overall_test <- survdiff(surv_obj ~ treatment_recommendation, data = combined_data)
  overall_pval <- 1 - pchisq(overall_test$chisq, length(overall_test$n) - 1)
  print(paste("Overall test across recommendation groups: p =", round(overall_pval, 4)))
  
  # Test 2: Treatment benefit in each recommendation group
  test_results <- data.frame()
  
  for (group in unique(combined_data$treatment_recommendation)) {
    group_data <- combined_data[combined_data$treatment_recommendation == group, ]
    
    if (nrow(group_data) > 10 && length(unique(group_data$Treatment)) > 1) {
      surv_group <- Surv(time = group_data$RFS_time, event = group_data$RFS_status)
      group_test <- survdiff(surv_group ~ Treatment, data = group_data)
      group_pval <- 1 - pchisq(group_test$chisq, length(group_test$n) - 1)
      
      test_results <- rbind(test_results, data.frame(
        Group = group,
        P_value = round(group_pval, 4),
        N_patients = nrow(group_data),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  print("Treatment effect within each recommendation group:")
  print(test_results)
  
  # Create Kaplan-Meier plots for the hierarchical strategy
  print("\n=== KAPLAN-MEIER VISUALIZATION ===")
  
  # Plot 1: Overall survival by treatment recommendation
  surv_fit_recommendation <- survfit(surv_obj ~ treatment_recommendation, data = combined_data)
  
  recommendation_colors <- c(
    "Chemotherapy_Sufficient" = "#7ca6dc",           # Blue - chemo OK
    "Need_Target_Therapy" = "#cd534c",               # Red - needs target therapy
    "Alternative_Combination_Needed" = "#ff9900"     # Orange - needs alternative
  )
  
  km_recommendation <- ggsurvplot(
    surv_fit_recommendation,
    data = combined_data,
    pval = TRUE,
    conf.int = TRUE,
    palette = unname(recommendation_colors),
    ggtheme = theme_bw(),
    title = "Hierarchical Biomarker Strategy",
    subtitle = "GLUT1 (Tier 1) → FASN (Tier 2) Decision Tree",
    xlab = "Time to Recurrence (months)",
    ylab = "Recurrence-Free Survival",
    legend.title = "Treatment Recommendation",
    risk.table = TRUE,
    risk.table.height = 0.3
  )
  
  pdf(file.path(figureDir,"Overall survival of treatment strategy.pdf"), width = 8, height = 6)
  print(km_recommendation)
  dev.off()

  # Plot 2: Treatment-stratified analysis for key groups
  key_groups <- c("Chemotherapy_Sufficient", "Need_Target_Therapy")
  key_group_data <- combined_data[combined_data$treatment_recommendation %in% key_groups, ]
  key_group_data$treatment_recommendation_actual <- paste(key_group_data$Treatment, key_group_data$treatment_recommendation, sep = "_")
  
  surv_obj_key <- Surv(time = key_group_data$RFS_time, event = key_group_data$RFS_status)
  surv_fit_key <- survfit(surv_obj_key ~ treatment_recommendation_actual, data = key_group_data)
  
  key_colors <- c(
    "Chemo_Chemotherapy_Sufficient" = "#7ca6dc",     # Blue - appropriate match
    "Combo_Chemotherapy_Sufficient" = "#b3d1ff",     # Light blue - over-treatment
    "Chemo_Need_Target_Therapy" = "#ffcccc",         # Light red - under-treatment  
    "Combo_Need_Target_Therapy" = "#cd534c"          # Red - appropriate match
  )
  
  km_key_groups <- ggsurvplot(
    surv_fit_key,
    data = key_group_data,
    pval = TRUE,
    conf.int = FALSE,
    palette = unname(key_colors),
    ggtheme = theme_bw(),
    title = "Treatment Appropriateness Analysis",
    subtitle = "Matching biomarker recommendation to actual treatment",
    xlab = "Time to Recurrence (months)",
    ylab = "Recurrence-Free Survival",
    legend.title = "Treatment + Recommendation",
    risk.table = TRUE,
    risk.table.height = 0.3
  )
  
  pdf(file.path(figureDir,"Overall survival of Treatment-stratified analysis for key groupsy.pdf"), width = 8, height = 6)
  print(km_key_groups)
  dev.off()
  
  # Create decision tree visualization
  print("\n=== CLINICAL DECISION TREE ===")
  
  decision_tree_summary <- combined_data %>%
    group_by(GLUT1_group, FASN_group, treatment_recommendation) %>%
    summarise(
      n_patients = n(),
      event_rate_chemo = round(100 * sum(RFS_status[Treatment == "Chemo"]) / sum(Treatment == "Chemo"), 1),
      event_rate_combo = round(100 * sum(RFS_status[Treatment == "Combo"]) / sum(Treatment == "Combo"), 1),
      .groups = 'drop'
    ) %>%
    mutate(
      combo_benefit = event_rate_chemo - event_rate_combo,  # Positive = combo better
      recommendation_logic = case_when(
        treatment_recommendation == "Chemotherapy_Sufficient" ~ "Low GLUT1 → Chemo OK",
        treatment_recommendation == "Need_Target_Therapy" ~ "High GLUT1 + High FASN → Add Target",
        treatment_recommendation == "Alternative_Combination_Needed" ~ "High GLUT1 + Low FASN → Try Alternative",
        TRUE ~ "Other"
      )
    )
  
  print("Decision tree with clinical outcomes:")
  print(decision_tree_summary[, c("GLUT1_group", "FASN_group", "recommendation_logic", "n_patients", "combo_benefit")])
  
  # Recreate the key summary data for visualization
  hierarchical_viz_data <- hierarchical_summary
  
  # Prepare data for enhanced visualization
  outcome_comparison_data <- hierarchical_viz_data %>%
    mutate(
      treatment_recommendation = factor(treatment_recommendation,
                                        levels = c("Chemotherapy_Sufficient", "Need_Target_Therapy", "Alternative_Combination_Needed"),
                                        labels = c("Chemotherapy Sufficient\n(Low GLUT1)", 
                                                   "Need Target Therapy\n(High GLUT1 + High FASN)", 
                                                   "Alternative Needed\n(High GLUT1 + Low FASN)")),
      Treatment = factor(Treatment, levels = c("Chemo", "Combo")),
      # Create enhanced label with both event rate and median survival
      enhanced_label = paste0(event_rate, "% events\n", 
                              median_rfs, " mo\n", 
                              "(n=", n_patients, ")")
    )
  
  # Create enhanced side-by-side bar plot
  enhanced_outcome_plot <- ggplot(outcome_comparison_data, aes(x = treatment_recommendation, y = event_rate, fill = Treatment)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = enhanced_label), 
              position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c("Chemo" = "#7ca6dc", "Combo" = "#cd534c"),
                      labels = c("Chemotherapy Only", "Combination Therapy")) +
    scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, 20)) +
    labs(
      title = "Hierarchical Biomarker Strategy: Treatment Outcomes",
      subtitle = "Event rates and median survival by biomarker recommendation and treatment",
      x = "Biomarker-Based Treatment Recommendation",
      y = "Early Recurrence Rate (%)",
      fill = "Treatment Received",
      caption = "Lower event rates and longer survival times indicate better outcomes"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0.5),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  pdf(file.path(figureDir,"Treatment-stratified comprasion analysis.pdf"), width = 8, height = 6)
  print(enhanced_outcome_plot)
  dev.off()
  
  # Save all results
  write.csv(combined_data, "Hierarchical_Biomarker_Data.csv", row.names = FALSE)
  write.csv(hierarchical_summary, "Hierarchical_Strategy_Summary.csv", row.names = FALSE)
  write.csv(treatment_appropriateness, "Treatment_Appropriateness_Analysis.csv", row.names = FALSE)
  write.csv(decision_tree_summary, "Clinical_Decision_Tree_Summary.csv", row.names = FALSE)

# ================================================================================
# Multivariable Cox regression of malignant gene expression
# ================================================================================
malignant_patient_expression <- malignant_roi_expression[,-match(c("sample_id","Tissue"),colnames(malignant_roi_expression))]
malignant_patient_expression <- malignant_patient_expression %>% 
  group_by(patient_id, Treatment, RFS_status, RFS_group) %>%
  summarise(across(c(1:(ncol(malignant_patient_expression) - 4)), mean, na.rm = TRUE))
malignant_patient_expression <- as.data.frame(malignant_patient_expression)

clinical_predictors <- c(
  "RFS_time","Age", "Gender", "fong_score", "TBS", 
  "CEA", "CA199", "CRLM_number", "CRLM_size", 
  "T_stage", "Differential_grade", "Lymph_positive", 
  "KRAS_mutation"
)

for(add_feature in clinical_predictors){
  malignant_patient_expression[[add_feature]] <- colData(spe)[match(malignant_patient_expression$patient_id,spe$patient_id),add_feature]
}

head(malignant_patient_expression)

malignant_patient_expression_scaled <- malignant_patient_expression %>% ## log for CEA and CA19-9
  mutate(
    CEA_log = log(CEA + 1),  # +1 to handle any zeros
    CA199_log = log(CA199 + 1)
  )

# Identify continuous variables to standardize
continuous_vars <- c("EpCAM", "Vimentin", "Ki67", "VEGF", "CA_IX", 
                     "GLUT1", "HK2", "FASN", "Age", "TBS", "CEA", "CA199",
                     "CEA_log", "CA199_log", "CRLM_size", "CRLM_number")

# Standardize continuous variables (mean=0, sd=1)
malignant_patient_expression_scaled <- malignant_patient_expression_scaled %>%
  mutate(across(all_of(continuous_vars), ~ scale(.)[,1]))

# Keep categorical variables as-is
categorical_vars <- c("Treatment", "Gender", "fong_score", 
                      "T_stage", "Differential_grade", "Lymph_positive", 
                      "KRAS_mutation")

clinical_cox_model <- coxph(Surv(RFS_time, RFS_status) ~ 
                              Treatment +
                              Age + fong_score + TBS + CEA + CA199 + CRLM_size + Lymph_positive +
                              GLUT1 + FASN,
                            data = malignant_patient_expression_scaled)

custom_names <- c("Treatment (Combo vs Chemo)", "Age", "Fong Score", "TBS",
                  "CEA", "CA199", "CRLM Size","Lymph Node Positive",
                  "GLUT1","FASN")

var_groups <- list(
  "Treatment" = 1,
  "Patient Factors" = 2:8,  # Age, fong_score, TBS, CEA_log
  "Tumor Characteristics" = 9:10
)
result <- createCoxForestPlot(clinical_cox_model,
                              var_display_names = custom_names,
                              var_groups = var_groups,
                              savePath = figureDir,
                              filename = "Multi-cox of clinical variables.pdf",
                              plot_title = "Clinical Variables Cox Regression")

multicox_with_interaction <- coxph(Surv(RFS_time, RFS_status) ~ 
                                     Age + fong_score + TBS + CEA + CA199 + 
                                     Treatment*FASN + Treatment*GLUT1,
                   data = malignant_patient_expression_scaled)

custom_names <- c("Age", "Fong Score", "TBS",
                  "CEA", "CA199", "Treatment (Combo vs Chemo)",
                  "FASN","GLUT1",
                  "Combo:FASN","Combo:GLUT1")

var_groups <- list(
  "Treatment" = 6,
  "Patient Factors" = 1:5,  # Age, fong_score, TBS, CEA_log
  "Tumor Characteristics" = 7:10
)
result <- createCoxForestPlot(multicox_with_interaction,
                              var_display_names = custom_names,
                              var_groups = var_groups,
                              savePath = figureDir,
                              filename = "Multi-cox of variables with interaction factor.pdf",
                              plot_title = "Clinical Variables Cox Regression")

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
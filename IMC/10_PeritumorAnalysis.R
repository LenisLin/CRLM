# For Peritumor cell analysis

# =============================================================================
# LOAD PACKAGES
# =============================================================================
library(SpatialExperiment)
library(imcRtools)
library(igraph)

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
library(ConsensusClusterPlus)
library(survival)
library(survminer)

library(dplyr)
library(tidyr)

# =============================================================================
# LOAD AND VALIDATE DATA
# =============================================================================
codeSpace <- "/home/lenislin/Experiment/projects/CRLM_2025/IMC"

source(file.path(codeSpace,"Peritumor_analysis_functions.R"))

# Set Working Directory and Source Utilities
workDir <- "/mnt/public/lyx/IMC_HE_Merge/CRLM"
setwd(workDir)

# Define directories for saving results and Steinbock outputs
saveDir <- file.path(workDir, "results")
steinResultDir <- file.path(workDir, "steinbock")
figureDir <- file.path(workDir, "figures","5_PTAnalysis")

if(!dir.exists(figureDir)){
  dir.create(figureDir,recursive = T)
}

# Load your SpatialExperiment object
date_time <- "0730"
spe <- readRDS(file.path(saveDir,paste0("spatial_spe_",date_time,".rds")))
img_id_ <- "sample_id"

# Defined color for CN
cn_column <- "CN_knn_20_cluster_10"  # CHANGE THIS - choose your CN parameter
CN_color <- setNames(metadata(spe)$color_vectors$color_10,unique(colData(spe)[[cn_column]]))
RFS_color <- setNames(c("#EFC000FF", "#0073C2FF"),unique(spe$RFS_status))
Tissue_color <- setNames(metadata(spe)$color_vectors$color_20[1:3],unique(spe$Tissue)[1:3])

# =============================================================================
# Step 0: Extract cells from PT tissues
# =============================================================================
spe_pt <- spe[,spe$Tissue %in% c("PT")]

# Extract metadata and expression data for easier handling
pt_meta <- as.data.frame(colData(spe_pt))
pt_expr <- as.data.frame(assay(spe_pt))

cat("PT Analysis: ", ncol(spe_pt), "cells from", 
    length(unique(pt_meta$patient_id)), "patients\n")

# =============================================================================
# Step 1: EXPRESSION DYNAMICS: Compare markers between clinical groups
# =============================================================================
cat("\n=== 1. Expression Dynamics Analysis ===\n")

# Calculate patient-level marker expression by cell type
expr_by_patient <- pt_meta %>%
  bind_cols(as.data.frame(t(pt_expr))) %>%
  group_by(patient_id, major_celltype, Treatment, RFS_status) %>%
  summarise(across(all_of(rownames(pt_expr)), mean, na.rm = TRUE), .groups = "drop")

# Test each marker in each cell type for R vs NR differences
expr_results <- data.frame()

for (celltype in unique(pt_meta$major_celltype)) {
  for (treatment in c("Chemo", "Combo")) {
    
    subset_data <- expr_by_patient %>%
      filter(major_celltype == celltype, Treatment == treatment)
    
    if (nrow(subset_data) < 4) next  # Skip if too few samples
    
    for (marker in rownames(pt_expr)) {
      
      r_values <- subset_data[[marker]][subset_data$RFS_status == 1]
      nr_values <- subset_data[[marker]][subset_data$RFS_status == 0]
      
      if (length(r_values) > 0 & length(nr_values) > 0) {
        test_result <- t.test(r_values, nr_values)
        
        expr_results <- rbind(expr_results, data.frame(
          celltype = celltype,
          treatment = treatment,
          marker = marker,
          R_mean = mean(r_values, na.rm = TRUE),
          NR_mean = mean(nr_values, na.rm = TRUE),
          fold_change = mean(r_values, na.rm = TRUE) / mean(nr_values, na.rm = TRUE),
          p_value = test_result$p.value
        ))
      }
    }
  }
}

expr_results$log2_fc <- log2(expr_results$fold_change)
expr_results$p_adjust <- p.adjust(expr_results$p_value,method = "BH")
expr_results$neg_log10_p <- -log10(expr_results$p_value)
expr_results$neg_log10_p_adj <- -log10(expr_results$p_adjust)

# =============================================================================
# Step 2: CELL ABUNDANCE: R vs NR differences
# =============================================================================

cat("\n=== 2. Cell Abundance Analysis ===\n")

# Calculate cell density per patient (Subtype)
cell_density <- pt_meta %>%
  group_by(patient_id, sub_celltype, Treatment, RFS_status) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(total_cells = sum(cell_count),
         density = cell_count / total_cells * 100) %>%
  ungroup()

# Test abundance differences
abundance_results <- cell_density %>%
  group_by(sub_celltype, Treatment) %>%
  summarise(
    R_mean = mean(density[RFS_status == 1], na.rm = TRUE),
    NR_mean = mean(density[RFS_status == 0], na.rm = TRUE),
    fold_change = R_mean / NR_mean,
    p_value = tryCatch({
      wilcox.test(density[RFS_status == 1], density[RFS_status == 0])$p.value
    }, error = function(e) NA),
    .groups = "drop"
  ) %>%
  filter(!is.na(p_value)) %>%
  mutate(
    log2_fc = log2(fold_change),
    neg_log10_p = -log10(p_value),
    significant = p_value < 0.05
  )

abundance_results$p_adjust <- p.adjust(abundance_results$p_value,method = "BH")
abundance_results$neg_log10_p_adj <- -log10(abundance_results$p_adjust)

# =============================================================================
# Step 3: CN ABUNDANCE: R vs NR differences
# =============================================================================

cat("\n=== 3. Cell Neighbors Analysis ===\n")

# Calculate cell density per patient
cn_density <- pt_meta %>%
  group_by(patient_id, CN_knn_20_cluster_10, Treatment, RFS_status) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(total_cells = sum(cell_count),
         density = cell_count / total_cells * 100) %>%
  ungroup()

# Test abundance differences
cn_abundance_results <- cn_density %>%
  group_by(CN_knn_20_cluster_10, Treatment) %>%
  summarise(
    R_mean = mean(density[RFS_status == 1], na.rm = TRUE),
    NR_mean = mean(density[RFS_status == 0], na.rm = TRUE),
    fold_change = R_mean / NR_mean,
    p_value = tryCatch({
      wilcox.test(density[RFS_status == 1], density[RFS_status == 0])$p.value
    }, error = function(e) NA),
    .groups = "drop"
  ) %>%
  filter(!is.na(p_value)) %>%
  mutate(
    log2_fc = log2(fold_change),
    neg_log10_p = -log10(p_value),
    significant = p_value < 0.05
  )

cn_abundance_results$p_adjust <- p.adjust(cn_abundance_results$p_value,method = "BH")
cn_abundance_results$neg_log10_p_adj <- -log10(cn_abundance_results$p_adjust)

# ============================================================================
# 4. BUBBLE PLOTS: Visualize results from 1-3
# ============================================================================
library(ggrepel)

cat("\n=== 4. Creating Bubble Plots ===\n")

# Plot 1: Expression dynamics (multi-panel by cell type)
expr_plot_data <- expr_results %>%
  filter(!is.na(p_value)) %>%
  prep_plot_data("marker")

p1 <- expr_plot_data %>%
  ggplot(aes(x = celltype, y = log2_fc)) +
  geom_jitter(aes(alpha = alpha_val, size = size_val, color = direction),width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf) +
  facet_wrap(~treatment) +
  scale_color_manual(values = c("R_higher" = "#EFC000FF", "NR_higher" = "#0073C2FF", "n.s." = "gray70")) +
  scale_alpha_identity() +
  scale_size_identity() +
  labs(title = "Expression Dynamics: R vs NR by Treatment",
       x = "Treatment", y = "Log2 Fold Change (R/NR)",
       color = "Higher in") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: Cell abundance 
abundance_plot_data <- abundance_results %>%
  prep_plot_data("sub_celltype")

p2 <- abundance_plot_data %>%
  ggplot(aes(x = Treatment, y = log2_fc)) +
  geom_jitter(aes(alpha = alpha_val, size = size_val, color = direction),width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf) +
  scale_color_manual(values = c("R_higher" = "#EFC000FF", "NR_higher" = "#0073C2FF", "n.s." = "gray70")) +
  scale_alpha_identity() +
  scale_size_identity() +
  labs(title = "Cell Abundance: R vs NR by Treatment",
       x = "Treatment", y = "Log2 Fold Change (R/NR)",
       color = "Higher in") +
  theme_minimal()

# Plot 3: CN abundance
cn_plot_data <- cn_abundance_results %>%
  prep_plot_data("CN_knn_20_cluster_10")
cn_plot_data$label <- ifelse(cn_plot_data$label=="","",paste0("CN",cn_plot_data$label)) 

p3 <- cn_plot_data %>%
  ggplot(aes(x = Treatment, y = log2_fc)) +
  geom_jitter(aes(alpha = alpha_val, size = size_val, color = direction),width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf) +
  scale_color_manual(values = c("R_higher" = "#EFC000FF", "NR_higher" = "#0073C2FF", "n.s." = "gray70")) +
  scale_alpha_identity() +
  scale_size_identity() +
  labs(title = "CN Cluster Abundance: R vs NR by Treatment", 
       x = "Treatment", y = "Log2 Fold Change (R/NR)",
       color = "Higher in") +
  theme_minimal()

ggsave(filename = file.path(figureDir,"Marker expression dynamic of major types.pdf"),p1,width = 10,height = 5)
ggsave(filename = file.path(figureDir,"Cell subpopulation abudance dynamic of major types.pdf"),p2,width = 6,height = 5)
ggsave(filename = file.path(figureDir,"CN types abudance dynamic of major types.pdf"),p3,width = 6,height = 5)

# =============================================================================
# Step 5: EC CO-LOCALIZATION ANALYSIS
# =============================================================================

cat("\n=== 5. EC Co-localization Analysis ===\n")

# Identify EC cells (EpCAM+ epithelial cells in PT)
ec_cells <- pt_meta %>%
  filter(grepl("^EC_", sub_celltype)) %>%
  pull(sub_celltype) %>%
  unique()

cat("Identified EC subtypes:", paste(ec_cells, collapse = ", "), "\n")

cat("\n--- Ripley's K Function Analysis ---\n")
# Analyze EC co-localization with major cell types
major_types <- unique(pt_meta$major_celltype)
major_types <- major_types[!is.na(major_types)]

ripleys_results <- list()

meta_all <- colData(spe_pt)
coor_all <- spatialCoords(spe_pt)

non_ec_subtypes <- unique(meta_all$sub_celltype)
non_ec_subtypes <- non_ec_subtypes[!startsWith(non_ec_subtypes,"EC")]

for (sample in unique(spe_pt$sample_id)) {
  index_ <- meta_all$sample_id == sample
  sample_meta <- meta_all[index_,]
  sample_coor <- coor_all[index_,]
  
  # Get patient info
  patient_id <- unique(sample_meta$patient_id)
  treatment <- unique(sample_meta$Treatment)
  rfs_status <- unique(sample_meta$RFS_status)
  
  cat("Analyzing sample:", sample, "Patient:", patient_id, "\n")
  
  # Test EC against each major cell type
  
  for (target_type in non_ec_subtypes) {
    # target_cells <- sample_meta$sub_celltype[sample_meta$major_celltype == target_type]
    # if (length(target_cells) < 5) next
    # # Use the most abundant subtype as representative
    # target_subtype <- names(sort(table(target_cells), decreasing = TRUE))[1]
    
    k_result <- perform_ripleys_k(coords = sample_coor, sample_meta = sample_meta, focal_celltype = "EC", target_type, max_dist = 150)
    
    if (!is.null(k_result)) {
      k_result$sample_id <- sample
      k_result$patient_id <- patient_id
      k_result$treatment <- treatment
      k_result$rfs_status <- rfs_status
      k_result$focal_type <- "EC"
      k_result$target_type <- target_type

      ripleys_results[[paste(sample, "EC", target_type, sep = "_")]] <- k_result
    }
  }
}

# Combine all results
ripleys_combined <- do.call(rbind, ripleys_results)
rownames(ripleys_combined) <- NULL

# =============================================================================
# Step 6: L-FUNCTION VISUALIZATION ANALYSIS
# =============================================================================

cat("\n=== 6. L-function Visualization Analysis ===\n")

# Step 6.1: Data Preparation
# Add some summary metrics for better analysis
ripleys_summary <- ripleys_combined %>%
  mutate(
    RFS_group = factor(rfs_status, levels = c(0, 1), labels = c("No Recurrence", "Recurrence")),
    Treatment_RFS = paste(treatment, RFS_group, sep = "_")
  )
ripleys_summary$log2_L <- sign(ripleys_summary$L) * log2(abs(ripleys_summary$L))

# Step 6.2: Plot 1 - Heatmaps for specific distances (20, 50, 100μm)
cat("\n--- 6.2: Creating Distance-specific Heatmaps ---\n")

heatmap_distances <- c(20, 50, 100)  # Specific distances for heatmaps
# Create heatmaps for each specific distance
for(dist in heatmap_distances) {
  cat("Creating heatmap for distance:", dist, "μm\n")
  
  # Prepare heatmap data for this distance
  heatmap_matrix <- ripleys_summary %>%
    filter(distance == dist) %>%
    select(sample_id, target_type, log2_L) %>%
    group_by(sample_id, target_type) %>%
    summarise(log2_L = mean(log2_L), .groups = "drop") %>%
    pivot_wider(names_from = target_type, values_from = log2_L) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()
  
  # Create the heatmap
  p_heatmap_dist <- pheatmap(
    heatmap_matrix,
    main = paste0("EC Co-localization Heatmap at ", dist, "μm"),
    color = colorRampPalette(c("blue", "white", "red"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    angle_col = 45,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 10,
    border_color = "white"
  )
  
  # Save immediately
  ggsave(file.path(figureDir, paste0("EC_colocalization_heatmap_", dist, "um.pdf")), 
         p_heatmap_dist, width = 12, height = 8)
  cat("Saved heatmap for", dist, "μm\n")
}

# Step 6.3: Plot 2 - Split Violin Plots with Statistical Testing
cat("\n--- 6.3: Creating Split Violin Plots with Statistics ---\n")

# Calculate key distance summaries for violin plots
key_distances <- c(20, 30, 50, 75, 100, 150)  # Focus on biologically relevant distances

ripleys_key_distances <- ripleys_summary %>%
  filter(distance %in% key_distances) %>%
  mutate(distance_label = paste0(distance, "μm"))

# Get key target types
key_target_types <- c("CD8T","B","Macro_CD163","SC_Collagen","SC_Collagen_Vimentin")

# Ensure proper factor levels for RFS
ripleys_key_distances$rfs_status <- factor(ripleys_key_distances$rfs_status)

# Plot 2A: Split violin plot for RFS comparison
ripleys_for_violin <- ripleys_key_distances %>%
  filter(target_type %in% key_target_types)

# Calculate statistical annotations
violin_stats <- ripleys_for_violin %>%
  group_by(target_type, distance_label) %>%
  group_modify(~ add_significance(.x, max(.x$log2_L, na.rm = TRUE) * 1.1)) %>%
  ungroup()

p_violin_rfs <- ripleys_for_violin %>%
  ggplot(aes(x = rfs_status, y = log2_L)) +
  geom_violin(aes(fill = rfs_status), alpha = 0.7, scale = "width", 
              position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(fill = rfs_status), width = 0.2, alpha = 0.8,
               position = position_dodge(width = 0.8)) +
  geom_text(data = violin_stats, aes(x = x, y = y, label = label), 
            size = 3, hjust = 0.5) +
  facet_grid(distance_label ~ target_type, scales = "free_y") +
  scale_fill_manual(values = RFS_color, name = "RFS Status") +
  labs(
    title = "Split Violin Plot: L-values by Recurrence Status",
    subtitle = "Distribution of clustering patterns with Wilcoxon test results",
    x = "RFS Status",
    y = "log2(L-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save immediately
ggsave(file.path(figureDir, "EC_colocalization_split_violin_rfs.pdf"), 
       p_violin_rfs, width = 14, height = 12)
cat("Saved split violin plot for RFS comparison\n")

# Plot 2B: Treatment-stratified split violin plot
ripleys_treatment_violin <- ripleys_key_distances %>%
  filter(target_type %in% key_target_types)

# Calculate treatment-specific statistics
treatment_violin_stats <- ripleys_treatment_violin %>%
  group_by(target_type, distance_label, treatment) %>%
  group_modify(~ add_significance(.x, max(.x$log2_L, na.rm = TRUE) * 1.1)) %>%
  ungroup()

p_violin_treatment <- ripleys_treatment_violin %>%
  ggplot(aes(x = rfs_status, y = log2_L)) +
  geom_violin(aes(fill = rfs_status), alpha = 0.7, scale = "width") +
  geom_boxplot(aes(fill = rfs_status), width = 0.2, alpha = 0.8) +
  geom_text(data = treatment_violin_stats, aes(x = x, y = y, label = label), 
            size = 2.5, hjust = 0.5) +
  facet_grid(distance_label ~ target_type + treatment, scales = "free_y") +
  scale_fill_manual(values = RFS_color, name = "RFS Status") +
  labs(
    title = "Treatment-Stratified Split Violin Plot",
    subtitle = "L-values by RFS Status within each treatment group",
    x = "RFS Status",
    y = "L-value"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  )

# Save immediately
ggsave(file.path(figureDir, "EC_colocalization_treatment_violin.pdf"), 
       p_violin_treatment, width = 16, height = 12)
cat("Saved treatment-stratified violin plot\n")

# Step 6.4: Plot 3 - Combined Line Plots (3 rows × n columns)
cat("\n--- 6.4: Creating Combined Line Plots (3×n layout) ---\n")

# Prepare data for line plots with more points
line_data_detailed <- ripleys_summary %>%
  filter(target_type %in% key_target_types) %>%  # Focus on top 4 target types
  mutate(
    treatment_category = case_when(
      TRUE ~ "All Patients",
      treatment == "Chemo" ~ "Chemotherapy Only", 
      treatment == "Combo" ~ "Combination Therapy"
    )
  )

# Create separate datasets for each row
line_data_all <- line_data_detailed %>%
  mutate(treatment_category = "All Patients") %>%
  group_by(distance, target_type, rfs_status, RFS_group, treatment_category) %>%
  summarise(
    mean_L = mean(log2_L, na.rm = TRUE),
    se_L = sd(log2_L, na.rm = TRUE) / sqrt(n()),
    n_samples = n(),
    .groups = "drop"
  )

line_data_chemo <- line_data_detailed %>%
  filter(treatment == "Chemo") %>%
  mutate(treatment_category = "Chemotherapy Only") %>%
  group_by(distance, target_type, rfs_status, RFS_group, treatment_category) %>%
  summarise(
    mean_L = mean(log2_L, na.rm = TRUE),
    se_L = sd(log2_L, na.rm = TRUE) / sqrt(n()),
    n_samples = n(),
    .groups = "drop"
  )

line_data_combo <- line_data_detailed %>%
  filter(treatment == "Combo") %>%
  mutate(treatment_category = "Combination Therapy") %>%
  group_by(distance, target_type, rfs_status, RFS_group, treatment_category) %>%
  summarise(
    mean_L = mean(log2_L, na.rm = TRUE),
    se_L = sd(log2_L, na.rm = TRUE) / sqrt(n()),
    n_samples = n(),
    .groups = "drop"
  )

# Combine all data
line_data_combined <- bind_rows(line_data_all, line_data_chemo, line_data_combo) %>%
  mutate(treatment_category = factor(treatment_category, 
                                     levels = c("All Patients", "Chemotherapy Only", "Combination Therapy")))

# Create the combined plot
line_data_combined$rfs_status <- as.factor(line_data_combined$rfs_status)
p_line_combined <- line_data_combined %>%
  ggplot(aes(x = distance, y = mean_L, color = rfs_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_point(size = 1, alpha = 0.8) +  # More points with smaller size
  geom_smooth(aes(group = rfs_status), method = "loess", se = FALSE, size = 1) +
  geom_ribbon(aes(ymin = mean_L - se_L, ymax = mean_L + se_L, fill = rfs_status), 
              alpha = 0.15, color = NA) +
  facet_grid(treatment_category ~ target_type, scales = "free_y") +
  scale_color_manual(values = RFS_color, name = "RFS Status") +
  scale_fill_manual(values = RFS_color, name = "RFS Status") +
  scale_x_continuous(breaks = seq(20, 200, 40)) +
  labs(
    title = "L-function Curves: Combined Analysis",
    subtitle = "EC co-localization patterns across treatment contexts",
    x = "Distance (μm)",
    y = "Mean L-value (±SE)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save immediately
ggsave(file.path(figureDir, "EC_colocalization_combined_lines.pdf"), 
       p_line_combined, width = 16, height = 12)
cat("Saved combined line plot (3×n layout)\n")

# =============================================================================
# Step 6.7: Save All Results to CSV
# =============================================================================

cat("\n--- 6.7: Saving Results to CSV ---\n")

# Save all analysis results
write.csv(ripleys_summary, file.path(saveDir, "EC_ripleys_complete_analysis.csv"), row.names = FALSE)
write.csv(violin_stats, file.path(saveDir, "EC_violin_plot_statistics.csv"), row.names = FALSE)
write.csv(treatment_violin_stats, file.path(saveDir, "EC_treatment_violin_statistics.csv"), row.names = FALSE)

# =============================================================================
# Interaction based analysis
# =============================================================================
library(scales)
library(BiocParallel)
out <- testInteractions(spe_pt, 
                        group_by = "sample_id",
                        label = "major_celltype", 
                        colPairName = "knn_20",
                        BPPARAM = SerialParam(RNGseed = 619))



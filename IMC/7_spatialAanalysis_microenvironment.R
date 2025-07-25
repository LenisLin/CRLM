# For spatial analysis

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
figureDir <- file.path(workDir, "figures","3_SpatialAnalysis","CN_Analysis")

if(!dir.exists(figureDir)){
  dir.create(figureDir,recursive = T)
}

# Load your SpatialExperiment object
date_time <- "0722"
spe <- readRDS(file.path(saveDir,paste0("spatial_spe_",date_time,".rds")))
img_id_ <- "sample_id"

# ============================================================================
# Task 1: CN-Celltype Composition Heatmap with Significance Testing
# ============================================================================

# Choose CN parameter to analyze
cn_column <- "CN_knn_20_cluster_10"  # CHANGE THIS - choose your CN parameter

# Focus tissues
tissue_focus <- c("TC", "IM", "PT")  # Tumor Core and Invasive Margin

print("Starting Task 1: CN-Celltype Composition Analysis")
print(paste("CN parameter:", cn_column))
print(paste("Focus tissues:", paste(tissue_focus, collapse = ", ")))

# Filter data for analysis
spe_filtered <- spe[, spe$Tissue %in% tissue_focus]
print(paste("Cells after tissue filtering:", ncol(spe_filtered)))

# Extract metadata
meta_df <- as.data.frame(colData(spe_filtered))
print(paste("Total cells in analysis:", nrow(meta_df)))

# Check if CN column exists
if (!cn_column %in% colnames(meta_df)) {
  print("Available CN columns:")
  cn_cols <- colnames(meta_df)[grepl("^CN_", colnames(meta_df))]
  print(cn_cols)
  stop(paste("CN column not found:", cn_column))
}

# Remove rows with missing CN or celltype data
meta_df <- meta_df[!is.na(meta_df[[cn_column]]) & !is.na(meta_df$sub_celltype), ]
print(paste("Cells after removing missing data:", nrow(meta_df)))

# Defined color for CN
CN_color <- setNames(metadata(spe)$color_vectors$color_10,unique(meta_df[[cn_column]]))
RFS_color <- setNames(c("#EFC000FF", "#0073C2FF"),unique(meta_df$RFS_status))
Tissue_color <- setNames(metadata(spe)$color_vectors$color_20[1:3],unique(meta_df$Tissue))

# CREATE CONTINGENCY TABLE AND CALCULATE PROPORTIONS
print("Creating contingency table...")

# Create contingency table: rows = CNs, columns = cell types
contingency_table <- table(meta_df[[cn_column]], meta_df$sub_celltype)
print(paste("Contingency table dimensions:", nrow(contingency_table), "CNs x", ncol(contingency_table), "cell types"))

# Show table summary
print("CN distribution:")
print(rowSums(contingency_table))
print("Cell type distribution:")
print(colSums(contingency_table))

# Calculate proportions within each CN (row proportions)
prop_table <- prop.table(contingency_table, margin = 1)

pvalue_adj <- perform_test_for_CN_enrichment(contingency_table)
sig_matrix <- creat_sig_anno_matrix(pvalue_adj)

# CREATE HEATMAPS
# Proportion heatmap with significance annotations
p <- pheatmap(prop_table,
         color = colorRampPalette(c("darkblue","white", "darkred"))(200),
         scale = "column",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = sig_matrix,
         number_color = "black",
         fontsize_number = 12,
         fontsize_row = 10,
         fontsize_col = 9,
         main = paste("CN-Celltype Composition Proportions (", cn_column, ")\n* p<0.05, ** p<0.01, *** p<0.001 (enrichment)"),
         border_color = "grey60")

pdf(file.path(figureDir, paste0("CN_celltype_composition_heatmap_", cn_column, ".pdf")), 
    width = 14, height = 8)
print(p)
dev.off()

# SAVE RESULTS TABLES
print("Saving results tables...")

write.csv(prop_table, file.path(figureDir, paste0("CN_celltype_proportions_", cn_column, ".csv")))
write.csv(pvalue_adj, file.path(figureDir, paste0("CN_celltype_significance_results_", cn_column, ".csv")), row.names = FALSE)

print("Task 1 completed successfully!")

# ============================================================================
# Task 2: CN Composition Stacked Barplots and Pie Charts
# ============================================================================
# Load additional libraries for Task 2
library(viridis)
library(scales)

print("Starting Task 2: CN Composition Visualization")

# DATA PREPARATION FOR VISUALIZATION
print("Preparing data for visualization...")

# Calculate CN composition per patient and ROI
cn_patient_roi <- meta_df %>%
  group_by(patient_id, sample_id, Tissue, Treatment, RFS_status, !!sym(cn_column)) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(patient_id, sample_id) %>%
  mutate(
    total_cells_roi = sum(cell_count),
    cn_fraction = cell_count / total_cells_roi
  ) %>%
  ungroup()

print(paste("Total patient-ROI combinations:", length(unique(paste(cn_patient_roi$patient_id, cn_patient_roi$sample_id)))))
print(paste("CNs found:", length(unique(cn_patient_roi[[cn_column]]))))

# Create patient-ROI identifier for plotting
cn_patient_roi$patient_roi <- cn_patient_roi$sample_id

# Get unique CNs for consistent coloring
unique_cns <- sort(unique(cn_patient_roi[[cn_column]]))
n_cns <- length(unique_cns)

print(paste("Unique CNs for visualization:", paste(unique_cns, collapse = ", ")))

# TASK 2A: STACKED BARPLOT - CN COMPOSITION ACROSS PATIENTS/ROIs
print("Creating stacked barplot for patient/ROI composition...")

# Create stacked barplot
p_stacked_tissue <- ggplot(cn_patient_roi, aes(x = patient_roi, y = cn_fraction, fill = factor(!!sym(cn_column)))) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  facet_wrap(~ Tissue, scales = "free_x", ncol = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    title = paste("CN Composition by Tissue (", cn_column, ")"),
    x = "Patient_ROI", 
    y = "CN Fraction", 
    fill = "CN"
  ) +
  scale_fill_manual(values = CN_color) +
  scale_y_continuous(labels = percent_format(accuracy = 1))

ggsave(file.path(figureDir, paste0("CN_composition_stacked_by_tissue_", cn_column, ".pdf")), 
       p_stacked_tissue, width = 18, height = 10)

# TASK 2B: PIE CHARTS - CN COMPOSITION BY TISSUE
print("Creating pie charts for tissue-specific composition...")

# Calculate CN composition by tissue
cn_tissue_summary <- meta_df %>%
  group_by(Tissue, !!sym(cn_column)) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(Tissue) %>%
  mutate(
    total_cells_tissue = sum(cell_count),
    cn_fraction = cell_count / total_cells_tissue,
    cn_percentage = round(cn_fraction * 100, 1)
  ) %>%
  ungroup()

print("Tissue-specific CN composition calculated")
print(cn_tissue_summary)

# Create pie charts for each tissue
p_pie_combined <- ggplot(cn_tissue_summary, aes(x = "", y = cn_fraction, fill = factor(!!sym(cn_column)))) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 0.5) +
  coord_polar("y", start = 0) +
  facet_wrap(~ Tissue) +
  theme_void() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    strip.text = element_text(size = 14, face = "bold")
  ) +
  labs(
    title = paste("CN Composition by Tissue (", cn_column, ")"),
    fill = "CN"
  ) +
  scale_fill_manual(values = CN_color) +
  geom_text(aes(label = ifelse(cn_percentage >= 8, paste0(cn_percentage, "%"), "")), 
            position = position_stack(vjust = 0.5), 
            color = "white", 
            fontface = "bold", 
            size = 3)

ggsave(file.path(figureDir, paste0("CN_composition_pie_combined_", cn_column, ".pdf")), 
       p_pie_combined, width = 12, height = 8)

# ADDITIONAL VISUALIZATIONS
print("Creating additional composition visualizations...")

# Patient-level CN diversity analysis
patient_cn_diversity <- cn_patient_roi %>%
  group_by(patient_id, Treatment, RFS_status) %>%
  summarise(
    n_cns_present = sum(cn_fraction > 0),
    cn_shannon_diversity = -sum(cn_fraction * log(cn_fraction + 1e-10)),
    dominant_cn = CN_knn_20_cluster_10[which.max(cn_fraction)], ## Replace column
    max_cn_fraction = max(cn_fraction),
    .groups = 'drop'
  )

# CN diversity boxplot
p_diversity <- ggplot(patient_cn_diversity, aes(x = factor(RFS_status, labels = c("No Early Relapse", "Early Relapse")), 
                                                y = cn_shannon_diversity, 
                                                fill = factor(RFS_status))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  labs(
    title = paste("CN Diversity by RFS Status (", cn_column, ")"),
    x = "RFS Status", 
    y = "Shannon Diversity Index"
  ) +
  scale_fill_manual(values = RFS_color) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4)

ggsave(file.path(figureDir, paste0("CN_diversity_by_RFS_", cn_column, ".pdf")), 
       p_diversity, width = 8, height = 6)

print("Saving composition summary data...")

# Save patient-ROI composition data
write.csv(cn_patient_roi, 
          file.path(figureDir, paste0("CN_composition_patient_ROI_", cn_column, ".csv")), 
          row.names = FALSE)

# Save tissue composition summary
write.csv(cn_tissue_summary, 
          file.path(figureDir, paste0("CN_composition_tissue_summary_", cn_column, ".csv")), 
          row.names = FALSE)

print("Task 2 completed successfully!")

# ============================================================================
# Task 3: CN Fraction Calculation Across ROIs
# ============================================================================
print("Starting Task 3: CN Fraction Calculation Across ROIs")

# CALCULATE CN FRACTIONS PER ROI
print("Calculating CN fractions per ROI...")

# Create comprehensive CN fraction data per ROI
cn_roi_fractions <- meta_df %>%
  group_by(patient_id, sample_id, Tissue, Treatment, RFS_status, !!sym(cn_column)) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(patient_id, sample_id) %>%
  mutate(
    total_cells_roi = sum(cell_count),
    cn_fraction = cell_count / total_cells_roi,
    cn_percentage = round(cn_fraction * 100, 2)
  ) %>%
  ungroup()

# Add additional metadata
cn_roi_fractions$patient_roi_id <- cn_roi_fractions$sample_id

print(paste("Total ROI records:", nrow(cn_roi_fractions)))
print(paste("Unique ROIs:", length(unique(cn_roi_fractions$patient_roi_id))))
print(paste("Unique patients:", length(unique(cn_roi_fractions$patient_id))))

# Check data completeness
print("=== DATA COMPLETENESS CHECK ===")
print("ROIs per patient:")
roi_per_patient <- cn_roi_fractions %>%
  group_by(patient_id) %>%
  summarise(n_rois = length(unique(sample_id)), .groups = 'drop')
print(table(roi_per_patient$n_rois))

print("CNs per ROI:")
cn_per_roi <- cn_roi_fractions %>%
  group_by(patient_roi_id) %>%
  summarise(n_cns = n(), .groups = 'drop')
print(table(cn_per_roi$n_cns))

# CREATE WIDE FORMAT FOR ANALYSIS
print("Creating wide format data for analysis...")

# Convert to wide format - each CN becomes a column
cn_fractions_wide <- cn_roi_fractions %>%
  select(patient_id, sample_id, patient_roi_id, Tissue, Treatment, RFS_status, !!sym(cn_column), cn_fraction) %>%
  pivot_wider(
    names_from = !!sym(cn_column),
    values_from = cn_fraction,
    values_fill = 0,
    names_prefix = "CN_"
  )

print(paste("Wide format dimensions:", nrow(cn_fractions_wide), "rows x", ncol(cn_fractions_wide), "columns"))

# Ensure we have all CNs represented
cn_columns <- colnames(cn_fractions_wide)[grepl("^CN_", colnames(cn_fractions_wide))]
print(paste("CN fraction columns:", paste(cn_columns, collapse = ", ")))


# PATIENT-LEVEL AGGREGATION
print("Creating patient-level aggregated data...")

# Aggregate to patient level (mean across ROIs per patient)
cn_patient_level <- cn_fractions_wide %>%
  group_by(patient_id, Treatment, RFS_status) %>%
  summarise(
    n_rois = n(),
    across(all_of(cn_columns), mean, na.rm = TRUE, .names = "{.col}"),
    .groups = 'drop'
  )

print(paste("Patient-level data:", nrow(cn_patient_level), "patients"))

# SUMMARY STATISTICS
print("Calculating summary statistics...")

# Overall CN fraction statistics
cn_summary_stats <- cn_roi_fractions %>%
  group_by(!!sym(cn_column)) %>%
  summarise(
    n_observations = n(),
    n_rois = length(unique(patient_roi_id)),
    n_patients = length(unique(patient_id)),
    mean_fraction = round(mean(cn_fraction), 4),
    median_fraction = round(median(cn_fraction), 4),
    sd_fraction = round(sd(cn_fraction), 4),
    min_fraction = round(min(cn_fraction), 4),
    max_fraction = round(max(cn_fraction), 4),
    q25_fraction = round(quantile(cn_fraction, 0.25), 4),
    q75_fraction = round(quantile(cn_fraction, 0.75), 4),
    .groups = 'drop'
  )

print("=== CN FRACTION SUMMARY STATISTICS ===")
print(cn_summary_stats)

# DISTRIBUTION ANALYSIS

print("Analyzing CN fraction distributions...")

# Create distribution plots for each CN
library(ggplot2)
library(patchwork)

# Prepare data for distribution plotting
cn_dist_data <- cn_roi_fractions %>%
  select(patient_roi_id, Tissue, Treatment, RFS_status, !!sym(cn_column), cn_fraction)

# Get top CNs by abundance for plotting
top_cns <- cn_summary_stats %>%
  arrange(desc(mean_fraction)) %>%
  slice_head(n = 8) %>%
  pull(!!sym(cn_column))

print(paste("Top CNs for distribution analysis:", paste(top_cns, collapse = ", ")))

# Create distribution plots
dist_plots <- list()
for (i in seq_along(top_cns)) {
  cn_val <- top_cns[i]
  cn_data <- cn_dist_data[cn_dist_data[[cn_column]] == cn_val, ]
  
  p <- ggplot(cn_data, aes(x = cn_fraction)) +
    geom_histogram(aes(fill = Tissue), alpha = 0.7, bins = 20, position = "identity") +
    theme_minimal() +
    labs(title = paste("CN", cn_val), x = "Fraction", y = "Count") +
    scale_fill_manual(values = Tissue_color) +
    theme(legend.position = "right", plot.title = element_text(size = 10))
  
  dist_plots[[i]] <- p
}

# Combine distribution plots
p_dist_combined <- wrap_plots(dist_plots, ncol = 4)

ggsave(file.path(figureDir, paste0("CN_fraction_distributions_", cn_column, ".pdf")), 
       p_dist_combined, width = 16, height = 8)

# CORRELATION ANALYSIS
print("Analyzing CN-CN correlations...")

# Calculate correlations between CNs at patient level
cn_correlation_matrix <- cor(cn_patient_level[, cn_columns], use = "complete.obs")

# Create correlation heatmap
library(corrplot)
p <- corrplot(cn_correlation_matrix, 
              method = "color", type = "upper",order = "hclust",
              tl.col = "black", tl.srt = 45, col = colorRampPalette(c("blue", "white", "red"))(100),
              title = paste("CN-CN Correlations (Patient Level) -", cn_column),
              mar = c(0,0,2,0),addCoef.col = "black",number.cex = 0.7)

pdf(file.path(figureDir, paste0("CN_correlation_matrix_", cn_column, ".pdf")), 
    width = 10, height = 8)
print(p)
dev.off()

# SAVE ALL RESULTS
# 1. Main ROI-level data (long format)
write.csv(cn_roi_fractions, 
          file.path(figureDir, paste0("CN_fractions_per_ROI_", cn_column, ".csv")), 
          row.names = FALSE)

# 2. Patient-level aggregated data
write.csv(cn_patient_level, 
          file.path(figureDir, paste0("CN_fractions_patient_level_", cn_column, ".csv")), 
          row.names = FALSE)

# 3. Summary statistics
write.csv(cn_summary_stats, 
          file.path(figureDir, paste0("CN_fraction_summary_statistics_", cn_column, ".csv")), 
          row.names = FALSE)

# 4. Correlation matrix
write.csv(cn_correlation_matrix, 
          file.path(figureDir, paste0("CN_correlation_matrix_", cn_column, ".csv")))

# PREPARE DATA FOR TASKS 4-5 (CLINICAL COMPARISONS)
print("Preparing data for clinical comparisons (Tasks 4-5)...")

clinical_analysis_data <- cn_patient_level

clinical_analysis_data$RFS_label <- factor(clinical_analysis_data$RFS_status, 
                                           levels = c(0, 1), 
                                           labels = c("No Early Relapse", "Early Relapse"))

clinical_analysis_data$Treatment_RFS <- paste(clinical_analysis_data$Treatment, 
                                              clinical_analysis_data$RFS_label, 
                                              sep = "_")

# Save for Tasks 4-5
write.csv(clinical_analysis_data, 
          file.path(figureDir, paste0("CN_data_for_clinical_analysis_", cn_column, ".csv")), 
          row.names = FALSE)

print("Task 3 completed successfully!")

# ============================================================================
# Task 4: Clinical Group Comparison of CN Composition
# ============================================================================
print("Starting Task 4: Clinical Group Comparison of CN Composition")

# LOAD DATA FROM TASK 3
print("Loading data from Task 3...")

if (!exists("clinical_analysis_data")) {
  # If data doesn't exist in environment, load from file
  clinical_analysis_data <- read.csv(file.path(figurePath, paste0("CN_data_for_clinical_analysis_", cn_column, ".csv")))
  print("Loaded clinical analysis data from file")
}

print(paste("Clinical analysis data:", nrow(clinical_analysis_data), "patients"))

# Get CN column names
cn_columns <- colnames(clinical_analysis_data)[grepl("^CN_", colnames(clinical_analysis_data))]
print(paste("CN columns for analysis:", paste(cn_columns, collapse = ", ")))

# Check clinical group distribution
print("=== CLINICAL GROUP DISTRIBUTION ===")
clinical_dist <- table(clinical_analysis_data$Treatment, clinical_analysis_data$RFS_status)
print("Treatment × RFS Status:")
print(clinical_dist)


# TASK 4A: OVERALL RFS STATUS COMPARISON (NOT CONSIDERING TREATMENT)
print("Task 4A: Overall RFS Status Comparison")

# Prepare data for statistical testing
rfs_comparison_results <- data.frame()

# Perform Wilcoxon test for each CN
for (cn_col in cn_columns) {
  cn_data <- clinical_analysis_data[, c("patient_id", "RFS_status", "RFS_label", cn_col)]
  colnames(cn_data)[4] <- "cn_value"
  
  # Remove any missing values
  cn_data <- cn_data[!is.na(cn_data$cn_value), ]
  
  if (nrow(cn_data) < 4) {
    print(paste("Skipping", cn_col, "- insufficient data"))
    next
  }
  
  # Check if we have both RFS groups
  rfs_groups <- unique(cn_data$RFS_status)
  if (length(rfs_groups) < 2) {
    print(paste("Skipping", cn_col, "- only one RFS group"))
    next
  }
  
  # Calculate group statistics
  group_stats <- cn_data %>%
    group_by(RFS_status, RFS_label) %>%
    summarise(
      n = n(),
      mean = mean(cn_value, na.rm = TRUE),
      median = median(cn_value, na.rm = TRUE),
      sd = sd(cn_value, na.rm = TRUE),
      q25 = quantile(cn_value, 0.25, na.rm = TRUE),
      q75 = quantile(cn_value, 0.75, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Perform Wilcoxon test
  wilcox_result <- wilcox.test(cn_value ~ RFS_status, data = cn_data)
  
  # Calculate effect size (rank-biserial correlation)
  n1 <- sum(cn_data$RFS_status == 0)
  n2 <- sum(cn_data$RFS_status == 1)
  U <- wilcox_result$statistic
  effect_size <- 1 - (2 * U) / (n1 * n2)
  
  # Store results
  rfs_comparison_results <- rbind(rfs_comparison_results, data.frame(
    CN = gsub("CN_", "", cn_col),
    comparison_type = "Overall_RFS",
    treatment_group = "All",
    n_no_relapse = group_stats$n[group_stats$RFS_status == 0],
    n_early_relapse = group_stats$n[group_stats$RFS_status == 1],
    mean_no_relapse = round(group_stats$mean[group_stats$RFS_status == 0], 4),
    mean_early_relapse = round(group_stats$mean[group_stats$RFS_status == 1], 4),
    median_no_relapse = round(group_stats$median[group_stats$RFS_status == 0], 4),
    median_early_relapse = round(group_stats$median[group_stats$RFS_status == 1], 4),
    p_value = wilcox_result$p.value,
    effect_size = effect_size,
    test_statistic = wilcox_result$statistic,
    stringsAsFactors = FALSE
  ))
}

print(paste("Overall RFS comparisons completed:", nrow(rfs_comparison_results)))


# TASK 4B: TREATMENT-STRATIFIED COMPARISON
print("Task 4B: Treatment-Stratified RFS Comparison")

# Get available treatments
treatments <- unique(clinical_analysis_data$Treatment)
print(paste("Treatments for analysis:", paste(treatments, collapse = ", ")))

treatment_comparison_results <- data.frame()

# Perform analysis within each treatment group
for (treatment in treatments) {
  print(paste("Analyzing treatment group:", treatment))
  
  treatment_data <- clinical_analysis_data[clinical_analysis_data$Treatment == treatment, ]
  
  # Check if we have both RFS groups in this treatment
  rfs_groups_treatment <- unique(treatment_data$RFS_status)
  if (length(rfs_groups_treatment) < 2) {
    print(paste("Skipping", treatment, "- only one RFS group"))
    next
  }
  
  print(paste("  Patients in", treatment, "group:", nrow(treatment_data)))
  print(paste("  RFS distribution:", paste(table(treatment_data$RFS_status), collapse = " / ")))
  
  # Test each CN within this treatment group
  for (cn_col in cn_columns) {
    cn_data <- treatment_data[, c("patient_id", "RFS_status", "RFS_label", cn_col)]
    colnames(cn_data)[4] <- "cn_value"
    
    # Remove missing values
    cn_data <- cn_data[!is.na(cn_data$cn_value), ]
    
    if (nrow(cn_data) < 4) next
    
    # Check RFS groups again after filtering
    if (length(unique(cn_data$RFS_status)) < 2) next
    
    # Calculate group statistics
    group_stats <- cn_data %>%
      group_by(RFS_status, RFS_label) %>%
      summarise(
        n = n(),
        mean = mean(cn_value, na.rm = TRUE),
        median = median(cn_value, na.rm = TRUE),
        sd = sd(cn_value, na.rm = TRUE),
        q25 = quantile(cn_value, 0.25, na.rm = TRUE),
        q75 = quantile(cn_value, 0.75, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Perform Wilcoxon test
    wilcox_result <- wilcox.test(cn_value ~ RFS_status, data = cn_data)
    
    # Calculate effect size
    n1 <- sum(cn_data$RFS_status == 0)
    n2 <- sum(cn_data$RFS_status == 1)
    if (n1 > 0 && n2 > 0) {
      U <- wilcox_result$statistic
      effect_size <- 1 - (2 * U) / (n1 * n2)
    } else {
      effect_size <- NA
    }
    
    # Store results
    treatment_comparison_results <- rbind(treatment_comparison_results, data.frame(
      CN = gsub("CN_", "", cn_col),
      comparison_type = "Treatment_Stratified",
      treatment_group = treatment,
      n_no_relapse = ifelse(0 %in% group_stats$RFS_status, group_stats$n[group_stats$RFS_status == 0], 0),
      n_early_relapse = ifelse(1 %in% group_stats$RFS_status, group_stats$n[group_stats$RFS_status == 1], 0),
      mean_no_relapse = ifelse(0 %in% group_stats$RFS_status, round(group_stats$mean[group_stats$RFS_status == 0], 4), NA),
      mean_early_relapse = ifelse(1 %in% group_stats$RFS_status, round(group_stats$mean[group_stats$RFS_status == 1], 4), NA),
      median_no_relapse = ifelse(0 %in% group_stats$RFS_status, round(group_stats$median[group_stats$RFS_status == 0], 4), NA),
      median_early_relapse = ifelse(1 %in% group_stats$RFS_status, round(group_stats$median[group_stats$RFS_status == 1], 4), NA),
      p_value = wilcox_result$p.value,
      effect_size = effect_size,
      test_statistic = wilcox_result$statistic,
      stringsAsFactors = FALSE
    ))
  }
}

print(paste("Treatment-stratified comparisons completed:", nrow(treatment_comparison_results)))

# COMBINE RESULTS AND APPLY MULTIPLE TESTING CORRECTION
print("Combining results and applying multiple testing correction...")

# Combine all results
all_comparison_results <- rbind(rfs_comparison_results, treatment_comparison_results)

# Apply Benjamini-Hochberg correction within each comparison type
all_comparison_results <- all_comparison_results %>%
  group_by(comparison_type) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    .groups = 'drop'
  )

# Add significance classifications
all_comparison_results$significance <- ifelse(all_comparison_results$p_adjusted < 0.001, "***",
                                              ifelse(all_comparison_results$p_adjusted < 0.01, "**",
                                                     ifelse(all_comparison_results$p_adjusted < 0.05, "*", "ns")))

all_comparison_results$significant <- all_comparison_results$p_adjusted <= 0.05
print(paste("Total comparisons performed:", nrow(all_comparison_results)))

# EFFECT SIZE INTERPRETATION
print("Adding effect size interpretation...")

# Add effect size interpretation (for rank-biserial correlation)
all_comparison_results$effect_size_interpretation <- ifelse(
  is.na(all_comparison_results$effect_size), "Unknown",
  ifelse(abs(all_comparison_results$effect_size) < 0.1, "Negligible",
         ifelse(abs(all_comparison_results$effect_size) < 0.3, "Small",
                ifelse(abs(all_comparison_results$effect_size) < 0.5, "Medium", "Large")))
)

# Add direction of effect
all_comparison_results$effect_direction <- ifelse(
  is.na(all_comparison_results$mean_no_relapse) | is.na(all_comparison_results$mean_early_relapse), "Unknown",
  ifelse(all_comparison_results$mean_early_relapse > all_comparison_results$mean_no_relapse, 
         "Higher_in_Early_Relapse", "Higher_in_No_Relapse")
)

# IDENTIFY SIGNIFICANT FINDINGS
print("Identifying significant findings...")

# Overall significant findings
overall_significant <- all_comparison_results[
  all_comparison_results$comparison_type == "Overall_RFS" & all_comparison_results$significant, 
]

if (nrow(overall_significant) > 0) {
  print("=== SIGNIFICANT OVERALL RFS FINDINGS ===")
  overall_sig_summary <- overall_significant[, c("CN", "p_value", "p_adjusted", "significance", 
                                                 "effect_size", "effect_direction")]
  print(overall_sig_summary)
} else {
  print("No significant overall RFS differences found")
}

# Treatment-specific significant findings
treatment_significant <- all_comparison_results[
  all_comparison_results$comparison_type == "Treatment_Stratified" & all_comparison_results$significant, 
]

if (nrow(treatment_significant) > 0) {
  print("=== SIGNIFICANT TREATMENT-STRATIFIED FINDINGS ===")
  treatment_sig_summary <- treatment_significant[, c("CN", "treatment_group", "p_value", "p_adjusted", 
                                                     "significance", "effect_size", "effect_direction")]
  print(treatment_sig_summary)
} else {
  print("No significant treatment-stratified differences found")
}

# CN-specific summary (which CNs show most differences)
cn_summary <- all_comparison_results %>%
  group_by(CN) %>%
  summarise(
    total_tests = n(),
    significant_tests = sum(significant),
    prop_significant = round(mean(significant), 3),
    min_p_value = round(min(p_adjusted), 4),
    max_effect_size = round(max(abs(effect_size), na.rm = TRUE), 3),
    .groups = 'drop'
  ) %>%
  arrange(desc(prop_significant), min_p_value)

print("=== CN-SPECIFIC SUMMARY (TOP DIFFERENTIAL CNs) ===")
print(head(cn_summary, 10))

# ============================================================================
# PREPARE DATA FOR TASK 5 (BOXPLOTS)
# ============================================================================

print("Preparing data for Task 5 boxplot visualization...")

# Create plotting data with proper factors
plotting_data <- clinical_analysis_data
plotting_data$RFS_label <- factor(plotting_data$RFS_status, 
                                  levels = c(0, 1), 
                                  labels = c("No Early Relapse", "Early Relapse"))

# Convert CN columns to long format for easier plotting
cn_long_data <- plotting_data %>%
  select(patient_id, RFS_status, RFS_label, Treatment, all_of(cn_columns)) %>%
  pivot_longer(cols = all_of(cn_columns), names_to = "CN", values_to = "CN_fraction") %>%
  mutate(CN = gsub("CN_", "", CN))

# Add significance information to plotting data
cn_long_data <- cn_long_data %>%
  left_join(
    all_comparison_results %>% 
      select(CN, comparison_type, treatment_group, p_adjusted, significance, significant),
    by = c("CN" = "CN"),
    relationship = "many-to-many"
  )

# Save plotting data for Task 5
write.csv(cn_long_data, 
          file.path(figurePath, paste0("CN_data_for_boxplots_", cn_column, ".csv")), 
          row.names = FALSE)

# ============================================================================
# Task 5: Boxplot Visualization with Significance Testing
# ============================================================================
print("Starting Task 5: Boxplot Visualization with Significance Testing")

# Load additional libraries for advanced plotting
library(ggpubr)
library(patchwork)

# LOAD DATA FROM TASK 4
print("Loading data from Task 4...")

# Load plotting data
if (!exists("cn_long_data")) {
  cn_long_data <- read.csv(file.path(figurePath, paste0("CN_data_for_boxplots_", cn_column, ".csv")))
  print("Loaded plotting data from file")
}

# Load statistical results
if (!exists("all_comparison_results")) {
  all_comparison_results <- read.csv(file.path(figurePath, paste0("CN_clinical_comparison_statistics_", cn_column, ".csv")))
  print("Loaded statistical results from file")
}

print(paste("Plotting data:", nrow(cn_long_data), "observations"))
print(paste("Statistical results:", nrow(all_comparison_results), "comparisons"))

# Create boxplot output directory
boxplot_dir <- file.path(figurePath, "boxplots")
if (!dir.exists(boxplot_dir)) {
  dir.create(boxplot_dir, recursive = TRUE)
}

# Get unique CNs
unique_cns <- sort(as.numeric(unique(cn_long_data$CN)))
print(paste("CNs for visualization:", paste(unique_cns, collapse = ", ")))

# TASK 5A: INDIVIDUAL CN BOXPLOTS - OVERALL RFS COMPARISON
print("Task 5A: Creating individual CN boxplots - Overall RFS comparison")

# Create individual boxplots for each CN (overall RFS comparison)
cn_long_data$RFS_status <- as.factor(cn_long_data$RFS_status)
for (cn_val in unique_cns) {
  print(paste("Creating overall RFS boxplot for CN", cn_val))
  
  # Filter data for this CN
  cn_data <- cn_long_data[cn_long_data$CN == cn_val, ]
  
  # Get statistical result for this CN (overall comparison)
  stat_result <- all_comparison_results[
    all_comparison_results$CN == cn_val & 
      all_comparison_results$comparison_type == "Overall_RFS", 
  ]
  
  if (nrow(stat_result) == 0) {
    print(paste("No statistical result found for CN", cn_val))
    next
  }
  
  # Prepare significance annotation
  p_val <- stat_result$p_adjusted[1]
  sig_label <- stat_result$significance[1]
  effect_direction <- stat_result$effect_direction[1]
  
  # Create subtitle with statistics
  subtitle_text <- paste0(
    "Wilcoxon p = ", format(p_val, scientific = TRUE, digits = 3),
    " (adjusted), Effect: ", effect_direction
  )
  
  # Create boxplot
  cn_data$RFS_status <- as.factor(cn_data$RFS_status)
  p_overall <- ggplot(cn_data, aes(x = RFS_label, y = CN_fraction, fill = RFS_status)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    scale_fill_manual(values = RFS_color) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "none",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste("CN", cn_val, "Fraction - Overall RFS Comparison"),
      subtitle = subtitle_text,
      x = "RFS Status", 
      y = "CN Fraction"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))
  
  # Add significance annotation
  y_max <- max(cn_data$CN_fraction, na.rm = TRUE)
  y_pos <- y_max * 1.1
  
  p_overall <- p_overall + 
    stat_compare_means(method = "wilcox.test", label = "p.signif", 
                       hide.ns = FALSE, size = 3)
  
  # Save individual boxplot
  ggsave(file.path(boxplot_dir, paste0("CN", cn_val, "_RFS_overall_boxplot_", cn_column, ".pdf")), 
         p_overall, width = 6, height = 5)
}

# ============================================================================
# TASK 5B: INDIVIDUAL CN BOXPLOTS - TREATMENT-STRATIFIED COMPARISON
# ============================================================================
print("Task 5B: Creating individual CN boxplots - Treatment-stratified comparison")

# Get available treatments
treatments <- unique(cn_long_data$Treatment)
print(paste("Treatments:", paste(treatments, collapse = ", ")))

# Create treatment-stratified boxplots for each CN
for (cn_val in unique_cns) {
  print(paste("Creating treatment-stratified boxplots for CN", cn_val))
  
  treatment_plots <- list()
  
  for (treatment in treatments) {
    # Filter data for this CN and treatment
    cn_treatment_data <- cn_long_data[
      cn_long_data$CN == cn_val & cn_long_data$Treatment == treatment, 
    ]
    
    if (nrow(cn_treatment_data) < 4) {
      print(paste("Insufficient data for CN", cn_val, "in", treatment))
      next
    }
    
    # Check if we have both RFS groups
    if (length(unique(cn_treatment_data$RFS_status)) < 2) {
      print(paste("Only one RFS group for CN", cn_val, "in", treatment))
      next
    }
    
    # Get statistical result
    stat_result <- all_comparison_results[
      all_comparison_results$CN == cn_val & 
        all_comparison_results$comparison_type == "Treatment_Stratified" &
        all_comparison_results$treatment_group == treatment, 
    ]
    
    if (nrow(stat_result) == 0) {
      print(paste("No statistical result for CN", cn_val, "in", treatment))
      next
    }
    
    # Prepare significance annotation
    p_val <- stat_result$p_adjusted[1]
    sig_label <- stat_result$significance[1]
    effect_direction <- stat_result$effect_direction[1]
    
    # Create subtitle
    subtitle_text <- paste0(
      "Wilcoxon p = ", format(p_val, scientific = TRUE, digits = 3),
      " (adj.)"
    )
    
    # Create boxplot
    p_treatment <- ggplot(cn_treatment_data, aes(x = RFS_label, y = CN_fraction, fill = RFS_status)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
      scale_fill_manual(values = RFS_color) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 10),
        legend.position = "none",
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = paste("CN", cn_val, "-", treatment, "Treatment"),
        subtitle = subtitle_text,
        x = "RFS Status", 
        y = "CN Fraction"
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))
    
    # Add significance annotation
    y_max <- max(cn_treatment_data$CN_fraction, na.rm = TRUE)
    y_pos <- y_max * 1.1
    
    p_treatment <- p_treatment +
      stat_compare_means(method = "wilcox.test", label = "p.signif", 
                         hide.ns = FALSE, size = 3)
    
    treatment_plots[[treatment]] <- p_treatment
  }
  
  # Combine treatment plots if we have both treatments
  if (length(treatment_plots) > 0) {
    if (length(treatment_plots) == 2) {
      combined_plot <- treatment_plots[[1]] + treatment_plots[[2]]
      width_val <- 12
    } else {
      combined_plot <- treatment_plots[[1]]
      width_val <- 6
    }
    
    ggsave(file.path(boxplot_dir, paste0("CN", cn_val, "_RFS_treatment_stratified_", cn_column, ".pdf")), 
           combined_plot, width = width_val, height = 5)
  }
}

# ============================================================================
# TASK 5C: SUMMARY BOXPLOTS - ALL CNs TOGETHER
# ============================================================================
print("Task 5C: Creating summary boxplots - All CNs together")

# Overall RFS comparison - all CNs in one plot
p_summary_overall <- ggplot(cn_long_data, aes(x = RFS_label, y = CN_fraction, fill = RFS_status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.5) +
  facet_wrap(~ paste("CN", CN), scales = "free_y", ncol = 4) +
  scale_fill_manual(values = RFS_color) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "CN Fraction Comparison - Overall RFS Status",
    x = "RFS Status", 
    y = "CN Fraction", 
    fill = "RFS Status"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

# Add significance annotations using stat_compare_means
p_summary_overall <- p_summary_overall +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     hide.ns = FALSE, size = 3)

ggsave(file.path(boxplot_dir, paste0("CN_all_RFS_overall_summary_", cn_column, ".pdf")), 
       p_summary_overall, width = 16, height = 12)

# Treatment-stratified comparison - all CNs
p_summary_treatment <- ggplot(cn_long_data, aes(x = RFS_label, y = CN_fraction, fill = RFS_status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.4) +
  facet_grid(Treatment ~ paste("CN", CN), scales = "free_y") +
  scale_fill_manual(values = RFS_color) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "CN Fraction Comparison - Treatment Stratified",
    x = "RFS Status", 
    y = "CN Fraction", 
    fill = "RFS Status"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

# Add significance annotations
p_summary_treatment <- p_summary_treatment +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     hide.ns = FALSE, size = 2.5)

ggsave(file.path(boxplot_dir, paste0("CN_all_RFS_treatment_summary_", cn_column, ".pdf")), 
       p_summary_treatment, width = 20, height = 10)

# ============================================================================
# TASK 5E: EFFECT SIZE VISUALIZATION
# ============================================================================

print("Task 5E: Creating effect size visualization")

# Create effect size plot
effect_size_data <- all_comparison_results %>%
  filter(!is.na(effect_size)) %>%
  mutate(
    abs_effect_size = abs(effect_size),
    comparison_label = paste(comparison_type, treatment_group, sep = " - ")
  )

if (nrow(effect_size_data) > 0) {
  p_effect_size <- ggplot(effect_size_data, aes(x = reorder(paste("CN", CN), abs_effect_size), 
                                                y = abs_effect_size, 
                                                fill = effect_size_interpretation)) +
    geom_col(alpha = 0.8) +
    facet_wrap(~ comparison_label, scales = "free_x") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold")
    ) +
    labs(
      title = "Effect Sizes of CN Differences",
      x = "CN", 
      y = "Absolute Effect Size (Rank-Biserial Correlation)", 
      fill = "Effect Size"
    ) +
    scale_fill_manual(values = c("Negligible" = "grey80", "Small" = "lightblue", 
                                 "Medium" = "orange", "Large" = "red")) +
    geom_hline(yintercept = c(0.1, 0.3, 0.5), linetype = "dashed", alpha = 0.5)
  
  ggsave(file.path(boxplot_dir, paste0("CN_effect_sizes_", cn_column, ".pdf")), 
         p_effect_size, width = 14, height = 8)
}










##### Cellular neighborhood analysis
# Choose K clusters
k_clusters <- c(10, 15, 20)
Pairnames <- colPairNames(spe)

for(pairname_ in Pairnames){
  
  ## aggregate neighbor
  spe <- aggregateNeighbors(
    spe,
    colPairName = pairname_,
    aggregate_by = "metadata",
    count_by = "sub_celltype"
  )
  
  ## Iterative on number of clusters
  for(k_cluster_ in k_clusters){
    ## create folder
    figureDir_ <- file.path(figureDir,paste0(pairname_," ",k_cluster_))
    if(!file.exists(figureDir_)){
      dir.create(figureDir_,recursive = T)
    }
    
    cn <- kmeans(spe$aggregatedNeighbors, centers = k_cluster_)
    
    # 创建动态列名避免覆盖
    cluster_colname <- paste0(pairname_,"_cluster_", k_cluster_)
    colData(spe)[[cluster_colname]] <- as.factor(cn$cluster)
    spe$cn_celltypes <- as.factor(colData(spe)[[cluster_colname]])
    
    # 绘制热图
    for_plot <- prop.table(table(as.character(spe$cn_celltypes), spe$sub_celltype), margin = 1)
    p <- pheatmap(for_plot, 
                  color = colorRampPalette(c("dark blue", "white", "dark red"))(100),
                  scale = "column")
    
    pdf(file.path(figureDir_, paste0("CN analysis of ",pairname_," in ",k_cluster_," cluster.pdf")),width = 15,height = 10)
    print(p)
    dev.off()
    
    # 绘制组成条形图 
    if(T){
      # get the number of celltypes within different CN
      countdf_ <- Transform_CellCountMat(
        spe_ = spe,
        clinicalFeatures = "cn_celltypes",
        img_id = img_id_,
        count_by = "sub_celltype",  # 这里保持使用临时列
        is.fraction = FALSE
      )
      
      # Identify cell subpopulation columns by excluding "PID" and "cn_celltypes"
      cell_subpop_cols <- setdiff(names(countdf_), c("PID", "cn_celltypes"))
      
      # Sum cell subpopulation counts for each cn_celltypes
      summed_df <- countdf_ %>%
        group_by(cn_celltypes) %>%
        summarise(across(all_of(cell_subpop_cols), sum, na.rm = TRUE)) %>%
        ungroup()
      
      # Convert cn_celltypes to factor for categorical plotting
      summed_df$cn_celltypes <- factor(summed_df$cn_celltypes)
      
      # Transform to long format for ggplot2
      long_df <- summed_df %>%
        pivot_longer(cols = all_of(cell_subpop_cols), 
                     names_to = "cell_subpopulation", 
                     values_to = "count")
      
      long_df <- as.data.frame(long_df)
      long_df$cell_subpopulation <- as.factor(long_df$cell_subpopulation)
      
      # Create the stacked barplot with horizontal bars
      p <- ggplot(long_df, aes(x = cn_celltypes, y = count, fill = cell_subpopulation)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(x = "CN Cell Types", 
             y = "Total Cell Count", 
             fill = "Cell Subpopulation") +
        theme_minimal() +
        scale_fill_manual(values = metadata(spe)$color_vectors[["sub_celltype"]]) +
        coord_flip()
      
      pdf(file.path(figureDir_, paste0("CNAtype_stacked_barplot_k", k_cluster_, ".pdf")),
          width = 10,
          height = 7.5)
      print(p)
      dev.off()
    }
    
    # Save Matrix
    TRG_countdf <- Transform_CellCountMat(spe_ = spe,clinicalFeatures = c(names(compare_groups)),img_id = img_id_,count_by = "cn_celltypes",is.fraction = TRUE)
    TRG_countdf <- as.data.frame(TRG_countdf)
    
    write.csv(TRG_countdf,file.path(figureDir_, paste0("CN analysis of ",pairname_," in ",k_cluster_," cluster.csv")))
    
    # Compare within different clinical matrix
    for(compare_group_ in names(compare_groups)){
      plotdf <- TRG_countdf[,c(1:k_cluster_,match(compare_group_,colnames(TRG_countdf)))]
      
      # Data transform
      plotdf <- pivot_longer(data = plotdf,cols = 1:k_cluster_,names_to = "CN_type",values_to = "Fraction")
      plotdf <- as.data.frame(plotdf)
      plotdf <- na.omit(plotdf)
      colnames(plotdf)[1] <- "group"
      
      plotdf$Fraction <- as.numeric(plotdf$Fraction)
      plotdf$group <- as.factor(plotdf$group)
      
      # Plot boxplot
      p1 <- ggplot(plotdf, aes(x = group, y = Fraction, fill = group)) +
        geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
        geom_jitter(color = "darkgrey", position = position_jitter(width = 0.2), size = 0.1, alpha = 0.3) +
        scale_fill_manual(values = ggsci::pal_jco("default")(length(unique(plotdf$group)))) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 12, face = "bold"),
          text = element_text(size = 12),
          axis.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(size = 12),
          strip.background = element_blank()
        ) +
        ggtitle(paste0("CN analysis of ", pairname_, " in ", compare_group_)) +
        facet_wrap(~CN_type,scales = "free_y") +
        stat_compare_means(
          comparisons = compare_groups[[clinical.info]], 
          method = "wilcox.test",        # More robust than t-test
          p.adjust.method = "fdr",       # Less conservative than Bonferroni
          label = "p.adj",               # Shows adjusted p-values
          hide.ns = FALSE                 # Cleaner visualization
        )
      
      pdf(file.path(figureDir_, paste0("CN analysis of ",compare_group_," in ",pairname_," with ",k_cluster_," cluster.pdf")),width = 12,height = 15)
      print(p1)
      dev.off()
    }
  }
}

## Visualize the CN types on images
pairname_ <- "knn_20"
all_images <- unique(spe$sample_id)

num_figures <- 10
length_ <- length(all_images) %/% num_figures
for(i in 1:num_figures){
  sample_image <- all_images[c((length_*(i-1)):(length_*(i)))]
  
  #print(sample_image)
  spe_subset <- spe[,spe$sample_id %in% sample_image]
  
  ## Plot
  p <- plotSpatial(spe_subset,  # spe
                   node_color_by = "knn_20_cluster_10", img_id = "sample_id",node_size_fix = 0.15) + 
    scale_color_brewer(palette = "Set3")
  
  pdf(file.path(figureDir, paste0("CN Types on select images of figure id ",i,".pdf")),width = 12,height = 9)
  print(p)
  dev.off()
}

## Cell-cell pairwise interaction analysis
out <- readRDS(file = file.path(saveDir, "Interaction_analysis_out.rds"))
out <- out[out$group_by %in% spe$sample_id, ] ## subset tumor interaction

out$group <- spe$patient_id[match(out$group_by, spe$sample_id)]
head(out)

## overall interaction
df_ <- out %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE))
p <- df_ %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  geom_text(aes(from_label, to_label, label = sum_sigval), size = 2) +
  scale_fill_gradient2(low = muted("#7aa6dcff"),
                       mid = "white",
                       high = muted("#cd534cff")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(figureDir, "celltype pair interaction heatmap in all images.pdf"),width = 12,height = 9)
print(p)
dev.off()

# Cellular interation between different clinical matrix
for(clinical_ in names(compare_groups)){
  clinical_groups <- unique(unlist(compare_groups[[clinical_]]) )
  
  for(clinical_group_ in clinical_groups){
    idx_ <- colData(spe)[,clinical_] %in% clinical_group_
    select_roi_ <- unique(spe[,idx_]$sample_id)
    
    out_ <- out[out$group_by %in% select_roi_, ] %>% 
      as_tibble() %>%
      group_by(from_label, to_label) %>%
      summarize(sum_sigval = sum(sigval, na.rm = TRUE))
    
    p <- out_ %>%
      ggplot() +
      geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
      geom_text(aes(from_label, to_label, label = sum_sigval), size = 2) +
      scale_fill_gradient2(low = muted("#7aa6dcff"),
                           mid = "white",
                           high = muted("#cd534cff")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    pdf(file.path(figureDir,paste0("Pair interaction in ",clinical_,"_group_",clinical_group_,".pdf")) ,width = 12,height = 9)
    print(p)
    dev.off()
  }
}

## Display specific cell-cell interaction
print(unique(out$from_label))
sourceType <- "YAP+ Fibro"
group_ <- "Group_A"

out_ <- out[out$from_label %in% sourceType,]
out_ <- as.data.frame(out_[out_$group %in% group_,]) 

if(T){
  # Step 1: Sum sigval for each to_label across all ROIs
  sum_sigval <- out_ %>%
    group_by(to_label) %>%
    summarise(total_sigval = sum(sigval, na.rm = TRUE)) %>%
    ungroup()
  
  # Step 2: Compute the fraction
  total_rois <- length(unique(out_$group_by))
  sum_sigval <- sum_sigval %>%
    mutate(fraction = total_sigval / total_rois)
  
  # Step 3: Sort by fraction in decreasing order
  sum_sigval <- sum_sigval %>%
    arrange(desc(fraction))
  
  # Step 4: Set factor levels and add sign column
  sum_sigval$to_label <- factor(sum_sigval$to_label, levels = sum_sigval$to_label)
  sum_sigval <- sum_sigval %>%
    mutate(sign = ifelse(fraction > 0, "Positive", "Negative"))
  
  # Step 5: Create the barplot
  p <- ggplot(sum_sigval, aes(x = to_label, y = fraction, fill = sign)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Positive" = "#cd534cff", "Negative" = "#7aa6dcff")) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = "To Label", y = "Fraction", title = "Fraction of Significant Cell-Cell Interactions")
  
  pdf(file.path(figureDir,paste0("celltype pair interaction barplot in ",group_,".pdf")) ,width = 8,height = 6)
  print(p)
  dev.off()
}

# ------------------------------------------------------------------------------------------------
##### Patch Analysis
colnames(colData(spe)) ## view column of data
patch_types <- colnames(colData(spe))[endsWith(colnames(colData(spe)),suffix = "_patch")] ## get all patch type
patch_names <- sapply(patch_types, function(x){
  return(strsplit(x,"_")[[1]][1])
})
all_images <- unique(spe$sample_id)

figureDir <- file.path(workDir, "figures", "4_Patch")
if(!file.exists(figureDir)){
  dir.create(figureDir,recursive = T)
}

## Plot different types of patches
custom_palette <- c(
  "EndoMT_patch" = "#A65628",
  "Endothelial_patch" = "#4DAF4A"
)

## Plot Patch on spatial
for(patch_type_ in patch_types){
  
  ## Sample for plotting
  sample_rois_ <- sample(all_images,12)
  spe_temp <- spe[,spe$sample_id %in% sample_rois_]
  
  patch_type_vec_ <- colData(spe_temp)[,patch_type_] ## get correspond patch type
  colData(spe_temp)[,patch_type_] <- ifelse(is.na(patch_type_vec_),NA,patch_type_)
  
  p <- plotSpatial(
    spe_temp,
    node_color_by = patch_type_,
    img_id = "sample_id",
    node_size_fix = 0.5
  ) +
    theme(legend.position = "right") +
    scale_color_manual(values = custom_palette)
  
  rm(spe_temp,sample_rois_,patch_type_vec_)
  gc()
  
  pdf(file.path(figureDir,paste0(patch_type_," in random ROIs.pdf")), width = 24,height = 18)
  print(p)
  dev.off()
}


### Analysis configuration in different patch types
analysisDF <- as.data.frame(colData(spe)[,c("sample_id","sub_celltype",names(compare_groups),patch_types,unname(patch_names))])
analysisDF_Back <- analysisDF

for(patch_type_ in patch_types){
  
  ## Get patch name
  patch_name_ <- strsplit(patch_type_,"_")[[1]]
  patch_name_ <- paste(patch_name_[1:length(patch_name_)-1],collapse =  "_")
  
  ## Get metadata
  analysisDF <- analysisDF[,match(c("sample_id","sub_celltype",names(compare_groups),patch_type_,patch_name_),
                                  colnames(analysisDF))]
  colnames(analysisDF)[c(11,12)] <- c("patch_type","patch_name")
  
  # --- 1. Patch Frequency and Size Analysis ---
  if(T){
    # Step 1: Calculate total number of cells per sample_id
    total_cells <- analysisDF %>%
      group_by(sample_id) %>%
      summarise(total_cells = n(), .groups = 'drop')
    
    # Step 2: Calculate number of cells with non-NA patch_type per sample_id
    patch_cells <- analysisDF %>%
      filter(!is.na(patch_type)) %>%
      group_by(sample_id) %>%
      summarise(patch_cells = n(), .groups = 'drop')
    
    # Step 3: Combine and calculate the ratio (patch_cells / total_cells)
    cell_ratio <- total_cells %>%
      left_join(patch_cells, by = "sample_id") %>%
      mutate(
        patch_cells = ifelse(is.na(patch_cells), 0, patch_cells),
        ratio = patch_cells / total_cells
      ) %>%
      arrange(desc(ratio))  # Sort in decreasing order
    
    # Print summary statistics
    cat("Summary of cell counts and ratios:\n")
    print(head(cell_ratio, 10))
    cat("\nTotal samples:", nrow(cell_ratio), "\n")
    cat("Samples with patch_type data:", sum(cell_ratio$patch_cells > 0), "\n")
    
    # Step 4 & 5: Create decreasing bar plot with top-5 sample labels
    # Identify top 5 samples with non-zero ratios for labeling
    top5_samples <- cell_ratio %>%
      filter(ratio > 0) %>%
      slice_head(n = 10) %>%
      pull(sample_id)
    
    # Create the lollipop plot
    p <- ggplot(cell_ratio, aes(x = reorder(sample_id, ratio), y = ratio)) +
      # Add the lollipop stick
      geom_segment(aes(xend = reorder(sample_id, ratio), yend = 0), 
                   color = "grey85", linewidth = 0.3) +
      # Add the lollipop circle
      geom_point(color = "#ffccccff", size = 4) +
      # Add text labels inside circles for top 5 samples
      geom_text(
        data = cell_ratio %>% filter(sample_id %in% top5_samples),
        #data = cell_ratio,
        aes(label = round(ratio,digits = 4)),
        color = "grey30",
        size = 2.5,
        fontface = "bold"
      ) +
      coord_flip() +  # Flip coordinates for better readability
      labs(
        title = "Proportion of Cells with Patch Type Information by Sample",
        x = "Sample ID",
        y = "Proportion of Cells with Patch Type",
        caption = paste("Total samples:", nrow(cell_ratio), "| Total cells:", sum(cell_ratio$total_cells))
      ) +
      theme_classic2() +
      theme(
        # axis.text.y = element_blank(),  # Hide y-axis labels 
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
        panel.grid.minor.y = element_blank()
      )
    
    # Display the plot
    pdf(file = file.path(figureDir,paste0(patch_type_," Patch size distribution.pdf")),width = 6,height = 12)
    print(p)
    dev.off()
    
  }
  patch_freq_roi <- analysisDF %>%
    filter(!is.na(patch_type)) %>%                  # Filter cells in patches
    group_by(sample_id) %>%
    summarise(n_patches = n_distinct(patch_type)) %>%  # Count unique patches
    left_join(analysisDF %>% select(sample_id, names(compare_groups)) %>% distinct(), by = "sample_id")
  
  patch_freq_roi$ratio <- cell_ratio[match(patch_freq_roi$sample_id,cell_ratio$sample_id),]$ratio 
  write.csv(patch_freq_roi,file.path(figureDir,paste0(patch_type_," Patch number distribution(ROI).csv")))
  
  # Boxplot by clinical group
  for(clinical.info in names(compare_groups)){
    patch_freq_roiTemp <- patch_freq_roi[,c("sample_id","n_patches","ratio",clinical.info)]
    patch_freq_roiTemp <- as.data.frame(na.omit(patch_freq_roiTemp)) 
    patch_freq_roiTemp$group <- as.factor(patch_freq_roiTemp[,clinical.info])
    
    p1 <- ggplot(patch_freq_roiTemp, aes(x = group, y = n_patches, fill = group)) +
      geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
      geom_jitter(color = "darkgrey", position = position_jitter(width = 0.2), size = 0.1, alpha = 0.3) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        strip.background = element_blank()
      ) +
      labs(title = paste0("Number of ",patch_type_," by ",clinical.info), x = "Group", y = "Number of Patches")+
      scale_fill_manual(values = ggsci::pal_jco("default")(length(unique(patch_freq_roiTemp$group)))) +
      stat_compare_means(
        comparisons = compare_groups[[clinical.info]], 
        method = "wilcox.test",        # More robust than t-test
        p.adjust.method = "fdr",       # Less conservative than Bonferroni
        label = "p.adj",               # Shows adjusted p-values
        hide.ns = FALSE                 # Cleaner visualization
      )
    
    p2 <- ggplot(patch_freq_roiTemp, aes(x = group, y = ratio, fill = group)) +
      geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
      geom_jitter(color = "darkgrey", position = position_jitter(width = 0.2), size = 0.1, alpha = 0.3) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        strip.background = element_blank()
      ) +
      labs(title = paste0("Fraction of ",patch_type_," by ",clinical.info), x = "Group", y = "Fraction of Patches")+
      scale_fill_manual(values = ggsci::pal_jco("default")(length(unique(patch_freq_roiTemp$group)))) +
      stat_compare_means(
        comparisons = compare_groups[[clinical.info]], 
        method = "wilcox.test",        # More robust than t-test
        p.adjust.method = "fdr",       # Less conservative than Bonferroni
        label = "p.adj",               # Shows adjusted p-values
        hide.ns = FALSE                 # Cleaner visualization
      )
    
    pdf(file = file.path(figureDir,paste0(patch_type_," boxplot for Patch number and size comparison in ",clinical.info,".pdf")),width = 10,height = 5)
    print(p1 + p2)
    dev.off()
  }
  
  # --- 2. Patch Cellular Composition ---  
  # Microenvironment composition
  microenv_comp <- analysisDF %>%
    filter(!is.na(patch_type) & !patch_name) %>%  # Non-CD45_TC cells in patches
    group_by(sample_id, patch_type, sub_celltype) %>%
    summarise(n_cells = n()) %>%
    group_by(sample_id, patch_type) %>%
    mutate(prop = n_cells / sum(n_cells)) %>%     # Proportion per patch
    group_by(sample_id, sub_celltype) %>%
    summarise(avg_prop = mean(prop)) %>%          # Average across patches
    left_join(analysisDF %>% select(sample_id, names(compare_groups)) %>% distinct(), by = "sample_id")
  
  for(clinical.info in names(compare_groups)){
    # Aggregate by clinical group
    microenv_compTemp <- microenv_comp[,c("sample_id","sub_celltype","avg_prop",clinical.info)]
    colnames(microenv_compTemp)[4] <- "group"
    microenv_compTemp <- microenv_compTemp[!is.na(microenv_compTemp$group),]
    
    microenv_group <- microenv_compTemp %>%
      group_by(group, sub_celltype) %>%
      summarise(mean_prop = mean(avg_prop))
    
    microenv_group$group <- as.factor(microenv_group$group)
    
    # Stacked bar plot for microenvironment composition
    p <- ggplot(microenv_group, aes(x = group, y = mean_prop, fill = sub_celltype)) +
      geom_bar(stat = "identity", position = "stack") +
      theme_classic2()+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
      labs(title = "Microenvironment Composition in ",patch_name_," Patches", x = "Group", y = "Mean Proportion")+
      scale_fill_manual(values = metadata(spe)$color_vectors[["sub_celltype"]])
    
    pdf(file = file.path(figureDir,paste0(patch_type_," of ",clinical.info," microenvironment composition.pdf")),width = 5,height = 8)
    print(p)
    dev.off()
    
    #write.csv(microenv_comp,file.path(figureDir,paste0(patch_type_," microenvironment composition.csv")))
    
    ## Compare across groups of subpopulation
    groups <- unique(microenv_compTemp$group)
    
    # Loop through each clinical group
    diff_results <- data.frame()
    for (grp in groups) {
      # Create a binary column: 1 if in the group, 0 otherwise
      microenv_compTemp$in_group <- ifelse(microenv_compTemp$group == grp, 1, 0)
      
      # Get unique cell subpopulations
      cell_types <- unique(microenv_compTemp$sub_celltype)
      
      # Perform Wilcoxon test for each cell type
      for (cell_type in cell_types) {
        # Filter data for the current cell type
        cell_data <- microenv_compTemp %>% filter(sub_celltype == cell_type)
        
        if(!(grp %in% unique(cell_data$group))){
          next;
        }
        if(length(unique(cell_data$group))==1){
          next;
        }
        # Perform Wilcoxon test (non-parametric test for comparing two groups)
        test_result <- wilcox.test(avg_prop ~ in_group, data = cell_data)
        
        # Calculate fold change (mean in group / mean in others)
        mean_in_group <- mean(cell_data$avg_prop[cell_data$in_group == 1], na.rm = TRUE)
        mean_others <- mean(cell_data$avg_prop[cell_data$in_group == 0], na.rm = TRUE)
        fold_change <- mean_in_group / mean_others
        log2_fc <- log2(fold_change)  # Log2 fold change for direction and magnitude
        
        # Store results
        diff_results <- rbind(diff_results, data.frame(
          group = grp,
          major_celltype = cell_type,
          p_value = test_result$p.value,
          fold_change = fold_change,
          log2_fc = log2_fc
        ))
      }
    }
    
    # Adjust p-values for multiple testing using Benjamini-Hochberg correction
    diff_results <- diff_results %>%
      group_by(group) %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      ungroup()
    
    # Determine significance and direction
    p_sig_cutoff <- 0.05
    logfc_cutoff <- 1
    
    # Create the bubble plot
    if(T){
      diff_results <- diff_results %>%
        mutate(
          significance = case_when(
            p_value <= p_sig_cutoff & log2_fc >= logfc_cutoff ~ "Up",
            p_value <= p_sig_cutoff & log2_fc <= logfc_cutoff ~ "Down",
            TRUE ~ "NS"  # Not significant
          )
        )
      
      p <- ggplot(diff_results, aes(x = log2_fc, y = -log10(p_adj))) +
        # Plot all points in grey as background
        geom_point(color = "grey", size = 0.8) +
        # Overlay significant points in color
        geom_point(data = subset(diff_results, significance %in% c("Up", "Down")), 
                   aes(color = significance), size = 0.8) +
        # Set color scale for significant points21
        scale_color_manual(values = c("Up" = "#cd534cff", "Down" = "#7ca6dcff", "NS" = "grey")) +
        # Facet by group, one row with multiple columns
        facet_grid(. ~ group) +
        # Add threshold lines: vertical for log2_fc, horizontal for p_adj = 0.05
        geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "grey50", size = 0.5) +
        geom_hline(yintercept = -log10(p_sig_cutoff), linetype = "dashed", color = "grey50", size = 0.5) +
        # Use a clean theme similar to the new code
        theme_bw() +
        theme(
          panel.grid = element_blank(),           # Remove grid lines
          axis.text = element_text(size = 10),    # Axis text size from new code
          strip.text.x = element_text(size = 10, face = "bold"),  # Bold facet titles
          legend.position = "bottom",             # Legend at bottom, inspired by your code
          axis.title = element_text(face = "bold", size = 14),    # Bold titles from your code
          text = element_text(size = 12)          # General text size from your code
        ) +
        # Add labels and title
        labs(
          title = paste0("Differential Abundance of Cell Subpopulations by ",clinical.info),
          x = "Log2 Fold Change",
          y = "-Log10 Adjusted P-value",
          color = "Significance"
        ) +
        # Add labels for significant points, matching your original approach
        geom_text_repel(
          data = subset(diff_results, significance %in% c("Up", "Down")),
          aes(label = major_celltype),
          size = 2,
          box.padding = 0.5,
          max.overlaps = Inf
        )
      
      pdf(file = file.path(figureDir,paste0(patch_type_," of ",clinical.info," cell subpopulation fraction difference.pdf")),width = 8,height = 6)
      print(p)
      dev.off()
      
      # write.csv(diff_results,file.path(figureDir,paste0(patch_type_," cell subpopulation fraction difference.csv") ))
    }
    
  }
    # --- 3. Patch Cellular Expression ---  
  if(T){
    # Step 1: Filter cells by patch_type
    spe_filtered <- spe[,rownames(analysisDF)]
    
    # Step 2: Get expression matrix and metadata
    expr_matrix_all <- assay(spe_filtered)
    
    # Step 3: Calculate mean expression per celltype per sample
    
    for(clinical.info in names(compare_groups)){
      
      ## Get subset expression matrix and metadata
      cell_meta <- analysisDF[,c("sample_id","sub_celltype",clinical.info,"patch_type","patch_name")]
      cell_meta <- cell_meta[!is.na(cell_meta$patch_type),]
      cell_meta <- cell_meta[!is.na(cell_meta[,clinical.info]),]
      
      expr_matrix <- expr_matrix_all[,match(rownames(cell_meta),colnames(expr_matrix_all))]
      
      # Calculate mean expression for each marker in each group
      mean_expr_list <- list()
      
      for (i in 1:nrow(expr_matrix)) {
        marker_data <- data.frame(
          expression = expr_matrix[i, ],
          sample_id = cell_meta$sample_id,
          sub_celltype = cell_meta$sub_celltype,
          clinical_group = cell_meta[[clinical.info]],
          stringsAsFactors = FALSE
        )
        
        # Calculate mean expression per celltype per sample per clinical group
        marker_summary <- marker_data %>%
          group_by(sample_id, sub_celltype, clinical_group) %>%
          summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = 'drop') %>%
          mutate(marker = rownames(expr_matrix)[i])
        
        mean_expr_list[[i]] <- marker_summary
      }
      
      # Combine all markers
      all_mean_expr <- do.call(rbind, mean_expr_list)
      group_vals <- unique(all_mean_expr$clinical_group)
      n_groups <- length(group_vals)
      
      # Step 4: Perform statistical testing and calculate fold change
      results_list <- list()
      
      celltypes <- unique(all_mean_expr$sub_celltype)
      markers <- unique(all_mean_expr$marker)
      markers <- markers[!startsWith(markers,prefix = "DNA")]
      
      for (celltype in celltypes) {
        for (marker in markers) {
          # Get data for this celltype and marker
          subset_data <- all_mean_expr %>%
            filter(sub_celltype == celltype, marker == !!marker)
          
          if (nrow(subset_data) < n_groups * 3) next  # Need at least 3 samples per group
          
          # Check if we have data for all groups
          groups_present <- unique(subset_data$clinical_group)
          if (length(groups_present) < 2) next
          
          if (n_groups == 2) {
            # Two groups: use t-test (original logic)
            group1_data <- subset_data %>% filter(clinical_group == group_vals[1]) %>% pull(mean_expr)
            group2_data <- subset_data %>% filter(clinical_group == group_vals[2]) %>% pull(mean_expr)
            
            if (length(group1_data) == 0 || length(group2_data) == 0) next
            
            # Calculate fold change (group1 vs group2)
            mean_group1 <- mean(group1_data, na.rm = TRUE)
            mean_group2 <- mean(group2_data, na.rm = TRUE)
            fold_change <- log2((mean_group1 + 0.001) / (mean_group2 + 0.001))
            
            # Perform t-test
            if (length(group1_data) > 1 && length(group2_data) > 1) {
              test_result <- t.test(group1_data, group2_data)
              p_value <- test_result$p.value
            } else {
              p_value <- 1
            }
            
            results_list[[paste(celltype, marker, sep = "_")]] <- data.frame(
              sub_celltype = celltype,
              marker = marker,
              fold_change = fold_change,
              p_value = p_value,
              mean_group1 = mean_group1,
              mean_group2 = mean_group2,
              n_group1 = length(group1_data),
              n_group2 = length(group2_data),
              comparison = paste(group_vals[1], "vs", group_vals[2])
            )
            
          } else {
            # Multiple groups: use "1 vs. remaining" approach
            
            # Perform each group vs. all others
            for (target_group in group_vals) {
              
              # Split data into target group vs. all others
              target_data <- subset_data %>% filter(clinical_group == target_group) %>% pull(mean_expr)
              other_data <- subset_data %>% filter(clinical_group != target_group) %>% pull(mean_expr)
              
              if (length(target_data) == 0 || length(other_data) == 0) next
              
              # Calculate fold change (target vs others)
              mean_target <- mean(target_data, na.rm = TRUE)
              mean_others <- mean(other_data, na.rm = TRUE)
              fold_change <- log2((mean_target + 0.001) / (mean_others + 0.001))
              
              # Perform t-test: target group vs. all others
              if (length(target_data) > 1 && length(other_data) > 1) {
                test_result <- t.test(target_data, other_data)
                p_value <- test_result$p.value
              } else {
                p_value <- 1        
              }
              
              comparison_name <- paste(celltype, marker, target_group, "vs_others", sep = "_")
              results_list[[comparison_name]] <- data.frame(
                sub_celltype = celltype,
                marker = marker,
                fold_change = fold_change,
                p_value = p_value,
                mean_group1 = mean_target,
                mean_group2 = mean_others,
                n_group1 = length(target_data),
                n_group2 = length(other_data),
                comparison = paste(target_group, "vs Others")
              )
            }
          }
        }
      }
      
      # Combine results
      if (length(results_list) == 0) {
        cat("No valid comparisons found\n")
        next
      }
      
      results_df <- do.call(rbind, results_list)
      results_df$p_adj <- p.adjust(results_df$p_value, method = "BH") # Adjust p-values
      
      # Classify changes        
      results_df$change_type <- ifelse(results_df$p_adj <= 0.05 & results_df$fold_change >= logfc_cutoff, "up",
                                       ifelse(results_df$p_adj <= 0.05 & results_df$fold_change <= -logfc_cutoff, "down", 
                                              "n.s."))
      
      # Calculate transparency (-log10(p_adj))
      results_df$neg_log_p <- -log10(results_df$p_adj + 1e-10)  # Add small value to avoid infinite
      
      # Step 5: Create bubble plot
      if (n_groups == 2) {
        # Simple plot for two groups
        p <- ggplot(results_df, aes(x = sub_celltype, y = marker, 
                                    size = abs(fold_change), 
                                    color = change_type,
                                    alpha = neg_log_p)) +
          geom_point() +
          scale_color_manual(values = c("up" = "#cd534cff", "down" = "#7ca6dcff", "n.s." = "gray"),name = "Change") +
          scale_size_continuous(name = "Abs(log2FC)", range = c(1, 8)) +
          scale_alpha_continuous(name = "-log10(p.adj)", range = c(0.3, 1)) +
          theme_classic2() +
          theme(
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "lightgray"),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(hjust = 0.5)
          ) +
          labs(title = paste("Differential Expression:",patch_name_ ),x = "Cell Type",y = "Marker") +
          guides(
            color = guide_legend(override.aes = list(size = 5, alpha = 1)),
            size = guide_legend(override.aes = list(alpha = 1)),
            alpha = guide_legend(override.aes = list(size = 5))
          )
        
        pdf(file = file.path(figureDir,paste0(patch_type_," of ",clinical.info," cellular expression change.pdf")),width = 10,height = 7.5)
        print(p)
        dev.off()
      } 
      else {
        # Multiple comparison plot - create faceted plot by comparison
        p <- ggplot(results_df, aes(x = sub_celltype, y = marker, 
                                    size = abs(fold_change), 
                                    color = change_type,
                                    alpha = neg_log_p)) +
          geom_point() +
          facet_wrap(~comparison, scales = "free") +
          scale_color_manual(values = c("up" = "#cd534cff", "down" = "#7ca6dcff", "n.s." = "gray"),name = "Change") +
          scale_size_continuous(name = "Abs(log2FC)", range = c(1, 6)) +
          scale_alpha_continuous(name = "-log10(p.adj)", range = c(0.3, 1)) +
          theme_classic2() +
          theme(
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "lightgray"),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
            axis.text.y = element_text(size = 7),
            plot.title = element_text(hjust = 0.5),
            strip.text = element_text(size = 8)
          ) +
          labs(title = paste("One vs Others Differential Expression:",patch_name_ ),x = "Cell Type",y = "Marker") +
          guides(
            color = guide_legend(override.aes = list(size = 5, alpha = 1)),
            size = guide_legend(override.aes = list(alpha = 1)),
            alpha = guide_legend(override.aes = list(size = 5)))
        
        pdf(file = file.path(figureDir,paste0(patch_type_," of ",clinical.info," one_vs_others cellular expression change.pdf")),width = 15,height = 7.5)
        print(p)
        dev.off()
      }
    }
  }

  analysisDF <- analysisDF_Back
}

# --- 3. Overlap and Network Analysis of Patches ---
# Assuming analysisDF is your data frame
# Identify patch columns
patch_cols <- grep("_patch$", names(analysisDF), value = TRUE)

# Define short labels for patch types
patch_labels <- c(
  "CD45_TC_patch" = "CD45",
  "YAP_TC_DTC_patch" = "YAP_DTC",
  "PD1_T_patch" = "PD1_T",
  "Treg_patch" = "Treg",
  "PDL1_Macro_patch" = "PDL1_M",
  "Endo_patch" = "Endo",
  "MC_2_patch" = "MC_2",
  "YAP_Fibro_patch" = "YAP_Fibro"
)
all_vertices <- unique(patch_labels)

# Get unique clinical groups
groups <- unique(analysisDF$group)

# Function to create a graph for a given group
create_group_graph <- function(df_grp, patch_cols, patch_labels) {
  # Get all pairwise combinations of patch columns
  pairs <- combn(patch_cols, 2, simplify = FALSE)
  
  # Calculate intersection sizes for each pair
  edge_list <- lapply(pairs, function(pair) {
    patch1 <- pair[1]
    patch2 <- pair[2]
    # Intersection size: count where both patches are non-NA
    weight <- sum(!is.na(df_grp[[patch1]]) & !is.na(df_grp[[patch2]]))
    c(patch_labels[patch1], patch_labels[patch2], weight)
  })
  
  # Convert to data frame
  edge_df <- do.call(rbind, edge_list)
  edge_df <- as.data.frame(edge_df)
  names(edge_df) <- c("from", "to", "weight")
  edge_df$weight <- as.numeric(edge_df$weight)
  
  # Filter out edges with zero weight
  edge_df <- edge_df %>% filter(weight > 0)
  
  # Create graph, ensuring all vertices are included (even isolates)
  g <- graph_from_data_frame(edge_df, directed = FALSE, 
                             vertices = data.frame(name = all_vertices))
  return(g)
}

# Create a list of graphs, one for each group
group_graphs <- lapply(groups, function(grp) {
  df_grp <- analysisDF %>% filter(group == grp)
  create_group_graph(df_grp, patch_cols, patch_labels)
})
names(group_graphs) <- groups
saveRDS(groups,file.path(figureDir,paste0("Network Analysis of Patches plot data.rds")))

# Visualize each graph using ggraph
plot_list <- lapply(groups, function(grp) {
  g <- group_graphs[[grp]]
  
  # Create the ggraph plot
  ggraph(g, layout = "circle") +
    # Edges with width proportional to intersection size
    geom_edge_link(aes(edge_width = weight), edge_colour = "grey50", alpha = 0.7) +
    # Nodes as points
    geom_node_point(size = 12, color = "lightblue", alpha = 0.9) +
    # Node labels
    geom_node_text(aes(label = name), repel = TRUE, size = 5, fontface = "bold") +
    # Scale edge width for better visualization
    scale_edge_width(range = c(0.5, 5), name = "Intersection Size") +
    # Clean theme
    theme_void() +
    # Titles and legend
    labs(title = paste("Group:", grp)) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )})

p <- wrap_plots(plot_list) + plot_layout(ncol = 2)

pdf(file = file.path(figureDir,"Patch Interaction.pdf"),width = 20,height = 16)
print(p)
dev.off()
# For spatial analysis

library(SpatialExperiment)
library(BiocParallel)

library(imcRtools)
library(lisaClust)

library(ggplot2)
library(ggridges)
library(ggpubr)
library(ggraph)
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
CN_color <- setNames(metadata(spe_filtered)$color_vectors$color_10,unique(meta_df[[cn_column]]))
RFS_color <- setNames(c("#EFC000FF", "#0073C2FF"),unique(meta_df$RFS_status))
Tissue_color <- setNames(metadata(spe_filtered)$color_vectors$color_20[1:3],unique(meta_df$Tissue))

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
pdf(file.path(figureDir, paste0("CN_correlation_matrix_", cn_column, ".pdf")), 
    width = 8, height = 6)
corrplot(cn_correlation_matrix, 
              method = "color", type = "upper",order = "hclust",
              tl.col = "black", tl.srt = 90, col = colorRampPalette(c("darkblue", "white", "darkred"))(100),
              title = paste("CN-CN Correlations (Patient Level) -", cn_column),
              mar = c(0,0,2,0),addCoef.col = "black",number.cex = 0.7)
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
  clinical_analysis_data <- read.csv(file.path(figureDir, paste0("CN_data_for_clinical_analysis_", cn_column, ".csv")))
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

# PREPARE DATA FOR TASK 5 (BOXPLOTS)
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

# Save plotting data for Task 5
write.csv(cn_long_data, 
          file.path(figureDir, paste0("CN_data_for_boxplots_", cn_column, ".csv")), 
          row.names = FALSE)

# ============================================================================
# Task 5: Boxplot Visualization with Significance Testing
# ============================================================================
print("Starting Task 5: Boxplot Visualization with Significance Testing")

# Load plotting data
if (!exists("cn_long_data")) {
  cn_long_data <- read.csv(file.path(figureDir, paste0("CN_data_for_boxplots_", cn_column, ".csv")))
  print("Loaded plotting data from file")
}

print(paste("Plotting data:", nrow(cn_long_data), "observations"))

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
  
  # Create boxplot
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
      x = "RFS Status", 
      y = "CN Fraction"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))
  
  # Add significance annotation
  y_max <- max(cn_data$CN_fraction, na.rm = TRUE)
  y_pos <- y_max * 1.1
  
  p_overall <- p_overall + 
    stat_compare_means(method = "wilcox.test", label = "p", 
                       hide.ns = FALSE, size = 3)
  
  # Save individual boxplot
  ggsave(file.path(figureDir, paste0("CN", cn_val, "_RFS_overall_boxplot_", cn_column, ".pdf")), 
         p_overall, width = 6, height = 5)
}

# TASK 5B: INDIVIDUAL CN BOXPLOTS - TREATMENT-STRATIFIED COMPARISON
print("Task 5B: Creating individual CN boxplots - Treatment-stratified comparison")

# Get available treatments
treatments <- unique(cn_long_data$Treatment)
print(paste("Treatments:", paste(treatments, collapse = ", ")))

# Create treatment-stratified boxplots for each CN
for (cn_val in unique_cns) {
  print(paste("Creating treatment-stratified boxplots for CN", cn_val))
  
    # Filter data for this CN and treatment
    cn_treatment_data <- cn_long_data[cn_long_data$CN == cn_val,]
    
    if (nrow(cn_treatment_data) < 4) {
      print(paste("Insufficient data for CN", cn_val, "in", treatment))
      next
    }
    
    # Check if we have both RFS groups
    if (length(unique(cn_treatment_data$RFS_status)) < 2) {
      print(paste("Only one RFS group for CN", cn_val, "in", treatment))
      next
    }
    
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
        x = "RFS Status", 
        y = "CN Fraction"
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))
    
    # Add significance annotation
    y_max <- max(cn_treatment_data$CN_fraction, na.rm = TRUE)
    y_pos <- y_max * 1.1
    
    p_treatment <- p_treatment +
      stat_compare_means(method = "wilcox.test", label = "p", 
                         hide.ns = FALSE, size = 3) + 
      facet_grid(~Treatment)
  
  # Combine treatment plots if we have both treatments
    ggsave(file.path(figureDir, paste0("CN", cn_val, "_RFS_treatment_stratified_", cn_column, ".pdf")), 
           p_treatment, width = 12, height = 5)
  
}

# TASK 5C: ALL CNs TOGETHER by Tissues
print("Task 5C: Creating summary boxplots - All CNs together - by tissue")

cn_long_data_roi <- cn_fractions_wide
cn_long_data_roi$RFS_label <- ifelse(cn_long_data_roi$RFS_status == 0, "No Early Relapse", "Early Relapse")
cn_long_data_roi <- cn_long_data_roi %>%
            select(sample_id, RFS_status, Tissue, RFS_label, Treatment, all_of(cn_columns)) %>%
            pivot_longer(cols = all_of(cn_columns), names_to = "CN", values_to = "CN_fraction") %>%
            mutate(CN = gsub("CN_", "", CN))
cn_long_data_roi$RFS_status <- as.factor(cn_long_data_roi$RFS_status)

# Overall RFS comparison - all CNs in one plot
p_summary_overall <- ggplot(cn_long_data_roi, aes(x = Tissue, y = CN_fraction, fill = RFS_status)) +
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
  stat_compare_means(method = "wilcox.test", label = "p", 
                     hide.ns = FALSE, size = 3)

ggsave(file.path(figureDir, paste0("CN_all_RFS_overall_summary_", cn_column, ".pdf")), 
       p_summary_overall, width = 16, height = 12)

# Treatment-stratified comparison - all CNs
p_summary_treatment <- ggplot(cn_long_data_roi, aes(x = Tissue, y = CN_fraction, fill = RFS_status)) +
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
  stat_compare_means(method = "wilcox.test", label = "p", 
                     hide.ns = FALSE, size = 3)

ggsave(file.path(figureDir, paste0("CN_all_RFS_treatment_summary_", cn_column, ".pdf")), 
       p_summary_treatment, width = 20, height = 8)

# ============================================================================
# Task 6: CN Functionality Analysis - Epithelial vs T cell CNs
# ============================================================================
# Based on proportion table analysis identifying epithelial vs T cell enriched CNs
print("Starting Task 6: Functionality Analysis - Epithelial vs T cell CNs")

# CN GROUP SELECTIONS BASED ON PROPORTION TABLE ANALYSIS
# Epithelial-enriched CNs (based on EC_* cell enrichment)
epithelial_cns <- c("3", "5", "7", "8")

# T cell-enriched CNs (based on CD4T, CD8T, Treg enrichment)
tcell_cns <- c("1", "5", "6", "8", "9")

# Note: CN5 and CN8 appear in both groups (mixed epithelial-immune neighborhoods)
print("=== CN GROUP SELECTIONS ===")
print(paste("Epithelial-enriched CNs:", paste(epithelial_cns, collapse = ", ")))
print(paste("T cell-enriched CNs:", paste(tcell_cns, collapse = ", ")))
print(paste("Mixed CNs (both groups):", paste(intersect(epithelial_cns, tcell_cns), collapse = ", ")))

# FUNCTIONALITY MARKER SELECTIONS
# Epithelial functionality markers (metabolism, hypoxia, proliferation, immune evasion)
epithelial_markers <- c(
  "Ki67",      # Proliferation
  "GLUT1",     # Glucose metabolism  
  "CA_IX",     # Hypoxia response
  "FASN",      # Fatty acid synthesis
  "HK2",       # Glycolysis
  "PRPS1",     # Nucleotide synthesis
  "VEGF",      # Angiogenesis
  "CD274"      # PD-L1 (immune evasion)
)

# T cell functionality markers (exhaustion, activation, memory)
tcell_markers <- c(
  "CD279",     # PD-1 (exhaustion)
  "TIGIT",     # Exhaustion
  "CD366",     # TIM-3 (exhaustion)
  "CD127",     # IL-7R (memory)
  "CD27",      # Costimulation
  "Ki67",      # T cell proliferation
  "CD274"      # PD-L1 (for T cell-tumor interactions)
)

print("=== FUNCTIONALITY MARKER SELECTIONS ===")
print(paste("Epithelial markers:", paste(epithelial_markers, collapse = ", ")))
print(paste("T cell markers:", paste(tcell_markers, collapse = ", ")))

# CELL TYPE SELECTIONS FOR ANALYSIS
# Epithelial cell types (from proportion table)
epithelial_celltypes <- c(
  "EC_CAIX", "EC_EpCAM", "EC_GLUT1", "EC_Ki67", 
  "EC_Ki67_CAIX", "EC_Vimentin"
)

# T cell types
tcell_celltypes <- c(
  "CD4T", "CD8T", "Treg", "NK"
)

print("=== CELL TYPE SELECTIONS ===")
print(paste("Epithelial cell types:", paste(epithelial_celltypes, collapse = ", ")))
print(paste("T cell types:", paste(tcell_celltypes, collapse = ", ")))

# DATA PREPARATION
print("Preparing data for functionality analysis...")

# Use the clinical analysis data from previous tasks
if (!exists("clinical_analysis_data")) {
  clinical_analysis_data <- read.csv(file.path(figureDir, paste0("CN_data_for_clinical_analysis_", cn_column, ".csv")))
}

# Also need the original metadata for cell-level analysis
meta_df_func <- meta_df

# Filter for cells in CNs of interest
all_cns_of_interest <- unique(c(epithelial_cns, tcell_cns))
meta_df_func <- meta_df_func[meta_df_func[[cn_column]] %in% all_cns_of_interest, ]

print(paste("Cells in CNs of interest:", nrow(meta_df_func)))

# Check available markers in data
available_epithelial_markers <- intersect(epithelial_markers, rownames(spe_filtered))
available_tcell_markers <- intersect(tcell_markers, rownames(spe_filtered))

print("=== AVAILABLE MARKERS CHECK ===")
print(paste("Available epithelial markers:", paste(available_epithelial_markers, collapse = ", ")))
print(paste("Available T cell markers:", paste(available_tcell_markers, collapse = ", ")))

# EPITHELIAL CN FUNCTIONALITY ANALYSIS
print("=== EPITHELIAL CN FUNCTIONALITY ANALYSIS ===")

# Filter for epithelial cells in epithelial CNs
epithelial_data <- meta_df_func[
  meta_df_func[[cn_column]] %in% epithelial_cns & 
    meta_df_func$sub_celltype %in% epithelial_celltypes, 
]

print(paste("Epithelial cells in epithelial CNs:", nrow(epithelial_data)))

if (nrow(epithelial_data) > 0) {
  # Get expression data for epithelial cells
  epithelial_expr <- assay(spe_filtered)[available_epithelial_markers, rownames(epithelial_data)]
  epithelial_expr_df <- as.data.frame(t(epithelial_expr))
  
  # Combine with metadata
  epithelial_func_df <- cbind(epithelial_data, epithelial_expr_df)
  
  # Calculate mean expression per marker per CN per cell type per patient
  epithelial_summary <- epithelial_func_df %>%
    group_by(patient_id, !!sym(cn_column), sub_celltype, RFS_status, Treatment) %>%
    summarise(across(all_of(available_epithelial_markers), mean, na.rm = TRUE), .groups = 'drop')
  
  print(paste("Epithelial summary records:", nrow(epithelial_summary)))
  
  # Save epithelial functionality data
  write.csv(epithelial_summary, 
            file.path(figureDir, paste0("epithelial_functionality_summary_", cn_column, ".csv")), 
            row.names = FALSE)
}

# T CELL CN FUNCTIONALITY ANALYSIS  
print("=== T CELL CN FUNCTIONALITY ANALYSIS ===")

# Filter for T cells in T cell CNs
tcell_data <- meta_df_func[
  meta_df_func[[cn_column]] %in% tcell_cns & 
    meta_df_func$sub_celltype %in% tcell_celltypes, 
]

print(paste("T cells in T cell CNs:", nrow(tcell_data)))

if (nrow(tcell_data) > 0) {
  # Get expression data for T cells
  tcell_expr <- assay(spe_filtered)[available_tcell_markers, rownames(tcell_data)]
  tcell_expr_df <- as.data.frame(t(tcell_expr))
  
  # Combine with metadata
  tcell_func_df <- cbind(tcell_data, tcell_expr_df)
  
  # Calculate mean expression per marker per CN per cell type per patient
  tcell_summary <- tcell_func_df %>%
    group_by(patient_id, !!sym(cn_column), sub_celltype, RFS_status, Treatment) %>%
    summarise(across(all_of(available_tcell_markers), mean, na.rm = TRUE), .groups = 'drop')
  
  print(paste("T cell summary records:", nrow(tcell_summary)))
  
  # Save T cell functionality data
  write.csv(tcell_summary, 
            file.path(figureDir, paste0("tcell_functionality_summary_", cn_column, ".csv")), 
            row.names = FALSE)
}

print("Task 6 completed successfully!")

# ============================================================================
# Task 7: Violin Plot Visualization for Functionality Analysis
# ============================================================================
# This continues from Task 6 - make sure you have run Task 6 first
print("Starting Task 7: Violin Plot Visualization for Functionality Analysis")

# LOAD DATA FROM TASK 6
print("Loading functionality data from Task 6...")

# Load functionality summaries
epithelial_summary <- read.csv(file.path(figureDir, paste0("epithelial_functionality_summary_", cn_column, ".csv")))
tcell_summary <- read.csv(file.path(figureDir, paste0("tcell_functionality_summary_", cn_column, ".csv")))

print(paste("Epithelial data records:", nrow(epithelial_summary)))
print(paste("T cell data records:", nrow(tcell_summary)))

# Get CN groups from Task 6
print("=== MARKERS AVAILABLE FOR VISUALIZATION ===")
print(paste("Epithelial markers:", paste(epithelial_markers, collapse = ", ")))
print(paste("T cell markers:", paste(tcell_markers, collapse = ", ")))


# TASK 7A: EPITHELIAL CN VIOLIN PLOTS
print("Task 7A: Creating epithelial CN violin plots")
if (nrow(epithelial_summary) > 0 && length(epithelial_markers) > 0) {
  
  # Create violin plots for epithelial markers across epithelial CNs
  # Reshape data
  plotdf <- epithelial_summary
  plotdf$RFS_status <- as.factor(plotdf$RFS_status)
  plotdf <- pivot_longer(epithelial_summary, 
                         cols = (ncol(epithelial_summary) - length(epithelial_markers) + 1):ncol(epithelial_summary),
                         values_to = "expression",names_to = "marker")
  
  plotdf$CN_factor <- factor(plotdf[[cn_column]], levels = epithelial_cns)
  plotdf$RFS_label <- factor(plotdf$RFS_status, 
                             levels = c(0, 1), labels = c("No Early Relapse", "Early Relapse"))
  

  print(paste("  Creating epithelial violin plot for epithelial CNs"))
  p_cn <- ggplot(plotdf, aes(x = CN_factor, y = expression, fill = CN_factor)) +
    geom_violin(alpha = 0.7, scale = "width", trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 0.75, color = "grey75") +
    scale_fill_manual(values = CN_color) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "top",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste(marker, "Expression in", celltype, "across CNs"),
      x = "Cellular Neighborhood", 
      y = paste(marker, "Expression")
    ) +
    stat_compare_means(comparisons = combn(as.character(unique(plotdf$CN_factor)), 2, simplify = FALSE),
                       method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    facet_grid(sub_celltype~marker)
  
    ggsave(file.path(figureDir, paste0("epithelial_", marker, "_CN_comparison_", cn_column, ".pdf")), 
         p_cn, width = 20, height = 16)
    
    # RFS comparison plot
    rfs_groups <- unique(plotdf$RFS_status)
    plotdf$RFS_status <- as.factor(plotdf$RFS_status)
    # CN-stratified RFS comparison (if we have multiple CNs)

      p_cn_rfs <- ggplot(plotdf, aes(x = CN_factor, y = expression, fill = RFS_status)) +
        geom_violin(alpha = 0.7, scale = "width", trim = FALSE) +
        facet_wrap(~ paste("CN", CN_factor), scales = "free_y") +
        scale_fill_manual(values = RFS_color) +
        theme_bw() +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 10),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 11, face = "bold"),
          legend.position = "bottom",
          panel.grid.minor = element_blank()
        ) +
        labs(
          title = paste(marker, "Expression - CN-Stratified RFS Comparison"),
          x = "RFS Status", 
          y = paste(marker, "Expression"),
          fill = "RFS Status"
        ) +
        stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3, hide.ns = TRUE) +
        facet_grid(marker~sub_celltype, scales = "free_y")
      
      ggsave(file.path(figureDir, paste0("epithelial_", marker, "_CN_stratified_RFS_", cn_column, ".pdf")), 
             p_cn_rfs, width = 24, height = 18)
  
}

# TASK 7B: T CELL CN VIOLIN PLOTS
print("Task 7B: Creating T cell CN violin plots")

if (nrow(tcell_summary) > 0 && length(tcell_markers) > 0) {
  
  # Reshape data
  plotdf <- tcell_summary
  plotdf$RFS_status <- as.factor(plotdf$RFS_status)
  plotdf <- pivot_longer(tcell_summary, 
                         cols = (ncol(tcell_summary) - length(tcell_markers) + 1):ncol(tcell_summary),
                         values_to = "expression",names_to = "marker")
  
  plotdf$CN_factor <- factor(plotdf[[cn_column]], levels = tcell_cns)
  plotdf$RFS_label <- factor(plotdf$RFS_status, 
                                    levels = c(0, 1), labels = c("No Early Relapse", "Early Relapse"))
  
  # Create violin plots for T cell markers across T cell CNs
  p_cn <- ggplot(plotdf, aes(x = CN_factor, y = expression, fill = CN_factor)) +
    geom_violin(alpha = 0.7, scale = "width", trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 0.75, color = "grey75") +
    scale_fill_manual(values = CN_color) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "top",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste(marker, "Expression in", celltype, "across CNs"),
      x = "Cellular Neighborhood", 
      y = paste(marker, "Expression")
    ) +
    stat_compare_means(comparisons = combn(as.character(unique(plotdf$CN_factor)), 2, simplify = FALSE),
                       method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    facet_grid(sub_celltype~marker)
  
    ggsave(file.path(figureDir, paste0("tcell_", marker, "_CN_comparison_", cn_column, ".pdf")), 
         p_cn, width = 20, height = 16)
      
    # RFS comparison for T cell types
    p_rfs <- ggplot(plotdf, aes(x = RFS_label, y = expression, fill = RFS_status)) +
      geom_violin(alpha = 0.7, scale = "width", trim = FALSE) +
      geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
      geom_jitter(width = 0.1, alpha = 0.5, size = 0.75, color = "grey75") +
      scale_fill_manual(values = RFS_color) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        legend.position = "top",
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = paste(marker, "in", celltype, "- RFS Comparison"),
        x = "RFS Status", 
        y = paste(marker, "Expression")
      ) +
      stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) + 
      facet_grid(sub_celltype~marker)
      
      ggsave(file.path(figureDir, paste0("tcell_", marker, "_", celltype, "_RFS_comparison_", cn_column, ".pdf")), 
             p_rfs, width = 16, height = 12)
  }

print("Task 7 completed successfully!")

# ============================================================================
# Task 8: CN Interaction Analysis using Canonical Correlation Analysis (CCA)
# ============================================================================
print("Starting Task 8: CN Interaction Analysis using CCA Networks")

# SETUP AND PARAMETERS FROM PREVIOUS TASKS
print("Setting up CCA analysis parameters...")

# CN groups
all_cns_analysis <- unique(meta_df[[cn_column]])

# Combined functional markers for CCA (remove duplicates)
# Functional markers from Tasks 6-7
epithelial_markers <- c("Ki67", "GLUT1", "CA_IX", "FASN", "HK2", "PRPS1", "VEGF", "CD274")
tcell_markers <- c("CD279", "TIGIT", "CD366", "CD127", "CD27", "Ki67", "CD274")
cca_functional_markers <- unique(c(epithelial_markers,tcell_markers))
# cca_functional_markers <- rownames(spe_filtered)[rowData(spe_filtered)$use_channel]

# Cell types for CCA
epithelial_celltypes <- c("EC_CAIX", "EC_EpCAM", "EC_GLUT1", "EC_Ki67", "EC_Ki67_CAIX", "EC_Vimentin")
tcell_celltypes <- c("CD4T", "CD8T", "Treg", "NK")
cca_celltypes <- c(epithelial_celltypes, tcell_celltypes)
# cca_celltypes <- unique(meta_df$sub_celltype)

print("=== CCA ANALYSIS SETUP ===")
print(paste("CNs for CCA:", paste(all_cns_analysis, collapse = ", ")))
print(paste("Functional markers:", paste(cca_functional_markers, collapse = ", ")))
print(paste("Cell types:", paste(cca_celltypes, collapse = ", ")))

# DATA PREPARATION FOR CCA
print("Preparing data for CCA analysis...")

# Filter metadata for CCA analysis
cca_meta <- meta_df[
  meta_df[[cn_column]] %in% all_cns_analysis & 
    meta_df$sub_celltype %in% cca_celltypes, 
]

print(paste("Cells for CCA analysis:", nrow(cca_meta)))

# Check available markers in the data
available_cca_markers <- intersect(cca_functional_markers, rownames(spe_filtered))
print(paste("Available CCA markers:", paste(available_cca_markers, collapse = ", ")))

if (length(available_cca_markers) == 0) {
  stop("No CCA markers found in expression data. Please check marker names.")
}

# Get expression data for CCA cells and markers
# cca_expr <- assay(spe_filtered)[available_cca_markers, rownames(cca_meta)]
# cca_expr_df <- as.data.frame(t(cca_expr))

# Combine with metadata
# cca_data <- cbind(cca_meta, cca_expr_df)
cca_data <- cca_meta
print("CCA data preparation completed")

# CALCULATE CN-SPECIFIC PROFILES FOR CCA
print("Calculating CN-specific profiles for CCA...")

# For each ROI and CN, calculate:
# 1. Cell type frequencies
# 2. Mean functional marker expression
cca_patient_profiles <- cca_data %>%
  group_by(patient_id, RFS_status, Treatment, !!sym(cn_column), sub_celltype) %>%
  summarise(cell_count = n()) %>%
  # summarise(
  #   cell_count = n(),
  #   across(all_of(available_cca_markers), mean, na.rm = TRUE, .names = "mean_{.col}"),
  #   .groups = 'drop'
  # ) %>%
  # Calculate cell type frequencies within each CN per patient
  group_by(patient_id, !!sym(cn_column)) %>%
  mutate(
    total_cn_cells = sum(cell_count),
    celltype_frequency = cell_count / total_cn_cells
  ) %>%
  ungroup()

print(paste("CN-specific profiles calculated:", nrow(cca_patient_profiles)))

# Create patient-level CN feature matrices
# For each patient and CN, combine cell type frequencies and marker expressions
# cca_features <- cca_patient_profiles %>%
#   group_by(patient_id, RFS_status, Treatment, !!sym(cn_column)) %>%
#   summarise(
#     # Aggregate cell type frequencies
#     total_cells = sum(cell_count),
#     epithelial_freq = sum(celltype_frequency[sub_celltype %in% epithelial_celltypes], na.rm = TRUE),
#     tcell_freq = sum(celltype_frequency[sub_celltype %in% tcell_celltypes], na.rm = TRUE),
#     # Aggregate functional marker expressions (weighted by cell count)
#     across(starts_with("mean_"), ~ weighted.mean(.x, cell_count, na.rm = TRUE)),
#     .groups = 'drop'
#   )

cca_features <- cca_patient_profiles %>%
  # Step 1: Remove the specified columns
  select(-cell_count, -total_cn_cells) %>%
  
  # Step 2: Pivot the remaining data
  pivot_wider(
    names_from = sub_celltype,
    values_from = celltype_frequency,
    values_fill = 0
  )

print(paste("Patient-level CN features:", nrow(cca_features)))
feature_cols <- colnames(cca_features)[5:ncol(cca_features)]

# CANONICAL CORRELATION ANALYSIS BETWEEN CN PAIRS
print("Performing Canonical Correlation Analysis between CN pairs...")

# Initialize results storage
cca_results <- list()

# OVERALL CCA ANALYSIS (ALL PATIENTS)
print("CCA Analysis - Overall (all patients)")

cca_results[["overall"]] <- list()

for (i in 1:(length(all_cns_analysis)-1)) {
  for (j in (i+1):length(all_cns_analysis)) {
    cn1 <- all_cns_analysis[i]
    cn2 <- all_cns_analysis[j]
    pair_name <- paste0("CN", cn1, "_CN", cn2)
    
    print(paste("  Analyzing:", pair_name))
    
    result <- perform_cca_analysis(data = cca_features, cn1, cn2, cn_column, feature_cols)
    if (!is.null(result)) {
      cca_results[["overall"]][[pair_name]] <- result
      print(paste("    Correlation:", round(result$observed_correlation, 3), 
                  "| Samples:", result$n_patients))
    }
  }
}

# RFS-STRATIFIED CCA ANALYSIS
rfs_groups <- c(0, 1)
rfs_labels <- c("no_relapse", "early_relapse")

for (g in 1:length(rfs_groups)) {
  rfs_group <- rfs_groups[g]
  rfs_label <- rfs_labels[g]
  
  print(paste("CCA Analysis -", rfs_label))
  
  cca_results[[rfs_label]] <- list()
  
  for (i in 1:(length(all_cns_analysis)-1)) {
    for (j in (i+1):length(all_cns_analysis)) {
      cn1 <- all_cns_analysis[i]
      cn2 <- all_cns_analysis[j]
      pair_name <- paste0("CN", cn1, "_CN", cn2)
      
      result <- perform_cca_analysis(cca_features, cn1, cn2, cn_column, feature_cols, rfs_group)
      if (!is.null(result)) {
        cca_results[[rfs_label]][[pair_name]] <- result
        print(paste("    ", pair_name, "Correlation:", round(result$observed_correlation, 3)))
      }
    }
  }
}
print("CCA analysis completed")

# CREATE CORRELATION MATRICES
print("Creating correlation matrices...")

# Create correlation matrices for each group
cor_matrices <- list()
for (group_name in names(cca_results)) {
  if (length(cca_results[[group_name]]) > 0) {
    cor_matrices[[group_name]] <- create_cca_correlation_matrix(cca_results[[group_name]], cn_list = all_cns_analysis, p.adjust = FALSE)
    print(paste("Correlation matrix created for", group_name))
  }
}

# VISUALIZATION: NETWORK GRAPHS
print("Creating network visualizations...")

# Create network plots for each group
correlation_threshold <- 0.4  # Adjust as needed
p_value_threshold <- 0.1
p_list <- list()

for (group_name in names(cor_matrices)) {
  if (!is.null(cor_matrices[[group_name]])) {
    
    if (group_name == "overall") {
      title_text <- "CN Interaction Network - Overall"
    } else if (group_name == "no_relapse") {
      title_text <- "CN Interaction Network - No Early Relapse"
    } else {
      title_text <- "CN Interaction Network - Early Relapse"
    }
    
    p_network <- create_cca_network(
        matrices_data = cor_matrices[[group_name]], 
        correlation_threshold = correlation_threshold,
        p_value_threshold = p_value_threshold,
        node_colors = CN_color,
        title = paste("CN Communication Network -", group_name)
      )
    
    p_list[[group_name]] <- p_network
    
    if (!is.null(p_network)) {
      ggsave(file.path(figureDir, paste0("CCA_network_", group_name, "_", cn_column, ".pdf")), 
             p_network, width = 10, height = 8)
    }
  }
}

# NETWORK COMPARISON BETWEEN RFS GROUPS
print("Creating network comparison between RFS groups...")
    if (!is.null(p_list)) {
      # Combine plots
      p_combined <- p_list[["early_relapse"]] + p_list[["no_relapse"]]
      
      ggsave(file.path(figureDir, paste0("CCA_network_comparison_", cn_column, ".pdf")), 
             p_combined, width = 20, height = 8)
    
}

# ============================================================================
# SUMMARY STATISTICS AND RESULTS
# ============================================================================

print("Creating summary statistics...")

# Compile all CCA results into summary table
cca_summary <- data.frame()

for (group_name in names(cca_results)) {
  for (pair_name in names(cca_results[[group_name]])) {
    result <- cca_results[[group_name]][[pair_name]]
    
    cca_summary <- rbind(cca_summary, data.frame(
      group = group_name,
      cn_pair = pair_name,
      canonical_correlation = result$observed_correlation,
      p_value = result$p_value,
      n_patients = result$n_patients,
      n_permutations = result$n_permutations
      ),
      stringsAsFactors = FALSE
    )
  }
}

# Save CCA summary
write.csv(cca_summary, 
          file.path(figureDir, paste0("CCA_summary_statistics_", cn_column, ".csv")), 
          row.names = FALSE)
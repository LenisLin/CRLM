# Spatial Cellular Interaction Analysis for CRLM

library(SpatialExperiment)
library(imcRtools)
library(igraph)

library(ggplot2)
library(ggridges)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggpubr)
library(ggraph)
library(gridExtra)
library(viridis)
library(patchwork)

library(dplyr)
library(tidyr)
library(tidyverse)

library(BiocParallel)
library(lisaClust)
library(scales)
library(circlize)

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
figureDir <- file.path(workDir, "figures","3_SpatialAnalysis","Interaction_Analysis")

if(!dir.exists(figureDir)){
  dir.create(figureDir,recursive = T)
}

# Load your SpatialExperiment object
date_time <- "0728"
spe <- readRDS(file.path(saveDir,paste0("spatial_spe_",date_time,".rds")))
img_id_ <- "sample_id"

RFS_color <<- setNames(c("#0073C2FF","#EFC000FF"),sort(unique(spe$RFS_status)))

# =============================================================================
# LOAD DATA
# =============================================================================
out <- readRDS(file.path(saveDir,"Interaction_analysis_out_of_knn_20.rds"))
sample_metadata <- colData(spe) %>% 
  as.data.frame() %>%
  select(sample_id, patient_id, Tissue, Treatment, RFS_status) %>%
  distinct() %>%
  rename(group_by = sample_id)

# First, let's explore the structure of our permutation test results
cat("Basic data overview:\n")
cat("Total interactions tested:", nrow(out), "\n")
cat("Number of samples:", length(unique(out$group_by)), "\n")
cat("Number of cell types:", length(unique(out$from_label)), "\n")
cat("Significant interactions:", sum(out$sig, na.rm = TRUE), "\n")

# Check the distribution of significance values
table(out$sigval, useNA = "ifany")

# ============================================================================
# MERGE WITH CLINICAL METADATA
# ============================================================================
# Merge permutation results with clinical data
out_clinical <- as.data.frame(out) %>%
  left_join(sample_metadata, by = "group_by")

# Check merge success
cat("Samples with clinical data:", sum(!is.na(out_clinical$Tissue)), "/", nrow(out_clinical), "\n")

# ============================================================================
# OVERALL INTERACTION HEATMAP BY TISSUE
# ============================================================================
# Create summary statistics for heatmap visualization
# Focus on significant interactions first
interaction_summary <- out_clinical %>%
  filter(!is.na(Tissue)) %>%
  group_by(Tissue, from_label, to_label) %>%
  summarise(
    mean_ct = mean(ct, na.rm = TRUE),
    median_ct = median(ct, na.rm = TRUE),
    prop_significant = mean(sig, na.rm = TRUE),
    prop_attractive = mean(interaction, na.rm = TRUE),
    mean_sigval = mean(sigval, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )


# Create matrices for each Tissue
tc_matrix <- create_interaction_matrix(interaction_summary, "TC", "mean_sigval")
im_matrix <- create_interaction_matrix(interaction_summary, "IM", "mean_sigval")
pt_matrix <- create_interaction_matrix(interaction_summary, "PT", "mean_sigval")

# Generate the combined appealing heatmap
combined_heatmap <- create_combined_interaction_heatmap(spe, tc_matrix, im_matrix, pt_matrix)

# Draw with enhanced layout
pdf(file.path(figureDir, paste0("ComplexHeatmap for multiple tissues.pdf")), width = 12, height = 5)
draw(combined_heatmap, 
     main_heatmap = "TC",  # Which heatmap controls the row clustering
     column_title = "Spatial Cellular Interaction Atlas Across CRLM Tissue Regions",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE,
     ht_gap = unit(0.5, "cm"))  # Gap between heatmaps
dev.off()

# ============================================================================
# COMPARATIVE ANALYSIS BY CLINICAL GROUPS
# ============================================================================

# Define clinically relevant cell type categories
epithelial_cells <- c("EC_CAIX", "EC_EpCAM", "EC_GLUT1", "EC_Ki67", "EC_Vimentin")

stromal_cells <- c("CAF", "SC_Collagen", "SC_Collagen_Vimentin", "SC_Vimentin", 
                   "SC_aSMA_Collagen", "SC_aSMA_Vimentin")

immune_cells <- c("B", "CD4T", "CD8T", "Treg", "NK", "DC", "Other_Immune",
                  "Macro_CD11b", "Macro_CD163", "Macro_HLADR", "Macro_Other")

# Perform tissue-specific, clinically-focused comparisons for RFS_status
cat("=== TISSUE-SPECIFIC CLINICAL INTERACTION ANALYSIS ===\n\n")

tissues <- c("TC", "IM")
tissue_results <- list()

for(tissue in tissues) {
  cat("Analyzing", tissue, "tissue...\n")
  tissue_results[[tissue]] <- focused_systematic_comparison(out_clinical, "RFS_status", tissue)
  
  if(!is.null(tissue_results[[tissue]])) {
    sig_results <- tissue_results[[tissue]]$results_table %>% filter(significant)
    cat("Significant interactions in", tissue, ":", nrow(sig_results), "\n")
    if(nrow(sig_results) > 0) {
      print(sig_results)
    }
    cat("\n")
  }
}

## Not significant results, directly to the treatment analysis

# TREATMENT-STRATIFIED CLINICAL INTERACTION ANALYSIS
# Perform treatment-stratified comparisons
cat("=== TREATMENT-STRATIFIED RFS INTERACTION ANALYSIS ===\n\n")

# For IM only
out_clinical_im <- out_clinical[out_clinical$Tissue == "IM",]
treatments <- c("Chemo", "Combo")
treatment_results <- list()

for(treatment in treatments) {
  treatment_results[[treatment]] <- treatment_stratified_comparison(out_clinical_im, treatment)
  
  if(!is.null(treatment_results[[treatment]])) {
    sig_results <- treatment_results[[treatment]]$results_table %>% filter(significant)
    cat("Significant interactions in", treatment, "group:", nrow(sig_results), "\n")
    if(nrow(sig_results) > 0) {
      print(sig_results)
    }
    cat("\n")
  }
}

# Combine results across treatments for comparison
combined_treatment_results <- bind_rows(
  lapply(treatment_results, function(x) if(!is.null(x)) x$results_table else NULL)
)

cat("=== SUMMARY: TREATMENT-SPECIFIC RFS PREDICTIVE INTERACTIONS ===\n")
significant_by_treatment <- combined_treatment_results %>% 
  filter(significant) %>%
  arrange(treatment, p_value)

# For TC only
out_clinical_tc <- out_clinical[out_clinical$Tissue == "TC",]
treatment_results <- list()

for(treatment in treatments) {
  treatment_results[[treatment]] <- treatment_stratified_comparison(out_clinical_tc, treatment)
  
  if(!is.null(treatment_results[[treatment]])) {
    sig_results <- treatment_results[[treatment]]$results_table %>% filter(significant)
    cat("Significant interactions in", treatment, "group:", nrow(sig_results), "\n")
    if(nrow(sig_results) > 0) {
      print(sig_results)
    }
    cat("\n")
  }
}

# Combine results across treatments for comparison
combined_treatment_results <- bind_rows(
  lapply(treatment_results, function(x) if(!is.null(x)) x$results_table else NULL)
)

cat("=== SUMMARY: TREATMENT-SPECIFIC RFS PREDICTIVE INTERACTIONS ===\n")
significant_by_treatment <- combined_treatment_results %>% 
  filter(significant) %>%
  arrange(treatment, p_value)

## Not significant results, directly to the treatment analysis

# ============================================================================
# FOCUSED ANALYSIS: EC_GLUT1 vs EC_EpCAM INTERACTION PATTERNS
# ============================================================================
## IM
# Perform EC interaction analysis
ec_results <- ec_interaction_analysis(interaction_summary)

# VISUALIZATION: BIDIRECTIONAL BARPLOTS FOR EC INTERACTIONS
# Create barplots for both EC types
glut1_plot <- create_ec_barplot(
  ec_results$glut1_data, 
  "EC_GLUT1", 
  "Associated with Early Recurrence (Chemotherapy)"
)

epcam_plot <- create_ec_barplot(
  ec_results$epcam_data,
  "EC_EpCAM", 
  "Associated with No Recurrence"
)

# Display the plots
pdf(file.path(figureDir,"EC_GLUT1 and EC_EpCAM interaction barplot.pdf"),width = 10,height = 12)
print(glut1_plot / epcam_plot)
dev.off()

# ============================================================================
# Display specific cell interaction
# ============================================================================
ec_subset <- c("EC_EpCAM", "EC_Vimentin", "EC_GLUT1", "EC_CAIX")
stromal_cells <- c("CAF", "SC_Collagen", "SC_Collagen_Vimentin", "SC_Vimentin", "SC_aSMA_Collagen", "SC_aSMA_Vimentin")
immune_cells <- c("B", "CD4T", "CD8T", "Treg", "NK", "Macro_CD11b", "Macro_CD163", "Macro_HLADR")

cat("Creating EC interaction chord plot...\n")

for(tissue_type in c("IM","TC")){
  for(ec_cell in ec_subset){
    # 1. Create chord plot for EC interactions
    chord_result <- create_ec_chord_plot(interaction_summary, tissue_type = tissue_type, 
                                         min_prop_significant = 0.3, 
                                         ec_cells = ec_cell, selected_targets = c(stromal_cells,immune_cells),
                                         savePath = file.path(figureDir,paste0(ec_cell," chord within ",tissue_type,".pdf")))
    
    # 2. Multiple comparisons by RFS status
    rfs_comparisons <- create_multiple_ec_comparisons(out_clinical,
                                                      ec_cell,c(stromal_cells,immune_cells),
                                                      tissue_type = tissue_type, comparison_var = "RFS_status")
    
    pdf(file.path(figureDir,paste0(ec_cell," interaction barplot within ",tissue_type,".pdf")),width = 12,height = 9)
    grid.arrange(grobs = rfs_comparisons$plots[1:length(rfs_comparisons$plots)], ncol = 4)
    dev.off()
  }
}


# ============================================================================
# SUMMARY STATISTICS TABLE
# ============================================================================

# Create summary table of interaction patterns by clinical variables
summary_table <- out_clinical %>%
  filter(!is.na(Tissue), !is.na(Treatment), !is.na(RFS_status)) %>%
  group_by(Tissue, Treatment, RFS_status) %>%
  summarise(
    n_samples = length(unique(group_by)),
    total_interactions = n(),
    significant_interactions = sum(sig, na.rm = TRUE),
    prop_significant = mean(sig, na.rm = TRUE),
    attractive_interactions = sum(interaction & sig, na.rm = TRUE),
    prop_attractive = mean(interaction[sig], na.rm = TRUE),
    mean_interaction_count = mean(ct, na.rm = TRUE),
    .groups = "drop"
  )

print("Summary of interaction patterns by clinical variables:")
print(summary_table)
write.csv(summary_table,file.path(figureDir,"interaction summary table.csv"))

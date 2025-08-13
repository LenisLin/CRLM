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
# 5. Interaction based analysis
# =============================================================================
library(scales)
library(BiocParallel)
library(dplyr)
library(circlize)
library(ggplot2)

out <- testInteractions(spe_pt, 
                        group_by = "sample_id",
                        label = "major_celltype", 
                        colPairName = "knn_20",
                        BPPARAM = SerialParam(RNGseed = 619))

saveRDS(out, file = file.path(figureDir,"PT_major_interaction_out.rds"))

# Step 1: Filter out "Other" cell type
out_filtered <- out %>% 
  as_tibble() %>%
  filter(from_label != "Other" & to_label != "Other")

# Check the filtered data
cat("Data after filtering:\n")
print(table(out_filtered$from_label))

# Step 2: Heatmap
p <- out_filtered %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("majortype_interaction_heatmap_PT.pdf", p, width = 8, height = 6)

# Step 3: Sum sigval entries across all images for each cell type pair
interaction_matrix <- out_filtered %>%
  group_by(from_label, to_label) %>%
  summarise(total_sigval = sum(sigval, na.rm = TRUE), .groups = 'drop')

# Create a square matrix for chord diagram
cell_types <- unique(c(interaction_matrix$from_label, interaction_matrix$to_label))
mat <- matrix(0, nrow = length(cell_types), ncol = length(cell_types))
rownames(mat) <- cell_types
colnames(mat) <- cell_types

# Fill the matrix
for(i in 1:nrow(interaction_matrix)) {
  from <- interaction_matrix$from_label[i]
  to <- interaction_matrix$to_label[i]
  mat[from, to] <- interaction_matrix$total_sigval[i]
}

# Step 4: Create chord plot
# Set up colors for different cell types
colors <- metadata(spe_pt)$color_vectors$major_celltype

# Create the chord diagram
pdf("major type interaction chord plot.pdf", width = 10, height = 10)
chordDiagram(mat, 
             grid.col = colors,
             transparency = 0.2,
             directional = 1,
             direction.type = c("arrows", "diffHeight"), 
             diffHeight = -0.04,
             annotationTrack = "grid",
             annotationTrackHeight = c(0.05, 0.1),
             link.arr.type = "big.arrow",
             link.sort = TRUE,
             link.largest.ontop = TRUE)

# Add cell type labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)

title("Cell-Cell Spatial Interactions (Chord Diagram)")
dev.off()

# For display in environment without PDF output
chordDiagram(mat, 
             grid.col = colors,
             transparency = 0.2,
             directional = 1,
             direction.type = c("arrows", "diffHeight"), 
             diffHeight = -0.04,
             annotationTrack = "grid",
             annotationTrackHeight = c(0.05, 0.1),
             link.arr.type = "big.arrow",
             link.sort = TRUE,
             link.largest.ontop = TRUE)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)

title("Cell-Cell Spatial Interactions (Chord Diagram)")
circos.clear()

# Step 5: Create barplot for epithelial interactions
epithelial_interactions <- out_filtered %>%
  filter(from_label == "Epithelial" | to_label == "Epithelial") %>%
  mutate(
    other_cell_type = ifelse(from_label == "Epithelial", to_label, from_label),
    interaction_direction = ifelse(from_label == "Epithelial", "from_epithelial", "to_epithelial")
  ) %>%
  group_by(other_cell_type) %>%
  summarise(total_interaction = sum(sigval, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(abs(total_interaction)))

# Create the barplot
p <- ggplot(epithelial_interactions, aes(x = total_interaction, 
                                         y = reorder(other_cell_type, total_interaction))) +
  geom_col(aes(fill = ifelse(total_interaction > 0, "Interaction", "Avoidance")),
           width = 0.7) +
  scale_fill_manual(values = c("Interaction" = "#2E86AB", "Avoidance" = "#F24236"),
                    name = "Type") +
  coord_flip() +
  labs(
    title = "Epithelial Cell Spatial Interactions",
    subtitle = "Positive values = Attraction, Negative values = Avoidance",
    x = "Cell Type",
    y = "Total Interaction Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

print(p)

# Save the plot
ggsave("epithelial_interactions_barplot.pdf", p, width = 8, height = 6)

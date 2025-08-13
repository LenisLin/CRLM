# Cholangiocyte Microenvironment Analysis for Cancer Recurrence
# Focus on PT (peritumor) tissue analysis

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
library(ComplexHeatmap)
library(circlize)

library(dplyr)
library(tidyr)

# =============================================================================
# LOAD AND VALIDATE DATA
# =============================================================================
codeSpace <- "/home/lenislin/Experiment/projects/CRLM_2025/IMC"

source(file.path(codeSpace,"BDME_analysis_functions.R"))

# Set Working Directory and Source Utilities
workDir <- "/mnt/public/lyx/IMC_HE_Merge/CRLM"
setwd(workDir)

# Define directories for saving results and Steinbock outputs
saveDir <- file.path(workDir, "results")
steinResultDir <- file.path(workDir, "steinbock")
output_dir <- file.path(workDir, "figures","6_PTME_Analysis")

if(!dir.exists(output_dir)){
  dir.create(output_dir,recursive = T)
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
# Extract metadata and expression data for easier handling
cat("Step 0: Preparing PT tissue data...\n")

# Filter for PT tissue only
pt_spe <- spe[, spe$Tissue == "PT"]

pt_cells <- as.data.frame(colData(pt_spe)) # Extract metadata
expr_data <- assay(pt_spe) # Get expression data for PT cells
coords <- spatialCoords(pt_spe) # Get spatial coordinates for PT cells  

pt_cells$Treatment_RFS <- paste(pt_cells$Treatment, # Create treatment groups
                                ifelse(pt_cells$RFS_status == 0, "NR", "R"), 
                                sep = "_")
pt_cells$sub_celltype <- ifelse(startsWith(pt_cells$sub_celltype,"EC_"),"EC",pt_cells$sub_celltype)


cat("Data prepared:\n")
cat("- Total PT cells:", nrow(pt_cells), "\n")
cat("- Cholangiocytes (Tumor = TRUE):", sum(pt_cells$Tumor, na.rm = TRUE), "\n")
cat("- Treatment groups:", table(pt_cells$Treatment_RFS), "\n")
cat("=== Starting Cholangiocyte Microenvironment Analysis ===\n\n")

rm()
gc()

# =============================================================================
# Step 1
# =============================================================================
# Step 1.1: Calculate zone sizes
cat("\nStep 1.1: Calculating zone sizes...\n")

# For each patient and zone, calculate area
zone_sizes <- calculate_zone_sizes(meta=pt_cells, coords = coords)
head(zone_sizes)

# =============================================================================
# Step 1.2: Calculate cellular composition
cat("\nStep 1.2: Calculating cellular composition...\n")

# For each zone, calculate cell type composition
zone_composition <- pt_cells %>%
  filter(!is.na(Tumor_patch) & Tumor_patch != "NA") %>%
  count(patient_id, sample_id, Tumor_patch, Treatment, RFS_status, Treatment_RFS, sub_celltype) %>%
  pivot_wider(names_from = sub_celltype, values_from = n, values_fill = 0) %>%
  rowwise() %>%
  mutate(total_cells = sum(c_across(where(is.numeric)))) %>%
  ungroup() %>%
  # Calculate proportions only for cell type columns (exclude metadata columns)
  mutate(across(where(is.numeric) & !any_of(c("patient_id", "sample_id", "Treatment", "RFS_status", "Treatment_RFS", "total_cells")), 
                ~ .x / total_cells, .names = "{.col}_prop"))

cat("Cellular composition calculated for", nrow(zone_composition), "zones\n")
head(zone_composition)

# =============================================================================
# Step 1.3: Calculate marker expression
cat("\nStep 1.3: Calculating marker expression...\n")

## select markers
markers <- c(
  "FASN","PRPS1","CA_IX","GLUT1","Ki67","HK2", ## Metabolism
  "CD127","CD80","CD274","TIGIT","CD279","CD366","CD27", ## Immune
  "Vimentin", "FAP", "VEGF" ## CAF
)

# Check which markers are available
available_markers <- intersect(markers, rownames(expr_data))
cat("Available markers:", paste(available_markers, collapse = ", "), "\n")

# Calculate mean expression for each zone
zone_expression <- as.data.frame(matrix(data = NA, nrow = 0,ncol = (6 + length(markers))))
patch_ids <- unique(na.omit(pt_cells$Tumor_patch))

for(patch_id in patch_ids){
  pt_cells_ <- pt_cells[pt_cells$Tumor_patch %in% patch_id,]
  expr_ <- expr_data[available_markers, rownames(pt_cells_)]
  expr_mean_ <- rowMeans(expr_)
  zone_expression <- rbind(zone_expression,c(
    pt_cells_$patient_id[1],pt_cells_$sample_id[1],
    patch_id,pt_cells_$Treatment[1],
    pt_cells_$RFS_status[1],pt_cells_$Treatment_RFS[1],
    expr_mean_))
}
colnames(zone_expression) <- c("patient_id","sample_id","Tumor_patch","Treatment","RFS_status","Treatment_RFS",
                               available_markers)

cat("Expression calculated for", length(unique(zone_expression$Tumor_patch)), "zones\n")
rm(patch_id, pt_cells_, expr_, expr_mean_)

# =============================================================================
# Step 3
# =============================================================================

# Step 3.1: Compare zone sizes
cat("\nStep 3.1: Comparing zone sizes (R vs NR)...\n")
comparisons <- create_comparison_groups(zone_sizes)

# Generate plots for each comparison
size_results  <- list()
size_results $all <- create_combined_size_boxplots(
  comparisons$all, 
  "All Patients", 
  file.path(output_dir, "zone_size_all_patients")
)

size_results $chemo <- create_combined_size_boxplots(
  comparisons$chemo, 
  "Chemotherapy Only", 
  file.path(output_dir, "zone_size_chemo_only")
)

size_results $combo <- create_combined_size_boxplots(
  comparisons$combo, 
  "Combination Therapy", 
  file.path(output_dir, "zone_size_combo_therapy")
)

# =============================================================================
# Step 3.2: Compare cellular composition
cat("\nStep 3.2: Comparing cellular composition (R vs NR)...\n")
comparisons <- create_comparison_groups(zone_composition)

# Calculate stats for all conditions
all_stats <- bind_rows(
  calculate_composition_stats(comparisons$all, "All Patients"),
  calculate_composition_stats(comparisons$chemo, "Chemotherapy Only"),
  calculate_composition_stats(comparisons$combo, "Combination Therapy")
) %>%
  filter(!is.na(p_value)) %>%
  group_by(condition) %>%
  mutate(
    # BH adjustment within each condition
    p_adj = p.adjust(p_value, method = "BH"),
    neg_log10_p = -log10(p_adj),
    # Cap the maximum -log10(p) value at 10
    neg_log10_p_capped = pmin(neg_log10_p, 10),
    significant = p_adj < 0.05,
    direction = case_when(
      difference > 0 & significant ~ "Higher in R",
      difference < 0 & significant ~ "Higher in NR", 
      TRUE ~ "Not significant"
    ),
    abs_difference = abs(difference),
    # Create labels for significant points
    label = ifelse(significant, cell_type, "")
  ) %>%
  ungroup()

# Create combined volcano plot
volcano_plot <- ggplot(all_stats, aes(x = difference, y = neg_log10_p_capped)) +
  geom_point(aes(color = direction, size = abs_difference), alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  ggrepel::geom_text_repel(aes(label = label),size = 3,max.overlaps = Inf,box.padding = 0.3,point.padding = 0.3) +
  scale_color_manual(
    values = c(
      "Higher in R" = "#EFC000FF",
      "Higher in NR" = "#0073C2FF", 
      "Not significant" = "gray60"
    )
  ) +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  facet_wrap(~ condition, ncol = 3) +
  labs(
    title = "Cell Composition Changes: Recurrence vs Non-recurrence",
    x = "Difference in Proportion log2(R/NR)",
    y = "-log10(adjusted p-value)",
    color = "Direction",
    caption = "p-values adjusted using Benjamini-Hochberg method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  )

# Save the plot
ggsave(file.path(output_dir, "composition_volcano_plot.pdf"), 
       volcano_plot, width = 15, height = 6, dpi = 300)


# =============================================================================
# Step 3.3: Compare marker expression
comparisons <- create_comparison_groups(zone_expression)

# Calculate stats for all conditions
expression_results <- bind_rows(
  calculate_expression_stats(comparisons$all, "All Patients"),
  calculate_expression_stats(comparisons$chemo, "Chemotherapy Only"),
  calculate_expression_stats(comparisons$combo, "Combination Therapy")
) %>%
  filter(!is.na(p_value)) %>%
  group_by(condition) %>%
  mutate(
    # BH adjustment within each condition
    p_adj = p.adjust(p_value, method = "BH"),
    neg_log10_p = -log10(p_adj),
    # Cap the maximum -log10(p) value at 10
    neg_log10_p_capped = pmin(neg_log10_p, 10),
    significant = p_adj < 0.05,
    direction = case_when(
      difference > 0 & significant ~ "Higher in R",
      difference < 0 & significant ~ "Higher in NR", 
      TRUE ~ "Not significant"
    ),
    abs_difference = abs(difference),
    # Create labels for significant points
    label = ifelse(significant, marker, "")
  ) %>%
  ungroup()

# Create combined volcano plot
volcano_plot <- ggplot(expression_results, aes(x = log2_fc, y = neg_log10_p_capped)) +
  geom_point(aes(color = direction, size = abs_difference), alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  ggrepel::geom_text_repel(aes(label = label),size = 3,max.overlaps = Inf,box.padding = 0.3,point.padding = 0.3) +
  scale_color_manual(
    values = c(
      "Higher in R" = "#EFC000FF",
      "Higher in NR" = "#0073C2FF", 
      "Not significant" = "gray60"
    )
  ) +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  facet_wrap(~ condition, ncol = 3) +
  labs(
    title = "Marker Expression Changes: Recurrence vs Non-recurrence",
    x = "Difference in Proportion log2(R/NR)",
    y = "-log10(adjusted p-value)",
    color = "Direction",
    caption = "p-values adjusted using Benjamini-Hochberg method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  )
ggsave(file.path(output_dir, "marker_expression_volcano_plot.pdf"), 
       volcano_plot, width = 10, height = 8, dpi = 300)


# ===============================================================================
# Step 4: Defined two peritumor cholangiocyte microenvironments (PCME)
# PCME-I: immuno-competent / PCME-S: stromal-dominant
# ===============================================================================
# Step 4.1: Calculate the Expression of cholangiocytes within different PCME

# Calculate CAIX expression in epithelial cells by zone
markers <- c(
  "FASN","PRPS1","CA_IX","GLUT1","Ki67","HK2", ## Metabolism
  "CD274", ## Immune
  "VEGF" ## CAF
)

epithelial_exp_mean <- calculate_ec_mean_expression(pt_cells,expr_data,markers)
# pcme_classification <-classify_by_caix_only(zone_composition, epithelial_exp_mean)
pcme_classification <- classify_by_dual_characteristics(zone_composition, epithelial_exp_mean)
validate_pcme_classification(pcme_classification)

# Step 4.2: Aggregate PCME Zones to Patient-Level Scores and Extract clinical features
# ===============================================================================
patient_scores <- aggregate_pcme_to_patient(pcme_classification, zone_sizes)

patient_data <- patient_scores

# Clinical features to extract
clinical_vars <- c("RFS_time","Age", "Gender", "KRAS_mutation", "TBS", "Differential_grade", 
                   "T_stage", "CRLM_number", "CRLM_size", "CEA", "CA199")

for(var_ in clinical_vars){
  patient_data[var_] <- pt_cells[match(patient_data$patient_id,pt_cells$patient_id),var_]
}

threshold_result  <- find_optimal_threshold(patient_data, score_var = "log_ratio_weighted", 
                       outcome_var = "RFS_status", time_var = "RFS_time", 
                       method = "youden")

analysis_data <- threshold_result$analysis_data

# Step 4.3: Draw KM curves
# ===============================================================================
# Overall survival analysis
fit_overall <- survfit(Surv(RFS_time, RFS_status) ~ risk_group, data = analysis_data)

# Create overall KM plot
km_plot_overall <- ggsurvplot(
  fit_overall,
  data = analysis_data,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  tables.height = 0.2,
  tables.theme = theme_cleantable(),
  palette = c("#E7B800", "#2E9FDF"),
  title = paste("Kaplan-Meier Curves by PCME Score"),
  subtitle = paste("Threshold:", round(threshold_result$threshold, 3)),
  xlab = "Time (months)",
  ylab = "Recurrence-free Survival",
  legend.title = "PCME Risk Group",
  legend.labs = c("High-risk (PCME-S dominant)", "Low-risk (PCME-I dominant)")
)

# Save overall plot
ggsave(file.path(output_dir, "km_curve_overall_pcme.pdf"), 
       km_plot_overall$plot, width = 10, height = 8, dpi = 300)

# Treatment-stratified analysis
Chemo_data <- 
  analysis_data %>% filter(Treatment == "Chemo")
fit_Chemo <- survfit(Surv(RFS_time, RFS_status) ~ risk_group, data = Chemo_data)
km_plot_Chemo <- ggsurvplot(
  fit_Chemo,
  data = Chemo_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  tables.height = 0.2,
  palette = c("#E7B800", "#2E9FDF"),
  title = paste("KM Curves - Chemotherapy"),
  xlab = "Time (months)",
  ylab = "Recurrence-free Survival"
)

Combo_data <- 
  analysis_data %>% filter(Treatment == "Combo")
fit_Combo <- survfit(Surv(RFS_time, RFS_status) ~ risk_group, data = Combo_data)
km_plot_Combo <- ggsurvplot(
  fit_Combo,
  data = Combo_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  tables.height = 0.2,
  palette = c("#E7B800", "#2E9FDF"),
  title = paste("KM Curves - Combination Therapy"),
  xlab = "Time (months)",
  ylab = "Recurrence-free Survival"
)

combined_plot <- arrange_ggsurvplots(list(km_plot_overall, km_plot_Chemo, km_plot_Combo), ncol = 3, nrow = 1)
pdf(file.path(output_dir, paste0("km_curve_pcme.pdf")),width = 20,height = 8)
print(combined_plot)
dev.off()

# Step 4.4: Draw AUC Curve
# ===============================================================================
# Prepare data for ROC analysis
clinical_vars <- clinical_vars[-1] ## Remove RFS_time
score_var <- "log_ratio_weighted"
roc_data <- patient_data %>%
  select(patient_id, RFS_status, Treatment, all_of(c(score_var, clinical_vars))) %>%
  filter(!is.na(RFS_status) & !is.na(.data[[score_var]]))

# Individual ROC curves
roc_results <- list()

# PCME score alone
roc_pcme <- roc(roc_data$RFS_status, roc_data[[score_var]], direction = "<", quiet = TRUE)
roc_results[["PCME_Score"]] <- roc_pcme

# Clinical variables (where possible)
for(var in clinical_vars) {
  if(var %in% colnames(roc_data) && is.numeric(roc_data[[var]])) {
    var_data <- roc_data %>% filter(!is.na(.data[[var]]))
    if(nrow(var_data) > 10) {
      roc_var <- roc(var_data$RFS_status, var_data[[var]], direction = "<", quiet = TRUE)
      roc_results[[var]] <- roc_var
    }
  }
}

# Create comparison plot
auc_values <- sapply(roc_results, auc)
auc_df <- data.frame(
  Variable = names(auc_values),
  AUC = as.numeric(auc_values)
) %>%
  arrange(desc(AUC))

# Plot individual AUCs
auc_plot <- ggplot(auc_df, aes(x = reorder(Variable, AUC), y = AUC)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_text(aes(label = round(AUC, 3)), hjust = -0.1) +
  coord_flip() +
  labs(
    title = "ROC Analysis: PCME Score vs Clinical Features",
    x = "Variables",
    y = "Area Under Curve (AUC)",
    subtitle = paste("Primary outcome: RFS_status")
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "roc_comparison_auc.pdf"), 
       auc_plot, width = 10, height = 8, dpi = 300)

# Combined ROC curve plot (top variables)
# top_vars <- c("PCME_Score","TBS","CRLM_size","Differential_grade","CEA","CA199")
top_vars <- auc_df$Variable

for(i in 1:length(top_vars)){
  pdf(file.path(output_dir, "AUC_each_feature", paste0(top_vars[i]," roc curves.pdf")), width = 6, height = 4.5)
  plot(roc_results[[top_vars[i]]], main = paste0(top_vars[i]),
       col = "blue", lwd = 4, 
       print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, 
       grid.col = "darkgray", grid.lwd = 1.5, auc.polygon.col = "skyblue")
  dev.off()
}


# Step 4.5: Draw Forest plot
# ===============================================================================
cat("\n=== Multivariate Analysis ===\n")

score_var = "log_ratio_weighted"

# Prepare complete case data
model_vars <- c("RFS_time", "RFS_status", "Treatment", score_var, clinical_vars)
model_data <- patient_data %>%
  select(all_of(model_vars)) %>%
  filter(complete.cases(.))

cat("Complete cases for multivariate analysis:", nrow(model_data), "\n")

# Cox proportional hazards models
clinical_formula <- paste(clinical_vars, collapse = " + ")

# Model 2: PCME score + clinical variables
formula2 <- as.formula(paste("Surv(RFS_time, RFS_status) ~", score_var, "+", clinical_formula))
model2 <- coxph(formula2, data = model_data)

# Model 3: PCME score + treatment + clinical variables
formula3 <- as.formula(paste("Surv(RFS_time, RFS_status) ~", score_var, 
                             "+ Treatment +", clinical_formula))
model3 <- coxph(formula3, data = model_data)

# Model comparison
model_comparison <- list(
  Model2 = list(
    formula = formula2,
    model = model2,
    aic = AIC(model2),
    concordance = summary(model2)$concordance[1]
  ),
  Model3 = list(
    formula = formula3,
    model = model3,
    aic = AIC(model3),
    concordance = summary(model3)$concordance[1]
  )
)

# Extract results for plotting
model_summary <- data.frame(
  Model = c("PCME + Treatment", "PCME + Treatment + Clinical"),
  AIC = sapply(model_comparison, function(x) x$aic),
  Concordance = sapply(model_comparison, function(x) x$concordance),
  stringsAsFactors = FALSE
)

# Create forest plot for final model
model3_summary <- summary(model3)
coef_df <- data.frame(
  Variable = rownames(model3_summary$coefficients),
  HR = exp(model3_summary$coefficients[, "coef"]),
  Lower_CI = exp(model3_summary$coefficients[, "coef"] - 1.96 * model3_summary$coefficients[, "se(coef)"]),
  Upper_CI = exp(model3_summary$coefficients[, "coef"] + 1.96 * model3_summary$coefficients[, "se(coef)"]),
  P_value = model3_summary$coefficients[, "Pr(>|z|)"]
) %>%
  mutate(
    Significant = P_value < 0.05,
    Variable = gsub(score_var, "PCME Score", Variable)
  )

generate_enhanced_forest_plot(cox_results = coef_df, tissue_name = "PT", 
                              save_path = file.path(output_dir,"Multi-variable COX of PCME score.pdf"), show_all = TRUE, 
                              max_variables = 20, alpha_level = 0.05)

cat("=== MULTIVARIATE COX REGRESSION RESULTS ===\n\n")

cat("\n\nMODEL 2: PCME Score + Treatment\n")
cat("AIC:", AIC(model2), "| C-index:", summary(model2)$concordance[1], "\n")
print(summary(model2))

cat("\n\nMODEL 3: Full Model\n")
cat("AIC:", AIC(model3), "| C-index:", summary(model3)$concordance[1], "\n")
print(summary(model3))


# ===============================================================================
# Step 5: Expression Correlation Analysis of PCME-S
# ===============================================================================

# 5.1: Cholangiocyte Marker Expression vs PCME-S Patch Size
# ===============================================================================

# Filter for PCME-S zones only
pcme_s_zones <- pcme_classification %>%
  filter(PCME_binary == "Stromal-fibrotic") %>%
  select(patient_id, sample_id, Tumor_patch)

# Merge with expression data
merged_data <- pcme_s_zones %>%
  left_join(epithelial_exp_mean, by = c("patient_id", "sample_id", "Tumor_patch"))
merged_data$log2_patch_size <- log2(merged_data$patch_size)

if(nrow(merged_data) == 0) {
  cat("No data available for size-expression correlation analysis\n")
  return(NULL)
}

# Create discrete size bins based on log2(patch_size)
merged_data <- merged_data %>%
  mutate(
    size_bin = cut(log2_patch_size, 
                   breaks = 5, 
                   labels = c("[2,3)", "[3,4)", "[4,5)", "[5,6)", "[6,7)"),
                   include.lowest = TRUE)
  ) %>%
  filter(!is.na(size_bin))

cat("Size bins created with distribution:\n")
print(table(merged_data$size_bin))

# Function to test trend across size bins
test_size_trend <- function(data, marker_name) {
  marker_data <- data %>% select("patient_id", "sample_id", "Tumor_patch","log2_patch_size",marker_name)
  
  if(nrow(marker_data) < 10) {
    return(data.frame(
      marker = marker_name,
      rho = NA,
      p_value = NA,
      n_observations = nrow(marker_data)
    ))
  }
  
  # Convert size_bin to numeric for correlation
  marker_data$log2_patch_size <- as.numeric(marker_data$log2_patch_size)
  cor_test <- cor.test(marker_data$expression, marker_data$log2_patch_size, method = "spearman")
  
  data.frame(
    marker = marker_name,
    rho = cor_test$estimate,
    p_value = cor_test$p.value,
    n_observations = nrow(marker_data)
  )
}

# Test trends for all markers
trend_results <- map_dfr(markers, ~ test_size_trend(merged_data, .x)) %>%
  filter(!is.na(rho)) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significant = p_adj <= 0.05,
    direction = case_when(
      significant & rho > 0 ~ "Positive",
      significant & rho < 0 ~ "Negative", 
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(desc(abs(rho)))

cat("Trend analysis results:\n")
print(trend_results)

# Create boxplots for each marker
plot_list <- list()

for(marker_name in markers) {
  marker_data <- merged_data %>% select("patient_id", "sample_id", "Tumor_patch","size_bin",marker_name)
  colnames(marker_data)[ncol(marker_data)] <- "expression"
  
  if(nrow(marker_data) < 5) next
  
  # Get trend info for this marker
  trend_info <- trend_results[trend_results$marker == marker_name, ]
  
  if(nrow(trend_info) == 0) {
    trend_info <- data.frame(rho = NA, p_adj = NA)
  }
  
  p <- ggplot(marker_data, aes(x = size_bin, y = expression)) +
    geom_boxplot(aes(fill = size_bin), alpha = 0.7, outlier.alpha = 0.5) +
    geom_jitter(aes(color = "grey75"), width = 0.2, alpha = 0.6, size = 1) +
    scale_fill_viridis_d(option = "plasma", alpha = 0.7) +
    labs(
      title = paste(marker_name, "in PCME-S Zones"),
      subtitle = paste("Spearman ρ =", round(trend_info$rho, 3), 
                       ", p =", format(trend_info$p_value, digits = 3)),
      x = "PCME-S Patch Size Category",
      y = paste(marker_name, "Expression"),
      fill = "Size Category",
      color = "Treatment & Status"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  plot_list[[marker_name]] <- p
}

# Combine all plots
if(length(plot_list) > 0) {
  combined_plot <- wrap_plots(plot_list, ncol = 3)
  
  # Save combined plot
  ggsave(file.path(output_dir, "cholangiocyte_expression_by_pcme_s_size.pdf"), 
         combined_plot, width = 18, height = ceiling(length(plot_list)/3) * 4, dpi = 300)
  
  # Create trend summary plot
  if(nrow(trend_results) > 0) {
    trend_plot <- ggplot(trend_results, aes(x = reorder(marker, rho), y = rho)) +
      geom_col(alpha = 0.7) +
      geom_text(aes(label = paste("p =", format(p_value, digits = 3))), 
                hjust = ifelse(trend_results$rho > 0, -0.1, 1.1), size = 3) +
      scale_fill_manual(values = c("Positive" = "#e74c3c", "Negative" = "#3498db", "Not significant" = "gray60")) +
      coord_flip() +
      labs(
        title = "Cholangiocyte Marker Trends with PCME-S Size",
        x = "Markers",
        y = "Spearman Correlation (ρ)",
        fill = "Trend Direction",
        caption = "p-values adjusted using Benjamini-Hochberg method"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave(file.path(output_dir, "cholangiocyte_size_trend_summary.pdf"), 
           trend_plot, width = 10, height = 8, dpi = 300)
  }
}

# Save trend results
write.csv(trend_results, file.path(output_dir, "patch_size_marker_trends.csv"), row.names = FALSE)

# 5.2: PCME-I vs PCME-S Marker Expression by Cell Subtypes
# Get PCME classification for each zone
pcme_zones <- pcme_classification %>%
  filter(PCME_binary %in% c("Immune-permissive", "Stromal-fibrotic")) %>%
  select(patient_id, sample_id, Tumor_patch, PCME_binary)

# Merge PT cells with PCME classification
cells_with_pcme <- pt_cells
cells_with_pcme <- cells_with_pcme[!is.na(cells_with_pcme$Tumor_patch) ,]
cells_with_pcme$PCME_binary <- pcme_zones$PCME_binary[match(cells_with_pcme$Tumor_patch,pcme_zones$Tumor_patch)]
cells_with_pcme <- cells_with_pcme[!is.na(cells_with_pcme$PCME_binary) ,]

cat("Cells with PCME classification:\n")
cat("- Total cells:", nrow(cells_with_pcme), "\n")
cat("- PCME distribution:\n")
print(table(cells_with_pcme$PCME_binary))
cat("- Cell subtype distribution:\n")
print(table(cells_with_pcme$sub_celltype))

# Get expression data for these cells
expr_data <- assay(pt_spe)[, match(rownames(cells_with_pcme),colnames(pt_spe))]

# Check marker availability
available_markers <- intersect(markers, rownames(expr_data))
cat("Available markers:", paste(available_markers, collapse = ", "), "\n")

if(length(available_markers) == 0) {
  cat("No markers available in expression data\n")
}

# Get all cell subtypes with sufficient representation
valid_subtypes <- cells_with_pcme %>%
  group_by(sub_celltype, PCME_binary) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = PCME_binary, values_from = n, values_fill = 0) %>%
  filter(`Immune-permissive` >= 5 & `Stromal-fibrotic` >= 5) %>%
  pull(sub_celltype)

cat("Valid subtypes for analysis:", length(valid_subtypes), "\n")
cat("Subtypes:", paste(valid_subtypes, collapse = ", "), "\n")

# Run comparison for all valid subtypes
available_markers <- rownames(spe)
all_results <- list()

for(celltype_ in valid_subtypes){
  result_ <- compare_pcme_markers(cells_with_pcme,expr_data,celltype_,available_markers) %>%
    filter(!is.na(p_value))
  
  all_results[[celltype_]] <- result_
}
all_results <- do.call(rbind, all_results)
rownames(all_results) <- 1:nrow(all_results)

if(nrow(all_results) == 0) {
  cat("No valid comparisons could be performed\n")
}

# Adjust p-values and classify results
final_results <- all_results %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    neg_log10_p = -log10(p_adj),
    neg_log10_p_capped = pmin(neg_log10_p, 10),  # Cap for visualization
    significant = p_adj <= 0.05,
    direction = case_when(
      significant & log2_fc > 0 ~ "Higher in PCME-S",
      significant & log2_fc < 0 ~ "Higher in PCME-I",
      TRUE ~ "Not significant"
    ),
    abs_log2_fc = abs(log2_fc)
  ) %>%
  arrange(p_adj)

cat("Significant comparisons:", sum(final_results$significant), "out of", nrow(final_results), "\n")

# Create heatmap data
heatmap_data <- final_results %>%
  select(cell_subtype, marker, direction, neg_log10_p_capped, log2_fc, significant)

# Create heatmap
create_pcme_heatmap(heatmap_data, output_dir)

# Save detailed results
write.csv(final_results, file.path(output_dir, "pcme_subtype_marker_comparison.csv"), row.names = FALSE)

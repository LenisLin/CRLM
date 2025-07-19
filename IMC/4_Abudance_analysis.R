# =============================================================================
# ABUNDANCE_ANALYSIS.R - Main Analysis Script (Line-by-Line Approach)
# =============================================================================

# =============================================================================
# LOAD AND VALIDATE DATA
# =============================================================================

codeSpace <- "/home/lenislin/Experiment/projects/CRLM_2025/IMC"
source(file.path(codeSpace,"Abudance_analysis_functions.R"))

# Set Working Directory and Source Utilities
workDir <- "/mnt/public/lyx/IMC_HE_Merge/CRLM"
setwd(workDir)

# Define directories for saving results and Steinbock outputs
saveDir <- file.path(workDir, "results")
steinResultDir <- file.path(workDir, "steinbock")
figureDir <- file.path(workDir, "figures","2_Abundance")

if(!dir.exists(figureDir)){
  dir.create(figureDir,recursive = T)
}

# Load your SpatialExperiment object
spe <- readRDS(file.path(saveDir,"subanno_spe_0714.rds"))

# Combine new clinical column
if(F){
  clinical_df <- read.csv("/home/lenislin/Experiment/projects/CRLM_2025/IMC/clinical.csv",row.names = 1)
    
  ## Add Treatment strategy
  spe$Treatment <- NA
  spe$Treatment <- clinical_df$Treatment[match(spe$patient_id,rownames(clinical_df))]
  
  ## Modify the Tissue column
  tissue_ <- as.character(spe$Tissue)
  tissue_ <- ifelse(tissue_ == "CT", "TC", tissue_)
  tissue_ <- ifelse(tissue_ == "TAT", "PT", tissue_)
  spe$Tissue <- tissue_
  
  saveRDS(spe, file.path(saveDir,"subanno_spe_0714.rds"))
}

# =============================================================================
# CONFIGURATION SECTION - MODIFY THESE PARAMETERS AS NEEDED
# =============================================================================

# Data structure configuration
celltype_col <- "sub_celltype"           # Cell type/subtype column
majortype_col <- "major_celltype"        # Cell type/subtype column
tissue_col <- "Tissue"                   # Tissue type column  
patient_col <- "patient_id"              # Patient ID column
roi_col <- "sample_id"                   # ROI/sample ID column
treat_col <- "Treatment"                 # Treatment Stratagy column

# Clinical variables configuration
discrete_clinical_vars <- c(
  "RFS_status",                          # Relapse-free survival status
  "Treatment",
  "Gender",                              # Patient gender
  "KRAS_mutation",                       # KRAS mutation status
  "BRAF_mutation",                       # BRAF mutation status
  "Differential_grade",                  # Differentiation grade
  "T_stage",                             # T stage
  "Lymph_positive"                       # Lymph node positivity
)

continuous_clinical_vars <- c(
  "RFS_time",                            # Relapse-free survival time
  "Age",                                 # Patient age
  "TBS",                                 # Tumor burden score
  "CRLM_number",                         # Number of CRLM
  "CRLM_size",                           # Size of CRLM
  "CEA",                                 # CEA levels
  "CA199"                                # CA19-9 levels
)

# Survival analysis variables
survival_time_col <- "RFS_time"          # Time to event
survival_event_col <- "RFS_status"       # Event indicator (0/1)

# Analysis scope settings
tissues_to_analyze <- c("TC", "IM", "PT")
include_unknown_celltypes <- FALSE       # Include unknown/unclassified cell types
min_fraction_threshold <- 0.00001          # Minimum fraction to consider (0.01%)
aggregate_level <- "roi"                 # "roi" or "patient" level analysis

# Statistical settings
p_adjust_method <- "BH"                  # "BH", "bonferroni", "holm", etc.
alpha_level <- 0.05                      # Significance threshold
test_method <- "wilcox.test"             # "wilcox.test" or "t.test"
min_cells_per_group <- 5                 # Minimum samples per group
fc_threshold <- 1.2                      # Fold change threshold for significance
correlation_method <- "spearman"         # "pearson", "spearman", "kendall"
min_observations <- 10                   # Minimum observations for correlation

# Plotting settings
color_palette <- metadata(spe)$color_vector
plot_dpi <- 300                          # Plot resolution
plot_format <- "pdf"                     # Output format: "pdf", "png", "svg"
text_size <- 12                          # Base text size
title_size <- 16                         # Title text size
axis_text_angle <- 45                    # X-axis text rotation

# Output settings
output_base_dir <- figureDir             # Base output directory
save_intermediate_results <- TRUE        # Save intermediate results
verbose_output <- TRUE                   # Print progress messages
create_subdirectories <- TRUE            # Create subdirectories by analysis type

# Validate input data (uncomment when you have data loaded)
validate_spe_input(spe, config = list(
  celltype = celltype_col,
  tissue = tissue_col,
  patient = patient_col,
  roi = roi_col,
  treatment = treat_col,
  verbose = verbose_output,
  discrete_clinical = discrete_clinical_vars,
  continuous_clinical = continuous_clinical_vars,
  survival = list(time = survival_time_col, event = survival_event_col)
))

# =============================================================================
# SETUP OUTPUT DIRECTORIES
# =============================================================================

# Create subdirectories
if (create_subdirectories) {
  subdirs <- c(
    "01_data_description",
    "02_discrete_clinical", 
    "03_continuous_clinical",
    "04_survival_analysis",
    "05_summary"
  )
  
  for (subdir in subdirs) {
    full_path <- file.path(output_base_dir, subdir)
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = TRUE)
    }
  }
}

if (verbose_output) {
  cat("✓ Output directories created in:", output_base_dir, "\n")
}

# =============================================================================
# DATA PROCESSING
# =============================================================================

if (verbose_output) {
  cat("\n", rep("=", 60), "\n")
  cat("STEP 0: DATA PROCESSING\n")
  cat(rep("=", 60), "\n")
}

# Create configuration object for functions
config <- list(
  columns = list(
    celltype = celltype_col,
    tissue = tissue_col,
    patient = patient_col,
    roi = roi_col,
    treatment = treat_col,
    discrete_clinical = discrete_clinical_vars,
    continuous_clinical = continuous_clinical_vars,
    survival = list(time = survival_time_col, event = survival_event_col)
  ),
  scope = list(
    tissues = tissues_to_analyze,
    include_unknown = include_unknown_celltypes,
    min_fraction_threshold = min_fraction_threshold
  ),
  stats = list(
    p_adjust_method = p_adjust_method,
    alpha_level = alpha_level,
    test_method = test_method,
    min_cells_per_group = min_cells_per_group,
    fc_threshold = fc_threshold,
    correlation_method = correlation_method,
    min_observations = min_observations
  ),
  plots = list(
    dpi = plot_dpi,
    format = plot_format,
    color_palette = color_palette,
    text_size = text_size,
    title_size = title_size,
    axis_text_angle = axis_text_angle
  ),
  output = list(
    base_dir = output_base_dir,
    save_intermediate = save_intermediate_results,
    verbose = verbose_output
  )
)

# Process SpatialExperiment data
analysis_data <- process_spe_data(spe, config, aggregate_level)

# Save processed data if requested
if (save_intermediate_results) {
  output_file <- file.path(output_base_dir, "processed_abundance_data.csv")
  write.csv(analysis_data, output_file, row.names = FALSE)
  if (verbose_output) {
    cat("✓ Processed data saved to:", output_file, "\n")
  }
}

# =============================================================================
# PART 1: DATA DESCRIPTION
# =============================================================================

if (verbose_output) {
  cat("\n", rep("=", 60), "\n")
  cat("PART 1: DATA DESCRIPTION\n")
  cat(rep("=", 60), "\n")
}

# Set output directory for this section
desc_output_dir <- file.path(output_base_dir, "01_data_description")

# 1.1 Cell type composition pie charts
if (verbose_output) cat("Generating cell type composition pie charts...\n")

pie_plot <- plot_celltype_pie_charts(
  fractions_df = analysis_data,
  config = config,
  save_path = file.path(desc_output_dir, "celltype_composition_pie_charts.pdf"),
  show_percentages = TRUE,
  min_percentage_label = 5,
  width = 20, height = 8
)

# 1.2 Stacked abundance bar plots
if (verbose_output) cat("Generating abundance stacked bar plots...\n")

bar_plot <- plot_stacked_abundance_bars(
  fractions_df = analysis_data,
  config = config,
  save_path = file.path(desc_output_dir, "abundance_stacked_bars.pdf"),
  order_by = "tissue",
  show_tissue_separator = TRUE,
  height = 7.5,width = 15
)

# 1.3 Calculate the entropy
# Add cellular diversity analysis
if (verbose_output) cat("Calculating cellular diversity metrics...\n")

# Calculate entropy and Simpson index
analysis_data <- calculate_cellular_diversity(
  analysis_data = analysis_data,
  config = config,
  add_scaled = TRUE,
  verbose = verbose_output
)

# Analyze diversity by clinical variables
for(clinical_var in c("RFS_status","Treatment")){
plot_cellular_diversity(
  analysis_data = analysis_data,
  clinical_var = clinical_var,
  config = config,
  save_path = file.path(desc_output_dir, paste0("diversity_", clinical_var, ".pdf")),
  use_scaled = TRUE
)
}
# Save enhanced data with diversity metrics
if (save_intermediate_results) {
  write.csv(analysis_data, 
            file.path(desc_output_dir, "analysis_data_with_diversity.csv"),
            row.names = FALSE)
}

# 1.4 Clinical features heatmap
if (verbose_output) cat("Generating clinical features heatmap...\n")

heatmap_plot <- plot_clinical_heatmap(
  fractions_df = analysis_data,
  config = config,
  save_path = file.path(desc_output_dir, "clinical_features_heatmap.pdf"),
  order_by = "event",              # Order patients by survival status
  cluster_rows = FALSE,            # Don't cluster patients
  cluster_columns = FALSE,          # Cluster cell types
  show_row_names = TRUE,           # Show patient IDs
  show_column_names = TRUE,        # Show cell type names
  width = 15, height = 9
)

# 1.5 KM plot for patients with different treatment
if (verbose_output) cat("Generating KM curve for different treatment...\n")

# Get information for plot
plotdf <- analysis_data[match(unique(analysis_data$patient_id),analysis_data$patient_id),]
plotdf <- plotdf[,c(config$columns$treatment, unname(unlist(config$columns$survival)))]
colnames(plotdf) <- c("treatment", names(unlist(config$columns$survival)))
plotdf <- na.omit(plotdf)

# Fit survival model with error checking
surv_object <- Surv(time = plotdf$time, event = plotdf$event)
fit <- survfit(surv_object ~ treatment, data = plotdf)

# Generate Kaplan-Meier plot
km_plot <- ggsurvplot(
  fit = fit,
  data = plotdf,
  pval = TRUE,                    # Add p-value
  conf.int = FALSE,                # Add confidence intervals
  risk.table = TRUE,              # Add risk table
  risk.table.col = "strata",      # Color risk table by groups
  linetype = "strata",            # Different line types for groups
  surv.median.line = "hv",        # Add median survival lines
  ggtheme = theme_bw(),           # Clean theme
  palette = c("#f39b7fff", "#8491b4ff"), # Custom colors
  title = "Recurrence-Free Survival by Cellular compostion Group",
  xlab = "Time (months)",
  ylab = "Survival probability",
  legend.title = "Treatment strategy",
  legend.labs = c("Chemo", "Combo")
)
pdf(file = file.path(desc_output_dir, "km_curve_for_treatment.pdf"), height = 6,width = 8)
print(km_plot)
dev.off()

# =============================================================================
# PART 2: DISCRETE CLINICAL ANALYSIS
# =============================================================================

if (verbose_output) {
  cat("\n", rep("=", 60), "\n")
  cat("PART 2: DISCRETE CLINICAL ANALYSIS\n")
  cat(rep("=", 60), "\n")
}

# Set output directory for this section
discrete_output_dir <- file.path(output_base_dir, "02_discrete_clinical")

# 2.1 Perform group comparisons
if (verbose_output) cat("Performing group comparisons for discrete variables...\n")

group_comparisons <- calculate_group_comparisons(
  fractions_df = analysis_data,
  config = config,
  test_method = test_method,
  effect_size = TRUE
)

# Save comparison results
if (!is.null(group_comparisons) && save_intermediate_results) {
  write.csv(
    group_comparisons,
    file.path(discrete_output_dir, "group_comparison_results.csv"),
    row.names = FALSE
  )
}

# 2.2 Generate bubble plot
if (!is.null(group_comparisons)) {
  if (verbose_output) cat("Generating abundance bubble plot...\n")

  bubble_plot <- plot_abundance_bubble_chart(
    stats_df = group_comparisons,
    config = config,
    save_path = file.path(discrete_output_dir, "abundance_bubble_plot.pdf"),
    fc_limits = c(0.5, 2), show_ns = TRUE,
    width = 18, height = 5
  )
}

# 2.3 Generate individual boxplots for each clinical variable
if (verbose_output) cat("Generating individual boxplots...\n")

# Loop through each discrete clinical variable
for (clinical_var in discrete_clinical_vars) {
  if (verbose_output) cat("  - Processing", clinical_var, "...\n")

  # Generate boxplots for this variable across all tissues
  plot_abundance_boxplots(
    fractions_df = analysis_data,
    clinical_var = clinical_var,
    config = config,
    save_path = file.path(discrete_output_dir, paste0("boxplots_", clinical_var, ".pdf")),
    test_method = test_method,
    show_significance = TRUE, show_points = TRUE,
    height = 8, width = 5
  )
}

# Add treatment information for compare
plot_abundance_boxplots_with_split(
  fractions_df = analysis_data,
  clinical_var = "RFS_status",
  split_var = "Treatment",
  config = config,
  save_path = file.path(discrete_output_dir, paste0("boxplots compare RFS_status with different Treatment.pdf")),
  test_method = test_method,
  show_significance = TRUE, show_points = TRUE,
  height = 8, width = 6
)

# =============================================================================
# PART 3: CONTINUOUS CLINICAL ANALYSIS
# =============================================================================

if (verbose_output) {
  cat("\n", rep("=", 60), "\n")
  cat("PART 3: CONTINUOUS CLINICAL ANALYSIS\n")
  cat(rep("=", 60), "\n")
}

# Set output directory for this section
continuous_output_dir <- file.path(output_base_dir, "03_continuous_clinical")

# 3.1 Calculate correlations
if (verbose_output) cat("Calculating correlations with continuous variables...\n")

correlations <- calculate_abundance_correlations(
  fractions_df = analysis_data,
  config = config,
  method = correlation_method,
  min_observations = min_observations
)

# 3.2 Generate correlation heatmaps for each tissue
if (length(correlations) > 0) {
  for (tissue in names(correlations)) {
    if (verbose_output) cat("  - Generating correlation heatmap for", tissue, "...\n")

    plot_correlation_heatmap(
      correlation_results = correlations[[tissue]], # Pass the list for the tissue
      tissue_name = tissue,
      config = config,
      save_path = file.path(continuous_output_dir, paste0("correlation_heatmap_", tissue, ".pdf")),
      height = 8, width = 8
    )
  }
}

# 3.3 Generate scatter plots for significant correlations
if (verbose_output) cat("Generating scatter plots for significant correlations...\n")

for (tissue in tissues_to_analyze) {
  if (tissue %in% names(correlations)) {
    if (verbose_output) cat("  - Processing", tissue, "tissue...\n")

    scatter_plots <- generate_correlation_scatter_plots(
      fractions_df = analysis_data,
      correlation_matrix = correlations[[tissue]],
      tissue_name = tissue,
      config = config,
      output_dir = continuous_output_dir,
      significance_threshold = alpha_level,
      width = 8,height = 6
    )
  }
}

# 3.4 Cell type correlations within tissues
if (verbose_output) cat("Analyzing cell type correlations within tissues...\n")

celltype_correlations <- analyze_celltype_correlations(
  fractions_df = analysis_data,
  config = config,
  output_dir = continuous_output_dir,
  method = correlation_method,
  min_correlation = 0.0
)

# 3.5 Generate correlation summary table
if (verbose_output) cat("Creating correlation summary table...\n")

correlation_summary <- create_correlation_summary(
  correlations = correlations,
  config = config,
  save_path = file.path(continuous_output_dir, "correlation_summary.csv")
)


# =============================================================================
# PART 4: SURVIVAL ANALYSIS
# =============================================================================

if (verbose_output) {
  cat("\n", rep("=", 60), "\n")
  cat("PART 4: SURVIVAL ANALYSIS\n")
  cat(rep("=", 60), "\n")
}

# Set output directory for this section
survival_output_dir <- file.path(output_base_dir, "04_survival_analysis")

# 4.1 Prepare survival data (aggregate to patient level if needed)
if (verbose_output) cat("Preparing survival data...\n")

survival_data <- prepare_survival_data(
  fractions_df = analysis_data,
  clinical_df = NULL,  # Already merged
  config = config,
  aggregation_method = "mean"
)

# 4.2 Univariate Cox regression for each cell type
if (verbose_output) cat("Performing univariate Cox regression...\n")

univariate_cox_results <- list()

for (tissue in tissues_to_analyze) {
  if (verbose_output) cat("  - Processing", tissue, "tissue...\n")

  tissue_survival_data <- survival_data[survival_data[[tissue_col]] == tissue, ]

  if (nrow(tissue_survival_data) >= 20) {  # Minimum for Cox regression
    univariate_results <- perform_univariate_cox(
      tissue_data = tissue_survival_data,
      config = config, scale_variables = TRUE
    )

    univariate_cox_results[[tissue]] <- univariate_results

    # Save individual results
    if (save_intermediate_results) {
      write.csv(
        univariate_results,
        file.path(survival_output_dir, paste0("univariate_cox_", tissue, ".csv")),
        row.names = FALSE
      )
    }
  }
}

# 4.3 Multivariate Cox regression (include clinical variables)
if (verbose_output) cat("Performing multivariate Cox regression...\n")

multivariate_cox_results <- list()

for (tissue in tissues_to_analyze) {
  if (verbose_output) cat("  - Processing", tissue, "tissue...\n")

  tissue_survival_data <- survival_data[survival_data[[tissue_col]] == tissue, ]

  if (nrow(tissue_survival_data) >= 20) {
    multivariate_results <- perform_multivariate_cox(
      tissue_data = tissue_survival_data,
      config = config,
      include_clinical = TRUE
    )

    multivariate_cox_results[[tissue]] <- multivariate_results

    # Save individual results
    if (save_intermediate_results) {
      write.csv(
        multivariate_results$results,
        file.path(survival_output_dir, paste0("multivariate_cox_", tissue, ".csv")),
        row.names = FALSE
      )
    }
  }
}

# 4.4 Generate forest plots
if (verbose_output) cat("Generating forest plots...\n")

for (tissue in names(univariate_cox_results)) {
  if (verbose_output) cat("  - Creating forest plot for", tissue, "...\n")

  forest_plot <- generate_forest_plot(
    cox_results = univariate_cox_results[[tissue]],
    tissue_name = tissue,
    config = config,
    save_path = file.path(survival_output_dir, paste0("forest_plot_", tissue, ".pdf"))
  )
}

# 4.5 Generate Kaplan-Meier curves for significant variables
if (verbose_output) cat("Generating Kaplan-Meier curves...\n")

for (tissue in names(univariate_cox_results)) {
  significant_vars <- univariate_cox_results[[tissue]][
    univariate_cox_results[[tissue]]$P_value <= alpha_level, "Variable"
  ]

  if (length(significant_vars) > 0) {
    for (var in significant_vars) {
      if (verbose_output) cat("  - Creating KM curve for", var, "in", tissue, "...\n")

      tissue_survival_data <- survival_data[survival_data[[tissue_col]] == tissue, ]

      analysis_data <- prepare_kaplan_meier_plot(
        survival_df = tissue_survival_data,
        variable = var,
        config = config,
        cutoff_method = "optimal"
      )
      
      # Fit survival model with error checking
      tryCatch({
        surv_object <- Surv(time = analysis_data$time, event = analysis_data$event)
        fit <- survfit(surv_object ~ Group, data = analysis_data)
      }, error = function(e) {
        stop("Error fitting survival model for ", var, ": ", e$message)
      })
      
      # Create clean variable name for plot
      clean_var_name <- clean_variable_names(var)
      tissue_name <- if ("Tissue" %in% colnames(survival_df)) unique(survival_df$Tissue)[1] else "Unknown"
      
      # Check the number of category
      category_ <- as.data.frame(table(analysis_data$Group))
      if(category_$Freq[1] * category_$Freq[2] == 0){
        cat(var, "in", tissue, "only has one category", "\n")
        next
      }
      
      # Generate Kaplan-Meier plot
      km_plot <- ggsurvplot(
        fit = fit,
        data = analysis_data,
        pval = TRUE,                    # Add p-value
        conf.int = TRUE,                # Add confidence intervals
        risk.table = TRUE,              # Add risk table
        risk.table.col = "strata",      # Color risk table by groups
        linetype = "strata",            # Different line types for groups
        surv.median.line = "hv",        # Add median survival lines
        ggtheme = theme_bw(),           # Clean theme
        palette = c("#bc3c29ff", "#0072b5ff"), # Custom colors
        title = "Recurrence-Free Survival by Cellular compostion Group",
        xlab = "Time (months)",
        ylab = "Survival probability",
        legend.title = "Group",
        legend.labs = c("High", "Low")
      )
      
      # Save plot if path provided
      km_output_dir <- file.path(survival_output_dir,"km")
      if(!file.exists(km_output_dir)){dir.create(km_output_dir,recursive = T)}
      pdf(file = file.path(km_output_dir, paste0("km_curve_", tissue, "_", var, ".pdf")), width = 8, height = 6)
      print(km_plot)
      dev.off()
      
      if (config$output$verbose) {
        cat("✓ Kaplan-Meier plot saved:", basename(km_output_dir), "\n")
        cat("  - Variable:", variable, "cutoff:", round(cutoff, 3), "\n")
        cat("  - Group sizes: Low =", group_sizes["Low"], ", High =", group_sizes["High"], "\n")
      }
    }
  }
}

# =============================================================================
# PART 5: SUMMARY AND REPORT GENERATION
# =============================================================================

# 5.2 Create combined plots summary (optional)
if (verbose_output) cat("Creating combined summary plots...\n")

# This section can include code to combine key plots into summary figures

# =============================================================================
# COMPLETION MESSAGE
# =============================================================================

cat("\n", rep("=", 70), "\n")
cat("CELL ABUNDANCE ANALYSIS SCRIPT LOADED SUCCESSFULLY\n")
cat(rep("=", 70), "\n")

cat("\nTo run the analysis:\n")
cat("1. Modify the configuration section at the top of this script\n")
cat("2. Load your SpatialExperiment object\n") 
cat("3. Uncomment the analysis sections you want to run\n")
cat("4. Execute the script line by line or in sections\n\n")

cat("Configuration summary:\n")
cat("- Tissues to analyze:", paste(tissues_to_analyze, collapse = ", "), "\n")
cat("- Output directory:", output_base_dir, "\n")
cat("- Statistical test:", test_method, "\n")
cat("- P-value adjustment:", p_adjust_method, "\n\n")

cat("Analysis sections available:\n")
cat("- Part 1: Data Description (pie charts, bar plots, heatmaps)\n")
cat("- Part 2: Discrete Clinical Analysis (group comparisons, bubble plots)\n") 
cat("- Part 3: Continuous Clinical Analysis (correlations, scatter plots)\n")
cat("- Part 4: Survival Analysis (Cox regression, forest plots, KM curves)\n")
cat("- Part 5: Summary and Report Generation\n\n")

cat(rep("=", 70), "\n")
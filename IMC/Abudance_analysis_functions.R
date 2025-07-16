# =============================================================================
# ABUNDANCE_FUNCTIONS.R - Modular Cell Abundance Analysis Functions
# =============================================================================

# Required libraries
required_packages <- c(
  "SpatialExperiment", "SingleCellExperiment", "dplyr", "tidyr",
  "ggplot2", "ggpubr", "pheatmap", "ComplexHeatmap", "circlize","igraph",
  "ggsci", "RColorBrewer", "patchwork", "survival", "survminer","ggcorrplot",
  "forestplot", "ggrepel", "stringr", "scales","tibble","tidyverse","Hmisc"
)  

# Load required packages
load_required_packages <- function(packages = required_packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
    library(pkg, character.only = TRUE)
  }
}

# =============================================================================
# 1. CONFIGURATION AND VALIDATION FUNCTIONS
# =============================================================================

#' Create Analysis Configuration
#' @param discrete_vars Character vector of discrete clinical variables
#' @param continuous_vars Character vector of continuous clinical variables
#' @param survival_vars Named list with 'time' and 'event' elements
#' @param tissues Character vector of tissue types to analyze
#' @param celltype_col Column name for cell type information
#' @param tissue_col Column name for tissue information
#' @param patient_col Column name for patient ID
#' @param roi_col Column name for ROI/ID information
#' @param plot_config List of plotting parameters
#' @param stats_config List of statistical parameters
create_analysis_config <- function(
  discrete_vars = c("RFS_status", "Gender", "KRAS_mutation"),
  continuous_vars = c("RFS_time", "Age", "TBS"),
  survival_vars = list(time = "RFS_time", event = "RFS_status"),
  tissues = c("CT", "IM", "TAT"),
  celltype_col = "SubType",
  tissue_col = "Tissue", 
  patient_col = "PID",
  roi_col = "ID",
  plot_config = list(
    width = 10, height = 8, dpi = 300, format = "pdf",
    color_palette = "Set2", point_size = 1.5, text_size = 12,
    title_size = 16, axis_text_angle = 45
  ),
  stats_config = list(
    p_adjust_method = "BH", correlation_method = "spearman",
    test_method = "wilcox.test", alpha_level = 0.05,
    fc_threshold = 1.2, min_cells_per_group = 5
  )
) {
  
  config <- list(
    # Data columns
    columns = list(
      discrete_clinical = discrete_vars,
      continuous_clinical = continuous_vars,
      survival = survival_vars,
      celltype = celltype_col,
      tissue = tissue_col,
      patient = patient_col,
      roi = roi_col
    ),
    
    # Analysis scope
    scope = list(
      tissues = tissues,
      include_unknown = FALSE,
      min_fraction_threshold = 0.001
    ),
    
    # Plotting parameters
    plots = plot_config,
    
    # Statistical parameters  
    stats = stats_config,
    
    # Output settings
    output = list(
      save_intermediate = TRUE,
      verbose = TRUE,
      create_report = TRUE
    )
  )
  
  # Validate configuration
  validate_config(config)
  return(config)
}

#' Validate Analysis Configuration
validate_config <- function(config) {
  required_fields <- c("columns", "scope", "plots", "stats", "output")
  missing_fields <- setdiff(required_fields, names(config))
  if (length(missing_fields) > 0) {
    stop("Missing required config fields: ", paste(missing_fields, collapse = ", "))
  }
  
  # Validate survival configuration
  if (!all(c("time", "event") %in% names(config$columns$survival))) {
    stop("Survival config must include 'time' and 'event' elements")
  }
}

#' Validate SpatialExperiment Input
#' @param spe SpatialExperiment object
#' @param config Analysis configuration
validate_spe_input <- function(spe, config) {
  
  # Check object type
  if (!is(spe, "SpatialExperiment") && !is(spe, "SingleCellExperiment")) {
    stop("Input must be a SpatialExperiment or SingleCellExperiment object")
  }
  
  # Get clinical table
  meta_table <- as.data.frame(colData(spe))
  
  # Check required columns
  coldata_cols <- colnames(meta_table)
  required_cols <- c(
    config$celltype, config$tissue,
    config$patient, config$roi,
    config$discrete_clinical, config$continuous_clinical,
    unlist(config$survival)
  )
  
  missing_cols <- setdiff(required_cols, coldata_cols)
  if (length(missing_cols) > 0) {
    stop("Missing required columns in colData: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check tissue types
  available_tissues <- unique(meta_table[[config$tissue]])
  missing_tissues <- setdiff(config$tissues, available_tissues)
  if (length(missing_tissues) > 0) {
    warning("Some specified tissues not found in data: ", paste(missing_tissues, collapse = ", "))
  }
  
  if (config$verbose) {
    cat("✓ Input validation passed\n")
    cat("  - Cells:", ncol(spe), "\n")
    cat("  - Tissues:", paste(available_tissues, collapse = ", "), "\n")
    cat("  - Cell types:", length(unique(meta_table[[config$celltype]])), "\n")
  }
}

# =============================================================================
# 2. DATA PROCESSING FUNCTIONS
# =============================================================================

#' Process SpatialExperiment for Abundance Analysis
#' @param spe SpatialExperiment object
#' @param config Analysis configuration
#' @param aggregate_level Character: "roi" or "patient" for aggregation level
process_spe_data <- function(spe, config, aggregate_level = "roi") {
  
  # Extract and filter metadata
  meta <- colData(spe) %>% as.data.frame()
  
  # Filter by specified tissues
  meta <- meta %>% 
    filter(!!sym(config$columns$tissue) %in% config$scope$tissues)
  
  # Remove unknown cell types if specified
  if (!config$scope$include_unknown) {
    meta <- meta %>% 
      filter(!grepl("unknown|Unknown|UNKNOWN|unlabelled|Unlabelled", !!sym(config$columns$celltype), ignore.case = TRUE))
  }
  
  # Calculate cell fractions
  if (aggregate_level == "roi") {
    fraction_df <- calculate_roi_fractions(meta, config)
  } else if (aggregate_level == "patient") {
    fraction_df <- calculate_patient_fractions(meta, config)
  } else {
    stop("aggregate_level must be 'roi' or 'patient'")
  }
  
  # Add clinical information
  clinical_df <- extract_clinical_data(meta, config)
  analysis_df <- merge_fractions_clinical(fraction_df, clinical_df, config)
  
  # Apply minimum fraction threshold
  if (config$scope$min_fraction_threshold > 0) {
    analysis_df <- apply_fraction_threshold(analysis_df, config)
  }
  
  if (config$output$verbose) {
    cat("✓ Data processing completed\n")
    cat("  - ROIs/Patients:", nrow(analysis_df), "\n")
    cat("  - Cell types:", sum(grepl("^[A-Z]", colnames(analysis_df))), "\n")
  }
  
  return(analysis_df)
}

#' Calculate Cell Type Fractions per ROI
#' @param meta Metadata data frame
#' @param config Analysis configuration
calculate_roi_fractions <- function(meta, config) {
  
  roi_fractions <- meta %>%
    group_by(!!sym(config$columns$roi), !!sym(config$columns$tissue)) %>%
    count(!!sym(config$columns$celltype), name = "count") %>%
    mutate(fraction = count / sum(count)) %>%
    select(-count) %>%
    pivot_wider(
      names_from = !!sym(config$columns$celltype),
      values_from = fraction,
      values_fill = 0
    ) %>%
    ungroup()
  
  return(roi_fractions)
}

#' Calculate Cell Type Fractions per Patient
#' @param meta Metadata data frame  
#' @param config Analysis configuration
calculate_patient_fractions <- function(meta, config) {
  
  patient_fractions <- meta %>%
    group_by(!!sym(config$columns$patient), !!sym(config$columns$tissue)) %>%
    count(!!sym(config$columns$celltype), name = "count") %>%
    mutate(fraction = count / sum(count)) %>%
    select(-count) %>%
    pivot_wider(
      names_from = !!sym(config$columns$celltype),
      values_from = fraction,
      values_fill = 0
    ) %>%
    ungroup()
  
  return(patient_fractions)
}

#' Extract Clinical Data
#' @param meta Metadata data frame
#' @param config Analysis configuration
extract_clinical_data <- function(meta, config) {
  
  clinical_cols <- c(
    config$columns$patient, config$columns$roi, config$columns$tissue,
    config$columns$discrete_clinical, config$columns$continuous_clinical,
    unname(unlist(config$columns$survival)) 
  )
  
  clinical_df <- meta %>%
    select(all_of(clinical_cols)) %>%
    distinct()
  
  return(clinical_df)
}

#' Merge Fractions with Clinical Data
#' @param fraction_df Data frame with cell fractions
#' @param clinical_df Data frame with clinical variables
#' @param config Analysis configuration
merge_fractions_clinical <- function(fraction_df, clinical_df, config) {
  
  # Determine merge key based on data structure
  if (config$columns$roi %in% colnames(fraction_df)) {
    merge_key <- config$columns$roi
  } else {
    merge_key <- config$columns$patient
  }
  
  analysis_df <- fraction_df %>%
    left_join(clinical_df, by = c(merge_key, config$columns$tissue))
  
  return(analysis_df)
}

# =============================================================================
# 3. DATA DESCRIPTION FUNCTIONS
# =============================================================================

#' Generate Pie Charts for Cell Type Composition
#' @param fractions_df Data frame with cell fractions
#' @param config Analysis configuration
#' @param save_path Path to save the plot
#' @param show_percentages Logical: show percentage labels
#' @param min_percentage_label Minimum percentage to show label
plot_celltype_pie_charts <- function(fractions_df, config, save_path = NULL,
                                   show_percentages = TRUE, min_percentage_label = 2,
                                   width = 8, height = 6) {
  
  tissues <- intersect(config$scope$tissues, unique(fractions_df[[config$columns$tissue]]))
  plot_list <- list()
  
  # Get cell type columns
  celltype_cols <- get_celltype_columns(fractions_df, config)
  
  for (tissue in tissues) {
    tissue_data <- fractions_df %>%
      filter(!!sym(config$columns$tissue) == tissue) %>%
      select(all_of(celltype_cols)) %>%
      summarise_all(mean, na.rm = TRUE) %>%
      pivot_longer(everything(), names_to = "CellType", values_to = "Fraction") %>%
      filter(Fraction > 0) %>%
      arrange(desc(Fraction))
    
    p <- ggplot(tissue_data, aes(x = "", y = Fraction, fill = CellType)) +
      geom_bar(width = 1, stat = "identity", color = "black", size = 0.5) +
      coord_polar("y", start = 0) +
      labs(title = tissue) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = config$plots$title_size, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = config$plots$text_size - 2)
      )
    
    # Add percentage labels
    if (show_percentages) {
      tissue_data <- tissue_data %>%
        mutate(
          percentage = Fraction * 100,
          label = ifelse(percentage >= min_percentage_label,
                        paste0(round(percentage, 1), "%"), "")
        )
      
      p <- p + geom_text(
        data = tissue_data,
        aes(label = label),
        position = position_stack(vjust = 0.5),
        color = "white", size = 3, fontface = "bold"
      )
    }
    
    # Apply color palette
    if (!is.null(config$plots$color_palette)) {
      palette_ <- config$plots$color_palette[[config$columns$celltype]]
      p <- p + scale_fill_manual(values = palette_)
    }
    
    plot_list[[tissue]] <- p
  }
  
  combined_plot <- wrap_plots(plot_list, nrow = 1)
  
  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, width = width, height = height, dpi = config$plots$dpi)}
  
  
  return(combined_plot)
}

#' Generate Stacked Bar Plot for Patient Composition
#' @param fractions_df Data frame with cell fractions
#' @param config Analysis configuration
#' @param save_path Path to save the plot
#' @param order_by Method to order patients ("tissue", "similarity", or column name)
#' @param show_tissue_separator Logical: show tissue separators
plot_stacked_abundance_bars <- function(fractions_df, config, save_path = NULL,
                                      order_by = "tissue", show_tissue_separator = TRUE,
                                      width = 8, height = 6) {
  
  # Get cell type columns
  celltype_cols <- get_celltype_columns(fractions_df, config)
  
  # Prepare data for plotting
  plot_data <- fractions_df %>%
    select(all_of(c(config$columns$roi, config$columns$tissue, celltype_cols))) %>%
    pivot_longer(
      cols = all_of(celltype_cols),
      names_to = "CellType", 
      values_to = "Fraction"
    )
  
  # Order ROIs/patients
  if (order_by == "tissue") {
    roi_order <- fractions_df %>%
      arrange(!!sym(config$columns$tissue), !!sym(config$columns$roi)) %>%
      pull(!!sym(config$columns$roi))
  } else if (order_by == "similarity") {
    # Order by hierarchical clustering
    roi_order <- order_rois_by_similarity(fractions_df, celltype_cols)
  } else if (order_by %in% colnames(fractions_df)) {
    roi_order <- fractions_df %>%
      arrange(!!sym(config$columns$tissue), !!sym(order_by)) %>%
      pull(!!sym(config$columns$roi))
  } else {
    roi_order <- unique(plot_data[[config$columns$roi]])
  }
  
  plot_data[[config$columns$roi]] <- factor(plot_data[[config$columns$roi]], levels = roi_order)
  
  # Create main plot
  p <- ggplot(plot_data, aes(x = !!sym(config$columns$roi), y = Fraction, fill = CellType)) +
    geom_bar(stat = "identity", color = "black", size = 0.1) +
    labs(
      x = "", y = "Cell Fraction",
      title = "Cell Type Composition Across Samples"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = config$plots$title_size, face = "bold"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = config$plots$text_size),
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid = element_blank()
    )
  
  # Apply color palette
  if (!is.null(config$plots$color_palette)) {
    palette_ <- config$plots$color_palette[[config$columns$celltype]]
    p <- p + scale_fill_manual(values = palette_)
    }
  
  # Add tissue separator if requested
  if (show_tissue_separator) {
    tissue_bar <- create_tissue_separator(fractions_df, config, roi_order)
    combined_plot <- p / tissue_bar + plot_layout(heights = c(10, 1))
  } else {
    combined_plot <- p
  }
  
  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot,dpi = config$plots$dpi, height = height, width = width)
  }
  
  return(combined_plot)
}

# =============================================================================
# 4. DISCRETE CLINICAL ANALYSIS FUNCTIONS
# =============================================================================

#' Perform Group Comparisons for Discrete Variables
#' @param fractions_df Data frame with cell fractions
#' @param config Analysis configuration
#' @param test_method Statistical test method
#' @param effect_size Logical: calculate effect sizes
calculate_group_comparisons <- function(fractions_df, config, 
                                      test_method = NULL, effect_size = TRUE) {
  
  if (is.null(test_method)) {
    test_method <- config$stats$test_method
  }
  
  results_list <- list()
  celltype_cols <- get_celltype_columns(fractions_df, config)
  
  for (var in config$columns$discrete_clinical) {
    for (tissue in config$scope$tissues) {
      
      tissue_data <- fractions_df %>% 
        filter(!!sym(config$columns$tissue) == tissue) %>%
        filter(!is.na(!!sym(var)))
      
      if (nrow(tissue_data) < config$stats$min_cells_per_group * 2) {
        if (config$output$verbose) {
          cat("Skipping", var, "in", tissue, "- insufficient data\n")
        }
        next
      }
      
      var_results <- perform_variable_comparison(
        tissue_data, var, celltype_cols, test_method, effect_size, config
      )
      
      if (!is.null(var_results)) {
        var_results$Variable <- var
        var_results$Tissue <- tissue
        results_list[[paste(var, tissue, sep = "_")]] <- var_results
      }
    }
  }
  
  if (length(results_list) == 0) {
    warning("No valid comparisons could be performed")
    return(NULL)
  }
  
  # Combine all results
  all_results <- bind_rows(results_list)
  
  # Apply multiple testing correction
  all_results$p_adjusted <- p.adjust(all_results$p_value, method = config$stats$p_adjust_method)
  
  # Add significance categories
  all_results <- all_results %>%
    mutate(
      significance = case_when(
        p_adjusted < 0.001 ~ "***",
        p_adjusted < 0.01 ~ "**",
        p_adjusted < 0.05 ~ "*", 
        TRUE ~ "ns"
      ),
      is_significant = p_adjusted < config$stats$alpha_level,
      passes_fc_threshold = abs(log2(fold_change)) > log2(config$stats$fc_threshold)
    )
  
  if (config$output$verbose) {
    n_sig <- sum(all_results$is_significant, na.rm = TRUE)
    cat("✓ Group comparisons completed\n")
    cat("  - Total tests:", nrow(all_results), "\n")
    cat("  - Significant results:", n_sig, "\n")
  }
  
  return(all_results)
}

#' Perform Statistical Comparison for Single Variable
#' @param data Data frame with fractions and grouping variable
#' @param var_name Name of grouping variable
#' @param celltype_cols Names of cell type columns
#' @param test_method Statistical test method
#' @param effect_size Logical: calculate effect sizes
#' @param config Analysis configuration
perform_variable_comparison <- function(data, var_name, celltype_cols, 
                                      test_method, effect_size, config) {
  
  groups <- unique(data[[var_name]])
  groups <- groups[!is.na(groups)]
  
  if (length(groups) < 2) {
    return(NULL)
  }
  
  # For now, handle binary comparisons
  if (length(groups) == 2) {
    groups <- sort(groups)
    group1_data <- data[data[[var_name]] == groups[1], ]
    group2_data <- data[data[[var_name]] == groups[2], ]
    
    # Check minimum group sizes
    if (nrow(group1_data) < config$stats$min_cells_per_group || 
        nrow(group2_data) < config$stats$min_cells_per_group) {
      return(NULL)
    }
    
    results <- data.frame(
      CellType = celltype_cols,
      fold_change = numeric(length(celltype_cols)),
      p_value = numeric(length(celltype_cols)),
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(celltype_cols)) {
      celltype <- celltype_cols[i]
      
      v1 <- group1_data[[celltype]]
      v2 <- group2_data[[celltype]]
      
      # Remove NA values
      v1 <- v1[!is.na(v1)]
      v2 <- v2[!is.na(v2)]
      
      if (length(v1) < 3 || length(v2) < 3) {
        results$fold_change[i] <- NA
        results$p_value[i] <- NA
        next
      }
      
      # Calculate fold change (group2 vs group1)
      mean1 <- mean(v1)
      mean2 <- mean(v2)
      results$fold_change[i] <- ifelse(mean1 > 0, mean2 / mean1, NA)
      
      # Perform statistical test
      if (test_method == "wilcox.test") {
        test_result <- wilcox.test(v2, v1)
      } else if (test_method == "t.test") {
        test_result <- t.test(v2, v1)
      } else {
        stop("Unsupported test method: ", test_method)
      }
      
      results$p_value[i] <- test_result$p.value
      
      # Add effect size if requested
      if (effect_size) {
        if (test_method == "wilcox.test") {
          # Calculate rank-biserial correlation
          results$effect_size[i] <- calculate_rank_biserial(v1, v2)
        } else {
          # Calculate Cohen's d
          results$effect_size[i] <- calculate_cohens_d(v1, v2)
        }
      }
    }
    
    return(results)
  }
  
  return(NULL)
}

#' Create Bubble Plot for Group Comparisons
#' @param stats_df Data frame with statistical results
#' @param config Analysis configuration
#' @param save_path Path to save the plot
#' @param fc_limits Limits for fold change color scale
#' @param show_ns Logical: show non-significant results
plot_abundance_bubble_chart <- function(stats_df, config, save_path = NULL,
                                       fc_limits = c(0.5, 2), show_ns = FALSE, width, height) {
  
  if (is.null(stats_df) || nrow(stats_df) == 0) {
    warning("No data available for bubble chart")
    return(NULL)
  }
  
  # Filter data if requested
  if (!show_ns) {
    stats_df <- stats_df %>% filter(is_significant)
  }
  
  if (nrow(stats_df) == 0) {
    warning("No significant results to plot")
    return(NULL)
  }
  
  # Prepare plotting data
  plot_data <- stats_df %>%
    mutate(
      log2_fc = log2(fold_change),
      log2_fc_capped = pmax(pmin(log2_fc, log2(fc_limits[2])), log2(fc_limits[1])),
      size_category = case_when(
        p_adjusted <= 0.001 ~ 4,
        p_adjusted <= 0.01 ~ 3,
        p_adjusted <= 0.05 ~ 2,
        TRUE ~ 1
      ),
      Variable_clean = clean_variable_names(Variable)
    ) %>%
    filter(!is.na(log2_fc))
  
  p <- ggplot(plot_data, aes(x = CellType, y = Variable_clean,
                            size = size_category, color = log2_fc_capped)) +
    geom_point(alpha = 0.8) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, name = "log2(FC)",
      limits = log2(fc_limits)
    ) +
    scale_size_continuous(
      range = c(2, 8), name = "Significance",
      breaks = c(1, 2, 3, 4),
      labels = c("ns", "*", "**", "***"),
      guide = guide_legend(override.aes = list(color = "black"))
    ) +
    facet_wrap(~Tissue, scales = "free_x", ncol = length(config$scope$tissues)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                size = config$plots$text_size - 2),
      axis.text.y = element_text(size = config$plots$text_size),
      axis.title = element_text(size = config$plots$text_size + 2, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = config$plots$title_size, face = "bold"),
      strip.text = element_text(size = config$plots$text_size + 1, face = "bold"),
      legend.position = "bottom"
    ) +
    labs(
      title = "Cell Type Abundance Differences Across Clinical Groups",
      x = "Cell Type", y = "Clinical Variable"
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = width, height = height,dpi = config$plots$dpi)
  }
  
  return(p)
}

#' Generate Individual Boxplots for Clinical Variables
#' @param fractions_df Data frame with cell fractions and clinical data
#' @param clinical_var Name of clinical variable to compare
#' @param config Analysis configuration
#' @param save_path Base path for saving plots (will create individual files)
#' @param test_method Statistical test method ("t.test" or "wilcox.test")
#' @param show_significance Logical: show statistical significance
#' @param width Plot width
#' @param height Plot height
#' @param show_points Logical: show individual data points
#' @param point_alpha Alpha value for points
#' @param box_alpha Alpha value for boxes
plot_abundance_boxplots <- function(fractions_df, clinical_var, config, 
                                    save_path = NULL, test_method = "t.test",
                                    show_significance = TRUE, width = 6, 
                                    height = 10, show_points = TRUE,
                                    point_alpha = 0.6, box_alpha = 0.7) {
  
  # Validate inputs
  if (!clinical_var %in% colnames(fractions_df)) {
    stop("Clinical variable '", clinical_var, "' not found in data")
  }
  
  # Get cell type columns
  celltype_cols <- get_celltype_columns(fractions_df, config)
  
  if (length(celltype_cols) == 0) {
    stop("No cell type columns found in data")
  }
  
  # Get tissue column name
  tissue_col <- config$columns$tissue
  if (!"Tissue" %in% colnames(fractions_df) && tissue_col %in% colnames(fractions_df)) {
    # Rename for consistency with reference code
    fractions_df$Tissue <- fractions_df[[tissue_col]]
  }
  
  # Check if we have Tissue column
  if (!"Tissue" %in% colnames(fractions_df)) {
    stop("Tissue column not found in data")
  }
  
  # Clean clinical variable values and convert to factor
  fractions_df[[clinical_var]] <- as.factor(fractions_df[[clinical_var]])
  
  # Remove rows with missing clinical variable
  fractions_df <- fractions_df[!is.na(fractions_df[[clinical_var]]), ]
  
  if (nrow(fractions_df) == 0) {
    stop("No data remaining after removing missing values for ", clinical_var)
  }
  
  # Check number of groups
  n_groups <- length(unique(fractions_df[[clinical_var]]))
  if (n_groups < 2) {
    warning("Less than 2 groups found for ", clinical_var, ". Skipping.")
    return(NULL)
  }
  
  # Set up output directory if save_path provided
  if (!is.null(save_path)) {
    output_dir <- dirname(save_path)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Create subdirectory for boxplots
    boxplot_dir <- file.path(output_dir, "individual_boxplots")
    if (!dir.exists(boxplot_dir)) {
      dir.create(boxplot_dir, recursive = TRUE)
    }
  }
  
  # Create clean variable name for titles
  clean_var_name <- clean_variable_names(clinical_var)
  
  # Get appropriate colors
  if (n_groups == 2) {
    fill_colors <- pal_jco("default")(2)
  } else {
    fill_colors <- pal_jco("default")(min(n_groups, 10))
  }
  
  # Store plots for return
  plot_list <- list()
  
  # Loop through each cell type
  for (celltype in celltype_cols) {
    
    if (config$output$verbose) {
      cat("  - Creating boxplot for", celltype, "...\n")
    }
    
    # Prepare data for this cell type
    plot_data <- fractions_df[, c(celltype, "Tissue", clinical_var)]
    colnames(plot_data) <- c("Fraction", "Tissue", "Clinical_Group")
    
    # Convert to appropriate types
    plot_data$Fraction <- as.numeric(plot_data$Fraction)
    plot_data$Clinical_Group <- as.factor(plot_data$Clinical_Group)
    plot_data$Tissue <- as.factor(plot_data$Tissue)
    
    # Remove rows with missing fractions
    plot_data <- plot_data[!is.na(plot_data$Fraction), ]
    
    if (nrow(plot_data) == 0) {
      warning("No data for cell type ", celltype, ". Skipping.")
      next
    }
    
    # Create the base plot
    p <- ggplot(plot_data, aes(x = Clinical_Group, y = Fraction, fill = Clinical_Group)) +
      geom_boxplot(alpha = box_alpha, color = "black", outlier.shape = NA) +
      scale_fill_manual(values = fill_colors) +
      theme_bw() +
      labs(
        x = clean_var_name,
        y = "Cell Fraction", 
        title = paste0("Abundance of ", celltype, " by ", clean_var_name)
      ) +
      theme(
        plot.title = element_text(size = config$plots$title_size, face = "bold", hjust = 0.5),
        text = element_text(size = config$plots$text_size),
        axis.title = element_text(face = "bold", size = config$plots$text_size + 2),
        axis.text.x = element_text(size = config$plots$text_size - 1, vjust = 0.5),
        axis.text.y = element_text(size = config$plots$text_size - 1),
        strip.background = element_blank(),
        strip.text = element_text(size = config$plots$text_size, face = "bold"),
        legend.position = "top"
      )
    
    # Add individual points if requested
    if (show_points) {
      p <- p + geom_jitter(
        color = "darkgrey", 
        position = position_jitter(width = 0.2), 
        size = 1, 
        alpha = point_alpha
      )
    }
    
    # Add statistical testing if requested
    if (show_significance && n_groups == 2) {
      p <- p + stat_compare_means(
        aes(group = Clinical_Group),
        method = test_method,
        hide.ns = FALSE,
        label = "p",
        label.y.npc = 0.75,
        size = 6
      )
    } else if (show_significance && n_groups > 2) {
      # For more than 2 groups, use ANOVA or Kruskal-Wallis
      if (test_method == "t.test") {
        p <- p + stat_compare_means(method = "anova", label.y.npc = "0.75")
      } else {
        p <- p + stat_compare_means(method = "kruskal.test", label.y.npc = "0.75")
      }
    }
    
    # Add faceting by tissue
    p <- p + facet_grid(Tissue ~ ., scales = "free_y")
    
    # Adjust plot margins
    p <- p + theme(plot.margin = margin(10, 10, 10, 10))
    
    # Store plot
    plot_list[[celltype]] <- p
    
    # Save individual plot if path provided
    if (!is.null(save_path)) {
      
      # Create filename
      safe_celltype <- gsub("[^A-Za-z0-9_]", "_", celltype)  # Make filename safe
      filename <- paste0("boxplot_", safe_celltype, "_by_", clinical_var, ".pdf")
      full_path <- file.path(boxplot_dir, filename)
      
      # Save plot
      ggsave(
        filename = full_path,
        plot = p,
        width = width,
        height = height,
        dpi = config$plots$dpi
      )
      
      if (config$output$verbose) {
        cat("    ✓ Saved:", filename, "\n")
      }
    }
  }
  
  if (config$output$verbose) {
    cat("✓ Generated", length(plot_list), "boxplots for", clinical_var, "\n")
  }
  
  return(NULL)
}

# =============================================================================
# 4. DISCRETE CLINICAL ANALYSIS FUNCTIONS
# =============================================================================

#' Perform Statistical Comparison for Single Variable
#' @param data Data frame with fractions and grouping variable
#' @param var_name Name of grouping variable
#' @param celltype_cols Names of cell type columns
#' @param test_method Statistical test method
#' @param effect_size Logical: calculate effect sizes
#' @param config Analysis configuration
perform_variable_comparison <- function(data, var_name, celltype_cols, 
                                        test_method, effect_size, config) {
  
  groups <- unique(data[[var_name]])
  groups <- groups[!is.na(groups)]
  
  if (length(groups) < 2) {
    return(NULL)
  }
  
  # For now, handle binary comparisons
  if (length(groups) == 2) {
    groups <- sort(groups)
    group1_data <- data[data[[var_name]] == groups[1], ]
    group2_data <- data[data[[var_name]] == groups[2], ]
    
    # Check minimum group sizes
    if (nrow(group1_data) < config$stats$min_cells_per_group || 
        nrow(group2_data) < config$stats$min_cells_per_group) {
      return(NULL)
    }
    
    results <- data.frame(
      CellType = celltype_cols,
      fold_change = numeric(length(celltype_cols)),
      p_value = numeric(length(celltype_cols)),
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(celltype_cols)) {
      celltype <- celltype_cols[i]
      
      v1 <- group1_data[[celltype]]
      v2 <- group2_data[[celltype]]
      
      # Remove NA values
      v1 <- v1[!is.na(v1)]
      v2 <- v2[!is.na(v2)]
      
      if (length(v1) < 3 || length(v2) < 3) {
        results$fold_change[i] <- NA
        results$p_value[i] <- NA
        next
      }
      
      # Calculate fold change (group2 vs group1)
      mean1 <- mean(v1)
      mean2 <- mean(v2)
      results$fold_change[i] <- ifelse(mean1 > 0, mean2 / mean1, NA)
      
      # Perform statistical test
      if (test_method == "wilcox.test") {
        test_result <- wilcox.test(v2, v1)
      } else if (test_method == "t.test") {
        test_result <- t.test(v2, v1)
      } else {
        stop("Unsupported test method: ", test_method)
      }
      
      results$p_value[i] <- test_result$p.value
      
      # Add effect size if requested
      if (effect_size) {
        if (test_method == "wilcox.test") {
          # Calculate rank-biserial correlation
          results$effect_size[i] <- calculate_rank_biserial(v1, v2)
        } else {
          # Calculate Cohen's d
          results$effect_size[i] <- calculate_cohens_d(v1, v2)
        }
      }
    }
    
    return(results)
  }
  
  return(NULL)
}

# =============================================================================
# 5. CONTINUOUS CLINICAL ANALYSIS FUNCTIONS
# =============================================================================

#' Calculate Correlations with Continuous Variables
#' @param fractions_df Data frame with cell fractions
#' @param config Analysis configuration
#' @param method Correlation method
#' @param min_observations Minimum observations for correlation
calculate_abundance_correlations <- function(fractions_df, config, 
                                             method = NULL, min_observations = 10) {
  
  if (is.null(method)) {
    method <- config$stats$correlation_method
  }
  
  correlation_results <- list()
  celltype_cols <- get_celltype_columns(fractions_df, config)
  
  for (tissue in config$scope$tissues) {
    tissue_data <- fractions_df %>% 
      filter(!!sym(config$columns$tissue) == tissue)
    
    if (nrow(tissue_data) < min_observations) {
      if (config$output$verbose) {
        cat("Skipping", tissue, "- insufficient observations for correlation\n")
      }
      next
    }
    
    # Calculate correlations
    cor_results <- calculate_tissue_correlations(
      tissue_data, celltype_cols, config$columns$continuous_clinical, 
      method, min_observations
    )
    
    if (!is.null(cor_results)) {
      correlation_results[[tissue]] <- cor_results
    }
  }
  
  if (config$output$verbose && length(correlation_results) > 0) {
    cat("✓ Correlation analysis completed\n")
    cat("  - Tissues analyzed:", paste(names(correlation_results), collapse = ", "), "\n")
  }
  
  return(correlation_results)
}

#' Plot Correlation Heatmap with Significance
#'
#' This function creates a heatmap of correlation coefficients and displays
#' significance levels as asterisks.
#'
#' @param correlation_results A list containing 'cor_matrix' and 'p_matrix'.
#' @param tissue_name The name of the tissue for the plot title.
#' @param config The analysis configuration object.
#' @param save_path The full file path to save the generated PDF.
#'
plot_correlation_heatmap <- function(correlation_results, tissue_name, config, save_path, width, height) {
  # Ensure pheatmap is available for plotting
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Please install it with install.packages('pheatmap')", call. = FALSE)
  }
  
  # Extract the correlation and p-value matrices from the input list
  cor_matrix <- correlation_results$cor_matrix
  p_matrix <- correlation_results$p_matrix
  
  # Exit if the correlation matrix is empty or invalid
  if (is.null(cor_matrix) || all(is.na(cor_matrix)) || nrow(cor_matrix) == 0 || ncol(cor_matrix) == 0) {
    if (config$output$verbose) {
      cat("Skipping heatmap for", tissue_name, "- no valid correlation data.\n")
    }
    return(invisible(NULL))
  }
  
  # Create a matrix of significance symbols based on p-values
  significance_matrix <- ifelse(p_matrix <= 0.001, "***",
                                ifelse(p_matrix <= 0.01, "**",
                                       ifelse(p_matrix <= 0.05, "*", "")))
  
  # Define a diverging color palette for correlations from -1 to 1
  color_palette <- colorRampPalette(c("#053061", "blue", "white", "red", "#B2182B"))(100)
  
  # Generate the heatmap and save it to the specified path
  pheatmap::pheatmap(
    mat               = cor_matrix,
    color             = color_palette,
    breaks            = seq(-1, 1, length.out = 101), # Ensure white is centered at 0
    display_numbers   = significance_matrix,
    fontsize_number   = 12,
    main              = paste("Correlations in", tissue_name),
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    filename          = save_path,
    width             = width, # Adjust width as needed for your plot
    height            = height  # Adjust height as needed for your plot
  )
  
  if (config$output$verbose) {
    cat("  - Correlation heatmap for", tissue_name, "saved to:", save_path, "\n")
  }
}

# =============================================================================
# 6. SURVIVAL ANALYSIS FUNCTIONS  
# =============================================================================

#' Prepare Data for Survival Analysis
#' @param fractions_df Data frame with cell fractions
#' @param clinical_df Data frame with clinical data (optional, can be NULL if already merged)
#' @param config Analysis configuration
#' @param aggregation_method Method for aggregating ROI data to patient level
prepare_survival_data <- function(fractions_df, clinical_df = NULL, config,
                                  aggregation_method = "mean") {
  
  # If clinical data is separate, merge it
  if (!is.null(clinical_df)) {
    fractions_df <- merge_fractions_clinical(fractions_df, clinical_df, config)
  }
  
  # Check if we need to aggregate by patient
  if (config$columns$roi %in% colnames(fractions_df)) {
    # Aggregate ROI-level data to patient level
    celltype_cols <- get_celltype_columns(fractions_df, config)
    
    if (aggregation_method == "mean") {
      agg_func <- mean
    } else if (aggregation_method == "median") {
      agg_func <- median
    } else {
      stop("Unsupported aggregation method: ", aggregation_method)
    }
    
    survival_df <- fractions_df %>%
      group_by(!!sym(config$columns$patient), !!sym(config$columns$tissue)) %>%
      summarise_at(vars(all_of(celltype_cols)), agg_func, na.rm = TRUE) %>%
      ungroup()
    
    # Add back clinical data (take first occurrence per patient)
    clinical_cols <- c(
      config$columns$discrete_clinical, config$columns$continuous_clinical,
      unname(unlist(config$columns$survival)) 
    )
    clinical_cols <- unique(clinical_cols)
    
    clinical_data <- fractions_df %>%
      select(all_of(c(config$columns$patient, clinical_cols))) %>%
      distinct()
    
    survival_df <- survival_df %>%
      left_join(clinical_data, by = config$columns$patient)
    
  } else {
    survival_df <- fractions_df
  }
  
  # Remove patients with missing survival data
  survival_cols <- unname(unlist(config$columns$survival))  
  survival_df <- survival_df %>%
    filter(!is.na(!!sym(survival_cols[1])) & !is.na(!!sym(survival_cols[2])))
  
  if (config$output$verbose) {
    cat("✓ Survival data prepared\n")
    cat("  - Patients:", length(unique(survival_df[[config$columns$patient]])), "\n")
    cat("  - Tissues:", paste(unique(survival_df[[config$columns$tissue]]), collapse = ", "), "\n")
  }
  
  return(survival_df)
}

#' Perform Cox Regression Analysis
#' @param survival_df Data frame prepared for survival analysis
#' @param config Analysis configuration
#' @param analysis_type Type of analysis: "univariate", "multivariate", or "both"
#' @param include_clinical Logical: include clinical variables in multivariate analysis
perform_cox_regression <- function(survival_df, config, analysis_type = "both",
                                   include_clinical = TRUE) {
  
  results <- list()
  celltype_cols <- get_celltype_columns(survival_df, config)
  survival_cols <- unlist(config$columns$survival)
  
  for (tissue in config$scope$tissues) {
    tissue_data <- survival_df %>% 
      filter(!!sym(config$columns$tissue) == tissue)
    
    if (nrow(tissue_data) < 20) {  # Minimum for Cox regression
      if (config$output$verbose) {
        cat("Skipping", tissue, "- insufficient patients for Cox regression\n")
      }
      next
    }
    
    tissue_results <- list()
    
    # Univariate analysis
    if (analysis_type %in% c("univariate", "both")) {
      uni_results <- perform_univariate_cox(tissue_data, celltype_cols, survival_cols, config)
      tissue_results$univariate <- uni_results
    }
    
    # Multivariate analysis
    if (analysis_type %in% c("multivariate", "both")) {
      multi_results <- perform_multivariate_cox(tissue_data, celltype_cols, 
                                                survival_cols, config, include_clinical)
      tissue_results$multivariate <- multi_results
    }
    
    results[[tissue]] <- tissue_results
  }
  
  if (config$output$verbose) {
    cat("✓ Cox regression analysis completed\n")
    cat("  - Analysis type:", analysis_type, "\n")
    cat("  - Tissues analyzed:", paste(names(results), collapse = ", "), "\n")
  }
  
  return(results)
}

# =============================================================================
# 7. UTILITY FUNCTIONS
# =============================================================================

#' Clean Variable Names for Display
#' @param variable_names Character vector of variable names
clean_variable_names <- function(variable_names) {
  cleaned <- variable_names
  
  # Define cleaning mappings
  name_mappings <- c(
    "RFS_status" = "Relapse vs Non-relapse",
    "Gender" = "Female vs Male", 
    "KRAS_mutation" = "KRAS Mut vs WT",
    "BRAF_mutation" = "BRAF Mut vs WT",
    "Differential_grade" = "Well vs Moderate Diff.",
    "Lymph_positive" = "LN+ vs LN-"
  )
  
  for (old_name in names(name_mappings)) {
    cleaned <- gsub(old_name, name_mappings[old_name], cleaned)
  }
  
  return(cleaned)
}

#' Calculate Effect Sizes
calculate_cohens_d <- function(x, y) {
  pooled_sd <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / 
                      (length(x) + length(y) - 2))
  (mean(y) - mean(x)) / pooled_sd
}

calculate_rank_biserial <- function(x, y) {
  # Simplified rank-biserial correlation
  n1 <- length(x)
  n2 <- length(y)
  U <- wilcox.test(x, y)$statistic
  (2 * U) / (n1 * n2) - 1
}

# =============================================================================
# 8. MISSING HELPER FUNCTIONS (NEED TO BE IMPLEMENTED)
# =============================================================================

#' Create Tissue Separator Bar for Stacked Plots
create_tissue_separator <- function(fractions_df, config, roi_order) {
  # Implementation needed for tissue separator in stacked bar plots
  tissue_data <- fractions_df %>%
    select(all_of(c(config$columns$roi, config$columns$tissue))) %>%
    distinct() %>%
    mutate(!!sym(config$columns$roi) := factor(!!sym(config$columns$roi), levels = roi_order))
  
  tissue_bar <- ggplot(tissue_data, aes(x = !!sym(config$columns$roi), y = 0.5, fill = !!sym(config$columns$tissue))) +
    geom_bar(stat = "identity", width = 1) +
    theme_void() +
    theme(legend.position = "bottom")
  
  return(tissue_bar)
}

#' Order ROIs by Similarity (Hierarchical Clustering)
order_rois_by_similarity <- function(fractions_df, celltype_cols) {
  # Implementation needed for ordering samples by similarity
  fraction_matrix <- fractions_df %>% 
    select(all_of(celltype_cols)) %>%
    as.matrix()
  
  if (nrow(fraction_matrix) > 1) {
    hc <- hclust(dist(fraction_matrix))
    return(fractions_df[[config$columns$roi]][hc$order])
  } else {
    return(fractions_df[[config$columns$roi]])
  }
}

#' Apply Fraction Threshold Filter
apply_fraction_threshold <- function(analysis_df, config) {
  # Implementation needed for filtering low-abundance cell types
  celltype_cols <- get_celltype_columns(analysis_df, config)
  
  for (col in celltype_cols) {
    analysis_df[[col]][analysis_df[[col]] < config$scope$min_fraction_threshold] <- 0
  }
  
  return(analysis_df)
}

#' Calculate Summary Statistics
calculate_summary_statistics <- function(analysis_data, config) {
  # Implementation needed for summary statistics
  celltype_cols <- get_celltype_columns(analysis_data, config)
  
  summary_stats <- list()
  
  for (tissue in config$scope$tissues) {
    tissue_data <- analysis_data %>% filter(!!sym(config$columns$tissue) == tissue)
    
    tissue_summary <- tissue_data %>%
      select(all_of(celltype_cols)) %>%
      summarise_all(list(mean = mean, median = median, sd = sd), na.rm = TRUE)
    
    summary_stats[[tissue]] <- tissue_summary
  }
  
  return(summary_stats)
}

#' Generate Individual Boxplots (wrapper function for backward compatibility)
#' @param analysis_data Data frame with cell fractions and clinical data
#' @param config Analysis configuration
#' @param output_dir Output directory for saving plots
generate_individual_boxplots <- function(analysis_data, config, output_dir) {
  
  if (config$output$verbose) {
    cat("Generating individual boxplots for discrete clinical variables...\n")
  }
  
  all_plots <- list()
  
  # Loop through each discrete clinical variable
  for (clinical_var in config$columns$discrete_clinical) {
    
    # Check if variable exists in data
    if (!clinical_var %in% colnames(analysis_data)) {
      if (config$output$verbose) {
        cat("  Warning: Variable", clinical_var, "not found in data. Skipping.\n")
      }
      next
    }
    
    if (config$output$verbose) {
      cat("  Processing", clinical_var, "...\n")
    }
    
    # Create boxplots for this variable
    var_plots <- plot_abundance_boxplots(
      fractions_df = analysis_data,
      clinical_var = clinical_var,
      config = config,
      save_path = file.path(output_dir, paste0("boxplots_", clinical_var, ".pdf")),
      test_method = config$stats$test_method,
      show_significance = TRUE,
      plot_width = 6,
      plot_height = 10
    )
    
    if (!is.null(var_plots)) {
      all_plots[[clinical_var]] <- var_plots
    }
  }
  
  return(all_plots)
}

# =============================================================================
# 5. CONTINUOUS CLINICAL ANALYSIS FUNCTIONS
# =============================================================================

#' Calculate Correlations with Continuous Variables
#' @param fractions_df Data frame with cell fractions
#' @param config Analysis configuration
#' @param method Correlation method
#' @param min_observations Minimum observations for correlation
calculate_abundance_correlations <- function(fractions_df, config, 
                                             method = NULL, min_observations = 10) {
  
  if (is.null(method)) {
    method <- config$stats$correlation_method
  }
  
  correlation_results <- list()
  celltype_cols <- get_celltype_columns(fractions_df, config)
  
  for (tissue in config$scope$tissues) {
    tissue_data <- fractions_df %>% 
      filter(!!sym(config$columns$tissue) == tissue)
    
    if (nrow(tissue_data) < min_observations) {
      if (config$output$verbose) {
        cat("Skipping", tissue, "- insufficient observations for correlation\n")
      }
      next
    }
    
    # Calculate correlations
    cor_results <- calculate_tissue_correlations(
      tissue_data, celltype_cols, config$columns$continuous_clinical, 
      method, min_observations
    )
    
    if (!is.null(cor_results)) {
      correlation_results[[tissue]] <- cor_results
    }
  }
  
  if (config$output$verbose && length(correlation_results) > 0) {
    cat("✓ Correlation analysis completed\n")
    cat("  - Tissues analyzed:", paste(names(correlation_results), collapse = ", "), "\n")
  }
  
  return(correlation_results)
}

# =============================================================================
# 6. SURVIVAL ANALYSIS FUNCTIONS  
# =============================================================================

#' Perform Cox Regression Analysis
#' @param survival_df Data frame prepared for survival analysis
#' @param config Analysis configuration
#' @param analysis_type Type of analysis: "univariate", "multivariate", or "both"
#' @param include_clinical Logical: include clinical variables in multivariate analysis
perform_cox_regression <- function(survival_df, config, analysis_type = "both",
                                 include_clinical = TRUE) {
  
  results <- list()
  celltype_cols <- get_celltype_columns(survival_df, config)
  survival_cols <- unlist(config$columns$survival)
  
  for (tissue in config$scope$tissues) {
    tissue_data <- survival_df %>% 
      filter(!!sym(config$columns$tissue) == tissue)
    
    if (nrow(tissue_data) < 20) {  # Minimum for Cox regression
      if (config$output$verbose) {
        cat("Skipping", tissue, "- insufficient patients for Cox regression\n")
      }
      next
    }
    
    tissue_results <- list()
    
    # Univariate analysis
    if (analysis_type %in% c("univariate", "both")) {
      uni_results <- perform_univariate_cox(tissue_data, celltype_cols, survival_cols, config)
      tissue_results$univariate <- uni_results
    }
    
    # Multivariate analysis
    if (analysis_type %in% c("multivariate", "both")) {
      multi_results <- perform_multivariate_cox(tissue_data, celltype_cols, 
                                              survival_cols, config, include_clinical)
      tissue_results$multivariate <- multi_results
    }
    
    results[[tissue]] <- tissue_results
  }
  
  if (config$output$verbose) {
    cat("✓ Cox regression analysis completed\n")
    cat("  - Analysis type:", analysis_type, "\n")
    cat("  - Tissues analyzed:", paste(names(results), collapse = ", "), "\n")
  }
  
  return(results)
}

# =============================================================================
# 7. UTILITY FUNCTIONS
# =============================================================================

#' Get Cell Type Column Names
#' @param data Data frame
#' @param config Analysis configuration
get_celltype_columns <- function(data, config) {
  # Identify cell type columns (typically start with capital letters or specific patterns)
  excluded_cols <- c(
    config$columns$roi, config$columns$tissue, config$columns$patient,
    config$columns$discrete_clinical, config$columns$continuous_clinical,
    unlist(config$columns$survival), "p_adjusted", "significance", 
    "is_significant", "fold_change", "p_value", "time", "event"
  )
  
  all_cols <- colnames(data)
  celltype_cols <- setdiff(all_cols, excluded_cols)
  
  # Further filter to remove obviously non-celltype columns
  celltype_cols <- celltype_cols[!grepl("^(p_|adj_|log2_|effect_)", celltype_cols)]
  
  return(celltype_cols)
}

#' Clean Variable Names for Display
#' @param variable_names Character vector of variable names
clean_variable_names <- function(variable_names) {
  cleaned <- variable_names
  
  # Define cleaning mappings
  name_mappings <- c(
    "RFS_status" = "Relapse vs Non-relapse",
    "Gender" = "Female vs Male", 
    "KRAS_mutation" = "KRAS Mut vs WT",
    "BRAF_mutation" = "BRAF Mut vs WT",
    "Differential_grade" = "Well vs Moderate Diff.",
    "Lymph_positive" = "LN+ vs LN-"
  )
  
  for (old_name in names(name_mappings)) {
    cleaned <- gsub(old_name, name_mappings[old_name], cleaned)
  }
  
  return(cleaned)
}

#' Calculate Effect Sizes
calculate_cohens_d <- function(x, y) {
  pooled_sd <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / 
                   (length(x) + length(y) - 2))
  (mean(y) - mean(x)) / pooled_sd
}

calculate_rank_biserial <- function(x, y) {
  # Simplified rank-biserial correlation
  n1 <- length(x)
  n2 <- length(y)
  U <- wilcox.test(x, y)$statistic
  (2 * U) / (n1 * n2) - 1
}

# =============================================================================
# 8. MISSING HELPER FUNCTIONS (NEED TO BE IMPLEMENTED)
# =============================================================================

#' Create Tissue Separator Bar for Stacked Plots
create_tissue_separator <- function(fractions_df, config, roi_order) {
  # Implementation needed for tissue separator in stacked bar plots
  tissue_data <- fractions_df %>%
    select(all_of(c(config$columns$roi, config$columns$tissue))) %>%
    distinct() %>%
    mutate(!!sym(config$columns$roi) := factor(!!sym(config$columns$roi), levels = roi_order))
  
  tissue_bar <- ggplot(tissue_data, aes(x = !!sym(config$columns$roi), y = 0.5, fill = !!sym(config$columns$tissue))) +
    geom_bar(stat = "identity", width = 1) +
    theme_void() +
    theme(legend.position = "bottom")
  
  return(tissue_bar)
}

#' Order ROIs by Similarity (Hierarchical Clustering)
order_rois_by_similarity <- function(fractions_df, celltype_cols) {
  # Implementation needed for ordering samples by similarity
  fraction_matrix <- fractions_df %>% 
    select(all_of(celltype_cols)) %>%
    as.matrix()
  
  if (nrow(fraction_matrix) > 1) {
    hc <- hclust(dist(fraction_matrix))
    return(fractions_df[[config$columns$roi]][hc$order])
  } else {
    return(fractions_df[[config$columns$roi]])
  }
}

#' Apply Fraction Threshold Filter
apply_fraction_threshold <- function(analysis_df, config) {
  # Implementation needed for filtering low-abundance cell types
  celltype_cols <- get_celltype_columns(analysis_df, config)
  
  for (col in celltype_cols) {
    analysis_df[[col]][analysis_df[[col]] < config$scope$min_fraction_threshold] <- 0
  }
  
  return(analysis_df)
}

#' Create Complex Clinical Heatmap
#' @param fractions_df Data frame with cell fractions and clinical data
#' @param config Analysis configuration
#' @param save_path Path to save the plot
#' @param order_by Variable to order patients by (default: survival event)
#' @param cluster_rows Logical: cluster rows (patients)
#' @param cluster_columns Logical: cluster columns (cell types)
#' @param show_row_names Logical: show patient names
#' @param show_column_names Logical: show cell type names
plot_clinical_heatmap <- function(fractions_df, config, save_path = NULL,
                                  order_by = NULL, cluster_rows = FALSE, 
                                  cluster_columns = TRUE, show_row_names = TRUE,
                                  show_column_names = TRUE, width, height) {
  
  # Set default ordering variable
  if (is.null(order_by)) {
    # Try to find survival event column
    if ("event" %in% colnames(fractions_df)) {
      order_by <- "event"
    } else if (config$columns$survival$event %in% colnames(fractions_df)) {
      order_by <- config$columns$survival$event
    } else {
      order_by <- config$columns$discrete_clinical[1]  # Use first discrete variable
    }
  }
  
  # Get cell type columns
  celltype_cols <- get_celltype_columns(fractions_df, config)
  
  # 1. Aggregate fractions by patient (mean across ROIs if needed)
  if ("sample_id" %in% colnames(fractions_df) && "patient_id" %in% colnames(fractions_df)) {
    # Data has sample-level fractions, need to aggregate by patient
    
    # Get clinical columns
    clinical_cols <- c(
      "patient_id", config$columns$discrete_clinical, 
      config$columns$continuous_clinical, unlist(config$columns$survival)
    )
    # Handle column name variations
    available_clinical_cols <- intersect(clinical_cols, colnames(fractions_df))[-1]
    
    # Add common alternative names
    alt_names <- c("event", "time", "Gender", "Age", "TBS", "CRLM_number", "CRLM_size", "CEA", "CA199")
    available_clinical_cols <- unique(c(available_clinical_cols, intersect(alt_names, colnames(fractions_df))))
    
    patient_fractions <- fractions_df %>%
      group_by(patient_id) %>%
      summarise(
        across(all_of(celltype_cols), mean, na.rm = TRUE),
        across(all_of(available_clinical_cols), first),  # Take first value for clinical vars
        .groups = "drop"
      )
  } else {
    # Data is already at patient level
    patient_fractions <- fractions_df
    # Ensure we have patient_id column
    if (!"patient_id" %in% colnames(patient_fractions) && config$columns$patient %in% colnames(patient_fractions)) {
      patient_fractions$patient_id <- patient_fractions[[config$columns$patient]]
    }
  }
  
  # 2. Create fraction matrix for heatmap
  fraction_matrix <- patient_fractions %>%
    select(patient_id, all_of(celltype_cols)) %>%
    column_to_rownames("patient_id") %>%
    as.matrix()
  
  # Remove any rows with all zeros
  fraction_matrix <- fraction_matrix[rowSums(fraction_matrix, na.rm = TRUE) > 0, ]
  
  # 3. Order patients by specified variable
  if (order_by %in% colnames(patient_fractions)) {
    patient_order <- patient_fractions %>%
      filter(patient_id %in% rownames(fraction_matrix)) %>%
      arrange(!!sym(order_by)) %>%
      pull(patient_id)
    
    fraction_matrix_ordered <- fraction_matrix[patient_order, ]
  } else {
    fraction_matrix_ordered <- fraction_matrix
    patient_order <- rownames(fraction_matrix)
  }
  
  # 4. Prepare clinical annotation data
  clinical_data <- patient_fractions %>%
    filter(patient_id %in% patient_order) %>%
    arrange(match(patient_id, patient_order)) %>%
    column_to_rownames("patient_id")
  
  # Identify discrete and continuous clinical variables
  discrete_vars <- intersect(config$columns$discrete_clinical,colnames(clinical_data))
  continuous_vars <- intersect(config$columns$continuous_clinical, colnames(clinical_data))
  
  # Remove duplicates and ensure we have the ordering variable
  discrete_vars <- unique(discrete_vars)
  continuous_vars <- unique(continuous_vars)
  
  # 5. Create color schemes for annotations
  
  # Define colors for discrete variables
  discrete_colors <- list()
  
  for (var in discrete_vars) {
    unique_vals <- unique(clinical_data[[var]])
    unique_vals <- unique_vals[!is.na(unique_vals)]
    
    if (var %in% c("event", "RFS_status")) {
      discrete_colors[[var]] <- c("0" = "#bc3c29ff", "1" = "#0072b5ff")
    } else if (var == "Gender") {
      discrete_colors[[var]] <- c("0" = "#7876b1ff", "1" = "#ee4c97ff")
    } else if (var %in% c("KRAS_mutation", "BRAF_mutation")) {
      discrete_colors[[var]] <- c("0" = "grey70", "1" = "#e18727ff")
    } else if (var == "Differential_grade") {
      discrete_colors[[var]] <- c("0" = "lightblue", "1" = "darkblue")
    } else if (var == "T_stage") {
      stage_colors <- c("lemonchiffon", "lemonchiffon2", "lemonchiffon3", "lemonchiffon4")
      names(stage_colors) <- as.character(sort(unique_vals))
      discrete_colors[[var]] <- stage_colors[1:length(unique_vals)]
    } else if (var == "Lymph_positive") {
      discrete_colors[[var]] <- c("0" = "white", "1" = "grey60")
    } else {
      # Default colors for other variables
      n_vals <- length(unique_vals)
      colors <- RColorBrewer::brewer.pal(min(max(3, n_vals), 8), "Set2")[1:n_vals]
      names(colors) <- as.character(unique_vals)
      discrete_colors[[var]] <- colors
    }
  }
  
  # Define colors for continuous variables
  continuous_colors <- list()
  
  for (var in continuous_vars) {
    var_range <- range(clinical_data[[var]], na.rm = TRUE)
    
    if (var %in% c("time", "RFS_time")) {
      continuous_colors[[var]] <- colorRamp2(var_range, c("white", "red4"))
    } else if (var == "Age") {
      continuous_colors[[var]] <- colorRamp2(var_range, c("white", "slategray"))
    } else if (var == "TBS") {
      continuous_colors[[var]] <- colorRamp2(var_range, c("white", "yellowgreen"))
    } else if (var == "CRLM_number") {
      continuous_colors[[var]] <- colorRamp2(var_range, c("white", "orangered1"))
    } else if (var == "CRLM_size") {
      continuous_colors[[var]] <- colorRamp2(var_range, c("white", "orangered3"))
    } else if (var == "CEA") {
      continuous_colors[[var]] <- colorRamp2(var_range, c("white", "tomato3"))
    } else if (var == "CA199") {
      continuous_colors[[var]] <- colorRamp2(var_range, c("white", "tomato1"))
    } else {
      # Default color scheme
      continuous_colors[[var]] <- colorRamp2(var_range, c("white", "darkblue"))
    }
  }
  
  # 6. Create heatmap annotations
  
  # Discrete annotations
  if (length(discrete_vars) > 0) {
    ha_discrete <- HeatmapAnnotation(
      df = clinical_data[, discrete_vars, drop = FALSE],
      col = discrete_colors[discrete_vars],
      which = "row",
      annotation_legend_param = list(title_position = "topcenter"),
      show_annotation_name = TRUE
    )
  } else {
    ha_discrete <- NULL
  }
  
  # Continuous annotations
  if (length(continuous_vars) > 0) {
    ha_continuous <- HeatmapAnnotation(
      df = clinical_data[, continuous_vars, drop = FALSE],
      col = continuous_colors[continuous_vars],
      which = "row",
      show_annotation_name = TRUE
    )
  } else {
    ha_continuous <- NULL
  }
  
  # Combine annotations
  if (!is.null(ha_discrete) && !is.null(ha_continuous)) {
    ha_combined <- ha_discrete + ha_continuous
  } else if (!is.null(ha_discrete)) {
    ha_combined <- ha_discrete
  } else if (!is.null(ha_continuous)) {
    ha_combined <- ha_continuous
  } else {
    ha_combined <- NULL
  }
  
  # 7. Create main heatmap
  
  # Color scheme for fractions
  max_fraction <- max(fraction_matrix_ordered, na.rm = TRUE)
  heatmap_colors <- colorRamp2(c(0, max_fraction), c("white", "navy"))
  
  ht_main <- Heatmap(
    fraction_matrix_ordered,
    name = "Cell Fraction",
    col = heatmap_colors,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    column_names_side = "bottom",
    row_names_side = "right",
    column_names_rot = 90,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(
      title = "Cell Fraction",
      title_position = "topcenter",
      legend_direction = "vertical"
    )
  )
  
  # 8. Combine heatmap with annotations
  if (!is.null(ha_combined)) {
    final_plot <- ha_combined + ht_main
  } else {
    final_plot <- ht_main
  }
  
  # 9. Save plot if path provided
  if (!is.null(save_path)) {
    pdf(save_path, width = width, height = height)
    draw(final_plot)
    dev.off()
    
    if (config$output$verbose) {
      cat("✓ Clinical heatmap saved to:", save_path, "\n")
    }
  }
  
  return(final_plot)
}


#' Calculate Cellular Diversity Metrics (Entropy and Simpson Index)
#' @param analysis_data Data frame with cell fractions per ROI
#' @param config Analysis configuration
#' @param add_scaled Logical: add scaled versions of metrics
#' @param verbose Logical: print progress messages
calculate_cellular_diversity <- function(analysis_data, config, add_scaled = TRUE, verbose = TRUE) {
  
  if (config$output$verbose) {
    cat("Calculating cellular diversity metrics...\n")
  }
  
  # Get cell type columns
  celltype_cols <- get_celltype_columns(analysis_data, config)
  
  if (length(celltype_cols) == 0) {
    stop("No cell type columns found in data")
  }
  
  # Extract cell fraction matrix
  fraction_matrix <- analysis_data[, celltype_cols]
  
  # Calculate Shannon Entropy for each ROI
  if (config$output$verbose) cat("  - Calculating Shannon entropy...\n")
  analysis_data$entropy <- apply(fraction_matrix, 1, function(x) {
    p <- x[x > 0]  # Keep only positive proportions
    if (length(p) == 0 || sum(p) == 0) {
      return(0)  # If all proportions are 0, entropy is 0
    } else {
      # Normalize to ensure sum = 1 (in case of small rounding errors)
      p <- p / sum(p)
      return(-sum(p * log(p)))  # Shannon entropy with natural log
    }
  })
  
  # Calculate Simpson Index for each ROI
  if (config$output$verbose) cat("  - Calculating Simpson index...\n")
  analysis_data$simpson <- apply(fraction_matrix, 1, function(x) {
    if (sum(x) == 0) {
      return(0)
    } else {
      # Normalize to ensure sum = 1
      p <- x / sum(x)
      return(sum(p^2))  # Simpson index D
    }
  })
  
  # Calculate Simpson Diversity Index (1 - D)
  analysis_data$simpson_diversity <- 1 - analysis_data$simpson
  
  # Calculate Inverse Simpson Index (1/D)
  analysis_data$inverse_simpson <- ifelse(analysis_data$simpson > 0, 
                                          1 / analysis_data$simpson, 
                                          0)
  
  # Add scaled versions if requested
  if (add_scaled) {
    if (config$output$verbose) cat("  - Creating scaled versions...\n")
    analysis_data$scale_entropy <- as.numeric(scale(analysis_data$entropy))
    analysis_data$scale_simpson <- as.numeric(scale(analysis_data$simpson))
    analysis_data$scale_simpson_diversity <- as.numeric(scale(analysis_data$simpson_diversity))
    analysis_data$scale_inverse_simpson <- as.numeric(scale(analysis_data$inverse_simpson))
  }
  
  # Add summary statistics
  if (config$output$verbose) {
    cat("✓ Diversity metrics calculated\n")
    cat("  - Entropy range:", round(range(analysis_data$entropy, na.rm = TRUE), 3), "\n")
    cat("  - Simpson index range:", round(range(analysis_data$simpson, na.rm = TRUE), 3), "\n")
    cat("  - Mean diversity per tissue:\n")
    
    diversity_summary <- analysis_data %>%
      group_by(!!sym(config$columns$tissue)) %>%
      summarise(
        mean_entropy = round(mean(entropy, na.rm = TRUE), 3),
        mean_simpson = round(mean(simpson, na.rm = TRUE), 3),
        .groups = "drop"
      )
    
    print(diversity_summary)
  }
  
  return(analysis_data)
}

#' Plot Cellular Diversity Metrics by Clinical Groups
#' @param analysis_data Data frame with diversity metrics
#' @param clinical_var Clinical variable to compare
#' @param config Analysis configuration
#' @param save_path Path to save the plot
#' @param use_scaled Logical: use scaled metrics for plotting
#' @param plot_width Plot width
#' @param plot_height Plot height
plot_cellular_diversity <- function(analysis_data, clinical_var, config, 
                                    save_path = NULL, use_scaled = TRUE,
                                    plot_width = 12, plot_height = 6) {
  
  # Validate clinical variable
  if (!clinical_var %in% colnames(analysis_data)) {
    stop("Clinical variable '", clinical_var, "' not found in data")
  }
  
  # Prepare data
  plot_data <- analysis_data %>%
    filter(!is.na(!!sym(clinical_var)))
  
  if (nrow(plot_data) == 0) {
    stop("No data available after removing missing values for ", clinical_var)
  }
  
  # Determine which metrics to plot
  if (use_scaled) {
    entropy_col <- "scale_entropy"
    simpson_col <- "scale_simpson"
    y_label_entropy <- "Scaled Shannon Entropy"
    y_label_simpson <- "Scaled Simpson Index"
  } else {
    entropy_col <- "entropy"
    simpson_col <- "simpson"
    y_label_entropy <- "Shannon Entropy"
    y_label_simpson <- "Simpson Index"
  }
  
  # Check if columns exist
  if (!entropy_col %in% colnames(plot_data)) {
    stop("Column '", entropy_col, "' not found. Did you calculate diversity metrics with add_scaled=TRUE?")
  }
  
  # Convert clinical variable to factor
  plot_data[[clinical_var]] <- as.factor(plot_data[[clinical_var]])
  n_groups <- length(unique(plot_data[[clinical_var]]))
  
  # Get appropriate colors
  colors <- pal_jco("default")(min(n_groups, 10))
  
  # Create entropy plot
  p1 <- ggplot(plot_data, aes(x = Tissue, y = !!sym(entropy_col), 
                              fill = !!sym(clinical_var))) +
    geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
    geom_jitter(color = "darkgrey", position = position_jitter(width = 0.2), 
                size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      text = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 11),
      strip.background = element_blank(),
      legend.position = "none"
    ) +
    labs(
      title = paste0("Shannon Entropy by ", clean_variable_names(clinical_var)),
      x = clean_variable_names(clinical_var),
      y = y_label_entropy
    )
  
  # Add statistical comparisons
  if (n_groups == 2) {
    p1 <- p1 + stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      label.x.npc = "center"
    )
  } else if (n_groups > 2) {
    # For multiple groups, add overall p-value
    p1 <- p1 + stat_compare_means(
      method = "kruskal.test",
      label = "p.format",
      label.y.npc = "top"
    )
  }
  
  # Create Simpson index plot
  p2 <- ggplot(plot_data, aes(x = Tissue, y = !!sym(simpson_col), 
                              fill = !!sym(clinical_var))) +
    geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
    geom_jitter(color = "darkgrey", position = position_jitter(width = 0.2), 
                size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      text = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 11),
      strip.background = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = paste0("Simpson Index by ", clean_variable_names(clinical_var)),
      x = clean_variable_names(clinical_var),
      y = y_label_simpson
    )

  # Add statistical comparisons
  if (n_groups == 2) {
    p2 <- p2 + stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      label.x.npc = "center"
    )
  } else if (n_groups > 2) {
    p2 <- p2 + stat_compare_means(
      method = "kruskal.test",
      label = "p.format",
      label.y.npc = "top"
    )
  }
  
  # Combine plots
  combined_plot <- p1 + p2
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, 
           width = plot_width, height = plot_height, dpi = config$plots$dpi)
    
    if (config$output$verbose) {
      cat("✓ Diversity plots saved:", basename(save_path), "\n")
    }
  }
  
  return(combined_plot)
}

#' Calculate Tissue Correlations with Significance Testing
#'
#' This function uses Hmisc::rcorr to calculate a correlation matrix and a
#' corresponding p-value matrix.
#'
#' @param tissue_data Data frame for a single tissue.
#' @param celltype_cols A character vector of cell type column names.
#' @param continuous_vars A character vector of continuous variable column names.
#' @param method Correlation method (e.g., "pearson", "spearman").
#' @return A list containing 'cor_matrix' (correlation coefficients) and
#'   'p_matrix' (p-values).
#'
calculate_tissue_correlations <- function(tissue_data, celltype_cols, continuous_vars, method, min_obs) {
  
  # Hmisc::rcorr requires a matrix and computes correlations for all column pairs
  data_for_cor <- as.matrix(tissue_data[, c(celltype_cols, continuous_vars)])
  
  # Perform correlation analysis
  # The result is a list containing matrices for r, n, and P (p-values)
  corr_test_result <- Hmisc::rcorr(data_for_cor, type = method)
  
  # Extract the correlation coefficients (r) and p-values (P)
  cor_matrix_full <- corr_test_result$r
  p_matrix_full <- corr_test_result$P
  
  # Subset the matrices to show correlations between cell types and clinical variables
  cor_matrix <- cor_matrix_full[celltype_cols, continuous_vars, drop = FALSE]
  p_matrix <- p_matrix_full[celltype_cols, continuous_vars, drop = FALSE]
  
  if (!is.null(config$stats$p_adjust_method)){
    p_matrix <- apply(p_matrix, MARGIN = 2, function(x){
      return(p.adjust(x, method = config$stats$p_adjust_method))
    })
  }
  
  # Return a list containing both the correlation and p-value matrices
  return(list(cor_matrix = cor_matrix, p_matrix = p_matrix))
}

#' Generate Correlation Scatter Plots for Significant Pairs
#'
#' This function iterates through significant correlations found in the analysis
#' and generates a scatter plot with a regression line for each pair.
#'
#' @param analysis_data The full data frame containing cell fractions and clinical data.
#' @param correlations The list object from `calculate_abundance_correlations`,
#'   which contains correlation and p-value matrices for each tissue.
#' @param config The analysis configuration object.
#' @param output_dir The base directory where plot folders will be created.
#' @param significance_threshold The significance threshold to select pairs for plotting.
#'   Defaults to 0.05. For robust analysis, you should use this threshold on
#'   p-values that have been adjusted for multiple comparisons (e.g., with the BH method).
#' @return This function is called for its side effect of saving plots and does
#'   not return a value.
#'
generate_correlation_scatter_plots <- function(fractions_df,
                                               correlation_matrix,
                                               tissue_name,
                                               config,
                                               output_dir,
                                               significance_threshold = 0.05,
                                               width, height) {
  
  if (config$output$verbose) {
    cat("📈 Generating scatter plots for significant correlations (p <", significance_threshold, ")...\n")
  }
  
  plot_count <- 0
  
  # 2. Loop through each tissue in the correlation results
    p_matrix <- correlation_matrix$p_matrix
    
    # Skip if there's no p-value matrix for the tissue
    if (is.null(p_matrix) || all(is.na(p_matrix))) {
      next
    }
    
    # 3. Find pairs with p-values below the threshold
    significant_pairs <- which(p_matrix <= significance_threshold, arr.ind = TRUE)
    
    if (nrow(significant_pairs) == 0) {
      next # No significant correlations for this tissue
    }
    
    # Filter the main dataframe to get data for the current tissue
    tissue_data <- analysis_data %>%
      dplyr::filter(!!sym(config$columns$tissue) == tissue)
    
    # Create a dedicated subdirectory for the tissue's plots
    tissue_output_dir <- file.path(output_dir, "scatter_plots", tissue)
    if (!dir.exists(tissue_output_dir)) {
      dir.create(tissue_output_dir, recursive = TRUE)
    }
    
    # 4. Loop through each significant pair to create a plot
    for (i in 1:nrow(significant_pairs)) {
      cell_type_name <- rownames(p_matrix)[significant_pairs[i, "row"]]
      clinical_var_name <- colnames(p_matrix)[significant_pairs[i, "col"]]
      
      # Generate the plot object
      p <- ggplot2::ggplot(tissue_data, ggplot2::aes(x = !!sym(clinical_var_name), y = !!sym(cell_type_name))) +
        ggplot2::geom_point(alpha = 0.7, color = "steelblue", size = 2) +
        ggplot2::geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
        ggpubr::stat_cor(
          method = config$stats$correlation_method,
          label.x.npc = "left",   # Position stats on the top-left
          label.y.npc = "top"
        ) +
        ggplot2::labs(
          title = paste(cell_type_name, "vs.", clinical_var_name),
          subtitle = paste("Tissue:", tissue),
          x = clinical_var_name,
          y = paste(cell_type_name, "Abundance")
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = ggplot2::element_text(hjust = 0.5)
        )
      
      # 5. Define a clean filename and save the plot
      safe_filename <- paste0("scatter_", gsub("[^a-zA-Z0-9_.-]", "_", cell_type_name), "_vs_", gsub("[^a-zA-Z0-9_.-]", "_", clinical_var_name, "_",tissue_name), ".pdf")
      save_path <- file.path(tissue_output_dir, safe_filename)
      
      ggplot2::ggsave(save_path, plot = p, width = width, height = height, device = "pdf")
      plot_count <- plot_count + 1
    }
  
  
  if (config$output$verbose) {
    if (plot_count > 0) {
      cat("✓ Successfully generated", plot_count, "scatter plot(s).\n")
    } else {
      cat("✓ No significant correlations found to generate plots.\n")
    }
  }
  
  return(invisible(NULL))
}

#' Analyze Cell Type Correlations Within Tissues
#' @param fractions_df Data frame with cell fractions
#' @param config Analysis configuration
#' @param output_dir Output directory for plots
#' @param method Correlation method ("spearman", "pearson")
#' @param min_correlation Minimum correlation to highlight
analyze_celltype_correlations <- function(fractions_df, config, output_dir, 
                                          method = "spearman", min_correlation = 0.3) {
  
  if (config$output$verbose) {
    cat("Analyzing cell type correlations within tissues...\n")
  }
  
  celltype_cols <- get_celltype_columns(fractions_df, config)
  correlation_results <- list()
  
  for (tissue in config$scope$tissues) {
    if (config$output$verbose) {
      cat("  - Processing", tissue, "tissue...\n")
    }
    
    tissue_data <- fractions_df %>% 
      filter(!!sym(config$columns$tissue) == tissue) %>%
      select(all_of(celltype_cols))
    
    tissue_data <- as.data.frame(tissue_data)
    
    if (nrow(tissue_data) < 10) {
      next
    }
    
    # Calculate correlations
    cor_matrix <- cor(tissue_data, method = method, use = "complete.obs")
    
    # Create p-value matrix
    p_matrix <- matrix(NA, nrow = ncol(tissue_data), ncol = ncol(tissue_data))
    rownames(p_matrix) <- colnames(p_matrix) <- colnames(tissue_data)
    
    for (i in 1:ncol(tissue_data)) {
      for (j in 1:ncol(tissue_data)) {
        if (i != j) {
          vec_1 <- as.numeric(tissue_data[,i])
          vec_2 <- as.numeric(tissue_data[,j])
          
          test_result <- cor.test(vec_1, vec_2, method = method)
          p_matrix[i,j] <- test_result$p.value
        }
      }
    }
    
    if (!is.null(config$stats$p_adjust_method)){
      p_vector <- p_matrix[upper.tri(p_matrix)] 
      p_adjusted_vector <- p.adjust(p_vector, method = config$stats$p_adjust_method)
      
      p_matrix[upper.tri(p_matrix)]  <- p_adjusted_vector
      p_matrix[lower.tri(p_matrix)]  <- p_adjusted_vector
  }
    
    # Store results
    correlation_results[[tissue]] <- list(
      correlations = cor_matrix,
      p_values = p_matrix,
      n_samples = nrow(tissue_data)
    )
    
    # Create heatmap
    create_celltype_correlation_heatmap(
      cor_matrix, p_matrix, tissue, method, config, 
      file.path(output_dir, paste0("celltype_correlations_", tissue, ".pdf")),
      height = 9,width = 12
    )
    
    # Find strong correlations
    strong_correlations <- find_strong_correlations(
      cor_matrix, p_matrix, tissue, min_correlation, config
    )
    
    if (nrow(strong_correlations) > 0) {
      write.csv(strong_correlations, 
                file.path(output_dir, paste0("strong_correlations_", tissue, ".csv")),
                row.names = FALSE)
    }
  }
  
  return(correlation_results)
}

#' Create a Cell Type Correlation Heatmap with ggcorrplot
#'
#' This function generates a modern correlation heatmap using the ggcorrplot package.
#' It displays the lower triangle of the correlation matrix, using colored circles
#' for visual representation and overlaying numeric coefficients. Non-significant
#' correlations (p > 0.05) are left blank.
#'
#' @param cor_matrix A numeric matrix of correlation coefficients.
#' @param p_matrix A numeric matrix of corresponding p-values.
#' @param tissue A character string for the plot title (e.g., "Brain").
#' @param method A character string for the correlation method (e.g., "pearson").
#' @param config A list containing configuration, including `config$output$verbose`.
#' @param save_path The full file path to save the PDF output.
#' @param width The width of the output PDF.
#' @param height The height of the output PDF.
#' @return A ggplot object representing the correlation heatmap.

create_celltype_correlation_heatmap <- function(cor_matrix, p_matrix, tissue,
                                                    method, config, save_path, width, height) {
  
  p_matrix[is.na(p_matrix)] <- 0
  
  # Create the correlation heatmap
  p <- ggcorrplot(
    corr = cor_matrix,
    p.mat = p_matrix,
    hc.order = FALSE,                   # Reorder based on hierarchical clustering
    type = "lower",                    # Display the lower triangle
    method = "circle",                 # Use circles for correlation strength
    insig = "blank",                   # Leave non-significant cells (p > 0.05) blank
    lab = TRUE,                        # Add correlation coefficients
    lab_size = 3.5,                    # Set the label size
    colors = c("#4477AA", "white", "#CC3311"), # Set a colorblind-friendly palette
    ggtheme = theme_minimal()          # Use a minimal theme
  ) +
    labs(
      title = paste("Cell Type Correlations in", tissue, "Tissue"),
      subtitle = paste(stringr::str_to_title(method), "Correlation"),
      x = "", # Remove x-axis label
      y = ""  # Remove y-axis label
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  
  # Save the plot using ggsave for better quality and control
  ggsave(save_path, plot = p, width = width, height = height, device = "pdf")
  
  if (config$output$verbose) {
    cat("✓ ggcorrplot heatmap saved:", basename(save_path), "\n")
  }
  
  return(p)
}

#' Find Strong Cell Type Correlations
find_strong_correlations <- function(cor_matrix, p_matrix, tissue, 
                                     min_correlation, config) {
  
  strong_cors <- data.frame(
    Tissue = character(),
    CellType1 = character(),
    CellType2 = character(),
    Correlation = numeric(),
    P_value = numeric(),
    Relationship = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:(nrow(cor_matrix)-1)) {
    for (j in (i+1):ncol(cor_matrix)) {
      cor_val <- cor_matrix[i,j]
      p_val <- p_matrix[i,j]
      
      if (!is.na(cor_val) && !is.na(p_val) && abs(cor_val) >= min_correlation && p_val < 0.05) {
        relationship <- ifelse(cor_val > 0, "Positive", "Negative")
        
        strong_cors <- rbind(strong_cors, data.frame(
          Tissue = tissue,
          CellType1 = rownames(cor_matrix)[i],
          CellType2 = colnames(cor_matrix)[j],
          Correlation = cor_val,
          P_value = p_val,
          Relationship = relationship,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Sort by absolute correlation
  if (nrow(strong_cors) > 0) {
    strong_cors <- strong_cors[order(abs(strong_cors$Correlation), decreasing = TRUE), ]
    
    if (config$output$verbose) {
      cat("    Strong correlations found:", nrow(strong_cors), "\n")
      if (nrow(strong_cors) > 0) {
        top_cor <- strong_cors[1, ]
        cat("    Strongest:", top_cor$CellType1, "↔", top_cor$CellType2, 
            "=", round(top_cor$Correlation, 3), "\n")
      }
    }
  }
  
  return(strong_cors)
}

#' Create a Tidy Summary of All Correlation Results
#'
#' This function consolidates correlation results from all tissues into a single
#' data frame, calculates adjusted p-values, and saves the output to a CSV file.
#'
#' @param correlations The list object from `calculate_abundance_correlations`,
#'   containing correlation and p-value matrices for each tissue.
#' @param config The analysis configuration object.
#' @param save_path The full file path to save the summary CSV file.
#' @return A tidy data frame summarizing all correlation results, returned invisibly.
#'
create_correlation_summary <- function(correlations, config, save_path) {
  
  # 1. Check for required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required. Please install it.", call. = FALSE)
  }
  
  # List to hold the tidy data frame for each tissue
  summary_list <- list()
  
  # 2. Loop through each tissue to process its results
  for (tissue in names(correlations)) {
    cor_matrix <- correlations[[tissue]]$cor_matrix
    p_matrix <- correlations[[tissue]]$p_matrix
    
    # Skip if data is missing for the tissue
    if (is.null(cor_matrix) || is.null(p_matrix)) {
      next
    }
    
    # Convert the correlation matrix to a long format
    cor_df <- as.data.frame(cor_matrix) %>%
      dplyr::mutate(CellType = rownames(cor_matrix), .before = 1) %>%
      tidyr::pivot_longer(
        cols = -CellType,
        names_to = "ClinicalVariable",
        values_to = "CorrelationCoefficient"
      )
    
    # Convert the p-value matrix to a long format
    p_df <- as.data.frame(p_matrix) %>%
      dplyr::mutate(CellType = rownames(p_matrix), .before = 1) %>%
      tidyr::pivot_longer(
        cols = -CellType,
        names_to = "ClinicalVariable",
        values_to = "PValue"
      )
    
    # Join the two data frames and add the tissue name
    tissue_summary <- dplyr::left_join(cor_df, p_df, by = c("CellType", "ClinicalVariable")) %>%
      dplyr::mutate(Tissue = tissue, .before = 1)
    
    summary_list[[tissue]] <- tissue_summary
  }
  
  # 3. Combine all tissue summaries into a single data frame
  if (length(summary_list) == 0) {
    if (config$output$verbose) {
      cat("No correlation results to summarize.\n")
    }
    return(invisible(NULL))
  }
  
  full_summary <- dplyr::bind_rows(summary_list)
  
  # 4. Calculate adjusted p-values and add a significance column
  # The adjustment is performed on all p-values across all tissues at once.
  final_summary <- full_summary %>%
    dplyr::mutate(
      PValueAdjusted = p.adjust(PValue, method = "BH"),
      Significance = ifelse(PValueAdjusted < 0.001, "***",
                            ifelse(PValueAdjusted < 0.01, "**",
                                   ifelse(PValueAdjusted < 0.05, "*", "ns")))
    ) %>%
    # Sort by significance for easy review
    dplyr::arrange(PValueAdjusted)
  
  # 5. Save the final summary to a CSV file
  readr::write_csv(final_summary, save_path)
  
  if (config$output$verbose) {
    cat("📊 Correlation summary saved to:", save_path, "\n")
  }
  
  return(invisible(final_summary))
}

# =============================================================================
# 4. DISCRETE CLINICAL ANALYSIS FUNCTIONS
# =============================================================================

#' Perform Statistical Comparison for Single Variable
#' @param data Data frame with fractions and grouping variable
#' @param var_name Name of grouping variable
#' @param celltype_cols Names of cell type columns
#' @param test_method Statistical test method
#' @param effect_size Logical: calculate effect sizes
#' @param config Analysis configuration
perform_variable_comparison <- function(data, var_name, celltype_cols, 
                                        test_method, effect_size, config) {
  
  groups <- unique(data[[var_name]])
  groups <- groups[!is.na(groups)]
  
  if (length(groups) < 2) {
    return(NULL)
  }
  
  # For now, handle binary comparisons
  if (length(groups) == 2) {
    groups <- sort(groups)
    group1_data <- data[data[[var_name]] == groups[1], ]
    group2_data <- data[data[[var_name]] == groups[2], ]
    
    # Check minimum group sizes
    if (nrow(group1_data) < config$stats$min_cells_per_group || 
        nrow(group2_data) < config$stats$min_cells_per_group) {
      return(NULL)
    }
    
    results <- data.frame(
      CellType = celltype_cols,
      fold_change = numeric(length(celltype_cols)),
      p_value = numeric(length(celltype_cols)),
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(celltype_cols)) {
      celltype <- celltype_cols[i]
      
      v1 <- group1_data[[celltype]]
      v2 <- group2_data[[celltype]]
      
      # Remove NA values
      v1 <- v1[!is.na(v1)]
      v2 <- v2[!is.na(v2)]
      
      if (length(v1) < 3 || length(v2) < 3) {
        results$fold_change[i] <- NA
        results$p_value[i] <- NA
        next
      }
      
      # Calculate fold change (group2 vs group1)
      mean1 <- mean(v1)
      mean2 <- mean(v2)
      results$fold_change[i] <- ifelse(mean1 > 0, mean2 / mean1, NA)
      
      # Perform statistical test
      if (test_method == "wilcox.test") {
        test_result <- wilcox.test(v2, v1)
      } else if (test_method == "t.test") {
        test_result <- t.test(v2, v1)
      } else {
        stop("Unsupported test method: ", test_method)
      }
      
      results$p_value[i] <- test_result$p.value
      
      # Add effect size if requested
      if (effect_size) {
        if (test_method == "wilcox.test") {
          # Calculate rank-biserial correlation
          results$effect_size[i] <- calculate_rank_biserial(v1, v2)
        } else {
          # Calculate Cohen's d
          results$effect_size[i] <- calculate_cohens_d(v1, v2)
        }
      }
    }
    
    return(results)
  }
  
  return(NULL)
}

# =============================================================================
# 5. CONTINUOUS CLINICAL ANALYSIS FUNCTIONS
# =============================================================================

#' Calculate Correlations with Continuous Variables
#' @param fractions_df Data frame with cell fractions
#' @param config Analysis configuration
#' @param method Correlation method
#' @param min_observations Minimum observations for correlation
calculate_abundance_correlations <- function(fractions_df, config, 
                                             method = NULL, min_observations = 10) {
  
  if (is.null(method)) {
    method <- config$stats$correlation_method
  }
  
  correlation_results <- list()
  celltype_cols <- get_celltype_columns(fractions_df, config)
  
  for (tissue in config$scope$tissues) {
    tissue_data <- fractions_df %>% 
      filter(!!sym(config$columns$tissue) == tissue)
    
    if (nrow(tissue_data) < min_observations) {
      if (config$output$verbose) {
        cat("Skipping", tissue, "- insufficient observations for correlation\n")
      }
      next
    }
    
    # Calculate correlations
    cor_results <- calculate_tissue_correlations(
      tissue_data, celltype_cols, config$columns$continuous_clinical, 
      method, min_observations
    )
    
    if (!is.null(cor_results)) {
      correlation_results[[tissue]] <- cor_results
    }
  }
  
  if (config$output$verbose && length(correlation_results) > 0) {
    cat("✓ Correlation analysis completed\n")
    cat("  - Tissues analyzed:", paste(names(correlation_results), collapse = ", "), "\n")
  }
  
  return(correlation_results)
}

# =============================================================================
# 6. SURVIVAL ANALYSIS FUNCTIONS  
# =============================================================================

#' Perform Cox Regression Analysis
#' @param survival_df Data frame prepared for survival analysis
#' @param config Analysis configuration
#' @param analysis_type Type of analysis: "univariate", "multivariate", or "both"
#' @param include_clinical Logical: include clinical variables in multivariate analysis
perform_cox_regression <- function(survival_df, config, analysis_type = "both",
                                   include_clinical = TRUE) {
  
  results <- list()
  celltype_cols <- get_celltype_columns(survival_df, config)
  survival_cols <- unlist(config$columns$survival)
  
  for (tissue in config$scope$tissues) {
    tissue_data <- survival_df %>% 
      filter(!!sym(config$columns$tissue) == tissue)
    
    if (nrow(tissue_data) < 20) {  # Minimum for Cox regression
      if (config$output$verbose) {
        cat("Skipping", tissue, "- insufficient patients for Cox regression\n")
      }
      next
    }
    
    tissue_results <- list()
    
    # Univariate analysis
    if (analysis_type %in% c("univariate", "both")) {
      uni_results <- perform_univariate_cox(tissue_data, celltype_cols, survival_cols, config)
      tissue_results$univariate <- uni_results
    }
    
    # Multivariate analysis
    if (analysis_type %in% c("multivariate", "both")) {
      multi_results <- perform_multivariate_cox(tissue_data, celltype_cols, 
                                                survival_cols, config, include_clinical)
      tissue_results$multivariate <- multi_results
    }
    
    results[[tissue]] <- tissue_results
  }
  
  if (config$output$verbose) {
    cat("✓ Cox regression analysis completed\n")
    cat("  - Analysis type:", analysis_type, "\n")
    cat("  - Tissues analyzed:", paste(names(results), collapse = ", "), "\n")
  }
  
  return(results)
}

# =============================================================================
# 7. UTILITY FUNCTIONS
# =============================================================================

#' Clean Variable Names for Display
#' @param variable_names Character vector of variable names
clean_variable_names <- function(variable_names) {
  cleaned <- variable_names
  
  # Define cleaning mappings
  name_mappings <- c(
    "RFS_status" = "Relapse vs Non-relapse",
    "Gender" = "Female vs Male", 
    "KRAS_mutation" = "KRAS Mut vs WT",
    "BRAF_mutation" = "BRAF Mut vs WT",
    "Differential_grade" = "Well vs Moderate Diff.",
    "Lymph_positive" = "LN+ vs LN-"
  )
  
  for (old_name in names(name_mappings)) {
    cleaned <- gsub(old_name, name_mappings[old_name], cleaned)
  }
  
  return(cleaned)
}

#' Calculate Effect Sizes
calculate_cohens_d <- function(x, y) {
  pooled_sd <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / 
                      (length(x) + length(y) - 2))
  (mean(y) - mean(x)) / pooled_sd
}

calculate_rank_biserial <- function(x, y) {
  # Simplified rank-biserial correlation
  n1 <- length(x)
  n2 <- length(y)
  U <- wilcox.test(x, y)$statistic
  (2 * U) / (n1 * n2) - 1
}

# =============================================================================
# 8. MISSING HELPER FUNCTIONS (NEED TO BE IMPLEMENTED)
# =============================================================================

#' Create Tissue Separator Bar for Stacked Plots
create_tissue_separator <- function(fractions_df, config, roi_order) {
  # Implementation needed for tissue separator in stacked bar plots
  tissue_data <- fractions_df %>%
    select(all_of(c(config$columns$roi, config$columns$tissue))) %>%
    distinct() %>%
    mutate(!!sym(config$columns$roi) := factor(!!sym(config$columns$roi), levels = roi_order))
  
  tissue_bar <- ggplot(tissue_data, aes(x = !!sym(config$columns$roi), y = 0.5, fill = !!sym(config$columns$tissue))) +
    geom_bar(stat = "identity", width = 1) +
    theme_void() +
    theme(legend.position = "bottom")
  
  return(tissue_bar)
}

#' Order ROIs by Similarity (Hierarchical Clustering)
order_rois_by_similarity <- function(fractions_df, celltype_cols) {
  # Implementation needed for ordering samples by similarity
  fraction_matrix <- fractions_df %>% 
    select(all_of(celltype_cols)) %>%
    as.matrix()
  
  if (nrow(fraction_matrix) > 1) {
    hc <- hclust(dist(fraction_matrix))
    return(fractions_df[[config$columns$roi]][hc$order])
  } else {
    return(fractions_df[[config$columns$roi]])
  }
}

#' Apply Fraction Threshold Filter
apply_fraction_threshold <- function(analysis_df, config) {
  # Implementation needed for filtering low-abundance cell types
  celltype_cols <- get_celltype_columns(analysis_df, config)
  
  for (col in celltype_cols) {
    analysis_df[[col]][analysis_df[[col]] < config$scope$min_fraction_threshold] <- 0
  }
  
  return(analysis_df)
}

#' Calculate Summary Statistics
calculate_summary_statistics <- function(analysis_data, config) {
  # Implementation needed for summary statistics
  celltype_cols <- get_celltype_columns(analysis_data, config)
  
  summary_stats <- list()
  
  for (tissue in config$scope$tissues) {
    tissue_data <- analysis_data %>% filter(!!sym(config$columns$tissue) == tissue)
    
    tissue_summary <- tissue_data %>%
      select(all_of(celltype_cols)) %>%
      summarise_all(list(mean = mean, median = median, sd = sd), na.rm = TRUE)
    
    summary_stats[[tissue]] <- tissue_summary
  }
  
  return(summary_stats)
}

#' Generate Volcano Plot for Abundance Differences
#' @param stats_df Data frame with statistical results (tissue_comparisons format)
#' @param config Analysis configuration
#' @param save_path Path to save the plot
#' @param fc_threshold Fold change threshold for significance
#' @param p_threshold P-value threshold for significance
#' @param use_adjusted_p Logical: use adjusted p-values instead of raw p-values
#' @param label_significant Logical: label significant points
#' @param max_labels Maximum number of labels to show
#' @param point_size Size of points
#' @param point_alpha Alpha transparency for points
plot_abundance_volcano <- function(stats_df, config, save_path = NULL,
                                   fc_threshold = 1.2, p_threshold = 0.05,
                                   use_adjusted_p = TRUE, label_significant = TRUE,
                                   max_labels = 20, point_size = 3, point_alpha = 0.7) {
  
  # Load required packages
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("ggrepel package is required for volcano plot labels")
  }
  
  # Validate input
  if (is.null(stats_df) || nrow(stats_df) == 0) {
    warning("No data provided for volcano plot")
    return(NULL)
  }
  
  # Required columns
  required_cols <- c("CellType", "fold_change", "p_value")
  missing_cols <- setdiff(required_cols, colnames(stats_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Prepare data for plotting
  plot_data <- stats_df %>%
    filter(!is.na(fold_change) & !is.na(p_value)) %>%
    filter(fold_change > 0 & p_value > 0)  # Remove invalid values
  
  if (nrow(plot_data) == 0) {
    warning("No valid data points for volcano plot")
    return(NULL)
  }
  
  # Choose p-value column
  if (use_adjusted_p && "p_adjusted" %in% colnames(plot_data)) {
    p_col <- "p_adjusted"
    p_label <- "Adjusted P-value"
  } else {
    p_col <- "p_value"
    p_label <- "P-value"
  }
  
  # Calculate log2 fold change and -log10 p-value
  plot_data <- plot_data %>%
    mutate(
      log2_fc = log2(fold_change),
      neg_log10_p = -log10(!!sym(p_col)),
      
      # Classify points based on significance
      change_category = case_when(
        !!sym(p_col) < p_threshold & log2_fc > log2(fc_threshold) ~ "Up-regulated",
        !!sym(p_col) < p_threshold & log2_fc < -log2(fc_threshold) ~ "Down-regulated",
        TRUE ~ "Non-significant"
      ),
      
      # Create labels for significant points
      point_label = case_when(
        !!sym(p_col) < p_threshold & abs(log2_fc) > log2(fc_threshold) ~ CellType,
        TRUE ~ ""
      )
    )
  
  # Limit number of labels if requested
  if (label_significant && max_labels > 0) {
    # Select top significant points by p-value
    top_significant <- plot_data %>%
      filter(point_label != "") %>%
      arrange(!!sym(p_col)) %>%
      head(max_labels)
    
    # Reset labels
    plot_data$point_label <- ifelse(plot_data$CellType %in% top_significant$CellType, 
                                    plot_data$CellType, "")
  }
  
  # Set up colors
  colors <- c(
    "Up-regulated" = "red",
    "Down-regulated" = "blue", 
    "Non-significant" = "grey60"
  )
  
  # Get plot title information
  tissue_name <- if ("Tissue" %in% colnames(plot_data)) unique(plot_data$Tissue)[1] else "Unknown"
  variable_name <- if ("Variable" %in% colnames(plot_data)) unique(plot_data$Variable)[1] else "Unknown"
  clean_var_name <- clean_variable_names(variable_name)
  
  # Create the volcano plot
  p <- ggplot(plot_data, aes(x = log2_fc, y = neg_log10_p, color = change_category)) +
    geom_point(size = point_size, alpha = point_alpha) +
    
    # Add threshold lines
    geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)), 
               linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_hline(yintercept = -log10(p_threshold), 
               linetype = "dashed", color = "black", linewidth = 0.8) +
    
    # Color scale
    scale_color_manual(values = colors, name = "Regulation") +
    
    # Labels and theme
    labs(
      x = "log2(Fold Change)",
      y = paste0("-log10(", p_label, ")"),
      title = paste0("Volcano Plot: ", clean_var_name, " in ", tissue_name, " Tissue"),
      subtitle = paste0("FC threshold: ", fc_threshold, ", P threshold: ", p_threshold)
    ) +
    
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = config$plots$title_size, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = config$plots$text_size, color = "grey40"),
      axis.text = element_text(size = config$plots$text_size - 1),
      axis.title = element_text(size = config$plots$text_size + 1, face = "bold"),
      legend.text = element_text(size = config$plots$text_size - 1),
      legend.title = element_text(size = config$plots$text_size, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # Add labels for significant points
  if (label_significant && sum(plot_data$point_label != "") > 0) {
    p <- p + geom_text_repel(
      aes(label = point_label),
      size = 3,
      box.padding = unit(0.6, "lines"),
      point.padding = unit(0.7, "lines"),
      segment.color = "black",
      segment.size = 0.5,
      show.legend = FALSE,
      max.overlaps = max_labels
    )
  }
  
  # Add summary statistics to the plot
  n_total <- nrow(plot_data)
  n_up <- sum(plot_data$change_category == "Up-regulated")
  n_down <- sum(plot_data$change_category == "Down-regulated")
  n_sig <- n_up + n_down
  
  # Add text annotation with summary
  p <- p + annotate(
    "text",
    x = Inf, y = Inf,
    hjust = 1.1, vjust = 1.1,
    label = paste0(
      "Total: ", n_total, "\n",
      "Significant: ", n_sig, "\n",
      "Up: ", n_up, "\n", 
      "Down: ", n_down
    ),
    size = 3,
    color = "black",
    fontface = "bold"
  )
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(
      filename = save_path,
      plot = p,
      width = config$plots$width,
      height = config$plots$height,
      dpi = config$plots$dpi
    )
    
    if (config$output$verbose) {
      cat("✓ Volcano plot saved:", basename(save_path), "\n")
      cat("  - Total points:", n_total, "\n")
      cat("  - Significant:", n_sig, "(", round(n_sig/n_total*100, 1), "%)\n")
    }
  }
  
  return(p)
}

#' Perform Univariate Cox Regression
#' @param tissue_data Data frame with survival data for specific tissue
#' @param config Analysis configuration
#' @param scale_variables Logical: scale variables before analysis
perform_univariate_cox <- function(tissue_data, config, scale_variables = TRUE) {
  
  # Get survival variable names
  time_var <- config$columns$survival$time
  event_var <- config$columns$survival$event
  
  # Map to actual column names if needed
  if ("RFS_time" %in% colnames(tissue_data)) time_var <- "RFS_time"
  if ("RFS_status" %in% colnames(tissue_data)) event_var <- "RFS_status"
  
  # Validate survival columns
  if (!time_var %in% colnames(tissue_data)) {
    stop("Time variable '", time_var, "' not found in data")
  }
  if (!event_var %in% colnames(tissue_data)) {
    stop("Event variable '", event_var, "' not found in data")
  }
  
  # Get cell type columns
  celltype_cols <- get_celltype_columns(tissue_data, config)
  
  if (length(celltype_cols) == 0) {
    stop("No cell type columns found in data")
  }
  
  # Get clinical feature columns from config
  if (is.null(config$columns)) {
    warning("`config$columns$clinical_features` is not defined. Only analyzing cell types.")
    clinical_cols <- c()
  } else {
    # Ensure the specified clinical columns actually exist in the data
    features_ <- config$columns$continuous_clinical
    features_ <- features_[-1]
    clinical_cols <- intersect(features_, colnames(tissue_data))
  }
  
  # Combine cell types and clinical features into one list of variables to test
  variables_to_test <- unique(c(celltype_cols, clinical_cols))
  
  if (length(variables_to_test) == 0) {
    stop("No cell type or clinical feature columns found to analyze.")
  }
  
  # Remove rows with missing survival data
  survival_complete <- complete.cases(tissue_data[[time_var]], tissue_data[[event_var]])
  tissue_data <- tissue_data[survival_complete, ]
  
  if (nrow(tissue_data) < 10) {
    warning("Insufficient patients with complete survival data (n=", nrow(tissue_data), ")")
    return(NULL)
  }
  
  # Initialize results data frame
  results <- data.frame(
    Variable = character(),
    HR = numeric(),
    Lower_CI = numeric(), 
    Upper_CI = numeric(),
    P_value = numeric(),
    N = integer(),
    Events = integer(),
    stringsAsFactors = FALSE
  )
  
  # Perform univariate Cox regression for each variable (cell type or clinical)
  for (variable in variables_to_test) {
    
    # Get complete cases for this specific predictor variable
    complete_idx <- complete.cases(tissue_data[[variable]], 
                                   tissue_data[[time_var]], 
                                   tissue_data[[event_var]])
    
    if (sum(complete_idx) < 10) {
      if (config$output$verbose) {
        cat("    Skipping", variable, "- insufficient complete observations\n")
      }
      next
    }
    
    analysis_data <- tissue_data[complete_idx, ]
    
    # Scale variable if requested. This works for continuous clinical variables (e.g., Age, CEA)
    # and binary variables. It correctly skips variables with no variance.
    if (scale_variables) {
      if (var(analysis_data[[variable]], na.rm = TRUE) > 1e-6) { # Use a small tolerance for variance
        analysis_data[[variable]] <- scale(analysis_data[[variable]])[, 1]
      } else {
        if (config$output$verbose) {
          cat("    Skipping", variable, "- no variance\n")
        }
        next
      }
    }
    
    tryCatch({
      surv_obj <- Surv(time = analysis_data[[time_var]], event = analysis_data[[event_var]])
      
      # Use the generic 'variable' in the formula
      cox_formula <- as.formula(paste("surv_obj ~", variable))
      cox_model <- coxph(cox_formula, data = analysis_data)
      cox_summary <- summary(cox_model)
      
      results <- rbind(results, data.frame(
        Variable = variable, # The variable name is stored here
        HR = cox_summary$conf.int[1, "exp(coef)"],
        Lower_CI = cox_summary$conf.int[1, "lower .95"],
        Upper_CI = cox_summary$conf.int[1, "upper .95"],
        P_value = cox_summary$coefficients[1, "Pr(>|z|)"],
        N = nrow(analysis_data),
        Events = sum(analysis_data[[event_var]]),
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      if (config$output$verbose) {
        cat("    Error with", variable, ":", e$message, "\n")
      }
    })
  }
  
  if (nrow(results) == 0) {
    warning("No successful Cox regressions were performed.")
    return(NULL)
  }
  
  # Add significance categories
  results$Significance <- case_when(
    results$P_value < 0.001 ~ "***",
    results$P_value < 0.01  ~ "**", 
    results$P_value < 0.05  ~ "*",
    TRUE ~ "ns"
  )
  
  # Add HR confidence interval as text
  results$HR_CI <- paste0(
    format(round(results$HR, 2), nsmall = 2), " (",
    format(round(results$Lower_CI, 2), nsmall = 2), "-",
    format(round(results$Upper_CI, 2), nsmall = 2), ")"
  )
  
  results <- results[order(results$P_value), ]
  rownames(results) <- NULL
  
  if (config$output$verbose) {
    n_sig <- sum(results$P_value < config$stats$alpha_level)
    # Updated log message
    cat("    ✓ Univariate Cox completed:", nrow(results), "variables analyzed,", n_sig, "significant\n")
  }
  
  return(results)
}

#' Perform Multivariate Cox Regression
#' @param tissue_data Data frame with survival data for specific tissue
#' @param config Analysis configuration
#' @param include_clinical Logical: include clinical variables
#' @param max_variables Maximum number of variables to include
perform_multivariate_cox <- function(tissue_data, config, include_clinical = TRUE, 
                                     max_variables = 10) {
  
  # Load required packages
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("survival package is required for Cox regression")
  }
  
  # Get survival variable names
  time_var <- config$columns$survival$time
  event_var <- config$columns$survival$event
  
  # Map to actual column names if needed
  if ("RFS_time" %in% colnames(tissue_data)) time_var <- "RFS_time"
  if ("RFS_status" %in% colnames(tissue_data)) event_var <- "RFS_status"
  
  # First perform univariate analysis to select variables
  univariate_results <- perform_univariate_cox(tissue_data, config, scale_variables = TRUE)
  
  if (is.null(univariate_results) || nrow(univariate_results) == 0) {
    warning("No univariate results available for multivariate analysis")
    return(NULL)
  }
  
  # Select significant variables from univariate analysis
  p_threshold <- 0.1  # Liberal threshold for multivariate inclusion
  significant_vars <- univariate_results$Variable[univariate_results$P_value < p_threshold]
  
  # Limit number of variables to prevent overfitting
  if (length(significant_vars) > max_variables) {
    significant_vars <- significant_vars[1:max_variables]
    if (config$output$verbose) {
      cat("    Limiting to top", max_variables, "variables to prevent overfitting\n")
    }
  }
  
  # Add clinical variables if requested
  variables_to_include <- significant_vars
  
  if (include_clinical) {
    # Add important clinical variables
    clinical_vars <- intersect(c("Age", "Gender", "KRAS_mutation", "TBS", "Differential_grade", "T_stage","CRLM_number","CRLM_size","CEA","CA199"), 
                               colnames(tissue_data))
    variables_to_include <- c(variables_to_include, clinical_vars)
    variables_to_include <- unique(variables_to_include)
  }
  
  if (length(variables_to_include) == 0) {
    warning("No variables selected for multivariate analysis")
    return(NULL)
  }
  
  # Prepare data with complete cases
  analysis_cols <- c(time_var, event_var, variables_to_include)
  complete_data <- tissue_data[complete.cases(tissue_data[, analysis_cols]), ]
  
  if (nrow(complete_data) < 20) {
    warning("Insufficient complete observations for multivariate analysis (n=", 
            nrow(complete_data), ")")
    return(NULL)
  }
  
  # Scale continuous variables
  for (var in variables_to_include) {
    if (is.numeric(complete_data[[var]]) && var %in% get_celltype_columns(complete_data, config)) {
      if (var(complete_data[[var]], na.rm = TRUE) > 0) {
        complete_data[[var]] <- scale(complete_data[[var]])[, 1]
      }
    }
  }
  
  # Build multivariate model
  tryCatch({
    # Create formula
    formula_str <- paste("Surv(", time_var, ",", event_var, ") ~", 
                         paste(variables_to_include, collapse = " + "))
    cox_formula <- as.formula(formula_str)
    
    # Fit multivariate Cox model
    multi_cox <- coxph(cox_formula, data = complete_data)
    
    # Extract results
    cox_summary <- summary(multi_cox)
    
    # Create results data frame
    results <- data.frame(
      Variable = rownames(cox_summary$coefficients),
      HR = cox_summary$conf.int[, "exp(coef)"],
      Lower_CI = cox_summary$conf.int[, "lower .95"],
      Upper_CI = cox_summary$conf.int[, "upper .95"],
      P_value = cox_summary$coefficients[, "Pr(>|z|)"],
      N = nrow(complete_data),
      Events = sum(complete_data[[event_var]]),
      stringsAsFactors = FALSE
    )
    
    # Add significance and confidence intervals
    results$Significance <- case_when(
      results$P_value < 0.001 ~ "***",
      results$P_value < 0.01 ~ "**",
      results$P_value < 0.05 ~ "*", 
      TRUE ~ "ns"
    )
    
    results$HR_CI <- paste0(
      round(results$HR, 2), " (",
      round(results$Lower_CI, 2), "-",
      round(results$Upper_CI, 2), ")"
    )
    
    # Add model statistics
    model_stats <- list(
      model = multi_cox,
      results = results,
      concordance = cox_summary$concordance[1],
      rsquare = cox_summary$rsq[1],
      likelihood_ratio_p = cox_summary$logtest[3],
      n_variables = length(variables_to_include),
      variables_included = variables_to_include
    )
    
    if (config$output$verbose) {
      n_sig <- sum(results$P_value < config$stats$alpha_level)
      cat("    ✓ Multivariate Cox completed:", nrow(results), "variables,", n_sig, "significant\n")
      cat("    ✓ Concordance:", round(model_stats$concordance, 3), "\n")
    }
    
    return(model_stats)
    
  }, error = function(e) {
    warning("Multivariate Cox regression failed: ", e$message)
    return(NULL)
  })
}

#' Generate Forest Plot for Cox Regression Results
#' @param cox_results Data frame with Cox regression results
#' @param tissue_name Name of the tissue
#' @param config Analysis configuration
#' @param save_path Path to save the plot
#' @param show_all Logical: show all variables or only significant ones
#' @param max_variables Maximum number of variables to show
generate_forest_plot <- function(cox_results, tissue_name = "Unknown", config, 
                                 save_path = NULL, show_all = TRUE, max_variables = 20) {
  
  # Validate input
  if (is.null(cox_results) || nrow(cox_results) == 0) {
    warning("No Cox regression results provided for forest plot")
    return(NULL)
  }
  
  # Required columns
  required_cols <- c("Variable", "HR", "Lower_CI", "Upper_CI", "P_value")
  missing_cols <- setdiff(required_cols, colnames(cox_results))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Prepare data
  plot_data <- cox_results
  
  # Filter data if needed
  if (!show_all) {
    plot_data <- plot_data[plot_data$P_value < config$stats$alpha_level, ]
  }
  
  # Limit number of variables
  if (nrow(plot_data) > max_variables) {
    plot_data <- plot_data[1:max_variables, ]
    if (config$output$verbose) {
      cat("  Limiting forest plot to top", max_variables, "variables\n")
    }
  }
  
  if (nrow(plot_data) == 0) {
    warning("No variables to display in forest plot")
    return(NULL)
  }
  
  # Clean variable names for display
  plot_data$Variable_clean <- clean_variable_names(plot_data$Variable)
  
  # Prepare data for forestplot
  # Create matrix for forestplot function
  tabletext <- cbind(
    c("Cell Type", plot_data$Variable_clean),
    c("HR (95% CI)", plot_data$HR_CI),
    c("P-value", sprintf("%.3f", plot_data$P_value))
  )
  
  # Prepare mean (HR) and confidence intervals
  mean_values <- c(NA, plot_data$HR)
  lower_values <- c(NA, plot_data$Lower_CI)
  upper_values <- c(NA, plot_data$Upper_CI)
  
  # Create forest plot
  forest_plot <- forestplot(
    tabletext,
    mean = mean_values,
    lower = lower_values, 
    upper = upper_values,
    zero = 1,
    boxsize = 0.3,
    graph.pos = 2,
    hrzl_lines = list(
      "1" = gpar(lty = 1, lwd = 2),
      "2" = gpar(lty = 2)
    ),
    graphwidth = unit(.25, "npc"),
    xlab = "Hazard Ratio",
    xticks = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
    is.summary = c(TRUE, rep(FALSE, nrow(plot_data))),
    txt_gp = fpTxtGp(
      label = gpar(cex = 1),
      ticks = gpar(cex = 1),
      xlab = gpar(cex = 1.2),
      title = gpar(cex = 1.4)
    ),
    lwd.zero = 2,
    lwd.ci = 1.5,
    lwd.xaxis = 2,
    lty.ci = 1,
    ci.vertices = TRUE,
    ci.vertices.height = 0.15,
    clip = c(0.05, 10),
    title = paste0("Forest Plot: ", tissue_name, " Tissue\nUnivariate Cox Regression"),
    col = fpColors(
      box = ifelse(c(FALSE, plot_data$P_value < config$stats$alpha_level), "red", "blue"),
      lines = ifelse(c(FALSE, plot_data$P_value < config$stats$alpha_level), "red", "blue"),
      zero = "black"
    )
  )
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    pdf(save_path, width = config$plots$width, height = max(6, nrow(plot_data) * 0.4))
    print(forest_plot)
    dev.off()
    
    if (config$output$verbose) {
      cat("✓ Forest plot saved:", basename(save_path), "\n")
      n_sig <- sum(plot_data$P_value < config$stats$alpha_level)
      cat("  - Variables shown:", nrow(plot_data), "\n")
      cat("  - Significant:", n_sig, "\n")
    }
  }
  
  return(forest_plot)
}

#' Generate Kaplan-Meier Plot for Single Variable
#' @param survival_df Data frame with survival data
#' @param variable Variable name to analyze
#' @param config Analysis configuration
#' @param cutoff_method Method to dichotomize variable ("median", "optimal", or numeric value)
prepare_kaplan_meier_plot <- function(survival_df, variable, config,
                                      cutoff_method = "optimal") {
  
  # Get survival variable names
  time_var <- config$columns$survival$time
  event_var <- config$columns$survival$event
  
  # Map to actual column names if needed
  if ("RFS_time" %in% colnames(survival_df)) time_var <- "RFS_time"
  if ("RFS_status" %in% colnames(survival_df)) event_var <- "RFS_status"
  
  # Validate inputs
  if (!variable %in% colnames(survival_df)) {
    stop("Variable '", variable, "' not found in survival data")
  }
  if (!time_var %in% colnames(survival_df)) {
    stop("Time variable '", time_var, "' not found in survival data")
  }
  if (!event_var %in% colnames(survival_df)) {
    stop("Event variable '", event_var, "' not found in survival data")
  }
  
  # Prepare data - safer column selection and missing value handling
  required_cols <- c(variable, time_var, event_var)
  col_indices <- match(required_cols, colnames(survival_df))
  
  # Check if all columns were found
  if (any(is.na(col_indices))) {
    missing_cols <- required_cols[is.na(col_indices)]
    stop("Required columns not found: ", paste(missing_cols, collapse = ", "))
  }
  
  # Select columns safely
  analysis_data <- survival_df[, col_indices, drop = FALSE]
  colnames(analysis_data) <- c("variable", "time", "event")
  
  # Remove rows with any missing values
  complete_rows <- complete.cases(analysis_data)
  analysis_data <- analysis_data[complete_rows, , drop = FALSE]
  
  if (nrow(analysis_data) < 10) {
    warning("Insufficient data for Kaplan-Meier plot (n=", nrow(analysis_data), ")")
    return(NULL)
  }
  
  # Dichotomize the variable
  if (cutoff_method == "median") {
    cutoff <- median(analysis_data$variable, na.rm = TRUE)
  } else if (cutoff_method == "optimal") {
    # Use survminer's surv_cutpoint for optimal cutoff with proper column names
    tryCatch({
      cutpoint_result <- surv_cutpoint(
        data = analysis_data,
        time = "time",        # Use renamed column
        event = "event",      # Use renamed column
        variables = "variable" # Use renamed column
      )
      cutoff <- summary(cutpoint_result)$cutpoint
    }, error = function(e) {
      warning("Optimal cutpoint failed for ", variable, ": ", e$message, 
              ". Using median instead.")
      cutoff <- median(analysis_data$variable, na.rm = TRUE)
    })
  } else if (is.numeric(cutoff_method)) {
    cutoff <- cutoff_method
  } else {
    cutoff <- median(analysis_data$variable, na.rm = TRUE)
  }
  
  # Create high/low groups
  analysis_data$Group <- ifelse(analysis_data$variable >= cutoff, "High", "Low")
  analysis_data$Group <- factor(analysis_data$Group, levels = c("Low", "High"))
  analysis_data$event <- as.numeric(as.character(analysis_data$event))
  
  # Check group sizes
  group_sizes <- table(analysis_data$Group)
  if (any(group_sizes < 5)) {
    warning("Small group sizes for ", variable, ": ", paste(group_sizes, collapse = ", "))
  }

  return(analysis_data)
}

# Initialize package loading when script is sourced
if (!exists("packages_loaded")) {
  load_required_packages()
  packages_loaded <- TRUE
  cat("✓ abundance_functions.R loaded successfully\n")
  cat("⚠️  Note: Some helper functions need implementation (see section 8)\n")
}
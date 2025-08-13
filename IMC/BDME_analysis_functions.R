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

library(pheatmap)
library(ComplexHeatmap)
library(forestplot)
library(corrplot)
library(RColorBrewer)
library(viridis)
library(patchwork)

library(pROC)
library(survival)
library(survminer)
library(circlize)
library(ConsensusClusterPlus)

library(dplyr)
library(tidyr)
library(purrr)

# ===============================================================================
# Calculate Cholangiocyte Zone Sizes
# ===============================================================================

calculate_zone_sizes <- function(meta, coords) {
  
  # Filter for valid zones
  valid_zones <- meta %>%
    filter(!is.na(Tumor_patch) & Tumor_patch != "NA")
  
  # Get unique combinations of patient, sample, and zone
  unique_zones <- valid_zones %>%
    select(patient_id, sample_id, Tumor_patch, Treatment, RFS_status, Treatment_RFS) %>%
    distinct()
  
  # Calculate total cells per sample for fraction calculation
  cells_per_sample <- meta %>%
    group_by(sample_id) %>%
    summarise(total_sample_cells = n(), .groups = "drop")
  
  # Initialize results dataframe
  zone_sizes <- data.frame()
  
  # Use for-loop to calculate zone sizes
  for(i in 1:nrow(unique_zones)) {
    
    current_zone <- unique_zones[i, ]
    
    # Get cells in this specific zone
    idx <- meta$Tumor_patch %in% current_zone$Tumor_patch
    zone_cells <- meta[idx,] %>%
      filter(patient_id == current_zone$patient_id,
             sample_id == current_zone$sample_id,
             Tumor_patch == current_zone$Tumor_patch)
    
    zone_cell_count <- nrow(zone_cells)
    
    # Calculate zone area if we have enough cells
    if(zone_cell_count >= 3) {
      zone_coords <- coords[idx, ]
      
      if(nrow(zone_coords) >= 3) {
        # Calculate area using bounding box approach
        x_range <- diff(range(zone_coords[,1]))
        y_range <- diff(range(zone_coords[,2]))
        zone_area <- x_range * y_range
      } else {
        zone_area <- NA
      }
    } else {
      zone_area <- NA
    }
    
    # Get total cells in this sample for fraction calculation
    total_sample_cells <- cells_per_sample$total_sample_cells[
      cells_per_sample$sample_id == current_zone$sample_id
    ]
    
    # Calculate zone cell fraction
    zone_cell_fraction <- zone_cell_count / total_sample_cells
    
    # Store results
    zone_result <- data.frame(
      patient_id = current_zone$patient_id,
      sample_id = current_zone$sample_id,
      Tumor_patch = current_zone$Tumor_patch,
      Treatment = current_zone$Treatment,
      RFS_status = current_zone$RFS_status,
      Treatment_RFS = current_zone$Treatment_RFS,
      zone_cell_count = zone_cell_count,
      zone_area = zone_area,
      zone_cell_fraction = zone_cell_fraction,
      total_sample_cells = total_sample_cells
    )
    
    zone_sizes <- rbind(zone_sizes, zone_result)
    
    # Progress indicator
    if(i %% 50 == 0 || i == nrow(unique_zones)) {
      cat("Processed", i, "of", nrow(unique_zones), "zones\n")
    }
  }
  
  # Filter out zones with invalid areas
  valid_zone_sizes <- zone_sizes %>%
    filter(!is.na(zone_area))
  
  cat("Zone sizes calculated for", nrow(valid_zone_sizes), "zones (from", nrow(zone_sizes), "total zones)\n")
  cat("Zone cell fraction range:", round(range(zone_sizes$zone_cell_fraction, na.rm = TRUE), 4), "\n")
  
  return(zone_sizes)  # Return all zones, including those with NA area
}

# ===============================================================================
# Create Treatment-specific Comparisons (R vs NR)
# ===============================================================================

create_comparison_groups <- function(zone_data) {
  
  # Ensure we have the RFS grouping
  zone_data$RFS_group <- ifelse(zone_data$RFS_status == 0, "NR", "R")
  
  comparisons <- list(
    # All patients: R vs NR
    all = zone_data,
    
    # Chemotherapy only: R vs NR  
    chemo = zone_data[zone_data$Treatment == "Chemo", ],
    
    # Combination therapy: R vs NR
    combo = zone_data[zone_data$Treatment == "Combo", ]
  )
  
  # Print sample sizes
  for(comp_name in names(comparisons)) {
    comp_data <- comparisons[[comp_name]]
    cat("\n", toupper(comp_name), "comparison:\n")
    print(table(comp_data$RFS_group))
  }
  
  return(comparisons)
}

# ===============================================================================
# Function to create boxplot with wilcox test
# ===============================================================================

create_combined_size_boxplots <- function(data, title, filename_prefix) {
  
  if(nrow(data) == 0 || length(unique(data$RFS_group)) < 2) {
    cat("Insufficient data for", title, "\n")
    return(NULL)
  }
  
  plots_list <- list()
  results_list <- list()
  
  # 1. Zone Area Plot (only for zones with valid area)
  area_data <- data %>% filter(!is.na(zone_area))
  
  if(nrow(area_data) > 5 && length(unique(area_data$RFS_group)) == 2) {
    wilcox_area <- wilcox.test(zone_area ~ RFS_group, data = area_data)
    p_value_area <- wilcox_area$p.value
    
    p_area <- ggplot(area_data, aes(x = RFS_group, y = zone_area, fill = RFS_group)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.6, color = "grey75") +
      scale_fill_manual(values = c("NR" = "#0073C2FF", "R" = "#EFC000FF")) +
      labs(
        title = "Zone Area",
        subtitle = paste("p =", format(p_value_area, digits = 3)),
        x = "Recurrence Status",
        y = "Zone Area (um²)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none"
      )
    
    plots_list[["area"]] <- p_area
    results_list[["area"]] <- list(p_value = p_value_area, n_zones = nrow(area_data))
  }
  
  # 2. Zone Cell Fraction Plot
  wilcox_fraction <- wilcox.test(zone_cell_fraction ~ RFS_group, data = data)
  p_value_fraction <- wilcox_fraction$p.value
  
  p_fraction <- ggplot(data, aes(x = RFS_group, y = zone_cell_fraction, fill = RFS_group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, color = "grey75") +
    scale_fill_manual(values = c("NR" = "#0073C2FF", "R" = "#EFC000FF")) +
    labs(
      title = "Zone Cell Fraction",
      subtitle = paste("p =", format(p_value_fraction, digits = 3)),
      x = "Recurrence Status",
      y = "Zone Cell Fraction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  plots_list[["fraction"]] <- p_fraction
  results_list[["fraction"]] <- list(p_value = p_value_fraction, n_zones = nrow(data))
  
  # 3. Zone Cell Count Plot
  wilcox_count <- wilcox.test(zone_cell_count ~ RFS_group, data = data)
  p_value_count <- wilcox_count$p.value
  
  p_count <- ggplot(data, aes(x = RFS_group, y = zone_cell_count, fill = RFS_group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, color = "grey75") +
    scale_fill_manual(values = c("NR" = "#0073C2FF", "R" = "#EFC000FF")) +
    labs(
      title = "Zone Cell Count",
      subtitle = paste("p =", format(p_value_count, digits = 3)),
      x = "Recurrence Status", 
      y = "Zone Cell Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  plots_list[["count"]] <- p_count
  results_list[["count"]] <- list(p_value = p_value_count, n_zones = nrow(data))
  
  # Combine plots using patchwork
  if(length(plots_list) == 3) {
    combined_plot <- plots_list[["area"]] + plots_list[["fraction"]] + plots_list[["count"]] +
      plot_annotation(
        title = paste("Cholangiocyte Zone Analysis -", title),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
      ) +
      plot_layout(ncol = 3)
  } else {
    # If no area data, just combine fraction and count
    combined_plot <- plots_list[["fraction"]] + plots_list[["count"]] +
      plot_annotation(
        title = paste("Cholangiocyte Zone Analysis -", title),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
      ) +
      plot_layout(ncol = 2)
  }
  
  # Save combined plot
  combined_filename <- paste0(filename_prefix, "_combined.pdf")
  ggsave(combined_filename, combined_plot, width = 15, height = 5, dpi = 300)
  
  return(list(
    combined_plot = combined_plot,
    individual_plots = plots_list,
    statistics = results_list
  ))
}

# ===============================================================================
# Function to calculate composition differences and stats
# ===============================================================================

calculate_composition_stats <- function(data, condition_name) {
  if(nrow(data) == 0) return(NULL)
  
  # Get proportion columns
  prop_cols <-  colnames(data)[endsWith(colnames(data),suffix = "_prop")]
  
  # Calculate differences and stats for each cell type
  results <- prop_cols %>%
    map_dfr(~ {
      col_name <- .x
      cell_name <- gsub("_prop$", "", .x)
      
      if(length(unique(data$RFS_group)) == 2) {
        # Calculate mean difference (R - NR)
        means <- data %>%
          group_by(RFS_group) %>%
          summarise(mean_prop = mean(.data[[col_name]], na.rm = TRUE), .groups = "drop")
        
        r_mean <- means$mean_prop[means$RFS_group == "R"]
        nr_mean <- means$mean_prop[means$RFS_group == "NR"]
        difference = (r_mean + 1e-6) / (nr_mean + 1e-6)
        
        # Statistical test
        test_result <- wilcox.test(data[[col_name]] ~ data$RFS_group)
        
        data.frame(
          condition = condition_name,
          cell_type = cell_name,
          difference = log2(difference),
          p_value = test_result$p.value,
          r_mean = r_mean,
          nr_mean = nr_mean
        )
      } else {
        NULL
      }
    })
  
  return(results)
}

calculate_expression_stats <- function(data, condition_name) {
  if(nrow(data) == 0) return(NULL)
  
  # Get expression marker columns (exclude metadata columns)
  metadata_cols <- c("patient_id", "sample_id", "Tumor_patch", "Treatment", 
                     "RFS_status", "Treatment_RFS", "RFS_group")
  expr_cols <- setdiff(colnames(data), metadata_cols)
  
  # Calculate differences and stats for each marker
  results <- expr_cols %>%
    map_dfr(~ {
      col_name <- .x
      marker_name <- .x
      
      if(length(unique(data$RFS_group)) == 2) {
        # Calculate mean expression (R vs NR)
        means <- data %>%
          group_by(RFS_group) %>%
          summarise(mean_expr = mean(as.numeric(.data[[col_name]]) , na.rm = TRUE), .groups = "drop")
        
        r_mean <- means$mean_expr[means$RFS_group == "R"]
        nr_mean <- means$mean_expr[means$RFS_group == "NR"]
        
        # Calculate log2 fold change (add small value to avoid log(0))
        log2_fc <- log2((r_mean + 1e-6) / (nr_mean + 1e-6))
        # Or simple difference if preferred
        difference <- r_mean - nr_mean
        
        # Statistical test
        test_result <- wilcox.test(as.numeric( data[[col_name]]) ~ data$RFS_group)
        
        data.frame(
          condition = condition_name,
          marker = marker_name,
          log2_fc = log2_fc,
          difference = difference,
          p_value = test_result$p.value,
          r_mean = r_mean,
          nr_mean = nr_mean
        )
      } else {
        NULL
      }
    })
  
  return(results)
}

# ===============================================================================
# Calculate the Expression of cholangiocytes within different PCME
# ===============================================================================
calculate_ec_mean_expression <- function(pt_cells,expr_data,markers){
  cat("\nStep 4.1: Calculating cholangiocytes mean expression ...\n")
  
  # Get CAIX expression data for epithelial cells
  # Check if CAIX is available
  if(!"CAIX" %in% rownames(expr_data)) {
    available_genes <- grep("CAIX|CA9|CA_IX", rownames(expr_data), value = TRUE, ignore.case = TRUE)
    cat("CAIX not found. Available related genes:", paste(available_genes, collapse = ", "), "\n")
    return(NULL)
  }
  
  # Get unique tumor patches
  tumor_patches <- pt_cells %>%
    filter(!is.na(Tumor_patch) & Tumor_patch != "NA") %>%
    distinct(patient_id, sample_id, Tumor_patch, Treatment, RFS_status, Treatment_RFS)
  
  # Initialize results list
  results_list <- list()
  
  # Loop through each tumor patch
  for(i in 1:nrow(tumor_patches)) {
    patch_info <- tumor_patches[i, ]
    
    # Get cells in this patch
    patch_cells <- pt_cells %>%
      filter(
        sub_celltype == "EC", 
        patient_id == patch_info$patient_id,
        sample_id == patch_info$sample_id,
        Tumor_patch == patch_info$Tumor_patch
      )
    
    # Get cell IDs and check which exist in expr_data
    cell_ids <- rownames(patch_cells)
    valid_cells <- cell_ids[cell_ids %in% colnames(expr_data)]
    
    # Calculate patch size and mean expressions
    if(length(valid_cells) > 0) {
      # Get expression data for valid cells
      patch_expr <- expr_data[markers, valid_cells, drop = FALSE]
      
      # Calculate means for each marker
      marker_means <- rowMeans(patch_expr, na.rm = TRUE)
      
      # Create result row
      result_row <- data.frame(
        patient_id = patch_info$patient_id,
        sample_id = patch_info$sample_id,
        Tumor_patch = patch_info$Tumor_patch,
        Treatment = patch_info$Treatment,
        RFS_status = patch_info$RFS_status,
        Treatment_RFS = patch_info$Treatment_RFS,
        patch_size = length(valid_cells),
        t(marker_means)  # transpose to make markers as columns
      )
    } else {
      # No valid cells
      result_row <- data.frame(
        patient_id = patch_info$patient_id,
        sample_id = patch_info$sample_id,
        Tumor_patch = patch_info$Tumor_patch,
        Treatment = patch_info$Treatment,
        RFS_status = patch_info$RFS_status,
        Treatment_RFS = patch_info$Treatment_RFS,
        patch_size = 0,
        matrix(NA, nrow = 1, ncol = length(markers), dimnames = list(NULL, markers))
      )
    }
    
    results_list[[i]] <- result_row
  }
  
  # Combine all results
  epithelial_caix_mean <- do.call(rbind, results_list)
  rownames(epithelial_caix_mean) <- NULL
  
  return(epithelial_caix_mean)
}

# ===============================================================================
# CAIX-only Classification (Simple)
# ===============================================================================

classify_by_caix_only <- function(zone_composition, epithelial_exp_mean) {
  
  # Extract CAIX expression for each zone
  caix_data <- epithelial_exp_mean %>%
    select(patient_id, sample_id, Tumor_patch, Treatment, RFS_status, Treatment_RFS) %>%
    mutate(
      caix_expression = as.numeric(epithelial_exp_mean$CA_IX)
    )
  
  # Define CAIX threshold (median split or quantile-based)
  caix_threshold <- median(caix_data$caix_expression, na.rm = TRUE)
  
  pcme_classification <- caix_data %>%
    mutate(
      PCME_type = ifelse(caix_expression >= caix_threshold, 
                         "Stromal-fibrotic", "Immune-permissive"),
      classification_method = "CAIX-only",
      caix_threshold = caix_threshold
    )
  cat("CAIX threshold:", round(caix_threshold, 3), "\n")
  cat("PCME distribution:\n")
  print(table(pcme_classification$PCME_type))
  pcme_classification$PCME_binary <- pcme_classification$PCME_type
  
  return(pcme_classification)
}

# ===============================================================================
# Dual Characteristics Classification (Comprehensive)
# ===============================================================================

classify_by_dual_characteristics <- function(zone_composition, epithelial_exp_mean) {
  
  # Step 1: Extract cholangiocyte characteristics
  cholangiocyte_features <- epithelial_exp_mean %>%
    select(patient_id, sample_id, Tumor_patch, Treatment, RFS_status, Treatment_RFS) %>%
    mutate(
      caix_expression = as.numeric(epithelial_exp_mean$CA_IX)
    )
  
  # Step 2: Extract microenvironment characteristics
  microenv_features <- zone_composition %>%
    select(patient_id, sample_id, Tumor_patch, Treatment, RFS_status, Treatment_RFS) %>%
    mutate(
      # Immune signature (CD8T, CD4T, Macrophages)
      immune_score = rowSums(select(zone_composition, matches("CD8T|CD4T|Macro_HLADR|Macro_CD11b")), na.rm = TRUE),
      
      # Stromal signature (Collagen, Vimentin-positive cells)
      stromal_score = rowSums(select(zone_composition, matches("SC_.*Collagen|SC_.*Vimentin|CAF")), na.rm = TRUE),
      
      # Calculate relative dominance
      immune_stromal_ratio = ifelse(stromal_score > 0, immune_score / (stromal_score), 1)
      # immune_stromal_ratio = ifelse(stromal_score > 0, immune_score / (immune_score + stromal_score), 1)
    )
  
  # Step 3: Define thresholds
  caix_threshold <- median(cholangiocyte_features$caix_expression, na.rm = TRUE)
  immune_threshold <- median(microenv_features$immune_stromal_ratio, na.rm = TRUE)
  
  # Step 4: Combine features and classify
  pcme_classification <- cholangiocyte_features %>%
    left_join(microenv_features, by = c("patient_id", "sample_id", "Tumor_patch", 
                                        "Treatment", "RFS_status", "Treatment_RFS")) %>%
    mutate(
      # Cholangiocyte state
      cholangiocyte_state = ifelse(caix_expression >= caix_threshold, "Abnormal", "Normal"),
      
      # Microenvironment state  
      microenv_state = ifelse(immune_stromal_ratio >= immune_threshold, "Immune-dominant", "Stromal-dominant"),
      
      # Combined PCME classification
      PCME_type = case_when(
        cholangiocyte_state == "Normal" & microenv_state == "Immune-dominant" ~ "Immune-permissive",
        cholangiocyte_state == "Abnormal" & microenv_state == "Stromal-dominant" ~ "Stromal-fibrotic", 
        cholangiocyte_state == "Normal" & microenv_state == "Stromal-dominant" ~ "Transitional-stromal",
        cholangiocyte_state == "Abnormal" & microenv_state == "Immune-dominant" ~ "Transitional-immune",
        TRUE ~ "Undefined"
      ),
      
      # Simplified binary classification (focus on clear types)
      PCME_binary = case_when(
        PCME_type == "Immune-permissive" ~ "Immune-permissive",
        PCME_type == "Stromal-fibrotic" ~ "Stromal-fibrotic",
        TRUE ~ "Mixed/Transitional"
      ),
      
      classification_method = "Dual-characteristics",
      caix_threshold = caix_threshold,
      immune_threshold = immune_threshold
    )
  
  # Print classification summary
  cat("Classification thresholds:\n")
  cat("- CAIX expression:", round(caix_threshold, 3), "\n")
  cat("- Immune/Stromal ratio:", round(immune_threshold, 3), "\n\n")
  
  cat("PCME type distribution (detailed):\n")
  print(table(pcme_classification$PCME_type))
  
  cat("\nPCME type distribution (binary):\n")
  print(table(pcme_classification$PCME_binary))
  
  cat("\nPCME type and treatment response:\n")
  print(table(pcme_classification$Treatment_RFS ,pcme_classification$PCME_binary))
  
  return(pcme_classification)
}

# ===============================================================================
# Step 2: Validate PCME Classification
# ===============================================================================

validate_pcme_classification <- function(pcme_classification) {
  
  # Test association with clinical outcomes
  survival_analysis <- pcme_classification %>%
    group_by(PCME_type, Treatment) %>%
    summarise(
      n_zones = n(),
      recurrence_rate = mean(RFS_status, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Treatment, desc(recurrence_rate))
  
  print("PCME types vs Recurrence rates:")
  print(survival_analysis)
  
  # Statistical testing
  if("PCME_binary" %in% colnames(pcme_classification)) {
    
    # Test for treatment-specific associations
    chemo_test <- pcme_classification %>%
      filter(Treatment == "Chemo" & PCME_binary %in% c("Immune-permissive", "Stromal-fibrotic"))
    
    combo_test <- pcme_classification %>%
      filter(Treatment == "Combo" & PCME_binary %in% c("Immune-permissive", "Stromal-fibrotic"))
    
    if(nrow(chemo_test) > 10) {
      chemo_fisher <- fisher.test(table(chemo_test$PCME_binary, chemo_test$RFS_status))
      cat("\nChemotherapy - PCME vs Recurrence (Fisher's test):\n")
      cat("p-value:", format(chemo_fisher$p.value, digits = 4), "\n")
    }
    
    if(nrow(combo_test) > 10) {
      combo_fisher <- fisher.test(table(combo_test$PCME_binary, combo_test$RFS_status))
      cat("\nCombination therapy - PCME vs Recurrence (Fisher's test):\n")
      cat("p-value:", format(combo_fisher$p.value, digits = 4), "\n")
    }
  }
  
  return(NULL)
}

# ===============================================================================
# Step 3: Visualize PCME Classification
# ===============================================================================

visualize_pcme_classification <- function(pcme_classification, zone_composition, 
                                          output_dir = ".") {
  
  # 1. PCME distribution by treatment and outcome
  p1 <- ggplot(pcme_classification, aes(x = Treatment_RFS, fill = PCME_binary)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    labs(
      title = "PCME Distribution by Treatment and Outcome",
      x = "Treatment & Recurrence Status",
      y = "Proportion of Zones",
      fill = "PCME Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(file.path(output_dir, "pcme_distribution_by_treatment.png"), 
         p1, width = 10, height = 6, dpi = 300)
  
  # 2. Feature space visualization (if dual classification)
  if("cholangiocyte_state" %in% colnames(pcme_classification)) {
    
    p2 <- ggplot(pcme_classification, aes(x = caix_expression, y = immune_stromal_ratio)) +
      geom_point(aes(color = PCME_binary, shape = Treatment_RFS), size = 3, alpha = 0.7) +
      geom_hline(yintercept = unique(pcme_classification$immune_threshold), 
                 linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = unique(pcme_classification$caix_threshold), 
                 linetype = "dashed", color = "gray50") +
      scale_color_brewer(type = "qual", palette = "Set1") +
      labs(
        title = "PCME Classification in Feature Space",
        x = "CAIX Expression Level",
        y = "Immune/Stromal Ratio",
        color = "PCME Type",
        shape = "Treatment & Outcome"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave(file.path(output_dir, "pcme_feature_space.png"), 
           p2, width = 12, height = 8, dpi = 300)
  }
  
  # 3. Heatmap of cellular composition by PCME type
  if("PCME_binary" %in% colnames(pcme_classification)) {
    
    # Calculate mean proportions by PCME type
    cell_props <- zone_composition %>%
      select(patient_id, sample_id, Tumor_patch, ends_with("_prop")) %>%
      left_join(pcme_classification %>% select(patient_id, sample_id, Tumor_patch, PCME_binary),
                by = c("patient_id", "sample_id", "Tumor_patch")) %>%
      filter(PCME_binary %in% c("Immune-permissive", "Stromal-fibrotic")) %>%
      select(PCME_binary, ends_with("_prop")) %>%
      group_by(PCME_binary) %>%
      summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
      tibble::column_to_rownames("PCME_binary") %>%
      as.matrix()
    
    # Create heatmap
    if(ncol(cell_props) > 0 && nrow(cell_props) > 1) {
      
      # Color palette
      col_fun <- colorRamp2(c(0,0.15,1), c("white", "red", "darkred"))
      
      pdf(file.path(output_dir, "pcme_cellular_composition_heatmap.pdf"), width = 12, height = 6)
      
      ht <- Heatmap(
        cell_props,
        name = "Proportion",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        row_names_side = "left",
        column_names_rot = 90,
        heatmap_legend_param = list(title = "Cell\nProportion"),
        row_title = "PCME Type",
        column_title = "Cell Types"
      )
      
      draw(ht)
      dev.off()
    }
  }
  
  return(list(distribution_plot = p1, feature_plot = p2))
}

# ===============================================================================
# Aggregate PCME Zones to Patient-Level Scores
# ===============================================================================

aggregate_pcme_to_patient <- function(pcme_classification, zone_sizes = NULL) {
  
  cat("=== Aggregating PCME to Patient Level ===\n")
  
  # Get patient-level clinical features (assuming they're consistent within patient)
  patient_clinical <- pcme_classification %>%
    group_by(patient_id) %>%
    summarise(
      Treatment = first(Treatment),
      RFS_status = first(RFS_status),
      .groups = "drop"
    )
  
  # Method 1: Simple zone counting approach
  patient_pcme_counts <- pcme_classification %>%
    filter(PCME_binary %in% c("Immune-permissive", "Stromal-fibrotic")) %>%
    group_by(patient_id) %>%
    summarise(
      total_zones = n(),
      PCME_I_zones = sum(PCME_binary == "Immune-permissive"),
      PCME_S_zones = sum(PCME_binary == "Stromal-fibrotic"),
      
      # Calculate ratios and proportions
      PCME_S_ratio = PCME_S_zones / (PCME_I_zones + PCME_S_zones),
      PCME_I_ratio = PCME_I_zones / (PCME_I_zones + PCME_S_zones),
      
      # Log ratio (better for statistical analysis)
      log_ratio = log((PCME_S_zones + 0.5) / (PCME_I_zones + 0.5)),
      
      # Dominant pattern
      dominant_PCME = ifelse(PCME_S_zones > PCME_I_zones, "PCME-S", "PCME-I"),
      
      .groups = "drop"
    )
  
  # Method 2: If zone sizes are available, use weighted approach
  if(!is.null(zone_sizes)) {
    
    # Merge with zone sizes
    pcme_with_sizes <- pcme_classification %>%
      left_join(zone_sizes %>% select(patient_id, sample_id, Tumor_patch, 
                                      zone_cell_count, zone_area, zone_cell_fraction),
                by = c("patient_id", "sample_id", "Tumor_patch")) %>%
      filter(PCME_binary %in% c("Immune-permissive", "Stromal-fibrotic"))
    
    # Weighted by cell count
    patient_pcme_weighted <- pcme_with_sizes %>%
      group_by(patient_id) %>%
      summarise(
        total_cells = sum(zone_cell_count, na.rm = TRUE),
        PCME_I_cells = sum(ifelse(PCME_binary == "Immune-permissive", zone_cell_count, 0), na.rm = TRUE),
        PCME_S_cells = sum(ifelse(PCME_binary == "Stromal-fibrotic", zone_cell_count, 0), na.rm = TRUE),
        
        # Cell-weighted ratios
        PCME_S_ratio_weighted = PCME_S_cells / (PCME_I_cells + PCME_S_cells),
        log_ratio_weighted = log((PCME_S_cells + 0.5) / (PCME_I_cells + 0.5)),
        
        .groups = "drop"
      )
    
    # Merge weighted scores
    patient_pcme_counts <- patient_pcme_counts %>%
      left_join(patient_pcme_weighted, by = "patient_id")
  }
  
  # Combine with clinical data
  patient_scores <- patient_clinical %>%
    left_join(patient_pcme_counts, by = "patient_id") %>%
    filter(!is.na(total_zones))
  
  cat("Patient-level aggregation completed:\n")
  cat("- Total patients with PCME data:", nrow(patient_scores), "\n")
  cat("- Patients by treatment:", table(patient_scores$Treatment), "\n")
  cat("- PCME-S ratio range:", round(range(patient_scores$PCME_S_ratio, na.rm = TRUE), 3), "\n")
  
  return(patient_scores)
}

# ===============================================================================
# Find Optimal Threshold for PCME Score
# ===============================================================================

find_optimal_threshold <- function(patient_data, score_var = "PCME_S_ratio", 
                                   outcome_var = "RFS_status", time_var = "RFS_time", 
                                   method = "youden") {
  
  cat("\n=== Finding Optimal Threshold ===\n")
  cat("Score variable:", score_var, "\n")
  cat("Method:", method, "\n")
  
  # Remove missing values
  if(method == "youden") {
    analysis_data <- patient_data %>%
      filter(!is.na(.data[[score_var]]) & !is.na(.data[[outcome_var]]))
  } else {
    analysis_data <- patient_data %>%
      filter(!is.na(.data[[score_var]]) & !is.na(.data[[outcome_var]]) & !is.na(.data[[time_var]]))
  }
  
  if(nrow(analysis_data) < 10) {
    cat("Insufficient data for threshold analysis\n")
    return(NULL)
  }
  
  if(method == "youden") {
    # Use ROC-based Youden index for AUC optimization
    roc_obj <- roc(analysis_data[[outcome_var]], analysis_data[[score_var]], 
                   direction = "<", quiet = TRUE)
    
    threshold <- coords(roc_obj, "best", ret = "threshold", 
                        best.method = "youden")[1,1]
    
    auc_value <- auc(roc_obj)
    surv_pvalue <- NULL
    
  } else if(method == "best") {
    # Use survminer surv_cutpoint for survival-based optimization
    
    # Prepare data for surv_cutpoint (needs specific column names)
    surv_data <- analysis_data %>%
      select(all_of(c(score_var, outcome_var, time_var))) %>%
      rename(
        score = !!score_var,
        status = !!outcome_var, 
        time = !!time_var
      )
    
    # Find optimal cutpoint using survminer
    cutpoint_result <- surv_cutpoint(
      data = surv_data,
      time = "time",
      event = "status", 
      variables = "score"
    )
    
    threshold <- cutpoint_result$cutpoint$cutpoint[1]
    surv_pvalue <- cutpoint_result$cutpoint$pvalue[1]
    auc_value <- NULL
    
  } 
else if(method == "median") {
  threshold <- 0 # median(analysis_data[[score_var]])
  auc_value <- NULL
  surv_pvalue <- NULL
  
}
  else {
    stop("Method must be either 'youden' or 'survminer'")
  }
  
  
  
  # Create risk groups
  analysis_data$risk_group <- ifelse(analysis_data[[score_var]] >= threshold, 
                                     "High-risk", "Low-risk")
  
  cat("Optimal threshold:", round(threshold, 4), "\n")
  cat("Risk group distribution:\n")
  print(table(analysis_data$risk_group))
  
  if(!is.null(auc_value)) {
    cat("AUC:", round(auc_value, 3), "\n")
  }
  
  if(!is.null(surv_pvalue)) {
    cat("Survival p-value:", format(surv_pvalue, digits = 4), "\n")
  }
  
  return(list(
    threshold = threshold,
    risk_groups = analysis_data$risk_group,
    auc = auc_value,
    surv_pvalue = surv_pvalue,
    analysis_data = analysis_data,
    method = method
  ))
}


# ===============================================================================
# Draw forest plot
# ===============================================================================

generate_forest_plot <- function(cox_results, tissue_name = "Unknown", 
                                 save_path = NULL, show_all = TRUE, 
                                 max_variables = 20, alpha_level = 0.05) {
  
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
    plot_data <- plot_data[plot_data$P_value < alpha_level, ]
  }
  
  # Limit number of variables and sort by p-value
  if (nrow(plot_data) > max_variables) {
    plot_data <- plot_data[order(plot_data$P_value), ]
    plot_data <- plot_data[1:max_variables, ]
    cat("Limiting forest plot to top", max_variables, "variables\n")
  }
  
  if (nrow(plot_data) == 0) {
    warning("No variables to display in forest plot")
    return(NULL)
  }
  
  # Clean variable names
  clean_variable_names <- function(var_names) {
    var_names %>%
      gsub("_prop$", "", .) %>%
      gsub("_", " ", .) %>%
      stringr::str_to_title(.)
  }
  
  # Prepare plot data
  plot_data <- plot_data %>%
    mutate(
      Variable_clean = clean_variable_names(Variable),
      Significant = P_value < alpha_level,
      HR_CI = paste0(round(HR, 2), " (", round(Lower_CI, 2), "-", round(Upper_CI, 2), ")"),
      P_text = case_when(
        P_value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", P_value)
      )
    ) %>%
    arrange(desc(HR))
  
  # Debug: Print significance status
  cat("\n=== Debug: Significance Status ===\n")
  print(data.frame(
    Variable = plot_data$Variable_clean,
    P_value = plot_data$P_value,
    Significant = plot_data$Significant
  ))
  
  # Create table text for forestplot
  tabletext <- cbind(
    c("Variable", plot_data$Variable_clean),
    c("HR (95% CI)", plot_data$HR_CI),
    c("P-value", plot_data$P_text)
  )
  
  # Prepare mean (HR) and confidence intervals
  mean_values <- c(NA, plot_data$HR)
  lower_values <- c(NA, plot_data$Lower_CI)
  upper_values <- c(NA, plot_data$Upper_CI)
  
  # Create colors and box sizes - DEBUG VERSION
  colors_box <- c(NA, ifelse(plot_data$Significant, "#e74c3c", "#3498db"))
  colors_lines <- c(NA, ifelse(plot_data$Significant, "#e74c3c", "#3498db"))
  box_sizes <- c(NA, ifelse(plot_data$Significant, 0.3, 0.2))
  
  # Debug: Print color assignments
  cat("\n=== Debug: Color Assignments ===\n")
  print(data.frame(
    Row = 1:length(colors_box),
    Variable = c("Header", plot_data$Variable_clean),
    Significant = c(NA, plot_data$Significant),
    Color = colors_box,
    BoxSize = box_sizes
  ))
  
  # Determine x-axis range
  x_min <- max(0.1, min(plot_data$Lower_CI, na.rm = TRUE) * 0.8)
  x_max <- min(10, max(plot_data$Upper_CI, na.rm = TRUE) * 1.2)
  
  # Create forest plot
  forest_plot <- forestplot(
    tabletext,
    mean = mean_values,
    lower = lower_values, 
    upper = upper_values,
    zero = 1,
    graph.pos = 2,
    hrzl_lines = list(
      "1" = gpar(lty = 1, lwd = 2),
      "2" = gpar(lty = 2, lwd = 1)
    ),
    graphwidth = unit(.3, "npc"),
    xlab = "Hazard Ratio (95% CI)",
    xticks = c(0.25, 0.5, 1, 2, 4),
    clip = c(x_min, x_max),
    is.summary = c(TRUE, rep(FALSE, nrow(plot_data))),
    
    txt_gp = fpTxtGp(
      label = gpar(cex = 0.9),
      ticks = gpar(cex = 0.8), 
      xlab = gpar(cex = 1.1, fontface = "bold"),
      title = gpar(cex = 1.3, fontface = "bold")
    ),
    
    lwd.zero = 2,
    lwd.ci = 1.5,
    lwd.xaxis = 1,
    lty.ci = 1,
    ci.vertices = TRUE,
    ci.vertices.height = 0.1,
    
    title = paste0("Forest Plot: ", tissue_name, "\nCox Regression Hazard Ratios"),
    
    # Try explicit color specification
    col = fpColors(
      box = colors_box,
      lines = colors_lines,
      zero = "black"
    ),
    
    boxsize = box_sizes
  )
  
  if (!is.null(save_path)) {
    dev.off()
    cat("Forest plot saved to:", save_path, "\n")
  }
  
  # Print summary
  cat("\n=== Forest Plot Summary ===\n")
  cat("Total variables:", nrow(plot_data), "\n")
  cat("Significant variables (p <", alpha_level, "):", sum(plot_data$Significant), "\n")
  cat("HR range:", round(min(plot_data$HR), 3), "to", round(max(plot_data$HR), 3), "\n")
  
  return (NULL)
}

# ===============================================================================
# Function to compare markers between PCME types for each cell subtype
# ===============================================================================

compare_pcme_markers <- function(cells_with_pcme,expr_data,subtype,available_markers) {
  
  idx <- (cells_with_pcme$sub_celltype == subtype)
  subtype_cells <- cells_with_pcme[idx,]
  
  if(nrow(subtype_cells) < 10) {
    return(NULL)  # Skip subtypes with too few cells
  }
  
  pcme_counts <- table(subtype_cells$PCME_binary)
  if(length(pcme_counts) < 2 || min(pcme_counts) < 5) {
    return(NULL)  # Skip if not enough cells in both PCME types
  }
  
  cat("Processing subtype:", subtype, "- cells:", nrow(subtype_cells), "\n")
  
  # Get expression data for this subtype
  subtype_expr <- expr_data[available_markers, idx]
  
  # Perform wilcoxon test for each marker
  marker_results <- map_dfr(available_markers, function(marker) {
    
    marker_expr <- as.numeric(subtype_expr[marker, ])
    pcme_groups <- subtype_cells$PCME_binary
    
    # Calculate means for each group
    mean_immune <- mean(marker_expr[pcme_groups == "Immune-permissive"], na.rm = TRUE)
    mean_stromal <- mean(marker_expr[pcme_groups == "Stromal-fibrotic"], na.rm = TRUE)
    
    # Wilcoxon test
    if(length(unique(pcme_groups)) == 2) {
      wilcox_result <- wilcox.test(marker_expr ~ pcme_groups)
      
      data.frame(
        cell_subtype = subtype,
        marker = marker,
        mean_immune_permissive = mean_immune,
        mean_stromal_fibrotic = mean_stromal,
        log2_fc = log2((mean_stromal + 1e-6) / (mean_immune + 1e-6)),
        p_value = wilcox_result$p.value,
        n_immune = sum(pcme_groups == "Immune-permissive"),
        n_stromal = sum(pcme_groups == "Stromal-fibrotic")
      )
    } else {
      data.frame(
        cell_subtype = subtype,
        marker = marker,
        mean_immune_permissive = mean_immune,
        mean_stromal_fibrotic = mean_stromal,
        log2_fc = NA,
        p_value = NA,
        n_immune = sum(pcme_groups == "Immune-permissive"),
        n_stromal = sum(pcme_groups == "Stromal-fibrotic")
      )
    }
  })
  
  return(marker_results)
}

# ===============================================================================
# Create PCME Comparison Heatmap
# ===============================================================================
# Function to filter biologically irrelevant combinations
filter_biological_associations <- function(matrix_data, associations) {
  filtered_matrix <- matrix_data
  
  for(cell_type in rownames(matrix_data)) {
    if(cell_type %in% names(associations)) {
      relevant_markers <- associations[[cell_type]]
      # Set irrelevant markers to NA
      irrelevant_markers <- setdiff(colnames(matrix_data), relevant_markers)
      filtered_matrix[cell_type, irrelevant_markers] <- NA
    }
  }
  
  return(filtered_matrix)
}

# Custom function to handle NA values in visualization
create_significance_annotation <- function(mat, size_mat, color_mat) {
  function(j, i, x, y, width, height, fill) {
    if(!is.na(color_mat[i, j]) && 
       !is.na(mat[i, j]) && 
       mat[i, j] == TRUE &&
       !is.na(size_mat[i, j])) {
      
      raw_size <- size_mat[i, j]
      
      # Much smaller points and more selective
      if(raw_size >= 4) {        # p < 0.0001
        point_size <- unit(0.0225, "npc")
        grid.points(x, y, pch = 16, size = point_size, gp = gpar(col = "black"))
      } else if(raw_size >= 3) { # p < 0.001  
        point_size <- unit(0.01875, "npc")
        grid.points(x, y, pch = 16, size = point_size, gp = gpar(col = "black"))
      } else if(raw_size >= 2) { # p < 0.01
        point_size <- unit(0.015, "npc")
        grid.points(x, y, pch = 16, size = point_size, gp = gpar(col = "black"))
      }
      else if(raw_size >= 1) { # p < 0.05
        point_size <- unit(0.0075, "npc")
        grid.points(x, y, pch = 1, size = point_size, gp = gpar(col = "black"))  # Hollow circle
      }    }
  }
}

create_pcme_heatmap <- function(heatmap_data, output_dir) {
  
  cat("Creating biologically filtered PCME heatmap...\n")
  
  # Define biologically meaningful marker-cell type associations
  # Based on canonical expression patterns and functional relevance
  # Irrelevant combinations (e.g., T cell markers on stromal cells) will be set to NA
  cell_marker_associations <- list(
    # Epithelial Cells (Cholangiocytes)
    "EC" = c("EpCAM", "CA_IX", "Ki67", "GLUT1", "FASN", "PRPS1", "HK2", "Vimentin", "CD274", "VEGF"),
    
    # Macrophages (all subtypes)
    "Macro_HLADR" = c("CD68", "CD45", "CD169", "CD80", "CD11b", "CD163", "HLA_DR", "CD14", "CD11c", "Ki67", "VEGF", "CD274", "CLEC9A"),
    "Macro_CD11b" = c("CD68", "CD45", "CD169", "CD80", "CD11b", "CD163", "HLA_DR", "CD14", "CD11c", "Ki67", "VEGF", "CD274"),
    "Macro_CD163" = c("CD68", "CD45", "CD169", "CD80", "CD11b", "CD163", "HLA_DR", "CD14", "CD11c", "Ki67", "VEGF", "CD274"),
    "Macro_Other" = c("CD68", "CD45", "CD169", "CD80", "CD11b", "CD163", "HLA_DR", "CD14", "CD11c", "Ki67", "VEGF", "CD274"),
    
    # Stromal Cells
    "SC_Collagen" = c("Vimentin", "Collagen_I", "Ki67", "VEGF", "GLUT1", "FASN", "PRPS1", "HK2", "CD274", "Alpha_SMA"),
    "SC_Collagen_Vimentin" = c("Vimentin", "Collagen_I", "Ki67", "VEGF", "GLUT1", "FASN", "PRPS1", "HK2", "CD274", "Alpha_SMA"),
    "SC_aSMA_Vimentin" = c("Vimentin", "Alpha_SMA", "Ki67", "VEGF", "GLUT1", "FASN", "PRPS1", "HK2", "CD274", "Collagen_I"),
    "SC_Vimentin" = c("Vimentin", "Ki67", "VEGF", "GLUT1", "FASN", "PRPS1", "HK2", "CD274", "Alpha_SMA", "Collagen_I"),
    "SC_aSMA_Collagen" = c("Alpha_SMA", "Collagen_I", "Vimentin", "Ki67", "VEGF", "GLUT1", "FASN", "PRPS1", "HK2", "CD274"),
    
    # Cancer Associated Fibroblasts
    "CAF" = c("FAP", "Vimentin", "Alpha_SMA", "Ki67", "VEGF", "GLUT1", "FASN", "PRPS1", "HK2", "CD274", "Collagen_I"),
    
    # B Cells
    "B" = c("CD20", "CD45", "Ki67", "CD80", "HLA_DR", "CD27", "CD274"),
    
    # Monocytes
    "Mono_CD14" = c("CD14", "CD45", "CD11b", "CD68", "Ki67", "HLA_DR", "CD80", "CD11c", "CD274"),
    
    # T Cells
    "CD4T" = c("CD3", "CD4", "CD45", "CD127", "TIGIT", "CD366", "CD279", "CD27", "Ki67", "CD80", "HLA_DR"),
    "CD8T" = c("CD3", "CD8a", "CD45", "CD127", "TIGIT", "CD366", "CD279", "CD27", "Ki67", "CD57"),
    "Treg" = c("CD3", "CD4", "CD45", "FoxP3", "TIGIT", "CD366", "CD279", "CD27", "Ki67", "CD127"),
    
    # NK Cells
    "NK" = c("CD16", "CD57", "CD45", "Ki67", "CD366", "TIGIT"),
    
    # Other Immune (broad category)
    "Other_Immune" = c("CD45", "Ki67", "CD3", "CD68", "CD80", "HLA_DR", "CD274", "CD11b", "CD11c", "CLEC9A")
  )
  
  # Prepare matrix for heatmap
  # Use log2_fc for color, neg_log10_p for size
  color_matrix <- heatmap_data %>%
    select(cell_subtype, marker, log2_fc) %>%
    pivot_wider(names_from = marker, values_from = log2_fc, values_fill = 0) %>%
    tibble::column_to_rownames("cell_subtype") %>%
    as.matrix()
  
  size_matrix <- heatmap_data %>%
    select(cell_subtype, marker, neg_log10_p_capped) %>%
    pivot_wider(names_from = marker, values_from = neg_log10_p_capped, values_fill = 0) %>%
    tibble::column_to_rownames("cell_subtype") %>%
    as.matrix()
  
  significance_matrix <- heatmap_data %>%
    select(cell_subtype, marker, significant) %>%
    pivot_wider(names_from = marker, values_from = significant, values_fill = FALSE) %>%
    tibble::column_to_rownames("cell_subtype") %>%
    as.matrix()
  
  # Ensure matrices have same dimensions
  common_rows <- intersect(rownames(color_matrix), rownames(size_matrix))
  common_cols <- intersect(colnames(color_matrix), colnames(size_matrix))
  
  color_matrix <- color_matrix[common_rows, common_cols]
  size_matrix <- size_matrix[common_rows, common_cols]
  significance_matrix <- significance_matrix[common_rows, common_cols]
  
  # Apply biological filtering
  color_matrix <- filter_biological_associations(color_matrix, cell_marker_associations)
  size_matrix <- filter_biological_associations(size_matrix, cell_marker_associations)
  significance_matrix <- filter_biological_associations(significance_matrix, cell_marker_associations)
  
  # Color function for log2 fold change (handle NA values)
  col_fun <- colorRamp2(c(-2, 0, 2), c("#3498db", "white", "#e74c3c"))
  
  
  # Group markers by functional categories for better visualization
  marker_groups <- list(
    "Pan_Immune" = c("CD45", "Ki67"),
    "T_Cell" = c("CD3", "CD4", "CD8a", "CD127", "CD27", "FoxP3"),
    "T_Exhaustion" = c("TIGIT", "CD366", "CD279"),
    "Myeloid" = c("CD68", "CD14", "CD11b", "CD11c", "CD163", "CD169", "HLA_DR"),
    "APC_Function" = c("CD80", "CLEC9A"),
    "B_NK" = c("CD20", "CD16", "CD57"),
    "Epithelial" = c("EpCAM", "CA_IX"),
    "Stromal" = c("Vimentin", "Alpha_SMA", "Collagen_I", "FAP"),
    "Metabolism" = c("GLUT1", "FASN", "PRPS1", "HK2"),
    "Signaling" = c("VEGF", "CD274"),
    "Other" = c("DNA1", "DNA2")
  )
  
  # Reorder markers by functional groups
  ordered_markers <- unlist(marker_groups)
  ordered_markers <- intersect(ordered_markers, colnames(color_matrix))
  remaining_markers <- setdiff(colnames(color_matrix), ordered_markers)
  final_marker_order <- c(ordered_markers, remaining_markers)
  
  # Reorder matrices
  color_matrix <- color_matrix[, final_marker_order]
  size_matrix <- size_matrix[, final_marker_order]
  significance_matrix <- significance_matrix[, final_marker_order]
  
  # Group cell types by function
  cell_groups <- list(
    "Epithelial" = c("EC"),
    "T_Cells" = c("CD4T", "CD8T", "Treg"),
    "Myeloid" = c("Macro_HLADR", "Macro_CD11b", "Macro_CD163", "Macro_Other", "Mono_CD14"),
    "Other_Immune" = c("B", "NK", "Other_Immune"),
    "Stromal" = c("SC_Collagen", "SC_Collagen_Vimentin", "SC_aSMA_Vimentin", 
                  "SC_Vimentin", "SC_aSMA_Collagen", "CAF")
  )
  
  # Reorder rows by cell type groups
  ordered_cells <- unlist(cell_groups)
  ordered_cells <- intersect(ordered_cells, rownames(color_matrix))
  remaining_cells <- setdiff(rownames(color_matrix), ordered_cells)
  final_cell_order <- c(ordered_cells, remaining_cells)
  
  # Reorder matrices
  color_matrix <- color_matrix[final_cell_order, ]
  size_matrix <- size_matrix[final_cell_order, ]
  significance_matrix <- significance_matrix[final_cell_order, ]
  
  # Create column annotations for marker groups
  marker_annotation <- rep(names(marker_groups), lengths(marker_groups))
  names(marker_annotation) <- unlist(marker_groups)
  marker_annotation <- marker_annotation[colnames(color_matrix)]
  marker_annotation[is.na(marker_annotation)] <- "Other"
  
  col_ha <- HeatmapAnnotation(
    Marker_Type = marker_annotation,
    col = list(Marker_Type = c("Pan_Immune" = "#FF6B6B", "T_Cell" = "#4ECDC4", 
                               "T_Exhaustion" = "#45B7D1", "Myeloid" = "#96CEB4",
                               "APC_Function" = "#FFEAA7", "B_NK" = "#DDA0DD",
                               "Epithelial" = "#FFA07A", "Stromal" = "#98D8C8",
                               "Metabolism" = "#F7DC6F", "Signaling" = "#BB8FCE",
                               "Other" = "lightgray")),
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  # Create row annotations for cell groups
  cell_annotation <- rep(names(cell_groups), lengths(cell_groups))
  names(cell_annotation) <- unlist(cell_groups)
  cell_annotation <- cell_annotation[rownames(color_matrix)]
  cell_annotation[is.na(cell_annotation)] <- "Other"
  
  row_ha <- rowAnnotation(
    Cell_Type = cell_annotation,
    col = list(Cell_Type = c("Epithelial" = "#FFB6C1", "T_Cells" = "#87CEEB", 
                             "Myeloid" = "#98FB98", "Other_Immune" = "#DDA0DD",
                             "Stromal" = "#F0E68C", "Other" = "lightgray")),
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  # Create heatmap with biological filtering and functional grouping
  pdf(file.path(output_dir, "pcme_subtype_marker_heatmap_filtered.pdf"), width = 20, height = 12)
  
  ht <- Heatmap(
    color_matrix,
    name = "log2(FC)",
    col = col_fun,
    na_col = "lightgray",
    
    cell_fun = create_significance_annotation(significance_matrix, size_matrix, color_matrix),
    
    top_annotation = col_ha,
    right_annotation = row_ha,
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    row_names_side = "left",
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 9),
    
    heatmap_legend_param = list(
      title = "log2(PCME-S/PCME-I)",
      title_position = "topcenter",
      at = c(-2, -1, 0, 1, 2),
      labels = c("-2", "-1", "0", "1", "2")
    ),
    
    width = unit(14, "cm"),
    height = unit(10, "cm"),
    border = TRUE
  )
  
  # Add custom legend for significance dots
  lgd_sig <- Legend(
    at = c("p < 0.05", "p < 0.01", "p < 0.001", "p < 0.0001"), 
    title = "Significance",
    legend_gp = gpar(fill = "black"),
    type = "points",
    pch = 16,
    size = unit(c(0.25, 0.4, 0.5, 0.7), "cm"),  # Match the scaling above
    grid_width = unit(1, "cm"),
    grid_height = unit(1, "cm")
  )
  
  # Add legend for NA values
  lgd_na <- Legend(
    at = "Not expressed",
    title = "Biological\nRelevance",
    legend_gp = gpar(fill = "lightgray"),
    type = "grid",
    grid_width = unit(0.8, "cm"),
    grid_height = unit(0.4, "cm")
  )
  
  draw(ht, annotation_legend_list = list(lgd_sig, lgd_na))
  
  dev.off()
  
  # Create a summary table of biologically relevant significant changes
  relevant_changes <- heatmap_data %>%
    filter(significant) %>%
    mutate(
      cell_type_group = case_when(
        cell_subtype == "EC" ~ "Epithelial",
        cell_subtype %in% c("CD4T", "CD8T", "Treg") ~ "T_Cells",
        cell_subtype %in% c("Macro_HLADR", "Macro_CD11b", "Macro_CD163", "Macro_Other", "Mono_CD14") ~ "Myeloid",
        cell_subtype %in% c("B", "NK", "Other_Immune") ~ "Other_Immune",
        TRUE ~ "Stromal"
      )
    )
  
  # Filter for biologically relevant combinations
  relevant_significant <- relevant_changes %>%
    filter(
      # Check if this combination should be included based on our associations
      pmap_lgl(list(cell_subtype, marker), function(ct, mk) {
        if(ct %in% names(cell_marker_associations)) {
          mk %in% cell_marker_associations[[ct]]
        } else {
          TRUE  # Include if cell type not in our predefined list
        }
      })
    ) %>%
    arrange(cell_type_group, desc(abs(log2_fc)))
  
  # Save filtered significant results
  write.csv(relevant_significant, 
            file.path(output_dir, "pcme_biologically_relevant_significant_changes.csv"), 
            row.names = FALSE)
  
  cat("Biologically relevant significant changes:", nrow(relevant_significant), 
      "out of", nrow(heatmap_data %>% filter(significant)), "total significant\n")
  
  # Create simplified dot plot version (biologically relevant only)
  dot_plot_data <- relevant_significant %>%
    mutate(
      marker = factor(marker, levels = final_marker_order),
      cell_subtype = factor(cell_subtype, levels = final_cell_order)
    )
  
  if(nrow(dot_plot_data) > 0) {
    # Create main dot plot
    dot_plot <- ggplot(dot_plot_data, aes(x = marker, y = cell_subtype)) +
      geom_point(aes(color = log2_fc, size = neg_log10_p_capped), alpha = 0.8) +
      scale_color_gradient2(low = "#3498db", mid = "white", high = "#e74c3c", 
                            midpoint = 0, name = "log2(FC)") +
      scale_size_continuous(range = c(2, 10), name = "-log10(p_adj)") +
      labs(
        title = "Biologically Relevant PCME Marker Differences",
        subtitle = "PCME-S vs PCME-I (significant & biologically relevant only)",
        x = "Markers",
        y = "Cell Subtypes"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right"
      ) +
      # Add functional group separators
      geom_vline(xintercept = cumsum(lengths(marker_groups))[1:length(marker_groups)-1] + 0.5, 
                 linetype = "dashed", alpha = 0.3) +
      geom_hline(yintercept = cumsum(lengths(cell_groups))[1:length(cell_groups)-1] + 0.5, 
                 linetype = "dashed", alpha = 0.3)
    
    ggsave(file.path(output_dir, "pcme_subtype_marker_dotplot_filtered.pdf"), 
           dot_plot, width = 16, height = 10, dpi = 300)
    
    # Create summary plot by cell type groups
    summary_plot_data <- relevant_significant %>%
      group_by(cell_type_group, direction) %>%
      summarise(
        count = n(),
        mean_abs_fc = mean(abs(log2_fc)),
        .groups = "drop"
      )
    
    summary_plot <- ggplot(summary_plot_data, aes(x = cell_type_group, y = count, fill = direction)) +
      geom_col(position = "dodge", alpha = 0.8) +
      geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.5) +
      scale_fill_manual(values = c("Higher in PCME-S" = "#e74c3c", "Higher in PCME-I" = "#3498db")) +
      labs(
        title = "Distribution of Significant PCME Differences by Cell Type",
        x = "Cell Type Groups",
        y = "Number of Significant Markers",
        fill = "Direction"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    ggsave(file.path(output_dir, "pcme_summary_by_celltype.pdf"), 
           summary_plot, width = 10, height = 6, dpi = 300)
  }
  
  # Print summary statistics
  cat("\nSummary of biologically relevant findings:\n")
  if(nrow(relevant_significant) > 0) {
    cat("Cell type groups with most differences:\n")
    print(relevant_significant %>% 
            count(cell_type_group, sort = TRUE))
    
    cat("\nMarkers with most differences across cell types:\n")
    print(relevant_significant %>% 
            count(marker, sort = TRUE) %>% 
            head(10))
    
    cat("\nDirection of changes:\n")
    print(table(relevant_significant$direction))
  }
}
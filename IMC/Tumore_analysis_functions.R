library(SpatialExperiment)
library(imcRtools)
library(igraph)

library(ggplot2)
library(ggridges)
library(ggalluvial)
library(RColorBrewer)
library(ggpubr)
library(ggraph)
library(viridis)
library(patchwork)
library(gridExtra)

library(dplyr)
library(tidyr)

library(survival)
library(survminer)
library(maxstat)
library(forestplot)

# Function to calculate distances between clusters
calculate_cluster_distances <- function(cluster1_coords, cluster2_coords) {
  # Calculate all pairwise distances between cells in two clusters
  dist_matrix <- as.matrix(dist(rbind(cluster1_coords, cluster2_coords)))
  
  # Extract the submatrix for cross-cluster distances
  n1 <- nrow(cluster1_coords)
  n2 <- nrow(cluster2_coords)
  cross_distances <- dist_matrix[1:n1, (n1+1):(n1+n2)]
  
  # Return different distance metrics
  list(
    min_distance = min(cross_distances),
    mean_distance = mean(cross_distances),
    median_distance = median(cross_distances),
    centroid_distance = sqrt(sum((colMeans(cluster1_coords) - colMeans(cluster2_coords))^2))
  )
}

# Function to perform Wilcoxon tests with BH adjustment
perform_wilcox_tests <- function(data, group_name) {
  # Perform tests for each tissue-subtype combination
  test_results <- data %>%
    group_by(Tissue, sub_celltype) %>%
    summarise(
      n0 = sum(RFS_status == 0),
      n1 = sum(RFS_status == 1),
      .groups = "drop"
    ) %>%
    filter(n0 >= 3 & n1 >= 3) %>%  # Only test if sufficient samples
    rowwise() %>%
    mutate(
      p_value = {
        current_tissue <- Tissue
        current_subtype <- sub_celltype
        tissue_data <- data %>% 
          filter(Tissue == current_tissue, sub_celltype == current_subtype)
        if(nrow(tissue_data) > 5) {
          test_result <- wilcox.test(fraction ~ RFS_status, data = tissue_data)
          test_result$p.value
        } else {
          NA_real_
        }
      }
    ) %>%
    ungroup() %>%
    filter(!is.na(p_value))
  
  # Apply BH adjustment
  if(nrow(test_results) > 0) {
    test_results$p_adjusted <- p.adjust(test_results$p_value, method = "BH")
    test_results$significance <- case_when(
      test_results$p_adjusted < 0.001 ~ "***",
      test_results$p_adjusted < 0.01 ~ "**", 
      test_results$p_adjusted < 0.05 ~ "*",
      TRUE ~ "ns"
    )
    test_results$group <- group_name
  }
  
  return(test_results)
}

# Function to calculate budding metrics for a single ROI
calculate_roi_budding_metrics <- function(roi_components) {
  # Check if roi_components is empty or NULL
  if (is.null(roi_components) || length(roi_components) == 0) {
    return(data.frame(
      budding_number = 0,
      total_size = 0,
      mean_size = 0
    ))
  }
  
  # Calculate budding number (number of components)
  budding_number <- length(roi_components)
  
  # Calculate size of each component
  component_sizes <- sapply(roi_components, length)
  
  # Calculate total size (sum of all elements across all components)
  total_size <- sum(component_sizes)
  
  # Calculate mean size (average elements per component)
  mean_size <- mean(component_sizes)
  
  return(data.frame(
    budding_number = budding_number,
    total_size = total_size,
    mean_size = round(mean_size, 2)
  ))
}

# Function to create boxplot with annotations
create_comparison_plot <- function(data, test_data, title_suffix, group_filter) {
  
  # Filter test results for this group
  filtered_tests <- test_data %>% filter(group == group_filter)
  
  # Create significance summary for subtitle
  sig_summary <- filtered_tests %>%
    summarise(
      total_tests = n(),
      significant = sum(p_adjusted < 0.05),
      .groups = "drop"
    )
  
  subtitle_text <- paste0("Wilcoxon tests with BH adjustment: ", 
                          sig_summary$significant, "/", sig_summary$total_tests, 
                          " comparisons significant (p < 0.05)")
  
  p <- ggplot(data, aes(x = factor(RFS_status), y = fraction, fill = factor(RFS_status))) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
    geom_point(position = position_jitter(width = 0.2, seed = 123), 
               alpha = 0.6, size = 1) +
    facet_grid(sub_celltype ~ Tissue, scales = "free_y") +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.format",
                       size = 3,
                       vjust = 0.5) +
    scale_fill_manual(values = c("0" = "#0073c2b2", "1" = "#efc000b2"),
                      labels = c("0" = "No Early Relapse", "1" = "Early Relapse")) +
    scale_x_discrete(labels = c("0" = "No\nRelapse", "1" = "Early\nRelapse")) +
    labs(
      title = paste0("Malignant Subtype Fractions by RFS Status - ", title_suffix),
      subtitle = subtitle_text,
      x = "RFS Status",
      y = "Fraction",
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
  
  return(p)
}

# Function to find optimal cutoff using maxstat
find_optimal_cutoff <- function(data, gene_name, time_var = "RFS_time", event_var = "RFS_status") {
  clean_data <- data[!is.na(data[[gene_name]]) & !is.na(data[[time_var]]) & !is.na(data[[event_var]]), ]
  
  if(nrow(clean_data) < 10) return(median(data[[gene_name]], na.rm = TRUE))
  
  maxstat_formula <- as.formula(paste("Surv(", time_var, ",", event_var, ") ~", gene_name))
  
  tryCatch({
    maxstat_result <- maxstat.test(maxstat_formula, data = clean_data, smethod = "LogRank")
    return(maxstat_result$estimate)
  }, error = function(e) {
    return(median(data[[gene_name]], na.rm = TRUE))
  })
}

# Function to create forest plot from Cox regression results
createCoxForestPlot <- function(cox_model, 
                                var_display_names = NULL, 
                                var_groups = NULL,
                                plot_title = "Cox Regression Forest Plot", 
                                savePath = NULL, 
                                filename = "Cox_forest_plot.pdf") {
  
  # Extract coefficients and statistics from Cox model
  coef_summary <- summary(cox_model)
  coef_table <- coef_summary$coefficients
  conf_int <- coef_summary$conf.int
  
  # Get variable names from the model
  var_names <- rownames(coef_table)
  n_vars <- length(var_names)
  
  # Create display names if not provided
  if(is.null(var_display_names)){
    var_display_names <- var_names
    # Try to clean up variable names for better display
    var_display_names <- gsub("_", " ", var_display_names)
    var_display_names <- tools::toTitleCase(var_display_names)
  }
  
  # Check if display names match number of variables
  if(length(var_display_names) != n_vars) {
    warning("Length of var_display_names doesn't match number of variables. Using original variable names.")
    var_display_names <- var_names
  }
  
  # Extract key statistics
  HR <- coef_table[, "exp(coef)"]
  lower_CI <- conf_int[, "lower .95"]
  upper_CI <- conf_int[, "upper .95"]
  p_values <- coef_table[, "Pr(>|z|)"]
  
  # Format HR(95%CI) for display
  HR_CI_text <- paste0(sprintf("%.2f", HR), " (", 
                       sprintf("%.2f", lower_CI), "-", 
                       sprintf("%.2f", upper_CI), ")")
  
  # Format p-values for display
  p_display <- ifelse(p_values < 0.001, "<0.001", 
                      ifelse(p_values < 0.01, sprintf("%.3f", p_values),
                             sprintf("%.3f", p_values)))
  
  # Helper function to insert section headers
  ins <- function(x) {
    c(as.character(x), NA, NA)
  }
  
  # Build result matrix
  if(is.null(var_groups)) {
    # No grouping - show all variables in one section
    result_df <- rbind(
      c("Variables", "HR(95%CI)", "p-value"),
      cbind(var_display_names, HR_CI_text, p_display),
      c(NA, NA, NA)
    )
    
    # Create vectors for forestplot function
    mean_values <- c(NA, HR, NA)
    lower_values <- c(NA, lower_CI, NA)
    upper_values <- c(NA, upper_CI, NA)
    
    # Create is_summary vector
    is_summary_vector <- c(TRUE, rep(FALSE, n_vars), TRUE)
    
  } else {
    # With grouping
    result_df <- c("Variables", "HR(95%CI)", "p-value")
    mean_values <- NA
    lower_values <- NA
    upper_values <- NA
    is_summary_vector <- TRUE
    
    # Process each group
    for(group_name in names(var_groups)) {
      group_indices <- var_groups[[group_name]]
      
      # Validate indices
      group_indices <- group_indices[group_indices <= n_vars & group_indices > 0]
      
      if(length(group_indices) > 0) {
        # Add group header
        result_df <- rbind(result_df, ins(group_name))
        mean_values <- c(mean_values, NA)
        lower_values <- c(lower_values, NA)
        upper_values <- c(upper_values, NA)
        is_summary_vector <- c(is_summary_vector, TRUE)
        
        # Add group variables
        result_df <- rbind(result_df, 
                           cbind(var_display_names[group_indices], 
                                 HR_CI_text[group_indices], 
                                 p_display[group_indices]))
        mean_values <- c(mean_values, HR[group_indices])
        lower_values <- c(lower_values, lower_CI[group_indices])
        upper_values <- c(upper_values, upper_CI[group_indices])
        is_summary_vector <- c(is_summary_vector, rep(FALSE, length(group_indices)))
      }
    }
    
    # Add final empty row
    result_df <- rbind(result_df, c(NA, NA, NA))
    mean_values <- c(mean_values, NA)
    lower_values <- c(lower_values, NA)
    upper_values <- c(upper_values, NA)
    is_summary_vector <- c(is_summary_vector, TRUE)
  }
  
  # Ensure all entries are characters
  result_df <- apply(result_df, 2, as.character)
  
  # Create the forest plot
  p <- forestplot(result_df,
                  mean = mean_values,
                  lower = lower_values,
                  upper = upper_values,
                  zero = 1,
                  boxsize = 0.4,
                  graph.pos = "right",
                  hrzl_lines = list(
                    "1" = gpar(lty = 1, lwd = 2),
                    "2" = gpar(lty = 2)
                  ),
                  graphwidth = unit(.3, "npc"),
                  xlab = "Hazard Ratio",
                  xticks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
                  is.summary = is_summary_vector,
                  txt_gp = fpTxtGp(
                    label = gpar(cex = 0.9),
                    ticks = gpar(cex = 0.8),
                    xlab = gpar(cex = 1.2),
                    title = gpar(cex = 1.4)
                  ),
                  lwd.zero = 2,
                  lwd.ci = 1.5,
                  lwd.xaxis = 2,
                  lty.ci = 1,
                  ci.vertices = TRUE,
                  ci.vertices.height = 0.15,
                  clip = c(0.05, 15),
                  lineheight = unit(8, "mm"),
                  line.margin = unit(5, "mm"),
                  colgap = unit(4, "mm"),
                  fn.ci_norm = "fpDrawNormalCI",
                  title = plot_title,
                  col = fpColors(
                    box = "#D32F2F",
                    lines = "#D32F2F",
                    zero = "black",
                    text = "black",
                    summary = "black"
                  )
  )
  
  # Save plot if path provided
  if (!is.null(savePath)) {
    pdf(file.path(savePath, filename), width = 10, height = 7.5)
    print(p)
    dev.off()
    cat("Forest plot saved to:", file.path(savePath, filename), "\n")
  }
  
  # Display plot
  print(p)
  
  # Return summary data
  summary_data <- data.frame(
    Variable = var_display_names,
    HR = HR,
    Lower_CI = lower_CI,
    Upper_CI = upper_CI,
    P_value = p_values,
    Significance = ifelse(p_values < 0.001, "***",
                          ifelse(p_values < 0.01, "**",
                                 ifelse(p_values < 0.05, "*", "")))
  )
  
  return(summary_data)
}

# Function to calculate correlation matrix with p-values
calculate_correlation_with_pval <- function(data, method = "spearman") {
  # Extract only marker columns
  marker_data <- data[, marker_columns]
  
  # Remove rows with any NA values
  marker_data <- marker_data[complete.cases(marker_data), ]
  
  if (nrow(marker_data) < 10) {
    return(NULL)  # Not enough data for reliable correlation
  }
  
  # Calculate correlation matrix and p-values using Hmisc
  cor_result <- rcorr(as.matrix(marker_data), type = method)
  
  # Return both correlation matrix and p-values
  return(list(
    correlations = cor_result$r,
    p_values = cor_result$P,
    n_samples = nrow(marker_data)
  ))
}

# Function to create significance stars
get_significance_stars <- function(p_values) {
  stars <- matrix("", nrow = nrow(p_values), ncol = ncol(p_values))
  stars[p_values <= 0.001] <- "***"
  stars[p_values > 0.001 & p_values <= 0.01] <- "**"
  stars[p_values > 0.01 & p_values <= 0.05] <- "*"
  stars[p_values > 0.05] <- ""
  
  # Set diagonal to empty (self-correlation)
  diag(stars) <- ""
  
  return(stars)
}

# Function to extract GLUT1 correlations for chord plot
extract_glut1_correlations <- function(data, tissue_name) {
  # Extract only marker columns
  marker_data <- data[, marker_columns]
  
  # Remove rows with any NA values
  marker_data <- marker_data[complete.cases(marker_data), ]
  
  if (nrow(marker_data) < 10) {
    return(NULL)  # Not enough data
  }
  
  # Calculate correlations with GLUT1
  glut1_cors <- cor(marker_data$GLUT1, marker_data, method = "spearman")
  
  # Calculate p-values for GLUT1 correlations
  glut1_pvals <- sapply(marker_columns, function(marker) {
    if (marker == "GLUT1") return(NA)  # Self-correlation
    cor.test(marker_data$GLUT1, marker_data[[marker]], method = "spearman")$p.value
  })
  
  # Create data frame for chord plot
  chord_data <- data.frame(
    from = "GLUT1",
    to = marker_columns[marker_columns != "GLUT1"],
    correlation = as.numeric(glut1_cors[marker_columns != "GLUT1"]),
    p_value = glut1_pvals[marker_columns != "GLUT1"],
    tissue = tissue_name,
    stringsAsFactors = FALSE
  )
  
  # Add transparency based on p-value (more significant = more opaque)
  chord_data$transparency <- pmax(0.2, 1 - chord_data$p_value)  # Min transparency of 0.2
  
  # Add chord width based on absolute correlation
  chord_data$width <- abs(chord_data$correlation) * 10  # Scale for visibility
  
  # Add direction for color coding
  chord_data$direction <- ifelse(chord_data$correlation > 0, "positive", "negative")
  
  return(list(
    data = chord_data,
    n_samples = nrow(marker_data)
  ))
}

# Function to create chord plot for GLUT1 correlations
create_glut1_chord_plot <- function(chord_result, tissue_name) {
  
  chord_data <- chord_result$data
  n_samples <- chord_result$n_samples
  
  # Prepare data matrix for circlize
  # Create adjacency matrix
  markers <- unique(c(chord_data$from, chord_data$to))
  n_markers <- length(markers)
  
  # Create matrix
  mat <- matrix(0, nrow = n_markers, ncol = n_markers)
  rownames(mat) <- markers
  colnames(mat) <- markers
  
  # Fill matrix with correlation values (use absolute values for chord thickness)
  for (i in 1:nrow(chord_data)) {
    from_idx <- which(markers == chord_data$from[i])
    to_idx <- which(markers == chord_data$to[i])
    mat[from_idx, to_idx] <- abs(chord_data$correlation[i])
  }
  
  # Create color matrix for chords based on correlation direction
  col_mat <- matrix(NA, nrow = n_markers, ncol = n_markers)
  rownames(col_mat) <- markers
  colnames(col_mat) <- markers
  
  # Fill color matrix
  for (i in 1:nrow(chord_data)) {
    from_idx <- which(markers == chord_data$from[i])
    to_idx <- which(markers == chord_data$to[i])
    
    # Color based on correlation direction
    if (chord_data$direction[i] == "positive") {
      col_mat[from_idx, to_idx] <- "#cd534c"  # Red for positive
    } else {
      col_mat[from_idx, to_idx] <- "#7ca6dc"  # Blue for negative
    }
  }
  
  # Set up colors for markers (sector colors)
  marker_colors <- c(
    "GLUT1" = "#ff6b35",      # Orange for GLUT1 (central)
    "EpCAM" = "#4a90e2",      # Blue for epithelial
    "Vimentin" = "#7b68ee",   # Purple for mesenchymal
    "Ki67" = "#32cd32",       # Green for proliferation
    "VEGF" = "#ff69b4",       # Pink for angiogenesis
    "CA_IX" = "#8b4513",      # Brown for hypoxia
    "HK2" = "#ffa500",        # Orange for metabolism
    "FASN" = "#ff8c00",       # Dark orange for metabolism
    "CD279" = "#dc143c"       # Red for immune
  )
  
  # Start plotting
  # Clear any existing plots
  circos.clear()
  
  # Set up circos parameters
  circos.par(start.degree = 90, gap.degree = 4, 
             track.margin = c(-0.1, 0.1), cell.padding = c(0.02, 0, 0.02, 0))
  
  # Create chord diagram with simplified coloring
  chordDiagram(mat, 
               grid.col = marker_colors[markers],
               col = col_mat,
               transparency = 0.3,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
  
  # Add marker labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                cex = 0.8, font = 2)
  }, bg.border = NA)
  
  # Add title
  title(main = paste("GLUT1 Correlations -", tissue_name), 
        sub = paste("n =", n_samples, "samples"), 
        cex.main = 1.2, cex.sub = 0.9)
  
  # Add legend for correlation direction
  legend("bottomleft", 
         legend = c("Positive correlation", "Negative correlation"),
         fill = c("#cd534c", "#7ca6dc"),
         cex = 0.7, bty = "n")
  
  # Add legend for chord thickness
  legend("bottomright",
         legend = c("Chord thickness = |correlation|"),
         bty = "n", cex = 0.7)
  
  return(chord_data)
}

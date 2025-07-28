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
createCoxForestPlot <- function(cox_model, var_display_names = NULL, plot_title = "Multi-variables Cox Regression Forest Plot", 
                                savePath = NULL, filename = "Cox_forest_plot.pdf") {
  
  # Extract coefficients and statistics from Cox model
  coef_summary <- summary(cox_model)
  coef_table <- coef_summary$coefficients
  conf_int <- coef_summary$conf.int
  
  # Get variable names (remove reference levels if factor)
  var_names <- rownames(coef_table)
  
  # Create formatted variable names for display
  if(is.null(var_display_names)){
    var_display_names <- c(
      "Treatment (Combo vs Chemo)",
      "FASN", "GLUT1", "EpCAM", 
      "Age", "Fong Score", "TBS", "CEA (log)", "CA199 (log)", "CRLM Size", "KRAS Mutation"
    )
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
  
  # Create result matrix similar to original function
  # Separate biomarkers from clinical variables
  biomarker_indices <- c(2, 3, 4)  # FASN, GLUT1, EpCAM
  clinical_indices <- c(1, 5:11)   # Treatment, Age, fong_score, etc.
  
  # Helper function to insert section headers
  ins <- function(x) {
    c(as.character(x), NA, NA)
  }
  
  # Build result matrix with sections
  result_df <- rbind(
    c("Features", "HR(95%CI)", "p-value"),
    ins("Biomarkers"),
    cbind(var_display_names[biomarker_indices], 
          HR_CI_text[biomarker_indices], 
          p_display[biomarker_indices]),
    ins("Clinical Variables"),
    cbind(var_display_names[clinical_indices], 
          HR_CI_text[clinical_indices], 
          p_display[clinical_indices]),
    c(NA, NA, NA)
  )
  
  # Ensure all entries are characters
  result_df <- apply(result_df, 2, as.character)
  
  # Create vectors for forestplot function
  mean_values <- c(NA, NA, HR[biomarker_indices], NA, HR[clinical_indices], NA)
  lower_values <- c(NA, NA, lower_CI[biomarker_indices], NA, lower_CI[clinical_indices], NA)
  upper_values <- c(NA, NA, upper_CI[biomarker_indices], NA, upper_CI[clinical_indices], NA)
  
  # Create is_summary vector (TRUE for section headers, FALSE for data rows)
  is_summary_vector <- c(TRUE, TRUE, rep(FALSE, length(biomarker_indices)), 
                         TRUE, rep(FALSE, length(clinical_indices)), TRUE)
  
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
                    box = "#D32F2F",        # Red for significant
                    lines = "#D32F2F",      # Red lines
                    zero = "black",         # Black reference line
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
  
  # Return summary data for reference
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

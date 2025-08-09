# Peritumor Analysis for CRLM Recurrence

library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(ConsensusClusterPlus)
library(pheatmap)
library(tidyr)

# Only pack highly reusable functions
get_patient_summary <- function(cell_data, group_vars, value_var) {
  cell_data %>%
    group_by(patient_id, !!!syms(group_vars)) %>%
    summarise(
      mean_value = mean(!!sym(value_var), na.rm = TRUE),
      .groups = "drop"
    )
}

# Prepare data for plots with proper direction colors and significance labels
prep_plot_data <- function(data, label_col) {
  data %>%
    mutate(
      direction = ifelse(p_value > 0.05, "n.s.", ifelse(log2_fc > 0, "R_higher", "NR_higher")),
      alpha_val = pmax(0.3, pmin(1, -log10(p_value) / 3)),  # transparency based on p-value
      size_val = pmax(1, pmin(6, -log10(p_value))),         # size based on p-value
      label = ifelse(p_value < 0.05 & abs(log2_fc) > 0.26, !!sym(label_col), "")
    )
}


# Function for Ripley's K analysis
# Helper function to select cells based on type pattern
get_cell_coords <- function(coords, cell_types, celltype_pattern) {
  if (celltype_pattern == "EC") {
    selected <- startsWith(cell_types, "EC")
  } else {
    selected <- cell_types == celltype_pattern
  }
  return(coords[selected, , drop = FALSE])
}

# Helper function to calculate study area
calculate_study_area <- function(coords) {
  x_range <- diff(range(coords[,1]))
  y_range <- diff(range(coords[,2]))
  return(x_range * y_range)
}

# Helper function to calculate distances between two sets of points
calculate_distance_matrix <- function(focal_cells, target_cells) {
  n_focal <- nrow(focal_cells)
  n_target <- nrow(target_cells)
  distances <- matrix(0, nrow = n_focal, ncol = n_target)
  
  for (i in 1:n_focal) {
    dx <- target_cells[,1] - focal_cells[i,1]
    dy <- target_cells[,2] - focal_cells[i,2]
    distances[i,] <- sqrt(dx^2 + dy^2)
  }
  
  return(distances)
}


# Helper function to select cells based on type pattern
get_cell_coords <- function(coords, cell_types, celltype_pattern) {
  if (celltype_pattern == "EC") {
    selected <- startsWith(cell_types, "EC")
  } else {
    selected <- cell_types == celltype_pattern
  }
  return(coords[selected, , drop = FALSE])
}

# Helper function to calculate study area
calculate_study_area <- function(coords) {
  x_range <- diff(range(coords[,1]))
  y_range <- diff(range(coords[,2]))
  return(x_range * y_range)
}

# Helper function to calculate distances between two sets of points
calculate_distance_matrix <- function(focal_cells, target_cells) {
  n_focal <- nrow(focal_cells)
  n_target <- nrow(target_cells)
  distances <- matrix(0, nrow = n_focal, ncol = n_target)
  
  for (i in 1:n_focal) {
    dx <- target_cells[,1] - focal_cells[i,1]
    dy <- target_cells[,2] - focal_cells[i,2]
    distances[i,] <- sqrt(dx^2 + dy^2)
  }
  
  return(distances)
}

# Main function for Ripley's K analysis
perform_ripleys_k <- function(coords, sample_meta, focal_celltype, target_celltype, max_dist = 200) {
  
  # Extract cell types
  cell_types <- sample_meta$sub_celltype
  
  # Get focal and target cell coordinates
  focal_cells <- get_cell_coords(coords, cell_types, focal_celltype)
  target_cells <- get_cell_coords(coords, cell_types, target_celltype)
  
  # Check minimum cell count requirement
  if (nrow(focal_cells) < 5 || nrow(target_cells) < 5) {
    return(NULL)
  }
  
  # Calculate distance matrix once
  distance_matrix <- calculate_distance_matrix(focal_cells, target_cells)
  
  # Set up distance sequence and results
  r_seq <- seq(10, max_dist, by = 10)
  k_values <- numeric(length(r_seq))
  
  # Calculate study area and target cell density
  area <- calculate_study_area(coords)
  lambda <- nrow(target_cells) / area
  n_focal <- nrow(focal_cells)
  
  # Calculate K-values for each distance
  for (i in seq_along(r_seq)) {
    r <- r_seq[i]
    
    # Count pairs within distance r
    within_distance <- distance_matrix <= r
    total_count <- sum(within_distance)
    
    # Ripley's K calculation
    k_values[i] <- total_count / (n_focal * lambda)
  }
  
  # L-function transformation
  l_values <- sqrt(k_values / pi) - r_seq
  
  return(data.frame(
    distance = r_seq,
    K = k_values,
    L = l_values
  ))
}

# Function to add significance annotations
add_significance <- function(data, y_position, p.adjust = TRUE) {
  # Perform Wilcoxon test
  test_result <- wilcox.test(log2_L ~ rfs_status, data = data)
  p_val <- test_result$p.value
  
  if(p.adjust){
    p.adjust <- p.adjust(p.adjust,"BH")
  }
  
  # Determine significance level
  if(p_val < 0.001) {
    sig_text <- "***"
  } else if(p_val < 0.01) {
    sig_text <- "**"
  } else if(p_val < 0.05) {
    sig_text <- "*"
  } else {
    sig_text <- "ns"
  }
  
  return(data.frame(
    x = 1.5,  # Center between the two groups
    y = y_position,
    label = paste0("p = ", format(p_val, digits = 3), " (", sig_text, ")")
  ))
}



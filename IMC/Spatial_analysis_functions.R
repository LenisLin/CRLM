# For spatial analysis

library(SpatialExperiment)
library(imcRtools)
library(igraph)

library(ggplot2)
library(ggridges)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(ggraph)
library(viridis)
library(patchwork)

library(dplyr)
library(tidyr)

library(BiocParallel)
library(lisaClust)
library(scales)

# =============================================================================
## CN analysis
# =============================================================================

# ============================================================================
# SIGNIFICANCE TESTING FOR ENRICHMENT
# ============================================================================
perform_test_for_CN_enrichment <- function(contingency_table, method = "fisher"){
  
  print("Performing significance testing for CN enrichment with method...")
  
  # Get dimensions
  cn_names <- rownames(contingency_table)
  celltype_names <- colnames(contingency_table)
  n_cn <- length(cn_names)
  n_celltype <- length(celltype_names)
  
  # Initialize matrices for results
  pvalue_matrix <- matrix(NA, nrow = n_cn, ncol = n_celltype)
  enrichment_matrix <- matrix(NA, nrow = n_cn, ncol = n_celltype)
  expected_matrix <- matrix(NA, nrow = n_cn, ncol = n_celltype)
  fold_change_matrix <- matrix(NA, nrow = n_cn, ncol = n_celltype)
  
  # Set row and column names
  matrices <- list(pvalue_matrix, enrichment_matrix, expected_matrix, fold_change_matrix)
  for (mat in matrices) {
    rownames(mat) <- cn_names
    colnames(mat) <- celltype_names
  }
  
  # Calculate totals
  cn_totals <- rowSums(contingency_table)
  celltype_totals <- colSums(contingency_table)
  grand_total <- sum(contingency_table)
  
  # Test each CN-celltype combination
  for (i in 1:n_cn) {
    for (j in 1:n_celltype) {
      
      # Skip if cell type is too rare globally
      if (celltype_totals[j] < 5) {next}
      
      # Observed count
      observed <- contingency_table[i, j]
      
      # Expected count under null hypothesis of no enrichment
      expected <- (cn_totals[i] * celltype_totals[j]) / grand_total
      expected_matrix[i, j] <- expected
      
      # Skip if expected count is too low
      if (expected < 1) {next}
      
      # Calculate fold change (observed/expected)
      fold_change_matrix[i, j] <- observed / expected
      
      if (method == "fisher") {
        # Fisher's exact test approach
        # Create 2x2 contingency table:
        # [celltype_in_CN, other_celltypes_in_CN]
        # [celltype_in_other_CNs, other_celltypes_in_other_CNs]
        
        celltype_in_cn <- observed
        other_celltypes_in_cn <- cn_totals[i] - observed
        celltype_in_other_cns <- celltype_totals[j] - observed
        other_celltypes_in_other_cns <- grand_total - cn_totals[i] - celltype_totals[j] + observed
        
        # Create contingency matrix
        cont_matrix <- matrix(c(celltype_in_cn, other_celltypes_in_cn,
                                celltype_in_other_cns, other_celltypes_in_other_cns),
                              nrow = 2, byrow = TRUE)
        
        # Perform Fisher's exact test (one-sided for enrichment)
        fisher_result <- fisher.test(cont_matrix, alternative = "greater")
        pvalue_matrix[i, j] <- fisher_result$p.value
        
      } else if (method == "hypergeometric") {
        # Hypergeometric test
        # H0: celltype is randomly distributed across CNs
        # H1: celltype is enriched in this CN
        
        # Parameters for hypergeometric distribution:
        # Population size: grand_total
        # Number of "successes" in population: celltype_totals[j]
        # Sample size: cn_totals[i]
        # Observed successes: observed
        
        pvalue_matrix[i, j] <- phyper(observed - 1, celltype_totals[j], 
                                      grand_total - celltype_totals[j], 
                                      cn_totals[i], lower.tail = FALSE)
        
      } else if (method == "binomial") {
        # Binomial test
        # H0: probability of this celltype in this CN equals global probability
        # H1: probability is higher (enrichment)
        
        global_prob <- celltype_totals[j] / grand_total
        binom_result <- binom.test(observed, cn_totals[i], global_prob, alternative = "greater")
        pvalue_matrix[i, j] <- binom_result$p.value
      }
      
      # Determine enrichment (binary)
      enrichment_matrix[i, j] <- as.numeric(observed > expected)
    }
  }
  
  # Adjust p-values for multiple testing
  pvalue_vector <- as.vector(pvalue_matrix)
  pvalue_vector_clean <- pvalue_vector[!is.na(pvalue_vector)]
  
  if (length(pvalue_vector_clean) > 0) {
    pvalue_adj_vector <- p.adjust(pvalue_vector_clean, method = "BH")
    
    # Reconstruct adjusted p-value matrix
    pvalue_adj_matrix <- matrix(NA, nrow = n_cn, ncol = n_celltype)
    rownames(pvalue_adj_matrix) <- cn_names
    colnames(pvalue_adj_matrix) <- celltype_names
    
    idx <- 1
    for (i in 1:n_cn) {
      for (j in 1:n_celltype) {
        if (!is.na(pvalue_matrix[i, j])) {
          pvalue_adj_matrix[i, j] <- pvalue_adj_vector[idx]
          idx <- idx + 1
        }
      }
    }
  } else {
    pvalue_adj_matrix <- pvalue_matrix
  }
  
  print("Significance testing completed")
  
  return(pvalue_adj_matrix)
}

# ============================================================================
# CREATE SIGNIFICANCE ANNOTATION MATRIX
# ============================================================================

creat_sig_anno_matrix  <- function(pvalue_adj){
  print("Creating significance annotations...")
  
  # Get celltype and CN names
  cn_names <- rownames(pvalue_adj)
  celltype_names <- colnames(pvalue_adj)
  
  # Create significance annotation matrix for heatmap
  sig_matrix <- matrix("", nrow = length(cn_names), ncol = length(celltype_names))
  rownames(sig_matrix) <- cn_names
  colnames(sig_matrix) <- celltype_names
  
  # Add significance stars
  sig_matrix[!is.na(pvalue_adj) & pvalue_adj < 0.001] <- "***"
  sig_matrix[!is.na(pvalue_adj) & pvalue_adj >= 0.001 & pvalue_adj < 0.01] <- "**"
  sig_matrix[!is.na(pvalue_adj) & pvalue_adj >= 0.01 & pvalue_adj < 0.05] <- "*"
  
  # Count significant enrichments
  n_sig_001 <- sum(sig_matrix == "***", na.rm = TRUE)
  n_sig_01 <- sum(sig_matrix == "**", na.rm = TRUE)
  n_sig_05 <- sum(sig_matrix == "*", na.rm = TRUE)
  
  print(paste("Significant enrichments found:"))
  print(paste("  p < 0.001:", n_sig_001))
  print(paste("  p < 0.01:", n_sig_01))
  print(paste("  p < 0.05:", n_sig_05))
  
  return(sig_matrix)
}

# ============================================================================
# Function to perform CCA between two CNs
# ============================================================================
perform_cca_analysis <- function(data, cn1, cn2, cn_column, feature_cols, rfs_group = NULL, n_permutations = 5000) {
  
  # 1. --- Data Filtering and Preparation ---
  # Filter data for the two CNs of interest
  if (!is.null(rfs_group)) {
    analysis_data <- data[data[[cn_column]] %in% c(cn1, cn2) & data$RFS_status == rfs_group, ]
  } else {
    analysis_data <- data[data[[cn_column]] %in% c(cn1, cn2), ]
  }
  
  if (nrow(analysis_data) == 0) return(NULL)
  
  # Split by CN
  cn1_data <- analysis_data[analysis_data[[cn_column]] == cn1, ]
  cn2_data <- analysis_data[analysis_data[[cn_column]] == cn2, ]
  
  # Find common patients to ensure data is paired
  common_patients <- intersect(cn1_data$patient_id, cn2_data$patient_id)
  
  # Need at least 3 common patients to proceed
  if (length(common_patients) < 3) return(NULL)
  
  # Filter data for common patients
  cn1_data <- cn1_data[cn1_data$patient_id %in% common_patients, ]
  cn2_data <- cn2_data[cn2_data$patient_id %in% common_patients, ]
  
  # 2. --- Log Transformation (NEW STEP) ---
  # Find the smallest non-zero frequency to use as a pseudocount
  pseudocount <- min(analysis_data[,feature_cols]) / 2
  
  # Apply log transformation to the feature columns for both datasets
  cn1_data <- cn1_data %>%
    mutate(across(all_of(feature_cols), ~ log(.x + pseudocount)))

  cn2_data <- cn2_data %>%
    mutate(across(all_of(feature_cols), ~ log(.x + pseudocount)))
  
  # cn1_data <- cn1_data %>%
  #   mutate(across(all_of(feature_cols), ~ .x + pseudocount)) %>%  # Add pseudocount first
  #   rowwise() %>%
  #   mutate(
  #     # Calculate geometric mean for this row
  #     geom_mean = exp(mean(log(c_across(all_of(feature_cols))))),
  #     # Apply CLR transformation to all feature columns
  #     across(all_of(feature_cols), ~ log(.x / geom_mean))
  #   ) %>%
  #   select(-geom_mean) %>%  # Remove the temporary geometric mean column
  #   ungroup()
  # 
  # cn2_data <- cn2_data %>%
  #   mutate(across(all_of(feature_cols), ~ .x + pseudocount)) %>%  # Add pseudocount first
  #   rowwise() %>%
  #   mutate(
  #     # Calculate geometric mean for this row
  #     geom_mean = exp(mean(log(c_across(all_of(feature_cols))))),
  #     # Apply CLR transformation to all feature columns
  #     across(all_of(feature_cols), ~ log(.x / geom_mean))
  #   ) %>%
  #   select(-geom_mean) %>%  # Remove the temporary geometric mean column
  #   ungroup()
  
  # 3. --- Prepare Matrices and Calculate Observed Correlation ---
  # Arrange by sample_id to ensure rows are perfectly aligned
  X <- cn1_data %>% arrange(patient_id) %>% select(all_of(feature_cols)) %>% as.matrix()
  Y <- cn2_data %>% arrange(patient_id) %>% select(all_of(feature_cols)) %>% as.matrix()
  
  # Remove zero-variance columns that may arise after filtering
  X_var <- apply(X, 2, var, na.rm = TRUE)
  Y_var <- apply(Y, 2, var, na.rm = TRUE)
  X <- X[, !is.na(X_var) & X_var > 0, drop = FALSE]
  Y <- Y[, !is.na(Y_var) & Y_var > 0, drop = FALSE]
  
  if (ncol(X) == 0 || ncol(Y) == 0) return(NULL)
  
  # Calculate the TRUE observed canonical correlation
  observed_cca <- tryCatch(cancor(X, Y), error = function(e) NULL)
  if (is.null(observed_cca)) return(NULL)
  observed_correlation <- observed_cca$cor[1]
  
  # 4. --- Permutation Test (NEW STEP) ---
  perm_correlations <- numeric(n_permutations)
  
  for (i in 1:n_permutations) {
    # Shuffle the rows of one matrix to break the patient-to-patient pairing
    Y_shuffled <- Y[sample(nrow(Y)), , drop = FALSE]
    
    # Calculate CCA on the permuted data
    perm_cca_result <- tryCatch(cancor(X, Y_shuffled), error = function(e) NULL)
    
    if (!is.null(perm_cca_result)) {
      perm_correlations[i] <- perm_cca_result$cor[1]
    } else {
      perm_correlations[i] <- NA # Store NA if CCA fails on a permutation
    }
  }
  
  # 5. --- Calculate P-Value ---
  # P-value is the proportion of permutations with a correlation >= observed
  # Add 1 to numerator and denominator to be conservative and avoid p=0
  p_value <- (sum(perm_correlations >= observed_correlation, na.rm = TRUE) + 1) / (n_permutations + 1)
  
  # Return the final results
  return(list(
    observed_correlation = observed_correlation,
    p_value = p_value,
    n_patients = length(common_patients),
    n_permutations = n_permutations
  ))
}

# ============================================================================
# Function to create correlation matrix from CCA results
# ============================================================================
create_cca_correlation_matrix <- function(cca_results_group, cn_list, p.adjust = TRUE) {
  n_cn <- length(cn_list)
  
  # Initialize correlation matrix
  cor_matrix <- matrix(0, nrow = n_cn, ncol = n_cn)
  rownames(cor_matrix) <- paste0("CN", cn_list)
  colnames(cor_matrix) <- paste0("CN", cn_list)
  diag(cor_matrix) <- 1  # Perfect correlation with itself
  
  # Also create p-value and sample size matrices for additional information
  pval_matrix <- matrix(1, nrow = n_cn, ncol = n_cn)  # Initialize with 1s
  rownames(pval_matrix) <- paste0("CN", cn_list)
  colnames(pval_matrix) <- paste0("CN", cn_list)
  diag(pval_matrix) <- 0  # p-value for self-correlation is 0
  
  n_patients_matrix <- matrix(0, nrow = n_cn, ncol = n_cn)
  rownames(n_patients_matrix) <- paste0("CN", cn_list)
  colnames(n_patients_matrix) <- paste0("CN", cn_list)
  
  # Fill in the matrices
  for (pair_name in names(cca_results_group)) {
    # Parse pair name - handle different possible formats
    if (grepl("_vs_", pair_name)) {
      # Format: "CN1_vs_CN2"
      cns <- gsub("CN", "", strsplit(pair_name, "_vs_")[[1]])
    } else {
      # Format: "CN1_CN2" (assuming underscore separated)
      parts <- strsplit(pair_name, "_")[[1]]
      # Find parts that start with "CN"
      cn_parts <- parts[grepl("^CN", parts)]
      if (length(cn_parts) >= 2) {
        cns <- gsub("CN", "", cn_parts[1:2])
      } else {
        # Try extracting numbers directly
        cns <- gsub("CN", "", parts)
        cns <- cns[cns != ""]
      }
    }
    
    if (length(cns) < 2) {
      warning(paste("Could not parse CN pair from:", pair_name))
      next
    }
    
    cn1 <- cns[1]
    cn2 <- cns[2]
    
    # Find indices in cn_list
    cn1_idx <- which(cn_list == cn1)
    cn2_idx <- which(cn_list == cn2)
    
    if (length(cn1_idx) > 0 && length(cn2_idx) > 0) {
      # Access the correct field name
      result <- cca_results_group[[pair_name]]
      
      # Handle nested structure if needed
      if (is.list(result) && length(result) == 1 && is.list(result[[1]])) {
        result <- result[[1]]
      }
      
      # Extract values with error checking
      if ("observed_correlation" %in% names(result)) {
        cor_val <- result$observed_correlation
        cor_matrix[cn1_idx, cn2_idx] <- cor_val
        cor_matrix[cn2_idx, cn1_idx] <- cor_val
      }
      
      if ("p_value" %in% names(result)) {
        pval <- result$p_value
        pval_matrix[cn1_idx, cn2_idx] <- pval
        pval_matrix[cn2_idx, cn1_idx] <- pval
      
      }
    } else {
      warning(paste("CN indices not found for:", cn1, "and", cn2))
    }
  }
  
  if(p.adjust){
    names <- colnames(pval_matrix)
    pval_matrix <- matrix(data = p.adjust(pval_matrix,method = "BH"), nrow = nrow(pval_matrix), ncol = ncol(pval_matrix))
    colnames(pval_matrix) <- names
    rownames(pval_matrix) <- names
  }
  
  return(list(
    correlation_matrix = cor_matrix,
    pvalue_matrix = pval_matrix
  ))
}

# ============================================================================
# Function to create network from correlation matrix
# ============================================================================
create_cca_network <- function(matrices_data, correlation_threshold = 0.3, 
                               p_value_threshold = 0.05, title = "", 
                               node_colors = NULL, show_all_edges = FALSE) {
  
  # Handle different input formats
  if (is.matrix(matrices_data)) {
    # Old format - just correlation matrix
    cor_matrix <- matrices_data
    pval_matrix <- NULL
  } else if (is.list(matrices_data)) {
    # New format - list with correlation_matrix and pvalue_matrix
    cor_matrix <- matrices_data$correlation_matrix
    pval_matrix <- matrices_data$pvalue_matrix
  } else {
    stop("Input must be either a correlation matrix or a list with correlation_matrix and pvalue_matrix")
  }
  
  if (is.null(cor_matrix)) {
    stop("Correlation matrix not found in input data")
  }
  
  # Create adjacency matrix
  adj_matrix <- cor_matrix
  
  # Apply correlation threshold
  adj_matrix[abs(adj_matrix) < correlation_threshold] <- 0
  
  # Apply p-value threshold if available
  if (!is.null(pval_matrix) && !show_all_edges) {
    # Only keep edges that are both above correlation threshold AND significant
    adj_matrix[pval_matrix >= p_value_threshold] <- 0
  }
  
  # Remove diagonal
  diag(adj_matrix) <- 0
  
  # Create igraph object
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
  
  if (vcount(g) == 0) {
    warning("No edges meet the specified thresholds")
    return(NULL)
  }
  
  # Set vertex attributes
  V(g)$name <- rownames(cor_matrix)
  V(g)$label <- V(g)$name
  
  # Default node coloring if not provided
  if (is.null(node_colors)) {
    # Create a color scheme based on CN numbers
    cn_numbers <- as.numeric(gsub("CN", "", V(g)$name))
    V(g)$color <- rainbow(length(unique(cn_numbers)))[as.numeric(factor(cn_numbers))]
  } else {
    V(g)$color <- node_colors
  }
  
  # Set edge attributes with enhanced styling
  E(g)$weight_abs <- abs(E(g)$weight)
  E(g)$color <- ifelse(E(g)$weight > 0, "darkgrey")  # Red for positive, blue for negative
  E(g)$width <- pmax(E(g)$weight_abs * 4, 0.5)  # Scale width, minimum 0.5
  
  # Add significance information to edges if p-values available
  if (!is.null(pval_matrix)) {
    # Get p-values for each edge
    edge_pvals <- numeric(ecount(g))
    for (i in 1:ecount(g)) {
      edge <- ends(g, i)
      v1_name <- edge[1]
      v2_name <- edge[2]
      edge_pvals[i] <- pval_matrix[v1_name, v2_name]
    }
    E(g)$p_value <- edge_pvals
    E(g)$significant <- edge_pvals < p_value_threshold
    
    # Style edges based on significance
    E(g)$alpha <- ifelse(E(g)$significant, 0.9, 0.4)
    E(g)$line_type <- ifelse(E(g)$significant, "solid", "dashed")
  } else {
    E(g)$alpha <- 0.7
    E(g)$line_type <- "solid"
  }
  
  # Create layout
  layout <- layout_with_fr(g, niter = 50)

  
  # Create the plot
  p <- g %>%
    as_tbl_graph() %>%
    ggraph(layout = layout) +
    geom_edge_link(aes(width = weight_abs, 
                       color = weight, 
                       alpha = alpha),
                   show.legend = FALSE) +
    geom_node_point(aes(color = I(color)), size = 15, alpha = 0.9) +
    geom_node_text(aes(label = label), color = "white", fontface = "bold", size = 4) +
    scale_edge_width(range = c(0.5, 4)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom",
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # Create subtitle with threshold information
  subtitle_text <- paste("Correlation ≥", correlation_threshold)
  if (!is.null(pval_matrix) && !show_all_edges) {
    subtitle_text <- paste(subtitle_text, "& p-value <", p_value_threshold)
  }
  subtitle_text <- paste(subtitle_text, "| Edges:", ecount(g), "of", 
                         choose(vcount(g), 2), "possible")
  
  p <- p + labs(
    title = title,
    subtitle = subtitle_text
  )
  
  return(p)
}

# =============================================================================
## Interaction analysis
# =============================================================================

# =============================================================================
# Create interaction matrices for each tissue
# =============================================================================
create_interaction_matrix <- function(data, tissue_type, value_col = "mean_sigval") {
  tissue_data <- data %>% filter(Tissue == tissue_type)
  
  # Get all unique cell types
  all_celltypes <- unique(c(tissue_data$from_label, tissue_data$to_label))
  
  # Create matrix
  interaction_mat <- matrix(0, nrow = length(all_celltypes), ncol = length(all_celltypes))
  rownames(interaction_mat) <- all_celltypes
  colnames(interaction_mat) <- all_celltypes
  
  # Fill matrix
  for(i in 1:nrow(tissue_data)) {
    from_idx <- which(all_celltypes == tissue_data$from_label[i])
    to_idx <- which(all_celltypes == tissue_data$to_label[i])
    interaction_mat[from_idx, to_idx] <- tissue_data[[value_col]][i]
  }
  
  return(interaction_mat)
}

# =============================================================================
# Enhanced function to create appealing combined heatmap
# =============================================================================
create_combined_interaction_heatmap <- function(spe, tc_mat, im_mat, pt_mat) {
  
  # Ensure all matrices have the same dimensions by getting union of all cell types
  all_celltypes <- unique(c(rownames(tc_mat), rownames(im_mat), rownames(pt_mat),
                            colnames(tc_mat), colnames(im_mat), colnames(pt_mat)))
  all_celltypes <- sort(all_celltypes)
  
  # Function to standardize matrix dimensions
  standardize_matrix <- function(mat, all_types) {
    std_mat <- matrix(0, nrow = length(all_types), ncol = length(all_types))
    rownames(std_mat) <- all_types
    colnames(std_mat) <- all_types
    
    # Fill in existing values
    common_rows <- intersect(rownames(mat), all_types)
    common_cols <- intersect(colnames(mat), all_types)
    std_mat[common_rows, common_cols] <- mat[common_rows, common_cols]
    
    return(std_mat)
  }
  
  # Standardize all matrices
  tc_std <- standardize_matrix(tc_matrix, all_celltypes)
  im_std <- standardize_matrix(im_matrix, all_celltypes)
  pt_std <- standardize_matrix(pt_matrix, all_celltypes)
  
  # Create cell type annotations
  cell_type_categories <- data.frame(
    CellType = all_celltypes,
    Category = case_when(
      str_detect(all_celltypes, "^EC_") ~ "Epithelial",
      all_celltypes %in% c("B") ~ "B_cells",
      all_celltypes %in% c("CD4T", "CD8T", "Treg","Other_Immune") ~ "T_cells",
      all_celltypes %in% c("NK") ~ "NK_cells",
      str_detect(all_celltypes, "^DC") ~ "DCs",
      str_detect(all_celltypes, "^Mono_") ~ "Monocytes",
      str_detect(all_celltypes, "^Macro_") ~ "Macrophages",
      str_detect(all_celltypes, "^SC_|^CAF") ~ "Stromal",
      TRUE ~ "Other"
    )
  )
  
  # Color scheme for categories
  category_colors <- metadata(spe)$color_vectors$major_celltype
  
  # Row and column annotations
  row_annotation <- rowAnnotation(
    Category = cell_type_categories$Category,
    col = list(Category = category_colors),
    annotation_name_gp = gpar(fontsize = 10),
    simple_anno_size = unit(0.3, "cm")
  )
  
  col_annotation <- HeatmapAnnotation(
    Category = cell_type_categories$Category,
    col = list(Category = category_colors),
    annotation_name_gp = gpar(fontsize = 10),
    simple_anno_size = unit(0.3, "cm")
  )
  
  # Enhanced color scheme with better contrast
  col_fun <- colorRamp2(
    c(-1, -0.5, 0, 0.5, 1), 
    c("#053061", "#67a9cf", "white", "#ef8a62", "#67001f")
  )
  
  # Create individual heatmaps with consistent styling
  h1 <- Heatmap(tc_std,
                name = "TC",
                col = col_fun,
                column_title = "Tumor Core (TC)",
                column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize = 7),
                column_names_gp = gpar(fontsize = 7),
                row_names_rot = 0,
                column_names_rot = 90,
                top_annotation = col_annotation,
                border = TRUE,
                rect_gp = gpar(col = "white", lwd = 0.5))
  
  h2 <- Heatmap(im_std,
                name = "IM", 
                col = col_fun,
                column_title = "Invasive Margin (IM)",
                column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = FALSE,  # Hide row names for middle heatmap
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 7),
                column_names_rot = 90,
                top_annotation = col_annotation,
                border = TRUE,
                rect_gp = gpar(col = "white", lwd = 0.5))
  
  h3 <- Heatmap(pt_std,
                name = "PT",
                col = col_fun,
                column_title = "Peri-Tumor (PT)", 
                column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = FALSE,  # Hide row names for right heatmap
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 7),
                column_names_rot = 90,
                top_annotation = col_annotation,
                border = TRUE,
                rect_gp = gpar(col = "white", lwd = 0.5),
                heatmap_legend_param = list(
                  title = "Interaction\nSignal",
                  title_gp = gpar(fontsize = 11, fontface = "bold"),
                  labels_gp = gpar(fontsize = 9),
                  legend_height = unit(4, "cm"),
                  at = c(-1, -0.5, 0, 0.5, 1),
                  labels = c("Strong\nAvoidance", "Weak\nAvoidance", "No\nInteraction", 
                             "Weak\nAttraction", "Strong\nAttraction")
                ))
  
  # Combine heatmaps horizontally
  combined_heatmap <- h1 + h2 + h3
  
  # Add row annotation only to the first heatmap
  final_heatmap <- row_annotation + combined_heatmap
  
  return(final_heatmap)
}

# =============================================================================
# Function to compare interaction counts between groups
# =============================================================================
compare_interactions <- function(data, group_var, cell_from, cell_to) {
  
  comparison_data <- data %>%
    filter(from_label == cell_from, to_label == cell_to, !is.na(!!sym(group_var))) %>%
    select(group_by, ct, !!sym(group_var))
  
  # Summary statistics
  summary_stats <- comparison_data %>%
    group_by(!!sym(group_var)) %>%
    summarise(
      n = n(),
      mean_ct = mean(ct, na.rm = TRUE),
      median_ct = median(ct, na.rm = TRUE),
      sd_ct = sd(ct, na.rm = TRUE),
      se_ct = sd_ct / sqrt(n),
      .groups = "drop"
    )
  
  # Statistical test (Wilcoxon rank-sum test for non-parametric comparison)
  if(length(unique(comparison_data[[group_var]])) == 2) {
    groups <- unique(comparison_data[[group_var]])
    group1_data <- comparison_data[comparison_data[[group_var]] == groups[1], "ct"]
    group2_data <- comparison_data[comparison_data[[group_var]] == groups[2], "ct"]
    
    if(sum(!is.na(group2_data)) == 0 | sum(!is.na(group1_data)) == 0){
      test_result <- NULL
      p_value <- 1
    }
    else{
      test_result <- wilcox.test(group1_data, group2_data)
      p_value <- test_result$p.value
    }
  }
  
  return(list(
    summary = summary_stats,
    p_value = p_value,
    test_result = test_result,
    raw_data = comparison_data
  ))
}

# =============================================================================
# Function for tissue-specific and clinically-focused systematic comparison
# =============================================================================
focused_systematic_comparison <- function(data, group_var, tissue_type, min_samples = 3) {
  
  # Filter for specific tissue
  tissue_data <- data %>% filter(Tissue == tissue_type, !is.na(!!sym(group_var)))
  
  # Define clinically relevant interaction pairs
  clinical_pairs <- bind_rows(
    # Tumor-Stromal interactions (barrier mechanisms)
    expand_grid(from = epithelial_cells, to = stromal_cells) %>%
      mutate(interaction_type = "Tumor_Stromal"),
    expand_grid(from = stromal_cells, to = epithelial_cells) %>%
      mutate(interaction_type = "Stromal_Tumor"),
    
    # Tumor-Immune interactions (immunotherapy potential)  
    expand_grid(from = epithelial_cells, to = immune_cells) %>%
      mutate(interaction_type = "Tumor_Immune"),
    expand_grid(from = immune_cells, to = epithelial_cells) %>%
      mutate(interaction_type = "Immune_Tumor"),
    
    # # Stromal-Immune interactions (microenvironment modulation)
    # expand_grid(from = stromal_cells, to = immune_cells) %>%
    #   mutate(interaction_type = "Stromal_Immune"),
    # expand_grid(from = immune_cells, to = stromal_cells) %>%
    #   mutate(interaction_type = "Immune_Stromal")
  )
  
  # Check which pairs have sufficient data
  valid_pairs <- tissue_data %>%
    group_by(from_label, to_label) %>%
    summarise(n_samples = n(), .groups = "drop") %>%
    filter(n_samples >= min_samples) %>%
    inner_join(clinical_pairs, by = c("from_label" = "from", "to_label" = "to"))
  
  if(nrow(valid_pairs) == 0) {
    cat("No valid clinical interaction pairs found for", tissue_type, "\n")
    return(NULL)
  }
  
  results <- list()
  
  for(i in 1:nrow(valid_pairs)) {
    from_cell <- valid_pairs$from_label[i]
    to_cell <- valid_pairs$to_label[i]
    int_type <- valid_pairs$interaction_type[i]
    
    comparison <- compare_interactions(tissue_data, group_var, from_cell, to_cell)
    
    results[[paste(from_cell, to_cell, sep = "_")]] <- list(
      from = from_cell,
      to = to_cell,
      interaction_type = int_type,
      tissue = tissue_type,
      p_value = comparison$p_value,
      summary = comparison$summary
    )
  }
  
  # Extract p-values and apply multiple testing correction
  p_values <- sapply(results, function(x) x$p_value)
  adj_p_values <- p.adjust(p_values, method = "BH")
  
  # Create results summary
  results_df <- data.frame(
    tissue = tissue_type,
    from_label = sapply(results, function(x) x$from),
    to_label = sapply(results, function(x) x$to),
    interaction_type = sapply(results, function(x) x$interaction_type),
    p_value = p_values,
    adj_p_value = adj_p_values,
    significant = adj_p_values < 0.05
  ) %>%
    arrange(p_value)
  
  return(list(
    results_table = results_df,
    detailed_results = results
  ))
}

# =============================================================================
# Function for treatment-stratified systematic comparison
# =============================================================================
treatment_stratified_comparison <- function(data, treatment_type, min_samples = 3) {
  
  # Filter for specific treatment (pooling across all tissues)
  treatment_data <- data %>% 
    filter(Treatment == treatment_type, !is.na(RFS_status)) %>%
    # Add tissue information to track which tissues contribute to each interaction
    mutate(interaction_id = paste(from_label, to_label, sep = "_"))
  
  cat("Analyzing", treatment_type, "treatment group...\n")
  cat("Total samples:", length(unique(treatment_data$group_by)), "\n")
  cat("RFS status distribution:\n")
  print(table(treatment_data$RFS_status))
  cat("\n")
  
  # Define clinically relevant interaction pairs
  clinical_pairs <- bind_rows(
    # Tumor-Stromal interactions (barrier mechanisms)
    expand_grid(from = epithelial_cells, to = stromal_cells) %>%
      mutate(interaction_type = "Tumor_Stromal"),
    expand_grid(from = stromal_cells, to = epithelial_cells) %>%
      mutate(interaction_type = "Stromal_Tumor"),
    
    # Tumor-Immune interactions (immunotherapy potential)  
    expand_grid(from = epithelial_cells, to = immune_cells) %>%
      mutate(interaction_type = "Tumor_Immune"),
    expand_grid(from = immune_cells, to = epithelial_cells) %>%
      mutate(interaction_type = "Immune_Tumor"),
    
    # Stromal-Immune interactions (microenvironment modulation)
    expand_grid(from = stromal_cells, to = immune_cells) %>%
      mutate(interaction_type = "Stromal_Immune"),
    expand_grid(from = immune_cells, to = stromal_cells) %>%
      mutate(interaction_type = "Immune_Stromal")
  )
  
  # Check which pairs have sufficient data
  valid_pairs <- treatment_data %>%
    group_by(from_label, to_label) %>%
    summarise(
      n_samples = n(),
      n_tissues = length(unique(Tissue)),
      tissues = paste(unique(Tissue), collapse = ","),
      rfs_0 = sum(RFS_status == 0),
      rfs_1 = sum(RFS_status == 1),
      .groups = "drop"
    ) %>%
    filter(n_samples >= min_samples, rfs_0 >= 1, rfs_1 >= 1) %>%  # Need both RFS groups
    inner_join(clinical_pairs, by = c("from_label" = "from", "to_label" = "to"))
  
  cat("Valid clinical interaction pairs:", nrow(valid_pairs), "\n")
  
  if(nrow(valid_pairs) == 0) {
    cat("No valid clinical interaction pairs found for", treatment_type, "\n")
    return(NULL)
  }
  
  results <- list()
  
  for(i in 1:nrow(valid_pairs)) {
    from_cell <- valid_pairs$from_label[i]
    to_cell <- valid_pairs$to_label[i]
    int_type <- valid_pairs$interaction_type[i]
    
    comparison <- compare_interactions(treatment_data, "RFS_status", from_cell, to_cell)
    
    results[[paste(from_cell, to_cell, sep = "_")]] <- list(
      from = from_cell,
      to = to_cell,
      interaction_type = int_type,
      treatment = treatment_type,
      tissues_present = valid_pairs$tissues[i],
      n_tissues = valid_pairs$n_tissues[i],
      p_value = comparison$p_value,
      summary = comparison$summary
    )
  }
  
  # Extract p-values and apply multiple testing correction
  p_values <- sapply(results, function(x) x$p_value)
  adj_p_values <- p.adjust(p_values, method = "BH")
  
  # Create results summary
  results_df <- data.frame(
    treatment = treatment_type,
    from_label = sapply(results, function(x) x$from),
    to_label = sapply(results, function(x) x$to),
    interaction_type = sapply(results, function(x) x$interaction_type),
    tissues_present = sapply(results, function(x) x$tissues_present),
    n_tissues = sapply(results, function(x) x$n_tissues),
    p_value = p_values,
    adj_p_value = adj_p_values,
    significant = adj_p_values < 0.05
  ) %>%
    arrange(p_value)
  
  return(list(
    results_table = results_df,
    detailed_results = results
  ))
}

# =============================================================================
# EC interaction analysis
# =============================================================================
ec_interaction_analysis <- function(data) {
  
  # Filter for IM regions only (as specified)
  tc_im_data <- data %>%
    filter(Tissue %in% c("IM"), !is.na(mean_sigval))
  
  cat("Analyzing EC interaction patterns in IM regions...\n")
  cat("Total interactions:", nrow(tc_im_data), "\n")
  
  # Get all other cell types (excluding the two ECs of interest)
  all_celltypes <- unique(c(tc_im_data$from_label, tc_im_data$to_label))
  other_celltypes <- setdiff(all_celltypes, c("EC_GLUT1", "EC_EpCAM"))
  
  cat("Other cell types:", length(other_celltypes), "\n")
  
  # Function to extract undirected interactions for a specific EC type
  extract_ec_interactions <- function(ec_type) {
    
    # Get both directions: EC ↔ Other
    ec_interactions <- tc_im_data %>%
      filter(
        (from_label == ec_type & to_label %in% other_celltypes) |
          (from_label %in% other_celltypes & to_label == ec_type)
      ) %>%
      mutate(
        # Standardize interaction pairs (always put EC first for consistency)
        partner = ifelse(from_label == ec_type, to_label, from_label),
        interaction_pair = paste(ec_type, partner, sep = "_")
      )
    
    # Combine bidirectional interactions by averaging
    combined_interactions <- ec_interactions %>%
      group_by(partner) %>%
      summarise(
        # Average the interaction values from both directions
        mean_sigval = mean(mean_sigval, na.rm = TRUE),
        # Count how many samples contributed (from both directions)
        n_samples = n(),
        # Track if we have both directions
        n_directions = length(unique(paste(from_label, to_label, sep = "→"))),
        # Get the raw values for transparency
        raw_values = paste(round(mean_sigval, 3), collapse = ", "),
        .groups = "drop"
      ) %>%
      mutate(
        ec_focus = ec_type,
        abs_sigval = abs(mean_sigval),
        interaction_type = paste(ec_type, "↔", partner),  # Use ↔ for undirected
        interaction_strength = case_when(
          mean_sigval > 0.5 ~ "Strong Attraction",
          mean_sigval > 0 ~ "Weak Attraction", 
          mean_sigval > -0.5 ~ "Weak Avoidance",
          TRUE ~ "Strong Avoidance"
        ),
        data_completeness = case_when(
          n_directions == 2 ~ "Bidirectional",
          n_directions == 1 ~ "Unidirectional"
        )
      )
    
    return(combined_interactions)
  }
  
  # Extract interactions for both EC types
  glut1_interactions <- extract_ec_interactions("EC_GLUT1")
  epcam_interactions <- extract_ec_interactions("EC_EpCAM")
  
  return(list(
    glut1_data = glut1_interactions,
    epcam_data = epcam_interactions,
    combined_data = bind_rows(glut1_interactions, epcam_interactions)
  ))
}

# =============================================================================
# Function to create ordered bidirectional barplot
# =============================================================================
create_ec_barplot <- function(ec_data, ec_name, clinical_context) {
  
  # Prepare data with proper ordering
  plot_data <- ec_data %>%
    arrange(desc(mean_sigval)) %>%  # Order by actual signal value (decreasing)
    mutate(
      # Create clean labels
      partner_clean = str_replace_all(partner, "_", " "),
      # Color based on interaction strength and direction
      fill_color = case_when(
        mean_sigval > 0.3 ~ "Strong Attraction",
        mean_sigval > 0 ~ "Weak Attraction", 
        mean_sigval > -0.3 ~ "Weak Avoidance",
        TRUE ~ "Strong Avoidance"
      ),
      # Add data quality indicator
      label_text = paste0(sprintf("%.2f", mean_sigval), 
                          ifelse(data_completeness == "Unidirectional", "*", ""))
    )
  
  # Define colors for interaction strength
  colors <- c(
    "Strong Attraction" = "#d62728",     # Red
    "Weak Attraction" = "#ff7f0e",       # Orange  
    "Weak Avoidance" = "#2ca02c",        # Green
    "Strong Avoidance" = "#1f77b4"       # Blue
  )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = reorder(partner_clean, mean_sigval), 
                             y = mean_sigval, fill = fill_color)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = label_text,
                  y = mean_sigval + ifelse(mean_sigval > 0, 0.05, -0.05)),
              size = 2.5, color = "black") +
    labs(
      title = paste0(str_replace(ec_name, "_", " "), " Undirected Interaction Profile"),
      subtitle = paste0("Clinical Context: ", clinical_context, " | IM Regions | *Unidirectional data"),
      x = "Interaction Partner Cell Type",
      y = "Mean Interaction Signal (-1: Avoidance, +1: Attraction)",
      fill = "Interaction Strength"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "darkblue"),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.title = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    ) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(-1.1, 1.1), breaks = seq(-1, 1, 0.5))
  
  return(p)
}

# ============================================================================
# CHORD PLOT FOR EC INTERACTIONS
# ============================================================================

create_ec_chord_plot <- function(interaction_data, tissue_type = "IM", 
                                 min_prop_significant = 0.3, ec_cells = NULL, selected_targets = NULL, savePath) {
  
  # Default target selection if not specified
  if(is.null(selected_targets)) {
    selected_targets <- c(stromal_cells, immune_cells[1:8])  # Select key immune cells
  }
  
  # Filter data for EC interactions with selected targets
  chord_data <- interaction_data %>%
    filter(
      Tissue == tissue_type,
      from_label %in% ec_cells,
      to_label %in% selected_targets,
      prop_significant >= min_prop_significant  # Filter for meaningful interactions
    ) %>%
    select(from_label, to_label, mean_sigval, prop_significant, mean_ct) %>%
    mutate(
      # Use interaction strength as weight
      interaction_weight = abs(mean_sigval) * mean_ct * 10,  # Scale for visualization
      interaction_direction = ifelse(mean_sigval > 0, "Attractive", "Avoidance")
    )
  
  if(nrow(chord_data) == 0) {
    cat("No significant EC interactions found for", tissue_type, "\n")
    return(NULL)
  }
  
  # Create adjacency matrix
  all_cells <- unique(c(chord_data$from_label, chord_data$to_label))
  adj_matrix <- matrix(0, nrow = length(all_cells), ncol = length(all_cells))
  rownames(adj_matrix) <- all_cells
  colnames(adj_matrix) <- all_cells
  
  # Fill matrix with interaction weights
  for(i in 1:nrow(chord_data)) {
    from_idx <- which(all_cells == chord_data$from_label[i])
    to_idx <- which(all_cells == chord_data$to_label[i])
    adj_matrix[from_idx, to_idx] <- chord_data$interaction_weight[i]
  }
  
  # Color scheme
  cell_colors <- metadata(spe)$color_vectors$sub_celltype
  
  # Set colors for present cells
  present_cells <- intersect(names(cell_colors), all_cells)
  colors <- cell_colors[present_cells]
  
  # Create chord diagram
  circos.clear()
  
  pdf(savePath,width = 12,height = 9)
  circos.par(start.degree = 90, gap.degree = 4)
  chordDiagram(t(adj_matrix),
               grid.col = colors,
               transparency = 0.3,
               directional = 1,
               direction.type = "arrows",
               link.arr.type = "big.arrow",
               annotationTrack = "grid",
               preAllocateTracks = 1,
               link.border = "grey40")
  
  # Add sector labels
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    
    # Color-code labels
    label_color <- if(sector.name %in% ec_cells) "darkred" else
      if(sector.name %in% immune_cells) "darkblue" else "darkorange"
    
    circos.text(mean(xlim), ylim[1] + 0.2, sector.name, 
                facing = "clockwise", niceFacing = TRUE, 
                adj = c(0, 0.5), cex = 0.8, col = label_color, font = 2)
  }, bg.border = NA)
  
  title(paste("EC Cell Interaction Network -", tissue_type), cex.main = 1.3)
  
  # Add legend
  legend("bottomright", 
         legend = c("EC Cells", "Immune Cells", "Stromal Cells"),
         col = c("darkred", "darkblue", "darkorange"),
         pch = 15, cex = 0.9, bty = "n")
  dev.off()
  circos.clear()
  
  return(chord_data)
}

# ============================================================================
# BARPLOT WITH ERROR BARS FOR GROUP COMPARISONS
# ============================================================================
create_interaction_comparison_barplot <- function(raw_data, ec_cell, target_cell, 
                                                  tissue_type = "IM", comparison_var = "RFS_status") {
  
  # Convert comparison variable to factor if needed
  raw_data <- raw_data %>%
    mutate(!!sym(comparison_var) := as.factor(!!sym(comparison_var)))
  
  # Filter and prepare data
  filtered_data <- raw_data %>%
    filter(
      Tissue == tissue_type,
      from_label == ec_cell,
      to_label == target_cell
    )
  
  # Check if we have data
  if(nrow(filtered_data) == 0) {
    warning(paste("No data found for", ec_cell, "->", target_cell, "in", tissue_type))
    return(NULL)
  }
  
  # Calculate summary statistics
  plot_data <- filtered_data %>%
    group_by(!!sym(comparison_var)) %>%
    summarise(
      mean_interaction = mean(sigval, na.rm = TRUE),
      se_interaction = sd(sigval, na.rm = TRUE) / sqrt(n()),
      mean_count = mean(ct, na.rm = TRUE),
      se_count = sd(ct, na.rm = TRUE) / sqrt(n()),
      prop_significant = mean(sig, na.rm = TRUE),
      prop_attractive = mean(interaction, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    # Handle cases with small sample sizes
    mutate(
      se_interaction = ifelse(is.na(se_interaction) | n_samples < 2, 0, se_interaction),
      se_count = ifelse(is.na(se_count) | n_samples < 2, 0, se_count)
    )
  
  # Perform statistical test between groups
  if(length(unique(plot_data[[comparison_var]])) == 2) {
    # wilcox_test for count
    group1_ct <- filtered_data %>% filter(!!sym(comparison_var) == unique(plot_data[[comparison_var]])[1]) %>% pull(ct)
    group2_ct <- filtered_data %>% filter(!!sym(comparison_var) == unique(plot_data[[comparison_var]])[2]) %>% pull(ct)
    
    if(length(group1_ct) > 1 && length(group2_ct) > 1) {
      wilcox_test_result_ct <- wilcox.test(group1_ct, group2_ct)
      p_value_count <- wilcox_test_result_ct$p.value
    } else {
      p_value_count <- NA
    }
  } else {
    p_value_interaction <- NA
    p_value_count <- NA
  }
  
  # Format p-values for display
  p_text_count <- if(!is.na(p_value_count)) {
    if(p_value_count < 0.001) "p < 0.001"
    else paste0("p = ", round(p_value_count, 3))
  } else "p = N/A"
  
  # Create barplot for interaction count
  p2 <- ggplot(plot_data, aes(x = !!sym(comparison_var), y = mean_count, 
                              fill = !!sym(comparison_var))) +
    geom_bar(stat = "identity", alpha = 0.7, width = 0.6) +
    geom_errorbar(aes(ymin = pmax(0, mean_count - se_count), 
                      ymax = mean_count + se_count), 
                  width = 0.2, color = "black") +
    # Add sample size labels on bars
    geom_text(aes(label = paste0("n=", n_samples)), 
              position = position_dodge(width = 0.6), 
              vjust = -0.5, size = 3, color = "black") +
    labs(
      title = paste0(ec_cell, " → ", target_cell, " (Interaction Count)"),
      subtitle = paste0(comparison_var, " comparison (", tissue_type, ") | ", p_text_count),
      x = comparison_var,
      y = "Mean Interaction Count",
      fill = comparison_var
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9, color = ifelse(!is.na(p_value_count) && p_value_count < 0.05, "red", "black")),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      legend.position = "none"
    ) +
    scale_fill_manual(values = RFS_color) +
    scale_y_continuous(limits = c(0, max(plot_data$mean_count + plot_data$se_count) * 1.1))
  
  # Create summary statistics table
  summary_stats <- plot_data %>%
    mutate(
      p_value_interaction = p_value_interaction,
      p_value_count = p_value_count
    )
  
  return(list(
    interaction_count = p2, 
    data = summary_stats,
    raw_data = filtered_data,
    statistical_tests = list(
      count_p = p_value_count
    )
  ))
}

# Function to create multiple comparison plots at once
create_multiple_ec_comparisons <- function(raw_data, ec_cells, target_cells, 
                                           tissue_type = "IM", comparison_var = "RFS_status") {
  
  all_plots <- list()
  significant_results <- data.frame()
  
  for(ec_cell in ec_cells) {
    for(target_cell in target_cells) {
      tryCatch({
        result <- create_interaction_comparison_barplot(raw_data, ec_cell, target_cell, 
                                                        tissue_type, comparison_var)
        
        if(!is.null(result)) {
          plot_name <- paste0(ec_cell, "_", target_cell)
          all_plots[[paste0(plot_name, "_count")]] <- result$interaction_count
          
          # Collect significant results
          if(!is.na(result$statistical_tests$count_p) && 
             result$statistical_tests$count_p < 0.05) {
            
            significant_results <- rbind(significant_results, 
                                         data.frame(
                                           ec_cell = ec_cell,
                                           target_cell = target_cell,
                                           comparison = comparison_var,
                                           tissue = tissue_type,
                                           p_value_count = result$statistical_tests$count_p,
                                           stringsAsFactors = FALSE
                                         ))
          }
        }
      }, error = function(e) {
        cat("Error processing", ec_cell, "->", target_cell, ":", e$message, "\n")
      })
    }
  }
  
  return(list(
    plots = all_plots,
    significant_results = significant_results
  ))
}

# =============================================================================
## Patch analysis
# =============================================================================
analysis_patch_characteristic <- function(data,patch_column_name,patche_ids,patch_type){
  # Initialize results dataframe
  patch_size_results <- data.frame()
  
  for (patch_id in patche_ids) {
    patch_cells_subset <- data[data[`patch_column_name`] == patch_id, ]
    
      # Basic metrics
      n_cells <- nrow(patch_cells_subset)
      sample_id <- patch_cells_subset$sample_id[1]
      patient_id <- patch_cells_subset$patient_id[1]
      RFS_status <- patch_cells_subset$RFS_status[1]
      Treatment <- patch_cells_subset$Treatment[1]
      
      # Spatial coordinates
      patch_coords <- patch_cells_subset[, c("Pos_X", "Pos_Y")]
      
      # Calculate spatial metrics
      # 1. Bounding box area
      x_range <- max(patch_coords$Pos_X) - min(patch_coords$Pos_X)
      y_range <- max(patch_coords$Pos_Y) - min(patch_coords$Pos_Y)
      bounding_box_area <- x_range * y_range
      
      # 2. Convex hull area (if >= 3 points)
      convex_hull_area <- NA
      if (nrow(patch_coords) >= 3) {
        tryCatch({
          hull_indices <- chull(patch_coords$Pos_X, patch_coords$Pos_Y)
          hull_coords <- patch_coords[hull_indices, ]
          # Calculate polygon area using shoelace formula
          n_hull <- nrow(hull_coords)
          convex_hull_area <- 0.5 * abs(sum(hull_coords$Pos_X * c(hull_coords$Pos_Y[-1], hull_coords$Pos_Y[1]) - 
                                              c(hull_coords$Pos_X[-1], hull_coords$Pos_X[1]) * hull_coords$Pos_Y))
        }, error = function(e) {
          convex_hull_area <<- NA
        })
      }
      
      # 3. Within-patch compactness metrics (cells per patch area)
      within_patch_compactness_bbox <- ifelse(bounding_box_area > 0, n_cells / bounding_box_area, NA)
      within_patch_compactness_convex <- ifelse(!is.na(convex_hull_area) & convex_hull_area > 0, n_cells / convex_hull_area, NA)
      
      # 4. Proper density metrics (relative to total tissue area)
      sample_tissue_area <- data[1,"width_px"] * data[1,"height_px"]
      patch_area_density <- ifelse(sample_tissue_area > 0, bounding_box_area / sample_tissue_area, NA)  # Patch coverage
      patch_convex_density <- ifelse(!is.na(convex_hull_area) & sample_tissue_area > 0, convex_hull_area / sample_tissue_area, NA)
      
      # 5. Patch shape metrics
      aspect_ratio <- ifelse(y_range > 0, x_range / y_range, NA)
      
      # 6. Patch compactness (convex hull area / bounding box area)
      shape_compactness <- ifelse(!is.na(convex_hull_area) & bounding_box_area > 0, convex_hull_area / bounding_box_area, NA)
      
      # 7. Average distance from centroid
      centroid_x <- mean(patch_coords$Pos_X)
      centroid_y <- mean(patch_coords$Pos_Y)
      distances_from_centroid <- sqrt((patch_coords$Pos_X - centroid_x)^2 + (patch_coords$Pos_Y - centroid_y)^2)
      avg_distance_from_centroid <- mean(distances_from_centroid)
      max_distance_from_centroid <- max(distances_from_centroid)
      
      # Store results
      patch_result <- data.frame(
        patch_id = paste0(patch_column_name,"_", patch_id),
        patch_type = patch_type,
        sample_id = sample_id,
        patient_id = patient_id,
        RFS_status = RFS_status,
        Treatment = Treatment,
        n_cells = n_cells,
        bounding_box_area = bounding_box_area,
        convex_hull_area = convex_hull_area,
        x_range = x_range,
        y_range = y_range,
        within_patch_compactness_bbox = within_patch_compactness_bbox,
        within_patch_compactness_convex = within_patch_compactness_convex,
        patch_area_density = patch_area_density,
        patch_convex_density = patch_convex_density,
        aspect_ratio = aspect_ratio,
        shape_compactness = shape_compactness,
        avg_distance_from_centroid = avg_distance_from_centroid,
        max_distance_from_centroid = max_distance_from_centroid,
        total_tissue_area = sample_tissue_area,
        stringsAsFactors = FALSE
      )
      
      patch_size_results <- rbind(patch_size_results, patch_result)
    
  }
  
  return(patch_size_results)
}

analysis_patch_characteristic_with_expression <- function(data,
                                                          expr_matrix_clean,analysis_markers,
                                                          patch_column_name,patche_ids,patch_type){
  # Initialize results dataframe
  patch_size_expression_data <- data.frame()
  
  for (patch_id in patche_ids) {
    patch_cells_subset <- data[data[`patch_column_name`] == patch_id, ]
    
      
      # Basic patch info
      sample_id <- patch_cells_subset$sample_id[1]
      patient_id <- patch_cells_subset$patient_id[1]
      RFS_status <- patch_cells_subset$RFS_status[1]
      Treatment <- patch_cells_subset$Treatment[1]
      n_cells <- nrow(patch_cells_subset)
      
      # Calculate spatial size metrics
      patch_coords <- patch_cells_subset[, c("Pos_X", "Pos_Y")]
      x_range <- max(patch_coords$Pos_X) - min(patch_coords$Pos_X)
      y_range <- max(patch_coords$Pos_Y) - min(patch_coords$Pos_Y)
      bounding_box_area <- x_range * y_range
      
      # Calculate convex hull area
      convex_hull_area <- NA
      if (nrow(patch_coords) >= 3) {
        tryCatch({
          hull_indices <- chull(patch_coords$Pos_X, patch_coords$Pos_Y)
          hull_coords <- patch_coords[hull_indices, ]
          convex_hull_area <- 0.5 * abs(sum(hull_coords$Pos_X * c(hull_coords$Pos_Y[-1], hull_coords$Pos_Y[1]) - 
                                              c(hull_coords$Pos_X[-1], hull_coords$Pos_X[1]) * hull_coords$Pos_Y))
        }, error = function(e) {
          convex_hull_area <<- NA
        })
      }
      
      # Get expression data for this patch
      patch_expr <- expr_matrix_clean[, match(rownames(patch_cells_subset), colnames(expr_matrix_clean))]
      
      # Calculate mean expression per marker for this patch
      patch_means <- apply(patch_expr, 1, mean, na.rm = TRUE)
      
      # Create patch data with size metrics and log-transformed expression
      patch_data <- data.frame(
        patch_id = paste0(patch_column_name,"_", patch_id),
        patch_type = patch_type,
        sample_id = sample_id,
        patient_id = patient_id,
        RFS_status = RFS_status,
        Treatment = Treatment,
        n_cells = n_cells,
        bounding_box_area = bounding_box_area,
        convex_hull_area = convex_hull_area,
        log_n_cells = log2(n_cells),
        log_bounding_box_area = log2(bounding_box_area + 0.001),
        log_convex_hull_area = log2(ifelse(is.na(convex_hull_area), 0.001, convex_hull_area) + 0.001),
        stringsAsFactors = FALSE
      )
      
      # Add log-transformed marker expression
      for (marker in analysis_markers) {
        patch_data[[paste0("mean_", marker)]] <- patch_means[marker]
      }
      
      patch_size_expression_data <- rbind(patch_size_expression_data, patch_data)
    
  }
  
  return(patch_size_expression_data)
}

# Function to process patches and calculate microenvironment composition
process_patches <- function(cell_data, patch_column, patch_type_label, patch_id_prefix) {
  
  # Get cells belonging to patches of this type
  patch_cells <- cell_data[!is.na(cell_data[[patch_column]]), ]
  patch_ids <- unique(patch_cells[[patch_column]])
  sub_celltypes <- unique(patch_cells$sub_celltype)
  
  print(paste("Processing", length(patch_ids), patch_type_label, "patches..."))
  
  patch_results <- data.frame()
  
  for (patch_id in patch_ids) {
    # Get all cells within this patch (including tumor and microenvironment)
    patch_all_cells <- patch_cells[patch_cells[[patch_column]] == patch_id, ]
      
      # Basic patch info
      sample_id <- patch_all_cells$sample_id[1]
      patient_id <- patch_all_cells$patient_id[1]
      RFS_status <- patch_all_cells$RFS_status[1]
      Treatment <- patch_all_cells$Treatment[1]
      
      # Skip if RFS_status is NA
      if (is.na(RFS_status)) next
      
      # Separate tumor cells from microenvironment cells
      # Tumor cells are typically epithelial cells (EpCAM+, cytokeratin+)
      # Microenvironment cells are immune, stromal, endothelial, etc.
      tumor_cell_types <- sub_celltypes[startsWith(sub_celltypes,prefix = "EC")]   # Adjust based on your data
      
      # Get microenvironment cells (non-tumor cells within the patch)
      microenv_cells <- patch_all_cells[!patch_all_cells$sub_celltype %in% tumor_cell_types, ]
      
      # Calculate total cells and microenvironment composition
      total_cells_in_patch <- nrow(patch_all_cells)
      tumor_cells_in_patch <- nrow(patch_all_cells[patch_all_cells$sub_celltype %in% tumor_cell_types, ])
      microenv_cells_in_patch <- nrow(microenv_cells)
      
      # Calculate microenvironment cell type proportions
      if (nrow(microenv_cells) > 0) {
        cell_type_counts <- table(microenv_cells$sub_celltype)
        cell_type_proportions <- cell_type_counts / total_cells_in_patch  # Proportion of total patch
        cell_type_microenv_proportions <- cell_type_counts / nrow(microenv_cells)  # Proportion of microenvironment
      } else {
        cell_type_proportions <- numeric(0)
        cell_type_microenv_proportions <- numeric(0)
      }
      
      # Calculate spatial metrics for the patch
      patch_coords <- patch_all_cells[, c("Pos_X", "Pos_Y")]
      x_range <- max(patch_coords$Pos_X) - min(patch_coords$Pos_X)
      y_range <- max(patch_coords$Pos_Y) - min(patch_coords$Pos_Y)
      patch_area <- x_range * y_range
      
      # Create base patch data
      patch_data <- data.frame(
        patch_id = paste0(patch_id_prefix, "_", patch_id),
        patch_type = patch_type_label,
        sample_id = sample_id,  
        patient_id = patient_id,
        RFS_status = RFS_status,
        RFS_group = ifelse(RFS_status == 0, "No_Recurrence", "Recurrence"),
        Treatment = Treatment,
        total_cells = total_cells_in_patch,
        tumor_cells = tumor_cells_in_patch,
        microenv_cells = microenv_cells_in_patch,
        microenv_fraction = microenv_cells_in_patch / total_cells_in_patch,
        patch_area = patch_area,
        microenv_density = microenv_cells_in_patch / patch_area,  # cells per area unit
        stringsAsFactors = FALSE
      )
      
      # Add cell type proportions (of total patch)
      all_cell_types <- unique(cell_data$sub_celltype)
      for (cell_type in all_cell_types) {
        if (cell_type %in% names(cell_type_proportions)) {
          patch_data[[paste0("prop_", cell_type)]] <- as.numeric(cell_type_proportions[cell_type])
          patch_data[[paste0("count_", cell_type)]] <- as.numeric(cell_type_counts[cell_type])
        } else {
          patch_data[[paste0("prop_", cell_type)]] <- 0
          patch_data[[paste0("count_", cell_type)]] <- 0
        }
      }
      
      # Add microenvironment-specific proportions (of microenvironment only)
      for (cell_type in all_cell_types) {
        if (cell_type %in% names(cell_type_microenv_proportions)) {
          patch_data[[paste0("microenv_prop_", cell_type)]] <- as.numeric(cell_type_microenv_proportions[cell_type])
        } else {
          patch_data[[paste0("microenv_prop_", cell_type)]] <- 0
        }
      }
      
      patch_results <- rbind(patch_results, patch_data)
    
  }
  
  return(patch_results)
}

# Function to extract microenvironment cells from patches with their expression data
extract_microenv_cells_with_expression <- function(cell_data, expr_matrix, patch_column, patch_type_label) {
  
  # Get cells belonging to patches of this type
  patch_cells <- cell_data[!is.na(cell_data[[patch_column]]), ]
  patch_ids <- unique(patch_cells[[patch_column]])
  
  print(paste("Extracting microenvironment cells from", length(patch_ids), patch_type_label, "patches..."))
  
  microenv_cells_data <- data.frame()
  
  for (patch_id in patch_ids) {
    # Get all cells within this patch
    patch_all_cells <- patch_cells[patch_cells[[patch_column]] == patch_id, ]
    
    # Require minimum 10 total cells per patch
    if (nrow(patch_all_cells) >= 10) {
      
      # Get microenvironment cells (non-tumor cells within the patch)
      microenv_cells <- patch_all_cells[!patch_all_cells$sub_celltype %in% tumor_cell_types, ]
      
      # Require at least some microenvironment cells
      if (nrow(microenv_cells) >= 3) {
        
        # Add patch information to microenvironment cells
        microenv_cells$patch_id <- paste0(patch_type_label, "_", patch_id)
        microenv_cells$patch_type <- patch_type_label
        microenv_cells$patch_numeric_id <- patch_id
        
        # Get expression data for these microenvironment cells
        cell_indices <- match(rownames(microenv_cells), colnames(expr_matrix))
        cell_expr <- expr_matrix[, cell_indices, drop = FALSE]
        
        # Add expression data to cell metadata
        for (marker in analysis_markers) {
          microenv_cells[[paste0("expr_", marker)]] <- as.numeric(cell_expr[marker, ])
        }
        
        microenv_cells_data <- rbind(microenv_cells_data, microenv_cells)
      }
    }
  }
  
  return(microenv_cells_data)
}

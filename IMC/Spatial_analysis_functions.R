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

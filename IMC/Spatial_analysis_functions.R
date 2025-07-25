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

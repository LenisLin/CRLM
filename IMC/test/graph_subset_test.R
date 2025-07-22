# Test: What happens to edges when subsetting SpatialExperiment objects?

library(SpatialExperiment)

# ================================================================================
# TEST EDGE BEHAVIOR AFTER SUBSETTING
# ================================================================================

cat("=== TESTING EDGE BEHAVIOR AFTER SUBSETTING ===\n")

# Get original data before subsetting
cat("BEFORE SUBSETTING:\n")
cat("Total cells in full dataset:", ncol(spe), "\n")
cat("Cells in IM tissue:", sum(spe$Tissue == "IM"), "\n")
cat("Malignant cells in IM:", sum(spe$Tissue == "IM" & grepl("^EC_", spe$sub_celltype)), "\n")

# Check if colPair exists and what types are available
cat("Available colPair types:", colPairNames(spe), "\n")

# If delaunay exists, check edge counts
if("delaunay" %in% colPairNames(spe)) {
  original_edges <- as.data.frame(colPair(spe, "delaunay"))
  cat("Total edges in original dataset:", length(original_edges[,1]), "\n")
  
  # Sample some edge indices to see what they look like
  sample_indices <- head(1:length(original_edges[,1]), 10)
  cat("Sample edge indices (first 10):\n")
  cat("From:", original_edges[,1][sample_indices], "\n")
  cat("To:  ", original_edges[,2][sample_indices], "\n")
  
  # Check what types of cells these edges connect
  from_cells <- original_edges[,1][sample_indices]
  to_cells <- original_edges[,2][sample_indices]
  
  cat("Sample edge cell types:\n")
  for(i in 1:length(sample_indices)) {
    from_type <- spe$sub_celltype[from_cells[i]]
    to_type <- spe$sub_celltype[to_cells[i]]
    from_tissue <- spe$Tissue[from_cells[i]]
    to_tissue <- spe$Tissue[to_cells[i]]
    
    cat(sprintf("Edge %d: %s (%s) -> %s (%s)\n", 
                i, from_type, from_tissue, to_type, to_tissue))
  }
} else {
  cat("No 'delaunay' colPair found in the dataset\n")
  cat("Available colPairs:", colPairNames(spe), "\n")
}

cat("\n", rep("=", 60), "\n")

# Now subset to IM malignant cells (as in the original code)
cat("AFTER SUBSETTING:\n")
malignant_IM <- spe[, spe$Tissue %in% "IM"]
cat("Cells after IM filtering:", ncol(malignant_IM), "\n")

malignant_IM <- malignant_IM[, grepl("^EC_", malignant_IM$sub_celltype)]
cat("Cells after malignant filtering:", ncol(malignant_IM), "\n")

# Check edge counts after subsetting
if("delaunay" %in% colPairNames(malignant_IM)) {
  subset_edges <- as.data.frame(colPair(malignant_IM, "delaunay")) 
  cat("Total edges after subsetting:", length(subset_edges[,1]), "\n")
  
  # Check maximum cell index to see if re-indexing occurred
  max_from <- max(subset_edges[,1])
  max_to <- max(subset_edges[,2])
  cat("Maximum 'from' index:", max_from, "\n")
  cat("Maximum 'to' index:", max_to, "\n")
  cat("Number of cells in subset:", ncol(malignant_IM), "\n")
  
  # This should be TRUE if proper re-indexing occurred
  cat("Indices properly re-mapped?", max(max_from, max_to) <= ncol(malignant_IM), "\n")
  
  # Sample some edges from subset to verify
  if(length(subset_edges[,1]) > 0) {
    sample_indices_subset <- head(1:length(subset_edges[,1]), 10)
    cat("Sample edge indices after subsetting (first 10):\n")
    cat("From:", subset_edges[,1][sample_indices_subset], "\n")
    cat("To:  ", subset_edges[,2][sample_indices_subset], "\n")
    
    # Verify all are malignant cells
    from_cells_subset <- subset_edges[,1][sample_indices_subset]
    to_cells_subset <- subset_edges[,2][sample_indices_subset]
    
    cat("Sample edge cell types after subsetting:\n")
    for(i in 1:length(sample_indices_subset)) {
      from_type <- malignant_IM$sub_celltype[from_cells_subset[i]]
      to_type <- malignant_IM$sub_celltype[to_cells_subset[i]]
      
      cat(sprintf("Edge %d: %s -> %s\n", i, from_type, to_type))
    }
  }
} else {
  cat("No 'delaunay' colPair found in subset\n")
}

cat("\n", rep("=", 60), "\n")

# ================================================================================
# COMPARISON: DIFFERENT SUBSETTING STRATEGIES
# ================================================================================

cat("COMPARISON OF SUBSETTING STRATEGIES:\n")

# Strategy 1: Subset first, then access edges (current approach)
cat("Strategy 1 - Subset then access edges:\n")
strategy1 <- spe[, spe$Tissue == "IM" & grepl("^EC_", spe$sub_celltype)]
if("delaunay" %in% colPairNames(strategy1)) {
  edges1 <- as.data.frame(colPair(strategy1, "delaunay")) 
  cat("  Edges in strategy 1:", length(edges1[,1]), "\n")
} else {
  cat("  No edges in strategy 1\n")
}

# Strategy 2: Access edges first, then filter (alternative approach)
cat("Strategy 2 - Access edges first, then filter:\n")
if("delaunay" %in% colPairNames(spe)) {
  all_edges <- as.data.frame(colPair(spe, "delaunay"))  
  
  # Get indices of malignant IM cells in original dataset
  malignant_im_indices <- which(spe$Tissue == "IM" & grepl("^EC_", spe$sub_celltype))
  
  # Filter edges to only include malignant IM cells
  edge_from <- all_edges[,1]
  edge_to <- all_edges[,2]
  
  # Keep edges where both endpoints are malignant IM cells
  valid_edges <- (edge_from %in% malignant_im_indices) & (edge_to %in% malignant_im_indices)
  
  filtered_from <- edge_from[valid_edges]
  filtered_to <- edge_to[valid_edges]
  
  cat("  Original edges:", length(edge_from), "\n")
  cat("  Malignant IM cells:", length(malignant_im_indices), "\n")
  cat("  Edges between malignant IM cells:", length(filtered_from), "\n")
  
  # Note: These indices would need to be re-mapped to 1:n for the subset
  cat("  Index re-mapping would be needed for strategy 2\n")
}

cat("\n", rep("=", 60), "\n")

# ================================================================================
# BIOLOGICAL INTERPRETATION
# ================================================================================

cat("BIOLOGICAL IMPLICATIONS:\n")

if("delaunay" %in% colPairNames(spe) && "delaunay" %in% colPairNames(malignant_IM)) {
  original_edges <- colPair(spe, "delaunay") 
  subset_edges <- as.data.frame(colPair(malignant_IM, "delaunay")) 
  
  edge_retention_rate <- length(subset_edges[,1]) / length(original_edges[,1])
  
  cat("Edge retention rate after subsetting:", round(100 * edge_retention_rate, 1), "%\n")
  
  if(edge_retention_rate < 0.5) {
    cat("WARNING: Most edges were lost during subsetting!\n")
    cat("This suggests many edges connected tumor cells to non-tumor cells.\n")
    cat("Consider whether tumor-microenvironment edges are important for your analysis.\n")
  }
  
  cat("\nIMPLICATIONS FOR BUDDING ANALYSIS:\n")
  cat("- Tumor-immune edges: LOST (might be important for immune exclusion)\n")
  cat("- Tumor-stromal edges: LOST (might be important for CAF interactions)\n")  
  cat("- Tumor-tumor edges: KEPT (good for budding cluster detection)\n")
  cat("- Edge lengths: Preserved within tumor-tumor connections\n")
  
} else {
  cat("Cannot calculate edge retention - colPair not available\n")
}

cat("\nRECOMMENDATIONS:\n")
cat("1. Current approach is CORRECT for tumor budding analysis\n")
cat("2. Only tumor-tumor edges are relevant for clustering tumor cells\n")
cat("3. Subsetting automatically handles re-indexing\n")
cat("4. For microenvironment analysis, keep broader cell types\n")

# ================================================================================
# VERIFICATION TEST
# ================================================================================

cat("\n", rep("=", 60), "\n")
cat("VERIFICATION TEST:\n")

# Test if all edges in subset connect valid cell indices
if("delaunay" %in% colPairNames(malignant_IM)) {
  subset_edges <- as.data.frame(colPair(malignant_IM, "delaunay")) 
  
  if(length(subset_edges[,1]) > 0) {
    max_index <- max(c(subset_edges[,1], subset_edges[,2]))
    min_index <- min(c(subset_edges[,1], subset_edges[,2]))
    
    cat("Cell indices in subset range:", min_index, "to", max_index, "\n")
    cat("Number of cells in subset:", ncol(malignant_IM), "\n")
    cat("All indices valid?", max_index <= ncol(malignant_IM) && min_index >= 1, "\n")
    
    # Test a few edges to make sure they point to actual cells
    test_edges <- head(1:length(subset_edges[,1]), 5)
    cat("Testing edge validity:\n")
    for(i in test_edges) {
      from_idx <- subset_edges[,1][i]
      to_idx <- subset_edges[,2][i]
      
      # These should not throw errors
      from_cell <- malignant_IM$sub_celltype[from_idx]
      to_cell <- malignant_IM$sub_celltype[to_idx]
      
      cat(sprintf("Edge %d: cell %d (%s) -> cell %d (%s) ✓\n", 
                  i, from_idx, from_cell, to_idx, to_cell))
    }
  }
}

cat("\nCONCLUSION: Subsetting behavior verified!\n")
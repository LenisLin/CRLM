## scRNA-seq for CRLM

library(Seurat)
library(DoubletFinder)
library(harmony)
library(Matrix)
library(glmGamPoi)
library(R.utils)

library(dplyr)
library(tidyr)

library(ggplot2)

library(RColorBrewer)
library(cowplot)
library(patchwork)
library(viridis)
library(dittoSeq)

## Function to handle duplicate feature names in expression matrices
## =============================================================================

#' Resolve Duplicate Feature Names in Expression Matrix
#' 
#' This function identifies duplicate feature names (genes) in an expression matrix
#' and averages their expression values. For count data, it applies ceiling operation
#' to maintain integer values.
#' 
#' @param expression_matrix A numeric matrix or sparse matrix (features x cells)
#' @param feature_names A character vector of feature names (same length as nrow(expression_matrix))
#' @param is_count_data Logical, whether the data represents raw counts (default: TRUE)
#' @param verbose Logical, whether to print processing information (default: TRUE)
#' 
#' @return A matrix with unique feature names and aggregated expression values
#' 
#' @examples
#' # Example with duplicate gene names
#' set.seed(123)
#' expr_matrix <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
#' gene_names <- c("GENE1", "GENE2", "GENE1", "GENE3", "GENE2", 
#'                 "GENE4", "GENE5", "GENE1", "GENE6", "GENE7")
#' 
#' cleaned_matrix <- resolve_duplicate_features(expr_matrix, gene_names)

resolve_duplicate_features <- function(expression_matrix, 
                                     feature_names, 
                                     is_count_data = TRUE, 
                                     verbose = TRUE) {
    
    # Input validation
    if (length(feature_names) != nrow(expression_matrix)) {
        stop("Length of feature_names must equal number of rows in expression_matrix")
    }
    
    # Check for duplicates
    duplicated_features <- duplicated(feature_names)
    unique_features <- unique(feature_names)
    
    if (verbose) {
        cat("Original number of features:", length(feature_names), "\n")
        cat("Number of unique features:", length(unique_features), "\n")
        cat("Number of duplicated features:", sum(duplicated_features), "\n")
    }
    
    # If no duplicates, return original matrix with proper rownames
    if (!any(duplicated_features)) {
        if (verbose) cat("No duplicate features found. Returning original matrix.\n")
        rownames(expression_matrix) <- feature_names
        return(expression_matrix)
    }
    
    # Handle sparse matrices
    is_sparse <- inherits(expression_matrix, "sparseMatrix")
    
    if (is_sparse && verbose) {
        cat("Processing sparse matrix...\n")
    }
    
    # Create new matrix for unique features
    if (is_sparse) {
        # For sparse matrices, convert to dense for processing
        dense_matrix <- as.matrix(expression_matrix)
        new_matrix <- matrix(0, nrow = length(unique_features), ncol = ncol(dense_matrix))
        colnames(new_matrix) <- colnames(dense_matrix)
    } else {
        new_matrix <- matrix(0, nrow = length(unique_features), ncol = ncol(expression_matrix))
        colnames(new_matrix) <- colnames(expression_matrix)
    }
    
    rownames(new_matrix) <- unique_features
    
    # Aggregate duplicate features
    if (verbose) cat("Aggregating duplicate features...\n")
    
    for (i in seq_along(unique_features)) {
        feature <- unique_features[i]
        
        # Find all rows with this feature name
        feature_indices <- which(feature_names == feature)
        
        if (length(feature_indices) == 1) {
            # No duplicates for this feature
            if (is_sparse) {
                new_matrix[i, ] <- dense_matrix[feature_indices, ]
            } else {
                new_matrix[i, ] <- expression_matrix[feature_indices, ]
            }
        } else {
            # Multiple instances - take mean
            if (is_sparse) {
                feature_data <- dense_matrix[feature_indices, , drop = FALSE]
            } else {
                feature_data <- expression_matrix[feature_indices, , drop = FALSE]
            }
            
            # Calculate mean across duplicate features
            mean_values <- colMeans(feature_data)
            
            # For count data, apply ceiling to get integers
            if (is_count_data) {
                mean_values <- ceiling(mean_values)
            }
            
            new_matrix[i, ] <- mean_values
            
            if (verbose) {
                cat("  Feature", feature, "had", length(feature_indices), "duplicates - averaged\n")
            }
        }
    }
    
    # Convert back to sparse if original was sparse
    if (is_sparse) {
        if (verbose) cat("Converting back to sparse matrix...\n")
        new_matrix <- Matrix::Matrix(new_matrix, sparse = TRUE)
    }
    
    # Ensure integer type for count data
    if (is_count_data && !is_sparse) {
        storage.mode(new_matrix) <- "integer"
    }
    
    if (verbose) {
        cat("Processing complete!\n")
        cat("Final matrix dimensions:", nrow(new_matrix), "x", ncol(new_matrix), "\n")
    }
    gc()
    return(new_matrix)
}


# Function to perform QC on individual Seurat object
# =============================================================================

#' Perform Quality Control on Seurat Object (V5 Compatible)
#' 
#' This function performs comprehensive quality control and preprocessing
#' on a Seurat object following Seurat V5 best practices
#' 
#' @param seurat_obj Seurat object
#' @param nFeature_RNA_min Minimum number of features per cell (default: 200)
#' @param nFeature_RNA_max Maximum number of features per cell (default: 5000)
#' @param percent.mt Maximum mitochondrial gene percentage (default: 20)
#' @param percent.ribo Maximum ribosomal gene percentage (default: NULL, no filtering)
#' @param nCount_RNA_min Minimum number of UMIs per cell (default: NULL, no filtering)
#' @param nCount_RNA_max Maximum number of UMIs per cell (default: NULL, no filtering)
#' @param n_hvgs Number of highly variable genes (default: 2000)
#' @param use_sct Logical, whether to use SCTransform instead of traditional workflow (default: FALSE)
#' @param pca_dims Number of PCA dimensions to use (default: 1:30)
#' @param umap_dims PCA dimensions for UMAP (default: 1:15)
#' @param cluster_resolution Clustering resolution (default: 0.5)
#' @param mt_pattern Pattern for mitochondrial genes (default: "^MT-")
#' @param ribo_pattern Pattern for ribosomal genes (default: "^RP[SL]")
#' @param verbose Logical, whether to print progress information (default: TRUE)
#' 
#' @return Processed Seurat object with QC metrics and preprocessing completed

perform_qc <- function(seurat_obj, 
                         nFeature_RNA_min = 200, 
                         nFeature_RNA_max = 5000, 
                         percent.mt = 20,
                         percent.ribo = NULL,
                         nCount_RNA_min = NULL,
                         nCount_RNA_max = NULL,
                         n_hvgs = 2000,
                         use_sct = FALSE,
                         pca_dims = 1:30,
                         umap_dims = 1:15,
                         cluster_resolution = 0.5,
                         mt_pattern = "^MT-",
                         ribo_pattern = "^RP[SL]",
                         verbose = TRUE) {
    
    # Load required packages
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat package is required")
    }
    
    if (verbose) {
        cat("=== Starting Quality Control Analysis (Seurat V5) ===\n")
        cat("Initial object dimensions:", nrow(seurat_obj), "features x", ncol(seurat_obj), "cells\n")
    }
    
    # Step 1: Ensure proper layer structure for V5
    if (verbose) cat("\n1. Checking and organizing data layers...\n")
    
    # Make sure we have the counts layer
    if (!"counts" %in% names(seurat_obj@assays$RNA@layers)) {
        if (verbose) cat("  Converting to V5 layer structure...\n")
        seurat_obj <- UpdateSeuratObject(seurat_obj)
    }
    
    # Join layers if they are split (common after integration)
    if (length(seurat_obj@assays$RNA@layers$counts) > 1) {
        if (verbose) cat("  Joining split layers...\n")
        seurat_obj <- JoinLayers(seurat_obj)
    }
    
    # Step 2: Calculate comprehensive QC metrics
    if (verbose) cat("\n2. Calculating QC metrics...\n")
    
    # Mitochondrial gene percentage
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
    
    # Ribosomal gene percentage (optional)
    if (!is.null(percent.ribo)) {
        seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = ribo_pattern)
    }
    
    # Cell complexity (log10(nFeature_RNA) / log10(nCount_RNA))
    seurat_obj[["complexity"]] <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
    
    if (verbose) {
        cat("  QC metrics calculated:\n")
        cat("    - Mitochondrial gene %\n")
        if (!is.null(percent.ribo)) cat("    - Ribosomal gene %\n")
        cat("    - Cell complexity\n")
    }
    
    # Step 3: Report QC statistics before filtering
    if (verbose) {
        cat("\n3. QC Statistics (before filtering):\n")
        cat("  nFeature_RNA - Min:", min(seurat_obj$nFeature_RNA), 
            "Max:", max(seurat_obj$nFeature_RNA), 
            "Median:", median(seurat_obj$nFeature_RNA), "\n")
        cat("  nCount_RNA - Min:", min(seurat_obj$nCount_RNA), 
            "Max:", max(seurat_obj$nCount_RNA), 
            "Median:", median(seurat_obj$nCount_RNA), "\n")
        cat("  percent.mt - Min:", round(min(seurat_obj$percent.mt), 2), 
            "Max:", round(max(seurat_obj$percent.mt), 2), 
            "Median:", round(median(seurat_obj$percent.mt), 2), "\n")
    }
    
    # Step 4: Filter low-quality cells
    if (verbose) cat("\n4. Filtering low-quality cells...\n")
    
    # Build filtering criteria
    filter_criteria <- paste0("nFeature_RNA > ", nFeature_RNA_min, 
                              " & nFeature_RNA < ", nFeature_RNA_max,
                              " & percent.mt < ", percent.mt)
    
    if (!is.null(nCount_RNA_min)) {
        filter_criteria <- paste0(filter_criteria, " & nCount_RNA > ", nCount_RNA_min)
    }
    if (!is.null(nCount_RNA_max)) {
        filter_criteria <- paste0(filter_criteria, " & nCount_RNA < ", nCount_RNA_max)
    }
    if (!is.null(percent.ribo)) {
        filter_criteria <- paste0(filter_criteria, " & percent.ribo < ", percent.ribo)
    }
    
    if (verbose) cat("  Filter criteria:", filter_criteria, "\n")
    
    cells_before <- ncol(seurat_obj)
    
    # Manual filtering to avoid Seurat subset bug
    # Create logical vector for cells that pass filtering
    metadata <- seurat_obj@meta.data
    
    keep_cells <- metadata$nFeature_RNA > nFeature_RNA_min & 
                  metadata$nFeature_RNA < nFeature_RNA_max & 
                  metadata$percent.mt < percent.mt
    
    # Apply additional filters if specified
    if (!is.null(nCount_RNA_min)) {
        keep_cells <- keep_cells & metadata$nCount_RNA > nCount_RNA_min
    }
    if (!is.null(nCount_RNA_max)) {
        keep_cells <- keep_cells & metadata$nCount_RNA < nCount_RNA_max
    }
    if (!is.null(percent.ribo)) {
        keep_cells <- keep_cells & metadata$percent.ribo < percent.ribo
    }
    
    # Remove NA values from logical vector
    keep_cells[is.na(keep_cells)] <- FALSE
    
    cells_after <- sum(keep_cells)
    
    if (verbose) {
        cat("  Cells before filtering:", cells_before, "\n")
        cat("  Cells after filtering:", cells_after, "\n")
        cat("  Cells removed:", cells_before - cells_after, 
            "(", round((cells_before - cells_after)/cells_before * 100, 2), "%)\n")
    }
    
    # Extract filtered data and create new object (workaround for Seurat subset bug)
    if (verbose) cat("  Creating new object with filtered data...\n")
    
    # Get the counts data from the appropriate layer
    if ("counts" %in% names(seurat_obj@assays$RNA@layers)) {
        counts_data <- seurat_obj@assays$RNA@layers$counts[, keep_cells]
    } else {
        # Fallback to counts slot if layers not available
        counts_data <- seurat_obj@assays$RNA@counts[, keep_cells]
    }
    
    # Filter metadata
    filtered_metadata <- metadata[keep_cells, , drop = FALSE]
    
    # Create new Seurat object with filtered data
    seurat_obj_filter <- CreateSeuratObject(
        counts = counts_data,
        meta.data = filtered_metadata,
        project = seurat_obj@project.name
    )
    rownames(seurat_obj_filter) <- rownames(seurat_obj)  # Ensure gene names are preserved
    colnames(seurat_obj_filter) <- rownames(filtered_metadata)  # Ensure cell barcodes are preserved

    if (verbose) cat("  New object created successfully!\n")
    
    # Step 5: Normalization and feature selection
    if (verbose) cat("\n5. Normalization and feature selection...\n")
    
    if (use_sct) {
        if (verbose) cat("  Using SCTransform workflow...\n")
        
        # SCTransform approach (recommended for V5)
        seurat_obj_filter <- SCTransform(seurat_obj_filter, 
                                 variable.features.n = n_hvgs,
                                 verbose = FALSE)
        
        # Set default assay to SCT
        DefaultAssay(seurat_obj_filter) <- "SCT"
        
    } else {
        if (verbose) cat("  Using traditional normalization workflow...\n")
        
        # Traditional approach
        seurat_obj_filter <- NormalizeData(seurat_obj_filter, verbose = FALSE)
        seurat_obj_filter <- FindVariableFeatures(seurat_obj_filter, 
                                          selection.method = "vst", 
                                          nfeatures = n_hvgs,
                                          verbose = FALSE)
        seurat_obj_filter <- ScaleData(seurat_obj_filter, verbose = FALSE)
    }
    
    # Step 6: Principal Component Analysis
    if (verbose) cat("\n6. Running Principal Component Analysis...\n")
    
    seurat_obj_filter <- RunPCA(seurat_obj_filter, 
                        npcs = max(pca_dims),
                        verbose = FALSE)
    
    # # Step 7: UMAP embedding
    # if (verbose) cat("\n7. Computing UMAP embedding...\n")
    
    # seurat_obj_filter <- RunUMAP(seurat_obj_filter, 
    #                      dims = umap_dims,
    #                      verbose = FALSE)
    
    # # Step 8: Graph-based clustering
    # if (verbose) cat("\n8. Performing graph-based clustering...\n")
    
    # seurat_obj_filter <- FindNeighbors(seurat_obj_filter, 
    #                            dims = umap_dims,
    #                            verbose = FALSE)
    # seurat_obj_filter <- FindClusters(seurat_obj_filter, 
    #                           resolution = cluster_resolution,
    #                           verbose = FALSE)
    
    # Step 9: Final summary
    if (verbose) {
        cat("\n=== QC Analysis Complete ===\n")
        cat("Final object dimensions:", nrow(seurat_obj_filter), "features x", ncol(seurat_obj_filter), "cells\n")
        cat("Number of clusters:", length(unique(Idents(seurat_obj_filter))), "\n")
        cat("Default assay:", DefaultAssay(seurat_obj_filter), "\n")
        
        if (use_sct) {
            cat("Workflow: SCTransform\n")
        } else {
            cat("Workflow: Traditional (NormalizeData + ScaleData)\n")
        }
    }
    
    return(seurat_obj_filter)
}

# =============================================================================
# Function to perform QC Doublet distinguish
perform_doublet_detect <- function(seurat_obj) {
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_kidney <- paramSweep(seurat_obj, PCs = 1:15, sct = FALSE)
    sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
    bcmvn_kidney <- find.pK(sweep.stats_kidney)
    pk_best <- bcmvn_kidney %>%
        dplyr::arrange(desc(BCmetric)) %>%
        dplyr::pull(pK) %>%
        .[1] %>%
        as.character() %>%
        as.numeric()

    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- seurat_obj@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(0.075 * nrow(seurat_obj@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    seurat_obj <- doubletFinder(seurat_obj, PCs = 1:15, pN = 0.25, pK = pk_best, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)

    return(seurat_obj)
}

# Function to align matrix to common gene set
align_matrix_to_genes <- function(matrix, target_genes) {
    current_genes <- rownames(matrix)
    
    # Create empty matrix with all genes
    aligned_matrix <- Matrix(0, 
                           nrow = length(target_genes), 
                           ncol = ncol(matrix),
                           sparse = TRUE)
    rownames(aligned_matrix) <- target_genes
    colnames(aligned_matrix) <- colnames(matrix)
    
    # Fill in values for genes present in current matrix
    common_genes <- intersect(current_genes, target_genes)
    aligned_matrix[common_genes, ] <- matrix[common_genes, ]
    
    return(aligned_matrix)
}

# Function to export Seurat object to files compatible with scanpy
export_seurat_to_scanpy <- function(seurat_obj, output_dir, prefix = "seurat_data", compress = TRUE) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract count matrix (use counts layer for raw data)
  if ("counts" %in% Layers(seurat_obj)) {
    count_matrix <- LayerData(seurat_obj, layer = "counts")
  } else {
    # Fallback to GetAssayData if no layers
    count_matrix <- GetAssayData(seurat_obj, slot = "counts")
  }
  
  # Ensure matrix is sparse
  if (!inherits(count_matrix, "sparseMatrix")) {
    count_matrix <- as(count_matrix, "sparseMatrix")
  }
  
  # File paths
  mtx_extension <- if (compress) ".mtx.gz" else ".mtx"
  features_extension <- if (compress) ".tsv.gz" else ".tsv"
  barcodes_extension <- if (compress) ".tsv.gz" else ".tsv"
  
  mtx_file <- file.path(output_dir, paste0(prefix, "_matrix", mtx_extension))
  features_file <- file.path(output_dir, paste0(prefix, "_features", features_extension))
  barcodes_file <- file.path(output_dir, paste0(prefix, "_barcodes", barcodes_extension))
  metadata_file <- file.path(output_dir, paste0(prefix, "_metadata.csv"))
  
  # 1. Save count matrix as .mtx file (with optional compression)
  if (compress) {
    # Write to temporary file first, then compress
    temp_mtx <- tempfile(fileext = ".mtx")
    Matrix::writeMM(count_matrix, temp_mtx)
    
    # Compress the file
    gzip(temp_mtx, destname = mtx_file, remove = TRUE)
  } else {
    Matrix::writeMM(count_matrix, mtx_file)
  }
  
  # 2. Save gene names (features)
  features_df <- data.frame(
    gene_id = rownames(count_matrix),
    gene_symbol = rownames(seurat_obj),
    gene_type = "Gene Expression"
  )
  
  if (compress) {
    write.table(features_df, gzfile(features_file), 
                sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
  } else {
    write.table(features_df, features_file, 
                sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  # 3. Save cell barcodes
  barcodes_df <- data.frame(barcode = colnames(count_matrix))
  
  if (compress) {
    write.table(barcodes_df, gzfile(barcodes_file), 
                sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
  } else {
    write.table(barcodes_df, barcodes_file, 
                sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  # 4. Save metadata as CSV (usually not compressed for readability)
  metadata <- seurat_obj@meta.data
  metadata$cell_barcode <- rownames(metadata)
  write.csv(metadata, metadata_file, row.names = FALSE)
  
  # Print summary
  compression_status <- if (compress) " (compressed)" else ""
  cat("Files saved to:", output_dir, "\n")
  cat("- Count matrix:", basename(mtx_file), compression_status, "\n")
  cat("- Features:", basename(features_file), compression_status, "\n") 
  cat("- Barcodes:", basename(barcodes_file), compression_status, "\n")
  cat("- Metadata:", basename(metadata_file), "\n")
  cat("Matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "cells\n")
  
  # Return file paths for convenience
  return(NULL)
}

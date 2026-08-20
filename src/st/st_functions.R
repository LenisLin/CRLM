# Shared spatial-analysis functions for Figure 7.
#
# Purpose: construct the Bin100 SpatialExperiment input, prepare the scRNA-seq
# reference for RCTD, map registered SpMap predictions to spots, attach and
# aggregate normalized RCTD weights, and run the prespecified GSEA stage.
#
# Callers: Figure7 scripts set their editable CONFIG st_functions path and
# source this module.
#
# Inputs: a Bin100 MatrixMarket quartet; annotated scRNA-seq raw counts;
# registered SpMap tile coordinates; normalized RCTD weights; cell-type
# aggregation mappings; ranked genes; and preloaded pathway collections.
#
# Outputs: in-memory SpatialExperiment objects, balanced scRNA-seq references,
# spot-level SpMap class assignments, RCTD composition tables, and GSEA result
# tables. Coordinates use the shared registered system. SpMap classes are
# 0 = TC, 1 = PIR-Niche, 2 = PSM-Niche, and 3 = other PT. The spatial analysis
# unit is one 50 x 50 um^2 Bin100 spot.
#
# Ordered use:
# 1. `read_anndata_to_spe()` loads the Bin100 input for each Figure 7 stage.
# 2. `balance_cell_types()` prepares the scRNA-seq reference for RCTD.
# 3. `.canonicalize_spot_ids()` aligns spot identifiers across inputs and RCTD
#    outputs.
# 4. `assign_predictions_to_spots()` maps registered SpMap tile classes to
#    Bin100 spots.
# 5. `combine_deconvolution()` adds normalized RCTD weights to retained spots.
# 6. `aggregate_celltype_proportions()` forms displayed composition groups.
# 7. `run_gsea_analysis()` evaluates ranked PIR-versus-PSM genes.

#' Canonicalize Bin100 spot identifiers
#'
#' Converts spot identifiers to trimmed character strings and expands numeric
#' scientific notation while retaining identifiers as character values. This allows
#' MatrixMarket metadata, Bin100 spot names, and RCTD outputs to
#' use the same representation.
#'
#' @param spot_ids A vector of spot identifiers coercible to character.
#'
#' @return A character vector of trimmed, canonicalized spot identifiers.
#'
#' @keywords internal
.canonicalize_spot_ids <- function(spot_ids) {
    # Normalize the shared spot key without converting identifiers to numeric values.
    spot_ids <- trimws(as.character(spot_ids))
    scientific_ids <- grepl(
        "^[0-9]+(?:\\.[0-9]+)?[eE][+-]?[0-9]+$",
        spot_ids,
        perl = TRUE
    )
    # Expand scientific notation explicitly so long IDs remain exact character keys.
    if (any(scientific_ids)) {
        spot_ids[scientific_ids] <- vapply(
            spot_ids[scientific_ids],
            function(spot_id) {
                parts <- regmatches(
                    spot_id,
                    regexec(
                        "^([0-9]+)(?:\\.([0-9]+))?[eE]([+-]?[0-9]+)$",
                        spot_id,
                        perl = TRUE
                    )
                )[[1]]
                digits <- paste0(parts[[2]], parts[[3]])
                decimal_position <- nchar(parts[[2]]) + as.integer(parts[[4]])
                if (decimal_position >= nchar(digits)) {
                    return(paste0(
                        digits,
                        strrep("0", decimal_position - nchar(digits))
                    ))
                }
                if (decimal_position > 0) {
                    return(paste0(
                        substr(digits, 1, decimal_position),
                        ".",
                        substr(digits, decimal_position + 1, nchar(digits))
                    ))
                }
                paste0("0.", strrep("0", -decimal_position), digits)
            },
            character(1)
        )
    }
    spot_ids
}

#' Read a processed Bin100 export as a SpatialExperiment
#'
#' Reads the FDZS-4 Bin100 MatrixMarket quartet, retains genes with nonempty
#' names, and constructs a spot-level SpatialExperiment with count data,
#' metadata, and registered coordinates.
#'
#' @param data_path Path to the directory containing `expression_profile.mtx`,
#'   `metadata.csv`, `row_names.csv`, and `column_names.csv`.
#' @param verbose Logical; whether to report input loading and removal of genes
#'   with empty names.
#'
#' @return A `SpatialExperiment` with a `counts` assay, `x` and `y` spatial
#'   coordinates, and metadata containing `sample_id` copied from `core_name`.
read_anndata_to_spe <- function(data_path, verbose = TRUE) {
    # Resolve the four MatrixMarket components as one spatial input bundle.
    required_files <- c(
        "expression_profile.mtx",
        "metadata.csv",
        "row_names.csv",
        "column_names.csv"
    )
    input_paths <- file.path(data_path, required_files)
    missing_files <- required_files[!file.exists(input_paths)]
    if (length(missing_files) > 0) {
        stop(
            "Missing required Bin100 input file(s): ",
            paste(missing_files, collapse = ", "),
            call. = FALSE
        )
    }

    if (verbose) message("Reading processed Bin100 MatrixMarket export")
    # The exported matrix is spot-by-gene; SpatialExperiment uses genes by spots.
    expression_matrix <- Matrix::t(Matrix::readMM(input_paths[[1]]))
    metadata_header <- names(utils::read.csv(
        input_paths[[2]],
        nrows = 0,
        check.names = FALSE
    ))
    metadata_col_classes <- rep(NA_character_, length(metadata_header))
    metadata_col_classes[[1]] <- "character"
    metadata <- utils::read.csv(
        input_paths[[2]],
        row.names = 1,
        colClasses = metadata_col_classes,
        check.names = FALSE
    )
    gene_names <- utils::read.csv(
        input_paths[[3]],
        header = FALSE,
        colClasses = "character",
        check.names = FALSE
    )[[1]]
    spot_names <- utils::read.csv(
        input_paths[[4]],
        header = FALSE,
        colClasses = "character",
        check.names = FALSE
    )[[1]]
    # Canonicalize both identifier sources before checking their positional alignment.
    rownames(metadata) <- .canonicalize_spot_ids(rownames(metadata))
    spot_names <- .canonicalize_spot_ids(spot_names)

    if (nrow(expression_matrix) != length(gene_names) ||
        ncol(expression_matrix) != length(spot_names)) {
        stop(
            "Matrix dimensions do not match row_names.csv and column_names.csv",
            call. = FALSE
        )
    }
    if (nrow(metadata) != length(spot_names) ||
        !identical(rownames(metadata), spot_names)) {
        stop(
            "metadata.csv rows are not aligned with column_names.csv",
            call. = FALSE
        )
    }

    required_metadata <- c("core_name", "x", "y")
    missing_metadata <- setdiff(required_metadata, colnames(metadata))
    if (length(missing_metadata) > 0) {
        stop(
            "metadata.csv is missing required column(s): ",
            paste(missing_metadata, collapse = ", "),
            call. = FALSE
        )
    }

    # Filter the feature axis and matrix rows together before assigning dimnames.
    valid_genes <- !is.na(gene_names) & gene_names != ""
    if (verbose && any(!valid_genes)) {
        message("Removing ", sum(!valid_genes), " genes with empty names")
    }
    expression_matrix <- expression_matrix[valid_genes, , drop = FALSE]
    gene_names <- gene_names[valid_genes]
    rownames(expression_matrix) <- gene_names
    colnames(expression_matrix) <- spot_names

    # Attach registered coordinates and metadata in the verified spot order.
    metadata$sample_id <- metadata$core_name
    spatial_coords <- as.matrix(metadata[, c("x", "y"), drop = FALSE])
    colnames(spatial_coords) <- c("x", "y")

    SpatialExperiment::SpatialExperiment(
        assays = list(counts = expression_matrix),
        spatialCoords = spatial_coords,
        colData = S4Vectors::DataFrame(metadata)
    )
}

#' Balance a single-cell reference across eligible cell types
#'
#' Retains cell types with at least `min_cells` cells and samples `k` cells per
#' retained type, with replacement when a retained type has fewer than `k`
#' cells. The result is used to construct the Figure 7 RCTD reference.
#'
#' @param sc_data A single-cell `SummarizedExperiment`-like object with cell
#'   identifiers in its column names and annotations in `colData()`.
#' @param cell_type_col Character scalar naming the cell-type column in
#'   `colData(sc_data)`.
#' @param k Number of cells sampled from each eligible cell type.
#' @param min_cells Minimum number of cells required for a cell type to be
#'   retained.
#' @param seed Integer random seed used for the per-cell-type sampling.
#'
#' @return A named list with `sc_data`, the balanced single-cell object, and
#'   `stats`, a list describing the original and balanced cell-type counts.
balance_cell_types <- function(sc_data, cell_type_col, k = 100,
                               min_cells = 10, seed = 619) {
    cell_types <- SummarizedExperiment::colData(sc_data)[[cell_type_col]]
    if (is.null(cell_types)) {
        stop("cell_type_col is not present in colData(sc_data)", call. = FALSE)
    }
    names(cell_types) <- colnames(sc_data)

    # Exclude sparse reference types before constructing equal-sized groups.
    cell_type_counts <- table(cell_types)
    valid_cell_types <- names(cell_type_counts)[cell_type_counts >= min_cells]
    if (length(valid_cell_types) == 0) {
        stop("No cell types meet min_cells", call. = FALSE)
    }

    set.seed(seed)
    sc_data_filtered <- sc_data[, cell_types %in% valid_cell_types]
    filtered_types <- SummarizedExperiment::colData(sc_data_filtered)[[cell_type_col]]
    names(filtered_types) <- colnames(sc_data_filtered)

    # Draw exactly k cells per type, using replacement only for smaller groups.
    sampled_cells <- character(0)
    for (cell_type in valid_cell_types) {
        candidates <- names(filtered_types)[filtered_types == cell_type]
        sampled_cells <- c(
            sampled_cells,
            sample(candidates, size = k, replace = length(candidates) < k)
        )
    }

    # Preserve repeated draws as distinct reference columns for RCTD attachment.
    sc_data_balanced <- sc_data_filtered[, sampled_cells]
    colnames(sc_data_balanced) <- make.unique(sampled_cells, sep = "_")

    balance_stats <- list(
        original_n_celltypes = length(cell_type_counts),
        filtered_n_celltypes = length(valid_cell_types),
        original_total_cells = ncol(sc_data),
        balanced_total_cells = ncol(sc_data_balanced),
        k = k,
        min_cells = min_cells,
        unique_barcodes = length(unique(colnames(sc_data_balanced)))
    )

    list(sc_data = sc_data_balanced, stats = balance_stats)
}

#' Assign registered SpMap tile predictions to Bin100 spots
#'
#' Assigns each spot to the first registered tile containing its coordinates.
#' For spots outside all tile bounds, assigns the class of the tile whose
#' top-left coordinate is nearest within the same sample.
#'
#' @param spatial_df A data frame with `x`, `y`, and `sample_id` columns for
#'   Bin100 spots in the shared registered coordinate system.
#' @param tile_anno A data frame with `sample_id`, tile bounds (`tile_left`,
#'   `tile_right`, `tile_top`, and `tile_bottom`), and `pred_class` columns.
#'
#' @return An integer vector of predicted SpMap classes in the row order of
#'   `spatial_df`.
assign_predictions_to_spots <- function(spatial_df, tile_anno) {
    # Confirm both coordinate tables expose the shared sample and geometry fields.
    required_spot_columns <- c("x", "y", "sample_id")
    required_tile_columns <- c(
        "sample_id", "tile_left", "tile_right", "tile_top", "tile_bottom",
        "pred_class"
    )
    missing_spot_columns <- setdiff(required_spot_columns, colnames(spatial_df))
    missing_tile_columns <- setdiff(required_tile_columns, colnames(tile_anno))
    if (length(missing_spot_columns) > 0) {
        stop(
            "spatial_df is missing required column(s): ",
            paste(missing_spot_columns, collapse = ", "),
            call. = FALSE
        )
    }
    if (length(missing_tile_columns) > 0) {
        stop(
            "tile_anno is missing required column(s): ",
            paste(missing_tile_columns, collapse = ", "),
            call. = FALSE
        )
    }

    unmatched <- !spatial_df$sample_id %in% tile_anno$sample_id
    if (any(unmatched)) {
        stop(
            "No tile annotations for sample(s): ",
            paste(unique(spatial_df$sample_id[unmatched]), collapse = ", "),
            call. = FALSE
        )
    }

    # Match within each sample, preferring tile containment in annotation order.
    pred_classes <- integer(nrow(spatial_df))
    for (i in seq_len(nrow(spatial_df))) {
        spot_x <- spatial_df$x[[i]]
        spot_y <- spatial_df$y[[i]]
        sample_tiles <- tile_anno[
            tile_anno$sample_id == spatial_df$sample_id[[i]],
            ,
            drop = FALSE
        ]

        in_tile <- which(
            spot_x >= sample_tiles$tile_left &
                spot_x <= sample_tiles$tile_right &
                spot_y >= sample_tiles$tile_top &
                spot_y <= sample_tiles$tile_bottom
        )
        if (length(in_tile) > 0) {
            selected_class <- sample_tiles$pred_class[[in_tile[[1]]]]
        } else {
            # For uncovered spots, inherit the nearest tile by top-left coordinate.
            distances <- sqrt(
                (sample_tiles$tile_left - spot_x)^2 +
                    (sample_tiles$tile_top - spot_y)^2
            )
            selected_class <- sample_tiles$pred_class[[which.min(distances)]]
        }
        pred_classes[[i]] <- as.integer(as.character(selected_class))
    }

    pred_classes
}

#' Attach normalized deconvolution weights to retained spatial spots
#'
#' Intersects canonical spot identifiers from a SpatialExperiment and a
#' deconvolution result table, retains the shared spots in spatial-object order,
#' and appends one `colData()` column per deconvolution cell type.
#'
#' @param spe A SpatialExperiment-like object with spot identifiers in
#'   `colnames(spe)`.
#' @param result_df A data frame or matrix-like object with spot identifiers in
#'   row names and normalized deconvolution weights in columns.
#' @param deconv_method Character prefix used to name appended weight columns.
#'
#' @return The subsetted spatial object with deconvolution columns named
#'   `<deconv_method>_<cell_type>` in `colData()`.
combine_deconvolution <- function(spe, result_df, deconv_method = "RCTD") {
    # Canonicalize both spot axes before intersecting them in spatial-object order.
    result_spots <- .canonicalize_spot_ids(rownames(result_df))
    if (anyDuplicated(result_spots)) {
        stop("result_df contains duplicate spot IDs", call. = FALSE)
    }

    spe_spots <- .canonicalize_spot_ids(colnames(spe))
    common_spots <- intersect(spe_spots, result_spots)
    spe_filtered <- spe[, match(common_spots, spe_spots)]
    colnames(spe_filtered) <- common_spots
    result_matrix <- as.matrix(
        result_df[match(common_spots, result_spots), , drop = FALSE]
    )
    rownames(result_matrix) <- common_spots

    # Append each aligned weight vector as a method-prefixed spot annotation.
    col_data <- SummarizedExperiment::colData(spe_filtered)
    for (cell_type in colnames(result_matrix)) {
        col_data[[paste0(deconv_method, "_", cell_type)]] <-
            result_matrix[, cell_type]
    }
    SummarizedExperiment::colData(spe_filtered) <- col_data

    spe_filtered
}

#' Aggregate deconvolution proportions into displayed cell-type groups
#'
#' Sums the specified deconvolution columns for each spot while preserving the
#' spot-level `sample_id` column.
#'
#' @param col_data A data frame or DataFrame containing `sample_id` and the
#'   deconvolution columns named in `mapping`.
#' @param mapping A named list whose names are output cell-type groups and whose
#'   values are one or more source column names in `col_data`.
#'
#' @return A data frame with `sample_id` followed by one summed proportion
#'   column for every named group in `mapping`.
aggregate_celltype_proportions <- function(col_data, mapping) {
    if (!"sample_id" %in% colnames(col_data)) {
        stop("col_data is missing sample_id", call. = FALSE)
    }

    # Sum mapped RCTD components per spot while retaining the sample grouping key.
    results <- data.frame(sample_id = col_data$sample_id, check.names = FALSE)
    for (cell_type in names(mapping)) {
        results[[cell_type]] <- rowSums(
            as.matrix(col_data[, mapping[[cell_type]], drop = FALSE])
        )
    }

    results
}

#' Run prespecified GSEA for one pathway collection
#'
#' Evaluates ranked genes against a preloaded pathway collection using the
#' Figure 7 fgsea settings, then records the source database and a lowercase,
#' space-separated pathway label for downstream selection and plotting.
#'
#' @param gene_ranks A named numeric vector of ranked gene statistics.
#' @param pathways A named list of gene vectors defining one pathway
#'   collection.
#' @param database Character label identifying the pathway collection.
#'
#' @return An `fgsea` result table with added `database` and `pathway_clean`
#'   columns.
run_gsea_analysis <- function(gene_ranks, pathways, database) {
    # Keep the supplied rank statistic and pathway collection in fgsea's input form.
    gsea_results <- fgsea::fgsea(
        pathways = pathways,
        stats = gene_ranks,
        minSize = 3,
        maxSize = 500,
        eps = 0
    )
    # Add collection provenance and plot-ready pathway labels to the result table.
    gsea_results$database <- database
    gsea_results$pathway_clean <- tolower(
        gsub("_", " ", gsea_results$pathway, fixed = TRUE)
    )

    gsea_results
}

# Target: Supplementary Figure 6A-C ROI interaction P-value matrices.
# Purpose: Generate the saved ClosePvalue.csv inputs consumed by the fixed-roster
# Supplementary Figure 6 plotting workflow.
# Inputs: A configured IMC SingleCellExperiment RDS with `Position`, current
# `SubType`, and ROI metadata in `colData()`; shared coordinate and spatial-
# interaction modules.
# Outputs: One symmetric ClosePvalue.csv matrix per ROI under
# permutation_<Tissue>/<ROI>/ClosePvalue. Significant enrichment uses
# `p_enrichment`, significant avoidance uses `1 - p_avoidance`, and other
# entries are recorded as 0.5.
# Ordered workflow: configure paths and seed; source the shared modules; load and
# filter the SCE; parse micrometre coordinates; run the shared 22-um, 1000-
# permutation interaction calculation; then write fixed-axis symmetric matrices.
# The configured analysis loads the canonical SCE.

# Edit these paths and seed before execution, then set RUN to TRUE.
RUN <- FALSE
figure_s6_precursor_config <- list(
    sce_path = "PATH/TO/FDZS1_IMC_processed.rds",
    permutation_root = "PATH/TO/saved_permutation_results",
    seed = 20260820L
)

#' Purpose: Resolve the absolute path to this Supplementary Figure 6 script.
#'
#' @return Character scalar containing the normalized script path.
.figure_s6_script_path <- function() {
    # Prefer the source-file context when this script is loaded from another R session.
    source_path <- tryCatch(sys.frame(1L)$ofile, error = function(error) NULL)
    if (is.null(source_path)) {
        # Fall back to the command-line file argument for direct Rscript execution.
        file_argument <- grep(
            "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
        )
        if (length(file_argument) != 1L) {
            stop("Cannot determine the Supplementary Figure 6 script path")
        }
        source_path <- sub("^--file=", "", file_argument)
    }
    normalizePath(source_path, mustWork = TRUE)
}

FIGURE_S6_SCRIPT_DIR <- dirname(.figure_s6_script_path())

DISTANCE_UM <- 22
N_PERMUTATIONS <- 1000L
SIGNIFICANCE_THRESHOLD <- 0.05
TISSUES <- c("IM", "PT", "TC")
CELLTYPE_ORDER <- c(
    "Macro_CD169", "Macro_HLADR", "Mono_CD11c", "SC_COLLAGEN", "CD8T", "CD4T", "B",
    "SC_FAP", "Macro_CD163", "UNKNOWN", "SC_Vimentin", "Mono_Intermediate", "Mono_Classic",
    "TC_EpCAM", "SC_aSMA", "TC_Ki67", "NK", "Treg", "Macro_CD11b", "TC_CAIX", "TC_VEGF"
)

#' Purpose: Source and verify the shared coordinate and spatial-interaction interfaces.
#'
#' @return Invisibly `NULL`; stops when a shared file or function is unavailable.
load_shared_interaction_modules <- function() {
    # Resolve the shared IMC interfaces relative to the public Figure S6 script.
    shared_dir <- normalizePath(
        file.path(FIGURE_S6_SCRIPT_DIR, "..", "..", "..", "..", "src", "imc"),
        mustWork = FALSE
    )
    shared_files <- file.path(shared_dir, c("coordinates.R", "spatial_interactions.R"))
    missing <- shared_files[!file.exists(shared_files)]
    if (length(missing) > 0L) {
        stop("Missing required src/imc file(s): ", paste(missing, collapse = ", "))
    }
    # Source both interfaces into the script environment used by the analysis.
    script_environment <- environment(load_shared_interaction_modules)
    invisible(lapply(shared_files, sys.source, envir = script_environment))

    # Require the coordinate parser and interaction engine exposed by the shared modules.
    required_functions <- c("position_to_xy", "calculate_spatial_interactions")
    unavailable <- required_functions[!vapply(
        required_functions, exists, logical(1L), mode = "function",
        envir = script_environment, inherits = FALSE
    )]
    if (length(unavailable) > 0L) {
        stop("Required src/imc interface is missing: ", paste(unavailable, collapse = ", "))
    }
    invisible(NULL)
}

#' Purpose: Load the configured SCE and prepare current Figure S6 cell metadata.
#'
#' @param sce_path Path to the configured IMC SingleCellExperiment RDS.
#' @param tissues Character vector of tissue labels retained for calculation.
#' @param subtype_order Complete subtype axis required by the Figure S6 consumer.
#' @return Cell metadata with numeric micrometre `x` and `y` columns.
prepare_interaction_metadata <- function(sce_path, tissues, subtype_order) {
    # Load cell metadata and require the identifiers, annotations, and Positions used by the engine.
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("SingleCellExperiment is required")
    }
    sce <- readRDS(sce_path)
    cell_metadata <- as.data.frame(SingleCellExperiment::colData(sce))
    required <- c("CellID", "PID", "ID", "Tissue", "SubType", "Position")
    missing <- setdiff(required, names(cell_metadata))
    if (length(missing) > 0L) {
        stop("SCE colData is missing: ", paste(missing, collapse = ", "))
    }

    # Restrict the cell-level input to the configured tissues and required metadata columns.
    cell_metadata <- cell_metadata[
        as.character(cell_metadata$Tissue) %in% tissues,
        required,
        drop = FALSE
    ]
    if (any(!stats::complete.cases(cell_metadata))) {
        stop("Required Figure S6 colData fields must be complete")
    }
    # Keep a complete fixed subtype axis and a one-to-one ROI-to-patient/tissue mapping.
    observed_subtypes <- sort(unique(as.character(cell_metadata$SubType)))
    if (!setequal(observed_subtypes, subtype_order)) {
        stop("Filtered SCE SubType values must match the Figure S6 cell-type roster")
    }
    roi_mapping <- unique(cell_metadata[, c("PID", "ID", "Tissue"), drop = FALSE])
    if (anyDuplicated(roi_mapping$ID)) {
        stop("Each ROI ID must map to one patient and tissue")
    }

    # Parse stored Positions into numeric micrometre coordinates for spatial calculations.
    coordinates <- position_to_xy(as.character(cell_metadata$Position))
    cell_metadata$x <- coordinates$x
    cell_metadata$y <- coordinates$y
    cell_metadata
}

#' Purpose: Convert one ROI's two permutation tails to a fixed-axis signed-tail matrix.
#'
#' @param roi_results Long-form results for one ROI returned by
#'   `calculate_spatial_interactions()`.
#' @param subtype_order Character vector defining both matrix axes.
#' @return Symmetric numeric ClosePvalue matrix in `subtype_order`; significant
#'   enrichment values are at or below 0.05, significant avoidance values are
#'   above 0.95, and other entries are 0.5.
close_pvalue_matrix <- function(roi_results, subtype_order) {
    # Map long-form subtype pairs to the fixed row and column axes.
    row_index <- match(roi_results$SubType_1, subtype_order)
    column_index <- match(roi_results$SubType_2, subtype_order)
    if (anyNA(row_index) || anyNA(column_index)) {
        stop("Interaction results contain a subtype outside the Figure S6 roster")
    }

    # Encode significant enrichment, significant avoidance, and neutral pairs on one tail scale.
    matrix <- matrix(
        NA_real_, nrow = length(subtype_order), ncol = length(subtype_order),
        dimnames = list(subtype_order, subtype_order)
    )
    close_pvalues <- ifelse(
        roi_results$p_enrichment <= SIGNIFICANCE_THRESHOLD,
        roi_results$p_enrichment,
        ifelse(
            roi_results$p_avoidance <= SIGNIFICANCE_THRESHOLD,
            1 - roi_results$p_avoidance,
            0.5
        )
    )
    # Mirror each unordered interaction pair to produce a symmetric consumer matrix.
    matrix[cbind(row_index, column_index)] <- close_pvalues
    matrix[cbind(column_index, row_index)] <- close_pvalues
    if (anyNA(matrix)) {
        stop("Interaction results do not contain every Figure S6 subtype pair")
    }
    matrix
}

#' Purpose: Write ROI signed-tail matrices in the Figure S6 saved-result layout.
#'
#' @param results Long-form output from `calculate_spatial_interactions()`.
#' @param permutation_root Root directory receiving `permutation_<Tissue>` directories.
#' @param tissues Character vector defining the tissues to write.
#' @param subtype_order Character vector defining matrix row and column order.
#' @return Invisibly returns the written ClosePvalue.csv paths.
write_close_pvalue_matrices <- function(
    results, permutation_root, tissues, subtype_order
) {
    # Write one fixed-axis matrix for every ROI, nested under its tissue result directory.
    written <- character()
    for (tissue in tissues) {
        # Preserve tissue as the outer analysis partition in the saved-result layout.
        tissue_results <- results[results$Tissue == tissue, , drop = FALSE]
        for (roi in unique(as.character(tissue_results$ID))) {
            # Convert and export the complete subtype-pair result for one ROI.
            roi_results <- tissue_results[tissue_results$ID == roi, , drop = FALSE]
            close_pvalues <- close_pvalue_matrix(roi_results, subtype_order)
            output_dir <- file.path(
                permutation_root, paste0("permutation_", tissue), roi, "ClosePvalue"
            )
            dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
            output_path <- file.path(output_dir, "ClosePvalue.csv")
            utils::write.csv(close_pvalues, output_path, quote = FALSE)
            written <- c(written, output_path)
        }
    }
    invisible(written)
}

#' Purpose: Run the configured Supplementary Figure 6 interaction precursor.
#'
#' @param config List with `sce_path`, `permutation_root`, and explicit integer `seed`
#'   entries.
#' @return Invisibly returns the written ClosePvalue.csv paths.
run_figure_s6_interaction_precursor <- function(config) {
    # Validate the input, output, and random-seed configuration before loading modules.
    required_config <- c("sce_path", "permutation_root", "seed")
    missing <- setdiff(required_config, names(config))
    if (length(missing) > 0L) {
        stop("Supplementary Figure 6 configuration is missing: ", paste(missing, collapse = ", "))
    }
    # Load shared interfaces and prepare the filtered cell-level spatial metadata.
    load_shared_interaction_modules()
    cell_metadata <- prepare_interaction_metadata(
        config$sce_path, TISSUES, CELLTYPE_ORDER
    )
    # Calculate ROI-level subtype-pair interactions with the configured distance and permutations.
    results <- calculate_spatial_interactions(
        cell_metadata,
        distance_um = DISTANCE_UM,
        n_permutations = N_PERMUTATIONS,
        seed = config$seed
    )
    # Export the long-form results as one symmetric ClosePvalue matrix per ROI.
    write_close_pvalue_matrices(
        results, config$permutation_root, TISSUES, CELLTYPE_ORDER
    )
}

if (sys.nframe() == 0L && isTRUE(RUN)) {
    run_figure_s6_interaction_precursor(figure_s6_precursor_config)
}

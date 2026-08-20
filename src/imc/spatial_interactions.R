# Module: IMC spatial-interaction helpers.
# Purpose: Calculate ROI-local subtype-pair permutation results or read and summarize saved
# interaction matrices.
# Callers: `figures/Figure3/supplementary/FigureS6/00_generate_interaction_pvalues.R`
# calls `calculate_spatial_interactions()` to generate the enrichment matrices consumed by
# the fixed-roster Supplementary Figure 6 workflow.
# Inputs: Cell metadata with numeric micrometre coordinates, or paths and layout metadata for
# saved interaction matrices.
# Outputs: Long-form ROI interaction results or tissue-level interaction-frequency summaries.
# Ordered use: Calculate permutation results from cell metadata and write the required saved
# matrices, or read existing matrices and then summarize them.

#' Purpose: Enumerate unordered subtype pairs, including same-subtype pairs.
#'
#' @param subtype_levels Character vector defining subtype order.
#' @return Data frame with `SubType_1` and `SubType_2` for every upper-triangular pair.
.interaction_subtype_pairs <- function(subtype_levels) {
    # Enumerate the upper triangle so every unordered pair appears exactly once.
    indices <- which(upper.tri(
        matrix(FALSE, length(subtype_levels), length(subtype_levels)), diag = TRUE
    ), arr.ind = TRUE)
    data.frame(
        SubType_1 = subtype_levels[indices[, 1L]],
        SubType_2 = subtype_levels[indices[, 2L]],
        stringsAsFactors = FALSE
    )
}

#' Purpose: Encode unordered subtype pairs by their positions in a subtype ordering.
#'
#' @param type_1 Character vector of first subtype labels.
#' @param type_2 Character vector of second subtype labels, parallel to `type_1`.
#' @param subtype_levels Character vector defining the reference subtype order.
#' @return Character vector of order-invariant `"rank_1::rank_2"` pair keys.
.interaction_pair_keys <- function(type_1, type_2, subtype_levels) {
    # Encode pairs by ordered subtype ranks to make endpoint order irrelevant.
    rank_1 <- match(type_1, subtype_levels)
    rank_2 <- match(type_2, subtype_levels)
    paste(
        pmin(rank_1, rank_2), pmax(rank_1, rank_2), sep = "::"
    )
}

#' Purpose: Test ROI-local subtype-pair proximity with label permutations.
#'
#' @param cell_metadata Data frame with unique `CellID`, `PID`, `ID`, `SubType`, `x`, and
#'   `y` columns; `Tissue` is retained when present.
#' @param distance_um Positive interaction-distance threshold in micrometres.
#' @param n_permutations Positive integer number of within-ROI subtype-label permutations.
#' @param seed Integer random seed supplied explicitly by the caller.
#' @return Data frame with observed subtype-pair counts and Laplace-smoothed enrichment and
#'   avoidance P values for each ROI.
calculate_spatial_interactions <- function(
    cell_metadata, distance_um = 22, n_permutations = 1000L, seed
) {
    # Fix the spatial threshold, permutation count, and reproducible randomization state.
    if (missing(seed)) {
        stop("seed must be supplied explicitly")
    }
    n_permutations <- as.integer(n_permutations)
    seed <- as.integer(seed)
    if (length(distance_um) != 1L || !is.finite(distance_um) || distance_um <= 0) {
        stop("distance_um must be one positive number")
    }
    if (length(n_permutations) != 1L || is.na(n_permutations) ||
        n_permutations < 1L) {
        stop("n_permutations must be one positive integer")
    }
    if (length(seed) != 1L || is.na(seed)) {
        stop("seed must be one integer")
    }

    # Shape cell metadata and enforce finite micrometre coordinates for distance calculations.
    required <- c("CellID", "PID", "ID", "SubType", "x", "y")
    missing <- setdiff(required, names(cell_metadata))
    if (length(missing) > 0L) {
        stop("Missing cell_metadata fields: ", paste(missing, collapse = ", ") )
    }
    if (anyDuplicated(cell_metadata$CellID)) {
        stop("CellID values must be unique")
    }
    coordinates <- as.matrix(cell_metadata[, c("x", "y"), drop = FALSE])
    storage.mode(coordinates) <- "double"
    if (any(!is.finite(coordinates))) {
        stop("x and y coordinates must be finite micrometre values")
    }

    # Define a common unordered subtype-pair basis across all ROI-specific tests.
    subtype_levels <- sort(unique(as.character(cell_metadata$SubType)))
    subtype_levels <- subtype_levels[!is.na(subtype_levels)]
    pairs <- .interaction_subtype_pairs(subtype_levels)
    pair_levels <- .interaction_pair_keys(
        pairs$SubType_1, pairs$SubType_2, subtype_levels
    )
    roi_key <- interaction(
        cell_metadata$PID, cell_metadata$ID, drop = TRUE, lex.order = TRUE
    )
    set.seed(seed)

    results <- lapply(split(seq_len(nrow(cell_metadata)), roi_key), function(indices) {
        # Construct ROI-local proximity edges from the prespecified distance threshold.
        roi_coordinates <- coordinates[indices, , drop = FALSE]
        if (length(indices) < 2L) {
            edge_indices <- matrix(integer(), ncol = 2L)
        } else {
            distances <- as.matrix(stats::dist(roi_coordinates))
            edge_indices <- which(
                upper.tri(distances) & distances <= distance_um, arr.ind = TRUE
            )
        }
        roi_subtypes <- as.character(cell_metadata$SubType[indices])

        #' Purpose: Count ROI edge endpoints by unordered subtype-pair level.
        #'
        #' @param types Character subtype labels parallel to the enclosing ROI cells.
        #' @return Integer vector of counts aligned to the enclosing `pair_levels` vector.
        count_pairs <- function(types) {
            # Count edge labels on the common pair basis, including zero-count pairs.
            if (nrow(edge_indices) == 0L) {
                return(integer(length(pair_levels)))
            }
            keys <- .interaction_pair_keys(
                types[edge_indices[, 1L]], types[edge_indices[, 2L]], subtype_levels
            )
            as.integer(table(factor(keys, levels = pair_levels)))
        }

        # Hold graph topology fixed while permuting subtype labels within the ROI.
        observed <- count_pairs(roi_subtypes)
        permuted <- matrix(
            replicate(
                n_permutations, count_pairs(sample(roi_subtypes, replace = FALSE))
            ),
            nrow = length(pair_levels), ncol = n_permutations
        )

        # Report observed counts and smoothed one-sided permutation probabilities.
        data.frame(
            Tissue = if ("Tissue" %in% names(cell_metadata)) {
                unique(cell_metadata$Tissue[indices])
            } else {
                NA_character_
            },
            PID = cell_metadata$PID[indices[[1L]]],
            ID = cell_metadata$ID[indices[[1L]]],
            pairs,
            observed_count = observed,
            p_enrichment = (rowSums(permuted >= observed) + 1) /
                (n_permutations + 1),
            p_avoidance = (rowSums(permuted <= observed) + 1) /
                (n_permutations + 1),
            stringsAsFactors = FALSE
        )
    })

    # Combine ROI results without pooling nested cell-level observations across ROIs.
    do.call(rbind, results)
}

#' Purpose: Read and validate one square saved interaction matrix.
#'
#' @param path CSV path with identical row and column subtype names.
#' @return Numeric square matrix with identical row and column names.
.read_interaction_matrix <- function(path) {
    # Read one numeric matrix while preserving subtype labels on both axes.
    if (!file.exists(path)) {
        stop("Missing expected interaction result: ", path)
    }
    matrix_data <- as.matrix(utils::read.csv(
        path, row.names = 1L, check.names = FALSE
    ))
    storage.mode(matrix_data) <- "double"
    if (nrow(matrix_data) != ncol(matrix_data) ||
        !identical(rownames(matrix_data), colnames(matrix_data))) {
        stop("Interaction result must be a square, identically named matrix: ", path)
    }
    matrix_data
}

#' Purpose: Read saved ROI interaction matrices into long-form permutation results.
#'
#' @param root Root directory containing tissue-specific interaction-result directories.
#' @param tissues Character vector of tissue names to read.
#' @param roi_ids Named list mapping each tissue to ROI identifiers.
#' @param tissue_directory Function mapping a tissue name to its directory name.
#' @param enrichment_file Relative path to an enrichment P-value matrix.
#' @param avoidance_file Relative path to an avoidance P-value matrix.
#' @param observed_file Relative path to an observed-count matrix.
#' @return Data frame with tissue, ROI, subtype-pair, observed-count, enrichment-P-value,
#'   and avoidance-P-value columns.
read_interaction_permutation_results <- function(
    root, tissues, roi_ids,
    tissue_directory = function(tissue) paste0("permutation_", tissue),
    enrichment_file = file.path("ClosePvalue", "ClosePvalue.csv"),
    avoidance_file = file.path("AvoidPvalue", "AvoidPvalue.csv"),
    observed_file = file.path("TureInteraction", "TureInteraction.csv")
) {
    # Require an explicit tissue-to-ROI roster for deterministic file discovery.
    if (!is.list(roi_ids) || is.null(names(roi_ids))) {
        stop("roi_ids must be a named list keyed by tissue")
    }
    missing_tissues <- setdiff(tissues, names(roi_ids))
    if (length(missing_tissues) > 0L) {
        stop("roi_ids is missing tissues: ", paste(missing_tissues, collapse = ", ") )
    }

    output <- list()
    output_index <- 1L
    for (tissue in tissues) {
        for (roi_id in roi_ids[[tissue]]) {
            # Load the three matrices for one ROI and require a shared subtype coordinate system.
            roi_root <- file.path(root, tissue_directory(tissue), roi_id)
            enrichment <- .read_interaction_matrix(file.path(roi_root, enrichment_file))
            avoidance <- .read_interaction_matrix(file.path(roi_root, avoidance_file))
            observed <- .read_interaction_matrix(file.path(roi_root, observed_file))
            if (!identical(dimnames(enrichment), dimnames(avoidance)) ||
                !identical(dimnames(enrichment), dimnames(observed))) {
                stop("Saved interaction matrices have different subtype axes for ", roi_id)
            }

            # Convert the upper triangle to one record per unordered subtype pair.
            indices <- which(upper.tri(enrichment, diag = TRUE), arr.ind = TRUE)
            output[[output_index]] <- data.frame(
                Tissue = tissue,
                ID = roi_id,
                SubType_1 = rownames(enrichment)[indices[, 1L]],
                SubType_2 = colnames(enrichment)[indices[, 2L]],
                observed_count = observed[indices],
                p_enrichment = enrichment[indices],
                p_avoidance = avoidance[indices],
                stringsAsFactors = FALSE
            )
            output_index <- output_index + 1L
        }
    }
    # Concatenate ROI records while retaining tissue and ROI provenance.
    do.call(rbind, output)
}

#' Purpose: Summarize significant enrichment and avoidance by tissue and subtype pair.
#'
#' @param results Long-form ROI interaction results with tissue, ROI, subtype-pair, and
#'   enrichment and avoidance P-value columns.
#' @param alpha Scalar significance threshold strictly between zero and one.
#' @return Data frame with one interaction and one avoidance row per tissue/subtype pair,
#'   including significant-ROI counts and signed normalized frequencies.
summarize_interaction_frequency <- function(results, alpha = 0.05) {
    # Require pair-level ROI results and a valid significance threshold.
    required <- c(
        "Tissue", "ID", "SubType_1", "SubType_2",
        "p_enrichment", "p_avoidance"
    )
    missing <- setdiff(required, names(results))
    if (length(missing) > 0L) {
        stop("Missing interaction result fields: ", paste(missing, collapse = ", ") )
    }
    if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
        stop("alpha must be between 0 and 1")
    }

    # Canonicalize pair orientation before enforcing one result per tissue, ROI, and pair.
    subtype_1 <- as.character(results$SubType_1)
    subtype_2 <- as.character(results$SubType_2)
    results$SubType_1 <- pmin(subtype_1, subtype_2)
    results$SubType_2 <- pmax(subtype_1, subtype_2)
    result_key <- results[c("Tissue", "ID", "SubType_1", "SubType_2")]
    if (anyDuplicated(result_key)) {
        stop("Interaction results must contain one row per tissue, ROI, and subtype pair")
    }

    # Use the number of distinct ROIs per tissue as the frequency denominator.
    roi_totals <- vapply(
        split(results$ID, results$Tissue),
        function(ids) length(unique(ids)),
        integer(1L)
    )
    pair_key <- interaction(
        results$Tissue, results$SubType_1, results$SubType_2,
        drop = TRUE, lex.order = TRUE
    )
    summaries <- lapply(split(results, pair_key), function(pair_results) {
        # Count significant enrichment and avoidance across ROIs for one tissue-pair unit.
        n_rois <- roi_totals[[as.character(pair_results$Tissue[[1L]])]]
        n_interaction <- sum(pair_results$p_enrichment <= alpha, na.rm = TRUE)
        n_avoidance <- sum(pair_results$p_avoidance <= alpha, na.rm = TRUE)
        data.frame(
            Tissue = pair_results$Tissue[[1L]],
            SubType_1 = pair_results$SubType_1[[1L]],
            SubType_2 = pair_results$SubType_2[[1L]],
            direction = c("interaction", "avoidance"),
            R_sig = c(n_interaction, n_avoidance),
            R_total = n_rois,
            normalized_frequency = c(n_interaction, n_avoidance) / n_rois,
            signed_frequency = c(n_interaction, -n_avoidance) / n_rois,
            stringsAsFactors = FALSE
        )
    })
    # Combine tissue-pair summaries with avoidance encoded on the negative axis.
    do.call(rbind, summaries)
}

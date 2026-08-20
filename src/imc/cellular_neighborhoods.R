# Module: IMC cellular-neighborhood composition.
# Purpose: Calculate ROI-local subtype-neighborhood compositions.
# Callers: `figures/Figure4/01_methods_aligned_orchestration.R` calls
# `calculate_cn_composition()`.
# Inputs: A filtered `SingleCellExperiment` with the required IMC `colData()` fields.
# Outputs: A cell-level subtype-composition data frame keyed by `CellID`, `PID`, and `ID`.
# Ordered use: Source `coordinates.R`, source this module, then call
# `calculate_cn_composition()`.

#' Purpose: Calculate unweighted subtype proportions among each cell's ROI-local nearest neighbors.
#'
#' @param sce Filtered `SingleCellExperiment` with `CellID`, `PID`, `ID`, `Tissue`,
#'   `SubType`, and `Position` in `colData()`.
#' @param k Positive integer number of nearest neighbors per cell.
#' @return Data frame with `CellID`, `PID`, `ID`, and one subtype-proportion column per
#'   observed subtype.
calculate_cn_composition <- function(sce, k = 20L) {
    # Establish a single fixed neighborhood size and the required spatial dependencies.
    k <- as.integer(k)
    if (length(k) != 1L || is.na(k) || k < 1L) {
        stop("k must be one positive integer")
    }
    if (!exists("position_to_xy", mode = "function")) {
        stop("source coordinates.R before calculating cellular neighborhoods")
    }
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
        !requireNamespace("FNN", quietly = TRUE)) {
        stop("SingleCellExperiment and FNN are required")
    }

    # Shape cell metadata around the canonical identifiers, tissue labels, and positions.
    metadata <- as.data.frame(SingleCellExperiment::colData(sce))
    required <- c("CellID", "PID", "ID", "Tissue", "SubType", "Position")
    missing <- setdiff(required, names(metadata))
    if (length(missing) > 0L) {
        stop("Missing colData fields: ", paste(missing, collapse = ", ") )
    }
    if (!all(metadata$Tissue %in% c("IM", "PT", "TC"))) {
        stop("sce must be filtered to IM, PT, and TC cells")
    }
    if (anyDuplicated(metadata$CellID)) {
        stop("CellID values must be unique")
    }

    # Derive numeric coordinates and a stable subtype basis shared across all ROIs.
    coordinates <- position_to_xy(as.character(metadata$Position))
    subtype_levels <- if (is.factor(metadata$SubType)) {
        levels(metadata$SubType)
    } else {
        sort(unique(as.character(metadata$SubType)))
    }
    subtype_levels <- subtype_levels[!is.na(subtype_levels)]
    if (length(subtype_levels) == 0L) {
        stop("SubType contains no values")
    }
    if (any(subtype_levels %in% c("CellID", "PID", "ID", "CN"))) {
        stop("SubType names conflict with identifier columns")
    }

    # Allocate one composition row per cell and keep neighborhoods within patient-ROI units.
    composition <- matrix(
        0, nrow = nrow(metadata), ncol = length(subtype_levels),
        dimnames = list(NULL, subtype_levels)
    )
    roi_key <- interaction(metadata$PID, metadata$ID, drop = TRUE, lex.order = TRUE)

    for (indices in split(seq_len(nrow(metadata)), roi_key)) {
        # Require enough ROI-local cells to exclude the index cell from a full k-neighbor set.
        if (length(indices) <= k) {
            stop(
                "ROI ", metadata$ID[indices[[1L]]], " for patient ",
                metadata$PID[indices[[1L]]], " has fewer than k + 1 cells"
            )
        }
        # Compute unweighted nearest-neighbor subtype fractions on the common subtype basis.
        neighbors <- FNN::get.knn(
            as.matrix(coordinates[indices, , drop = FALSE]), k = k
        )$nn.index
        roi_subtypes <- as.character(metadata$SubType[indices])

        for (row in seq_along(indices)) {
            counts <- table(factor(
                roi_subtypes[neighbors[row, ]], levels = subtype_levels
            ))
            composition[indices[[row]], ] <- as.numeric(counts) / k
        }
    }

    # Attach compositions to canonical cell, patient, and ROI identifiers without reordering.
    data.frame(
        CellID = metadata$CellID,
        PID = metadata$PID,
        ID = metadata$ID,
        composition,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

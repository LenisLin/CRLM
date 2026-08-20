# Module: Cholangiocyte-associated peritumoral microenvironment (PTME) helpers.
# Purpose: Delineate PTMEs from PT-cell spatial graphs, classify their niche composition,
# and summarize IMC ISR across eligible ROI or patient units.
# Callers: `figures/Figure4/01_methods_aligned_orchestration.R` calls the three public
# functions; the private helpers are called internally during PTME delineation.
# Inputs: PT-filtered IMC `SingleCellExperiment` objects, logical cholangiocyte masks,
# PTME membership tables, subtype metadata, and complete eligible-unit rosters.
# Outputs: PTME membership and classification tables plus ROI-level or patient-level ISR tables.
# Ordered use: Source `coordinates.R`, source this module, delineate PTMEs, classify them,
# then summarize ISR against the caller-supplied eligible roster.

#' Purpose: Derive unique undirected Delaunay-triangulation edges for coordinate rows.
#'
#' @param coordinates Numeric two-column coordinate matrix.
#' @return Integer matrix of unique sorted row-index pairs defining Delaunay edges.
.ptme_delaunay_edges <- function(coordinates) {
    # Connect all three vertices directly when exactly three coordinate rows are supplied.
    if (nrow(coordinates) == 3L) {
        return(matrix(c(1L, 2L, 1L, 3L, 2L, 3L), ncol = 2L, byrow = TRUE))
    }

    # Convert triangle faces to a unique undirected edge set for graph construction.
    triangles <- geometry::delaunayn(coordinates, options = "QJ")
    edges <- rbind(
        triangles[, c(1L, 2L), drop = FALSE],
        triangles[, c(1L, 3L), drop = FALSE],
        triangles[, c(2L, 3L), drop = FALSE]
    )
    unique(t(apply(edges, 1L, sort)))
}

#' Purpose: Identify points inside a convex hull or within its expansion distance.
#'
#' @param points Numeric two-column matrix of points to classify.
#' @param hull_points Numeric two-column matrix of convex-hull vertices in boundary order.
#' @param expansion_um Non-negative distance extending the hull boundary.
#' @return Logical vector indicating whether each row of `points` is inside or within the
#'   expanded hull.
.points_in_expanded_hull <- function(points, hull_points, expansion_um) {
    # Orient each point against every hull edge to identify points inside the convex polygon.
    next_vertex <- c(seq.int(2L, nrow(hull_points)), 1L)
    cross_products <- vapply(seq_len(nrow(hull_points)), function(index) {
        edge <- hull_points[next_vertex[[index]], ] - hull_points[index, ]
        relative <- sweep(points, 2L, hull_points[index, ], "-")
        edge[[1L]] * relative[, 2L] - edge[[2L]] * relative[, 1L]
    }, numeric(nrow(points)))
    if (is.null(dim(cross_products))) {
        cross_products <- matrix(cross_products, nrow = nrow(points))
    }
    inside <- apply(cross_products, 1L, function(value) {
        all(value >= -sqrt(.Machine$double.eps)) ||
            all(value <= sqrt(.Machine$double.eps))
    })

    #' Purpose: Calculate each point's Euclidean distance to one line segment.
    #'
    #' @param start Numeric length-two vector for the segment start.
    #' @param end Numeric length-two vector for the segment end.
    #' @return Numeric vector of distances, one per row of the enclosing `points` matrix.
    distance_to_segment <- function(start, end) {
        # Project onto the finite segment so expansion is measured from the hull boundary.
        segment <- end - start
        relative <- sweep(points, 2L, start, "-")
        projection <- (relative[, 1L] * segment[[1L]] +
            relative[, 2L] * segment[[2L]]) / sum(segment^2)
        projection <- pmin(1, pmax(0, projection))
        closest <- sweep(projection * matrix(
            segment, nrow = nrow(points), ncol = 2L, byrow = TRUE
        ), 2L, start, "+")
        sqrt(rowSums((points - closest)^2))
    }
    # Use the nearest boundary segment as the Euclidean hull-expansion distance.
    boundary_distance <- do.call(pmin, lapply(seq_len(nrow(hull_points)), function(index) {
        distance_to_segment(hull_points[index, ], hull_points[next_vertex[[index]], ])
    }))

    inside | boundary_distance <= expansion_um
}

#' Purpose: Delineate ROI-local PTMEs from cholangiocyte connectivity and expanded hulls.
#'
#' @param sce PT-only `SingleCellExperiment` with `CellID`, `PID`, `ID`, `Tissue`, and
#'   `Position` in `colData()`.
#' @param cholangiocyte_mask Complete logical vector identifying cholangiocyte cells in `sce`.
#' @param min_cholangiocytes Positive integer minimum connected-component size to retain.
#' @param expansion_um Non-negative convex-hull expansion in micrometres.
#' @return Data frame of PTME member cells with `CellID`, `PID`, `ID`, and `ptme_id`;
#'   returns a typed zero-row table when no PTME is retained.
delineate_cholangiocyte_ptmes <- function(
    sce, cholangiocyte_mask, min_cholangiocytes = 5L, expansion_um = 22
) {
    # Establish spatial dependencies and biologically meaningful component-size parameters.
    if (!exists("position_to_xy", mode = "function")) {
        stop("source coordinates.R before delineating PTMEs")
    }
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
        !requireNamespace("igraph", quietly = TRUE) ||
        !requireNamespace("geometry", quietly = TRUE)) {
        stop("SingleCellExperiment, igraph, and geometry are required")
    }
    min_cholangiocytes <- as.integer(min_cholangiocytes)
    if (length(min_cholangiocytes) != 1L || is.na(min_cholangiocytes) ||
        min_cholangiocytes < 3L) {
        stop("min_cholangiocytes must be one integer of at least 3")
    }
    if (length(expansion_um) != 1L || !is.finite(expansion_um) || expansion_um < 0) {
        stop("expansion_um must be one non-negative number")
    }

    # Shape PT-only cell metadata and align the cholangiocyte mask to unique cells.
    metadata <- as.data.frame(SingleCellExperiment::colData(sce))
    required <- c("CellID", "PID", "ID", "Tissue", "Position")
    missing <- setdiff(required, names(metadata))
    if (length(missing) > 0L) {
        stop("Missing colData fields: ", paste(missing, collapse = ", ") )
    }
    if (!all(metadata$Tissue == "PT")) {
        stop("sce must contain PT cells only")
    }
    if (!is.logical(cholangiocyte_mask) || length(cholangiocyte_mask) != nrow(metadata) ||
        anyNA(cholangiocyte_mask)) {
        stop("cholangiocyte_mask must be a complete logical vector with one value per cell")
    }
    if (anyDuplicated(metadata$CellID)) {
        stop("CellID values must be unique")
    }

    # Derive micrometre coordinates and isolate graph construction within patient-ROI units.
    coordinates <- as.matrix(position_to_xy(as.character(metadata$Position)))
    roi_key <- interaction(metadata$PID, metadata$ID, drop = TRUE, lex.order = TRUE)
    memberships <- list()
    membership_index <- 1L

    for (indices in split(seq_len(nrow(metadata)), roi_key)) {
        # Identify the cholangiocyte core eligible to define PTME components in this ROI.
        roi_coordinates <- coordinates[indices, , drop = FALSE]
        roi_core <- which(cholangiocyte_mask[indices])
        if (length(roi_core) < min_cholangiocytes) {
            next
        }
        core_coordinates <- roi_coordinates[roi_core, , drop = FALSE]
        if (nrow(unique(core_coordinates)) < 3L) {
            stop("Delaunay triangulation requires at least three unique cholangiocyte coordinates")
        }

        # Build the full-cell Delaunay graph, then retain connected cholangiocyte components.
        edges <- .ptme_delaunay_edges(roi_coordinates)
        graph <- igraph::make_empty_graph(n = nrow(roi_coordinates), directed = FALSE)
        igraph::V(graph)$name <- as.character(seq_len(nrow(roi_coordinates)))
        graph <- igraph::add_edges(graph, as.vector(t(edges)))
        core_graph <- igraph::induced_subgraph(graph, vids = as.character(roi_core))
        components <- igraph::components(core_graph)$membership
        retained <- names(which(table(components) >= min_cholangiocytes))

        for (component in retained) {
            # Expand each retained component hull to assign all spatially enclosed PT cells.
            component_core <- which(components == as.integer(component))
            core_local <- as.integer(igraph::V(core_graph)$name[component_core])
            component_coordinates <- roi_coordinates[core_local, , drop = FALSE]
            hull_index <- grDevices::chull(component_coordinates)
            if (length(hull_index) < 3L) {
                stop("A retained cholangiocyte component has no two-dimensional convex hull")
            }
            member_local <- which(.points_in_expanded_hull(
                roi_coordinates,
                component_coordinates[hull_index, , drop = FALSE],
                expansion_um
            ))
            global_indices <- indices[member_local]
            ptme_id <- paste(
                metadata$PID[indices[[1L]]], metadata$ID[indices[[1L]]],
                sprintf("PTME%03d", match(component, retained)), sep = "::"
            )
            memberships[[membership_index]] <- data.frame(
                CellID = metadata$CellID[global_indices],
                PID = metadata$PID[global_indices],
                ID = metadata$ID[global_indices],
                ptme_id = ptme_id,
                stringsAsFactors = FALSE
            )
            membership_index <- membership_index + 1L
        }
    }

    # Preserve the membership schema even when no ROI contains a retained PTME.
    if (length(memberships) == 0L) {
        return(data.frame(
            CellID = metadata$CellID[FALSE], PID = metadata$PID[FALSE],
            ID = metadata$ID[FALSE], ptme_id = character(),
            stringsAsFactors = FALSE
        ))
    }
    # Concatenate PTME memberships while retaining cell, patient, and ROI provenance.
    rownames_result <- do.call(rbind, memberships)
    rownames(rownames_result) <- NULL
    rownames_result
}

#' Purpose: Classify PTMEs by the stromal-to-immune niche ratio.
#'
#' @param membership PTME membership data frame with `CellID`, `PID`, `ID`, and `ptme_id`.
#' @param cell_metadata Cell metadata data frame with `CellID`, `PID`, `ID`, and `SubType`.
#' @param threshold Positive scalar ratio separating PIR and PSM classifications.
#' @param immune_subtypes Character vector of subtype labels counted as immune cells.
#' @param stromal_subtypes Character vector of subtype labels counted as stromal cells.
#' @return Data frame with one row per PTME, subtype counts and fractions, `R_niche`, and
#'   `niche_class`; equal or undefined ratios are labeled `Unclassified`.
classify_ptmes <- function(
    membership, cell_metadata, threshold = 1,
    immune_subtypes = c("B", "CD8T"),
    stromal_subtypes = c("Macro_CD163", "SC_COLLAGEN")
) {
    # Require the identifiers and subtype fields needed to align cells with PTME membership.
    required_membership <- c("CellID", "PID", "ID", "ptme_id")
    missing <- setdiff(required_membership, names(membership))
    if (length(missing) > 0L) {
        stop("Missing membership fields: ", paste(missing, collapse = ", ") )
    }
    required_metadata <- c("CellID", "PID", "ID", "SubType")
    missing <- setdiff(required_metadata, names(cell_metadata))
    if (length(missing) > 0L) {
        stop("Missing cell_metadata fields: ", paste(missing, collapse = ", ") )
    }
    if (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0) {
        stop("threshold must be one positive number")
    }
    if (anyDuplicated(cell_metadata$CellID)) {
        stop("cell_metadata CellID values must be unique")
    }
    if (anyDuplicated(membership[c("CellID", "ptme_id")])) {
        stop("Each cell may occur at most once within a PTME")
    }

    # Transfer subtype labels by CellID while preserving patient and ROI consistency.
    metadata_index <- match(membership$CellID, cell_metadata$CellID)
    if (anyNA(metadata_index)) {
        stop("membership contains CellID values absent from cell_metadata")
    }
    if (any(membership$PID != cell_metadata$PID[metadata_index]) ||
        any(membership$ID != cell_metadata$ID[metadata_index])) {
        stop("membership and cell_metadata disagree on PID or ID")
    }
    membership$SubType <- as.character(cell_metadata$SubType[metadata_index])

    # Preserve the one-row-per-PTME output contract for an empty membership table.
    if (nrow(membership) == 0L) {
        return(data.frame(
            ptme_id = character(), PID = character(), ID = character(),
            n_cells = integer(), n_immune = integer(), n_stromal = integer(),
            immune_fraction = numeric(), stromal_fraction = numeric(),
            R_niche = numeric(), niche_class = character(),
            stringsAsFactors = FALSE
        ))
    }

    classified <- lapply(split(membership, membership$ptme_id), function(ptme) {
        # Summarize immune and stromal representation within one delineated PTME.
        n_total <- nrow(ptme)
        n_immune <- sum(ptme$SubType %in% immune_subtypes)
        n_stromal <- sum(ptme$SubType %in% stromal_subtypes)
        ratio <- if (n_immune > 0L) {
            n_stromal / n_immune
        } else if (n_stromal > 0L) {
            Inf
        } else {
            NA_real_
        }
        # Classify only ratios strictly below or above the prespecified niche threshold.
        niche_class <- if (is.na(ratio) || ratio == threshold) {
            "Unclassified"
        } else if (ratio < threshold) {
            "PIR"
        } else {
            "PSM"
        }

        # Emit counts, fractions, and the niche ratio as one PTME-level record.
        data.frame(
            ptme_id = ptme$ptme_id[[1L]],
            PID = ptme$PID[[1L]],
            ID = ptme$ID[[1L]],
            n_cells = n_total,
            n_immune = n_immune,
            n_stromal = n_stromal,
            immune_fraction = n_immune / n_total,
            stromal_fraction = n_stromal / n_total,
            R_niche = ratio,
            niche_class = niche_class,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, classified)
}

#' Purpose: Summarize PIR and PSM counts and clipped IMC ISR for eligible units.
#'
#' @param ptmes One-row-per-PTME classification data frame containing `ptme_id`, `PID`,
#'   `ID`, and `niche_class`.
#' @param roster Data frame of unique complete eligible ROI or patient units.
#' @param level Unit level, either `"ROI"` or `"patient"`.
#' @param lower Positive lower bound for clipped ISR.
#' @param upper Upper bound for clipped ISR, greater than or equal to `lower`.
#' @return Data frame containing each roster unit, `N_PIR`, `N_PSM`, `ISR_raw`, and
#'   clipped `ISR`.
summarize_imc_isr <- function(
    ptmes, roster, level = c("ROI", "patient"), lower = 0.1, upper = 10
) {
    # Fix the summarization unit and require one classified record per PTME.
    level <- match.arg(level)
    required <- c("ptme_id", "PID", "ID", "niche_class")
    missing <- setdiff(required, names(ptmes))
    if (length(missing) > 0L) {
        stop("Missing PTME fields: ", paste(missing, collapse = ", ") )
    }
    if (anyDuplicated(ptmes$ptme_id)) {
        stop("ptmes must contain one row per ptme_id")
    }
    if (length(lower) != 1L || length(upper) != 1L ||
        !is.finite(lower) || !is.finite(upper) || lower <= 0 || lower > upper) {
        stop("lower and upper must define a positive ordered interval")
    }

    # Define the complete eligible-unit denominator independently of observed PTMEs.
    roster_fields <- if (level == "ROI") c("PID", "ID") else "PID"
    missing <- setdiff(roster_fields, names(roster))
    if (length(missing) > 0L) {
        stop("Missing roster fields: ", paste(missing, collapse = ", ") )
    }
    roster <- roster[, roster_fields, drop = FALSE]
    if (nrow(roster) == 0L ||
        any(!stats::complete.cases(roster)) || anyDuplicated(roster)) {
        stop("roster must contain unique, complete eligible units")
    }

    # Map every PTME to its eligible roster unit using character-stable composite keys.
    roster_units <- data.frame(
        lapply(roster, as.character), check.names = FALSE,
        stringsAsFactors = FALSE
    )
    ptme_units <- data.frame(
        lapply(ptmes[, roster_fields, drop = FALSE], as.character),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    combined_units <- rbind(roster_units, ptme_units)
    combined_key <- do.call(interaction, c(
        combined_units, list(drop = TRUE, lex.order = TRUE)
    ))
    roster_key <- combined_key[seq_len(nrow(roster))]
    ptme_key <- combined_key[nrow(roster) + seq_len(nrow(ptmes))]
    roster_index <- match(ptme_key, roster_key)
    if (anyNA(roster_index)) {
        stop("ptmes contains a unit absent from the eligible roster")
    }

    # Count classified niches per unit and calculate the pseudocount-adjusted ISR.
    n_pir <- tabulate(
        roster_index[ptmes$niche_class == "PIR"], nbins = nrow(roster)
    )
    n_psm <- tabulate(
        roster_index[ptmes$niche_class == "PSM"], nbins = nrow(roster)
    )
    raw_isr <- (n_pir + 1) / (n_psm + 1)
    # Return the requested unit identifiers with ISR constrained to the analysis bounds.
    result <- if (level == "ROI") {
        data.frame(
            PID = as.character(roster$PID), ID = as.character(roster$ID),
            N_PIR = n_pir, N_PSM = n_psm, ISR_raw = raw_isr,
            ISR = pmin(upper, pmax(lower, raw_isr)),
            stringsAsFactors = FALSE
        )
    } else {
        data.frame(
            PID = as.character(roster$PID), N_PIR = n_pir, N_PSM = n_psm,
            ISR_raw = raw_isr, ISR = pmin(upper, pmax(lower, raw_isr)),
            stringsAsFactors = FALSE
        )
    }
    rownames(result) <- NULL
    result
}

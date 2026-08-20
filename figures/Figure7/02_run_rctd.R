# Target: Figure 7 RCTD deconvolution outputs.
# Purpose: deconvolve FDZS-4 Bin100 spots with the annotated FDZS-3 reference.
# Inputs: edit the paths in CONFIG for the scRNA bundle, ST bundle, helpers, and output.
# Outputs: normalized RCTD weights, spot-status records, reference counts, and run summary TSVs.
# Ordered workflow: load helpers and bundles, validate spots, balance reference, run RCTD, and export tables.

EXPECTED_ALL_SPOTS <- 18186L
EXPECTED_RCTD_SPOTS <- 17902L

# Edit these paths for the Figure 7 analysis.
CONFIG <- list(
    sc_bundle = "/path/to/figure7_scrna_output/r_bundle",
    st_bundle = "/path/to/fdzs4_bin100_bundle",
    output_dir = "/path/to/figure7_rctd_output",
    st_functions = "src/st/st_functions.R",
    scrna_functions = "src/scrna/cellchat.R"
)

#' Purpose: write a data frame as a tab-separated output table.
#'
#' @param data Data frame to write.
#' @param path Character scalar destination path for the TSV file.
#' @return Invisibly returns `NULL` after writing the table.
write_tsv <- function(data, path) {
    utils::write.table(
        data,
        file = path,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
}

#' Purpose: execute the configured Figure 7 subtype RCTD workflow.
#'
#' @param config List containing scRNA, spatial, helper, and output paths.
#' @return Invisibly returns `NULL` after writing the RCTD output tables.
main <- function(config = CONFIG) {
    # Initialize helper functions and destinations before loading reference cells or spatial spots.
    output_dir <- config[["output_dir"]]
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    source(config[["st_functions"]])
    source(config[["scrna_functions"]])

    # Load the annotated scRNA reference; columns remain individual cells keyed by Sub_type.
    sc_object <- read_anndata_to_seurat(
        config[["sc_bundle"]],
        verbose = TRUE
    )
    if (!"Sub_type" %in% colnames(sc_object[[]])) {
        stop("The scRNA R bundle metadata must contain Sub_type", call. = FALSE)
    }
    sc_data <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            counts = Seurat::GetAssayData(
                sc_object,
                assay = "RNA",
                layer = "counts"
            )
        ),
        colData = S4Vectors::DataFrame(sc_object[[]])
    )

    # Load all Bin100 spots and preserve their count columns and spatial-coordinate row correspondence.
    st_data <- read_anndata_to_spe(config[["st_bundle"]], verbose = TRUE)
    if (ncol(st_data) != EXPECTED_ALL_SPOTS) {
        stop(
            "The all-spot Bin100 input contains ",
            EXPECTED_ALL_SPOTS, " spots; observed ", ncol(st_data),
            call. = FALSE
        )
    }

    # Balance cell-level subtype representation before constructing the RCTD reference matrix.
    balanced <- balance_cell_types(
        sc_data,
        cell_type_col = "Sub_type",
        k = 100,
        min_cells = 10,
        seed = 619
    )
    balanced_sc <- balanced$sc_data
    cell_types <- SummarizedExperiment::colData(balanced_sc)[["Sub_type"]]
    names(cell_types) <- colnames(balanced_sc)
    reference <- spacexr::Reference(
        SummarizedExperiment::assay(balanced_sc),
        as.factor(cell_types)
    )

    # Prepare integer spot counts and coordinates with identical spot ordering for SpatialRNA.
    counts <- Matrix::as(SummarizedExperiment::assay(st_data), "dgCMatrix")
    counts@x <- as.numeric(round(counts@x))
    coordinates <- as.data.frame(SpatialExperiment::spatialCoords(st_data))
    rownames(coordinates) <- colnames(st_data)
    puck <- spacexr::SpatialRNA(coordinates, counts)

    # Estimate subtype mixtures per Bin100 spot and normalize weights within retained spots.
    rctd <- spacexr::create.RCTD(puck, reference, max_cores = 8)
    rctd <- spacexr::run.RCTD(rctd, doublet_mode = "full")
    normalized_weights <- spacexr::normalize_weights(rctd@results$weights)

    # Canonicalize spot identities before alignment; each retained spot must occur exactly once.
    weight_ids <- .canonicalize_spot_ids(rownames(normalized_weights))
    spot_ids <- .canonicalize_spot_ids(colnames(st_data))
    if (length(weight_ids) != EXPECTED_RCTD_SPOTS) {
        stop(
            "The normalized RCTD result contains ",
            EXPECTED_RCTD_SPOTS, " spots; observed ", length(weight_ids),
            call. = FALSE
        )
    }
    if (anyDuplicated(weight_ids)) {
        stop(
            "Normalized RCTD weights contain duplicate canonical spot IDs",
            call. = FALSE
        )
    }
    if (anyDuplicated(spot_ids)) {
        stop(
            "The Bin100 input contains duplicate canonical spot IDs",
            call. = FALSE
        )
    }

    # Export one normalized subtype-composition row per retained Bin100 spot.
    weight_table <- data.frame(
        spot_id = weight_ids,
        as.data.frame(as.matrix(normalized_weights), check.names = FALSE),
        check.names = FALSE
    )
    write_tsv(
        weight_table,
        file.path(output_dir, "rctd_normalized_weights.tsv")
    )

    # Classify every input spot by RCTD retention while retaining the all-spot denominator.
    retained <- spot_ids %in% weight_ids
    spot_audit <- data.frame(
        spot_id = spot_ids,
        rctd_status = ifelse(retained, "retained", "excluded"),
        reason = ifelse(
            retained,
            "present_in_normalized_weights",
            "absent_from_normalized_weights"
        ),
        stringsAsFactors = FALSE
    )
    retained_audit_rows <- sum(spot_audit$rctd_status == "retained")
    excluded_audit_rows <- sum(spot_audit$rctd_status == "excluded")
    expected_excluded_rows <- EXPECTED_ALL_SPOTS - EXPECTED_RCTD_SPOTS
    if (retained_audit_rows != EXPECTED_RCTD_SPOTS ||
        excluded_audit_rows != expected_excluded_rows) {
        stop(
            "The RCTD spot audit must contain ", EXPECTED_RCTD_SPOTS,
            " retained and ", expected_excluded_rows,
            " excluded rows; observed ", retained_audit_rows,
            " retained and ", excluded_audit_rows, " excluded",
            call. = FALSE
        )
    }
    write_tsv(spot_audit, file.path(output_dir, "rctd_spot_audit.tsv"))

    # Summarize reference-cell balancing by subtype for interpretation of the deconvolution input.
    original_counts <- table(
        SummarizedExperiment::colData(sc_data)[["Sub_type"]]
    )
    balanced_counts <- table(
        SummarizedExperiment::colData(balanced_sc)[["Sub_type"]]
    )
    reference_counts <- data.frame(
        subtype = names(original_counts),
        input_cells = as.integer(original_counts),
        balanced_cells = as.integer(balanced_counts[names(original_counts)]),
        stringsAsFactors = FALSE
    )
    reference_counts$balanced_cells[is.na(reference_counts$balanced_cells)] <- 0L
    write_tsv(
        reference_counts,
        file.path(output_dir, "rctd_reference_subtype_counts.tsv")
    )

    # Export fixed population sizes and RCTD settings used by this run.
    run_summary <- data.frame(
        metric = c(
            "input_st_spots",
            "normalized_weight_spots",
            "excluded_st_spots",
            "reference_input_cells",
            "reference_balanced_cells",
            "reference_subtypes_retained",
            "balance_k",
            "balance_min_cells",
            "balance_seed",
            "rctd_cores",
            "rctd_mode"
        ),
        value = c(
            length(spot_ids),
            length(weight_ids),
            length(spot_ids) - length(weight_ids),
            ncol(sc_data),
            ncol(balanced_sc),
            sum(reference_counts$balanced_cells > 0),
            100,
            10,
            619,
            8,
            "full"
        ),
        stringsAsFactors = FALSE
    )
    write_tsv(run_summary, file.path(output_dir, "rctd_run_summary.tsv"))
}

if (sys.nframe() == 0L) {
    main()
}

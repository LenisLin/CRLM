# Target: Figure 3B-F cellular-neighborhood analysis.
# Purpose: derive inverse-distance weighted, within-ROI cellular neighborhoods;
# summarize their patient-level abundances; and produce the figure-source tables
# and plots. The configured analysis loads the canonical SCE.
# Inputs: an FDZS1_IMC_processed.rds SingleCellExperiment and the shared coordinate helper.
# Outputs: cellular-neighborhood assignments, neighbor fractions, composition
# matrices/heatmaps, patient-level fractions, and tissue-specific summaries/plots.
# Ordered workflow: validate dependencies and SCE metadata; convert Position to
# coordinates; derive neighborhoods; save assignments and composition; summarize
# patient fractions; then write IM and TC statistical summaries and plots.

# Edit these paths and set run to TRUE for a sequential analysis invocation.
figure3_config <- list(
    run = FALSE,
    sce_path = "PATH/TO/FDZS1_IMC_processed.rds",
    output_dir = "PATH/TO/Figure3_output",
    script_dir = "figures/Figure3"
)

#' Verify that required R packages are installed.
#'
#' @param packages Character vector of package names.
#' @return Invisibly `NULL`; stops when a required package is unavailable.
require_packages <- function(packages) {
    # Resolve all analysis dependencies before any data object is loaded.
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0L) stop("Missing required packages: ", paste(missing, collapse = ", "))
}

#' Load the coordinate-conversion interface.
#'
#' @param script_dir Directory containing this Figure 3 script.
#' @return Invisibly `NULL`; sources `coordinates.R` and requires
#'   `position_to_xy()` to be defined.
load_shared_coordinates <- function(script_dir) {
    # Source the shared Position parser relative to the public Figure 3 directory.
    helper <- normalizePath(file.path(script_dir, "..", "..", "src", "imc", "coordinates.R"), mustWork = FALSE)
    if (!file.exists(helper)) stop("Required src/imc coordinate interface is missing: ", helper)
    source(helper)
    if (!exists("position_to_xy", mode = "function")) {
        stop("src/imc coordinates.R must define position_to_xy().")
    }
}

required_cell_columns <- c("CellID", "ID", "PID", "Tissue", "Position", "SubType", "Treatment", "RFS_time", "RFS_event")

#' Validate required cell-level metadata at the SCE boundary.
#'
#' @param cell_data Data frame created from SCE `colData()`.
#' @return Invisibly `NULL`; stops for missing required columns or non-unique
#'   `CellID` values.
validate_cell_metadata <- function(cell_data) {
    # Enforce the cell-level identifier and clinical-field contract at the SCE boundary.
    missing <- setdiff(required_cell_columns, names(cell_data))
    if (length(missing) > 0L) stop("SCE colData is missing: ", paste(missing, collapse = ", "))
    if (anyDuplicated(cell_data$CellID)) stop("CellID must be unique at the SCE input boundary.")
}

#' Calculate inverse-distance weighted subtype fractions within each ROI.
#'
#' @param cell_data Cell metadata with `CellID`, `ID`, `SubType`, `x`, and `y`.
#' @param k_neighbors Number of within-ROI nearest neighbors per cell.
#' @param n_cores Number of parallel ROI workers.
#' @return Data frame of per-cell neighborhood subtype fractions and coordinates.
inverse_distance_fractions <- function(cell_data, k_neighbors = 20L, n_cores = 12L) {
    # Fix a common subtype feature axis for every cell and ROI.
    subtype_order <- sort(unique(as.character(cell_data$SubType)))
    # Preserve the patient-to-ROI hierarchy even if an ROI label is reused.
    by_roi <- split(
        cell_data,
        interaction(cell_data$PID, cell_data$ID, drop = TRUE, lex.order = TRUE)
    )

    #' Calculate weighted subtype fractions for one ROI.
    #'
    #' @param roi_data Cell metadata for one ROI.
    #' @return A per-cell fraction data frame, or `NULL` when the ROI has fewer
    #'   than `k_neighbors + 1` cells.
    worker <- function(roi_data) {
        # Retain only ROIs large enough to define the requested within-ROI neighborhood.
        if (nrow(roi_data) < k_neighbors + 1L) return(NULL)
        # Find spatial neighbors in micrometre coordinates for every cell in this ROI.
        knn <- FNN::get.knn(as.matrix(roi_data[, c("x", "y")]), k = k_neighbors)
        fractions <- matrix(0, nrow = nrow(roi_data), ncol = length(subtype_order),
            dimnames = list(NULL, subtype_order))
        # Convert neighbor labels to inverse-distance weighted subtype fractions per index cell.
        for (i in seq_len(nrow(roi_data))) {
            weights <- 1 / (knn$nn.dist[i, ] + 1e-9)
            weighted <- tapply(
                weights,
                factor(as.character(roi_data$SubType[knn$nn.index[i, ]]), levels = subtype_order),
                sum
            )
            weighted[is.na(weighted)] <- 0
            fractions[i, ] <- as.numeric(weighted) / sum(weights)
        }
        # Return cell identifiers and coordinates aligned row-wise with the fraction matrix.
        data.frame(CellID = roi_data$CellID, ID = roi_data$ID, fractions,
            x = roi_data$x, y = roi_data$y, check.names = FALSE)
    }
    # Calculate each ROI independently, then restore a single per-cell table.
    fractions <- parallel::mclapply(by_roi, worker, mc.cores = n_cores)
    do.call(rbind, Filter(Negate(is.null), fractions))
}

#' Cluster inverse-distance neighborhood profiles into cellular neighborhoods.
#'
#' @param cell_data Cell metadata with subtype and coordinate columns.
#' @param k_neighbors Number of within-ROI nearest neighbors per cell.
#' @param n_clusters Number of k-means neighborhood clusters.
#' @param n_cores Number of parallel ROI workers.
#' @return List containing `assignments` and `neighbor_fractions` data frames.
detect_cellular_neighborhoods <- function(cell_data, k_neighbors = 20L, n_clusters = 10L,
                                          n_cores = 12L) {
    # Build the per-cell neighborhood feature matrix on the shared subtype axis.
    fractions <- inverse_distance_fractions(cell_data, k_neighbors, n_cores)
    if (is.null(fractions) || nrow(fractions) < n_clusters) {
        stop("Fewer than ", n_clusters, " cells have valid within-ROI neighborhoods.")
    }
    # Fit the specified k-means model using only subtype-fraction features.
    feature_columns <- setdiff(names(fractions), c("CellID", "ID", "x", "y"))
    set.seed(42)
    cluster <- stats::kmeans(as.matrix(fractions[, feature_columns, drop = FALSE]),
        centers = n_clusters, nstart = 20, iter.max = 100)$cluster
    # Construct per-cell labels and merge clusters below the fixed cell-count threshold.
    assignments <- fractions[, c("CellID", "ID"), drop = FALSE]
    assignments$CNP20 <- paste0("CN", cluster)
    counts <- table(assignments$CNP20)
    assignments$CNP20[assignments$CNP20 %in% names(counts)[counts < 5000L]] <- "CN_Other"
    list(assignments = assignments, neighbor_fractions = fractions)
}

#' Attach cellular-neighborhood assignments to cell metadata.
#'
#' @param cell_data Cell metadata with `CellID`.
#' @param assignments Per-cell neighborhood assignments with `CellID` and `CNP20`.
#' @return `cell_data` with a `CNP20` column.
attach_cn_assignments <- function(cell_data, assignments) {
    # Match labels by unique CellID so the original cell metadata order is preserved.
    if (anyDuplicated(assignments$CellID)) stop("CN assignment CellID values must be unique.")
    index <- match(cell_data$CellID, assignments$CellID)
    cell_data$CNP20 <- assignments$CNP20[index]
    cell_data
}

#' Extract one clinical record per patient from cell metadata.
#'
#' @param cell_data Cell metadata containing patient-level clinical variables.
#' @return Data frame with `PID`, `Treatment`, `RFS_time`, and `RFS_event`.
patient_clinical_data <- function(cell_data) {
    # Collapse repeated cell metadata to one invariant clinical record per patient.
    clinical <- unique(cell_data[, c("PID", "Treatment", "RFS_time", "RFS_event"), drop = FALSE])
    duplicates <- duplicated(clinical$PID) | duplicated(clinical$PID, fromLast = TRUE)
    if (any(duplicates)) stop("Clinical variables must be unique within PID.")
    clinical
}

#' Calculate cellular-neighborhood fractions for each patient and tissue.
#'
#' @param cell_data Cell metadata with `CNP20` assignments.
#' @return Patient-tissue-neighborhood fraction data joined to clinical variables.
patient_cn_fractions <- function(cell_data) {
    # Retain assigned cells in the three tissue compartments used for patient summaries.
    selected <- cell_data[cell_data$Tissue %in% c("IM", "PT", "TC") & !is.na(cell_data$CNP20), , drop = FALSE]
    # Count cells at the patient-tissue-neighborhood analysis unit.
    counts <- dplyr::summarise(
        dplyr::group_by(
            dplyr::mutate(selected, PID = as.character(PID), CNP20 = as.character(CNP20)),
            PID, Tissue, CNP20
        ),
        cell_count_in_cn = dplyr::n(), .groups = "drop"
    )
    # Complete absent patient-neighborhood combinations within each tissue with zero counts.
    counts <- dplyr::group_modify(
        dplyr::group_by(counts, Tissue),
        ~ tidyr::complete(.x, PID, CNP20, fill = list(cell_count_in_cn = 0L))
    )
    # Normalize neighborhood counts within each patient-tissue total.
    fractions <- dplyr::mutate(
        dplyr::group_by(counts, PID, Tissue),
        total_cells_in_group = sum(cell_count_in_cn),
        fraction = cell_count_in_cn / total_cells_in_group
    )
    # Attach one clinical record to every patient-tissue-neighborhood fraction.
    dplyr::left_join(dplyr::ungroup(fractions), patient_clinical_data(selected), by = "PID")
}

#' Perform nominal Wilcoxon comparisons by treatment and neighborhood.
#'
#' @param tissue_data Patient-level cellular-neighborhood fractions for one tissue.
#' @return Data frame of eligible tests with nominal and BH-adjusted P values.
nominal_wilcoxon <- function(tissue_data) {
    # Partition patient fractions into independent treatment-neighborhood test groups.
    groups <- split(tissue_data, interaction(tissue_data$Treatment, tissue_data$CNP20, drop = TRUE))

    #' Test one treatment-neighborhood group when both recurrence groups are eligible.
    #'
    #' @param group_data Patient-level fractions for one treatment-neighborhood group.
    #' @return One-row nominal Wilcoxon result, or `NULL` when eligibility criteria
    #'   are not met.
    test_group <- function(group_data) {
        # Require both recurrence groups and the prespecified minimum patient representation.
        status_counts <- table(group_data$RFS_event)
        if (nrow(group_data) < 6L || length(status_counts) != 2L || min(status_counts) < 3L) return(NULL)
        # Compare patient-level neighborhood fractions between recurrence groups.
        test <- stats::wilcox.test(fraction ~ RFS_event, data = group_data)
        data.frame(
            Treatment = as.character(group_data$Treatment[1]), CNP20 = as.character(group_data$CNP20[1]),
            n = nrow(group_data), n_RFS_event0 = unname(status_counts["0"]), n_RFS_event1 = unname(status_counts["1"]),
            p_value = unname(test$p.value), stringsAsFactors = FALSE
        )
    }
    # Combine eligible tests and adjust across the tissue-level test family.
    rows <- lapply(groups, test_group)
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0L) return(data.frame())
    result <- do.call(rbind, rows)
    result$p_adjusted_bh <- stats::p.adjust(result$p_value, method = "BH")
    result
}

#' Plot patient cellular-neighborhood fractions by treatment and recurrence status.
#'
#' @param tissue_data Patient-level cellular-neighborhood fractions for one tissue.
#' @param wilcoxon_results Results returned by `nominal_wilcoxon()`.
#' @param tissue Tissue label used in the plot subtitle.
#' @param output_file PDF path for the saved plot.
#' @return Invisibly returns `output_file`, as returned by `ggsave()`, after writing one PDF.
plot_cn_frequency <- function(tissue_data, wilcoxon_results, tissue, output_file) {
    # Construct patient-level recurrence comparisons faceted by treatment and neighborhood.
    plot <- ggplot2::ggplot(tissue_data, ggplot2::aes(x = factor(RFS_event), y = fraction)) +
        ggplot2::geom_jitter(ggplot2::aes(color = factor(RFS_event)), width = 0.1, size = 2, alpha = 0.6) +
        ggplot2::geom_boxplot(ggplot2::aes(fill = factor(RFS_event)), alpha = 0.5, outlier.shape = NA) +
        ggplot2::facet_grid(Treatment ~ CNP20, scales = "free_y") +
        ggplot2::scale_fill_manual(values = c("0" = "#0072B2", "1" = "#D55E00")) +
        ggplot2::scale_color_manual(values = c("0" = "#0072B2", "1" = "#D55E00")) +
        ggplot2::labs(
            subtitle = paste0(tissue, " - CN Abundance by RFS Event, Grouped by Treatment"),
            x = "RFS Event (0 = No Recurrence, 1 = Recurrence)", y = "Fraction of Cells in Neighborhood"
        ) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(legend.position = "none", strip.text.y = ggplot2::element_text(angle = 0))
    # Place nominal test annotations above their corresponding treatment-neighborhood facets.
    if (nrow(wilcoxon_results) > 0L) {
        labels <- dplyr::mutate(wilcoxon_results,
            label = paste0("p = ", format.pval(p_value, digits = 2, eps = 0.001)),
            x = 1.5,
            y = vapply(split(tissue_data$fraction, interaction(tissue_data$Treatment, tissue_data$CNP20, drop = TRUE)), max, numeric(1))[interaction(wilcoxon_results$Treatment, wilcoxon_results$CNP20, drop = TRUE)] * 1.08
        )
        plot <- plot + ggplot2::geom_text(data = labels, ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3)
    }
    # Export the tissue-specific abundance panel.
    ggplot2::ggsave(output_file, plot, width = 12, height = 7.5)
}

#' Fit cellular-neighborhood fraction by treatment Cox interaction models.
#'
#' @param tissue_data Patient-level cellular-neighborhood fractions for one tissue.
#' @return Data frame of baseline and treatment-interaction hazard-ratio terms.
cox_summary <- function(tissue_data) {
    # Reshape patient-level neighborhood fractions to one row per patient for modeling.
    wide <- tidyr::pivot_wider(tissue_data[, c("PID", "CNP20", "fraction", "Treatment", "RFS_time", "RFS_event")],
        names_from = "CNP20", values_from = "fraction")
    features <- setdiff(names(wide), c("PID", "Treatment", "RFS_time", "RFS_event"))

    #' Fit and extract the prespecified Cox model for one neighborhood feature.
    #'
    #' @param feature Name of one wide-format neighborhood-fraction column.
    #' @return Two-row data frame containing the baseline and treatment-interaction
    #'   model terms.
    fit_feature <- function(feature) {
        # Fit the neighborhood main effect and its treatment interaction for one feature.
        model <- survival::coxph(stats::as.formula(paste0("survival::Surv(RFS_time, RFS_event) ~ scale(`", feature, "`) * Treatment")), data = wide)
        # Extract the baseline neighborhood and treatment-interaction estimates with confidence limits.
        summary <- summary(model)
        terms <- summary$coefficients[, "Pr(>|z|)"]
        keep <- c(1L, 3L)
        data.frame(
            CNP20 = feature, term = rownames(summary$coefficients)[keep],
            HR = summary$conf.int[keep, "exp(coef)"], Lower = summary$conf.int[keep, "lower .95"],
            Upper = summary$conf.int[keep, "upper .95"], p_value = terms[keep], row.names = NULL
        )
    }
    # Stack feature-specific model terms into a single tissue summary.
    rows <- lapply(features, fit_feature)
    do.call(rbind, rows)
}

#' Plot a forest summary of cellular-neighborhood Cox model terms.
#'
#' @param cox_results Results returned by `cox_summary()`.
#' @param tissue Tissue label used in the plot title.
#' @param output_file PDF path for the saved forest plot.
#' @return Invisibly returns `output_file`, as returned by `ggsave()`, after writing one PDF.
plot_cox_summary <- function(cox_results, tissue, output_file) {
    # Label baseline and interaction terms before arranging hazard ratios for display.
    cox_results$term_label <- ifelse(grepl(":Treatment", cox_results$term), "Interaction", "Baseline hazard")
    cox_results$label <- paste(cox_results$CNP20, cox_results$term_label, sep = " - ")
    plot <- ggplot2::ggplot(cox_results, ggplot2::aes(x = HR, y = reorder(label, HR), color = term_label)) +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower, xmax = Upper), height = 0.2) +
        ggplot2::geom_point() + ggplot2::scale_x_log10() +
        ggplot2::labs(title = paste0(tissue, " - Univariable Interaction Model Results"), x = "Hazard ratio", y = NULL, color = NULL) +
        ggplot2::theme_bw(base_size = 12)
    # Export a height-adjusted forest plot containing all tissue model terms.
    ggplot2::ggsave(output_file, plot, width = 10, height = max(8, 0.3 * nrow(cox_results) + 3))
}

#' Run the configured Figure 3B-F cellular-neighborhood workflow.
#'
#' @param config List with `sce_path`, `output_dir`, and `script_dir` entries.
#' @return Invisibly `NULL`; writes the configured Figure 3 output files.
run_figure3_cellular_neighborhoods <- function(config) {
    # Require the analysis configuration and package interfaces before loading the SCE.
    required_config <- c("sce_path", "output_dir", "script_dir")
    missing_config <- setdiff(required_config, names(config))
    if (length(missing_config) > 0L) {
        stop("Figure 3 configuration is missing: ", paste(missing_config, collapse = ", "))
    }
    require_packages(c("SingleCellExperiment", "dplyr", "tidyr", "FNN", "pheatmap", "ggplot2", "survival"))
    # Load cell metadata and convert stored Positions to numeric micrometre coordinates.
    set.seed(619)
    load_shared_coordinates(config$script_dir)
    sce <- readRDS(config$sce_path)
    cell_data <- as.data.frame(SummarizedExperiment::colData(sce))
    validate_cell_metadata(cell_data)
    coordinates <- position_to_xy(cell_data$Position)
    cell_data$x <- coordinates$x
    cell_data$y <- coordinates$y
    # Prepare the panel and intermediate-results output directories.
    dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
    results_dir <- file.path(config$output_dir, "cn_results")
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

    # Detect per-cell neighborhoods and export assignments plus their feature fractions.
    detected <- detect_cellular_neighborhoods(cell_data)
    utils::write.csv(detected$assignments, file.path(results_dir, "peritumor_cn_results.csv"), row.names = FALSE)
    utils::write.csv(detected$neighbor_fractions, file.path(results_dir, "Spatial neighbors of All with k=20 neighbors.csv"), row.names = FALSE)
    cells_with_cn <- attach_cn_assignments(cell_data, detected$assignments)
    selected_cells <- cells_with_cn[cells_with_cn$Tissue %in% c("IM", "PT", "TC") & !is.na(cells_with_cn$CNP20), , drop = FALSE]

    # Aggregate assigned cells into a neighborhood-by-subtype composition matrix.
    composition <- dplyr::summarise(dplyr::group_by(selected_cells, CNP20, SubType), n = dplyr::n(), .groups = "drop")
    composition <- dplyr::mutate(dplyr::group_by(composition, CNP20), subtype_fraction = n / sum(n))
    composition_wide <- tidyr::pivot_wider(composition[, c("CNP20", "SubType", "subtype_fraction")], names_from = "SubType", values_from = "subtype_fraction", values_fill = 0)
    composition_matrix <- as.matrix(composition_wide[, -1])
    rownames(composition_matrix) <- composition_wide$CNP20
    utils::write.table(composition_matrix, file.path(config$output_dir, "Cellular Neighbors celltype fraction matrix.tsv"), sep = "\t", quote = FALSE)
    # Construct and export the column-scaled neighborhood composition heatmap.
    heatmap <- pheatmap::pheatmap(composition_matrix, color = grDevices::colorRampPalette(c("#436eee", "white", "#EE0000"))(100), scale = "column", cluster_rows = FALSE, cluster_cols = TRUE, silent = TRUE)
    grDevices::pdf(file.path(config$output_dir, "Cellular Neighbors celltype fraction heatmap.pdf"), width = 8, height = 6)
    print(heatmap)
    grDevices::dev.off()

    # Derive and export patient-level neighborhood fractions for downstream inference.
    fractions <- patient_cn_fractions(selected_cells)
    utils::write.csv(fractions, file.path(results_dir, "patient_level_cn_fractions.csv"), row.names = FALSE)
    # Run recurrence comparisons and Cox interaction models independently for IM and TC.
    for (tissue in c("IM", "TC")) {
        tissue_data <- fractions[fractions$Tissue == tissue, , drop = FALSE]
        wilcoxon_results <- nominal_wilcoxon(tissue_data)
        # Export the statistical summaries and their corresponding figure panels.
        utils::write.csv(wilcoxon_results, file.path(results_dir, paste0(tissue, "_CN_nominal_wilcoxon.csv")), row.names = FALSE)
        plot_cn_frequency(tissue_data, wilcoxon_results, tissue, file.path(config$output_dir, paste0(tissue, "_CN_Frequency_by_treatment_group.pdf")))
        cox_results <- cox_summary(tissue_data)
        utils::write.csv(cox_results, file.path(results_dir, paste0(tissue, "_CN_cox_results.csv")), row.names = FALSE)
        plot_cox_summary(cox_results, tissue, file.path(config$output_dir, paste0("ForestPlot_Univariable_CN_Fraction_", tissue, ".pdf")))
    }
}

if (isTRUE(figure3_config$run)) {
    run_figure3_cellular_neighborhoods(figure3_config)
}

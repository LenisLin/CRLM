# Target: Supplementary Figure 6A-C pairwise cellular-interaction frequency panels.
# Purpose: aggregate the specified saved ClosePvalue.csv matrices into
# tissue-specific interaction-frequency summaries.
# Inputs: saved permutation-result directories containing the explicit IM, PT,
# and TC ROI rosters and their ClosePvalue.csv matrices.
# Outputs: a saved-result roster table plus one frequency matrix TSV and heatmap
# PDF for each tissue panel.
# Ordered workflow: validate the fixed rosters and dependencies; preflight every
# saved p-value matrix; record the roster; compute signed frequencies; then write
# the tissue-specific matrices and heatmaps.

# Edit these paths and set run to TRUE for a sequential analysis invocation.
figure_s6_config <- list(
    run = FALSE,
    permutation_root = "PATH/TO/saved_permutation_results",
    output_dir = "PATH/TO/FigureS6_output"
)

celltype_order <- c(
    "Macro_CD169", "Macro_HLADR", "Mono_CD11c", "SC_COLLAGEN", "CD8T", "CD4T", "B",
    "SC_FAP", "Macro_CD163", "UNKNOWN", "SC_Vimentin", "Mono_Intermediate", "Mono_Classic",
    "TC_EpCAM", "SC_aSMA", "TC_Ki67", "NK", "Treg", "Macro_CD11b", "TC_CAIX", "TC_VEGF"
)

# The specified saved-result roster fixes the denominator for each panel.
# A=IM, B=PT, and C=TC.
roi_rosters <- list(
    IM = c(
    "B10_ROI10",
    "B10_ROI15",
    "B10_ROI5",
    "B11_ROI11",
    "B11_ROI13",
    "B11_ROI6",
    "B12_ROI11",
    "B12_ROI12",
    "B12_ROI6",
    "B12_ROI9",
    "B13_ROI7",
    "B13_ROI8",
    "B13_ROI9",
    "B14_ROI10",
    "B14_ROI11",
    "B14_ROI6",
    "B14_ROI9",
    "B15_ROI13",
    "B15_ROI7",
    "B15_ROI8",
    "B16_ROI12",
    "B16_ROI13",
    "B16_ROI5",
    "B16_ROI9",
    "B18_ROI11",
    "B18_ROI12",
    "B18_ROI6",
    "B18_ROI9",
    "B20_ROI10",
    "B20_ROI11",
    "B20_ROI12",
    "B20_ROI9",
    "B21_ROI12",
    "B21_ROI6",
    "B21_ROI7",
    "B21_ROI9",
    "B22_ROI10",
    "B22_ROI13",
    "B22_ROI7",
    "B22_ROI9",
    "B2_ROI4",
    "B2_ROI6",
    "B2_ROI9",
    "B3_ROI11",
    "B3_ROI12",
    "B3_ROI9",
    "B4_ROI10",
    "B4_ROI6",
    "B4_ROI8",
    "B5_ROI14",
    "B5_ROI15",
    "B5_ROI8",
    "B6_ROI12",
    "B6_ROI7",
    "B6_ROI8",
    "B8_ROI12",
    "B8_ROI5",
    "B8_ROI9",
    "W10_ROI10",
    "W10_ROI11",
    "W10_ROI7",
    "W12_ROI12",
    "W12_ROI13",
    "W12_ROI15",
    "W15_ROI10",
    "W15_ROI12",
    "W15_ROI7",
    "W16_ROI12",
    "W16_ROI5",
    "W16_ROI8",
    "W16_ROI9",
    "W17_ROI10",
    "W17_ROI11",
    "W17_ROI6",
    "W17_ROI8",
    "W18_ROI11",
    "W18_ROI13",
    "W18_ROI14",
    "W19_ROI10",
    "W19_ROI12",
    "W19_ROI13",
    "W19_ROI8",
    "W19_ROI9",
    "W1_ROI5",
    "W1_ROI6",
    "W1_ROI7",
    "W20_ROI10",
    "W20_ROI12",
    "W20_ROI7",
    "W21_ROI10",
    "W21_ROI11",
    "W21_ROI8",
    "W22_ROI10",
    "W22_ROI13",
    "W22_ROI5",
    "W23_ROI10",
    "W23_ROI7",
    "W23_ROI8",
    "W25_ROI14",
    "W25_ROI15",
    "W25_ROI6",
    "W25_ROI8",
    "W2_ROI12",
    "W2_ROI5",
    "W2_ROI7",
    "W3_ROI11",
    "W3_ROI6",
    "W3_ROI7",
    "W4_ROI11",
    "W4_ROI5",
    "W4_ROI8",
    "W6_ROI13",
    "W6_ROI5",
    "W6_ROI8",
    "W7_ROI13",
    "W7_ROI7",
    "W7_ROI9"
),
    PT = c(
    "B10_ROI7",
    "B10_ROI8",
    "B10_ROI9",
    "B11_ROI10",
    "B11_ROI5",
    "B11_ROI9",
    "B12_ROI10",
    "B12_ROI5",
    "B12_ROI8",
    "B13_ROI10",
    "B13_ROI12",
    "B13_ROI13",
    "B14_ROI5",
    "B14_ROI7",
    "B14_ROI8",
    "B15_ROI5",
    "B15_ROI6",
    "B15_ROI9",
    "B16_ROI10",
    "B16_ROI11",
    "B16_ROI6",
    "B18_ROI10",
    "B18_ROI7",
    "B18_ROI8",
    "B20_ROI5",
    "B20_ROI6",
    "B20_ROI7",
    "B21_ROI10",
    "B21_ROI5",
    "B21_ROI8",
    "B22_ROI5",
    "B22_ROI6",
    "B22_ROI8",
    "B2_ROI1",
    "B2_ROI2",
    "B2_ROI3",
    "B3_ROI5",
    "B3_ROI6",
    "B3_ROI7",
    "B3_ROI8",
    "B4_ROI11",
    "B4_ROI12",
    "B4_ROI5",
    "B5_ROI12",
    "B5_ROI5",
    "B5_ROI6",
    "B6_ROI10",
    "B6_ROI11",
    "B6_ROI9",
    "B8_ROI11",
    "B8_ROI7",
    "B8_ROI8",
    "W10_ROI5",
    "W10_ROI8",
    "W10_ROI9",
    "W12_ROI11",
    "W12_ROI16",
    "W12_ROI6",
    "W15_ROI13",
    "W15_ROI5",
    "W15_ROI6",
    "W16_ROI11",
    "W16_ROI13",
    "W16_ROI6",
    "W17_ROI12",
    "W17_ROI13",
    "W17_ROI7",
    "W18_ROI5",
    "W18_ROI9",
    "W19_ROI11",
    "W19_ROI5",
    "W19_ROI6",
    "W20_ROI11",
    "W20_ROI5",
    "W20_ROI6",
    "W21_ROI12",
    "W21_ROI13",
    "W21_ROI9",
    "W22_ROI12",
    "W22_ROI14",
    "W22_ROI7",
    "W23_ROI11",
    "W23_ROI13",
    "W23_ROI5",
    "W25_ROI10",
    "W25_ROI5",
    "W25_ROI9",
    "W2_ROI11",
    "W2_ROI13",
    "W2_ROI9",
    "W3_ROI12",
    "W3_ROI5",
    "W3_ROI8",
    "W4_ROI12",
    "W4_ROI13",
    "W4_ROI7",
    "W6_ROI10",
    "W6_ROI15",
    "W6_ROI9",
    "W7_ROI5",
    "W7_ROI6",
    "W7_ROI8"
),
    TC = c(
    "B10_ROI11",
    "B11_ROI12",
    "B11_ROI7",
    "B13_ROI11",
    "B13_ROI6",
    "B15_ROI11",
    "B15_ROI12",
    "B15_ROI15",
    "B16_ROI8",
    "B18_ROI13",
    "B20_ROI13",
    "B20_ROI14",
    "B21_ROI13",
    "B22_ROI11",
    "B2_ROI7",
    "B2_ROI8",
    "B3_ROI10",
    "B3_ROI13",
    "B3_ROI14",
    "B4_ROI13",
    "B4_ROI7",
    "B5_ROI10",
    "B5_ROI11",
    "B5_ROI9",
    "B6_ROI4",
    "B6_ROI6",
    "B8_ROI10",
    "B8_ROI13",
    "W10_ROI12",
    "W10_ROI13",
    "W10_ROI6",
    "W12_ROI14",
    "W12_ROI8",
    "W12_ROI9",
    "W15_ROI8",
    "W16_ROI10",
    "W16_ROI7",
    "W17_ROI9",
    "W18_ROI16",
    "W19_ROI7",
    "W1_ROI10",
    "W1_ROI8",
    "W1_ROI9",
    "W20_ROI13",
    "W20_ROI8",
    "W20_ROI9",
    "W21_ROI5",
    "W21_ROI6",
    "W21_ROI7",
    "W22_ROI11",
    "W22_ROI9",
    "W23_ROI12",
    "W23_ROI6",
    "W23_ROI9",
    "W25_ROI12",
    "W2_ROI10",
    "W2_ROI6",
    "W2_ROI8",
    "W3_ROI10",
    "W3_ROI13",
    "W3_ROI9",
    "W4_ROI10",
    "W4_ROI6",
    "W4_ROI9",
    "W6_ROI14",
    "W6_ROI6",
    "W7_ROI10",
    "W7_ROI11",
    "W7_ROI12"
)
)

if (!identical(vapply(roi_rosters, length, integer(1)), c(IM = 117L, PT = 102L, TC = 69L))) {
    stop("Supplementary Figure 6 roster lengths must be IM=117, PT=102, and TC=69.")
}

#' Verify that required R packages are installed.
#'
#' @param packages Character vector of package names.
#' @return Invisibly `NULL`; stops when a required package is unavailable.
require_packages <- function(packages) {
    # Resolve plotting dependencies before loading the saved ROI matrices.
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0L) stop("Missing required packages: ", paste(missing, collapse = ", "))
}

#' Read one saved ClosePvalue matrix and validate its subtype axes.
#'
#' @param permutation_root Tissue-specific saved-result root directory.
#' @param roi ROI identifier in the specified tissue roster.
#' @return Numeric p-value matrix in the required cell-type order.
read_pvalue_table <- function(permutation_root, roi) {
    # Read the saved matrix for one ROI from the established directory layout.
    table_path <- file.path(permutation_root, roi, "ClosePvalue", "ClosePvalue.csv")
    if (!file.exists(table_path)) stop("Missing saved p-value table for ROI ", roi, ": ", table_path)
    pvalues <- as.matrix(utils::read.csv(table_path, row.names = 1, check.names = FALSE))
    # Preserve an identical cell-type axis across all matrices before aggregation.
    if (!identical(rownames(pvalues), celltype_order) || !identical(colnames(pvalues), celltype_order)) {
        stop("Unexpected cell-type order in saved p-value table for ROI ", roi, ": ", table_path)
    }
    storage.mode(pvalues) <- "numeric"
    pvalues
}

#' Verify and load every saved ClosePvalue matrix in the fixed rosters.
#'
#' @param permutation_root Root directory containing `permutation_IM`,
#'   `permutation_PT`, and `permutation_TC` result directories.
#' @return Named list of tissue-specific, named ROI p-value matrices.
preflight_saved_pvalue_tables <- function(permutation_root) {
    #' List all required ClosePvalue paths for one tissue roster.
    #'
    #' @param tissue Tissue key in `roi_rosters`.
    #' @return Named character vector of expected saved p-value file paths.
    expected_paths_for_tissue <- function(tissue) {
        # Expand the fixed tissue roster into the complete set of required matrix paths.
        rosters <- roi_rosters[[tissue]]
        stats::setNames(
            file.path(permutation_root, paste0("permutation_", tissue), rosters, "ClosePvalue", "ClosePvalue.csv"),
            paste(tissue, rosters, sep = "/")
        )
    }
    # Require every roster member before any tissue-level frequency is calculated.
    table_paths <- unlist(lapply(names(roi_rosters), expected_paths_for_tissue), use.names = TRUE)
    missing_paths <- table_paths[!file.exists(table_paths)]
    if (length(missing_paths) > 0L) {
        stop("Missing required saved p-value table(s):\n", paste(paste(names(missing_paths), missing_paths, sep = ": "), collapse = "\n"))
    }

    #' Load and name every saved ClosePvalue matrix for one tissue roster.
    #'
    #' @param tissue Tissue key in `roi_rosters`.
    #' @return Named list of p-value matrices for the specified tissue roster.
    read_tissue_roster <- function(tissue) {
        # Load matrices within one tissue while preserving the fixed ROI roster names.
        permutation_tissue_root <- file.path(permutation_root, paste0("permutation_", tissue))

        #' Read one specified ROI matrix from the enclosing tissue directory.
        #'
        #' @param roi ROI identifier in the specified tissue roster.
        #' @return Numeric p-value matrix validated by `read_pvalue_table()`.
        read_roi <- function(roi) {
            # Resolve one ROI relative to its enclosing tissue directory.
            read_pvalue_table(permutation_tissue_root, roi)
        }
        stats::setNames(lapply(roi_rosters[[tissue]], read_roi), roi_rosters[[tissue]])
    }
    # Return a tissue-keyed collection with ROI-keyed matrices at the inner level.
    pvalue_tables <- lapply(names(roi_rosters), read_tissue_roster)
    stats::setNames(pvalue_tables, names(roi_rosters))
}

#' Calculate the signed saved-result interaction frequency matrix for one tissue.
#'
#' @param pvalue_tables Named list of saved ClosePvalue matrices.
#' @param roster Character vector of ROI identifiers used as the denominator.
#' @return Matrix of the existing `<= 0.05` enrichment and `> 0.95` avoidance
#'   coding divided by the fixed roster length.
frequency_matrix <- function(pvalue_tables, roster) {
    #' Encode one saved p-value matrix with the established signed thresholds.
    #'
    #' @param pvalues Numeric ClosePvalue matrix for one ROI.
    #' @return Matrix with enrichment `1`, avoidance `-1`, and neutral `0` values.
    encode_roi_pvalues <- function(pvalues) {
        # Convert each ROI matrix to the established signed enrichment/avoidance coding.
        ifelse(pvalues <= 0.05, 1, ifelse(pvalues > 0.95, -1, 0))
    }
    # Average signed indicators over the complete fixed ROI denominator for the tissue.
    roi_values <- lapply(pvalue_tables[roster], encode_roi_pvalues)
    Reduce("+", roi_values) / length(roster)
}

#' Write the fixed ROI roster used by the Supplementary Figure 6 panels.
#'
#' @param output_dir Directory receiving `FigureS6_saved_pvalue_roster.tsv`.
#' @return Invisibly `NULL`; writes the roster status table.
write_roster_status <- function(output_dir) {
    #' Create roster-status rows for one tissue panel.
    #'
    #' @param tissue Tissue key in `roi_rosters`.
    #' @return Data frame of panel index, tissue, and ROI values.
    roster_status_for_tissue <- function(tissue) {
        # Record each ROI with the tissue panel that contributes to its denominator.
        data.frame(panel_index = match(tissue, c("IM", "PT", "TC")), tissue = tissue,
            roi = roi_rosters[[tissue]], stringsAsFactors = FALSE)
    }
    # Consolidate and export the exact saved-result roster used by all three panels.
    status <- do.call(rbind, lapply(names(roi_rosters), roster_status_for_tissue))
    status$panel <- c("A", "B", "C")[status$panel_index]
    status$expected_pvalue_table <- file.path(status$tissue, status$roi, "ClosePvalue", "ClosePvalue.csv")
    utils::write.table(status, file.path(output_dir, "FigureS6_saved_pvalue_roster.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

#' Plot one Supplementary Figure 6 tissue interaction-frequency heatmap.
#'
#' @param matrix Signed interaction-frequency matrix for one tissue.
#' @param tissue Tissue label used in the panel title.
#' @param panel Supplementary Figure 6 panel label.
#' @param output_file PDF path for the saved heatmap.
#' @return Integer device number returned by `dev.off()` after writing one PDF.
plot_frequency_matrix <- function(matrix, tissue, panel, output_file) {
    # Map signed interaction frequencies to a symmetric diverging color scale.
    colour_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#2166ac", "#f7f7f7", "#b2182b"))
    # Construct the clustered tissue matrix with the panel identity and fixed legend scale.
    heatmap <- ComplexHeatmap::Heatmap(
        matrix, name = "Interaction\nFrequency", col = colour_fun,
        cluster_rows = TRUE, cluster_columns = TRUE, row_dend_reorder = TRUE, column_dend_reorder = TRUE,
        row_dend_width = grid::unit(25, "mm"), column_dend_height = grid::unit(25, "mm"),
        rect_gp = grid::gpar(col = "white", lwd = 1), show_row_names = TRUE, show_column_names = TRUE,
        column_title = paste0("Panel ", panel, ": Cell-Cell Interaction Frequency in ", tissue),
        column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
        row_names_gp = grid::gpar(fontsize = 10), column_names_gp = grid::gpar(fontsize = 10),
        column_names_rot = 45,
        heatmap_legend_param = list(title = "Frequency", at = c(-1, -0.5, 0, 0.5, 1),
            legend_height = grid::unit(6, "cm"))
    )
    # Export the tissue-specific interaction-frequency heatmap.
    grDevices::pdf(output_file, width = 10, height = 9)
    ComplexHeatmap::draw(heatmap)
    grDevices::dev.off()
}

#' Run the configured Supplementary Figure 6 saved-result workflow.
#'
#' @param config List with `permutation_root` and `output_dir` entries.
#' @return Invisibly `NULL`; writes the configured roster, matrices, and PDFs.
run_figure_s6_pairwise_interactions <- function(config) {
    # Validate saved-result and output locations before resolving plotting packages.
    required_config <- c("permutation_root", "output_dir")
    missing_config <- setdiff(required_config, names(config))
    if (length(missing_config) > 0L) {
        stop("Supplementary Figure 6 configuration is missing: ", paste(missing_config, collapse = ", "))
    }
    require_packages(c("ComplexHeatmap", "circlize"))
    permutation_root <- config$permutation_root
    output_dir <- config$output_dir
    # Load the complete fixed roster of ROI matrices before producing any outputs.
    pvalue_tables <- preflight_saved_pvalue_tables(permutation_root)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    # Export the roster that defines each tissue panel denominator.
    write_roster_status(output_dir)
    for (tissue in c("IM", "PT", "TC")) {
        # Aggregate, export, and plot one signed frequency matrix per tissue panel.
        matrix <- frequency_matrix(pvalue_tables[[tissue]], roi_rosters[[tissue]])
        utils::write.table(matrix, file.path(output_dir, paste0("Cell-Cell Interaction Frequency in ", tissue, ".tsv")),
            sep = "\t", quote = FALSE, col.names = NA)
        panel <- c(IM = "A", PT = "B", TC = "C")[[tissue]]
        plot_frequency_matrix(matrix, tissue, panel,
            file.path(output_dir, paste0("Cell-Cell Interaction Frequency in ", tissue, ".pdf")))
    }
}

if (isTRUE(figure_s6_config$run)) {
    run_figure_s6_pairwise_interactions(figure_s6_config)
}

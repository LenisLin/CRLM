#!/usr/bin/env Rscript

# Figure 4 PTME/ISR plot consumer.
#
# Purpose: Render only plot panels whose full input and statistical contracts
# are present in the TSV outputs of 01_methods_aligned_orchestration.R.
# Current coverage: Supplementary Figure 9C.

RUN_FIGURE4_PLOT_CONSUMERS <- FALSE

TABLE_DIRECTORY <- file.path("outputs", "Figure4")
OUTPUT_DIRECTORY <- file.path("outputs", "Figure4", "plots")

#' Purpose: Stop if a plotting dependency is unavailable.
#'
#' @param packages Character vector of package names.
#' @return Invisibly returns NULL; otherwise stops with missing package names.
require_packages <- function(packages) {
    missing <- packages[!vapply(
        packages, requireNamespace, logical(1L), quietly = TRUE
    )]
    if (length(missing) > 0L) {
        stop("Missing required packages: ", paste(missing, collapse = ", "))
    }
}

#' Purpose: Read and validate the per-PTME table required for Supplementary Figure 9C.
#'
#' @param table_directory Directory written by 01_methods_aligned_orchestration.R.
#' @return One row per classified PTME with a positive integer n_cells value.
read_ptme_classification <- function(table_directory) {
    path <- file.path(table_directory, "figure4_ptme_classification.tsv")
    if (!file.exists(path)) {
        stop("Missing Figure 4 PTME classification table: ", path)
    }
    ptmes <- utils::read.delim(
        path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
        check.names = FALSE
    )
    required <- c("ptme_id", "PID", "ID", "n_cells")
    missing <- setdiff(required, names(ptmes))
    if (length(missing) > 0L) {
        stop("PTME classification table is missing: ", paste(missing, collapse = ", "))
    }
    if (nrow(ptmes) == 0L || anyDuplicated(ptmes$ptme_id) ||
        anyNA(ptmes[, required, drop = FALSE])) {
        stop("PTME classification table must contain complete unique PTME rows")
    }
    ptmes$n_cells <- suppressWarnings(as.numeric(ptmes$n_cells))
    if (any(!is.finite(ptmes$n_cells)) || any(ptmes$n_cells < 1) ||
        any(ptmes$n_cells != floor(ptmes$n_cells))) {
        stop("PTME n_cells must contain positive integers")
    }
    ptmes$n_cells <- as.integer(ptmes$n_cells)
    ptmes
}

#' Purpose: Render the PTME structure-size histogram for Supplementary Figure 9C.
#'
#' @param ptmes One-row-per-PTME classification data frame.
#' @return A ggplot object with PTME size measured in cells.
plot_supplementary_figure9c <- function(ptmes) {
    ggplot2::ggplot(ptmes, ggplot2::aes(x = n_cells)) +
        ggplot2::geom_histogram(
            bins = 30L, fill = "#C73E1D", alpha = 0.8, color = "white"
        ) +
        ggplot2::labs(
            x = "PTME size (number of cells)",
            y = "Count"
        ) +
        ggplot2::theme_minimal(base_size = 12)
}

#' Purpose: Render all currently supported Figure 4 PTME/ISR plot consumers.
#'
#' @param table_directory Directory written by 01_methods_aligned_orchestration.R.
#' @param output_directory Directory for consumer PDFs.
#' @return Invisibly returns a named list of rendered ggplot objects.
run_figure4_plot_consumers <- function(table_directory, output_directory) {
    require_packages("ggplot2")
    ptmes <- read_ptme_classification(table_directory)
    supplementary_figure9c <- plot_supplementary_figure9c(ptmes)

    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
        filename = file.path(output_directory, "SupplementaryFigure9C_ptme_size_distribution.pdf"),
        plot = supplementary_figure9c, width = 8, height = 6
    )
    invisible(list(supplementary_figure9c = supplementary_figure9c))
}

if (isTRUE(RUN_FIGURE4_PLOT_CONSUMERS)) {
    run_figure4_plot_consumers(TABLE_DIRECTORY, OUTPUT_DIRECTORY)
}

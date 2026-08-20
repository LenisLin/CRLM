#!/usr/bin/env Rscript

# Target: Figure 2B cell-subpopulation marker heatmap.
# Purpose: Summarize IMC cell subpopulations and visualize scaled marker expression.
# Inputs: A SingleCellExperiment RDS with assay data, current annotations, and palettes.
# Outputs: A subtype summary TSV and Figure 2B marker-heatmap PDF in `output_dir`.
# Workflow: 1. Load and validate the SCE. 2. Retain IM/PT/TC cells. 3. Summarize
# subtype abundance and expression. 4. Write the summary and render the heatmap.

# Directly editable interactive configuration.
sce_path <- file.path("data", "FDZS1_IMC_processed.rds")
output_dir <- file.path("results", "Figure2", "figure2B")

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

#' Create the Figure 2B marker-expression heatmap and subtype summary.
#'
#' @param sce_path Character scalar path to the input SingleCellExperiment RDS.
#' @param output_dir Character scalar path to the directory for Figure 2B outputs.
#' @return Integer device number returned by `dev.off()` after writing the TSV summary and PDF heatmap.
run_figure2b <- function(sce_path, output_dir) {
  # Load the cell-level object and require the annotation fields used throughout the panel.
  if (!file.exists(sce_path)) stop("SCE input does not exist: ", sce_path)
  sce <- readRDS(sce_path)
  required_fields <- c("Tissue", "SubType", "MajorType2")
  missing_fields <- setdiff(required_fields, colnames(colData(sce)))
  if (length(missing_fields) > 0L) stop("Missing colData field(s): ", paste(missing_fields, collapse = ", "))

  # Restrict the analysis unit to cells from the three tissue compartments shown in Figure 2.
  retained <- colData(sce)$Tissue %in% c("IM", "PT", "TC")
  if (!any(retained)) stop("No IM, PT, or TC cells are available.")
  sce <- sce[, retained]
  meta <- as.data.frame(colData(sce))

  # Separate lineage-annotation markers from the remaining status-associated marker columns.
  annotation_markers <- c(
    "CD45", "CD20", "EpCAM", "CD11c", "CD14", "CD16", "HLADR", "CD68",
    "CD11b", "CD163", "CD169", "CD57", "CollagenI", "AlphaSMA", "FAP",
    "Vimentin", "CD3", "CD4", "CD8a", "FoxP3"
  )
  markers <- rownames(sce)
  missing_markers <- setdiff(annotation_markers, markers)
  if (length(missing_markers) > 0L) stop("Missing annotation marker(s): ", paste(missing_markers, collapse = ", "))

  # Use the SCE palette so subtype-to-major-type annotations remain consistent across figures.
  color_vectors <- metadata(sce)$color_vectors
  major_palette <- color_vectors$MajorType2
  if (is.null(major_palette)) stop("metadata(sce)$color_vectors$MajorType2 is required.")

  subtypes <- sort(unique(as.character(meta$SubType)))
  #' Map one subtype to its observed MajorType2 label.
  #'
  #' @param subtype Character scalar subtype label present in `meta$SubType`.
  #' @return Character scalar MajorType2 label for the subtype.
  major_types <- vapply(subtypes, function(subtype) {
    # Each subtype inherits the major-type label observed for its cells.
    as.character(meta$MajorType2[match(subtype, meta$SubType)])
  }, character(1L))
  if (any(!major_types %in% names(major_palette))) stop("MajorType2 palette is incomplete.")

  # Aggregate expression at the subtype level before marker-wise scaling for display.
  expression <- assay(sce)
  #' Calculate marker-wise mean expression for one subtype.
  #'
  #' @param subtype Character scalar subtype label present in `meta$SubType`.
  #' @return Numeric vector of mean expression values, one per SCE marker.
  mean_expression <- vapply(subtypes, function(subtype) {
    # Retain one mean per marker using all retained cells assigned to this subtype.
    rowMeans(expression[, meta$SubType == subtype, drop = FALSE])
  }, numeric(nrow(expression)))
  mean_expression <- t(mean_expression)
  colnames(mean_expression) <- markers

  #' Count retained cells assigned to one subtype.
  #'
  #' @param subtype Character scalar subtype label present in `meta$SubType`.
  #' @return Numeric scalar number of retained cells with the subtype label.
  # Pair expression summaries with retained-cell counts and cohort-wide cell fractions.
  cell_count <- vapply(subtypes, function(subtype) sum(meta$SubType == subtype), numeric(1L))
  summary_table <- data.frame(
    subtype = subtypes,
    major_type = major_types,
    cell_count = cell_count,
    cell_fraction = cell_count / ncol(sce),
    check.names = FALSE
  )
  summary_table <- cbind(summary_table, as.data.frame(mean_expression, check.names = FALSE))

  # Export the unscaled subtype summary as the numerical source for the panel.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    summary_table,
    file.path(output_dir, "figure2B_cell_subpopulation_summary.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  # Construct row labels, scaled marker values, and aligned annotations for the heatmap.
  rownames(mean_expression) <- paste0(subtypes, " (", format(round(summary_table$cell_fraction * 100, 1), nsmall = 1), "%)")
  scaled_expression <- scale(mean_expression)
  status_markers <- setdiff(markers, annotation_markers)
  color_function <- circlize::colorRamp2(c(-2, 0, 2), c("#2166ac", "#f7f7f7", "#b2182b"))
  cell_count_colors <- unname(major_palette[major_types])

  # Assemble abundance and major-lineage annotations beside the two marker blocks.
  cell_count_annotation <- rowAnnotation(
    CellCount = anno_barplot(
      log10(summary_table$cell_count + 1), bar_width = 0.8,
      gp = gpar(fill = cell_count_colors, col = "white"),
      axis_param = list(side = "bottom", at = c(0, 2, 4, 6), labels = c("1", "100", "10K", "1M"), gp = gpar(fontsize = 8)),
      width = unit(3, "cm")
    ),
    width = unit(3, "cm")
  )
  major_annotation <- rowAnnotation(
    MajorType = major_types,
    col = list(MajorType = major_palette),
    width = unit(0.5, "cm")
  )
  annotation_heatmap <- Heatmap(
    scaled_expression[, annotation_markers, drop = FALSE], name = "Z-score", col = color_function,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
    column_title = "Annotation Markers", width = unit(10, "cm"), rect_gp = gpar(col = "white", lwd = 0.5)
  )
  status_heatmap <- Heatmap(
    scaled_expression[, status_markers, drop = FALSE], name = "Z-score2", col = color_function,
    cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
    column_title = "Status-Associated Markers", width = unit(8, "cm"), rect_gp = gpar(col = "white", lwd = 0.5)
  )

  # Render the combined annotations and heatmaps to the Figure 2B PDF.
  pdf(file.path(output_dir, "figure2B_cell_subpopulation_marker_heatmap.pdf"), width = 16, height = 6)
  draw(
    cell_count_annotation + major_annotation + annotation_heatmap + status_heatmap,
    column_title = "Cell Type Marker Expression Profile", heatmap_legend_side = "right"
  )
  dev.off()
}

if (sys.nframe() == 0L) {
  run_figure2b(sce_path, output_dir)
}

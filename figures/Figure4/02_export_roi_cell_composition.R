#!/usr/bin/env Rscript

# Figure 4 public ROI-level cell-composition export.
#
# Purpose: Aggregate the public FDZS-1 IMC cell annotations to one wide row per
# patient, ROI, and tissue region for distribution through Zenodo.
# Inputs: The public FDZS-1 SingleCellExperiment RDS selected in EXPORT_CONFIG.
# Output: A gzip-compressed TSV containing ROI identifiers, total cell counts,
# and SubType fractions.
# Runtime: Host Conda environment Spatial with R 4.2.2.

RUN_ROI_EXPORT <- FALSE

EXPORT_CONFIG <- list(
  sce_input = "path/to/FDZS1_IMC_processed.rds",
  output = "path/to/FDZS1_ROI_cell_composition.tsv.gz",
  expected_rois = 311L
)

required_packages <- c("SingleCellExperiment", "SummarizedExperiment")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

#' Read the public FDZS-1 IMC object and select cell-level annotation fields.
#'
#' @param path Path to the public SingleCellExperiment RDS.
#' @return A data frame containing patient, ROI, tissue, and cell-subtype fields.
read_cell_annotations <- function(path) {
  # The export reads colData only; expression assays are not transformed.
  if (!file.exists(path)) stop("Input does not exist: ", path)
  sce <- readRDS(path)
  required <- c("PID", "ID", "Tissue", "SubType")
  missing <- setdiff(required, names(SummarizedExperiment::colData(sce)))
  if (length(missing) > 0L) {
    stop("FDZS-1 object is missing fields: ", paste(missing, collapse = ", "))
  }
  annotations <- as.data.frame(
    SummarizedExperiment::colData(sce)[, required],
    stringsAsFactors = FALSE
  )
  names(annotations) <- c("patient_id", "roi_id", "tissue", "cell_subtype")
  if (anyNA(annotations) || any(!nzchar(as.matrix(annotations)))) {
    stop("ROI composition fields must be complete and non-empty.")
  }
  annotations
}

#' Aggregate cell annotations to a wide ROI-level fraction matrix.
#'
#' @param annotations Cell-level annotation data frame from read_cell_annotations().
#' @param expected_rois Expected number of unique patient-ROI-tissue groups.
#' @return One row per ROI with total_cells and one fraction column per SubType.
build_roi_matrix <- function(annotations, expected_rois) {
  # A composite key preserves the nesting of ROI identifiers within patients.
  group_key <- interaction(
    annotations$patient_id,
    annotations$roi_id,
    annotations$tissue,
    drop = TRUE,
    lex.order = TRUE
  )
  subtype <- factor(
    annotations$cell_subtype,
    levels = sort(unique(annotations$cell_subtype))
  )
  counts <- unclass(table(group_key, subtype))
  total_cells <- rowSums(counts)
  fractions <- sweep(counts, 1L, total_cells, "/")

  # Recover the explicit group fields in the same order as the table rows.
  group_rows <- unique(data.frame(
    group_key = as.character(group_key),
    patient_id = annotations$patient_id,
    roi_id = annotations$roi_id,
    tissue = annotations$tissue,
    stringsAsFactors = FALSE
  ))
  group_rows <- group_rows[match(rownames(counts), group_rows$group_key), , drop = FALSE]
  if (anyNA(group_rows$group_key) || nrow(group_rows) != expected_rois) {
    stop("Expected ", expected_rois, " ROI groups; found ", nrow(group_rows), ".")
  }
  if (any(abs(rowSums(fractions) - 1) > 1e-12)) {
    stop("ROI cell fractions must sum to one.")
  }

  colnames(fractions) <- paste0("cell_fraction__", colnames(fractions))
  result <- cbind(
    group_rows[, c("patient_id", "roi_id", "tissue"), drop = FALSE],
    total_cells = as.integer(total_cells),
    as.data.frame(fractions, check.names = FALSE)
  )
  result[order(result$patient_id, result$roi_id, result$tissue), , drop = FALSE]
}

#' Write the public ROI-level matrix as a compressed TSV.
#'
#' @param data ROI-level cell-composition data frame.
#' @param path Destination .tsv.gz path.
#' @return An invisible NULL after the table is written.
write_roi_matrix <- function(data, path) {
  # Use gzip compression while retaining a directly readable tabular format.
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.table(data, connection, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(NULL)
}

#' Run the public FDZS-1 ROI composition export.
#'
#' @param config List containing sce_input, output, and expected_rois.
#' @return An invisible ROI-level data frame after writing the configured output.
run_roi_export <- function(config) {
  # Execute the read, patient-ROI aggregation, and compressed export in order.
  annotations <- read_cell_annotations(config$sce_input)
  result <- build_roi_matrix(annotations, config$expected_rois)
  write_roi_matrix(result, config$output)
  message(config$output)
  invisible(result)
}

if (isTRUE(RUN_ROI_EXPORT)) {
  run_roi_export(EXPORT_CONFIG)
}

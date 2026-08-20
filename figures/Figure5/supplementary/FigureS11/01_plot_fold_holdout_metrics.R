#!/usr/bin/env Rscript

# Target/purpose: plot Supplementary Figure 11A fold-validation and common-
# holdout metrics. Inputs: two five-row CSVs selected in CONFIG. Outputs: one
# PDF plus plotted-data and summary TSVs in output_dir. Ordered workflow:
# validate both tables, construct long-form and summary data, render the panel,
# then write the TSVs and PDF.

suppressPackageStartupMessages(library(ggplot2))

# Edit these paths before running this interactive supplementary-figure entry point.
CONFIG <- list(
  validation = "/path/to/fold_validation_metrics.csv",
  holdout = "/path/to/common_holdout_metrics.csv",
  output_dir = "/path/to/FigureS11"
)

#' Purpose: Validate one five-fold metric table against the Figure S11A contract.
#'
#' @param frame Data frame loaded from a fold-metric CSV.
#' @param path Character scalar identifying the source file in error messages.
#' @return Data frame containing the required metric columns in fold order.
require_metrics <- function(frame, path) {
  # Select the fixed fold and metric schema before combining evaluation sources.
  required <- c("fold", "accuracy", "precision_macro", "recall_macro", "f1_macro")
  missing <- setdiff(required, names(frame))
  if (length(missing) > 0) {
    stop(paste(path, "is missing columns:", paste(missing, collapse = ", ")))
  }
  frame <- frame[, required]
  if (nrow(frame) != 5 || anyDuplicated(frame$fold) || !setequal(frame$fold, 1:5)) {
    stop(paste(path, "must contain one row for each fold 1 through 5."))
  }
  numeric_values <- as.matrix(frame[, setdiff(required, "fold")])
  if (any(!is.finite(numeric_values))) {
    stop(paste(path, "contains non-finite metric values."))
  }
  frame
}

#' Purpose: Convert validation and holdout metrics to long format for plotting.
#'
#' @param validation Validated five-fold metric data frame for fold validation.
#' @param holdout Validated five-fold metric data frame for common-holdout evaluation.
#' @return Data frame with one row per evaluation, fold, and displayed metric.
build_plot_data <- function(validation, holdout) {
  # Preserve evaluation identity while mapping wide metric columns to plotted rows.
  validation$evaluation <- "Fold validation"
  holdout$evaluation <- "Common holdout"
  results <- rbind(validation, holdout)
  metric_columns <- c("accuracy", "precision_macro", "recall_macro", "f1_macro")
  metric_labels <- c("Accuracy", "Macro-precision", "Macro-recall", "Macro-F1")
  rows <- vector("list", length(metric_columns))
  for (index in seq_along(metric_columns)) {
    rows[[index]] <- data.frame(
      metric = metric_labels[[index]],
      evaluation = results$evaluation,
      fold = results$fold,
      value = results[[metric_columns[[index]]]],
      stringsAsFactors = FALSE
    )
  }
  plot_data <- do.call(rbind, rows)
  plot_data$metric <- factor(plot_data$metric, levels = metric_labels)
  plot_data$evaluation <- factor(
    plot_data$evaluation,
    levels = c("Fold validation", "Common holdout")
  )
  plot_data
}

#' Purpose: Summarize each metric and evaluation combination by mean and SD.
#'
#' @param plot_data Long-form data frame returned by `build_plot_data()`.
#' @return Data frame with one summary row per displayed metric and evaluation.
build_summary_data <- function(plot_data) {
  # Summaries retain one row for every displayed metric-by-evaluation combination.
  keys <- unique(plot_data[c("metric", "evaluation")])
  rows <- vector("list", nrow(keys))
  for (index in seq_len(nrow(keys))) {
    values <- plot_data$value[
      plot_data$metric == keys$metric[[index]] &
        plot_data$evaluation == keys$evaluation[[index]]
    ]
    rows[[index]] <- data.frame(
      metric = keys$metric[[index]],
      evaluation = keys$evaluation[[index]],
      mean = mean(values),
      sd = sd(values)
    )
  }
  summary_data <- do.call(rbind, rows)
  summary_data$metric <- factor(summary_data$metric, levels = levels(plot_data$metric))
  summary_data$evaluation <- factor(
    summary_data$evaluation,
    levels = levels(plot_data$evaluation)
  )
  summary_data
}

#' Purpose: Construct the Supplementary Figure 11A metric panel.
#'
#' @param plot_data Long-form data frame containing fold-level metric values.
#' @param summary_data Data frame containing metric means and standard deviations.
#' @return A ggplot object ready for PDF export.
build_panel <- function(plot_data, summary_data) {
  # Overlay fold observations on mean-SD bars using a shared evaluation encoding.
  colors <- c("Fold validation" = "#76AACB", "Common holdout" = "#F5B079")
  bar_dodge <- position_dodge(width = 0.72)
  point_dodge <- position_jitterdodge(
    jitter.width = 0.035, jitter.height = 0, dodge.width = 0.72, seed = 24
  )
  ggplot(summary_data, aes(x = metric, y = mean, fill = evaluation)) +
    geom_col(position = bar_dodge, width = 0.64, color = "#333333", linewidth = 0.35) +
    geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = pmin(1, mean + sd)), position = bar_dodge, width = 0.14, linewidth = 0.55, color = "#333333") +
    geom_point(data = plot_data, aes(y = value), position = point_dodge, size = 2.5, stroke = 0.55, shape = 21, color = "#252525", show.legend = FALSE) +
    scale_fill_manual(values = colors, name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), labels = sprintf("%.1f", seq(0, 1, by = 0.2)), expand = expansion(mult = c(0, 0.015))) +
    labs(x = NULL, y = "Performance") +
    theme_classic(base_size = 10, base_family = "Helvetica") +
    theme(axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 8)), axis.text.x = element_text(size = 9, color = "#222222", margin = margin(t = 5)), axis.text.y = element_text(size = 9, color = "#222222"), axis.ticks = element_line(linewidth = 0.4, color = "#333333"), axis.line = element_line(linewidth = 0.55, color = "#333333"), panel.grid.major.y = element_line(linewidth = 0.3, color = "#E5E5E5"), panel.grid.minor = element_blank(), legend.position = "right", legend.text = element_text(size = 9), plot.margin = margin(8, 10, 8, 8))
}

#' Purpose: Execute the Supplementary Figure 11A workflow.
#'
#' @param config List containing `validation`, `holdout`, and `output_dir` paths.
#' @return Invisibly returns `NULL` after writing two TSV files and one PDF.
run <- function(config) {
  # Orchestrate metric preparation, panel construction, and matched data exports.
  validation <- require_metrics(read.csv(config$validation, check.names = FALSE), config$validation)
  holdout <- require_metrics(read.csv(config$holdout, check.names = FALSE), config$holdout)
  plot_data <- build_plot_data(validation, holdout)
  summary_data <- build_summary_data(plot_data)
  panel <- build_panel(plot_data, summary_data)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(plot_data, file.path(config$output_dir, "FigureS11A_plotted_data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(summary_data, file.path(config$output_dir, "FigureS11A_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  ggsave(file.path(config$output_dir, "FigureS11A_fold_holdout_metrics.pdf"), panel, device = cairo_pdf, width = 7.2, height = 4.3, units = "in")
  invisible(NULL)
}

if (sys.nframe() == 0) {
  run(CONFIG)
}

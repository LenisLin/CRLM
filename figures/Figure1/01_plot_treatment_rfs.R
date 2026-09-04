#!/usr/bin/env Rscript

# Figure 1B treatment-stratified recurrence-free survival curves.
#
# Input: the current Supplementary Tables workbook. Supplementary Tables 1 and
# 5 provide the FDZS-1 and FDZS-2 patient-level RFS records, respectively.
# Outputs: one two-panel Kaplan-Meier PDF and a cohort-level summary TSV.
# Runtime: host Conda environment Spatial with R 4.2.2.

# Edit these paths and set RUN_FIGURE1B to TRUE before direct execution.
RUN_FIGURE1B <- FALSE

FIGURE1B_CONFIG <- list(
  supplementary_workbook = "path/to/Supplementary Tables.xlsx",
  discovery_sheet = "Supplementary Table 1",
  test_sheet = "Supplementary Table 5",
  output_dir = "path/to/figure1_outputs"
)

required_fields <- c("patient_id", "rfs_time_months", "rfs_event", "treatment")
table1_fields <- c(
  "Patient ID",
  "Recurrence status during follow-up",
  "RFS time from completion of primary and liver metastasis resection (months)",
  "Postoperative adjuvant treatment group"
)
table5_fields <- c(
  "Patient ID",
  "RFS event status (1 = recurrence; 0 = censored)",
  "RFS time from completion of primary and liver metastasis resection (months)",
  "Treatment"
)

#' Read one patient-level supplementary table.
#'
#' @param path Path to the current Supplementary Tables workbook.
#' @param sheet Exact worksheet name.
#' @return A data frame retaining the displayed field names.
read_supplementary_table <- function(path, sheet) {
  if (!file.exists(path)) stop("Input does not exist: ", path)
  as.data.frame(
    readxl::read_excel(path, sheet = sheet, skip = 1, .name_repair = "minimal"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Require the displayed source fields for one worksheet.
#'
#' @param data Raw supplementary-table data.
#' @param fields Required displayed column names.
#' @param sheet Worksheet label used in errors.
#' @return The input data invisibly.
require_source_fields <- function(data, fields, sheet) {
  missing_fields <- setdiff(fields, names(data))
  if (length(missing_fields) > 0L) {
    stop(sheet, " is missing: ", paste(missing_fields, collapse = ", "))
  }
  invisible(data)
}

#' Normalize Supplementary Table 1 to the Figure 1B contract.
#'
#' @param data Current FDZS-1 supplementary table.
#' @return A patient-level data frame with the common analysis fields.
normalize_discovery_table <- function(data) {
  require_source_fields(data, table1_fields, "Supplementary Table 1")
  event_map <- c("No recurrence" = 0L, "Recurrence" = 1L)
  treatment_map <- c(
    "FOLFOX alone" = "Chemo",
    "FOLFOX plus targeted therapy" = "Combo"
  )
  data.frame(
    patient_id = trimws(as.character(data[["Patient ID"]])),
    rfs_time_months = data[["RFS time from completion of primary and liver metastasis resection (months)"]],
    rfs_event = unname(event_map[as.character(data[["Recurrence status during follow-up"]])]),
    treatment = unname(treatment_map[as.character(data[["Postoperative adjuvant treatment group"]])]),
    stringsAsFactors = FALSE
  )
}

#' Normalize Supplementary Table 5 to the Figure 1B contract.
#'
#' @param data Current FDZS-2 supplementary table.
#' @return A patient-level data frame with the common analysis fields.
normalize_test_table <- function(data) {
  require_source_fields(data, table5_fields, "Supplementary Table 5")
  data.frame(
    patient_id = trimws(as.character(data[["Patient ID"]])),
    rfs_time_months = data[["RFS time from completion of primary and liver metastasis resection (months)"]],
    rfs_event = data[["RFS event status (1 = recurrence; 0 = censored)"]],
    treatment = as.character(data[["Treatment"]]),
    stringsAsFactors = FALSE
  )
}

#' Validate one cohort's treatment/RFS input contract.
#'
#' @param data Patient-level public table.
#' @param cohort Cohort label used in validation errors.
#' @param expected_rows Required number of unique patient rows.
#' @return A validated data frame with analysis fields normalized for plotting.
validate_patient_table <- function(data, cohort, expected_rows) {
  missing_fields <- setdiff(required_fields, names(data))
  if (length(missing_fields) > 0L) {
    stop(cohort, " table is missing: ", paste(missing_fields, collapse = ", "))
  }
  if (nrow(data) != expected_rows) {
    stop(cohort, " table must contain exactly ", expected_rows, " rows; found ", nrow(data), ".")
  }

  data$patient_id <- as.character(data$patient_id)
  if (anyNA(data$patient_id) || any(!nzchar(data$patient_id))) {
    stop(cohort, " table has missing patient_id values.")
  }
  if (anyDuplicated(data$patient_id)) {
    stop(cohort, " table has duplicate patient_id values; no deduplication is performed.")
  }

  for (field in c("rfs_time_months", "rfs_event")) {
    data[[field]] <- suppressWarnings(as.numeric(data[[field]]))
    if (any(!is.finite(data[[field]]))) stop(cohort, " table has non-finite ", field, " values.")
  }
  if (any(data$rfs_time_months <= 0)) stop(cohort, " table has non-positive rfs_time_months values.")
  if (any(!data$rfs_event %in% c(0, 1))) stop(cohort, " table rfs_event must use 1=recurrence and 0=censored.")
  if (!all(as.character(data$treatment) %in% c("Chemo", "Combo"))) {
    stop(cohort, " table treatment must contain only Chemo or Combo.")
  }
  data$treatment <- factor(as.character(data$treatment), levels = c("Chemo", "Combo"))
  if (any(table(data$treatment) == 0L)) stop(cohort, " table must contain both treatment groups.")
  data
}

#' Calculate cohort-level counts and a nominal two-sided log-rank P value.
#'
#' @param data Validated patient-level table.
#' @param cohort Cohort label written to the output summary.
#' @return A one-row data frame.
summarize_cohort <- function(data, cohort) {
  fit <- survival::survdiff(survival::Surv(rfs_time_months, rfs_event) ~ treatment, data = data)
  data.frame(
    cohort = cohort,
    n_total = nrow(data),
    n_chemo = sum(data$treatment == "Chemo"),
    n_combo = sum(data$treatment == "Combo"),
    recurrence_events = sum(data$rfs_event),
    nominal_logrank_p_value = stats::pchisq(fit$chisq, df = length(fit$n) - 1L, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

#' Build one treatment-stratified Kaplan-Meier panel.
#'
#' @param data Validated patient-level table.
#' @param cohort_title Plot title.
#' @param p_value Nominal two-sided log-rank P value from the cohort summary.
#' @param show_y_axis Whether to show the y-axis label and tick labels.
#' @return A ggplot object.
plot_treatment_rfs <- function(data, cohort_title, p_value, show_y_axis = TRUE) {
  fit <- survival::survfit(survival::Surv(rfs_time_months, rfs_event) ~ treatment, data = data)
  p_label <- if (p_value >= 0.1) {
    sprintf("log-rank p = %.2f", p_value)
  } else if (p_value >= 0.001) {
    sprintf("log-rank p = %.3f", p_value)
  } else {
    "log-rank p < 0.001"
  }
  max_time <- max(data$rfs_time_months)

  km_plot <- survminer::ggsurvplot(
    fit = fit,
    data = data,
    pval = FALSE,
    conf.int = FALSE,
    risk.table = FALSE,
    linetype = "strata",
    surv.median.line = "hv",
    censor.shape = 3,
    censor.size = 2.2,
    size = 0.8,
    break.time.by = 10,
    break.y.by = 0.25,
    xlim = c(0, max_time),
    ylim = c(0, 1),
    axes.offset = FALSE,
    ggtheme = ggplot2::theme_classic(),
    palette = c("#0073C2FF", "#EFC000FF"),
    title = cohort_title,
    xlab = "RFS Time (months)",
    ylab = "Survival probability",
    legend = "none"
  )
  panel <- km_plot$plot +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = 0.92,
      label = p_label,
      hjust = 1.05,
      vjust = 0,
      size = 3.5,
      fontface = "bold"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 9, face = "bold"),
      axis.text = ggplot2::element_text(size = 8, colour = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.45, colour = "black"),
      plot.margin = ggplot2::margin(4, 8, 4, 6)
    )
  if (!show_y_axis) {
    panel <- panel + ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  }
  panel
}

#' Render the two-panel Figure 1B replacement.
#'
#' @param discovery Validated FDZS-1 patient-level table.
#' @param test Validated FDZS-2 patient-level table.
#' @param p_values Named vector of cohort-level nominal log-rank P values.
#' @param path Destination PDF path.
#' @return An invisible NULL after writing the PDF.
write_figure1b_pdf <- function(discovery, test, p_values, path) {
  legend_data <- data.frame(
    x = rep(c(1, 2), times = 2),
    y = rep(c(1, 2), each = 2),
    treatment = factor(rep(c("Chemo", "Combo"), each = 2), levels = c("Chemo", "Combo"))
  )
  legend_plot <- ggplot2::ggplot(
    legend_data,
    ggplot2::aes(x = x, y = y, colour = treatment, linetype = treatment, shape = treatment)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(
      name = "Treatment strategy",
      values = c(Chemo = "#0073C2FF", Combo = "#EFC000FF"),
      labels = c(Chemo = "Chemotherapy Only", Combo = "Combination Therapy")
    ) +
    ggplot2::scale_linetype_manual(
      name = "Treatment strategy",
      values = c(Chemo = "solid", Combo = "dashed"),
      labels = c(Chemo = "Chemotherapy Only", Combo = "Combination Therapy")
    ) +
    ggplot2::scale_shape_manual(
      name = "Treatment strategy",
      values = c(Chemo = 3, Combo = 3),
      labels = c(Chemo = "Chemotherapy Only", Combo = "Combination Therapy")
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = 9, face = "bold"),
      legend.text = ggplot2::element_text(size = 9),
      legend.key.width = grid::unit(1.1, "cm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    )
  legend <- cowplot::get_legend(legend_plot)
  panels <- cowplot::plot_grid(
    plot_treatment_rfs(discovery, "Discovery cohort (n=35)", p_values[["FDZS-1"]]),
    plot_treatment_rfs(test, "Validation cohort (n=95)", p_values[["FDZS-2"]], show_y_axis = FALSE),
    ncol = 2,
    align = "hv",
    axis = "tblr",
    rel_widths = c(1, 1)
  )
  header <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label("B", x = 0, hjust = 0, fontface = "bold", size = 12),
    legend,
    ncol = 2,
    rel_widths = c(0.06, 0.94)
  )
  figure <- cowplot::plot_grid(header, panels, ncol = 1, rel_heights = c(0.16, 0.84))

  grDevices::cairo_pdf(path, width = 8.0, height = 3.3, family = "sans")
  print(figure)
  grDevices::dev.off()
  invisible(NULL)
}

#' Generate the Figure 1B public treatment/RFS outputs.
#'
#' @param config Editable input and output paths.
#' @return An invisible list containing the two cohort summaries.
run_figure1b <- function(config) {
  required_packages <- c("readxl", "survival", "survminer", "ggplot2", "cowplot")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0L) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
  }

  discovery <- validate_patient_table(
    normalize_discovery_table(read_supplementary_table(config$supplementary_workbook, config$discovery_sheet)),
    "FDZS-1",
    35L
  )
  test <- validate_patient_table(
    normalize_test_table(read_supplementary_table(config$supplementary_workbook, config$test_sheet)),
    "FDZS-2",
    95L
  )
  summary <- rbind(
    summarize_cohort(discovery, "FDZS-1"),
    summarize_cohort(test, "FDZS-2")
  )
  p_values <- stats::setNames(summary$nominal_logrank_p_value, summary$cohort)

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  write_figure1b_pdf(
    discovery,
    test,
    p_values,
    file.path(config$output_dir, "Figure1B_treatment_RFS_revised.pdf")
  )
  write.table(summary, file.path(config$output_dir, "Figure1B_treatment_RFS_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(summary)
}

if (isTRUE(RUN_FIGURE1B)) {
  run_figure1b(FIGURE1B_CONFIG)
}

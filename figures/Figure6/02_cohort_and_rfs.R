#!/usr/bin/env Rscript

# Figure 6B-D clinical distribution and recurrence-free survival analysis.
#
# Purpose: Generate the Figure 6B test-cohort clinical-distribution component
# and the Figure 6C-D treatment-by-ISR Kaplan-Meier components from validated
# prepared patient tables.
#
# Inputs: A 34-row discovery RFS TSV and a 95-row test RFS TSV with the schema
# validated by validate_prepared_rfs().
# Outputs: Figure 6B-D PDF components and cutoff, distribution, group-count,
# risk-count, and pairwise log-rank TSV files in FIGURE6_CONFIG$output_dir.
# Runtime: Host Conda environment Spatial with R 4.2.2.
#
# Set RUN_FIGURE6 to TRUE after editing all paths in FIGURE6_CONFIG.

RUN_FIGURE6 <- FALSE

FIGURE6_CONFIG <- list(
  discovery_rfs = "path/to/discovery_rfs.tsv",
  test_rfs = "path/to/test_rfs.tsv",
  output_dir = "path/to/figure6_outputs"
)

required_packages <- c("survival", "survminer", "ggplot2", "gridExtra")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(gridExtra)
})

#' Read a prepared Figure 6 patient table.
#'
#' @param path Path to a tab-delimited prepared patient table.
#' @return A data frame with original column names retained.
read_prepared_table <- function(path) {
  # Read the prepared patient table without altering its published field names.
  if (!file.exists(path)) stop("Input does not exist: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

#' Convert required table fields to finite numeric values.
#'
#' @param data Prepared patient data frame.
#' @param fields Character vector of field names to convert.
#' @param label Cohort label used in validation errors.
#' @return The input data frame with the requested fields stored as numeric.
check_numeric <- function(data, fields, label) {
  # Normalize analysis fields to one finite numeric representation per cohort.
  for (field in fields) {
    value <- suppressWarnings(as.numeric(data[[field]]))
    if (any(!is.finite(value))) stop(label, " has non-finite ", field, " values.")
    data[[field]] <- value
  }
  data
}

#' Validate a prepared discovery or test RFS table.
#'
#' @param data Prepared patient data frame.
#' @param cohort Either "discovery" or "test".
#' @return A validated data frame with numeric analysis fields and ordered
#' treatment levels.
validate_prepared_rfs <- function(data, cohort) {
  # Map each cohort role to its required patient-level analysis schema.
  common_fields <- c("patient_id", "wsi_id", "isr", "rfs_time_months", "rfs_event", "treatment")
  required_fields <- if (identical(cohort, "test")) {
    c(common_fields, "TBS", "CRLM_number", "CRLM_size", "fong_score", "KRAS_mutation")
  } else {
    common_fields
  }
  missing_fields <- setdiff(required_fields, names(data))
  if (length(missing_fields) > 0L) {
    stop(cohort, " table is missing: ", paste(missing_fields, collapse = ", "))
  }
  expected_rows <- if (identical(cohort, "discovery")) 34L else 95L
  if (nrow(data) != expected_rows) {
    stop(cohort, " table must contain exactly ", expected_rows, " rows; found ", nrow(data), ".")
  }

  # Patient identifiers remain unique rows and treatment uses a fixed display order.
  data$patient_id <- as.character(data$patient_id)
  data$wsi_id <- as.character(data$wsi_id)
  if (anyNA(data$patient_id) || any(!nzchar(data$patient_id))) stop(cohort, " table has missing patient_id values.")
  if (anyDuplicated(data$patient_id)) stop(cohort, " table has duplicate patient_id values; no deduplication is performed.")
  if (anyNA(data$wsi_id) || any(!nzchar(data$wsi_id))) stop(cohort, " table has missing wsi_id values.")

  data <- check_numeric(data, c("isr", "rfs_time_months", "rfs_event"), cohort)
  if (any(data$rfs_time_months <= 0)) stop(cohort, " table has non-positive rfs_time_months values.")
  if (any(!data$rfs_event %in% c(0, 1))) stop(cohort, " table rfs_event must use 0 or 1.")
  if (!all(as.character(data$treatment) %in% c("Chemo", "Combo"))) {
    stop(cohort, " table treatment must contain only Chemo or Combo.")
  }
  data$treatment <- factor(as.character(data$treatment), levels = c("Combo", "Chemo"))

  if (identical(cohort, "test")) {
    data <- check_numeric(data, c("TBS", "CRLM_number", "CRLM_size", "fong_score", "KRAS_mutation"), cohort)
  }
  data
}

#' Write a data frame as a tab-delimited table.
#'
#' @param data Data frame to write.
#' @param path Destination TSV path.
#' @return An invisible NULL value after writing the table.
write_tsv <- function(data, path) {
  # Use the same tabular representation for all Figure 6 numerical exports.
  write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE)
}

#' Assign ISR and treatment-by-ISR groups.
#'
#' @param data Validated prepared patient data frame.
#' @param cutoff Discovery-cohort median ISR cutoff.
#' @return The input data frame with isr_group and combined_group columns.
prepare_groups <- function(data, cutoff) {
  # Apply the discovery-derived cutoff identically to both cohort tables.
  data$isr_group <- ifelse(data$isr > cutoff, "High", "Low")
  data$combined_group <- factor(
    paste(data$treatment, data$isr_group, sep = "_"),
    levels = c("Combo_High", "Combo_Low", "Chemo_High", "Chemo_Low")
  )
  if (any(table(data$combined_group) == 0L)) {
    stop("All four treatment-by-ISR groups must be represented in the prepared table.")
  }
  data
}

#' Calculate all unordered pairwise log-rank comparisons.
#'
#' @param data Prepared patient data with combined_group assignments.
#' @param cohort Cohort label written to the result table.
#' @return A data frame containing six group-pair comparison results.
pairwise_logrank <- function(data, cohort) {
  # Evaluate every unordered pair on its own two-group survival subset.
  pairs <- combn(levels(data$combined_group), 2L, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    # Retain patient and event counts beside each nominal log-rank result.
    subset <- data[data$combined_group %in% pair, , drop = FALSE]
    subset$combined_group <- droplevels(subset$combined_group)
    fit <- survdiff(Surv(rfs_time_months, rfs_event) ~ combined_group, data = subset)
    p_value <- stats::pchisq(fit$chisq, df = length(fit$n) - 1L, lower.tail = FALSE)
    data.frame(
      cohort = cohort,
      group_1 = pair[[1L]],
      group_2 = pair[[2L]],
      n_group_1 = sum(data$combined_group == pair[[1L]]),
      events_group_1 = sum(data$rfs_event[data$combined_group == pair[[1L]]]),
      n_group_2 = sum(data$combined_group == pair[[2L]]),
      events_group_2 = sum(data$rfs_event[data$combined_group == pair[[2L]]]),
      chisq = unname(fit$chisq),
      nominal_p_value = p_value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

displayed_comparison_keys <- c(
  "Combo_High::Chemo_High",
  "Combo_High::Combo_Low",
  "Chemo_High::Chemo_Low",
  "Combo_Low::Chemo_Low"
)

#' Count patients and RFS events in each treatment-by-ISR group.
#'
#' @param data Prepared patient data with combined_group assignments.
#' @param cohort Cohort label written to the result table.
#' @param cutoff Discovery-cohort median ISR cutoff.
#' @return A four-row data frame of group counts, event counts, and cutoff.
group_counts <- function(data, cohort, cutoff) {
  # Keep counts aligned to the fixed treatment-by-ISR factor levels.
  groups <- levels(data$combined_group)
  data.frame(
    cohort = cohort,
    group = groups,
    n = as.integer(table(data$combined_group)[groups]),
    events = vapply(groups, function(group) sum(data$rfs_event[data$combined_group == group]), numeric(1L)),
    isr_cutoff = cutoff,
    stringsAsFactors = FALSE
  )
}

#' Calculate six-month interval numbers at risk by group.
#'
#' @param data Prepared patient data with combined_group assignments.
#' @param cohort Cohort label written to the result table.
#' @return A data frame of group-specific risk counts at six-month intervals.
risk_counts <- function(data, cohort) {
  # Use a common six-month grid across every group within a cohort.
  times <- seq.int(0L, ceiling(max(data$rfs_time_months) / 6) * 6L, by = 6L)
  do.call(rbind, lapply(levels(data$combined_group), function(group) {
    # Count patients whose observed time reaches each displayed risk-table point.
    subset <- data[data$combined_group == group, , drop = FALSE]
    data.frame(
      cohort = cohort,
      group = group,
      time_months = times,
      n_at_risk = vapply(times, function(time) sum(subset$rfs_time_months >= time), numeric(1L)),
      stringsAsFactors = FALSE
    )
  }))
}

#' Render a treatment-by-ISR Kaplan-Meier PDF component.
#'
#' @param data Prepared patient data with combined_group assignments.
#' @param pairwise_results Four ordered pairwise comparison rows for annotation.
#' @param path Destination PDF path.
#' @return The graphics device number returned after writing the PDF.
plot_km <- function(data, pairwise_results, path) {
  # Fit patient-level Kaplan-Meier curves for the four fixed combined groups.
  fit <- survfit(Surv(rfs_time_months, rfs_event) ~ combined_group, data = data)
  # Prepare the four prespecified pairwise results as a compact panel annotation.
  annotation <- paste(
    paste(
      gsub("_", " ", pairwise_results$group_1), "vs",
      gsub("_", " ", pairwise_results$group_2), ":",
      ifelse(pairwise_results$nominal_p_value < 0.001, "p < 0.001", sprintf("p = %.3f", pairwise_results$nominal_p_value))
    ),
    collapse = "\n"
  )
  plot <- ggsurvplot(
    fit,
    data = data,
    palette = c("#1f77b4", "#aec7e8", "#d62728", "#ff9896"),
    size = 1.1,
    conf.int = TRUE,
    conf.int.style = "ribbon",
    conf.int.alpha = 0.2,
    pval = TRUE,
    pval.size = 6,
    pval.coord = c(0, 0.1),
    xlab = "Time (Months)",
    ylab = "Recurrence-Free Survival",
    title = "Treatment and Risk Score Stratified Survival Analysis",
    break.time.by = 6,
    xlim = c(0, max(data$rfs_time_months) + 2),
    legend.title = "Patient Group",
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    risk.table.title = "Number at Risk",
    tables.theme = theme_survminer(font.main = 12),
    ggtheme = theme_classic2(base_size = 14) + theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = c(0.8, 0.35),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.title = element_text(face = "bold")
    )
  )
  plot$plot <- plot$plot + annotate(
    geom = "label",
    x = max(data$rfs_time_months) * 0.6,
    y = 0.95,
    label = annotation,
    hjust = 0,
    vjust = 1,
    size = 4.5,
    fill = "white",
    label.padding = grid::unit(0.4, "lines"),
    label.r = grid::unit(0.1, "lines"),
    color = "black"
  )
  # Export the survival curve and number-at-risk table as one PDF component.
  grDevices::pdf(path, width = 12, height = 9)
  print(plot)
  grDevices::dev.off()
}

#' Define the compact publication theme for Figure 6B strips.
#'
#' @return A ggplot2 theme object.
theme_pub_strip <- function() {
  # Share compact typography and margins across all five clinical strips.
  theme_classic(base_size = 9) + theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )
}

#' Create a continuous-variable distribution strip.
#'
#' @param data Test-cohort prepared patient data.
#' @param variable Name of the numeric variable to display.
#' @return A ggplot object with histogram and density layers.
plot_continuous <- function(data, variable) {
  # Pair histogram density with a smooth density trace for one numeric field.
  ggplot(data, aes(x = .data[[variable]])) +
    geom_histogram(aes(y = after_stat(density)), fill = "#4E79A7", color = "white", alpha = 0.4, bins = 20) +
    geom_density(color = "#4E79A7", linewidth = 0.8) +
    labs(title = variable, x = NULL, y = "Density") +
    theme_pub_strip()
}

#' Create a categorical-variable distribution strip.
#'
#' @param data Test-cohort prepared patient data.
#' @param variable Name of the categorical variable to display.
#' @return A ggplot object with category counts.
plot_categorical <- function(data, variable) {
  # Display category counts directly above bars for one discrete field.
  ggplot(data, aes(x = factor(.data[[variable]]))) +
    geom_bar(fill = "#4E79A7", width = 0.7) +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 2.5) +
    labs(title = variable, x = NULL, y = "N") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme_pub_strip()
}

#' Summarize selected Figure 6B clinical variables.
#'
#' @param data Test-cohort prepared patient data.
#' @return A data frame of continuous summaries and categorical counts.
distribution_summary <- function(data) {
  # Construct one export table spanning continuous summaries and category counts.
  continuous <- c("TBS", "CRLM_number", "CRLM_size")
  categorical <- c("fong_score", "KRAS_mutation")
  continuous_rows <- do.call(rbind, lapply(continuous, function(variable) {
    value <- data[[variable]]
    data.frame(
      variable = variable,
      data_type = "continuous",
      level = NA_character_,
      statistic = c("n", "min", "q1", "median", "q3", "max", "mean", "sd"),
      value = c(length(value), min(value), unname(quantile(value, 0.25)), median(value), unname(quantile(value, 0.75)), max(value), mean(value), sd(value)),
      stringsAsFactors = FALSE
    )
  }))
  categorical_rows <- do.call(rbind, lapply(categorical, function(variable) {
    counts <- table(data[[variable]], useNA = "no")
    data.frame(
      variable = variable,
      data_type = "categorical",
      level = names(counts),
      statistic = "count",
      value = as.numeric(counts),
      stringsAsFactors = FALSE
    )
  }))
  rbind(continuous_rows, categorical_rows)
}

#' Render the five-panel Figure 6B clinical-distribution PDF component.
#'
#' @param data Test-cohort prepared patient data.
#' @param path Destination PDF path.
#' @return The graphics device number returned after writing the PDF.
plot_figure6b <- function(data, path) {
  # Assemble the five single-variable plots in the manuscript panel order.
  plots <- list(
    plot_continuous(data, "TBS"),
    plot_continuous(data, "CRLM_number"),
    plot_continuous(data, "CRLM_size"),
    plot_categorical(data, "fong_score"),
    plot_categorical(data, "KRAS_mutation")
  )
  grDevices::pdf(path, width = 10, height = 3)
  gridExtra::grid.arrange(grobs = plots, ncol = 5)
  grDevices::dev.off()
}

#' Generate all Figure 6B-D output components.
#'
#' @param discovery_path Path to the prepared 34-row discovery RFS TSV.
#' @param test_path Path to the prepared 95-row test RFS TSV.
#' @param output_dir Directory for all Figure 6B-D output files.
#' @return An invisible NULL value after writing the Figure 6B-D components.
run_figure6bcd <- function(discovery_path, test_path, output_dir) {
  # Load both patient cohorts before deriving the shared discovery ISR cutoff.
  discovery <- validate_prepared_rfs(read_prepared_table(discovery_path), "discovery")
  test <- validate_prepared_rfs(read_prepared_table(test_path), "test")
  cutoff <- stats::median(discovery$isr)
  if (!is.finite(cutoff)) stop("Discovery ISR median is not finite.")

  # Assign comparable groups, then create the common output location.
  discovery <- prepare_groups(discovery, cutoff)
  test <- prepare_groups(test, cutoff)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Export the cutoff and Figure 6B distribution inputs before survival outputs.
  write_tsv(data.frame(
    cutoff_name = "discovery_median_isr",
    cutoff_value = cutoff,
    discovery_rows = nrow(discovery),
    test_rows = nrow(test),
    high_definition = "isr > discovery_median_isr",
    low_definition = "isr <= discovery_median_isr",
    stringsAsFactors = FALSE
  ), file.path(output_dir, "Figure6_isr_cutoff.tsv"))
  write_tsv(distribution_summary(test), file.path(output_dir, "Figure6B_distribution_summary.tsv"))
  plot_figure6b(test, file.path(output_dir, "Figure6B_test_clinical_distributions.pdf"))

  # Generate matched statistical tables and Kaplan-Meier components for each cohort.
  for (cohort in c("discovery", "test")) {
    data <- if (identical(cohort, "discovery")) discovery else test
    panel <- if (identical(cohort, "discovery")) "Figure6C" else "Figure6D"
    comparisons <- pairwise_logrank(data, cohort)
    write_tsv(comparisons, file.path(output_dir, paste0(panel, "_pairwise_logrank.tsv")))
    comparison_keys <- paste(comparisons$group_1, comparisons$group_2, sep = "::")
    displayed_comparisons <- comparisons[match(displayed_comparison_keys, comparison_keys), , drop = FALSE]
    write_tsv(group_counts(data, cohort, cutoff), file.path(output_dir, paste0(panel, "_group_counts.tsv")))
    write_tsv(risk_counts(data, cohort), file.path(output_dir, paste0(panel, "_risk_counts.tsv")))
    plot_km(data, displayed_comparisons, file.path(output_dir, paste0(panel, if (identical(cohort, "discovery")) "_discovery_rfs.pdf" else "_test_rfs.pdf")))
  }
}

if (sys.nframe() == 0L && isTRUE(RUN_FIGURE6)) {
  run_figure6bcd(
    FIGURE6_CONFIG$discovery_rfs,
    FIGURE6_CONFIG$test_rfs,
    FIGURE6_CONFIG$output_dir
  )
}

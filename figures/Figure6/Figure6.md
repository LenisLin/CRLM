Figure 6. SpMap-derived ISR stratifies recurrence risk after adjuvant
therapy in CRLM.
================

## Purpose

Figure 6 examines the association between H&E-derived immune-to-stromal
ratio (ISR), postoperative adjuvant treatment, and recurrence-free
survival (RFS) in colorectal liver metastases (CRLM). The report
displays the complete publication Figure 6 and its survival supplements
alongside the clinical-distribution, Kaplan-Meier, continuous-ISR, and
binary-ISR code. Supplementary Figures 13 and 14 provide
overall-survival (OS) and pathology-model comparisons.

## Inputs and analytical units

The patient is the statistical unit. The [SpMap
workflow](../Figure5/README.md) supplies patient-linked H&E scores,
calculated as the number of PIR tiles divided by the number of PSM
tiles. The [Python export entry point](03_prepare_public_isr_tables.py)
joins saved scores to the clinical workbook, retains the 34 FDZS-1
patients with available H&E ISR, assigns public FDZS-2 identifiers, and
exports patient ISR/RFS, cohort-assignment, and score-definition tables.
Its IMC-score field uses `clip((N_PIR + 1)/(N_PSM + 1), 0.1, 10)`.

The [clinical and RFS entry point](02_cohort_and_rfs.R) requires one
complete row per patient: 34 discovery rows and 95 independent-test
rows. Both prepared TSVs contain `patient_id`, `wsi_id`, `isr`,
`rfs_time_months`, `rfs_event`, and `treatment`; the test table also
contains `TBS`, `CRLM_number`, `CRLM_size`, `fong_score`, and
`KRAS_mutation`. A WSI identifier may be shared across patient rows. RFS
time is measured in months from completion of primary and liver
metastasis resection; `rfs_event` is coded `1 = recurrence` and
`0 = censored`. Treatment is `Chemo` (FOLFOX alone) or `Combo` (FOLFOX
plus targeted therapy).

The [adjusted continuous-ISR entry
point](04_continuous_isr_multivariable_cox.R) reads the Supplementary
Tables workbook with a title row followed by column headers.
Supplementary Table 1 supplies 35 discovery clinical rows, joined to the
34 rows with `figure6_included == 1` in `FDZS1_patient_level_ISR.tsv`;
Supplementary Table 5 supplies 95 test rows. Required adjustment fields
are age, gender, KRAS status, Fong score, CEA, CA19-9, TBS, CRLM number
and largest diameter, differentiation grade, pathological T stage, and
lymph-node status. Source column labels and category mappings are given
in the code below. The [data contract](../../docs/data_contract.md)
describes the released tables and the supplementary-workbook
requirements.

## Figure 6B-D: clinical distributions and treatment-stratified RFS

### Configuration and prepared-table validation

The source defines paths and packages, reads the two prepared tables,
and checks cohort size, patient identifiers, numerical fields, event
coding, and treatment labels.

<!-- source: 02_cohort_and_rfs.R; lines: 1-37 -->

``` r
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
```

<!-- source: 02_cohort_and_rfs.R; lines: 38-117 -->

``` r
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
```

### ISR groups and survival comparisons

High ISR is defined by a score strictly above the discovery median,
0.581262877; the same cutoff is applied to the test cohort. Values at or
below the cutoff are Low ISR. Kaplan-Meier estimates compare four
treatment-by-ISR strata. The code calculates all six pairwise nominal
two-sided log-rank tests and selects four comparisons for panel
annotations; it also exports group/event counts and six-month numbers at
risk.

<!-- source: 02_cohort_and_rfs.R; lines: 118-172 -->

``` r
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
```

<!-- source: 02_cohort_and_rfs.R; lines: 173-212 -->

``` r
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
```

### Kaplan-Meier display

Curves include confidence bands, censoring marks, an overall log-rank P
value, four pairwise annotations, and a number-at-risk table.

<!-- source: 02_cohort_and_rfs.R; lines: 213-278 -->

``` r
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
```

The publication curves are shown in the complete Figure 6 below.
High-ISR treatment comparisons have nominal log-rank P = 0.036 in
discovery and P = 0.023 in test. The High-ISR Combo discovery stratum
contains six patients and zero recurrence events. The saved analysis
PDFs are available for
[discovery](report_assets/figure6c_discovery_rfs.pdf) and
[test](report_assets/figure6d_test_rfs.pdf).

### Test-cohort clinical distributions

The source displays TBS, CRLM number, and largest CRLM diameter with
20-bin density-scaled histograms and kernel-density curves. Fong score
and KRAS status are displayed as category counts. The summary export
contains continuous-variable summaries and categorical counts.

<!-- source: 02_cohort_and_rfs.R; lines: 279-321 -->

``` r
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
```

<!-- source: 02_cohort_and_rfs.R; lines: 322-373 -->

``` r
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
```

Panel B of the complete publication sheet below contains the three
continuous distributions and the Fong-score and KRAS bars. The saved
source sheets are available as [continuous
distributions](report_assets/figure6b_continuous.pdf) and [categorical
distributions](report_assets/figure6b_categorical.pdf).

### Output assembly

The driver derives the discovery cutoff, applies it to both cohorts, and
writes the distribution, Kaplan-Meier, group-count, risk-count, and
log-rank outputs.

<!-- source: 02_cohort_and_rfs.R; lines: 374-425 -->

``` r
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
```

## Published Figure 6

![Complete publication Figure 6: cohort flow, clinical distributions,
RFS, and Cox forests](report_assets/figure6_publication.png)

**Figure 6. SpMap-derived ISR stratifies recurrence risk after adjuvant
therapy in CRLM.** (A) FDZS-2 patient inclusion and treatment groups.
(B) Test-cohort TBS, metastases count, largest metastasis size, Fong
score, and KRAS status. (C-D) RFS for 34 discovery and 95 test patients
in four ISR-by-treatment groups; red and blue denote Combo and Chemo,
with darker shades indicating High ISR. Shading shows 95% confidence
intervals, crosses indicate censoring, and risk tables give patients at
risk. (E-F) Cox forest plots for the discovery and test cohorts. Forest
coordinates are log10(HR), with HR = 1 as reference; numerical columns
report HRs, 95% confidence intervals, and nominal P values. [Download
the original publication PDF](report_assets/figure6_publication.pdf).

## Adjusted continuous-ISR models

### Model specification and input normalization

The same Cox model is specified separately for the discovery and test
cohorts: continuous ISR, treatment, their interaction, and 12
clinicopathological covariates. Chemo is the treatment reference. Age,
TBS, CRLM number, largest CRLM diameter, CEA, and CA19-9 are
standardized within cohort; Fong score and T stage remain numeric on
their original scales. Gender uses the cohort-specific factor levels
shown in the source. Each fit contains 15 coefficients.

<!-- source: 04_continuous_isr_multivariable_cox.R; lines: 1-75 -->

``` r
#!/usr/bin/env Rscript

# Figure 6E-F cohort-specific multivariable Cox models with continuous ISR.
#
# Inputs: current Supplementary Tables 1 and 5 plus the FDZS-1 patient-level
# H&E-derived ISR table. Outputs: one coefficient table and two forest-plot PDFs.
# Runtime: host Conda environment Spatial with R 4.2.2.

# Edit these paths and set RUN_FIGURE6_COX to TRUE before direct execution.
RUN_FIGURE6_COX <- FALSE

FIGURE6_COX_CONFIG <- list(
  supplementary_workbook = "path/to/Supplementary Tables.xlsx",
  discovery_isr = "path/to/FDZS1_patient_level_ISR.tsv",
  discovery_sheet = "Supplementary Table 1",
  test_sheet = "Supplementary Table 5",
  output_dir = "path/to/figure6_outputs"
)

discovery_fields <- c(
  "Patient ID",
  "Recurrence status during follow-up",
  "RFS time from completion of primary and liver metastasis resection (months)",
  "Postoperative adjuvant treatment group",
  "Age (years)",
  "Gender",
  "Fong clinical risk score",
  "KRAS mutation status",
  "Tumor burden score (TBS)",
  "Number of CRLM",
  "Largest CRLM diameter (cm)",
  "Preoperative CEA (ng/mL)",
  "Preoperative CA19-9 (U/mL)",
  "Differentiation grade",
  "Pathological T stage",
  "Lymph node metastasis status"
)

test_fields <- c(
  "Patient ID",
  "RFS event status (1 = recurrence; 0 = censored)",
  "RFS time from completion of primary and liver metastasis resection (months)",
  "Treatment",
  "ISR",
  "Age",
  "Gender",
  "Fong score",
  "KRAS mutation",
  "TBS",
  "CRLM number",
  "CRLM_size",
  "CEA",
  "CA199",
  "Differential_grade",
  "T_stage",
  "Lymph_positive"
)

model_fields <- c(
  "patient_id", "rfs_time_months", "rfs_event", "isr", "treatment",
  "age", "gender", "kras", "fong", "cea", "ca199", "tbs",
  "crlm_number", "crlm_size", "differentiation_grade", "t_stage",
  "lymph_positive"
)

scaled_covariates <- c("age", "tbs", "crlm_number", "crlm_size", "cea", "ca199")

continuous_isr_formula <- stats::as.formula(
  paste(
    "survival::Surv(rfs_time_months, rfs_event) ~",
    "isr * treatment + age + gender + kras + fong + cea + ca199 +",
    "tbs + crlm_number + crlm_size + differentiation_grade + t_stage + lymph_positive"
  )
)
```

<!-- source: 04_continuous_isr_multivariable_cox.R; lines: 76-103 -->

``` r
#' Read one worksheet from the current supplementary workbook.
#'
#' @param path Workbook path.
#' @param sheet Exact worksheet name.
#' @return A data frame with displayed column names retained.
read_supplementary_table <- function(path, sheet) {
  if (!file.exists(path)) stop("Input does not exist: ", path)
  as.data.frame(
    readxl::read_excel(path, sheet = sheet, skip = 1, .name_repair = "minimal"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Require source columns at the workbook or TSV boundary.
#'
#' @param data Input data frame.
#' @param fields Required columns.
#' @param source Source label used in an error.
#' @return The input invisibly.
require_fields <- function(data, fields, source) {
  missing_fields <- setdiff(fields, names(data))
  if (length(missing_fields) > 0L) {
    stop(source, " is missing: ", paste(missing_fields, collapse = ", "))
  }
  invisible(data)
}
```

The discovery normalization maps displayed clinical categories to model
fields and joins H&E ISR by patient. The test normalization reads ISR
and clinical covariates directly from Supplementary Table 5. Both model
frames require complete fields, unique patient rows, positive RFS times,
and conventional event coding.

<!-- source: 04_continuous_isr_multivariable_cox.R; lines: 104-148 -->

``` r
#' Normalize the FDZS-1 clinical and H&E-derived ISR inputs.
#'
#' @param clinical Supplementary Table 1.
#' @param scores Public FDZS-1 patient-level ISR table.
#' @return The 34-patient discovery model frame.
normalize_discovery_cohort <- function(clinical, scores) {
  require_fields(clinical, discovery_fields, "Supplementary Table 1")
  require_fields(scores, c("patient_id", "figure6_included", "isr"), "FDZS-1 ISR table")

  event_map <- c("No recurrence" = 0, "Recurrence" = 1)
  treatment_map <- c(
    "FOLFOX alone" = "Chemo",
    "FOLFOX plus targeted therapy" = "Combo"
  )
  gender_map <- c("Male" = "Male", "Female" = "Female")
  kras_map <- c("Wild-type" = 0, "Mutant" = 1)
  grade_map <- c("Moderately differentiated" = 0, "Poorly differentiated" = 1)
  lymph_map <- c("Negative" = 0, "Positive" = 1)

  base <- data.frame(
    patient_id = trimws(as.character(clinical[["Patient ID"]])),
    rfs_time_months = clinical[["RFS time from completion of primary and liver metastasis resection (months)"]],
    rfs_event = unname(event_map[as.character(clinical[["Recurrence status during follow-up"]])]),
    treatment = unname(treatment_map[as.character(clinical[["Postoperative adjuvant treatment group"]])]),
    age = clinical[["Age (years)"]],
    gender = unname(gender_map[as.character(clinical[["Gender"]])]),
    kras = unname(kras_map[as.character(clinical[["KRAS mutation status"]])]),
    fong = clinical[["Fong clinical risk score"]],
    cea = clinical[["Preoperative CEA (ng/mL)"]],
    ca199 = clinical[["Preoperative CA19-9 (U/mL)"]],
    tbs = clinical[["Tumor burden score (TBS)"]],
    crlm_number = clinical[["Number of CRLM"]],
    crlm_size = clinical[["Largest CRLM diameter (cm)"]],
    differentiation_grade = unname(grade_map[as.character(clinical[["Differentiation grade"]])]),
    t_stage = sub("^T", "", as.character(clinical[["Pathological T stage"]])),
    lymph_positive = unname(lymph_map[as.character(clinical[["Lymph node metastasis status"]])]),
    stringsAsFactors = FALSE
  )
  score_rows <- scores[scores$figure6_included == 1, c("patient_id", "isr"), drop = FALSE]
  score_rows$patient_id <- trimws(as.character(score_rows$patient_id))
  result <- merge(base, score_rows, by = "patient_id", all = FALSE, sort = FALSE)
  result$gender <- factor(result$gender, levels = c("Male", "Female"))
  validate_model_frame(result, "FDZS-1", 34L)
}
```

<!-- source: 04_continuous_isr_multivariable_cox.R; lines: 149-206 -->

``` r
#' Normalize Supplementary Table 5 to the common Cox-model fields.
#'
#' @param clinical Supplementary Table 5.
#' @return The 95-patient test model frame.
normalize_test_cohort <- function(clinical) {
  require_fields(clinical, test_fields, "Supplementary Table 5")
  result <- data.frame(
    patient_id = trimws(as.character(clinical[["Patient ID"]])),
    rfs_time_months = clinical[["RFS time from completion of primary and liver metastasis resection (months)"]],
    rfs_event = clinical[["RFS event status (1 = recurrence; 0 = censored)"]],
    isr = clinical[["ISR"]],
    treatment = as.character(clinical[["Treatment"]]),
    age = clinical[["Age"]],
    gender = factor(
      as.character(clinical[["Gender"]]),
      levels = c("1", "2")
    ),
    kras = clinical[["KRAS mutation"]],
    fong = clinical[["Fong score"]],
    cea = clinical[["CEA"]],
    ca199 = clinical[["CA199"]],
    tbs = clinical[["TBS"]],
    crlm_number = clinical[["CRLM number"]],
    crlm_size = clinical[["CRLM_size"]],
    differentiation_grade = clinical[["Differential_grade"]],
    t_stage = clinical[["T_stage"]],
    lymph_positive = clinical[["Lymph_positive"]],
    stringsAsFactors = FALSE
  )
  validate_model_frame(result, "FDZS-2", 95L)
}

#' Validate and type one cohort-level Cox model frame.
#'
#' @param data Normalized patient rows.
#' @param cohort Cohort label.
#' @param expected_rows Expected patient count.
#' @return A typed model frame.
validate_model_frame <- function(data, cohort, expected_rows) {
  require_fields(data, model_fields, cohort)
  if (nrow(data) != expected_rows || anyDuplicated(data$patient_id)) {
    stop(cohort, " must contain exactly ", expected_rows, " unique patient rows.")
  }
  numeric_fields <- setdiff(model_fields, c("patient_id", "treatment", "gender"))
  for (field in numeric_fields) data[[field]] <- suppressWarnings(as.numeric(data[[field]]))
  if (any(!stats::complete.cases(data[, model_fields, drop = FALSE]))) {
    stop(cohort, " has missing or non-numeric model fields.")
  }
  if (any(data$rfs_time_months <= 0) || any(!data$rfs_event %in% c(0, 1))) {
    stop(cohort, " has invalid RFS time or event coding.")
  }
  if (!setequal(as.character(data$treatment), c("Chemo", "Combo"))) {
    stop(cohort, " treatment must contain Chemo and Combo.")
  }
  data$treatment <- factor(as.character(data$treatment), levels = c("Chemo", "Combo"))
  data
}
```

### Estimation and forest-plot output

The model uses Efron handling of tied event times and exports
coefficients, standard errors, hazard ratios (HRs), 95% confidence
intervals, nominal Wald P values, patient/event counts, formula, and
scaling metadata. The interaction HR describes the multiplicative change
in the Combo-versus-Chemo HR per one-unit increase in ISR; the treatment
main effect is the contrast at ISR = 0.

<!-- source: 04_continuous_isr_multivariable_cox.R; lines: 207-265 -->

``` r
#' Fit the fixed continuous-ISR multivariable Cox model.
#'
#' @param data Validated cohort model frame.
#' @param cohort Cohort label written to results.
#' @return A list containing the fit and tidy coefficient table.
fit_continuous_isr_cox <- function(data, cohort) {
  model_data <- data
  for (field in scaled_covariates) model_data[[field]] <- as.numeric(scale(model_data[[field]]))
  fit <- survival::coxph(
    continuous_isr_formula,
    data = model_data,
    ties = "efron",
    model = TRUE,
    x = TRUE
  )
  fit_summary <- summary(fit)
  coefficients <- fit_summary$coefficients
  intervals <- fit_summary$conf.int
  terms <- rownames(coefficients)
  labels <- c(
    isr = "Continuous ISR",
    treatmentCombo = "Treatment: Combo vs Chemo",
    age = "Age (per SD)",
    genderFemale = "Gender: Female vs Male",
    gender2 = "Gender: code 2 vs 1",
    kras = "KRAS mutation",
    fong = "Fong score",
    cea = "CEA (per SD)",
    ca199 = "CA19-9 (per SD)",
    tbs = "TBS (per SD)",
    crlm_number = "CRLM number (per SD)",
    crlm_size = "Largest CRLM size (per SD)",
    differentiation_grade = "Differentiation grade",
    t_stage = "Pathological T stage",
    lymph_positive = "Lymph-node positive",
    `isr:treatmentCombo` = "Continuous ISR x treatment"
  )
  result <- data.frame(
    cohort = cohort,
    model_term = terms,
    display_label = unname(labels[terms]),
    coefficient = coefficients[, "coef"],
    standard_error = coefficients[, "se(coef)"],
    hazard_ratio = intervals[, "exp(coef)"],
    lower_95_ci = intervals[, "lower .95"],
    upper_95_ci = intervals[, "upper .95"],
    nominal_p_value = coefficients[, "Pr(>|z|)"],
    n_patients = fit$n,
    n_events = fit$nevent,
    model_formula = paste(deparse(continuous_isr_formula), collapse = " "),
    scaled_covariates = paste(scaled_covariates, collapse = ";"),
    stringsAsFactors = FALSE
  )
  if (nrow(result) != 15L || anyNA(result$display_label)) {
    stop(cohort, " model did not return the expected 15 named coefficients.")
  }
  list(fit = fit, coefficients = result)
}
```

Forest plots place HRs and their 95% confidence intervals on a
logarithmic axis with HR = 1 as the reference. The driver writes a
combined coefficient table and one PDF for each cohort.

<!-- source: 04_continuous_isr_multivariable_cox.R; lines: 266-298 -->

``` r
#' Write one forest plot from a fitted coefficient table.
#'
#' @param results Fifteen-row coefficient table.
#' @param title Panel title.
#' @param path Output PDF path.
#' @return Invisible NULL after writing the PDF.
write_forest_plot <- function(results, title, path) {
  format_p <- function(value) if (value < 0.001) "<0.001" else sprintf("%.3f", value)
  label_text <- cbind(
    c("Model term", results$display_label),
    c("HR (95% CI)", sprintf("%.2g (%.2g-%.2g)", results$hazard_ratio, results$lower_95_ci, results$upper_95_ci)),
    c("P value", vapply(results$nominal_p_value, format_p, character(1)))
  )
  plot_object <- forestplot::forestplot(
    labeltext = label_text,
    mean = c(NA, results$hazard_ratio),
    lower = c(NA, results$lower_95_ci),
    upper = c(NA, results$upper_95_ci),
    is.summary = c(TRUE, rep(FALSE, nrow(results))),
    zero = 1,
    xlog = TRUE,
    boxsize = 0.18,
    graph.pos = 2,
    title = title,
    xlab = "Hazard ratio",
    new_page = FALSE
  )
  grDevices::pdf(path, width = 10, height = 10)
  print(plot_object)
  grDevices::dev.off()
  invisible(NULL)
}
```

<!-- source: 04_continuous_isr_multivariable_cox.R; lines: 299-344 -->

``` r
#' Generate the fixed Figure 6E-F model outputs.
#'
#' @param config Editable input and output paths.
#' @return Invisible combined coefficient table.
run_figure6_continuous_isr_cox <- function(config) {
  required_packages <- c("readxl", "survival", "forestplot")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0L) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
  }
  discovery_scores <- utils::read.delim(config$discovery_isr, check.names = FALSE)
  discovery <- normalize_discovery_cohort(
    read_supplementary_table(config$supplementary_workbook, config$discovery_sheet),
    discovery_scores
  )
  test <- normalize_test_cohort(
    read_supplementary_table(config$supplementary_workbook, config$test_sheet)
  )
  discovery_fit <- fit_continuous_isr_cox(discovery, "FDZS-1")
  test_fit <- fit_continuous_isr_cox(test, "FDZS-2")
  results <- rbind(discovery_fit$coefficients, test_fit$coefficients)

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    results,
    file.path(config$output_dir, "Figure6EF_continuous_ISR_multivariable_Cox.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write_forest_plot(
    discovery_fit$coefficients,
    "Discovery set: multivariable Cox model",
    file.path(config$output_dir, "Figure6E_continuous_ISR_multivariable_Cox.pdf")
  )
  write_forest_plot(
    test_fit$coefficients,
    "Test set: multivariable Cox model",
    file.path(config$output_dir, "Figure6F_continuous_ISR_multivariable_Cox.pdf")
  )
  invisible(results)
}

if (isTRUE(RUN_FIGURE6_COX)) {
  run_figure6_continuous_isr_cox(FIGURE6_COX_CONFIG)
}
```

## Continuous and binary ISR models for both cohorts

The [model entry point](05_continuous_and_binary_isr_cox.R) fits
continuous ISR and binary ISR models separately in both cohorts. Each
model includes treatment, the ISR-by-treatment interaction, and the same
12 clinicopathological covariates. Binary High ISR is strictly above
0.581262877; High and Chemo are the reference levels. The six continuous
clinical covariates are standardized within cohort. Coefficient tables
retain the fitted terms, and binary forest tables express the
interaction as the High/Low treatment-HR ratio.

The entry reads the current supplementary workbook and the public FDZS-1
ISR table. Set the paths and sheet names in `ISR_MODEL_CONFIG`, enable
`RUN_ISR_MODELS`, and run the entry from the repository root. The
[README](README.md) lists inputs and the four forest outputs.

<!-- source: 05_continuous_and_binary_isr_cox.R; lines: 1-37 -->

``` r
#!/usr/bin/env Rscript

# Figure 6: continuous and binary ISR interaction models in both cohorts.
# Inputs: current Supplementary Tables 1/5 and FDZS1_patient_level_ISR.tsv.
# Outputs: coefficient/contrast tables and four cohort-by-model forest PDFs.
# Run from the repository root in Spatial after editing paths and enabling the switch.
# Entry 04 supplies discovery normalization and the shared forest display function.

source("figures/Figure6/04_continuous_isr_multivariable_cox.R")

RUN_ISR_MODELS <- FALSE
ISR_MODEL_CONFIG <- list(
  supplementary_workbook = "path/to/Supplementary Tables.xlsx",
  discovery_isr = "path/to/FDZS1_patient_level_ISR.tsv",
  discovery_sheet = "Supplementary Table 1",
  test_sheet = "Supplementary Table 5",
  output_dir = "path/to/figure6_cox_outputs"
)
ISR_CUTOFF <- 0.581262877
clinical_covariates <- c(
  "age", "gender", "kras", "fong", "cea", "ca199", "tbs",
  "crlm_number", "crlm_size", "differentiation_grade", "t_stage", "lymph_positive"
)
term_labels <- c(
  isr = "Continuous ISR (per unit)",
  treatmentCombo = "Treatment: Combo vs Chemo at ISR = 0",
  `isr:treatmentCombo` = "Continuous ISR x treatment",
  isr_groupLow = "ISR group: Low vs High in Chemo",
  `isr_groupLow:treatmentCombo` = "Treatment HR ratio: Low / High",
  age = "Age (per SD)", genderFemale = "Gender: Female vs Male",
  kras = "KRAS mutation", fong = "Fong score", cea = "CEA (per SD)",
  ca199 = "CA19-9 (per SD)", tbs = "TBS (per SD)",
  crlm_number = "Metastases count (per SD)", crlm_size = "Largest metastasis size (per SD)",
  differentiation_grade = "Differentiation grade", t_stage = "Pathological T stage",
  lymph_positive = "Lymph-node positive"
)
```

<!-- source: 05_continuous_and_binary_isr_cox.R; lines: 38-65 -->

``` r
# Normalize the current reader-facing Table 5 headers and clinical categories.
# Returns the same 95-patient model-frame schema used by the continuous entry.
normalize_current_test <- function(clinical) {
  fields <- c(
    patient_id = "Patient ID", rfs_event = "Recurrence status during follow-up",
    rfs_time_months = "RFS time from completion of primary and liver metastasis resection (months)",
    treatment = "Postoperative adjuvant treatment group",
    isr = "SpMap-derived immune-stromal ratio (ISR)", age = "Age (years)",
    gender = "Gender", fong = "Fong clinical risk score", kras = "KRAS mutation status",
    cea = "Preoperative CEA (ng/mL)", ca199 = "Preoperative CA19-9 (U/mL)",
    tbs = "Tumor burden score (TBS)", crlm_number = "Number of CRLM",
    crlm_size = "Largest CRLM diameter (cm)", differentiation_grade = "Differentiation grade",
    t_stage = "Pathological T stage", lymph_positive = "Lymph node metastasis status"
  )
  require_fields(clinical, unname(fields), "Supplementary Table 5")
  # Current proof workbooks include a footnote row after the patient observations.
  data <- clinical[!is.na(clinical[[fields[["rfs_time_months"]]]]), unname(fields), drop = FALSE]
  names(data) <- names(fields)
  data$rfs_event <- unname(c("No recurrence" = 0, "Recurrence" = 1)[data$rfs_event])
  data$treatment <- unname(c("FOLFOX alone" = "Chemo", "FOLFOX plus targeted therapy" = "Combo")[data$treatment])
  data$gender <- factor(data$gender, levels = c("Male", "Female"))
  data$kras <- unname(c("Wild-type" = 0, "Mutant" = 1)[data$kras])
  data$differentiation_grade <- unname(c("Moderately differentiated" = 0, "Poorly differentiated" = 1)[data$differentiation_grade])
  data$lymph_positive <- unname(c("Negative" = 0, "Positive" = 1)[data$lymph_positive])
  data$t_stage <- sub("^T", "", data$t_stage)
  validate_model_frame(data, "FDZS-2", 95L)
}
```

<!-- source: 05_continuous_and_binary_isr_cox.R; lines: 66-88 -->

``` r
# Fit one specified Cox model and retain its actual model identity and fitted coefficients.
# data contains one row per patient; formula defines adjustment and ISR representation.
# Returns an HR/CI table compatible with the shared forest display function.
estimate_isr_cox <- function(formula, data, cohort, specification) {
  fit <- survival::coxph(
    formula, data = data, ties = "efron", na.action = stats::na.fail,
    control = survival::coxph.control(iter.max = 100L, eps = 1e-9)
  )
  estimates <- summary(fit)$coefficients
  beta <- estimates[, "coef"]
  se <- estimates[, "se(coef)"]
  data.frame(
    cohort = cohort, specification = specification,
    model_formula = paste(deparse(formula), collapse = " "),
    model_term = names(beta), display_label = unname(term_labels[names(beta)]),
    coefficient = unname(beta), standard_error = unname(se),
    hazard_ratio = exp(beta), lower_95_ci = exp(beta - stats::qnorm(0.975) * se),
    upper_95_ci = exp(beta + stats::qnorm(0.975) * se),
    nominal_p_value = estimates[, "Pr(>|z|)"],
    n_patients = fit$n, n_events = fit$nevent, row.names = NULL
  )
}
```

<!-- source: 05_continuous_and_binary_isr_cox.R; lines: 89-155 -->

``` r
if (isTRUE(RUN_ISR_MODELS)) {
  config <- ISR_MODEL_CONFIG
  cohorts <- list(
    `FDZS-1` = normalize_discovery_cohort(
      read_supplementary_table(config$supplementary_workbook, config$discovery_sheet),
      utils::read.delim(config$discovery_isr, check.names = FALSE)
    ),
    `FDZS-2` = normalize_current_test(
      read_supplementary_table(config$supplementary_workbook, config$test_sheet)
    )
  )
  coefficient_parts <- list()
  display_parts <- list()
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  for (cohort in names(cohorts)) {
    data <- cohorts[[cohort]]
    data$isr_group <- factor(ifelse(data$isr > ISR_CUTOFF, "High", "Low"),
                             levels = c("High", "Low"))
    for (field in scaled_covariates) data[[field]] <- as.numeric(scale(data[[field]]))

    # Fit both ISR representations with the same clinical adjustment set.
    for (representation in c("continuous", "binary")) {
      score_term <- if (representation == "continuous") "isr" else "isr_group"
      formula <- stats::as.formula(paste(
        "survival::Surv(rfs_time_months, rfs_event) ~", score_term, "* treatment +",
        paste(clinical_covariates, collapse = " + ")
      ))
      estimates <- estimate_isr_cox(formula, data, cohort,
                                    paste("Adjusted", representation, "ISR-by-treatment model"))
      estimates$isr_representation <- representation
      if (representation == "binary") {
        estimates$display_label[estimates$model_term == "treatmentCombo"] <-
          "Treatment: Combo vs Chemo within High ISR"
      }
      key <- paste(cohort, representation, sep = "_")
      coefficient_parts[[key]] <- estimates
      displayed <- estimates

      # Express the binary interaction as the High/Low treatment-HR ratio.
      if (representation == "binary") {
        interaction <- displayed$model_term == "isr_groupLow:treatmentCombo"
        displayed$coefficient[interaction] <- -estimates$coefficient[interaction]
        displayed$hazard_ratio[interaction] <- 1 / estimates$hazard_ratio[interaction]
        displayed$lower_95_ci[interaction] <- 1 / estimates$upper_95_ci[interaction]
        displayed$upper_95_ci[interaction] <- 1 / estimates$lower_95_ci[interaction]
        displayed$model_term[interaction] <- "interaction_high_vs_low"
        displayed$display_label[interaction] <- "Treatment HR ratio: High / Low"
      }
      display_parts[[key]] <- displayed
    }
  }

  coefficients <- do.call(rbind, coefficient_parts)
  displayed <- do.call(rbind, display_parts)
  utils::write.table(coefficients, file.path(config$output_dir, "Figure6_continuous_binary_Cox_coefficients.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(displayed, file.path(config$output_dir, "Figure6_continuous_binary_Cox_displayed_contrasts.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  for (key in names(display_parts)) {
    result <- display_parts[[key]]
    panel <- if (result$cohort[1] == "FDZS-1") "Figure6E" else "Figure6F"
    representation <- result$isr_representation[1]
    write_forest_plot(result, paste(result$cohort[1], representation, "ISR: adjusted Cox model"),
                      file.path(config$output_dir, paste0(panel, "_", representation, "_ISR_Cox.pdf")))
  }
}
```

## Supplementary Figure 13: overall survival

The saved OS analysis uses 34 discovery patients with nine deaths and 62
test patients with 11 deaths. Its source selects a separate ISR cutoff
within each cohort’s Combo subgroup using maximally selected rank
statistics (`minprop = 0.3`). Four-group Kaplan-Meier curves and nominal
pairwise log-rank tests summarize OS.

![Complete publication Figure S13: discovery and test overall
survival](report_assets/figure_s13_publication.png)

**Figure S13. SpMap-derived ISR stratifies overall survival in CRLM.**
OS in discovery (34 patients, nine deaths; median follow-up, 34.6
months) and test (62 patients, 11 deaths), stratified by ISR and
adjuvant regimen. Red and blue denote Combo and Chemo; darker shades
indicate High ISR. Shading shows 95% confidence intervals, crosses
indicate censoring, and the lower panels give numbers at risk. Nominal
log-rank P values appear in each panel. Saved analysis PDFs:
[discovery](report_assets/figure_s13_discovery_os.pdf) and
[test](report_assets/figure_s13_test_os.pdf).

Among High-ISR patients, the treatment comparisons have nominal log-rank
P = 0.013 in discovery and P = 0.089 in test. These OS comparisons use
cohort-specific, outcome-selected ISR groups.

## Supplementary Figure 14: treatment-stratified pathology-model comparison

The saved analysis compares ISR with CHIEF, CONCH, CTransPath, and UNI
prediction groups. ISR uses the discovery-derived RFS cutoff; each
comparator uses the median of its prediction table. A
grouped-score-by-treatment Cox model supplies the interaction test, and
separate treatment-only Cox models estimate Combo-versus-Chemo HRs
within High and Low groups. The display transforms HRs to
`100 * (1 - HR)` and transforms their 95% confidence limits accordingly.

![Complete publication Figure S14: treatment-stratified ISR and
pathology-model comparison](report_assets/figure_s14_publication.png)

**Figure S14. Treatment-stratified recurrence associations of ISR and
pathology-based models.** Estimated percentage reduction in RFS hazard
associated with Combo versus Chemo within each predicted-risk group. Red
triangles and blue circles denote High and Low groups. Horizontal
intervals show transformed 95% confidence intervals; the dashed
reference marks HR = 1. Numerical columns give HR reductions, nominal
subgroup P values, and sample sizes; grouped interaction P values appear
at left. [Download the saved analysis
PDF](report_assets/figure_s14_pathology_comparison.pdf).

The High-ISR subgroup has a displayed HR reduction of 74.6% (nominal
subgroup P = 0.034; n = 24); the ISR-group-by-treatment interaction has
P = 0.116. Comparator estimates use 90 patients per model, and ISR uses
95 patients.

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

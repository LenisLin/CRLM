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

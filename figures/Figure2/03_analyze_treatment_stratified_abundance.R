#!/usr/bin/env Rscript

# Target: Figure 2D and treatment-stratified abundance analyses.
# Purpose: Relate patient-level tissue subtype fractions to clinical variables,
# recurrence, and treatment strata.
# Inputs: A SingleCellExperiment RDS with subtype, tissue, treatment, recurrence,
# and clinical metadata.
# Outputs: Statistical-result TSVs plus Figure 2D and Supplementary Figure 4/5 PDFs.
# Workflow: 1. Load and validate the SCE. 2. Derive patient tissue fractions. 3. Run
# prespecified stratum-specific tests. 4. Write results and render figure panels.

# Directly editable interactive configuration.
sce_path <- file.path("data", "FDZS1_IMC_processed.rds")
output_dir <- file.path("results", "Figure2", "treatment_stratified_abundance")

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(forestplot)
  library(grid)
})

#' Select tissue-specific subtype labels for downstream analyses.
#'
#' @param tissue Character scalar tissue code: `TC`, `IM`, or `PT`.
#' @param all_cell_types Character vector of observed subtype labels.
#' @return Character vector of unique subtype labels selected for the tissue.
cell_types_for_tissue <- function(tissue, all_cell_types) {
  # Partition observed subtype labels into biological lineages used by tissue-specific panels.
  tumor_cells <- all_cell_types[grepl("^TC_", all_cell_types)]
  stromal_cells <- all_cell_types[grepl("^SC_", all_cell_types)]
  mono_cells <- all_cell_types[grepl("^Mono_", all_cell_types)]
  macro_cells <- all_cell_types[grepl("^Macro_", all_cell_types)]
  immune_cells <- intersect(c("B", "CD4T", "CD8T", "NK", "Treg"), all_cell_types)
  # Retain only the lineage combinations defined for the requested tissue compartment.
  selected <- switch(
    tissue,
    TC = c(tumor_cells, stromal_cells, intersect(c("CD8T", "NK", "Treg", "Mono_CD11c"), all_cell_types)),
    IM = c(tumor_cells, stromal_cells, immune_cells, intersect(c("Mono_CD11c", "Macro_CD163", "Macro_CD169"), all_cell_types)),
    PT = c(stromal_cells, immune_cells, mono_cells, macro_cells)
  )
  unique(selected)
}

#' Calculate patient-level mean ROI fractions by tissue and subtype.
#'
#' @param meta Data frame of retained-cell metadata containing PID, ID, Tissue, and SubType.
#' @return Tibble with one mean ROI fraction for each patient, tissue, and subtype.
calculate_patient_tissue_fractions <- function(meta) {
  # Complete subtype counts within each ROI, then average ROI fractions at the patient-tissue level.
  meta %>%
    count(PID, ID, Tissue, SubType, name = "cell_count") %>%
    group_by(Tissue) %>%
    complete(nesting(PID, ID), SubType, fill = list(cell_count = 0L)) %>%
    ungroup() %>%
    group_by(PID, ID, Tissue) %>%
    mutate(roi_fraction = cell_count / sum(cell_count)) %>%
    ungroup() %>%
    group_by(PID, Tissue, SubType) %>%
    summarise(mean_roi_fraction = mean(roi_fraction), .groups = "drop")
}

#' Create a clinical-correlation bubble plot for selected cell subtypes.
#'
#' @param correlation_results Data frame of subtype-clinical Spearman test results.
#' @param cell_type_order Character vector specifying subtype display order.
#' @return A ggplot object encoding correlation coefficients and adjusted P-value categories.
create_correlation_plot <- function(correlation_results, cell_type_order) {
  # Convert adjusted P values to display categories while preserving the requested subtype order.
  plot_data <- correlation_results %>%
    mutate(
      p_level = case_when(
        p_adjusted <= 0.01 ~ "p.adj <= 0.01",
        p_adjusted <= 0.05 ~ "0.01 < p.adj <= 0.05",
        p_adjusted <= 0.25 ~ "0.05 < p.adj <= 0.25",
        TRUE ~ "p.adj > 0.25"
      ),
      subtype = factor(subtype, levels = rev(cell_type_order)),
      clinical_variable = factor(clinical_variable, levels = unique(clinical_variable)),
      p_level = factor(p_level, levels = c("p.adj <= 0.01", "0.01 < p.adj <= 0.05", "0.05 < p.adj <= 0.25", "p.adj > 0.25"))
    )
  # Encode association direction by color and adjusted-significance category by point size.
  ggplot(plot_data, aes(x = clinical_variable, y = subtype, color = spearman_rho, size = p_level)) +
    geom_point(alpha = 0.8) +
    scale_color_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-1, 1), name = "Correlation") +
    scale_size_manual(values = c("p.adj <= 0.01" = 8, "0.01 < p.adj <= 0.05" = 6, "0.05 < p.adj <= 0.25" = 4, "p.adj > 0.25" = 2), name = "Significance") +
    labs(title = "Correlation Between Cell Abundance and Clinical Factors", x = NULL, y = NULL) +
    theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
}

#' Create a univariable Cox-model forest plot.
#'
#' @param cox_results Data frame containing subtype hazard ratios, confidence intervals, and P values.
#' @return A forestplot graphical object drawn on the active graphics device.
create_forest_plot <- function(cox_results) {
  # Order model terms by nominal evidence and format effect estimates for the label table.
  plot_data <- cox_results %>% arrange(p_raw) %>%
    mutate(
      hr_ci = sprintf("%.2f (%.2f-%.2f)", hazard_ratio, ci_lower, ci_upper),
      p_text = scales::pvalue(p_raw, accuracy = 0.001)
    )
  labeltext <- rbind(c("Features", "HR (95% CI)", "P value"), as.matrix(plot_data[, c("subtype", "hr_ci", "p_text")]))
  # Draw hazard ratios and confidence intervals relative to the null value of one.
  forestplot::forestplot(
    labeltext = labeltext,
    mean = c(NA, plot_data$hazard_ratio), lower = c(NA, plot_data$ci_lower), upper = c(NA, plot_data$ci_upper),
    is.summary = c(TRUE, rep(FALSE, nrow(plot_data))), zero = 1, xlab = "HR (exp(coef))", title = "Univariable Cox Analysis",
    graph.pos = "right", boxsize = 0.35, col = fpColors(box = "#2E86AB", lines = "#2E86AB", summary = "#2E86AB")
  )
}

#' Run treatment-stratified abundance, association, and recurrence analyses.
#'
#' @param sce_path Character scalar path to the input SingleCellExperiment RDS.
#' @param output_dir Character scalar path to the directory for Figure 2 and supplementary outputs.
#' @return Invisibly returns the final PDF path from `ggsave()` after writing statistical tables and PDF figures.
run_treatment_stratified_abundance <- function(sce_path, output_dir) {
  # Load the SCE and require all identifiers, outcomes, treatments, and covariates used below.
  if (!file.exists(sce_path)) stop("SCE input does not exist: ", sce_path)
  sce <- readRDS(sce_path)
  continuous_vars <- c("fong_score", "TBS", "CRLM_number", "CEA", "CA199")
  discrete_vars <- c("RFS_event", "Treatment", "Gender", "KRAS_mutation", "Differential_grade", "Lymph_positive")
  required_fields <- c("PID", "ID", "Tissue", "SubType", "RFS_time", unique(c(continuous_vars, discrete_vars)))
  missing_fields <- setdiff(required_fields, colnames(colData(sce)))
  if (length(missing_fields) > 0L) stop("Missing colData field(s): ", paste(missing_fields, collapse = ", "))

  # Restrict the cell input to the three analyzed tissue compartments.
  retained <- colData(sce)$Tissue %in% c("IM", "PT", "TC")
  if (!any(retained)) stop("No IM, PT, or TC cells are available.")
  meta <- as.data.frame(colData(sce[, retained]))
  if (!all(c("Chemo", "Combo") %in% unique(meta$Treatment))) stop("Treatment must include Chemo and Combo.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Derive patient-level mean ROI fractions and attach one clinical record per patient.
  patient_tissue_fractions <- calculate_patient_tissue_fractions(meta)
  clinical_data <- meta %>% select(PID, all_of(c(continuous_vars, discrete_vars, "RFS_time"))) %>% distinct(PID, .keep_all = TRUE)
  patient_tissue_fractions <- left_join(patient_tissue_fractions, clinical_data, by = "PID")
  write.table(patient_tissue_fractions, file.path(output_dir, "patient_tissue_subtype_fractions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # Create one row per patient-tissue combination for subtype-wise tests and models.
  patient_wide <- patient_tissue_fractions %>%
    select(PID, Tissue, SubType, mean_roi_fraction) %>%
    pivot_wider(names_from = SubType, values_from = mean_roi_fraction, values_fill = 0) %>%
    left_join(clinical_data, by = "PID")
  all_cell_types <- setdiff(sort(unique(as.character(meta$SubType))), "UNKNOWN")
  treatment_groups <- list(Overall = c("Chemo", "Combo"), Chemo = "Chemo", Combo = "Combo")
  analysis_tissues <- c("IM", "TC")
  correlation_rows <- list()
  discrete_rows <- list()
  cox_rows <- list()

  # Evaluate each tissue separately within the overall cohort and each treatment stratum.
  for (treatment_name in names(treatment_groups)) {
    for (tissue_name in analysis_tissues) {
      # Select the patient-tissue analysis rows and subtype roster for this stratum.
      data_subset <- patient_wide %>% filter(Tissue == tissue_name, Treatment %in% treatment_groups[[treatment_name]])
      if (nrow(data_subset) < 10L) next
      cell_types <- cell_types_for_tissue(tissue_name, all_cell_types)
      cell_types <- intersect(cell_types, names(data_subset))

      # Test patient-level subtype fractions against each continuous clinical variable.
      correlation_stratum <- list()
      for (subtype in cell_types) {
        for (clinical_variable in continuous_vars) {
          complete <- stats::complete.cases(data_subset[, c(subtype, clinical_variable)])
          test <- cor.test(data_subset[[subtype]][complete], data_subset[[clinical_variable]][complete], method = "spearman", exact = FALSE)
          correlation_stratum[[length(correlation_stratum) + 1L]] <- data.frame(
            tissue = tissue_name, treatment_stratum = treatment_name, subtype = subtype, clinical_variable = clinical_variable,
            spearman_rho = unname(test$estimate), p_raw = test$p.value, n = sum(complete), stringsAsFactors = FALSE
          )
        }
      }
      # Control the false-discovery rate within this tissue-treatment correlation family.
      correlation_stratum <- bind_rows(correlation_stratum) %>% mutate(p_adjusted = p.adjust(p_raw, method = "BH"))
      correlation_rows[[length(correlation_rows) + 1L]] <- correlation_stratum

      # Compare subtype fractions across eligible binary clinical groups at the patient level.
      discrete_stratum <- list()
      for (clinical_variable in discrete_vars) {
        observed_levels <- sort(unique(stats::na.omit(data_subset[[clinical_variable]])))
        if (length(observed_levels) != 2L) next
        variable_rows <- list()
        for (subtype in cell_types) {
          complete <- stats::complete.cases(data_subset[, c(subtype, clinical_variable)])
          test_data <- data_subset[complete, ]
          test <- wilcox.test(test_data[[subtype]] ~ test_data[[clinical_variable]], alternative = "two.sided")
          mean_group_1 <- mean(test_data[[subtype]][test_data[[clinical_variable]] == observed_levels[[1L]]])
          mean_group_2 <- mean(test_data[[subtype]][test_data[[clinical_variable]] == observed_levels[[2L]]])
          variable_rows[[length(variable_rows) + 1L]] <- data.frame(
            tissue = tissue_name, treatment_stratum = treatment_name, clinical_variable = clinical_variable, subtype = subtype,
            group_1 = as.character(observed_levels[[1L]]), group_2 = as.character(observed_levels[[2L]]), fold_change_group2_vs_group1 = mean_group_2 / (mean_group_1 + 1e-9),
            p_raw = test$p.value, n = nrow(test_data), recurrence_count = sum(test_data$RFS_event == 1), stringsAsFactors = FALSE
          )
        }
        # Adjust subtype tests together within each binary clinical variable.
        discrete_stratum[[length(discrete_stratum) + 1L]] <- bind_rows(variable_rows) %>% mutate(p_adjusted = p.adjust(p_raw, method = "BH"))
      }
      discrete_stratum <- bind_rows(discrete_stratum)
      discrete_rows[[length(discrete_rows) + 1L]] <- discrete_stratum

      # Fit one patient-level univariable recurrence model per selected subtype.
      cox_stratum <- list()
      for (subtype in cell_types) {
        complete <- stats::complete.cases(data_subset[, c(subtype, "RFS_time", "RFS_event")])
        cox_data <- data_subset[complete, ]
        model <- coxph(Surv(cox_data$RFS_time, cox_data$RFS_event) ~ scale(cox_data[[subtype]]))
        model_summary <- summary(model)
        cox_stratum[[length(cox_stratum) + 1L]] <- data.frame(
          tissue = tissue_name, treatment_stratum = treatment_name, subtype = subtype,
          hazard_ratio = model_summary$conf.int[1L, "exp(coef)"], ci_lower = model_summary$conf.int[1L, "lower .95"], ci_upper = model_summary$conf.int[1L, "upper .95"],
          p_raw = model_summary$logtest["pvalue"], n = nrow(cox_data), recurrence_count = sum(cox_data$RFS_event == 1), stringsAsFactors = FALSE
        )
      }
      # Adjust model P values across the subtype family for this tissue-treatment stratum.
      cox_stratum <- bind_rows(cox_stratum) %>% mutate(p_adjusted = p.adjust(p_raw, method = "BH"))
      cox_rows[[length(cox_rows) + 1L]] <- cox_stratum

      # Construct the treatment-specific discrete, correlation, and Cox supplementary panels.
      if (treatment_name != "Overall") {
        figure_prefix <- if (tissue_name == "IM") "figureS4" else "figureS5"
        discrete_plot <- discrete_stratum %>%
          mutate(log2_fold_change = pmax(-1.5, pmin(1.5, log2(fold_change_group2_vs_group1))),
                 p_level = case_when(p_adjusted <= 0.01 ~ "p.adj <= 0.01", p_adjusted <= 0.05 ~ "0.01 < p.adj <= 0.05", p_adjusted <= 0.25 ~ "0.05 < p.adj <= 0.25", TRUE ~ "p.adj > 0.25"))
        p_discrete <- ggplot(discrete_plot, aes(x = subtype, y = clinical_variable, size = p_level, colour = log2_fold_change)) +
          geom_point() + scale_color_gradient2(low = "royalblue", mid = "white", high = "firebrick3", name = "Log2 FC") +
          labs(title = paste("Discrete Variable Analysis:", tissue_name, treatment_name), x = "Cell Subtype", y = "Clinical Feature") +
          theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        # Export the three aligned evidence views for this treatment and tissue.
        ggsave(file.path(output_dir, paste0(figure_prefix, "A_discrete_", tissue_name, "_", treatment_name, ".pdf")), p_discrete, width = 10, height = 6)
        ggsave(file.path(output_dir, paste0(figure_prefix, "B_correlation_", tissue_name, "_", treatment_name, ".pdf")), create_correlation_plot(correlation_stratum, cell_types), width = 8, height = 6)
        pdf(file.path(output_dir, paste0(figure_prefix, "C_univariable_cox_", tissue_name, "_", treatment_name, ".pdf")), width = 10, height = 8)
        create_forest_plot(cox_stratum)
        dev.off()
      }
    }
  }

  # Consolidate and export all treatment-tissue statistical results.
  correlation_results <- bind_rows(correlation_rows)
  discrete_results <- bind_rows(discrete_rows)
  cox_results <- bind_rows(cox_rows)
  write.table(correlation_results, file.path(output_dir, "treatment_stratified_spearman_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(discrete_results, file.path(output_dir, "treatment_stratified_wilcoxon_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(cox_results, file.path(output_dir, "treatment_stratified_univariable_cox_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # Select the predefined IM subtypes and reshape patient fractions for Figure 2D.
  figure2d_subtypes <- c("CD4T", "CD8T", "Macro_CD169", "Macro_CD163")
  missing_figure2d_subtypes <- setdiff(figure2d_subtypes, names(patient_wide))
  if (length(missing_figure2d_subtypes) > 0L) stop("Figure 2D subtype(s) are absent: ", paste(missing_figure2d_subtypes, collapse = ", "))
  figure2d_data <- patient_wide %>% filter(Tissue == "IM", Treatment %in% c("Chemo", "Combo")) %>%
    select(PID, Treatment, RFS_event, all_of(figure2d_subtypes)) %>%
    pivot_longer(cols = all_of(figure2d_subtypes), names_to = "subtype", values_to = "mean_roi_fraction")
  figure2d_results <- figure2d_data %>%
    group_by(Treatment, subtype) %>%
    #' Test recurrence-group abundance differences within one treatment-subtype group.
    #'
    #' @param data Tibble containing one treatment and subtype group of patient fractions.
    #' @param key One-row tibble identifying the current group.
    #' @return One-row data frame of Wilcoxon results, or a zero-row data frame when the group is insufficient.
    group_modify(function(data, key) {
      # Test only treatment-subtype groups with both recurrence states sufficiently represented.
      group_counts <- table(data$RFS_event)
      if (nrow(data) < 6L || length(group_counts) != 2L || min(group_counts) < 3L) return(data.frame())
      test <- wilcox.test(mean_roi_fraction ~ RFS_event, data = data, alternative = "two.sided")
      data.frame(p_raw = test$p.value, n = nrow(data), recurrence_count = sum(data$RFS_event == 1))
    }) %>% ungroup()
  if (nrow(figure2d_results) == 0L) stop("Figure 2D requires at least six patients and three patients per RFS_event group in each treatment-subtype comparison.")
  # Export Figure 2D test results and construct their panel-specific annotations.
  write.table(figure2d_results, file.path(output_dir, "figure2D_IM_recurrence_wilcoxon.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  figure2d_labels <- figure2d_results %>% mutate(label = paste0("p = ", format.pval(p_raw, digits = 3, eps = 0.001)), y = max(figure2d_data$mean_roi_fraction) * 1.05)
  # Plot patient observations and recurrence-group distributions within each treatment facet.
  p_figure2d <- ggplot(figure2d_data, aes(x = factor(RFS_event), y = mean_roi_fraction)) +
    geom_jitter(aes(color = factor(RFS_event)), width = 0.1, size = 2, alpha = 0.6) +
    geom_boxplot(aes(fill = factor(RFS_event)), alpha = 0.5, outlier.shape = NA) +
    facet_grid(subtype ~ Treatment, scales = "free_y") +
    geom_text(data = figure2d_labels, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 3) +
    labs(x = "RFS_event (0 = no recurrence, 1 = recurrence)", y = "Mean ROI cell fraction", fill = "RFS_event", color = "RFS_event") +
    theme_bw(base_size = 14) + theme(legend.position = "top", strip.text.y = element_text(angle = 0))
  # Export the completed Figure 2D recurrence-abundance panel.
  ggsave(file.path(output_dir, "figure2D_IM_recurrence_boxplots.pdf"), p_figure2d, width = 6, height = 12)
}

if (sys.nframe() == 0L) {
  run_treatment_stratified_abundance(sce_path, output_dir)
}

#!/usr/bin/env Rscript

# Target: Supplementary Figure 3 tissue-composition panels associated with Figure 2.
# Purpose: Quantify subtype composition across tissues, ROIs, and patients.
# Inputs: A SingleCellExperiment RDS with current annotations and clinical metadata.
# Outputs: Three composition TSV files and three Supplementary Figure 3 PDF panels.
# Workflow: 1. Load and validate the SCE. 2. Retain IM/PT/TC cells. 3. Calculate
# pooled, ROI, and patient fractions. 4. Write tables and render composition plots.

# Directly editable interactive configuration.
sce_path <- file.path("data", "FDZS1_IMC_processed.rds")
output_dir <- file.path("results", "Figure2", "supplementary_figureS3")

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})

#' Create tissue, ROI, and patient cell-subtype composition outputs.
#'
#' @param sce_path Character scalar path to the input SingleCellExperiment RDS.
#' @param output_dir Character scalar path to the directory for Supplementary Figure 3 outputs.
#' @return Integer device number returned by `dev.off()` after writing composition tables and PDF plots.
run_tissue_composition <- function(sce_path, output_dir) {
  # Load the SCE and require the identifiers and clinical fields used by all three panels.
  if (!file.exists(sce_path)) stop("SCE input does not exist: ", sce_path)
  sce <- readRDS(sce_path)
  clinical_fields <- c(
    "PID", "ID", "Tissue", "SubType", "RFS_event", "Gender", "KRAS_mutation",
    "CRC_site", "Differential_grade", "T_stage", "Lymph_positive",
    "RFS_time", "Age", "TBS", "CRLM_number", "CRLM_size", "CEA", "CA199"
  )
  missing_fields <- setdiff(clinical_fields, colnames(colData(sce)))
  if (length(missing_fields) > 0L) stop("Missing colData field(s): ", paste(missing_fields, collapse = ", "))

  # Restrict composition summaries to cells from the displayed TC, IM, and PT compartments.
  retained <- colData(sce)$Tissue %in% c("IM", "PT", "TC")
  if (!any(retained)) stop("No IM, PT, or TC cells are available.")
  sce <- sce[, retained]
  meta <- as.data.frame(colData(sce))
  tissue_order <- c("TC", "IM", "PT")
  if (!all(tissue_order %in% meta$Tissue)) stop("All TC, IM, and PT tissues are required.")

  # Preserve the shared subtype color mapping across pooled, ROI, and patient views.
  subtype_palette <- metadata(sce)$color_vectors$SubType
  if (is.null(subtype_palette)) stop("metadata(sce)$color_vectors$SubType is required.")
  observed_subtypes <- sort(unique(as.character(meta$SubType)))
  if (any(!observed_subtypes %in% names(subtype_palette))) stop("SubType palette is incomplete.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Aggregate all retained cells within tissue to obtain the pooled panel fractions.
  pooled_fractions <- meta %>%
    count(Tissue, SubType, name = "cell_count") %>%
    group_by(Tissue) %>%
    mutate(cell_fraction = cell_count / sum(cell_count)) %>%
    ungroup() %>%
    mutate(Tissue = factor(Tissue, levels = tissue_order), SubType = factor(SubType, levels = observed_subtypes)) %>%
    arrange(Tissue, SubType)
  write.table(pooled_fractions, file.path(output_dir, "figureS3A_pooled_tissue_subtype_fractions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  #' Build one pooled tissue-composition pie chart.
  #'
  #' @param tissue Character scalar tissue code included in `tissue_order`.
  #' @return A ggplot object showing subtype fractions for the requested tissue.
  pie_plots <- lapply(tissue_order, function(tissue) {
    # Build one composition panel from the pooled fractions for a single tissue.
    plot_data <- filter(pooled_fractions, Tissue == tissue)
    ggplot(plot_data, aes(x = "", y = cell_fraction, fill = SubType)) +
      geom_col(width = 1, color = "black", show.legend = TRUE) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = subtype_palette) +
      geom_text(aes(label = ifelse(cell_fraction > 0.05, scales::percent(cell_fraction, accuracy = 1), "")), position = position_stack(vjust = 0.5), size = 3) +
      labs(title = tissue) + theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), legend.position = "top", legend.title = element_blank())
  })
  # Export the aligned tissue-level pie charts as Supplementary Figure 3A.
  pdf(file.path(output_dir, "figureS3A_tissue_composition_pie_charts.pdf"), width = 20, height = 10)
  print(cowplot::plot_grid(plotlist = pie_plots, ncol = 3))
  dev.off()

  # Calculate subtype fractions independently within each patient ROI.
  roi_fractions <- meta %>%
    count(PID, ID, Tissue, SubType, name = "cell_count") %>%
    group_by(PID, ID, Tissue) %>%
    mutate(cell_fraction = cell_count / sum(cell_count)) %>%
    ungroup() %>%
    mutate(Tissue = factor(Tissue, levels = tissue_order), SubType = factor(SubType, levels = observed_subtypes)) %>%
    arrange(Tissue, ID, SubType)
  write.table(roi_fractions, file.path(output_dir, "figureS3B_roi_subtype_fractions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # Construct aligned ROI composition bars and a tissue strip on the same ROI order.
  roi_order <- roi_fractions %>% distinct(ID, Tissue) %>% arrange(Tissue, ID)
  stacked_bars <- ggplot(roi_fractions, aes(x = factor(ID, levels = roi_order$ID), y = cell_fraction, fill = SubType)) +
    geom_col(color = "black", width = 1) + scale_fill_manual(values = subtype_palette) +
    labs(x = "", y = "Cell Fraction", title = "Cell Fractions in Each ROI by Tissue") +
    theme_minimal() + theme(axis.text.x = element_blank(), panel.grid = element_blank(), legend.position = "top", legend.title = element_blank())
  tissue_bar <- ggplot(roi_order, aes(x = factor(ID, levels = ID), y = 0.5, fill = Tissue)) +
    geom_col(width = 1) + scale_fill_brewer(palette = "Set2") + theme_minimal() +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom")
  # Export the stacked ROI composition view as Supplementary Figure 3B.
  pdf(file.path(output_dir, "figureS3B_roi_subtype_fractions.pdf"), width = 15, height = 8)
  print(stacked_bars / tissue_bar + plot_layout(heights = c(14, 1)))
  dev.off()

  # Aggregate cells by patient and complete absent subtype combinations with zero fractions.
  patient_order <- meta %>% distinct(PID, RFS_event) %>% arrange(RFS_event) %>% pull(PID)
  patient_fractions <- meta %>%
    count(PID, SubType, name = "cell_count") %>%
    group_by(PID) %>%
    mutate(cell_fraction = cell_count / sum(cell_count)) %>%
    ungroup() %>%
    complete(PID = patient_order, SubType = observed_subtypes, fill = list(cell_count = 0, cell_fraction = 0)) %>%
    mutate(PID = factor(PID, levels = patient_order), SubType = factor(SubType, levels = observed_subtypes)) %>%
    arrange(PID, SubType)
  write.table(patient_fractions, file.path(output_dir, "figureS3C_patient_subtype_fractions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # Align the patient-by-subtype matrix with one clinical annotation record per patient.
  fraction_table <- patient_fractions %>%
    select(PID, SubType, cell_fraction) %>%
    pivot_wider(names_from = SubType, values_from = cell_fraction)
  rownames(fraction_table) <- NULL
  fraction_matrix <- fraction_table %>% column_to_rownames("PID") %>% as.matrix()
  clinical_annotation <- meta %>%
    select(all_of(clinical_fields[-c(2L, 3L, 4L)])) %>%
    distinct(PID, .keep_all = TRUE) %>%
    slice(match(patient_order, PID))
  rownames(clinical_annotation) <- NULL
  clinical_annotation <- clinical_annotation %>% column_to_rownames("PID")
  discrete_vars <- c("RFS_event", "Gender", "KRAS_mutation", "CRC_site", "Differential_grade", "T_stage", "Lymph_positive")
  continuous_vars <- c("RFS_time", "Age", "TBS", "CRLM_number", "CRLM_size", "CEA", "CA199")
  discrete_colors <- list(
    RFS_event = c("0" = "#bc3c29ff", "1" = "#0072b5ff"), Gender = c("0" = "#7876b1ff", "1" = "#ee4c97ff"),
    KRAS_mutation = c("0" = "grey", "1" = "#e18727ff"),
    CRC_site = c("1" = "lightgreen", "2" = "forestgreen", "3" = "darkgreen"), Differential_grade = c("0" = "lightblue", "1" = "darkblue"),
    T_stage = c("1" = "lemonchiffon", "2" = "lemonchiffon2", "3" = "lemonchiffon3", "4" = "lemonchiffon4"), Lymph_positive = c("0" = "white", "1" = "grey60")
  )
  #' Create a continuous clinical-variable color function.
  #'
  #' @param variable Character scalar naming a continuous clinical annotation column.
  #' @return A circlize color-mapping function spanning the observed annotation range.
  # Map each continuous variable across its observed patient-level range.
  continuous_colors <- lapply(continuous_vars, function(variable) colorRamp2(range(clinical_annotation[[variable]], na.rm = TRUE), c("white", "red4")))
  names(continuous_colors) <- continuous_vars
  # Construct clinical annotations and the composition heatmap on the shared patient axis.
  clinical_ha <- HeatmapAnnotation(df = clinical_annotation[, discrete_vars, drop = FALSE], col = discrete_colors, which = "row") +
    HeatmapAnnotation(df = clinical_annotation[, continuous_vars, drop = FALSE], col = continuous_colors, which = "row")
  heatmap_colors <- colorRamp2(c(0, max(fraction_matrix)), c("white", "navy"))
  clinical_heatmap <- Heatmap(fraction_matrix, name = "Fraction", col = heatmap_colors, show_row_names = TRUE, show_column_names = TRUE, cluster_rows = FALSE, cluster_columns = TRUE, column_names_side = "bottom")
  # Export the patient-level composition and clinical annotation panel.
  pdf(file.path(output_dir, "figureS3C_patient_clinical_composition_heatmap.pdf"), width = 15, height = 15)
  draw(clinical_ha + clinical_heatmap)
  dev.off()
}

if (sys.nframe() == 0L) {
  run_tissue_composition(sce_path, output_dir)
}

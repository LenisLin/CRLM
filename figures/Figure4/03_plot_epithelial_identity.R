#!/usr/bin/env Rscript

# Figure 4D and Supplementary Figure 9A-B epithelial identity panels.
#
# Input contract: one SingleCellExperiment RDS whose expression assay contains
# PRPS1, FASN, GLUT1, HK2, Ki67, VEGF, and CAIX, and whose colData contains
# unique CellID plus complete PID, ID, Tissue, SubType, and MajorType2 fields.
# The analyzed tissues are TC, IM, and PT; the epithelial mask is MajorType2 ==
# EC with the four configured TargetSubTypes. Output contract: four TSVs named
# figureS9A_epithelial_tissue_distribution.tsv (tissue),
# epithelial_roi_marker_means.tsv (patient-ROI-tissue-marker),
# epithelial_patient_tissue_marker_means.tsv (patient-tissue-marker), and
# epithelial_patient_tissue_wilcoxon.tsv (marker-tissue-pair), plus
# FigureS9A_epithelial_tissue_distribution.pdf,
# Figure4D_CAIX_Ki67_expression.pdf, and
# FigureS9B_epithelial_metabolic_marker_expression.pdf. The inferential unit is
# one patient-tissue marker mean formed by averaging that patient's ROI means.
#
# Edit these paths and set RUN_EPITHELIAL_IDENTITY to TRUE before direct
# execution. This no-CLI entry point intentionally does not read data by
# default, because the configured SCE is a large project object.
RUN_EPITHELIAL_IDENTITY <- FALSE

EPITHELIAL_IDENTITY_CONFIG <- list(
    sce_input = "path/to/FDZS1_IMC_processed.rds",
    output_dir = "path/to/Figure4_epithelial_identity_outputs",
    tissues = c("TC", "IM", "PT"),
    TargetSubTypes = c("TC_CAIX", "TC_EpCAM", "TC_Ki67", "TC_VEGF"),
    MetaMarkers = c("PRPS1", "FASN", "GLUT1", "HK2", "Ki67", "VEGF", "CAIX"),
    figure4d_markers = c("CAIX", "Ki67"),
    supplementary_figure9b_markers = c("VEGF", "FASN", "GLUT1", "HK2", "PRPS1")
)

TISSUE_COLORS <- c("TC" = "#1B9E77", "IM" = "#D95F02", "PT" = "#7570B3")
PAIRWISE_TISSUE_COMPARISONS <- list(c("TC", "IM"), c("TC", "PT"), c("IM", "PT"))

#' Verify packages used by the Figure 4D and Supplementary Figure 9 panels.
#'
#' @return Invisibly returns NULL or stops before the SCE is read.
require_epithelial_identity_packages <- function() {
    packages <- c("SingleCellExperiment", "SummarizedExperiment", "ggplot2", "scales")
    missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
    if (length(missing) > 0L) {
        stop("Missing required packages: ", paste(missing, collapse = ", "))
    }
}

#' Validate the configured Figure 4 epithelial identity contract.
#'
#' @param sce SingleCellExperiment loaded from the configured RDS.
#' @param config Editable Figure 4 epithelial identity configuration.
#' @return Cell metadata aligned to the SCE columns.
validate_epithelial_identity_input <- function(sce, config) {
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("Configured input must be a SingleCellExperiment")
    }
    metadata <- as.data.frame(SummarizedExperiment::colData(sce))
    required_fields <- c("CellID", "PID", "ID", "Tissue", "SubType", "MajorType2")
    missing_fields <- setdiff(required_fields, names(metadata))
    if (length(missing_fields) > 0L) {
        stop("SCE colData is missing: ", paste(missing_fields, collapse = ", "))
    }
    if (anyDuplicated(metadata$CellID)) {
        stop("CellID must be unique")
    }
    if (any(!stats::complete.cases(metadata[, required_fields, drop = FALSE]))) {
        stop("Required epithelial identity metadata fields must be complete")
    }
    if (!all(config$tissues %in% as.character(metadata$Tissue))) {
        stop("Configured tissues are not all present in the SCE")
    }
    if (!all(config$TargetSubTypes %in% as.character(metadata$SubType))) {
        stop("Configured TargetSubTypes are not all present in the SCE")
    }
    available_markers <- rownames(sce)
    required_markers <- unique(c(
        config$MetaMarkers,
        config$figure4d_markers,
        config$supplementary_figure9b_markers
    ))
    missing_markers <- setdiff(required_markers, available_markers)
    if (length(missing_markers) > 0L) {
        stop("Configured markers are not present in the SCE assay: ",
             paste(missing_markers, collapse = ", "))
    }
    metadata
}

#' Summarize the epithelial cell distribution used by Supplementary Figure 9A.
#'
#' @param metadata Validated SCE cell metadata.
#' @param config Editable Figure 4 epithelial identity configuration.
#' @return One row per displayed tissue with epithelial cell counts and fractions.
summarize_epithelial_distribution <- function(metadata, config) {
    is_epithelial <- as.character(metadata$MajorType2) == "EC" &
        as.character(metadata$SubType) %in% config$TargetSubTypes &
        as.character(metadata$Tissue) %in% config$tissues
    counts <- table(factor(as.character(metadata$Tissue[is_epithelial]), levels = config$tissues))
    if (sum(counts) == 0L) {
        stop("No configured epithelial cells are available in the configured tissues")
    }
    data.frame(
        Tissue = config$tissues,
        epithelial_cell_count = as.integer(counts),
        epithelial_cell_fraction = as.numeric(counts) / sum(counts),
        stringsAsFactors = FALSE
    )
}

#' Derive expression means at the historical ROI boundary without testing ROIs.
#'
#' @param sce SingleCellExperiment with validated feature names.
#' @param metadata Validated SCE cell metadata.
#' @param config Editable Figure 4 epithelial identity configuration.
#' @return One row per patient, ROI, tissue, and marker.
summarize_epithelial_roi_expression <- function(sce, metadata, config) {
    keep <- as.character(metadata$MajorType2) == "EC" &
        as.character(metadata$SubType) %in% config$TargetSubTypes &
        as.character(metadata$Tissue) %in% config$tissues
    if (!any(keep)) {
        stop("No cells meet the configured epithelial identity mask")
    }
    expression <- as.matrix(SummarizedExperiment::assay(sce)[config$MetaMarkers, keep, drop = FALSE])
    if (any(!is.finite(expression))) {
        stop("Configured epithelial marker expression contains non-finite values")
    }
    selected_metadata <- metadata[keep, c("PID", "ID", "Tissue"), drop = FALSE]
    roi_key <- interaction(
        selected_metadata$PID, selected_metadata$ID, selected_metadata$Tissue,
        drop = TRUE, lex.order = TRUE
    )
    roi_metadata <- unique(data.frame(
        roi_key = as.character(roi_key),
        PID = as.character(selected_metadata$PID),
        ID = as.character(selected_metadata$ID),
        Tissue = as.character(selected_metadata$Tissue),
        stringsAsFactors = FALSE
    ))
    roi_metadata <- roi_metadata[match(levels(roi_key), roi_metadata$roi_key), , drop = FALSE]
    marker_means <- t(vapply(
        levels(roi_key),
        function(key) rowMeans(expression[, roi_key == key, drop = FALSE]),
        numeric(length(config$MetaMarkers))
    ))
    colnames(marker_means) <- config$MetaMarkers
    wide <- cbind(roi_metadata[, c("PID", "ID", "Tissue"), drop = FALSE], marker_means)
    long <- stats::reshape(
        wide,
        varying = config$MetaMarkers,
        v.names = "expression_mean",
        timevar = "marker",
        times = config$MetaMarkers,
        direction = "long"
    )
    rownames(long) <- NULL
    long$Tissue <- factor(as.character(long$Tissue), levels = config$tissues)
    long[order(long$PID, long$ID, long$Tissue, long$marker), , drop = FALSE]
}

#' Average ROI means to one patient-tissue-marker value for statistical testing.
#'
#' @param roi_expression One row per patient, ROI, tissue, and marker.
#' @param config Editable Figure 4 epithelial identity configuration.
#' @return One row per patient, tissue, and marker with ROI count and mean expression.
summarize_epithelial_patient_expression <- function(roi_expression, config) {
    groups <- interaction(roi_expression$PID, roi_expression$Tissue, roi_expression$marker,
                          drop = TRUE, lex.order = TRUE)
    grouped <- split(roi_expression, groups)
    patient_expression <- do.call(rbind, lapply(grouped, function(data) {
        data.frame(
            PID = as.character(data$PID[[1L]]),
            Tissue = as.character(data$Tissue[[1L]]),
            marker = as.character(data$marker[[1L]]),
            roi_count = nrow(data),
            expression_mean = mean(data$expression_mean),
            stringsAsFactors = FALSE
        )
    }))
    rownames(patient_expression) <- NULL
    patient_expression$Tissue <- factor(patient_expression$Tissue, levels = config$tissues)
    patient_expression[order(patient_expression$PID, patient_expression$Tissue,
                            patient_expression$marker), , drop = FALSE]
}

#' Calculate nominal two-sided Wilcoxon tests for each marker and tissue pair.
#'
#' @param patient_expression One row per patient, tissue, and marker.
#' @param config Editable Figure 4 epithelial identity configuration.
#' @return One row per marker and pairwise tissue comparison; P values are nominal.
calculate_epithelial_wilcoxon_tests <- function(patient_expression, config) {
    rows <- list()
    for (marker in config$MetaMarkers) {
        marker_data <- patient_expression[patient_expression$marker == marker, , drop = FALSE]
        for (comparison in PAIRWISE_TISSUE_COMPARISONS) {
            group_1 <- marker_data$expression_mean[marker_data$Tissue == comparison[[1L]]]
            group_2 <- marker_data$expression_mean[marker_data$Tissue == comparison[[2L]]]
            if (length(group_1) == 0L || length(group_2) == 0L) {
                stop("Each marker comparison requires both configured tissue groups")
            }
            test <- stats::wilcox.test(group_1, group_2, alternative = "two.sided", exact = FALSE)
            rows[[length(rows) + 1L]] <- data.frame(
                marker = marker,
                tissue_1 = comparison[[1L]],
                tissue_2 = comparison[[2L]],
                n_patients_tissue_1 = length(group_1),
                n_patients_tissue_2 = length(group_2),
                median_tissue_1 = stats::median(group_1),
                median_tissue_2 = stats::median(group_2),
                wilcoxon_w = unname(test$statistic),
                nominal_p_value = test$p.value,
                stringsAsFactors = FALSE
            )
        }
    }
    do.call(rbind, rows)
}

#' Format a nominal P value for a panel annotation.
#'
#' @param p_value Numeric nominal P value.
#' @return Character display label.
format_nominal_p_value <- function(p_value) {
    if (p_value <= 0.0001) return("****")
    if (p_value <= 0.001) return("***")
    if (p_value <= 0.01) return("**")
    if (p_value <= 0.05) return("*")
    "n.s."
}

#' Plot one marker subset using patient-tissue values and nominal Wilcoxon labels.
#'
#' @param patient_expression One row per patient, tissue, and marker.
#' @param tests Nominal Wilcoxon table from `calculate_epithelial_wilcoxon_tests()`.
#' @param markers Marker subset assigned to one PDF.
#' @param config Editable Figure 4 epithelial identity configuration.
#' @return A ggplot object with one violin panel per marker.
plot_epithelial_markers <- function(patient_expression, tests, markers, config) {
    plot_data <- patient_expression[patient_expression$marker %in% markers, , drop = FALSE]
    plot_data$marker <- factor(plot_data$marker, levels = markers)
    labels <- tests[tests$marker %in% markers, , drop = FALSE]
    labels$marker <- factor(labels$marker, levels = markers)
    marker_maximum <- tapply(plot_data$expression_mean, plot_data$marker, max)
    marker_range <- tapply(plot_data$expression_mean, plot_data$marker, range)
    labels$comparison_index <- ave(
        seq_len(nrow(labels)), labels$marker, FUN = seq_along
    )
    labels$y_step <- vapply(as.character(labels$marker), function(marker) {
        values <- marker_range[[marker]]
        maximum <- marker_maximum[[marker]]
        max(diff(values) * 0.10, abs(maximum) * 0.05, 0.05)
    }, numeric(1L))
    labels$y <- vapply(seq_len(nrow(labels)), function(index) {
        marker <- as.character(labels$marker[[index]])
        marker_maximum[[marker]] + labels$y_step[[index]] *
            labels$comparison_index[[index]]
    }, numeric(1L))
    labels$x_start <- match(labels$tissue_1, config$tissues)
    labels$x_end <- match(labels$tissue_2, config$tissues)
    labels$x_label <- (labels$x_start + labels$x_end) / 2
    labels$label <- vapply(labels$nominal_p_value, format_nominal_p_value, character(1L))

    ggplot2::ggplot(plot_data, ggplot2::aes(x = Tissue, y = expression_mean, fill = Tissue)) +
        ggplot2::geom_violin(alpha = 0.70, trim = TRUE, color = "black") +
        ggplot2::geom_jitter(width = 0.08, height = 0, size = 1.4, alpha = 0.65) +
        ggplot2::geom_segment(
            data = labels,
            ggplot2::aes(x = x_start, xend = x_end, y = y, yend = y),
            inherit.aes = FALSE
        ) +
        ggplot2::geom_text(
            data = labels,
            ggplot2::aes(x = x_label, y = y, label = label),
            inherit.aes = FALSE,
            size = 3,
            vjust = -0.3
        ) +
        ggplot2::facet_wrap(~marker, scales = "free_y", ncol = length(markers)) +
        ggplot2::scale_fill_manual(values = TISSUE_COLORS) +
        ggplot2::labs(x = NULL, y = "Patient mean of ROI mean expression") +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(legend.position = "none", strip.background = ggplot2::element_blank())
}

#' Plot the epithelial regional distribution for Supplementary Figure 9A.
#'
#' @param distribution One row per tissue with epithelial cell fractions.
#' @return A ggplot object.
plot_epithelial_distribution <- function(distribution) {
    ggplot2::ggplot(distribution, ggplot2::aes(x = "", y = epithelial_cell_fraction, fill = Tissue)) +
        ggplot2::geom_col(width = 1, color = "black") +
        ggplot2::coord_polar("y", start = 0) +
        ggplot2::geom_text(
            ggplot2::aes(label = scales::percent(epithelial_cell_fraction, accuracy = 1)),
            position = ggplot2::position_stack(vjust = 0.5),
            size = 3
        ) +
        ggplot2::scale_fill_manual(values = TISSUE_COLORS) +
        ggplot2::theme_void() +
        ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "top")
}

#' Write a tab-separated output table.
#'
#' @param data Data frame to write.
#' @param path Destination path.
#' @return Invisibly returns NULL after writing.
write_epithelial_identity_tsv <- function(data, path) {
    utils::write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
    invisible(NULL)
}

#' Generate Figure 4D and the linked Supplementary Figure 9A-B outputs.
#'
#' @param config Editable input, output, marker, and tissue configuration.
#' @return Invisibly returns the distribution, ROI, patient, and test tables.
run_epithelial_identity <- function(config = EPITHELIAL_IDENTITY_CONFIG) {
    require_epithelial_identity_packages()
    if (!file.exists(config$sce_input)) {
        stop("SCE input does not exist: ", config$sce_input)
    }
    sce <- readRDS(config$sce_input)
    metadata <- validate_epithelial_identity_input(sce, config)
    distribution <- summarize_epithelial_distribution(metadata, config)
    roi_expression <- summarize_epithelial_roi_expression(sce, metadata, config)
    patient_expression <- summarize_epithelial_patient_expression(roi_expression, config)
    tests <- calculate_epithelial_wilcoxon_tests(patient_expression, config)

    dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
    write_epithelial_identity_tsv(distribution, file.path(
        config$output_dir, "figureS9A_epithelial_tissue_distribution.tsv"
    ))
    write_epithelial_identity_tsv(roi_expression, file.path(
        config$output_dir, "epithelial_roi_marker_means.tsv"
    ))
    write_epithelial_identity_tsv(patient_expression, file.path(
        config$output_dir, "epithelial_patient_tissue_marker_means.tsv"
    ))
    write_epithelial_identity_tsv(tests, file.path(
        config$output_dir, "epithelial_patient_tissue_wilcoxon.tsv"
    ))

    ggplot2::ggsave(
        file.path(config$output_dir, "FigureS9A_epithelial_tissue_distribution.pdf"),
        plot_epithelial_distribution(distribution), width = 6, height = 6
    )
    ggplot2::ggsave(
        file.path(config$output_dir, "Figure4D_CAIX_Ki67_expression.pdf"),
        plot_epithelial_markers(patient_expression, tests, config$figure4d_markers, config),
        width = 7, height = 4.5
    )
    ggplot2::ggsave(
        file.path(config$output_dir, "FigureS9B_epithelial_metabolic_marker_expression.pdf"),
        plot_epithelial_markers(
            patient_expression, tests, config$supplementary_figure9b_markers, config
        ),
        width = 13, height = 4.5
    )
    invisible(list(
        distribution = distribution,
        roi_expression = roi_expression,
        patient_expression = patient_expression,
        tests = tests
    ))
}

if (isTRUE(RUN_EPITHELIAL_IDENTITY)) {
    run_epithelial_identity()
}

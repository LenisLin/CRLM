# Target: Figure 7C-E and Supplementary Figures 17-18 spatial composition outputs.
# Purpose: map registered SpMap tiles to Bin100 spots and compare RCTD compositions.
# Inputs: edit the paths in CONFIG for the ST bundle, tile annotations, RCTD weights, helpers, and output.
# Outputs: mapping records, composition and test TSVs, and Figure 7/Supplementary PDF panels.
# Ordered workflow: load ST and tiles, map classes, intersect RCTD weights, summarize composition, test, plot.

EXPECTED_ALL_SPOTS <- 18186L
EXPECTED_RCTD_SPOTS <- 17902L
EXPECTED_ALL_COUNTS <- c("0" = 463L, "1" = 382L, "2" = 1526L, "3" = 15815L)
EXPECTED_RCTD_COUNTS <- c("0" = 449L, "1" = 377L, "2" = 1522L, "3" = 15554L)

CLASS_LABELS <- c(
    "0" = "TC", "1" = "PIR-Niche", "2" = "PSM-Niche", "3" = "Other PT"
)
CLASS_COLORS <- c(
    "0" = "#E41A1C", "1" = "#1F78B4", "2" = "#FF7F00", "3" = "grey70"
)

MAIN_CELLTYPE_MAPPING <- list(
    B = c("RCTD_B_AIM2", "RCTD_B_HSP", "RCTD_B_IGHM"),
    Plasma = "RCTD_Plasma",
    Epithelial = c("RCTD_Epithelial", "RCTD_Malignant"),
    Hepatocyte = c("RCTD_Hepatocyte", "RCTD_Cholangiocyte"),
    Fibroblast = "RCTD_Fibroblast",
    Endothelial = "RCTD_Endothelial",
    Myeloid = c(
        "RCTD_Macro_CD163", "RCTD_Macro_CD169", "RCTD_Macro_SPP1",
        "RCTD_Monocyte", "RCTD_cDC_CD1C", "RCTD_cDC_CLEC9A",
        "RCTD_Neutrophil", "RCTD_Myeloid_Other", "RCTD_pDC_LILRA4"
    ),
    NK = c("RCTD_NK_FCGR3A", "RCTD_NK_KLRB1"),
    T = c(
        "RCTD_CD4T_IL7R", "RCTD_CD8T_CXCL13", "RCTD_CD8T_GZMK",
        "RCTD_MAIT_KLRB1", "RCTD_T_HSP", "RCTD_T_MKI67",
        "RCTD_T_TCF7", "RCTD_Treg_FOXP3"
    )
)

IMMUNE_SUBTYPE_MAPPING <- list(
    B_AIM2 = "RCTD_B_AIM2",
    B_IGHM = "RCTD_B_IGHM",
    CD4_IL7R = "RCTD_CD4T_IL7R",
    CD8_CXCL13 = "RCTD_CD8T_CXCL13",
    CD8_GZMK = "RCTD_CD8T_GZMK",
    MAIT_KLRB1 = "RCTD_MAIT_KLRB1",
    Macro_CD163 = "RCTD_Macro_CD163",
    Macro_CD169 = "RCTD_Macro_CD169",
    Macro_SPP1 = "RCTD_Macro_SPP1",
    Plasma = "RCTD_Plasma"
)

MAIN_COLORS <- c(
    B = "#FFD700", Plasma = "#FFFF99", Epithelial = "#33A02C",
    Hepatocyte = "#B2DF8A", Fibroblast = "#FF7F00",
    Endothelial = "#E31A1C", Myeloid = "#6A3D9A",
    NK = "#A6CEE3", T = "#1F78B4"
)

IMMUNE_COLORS <- c(
    B_AIM2 = "#FFD700", B_IGHM = "#FFDAB9", CD4_IL7R = "#1F78B4",
    CD8_CXCL13 = "#3C8DBC", CD8_GZMK = "#5DA5D8",
    MAIT_KLRB1 = "#7EBDEF", Macro_CD163 = "#6A3D9A",
    Macro_CD169 = "#8A63A9", Macro_SPP1 = "#AA89B8", Plasma = "#FFFF00"
)

FOCUSED_CELLTYPES <- c(
    "B_IGHM", "CD8T_GZMK", "Fibroblast", "Macro_CD163", "Macro_SPP1"
)

# Edit these paths for the Figure 7 analysis.
CONFIG <- list(
    st_bundle = "/path/to/fdzs4_bin100_bundle",
    tile_annotations = "/path/to/registered_spmap_tiles.tsv",
    rctd_weights = "/path/to/figure7_rctd_output/rctd_normalized_weights.tsv",
    output_dir = "/path/to/figure7_composition_output",
    st_functions = "src/st/st_functions.R"
)

#' Purpose: read a CSV, TSV, or TXT table while retaining an optional string identifier.
#'
#' @param path Character scalar path to the input table.
#' @param first_column_character Logical scalar indicating whether to read the first
#'   column as character data.
#' @return A data frame containing the parsed table.
read_delimited <- function(path, first_column_character = FALSE) {
    # Infer the delimiter once and preserve source column names for downstream identity matching.
    separator <- if (tolower(tools::file_ext(path)) %in% c("tsv", "txt")) {
        "\t"
    } else {
        ","
    }
    header <- utils::read.table(
        path,
        header = TRUE,
        sep = separator,
        nrows = 0,
        quote = "\"",
        comment.char = "",
        check.names = FALSE
    )
    column_classes <- rep(NA_character_, ncol(header))
    if (first_column_character) {
        column_classes[[1]] <- "character"
    }
    utils::read.table(
        path,
        header = TRUE,
        sep = separator,
        quote = "\"",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE,
        colClasses = column_classes
    )
}

#' Purpose: write a data frame as a tab-separated output table.
#'
#' @param data Data frame to write.
#' @param path Character scalar destination path for the TSV file.
#' @return Invisibly returns `NULL` after writing the table.
write_tsv <- function(data, path) {
    utils::write.table(
        data,
        file = path,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
}

#' Purpose: verify class counts against the specified Figure 7 population contract.
#'
#' @param classes Vector of mapped class labels.
#' @param expected Named vector of expected class counts.
#' @param expected_total Integer scalar expected population size.
#' @param population Character scalar naming the checked population.
#' @return Invisibly returns `NULL` when the observed counts match the contract.
assert_class_counts <- function(classes, expected, expected_total, population) {
    # Count mapped Bin100 spots on a fixed four-class axis, including classes with zero observations.
    observed <- table(factor(as.character(classes), levels = names(expected)))
    if (length(classes) != expected_total ||
        !identical(as.integer(observed), as.integer(expected))) {
        stop(
            population, " does not match the fixed Figure 7 class contract. ",
            "Observed total/classes: ", length(classes), "/",
            paste(as.integer(observed), collapse = ","),
            call. = FALSE
        )
    }
}

#' Purpose: record the tile-selection route for every mapped spatial spot.
#'
#' @param spatial_df Data frame containing spatial spot identifiers and coordinates.
#' @param tile_anno Data frame containing tile boundaries and predicted classes.
#' @param assigned_classes Vector of mapped class labels, one per spatial spot.
#' @return A spot-level mapping-record data frame.
build_mapping_audit <- function(spatial_df, tile_anno, assigned_classes) {
    # Index candidate tiles by sample so spatial containment is evaluated only within each specimen.
    tile_rows_by_sample <- split(seq_len(nrow(tile_anno)), tile_anno$sample_id)
    records <- vector("list", nrow(spatial_df))

    # Record one selected registered tile per Bin100 spot, using containment before nearest-top-left.
    for (i in seq_len(nrow(spatial_df))) {
        tile_rows <- tile_rows_by_sample[[as.character(spatial_df$sample_id[[i]])]]
        sample_tiles <- tile_anno[tile_rows, , drop = FALSE]
        in_tile <- which(
            spatial_df$x[[i]] >= sample_tiles$tile_left &
                spatial_df$x[[i]] <= sample_tiles$tile_right &
                spatial_df$y[[i]] >= sample_tiles$tile_top &
                spatial_df$y[[i]] <= sample_tiles$tile_bottom
        )
        if (length(in_tile) > 0) {
            selected_local <- in_tile[[1]]
            assignment_mode <- "contained"
        } else {
            distances <- sqrt(
                (sample_tiles$tile_left - spatial_df$x[[i]])^2 +
                    (sample_tiles$tile_top - spatial_df$y[[i]])^2
            )
            selected_local <- which.min(distances)
            assignment_mode <- "nearest_top_left"
        }
        selected_row <- tile_rows[[selected_local]]
        selected_distance <- sqrt(
            (tile_anno$tile_left[[selected_row]] - spatial_df$x[[i]])^2 +
                (tile_anno$tile_top[[selected_row]] - spatial_df$y[[i]])^2
        )
        records[[i]] <- data.frame(
            assignment_mode = assignment_mode,
            containing_tile_count = length(in_tile),
            selected_tile_source_row = selected_row,
            distance_to_selected_top_left = selected_distance,
            stringsAsFactors = FALSE
        )
    }

    # Reconcile the independently recorded tile route with the assigned class for every spot.
    audit <- cbind(spatial_df, do.call(rbind, records))
    audit$pred_class <- as.integer(assigned_classes)
    selected_classes <- as.integer(as.character(
        tile_anno$pred_class[audit$selected_tile_source_row]
    ))
    if (!identical(selected_classes, audit$pred_class)) {
        stop("Mapping audit disagrees with assign_predictions_to_spots", call. = FALSE)
    }
    audit$class_label <- unname(CLASS_LABELS[as.character(audit$pred_class)])
    audit
}

#' Purpose: summarize mapped SpMap classes and assignment modes by sample.
#'
#' @param audit Spot-level mapping records returned by `build_mapping_audit()`.
#' @return A data frame containing per-sample and total mapping summaries.
summarize_mapping <- function(audit) {
    # Aggregate spot-level mapping outcomes by sample while retaining the all-sample denominator.
    grouped <- split(audit, audit$sample_id)
    rows <- lapply(names(grouped), function(sample_id) {
        data <- grouped[[sample_id]]
        data.frame(
            sample_id = sample_id,
            all_spots = nrow(data),
            TC = sum(data$pred_class == 0L),
            PIR = sum(data$pred_class == 1L),
            PSM = sum(data$pred_class == 2L),
            other_PT = sum(data$pred_class == 3L),
            contained = sum(data$assignment_mode == "contained"),
            nearest_top_left = sum(data$assignment_mode == "nearest_top_left"),
            multiple_containing_tiles = sum(data$containing_tile_count > 1L),
            stringsAsFactors = FALSE
        )
    })
    total <- data.frame(
        sample_id = "Total",
        all_spots = nrow(audit),
        TC = sum(audit$pred_class == 0L),
        PIR = sum(audit$pred_class == 1L),
        PSM = sum(audit$pred_class == 2L),
        other_PT = sum(audit$pred_class == 3L),
        contained = sum(audit$assignment_mode == "contained"),
        nearest_top_left = sum(audit$assignment_mode == "nearest_top_left"),
        multiple_containing_tiles = sum(audit$containing_tile_count > 1L),
        stringsAsFactors = FALSE
    )
    do.call(rbind, c(rows, list(total)))
}

#' Purpose: calculate per-sample cell-type composition fractions from RCTD proportions.
#'
#' @param proportions Data frame containing sample IDs and cell-type proportions.
#' @return A long-format composition data frame with normalized display fractions.
summarize_composition <- function(proportions) {
    # Average spot-level RCTD proportions within each sample before display normalization.
    value_columns <- setdiff(colnames(proportions), "sample_id")
    means <- stats::aggregate(
        proportions[, value_columns, drop = FALSE],
        by = list(sample_id = proportions$sample_id),
        FUN = mean
    )
    long <- tidyr::pivot_longer(
        means,
        cols = -sample_id,
        names_to = "cell_type",
        values_to = "mean_proportion"
    )
    # Normalize only across displayed cell types within a sample so each stacked bar sums to one.
    totals <- ave(long$mean_proportion, long$sample_id, FUN = sum)
    long$display_fraction <- long$mean_proportion / totals
    long
}

#' Purpose: construct a stacked composition plot.
#'
#' @param data Long-format composition data returned by `summarize_composition()`.
#' @param colors Named character vector mapping cell types to colors.
#' @param title Character scalar plot title.
#' @param legend_title Character scalar legend title.
#' @return A ggplot composition object.
make_composition_plot <- function(data, colors, title, legend_title) {
    # Keep cell-type factor order synchronized with the named palette across samples.
    data$cell_type <- factor(data$cell_type, levels = names(colors))
    ggplot2::ggplot(
        data,
        ggplot2::aes(x = display_fraction, y = sample_id, fill = cell_type)
    ) +
        ggplot2::geom_col(width = 0.8) +
        ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
        ggplot2::scale_x_continuous(
            labels = scales::percent,
            expand = c(0, 0)
        ) +
        ggplot2::labs(
            title = title,
            x = "Proportion",
            y = "Sample ID",
            fill = legend_title
        ) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            legend.position = "right"
        )
}

#' Purpose: test spot-level RCTD abundance between PIR and PSM niches.
#'
#' @param long_data Data frame containing cell type, niche, and abundance columns.
#' @return A BH-adjusted Wilcoxon-test data frame with annotation positions.
run_spot_wilcoxon <- function(long_data) {
    # Compare PIR and PSM at the Bin100-spot level separately for each RCTD cell type.
    grouped_data <- dplyr::group_by(long_data, cell_type)
    tests <- rstatix::wilcox_test(
        grouped_data,
        abundance ~ niche,
        paired = FALSE,
        alternative = "two.sided"
    )
    tests <- rstatix::adjust_pvalue(tests, method = "BH")
    tests <- rstatix::add_significance(tests, p.col = "p.adj")
    # Derive annotation heights from each cell type independently of the test statistics.
    maxima <- stats::aggregate(
        abundance ~ cell_type,
        data = long_data,
        FUN = max
    )
    names(maxima)[[2]] <- "maximum"
    tests <- merge(tests, maxima, by = "cell_type", sort = FALSE)
    tests$y.position <- pmax(0.02, tests$maximum * 1.12)
    tests$analysis_unit <- "Bin100_spot"
    tests
}

#' Purpose: construct faceted spot-level abundance violin plots.
#'
#' @param data Data frame containing spot-level cell-type abundances.
#' @param tests Data frame containing P-value annotations and positions.
#' @param title Character scalar plot title.
#' @param columns Integer scalar number of facet columns.
#' @return A ggplot abundance-comparison object.
make_spot_plot <- function(data, tests, title, columns) {
    # Display spot-level abundance distributions with their matching per-cell-type tests.
    ggplot2::ggplot(
        data,
        ggplot2::aes(x = niche, y = abundance, fill = niche)
    ) +
        ggplot2::geom_violin(alpha = 0.7, width = 0.8, trim = TRUE) +
        ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
        ggpubr::stat_pvalue_manual(
            tests,
            label = "p",
            tip.length = 0.01,
            hide.ns = FALSE,
            inherit.aes = FALSE
        ) +
        ggplot2::facet_wrap(~cell_type, scales = "free_y", ncol = columns) +
        ggplot2::scale_fill_manual(
            values = c("PIR-Niche" = "#377EB8", "PSM-Niche" = "#FDBF6F")
        ) +
        ggplot2::labs(
            title = title,
            x = "Spatial niche",
            y = "RCTD proportion",
            fill = "Spatial niche"
        ) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            legend.position = "bottom"
        ) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18)))
}

#' Purpose: execute the configured SpMap mapping and RCTD composition workflow.
#'
#' @param config List containing ST, tile, RCTD, helper, and output paths.
#' @return Invisibly returns the path to the final saved PDF after writing all outputs.
main <- function(config = CONFIG) {
    # Load the complete Bin100 spatial object; its columns define the all-spot mapping population.
    output_dir <- config[["output_dir"]]
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    source(config[["st_functions"]])
    spe <- read_anndata_to_spe(config[["st_bundle"]], verbose = TRUE)
    if (ncol(spe) != EXPECTED_ALL_SPOTS) {
        stop(
            "The all-spot Bin100 input contains ", EXPECTED_ALL_SPOTS,
            " spots; observed ", ncol(spe), call. = FALSE
        )
    }

    # Map registered SpMap tiles to spot coordinates within samples and retain each selected tile.
    tile_anno <- read_delimited(config[["tile_annotations"]])
    spatial_coordinates <- SpatialExperiment::spatialCoords(spe)
    spatial_df <- data.frame(
        spot_id = as.character(colnames(spe)),
        x = spatial_coordinates[, "x"],
        y = spatial_coordinates[, "y"],
        sample_id = as.character(
            SummarizedExperiment::colData(spe)[["sample_id"]]
        ),
        stringsAsFactors = FALSE
    )
    assigned_classes <- assign_predictions_to_spots(spatial_df, tile_anno)
    assert_class_counts(
        assigned_classes,
        EXPECTED_ALL_COUNTS,
        EXPECTED_ALL_SPOTS,
        "All-spot SpMap mapping"
    )
    mapping_audit <- build_mapping_audit(
        spatial_df,
        tile_anno,
        assigned_classes
    )
    spe$pred_class <- assigned_classes
    spe$class_label <- unname(CLASS_LABELS[as.character(assigned_classes)])

    # Align normalized RCTD rows to spatial spot identities, yielding the retained-spot intersection.
    weight_table <- read_delimited(
        config[["rctd_weights"]],
        first_column_character = TRUE
    )
    if (!"spot_id" %in% colnames(weight_table)) {
        stop("RCTD weights must contain a spot_id column", call. = FALSE)
    }
    if (anyDuplicated(weight_table$spot_id)) {
        stop("RCTD weights contain duplicate spot_id values", call. = FALSE)
    }
    weight_columns <- setdiff(colnames(weight_table), "spot_id")
    weights <- data.frame(
        lapply(weight_table[, weight_columns, drop = FALSE], as.numeric),
        check.names = FALSE
    )
    rownames(weights) <- weight_table$spot_id
    spe_rctd <- combine_deconvolution(spe, weights, deconv_method = "RCTD")
    assert_class_counts(
        spe_rctd$pred_class,
        EXPECTED_RCTD_COUNTS,
        EXPECTED_RCTD_SPOTS,
        "RCTD-intersection SpMap mapping"
    )

    # Export both all-spot mapping records and the fixed RCTD-intersection analysis population.
    retained_ids <- as.character(colnames(spe_rctd))
    mapping_audit$rctd_status <- ifelse(
        mapping_audit$spot_id %in% retained_ids,
        "retained",
        "excluded"
    )
    write_tsv(
        mapping_audit,
        file.path(output_dir, "spmap_spot_assignments.tsv")
    )
    write_tsv(
        summarize_mapping(mapping_audit),
        file.path(output_dir, "spmap_mapping_summary.tsv")
    )

    intersection_spots <- data.frame(
        spot_id = retained_ids,
        sample_id = as.character(spe_rctd$sample_id),
        pred_class = as.integer(spe_rctd$pred_class),
        class_label = as.character(spe_rctd$class_label),
        stringsAsFactors = FALSE
    )
    write_tsv(
        intersection_spots,
        file.path(output_dir, "rctd_intersection_spots.tsv")
    )
    denominator_summary <- data.frame(
        population = c("all_mapped_spots", "rctd_intersection"),
        total_spots = c(EXPECTED_ALL_SPOTS, EXPECTED_RCTD_SPOTS),
        TC = c(EXPECTED_ALL_COUNTS[["0"]], EXPECTED_RCTD_COUNTS[["0"]]),
        PIR = c(EXPECTED_ALL_COUNTS[["1"]], EXPECTED_RCTD_COUNTS[["1"]]),
        PSM = c(EXPECTED_ALL_COUNTS[["2"]], EXPECTED_RCTD_COUNTS[["2"]]),
        other_PT = c(EXPECTED_ALL_COUNTS[["3"]], EXPECTED_RCTD_COUNTS[["3"]])
    )
    write_tsv(
        denominator_summary,
        file.path(output_dir, "rctd_intersection_summary.tsv")
    )

    # Plot mapped classes from all spots, preserving each sample's native spatial coordinates.
    map_plot <- ggplot2::ggplot(
        mapping_audit,
        ggplot2::aes(x = x, y = y, color = factor(pred_class))
    ) +
        ggplot2::geom_point(size = 0.3, alpha = 0.8) +
        ggplot2::scale_color_manual(
            values = CLASS_COLORS,
            labels = CLASS_LABELS,
            name = "Region type"
        ) +
        ggplot2::scale_y_reverse() +
        ggplot2::facet_wrap(~sample_id, ncol = 3, scales = "free") +
        ggplot2::labs(x = "X coordinate", y = "Y coordinate") +
        ggplot2::theme_minimal(base_size = 9) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(face = "bold"),
            aspect.ratio = 1,
            legend.position = "bottom"
        )
    ggplot2::ggsave(
        file.path(output_dir, "Figure7C_FigureS17_spmap_class_maps.pdf"),
        map_plot,
        width = 12,
        height = 12
    )

    # Aggregate subtype weights only after restricting to spots retained by RCTD.
    col_data <- as.data.frame(SummarizedExperiment::colData(spe_rctd))
    required_rctd_columns <- unique(c(
        unlist(MAIN_CELLTYPE_MAPPING),
        unlist(IMMUNE_SUBTYPE_MAPPING)
    ))
    missing_rctd_columns <- setdiff(required_rctd_columns, colnames(col_data))
    if (length(missing_rctd_columns) > 0) {
        stop(
            "RCTD weights are missing required subtype column(s): ",
            paste(missing_rctd_columns, collapse = ", "),
            call. = FALSE
        )
    }

    main_proportions <- aggregate_celltype_proportions(
        col_data,
        MAIN_CELLTYPE_MAPPING
    )
    immune_proportions <- aggregate_celltype_proportions(
        col_data,
        IMMUNE_SUBTYPE_MAPPING
    )
    main_summary <- summarize_composition(main_proportions)
    immune_summary <- summarize_composition(immune_proportions)
    write_tsv(
        main_summary,
        file.path(output_dir, "major_celltype_composition_by_sample.tsv")
    )
    write_tsv(
        immune_summary,
        file.path(output_dir, "immune_subtype_composition_by_sample.tsv")
    )

    # Export sample-level major and immune composition panels from the same retained spots.
    major_plot <- make_composition_plot(
        main_summary,
        MAIN_COLORS,
        "Cellular composition across spatial transcriptomics samples",
        "Cell type"
    )
    immune_plot <- make_composition_plot(
        immune_summary,
        IMMUNE_COLORS,
        "Immune-cell subtype composition across spatial transcriptomics samples",
        "Immune subtype"
    )
    ggplot2::ggsave(
        file.path(output_dir, "FigureS18A_major_cell_composition.pdf"),
        major_plot,
        width = 8,
        height = 6
    )
    ggplot2::ggsave(
        file.path(output_dir, "Figure7D_immune_subtype_composition.pdf"),
        immune_plot,
        width = 9,
        height = 6
    )

    # Reshape PIR and PSM spot weights for cell-type-specific spot-level comparisons.
    rctd_columns <- grep("^RCTD_", colnames(col_data), value = TRUE)
    niche_data <- col_data[col_data$pred_class %in% c(1L, 2L), , drop = FALSE]
    abundance_long <- tidyr::pivot_longer(
        niche_data[, c("sample_id", "pred_class", rctd_columns), drop = FALSE],
        cols = tidyselect::all_of(rctd_columns),
        names_to = "cell_type",
        values_to = "abundance"
    )
    abundance_long$cell_type <- sub("^RCTD_", "", abundance_long$cell_type)
    abundance_long$niche <- factor(
        abundance_long$pred_class,
        levels = c(1L, 2L),
        labels = c("PIR-Niche", "PSM-Niche")
    )
    spot_tests <- run_spot_wilcoxon(abundance_long)
    write_tsv(
        spot_tests,
        file.path(output_dir, "celltype_abundance_spot_wilcoxon.tsv")
    )

    # Produce the prespecified focused panel and the complete subtype companion panel.
    focused_data <- abundance_long[
        abundance_long$cell_type %in% FOCUSED_CELLTYPES,
        ,
        drop = FALSE
    ]
    focused_tests <- spot_tests[
        spot_tests$cell_type %in% FOCUSED_CELLTYPES,
        ,
        drop = FALSE
    ]
    focused_data$cell_type <- factor(
        focused_data$cell_type,
        levels = FOCUSED_CELLTYPES
    )
    focused_tests$cell_type <- factor(
        focused_tests$cell_type,
        levels = FOCUSED_CELLTYPES
    )
    focused_plot <- make_spot_plot(
        focused_data,
        focused_tests,
        "Selected RCTD proportions in PIR-Niche and PSM-Niche spots",
        5
    )
    all_subtype_plot <- make_spot_plot(
        abundance_long,
        spot_tests,
        "RCTD subtype proportions in PIR-Niche and PSM-Niche spots",
        6
    )
    ggplot2::ggsave(
        file.path(output_dir, "Figure7E_focused_celltypes_spot_level.pdf"),
        focused_plot,
        width = 12,
        height = 6
    )
    ggplot2::ggsave(
        file.path(output_dir, "FigureS18B_all_subtypes_spot_level.pdf"),
        all_subtype_plot,
        width = 16,
        height = 12
    )
}

if (sys.nframe() == 0L) {
    main()
}

# Target: Figure 7F-H spot-level DEG and GSEA outputs.
# Purpose: compare PIR-Niche and PSM-Niche Bin100 spots and display pathway enrichment.
# Inputs: edit the paths in CONFIG for the ST bundle, labels, gene sets, manifest, helpers, and output.
# Outputs: DEG, GSEA, resource, and summary TSVs plus Figure 7F-H PDF panels.
# Ordered workflow: validate labels, select niches, run DEG, build ranks, run GSEA, export and plot.

EXPECTED_ALL_SPOTS <- 18186L
EXPECTED_RCTD_SPOTS <- 17902L
EXPECTED_PIR_SPOTS <- 377L
EXPECTED_PSM_SPOTS <- 1522L
EXPECTED_GSEA_RANK_GENES <- 2559L
GSEA_COLLECTIONS <- c("Hallmark", "KEGG", "Reactome")

FIGURE7H_PATHWAYS <- c(
    Hallmark = "HALLMARK_ANGIOGENESIS",
    Reactome = "REACTOME_MTORC1_MEDIATED_SIGNALLING"
)

PATHWAY_KEYWORDS <- c(
    "bile", "cholangiocyte", "cholestasis", "liver", "hepatic",
    "jaundice", "cholesterol", "lipid metabolism",
    "immune", "immunity", "inflammatory", "interferon", "cytokine",
    "interleukin", "chemokine", "t cell", "b cell", "nk cell",
    "dendritic", "macrophage", "leukocyte", "lymphocyte",
    "tnf", "nfkb", "nf kappa b", "jak stat", "pd1", "pd l1",
    "ctla4", "checkpoint", "immunotherapy", "antigen presentation",
    "extracellular matrix", "ecm", "collagen", "fibrillin", "elastin",
    "fibronectin", "laminin", "integrin", "focal adhesion",
    "cell adhesion", "cell migration", "wound healing",
    "fibrosis", "scarring", "stroma", "mesenchymal",
    "epithelial mesenchymal transition", "emt",
    "drug", "xenobiotic", "detoxification", "cyp450", "cytochrome",
    "metabolism", "metabolic", "chemotherapy", "5 fu", "fluorouracil",
    "oxaliplatin", "irinotecan", "bevacizumab", "cetuximab",
    "targeted therapy", "egfr", "vegf", "angiogenesis",
    "apoptosis", "cell death", "autophagy", "senescence",
    "cell cycle", "proliferation", "dna repair", "dna damage",
    "hypoxia", "oxygen", "reactive oxygen", "oxidative stress",
    "wnt", "notch", "hedgehog", "tgf beta", "tgfb", "bmp",
    "pi3k", "akt", "mtor", "ras", "mapk", "erk", "myc",
    "metastasis", "invasion", "migration", "motility"
)

# Edit these paths for the Figure 7 analysis.
CONFIG <- list(
    st_bundle = "/path/to/fdzs4_bin100_bundle",
    spot_labels = "/path/to/figure7_composition_output/rctd_intersection_spots.tsv",
    gene_sets = "/path/to/offline_gene_sets.rds",
    resource_manifest = "/path/to/offline_gene_set_resource_manifest.tsv",
    output_dir = "/path/to/figure7_deg_gsea_output",
    st_functions = "src/st/st_functions.R"
)

#' Purpose: read a TSV while retaining the spot identifier as a character value.
#'
#' @param path Character scalar path to the input TSV file.
#' @return A data frame containing the parsed table.
read_tsv_with_character_id <- function(path) {
    # Preserve the leading spot identifier as text so it can align exactly to ST columns.
    header <- utils::read.delim(
        path,
        nrows = 0,
        check.names = FALSE
    )
    column_classes <- rep(NA_character_, ncol(header))
    column_classes[[1]] <- "character"
    utils::read.delim(
        path,
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

#' Purpose: construct the specified filtered and ranked GSEA gene table.
#'
#' @param de_results Data frame containing Seurat differential-expression results.
#' @param expected_rows Integer scalar expected number of retained genes.
#' @return The descending log2-fold-change rank table with filter fields.
build_gsea_rank_table <- function(de_results, expected_rows = EXPECTED_GSEA_RANK_GENES) {
    # Require the gene-level effect and detection fractions used by the prespecified rank filter.
    required <- c("gene", "avg_log2FC", "pct.1", "pct.2")
    missing <- setdiff(required, colnames(de_results))
    if (length(missing) > 0) {
        stop(
            "DEG results are missing GSEA filter column(s): ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
    # Filter genes by spot detection contrast, then rank the retained genes by PIR-versus-PSM effect.
    rank_table <- de_results
    rank_table$delta_pct <- abs(rank_table$pct.1 - rank_table$pct.2)
    rank_table$filter_pct <-
        rank_table$delta_pct < rank_table$pct.1 &
        rank_table$delta_pct < rank_table$pct.2
    rank_table <- rank_table[!rank_table$filter_pct, , drop = FALSE]
    rank_table <- rank_table[!is.na(rank_table$avg_log2FC), , drop = FALSE]
    rank_table <- rank_table[
        order(rank_table$avg_log2FC, decreasing = TRUE),
        ,
        drop = FALSE
    ]
    if (nrow(rank_table) != expected_rows) {
        stop(
            "The pre-GSEA detection-fraction filter retains ",
            expected_rows, " genes; observed ", nrow(rank_table),
            call. = FALSE
        )
    }
    rank_table$rank_position <- seq_len(nrow(rank_table))
    rank_table
}

#' Purpose: validate the configured offline pathway collections and version manifest.
#'
#' @param pathways Named list of pathway collections.
#' @param manifest Data frame recording the version of each pathway collection.
#' @return Invisibly returns `TRUE` when the resource contract is satisfied.
validate_gene_set_resource <- function(pathways, manifest) {
    # Keep pathway collections and their version rows in a one-to-one named correspondence.
    if (!is.list(pathways) || length(pathways) != length(GSEA_COLLECTIONS) ||
        anyDuplicated(names(pathways)) ||
        !setequal(names(pathways), GSEA_COLLECTIONS)) {
        stop(
            "The offline gene-set RDS must contain exactly Hallmark, KEGG, and Reactome",
            call. = FALSE
        )
    }
    required_manifest <- c("collection", "resource_version")
    missing <- setdiff(required_manifest, colnames(manifest))
    if (length(missing) > 0) {
        stop(
            "Resource manifest is missing column(s): ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
    if (anyDuplicated(manifest$collection) ||
        !setequal(manifest$collection, GSEA_COLLECTIONS) ||
        any(is.na(manifest$resource_version) | manifest$resource_version == "")) {
        stop(
            "Resource manifest must provide one versioned row for each configured collection",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

#' Purpose: convert GSEA list columns into TSV-compatible fields.
#'
#' @param results An fgsea-compatible result table.
#' @return A data frame suitable for TSV export.
export_gsea_results <- function(results) {
    # Flatten each gene-level leading edge into one field while preserving one row per pathway.
    exported <- as.data.frame(results)
    if ("leadingEdge" %in% colnames(exported)) {
        exported$leadingEdge <- vapply(
            exported$leadingEdge,
            paste,
            collapse = ";",
            FUN.VALUE = character(1)
        )
    }
    exported
}

#' Purpose: execute the configured Figure 7 PIR-versus-PSM DEG and GSEA workflow.
#'
#' @param config List containing ST, label, resource, helper, and output paths.
#' @return Invisibly returns `NULL` after writing DEG/GSEA tables and PDF panels.
main <- function(config = CONFIG) {
    # Load the retained-spot label table; each row is one unique Bin100 analysis unit.
    output_dir <- config[["output_dir"]]
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    source(config[["st_functions"]])
    spot_labels <- read_tsv_with_character_id(config[["spot_labels"]])
    required_labels <- c("spot_id", "sample_id", "pred_class", "class_label")
    missing_labels <- setdiff(required_labels, colnames(spot_labels))
    if (length(missing_labels) > 0) {
        stop(
            "Spot-label table is missing column(s): ",
            paste(missing_labels, collapse = ", "),
            call. = FALSE
        )
    }
    if (nrow(spot_labels) != EXPECTED_RCTD_SPOTS ||
        anyDuplicated(spot_labels$spot_id)) {
        stop(
            "The DEG/GSEA spot-label input must contain 17,902 unique RCTD-intersection spots",
            call. = FALSE
        )
    }
    niche_counts <- table(factor(spot_labels$pred_class, levels = c(1L, 2L)))
    if (!identical(
        as.integer(niche_counts),
        c(EXPECTED_PIR_SPOTS, EXPECTED_PSM_SPOTS)
    )) {
        stop(
            "The DEG/GSEA intersection contains 377 PIR and 1,522 PSM spots",
            call. = FALSE
        )
    }

    # Align the fixed RCTD-intersection labels to the complete ST object by spot identity.
    spe <- read_anndata_to_spe(config[["st_bundle"]], verbose = TRUE)
    if (ncol(spe) != EXPECTED_ALL_SPOTS) {
        stop(
            "The all-spot Bin100 input contains ", EXPECTED_ALL_SPOTS,
            " spots; observed ", ncol(spe), call. = FALSE
        )
    }
    label_positions <- match(spot_labels$spot_id, colnames(spe))
    if (anyNA(label_positions)) {
        stop("Some RCTD-intersection spot IDs are absent from the Bin100 input", call. = FALSE)
    }
    spe_intersection <- spe[, label_positions]
    spe_intersection$pred_class <- as.integer(spot_labels$pred_class)
    spe_intersection$class_label <- as.character(spot_labels$class_label)
    niche_mask <- spe_intersection$pred_class %in% c(1L, 2L)
    spe_niche <- spe_intersection[, niche_mask]

    # Compare PIR and PSM spots at the Bin100 level using the shared gene-count matrix.
    seurat_object <- Seurat::as.Seurat(spe_niche, counts = "counts", data = NULL)
    seurat_object$BDME_class <- factor(
        ifelse(spe_niche$pred_class == 1L, "PIR-Niche", "PSM-Niche"),
        levels = c("PIR-Niche", "PSM-Niche")
    )
    Seurat::Idents(seurat_object) <- "BDME_class"
    de_results <- Seurat::FindMarkers(
        seurat_object,
        ident.1 = "PIR-Niche",
        ident.2 = "PSM-Niche",
        test.use = "wilcox",
        min.pct = 0.01,
        logfc.threshold = 0,
        only.pos = FALSE,
        verbose = TRUE
    )
    de_results$gene <- rownames(de_results)
    de_results$de_status <- ifelse(
        de_results$p_val < 0.05 & de_results$avg_log2FC > 0,
        "Up in PIR",
        ifelse(
            de_results$p_val < 0.05 & de_results$avg_log2FC < 0,
            "Up in PSM",
            "Not significant"
        )
    )
    write_tsv(
        de_results,
        file.path(output_dir, "Task2_Differential_Expression_All_Genes.tsv")
    )

    # Export DEG population sizes, test settings, and nominal direction counts.
    de_summary <- data.frame(
        metric = c(
            "analysis_unit", "rctd_intersection_spots", "pir_spots", "psm_spots",
            "tested_genes", "nominal_p_lt_0.05_pir_up",
            "nominal_p_lt_0.05_psm_up", "test", "min_pct", "logfc_threshold"
        ),
        value = c(
            "Bin100_spot", EXPECTED_RCTD_SPOTS, EXPECTED_PIR_SPOTS,
            EXPECTED_PSM_SPOTS, nrow(de_results),
            sum(de_results$de_status == "Up in PIR"),
            sum(de_results$de_status == "Up in PSM"),
            "Seurat Wilcoxon rank-sum", 0.01, 0
        ),
        stringsAsFactors = FALSE
    )
    write_tsv(de_summary, file.path(output_dir, "deg_run_summary.tsv"))

    # Display the exported gene-level fold changes and nominal P values on bounded plot axes.
    volcano_data <- de_results
    volcano_data$plot_log2fc <- pmax(-8, pmin(8, volcano_data$avg_log2FC))
    top_genes <- volcano_data[volcano_data$p_val < 0.05, , drop = FALSE]
    top_genes <- head(
        top_genes[order(abs(top_genes$plot_log2fc), decreasing = TRUE), , drop = FALSE],
        20
    )
    volcano_plot <- ggplot2::ggplot(
        volcano_data,
        ggplot2::aes(x = plot_log2fc, y = -log10(p_val), color = de_status)
    ) +
        ggplot2::geom_point(alpha = 0.6, size = 1) +
        ggplot2::scale_color_manual(values = c(
            "Up in PIR" = "#377EB8",
            "Up in PSM" = "#FDBF6F",
            "Not significant" = "grey70"
        )) +
        ggplot2::geom_hline(
            yintercept = -log10(0.05),
            linetype = "dashed",
            color = "red"
        ) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
        ggplot2::geom_text(
            data = top_genes,
            ggplot2::aes(label = gene),
            size = 3,
            vjust = -0.5,
            check_overlap = TRUE
        ) +
        ggplot2::coord_cartesian(xlim = c(-8, 8)) +
        ggplot2::labs(
            title = "Differential expression: PIR-Niche versus PSM-Niche",
            x = "Average log2 fold change (PIR/PSM)",
            y = "-Log10 nominal P value",
            color = "Expression"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            legend.position = "bottom",
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )
    ggplot2::ggsave(
        file.path(output_dir, "Figure7F_DEG_volcano.pdf"),
        volcano_plot,
        width = 10,
        height = 8
    )

    # Build the single ordered gene statistic shared by all configured pathway collections.
    rank_table <- build_gsea_rank_table(de_results)
    write_tsv(rank_table, file.path(output_dir, "gsea_rank_table_2559_genes.tsv"))
    gene_ranks <- stats::setNames(rank_table$avg_log2FC, rank_table$gene)

    # Load offline gene sets and require one versioned manifest row per collection.
    pathway_databases <- readRDS(config[["gene_sets"]])
    resource_manifest <- utils::read.delim(
        config[["resource_manifest"]],
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    validate_gene_set_resource(pathway_databases, resource_manifest)

    # Run pathway enrichment independently by collection against the same Bin100-derived rank vector.
    gsea_results <- lapply(GSEA_COLLECTIONS, function(database) {
        run_gsea_analysis(
            gene_ranks,
            pathway_databases[[database]],
            database
        )
    })
    names(gsea_results) <- GSEA_COLLECTIONS
    combined_gsea <- as.data.frame(do.call(rbind, gsea_results))
    combined_gsea_export <- export_gsea_results(combined_gsea)
    write_tsv(
        combined_gsea_export,
        file.path(output_dir, "Task3_GSEA_PIR_vs_PSM_All_Results.tsv")
    )

    # Select nominally significant, topic-matched pathways for display while retaining all results.
    significant <- combined_gsea[combined_gsea$pval <= 0.05, , drop = FALSE]
    keyword_pattern <- paste(PATHWAY_KEYWORDS, collapse = "|")
    selected <- significant[
        grepl(keyword_pattern, significant$pathway_clean, ignore.case = TRUE),
        ,
        drop = FALSE
    ]
    selected <- selected[order(selected$NES, decreasing = TRUE), , drop = FALSE]
    if (nrow(selected) == 0) {
        stop("No nominally significant pathways match the display filter", call. = FALSE)
    }
    write_tsv(
        export_gsea_results(selected),
        file.path(output_dir, "Figure7G_selected_pathways.tsv")
    )

    # Encode selected pathway direction, nominal significance, and collection in one bubble panel.
    bubble_data <- as.data.frame(selected)
    bubble_data$neg_log_pval <- -log10(bubble_data$pval)
    bubble_data$pathway_display <- tools::toTitleCase(
        gsub("_", " ", bubble_data$pathway)
    )
    bubble_data$pathway_display <- factor(
        bubble_data$pathway_display,
        levels = bubble_data$pathway_display[order(bubble_data$NES)]
    )
    bubble_plot <- ggplot2::ggplot(
        bubble_data,
        ggplot2::aes(x = NES, y = pathway_display)
    ) +
        ggplot2::geom_point(
            ggplot2::aes(size = neg_log_pval, fill = database),
            shape = 21,
            alpha = 0.8
        ) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
        ggplot2::scale_fill_brewer(palette = "Set2", name = "Database") +
        ggplot2::labs(
            title = "Selected pathway enrichment in PIR-Niche versus PSM-Niche",
            x = "Normalized enrichment score",
            y = "Pathway",
            size = "-Log10 nominal P"
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )
    ggplot2::ggsave(
        file.path(output_dir, "Figure7G_GSEA_bubble_plot.pdf"),
        bubble_plot,
        width = 14,
        height = 10
    )

    # Require both prespecified Figure 7H pathways in the resources and computed results.
    for (database in names(FIGURE7H_PATHWAYS)) {
        pathway <- FIGURE7H_PATHWAYS[[database]]
        if (!pathway %in% names(pathway_databases[[database]])) {
            stop(
                "Figure 7H pathway is absent from the offline resource: ", pathway,
                call. = FALSE
            )
        }
        if (!pathway %in% combined_gsea$pathway) {
            stop("Figure 7H pathway is absent from GSEA results: ", pathway, call. = FALSE)
        }
    }
    #' Purpose: construct an enrichment curve for one prespecified Figure 7H pathway.
    #'
    #' @param database Character scalar naming the pathway collection.
    #' @param pathway Character scalar naming the gene set.
    #' @return A ggplot enrichment-curve object.
    make_enrichment_plot <- function(database, pathway) {
        # Pair one pathway's running enrichment curve with its result from the shared rank vector.
        result <- combined_gsea[combined_gsea$pathway == pathway, , drop = FALSE]
        fgsea::plotEnrichment(pathway_databases[[database]][[pathway]], gene_ranks) +
            ggplot2::labs(
                title = tools::toTitleCase(gsub("_", " ", pathway)),
                subtitle = sprintf("NES = %.3f, nominal P = %.3g", result$NES, result$pval)
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = ggplot2::element_text(hjust = 0.5)
            )
    }
    angiogenesis_plot <- make_enrichment_plot(
        "Hallmark",
        FIGURE7H_PATHWAYS[["Hallmark"]]
    )
    mtor_plot <- make_enrichment_plot(
        "Reactome",
        FIGURE7H_PATHWAYS[["Reactome"]]
    )
    figure7h <- patchwork::wrap_plots(angiogenesis_plot, mtor_plot, ncol = 2)
    ggplot2::ggsave(
        file.path(output_dir, "Figure7H_enrichment_curves.pdf"),
        figure7h,
        width = 12,
        height = 5
    )

    # Export pathway-resource coverage and run dimensions after all DEG/GSEA panels are complete.
    resource_summary <- merge(
        resource_manifest,
        data.frame(
            collection = GSEA_COLLECTIONS,
            pathways_loaded = vapply(
                pathway_databases[GSEA_COLLECTIONS],
                length,
                integer(1)
            )
        ),
        by = "collection",
        sort = FALSE
    )
    write_tsv(
        resource_summary,
        file.path(output_dir, "gsea_resource_summary.tsv")
    )
    gsea_summary <- data.frame(
        metric = c(
            "analysis_unit", "input_deg_genes", "rank_genes_after_filter",
            "filter_excluded_genes", "collections", "gsea_min_size",
            "gsea_max_size", "gsea_eps", "nominal_significant_pathways",
            "selected_display_pathways"
        ),
        value = c(
            "gene_rank_from_Bin100_spot_comparison", nrow(de_results),
            nrow(rank_table), nrow(de_results) - nrow(rank_table),
            paste(GSEA_COLLECTIONS, collapse = ";"), 3, 500, 0,
            nrow(significant), nrow(selected)
        ),
        stringsAsFactors = FALSE
    )
    write_tsv(gsea_summary, file.path(output_dir, "gsea_run_summary.tsv"))
}

if (sys.nframe() == 0L) {
    main()
}

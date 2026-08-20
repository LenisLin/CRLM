# Target: Supplementary Figure 16 CellChat outputs.
# Purpose: define cholangiocyte states and generate CellChat interaction summaries and panels.
# Inputs: edit the paths in CONFIG for the annotated scRNA R bundle, helpers, and output.
# Outputs: CellChat RDS, interaction and parameter TSVs, and Supplementary Figure 16 PDF panels.
# Ordered workflow: load reference, assign states, compute CellChat, summarize interactions, plot panels.

FIBROBLAST_PATHWAYS <- c("TGFb", "SPP1", "PDGF", "TNF", "COL1A1", "FN1")
ENDOTHELIAL_PATHWAYS <- c("VEGF", "SPP1", "TGFb", "ANGPT")
VASCULAR_PATHWAYS <- c("ANGPT", "FGF", "CXCL", "PDGF", "VEGF")
RECRUITMENT_PATHWAYS <- c("CXCL", "CCL")
SUPPRESSION_PATHWAYS <- c("TGFb", "SPP1", "LGALS9")

CHOLANGIOCYTE_GROUPS <- c("Chol_HypoxiaHigh", "Chol_HypoxiaLow")
VASCULAR_SUPPORTERS <- c(
    "Fibroblast", "Macro_CD163", "Macro_SPP1", "B_IGHM", "CD8T_GZMK",
    "CD8T_CXCL13", "MAIT_KLRB1"
)
REGULATORY_SENDERS <- c(
    "Fibroblast", "Macro_SPP1", "Macro_CD163", "cDC_CD1C"
)
LYMPHOID_RECEIVERS <- c("CD8T_GZMK", "CD8T_CXCL13", "Treg_FOXP3")

# Edit these paths for the Figure 7 analysis.
CONFIG <- list(
    sc_bundle = "/path/to/figure7_scrna_output/r_bundle",
    output_dir = "/path/to/figureS16_cellchat_output",
    scrna_functions = "src/scrna/cellchat.R"
)

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

#' Purpose: retain requested cell groups present in a CellChat result.
#'
#' @param requested Character vector of requested group names.
#' @param available Character vector of available group names.
#' @return Character vector of requested groups present in `available`.
existing_groups <- function(requested, available) {
    # Preserve requested plot order while limiting selections to modeled cell groups.
    requested[requested %in% available]
}

#' Purpose: retain requested signaling pathways present in a CellChat result.
#'
#' @param requested Character vector of requested signaling pathways.
#' @param cellchat A computed CellChat object.
#' @return Character vector of requested pathways present in `cellchat@netP$pathways`.
existing_pathways <- function(requested, cellchat) {
    # Preserve requested plot order while limiting selections to inferred pathways.
    requested[requested %in% cellchat@netP$pathways]
}

#' Purpose: generate and save one CellChat bubble-plot PDF.
#'
#' @param cellchat A computed CellChat object.
#' @param path Character scalar destination path for the PDF.
#' @param sources Character vector of source group names.
#' @param targets Character vector of target group names.
#' @param pathways Character vector of signaling pathways to display.
#' @param title Character scalar plot title.
#' @param width Numeric scalar PDF width in inches.
#' @param height Numeric scalar PDF height in inches.
#' @return The CellChat bubble-plot object after writing the PDF.
plot_bubble_pdf <- function(
    cellchat,
    path,
    sources,
    targets,
    pathways,
    title,
    width,
    height
) {
    # Plot inferred group-to-group communication only for pathways present in this CellChat result.
    plot <- CellChat::netVisual_bubble(
        cellchat,
        sources.use = sources,
        targets.use = targets,
        signaling = existing_pathways(pathways, cellchat),
        remove.isolate = FALSE,
        angle.x = 45,
        title.name = title
    )
    grDevices::pdf(path, width = width, height = height)
    print(plot)
    grDevices::dev.off()
    plot
}

#' Purpose: define the group and pathway selections for Supplementary Figure 16 panels.
#'
#' @param available_groups Character vector of groups retained by CellChat.
#' @return A named list of plot specifications.
build_plot_specs <- function(available_groups) {
    # Define each panel's sender, receiver, and pathway subset on the modeled group axis.
    list(
        FigureS16B_fibroblast = list(
            source = CHOLANGIOCYTE_GROUPS,
            target = "Fibroblast",
            pathway = FIBROBLAST_PATHWAYS
        ),
        FigureS16B_endothelial = list(
            source = CHOLANGIOCYTE_GROUPS,
            target = "Endothelial",
            pathway = ENDOTHELIAL_PATHWAYS
        ),
        FigureS16C_vascular = list(
            source = existing_groups(VASCULAR_SUPPORTERS, available_groups),
            target = "Endothelial",
            pathway = VASCULAR_PATHWAYS
        ),
        FigureS16D_recruitment = list(
            source = existing_groups(REGULATORY_SENDERS, available_groups),
            target = existing_groups(LYMPHOID_RECEIVERS, available_groups),
            pathway = RECRUITMENT_PATHWAYS
        ),
        FigureS16D_suppression = list(
            source = existing_groups(REGULATORY_SENDERS, available_groups),
            target = existing_groups(LYMPHOID_RECEIVERS, available_groups),
            pathway = SUPPRESSION_PATHWAYS
        )
    )
}

#' Purpose: retain interaction rows represented by Supplementary Figure 16 plot specifications.
#'
#' @param interactions Data frame of inferred CellChat interactions.
#' @param plot_specs Named list defining sources, targets, and pathways for each plot.
#' @return A data frame of plotted interactions with plot identifiers.
select_plotted_interactions <- function(interactions, plot_specs) {
    # Filter the interaction table with the same group and pathway selections used by each panel.
    selected <- lapply(names(plot_specs), function(plot_id) {
        specification <- plot_specs[[plot_id]]
        rows <- interactions[
            interactions$source %in% specification$source &
                interactions$target %in% specification$target &
                interactions$pathway_name %in% specification$pathway,
            ,
            drop = FALSE
        ]
        rows$plot_id <- rep(plot_id, nrow(rows))
        rows
    })
    # Retain one row per inferred interaction-panel membership, including shared interactions.
    do.call(rbind, selected)
}

#' Purpose: execute the configured Supplementary Figure 16 CellChat workflow.
#'
#' @param config List containing the scRNA bundle, helper, and output paths.
#' @return Invisibly returns the path to the final saved PDF after writing all outputs.
main <- function(config = CONFIG) {
    # Initialize the CellChat runtime and output locations before loading cell-level expression.
    output_dir <- config[["output_dir"]]
    panel_dir <- file.path(output_dir, "panels")
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

    suppressPackageStartupMessages({
        library(Seurat)
        library(CellChat)
        library(ggplot2)
        library(patchwork)
    })
    source(config[["scrna_functions"]])

    # Load the annotated scRNA reference and normalize expression on the existing cell identities.
    object <- read_anndata_to_seurat(config[["sc_bundle"]], verbose = TRUE)
    if (utils::packageVersion("Seurat") >= "5.0.0") {
        object <- SeuratObject::JoinLayers(object)
    }
    object <- Seurat::NormalizeData(object, verbose = FALSE)
    if (!all(CHOLANGIOCYTE_SCORE_GENES %in% rownames(object))) {
        stop(
            "The scRNA bundle must contain all six required cholangiocyte score genes",
            call. = FALSE
        )
    }
    # Assign cholangiocyte states from cell-level scores while retaining other subtype labels.
    object <- assign_cholangiocyte_states(
        object,
        genes = CHOLANGIOCYTE_SCORE_GENES,
        lower_quantile = 0.30,
        upper_quantile = 0.70,
        subtype_column = "Sub_type",
        score_column = "Metabolic_Score",
        output_column = "Sub_type_Revised"
    )

    # Export state membership and score cutoffs derived from the analyzed cholangiocyte cells.
    metadata <- object[[]]
    state_summary <- as.data.frame(table(
        original_subtype = metadata$Sub_type,
        revised_subtype = metadata$Sub_type_Revised
    ))
    state_summary <- state_summary[state_summary$Freq > 0, , drop = FALSE]
    write_tsv(
        state_summary,
        file.path(output_dir, "cellchat_state_summary.tsv")
    )
    cholangiocyte_scores <- metadata$Metabolic_Score[
        metadata$Sub_type == "Cholangiocyte"
    ]
    parameters <- data.frame(
        parameter = c(
            "input_subtypes", "cholangiocyte_score_genes", "lower_quantile",
            "upper_quantile", "lower_score_cutoff", "upper_score_cutoff",
            "database", "raw_use", "minimum_cells"
        ),
        value = c(
            paste(CELLCHAT_SUBTYPES, collapse = ";"),
            paste(CHOLANGIOCYTE_SCORE_GENES, collapse = ";"),
            0.30,
            0.70,
            unname(stats::quantile(cholangiocyte_scores, 0.30)),
            unname(stats::quantile(cholangiocyte_scores, 0.70)),
            "CellChatDB.human",
            TRUE,
            10
        ),
        stringsAsFactors = FALSE
    )
    write_tsv(parameters, file.path(output_dir, "cellchat_parameters.tsv"))

    # Infer communication between revised cell groups, requiring the configured minimum cells per group.
    cellchat <- compute_cellchat(
        object,
        group_by = "Sub_type_Revised",
        assay = "RNA",
        layer = "data",
        database = CellChat::CellChatDB.human,
        min_cells = 10
    )
    saveRDS(cellchat, file.path(output_dir, "cellchat_final.rds"))

    # Export all inferred interactions, then the subset represented by the panel specifications.
    interactions <- CellChat::subsetCommunication(cellchat)
    write_tsv(
        interactions,
        file.path(output_dir, "cellchat_interactions.tsv")
    )
    available_groups <- levels(cellchat@idents)
    plot_specs <- build_plot_specs(available_groups)
    plotted_interactions <- select_plotted_interactions(interactions, plot_specs)
    write_tsv(
        plotted_interactions,
        file.path(output_dir, "cellchat_plotted_interactions.tsv")
    )

    # Summarize incoming and outgoing interaction strength on the same modeled group axis.
    network_weights <- cellchat@net$weight
    activity <- data.frame(
        group = rep(rownames(network_weights), 2),
        direction = rep(c("Outgoing", "Incoming"), each = nrow(network_weights)),
        interaction_strength = c(
            rowSums(network_weights),
            colSums(network_weights)
        ),
        stringsAsFactors = FALSE
    )
    outgoing_order <- rownames(network_weights)[
        order(rowSums(network_weights), decreasing = TRUE)
    ]
    activity$group <- factor(activity$group, levels = outgoing_order)
    write_tsv(activity, file.path(output_dir, "cellchat_group_activity.tsv"))
    activity_plot <- ggplot2::ggplot(
        activity,
        ggplot2::aes(x = group, y = interaction_strength, fill = direction)
    ) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::scale_fill_manual(
            values = c("Outgoing" = "#E41A1C", "Incoming" = "#377EB8")
        ) +
        ggplot2::labs(
            title = "Global signaling landscape",
            x = NULL,
            y = "Total interaction strength"
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    ggplot2::ggsave(
        file.path(panel_dir, "01_Global_Activity_Landscape.pdf"),
        activity_plot,
        width = 10,
        height = 6
    )

    # Generate cholangiocyte-state and vascular panels from prespecified signaling selections.
    plot_bubble_pdf(
        cellchat,
        file.path(panel_dir, "02A_Chol_to_Fibroblast_Activation.pdf"),
        CHOLANGIOCYTE_GROUPS,
        "Fibroblast",
        FIBROBLAST_PATHWAYS,
        "Cholangiocyte states to fibroblast signaling",
        9,
        7
    )
    plot_bubble_pdf(
        cellchat,
        file.path(panel_dir, "02B_Chol_to_Endothelial_Interference.pdf"),
        CHOLANGIOCYTE_GROUPS,
        "Endothelial",
        ENDOTHELIAL_PATHWAYS,
        "Cholangiocyte states to endothelial signaling",
        9,
        7
    )
    plot_bubble_pdf(
        cellchat,
        file.path(panel_dir, "03_Ecosystem_Vascular_Maintenance.pdf"),
        existing_groups(VASCULAR_SUPPORTERS, available_groups),
        "Endothelial",
        VASCULAR_PATHWAYS,
        "Predicted vascular signaling to endothelial cells",
        11,
        7
    )

    # Pair recruitment and suppression views using identical regulatory sender and receiver sets.
    recruitment_plot <- CellChat::netVisual_bubble(
        cellchat,
        sources.use = existing_groups(REGULATORY_SENDERS, available_groups),
        targets.use = existing_groups(LYMPHOID_RECEIVERS, available_groups),
        signaling = existing_pathways(RECRUITMENT_PATHWAYS, cellchat),
        remove.isolate = FALSE,
        angle.x = 45,
        title.name = "A. Recruitment"
    )
    suppression_plot <- CellChat::netVisual_bubble(
        cellchat,
        sources.use = existing_groups(REGULATORY_SENDERS, available_groups),
        targets.use = existing_groups(LYMPHOID_RECEIVERS, available_groups),
        signaling = existing_pathways(SUPPRESSION_PATHWAYS, cellchat),
        remove.isolate = FALSE,
        angle.x = 45,
        title.name = "B. Suppression"
    )
    immune_plot <- patchwork::wrap_plots(
        recruitment_plot,
        suppression_plot,
        ncol = 2
    )
    ggplot2::ggsave(
        file.path(panel_dir, "04_Immune_Regulation_Combined.pdf"),
        immune_plot,
        width = 16,
        height = 8
    )
}

if (sys.nframe() == 0L) {
    main()
}

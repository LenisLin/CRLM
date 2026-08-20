# Purpose: Read the Figure 7 R bundle, assign cholangiocyte states, and infer
# CellChat communication networks.
# Callers: figures/Figure7/02_run_rctd.R and
# figures/Figure7/supplementary/FigureS16/01_run_plot_cellchat.R.
# Inputs: R-readable AnnData export bundle, Seurat metadata, subtype labels, and
# the CellChat human interaction database.
# Outputs: A Seurat object with revised cholangiocyte states and a CellChat
# object containing inferred communication probabilities and aggregated networks.
# Ordered use: Read the R bundle, assign cholangiocyte states, then compute
# CellChat using the revised subtype column.

CHOLANGIOCYTE_SCORE_GENES <- c("CA9", "SLC2A1", "LDHA", "PGK1", "ENO1", "ALDOA")
CELLCHAT_SUBTYPES <- c(
  "Cholangiocyte", "Endothelial", "Fibroblast", "Macro_CD163", "Macro_SPP1",
  "Malignant", "CD8T_GZMK", "B_IGHM", "CD8T_CXCL13", "MAIT_KLRB1",
  "cDC_CD1C", "Treg_FOXP3"
)

#' Read an exported AnnData MatrixMarket bundle into a Seurat object.
#'
#' @param data_path Path to the directory created by
#'   `export_anndata_for_r()`.
#' @param verbose Whether to report bundle-reading progress to the console.
#'
#' @return A Seurat object with the exported counts and cell metadata.
read_anndata_to_seurat <- function(data_path, verbose = TRUE) {
  if (verbose) cat("Reading expression matrix and metadata\n")

  # The Python export is cell-by-gene; Seurat requires genes by cells.
  expression_matrix <- Matrix::readMM(file.path(data_path, "expression_profile.mtx"))
  expression_matrix <- Matrix::t(expression_matrix)
  metadata <- read.csv(file.path(data_path, "metadata.csv"), row.names = 1)
  gene_names <- read.csv(
    file.path(data_path, "row_names.csv"), header = FALSE
  )$V1
  cell_names <- read.csv(
    file.path(data_path, "column_names.csv"), header = FALSE
  )$V1

  # Apply exported axis names before filtering unnamed genes and attaching metadata.
  rownames(expression_matrix) <- gene_names
  colnames(expression_matrix) <- cell_names
  expression_matrix <- expression_matrix[rownames(expression_matrix) != "", , drop = FALSE]
  object <- Seurat::CreateSeuratObject(
    counts = expression_matrix,
    meta.data = metadata,
    project = "FDZS3_scRNA"
  )
  object[rownames(object) != "", ]
}


#' Retain cells whose subtype is included in the CellChat allow-list.
#'
#' @param object Seurat object containing subtype annotations.
#' @param subtype_column Metadata column holding subtype labels.
#' @param allowed_subtypes Character vector of subtype labels retained for
#'   CellChat analysis.
#'
#' @return A Seurat object restricted to the allowed subtypes.
#'
#' @details Stops when no cells match the allow-list.
filter_cellchat_subtypes <- function(
  object,
  subtype_column = "Sub_type",
  allowed_subtypes = CELLCHAT_SUBTYPES
) {
  # Resolve the allow-list against metadata while preserving Seurat cell order.
  subtype <- as.character(object[[subtype_column, drop = TRUE]])
  cells <- colnames(object)[subtype %in% allowed_subtypes]
  if (!length(cells)) {
    stop("No cells match the CellChat subtype allow-list")
  }
  subset(object, cells = cells)
}


#' Assign hypoxia-associated cholangiocyte states from expression scores.
#'
#' @param object Seurat object containing normalized RNA expression and subtype
#'   annotations.
#' @param genes Gene symbols used to compute the metabolic score.
#' @param lower_quantile Quantile defining the low-score cholangiocyte state.
#' @param upper_quantile Quantile defining the high-score cholangiocyte state.
#' @param subtype_column Metadata column holding input subtype labels.
#' @param score_column Metadata column receiving per-cell metabolic scores.
#' @param output_column Metadata column receiving revised subtype labels.
#' @param cholangiocyte_label Input subtype label for cholangiocytes.
#' @param high_label Revised label for high-score cholangiocytes.
#' @param low_label Revised label for low-score cholangiocytes.
#'
#' @return A subtype-filtered Seurat object with metabolic scores and revised
#'   subtype labels.
assign_cholangiocyte_states <- function(
  object,
  genes = CHOLANGIOCYTE_SCORE_GENES,
  lower_quantile = 0.30,
  upper_quantile = 0.70,
  subtype_column = "Sub_type",
  score_column = "Metabolic_Score",
  output_column = "Sub_type_Revised",
  cholangiocyte_label = "Cholangiocyte",
  high_label = "Chol_HypoxiaHigh",
  low_label = "Chol_HypoxiaLow"
) {
  # Restrict scoring to the communication-analysis populations and available genes.
  object <- filter_cellchat_subtypes(
    object,
    subtype_column = subtype_column
  )
  valid_genes <- genes[genes %in% rownames(object)]
  gene_data <- Seurat::FetchData(object, vars = valid_genes, layer = "data")
  scores <- rowMeans(gene_data)
  object[[score_column]] <- scores

  # Derive the state thresholds only from cholangiocyte scores.
  subtype <- object[[subtype_column, drop = TRUE]]
  cholangiocyte_cells <- colnames(object)[subtype == cholangiocyte_label]
  cholangiocyte_scores <- scores[cholangiocyte_cells]
  high_cutoff <- quantile(cholangiocyte_scores, upper_quantile)
  low_cutoff <- quantile(cholangiocyte_scores, lower_quantile)

  # Relabel the lower and upper tails while retaining all other subtype labels.
  revised_subtype <- as.character(subtype)
  high_cells <- names(cholangiocyte_scores)[cholangiocyte_scores >= high_cutoff]
  low_cells <- names(cholangiocyte_scores)[cholangiocyte_scores <= low_cutoff]
  revised_subtype[match(high_cells, colnames(object))] <- high_label
  revised_subtype[match(low_cells, colnames(object))] <- low_label
  object[[output_column]] <- revised_subtype
  object
}


#' Infer and aggregate CellChat communication networks.
#'
#' @param object Seurat object containing expression data and group labels.
#' @param group_by Metadata column defining CellChat groups.
#' @param assay Seurat assay supplying expression data.
#' @param layer Assay layer supplied to `GetAssayData()`.
#' @param database CellChat interaction database assigned to the analysis.
#' @param min_cells Minimum cells per group retained for communication filtering.
#'
#' @return A CellChat object with inferred communication probabilities,
#' pathways, and aggregated network summaries.
compute_cellchat <- function(
  object,
  group_by = "Sub_type_Revised",
  assay = "RNA",
  layer = "data",
  database = CellChat::CellChatDB.human,
  min_cells = 10
) {
  # Supply expression columns and metadata rows from the same retained cells.
  object <- filter_cellchat_subtypes(object)
  data_use <- Seurat::GetAssayData(object, assay = assay, layer = layer)
  metadata <- object@meta.data
  cellchat <- CellChat::createCellChat(
    object = data_use,
    meta = metadata,
    group.by = group_by
  )
  # Restrict the database to expressed interactions before estimating cell-group
  # probabilities and aggregating them into pathway-level networks.
  cellchat@DB <- database
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  CellChat::aggregateNet(cellchat)
}

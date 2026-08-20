#!/usr/bin/env Rscript

# Target: Figure 4 methods-aligned table generation.
# Purpose: Derive peritumoral cellular-neighborhood, cholangiocyte-associated
# PTME, and ISR tables using the shared IMC implementations.
# Inputs: An SCE RDS with the required Figure 4 cell metadata; src/imc source
# files; and the fixed parameters below.
# Outputs: Parameter, cellular-neighborhood, cholangiocyte-mask, PTME,
# ROI-level ISR, and patient-level ISR TSV files in OUTPUT_DIRECTORY.
# Ordered workflow:
# 1. Configure SCE_INPUT, OUTPUT_DIRECTORY, and SCRIPT_DIRECTORY.
# 2. Verify required R packages and src/imc interfaces.
# 3. Load the configured SCE and validate its required cell metadata.
# 4. Restrict the SCE to peritumoral cells and derive patient treatment labels.
# 5. Calculate cellular-neighborhood composition and define cholangiocyte masks.
# 6. Delineate and classify cholangiocyte-associated PTMEs.
# 7. Summarize ROI-level and patient-level ISR.
# 8. Write all Figure 4 tables and the fixed-parameter record.

# Configure these paths before direct execution from the repository root.
SCE_INPUT <- "FDZS1_IMC_processed.rds"
OUTPUT_DIRECTORY <- file.path("outputs", "Figure4")
SCRIPT_DIRECTORY <- file.path("figures", "Figure4")

CN_NEIGHBORS <- 20L
PTME_MIN_CHOLANGIOCYTES <- 5L
PTME_EXPANSION_UM <- 22
NICHE_RATIO_THRESHOLD <- 1
ISR_LOWER <- 0.1
ISR_UPPER <- 10
IMMUNE_SUBTYPES <- c("B", "CD8T")
STROMAL_SUBTYPES <- c("SC_COLLAGEN", "Macro_CD163")

#' Purpose: Stop if any required package is unavailable.
#'
#' @param packages Character vector of package names.
#' @return Invisibly returns `NULL`; otherwise stops with missing package names.
require_packages <- function(packages) {
    # Resolve all namespaces before the workflow reads data or creates outputs.
    missing <- packages[!vapply(
        packages, requireNamespace, logical(1L), quietly = TRUE
    )]
    if (length(missing) > 0L) {
        stop("Missing required packages: ", paste(missing, collapse = ", "))
    }
}

#' Purpose: Source the shared IMC implementations and verify their interface.
#'
#' @param script_directory Directory containing this Figure 4 script.
#' @return Invisibly returns `NULL`; otherwise stops when a required file or
#'   function is absent.
load_shared_imc <- function(script_directory) {
    # Map the figure entry point to the shared IMC implementations it orchestrates.
    shared_directory <- normalizePath(
        file.path(script_directory, "..", "..", "src", "imc"),
        mustWork = FALSE
    )
    shared_files <- file.path(shared_directory, c(
        "coordinates.R", "cellular_neighborhoods.R", "ptme.R"
    ))
    missing <- shared_files[!file.exists(shared_files)]
    if (length(missing) > 0L) {
        stop("Missing required src/imc file(s): ", paste(missing, collapse = ", "))
    }
    invisible(lapply(shared_files, source))

    # Keep the expected shared-function interface explicit at the call boundary.
    required_functions <- c(
        "calculate_cn_composition", "delineate_cholangiocyte_ptmes",
        "classify_ptmes", "summarize_imc_isr"
    )
    unavailable <- required_functions[!vapply(
        required_functions, exists, logical(1L), mode = "function"
    )]
    if (length(unavailable) > 0L) {
        stop("Required src/imc interface is missing: ", paste(unavailable, collapse = ", "))
    }
}

#' Purpose: Verify that the SCE cell metadata satisfies the Figure 4 contract.
#'
#' @param cell_metadata Data frame converted from `SingleCellExperiment::colData()`.
#' @return Invisibly returns `NULL`; otherwise stops when required fields or
#'   treatment assignments are invalid.
validate_cell_metadata <- function(cell_metadata) {
    # Establish the cell, ROI, patient, tissue, and treatment fields used downstream.
    required <- c(
        "CellID", "PID", "ID", "Tissue", "SubType", "Position",
        "MajorType2", "Treatment"
    )
    missing <- setdiff(required, names(cell_metadata))
    if (length(missing) > 0L) {
        stop("SCE colData is missing: ", paste(missing, collapse = ", "))
    }
    if (anyDuplicated(cell_metadata$CellID)) {
        stop("CellID must be unique")
    }
    complete_fields <- c(
        "CellID", "PID", "ID", "Tissue", "SubType", "Position",
        "MajorType2", "Treatment"
    )
    if (any(!stats::complete.cases(cell_metadata[, complete_fields, drop = FALSE]))) {
        stop("Required Figure 4 colData fields must be complete")
    }
    if (!all(as.character(cell_metadata$Treatment) %in% c("Chemo", "Combo"))) {
        stop("Treatment must be the explicit clinical field Chemo or Combo")
    }

    # Each patient must map to exactly one explicit treatment assignment.
    patient_treatment <- unique(cell_metadata[, c("PID", "Treatment"), drop = FALSE])
    if (anyDuplicated(patient_treatment$PID)) {
        stop("Treatment must be constant within each PID")
    }
}

#' Purpose: Produce one explicit treatment assignment per patient.
#'
#' @param cell_metadata Data frame containing `PID` and `Treatment` columns.
#' @return Data frame with character `PID` and `Treatment` columns and one row
#'   per patient.
patient_treatment_table <- function(cell_metadata) {
    # Collapse repeated cell metadata to the one-row-per-patient treatment mapping.
    treatment <- unique(cell_metadata[, c("PID", "Treatment"), drop = FALSE])
    treatment$PID <- as.character(treatment$PID)
    treatment$Treatment <- as.character(treatment$Treatment)
    rownames(treatment) <- NULL
    treatment
}

#' Purpose: Add the explicit patient treatment assignment to a derived table.
#'
#' @param data Derived data frame containing a `PID` column.
#' @param patient_treatment One-row-per-patient treatment data frame.
#' @return The input data frame with a character `Treatment` column; otherwise
#'   stops if a derived row lacks a treatment assignment.
attach_treatment <- function(data, patient_treatment) {
    # Match by PID while preserving the row order of each derived analysis table.
    index <- match(as.character(data$PID), patient_treatment$PID)
    if (anyNA(index)) {
        stop("A derived row has no explicit Treatment value")
    }
    data$Treatment <- patient_treatment$Treatment[index]
    data
}

#' Purpose: Record the fixed scientific parameters used for Figure 4 tables.
#'
#' This function has no parameters.
#' @return Data frame of parameter names, configured values, and units.
methods_parameters <- function() {
    # Encode the fixed analysis choices in the same tabular form as other outputs.
    data.frame(
        parameter = c(
            "cn_neighbors", "cn_weighting", "cholangiocyte_mask",
            "ptme_grouping", "ptme_min_cholangiocytes", "ptme_expansion",
            "niche_numerator_subtypes", "niche_denominator_subtypes",
            "niche_threshold", "niche_equality", "isr_definition",
            "isr_roster", "isr_clipping", "treatment_source"
        ),
        value = c(
            as.character(CN_NEIGHBORS), "unweighted",
            "Tissue == PT and MajorType2 == EC",
            "all-PT-cell within-ROI Delaunay graph; induced cholangiocyte components",
            as.character(PTME_MIN_CHOLANGIOCYTES),
            as.character(PTME_EXPANSION_UM),
            paste(STROMAL_SUBTYPES, collapse = ";"),
            paste(IMMUNE_SUBTYPES, collapse = ";"),
            as.character(NICHE_RATIO_THRESHOLD), "Unclassified",
            "(N_PIR + 1) / (N_PSM + 1)",
            "all eligible PT ROIs and patients, including zero-structure units",
            paste0("[", ISR_LOWER, ", ", ISR_UPPER, "]"),
            "colData.Treatment"
        ),
        unit = c(
            "cells", NA, NA, NA, "cells", "micrometres", NA, NA,
            "ratio", NA, "structures", NA, "ratio", NA
        ),
        stringsAsFactors = FALSE
    )
}

#' Purpose: Write a data frame as a tab-separated table with explicit missing values.
#'
#' @param data Data frame to export.
#' @param path Destination TSV file path.
#' @return Invisibly returns `NULL` after writing the table.
write_tsv <- function(data, path) {
    # Use a stable, unquoted representation shared by every exported Figure 4 table.
    utils::write.table(
        data, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
    )
}

#' Purpose: Generate all Figure 4 methods-aligned tables from an IMC SCE.
#'
#' @param sce `SingleCellExperiment` containing the Figure 4 IMC data.
#' @param output_directory Directory for the generated TSV files.
#' @return Invisibly returns a named list of cellular-neighborhood, mask, PTME,
#'   ROI-level ISR, and patient-level ISR data frames.
run_figure4_methods <- function(sce, output_directory) {
    # Map authoritative SCE cell metadata into the peritumoral analysis population.
    cell_metadata <- as.data.frame(SingleCellExperiment::colData(sce))
    validate_cell_metadata(cell_metadata)

    pt_index <- as.character(cell_metadata$Tissue) == "PT"
    if (!any(pt_index)) {
        stop("The SCE contains no PT cells")
    }
    pt_sce <- sce[, pt_index]
    pt_metadata <- as.data.frame(SingleCellExperiment::colData(pt_sce))
    patient_treatment <- patient_treatment_table(pt_metadata)

    # Derive neighborhood composition before constructing the cholangiocyte mask.
    cn_composition <- calculate_cn_composition(pt_sce, k = CN_NEIGHBORS)
    cn_composition <- attach_treatment(cn_composition, patient_treatment)

    epcam_positive_ec_mask <- as.character(pt_metadata$MajorType2) == "EC"
    cholangiocyte_mask <- data.frame(
        CellID = pt_metadata$CellID,
        PID = pt_metadata$PID,
        ID = pt_metadata$ID,
        MajorType2 = as.character(pt_metadata$MajorType2),
        epcam_positive_ec_cholangiocyte = epcam_positive_ec_mask,
        stringsAsFactors = FALSE
    )
    cholangiocyte_mask <- attach_treatment(cholangiocyte_mask, patient_treatment)

    # Delineate spatial memberships, then classify each retained PTME from cell content.
    membership <- delineate_cholangiocyte_ptmes(
        pt_sce,
        cholangiocyte_mask = epcam_positive_ec_mask,
        min_cholangiocytes = PTME_MIN_CHOLANGIOCYTES,
        expansion_um = PTME_EXPANSION_UM
    )
    membership <- attach_treatment(membership, patient_treatment)

    ptmes <- classify_ptmes(
        membership,
        cell_metadata = pt_metadata,
        threshold = NICHE_RATIO_THRESHOLD,
        immune_subtypes = IMMUNE_SUBTYPES,
        stromal_subtypes = STROMAL_SUBTYPES
    )
    ptmes <- attach_treatment(ptmes, patient_treatment)

    # Preserve zero-structure ROIs and patients through explicit summary rosters.
    roi_roster <- unique(pt_metadata[, c("PID", "ID"), drop = FALSE])
    patient_roster <- unique(pt_metadata[, "PID", drop = FALSE])
    roi_isr <- summarize_imc_isr(
        ptmes, roster = roi_roster, level = "ROI",
        lower = ISR_LOWER, upper = ISR_UPPER
    )
    roi_isr <- attach_treatment(roi_isr, patient_treatment)
    patient_isr <- summarize_imc_isr(
        ptmes, roster = patient_roster, level = "patient",
        lower = ISR_LOWER, upper = ISR_UPPER
    )
    patient_isr <- attach_treatment(patient_isr, patient_treatment)

    # Export parameters and every intermediate table needed to reproduce the figure inputs.
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    write_tsv(methods_parameters(), file.path(
        output_directory, "figure4_methods_parameters.tsv"
    ))
    write_tsv(cn_composition, file.path(
        output_directory, "figure4_pt_cn_composition.tsv"
    ))
    write_tsv(cholangiocyte_mask, file.path(
        output_directory, "figure4_pt_cholangiocyte_mask.tsv"
    ))
    write_tsv(membership, file.path(
        output_directory, "figure4_ptme_membership.tsv"
    ))
    write_tsv(ptmes, file.path(
        output_directory, "figure4_ptme_classification.tsv"
    ))
    write_tsv(roi_isr, file.path(
        output_directory, "figure4_roi_isr.tsv"
    ))
    write_tsv(patient_isr, file.path(
        output_directory, "figure4_patient_isr.tsv"
    ))

    invisible(list(
        cn_composition = cn_composition,
        cholangiocyte_mask = cholangiocyte_mask,
        membership = membership,
        ptmes = ptmes,
        roi_isr = roi_isr,
        patient_isr = patient_isr
    ))
}

#' Purpose: Run the configured Figure 4 table-generation workflow.
#'
#' This function reads the editable in-script path configuration and has no parameters.
#' @return Invisibly returns the named output list from `run_figure4_methods()`;
#'   otherwise stops when configured inputs or dependencies are unavailable.
main <- function() {
    # Resolve dependencies and shared code before loading the configured SCE.
    require_packages(c(
        "SingleCellExperiment", "SummarizedExperiment", "FNN", "igraph",
        "geometry"
    ))
    load_shared_imc(SCRIPT_DIRECTORY)
    if (!file.exists(SCE_INPUT)) {
        stop("Configured SCE input does not exist: ", SCE_INPUT)
    }
    sce <- readRDS(SCE_INPUT)
    run_figure4_methods(sce, OUTPUT_DIRECTORY)
}

if (sys.nframe() == 0L) {
    main()
}

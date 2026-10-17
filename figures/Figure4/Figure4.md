Figure 4. Associations between cholangiocyte-associated peritumor niches
and treatment-stratified RFS in CRLM.
================

## Purpose

Figure 4 examines peritumoral cellular neighborhoods, epithelial
identity, and cholangiocyte-associated peritumor microenvironments
(PTMEs) in colorectal liver metastases (CRLM). The workflow derives
immune-rich (PIR) and stromal-rich (PSM) niche classes and the
immune-to-stromal ratio (ISR) used to characterize patients within
postoperative treatment groups.

## Inputs and analytical units

The [Figure 4 entry points](README.md) use the processed FDZS-1 IMC
`SingleCellExperiment` described in the [public data
contract](../../docs/data_contract.md). Required cell metadata are
unique `CellID` and complete `PID`, `ID`, `Tissue`, `SubType`,
`Position`, `MajorType2`, and `Treatment`. `PID` identifies patients,
`ID` identifies nested ROIs, and `Position` contains coordinates in
micrometres. The clinical `Treatment` field is constant within each
patient: `Chemo` denotes FOLFOX alone and `Combo` denotes FOLFOX plus
targeted therapy.

In the code below, spatial analyses use PT cells. Epithelial comparisons
use TC, IM, and PT, the four configured epithelial subtypes (`TC_CAIX`,
`TC_EpCAM`, `TC_Ki67`, `TC_VEGF`), and PRPS1, FASN, GLUT1, HK2, Ki67,
VEGF, and CAIX expression from the existing assay. Cells are the units
for neighborhood composition; structures are the units for PTME
classification and size; ROI and patient summaries retain their
respective grouping levels. Marker comparisons use one patient-tissue
mean obtained by equally averaging ROI means. The ROI export contains
one row per patient–ROI–tissue group.

## Peritumoral neighborhoods, PTME classification, and ISR

The [table-generation entry point](01_methods_aligned_orchestration.R)
calculates unweighted subtype proportions among each PT cell’s 20
nearest neighbors within its patient and ROI. Shared implementations are
linked in [coordinates.R](../../src/imc/coordinates.R),
[cellular_neighborhoods.R](../../src/imc/cellular_neighborhoods.R), and
[ptme.R](../../src/imc/ptme.R).

### Published peritumoral niches and treatment-stratified RFS

![Complete published Figure 4 showing peritumoral neighborhoods,
epithelial identity, PTMEs, ISR, and
RFS](report_assets/figure4_publication.png)

**Figure 4. Associations between cholangiocyte-associated peritumor
niches and treatment-stratified RFS in CRLM.** (A) Peritumoral CN
composition and observed-to-expected enrichment. (B) ROI-level CN1/CN4
fractions by recurrence status and regimen, with BH-adjusted P values.
(C) H&E and EpCAM/CK20 multiplex immunofluorescence. (D) ROI-mean
epithelial CAIX and Ki67 distributions across TC, IM, and PT, with the
published nominal Welch-test annotations. (E) PTME delineation and ISR
schematic. (F) Cell fractions within individual PTMEs by recurrence
status and treatment. (G) Representative PIR/PSM IMC, H&E, and cell
maps. (H) ROI-level ISR above and patient-level ISR below. (I)
Patient-level RFS curves for four ISR/treatment groups, with 95%
confidence bands, numbers at risk, and nominal log-rank P values. (J)
Niche composition and ISR concept. Imaging scale bars are 50 µm;
blue/red denote no recurrence/recurrence in the comparison panels.

[Download Figure 4 as PDF](report_assets/figure4_publication.pdf)

The published cetuximab-stratum ROI comparisons show higher CN1
fractions with recurrence and higher CN4 fractions without recurrence
(BH-adjusted P = 0.0049 and 0.0156). In the combination-therapy group,
the published ISR comparisons give nominal P = 0.013 at ROI level and P
= 0.045 at patient level. Additional peritumoral CN comparisons are
described in the [Figure 3 report](../Figure3/Figure3.md).

### Configuration and input validation

The configuration fixes the spatial and classification parameters. The
functions check package availability, shared interfaces, complete cell
metadata, and one explicit treatment assignment per patient.

<!-- source: 01_methods_aligned_orchestration.R; lines: 1-80 -->

``` r
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
```

<!-- source: 01_methods_aligned_orchestration.R; lines: 81-146 -->

``` r
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
```

### Spatial definitions and structure counts

PT `MajorType2 == "EC"` cells define the cholangiocyte mask. The shared
PTME function builds an all-PT-cell Delaunay graph within each ROI,
retains its induced cholangiocyte components containing at least five
cholangiocytes, and expands each convex hull by 22 µm. The niche ratio
is `(SC_COLLAGEN + Macro_CD163) / (B + CD8T)`: values below 1 define PIR
and values above 1 define PSM. Equality and an empty numerator and
denominator yield an unclassified structure; a positive numerator with a
zero denominator yields PSM.

ISR is `(N_PIR + 1) / (N_PSM + 1)`, clipped to \[0.1, 10\]. Explicit PT
ROI and patient rosters retain units with zero PIR and PSM counts, for
which ISR equals 1. The driver writes cell-level compositions, the
cholangiocyte mask, PTME membership and classification, and ROI- and
patient-level ISR tables.

<!-- source: 01_methods_aligned_orchestration.R; lines: 147-194 -->

``` r
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
```

<!-- source: 01_methods_aligned_orchestration.R; lines: 195-294 -->

``` r
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
```

<!-- source: 01_methods_aligned_orchestration.R; lines: 295-316 -->

``` r
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
```

## ROI-level cell-composition export

The [ROI export](02_export_roi_cell_composition.R) groups the public
cell annotations by patient, ROI, and tissue, then divides subtype
counts by the total cells in each group. Its configured output contains
311 groups and a fraction column for every observed subtype. The
exported matrix supports ROI-level composition summaries.

### Annotation fields and fraction matrix

<!-- source: 02_export_roi_cell_composition.R; lines: 1-49 -->

``` r
#!/usr/bin/env Rscript

# Figure 4 public ROI-level cell-composition export.
#
# Purpose: Aggregate the public FDZS-1 IMC cell annotations to one wide row per
# patient, ROI, and tissue region for distribution through Zenodo.
# Inputs: The public FDZS-1 SingleCellExperiment RDS selected in EXPORT_CONFIG.
# Output: A gzip-compressed TSV containing ROI identifiers, total cell counts,
# and SubType fractions.
# Runtime: Host Conda environment Spatial with R 4.2.2.

RUN_ROI_EXPORT <- FALSE

EXPORT_CONFIG <- list(
  sce_input = "path/to/FDZS1_IMC_processed.rds",
  output = "path/to/FDZS1_ROI_cell_composition.tsv.gz",
  expected_rois = 311L
)

required_packages <- c("SingleCellExperiment", "SummarizedExperiment")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

#' Read the public FDZS-1 IMC object and select cell-level annotation fields.
#'
#' @param path Path to the public SingleCellExperiment RDS.
#' @return A data frame containing patient, ROI, tissue, and cell-subtype fields.
read_cell_annotations <- function(path) {
  # The export reads colData only; expression assays are not transformed.
  if (!file.exists(path)) stop("Input does not exist: ", path)
  sce <- readRDS(path)
  required <- c("PID", "ID", "Tissue", "SubType")
  missing <- setdiff(required, names(SummarizedExperiment::colData(sce)))
  if (length(missing) > 0L) {
    stop("FDZS-1 object is missing fields: ", paste(missing, collapse = ", "))
  }
  annotations <- as.data.frame(
    SummarizedExperiment::colData(sce)[, required],
    stringsAsFactors = FALSE
  )
  names(annotations) <- c("patient_id", "roi_id", "tissue", "cell_subtype")
  if (anyNA(annotations) || any(!nzchar(as.matrix(annotations)))) {
    stop("ROI composition fields must be complete and non-empty.")
  }
  annotations
}
```

<!-- source: 02_export_roi_cell_composition.R; lines: 50-96 -->

``` r
#' Aggregate cell annotations to a wide ROI-level fraction matrix.
#'
#' @param annotations Cell-level annotation data frame from read_cell_annotations().
#' @param expected_rois Expected number of unique patient-ROI-tissue groups.
#' @return One row per ROI with total_cells and one fraction column per SubType.
build_roi_matrix <- function(annotations, expected_rois) {
  # A composite key preserves the nesting of ROI identifiers within patients.
  group_key <- interaction(
    annotations$patient_id,
    annotations$roi_id,
    annotations$tissue,
    drop = TRUE,
    lex.order = TRUE
  )
  subtype <- factor(
    annotations$cell_subtype,
    levels = sort(unique(annotations$cell_subtype))
  )
  counts <- unclass(table(group_key, subtype))
  total_cells <- rowSums(counts)
  fractions <- sweep(counts, 1L, total_cells, "/")

  # Recover the explicit group fields in the same order as the table rows.
  group_rows <- unique(data.frame(
    group_key = as.character(group_key),
    patient_id = annotations$patient_id,
    roi_id = annotations$roi_id,
    tissue = annotations$tissue,
    stringsAsFactors = FALSE
  ))
  group_rows <- group_rows[match(rownames(counts), group_rows$group_key), , drop = FALSE]
  if (anyNA(group_rows$group_key) || nrow(group_rows) != expected_rois) {
    stop("Expected ", expected_rois, " ROI groups; found ", nrow(group_rows), ".")
  }
  if (any(abs(rowSums(fractions) - 1) > 1e-12)) {
    stop("ROI cell fractions must sum to one.")
  }

  colnames(fractions) <- paste0("cell_fraction__", colnames(fractions))
  result <- cbind(
    group_rows[, c("patient_id", "roi_id", "tissue"), drop = FALSE],
    total_cells = as.integer(total_cells),
    as.data.frame(fractions, check.names = FALSE)
  )
  result[order(result$patient_id, result$roi_id, result$tissue), , drop = FALSE]
}
```

### Compressed output

The writer saves the matrix as a gzip-compressed TSV; the driver calls
the annotation reader and aggregation function in sequence.

<!-- source: 02_export_roi_cell_composition.R; lines: 97-126 -->

``` r
#' Write the public ROI-level matrix as a compressed TSV.
#'
#' @param data ROI-level cell-composition data frame.
#' @param path Destination .tsv.gz path.
#' @return An invisible NULL after the table is written.
write_roi_matrix <- function(data, path) {
  # Use gzip compression while retaining a directly readable tabular format.
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.table(data, connection, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(NULL)
}

#' Run the public FDZS-1 ROI composition export.
#'
#' @param config List containing sce_input, output, and expected_rois.
#' @return An invisible ROI-level data frame after writing the configured output.
run_roi_export <- function(config) {
  # Execute the read, patient-ROI aggregation, and compressed export in order.
  annotations <- read_cell_annotations(config$sce_input)
  result <- build_roi_matrix(annotations, config$expected_rois)
  write_roi_matrix(result, config$output)
  message(config$output)
  invisible(result)
}

if (isTRUE(RUN_ROI_EXPORT)) {
  run_roi_export(EXPORT_CONFIG)
}
```

## Epithelial identity: Figure 4D and Supplementary Figure 9A–B

The [epithelial-identity entry point](03_plot_epithelial_identity.R)
summarizes the regional distribution of epithelial cells and compares
marker expression across TC, IM, and PT.

### Epithelial population and regional distribution

The population mask combines `MajorType2 == "EC"`, the four configured
epithelial subtypes, and the three tissue regions. Regional fractions
divide each tissue’s epithelial-cell count by the total epithelial-cell
count across these regions.

<!-- source: 03_plot_epithelial_identity.R; lines: 1-88 -->

``` r
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
```

<!-- source: 03_plot_epithelial_identity.R; lines: 89-109 -->

``` r
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
```

### ROI and patient marker summaries

The code first calculates an epithelial-cell expression mean for each
patient–ROI–tissue–marker group. It then averages the ROI means within
each patient and tissue, giving each ROI equal weight.

<!-- source: 03_plot_epithelial_identity.R; lines: 110-159 -->

``` r
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
```

<!-- source: 03_plot_epithelial_identity.R; lines: 160-184 -->

``` r
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
```

### Tissue comparisons and marker panels

For each marker, the public entry point applies nominal two-sided
Wilcoxon rank-sum tests to the three tissue pairs, with `exact = FALSE`.
The code draws patient-tissue means and individual points. Its Figure 4D
output uses CAIX and Ki67; its Supplementary Figure 9B output uses VEGF,
FASN, GLUT1, HK2, and PRPS1. Stars encode nominal P-value thresholds of
0.05, 0.01, 0.001, and 0.0001.

<!-- source: 03_plot_epithelial_identity.R; lines: 185-229 -->

``` r
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
```

<!-- source: 03_plot_epithelial_identity.R; lines: 230-283 -->

``` r
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
```

### Distribution panel and output generation

The regional-distribution function draws the epithelial-cell fractions
as a pie chart. The driver saves the distribution, ROI means,
patient-tissue means, tissue-pair tests, and the three panel PDFs.

<!-- source: 03_plot_epithelial_identity.R; lines: 284-311 -->

``` r
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
```

<!-- source: 03_plot_epithelial_identity.R; lines: 312-368 -->

``` r
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
```

### Published quantitative PTME features

![Complete published Figure S9 showing epithelial distribution and
marker expression, PTME size and composition, ROI summaries, and niche
protein differences](report_assets/figureS9_publication.png)

**Figure S9. Quantitative features of cholangiocyte-centered peritumor
microenvironments.** (A) Pooled EpCAM-positive cell distribution: TC
66%, IM 32%, and PT 2%. (B) ROI-mean epithelial VEGF, FASN, GLUT1, HK2,
and PRPS1 distributions with the published nominal Welch-test
annotations. (C) Number of cells per PTME. (D) Cell-type fractions
within PTMEs and the overall PT region. (E) PTME counts per ROI and mean
PTME size, overall and by treatment. (F) ROI-mean marker expression
within PTMEs, by recurrence status and treatment. (G) Cell-type-specific
protein differences between PSM and PIR niches; color shows log2 fold
change and size shows −log10(BH-adjusted P). Comparison annotations
retain the publication values.

The pooled epithelial distribution places only 2% of epithelial cells in
PT. The publication marker panels show higher CAIX and lower Ki67 in
peritumoral epithelial cells than in tumor-core epithelial cells (Figure
4D), alongside lower VEGF, FASN, GLUT1, HK2, and PRPS1 (Figure S9B).

[Download the standalone epithelial-distribution
PDF](report_assets/supplementary9a_epithelial_distribution.pdf). Its CT
and TAT labels correspond to TC and PT.

## PTME structure size: Supplementary Figure 9C

The [PTME plot consumer](04_plot_ptme_isr_panels.R) reads
`figure4_ptme_classification.tsv` from the table-generation workflow.
Each row represents one PTME, identified by `ptme_id`, patient, and ROI;
`n_cells` measures structure size as the number of member cells.

### Structure table and size histogram

The reader requires complete, unique PTME rows and positive integer cell
counts. The histogram uses 30 bins and displays the number of PTMEs at
each structure size.

<!-- source: 04_plot_ptme_isr_panels.R; lines: 1-57 -->

``` r
#!/usr/bin/env Rscript

# Figure 4 PTME/ISR plot consumer.
#
# Purpose: Render only plot panels whose full input and statistical contracts
# are present in the TSV outputs of 01_methods_aligned_orchestration.R.
# Current coverage: Supplementary Figure 9C.

RUN_FIGURE4_PLOT_CONSUMERS <- FALSE

TABLE_DIRECTORY <- file.path("outputs", "Figure4")
OUTPUT_DIRECTORY <- file.path("outputs", "Figure4", "plots")

#' Purpose: Stop if a plotting dependency is unavailable.
#'
#' @param packages Character vector of package names.
#' @return Invisibly returns NULL; otherwise stops with missing package names.
require_packages <- function(packages) {
    missing <- packages[!vapply(
        packages, requireNamespace, logical(1L), quietly = TRUE
    )]
    if (length(missing) > 0L) {
        stop("Missing required packages: ", paste(missing, collapse = ", "))
    }
}

#' Purpose: Read and validate the per-PTME table required for Supplementary Figure 9C.
#'
#' @param table_directory Directory written by 01_methods_aligned_orchestration.R.
#' @return One row per classified PTME with a positive integer n_cells value.
read_ptme_classification <- function(table_directory) {
    path <- file.path(table_directory, "figure4_ptme_classification.tsv")
    if (!file.exists(path)) {
        stop("Missing Figure 4 PTME classification table: ", path)
    }
    ptmes <- utils::read.delim(
        path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
        check.names = FALSE
    )
    required <- c("ptme_id", "PID", "ID", "n_cells")
    missing <- setdiff(required, names(ptmes))
    if (length(missing) > 0L) {
        stop("PTME classification table is missing: ", paste(missing, collapse = ", "))
    }
    if (nrow(ptmes) == 0L || anyDuplicated(ptmes$ptme_id) ||
        anyNA(ptmes[, required, drop = FALSE])) {
        stop("PTME classification table must contain complete unique PTME rows")
    }
    ptmes$n_cells <- suppressWarnings(as.numeric(ptmes$n_cells))
    if (any(!is.finite(ptmes$n_cells)) || any(ptmes$n_cells < 1) ||
        any(ptmes$n_cells != floor(ptmes$n_cells))) {
        stop("PTME n_cells must contain positive integers")
    }
    ptmes$n_cells <- as.integer(ptmes$n_cells)
    ptmes
}
```

<!-- source: 04_plot_ptme_isr_panels.R; lines: 58-73 -->

``` r
#' Purpose: Render the PTME structure-size histogram for Supplementary Figure 9C.
#'
#' @param ptmes One-row-per-PTME classification data frame.
#' @return A ggplot object with PTME size measured in cells.
plot_supplementary_figure9c <- function(ptmes) {
    ggplot2::ggplot(ptmes, ggplot2::aes(x = n_cells)) +
        ggplot2::geom_histogram(
            bins = 30L, fill = "#C73E1D", alpha = 0.8, color = "white"
        ) +
        ggplot2::labs(
            x = "PTME size (number of cells)",
            y = "Count"
        ) +
        ggplot2::theme_minimal(base_size = 12)
}
```

<!-- source: 04_plot_ptme_isr_panels.R; lines: 74-94 -->

``` r
#' Purpose: Render all currently supported Figure 4 PTME/ISR plot consumers.
#'
#' @param table_directory Directory written by 01_methods_aligned_orchestration.R.
#' @param output_directory Directory for consumer PDFs.
#' @return Invisibly returns a named list of rendered ggplot objects.
run_figure4_plot_consumers <- function(table_directory, output_directory) {
    require_packages("ggplot2")
    ptmes <- read_ptme_classification(table_directory)
    supplementary_figure9c <- plot_supplementary_figure9c(ptmes)

    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
        filename = file.path(output_directory, "SupplementaryFigure9C_ptme_size_distribution.pdf"),
        plot = supplementary_figure9c, width = 8, height = 6
    )
    invisible(list(supplementary_figure9c = supplementary_figure9c))
}

if (isTRUE(RUN_FIGURE4_PLOT_CONSUMERS)) {
    run_figure4_plot_consumers(TABLE_DIRECTORY, OUTPUT_DIRECTORY)
}
```

## Representative PIR and PSM niches

![Complete published Figure S10 with four PIR and four PSM examples,
each showing IMC, H&E, and a cell-type
map](report_assets/figureS10_publication.png)

**Figure S10. Representative PIR- and PSM-Niches in peritumoral region
of CRLM.** Four PIR examples (left) and four PSM examples (right) pair
IMC overlays, matched H&E, and cell-type maps. PIR overlays show DNA,
EpCAM, CD20, and CD8a; PSM overlays show DNA, EpCAM, CD163, and
Collagen-I. Cell maps identify putative cholangiocytes, CD8-positive T
cells, CD20-positive B cells, CD163-positive macrophages, and
collagen-I-positive stromal cells. Scale bars, 50 µm.

These examples complement Figure 4C/G by showing the morphology and
surrounding immune/stromal composition used to characterize the two
niche states. The [panel manifest](../../docs/figure_table_manifest.tsv)
connects the publication panels to their analysis entry points.

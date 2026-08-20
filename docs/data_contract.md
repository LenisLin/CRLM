# Public data contract

The public data package contains four processed-data payloads with independent
cohort identifiers. The Zenodo package includes a field-level data dictionary
for each object.

## FDZS-1 IMC

`FDZS1_IMC_processed.rds` is a `SingleCellExperiment` with 35 markers x
2,018,260 cells from 35 B/W research-pseudonym patients and 311 ROIs: 102 PT,
92 TC, and 117 IM. `Position` stores cell coordinates in micrometres.

Public analysis fields include `CellID`, ROI `ID`, patient `PID`, `Tissue`,
`MajorType2`, `SubType`, `RFS_event`, `RFS_time`, `Treatment`, `Treatment2`,
and clinicopathological covariates. `RFS_event = 1` denotes recurrence and
`RFS_event = 0` denotes censored/no recurrence. `RFS_time` is measured in
months from completion of resection of the primary colorectal tumor and liver
metastases.

## FDZS-3 scRNA-seq

`FDZS3_scRNA_processed.h5ad` contains 35,900 cells. `X` is 35,900 x 2,000,
and embedded `raw.X` is a 35,900 x 20,237 count matrix used by downstream
analysis. Public samples are `FDZS3-P001` through `FDZS3-P006`. Metadata
includes `patient`, `batch`, `doublet`, `Major_type`, and `Sub_type`.
Technical cell identifiers retain their analytical suffixes.

## FDZS-4 spatial transcriptomics

`FDZS4_ST_processed.h5ad` contains 17,902 Bin100 spots x 56,695 features.
Public samples are `FDZS4-P001` through `FDZS4-P009`. `obs` contains
`sample_id`, `x`, and `y`; `obsm["spatial"]` contains the same native
Stereo-seq coordinates; and `obsm["RCTD_cell_type_abundances"]` contains 30
normalized subtype weights.

## Supplementary tables

`Supplementary_Tables_public.xlsx` contains Supplementary Tables 1-4 in
their manuscript numbering and order. The tables provide patient, ROI,
marker/channel, and reagent records.

## Script interfaces

The IMC entry points consume `FDZS1_IMC_processed.rds`. Figure 7 entry points
consume processed FDZS-3 and FDZS-4 MatrixMarket bundles, registered SpMap
tile annotations, RCTD weights, and offline pathway resources. Each entry
point reads its editable in-script configuration and uses no-argument direct
execution. `config/paths.example.yml` is a path-reference worksheet, and the
required file schemas are specified in `figures/Figure7/README.md`.

## Privacy model

The public objects use cohort-specific pseudonyms and retain technical
cell, ROI, and spot identifiers and spatial coordinates required by the
analytical interfaces. Each payload's public fields are listed above and in
the accompanying Zenodo data dictionary.

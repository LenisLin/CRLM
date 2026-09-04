# Public data contract

The public data package contains four processed-data payloads and twelve
derived analysis artifacts with cohort-specific identifiers. The Zenodo package
includes a field-level data dictionary for every object and table.

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

## Patient-level ISR and cohort tables

`FDZS1_patient_level_ISR.tsv` contains 35 FDZS-1 patient rows with B/W research
pseudonyms, recurrence-free-survival fields, treatment group, saved IMC-derived
ISR components, and available H&E-derived ISR components. The
`figure6_included = 1` identifies the 34 patients with H&E-derived ISR used in
the discovery analysis; `wsi_id` gives the corresponding H&E slide identifier.

`FDZS2_patient_level_ISR_RFS.tsv` contains 95 independent-test patient rows with
`FDZS2-P*` patient IDs and `FDZS2-WSI*` slide IDs. It provides the patient-level
fields consumed by the Figure 6 workflow: ISR, ISR group, RFS time/event,
treatment, TBS, CRLM number and size, Fong score, KRAS mutation code, age, and
gender code. `SpMap_patient_cohort_assignment.tsv` lists the 34 discovery and 95
independent-test analysis assignments.

`ISR_definition_and_cutoff.tsv` records the operational IMC- and H&E-derived ISR
formulas and the discovery-cohort median H&E-derived ISR cutoff applied to both
cohorts.

## ROI-level cell composition

`FDZS1_ROI_cell_composition.tsv.gz` contains 311 rows, one per FDZS-1
patient-ROI-tissue combination. Columns provide the B/W patient and ROI IDs,
tissue region, total cell count, and the fraction of every public `SubType`.

## SpMap artifacts

`SpMap_reference_tile_labels.tsv.gz` contains the canonical 61,628 four-class
reference tiles. `SpMap_tile_partitions.tsv.gz` adds the development or common
internal-holdout assignment and the five-fold validation membership for each
development tile.

`SpMap_model_performance_5seeds.tsv` contains class-specific and overall
performance values for seeds 24, 101, 202, 303, and 404.
`SpMap_confusion_matrix.tsv` contains the seed-24 pooled out-of-fold four-class
confusion counts and true-class fractions. The OOF archive contains the 25
fold-specific validation-prediction tables. The final-model archive contains
exactly five `primary_1280` fold checkpoints, one per fold, with portable model
manifests. The model uses paired 768- and 512-dimensional CONCH
representations (`768 + 512 = 1280` input) in a `1280 -> 450 -> 4` MLP with
ReLU, dropout 0.25, no BatchNorm, and validation-accuracy checkpoint selection.
The CONCH archive contains the canonical feature tar shards and their member
inventories.

## Script interfaces

The IMC entry points consume `FDZS1_IMC_processed.rds`. Figure 1B consumes
Supplementary Tables 1 and 5 and writes a treatment-stratified RFS PDF and
summary table. Figure 6E-F consumes those tables plus
`FDZS1_patient_level_ISR.tsv` and writes the cohort-specific continuous-ISR Cox
table and two PDFs. Figure 7 entry points consume processed FDZS-3 and FDZS-4
MatrixMarket bundles, registered SpMap tile annotations, RCTD weights, and
offline pathway resources. Each entry point reads its editable in-script
configuration and uses no-argument direct execution. The required file schemas
are specified in `figures/Figure7/README.md`. Figure 4, Figure 5, and Figure 6
provide dedicated entry points for generating public derived artifacts from
their documented source tables and objects.

## Privacy model

The public objects and derived tables use cohort-specific pseudonyms and retain technical
cell, ROI, and spot identifiers and spatial coordinates required by the
analytical interfaces. Each payload's public fields are listed above and in
the accompanying Zenodo data dictionary.

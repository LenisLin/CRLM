# Figure 4

This directory provides `01_methods_aligned_orchestration.R` for **Figure 4.
Associations between cholangiocyte-associated peritumor niches and treatment-
stratified RFS in CRLM.** The script generates cellular-neighborhood, PTME,
and ISR tables from an IMC SCE. `02_export_roi_cell_composition.R` exports the
public ROI-level cell-composition matrix from the processed FDZS-1 object.

## Run

```text
Rscript figures/Figure4/01_methods_aligned_orchestration.R
Rscript figures/Figure4/02_export_roi_cell_composition.R
```

Before running from the repository root, edit `SCE_INPUT`,
`OUTPUT_DIRECTORY`, and `SCRIPT_DIRECTORY` near the top of the script. Invoke
the entry point directly with the no-argument command above. For the public
ROI export, edit `EXPORT_CONFIG` and set `RUN_ROI_EXPORT` to `TRUE`.

The script requires `SingleCellExperiment`, `SummarizedExperiment`, `FNN`,
`igraph`, and `geometry`. It uses the shared interfaces in
`../../src/imc/coordinates.R`, `../../src/imc/cellular_neighborhoods.R`, and
`../../src/imc/ptme.R`.

## Input schema

The SCE `colData()` must contain complete `CellID`, `PID`, `ID`, `Tissue`,
`SubType`, `Position`, `MajorType2`, and `Treatment` fields. `CellID` must be
unique. `Treatment` is the explicit clinical field, uses `Chemo` or `Combo`,
and is constant within each patient. Execution requires the listed packages,
the three shared IMC interface files, complete required fields, and PT cells.

## Fixed analysis parameters

The analysis is restricted to PT cells. Neighborhood composition uses the
unweighted subtype composition of the 20 nearest cells within each patient and
ROI. `MajorType2 == "EC"` identifies the cholangiocyte mask. PTMEs are derived
from the induced cholangiocyte subgraph of an all-PT-cell, within-ROI Delaunay
graph; components containing at least 5 cholangiocytes receive a 22-micrometre
expanded convex hull.

PTMEs are classified from the ratio
`(SC_COLLAGEN + Macro_CD163) / (B + CD8T)`: PIR is below 1, PSM is above 1,
and a ratio equal to 1 is unclassified. ISR is calculated as
`(N_PIR + 1) / (N_PSM + 1)` at patient-ROI and patient levels, then clipped to
`[0.1, 10]`. All eligible PT ROIs and patients are included. Units with no PIR
or PSM structure have zero counts and ISR 1. ISR counts use PIR and PSM
structures.

Analysis units are PT cells for neighborhood composition, the cholangiocyte
mask, and PTME membership; PTMEs for niche classification; and PT
patient-ROIs or patients for ISR summaries. Each derived table carries the
explicit patient treatment value.

## Outputs

| Output | Analysis unit and content |
|---|---|
| `figure4_methods_parameters.tsv` | Fixed parameter record. |
| `figure4_pt_cn_composition.tsv` | One row per PT cell with 20-neighbor subtype proportions. |
| `figure4_pt_cholangiocyte_mask.tsv` | One row per PT cell with the `MajorType2 == "EC"` cholangiocyte mask. |
| `figure4_ptme_membership.tsv` | PT-cell membership in expanded-hull PTMEs. |
| `figure4_ptme_classification.tsv` | One row per PTME with subtype counts, fractions, ratio, niche class, and treatment. |
| `figure4_roi_isr.tsv` | One row per eligible PT patient-ROI with PIR/PSM counts, raw and clipped ISR, and treatment. |
| `figure4_patient_isr.tsv` | One row per eligible PT patient with PIR/PSM counts, raw and clipped ISR, and treatment. |
| `FDZS1_ROI_cell_composition.tsv.gz` | One row per patient-ROI-tissue group with total cells and public SubType fractions. |

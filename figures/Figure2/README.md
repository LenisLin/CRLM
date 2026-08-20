# Figure 2: Cellular composition of CRLM is associated with recurrence

This directory contains the interactive R entry points for Figure 2B, Figure
2D, and Supplementary Figures 3-5. The `imc` package family and observed
versions are listed in [`../../envs/package_versions.yml`](../../envs/package_versions.yml).

## Shared input contract

Each entry point reads the directly editable `sce_path` and `output_dir`
variables near the top of the script. Set those variables, then invoke the
script from the repository root without arguments:

```text
Rscript figures/Figure2/01_plot_cell_subpopulation_heatmap.R
Rscript figures/Figure2/02_plot_tissue_composition.R
Rscript figures/Figure2/03_analyze_treatment_stratified_abundance.R
```

The analyses use cells from the `IM`, `PT`, and `TC` tissue compartments.
`colData()` provides the listed cell and patient fields, and `metadata(sce)`
provides the named colour vectors. Each command creates its selected output
directory.

## Ordered entry points

### 1. Figure 2B cell-subpopulation heatmap

```text
Rscript figures/Figure2/01_plot_cell_subpopulation_heatmap.R
```

**Required SCE contents.** `colData()` fields `Tissue`, `SubType`, and
`MajorType2`; `metadata(sce)$color_vectors$MajorType2`; and the annotation
markers `CD45`, `CD20`, `EpCAM`, `CD11c`, `CD14`, `CD16`, `HLADR`, `CD68`,
`CD11b`, `CD163`, `CD169`, `CD57`, `CollagenI`, `AlphaSMA`, `FAP`,
`Vimentin`, `CD3`, `CD4`, `CD8a`, and `FoxP3`.

**Outputs.** `figure2B_cell_subpopulation_marker_heatmap.pdf` and
`figure2B_cell_subpopulation_summary.tsv`.

**Analysis contract.** The script calculates mean marker expression for each
`SubType`, standardizes each marker across subtypes, and reports subtype cell
counts and fractions using all retained cells.

### 2. Supplementary Figure 3 tissue composition

```text
Rscript figures/Figure2/02_plot_tissue_composition.R
```

**Required SCE contents.** `colData()` fields `PID`, `ID`, `Tissue`,
`SubType`, `RFS_event`, `Gender`, `KRAS_mutation`, `CRC_site`,
`Differential_grade`, `T_stage`, `Lymph_positive`, `RFS_time`, `Age`, `TBS`,
`CRLM_number`, `CRLM_size`, `CEA`, and `CA199`; and
`metadata(sce)$color_vectors$SubType`.

**Outputs.** `figureS3A_tissue_composition_pie_charts.pdf`,
`figureS3B_roi_subtype_fractions.pdf`,
`figureS3C_patient_clinical_composition_heatmap.pdf`,
`figureS3A_pooled_tissue_subtype_fractions.tsv`,
`figureS3B_roi_subtype_fractions.tsv`, and
`figureS3C_patient_subtype_fractions.tsv`.

**Analysis contract.** S3A pools cells within tissue, S3B calculates subtype
fractions within ROIs, and S3C pools cells within patients and displays the
listed clinical annotations.

### 3. Figure 2D and Supplementary Figures 4-5 treatment-stratified abundance

```text
Rscript figures/Figure2/03_analyze_treatment_stratified_abundance.R
```

**Required SCE contents.** `colData()` fields `PID`, `ID`, `Tissue`,
`SubType`, `RFS_time`, `RFS_event`, `Treatment`, `Gender`, `KRAS_mutation`,
`Differential_grade`, `Lymph_positive`, `fong_score`, `TBS`, `CRLM_number`,
`CEA`, and `CA199`. Treatment strata are `Chemo` and `Combo`; the output
tables also include `Overall`.

**Outputs.** `patient_tissue_subtype_fractions.tsv`,
`treatment_stratified_spearman_results.tsv`,
`treatment_stratified_wilcoxon_results.tsv`,
`treatment_stratified_univariable_cox_results.tsv`,
`figure2D_IM_recurrence_wilcoxon.tsv`, and
`figure2D_IM_recurrence_boxplots.pdf`. Supplementary Figure 4 outputs use
the `figureS4[A-C]_<analysis>_IM_<treatment>.pdf` pattern, and Supplementary
Figure 5 outputs use the `figureS5[A-C]_<analysis>_TC_<treatment>.pdf`
pattern for `Chemo` and `Combo`.

**Analysis contract.** The script calculates subtype fractions within each
ROI and averages them with equal ROI weight within each patient, tissue, and
subtype. It reports Spearman associations for `fong_score`, `TBS`,
`CRLM_number`, `CEA`, and `CA199`; two-sided Wilcoxon rank-sum comparisons
for `RFS_event`, `Treatment`, `Gender`, `KRAS_mutation`,
`Differential_grade`, and `Lymph_positive`; and univariable Cox models using
`RFS_time` and `RFS_event`. Spearman and Wilcoxon tables apply Benjamini-
Hochberg adjustment within their tissue-treatment analysis strata; Cox
predictors are scaled by one sample SD and report likelihood-ratio-test P
values with Benjamini-Hochberg adjustment. `RFS_event = 1` denotes
recurrence.

Figure 2D uses IM patient-level mean ROI fractions for `CD4T`, `CD8T`,
`Macro_CD169`, and `Macro_CD163` in the `Chemo` and `Combo` strata. Each
treatment-subtype comparison uses a two-sided Wilcoxon rank-sum test with at
least six patients and at least three patients in each `RFS_event` group.

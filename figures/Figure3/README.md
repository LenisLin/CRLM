# Figure 3

This directory provides the analysis scripts for **Figure 3. Cellular
neighborhood patterns in tumor core and invasive margin.** It also contains the
Supplementary Figure 6 workflow for pairwise cell-cell interaction frequency
matrices across CRLM regions.

## Figure 3B-F cellular neighborhoods

Edit `figure3_config$sce_path`, `output_dir`, and `script_dir`, set
`figure3_config$run` to `TRUE`, and run from the repository root:

```text
Rscript figures/Figure3/01_cellular_neighborhoods.R
```

The script requires `SingleCellExperiment`, `dplyr`, `tidyr`, `FNN`,
`pheatmap`, `ggplot2`, and `survival`, and uses
`../../src/imc/coordinates.R` for `position_to_xy()`.

The SCE `colData()` must contain `CellID`, `ID`, `PID`, `Tissue`, `Position`,
`SubType`, `Treatment`, `RFS_time`, and `RFS_event`. `CellID` must be unique,
and the clinical fields must define one `Treatment`, `RFS_time`, and
`RFS_event` value per patient. `Position` is parsed as micrometre coordinates.
The input contains complete listed fields, unique `CellID` values, one
treatment/RFS record per patient, and at least 10 cells with valid within-ROI
neighborhoods.

Neighborhoods are calculated separately within ROI using inverse-distance-
weighted fractions of the 20 nearest cells. The fixed settings are a top-level
seed of 619, k-means seed of 42, 10 clusters, `nstart = 20`, `iter.max = 100`,
12 workers, and reassignment of clusters with fewer than 5,000 cells to
`CN_Other`. Analysis units are cells for neighborhood assignment, patients for
CN fractions and clinical summaries, and IM and TC tissues for the
treatment-stratified summaries. The nominal Wilcoxon tests compare patient CN
fractions by `RFS_event` within treatment and CN; the Cox models are
`Surv(RFS_time, RFS_event) ~ scale(CN fraction) * Treatment`.

The output directory contains:

| Output | Content and analysis unit |
|---|---|
| `cn_results/peritumor_cn_results.csv` | One cellular-neighborhood assignment per cell. |
| `cn_results/Spatial neighbors of All with k=20 neighbors.csv` | One within-ROI 20-neighbor subtype-fraction profile per cell. |
| `Cellular Neighbors celltype fraction matrix.tsv` | CN-by-subtype composition matrix. |
| `Cellular Neighbors celltype fraction heatmap.pdf` | Heatmap of the CN-by-subtype composition matrix. |
| `cn_results/patient_level_cn_fractions.csv` | Patient-by-tissue-by-CN cell counts and fractions. |
| `cn_results/IM_CN_nominal_wilcoxon.csv`, `cn_results/TC_CN_nominal_wilcoxon.csv` | Treatment- and CN-specific patient-level Wilcoxon results. |
| `IM_CN_Frequency_by_treatment_group.pdf`, `TC_CN_Frequency_by_treatment_group.pdf` | Patient CN fractions by treatment and RFS event. |
| `cn_results/IM_CN_cox_results.csv`, `cn_results/TC_CN_cox_results.csv` | CN fraction and CN-by-treatment Cox-model estimates. |
| `ForestPlot_Univariable_CN_Fraction_IM.pdf`, `ForestPlot_Univariable_CN_Fraction_TC.pdf` | Forest plots of the Cox-model estimates. |

## Supplementary Figure 6A-C pairwise interactions

Supplementary Figure 6 uses an ordered producer-consumer workflow. First edit
`figure_s6_precursor_config`, set `RUN` to `TRUE`, and generate ROI interaction
P-value matrices:

```text
Rscript figures/Figure3/supplementary/FigureS6/00_generate_interaction_pvalues.R
```

The precursor consumes the FDZS-1 SCE fields `CellID`, `PID`, `ID`, `Tissue`,
`SubType`, and `Position`, sources `src/imc/coordinates.R` and
`src/imc/spatial_interactions.R`, and applies a 22-micrometre distance,
1,000 permutations, tissues IM/PT/TC, and the fixed 21-subtype order. Its
analysis unit is an ROI-specific cell-pair permutation, and it writes one
`ClosePvalue.csv` matrix under
`permutation_<Tissue>/<ROI>/ClosePvalue/`.

Then edit `figure_s6_config$permutation_root` and `output_dir`, set
`figure_s6_config$run` to `TRUE`, and run the consumer:

```text
Rscript figures/Figure3/supplementary/FigureS6/01_pairwise_interactions.R
```

The script requires `ComplexHeatmap` and `circlize`. For every ROI in the fixed
roster, it reads `ClosePvalue/ClosePvalue.csv` from the corresponding
`permutation_<tissue>/<ROI>` directory. The fixed panel mapping is A = IM,
B = PT, and C = TC. The roster contains 117 IM ROIs, 102 PT ROIs, and 69 TC
ROIs; these counts are checked before execution.

Every saved table must have the following row and column order:

```text
Macro_CD169, Macro_HLADR, Mono_CD11c, SC_COLLAGEN, CD8T, CD4T, B,
SC_FAP, Macro_CD163, UNKNOWN, SC_Vimentin, Mono_Intermediate, Mono_Classic,
TC_EpCAM, SC_aSMA, TC_Ki67, NK, Treg, Macro_CD11b, TC_CAIX, TC_VEGF
```

The consumer requires one schema-matched table for every rostered ROI. For each
cell-type pair and ROI, p-values at or below 0.05 are scored as 1, p-values
above 0.95 as -1, and all other p-values as 0. Tissue-specific matrices are
the ROI-average signed scores, using the fixed tissue roster as denominator.
The analysis unit is an ROI-level saved p-value matrix.

The output directory contains `FigureS6_saved_pvalue_roster.tsv`, three signed
frequency matrices named `Cell-Cell Interaction Frequency in <tissue>.tsv`,
and the matching PDF heatmaps `Cell-Cell Interaction Frequency in <tissue>.pdf`
for IM, PT, and TC.

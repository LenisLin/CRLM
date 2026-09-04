# Figure 1

This directory provides the patient-level treatment/RFS entry point for Figure
1B and documents the assembled publication panels for Figure 1 and
Supplementary Figures 1-2.

## Figure 1B entry point

Edit `FIGURE1B_CONFIG$supplementary_workbook`, set `RUN_FIGURE1B <- TRUE`, and
run from the repository root:

```text
Rscript figures/Figure1/01_plot_treatment_rfs.R
```

The script reads `Supplementary Table 1` and `Supplementary Table 5` directly
from the Supplementary Tables workbook. The first table supplies the
35-patient FDZS-1 discovery cohort and the second supplies the 95-patient
FDZS-2 test cohort. Displayed recurrence and treatment fields are normalized
to `rfs_event` (`1=recurrence`, `0=censored`) and `treatment` (`Chemo` or
`Combo`) before analysis. The script requires one complete row per unique
patient and both treatment groups in each cohort.

The entry point writes the two-panel vector PDF
`Figure1B_treatment_RFS_revised.pdf` and
`Figure1B_treatment_RFS_summary.tsv`, which records cohort size, group counts,
event count, and the nominal two-sided log-rank P value. The PDF uses the
publication palette (`#0073C2FF` and `#EFC000FF`).

## Assembled publication panels

- **Figure 1A. Study workflow.** Study-design schematic.
- **Figure 1C. Custom 35-plex IMC antibody-panel overview.** Categorical
  overview of the custom 35-plex IMC panel, with supporting antibody
  information in Supplementary Table 3.
- **Figure 1D. Representative IMC region with cell-type map and marker
  overlays.** Representative IMC cell-type map with marker overlays.
- **Supplementary Figure 1. Spatial characterization of CRLM across tissue
  regions.** WSI, IMC, and H&E composite.
- **Supplementary Figure 2. Representative IMC marker images for cell type
  annotation in CRLM.** H&E and IMC marker-overlay composite.

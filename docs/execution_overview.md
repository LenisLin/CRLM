# Execution overview

This repository uses 44 figure-owned and shared R/Python scripts. Each Figure
README specifies script order, editable configuration, required inputs,
outputs, analysis units, fixed parameters, and package-version family.

All entry points are invoked from the repository root without command-line
flags or positional arguments. Edit the top-level path configuration in each
script and set its run switch where present. `config/paths.example.yml` provides
shared path references.

## IMC branch

1. Provide `FDZS1_IMC_processed.rds` and use the `imc` package family.
2. Run the three Figure 2 scripts to generate cell-subpopulation, tissue-composition,
   and treatment-stratified abundance outputs.
3. Run Figure 3 to construct inverse-distance-weighted 20-neighbor cellular
   neighborhoods. For Supplementary Figure 6, run
   `00_generate_interaction_pvalues.R` to produce ROI `ClosePvalue.csv`
   matrices before `01_pairwise_interactions.R` aggregates the fixed rosters.
4. Run Figure 4 to generate CN, PTME, ISR, epithelial-identity, and supported
   Supplementary Figure 9C outputs from the IMC object and emitted tables.
5. Run `figures/Figure4/02_export_roi_cell_composition.R` to generate the public
   311-ROI cell-composition matrix from the processed FDZS-1 object.

## SpMap and clinical branch

1. `src/spmap` provides reference-tile construction, preprocessing, paired
   CONCH feature extraction, grouped splits, CUDA MLP training, WSI inference,
   and H&E-derived ISR aggregation.
2. Configure `figures/Figure5/02_run_spmap_workflow.py` for the reference-
   manifest, feature, training/evaluation, and WSI-inference stages.
3. Run the Figure 5 count producer and model-development plot consumer,
   then the Supplementary Figure 11 consumers with their configured
   panel-specific input tables.
4. Run Figure 1B directly from Supplementary Tables 1 and 5. Run Figure 6B-D
   with its prepared 34-patient discovery and 95-patient test tables. Run the
   Figure 6E-F continuous-ISR Cox entry point with Supplementary Tables 1 and
   5 plus the FDZS-1 patient-level ISR table.
5. Run `figures/Figure5/03_prepare_public_spmap_artifacts.py` to prepare the
   public SpMap labels, partitions, metrics, OOF predictions, feature shards,
   and five `primary_1280` fold checkpoints.
6. Run `figures/Figure6/03_prepare_public_isr_tables.py` to generate the public
   patient-level ISR/RFS, cohort-assignment, and cutoff tables.

## scRNA-seq and spatial-transcriptomics branch

1. Run `figures/Figure7/01_integrate_annotate_scrna.py` with the processed
   FDZS-3 MatrixMarket bundle to generate annotated scRNA objects and an R
   interoperability bundle.
   Use `figures/Figure7/05_plot_public_scrna_annotations.py` to render Figure
   7B and Supplementary Figure 15 annotation panels directly from the released
   H5AD.
2. Run the Supplementary Figure 16 CellChat script with the R bundle.
3. Run `02_run_rctd.R` with the processed FDZS-4 Bin100 bundle and scRNA
   reference to generate normalized RCTD weights.
4. Run `03_map_spmap_plot_composition.R` with registered SpMap tiles and RCTD
   weights to generate spatial-composition outputs.
5. Run `04_run_deg_gsea.R` with RCTD-intersection spot labels and the offline
   gene-set resource to generate DEG and GSEA outputs.

`config/analysis_parameters.yml` and the Figure READMEs define fixed parameters
and analysis units. The Figure READMEs provide the exact no-argument commands.

## Manifest statuses

The figure/table manifest uses the following factual statuses:

- `executable`: a listed entry point produces the analytical item from its
  documented inputs;
- `partial`: a listed entry point produces a documented analytical component
  using additional processed inputs;
- `assembled`: an assembled publication panel or descriptive item;
- `documented`: a manuscript item described by its Figure README; and
- `data_only`: a public data or table artifact.

The provenance manifest contains all 44 R/Python scripts and the native and
derived public data objects. It uses `current` for analysis code,
`implementation` for figure-specific table-generation code, `interface` for
shared computational modules, `processed_data` for public native objects, and
`derived_data` for public tables, predictions, features, and model artifacts.
The `decision` field records `release` for all public entries.

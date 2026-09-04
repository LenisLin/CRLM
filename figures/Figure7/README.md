# Figure 7

This directory provides the Figure 7 workflow for processed FDZS-3 scRNA-seq
and FDZS-4 Stereo-seq data. The workflow constructs the annotated scRNA
reference, estimates Bin100 spot composition with RCTD, maps registered SpMap
classes, compares PIR-Niche and PSM-Niche spots, and generates the Figure 7
and Supplementary Figures 15-18 analytical outputs.

Key package versions are listed in `envs/package_versions.yml`.

## Workflow order

Edit each script's top-level `CONFIG`, then run the following no-argument
commands from the repository root in order:

1. `python figures/Figure7/01_integrate_annotate_scrna.py` creates the annotated FDZS-3 reference,
   Figure 7B and Supplementary Figure 15 PDFs, and the R MatrixMarket bundle.
2. `Rscript figures/Figure7/02_run_rctd.R` uses the scRNA R bundle and the FDZS-4 Bin100 bundle to
   generate normalized subtype RCTD weights.
3. `Rscript figures/Figure7/03_map_spmap_plot_composition.R` combines the Bin100 bundle, registered
   SpMap tile annotations, and RCTD weights to generate Figure 7C-E and
   Supplementary Figures 17-18 PDFs.
4. `Rscript figures/Figure7/04_run_deg_gsea.R` uses the RCTD-intersection spot table with the Bin100
   bundle and versioned pathway resources to generate Figure 7F-H PDFs.

For direct plotting from the released processed object, run
`python figures/Figure7/05_plot_public_scrna_annotations.py` after editing its
`CONFIG`. It reads `FDZS3_scRNA_processed.h5ad` directly and renders Figure 7B
and Supplementary Figure 15.

`Rscript figures/Figure7/supplementary/FigureS16/01_run_plot_cellchat.R`
uses the annotated scRNA R bundle after step 1 to generate the five
Supplementary Figure 16 PDFs.

## Required inputs

For `01_integrate_annotate_scrna.py`, set `CONFIG["input_dir"]`, `gtf`,
`output_dir`, `prefix`, and `excluded_epithelial_samples`. The fixed prefix is
`processed`, and the documented epithelial selection is
`excluded_epithelial_samples=("FDZS3-P006",)`. The processed FDZS-3 MatrixMarket bundle
contains `processed_matrix.mtx[.gz]`, `processed_features.tsv[.gz]`,
`processed_barcodes.tsv[.gz]`, and `processed_metadata.csv`. The on-disk
matrix has genes as rows and cells as columns. Metadata contains `patient`,
`doublet`, and `sample_id`.

For `05_plot_public_scrna_annotations.py`, set `CONFIG["input_h5ad"]` and
`CONFIG["output_dir"]`. The H5AD must provide `Major_type`, `Sub_type`, the
two-dimensional `X_umap` embedding, and `raw.X` with the configured marker
genes. The public object contains 35,900 cells. This consumer uses stored
annotations and coordinates; lineage-specific views subset the stored global
UMAP.

For `02_run_rctd.R`, set `CONFIG$sc_bundle`, `st_bundle`, `output_dir`,
`st_functions`, and `scrna_functions`. For repository-root execution, set the
helper paths to `src/st/st_functions.R` and
`src/scrna/cellchat.R`. The FDZS-4 Bin100 bundle contains
`expression_profile.mtx`, `metadata.csv`, `row_names.csv`, and
`column_names.csv`. The on-disk matrix has spots as rows and genes as columns;
`src/st/st_functions.R` transposes it for the spatial object. Metadata aligns
with spot identifiers and provides `core_name`, `x`, and `y`.

For `03_map_spmap_plot_composition.R`, set `CONFIG$st_bundle`,
`tile_annotations`, `rctd_weights`, `output_dir`, and `st_functions`.
Set `st_functions` to `src/st/st_functions.R` for repository-root
execution. Registered SpMap tile annotations provide `sample_id`,
`tile_left`, `tile_right`, `tile_top`, `tile_bottom`, and `pred_class` in the
Bin100 coordinate system. The RCTD input is
`rctd_normalized_weights.tsv`, with a string `spot_id` column followed by
subtype proportions.

For `04_run_deg_gsea.R`, set `CONFIG$st_bundle`, `spot_labels`, `gene_sets`,
`resource_manifest`, `output_dir`, and `st_functions`. Set `st_functions` to
`src/st/st_functions.R` for repository-root execution. `spot_labels` is
`rctd_intersection_spots.tsv`. `gene_sets` is an RDS list named `Hallmark`,
`KEGG`, and `Reactome`; `resource_manifest` provides
one row per collection with `collection` and nonempty `resource_version`.

For `supplementary/FigureS16/01_run_plot_cellchat.R`, set
`CONFIG$sc_bundle`, `output_dir`, and `scrna_functions`. Set
`scrna_functions` to `src/scrna/cellchat.R` for repository-root
execution.

## Fixed settings and analysis units

- The scRNA reference contains 35,900 singlet cells. It uses normalization to
  10,000, log1p transformation, 2,000 patient-batched highly variable genes,
  scale maximum 10, 50 PCs, Harmony integration, and major clustering with 10
  neighbors, 40 Harmony PCs, and Leiden resolution 0.8. B/plasma
  subclustering uses resolution 0.8; T/NK and myeloid subclustering use 1.0.
  The reference uses the inherited Harmony representation, major and subtype
  annotation maps, inferCNV window 250, and specified cholangiocyte CNV clusters.
- RCTD balances 100 cells per eligible subtype, uses a minimum of 10 cells,
  seed 619, eight cores, `full` mode, and normalized weights.
- SpMap classes are `0=TC`, `1=PIR-Niche`, `2=PSM-Niche`, and `3=Other PT`.
  Registered-tile assignment uses inclusive containment, the first source row
  for multiple containing tiles, and nearest tile top-left selection within
  each sample.
- The mapped-spots population contains 18,186 Bin100 spots: 463 TC, 382 PIR,
  1,526 PSM, and 15,815 Other PT. The RCTD-intersection population contains
  17,902 Bin100 spots: 449 TC, 377 PIR, 1,522 PSM, and 15,554 Other PT.
  Figure 7D-E, Supplementary Figure 18, DEG, and GSEA use the
  RCTD-intersection population. Figure 7E and Supplementary Figure 18 use
  Bin100 spots as observations, two-sided Wilcoxon rank-sum tests, and nominal
  P-value labels. Figure 7E displays `B_IGHM`, `CD8T_GZMK`, `Fibroblast`,
  `Macro_CD163`, and `Macro_SPP1`.
- DEG compares 377 PIR-Niche and 1,522 PSM-Niche Bin100 spots with the Seurat
  Wilcoxon test, `min.pct = 0.01`, `logfc.threshold = 0`, and both expression
  directions. The GSEA rank table contains 2,559 genes after the
  detection-fraction filter. GSEA uses Hallmark, KEGG, and Reactome with
  `minSize = 3`, `maxSize = 500`, `eps = 0`, and nominal pathway P values.
  Figure 7H displays `HALLMARK_ANGIOGENESIS` and the exact Reactome pathway
  `REACTOME_MTORC1_MEDIATED_SIGNALLING`.
- CellChat uses the 12-subtype input universe, the cholangiocyte score genes
  `CA9`, `SLC2A1`, `LDHA`, `PGK1`, `ENO1`, and `ALDOA`, the 0.30 and 0.70
  quantiles, `CellChatDB.human`, `raw.use = TRUE`, a minimum of 10 cells,
  pathway probability, and network aggregation.

## Outputs

`01_integrate_annotate_scrna.py` writes `integrated_scRNA_stereoseq.h5ad`,
`major_anno_all.h5ad`, `final_annotated_scRNA.h5ad`,
`scrna_run_summary.tsv`, Figure 7B and Supplementary Figure 15 PDFs, and the
`r_bundle` MatrixMarket directory.

`05_plot_public_scrna_annotations.py` writes the Figure 7B annotation UMAP and
eight Supplementary Figure 15 major-type, myeloid, T/NK, and B/plasma UMAP or
marker-dot-plot PDFs from the released H5AD.

`02_run_rctd.R` writes `rctd_normalized_weights.tsv`,
`rctd_spot_audit.tsv`, `rctd_reference_subtype_counts.tsv`, and
`rctd_run_summary.tsv`.

`03_map_spmap_plot_composition.R` writes `spmap_spot_assignments.tsv`,
`spmap_mapping_summary.tsv`, `rctd_intersection_spots.tsv`,
`rctd_intersection_summary.tsv`, composition tables,
`celltype_abundance_spot_wilcoxon.tsv`, and Figure 7C-E and Supplementary
Figure 17-18 PDFs.

`04_run_deg_gsea.R` writes complete DEG and GSEA TSVs,
`gsea_rank_table_2559_genes.tsv`, resource and run summaries,
`Figure7G_selected_pathways.tsv`, and Figure 7F-H PDFs.

`supplementary/FigureS16/01_run_plot_cellchat.R` writes `cellchat_final.rds`,
state, parameter, interaction, and activity tables, and the five
Supplementary Figure 16 PDFs.

# Figure 5

This directory provides the SpMap workflow, Figure 5B count
preparation, and Supplementary Figure 11 evaluation entry points.

## Entry points

| Panel | Entry point | Required inputs | Outputs | Analysis unit |
|---|---|---|---|---|
| Figure 5A/C workflow | `02_run_spmap_workflow.py` | Stage-specific paths and adapter callables in `CONFIG` | Selected reference manifests, paired features, fold checkpoints/metrics, WSI predictions, and WSI ISR tables | structure/ROI tile, parent group, model evaluation, WSI tile, and WSI |
| Figure 5B count table | `01_prepare_panel_b_counts.py` | Canonical tile-label manifest; development and internal-holdout split manifests | `Figure5B_training_count_provenance.tsv` | 256 x 256 tile nested in `parent_id` |
| Figure 5B count plot | `04_plot_panel_b_counts.py` | `Figure5B_training_count_provenance.tsv` | Nested-donut PDF and plotted-data TSV | Tile and parent-group split |
| Public SpMap artifacts | `03_prepare_public_spmap_artifacts.py` | Canonical tile/split manifests; selected evaluation runs; final `primary_1280` run; canonical CONCH feature shards | Public label, partition, metric, OOF-prediction, feature, and five-fold final-model artifacts | Tile, parent group, seed/fold model evaluation, and final fold model |
| Figure S11A | `supplementary/FigureS11/01_plot_fold_holdout_metrics.R` | Five-row fold-validation and five-row common-internal-holdout metric CSVs | `FigureS11A_fold_holdout_metrics.pdf`, plotted-data TSV, summary TSV | Model evaluation |
| Figure S11B-C | `supplementary/FigureS11/02_plot_oof_class_metrics.py` | Four-class pooled out-of-fold confusion count and row-proportion tables; five-seed class and macro-F1 metric tables | Two panel PDFs with plotted-data and statistics TSVs | Tile for S11B; seed-by-class metric for S11C |
| Figure S11D | `supplementary/FigureS11/03_plot_isr_correlation.py` | One-row-per-patient IMC ISR table; discovery H&E ISR table | Matched-patient TSV, correlation-statistics TSV, `FigureS11D_isr_correlation.pdf` | Matched patient |

Edit each script's top-level `CONFIG` and run it from the repository root
without arguments. For a complete upstream SpMap run, first configure
`02_run_spmap_workflow.py`, supply the required adapter callables, and enable
the selected stages. Its fixed stage order is reference-manifest construction,
paired CONCH extraction, grouped training/evaluation, then WSI inference and
ISR aggregation. Run the count and Supplementary Figure 11 consumers after
their input tables are available. Python tabular preparation and plotting use
`Spatial_py`; SpMap feature/training/inference stages use the CUDA `OSMOSIS_V2`
environment; S11A uses `Spatial` with R 4.2.2. Package versions are recorded
in `envs/package_versions.yml`.

```text
python figures/Figure5/02_run_spmap_workflow.py
python figures/Figure5/01_prepare_panel_b_counts.py
python figures/Figure5/04_plot_panel_b_counts.py
python figures/Figure5/03_prepare_public_spmap_artifacts.py
Rscript figures/Figure5/supplementary/FigureS11/01_plot_fold_holdout_metrics.R
python figures/Figure5/supplementary/FigureS11/02_plot_oof_class_metrics.py
python figures/Figure5/supplementary/FigureS11/03_plot_isr_correlation.py
```

## SpMap workflow configuration

`02_run_spmap_workflow.py` exposes `CONFIG["enabled"]` switches for
`build_reference_manifest`, `extract_reference_conch`, `train_evaluate`, and
`infer_wsi`. `CONFIG["reference"]` holds registered structure/ROI inputs,
image sizes, image loading and CONCH adapters, and manifest/feature outputs.
`CONFIG["training"]` holds the canonical tile manifest, canonical CONCH tar
root, and training output directory. `CONFIG["inference"]` holds WSI paths,
checkpoint and `model_id`, WSI level, normalization/CONCH adapters, and output
directory. CONCH extraction, training, and inference use the configured CUDA
GPU. Training loads the canonical tar-shard features once into RAM before the
fold and epoch loops.

## Figure 5B training-count contract

`01_prepare_panel_b_counts.py` requires three CSV files:

- Canonical tile manifest: `tile_id`, `class_id`, `class_name`, `parent_id`.
- Development manifest: `tile_id`, `class_id`, `parent_id`.
- Internal-holdout manifest: `tile_id`, `class_id`, `parent_id`.

The four class IDs are `0=TC`, `1=PIR`, `2=PSM`, and `3=OTHER_PT`. The
canonical manifest contains 61,628 tiles: TC 23,894; PIR 9,690; PSM 9,060;
and OTHER_PT 18,984. The development and internal-holdout manifests contain
49,001 and 12,627 tiles, respectively, with disjoint `parent_id` groups.
The output columns are `count_group`, `category`, `count`, `analysis_unit`,
`source_path`, and `selection_rule`. Set `CONFIG["tile_manifest"]`,
`CONFIG["development"]`, `CONFIG["holdout"]`, and `CONFIG["output"]` before
running the no-argument command above.

Set `CONFIG["counts"]` and `CONFIG["output_dir"]` in
`04_plot_panel_b_counts.py` after generating the count table. This consumer
renders the four-class tile counts and grouped development/internal-holdout
split.

For S11A, set `CONFIG$validation`, `holdout`, and `output_dir`. For S11B-C,
set `CONFIG["confusion_counts"]`, `confusion_row_proportions`,
`seed_class_metrics`, `seed_overall_metrics`, and `output_dir`. For S11D, set
`CONFIG["imc_scores"]`, `CONFIG["he_scores"]`, and `CONFIG["output_dir"]`.

For the public SpMap export, set the canonical label, primary split,
five-fold-validation, five selected evaluation-run, final `primary_1280`,
CONCH-feature, and output paths in `03_prepare_public_spmap_artifacts.py` and run
it in `Spatial_py`. The script writes 61,628 reference labels, the combined
primary/five-fold partition table, five-seed performance, the seed-24 confusion
matrix, 25 fold OOF prediction files, the five final `primary_1280` fold
checkpoints with portable manifests, and the canonical feature shards.

## Final SpMap model contract

The final SpMap model is the executed `primary_1280` configuration. It
uses paired 768- and 512-dimensional CONCH representations (`768 + 512 =
1280` input) in a `1280 -> 450 -> 4` MLP with ReLU, dropout 0.25, and no
BatchNorm. Checkpoints are selected by validation accuracy. The public
final-model archive contains exactly five `primary_1280` fold checkpoints, one
per fold, with portable manifests.

## Supplementary Figure 11 contracts

- S11A contains five fold-validation and five common-internal-holdout model
  evaluations. Bars show mean plus/minus SD and points show individual model
  evaluations; all common-holdout evaluations use the same reserved tiles.
- S11B summarizes the 49,001-tile pooled out-of-fold confusion table. Each
  cell reports its count and true-class row percentage.
- S11C reports precision, recall, one-versus-rest AUROC, one-versus-rest
  AUPRC, and macro-F1 for TC, PIR, PSM, and OTHER_PT. Two-sided 95% t intervals
  use training seeds 24, 101, 202, 303, and 404.
- S11D performs a one-to-one inner join of the discovery-cohort source tables. IMC
  ISR is `(n_bdme_i + 1) / (n_bdme_s + 1)`, clipped to `[0.1, 10]`; H&E ISR is
  `bile_1_count / bile_2_count`. The panel reports Pearson correlation, its
  two-sided P value, and a Fisher-z 95% CI for 33 matched patients.

Shared reference-tile preparation, feature extraction, training, inference,
and WSI-level ISR aggregation are provided by `src/spmap`.

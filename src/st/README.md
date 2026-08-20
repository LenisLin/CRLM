# Shared Stereo-seq Functions

`st_functions.R` provides reusable in-memory R operations for processed FDZS-4
Bin100 spatial-transcriptomics data, subtype-level RCTD reference preparation,
normalized RCTD-weight attachment, canonical SpMap tile-to-spot mapping,
spot-level composition aggregation, and GSEA result formatting. The spatial
analysis unit is one 50 x 50 um^2 Bin100 spot; the RCTD reference unit is one
annotated scRNA-seq cell.

## Function Reference

| Function | Required input | Return contract |
|---|---|---|
| `read_anndata_to_spe(data_path, verbose = TRUE)` | A directory containing the four processed MatrixMarket files listed below | A counts-only `SpatialExperiment` with genes x spots counts, `x`/`y` coordinates, source metadata, and `sample_id = core_name` |
| `balance_cell_types(sc_data, cell_type_col, k = 100, min_cells = 10, seed = 619)` | Annotated scRNA-seq raw counts with `cell_type_col` in `colData(sc_data)` | A list with balanced raw-count `sc_data` and sampling `stats`; every retained subtype has `k` columns |
| `assign_predictions_to_spots(spatial_df, tile_anno)` | Registered spot and tile coordinates with all required sample IDs | One integer SpMap class per input spot in input-row order |
| `combine_deconvolution(spe, result_df, deconv_method = "RCTD")` | A labelled ST object and normalized weights indexed by unique spot ID | The spot-ID intersection in ST-object order with one `<deconv_method>_<subtype>` column per weight |
| `aggregate_celltype_proportions(col_data, mapping)` | Spot metadata with `sample_id` and a named list of subtype-column groups | A data frame with `sample_id` and one row-sum column per mapping entry |
| `run_gsea_analysis(gene_ranks, pathways, database)` | Named gene ranks and a preloaded versioned pathway list | A complete `fgsea` result table with `database` and cleaned pathway names |

## Processed Bin100 Input

`read_anndata_to_spe` reads these files from one directory:

- `expression_profile.mtx`: spots x genes before import; transposed to genes x spots.
- `metadata.csv`: spot-indexed rows with `core_name`, `x`, and `y` columns.
- `row_names.csv`: gene names in matrix-column order; the first column supplies names.
- `column_names.csv`: spot IDs in matrix-row order; the first column supplies names.

The returned counts assay has gene row names and spot column names;
`spatialCoords` contains `x` and `y`; and `colData` retains the source metadata
with `sample_id = core_name`. Input construction requires all four files,
aligned dimensions and spot-ID order, the listed metadata fields, and nonempty
gene names.

## RCTD Reference And Weights

`balance_cell_types` uses `k = 100`, `min_cells = 10`, and `seed = 619` by
default. Each subtype with at least 10 cells contributes exactly 100 sampled
columns: sampling is without replacement at 100 or more available cells and
with replacement below 100 cells. Counts remain raw, duplicated sampled cells
receive unique column names, and the input supplies `cell_type_col` with at
least one subtype meeting `min_cells`.

`combine_deconvolution` matches normalized weights to ST spots by canonicalized
spot ID, preserves ST-object order, and appends `RCTD_<subtype>` columns by
default. The SpMap-labelled mapping contains 18,186 spots, and the normalized
RCTD table contains 17,902 uniquely indexed matched spots, including
377 PIR-Niche and 1,522 PSM-Niche spots.

## SpMap Tile Mapping

Spot and tile coordinates use one StereoMap-registered coordinate system.
`spatial_df` requires `x`, `y`, and `sample_id`; `tile_anno` requires
`sample_id`, `tile_left`, `tile_right`, `tile_top`, `tile_bottom`, and
`pred_class`. Mapping operates within each sample. Inclusive tile bounds select
the first containing tile in source-row order. Remaining spots receive the class
of the nearest tile measured from (`tile_left`, `tile_top`). A spot sample
requires at least one corresponding tile annotation. Class IDs are `0=TC`,
`1=PIR-Niche`, `2=PSM-Niche`, and `3=other PT`.

This mapping yields 463 TC, 382 PIR-Niche, 1,526 PSM-Niche, and 15,815 other-PT
spots, with 1,319 nearest-tile assignments.

## GSEA Parameters And Packages

`run_gsea_analysis` accepts caller-preloaded pathway collections and runs
`fgsea` with `minSize = 3`, `maxSize = 500`, and `eps = 0`. It returns every
tested pathway with the supplied database label and a lowercase,
underscore-to-space `pathway_clean` field. Key package versions are listed as
`st_rctd`, `st_spatial`, and `st_gsea` in `envs/package_versions.yml`.

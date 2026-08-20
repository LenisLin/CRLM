# FDZS-3 single-cell RNA-seq shared source

This directory provides the shared FDZS-3 operations for processed MatrixMarket
loading, Harmony integration, cell annotation, inferCNV, R-reference export,
and CellChat inference. The Figure 7 entry point supplies figure-specific
cluster labels, output locations, and plotting.

## Processed MatrixMarket input

The workflow starts from the processed six-sample FDZS-3 MatrixMarket bundle:

| File | Contents | Orientation or required fields |
|---|---|---|
| `processed_matrix.mtx` | count matrix | genes by cells |
| `processed_features.tsv` | gene annotation | `gene_id`, `gene_symbol`, `gene_type` |
| `processed_barcodes.tsv` | cell identifiers | one unique barcode per matrix column |
| `processed_metadata.csv` | cell metadata | unique `cell_barcode` index with exact barcode coverage |

`load_processed_matrixmarket()` transposes the matrix into a cell-by-gene
`AnnData` object, preserves the stable gene ID in `adata.var["gene_id"]`,
makes display gene symbols unique, and stores the count matrix in `adata.raw`.
The Figure 7 entry point requires `patient`, `doublet`, and `sample_id` in
`adata.obs`, then creates `batch` directly from `patient`.

## Pipeline order and parameters

1. **Integration.** `normalize_log1p()` normalizes each cell to 10,000 counts
   and applies `log1p`. `prepare_harmony_input()` selects 2,000 highly variable
   genes with `batch` as the batch field, scales values to a maximum of 10, and
   calculates 50 principal components. `integrate_harmony()` writes
   `X_pca_harmony` using `batch`.
2. **Major-cell annotation.** `remove_doublets()` retains cells whose
   `doublet` value differs from `Doublet`. `cluster_major_cells()` builds a
   10-neighbor graph from the first 40 `X_pca_harmony` components, computes a
   UMAP, and applies Leiden clustering at resolution 0.8. Marker ranking uses
   Wilcoxon tests, and the Figure 7 entry point writes `Major_type` from its
   explicit cluster-label map.
3. **Lineage annotation.** `preprocess_lineage()` creates B/plasma, T/NK, and
   myeloid subsets from `Major_type`; each subset is normalized to 10,000,
   log-transformed, processed with 2,000 highly variable genes, scaled,
   projected by PCA, and graphed with 10 neighbors over 50 inherited
   `X_pca_harmony` components. B/plasma uses Leiden resolution 0.8; T/NK and
   myeloid use resolution 1.0. The Figure 7 entry point writes `Sub_type` from
   its explicit lineage cluster-label maps.
4. **inferCNV and epithelial annotation.** `load_gene_gtf()` and
   `add_gtf_annotations()` align GRCh38-2024-A coordinates through stable
   Ensembl gene IDs. `run_infercnv()` uses `Major_type` reference populations
   `T`, `Myeloid`, `NK`, `B`, `Plasma`, and `Hepatocyte` with a 250-gene window.
   `prepare_epithelial_subset()` applies the FDZS3-P006 sample filter through
   `excluded_samples=("FDZS3-P006",)`, selects 2,000 highly variable genes,
   integrates on `batch`, constructs a 10-neighbor graph from 50 Harmony
   components, and uses expression Leiden resolution 0.8. `apply_cnv_status()`
   assigns cholangiocyte status to CNV clusters `5`, `25`, `17`, `26`, `20`,
   `28`, and `29`.
5. **R-reference export.** `export_anndata_for_r()` writes the final annotated
   count reference directly from `adata.raw`. The Figure 7 entry point exports
   cells after the FDZS3-P006 epithelial filter and retains cells with assigned
   `Major_type` values.
6. **CellChat.** `read_anndata_to_seurat()` reads the R reference bundle.
   `filter_cellchat_subtypes()` retains the 12 FDZS-3 cell populations,
   `assign_cholangiocyte_states()` scores cholangiocytes with `CA9`, `SLC2A1`,
   `LDHA`, `PGK1`, `ENO1`, and `ALDOA`, and assigns the lower and upper 30%
   quantiles to `Chol_HypoxiaLow` and `Chol_HypoxiaHigh`. `compute_cellchat()`
   uses `CellChatDB.human`, `raw.use = TRUE`, `min.cells = 10`, pathway
   probability calculation, and network aggregation.

## Required fields and analysis units

| Stage | Required and generated fields | Analysis unit |
|---|---|---|
| Integration | input: `cell_barcode`, `patient`, `doublet`, `sample_id`; generated: `batch` | cell; `patient` supplies the integration batch |
| Major annotation | input: `doublet`, `X_pca_harmony`; generated: `Major_type` | cell |
| Lineage annotation | input: `Major_type`, `X_pca_harmony`; generated: `Sub_type` | cell |
| inferCNV and epithelial annotation | `Major_type`, `Sub_type`, `sample_id`, `batch`, `gene_id` | cell |
| R export | `adata.raw`, `adata.obs_names`, `adata.raw.var_names`, `adata.obs` | cell |
| CellChat | `Sub_type`, `Sub_type_Revised` | cell-state/type group |

CellChat network values summarize inferred communication probabilities between
cell-state/type groups from cell-level normalized expression.

## R interoperability contract

`interop.py` writes a cell-by-gene MatrixMarket bundle for the R stages:

| File | Contents |
|---|---|
| `expression_profile.mtx` | cells by genes count matrix |
| `row_names.csv` | gene names in matrix-column order |
| `column_names.csv` | cell names in matrix-row order |
| `metadata.csv` | one metadata row per cell in matrix-row order |
| `<embedding>.csv` | optional cell-by-component embedding |

`cellchat.R` transposes `expression_profile.mtx` to genes by cells, assigns
the listed gene and cell names, and creates the Seurat object with aligned
metadata rows.

## Outputs and package reference

The Figure 7 entry point writes `integrated_scRNA_stereoseq.h5ad`,
`major_anno_all.h5ad`, `final_annotated_scRNA.h5ad`, `scrna_run_summary.tsv`,
and the `r_bundle/` MatrixMarket reference. The Supplementary Figure 16
entry point writes the CellChat object, state summary, parameter table,
interaction tables, group activity table, and figure files.

The Python stages use the `scrna` analysis family and the R CellChat stage uses
the `cellchat` analysis family. Their observed runtimes and key package
versions are listed in `envs/package_versions.yml`.

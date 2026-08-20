"""Load and integrate processed FDZS-3 scRNA MatrixMarket inputs.

Purpose:
    Reconstruct the processed cell-by-gene AnnData object and run the
    normalization, feature-selection, PCA, and Harmony steps for Figure 7.
Callers:
    ``figures/Figure7/01_integrate_annotate_scrna.py``.
Inputs:
    The processed MatrixMarket matrix, features, barcodes, and metadata files.
Outputs:
    An AnnData object updated in place with normalized expression, PCA, and
    Harmony coordinates.
Ordered use:
    Call ``load_processed_matrixmarket``, ``normalize_log1p``,
    ``prepare_harmony_input``, and ``integrate_harmony`` in that order.
"""

from pathlib import Path
from typing import Union

import pandas as pd
import scanpy as sc


def _find_matrixmarket_file(data_dir: Path, name: str) -> Path:
    """Locate one uncompressed or gzip-compressed MatrixMarket bundle file.

    Args:
        data_dir: Directory containing the processed MatrixMarket bundle.
        name: Uncompressed filename to resolve.

    Returns:
        Path to ``name`` when present, otherwise to ``name`` with ``.gz``
        appended.

    Raises:
        FileNotFoundError: Neither uncompressed nor gzip-compressed file is
            present in ``data_dir``.
    """
    # Prefer the compressed bundle member when both representations are present.
    compressed = data_dir / f"{name}.gz"
    if compressed.exists():
        return compressed

    uncompressed = data_dir / name
    if uncompressed.exists():
        return uncompressed

    raise FileNotFoundError(f"Neither {compressed} nor {uncompressed} exists")


def load_processed_matrixmarket(
    data_dir: Union[str, Path],
    *,
    prefix: str = "processed",
):
    """Load the processed MatrixMarket quartet as a cell-by-gene AnnData object.

    Args:
        data_dir: Directory containing the prefixed matrix, feature, barcode,
            and metadata files.
        prefix: Shared filename prefix for the processed input quartet.

    Returns:
        AnnData object with cell-by-gene counts in ``X``, stable gene IDs and
        gene types in ``var``, observation metadata in ``obs``, and the
        unmodified counts assigned to ``raw``.

    Raises:
        FileNotFoundError: A required MatrixMarket bundle file is absent.
        ValueError: Barcode or metadata identifiers are nonunique or do not
            match exactly.
    """
    # Resolve the three MatrixMarket companions from the same named bundle.
    data_dir = Path(data_dir)
    matrix_path = _find_matrixmarket_file(data_dir, f"{prefix}_matrix.mtx")
    features_path = _find_matrixmarket_file(data_dir, f"{prefix}_features.tsv")
    barcodes_path = _find_matrixmarket_file(data_dir, f"{prefix}_barcodes.tsv")
    metadata_path = data_dir / f"{prefix}_metadata.csv"

    # MatrixMarket is gene-by-cell on disk; AnnData stores cells by genes.
    adata = sc.read_mtx(matrix_path).T
    features = pd.read_csv(
        features_path,
        sep="\t",
        header=None,
        names=["gene_id", "gene_symbol", "gene_type"],
    )
    barcodes = pd.read_csv(
        barcodes_path,
        sep="\t",
        header=None,
        names=["barcode"],
    )
    metadata = pd.read_csv(metadata_path, index_col="cell_barcode")

    # Establish a one-to-one barcode key before transferring observation metadata.
    barcode_index = pd.Index(barcodes["barcode"])
    if not barcode_index.is_unique:
        raise ValueError("processed barcodes must be unique")
    if not metadata.index.is_unique:
        raise ValueError("processed metadata cell_barcode values must be unique")

    missing_barcodes = barcode_index.difference(metadata.index)
    unexpected_barcodes = metadata.index.difference(barcode_index)
    if len(missing_barcodes) or len(unexpected_barcodes):
        raise ValueError(
            "processed metadata barcodes must exactly match processed barcodes; "
            f"missing={missing_barcodes.tolist()}, "
            f"unexpected={unexpected_barcodes.tolist()}"
        )

    # Attach feature fields and cell names in their respective matrix-axis order.
    adata.var_names = features["gene_symbol"].to_numpy()
    adata.var["gene_id"] = features["gene_id"].to_numpy()
    adata.var["gene_type"] = features["gene_type"].to_numpy()
    adata.obs_names = barcodes["barcode"].to_numpy()

    # Reindex metadata to the matrix barcode order rather than CSV row order.
    aligned_metadata = metadata.reindex(adata.obs_names)
    for column in aligned_metadata.columns:
        adata.obs[column] = aligned_metadata[column].to_numpy()

    # Preserve imported counts before subsequent normalization and HVG subsetting.
    adata.var_names_make_unique()
    adata.raw = adata
    return adata


def normalize_log1p(adata, *, target_sum: float = 1e4):
    """Store raw counts, normalize to the target sum, and apply log1p in place.

    Args:
        adata: AnnData object with count data in ``X``.
        target_sum: Per-cell total used by Scanpy normalization.

    Returns:
        The same AnnData object after normalization and log transformation.
    """
    # Capture counts once, then normalize each cell to a common library size.
    if adata.raw is None:
        adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata


def prepare_harmony_input(
    adata,
    *,
    batch_key: str = "batch",
    n_top_genes: int = 2000,
    scale_max_value: float = 10,
    n_comps: int = 50,
):
    """Select batch-aware HVGs, scale expression, and calculate PCA.

    Args:
        adata: Normalized and log-transformed AnnData object.
        batch_key: Observation column identifying batches for HVG selection.
        n_top_genes: Number of highly variable genes retained.
        scale_max_value: Maximum scaled value passed to Scanpy.
        n_comps: Number of principal components to calculate.

    Returns:
        The same AnnData object after HVG subsetting, scaling, and PCA.
    """
    # Select features across batches before scaling and projecting the shared space.
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        batch_key=batch_key,
        subset=True,
    )
    sc.pp.scale(adata, max_value=scale_max_value)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=n_comps)
    return adata


def integrate_harmony(
    adata,
    *,
    batch_key: str = "batch",
    basis: str = "X_pca",
    adjusted_basis: str = "X_pca_harmony",
):
    """Apply Harmony using the supplied observation-level batch field.

    Args:
        adata: PCA-prepared AnnData object.
        batch_key: Observation column specifying Harmony batch assignments.
        basis: ``obsm`` key containing PCA coordinates to integrate.
        adjusted_basis: ``obsm`` key receiving Harmony-adjusted coordinates.

    Returns:
        The same AnnData object with Harmony coordinates stored at
        ``adjusted_basis``.
    """
    # Correct the PCA coordinates by the observation-level batch assignments.
    sc.external.pp.harmony_integrate(
        adata,
        key=batch_key,
        basis=basis,
        adjusted_basis=adjusted_basis,
    )
    return adata

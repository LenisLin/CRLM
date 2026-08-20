"""Export annotated FDZS-3 scRNA data as an R-readable MatrixMarket bundle.

Purpose:
    Serialize the annotated AnnData reference used by the R CellChat stage.
Callers:
    ``figures/Figure7/01_integrate_annotate_scrna.py``.
Inputs:
    An AnnData object with ``raw`` expression data and a destination directory.
Outputs:
    A MatrixMarket expression matrix, cell and gene name files, cell metadata,
    and optionally one CSV file per observation-matrix entry.
Ordered use:
    Call after Figure 7 integration and annotation, then pass the output bundle
    to ``src/scrna/cellchat.R::read_anndata_to_seurat``.
"""

from pathlib import Path
from typing import Union

import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def export_anndata_for_r(
    adata, output_dir: Union[str, Path], *, save_obsm: bool = False
):
    """Write a cell-by-gene MatrixMarket bundle for the R analysis stages.

    Args:
        adata: Annotated AnnData object whose ``raw`` slot contains the
            cell-by-gene expression matrix and whose observations identify
            cells for metadata export.
        output_dir: Directory that will receive the R-readable bundle.
        save_obsm: Whether to export each entry in ``adata.obsm`` as a
            cell-indexed CSV file.

    Returns:
        None. Writes ``expression_profile.mtx``, ``metadata.csv``,
        ``row_names.csv``, and ``column_names.csv`` to ``output_dir``; when
        requested, also writes one CSV file for each observation-matrix entry.
    """
    if adata.raw is None:
        raise ValueError("adata.raw is required for the R reference export")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Export raw cells-by-genes counts; the R importer transposes this orientation.
    matrix = sparse.csr_matrix(adata.raw.X, dtype="float64")
    mmwrite(output_dir / "expression_profile.mtx", matrix)
    adata.obs.to_csv(output_dir / "metadata.csv")
    pd.Series(adata.raw.var_names).to_csv(
        output_dir / "row_names.csv", header=False, index=False
    )
    pd.Series(adata.obs_names).to_csv(
        output_dir / "column_names.csv", header=False, index=False
    )

    # Keep every optional embedding row aligned to the exported cell-name vector.
    if save_obsm:
        for key, values in adata.obsm.items():
            embedding_name = key[2:] if key.startswith("X_") else key
            pd.DataFrame(values, index=adata.obs_names).to_csv(
                output_dir / f"{embedding_name}.csv"
            )

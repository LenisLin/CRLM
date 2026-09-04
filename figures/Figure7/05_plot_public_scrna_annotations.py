"""Render Figure 7B and Supplementary Figure 15 from the public FDZS-3 H5AD.

This consumer uses the stored ``X_umap``, ``Major_type``, ``Sub_type``, and
``raw`` count layer. It does not rerun integration, clustering, CNV inference,
or annotation. Lineage panels are subsets of the stored global UMAP; the
lineage-specific UMAP coordinates from the private intermediate objects are not
part of the public-object contract.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


CONFIG = {
    "input_h5ad": Path("/path/to/FDZS3_scRNA_processed.h5ad"),
    "output_dir": Path("/path/to/figure7_public_scrna_annotations"),
}

MAJOR_MARKERS = {
    "Malignant": ["EPCAM", "CDH1", "KRT20", "KRT19", "KRT7"],
    "Cholangiocytes": ["SOX9"],
    "Hepatocytes": ["ALB"],
    "T": ["CD3E"],
    "NKT": ["GNLY", "NKG7"],
    "B": ["CD19", "MS4A1"],
    "Plasma": ["IGKC"],
    "Macrophage": ["CD68", "LYZ"],
    "Monocyte": ["CD14"],
    "cDC": ["ITGAX", "S100A9"],
    "Stromal": ["COL1A2", "FAP", "POSTN"],
    "Endothelial": ["PECAM1", "VWF"],
}

B_PLASMA_MARKERS = {
    "Pan_B_Plasma": ["CD19", "MS4A1", "IGKC"],
    "Memory": ["AIM2"],
    "HSP": ["HSPA1A", "HSP90AA1", "HSPE1"],
    "Activate": ["CD27"],
    "IgA": ["IGHA1", "IGHA2"],
    "IgG": ["IGHG1", "IGHG2", "IGHG3", "IGHG4"],
    "IgM": ["IGHM"],
}

NKT_MARKERS = {
    "Pan_T": ["CD3E", "CD8A", "CD8B", "CD4"],
    "NK": ["NKG7", "KLRD1", "FCGR3A", "GNLY"],
    "Naive": ["SELL", "CCR7", "TCF7"],
    "Exhausted": ["CTLA4", "HAVCR2", "TIGIT", "PDCD1"],
    "HSP": ["HSPA6", "HSPA1A", "HSP90AA1"],
    "CD4T_IL7R": ["IL7R"],
    "Treg": ["FOXP3", "TIGIT", "CTLA4", "IL2RA"],
    "CD8_GZMK": ["GZMK", "NKG7"],
    "CD8_Effector": ["FGFBP2", "GZMB", "PRF1"],
    "CD8_MKI67": ["MKI67"],
    "CD8_CXCL13": ["CXCL13"],
    "MAIT": ["SLC4A10", "DPP4", "KLRB1"],
    "GammaDelta_T": ["TRGC1", "TRGC2"],
}

MYELOID_MARKERS = {
    "Macrophage": ["CD68"],
    "Monocyte": ["CD14"],
    "Macro_APOE": ["APOE"],
    "Macro_FOLR2": ["FOLR2"],
    "Macro_SPP1": ["SPP1", "ITGAM"],
    "Macro_CD169": ["SIGLEC1"],
    "Macro_CD163": ["CD163"],
    "Kupffer": ["CLEC4F", "VSIG4"],
    "Neutrophil": ["FCGR3B", "CXCR2"],
    "cDC_CD1C": ["CLEC10A", "CD1C"],
    "cDC_CLEC9A": ["CLEC9A", "XCR1"],
    "pDC_LILRA4": ["LILRA4"],
}

REQUIRED_OBS_COLUMNS = ("Major_type", "Sub_type")
REQUIRED_OBSM_KEY = "X_umap"
REQUIRED_LAYER_KEYS = ()


def marker_names(marker_groups: dict[str, list[str]]) -> list[str]:
    """Return marker names once, retaining their configured display order."""
    return list(dict.fromkeys(gene for genes in marker_groups.values() for gene in genes))


def require_public_plot_contract(adata) -> None:
    """Raise when the public object lacks fields required by this consumer."""
    missing_obs = [name for name in REQUIRED_OBS_COLUMNS if name not in adata.obs]
    if missing_obs:
        raise ValueError("Missing required obs field(s): " + ", ".join(missing_obs))
    if REQUIRED_OBSM_KEY not in adata.obsm:
        raise ValueError("Missing required embedding: " + REQUIRED_OBSM_KEY)
    if adata.obsm[REQUIRED_OBSM_KEY].shape != (adata.n_obs, 2):
        raise ValueError("X_umap must contain one two-dimensional coordinate per cell")
    if adata.raw is None:
        raise ValueError("Missing required raw count layer")
    missing_layers = [name for name in REQUIRED_LAYER_KEYS if name not in adata.layers]
    if missing_layers:
        raise ValueError("Missing required layer(s): " + ", ".join(missing_layers))

    required_markers = set().union(
        marker_names(MAJOR_MARKERS),
        marker_names(B_PLASMA_MARKERS),
        marker_names(NKT_MARKERS),
        marker_names(MYELOID_MARKERS),
    )
    missing_markers = sorted(required_markers.difference(adata.raw.var_names))
    if missing_markers:
        raise ValueError("Missing required marker gene(s): " + ", ".join(missing_markers))


def ordered_categories(adata, obs_key: str) -> list[str]:
    """Return the stored categorical order, or observed values for plain strings."""
    values = adata.obs[obs_key]
    if hasattr(values, "cat"):
        return [str(value) for value in values.cat.categories]
    return sorted(str(value) for value in values.dropna().unique())


def stored_colors(adata, obs_key: str, categories: list[str]) -> dict[str, str]:
    """Map stored AnnData category colors to their matching labels."""
    colors = adata.uns.get(f"{obs_key}_colors")
    if colors is None or len(colors) != len(categories):
        return {}
    return dict(zip(categories, colors))


def save_umap(adata, obs_key: str, output_path: Path, *, mask=None) -> None:
    """Render stored global UMAP coordinates, optionally restricted to a lineage."""
    import matplotlib.pyplot as plt

    coordinates = np.asarray(adata.obsm[REQUIRED_OBSM_KEY])
    if mask is None:
        mask = np.ones(adata.n_obs, dtype=bool)
    mask = np.asarray(mask, dtype=bool)
    labels = adata.obs[obs_key].astype(str).to_numpy()
    categories = ordered_categories(adata, obs_key)
    colors = stored_colors(adata, obs_key, categories)

    figure, axis = plt.subplots(figsize=(5.5, 4.5), dpi=300)
    for category in categories:
        category_mask = mask & (labels == category)
        if category_mask.any():
            axis.scatter(
                coordinates[category_mask, 0],
                coordinates[category_mask, 1],
                s=1.5,
                alpha=0.8,
                linewidths=0,
                color=colors.get(category),
                label=category,
            )
    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    axis.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7)
    figure.tight_layout()
    figure.savefig(output_path, format="pdf", dpi=300, bbox_inches="tight")
    plt.close(figure)


def save_marker_dotplot(adata, obs_key: str, marker_groups: dict[str, list[str]], output_path: Path, *, mask) -> None:
    """Render raw-count marker summaries after per-cell library normalization."""
    import matplotlib.pyplot as plt

    mask = np.asarray(mask, dtype=bool)
    labels = adata.obs[obs_key].astype(str).to_numpy()[mask]
    categories = [category for category in ordered_categories(adata, obs_key) if (labels == category).any()]
    genes = marker_names(marker_groups)
    gene_indices = adata.raw.var_names.get_indexer(genes)
    raw_counts = adata.raw.X[:, gene_indices].tocsr()
    library_sizes = np.asarray(adata.raw.X.sum(axis=1)).ravel()
    scale = 1e4 / np.maximum(library_sizes, 1)
    expression = raw_counts.multiply(scale[:, None]).tocsr()
    expression.data = np.log1p(expression.data)
    expression = expression[mask]

    mean_expression = np.empty((len(categories), len(genes)))
    fraction_expressing = np.empty_like(mean_expression)
    for index, category in enumerate(categories):
        group_expression = expression[labels == category]
        mean_expression[index] = np.asarray(group_expression.mean(axis=0)).ravel()
        fraction_expressing[index] = group_expression.getnnz(axis=0) / group_expression.shape[0]

    minima = mean_expression.min(axis=0)
    spans = mean_expression.max(axis=0) - minima
    scaled_expression = np.divide(
        mean_expression - minima,
        spans,
        out=np.zeros_like(mean_expression),
        where=spans > 0,
    )

    figure_width = max(6.0, 0.32 * len(genes))
    figure_height = max(3.0, 0.28 * len(categories) + 1.5)
    figure, axis = plt.subplots(figsize=(figure_width, figure_height), dpi=300)
    for row in range(len(categories)):
        axis.scatter(
            np.arange(len(genes)),
            np.full(len(genes), row),
            s=16 + 130 * fraction_expressing[row],
            c=scaled_expression[row],
            cmap="Reds",
            vmin=0,
            vmax=1,
            linewidths=0,
        )
    axis.set_xticks(np.arange(len(genes)))
    axis.set_xticklabels(genes, rotation=90, fontsize=7)
    axis.set_yticks(np.arange(len(categories)))
    axis.set_yticklabels(categories, fontsize=7)
    axis.invert_yaxis()
    axis.set_xlim(-0.7, len(genes) - 0.3)
    axis.set_xlabel("Canonical marker genes")
    axis.set_ylabel(obs_key)
    figure.tight_layout()
    figure.savefig(output_path, format="pdf", dpi=300, bbox_inches="tight")
    plt.close(figure)


def main(config: dict = CONFIG) -> None:
    """Write the public-H5AD Figure 7B and Supplementary Figure 15 panel files."""
    import anndata as ad

    input_h5ad = Path(config["input_h5ad"])
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(input_h5ad)
    require_public_plot_contract(adata)

    major_mask = np.ones(adata.n_obs, dtype=bool)
    myeloid_mask = adata.obs["Major_type"].astype(str).to_numpy() == "Myeloid"
    nkt_mask = adata.obs["Major_type"].astype(str).isin(("T", "NK")).to_numpy()
    b_plasma_mask = adata.obs["Major_type"].astype(str).isin(("B", "Plasma")).to_numpy()

    save_umap(adata, "Major_type", output_dir / "Figure7B_annotation_umap.pdf", mask=major_mask)
    save_umap(adata, "Major_type", output_dir / "FigureS15_major_umap.pdf", mask=major_mask)
    save_marker_dotplot(adata, "Major_type", MAJOR_MARKERS, output_dir / "FigureS15_major_marker_dotplot.pdf", mask=major_mask)
    save_umap(adata, "Sub_type", output_dir / "FigureS15_myeloid_umap.pdf", mask=myeloid_mask)
    save_marker_dotplot(adata, "Sub_type", MYELOID_MARKERS, output_dir / "FigureS15_myeloid_marker_dotplot.pdf", mask=myeloid_mask)
    save_umap(adata, "Sub_type", output_dir / "FigureS15_nkt_umap.pdf", mask=nkt_mask)
    save_marker_dotplot(adata, "Sub_type", NKT_MARKERS, output_dir / "FigureS15_nkt_marker_dotplot.pdf", mask=nkt_mask)
    save_umap(adata, "Sub_type", output_dir / "FigureS15_b_plasma_umap.pdf", mask=b_plasma_mask)
    save_marker_dotplot(adata, "Sub_type", B_PLASMA_MARKERS, output_dir / "FigureS15_b_plasma_marker_dotplot.pdf", mask=b_plasma_mask)


if __name__ == "__main__":
    main()

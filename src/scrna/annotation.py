"""Annotate FDZS-3 scRNA cells and prepare Figure 7 lineage analyses.

Purpose:
    Define markers and label mappings, cluster major and lineage-restricted
    populations, and prepare epithelial inferCNV inputs for Figure 7.
Callers:
    ``figures/Figure7/01_integrate_annotate_scrna.py``.
Inputs:
    Harmony-integrated AnnData objects, curated cluster mappings, and a GTF
    gene-coordinate table for inferCNV.
Outputs:
    AnnData objects updated in place with clusters, labels, embeddings, and
    inferCNV fields; a parsed GTF coordinate DataFrame.
Ordered use:
    Remove doublets, cluster and label major cells, analyze lineage subsets,
    then add GTF annotations, run inferCNV, and classify epithelial cells.
"""

from __future__ import annotations

import gzip
import re

import anndata as ad
import infercnvpy as cnv
import pandas as pd
import scanpy as sc


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

CURRENT_MAJOR_CLUSTER_LABELS = {
    str(index): label
    for index, label in enumerate(
        [
            "T", "Epithelial", "T", "NK", "Epithelial", "Myeloid", "B",
            "Epithelial", "T", "T", "T", "Myeloid", "Myeloid",
            "Endothelial", "Epithelial", "Epithelial", "Fibroblast", "T",
            "Plasma", "T", "Hepatocyte", "B",
        ]
    )
}

CURRENT_SUBTYPE_CLUSTER_LABELS = {
    "B_Plasma": {
        str(index): label
        for index, label in enumerate(
            [
                "B_AIM2", "B_HSP", "B_IGHM", "B_IGHM", "Plasma", "B_HSP",
                "Unknown", "B_HSP", "T", "pDC_LILRA4", "B_CD27",
            ]
        )
    },
    "NKT": {
        str(index): label
        for index, label in enumerate(
            [
                "CD4T_IL7R", "CD8T_GZMK", "NK_FCGR3A", "MAIT_KLRB1",
                "T_TCF7", "CD4T_IL7R", "Treg_FOXP3", "NK_KLRB1",
                "CD8T_CXCL13", "T_HSP", "T_MKI67", "B_IGHM", "CD4T_IL7R",
            ]
        )
    },
    "Myeloid": {
        str(index): label
        for index, label in enumerate(
            [
                "cDC_CD1C", "Monocyte", "Macro_CD169", "Macro_SPP1",
                "Macro_CD163", "Macro_CD163", "Myeloid_Other", "Monocyte",
                "Macro_SPP1", "Myeloid_Other", "cDC_CLEC9A", "Macro_CD169",
                "Neutrophil", "pDC_LILRA4",
            ]
        )
    },
}


def remove_doublets(
    adata,
    *,
    doublet_key: str = "doublet",
    doublet_label: str = "Doublet",
):
    """Return the singlet cells for major-cell annotation.

    Args:
        adata: AnnData object containing a doublet classification column.
        doublet_key: Observation column holding doublet classifications.
        doublet_label: Value in ``doublet_key`` identifying doublets.

    Returns:
        Copy of ``adata`` restricted to cells not labeled ``doublet_label``.
    """
    # Remove doublet calls before graph construction and curated annotation.
    return adata[adata.obs[doublet_key] != doublet_label].copy()


def cluster_major_cells(
    adata,
    *,
    resolution: float = 0.8,
    n_neighbors: int = 10,
    n_pcs: int = 40,
    use_rep: str = "X_pca_harmony",
    cluster_key: str = "leiden",
):
    """Build the major-cell graph, UMAP, and Leiden clustering.

    Args:
        adata: Harmony-integrated AnnData object.
        resolution: Leiden clustering resolution.
        n_neighbors: Number of neighbors in the graph.
        n_pcs: Number of components used when constructing neighbors.
        use_rep: ``obsm`` key containing the representation for neighbors.
        cluster_key: Observation column receiving Leiden cluster labels.

    Returns:
        The same AnnData object with neighbor graph, UMAP, and cluster labels.
    """
    # Derive UMAP and Leiden labels from the same Harmony-based neighbor graph.
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key)
    return adata


def rank_cluster_markers(adata, *, cluster_key: str = "leiden"):
    """Rank cluster markers with Wilcoxon tests.

    Args:
        adata: AnnData object containing cluster labels.
        cluster_key: Observation column defining groups for marker ranking.

    Returns:
        The same AnnData object with ranked genes stored in ``uns``.
    """
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon")
    return adata


def apply_cluster_labels(
    adata,
    cluster_to_label: dict[str, str],
    *,
    cluster_key: str = "leiden",
    output_key: str = "Major_type",
):
    """Apply the supplied cluster-to-label mapping.

    Args:
        adata: AnnData object containing cluster labels.
        cluster_to_label: Mapping from string cluster IDs to annotation labels.
        cluster_key: Observation column holding cluster IDs.
        output_key: Observation column receiving mapped labels.

    Returns:
        The same AnnData object with mapped labels in ``output_key``.
    """
    # Match cluster IDs as strings to the curated annotation mapping.
    adata.obs[output_key] = adata.obs[cluster_key].astype(str).map(cluster_to_label)
    return adata


def preprocess_lineage(
    adata,
    major_types: tuple[str, ...],
    *,
    resolution: float,
    major_key: str = "Major_type",
    n_top_genes: int = 2000,
    n_neighbors: int = 10,
    n_pcs: int = 50,
    use_rep: str = "X_pca_harmony",
    cluster_key: str = "leiden",
):
    """Prepare a lineage subset with inherited Harmony coordinates.

    Args:
        adata: Major-cell annotated AnnData object with ``raw`` counts and
            Harmony coordinates.
        major_types: Major-cell labels retained in the lineage subset.
        resolution: Leiden clustering resolution for the lineage subset.
        major_key: Observation column holding major-cell labels.
        n_top_genes: Number of highly variable genes retained.
        n_neighbors: Number of neighbors in the lineage graph.
        n_pcs: Number of components used when constructing neighbors.
        use_rep: ``obsm`` key containing inherited Harmony coordinates.
        cluster_key: Observation column receiving lineage cluster labels.

    Returns:
        New AnnData lineage subset with normalized expression, PCA, graph,
        UMAP, Leiden labels, and Wilcoxon marker rankings.

    Notes:
        The Figure 7 entry point uses resolution 0.8 for B/plasma and 1.0 for
        NKT and myeloid lineages.
    """
    # Restore raw counts for the selected lineage while retaining cell metadata
    # and the parent object's integrated coordinates.
    cell_mask = adata.obs[major_key].isin(major_types)
    subset = ad.AnnData(
        X=adata[cell_mask].raw.X,
        obs=adata[cell_mask].obs.copy(),
        var=adata.raw.var.copy(),
        obsm=adata[cell_mask].obsm.copy(),
    )
    # Recompute lineage-specific normalized expression, variable genes, and PCA.
    sc.pp.normalize_total(subset, target_sum=1e4)
    sc.pp.log1p(subset)
    sc.pp.highly_variable_genes(subset, n_top_genes=n_top_genes)
    sc.pp.scale(subset)
    sc.tl.pca(subset)
    # Cluster the lineage in the inherited Harmony space and rank its markers.
    sc.pp.neighbors(subset, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)
    sc.tl.umap(subset)
    sc.tl.leiden(subset, resolution=resolution, key_added=cluster_key)
    rank_cluster_markers(subset, cluster_key=cluster_key)
    return subset


def parse_gtf_attributes(attributes: str) -> dict[str, str]:
    """Parse quoted key-value fields from one GTF attributes column.

    Args:
        attributes: Ninth-column GTF attribute string.

    Returns:
        Mapping of parsed attribute keys to their quoted values.
    """
    # GTF attributes are semicolon-delimited quoted key-value fields.
    parsed = {}
    for item in attributes.strip().split(";"):
        match = re.match(r'^\s*(\w+)\s+"([^"]+)"', item.strip())
        if match:
            key, value = match.groups()
            parsed[key] = value
    return parsed


def _canonical_gene_ids(values) -> pd.Index:
    """Normalize Ensembl gene IDs by removing an optional version suffix.

    Args:
        values: Sequence of Ensembl-style gene identifiers.

    Returns:
        String Index containing gene identifiers with ``.version`` suffixes removed.
    """
    return pd.Index(values, dtype="string").str.replace(r"\.\d+$", "", regex=True)


def load_gene_gtf(gtf_path: str) -> pd.DataFrame:
    """Read gene coordinates from an uncompressed or gzip-compressed GTF.

    Args:
        gtf_path: Path to a GTF file, optionally gzip-compressed.

    Returns:
        DataFrame indexed by stable Ensembl gene ID with gene name, chromosome,
        start, end, and strand columns for gene features.

    Raises:
        ValueError: Stable Ensembl gene IDs are duplicated in the GTF.
    """
    # Retain only gene features needed for inferCNV genomic ordering.
    rows = []
    open_file = gzip.open if str(gtf_path).endswith(".gz") else open
    with open_file(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = parse_gtf_attributes(fields[8])
            rows.append(
                {
                    "gene_name": attributes.get("gene_name", "Unknown"),
                    "gene_id": attributes.get("gene_id", "Unknown"),
                    "chromosome": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                }
            )
    # Version-stripped Ensembl IDs provide the stable alignment key to AnnData.
    gtf = pd.DataFrame(rows)
    gtf["ensg"] = _canonical_gene_ids(gtf["gene_id"])
    if gtf["ensg"].duplicated().any():
        duplicated = gtf.loc[gtf["ensg"].duplicated(keep=False), "ensg"].unique()
        raise ValueError(f"GTF contains duplicate stable gene IDs: {duplicated.tolist()}")
    return gtf.set_index("ensg").drop(columns="gene_id")


def add_gtf_annotations(adata, gtf: pd.DataFrame):
    """Attach infercnvpy gene-coordinate fields to ``adata.var``.

    Args:
        adata: AnnData object with stable gene IDs in ``var['gene_id']``.
        gtf: Gene-coordinate DataFrame returned by ``load_gene_gtf``.

    Returns:
        The same AnnData object with inferCNV coordinate fields in ``var``.

    Raises:
        ValueError: ``adata.var`` does not contain a ``gene_id`` column.
    """
    if "gene_id" not in adata.var:
        raise ValueError("adata.var must contain stable gene IDs in 'gene_id'")
    # Reindex coordinates in AnnData gene order so every annotation row stays aligned.
    gene_ids = _canonical_gene_ids(adata.var["gene_id"])
    annotations = gtf.reindex(gene_ids)
    adata.var["ensg"] = gene_ids.fillna("Unknown").to_numpy()
    adata.var["chromosome"] = annotations["chromosome"].fillna("Unknown").to_numpy()
    adata.var["start"] = annotations["start"].fillna(0).astype(int).to_numpy()
    adata.var["end"] = annotations["end"].fillna(0).astype(int).to_numpy()
    adata.var["strand"] = annotations["strand"].fillna("Unknown").to_numpy()
    return adata


def run_infercnv(
    adata,
    reference_categories: tuple[str, ...],
    *,
    reference_key: str = "Major_type",
    window_size: int = 250,
):
    """Run inferCNV with the supplied reference populations.

    Args:
        adata: Gene-coordinate annotated AnnData object.
        reference_categories: Labels identifying reference cell populations.
        reference_key: Observation column containing reference labels.
        window_size: Number of genes in each inferCNV smoothing window.

    Returns:
        The same AnnData object with inferCNV results added by infercnvpy.
    """
    # Use the supplied non-epithelial populations as the expression reference.
    cnv.tl.infercnv(
        adata,
        reference_key=reference_key,
        reference_cat=list(reference_categories),
        window_size=window_size,
    )
    return adata


def prepare_epithelial_subset(
    adata,
    *,
    excluded_samples: tuple[str, ...],
    sample_key: str = "sample_id",
    batch_key: str = "batch",
    major_key: str = "Major_type",
    epithelial_label: str = "Epithelial",
    expression_resolution: float = 0.8,
):
    """Prepare the epithelial subset after applying the supplied sample filter.

    Args:
        adata: Major-cell annotated, inferCNV-processed AnnData object.
        excluded_samples: Sample IDs designated for the epithelial-sample filter.
        sample_key: Observation column containing sample IDs.
        batch_key: Observation column defining Harmony batches.
        major_key: Observation column containing major-cell labels.
        epithelial_label: Value in ``major_key`` identifying epithelial cells.
        expression_resolution: Leiden resolution for expression clustering.

    Returns:
        New epithelial AnnData subset with expression and CNV embeddings,
        clusters, and CNV scores.

    Notes:
        The Figure 7 entry point supplies ``excluded_samples=(\"FDZS3-P006\",)``.
    """
    # Apply the sample exclusion within epithelial cells, then restore raw counts
    # while carrying the full-object inferCNV representation into the subset.
    cell_mask = (adata.obs[major_key] == epithelial_label) & ~adata.obs[
        sample_key
    ].isin(excluded_samples)
    subset = ad.AnnData(
        X=adata[cell_mask].raw.X,
        obs=adata[cell_mask].obs.copy(),
        var=adata.raw.var.copy(),
        obsm=adata[cell_mask].obsm.copy(),
        uns=adata.uns.copy(),
    )
    # Build a batch-integrated expression space for epithelial clustering.
    sc.pp.normalize_total(subset, target_sum=1e4)
    sc.pp.log1p(subset)
    sc.pp.highly_variable_genes(subset, n_top_genes=2000)
    sc.pp.scale(subset)
    sc.tl.pca(subset)
    sc.external.pp.harmony_integrate(
        subset,
        key=batch_key,
        basis="X_pca",
        adjusted_basis="X_pca_harmony",
    )
    sc.pp.neighbors(subset, n_neighbors=10, n_pcs=50, use_rep="X_pca_harmony")
    sc.tl.umap(subset)
    sc.tl.leiden(subset, resolution=expression_resolution, key_added="leiden")

    # Independently embed and cluster the inherited inferCNV profiles, then score CNV.
    cnv.tl.pca(subset)
    cnv.pp.neighbors(subset)
    cnv.tl.leiden(subset)
    cnv.tl.umap(subset)
    cnv.tl.cnv_score(subset)
    return subset


def apply_cnv_status(
    epithelial_adata,
    cholangiocyte_clusters: tuple[str, ...],
    *,
    cluster_key: str = "cnv_leiden",
    output_key: str = "Sub_type",
):
    """Apply the supplied CNV cluster classification.

    Args:
        epithelial_adata: Epithelial AnnData object with CNV Leiden labels.
        cholangiocyte_clusters: CNV cluster IDs classified as cholangiocytes.
        cluster_key: Observation column holding CNV Leiden labels.
        output_key: Observation column receiving malignant or cholangiocyte
            subtype labels.

    Returns:
        The same epithelial AnnData object with subtype labels in ``output_key``.
    """
    # Initialize the epithelial lineage as malignant and relabel curated CNV clusters.
    epithelial_adata.obs[output_key] = "Malignant"
    epithelial_adata.obs.loc[
        epithelial_adata.obs[cluster_key].astype(str).isin(cholangiocyte_clusters),
        output_key,
    ] = "Cholangiocyte"
    return epithelial_adata

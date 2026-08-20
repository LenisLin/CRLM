"""Target: Figure 7B and Supplementary Figure 15 FDZS-3 scRNA outputs.

Purpose: integrate and annotate the FDZS-3 scRNA reference.
Inputs: editable paths and sample settings in ``CONFIG``.
Outputs: annotated AnnData files, R bundle, run summary, and panel PDFs.
Ordered workflow: load, integrate, annotate, infer CNV, export, summarize.
"""

# Target: Figure 7B and Supplementary Figure 15 FDZS-3 scRNA outputs.
# Purpose: integrate and annotate the FDZS-3 scRNA reference.
# Inputs: edit paths and sample settings in CONFIG.
# Outputs: annotated AnnData files, R bundle, run summary, and panel PDFs.
# Ordered workflow: load, integrate, annotate, infer CNV, export, summarize.

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Iterable, Tuple

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

EXPECTED_ANNOTATED_CELLS = 35_900

# Edit these paths for the Figure 7 analysis.
CONFIG = {
    "input_dir": Path("/path/to/fdzs3_processed_matrixmarket_bundle"),
    "gtf": Path("/path/to/genes.gtf"),
    "output_dir": Path("/path/to/figure7_scrna_output"),
    "prefix": "processed",
    "excluded_epithelial_samples": ("FDZS3-P006",),
}

MAJOR_COLORS = {
    "B": "#FFD700",
    "Plasma": "#FFFF99",
    "Epithelial": "#33A02C",
    "Hepatocyte": "#B2DF8A",
    "Fibroblast": "#FF7F00",
    "Endothelial": "#E31A1C",
    "Myeloid": "#6A3D9A",
    "NK": "#A6CEE3",
    "T": "#1F78B4",
    "Unknown": "#BDBDBD",
}

B_PLASMA_TO_MAJOR = {
    "B_AIM2": "B",
    "B_HSP": "B",
    "B_IGHM": "B",
    "B_CD27": "Unknown",
    "Plasma": "Plasma",
    "T": "T",
    "pDC_LILRA4": "Myeloid",
    "Unknown": "Unknown",
}

NKT_TO_MAJOR = {
    "T_HSP": "T",
    "T_MKI67": "T",
    "T_TCF7": "T",
    "CD4T_IL7R": "T",
    "Treg_FOXP3": "T",
    "CD8T_GZMK": "T",
    "CD8T_CXCL13": "T",
    "MAIT_KLRB1": "T",
    "NK_FCGR3A": "NK",
    "NK_KLRB1": "NK",
    "B_IGHM": "B",
}


def load_analysis_dependencies() -> None:
    """Purpose: load Figure 7 scRNA dependencies when analysis execution begins.

    Parameters: none.
    Returns: ``None`` after binding the analysis dependencies used by this module.
    """
    # Bind heavy analysis libraries only for an executed workflow; module import remains lightweight.
    global plt, pd, sc
    global B_PLASMA_MARKERS, CURRENT_MAJOR_CLUSTER_LABELS, CURRENT_SUBTYPE_CLUSTER_LABELS
    global MAJOR_MARKERS, MYELOID_MARKERS, NKT_MARKERS
    global add_gtf_annotations, apply_cluster_labels, apply_cnv_status
    global cluster_major_cells, load_gene_gtf, prepare_epithelial_subset
    global preprocess_lineage, rank_cluster_markers, remove_doublets, run_infercnv
    global integrate_harmony, load_processed_matrixmarket, normalize_log1p
    global prepare_harmony_input, export_anndata_for_r

    import matplotlib.pyplot as plt
    import pandas as pd
    import scanpy as sc

    from src.scrna.annotation import (
        B_PLASMA_MARKERS,
        CURRENT_MAJOR_CLUSTER_LABELS,
        CURRENT_SUBTYPE_CLUSTER_LABELS,
        MAJOR_MARKERS,
        MYELOID_MARKERS,
        NKT_MARKERS,
        add_gtf_annotations,
        apply_cluster_labels,
        apply_cnv_status,
        cluster_major_cells,
        load_gene_gtf,
        prepare_epithelial_subset,
        preprocess_lineage,
        rank_cluster_markers,
        remove_doublets,
        run_infercnv,
    )
    from src.scrna.integration import (
        integrate_harmony,
        load_processed_matrixmarket,
        normalize_log1p,
        prepare_harmony_input,
    )
    from src.scrna.interop import export_anndata_for_r


def require_obs_columns(adata, columns: Iterable[str]) -> None:
    """Purpose: verify required observation metadata before analysis.

    Parameters: ``adata`` is the scRNA AnnData object; ``columns`` are required names.
    Returns: ``None`` after validation, or raises ``ValueError`` for missing columns.
    """
    # Observation columns are cell-level metadata and must remain aligned to AnnData row identities.
    missing = [column for column in columns if column not in adata.obs.columns]
    if missing:
        raise ValueError("Missing required scRNA metadata column(s): " + ", ".join(missing))


def save_umap(adata, color: str, output_path: Path, *, palette=None) -> None:
    """Purpose: render a UMAP panel with the Figure 7 settings.

    Parameters: ``adata`` is annotated AnnData; ``color`` is the observation field;
    ``output_path`` is the PDF destination; ``palette`` optionally defines colors.
    Returns: ``None`` after writing the PDF.
    """
    # Render one point per cell using the existing embedding coordinates.
    figure, axis = plt.subplots(figsize=(4, 4), dpi=300)
    sc.pl.umap(
        adata,
        color=color,
        palette=palette,
        size=1.5,
        alpha=0.8,
        frameon=True,
        legend_loc="right margin",
        legend_fontsize=7,
        show=False,
        ax=axis,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight", format="pdf")
    plt.close(figure)


def save_dotplot(adata, markers: Dict[str, list[str]], groupby: str, output_path: Path) -> None:
    """Purpose: render a marker dot plot for an annotated cell population.

    Parameters: ``adata`` is AnnData; ``markers`` maps groups to genes; ``groupby``
    identifies the observation field; ``output_path`` is the output PDF.
    Returns: ``None`` after writing the PDF.
    """
    # Summarize cell-level marker expression by the supplied annotation groups.
    plot = sc.pl.dotplot(
        adata,
        var_names=markers,
        groupby=groupby,
        standard_scale="var",
        color_map="Reds",
        dot_max=0.8,
        dot_min=0.1,
        show=False,
    )
    plot.savefig(output_path)
    plt.close("all")


def merge_lineage_labels(adata, lineage, lineage_name: str, major_mapping: Dict[str, str]) -> None:
    """Purpose: transfer lineage subtype and major labels to the full reference.

    Parameters: ``adata`` is the full AnnData object; ``lineage`` is its subset;
    ``lineage_name`` selects specified labels; ``major_mapping`` maps subtypes.
    Returns: ``None`` after updating ``adata.obs`` in place.
    """
    # Assign subtype labels within the lineage subset before deriving its broad annotation.
    apply_cluster_labels(
        lineage,
        CURRENT_SUBTYPE_CLUSTER_LABELS[lineage_name],
        output_key="Sub_type",
    )
    lineage.obs["Major_type"] = lineage.obs["Sub_type"].map(major_mapping).fillna(
        "Unknown"
    )
    # Transfer labels by cell identity so the full reference retains its original observation order.
    shared_cells = adata.obs_names[adata.obs_names.isin(lineage.obs_names)]
    adata.obs.loc[shared_cells, "Sub_type"] = lineage.obs.loc[shared_cells, "Sub_type"]
    adata.obs.loc[shared_cells, "Major_type"] = lineage.obs.loc[
        shared_cells, "Major_type"
    ]


def prepare_output_paths(output_dir: Path) -> Tuple[Path, Path]:
    """Purpose: create panel directories for Figure 7 scRNA outputs.

    Parameters: ``output_dir`` is the configured analysis-output directory.
    Returns: panel and Supplementary Figure 15 directory paths.
    """
    # Keep primary and supplementary panels under the configured analysis output root.
    panel_dir = output_dir / "panels"
    supplementary_dir = panel_dir / "supplementary" / "FigureS15"
    panel_dir.mkdir(parents=True, exist_ok=True)
    supplementary_dir.mkdir(parents=True, exist_ok=True)
    return panel_dir, supplementary_dir


def main(config: dict = CONFIG) -> None:
    """Purpose: execute the configured Figure 7 FDZS-3 integration workflow.

    Parameters: ``config`` supplies paths, the MatrixMarket prefix, and the epithelial-sample filter.
    Returns: ``None`` after writing the configured AnnData, R-bundle, TSV, and PDF outputs.
    """
    # Resolve runtime dependencies and output locations before loading cell-level matrices.
    load_analysis_dependencies()
    input_dir = Path(config["input_dir"])
    gtf_path = Path(config["gtf"])
    output_dir = Path(config["output_dir"])
    prefix = config["prefix"]
    excluded_epithelial_samples = tuple(config["excluded_epithelial_samples"])
    output_dir.mkdir(parents=True, exist_ok=True)
    panel_dir, supplementary_dir = prepare_output_paths(output_dir)

    # Load the processed MatrixMarket bundle as one cell-by-gene object; obs rows define cell identity.
    adata = load_processed_matrixmarket(input_dir, prefix=prefix)
    require_obs_columns(adata, ("patient", "doublet", "sample_id"))
    input_cells = adata.n_obs
    adata.obs["batch"] = adata.obs["patient"].astype(str)

    # Integrate patient batches on normalized cell profiles while preserving the shared gene axis.
    normalize_log1p(adata, target_sum=1e4)
    prepare_harmony_input(
        adata,
        batch_key="batch",
        n_top_genes=2000,
        scale_max_value=10,
        n_comps=50,
    )
    integrate_harmony(adata, batch_key="batch")
    adata.write(output_dir / "integrated_scRNA_stereoseq.h5ad")

    # Restrict annotation to singlet cells, whose fixed count anchors all downstream cell-level labels.
    adata = remove_doublets(adata, doublet_key="doublet", doublet_label="Doublet")
    if adata.n_obs != EXPECTED_ANNOTATED_CELLS:
        raise ValueError(
            "The FDZS-3 singlet reference contains "
            f"{EXPECTED_ANNOTATED_CELLS} cells; observed {adata.n_obs}"
        )
    cluster_major_cells(
        adata,
        resolution=0.8,
        n_neighbors=10,
        n_pcs=40,
        use_rep="X_pca_harmony",
    )
    rank_cluster_markers(adata)
    apply_cluster_labels(adata, CURRENT_MAJOR_CLUSTER_LABELS, output_key="Major_type")
    adata.uns["Major_type_colors"] = MAJOR_COLORS
    save_umap(adata, "Major_type", supplementary_dir / "FigureS15_major_umap.pdf", palette=MAJOR_COLORS)
    save_dotplot(adata, MAJOR_MARKERS, "Major_type", supplementary_dir / "FigureS15_major_marker_dotplot.pdf")
    adata.write(output_dir / "major_anno_all.h5ad")

    # Resolve immune subtypes within lineage-specific cell subsets, then map them to the full reference.
    adata.obs["Sub_type"] = "Unknown"
    b_plasma = preprocess_lineage(
        adata,
        ("B", "Plasma"),
        resolution=0.8,
        n_neighbors=10,
        n_pcs=50,
        use_rep="X_pca_harmony",
    )
    merge_lineage_labels(adata, b_plasma, "B_Plasma", B_PLASMA_TO_MAJOR)
    save_umap(b_plasma, "Sub_type", supplementary_dir / "FigureS15_b_plasma_umap.pdf")
    save_dotplot(b_plasma, B_PLASMA_MARKERS, "Sub_type", supplementary_dir / "FigureS15_b_plasma_marker_dotplot.pdf")

    nkt = preprocess_lineage(
        adata,
        ("T", "NK"),
        resolution=1.0,
        n_neighbors=10,
        n_pcs=50,
        use_rep="X_pca_harmony",
    )
    merge_lineage_labels(adata, nkt, "NKT", NKT_TO_MAJOR)
    save_umap(nkt, "Sub_type", supplementary_dir / "FigureS15_nkt_umap.pdf")
    save_dotplot(nkt, NKT_MARKERS, "Sub_type", supplementary_dir / "FigureS15_nkt_marker_dotplot.pdf")

    myeloid = preprocess_lineage(
        adata,
        ("Myeloid",),
        resolution=1.0,
        n_neighbors=10,
        n_pcs=50,
        use_rep="X_pca_harmony",
    )
    merge_lineage_labels(
        adata,
        myeloid,
        "Myeloid",
        {label: "Myeloid" for label in CURRENT_SUBTYPE_CLUSTER_LABELS["Myeloid"].values()},
    )
    save_umap(myeloid, "Sub_type", supplementary_dir / "FigureS15_myeloid_umap.pdf")
    save_dotplot(myeloid, MYELOID_MARKERS, "Sub_type", supplementary_dir / "FigureS15_myeloid_marker_dotplot.pdf")

    # Non-immune lineages retain their broad annotation as subtype unless epithelial CNV refines it.
    non_lineage_mask = ~adata.obs["Major_type"].isin(("T", "B", "NK", "Myeloid", "Plasma"))
    adata.obs.loc[non_lineage_mask, "Sub_type"] = adata.obs.loc[
        non_lineage_mask, "Major_type"
    ]

    # Align genes to genomic positions and infer epithelial CNV states against non-epithelial references.
    gtf = load_gene_gtf(str(gtf_path))
    add_gtf_annotations(adata, gtf)
    run_infercnv(
        adata,
        ("T", "Myeloid", "NK", "B", "Plasma", "Hepatocyte"),
        reference_key="Major_type",
        window_size=250,
    )
    epithelial = prepare_epithelial_subset(
        adata,
        excluded_samples=excluded_epithelial_samples,
        sample_key="sample_id",
        batch_key="batch",
        major_key="Major_type",
        epithelial_label="Epithelial",
        expression_resolution=0.8,
    )
    apply_cnv_status(
        epithelial,
        ("5", "25", "17", "26", "20", "28", "29"),
        cluster_key="cnv_leiden",
        output_key="Sub_type",
    )
    epithelial.obs["Major_type"] = "Epithelial"
    epithelial_cells = adata.obs_names[adata.obs_names.isin(epithelial.obs_names)]
    adata.obs.loc[epithelial_cells, "Sub_type"] = epithelial.obs.loc[
        epithelial_cells, "Sub_type"
    ]

    # Export final cell annotations and panels from the same in-memory reference state.
    save_umap(adata, "Major_type", panel_dir / "Figure7B_annotation_umap.pdf", palette=MAJOR_COLORS)
    save_umap(adata, "Sub_type", supplementary_dir / "FigureS15_subtype_umap.pdf")
    adata.write(output_dir / "final_annotated_scRNA.h5ad")

    # Build the R reference from identified cells, excluding only the configured epithelial sample cells.
    export_mask = ~(
        adata.obs["sample_id"].isin(excluded_epithelial_samples)
        & (adata.obs["Major_type"] == "Epithelial")
    ) & (adata.obs["Major_type"] != "Unknown")
    export_anndata_for_r(adata[export_mask].copy(), output_dir / "r_bundle")

    # Export the cell-population dimensions and analysis settings represented by these outputs.
    summary = pd.DataFrame(
        {
            "metric": [
                "input_cells",
                "annotated_cells",
                "exported_reference_cells",
                "major_clusters",
                "batch_key",
                "major_resolution",
                "b_plasma_resolution",
                "nkt_myeloid_resolution",
                "cnv_window_size",
                "excluded_epithelial_samples",
            ],
            "value": [
                input_cells,
                adata.n_obs,
                int(export_mask.sum()),
                adata.obs["leiden"].nunique(),
                "batch=patient",
                0.8,
                0.8,
                1.0,
                250,
                ";".join(excluded_epithelial_samples),
            ],
        }
    )
    summary.to_csv(output_dir / "scrna_run_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()

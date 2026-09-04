#!/usr/bin/env python3
"""Plot the public Figure 5B model-development count component.

Target/purpose: render the nested donut plot for validated four-class tile
counts and the grouped development/internal-holdout split. Inputs: the
``Figure5B_training_count_provenance.tsv`` table produced by
``01_prepare_panel_b_counts.py``. Outputs: one PDF and its plotted-data TSV in
the configured directory. The held discovery/test inference-tile component is
not represented because its denominators are not public reproducible inputs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


CLASS_ORDER = ["PIR", "PSM", "OTHER_PT", "TC"]
CLASS_LABELS = {
    "PIR": "PIR-Niche",
    "PSM": "PSM-Niche",
    "OTHER_PT": "Other PT",
    "TC": "TC",
}
CLASS_COLORS = {
    "PIR": "#9EC3D5",
    "PSM": "#2176A8",
    "OTHER_PT": "#AED581",
    "TC": "#399C39",
}
SPLIT_ORDER = ["development", "internal_holdout"]
SPLIT_LABELS = {
    "development": "Model-development subset (80%)",
    "internal_holdout": "Common internal holdout subset (20%)",
}
SPLIT_COLORS = {"development": "#F0C7C9", "internal_holdout": "#9DB6DE"}

# Edit these paths before running this interactive Figure 5 entry point.
CONFIG = {
    "counts": Path("/path/to/Figure5B_training_count_provenance.tsv"),
    "output_dir": Path("/path/to/Figure5B"),
}


def _require_columns(frame: pd.DataFrame, path: Path) -> None:
    """Require the published Figure 5B count-table schema.

    Parameters: ``frame`` is the source table and ``path`` identifies it in an
    error message. Returns: no value; raises ``ValueError`` on schema mismatch.
    """
    required = {
        "count_group",
        "category",
        "count",
        "analysis_unit",
        "source_path",
        "selection_rule",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")


def _count_rows(
    frame: pd.DataFrame, count_group: str, categories: list[str]
) -> pd.DataFrame:
    """Select one ordered, complete count-group subset.

    Parameters: ``frame`` is the input table, ``count_group`` identifies a
    table role, and ``categories`` defines its required order. Returns: the
    selected rows in that order.
    """
    subset = frame.loc[frame["count_group"].eq(count_group)].copy()
    if subset["category"].duplicated().any() or set(subset["category"]) != set(categories):
        raise ValueError(
            f"{count_group} must contain exactly these categories: {categories}"
        )
    return subset.set_index("category").loc[categories].reset_index()


def read_counts(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read and reconcile the public Figure 5B training-count table.

    Parameters: ``path`` identifies the TSV produced by the panel-B count
    preparation entry point. Returns: ordered class and split count rows.
    """
    frame = pd.read_csv(path, sep="\t")
    _require_columns(frame, path)
    frame = frame.copy()
    frame["count"] = pd.to_numeric(frame["count"], errors="raise")
    if not np.isfinite(frame["count"].to_numpy(dtype=float)).all():
        raise ValueError("Figure 5B counts must be finite")
    if not np.equal(frame["count"], np.floor(frame["count"])).all() or (frame["count"] < 0).any():
        raise ValueError("Figure 5B counts must be non-negative integers")
    frame["count"] = frame["count"].astype(int)

    class_rows = _count_rows(frame, "training_class_tiles", CLASS_ORDER)
    split_rows = _count_rows(frame, "split_tiles", SPLIT_ORDER)
    total_rows = frame.loc[
        frame["count_group"].eq("training_total_tiles")
        & frame["category"].eq("all_classes")
    ]
    if len(total_rows) != 1:
        raise ValueError("Figure 5B count table must contain one training_total_tiles row")
    total = int(total_rows.iloc[0]["count"])
    if int(class_rows["count"].sum()) != total:
        raise ValueError("Training class counts do not sum to the training total")
    if int(split_rows["count"].sum()) != total:
        raise ValueError("Split tile counts do not sum to the training total")
    return class_rows, split_rows


def build_plotted_data(class_rows: pd.DataFrame, split_rows: pd.DataFrame) -> pd.DataFrame:
    """Create the exact normalized values displayed in the nested donut plot.

    Parameters: ``class_rows`` and ``split_rows`` are validated ordered table
    subsets. Returns: a plotting table with count fractions and display labels.
    """
    total = int(class_rows["count"].sum())
    outer = class_rows.assign(
        ring="class_tiles",
        display_label=class_rows["category"].map(CLASS_LABELS),
        display_color=class_rows["category"].map(CLASS_COLORS),
        fraction=class_rows["count"] / total,
    )
    inner = split_rows.assign(
        ring="grouped_split_tiles",
        display_label=split_rows["category"].map(SPLIT_LABELS),
        display_color=split_rows["category"].map(SPLIT_COLORS),
        fraction=split_rows["count"] / total,
    )
    columns = [
        "ring",
        "count_group",
        "category",
        "display_label",
        "count",
        "fraction",
        "analysis_unit",
        "source_path",
        "selection_rule",
        "display_color",
    ]
    return pd.concat([outer[columns], inner[columns]], ignore_index=True)


def plot_panel(plotted_data: pd.DataFrame, output: Path) -> None:
    """Render the supported Figure 5B nested training-count donut plot.

    Parameters: ``plotted_data`` contains validated class and split plot values;
    ``output`` identifies the PDF to write. Returns: no value.
    """
    outer = plotted_data.loc[plotted_data["ring"].eq("class_tiles")]
    inner = plotted_data.loc[plotted_data["ring"].eq("grouped_split_tiles")]
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure, axis = plt.subplots(figsize=(6.2, 4.0))
    axis.pie(
        outer["count"],
        colors=outer["display_color"],
        startangle=90,
        counterclock=False,
        radius=1.0,
        wedgeprops={"width": 0.31, "edgecolor": "white", "linewidth": 1.2},
    )
    axis.pie(
        inner["count"],
        colors=inner["display_color"],
        startangle=90,
        counterclock=False,
        radius=0.67,
        wedgeprops={"width": 0.31, "edgecolor": "white", "linewidth": 1.2},
    )
    total = int(outer["count"].sum())
    axis.text(0, 0.05, f"{total:,}", ha="center", va="center", fontsize=15, fontweight="bold")
    axis.text(0, -0.15, "tile records", ha="center", va="center", fontsize=8, color="#444444")
    axis.set_title("Model-development dataset composition", fontsize=10, fontweight="bold", pad=14)
    axis.text(-1.23, 1.12, "B", fontsize=12, fontweight="bold")
    class_handles = [
        Patch(facecolor=color, edgecolor="white", label=f"{CLASS_LABELS[name]} ({int(outer.loc[outer['category'].eq(name), 'count'].iloc[0]):,})")
        for name, color in CLASS_COLORS.items()
    ]
    split_handles = [
        Patch(facecolor=color, edgecolor="white", label=SPLIT_LABELS[name])
        for name, color in SPLIT_COLORS.items()
    ]
    class_legend = axis.legend(
        handles=class_handles,
        title="Tile class",
        loc="upper left",
        bbox_to_anchor=(1.00, 0.98),
        frameon=False,
        fontsize=8,
        title_fontsize=8.5,
        handlelength=1.4,
        handletextpad=0.55,
    )
    axis.add_artist(class_legend)
    axis.legend(
        handles=split_handles,
        title="Grouped split",
        loc="lower left",
        bbox_to_anchor=(1.00, 0.02),
        frameon=False,
        fontsize=8,
        title_fontsize=8.5,
        handlelength=1.4,
        handletextpad=0.55,
    )
    axis.set_aspect("equal")
    figure.subplots_adjust(left=0.08, right=0.63, bottom=0.08, top=0.88)
    figure.savefig(
        output,
        format="pdf",
        bbox_inches="tight",
        metadata={"Title": "Figure 5B model-development dataset composition"},
    )
    plt.close(figure)


def run(config: dict[str, Path]) -> None:
    """Write Figure 5B plotted data and the supported nested-donut PDF.

    Parameters: ``config`` supplies the public count-table path and output
    directory. Returns: no value; writes the TSV and PDF.
    """
    class_rows, split_rows = read_counts(config["counts"])
    plotted_data = build_plotted_data(class_rows, split_rows)
    output_dir = config["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    plotted_data.to_csv(
        output_dir / "Figure5B_model_development_plotted_data.tsv",
        sep="\t",
        index=False,
    )
    plot_panel(plotted_data, output_dir / "Figure5B_model_development_counts.pdf")


if __name__ == "__main__":
    run(CONFIG)

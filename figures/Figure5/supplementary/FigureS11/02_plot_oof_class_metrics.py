#!/usr/bin/env python3
"""Plot Supplementary Figure 11B-C SpMap performance summaries.

Target/purpose: render the pooled OOF and five-seed summaries.
Inputs: the four small summary tables selected in ``CONFIG``. Outputs: two
PDFs and their plotted-data/statistics TSVs in the configured directory.
Ordered workflow: validate tables, calculate intervals, write normalized
tables, render the confusion panel, then render the class-performance panel.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy.stats import t


CLASS_ORDER = ["TC", "PIR", "PSM", "OTHER_PT"]
CLASS_LABELS = {"TC": "TC", "PIR": "PIR", "PSM": "PSM", "OTHER_PT": "Other PT"}
CLASS_COLORS = {"TC": "#555555", "PIR": "#3F8FB5", "PSM": "#E79A55", "OTHER_PT": "#A5A5A5"}
SEEDS = [24, 101, 202, 303, 404]
METRICS = [("precision", "Precision"), ("recall", "Recall"), ("auroc_ovr", "AUROC"), ("auprc_ovr", "AUPRC")]

# Edit these paths before running this interactive supplementary-figure entry point.
CONFIG = {
    "confusion_counts": Path("/path/to/confusion_counts.csv"),
    "confusion_row_proportions": Path("/path/to/confusion_row_proportions.csv"),
    "seed_class_metrics": Path("/path/to/seed_class_metrics.csv"),
    "seed_overall_metrics": Path("/path/to/seed_overall_metrics.csv"),
    "output_dir": Path("/path/to/FigureS11"),
}


def _read_matrix(path: Path, name: str) -> pd.DataFrame:
    """Read one class-indexed confusion matrix.

    Parameters: ``path`` is the CSV matrix and ``name`` identifies it in errors.
    Returns: a finite floating-point matrix in the fixed four-class order.
    """
    # Preserve the fixed true-class and predicted-class axes across both matrices.
    matrix = pd.read_csv(path, index_col=0)
    if list(matrix.index) != CLASS_ORDER or list(matrix.columns) != CLASS_ORDER:
        raise ValueError(f"{name} must use the exact class order {CLASS_ORDER}")
    matrix = matrix.astype(float)
    if not np.isfinite(matrix.to_numpy()).all():
        raise ValueError(f"{name} contains non-finite values")
    return matrix


def read_confusion(count_path: Path, proportion_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read and validate the fixed 49,001-tile S11B confusion contract.

    Parameters: ``count_path`` and ``proportion_path`` identify the two CSVs.
    Returns: integer confusion counts and matching row-normalized proportions.
    """
    # Reconcile integer counts with their row-normalized plotting representation.
    counts = _read_matrix(count_path, "Confusion counts")
    if not np.equal(counts.to_numpy(), np.floor(counts.to_numpy())).all() or (counts.to_numpy() < 0).any():
        raise ValueError("Confusion counts must be non-negative integers")
    counts = counts.astype(int)
    if int(counts.to_numpy().sum()) != 49001:
        raise ValueError("Confusion counts must sum to 49,001 pooled OOF tiles")
    proportions = _read_matrix(proportion_path, "Confusion row proportions")
    if not np.allclose(proportions.sum(axis=1), 1.0, rtol=0, atol=1e-8):
        raise ValueError("Confusion row proportions must sum to one")
    expected = counts.div(counts.sum(axis=1), axis=0)
    if not np.allclose(proportions.to_numpy(), expected.to_numpy(), rtol=0, atol=1e-8):
        raise ValueError("Confusion row proportions do not match the count matrix")
    return counts, proportions


def read_seed_metrics(class_path: Path, overall_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read the complete five-seed S11C metric tables.

    Parameters: ``class_path`` identifies class-level metrics and ``overall_path``
    identifies macro-F1 metrics.
    Returns: ordered class-level and overall five-seed dataframes.
    """
    # Map class-level metrics to the complete seed-by-class grid in display order.
    class_metrics = pd.read_csv(class_path)
    required_class = {"seed", "class_name", *(metric for metric, _ in METRICS)}
    missing_class = required_class - set(class_metrics.columns)
    if missing_class:
        raise ValueError(f"Seed class metrics are missing columns: {sorted(missing_class)}")
    class_metrics = class_metrics.copy()
    class_metrics["seed"] = class_metrics["seed"].astype(int)
    observed = set(zip(class_metrics["seed"], class_metrics["class_name"]))
    expected = {(seed, class_name) for seed in SEEDS for class_name in CLASS_ORDER}
    if observed != expected or len(class_metrics) != len(expected):
        raise ValueError("Seed class metrics must contain exactly one row per configured seed and class")
    class_metrics = class_metrics.set_index(["seed", "class_name"]).loc[
        pd.MultiIndex.from_product([SEEDS, CLASS_ORDER], names=["seed", "class_name"])
    ].reset_index()
    metric_columns = [metric for metric, _ in METRICS]
    if not np.isfinite(class_metrics[metric_columns].to_numpy()).all() or not ((class_metrics[metric_columns] >= 0) & (class_metrics[metric_columns] <= 1)).all().all():
        raise ValueError("Seed class metrics must be finite proportions")

    # Keep overall macro-F1 on the same ordered five-seed axis.
    overall_metrics = pd.read_csv(overall_path)
    required_overall = {"seed", "f1_macro"}
    missing_overall = required_overall - set(overall_metrics.columns)
    if missing_overall:
        raise ValueError(f"Seed overall metrics are missing columns: {sorted(missing_overall)}")
    overall_metrics = overall_metrics[["seed", "f1_macro"]].copy()
    overall_metrics["seed"] = overall_metrics["seed"].astype(int)
    overall_metrics = overall_metrics.sort_values("seed").reset_index(drop=True)
    if overall_metrics["seed"].tolist() != SEEDS or overall_metrics["seed"].duplicated().any():
        raise ValueError("Seed overall metrics must contain seeds 24, 101, 202, 303, and 404 once")
    if not np.isfinite(overall_metrics["f1_macro"]).all() or not overall_metrics["f1_macro"].between(0, 1).all():
        raise ValueError("Macro-F1 values must be finite proportions")
    return class_metrics, overall_metrics


def interval_table(class_metrics: pd.DataFrame, overall_metrics: pd.DataFrame) -> pd.DataFrame:
    """Calculate two-sided 95% t intervals across the five specified seeds.

    Parameters: ``class_metrics`` and ``overall_metrics`` are validated seed tables.
    Returns: one interval row per class metric plus overall macro-F1.
    """
    # Apply one interval definition to every class metric and overall macro-F1.
    rows: list[dict[str, object]] = []
    critical = t.ppf(0.975, df=len(SEEDS) - 1)
    for metric, label in METRICS:
        for class_name in CLASS_ORDER:
            values = class_metrics.loc[class_metrics["class_name"] == class_name, metric].to_numpy(dtype=float)
            mean = float(values.mean())
            standard_deviation = float(values.std(ddof=1))
            margin = float(critical * standard_deviation / np.sqrt(len(values)))
            rows.append({"metric": label, "class_name": class_name, "seeds": len(values), "mean": mean, "standard_deviation": standard_deviation, "ci95_lower": mean - margin, "ci95_upper": mean + margin, "interval_definition": "two-sided t interval across five training seeds"})
    values = overall_metrics["f1_macro"].to_numpy(dtype=float)
    mean = float(values.mean())
    standard_deviation = float(values.std(ddof=1))
    margin = float(critical * standard_deviation / np.sqrt(len(values)))
    rows.append({"metric": "Macro-F1", "class_name": "overall", "seeds": len(values), "mean": mean, "standard_deviation": standard_deviation, "ci95_lower": mean - margin, "ci95_upper": mean + margin, "interval_definition": "two-sided t interval across five training seeds"})
    return pd.DataFrame(rows)


def plot_confusion(counts: pd.DataFrame, proportions: pd.DataFrame, output: Path) -> None:
    """Render the S11B four-class pooled OOF confusion PDF.

    Parameters: ``counts`` and ``proportions`` are validated confusion matrices;
    ``output`` is the target PDF path.
    Returns: no value; writes the PDF.
    """
    # Encode row percentages by color while retaining counts in each matrix cell.
    cmap = LinearSegmentedColormap.from_list("spmap_confusion", ["#F7F9FA", "#B8D5E5", "#3F8FB5", "#174A66"])
    figure, axis = plt.subplots(figsize=(4.45, 3.95))
    image = axis.imshow(proportions.to_numpy() * 100, cmap=cmap, vmin=0, vmax=100)
    labels = [CLASS_LABELS[name] for name in CLASS_ORDER]
    axis.set_xticks(range(4), labels=labels)
    axis.set_yticks(range(4), labels=labels)
    axis.set_xlabel("Predicted class", labelpad=7)
    axis.set_ylabel("True class", labelpad=7)
    axis.set_title("Four-class pooled out-of-fold confusion matrix", fontweight="bold", pad=10)
    axis.tick_params(axis="both", length=0, pad=5)
    axis.set_xticks(np.arange(-0.5, 4, 1), minor=True)
    axis.set_yticks(np.arange(-0.5, 4, 1), minor=True)
    axis.grid(which="minor", color="white", linewidth=1.2)
    axis.tick_params(which="minor", bottom=False, left=False)
    for ticks in (axis.get_xticklabels(), axis.get_yticklabels()):
        for tick, class_name in zip(ticks, CLASS_ORDER):
            tick.set_color(CLASS_COLORS[class_name])
            if class_name in ("PIR", "PSM"):
                tick.set_fontweight("bold")
    for row in range(4):
        for column in range(4):
            percentage = proportions.iloc[row, column] * 100
            axis.text(column, row, f"{counts.iloc[row, column]:,}\n({percentage:.1f}%)", ha="center", va="center", color="white" if percentage >= 42 else "#222222", fontsize=7.5, linespacing=1.25)
    axis.add_patch(Rectangle((0.5, 0.5), 2, 2, fill=False, edgecolor="#202020", linewidth=1.1, linestyle=(0, (3, 2))))
    colorbar = figure.colorbar(image, ax=axis, fraction=0.046, pad=0.045)
    colorbar.set_label("True-class row percentage (%)", rotation=270, labelpad=13)
    colorbar.set_ticks([0, 25, 50, 75, 100])
    colorbar.outline.set_linewidth(0.5)
    axis.text(-0.20, 1.11, "B", transform=axis.transAxes, fontsize=11, fontweight="bold")
    figure.subplots_adjust(left=0.19, right=0.86, bottom=0.18, top=0.84)
    figure.savefig(output, format="pdf", bbox_inches="tight", metadata={"Title": "Supplementary Figure 11B four-class pooled OOF confusion matrix"})
    plt.close(figure)


def plot_performance(class_metrics: pd.DataFrame, overall_metrics: pd.DataFrame, intervals: pd.DataFrame, output: Path) -> None:
    """Render the S11C class-specific five-seed performance PDF.

    Parameters: ``class_metrics`` and ``overall_metrics`` supply observed values,
    ``intervals`` supplies t intervals, and ``output`` is the PDF path.
    Returns: no value; writes the PDF.
    """
    # Allocate stable metric and class offsets before adding seed values and intervals.
    metric_positions = {metric: position for position, (metric, _) in enumerate(METRICS)}
    class_offsets = dict(zip(CLASS_ORDER, [-0.24, -0.08, 0.08, 0.24]))
    seed_jitter = dict(zip(SEEDS, np.linspace(-0.025, 0.025, len(SEEDS))))
    figure, axis = plt.subplots(figsize=(7.2, 4.35))
    for metric, label in METRICS:
        for class_name in CLASS_ORDER:
            values = class_metrics.loc[class_metrics["class_name"] == class_name].set_index("seed").loc[SEEDS, metric].astype(float)
            x_position = metric_positions[metric] + class_offsets[class_name]
            for seed, value in values.items():
                axis.scatter(x_position + seed_jitter[seed], value, marker="x", s=19, color=CLASS_COLORS[class_name], linewidth=0.85, alpha=0.65, zorder=2)
            summary = intervals.loc[(intervals["metric"] == label) & (intervals["class_name"] == class_name)].iloc[0]
            axis.errorbar(x_position, summary["mean"], yerr=[[summary["mean"] - summary["ci95_lower"]], [summary["ci95_upper"] - summary["mean"]]], fmt="o", markersize=4.8, markerfacecolor=CLASS_COLORS[class_name], markeredgecolor="#202020", markeredgewidth=0.45, ecolor=CLASS_COLORS[class_name], elinewidth=1.35, capsize=2.8, capthick=1.0, zorder=3)
    # Plot overall macro-F1 separately from the four class-specific metric groups.
    macro_values = overall_metrics.set_index("seed").loc[SEEDS, "f1_macro"].astype(float)
    for seed, value in macro_values.items():
        axis.scatter(4 + seed_jitter[seed], value, marker="x", s=19, color="#222222", linewidth=0.85, alpha=0.65, zorder=2)
    macro = intervals.loc[intervals["metric"] == "Macro-F1"].iloc[0]
    axis.errorbar(4, macro["mean"], yerr=[[macro["mean"] - macro["ci95_lower"]], [macro["ci95_upper"] - macro["mean"]]], fmt="D", markersize=5.0, markerfacecolor="#222222", markeredgecolor="#222222", ecolor="#222222", elinewidth=1.35, capsize=2.8, capthick=1.0, zorder=3)
    axis.hlines(0.5, 2 - 0.44, 2 + 0.44, color="#777777", linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
    axis.text(2, 0.515, "Random AUROC = 0.5", color="#666666", fontsize=6.8, ha="center", va="bottom")
    axis.set_xlim(-0.55, 4.45)
    axis.set_ylim(0, 1.025)
    axis.set_xticks(range(5), [label for _, label in METRICS] + ["Macro-F1"])
    axis.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axis.set_yticklabels(["0", ".25", ".50", ".75", "1"])
    axis.set_ylabel("Performance", labelpad=7)
    axis.grid(axis="y", color="#E2E2E2", linewidth=0.55)
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(axis="x", length=0, pad=7)
    axis.tick_params(axis="y", length=3, color="#444444")
    handles = [Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=CLASS_COLORS[name], markeredgecolor="#222222", markeredgewidth=0.4, markersize=5, label=CLASS_LABELS[name]) for name in CLASS_ORDER]
    handles.extend([Line2D([0], [0], marker="D", linestyle="none", markerfacecolor="#222222", markeredgecolor="#222222", markersize=5, label="Overall macro-F1"), Line2D([0], [0], linestyle=(0, (4, 3)), color="#777777", linewidth=1.0, label="Random AUROC")])
    axis.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.20), frameon=False, ncol=6, handlelength=1.4, handletextpad=0.45, columnspacing=1.0, fontsize=7.2)
    axis.set_title("Class-specific performance across five training seeds", fontsize=9, fontweight="bold", pad=45)
    axis.text(-0.09, 1.18, "C", transform=axis.transAxes, fontsize=11, fontweight="bold")
    figure.text(0.985, 0.025, "Crosses show individual seeds; circles/diamond show means; error bars are 95% t intervals across seeds.", fontsize=6.8, color="#555555", ha="right")
    figure.subplots_adjust(left=0.10, right=0.985, bottom=0.17, top=0.78)
    figure.savefig(output, format="pdf", bbox_inches="tight", metadata={"Title": "Supplementary Figure 11C class-specific performance across training seeds"})
    plt.close(figure)


def run(config: dict[str, Path]) -> None:
    """Write and plot S11B-C outputs from editable in-script configuration.

    Parameters: ``config`` supplies the four input CSV paths and ``output_dir``.
    Returns: no value; writes the configured TSV and PDF artifacts.
    """
    # Orchestrate metric preparation, normalized table export, and both figure panels.
    counts, proportions = read_confusion(
        config["confusion_counts"], config["confusion_row_proportions"]
    )
    class_metrics, overall_metrics = read_seed_metrics(
        config["seed_class_metrics"], config["seed_overall_metrics"]
    )
    intervals = interval_table(class_metrics, overall_metrics)
    output_dir = config["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    counts.rename_axis("true_class").reset_index().to_csv(output_dir / "FigureS11B_confusion_counts.tsv", sep="\t", index=False)
    proportions.rename_axis("true_class").reset_index().to_csv(output_dir / "FigureS11B_confusion_row_proportions.tsv", sep="\t", index=False)
    class_metrics.to_csv(output_dir / "FigureS11C_seed_class_metrics.tsv", sep="\t", index=False)
    overall_metrics.to_csv(output_dir / "FigureS11C_seed_macro_f1.tsv", sep="\t", index=False)
    intervals.to_csv(output_dir / "FigureS11C_metric_t_intervals.tsv", sep="\t", index=False)
    plot_confusion(counts, proportions, output_dir / "FigureS11B_confusion_matrix.pdf")
    plot_performance(class_metrics, overall_metrics, intervals, output_dir / "FigureS11C_class_metrics.pdf")


if __name__ == "__main__":
    run(CONFIG)

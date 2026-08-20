#!/usr/bin/env python3
"""Plot Supplementary Figure 11D from specified patient-level ISR score tables.

Target/purpose: correlate matched IMC-derived and H&E-derived ISR values.
Inputs: unique-PID IMC and unique-WSI H&E score CSVs set in ``CONFIG``.
Outputs: matched-patient/statistics TSVs and one PDF in ``output_dir``. Ordered
workflow: validate ISR arithmetic, perform the one-to-one join, calculate Pearson
statistics, write tables, and render the correlation panel.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


EXPECTED_MATCHED_PATIENTS = 33

# Edit these paths before running this interactive supplementary-figure entry point.
CONFIG = {
    "imc_scores": Path("/path/to/imc_scores.csv"),
    "he_scores": Path("/path/to/he_scores.csv"),
    "output_dir": Path("/path/to/FigureS11"),
}


def _require_columns(frame: pd.DataFrame, columns: set[str], source: str) -> None:
    """Fail when an ISR input table lacks a required column.

    Parameters: ``frame`` is the loaded table, ``columns`` is the required set,
    and ``source`` identifies the table in errors.
    Returns: no value; raises ``ValueError`` when a column is absent.
    """
    # Establish the identifier, count, and saved-score fields used for matching.
    missing = columns - set(frame.columns)
    if missing:
        raise ValueError(f"{source} is missing columns: {sorted(missing)}")


def fisher_z_interval(r_value: float, n: int) -> tuple[float, float]:
    """Calculate a two-sided 95% Fisher-z interval for Pearson correlation.

    Parameters: ``r_value`` is the observed Pearson correlation and ``n`` is the
    matched-patient count.
    Returns: lower and upper correlation-scale confidence limits.
    """
    # The Fisher transform requires a non-perfect correlation and residual degrees of freedom.
    if n <= 3 or not -1 < r_value < 1:
        raise ValueError("Fisher-z confidence interval requires n > 3 and a non-perfect correlation")
    z_value = np.arctanh(r_value)
    margin = stats.norm.ppf(0.975) / np.sqrt(n - 3)
    return tuple(np.tanh([z_value - margin, z_value + margin]))


def load_and_match(imc_path: Path, he_path: Path) -> pd.DataFrame:
    """Validate source ISR arithmetic and create the specified one-to-one join.

    Parameters: ``imc_path`` and ``he_path`` identify the specified source CSVs.
    Returns: a sorted 33-patient dataframe containing both derived ISR values.
    """
    # Map each source to one identifier and its operational ISR definition.
    imc = pd.read_csv(imc_path)
    he = pd.read_csv(he_path)
    _require_columns(imc, {"PID", "n_bdme_i", "n_bdme_s", "easy_bdme_score"}, "IMC scores")
    _require_columns(he, {"wsi_id", "bile_1_count", "bile_2_count", "ratio_bile1_to_bile2"}, "H&E scores")
    if imc["PID"].duplicated().any() or he["wsi_id"].duplicated().any():
        raise ValueError("IMC PID and H&E wsi_id values must each be unique")
    imc = imc[["PID", "n_bdme_i", "n_bdme_s", "easy_bdme_score"]].copy()
    he = he[["wsi_id", "bile_1_count", "bile_2_count", "ratio_bile1_to_bile2"]].copy()
    numeric_imc = ["n_bdme_i", "n_bdme_s", "easy_bdme_score"]
    numeric_he = ["bile_1_count", "bile_2_count", "ratio_bile1_to_bile2"]
    if not np.isfinite(imc[numeric_imc].to_numpy(dtype=float)).all() or not np.isfinite(he[numeric_he].to_numpy(dtype=float)).all():
        raise ValueError("ISR score inputs must be finite")
    expected_imc = np.clip((imc["n_bdme_i"] + 1) / (imc["n_bdme_s"] + 1), 0.1, 10)
    if not np.allclose(imc["easy_bdme_score"], expected_imc, rtol=0, atol=1e-8):
        raise ValueError("Saved IMC ISR does not match the operational clipped pseudocount definition")
    if (he["bile_2_count"] == 0).any():
        raise ZeroDivisionError("H&E PSM tile count is zero")
    expected_he = he["bile_1_count"] / he["bile_2_count"]
    if not np.allclose(he["ratio_bile1_to_bile2"], expected_he, rtol=0, atol=1e-8):
        raise ValueError("Saved H&E ISR does not match the literal PIR/PSM ratio")
    # Preserve patient-level inference through an explicitly one-to-one inner match.
    matched = imc.merge(he, left_on="PID", right_on="wsi_id", how="inner", validate="one_to_one")
    if len(matched) != EXPECTED_MATCHED_PATIENTS:
        raise ValueError(f"Expected {EXPECTED_MATCHED_PATIENTS} matched discovery-cohort patients, found {len(matched)}")
    matched = matched.rename(columns={"easy_bdme_score": "imc_derived_isr", "ratio_bile1_to_bile2": "he_derived_isr"}).drop(columns="wsi_id")
    return matched.sort_values("PID").reset_index(drop=True)


def correlation_statistics(matched: pd.DataFrame) -> dict[str, object]:
    """Calculate the reported Pearson correlation and Fisher-z interval.

    Parameters: ``matched`` is the validated one-row-per-patient ISR dataframe.
    Returns: a statistics dictionary written to the S11D results TSV.
    """
    # Prepare paired patient vectors once for the correlation and confidence interval.
    x_values = matched["imc_derived_isr"].to_numpy(dtype=float)
    y_values = matched["he_derived_isr"].to_numpy(dtype=float)
    r_value, p_value = stats.pearsonr(x_values, y_values)
    ci_low, ci_high = fisher_z_interval(float(r_value), len(matched))
    return {"n": len(matched), "pearson_r": r_value, "pearson_p_two_sided": p_value, "pearson_ci95_low_fisher_z": ci_low, "pearson_ci95_high_fisher_z": ci_high, "imc_isr_definition": "clip((n_bdme_i + 1) / (n_bdme_s + 1), 0.1, 10)", "he_isr_definition": "bile_1_count / bile_2_count", "ci_definition": "two-sided 95% Fisher-z confidence interval"}


def plot_panel(matched: pd.DataFrame, summary: dict[str, object], output: Path) -> None:
    """Render the S11D patient-level ISR correlation PDF.

    Parameters: ``matched`` contains paired ISR values, ``summary`` contains
    correlation statistics, and ``output`` is the target PDF path.
    Returns: no value; writes the PDF.
    """
    # Fit the displayed ordinary least-squares mean and its pointwise confidence band.
    x_values = matched["imc_derived_isr"].to_numpy(dtype=float)
    y_values = matched["he_derived_isr"].to_numpy(dtype=float)
    n = len(matched)
    design = np.column_stack([np.ones(n), x_values])
    intercept, slope = np.linalg.lstsq(design, y_values, rcond=None)[0]
    x_grid = np.linspace(x_values.min(), x_values.max(), 300)
    fitted = intercept + slope * x_grid
    residuals = y_values - (intercept + slope * x_values)
    mse = np.sum(residuals ** 2) / (n - 2)
    sxx = np.sum((x_values - x_values.mean()) ** 2)
    if sxx == 0:
        raise ValueError("IMC ISR values must vary to plot the fitted mean confidence band")
    fit_se = np.sqrt(mse * (1 / n + (x_grid - x_values.mean()) ** 2 / sxx))
    critical = stats.t.ppf(0.975, df=n - 2)
    # Render paired patients, the fitted mean, and the reported correlation together.
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 8, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "pdf.fonttype": 42, "ps.fonttype": 42})
    figure, axis = plt.subplots(figsize=(3.5, 3.15))
    axis.fill_between(x_grid, fitted - critical * fit_se, fitted + critical * fit_se, color="#B3B3B3", alpha=0.38, linewidth=0)
    axis.plot(x_grid, fitted, color="#333333", linewidth=1.25, zorder=2)
    axis.scatter(x_values, y_values, s=25, color="#0072B2", edgecolor="white", linewidth=0.5, alpha=0.9, zorder=3)
    annotation = f"Pearson r = {summary['pearson_r']:.3f}\n95% CI: {summary['pearson_ci95_low_fisher_z']:.3f} to {summary['pearson_ci95_high_fisher_z']:.3f}\nP = {summary['pearson_p_two_sided']:.2e}\nn = {summary['n']}"
    axis.text(0.04, 0.96, annotation, transform=axis.transAxes, ha="left", va="top", fontsize=7.5, linespacing=1.35)
    axis.set_xlabel("IMC-derived ISR")
    axis.set_ylabel("H&E-derived ISR")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(direction="out", length=3, width=0.7)
    axis.margins(x=0.05, y=0.08)
    figure.tight_layout(pad=0.8)
    figure.savefig(output, format="pdf", metadata={"Title": "Supplementary Figure 11D patient-level IMC and H&E ISR correlation"})
    plt.close(figure)


def run(config: dict[str, Path]) -> None:
    """Write and plot S11D outputs from editable in-script configuration.

    Parameters: ``config`` supplies IMC/H&E score paths and ``output_dir``.
    Returns: no value; writes the matched table, statistics table, and PDF.
    """
    # Orchestrate matching, patient-level statistics, tabular export, and plotting.
    matched = load_and_match(config["imc_scores"], config["he_scores"])
    summary = correlation_statistics(matched)
    output_dir = config["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    matched.to_csv(output_dir / "FigureS11D_matched_patient_isr.tsv", sep="\t", index=False)
    pd.DataFrame([summary]).to_csv(output_dir / "FigureS11D_correlation_statistics.tsv", sep="\t", index=False)
    plot_panel(matched, summary, output_dir / "FigureS11D_isr_correlation.pdf")


if __name__ == "__main__":
    run(CONFIG)

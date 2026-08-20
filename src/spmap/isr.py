"""Aggregate four-class SpMap tile predictions into WSI-level H&E ISR.

Purpose:
    Count predicted tile classes per WSI and calculate the literal PIR/PSM ISR.
Figure 5 callers:
    Figure 5 SpMap workflows use ``aggregate_isr`` after inference to prepare
    H&E ISR inputs for Supplementary Figure S11D.
Inputs:
    A prediction DataFrame containing a WSI grouping column and four-class
    integer predictions.
Outputs:
    One row per WSI with class counts, total 256-pixel tiles, and ISR.
Ordered use:
    Run ``wsi.py`` and ``inference.py`` first, then aggregate their prediction
    table with ``aggregate_isr``.
"""

from __future__ import annotations

import pandas as pd


CLASS_COUNT_COLUMNS = ("tc_count", "pir_count", "psm_count", "other_pt_count")


def aggregate_isr(
    predictions: pd.DataFrame,
    *,
    group_col: str = "wsi_id",
    prediction_col: str = "pred_class",
) -> pd.DataFrame:
    """Count predicted classes and compute literal PIR/PSM ISR per group.

    Parameters:
        predictions: Tile prediction rows containing ``group_col`` and
            ``prediction_col``.
        group_col: Column defining one WSI-level aggregation group.
        prediction_col: Integer four-class prediction column.

    Returns:
        A DataFrame with the grouping column, all fixed class-count columns,
        ``total_256_tiles``, and ``isr``.
    """
    # Require the WSI grouping key and canonical four-class prediction field.
    required = {group_col, prediction_col}
    missing = required - set(predictions.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    # Preserve the full aggregation schema when no prediction rows are supplied.
    if predictions.empty:
        return pd.DataFrame(
            columns=[group_col, *CLASS_COUNT_COLUMNS, "total_256_tiles", "isr"]
        )
    # Normalize and validate class IDs before WSI-level counting.
    classes = predictions[prediction_col].astype(int)
    invalid = sorted(set(classes) - set(range(4)))
    if invalid:
        raise ValueError(f"Invalid predicted class IDs: {invalid}")

    # Count every canonical class per WSI and materialize absent classes as zero.
    counts = (
        predictions.assign(**{prediction_col: classes})
        .groupby([group_col, prediction_col], sort=False)
        .size()
        .unstack(fill_value=0)
        .reindex(columns=range(4), fill_value=0)
    )
    counts.columns = CLASS_COUNT_COLUMNS
    # Aggregate the complete-tile denominator across all four predicted classes.
    counts["total_256_tiles"] = counts[list(CLASS_COUNT_COLUMNS)].sum(axis=1)
    zero_psm = counts["psm_count"] == 0
    if zero_psm.any():
        groups = counts.index[zero_psm].astype(str).tolist()
        raise ZeroDivisionError(f"PSM tile count is zero for: {groups[0]}")
    # Compute the literal PIR-to-PSM count ratio after enforcing a nonzero divisor.
    counts["isr"] = counts["pir_count"] / counts["psm_count"]
    return counts.reset_index()

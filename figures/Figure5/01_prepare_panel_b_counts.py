#!/usr/bin/env python3
"""Prepare the Figure 5B training-count table.

Target/purpose: report canonical four-class tile and grouped-split counts.
Inputs: canonical tile, development, and internal-holdout CSV manifests selected
in ``CONFIG``. Outputs: one Figure5B TSV selected in ``CONFIG``. Ordered
workflow: read and validate manifests, verify split coverage and group
disjointness, build the training-count table, then write the TSV.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


CODE_ROOT = Path(__file__).resolve().parents[2]
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from src.spmap.splits import assert_figure5_training_contract, assert_group_disjoint


CLASS_NAMES = {0: "TC", 1: "PIR", 2: "PSM", 3: "OTHER_PT"}

# Edit these paths before running this interactive Figure 5 entry point.
CONFIG = {
    "tile_manifest": Path("/path/to/tile_labels.csv"),
    "development": Path("/path/to/development.csv"),
    "holdout": Path("/path/to/internal_holdout.csv"),
    "output": Path("/path/to/Figure5B_training_count_provenance.tsv"),
}


def _require_columns(frame: pd.DataFrame, columns: set[str], path: Path) -> None:
    """Fail when a manifest lacks required columns.

    Parameters: ``frame`` is the loaded table, ``columns`` is the required
    column-name set, and ``path`` identifies the input in an error message.
    Returns: no value; raises ``ValueError`` for a missing column.
    """
    # Treat the declared columns as the input boundary for downstream counting.
    missing = columns - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")


def _normalize_manifest(path: Path) -> pd.DataFrame:
    """Load and validate the canonical tile-label manifest.

    Parameters: ``path`` is the canonical tile-label CSV path.
    Returns: a normalized dataframe with tile ID, class ID/name, and parent ID.
    """
    # Map the canonical manifest to the four fields used for split reconciliation.
    frame = pd.read_csv(path)
    _require_columns(frame, {"tile_id", "class_id", "class_name", "parent_id"}, path)
    frame = frame[["tile_id", "class_id", "class_name", "parent_id"]].copy()
    frame["tile_id"] = frame["tile_id"].astype(str)
    frame["parent_id"] = frame["parent_id"].astype(str)
    frame["class_id"] = frame["class_id"].astype(int)
    if frame["tile_id"].duplicated().any():
        raise ValueError("Canonical tile_id values must be unique")
    if frame["parent_id"].eq("").any():
        raise ValueError("Canonical parent_id values must be non-empty")
    expected_names = frame["class_id"].map(CLASS_NAMES)
    if expected_names.isna().any() or not expected_names.eq(frame["class_name"]).all():
        raise ValueError("Canonical class_id and class_name values do not match the four-class contract")
    return frame


def _normalize_split(path: Path, label: str) -> pd.DataFrame:
    """Load and normalize one specified grouped split manifest.

    Parameters: ``path`` is a split CSV and ``label`` identifies its role.
    Returns: a dataframe standardized to tile ID, class ID, and parent ID.
    """
    # Standardize both split roles to the same tile, class, and parent-group schema.
    frame = pd.read_csv(path)
    _require_columns(frame, {"tile_id", "class_id", "parent_id"}, path)
    frame = frame[["tile_id", "class_id", "parent_id"]].copy()
    frame["tile_id"] = frame["tile_id"].astype(str)
    frame["parent_id"] = frame["parent_id"].astype(str)
    frame["class_id"] = frame["class_id"].astype(int)
    if frame["tile_id"].duplicated().any():
        raise ValueError(f"{label} split contains duplicate tile IDs")
    if set(frame["class_id"]) - set(CLASS_NAMES):
        raise ValueError(f"{label} split contains invalid class IDs")
    return frame


def read_training_counts(
    tile_manifest: Path, development: Path, holdout: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read and cross-check the canonical and grouped split manifests.

    Parameters: ``tile_manifest``, ``development``, and ``holdout`` identify
    the three Figure 5B CSV inputs.
    Returns: normalized canonical, development, and holdout dataframes.
    """
    # Load each source independently before reconciling split membership to the manifest.
    manifest = _normalize_manifest(tile_manifest)
    development_frame = _normalize_split(development, "development")
    holdout_frame = _normalize_split(holdout, "internal_holdout")
    assert_group_disjoint(development_frame, holdout_frame, group_col="parent_id")

    # Development and holdout must partition the canonical tiles without group overlap.
    combined = pd.concat([development_frame, holdout_frame], ignore_index=True)
    if combined["tile_id"].duplicated().any():
        raise ValueError("Development and internal-holdout tiles must be disjoint")
    if set(combined["tile_id"]) != set(manifest["tile_id"]):
        raise ValueError("Split tile IDs must cover the canonical manifest exactly")

    expected = manifest.set_index("tile_id")[["class_id", "parent_id"]]
    observed = combined.set_index("tile_id")[["class_id", "parent_id"]]
    if not observed.sort_index().equals(expected.sort_index()):
        raise ValueError("Split class_id or parent_id values do not match the canonical manifest")

    assert_figure5_training_contract(manifest, development_frame, holdout_frame)
    return manifest, development_frame, holdout_frame


def build_provenance(
    manifest: pd.DataFrame,
    development: pd.DataFrame,
    holdout: pd.DataFrame,
    *,
    tile_manifest_path: Path,
    development_path: Path,
    holdout_path: Path,
) -> pd.DataFrame:
    """Create the Figure 5B training-count table.

    Parameters: ``manifest``, ``development``, and ``holdout`` are validated
    dataframes; the three ``*_path`` values preserve their source locations.
    Returns: the table written as the Figure 5B TSV.
    """
    # Build one table that distinguishes class counts from split tile and group counts.
    rows: list[dict[str, object]] = []
    for class_name in CLASS_NAMES.values():
        rows.append(
            {
                "count_group": "training_class_tiles",
                "category": class_name,
                "count": int((manifest["class_name"] == class_name).sum()),
                "analysis_unit": "tile",
                "source_path": str(tile_manifest_path),
                "selection_rule": "canonical tile manifest; fixed four-class label",
            }
        )
    rows.append(
        {
            "count_group": "training_total_tiles",
            "category": "all_classes",
            "count": len(manifest),
            "analysis_unit": "tile",
            "source_path": str(tile_manifest_path),
            "selection_rule": "all canonical tile rows",
        }
    )
    for split_name, frame, path in (
        ("development", development, development_path),
        ("internal_holdout", holdout, holdout_path),
    ):
        rows.extend(
            [
                {
                    "count_group": "split_tiles",
                    "category": split_name,
                    "count": len(frame),
                    "analysis_unit": "tile",
                    "source_path": str(path),
                    "selection_rule": "stratified grouped split manifest",
                },
                {
                    "count_group": "split_parent_groups",
                    "category": split_name,
                    "count": frame["parent_id"].nunique(),
                    "analysis_unit": "parent_id group",
                    "source_path": str(path),
                    "selection_rule": "unique explicit parent_id values; disjoint across splits",
                },
            ]
        )
    return pd.DataFrame(rows)


def run(config: dict[str, Path]) -> None:
    """Generate the Figure 5B training-count TSV from editable configuration.

    Parameters: ``config`` supplies ``tile_manifest``, ``development``,
    ``holdout``, and ``output`` paths.
    Returns: no value; writes the configured TSV after all contract checks pass.
    """
    # Orchestrate input reconciliation, table construction, and final TSV export.
    manifest, development, holdout = read_training_counts(
        config["tile_manifest"], config["development"], config["holdout"]
    )
    output = build_provenance(
        manifest,
        development,
        holdout,
        tile_manifest_path=config["tile_manifest"],
        development_path=config["development"],
        holdout_path=config["holdout"],
    )
    config["output"].parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(config["output"], sep="\t", index=False)
    print(config["output"])


if __name__ == "__main__":
    run(CONFIG)

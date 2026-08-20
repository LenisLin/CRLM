"""Build explicit structure-grouped SpMap development and fold splits.

Purpose:
    Validate a single-class structure manifest and create the canonical 80/20
    development/common-holdout split with five grouped development folds.
Figure 5 callers:
    Figure 5 training preparation uses ``build_grouped_splits`` and
    ``assert_figure5_training_contract``; the Figure 5B count-preparation entry
    point uses the same fixed-count contract and ``assert_group_disjoint``.
Inputs:
    Tile manifests containing unique tile IDs, canonical class IDs, and
    structure-level parent IDs.
Outputs:
    Validated development and common-holdout DataFrames plus ordered fold
    train/validation DataFrames.
Ordered use:
    Build the splits before feature-store construction, then train models from
    the returned fold manifests while retaining the common holdout separately.
"""

from __future__ import annotations

import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold, StratifiedShuffleSplit


FIGURE5_CLASS_COUNTS = {0: 23894, 1: 9690, 2: 9060, 3: 18984}
FIGURE5_TOTAL_TILES = 61628
FIGURE5_DEVELOPMENT_TILES = 49001
FIGURE5_INTERNAL_HOLDOUT_TILES = 12627


def _validate_manifest(
    frame: pd.DataFrame,
    *,
    tile_col: str,
    class_col: str,
    group_col: str,
) -> pd.DataFrame:
    """Validate and type-normalize the manifest used for grouped splitting.

    Parameters:
        frame: Tile manifest.
        tile_col: Column containing unique tile identifiers.
        class_col: Column containing canonical integer class IDs.
        group_col: Column defining the structure-level split groups.

    Returns:
        A copied, normalized manifest ready for grouped split construction.
    """
    # Require the identifiers needed to keep tiles and structures traceable.
    required = {tile_col, class_col, group_col}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    # Normalize identifier and class types on a copy used only for splitting.
    manifest = frame.copy()
    manifest[tile_col] = manifest[tile_col].astype(str)
    manifest[group_col] = manifest[group_col].astype(str)
    manifest[class_col] = manifest[class_col].astype(int)
    # Enforce unique tiles, canonical classes, and one class per split group.
    if manifest[tile_col].duplicated().any():
        raise ValueError(f"{tile_col} values must be unique")
    invalid = sorted(set(manifest[class_col]) - set(range(4)))
    if invalid:
        raise ValueError(f"Invalid class IDs: {invalid}")
    group_class_counts = manifest.groupby(group_col)[class_col].nunique()
    mixed = group_class_counts[group_class_counts != 1]
    if not mixed.empty:
        raise ValueError(f"Groups must have one class; mixed group: {mixed.index[0]}")
    return manifest


def assert_group_disjoint(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    group_col: str = "parent_id",
) -> None:
    """Raise when two manifest partitions share a structure-level group.

    Parameters:
        left: First partition to compare.
        right: Second partition to compare.
        group_col: Column defining structure-level groups.

    Returns:
        ``None`` when the partitions have no overlapping group values.
    """
    # Compare structure identities rather than tile rows across partitions.
    overlap = set(left[group_col]) & set(right[group_col])
    if overlap:
        raise ValueError(f"Found {len(overlap)} overlapping {group_col} values")


def assert_figure5_training_contract(
    manifest: pd.DataFrame,
    development: pd.DataFrame,
    internal_holdout: pd.DataFrame,
    *,
    class_col: str = "class_id",
) -> None:
    """Enforce the canonical Figure 5 tile and partition counts.

    Parameters:
        manifest: Complete canonical four-class tile manifest.
        development: Grouped development partition.
        internal_holdout: Grouped common internal-holdout partition.
        class_col: Column containing canonical integer class IDs.

    Returns:
        ``None`` when total, class-specific, development, and holdout counts
        match the fixed Figure 5 contract.
    """
    # Match the full manifest to the fixed total and per-class tile counts.
    if class_col not in manifest.columns:
        raise ValueError(f"Manifest is missing column: {class_col}")
    class_counts = manifest[class_col].astype(int).value_counts().sort_index().to_dict()
    if len(manifest) != FIGURE5_TOTAL_TILES or class_counts != FIGURE5_CLASS_COUNTS:
        raise ValueError(
            "Canonical training counts do not match the fixed Figure 5 contract: "
            f"expected total={FIGURE5_TOTAL_TILES}, classes={FIGURE5_CLASS_COUNTS}; "
            f"observed total={len(manifest)}, classes={class_counts}"
        )
    # Match the grouped development/holdout sizes to the Figure 5 partition.
    observed_partitions = {
        "development": len(development),
        "internal_holdout": len(internal_holdout),
    }
    expected_partitions = {
        "development": FIGURE5_DEVELOPMENT_TILES,
        "internal_holdout": FIGURE5_INTERNAL_HOLDOUT_TILES,
    }
    if observed_partitions != expected_partitions:
        raise ValueError(
            "Grouped split tile counts do not match the fixed Figure 5 contract: "
            f"expected {expected_partitions}; observed {observed_partitions}"
        )


def build_grouped_splits(
    frame: pd.DataFrame,
    *,
    tile_col: str = "tile_id",
    class_col: str = "class_id",
    group_col: str = "parent_id",
    holdout_fraction: float = 0.20,
    n_folds: int = 5,
    seed: int = 42,
) -> dict[str, object]:
    """Build the canonical grouped 80/20 split and five development folds.

    Parameters:
        frame: Tile manifest with unique tiles and single-class parent groups.
        tile_col: Column containing unique tile identifiers.
        class_col: Column containing canonical integer class IDs.
        group_col: Column defining structure-level split groups.
        holdout_fraction: Fraction of groups reserved for common internal holdout.
        n_folds: Number of grouped folds constructed within development.
        seed: Random state used by both splitters.

    Returns:
        A mapping containing ``development``, ``internal_holdout``, and ordered
        fold mappings with ``fold``, ``train``, and ``validation`` entries.
    """
    # Normalize the tile manifest before deriving structure-level split rows.
    manifest = _validate_manifest(
        frame, tile_col=tile_col, class_col=class_col, group_col=group_col
    )
    groups = manifest[[group_col, class_col]].drop_duplicates().reset_index(drop=True)
    # Reserve the common holdout by stratifying unique single-class structures.
    splitter = StratifiedShuffleSplit(
        n_splits=1,
        test_size=holdout_fraction,
        random_state=seed,
    )
    development_index, holdout_index = next(
        splitter.split(groups[group_col], groups[class_col])
    )
    # Project selected structure IDs back to complete tile-level partitions.
    development_groups = set(groups.iloc[development_index][group_col])
    holdout_groups = set(groups.iloc[holdout_index][group_col])
    development = manifest[manifest[group_col].isin(development_groups)].reset_index(drop=True)
    holdout = manifest[manifest[group_col].isin(holdout_groups)].reset_index(drop=True)
    assert_group_disjoint(development, holdout, group_col=group_col)

    # Construct stratified, structure-disjoint folds only within development.
    fold_splitter = StratifiedGroupKFold(
        n_splits=n_folds,
        shuffle=True,
        random_state=seed,
    )
    folds = []
    for fold, (train_index, validation_index) in enumerate(
        fold_splitter.split(
            development[tile_col],
            development[class_col],
            groups=development[group_col],
        ),
        start=1,
    ):
        # Preserve fold order while materializing disjoint train/validation rows.
        train = development.iloc[train_index].reset_index(drop=True)
        validation = development.iloc[validation_index].reset_index(drop=True)
        assert_group_disjoint(train, validation, group_col=group_col)
        folds.append({"fold": fold, "train": train, "validation": validation})
    return {"development": development, "internal_holdout": holdout, "folds": folds}

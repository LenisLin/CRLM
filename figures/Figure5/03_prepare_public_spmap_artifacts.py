#!/usr/bin/env python3
"""Prepare SpMap labels, partitions, metrics, predictions, and model archives.

Target/purpose: convert the selected SpMap C10 five-seed/five-fold artifacts to
portable public files. Inputs: canonical tile and split manifests, the five
selected run directories, and canonical CONCH feature shards set in ``CONFIG``.
Outputs: two compressed TSVs, two metric TSVs, and three tar.gz archives in
``output_dir``. Ordered workflow: reconcile labels and grouped partitions,
collect saved metrics, validate each seed's out-of-fold coverage, write portable
model metadata, and package the selected artifacts.
"""

from __future__ import annotations

import io
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd


CLASS_NAMES = {0: "TC", 1: "PIR", 2: "PSM", 3: "OTHER_PT"}
SEEDS = (24, 101, 202, 303, 404)
FOLDS = (1, 2, 3, 4, 5)

# Edit these paths before running this interactive Figure 5 entry point.
CONFIG = {
    "tile_manifest": Path("/path/to/tile_labels_canonical.csv"),
    "development": Path("/path/to/development.csv"),
    "holdout": Path("/path/to/internal_holdout.csv"),
    "validation_fold_dir": Path("/path/to/five_fold"),
    "seed_runs": {
        24: Path("/path/to/C10"),
        101: Path("/path/to/seed101"),
        202: Path("/path/to/seed202"),
        303: Path("/path/to/seed303"),
        404: Path("/path/to/seed404"),
    },
    "conch_features": Path("/path/to/03_conch_features"),
    "output_dir": Path("/path/to/derived_results/SpMap"),
}


def _require_columns(frame: pd.DataFrame, columns: set[str], path: Path) -> None:
    """Fail when an input table lacks required columns.

    Parameters: ``frame`` is the loaded table, ``columns`` is the required set,
    and ``path`` identifies the source in an error message.
    Returns: no value; raises ``ValueError`` for missing columns.
    """
    # Validate the source boundary before normalizing legacy field names.
    missing = columns - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")


def read_tile_manifest(path: Path) -> pd.DataFrame:
    """Read the canonical four-class tile-label manifest.

    Parameters: ``path`` identifies the canonical label CSV.
    Returns: a validated dataframe with tile, class, and parent identifiers.
    """
    # Retain only the public fields required for split and prediction matching.
    frame = pd.read_csv(path)
    _require_columns(frame, {"tile_id", "class_id", "class_name", "parent_id"}, path)
    frame = frame[["tile_id", "class_id", "class_name", "parent_id"]].copy()
    frame["tile_id"] = frame["tile_id"].astype(str)
    frame["parent_id"] = frame["parent_id"].astype(str)
    frame["class_id"] = frame["class_id"].astype(int)
    if len(frame) != 61628 or frame["tile_id"].duplicated().any():
        raise ValueError("Canonical tile manifest must contain 61,628 unique tiles")
    expected_names = frame["class_id"].map(CLASS_NAMES)
    if expected_names.isna().any() or not expected_names.eq(frame["class_name"]).all():
        raise ValueError("Tile class IDs and names do not match the four-class contract")
    return frame


def read_split(path: Path) -> pd.DataFrame:
    """Read one legacy or normalized SpMap split manifest.

    Parameters: ``path`` identifies a development, holdout, or validation CSV.
    Returns: a dataframe normalized to tile_id, class_id, and parent_id.
    """
    # Current saved splits use fileName/type_num; public exports use canonical names.
    frame = pd.read_csv(path)
    if {"fileName", "type_num", "parent_id"}.issubset(frame.columns):
        frame = frame.rename(columns={"fileName": "tile_id", "type_num": "class_id"})
    _require_columns(frame, {"tile_id", "class_id", "parent_id"}, path)
    frame = frame[["tile_id", "class_id", "parent_id"]].copy()
    frame["tile_id"] = frame["tile_id"].astype(str)
    frame["parent_id"] = frame["parent_id"].astype(str)
    frame["class_id"] = frame["class_id"].astype(int)
    if frame["tile_id"].duplicated().any():
        raise ValueError(f"{path} contains duplicate tile IDs")
    return frame


def build_partitions(
    manifest: pd.DataFrame,
    development_path: Path,
    holdout_path: Path,
    validation_fold_dir: Path,
) -> pd.DataFrame:
    """Combine primary and five-fold memberships in one tile partition table.

    Parameters: ``manifest`` is canonical; the remaining paths identify saved
    development, holdout, and validation-fold manifests.
    Returns: all 61,628 tiles with primary_partition and validation_fold.
    """
    # The primary split partitions every canonical tile exactly once.
    development = read_split(development_path)
    holdout = read_split(holdout_path)
    combined = pd.concat(
        [
            development.assign(primary_partition="development"),
            holdout.assign(primary_partition="internal_holdout"),
        ],
        ignore_index=True,
    )
    if len(development) != 49001 or len(holdout) != 12627:
        raise ValueError("Expected 49,001 development and 12,627 holdout tiles")
    if combined["tile_id"].duplicated().any() or set(combined["tile_id"]) != set(manifest["tile_id"]):
        raise ValueError("Development and holdout must partition the canonical tiles")
    if set(development["parent_id"]) & set(holdout["parent_id"]):
        raise ValueError("Development and holdout parent groups must be disjoint")

    # Each development tile appears in one saved validation fold.
    fold_rows = []
    for fold in FOLDS:
        frame = read_split(validation_fold_dir / f"val_fold{fold}.csv")
        fold_rows.append(frame.assign(validation_fold=fold))
    folds = pd.concat(fold_rows, ignore_index=True)
    if folds["tile_id"].duplicated().any() or set(folds["tile_id"]) != set(development["tile_id"]):
        raise ValueError("Five validation folds must partition the development tiles")

    result = manifest.merge(
        combined[["tile_id", "class_id", "parent_id", "primary_partition"]],
        on="tile_id",
        how="left",
        suffixes=("", "_split"),
        validate="one_to_one",
    )
    if not result["class_id"].eq(result["class_id_split"]).all() or not result["parent_id"].eq(result["parent_id_split"]).all():
        raise ValueError("Primary split labels do not match the canonical manifest")
    result = result.drop(columns=["class_id_split", "parent_id_split"]).merge(
        folds[["tile_id", "validation_fold"]],
        on="tile_id",
        how="left",
        validate="one_to_one",
    )
    result["validation_fold"] = result["validation_fold"].astype("Int64")
    if result.loc[result["primary_partition"] == "development", "validation_fold"].isna().any():
        raise ValueError("Every development tile must have a validation-fold assignment")
    if result.loc[result["primary_partition"] == "internal_holdout", "validation_fold"].notna().any():
        raise ValueError("Internal-holdout tiles must not have a validation fold")
    return result


def read_seed_metrics(seed_runs: dict[int, Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collect saved class-specific and overall metrics for five seeds.

    Parameters: ``seed_runs`` maps the specified seeds to selected run roots.
    Returns: tidy performance rows and the seed-24 long confusion table.
    """
    # Convert saved wide summaries to one inspectable metric-value table.
    metric_rows: list[dict[str, object]] = []
    seed24_counts: pd.DataFrame | None = None
    seed24_rows: pd.DataFrame | None = None
    for seed in SEEDS:
        root = seed_runs[seed] / "comment12_summary"
        class_frame = pd.read_csv(root / "pooled_oof_per_class_metrics.csv")
        overall = pd.read_csv(root / "pooled_oof_metrics.csv")
        _require_columns(
            class_frame,
            {"class_id", "class_name", "support", "precision", "recall", "f1", "auroc_ovr", "auprc_ovr"},
            root / "pooled_oof_per_class_metrics.csv",
        )
        if len(class_frame) != 4 or len(overall) != 1:
            raise ValueError(f"Seed {seed} metric summaries have unexpected row counts")
        for row in class_frame.itertuples(index=False):
            for metric in ("precision", "recall", "f1", "auroc_ovr", "auprc_ovr"):
                metric_rows.append(
                    {
                        "seed": seed,
                        "scope": "class",
                        "class_name": row.class_name,
                        "support": int(row.support),
                        "metric": metric,
                        "value": float(getattr(row, metric)),
                    }
                )
        overall_metrics = {
            "accuracy": "accuracy",
            "precision_macro": "macro_precision",
            "recall_macro": "macro_recall",
            "f1_macro": "macro_f1",
            "auroc_macro_ovr": "macro_auroc_ovr",
            "auprc_macro_ovr": "macro_auprc_ovr",
            "auprc_pir_psm_mean": "pir_psm_mean_auprc",
        }
        for source, public_name in overall_metrics.items():
            metric_rows.append(
                {
                    "seed": seed,
                    "scope": "overall",
                    "class_name": "ALL",
                    "support": int(overall.iloc[0]["rows"]),
                    "metric": public_name,
                    "value": float(overall.iloc[0][source]),
                }
            )
        if seed == 24:
            seed24_counts = pd.read_csv(root / "pooled_oof_confusion_counts.csv", index_col=0)
            seed24_rows = pd.read_csv(root / "pooled_oof_confusion_row_normalized.csv", index_col=0)

    if seed24_counts is None or seed24_rows is None:
        raise ValueError("Seed 24 confusion summaries were not loaded")
    if list(seed24_counts.index) != list(CLASS_NAMES.values()) or list(seed24_counts.columns) != list(CLASS_NAMES.values()):
        raise ValueError("Seed 24 confusion table does not use the canonical class order")
    confusion_rows = []
    for true_class in CLASS_NAMES.values():
        for predicted_class in CLASS_NAMES.values():
            confusion_rows.append(
                {
                    "seed": 24,
                    "true_class": true_class,
                    "predicted_class": predicted_class,
                    "count": int(seed24_counts.loc[true_class, predicted_class]),
                    "true_class_fraction": float(seed24_rows.loc[true_class, predicted_class]),
                }
            )
    return pd.DataFrame(metric_rows), pd.DataFrame(confusion_rows)


def validate_oof_predictions(seed_runs: dict[int, Path], development: pd.DataFrame) -> None:
    """Validate complete out-of-fold coverage for every selected seed.

    Parameters: ``seed_runs`` maps seeds to run roots and ``development`` is the
    public partition table restricted to development tiles.
    Returns: no value; raises ``ValueError`` for incomplete or inconsistent OOF data.
    """
    # Each seed's five validation files must cover every development tile once.
    expected = development.set_index("tile_id")[["class_id", "parent_id"]].sort_index()
    for seed in SEEDS:
        frames = []
        for fold in FOLDS:
            path = seed_runs[seed] / f"fold{fold}" / "validation_predictions.csv"
            frame = pd.read_csv(path)
            _require_columns(
                frame,
                {"fileName", "type_num", "parent_id", "pred_class", "prob_tc", "prob_pir", "prob_psm", "prob_other_pt"},
                path,
            )
            frames.append(frame.rename(columns={"fileName": "tile_id", "type_num": "class_id"}))
        pooled = pd.concat(frames, ignore_index=True)
        if len(pooled) != 49001 or pooled["tile_id"].duplicated().any():
            raise ValueError(f"Seed {seed} must provide 49,001 unique OOF predictions")
        observed = pooled.set_index("tile_id")[["class_id", "parent_id"]].sort_index()
        if not observed.equals(expected):
            raise ValueError(f"Seed {seed} OOF labels do not match the development partition")
        probabilities = pooled[["prob_tc", "prob_pir", "prob_psm", "prob_other_pt"]].to_numpy(dtype=float)
        if not np.isfinite(probabilities).all() or not np.allclose(probabilities.sum(axis=1), 1, atol=1e-6):
            raise ValueError(f"Seed {seed} OOF probabilities must be finite and sum to one")


def _add_bytes(archive: tarfile.TarFile, arcname: str, content: bytes) -> None:
    """Add generated bytes to a tar archive under a portable name.

    Parameters: ``archive`` is open for writing, ``arcname`` is the member name,
    and ``content`` is the serialized payload.
    Returns: no value after adding one archive member.
    """
    # Generated manifests are written directly without temporary filesystem files.
    member = tarfile.TarInfo(arcname)
    member.size = len(content)
    archive.addfile(member, io.BytesIO(content))


def public_model_manifest(path: Path, seed: int, fold: int) -> dict[str, object]:
    """Select portable training metadata from one saved run manifest.

    Parameters: ``path`` identifies a private saved manifest; ``seed`` and
    ``fold`` identify the public model member.
    Returns: a dictionary containing portable model and training fields.
    """
    # Preserve model, optimization, performance, and environment metadata only.
    source = json.loads(path.read_text(encoding="utf-8"))
    data_summary = {
        role: {key: value for key, value in details.items() if key in {"rows", "parent_ids"}}
        for role, details in source.get("data", {}).items()
        if isinstance(details, dict)
    }
    keys = (
        "feature_mode",
        "input_dim",
        "class_map",
        "architecture",
        "activation",
        "dropout",
        "depth",
        "batch_norm",
        "input_normalization",
        "loss",
        "class_weight_mode",
        "target_multiplier",
        "aux_pir_psm_binary_weight",
        "aux_target_ovr_weight",
        "sampler",
        "focal_gamma",
        "optimizer",
        "parameters",
        "selected_epoch",
        "best_selection_value",
        "validation_metrics",
        "runtime_seconds",
        "environment",
        "feature_preload_seconds",
    )
    result = {key: source.get(key) for key in keys}
    result.update({"seed": seed, "fold": fold, "data": data_summary})
    return result


def write_prediction_archive(seed_runs: dict[int, Path], output: Path) -> None:
    """Archive the 25 selected validation-prediction tables.

    Parameters: ``seed_runs`` maps seeds to selected run roots and ``output`` is
    the target tar.gz path.
    Returns: no value after writing the archive.
    """
    # Seed/fold member paths retain the public evaluation hierarchy.
    with tarfile.open(output, "w:gz") as archive:
        for seed in SEEDS:
            for fold in FOLDS:
                source = seed_runs[seed] / f"fold{fold}" / "validation_predictions.csv"
                archive.add(source, arcname=f"seed{seed}/fold{fold}/validation_predictions.csv")


def write_weight_archive(seed_runs: dict[int, Path], output: Path) -> None:
    """Archive 25 selected checkpoints with portable model manifests.

    Parameters: ``seed_runs`` maps seeds to selected run roots and ``output`` is
    the target tar.gz path.
    Returns: no value after writing the archive.
    """
    # Every checkpoint is paired with its portable model and training metadata.
    with tarfile.open(output, "w:gz") as archive:
        for seed in SEEDS:
            for fold in FOLDS:
                fold_root = seed_runs[seed] / f"fold{fold}"
                prefix = f"seed{seed}/fold{fold}"
                archive.add(fold_root / "best.pt", arcname=f"{prefix}/best.pt")
                manifest = public_model_manifest(fold_root / "run_manifest.json", seed, fold)
                content = json.dumps(manifest, indent=2, sort_keys=True).encode("utf-8")
                _add_bytes(archive, f"{prefix}/model_manifest.json", content)


def write_feature_archive(source_root: Path, output: Path) -> None:
    """Archive canonical CONCH tar shards and member inventories.

    Parameters: ``source_root`` is the canonical feature directory and
    ``output`` is the target tar.gz path.
    Returns: no value after writing feature files under portable member paths.
    """
    # Package the actual feature shards and text member lists in sorted order.
    sources = sorted(
        path
        for path in source_root.rglob("*")
        if path.is_file() and (path.suffix == ".tar" or path.name == "archive_members.txt")
    )
    if not sources:
        raise ValueError("No canonical CONCH feature shards were found")
    with tarfile.open(output, "w:gz") as archive:
        for source in sources:
            archive.add(source, arcname=str(Path("SpMap_CONCH_features") / source.relative_to(source_root)))


def run(config: dict[str, object]) -> None:
    """Write all public SpMap artifacts from editable configuration.

    Parameters: ``config`` supplies label, split, run, feature, and output paths.
    Returns: no value; writes seven public release artifacts.
    """
    # Reconcile public labels and split membership before packaging model outputs.
    seed_runs = config["seed_runs"]
    if set(seed_runs) != set(SEEDS):
        raise ValueError("seed_runs must define seeds 24, 101, 202, 303, and 404")
    manifest = read_tile_manifest(config["tile_manifest"])
    partitions = build_partitions(
        manifest,
        config["development"],
        config["holdout"],
        config["validation_fold_dir"],
    )
    performance, confusion = read_seed_metrics(seed_runs)
    development = partitions.loc[partitions["primary_partition"] == "development"]
    validate_oof_predictions(seed_runs, development)

    output_dir = config["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(
        output_dir / "SpMap_reference_tile_labels.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    partitions.to_csv(
        output_dir / "SpMap_tile_partitions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    performance.to_csv(output_dir / "SpMap_model_performance_5seeds.tsv", sep="\t", index=False)
    confusion.to_csv(output_dir / "SpMap_confusion_matrix.tsv", sep="\t", index=False)
    write_prediction_archive(
        seed_runs, output_dir / "SpMap_OOF_predictions_C10_5seeds_5folds.tar.gz"
    )
    write_weight_archive(
        seed_runs, output_dir / "SpMap_model_weights_C10_5seeds_5folds.tar.gz"
    )
    write_feature_archive(
        config["conch_features"], output_dir / "SpMap_CONCH_features.tar.gz"
    )
    for path in sorted(output_dir.iterdir()):
        print(path)


if __name__ == "__main__":
    run(CONFIG)

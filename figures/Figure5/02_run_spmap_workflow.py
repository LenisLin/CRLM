#!/usr/bin/env python3
"""Run the Figure 5 SpMap workflow from editable configuration.

Target/purpose: call the shared SpMap reference, preprocessing, CONCH,
grouped split/training/evaluation, WSI-inference, and ISR APIs. Inputs: the
paths and adapter callables set in ``CONFIG``. Outputs: selected manifests,
paired reference features, fold checkpoints/metrics, and WSI prediction/ISR
tables. Ordered workflow: build reference tiles, preprocess and extract CONCH
features, build grouped splits and train/evaluate folds, then infer WSI tiles
and aggregate ISR. The ``enabled`` mapping selects the stages for each run.
"""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
import torch


CODE_ROOT = Path(__file__).resolve().parents[2]
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from src.spmap.conch_features import ConchFeatureExtractor
from src.spmap.inference import load_checkpoint, predict_wsi_tiles
from src.spmap.isr import aggregate_isr
from src.spmap.mlp import (
    FeatureDataset,
    InMemoryFeatureStore,
    checkpoint_manifest,
    evaluate_predictions,
    prediction_frame,
    require_cuda,
    set_seed,
    summarize_oof_predictions,
    summarize_training_seed_variability,
    train_fold,
)
from src.spmap.preprocessing import VahadaneTrainingNormalizer, preprocess_training_tile
from src.spmap.reference_tiles import (
    build_complete_tile_manifest,
    build_structure_bboxes,
    extract_tile,
)
from src.spmap.splits import assert_figure5_training_contract, build_grouped_splits
from src.spmap.wsi import iter_wsi_tiles, open_wsi


SPLIT_SEED = 42
TRAINING_SEEDS = (24, 101, 202, 303, 404)
FOLDS = 5
EPOCHS = 100
BATCH_SIZE = 256
NUM_WORKERS = 6
GPU_INDEX = 0

# Set paths and adapter callables, then select workflow stages with ``enabled``.
CONFIG: dict[str, Any] = {
    "enabled": {
        "build_reference_manifest": False,
        "extract_reference_conch": False,
        "train_evaluate": False,
        "infer_wsi": False,
    },
    "reference": {
        "structure_cells_csv": None,
        "selected_regions_csv": None,
        "image_sizes": {},
        "manifest_output": Path("/path/to/reference_tile_manifest.csv"),
        "image_loader": None,
        "conch_checkpoint": None,
        "conch_model_loader": None,
        "conch_preprocess": None,
        "conch_forward_representations": None,
        "feature_output": Path("/path/to/reference_conch_features.npz"),
    },
    "training": {
        "tile_manifest_csv": Path("/path/to/canonical_tile_manifest.csv"),
        "feature_tar_root": Path("/path/to/03_conch_features"),
        "output_dir": Path("/path/to/spmap_training"),
    },
    "inference": {
        "wsi_paths": [],
        "checkpoint": Path("/path/to/checkpoint.pt"),
        "model_id": "replace_with_checkpoint_identifier",
        "output_dir": Path("/path/to/spmap_wsi_inference"),
        "wsi_level": 0,
        "macenko_normalizer": None,
        "conch_checkpoint": None,
        "conch_model_loader": None,
        "conch_preprocess": None,
        "conch_forward_representations": None,
    },
}


def resolve_callable(specification: str | Callable[..., Any] | None, name: str) -> Callable[..., Any]:
    """Resolve a configured adapter callable or fail with its configuration name.

    Parameters: ``specification`` is either a callable or a ``module:attribute``
    string, and ``name`` identifies the required configuration field.
    Returns: the configured callable.
    """
    # Normalize direct callables and import specifications to one adapter interface.
    if callable(specification):
        return specification
    if isinstance(specification, str) and ":" in specification:
        module_name, attribute_name = specification.split(":", maxsplit=1)
        resolved = getattr(importlib.import_module(module_name), attribute_name)
        if callable(resolved):
            return resolved
    raise ValueError(f"{name} must be a callable or 'module:attribute' string")


def require_file(path: str | Path | None, name: str) -> Path:
    """Return an existing configured file path.

    Parameters: ``path`` is a configured path-like value and ``name`` labels it.
    Returns: the existing path as a ``Path`` instance.
    """
    # Convert configured file values once so downstream stages receive concrete paths.
    if path is None:
        raise ValueError(f"{name} is required")
    resolved = Path(path)
    if not resolved.is_file():
        raise FileNotFoundError(f"{name} does not exist: {resolved}")
    return resolved


def require_directory(path: str | Path | None, name: str) -> Path:
    """Return an existing configured directory path.

    Parameters: ``path`` is a configured path-like value and ``name`` labels it.
    Returns: the existing directory as a ``Path`` instance.
    """
    # Convert configured directory values once at the workflow boundary.
    if path is None:
        raise ValueError(f"{name} is required")
    resolved = Path(path)
    if not resolved.is_dir():
        raise NotADirectoryError(f"{name} is not a directory: {resolved}")
    return resolved


def prepare_empty_output_directory(path: str | Path, name: str) -> Path:
    """Create or validate an empty output directory for a new run.

    Parameters: ``path`` is the configured output location and ``name`` labels
    the directory in error messages.
    Returns: the empty directory as a ``Path`` instance.
    """
    # Each training or inference invocation owns one initially empty output directory.
    output = Path(path)
    if output.exists():
        if not output.is_dir():
            raise NotADirectoryError(f"{name} is not a directory: {output}")
        if next(output.iterdir(), None) is not None:
            raise FileExistsError(f"{name} must be empty for a new run: {output}")
    else:
        output.mkdir(parents=True)
    return output


def write_prediction_artifacts(
    predictions: pd.DataFrame,
    output_dir: Path,
    prefix: str,
) -> dict[str, object]:
    """Write classification, per-class, and confusion artifacts for predictions.

    Parameters: ``predictions`` is one fold's validated prediction table,
    ``output_dir`` is its seed/fold directory, and ``prefix`` identifies the
    validation or common-holdout evaluation.
    Returns: the classification metric mapping written for the evaluation.
    """
    # Derive all metric and confusion exports from the same prediction rows.
    summary = summarize_oof_predictions([predictions])
    pd.DataFrame([summary["metrics"]]).to_csv(
        output_dir / f"{prefix}_classification_metrics.csv", index=False
    )
    summary["per_class_metrics"].to_csv(
        output_dir / f"{prefix}_per_class_metrics.csv", index=False
    )
    summary["confusion_counts"].to_csv(
        output_dir / f"{prefix}_confusion_counts.csv", index_label="true_class"
    )
    summary["confusion_row_normalized"].to_csv(
        output_dir / f"{prefix}_confusion_row_proportions.csv",
        index_label="true_class",
    )
    return summary["metrics"]


def build_reference_manifest(config: dict[str, Any]) -> pd.DataFrame:
    """Build the caller-selected SpMap reference-tile manifest.

    Parameters: ``config`` provides optional registered-cell and selected-region
    CSVs, image dimensions by source image ID, and a manifest output path.
    Returns: the complete-tile manifest written to the configured CSV.
    """
    # Map registered structures and selected regions into one reference-region stream.
    reference = config["reference"]
    region_frames: list[pd.DataFrame] = []
    cells_path = reference["structure_cells_csv"]
    if cells_path is not None:
        cells = pd.read_csv(require_file(cells_path, "structure_cells_csv"))
        image_sizes = reference["image_sizes"]
        for source_image_id, group in cells.groupby("source_image_id", sort=False):
            if str(source_image_id) not in image_sizes:
                raise ValueError(f"image_sizes lacks source image {source_image_id}")
            region_frames.append(
                build_structure_bboxes(
                    group,
                    image_size=tuple(image_sizes[str(source_image_id)]),
                )
            )
    regions_path = reference["selected_regions_csv"]
    if regions_path is not None:
        region_frames.append(pd.read_csv(require_file(regions_path, "selected_regions_csv")))
    if not region_frames:
        raise ValueError("Provide structure_cells_csv and/or selected_regions_csv")
    # Complete tiling is applied after all configured region sources share one schema.
    manifest = build_complete_tile_manifest(pd.concat(region_frames, ignore_index=True))
    output = Path(reference["manifest_output"])
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output, index=False)
    return manifest


def create_conch_extractor(config: dict[str, Any]) -> ConchFeatureExtractor:
    """Instantiate the CUDA-only CONCH adapter from editable configuration.

    Parameters: ``config`` provides a checkpoint and model/preprocess/forward
    adapter callables for a compatible CONCH v1 implementation.
    Returns: a ready CUDA ``ConchFeatureExtractor``.
    """
    # Bind the checkpoint and three CONCH adapters to the fixed workflow GPU.
    return ConchFeatureExtractor(
        model_loader=resolve_callable(config["conch_model_loader"], "conch_model_loader"),
        checkpoint=str(require_file(config["conch_checkpoint"], "conch_checkpoint")),
        preprocess=resolve_callable(config["conch_preprocess"], "conch_preprocess"),
        forward_representations=resolve_callable(
            config["conch_forward_representations"],
            "conch_forward_representations",
        ),
        gpu_index=GPU_INDEX,
    )


def extract_reference_conch(config: dict[str, Any]) -> None:
    """Preprocess selected reference tiles and save paired CONCH arrays.

    Parameters: ``config`` provides the reference manifest, RGB image loader,
    CONCH adapters, and paired-feature NPZ output path.
    Returns: no value; writes a compressed NPZ with tile IDs and 768/512 arrays.
    """
    # Initialize preprocessing and feature extraction once for the full manifest.
    reference = config["reference"]
    manifest = pd.read_csv(require_file(reference["manifest_output"], "manifest_output"))
    image_loader = resolve_callable(reference["image_loader"], "image_loader")
    normalizer = VahadaneTrainingNormalizer()
    extractor = create_conch_extractor(reference)
    images: dict[str, np.ndarray] = {}
    kept_ids: list[str] = []
    kept_images: list[np.ndarray] = []
    # Cache each source image and retain tile IDs in exact processed-image order.
    for tile in manifest.itertuples(index=False):
        source_image_id = str(tile.source_image_id)
        if source_image_id not in images:
            images[source_image_id] = np.asarray(image_loader(source_image_id))
        processed = preprocess_training_tile(
            extract_tile(images[source_image_id], pd.Series(tile._asdict())),
            normalizer=normalizer,
        )
        if processed.kept:
            kept_ids.append(str(tile.tile_id))
            kept_images.append(processed.image)
    if not kept_ids:
        raise ValueError("Reference preprocessing retained no tiles")

    # Batch extraction preserves alignment between retained IDs and both representations.
    semantic_batches, alignment_batches = [], []
    for start in range(0, len(kept_ids), BATCH_SIZE):
        features = extractor.extract(
            kept_ids[start : start + BATCH_SIZE],
            kept_images[start : start + BATCH_SIZE],
        )
        semantic_batches.append(features.semantic_768)
        alignment_batches.append(features.alignment_512)
    # Export the paired feature matrices with their shared tile-ID axis.
    output = Path(reference["feature_output"])
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        tile_ids=np.asarray(kept_ids, dtype=str),
        semantic_768=np.concatenate(semantic_batches, axis=0),
        alignment_512=np.concatenate(alignment_batches, axis=0),
    )


def train_and_evaluate(config: dict[str, Any]) -> None:
    """Run five ordered fold sequences using one RAM-resident feature store.

    Parameters: ``config`` provides the canonical manifest, CONCH tar root, and
    output directory. Fixed split/training parameters are module constants.
    Returns: no value; writes seed/fold checkpoints and evaluation artifacts,
    seed-24 S11A/B inputs, and five-seed S11C summaries.
    """
    # Create the grouped development, fold, and common-holdout evaluation structure.
    training = config["training"]
    output_dir = prepare_empty_output_directory(
        training["output_dir"], "training.output_dir"
    )
    manifest = pd.read_csv(require_file(training["tile_manifest_csv"], "tile_manifest_csv"))
    splits = build_grouped_splits(
        manifest,
        n_folds=FOLDS,
        seed=SPLIT_SEED,
    )
    assert_figure5_training_contract(
        manifest,
        splits["development"],
        splits["internal_holdout"],
    )
    splits["development"].to_csv(output_dir / "development.csv", index=False)
    splits["internal_holdout"].to_csv(output_dir / "internal_holdout.csv", index=False)

    # Load all required tar-shard features once before any seed or fold starts.
    feature_store = InMemoryFeatureStore(
        require_directory(training["feature_tar_root"], "feature_tar_root"),
        [splits["development"], splits["internal_holdout"]],
    )
    # Repeat the fixed folds across seeds while sharing the RAM-resident feature store.
    seed_summaries: dict[int, dict[str, object]] = {}
    for seed in TRAINING_SEEDS:
        seed_dir = output_dir / f"seed_{seed}"
        seed_dir.mkdir()
        set_seed(seed)
        validation_predictions = []
        validation_metrics, holdout_metrics = [], []
        # Each fold trains once, freezes its selected model, and evaluates two disjoint sets.
        for fold_data in splits["folds"]:
            fold = int(fold_data["fold"])
            fold_dir = seed_dir / f"fold_{fold}"
            fold_dir.mkdir()
            result = train_fold(
                fold_data["train"],
                fold_data["validation"],
                feature_store,
                epochs=EPOCHS,
                batch_size=BATCH_SIZE,
                seed=seed,
                gpu_index=GPU_INDEX,
                num_workers=NUM_WORKERS,
                reset_seed=False,
            )
            model_id = f"spmap_seed{seed}_fold{fold}"
            torch.save(result.model.state_dict(), fold_dir / "best_checkpoint.pt")
            result.history.to_csv(fold_dir / "training_history.csv", index=False)
            with (fold_dir / "checkpoint_manifest.json").open("x", encoding="utf-8") as handle:
                json.dump(
                    checkpoint_manifest(
                        model_id=model_id,
                        selected_epoch=result.selected_epoch,
                        fold=fold,
                        seed=seed,
                    ),
                    handle,
                    indent=2,
                    sort_keys=True,
                )
                handle.write("\n")

            # Prepare validation and common-holdout predictions from the same fold checkpoint.
            validation = prediction_frame(
                result.model,
                FeatureDataset(fold_data["validation"], feature_store),
                batch_size=BATCH_SIZE,
                gpu_index=GPU_INDEX,
                num_workers=NUM_WORKERS,
            )
            holdout = prediction_frame(
                result.model,
                FeatureDataset(splits["internal_holdout"], feature_store),
                batch_size=BATCH_SIZE,
                gpu_index=GPU_INDEX,
                num_workers=NUM_WORKERS,
            )
            validation.to_csv(fold_dir / "validation_predictions.csv", index=False)
            holdout.to_csv(fold_dir / "common_holdout_predictions.csv", index=False)
            write_prediction_artifacts(validation, fold_dir, "validation")
            write_prediction_artifacts(holdout, fold_dir, "common_holdout")
            validation_predictions.append(validation)
            validation_metrics.append({"fold": fold, **evaluate_predictions(validation)})
            holdout_metrics.append({"fold": fold, **evaluate_predictions(holdout)})

        # Consolidate fold metrics and pooled out-of-fold predictions for this seed.
        fold_summary = pd.DataFrame(validation_metrics).drop(columns="confusion_matrix")
        holdout_summary = pd.DataFrame(holdout_metrics).drop(columns="confusion_matrix")
        fold_summary.to_csv(seed_dir / "fold_validation_metrics.csv", index=False)
        holdout_summary.to_csv(seed_dir / "common_holdout_metrics.csv", index=False)
        if seed == TRAINING_SEEDS[0]:
            fold_summary.to_csv(output_dir / "fold_validation_metrics.csv", index=False)
            holdout_summary.to_csv(output_dir / "common_holdout_metrics.csv", index=False)

        oof = summarize_oof_predictions(validation_predictions)
        oof["pooled_predictions"].to_csv(seed_dir / "oof_predictions.csv", index=False)
        pd.DataFrame([oof["metrics"]]).to_csv(
            seed_dir / "oof_classification_metrics.csv", index=False
        )
        oof["per_class_metrics"].to_csv(
            seed_dir / "oof_per_class_metrics.csv", index=False
        )
        oof["confusion_counts"].to_csv(
            seed_dir / "oof_confusion_counts.csv", index_label="true_class"
        )
        oof["confusion_row_normalized"].to_csv(
            seed_dir / "oof_confusion_row_proportions.csv", index_label="true_class"
        )
        seed_summaries[seed] = oof

    # Export the prespecified seed panels and cross-seed variability tables.
    seed_24 = seed_summaries[TRAINING_SEEDS[0]]
    seed_24["confusion_counts"].to_csv(
        output_dir / "confusion_counts.csv", index_label="true_class"
    )
    seed_24["confusion_row_normalized"].to_csv(
        output_dir / "confusion_row_proportions.csv", index_label="true_class"
    )
    seed_class_metrics = pd.concat(
        [
            summary["per_class_metrics"].assign(seed=seed)
            for seed, summary in seed_summaries.items()
        ],
        ignore_index=True,
    )
    seed_class_metrics = seed_class_metrics[
        ["seed", *[column for column in seed_class_metrics.columns if column != "seed"]]
    ]
    seed_overall_metrics = pd.DataFrame(
        [{"seed": seed, **summary["metrics"]} for seed, summary in seed_summaries.items()]
    )
    seed_class_metrics.to_csv(output_dir / "seed_class_metrics.csv", index=False)
    seed_overall_metrics.to_csv(output_dir / "seed_overall_metrics.csv", index=False)
    summarize_training_seed_variability(seed_overall_metrics).to_csv(
        output_dir / "training_seed_variability.csv", index=False
    )


def infer_wsi_and_aggregate_isr(config: dict[str, Any]) -> None:
    """Predict configured WSIs and aggregate their literal PIR/PSM ISR.

    Parameters: ``config`` provides WSI paths, one checkpoint/model ID, Macenko
    and CONCH adapters, and an output directory.
    Returns: no value; writes retained-tile prediction and WSI ISR CSVs.
    """
    # Initialize one normalizer, extractor, and checkpoint for all configured WSIs.
    inference = config["inference"]
    output_dir = prepare_empty_output_directory(
        inference["output_dir"], "inference.output_dir"
    )
    paths = [require_file(path, "wsi_paths item") for path in inference["wsi_paths"]]
    if not paths:
        raise ValueError("inference.wsi_paths must contain at least one WSI")
    normalizer = resolve_callable(inference["macenko_normalizer"], "macenko_normalizer")
    extractor = create_conch_extractor(inference)
    model = load_checkpoint(require_file(inference["checkpoint"], "checkpoint"), gpu_index=GPU_INDEX)
    # Keep each slide open only while its retained tiles are generated and predicted.
    prediction_frames = []
    for path in paths:
        slide = open_wsi(path)
        try:
            prediction_frames.append(
                predict_wsi_tiles(
                    iter_wsi_tiles(
                        slide,
                        wsi_id=path.stem,
                        level=int(inference["wsi_level"]),
                        macenko_normalizer=normalizer,
                    ),
                    feature_extractor=extractor,
                    model=model,
                    model_id=str(inference["model_id"]),
                    batch_size=BATCH_SIZE,
                    gpu_index=GPU_INDEX,
                )
            )
        finally:
            slide.close()
    # Aggregate WSI ISR from the same concatenated tile predictions that are exported.
    predictions = pd.concat(prediction_frames, ignore_index=True)
    predictions.to_csv(output_dir / "wsi_tile_predictions.csv", index=False)
    aggregate_isr(predictions).to_csv(output_dir / "wsi_isr.csv", index=False)


def run(config: dict[str, Any]) -> None:
    """Execute each explicitly enabled Figure 5 SpMap stage in workflow order.

    Parameters: ``config`` is the editable ``CONFIG`` mapping at module scope.
    Returns: no value; creates outputs only for enabled stages and fails fast on
    missing CUDA, invalid adapters, missing inputs, or retained contract checks.
    """
    # Enforce CUDA before dispatching enabled stages in their scientific workflow order.
    enabled = config["enabled"]
    cuda_stages = ("extract_reference_conch", "train_evaluate", "infer_wsi")
    if any(enabled[stage] for stage in cuda_stages):
        require_cuda(GPU_INDEX)
    if enabled["build_reference_manifest"]:
        build_reference_manifest(config)
    if enabled["extract_reference_conch"]:
        extract_reference_conch(config)
    if enabled["train_evaluate"]:
        train_and_evaluate(config)
    if enabled["infer_wsi"]:
        infer_wsi_and_aggregate_isr(config)


if __name__ == "__main__":
    run(CONFIG)

"""Run CUDA-only CONCH-to-MLP inference for retained SpMap WSI tiles.

Purpose:
    Generate the canonical probability-bearing tile prediction table from WSI
    tiles, paired CONCH representations, and a selected SpMap checkpoint.
Figure 5 callers:
    Figure 5 SpMap WSI-inference workflows use ``load_checkpoint`` and
    ``predict_wsi_tiles`` before WSI-level ISR aggregation.
Inputs:
    Retained ``WSITile`` records, a configured ``ConchFeatureExtractor``, an
    MLP checkpoint/model, and a stable model identifier.
Outputs:
    A DataFrame with tile coordinates, model identity, four probabilities, and
    the argmax class under ``PREDICTION_COLUMNS``.
Ordered use:
    Open and tile the WSI with ``wsi.py``, load the checkpoint, call
    ``predict_wsi_tiles``, and supply the resulting table to ``isr.py``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd
import torch

from .conch_features import ConchFeatureExtractor
from .mlp import PROBABILITY_COLUMNS, SpMapMLP, require_cuda
from .wsi import WSITile


PREDICTION_COLUMNS = (
    "tile_id",
    "wsi_id",
    "level",
    "x",
    "y",
    "x_level0",
    "y_level0",
    "w_real",
    "h_real",
    "pad_right",
    "pad_bottom",
    "model_id",
    "pred_class",
    *PROBABILITY_COLUMNS,
)


def load_checkpoint(
    checkpoint_path: str | Path,
    *,
    gpu_index: int = 0,
) -> SpMapMLP:
    """Load one canonical-architecture checkpoint onto the selected CUDA device.

    Parameters:
        checkpoint_path: Path to a state dictionary emitted by the SpMap trainer.
        gpu_index: Zero-based CUDA device index.

    Returns:
        An evaluation-mode ``SpMapMLP`` with the strict checkpoint state loaded.
    """
    # Recreate the canonical architecture on CUDA before loading its strict state.
    device = require_cuda(gpu_index)
    model = SpMapMLP().to(device)
    state = torch.load(checkpoint_path, map_location=device, weights_only=True)
    model.load_state_dict(state, strict=True)
    # Freeze dropout behavior for all downstream tile predictions.
    model.eval()
    return model


def _predict_batch(
    records: list[WSITile],
    *,
    feature_extractor: ConchFeatureExtractor,
    model: SpMapMLP,
    model_id: str,
    device: torch.device,
) -> list[dict[str, object]]:
    """Predict one complete retained-tile batch and build output row mappings.

    Parameters:
        records: Retained WSI tile records for one inference batch.
        feature_extractor: CUDA CONCH extractor configured for the same run.
        model: Evaluation-mode SpMap MLP.
        model_id: Stable identifier written to every prediction row.
        device: CUDA device receiving concatenated feature tensors.

    Returns:
        Prediction-row dictionaries following ``PREDICTION_COLUMNS`` order.
    """
    # Extract paired CONCH vectors in the same order as the WSI tile records.
    features = feature_extractor.extract(
        [record.tile_id for record in records],
        [record.image for record in records],
    )
    # Concatenate the canonical representations and transfer one batch to CUDA.
    feature_tensor = torch.from_numpy(features.concatenated()).to(
        device, non_blocking=True
    )
    # Convert logits to four-class probabilities without retaining gradients.
    with torch.no_grad():
        probabilities = torch.softmax(model(feature_tensor), dim=1).cpu().numpy()
    # Join predictions back to complete spatial metadata by batch order.
    rows = []
    for record, probability in zip(records, probabilities):
        row = {
            "tile_id": record.tile_id,
            "wsi_id": record.wsi_id,
            "level": record.level,
            "x": record.x,
            "y": record.y,
            "x_level0": record.x_level0,
            "y_level0": record.y_level0,
            "w_real": record.w_real,
            "h_real": record.h_real,
            "pad_right": record.pad_right,
            "pad_bottom": record.pad_bottom,
            "model_id": model_id,
            "pred_class": int(probability.argmax()),
        }
        # Store probabilities in the canonical class-column order.
        for class_id, column in enumerate(PROBABILITY_COLUMNS):
            row[column] = float(probability[class_id])
        rows.append(row)
    return rows


def predict_wsi_tiles(
    tiles: Iterable[WSITile],
    *,
    feature_extractor: ConchFeatureExtractor,
    model: SpMapMLP,
    model_id: str,
    batch_size: int = 256,
    gpu_index: int = 0,
) -> pd.DataFrame:
    """Predict retained WSI tiles from preloaded feature arrays.

    Parameters:
        tiles: WSI tile iterator that includes keep/skip metadata.
        feature_extractor: Configured CUDA CONCH feature extractor.
        model: Selected canonical SpMap MLP.
        model_id: Non-empty model identifier stored in the output table.
        batch_size: Number of retained tiles evaluated per CONCH/MLP batch.
        gpu_index: Zero-based CUDA device index.

    Returns:
        A DataFrame with the fixed ``PREDICTION_COLUMNS`` schema.
    """
    # Validate caller-visible output identity and the retained-tile batch size.
    if not model_id:
        raise ValueError("model_id is required")
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    # Place the selected MLP in CUDA evaluation mode for the full WSI stream.
    device = require_cuda(gpu_index)
    model = model.to(device)
    model.eval()
    # Skip rejected grid positions and collect retained records into fixed batches.
    rows = []
    batch = []
    for tile in tiles:
        if tile.status != "kept":
            continue
        if tile.image is None:
            raise ValueError(f"{tile.tile_id}: retained tile has no image")
        batch.append(tile)
        # Run each complete batch through the shared CONCH-to-MLP path.
        if len(batch) == batch_size:
            rows.extend(
                _predict_batch(
                    batch,
                    feature_extractor=feature_extractor,
                    model=model,
                    model_id=model_id,
                    device=device,
                )
            )
            batch = []
    # Evaluate the final partial batch without dropping retained edge records.
    if batch:
        rows.extend(
            _predict_batch(
                batch,
                feature_extractor=feature_extractor,
                model=model,
                model_id=model_id,
                device=device,
            )
        )
    return pd.DataFrame(rows, columns=PREDICTION_COLUMNS)

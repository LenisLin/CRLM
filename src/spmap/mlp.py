"""Train and evaluate the canonical CUDA-only SpMap MLP.

Purpose:
    Preload paired CONCH tar-shard features once into RAM, train the fixed
    four-class MLP by grouped fold, and produce probability, metric, and
    checkpoint metadata contracts.
Figure 5 callers:
    Figure 5 SpMap training workflows use these helpers for five-fold model
    fitting and Supplementary Figure S11 evaluation summaries.
Inputs:
    Split manifests with tile/class keys, canonical 768- and 512-dimensional
    tar shards, and CUDA training parameters.
Outputs:
    Fitted MLPs, tile-level probability DataFrames, metric mappings/tables,
    seed-variability summaries, and checkpoint manifest mappings.
Ordered use:
    Build splits, construct one ``InMemoryFeatureStore`` from every required
    manifest, train each fold, generate predictions, then summarize outputs.
"""

from __future__ import annotations

import glob
import io
import random
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    classification_report,
    confusion_matrix,
    precision_recall_fscore_support,
    roc_auc_score,
)
from torch import nn
from torch.utils.data import DataLoader, Dataset


CLASS_NAMES = ("TC", "PIR", "PSM", "OTHER_PT")
PROBABILITY_COLUMNS = tuple(f"prob_{name.lower()}" for name in CLASS_NAMES)
INPUT_DIM = 1280
HIDDEN_DIM = 450
OUTPUT_DIM = 4


def require_cuda(gpu_index: int = 0) -> torch.device:
    """Resolve the CUDA device used for MLP training and evaluation.

    Parameters:
        gpu_index: Zero-based index of the CUDA device.

    Returns:
        The selected ``torch.device``.
    """
    # Require an available in-range CUDA device; training has no CPU path.
    if not torch.cuda.is_available():
        raise RuntimeError("SpMap MLP training requires CUDA")
    if gpu_index < 0 or gpu_index >= torch.cuda.device_count():
        raise ValueError(
            f"Requested GPU {gpu_index}; available device count is {torch.cuda.device_count()}"
        )
    return torch.device(f"cuda:{gpu_index}")


def set_seed(seed: int) -> None:
    """Set the retained random seeds for one MLP training run.

    Parameters:
        seed: Integer applied to Python, NumPy, CPU PyTorch, and CUDA PyTorch.

    Returns:
        ``None``. Random-number-generator state is updated in place.
    """
    # Synchronize Python, NumPy, and PyTorch RNG state for one training run.
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def _archive_group(class_id: int) -> str:
    """Map a fixed SpMap class identifier to its tar-archive group.

    Parameters:
        class_id: Canonical class ID in the range 0 through 3.

    Returns:
        The archive directory name containing the class's feature shards.
    """
    # Preserve the canonical class-to-shard layout shared by both feature widths.
    if class_id == 0:
        return "tc"
    if class_id in (1, 2):
        return "pir_psm"
    if class_id == 3:
        return "other_pt"
    raise ValueError(f"Invalid class_id {class_id}")


class InMemoryFeatureStore:
    """Store required canonical tar-shard features in one contiguous RAM matrix.

    Attributes:
        index: Mapping from ``(archive_group, tile_id)`` to matrix row index.
        matrix: Float32 matrix with one 1,280-dimensional feature row per tile.
    """

    def __init__(
        self,
        archive_root: str | Path,
        frames: Iterable[pd.DataFrame],
        *,
        tile_col: str = "tile_id",
        class_col: str = "class_id",
        npz_key: str = "arr_0",
    ) -> None:
        """Load each required paired tar-shard feature once before fold training.

        Parameters:
            archive_root: Root directory containing class/dimension tar shards.
            frames: All manifests whose tile features will be accessed.
            tile_col: Manifest column containing canonical tile IDs.
            class_col: Manifest column containing canonical class IDs.
            npz_key: Key of the feature vector in every NPZ archive member.

        Returns:
            ``None``. The instance retains the tile index and contiguous feature
            matrix for dataset construction.
        """
        # Collect the unique archive-group/tile keys required by all supplied manifests.
        keys = set()
        for frame in frames:
            missing = {tile_col, class_col} - set(frame.columns)
            if missing:
                raise ValueError(f"Feature manifest is missing columns: {sorted(missing)}")
            keys.update(
                (_archive_group(int(class_id)), str(tile_id))
                for tile_id, class_id in frame[[tile_col, class_col]].itertuples(index=False)
            )
        # Allocate one stable row per tile for the complete concatenated feature vector.
        ordered_keys = sorted(keys)
        self.index = {key: index for index, key in enumerate(ordered_keys)}
        self.matrix = np.empty((len(ordered_keys), INPUT_DIM), dtype=np.float32)
        archive_root = Path(archive_root)

        # Scan each canonical shard set once and write vectors into their fixed offsets.
        for dimension, offset in ((768, 0), (512, 768)):
            seen = set()
            for group in ("tc", "pir_psm", "other_pt"):
                shard_pattern = str(archive_root / group / str(dimension) / "part*.tar")
                shard_paths = sorted(glob.glob(shard_pattern))
                if not shard_paths:
                    raise FileNotFoundError(shard_pattern)
                # Stream archive members while retaining only manifest-requested tiles.
                for shard_path in shard_paths:
                    with tarfile.open(shard_path, "r:") as archive:
                        for member in archive:
                            if not member.isfile():
                                continue
                            tile_id = Path(member.name).stem
                            key = (group, tile_id)
                            row_index = self.index.get(key)
                            if row_index is None:
                                continue
                            handle = archive.extractfile(member)
                            # Decode in memory and enforce the feature-width contract.
                            with np.load(io.BytesIO(handle.read())) as npz:
                                vector = np.asarray(npz[npz_key], dtype=np.float32)
                            if vector.shape != (dimension,):
                                raise ValueError(
                                    f"{tile_id}: expected ({dimension},), observed {vector.shape}"
                                )
                            self.matrix[row_index, offset : offset + dimension] = vector
                            seen.add(key)
            # Require every requested tile before moving to training from RAM.
            missing_keys = set(ordered_keys) - seen
            if missing_keys:
                first = sorted(missing_keys)[0]
                raise RuntimeError(
                    f"Missing {dimension}-dimensional features for "
                    f"{len(missing_keys)} tiles; first: {first}"
                )

    def row_indices(
        self,
        frame: pd.DataFrame,
        *,
        tile_col: str = "tile_id",
        class_col: str = "class_id",
    ) -> np.ndarray:
        """Translate ordered manifest rows to preloaded feature-matrix indices.

        Parameters:
            frame: Manifest rows requested by a dataset or prediction operation.
            tile_col: Column containing canonical tile IDs.
            class_col: Column containing canonical class IDs.

        Returns:
            Int64 matrix-row indices preserving the input frame order.
        """
        # Resolve every manifest key against the archive index while preserving row order.
        return np.asarray(
            [
                self.index[(_archive_group(int(class_id)), str(tile_id))]
                for tile_id, class_id in frame[[tile_col, class_col]].itertuples(index=False)
            ],
            dtype=np.int64,
        )


class FeatureDataset(Dataset):
    """Expose manifest-aligned RAM feature rows through the PyTorch dataset API."""

    def __init__(
        self,
        frame: pd.DataFrame,
        feature_store: InMemoryFeatureStore,
        *,
        tile_col: str = "tile_id",
        class_col: str = "class_id",
    ) -> None:
        """Construct a dataset from one fold or evaluation manifest.

        Parameters:
            frame: Tile manifest with feature keys and the target class column.
            feature_store: Preloaded paired CONCH feature store.
            tile_col: Column containing canonical tile IDs.
            class_col: Column containing integer class labels.

        Returns:
            ``None``. A reset-index manifest and aligned feature indices are kept
            for PyTorch data loading.
        """
        # Freeze manifest row order before translating rows to RAM feature indices.
        self.frame = frame.reset_index(drop=True).copy()
        self.class_col = class_col
        self.feature_store = feature_store
        self.feature_indices = feature_store.row_indices(
            self.frame, tile_col=tile_col, class_col=class_col
        )

    def __len__(self) -> int:
        """Return the number of manifest rows available to the data loader.

        Returns:
            The dataset length.
        """
        return len(self.frame)

    def __getitem__(self, index: int):
        """Return one feature vector, class label, and manifest row index.

        Parameters:
            index: Zero-based row index in the reset dataset manifest.

        Returns:
            A tuple containing the feature tensor, integer-label tensor, and
            original dataset-row-index tensor.
        """
        # Package one preloaded feature row with its aligned class and dataset index.
        vector = self.feature_store.matrix[self.feature_indices[index]]
        label = int(self.frame.iloc[index][self.class_col])
        return torch.from_numpy(vector), torch.tensor(label), torch.tensor(index)


class SpMapMLP(nn.Module):
    """Implement the fixed 1280 -> 450 -> ReLU -> dropout -> 4 classifier."""

    def __init__(self, *, dropout: float = 0.25) -> None:
        """Initialize the retained fully connected SpMap architecture.

        Parameters:
            dropout: Dropout probability between the hidden ReLU and output layer.

        Returns:
            ``None``. The sequential model layers are assigned to ``self.net``.
        """
        # Instantiate the fixed hidden width, activation, dropout, and class head.
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(INPUT_DIM, HIDDEN_DIM),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(HIDDEN_DIM, OUTPUT_DIM),
        )

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        """Calculate class logits for a batch of concatenated CONCH features.

        Parameters:
            features: Tensor with final dimension ``INPUT_DIM``.

        Returns:
            Unnormalized four-class logits with one row per input feature vector.
        """
        return self.net(features)


def inverse_frequency_weights(
    frame: pd.DataFrame,
    *,
    class_col: str = "class_id",
    device: torch.device,
) -> torch.Tensor:
    """Calculate mean-normalized inverse-frequency training-class weights.

    Parameters:
        frame: Training manifest containing the complete four-class target set.
        class_col: Column containing canonical integer class IDs.
        device: CUDA device receiving the returned tensor.

    Returns:
        Four float32 class weights ordered by canonical class ID.
    """
    # Count canonical training classes and require every output class to be represented.
    counts = frame[class_col].astype(int).value_counts()
    missing = [class_id for class_id in range(4) if class_id not in counts]
    if missing:
        raise ValueError(f"Training data omit classes: {missing}")
    # Order inverse frequencies by class ID and normalize their mean to one.
    frequencies = np.asarray([counts[class_id] for class_id in range(4)], dtype=float)
    weights = (len(frame) / frequencies).astype(np.float32)
    weights /= weights.mean()
    return torch.tensor(weights, device=device)


def _loader(
    dataset: FeatureDataset,
    *,
    batch_size: int,
    shuffle: bool,
    num_workers: int = 6,
) -> DataLoader:
    """Build the fixed pinned-memory data loader for an in-memory dataset.

    Parameters:
        dataset: Manifest-aligned RAM feature dataset.
        batch_size: Number of rows yielded in each batch.
        shuffle: Whether to shuffle row order for this loader.
        num_workers: Number of worker processes used to load RAM-resident rows.

    Returns:
        A ``DataLoader`` with the retained worker and memory settings.
    """
    # Use pinned batches and configured workers for RAM-resident feature transfer.
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        pin_memory=True,
    )


def _predict_loader(
    model: SpMapMLP,
    loader: DataLoader,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate a loader and collect labels, probabilities, and row indices.

    Parameters:
        model: SpMap MLP evaluated in inference mode.
        loader: Data loader yielding features, labels, and dataset row indices.
        device: CUDA device used for feature tensors and model execution.

    Returns:
        A tuple of labels, four-class probabilities, and dataset row indices.
    """
    # Disable dropout and gradient tracking during loader evaluation.
    model.eval()
    labels, probabilities, indices = [], [], []
    with torch.no_grad():
        # Move only features to CUDA and retain labels plus row indices on the host.
        for features, targets, row_indices in loader:
            logits = model(features.to(device, non_blocking=True))
            labels.append(targets.numpy())
            probabilities.append(torch.softmax(logits, dim=1).cpu().numpy())
            indices.append(row_indices.numpy())
    # Preserve fixed output ranks when the loader contains no rows.
    if not labels:
        return (
            np.empty(0, dtype=int),
            np.empty((0, OUTPUT_DIM), dtype=float),
            np.empty(0, dtype=int),
        )
    # Concatenate batches without changing their loader-defined row order.
    return np.concatenate(labels), np.concatenate(probabilities), np.concatenate(indices)


def prediction_frame(
    model: SpMapMLP,
    dataset: FeatureDataset,
    *,
    batch_size: int = 256,
    gpu_index: int = 0,
    num_workers: int = 6,
) -> pd.DataFrame:
    """Generate one probability-bearing prediction row for every dataset tile.

    Parameters:
        model: Trained SpMap MLP to evaluate.
        dataset: Dataset supplying metadata, features, and class labels.
        batch_size: Number of dataset rows evaluated per MLP batch.
        gpu_index: Zero-based CUDA device index.
        num_workers: Number of worker processes used by the evaluation loader.

    Returns:
        The dataset manifest enriched with ``pred_class`` and all named
        probability columns.
    """
    # Evaluate the model on CUDA with a nonshuffled manifest-aligned loader.
    device = require_cuda(gpu_index)
    model = model.to(device)
    _, probabilities, indices = _predict_loader(
        model,
        _loader(
            dataset,
            batch_size=batch_size,
            shuffle=False,
            num_workers=num_workers,
        ),
        device,
    )
    # Restore manifest metadata by returned row index before adding predictions.
    output = dataset.frame.iloc[indices].reset_index(drop=True)
    output["pred_class"] = probabilities.argmax(axis=1)
    # Store probability columns in canonical class-ID order.
    for class_id, column in enumerate(PROBABILITY_COLUMNS):
        output[column] = probabilities[:, class_id]
    return output


@dataclass(frozen=True)
class TrainingResult:
    """Describe the selected state from one fold-training run.

    Attributes:
        model: MLP restored to the best validation-accuracy state.
        selected_epoch: One-based epoch at which the selected state was observed.
        validation_accuracy: Accuracy of the selected validation state.
        history: Per-epoch training loss and validation-accuracy DataFrame.
    """

    model: SpMapMLP
    selected_epoch: int
    validation_accuracy: float
    history: pd.DataFrame


def train_fold(
    train_frame: pd.DataFrame,
    validation_frame: pd.DataFrame,
    feature_store: InMemoryFeatureStore,
    *,
    epochs: int = 100,
    batch_size: int = 256,
    learning_rate: float = 1e-3,
    weight_decay: float = 1e-4,
    dropout: float = 0.25,
    seed: int = 24,
    gpu_index: int = 0,
    num_workers: int = 6,
    reset_seed: bool = True,
) -> TrainingResult:
    """Train one fold and select its state by validation accuracy.

    Parameters:
        train_frame: Training-manifest rows for the current fold.
        validation_frame: Validation-manifest rows for the current fold.
        feature_store: One preconstructed RAM feature store shared by all folds.
        epochs: Number of complete training epochs.
        batch_size: Number of tile rows in each training/evaluation batch.
        learning_rate: AdamW learning rate.
        weight_decay: AdamW weight-decay coefficient.
        dropout: Hidden-layer dropout probability for the MLP.
        seed: Integer training seed.
        gpu_index: Zero-based CUDA device index.
        num_workers: Number of worker processes used by each data loader.
        reset_seed: Whether to reset RNG state before this fold. The default
            preserves standalone library behavior; ordered multi-fold callers
            can seed once before their fold sequence and pass ``False``.

    Returns:
        A ``TrainingResult`` containing the selected model, epoch, validation
        accuracy, and epoch history.
    """
    # Validate the epoch count and resolve the mandatory CUDA device.
    if epochs < 1:
        raise ValueError("epochs must be positive")
    device = require_cuda(gpu_index)
    # Reset RNGs for standalone runs or retain caller-managed fold sequencing.
    if reset_seed:
        set_seed(seed)
    # Bind both manifests to the shared in-memory feature matrix.
    train_dataset = FeatureDataset(train_frame, feature_store)
    validation_dataset = FeatureDataset(validation_frame, feature_store)
    # Shuffle only training rows; keep validation order stable across epochs.
    train_loader = _loader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
    )
    validation_loader = _loader(
        validation_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
    )
    # Configure the fixed MLP, class-balanced loss, and AdamW optimizer on CUDA.
    model = SpMapMLP(dropout=dropout).to(device)
    criterion = nn.CrossEntropyLoss(
        weight=inverse_frequency_weights(train_frame, device=device)
    )
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=learning_rate, weight_decay=weight_decay
    )
    # Track the single validation-accuracy-selected state and per-epoch history.
    best_accuracy = -1.0
    selected_epoch = 0
    best_state = None
    history = []

    # Train for the prespecified number of complete epochs.
    for epoch in range(1, epochs + 1):
        model.train()
        loss_total = 0.0
        sample_total = 0
        # Transfer each RAM batch nonblockingly and update model parameters once.
        for features, targets, _ in train_loader:
            features = features.to(device, non_blocking=True)
            targets = targets.to(device, non_blocking=True)
            optimizer.zero_grad()
            loss = criterion(model(features), targets)
            if not torch.isfinite(loss):
                raise RuntimeError(f"Epoch {epoch}: non-finite loss")
            loss.backward()
            optimizer.step()
            loss_total += float(loss.item()) * len(features)
            sample_total += len(features)

        # Evaluate the complete validation fold after the epoch update phase.
        labels, probabilities, _ = _predict_loader(model, validation_loader, device)
        predicted = probabilities.argmax(axis=1)
        accuracy = float(accuracy_score(labels, predicted))
        # Snapshot the first strictly improved validation-accuracy state on the host.
        if accuracy > best_accuracy:
            best_accuracy = accuracy
            selected_epoch = epoch
            best_state = {
                name: value.detach().cpu().clone()
                for name, value in model.state_dict().items()
            }
        # Record sample-weighted training loss and the epoch selection indicator.
        history.append(
            {
                "epoch": epoch,
                "train_loss": loss_total / sample_total,
                "validation_accuracy": accuracy,
                "selected": accuracy > max(
                    (row["validation_accuracy"] for row in history), default=-1.0
                ),
            }
        )

    # Restore the selected state before returning the fold result.
    model.load_state_dict(best_state)
    return TrainingResult(
        model=model,
        selected_epoch=selected_epoch,
        validation_accuracy=best_accuracy,
        history=pd.DataFrame(history),
    )


def evaluate_predictions(
    predictions: pd.DataFrame,
    *,
    true_col: str = "class_id",
) -> dict[str, object]:
    """Calculate tile-level macro metrics and the fixed four-class matrix.

    Parameters:
        predictions: Probability-bearing tile prediction table.
        true_col: Column containing integer ground-truth class IDs.

    Returns:
        Accuracy, macro precision/recall/F1, and a fixed-order confusion matrix.
    """
    # Require the complete prediction schema used by the four-class metric contract.
    required = {true_col, "pred_class", *PROBABILITY_COLUMNS}
    missing = required - set(predictions.columns)
    if missing:
        raise ValueError(f"Prediction table is missing columns: {sorted(missing)}")
    # Calculate fixed-class macro metrics independently of observed class order.
    y_true = predictions[true_col].astype(int)
    y_pred = predictions["pred_class"].astype(int)
    precision, recall, f1, _ = precision_recall_fscore_support(
        y_true, y_pred, labels=range(4), average="macro", zero_division=0
    )
    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "precision_macro": float(precision),
        "recall_macro": float(recall),
        "f1_macro": float(f1),
        "confusion_matrix": confusion_matrix(y_true, y_pred, labels=range(4)),
    }


def validate_prediction_table(
    predictions: pd.DataFrame,
    *,
    tile_col: str = "tile_id",
    true_col: str = "class_id",
) -> None:
    """Validate the probability-bearing four-class tile output contract.

    Parameters:
        predictions: Prediction table to validate.
        tile_col: Column containing unique tile identifiers.
        true_col: Column containing integer ground-truth class IDs.

    Returns:
        ``None`` when the required schema, class values, and probability
        relationships satisfy the contract.
    """
    # Require identifiers, labels, argmax predictions, and all probability columns.
    required = {tile_col, true_col, "pred_class", *PROBABILITY_COLUMNS}
    missing = required - set(predictions.columns)
    if missing:
        raise ValueError(f"Prediction table is missing columns: {sorted(missing)}")
    # Enforce unique tile rows and canonical true/predicted class identifiers.
    if predictions[tile_col].astype(str).duplicated().any():
        raise ValueError(f"{tile_col} values must be unique")
    for column in (true_col, "pred_class"):
        invalid = sorted(set(predictions[column].astype(int)) - set(range(4)))
        if invalid:
            raise ValueError(f"Invalid {column} values: {invalid}")
    # Check finite normalized probabilities and their agreement with pred_class.
    probabilities = predictions.loc[:, PROBABILITY_COLUMNS].to_numpy(dtype=float)
    if not np.isfinite(probabilities).all():
        raise ValueError("Prediction probabilities must be finite")
    if not np.allclose(probabilities.sum(axis=1), 1.0, atol=1e-6, rtol=0):
        raise ValueError("Prediction probability rows must sum to one")
    if not np.array_equal(
        probabilities.argmax(axis=1), predictions["pred_class"].to_numpy(dtype=int)
    ):
        raise ValueError("pred_class must equal the probability argmax")


def summarize_oof_predictions(
    fold_predictions: Iterable[pd.DataFrame],
    *,
    tile_col: str = "tile_id",
    true_col: str = "class_id",
) -> dict[str, object]:
    """Pool fold-validation rows and calculate S11B/C summary artifacts.

    Parameters:
        fold_predictions: One validated prediction DataFrame per fold.
        tile_col: Column containing unique tile identifiers.
        true_col: Column containing integer ground-truth class IDs.

    Returns:
        A mapping with pooled predictions, summary metrics, per-class metrics,
        confusion counts, and row-normalized confusion proportions.
    """
    # Validate each fold independently and attach its one-based fold identity.
    validated = []
    for fold, frame in enumerate(fold_predictions, start=1):
        validate_prediction_table(frame, tile_col=tile_col, true_col=true_col)
        validated.append(frame.assign(fold=fold))
    # Pool each validation tile exactly once across the OOF prediction table.
    pooled = pd.concat(validated, ignore_index=True)
    if pooled[tile_col].astype(str).duplicated().any():
        raise ValueError("Pooled OOF predictions contain repeated tile IDs")

    # Materialize aligned labels and canonical-order probabilities for all summaries.
    y_true = pooled[true_col].to_numpy(dtype=int)
    y_pred = pooled["pred_class"].to_numpy(dtype=int)
    probabilities = pooled.loc[:, PROBABILITY_COLUMNS].to_numpy(dtype=float)
    # Build the fixed four-class precision, recall, F1, and support table.
    report = pd.DataFrame(
        classification_report(
            y_true,
            y_pred,
            labels=range(4),
            target_names=CLASS_NAMES,
            output_dict=True,
            zero_division=0,
        )
    ).transpose()
    # Add one-versus-rest ranking metrics for each canonical class.
    per_class_rows = []
    for class_id, class_name in enumerate(CLASS_NAMES):
        binary = (y_true == class_id).astype(int)
        per_class_rows.append(
            {
                "class_id": class_id,
                "class_name": class_name,
                "support": int(binary.sum()),
                "precision": float(report.loc[class_name, "precision"]),
                "recall": float(report.loc[class_name, "recall"]),
                "f1": float(report.loc[class_name, "f1-score"]),
                "auroc_ovr": float(
                    roc_auc_score(binary, probabilities[:, class_id])
                ),
                "auprc_ovr": float(
                    average_precision_score(binary, probabilities[:, class_id])
                ),
            }
        )
    per_class = pd.DataFrame(per_class_rows)
    # Preserve canonical axes for count and row-normalized confusion matrices.
    counts = confusion_matrix(y_true, y_pred, labels=range(4))
    row_normalized = counts / counts.sum(axis=1, keepdims=True)
    # Aggregate global, macro, and PIR/PSM-specific scalar metrics.
    metrics = {
        "rows": len(pooled),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "precision_macro": float(report.loc["macro avg", "precision"]),
        "recall_macro": float(report.loc["macro avg", "recall"]),
        "f1_macro": float(report.loc["macro avg", "f1-score"]),
        "auroc_macro_ovr": float(per_class["auroc_ovr"].mean()),
        "auprc_macro_ovr": float(per_class["auprc_ovr"].mean()),
        "auroc_pir": float(per_class.loc[1, "auroc_ovr"]),
        "auroc_psm": float(per_class.loc[2, "auroc_ovr"]),
        "auprc_pir": float(per_class.loc[1, "auprc_ovr"]),
        "auprc_psm": float(per_class.loc[2, "auprc_ovr"]),
    }
    metrics["auprc_pir_psm_mean"] = (
        metrics["auprc_pir"] + metrics["auprc_psm"]
    ) / 2
    return {
        "pooled_predictions": pooled,
        "metrics": metrics,
        "per_class_metrics": per_class,
        "confusion_counts": pd.DataFrame(
            counts, index=CLASS_NAMES, columns=CLASS_NAMES
        ),
        "confusion_row_normalized": pd.DataFrame(
            row_normalized, index=CLASS_NAMES, columns=CLASS_NAMES
        ),
    }


def summarize_training_seed_variability(
    seed_metrics: pd.DataFrame,
    *,
    metric_columns: tuple[str, ...] = (
        "accuracy",
        "f1_macro",
        "auroc_pir",
        "auroc_psm",
        "auprc_pir",
        "auprc_psm",
        "auprc_pir_psm_mean",
    ),
) -> pd.DataFrame:
    """Describe prespecified training-seed variability with t intervals.

    Parameters:
        seed_metrics: One row per training seed with requested metric columns.
        metric_columns: Metrics summarized across the seed rows.

    Returns:
        A table of means, standard deviations, t intervals, minima, and maxima.
    """
    # Import the t distribution only when seed-variability intervals are requested.
    try:
        from scipy.stats import t
    except ImportError as error:
        raise ImportError("Training-seed intervals require scipy") from error
    # Require one unique row per seed and every requested metric column.
    missing = {"seed", *metric_columns} - set(seed_metrics.columns)
    if missing:
        raise ValueError(f"Seed metric table is missing columns: {sorted(missing)}")
    if seed_metrics["seed"].duplicated().any() or len(seed_metrics) < 2:
        raise ValueError("Seed metrics require at least two unique seeds")
    # Use one shared two-sided 95% critical value across metric summaries.
    critical_value = float(t.ppf(0.975, df=len(seed_metrics) - 1))
    rows = []
    # Summarize each metric over the same prespecified seed set.
    for metric in metric_columns:
        values = seed_metrics[metric].astype(float)
        mean = float(values.mean())
        standard_deviation = float(values.std(ddof=1))
        margin = critical_value * standard_deviation / np.sqrt(len(values))
        rows.append(
            {
                "metric": metric,
                "seeds": len(values),
                "mean": mean,
                "standard_deviation": standard_deviation,
                "training_seed_95ci_lower": mean - margin,
                "training_seed_95ci_upper": mean + margin,
                "minimum": float(values.min()),
                "maximum": float(values.max()),
            }
        )
    return pd.DataFrame(rows)


def checkpoint_manifest(
    *,
    model_id: str,
    selected_epoch: int,
    fold: int,
    seed: int,
) -> dict[str, object]:
    """Build metadata describing one selected checkpoint.

    Parameters:
        model_id: Stable caller-defined checkpoint identifier.
        selected_epoch: One-based epoch selected by validation accuracy.
        fold: One-based cross-validation fold number.
        seed: Integer training seed for the fold.

    Returns:
        A mapping with model identity, architecture, selection, and class-map
        metadata.
    """
    # Record the architecture and selection metadata needed to identify the checkpoint.
    return {
        "model_id": model_id,
        "fold": fold,
        "seed": seed,
        "selected_epoch": selected_epoch,
        "architecture": [INPUT_DIM, HIDDEN_DIM, OUTPUT_DIM],
        "activation": "ReLU",
        "dropout": 0.25,
        "selection_metric": "validation_accuracy",
        "class_map": {str(index): name for index, name in enumerate(CLASS_NAMES)},
    }

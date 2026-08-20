"""Extract paired CONCH v1 representations for SpMap tiles on CUDA.

Purpose:
    Convert caller-supplied RGB tiles into the fixed 768-dimensional semantic
    and 512-dimensional alignment representations used by the SpMap MLP.
Figure 5 callers:
    Figure 5 SpMap training and WSI-inference workflows use this shared module
    through ``ConchFeatureExtractor``.
Inputs:
    A CONCH checkpoint identifier, model-loader and representation-adapter
    callables, an RGB preprocessing callable, and ordered tile IDs and images.
Outputs:
    ``FeatureBatch`` instances retaining tile order and the paired feature
    arrays.
Ordered use:
    Construct ``ConchFeatureExtractor`` on the selected CUDA device, call
    ``extract``, then pass ``FeatureBatch.concatenated()`` to the MLP.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np
import torch


SEMANTIC_DIM = 768
ALIGNMENT_DIM = 512
CONCATENATED_DIM = SEMANTIC_DIM + ALIGNMENT_DIM


def require_cuda(gpu_index: int = 0) -> torch.device:
    """Resolve the CUDA device used for CONCH feature extraction.

    Parameters:
        gpu_index: Zero-based index of the CUDA device.

    Returns:
        The selected ``torch.device``.
    """
    # Require an available in-range CUDA device; extraction has no CPU path.
    if not torch.cuda.is_available():
        raise RuntimeError("SpMap CONCH extraction requires CUDA")
    if gpu_index < 0 or gpu_index >= torch.cuda.device_count():
        raise ValueError(
            f"Requested GPU {gpu_index}; available device count is {torch.cuda.device_count()}"
        )
    return torch.device(f"cuda:{gpu_index}")


@dataclass(frozen=True)
class FeatureBatch:
    """Store paired CONCH representations aligned to input tile IDs.

    Attributes:
        tile_ids: Ordered identifiers for the extracted tiles.
        semantic_768: Semantic representation array with one 768-vector per tile.
        alignment_512: Alignment representation array with one 512-vector per tile.
    """

    tile_ids: tuple[str, ...]
    semantic_768: np.ndarray
    alignment_512: np.ndarray

    def concatenated(self) -> np.ndarray:
        """Concatenate semantic and alignment features in canonical order.

        Returns:
            A row-aligned float array with 1,280 columns: semantic features
            followed by alignment features.
        """
        return np.concatenate((self.semantic_768, self.alignment_512), axis=1)


class ConchFeatureExtractor:
    """Hold the CUDA-resident CONCH model and its caller-provided adapters.

    The extractor serves Figure 5 tile-feature preparation and WSI inference.
    """

    def __init__(
        self,
        *,
        model_loader: Callable[[str, torch.device], torch.nn.Module],
        checkpoint: str,
        preprocess: Callable[[np.ndarray], torch.Tensor],
        forward_representations: Callable[
            [torch.nn.Module, torch.Tensor], tuple[torch.Tensor, torch.Tensor]
        ],
        gpu_index: int = 0,
    ) -> None:
        """Load the supplied CONCH model and place it in evaluation mode.

        Parameters:
            model_loader: Builds a CONCH ``torch.nn.Module`` from the checkpoint.
            checkpoint: Checkpoint identifier accepted by ``model_loader``.
            preprocess: Converts one RGB tile into a model input tensor.
            forward_representations: Returns semantic and alignment tensors for a batch.
            gpu_index: Zero-based CUDA device index.

        Returns:
            ``None``. The configured model, device, and adapters are retained on
            the instance.
        """
        # Resolve CUDA first, then retain the caller-provided preprocessing adapters.
        self.device = require_cuda(gpu_index)
        self.preprocess = preprocess
        self.forward_representations = forward_representations
        # Keep the CONCH model CUDA-resident with training-only layer behavior disabled.
        self.model = model_loader(checkpoint, self.device).to(self.device)
        self.model.eval()

    def extract(
        self,
        tile_ids: Iterable[str],
        images: Iterable[np.ndarray],
    ) -> FeatureBatch:
        """Extract paired CONCH representations for one ordered tile batch.

        Parameters:
            tile_ids: Unique tile identifiers in output-row order.
            images: RGB arrays aligned one-to-one with ``tile_ids``.

        Returns:
            A ``FeatureBatch`` whose arrays preserve the supplied tile order.
        """
        # Materialize aligned inputs once so order, length, and identity are fixed.
        ids = tuple(str(tile_id) for tile_id in tile_ids)
        image_list = list(images)
        if len(ids) != len(image_list):
            raise ValueError("tile_ids and images must have equal length")
        if len(set(ids)) != len(ids):
            raise ValueError("tile_ids must be unique within a batch")
        # Preserve both feature-width contracts for an empty input batch.
        if not ids:
            return FeatureBatch(
                ids,
                np.empty((0, SEMANTIC_DIM), dtype=np.float32),
                np.empty((0, ALIGNMENT_DIM), dtype=np.float32),
            )

        # Preprocess in tile order and transfer one stacked batch to CUDA.
        tensors = [self.preprocess(np.asarray(image)) for image in image_list]
        batch = torch.stack(tensors).to(self.device, non_blocking=True)
        # Extract paired representations without building an autograd graph.
        with torch.no_grad():
            semantic, alignment = self.forward_representations(self.model, batch)
        semantic_array = semantic.detach().float().cpu().numpy()
        alignment_array = alignment.detach().float().cpu().numpy()
        # Require one canonical-width vector per supplied tile in both outputs.
        expected_semantic = (len(ids), SEMANTIC_DIM)
        expected_alignment = (len(ids), ALIGNMENT_DIM)
        if semantic_array.shape != expected_semantic:
            raise ValueError(
                f"Expected semantic feature shape {expected_semantic}, "
                f"observed {semantic_array.shape}"
            )
        if alignment_array.shape != expected_alignment:
            raise ValueError(
                f"Expected alignment feature shape {expected_alignment}, "
                f"observed {alignment_array.shape}"
            )
        return FeatureBatch(ids, semantic_array, alignment_array)

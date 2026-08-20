"""Preprocess SpMap training tiles and generate named augmentations.

Purpose:
    Apply the retained tile-quality, luminosity-standardization, Vahadane, and
    offline-augmentation procedures used to prepare SpMap training images.
Figure 5 callers:
    Figure 5 SpMap training-tile preparation workflows use these shared
    preprocessing functions before CONCH feature extraction.
Inputs:
    RGB uint8 tile arrays, the specified Vahadane normalizer, and a named
    augmentation policy where augmentation is required.
Outputs:
    Processed tile records or augmented tile records preserving tile lineage.
Ordered use:
    Preprocess retained reference tiles, apply the selected augmentation policy,
    then submit resulting RGB arrays to the CONCH feature extractor.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


VAHADANE_TARGET_MAX_CONCENTRATIONS = np.array(
    [[1.92984228, 1.13784735]], dtype=np.float64
)
VAHADANE_TARGET_STAIN_MATRIX = np.array(
    [
        [0.48011481, 0.75784925, 0.44176271],
        [0.20473679, 0.80759910, 0.55305202],
    ],
    dtype=np.float64,
)
AUGMENTATION_POLICIES = {
    "pir_psm_geometry",
    "pir_psm_stochastic",
    "other_pt_geometry",
}


@dataclass(frozen=True)
class ProcessedTile:
    """Record one training tile's quality-filter and normalization result.

    Attributes:
        image: Normalized RGB image for a retained tile, otherwise ``None``.
        kept: Whether the tile satisfied the white-area criterion.
        skip_reason: Stable reason assigned when the tile is skipped.
    """

    image: np.ndarray | None
    kept: bool
    skip_reason: str | None


@dataclass(frozen=True)
class AugmentedTile:
    """Record one offline augmentation and its source-tile lineage.

    Attributes:
        tile_id: Identifier assigned to the augmented image.
        parent_tile_id: Identifier of the unaugmented source tile.
        policy: Name of the augmentation policy that generated the variant.
        image: Contiguous RGB uint8 variant array.
    """

    tile_id: str
    parent_tile_id: str
    policy: str
    image: np.ndarray


def _as_rgb_uint8(image: np.ndarray) -> np.ndarray:
    """Validate and make contiguous one RGB uint8 tile array.

    Parameters:
        image: H x W x 3 image array.

    Returns:
        A contiguous H x W x 3 ``uint8`` RGB array.
    """
    # Normalize array-like input without changing the required RGB byte contract.
    array = np.asarray(image)
    if array.ndim != 3 or array.shape[2] != 3 or array.dtype != np.uint8:
        raise ValueError("image must be an H x W x 3 uint8 RGB array")
    return np.ascontiguousarray(array)


def white_fraction(image: np.ndarray, *, white_threshold: int = 230) -> float:
    """Calculate the fraction of pixels white in all three RGB channels.

    Parameters:
        image: H x W x 3 uint8 RGB tile array.
        white_threshold: Per-channel inclusive RGB threshold for white pixels.

    Returns:
        The tile fraction whose three RGB channels meet the threshold.
    """
    rgb = _as_rgb_uint8(image)
    return float(np.all(rgb >= white_threshold, axis=2).mean())


def luminosity_standardize(
    image: np.ndarray,
    *,
    percentile: float = 95.0,
) -> np.ndarray:
    """Apply retained LAB-luminosity standardization to one RGB tile.

    Parameters:
        image: H x W x 3 uint8 RGB tile array.
        percentile: LAB luminosity percentile used as the scaling reference.

    Returns:
        Luminosity-standardized H x W x 3 uint8 RGB tile array.
    """
    # Import OpenCV only for callers that request luminosity standardization.
    try:
        import cv2
    except ImportError as error:
        raise ImportError("luminosity_standardize requires opencv-python") from error

    # Scale the LAB luminosity channel to the retained percentile reference.
    rgb = _as_rgb_uint8(image)
    lab = cv2.cvtColor(rgb, cv2.COLOR_RGB2LAB)
    luminosity = lab[:, :, 0].astype(float)
    scale = float(np.percentile(luminosity, percentile))
    if scale <= 0:
        raise ValueError("luminosity percentile must be positive")
    lab[:, :, 0] = np.clip(255 * luminosity / scale, 0, 255).astype(np.uint8)
    return cv2.cvtColor(lab, cv2.COLOR_LAB2RGB)


class VahadaneTrainingNormalizer:
    """Apply the retained training-side Vahadane transform with fixed targets."""

    def __init__(self) -> None:
        """Load the Vahadane extraction functions from staintools.

        Returns:
            ``None``. The stain-matrix and concentration functions are retained
            for calls to this normalizer.
        """
        # Resolve the two staintools operations used by every normalization call.
        try:
            from staintools.miscellaneous.get_concentrations import (
                get_concentrations,
            )
            from staintools.stain_extraction.vahadane_stain_extractor import (
                VahadaneStainExtractor,
            )
        except ImportError as error:
            raise ImportError(
                "VahadaneTrainingNormalizer requires the staintools stain-matrix API"
            ) from error
        self._extractor = VahadaneStainExtractor
        self._get_concentrations = get_concentrations

    def __call__(self, image: np.ndarray) -> np.ndarray:
        """Normalize one luminosity-standardized training tile by Vahadane.

        Parameters:
            image: H x W x 3 uint8 RGB tile array.

        Returns:
            H x W x 3 uint8 RGB tile transformed to the fixed target stain basis.
        """
        # Decompose the standardized source tile into its stain concentrations.
        standardized = luminosity_standardize(image)
        stain_matrix = self._extractor.get_stain_matrix(standardized)
        concentrations = self._get_concentrations(standardized, stain_matrix)
        # Match source concentration scales to the fixed training target basis.
        source_max = np.percentile(concentrations, 99, axis=0).reshape((1, 2))
        if np.any(source_max <= 0):
            raise ValueError("Vahadane source concentrations must be positive")
        concentrations *= VAHADANE_TARGET_MAX_CONCENTRATIONS / source_max
        # Reconstruct RGB pixels from the rescaled concentrations and target stains.
        transformed = 255 * np.exp(
            -np.dot(concentrations, VAHADANE_TARGET_STAIN_MATRIX)
        )
        return transformed.reshape(standardized.shape).astype(np.uint8)


def preprocess_training_tile(
    image: np.ndarray,
    *,
    normalizer: VahadaneTrainingNormalizer,
    white_threshold: int = 230,
    maximum_white_fraction: float = 0.80,
) -> ProcessedTile:
    """Apply the retained quality filter and normalization to one training tile.

    Parameters:
        image: H x W x 3 uint8 RGB tile array.
        normalizer: Initialized Vahadane normalizer for retained tiles.
        white_threshold: Per-channel threshold defining a white pixel.
        maximum_white_fraction: Largest retained fraction of white pixels.

    Returns:
        A ``ProcessedTile`` with either the normalized image or skip metadata.
    """
    # Apply the white-area gate to a validated RGB tile before stain processing.
    rgb = _as_rgb_uint8(image)
    if white_fraction(rgb, white_threshold=white_threshold) > maximum_white_fraction:
        return ProcessedTile(None, False, "more_than_80_percent_white")
    # Normalize only retained tiles and preserve their explicit kept state.
    return ProcessedTile(normalizer(rgb), True, None)


def _geometry_variants(image: np.ndarray, policy: str) -> list[tuple[int, np.ndarray]]:
    """Create deterministic geometric variants for one retained named policy.

    Parameters:
        image: Validated RGB source image.
        policy: ``pir_psm_geometry`` or ``other_pt_geometry``.

    Returns:
        Ordered ``(suffix, image)`` geometric variants, including the source image.
    """
    # Both geometry policies retain the source image and horizontal reflection.
    variants = [(0, image), (1, np.fliplr(image))]
    # Add the policy-specific rotations and vertical reflection in stable suffix order.
    if policy == "pir_psm_geometry":
        variants.extend(
            [(2, np.flipud(image)), (3, np.rot90(image, -1)), (4, np.rot90(image, 1))]
        )
    else:
        variants.append((3, np.rot90(image, -1)))
    return variants


def _stochastic_variants(
    image: np.ndarray,
    *,
    seed: int,
) -> list[tuple[int, np.ndarray]]:
    """Create deterministic-seed stochastic PIR/PSM variants.

    Parameters:
        image: Validated RGB source image.
        seed: Base seed incremented for each generated variant.

    Returns:
        Three ordered ``(suffix, image)`` stochastic variants.
    """
    # Import augmentation dependencies only for the stochastic policy.
    try:
        import albumentations as A
        from skimage.color import hed2rgb, rgb2hed
    except ImportError as error:
        raise ImportError(
            "pir_psm_stochastic requires albumentations and scikit-image"
        ) from error

    class HEDPerturb(A.ImageOnlyTransform):
        """Apply bounded random HED stain-space scaling to one RGB image."""

        def apply(self, img, **params):
            """Transform one image through randomized HED channel scaling.

            Parameters:
                img: RGB uint8 image accepted by Albumentations.
                **params: Additional Albumentations transform parameters.

            Returns:
                RGB uint8 image after HED-space perturbation.
            """
            # Perturb hematoxylin and eosin channels independently within fixed bounds.
            rgb = img.astype(np.float32) / 255.0
            hed = rgb2hed(rgb)
            hed[..., 0] *= np.random.uniform(0.9, 1.1)
            hed[..., 1] *= np.random.uniform(0.9, 1.1)
            return (np.clip(hed2rgb(hed), 0.0, 1.0) * 255.0).astype(np.uint8)

    # Compose stain, color, and low-probability blur/noise perturbations.
    augmenter = A.Compose(
        [
            HEDPerturb(p=0.5),
            A.ColorJitter(
                brightness=0.05,
                contrast=0.05,
                saturation=0.05,
                hue=0.05,
                p=0.3,
            ),
            A.OneOf(
                [
                    A.GaussianBlur(blur_limit=(3, 3), sigma_limit=(0.0, 1.0), p=1.0),
                    A.GaussNoise(var_limit=(1.0, 10.0), mean=0, per_channel=True, p=1.0),
                ],
                p=0.1,
            ),
        ]
    )
    # Reset the NumPy seed used by the retained augmentation sequence for each suffix.
    variants = []
    for suffix in range(1, 4):
        np.random.seed(seed + suffix)
        variants.append((suffix, augmenter(image=image)["image"]))
    return variants


def augment_training_tile(
    image: np.ndarray,
    *,
    tile_id: str,
    policy: str,
    seed: int = 24,
) -> list[AugmentedTile]:
    """Apply one retained named offline augmentation policy.

    Parameters:
        image: H x W x 3 uint8 RGB source tile.
        tile_id: Identifier of the unaugmented source tile.
        policy: One member of ``AUGMENTATION_POLICIES``.
        seed: Base random seed for the stochastic policy.

    Returns:
        Augmented tile records retaining source ID, policy, suffix, and image.
    """
    # Normalize input and require an explicitly retained augmentation policy.
    rgb = _as_rgb_uint8(image)
    if policy not in AUGMENTATION_POLICIES:
        raise ValueError(f"Unknown augmentation policy: {policy}")
    # Dispatch to the deterministic geometry or seeded stochastic transform set.
    variants = (
        _stochastic_variants(rgb, seed=seed)
        if policy == "pir_psm_stochastic"
        else _geometry_variants(rgb, policy)
    )
    # Preserve parent lineage while assigning stable suffix-based variant IDs.
    return [
        AugmentedTile(
            tile_id=f"{tile_id}_{suffix}",
            parent_tile_id=tile_id,
            policy=policy,
            image=np.ascontiguousarray(variant),
        )
        for suffix, variant in variants
    ]

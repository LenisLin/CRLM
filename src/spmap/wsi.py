"""Iterate SpMap WSI tiles on a complete non-overlapping grid.

Purpose:
    Read one OpenSlide level, classify each grid location by white-area status,
    normalize retained observed pixels, and right/bottom-pad border tiles.
Figure 5 callers:
    Figure 5 SpMap WSI-inference workflows use ``open_wsi`` and
    ``iter_wsi_tiles`` before CUDA CONCH-to-MLP inference.
Inputs:
    An OpenSlide WSI, a caller-selected level, a WSI identifier, and a Macenko
    normalization callable.
Outputs:
    Ordered ``WSITile`` records for every complete-grid position, including
    spatial metadata and kept/skip status.
Ordered use:
    Open the slide, iterate tiles at the selected level, send retained records
    to ``inference.predict_wsi_tiles``, then aggregate ISR.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterator

import numpy as np

from .preprocessing import white_fraction


TILE_SIZE = 256


@dataclass(frozen=True)
class WSITile:
    """Represent one WSI grid location and its retained-image state.

    Attributes:
        tile_id: Stable identifier for the WSI level and grid coordinates.
        wsi_id: Caller-defined WSI identifier.
        level: OpenSlide level used for the read.
        x: Tile x-coordinate at the selected level.
        y: Tile y-coordinate at the selected level.
        x_level0: Corresponding x-coordinate at level zero.
        y_level0: Corresponding y-coordinate at level zero.
        w_real: Observed tile width before edge padding.
        h_real: Observed tile height before edge padding.
        pad_right: Right-padding width in pixels.
        pad_bottom: Bottom-padding height in pixels.
        status: ``kept`` or the assigned white-area skip reason.
        image: Padded normalized RGB tile for kept records, otherwise ``None``.
    """

    tile_id: str
    wsi_id: str
    level: int
    x: int
    y: int
    x_level0: int
    y_level0: int
    w_real: int
    h_real: int
    pad_right: int
    pad_bottom: int
    status: str
    image: np.ndarray | None


def open_wsi(path: str | Path):
    """Open one whole-slide image through OpenSlide.

    Parameters:
        path: Filesystem path to the WSI.

    Returns:
        The opened OpenSlide handle.
    """
    # Import OpenSlide only when a caller opens a whole-slide image.
    try:
        import openslide
    except ImportError as error:
        raise ImportError("open_wsi requires openslide-python") from error
    return openslide.OpenSlide(str(path))


def iter_wsi_tiles(
    slide,
    *,
    wsi_id: str,
    level: int,
    macenko_normalizer: Callable[[np.ndarray], np.ndarray],
    tile_size: int = TILE_SIZE,
    white_threshold: int = 230,
    maximum_white_fraction: float = 0.80,
) -> Iterator[WSITile]:
    """Yield every selected-level grid position with keep/skip metadata.

    Parameters:
        slide: OpenSlide handle returned by ``open_wsi``.
        wsi_id: Stable identifier recorded in every yielded tile.
        level: OpenSlide level selected for the inference magnification.
        macenko_normalizer: Callable returning a same-shape uint8 RGB array.
        tile_size: Required side length for SpMap inference tiles.
        white_threshold: Per-channel inclusive RGB threshold for white pixels.
        maximum_white_fraction: Largest retained fraction of white pixels.

    Yields:
        ``WSITile`` records for all grid locations in row-major order.
    """
    # Validate the fixed grid size and caller-selected OpenSlide level.
    if tile_size != TILE_SIZE:
        raise ValueError(f"SpMap WSI tiles must be {TILE_SIZE} x {TILE_SIZE}")
    if level < 0 or level >= slide.level_count:
        raise ValueError(f"Requested level {level}; slide has {slide.level_count} levels")
    # Read selected-level geometry and its mapping back to level-zero coordinates.
    width, height = map(int, slide.level_dimensions[level])
    downsample = float(slide.level_downsamples[level])

    # Traverse every grid position in row-major order, including right/bottom edges.
    for y in range(0, height, tile_size):
        for x in range(0, width, tile_size):
            w_real = min(tile_size, width - x)
            h_real = min(tile_size, height - y)
            x_level0 = int(round(x * downsample))
            y_level0 = int(round(y * downsample))
            tile_id = f"{wsi_id}:L{level}:x{x}:y{y}"
            # Read observed edge dimensions while addressing OpenSlide at level zero.
            observed = np.asarray(
                slide.read_region(
                    (x_level0, y_level0), level, (w_real, h_real)
                ).convert("RGB"),
                dtype=np.uint8,
            )
            pad_right = tile_size - w_real
            pad_bottom = tile_size - h_real
            # Classify whiteness on observed pixels before normalization or padding.
            if (
                white_fraction(observed, white_threshold=white_threshold)
                > maximum_white_fraction
            ):
                yield WSITile(
                    tile_id,
                    str(wsi_id),
                    level,
                    x,
                    y,
                    x_level0,
                    y_level0,
                    w_real,
                    h_real,
                    pad_right,
                    pad_bottom,
                    "more_than_80_percent_white",
                    None,
                )
                continue

            # Normalize retained tissue without changing its observed edge geometry.
            normalized = np.asarray(macenko_normalizer(observed))
            if normalized.shape != observed.shape or normalized.dtype != np.uint8:
                raise ValueError(
                    "macenko_normalizer must return an RGB uint8 array with unchanged shape"
                )
            # Right/bottom-pad retained edge tiles to the fixed inference input size.
            padded = np.zeros((tile_size, tile_size, 3), dtype=np.uint8)
            padded[:h_real, :w_real] = normalized
            yield WSITile(
                tile_id,
                str(wsi_id),
                level,
                x,
                y,
                x_level0,
                y_level0,
                w_real,
                h_real,
                pad_right,
                pad_bottom,
                "kept",
                padded,
            )

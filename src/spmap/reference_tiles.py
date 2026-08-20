"""Construct SpMap reference regions and complete non-overlapping tiles.

Purpose:
    Convert caller-selected registered structures or ROIs into clipped H&E
    regions, then enumerate and extract complete 256 x 256 training tiles.
Figure 5 callers:
    Figure 5 SpMap reference-tile preparation workflows use these helpers
    before preprocessing, feature extraction, and grouped splitting.
Inputs:
    Registered-cell tables or selected regions, source-image dimensions, and
    in-memory RGB images for extraction.
Outputs:
    Bounding-box and tile-manifest DataFrames plus contiguous RGB tile arrays.
Ordered use:
    Build structure boxes where needed, combine caller-selected regions, build
    the complete manifest, and extract tiles from the source image.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


TILE_SIZE = 256
CLASS_NAMES = {0: "TC", 1: "PIR", 2: "PSM", 3: "OTHER_PT"}
BBOX_COLUMNS = (
    "parent_id",
    "source_image_id",
    "class_id",
    "class_name",
    "x",
    "y",
    "width",
    "height",
)
TILE_COLUMNS = (
    "tile_id",
    "class_id",
    "class_name",
    "parent_id",
    "source_image_id",
    "x",
    "y",
    "width",
    "height",
)


def _require_columns(frame: pd.DataFrame, columns: set[str]) -> None:
    """Require the named columns in a caller-supplied DataFrame.

    Parameters:
        frame: DataFrame validated before a reference-tile operation.
        columns: Required column-name set.

    Returns:
        ``None`` when every required column is present.
    """
    missing = columns - set(frame.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")


def _single_value(group: pd.DataFrame, column: str, parent_id: str):
    """Read the sole allowed group-level value for one metadata column.

    Parameters:
        group: One structure's registered-cell rows.
        column: Metadata column expected to have one distinct value.
        parent_id: Structure identifier included in validation errors.

    Returns:
        The group's unique value for ``column``.
    """
    # Preserve the invariant that one structure has one source and one class.
    values = group[column].drop_duplicates()
    if len(values) != 1:
        raise ValueError(f"{parent_id}: expected one {column}, observed {len(values)}")
    return values.iloc[0]


def _clip_preserving_size(
    x: int,
    y: int,
    width: int,
    height: int,
    image_width: int,
    image_height: int,
) -> tuple[int, int, int, int]:
    """Clip one rectangle to image bounds while retaining its clipped size.

    Parameters:
        x: Rectangle left coordinate.
        y: Rectangle top coordinate.
        width: Rectangle width.
        height: Rectangle height.
        image_width: Source image width.
        image_height: Source image height.

    Returns:
        Clipped ``(x, y, width, height)`` integer rectangle.
    """
    # Limit box dimensions first so subsequent coordinate clipping retains them.
    width = min(width, image_width)
    height = min(height, image_height)
    # Shift the fixed-size box inside the source-image coordinate bounds.
    x = min(max(x, 0), image_width - width)
    y = min(max(y, 0), image_height - height)
    return x, y, width, height


def build_structure_bboxes(
    cells: pd.DataFrame,
    *,
    image_size: tuple[int, int],
    structure_col: str = "structure_id",
    source_image_col: str = "source_image_id",
    class_col: str = "class_id",
    x_col: str = "x",
    y_col: str = "y",
    coordinate_scale: float = 2.0,
    padding: float = 10.0,
    minimum_size: tuple[int, int] = (256, 256),
    snap: int = 32,
) -> pd.DataFrame:
    """Convert registered-cell groups into clipped H&E bounding boxes.

    Parameters:
        cells: Registered-cell table with structure, image, class, and position columns.
        image_size: Target H&E image width and height in pixels.
        structure_col: Column grouping cells into caller-selected structures.
        source_image_col: Column identifying each H&E source image.
        class_col: Column containing canonical class IDs.
        x_col: Registered x-coordinate column.
        y_col: Registered y-coordinate column.
        coordinate_scale: Factor converting registered coordinates to H&E pixels.
        padding: Registered-coordinate padding around each structure.
        minimum_size: Minimum width and height before coordinate scaling.
        snap: Registered-coordinate grid spacing for box edges and sizes.

    Returns:
        Bounding-box DataFrame with the fixed ``BBOX_COLUMNS`` schema.
    """
    # Validate the registered-cell schema and target image geometry once up front.
    required = {structure_col, source_image_col, class_col, x_col, y_col}
    _require_columns(cells, required)
    image_width, image_height = map(int, image_size)
    if image_width <= 0 or image_height <= 0:
        raise ValueError("image_size values must be positive")
    if coordinate_scale <= 0 or snap <= 0:
        raise ValueError("coordinate_scale and snap must be positive")

    # Build one class-preserving region from each registered structure group.
    rows = []
    for structure_id, group in cells.groupby(structure_col, sort=False):
        parent_id = str(structure_id)
        source_image_id = str(_single_value(group, source_image_col, parent_id))
        class_id = int(_single_value(group, class_col, parent_id))
        if class_id not in CLASS_NAMES:
            raise ValueError(f"{parent_id}: invalid class_id {class_id}")

        # Expand cell extrema, snap edges and dimensions, and enforce minimum size.
        x_min = float(group[x_col].min()) - padding
        x_max = float(group[x_col].max()) + padding
        y_min = float(group[y_col].min()) - padding
        y_max = float(group[y_col].max()) + padding
        x = math.floor(x_min / snap) * snap
        y = math.floor(y_min / snap) * snap
        width = max(math.ceil((x_max - x) / snap) * snap, int(minimum_size[0]))
        height = max(math.ceil((y_max - y) / snap) * snap, int(minimum_size[1]))

        # Convert registered coordinates to H&E pixels before image-bound clipping.
        x = int(round(x * coordinate_scale))
        y = int(round(y * coordinate_scale))
        width = int(round(width * coordinate_scale))
        height = int(round(height * coordinate_scale))
        x, y, width, height = _clip_preserving_size(
            x, y, width, height, image_width, image_height
        )
        # Emit the fixed bbox schema with the canonical class label.
        rows.append(
            {
                "parent_id": parent_id,
                "source_image_id": source_image_id,
                "class_id": class_id,
                "class_name": CLASS_NAMES[class_id],
                "x": x,
                "y": y,
                "width": width,
                "height": height,
            }
        )
    return pd.DataFrame(rows, columns=BBOX_COLUMNS)


def build_complete_tile_manifest(
    regions: pd.DataFrame,
    *,
    tile_size: int = TILE_SIZE,
) -> pd.DataFrame:
    """Enumerate complete tiles from caller-selected structures or ROIs.

    Parameters:
        regions: Region table with parent, source, class, and rectangle columns.
        tile_size: Required side length of each complete SpMap tile.

    Returns:
        Tile-manifest DataFrame following the fixed ``TILE_COLUMNS`` schema.
    """
    # Validate region geometry and retain the fixed 256-pixel SpMap tile size.
    required = {
        "parent_id",
        "source_image_id",
        "class_id",
        "x",
        "y",
        "width",
        "height",
    }
    _require_columns(regions, required)
    if tile_size != TILE_SIZE:
        raise ValueError(f"SpMap reference tiles must be {TILE_SIZE} x {TILE_SIZE}")

    # Enumerate only complete non-overlapping tiles within each selected region.
    rows = []
    for region in regions.itertuples(index=False):
        values = region._asdict()
        parent_id = str(values["parent_id"])
        source_image_id = str(values["source_image_id"])
        class_id = int(values["class_id"])
        if class_id not in CLASS_NAMES:
            raise ValueError(f"{parent_id}: invalid class_id {class_id}")
        x0, y0 = int(values["x"]), int(values["y"])
        width, height = int(values["width"]), int(values["height"])
        # Use absolute image coordinates in stable row-major tile identifiers.
        for dy in range(0, height - tile_size + 1, tile_size):
            for dx in range(0, width - tile_size + 1, tile_size):
                x, y = x0 + dx, y0 + dy
                tile_id = f"{source_image_id}:{parent_id}:x{x}:y{y}"
                rows.append(
                    {
                        "tile_id": tile_id,
                        "class_id": class_id,
                        "class_name": CLASS_NAMES[class_id],
                        "parent_id": parent_id,
                        "source_image_id": source_image_id,
                        "x": x,
                        "y": y,
                        "width": tile_size,
                        "height": tile_size,
                    }
                )
    # Enforce globally unique tile identities across all parent regions.
    manifest = pd.DataFrame(rows, columns=TILE_COLUMNS)
    if not manifest.empty and manifest["tile_id"].duplicated().any():
        duplicates = manifest.loc[manifest["tile_id"].duplicated(), "tile_id"]
        raise ValueError(f"Duplicate tile_id: {duplicates.iloc[0]}")
    return manifest


def extract_tile(rgb_image: np.ndarray, tile_row: pd.Series) -> np.ndarray:
    """Extract one manifest-defined complete tile from an in-memory RGB image.

    Parameters:
        rgb_image: Source H x W x 3 uint8 RGB image array.
        tile_row: Manifest row containing tile coordinates and dimensions.

    Returns:
        Contiguous H x W x 3 uint8 tile array with manifest dimensions.
    """
    # Normalize the source array while preserving the RGB byte-image contract.
    image = np.asarray(rgb_image)
    if image.ndim != 3 or image.shape[2] != 3 or image.dtype != np.uint8:
        raise ValueError("rgb_image must be an H x W x 3 uint8 array")
    # Slice the manifest rectangle and require its complete requested extent.
    x, y = int(tile_row["x"]), int(tile_row["y"])
    width, height = int(tile_row["width"]), int(tile_row["height"])
    tile = image[y : y + height, x : x + width]
    if tile.shape != (height, width, 3):
        raise ValueError(f"{tile_row['tile_id']}: tile extends beyond the source image")
    return np.ascontiguousarray(tile)

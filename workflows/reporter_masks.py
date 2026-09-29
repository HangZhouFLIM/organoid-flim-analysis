"""Selected original reporter-mask functions; explicit parameters, no workstation paths.

Upstream SHA256: f7140a6897fbfd21bb74905282b3c994a27271fc640a225f4e64ca36229c918e
"""

from __future__ import annotations

from dataclasses import dataclass

from pathlib import Path

import numpy as np

import pandas as pd

from scipy import ndimage

from skimage import filters,measure,morphology

IMAGE_SHAPE=(246,246)

@dataclass
class GfpGroupResult:
    method: str
    low_mask: np.ndarray
    intermediate_mask: np.ndarray
    high_mask: np.ndarray
    p25: float
    p75: float
    raw_low_pixels: int
    raw_high_pixels: int
    low_pixels: int
    intermediate_pixels: int
    high_pixels: int
    high_initial_objects: int
    high_final_objects: int
    high_removed_small_objects: int
    high_removed_small_pixels: int
    high_filled_hole_pixels: int

def read_channel(path: Path, required_metrics: list[str]) -> pd.DataFrame:
    required = set(required_metrics) | {"Ycoord", "Xcoord"}
    df = pd.read_csv(path)
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")
    if df[["Ycoord", "Xcoord"]].duplicated().any():
        raise ValueError(f"{path} contains duplicate X/Y coordinates")
    return df[list(required)]

def reconstruct_images(
    df: pd.DataFrame,
    columns: list[str],
    shape: tuple[int, int] = IMAGE_SHAPE,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    y = df["Ycoord"].to_numpy(dtype=int) - 1
    x = df["Xcoord"].to_numpy(dtype=int) - 1
    if y.min() < 0 or x.min() < 0 or y.max() >= shape[0] or x.max() >= shape[1]:
        raise ValueError(
            f"Coordinates exceed image shape {shape}: "
            f"X={x.min() + 1}..{x.max() + 1}, Y={y.min() + 1}..{y.max() + 1}"
        )
    valid = np.zeros(shape, dtype=bool)
    valid[y, x] = True
    images: dict[str, np.ndarray] = {}
    for column in columns:
        image = np.full(shape, np.nan, dtype=float)
        image[y, x] = df[column].to_numpy(dtype=float)
        images[column] = image
    return images, valid

def object_summary(mask: np.ndarray) -> tuple[int, np.ndarray]:
    labeled = measure.label(mask, connectivity=2)
    sizes = np.bincount(labeled.ravel())
    object_sizes = sizes[1:] if len(sizes) > 1 else np.array([], dtype=int)
    return int(object_sizes.size), object_sizes

def cleanup_mask(
    initial: np.ndarray,
    valid: np.ndarray,
    min_object_size: int,
    hole_size: int,
    closing_radius: int = 1,
) -> tuple[np.ndarray, dict[str, int]]:
    initial = (initial & valid).astype(bool)
    initial_objects, object_sizes = object_summary(initial)
    small_sizes = object_sizes[object_sizes < min_object_size]
    without_small = morphology.remove_small_objects(
        initial, min_size=min_object_size, connectivity=2
    )
    without_small &= valid
    filled = morphology.remove_small_holes(
        without_small, area_threshold=hole_size, connectivity=2
    )
    filled &= valid
    filled_hole_pixels = int(filled.sum() - without_small.sum())
    if closing_radius > 0:
        closed = morphology.binary_closing(filled, morphology.disk(closing_radius))
    else:
        closed = filled
    closed &= valid
    final_objects, _ = object_summary(closed)
    return closed, {
        "initial_objects": int(initial_objects),
        "final_objects": int(final_objects),
        "removed_small_objects": int(small_sizes.size),
        "removed_small_pixels": int(small_sizes.sum()) if small_sizes.size else 0,
        "filled_hole_pixels": max(filled_hole_pixels, 0),
    }

def normalize_valid(values: np.ndarray, valid: np.ndarray, p_low: float = 1.0, p_high: float = 99.8) -> np.ndarray:
    finite = valid & np.isfinite(values)
    out = np.zeros(values.shape, dtype=float)
    if not np.any(finite):
        return out
    lo, hi = np.nanpercentile(values[finite], [p_low, p_high])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return out
    out[finite] = np.clip((values[finite] - lo) / (hi - lo), 0.0, 1.0)
    return out

def make_infer_subc_mito_mask(
    tmrm_image: np.ndarray,
    valid: np.ndarray,
    median_size: int,
    gaussian_sigma: float,
    dot_sigma: float,
    dot_cutoff: float,
    filament_sigma: float,
    filament_cutoff: float,
    min_object_size: int,
    hole_size: int,
    closing_radius: int,
) -> tuple[np.ndarray, dict[str, int], str, int]:
    normalized = normalize_valid(tmrm_image, valid)
    if median_size > 1:
        normalized = ndimage.median_filter(normalized, size=median_size, mode="nearest")
    smooth = ndimage.gaussian_filter(normalized, sigma=gaussian_sigma, mode="nearest")
    dot_response = -ndimage.gaussian_laplace(smooth, sigma=dot_sigma, mode="nearest")
    dot_response = normalize_valid(np.clip(dot_response, 0.0, None), valid, p_low=1.0, p_high=99.5)
    filament_response = filters.sato(
        smooth,
        sigmas=[filament_sigma],
        black_ridges=False,
        mode="nearest",
    )
    filament_response = normalize_valid(filament_response, valid, p_low=1.0, p_high=99.5)
    initial = valid & ((dot_response >= dot_cutoff) | (filament_response >= filament_cutoff))
    mask, stats_dict = cleanup_mask(
        initial,
        valid,
        min_object_size=min_object_size,
        hole_size=hole_size,
        closing_radius=closing_radius,
    )
    threshold = (
        "infer-subc 2.0.1-style mitochondria; "
        f"median={median_size}, gaussian={gaussian_sigma:g}, "
        f"dot_sigma={dot_sigma:g}, dot_cutoff={dot_cutoff:g}, "
        f"filament_sigma={filament_sigma:g}, filament_cutoff={filament_cutoff:g}"
    )
    return mask, stats_dict, threshold, int(initial.sum())

def make_gfp_groups(
    gfp_image: np.ndarray,
    valid: np.ndarray,
    min_object_size: int,
    hole_size: int,
) -> GfpGroupResult:
    values = gfp_image[valid]
    p25 = float(np.nanpercentile(values, 25))
    p75 = float(np.nanpercentile(values, 75))
    raw_low = valid & (gfp_image <= p25)
    raw_high = valid & (gfp_image >= p75)
    high_mask, high_stats = cleanup_mask(
        raw_high, valid, min_object_size=min_object_size, hole_size=hole_size
    )
    low_mask = raw_low & ~high_mask
    intermediate_mask = valid & ~(low_mask | high_mask)
    return GfpGroupResult(
        method="gfp_low_medium_high_p25_p75_high_clean",
        low_mask=low_mask,
        intermediate_mask=intermediate_mask,
        high_mask=high_mask,
        p25=p25,
        p75=p75,
        raw_low_pixels=int(raw_low.sum()),
        raw_high_pixels=int(raw_high.sum()),
        low_pixels=int(low_mask.sum()),
        intermediate_pixels=int(intermediate_mask.sum()),
        high_pixels=int(high_mask.sum()),
        high_initial_objects=high_stats["initial_objects"],
        high_final_objects=high_stats["final_objects"],
        high_removed_small_objects=high_stats["removed_small_objects"],
        high_removed_small_pixels=high_stats["removed_small_pixels"],
        high_filled_hole_pixels=high_stats["filled_hole_pixels"],
    )

"""Initial geometry and positive-fluorescence constraints (NumPy/SciPy only).

All intensities and thresholds are camera counts on the original 0..255 scale.
Functions do not modify their input arrays. Coordinates are source-pixel (x, y),
while ``yx`` indexes arrays as (row, column). Defaults reproduce the reviewed
experiment settings; they are parameters, not universal calibration values.

Configuration may contain ``grain``, ``wall`` and ``fluorescence`` dictionaries.
Flat keys are also accepted. Geometry: ``bed_polygon_xy`` or ``roi_xyxy``;
otherwise the complete image is the bed. Grain defaults: window=101, k=.25,
r=.5, p=2, q=10, sigma=.8, min_area=25, intensity_scale=255. Wall defaults:
normal_sigma=1.1, channel_sigma=.9, bed_clearance=10, profile_chunk_size=200000,
alpha=1.1; lit thresholds RB=14, D=12 with RB>6, B=10 with RB>5; unlit
thresholds RB=6, B=5, D=6. Fluorescence defaults: threshold=55, sigma=.9,
shiny_distance=4.5, shiny_contrast=10, outer_margin=10. No grain dilation is
applied to the excluded mask; an immediate 8-neighbor rim samples wall profiles.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def _options(cfg, section):
    out = dict(cfg or {})
    nested = out.get(section, {})
    if isinstance(nested, dict):
        out.update(nested)
    return out


def _image(a, name):
    a = np.asarray(a)
    if a.ndim != 2 or not a.size:
        raise ValueError(f"{name} must be a nonempty 2-D grayscale image")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} contains nonfinite intensities")
    return a


def _finite(value, name, positive=False, nonnegative=False):
    value = float(value)
    if not np.isfinite(value) or (positive and value <= 0) or (nonnegative and value < 0):
        raise ValueError(f"Invalid {name}: {value}")
    return value


def _bed_mask(shape, cfg):
    """Rasterize a polygon without an additional graphics dependency."""
    h, w = shape
    polygon = cfg.get("bed_polygon_xy")
    roi = cfg.get("roi_xyxy")
    if polygon is not None and roi is not None:
        raise ValueError("Choose either bed_polygon_xy or roi_xyxy")
    if roi is not None:
        if len(roi) != 4:
            raise ValueError("roi_xyxy needs [x0, y0, x1, y1]")
        if any(int(v) != v for v in roi):
            raise ValueError("roi_xyxy coordinates must be integers")
        x0, y0, x1, y1 = map(int, roi)
        if not (0 <= x0 < x1 <= w and 0 <= y0 < y1 <= h):
            raise ValueError("roi_xyxy lies outside the image")
        bed = np.zeros(shape, bool)
        bed[y0:y1, x0:x1] = True
        return bed, {"method": "rectangle", "roi_xyxy": [x0, y0, x1, y1]}
    if polygon is None:
        return np.ones(shape, bool), {"method": "full image"}
    points = np.asarray(polygon, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 3 or not np.all(np.isfinite(points)):
        raise ValueError("bed_polygon_xy needs at least three finite [x, y] vertices")
    if np.any(points[:, 0] < 0) or np.any(points[:, 0] > w) or np.any(points[:, 1] < 0) or np.any(points[:, 1] > h):
        raise ValueError("bed_polygon_xy lies outside the image")
    yy = np.arange(h)[:, None]
    xx = np.arange(w)[None, :]
    bed = np.zeros(shape, bool)
    for start, end in zip(points, np.roll(points, -1, axis=0)):
        x0, y0 = start
        x1, y1 = end
        if y1 == y0:
            continue
        crossing_rows = (y0 > yy) != (y1 > yy)
        crossing_x = x0 + (yy - y0) * ((x1 - x0) / (y1 - y0))
        bed ^= crossing_rows & (xx < crossing_x)
    if not bed.any():
        raise ValueError("bed_polygon_xy contains no image pixels")
    return bed, {"method": "polygon", "bed_polygon_xy": points.tolist()}


def make_masks(initial, cfg=None):
    """Return undilated PDMS, pore and bed masks plus JSON-safe geometry details.

    Phansalkar: t=m*(1+p*exp(-q*m)+k*(s/r-1)), with normalized intensities.
    Fill enclosed grain holes, remove components smaller than min_area, and
    retain pore space immediately next to grains for interface detection.
    """
    initial = _image(initial, "initial")
    opts = _options(cfg, "grain")
    window = opts.get("window", 101)
    if int(window) != window or int(window) < 3 or int(window) % 2 != 1:
        raise ValueError("grain.window must be an odd integer >=3")
    minimum = opts.get("min_area", opts.get("min_grain_area", 25))
    if int(minimum) != minimum or int(minimum) < 1:
        raise ValueError("grain.min_area must be an integer >=1")
    k = _finite(opts.get("k", .25), "grain.k")
    r = _finite(opts.get("r", .5), "grain.r", positive=True)
    p = _finite(opts.get("p", 2), "grain.p", nonnegative=True)
    q = _finite(opts.get("q", 10), "grain.q", nonnegative=True)
    sigma = _finite(opts.get("sigma", .8), "grain.sigma", nonnegative=True)
    scale = _finite(opts.get("intensity_scale", 255), "grain.intensity_scale", positive=True)
    bed, info = _bed_mask(initial.shape, opts)
    normalized = ndi.gaussian_filter(initial.astype(np.float32), sigma) / scale
    mean = ndi.uniform_filter(normalized, size=int(window), mode="reflect")
    variance = ndi.uniform_filter(normalized * normalized, size=int(window), mode="reflect")
    variance -= mean * mean
    np.maximum(variance, 0, out=variance)
    np.sqrt(variance, out=variance)
    threshold = mean * (1 + p * np.exp(-q * mean) + k * (variance / r - 1))
    grain = ndi.binary_fill_holes((normalized < threshold) & bed) & bed
    labels, _ = ndi.label(grain)
    sizes = np.bincount(labels.ravel())
    grain = (labels > 0) & (sizes[labels] >= int(minimum))
    pore = bed & ~grain
    info.update({"grain_method": "Phansalkar", "window": int(window), "k": k,
                 "r": r, "p": p, "q": q, "sigma": sigma,
                 "min_area": int(minimum), "intensity_scale": scale,
                 "grain_dilation_px": 0, "bed_pixels": int(bed.sum()),
                 "grain_pixels": int(grain.sum()), "pore_pixels": int(pore.sum()),
                 "warnings": []})
    if not grain.any():
        info["warnings"].append("No grains detected; grain-wall evidence is unavailable.")
    if not pore.any():
        info["warnings"].append("No pore pixels remain in the selected bed.")
    return grain, pore, bed, info


def _masked_column_median(values, mask):
    """Small (2 or 3)-sample column medians without all-NaN warnings."""
    count = np.sum(mask, axis=0)
    ordered = np.sort(np.where(mask, values, np.inf), axis=0)
    index = np.arange(values.shape[1])
    lo = ordered[np.maximum(count - 1, 0) // 2, index]
    hi = ordered[count // 2, index]
    lo = np.where(count > 0, lo, 0)
    hi = np.where(count > 0, hi, 0)
    return (lo + hi) * .5


def wall_profiles(B, RB, grain, pore, cfg=None):
    """Sample original channels across initial PDMS walls.

    Positive normal direction points from PDMS into pore space. Profile
    validity needs >=1 grain-side and >=2 outer-pore samples, normal norm>.9,
    and bed clearance>10px. D is sampled RB-alpha*B; subtraction is never
    unsigned. Returned yx/features are vectors with identical ordering.
    """
    B, RB = _image(B, "B"), _image(RB, "RB")
    grain, pore = np.asarray(grain, bool), np.asarray(pore, bool)
    if B.shape != RB.shape or B.shape != grain.shape or B.shape != pore.shape:
        raise ValueError("B, RB, grain and pore must have identical shapes")
    if np.any(grain & pore):
        raise ValueError("grain and pore masks overlap")
    opts = _options(cfg, "wall")
    alpha = _finite(opts.get("alpha", opts.get("blue_scale", 1.1)), "alpha", nonnegative=True)
    normal_sigma = _finite(opts.get("normal_sigma", 1.1), "wall.normal_sigma", nonnegative=True)
    channel_sigma = _finite(opts.get("channel_sigma", .9), "wall.channel_sigma", nonnegative=True)
    clearance = _finite(opts.get("bed_clearance", 10), "wall.bed_clearance", nonnegative=True)
    chunk_size = int(opts.get("profile_chunk_size", 200000))
    if chunk_size < 1:
        raise ValueError("wall.profile_chunk_size must be >=1")
    rim = pore & ndi.binary_dilation(grain, structure=np.ones((3, 3), bool))
    yy, xx = np.where(rim)
    offsets = np.array([-4, -2, 0, 1, 2, 3, 4, 6, 8], np.float32)
    n = len(yy)
    result = {"yx": np.column_stack((yy, xx)).astype(np.int32),
              "offsets": offsets, "valid": np.zeros(n, bool),
              "lit": np.zeros(n, bool), "unlit": np.zeros(n, bool),
              "normals_xy": np.zeros((n, 2), np.float32), "warnings": []}
    for channel in ("B", "RB", "D"):
        for feature in ("contrast", "near_peak", "outer", "inner_grain"):
            result[channel + "_" + feature] = np.zeros(n, np.float32)
    if n:
        sd = ndi.distance_transform_edt(~grain).astype(np.float32)
        sd -= ndi.distance_transform_edt(grain).astype(np.float32)
        sd = ndi.gaussian_filter(sd, normal_sigma)
        ny, nx = np.gradient(sd)
        del sd
        norm = np.maximum(np.hypot(nx, ny), 1e-6)
        vx, vy = nx[yy, xx] / norm[yy, xx], ny[yy, xx] / norm[yy, xx]
        del nx, ny, norm
        result["normals_xy"] = np.column_stack((vx, vy)).astype(np.float32)
        bed_distance = ndi.distance_transform_edt(grain | pore)
        valid_clearance = bed_distance[yy, xx] > clearance
        del bed_distance
        # A full-image bed has no exterior zeros; explicitly protect its edges.
        image_clearance = np.minimum.reduce([yy + 1, xx + 1, B.shape[0] - yy, B.shape[1] - xx])
        valid_clearance &= image_clearance > clearance
        bs = ndi.gaussian_filter(B.astype(np.float32), channel_sigma)
        rs = ndi.gaussian_filter(RB.astype(np.float32), channel_sigma)
        grain_u8, pore_u8 = grain.astype(np.uint8), pore.astype(np.uint8)
        for start in range(0, n, chunk_size):
            end = min(n, start + chunk_size)
            sl = slice(start, end)
            cy = yy[None, sl] + offsets[:, None] * vy[None, sl]
            cx = xx[None, sl] + offsets[:, None] * vx[None, sl]
            gp = ndi.map_coordinates(grain_u8, [cy, cx], order=0, mode="nearest") > 0
            pp = ndi.map_coordinates(pore_u8, [cy, cx], order=0, mode="nearest") > 0
            result["valid"][sl] = ((gp[:2].sum(0) >= 1) & (pp[6:].sum(0) >= 2)
                                  & (np.hypot(vx[sl], vy[sl]) > .9) & valid_clearance[sl])
            bp = ndi.map_coordinates(bs, [cy, cx], order=1, mode="nearest")
            rp = ndi.map_coordinates(rs, [cy, cx], order=1, mode="nearest")
            for name, values in (("B", bp), ("RB", rp), ("D", rp - alpha * bp)):
                peak = values[2:6].max(axis=0)
                outer = _masked_column_median(values[6:], pp[6:])
                inner = _masked_column_median(values[:2], gp[:2])
                baseline_exists = gp[:2].any(axis=0) & pp[6:].any(axis=0)
                result[name + "_contrast"][sl] = np.where(
                    baseline_exists, peak - np.maximum(outer, inner), 0)
                result[name + "_near_peak"][sl] = peak
                result[name + "_outer"][sl] = outer
                result[name + "_inner_grain"][sl] = inner
        rb, b, d = result["RB_contrast"], result["B_contrast"], result["D_contrast"]
        result["lit"] = result["valid"] & ((rb > opts.get("lit_rb", 14))
            | ((d > opts.get("lit_d", 12)) & (rb > opts.get("lit_d_min_rb", 6)))
            | ((b > opts.get("lit_b", 10)) & (rb > opts.get("lit_b_min_rb", 5))))
        result["unlit"] = (result["valid"] & (rb < opts.get("unlit_rb", 6))
                           & (b < opts.get("unlit_b", 5)) & (d < opts.get("unlit_d", 6)))
        score = np.clip(np.maximum.reduce([rb / 14, d / 12, b / 10]), 0, 3) / 3
    else:
        score = np.zeros(0, np.float32)
        result["warnings"].append("No grain/pore rim; wall evidence and shiny-edge exceptions are unavailable.")
    if n and not result["valid"].any():
        result["warnings"].append("No valid wall profiles; shiny-edge exceptions are unavailable.")
    for map_name, vector in (("wall_lit", result["lit"]), ("wall_unlit", result["unlit"]),
                             ("wall_valid", result["valid"]), ("wall_score", score)):
        full = np.zeros(B.shape, np.float32 if map_name == "wall_score" else bool)
        full[yy, xx] = vector
        result[map_name] = full
    return result


def fluorescence_water(B, grain, pore, profiles, cfg=None, threshold=None):
    """Return direct-water, shiny-exception, smoothed-B arrays.

    A positive B observation is water, except for a narrow peak no more than
    4.5px from a *valid* wall profile with contrast>10 and outer-pore B<T-10.
    This exception only withholds a water seed; it does not classify air.
    Unassigned dark pixels can be either diluted water or air.
    """
    B = _image(B, "B")
    grain, pore = np.asarray(grain, bool), np.asarray(pore, bool)
    if B.shape != grain.shape or B.shape != pore.shape:
        raise ValueError("B, grain and pore must have identical shapes")
    if np.any(grain & pore):
        raise ValueError("grain and pore masks overlap")
    opts = _options(cfg, "fluorescence")
    T = _finite(opts.get("threshold", 55) if threshold is None else threshold, "fluorescence.threshold")
    sigma = _finite(opts.get("sigma", .9), "fluorescence.sigma", nonnegative=True)
    distance = _finite(opts.get("shiny_distance", 4.5), "fluorescence.shiny_distance", nonnegative=True)
    contrast = _finite(opts.get("shiny_contrast", 10), "fluorescence.shiny_contrast", nonnegative=True)
    margin = _finite(opts.get("outer_margin", 10), "fluorescence.outer_margin", nonnegative=True)
    bs = ndi.gaussian_filter(B.astype(np.float32), sigma)
    bright = pore & (bs > T)
    shiny = np.zeros(B.shape, bool)
    yx = np.asarray(profiles["yx"], dtype=np.int32)
    valid = np.asarray(profiles["valid"], bool)
    if len(yx) != len(valid):
        raise ValueError("Profile yx and valid vectors have different lengths")
    if valid.any() and bright.any():
        wall = np.zeros(B.shape, bool)
        wall[yx[valid, 0], yx[valid, 1]] = True
        nearest = ndi.distance_transform_edt(~wall, return_distances=False, return_indices=True)
        # Evaluate only bright pixels to avoid full-frame nearest-feature arrays.
        yy, xx = np.where(bright)
        ny, nx = nearest[0, yy, xx], nearest[1, yy, xx]
        close = ((yy - ny) ** 2 + (xx - nx) ** 2) <= distance * distance
        del nearest
        feature_index = np.full(B.shape, -1, np.int32)
        feature_index[yx[valid, 0], yx[valid, 1]] = np.flatnonzero(valid)
        idx = feature_index[ny, nx]
        del feature_index
        is_shiny = (close & (np.asarray(profiles["B_contrast"])[idx] > contrast)
                    & (np.asarray(profiles["B_outer"])[idx] < T - margin))
        shiny[yy[is_shiny], xx[is_shiny]] = True
    return bright & ~shiny, shiny, bs

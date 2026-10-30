"""Read a fixed grain reference once, independently of the image-pair series.

``mask_input.mode`` is either ``binarize`` (Phansalkar on a separate camera
image) or ``load_binary`` (use a supplied grain mask exactly). Binary inputs
are never resized, interpolated, thresholded, filled, or cleaned. Paths may be
absolute or relative to ``input_folder``. TIFF and single-frame PNG are read
without changing original files.
"""
from __future__ import annotations

from pathlib import Path
import hashlib

import numpy as np
from scipy import ndimage as ndi
import tifffile

from .geometry import _bed_mask, make_masks


def _resolve_path(folder, value, label):
    if value is None or not str(value).strip():
        raise ValueError(f"Set mask_input.{label} to a file path")
    path = Path(value).expanduser()
    candidate = path if path.is_absolute() else folder / path
    if candidate.is_file():
        return candidate.resolve()
    # Preserve the original package's unique TIFF-suffix convenience.
    if not path.suffix and folder.is_dir():
        found = sorted(p for p in folder.iterdir() if p.is_file()
                       and p.suffix.lower() in {".tif", ".tiff"}
                       and p.stem.endswith("_" + str(value)))
        if len(found) == 1:
            return found[0].resolve()
        if len(found) > 1:
            raise ValueError(f"Ambiguous {label}: {len(found)} TIFFs end _{value}; use a complete filename")
    raise FileNotFoundError(f"Reference file does not exist: {candidate}")


def _read_2d(path):
    if path.suffix.lower() in {".tif", ".tiff"}:
        image = tifffile.imread(path)
    elif path.suffix.lower() == ".png":
        from PIL import Image
        with Image.open(path) as source:
            if getattr(source, "n_frames", 1) != 1:
                raise ValueError(f"Expected a single image, not a sequence: {path}")
            image = np.asarray(source)
    else:
        raise ValueError(f"Reference files must be TIFF or PNG: {path}")
    if image.ndim != 2 or not image.size or min(image.shape) < 3:
        raise ValueError(f"Expected one 2-D grayscale reference: {path}; got shape {image.shape}")
    if not np.all(np.isfinite(image)):
        raise ValueError(f"Nonfinite reference values: {path}")
    return image


def _camera_image(path):
    image = _read_2d(path)
    if image.min() < 0 or image.max() > 255:
        raise ValueError(f"{path.name} is outside the configured 0–255 camera-count scale; "
                         "do not silently rescale a different acquisition")
    return image


def _binary_mask(path):
    """Accept exact binary levels, including valid uniform empty/full masks."""
    image = _read_2d(path)
    levels = np.unique(image)
    nonzero = levels[levels != 0]
    if len(nonzero) > 1 or (len(nonzero) == 1 and nonzero[0] not in (1, 255, 65535)):
        raise ValueError(f"{path.name} must already be binary (0/1, 0/255, or 0/65535); "
                         f"found {len(levels)} level(s), starting {levels[:6].tolist()}. "
                         "Binarize it externally or choose mask_input.mode='binarize'.")
    return image != 0, [int(v) for v in levels]


def _shape(shape, label):
    values = np.asarray(shape)
    if values.shape != (2,) or not np.all(np.isfinite(values)) or np.any(values != values.astype(int)) \
            or np.any(values < 3):
        raise ValueError(f"{label} must be a (height, width) pair of integers >=3")
    return tuple(int(v) for v in values)


def _check_shape(actual, target, label):
    if target is not None and tuple(actual) != tuple(target):
        raise ValueError(f"{label} shape {tuple(actual)} differs from expected {tuple(target)}. "
                         "References must already share the series pixel grid; no automatic resizing is performed.")


def _source(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path), "sha256": digest.hexdigest()}


def _shift(image, shift_yx):
    if not np.any(shift_yx):
        return image
    return ndi.shift(image.astype(np.float32), shift_yx, order=1,
                     mode="nearest", prefilter=False)


def prepare_reference(cfg, shape_hint=None):
    """Return grain/pore/bed geometry and an explicitly identified display image.

    Required in ``mask_input``:
      * ``mode='binarize'``: ``image_path`` is a separate raw initial image.
      * ``mode='load_binary'``: ``grain_mask_path`` plus ``white_is_grain``.

    Optional ``bed_mask_path`` (nonzero=inside) intersects the configured
    ``bed_polygon_xy`` / ``roi_xyxy``; absent all three, the bed is the full
    frame. Optional ``display_reference_path`` supplies a camera image only
    for figures. In load-binary mode, absent a display image, ``initial`` is
    synthetic dark-grain/bright-pore geometry and ``reference_kind`` says so.
    A nonzero ``initial_shift_yx`` is rejected for external masks: align the
    supplied mask beforehand rather than silently altering its labels.
    """
    options = cfg.get("mask_input", {})
    if not isinstance(options, dict):
        raise ValueError("mask_input must be a dictionary")
    mode = str(options.get("mode", "binarize")).lower()
    if mode not in {"binarize", "load_binary"}:
        raise ValueError("mask_input.mode must be 'binarize' or 'load_binary'")
    folder = Path(cfg.get("input_folder", ".")).expanduser().resolve()
    expected = cfg.get("expected_shape_hw")
    expected = None if expected is None else _shape(expected, "expected_shape_hw")
    hint = None if shape_hint is None else _shape(shape_hint, "shape_hint")
    if expected is not None and hint is not None:
        _check_shape(hint, expected, "Series image")
    target = hint if hint is not None else expected
    shift = np.asarray(cfg.get("initial_shift_yx", (0, 0)), dtype=float)
    if shift.shape != (2,) or not np.all(np.isfinite(shift)):
        raise ValueError("initial_shift_yx must be a finite (dy, dx) pair")

    sources = {}
    if mode == "binarize":
        path = _resolve_path(folder, options.get("image_path", cfg.get("initial_image")), "image_path")
        raw = _camera_image(path)
        _check_shape(raw.shape, target, "Initial image")
        initial = _shift(raw, shift)
        grain, _, bed, info = make_masks(initial, cfg)
        sources["initial_image"] = _source(path)
        kind = "raw_initial"
        info["initial_shift_yx"] = shift.tolist()
        info["binarization_performed"] = True
        info["binary_input_levels"] = None
    else:
        if np.any(shift):
            raise ValueError("initial_shift_yx must be (0,0) for load_binary mode. "
                             "Provide a mask already aligned to the series; external masks are not shifted or resized.")
        white = options.get("white_is_grain", True)
        if not isinstance(white, (bool, np.bool_)):
            raise ValueError("mask_input.white_is_grain must be True or False")
        path = _resolve_path(folder, options.get("grain_mask_path"), "grain_mask_path")
        nonzero, levels = _binary_mask(path)
        _check_shape(nonzero.shape, target, "External grain mask")
        grain = nonzero if white else ~nonzero
        geometry_cfg = dict(cfg)
        if isinstance(cfg.get("grain"), dict):
            geometry_cfg.update(cfg["grain"])
        bed, info = _bed_mask(grain.shape, geometry_cfg)
        sources["grain_mask"] = _source(path)
        info.update({"grain_method": "external binary mask", "binarization_performed": False,
                     "binary_input_levels": levels, "white_is_grain": bool(white),
                     "grain_dilation_px": 0, "initial_shift_yx": [0.0, 0.0], "warnings": []})
        initial = None
        kind = "binary_reference"

    if options.get("bed_mask_path") is not None:
        path = _resolve_path(folder, options["bed_mask_path"], "bed_mask_path")
        supplied_bed, levels = _binary_mask(path)
        _check_shape(supplied_bed.shape, grain.shape, "External bed mask")
        bed &= supplied_bed
        sources["bed_mask"] = _source(path)
        info["bed_mask_levels"] = levels
        info["bed_mask_combination"] = "intersection with configured polygon/ROI, or full frame"
    if not np.any(bed):
        raise ValueError("The combined bed mask contains no pixels")
    outside = int(np.count_nonzero(grain & ~bed))
    grain = grain & bed
    pore = bed & ~grain
    if not np.any(pore):
        raise ValueError("No pore pixels remain in the selected bed; check grain polarity and bed geometry")
    warnings = info.setdefault("warnings", [])
    if outside:
        warnings.append(f"Excluded {outside} grain pixels outside the combined bed mask.")
    if not np.any(grain) and not any("No grains" in warning for warning in warnings):
        warnings.append("No grains present; grain-wall evidence is unavailable.")
    if options.get("display_reference_path") is not None:
        path = _resolve_path(folder, options["display_reference_path"], "display_reference_path")
        initial = _camera_image(path)
        _check_shape(initial.shape, grain.shape, "Display reference")
        sources["display_reference"] = _source(path)
        kind = "display_reference"
        info["display_reference_alignment"] = "provided on the analysis grid; no shift applied"
    elif initial is None:
        initial = np.where(grain, 0, 255).astype(np.uint8)
        warnings.append("The reference display is synthetic binary geometry, not a measured camera image.")
    info.update({"mask_input_mode": mode, "reference_kind": kind,
                 "shape_hw": list(grain.shape), "bed_pixels": int(bed.sum()),
                 "grain_pixels": int(grain.sum()), "pore_pixels": int(pore.sum()),
                 "grain_pixels_excluded_outside_bed": outside,
                 "external_masks_resized_or_cleaned": False})
    return {"initial": initial, "grain": grain, "pore": pore, "bed": bed,
            "geometry_info": info, "reference_sources": sources, "reference_kind": kind}

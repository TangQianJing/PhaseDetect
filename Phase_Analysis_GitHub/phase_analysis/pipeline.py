"""Portable orchestration, file validation, review flags and exports."""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime
from pathlib import Path
import hashlib
import json
import platform
import time

import numpy as np
from scipy import ndimage as ndi
import scipy
import skimage
import tifffile

from .geometry import make_masks, wall_profiles, fluorescence_water
from .segmentation import infer_phases
from .reference import prepare_reference
from .saturation import saturation_metrics, save_csv_atomic, write_saturation_plot


PHASE_CODES = {0: "outside bed", 1: "initial PDMS", 2: "liquid, dark or bright", 3: "air"}
REVIEW_BITS = {1: "weak or fallback region", 2: "watershed-parameter disagreement",
               4: "fluorescence-threshold disagreement", 8: "free-interface band",
               16: "possible shiny grain fringe", 32: "clipping near interface",
               64: "automatic liquid-inclusion boundary"}


def jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(v) for v in value]
    return value


def resolve_image(folder, selector):
    path = Path(str(selector)).expanduser()
    direct = path if path.is_absolute() else folder / path
    if direct.is_file():
        return direct.resolve()
    if path.suffix:
        raise FileNotFoundError(f"Image does not exist: {direct}")
    found = sorted(p for p in folder.iterdir() if p.is_file()
                   and p.suffix.lower() in {".tif", ".tiff"}
                   and p.stem.endswith("_" + str(selector)))
    if len(found) != 1:
        raise ValueError(f"Expected one TIFF ending _{selector} in {folder}; found {len(found)}. "
                         "Use its complete filename in settings.py.")
    return found[0].resolve()


def load_image(path):
    image = tifffile.imread(path)
    if image.ndim != 2 or not image.size or min(image.shape) < 3:
        raise ValueError(f"Expected one 2-D grayscale image: {path}; got shape {image.shape}")
    if not np.all(np.isfinite(image)):
        raise ValueError(f"Nonfinite intensities in {path}")
    if image.min() < 0 or image.max() > 255:
        raise ValueError(f"{path.name} is outside the 0–255 camera-count scale. "
                         "This configuration is calibrated for the supplied 8-bit images; "
                         "do not silently rescale a different acquisition.")
    return image


def shifted(image, shift_yx):
    shift_yx = np.asarray(shift_yx, float)
    if shift_yx.shape != (2,) or not np.all(np.isfinite(shift_yx)):
        raise ValueError("Manual shifts must be finite (dy,dx) pairs")
    if not np.any(shift_yx):
        return image
    return ndi.shift(image.astype(np.float32), shift_yx, order=1, mode="nearest", prefilter=False)


def _hash(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _output_dir(root):
    root = Path(root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    directory = root / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
    directory.mkdir(exist_ok=False)
    return directory


def process_pair(settings, reference, blue_path, red_blue_path, out, progress=None):
    """Process one pair using a prepared, fixed reference; keep arrays for one pair only."""
    cfg = deepcopy(settings)
    started = time.perf_counter()

    if progress is None:
        def progress(message):
            print(f"[{time.perf_counter() - started:6.1f}s] {message}", flush=True)
    progress("Reading B and R&B; reusing the fixed reference")
    paths = {"blue_image": Path(blue_path), "red_blue_image": Path(red_blue_path)}
    B, RB = load_image(paths["blue_image"]), load_image(paths["red_blue_image"])
    initial = reference["initial"]
    grain, pore, bed = (reference[k] for k in ("grain", "pore", "bed"))
    geometry_info = reference["geometry_info"]
    if not (initial.shape == B.shape == RB.shape == grain.shape):
        raise ValueError(f"Grid mismatch: reference={grain.shape}, B={B.shape}, R&B={RB.shape}")
    B = shifted(B, cfg.get("blue_shift_yx", (0, 0)))
    alpha = float(cfg.get("blue_scale", 1.1))
    if not np.isfinite(alpha) or alpha < 0:
        raise ValueError("blue_scale must be finite and nonnegative")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    progress(f"Output folder: {out}")

    progress("Step 1: using prepared grain and pore masks")
    if not pore.any():
        raise ValueError("The selected bed contains no pore pixels; inspect grain settings and ROI")
    progress("Step 2: signed subtraction and initial-wall normal profiles")
    difference = RB.astype(np.float32) - alpha * B.astype(np.float32)
    profiles = wall_profiles(B, RB, grain, pore, cfg)
    for warning in geometry_info.get("warnings", []) + profiles.get("warnings", []):
        progress("Warning: " + warning)

    progress("Step 3: fix positive fluorescence as liquid; withhold local shiny-wall peaks")
    direct, shiny, smoothed_B = fluorescence_water(B, grain, pore, profiles, cfg)
    progress(f"Fixed-water pixels: {int(direct.sum()):,}; possible shiny fringes: {int(shiny.sum()):,}")
    progress("Step 4: infer remaining pore phases while preserving fixed water")
    result = infer_phases(initial, B, RB, grain, pore, profiles, direct, cfg, progress=progress)

    sensitivity = {}
    main_threshold = float(cfg.get("fluorescence", {}).get("threshold", 55))
    for threshold in cfg.get("fluorescence_sensitivity", []):
        threshold = float(threshold)
        if threshold == main_threshold:
            continue
        progress(f"Fluorescence sensitivity: cutoff {threshold:g} camera counts")
        fixed, exception, _ = fluorescence_water(B, grain, pore, profiles, cfg, threshold=threshold)
        local_cfg = deepcopy(cfg)
        local_cfg.setdefault("segmentation", {})["watershed_sensitivity"] = False
        variant = infer_phases(initial, B, RB, grain, pore, profiles, fixed, local_cfg)
        sensitivity[threshold] = {"phase": variant["phase"], "direct_water": fixed,
                                  "shiny_exception": exception,
                                  "liquid_fraction": float((variant["phase"][pore] == 2).mean()),
                                  "air_fraction": float((variant["phase"][pore] == 3).mean())}
        if not np.all(variant["phase"][fixed] == 2):
            raise AssertionError("Sensitivity run changed a fixed water pixel")
    phase = result["phase"]
    progress("Step 5: optical-boundary, parameter and evidence review flags")
    flags = np.zeros(phase.shape, np.uint16)
    flags[result["weak"] & pore] |= 1
    flags[pore & ~result["stable"]] |= 2
    fluorescence_disagree = np.zeros(phase.shape, bool)
    for variant in sensitivity.values():
        fluorescence_disagree |= pore & (variant["phase"] != phase)
    flags[fluorescence_disagree] |= 4
    width = int(cfg.get("review_boundary_width_px", 2))
    if width < 0:
        raise ValueError("review_boundary_width_px cannot be negative")
    free = ((ndi.binary_dilation(phase == 2) & (phase == 3))
            | (ndi.binary_dilation(phase == 3) & (phase == 2)))
    band = (ndi.binary_dilation(free, iterations=width) if width else free) & pore
    flags[band] |= 8
    flags[shiny] |= 16
    clipped = pore & ((B >= 255) | (RB >= 255))
    flags[clipped & band] |= 32
    inclusion = result["inclusion_core"]
    inc_border = ndi.binary_dilation(inclusion) & ~ndi.binary_erosion(inclusion) & pore
    flags[inc_border] |= 64
    # Strong broad fluorescence resolves old dark-phase ambiguity, but not an
    # optical boundary or a threshold-dependent assignment.
    strong = direct & ~fluorescence_disagree & ~band & ~shiny
    flags[strong] &= np.uint16(0xffff ^ (1 | 2))
    review_phase = phase.copy()
    review_phase[flags > 0] = 4
    support = result["support"].copy()
    support[strong] = np.maximum(support[strong], .95)
    support[flags > 0] = np.minimum(support[flags > 0], .45)

    invariants = {"fixed_water_preserved": bool(np.all(phase[direct] == 2)),
                  "grains_preserved": bool(np.all(phase[grain] == 1)),
                  "outside_bed_excluded": bool(np.all(phase[~bed] == 0)),
                  "pore_partition_complete": bool(np.all(np.isin(phase[pore], (2, 3)))),
                  "no_fluorescence_on_grain": bool(not np.any(direct & grain))}
    if not all(invariants.values()):
        raise AssertionError(f"Physical/mask consistency check failed: {invariants}")
    versions = {"python": platform.python_version(), "numpy": np.__version__,
                "scipy": scipy.__version__, "scikit_image": skimage.__version__,
                "tifffile": tifffile.__version__}
    saturation = saturation_metrics(phase, pore, flags, B=B, RB=RB)
    summary = {"shape_hw": list(phase.shape), "pore_pixels": int(pore.sum()),
               "liquid_fraction_of_pore": float((phase[pore] == 2).mean()),
               "air_fraction_of_pore": float((phase[pore] == 3).mean()),
               "review_fraction_of_pore": float((flags[pore] > 0).mean()),
               "direct_fluorescence_water_pixels": int(direct.sum()),
               "shiny_exception_pixels": int(shiny.sum()),
               "blue_scale": alpha, "fluorescence_cutoff_counts": main_threshold,
               "fluorescence_sensitivity": [
                   {"threshold_counts": float(t),
                    "liquid_fraction_of_pore": v["liquid_fraction"],
                    "air_fraction_of_pore": v["air_fraction"],
                    "direct_water_fraction_of_pore": float(v["direct_water"][pore].mean())}
                   for t, v in sensitivity.items()],
               "watershed_sensitivity": result["records"].get("sensitivity", []),
               "grain_geometry": geometry_info, "checks": invariants,
               "reference_kind": reference["reference_kind"],
               "saturation": saturation,
               "runtime_versions": versions,
               "interpretation": "Projected pore-pixel estimates; no physical volume or accuracy calibration.",
               "annotations_used_as_seeds": False,
               "water_identified_before_remaining_phase_inference": True,
               "stage_order": ["initial geometry", "signed subtraction and wall profiles",
                               "fixed positive-fluorescence water", "remaining phase inference",
                               "sensitivity and review", "exports and illustrated guide"]}
    progress("Step 6: write full-resolution masks and provenance")

    def save(name, array, description):
        level = cfg.get("export_level", "full")
        if level not in {"full", "minimal", "summary_only"}:
            raise ValueError("export_level must be full, minimal, or summary_only")
        if level == "summary_only":
            return
        if level == "minimal" and name not in {"06_phase_distribution.tiff", "07_phase_with_review.tiff", "08_review_flags.tiff"}:
            return
        tifffile.imwrite(out / name, array, compression="deflate", photometric="minisblack",
                         metadata=None, description=description)

    save("01_initial_grain_mask.tiff", grain.astype(np.uint8) * 255, "255initialPDMS;0other;no dilation")
    save("02_analysed_pore_mask.tiff", pore.astype(np.uint8) * 255, "255analysedinitialpore;0grain/outside")
    save("03_RB_minus_scaled_B_float32.tiff", difference, "Signedfloat32RB-alpha*B;alpha in run_settings.json")
    save("04_fixed_fluorescence_water.tiff", direct.astype(np.uint8) * 255, "255fixedpositivefluorescencewater;0doesnotmean air")
    save("05_possible_shiny_grain_fringe.tiff", shiny.astype(np.uint8) * 255, "255localshinywallpeakexception;notprovenair")
    save("06_phase_distribution.tiff", phase, "0outsidebed;1PDMS;2liquiddarkorbright;3air;nominalproposal")
    save("07_phase_with_review.tiff", review_phase, "0outside;1PDMS;2liquid;3air;4review/unresolved")
    save("08_review_flags.tiff", flags, "uint16bitmask;see review_flag_codes.json")
    save("09_support_index_float32.tiff", support.astype(np.float32), "0:1heuristicsupport;notprobability")
    save("10_optical_region_labels.tiff", result["region_labels"].astype(np.int32), "AutomaticwatershedregionIDs;notphysicalphaseIDs")
    save("11_liquid_inclusion_cores.tiff", inclusion.astype(np.uint8) * 255, "255autobrightliquidcoreproposal;notfullphysicaldroplet")
    save("12_liquid_mask.tiff", (phase == 2).astype(np.uint8) * 255, "255nominalliquidphase;see reviewmask")
    save("13_air_mask.tiff", (phase == 3).astype(np.uint8) * 255, "255nominalairphase;see reviewmask")
    for threshold, variant in sensitivity.items():
        save(f"sensitivity_B_cutoff_{threshold:g}.tiff", variant["phase"], "NominalphasewithalternateBcutoff;0outside1grain2liquid3air")
    files = {key: {"path": str(path), "sha256": _hash(path) if cfg.get("hash_pair_sources", True) else None}
             for key, path in paths.items()}
    files["reference_sources"] = reference["reference_sources"]
    for name, value in [("run_settings.json", cfg), ("summary.json", summary),
                        ("source_images.json", files), ("phase_codes.json", PHASE_CODES),
                        ("review_flag_codes.json", REVIEW_BITS),
                        ("region_evidence.json", result["records"]),
                        ("algorithm_parameters.json", result["parameters"])]:
        (out / name).write_text(json.dumps(jsonable(value), indent=2), encoding="utf-8")
    (out / "saturation.json").write_text(json.dumps(jsonable(saturation), indent=2), encoding="utf-8")
    if cfg.get("save_intermediate_arrays", False):
        np.savez_compressed(out / "analysis_arrays.npz", initial=initial, B=B, RB=RB,
                            grain=grain, pore=pore, phase=phase, direct_water=direct,
                            shiny_exception=shiny, region_labels=result["region_labels"],
                            support=support, review_flags=flags)
    data = {"initial": initial, "B": B, "RB": RB, "grain": grain, "pore": pore,
            "bed": bed, "difference": difference, "direct_water": direct,
            "shiny_exception": shiny, "smoothed_B": smoothed_B, "phase": phase,
            "review_phase": review_phase, "review_flags": flags, "support": support,
            "region_labels": result["region_labels"], "profiles": profiles,
            "threshold_sensitivity": sensitivity, "summary": summary,
            "reference_kind": reference["reference_kind"],
            "watershed_sensitivity": result["records"].get("sensitivity", []),
            "run_directory": out}
    if cfg.get("make_report", True):
        progress("Drawing example figures and the detailed illustrated guide")
        from .reporting import write_report
        report = write_report(data, cfg, out)
        data["report"] = report
    summary["elapsed_seconds"] = round(time.perf_counter() - started, 2)
    (out / "summary.json").write_text(json.dumps(jsonable(summary), indent=2), encoding="utf-8")
    progress(f"Finished. Liquid {100 * summary['liquid_fraction_of_pore']:.1f}% / "
             f"air {100 * summary['air_fraction_of_pore']:.1f}% of pore pixels; "
             f"review flagged {100 * summary['review_fraction_of_pore']:.1f}%.")
    if cfg.get("make_report", True):
        progress(f"Open {out / 'GUIDE.html'} for the illustrated result")
    if cfg.get("show_figures", True):
        import matplotlib
        import matplotlib.pyplot as plt
        overview = out / "figures" / "01_input_images.png"
        # Report figures are saved regardless of Spyder's plotting backend.
        if overview.exists() and matplotlib.get_backend().lower() != "agg":
            fig, ax = plt.subplots(figsize=(12, 5))
            ax.imshow(plt.imread(overview)); ax.axis("off")
            fig.suptitle("Analysis complete — open GUIDE.html for all figures")
            plt.show(block=False)
        elif overview.exists():
            progress("Noninteractive plotting backend: all figures are saved in the report folder.")
    return data


def run_single(settings):
    """Backwards-compatible one-pair mode with separate reference input."""
    cfg = deepcopy(settings)
    folder = Path(cfg["input_folder"]).expanduser().resolve()
    blue = resolve_image(folder, cfg["blue_image"])
    rb = resolve_image(folder, cfg["red_blue_image"])
    with tifffile.TiffFile(blue) as tf:
        shape = tf.series[0].shape
    reference = prepare_reference(cfg, shape_hint=shape)
    out = _output_dir(cfg["output_folder"])
    data = process_pair(cfg, reference, blue, rb, out)
    row = {"pair_number": 0, "status": "completed", "blue_path": str(blue),
           "red_blue_path": str(rb), **data["summary"]["saturation"]}
    save_csv_atomic([row], out / "saturation.csv")
    write_saturation_plot([row], out / "saturation_vs_pair.png")
    return data


def run(settings):
    mode = str(settings.get("run_mode", "single")).lower()
    if mode == "single":
        return run_single(settings)
    if mode == "series":
        from .batch import run_series
        return run_series(settings)
    raise ValueError("run_mode must be 'single' or 'series'")

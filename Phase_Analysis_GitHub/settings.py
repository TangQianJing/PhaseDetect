"""Edit the input, mask and series blocks, then press F5 in RUN_ANALYSIS.py.

All pixel thresholds are on the original 8-bit camera-count scale, NOT the
displayed image brightness. The supplied paths and bed polygon fit this image
triplet. For a different field of view, update the expected shape and polygon.
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent

SETTINGS = {
    "run_mode": "series",  # "series" or "single" for one explicitly named pair.
    "input_folder": HERE / "input_images",  # Folder containing the alternating series.
    "output_folder": HERE / "results",  # New timestamped folder for each fresh run.

    # Independent geometry input. Absolute paths can point outside input_folder.
    "mask_input": {
        "mode": "binarize",  # "binarize" raw reference, OR "load_binary" exact mask.
        "image_path": HERE / "input_images" / "Basler_acA4112-20um__40402042__20261023_163522172_0001.tiff",
        "grain_mask_path": None,  # Required only for load_binary; TIFF or grayscale PNG.
        "white_is_grain": True,   # False means black pixels are grains.
        "bed_mask_path": None,   # Optional binary nonzero=inside; intersects polygon/ROI.
        "display_reference_path": None,  # Optional raw image for figures in mask mode.
    },

    # Inclusive source indices, not number of pairs. Demo uses the provided pair.
    # For future even-B series: start_index=0, end_index=1999, blue_parity="even".
    # For odd-B series change only blue_parity="odd"; channel order is never guessed.
    "series": {
        "start_index": 6544,
        "end_index": 6545,
        "blue_parity": "odd",  # ABSOLUTE index parity: even 0,2,... OR odd 1,3,...
        "file_pattern": None,  # Optional acquisition glob, e.g. "Basler_*20261023*.tiff".
        "index_regex": r"(?P<index>\d+)$",  # Final numeric filename suffix, ignoring padding.
        "unpaired_policy": "error",  # "drop_last" explicitly omits the final lone frame.
        "missing_pair_policy": "error",  # "skip" omits WHOLE pairs; never shifts pairing.
        "pair_interval_seconds": None,  # Actual time per pair, if known; no guessed timing.
    },
    "batch": {
        "report_policy": "first",  # "first", "all", "none": illustrated reports per pair.
        "sensitivity_policy": "none",  # "all" adds slow parameter sweeps to EVERY pair.
        "export_level": "minimal",  # "minimal": phase/review/flags; "full"; "summary_only".
        "pair_error_policy": "stop",  # "continue" saves the failed row and moves on.
        "plot_every": 10,  # CSV saved EVERY pair; plot refreshed this often and at end.
        "preflight_shapes": True,  # Check all TIFF headers before the first analysis.
        "resume_from": None,  # Existing run folder to resume; same science/range/files required.
        "retain_last_pair": True,  # Keep only last pair's arrays in Spyder's RESULT.
    },

    # Used only in run_mode="single". Unique suffix, full filename or absolute path.
    "red_blue_image": "6544",
    "blue_image": "6545",
    "expected_shape_hw": (3000, 4096),  # Set None only after checking a new field.
    "show_figures": True,
    "make_report": True,
    "save_intermediate_arrays": False,  # Optional compressed NumPy archive.
    "hash_pair_sources": True,  # SHA256 of each analysed frame in per-pair provenance.

    # Analysis area: source-pixel coordinates [x,y], x right, y down.
    # These slightly sloping bed edges were reviewed for the supplied experiment.
    # For a rectangular crop instead, set bed_polygon_xy=None and roi_xyxy below.
    "bed_polygon_xy": [[55, 576], [4030, 590], [4030, 2580], [55, 2570]],
    "roi_xyxy": None,  # e.g. [100,600,3900,2500], upper limits exclusive.

    "blue_scale": 1.1,  # Signed difference: R&B - blue_scale * B.
    "initial_shift_yx": (0.0, 0.0),  # Optional manual registration; no shift now.
    "blue_shift_yx": (0.0, 0.0),

    "grain": {
        "window": 101, "k": 0.25, "r": 0.5, "p": 2.0, "q": 10.0,
        "sigma": 0.8, "min_area": 25, "intensity_scale": 255.0,
    },
    "wall": {
        "normal_sigma": 1.1, "channel_sigma": 0.9,
        "bed_clearance": 10.0, "profile_chunk_size": 200000,
        "lit_rb": 14.0, "lit_d": 12.0, "lit_d_min_rb": 6.0,
        "lit_b": 10.0, "lit_b_min_rb": 5.0,
        "unlit_rb": 6.0, "unlit_b": 5.0, "unlit_d": 6.0,
    },
    "fluorescence": {
        "threshold": 55.0, "sigma": 0.9,
        "shiny_distance": 4.5, "shiny_contrast": 10.0, "outer_margin": 10.0,
    },
    "segmentation": {
        "watershed_sigma_px": 1.0, "watershed_h": 8.0,
        "watershed_sensitivity": [5.0, 8.0, 12.0],
        "air_wall_lit_fraction": 0.27,
        "air_wall_min_lit_samples": 12, "air_wall_min_valid_samples": 20,
        "air_wall_max_unlit_fraction": 0.65,
        "broad_water_region_fraction": 0.60, "broad_water_min_pixels": 20,
        "use_fringe_graph": True, "detect_liquid_inclusions": True,
    },
    # Each threshold reruns phase assignment with its own fixed-water mask.
    # Set [] for a faster run; these variants are sensitivity, not error bars.
    "fluorescence_sensitivity": [45.0, 65.0],
    "review_boundary_width_px": 2,
    "gradient_sigma_px": 0.9,

    # Example boxes only control report figures, never phase classification.
    "examples": [
        {"name": "Broad fluorescent water", "xyxy": [760, 1480, 1040, 1760]},
        {"name": "Dark water / air and liquid inclusions", "xyxy": [1155, 1025, 1335, 1205]},
        {"name": "Dry-grain fringes", "xyxy": [800, 650, 1100, 950]},
        {"name": "Dark pore region", "xyxy": [2800, 1180, 3040, 1420]},
    ],
}

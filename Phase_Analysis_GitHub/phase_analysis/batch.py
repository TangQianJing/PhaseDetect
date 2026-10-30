"""Sequential alternating-frame analysis with durable, resumable per-pair rows."""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime
from pathlib import Path
import csv
import gc
import hashlib
import json
import os
import tempfile
import time

import numpy as np
import tifffile

from .pipeline import process_pair, jsonable, _output_dir
from .reference import prepare_reference
from .series_io import discover_indexed_images, build_pair_plan
from .saturation import save_csv_atomic, write_saturation_plot
from . import __version__


def _write_json(path, value):
    """Replace metadata atomically, including on an interrupted batch."""
    path = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", delete=False,
                                         dir=path.parent, suffix=".tmp") as stream:
            temporary = stream.name
            json.dump(jsonable(value), stream, indent=2, allow_nan=False)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            Path(temporary).unlink(missing_ok=True)


def _choice(options, key, default, values):
    value = options.get(key, default)
    if value not in values:
        raise ValueError(f"batch.{key} must be one of {sorted(values)}")
    return value


def _shape(path):
    with tifffile.TiffFile(path) as source:
        shape = tuple(source.series[0].shape)
    if len(shape) != 2 or min(shape) < 3:
        raise ValueError(f"Expected a single 2-D grayscale TIFF: {path}; shape={shape}")
    return shape


def _stat(path):
    stat = Path(path).stat()
    return {"path": str(path), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _fingerprint(cfg, plan, reference, sensitivity_policy):
    # Runtime/report/export choices can change on resume; science cannot.
    runtime = {"batch", "input_folder", "output_folder", "run_mode", "show_figures",
               "make_report", "save_intermediate_arrays", "export_level", "hash_pair_sources",
               "examples", "initial_image", "blue_image", "red_blue_image", "mask_input", "series"}
    science = {k: v for k, v in cfg.items() if k not in runtime}
    science["sensitivity_policy"] = sensitivity_policy
    science["mask_mode"] = cfg.get("mask_input", {}).get("mode", "binarize")
    science["white_is_grain"] = cfg.get("mask_input", {}).get("white_is_grain", True)
    frame_stats = {str(p[key]): _stat(p[key]) for p in plan["pairs"]
                   for key in ("blue_path", "red_blue_path")}
    payload = {"schema": 1, "science_settings": science, "pair_plan": plan,
               "reference_sources": reference["reference_sources"],
               "reference_kind": reference["reference_kind"], "frame_stats": frame_stats,
               "algorithm_version": __version__,
               "algorithm_sources_sha256": {
                   name: hashlib.sha256((Path(__file__).parent / name).read_bytes()).hexdigest()
                   for name in ("geometry.py", "segmentation.py", "reference.py", "pipeline.py",
                                "batch.py", "series_io.py", "saturation.py")}}
    serial = json.dumps(jsonable(payload), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(serial.encode("utf-8")).hexdigest(), payload


def _pair_base(pair, out):
    row = {key: pair.get(key) for key in ("pair_number", "first_index", "second_index",
                                         "blue_index", "red_blue_index", "time_seconds",
                                         "blue_path", "red_blue_path")}
    row["pair_directory"] = str(out)
    return row


def run_series(settings):
    """Prepare geometry once, then save saturation immediately after every pair.

    Memory use is bounded by the reference and one pair. The returned last_pair
    is optional; all earlier arrays live only in their exported files.
    """
    cfg = deepcopy(settings)
    options = cfg.get("batch", {})
    report_policy = _choice(options, "report_policy", "first", {"none", "first", "all"})
    sensitivity_policy = _choice(options, "sensitivity_policy", "none", {"none", "first", "all"})
    export_level = _choice(options, "export_level", "minimal", {"full", "minimal", "summary_only"})
    error_policy = _choice(options, "pair_error_policy", "stop", {"stop", "continue"})
    plot_every = options.get("plot_every", 10)
    if isinstance(plot_every, bool) or not isinstance(plot_every, int) or plot_every < 1:
        raise ValueError("batch.plot_every must be a positive integer")
    started = time.perf_counter()

    def progress(message):
        print(f"[batch {time.perf_counter() - started:7.1f}s] {message}", flush=True)

    series = cfg["series"]
    mapping = discover_indexed_images(cfg["input_folder"], series.get("file_pattern"),
                                     series.get("index_regex", r"(?P<index>\d+)$"))
    plan = build_pair_plan(cfg, mapping)
    target_shape = _shape(plan["pairs"][0]["blue_path"])
    if options.get("preflight_shapes", True):
        progress(f"Checking image headers for {len(plan['pairs'])} pairs")
        for pair in plan["pairs"]:
            for key in ("blue_path", "red_blue_path"):
                if _shape(pair[key]) != target_shape:
                    raise ValueError(f"Grid mismatch before processing: {pair[key]}")
    progress("Preparing the separate grain reference once")
    reference = prepare_reference(cfg, shape_hint=target_shape)
    fingerprint, provenance = _fingerprint(cfg, plan, reference, sensitivity_policy)
    resume = options.get("resume_from")
    if resume:
        out = Path(resume).expanduser().resolve()
        saved = json.loads((out / "batch_manifest.json").read_text(encoding="utf-8"))
        if saved.get("fingerprint") != fingerprint:
            raise ValueError("Resume rejected: science settings, references, series plan or source-file "
                             "size/modification time changed. Start a fresh batch to avoid mixing results.")
        with (out / "saturation.csv").open(newline="", encoding="utf-8") as stream:
            rows = {int(row["pair_number"]): row for row in csv.DictReader(stream)}
        progress(f"Resuming {out}")
    else:
        out = _output_dir(cfg["output_folder"])
        rows = {}
        _write_json(out / "batch_manifest.json", {"fingerprint": fingerprint, **provenance})
        _write_json(out / "batch_settings.json", cfg)
        _write_json(out / "pair_plan.json", plan)
        pair_rows = [dict(p, planned_status="ready") for p in plan["pairs"]]
        pair_rows += [dict(p, planned_status="missing_pair_skipped") for p in plan["skipped_pairs"]]
        save_csv_atomic(sorted(pair_rows, key=lambda p: p["pair_number"]), out / "pair_plan.csv")
        ref_dir = out / "reference"
        ref_dir.mkdir()
        for key in ("grain", "pore", "bed"):
            tifffile.imwrite(ref_dir / f"{key}_mask.tiff", reference[key].astype(np.uint8) * 255,
                             compression="deflate", photometric="minisblack", metadata=None)
        _write_json(ref_dir / "reference_info.json", {
            "reference_kind": reference["reference_kind"],
            "geometry_info": reference["geometry_info"],
            "sources": reference["reference_sources"],
            "pore_pixel_count": int(reference["pore"].sum())})
    for skipped in plan["skipped_pairs"]:
        number = skipped["pair_number"]
        rows[number] = {**skipped, "status": "missing_pair_skipped",
                        "time_seconds": None if plan["pair_interval_seconds"] is None
                        else number * plan["pair_interval_seconds"],
                        "error": "Missing frame index/indices " + str(skipped["missing_indices"])}
    last_pair = None
    attempts = 0
    xaxis = "time_seconds" if plan["pair_interval_seconds"] is not None else "pair_number"

    def snapshot(status, plot=False):
        ordered = [rows[k] for k in sorted(rows)]
        # Header-only CSV is valid before the first pair finishes.
        save_csv_atomic(ordered, out / "saturation.csv",
                        fieldnames=None if ordered else ["pair_number", "status"])
        counts = {s: sum(row.get("status") == s for row in ordered)
                  for s in ("completed", "failed", "interrupted", "missing_pair_skipped")}
        _write_json(out / "batch_summary.json", {
            "status": status, "updated_local": datetime.now().isoformat(timespec="seconds"),
            "run_directory": out, "fingerprint": fingerprint, "counts": counts,
            "available_pairs": len(plan["pairs"]), "unpaired_indices": plan["unpaired_indices"],
            "pore_pixel_count": int(reference["pore"].sum()),
            "elapsed_this_invocation_seconds": round(time.perf_counter() - started, 2),
            "saturation_definition": "2D projected liquid/air pixel area divided by fixed pore pixel area",
            "review_range": "All flagged pixels may change phase; unflagged labels held fixed. Not a CI.",
            "source_check": "Reference SHA256; series file path, size and mtime_ns; per-pair hashes if enabled",
            "sensitivity_policy": sensitivity_policy})
        if plot:
            try:
                write_saturation_plot(ordered, out / "saturation_vs_pair.png", xaxis=xaxis)
            except Exception as exc:
                progress(f"Plot could not be updated; saturation CSV is saved: {exc}")

    snapshot("running")
    progress(f"Output: {out}; available pairs={len(plan['pairs'])}, "
             f"missing pairs skipped={len(plan['skipped_pairs'])}, dropped frames={plan['unpaired_indices']}")
    for position, pair in enumerate(plan["pairs"]):
        number = pair["pair_number"]
        pair_dir = out / "pairs" / f"pair_{number:06d}_B{pair['blue_index']}_RB{pair['red_blue_index']}"
        if rows.get(number, {}).get("status") == "completed":
            marker = pair_dir / "pair_complete.json"
            if marker.exists() and json.loads(marker.read_text(encoding="utf-8")).get("fingerprint") == fingerprint:
                progress(f"Pair {number}: already completed; skipping")
                continue
            progress(f"Pair {number}: completion marker missing; recalculating")
        # Drop previous arrays BEFORE allocating the next pair (important at full resolution).
        last_pair = None
        gc.collect()
        local = deepcopy(cfg)
        local["make_report"] = report_policy == "all" or (report_policy == "first" and position == 0)
        local["show_figures"] = False
        local["export_level"] = export_level
        check_sensitivity = sensitivity_policy == "all" or (sensitivity_policy == "first" and position == 0)
        if not check_sensitivity:
            local["fluorescence_sensitivity"] = []
            local.setdefault("segmentation", {})["watershed_sensitivity"] = False
        check_sensitivity = bool(local.get("fluorescence_sensitivity")
                                 or local.get("segmentation", {}).get("watershed_sensitivity"))
        row = {**_pair_base(pair, pair_dir), "status": "running", "sensitivity_checked": check_sensitivity,
               "started_local": datetime.now().isoformat(timespec="seconds"), "error": ""}
        tick = time.perf_counter()
        progress(f"Pair {number}: B={pair['blue_index']}, R&B={pair['red_blue_index']}")
        try:
            for key in ("blue_path", "red_blue_path"):
                if _stat(pair[key]) != provenance["frame_stats"][str(pair[key])]:
                    raise ValueError(f"Source changed after planning: {pair[key]}")
            data = process_pair(local, reference, pair["blue_path"], pair["red_blue_path"], pair_dir,
                                progress=lambda msg: progress(f"Pair {number}: {msg}"))
            metrics = data["summary"]["saturation"]
            if metrics["pore_pixel_count"] != int(reference["pore"].sum()):
                raise AssertionError("Pore denominator changed between pairs")
            row.update(metrics)
            row["status"] = "completed"
            row["elapsed_seconds"] = round(time.perf_counter() - tick, 3)
            _write_json(pair_dir / "pair_complete.json", {"fingerprint": fingerprint, "row": row})
            if options.get("retain_last_pair", True):
                last_pair = data
            del data
        except KeyboardInterrupt:
            row.update(status="interrupted", elapsed_seconds=round(time.perf_counter() - tick, 3),
                       error="KeyboardInterrupt; this pair will be retried on resume")
            rows[number] = row
            snapshot("interrupted", plot=True)
            progress("Interrupted. Earlier pairs and the saturation table are saved.")
            raise
        except Exception as exc:
            row.update(status="failed", elapsed_seconds=round(time.perf_counter() - tick, 3),
                       error=f"{type(exc).__name__}: {exc}")
            rows[number] = row
            snapshot("stopped_on_error" if error_policy == "stop" else "running", plot=True)
            progress(f"Pair {number} failed; saved its error in saturation.csv: {exc}")
            if error_policy == "stop":
                raise
        rows[number] = row
        attempts += 1
        snapshot("running", plot=attempts % plot_every == 0)
    status = "complete_with_errors" if any(r.get("status") == "failed" for r in rows.values()) else "complete"
    snapshot(status, plot=True)
    progress(f"Finished: {out / 'saturation.csv'}")
    return {"run_directory": out, "saturation_rows": [rows[k] for k in sorted(rows)],
            "pair_plan": plan, "last_pair": last_pair, "status": status}

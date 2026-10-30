"""Discover indexed TIFFs and plan alternating illumination pairs.

This module does not read image pixels. Illumination is supplied explicitly
through ``blue_parity`` and is never inferred from image brightness. Pairing is
anchored at ``start_index``; a missing image never shifts the following pairs.
Only the standard library is required, so planning can run before the heavier
image-analysis imports.
"""
from __future__ import annotations

from collections.abc import Mapping
from fnmatch import fnmatchcase
from numbers import Integral, Real
from pathlib import Path
import math
import os
import re


def _integer(value, name: str) -> int:
    """Accept genuine finite integers; reject bool, floats and numeric strings."""
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer, not {value!r}.")
    return int(value)


def _choice(value, name: str, choices: tuple[str, ...]) -> str:
    if not isinstance(value, str) or value.strip().lower() not in choices:
        raise ValueError(f"{name} must be one of {choices}; got {value!r}.")
    return value.strip().lower()


def discover_indexed_images(
    folder,
    pattern: str | None = None,
    index_regex: str = r"(?P<index>\d+)$",
) -> dict[int, Path]:
    """Return ``{frame_index: absolute_path}`` for TIFFs in one folder.

    Both .tif and .tiff are accepted, ignoring extension case. The default
    expression extracts the LAST numeric suffix from each filename stem, so
    camera serial numbers and acquisition timestamps earlier in the name are
    ignored. Leading zeroes are allowed. Nonmatching files are ignored.

    ``pattern`` is an optional, case-insensitive filename glob such as
    ``"Basler_acA4112*20261023*.tiff"``. It selects an acquisition before index
    extraction. Discovery is not recursive. A custom expression must contain
    the named capture group ``(?P<index>...)``. Duplicate numeric indices are
    rejected; use a narrower folder or glob rather than silently choosing one.
    """
    root = Path(folder).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Image folder does not exist: {root}")
    if pattern is not None and (not isinstance(pattern, str) or not pattern):
        raise ValueError("series.file_pattern must be None or a nonempty filename glob.")
    if not isinstance(index_regex, str):
        raise ValueError("series.index_regex must be a regular-expression string.")
    try:
        expression = re.compile(index_regex)
    except re.error as exc:
        raise ValueError(f"Invalid series.index_regex: {exc}") from exc
    if "index" not in expression.groupindex:
        raise ValueError("series.index_regex needs a named group: (?P<index>...).")

    found: dict[int, Path] = {}
    for path in sorted(root.iterdir(), key=lambda p: p.name.casefold()):
        if not path.is_file() or path.suffix.casefold() not in (".tif", ".tiff"):
            continue
        if pattern is not None and not fnmatchcase(path.name.casefold(), pattern.casefold()):
            continue
        match = expression.search(path.stem)
        if match is None:
            continue
        raw_index = match.group("index")
        if raw_index is None or re.fullmatch(r"[0-9]+", raw_index) is None:
            raise ValueError(f"The index group must contain nonnegative digits: {path.name}")
        index = int(raw_index)
        if index in found:
            raise ValueError(
                f"Duplicate frame index {index}: '{found[index].name}' and '{path.name}'. "
                "Choose one acquisition with series.file_pattern or a separate folder."
            )
        found[index] = path.resolve()
    if not found:
        selector = f" matching {pattern!r}" if pattern is not None else ""
        raise FileNotFoundError(f"No indexed .tif/.tiff images found in {root}{selector}.")
    return dict(sorted(found.items()))


def build_pair_plan(cfg: Mapping, mapping: Mapping[int, Path]) -> dict:
    """Validate every requested pair before image processing.

    Accept either the full SETTINGS mapping or its ``series`` mapping. Required
    settings are ``start_index``, ``end_index`` (both inclusive), and
    ``blue_parity`` (``'even'`` or ``'odd'`` in ABSOLUTE frame numbering).
    Optional policies default to ``'error'``: ``unpaired_policy`` may be
    ``'drop_last'``; ``missing_pair_policy`` may be ``'skip'``. A skip omits the
    whole requested pair without re-pairing other frames.

    Pairs are (start, start+1), (start+2, start+3), ... . ``pair_number`` is the
    zero-based position in that requested sequence; skipped positions remain
    gaps. Optional positive ``pair_interval_seconds`` gives elapsed time from
    the FIRST REQUESTED pair, preserving the same gaps.

    The result contains ``pairs`` plus explicit ``skipped_pairs`` and
    ``unpaired_indices`` for saving with the run. ``time_seconds`` is None when
    no acquisition interval has been supplied. Paths are pathlib.Path values.
    """
    if not isinstance(cfg, Mapping):
        raise ValueError("Settings must be a mapping.")
    series = cfg.get("series", cfg)
    if not isinstance(series, Mapping):
        raise ValueError("settings['series'] must be a mapping.")
    for required in ("start_index", "end_index", "blue_parity"):
        if required not in series:
            raise ValueError(f"Missing series setting: {required}")
    start = _integer(series["start_index"], "series.start_index")
    end = _integer(series["end_index"], "series.end_index")
    if start < 0 or end < 0:
        raise ValueError("Frame indices must be nonnegative.")
    if end < start:
        raise ValueError("series.end_index must be >= series.start_index (inclusive bounds).")
    blue_parity = _choice(series["blue_parity"], "series.blue_parity", ("even", "odd"))
    unpaired_policy = _choice(
        series.get("unpaired_policy", "error"), "series.unpaired_policy", ("error", "drop_last")
    )
    missing_policy = _choice(
        series.get("missing_pair_policy", "error"), "series.missing_pair_policy", ("error", "skip")
    )
    interval = series.get("pair_interval_seconds")
    if interval is not None:
        if (isinstance(interval, bool) or not isinstance(interval, Real)
                or not math.isfinite(float(interval)) or interval <= 0):
            raise ValueError("series.pair_interval_seconds must be a positive finite number or None.")
        interval = float(interval)

    if not isinstance(mapping, Mapping):
        raise ValueError("Image index mapping must be a mapping of integer indices to paths.")
    images: dict[int, Path] = {}
    for key, value in mapping.items():
        index = _integer(key, "Image mapping index")
        if index < 0:
            raise ValueError("Image mapping indices must be nonnegative.")
        try:
            images[index] = Path(value).expanduser().resolve()
        except TypeError as exc:
            raise ValueError(f"Invalid image path for index {index}: {value!r}") from exc

    count = end - start + 1
    unpaired = [end] if count % 2 else []
    if unpaired and unpaired_policy == "error":
        raise ValueError(
            f"Inclusive range {start}..{end} contains {count} frames; frame {end} is unpaired. "
            "Choose an even number of frames (for example 0..1999 or 0..2001), "
            "or explicitly set series.unpaired_policy='drop_last'."
        )
    last_paired = end - len(unpaired)
    pairs: list[dict] = []
    skipped: list[dict] = []
    seen_paths: dict[str, int] = {}
    for pair_number, first in enumerate(range(start, last_paired + 1, 2)):
        second = first + 1
        missing = [i for i in (first, second) if i not in images]
        if missing:
            skipped.append({
                "pair_number": pair_number,
                "first_index": first,
                "second_index": second,
                "missing_indices": missing,
            })
            continue
        blue = first if first % 2 == (0 if blue_parity == "even" else 1) else second
        red_blue = second if blue == first else first
        for index in (first, second):
            identity = os.path.normcase(str(images[index]))
            if identity in seen_paths:
                raise ValueError(
                    f"Image file reused for frame indices {seen_paths[identity]} and {index}: {images[index]}"
                )
            seen_paths[identity] = index
        elapsed = None if interval is None else pair_number * interval
        if elapsed is not None and not math.isfinite(elapsed):
            raise ValueError("Elapsed time overflow: reduce the range or pair_interval_seconds.")
        pairs.append({
            "pair_number": pair_number,
            "first_index": first,
            "second_index": second,
            "blue_index": blue,
            "red_blue_index": red_blue,
            "blue_path": images[blue],
            "red_blue_path": images[red_blue],
            "time_seconds": elapsed,
        })
    if skipped and missing_policy == "error":
        missing_indices = [i for pair in skipped for i in pair["missing_indices"]]
        preview = ", ".join(map(str, missing_indices[:30]))
        if len(missing_indices) > 30:
            preview += f", ... ({len(missing_indices)} missing total)"
        raise ValueError(
            f"Missing requested frame indices: {preview}. No processing started. "
            "Fix the range/files or explicitly set series.missing_pair_policy='skip' "
            "to omit whole pairs without shifting pairing."
        )
    if not pairs:
        raise ValueError("No complete image pairs remain after applying the requested range and policies.")
    return {
        "start_index": start,
        "end_index": end,
        "blue_parity": blue_parity,
        "unpaired_policy": unpaired_policy,
        "missing_pair_policy": missing_policy,
        "pair_interval_seconds": interval,
        "pairs": pairs,
        "skipped_pairs": skipped,
        "unpaired_indices": unpaired,
    }


def plan_pairs(cfg: Mapping, mapping: Mapping[int, Path]) -> list[dict]:
    """Convenience wrapper returning only the validated complete pairs."""
    return build_pair_plan(cfg, mapping)["pairs"]

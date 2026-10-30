"""Projected pore-area saturation, atomic tabular output, and batch plots.

Saturation here means the fraction of the initial 2-D pore mask occupied by a
phase. It is not a 3-D volume measurement and is unrelated to camera clipping.
The review range is the range obtained by assigning every review-flagged pore
pixel first to air and then to liquid, while holding unflagged labels fixed.
It is a sensitivity range, not a statistical confidence interval; unflagged
pixels have not thereby been independently validated.
"""
from __future__ import annotations

import csv
import math
import os
from pathlib import Path
import tempfile

import numpy as np


def _checked_array(value, name, shape=None):
    array = np.asarray(value)
    if array.ndim != 2 or not array.size:
        raise ValueError(f"{name} must be a nonempty 2-D array.")
    if shape is not None and array.shape != shape:
        raise ValueError(f"{name} shape {array.shape} does not match phase shape {shape}.")
    if not np.issubdtype(array.dtype, np.number) and array.dtype != np.bool_:
        raise ValueError(f"{name} must contain numeric or boolean values.")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains nonfinite values.")
    return array


def saturation_metrics(phase, pore, review_flags=None, B=None, RB=None):
    """Return scalar counts and fractions for one image pair.

    ``phase`` codes are 0=outside, 1=grain, 2=liquid and 3=air. All and only
    pore pixels must carry fluid codes 2 or 3. ``pore`` may be boolean or a
    nonzero-valued binary mask. Any nonzero review flag marks a pixel as
    reviewed/uncertain; flags outside the pore mask do not enter the metrics.

    All returned fractions use 0..1, not percentages. Unreviewed saturation
    uses only unflagged pore pixels as its denominator, returning NaN when
    none remain. Optional raw 8-bit-scale B/RB arrays give separate camera
    clipping fractions at >=255, again within pore space. These are NaN
    when the corresponding channel was not supplied.
    """
    phase = _checked_array(phase, "phase")
    pore = _checked_array(pore, "pore", phase.shape)
    if np.any(pore < 0):
        raise ValueError("pore mask cannot contain negative values.")
    pore = pore.astype(bool)
    if np.any(~np.isin(phase, [0, 1, 2, 3])):
        raise ValueError("phase labels must be 0, 1, 2 or 3.")
    if np.any(pore & ~np.isin(phase, [2, 3])):
        raise ValueError("Every pore pixel must be labeled liquid (2) or air (3).")
    if np.any(~pore & np.isin(phase, [2, 3])):
        raise ValueError("Fluid phase labels occur outside the supplied pore mask.")
    pore_count = int(np.count_nonzero(pore))
    if not pore_count:
        raise ValueError("Pore mask is empty; saturation has no denominator.")

    if review_flags is None:
        reviewed = np.zeros(phase.shape, dtype=bool)
    else:
        flags = _checked_array(review_flags, "review_flags", phase.shape)
        if np.any(flags < 0):
            raise ValueError("review_flags cannot contain negative values.")
        reviewed = pore & (flags != 0)
    unreviewed = pore & ~reviewed
    liquid = pore & (phase == 2)
    air = pore & (phase == 3)
    n_liquid = int(np.count_nonzero(liquid))
    n_air = int(np.count_nonzero(air))
    n_reviewed = int(np.count_nonzero(reviewed))
    n_unreviewed = pore_count - n_reviewed
    n_liquid_unreviewed = int(np.count_nonzero(liquid & unreviewed))
    n_air_unreviewed = int(np.count_nonzero(air & unreviewed))
    sw = n_liquid / pore_count
    sa = n_air / pore_count

    result = {
        "pore_pixel_count": pore_count,
        "liquid_pixel_count": n_liquid,
        "air_pixel_count": n_air,
        "liquid_saturation": sw,
        "air_saturation": sa,
        "nominal_liquid_fraction_of_pore": sw,
        "nominal_air_fraction_of_pore": sa,
        "reviewed_pore_pixels": n_reviewed,
        "unreviewed_pore_pixels": n_unreviewed,
        "unreviewed_liquid_pixels": n_liquid_unreviewed,
        "unreviewed_air_pixels": n_air_unreviewed,
        "review_fraction": n_reviewed / pore_count,
        "liquid_saturation_unreviewed": (
            n_liquid_unreviewed / n_unreviewed if n_unreviewed else float("nan")),
        "air_saturation_unreviewed": (
            n_air_unreviewed / n_unreviewed if n_unreviewed else float("nan")),
        "liquid_saturation_lower_possible": n_liquid_unreviewed / pore_count,
        "liquid_saturation_upper_possible": (
            n_liquid_unreviewed + n_reviewed) / pore_count,
        "air_saturation_lower_possible": n_air_unreviewed / pore_count,
        "air_saturation_upper_possible": (n_air_unreviewed + n_reviewed) / pore_count,
    }
    for name, channel in (("blue", B), ("RB", RB)):
        clipped_count = None
        fraction = float("nan")
        if channel is not None:
            channel = _checked_array(channel, name, phase.shape)
            if channel.min() < 0 or channel.max() > 255:
                raise ValueError(f"{name} camera values must be on the 0..255 scale.")
            clipped_count = int(np.count_nonzero(pore & (channel >= 255)))
            fraction = clipped_count / pore_count
        result[f"{name}_clipped_pore_pixels"] = clipped_count
        result[f"{name}_clipped_fraction"] = fraction
    return result


def save_csv_atomic(rows, path, fieldnames=None):
    """Atomically save the complete accumulated list of per-pair rows.

    Call after each pair, supplying *all* prior rows plus the new row, including
    failures. A write/replace failure preserves the previous complete CSV.
    The temporary file stays in the destination directory so replacement
    cannot cross filesystems. Extra row keys are never silently discarded.
    Missing values and nonfinite numeric values are stored as blank cells.
    """
    rows = list(rows)
    if fieldnames is None:
        columns = list(dict.fromkeys(key for row in rows for key in row))
    else:
        columns = list(fieldnames)
    if not columns:
        raise ValueError("CSV requires at least one column.")
    if len(set(columns)) != len(columns):
        raise ValueError("CSV fieldnames contain duplicate columns.")
    extras = {key for row in rows for key in row} - set(columns)
    if extras:
        raise ValueError(f"CSV fieldnames omit row keys: {sorted(extras)}")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_name = None
    try:
        with tempfile.NamedTemporaryFile(
                mode="w", encoding="utf-8", newline="", delete=False,
                prefix="." + path.name + ".", suffix=".tmp", dir=path.parent) as handle:
            temp_name = handle.name
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for row in rows:
                clean = {}
                for key, value in row.items():
                    if isinstance(value, np.generic):
                        value = value.item()
                    if value is None or (isinstance(value, float) and not math.isfinite(value)):
                        value = ""
                    clean[key] = value
                writer.writerow(clean)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_name, path)
        temp_name = None
    finally:
        if temp_name is not None:
            try:
                os.unlink(temp_name)
            except FileNotFoundError:
                pass
    return path


def _finite_float(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def write_saturation_plot(rows, path, xaxis="pair_number"):
    """Plot nominal phase saturation and review fraction for finite rows.

    ``xaxis`` is ``pair_number`` or ``time_seconds``. No status spelling is
    required: rows without finite saturation/x values are omitted, allowing
    failed, missing and incomplete pair rows to stay in the source CSV.
    Existing Matplotlib backend selection is respected (including Spyder).
    """
    if xaxis not in {"pair_number", "time_seconds"}:
        raise ValueError("xaxis must be 'pair_number' or 'time_seconds'.")
    import matplotlib.pyplot as plt

    records = []
    for row in rows:
        x = _finite_float(row.get(xaxis))
        sw = _finite_float(row.get("liquid_saturation"))
        sa = _finite_float(row.get("air_saturation"))
        if np.isfinite(x) and np.isfinite(sw) and np.isfinite(sa):
            records.append((x, sw, sa,
                            _finite_float(row.get("liquid_saturation_lower_possible")),
                            _finite_float(row.get("liquid_saturation_upper_possible")),
                            _finite_float(row.get("review_fraction"))))
    records.sort(key=lambda values: values[0])
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1]}, constrained_layout=True)
    main, review = axes
    if records:
        x, sw, sa, lower, upper, flagged = np.asarray(records, dtype=float).T
        main.plot(x, 100 * sw, "o-", color="#1488cc", markersize=3, label="Liquid (nominal)")
        main.plot(x, 100 * sa, "o-", color="#e59631", markersize=3, label="Air (nominal)")
        if np.any(np.isfinite(lower) & np.isfinite(upper)):
            main.fill_between(x, 100 * lower, 100 * upper,
                              color="#1488cc", alpha=0.15,
                              label="Liquid review range (not a confidence interval)")
        review.plot(x, 100 * flagged, "o-", color="#7c578b", markersize=3)
        if x.size == 1:
            review.set_xlim(x[0] - 0.5, x[0] + 0.5)
        main.legend(loc="best", fontsize=8)
    else:
        main.text(0.5, 0.5, "No completed pairs with finite saturation values",
                  ha="center", va="center", transform=main.transAxes)
    main.set(title="Projected pore-area phase saturation", ylabel="Pore area (%)", ylim=(0, 100))
    review.set(ylabel="For review (%)", ylim=(0, 100),
               xlabel="Elapsed time (s)" if xaxis == "time_seconds" else "Pair number")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.suptitle("Review range holds unflagged labels fixed; unflagged does not mean validated.",
                 fontsize=9)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.savefig(path, dpi=160)
    finally:
        plt.close(fig)
    return path

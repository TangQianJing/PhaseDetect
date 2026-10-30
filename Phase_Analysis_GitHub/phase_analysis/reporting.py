"""Portable, illustrated reports for the fluorescence-first phase workflow.

Only numerical arrays, ordinary configuration dictionaries and local output
paths are used. Report images are scientific visualizations of the measured
arrays; drawings are explicitly identified as conceptual schematics.
"""
from __future__ import annotations

import html
import json
import gc
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle, Patch


COLORS = ["#e6e7e9", "#343b45", "#3199cf", "#ee9c39", "#b749a3"]
PHASE_CMAP = ListedColormap(COLORS)
PHASE_NORM = BoundaryNorm(np.arange(-.5, 5.5), len(COLORS))
EXAMPLES = [
    {"name": "Candidate 2 / dark-liquid and inclusion example", "xyxy": [1155, 1025, 1335, 1205]},
    {"name": "Broader fluorescence example", "xyxy": [760, 1480, 1040, 1760]},
    {"name": "Earlier capsule example", "xyxy": [1469, 709, 1579, 819]},
]


def _plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return value


def _save(fig, directory, stem):
    path = directory / (stem + ".png")
    fig.savefig(path, dpi=145, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    gc.collect()
    return "figures/" + path.name


def _display(array, extent=None):
    """Limit figure memory; exported analytical TIFFs retain the original grid."""
    h, w = array.shape[:2]
    if extent is None:
        extent = (0, w, h, 0)
    step = max(1, int(np.ceil(max(h, w)/1800)))
    return array[::step, ::step], extent


def _cutoff(cfg):
    nested = cfg.get("fluorescence", {})
    return float(nested.get("threshold", cfg.get("fluorescence_threshold", 55)))


def _reference_kind(data, cfg):
    return data.get("reference_kind", "binary_reference" if
                    cfg.get("mask_input", {}).get("mode") == "load_binary"
                    and not cfg.get("mask_input", {}).get("display_reference_path")
                    else "raw_initial")


def _reference_labels(data, cfg):
    kind = _reference_kind(data, cfg)
    if kind == "binary_reference":
        return ("Binary geometry reference (synthetic)", "Binary geometry reference + grain overlay",
                "Reference: synthetic binary geometry (0/255), not a camera image.")
    if kind == "display_reference":
        return ("Measured display reference / PDMS geometry", "Grain mask on measured display reference",
                "Reference: separately supplied display image, stretched to its 2–99th percentiles.")
    return ("Initial camera image / PDMS geometry", "Initial PDMS mask overlay",
            "Reference: initial camera image, stretched to its 2–99th percentiles.")


def _sensitivity_status(data, cfg):
    """Describe actual alternatives, without treating a nominal-only run as a test."""
    threshold = data.get("threshold_sensitivity", {})
    if isinstance(threshold, dict):
        records = threshold.get("records", threshold.get("metrics"))
        if records is None:
            values = list(threshold)
        else:
            values = [r.get("threshold", r.get("B_threshold", r.get("fluorescence_threshold")))
                      for r in records if isinstance(r, dict)]
    else:
        values = [r.get("threshold", r.get("B_threshold", r.get("fluorescence_threshold")))
                  for r in threshold if isinstance(r, dict)]
    cutoffs = sorted({float(v) for v in values if v is not None and float(v) != _cutoff(cfg)})
    if cutoffs:
        fluorescence = ("Fluorescence-threshold sensitivity was calculated at "
                        + ", ".join(f"{v:g}" for v in cutoffs)
                        + f" counts in addition to the nominal {_cutoff(cfg):g}-count cutoff.")
    elif not cfg.get("fluorescence_sensitivity", []):
        fluorescence = ("Fluorescence-threshold sensitivity was disabled for this run; "
                        f"only the nominal {_cutoff(cfg):g}-count cutoff was used. "
                        "The review flags therefore do not assess alternative fluorescence cutoffs.")
    else:
        fluorescence = "No alternate fluorescence-cutoff results are available in this report."
    records = data.get("watershed_sensitivity", data.get("summary", {}).get("watershed_sensitivity", []))
    if isinstance(records, dict):
        records = records.get("records", records.get("metrics", []))
    prominences = sorted({float(r["watershed_h"]) for r in records
                          if isinstance(r, dict) and "watershed_h" in r})
    if len(prominences) > 1:
        watershed = ("Watershed sensitivity was calculated at h="
                     + ", ".join(f"{v:g}" for v in prominences)
                     + " counts. Agreement measures parameter stability, not experimental accuracy.")
    elif records or cfg.get("segmentation", {}).get("watershed_sensitivity", True) is False \
            or cfg.get("segmentation", {}).get("watershed_sensitivity") == []:
        watershed = ("Watershed-parameter sensitivity was disabled or no alternate prominence was run. "
                     "A nominal-only run does not establish watershed stability; other review flags still apply.")
    else:
        watershed = "No alternate watershed-prominence results are available in this report."
    return fluorescence, watershed


def _image(ax, array, title="", limits=None, extent=None, cmap="gray"):
    array, extent = _display(array, extent)
    kws = {"cmap": cmap, "origin": "upper", "interpolation": "nearest"}
    if limits is not None:
        kws.update(vmin=limits[0], vmax=limits[1])
    if extent is not None:
        kws["extent"] = extent
    im = ax.imshow(array, **kws)
    ax.set_title(title, fontsize=11, pad=7)
    ax.tick_params(labelsize=8)
    return im


def _phase(ax, array, title="", extent=None):
    array, extent = _display(array, extent)
    im = ax.imshow(array, origin="upper", interpolation="nearest", extent=extent,
                   cmap=PHASE_CMAP, norm=PHASE_NORM)
    ax.set_title(title, fontsize=11)
    ax.tick_params(labelsize=8)
    return im


def _legend(fig, review=False, y=.01):
    labels = ["Outside bed", "Fixed PDMS", "Liquid", "Air"]
    if review:
        labels.append("Review / unresolved")
    fig.legend(handles=[Patch(color=c, label=l) for c, l in zip(COLORS, labels)],
               loc="lower center", ncol=len(labels), frameon=False,
               fontsize=9, bbox_to_anchor=(.5, y))


def _overlay(ax, mask, color, alpha=.6, extent=None):
    mask, extent = _display(mask, extent)
    rgba = np.zeros((*mask.shape, 4), dtype=np.float32)
    rgba[..., :3] = matplotlib.colors.to_rgb(color)
    rgba[..., 3] = np.asarray(mask, bool) * alpha
    ax.imshow(rgba, origin="upper", interpolation="nearest", extent=extent)


def _grain_contour(ax, grain, xyxy, color="#ad72b8", linewidth=.65):
    x0, y0, x1, y1 = xyxy
    local = grain[y0:y1, x0:x1]
    if local.any() and not local.all():
        ax.contour(np.arange(x0, x1) + .5, np.arange(y0, y1) + .5,
                   local.astype(float), levels=[.5], colors=[color], linewidths=linewidth)


def _crop(array, xyxy):
    x0, y0, x1, y1 = xyxy
    return array[y0:y1, x0:x1]


def _examples(cfg, shape):
    height, width = shape
    examples = []
    for item in cfg.get("examples", EXAMPLES):
        x0, y0, x1, y1 = map(int, item["xyxy"])
        x0, x1 = np.clip([x0, x1], 0, width)
        y0, y1 = np.clip([y0, y1], 0, height)
        if x1 > x0 + 8 and y1 > y0 + 8:
            examples.append({"name": str(item.get("name", "Example")),
                             "xyxy": [int(x0), int(y0), int(x1), int(y1)]})
    if not examples:
        examples = [{"name": "Central example", "xyxy": [width//4, height//4,
                                                              3*width//4, 3*height//4]}]
    return examples


def _topology_example(examples):
    for item in examples:
        name = item["name"].lower()
        if (item["xyxy"] == [1155, 1025, 1335, 1205]
                or "inclusion" in name or "candidate 2" in name or "dark water" in name):
            return item
    return examples[0]


def _profile_figure(data, cfg, directory):
    profiles = data.get("profiles", {})
    yx = np.asarray(profiles.get("yx", []))
    normals = np.asarray(profiles.get("normals_xy", []))
    if yx.ndim != 2 or normals.shape != yx.shape or not len(yx):
        fig, ax = plt.subplots(figsize=(9, 4))
        ax.text(.5, .5, "No wall-profile samples available for these inputs.", ha="center")
        ax.set_axis_off()
        return _save(fig, directory, "03_measured_wall_profiles")
    n = len(yx)
    valid = np.asarray(profiles.get("valid", np.ones(n, bool)), bool)
    outer = np.asarray(profiles.get("B_outer", np.zeros(n)), float)
    contrast = np.asarray(profiles.get("B_contrast", np.zeros(n)), float)
    lit = np.asarray(profiles.get("lit", np.zeros(n, bool)), bool)
    # Select actual measured examples, without requiring all categories to exist.
    broad = np.flatnonzero(valid & (outer > _cutoff(cfg)))
    narrow = np.flatnonzero(valid & lit & (outer < 45) & (contrast > 10))
    weak = np.flatnonzero(valid & ~lit)
    choices = []
    for title, options, score in [("Broad fluorescent water", broad, outer),
                                   ("Narrow bright wall fringe", narrow, contrast),
                                   ("Weak / unlit wall sample", weak, -contrast)]:
        if len(options):
            ordered = options[np.argsort(score[options])]
            # A high, but not extreme, measured example avoids isolated outliers.
            index = int(ordered[min(len(ordered)-1, int(.80*len(ordered)))])
            choices.append((title, index, True))
        else:
            fallback = np.flatnonzero(valid)
            index = int(fallback[0]) if len(fallback) else 0
            choices.append((title + "\n(category absent; alternate sample)", index, False))
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.2))
    offsets = np.linspace(-7, 18, 101)
    channel_sigma = float(cfg.get("wall", {}).get("channel_sigma", .9))
    smooth_channels = {name: ndi.gaussian_filter(np.asarray(data[name], np.float32), channel_sigma)
                       for name in ("B", "RB")}
    for col, (title, index, present) in enumerate(choices):
        y, x = yx[index]
        nx, ny = normals[index]
        height, width = data["B"].shape
        box = [max(int(x)-25, 0), max(int(y)-25, 0),
               min(int(x)+26, width), min(int(y)+26, height)]
        extent = (box[0], box[2], box[3], box[1])
        ax = axes[0, col]
        _image(ax, _crop(data["B"], box), title, (6, 120), extent)
        _grain_contour(ax, data["grain"], box)
        # Imshow extent puts row/column i at display center i+0.5.
        # Interpolation below still uses the true NumPy row/column coordinates.
        ax.plot(x + .5 + offsets*nx, y + .5 + offsets*ny, color="#ffe16d", lw=1.2)
        ax.plot(x + .5, y + .5, "o", color="#e53e61", ms=4)
        ax.set_xlim(box[0], box[2]); ax.set_ylim(box[3], box[1])
        ax.set_xlabel(f"Source x; sample ({x:.1f}, {y:.1f})", fontsize=9)
        values = {}
        for name, color in [("B", "#2675b7"), ("RB", "#bd5262")]:
            vals = ndi.map_coordinates(smooth_channels[name],
                                      [y+offsets*ny, x+offsets*nx], order=1, mode="nearest")
            values[name] = vals
            axes[1, col].plot(offsets, vals, color=color, label="R&B" if name == "RB" else "B")
        axes[1, col].axvline(0, color="#ad72b8", ls="--", lw=1,
                            label="Sampled wall location")
        axes[1, col].axhline(_cutoff(cfg),
                            color="#2675b7", ls=":", lw=1, label="B cutoff")
        axes[1, col].set_xlabel("Distance along outward PDMS normal (pixels)", fontsize=9)
        axes[1, col].set_ylabel("Smoothed camera counts", fontsize=9)
        axes[1, col].grid(alpha=.18)
        axes[1, col].tick_params(labelsize=8)
        axes[1, col].set_title(f"B crest contrast {contrast[index]:.1f}; outer B {outer[index]:.1f}", fontsize=10)
    axes[1, 0].legend(fontsize=7, loc="best")
    fig.suptitle("Measured wall-normal profiles: a sustained bright interior differs from a narrow fringe", fontsize=13)
    fig.text(.5, .015, f"Purple = fixed grain boundary. Yellow = measured profile. Curves use σ={channel_sigma:g} px channel smoothing. Positive normal: PDMS → pore; not automatically air.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .045, 1, .96))
    return _save(fig, directory, "03_measured_wall_profiles")


def _gradient_figure(data, examples, directory, cfg):
    example = _topology_example(examples)
    box = example["xyxy"]
    extent = (box[0], box[2], box[3], box[1])
    b = _crop(data["B"], box).astype(np.float32)
    rb = _crop(data["RB"], box)
    sigma = float(cfg.get("gradient_sigma_px", .9))
    smooth = ndi.gaussian_filter(b, sigma)
    gy, gx = np.gradient(smooth)
    magnitude = np.hypot(gx, gy)
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 10.2))
    _image(axes[0, 0], b, "Measured B + fixed PDMS boundary", (6, 120), extent)
    _image(axes[0, 1], rb, "Measured R&B + fixed PDMS boundary", (15, 255), extent)
    for ax in axes[0]:
        _grain_contour(ax, data["grain"], box)
    ax = axes[1, 0]
    _image(ax, b, "Measured −∇B: bright → dark on BOTH flanks", (6, 120), extent)
    _grain_contour(ax, data["grain"], box)
    step = max(4, int(max(b.shape)/34))
    ys, xs = np.mgrid[step//2:b.shape[0]:step, step//2:b.shape[1]:step]
    local_mag = magnitude[ys, xs]
    threshold = max(1.5, float(np.percentile(magnitude, 70)))
    good = local_mag > threshold
    length = max(2.5, step*.68)
    u = -gx[ys, xs]/np.maximum(local_mag, 1e-6)*length
    v = -gy[ys, xs]/np.maximum(local_mag, 1e-6)*length
    ax.quiver(xs[good]+box[0]+.5, ys[good]+box[1]+.5,
              u[good], v[good], color="#1dd5d1", angles="xy", scale_units="xy",
              scale=1, width=.0045, headwidth=3.3)
    ax.set_xlim(box[0], box[2]); ax.set_ylim(box[3], box[1])
    ax = axes[1, 1]
    s = np.linspace(-7, 7, 280)
    line = 24 + 85*np.exp(-s*s/(2*.95**2))
    ax.plot(s, line, color="#39495a", lw=2)
    ax.axvline(0, ls="--", lw=1, color="#d38a2a", label="Bright-fringe centerline")
    for start, end in [(-.7, -3.0), (.7, 3.0)]:
        ax.annotate("", (end, 44), (start, 44), arrowprops={"arrowstyle":"->", "color":"#00a6a2", "lw":2})
    ax.text(-3.1, 51, "−∇I", color="#008d8a", ha="center")
    ax.text(3.1, 51, "−∇I", color="#008d8a", ha="center")
    ax.text(0, 120, "Conceptual profile — not measured data", ha="center", fontsize=10)
    ax.text(-4, 18, "Dark side A", ha="center", fontsize=10)
    ax.text(4, 18, "Dark side B", ha="center", fontsize=10)
    ax.set_ylim(0, 135); ax.set_xlabel("Position across a thin fringe (arbitrary units)")
    ax.set_ylabel("Brightness (arbitrary units)")
    ax.set_title("A bright line alone does not assign the air side", fontsize=11)
    ax.legend(fontsize=8, loc="upper right")
    fig.suptitle(example["name"] + ": gradient topology", fontsize=13)
    fig.text(.5, .015, f"Arrows are normalized directions after σ={sigma:g} px smoothing. Low-gradient samples are omitted; arrow length does not represent gradient strength.",
             ha="center", fontsize=8.5)
    fig.tight_layout(rect=(0, .04, 1, .96))
    return _save(fig, directory, "04_gradients_and_two_sided_fringe")


def _junction_figure(directory):
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.9))
    ax = axes[0]
    ax.add_patch(Rectangle((-2, -2), 2, 4, color=COLORS[1]))
    ax.plot([0, 0], [-2, 2], color="#ad72b8", lw=3)
    ax.plot([.035, .035], [-2, 0], color="#ffe26c", lw=6)
    x = np.linspace(0, 3.7, 100); y = .32*np.sin(x*.7)
    ax.fill_between(x, -2, y, color=COLORS[3], alpha=.18)
    ax.fill_between(x, y, 2, color=COLORS[2], alpha=.18)
    ax.plot(x, y, color="#ffe26c", lw=6)
    ax.scatter([0], [0], c="#d0397a", s=55, zorder=4)
    ax.text(-1, .8, "PDMS", color="white", ha="center", fontsize=12)
    ax.text(2.4, -.95, "Air candidate", fontsize=12, color="#9d5a13")
    ax.text(2.4, 1, "Liquid candidate", fontsize=12, color="#206a97")
    ax.annotate("Lit wall branch", (.04, -1), (.65, -1.65), fontsize=10,
                arrowprops={"arrowstyle":"->", "color":"#555"})
    ax.annotate("Free bright fringe", (2.5, .31), (1.1, 1.55), fontsize=10,
                arrowprops={"arrowstyle":"->", "color":"#555"})
    ax.set_xlim(-2, 4); ax.set_ylim(-2, 2); ax.set_aspect("equal"); ax.set_axis_off()
    ax.set_title("Concept: connected wall and free fringe establish a side", fontsize=11)
    ax = axes[1]
    ax.set_axis_off()
    ax.text(.04, .94, "Implemented approximation", weight="bold", fontsize=13, va="top")
    text = ("1  Keep reference PDMS geometry fixed.\n\n"
            "2  Sample wall normals and classify lit / unlit profiles.\n\n"
            "3  Aggregate wall evidence within raw R&B basins.\n\n"
            "4  Connect neighboring basins across supported free fringes.\n\n"
            "5  Use wall-anchored air neighbors for cautious liquid votes.\n\n"
            "Fluorescence-confirmed water remains fixed throughout.\n"
            "A precise contact-angle or triple-junction fit is not performed.")
    ax.text(.04, .83, text, fontsize=11, va="top", linespacing=1.3)
    fig.text(.5, .01, "Conceptual illustration only: the colored sectors are not a traced junction from the experiment. Connection need not be tangent; a contact angle is allowed.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .05, 1, 1))
    return _save(fig, directory, "05_wall_junction_principle")


def _basin_figure(data, examples, directory):
    box = _topology_example(examples)["xyxy"]
    extent = (box[0], box[2], box[3], box[1])
    labels = _crop(np.asarray(data.get("region_labels", np.zeros(data["B"].shape))), box)
    edge = np.zeros(labels.shape, bool)
    edge[1:, :] |= labels[1:, :] != labels[:-1, :]
    edge[:, 1:] |= labels[:, 1:] != labels[:, :-1]
    edge &= _crop(data["pore"], box)
    fig, axes = plt.subplots(2, 2, figsize=(10.9, 9.7))
    _image(axes[0, 0], _crop(data["RB"], box), "Raw R&B + watershed boundaries", (15, 255), extent)
    _overlay(axes[0, 0], edge, "#36d7c4", .85, extent)
    _grain_contour(axes[0, 0], data["grain"], box)
    # Reproducible pseudo-colors indicate optical basins, without phase meaning.
    pseudo = ((labels.astype(np.int64)*2654435761) % 997)/997
    shown = np.ma.masked_where(labels == 0, pseudo)
    _image(axes[0, 1], shown, "Optical basin identities (colors have no phase meaning)", (0, 1), extent, cmap="turbo")
    _image(axes[1, 0], _crop(data["B"], box), "Fixed water + sampled wall evidence", (6, 120), extent)
    _overlay(axes[1, 0], _crop(data["direct_water"], box), COLORS[2], .58, extent)
    _grain_contour(axes[1, 0], data["grain"], box)
    profiles = data.get("profiles", {})
    yx = np.asarray(profiles.get("yx", []))
    if yx.ndim == 2 and len(yx):
        inside = (yx[:, 1] >= box[0]) & (yx[:, 1] < box[2]) & (yx[:, 0] >= box[1]) & (yx[:, 0] < box[3])
        for key, color, label in [("lit", "#f6cc3c", "Lit-wall sample"), ("unlit", "#42d6a8", "Unlit-wall sample")]:
            use = inside & np.asarray(profiles.get(key, np.zeros(len(yx), bool)), bool)
            axes[1, 0].scatter(yx[use, 1]+.5, yx[use, 0]+.5, s=4, c=color, label=label, linewidths=0)
        axes[1, 0].legend(fontsize=8, loc="lower right")
    _phase(axes[1, 1], _crop(data["phase"], box), "Nominal phase after fixed seeds and wall/graph votes", extent)
    for ax in axes.ravel():
        ax.set_xlim(box[0], box[2]); ax.set_ylim(box[3], box[1])
        ax.set_xlabel("Source x (pixels)", fontsize=9); ax.set_ylabel("Source y (pixels)", fontsize=9)
    _legend(fig, y=.01)
    fig.tight_layout(rect=(0, .045, 1, 1))
    return _save(fig, directory, "10_local_watershed_and_wall_evidence")


def _sensitivity_figure(data, cfg, directory):
    fig, axes = plt.subplots(1, 2, figsize=(11.7, 4.6))
    fluorescence_status, watershed_status = _sensitivity_status(data, cfg)
    threshold = data.get("threshold_sensitivity", {})
    # Accept a list of records, a keyed dictionary of records, or a records field.
    if isinstance(threshold, dict):
        if "records" in threshold or "metrics" in threshold:
            records = threshold.get("records", threshold.get("metrics", []))
        else:
            records = [dict(row, threshold=key) for key, row in threshold.items()
                       if isinstance(row, dict)]
    else:
        records = threshold
    tx, ty, total_x, total_y = [], [], [], []
    if isinstance(records, list):
        for row in records:
            if not isinstance(row, dict):
                continue
            t = row.get("threshold", row.get("B_threshold", row.get("fluorescence_threshold")))
            y = row.get("direct_water_fraction_of_pore", row.get("water_seed_fraction_of_pore", row.get("direct_fraction")))
            if y is None and "direct_water" in row:
                mask = np.asarray(row["direct_water"], bool)
                pore = np.asarray(data["pore"], bool)
                if mask.shape == pore.shape and pore.any():
                    y = float(mask[pore].mean())
            if t is not None and y is not None:
                tx.append(float(t)); ty.append(float(y)*100)
            total = row.get("liquid_fraction", row.get("liquid_fraction_of_pore"))
            if t is not None and total is not None:
                total_x.append(float(t)); total_y.append(float(total)*100)
    if tx and _cutoff(cfg) not in tx:
        pore = np.asarray(data["pore"], bool)
        tx.append(_cutoff(cfg)); ty.append(100*float(np.asarray(data["direct_water"])[pore].mean()))
    if total_x and _cutoff(cfg) not in total_x:
        pore = np.asarray(data["pore"], bool)
        total_x.append(_cutoff(cfg)); total_y.append(100*float((np.asarray(data["phase"])[pore] == 2).mean()))
    if tx:
        order = np.argsort(tx)
        axes[0].plot(np.array(tx)[order], np.array(ty)[order], "o-", color=COLORS[2], label="Fixed fluorescent water")
        if total_x:
            ti = np.argsort(total_x)
            axes[0].plot(np.array(total_x)[ti], np.array(total_y)[ti], "s--", color="#214b69", label="Total nominal liquid")
            axes[0].legend(fontsize=8)
        axes[0].set_xlabel("B fluorescence cutoff (raw camera counts)")
        axes[0].set_ylabel("Fraction of pore pixels (%)")
        axes[0].grid(alpha=.2)
    else:
        axes[0].set_axis_off()
        import textwrap
        axes[0].text(.5, .55, textwrap.fill(fluorescence_status, 49),
                      ha="center", va="center", fontsize=10, transform=axes[0].transAxes)
    axes[0].set_title("Fluorescence-threshold sensitivity", fontsize=12)
    summary = data.get("summary", {})
    records = data.get("watershed_sensitivity", summary.get("watershed_sensitivity", []))
    if isinstance(records, dict):
        records = records.get("records", records.get("metrics", []))
    x, water, air = [], [], []
    if isinstance(records, list):
        for row in records:
            if not isinstance(row, dict) or "watershed_h" not in row:
                continue
            if "nominal_liquid_fraction_of_pore" in row:
                x.append(row["watershed_h"])
                water.append(100*row["nominal_liquid_fraction_of_pore"])
                air.append(100*row["nominal_air_fraction_of_pore"])
    if len(set(x)) > 1:
        order = np.argsort(x)
        axes[1].plot(np.array(x)[order], np.array(water)[order], "o-", color=COLORS[2], label="Nominal liquid")
        axes[1].plot(np.array(x)[order], np.array(air)[order], "o-", color=COLORS[3], label="Nominal air")
        axes[1].set_xlabel("Watershed h-minimum prominence (counts)")
        axes[1].set_ylabel("Fraction of pore pixels (%)")
        axes[1].legend(fontsize=9); axes[1].grid(alpha=.2)
    else:
        axes[1].set_axis_off()
        import textwrap
        axes[1].text(.5, .55, textwrap.fill(watershed_status, 49),
                      ha="center", va="center", fontsize=10, transform=axes[1].transAxes)
    axes[1].set_title("Watershed-parameter sensitivity", fontsize=12)
    fig.text(.5, .015, "Sensitivity measures parameter dependence; agreement does not establish experimental accuracy or a statistical confidence interval.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .06, 1, 1))
    return _save(fig, directory, "09_parameter_sensitivity")


PARAMETER_NOTES = [
    ("mask_input.mode", "binarize / load_binary", "Create geometry from a separate camera image, or use an already binary grain mask exactly."),
    ("mask_input.white_is_grain", "True / False", "For external masks, choose whether nonzero or zero means PDMS."),
    ("mask_input.bed_mask_path", "optional TIFF / PNG path", "Nonzero pixels are inside the bed; intersected with the configured polygon/ROI."),
    ("mask_input.display_reference_path", "optional TIFF / PNG path", "Measured reference for figures; does not alter an external grain mask."),
    ("grain.window", "pixels; odd integer", "Local Phansalkar neighborhood; must span grain and neighboring pore structure."),
    ("grain.sigma", "pixels", "Smoothing before initial-grain thresholding."),
    ("fluorescence.threshold", "camera counts", "Positive B cutoff for fixed liquid pixels."),
    ("fluorescence.sigma", "pixels", "B smoothing used for the fluorescence cutoff."),
    ("fluorescence.shiny_distance", "pixels", "Near-wall range in which a valid narrow-fringe exception may apply."),
    ("fluorescence.shiny_contrast", "camera counts", "Minimum narrow B crest contrast for the shiny exception."),
    ("fluorescence.outer_margin", "camera counts", "Outer B must be below fluorescence cutoff minus this margin."),
    ("blue_scale", "dimensionless", "B multiplier in the signed R&B minus scaled-B difference."),
    ("wall.normal_sigma", "pixels", "Smoothing of the signed grain-distance field before normal estimation."),
    ("wall.channel_sigma", "pixels", "Channel smoothing before measured wall-profile features."),
    ("wall.lit_rb", "camera counts", "One criterion for a lit R&B wall crest; alternate combined criteria also exist."),
    ("segmentation.watershed_h", "camera counts", "Nominal h-minimum prominence; affects basin splitting and merging."),
    ("segmentation.watershed_sensitivity", "camera counts", "Additional watershed prominences to compare with the nominal labels."),
    ("segmentation.air_wall_lit_fraction", "fraction of valid wall samples", "One of several requirements for an accepted air-region proposal."),
    ("segmentation.broad_water_region_fraction", "fraction of basin pixels", "Minimum substantial direct-water coverage for extending a liquid-region hypothesis."),
    ("fluorescence_sensitivity", "camera counts", "Additional B cutoffs; each reruns fixed-water and phase assignment."),
    ("gradient_sigma_px", "pixels", "Smoothing for the illustrative measured gradient arrows."),
]


GUIDE_SECTIONS = [
    ("Purpose and what the map means", """This workflow estimates liquid–air distribution in the initial pore space of the porous bed. Liquid includes both fluorescent water and dark diluted or fresh water. The initial PDMS grains are a separate fixed class. It uses the initial frame, blue-only fluorescence and the red-plus-blue interface image together. It does not assume that darkness means air. The nominal map supplies a complete display partition; the review map retains unresolved or weakly supported areas. Both should accompany any reported phase fractions. These are optical proposals from two-dimensional images, without manual capsule outlines or phase training labels.""", "06_whole_image_phase_and_review"),
    ("1. Read the original camera images", """Provide three two-dimensional TIFFs of identical shape: the initial frame 0001, R&B frame 6544, and B frame 6545. Edit input paths and output location in settings.py; open RUN_ANALYSIS.py in Spyder and press F5 to execute the whole workflow. The program reads camera values without normalizing each image independently. Thresholds therefore remain in the original camera-count units; the supplied settings were developed for these 8-bit images. Different gain, exposure, bit depth or illumination requires revisiting the settings. Display contrast is only for figures. B and R&B must refer to the same scene with sufficiently small motion; this implementation does not automatically register images or correct deformation and changing liquid geometry. Optional initial_shift_yx and blue_shift_yx apply a user-specified translation; verify any shift against grain landmarks. A shifted grain edge can look like an interface. The report records the settings and input provenance so the run can be reproduced.""", "01_input_images"),
    ("2. Establish the bed and initial PDMS reference", """The initial frame supplies grain geometry before the later phase interpretation. Uneven illumination is handled by a local Phansalkar threshold, using the configured window (nominally 101 pixels), rather than one global gray-level cutoff. For normalized initial intensity, the threshold is t=m[1+p exp(−qm)+k(s/r−1)], where m and s are the local mean and standard deviation. Nominal values are k=0.25, r=0.5, p=2 and q=10. Local statistics allow a dark grain to be separated from brighter pore space despite slow brightness variation. Small disconnected objects are removed and enclosed grain holes are filled. The bed mask restricts analysis to the experimental porous area. Pixels inside the bed are divided into PDMS and pores; grains are not expanded by a blanket dilation. This matters because real liquid films and menisci can be only a few pixels from a wall. Verify the purple overlay against the original frame, especially around small grains, dim edges and the bed perimeter. An incorrect initial mask propagates to every subsequent wall-normal and phase estimate. Distances and normals are pixel-based: there is no physical length calibration unless one is supplied separately.""", "02_geometry_fluorescence_and_difference"),
    ("3. Fix fluorescent water first", """B is Gaussian-smoothed gently (nominal σ=0.9 pixels), then compared with the fluorescence cutoff (nominally 55 camera counts). Bright pore pixels become fixed liquid seeds before the interface-based classification. This is positive, one-sided evidence: the dye can fluoresce only in water, but absence of fluorescence does not exclude water. Every later hypothesis and watershed sensitivity run respects these fixed seeds. A bright pixel fixes that pixel, rather than automatically converting its entire optical basin to liquid. A basin with a substantial fluorescent interior can provide additional, explicitly weaker liquid evidence for its remaining pixels. The raw B and R&B images remain available for fringe and geometry measurements after the seeds are established; replacing them with a binary seed mask would discard useful structure.""", None),
    ("4. Withhold narrow shiny grain fringes", """Some dry-grain outlines are bright even in B. To avoid forcing all such outlines to liquid, the program samples brightness along outward normals to the initial PDMS boundary. The signed distance to the initial grain supplies a smoothed outward normal. Samples at offsets −4, −2, 0, 1, 2, 3, 4, 6 and 8 pixels measure grain-side baseline, near-wall crest and outer pore brightness. A valid profile needs grain-side and outer-pore samples and adequate clearance from the bed edge. Crest contrast is the near peak minus the larger of the grain and outer-pore baselines. A narrow B peak above a darker neighboring pore is treated differently from brightness sustained across a broad water-filled pore. The nominal exception uses wall proximity of 4.5 pixels, B crest contrast above 10 counts, and adjacent pore B below 45 counts, together with a valid wall profile. Only the corresponding near-wall pixels are withheld from fixed fluorescence labeling; they are not automatically classified as air. Broad fluorescent water remains liquid up to the wall. Purple boundaries, sampled lines, measured B/R&B curves and the actual selected profile coordinates are shown in the report. The initial edge and optical crest need not coincide exactly. Thin liquid films, overlapping fringes, unusual grain shapes and inaccurate normals can remain ambiguous, so inspect the exception overlay.""", "03_measured_wall_profiles"),
    ("5. Prepare the subtraction without replacing raw evidence", """The signed difference is D=R&B−αB, with α initially 1.1 and adjustable. It reduces blue-associated signal and helps identify red-associated fringe evidence, but is not a pure interface image. Changes in exposure, gain or blue intensity alter the required scale. Negative values are retained for analysis rather than silently clipped. Saturated R&B pixels cannot reveal their original crest height, and subtraction can create misleading dark basins in bright water. For that reason raw R&B defines the optical basins and raw B constrains liquid; difference contrast is supplementary. The configured image pair has no independent red-only exposure or calibrated spectral unmixing, so D should not be interpreted as a quantitative red-channel measurement.""", None),
    ("6. Understand the two-sided gradient", """For image intensity I, ∇I points from darker to brighter pixels; −∇I points from brighter to darker pixels. A thin bright water–air fringe surrounded by dark liquid on one side and dark air on the other produces gradients on both flanks. A gradient arrow therefore cannot independently identify the air side. The measured example shows −∇B after mild smoothing, with normalized arrow lengths and low gradients omitted. No arrow direction is forced toward a proposed air interior. The companion curve is explicitly a conceptual bright-line profile. The useful geometric object is the fringe centerline and its connection to the initial grain wall. A free fringe joining a lit air–PDMS wall branch can anchor the surrounding air sector, while the opposite sector is a liquid candidate. Such a join may have a contact angle; straightest continuation is not a sufficient criterion.""", "04_gradients_and_two_sided_fringe"),
    ("7. Infer remaining regions from wall and fringe topology", """Smoothed raw R&B is partitioned by a marker-controlled watershed. H-minima suppress shallow minima; nominal prominence h=8 counts is compared with h=5 and h=12. Markers are kept away from the grain boundary where practical. Each resulting basin is an optical region, not automatically a pore or one physical phase. Valid PDMS-normal profiles are grouped by basin and counted as lit or unlit. Sufficient lit-wall coverage supports an air-side proposal, while broad fluorescence and suitable unlit-wall evidence support liquid. Neighboring basins are connected in a graph when their shared free boundary carries enough R&B and difference crest evidence. Cautious additional liquid votes can then come from the side opposite strongly wall-anchored air neighbors. This is an aggregate wall-side and region-adjacency approximation: it does not explicitly fit every fringe centerline, reconstruct every triple junction, or measure contact angles. Weak regions receive a nominal display label and remain flagged; completely unresolved regions display as liquid while retaining low support and review status. That fallback is not experimental confirmation of water.""", "05_wall_junction_principle"),
    ("7a. Inspect the actual optical partition", """The local diagnostic separates four things that are easy to confuse: measured R&B brightness, the optical watershed basins, sampled wall evidence, and the resulting nominal phase. Random-looking basin colors are identifiers and carry no liquid or air meaning. Cyan lines in the raw-image panel show actual watershed boundaries, which can follow a bright crest but can also split or merge a physical capsule. Yellow and green points indicate lit and unlit initial-wall samples. The blue fixed-water overlay is imposed before inference. A physical phase can span several basins, while one basin can contain both a fixed liquid pixel and an inferred air pixel. This is why the final mask is not simply one phase label per watershed region.""", "10_local_watershed_and_wall_evidence"),
    ("8. Keep possible liquid inclusions and inspect examples", """Fixed fluorescent liquid is preserved even when surrounded by an air proposal, allowing holes within a capsule rather than a compulsory solid air blob. An additional conservative detector proposes compact bright liquid cores using local B and R&B prominence, separation from PDMS and an air-dominated surrounding ring. Those extra cores are flagged because an optical peak alone does not establish composition. This detector estimates bright cores, not complete droplet perimeters, and it may miss weak or dark inclusions. The example panels show raw B, raw R&B and nominal phase on the same source coordinates so residual water, free fringes and grain contacts can be compared. Apparent enclosure is only two-dimensional: the image cannot prove that a droplet is disconnected from liquid above, below or outside the focal plane.""", "07_local_B_RB_phase_examples"),
    ("9. Interpret uncertainty and sensitivity", """Review status can arise from weak or conflicting wall evidence, parameter-dependent labels, clipping, a questionable grain fringe, inclusion-only evidence or other configured warnings. Consult the saved flag definitions: bit values are combinable, not mutually exclusive phase codes. The support image is a heuristic evidence index, not a probability of correctness. The nominal phase map uses 0 outside, 1 PDMS, 2 liquid and 3 air; the separate review map can replace uncertain pore labels by 4. Compare these maps before using a global fraction. Watershed agreement measures stability to selected h values; it cannot detect a shared systematic error. Fluorescence cutoffs such as 45/55/65 counts test how much of the liquid labeling depends on the brightness choice. Look for local changes around thin films and shiny edges as well as changes in the global fraction. Do not treat a stable fraction as proof that the local interfaces are correct.""", "08_uncertainty_and_support"),
    ("10. Run, tune and preserve the evidence", """Select a Spyder Python environment containing the packages in requirements.txt. Edit settings.py, then open RUN_ANALYSIS.py and press F5. The launcher reloads settings for each run, and saves a dated result folder. The supplied input suffixes select the unique TIFFs ending in 0001, 6544 and 6545; complete filenames are also accepted. Start with the supplied configuration and inspect input alignment, grain overlays, water seeds and withheld shiny fringes before adjusting phase inference. If the data are brighter or dimmer, tune the fluorescence cutoff in original units and rerun its sensitivity comparison. If grain boundaries are wrong, repair the initial mask rather than compensating by changing air thresholds. Change blue scale only with attention to exposure and illumination. Retain TIFF masks, raw-value difference, JSON settings and summary, flag definitions, coordinates and this illustrated report together. Original TIFFs are read only. All figure paths in the HTML are relative, so the result folder can be copied or zipped as a unit. This run is a reproducible exploratory estimate; manual review or independent phase reference data are needed for validation.""", "09_parameter_sensitivity"),
]


def _guide_sections(data, cfg):
    """Customize reference and sensitivity descriptions to the run actually shown."""
    sections = list(GUIDE_SECTIONS)
    kind = _reference_kind(data, cfg)
    mode = cfg.get("mask_input", {}).get("mode", "load_binary" if kind == "binary_reference" else "binarize")
    fluorescence_status, watershed_status = _sensitivity_status(data, cfg)
    title, paragraph, stem = sections[0]
    paragraph = paragraph.replace("It uses the initial frame, blue-only fluorescence and the red-plus-blue interface image together.",
                                  "It uses the fixed grain reference, blue-only fluorescence and the red-plus-blue interface image together.")
    sections[0] = title, paragraph, stem
    read = ("This report describes one blue-only (B) and red-plus-blue (R&B) image pair on a common pixel grid. "
            "The fixed grain reference is prepared separately and reused for the series. Set the series limits and even/odd "
            "channel assignment manually in settings.py; the program does not infer illumination from index parity. "
            "Camera images are read without independent normalization, so thresholds stay on the original 0–255 camera-count scale. "
            "Different exposure, gain, illumination or bit depth requires revisiting these settings. B and R&B must describe the "
            "same scene with sufficiently little motion; automatic registration or deformation correction is not performed. "
            "The configured blue_shift_yx applies a manual translation to B. Every reference and pair must already have matching "
            "shape and field of view; masks are never resized to make them fit. Sources and settings are saved for reproducibility.")
    if kind == "binary_reference":
        read += (" The first panel is a synthetic binary geometry reference (dark grains, bright other pixels), "
                 "not a measured camera image. Its brightness conveys geometry only and is never fluorescence evidence.")
    elif kind == "display_reference":
        read += " The first panel is a separately supplied measured display reference; the grain mask is a distinct geometry input."
    else:
        read += " The first panel is the separate initial camera image used for grain binarization."
    sections[1] = ("1. Read the pair and separate reference", read, sections[1][2])
    if mode == "load_binary":
        geometry = ("This run uses an externally supplied binary grain mask. No Phansalkar thresholding, hole filling, "
                    "object removal, dilation, interpolation or resizing is performed on that input. Accepted exact levels are "
                    "0/1, 0/255 or 0/65535; grayscale masks are rejected instead of silently thresholded. "
                    "mask_input.white_is_grain specifies its polarity. An optional external bed mask uses nonzero=inside, "
                    "and is intersected with the selected bed polygon or rectangular ROI; without any bed restriction, the "
                    "whole frame is analysed. Grain pixels outside the combined bed are excluded and counted in the provenance. "
                    "The fixed pore mask is bed minus grains. A nonzero initial_shift_yx is rejected in this mode: provide a "
                    "mask already aligned to the series. Verify the purple boundaries against camera landmarks when a measured "
                    "display reference is available; the synthetic geometry panel alone cannot validate alignment. An incorrect "
                    "mask propagates to all wall-normal and phase estimates. Distances are in pixels, without a physical length calibration.")
        sections[2] = ("2. Load the fixed binary PDMS reference and bed", geometry, sections[2][2])
    else:
        title, paragraph, stem = sections[2]
        paragraph += (" The binarization image may be stored separately from the B/R&B series. The reference is prepared "
                      "once, then kept fixed across pairs. Optional initial_shift_yx translates this raw initial image before "
                      "binarization. An optional external bed mask intersects the polygon/ROI restriction. Any separately supplied "
                      "display reference is shown on the analysis grid without an additional shift.")
        sections[2] = title, paragraph, stem
    title, paragraph, stem = sections[7]
    paragraph = paragraph.replace("nominal prominence h=8 counts is compared with h=5 and h=12.",
                                  f"the nominal prominence is h={float(cfg.get('segmentation', {}).get('watershed_h', 8)):g} counts.")
    paragraph += " " + watershed_status
    sections[7] = title, paragraph, stem
    title, paragraph, stem = sections[10]
    paragraph = ("Review status can arise from weak or conflicting wall evidence, parameter-dependent labels when the relevant "
                 "sensitivity tests were run, clipping, a questionable grain fringe or inclusion-only evidence. Flag bits can overlap; "
                 "consult the saved definitions. The support image is a heuristic evidence index, not a probability of correctness. "
                 "The nominal phase map uses 0 outside, 1 PDMS, 2 liquid and 3 air; the review map can replace uncertain pore labels "
                 "by 4. Compare both maps before using a global fraction. " + fluorescence_status + " " + watershed_status +
                 " Look for local changes around thin films and shiny edges as well as global fraction changes. Parameter agreement "
                 "cannot detect a shared systematic error and does not prove that local interfaces are correct.")
    sections[10] = title, paragraph, stem
    title, paragraph, stem = sections[11]
    paragraph = ("Select a Spyder environment containing requirements.txt, edit settings.py, then open RUN_ANALYSIS.py and press F5. "
                 "Choose the separate reference mode, series index range and manual channel parity before processing new data. "
                 "The fixed bed/pore geometry is reused across pairs, and each pair receives its own phase results and recorded saturation "
                 "fractions. Inspect alignment, fixed grain boundaries, fluorescent-water seeds and withheld shiny fringes before adjusting "
                 "phase inference. Change intensity cutoffs only in the original camera-count units, and enable additional parameter "
                 "sensitivity runs when needed. Correct a bad grain mask at its source. Retain TIFF masks, signed differences, settings, "
                 "pair and series summaries, source provenance and review flags together. Original inputs are read only, and report "
                 "figure paths are relative so the result folder remains portable. These are reproducible exploratory estimates; manual "
                 "review or independent phase references are needed for validation.")
    sections[11] = title, paragraph, stem
    return sections


def _guide(data, cfg, out_dir, figures):
    pore = np.asarray(data["pore"], bool)
    phase = np.asarray(data["phase"])
    review = np.asarray(data.get("review_phase", phase))
    n = int(pore.sum())
    water = int(np.sum((phase == 2) & pore))
    air = int(np.sum((phase == 3) & pore))
    uncertain = int(np.sum((review == 4) & pore))
    seeds = int(np.asarray(data["direct_water"], bool).sum())
    overview = (f"Image: {phase.shape[1]} × {phase.shape[0]} pixels. Analysed pore pixels: {n:,}. "
                f"Nominal liquid: {100*water/max(n,1):.2f}%; nominal air: {100*air/max(n,1):.2f}%. "
                f"Fixed fluorescence-water pixels: {100*seeds/max(n,1):.2f}% of pores. "
                f"Review-map unresolved/flagged pixels: {100*uncertain/max(n,1):.2f}% of pores. "
                "Nominal fractions include weak proposals and display fallbacks; use the review map with them.")
    sections = _guide_sections(data, cfg)
    reference_note = _reference_labels(data, cfg)[2]
    md = ["# Fluorescence-first liquid–air analysis\n", overview + "\n", reference_note + "\n"]
    body = ["<h1>Fluorescence-first liquid–air analysis</h1>",
            "<p class='lead'>A reproducible local workflow with measured examples and explicit uncertainty.</p>",
            "<div class='run'><strong>This run.</strong> " + html.escape(overview) + "</div>",
            "<p>" + html.escape(reference_note) + "</p>",
            "<p class='nav'>" + " · ".join(f"<a href='#s{i}'>{html.escape(title.split('.')[0])}</a>"
                                          for i,(title,_,_) in enumerate(sections)) + "</p>"]
    for i, (title, paragraph, stem) in enumerate(sections):
        md.extend(["## " + title + "\n", paragraph + "\n"])
        body.extend([f"<section id='s{i}'><h2>{html.escape(title)}</h2>",
                     "<p>" + html.escape(paragraph) + "</p>"])
        if stem and stem in figures:
            path = figures[stem]
            md.append(f"![{title}]({path})\n")
            body.append(f"<figure><a href='{html.escape(path, quote=True)}'><img src='{html.escape(path, quote=True)}' alt='{html.escape(title, quote=True)}'></a>"
                        f"<figcaption>{html.escape(title)}. Open the figure for its full resolution.</figcaption></figure>")
        body.append("</section>")
    config_json = json.dumps(_plain(cfg), indent=2, ensure_ascii=False)
    summary_json = json.dumps(_plain(data.get("summary", {})), indent=2, ensure_ascii=False)
    md.extend(["## Main settings and units\n", "| Setting in settings.py | Units | Function |",
               "| --- | --- | --- |"])
    body.append("<h2>Main settings and units</h2><table><thead><tr><th>Setting in settings.py</th><th>Units</th><th>Function</th></tr></thead><tbody>")
    for key, units, meaning in PARAMETER_NOTES:
        md.append(f"| `{key}` | {units} | {meaning} |")
        body.append("<tr><td><code>"+html.escape(key)+"</code></td><td>"+html.escape(units)+"</td><td>"+html.escape(meaning)+"</td></tr>")
    body.append("</tbody></table>")
    md.append("\n")
    md.extend(["## Settings recorded for this run\n", "```json\n" + config_json + "\n```\n",
               "## Summary recorded for this run\n", "```json\n" + summary_json + "\n```\n"])
    body.extend(["<h2>Settings recorded for this run</h2><details><summary>Expand complete configuration</summary><pre>" + html.escape(config_json) + "</pre></details>",
                 "<h2>Summary recorded for this run</h2><details><summary>Expand numerical summary</summary><pre>" + html.escape(summary_json) + "</pre></details>"])
    css = """body{font-family:Segoe UI,Arial,sans-serif;color:#243344;background:#f4f6f8;line-height:1.65;margin:0}main{max-width:1080px;margin:auto;background:white;padding:35px 45px}h1{line-height:1.2;font-size:34px}h2{margin-top:38px;line-height:1.3;font-size:24px}p{font-size:16px}.lead{font-size:19px;color:#58697a}.run{background:#eaf3fa;border-left:5px solid #3199cf;padding:16px;font-size:15px}.nav{font-size:13px}a{color:#186b9f}figure{margin:25px 0 32px}img{width:100%;height:auto;border:1px solid #dce2e8}figcaption{font-size:13px;color:#667480}pre{overflow:auto;background:#f1f3f6;padding:16px;font-size:12px}details{padding:12px;background:#f7f8fa}summary{cursor:pointer}table{width:100%;border-collapse:collapse;font-size:13px}td,th{padding:9px;border:1px solid #dce2e8;text-align:left;vertical-align:top}th{background:#eaf3fa}code{overflow-wrap:anywhere}@media(max-width:700px){main{padding:20px}h1{font-size:28px}}@media print{body{background:white}main{padding:0}section{break-inside:avoid}details{display:none}.nav{display:none}}"""
    document = "<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Liquid–air analysis — illustrated guide</title><style>" + css + "</style></head><body><main>" + "\n".join(body) + "</main></body></html>"
    (out_dir / "GUIDE.html").write_text(document, encoding="utf-8")
    (out_dir / "METHOD_GUIDE.md").write_text("\n".join(md), encoding="utf-8")


def write_report(data: dict, cfg: dict, out_dir: Path) -> dict:
    """Save illustrated HTML/Markdown guides and ten explanatory figure PNGs.

    Required arrays follow the pipeline's original pixel grid. ``profiles``
    contains N×2 ``yx`` and outward ``normals_xy`` plus scalar profile fields.
    Figures preserve source coordinates; display stretches never alter analysis.
    Return paths are relative to ``out_dir`` so the folder remains portable.
    """
    out_dir = Path(out_dir)
    directory = out_dir / "figures"
    directory.mkdir(parents=True, exist_ok=True)
    cfg = dict(cfg or {})
    initial, b, rb = (np.asarray(data[key]) for key in ("initial", "B", "RB"))
    grain, pore = (np.asarray(data[key], bool) for key in ("grain", "pore"))
    examples = _examples(cfg, b.shape)
    figures = {}
    plt.rcParams.update({"font.family":"DejaVu Sans", "axes.spines.top":False,
                         "axes.spines.right":False, "font.size":10})

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.1))
    reference_title, overlay_title, reference_note = _reference_labels(data, cfg)
    lo, hi = (0, 255) if _reference_kind(data, cfg) == "binary_reference" else np.percentile(initial, [2, 99])
    for ax, array, title, limits in zip(axes, [initial, b, rb],
        [reference_title, "Blue only / fluorescence", "Red + blue / fluorescence and fringes"],
        [(lo, hi), (6, 120), (15, 255)]):
        _image(ax, array, title, limits)
        ax.set_xlabel("Source x (pixels)"); ax.set_ylabel("Source y (pixels)")
        for i, example in enumerate(examples):
            x0,y0,x1,y1=example["xyxy"]
            ax.add_patch(Rectangle((x0,y0),x1-x0,y1-y0,fill=False,color="#f1ba4d",lw=.9))
            ax.text(x0,y0,str(i+1),color="#ffe58d",fontsize=9,va="bottom")
    fig.text(.5,.01,reference_note + " B 6–120; R&B 15–255 camera counts. Rectangles mark examples.",ha="center",fontsize=8.5)
    fig.tight_layout(rect=(0,.035,1,1))
    figures["01_input_images"]=_save(fig,directory,"01_input_images")

    fig, axes=plt.subplots(1,3,figsize=(15,5.2))
    _image(axes[0],initial,overlay_title,(lo,hi))
    _overlay(axes[0],grain,"#ad72b8",.65)
    _image(axes[1],b,"Fixed fluorescent water and withheld shiny fringe",(6,120))
    _overlay(axes[1],np.asarray(data["direct_water"],bool),COLORS[2],.75)
    _overlay(axes[1],np.asarray(data["shiny_exception"],bool),"#e25fbb",.9)
    difference=np.asarray(data.get("difference",rb-float(cfg.get("blue_scale",1.1))*b),float)
    scale=max(1.,float(np.percentile(np.abs(difference[pore]),98))) if pore.any() else 1.
    im=_image(axes[2],difference,f"Signed difference: R&B − {float(cfg.get('blue_scale',1.1)):g} B",(-scale,scale),cmap="RdBu_r")
    fig.colorbar(im,ax=axes[2],shrink=.6,label="Camera-count difference")
    fig.legend(handles=[Patch(color="#ad72b8",label="Fixed PDMS"),Patch(color=COLORS[2],label="Fixed water seed"),Patch(color="#e25fbb",label="Possible shiny grain fringe")],loc="lower center",ncol=3,frameon=False,fontsize=9)
    fig.tight_layout(rect=(0,.06,1,1))
    figures["02_geometry_fluorescence_and_difference"]=_save(fig,directory,"02_geometry_fluorescence_and_difference")
    figures["03_measured_wall_profiles"]=_profile_figure(data,cfg,directory)
    figures["04_gradients_and_two_sided_fringe"]=_gradient_figure(data,examples,directory,cfg)
    figures["05_wall_junction_principle"]=_junction_figure(directory)

    phase=np.asarray(data["phase"])
    review=np.asarray(data.get("review_phase",phase))
    fig,axes=plt.subplots(1,2,figsize=(13.8,6.2))
    _phase(axes[0],phase,"Nominal liquid–air distribution")
    _phase(axes[1],review,"Review map: uncertain pore pixels shown separately")
    _legend(fig,review=True,y=.015)
    fig.text(.5,.07,"Nominal labels include weak proposals and unresolved liquid display fallbacks. Inspect the review map before interpreting phase fractions.",ha="center",fontsize=9)
    fig.tight_layout(rect=(0,.10,1,1))
    figures["06_whole_image_phase_and_review"]=_save(fig,directory,"06_whole_image_phase_and_review")

    fig,axes=plt.subplots(len(examples),3,figsize=(12.5,4.1*len(examples)),squeeze=False)
    for row,example in enumerate(examples):
        box=example["xyxy"]; extent=(box[0],box[2],box[3],box[1])
        _image(axes[row,0],_crop(b,box),"B: "+example["name"],(6,120),extent)
        _image(axes[row,1],_crop(rb,box),"R&B",(15,255),extent)
        _phase(axes[row,2],_crop(phase,box),"Nominal phase",extent)
        for ax in axes[row]:
            _grain_contour(ax,grain,box)
            ax.set_xlim(box[0],box[2]);ax.set_ylim(box[3],box[1])
        for col in range(3):
            axes[row,col].set_xlabel("Source x (pixels)",fontsize=9)
            axes[row,col].set_ylabel("Source y (pixels)",fontsize=9)
    _legend(fig,y=.005)
    fig.tight_layout(rect=(0,.035,1,1))
    figures["07_local_B_RB_phase_examples"]=_save(fig,directory,"07_local_B_RB_phase_examples")

    fig,axes=plt.subplots(1,2,figsize=(13.8,6.2))
    support=np.asarray(data.get("support",np.zeros(phase.shape)),float)
    shown=np.ma.masked_where(~pore,support)
    im=_image(axes[0],shown,"Heuristic evidence support (not probability)",(0,1),cmap="viridis")
    fig.colorbar(im,ax=axes[0],shrink=.7,label="Support index")
    flags=np.asarray(data.get("review_flags",review==4))
    any_flag=(flags!=0)&pore
    _image(axes[1],b,"Any review flag over measured B",(6,120))
    _overlay(axes[1],any_flag,COLORS[4],.75)
    fig.text(.5,.02,"Review flags can overlap. Use the saved definitions to interpret individual bits; parameter stability is not independent phase validation.",ha="center",fontsize=9)
    fig.tight_layout(rect=(0,.05,1,1))
    figures["08_uncertainty_and_support"]=_save(fig,directory,"08_uncertainty_and_support")
    figures["09_parameter_sensitivity"]=_sensitivity_figure(data,cfg,directory)
    figures["10_local_watershed_and_wall_evidence"]=_basin_figure(data,examples,directory)
    _guide(data,cfg,out_dir,figures)
    return {"files":list(figures.values())+["GUIDE.html","METHOD_GUIDE.md"],
            "html":"GUIDE.html","markdown":"METHOD_GUIDE.md","figures":figures}

"""Fluorescence-first liquid/air proposals on the original camera grid.

This is a portable implementation of the exploratory wall/fringe method. It
does not reconstruct contact angles or explicitly trace every triple junction.
Positive fluorescence is a one-sided liquid constraint; a dark B pixel is never
positive air evidence. Keep the accompanying review mask when using the result.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy import ndimage as ndi
from skimage.morphology import h_minima
from skimage.segmentation import watershed


def _say(progress: Callable | None, message: str) -> None:
    if progress is not None:
        progress(message)


def _seeded_watershed(surface, pore, distance, prominence, clearance):
    """Find optical basins, giving each disconnected pore a usable marker."""
    minima = h_minima(surface, prominence) & pore & (distance > clearance)
    markers, number = ndi.label(minima)
    components, count = ndi.label(pore)
    have_marker = np.zeros(count + 1, dtype=bool)
    have_marker[np.unique(components[markers > 0])] = True
    slices = ndi.find_objects(components)
    for index in np.flatnonzero(~have_marker[1:]) + 1:
        sl = slices[index - 1]
        if sl is None:
            continue
        local = components[sl] == index
        values = np.where(local, surface[sl], np.inf)
        iy, ix = np.unravel_index(np.argmin(values), values.shape)
        number += 1
        markers[iy + sl[0].start, ix + sl[1].start] = number
    del components
    return watershed(surface, markers, mask=pore).astype(np.int32)


def _mean_by_label(labels, values, size, mask=None):
    if mask is None:
        ids, vals = labels.ravel(), values.ravel()
    else:
        ids, vals = labels[mask], values[mask]
    count = np.bincount(ids, minlength=size)
    total = np.bincount(ids, weights=vals, minlength=size)
    return total / np.maximum(count, 1), count


def _free_fringe_graph(labels, rb, difference, distance, rb_core, d_core, cfg):
    """Pair adjacent raw-image basins across a sufficiently positive free crest.

    A graph edge has no intrinsic air side. The side must already be anchored
    by lit initial PDMS walls before an opposite-side liquid vote is allowed.
    """
    number = int(labels.max()) + 1
    accumulated = {}
    clearance = float(cfg.get("graph_wall_clearance_px", 3.0))
    rb_threshold = float(cfg.get("graph_RB_crest_counts", 12.0))
    d_threshold = float(cfg.get("graph_difference_crest_counts", 5.0))
    for axis in (0, 1):
        if axis == 0:
            a, b = labels[:-1, :], labels[1:, :]
            r = np.maximum(rb[:-1, :], rb[1:, :])
            d = np.maximum(difference[:-1, :], difference[1:, :])
            clear = np.minimum(distance[:-1, :], distance[1:, :])
        else:
            a, b = labels[:, :-1], labels[:, 1:]
            r = np.maximum(rb[:, :-1], rb[:, 1:])
            d = np.maximum(difference[:, :-1], difference[:, 1:])
            clear = np.minimum(distance[:, :-1], distance[:, 1:])
        use = (a > 0) & (b > 0) & (a != b) & (clear > clearance)
        if not np.any(use):
            continue
        low = np.minimum(a[use], b[use])
        high = np.maximum(a[use], b[use])
        keys = low.astype(np.int64) * number + high
        crest = (r[use] - np.maximum(rb_core[low], rb_core[high]) > rb_threshold)
        crest &= d[use] - np.maximum(d_core[low], d_core[high]) > d_threshold
        unique, inverse, count = np.unique(keys, return_inverse=True, return_counts=True)
        supported = np.bincount(inverse, weights=crest, minlength=len(unique))
        for key, n, ns in zip(unique, count, supported):
            key = int(key)
            if key in accumulated:
                accumulated[key][0] += int(n)
                accumulated[key][1] += float(ns)
            else:
                accumulated[key] = [int(n), float(ns)]
    minimum = int(cfg.get("graph_min_boundary_samples", 4))
    return [
        {"region_a": int(key // number), "region_b": int(key % number),
         "samples": int(count), "crest_support": float(supported / count)}
        for key, (count, supported) in accumulated.items() if count >= minimum
    ]


def _classify_regions(labels, B, RB, difference, distance, pore, profiles,
                      direct_water, cfg):
    """Estimate nonfluorescent pixels from wall-side and free-fringe adjacency."""
    size = int(labels.max()) + 1
    area = np.bincount(labels.ravel(), minlength=size)
    yx = np.asarray(profiles["yx"], dtype=np.intp)
    valid = np.asarray(profiles["valid"], dtype=bool).copy()
    if yx.ndim != 2 or yx.shape[1] != 2 or len(yx) != len(valid):
        raise ValueError("profiles['yx'] must have shape (n,2), matching valid/lit/unlit.")
    if len(yx):
        wall_ids = labels[yx[:, 0], yx[:, 1]]
    else:
        wall_ids = np.zeros(0, dtype=np.int32)
    valid &= wall_ids > 0
    lit = np.asarray(profiles["lit"], dtype=bool) & valid
    unlit = np.asarray(profiles["unlit"], dtype=bool) & valid
    count = np.bincount(wall_ids[valid], minlength=size)
    nlit = np.bincount(wall_ids[lit], minlength=size)
    nunlit = np.bincount(wall_ids[unlit], minlength=size)
    lit_fraction = nlit / np.maximum(count, 1)
    unlit_fraction = nunlit / np.maximum(count, 1)
    water_pixels = np.bincount(labels[direct_water], minlength=size)
    water_fraction = water_pixels / np.maximum(area, 1)
    core = pore & (distance > float(cfg.get("region_core_clearance_px", 3.0)))
    rb_core, core_count = _mean_by_label(labels, RB, size, core)
    d_core, _ = _mean_by_label(labels, difference, size, core)
    whole_rb, _ = _mean_by_label(labels, RB, size, pore)
    whole_d, _ = _mean_by_label(labels, difference, size, pore)
    rb_core = np.where(core_count >= 4, rb_core, whole_rb)
    d_core = np.where(core_count >= 4, d_core, whole_d)

    air = ((lit_fraction > float(cfg.get("air_wall_lit_fraction", .27)))
           & (nlit >= int(cfg.get("air_wall_min_lit_samples", 12)))
           & (count >= int(cfg.get("air_wall_min_valid_samples", 20)))
           & (unlit_fraction < float(cfg.get("air_wall_max_unlit_fraction", .65))))
    strong_air = (air & (lit_fraction > float(cfg.get("strong_air_lit_fraction", .43)))
                  & (nlit >= int(cfg.get("strong_air_min_lit_samples", 20)))
                  & (count >= int(cfg.get("strong_air_min_valid_samples", 30)))
                  & (unlit_fraction < float(cfg.get("strong_air_max_unlit_fraction", .40))))

    # These are geometric/optical region hypotheses. A single bright pixel does
    # not relabel an entire basin, and the pixelwise water constraint below wins.
    accepted = np.zeros(size, dtype=np.uint8)
    accepted[air] = 3
    support = np.zeros(size, dtype=np.float32)
    support[air] = np.clip(.40 + .60 * lit_fraction[air], 0, .95)
    broad_water = water_fraction >= float(cfg.get("broad_water_region_fraction", .60))
    # A substantial region interior already known to be water supports extending
    # liquid locally. This extension remains tentative unless wall evidence agrees.
    broad_water &= water_pixels >= int(cfg.get("broad_water_min_pixels", 20))
    accepted[broad_water] = 2
    support[broad_water] = np.clip(.45 + .40 * water_fraction[broad_water], 0, .85)
    strong_air[broad_water] = False
    reason = np.full(size, "unresolved", dtype=object)
    reason[air] = "lit_PDMS_wall_side"
    reason[broad_water] = "substantial_fluorescence_within_basin"

    edges = []
    inferred = np.zeros(size, dtype=bool)
    if cfg.get("use_fringe_graph", True):
        edges = _free_fringe_graph(labels, RB, difference, distance, rb_core, d_core, cfg)
        votes = np.zeros(size, dtype=np.float32)
        neighbours = np.zeros(size, dtype=np.int32)
        required = float(cfg.get("graph_min_crest_support", .50))
        for edge in edges:
            if edge["crest_support"] < required:
                continue
            i, j = edge["region_a"], edge["region_b"]
            for src, dst in ((i, j), (j, i)):
                if (strong_air[src] and accepted[dst] == 0 and count[dst] >= 16
                        and lit_fraction[dst] < .25 and unlit_fraction[dst] > .40
                        and area[dst] >= 10):
                    votes[dst] += edge["samples"] * edge["crest_support"]
                    neighbours[dst] += 1
        inferred = ((votes >= float(cfg.get("graph_min_votes", 6)))
                    & (neighbours >= int(cfg.get("graph_min_air_neighbours", 2)))
                    & (accepted == 0))
        accepted[inferred] = 2
        support[inferred] = np.clip(.44 + .03 * np.minimum(neighbours[inferred], 4), 0, .60)
        reason[inferred] = "free_fringe_opposite_wall_anchored_air"

    nominal = accepted.copy()
    unlit_liquid = ((nominal == 0) & (count >= 12)
                    & (unlit_fraction > .50) & (unlit_fraction > lit_fraction))
    nominal[unlit_liquid] = 2
    support[unlit_liquid] = .25
    reason[unlit_liquid] = "tentative_unlit_wall_liquid"
    # A weak lit-wall air proposal reproduces the exploratory nominal map, but
    # is never promoted to accepted air and is always exposed in the review mask.
    weak_air = ((nominal == 0) & (count >= 20) & (lit_fraction > .18) & (nlit >= 10))
    nominal[weak_air] = 3
    support[weak_air] = .22
    reason[weak_air] = "tentative_lit_wall_air"
    unresolved = nominal == 0
    nominal[unresolved] = 2
    support[unresolved] = .10
    reason[unresolved] = "liquid_display_fallback_unresolved"
    # Index zero is outside a basin, never a physical phase inference.
    accepted[0] = nominal[0] = 0
    support[0] = 0
    reason[0] = "outside"
    return {"accepted": accepted, "nominal": nominal, "support": support,
            "area": area, "wall_samples": count, "lit_samples": nlit,
            "unlit_samples": nunlit, "lit_fraction": lit_fraction,
            "unlit_fraction": unlit_fraction, "water_fraction": water_fraction,
            "direct_water_pixels": water_pixels, "strong_air": strong_air,
            "graph_inferred": inferred, "reason": reason, "edges": edges}


def _liquid_inclusions(nominal, pore, distance, highpass_b, highpass_rb, cfg):
    """Conservative compact bright liquid cores surrounded by nominal air.

    This is independent of user drawings. Optical cores are not complete droplet
    perimeters, and a positive local peak is not proof of water on its own.
    """
    air = (nominal == 3) & pore
    cores = np.zeros(pore.shape, dtype=bool)
    records = []
    if not cfg.get("detect_liquid_inclusions", True) or not np.any(air):
        return cores, records
    bright = ((highpass_b > float(cfg.get("inclusion_B_prominence_counts", 13)))
              & (highpass_rb > float(cfg.get("inclusion_RB_prominence_counts", 17)))
              & pore)
    blobs, _ = ndi.label(bright)
    for i, sl in enumerate(ndi.find_objects(blobs), 1):
        if sl is None:
            continue
        local = blobs[sl] == i
        area = int(local.sum())
        if not int(cfg.get("inclusion_min_area_px", 8)) <= area <= int(cfg.get("inclusion_max_area_px", 450)):
            continue
        height, width = local.shape
        if max(height, width) > int(cfg.get("inclusion_max_extent_px", 40)):
            continue
        y0, x0 = sl[0].start, sl[1].start
        y1, x1 = sl[0].stop, sl[1].stop
        # A cropped component cannot establish a closed inclusion in this image.
        if y0 < 4 or x0 < 4 or y1 > pore.shape[0] - 4 or x1 > pore.shape[1] - 4:
            continue
        context = (slice(max(y0 - 5, 0), min(y1 + 5, pore.shape[0])),
                   slice(max(x0 - 5, 0), min(x1 + 5, pore.shape[1])))
        cmask = blobs[context] == i
        ring = ndi.binary_dilation(cmask, iterations=4) & ~ndi.binary_dilation(cmask)
        if not np.any(ring) or pore[context][ring].mean() < .90:
            continue
        ring &= pore[context]
        if not np.any(ring):
            continue
        surround = float(air[context][ring].mean())
        median_clearance = float(np.median(distance[sl][local]))
        near_fraction = float((distance[sl][local] < 2).mean())
        if (surround < float(cfg.get("inclusion_air_surround_fraction", .65))
                or median_clearance < float(cfg.get("inclusion_median_grain_clearance_px", 6.0))
                or near_fraction >= .20):
            continue
        yy, xx = np.where(local)
        cores[sl][local] = True
        records.append({"id": len(records) + 1, "area_px": area,
                        "center_xy": [float(xx.mean() + x0), float(yy.mean() + y0)],
                        "air_surround_fraction": surround,
                        "median_PDMS_clearance_px": median_clearance,
                        "near_PDMS_fraction": near_fraction})
    return cores, records


def infer_phases(initial, B, RB, grain, pore, profiles, direct_water, cfg,
                 progress=None):
    """Return nominal phases, evidence support and reproducibility warnings.

    Parameters are ordinary NumPy arrays; ``cfg`` is a dict. Fluorescent water
    must have been identified before this function is called. Those fixed pixels
    remain liquid during each hypothesis and sensitivity run. The raw B/R&B
    images still determine all optical basins and wall/fringe geometry.

    Codes: 0 outside, 1 initial PDMS, 2 liquid (dark or fluorescent), 3 air.
    ``weak`` exposes nominal guesses; ``stable`` compares requested watershed
    prominence settings, and is not statistical confidence. ``support`` is a
    heuristic 0--1 index, not a probability. Initial grain geometry is immutable.
    """
    supplied_cfg = dict(cfg or {})
    cfg = dict(supplied_cfg.get("segmentation", supplied_cfg))
    if "blue_scale" not in cfg and "blue_scale" in supplied_cfg:
        cfg["blue_scale"] = supplied_cfg["blue_scale"]
    B = np.asarray(B, dtype=np.float32)
    RB = np.asarray(RB, dtype=np.float32)
    grain = np.asarray(grain, dtype=bool)
    pore = np.asarray(pore, dtype=bool)
    direct_water = np.asarray(direct_water, dtype=bool)
    if any(a.shape != B.shape for a in (initial, RB, grain, pore, direct_water)):
        raise ValueError("All input images and masks must have identical shapes.")
    if B.ndim != 2 or not np.all(np.isfinite(B)) or not np.all(np.isfinite(RB)):
        raise ValueError("B and R&B must be finite two-dimensional grayscale arrays.")
    if np.any(grain & pore) or np.any(direct_water & ~pore):
        raise ValueError("Grain and pore must be disjoint; direct_water must lie in pore space.")
    if not np.any(pore):
        raise ValueError("The analysed pore mask is empty; check the initial-frame mask.")
    alpha = float(cfg.get("blue_scale", 1.1))
    sigma = float(cfg.get("watershed_sigma_px", 1.0))
    nominal_h = float(cfg.get("watershed_h", 8.0))
    settings = cfg.get("watershed_sensitivity", [5.0, 8.0, 12.0])
    if isinstance(settings, bool):
        settings = [5.0, nominal_h, 12.0] if settings else [nominal_h]
    settings = [float(h) for h in settings]
    settings = [nominal_h] + [h for h in dict.fromkeys(settings) if h != nominal_h]
    if sigma < 0 or any(h <= 0 for h in settings):
        raise ValueError("Smoothing must be nonnegative; watershed prominence must be positive.")
    _say(progress, "Preparing raw-image basins; fluorescence water is already fixed.")
    distance = ndi.distance_transform_edt(~grain).astype(np.float32)
    smooth_rb = ndi.gaussian_filter(RB, sigma)
    # High outside level prevents external background minima from becoming seeds.
    surface = smooth_rb.copy()
    surface[~pore] = float(np.max(smooth_rb)) + max(settings) + 1.0
    difference = RB - alpha * B
    highpass_b = ndi.gaussian_filter(B, .8) - ndi.gaussian_filter(B, 8.0)
    highpass_rb = ndi.gaussian_filter(RB, .8) - ndi.gaussian_filter(RB, 8.0)
    stable = pore.copy()
    metrics = []
    result = None
    for run, h in enumerate(settings):
        _say(progress, f"Watershed prominence {h:g}: wall-side and free-fringe inference.")
        # Known fluorescence is established before watershed/region inference.
        # None of the subsequent region votes can write into these fixed seeds.
        phase = np.zeros(B.shape, dtype=np.uint8)
        phase[grain] = 1
        phase[direct_water] = 2
        unassigned = pore & ~direct_water
        labels = _seeded_watershed(surface, pore, distance, h,
                                   float(cfg.get("watershed_marker_clearance_px", 2.0)))
        votes = _classify_regions(labels, B, RB, difference, distance, pore,
                                  profiles, direct_water, cfg)
        # Start from known liquid before assigning the remaining dark/ambiguous
        # pixels, so no wall hypothesis can ever relabel a fluorescence seed.
        phase[unassigned] = votes["nominal"][labels[unassigned]]
        # Inclusions are evaluated against the surrounding wall-derived phase;
        # fixed fluorescence cores remain visible in raw evidence and topology.
        wall_phase = votes["nominal"][labels]
        inclusion, inclusions = _liquid_inclusions(wall_phase, pore, distance,
                                                   highpass_b, highpass_rb, cfg)
        phase[inclusion] = 2
        support = votes["support"][labels]
        support[direct_water] = .95
        support[inclusion & ~direct_water] = .60
        support[~pore] = 0
        weak = pore & ((votes["accepted"][labels] == 0) | (support < .50))
        # A fluorescence-containing basin can have different phases locally.
        # Conflict flags concern only the remaining inferred pixels.
        mixed_conflict = ((votes["water_fraction"] > .10)
                          & (votes["nominal"] == 3))[labels]
        weak |= pore & mixed_conflict & ~direct_water
        weak[direct_water] = False
        weak[inclusion & ~direct_water] = True
        m = {"watershed_h": h, "regions": int(labels.max()),
             "nominal_liquid_fraction_of_pore": float((phase[pore] == 2).mean()),
             "nominal_air_fraction_of_pore": float((phase[pore] == 3).mean()),
             "weak_fraction_of_pore": float(weak[pore].mean()),
             "automatic_inclusion_core_count": len(inclusions),
             "automatic_inclusion_core_pixels": int(inclusion.sum())}
        metrics.append(m)
        if run == 0:
            records = []
            for i in range(1, len(votes["area"])):
                if not votes["area"][i]:
                    continue
                records.append({"id": i, "area_px": int(votes["area"][i]),
                                "nominal_region_phase": int(votes["nominal"][i]),
                                "accepted_region_phase": int(votes["accepted"][i]),
                                "support_index": float(votes["support"][i]),
                                "valid_wall_samples": int(votes["wall_samples"][i]),
                                "lit_wall_samples": int(votes["lit_samples"][i]),
                                "unlit_wall_samples": int(votes["unlit_samples"][i]),
                                "lit_wall_fraction": float(votes["lit_fraction"][i]),
                                "unlit_wall_fraction": float(votes["unlit_fraction"][i]),
                                "direct_water_fraction": float(votes["water_fraction"][i]),
                                "direct_water_pixels": int(votes["direct_water_pixels"][i]),
                                "reason": str(votes["reason"][i]),
                                "graph_inferred_liquid": bool(votes["graph_inferred"][i])})
            result = {"phase": phase, "region_labels": labels, "support": support,
                      "weak": weak, "inclusion_core": inclusion,
                      "records": {"regions": records, "inclusions": inclusions,
                                  "free_fringe_edges": votes["edges"]}}
        else:
            stable &= phase == result["phase"]
        _say(progress, f"  {m['regions']:,} basins; {100*m['nominal_liquid_fraction_of_pore']:.1f}% liquid.")
    stable &= pore
    result["stable"] = stable
    result["records"]["sensitivity"] = metrics
    result["parameters"] = {
        "blue_scale": alpha, "watershed_sigma_px": sigma,
        "nominal_watershed_h": nominal_h, "watershed_sensitivity": settings,
        "watershed_marker_clearance_px": float(cfg.get("watershed_marker_clearance_px", 2.0)),
        "air_wall_lit_fraction": float(cfg.get("air_wall_lit_fraction", .27)),
        "air_wall_min_lit_samples": int(cfg.get("air_wall_min_lit_samples", 12)),
        "air_wall_min_valid_samples": int(cfg.get("air_wall_min_valid_samples", 20)),
        "air_wall_max_unlit_fraction": float(cfg.get("air_wall_max_unlit_fraction", .65)),
        "broad_water_region_fraction": float(cfg.get("broad_water_region_fraction", .60)),
        "use_fringe_graph": bool(cfg.get("use_fringe_graph", True)),
        "graph_min_air_neighbours": int(cfg.get("graph_min_air_neighbours", 2)),
        "detect_liquid_inclusions": bool(cfg.get("detect_liquid_inclusions", True)),
        "fluorescence_priority": "Direct-water pixels are liquid before and throughout every inference run.",
        "basin_scope": "Nominal region phase describes only pixels not constrained by direct fluorescence or inclusion cores.",
        "contact_topology": "Initial-wall-side/region adjacency approximation, not explicit triple-junction reconstruction.",
        "nominal_fallback": "Unresolved pixels display as liquid and remain flagged weak; dark B alone never means air.",
        "inclusion_scope": "Automatic compact optical liquid-core proposals; no manual outlines, no complete droplet reconstruction.",
        "support_meaning": "Heuristic evidence index, not a probability.",
        "sensitivity_meaning": "Phase agreement under watershed prominence changes, not a confidence interval.",
        "stage_order": ["initial grain geometry", "fixed fluorescence water seeds",
                        "raw-image basin geometry", "unknown-pixel wall/fringe inference",
                        "compact liquid inclusion proposals", "sensitivity checks"],
    }
    if not np.all(result["phase"][direct_water] == 2):
        raise RuntimeError("Internal error: a fixed fluorescence water pixel changed phase.")
    if not np.all(np.isin(result["phase"][pore], [2, 3])):
        raise RuntimeError("Internal error: nominal phases do not partition the analysed pores.")
    result["records"]["water_preserved_audit"] = {
        "fixed_fluorescence_water_pixels": int(direct_water.sum()),
        "fixed_water_pixels_changed_from_liquid": int(np.count_nonzero(result["phase"][direct_water] != 2)),
        "fixed_water_phase_stable_across_runs": bool(np.all(stable[direct_water])),
    }
    return result

# Run an alternating B / R&B image series in Spyder

This folder analyses a fixed grain geometry and a sequence of image pairs. Each pair produces a nominal liquid–air distribution and a row of pore-area saturation measurements. Liquid includes fluorescent water and dark diluted or fresh water. Dark B pixels are not automatically air.

The underlying optical method, illustrated wall profiles, two-sided gradients and limitations are explained in [METHOD_GUIDE.md](METHOD_GUIDE.md). This guide covers inputs, pairing, batch controls and saved measurements. The method remains an optical estimate with review flags, rather than a validated phase measurement.

## 1. Open and run

1. Keep `RUN_ANALYSIS.py`, `settings.py`, `phase_analysis/` and the other package files together.
2. Install `requirements.txt` into the environment used by the Spyder kernel if required. See [START_HERE.md](START_HERE.md).
3. Edit `settings.py`. Set the image folder, mask source, inclusive index range and B parity explicitly.
4. Open `RUN_ANALYSIS.py` in Spyder and press **F5**. Inputs are read only. Progress appears in the console.
5. Open the new results directory. Check `pair_plan.csv` before interpreting `saturation.csv`, and inspect the first pair's phase and review maps.

The supplied demo uses indices **6544–6545**, with **odd indices assigned to B**. Thus frame 6545 is B and frame 6544 is R&B. It uses initial frame 0001 as the separate grain reference; that reference is outside the demo's requested pair range.

## 2. Supply the grain geometry independently

The grain reference does not need to be a member of the analysed series. Absolute paths may point to another folder. Relative mask/reference paths are resolved against `input_folder`. Use a full filename when a short numeric suffix could match more than one image.

### A. Binarize a separate raw initial image here

```python
"mask_input": {
    "mode": "binarize",
    "image_path": r"D:\Experiment\reference\initial.tiff",
    "grain_mask_path": None,
    "white_is_grain": True,
    "bed_mask_path": None,
    "display_reference_path": None,
},
```

This runs the configured Phansalkar threshold on the raw reference. Its local window handles slow illumination variation; small-object cleanup and enclosed-grain-hole filling follow. The grain mask is not dilated. Parameters under `grain`, such as `window` and `k`, apply in this mode. The fixed grain geometry is prepared once for the batch.

### B. Use a binary mask made elsewhere

```python
"mask_input": {
    "mode": "load_binary",
    "image_path": None,
    "grain_mask_path": r"D:\Experiment\reference\grain_mask.tiff",
    "white_is_grain": True,
    "bed_mask_path": r"D:\Experiment\reference\bed_mask.tiff",
    "display_reference_path": r"D:\Experiment\reference\initial.tiff",
},
```

Set `white_is_grain=True` when nonzero/white pixels are grains and zero/black pixels are not. Set it to `False` for the opposite polarity. Check this carefully: inverted grain polarity changes the meaning of every phase fraction.

External masks must already be binary grayscale TIFF or single-image PNG, using **0/1, 0/255 or 0/65535**. They must share the series image dimensions and pixel alignment. The program does not threshold, resize, shift, fill or clean these masks. Set `initial_shift_yx=(0, 0)` in this mode; align an external mask beforehand if needed. Grain labels outside the selected bed are excluded from the analysis.

`display_reference_path` is optional and supplies an image for figures only. It does not alter the binary grain mask. If omitted, figures use clearly identified synthetic dark-grain/bright-pore geometry instead of a measured initial image.

### Define the analysed bed

The bed contains the porous experimental area; pore space is **bed minus grains**. The bed boundary is essential when the camera view includes regions outside the experiment.

- Supply an optional binary `bed_mask_path`, with **nonzero = inside the bed**. Its polarity is independent of `white_is_grain`.
- Or set `bed_polygon_xy` using `[x, y]` vertices, or `roi_xyxy=[x0, y0, x1, y1]` for a rectangle. Choose a polygon or rectangle, not both. Rectangle upper bounds are exclusive.
- A supplied bed mask intersects the configured polygon/rectangle. To use the supplied mask alone, set both `bed_polygon_xy` and `roi_xyxy` to `None`.
- If all bed restrictions are absent, the full image is the bed. Do this only when that is appropriate.

The supplied polygon fits the current experiment. Check it for another field of view. `expected_shape_hw` is `(height, width)`, for example `(3000, 4096)`; it is a check, not a resizing instruction. The coordinate origin is the image's upper left, with x to the right and y downward.

## 3. Select the files and illumination order manually

```python
"input_folder": Path(r"D:\Experiment\series"),
"series": {
    "start_index": 0,
    "end_index": 1999,
    "blue_parity": "even",
    "file_pattern": "*.tiff",
    "index_regex": r"(?P<index>\d+)$",
    "unpaired_policy": "error",
    "missing_pair_policy": "error",
    "pair_interval_seconds": None,
},
```

`start_index` and `end_index` are both **inclusive**. Files are paired in consecutive slots anchored at the start: `(start, start+1)`, `(start+2, start+3)`, and so on. Pair numbers start at zero. The program does not infer illumination from brightness.

`blue_parity` refers to the **absolute numeric frame index**, not its position in the selected range:

| Requested range | `blue_parity` | First pair | Last pair | Number of pairs |
|---|---|---|---|---|
| 0–1999 | `"even"` | B=0000, R&B=0001 | B=1998, R&B=1999 | 1000 |
| 0–1999 | `"odd"` | B=0001, R&B=0000 | B=1999, R&B=1998 | 1000 |
| 6544–6545 | `"odd"` | B=6545, R&B=6544 | Same pair | 1 |

For example, a range beginning at 101 still assigns frame 101 to B when parity is `"odd"`; it assigns frame 102 to B when parity is `"even"`.

### The 0–2000 case

The inclusive range 0–2000 contains **2001 images**, leaving frame 2000 without a partner. The default `unpaired_policy="error"` stops before image processing. Either select 0–1999, include frame 2001 and select 0–2001, or explicitly set `unpaired_policy="drop_last"`. With `drop_last`, 0–2000 processes the same 1000 pairs as 0–1999 and records 2000 as unpaired; it does not invent a partner.

### Filenames and missing images

By default the index is the final run of digits in the filename stem. Thus `Basler_..._0007.tiff` becomes index 7; camera serial numbers earlier in the name do not affect it. Both `.tif` and `.tiff` are supported. `file_pattern=None` includes both; a pattern such as `"*.tiff"` deliberately selects only matching filenames. Discovery uses one folder, not its subfolders.

Use `file_pattern` to select one acquisition when a folder contains multiple series, for example `"Basler*20261023*.tiff"`. Duplicate numeric indices are rejected rather than chosen arbitrarily. For another filename convention, change `index_regex` while retaining its named `(?P<index>...)` group.

The default `missing_pair_policy="error"` stops when a requested frame is absent. To continue past incomplete pairs, set `"skip"`. **An entire incomplete pair is skipped, and later pairs never shift.** If frame 0003 is missing, `(0002, 0003)` is omitted; `(0004, 0005)` stays together. A missing B image is never replaced by the next R&B image. Skips remain visible in the plan and summary.

If known, set `pair_interval_seconds` to the time **between successive pairs**, not the B-to-R&B delay. `time_seconds` then starts at zero for the first requested pair and preserves skipped-pair gaps. Leave it `None` if the timing is unknown; frame filenames alone do not provide a calibrated time axis.

## 4. Choose batch speed, exports and recovery

```python
"batch": {
    "resume_from": None,
    "report_policy": "first",
    "sensitivity_policy": "none",
    "export_level": "minimal",
    "pair_error_policy": "stop",
    "plot_every": 10,
    "preflight_shapes": True,
    "retain_last_pair": True,
},
```

| Setting | Practical meaning |
|---|---|
| `resume_from` | `None` creates a new run. Set an existing batch run directory to resume that run. |
| `report_policy` | `"none"`, `"first"` or `"all"`: illustrated per-pair reports for none, the first pair or every pair. `"first"` is a useful starting point for a long series. |
| `sensitivity_policy` | `"none"`, `"first"` or `"all"`: run the configured fluorescence/watershed parameter variants for none, the first pair or every pair. These extra segmentations cost time. |
| `export_level` | `"minimal"`, `"full"` or `"summary_only"`; see below. |
| `pair_error_policy` | `"stop"` stops processing after a pair error; `"continue"` records the failure and tries the later planned pairs. This is separate from the policy for missing input files. |
| `plot_every` | Refresh the accumulated saturation plot after this many pairs; the final plot is also saved. The CSV is saved after each completed/failed pair, independently of plot frequency. |
| `preflight_shapes` | `True` checks planned series dimensions before the expensive analysis. Keep it enabled unless there is a specific reason to defer checks. |
| `retain_last_pair` | Keep the last analysed pair's arrays in the returned Spyder result for inspection. Earlier pair arrays are not accumulated in memory. |

Export levels control scientific image files; reporting is a separate choice:

- **`minimal`** saves nominal phase labels, phase labels with review pixels, and review flag bits for each pair, plus metadata and saturation. This is the normal batch choice.
- **`full`** additionally saves grain/pore masks, the signed subtraction, fluorescence/shiny-fringe masks, region labels, liquid/air masks and other diagnostic images. Configured sensitivity maps are included when those variants are run.
- **`summary_only`** omits per-pair scientific TIFFs while retaining measurements and metadata. Select reports separately with `report_policy`; requesting a report can still create figures. Use this mode when a compact measurement run is sufficient, because omitted phase TIFFs cannot be inspected later without rerunning.

The nominal scientific analysis retains the configured fluorescence priority and interface reasoning at every export level. Report/export choices do not turn dark B into air. The sensitivity policy controls additional parameter trials; choosing `"none"` does **not** establish that the nominal map is certain. In particular, disagreement flags from variants that were not run cannot appear. Review fractions from runs with different sensitivity policies therefore need care when comparing them.

For a new acquisition, start with a short range and `report_policy="first"`, preferably with `sensitivity_policy="first"`, to inspect grain alignment, fluorescence cutoff and shiny-edge handling. After that check, extend the range. No automatic setting replaces checking a few representative bright, diluted and air-dominated regions.

## 5. Read the saved liquid–air saturation correctly

For every successfully analysed pair, the nominal pore pixels have exactly two fluid labels:

```text
Npore   = number of pixels inside the bed and outside the fixed grains
Nliquid = pore pixels labeled liquid (2)
Nair    = pore pixels labeled air (3)

Sw = Nliquid / Npore
Sa = Nair / Npore
Sw + Sa = 1
```

Grains and pixels outside the bed are excluded from the denominator. This is a **projected two-dimensional pore-area fraction**, not a measured three-dimensional pore-volume saturation. It does not require pixel-size calibration to express the fraction. Different depths, unresolved films and optical fringe widths can affect its relation to physical volume.

Fractions in CSV/JSON use **0–1**; plots show percent. For example `liquid_saturation=0.63` means 63% of analysed pore area has a nominal liquid label. These nominal totals include pixels flagged for review, so read the review information alongside them.

Important `saturation.csv` columns include:

| Column | Meaning |
|---|---|
| `pair_number`, frame indices and source paths | Identify the actual requested pair and its illumination assignment. |
| `time_seconds` | Optional elapsed pair time; blank when no interval was supplied. |
| `status`, `error` | Distinguish completed measurements from failures or other recorded outcomes. Missing measurements are not zero saturation. |
| `pore_pixel_count`, `liquid_pixel_count`, `air_pixel_count` | Counts behind the nominal fractions. |
| `liquid_saturation`, `air_saturation` | Nominal Sw and Sa, each divided by all pore pixels. |
| `reviewed_pore_pixels`, `review_fraction` | Count/fraction of pore pixels carrying any review flag. |
| `unreviewed_pore_pixels`, `unreviewed_liquid_pixels`, `unreviewed_air_pixels` | Counts after excluding all flagged pixels. |
| `liquid_saturation_unreviewed` | Unflagged liquid / all unflagged pore pixels. The denominator differs from nominal Sw; blank if no unflagged pore pixels remain. |
| `liquid_saturation_lower_possible`, `liquid_saturation_upper_possible` | Review sensitivity range, described below. |
| `blue_clipped_fraction`, `RB_clipped_fraction` | Pore fraction clipped at the camera's 255-count ceiling in each raw channel; separate from fluid occupation. |

The review range holds all unflagged labels fixed and reassigns every flagged pore pixel:

```text
Sw_lower = unflagged liquid pixels / all pore pixels
Sw_upper = (unflagged liquid pixels + flagged pore pixels) / all pore pixels
```

It is a **review sensitivity range, not a confidence interval** and not a calibrated accuracy bound. Unflagged pixels have not been independently validated and may also be wrong. Air has corresponding columns. The plot's shaded liquid range uses these assumptions; it does not represent experimental repeatability or uncertainty in the grain mask.

“Camera saturation” means something else: clipping at 255 counts loses intensity information. Its B and R&B columns report that separately within pore space. A fluorescent liquid pixel can be clipped and still labeled liquid; clipping fraction is not air fraction or water saturation.

## 6. Find the results and resume an interrupted run

A batch run uses this structure:

```text
run_directory/
    reference/                 fixed grain / bed / pore geometry and provenance
    pairs/                     one directory per analysed pair
    pair_plan.csv              requested pairing, illumination and skipped slots
    saturation.csv             accumulated measurements and pair outcomes
    saturation_vs_pair.png     nominal phase fractions and review information
    batch_summary.json         run-level progress, settings/provenance and outcomes
```

Each pair directory also retains settings, source information and summary metadata according to the batch implementation; its TIFFs and figures depend on the export/report choices above. The per-pair nominal phase codes are **0 outside, 1 PDMS, 2 liquid, 3 air**. In the review phase image, **4 means review/unresolved**. `08_review_flags.tiff` keeps the individual reasons as a bit mask, described by `review_flag_codes.json`.

The accumulated CSV is saved after each pair by writing a temporary file in the same directory and replacing the old CSV only after the new file has been flushed. A failure during that save preserves the previous complete CSV. Avoid keeping the CSV open in software that locks it while analysis is writing; read a copy if needed. A sudden interruption can leave the currently running pair incomplete, while earlier persisted rows remain available.

To resume:

1. Keep the same scientific settings, reference masks/images and exact source/pair plan.
2. Set `SETTINGS["batch"]["resume_from"]` to the existing batch run directory.
3. Run `RUN_ANALYSIS.py` again. Completed pairs are retained; incomplete or failed work can be retried as supported by the recorded state.

Resume verifies the scientific configuration, reference SHA256 hashes, algorithm-module SHA256 hashes, and each source frame path, file size and modification time. Per-pair source SHA256 hashes are also stored when hash_pair_sources is enabled. This prevents silently mixing different experiments or code versions. To change the requested range, illumination parity, grain geometry, input sources or scientific thresholds, start a **new run** with `resume_from=None`.

Runtime choices such as report/export policies can change for remaining work. Already completed pairs keep their earlier exports and recorded settings. For example, resuming a `summary_only` run with `export_level="full"` does not retrospectively create full masks for pairs already marked completed. Start a fresh run if you need uniform exports or parameter trials across every pair.

## 7. Checks before interpreting a long series

The packaged workflow has been checked with the supplied **one real experimental pair** and synthetic tests for pairing, masks, saturation arithmetic, persistence and related batch behavior. This is not validation on an acquired long experimental series.

For future recordings, check that the illumination assignment really stays alternating, that B/R&B exposure and gain remain appropriate for the thresholds, and that the grain reference stays aligned throughout the sequence. The supplied camera scale is 0–255; different acquisition scales are rejected rather than silently normalized. The code does not infer parity or automatically correct scene motion, changing geometry or illumination drift. A fixed grain mask can remain numerically consistent while becoming physically misaligned.

Inspect representative early, middle and late pairs, including diluted-water regions. Weak interfaces, grain-edge glare, clipping and bright-fringe topology can still make phase labels ambiguous. Keep the nominal and review maps together when using the saturation curves for scientific interpretation.

The first report/sensitivity slot means the first complete pair in the validated plan. Missing-pair skips do not take that slot; a failed first available pair does not automatically move the slot to another pair.

# Run the batch analysis in Spyder

1. Extract the entire folder; keep phase_analysis/ next to RUN_ANALYSIS.py.
2. Open settings.py. The default demo already points to the bundled 6544/6545 pair and separate initial image 0001.
3. For your own experiment, set input_folder, mask_input, and the series start_index, end_index and blue_parity. End index is inclusive. See BATCH_GUIDE.md for examples.
4. Open RUN_ANALYSIS.py in Spyder and press F5. Choose the current console if asked. Paths are resolved from this folder, so Spyder's working directory can be elsewhere.
5. Find the timestamped run in results/. Open saturation.csv and saturation_vs_pair.png. The first pair's folder contains GUIDE.html and ten explanatory figures by default.

Open the root GUIDE.html to inspect the supplied tested example immediately.

## Dependencies

Use Python 3.10 or newer. Install NumPy, SciPy, scikit-image, tifffile, Matplotlib and Pillow in the same environment used by the Spyder kernel. If needed, activate that environment in Anaconda Prompt and run:

```text
python -m pip install -r "C:\path\to\Spyder_Phase_Analysis_Batch\requirements.txt"
```

INSTALL_REQUIREMENTS.bat uses the python on PATH. Check Spyder's interpreter with:

```python
import sys
print(sys.executable)
```

Restart the kernel after installing packages. The folder needs no Codex runtime, chat cache, work/ directory or internet connection during analysis. Source images are read only.

## Choose the inputs

- mask_input.mode="binarize": supply a separate raw image in image_path; local Phansalkar binarization is done once.
- mask_input.mode="load_binary": supply grain_mask_path and white_is_grain; no re-thresholding, cleanup or resizing is applied. Optional display_reference_path is only for figures.
- input_folder contains alternating TIFFs. Set absolute B parity manually: "even" means 0000,0002,... are B; "odd" means 0001,0003,... are B. The other parity is R&B.
- Pairs are anchored at start_index: (start,start+1), (start+2,start+3), etc. Missing frames never shift later pairs.
- 0..1999 contains 1000 complete pairs. 0..2000 contains a lone last frame; choose a complete range or explicitly set unpaired_policy="drop_last".

The supplied bed polygon and expected_shape_hw fit this particular field of view. Review them for a new acquisition. All masks must already share the series pixel grid. The supplied analysis uses 0–255 camera counts and blue_scale=1.1; different gain/exposure or bit depth needs parameter review. Dark B alone does not imply air.

## Outputs and speed

The reference/ directory stores common grain, pore and bed masks. Each pairs/pair_.../ directory stores its own settings, source hashes, evidence, summary and saturation. Default export_level="minimal" saves three TIFFs: nominal phase, review phase and review flags. Choose "full" for all diagnostics or "summary_only" to omit TIFF arrays. Reports are controlled independently by report_policy="first", "all" or "none".

saturation.csv is replaced safely after every completed or failed pair. Saturation is the liquid or air **2D pixel area divided by the fixed analysed pore area**; it is not a measured 3D volume fraction. CSV fractions use 0–1, plots use percent. Review fraction and review range accompany nominal estimates. The range holds unflagged labels fixed and is not a confidence interval. Camera clipping has separate columns.

For speed, batch sensitivity_policy defaults to "none". "first" or "all" runs the configured extra fluorescence/watershed variants. A disabled sensitivity test supplies no evidence of stability; resulting review fractions can differ between policies.

To resume, set batch.resume_from to the existing timestamped run directory. Completed pairs are skipped; failed/interrupted pairs are retried. References, scientific settings, code, exact range and source-file plan must match. Runtime export/report choices can change for remaining work, but earlier pairs retain their original outputs. Start a fresh run to extend the range or change thresholds/parity.

## Spyder variables

RESULT['run_directory'] gives the output folder; RESULT['saturation_rows'] contains the scalar table. RESULT['last_pair'] contains only the last processed pair's arrays, for example RESULT['last_pair']['phase']; it can be None if retention is disabled or every pair was already completed on resume. Earlier arrays are read from their saved files. run_mode="single" remains available for an explicitly named pair and returns its arrays directly.

Nominal labels: 0 outside, 1 PDMS, 2 liquid (dark or fluorescent), 3 air. Review labels additionally use 4 for flagged pore pixels. Keep nominal and review results together. The inference remains exploratory; see the illustrated method and TESTED_RUN.md for its scope.

# SDM pipeline changes

For the Shiny app maintainer: what changed in the inputs, outputs and
behaviour of the `sdm_*` scripts, and what the app needs to handle.

- **Compared against:** scaleosdm_gee commit `51d2251`, the version these
  scripts were copied from before the first sync (PR #1).
- **Current as of:** scaleosdm_gee `0c86769` (21 Oct 2026). The newest
  changes are in the first section below.

The full history and the measurements behind each change are in `CLAUDE.md`
of the scaleosdm_gee repository (github.com/forest-rev/scaleosdm_gee,
private).

## Update, 21 Oct 2026: fitting and prediction fail loudly

These changes respond to the `m4` run of 20 Oct. Wall-to-wall warned that
the model had no `predictor_source`, `fill_values` or `gap_years`, and its
feature list included `yyyy`.

### Refitting under an existing model ID now stops at the start

- **What went wrong.** GEE export tasks can't overwrite an existing asset:
  the task fails and the old asset stays. Fitting started its three exports
  (classifier, tree strings, selected features) without checking them. So
  refitting `m4` left the previous fit's assets in GEE, and wall-to-wall
  predicted with the old model, while the new `.pkl` sat on disk.
- **Now:** before fitting, the script checks for this model's outputs: the
  three GEE assets and the local `_rf_final.pkl`, `_selected_features.csv`,
  `_rf_sizing.csv` and `_rf_summary.csv`. If any exist, it stops straight
  away with a list of them. The cached `_training.csv` doesn't count,
  because fitting reuses it.
- **To refit deliberately:** pass `--overwrite` (`overwrite = TRUE` in
  `sdm_model_fitting_py()`). It deletes the old GEE assets just before
  exporting the new ones.
- **Fitting waits for its exports** and prints `OK` or `FAIL` for each. If
  any fails, the script exits with an error after saving the local outputs,
  so the GUI can tell the fit isn't usable for prediction yet.

For the GUI: a user who re-runs a fit with the same model ID now gets the
"already been fit" message. You may want an "overwrite" checkbox on the
fitting panels, passed as `overwrite = TRUE`.

### Prediction checks the model's predictors before starting

`m4` listed `yyyy` as a predictor. That's a column of the occurrence data
which the extraction of the time copied into the training tables. The
predictor stack can't build it, so every tile would have failed.
Wall-to-wall now checks this once, before the cost estimate (so
`estimate_only` shows it too). It stops with a message to re-extract with
the current scripts and refit.

### Your suggestions

- **Keep the gap-fill layers:** `--keep_intermediate` on
  `sdm_modis_export.py` (`keep_intermediate = TRUE` in
  `sdm_modis_export_py()`). It keeps the annual, merged and gap-filled
  assets in `{BASE_FOLDER}/temp` instead of deleting each once the next
  step has used it. It's off by default because of storage: on the 2° test
  tile these took 173 MB against 7 MB for the metrics, which would be about
  115 GB or more for the puma range. Re-running a stage still works with
  them kept (tested).
- **Don't train if the outputs exist:** covered by the check above, for both
  the GEE assets and the local `.pkl`/`.csv` files.
- **`tile_degrees = 0`** in export and prediction (and
  `sdm_cost_estimate.py`) means "choose from the scale": 2° at 250 m, 0.25°
  at 30 m, 0.1° at 10 m. In prediction, `NULL` works the same way. When
  prediction uses exported metrics, export and prediction must still end up
  with the same tile size; passing 0 to both does that.

### Also fixed

A relative `--working_dir` made fitting write its outputs one folder too
deep (`out/out/...`), because it changes into that folder after building
the paths. The GUI passes absolute paths, so it wasn't affected.

## New default run sequence

The MODIS export step is now optional. Wall-to-wall computes the imagery
predictors itself for every source; in testing that was cheaper and gave the
same predictions. Most runs only need three steps:

```
python sdm_extraction.py    --species S --model_id M --run_mode full --working_dir DIR [--column_train presence] [--predictor_source modis] [--scale 250]
python sdm_model_fitting.py --species S --model_id M --working_dir DIR [--predictor_source modis] [--scale 250]
python sdm_wall_to_wall.py  --species S --model_id M --target_year Y --run_mode full
```

Run the export (`sdm_modis_export.py`) only when two or more models will map
the same species and year, and then pass `--use_exported_metrics` to
wall-to-wall. For puma range-wide with MODIS, the estimates are:

- **Without export:** ~46 EECU-hours, ~47 h wall-clock.
- **With export:** ~66 EECU-hours, ~3.5 days wall-clock.

The GUI currently still runs the export steps and passes
`use_exported_metrics = TRUE`.

## Changes the app must handle

### Two scripts renamed

Both now handle every imagery source, so "modis" was dropped from their
names. The export script keeps its name because it only exports MODIS.

| Before | Now | Where cola refers to it |
|---|---|---|
| `sdm_modis_extraction.py` | `sdm_extraction.py` | default `pyscript` of `sdm_modis_extract_py()`, `R/cola_fun.R` |
| `sdm_modis_wall_to_wall.py` | `sdm_wall_to_wall.py` | default `pyscript` of `sdm_modis_prediction_py()`, `R/cola_fun.R`; path in `inst/app/app.R` |
| `sdm_modis_export.py` | unchanged | — |

The scripts now import three helper modules, `predictors.py`, `grid.py` and
`cost.py`, which must sit in the same folder as the scripts.

### Output lines containing `%` broke `catandcapt()` (fixed in PR #1)

`catandcapt()` echoed each Python output line with `sprintf(line)`. R's
`sprintf()` errors on a bare `%` ("too few arguments" or "unrecognised format
specification"; checked in R 4.5.3), so console output stopped at the first
such line. The scripts print percentages in several places, including every
cost estimate, so they keep their `%` signs. The fix is on the R side:

```r
msg <- line      # was: msg <- sprintf( line )
```

### Removed command-line arguments

Passing any of these to a script stops it with "unrecognized arguments". The
R wrappers still accept them, so existing calls keep working: they report
them as unused and don't pass them on.

| Script | Removed | Instead |
|---|---|---|
| `sdm_extraction.py` | `--min_year` `--max_year` `--target_scale` `--point_buffer` | Use `--scale` (the wrapper maps `target_scale` to it). Years are checked against the source's data coverage. |
| `sdm_modis_export.py` | `--min_year` `--max_year` `--crs` | Nothing; fixed internally. |
| `sdm_wall_to_wall.py` | `--min_year` `--max_year` `--crs` `--scale` | Source and scale are read from the trained model. |

`--working_dir` and `--column_train` on extraction still work as the app
passes them (see below).

### GEE folder and log names include source and scale

For MODIS at the default 250 m the tag is `modis_250m`. The new scripts don't
find existing folders under the old names. Those folders predate the October
2026 consistency fixes, so regenerate them rather than renaming them.

| What | Before | Now |
|---|---|---|
| Training data (extraction output, fitting input) | `{species}_SDM_modis_exports_{model_id}` | `{species}_SDM_{source}_{scale}m_exports_{model_id}` |
| Exported MODIS metrics (optional step) | `{species}_SDM_modis_{year}/metrics` | `{species}_SDM_modis_250m_{year}/metrics` |
| Extraction resume log | `{working_dir}/completed_batches_{species}_{model_id}.txt` | `{working_dir}/completed_batches_{species}_{source}_{scale}m_{model_id}.txt` |
| Extraction task names | `SDM_MODIS_{species}_{year}_bNN` | `SDM_{SOURCE}_{species}_{year}_bNN` (unchanged for MODIS) |

These are unchanged:

- the GEE model assets `{model_id}_RF_classifier`, `_RF_classifier_strings`
  and `_selected_features`
- the prediction folder `{model_id}_prediction_{year}` and the prediction
  tile names
- the local `.pkl` and `.csv` file names from fitting

### Wall-to-wall no longer reads exported metrics by default

If the app runs the export and wants wall-to-wall to use it, pass
`--use_exported_metrics` (`use_exported_metrics = TRUE` in the wrapper).
`--metrics_folder` now works only together with that flag. The flag is
refused for HLS and embedding models, which have no export step.

### Extraction stops early on bad occurrence data

In each of these cases extraction exits with code 1 before any task starts,
with a message the app should show to the user:

- **Years outside the source's coverage.** Previously years before
  `--min_year` were silently dropped and years after `--max_year` were
  clamped. Now the out-of-range years and their point counts are listed.
- **Missing response column.** If the `--column_train` column isn't on the
  occurrence asset, the message lists the properties the asset does have.
- **Years that aren't whole numbers, or are missing.** Years stored as
  decimals such as `2019.0` are fine and are converted automatically.

### `--run_mode single` picks points differently

It used to select points by the occurrence asset's `batch` property. It now
takes the n-th group of `--batch_size` points in `--batch_year`, which is the
same group full mode builds for that batch number.

### Wall-to-wall test mode runs one tile

It used to run `test_costa_rica_{d}deg` and `test_costa_rica_south_{d}deg`.
Now it runs only `test_costa_rica_{d}deg`; tiles smaller than 1° use
`test_costa_rica_inland_{d}deg` instead.

### Console messages changed

This only matters if the app parses script output.

- Emoji markers are plain words: ✅ → `OK`, ❌ → `FAIL`, ⚠️ → `WARNING`.
- Extraction prints `OK <task> -- N rows, X EECU-s`. When a batch wrote
  fewer rows than expected, it prints `SHORT <task> -- n/N`.
- Every launching script prints a cost-estimate table before it starts any
  task.
- Wall-to-wall prints a `Gap-fill:` line. It warns when a model predates
  the `gap_years` or `fill_values` records.

## Ivan's edits, now in the scripts

The edits made to cola's copies are now part of the scripts themselves, so
updating the copies keeps them.

- **`--working_dir` on extraction.** The resume log goes in this folder,
  which is created if missing. It's optional and defaults to the current
  directory. Fitting also creates its `--working_dir` if it doesn't exist.
- **`--column_train` on extraction.** This is the response column; it can
  also be given as `--target_col`, matching fitting, and defaults to
  `presence`. Only `year` and this column are copied from the occurrence
  asset into the training data. Fitting treats every other numeric column
  as a candidate predictor, so this keeps user properties such as animal
  IDs or counts out of the model.
- **Whole-number years.** The `int()` casts are now a check on every
  occurrence year, described under "Extraction stops early" above.
- **`%` signs in printed output.** These were not removed from the scripts;
  `catandcapt()` was fixed instead (see above).

## New options

All are optional. Apart from the in-memory default above, the defaults
reproduce the old MODIS 250 m behaviour.

- **`--predictor_source`** (extraction and fitting; `predictor_source` in
  the wrappers): one imagery source per model. Static predictors are always
  included. The options are:
  - `modis` (default): 250 m, 2001+
  - `hls_l30`: 30 m, 2013+
  - `hls_s30`: 30 m, 2016+
  - `embed` (Google Satellite Embedding): 10 m, 2017–2025
- **`--scale`** (extraction, fitting and export): defaults to the source's
  native resolution and can only be coarser. Source and scale must match
  between extraction and fitting; wall-to-wall reads both from the model.
- **`--estimate_only`** (extraction, export and wall-to-wall;
  `estimate_only` in the wrappers): prints the cost and exits without
  starting tasks.
- **`sdm_cost_estimate.py`**: estimates the whole pipeline before anything
  runs: EECU-hours, task count, wall-clock and share of monthly quota. It
  also takes `--use_exported_metrics`.
- **`--tile_degrees`** (export and wall-to-wall): accepts fractions. By
  default it's chosen from the scale (2° at 250 m, as before).

## Output contents

| Output | Change |
|---|---|
| Training data (GEE tables, local CSV) | Only `year` and the `--column_train` column are carried over from the occurrence asset. New columns `predictor_source`, `scale` and `gap_years`. HLS and embeddings have different predictor columns: monthly HLS indices, bands A00–A63, and suffixes such as `_30m` and `_10m`. MODIS 250 m column names are unchanged. |
| `{model_id}_selected_features.csv` | New `fill_value` column. |
| `{model_id}_rf_summary.csv` | New rows `predictor_source`, `scale` and `gap_years`. The SEDI sign is corrected, so older summaries report SEDI with the wrong sign. |
| `_selected_features` GEE asset | New properties `predictor_source`, `scale`, `gap_years` and `fill_values` (JSON). |
| Prediction tiles | Unchanged format: one float band, `suitability`, with ocean and water masked. |

## Values changed

Predictor values differ from older runs, because:

- extraction gap-fill bugs were fixed
- the Aqua product is now 250 m
- training samples on the prediction grid
- tile edges are handled better

Regenerate models and predictions made before 6 October 2026 rather than
mixing them with new ones.

## For the GUI, later (listed with PR #1)

- **Unused inputs:** the year ranges (`*_yy`), CRS (`*_crs`), point buffer
  (`*_pbuffer`) and prediction scale (`in_eepre_scale`) no longer do
  anything.
- **Scale and source:** the scale inputs allow 90 m, but MODIS can't go
  below 250 m (the script stops with a message). Using HLS or embeddings
  would need a `predictor_source` selector.
- **Export steps:** these are optional now. Prediction can compute the
  metrics itself (`use_exported_metrics = FALSE`), which costs less for a
  single model: about 41 vs 62 EECU-h for puma range-wide. If the GUI drops
  the export steps from its default flow, it should pass `FALSE`.
- **Standalone fitting panel:** it has no scale input, so it assumes the
  native 250 m. That matches the extraction panel's default.

# 1 · Landscape metric extraction

This part turns raw forest-inventory and human-footprint layers into per-site
**habitat amount** and **fragmentation** metrics. It is the upstream data
processing for the occupancy analysis in `../2_occupancy_analysis/`.

**No data are provided here — this is code plus documentation.** The raw layers
we used are licensed and not redistributable. The two scripts are documented so
you can run the same workflow on *your own* vegetation/disturbance data (for a
different study area, or different species with different habitat definitions).
If you only want to reproduce the published occupancy models, skip this part —
the metric tables it produces are provided on Dryad (see the top-level README).

The workflow has two steps that run in different software:

```
1_gis_habitat_processing.ipynb   ArcGIS Pro / arcpy   ->  habitat patches (shapefiles) + point tables
2_landscape_metrics.R            R / landscapemetrics ->  hyp{1,2,3}_landscape_metrics_2m_res.csv
```

The three output CSVs are the hand-off to the occupancy analysis: copy them into
`../2_occupancy_analysis/data/`.

---

## Step 1 — GIS habitat processing (`1_gis_habitat_processing.ipynb`)

**Software:** ArcGIS Pro with the `arcpy` Python environment (Python 3.11). Runs
inside an open ArcGIS Pro project (`arcpy.mp.ArcGISProject("CURRENT")`).

**Inputs you supply** (as layers in the ArcGIS Pro project):
- A combined vegetation + disturbance layer per survey year. Ours was named
  `veg_{year}_wfp_hf_ght_final` and carried these attribute fields, which the
  code queries — rename yours to match, or edit the queries in the notebook:
  - `Combined_v7` — vegetation class (we keep `Decid`, `Mixedwood`, `Spruce`)
  - `Origin_Year` — stand origin year (`9999` = natural/old forest)
  - `Origin_Source_NatDist` — natural-disturbance source (used to find fire)
  - `YEAR` — event year (for fire and harvest recency)
  - `FEATURE_TY` — disturbance/footprint feature type (drives the definitions)
- A survey-points layer with a `survey_year` field. Ours was
  `landscape_summary_camaru_pts_2010_2023`.

**What it does**, per year × fragmentation definition (hyp1/2/3):
1. Selects suitable habitat (`Combined_v7` in your habitat classes AND
   `Origin_Year = 9999` OR older than your age threshold — we used ≥ 80 yr).
2. Erases recent fire and harvest (within your age threshold).
3. Erases disturbances for the definition — the `FEATURE_TY` lists in
   `get_disturbance_query()` are increasingly inclusive:
   - **hyp1 (Def A):** polygonal disturbances only
   - **hyp2 (Def B):** + wide linear features (roads, pipelines, transmission)
   - **hyp3 (Def C):** + narrow linear features (seismic lines, trails)
4. Dissolves the remainder into single-part patches and computes each point's
   distance to the nearest patch.

**Outputs** (per definition h and year — these are the inputs to Step 2):
- `dissolved_habitat_hyp{h}_{year}.shp` — suitable-habitat patches
- `patch_metrics_hyp{h}_{year}.csv` — survey points + metadata + distance

**To adapt for another species / system:** change the habitat classes and age
threshold in the suitable-habitat query, and (if needed) the `FEATURE_TY` lists
that define each fragmentation level. Everything downstream is generic.

---

## Step 2 — Landscape metrics (`2_landscape_metrics.R`)

**Software:** R. Packages: `sf`, `landscapemetrics`, `raster`, `fasterize`,
`tidyverse`.

**Inputs:** the Step 1 outputs, in one folder (`input_dir` at the bottom of the
script):
- `dissolved_habitat_hyp{h}_{year}.shp`
- `patch_metrics_hyp{h}_{year}.csv` (needs a coordinate pair the script can find:
  `Easting`/`Northing`, `POINT_X`/`POINT_Y`, or `X`/`Y`)

**What it does:** for each point, buffers at 150 / 500 / 1000 m, rasterises the
clipped habitat at **2 m** (fine enough to preserve 2 m seismic lines), and
computes the class-level metrics for the habitat class:
`prop_habitat`, `edge_density`, `aggregation_index`, `cohesion`, `n_patches`,
`mean_patch_area`, `contagion`. Metrics undefined without habitat are returned
as `NA`.

**Output** (one per definition — the hand-off to the occupancy analysis):
- `outputs/landscape_metrics/2m_resolution/hyp{h}_landscape_metrics_2m_res.csv`

Copy the three `hyp*_landscape_metrics_2m_res.csv` files into
`../2_occupancy_analysis/data/`.

**Settings** are at the top (`RESOLUTION`, `BUFFER_SIZES`) and the bottom
(`input_dir`, `output_dir`). The main loop runs all three definitions
(`for (hyp in 1:3)`); uncomment the `test_year` call to trial one year first.

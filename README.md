# Beyond habitat amount: A species-specific approach reveals narrow linear features fragment old-forest songbird habitat

Dataset DOI: [10.5061/dryad.05qfttfhk](https://doi.org/10.5061/dryad.05qfttfhk)

## Description of the data

Code and data used to assess the response of Black-throated Green Warblers (BTNW; *Setophaga virens*) to anthropogenic fragmentation and habitat loss in the boreal and foothills regions of Alberta, Canada. BTNW detection/non-detection histories were compiled from autonomous recording units (ARUs) and human point counts at 508 survey locations, each visited up to three times in each of two years (3,048 visits total). 

For every survey location we calculated the amount of suitable BTNW habitat (defined as deciduous, white spruce, or mixedwood forest aged 80-140 years) and the degree of fragmentation (via edge density and other landscape metrics) within 150, 500, 1000 m buffers (7, 78, 314 ha), using three increasingly inclusive defintions of antrhopogenic fragmentation (A = large polygonal disturbnaces only; B = + wide linear features; C = + narrow linear features such as seismic lines). These data were used to fit multi-season spatial occupancy models to assess the impacts of habitat amount and fragmentation while accounting for imperfect detection and spatial autocorrelation. 

Missing values are recorded as **NA** throughout, including configuration metrics (aggregation_index, cohesion, contagion) that are undefined where a buffer contains no habitat. Coordinates are NAD83 / Alberta 10-TM (Forest), EPSG:3400 (metres).

### Files and variables

**2yearvisitmatrix.csv** — BTNW detection/non-detection history; one row per survey location.

* gisid: unique survey-location identifier
* 1_1, 1_2, 1_3: year 1, visits 1–3 — BTNW detected (1) or not detected (0); NA = no survey
* 2_1, 2_2, 2_3: year 2, visits 1–3 — BTNW detected (1) or not detected (0); NA = no survey

**visitsfordetcovs.csv** — per-visit detection covariates and survey metadata; one row per visit.

* surveyid: unique survey (visit-level) identifier
* gisid: survey-location identifier (links to the other files)
* BTNW: BTNW count at the visit
* BTNW_binary: BTNW detected (1) or not detected (0)
* offset: QPAD statistical offset for detectability and effective area (log scale)
* project_id, project, source, organization: data-provenance fields
* location: location name or label
* sensor: survey method — ARU or HPC (human point count)
* buffer: point-count / ARU detection-radius setting (m)
* latitude, longitude: location in geographic coordinates, NAD83 (decimal degrees)
* Easting, Northing: location in projected coordinates, EPSG:3400 (m)
* year: calendar year of survey
* year_id: survey-year index (1 or 2)
* visit_num: visit number within the year (1–3)
* visit_id: unique visit identifier
* date_time: date and time of the survey (YYYY-MM-DD HH:MM)
* task_method: survey protocol / task method
* duration: survey / recording duration category (min)
* distance: detection distance, where recorded (m)
* topsecret, useNorth, useSouth: internal data-handling flags

**hyp1_landscape_metrics_2m_res.csv**, **hyp2_landscape_metrics_2m_res.csv**, **hyp3_landscape_metrics_2m_res.csv** — landscape metrics per survey point and buffer size. hyp1 = fragmentation Definition A, hyp2 = B, hyp3 = C. The three files share the same columns. Metrics were computed on habitat rasterised at 2 m resolution.

* year: calendar year the metrics correspond to
* surveyid: survey-location identifier
* survey_year: survey year for that location
* veg_type: dominant suitable-habitat vegetation class (Decid, Mixedwood, or Spruce)
* origin_year: stand origin year (9999 = natural / old forest)
* buffer_size: buffer radius around the point (m; 150, 500, or 1000)
* resolution: raster resolution used for the metrics (m; 2)
* dist_to_habitat: distance to the nearest habitat patch (m; 0 = point is in habitat)
* prop_habitat: proportion of the buffer that is suitable habitat (0–1)
* edge_density: habitat edge density (m/ha)
* aggregation_index: aggregation index of habitat (%, 0–100; NA if no habitat)
* cohesion: patch cohesion index (0–100; NA if no habitat)
* n_patches: number of habitat patches (count)
* mean_patch_area: mean habitat patch area (ha)
* contagion: contagion index (%, 0–100; NA if no habitat)

## Code/software

Analyses were run in R version 4.4.2 (open source; [https://www.r-project.org](https://www.r-project.org)). Each script loads the packages it needs:

* data formatting (01): tidyr, dplyr
* occupancy models and diagnostics (02, 03, 05, 06, 07, 08): spOccupancy, tidyverse, spdep, pROC
* predictor coupling (04): mgcv, tidyverse
* figures (09): ggplot2, dplyr, gridExtra
* upstream landscape metrics (2_landscape_metrics.R): sf, landscapemetrics, raster, fasterize, tidyverse
* upstream GIS step (1_gis_habitat_processing.ipynb): ArcGIS Pro with the arcpy Python environment (Python 3.11); proprietary, and needed only to regenerate the habitat layers from licensed data.

**Workflow:** the three hyp landscape-metric files and the two detection files are the starting point for the occupancy analysis. Run the analysis scripts in order 01 through 09; run 03 before the diagnostic and validation scripts (04–08), which read the fitted models. The upstream habitat processing code shows how the hyp landscape-metric files were produced and can be adapted to other study areas or species by changing the habitat definition and disturbance queries.

### 1_landscape_metric_extraction

**Upstream habitat processing code:** in the repository but not required to reproduce the models, because its outputs (the three hyp landscape-metric files) are provided above.

* 1_gis_habitat_processing.ipynb: ArcGIS Pro / arcpy. Selects suitable habitat, removes fire/harvest and definition-specific disturbances, and exports dissolved habitat patches and survey-point tables. Requires licensed vegetation and human-footprint layers, which are not shared
* 2_landscape_metrics.R: R / landscapemetrics. Buffers each survey point at 150, 500, and 1000 m, rasterises the clipped habitat at 2 m, and computes the metrics in the hyp landscape-metric files

### 2_occupancy_analysis

**Occupancy analysis code (R):**. Run in order; set the working directory to the folder holding the scripts and data. Script 01 regenerates the nine def{A,B,C}{150,500,1000}_glmmdata.csv modelling tables from the files above, so those are not deposited separately.

* 01_format_data.R: joins the detection history and detection covariates to the landscape metrics for each fragmentation definition and spatial scale, standardises the covariates, and writes the nine modelling tables.
* 02_spatial_vs_nonspatial.R: fits a non-spatial occupancy model and tests its residuals for spatial autocorrelation (Moran's I at several distances, plus a variogram) to justify using a spatial model.
* 03_fit_spoccupancy_models.R: fits the four candidate spatial occupancy models (habitat only, edge only, additive, interaction) for each definition and scale and compares them with WAIC (Table S2).
* 04_predictor_coupling.R: measures how strongly each fragmentation metric is coupled to habitat amount, comparing linear (Pearson) with non-linear (GAM) diagnostics (Table S1).
* 05_structure_coefficients.R: compares each predictor's regression slope with its structure coefficient to identify effects driven by collinearity rather than a direct association (Table S3).
* 06_model_diagnostics.R: convergence (trace plots), model fit (Freeman-Tukey posterior predictive check), residual analysis, and Moran's I on residuals for the top-supported model.
* 07_spatial_block_cv.R: 5-fold spatial block cross-validation of the top-supported model, reporting AUC on spatially independent held-out sites.
* 08_alternative_frag_metrics.R: repeats the model comparison using aggregation index and mean patch area in place of edge density (Tables S4 and S5).
* 09_figures.R: produces the main-text figures.

### Access information

Code (all scripts above, including the upstream processing) is available at: [https://github.com/harttaylor/btnw-fragmentation-occupancy.git](https://github.com/harttaylor/btnw-fragmentation-occupancy.git)

Data were derived from:

\- Alberta Biodiversity Monitoring Institute (ABMI). Alberta Vegetation Inventory

and Human Footprint Inventory, 2024. [https://abmi.ca/](https://abmi.ca/)

\- Habitat-suitability age criteria: ABMI and Boreal Avian Modelling Project

(BAM), 2023.

\- BTNW detection data: compiled via WildTrax

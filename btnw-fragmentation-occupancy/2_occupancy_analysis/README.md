# 2 · Occupancy analysis

Fits the multi-season spatial occupancy models and produces every model result,
diagnostic, and figure in the paper. Unlike the upstream processing, **this part
is fully reproducible from the data on Dryad**.

## Data needed (from Dryad → put in `data/`)

| File | Role |
|------|------|
| `2yearvisitmatrix.csv` | detection/non-detection history (508 sites × 2 yr × 3 visits) |
| `visitsfordetcovs.csv` | per-visit detection covariates |
| `hyp1_landscape_metrics_2m_res.csv` | landscape metrics, fragmentation Definition A |
| `hyp2_landscape_metrics_2m_res.csv` | landscape metrics, Definition B |
| `hyp3_landscape_metrics_2m_res.csv` | landscape metrics, Definition C |

Script 01 joins these into the nine `def{A,B,C}{150,500,1000}_glmmdata.csv`
modelling files (written back into `data/`).

## Software

R ≥ 4.4.2. Each script loads only what it needs, so you can run a step months
later without installing the whole stack:

| Script | Packages |
|--------|----------|
| 01_format_data.R | tidyr, dplyr |
| 02_spatial_vs_nonspatial.R | spOccupancy, tidyverse, spdep |
| 03_fit_occupancy_models.R | spOccupancy, tidyverse |
| 04_predictor_coupling.R | mgcv, tidyverse |
| 05_structure_coefficients.R | spOccupancy, tidyverse |
| 06_model_diagnostics.R | spOccupancy, spdep |
| 07_spatial_block_cv.R | spOccupancy, pROC |
| 08_alternative_frag_metrics.R | spOccupancy, tidyverse |
| 09_figures.R | ggplot2, dplyr, gridExtra |

```r
install.packages(c("spOccupancy","tidyverse","spdep","mgcv","pROC","gridExtra"))
```

Set the working directory to this folder (open the `.Rproj`, or `setwd()` to
`2_occupancy_analysis/`). No script sets an absolute path.

## Run order

Run 03 before the diagnostics (04–08): they read its fitted models or the
objects it builds. Within a run, set `DEFINITION` and `SCALE` at the top of a
script and re-run for each combination you need.

1. **01_format_data.R** — build the nine modelling datasets. *(Skip if `data/`
   already has them.)*
2. **02_spatial_vs_nonspatial.R** — Moran's I on non-spatial residuals justifies
   the spatial model. → `outputs/tables/spatial_autocorrelation_results.csv`
3. **03_fit_occupancy_models.R** — fit models 1–4 + WAIC, per definition × scale.
   → `outputs/models/`, `outputs/tables/waic_comparison_*` (**Table S2**)
4. **04_predictor_coupling.R** — habitat vs fragmentation coupling (GAMs).
   (**Table S1**)
5. **05_structure_coefficients.R** — slope vs structure coefficient; flags terms
   that are not interpretable under collinearity. (**Table S3**)
6. **06_model_diagnostics.R** — trace plots, Freeman-Tukey PPC, residuals,
   Moran's I on residuals. → `outputs/figures/diagnostics_*.png`
7. **07_spatial_block_cv.R** — 5-fold spatial block CV of the top model (AUC).
   *(run script 03's data-prep block first)*
8. **08_alternative_frag_metrics.R** — aggregation index & mean patch area.
   (**Tables S4/S5**) *(run script 03's data-prep block first)*
9. **09_figures.R** — main-text figures. → `outputs/figures/`

## Manuscript ↔ script map

| Output | Script |
|--------|--------|
| Main results, WAIC (Fig. 3–5) | 03, 09 |
| Table S1 — predictor coupling | 04 |
| Table S2 — WAIC across definitions/scales | 03 |
| Table S3 — structure coefficients | 05 |
| Tables S4/S5 — alternative metrics | 08 |
| Convergence / fit / residual diagnostics | 06 |
| Spatial-model justification (Moran's I) | 02 |
| Cross-validation (AUC) | 07 |

## Notes

- Missing values are `NA` (including configuration metrics that are undefined
  where a buffer has no habitat).
- CV fold assignment is seeded (`CV_SEED <- 100`). Set a seed before fitting if
  you want bit-reproducible MCMC chains.

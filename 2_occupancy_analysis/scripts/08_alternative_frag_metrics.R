# =============================================================================
# 08  ALTERNATIVE FRAGMENTATION METRICS  (Supplementary Tables S4 / S5)
# =============================================================================
# Run the models using two alternative fragmentation metrics
#     - mean patch area (ha, log-transformed)   
#     - aggregation index                        
#
# IMPORTANT - fair WAIC comparisons:
#   * Mean patch area is defined at every site, so its models use all 508 sites
#     and are compared against the habitat-only model from script 03.
#   * Aggregation index is UNDEFINED where a buffer contains no habitat. Those
#     sites are dropped, and BOTH the alternative models AND a matched
#     habitat-only baseline are refit on that same (smaller) set of sites, so
#     the WAIC values are comparable.
#
# PREREQUISITE: run script 03 first (down through "MCMC SETTINGS") so these
#   objects exist: y, coords, site_covs, det_covs_final, det_formula, priors,
#   tuning, n_batch, batch_length, n_burn, n_thin, n_chains.
#
# OUTPUTS: outputs/models/altmetrics/*.rds
#          outputs/tables/waic_{patcharea,aggregation}_defC_{scale}m.csv
# =============================================================================

library(spOccupancy)
library(tidyverse)

DEFINITION <- "C"
SCALE      <- "150"

MODEL_DIR <- "outputs/models"
ALT_DIR   <- "outputs/models/altmetrics"
TABLE_DIR <- "outputs/tables"
dir.create(ALT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- helper: fit one model (optionally on a subset of sites) -----------------
# Reuses the priors, tuning, and MCMC settings created in script 03.
fit_alt <- function(occ_covs, occ_formula, n_beta,
                    y_use = y, coords_use = coords, det_use = det_covs_final) {
  stPGOcc(
    occ.formula  = occ_formula,
    det.formula  = det_formula,
    data         = list(y = y_use, coords = coords_use,
                        occ.covs = occ_covs, det.covs = det_use),
    inits        = list(beta = rnorm(n_beta, 0, 0.5), alpha = rnorm(4, 0, 0.5),
                        sigma.sq.psi = runif(1, 0.5, 2), phi = runif(1, 0.001, 0.004)),
    priors       = priors, tuning = tuning, cov.model = "exponential",
    NNGP = TRUE, n.neighbors = 5,
    n.batch = n_batch, batch.length = batch_length,
    n.burn  = n_burn,  n.thin = n_thin, n.chains = n_chains,
    n.report = 200, verbose = TRUE
  )
}

# ---- helper: build a WAIC comparison table -----------------------------------
make_waic <- function(models, labels, outfile) {
  w <- sapply(models, function(m) waicOcc(m, by.sp = FALSE))
  tab <- data.frame(Model = labels,
                    elpd = w["elpd", ], pD = w["pD", ], WAIC = w["WAIC", ]) %>%
    mutate(delta_WAIC = WAIC - min(WAIC),
           weight     = exp(-0.5 * delta_WAIC) / sum(exp(-0.5 * delta_WAIC))) %>%
    arrange(WAIC)
  print(tab, digits = 3)
  write.csv(tab, file.path(TABLE_DIR, outfile), row.names = FALSE)
  tab
}

ph <- site_covs$prop_habitat_std   # habitat amount (all sites)

# =============================================================================
# MEAN PATCH AREA  (all 508 sites)
# =============================================================================
# Baseline is the habitat-only model already fit in script 03 (same sites).
model_hab <- readRDS(file.path(MODEL_DIR,
              sprintf("model1_habitat_only_def%s_%sm.rds", DEFINITION, SCALE)))

mpa <- site_covs$mean_patch_area_std

mpa2 <- fit_alt(list(mean_patch_area = mpa), ~ mean_patch_area, 2)
mpa3 <- fit_alt(list(prop_habitat = ph, mean_patch_area = mpa),
                ~ prop_habitat + mean_patch_area, 3)
mpa4 <- fit_alt(list(prop_habitat = ph, mean_patch_area = mpa),
                ~ prop_habitat * mean_patch_area, 4)

summary(mpa2); summary(mpa3); summary(mpa4)

saveRDS(mpa2, file.path(ALT_DIR, sprintf("model2_patch_area_only_def%s_%sm.rds",        DEFINITION, SCALE)))
saveRDS(mpa3, file.path(ALT_DIR, sprintf("model3_habitat_patch_area_add_def%s_%sm.rds", DEFINITION, SCALE)))
saveRDS(mpa4, file.path(ALT_DIR, sprintf("model4_habitat_patch_area_interact_def%s_%sm.rds", DEFINITION, SCALE)))

make_waic(list(model_hab, mpa2, mpa3, mpa4),
          c("M1: Habitat only", "M2: Patch area only",
            "M3: Habitat + Patch area (add)", "M4: Habitat x Patch area (int)"),
          sprintf("waic_patcharea_def%s_%sm.csv", DEFINITION, SCALE))

# =============================================================================
# AGGREGATION INDEX  (habitat-present subset, with a matched baseline)
# =============================================================================
keep <- !is.na(site_covs$aggregation_std)
cat(sprintf("\nAggregation: %d of %d sites used (%d dropped: no habitat in buffer)\n",
            sum(keep), length(keep), sum(!keep)))

# Subset every data object to the same sites.
y_k     <- y[keep, , , drop = FALSE]
coords_k <- coords[keep, , drop = FALSE]
det_k   <- lapply(det_covs_final, function(a) a[keep, , , drop = FALSE])
ph_k    <- site_covs$prop_habitat_std[keep]
agg_k   <- site_covs$aggregation_std[keep]

# Refit the habitat-only baseline on the SAME subset (do not reuse the 508-site
# model above - WAIC must be computed on identical data).
agg1 <- fit_alt(list(prop_habitat = ph_k), ~ prop_habitat, 2,
                y_use = y_k, coords_use = coords_k, det_use = det_k)
agg2 <- fit_alt(list(aggregation = agg_k), ~ aggregation, 2,
                y_use = y_k, coords_use = coords_k, det_use = det_k)
agg3 <- fit_alt(list(prop_habitat = ph_k, aggregation = agg_k),
                ~ prop_habitat + aggregation, 3,
                y_use = y_k, coords_use = coords_k, det_use = det_k)
agg4 <- fit_alt(list(prop_habitat = ph_k, aggregation = agg_k),
                ~ prop_habitat * aggregation, 4,
                y_use = y_k, coords_use = coords_k, det_use = det_k)

summary(agg2); summary(agg3); summary(agg4)

saveRDS(agg1, file.path(ALT_DIR, sprintf("model1_habitat_only_aggSUBSET_def%s_%sm.rds",  DEFINITION, SCALE)))
saveRDS(agg2, file.path(ALT_DIR, sprintf("model2_aggreg_only_def%s_%sm.rds",             DEFINITION, SCALE)))
saveRDS(agg3, file.path(ALT_DIR, sprintf("model3_habitat_aggreg_add_def%s_%sm.rds",      DEFINITION, SCALE)))
saveRDS(agg4, file.path(ALT_DIR, sprintf("model4_habitat_aggreg_interact_def%s_%sm.rds", DEFINITION, SCALE)))

make_waic(list(agg1, agg2, agg3, agg4),
          c("M1: Habitat only (subset)", "M2: Aggregation only",
            "M3: Habitat + Aggregation (add)", "M4: Habitat x Aggregation (int)"),
          sprintf("waic_aggregation_def%s_%sm.csv", DEFINITION, SCALE))

# The structure-coefficient check for these alternative-metric models is done by
# script 05 - point it at these fitted objects to test interpretability.

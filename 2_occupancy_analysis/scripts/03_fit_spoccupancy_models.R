# =============================================================================
# 03  FIT OCCUPANCY MODELS  (main analysis)
# =============================================================================
# Fits the four candidate spatial occupancy models for ONE fragmentation
# definition x spatial scale, and compares them with WAIC:
#     Model 1  habitat amount only
#     Model 2  fragmentation (edge density) only
#     Model 3  additive     (habitat + edge)
#     Model 4  interactive  (habitat * edge)
#
# To reproduce the whole analysis, set DEFINITION and SCALE below and re-run for
# each of the 9 combinations (A/B/C x 150/500/1000)
#
# This script also builds the objects (y, coords, site_covs, det_covs_final,
# det_formula, priors, tuning, MCMC settings) that scripts 04-08 reuse
#
# INPUTS : data/def{D}{scale}_glmmdata.csv, data/2yearvisitmatrix.csv,
#          data/visitsfordetcovs.csv
# OUTPUTS: outputs/models/model{1-4}_..._def{D}_{scale}m.rds
#          outputs/tables/waic_comparison_def{D}_{scale}m.csv
# =============================================================================

library(spOccupancy)
library(tidyverse)

# ---- Choose the definition and scale to fit ---------------------------------
DEFINITION <- "A"     # "A", "B", or "C"
SCALE      <- "150"   # "150" (7 ha), "500" (78 ha), or "1000" (314 ha)

MODEL_DIR <- "outputs/models"
TABLE_DIR <- "outputs/tables"
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD DATA
# =============================================================================
data_file <- sprintf("data/def%s%s_glmmdata.csv", DEFINITION, SCALE)
cat("Loading:", data_file, "\n")

full_data    <- read.csv(data_file)
visit_matrix <- read.csv("data/2yearvisitmatrix.csv")
det_covs     <- read.csv("data/visitsfordetcovs.csv")

# =============================================================================
# DETECTION HISTORY ARRAY  (site x year x visit)
# =============================================================================
# The visit matrix has one row per site and six detection columns: visits 1-3 in
# year 1, then visits 1-3 in year 2. We reshape this into the 3D array spOccupancy
# expects.

y_matrix   <- as.matrix(visit_matrix[, -1])
n_sites    <- nrow(y_matrix)
n_years    <- 2
max_visits <- 3

y <- array(NA, dim = c(n_sites, n_years, max_visits))
for (i in 1:n_sites) {
  y[i, 1, ] <- y_matrix[i, 1:3]   # year 1 visits
  y[i, 2, ] <- y_matrix[i, 4:6]   # year 2 visits
}
cat(sprintf("Detection array: %d sites x %d years x %d visits\n",
            dim(y)[1], dim(y)[2], dim(y)[3]))

# =============================================================================
# SITE-LEVEL (OCCUPANCY) COVARIATES
# =============================================================================
# One row per site. Site covariates do not change between years, so we take the
# year-1 values, order them to match the visit matrix, then standardise. Mean
# patch area is log-transformed first because it is strongly right-skewed.

site_covs <- full_data %>%
  distinct(gisid, year_id, .keep_all = TRUE) %>%
  filter(year_id == 1) %>%
  select(gisid, sensor, Easting, Northing,
         prop_habitat, edge_density, mean_patch_area,
         aggregation_index, n_patches, contagion)

site_covs <- site_covs[match(visit_matrix$gisid, site_covs$gisid), ]

site_covs <- site_covs %>%
  mutate(
    prop_habitat_std    = scale(prop_habitat)[, 1],
    edge_density_std    = scale(edge_density)[, 1],
    mean_patch_area_std = scale(log(mean_patch_area + 0.01))[, 1],
    aggregation_std     = scale(aggregation_index)[, 1],
    contagion_std       = scale(contagion)[, 1]
  )

coords <- as.matrix(site_covs[, c("Easting", "Northing")])

# =============================================================================
# DETECTION COVARIATES  (site x year x visit)
# =============================================================================
# Derive day-of-year and time-of-day from the visit timestamp, code sensor as a
# 0/1 indicator (ARU = 1), standardise the two continuous covariates, then place
# each visit's values into the matching cell of a site x year x visit array.

det_covs_list <- det_covs %>%
  mutate(
    date_time     = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M"),
    day_of_year   = as.numeric(format(date_time, "%j")),
    time_of_day   = as.numeric(format(date_time, "%H")) +
                    as.numeric(format(date_time, "%M")) / 60,
    sensor_binary = ifelse(sensor == "ARU", 1, 0)
  ) %>%
  mutate(
    day_of_year_std = scale(day_of_year)[, 1],
    time_of_day_std = scale(time_of_day)[, 1]
  ) %>%
  select(gisid, year_id, visit_num,
         day_of_year_std, time_of_day_std, sensor_binary) %>%
  arrange(gisid, year_id, visit_num)

doy_array    <- array(NA, dim = c(n_sites, n_years, max_visits))
tod_array    <- array(NA, dim = c(n_sites, n_years, max_visits))
sensor_array <- array(NA, dim = c(n_sites, n_years, max_visits))

for (i in 1:nrow(det_covs_list)) {
  site_idx  <- which(visit_matrix$gisid == det_covs_list$gisid[i])
  year_idx  <- det_covs_list$year_id[i]
  visit_idx <- det_covs_list$visit_num[i]
  if (length(site_idx) > 0) {
    doy_array   [site_idx, year_idx, visit_idx] <- det_covs_list$day_of_year_std[i]
    tod_array   [site_idx, year_idx, visit_idx] <- det_covs_list$time_of_day_std[i]
    sensor_array[site_idx, year_idx, visit_idx] <- det_covs_list$sensor_binary[i]
  }
}

det_covs_final <- list(
  day_of_year = doy_array,
  time_of_day = tod_array,
  sensor      = sensor_array
)

# =============================================================================
# MCMC SETTINGS, PRIORS, AND SHARED FORMULAS
# =============================================================================
# 3 chains x 40,000 iterations (1,600 batches x 25), 20,000 burn-in, thin 5 -> 12,000 posterior samples for inference.

n_batch      <- 1600
batch_length <- 25
n_burn       <- 20000
n_thin       <- 5
n_chains     <- 3

priors <- list(
  beta.normal  = list(mean = 0, var = 2.72),   # occupancy coefficients
  alpha.normal = list(mean = 0, var = 2.72),   # detection coefficients
  sigma.sq.psi = c(2, 1),                       # inverse-gamma(2, 1)
  phi.unif     = c(0.0006, 0.006)               # 3/5000 to 3/500 (effective range)
)

tuning <- list(phi = 0.5, sigma.sq.psi = 0.5)

# Detection formula is the same for every model.
det_formula <- ~ day_of_year + time_of_day + sensor

# helper so all four models are fit with identical settings; only the
# occupancy covariates, formula, and number of beta parameters change
fit_model <- function(occ_covs, occ_formula, n_beta) {
  stPGOcc(
    occ.formula  = occ_formula,
    det.formula  = det_formula,
    data         = list(y = y, coords = coords,
                        occ.covs = occ_covs, det.covs = det_covs_final),
    inits        = list(beta = rnorm(n_beta, 0, 0.5), alpha = rnorm(4, 0, 0.5),
                        sigma.sq.psi = runif(1, 0.5, 2), phi = runif(1, 0.001, 0.004)),
    priors       = priors,
    tuning       = tuning,
    cov.model    = "exponential",
    NNGP         = TRUE,
    n.neighbors  = 5,
    n.batch      = n_batch,
    batch.length = batch_length,
    n.burn       = n_burn,
    n.thin       = n_thin,
    n.chains     = n_chains,
    n.report     = 200,
    verbose      = TRUE
  )
}

# =============================================================================
# FIT THE FOUR MODELS
# =============================================================================

cat("\nModel 1: habitat amount only\n")
model_1 <- fit_model(list(prop_habitat = site_covs$prop_habitat_std),
                     ~ prop_habitat, n_beta = 2)
saveRDS(model_1, file.path(MODEL_DIR,
        sprintf("model1_habitat_only_def%s_%sm.rds", DEFINITION, SCALE)))

cat("\nModel 2: edge density only\n")
model_2 <- fit_model(list(edge_density = site_covs$edge_density_std),
                     ~ edge_density, n_beta = 2)
saveRDS(model_2, file.path(MODEL_DIR,
        sprintf("model2_edge_only_def%s_%sm.rds", DEFINITION, SCALE)))

cat("\nModel 3: habitat + edge (additive)\n")
model_3 <- fit_model(list(prop_habitat = site_covs$prop_habitat_std,
                          edge_density = site_covs$edge_density_std),
                     ~ prop_habitat + edge_density, n_beta = 3)
saveRDS(model_3, file.path(MODEL_DIR,
        sprintf("model3_habitat_edge_add_def%s_%sm.rds", DEFINITION, SCALE)))

cat("\nModel 4: habitat * edge (interaction)\n")
model_4 <- fit_model(list(prop_habitat = site_covs$prop_habitat_std,
                          edge_density = site_covs$edge_density_std),
                     ~ prop_habitat * edge_density, n_beta = 4)
saveRDS(model_4, file.path(MODEL_DIR,
        sprintf("model4_habitat_edge_interact_def%s_%sm.rds", DEFINITION, SCALE)))

# =============================================================================
# MODEL COMPARISON (WAIC)
# =============================================================================
# Lower WAIC = better expected predictive accuracy. delta_WAIC is the difference
# from the best model; models within 2 units are treated as equivalent. The
# weight is the relative support (Akaike-style weights on the WAIC scale).

waic_1 <- waicOcc(model_1, by.sp = FALSE)
waic_2 <- waicOcc(model_2, by.sp = FALSE)
waic_3 <- waicOcc(model_3, by.sp = FALSE)
waic_4 <- waicOcc(model_4, by.sp = FALSE)

waic_table <- data.frame(
  Model = c("M1: Habitat only", "M2: Edge only",
            "M3: Habitat + Edge (add)", "M4: Habitat x Edge (int)"),
  elpd  = c(waic_1[1], waic_2[1], waic_3[1], waic_4[1]),
  pD    = c(waic_1[2], waic_2[2], waic_3[2], waic_4[2]),
  WAIC  = c(waic_1[3], waic_2[3], waic_3[3], waic_4[3])
) %>%
  mutate(delta_WAIC = WAIC - min(WAIC),
         weight     = exp(-0.5 * delta_WAIC) / sum(exp(-0.5 * delta_WAIC))) %>%
  arrange(WAIC)

cat("\nWAIC comparison:\n")
print(waic_table, digits = 3)

write.csv(waic_table,
          file.path(TABLE_DIR,
                    sprintf("waic_comparison_def%s_%sm.csv", DEFINITION, SCALE)),
          row.names = FALSE)

# Posterior summaries (means, 95% CRIs, R-hat) for each model:
summary(model_1)
summary(model_2)
summary(model_3)
summary(model_4)

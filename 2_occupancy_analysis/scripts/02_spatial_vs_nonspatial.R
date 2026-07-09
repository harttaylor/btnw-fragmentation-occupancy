# =============================================================================
# 02  DO WE NEED A SPATIAL MODEL?
# =============================================================================
# Fits a NON-spatial occupancy model (tPGOcc), extracts its site-level residuals,
# and tests them for spatial autocorrelation with Moran's I (at several distance
# thresholds) and an empirical variogram. Significant residual autocorrelation is
# what justifies using the spatial model (stPGOcc) in scripts 03 onwards.
#
# Uses the primary analysis dataset (Definition C, 500 m) as the test case.
#
# INPUTS : data/defC500_glmmdata.csv, data/2yearvisitmatrix.csv,
#          data/visitsfordetcovs.csv
# OUTPUTS: outputs/tables/spatial_autocorrelation_results.csv
# =============================================================================

library(spOccupancy)
library(tidyverse)
library(spdep)     # Moran's I

TABLE_DIR <- "outputs/tables"
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD DATA
# =============================================================================
full_data    <- read.csv("data/defC500_glmmdata.csv")
visit_matrix <- read.csv("data/2yearvisitmatrix.csv")
det_covs     <- read.csv("data/visitsfordetcovs.csv")

# =============================================================================
# PREPARE ARRAYS  (same structure as script 03)
# =============================================================================
y_matrix   <- as.matrix(visit_matrix[, -1])
n_sites    <- nrow(y_matrix)
n_years    <- 2
max_visits <- 3

y <- array(NA, dim = c(n_sites, n_years, max_visits))
for (i in 1:n_sites) {
  y[i, 1, ] <- y_matrix[i, 1:3]
  y[i, 2, ] <- y_matrix[i, 4:6]
}

# Site covariates (year-1 values, ordered to match the visit matrix)
site_covs <- full_data %>%
  distinct(gisid, year_id, .keep_all = TRUE) %>%
  filter(year_id == 1) %>%
  select(gisid, sensor, Easting, Northing,
         edge_density_std, prop_habitat_std)
site_covs <- site_covs[match(visit_matrix$gisid, site_covs$gisid), ]
coords    <- as.matrix(site_covs[, c("Easting", "Northing")])

occ_covs <- list(edge_density = site_covs$edge_density_std,
                 prop_habitat = site_covs$prop_habitat_std)

# Detection covariates
det_covs_list <- det_covs %>%
  mutate(
    date_time     = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M"),
    day_of_year   = as.numeric(format(date_time, "%j")),
    time_of_day   = as.numeric(format(date_time, "%H")) +
                    as.numeric(format(date_time, "%M")) / 60,
    sensor_binary = ifelse(sensor == "ARU", 1, 0),
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
  s <- which(visit_matrix$gisid == det_covs_list$gisid[i])
  t <- det_covs_list$year_id[i]; k <- det_covs_list$visit_num[i]
  if (length(s) > 0) {
    doy_array   [s, t, k] <- det_covs_list$day_of_year_std[i]
    tod_array   [s, t, k] <- det_covs_list$time_of_day_std[i]
    sensor_array[s, t, k] <- det_covs_list$sensor_binary[i]
  }
}
det_covs_final <- list(day_of_year = doy_array,
                       time_of_day = tod_array,
                       sensor      = sensor_array)

# =============================================================================
# FIT THE NON-SPATIAL MODEL  (tPGOcc)
# =============================================================================
# Reduced MCMC is fine here - this is only a diagnostic to detect leftover
# spatial structure, not a model we report.

priors <- list(beta.normal  = list(mean = 0, var = 2.72),
               alpha.normal = list(mean = 0, var = 2.72))

model_nonspatial <- tPGOcc(
  occ.formula  = ~ edge_density + prop_habitat,
  det.formula  = ~ day_of_year + time_of_day + sensor,
  data         = list(y = y, occ.covs = occ_covs, det.covs = det_covs_final),
  n.batch      = 800, batch.length = 25,
  inits        = list(beta = rnorm(3, 0, 0.5), alpha = rnorm(4, 0, 0.5)),
  priors       = priors,
  n.burn       = 10000, n.thin = 5, n.chains = 3,
  n.report     = 200, verbose = TRUE
)
summary(model_nonspatial)

# =============================================================================
# SITE-LEVEL RESIDUALS
# =============================================================================
# Residual = observed occupancy (ever detected at the site?) minus the posterior
# mean occupancy probability, averaged across years.

psi_mean <- apply(model_nonspatial$psi.samples, c(2, 3), mean)  # site x year
psi_site <- rowMeans(psi_mean)

y_site <- apply(y, 1, function(x) max(x, na.rm = TRUE))
y_site[is.infinite(y_site)] <- NA

residuals_raw <- y_site - psi_site

# =============================================================================
# MORAN'S I ACROSS DISTANCE THRESHOLDS
# =============================================================================
# Test whether residuals from nearby sites are more similar than expected. A
# significant positive Moran's I means spatial structure the non-spatial model
# failed to capture -> a spatial model is warranted.

valid            <- !is.na(residuals_raw)
coords_valid     <- coords[valid, ]
residuals_valid  <- residuals_raw[valid]
dist_thresholds  <- c(500, 1000, 2500, 5000)   # metres

moran_results <- data.frame()
for (d in dist_thresholds) {
  nb <- dnearneigh(coords_valid, 0, d)
  if (sum(card(nb) > 0) < 10) {
    cat(sprintf("Distance %d m: too few neighbours, skipping\n", d)); next
  }
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  mt <- moran.test(residuals_valid, lw, zero.policy = TRUE)

  moran_results <- rbind(moran_results, data.frame(
    threshold_m = d,
    moran_I     = unname(mt$estimate[1]),
    expected_I  = unname(mt$estimate[2]),
    p_value     = mt$p.value
  ))
  cat(sprintf("Distance %d m: Moran's I = %.3f, p = %.4f\n",
              d, mt$estimate[1], mt$p.value))
}

# =============================================================================
# EMPIRICAL VARIOGRAM
# =============================================================================
# Semivariance rising with distance is another signature of positive spatial
# autocorrelation in the residuals.

dist_matrix <- as.matrix(dist(coords_valid))
diff_matrix <- outer(residuals_valid, residuals_valid, function(x, y) (x - y)^2)
upper       <- upper.tri(dist_matrix)

breaks <- seq(0, 50000, by = 2500)
bins   <- cut(dist_matrix[upper], breaks = breaks, labels = FALSE)

variogram_data <- data.frame(
  distance_km  = (breaks[-length(breaks)] + 1250) / 1000,
  semivariance = tapply(diff_matrix[upper], bins, function(x) mean(x, na.rm = TRUE) / 2),
  n_pairs      = tapply(diff_matrix[upper], bins, length)
)
variogram_data <- variogram_data[complete.cases(variogram_data), ]
print(head(variogram_data, 10), digits = 3)

# =============================================================================
# SAVE
# =============================================================================
write.csv(moran_results,
          file.path(TABLE_DIR, "spatial_autocorrelation_results.csv"),
          row.names = FALSE)
cat("\nSaved:", file.path(TABLE_DIR, "spatial_autocorrelation_results.csv"), "\n")

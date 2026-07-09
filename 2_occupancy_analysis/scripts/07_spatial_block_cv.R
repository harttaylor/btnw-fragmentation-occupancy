# =============================================================================
# 07  SPATIAL BLOCK CROSS-VALIDATION  (out-of-sample validation)
# =============================================================================
# Validates the top-supported model with 5-fold spatial block cross-validation
# Sites are grouped into 20 x 20 km blocks and whole blocks are assigned to folds
# so training and test sites are spatially independent 
#
# PREREQUISITE: run script 03 first (down through "MCMC SETTINGS") so these
#   objects exist: y, coords, site_covs, det_covs_final, det_formula, priors,
#   tuning, n_years, max_visits.
#
# OUTPUTS: outputs/tables/cv_results_<MODEL_TAG>.csv
#          outputs/tables/fold_assignment_<MODEL_TAG>.csv
# =============================================================================

library(spOccupancy)
library(pROC)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
BLOCK_SIZE_M <- 20000      # 20 km blocks (exceeds the ~4 km spatial range)
K            <- 5          # number of folds
CV_SEED      <- 100        # reproducible fold assignment

# Top-supported model: habitat x edge interaction at 7 ha (150 m), Definition C.
cv_occ_formula <- ~ prop_habitat * edge_density
cv_occ_covs    <- list(prop_habitat = site_covs$prop_habitat_std,
                       edge_density = site_covs$edge_density_std)

# Reduced MCMC for CV (the model is fit K times); convergence checked per fold.
cv_n_batch      <- 1200
cv_batch_length <- 25
cv_n_burn       <- 15000
cv_n_thin       <- 5
cv_n_chains     <- 3

TABLE_DIR <- "outputs/tables"
MODEL_TAG <- "top_model_habitat_x_edge_defC_7ha"
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Assign sites to spatial blocks
# -----------------------------------------------------------------------------
block_x  <- floor(coords[, 1] / BLOCK_SIZE_M)
block_y  <- floor(coords[, 2] / BLOCK_SIZE_M)
block_id <- paste(block_x, block_y, sep = "_")

unique_blocks <- unique(block_id)
cat("Assigned", nrow(coords), "sites to", length(unique_blocks),
    "blocks of", BLOCK_SIZE_M / 1000, "km\n")

# -----------------------------------------------------------------------------
# 2. Assign blocks to folds (size-balanced greedy assignment)
# -----------------------------------------------------------------------------
# Sampling is strongly clustered (one block has >100 sites), so random
# assignment gives very unbalanced folds. We place each block (largest first)
# into the fold with the fewest sites so far. Whole blocks stay in one fold
# (preserving spatial independence) while fold sizes stay roughly balanced

set.seed(CV_SEED)
block_sizes <- table(block_id)
block_order <- names(block_sizes)[order(-as.numeric(block_sizes),
                                        sample(seq_along(block_sizes)))]

fold_sites_count <- rep(0, K)
block_fold       <- setNames(integer(length(unique_blocks)), unique_blocks)
for (b in block_order) {
  target <- which.min(fold_sites_count)
  block_fold[b] <- target
  fold_sites_count[target] <- fold_sites_count[target] + as.numeric(block_sizes[b])
}
site_fold <- block_fold[block_id]

cat("Sites per fold:\n"); print(table(site_fold))

write.csv(
  data.frame(site_index = seq_len(nrow(coords)),
             Easting = coords[, 1], Northing = coords[, 2],
             block_id = block_id, fold = as.integer(site_fold)),
  file.path(TABLE_DIR, paste0("fold_assignment_", MODEL_TAG, ".csv")),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 3. Helper: deviance + AUC on held-out visits
# -----------------------------------------------------------------------------
# For each visit the predicted detection probability is psi * p. Deviance is
# -2 * sum(log-likelihood) over all observed visits; AUC measures discrimination
compute_test_metrics <- function(pred_psi_samples, pred_p_samples, y_test) {
  psi_mean <- apply(pred_psi_samples, c(2, 3),    mean)   # site x year
  p_mean   <- apply(pred_p_samples,   c(2, 3, 4), mean)   # site x year x visit

  pred_y <- array(NA, dim = dim(y_test))
  for (k in seq_len(dim(y_test)[3])) pred_y[, , k] <- psi_mean * p_mean[, , k]

  y_vec    <- as.vector(y_test)
  pred_vec <- as.vector(pred_y)
  keep     <- !is.na(y_vec) & !is.na(pred_vec)
  y_vec    <- y_vec[keep]
  pred_vec <- pmin(pmax(pred_vec[keep], 1e-10), 1 - 1e-10)

  log_lik  <- y_vec * log(pred_vec) + (1 - y_vec) * log(1 - pred_vec)
  list(deviance         = -2 * sum(log_lik),
       n_obs            = length(y_vec),
       deviance_per_obs = -2 * sum(log_lik) / length(y_vec),
       auc              = as.numeric(pROC::auc(y_vec, pred_vec, quiet = TRUE)))
}

# -----------------------------------------------------------------------------
# 4. Cross-validation loop
# -----------------------------------------------------------------------------
n_beta_cv  <- ncol(model.matrix(cv_occ_formula, data = as.data.frame(cv_occ_covs)))
n_alpha_cv <- 4   # intercept + day + time + sensor

make_inits <- function() list(beta  = rnorm(n_beta_cv, 0, 0.5),
                              alpha = rnorm(n_alpha_cv, 0, 0.5),
                              sigma.sq.psi = runif(1, 0.5, 2),
                              phi = runif(1, 0.001, 0.004))

cv_results <- vector("list", K)

for (fold in seq_len(K)) {
  cat("\n---- Fold", fold, "of", K, "----\n")
  test_idx  <- which(site_fold == fold)
  train_idx <- which(site_fold != fold)
  cat("Train:", length(train_idx), " | Test:", length(test_idx), "\n")

  data_train <- list(
    y        = y[train_idx, , , drop = FALSE],
    coords   = coords[train_idx, , drop = FALSE],
    occ.covs = lapply(cv_occ_covs, function(x) x[train_idx]),
    det.covs = lapply(det_covs_final, function(a) a[train_idx, , , drop = FALSE])
  )

  set.seed(CV_SEED + fold)
  fold_model <- stPGOcc(
    occ.formula = cv_occ_formula, det.formula = det_formula, data = data_train,
    n.batch = cv_n_batch, batch.length = cv_batch_length, inits = make_inits(),
    priors = priors, tuning = tuning, cov.model = "exponential",
    NNGP = TRUE, n.neighbors = 5,
    n.burn = cv_n_burn, n.thin = cv_n_thin, n.chains = cv_n_chains,
    n.report = 200, verbose = FALSE
  )

  rhat_max <- max(c(fold_model$rhat$beta, fold_model$rhat$alpha), na.rm = TRUE)
  cat("Max Rhat:", round(rhat_max, 3), "\n")
  if (rhat_max > 1.1) warning("Fold ", fold, ": Rhat > 1.1.")

  # Occupancy design matrix at test sites, replicated across years
  X.0       <- model.matrix(cv_occ_formula,
                            data = as.data.frame(lapply(cv_occ_covs, function(x) x[test_idx])))
  X.0_array <- array(NA, dim = c(length(test_idx), n_years, ncol(X.0)))
  for (t in seq_len(n_years)) X.0_array[, t, ] <- X.0

  psi_pred <- predict(fold_model, X.0 = X.0_array,
                      coords.0 = coords[test_idx, , drop = FALSE],
                      t.cols = 1:n_years, type = "occupancy", verbose = FALSE)
  pred_psi_samples <- psi_pred$psi.0.samples

  # Detection probabilities at test sites from the posterior alpha draws
  alpha_samples <- fold_model$alpha.samples
  n_samples     <- nrow(alpha_samples)
  doy    <- det_covs_final$day_of_year[test_idx, , , drop = FALSE]
  tod    <- det_covs_final$time_of_day[test_idx, , , drop = FALSE]
  sensor <- det_covs_final$sensor     [test_idx, , , drop = FALSE]

  pred_p_samples <- array(NA, dim = c(n_samples, length(test_idx), n_years, max_visits))
  for (s in seq_len(n_samples))
    for (t in seq_len(n_years))
      for (k in seq_len(max_visits)) {
        eta <- alpha_samples[s, 1] + alpha_samples[s, 2] * doy[, t, k] +
               alpha_samples[s, 3] * tod[, t, k] + alpha_samples[s, 4] * sensor[, t, k]
        pred_p_samples[s, , t, k] <- plogis(eta)
      }

  metrics <- compute_test_metrics(pred_psi_samples, pred_p_samples,
                                  y[test_idx, , , drop = FALSE])
  cv_results[[fold]] <- c(list(fold = fold,
                               n_train = length(train_idx),
                               n_test  = length(test_idx)), metrics)
  cat("Fold", fold, "AUC:", round(metrics$auc, 3),
      " | Deviance/obs:", round(metrics$deviance_per_obs, 3), "\n")

  rm(fold_model, pred_psi_samples, pred_p_samples); gc(verbose = FALSE)
}

# -----------------------------------------------------------------------------
# Summarise and save
# -----------------------------------------------------------------------------
cv_summary <- do.call(rbind, lapply(cv_results, function(x)
  data.frame(fold = x$fold, n_train = x$n_train, n_test = x$n_test,
             n_obs = x$n_obs, deviance = x$deviance,
             deviance_per_obs = x$deviance_per_obs, auc = x$auc)))

print(cv_summary, digits = 4)
cat(sprintf("\nMean AUC: %.3f (SD %.3f)\n",
            mean(cv_summary$auc), sd(cv_summary$auc)))

write.csv(cv_summary,
          file.path(TABLE_DIR, paste0("cv_results_", MODEL_TAG, ".csv")),
          row.names = FALSE)
cat("Saved:", file.path(TABLE_DIR, paste0("cv_results_", MODEL_TAG, ".csv")), "\n")

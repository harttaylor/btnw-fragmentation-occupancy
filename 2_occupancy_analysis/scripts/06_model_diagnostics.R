# =============================================================================
# 06  MODEL DIAGNOSTICS  (for the top-supported model)
# =============================================================================
#
# INPUTS  : one fitted model from script 03.
# OUTPUTS : outputs/figures/diagnostics_traceplots.png
#           outputs/figures/diagnostics_residuals.png
#           (printed to console: PPC Bayesian p-value, Moran's I, coefficients)
# =============================================================================

library(spOccupancy)
library(spdep)     # Moran's I

DEFINITION <- "A"
SCALE      <- "150"
MODEL_DIR  <- "outputs/models"
FIG_DIR    <- "outputs/figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Load the additive (habitat + edge) model - the one reported at this scale
model <- readRDS(file.path(MODEL_DIR,
          sprintf("model3_habitat_edge_add_def%s_%sm.rds", DEFINITION, SCALE)))

summary(model)

# =============================================================================
# 1. TRACE PLOTS
# =============================================================================

png(file.path(FIG_DIR, "diagnostics_traceplots.png"),
    width = 12, height = 10, units = "in", res = 150)
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

plot(model$beta.samples[, 1],  type = "l", main = "Intercept (occupancy)",
     xlab = "Iteration", ylab = "Value", col = rgb(0, 0, 0, 0.5))
plot(model$beta.samples[, 2],  type = "l", main = "prop_habitat",
     xlab = "Iteration", ylab = "Value", col = rgb(0, 0, 0, 0.5))
plot(model$beta.samples[, 3],  type = "l", main = "edge_density",
     xlab = "Iteration", ylab = "Value", col = rgb(0, 0, 0, 0.5))

plot(model$alpha.samples[, 1], type = "l", main = "Intercept (detection)",
     xlab = "Iteration", ylab = "Value", col = rgb(0, 0, 0, 0.5))
plot(model$alpha.samples[, 2], type = "l", main = "day_of_year",
     xlab = "Iteration", ylab = "Value", col = rgb(0, 0, 0, 0.5))
plot(model$alpha.samples[, 3], type = "l", main = "time_of_day",
     xlab = "Iteration", ylab = "Value", col = rgb(0, 0, 0, 0.5))

plot(model$theta.samples[, 1], type = "l", main = "sigma.sq (spatial variance)",
     xlab = "Iteration", ylab = "Value", col = rgb(0, 0, 0, 0.5))
plot(model$theta.samples[, 2], type = "l", main = "phi (spatial decay)",
     xlab = "Iteration", ylab = "Value", col = rgb(0, 0, 0, 0.5))
plot(3 / model$theta.samples[, 2], type = "l", main = "Effective range (m)",
     xlab = "Iteration", ylab = "Range (m)", col = rgb(0, 0, 0, 0.5))

dev.off()
cat("Trace plots saved to", file.path(FIG_DIR, "diagnostics_traceplots.png"), "\n")

# =============================================================================
# 2. POSTERIOR PREDICTIVE CHECK  (Freeman-Tukey)
# =============================================================================
# ppcOcc simulates data from the fitted model and compares a discrepancy
# statistic between observed and simulated data
# summary() reports a Bayesian p-value where values near 0.5 indicate the model reproduces the data well
# values near 0 or 1 indicate poor fit

ppc <- ppcOcc(model, fit.stat = "freeman-tukey", group = 1)
summary(ppc)

# =============================================================================
# 3. RESIDUAL ANALYSIS
# =============================================================================
# Residual = observed occupancy (was the species ever detected at the site?)
# minus the posterior mean occupancy probability, computed per year and pooled.

psi_mean <- apply(model$psi.samples, c(2, 3), mean)   # site x year
y        <- model$y
n_years  <- dim(y)[2]

resid_by_year <- lapply(seq_len(n_years), function(t) {
  y_site <- apply(y[, t, , drop = FALSE], 1, function(x)
                  if (all(is.na(x))) NA else max(x, na.rm = TRUE))
  y_site - psi_mean[, t]
})

all_resid <- unlist(resid_by_year)
all_resid <- all_resid[!is.na(all_resid)]

cat(sprintf("\nResiduals: mean = %.4f (want ~0), SD = %.4f, range [%.3f, %.3f]\n",
            mean(all_resid), sd(all_resid), min(all_resid), max(all_resid)))

png(file.path(FIG_DIR, "diagnostics_residuals.png"),
    width = 12, height = 5, units = "in", res = 150)
par(mfrow = c(1, 2))

plot(as.vector(psi_mean), all_resid,
     xlab = "Fitted occupancy probability", ylab = "Residual (observed - fitted)",
     main = "Residuals vs fitted", pch = 19, col = rgb(0, 0, 0, 0.3))
abline(h = 0, col = "red", lwd = 2)
lines(lowess(as.vector(psi_mean), all_resid, f = 0.5), col = "blue", lwd = 2)

hist(all_resid, breaks = 30, freq = FALSE, col = "grey80", border = "white",
     main = "Distribution of residuals", xlab = "Residual")
curve(dnorm(x, mean(all_resid), sd(all_resid)), add = TRUE,
      col = "blue", lwd = 2, lty = 2)

dev.off()
cat("Residual plots saved to", file.path(FIG_DIR, "diagnostics_residuals.png"), "\n")

# =============================================================================
# 4. SPATIAL RESIDUAL CHECK  (Moran's I)
# =============================================================================
# If the spatial random effect did its job, the residuals should show no
# significant spatial autocorrelation (Moran's I p-value > 0.05).

coords     <- model$coords
resid_site <- rowMeans(do.call(cbind, resid_by_year), na.rm = TRUE)
ok         <- !is.na(resid_site)

nb <- dnearneigh(coords[ok, ], 0, 2500)             # neighbours within 2.5 km
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
mi <- moran.test(resid_site[ok], lw, zero.policy = TRUE)

cat(sprintf("\nMoran's I on residuals (2.5 km): I = %.4f, p = %.4f\n",
            mi$estimate[1], mi$p.value))
cat(if (mi$p.value > 0.05)
      "  -> No significant residual autocorrelation (spatial model is adequate).\n"
    else
      "  -> Residual autocorrelation remains; consider more neighbours / wider phi.\n")

# =============================================================================
# 5. COEFFICIENTS IN PLAIN LANGUAGE
# =============================================================================
# We summarise straight from the posterior sample matrices (which always carry
# the parameter names as column names). A term is "well supported" when its 95%
# credible interval excludes zero (flagged ***). For those we also give the odds
# ratio, exp(beta): a 1-SD change multiplies the odds of occupancy (or detection)
# by that factor.

report_samples <- function(samples, kind, prob = FALSE) {
  cat("\n", kind, " (log-odds scale)\n", sep = "")
  cat(strrep("-", 50), "\n")
  for (nm in colnames(samples)) {
    draws <- samples[, nm]
    m  <- mean(draws)
    lo <- quantile(draws, 0.025); hi <- quantile(draws, 0.975)
    star <- if (sign(lo) == sign(hi)) " ***" else ""
    cat(sprintf("  %-28s %6.3f  (95%% CRI: %6.3f, %6.3f)%s\n", nm, m, lo, hi, star))
    if (prob && grepl("Intercept", nm))
      cat(sprintf("      -> baseline probability: %.1f%%\n", plogis(m) * 100))
    else if (star != "" && !grepl("Intercept", nm))
      cat(sprintf("      -> odds ratio: %.2f\n", exp(m)))
  }
}

report_samples(model$beta.samples,  "OCCUPANCY parameters", prob = TRUE)
report_samples(model$alpha.samples, "DETECTION parameters", prob = TRUE)

# ---- spatial parameters + effective range -----------------------------------
# For an exponential NNGP, theta.samples holds sigma.sq (column 1) and phi
# (column 2). The effective range (distance at which spatial correlation decays
# to ~5%) is 3/phi. (Confirm the columns with colnames(model$theta.samples).)
sig <- model$theta.samples[, 1]   # sigma.sq (spatial variance)
phi <- model$theta.samples[, 2]   # phi      (spatial decay)
rng <- 3 / phi

cat("\nSPATIAL parameters\n"); cat(strrep("-", 50), "\n")
cat(sprintf("  sigma.sq (variance):   %.2f  (95%% CRI: %.2f, %.2f)\n",
            mean(sig), quantile(sig, 0.025), quantile(sig, 0.975)))
cat(sprintf("  phi (decay):           %.4f (95%% CRI: %.4f, %.4f)\n",
            mean(phi), quantile(phi, 0.025), quantile(phi, 0.975)))
cat(sprintf("  effective range 3/phi: %.0f m  (95%% CRI: %.0f, %.0f m)\n",
            mean(rng), quantile(rng, 0.025), quantile(rng, 0.975)))
cat("      -> distance at which spatial correlation decays to ~5%\n")

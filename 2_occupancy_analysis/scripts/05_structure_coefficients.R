# =============================================================================
# 05  STRUCTURE COEFFICIENTS  (Supplementary Table S3)
# =============================================================================
# Compare the posterior regression SLOPE (beta), and the STRUCTURE COEFFICIENT (SC) 
# for every occupancy predictor in every candidate model
#
# WHY WE DO THIS
#   Habitat amount and fragmentation are geometrically coupled, so a partial
#   slope can flip sign purely because of collinearity (a "suppression" effect)
#   rather than because of a real association with occupancy. A structure
#   coefficient is the correlation between a predictor and the model's fitted
#   values (Ray-Mukherjee et al. 2014). Its sign reflects the predictor's plain
#   bivariate relationship with the response
#       slope sign  ==  SC sign   ->  effect is interpretable
#       slope sign  !=  SC sign   ->  slope is driven by collinearity; DO NOT interpret that term
#
# Bayesian version of a structure coefficient:
#   For each posterior draw d,
#       eta_d = X %*% beta_d          # the fitted linear predictor (logit scale)
#       SC_j  = cor( eta_d , X[, j] ) # correlation of predictor j with the fit
#   then summarise SC_j across all draws (posterior mean + 95% interval),
#   and do the same for the slope beta_j
#
# INPUTS  : the four fitted models from script 03, for this stage.
# OUTPUTS : outputs/tables/structure_coefs_defC_150m.csv  (one row per term)
# =============================================================================

library(spOccupancy) 
library(tidyverse)

DEFINITION <- "A"
SCALE      <- "150"    

MODEL_DIR  <- "outputs/models"
OUT_DIR    <- "outputs/tables"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD FITTED MODELS
# =============================================================================
# Model 4 (the interaction) is the one whose interpretability we are most concerned about

model_1 <- readRDS(file.path(MODEL_DIR,
            sprintf("model1_habitat_only_def%s_%sm.rds",         DEFINITION, SCALE)))
model_2 <- readRDS(file.path(MODEL_DIR,
            sprintf("model2_edge_only_def%s_%sm.rds",            DEFINITION, SCALE)))
model_3 <- readRDS(file.path(MODEL_DIR,
            sprintf("model3_habitat_edge_add_def%s_%sm.rds",     DEFINITION, SCALE)))
model_4 <- readRDS(file.path(MODEL_DIR,
            sprintf("model4_habitat_edge_interact_def%s_%sm.rds", DEFINITION, SCALE)))

# =============================================================================
# FLATTEN THE DESIGN MATRIX
# =============================================================================
# spOccupancy stores the occupancy design matrix as a 3D array with dimensions
# [site, year, predictor]. To correlate a predictor column with the fitted
# values we need a 2D matrix [ (site*year) , predictor ]; stack
# the year-1 slice on top of the year-2 slice.

flatten_X <- function(X) {
  n_years <- dim(X)[2]
  X_flat  <- do.call(rbind, lapply(seq_len(n_years), function(t) X[, t, ]))
  colnames(X_flat) <- dimnames(X)[[3]]   # keep the predictor names
  X_flat
}

# ---- (optional) sanity check: flattened vs year-1-only give the same SC ------
X_full  <- model_1$X
X_flat  <- flatten_X(X_full)                     # stacked years  (1016 x 2)
X_year1 <- X_full[, 1, ]                          # year 1 only     (508 x 2)
colnames(X_year1) <- dimnames(X_full)[[3]]

beta1   <- model_1$beta.samples                   # draws x predictors
j_hab   <- which(colnames(X_flat) == "prop_habitat")

sc_flat  <- apply(beta1 %*% t(X_flat),  1, function(e) cor(e, X_flat[,  j_hab]))
sc_year1 <- apply(beta1 %*% t(X_year1), 1, function(e) cor(e, X_year1[, j_hab]))
cat("Sanity check (should match):  flattened =", round(mean(sc_flat), 4),
    " | year 1 only =", round(mean(sc_year1), 4), "\n\n")

# =============================================================================
# HELPER THAT DOES ONE TERM
# =============================================================================
# Given a flattened design matrix, the posterior slope draws, and a column
# index j, this returns the slope summary, the SC summary, and whether their
# signs agree

sc_for_term <- function(X_flat, beta_samples, j) {
  eta <- beta_samples %*% t(X_flat)                    # draws x (site*year)
  sc  <- apply(eta, 1, function(eta_d) cor(eta_d, X_flat[, j]))
  b   <- beta_samples[, j]
  pd  <- max(mean(b > 0), mean(b < 0))                 # probability of direction
  data.frame(
    term      = colnames(X_flat)[j],
    beta_mean = mean(b),
    beta_lwr  = unname(quantile(b, 0.025)),
    beta_upr  = unname(quantile(b, 0.975)),
    pd        = pd,
    sc_mean   = mean(sc),
    sc_lwr    = unname(quantile(sc, 0.025)),
    sc_upr    = unname(quantile(sc, 0.975)),
    row.names = NULL
  )
}

# =============================================================================
# RUN EACH TERM FOR EACH MODEL
# =============================================================================

# ---- MODEL 1: habitat amount only -------------------------------------------
X1 <- flatten_X(model_1$X)
b1 <- model_1$beta.samples
cat("Predictors in Model 1:", paste(colnames(X1), collapse = ", "), "\n")

m1_habitat <- sc_for_term(X1, b1, which(colnames(X1) == "prop_habitat"))

# ---- MODEL 2: edge density only ---------------------------------------------
X2 <- flatten_X(model_2$X)
b2 <- model_2$beta.samples

m2_edge <- sc_for_term(X2, b2, which(colnames(X2) == "edge_density"))

# ---- MODEL 3: habitat + edge (additive) -------------------------------------
X3 <- flatten_X(model_3$X)
b3 <- model_3$beta.samples

m3_habitat <- sc_for_term(X3, b3, which(colnames(X3) == "prop_habitat"))
m3_edge    <- sc_for_term(X3, b3, which(colnames(X3) == "edge_density"))

# ---- MODEL 4: habitat * edge (interaction) ----------------------------------
# The flattened design matrix already contains the interaction product column
# "prop_habitat:edge_density", so its SC is computed exactly like a main effect.
X4 <- flatten_X(model_4$X)
b4 <- model_4$beta.samples
cat("Predictors in Model 4:", paste(colnames(X4), collapse = ", "), "\n\n")

m4_habitat     <- sc_for_term(X4, b4, which(colnames(X4) == "prop_habitat"))
m4_edge        <- sc_for_term(X4, b4, which(colnames(X4) == "edge_density"))
m4_interaction <- sc_for_term(X4, b4, which(colnames(X4) == "prop_habitat:edge_density"))

# =============================================================================
# ASSEMBLE RESULTS TABLE
# =============================================================================
# Tag each row with its model, then add the interpretation flags: do the signs
# agree, and does the 95% credible interval of the slope exclude zero?

results <- bind_rows(
  cbind(model = "M1 habitat only",   m1_habitat),
  cbind(model = "M2 edge only",      m2_edge),
  cbind(model = "M3 habitat + edge", m3_habitat),
  cbind(model = "M3 habitat + edge", m3_edge),
  cbind(model = "M4 habitat * edge", m4_habitat),
  cbind(model = "M4 habitat * edge", m4_edge),
  cbind(model = "M4 habitat * edge", m4_interaction)
) %>%
  mutate(
    definition         = DEFINITION,
    scale              = SCALE,
    slope_sign         = ifelse(beta_mean > 0, "+", "-"),
    sc_sign            = ifelse(sc_mean   > 0, "+", "-"),
    signs_agree        = slope_sign == sc_sign,   # TRUE  = interpretable
    beta_excludes_zero = !(beta_lwr < 0 & beta_upr > 0)
  ) %>%
  select(definition, scale, model, term,
         beta_mean, beta_lwr, beta_upr, slope_sign, beta_excludes_zero,
         sc_mean, sc_lwr, sc_upr, sc_sign,
         signs_agree, pd)

print(results, digits = 3)


# Any term whose slope and SC disagree in sign is NOT interpretable
not_interpretable <- results %>% filter(!signs_agree)

if (nrow(not_interpretable) == 0) {
  cat("  none - every term is sign-concordant.\n")
} else {
  not_interpretable %>%
    transmute(model, term,
              beta = round(beta_mean, 2), sc = round(sc_mean, 2)) %>%
    print(row.names = FALSE)
}

# ---- Save --------------------------------------------------------------------
write.csv(results,
          file.path(OUT_DIR, sprintf("structure_coefs_def%s_%sm.csv",
                                     DEFINITION, SCALE)),
          row.names = FALSE)


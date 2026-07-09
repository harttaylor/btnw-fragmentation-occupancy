# =============================================================================
# 04  PREDICTOR COUPLING DIAGNOSTIC  (Supplementary Table S1)
# =============================================================================
# This script measures how tightly each configuration metric (edge density, mean patch area, aggregation index)
# is coupled to habitat amount, using both a linear and a non-linear measure.
#
# WHY WE DO THIS
#   Habitat amount and fragmentation are often geometrically coupled, and that
#   coupling is frequently NON-LINEAR. Linear tools (Pearson r, VIF) can miss
#   non-linear coupling entirely. A generalised additive model (GAM) fits a
#   smooth curve, so its "deviance explained" captures coupling of any shape.
#   If the GAM deviance explained is much larger than the squared Pearson r,
#   the two variables are non-linearly entangled.
#
# We fit the smooth in BOTH directions for each metric:
#       metric ~ s(amount) ->  does habitat amount shape the metric?
#       amount ~ s(metric) ->  how much of habitat amount does the metric
#                                reveal?  (this is the entanglement we care
#                                about: LOWER = the metric is MORE separable
#                                from amount, and therefore a cleaner
#                                fragmentation predictor)
#
# NOTE ON UNDEFINED VALUES
#   Aggregation index and contagion are undefined at sites with no habitat.
#   We assess each metric on COMPLETE CASES (sites where it is defined), i.e.
#   the domain over which configuration actually exists.
#
# INPUTS  : data/defC150_glmmdata.csv   (from script 01)
# OUTPUTS : outputs/tables/coupling_defC_150m.csv   (one row per metric)
# =============================================================================

library(mgcv) # gam()
library(tidyverse)

# Change and re-run for each DEFINITION / SCALE combination.
DEFINITION <- "A"
SCALE <- "150"  # 150 / 500 / 1000

OUT_DIR <- "outputs/tables"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# STEP 1 - LOAD AND STANDARDISE ONE ROW PER SITE
# =============================================================================
# The modeling data has one row per visit but we need one row per site
# Also standardise every predictor (mean 0, SD 1) so the smooths and the Pearson
# r are all on the same scale. mean_patch_area is log-transformed first, matching
# how it enters the occupancy models in script 03.

dat <- read.csv(sprintf("data/def%s%s_glmmdata.csv", DEFINITION, SCALE))

site_dat <- dat %>%
  distinct(gisid, year_id, .keep_all = TRUE) %>%
  filter(year_id == 1) %>%
  mutate(
    prop_habitat_std = scale(prop_habitat)[, 1],
    edge_density_std = scale(edge_density)[, 1],
    mean_patch_area_std = scale(log(mean_patch_area + 0.01))[, 1],
    aggregation_std = scale(aggregation_index)[, 1]   # NA where no habitat
  )

# =============================================================================
# STEP 2 - A SMALL HELPER FOR ONE METRIC
# =============================================================================
# For a single metric it fits both GAM directions on complete cases and returns
# the Pearson r plus the effective degrees of freedom (edf; ~1 means linear,
# >1 means curved) and deviance explained for each direction.

coupling <- function(d, metric_std, amount_std = "prop_habitat_std") {
  cc <- d[complete.cases(d[, c(amount_std, metric_std)]),
          c(amount_std, metric_std)]
  g_metric <- gam(as.formula(sprintf("%s ~ s(%s)", metric_std, amount_std)),
                  data = cc, method = "REML")   # metric ~ s(amount)
  g_amount <- gam(as.formula(sprintf("%s ~ s(%s)", amount_std, metric_std)),
                  data = cc, method = "REML")   # amount ~ s(metric)
  data.frame(
    metric = metric_std,
    n = nrow(cc),
    pearson_r = round(cor(cc[[amount_std]], cc[[metric_std]]), 3),
    edf_metric_on_amount = round(summary(g_metric)$s.table[1, "edf"], 2),
    dev_metric_on_amount = round(summary(g_metric)$dev.expl, 3),
    edf_amount_on_metric = round(summary(g_amount)$s.table[1, "edf"], 2),
    dev_amount_on_metric = round(summary(g_amount)$dev.expl, 3)
  )
}

# =============================================================================
# STEP 3 - RUN FOR EACH METRIC 
# =============================================================================

edge <- coupling(site_dat, "edge_density_std")
mpa  <- coupling(site_dat, "mean_patch_area_std")
agg  <- coupling(site_dat, "aggregation_std")

coupling_tab <- bind_rows(edge, mpa, agg) %>%
  mutate(definition = DEFINITION, scale = SCALE) %>%
  arrange(dev_amount_on_metric) # most separable (lowest coupling) at the top

print(coupling_tab)

# =============================================================================
# STEP 4 - LOOK AT THE KEY PAIR (habitat amount vs edge density)
# =============================================================================
# A picture of the pair we care about most. The GAM smooth (red) shows the shape
# of the relationship; compare it against how a straight line (Pearson r) would
# describe the same cloud.

ggplot(site_dat, aes(prop_habitat_std, edge_density_std)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x),
              colour = "firebrick", fill = "firebrick", alpha = 0.2) +
  labs(x = "Habitat amount (standardized)",
       y = "Edge density (standardized)",
       title = sprintf("Amount vs edge, Def %s / %s m", DEFINITION, SCALE),
       subtitle = sprintf("Pearson r = %.2f   |   GAM edf = %.2f, dev.expl = %.2f",
                          edge$pearson_r,
                          edge$edf_amount_on_metric,
                          edge$dev_amount_on_metric)) +
  theme_bw()

# ---- Save --------------------------------------------------------------------
write.csv(coupling_tab,
          file.path(OUT_DIR, sprintf("coupling_def%s_%sm.csv",
                                     DEFINITION, SCALE)),
          row.names = FALSE)



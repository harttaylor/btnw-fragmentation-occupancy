# =============================================================================
# 09  MAIN-TEXT FIGURES
# =============================================================================
# Reproduces the three R-generated main-text figures:
#   Figure 3  dWAIC by scale and fragmentation definition (bar plot)
#   Figure 4  occupancy coefficient forest plot (habitat, edge, interaction)
#   Figure 5  predicted-occupancy partial effects for habitat and edge, by scale
#
# The plotted values (WAIC deltas, coefficient means and 95% CRIs) are the
# results reported in the manuscript - i.e. the outputs of scripts 03 and 05.
# They are entered directly here so the figures can be regenerated without
# reloading the fitted models. If you re-fit the models, update the numbers in
# the labelled data blocks below.
#
# OUTPUTS: outputs/figures/fig3_dWAIC_by_scale.png
#          outputs/figures/fig4_forest_coefficients.png
#          outputs/figures/fig5_habitat_edge_partial_effects.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(gridExtra)

FIG_DIR <- "outputs/figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SHARED STYLE
# =============================================================================
scale_ar   <- c("150m" = "7 ha", "500m" = "78 ha", "1000m" = "314 ha")
def_cols   <- c("A" = "#F4A07A", "B" = "#A88FC4", "C" = "#56B4AC")
def_labs   <- c("A" = "A: Polygonal", "B" = "B: + Wide linear", "C" = "C: + Narrow linear")
scale_cols <- c("7 ha" = "#9A031E", "78 ha" = "#E69F00", "314 ha" = "#0F4C5C")

theme_clean <- function(b = 12) theme_bw(b) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border     = element_rect(color = "grey60", fill = NA),
        strip.background  = element_blank(),
        strip.text        = element_text(face = "bold"),
        legend.key        = element_blank(),
        plot.tag          = element_text(face = "bold", size = 12))

# Pull the legend grob out of a plot (used to share one legend across panels).
get_legend <- function(p) {
  g <- ggplotGrob(p + theme(legend.position = "bottom"))
  k <- which(grepl("guide-box", vapply(g$grobs, function(x) x$name, "")))
  if (length(k)) g$grobs[[k[1]]] else grid::nullGrob()
}

# =============================================================================
# FIGURE 3 - dWAIC by scale and fragmentation definition
# =============================================================================
# Definition C is the reference (dWAIC = 0) at each scale; bars show how much
# worse Definitions A and B are. The grey band marks dWAIC < 2 (equivalent
# support).

fig_waic <- tibble::tribble(
  ~scale,   ~def, ~dWAIC,
  "7 ha",   "A", 6.2,   "7 ha",   "B", 2.3,   "7 ha",   "C", 0.0,
  "78 ha",  "A", 2.8,   "78 ha",  "B", 2.6,   "78 ha",  "C", 0.0,
  "314 ha", "A", 1.0,   "314 ha", "B", 0.9,   "314 ha", "C", 0.0
) %>%
  mutate(scale  = factor(scale, levels = c("7 ha", "78 ha", "314 ha")),
         plot_h = ifelse(def == "C", 0.15, dWAIC),   # small visible bar for C
         lab    = ifelse(def == "C", "0", sprintf("%.1f", dWAIC)))

fig3 <- ggplot(fig_waic, aes(scale, plot_h, fill = def)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 2,
           fill = "grey88", alpha = 0.6) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72,
           colour = "grey25", linewidth = 0.4) +
  geom_text(aes(label = lab), position = position_dodge(width = 0.8),
            vjust = -0.45, size = 3.4, colour = "grey15") +
  scale_fill_manual(values = def_cols, labels = def_labs,
                    name = "Fragmentation definition") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     breaks = seq(0, 6, 2)) +
  labs(x = "Spatial scale", y = expression(Delta * "WAIC")) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme_classic(base_size = 13) +
  theme(axis.title      = element_text(face = "bold"),
        legend.position = "bottom",
        legend.title    = element_text(face = "bold"),
        plot.margin     = margin(t = 12, r = 22, b = 8, l = 12))

ggsave(file.path(FIG_DIR, "fig3_dWAIC_by_scale.png"), fig3,
       width = 9.5, height = 5.8, dpi = 600, bg = "white")

# =============================================================================
# FIGURE 4 - occupancy coefficient forest plot
# =============================================================================
# Habitat amount and edge density are from the additive model (M3); the
# interaction term is from the interaction model (M4). Hollow points overlap 0;
# the crossed point is the 7 ha Def C interaction, flagged not interpretable by
# the structure-coefficient check (script 05).

forest_data <- tibble::tribble(
  ~definition, ~scale,  ~parameter,     ~mean,   ~lower,  ~upper,
  # Habitat amount - additive (M3)
  "A","150m","prop_habitat",  0.7159,  0.2640,  1.2051,
  "B","150m","prop_habitat",  0.7193,  0.2707,  1.2073,
  "C","150m","prop_habitat",  0.8968,  0.4791,  1.3678,
  "A","500m","prop_habitat",  0.8912,  0.3698,  1.4916,
  "B","500m","prop_habitat",  0.9731,  0.4290,  1.6000,
  "C","500m","prop_habitat",  1.1558,  0.6220,  1.7454,
  "A","1000m","prop_habitat", 1.2466,  0.6643,  1.9495,
  "B","1000m","prop_habitat", 1.2799,  0.7017,  1.9621,
  "C","1000m","prop_habitat", 1.3717,  0.7644,  2.0506,
  # Edge density - additive (M3)
  "A","150m","edge_density", -0.5084, -0.9974, -0.0630,
  "B","150m","edge_density", -0.5753, -1.0244, -0.1716,
  "C","150m","edge_density", -0.6006, -1.0175, -0.2213,
  "A","500m","edge_density", -0.5718, -1.1237, -0.0600,
  "B","500m","edge_density", -0.5478, -1.0733, -0.0730,
  "C","500m","edge_density", -0.8086, -1.3858, -0.2833,
  "A","1000m","edge_density",-0.2979, -0.8641,  0.2453,
  "B","1000m","edge_density",-0.2499, -0.7634,  0.2560,
  "C","1000m","edge_density",-0.4511, -1.0577,  0.1132,
  # Habitat x Edge - interaction model (M4)
  "A","150m","interaction",   0.4694, -0.0667,  1.0304,
  "B","150m","interaction",   0.6111,  0.0510,  1.2173,
  "C","150m","interaction",   0.5465,  0.0215,  1.1111,
  "A","500m","interaction",  -0.1410, -0.6418,  0.3517,
  "B","500m","interaction",  -0.0750, -0.5505,  0.3756,
  "C","500m","interaction",  -0.1030, -0.6436,  0.3948,
  "A","1000m","interaction",  0.0378, -0.4932,  0.5562,
  "B","1000m","interaction",  0.0861, -0.4155,  0.5665,
  "C","1000m","interaction", -0.0786, -0.6349,  0.4560
) %>%
  mutate(
    status = case_when(
      parameter == "interaction" & definition == "C" & scale == "150m" ~ "Not interpretable",
      sign(lower) == sign(upper)                                        ~ "Excludes 0",
      TRUE                                                              ~ "Overlaps 0"),
    status = factor(status, levels = c("Excludes 0", "Overlaps 0", "Not interpretable")),
    param_label = factor(recode(parameter,
                                prop_habitat = "Habitat amount",
                                edge_density = "Edge density",
                                interaction  = "Habitat \u00d7 Edge"),
                         levels = c("Habitat \u00d7 Edge", "Edge density", "Habitat amount")),
    scale_area = factor(scale, levels = c("150m", "500m", "1000m"),
                        labels = c("7 ha", "78 ha", "314 ha")),
    definition = factor(definition, levels = c("A", "B", "C")))

dodge <- position_dodge(width = 0.6)

fig4 <- ggplot(forest_data, aes(mean, param_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = lower, xmax = upper, colour = definition),
                 height = 0.25, linewidth = 0.8, position = dodge) +
  geom_point(aes(fill = definition, shape = status),
             size = 3, stroke = 0.7, position = dodge) +
  scale_colour_manual(values = def_cols, labels = def_labs, name = "Definition") +
  scale_fill_manual(values = def_cols,  labels = def_labs, name = "Definition") +
  scale_shape_manual(values = c("Excludes 0" = 16, "Overlaps 0" = 1,
                                "Not interpretable" = 13),
                     name = "95% CRI / Interpretation") +
  facet_wrap(~scale_area) +
  labs(x = "Occupancy coefficient (\u03b2, 95% CRI)", y = NULL) +
  theme_clean(13) +
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave(file.path(FIG_DIR, "fig4_forest_coefficients.png"), fig4,
       width = 10, height = 4.5, dpi = 600, bg = "white")

# =============================================================================
# FIGURE 5 - predicted-occupancy partial effects (Def C additive model, M3)
# =============================================================================
# Left column: coefficient by scale. Right column: predicted occupancy across
# the observed range of the covariate, holding the other at its mean. Top row =
# habitat amount, bottom row = edge density.

# Posterior means / SDs for the Def C additive model, per scale, used to build
# the prediction curves (b0 = intercept, bh/be = habitat/edge slopes; hm/hsd and
# em/esd = mean and SD used to standardise habitat % and edge density).
m3C <- list(
  "150m"  = list(b0 = -3.2240, b0sd = 0.4745, bh = 0.8968, bhsd = 0.2273,
                 be = -0.6006, besd = 0.2013, hm = 0.6974, hsd = 0.3318, em = 91.55, esd = 73.90),
  "500m"  = list(b0 = -3.1869, b0sd = 0.4750, bh = 1.1558, bhsd = 0.2865,
                 be = -0.8086, besd = 0.2809, hm = 0.5400, hsd = 0.2592, em = 88.12, esd = 42.99),
  "1000m" = list(b0 = -3.2192, b0sd = 0.4809, bh = 1.3717, bhsd = 0.3293,
                 be = -0.4511, besd = 0.2993, hm = 0.4552, hsd = 0.2129, em = 83.18, esd = 35.00))

# Build the predicted-occupancy curve for one term ("hab" or "edge").
make_curve <- function(kind) {
  out <- list()
  for (sc in names(m3C)) {
    p <- m3C[[sc]]
    if (kind == "hab") {
      x  <- seq(0.05, 0.99, length.out = 120)
      xs <- (x - p$hm) / p$hsd
      lp <- p$b0 + p$bh * xs
      se <- sqrt(p$b0sd^2 + xs^2 * p$bhsd^2)
      xx <- x * 100
    } else {
      x  <- seq(0, 250, length.out = 120)
      xs <- (x - p$em) / p$esd
      lp <- p$b0 + p$be * xs
      se <- sqrt(p$b0sd^2 + xs^2 * p$besd^2)
      xx <- x
    }
    out[[sc]] <- data.frame(scale = sc, x = xx, psi = plogis(lp),
                            lo = plogis(lp - 1.96 * se), hi = plogis(lp + 1.96 * se))
  }
  do.call(rbind, out) %>%
    mutate(scale_area = factor(scale_ar[scale], levels = c("7 ha", "78 ha", "314 ha")))
}

# Coefficient values for the left-column panels.
coefC <- tibble::tribble(
  ~kind,  ~scale,   ~mean,    ~lower,   ~upper,
  "hab",  "150m",   0.8968,   0.4791,   1.3678,
  "hab",  "500m",   1.1558,   0.6220,   1.7454,
  "hab",  "1000m",  1.3717,   0.7644,   2.0506,
  "edge", "150m",  -0.6006,  -1.0175,  -0.2213,
  "edge", "500m",  -0.8086,  -1.3858,  -0.2833,
  "edge", "1000m", -0.4511,  -1.0577,   0.1132
) %>%
  mutate(scale_area = factor(scale_ar[scale], levels = c("7 ha", "78 ha", "314 ha")))

coef_panel <- function(kind, ylab, tag)
  ggplot(filter(coefC, kind == !!kind), aes(scale_area, mean, color = scale_area)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, linewidth = 0.9) +
  geom_point(size = 3.4) +
  scale_color_manual(values = scale_cols, name = "Scale") +
  labs(x = "Spatial scale", y = ylab, tag = tag) + theme_clean(12)

occ_panel <- function(df, xlab, ymax, tag)
  ggplot(df, aes(x, psi, color = scale_area, fill = scale_area)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10, color = NA) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = scale_cols, name = "Scale") +
  scale_fill_manual(values = scale_cols, name = "Scale") +
  labs(x = xlab, y = "Predicted occupancy", tag = tag) +
  coord_cartesian(ylim = c(0, ymax)) + theme_clean(12)

no_legend <- function(p) p + theme(legend.position = "none")

Hc  <- coef_panel("hab",  expression("Habitat " * beta), "a")
Ho  <- occ_panel(make_curve("hab"),  "Habitat amount (%)",  0.6, "b")
Ec  <- coef_panel("edge", expression("Edge " * beta),    "c")
Eo  <- occ_panel(make_curve("edge"), "Edge density (m/ha)", 0.5, "d")
leg <- get_legend(Ho)

fig5 <- arrangeGrob(
  arrangeGrob(no_legend(Hc), no_legend(Ho), ncol = 2, widths = c(0.8, 1)),
  arrangeGrob(no_legend(Ec), no_legend(Eo), ncol = 2, widths = c(0.8, 1)),
  leg, ncol = 1, heights = c(1, 1, 0.12))

ggsave(file.path(FIG_DIR, "fig5_habitat_edge_partial_effects.png"), fig5,
       width = 9.5, height = 7.6, dpi = 600, bg = "white")

cat("Saved 3 figures to", FIG_DIR, "\n")

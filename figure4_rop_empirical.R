###############################################################################
# figure4_rop_empirical.R
#
# Figure 4. Empirical illustration: AUROC by predictor configuration,
# CV strategy, and learner (primary outcome, any ROP).
#
# Design: two side-by-side panels (LR | XGBoost). Within each panel,
# three predictor configurations on the x-axis; three point shapes for
# naive 10-fold, subject-level 10-fold, and LOCO CV. LOCO shown with
# DeLong 95% CI error bars; 50-replicate 95% prediction intervals shown
# for naive and subject-level CV. Shared y-axis (0.4 to 1.0).
#
# All values hard-coded from verified analysis output:
#   rop_decomposition_xgb_2026-07-02.txt  (Tier 1, Ridge + XGBoost)
# Source file location:
#   rop_decomposition_xgb_2026-07-02.txt should be in the same directory
#   as this script, or update OUT_DIR below to point to your results folder.
#
# Do not re-derive values from memory. If any number is in doubt,
# re-pull from the source .txt file above.
#
# Author: Joseph L. Hagan, ScD, MSPH
# Section of Neonatology, Baylor College of Medicine
###############################################################################

library(ggplot2)
library(dplyr)
library(tibble)

# =============================================================================
# 1. DATA  (primary outcome, any ROP, Tier 1 full predictor set)
# =============================================================================
# Columns:
#   learner   : "LR" or "XGBoost"
#   config    : predictor configuration (factor, ordered left to right)
#   strategy  : CV strategy (factor, controls shape)
#   auroc     : point estimate
#   lo        : lower bound (2.5th percentile of 50 reps for naive/subject;
#               DeLong 95% CI lower for LOCO)
#   hi        : upper bound (97.5th percentile / DeLong 95% CI upper for LOCO)

fig4_data <- tribble(
  ~learner,  ~config,             ~strategy,    ~auroc,  ~lo,    ~hi,

  # ---- LR (penalized logistic regression) ----
  # Fixed-only
  "LR", "Fixed-only",       "Naive 10-fold",   0.792,  0.792,  0.793,
  "LR", "Fixed-only",       "Subject 10-fold", 0.703,  0.679,  0.727,
  "LR", "Fixed-only",       "LOCO",            0.706,  0.603,  0.810,
  # Longitudinal-only
  "LR", "Longitudinal-only","Naive 10-fold",   0.682,  0.679,  0.684,
  "LR", "Longitudinal-only","Subject 10-fold", 0.593,  0.555,  0.622,
  "LR", "Longitudinal-only","LOCO",            0.606,  0.494,  0.717,
  # Mixed
  "LR", "Mixed",            "Naive 10-fold",   0.818,  0.817,  0.819,
  "LR", "Mixed",            "Subject 10-fold", 0.715,  0.692,  0.740,
  "LR", "Mixed",            "LOCO",            0.718,  0.616,  0.820,

  # ---- XGBoost ----
  # Fixed-only
  "XGBoost", "Fixed-only",       "Naive 10-fold",   1.000,  1.000,  1.000,
  "XGBoost", "Fixed-only",       "Subject 10-fold", 0.644,  0.606,  0.675,
  "XGBoost", "Fixed-only",       "LOCO",            0.653,  0.544,  0.762,
  # Longitudinal-only
  "XGBoost", "Longitudinal-only","Naive 10-fold",   0.742,  0.736,  0.747,
  "XGBoost", "Longitudinal-only","Subject 10-fold", 0.591,  0.547,  0.622,
  "XGBoost", "Longitudinal-only","LOCO",            0.590,  0.477,  0.703,
  # Mixed
  "XGBoost", "Mixed",            "Naive 10-fold",   1.000,  1.000,  1.000,
  "XGBoost", "Mixed",            "Subject 10-fold", 0.642,  0.602,  0.681,
  "XGBoost", "Mixed",            "LOCO",            0.637,  0.527,  0.748
)

# Factor ordering
fig4_data <- fig4_data %>%
  mutate(
    learner  = factor(learner,   levels = c("LR", "XGBoost")),
    config   = factor(config,    levels = c("Fixed-only",
                                            "Longitudinal-only",
                                            "Mixed")),
    strategy = factor(strategy,  levels = c("Naive 10-fold",
                                            "Subject 10-fold",
                                            "LOCO"))
  )

# =============================================================================
# 2. AESTHETICS
# =============================================================================
# Shapes: filled circle (naive), open circle (subject), filled diamond (LOCO)
shape_vals  <- c("Naive 10-fold" = 16, "Subject 10-fold" = 1, "LOCO" = 18)
size_vals   <- c("Naive 10-fold" = 2.5, "Subject 10-fold" = 2.5, "LOCO" = 3.0)

# Horizontal dodge so the three strategy points don't overlap within config
dodge_w <- 0.45

# Panel labels
learner_labels <- c("LR" = "LR", "XGBoost" = "XGBoost")

# =============================================================================
# 3. PLOT
# =============================================================================
fig4 <- ggplot(fig4_data,
               aes(x      = config,
                   y      = auroc,
                   shape  = strategy,
                   size   = strategy,
                   group  = strategy)) +
  # Error bars first (drawn behind points)
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    position = position_dodge(width = dodge_w),
    width    = 0.12,
    linewidth = 0.45,
    color    = "grey40"
  ) +
  # Points on top
  geom_point(
    position = position_dodge(width = dodge_w),
    color    = "black",
    fill     = "black"
  ) +
  # Horizontal reference line at AUROC = 0.5 (chance)
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey60",
             linewidth  = 0.4) +
  # Panel per learner
  facet_wrap(~ learner, nrow = 1, labeller = labeller(learner = learner_labels)) +
  # Scales
  scale_shape_manual(values = shape_vals,  name = "CV strategy") +
  scale_size_manual( values = size_vals,   name = "CV strategy") +
  scale_y_continuous(
    limits = c(0.40, 1.02),
    breaks = seq(0.4, 1.0, by = 0.1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  # Axis labels (no figure title per Joe's convention)
  labs(
    x = "Predictor configuration",
    y = "AUROC"
  ) +
  # Theme: clean, minimal grid, consistent with Figures 2 and 3
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor    = element_blank(),
    panel.grid.major.x  = element_blank(),
    strip.background    = element_rect(fill = "grey92", color = NA),
    strip.text          = element_text(face = "bold", size = 9),
    axis.text.x         = element_text(angle = 20, hjust = 1, size = 8.5),
    axis.title          = element_text(size = 9),
    legend.position     = "bottom",
    legend.title        = element_text(size = 8.5),
    legend.text         = element_text(size = 8),
    legend.key.size     = unit(0.5, "lines")
  )

print(fig4)

# =============================================================================
# 4. EXPORT
# =============================================================================
# Set OUT_DIR to the folder where you want the figure saved.
# "." saves to the current working directory (getwd()).
# Change to an absolute path if preferred, e.g. "C:/my/path/to/results".
OUT_DIR <- "."

ggsave(
  file.path(OUT_DIR, "Hagan_Figure4.tiff"),
  plot        = fig4,
  width       = 17.4,   # cm; fits a two-panel figure at Wiley full-width
  height      = 10.0,   # cm
  units       = "cm",
  dpi         = 300,
  compression = "lzw"
)

ggsave(
  file.path(OUT_DIR, "Hagan_Figure4.eps"),
  plot   = fig4,
  width  = 17.4,
  height = 10.0,
  units  = "cm",
  device = cairo_ps
)

cat("Figure 4 saved to:", OUT_DIR, "\n")

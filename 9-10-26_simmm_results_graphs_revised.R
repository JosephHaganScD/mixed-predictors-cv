# =============================================================================
# All figures: Mixed-Predictor CV Optimism Decomposition Manuscript
# Target journal: Statistics in Medicine
# Author: Joseph L. Hagan, ScD, MSPH
#
# Figures produced:
#   Figure 1  -- figures/fig1_mechanism_schematic.pdf / .tiff
#   Figure 2  -- figures/figure2_icc_bridge.pdf / .tiff
#   Figure 3  -- figures/figure3_composition.tiff
#
# All input CSVs are read from the working directory set below.
# All output figures are written to a "figures" subfolder within that directory.
#
# Revision history (9-10-26):
#   Figure 2: LR recovery percentages and LR dashed reference line removed.
#             Only XGBoost percentages and reference line retained. LR
#             convergence documented numerically in Table 2 per Section 4.2.
#             Caption updated accordingly.
#   Figure 3: INPUT_CSV path corrected; all paths now use bare filenames
#             relative to the single setwd() call below.
#   setwd(): corrected to actual OneDrive file location.
# =============================================================================

# ── Packages ------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)

# ── Working directory ---------------------------------------------------------
# All input CSVs and output figures resolve relative to this path.
# Set working directory to the folder containing the summary CSV files.
# If running from the cloned repository, this line can be removed.
# setwd("path/to/ALL FILES")   # uncomment and update if needed
dir.create("figures", showWarnings = FALSE)


# =============================================================================
# Table 1 Brier column (printed to console; not a figure)
# =============================================================================

sim_results_v8 <- read.csv("sim_results_v8.csv")

table1_brier <- sim_results_v8 %>%
  mutate(brier_delta_naive = brier_true - brier_naive) %>%
  group_by(condition_id, arm, learner) %>%
  summarise(cond_mean = mean(brier_delta_naive, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(arm, learner) %>%
  summarise(
    n_conditions = n(),
    mean_brier   = mean(cond_mean),
    min_brier    = min(cond_mean),
    max_brier    = max(cond_mean),
    .groups = "drop"
  ) %>%
  arrange(arm, learner)

print(table1_brier, width = Inf)

# Confirm all rows are valid
table(sim_results_v8$warning_flags == "" | is.na(sim_results_v8$warning_flags))


# =============================================================================
# Figure 1. Mechanism of identity-mediated leakage under naive partitioning
#
# v3 changes from v2:
#   - Panel A fold assignment is now random per observation (seeded), not
#     aligned by time column, matching how naive record-level CV actually
#     splits. Highlighted subject is guaranteed a train/test mix so the
#     leakage point still reads clearly.
#   - Panel B logic unchanged (5 subjects train, 1 test) -- confirmed correct
#     by sampling pixel colors from the v2 render.
#   - Added S1-S6 (subject) and T1-T4 (time) axis labels to panels A and B.
#   - Added S1-S3 axis labels to panel C.
#   - Legend is now collected once via patchwork and placed under the full
#     A/B stack, instead of embedded only in panel A.
#   - Panel C bottom annotation shortened; interpretation belongs in the
#     figure legend text, not baked into the image.
#   - Fixed legend bug: scale_colour_manual() was pairing "train fold" /
#     "test fold" labels to alphabetically-sorted factor levels, silently
#     swapping legend text. Added breaks = c("train","test") to force correct
#     pairing. Tile colors themselves were never wrong.
#   - Panel B annotation reworded from "all rows of a subject share one fold"
#     to "all repeated observations from a subject are assigned to the same
#     fold" since rows are now explicitly labeled as subjects (S1-S6).
#
# Three-panel conceptual schematic; no data plotted. Referenced in Section 1.
# Layout: (A over B) beside C, via patchwork.
# =============================================================================

# ---- shared constants -------------------------------------------------------
n_subj    <- 6
n_obs     <- 4
train_col <- "#2b6cb0"
test_col  <- "#dd6b20"

subj_symbol <- tibble(subject = 1:n_subj,
                      symbol  = c("/", "\\", "x", ".", "o", "+"))

long_val <- function(s, t) 0.5 + 0.35 * sin(1.3 * (s - 1) + 1.1 * (t - 1))

symbol_offsets <- expand.grid(dx = c(-0.28, 0, 0.28), dy = c(-0.18, 0.18))

build_symbol_layer <- function(tiles_df) {
  do.call(rbind, lapply(seq_len(nrow(tiles_df)), function(i) {
    row <- tiles_df[i, ]
    data.frame(x     = row$obs + symbol_offsets$dx,
               y     = row$plot_y + symbol_offsets$dy,
               label = row$symbol)
  }))
}

# Subject rows read top-to-bottom as S1..S6; plot_y is flipped.
schematic_axes <- function(hl_left_margin   = 25,
                           hl_bottom_margin = 20,
                           hl_right_margin  = 95) {
  list(
    scale_x_continuous(breaks = 1:n_obs,  labels = paste0("T", 1:n_obs)),
    scale_y_continuous(breaks = 1:n_subj, labels = paste0("S", n_subj:1)),
    theme(axis.text.x = element_text(size = 8, margin = margin(t = 2)),
          axis.text.y = element_text(size = 8, margin = margin(r = 2)),
          axis.ticks  = element_blank(),
          plot.margin = margin(5, hl_right_margin,
                               hl_bottom_margin, hl_left_margin))
  )
}

# ---- panel A: naive partitioning --------------------------------------------
set.seed(42)
fold_matrix <- matrix(
  sample(c("train", "test"), n_subj * n_obs, replace = TRUE,
         prob = c(0.75, 0.25)),
  nrow = n_subj, ncol = n_obs)

hl_subject <- 3
if (length(unique(fold_matrix[hl_subject, ])) == 1) {
  fold_matrix[hl_subject, 1] <- "train"
  fold_matrix[hl_subject, 2] <- "test"
}

tiles_A <- expand.grid(subject = 1:n_subj, obs = 1:n_obs)
tiles_A$fold <- fold_matrix[cbind(tiles_A$subject, tiles_A$obs)]
tiles_A <- tiles_A %>%
  mutate(value  = mapply(long_val, subject, obs),
         plot_y = n_subj - subject + 1) %>%
  left_join(subj_symbol, by = "subject")

symbols_A <- build_symbol_layer(tiles_A)
hl_A      <- n_subj - hl_subject + 1

panel_A <- ggplot(tiles_A, aes(x = obs, y = plot_y)) +
  geom_tile(aes(fill = value, colour = fold),
            linewidth = 1.1, width = 0.92, height = 0.85) +
  geom_text(data = symbols_A, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.0, colour = "grey20", alpha = 0.8) +
  scale_fill_gradient(low = "#c6dbef", high = "#2b6cb0", guide = "none") +
  scale_colour_manual(values = c(train = train_col, test = test_col),
                      breaks = c("train", "test"),
                      name = NULL, labels = c("train fold", "test fold")) +
  annotate("rect",
           xmin = 0.5, xmax = n_obs + 0.5,
           ymin = hl_A - 0.45, ymax = hl_A + 0.45,
           fill = NA, colour = "black", linetype = "dashed", linewidth = 1.0) +
  annotate("segment",
           x = n_obs + 0.6, xend = n_obs + 1.1,
           y = hl_A, yend = hl_A + 0.9,
           arrow = arrow(length = unit(0.15, "cm"))) +
  annotate("text",
           x = n_obs + 1.15, y = hl_A + 1.2, hjust = 0, size = 3.0,
           label = "same fingerprint appears\nin both folds\n\u2192 constant outcome recoverable") +
  coord_fixed(clip = "off") +
  labs(title = "A. Naive partitioning (record-level folds)") +
  theme_void() +
  schematic_axes() +
  theme(plot.title = element_text(face = "bold", size = 11,
                                  hjust = 0, margin = margin(b = 8)))

# ---- panel B: subject-level partitioning ------------------------------------
tiles_B <- expand.grid(subject = 1:n_subj, obs = 1:n_obs) %>%
  mutate(value  = mapply(long_val, subject, obs),
         fold   = if_else(subject == 5, "test", "train"),
         plot_y = n_subj - subject + 1) %>%
  left_join(subj_symbol, by = "subject")

symbols_B <- build_symbol_layer(tiles_B)
hl_B      <- n_subj - 5 + 1

panel_B <- ggplot(tiles_B, aes(x = obs, y = plot_y)) +
  geom_tile(aes(fill = value, colour = fold),
            linewidth = 1.1, width = 0.92, height = 0.85) +
  geom_text(data = symbols_B, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.0, colour = "grey20", alpha = 0.8) +
  scale_fill_gradient(low = "#c6dbef", high = "#2b6cb0", guide = "none") +
  scale_colour_manual(values = c(train = train_col, test = test_col),
                      guide = "none") +
  annotate("rect",
           xmin = 0.5, xmax = n_obs + 0.5,
           ymin = hl_B - 0.45, ymax = hl_B + 0.45,
           fill = NA, colour = "black", linetype = "dashed", linewidth = 1.0) +
  annotate("segment",
           x = n_obs + 0.6, xend = n_obs + 1.1,
           y = hl_B, yend = hl_B + 0.9,
           arrow = arrow(length = unit(0.15, "cm"))) +
  annotate("text",
           x = n_obs + 1.15, y = hl_B + 1.2, hjust = 0, size = 3.0,
           label = "all repeated observations from a\nsubject are assigned to the same fold\n\u2192 no shared fingerprint") +
  coord_fixed(clip = "off") +
  labs(title = "B. Subject-level partitioning (subject-level folds)") +
  theme_void() +
  schematic_axes() +
  theme(plot.title = element_text(face = "bold", size = 11,
                                  hjust = 0, margin = margin(b = 8)))

# ---- panel C: ICC continuum -------------------------------------------------
subj_mean_base <- c(-1, 0, 1)
resid_draws <- matrix(c(-0.6,  0.3, -0.9,  0.7,
                         0.8, -0.4,  0.5, -0.6,
                        -0.3,  0.9, -0.7,  0.2),
                      nrow = 3, byrow = TRUE)
subj_col <- c("#2b6cb0", "#2f855a", "#b7791f")

build_icc_panel <- function(icc, title_lab, show_annotation = FALSE) {
  means    <- subj_mean_base * sqrt(icc)
  resid    <- resid_draws    * sqrt(1 - icc)
  df       <- expand.grid(subject = 1:3, obs = 1:4)
  df$value <- mapply(function(s, t) means[s] + resid[s, t],
                     df$subject, df$obs)
  means_df <- tibble(subject = 1:3, mean = means)

  p <- ggplot(df, aes(x = subject)) +
    geom_segment(data = means_df,
                 aes(x = subject - 0.32, xend = subject + 0.32,
                     y = mean, yend = mean, colour = factor(subject)),
                 linewidth = 1.0, linetype = "dashed") +
    geom_jitter(aes(y = value, fill = factor(subject)),
                width = 0.12, height = 0,
                shape = 21, size = 2.6, colour = "black", stroke = 0.5) +
    scale_fill_manual(values   = subj_col, guide = "none") +
    scale_colour_manual(values = subj_col, guide = "none") +
    scale_x_continuous(breaks = 1:3, labels = c("S1", "S2", "S3"),
                       limits = c(0.4, 3.6)) +
    coord_cartesian(ylim = c(-1.8, 1.8), clip = "off") +
    labs(title = title_lab) +
    theme_void() +
    theme(plot.title  = element_text(size = 9, hjust = 0,
                                     margin = margin(b = 4)),
          axis.text.x = element_text(size = 8, margin = margin(t = 2)),
          plot.margin = margin(5, 5, if (show_annotation) 32 else 5, 5))

  if (show_annotation) {
    p <- p + annotate("text", x = 2, y = -2.5, size = 2.7, hjust = 0.5,
                      label = "\u2192 same constant-fingerprint condition as panel B")
  }
  p
}

panel_C1 <- build_icc_panel(0.3, "ICC = 0.3")
panel_C2 <- build_icc_panel(0.7, "ICC = 0.7")
panel_C3 <- build_icc_panel(1.0, "ICC = 1.0", show_annotation = TRUE)
panel_C  <- panel_C1 / panel_C2 / panel_C3

# ---- combine and save -------------------------------------------------------
left_col <- (panel_A / panel_B) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank())

fig1 <- (left_col | panel_C) + plot_layout(widths = c(2, 1))

print(fig1)

ggsave("figures/fig1_mechanism_schematic.pdf",  fig1,
       width = 11, height = 6.5, units = "in")
ggsave("figures/fig1_mechanism_schematic.tiff", fig1,
       width = 11, height = 6.5, units = "in",
       dpi = 300, compression = "lzw")

message("Figure 1 saved.")


# =============================================================================
# Figure 2. The Intraclass Correlation Bridge  (Sections 3.2, 4.2)
#
# 9-10-26 revision (Item 1):
#   - LR recovery percentages removed; LR dashed reference line removed.
#   - Only the XGBoost dashed reference line and XGBoost recovery percentages
#     are shown. LR convergence is documented numerically in Table 2.
#   - ref_lines and recovery_df both filtered to XGBoost.
#   - Caption updated.
#
# Input:  sim_summary_v8.csv
# Output: figures/figure2_icc_bridge.pdf / .tiff
# =============================================================================

LEARNER_LABELS_F2 <- c(ridge_linear = "LR", xgboost = "XGBoost")
LEARNER_COLORS_F2 <- c(LR = "#0072B2", XGBoost = "#D55E00")  # Okabe-Ito

OUTPUT_PDF_F2 <- "figures/figure2_icc_bridge.pdf"
OUTPUT_TIF_F2 <- "figures/figure2_icc_bridge.tiff"

# recovery_df is XGBoost-only, so TRUE is safe.
SHOW_RECOVERY_LABELS <- TRUE

# ── 1. Arm B data (averaged over rho and R2_total) ----------------------------
summ <- read.csv("sim_summary_v8.csv")

armB <- summ %>%
  filter(arm == "B_long") %>%
  mutate(learner = factor(unname(LEARNER_LABELS_F2[learner]),
                          levels = c("LR", "XGBoost"))) %>%
  group_by(ICC, p_L, n, learner) %>%
  summarise(mean_delta_naive = mean(mean_delta_naive), .groups = "drop") %>%
  mutate(
    pL_label = factor(paste(p_L, "predictors"),
                      levels = paste(sort(unique(p_L)), "predictors")),
    n_label  = factor(paste("n =", n),
                      levels = paste("n =", sort(unique(n))))
  )

# ── 2. Arm A reference values (ICC = 1 limit; p_F in {2, 5}) ------------------
armA_ref <- summ %>%
  filter(arm == "A_fixed", p_F %in% c(2, 5)) %>%
  mutate(learner = factor(unname(LEARNER_LABELS_F2[learner]),
                          levels = c("LR", "XGBoost"))) %>%
  group_by(p_F, n, learner) %>%
  summarise(ref_value = mean(mean_delta_naive), .groups = "drop") %>%
  rename(p_L = p_F)

armB <- armB %>%
  left_join(armA_ref, by = c("p_L", "n", "learner")) %>%
  mutate(recovery_pct = round(100 * mean_delta_naive / ref_value))

# ── 3. XGBoost-only reference line data --------------------------------------
# LR's Arm A limit is not drawn; convergence documented in Table 2.
ref_lines <- armB %>%
  distinct(pL_label, n_label, learner, ref_value) %>%
  filter(learner == "XGBoost")

# ── 4. Reference line label (first facet only) --------------------------------
# pL_label/n_label wrapped in factor() to preserve level order across layers;
# plain character would cause ggplot2 to sort facets alphabetically.
label_df <- ref_lines %>%
  filter(pL_label == levels(armB$pL_label)[1],
         n_label  == levels(armB$n_label)[1]) %>%
  summarise(y = max(ref_value) + 0.03) %>%
  mutate(pL_label = factor(levels(armB$pL_label)[1],
                            levels = levels(armB$pL_label)),
         n_label  = factor(levels(armB$n_label)[1],
                            levels = levels(armB$n_label)),
         x = 0.5, label = "fixed covariates\n(ICC = 1)")

# ── 5. XGBoost recovery percentages at ICC = 0.9 -----------------------------
recovery_df <- armB %>%
  filter(ICC == max(ICC), learner == "XGBoost")

# ── 6. Plot -------------------------------------------------------------------
p2 <- ggplot(armB, aes(x = ICC, y = mean_delta_naive, color = learner)) +
  # Single dashed line per facet: XGBoost Arm A limit (orange).
  geom_hline(data = ref_lines,
             aes(yintercept = ref_value, color = learner),
             linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_grid(pL_label ~ n_label) +
  scale_x_continuous(breaks = c(0.3, 0.7, 0.9), limits = c(0.25, 0.95)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  scale_color_manual(values = LEARNER_COLORS_F2, name = NULL) +
  labs(x = "Intraclass correlation",
       y = "Naive cross-validation optimism (AUROC units)",
       caption = paste(
         "Percentages: XGBoost naive optimism at ICC = 0.9 as a percentage",
         "of the fixed-covariate (ICC = 1) reference.",
         "LR convergence to its ICC = 1 limit is documented numerically in Table 2.")) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(legend.position  = "top",
        strip.background = element_rect(fill = "grey95", color = NA),
        plot.margin      = margin(5.5, 12, 5.5, 5.5),
        plot.caption     = element_text(hjust = 0, size = 8, color = "grey30"))

# vjust = 1 anchors the TOP of the text block at y so it grows downward from
# a known point; combined with scale headroom this avoids clipping.
p2 <- p2 +
  geom_text(data = label_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1,
            size = 2.6, lineheight = 0.9, color = "grey30")

if (SHOW_RECOVERY_LABELS) {
  p2 <- p2 +
    geom_text(data = recovery_df,
              aes(label = paste0(recovery_pct, "%")),
              vjust = -0.8, size = 2.6, fontface = "bold",
              show.legend = FALSE)
}

print(p2)

# cairo_pdf: closer to WYSIWYG than base pdf(), which sometimes clips
# in-panel text due to different text-width metrics.
ggsave(OUTPUT_PDF_F2, p2, width = 7.5, height = 5,   units = "in",
       device = cairo_pdf)
ggsave(OUTPUT_TIF_F2, p2, width = 190, height = 127, units = "mm",
       dpi = 300, compression = "lzw", device = "tiff", type = "cairo")

message("Figure 2 saved.")


# =============================================================================
# Figure 3. Leakage Against Fixed-Predictor Composition, by Learner
#                                                         (Sections 3.3, 4.3)
#
# 9-10-26 revision:
#   - INPUT_CSV corrected to bare filename; categorical_results.csv is in
#     the same ALL FILES folder as the other CSVs.
#
# Input:  categorical_results.csv (90 rows)
# Output: figures/figure3_composition.tiff
# =============================================================================

OUTPUT_TIF_F3 <- "figures/figure3_composition.tiff"

LEARNER_LABELS_F3 <- c(logistic = "LR", xgboost = "XGBoost")
N_COLORS_F3 <- c("50" = "#9ECAE1", "100" = "#4292C6", "200" = "#08519C")

COMPOSITION_LEVELS <- c("5 continuous",
                         "3 continuous + 2 binary",
                         "2 continuous + 3 binary",
                         "1 continuous + 4 binary",
                         "5 binary")

# ── 1. Load and prepare ------------------------------------------------------
raw_df <- read.csv("categorical_results.csv")

raw_df <- raw_df %>%
  mutate(
    delta_naive = auroc_naive.auroc - auroc_true.auroc,
    learner     = factor(unname(LEARNER_LABELS_F3[learner]),
                         levels = c("LR", "XGBoost")),
    composition = factor(composition, levels = COMPOSITION_LEVELS)
  )

# Average over R2 (3 levels) within each composition x n x learner cell.
plot_df <- raw_df %>%
  group_by(composition, n, learner) %>%
  summarise(delta_naive = mean(delta_naive),
            u_raw       = mean(u_raw),
            .groups = "drop")

# ── 2. Uniqueness annotations -------------------------------------------------
# u = 1.000 for first four compositions (one label per composition).
# u varies by n only for "5 binary" (three per-point labels to avoid the
# long combined string clipping at the panel edge).
uniqueness_binary <- plot_df %>%
  filter(learner == "LR", composition == "5 binary") %>%
  mutate(label   = sprintf("u = %.3f", u_raw),
         y       = delta_naive + 0.015,
         learner = factor("LR", levels = c("LR", "XGBoost")))

uniqueness_const <- plot_df %>%
  filter(learner == "LR", composition != "5 binary") %>%
  group_by(composition) %>%
  summarise(label = sprintf("u = %.3f", u_raw[1]), .groups = "drop") %>%
  mutate(learner = factor("LR", levels = c("LR", "XGBoost")),
         y = max(uniqueness_binary$y) + 0.03)

# ── 3. Ceiling-percentage annotations (XGBoost only) -------------------------
ceiling_df <- raw_df %>%
  filter(learner == "XGBoost") %>%
  group_by(composition) %>%
  summarise(pct_ceiling = round(100 * mean(auroc_naive.auroc >= 0.999)),
            .groups = "drop") %>%
  mutate(label   = paste0(pct_ceiling, "%"),
         learner = factor("XGBoost", levels = c("LR", "XGBoost")),
         y       = max(plot_df$delta_naive[plot_df$learner == "XGBoost"]) * 1.05)

# ── 4. Caption (pre-wrapped; ggplot silently clips long single-line captions) -
CAPTION_TEXT_F3 <- paste(
  paste(strwrap(
    "u (LR panel): mean raw subject uniqueness for that composition and sample size; a property of the cohort, not the learner.",
    width = 95), collapse = "\n"),
  paste(strwrap(
    "Percentages (XGBoost panel): share of the 9 underlying conditions (3 n x 3 signal-strength levels) reaching the AUROC 0.999 ceiling.",
    width = 95), collapse = "\n"),
  sep = "\n")

# ── 5. Plot -------------------------------------------------------------------
p3 <- ggplot(plot_df, aes(x = composition, y = delta_naive,
                           color = factor(n), group = factor(n))) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_wrap(~ learner, nrow = 1) +
  scale_color_manual(values = N_COLORS_F3, name = "n") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x       = NULL,
       y       = "Naive cross-validation optimism (AUROC units)",
       caption = CAPTION_TEXT_F3) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(legend.position  = "top",
        strip.background = element_rect(fill = "grey95", color = NA),
        axis.text.x      = element_text(angle = 25, hjust = 1),
        plot.margin      = margin(5.5, 12, 5.5, 5.5),
        plot.caption     = element_text(hjust = 0, size = 8, color = "grey30"))

p3 <- p3 +
  geom_text(data = uniqueness_const,
            aes(x = composition, y = y, label = label),
            inherit.aes = FALSE, size = 2.4, color = "grey30") +
  geom_text(data = uniqueness_binary,
            aes(x = composition, y = y, label = label),
            inherit.aes = FALSE, size = 2.2, color = "grey30") +
  geom_text(data = ceiling_df,
            aes(x = composition, y = y, label = label),
            inherit.aes = FALSE, size = 2.6, color = "grey30",
            fontface = "bold")

print(p3)

ggsave(OUTPUT_TIF_F3, p3, width = 180, height = 100, units = "mm",
       dpi = 300, compression = "lzw", device = "tiff", type = "cairo")

message("Figure 3 saved.")
message("All figures complete. Output folder: ",
        normalizePath("figures", mustWork = FALSE))


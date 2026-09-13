###############################################################################
# metrics.R — Performance Metrics and Cohort Uniqueness Statistics
# Mixed-Predictor CV Optimism Decomposition Simulation
#
# Functions:
#   compute_auroc()           — subject-level AUROC via Wilcoxon-Mann-Whitney
#   compute_brier()           — mean squared error of predicted probabilities
#   compute_calibration()     — calibration intercept and slope (logistic)
#   compute_uniqueness_stats() — equivalence-class structure of fixed predictors
#
# Calibration outputs are stored in sim_results but are NOT analyzed in
# this manuscript; they are reserved for the planned standalone calibration
# follow-up paper (spec memo Section 8).
#
# Author: Joseph L. Hagan, ScD, MSPH
###############################################################################

suppressPackageStartupMessages(library(pROC))


# ---------------------------------------------------------------------------
# compute_auroc()
#
# Computes the area under the ROC curve at the subject level.
# Uses pROC::auc() with direction = "<" (higher score -> higher P(Y=1)).
# Returns NA with a warning if there are fewer than 2 unique outcome values.
# ---------------------------------------------------------------------------
compute_auroc <- function(prob_pred, Y) {

  if (length(unique(Y)) < 2L) {
    warning("compute_auroc: fewer than 2 outcome classes; returning NA.")
    return(NA_real_)
  }

  suppressMessages(
    as.numeric(pROC::auc(Y, prob_pred, direction = "<", quiet = TRUE))
  )
}


# ---------------------------------------------------------------------------
# compute_brier()
#
# Mean squared error of predicted probabilities against observed outcomes.
# Brier score = mean((p_hat - Y)^2); lower is better.
# ---------------------------------------------------------------------------
compute_brier <- function(prob_pred, Y) {
  mean((prob_pred - Y)^2)
}


# ---------------------------------------------------------------------------
# compute_calibration()
#
# Estimates calibration intercept and slope via logistic regression of Y on
# the logit of predicted probabilities (the "calibration-in-the-large" and
# "calibration slope" of Van Calster et al.).
#
#   logit(E[Y]) = a + b * logit(p_hat)
#
# Intercept a = 0 and slope b = 1 indicate perfect calibration.
# Returns list(intercept, slope); NA for both if fitting fails.
# ---------------------------------------------------------------------------
compute_calibration <- function(prob_pred, Y) {

  # Guard: clamp predictions away from 0/1 to avoid infinite logits
  p_clamped <- pmax(pmin(prob_pred, 1 - 1e-6), 1e-6)
  logit_p   <- log(p_clamped / (1 - p_clamped))

  result <- tryCatch({
    fit <- glm(Y ~ logit_p, family = binomial(link = "logit"))
    list(
      intercept = as.numeric(coef(fit)[1]),
      slope     = as.numeric(coef(fit)[2])
    )
  }, error = function(e) {
    list(intercept = NA_real_, slope = NA_real_)
  })

  result
}


# ---------------------------------------------------------------------------
# compute_fold_recognition_stats()
#
# DESCRIPTIVE DIAGNOSTIC: Fold-specific subject recognizability (v2,
# patched 2026-07-10).
#
# Motivation: characterize, per fold, the proportion of subjects whose
# predictor fingerprint is unique in the training fold, separately for
# fixed-only, longitudinal-only, and joint fingerprints, and, for Arm C
# (mixed predictors), the overlap structure between the fixed-only and
# longitudinal-only unique subject sets.
#
# CHANGE FROM v1: this diagnostic no longer compares joint recognizability
# against any formula-based reference (neither the simple sum u_F + u_L nor
# the independence-adjusted u_F + u_L - u_F*u_L). Both were found to be
# invalid null references for this statistic: joint fingerprint uniqueness
# is not the union of two independent success events, because combining two
# individually non-discriminating predictor sets can create discrimination
# through complementary disambiguation (e.g. fixed bins A,A,B,B and
# longitudinal bins C,D,C,D give u_F = u_L = 0 but u_joint = 1, since the
# joint bins AC,AD,BC,BD are all distinct). This is a structural property of
# combining equivalence classes, unrelated to whether the two leakage
# channels share an identification mechanism, so no formula-based reference
# is valid here (manuscript Section 2.6). Instead, this function now
# additionally returns the OBSERVED overlap, union, and newly-unique-only-
# jointly subject proportions, computed directly from the per-subject
# unique/non-unique indicators rather than inferred from u_F and u_L alone.
# These are purely descriptive; the mechanistic test is the subject-level
# leakage-credit diagnostic (compute_leakage_credits(), below), restricted
# to the paired additivity sub-study.
#
# Design (fold-specific / exact):
#   For each naive CV fold k, training subjects are all n subjects (since
#   under naive fold assignment, every subject has rows in both folds). For
#   each training subject i, the fold-specific X̄_L^k is the mean of
#   longitudinal predictor values over i's training-fold rows only (~(K-1)/K
#   of their T observations). This differs from the full-cohort X̄_L by
#   within-subject noise (SD ≈ σ_e / sqrt(T_train)) and is therefore
#   slightly less unique than the full mean — the approximation (full-cohort
#   X̄_L) would overstate long-predictor uniqueness. For T_const = 30 and
#   K = 10, each training fold retains ~27 observations per subject;
#   per-subject noise in the partial mean is small but non-negligible at
#   high ρ, making the exact version preferable.
#
# Arguments:
#   data_long     : long-format data frame
#   folds_naive   : naive fold vector (from make_naive_folds)
#   predictor_set : "fixed_only", "longitudinal_only", or "mixed"
#   p_F, p_L      : predictor counts
#   K             : number of folds (default 10)
#
# Returns a named list (NA where not applicable to the arm):
#   fold_recog_fixed_q4/q10  — mean across K folds of proportion of subjects
#     with unique X_F fingerprint in the training fold (Arm A and C only)
#   fold_recog_long_q4/q10   — same for fold-specific X̄_L fingerprint
#     (Arm B and C only)
#   fold_recog_joint_q4/q10  — same for joint (X_F + X̄_L) fingerprint
#     (Arm C only)
#   fold_recog_overlap_q4/q10 — mean proportion of subjects uniquely
#     identified under BOTH the fixed-only and longitudinal-only
#     fingerprints (Arm C only)
#   fold_recog_union_q4/q10   — mean proportion of subjects uniquely
#     identified under EITHER the fixed-only or longitudinal-only
#     fingerprint (Arm C only)
#   fold_recog_newly_joint_q4/q10 — mean proportion of subjects that become
#     uniquely identified only once fixed and longitudinal fingerprints are
#     combined, i.e. unique jointly but not unique under either component
#     alone (Arm C only); reflects complementary disambiguation
# ---------------------------------------------------------------------------
compute_fold_recognition_stats <- function(data_long, folds_naive,
                                            predictor_set, p_F, p_L,
                                            K = 10L) {

  has_fixed <- predictor_set %in% c("fixed_only", "mixed")
  has_long  <- predictor_set %in% c("longitudinal_only", "mixed")
  has_joint <- predictor_set == "mixed"

  cols_F <- if (p_F > 0L) paste0("X_F_", seq_len(p_F)) else character(0)
  cols_L <- if (p_L > 0L) paste0("X_L_", seq_len(p_L)) else character(0)

  # Returns a logical vector, one entry per subject (in fp's order), TRUE
  # where that subject's fingerprint is a singleton (unique) in this fold.
  .unique_indicator <- function(fp) {
    sz <- table(fp)
    as.integer(sz[fp]) == 1L
  }

  # Bin a data.frame of continuous columns into n_bins quantile groups.
  # Returns a data.frame of integer bin assignments, one column per input.
  .bin_df <- function(df, n_bins) {
    as.data.frame(lapply(df, function(col) {
      cuts <- quantile(col, probs = seq(0, 1, length.out = n_bins + 1L),
                       type = 7)
      cuts[1L]            <- -Inf
      cuts[length(cuts)]  <-  Inf
      as.integer(cut(col, breaks = cuts, labels = FALSE,
                     include.lowest = TRUE))
    }))
  }

  # Accumulators: one value per fold
  rf_q4  <- numeric(K); rf_q10 <- numeric(K)
  rl_q4  <- numeric(K); rl_q10 <- numeric(K)
  rj_q4  <- numeric(K); rj_q10 <- numeric(K)
  ov_q4  <- numeric(K); ov_q10 <- numeric(K)
  un_q4  <- numeric(K); un_q10 <- numeric(K)
  nj_q4  <- numeric(K); nj_q10 <- numeric(K)

  for (k in seq_len(K)) {

    train_rows <- folds_naive != k
    train_long <- data_long[train_rows, , drop = FALSE]

    # One row per subject (fixed covariates are constant within subject)
    subj_rows <- !duplicated(train_long$subject_id)
    XF_subj   <- if (has_fixed && length(cols_F) > 0L)
                   train_long[subj_rows, cols_F, drop = FALSE]
                 else NULL

    # Fold-specific X̄_L: mean over each subject's TRAINING-FOLD rows only.
    # Using tapply (single pass) rather than looping over subjects.
    if (has_long && length(cols_L) > 0L) {
      Lbar_list <- lapply(cols_L, function(v)
        tapply(train_long[[v]], train_long$subject_id, mean))
      # Ensure subject order matches XF_subj (alphabetical from tapply)
      subj_order <- as.character(train_long$subject_id[subj_rows])
      XL_subj <- as.data.frame(lapply(Lbar_list, function(v) v[subj_order]))
      names(XL_subj) <- cols_L
    } else {
      XL_subj <- NULL
    }

    # Compute fingerprints and unique-indicators at Q4 and Q10. uF4/uF10 and
    # uL4/uL10 (subject-order-aligned logical vectors) are retained, not just
    # their means, so overlap/union/newly-unique can be computed below.
    uF4 <- uF10 <- uL4 <- uL10 <- NULL

    if (has_fixed && !is.null(XF_subj)) {
      bf4  <- .bin_df(XF_subj, 4L);  fp <- do.call(paste, c(bf4,  sep="_"))
      bf10 <- .bin_df(XF_subj, 10L); fp10 <- do.call(paste, c(bf10, sep="_"))
      uF4  <- .unique_indicator(fp)
      uF10 <- .unique_indicator(fp10)
      rf_q4[k]  <- mean(uF4)
      rf_q10[k] <- mean(uF10)
    }

    if (has_long && !is.null(XL_subj)) {
      bl4  <- .bin_df(XL_subj, 4L);  fp <- do.call(paste, c(bl4,  sep="_"))
      bl10 <- .bin_df(XL_subj, 10L); fp10 <- do.call(paste, c(bl10, sep="_"))
      uL4  <- .unique_indicator(fp)
      uL10 <- .unique_indicator(fp10)
      rl_q4[k]  <- mean(uL4)
      rl_q10[k] <- mean(uL10)
    }

    if (has_joint && !is.null(XF_subj) && !is.null(XL_subj)) {
      # Re-bin jointly (not reusing above to avoid bind ordering issues)
      bf4  <- .bin_df(XF_subj, 4L);  bl4  <- .bin_df(XL_subj, 4L)
      bf10 <- .bin_df(XF_subj, 10L); bl10 <- .bin_df(XL_subj, 10L)
      fp4  <- do.call(paste, c(cbind(bf4,  bl4),  sep = "_"))
      fp10 <- do.call(paste, c(cbind(bf10, bl10), sep = "_"))
      uJ4  <- .unique_indicator(fp4)
      uJ10 <- .unique_indicator(fp10)
      rj_q4[k]  <- mean(uJ4)
      rj_q10[k] <- mean(uJ10)

      # Observed overlap/union/newly-unique-only-jointly, computed directly
      # from the per-subject indicators (uF4/uL4/uJ4 are all in the same
      # subject order within this fold), not inferred from a formula.
      ov_q4[k] <- mean(uF4 & uL4)
      un_q4[k] <- mean(uF4 | uL4)
      nj_q4[k] <- mean(uJ4 & !(uF4 | uL4))

      ov_q10[k] <- mean(uF10 & uL10)
      un_q10[k] <- mean(uF10 | uL10)
      nj_q10[k] <- mean(uJ10 & !(uF10 | uL10))
    }
  }

  list(
    fold_recog_fixed_q4       = if (has_fixed) mean(rf_q4)  else NA_real_,
    fold_recog_fixed_q10      = if (has_fixed) mean(rf_q10) else NA_real_,
    fold_recog_long_q4        = if (has_long)  mean(rl_q4)  else NA_real_,
    fold_recog_long_q10       = if (has_long)  mean(rl_q10) else NA_real_,
    fold_recog_joint_q4       = if (has_joint) mean(rj_q4)  else NA_real_,
    fold_recog_joint_q10      = if (has_joint) mean(rj_q10) else NA_real_,
    fold_recog_overlap_q4     = if (has_joint) mean(ov_q4)  else NA_real_,
    fold_recog_overlap_q10    = if (has_joint) mean(ov_q10) else NA_real_,
    fold_recog_union_q4       = if (has_joint) mean(un_q4)  else NA_real_,
    fold_recog_union_q10      = if (has_joint) mean(un_q10) else NA_real_,
    fold_recog_newly_joint_q4  = if (has_joint) mean(nj_q4)  else NA_real_,
    fold_recog_newly_joint_q10 = if (has_joint) mean(nj_q10) else NA_real_
  )
}
#
# Computes the equivalence-class structure of the fixed predictors under
# two levels of discretization: 4-quantile (quartile) and 10-quantile
# (decile) binning of each continuous fixed predictor.
#
# Arguments:
#   data_wide — wide-format dataset (must contain X_F_* columns)
#   p_F       — number of fixed predictors
#
# Returns: list with
#   pct_unique_q4   — proportion of subjects uniquely identified under Q4 bins
#   pct_unique_q10  — proportion of subjects uniquely identified under Q10 bins
#   mean_class_size — mean equivalence-class size under Q4 bins
#   max_class_size  — maximum equivalence-class size under Q4 bins
#
# If p_F == 0, all values are returned as NA (no fixed predictors to analyze).
# ---------------------------------------------------------------------------
compute_uniqueness_stats <- function(data_wide, p_F) {

  if (p_F == 0L) {
    return(list(
      pct_unique_q4   = NA_real_,
      pct_unique_q10  = NA_real_,
      mean_class_size = NA_real_,
      max_class_size  = NA_real_
    ))
  }

  xF_cols <- paste0("X_F_", seq_len(p_F))
  X_F     <- data_wide[, xF_cols, drop = FALSE]
  n       <- nrow(X_F)

  .bin_and_count <- function(X, n_bins) {
    # Discretize each column into n_bins quantile groups
    X_binned <- as.data.frame(lapply(X, function(col) {
      cuts <- quantile(col, probs = seq(0, 1, length.out = n_bins + 1),
                       type = 7)
      cuts[1]              <- -Inf
      cuts[length(cuts)]   <-  Inf
      as.integer(cut(col, breaks = cuts, labels = FALSE, include.lowest = TRUE))
    }))

    # Equivalence class: unique combination of bin values
    class_key   <- do.call(paste, c(X_binned, sep = "_"))
    class_sizes <- table(class_key)

    # pct_unique must be the proportion of SUBJECTS in a singleton class,
    # not the proportion of distinct classes that happen to be singletons.
    # table(class_key) has one entry per distinct class; mean(class_sizes==1)
    # therefore averages over classes, which over- or under-states the
    # subject-level proportion whenever class sizes are unequal. Looking up
    # each subject's own class size first and averaging over subjects fixes
    # this. mean_class_size is unaffected by this distinction: mean(table(.))
    # already equals n / n_classes, the standard k-anonymity definition,
    # regardless of weighting.
    subj_class_size <- as.integer(class_sizes[class_key])

    list(
      pct_unique      = mean(subj_class_size == 1L),
      mean_class_size = mean(class_sizes),
      max_class_size  = max(class_sizes)
    )
  }

  stats_q4  <- .bin_and_count(X_F, 4L)
  stats_q10 <- .bin_and_count(X_F, 10L)

  list(
    pct_unique_q4   = stats_q4$pct_unique,
    pct_unique_q10  = stats_q10$pct_unique,
    mean_class_size = stats_q4$mean_class_size,
    max_class_size  = stats_q4$max_class_size
  )
}


# ---------------------------------------------------------------------------
# compute_leakage_credits()
#
# NEW 2026-07-10: subject-level leakage-credit diagnostic, restricted to the
# paired additivity sub-study. Replaces a nearest-centroid subject-
# identification classifier considered during design, which was rejected as
# tautological: fixed predictors have ICC = 1 by construction (Section 2.1
# of the manuscript), so a held-out subject's fixed-predictor vector is
# identical to that subject's own training-fold centroid, meaning
# nearest-centroid identification using fixed predictors alone would
# trivially achieve ~100% accuracy by construction of the data-generating
# model, not as an empirical finding, and would not test whether the FITTED
# MODELS actually produce leakage for the same subjects.
#
# For each subject i, the leakage credit is
#   C_i = l(Y_i, p_hat_subject_i) - l(Y_i, p_hat_naive_i)
# where l(y,p) is log-loss, l(y,p) = -[y*log(p) + (1-y)*log(1-p)]. A
# POSITIVE credit means the subject's naive-CV prediction had LOWER
# (better) log-loss than their subject-level-CV prediction, i.e. that
# subject benefited from leakage under naive partitioning. This is
# computed separately for the fixed-only model and the longitudinal-only
# model (both fit within the same paired-sub-study iteration, on the same
# generated dataset), then compared for correlation and subject-set
# overlap via compare_leakage_credits(), below.
#
# Arguments:
#   Y            — subject-level binary outcome vector
#   prob_naive   — subject-level predicted probability under naive CV
#   prob_subject — subject-level predicted probability under subject-level CV
#   subject_id   — subject identifiers, aligned with Y/prob_naive/prob_subject
#   eps          — clamp bound to avoid infinite log-loss at p in {0,1}
#
# Returns: named numeric vector of credits, one per subject, names = subject_id
# ---------------------------------------------------------------------------
compute_leakage_credits <- function(Y, prob_naive, prob_subject, subject_id,
                                     eps = 1e-6) {

  p_naive_c   <- pmax(pmin(prob_naive,   1 - eps), eps)
  p_subject_c <- pmax(pmin(prob_subject, 1 - eps), eps)

  .logloss <- function(y, p) -(y * log(p) + (1 - y) * log(1 - p))

  loss_naive   <- .logloss(Y, p_naive_c)
  loss_subject <- .logloss(Y, p_subject_c)

  credits <- loss_subject - loss_naive
  names(credits) <- as.character(subject_id)
  credits
}


# ---------------------------------------------------------------------------
# compare_leakage_credits()
#
# Compares leakage credits (compute_leakage_credits(), above) between two
# models fit to the SAME subjects (e.g. fixed-only vs. longitudinal-only,
# within one paired-sub-study iteration). Tests whether the same subjects
# benefit from leakage under both models, the direct evidence for the
# single-mechanism account referenced in manuscript Section 2.9.
#
# "High credit" subjects are defined as credit > 0 (any degree of leakage
# benefit under naive CV for that model), a natural zero threshold rather
# than an arbitrary cutoff, since 0 is exactly the boundary between naive
# CV being more vs. less accurate than subject-level CV for that subject.
#
# Arguments:
#   credits_a, credits_b — named numeric vectors from compute_leakage_credits(),
#                          for two different models (e.g. fixed-only,
#                          longitudinal-only)
#
# Returns a named list:
#   n_subjects        — number of subjects common to both credit vectors
#   credit_correlation — Pearson correlation between credits_a and credits_b
#                        across common subjects (NA if either has zero
#                        variance or fewer than 3 common subjects)
#   n_high_a, n_high_b — count of subjects with credit > 0 under each model
#   n_high_overlap     — count of subjects with credit > 0 under BOTH models
#   n_high_union       — count of subjects with credit > 0 under EITHER model
#   jaccard_high       — n_high_overlap / n_high_union (NA if union is empty)
# ---------------------------------------------------------------------------
compare_leakage_credits <- function(credits_a, credits_b) {

  common <- intersect(names(credits_a), names(credits_b))
  ca <- credits_a[common]
  cb <- credits_b[common]

  credit_correlation <- if (length(common) >= 3L && sd(ca) > 0 && sd(cb) > 0) {
    cor(ca, cb)
  } else {
    NA_real_
  }

  high_a <- ca > 0
  high_b <- cb > 0
  n_high_overlap <- sum(high_a & high_b)
  n_high_union   <- sum(high_a | high_b)

  list(
    n_subjects         = length(common),
    credit_correlation = credit_correlation,
    n_high_a           = sum(high_a),
    n_high_b           = sum(high_b),
    n_high_overlap     = n_high_overlap,
    n_high_union       = n_high_union,
    jaccard_high       = if (n_high_union > 0L) n_high_overlap / n_high_union
                          else NA_real_
  )
}


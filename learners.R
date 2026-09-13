###############################################################################
# learners.R — Learner Flexibility Levels (v5)
# Mixed-Predictor CV Optimism Decomposition Simulation
#
# Two learner types:
#   "ridge_linear" — ridge logistic regression, linear main effects
#   "xgboost"      — gradient boosted trees
#
# Changes from v4 (patched 2026-07-10):
#   - Hyperparameters (ridge lambda; XGBoost min_child_weight) are no longer
#     computed inside .fit_ridge_linear()/.fit_xgboost() from whatever
#     data_long they are handed. Previously, since these functions were
#     called separately for the full-model fit and for each cross-validation
#     fold (cv_strategies.R), a fold's own subject count (~90% of the full
#     condition's n for subject-level folds) silently changed the
#     regularization strength between the full model and subject-level
#     folds, partially confounding the fold-partitioning-strategy comparison
#     with a regularization difference. Hyperparameters are now computed
#     once via compute_hyperparams(), from the full condition's subject
#     count, and passed explicitly into fit_learner() for the full model and
#     every fold of both cross-validation strategies (see cv_strategies.R).
#   - nthread = 1 added to XGBoost params, so multisession parallelization
#     (run_simulation.R) does not oversubscribe cores via XGBoost's own
#     internal threading.
#
# Changes from v3 (retained for history):
#   - BOTH learners now fit on row-level (data_long) data and predict at the
#     ROW level. v3 fit XGBoost on subject-level wide-aggregated data
#     unconditionally (even under naive/row-level CV), so XGBoost never saw
#     the row-duplicated fixed-covariate fingerprint that ridge experienced
#     during naive CV. That made the ridge-vs-XGBoost optimism comparison
#     confounded by training-data granularity in addition to learner
#     flexibility. Both learners now receive identical training data under
#     identical fold strategies; only the learner differs.
#   - predict_prob() now returns ONE PREDICTION PER ROW (not per subject).
#     Subject-level aggregation (mean of row-level predicted probabilities)
#     is handled centrally in cv_strategies.R via predict_subject_level(),
#     applied identically to both learners.
#
# Author: Joseph L. Hagan, ScD, MSPH
# Version: 5 (2026-07-10)
###############################################################################

suppressPackageStartupMessages({
  library(glmnet)
  library(xgboost)
  library(splines)
})


# ---------------------------------------------------------------------------
# compute_hyperparams()
#
# Computes each learner's hyperparameter(s) once from the full condition's
# subject count. Called once per condition (cv_strategies.R), and the
# returned value is passed unchanged into every subsequent fit_learner()
# call for that condition, full model and every fold of both CV strategies.
# ---------------------------------------------------------------------------
compute_hyperparams <- function(learner_type, n_subj_full) {
  switch(learner_type,
    "ridge_linear" = list(lambda = log(1 + 1 / n_subj_full) / n_subj_full),
    "xgboost"      = list(min_child_weight = max(1L, floor(n_subj_full / 20L))),
    stop("Unknown learner_type: ", learner_type)
  )
}


# ---------------------------------------------------------------------------
# fit_learner() — main entry point
# hyperparams is the list returned by compute_hyperparams(), computed once
# per condition from the full subject count and passed through unchanged
# regardless of which fold (or the full model) is being fit.
# ---------------------------------------------------------------------------
fit_learner <- function(data_long, learner_type, predictor_set, p_F, p_L,
                         hyperparams) {

  switch(learner_type,
    "ridge_linear" = .fit_ridge_linear(data_long, predictor_set, p_F, p_L,
                                        lambda = hyperparams$lambda),
    "xgboost"      = .fit_xgboost(data_long, predictor_set, p_F, p_L,
                                   min_child_weight = hyperparams$min_child_weight),
    stop("Unknown learner_type: ", learner_type)
  )
}


# ---------------------------------------------------------------------------
# Column name helpers
# Row-level predictors only: X_F_* (replicated per row) and X_L_* (raw
# per-timepoint values). No wide/X_L_bar helpers needed; both training and
# prediction now operate on the same row-level representation.
# ---------------------------------------------------------------------------
.cols_F <- function(p_F) if (p_F > 0) paste0("X_F_", seq_len(p_F)) else character(0)
.cols_L <- function(p_L) if (p_L > 0) paste0("X_L_", seq_len(p_L)) else character(0)

.use_cols <- function(predictor_set, p_F, p_L) {
  switch(predictor_set,
    "fixed_only"        = .cols_F(p_F),
    "longitudinal_only" = .cols_L(p_L),
    "mixed"             = c(.cols_F(p_F), .cols_L(p_L)),
    stop("Unknown predictor_set: ", predictor_set)
  )
}

.build_X <- function(data_long, predictor_set, p_F, p_L)
  as.matrix(data_long[, .use_cols(predictor_set, p_F, p_L), drop = FALSE])


# ---------------------------------------------------------------------------
# .fit_ridge_linear()
#
# Ridge logistic regression with a pre-specified fixed lambda, fit on
# row-level data. lambda is now supplied by the caller (computed once per
# condition via compute_hyperparams(), from the full condition's subject
# count), not recomputed here from whichever data_long this call happens to
# receive.
# ---------------------------------------------------------------------------
.fit_ridge_linear <- function(data_long, predictor_set, p_F, p_L, lambda) {

  X <- .build_X(data_long, predictor_set, p_F, p_L)
  y <- data_long$Y

  fit_full <- glmnet(X, y, family = "binomial", alpha = 0, lambda = lambda)

  predict_prob <- function(newdata_long) {
    Xnew <- .build_X(newdata_long, predictor_set, p_F, p_L)
    as.numeric(predict(fit_full, newx = Xnew, s = lambda, type = "response"))
  }

  list(model = fit_full, lambda_used = lambda, predict_prob = predict_prob)
}


# ---------------------------------------------------------------------------
# .fit_xgboost()
#
# Gradient boosted trees, fit on row-level data (one row per observation, Y
# replicated within subject), matching ridge's training granularity.
# min_child_weight is now supplied by the caller (computed once per
# condition via compute_hyperparams(), from the full condition's subject
# count), not recomputed here from whichever data_long this call happens to
# receive; as an integer-valued hyperparameter, recomputing it per fold
# previously produced discrete jumps between the full model and
# subject-level folds (~90% of subjects), not just a small continuous
# shift as with ridge's lambda.
#
# nthread = 1: prevents XGBoost's internal threading from oversubscribing
# cores when many conditions are fit concurrently under multisession
# parallelization (run_simulation.R).
#
# Other hyperparameters unchanged from v3:
#   max_depth = 4, eta = 0.1, nrounds = 50 (fixed), subsample = 0.8,
#   colsample_bytree = 0.8
# ---------------------------------------------------------------------------
.fit_xgboost <- function(data_long, predictor_set, p_F, p_L, min_child_weight) {

  X <- .build_X(data_long, predictor_set, p_F, p_L)
  y <- data_long$Y

  dfull <- xgb.DMatrix(data = X, label = y)

  params <- list(
    objective        = "binary:logistic",
    eval_metric      = "logloss",
    max_depth        = 4L,
    eta              = 0.1,
    subsample        = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = min_child_weight,
    nthread          = 1L
  )

  fit_full <- xgb.train(
    params  = params,
    data    = dfull,
    nrounds = 50L,
    verbose = 0
  )

  predict_prob <- function(newdata_long) {
    Xnew <- .build_X(newdata_long, predictor_set, p_F, p_L)
    as.numeric(predict(fit_full, xgb.DMatrix(data = Xnew)))
  }

  list(model = fit_full, nrounds_used = 50L, predict_prob = predict_prob)
}


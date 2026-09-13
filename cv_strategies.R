###############################################################################
# cv_strategies.R — CV Fold Assignment and Evaluation (v4)
# Mixed-Predictor CV Optimism Decomposition Simulation
#
# Changes from v3 (patched 2026-07-10):
#   - run_cv() now takes a `hyperparams` argument (from
#     learners.R::compute_hyperparams(), computed once per condition from
#     the full subject count) and passes it unchanged into every
#     fit_learner() call inside the fold loop, rather than letting
#     fit_learner() recompute hyperparameters from each fold's own,
#     smaller subject count. See learners.R v5 changelog for the rationale.
#
# Changes from v2 (retained for history):
#   - run_cv() rewritten around row-level prediction + subject-level
#     aggregation (predict_subject_level()), used identically for both fold
#     strategies and both learners now that fit_learner() takes data_long
#     for everyone (see learners.R v4). The data_wide argument and the
#     .long_to_wide_means() pre-aggregation helper are removed: there is no
#     longer a point in the pipeline where predictors are aggregated BEFORE
#     model fitting or prediction, only after, which is what makes the
#     ridge/XGBoost comparison and the naive/subject-level CV comparison
#     consistent across learners.
#   - No functional changes to make_naive_folds(), make_subject_folds(), or
#     make_loco_folds().
#
# Author: Joseph L. Hagan, ScD, MSPH
# Version: 4 (2026-07-10)
###############################################################################


# ---------------------------------------------------------------------------
# make_naive_folds()
# Fold assignment at the observation (row) level — naive strategy.
# ---------------------------------------------------------------------------
make_naive_folds <- function(data_long, K = 10, seed) {
  set.seed(seed)
  n_rows   <- nrow(data_long)
  fold_vec <- integer(n_rows)
  fold_vec[sample(n_rows)] <- rep(seq_len(K), length.out = n_rows)
  fold_vec
}


# ---------------------------------------------------------------------------
# make_subject_folds()
# Fold assignment at the subject level — cluster-aware strategy.
# ---------------------------------------------------------------------------
make_subject_folds <- function(data_long, K = 10, seed) {
  set.seed(seed)
  subjects   <- unique(data_long$subject_id)
  n_subjects <- length(subjects)
  subj_fold  <- integer(n_subjects)
  subj_fold[sample(n_subjects)] <- rep(seq_len(K), length.out = n_subjects)
  names(subj_fold) <- subjects
  subj_fold[as.character(data_long$subject_id)]
}


# ---------------------------------------------------------------------------
# make_loco_folds()
# EMPIRICAL ILLUSTRATION ONLY — not used in simulation factorial.
# Leave-one-cluster-out: one fold per subject.
# ---------------------------------------------------------------------------
make_loco_folds <- function(data_long) {
  data_long$subject_id
}


# ---------------------------------------------------------------------------
# predict_subject_level()  [internal helper]
#
# Generates row-level predicted probabilities from a fitted learner and
# averages them to subject level. Used identically for naive CV, subject-
# level CV, and the full-cohort/test-cohort evaluation, so the aggregation
# method is held fixed across every comparison in the manuscript.
#
# Aggregation uses tapply(), not vapply() over subject IDs with a full-vector
# logical scan per subject. The latter is O(n_subjects x n_rows): for the
# n_test=5000, T_const=30 test cohort (150,000 rows), that is 5000 separate
# 150,000-element scans per call, measured at ~2.2 sec/call versus ~0.03
# sec/call for the tapply form (verified identical output). This call fires
# twice per iteration (once per learner) for every iteration of every
# condition, so the original form alone was on the order of a week of
# wall-clock time across the full factorial. tapply()'s single grouped pass
# avoids it entirely with no change to the computed values.
# ---------------------------------------------------------------------------
.predict_subject_level <- function(fit, newdata_long) {
  row_pred  <- fit$predict_prob(newdata_long)
  prob_subj <- tapply(row_pred, newdata_long$subject_id, mean)
  Y_subj    <- tapply(newdata_long$Y, newdata_long$subject_id, function(y) y[1L])
  subj_ids  <- names(prob_subj)
  data.frame(subject_id = subj_ids,
            Y = as.numeric(Y_subj[subj_ids]),
            prob_pred = as.numeric(prob_subj[subj_ids]),
            stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------------
# run_cv()
#
# Orchestrates cross-validation for a single learner and fold assignment.
# Returns subject-level predicted probabilities, aggregated from row-level
# predictions identically regardless of fold_level or learner_type.
#
# hyperparams: computed once per condition via learners.R::compute_hyperparams()
#              from the full condition's subject count, and passed unchanged
#              into every fold's fit_learner() call (see learners.R v5).
# fold_level = "row"     -> naive CV (folds assigned at observation level)
# fold_level = "subject" -> subject-level CV or LOCO
# ---------------------------------------------------------------------------
run_cv <- function(data_long,
                   learner_type, predictor_set,
                   p_F, p_L,
                   fold_vector,
                   hyperparams,
                   fold_level = c("row", "subject")) {

  fold_level <- match.arg(fold_level)
  folds      <- sort(unique(fold_vector))
  row_preds  <- numeric(nrow(data_long))

  for (k in folds) {

    if (fold_level == "row") {
      train_rows <- fold_vector != k
      test_rows  <- fold_vector == k

      train_long <- data_long[train_rows, , drop = FALSE]
      test_long  <- data_long[test_rows,  , drop = FALSE]

    } else {
      test_subjects  <- unique(data_long$subject_id[fold_vector == k])
      train_subjects <- setdiff(unique(data_long$subject_id), test_subjects)

      train_long <- data_long[data_long$subject_id %in% train_subjects, , drop = FALSE]
      test_long  <- data_long[data_long$subject_id %in% test_subjects,  , drop = FALSE]
      test_rows  <- data_long$subject_id %in% test_subjects
    }

    fit <- fit_learner(train_long, learner_type, predictor_set, p_F, p_L,
                        hyperparams = hyperparams)
    row_preds[test_rows] <- fit$predict_prob(test_long)
  }

  # Aggregate row-level predictions (every row visited exactly once as test
  # across the K folds) to subject level. See .predict_subject_level() above
  # for why tapply() is used instead of vapply()+full-vector-scan.
  prob_subj <- tapply(row_preds, data_long$subject_id, mean)
  Y_subj    <- tapply(data_long$Y, data_long$subject_id, function(y) y[1L])
  subj_ids  <- names(prob_subj)

  data.frame(subject_id = subj_ids,
            Y = as.numeric(Y_subj[subj_ids]),
            prob_pred = as.numeric(prob_subj[subj_ids]),
            stringsAsFactors = FALSE)
}

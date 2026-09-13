###############################################################################
# simulate.R — Single-Iteration Orchestrator (v4)
# Mixed-Predictor CV Optimism Decomposition Simulation
#
# Changes from v3 (patched 2026-07-10):
#   - target_prev now read from condition_row$target_prevalence instead of
#     hardcoded (was 0.62, incorrect; corrected value 0.48 for the main
#     factorial, and variable per-condition for the supplemental grids).
#   - n_subj_full computed once per iteration (before the per-learner loop)
#     and used to compute each learner's hyperparameters via
#     learners.R::compute_hyperparams(), then passed unchanged into the
#     full-model fit and every fold of both run_cv() calls. Previously,
#     fit_learner() recomputed ridge's lambda and XGBoost's min_child_weight
#     internally from whichever data it was handed, so subject-level folds
#     (~90% of n_subj_full) received systematically different regularization
#     than the full model and naive folds, partially confounding the
#     fold-partitioning-strategy comparison with a regularization
#     difference. See learners.R v5 and cv_strategies.R v4 changelogs.
#
# Changes from v2 (retained for history):
#   - iter_seed redesigned to eliminate seed collisions. v2 used
#     master_seed + condition_id*10000 + iter_id, with internal sub-step
#     offsets (+0..+4) that overlapped the iter_id spacing of 1, so e.g. the
#     fold-assignment seed of iteration i was literally identical to the
#     beta-generation seed of iteration i+4 within the same condition
#     (verified: ~96% of seed-uses collided with some other sub-step in
#     steady state). Iteration spacing is now 100 and condition spacing is
#     1e6, leaving headroom far in excess of the handful of sub-step offsets
#     used per iteration, so no two distinct (condition, iteration, sub-step)
#     triples can ever produce the same literal seed.
#   - test_wide replaced by test cohort's data_long throughout: full_fit
#     evaluation now uses row-level prediction + subject-level aggregation
#     via .predict_subject_level() (cv_strategies.R), matching how naive and
#     subject-level CV predictions are produced. This keeps the prediction
#     method identical across all three (true/naive/subject) evaluations and
#     across both learners (see learners.R v4, cv_strategies.R v3).
#   - fit_learner() and run_cv() calls updated to the new (no data_wide)
#     signatures.
#
# Author: Joseph L. Hagan, ScD, MSPH
# Version: 4 (2026-07-10)
###############################################################################


# ---------------------------------------------------------------------------
# run_one_iteration()
# ---------------------------------------------------------------------------
run_one_iteration <- function(condition_row,
                              iter_id,
                              master_seed,
                              T_const = 30L,
                              n_test  = 5000L,
                              K       = 10L) {

  # Condition spacing 1e6, iteration spacing 100. Sub-step offsets used below
  # (+0..+4) stay far inside the 100-wide gap between consecutive iterations,
  # and the +1e9 offset inside generate_test_cohort() stays far outside the
  # ~5.5e8 ceiling of the largest base seed used directly anywhere in the
  # factorial. No two (condition_id, iter_id, sub-step) triples can collide.
  iter_seed <- master_seed + condition_row$condition_id * 1000000L + iter_id * 100L

  arm          <- condition_row$arm
  n            <- condition_row$n
  p_F          <- condition_row$p_F
  p_L          <- condition_row$p_L
  ICC          <- condition_row$ICC
  rho          <- condition_row$rho
  R2_total     <- condition_row$R2_total
  signal_alloc <- condition_row$signal_allocation

  # Patched 2026-07-09: was hardcoded `target_prev <- 0.62` (incorrect; the
  # actual any-ROP primary-outcome prevalence in the motivating dataset is
  # 48.0%, not 62%; see manuscript Section 2.1 for the correction history).
  # Now read from the condition grid so the same iteration function serves
  # both the corrected main factorial (target_prevalence fixed at 0.480 for
  # every condition) and the supplemental predictor-imbalance/prevalence
  # grid (target_prevalence varies by condition: 0.480 or 0.150).
  target_prev  <- condition_row$target_prevalence
  if (is.null(target_prev) || is.na(target_prev))
    stop("condition_row$target_prevalence is missing or NA for condition_id ",
         condition_row$condition_id,
         ". Every row of the condition grid must specify a target prevalence.")

  learner_types <- c("ridge_linear", "xgboost")

  predictor_set <- switch(arm,
    "A_fixed" = "fixed_only",
    "B_long"  = "longitudinal_only",
    "C_mixed" = "mixed",
    stop("Unknown arm: ", arm)
  )

  ICC_dgm <- if (is.na(ICC)) 0 else ICC
  rho_dgm <- if (is.na(rho)) 0 else rho

  result <- tryCatch({

    # --- 1. Beta vectors and intercept --------------------------------------
    betas <- compute_signal_allocation(
      R2_total   = R2_total,
      allocation = signal_alloc,
      p_F = p_F, p_L = p_L,
      ICC = ICC_dgm, rho = rho_dgm,
      T_const = T_const,
      seed = iter_seed
    )
    beta_F <- betas$beta_F
    beta_L <- betas$beta_L

    alpha <- calibrate_intercept(
      beta_F = beta_F, beta_L = beta_L,
      p_F = p_F, p_L = p_L,
      ICC = ICC_dgm, rho = rho_dgm,
      T_const = T_const,
      target_prev = target_prev,
      seed = iter_seed + 1L
    )

    # --- 2. Training cohort -------------------------------------------------
    cohort <- generate_cohort(
      n = n, T_const = T_const,
      p_F = p_F, p_L = p_L,
      beta_F = beta_F, beta_L = beta_L,
      ICC = ICC_dgm, rho = rho_dgm,
      alpha = alpha,
      seed = iter_seed + 2L
    )

    n_events <- sum(cohort$data_wide$Y)
    if (n_events < 3L || (n - n_events) < 3L)
      stop("Degenerate dataset: n_events = ", n_events)

    data_long   <- cohort$data_long
    data_wide   <- cohort$data_wide
    prev_actual <- cohort$prevalence_actual

    # --- 3. Test cohort -------------------------------------------------------
    test_cohort <- generate_test_cohort(
      n_test = n_test, T_const = T_const,
      p_F = p_F, p_L = p_L,
      beta_F = beta_F, beta_L = beta_L,
      ICC = ICC_dgm, rho = rho_dgm,
      alpha = alpha,
      seed = iter_seed + 3L
    )
    test_long <- test_cohort$data_long

    # --- 4. Uniqueness statistics --------------------------------------------
    uniq_stats <- compute_uniqueness_stats(data_wide, p_F)

    # --- 5. Fold assignments (naive and subject only; LOCO dropped) ---------
    fold_seed     <- iter_seed + 4L
    folds_naive   <- make_naive_folds(data_long,   K = K, seed = fold_seed)
    folds_subject <- make_subject_folds(data_long, K = K, seed = fold_seed)

    # --- 5b. Mechanistic diagnostic: fold-specific subject recognizability --
    # Computed once per iteration (not per learner), before the model-fitting
    # loop. The diagnostic characterises the TRAINING FOLD fingerprint
    # structure that drives leakage, independently of which learner is used.
    # See compute_fold_recognition_stats() in metrics.R for rationale.
    recog_stats <- compute_fold_recognition_stats(
      data_long     = data_long,
      folds_naive   = folds_naive,
      predictor_set = predictor_set,
      p_F = p_F, p_L = p_L,
      K = K
    )

    # --- 6. Per-learner evaluation ------------------------------------------
    # n_subj_full is the full condition's subject count, computed once here
    # (identical to `n`, computed defensively from the actual generated data)
    # and used for every learner's hyperparameter computation below, so that
    # ridge's lambda and XGBoost's min_child_weight are identical across the
    # full-model fit and every fold of both CV strategies (patched 2026-07-10;
    # see learners.R v5 changelog). Previously these were recomputed inside
    # fit_learner() from whichever fold's data it was handed, so subject-level
    # folds (~90% of n_subj_full) received systematically different
    # regularization than the full model and naive folds.
    n_subj_full <- length(unique(data_long$subject_id))

    learner_rows <- lapply(learner_types, function(ltype) {

      t_start    <- proc.time()["elapsed"]
      warn_flags <- character(0)

      hyperparams <- compute_hyperparams(ltype, n_subj_full)

      # Full fit -> AUROC_true on held-out test cohort. Row-level prediction
      # + subject-level aggregation, identical method to naive/subject CV.
      full_fit <- withCallingHandlers(
        fit_learner(data_long, ltype, predictor_set, p_F, p_L,
                    hyperparams = hyperparams),
        warning = function(w) {
          warn_flags <<- c(warn_flags, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      pred_true  <- .predict_subject_level(full_fit, test_long)
      auroc_true <- compute_auroc(pred_true$prob_pred, pred_true$Y)
      brier_true <- compute_brier(pred_true$prob_pred, pred_true$Y)
      calib_true <- compute_calibration(pred_true$prob_pred, pred_true$Y)

      # Naive 10-fold CV
      cv_naive <- withCallingHandlers(
        run_cv(data_long, ltype, predictor_set, p_F, p_L,
               folds_naive, hyperparams = hyperparams, fold_level = "row"),
        warning = function(w) {
          warn_flags <<- c(warn_flags, paste0("naive:", conditionMessage(w)))
          invokeRestart("muffleWarning")
        }
      )
      auroc_naive <- compute_auroc(cv_naive$prob_pred, cv_naive$Y)
      brier_naive <- compute_brier(cv_naive$prob_pred, cv_naive$Y)
      calib_naive <- compute_calibration(cv_naive$prob_pred, cv_naive$Y)

      # Subject-level 10-fold CV
      cv_subject <- withCallingHandlers(
        run_cv(data_long, ltype, predictor_set, p_F, p_L,
               folds_subject, hyperparams = hyperparams, fold_level = "subject"),
        warning = function(w) {
          warn_flags <<- c(warn_flags, paste0("subject:", conditionMessage(w)))
          invokeRestart("muffleWarning")
        }
      )
      auroc_subject <- compute_auroc(cv_subject$prob_pred, cv_subject$Y)
      brier_subject <- compute_brier(cv_subject$prob_pred, cv_subject$Y)
      calib_subject <- compute_calibration(cv_subject$prob_pred, cv_subject$Y)

      runtime_sec <- proc.time()["elapsed"] - t_start

      data.frame(
        # Identity
        arm              = arm,
        condition_id     = condition_row$condition_id,
        iteration        = iter_id,
        learner          = ltype,
        # DGM parameters
        n                = n,
        T_const          = T_const,
        p_F              = p_F,
        p_L              = p_L,
        ICC              = ICC,
        rho              = rho,
        R2_total         = R2_total,
        signal_allocation = signal_alloc,
        prevalence_target = target_prev,
        # DGM realizations
        prevalence_actual = prev_actual,
        pct_unique_q4    = uniq_stats$pct_unique_q4,
        pct_unique_q10   = uniq_stats$pct_unique_q10,
        mean_class_size  = uniq_stats$mean_class_size,
        # Fold-specific descriptive diagnostic (computed once, shared across
        # both learners per iteration). u_F, u_L, and u_joint are reported as
        # descriptive fingerprint statistics; overlap/union/newly_joint are
        # computed directly from per-subject unique-indicators (Arm C only).
        # Patched 2026-07-10: no formula-based reference (neither the simple
        # sum nor the independence-adjusted u_F+u_L-u_F*u_L) is valid for
        # testing the single-mechanism account with this statistic, since
        # joint fingerprint uniqueness reflects complementary disambiguation
        # as much as shared identification capacity (manuscript Section 2.6).
        # The mechanistic test is compute_leakage_credits() (metrics.R),
        # used in the paired additivity sub-study only.
        fold_recog_fixed_q4       = recog_stats$fold_recog_fixed_q4,
        fold_recog_fixed_q10      = recog_stats$fold_recog_fixed_q10,
        fold_recog_long_q4        = recog_stats$fold_recog_long_q4,
        fold_recog_long_q10       = recog_stats$fold_recog_long_q10,
        fold_recog_joint_q4       = recog_stats$fold_recog_joint_q4,
        fold_recog_joint_q10      = recog_stats$fold_recog_joint_q10,
        fold_recog_overlap_q4     = recog_stats$fold_recog_overlap_q4,
        fold_recog_overlap_q10    = recog_stats$fold_recog_overlap_q10,
        fold_recog_union_q4       = recog_stats$fold_recog_union_q4,
        fold_recog_union_q10      = recog_stats$fold_recog_union_q10,
        fold_recog_newly_joint_q4  = recog_stats$fold_recog_newly_joint_q4,
        fold_recog_newly_joint_q10 = recog_stats$fold_recog_newly_joint_q10,
        # True performance (held-out test cohort)
        auroc_true       = auroc_true,
        # CV-based AUROC estimates
        auroc_naive      = auroc_naive,
        auroc_subject    = auroc_subject,
        # Calibration outputs (kept for descriptive/sensitivity use only;
        # NOT used for the calibration follow-up paper, which requires a
        # dedicated simulation that generates genuine per-timepoint
        # predictions and propagates their within-subject correlation into
        # p_bar_i -- see calibration_w_repeats_idea.docx)
        brier_naive      = brier_naive,
        brier_subject    = brier_subject,
        brier_true       = brier_true,
        calib_int_naive  = calib_naive$intercept,
        calib_slope_naive = calib_naive$slope,
        calib_int_subject = calib_subject$intercept,
        calib_slope_subject = calib_subject$slope,
        calib_int_true   = calib_true$intercept,
        calib_slope_true = calib_true$slope,
        # Bookkeeping
        seed             = iter_seed,
        runtime_sec      = runtime_sec,
        warning_flags    = paste(warn_flags, collapse = ";"),
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, learner_rows)

  }, error = function(e) {

    data.frame(
      arm              = arm,
      condition_id     = condition_row$condition_id,
      iteration        = iter_id,
      learner          = "ERROR",
      n = n, T_const = T_const, p_F = p_F, p_L = p_L,
      ICC = ICC, rho = rho, R2_total = R2_total,
      signal_allocation = signal_alloc,
      prevalence_target = target_prev,
      prevalence_actual = NA_real_,
      pct_unique_q4 = NA_real_, pct_unique_q10 = NA_real_,
      mean_class_size = NA_real_,
      fold_recog_fixed_q4  = NA_real_, fold_recog_fixed_q10  = NA_real_,
      fold_recog_long_q4   = NA_real_, fold_recog_long_q10   = NA_real_,
      fold_recog_joint_q4  = NA_real_, fold_recog_joint_q10  = NA_real_,
      fold_recog_overlap_q4      = NA_real_, fold_recog_overlap_q10      = NA_real_,
      fold_recog_union_q4        = NA_real_, fold_recog_union_q10        = NA_real_,
      fold_recog_newly_joint_q4  = NA_real_, fold_recog_newly_joint_q10  = NA_real_,
      auroc_true    = NA_real_,
      auroc_naive   = NA_real_, auroc_subject = NA_real_,
      brier_naive   = NA_real_, brier_subject = NA_real_, brier_true = NA_real_,
      calib_int_naive = NA_real_, calib_slope_naive = NA_real_,
      calib_int_subject = NA_real_, calib_slope_subject = NA_real_,
      calib_int_true = NA_real_, calib_slope_true = NA_real_,
      seed = iter_seed, runtime_sec = NA_real_,
      warning_flags = paste0("ERROR: ", conditionMessage(e)),
      stringsAsFactors = FALSE
    )
  })

  result
}

###############################################################################
# dgm.R — Data-Generating Model
# Mixed-Predictor CV Optimism Decomposition Simulation
#
# Functions:
#   generate_cohort()         — training cohort (long + wide format)
#   generate_test_cohort()    — held-out test cohort (wide format)
#   compute_signal_allocation() — translate R2 + allocation into beta vectors
#   calibrate_intercept()     — bisection to hit target prevalence
#
# Predictor structure:
#   Fixed (time-invariant): X_F_1 ... X_F_pF ~ N(0,1) iid across subjects
#   Longitudinal (time-varying): ICC + AR(1) structure, marginal variance = 1
#
# Outcome:
#   logit P(Y_i=1) = alpha + X_F_i %*% beta_F + X_L_bar_i %*% beta_L
#   Y_i ~ Bernoulli(p_i)  [subject-level binary outcome]
#
# Author: Joseph L. Hagan, ScD, MSPH
# Spec:   5-5-26 mixed predictors spec memo
###############################################################################


# ---------------------------------------------------------------------------
# compute_signal_allocation()
#
# Translates a target R2_total and allocation type into beta_F and beta_L
# in closed form (patched 2026-07-10; see var_longitudinal_mean() below and
# manuscript Section 2.1 for the derivation). Replaces the v3 approach,
# which solved for beta numerically on a simulated reference sample using
# an approximate longitudinal subject-mean variance (ICC + (1-ICC)/T) that
# omitted the AR(1) parameter rho entirely, understating the true variance
# whenever rho > 0 and mislabeling the achieved fixed/longitudinal signal
# split as a result (see var_longitudinal_mean() for the magnitude).
#
# Derivation:
#   V_eta = (pi^2/3) * R2_total / (1 - R2_total)   [target variance of the
#           linear predictor on the logistic latent-variable scale]
#   V_F   = q * V_eta,  V_L = (1 - q) * V_eta        [q = fixed-signal share]
#   beta_F = sqrt(V_F) * d_F                          [fixed predictors have
#                                                       unit variance]
#   beta_L = sqrt(V_L / Var(X_L_bar)) * d_L            [longitudinal subject
#                                                       means do not have
#                                                       unit variance; scaled
#                                                       by the exact AR(1)
#                                                       mean variance]
# where d_F, d_L are unit-length, equal-magnitude direction vectors within
# each predictor type. This guarantees the population latent R2 and the
# population signal allocation equal their nominal targets directly;
# realized values in any finite simulated cohort will vary around these
# targets due to ordinary sampling variability.
#
# allocation: "fixed_dominant"   -> 70% of signal from fixed predictors
#             "long_dominant"    -> 30% of signal from fixed predictors
#             "balanced"         -> 50% of signal from fixed predictors
#             "NA"               -> used for single-arm cases (all signal
#                                   goes to whichever predictor type exists)
#
# Returns: list(beta_F, beta_L)
# ---------------------------------------------------------------------------

# Exact variance of a subject's mean over T equally spaced AR(1)
# observations, with between-subject variance ICC and within-subject
# variance (1 - ICC) at correlation rho. Reduces to ICC + (1-ICC)/T only
# when rho = 0; otherwise larger, since positive autocorrelation slows the
# rate at which the within-subject mean converges as T grows.
var_longitudinal_mean <- function(ICC, rho, T_const) {
  if (T_const <= 1L) return(ICC + (1 - ICC))  # degenerate case, T=1
  h <- seq_len(T_const - 1L)
  within_term <- (1 - ICC) / T_const^2 *
    (T_const + 2 * sum((T_const - h) * rho^h))
  ICC + within_term
}

compute_signal_allocation <- function(R2_total,
                                      allocation,
                                      p_F,
                                      p_L,
                                      ICC  = 0.5,
                                      rho  = 0.5,
                                      T_const = 30,
                                      n_ref = 5000,   # retained for call-site
                                      seed  = 42) {   # compatibility; unused,
                                                       # closed-form has no
                                                       # Monte Carlo step

  # Proportion of total signal assigned to fixed predictors
  prop_fixed <- switch(allocation,
    "fixed_dominant"  = 0.70,
    "long_dominant"   = 0.30,
    "balanced"        = 0.50,
    "NA"              = {
      if (p_F > 0 && p_L == 0) 1.0
      else if (p_F == 0 && p_L > 0) 0.0
      else 0.5   # fallback (should not occur with well-formed condition grid)
    }
  )

  logistic_error_var <- pi^2 / 3
  V_eta <- logistic_error_var * R2_total / (1 - R2_total)
  V_F   <- prop_fixed * V_eta
  V_L   <- (1 - prop_fixed) * V_eta

  # Unit-length, equal-magnitude direction vectors within each predictor type
  dir_F <- if (p_F > 0) rep(1, p_F) / sqrt(p_F) else numeric(0)
  dir_L <- if (p_L > 0) rep(1, p_L) / sqrt(p_L) else numeric(0)

  # Fixed predictors have unit variance by construction (Section 2.1)
  beta_F <- if (p_F > 0) sqrt(V_F) * dir_F else numeric(0)

  # Longitudinal subject means do not have unit variance; scale by the
  # exact AR(1) mean variance so that the realized variance contribution
  # matches V_L directly
  beta_L <- if (p_L > 0) {
    var_xL_bar <- var_longitudinal_mean(ICC, rho, T_const)
    sqrt(V_L / var_xL_bar) * dir_L
  } else {
    numeric(0)
  }

  list(beta_F = beta_F, beta_L = beta_L)
}


# ---------------------------------------------------------------------------
# calibrate_intercept()
#
# Bisection search for the logistic intercept alpha that achieves the
# target outcome prevalence on a large reference dataset.
# ---------------------------------------------------------------------------
calibrate_intercept <- function(beta_F, beta_L,
                                p_F, p_L,
                                ICC, rho, T_const,
                                target_prev = 0.48,  # patched 2026-07-09: was 0.62 (incorrect).
                                                      # Not actually used when called from
                                                      # simulate.R, which always passes
                                                      # target_prev explicitly from
                                                      # condition_row$target_prevalence, but
                                                      # corrected here for direct callers.
                                n_ref = 5000,
                                seed  = 99) {

  set.seed(seed)

  # Generate reference wide-format data
  ref <- .generate_wide_predictors(n     = n_ref,
                                   T_const = T_const,
                                   p_F   = p_F,
                                   p_L   = p_L,
                                   ICC   = ICC,
                                   rho   = rho,
                                   seed  = seed)

  compute_prev <- function(alpha) {
    eta <- alpha
    if (p_F > 0) eta <- eta + as.matrix(ref$X_F) %*% beta_F
    if (p_L > 0) eta <- eta + as.matrix(ref$X_L_bar) %*% beta_L
    mean(plogis(eta))
  }

  lo <- -20; hi <- 20
  for (i in seq_len(80)) {
    mid <- (lo + hi) / 2
    if (compute_prev(mid) < target_prev) lo <- mid else hi <- mid
    if (hi - lo < 1e-7) break
  }
  (lo + hi) / 2
}


# ---------------------------------------------------------------------------
# .generate_wide_predictors()  [internal helper]
#
# Generates subject-level fixed predictors and subject means of longitudinal
# predictors. Used by both calibrate_intercept() and generate_cohort().
# ---------------------------------------------------------------------------
.generate_wide_predictors <- function(n, T_const, p_F, p_L, ICC, rho, seed) {

  set.seed(seed)

  # Fixed predictors: iid N(0,1) across subjects
  X_F <- if (p_F > 0) {
    as.data.frame(matrix(rnorm(n * p_F), nrow = n,
                         dimnames = list(NULL, paste0("X_F_", seq_len(p_F)))))
  } else {
    data.frame(row.names = seq_len(n))
  }

  # Longitudinal predictors: ICC + AR(1) structure
  # Between-subject variance: sigma2_b = ICC  (marginal variance = 1)
  # Within-subject innovation variance: sigma2_e = 1 - ICC
  # AR(1) innovation sd: sigma_innov = sqrt(sigma2_e * (1 - rho^2))
  X_L_bar <- if (p_L > 0) {
    sigma2_b    <- ICC
    sigma2_e    <- 1 - ICC
    sigma_innov <- sqrt(sigma2_e * (1 - rho^2))

    xL_bar_mat <- matrix(NA_real_, nrow = n, ncol = p_L)
    for (j in seq_len(p_L)) {
      b_i <- rnorm(n, mean = 0, sd = sqrt(sigma2_b))  # subject random effect
      xL_sum <- numeric(n)
      e_prev <- rnorm(n, mean = 0, sd = sqrt(sigma2_e))  # t=1 error
      xL_sum <- xL_sum + (b_i + e_prev)
      for (t in 2:T_const) {
        e_t    <- rho * e_prev + rnorm(n, 0, sigma_innov)
        xL_sum <- xL_sum + (b_i + e_t)
        e_prev <- e_t
      }
      xL_bar_mat[, j] <- xL_sum / T_const
    }
    as.data.frame(xL_bar_mat,
                  col.names = paste0("X_L_bar_", seq_len(p_L)))
  } else {
    data.frame(row.names = seq_len(n))
  }

  list(X_F = X_F, X_L_bar = X_L_bar)
}


# ---------------------------------------------------------------------------
# generate_cohort()
#
# Generates the training cohort in both long and wide format.
#
# Arguments:
#   n         — number of subjects
#   T_const   — time points per subject (constant)
#   p_F, p_L  — number of fixed / longitudinal predictors
#   beta_F, beta_L — coefficient vectors (from compute_signal_allocation)
#   ICC, rho  — longitudinal predictor dependence parameters
#   alpha     — logistic intercept (from calibrate_intercept)
#   seed      — integer for reproducibility
#
# Returns: list(data_long, data_wide, prevalence_actual)
# ---------------------------------------------------------------------------
generate_cohort <- function(n, T_const, p_F, p_L,
                            beta_F, beta_L,
                            ICC, rho, alpha,
                            seed) {

  set.seed(seed)

  sigma2_b    <- if (p_L > 0) ICC       else NA_real_
  sigma2_e    <- if (p_L > 0) 1 - ICC   else NA_real_
  sigma_innov <- if (p_L > 0) sqrt(sigma2_e * (1 - rho^2)) else NA_real_

  # ---- Fixed predictors (subject level) -----------------------------------
  X_F <- if (p_F > 0) {
    matrix(rnorm(n * p_F), nrow = n,
           dimnames = list(NULL, paste0("X_F_", seq_len(p_F))))
  } else {
    matrix(nrow = n, ncol = 0)
  }

  # ---- Subject-level linear predictor contribution from fixed predictors --
  eta_F <- if (p_F > 0) as.numeric(X_F %*% beta_F) else rep(0, n)

  # ---- Longitudinal predictors + outcome ----------------------------------
  # We generate the full long-format dataset row by row over subjects.
  # Subject means (X_L_bar) are accumulated for the wide format.

  long_rows  <- vector("list", n)
  X_L_bar    <- if (p_L > 0) matrix(NA_real_, nrow = n, ncol = p_L) else
                              matrix(nrow = n, ncol = 0)

  # Subject-level random effects for longitudinal predictors
  b_mat <- if (p_L > 0) {
    matrix(rnorm(n * p_L, 0, sqrt(sigma2_b)), nrow = n, ncol = p_L)
  } else NULL

  for (i in seq_len(n)) {

    if (p_L > 0) {
      # AR(1) process for each longitudinal predictor
      xL_it <- matrix(NA_real_, nrow = T_const, ncol = p_L)
      e_prev <- rnorm(p_L, 0, sqrt(sigma2_e))
      xL_it[1, ] <- b_mat[i, ] + e_prev
      for (t in 2:T_const) {
        e_t        <- rho * e_prev + rnorm(p_L, 0, sigma_innov)
        xL_it[t, ] <- b_mat[i, ] + e_t
        e_prev     <- e_t
      }
      X_L_bar[i, ] <- colMeans(xL_it)
      eta_L_i <- sum(X_L_bar[i, ] * beta_L)
    } else {
      xL_it   <- matrix(nrow = T_const, ncol = 0)
      eta_L_i <- 0
    }

    # Subject-level outcome
    p_i <- plogis(alpha + eta_F[i] + eta_L_i)
    Y_i <- rbinom(1, 1, p_i)

    # Long-format rows for subject i
    row_df <- data.frame(
      subject_id = i,
      time       = seq_len(T_const),
      Y          = Y_i   # replicated across time points
    )
    if (p_F > 0) {
      xF_rep <- matrix(rep(X_F[i, ], each = T_const),
                       nrow = T_const,
                       dimnames = list(NULL, paste0("X_F_", seq_len(p_F))))
      row_df <- cbind(row_df, xF_rep)
    }
    if (p_L > 0) {
      colnames(xL_it) <- paste0("X_L_", seq_len(p_L))
      row_df <- cbind(row_df, xL_it)
    }

    long_rows[[i]] <- row_df
  }

  data_long <- do.call(rbind, long_rows)
  rownames(data_long) <- NULL

  # ---- Wide format (one row per subject) ----------------------------------
  data_wide <- data.frame(subject_id = seq_len(n))
  data_wide$Y <- vapply(seq_len(n), function(i)
    data_long$Y[data_long$subject_id == i][1L], numeric(1))

  if (p_F > 0) {
    xF_df <- as.data.frame(X_F)
    colnames(xF_df) <- paste0("X_F_", seq_len(p_F))
    data_wide <- cbind(data_wide, xF_df)
  }
  if (p_L > 0) {
    xLbar_df <- as.data.frame(X_L_bar)
    colnames(xLbar_df) <- paste0("X_L_bar_", seq_len(p_L))
    data_wide <- cbind(data_wide, xLbar_df)
  }

  list(
    data_long        = data_long,
    data_wide        = data_wide,
    prevalence_actual = mean(data_wide$Y)
  )
}


# ---------------------------------------------------------------------------
# generate_test_cohort()
#
# Generates the held-out test cohort used to estimate AUROC_true.
# Returns BOTH long and wide format: row-level predictions (used by
# fit$predict_prob()) require data_long; data_wide is retained for any
# descriptive use.
# n_test = 5000 as specified in the memo.
#
# Seed offset: +1e9 (was +1e7). With the iter_seed redesign in simulate.R
# (condition spacing 1e6, iteration spacing 100), the largest base seed any
# single iteration can use directly is on the order of 5.5e8. A +1e9 offset
# guarantees this internal call can never land on a seed value used directly
# by ANY (condition_id, iter_id) combination in the factorial, eliminating
# the cross-purpose seed collisions present in the prior design.
# ---------------------------------------------------------------------------
generate_test_cohort <- function(n_test = 5000,
                                 T_const, p_F, p_L,
                                 beta_F, beta_L,
                                 ICC, rho, alpha,
                                 seed) {

  # Offset seed so test cohort is independent of training cohort
  cohort <- generate_cohort(n       = n_test,
                            T_const = T_const,
                            p_F     = p_F,
                            p_L     = p_L,
                            beta_F  = beta_F,
                            beta_L  = beta_L,
                            ICC     = ICC,
                            rho     = rho,
                            alpha   = alpha,
                            seed    = seed + 1000000000L)
  list(data_long = cohort$data_long, data_wide = cohort$data_wide)
}

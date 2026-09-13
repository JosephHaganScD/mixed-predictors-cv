###############################################################################
# Penalty sensitivity analysis  (manuscript Section 3.9)
#
# The main factorial used lambda = log(1 + 1/n)/n, which is of order n^-2 and
# imposes negligible shrinkage. Pre-submission review objected that results
# attributed to "ridge logistic regression" are effectively unpenalized.
#
# Question: does naive-CV optimism depend on the degree of shrinkage?
# Prediction: shrinkage toward zero reduces a linear learner's capacity to fit
# subject-specific fingerprints, so naive optimism should DECREASE monotonically
# in lambda, while subject-level optimism should be comparatively unaffected.
#
# Only the linear learner is evaluated; the penalty is a property of that
# learner and XGBoost is unaffected by it.
#
# EFFICIENCY NOTE: glmnet fits an entire lambda path in a single call, so all
# penalty levels are obtained from one fit per fold rather than one fit per
# (fold, lambda). This is roughly a fivefold saving over the naive loop.
#
# Author: J. Hagan   Date: 2026-07-19
###############################################################################

library(glmnet); library(pROC); library(dplyr)

## Warnings are collected and tabulated at the end rather than deferred to
## warnings(), which truncates at 50 and loses context.
##
## HISTORY: an earlier version of this script fitted a dense lambda path and
## extracted predictions with predict(s = lam_target). That invoked glmnet's
## lambda.interp, which raised one approx() warning per prediction and resolved
## all six target lambdas to a single effective value, so every penalty level
## returned identical coefficients and the analysis measured nothing. Predictions
## are now formed from coef() at the fitted lambdas directly, with no
## interpolation path of any kind, and a coefficient-norm diagnostic verifies
## that shrinkage is actually taking effect.
WARN <- new.env(); WARN$msg <- character(0)
with_warn_log <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    WARN$msg <- c(WARN$msg, conditionMessage(w)); invokeRestart("muffleWarning")
  })
}

# Set OUT_DIR to the folder where results should be saved.
OUT_DIR <- "."   # change to your preferred output directory
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
RESUME      <- TRUE
MASTER_SEED <- 20260719
N_ITER      <- 200
K_FOLDS     <- 10
T_CONST     <- 30
N_TEST      <- 5000
PREV_TARGET <- 0.48

## --- fixed penalty grid, spanning negligible to substantial ------------------
## The first element is the rule used in the main factorial and is condition
## dependent, so it is appended per condition rather than stored here.
## A common near-zero anchor (1e-4) is included so that every condition shares
## the SAME five penalty levels. Without it the main-rule lambda differs by n
## (3.96e-4 at n = 50, 2.49e-5 at n = 200), and any marginal summary over lambda
## averages different condition sets in different rows.
LAMBDA_FIXED <- c(1e-4, 0.001, 0.01, 0.05, 0.20)

## ---------------------------------------------------------------------------
## Data-generating model (identical to the v8 main factorial, Section 2.1)
## ---------------------------------------------------------------------------
var_xbar <- function(T_const, ICC, rho) {
  if (ICC >= 1) return(1)
  h <- seq_len(T_const - 1)
  ICC + (1 - ICC) / T_const^2 * (T_const + 2 * sum((T_const - h) * rho^h))
}

gen_cohort <- function(n, p_F, p_L, T_const, ICC, rho, R2, q, alpha = NULL) {
  XF <- if (p_F > 0) matrix(rnorm(n * p_F), n, p_F) else matrix(0, n, 0)
  XL <- array(0, dim = c(n, T_const, max(p_L, 1)))
  if (p_L > 0) {
    s2 <- (1 - ICC) * (1 - rho^2)
    for (l in seq_len(p_L)) {
      b <- rnorm(n, 0, sqrt(ICC))
      e <- matrix(0, n, T_const); e[, 1] <- rnorm(n, 0, sqrt(1 - ICC))
      for (t in 2:T_const) e[, t] <- rho * e[, t - 1] + rnorm(n, 0, sqrt(s2))
      XL[, , l] <- b + e
    }
  }
  V_eta <- (pi^2 / 3) * R2 / (1 - R2)
  bF <- if (p_F > 0) rep(sqrt(q * V_eta / p_F), p_F) else numeric(0)
  bL <- if (p_L > 0) rep(sqrt((1 - q) * V_eta / (p_L * var_xbar(T_const, ICC, rho))), p_L) else numeric(0)

  # Built additively rather than as a single matrix product, because one of the
  # two blocks is always empty in this design (Arm A has p_L = 0, Arm B has
  # p_F = 0) and zero-column matrix products are a needless failure mode.
  eta0 <- rep(0, n)
  if (p_F > 0) eta0 <- eta0 + as.vector(XF %*% bF)
  if (p_L > 0) {
    xbar <- apply(XL[, , seq_len(p_L), drop = FALSE], c(1, 3), mean)
    eta0 <- eta0 + as.vector(xbar %*% bL)
  }

  if (is.null(alpha)) {                              # bisection on the intercept
    lo <- -20; hi <- 20
    for (i in 1:60) {
      mid <- (lo + hi) / 2
      if (mean(plogis(mid + eta0)) < PREV_TARGET) lo <- mid else hi <- mid
    }
    alpha <- (lo + hi) / 2
  }
  Y <- rbinom(n, 1, plogis(alpha + eta0))

  # expand to one row per (subject, time)
  idx <- rep(seq_len(n), each = T_const)
  Xrow <- cbind(
    if (p_F > 0) XF[idx, , drop = FALSE] else NULL,
    if (p_L > 0) matrix(aperm(XL[, , seq_len(p_L), drop = FALSE], c(2, 1, 3)),
                        nrow = n * T_const, ncol = p_L) else NULL
  )
  colnames(Xrow) <- paste0("v", seq_len(ncol(Xrow)))
  list(X = Xrow, subj = idx, Y = Y, Yrow = Y[idx], alpha = alpha, n = n)
}

## ---------------------------------------------------------------------------
## Fit at EXACTLY the target lambdas and form predictions from coef() directly.
## predict(s = ...) is never called, so lambda.interp is never invoked. Columns
## of coef() are aligned to f$lambda, which is sorted decreasing; the assertion
## confirms the alignment rather than assuming it.
## Returns subject-level predictions and the L2 norm of the slope coefficients
## at each lambda, the latter purely as a check that shrinkage is taking effect.
## ---------------------------------------------------------------------------
fit_path <- function(Xtr, ytr, Xte, subj_te, lam_target) {
  lam <- sort(lam_target, decreasing = TRUE)
  f <- glmnet(Xtr, ytr, family = "binomial", alpha = 0, lambda = lam,
              standardize = TRUE)
  stopifnot(length(f$lambda) == length(lam),
            max(abs(f$lambda - lam)) < 1e-12 * max(lam))
  B   <- as.matrix(coef(f))                       # (p + 1) x nlambda
  eta <- cbind(1, Xte) %*% B
  p   <- 1 / (1 + exp(-eta))
  ord <- match(lam_target, lam)                   # restore caller's order
  list(pred = sapply(ord, function(j) tapply(p[, j], subj_te, mean)),
       l2   = sqrt(colSums(B[-1, , drop = FALSE]^2))[ord])
}

## ---------------------------------------------------------------------------
## AUROC with an EXPLICITLY FIXED direction, matching metrics.R::compute_auroc()
## in the v8 main factorial.
##
## Do not revert this to pROC's default. direction = "auto" chooses the
## direction from the data, which folds estimates below 0.5 upward and biases
## AUROC upward by roughly E|AUC - 0.5| whenever true performance is near
## chance. Subject-level CV at n = 50 sits in exactly that regime: the bias was
## measured at +0.09 AUROC units there, and it silently reversed the sign of
## subject-level optimism.
## ---------------------------------------------------------------------------
metrics <- function(pred, ytrue) {
  ok <- length(unique(ytrue)) == 2
  c(auroc = if (ok) suppressMessages(as.numeric(
              pROC::auc(ytrue, pred, direction = "<", quiet = TRUE))) else NA_real_,
    brier = mean((pred - ytrue)^2))
}

## ---------------------------------------------------------------------------
## One iteration: naive CV, subject CV, and true performance, for all lambdas
## ---------------------------------------------------------------------------
one_iter <- function(cond, lambdas) {
  d  <- with(cond, gen_cohort(n, p_F, p_L, T_CONST, ICC, rho, R2, q))
  dt <- with(cond, gen_cohort(N_TEST, p_F, p_L, T_CONST, ICC, rho, R2, q, alpha = d$alpha))
  nl <- length(lambdas)
  sid <- sort(unique(d$subj)); Ysub <- d$Y[sid]

  # --- true out-of-sample on independent test cohort
  ft   <- fit_path(d$X, d$Yrow, dt$X, dt$subj, lambdas)
  pt   <- ft$pred
  mtrue <- t(sapply(seq_len(nl), function(j) metrics(pt[, j], dt$Y)))

  # --- naive CV: folds assigned to rows
  fr <- sample(rep_len(seq_len(K_FOLDS), nrow(d$X)))
  pn <- matrix(NA_real_, length(sid), nl)
  acc <- vector("list", K_FOLDS)
  for (k in seq_len(K_FOLDS)) {
    tr <- fr != k
    acc[[k]] <- list(p = fit_path(d$X[tr, , drop = FALSE], d$Yrow[tr],
                                  d$X[!tr, , drop = FALSE], d$subj[!tr], lambdas)$pred,
                     s = sort(unique(d$subj[!tr])))
  }
  for (j in seq_len(nl)) {                       # average across folds by subject
    num <- rep(0, length(sid)); den <- rep(0, length(sid))
    for (k in seq_len(K_FOLDS)) {
      m <- match(acc[[k]]$s, sid)
      num[m] <- num[m] + acc[[k]]$p[, j]; den[m] <- den[m] + 1
    }
    pn[, j] <- num / den
  }
  mnaive <- t(sapply(seq_len(nl), function(j) metrics(pn[, j], Ysub)))

  # --- subject-level CV: folds assigned to subjects
  fs <- sample(rep_len(seq_len(K_FOLDS), length(sid)))
  ps <- matrix(NA_real_, length(sid), nl)
  for (k in seq_len(K_FOLDS)) {
    hold <- sid[fs == k]; tr <- !(d$subj %in% hold)
    pk <- fit_path(d$X[tr, , drop = FALSE], d$Yrow[tr],
                   d$X[!tr, , drop = FALSE], d$subj[!tr], lambdas)$pred
    ps[match(sort(hold), sid), ] <- pk
  }
  msub <- t(sapply(seq_len(nl), function(j) metrics(ps[, j], Ysub)))

  data.frame(lambda = lambdas, coef_l2 = ft$l2,
             auroc_true = mtrue[, 1], auroc_naive = mnaive[, 1], auroc_subject = msub[, 1],
             brier_true = mtrue[, 2], brier_naive = mnaive[, 2], brier_subject = msub[, 2])
}

## ---------------------------------------------------------------------------
## Design: 8 fixed-only conditions spanning the Arm A extremes, plus 4
## longitudinal-only conditions to confirm the result is not specific to Arm A.
## ---------------------------------------------------------------------------
condA <- expand.grid(p_F = c(2, 10), n = c(50, 200), R2 = c(0.05, 0.30),
                     stringsAsFactors = FALSE)
condA <- transform(condA, arm = "A_fixed", p_L = 0, ICC = 0.7, rho = 0.3, q = 1)
condB <- expand.grid(p_L = c(2, 5), n = c(50, 200), stringsAsFactors = FALSE)
condB <- transform(condB, arm = "B_long", p_F = 0, R2 = 0.15, ICC = 0.7, rho = 0.3, q = 0)
cond  <- rbind(condA[, c("arm","p_F","p_L","n","R2","ICC","rho","q")],
               condB[, c("arm","p_F","p_L","n","R2","ICC","rho","q")])
cond$condition_id <- seq_len(nrow(cond))
cat("Conditions:", nrow(cond), " Iterations:", N_ITER,
    " Penalty levels:", length(LAMBDA_FIXED) + 1, "\n")

ck  <- file.path(OUT_DIR, "checkpoint.rds")
res_file <- file.path(OUT_DIR, "penalty_results.csv")
start <- 1
if (RESUME && file.exists(ck)) { start <- readRDS(ck) + 1; cat("Resuming at", start, "\n") }

t0 <- Sys.time()
for (i in seq(start, nrow(cond))) {
  g <- cond[i, ]
  lam <- sort(unique(c(log(1 + 1 / g$n) / g$n, LAMBDA_FIXED)), decreasing = TRUE)
  acc <- vector("list", N_ITER)
  for (it in seq_len(N_ITER)) {
    set.seed(MASTER_SEED + 1000L * g$condition_id + it)
    acc[[it]] <- cbind(iter = it, with_warn_log(one_iter(g, lam)))
  }
  out <- bind_rows(acc) |>
    mutate(condition_id = g$condition_id, arm = g$arm, p_F = g$p_F, p_L = g$p_L,
           n = g$n, R2 = g$R2,
           is_main_rule = abs(lambda - log(1 + 1 / g$n) / g$n) < 1e-12)
  write.table(out, res_file, sep = ",", row.names = FALSE,
              col.names = !file.exists(res_file), append = file.exists(res_file))
  saveRDS(i, ck)
  cat(sprintf("  %d/%d  elapsed %.1f min\n", i, nrow(cond),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## ---------------------------------------------------------------------------
## Summary against the pre-specified prediction
## ---------------------------------------------------------------------------
cat("\n=== Warning summary ===\n")
if (length(WARN$msg) == 0) cat("No warnings raised.\n") else {
  cat("Total warnings:", length(WARN$msg), "across",
      nrow(cond) * N_ITER, "iterations\n")
  print(head(sort(table(sub("[0-9.eE+-]{3,}", "<num>", WARN$msg)),
                  decreasing = TRUE), 10))
}

## ---------------------------------------------------------------------------
## SELF-VALIDATION: the main-rule lambda reproduces the main factorial's
## procedure, so those cells should reproduce the published v8 Arm A and Arm B
## values. A systematic discrepancy would indicate that fitting the main-rule
## lambda as one point on a warm-started path differs from fitting it alone.
## ---------------------------------------------------------------------------
published <- data.frame(
  arm = c(rep("A_fixed", 8), rep("B_long", 4)),
  p_F = c(2,2,2,2,10,10,10,10, 0,0,0,0),
  p_L = c(rep(0, 8), 2,2,5,5),
  R2  = c(0.05,0.30,0.05,0.30,0.05,0.30,0.05,0.30, rep(0.15, 4)),
  n   = c(50,50,200,200,50,50,200,200, 50,200,50,200),
  pub_naive   = c(0.0626,0.0161,0.0122,0.0041,0.2281,0.1636,0.0884,0.0438,
                  0.0301,0.0091,0.1029,0.0257),
  pub_subject = c(-0.0697,-0.0456,-0.0270,-0.0104,-0.0395,-0.0414,-0.0155,-0.0118,
                  -0.0708,-0.0141,-0.0510,-0.0205))

r <- read.csv(res_file) |>
  mutate(d_auroc_naive = auroc_naive - auroc_true,
         d_auroc_subj  = auroc_subject - auroc_true,
         d_brier_naive = brier_true - brier_naive)

cat("\n=== GATE CHECK: is the penalty actually taking effect? ===\n")
cat("Mean L2 norm of slope coefficients by lambda. MUST decrease monotonically.\n")
cat("If these are flat, the run is invalid and nothing below should be read.\n")
gate <- r |> group_by(arm, p_F, p_L, n, R2, lambda) |>
  summarise(coef_l2 = mean(coef_l2), .groups = "drop") |>
  group_by(condition_id = paste(arm, p_F, p_L, n, R2)) |>
  arrange(lambda, .by_group = TRUE) |>
  summarise(l2_min_lambda = round(first(coef_l2), 4),
            l2_max_lambda = round(last(coef_l2), 6),
            ratio = round(last(coef_l2) / first(coef_l2), 5), .groups = "drop")
print(as.data.frame(gate))
cat("\nAll ratios should be far below 1. Max ratio observed:",
    round(max(gate$ratio), 5), "\n")

cat("\n=== Common penalty grid only: balanced marginal summary ===\n")
cat("Restricted to the levels every condition shares, so rows are comparable.\n")
print(r |> filter(lambda %in% LAMBDA_FIXED) |> group_by(arm, lambda) |>
        summarise(n_cells = n(),
                  d_auroc_naive = round(mean(d_auroc_naive, na.rm = TRUE), 4),
                  d_auroc_subj  = round(mean(d_auroc_subj, na.rm = TRUE), 4),
                  d_brier_naive = round(mean(d_brier_naive), 4), .groups = "drop") |>
        arrange(arm, lambda) |> as.data.frame())

cat("\n=== Main-rule lambda reported separately (differs by n, not comparable above) ===\n")
print(r |> filter(is_main_rule) |> group_by(arm, n, lambda) |>
        summarise(n_cells = n(),
                  d_auroc_naive = round(mean(d_auroc_naive, na.rm = TRUE), 4),
                  d_brier_naive = round(mean(d_brier_naive), 4), .groups = "drop") |>
        as.data.frame())

cat("\n=== PRIMARY TEST: within-condition change relative to the main rule ===\n")
cat("Prediction: negative and increasingly negative as lambda rises.\n")
base <- r |> filter(is_main_rule) |>
  group_by(condition_id) |>
  summarise(base_naive = mean(d_auroc_naive, na.rm = TRUE),
            base_subj  = mean(d_auroc_subj, na.rm = TRUE), .groups = "drop")
print(r |> filter(lambda %in% LAMBDA_FIXED) |>
        group_by(condition_id, arm, p_F, p_L, n, R2, lambda) |>
        summarise(d_auroc_naive = mean(d_auroc_naive, na.rm = TRUE),
                  d_auroc_subj  = mean(d_auroc_subj, na.rm = TRUE), .groups = "drop") |>
        left_join(base, by = "condition_id") |>
        mutate(chg_naive = d_auroc_naive - base_naive,
               chg_subj  = d_auroc_subj - base_subj) |>
        group_by(arm, lambda) |>
        summarise(mean_chg_naive = round(mean(chg_naive), 4),
                  pct_conditions_decreasing = round(mean(chg_naive < 0), 3),
                  mean_chg_subject = round(mean(chg_subj), 4), .groups = "drop") |>
        as.data.frame())

cat("\n=== Per-condition detail: main rule versus heaviest penalty ===\n")
print(r |> filter(is_main_rule | lambda == max(LAMBDA_FIXED)) |>
        group_by(arm, p_F, p_L, n, R2, lambda, is_main_rule) |>
        summarise(d_auroc_naive = round(mean(d_auroc_naive, na.rm = TRUE), 4),
                  d_auroc_subj = round(mean(d_auroc_subj, na.rm = TRUE), 4),
                  .groups = "drop") |> as.data.frame())

cat("\n=== SELF-VALIDATION against published v8 values (main-rule lambda) ===\n")
cat("z = (observed - published) / SE(observed). |z| > 3 warrants investigation.\n")
chk <- r |> filter(is_main_rule) |>
  group_by(arm, p_F, p_L, n, R2) |>
  summarise(obs_naive = mean(d_auroc_naive, na.rm = TRUE),
            se_naive  = sd(d_auroc_naive, na.rm = TRUE) / sqrt(sum(!is.na(d_auroc_naive))),
            obs_subj  = mean(d_auroc_subj, na.rm = TRUE),
            se_subj   = sd(d_auroc_subj, na.rm = TRUE) / sqrt(sum(!is.na(d_auroc_subj))),
            .groups = "drop") |>
  left_join(published, by = c("arm", "p_F", "p_L", "n", "R2")) |>
  mutate(z_naive = round((obs_naive - pub_naive) / se_naive, 2),
         z_subj  = round((obs_subj - pub_subject) / se_subj, 2),
         across(c(obs_naive, pub_naive, obs_subj, pub_subject), \(x) round(x, 4)))
print(as.data.frame(chk[, c("arm","p_F","p_L","n","R2","obs_naive","pub_naive","z_naive",
                            "obs_subj","pub_subject","z_subj")]))
cat("\nMean z (naive):", round(mean(chk$z_naive, na.rm = TRUE), 3),
    " SE:", round(1 / sqrt(sum(!is.na(chk$z_naive))), 3),
    " Max |z|:", round(max(abs(chk$z_naive), na.rm = TRUE), 2), "\n")

cat("\nWritten to:", res_file, "\n")
cat("Total runtime:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

###############################################################################
# Categorical fixed-predictor sensitivity factorial  (manuscript Section 3.11)
#
# The main factorial evaluated continuous fixed predictors only. Many clinical
# baseline covariates are binary or categorical, and k binary predictors
# generate at most 2^k equivalence classes, so they cannot isolate individual
# subjects in cohorts of the sizes evaluated here.
#
# Prediction: identity-mediated leakage declines sharply as the fixed-predictor
# set becomes categorical, and the decline is steepest for the flexible learner,
# which has no individual fingerprint left to partition toward. If confirmed,
# the manuscript's claim sharpens from "baseline covariates leak" to
# "CONTINUOUS baseline covariates leak", which is more specific, more
# actionable, and directly consistent with the empirical illustration in which
# birth weight and gestational age alone identified every subject.
#
# Binary predictors are generated at prevalence 0.5, which maximizes the number
# of distinct equivalence classes and is therefore the CONSERVATIVE choice: it
# gives categorical predictors their best chance of producing leakage.
#
# Total predictor count is held at p_F = 5 across all compositions so that
# composition is not confounded with dimensionality.
#
# Author: J. Hagan   Date: 2026-07-19
###############################################################################

library(glmnet); library(xgboost); library(pROC); library(dplyr)
library(future); library(furrr)

# Set OUT_DIR to the folder where results should be saved.
OUT_DIR <- "."   # change to your preferred output directory
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
RESUME      <- TRUE
MASTER_SEED <- 20260719
N_ITER      <- 500
K_FOLDS     <- 10
T_CONST     <- 30
N_TEST      <- 5000
PREV_TARGET <- 0.48
N_WORKERS   <- 4     # keep low if another simulation is running concurrently

P_F <- 5
COMPOSITIONS <- list(
  "5 continuous"              = c(rep("cont", 5)),
  "3 continuous + 2 binary"   = c(rep("cont", 3), rep("bin", 2)),
  "2 continuous + 3 binary"   = c(rep("cont", 2), rep("bin", 3)),   # ROP-like
  "1 continuous + 4 binary"   = c("cont", rep("bin", 4)),
  "5 binary"                  = c(rep("bin", 5))
)

## ---------------------------------------------------------------------------
## Generate fixed predictors of specified types, each scaled to unit variance so
## that the closed-form coefficient derivation of Section 2.1 applies unchanged.
## ---------------------------------------------------------------------------
gen_fixed <- function(n, types, p_bin = 0.5) {
  X <- matrix(0, n, length(types))
  for (j in seq_along(types)) {
    X[, j] <- if (types[j] == "cont") rnorm(n) else
      (rbinom(n, 1, p_bin) - p_bin) / sqrt(p_bin * (1 - p_bin))
  }
  X
}

gen_cohort <- function(n, types, R2, alpha = NULL) {
  XF <- gen_fixed(n, types)
  V_eta <- (pi^2 / 3) * R2 / (1 - R2)
  bF <- rep(sqrt(V_eta / length(types)), length(types))    # unit-variance columns
  eta0 <- as.vector(XF %*% bF)
  if (is.null(alpha)) {
    lo <- -20; hi <- 20
    for (i in 1:60) {
      mid <- (lo + hi) / 2
      if (mean(plogis(mid + eta0)) < PREV_TARGET) lo <- mid else hi <- mid
    }
    alpha <- (lo + hi) / 2
  }
  Y <- rbinom(n, 1, plogis(alpha + eta0))
  idx <- rep(seq_len(n), each = T_CONST)
  Xrow <- XF[idx, , drop = FALSE]; colnames(Xrow) <- paste0("v", seq_len(ncol(Xrow)))
  list(X = Xrow, XF = XF, subj = idx, Y = Y, Yrow = Y[idx], alpha = alpha, n = n)
}

## ---------------------------------------------------------------------------
## Equivalence-class structure of the fixed-predictor set.
## For categorical compositions the raw values are already discrete, so raw
## uniqueness is the quantity of primary interest; Q4 binning is retained for
## comparability with the main factorial.
## ---------------------------------------------------------------------------
uniq_stats <- function(XF) {
  key_raw <- apply(round(XF, 10), 1, paste, collapse = "|")
  q4 <- apply(XF, 2, function(z) {
    br <- unique(quantile(z, seq(0, 1, 0.25))); br[1] <- -Inf; br[length(br)] <- Inf
    if (length(br) < 3) rep(1L, length(z)) else
      as.integer(cut(z, breaks = br, labels = FALSE, include.lowest = TRUE))
  })
  key_q4 <- apply(q4, 1, paste, collapse = "|")
  c(u_raw     = mean(table(key_raw)[key_raw] == 1),
    u_q4      = mean(table(key_q4)[key_q4] == 1),
    n_classes = length(unique(key_raw)),
    max_class = max(table(key_raw)))
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
## Learners. Hyperparameters computed once from the full condition's subject
## count and held fixed across the full fit and every fold (Section 2.3).
## ---------------------------------------------------------------------------
## The XGBoost call uses xgb.train() with xgb.DMatrix, copied from
## learners.R::.fit_xgboost() in the v8 main factorial. Do not substitute the
## xgboost() convenience wrapper: recent versions of the package renamed data,
## label, eta and verbose, and require a factor y for classification objectives,
## which is what caused the 19 July run to abort at the first XGBoost condition.
## The parameter list must also match v8 exactly, including subsample and
## colsample_bytree, or the "5 continuous" composition will not reproduce the
## published Arm A values and the self-validation below will fail.
fit_predict <- function(learner, Xtr, ytr, Xte, subj_te, lam, mcw) {
  if (learner == "logistic") {
    f <- glmnet(Xtr, ytr, family = "binomial", alpha = 0, lambda = lam)
    p <- as.numeric(predict(f, newx = Xte, s = lam, type = "response"))
  } else {
    params <- list(objective = "binary:logistic", eval_metric = "logloss",
                   max_depth = 4L, eta = 0.1, subsample = 0.8,
                   colsample_bytree = 0.8, min_child_weight = mcw, nthread = 1L)
    f <- xgb.train(params = params, data = xgb.DMatrix(data = Xtr, label = ytr),
                   nrounds = 50L, verbose = 0)
    p <- as.numeric(predict(f, xgb.DMatrix(data = Xte)))
  }
  tapply(p, subj_te, mean)
}

one_iter <- function(types, n, R2, learner) {
  d  <- gen_cohort(n, types, R2)
  dt <- gen_cohort(N_TEST, types, R2, alpha = d$alpha)
  lam <- log(1 + 1 / n) / n; mcw <- max(1, floor(n / 20))
  sid <- sort(unique(d$subj)); Ysub <- d$Y[sid]

  pt <- fit_predict(learner, d$X, d$Yrow, dt$X, dt$subj, lam, mcw)
  mt <- metrics(as.numeric(pt), dt$Y)

  fr <- sample(rep_len(seq_len(K_FOLDS), nrow(d$X)))
  num <- rep(0, length(sid)); den <- rep(0, length(sid))
  for (k in seq_len(K_FOLDS)) {
    tr <- fr != k
    pk <- fit_predict(learner, d$X[tr, , drop = FALSE], d$Yrow[tr],
                      d$X[!tr, , drop = FALSE], d$subj[!tr], lam, mcw)
    m <- match(as.integer(names(pk)), sid)
    num[m] <- num[m] + as.numeric(pk); den[m] <- den[m] + 1
  }
  mn <- metrics(num / den, Ysub)

  fs <- sample(rep_len(seq_len(K_FOLDS), length(sid)))
  ps <- rep(NA_real_, length(sid))
  for (k in seq_len(K_FOLDS)) {
    hold <- sid[fs == k]; tr <- !(d$subj %in% hold)
    pk <- fit_predict(learner, d$X[tr, , drop = FALSE], d$Yrow[tr],
                      d$X[!tr, , drop = FALSE], d$subj[!tr], lam, mcw)
    ps[match(as.integer(names(pk)), sid)] <- as.numeric(pk)
  }
  msu <- metrics(ps, Ysub)

  c(auroc_true = mt[1], auroc_naive = mn[1], auroc_subject = msu[1],
    brier_true = mt[2], brier_naive = mn[2], brier_subject = msu[2],
    uniq_stats(d$XF), prevalence = mean(d$Y))
}

## ---------------------------------------------------------------------------
## Design grid
## ---------------------------------------------------------------------------
grid <- expand.grid(composition = names(COMPOSITIONS), n = c(50, 100, 200),
                    R2 = c(0.05, 0.15, 0.30), learner = c("logistic", "xgboost"),
                    stringsAsFactors = FALSE)
grid$condition_id <- seq_len(nrow(grid))
cat("Condition-by-learner cells:", nrow(grid), " Iterations each:", N_ITER, "\n")

plan(multisession, workers = N_WORKERS)
ck <- file.path(OUT_DIR, "checkpoint.rds")
res_file <- file.path(OUT_DIR, "categorical_results.csv")
start <- 1
if (RESUME && file.exists(ck)) { start <- readRDS(ck) + 1; cat("Resuming at", start, "\n") }

t0 <- Sys.time()
for (i in seq(start, nrow(grid))) {
  g <- grid[i, ]; ty <- COMPOSITIONS[[g$composition]]
  M <- future_map_dfr(seq_len(N_ITER), function(it) {
    set.seed(MASTER_SEED + 1000L * g$condition_id + it)
    as.data.frame(t(one_iter(ty, g$n, g$R2, g$learner)))
  }, .options = furrr_options(seed = TRUE))
  out <- cbind(g[rep(1, 1), c("condition_id","composition","n","R2","learner")],
               as.data.frame(t(colMeans(M, na.rm = TRUE))),
               n_valid = sum(!is.na(M[[1]])),
               mcse_auroc_naive = sd(M$auroc_naive - M$auroc_true, na.rm = TRUE) / sqrt(N_ITER))
  write.table(out, res_file, sep = ",", row.names = FALSE,
              col.names = !file.exists(res_file), append = file.exists(res_file))
  saveRDS(i, ck)
  if (i %% 5 == 0 || i == nrow(grid))
    cat(sprintf("  %d/%d  elapsed %.1f min\n", i, nrow(grid),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
plan(sequential)

## ---------------------------------------------------------------------------
## Summary against the pre-specified prediction
## ---------------------------------------------------------------------------
## metrics() returns a named vector, so c(auroc_true = mt[1]) produced column
## names of the form "auroc_true.auroc". The names are stripped here rather than
## at source, so that rows written before this fix remain append-compatible.
r <- read.csv(res_file)
names(r) <- sub("\\..*$", "", names(r))
r <- r |>
  mutate(d_auroc_naive = auroc_naive - auroc_true,
         d_auroc_subj  = auroc_subject - auroc_true,
         d_brier_naive = brier_true - brier_naive)

cat("\n=== Equivalence-class structure by composition (n = 50, 100, 200) ===\n")
print(r |> filter(learner == "logistic") |> group_by(composition, n) |>
        summarise(u_raw = round(mean(u_raw), 3), u_q4 = round(mean(u_q4), 3),
                  n_classes = round(mean(n_classes), 1),
                  max_class = round(mean(max_class), 1), .groups = "drop") |>
        as.data.frame())

cat("\n=== Naive optimism by composition and learner (prediction: declines) ===\n")
print(r |> group_by(learner, composition) |>
        summarise(auroc_naive = round(mean(auroc_naive), 3),
                  d_auroc_naive = round(mean(d_auroc_naive), 4),
                  d_brier_naive = round(mean(d_brier_naive), 4),
                  d_auroc_subj = round(mean(d_auroc_subj), 4), .groups = "drop") |>
        as.data.frame())

cat("\n=== Ceiling check: proportion of cells with naive AUROC >= 0.999 ===\n")
print(r |> group_by(learner, composition) |>
        summarise(pct_ceiling = round(mean(auroc_naive >= 0.999), 3), .groups = "drop") |>
        as.data.frame())

## ---------------------------------------------------------------------------
## SELF-VALIDATION: the "5 continuous" composition at p_F = 5 reproduces the
## main factorial's Arm A conditions at the same predictor count, so those nine
## cells per learner should match the published values.
## ---------------------------------------------------------------------------
pub <- data.frame(
  learner = c(rep("logistic", 9), rep("xgboost", 9)),
  R2 = rep(c(0.05, 0.15, 0.30), 6)[1:18],
  n  = rep(c(50, 50, 50, 100, 100, 100, 200, 200, 200), 2),
  pub_naive = c(0.1424, 0.1078, 0.0756, 0.0832, 0.0566, 0.0430, 0.0466, 0.0295, 0.0195,
                0.4718, 0.4230, 0.3612, 0.4656, 0.4042, 0.3303, 0.4568, 0.3864, 0.3062))
cat("\n=== SELF-VALIDATION: '5 continuous' against published v8 Arm A, p_F = 5 ===\n")
chk <- r |> filter(composition == "5 continuous") |>
  select(learner, n, R2, obs_naive = d_auroc_naive, mcse = mcse_auroc_naive) |>
  left_join(pub, by = c("learner", "n", "R2")) |>
  mutate(z = round((obs_naive - pub_naive) / mcse, 2),
         across(c(obs_naive, pub_naive), \(x) round(x, 4)))
print(as.data.frame(chk))
cat("\nMean z:", round(mean(chk$z, na.rm = TRUE), 3),
    " Max |z|:", round(max(abs(chk$z), na.rm = TRUE), 2), "\n")

cat("\n=== PRIMARY TEST: does leakage depend on subject uniqueness? ===\n")
print(r |> group_by(learner, composition) |>
        summarise(u_raw = round(mean(u_raw), 3),
                  n_classes = round(mean(n_classes), 1),
                  d_auroc_naive = round(mean(d_auroc_naive), 4),
                  d_brier_naive = round(mean(d_brier_naive), 4), .groups = "drop") |>
        as.data.frame())
for (lr in unique(r$learner)) {
  s <- r[r$learner == lr, ]
  cat(sprintf("  %s: cor(u_raw, naive optimism) = %.3f (p = %.3g)\n", lr,
              cor(s$u_raw, s$d_auroc_naive), cor.test(s$u_raw, s$d_auroc_naive)$p.value))
}

cat("\nWritten to:", res_file, "\n")
cat("Total runtime:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

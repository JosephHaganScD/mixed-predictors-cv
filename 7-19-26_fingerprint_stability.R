###############################################################################
# Fingerprint stability diagnostic for naive cross-validation
#
# Purpose: the existing fold-specific recognizability statistic (u_L) measures
# whether a subject occupies a singleton bin within a training fold. It is
# invariant to ICC and rho because quantile binning is scale-invariant, so it
# cannot explain the ICC and rho gradients in naive-CV optimism (Section 3.4).
#
# This script computes a complementary STABILITY statistic: whether a subject's
# fingerprint computed from their training-fold rows still matches the same
# subject when computed from their held-out rows. Uniqueness without stability
# does not produce recognition.
#
# Predictions to be tested:
#   (1) stability = 1.000 exactly at ICC = 1 (fixed covariates), by construction
#   (2) stability increases monotonically with ICC
#   (3) stability increases with rho (held-out rows are interleaved in time
#       among training rows, so temporal smoothness aids matching)
#   (4) stability, unlike u_L, tracks the optimism gradients in Arm B
#
# No model fitting is performed. Runtime is a few minutes, sequential.
# Deliberately NOT parallelized to avoid core contention with any concurrent
# simulation run.
#
# Author: J. Hagan   Date: 2026-07-19
###############################################################################

library(dplyr)
library(tidyr)

# Set OUT_DIR to the folder where results should be saved.
OUT_DIR <- "."   # change to your preferred output directory
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

MASTER_SEED <- 20260719
N_ITER      <- 200      # proportions averaged over n subjects x 10 folds; ample
K_FOLDS     <- 10
T_CONST     <- 30       # matches main factorial

## ---------------------------------------------------------------------------
## Data generation: matches the v8 longitudinal generator exactly.
##   x_ilt = b_il + e_ilt,  b_il ~ N(0, ICC),
##   e_il1 ~ N(0, 1 - ICC), e_ilt = rho * e_il,t-1 + eps,
##   var(eps) = (1 - ICC)(1 - rho^2)   ->  marginal var = 1
## At ICC = 1 the within-subject error vanishes and the predictor is constant
## within subject; rho is then irrelevant (undefined, not 1), which is the
## fixed-covariate limiting case.
## ---------------------------------------------------------------------------
gen_long <- function(n, p_L, T_const, ICC, rho) {
  b <- matrix(rnorm(n * p_L, 0, sqrt(ICC)), nrow = n, ncol = p_L)
  X <- array(0, dim = c(n, T_const, p_L))
  if (ICC < 1) {
    s2_eps <- (1 - ICC) * (1 - rho^2)
    for (l in seq_len(p_L)) {
      e <- matrix(0, nrow = n, ncol = T_const)
      e[, 1] <- rnorm(n, 0, sqrt(1 - ICC))               # stationary start
      for (t in 2:T_const) e[, t] <- rho * e[, t - 1] + rnorm(n, 0, sqrt(s2_eps))
      X[, , l] <- b[, l] + e
    }
  } else {
    for (l in seq_len(p_L)) X[, , l] <- b[, l]           # constant within subject
  }
  # long format: one row per (subject, time)
  data.frame(
    subj = rep(seq_len(n), each = T_const),
    time = rep(seq_len(T_const), times = n)
  ) |>
    cbind(matrix(aperm(X, c(2, 1, 3)), nrow = n * T_const, ncol = p_L,
                 dimnames = list(NULL, paste0("x", seq_len(p_L)))))
}

## ---------------------------------------------------------------------------
## Bin a matrix of subject-level fingerprints into equivalence classes using
## quantile cut points derived from the TRAINING-fold subject means, so that
## held-out means are binned on the same scale the learner would have seen.
## ---------------------------------------------------------------------------
bin_codes <- function(train_mat, test_mat, n_bins) {
  codes_tr <- matrix("", nrow = nrow(train_mat), ncol = ncol(train_mat))
  codes_te <- matrix("", nrow = nrow(test_mat),  ncol = ncol(test_mat))
  for (j in seq_len(ncol(train_mat))) {
    br <- unique(quantile(train_mat[, j], probs = seq(0, 1, length.out = n_bins + 1),
                          na.rm = TRUE))
    br[1] <- -Inf; br[length(br)] <- Inf
    codes_tr[, j] <- as.character(cut(train_mat[, j], breaks = br, labels = FALSE,
                                      include.lowest = TRUE))
    codes_te[, j] <- as.character(cut(test_mat[, j],  breaks = br, labels = FALSE,
                                      include.lowest = TRUE))
  }
  list(train = apply(codes_tr, 1, paste, collapse = "|"),
       test  = apply(codes_te, 1, paste, collapse = "|"))
}

## ---------------------------------------------------------------------------
## Core: one simulated dataset, naive K-fold at the ROW level, then for each
## fold compute uniqueness and stability of the subject fingerprint.
## ---------------------------------------------------------------------------
one_iteration <- function(n, p_L, T_const, ICC, rho, K = K_FOLDS) {
  dat  <- gen_long(n, p_L, T_const, ICC, rho)
  xcol <- paste0("x", seq_len(p_L))
  dat$fold <- sample(rep_len(seq_len(K), nrow(dat)))    # naive: rows, not subjects

  res <- matrix(NA_real_, nrow = K, ncol = 6)
  colnames(res) <- c("u_q4", "u_q10", "stab_q4", "stab_q10",
                     "nn_subj", "nn_row")

  for (k in seq_len(K)) {
    tr <- dat[dat$fold != k, , drop = FALSE]
    te <- dat[dat$fold == k, , drop = FALSE]

    # subjects present in both parts of this fold (nearly all under naive CV)
    keep <- intersect(unique(tr$subj), unique(te$subj))
    if (length(keep) < 3) next
    tr <- tr[tr$subj %in% keep, , drop = FALSE]
    te <- te[te$subj %in% keep, , drop = FALSE]

    # fold-specific fingerprints: subject means over training rows and over
    # held-out rows, respectively
    m_tr <- as.matrix(aggregate(tr[xcol], by = list(subj = tr$subj), FUN = mean)[, -1, drop = FALSE])
    m_te <- as.matrix(aggregate(te[xcol], by = list(subj = te$subj), FUN = mean)[, -1, drop = FALSE])
    sid  <- sort(keep)

    for (nb in c(4, 10)) {
      bc <- bin_codes(m_tr, m_te, nb)
      # (a) UNIQUENESS: existing diagnostic, proportion of subjects in a
      #     singleton training-fold equivalence class
      u <- mean(table(bc$train)[bc$train] == 1)
      # (b) STABILITY: proportion of subjects whose held-out fingerprint falls
      #     in the SAME equivalence class as their own training fingerprint
      s <- mean(bc$test == bc$train)
      if (nb == 4)  { res[k, "u_q4"]  <- u; res[k, "stab_q4"]  <- s }
      if (nb == 10) { res[k, "u_q10"] <- u; res[k, "stab_q10"] <- s }
    }

    # (c) NEAREST-CENTROID MATCH RATE, subject level: is the training-fold
    #     centroid closest to a subject's held-out mean that subject's own?
    D <- as.matrix(dist(rbind(m_te, m_tr)))[seq_len(nrow(m_te)),
                                            nrow(m_te) + seq_len(nrow(m_tr)), drop = FALSE]
    res[k, "nn_subj"] <- mean(sid[apply(D, 1, which.min)] == sid)

    # (d) NEAREST-CENTROID MATCH RATE, row level: closest to what a partitioning
    #     learner actually does with an individual held-out observation
    Xte <- as.matrix(te[xcol])
    Dr  <- as.matrix(dist(rbind(Xte, m_tr)))[seq_len(nrow(Xte)),
                                             nrow(Xte) + seq_len(nrow(m_tr)), drop = FALSE]
    res[k, "nn_row"] <- mean(sid[apply(Dr, 1, which.min)] == te$subj)
  }
  colMeans(res, na.rm = TRUE)
}

## ---------------------------------------------------------------------------
## Design grid. ICC = 1 is the fixed-covariate limiting case and is included
## as a built-in anchor; rho is irrelevant there and is retained only so the
## grid stays rectangular (the two rho rows should agree to Monte Carlo error).
## ---------------------------------------------------------------------------
grid <- expand.grid(
  ICC = c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 1),
  rho = c(0.3, 0.7),
  p_L = c(2, 5),
  n   = c(50, 100, 200),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
grid$condition_id <- seq_len(nrow(grid))
cat("Conditions:", nrow(grid), " Iterations each:", N_ITER, "\n")

t0 <- Sys.time()
out <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  M <- matrix(NA_real_, nrow = N_ITER, ncol = 6)
  for (it in seq_len(N_ITER)) {
    set.seed(MASTER_SEED + 1000L * g$condition_id + it)   # reproducible, resumable
    M[it, ] <- one_iteration(g$n, g$p_L, T_CONST, g$ICC, g$rho)
  }
  colnames(M) <- c("u_q4", "u_q10", "stab_q4", "stab_q10", "nn_subj", "nn_row")
  out[[i]] <- cbind(
    g[rep(1, 1), c("condition_id", "n", "p_L", "ICC", "rho")],
    as.data.frame(t(colMeans(M, na.rm = TRUE))),
    as.data.frame(t(apply(M, 2, function(z) sd(z, na.rm = TRUE) / sqrt(sum(!is.na(z)))))) |>
      setNames(paste0("mcse_", colnames(M)))
  )
  if (i %% 10 == 0 || i == nrow(grid))
    cat(sprintf("  %d/%d  elapsed %.1f min\n", i, nrow(grid),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

res <- bind_rows(out)
write.csv(res, file.path(OUT_DIR, "fingerprint_stability_summary.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## Immediate checks against the four predictions
## ---------------------------------------------------------------------------
cat("\n=== Prediction 1: stability = 1.000 at ICC = 1 (fixed-covariate limit) ===\n")
print(res |> filter(ICC == 1) |>
        select(n, p_L, rho, stab_q4, stab_q10, nn_subj, nn_row) |>
        mutate(across(where(is.numeric), \(z) round(z, 4))))

cat("\n=== Predictions 2-3: stability by ICC and rho (uniqueness shown for contrast) ===\n")
print(res |> group_by(p_L, ICC, rho) |>
        summarise(u_q4 = mean(u_q4), stab_q4 = mean(stab_q4),
                  nn_subj = mean(nn_subj), nn_row = mean(nn_row), .groups = "drop") |>
        mutate(across(where(is.numeric), \(z) round(z, 4))) |>
        as.data.frame())

cat("\n=== Prediction 4: does stability track Arm B optimism where uniqueness did not? ===\n")
cat("Compare the ICC gradient below against Arm B naive optimism\n")
cat("  (XGBoost AUROC: 0.080 / 0.143 / 0.209 at ICC 0.3 / 0.7 / 0.9;\n")
cat("   XGBoost Brier: 0.017 / 0.045 / 0.084 at the same levels).\n")
print(res |> filter(ICC %in% c(0.3, 0.7, 0.9)) |> group_by(ICC) |>
        summarise(u_q4 = mean(u_q4), stab_q4 = mean(stab_q4),
                  nn_subj = mean(nn_subj), nn_row = mean(nn_row), .groups = "drop") |>
        mutate(across(where(is.numeric), \(z) round(z, 4))) |>
        as.data.frame())

cat("\nWritten to:", file.path(OUT_DIR, "fingerprint_stability_summary.csv"), "\n")
cat("Total runtime:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

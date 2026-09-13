###############################################################################
# Mixed-fingerprint stability diagnostic  (manuscript Sections 3.10 and 4.5)
#
# The 19 July stability run covered longitudinal fingerprints and the
# fixed-covariate limit, but not JOINT fixed-and-longitudinal fingerprints,
# which is the configuration the manuscript is actually about. This script
# closes that gap and is intended to strengthen Section 4.5, where the
# redundancy evidence is currently thinnest.
#
# THE ARGUMENT THIS TESTS
# Let s_F, s_L, s_joint denote the proportion of subjects whose held-out
# fingerprint has that same subject's training fingerprint as nearest neighbour,
# computed from fixed predictors only, longitudinal predictors only, and both
# combined. For CONTINUOUS fixed covariates, s_F = 1.000 exactly and by
# construction, since a fixed covariate is identical across all of a subject's
# rows. Because s_joint cannot exceed 1.000, the longitudinal channel has zero
# marginal identification capacity to contribute once continuous fixed
# covariates are present. That is the redundancy claim of Section 4.5 stated as
# a model-free property of the design matrix, requiring no fitted model, no
# AUROC, and no cross-arm comparison.
#
# The claim is clean here precisely BECAUSE s_F = 1.000 caps s_joint. It would
# not be clean if s_F were below one, since two individually non-discriminating
# fingerprints can jointly discriminate (the complementary disambiguation
# problem of Section 2.6). Categorical fixed compositions, for which s_F < 1,
# are therefore included so that this boundary is characterised rather than
# assumed, and they connect this analysis to the categorical arm of Section 3.11.
#
# PRE-SPECIFIED PREDICTIONS
#  (1) s_F = 1.000 exactly for all-continuous fixed sets, at every n.
#  (2) s_joint = 1.000, or very close, whenever at least two continuous fixed
#      covariates are present, irrespective of ICC, rho, and p_L.
#  (3) The marginal increment s_joint - s_F is approximately zero for continuous
#      fixed sets: the longitudinal channel adds no identification capacity.
#  (4) Adding noisy longitudinal dimensions may REDUCE s_joint below s_F when
#      p_L is large relative to p_F, since within-subject noise contributes to
#      the joint distance. If observed, this is a genuine finding and not a
#      failure: it would mean a mixed fingerprint can be less stable than the
#      fixed fingerprint alone, which bounds redundancy from the other side.
#  (5) For categorical fixed sets, s_F < 1 and s_joint may exceed both s_F and
#      s_L through complementary disambiguation, in which case redundancy is
#      partial rather than complete and Section 4.5 must say so.
#
# PERFORMANCE NOTE
# The 19 July script used as.matrix(dist(rbind(Xte, m_tr))), which computes the
# full (n_te + n_tr)^2 distance matrix when only the n_te x n_tr cross-block is
# needed. At n = 200 that is roughly 38 million entries per fold against 1.2
# million required. This script computes cross-distances directly from the
# identity ||a - b||^2 = ||a||^2 + ||b||^2 - 2 a.b, which is the reason the
# expected runtime here is far shorter than the 115 minutes that run took.
#
# No model fitting is performed.
#
# Author: J. Hagan   Date: 2026-07-19
###############################################################################

library(dplyr)

# Set OUT_DIR to the folder where results should be saved.
OUT_DIR <- "."   # change to your preferred output directory
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

MASTER_SEED <- 20260719
N_ITER      <- 200
K_FOLDS     <- 10
T_CONST     <- 30
ROW_SUBSAMPLE <- 300   # held-out rows sampled per fold for nn_row; NA = use all

## ---------------------------------------------------------------------------
## Fixed-predictor compositions. "cont" is continuous; "bin" is Bernoulli(0.5)
## standardised to unit variance, which maximises the number of distinct
## equivalence classes and is therefore the conservative choice.
## ---------------------------------------------------------------------------
FIXED_COMP <- list(
  "none"                    = character(0),
  "1 continuous"            = c("cont"),
  "2 continuous"            = c("cont", "cont"),
  "5 continuous"            = rep("cont", 5),
  "2 continuous + 3 binary" = c("cont", "cont", "bin", "bin", "bin"),
  "5 binary"                = rep("bin", 5)
)

gen_fixed <- function(n, types, p_bin = 0.5) {
  if (length(types) == 0) return(matrix(0, n, 0))
  X <- matrix(0, n, length(types))
  for (j in seq_along(types))
    X[, j] <- if (types[j] == "cont") rnorm(n) else
      (rbinom(n, 1, p_bin) - p_bin) / sqrt(p_bin * (1 - p_bin))
  X
}

gen_long <- function(n, p_L, T_const, ICC, rho) {
  if (p_L == 0) return(array(0, dim = c(n, T_const, 0)))
  X <- array(0, dim = c(n, T_const, p_L))
  s2 <- (1 - ICC) * (1 - rho^2)
  for (l in seq_len(p_L)) {
    b <- rnorm(n, 0, sqrt(ICC))
    e <- matrix(0, n, T_const); e[, 1] <- rnorm(n, 0, sqrt(1 - ICC))
    for (t in 2:T_const) e[, t] <- rho * e[, t - 1] + rnorm(n, 0, sqrt(s2))
    X[, , l] <- b + e
  }
  X
}

## ---------------------------------------------------------------------------
## Cross-distance nearest-neighbour match rate.
## A: query rows (held-out fingerprints); B: reference rows (training centroids,
## one per subject, ordered by sid). Returns the proportion of queries whose
## nearest reference is the subject given in truth_idx (indices into sid).
## ---------------------------------------------------------------------------
nn_match <- function(A, B, truth_idx) {
  if (ncol(A) == 0) return(NA_real_)                 # no fingerprint to match on
  cross <- A %*% t(B)
  d2 <- outer(rowSums(A^2), rowSums(B^2), "+") - 2 * cross
  mean(max.col(-d2, ties.method = "first") == truth_idx)
}

one_iteration <- function(n, fixed_types, p_L, ICC, rho) {
  p_F <- length(fixed_types)
  XF  <- gen_fixed(n, fixed_types)
  XL  <- gen_long(n, p_L, T_CONST, ICC, rho)
  idx <- rep(seq_len(n), each = T_CONST)
  XLrow <- if (p_L > 0) matrix(aperm(XL, c(2, 1, 3)), nrow = n * T_CONST, ncol = p_L) else
    matrix(0, n * T_CONST, 0)
  XFrow <- if (p_F > 0) XF[idx, , drop = FALSE] else matrix(0, n * T_CONST, 0)

  fold <- sample(rep_len(seq_len(K_FOLDS), n * T_CONST))
  res  <- matrix(NA_real_, K_FOLDS, 6)
  colnames(res) <- c("s_F", "s_L", "s_joint", "nnrow_F", "nnrow_L", "nnrow_joint")

  for (k in seq_len(K_FOLDS)) {
    tr <- fold != k
    sid <- sort(intersect(unique(idx[tr]), unique(idx[!tr])))
    if (length(sid) < 3) next
    ktr <- tr & (idx %in% sid); kte <- (!tr) & (idx %in% sid)

    # training centroids and held-out centroids, one row per subject in sid
    cen <- function(M, sel) {
      if (ncol(M) == 0) return(matrix(0, length(sid), 0))
      do.call(rbind, lapply(split(seq_len(nrow(M))[sel], idx[sel]),
                            function(r) colMeans(M[r, , drop = FALSE])))
    }
    trF <- cen(XFrow, ktr); teF <- cen(XFrow, kte)
    trL <- cen(XLrow, ktr); teL <- cen(XLrow, kte)
    truth <- seq_along(sid)

    res[k, "s_F"]     <- nn_match(teF, trF, truth)
    res[k, "s_L"]     <- nn_match(teL, trL, truth)
    res[k, "s_joint"] <- nn_match(cbind(teF, teL), cbind(trF, trL), truth)

    # row-level match rate, optionally on a subsample of held-out rows
    ridx <- which(kte)
    if (!is.na(ROW_SUBSAMPLE) && length(ridx) > ROW_SUBSAMPLE)
      ridx <- sample(ridx, ROW_SUBSAMPLE)
    rt <- match(idx[ridx], sid)
    res[k, "nnrow_F"]     <- nn_match(XFrow[ridx, , drop = FALSE], trF, rt)
    res[k, "nnrow_L"]     <- nn_match(XLrow[ridx, , drop = FALSE], trL, rt)
    res[k, "nnrow_joint"] <- nn_match(cbind(XFrow[ridx, , drop = FALSE],
                                            XLrow[ridx, , drop = FALSE]),
                                      cbind(trF, trL), rt)
  }
  colMeans(res, na.rm = TRUE)
}

## ---------------------------------------------------------------------------
## Design grid. ICC and rho are crossed at the main-factorial levels; the
## "none" fixed composition reproduces the pure-longitudinal reference from the
## 19 July run and serves as a consistency check against it.
## ---------------------------------------------------------------------------
grid <- expand.grid(fixed_comp = names(FIXED_COMP), p_L = c(2, 5),
                    ICC = c(0.3, 0.7, 0.9), rho = c(0.3, 0.7), n = c(50, 100, 200),
                    stringsAsFactors = FALSE)
grid$condition_id <- seq_len(nrow(grid))
cat("Conditions:", nrow(grid), " Iterations each:", N_ITER, "\n")

t0 <- Sys.time(); out <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]; ty <- FIXED_COMP[[g$fixed_comp]]
  M <- matrix(NA_real_, N_ITER, 6)
  for (it in seq_len(N_ITER)) {
    set.seed(MASTER_SEED + 1000L * g$condition_id + it)
    M[it, ] <- one_iteration(g$n, ty, g$p_L, g$ICC, g$rho)
  }
  colnames(M) <- c("s_F", "s_L", "s_joint", "nnrow_F", "nnrow_L", "nnrow_joint")
  out[[i]] <- cbind(g[rep(1, 1), c("condition_id","fixed_comp","p_L","ICC","rho","n")],
                    p_F = length(ty),
                    as.data.frame(t(colMeans(M, na.rm = TRUE))),
                    as.data.frame(t(apply(M, 2, function(z)
                      sd(z, na.rm = TRUE) / sqrt(sum(!is.na(z)))))) |>
                      setNames(paste0("mcse_", colnames(M))))
  if (i %% 20 == 0 || i == nrow(grid))
    cat(sprintf("  %d/%d  elapsed %.1f min\n", i, nrow(grid),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- bind_rows(out)
res$increment <- res$s_joint - pmax(res$s_F, res$s_L, na.rm = TRUE)
write.csv(res, file.path(OUT_DIR, "mixed_fingerprint_stability.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## Checks against the pre-specified predictions
## ---------------------------------------------------------------------------
cat("\n=== Predictions 1-2: s_F and s_joint by fixed composition ===\n")
print(res |> group_by(fixed_comp, p_F) |>
        summarise(s_F = mean(s_F, na.rm = TRUE), s_L = mean(s_L),
                  s_joint = mean(s_joint), .groups = "drop") |>
        mutate(across(where(is.numeric), \(z) round(z, 4))) |> as.data.frame())

cat("\n=== Prediction 3-4: marginal increment of the longitudinal channel ===\n")
cat("increment = s_joint - max(s_F, s_L). Near zero means no marginal capacity;\n")
cat("negative means noisy longitudinal dimensions degrade the fixed fingerprint.\n")
print(res |> group_by(fixed_comp, p_L) |>
        summarise(s_F = round(mean(s_F, na.rm = TRUE), 4),
                  s_joint = round(mean(s_joint), 4),
                  increment = round(mean(increment), 4), .groups = "drop") |>
        as.data.frame())

cat("\n=== Prediction 2: is s_joint invariant to ICC and rho once fixed present? ===\n")
print(res |> group_by(fixed_comp, ICC, rho) |>
        summarise(s_joint = round(mean(s_joint), 4), .groups = "drop") |>
        as.data.frame())

cat("\n=== Prediction 5: complementary disambiguation for categorical fixed sets ===\n")
print(res |> filter(fixed_comp %in% c("5 binary", "2 continuous + 3 binary")) |>
        group_by(fixed_comp, n, p_L) |>
        summarise(s_F = round(mean(s_F), 4), s_L = round(mean(s_L), 4),
                  s_joint = round(mean(s_joint), 4),
                  exceeds_both = round(mean(s_joint > pmax(s_F, s_L) + 0.01), 3),
                  .groups = "drop") |> as.data.frame())

cat("\n=== Consistency check against the 19 July run (fixed_comp = none) ===\n")
cat("These should reproduce nn_subj and nn_row from fingerprint_stability_summary.csv\n")
print(res |> filter(fixed_comp == "none") |> group_by(p_L, ICC, rho) |>
        summarise(s_L = round(mean(s_L), 4), nnrow_L = round(mean(nnrow_L), 4),
                  .groups = "drop") |> as.data.frame())

cat("\nWritten to:", file.path(OUT_DIR, "mixed_fingerprint_stability.csv"), "\n")
cat("Total runtime:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

###############################################################################
# run_supplemental.R — Driver Script, Combined Supplemental Analysis (v2)
# Mixed-Predictor CV Optimism Decomposition Simulation
#
# Manuscript Section 3.8. Addresses two gaps not covered by the main
# 567-condition factorial or the 48-condition paired sub-study
# (run_paired_substudy.R):
#
#   1. Arm C predictor-count coverage. Neither the main factorial nor the
#      paired sub-study evaluates the "crossed extremes" (p_F, p_L) =
#      (2,5) and (10,2), maximally imbalanced in opposite directions.
#
#   2. Outcome prevalence. Both the main factorial and the paired sub-study
#      hold prevalence fixed at 0.480 throughout, and cannot speak to
#      whether the additivity departure depends on prevalence, which the
#      empirical illustration suggests it does (more pronounced departure
#      for the low-prevalence secondary outcome, severe ROP).
#
# CHANGE FROM v1 (patched 2026-07-10): because this analysis makes a
# quantitative additivity claim, using the same criterion applied in the
# paired sub-study, it now uses the same within-iteration paired design
# (simulate_paired.R::run_one_paired_iteration()) rather than fitting only
# the mixed-predictor model per iteration and comparing against separately
# generated Arm A/B conditions. v1's design would have reintroduced exactly
# the unpaired-comparison flaw identified for the main factorial (manuscript
# Section 2.2), just relocated to this supplemental analysis.
#
# Grid (48 conditions, manuscript Section 3.8):
#   (p_F, p_L)         in {(2,5), (10,2)}                      [2 levels]
#   target_prevalence  in {0.480, 0.150}                       [2 levels]
#   ICC                = 0.7 (fixed; clearest ICC-by-rho interaction in
#                        Arm B, manuscript Section 3.2)
#   rho                in {0.3, 0.7}                           [2 levels]
#   R2_total           in {0.05, 0.15}                         [2 levels;
#                        swept rather than fixed, since neither value has a
#                        clean independent justification for this specific
#                        grid -- see build_supplemental_grid() below]
#   n                  in {50, 100, 200}                       [3 levels]
#   signal_allocation  = "balanced" (fixed)
#   Total: 2 x 2 x 2 x 2 x 3 = 48 conditions x 500 iterations x 2 learners.
#
# This is a SEPARATE script with a SEPARATE checkpoint/output directory from
# run_simulation.R and run_paired_substudy.R, run after the main factorial's
# checkpoint is stable (not concurrently -- see run_simulation.R header for
# the core-oversubscription reasoning).
#
# Author: Joseph L. Hagan, ScD, MSPH
# Version: 3 (2026-07-10)
###############################################################################

# =============================================================================
# USER SETTINGS
# =============================================================================

PILOT         <- FALSE       # TRUE = quick sanity check; FALSE = full run
N_ITER        <- 500L
N_ITER_PILOT  <- 10L
MASTER_SEED   <- 20260710L   # distinct from main factorial (20260505L) and
                              # from run_paired_substudy.R (20260709L)
N_CORES       <- parallel::detectCores() - 1L
T_CONST       <- 30L
N_TEST        <- 5000L
K_FOLDS       <- 10L
RESUME        <- TRUE

ICC_FIXED <- 0.7   # fixed for this grid; see header rationale

# Same local, non-OneDrive-synced root as the main v7 run, sibling folder
# Set LOCAL_ROOT to the folder containing the simulation output CSVs.
# If running from the cloned repository, use "." for current directory.
LOCAL_ROOT <- "."   # change to your preferred data directory
OUT_DIR    <- file.path(LOCAL_ROOT,
                         if (PILOT) "results_supplemental_v3_pilot"
                         else "results_supplemental_v3")

script_dir <- LOCAL_ROOT
setwd(LOCAL_ROOT)

# =============================================================================
# SETUP
# =============================================================================

cat(sprintf("Combined Supplemental Analysis (v2, paired design)\nStarted: %s\n\n",
            format(Sys.time())))

suppressPackageStartupMessages({
  library(glmnet)
  library(xgboost)
  library(splines)
  library(pROC)
  library(future)
  library(furrr)
})

source(file.path(script_dir, "io_retry.R"))
source(file.path(script_dir, "dgm.R"))
source(file.path(script_dir, "learners.R"))
source(file.path(script_dir, "cv_strategies.R"))
source(file.path(script_dir, "metrics.R"))
source(file.path(script_dir, "simulate_paired.R"))

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

n_iter_use <- if (PILOT) N_ITER_PILOT else N_ITER
cat(sprintf("Mode: %s | Iterations: %d | Cores: %d\n\n",
            if (PILOT) "PILOT" else "FULL", n_iter_use, N_CORES))

# =============================================================================
# CONDITION GRID
# =============================================================================

build_supplemental_grid <- function() {

  pF_pL_pairs <- data.frame(p_F = c(2L, 10L), p_L = c(5L, 2L))  # crossed extremes

  base <- expand.grid(
    pF_pL_idx         = seq_len(nrow(pF_pL_pairs)),
    target_prevalence = c(0.480, 0.150),
    rho               = c(0.3, 0.7),
    R2_total          = c(0.05, 0.15),  # patched 2026-07-10: was fixed at
                                         # 0.15; now swept, since neither
                                         # value has a clean, non-arbitrary
                                         # justification on its own for this
                                         # grid specifically -- 0.05 makes
                                         # the imbalance effect easiest to
                                         # see but risks XGBoost's naive
                                         # AUROC saturating at 1.000 for
                                         # small n crossed with the low
                                         # (0.150) prevalence condition,
                                         # while 0.15 avoids that risk but
                                         # is otherwise no better motivated
                                         # than 0.05. Sweeping both removes
                                         # the anchor-point vulnerability
                                         # entirely and lets the ceiling-risk
                                         # concern show up in the data
                                         # directly rather than being argued
                                         # around.
    n                 = c(50L, 100L, 200L),
    stringsAsFactors  = FALSE)

  grid <- cbind(
    p_F               = pF_pL_pairs$p_F[base$pF_pL_idx],
    p_L               = pF_pL_pairs$p_L[base$pF_pL_idx],
    ICC               = ICC_FIXED,
    signal_allocation = "balanced",
    base[, setdiff(names(base), "pF_pL_idx")],
    stringsAsFactors  = FALSE)

  grid$condition_id <- seq_len(nrow(grid))
  id_cols <- c("condition_id", "n", "p_F", "p_L", "ICC", "rho",
               "R2_total", "signal_allocation", "target_prevalence")
  grid[, id_cols]
}

condition_grid <- build_supplemental_grid()
cat(sprintf("Supplemental condition grid: %d conditions (expect 48)\n", nrow(condition_grid)))
cat(sprintf("  (p_F,p_L) pairs: %s\n",
            paste(unique(paste0("(", condition_grid$p_F, ",", condition_grid$p_L, ")")),
                  collapse = ", ")))
cat(sprintf("  target_prevalence levels: %s\n",
            paste(unique(condition_grid$target_prevalence), collapse = ", ")))
cat(sprintf("  R2_total levels: %s\n\n",
            paste(unique(condition_grid$R2_total), collapse = ", ")))

# =============================================================================
# CHECKPOINT
# =============================================================================

progress_file <- file.path(OUT_DIR, "sim_progress_supplemental_v3.rds")
results_file  <- file.path(OUT_DIR, "sim_results_supplemental_v3.csv")
failures_file <- file.path(OUT_DIR, "sim_failures_supplemental_v3.csv")

if (RESUME && file.exists(progress_file)) {
  last_completed <- readRDS(progress_file)
  start_cond     <- last_completed + 1L
  cat(sprintf("Resuming from condition %d (last completed: %d)\n\n",
              start_cond, last_completed))
} else {
  start_cond <- 1L
  if (!RESUME) {
    to_remove <- c(results_file, failures_file, progress_file)
    file.remove(to_remove[file.exists(to_remove)])
    cat("RESUME = FALSE: starting clean.\n\n")
  }
}

remaining <- condition_grid[condition_grid$condition_id >= start_cond, , drop = FALSE]
cat(sprintf("Conditions remaining: %d\n\n", nrow(remaining)))

# =============================================================================
# PARALLEL BACKEND
# =============================================================================

future::plan(future::multisession, workers = N_CORES)
cat(sprintf("Parallel backend: multisession, %d workers\n\n", N_CORES))

# =============================================================================
# MAIN LOOP
# =============================================================================

for (ci in seq_len(nrow(remaining))) {

  cond <- remaining[ci, , drop = FALSE]
  cid  <- cond$condition_id

  cat(sprintf("[%d/%d | ID %d | n=%d pF=%d pL=%d ICC=%s rho=%s R2=%.2f alloc=%s prev=%.3f]\n",
              ci, nrow(remaining), cid,
              cond$n, cond$p_F, cond$p_L,
              cond$ICC, cond$rho, cond$R2_total, cond$signal_allocation,
              cond$target_prevalence))

  t_start <- proc.time()["elapsed"]

  iter_results <- furrr::future_map_dfr(
    seq_len(n_iter_use),
    function(it) {
      run_one_paired_iteration(
        condition_row = cond, iter_id = it,
        master_seed   = MASTER_SEED,
        T_const = T_CONST, n_test = N_TEST, K = K_FOLDS
      )
    },
    .options = furrr::furrr_options(seed = TRUE)
  )

  elapsed <- proc.time()["elapsed"] - t_start

  is_error <- iter_results$learner == "ERROR"
  valid    <- iter_results[!is_error, , drop = FALSE]
  failures <- iter_results[is_error,  , drop = FALSE]

  write_table_retry(valid,
              file = results_file, sep = ",",
              col.names = !file.exists(results_file),
              row.names = FALSE, append = file.exists(results_file))

  if (nrow(failures) > 0L)
    write_table_retry(failures,
                file = failures_file, sep = ",",
                col.names = !file.exists(failures_file),
                row.names = FALSE, append = file.exists(failures_file))

  saveRDS_retry(cid, file = progress_file)

  cat(sprintf("  -> %d valid, %d failures | %.1f sec\n\n",
              nrow(valid), nrow(failures), elapsed))
}

# =============================================================================
# SUMMARY
# =============================================================================

last_progress <- if (file.exists(progress_file)) readRDS(progress_file) else 0L
if (last_progress < max(condition_grid$condition_id)) {
  stop(sprintf(
    "Incomplete run: last completed condition_id = %d, expected %d. ",
    last_progress, max(condition_grid$condition_id)),
    "Not computing summary. Re-run this script (RESUME = TRUE) to continue.")
}

cat("All conditions complete (verified against checkpoint). Computing summary...\n")
results_full <- read.csv(results_file, stringsAsFactors = FALSE)

missing_conditions <- setdiff(condition_grid$condition_id, unique(results_full$condition_id))
if (length(missing_conditions) > 0L) {
  stop(sprintf(
    "Checkpoint claims completion but %d condition(s) are absent from %s: %s",
    length(missing_conditions), results_file,
    paste(missing_conditions, collapse = ", ")))
}

id_cols_required <- c("condition_id", "iteration", "learner")
malformed <- !stats::complete.cases(results_full[, id_cols_required])
if (any(malformed)) {
  cat(sprintf("Warning: %d malformed row(s) dropped.\n", sum(malformed)))
  results_full <- results_full[!malformed, , drop = FALSE]
}
dup_key <- paste(results_full$condition_id, results_full$iteration,
                  results_full$learner, sep = "_")
is_dup  <- duplicated(dup_key)
if (any(is_dup)) {
  cat(sprintf("Warning: %d duplicate row(s) dropped, keeping first occurrence.\n",
              sum(is_dup)))
  results_full <- results_full[!is_dup, , drop = FALSE]
}
still_missing <- setdiff(condition_grid$condition_id, unique(results_full$condition_id))
if (length(still_missing) > 0L) {
  stop(sprintf(
    "After removing malformed/duplicate rows, %d condition(s) are no longer present: %s. ",
    length(still_missing), paste(head(still_missing, 10), collapse = ", ")),
    "Re-run this script (RESUME = TRUE) to regenerate them.")
}

summary_df <- do.call(rbind, lapply(
  split(results_full,
        interaction(results_full$condition_id, results_full$learner, drop = TRUE)),
  function(df) {
    data.frame(
      condition_id      = df$condition_id[1],
      learner           = df$learner[1],
      n = df$n[1], p_F = df$p_F[1], p_L = df$p_L[1],
      ICC = df$ICC[1], rho = df$rho[1],
      R2_total          = df$R2_total[1],
      signal_allocation = df$signal_allocation[1],
      target_prevalence = df$target_prevalence[1],
      n_valid                  = sum(!is.na(df$additivity_departure)),
      mean_additivity_departure = mean(df$additivity_departure, na.rm = TRUE),
      mcse_additivity_departure = sd(df$additivity_departure, na.rm = TRUE) /
                                  sqrt(sum(!is.na(df$additivity_departure))),
      mean_delta_naive_fixed   = mean(df$delta_naive_fixed,  na.rm = TRUE),
      mean_delta_naive_long    = mean(df$delta_naive_long,   na.rm = TRUE),
      mean_delta_naive_mixed   = mean(df$delta_naive_mixed,  na.rm = TRUE),
      mean_credit_correlation  = mean(df$credit_correlation, na.rm = TRUE),
      mean_credit_jaccard_high = mean(df$credit_jaccard_high, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))

write_csv_retry(summary_df,
          file = file.path(OUT_DIR, "sim_summary_supplemental_v3.csv"),
          row.names = FALSE)

cat(sprintf("\nSupplemental analysis complete: %s\n", format(Sys.time())))
cat(sprintf("Results:  %s\n", results_file))
cat(sprintf("Summary:  %s\n", file.path(OUT_DIR, "sim_summary_supplemental_v3.csv")))
cat(sprintf("Failures: %s\n", failures_file))

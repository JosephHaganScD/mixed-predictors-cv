###############################################################################
# run_simulation.R — Driver Script (v8)
# Mixed-Predictor CV Optimism Decomposition Simulation
#
# Changes from v6 (patched 2026-07-09):
#   - Corrected outcome prevalence bug: target_prevalence was hardcoded to an
#     incorrect 0.62 inside simulate.R, independent of any setting here. Now
#     added explicitly to the condition grid (fixed at 0.480 for every
#     condition in this main factorial) and read from condition_row inside
#     simulate.R. See dgm.R and simulate.R patch notes for detail.
#   - I/O made robust to transient write failures (e.g. the OneDrive-sync
#     collision that truncated the v6 run at condition 333/567): all
#     checkpoint writes now go through write_table_retry()/saveRDS_retry()
#     (io_retry.R), which retry with backoff and then stop() the script
#     outright on repeated failure rather than allowing silent continuation.
#   - Output moved off the OneDrive-synced folder entirely, to a local path
#     confirmed writable by the user (local path set in LOCAL_ROOT below).
#   - Output versioned as v7 (not v6) specifically so that a resumed run can
#     never silently mix conditions generated under the old 62% target with
#     conditions generated under the corrected 48.0% target: a new version
#     has no pre-existing progress file, so RESUME = TRUE still forces a
#     clean start by construction.
#   - Added a completeness assertion before the SUMMARY section, so the
#     script cannot print "complete" on a truncated results file again.
#   - IMPORTANT: run this via `Rscript run_simulation.R` from the command
#     line (or R CMD BATCH), not by pasting into an open interactive
#     console. If pasted, an uncaught error (including the stop() calls
#     above) will not halt subsequent pasted lines, which is what allowed
#     the v6 run's SUMMARY block to execute on incomplete data after the
#     write failure.
#
# Changes from v2 to v6 (unchanged, retained for history):
#   - n=500 removed from factorial. Rationale: identity-mediated leakage
#     is at ceiling at n=500 (near-perfect subject uniqueness under any
#     reasonable discretization); the interesting variation in optimism
#     occurs at n in {50, 100, 200}. Removes ~40% of compute cost.
#   - ridge_spline removed. Two-learner design (ridge_linear + xgboost)
#     represents the full flexibility gradient adequately; ridge_spline
#     dominated runtime at large n without adding qualitatively new findings.
#   - Ridge lambda now truly fixed (no cv.glmnet); see learners.R v3.
#   - 567 conditions x 500 iterations x 2 learners.
#
# CV strategies: naive 10-fold, subject-level 10-fold
# Learners: ridge_linear (fixed lambda), xgboost
#
# Author: Joseph L. Hagan, ScD, MSPH
# Version: 8 (2026-07-11)
#
# v8 note: output renamed from v7 to v8 (results_v7 -> results_v8, etc.)
# purely to force a clean restart. dgm.R, learners.R, cv_strategies.R, and
# simulate.R all changed after some v7 conditions had already completed
# (closed-form R2 calibration and the hyperparameter fix, both patched
# 2026-07-10), so any conditions already in results_v7 were generated under
# pre-fix code. A new version directory has no existing progress file, so
# RESUME = TRUE is forced to start clean by construction, the same
# reasoning as the v6 -> v7 bump for the prevalence correction. The old
# results_v7 folder is not touched by this script and can be archived or
# deleted once you've confirmed you no longer need it.
###############################################################################

# =============================================================================
# USER SETTINGS
# =============================================================================

PILOT             <- FALSE       # TRUE = pilot run; FALSE = full factorial
N_ITER            <- 500L       # iterations per condition (full run)
N_ITER_PILOT      <- 10L        # iterations per condition (pilot)
MASTER_SEED       <- 20260505L
N_CORES           <- parallel::detectCores() - 1L
T_CONST           <- 30L
N_TEST            <- 5000L
K_FOLDS           <- 10L
RESUME            <- TRUE      # FALSE = clean start; TRUE = resume from checkpoint
                                 # (safe either way under v8 --- see header note)
TARGET_PREVALENCE <- 0.480      # corrected 2026-07-09; was erroneously 0.62

# Local, non-OneDrive-synced output root (confirmed writable 2026-07-09)
# Set LOCAL_ROOT to the folder where simulation outputs should be written.
# This should be a local (non-cloud-synced) path for best I/O performance.
LOCAL_ROOT <- "."   # change to your preferred output directory
OUT_DIR           <- file.path(LOCAL_ROOT, if (PILOT) "results_v8_pilot" else "results_v8")

# All six source files (io_retry.R, dgm.R, learners.R, cv_strategies.R,
# metrics.R, simulate.R) must live directly in LOCAL_ROOT. script_dir is
# hardcoded here (not inferred from interactive()/sys.frame()), and setwd()
# targets LOCAL_ROOT rather than the old OneDrive folder, so nothing in this
# session touches the synced folder at all. This is the same fix already
# applied to run_paired_substudy.R and run_supplemental.R; it was
# inadvertently lost from this file when it was rebuilt for the v8 update
# (2026-07-10) starting from an earlier, pre-fix copy. If you see a
# "cannot change working directory" error referencing the old OneDrive path,
# it means you're running a stale copy of this file from before this note,
# not a problem with your folder setup.
script_dir <- LOCAL_ROOT
setwd(LOCAL_ROOT)

# =============================================================================
# SETUP
# =============================================================================

cat(sprintf("Mixed-Predictor CV Optimism Simulation (v8)\nStarted: %s\n\n",
            format(Sys.time())))

suppressPackageStartupMessages({
  library(glmnet)
  library(xgboost)   # replaces ranger
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
source(file.path(script_dir, "simulate.R"))

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

n_iter_use <- if (PILOT) N_ITER_PILOT else N_ITER
cat(sprintf("Mode: %s | Iterations: %d | Cores: %d | Target prevalence: %.3f\n\n",
            if (PILOT) "PILOT" else "FULL", n_iter_use, N_CORES, TARGET_PREVALENCE))

# =============================================================================
# CONDITION GRID (target_prevalence column added 2026-07-09; grid values
# and row order otherwise unchanged from v1-v6)
# =============================================================================

build_condition_grid <- function() {

  arm_A <- expand.grid(
    arm = "A_fixed", p_F = c(2L, 5L, 10L), p_L = 0L,
    ICC = NA_real_, rho = NA_real_,
    R2_total = c(0.05, 0.15, 0.30), signal_allocation = "NA",
    n = c(50L, 100L, 200L), target_prevalence = TARGET_PREVALENCE,
    stringsAsFactors = FALSE)  # n=500 removed

  arm_B <- expand.grid(
    arm = "B_long", p_F = 0L, p_L = c(2L, 5L),
    ICC = c(0.3, 0.7, 0.9), rho = c(0.3, 0.7),
    R2_total = c(0.05, 0.15, 0.30), signal_allocation = "NA",
    n = c(50L, 100L, 200L), target_prevalence = TARGET_PREVALENCE,
    stringsAsFactors = FALSE)  # n=500 removed; ICC=0.9 added

  pF_pL_pairs <- data.frame(p_F = c(2L, 5L, 5L, 10L),
                             p_L = c(2L, 2L, 5L,  5L))
  arm_C_base <- expand.grid(
    pF_pL_idx = seq_len(nrow(pF_pL_pairs)),
    ICC = c(0.3, 0.7), rho = c(0.3, 0.7),
    R2_total = c(0.05, 0.15, 0.30),
    signal_allocation = c("fixed_dominant", "long_dominant", "balanced"),
    n = c(50L, 100L, 200L), target_prevalence = TARGET_PREVALENCE,
    stringsAsFactors = FALSE)  # n=500 removed
  arm_C <- cbind(
    arm = "C_mixed",
    p_F = pF_pL_pairs$p_F[arm_C_base$pF_pL_idx],
    p_L = pF_pL_pairs$p_L[arm_C_base$pF_pL_idx],
    arm_C_base[, setdiff(names(arm_C_base), "pF_pL_idx")],
    stringsAsFactors = FALSE)

  grid <- rbind(arm_A, arm_B, arm_C)
  grid$condition_id <- seq_len(nrow(grid))
  id_cols <- c("condition_id", "arm", "n", "p_F", "p_L", "ICC", "rho",
               "R2_total", "signal_allocation", "target_prevalence")
  grid[, id_cols]
}

condition_grid <- build_condition_grid()
cat(sprintf("Condition grid: %d conditions\n", nrow(condition_grid)))
cat(sprintf("  Arm A: %d | Arm B: %d | Arm C: %d\n\n",
            sum(condition_grid$arm == "A_fixed"),
            sum(condition_grid$arm == "B_long"),
            sum(condition_grid$arm == "C_mixed")))

# =============================================================================
# CHECKPOINT
# =============================================================================

progress_file <- file.path(OUT_DIR, "sim_progress_v8.rds")
results_file  <- file.path(OUT_DIR, "sim_results_v8.csv")
failures_file <- file.path(OUT_DIR, "sim_failures_v8.csv")

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

  cat(sprintf("[%d/%d | ID %d | %s | n=%d pF=%d pL=%d ICC=%s rho=%s R2=%.2f alloc=%s]\n",
              ci, nrow(remaining), cid,
              cond$arm, cond$n, cond$p_F, cond$p_L,
              ifelse(is.na(cond$ICC), "NA", cond$ICC),
              ifelse(is.na(cond$rho), "NA", cond$rho),
              cond$R2_total, cond$signal_allocation))

  t_start <- proc.time()["elapsed"]

  iter_results <- furrr::future_map_dfr(
    seq_len(n_iter_use),
    function(it) {
      run_one_iteration(
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

# Completeness check (added 2026-07-09): verify every condition actually
# completed before claiming so. This is what the v6 run's console output
# failed to do -- it printed "All conditions complete" after the loop was
# aborted by a write failure partway through, because nothing checked.
last_progress <- if (file.exists(progress_file)) readRDS(progress_file) else 0L
if (last_progress < max(condition_grid$condition_id)) {
  stop(sprintf(
    "Incomplete run: last completed condition_id = %d, expected %d. ",
    last_progress, max(condition_grid$condition_id)),
    "Not computing summary. Re-run this script (RESUME = TRUE) to continue ",
    "from the checkpoint.")
}

cat("All conditions complete (verified against checkpoint). Computing summary...\n")
results_full <- read.csv(results_file, stringsAsFactors = FALSE)

# Second check: the results file itself should also cover every condition,
# independent of the checkpoint value.
missing_conditions <- setdiff(condition_grid$condition_id, unique(results_full$condition_id))
if (length(missing_conditions) > 0L) {
  stop(sprintf(
    "Checkpoint claims completion but %d condition(s) are absent from %s: %s",
    length(missing_conditions), results_file,
    paste(head(missing_conditions, 10), collapse = ", ")))
}

# Malformed-row check (added 2026-07-10): drop any row with a missing value
# in a key identifier column. This is the failure mode a hard process kill
# (power loss, forced shutdown) could produce mid-write, as opposed to a
# graceful interruption, which write_table_retry()/saveRDS_retry() already
# handle. A dropped row here just means that iteration is regenerated on the
# next resume (same seed, so identical output), not a repair attempt.
id_cols_required <- c("condition_id", "iteration", "learner")
malformed <- !stats::complete.cases(results_full[, id_cols_required])
if (any(malformed)) {
  cat(sprintf(
    "Warning: %d malformed row(s) found (missing identifier field(s)) and dropped: %s\n",
    sum(malformed), results_file))
  results_full <- results_full[!malformed, , drop = FALSE]
}

# Duplicate-row check (added 2026-07-10): the checkpoint write order is
# results-file-append THEN progress-file-save, so a crash between those two
# steps could cause a condition to be regenerated and appended a second time
# on resume. Because every iteration's random draws are seeded
# deterministically from master_seed + condition_id + iteration (simulate.R),
# a duplicate (condition_id, iteration, learner) key is a byte-identical
# duplicate, not a second independent draw, so it is safe to keep the first
# occurrence and drop the rest rather than needing to reconcile two different
# values.
dup_key <- paste(results_full$condition_id, results_full$iteration,
                  results_full$learner, sep = "_")
is_dup  <- duplicated(dup_key)
if (any(is_dup)) {
  cat(sprintf(
    "Warning: %d duplicate (condition_id, iteration, learner) row(s) found and dropped, keeping first occurrence: %s\n",
    sum(is_dup), results_file))
  results_full <- results_full[!is_dup, , drop = FALSE]
}

# Re-verify completeness after dropping malformed/duplicate rows, in case
# dropping malformed rows uncovered a condition that is now actually
# incomplete (was previously masked by a malformed duplicate).
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
      arm               = df$arm[1],
      learner           = df$learner[1],
      n = df$n[1], p_F = df$p_F[1], p_L = df$p_L[1],
      ICC = df$ICC[1], rho = df$rho[1],
      R2_total          = df$R2_total[1],
      signal_allocation = df$signal_allocation[1],
      target_prevalence = df$prevalence_target[1],
      n_valid           = sum(!is.na(df$auroc_true)),
      mean_auroc_true    = mean(df$auroc_true,    na.rm = TRUE),
      mean_auroc_naive   = mean(df$auroc_naive,   na.rm = TRUE),
      mean_auroc_subject = mean(df$auroc_subject, na.rm = TRUE),
      mean_delta_naive   = mean(df$auroc_naive   - df$auroc_true, na.rm = TRUE),
      mean_delta_subject = mean(df$auroc_subject - df$auroc_true, na.rm = TRUE),
      mcse_delta_naive   = sd(df$auroc_naive   - df$auroc_true, na.rm = TRUE) /
                           sqrt(sum(!is.na(df$auroc_naive))),
      mcse_delta_subject = sd(df$auroc_subject - df$auroc_true, na.rm = TRUE) /
                           sqrt(sum(!is.na(df$auroc_subject))),
      mean_pct_unique_q4  = mean(df$pct_unique_q4,  na.rm = TRUE),
      mean_pct_unique_q10 = mean(df$pct_unique_q10, na.rm = TRUE),
      mean_prevalence_actual = mean(df$prevalence_actual, na.rm = TRUE),
      # Fold-specific mechanistic diagnostic (mean across iterations)
      mean_fold_recog_fixed_q4  = mean(df$fold_recog_fixed_q4,  na.rm = TRUE),
      mean_fold_recog_fixed_q10 = mean(df$fold_recog_fixed_q10, na.rm = TRUE),
      mean_fold_recog_long_q4   = mean(df$fold_recog_long_q4,   na.rm = TRUE),
      mean_fold_recog_long_q10  = mean(df$fold_recog_long_q10,  na.rm = TRUE),
      mean_fold_recog_joint_q4  = mean(df$fold_recog_joint_q4,  na.rm = TRUE),
      mean_fold_recog_joint_q10 = mean(df$fold_recog_joint_q10, na.rm = TRUE),
      mean_fold_recog_overlap_q4      = mean(df$fold_recog_overlap_q4,      na.rm = TRUE),
      mean_fold_recog_overlap_q10     = mean(df$fold_recog_overlap_q10,     na.rm = TRUE),
      mean_fold_recog_union_q4        = mean(df$fold_recog_union_q4,        na.rm = TRUE),
      mean_fold_recog_union_q10       = mean(df$fold_recog_union_q10,       na.rm = TRUE),
      mean_fold_recog_newly_joint_q4  = mean(df$fold_recog_newly_joint_q4,  na.rm = TRUE),
      mean_fold_recog_newly_joint_q10 = mean(df$fold_recog_newly_joint_q10, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))

write_csv_retry(summary_df,
          file = file.path(OUT_DIR, "sim_summary_v8.csv"),
          row.names = FALSE)

cat(sprintf("\nSimulation complete: %s\n", format(Sys.time())))
cat(sprintf("Results:  %s\n", results_file))
cat(sprintf("Summary:  %s\n", file.path(OUT_DIR, "sim_summary_v8.csv")))
cat(sprintf("Failures: %s\n", failures_file))








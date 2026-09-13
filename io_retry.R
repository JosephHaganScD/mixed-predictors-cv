###############################################################################
# io_retry.R — Robust I/O helpers (patch 2026-07-09)
#
# Retries a write a few times with linear backoff before giving up, so a
# transient file-lock (e.g. sync-related, antivirus scan) doesn't kill an
# unattended run. On repeated failure it stops the script outright rather
# than letting execution continue on incomplete data.
#
# Sourced by both run_simulation.R and run_supplemental.R.
###############################################################################

write_table_retry <- function(x, file, max_attempts = 5, wait_sec = 5, ...) {
  for (attempt in seq_len(max_attempts)) {
    ok <- tryCatch({ write.table(x, file = file, ...); TRUE },
      error = function(e) {
        cat(sprintf("  [write_table_retry] attempt %d/%d failed: %s\n",
                     attempt, max_attempts, conditionMessage(e)))
        FALSE
      })
    if (isTRUE(ok)) return(invisible(TRUE))
    if (attempt < max_attempts) Sys.sleep(wait_sec * attempt)  # 5, 10, 15, 20 sec
  }
  stop(sprintf("write_table_retry: failed to write %s after %d attempts",
               file, max_attempts))
}

saveRDS_retry <- function(object, file, max_attempts = 5, wait_sec = 5) {
  for (attempt in seq_len(max_attempts)) {
    ok <- tryCatch({ saveRDS(object, file = file); TRUE },
      error = function(e) {
        cat(sprintf("  [saveRDS_retry] attempt %d/%d failed: %s\n",
                     attempt, max_attempts, conditionMessage(e)))
        FALSE
      })
    if (isTRUE(ok)) return(invisible(TRUE))
    if (attempt < max_attempts) Sys.sleep(wait_sec * attempt)
  }
  stop(sprintf("saveRDS_retry: failed to write %s after %d attempts",
               file, max_attempts))
}

write_csv_retry <- function(x, file, max_attempts = 5, wait_sec = 5, ...) {
  for (attempt in seq_len(max_attempts)) {
    ok <- tryCatch({ write.csv(x, file = file, ...); TRUE },
      error = function(e) {
        cat(sprintf("  [write_csv_retry] attempt %d/%d failed: %s\n",
                     attempt, max_attempts, conditionMessage(e)))
        FALSE
      })
    if (isTRUE(ok)) return(invisible(TRUE))
    if (attempt < max_attempts) Sys.sleep(wait_sec * attempt)
  }
  stop(sprintf("write_csv_retry: failed to write %s after %d attempts",
               file, max_attempts))
}

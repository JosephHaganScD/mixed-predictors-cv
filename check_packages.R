###############################################################################
# check_packages.R — Verify required packages after R reinstall
# Mixed-Predictor CV Optimism Decomposition Simulation
###############################################################################

cat(sprintf("R version: %s\n\n", R.version.string))

# Every package actually library()'d somewhere across dgm.R, learners.R,
# cv_strategies.R, metrics.R, simulate.R, simulate_paired.R, and the three
# driver scripts, PLUS dplyr, which nothing calls directly but which
# furrr::future_map_dfr() requires internally -- it fails with an unrelated-
# looking error ("future_map_dfr() requires dplyr") if it's missing, not an
# obvious "package not found" message, so it's easy to miss until mid-run.
required_pkgs <- c("glmnet", "xgboost", "pROC", "future", "furrr", "dplyr",
                    "splines")

# splines ships with base R and will always be found; included only for
# completeness against the actual library() calls in the scripts.
check_one <- function(pkg) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  ver <- if (ok) as.character(utils::packageVersion(pkg)) else NA_character_
  data.frame(package = pkg, installed = ok, version = ver, stringsAsFactors = FALSE)
}

status <- do.call(rbind, lapply(required_pkgs, check_one))
print(status, row.names = FALSE)

missing_pkgs <- status$package[!status$installed]

if (length(missing_pkgs) > 0) {
  cat(sprintf("\nMissing: %s\n", paste(missing_pkgs, collapse = ", ")))
  cat("Installing now (you may be prompted to choose a CRAN mirror)...\n\n")
  install.packages(setdiff(missing_pkgs, "splines"))  # splines is base R;
                                                        # install.packages()
                                                        # on it would error

  cat("\nRe-checking after install:\n")
  status2 <- do.call(rbind, lapply(required_pkgs, check_one))
  print(status2, row.names = FALSE)

  if (any(!status2$installed)) {
    still_missing <- status2$package[!status2$installed]
    stop("Still missing after install attempt: ", paste(still_missing, collapse = ", "),
         ". Install manually before running any driver script.")
  }
  cat("\nAll required packages now installed.\n")
} else {
  cat("\nAll required packages already installed.\n")
}

# Quick load test -- confirms each package not only installs but actually
# attaches without error (catches version-mismatch / dependency issues that
# requireNamespace() alone can miss).
cat("\nLoad test:\n")
for (pkg in required_pkgs) {
  ok <- tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    TRUE
  }, error = function(e) {
    cat(sprintf("  %s: FAILED TO LOAD -- %s\n", pkg, conditionMessage(e)))
    FALSE
  })
  if (ok) cat(sprintf("  %s: OK\n", pkg))
}

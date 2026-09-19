# s11_between_iteration_intervals.R
# Supplementary Table S11, mixed predictors manuscript.
# Between-iteration variability of subject-level 10-fold cross-validated AUROC in the
# simulation: Arm A (fixed-only predictors), penalized logistic regression (learner label
# "ridge_linear"), n = 100 subjects; 9 conditions (p_F x R2_total) x 500 iterations.
#
# Input : sim_results_v8.csv, iteration-level simulation output (567,000 rows)
# Output: S11_between_iteration_intervals.csv and Table_S11.docx, written to out_dir
# Packages: data.table, flextable, officer
# No random number generation; results are deterministic given the input file.

library(data.table)
library(flextable)
library(officer)

# ---- Paths: edit these two lines; run with the working directory set to the project folder ----
results_file <- "sim_results_v8.csv"
out_dir      <- "output"
stopifnot(file.exists(results_file))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Empirical comparators (Table 4 of the manuscript, primary outcome, LR fixed-only) ----
emp_auroc <- 0.703               # subject-level 10-fold AUROC
delong_ci <- c(0.603, 0.810)     # DeLong 95% CI of the leave-one-cluster-out AUROC (0.706)

# ---- Read needed columns; subset to Arm A, penalized logistic regression, n = 100 ----
dt  <- fread(results_file, select = c("arm", "condition_id", "iteration", "learner",
                                      "n", "p_F", "R2_total", "auroc_subject"))
sub <- dt[arm == "A_fixed" & learner == "ridge_linear" & n == 100]

# Integrity checks: 9 conditions x 500 iterations
stopifnot(uniqueN(sub$condition_id) == 9, nrow(sub) == 4500,
          all(sub[, .N, by = condition_id]$N == 500))

# ---- Between-iteration summary of subject-level cross-validated AUROC per condition ----
tab <- sub[, .(median = median(auroc_subject),
               lo     = quantile(auroc_subject, 0.025),   # default quantile type 7
               hi     = quantile(auroc_subject, 0.975),
               sd     = sd(auroc_subject)),
           by = .(p_F, R2_total)][order(p_F, R2_total)]
tab[, width := hi - lo]
print(tab[, lapply(.SD, round, 3)])
fwrite(tab, file.path(out_dir, "S11_between_iteration_intervals.csv"))

# ---- Display table (3 decimals); plain column names here, labels applied afterwards ----
f3   <- function(x) formatC(x, format = "f", digits = 3)
disp <- data.frame(p_F = tab$p_F,
                   R2  = formatC(tab$R2_total, format = "f", digits = 2),
                   med = f3(tab$median),
                   lo  = f3(tab$lo),
                   hi  = f3(tab$hi),
                   wid = f3(tab$width),
                   sd  = f3(tab$sd))

fn <- paste0("Arm A (fixed-only predictors), penalized logistic regression, n = 100 subjects, ",
             "outcome prevalence targeted at 48.0%; 500 iterations per condition. Percentiles are of ",
             "subject-level 10-fold cross-validated AUROC across iterations, and width is the difference ",
             "between the 97.5th and 2.5th percentiles. Subject-level cross-validated AUROC is used as a ",
             "proxy for the leave-one-cluster-out estimate, which was not computed in the simulation. ",
             "For comparison, the empirical subject-level AUROC for the fixed-only configuration under ",
             "penalized logistic regression was ", sprintf("%.3f", emp_auroc),
             ", and the width of the DeLong 95% confidence interval for the corresponding ",
             "leave-one-cluster-out estimate (", sprintf("%.3f", delong_ci[1]), " to ",
             sprintf("%.3f", delong_ci[2]), ") was ", sprintf("%.3f", diff(delong_ci)), ".")

ft <- flextable(disp) |>
  set_header_labels(values = c(p_F = "Fixed predictors", R2 = "R\u00B2",
                               med = "Median AUROC", lo = "2.5th percentile",
                               hi = "97.5th percentile", wid = "Width", sd = "SD")) |>
  theme_booktabs() |>
  align(align = "center", part = "all") |>
  add_footer_lines(values = fn) |>
  align(align = "left", part = "footer") |>
  fontsize(size = 9, part = "footer") |>
  autofit() |>
  fit_to_width(max_width = 9)                    # landscape letter, 1-inch margins

doc <- read_docx() |>
  body_add_par("Table S11. Between-iteration variability of subject-level cross-validated AUROC in the simulation, by condition.",
               style = "Normal") |>
  body_add_flextable(ft) |>
  body_end_section_landscape(w = 8.5, h = 11)    # after the content, US Letter
print(doc, target = file.path(out_dir, "Table_S11.docx"))

print(sessionInfo())

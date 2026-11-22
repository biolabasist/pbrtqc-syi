data_dir <- Sys.getenv("PBRTQC_DATA", "../data")
out_dir <- Sys.getenv("PBRTQC_OUT", "../output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source("R/pbrtqc.R")
source("R/analysis.R")

analytes <- utils::read.csv("config/analytes.csv", encoding = "UTF-8", stringsAsFactors = FALSE)

# Allowable total error (%). CLIA percentage criterion; analytes with an absolute
# criterion only are converted at the mean of the truncated stream. Alternatives
# are used in sensitivity analyses.
tea_percent <- function(a, mean_kept, rule = c("primary", "greater", "eflm")) {
  rule <- match.arg(rule)
  abs_pct <- if (is.na(a$clia_abs)) NA else 100 * a$clia_abs / mean_kept
  switch(rule,
         primary = if (!is.na(a$clia_pct)) a$clia_pct else abs_pct,
         greater = max(a$clia_pct, abs_pct, na.rm = TRUE),
         eflm = 1.65 * 0.5 * a$cvi + 0.25 * sqrt(a$cvi^2 + a$cvg^2))
}

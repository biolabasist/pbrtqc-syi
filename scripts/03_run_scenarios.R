# Error-detection simulations for the primary analysis and sensitivity analyses.
# Usage (from the repository root): Rscript scripts/03_run_scenarios.R <scenario>
#
# primary      combined stream of both analyzers, full year
# analyzer1    error confined to analyzer 1, monitored on the combined stream
# stream_a1    analyzer 1 monitored as a separate stream
# stream_a2    analyzer 2 monitored as a separate stream
# derivation   January-June (parameters and selection)
# validation   July-December, monitored with the January-June parameters
# cva_measured truncation with the measured analytical CV
# median       truncation limits centred on the median
# logsym       median-centred limits, symmetric on the logarithmic scale
# tea_greater  CLIA criterion as the greater of the percentage and absolute limit
# tea_eflm     desirable total allowable error from biological variation
# permuted     randomly permuted stream (order effect)

source("R/setup.R")

scenario <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(scenario)) scenario <- "primary"

stream <- readRDS(file.path(data_dir, "stream_2023.rds"))
cva <- readRDS(file.path(out_dir, "iqc_cv.rds"))
split_date <- as.POSIXct("2023-07-01", tz = "UTC")

primary_model <- function(a, x) stream_model(x, rcv_g(a$cvi, a$cvg))

run_one <- function(a) {
  d <- stream[stream$analyte == a$analyte, ]
  base <- primary_model(a, d$value)
  tea <- tea_percent(a, base$mean_kept, "primary")
  x <- d$value; time <- d$time; model <- base; affected <- NULL; starts <- NULL

  if (scenario == "analyzer1") {
    affected <- d$analyzer %in% 1L
    starts <- readRDS(file.path(out_dir, "scenario_primary.rds"))[[a$analyte]]$starts
  } else if (scenario %in% c("stream_a1", "stream_a2")) {
    k <- as.integer(substr(scenario, 9, 9))
    sel <- d$analyzer %in% k
    x <- x[sel]; time <- time[sel]; model <- primary_model(a, x)
  } else if (scenario %in% c("derivation", "validation")) {
    h1 <- d$time < split_date
    model <- primary_model(a, d$value[h1])
    sel <- if (scenario == "derivation") h1 else !h1
    x <- x[sel]; time <- time[sel]
  } else if (scenario == "cva_measured") {
    cv <- cva$cv_pooled[cva$analyte == a$analyte]
    model <- stream_model(x, rcv_g(a$cvi, a$cvg, cv))
  } else if (scenario == "median") {
    model <- stream_model(x, rcv_g(a$cvi, a$cvg), centre = "median")
  } else if (scenario == "logsym") {
    model <- stream_model(x, rcv_g(a$cvi, a$cvg), centre = "median", form = "multiplicative")
  } else if (scenario == "tea_greater") {
    tea <- tea_percent(a, base$mean_kept, "greater")
  } else if (scenario == "tea_eflm") {
    tea <- tea_percent(a, base$mean_kept, "eflm")
  } else if (scenario == "permuted") {
    set.seed(1)
    x <- x[sample.int(length(x))]
  }

  mult <- if (scenario %in% c("tea_greater", "tea_eflm")) bias_multiples[abs(bias_multiples) <= 1] else bias_multiples
  t0 <- Sys.time()
  res <- analyse_stream(x, time, model, tea, multiples = mult, affected = affected, starts = starts)
  message(sprintf("%s %s: %.0f s", scenario, a$analyte, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  list(ic = res$ic, sim = res$sim, starts = res$starts, tea = tea, n = length(x),
       model = model[c("lim", "offset", "mu0", "s0", "mean_kept", "excluded")])
}

out <- lapply(split(analytes, seq_len(nrow(analytes))), run_one)
names(out) <- analytes$analyte
saveRDS(out, file.path(out_dir, sprintf("scenario_%s.rds", scenario)))

# Sensitivity analysis: an error counts as detected only at the first signal in
# the direction of the error (signals in the opposite direction are ignored).
# Uses the onsets, models and in-control results of the primary analysis.
# Also records, for every configuration, the probability of at least one alarm
# within B results after the onsets when no error is introduced.

source("R/setup.R")
source("scripts/sim_direction.R")

stream <- readRDS(file.path(data_dir, "stream_2023.rds"))
prim <- readRDS(file.path(out_dir, "scenario_primary.rds"))
grid <- chart_grid()

out <- list()
for (k in analytes$analyte) {
  t0 <- Sys.time()
  p <- prim[[k]]
  x <- stream$value[stream$analyte == k]
  mults <- c(0, bias_multiples)
  sims <- lapply(seq_len(nrow(grid)), function(i) {
    cfg <- grid[i, ]
    run <- run_chart(x, p$model, cfg)
    r <- sim_direction(x, p$model, cfg, run, p$starts, 1 + mults * p$tea / 100, 1000)
    r$multiple <- rep(mults, each = length(p$starts))
    r$config <- cfg$config
    r
  })
  s <- do.call(rbind, sims)
  s$any <- pmin(s$up, s$dn, na.rm = TRUE)
  null <- s[s$multiple == 0, c("config", "start", "any")]
  s <- s[s$multiple != 0, ]
  s$first <- ifelse(s$multiple > 0, s$up, s$dn)
  s$bias_pct <- s$multiple * p$tea
  out[[k]] <- list(ic = p$ic, sim = s[, c("start", "mult", "first", "config", "multiple", "bias_pct")],
                   null = null, starts = p$starts, tea = p$tea, n = p$n, model = p$model)
  message(sprintf("direction %s: %.0f s", k, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
saveRDS(out, file.path(out_dir, "scenario_direction.rds"))

# Grid extension: configurations beyond the edges of the primary grid
# (EWMA lambda 0.01 and 0.30; MA w = 5; MA k = 2.5), evaluated with the onsets
# and models of the primary analysis so that they can be compared directly with
# the 53 primary configurations.

source("R/setup.R")

stream <- readRDS(file.path(data_dir, "stream_2023.rds"))
prim <- readRDS(file.path(out_dir, "scenario_primary.rds"))

label <- function(g) ifelse(g$method == "EWMA", sprintf("EWMA λ=%.2f L=%.1f", g$lambda, g$L),
                            sprintf("MA w=%d k=%.1f", as.integer(g$w), g$L))
e <- expand.grid(lambda = c(0.01, 0.30), L = c(2.5, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0))
e <- data.frame(method = "EWMA", lambda = e$lambda, L = e$L, w = NA_real_)
m1 <- data.frame(method = "MA", lambda = NA_real_, L = c(2.5, 3.0, 3.5, 4.0, 4.5, 5.0), w = 5)
m2 <- data.frame(method = "MA", lambda = NA_real_, L = 2.5, w = c(10, 25, 50, 75, 100))
extra <- rbind(e, m1, m2)
extra$config <- label(extra)

out <- list()
for (k in analytes$analyte) {
  t0 <- Sys.time()
  p <- prim[[k]]
  d <- stream[stream$analyte == k, ]
  runs <- lapply(seq_len(nrow(extra)), function(i) run_chart(d$value, p$model, extra[i, ]))
  ic <- do.call(rbind, lapply(seq_len(nrow(extra)), function(i) {
    st <- runs[[i]]$status
    ep <- alarm_episodes(st, d$time)
    data.frame(config = extra$config[i], spec = 1 - mean(st), far = mean(st),
               episodes = ep$n, episodes_per_day = ep$per_day)
  }))
  mult <- 1 + bias_multiples * p$tea / 100
  sim <- do.call(rbind, lapply(seq_len(nrow(extra)), function(i) {
    r <- simulate_onsets(d$value, p$model, extra[i, ], runs[[i]], p$starts, mult, 1000)
    r$config <- extra$config[i]
    r$multiple <- rep(bias_multiples, each = length(p$starts))
    r
  }))
  sim$bias_pct <- sim$multiple * p$tea
  # Onsets of the primary analysis were chosen where the primary configurations
  # were in control; record whether each extra configuration was in control too.
  inctl <- vapply(seq_len(nrow(extra)), function(i) mean(!runs[[i]]$status[p$starts - 1]), numeric(1))
  ic$in_control_at_onsets <- inctl
  out[[k]] <- list(ic = ic, sim = sim)
  message(sprintf("grid %s: %.0f s", k, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
saveRDS(list(grid = extra, res = out), file.path(out_dir, "scenario_grid_extension.rds"))

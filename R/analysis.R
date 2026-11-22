# Stream-level analysis and summary metrics.

bias_multiples <- c(-3, -2, -1.5, -1, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1, 1.5, 2, 3)

analyse_stream <- function(x, time, model, tea_pct, grid = chart_grid(),
                           multiples = bias_multiples, B = 1000, K = 100,
                           affected = NULL, starts = NULL, seed = 2023) {
  runs <- lapply(seq_len(nrow(grid)), function(i) run_chart(x, model, grid[i, ]))
  status <- vapply(runs, `[[`, logical(length(x)), "status")

  ic <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    ep <- alarm_episodes(status[, i], time)
    data.frame(config = grid$config[i], spec = 1 - mean(status[, i]),
               far = mean(status[, i]), episodes = ep$n, episodes_per_day = ep$per_day)
  }))

  if (is.null(starts)) starts <- select_onsets(status, K = K, seed = seed)
  mult <- 1 + multiples * tea_pct / 100
  sim <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    r <- simulate_onsets(x, model, grid[i, ], runs[[i]], starts, mult, B, affected)
    r$config <- grid$config[i]
    r$multiple <- rep(multiples, each = length(starts))
    r
  }))
  sim$bias_pct <- sim$multiple * tea_pct
  list(ic = ic, sim = sim, starts = starts, status = status, runs = runs)
}

# Youden index per configuration and bias level for block length B.
youden_table <- function(res, B) {
  s <- res$sim
  s$sens <- block_sensitivity(s$first, B)
  a <- stats::aggregate(sens ~ config + multiple + bias_pct, s, mean)
  a <- merge(a, res$ic[, c("config", "spec")], by = "config")
  a$J <- a$sens + a$spec - 1
  a[order(a$config, a$multiple), ]
}

# Summed Youden index over |bias| <= 1 x TEa with a paired bootstrap over onset
# positions. The reference is the best configuration among those whose alarm
# burden in the in-control stream does not exceed `max_episodes` per day;
# eligible configurations whose SYI difference from the reference has a 95%
# interval including zero are flagged as equivalent.
syi_table <- function(res, B, grid = chart_grid(), nboot = 2000, seed = 2023,
                      within = 1, max_episodes = Inf) {
  s <- res$sim[abs(res$sim$multiple) <= within, ]
  s$sens <- block_sensitivity(s$first, B)
  nlev <- length(unique(s$multiple))
  per <- stats::aggregate(sens ~ config + start, s, sum)
  per <- merge(per, res$ic[, c("config", "spec")], by = "config")
  per$syi <- per$sens + nlev * (per$spec - 1)
  wide <- stats::reshape(per[, c("config", "start", "syi")], idvar = "start",
                         timevar = "config", direction = "wide")
  m <- as.matrix(wide[, -1]); colnames(m) <- sub("^syi\\.", "", colnames(m))
  est <- colMeans(m)
  burden <- res$ic$episodes_per_day[match(names(est), res$ic$config)]
  eligible <- burden <= max_episodes
  if (!any(eligible)) eligible <- burden == min(burden)
  best <- names(est)[eligible][which.max(est[eligible])]
  set.seed(seed)
  bt <- replicate(nboot, colMeans(m[sample.int(nrow(m), replace = TRUE), , drop = FALSE]))
  diff <- bt - matrix(bt[best, ], nrow = nrow(bt), ncol = nboot, byrow = TRUE)
  out <- data.frame(config = names(est), syi = est,
                    lo = apply(bt, 1, stats::quantile, 0.025),
                    hi = apply(bt, 1, stats::quantile, 0.975),
                    diff_lo = apply(diff, 1, stats::quantile, 0.025),
                    diff_hi = apply(diff, 1, stats::quantile, 0.975),
                    eligible = eligible)
  out$equivalent <- out$eligible & out$diff_hi >= 0
  out$best <- out$config == best
  out <- merge(out, res$ic, by = "config")
  out <- merge(grid[, c("config", "method", "lambda", "L", "w")], out, by = "config")
  out[order(-out$syi), ]
}

# Detection delay for large errors: median number of results affected until the
# first signal (censored at B) and probability of detection within B.
delay_table <- function(res, B) {
  s <- res$sim
  s$d <- ifelse(is.na(s$first) | s$first > B, B, s$first)
  s$det <- !is.na(s$first) & s$first <= B
  a <- stats::aggregate(cbind(d, det, trunc_frac) ~ config + multiple + bias_pct, s,
                        function(v) c(med = stats::median(v), mean = mean(v),
                                      p95 = stats::quantile(v, 0.95, names = FALSE)))
  data.frame(config = a$config, multiple = a$multiple, bias_pct = a$bias_pct,
             mnped = a$d[, "med"], anped = a$d[, "mean"], p95 = a$d[, "p95"],
             p_detect = a$det[, "mean"], trunc_frac = a$trunc_frac[, "mean"])
}

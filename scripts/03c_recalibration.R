# Temporal validation with monthly re-estimation of the in-control mean and SD.
# Truncation limits and the configuration are derived from January-June, as in
# the "validation" scenario; in July-December, mu0 and sigma0 of each calendar
# month are the mean and SD of the transformed, truncated results of the
# preceding month. All 53 configurations are evaluated as in the primary
# analysis (in-control run, onsets, error simulation).

source("R/setup.R")

stream <- readRDS(file.path(data_dir, "stream_2023.rds"))
split_date <- as.POSIXct("2023-07-01", tz = "UTC")
grid <- chart_grid()

# Chart with a mean and SD that may change at every result.
run_chart_tv <- function(x, model, cfg, mu, sd) {
  ins <- inside_limits(x, model$lim)
  y <- log(x[ins] + model$offset)
  mu_m <- mu[ins]; sd_m <- sd[ins]
  n <- seq_along(y)
  if (cfg$method == "EWMA") {
    lam <- cfg$lambda
    stat <- as.numeric(stats::filter(lam * y, 1 - lam, method = "recursive", init = mu[1]))
    hw <- cfg$L * sd_m * sqrt(lam / (2 - lam) * (1 - (1 - lam)^(2 * n)))
  } else {
    w <- cfg$w
    cs <- cumsum(y)
    k <- pmin(n, w)
    stat <- (cs - c(rep(0, w), cs)[seq_along(cs)]) / k
    hw <- cfg$L * sd_m / sqrt(k)
  }
  alarm_mon <- abs(stat - mu_m) > hw
  k_in <- cumsum(ins)
  status <- logical(length(x))
  status[k_in > 0] <- alarm_mon[k_in[k_in > 0]]
  list(status = status, stat = stat, y = y, k_in = k_in, n_mon = length(y))
}

simulate_tv <- function(x, model, cfg, run, starts, mult, B, mu, sd) {
  K <- length(starts); nb <- length(mult)
  idx <- outer(0:(B - 1), starts, "+")
  raw <- matrix(x[idx], nrow = B)[, rep(seq_len(K), times = nb), drop = FALSE]
  mum <- matrix(mu[idx], nrow = B)[, rep(seq_len(K), times = nb), drop = FALSE]
  sdm <- matrix(sd[idx], nrow = B)[, rep(seq_len(K), times = nb), drop = FALSE]
  xb <- raw * matrix(rep(mult, each = K * B), nrow = B)
  ins <- xb >= model$lim[1] & xb <= model$lim[2]
  y <- log(pmax(xb, 0) + model$offset)
  ncol <- K * nb
  first <- rep(NA_integer_, ncol)
  st0 <- lapply(starts, function(s) chart_state(run, s, cfg))
  n <- rep(vapply(st0, `[[`, numeric(1), "n"), times = nb)
  if (cfg$method == "EWMA") {
    lam <- cfg$lambda
    z <- rep(vapply(st0, `[[`, numeric(1), "z"), times = nb)
    for (t in seq_len(B)) {
      m <- ins[t, ]
      z[m] <- lam * y[t, m] + (1 - lam) * z[m]
      n[m] <- n[m] + 1
      hw <- cfg$L * sdm[t, ] * sqrt(lam / (2 - lam) * (1 - (1 - lam)^(2 * n)))
      hit <- is.na(first) & abs(z - mum[t, ]) > hw
      first[hit] <- t
    }
  } else {
    w <- cfg$w
    buf <- matrix(unlist(lapply(st0, `[[`, "buf")), nrow = w)[, rep(seq_len(K), times = nb), drop = FALSE]
    pos <- rep(1L, ncol)
    s <- colSums(buf)
    for (t in seq_len(B)) {
      m <- which(ins[t, ])
      if (length(m)) {
        old <- buf[cbind(pos[m], m)]
        buf[cbind(pos[m], m)] <- y[t, m]
        s[m] <- s[m] - old + y[t, m]
        pos[m] <- pos[m] %% w + 1L
      }
      hit <- is.na(first) & abs(s / w - mum[t, ]) > cfg$L * sdm[t, ] / sqrt(w)
      first[hit] <- t
    }
  }
  data.frame(start = rep(starts, times = nb), mult = rep(mult, each = K), first = first)
}

out <- list()
for (i in seq_len(nrow(analytes))) {
  a <- analytes[i, ]
  t0 <- Sys.time()
  d <- stream[stream$analyte == a$analyte, ]
  base <- stream_model(d$value, rcv_g(a$cvi, a$cvg))
  tea <- tea_percent(a, base$mean_kept, "primary")
  h1 <- d$time < split_date
  model <- stream_model(d$value[h1], rcv_g(a$cvi, a$cvg))

  # Monthly mean and SD of the transformed, truncated results.
  ins <- inside_limits(d$value, model$lim)
  mon <- format(d$time, "%Y-%m")
  y <- log(d$value + model$offset)
  ms <- stats::aggregate(y[ins], list(month = mon[ins]), function(v) c(m = mean(v), s = stats::sd(v)))
  ms <- data.frame(month = ms$month, mu = ms$x[, "m"], sd = ms$x[, "s"])
  prev <- function(mo) format(seq(as.Date(paste0(mo, "-01")), by = "-1 month", length.out = 2)[2], "%Y-%m")

  x <- d$value[!h1]; time <- d$time[!h1]; mo2 <- mon[!h1]
  um <- unique(mo2)
  pm <- setNames(vapply(um, prev, character(1)), um)[mo2]
  mu <- ms$mu[match(pm, ms$month)]
  sd <- ms$sd[match(pm, ms$month)]

  runs <- lapply(seq_len(nrow(grid)), function(j) run_chart_tv(x, model, grid[j, ], mu, sd))
  status <- vapply(runs, `[[`, logical(length(x)), "status")
  ic <- do.call(rbind, lapply(seq_len(nrow(grid)), function(j) {
    ep <- alarm_episodes(status[, j], time)
    data.frame(config = grid$config[j], spec = 1 - mean(status[, j]),
               far = mean(status[, j]), episodes = ep$n, episodes_per_day = ep$per_day)
  }))
  starts <- select_onsets(status, K = 100)
  mult <- 1 + bias_multiples * tea / 100
  sim <- do.call(rbind, lapply(seq_len(nrow(grid)), function(j) {
    r <- simulate_tv(x, model, grid[j, ], runs[[j]], starts, mult, 1000, mu, sd)
    r$config <- grid$config[j]
    r$multiple <- rep(bias_multiples, each = length(starts))
    r
  }))
  sim$bias_pct <- sim$multiple * tea
  sim$trunc_frac <- NA_real_
  out[[a$analyte]] <- list(ic = ic, sim = sim, starts = starts, tea = tea, n = length(x),
                           model = model[c("lim", "offset", "mu0", "s0", "mean_kept", "excluded")],
                           monthly = ms)
  message(sprintf("recalibration %s: %.0f s", a$analyte, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
saveRDS(out, file.path(out_dir, "scenario_recalibration.rds"))

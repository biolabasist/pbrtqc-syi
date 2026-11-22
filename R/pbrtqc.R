# Core functions: truncation, transformation, EWMA/MA charts and error simulation.

# Truncation -----------------------------------------------------------------

# Between-subject reference change value (Coskun et al., 2021), in percent.
# The analytical term defaults to the desirable specification CVa = 0.5 * CVi.
rcv_g <- function(cvi, cvg, cva = 0.5 * cvi) {
  sqrt(2) * 1.96 * sqrt(cvg^2 + cvi^2 + cva^2)
}

truncation_limits <- function(x, rcv, centre = c("mean", "median"),
                              form = c("additive", "multiplicative")) {
  centre <- match.arg(centre)
  form <- match.arg(form)
  c0 <- if (centre == "mean") mean(x) else stats::median(x)
  f <- rcv / 100
  lim <- if (form == "additive") c0 * c(1 - f, 1 + f) else c(c0 / (1 + f), c0 * (1 + f))
  lim[1] <- max(lim[1], 0)
  lim
}

inside_limits <- function(x, lim) x >= lim[1] & x <= lim[2]

# ln(x + 1) is used when any monitored value is below 1 in the working units.
log_offset <- function(x) if (any(x < 1)) 1 else 0

# In-control model of a stream: truncation limits, offset, mean and SD of
# the transformed, truncated results.
stream_model <- function(x, rcv, centre = "mean", form = "additive") {
  lim <- truncation_limits(x, rcv, centre, form)
  kept <- x[inside_limits(x, lim)]
  off <- log_offset(kept)
  y <- log(kept + off)
  list(lim = lim, offset = off, mu0 = mean(y), s0 = stats::sd(y),
       mean_kept = mean(kept), excluded = 1 - length(kept) / length(x))
}

# Configurations ---------------------------------------------------------------

chart_grid <- function(lambda = c(0.02, 0.05, 0.10, 0.20),
                       L = c(2.5, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0),
                       w = c(10, 25, 50, 75, 100), k = c(3.0, 3.5, 4.0, 4.5, 5.0)) {
  e <- expand.grid(lambda = lambda, L = L)
  e <- data.frame(method = "EWMA", lambda = e$lambda, L = e$L, w = NA_real_)
  m <- expand.grid(w = w, L = k)
  m <- data.frame(method = "MA", lambda = NA_real_, L = m$L, w = m$w)
  g <- rbind(e, m)
  g$config <- ifelse(g$method == "EWMA",
                     sprintf("EWMA λ=%.2f L=%.1f", g$lambda, g$L),
                     sprintf("MA w=%d k=%.1f", as.integer(g$w), g$L))
  g
}

# Continuous monitoring of an in-control stream ------------------------------

# Runs one configuration over an ordered stream of raw results. Results outside
# the truncation limits do not update the statistic; the chart status at such a
# result is the status after the previous monitored result.
run_chart <- function(x, model, cfg) {
  ins <- inside_limits(x, model$lim)
  y <- log(x[ins] + model$offset)
  n <- seq_along(y)
  if (cfg$method == "EWMA") {
    lam <- cfg$lambda
    stat <- as.numeric(stats::filter(lam * y, 1 - lam, method = "recursive",
                                     init = model$mu0))
    hw <- cfg$L * model$s0 * sqrt(lam / (2 - lam) * (1 - (1 - lam)^(2 * n)))
  } else {
    w <- cfg$w
    cs <- cumsum(y)
    k <- pmin(n, w)
    stat <- (cs - c(rep(0, w), cs)[seq_along(cs)]) / k
    hw <- cfg$L * model$s0 / sqrt(k)
  }
  alarm_mon <- abs(stat - model$mu0) > hw
  k_in <- cumsum(ins)
  status <- logical(length(x))
  status[k_in > 0] <- alarm_mon[k_in[k_in > 0]]
  list(status = status, stat = stat, y = y, k_in = k_in, n_mon = length(y))
}

# Chart state just before result s (used to start an error on a running chart).
chart_state <- function(run, s, cfg) {
  k <- run$k_in[s - 1]
  if (cfg$method == "EWMA") {
    list(z = run$stat[k], n = k)
  } else {
    w <- cfg$w
    list(buf = run$y[(k - w + 1):k], n = k)
  }
}

# Error simulation -------------------------------------------------------------

# Introduces proportional systematic errors at onset positions `starts` and
# follows the chart for B results. `mult` is a vector of multipliers (1 + bias).
# `affected` optionally restricts the error to a subset of results (e.g. one
# analyzer). Returns the index (1..B) of the first out-of-control status or NA,
# and the fraction of affected results removed by truncation.
simulate_onsets <- function(x, model, cfg, run, starts, mult, B,
                            affected = NULL) {
  K <- length(starts)
  nb <- length(mult)
  idx <- outer(0:(B - 1), starts, "+")
  raw <- matrix(x[idx], nrow = B)
  aff <- if (is.null(affected)) matrix(TRUE, B, K) else matrix(affected[idx], nrow = B)
  raw <- raw[, rep(seq_len(K), times = nb), drop = FALSE]
  aff <- aff[, rep(seq_len(K), times = nb), drop = FALSE]
  fac <- matrix(rep(mult, each = K * B), nrow = B)
  xb <- ifelse(aff, raw * fac, raw)
  ins <- xb >= model$lim[1] & xb <= model$lim[2]
  y <- log(pmax(xb, 0) + model$offset)
  trunc_frac <- colSums(!ins & aff) / pmax(colSums(aff), 1)

  ncol <- K * nb
  first <- rep(NA_integer_, ncol)
  st0 <- lapply(starts, function(s) chart_state(run, s, cfg))
  n <- rep(vapply(st0, `[[`, numeric(1), "n"), times = nb)
  mu0 <- model$mu0
  s0 <- model$s0

  if (cfg$method == "EWMA") {
    lam <- cfg$lambda
    z <- rep(vapply(st0, `[[`, numeric(1), "z"), times = nb)
    for (t in seq_len(B)) {
      m <- ins[t, ]
      z[m] <- lam * y[t, m] + (1 - lam) * z[m]
      n[m] <- n[m] + 1
      hw <- cfg$L * s0 * sqrt(lam / (2 - lam) * (1 - (1 - lam)^(2 * n)))
      hit <- is.na(first) & abs(z - mu0) > hw
      first[hit] <- t
    }
  } else {
    w <- cfg$w
    buf <- matrix(unlist(lapply(st0, `[[`, "buf")), nrow = w)
    buf <- buf[, rep(seq_len(K), times = nb), drop = FALSE]
    pos <- rep(1L, ncol)
    s <- colSums(buf)
    hw <- cfg$L * s0 / sqrt(w)
    for (t in seq_len(B)) {
      m <- which(ins[t, ])
      if (length(m)) {
        old <- buf[cbind(pos[m], m)]
        buf[cbind(pos[m], m)] <- y[t, m]
        s[m] <- s[m] - old + y[t, m]
        pos[m] <- pos[m] %% w + 1L
      }
      hit <- is.na(first) & abs(s / w - mu0) > hw
      first[hit] <- t
    }
  }
  data.frame(start = rep(starts, times = nb), mult = rep(mult, each = K),
             first = first, trunc_frac = trunc_frac)
}

# Onset positions ---------------------------------------------------------------

# Stratified selection of K onset positions across the stream at which every
# configuration is in control at the preceding result.
select_onsets <- function(status_matrix, K = 100, margin = 1000, seed = 2023) {
  N <- nrow(status_matrix)
  ok <- rowSums(status_matrix) == 0
  cand <- which(c(FALSE, ok[-N]))
  cand <- cand[cand > margin & cand <= N - margin]
  set.seed(seed)
  bins <- cut(seq_along(cand), K, labels = FALSE)
  sort(vapply(split(cand, bins), function(v) v[sample.int(length(v), 1)], numeric(1)))
}

# Metrics ---------------------------------------------------------------------

# Per-result sensitivity within a block of B results when the laboratory acts at
# the first signal: (B - d + 1) / B, zero when no signal occurs within B.
block_sensitivity <- function(first, B) ifelse(!is.na(first) & first <= B, (B - first + 1) / B, 0)

alarm_episodes <- function(status, time) {
  onset <- which(status & !c(FALSE, status[-length(status)]))
  days <- length(unique(as.Date(time)))
  list(n = length(onset), per_day = length(onset) / days, onset_time = time[onset])
}

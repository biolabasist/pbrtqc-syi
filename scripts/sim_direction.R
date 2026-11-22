# Error simulation recording the first upward and the first downward signal.

sim_direction <- function(x, model, cfg, run, starts, mult, B) {
  K <- length(starts); nb <- length(mult)
  idx <- outer(0:(B - 1), starts, "+")
  raw <- matrix(x[idx], nrow = B)
  raw <- raw[, rep(seq_len(K), times = nb), drop = FALSE]
  fac <- matrix(rep(mult, each = K * B), nrow = B)
  xb <- raw * fac
  ins <- xb >= model$lim[1] & xb <= model$lim[2]
  y <- log(pmax(xb, 0) + model$offset)
  ncol <- K * nb
  up <- dn <- rep(NA_integer_, ncol)
  st0 <- lapply(starts, function(s) chart_state(run, s, cfg))
  n <- rep(vapply(st0, `[[`, numeric(1), "n"), times = nb)
  mu0 <- model$mu0; s0 <- model$s0
  if (cfg$method == "EWMA") {
    lam <- cfg$lambda
    z <- rep(vapply(st0, `[[`, numeric(1), "z"), times = nb)
    for (t in seq_len(B)) {
      m <- ins[t, ]
      z[m] <- lam * y[t, m] + (1 - lam) * z[m]
      n[m] <- n[m] + 1
      hw <- cfg$L * s0 * sqrt(lam / (2 - lam) * (1 - (1 - lam)^(2 * n)))
      h <- is.na(up) & (z - mu0) > hw; up[h] <- t
      h <- is.na(dn) & (mu0 - z) > hw; dn[h] <- t
    }
  } else {
    w <- cfg$w
    buf <- matrix(unlist(lapply(st0, `[[`, "buf")), nrow = w)
    buf <- buf[, rep(seq_len(K), times = nb), drop = FALSE]
    pos <- rep(1L, ncol); s <- colSums(buf); hw <- cfg$L * s0 / sqrt(w)
    for (t in seq_len(B)) {
      m <- which(ins[t, ])
      if (length(m)) {
        old <- buf[cbind(pos[m], m)]
        buf[cbind(pos[m], m)] <- y[t, m]
        s[m] <- s[m] - old + y[t, m]
        pos[m] <- pos[m] %% w + 1L
      }
      h <- is.na(up) & (s / w - mu0) > hw; up[h] <- t
      h <- is.na(dn) & (mu0 - s / w) > hw; dn[h] <- t
    }
  }
  data.frame(start = rep(starts, times = nb), mult = rep(mult, each = K), up = up, dn = dn)
}


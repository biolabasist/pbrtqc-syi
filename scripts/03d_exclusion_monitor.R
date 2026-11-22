# Exclusion-rate monitor: number of results outside the truncation limits among
# the last n = 50 results. The limit h is the smallest count for which the
# monitor alone produces no more than 0.1 alarm episodes per day in the 2023
# stream. The monitor is evaluated alone and combined with the recommended
# PBRTQC configuration (alarm when either signals), using the onsets and error
# sizes of the primary analysis.

source("R/setup.R")

n_win <- 50
max_ep <- 0.1

stream <- readRDS(file.path(data_dir, "stream_2023.rds"))
prim <- readRDS(file.path(out_dir, "scenario_primary.rds"))
rec <- read.csv(file.path(out_dir, "tables", "recommended.csv"))
iqc_cv <- readRDS(file.path(out_dir, "iqc_cv.rds"))
b_iqc <- setNames(iqc_cv$b_iqc, iqc_cv$analyte)
grid <- chart_grid()

moving_count <- function(flag, n) {
  cs <- cumsum(flag)
  cs - c(rep(0, n), cs)[seq_along(cs)]
}

episodes <- function(status, time) alarm_episodes(status, time)$per_day

summ <- list(); det <- list()
for (k in analytes$analyte) {
  d <- stream[stream$analyte == k, ]
  p <- prim[[k]]
  m <- p$model
  B <- b_iqc[[k]]
  cfgname <- rec$config[rec$analyte == k]
  cfg <- grid[grid$config == cfgname, ]
  out <- !inside_limits(d$value, m$lim)
  cnt <- moving_count(out, n_win)

  hs <- seq_len(n_win)
  ep <- vapply(hs, function(h) episodes(cnt > h, d$time), numeric(1))
  h <- hs[which(ep <= max_ep)[1]]
  mon <- cnt > h
  chart <- run_chart(d$value, m, cfg)$status
  summ[[k]] <- data.frame(analyte = k, config = cfgname, p0 = mean(out), h = h,
                          monitor_episodes = episodes(mon, d$time),
                          chart_episodes = episodes(chart, d$time),
                          combined_episodes = episodes(mon | chart, d$time))

  s0 <- p$sim[p$sim$config == cfgname, c("start", "multiple", "first")]
  mults <- c(0, bias_multiples)
  rows <- list()
  for (mlt in mults) {
    f_mon <- vapply(p$starts, function(s) {
      xb <- d$value[s:(s + B - 1)] * (1 + mlt * p$tea / 100)
      v <- c(out[(s - n_win + 1):(s - 1)], !inside_limits(xb, m$lim))
      c2 <- moving_count(v, n_win)[n_win:(n_win + B - 1)]
      i <- which(c2 > h)[1]
      if (is.na(i)) NA_integer_ else i
    }, integer(1))
    if (mlt == 0) {
      rows[[length(rows) + 1]] <- data.frame(multiple = 0, p_chart = NA, p_combined = NA,
                                             mnped_chart = NA, mnped_combined = NA,
                                             p_monitor = mean(!is.na(f_mon)))
      next
    }
    f_ch <- s0$first[s0$multiple == mlt][match(p$starts, s0$start[s0$multiple == mlt])]
    f_ch[!is.na(f_ch) & f_ch > B] <- NA
    f_cb <- pmin(f_ch, f_mon, na.rm = TRUE)
    dd <- function(f) stats::median(ifelse(is.na(f), B, f))
    rows[[length(rows) + 1]] <- data.frame(multiple = mlt,
                                           p_chart = mean(!is.na(f_ch)), p_combined = mean(!is.na(f_cb)),
                                           mnped_chart = dd(f_ch), mnped_combined = dd(f_cb),
                                           p_monitor = mean(!is.na(f_mon)))
  }
  r <- do.call(rbind, rows); r$analyte <- k
  det[[k]] <- r
  message(k)
}
summ <- do.call(rbind, summ)
det <- do.call(rbind, det)
write.csv(summ, file.path(out_dir, "tables", "exclusion_monitor_summary.csv"), row.names = FALSE)
write.csv(det, file.path(out_dir, "tables", "exclusion_monitor_detection.csv"), row.names = FALSE)
print(summ, digits = 3, row.names = FALSE)

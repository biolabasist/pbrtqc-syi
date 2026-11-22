# Recommended configurations applied to the full 2023 stream: alarm burden and
# agreement with the internal quality control (IQC) record.

source("R/setup.R")

tab_dir <- file.path(out_dir, "tables")
stream <- readRDS(file.path(data_dir, "stream_2023.rds"))
iqc <- readRDS(file.path(data_dir, "iqc_2023.rds"))
rec <- utils::read.csv(file.path(tab_dir, "recommended.csv"), encoding = "UTF-8")
prim <- readRDS(file.path(out_dir, "scenario_primary.rds"))
grid <- chart_grid()

# IQC deviation (%) from the median of each analyzer module, level and control lot
iqc$dev <- 100 * (iqc$value / stats::ave(iqc$value, iqc$iqc_name, iqc$analyzer, iqc$module,
                                         iqc$level, iqc$lot, FUN = stats::median) - 1)
iqc$day <- as.Date(iqc$time)
reject <- c("1:3S", "2:2S", "R:4S", "4:1S", "3:1S")

out <- lapply(seq_len(nrow(analytes)), function(i) {
  a <- analytes[i, ]
  d <- stream[stream$analyte == a$analyte, ]
  m <- prim[[a$analyte]]$model
  r <- rec[rec$analyte == a$analyte, ]
  cfg <- grid[grid$config == r$config, ]
  run <- run_chart(d$value, m, cfg)
  ep <- alarm_episodes(run$status, d$time)
  onset <- which(run$status & !c(FALSE, run$status[-length(run$status)]))
  k <- run$k_in[onset]
  dir_pbrtqc <- sign(run$stat[k] - m$mu0)

  q <- iqc[iqc$iqc_name == a$iqc_name, ]
  qd <- tapply(q$dev, q$day, mean)
  qabs <- tapply(abs(q$dev), q$day, mean)
  flag_days <- unique(q$day[q$westgard %in% reject])

  alarm_days <- unique(as.Date(d$time[onset]))
  iqc_days <- as.Date(names(qd))
  on_alarm <- iqc_days %in% alarm_days
  w <- if (sum(on_alarm) > 2 && sum(!on_alarm) > 2)
    stats::wilcox.test(qabs[on_alarm], qabs[!on_alarm], exact = FALSE)$p.value else NA

  # direction of the IQC deviation on the day of each alarm onset
  od <- as.character(as.Date(d$time[onset]))
  iqc_dir <- sign(qd[od])
  concord <- mean(iqc_dir == dir_pbrtqc, na.rm = TRUE)

  # daily patient mean against daily IQC deviation, per analyzer
  ins <- inside_limits(d$value, m$lim)
  y <- log(d$value + m$offset)
  rho <- sapply(1:2, function(an) {
    sel <- ins & d$analyzer %in% an
    pm <- tapply(y[sel], as.Date(d$time[sel]), mean)
    pm <- 100 * (exp(pm - m$mu0) - 1)
    qa <- tapply(q$dev[q$analyzer == an], q$day[q$analyzer == an], mean)
    cm <- intersect(names(pm), names(qa))
    ct <- stats::cor.test(pm[cm], qa[cm], method = "spearman", exact = FALSE)
    c(ct$estimate, ct$p.value)
  })

  data.frame(analyte = a$analyte, config = r$config, episodes = ep$n,
             episodes_per_day = ep$per_day, far = mean(run$status),
             alarm_days = length(alarm_days), iqc_days = length(iqc_days),
             iqc_flag_days = length(flag_days),
             flag_days_with_alarm = sum(flag_days %in% alarm_days),
             iqc_absdev_alarm = stats::median(qabs[on_alarm]),
             iqc_absdev_noalarm = stats::median(qabs[!on_alarm]), p_absdev = w,
             direction_concordance = concord,
             rho_a1 = rho[1, 1], p_a1 = rho[2, 1], rho_a2 = rho[1, 2], p_a2 = rho[2, 2])
})
out <- do.call(rbind, out)
utils::write.csv(out, file.path(tab_dir, "real_world.csv"), row.names = FALSE)
print(out, digits = 3, row.names = FALSE)

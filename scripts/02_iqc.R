# Analytical imprecision, sigma metrics and IQC intervals from the 2023 IQC data.

source("R/setup.R")

stream <- readRDS(file.path(data_dir, "stream_2023.rds"))
iqc <- readRDS(file.path(data_dir, "iqc_2023.rds"))
iqc$month <- format(iqc$time, "%m")

rows <- lapply(seq_len(nrow(analytes)), function(i) {
  a <- analytes[i, ]
  q <- iqc[iqc$iqc_name == a$iqc_name, ]
  d <- stream[stream$analyte == a$analyte, ]
  base <- stream_model(d$value, rcv_g(a$cvi, a$cvg))
  tea <- tea_percent(a, base$mean_kept, "primary")

  # Monthly CV for each analyzer module and control level (at least 5 results, non-zero SD)
  g <- split(q$value, list(q$analyzer, q$module, q$level, q$month), drop = TRUE)
  g <- g[lengths(g) >= 5 & vapply(g, stats::sd, numeric(1)) > 0]
  cv <- vapply(g, function(v) 100 * stats::sd(v) / mean(v), numeric(1))
  lev <- sub("^[0-9]+\\.[0-9]+\\.([0-9]+)\\..*$", "\\1", names(cv))

  # Patient results between consecutive IQC events (first IQC result of each day)
  ev <- sort(as.POSIXct(tapply(as.numeric(q$time), as.Date(q$time), min),
                        origin = "1970-01-01", tz = "UTC"))
  between <- as.numeric(table(cut(as.numeric(d$time), breaks = as.numeric(ev))))

  data.frame(analyte = a$analyte, tea_pct = tea, n_iqc = nrow(q), n_groups = length(cv),
             cv_mean = mean(cv), cv_pooled = sqrt(mean(cv^2)),
             cv_level1 = sqrt(mean(cv[lev == "1"]^2)), cv_level2 = sqrt(mean(cv[lev == "2"]^2)),
             sigma_mean = mean(tea / cv), sigma_pooled = tea / sqrt(mean(cv^2)),
             sigma_min = min(tea / cv),
             iqc_days = length(ev), between_median = stats::median(between),
             between_q1 = stats::quantile(between, 0.25, names = FALSE),
             between_q3 = stats::quantile(between, 0.75, names = FALSE),
             b_iqc = 10 * round(stats::median(between) / 10))
})
iqc_cv <- do.call(rbind, rows)
saveRDS(iqc_cv, file.path(out_dir, "iqc_cv.rds"))
utils::write.csv(iqc_cv, file.path(out_dir, "iqc_summary.csv"), row.names = FALSE)
print(iqc_cv[, c("analyte", "tea_pct", "cv_pooled", "sigma_mean", "sigma_pooled", "sigma_min", "b_iqc")],
      digits = 3, row.names = FALSE)

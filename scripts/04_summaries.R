# Summary tables for the primary analysis and the sensitivity analyses.

source("R/setup.R")

tab_dir <- file.path(out_dir, "tables")
dir.create(tab_dir, showWarnings = FALSE)
max_episodes <- 1
grid <- chart_grid()

iqc_cv <- readRDS(file.path(out_dir, "iqc_cv.rds"))
b_iqc <- setNames(iqc_cv$b_iqc, iqc_cv$analyte)
load_scenario <- function(s) readRDS(file.path(out_dir, sprintf("scenario_%s.rds", s)))
prim <- load_scenario("primary")

select <- function(res, B, max_ep = max_episodes) syi_table(res, B, grid, max_episodes = max_ep)

# Primary analysis --------------------------------------------------------------

syi_all <- list(); rec <- list(); yt_all <- list(); dl_all <- list()
for (k in analytes$analyte) {
  res <- prim[[k]]
  B <- b_iqc[[k]]
  con <- select(res, B)
  unc <- select(res, B, Inf)
  con$analyte <- k
  con$best_unconstrained <- con$config == unc$config[unc$best]
  syi_all[[k]] <- con
  b1 <- con[con$best, ]
  b0 <- unc[unc$best, ]
  eq <- con$config[con$equivalent & !con$best]
  rec[[k]] <- data.frame(
    analyte = k, tea_pct = res$tea, B = B,
    config = b1$config, method = b1$method, lambda = b1$lambda, L = b1$L, w = b1$w,
    syi = b1$syi, syi_lo = b1$lo, syi_hi = b1$hi, far = b1$far,
    episodes_per_day = b1$episodes_per_day, n_equivalent = length(eq),
    equivalent = paste(eq, collapse = "; "),
    unc_config = b0$config, unc_syi = b0$syi, unc_far = b0$far,
    unc_episodes_per_day = b0$episodes_per_day,
    n_eligible = sum(con$eligible))
  yt <- youden_table(res, B); yt$analyte <- k; yt_all[[k]] <- yt
  dl <- delay_table(res, B); dl$analyte <- k; dl_all[[k]] <- dl
}
syi_all <- do.call(rbind, syi_all)
rec <- do.call(rbind, rec)
yt_all <- do.call(rbind, yt_all)
dl_all <- do.call(rbind, dl_all)
write.csv(syi_all, file.path(tab_dir, "syi_primary_all.csv"), row.names = FALSE)
write.csv(rec, file.path(tab_dir, "recommended.csv"), row.names = FALSE)
write.csv(yt_all, file.path(tab_dir, "youden_primary.csv"), row.names = FALSE)
write.csv(dl_all, file.path(tab_dir, "delay_primary.csv"), row.names = FALSE)

# Truncation and stream description ------------------------------------------------

stream <- readRDS(file.path(data_dir, "stream_2023.rds"))
trunc <- do.call(rbind, lapply(seq_len(nrow(analytes)), function(i) {
  a <- analytes[i, ]
  x <- stream$value[stream$analyte == a$analyte]
  tm <- stream$time[stream$analyte == a$analyte]
  m <- prim[[a$analyte]]$model
  kept <- x[inside_limits(x, m$lim)]
  f <- a$si_factor
  data.frame(analyte = a$analyte, rcv = rcv_g(a$cvi, a$cvg), n = length(x),
             mean = f * mean(x), sd = f * stats::sd(x), median = f * stats::median(x),
             n_kept = length(kept), mean_kept = f * mean(kept), sd_kept = f * stats::sd(kept),
             median_kept = f * stats::median(kept), excluded_pct = 100 * m$excluded,
             ltl = f * m$lim[1], utl = f * m$lim[2], offset = m$offset, mu0 = m$mu0, s0 = m$s0,
             days = length(unique(as.Date(tm))),
             per_day_median = stats::median(table(as.Date(tm[inside_limits(x, m$lim)]))),
             tea_pct = prim[[a$analyte]]$tea)
}))
write.csv(trunc, file.path(tab_dir, "truncation.csv"), row.names = FALSE)

# Block length --------------------------------------------------------------------

bsens <- do.call(rbind, lapply(analytes$analyte, function(k) {
  do.call(rbind, lapply(c(100, 250, 500, 1000), function(B) {
    st <- select(prim[[k]], B)
    b <- st[st$best, ]
    r <- rec[rec$analyte == k, ]
    data.frame(analyte = k, B = B, config = b$config, syi = b$syi,
               rec_config = r$config, rec_syi = st$syi[st$config == r$config],
               rec_equivalent = st$equivalent[st$config == r$config])
  }))
}))
write.csv(bsens, file.path(tab_dir, "b_sensitivity.csv"), row.names = FALSE)

# Sensitivity analyses -----------------------------------------------------------

compare_scenario <- function(s) {
  sc <- load_scenario(s)
  do.call(rbind, lapply(analytes$analyte, function(k) {
    st <- select(sc[[k]], b_iqc[[k]])
    b <- st[st$best, ]
    r <- rec[rec$analyte == k, ]
    p <- st[st$config == r$config, ]
    data.frame(scenario = s, analyte = k, tea_pct = sc[[k]]$tea, n = sc[[k]]$n,
               excluded_pct = 100 * sc[[k]]$model$excluded,
               config = b$config, syi = b$syi, episodes_per_day = b$episodes_per_day,
               rec_config = r$config, rec_syi = p$syi, rec_episodes_per_day = p$episodes_per_day,
               rec_eligible = p$eligible, rec_equivalent = p$equivalent | p$best,
               same = b$config == r$config, rec_rank = which(st$config == r$config))
  }))
}
scen <- c("stream_a1", "stream_a2", "derivation", "validation", "cva_measured", "median",
          "logsym", "tea_greater", "tea_eflm", "permuted", "analyzer1")
scen <- scen[file.exists(file.path(out_dir, sprintf("scenario_%s.rds", scen)))]
cmp <- do.call(rbind, lapply(scen, compare_scenario))
write.csv(cmp, file.path(tab_dir, "scenarios.csv"), row.names = FALSE)

# Temporal validation: configuration selected in January-June, evaluated in
# July-December with the January-June parameters.
if (all(c("derivation", "validation") %in% scen)) {
  der <- load_scenario("derivation"); val <- load_scenario("validation")
  tv <- do.call(rbind, lapply(analytes$analyte, function(k) {
    sd1 <- select(der[[k]], b_iqc[[k]]); sv <- select(val[[k]], b_iqc[[k]])
    c1 <- sd1$config[sd1$best]
    v <- sv[sv$config == c1, ]
    data.frame(analyte = k, derivation_config = c1, derivation_syi = sd1$syi[sd1$best],
               derivation_episodes = sd1$episodes_per_day[sd1$best],
               validation_syi = v$syi, validation_episodes = v$episodes_per_day,
               validation_best = sv$config[sv$best], validation_best_syi = sv$syi[sv$best],
               within_equivalence = v$equivalent | v$best, validation_eligible = v$eligible,
               rank = which(sv$config == c1),
               rho = stats::cor(sd1$syi[order(sd1$config)], sv$syi[order(sv$config)],
                                method = "spearman"))
  }))
  write.csv(tv, file.path(tab_dir, "temporal_validation.csv"), row.names = FALSE)

  # Change of the in-control mean between the half-years
  split_date <- as.POSIXct("2023-07-01", tz = "UTC")
  shift <- do.call(rbind, lapply(analytes$analyte, function(k) {
    d <- stream[stream$analyte == k & stream$time >= split_date, ]
    m <- der[[k]]$model
    y <- log(d$value[inside_limits(d$value, m$lim)] + m$offset)
    data.frame(analyte = k, shift_pct = 100 * (exp(mean(y) - m$mu0) - 1), sd_ratio = stats::sd(y) / m$s0)
  }))
  write.csv(shift, file.path(tab_dir, "halfyear_shift.csv"), row.names = FALSE)
}

# Error confined to one analyzer, recommended configurations.
if ("analyzer1" %in% scen) {
  a1 <- load_scenario("analyzer1")
  one <- do.call(rbind, lapply(analytes$analyte, function(k) {
    r <- rec$config[rec$analyte == k]
    d0 <- dl_all[dl_all$analyte == k & dl_all$config == r, ]
    d1 <- delay_table(a1[[k]], b_iqc[[k]]); d1 <- d1[d1$config == r, ]
    s0 <- syi_all$syi[syi_all$analyte == k & syi_all$config == r]
    s1 <- select(a1[[k]], b_iqc[[k]]); s1 <- s1$syi[s1$config == r]
    data.frame(analyte = k, config = r, syi_both = s0, syi_one = s1,
               mnped_both_1tea = mean(d0$mnped[abs(d0$multiple) == 1]),
               mnped_one_1tea = mean(d1$mnped[abs(d1$multiple) == 1]),
               pdet_both_1tea = mean(d0$p_detect[abs(d0$multiple) == 1]),
               pdet_one_1tea = mean(d1$p_detect[abs(d1$multiple) == 1]))
  }))
  write.csv(one, file.path(tab_dir, "single_analyzer_error.csv"), row.names = FALSE)
}

print(rec[, c("analyte", "tea_pct", "B", "config", "syi", "episodes_per_day", "n_equivalent",
              "unc_config", "unc_syi", "unc_episodes_per_day", "n_eligible")], digits = 3, row.names = FALSE)

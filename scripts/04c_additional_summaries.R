# Summaries of the additional analyses: alarm-burden limit, monthly
# re-estimation of mu0 and sigma0, and grid extension.

source("R/setup.R")

tab_dir <- file.path(out_dir, "tables")
grid <- chart_grid()
iqc_cv <- readRDS(file.path(out_dir, "iqc_cv.rds"))
b_iqc <- setNames(iqc_cv$b_iqc, iqc_cv$analyte)
prim <- readRDS(file.path(out_dir, "scenario_primary.rds"))
rec <- read.csv(file.path(tab_dir, "recommended.csv"))

# Alarm-burden limit ------------------------------------------------------------

lim <- do.call(rbind, lapply(analytes$analyte, function(k) {
  r <- rec$config[rec$analyte == k]
  do.call(rbind, lapply(c(0.5, 1, 2), function(L) {
    st <- syi_table(prim[[k]], b_iqc[[k]], grid, max_episodes = L)
    b <- st[st$best, ]
    p <- st[st$config == r, ]
    data.frame(analyte = k, limit = L, config = b$config, syi = b$syi,
               episodes_per_day = b$episodes_per_day,
               rec_eligible = p$eligible, rec_equivalent = p$best | p$equivalent)
  }))
}))
write.csv(lim, file.path(tab_dir, "alarm_limit.csv"), row.names = FALSE)
print(aggregate(cbind(syi, episodes_per_day) ~ limit, lim, sum))

# Monthly re-estimation ---------------------------------------------------------

f <- file.path(out_dir, "scenario_recalibration.rds")
if (file.exists(f)) {
  rc <- readRDS(f)
  tv <- read.csv(file.path(tab_dir, "temporal_validation.csv"))
  recal <- do.call(rbind, lapply(analytes$analyte, function(k) {
    c1 <- tv$derivation_config[tv$analyte == k]
    st <- syi_table(rc[[k]], b_iqc[[k]], grid, max_episodes = 1)
    v <- st[st$config == c1, ]
    data.frame(analyte = k, derivation_config = c1,
               fixed_episodes = tv$validation_episodes[tv$analyte == k],
               fixed_syi = tv$validation_syi[tv$analyte == k],
               fixed_within = tv$within_equivalence[tv$analyte == k],
               recal_episodes = v$episodes_per_day, recal_syi = v$syi,
               recal_eligible = v$eligible, recal_within = v$best | v$equivalent,
               recal_best = st$config[st$best], recal_best_syi = st$syi[st$best],
               n_eligible = sum(st$eligible))
  }))
  write.csv(recal, file.path(tab_dir, "recalibration.csv"), row.names = FALSE)
  print(recal, digits = 3, row.names = FALSE)
}

# Grid extension ------------------------------------------------------------------

f <- file.path(out_dir, "scenario_grid_extension.rds")
if (file.exists(f)) {
  ge <- readRDS(f)
  g2 <- rbind(grid, ge$grid)
  ext <- do.call(rbind, lapply(analytes$analyte, function(k) {
    res <- prim[[k]]
    res$ic <- rbind(res$ic, ge$res[[k]]$ic[, names(res$ic)])
    res$sim <- rbind(res$sim, ge$res[[k]]$sim[, names(res$sim)])
    st <- syi_table(res, b_iqc[[k]], g2, max_episodes = 1)
    r <- rec$config[rec$analyte == k]
    p <- st[st$config == r, ]
    b <- st[st$best, ]
    data.frame(analyte = k, rec_config = r, rec_syi = p$syi, ext_best = b$config,
               ext_best_syi = b$syi, ext_best_episodes = b$episodes_per_day,
               best_is_extension = b$config %in% ge$grid$config,
               rec_within = p$best | p$equivalent,
               min_in_control_extra = min(ge$res[[k]]$ic$in_control_at_onsets))
  }))
  write.csv(ext, file.path(tab_dir, "grid_extension.csv"), row.names = FALSE)
  print(ext, digits = 3, row.names = FALSE)
}

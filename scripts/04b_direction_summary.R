# Summary of the direction-specific sensitivity analysis and of chance signals.

source("R/setup.R")

tab_dir <- file.path(out_dir, "tables")
grid <- chart_grid()
iqc_cv <- readRDS(file.path(out_dir, "iqc_cv.rds"))
b_iqc <- setNames(iqc_cv$b_iqc, iqc_cv$analyte)
dir <- readRDS(file.path(out_dir, "scenario_direction.rds"))
rec <- read.csv(file.path(tab_dir, "recommended.csv"))

cmp <- do.call(rbind, lapply(analytes$analyte, function(k) {
  B <- b_iqc[[k]]
  st <- syi_table(dir[[k]], B, grid, max_episodes = 1)
  b <- st[st$best, ]
  r <- rec[rec$analyte == k, ]
  p <- st[st$config == r$config, ]
  nl <- dir[[k]]$null
  nl <- nl[nl$config == r$config, ]
  data.frame(analyte = k, B = B, rec_config = r$config, rec_syi_any = r$syi,
             rec_syi_dir = p$syi, dir_config = b$config, dir_syi = b$syi, same = b$config == r$config,
             rec_equivalent = p$equivalent | p$best,
             null_p_alarm = mean(!is.na(nl$any) & nl$any <= B))
}))
write.csv(cmp, file.path(tab_dir, "direction_sensitivity.csv"), row.names = FALSE)

# Chance signals for every configuration: probability of at least one alarm
# within B results after the onsets without an error.
null_all <- do.call(rbind, lapply(analytes$analyte, function(k) {
  B <- b_iqc[[k]]
  nl <- dir[[k]]$null
  a <- stats::aggregate(any ~ config, transform(nl, any = !is.na(any) & any <= B), mean)
  names(a)[2] <- "p_alarm_no_error"
  a$analyte <- k
  a
}))
write.csv(null_all, file.path(tab_dir, "null_alarm.csv"), row.names = FALSE)

print(cmp, digits = 3, row.names = FALSE)
cat("same:", sum(cmp$same), " equivalent:", sum(cmp$same | cmp$rec_equivalent), "\n")

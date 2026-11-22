# Figure 2 and interactive control charts for the recommended configurations.

source("R/setup.R")
source("R/plots.R")
library(patchwork)

tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
html_dir <- file.path(out_dir, "interactive")
dir.create(html_dir, showWarnings = FALSE)

stream <- readRDS(file.path(data_dir, "stream_2023.rds"))
iqc <- readRDS(file.path(data_dir, "iqc_2023.rds"))
prim <- readRDS(file.path(out_dir, "scenario_primary.rds"))
rec <- utils::read.csv(file.path(tab_dir, "recommended.csv"), encoding = "UTF-8")
grid <- chart_grid()

chart_data <- function(k) {
  a <- analytes[analytes$analyte == k, ]
  d <- stream[stream$analyte == k, ]
  m <- prim[[k]]$model
  cfg <- grid[grid$config == rec$config[rec$analyte == k], ]
  run <- run_chart(d$value, m, cfg)
  ins <- inside_limits(d$value, m$lim)
  n <- seq_len(run$n_mon)
  hw <- if (cfg$method == "EWMA") {
    cfg$L * m$s0 * sqrt(cfg$lambda / (2 - cfg$lambda) * (1 - (1 - cfg$lambda)^(2 * n)))
  } else cfg$L * m$s0 / sqrt(pmin(n, cfg$w))
  data.frame(time = d$time[ins], raw = d$value[ins] * a$si_factor, analyzer = d$analyzer[ins],
             dev = 100 * (exp(run$stat - m$mu0) - 1),
             ucl = 100 * (exp(hw) - 1), lcl = 100 * (exp(-hw) - 1),
             alarm = abs(run$stat - m$mu0) > hw, config = cfg$config, si_unit = a$si_unit)
}

iqc_dev <- function(k) {
  a <- analytes[analytes$analyte == k, ]
  q <- iqc[iqc$iqc_name == a$iqc_name, ]
  q$dev <- 100 * (q$value / stats::ave(q$value, q$analyzer, q$module, q$level, q$lot,
                                       FUN = stats::median) - 1)
  q
}

# Figure 2: one working week of the creatinine stream with the IQC record -------------

k <- "creatinine"
cd <- chart_data(k)
q <- iqc_dev(k)
from <- as.POSIXct("2023-03-13", tz = "UTC"); to <- as.POSIXct("2023-03-18", tz = "UTC")
w <- cd[cd$time >= from & cd$time < to, ]
w$i <- seq_len(nrow(w))
w$day <- as.Date(w$time)
day_start <- tapply(w$i, w$day, min)
qw <- q[q$time >= from & q$time < to, ]
qw$i <- vapply(qw$time, function(t) { j <- which(w$time >= t)[1]; if (is.na(j)) NA_integer_ else w$i[j] }, integer(1))
qw$series <- sprintf("AU5800-%d%s, level %d", qw$analyzer,
                     ifelse(qw$analyzer == 1, sprintf(" module %d", qw$module), ""), qw$level)
iqc_i <- unique(tapply(qw$i, as.Date(qw$time), min))
brk <- as.numeric(day_start)
lab <- format(as.Date(names(day_start)), "%d %b")

pa <- ggplot(w, aes(i, dev)) +
  geom_vline(xintercept = iqc_i, linetype = "dotted", colour = "grey45") +
  geom_line(aes(y = ucl), colour = "#E8112D", linewidth = 0.4) +
  geom_line(aes(y = lcl), colour = "#E8112D", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_line(colour = "#08519c", linewidth = 0.3) +
  geom_point(data = w[w$alarm, ], colour = "#E8112D", size = 0.8) +
  scale_x_continuous(breaks = brk, labels = lab) +
  labs(x = "Patient results in order of reporting (labels mark the first result of each day)",
       y = "Control statistic
(% from in-control mean)",
       subtitle = sprintf("A  Creatinine, %s; dotted lines: daily IQC", cd$config[1])) +
  theme_pub(8)
pb <- ggplot(qw[!is.na(qw$i), ], aes(i, dev, colour = series)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_point(size = 1.2) + geom_line(linewidth = 0.3) +
  scale_x_continuous(breaks = brk, labels = lab, limits = range(w$i)) +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(x = NULL, y = "IQC deviation (%)", subtitle = "B  IQC results, deviation from the lot median") +
  theme_pub(8) + theme(legend.position = "bottom", legend.text = element_text(size = 6.5)) +
  guides(colour = guide_legend(nrow = 2))
p2 <- pa / pb + plot_layout(heights = c(1.4, 1))
save_figure(p2, file.path(fig_dir, "Figure2"), width = 17.4, height = 13)

# Data for the interactive full-year charts (scripts/09_interactive_html.py).
# Only the control statistic is exported (rounded, time to the minute); individual
# patient results are not included.

for (k in c("creatinine", "sodium")) {
  cd <- chart_data(k)
  out <- data.frame(time = format(trunc(cd$time, "mins"), "%Y-%m-%d %H:%M"), stat = round(cd$dev, 2),
                    alarm = as.integer(cd$alarm))
  utils::write.csv(out, file.path(html_dir, sprintf("chart_%s.csv", k)), row.names = FALSE)
  writeLines(c(analytes$label[analytes$analyte == k], cd$config[1],
               sprintf("%.4f", utils::tail(cd$ucl, 1)), sprintf("%.4f", utils::tail(cd$lcl, 1))),
             file.path(html_dir, sprintf("chart_%s.txt", k)))
}

# Figure 1: study workflow.

source("R/setup.R")
source("R/plots.R")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

boxes <- data.frame(
  id = 1:8,
  x = c(1, 1, 1, 1, 1, 1, 1, 1),
  y = c(8, 7, 6, 5, 4, 3, 2, 1),
  title = c("Data",
            "Eligibility and quality goals",
            "Truncation and transformation",
            "Candidate charts (53 configurations)",
            "In-control performance",
            "Error simulation",
            "Youden index and SYI",
            "Selection and robustness"),
  body = c("LIS, 2023: 23 analytes, two AU5800 analyzers, adults, non-acute units\n2,876,239 results kept in chronological order",
           "Sigma metric from IQC (pooled CV) ≥ 2; CLIA 2024 TEa; EFLM CVi and CVg",
           "RCVg = √2 × 1.96 × √(CVg² + CVi² + CVa²), CVa = 0.5 CVi; limits around the annual mean\nNatural logarithm of retained results",
           "EWMA: λ 0.02–0.20 × L 2.5–5.0 (28)    MA: w 10–100 × k 3.0–5.0 (25)",
           "Each chart run over the full year: false-alarm proportion and alarm episodes per day",
           "100 error onsets spread across the year; proportional bias of ±0.25 to ±3 × TEa added to raw results\nbefore truncation; follow-up of B results (B = median number of results between daily IQC events)",
           "Se = (B − d + 1)/B, d = first signal; J = Se + Sp − 1\nSYI = sum of J over |bias| ≤ TEa (8 levels); bootstrap 95% CI over onsets",
           "Highest SYI with ≤ 1 alarm episode per day; equivalent alternatives\nB, analyzer, half-year validation, truncation, TEa and order sensitivity; IQC concordance"))

p <- ggplot(boxes) +
  geom_segment(data = boxes[-nrow(boxes), ], aes(x = x, xend = x, y = y - 0.33, yend = y - 0.62),
               arrow = arrow(length = unit(1.6, "mm"), type = "closed"), colour = "grey30") +
  geom_label(aes(x = x, y = y, label = paste0(title, "\n", body)), size = 2.6, hjust = 0.5,
             label.padding = unit(1.6, "mm"), label.r = unit(1, "mm"), fill = "grey97",
             lineheight = 1.05, label.size = 0.3) +
  scale_x_continuous(limits = c(0.4, 1.6)) +
  scale_y_continuous(limits = c(0.5, 8.5)) +
  theme_void()

save_figure(p, file.path(fig_dir, "Figure1"), width = 17, height = 19)

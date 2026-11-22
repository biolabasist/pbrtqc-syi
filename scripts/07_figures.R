# Figures 2-5 and supplementary figures.

source("R/setup.R")
source("R/plots.R")
library(patchwork)

tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

grid <- chart_grid()
prim <- readRDS(file.path(out_dir, "scenario_primary.rds"))
syi <- utils::read.csv(file.path(tab_dir, "syi_primary_all.csv"), encoding = "UTF-8")
rec <- utils::read.csv(file.path(tab_dir, "recommended.csv"), encoding = "UTF-8")
yt <- utils::read.csv(file.path(tab_dir, "youden_primary.csv"), encoding = "UTF-8")
iqc_cv <- readRDS(file.path(out_dir, "iqc_cv.rds"))
b_iqc <- setNames(iqc_cv$b_iqc, iqc_cv$analyte)
labels <- setNames(analytes$label, analytes$analyte)
short <- setNames(analytes$label, analytes$analyte)
short[c("alt", "alp", "ast", "ggt", "ldh")] <- c("ALT", "ALP", "AST", "GGT", "LDH")
mult_lab <- function(m) sub("-", "−", sprintf("%g", m))

row_label <- function(method, lambda, w) {
  ifelse(method == "EWMA", sprintf("λ %.2f", lambda), sprintf("w %d", as.integer(w)))
}
row_levels <- c(sprintf("λ %.2f", c(0.02, 0.05, 0.10, 0.20)), sprintf("w %d", c(10, 25, 50, 75, 100)))
mult_levels <- c("2.5", "2.7", "3.0", "3.5", "4.0", "4.5", "5.0")

# Figure 3: SYI landscape for all analytes -----------------------------------------

d3 <- syi
d3$row <- factor(row_label(d3$method, d3$lambda, d3$w), levels = rev(row_levels))
d3$col <- factor(sprintf("%.1f", d3$L), levels = mult_levels)
d3$mj <- d3$syi / 8
d3$panel <- factor(short[d3$analyte], levels = sort(unname(short)))
p3 <- ggplot(d3, aes(col, row)) +
  geom_tile(aes(fill = mj), colour = "white", linewidth = 0.25) +
  geom_point(data = d3[!d3$eligible, ], shape = 4, size = 0.9, colour = "grey85", stroke = 0.35) +
  geom_point(data = d3[d3$best_unconstrained, ], shape = 1, size = 2.1, colour = "white", stroke = 0.7) +
  geom_point(data = d3[d3$best, ], shape = 8, size = 1.9, colour = "#E8112D", stroke = 0.8) +
  geom_hline(yintercept = 5.5, colour = "white", linewidth = 0.8) +
  scale_fill_viridis_c(option = "viridis", limits = c(0, 1), oob = scales::squish,
                       name = "Mean Youden index within ±TEa (SYI/8)") +
  facet_wrap(~panel, ncol = 5) +
  labs(x = "Control limit multiplier (EWMA L; MA k)", y = NULL) +
  theme_pub(8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 6), panel.grid = element_blank(),
        legend.key.width = unit(1.4, "cm"), legend.key.height = unit(0.25, "cm"),
        strip.text = element_text(size = 6.5))
save_figure(p3, file.path(fig_dir, "Figure3"), width = 17.4, height = 21)

# Analyte detail panels ----------------------------------------------------------

detail_panels <- function(k) {
  r <- rec[rec$analyte == k, ]
  res <- prim[[k]]
  B <- b_iqc[[k]]
  s <- syi[syi$analyte == k, ]
  e1 <- s[s$method == "EWMA" & s$eligible, ]
  m1 <- s[s$method == "MA" & s$eligible, ]
  Lstar <- if (r$method == "EWMA") r$L else e1$L[which.max(e1$syi)]
  kstar <- if (r$method == "MA") r$L else m1$L[which.max(m1$syi)]

  # (A) Youden index by configuration and bias
  y <- yt[yt$analyte == k, ]
  y <- merge(y, grid, by = "config")
  y <- y[(y$method == "EWMA" & y$L == Lstar) | (y$method == "MA" & y$L == kstar), ]
  y$row <- factor(row_label(y$method, y$lambda, y$w), levels = rev(row_levels))
  y$colm <- factor(mult_lab(y$multiple), levels = mult_lab(sort(unique(y$multiple))))
  y$fill <- fill_colours(y$J, c(0, 1))
  y$txt <- contrast_text(y$fill)
  y$lab <- sprintf("%.2f", y$J)
  inner <- which(abs(sort(unique(y$multiple))) <= 1)
  pa <- ggplot(y, aes(colm, row)) +
    geom_tile(aes(fill = J), colour = "white", linewidth = 0.3) +
    geom_text(aes(label = lab, colour = txt), size = 2.3) +
    scale_colour_identity() +
    annotate("rect", xmin = min(inner) - 0.5, xmax = max(inner) + 0.5, ymin = 0.5,
             ymax = length(levels(droplevels(y$row))) + 0.5, fill = NA, colour = "#E8112D",
             linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = 5.5, colour = "white", linewidth = 1) +
    scale_fill_viridis_c(limits = c(0, 1), oob = scales::squish, name = "J") +
    labs(x = "Bias (multiples of TEa)", y = NULL,
         subtitle = sprintf("A  Youden index (EWMA L = %.1f; MA k = %.1f)", Lstar, kstar)) +
    theme_pub(8) + theme(panel.grid = element_blank(), legend.position = "right",
                         legend.key.height = unit(0.5, "cm"))

  # (B) Results to first signal for the recommended configuration
  sm <- res$sim[res$sim$config == r$config, ]
  sm$det <- !is.na(sm$first) & sm$first <= B
  sm$colm <- factor(mult_lab(sm$multiple), levels = mult_lab(sort(unique(sm$multiple))))
  pd <- stats::aggregate(det ~ colm, sm, mean)
  ymax <- B
  pb <- ggplot(sm[sm$det, ], aes(colm, first)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.3, fill = "#9ecae1") +
    geom_text(data = pd, aes(colm, ymax * 1.6, label = sprintf("%.0f", 100 * det)), size = 2,
              inherit.aes = FALSE) +
    scale_y_log10(limits = c(1, ymax * 2)) +
    scale_x_discrete(drop = FALSE) +
    labs(x = "Bias (multiples of TEa)", y = "Results to first signal",
         subtitle = sprintf("B  %s: results to first signal (B = %d)", r$config, B)) +
    theme_pub(8) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  # (C) Probability of detection within B
  pw <- res$sim
  pw <- merge(pw, grid, by = "config")
  pw <- pw[(pw$method == "EWMA" & pw$L == Lstar) | (pw$method == "MA" & pw$L == kstar), ]
  pw$det <- !is.na(pw$first) & pw$first <= B
  pw <- stats::aggregate(det ~ config + method + multiple, pw, mean)
  pw$rec <- pw$config == r$config
  pw$config <- factor(pw$config, levels = grid$config[grid$config %in% pw$config])
  pc <- ggplot(pw, aes(multiple, 100 * det, group = config, colour = config)) +
    geom_line(aes(linewidth = rec)) + geom_point(size = 0.6) +
    scale_linewidth_manual(values = c(`FALSE` = 0.35, `TRUE` = 1.1), guide = "none") +
    scale_colour_viridis_d(option = "turbo", name = NULL) +
    scale_x_continuous(breaks = bias_multiples, labels = mult_lab(bias_multiples)) +
    labs(x = "Bias (multiples of TEa)", y = "Detected within B (%)",
         subtitle = "D  Power function") +
    theme_pub(8) + theme(legend.position = "right", legend.key.size = unit(0.3, "cm"),
                         legend.text = element_text(size = 6.5),
                         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  # (D) SYI against alarm burden
  s$tag <- ifelse(s$best, "Recommended", ifelse(s$best_unconstrained, "Highest SYI", "Other"))
  pdd <- ggplot(s, aes(pmax(episodes_per_day, 0.01), syi)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_point(aes(colour = method, shape = tag, size = tag)) +
    scale_shape_manual(values = c(Recommended = 8, `Highest SYI` = 1, Other = 16), name = NULL) +
    scale_size_manual(values = c(Recommended = 2.6, `Highest SYI` = 2.6, Other = 1.1), name = NULL) +
    scale_colour_manual(values = c(EWMA = "#3182bd", MA = "#e6550d"), name = NULL) +
    scale_x_log10() +
    labs(x = "Alarm episodes per day (2023 stream)", y = "SYI", subtitle = "C  Detection versus alarm burden") +
    theme_pub(8) + theme(legend.position = "right")

  pa / (pb | pdd) / pc + plot_layout(heights = c(1.15, 1, 0.9)) +
    plot_annotation(title = labels[[k]], theme = theme(plot.title = element_text(face = "bold", size = 10)))
}

save_figure(detail_panels("creatinine"), file.path(fig_dir, "Figure4"), width = 17.4, height = 22)
save_figure(detail_panels("sodium"), file.path(fig_dir, "Figure5"), width = 17.4, height = 22)

supp_dir <- file.path(fig_dir, "supplementary")
dir.create(supp_dir, showWarnings = FALSE)
ord <- analytes$analyte[order(short[analytes$analyte])]
for (i in seq_along(ord)) {
  k <- ord[i]
  ggsave(file.path(supp_dir, sprintf("FigureS%02d_%s.png", i, k)), detail_panels(k),
         width = 17.4, height = 22, units = "cm", dpi = 300, bg = "white")
}

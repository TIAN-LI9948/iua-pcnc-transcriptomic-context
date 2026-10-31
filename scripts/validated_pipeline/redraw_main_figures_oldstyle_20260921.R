#!/usr/bin/env Rscript

# Restyles the locked, reframed main figures using the visual language of the
# earlier submission figures. This script performs no statistical reanalysis:
# all estimates, intervals, q values, eligibility decisions, and cohort sizes
# are read from the locked figure source-data tables.

options(stringsAsFactors = FALSE, warn = 1)

required <- c("ggplot2", "cowplot", "ragg", "svglite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 1L) {
  stop("Usage: Rscript redraw_main_figures_oldstyle_20260921.R <output_directory>", call. = FALSE)
}

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
script_file <- normalizePath(script_arg[[1L]], winslash = "/", mustWork = TRUE)
base_root <- normalizePath(dirname(script_file), winslash = "/", mustWork = TRUE)
source_root <- normalizePath(file.path(base_root, "revised_r2", "source_data"), winslash = "/", mustWork = TRUE)
output_root <- normalizePath(argv[[1L]], winslash = "/", mustWork = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)

read_tsv <- function(name) {
  path <- file.path(source_root, name)
  if (!file.exists(path)) stop("Missing locked source table: ", path, call. = FALSE)
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

donor <- read_tsv("Figure_1_donor_pseudobulks.tsv")
estimates <- read_tsv("Figure_1_locked_endpoint_estimates.tsv")
effects <- read_tsv("Figure_2_cross_context_effect_vectors.tsv")
concordance <- read_tsv("Figure_2_concordance_summary.tsv")
bulk <- read_tsv("Figure_3_independent_bulk_endpoint_tests.tsv")

# Frozen source contract. Fail closed if a staging table no longer matches the
# source-to-figure audit that underlies the current manuscript.
stopifnot(
  nrow(donor) == 16L,
  sum(donor$CycleType == "Controlled") == 7L,
  sum(donor$CycleType == "Natural") == 9L,
  nrow(estimates) == 4L,
  nrow(effects) == 10465L,
  nrow(concordance) == 1L,
  concordance$genes[[1L]] == 10465L,
  nrow(bulk) == 4L
)

for (stem in c("Figure_1", "Figure_2", "Figure_3")) {
  targets <- file.path(output_root, paste0(stem, c(".pdf", ".png", ".svg", ".tiff")))
  if (any(file.exists(targets))) {
    stop("Refusing to overwrite an existing preview bundle for ", stem, call. = FALSE)
  }
}

# Style-only inheritance from the earlier figures: orange circles and blue
# squares for conditions, classic black axes, restrained grey intervals, and
# large uppercase panel labels. Scientific mappings are defined afresh below.
PC_ORANGE <- "#D55E00"
NC_BLUE <- "#0072B2"
PGR_ORANGE <- "#D55E00"
ESR1_BLUE <- "#0072B2"
WOI_PURPLE <- "#6A5ACD"
CORE_GREEN <- "#009E73"
NEUTRAL_DARK <- "#333333"
NEUTRAL_MID <- "#7A7A7A"
NEUTRAL_LIGHT <- "#D9D9D9"

theme_submission <- function(base_size = 9.2) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.45, colour = "black"),
      axis.ticks.length = grid::unit(1.35, "mm"),
      axis.text = ggplot2::element_text(size = base_size - 0.55, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size, colour = "black"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(7, 6, 5, 7)
    )
}

format_q_expr <- function(x) {
  if (!is.finite(x)) return("italic(q)~'not estimated'")
  if (x < 0.001) return("italic(q) < 0.001")
  paste0("italic(q) == ", formatC(x, format = "f", digits = 4L))
}

compose_panels <- function(..., labels, ncol = 2L, rel_widths = 1, rel_heights = 1) {
  cowplot::plot_grid(
    ...,
    labels = labels,
    ncol = ncol,
    rel_widths = rel_widths,
    rel_heights = rel_heights,
    align = "none",
    label_x = 0.012,
    label_y = 0.995,
    hjust = 0,
    vjust = 1,
    label_size = 13,
    label_fontface = "bold",
    label_fontfamily = "Helvetica"
  )
}

save_preview_bundle <- function(plot, stem, width_mm, height_mm) {
  pdf_path <- file.path(output_root, paste0(stem, ".pdf"))
  png_path <- file.path(output_root, paste0(stem, ".png"))
  svg_path <- file.path(output_root, paste0(stem, ".svg"))
  tiff_path <- file.path(output_root, paste0(stem, ".tiff"))
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in,
                       family = "Helvetica", bg = "white")
  print(plot)
  grDevices::dev.off()

  svglite::svglite(svg_path, width = width_in, height = height_in, bg = "white")
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(png_path, width = width_in, height = height_in, units = "in",
                res = 600, background = "white", scaling = 1)
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(tiff_path, width = width_in, height = height_in, units = "in",
                 res = 600, background = "white", compression = "lzw")
  print(plot)
  grDevices::dev.off()
}

# Figure 1: donor-level PC/NC endpoints, rendered as the small-n box-and-point
# panels used in the earlier figures. Points remain the biological donors and
# the thick black crossbar is the group mean.
donor$context <- factor(
  ifelse(donor$CycleType == "Controlled", "Programmed cycle", "Natural cycle"),
  levels = c("Programmed cycle", "Natural cycle")
)

endpoint_q <- setNames(estimates$exact_BH_q_within_family, estimates$endpoint)

small_n_panel <- function(value_column, y_label, seed) {
  d <- data.frame(context = donor$context, value = donor[[value_column]])
  means <- stats::aggregate(value ~ context, data = d, FUN = mean)
  y_span <- diff(range(d$value, na.rm = TRUE))
  if (!is.finite(y_span) || y_span <= 0) y_span <- 1
  y_annot <- max(d$value, na.rm = TRUE) + 0.15 * y_span

  ggplot2::ggplot(d, ggplot2::aes(x = context, y = value, fill = context, shape = context)) +
    ggplot2::geom_boxplot(
      width = 0.45, outlier.shape = NA, alpha = 0.15,
      colour = NEUTRAL_DARK, linewidth = 0.42, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.075, height = 0, seed = seed),
      size = 2.45, colour = NEUTRAL_DARK, stroke = 0.45, show.legend = FALSE
    ) +
    ggplot2::geom_crossbar(
      data = means,
      ggplot2::aes(x = context, y = value, ymin = value, ymax = value),
      inherit.aes = FALSE, width = 0.28, linewidth = 0.72,
      colour = "black", show.legend = FALSE
    ) +
    ggplot2::annotate(
      "text", x = 1.5, y = y_annot,
      label = format_q_expr(endpoint_q[[value_column]]),
      parse = TRUE, size = 2.9, family = "Helvetica"
    ) +
    ggplot2::scale_fill_manual(values = c("Programmed cycle" = PC_ORANGE, "Natural cycle" = NC_BLUE)) +
    ggplot2::scale_shape_manual(values = c("Programmed cycle" = 21, "Natural cycle" = 22)) +
    ggplot2::scale_x_discrete(labels = c(
      "Programmed cycle" = "Programmed\n(n = 7)",
      "Natural cycle" = "Natural\n(n = 9)"
    )) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.07, 0.21))) +
    ggplot2::labs(x = NULL, y = y_label) +
    theme_submission() +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 7, 5, 7))
}

f1a <- small_n_panel("PGR_logCPM", "PGR TMM logCPM", 202609211L)
f1b <- small_n_panel("frozen_candidate25_score", "Frozen 22-gene WOI score", 202609212L)
f1c <- small_n_panel("ESR1_logCPM", "ESR1 TMM logCPM", 202609213L)
f1d <- small_n_panel("stable_epithelial_core_score", "Stable 9-gene epithelial-core score", 202609214L)
figure_1 <- compose_panels(f1a, f1b, f1c, f1d, labels = c("A", "B", "C", "D"), ncol = 2)
save_preview_bundle(figure_1, "Figure_1", width_mm = 183, height_mm = 132)

# Figure 2: one large gene-effect cloud plus a compact interval summary. The
# second panel visualizes the two descriptive quantities already reported in
# the locked concordance table; it does not add a new model or threshold.
effects <- effects[is.finite(effects$pc_nc) & is.finite(effects$as_hrt), ]
stopifnot(nrow(effects) == 10465L)
highlight <- effects[effects$gene %in% c("PGR", "ESR1"), ]
stopifnot(nrow(highlight) == 2L)

axis_limit <- ceiling(max(abs(c(effects$pc_nc, effects$as_hrt))))
highlight$label_x <- highlight$pc_nc + 0.33
highlight$label_y <- highlight$as_hrt + ifelse(highlight$gene == "PGR", 0.25, -0.25)

f2a <- ggplot2::ggplot(effects, ggplot2::aes(x = pc_nc, y = as_hrt)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.38, colour = NEUTRAL_MID) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.38, colour = NEUTRAL_MID) +
  ggplot2::geom_point(colour = NEUTRAL_MID, alpha = 0.20, size = 0.48) +
  ggplot2::geom_point(
    data = highlight,
    ggplot2::aes(fill = gene),
    shape = 21, colour = NEUTRAL_DARK, stroke = 0.45, size = 3.0,
    inherit.aes = TRUE
  ) +
  ggplot2::geom_label(
    data = highlight,
    ggplot2::aes(x = label_x, y = label_y, label = gene, colour = gene),
    hjust = 0, size = 3.2, fontface = "bold", family = "Helvetica",
    fill = "white", linewidth = 0.15, label.padding = grid::unit(0.10, "lines"),
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(values = c("PGR" = PGR_ORANGE, "ESR1" = ESR1_BLUE)) +
  ggplot2::scale_colour_manual(values = c("PGR" = PGR_ORANGE, "ESR1" = ESR1_BLUE)) +
  ggplot2::coord_equal(
    xlim = c(-axis_limit, axis_limit), ylim = c(-axis_limit, axis_limit),
    expand = FALSE, clip = "off"
  ) +
  ggplot2::labs(
    x = "PC minus NC epithelial effect\n(logCPM difference)",
    y = "AS/HRT minus local-control epithelial effect\n(logCPM difference)"
  ) +
  theme_submission() +
  ggplot2::theme(plot.margin = ggplot2::margin(10, 6, 6, 8))

summary_df <- data.frame(
  metric = factor(c("Spearman correlation", "Sign concordance"),
                  levels = c("Sign concordance", "Spearman correlation")),
  estimate = c(concordance$spearman, concordance$sign_concordance),
  low = c(concordance$spearman_bootstrap_ci95_low, concordance$sign_concordance_bootstrap_ci95_low),
  high = c(concordance$spearman_bootstrap_ci95_high, concordance$sign_concordance_bootstrap_ci95_high),
  colour = c(WOI_PURPLE, CORE_GREEN)
)
summary_df$value_label <- sprintf("%.3f", summary_df$estimate)
# Reserve a dedicated right-hand value column beyond the longest interval.
summary_df$label_x <- 0.955

f2b <- ggplot2::ggplot(summary_df, ggplot2::aes(y = metric)) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.38, colour = NEUTRAL_MID) +
  ggplot2::geom_segment(
    ggplot2::aes(x = low, xend = high, yend = metric),
    linewidth = 0.85, colour = NEUTRAL_DARK
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = estimate, fill = metric),
    shape = 21, colour = NEUTRAL_DARK, stroke = 0.45, size = 3.0
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = label_x, label = value_label),
    hjust = 1, size = 3.0, family = "Helvetica"
  ) +
  ggplot2::scale_fill_manual(values = c(
    "Spearman correlation" = WOI_PURPLE,
    "Sign concordance" = CORE_GREEN
  )) +
  ggplot2::scale_x_continuous(
    limits = c(0, 0.98), breaks = c(0, 0.2, 0.4, 0.6, 0.8),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::coord_cartesian(ylim = c(0.55, 2.45), clip = "off") +
  ggplot2::labs(x = "Estimate (bootstrap 95% CI)", y = NULL) +
  theme_submission() +
  ggplot2::theme(
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    # Keep the long x-axis title inside the final PDF page boundary.
    plot.margin = ggplot2::margin(10, 22, 6, 9)
  )

figure_2 <- compose_panels(
  f2a, f2b, labels = c("A", "B"), ncol = 2,
  rel_widths = c(1.55, 1.0)
)
save_preview_bundle(figure_2, "Figure_2", width_mm = 183, height_mm = 112)

# Figure 3: scale-specific forest panels with the classic visual language from
# zero-line and coloured-point language. Intervals and q values are unchanged.
endpoint_fill <- c(
  PGR_log2_FPKM = PGR_ORANGE,
  ESR1_log2_FPKM = ESR1_BLUE,
  frozen_WOI25_zscore = WOI_PURPLE,
  stable_epithelial_core_zscore = CORE_GREEN
)

make_forest_panel <- function(d, order_top, x_label) {
  d$label <- factor(d$label, levels = rev(order_top))
  d$q_label <- vapply(d$exact_BH_q_within_external_family, format_q_expr, character(1))
  x_values <- c(d$bootstrap_ci95_low, d$bootstrap_ci95_high, 0)
  span <- diff(range(x_values, na.rm = TRUE))
  if (!is.finite(span) || span <= 0) span <- 1
  q_x <- max(x_values, na.rm = TRUE) + 0.14 * span
  upper <- q_x + 0.43 * span
  lower <- min(x_values, na.rm = TRUE) - 0.10 * span

  ggplot2::ggplot(d, ggplot2::aes(y = label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.38, colour = NEUTRAL_MID) +
    ggplot2::geom_segment(
      ggplot2::aes(x = bootstrap_ci95_low, xend = bootstrap_ci95_high, yend = label),
      linewidth = 0.85, colour = NEUTRAL_DARK
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = difference_severe_IUA_minus_control, fill = endpoint),
      shape = 21, colour = NEUTRAL_DARK, stroke = 0.45, size = 3.0
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = q_x, label = q_label),
      hjust = 0, parse = TRUE, size = 2.7, family = "Helvetica", colour = NEUTRAL_DARK
    ) +
    ggplot2::scale_fill_manual(values = endpoint_fill, guide = "none") +
    ggplot2::coord_cartesian(xlim = c(lower, upper), clip = "off") +
    ggplot2::labs(x = x_label, y = NULL) +
    theme_submission() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 6, 8)
    )
}

f3a <- make_forest_panel(
  bulk[bulk$endpoint %in% c("PGR_log2_FPKM", "ESR1_log2_FPKM"), ],
  order_top = c("PGR", "ESR1"),
  x_label = "Severe IUA minus control\n(log2[FPKM + 1])"
)
f3b <- make_forest_panel(
  bulk[bulk$endpoint %in% c("frozen_WOI25_zscore", "stable_epithelial_core_zscore"), ],
  order_top = c("Frozen WOI score", "Epithelial-core score"),
  x_label = "Severe IUA minus control\n(score units)"
)

figure_3 <- compose_panels(f3a, f3b, labels = c("A", "B"), ncol = 2)
save_preview_bundle(figure_3, "Figure_3", width_mm = 183, height_mm = 84)

message("Rendered style-adapted main figures from locked source data: ", output_root)

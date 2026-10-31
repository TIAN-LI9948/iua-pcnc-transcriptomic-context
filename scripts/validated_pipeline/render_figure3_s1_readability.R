#!/usr/bin/env Rscript

# Readability repair for Figure 3 and Supplementary Figure S1.
# The script consumes only locked source-data tables and does not refit models.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggtext)
  library(cowplot)
  library(ragg)
  library(svglite)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript render_figure3_s1_readability.R SOURCE_DATA_DIR DESIGN_DIR OUTPUT_DIR", call. = FALSE)
}

# Keep caller-supplied ASCII junctions unresolved. Resolving them on Windows can
# reintroduce non-ASCII path components that some R file APIs cannot convert.
source_data_dir <- gsub("\\\\", "/", args[[1]])
design_dir <- gsub("\\\\", "/", args[[2]])
output_dir <- gsub("\\\\", "/", args[[3]])
if (!dir.exists(source_data_dir)) stop("Source-data directory does not exist: ", source_data_dir, call. = FALSE)
if (!dir.exists(design_dir)) stop("Design directory does not exist: ", design_dir, call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
layout_file <- normalizePath("C:/Users/zh/.codex/skills/zhplot/scripts/zhplot_layout.R", winslash = "/", mustWork = TRUE)
source(layout_file)

read_locked <- function(name) {
  path <- file.path(source_data_dir, name)
  if (!file.exists(path)) stop("Locked source table is missing: ", path, call. = FALSE)
  fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
}

f3 <- read_locked("Figure_3_endpoint_estimates.tsv")
f1 <- read_locked("Figure_1_endpoint_estimates.tsv")
s1_threshold <- read_locked("Supplementary_Figure_S1_threshold_eligibility.tsv")
s1_compartment <- read_locked("Supplementary_Figure_S1_compartment_sensitivity.tsv")
s1_lodo <- read_locked("Supplementary_Figure_S1_leave_one_donor_out.tsv")

stopifnot(
  nrow(f3) == 4L,
  setequal(f3$endpoint, c("PGR_log2_FPKM", "ESR1_log2_FPKM", "frozen_WOI25_zscore", "stable_epithelial_core_zscore")),
  nrow(f1) == 4L,
  setequal(s1_threshold$threshold, c(20, 50, 100)),
  all(s1_threshold$Total == 16L),
  all(s1_threshold$`Same donor IDs` == "Yes"),
  nrow(s1_compartment) == 6L,
  nrow(s1_lodo) == 64L,
  uniqueN(s1_lodo$left_out_donor) == 16L,
  uniqueN(s1_lodo$endpoint) == 4L,
  isTRUE(all.equal(f3[endpoint == "PGR_log2_FPKM", difference_severe_IUA_minus_control], -0.298876377634213, tolerance = 1e-12)),
  isTRUE(all.equal(f1[endpoint == "PGR_logCPM", difference_PC_minus_NC], 1.20080706566026, tolerance = 1e-12))
)

fmt_num <- function(x, digits = 2L) formatC(x, format = "f", digits = digits)
fmt_prob <- function(x) {
  vapply(x, function(z) {
    if (!is.finite(z)) return("NA")
    if (z < 0.001) return("<0.001")
    formatC(z, format = "f", digits = 3)
  }, character(1))
}
fmt_stat_html <- function(symbol, x) {
  prob <- fmt_prob(x)
  ifelse(
    startsWith(prob, "<"),
    paste0("<i>", symbol, "</i> &lt; ", substring(prob, 2L)),
    paste0("<i>", symbol, "</i> = ", prob)
  )
}

theme_pub <- function(base_size = 8.8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_blank(), plot.subtitle = element_blank(),
      axis.title = element_text(size = 9.0, colour = "#202428"),
      axis.text = element_text(size = 8.0, colour = "#202428"),
      axis.line = element_line(linewidth = 0.34, colour = "#34383C"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#34383C"),
      strip.background = element_blank(),
      strip.text = element_text(size = 8.1, face = "bold", colour = "#202428"),
      panel.grid = element_blank(), legend.position = "none",
      plot.margin = margin(4, 5, 4, 5, unit = "pt")
    )
}

manifest_for <- function(id, roles, colours, brief, concept, first, second, third, use_forest_reference = FALSE) {
  list(
    figureId = id,
    journalTarget = "International Journal of Women's Health candidate",
    panelRoles = roles,
    semanticColors = as.list(colours),
    sharedLegends = list(),
    fontSizes = list(panelLabelPt = 10, axisTitlePt = 8.5, axisTickPt = 7.5, legendPt = 7.3, annotationMinPt = 7.5),
    gridColumns = 12L,
    outerMarginMm = 3.0,
    designProvenance = list(
      plotcase_status = if (use_forest_reference) "USED" else "NO_COMPATIBLE_RESULT",
      plotcase_queries = if (use_forest_reference) list(
        "effect-size forest plot grouped comparison confidence interval annotation",
        "forest plot effect estimate confidence interval aligned numeric columns q value"
      ) else list("eligibility matrix and influence range plot"),
      plotcase_candidates = if (use_forest_reference) list("PlotCase effect-size forest plot") else list(),
      plotcase_reuse_level = if (use_forest_reference) "STYLE_ONLY" else "NOT_USED",
      plotcase_reason = if (use_forest_reference) {
        "A compatible forest case informed interval hierarchy and aligned numerical columns; no statistical code or data were reused."
      } else {
        "No compatible mixed eligibility-matrix and donor-influence case was found; the locked analysis contract determined the design."
      },
      plotcase_candidate_sufficiency = if (use_forest_reference) "LIMITED" else "NONE",
      plotcase_compatible_candidate_count = if (use_forest_reference) 1L else 0L,
      user_reference_status = if (use_forest_reference) "FOUND" else "NOT_FOUND",
      user_reference_ids = if (use_forest_reference) "user-subgroup-hr-ci-forest-001" else character(),
      user_reference_id = if (use_forest_reference) "user-subgroup-hr-ci-forest-001" else NULL,
      user_reference_transferability = if (use_forest_reference) "VISUAL_REFERENCE_ONLY" else NULL,
      reference_features_used = if (use_forest_reference) "lightweight zero line, interval hierarchy, aligned numerical columns" else character(),
      reference_adaptation_notes = if (use_forest_reference) "Presentation geometry only; the locked estimands, intervals, and q values were retained." else NULL,
      scidraw_status = "NOT_REQUIRED_WITH_REASON",
      scidraw_reason = "All panels are quantitative data displays.",
      design_brief = file.path(design_dir, brief),
      selected_concept = concept,
      visual_intent = list(first_read = first, second_read = second, third_read = third),
      panel_visual_review = "MANUAL_PANEL_VISUAL_REVIEW_REQUIRED"
    )
  )
}

make_table_forest <- function(d, endpoint_order, x_label) {
  d <- copy(d[endpoint %in% endpoint_order])
  d[, order__ := match(endpoint, endpoint_order)]
  setorder(d, order__)
  d[, order__ := NULL]
  d[, y := rev(seq_len(.N))]
  d[, endpoint_label := fifelse(
    endpoint == "PGR_log2_FPKM", "PGR",
    fifelse(endpoint == "ESR1_log2_FPKM", "ESR1",
      fifelse(endpoint == "frozen_WOI25_zscore", "Frozen WOI\nprogramme", "Epithelial-core\nprogramme"))
  )]
  d[, value_label := paste0(
    fmt_num(difference_severe_IUA_minus_control), " [",
    fmt_num(bootstrap_ci95_low), ", ", fmt_num(bootstrap_ci95_high), "]"
  )]
  d[, stat_label := paste0(
    fmt_stat_html("P", exact_permutation_p_two_sided),
    "; ", fmt_stat_html("q", exact_BH_q_within_external_family)
  )]

  xr <- range(c(0, d$bootstrap_ci95_low, d$bootstrap_ci95_high), finite = TRUE)
  pad <- max(diff(xr) * 0.10, 0.05)
  forest_lo <- xr[[1]] - pad
  forest_hi <- xr[[2]] + pad
  span <- forest_hi - forest_lo
  value_x <- forest_hi + 0.20 * span
  stat_x <- forest_hi + 1.78 * span
  full_hi <- forest_hi + 2.95 * span
  axis_y <- 0.48
  axis_breaks <- pretty(c(forest_lo, forest_hi), n = 3)
  axis_breaks <- axis_breaks[axis_breaks >= forest_lo & axis_breaks <= forest_hi]

  ggplot(d) +
    annotate("rect", xmin = forest_lo, xmax = full_hi, ymin = 2.47, ymax = 2.92, fill = "#EDF2F5", colour = NA) +
    annotate("rect", xmin = forest_lo, xmax = full_hi, ymin = 1.60, ymax = 2.40, fill = "#F7F9FA", colour = NA) +
    annotate("rect", xmin = forest_lo, xmax = full_hi, ymin = 0.60, ymax = 1.40, fill = "#FFFFFF", colour = NA) +
    annotate("segment", x = 0, xend = 0, y = 0.57, yend = 2.42, linetype = "22", linewidth = 0.40, colour = "#747B82") +
    geom_segment(aes(x = bootstrap_ci95_low, xend = bootstrap_ci95_high, y = y, yend = y),
      linewidth = 0.98, lineend = "round", colour = "#858B91") +
    geom_point(aes(x = difference_severe_IUA_minus_control, y = y), shape = 21, size = 3.0,
      fill = "#3F6E8C", colour = "#252A2E", stroke = 0.42) +
    geom_text(aes(x = value_x, y = y, label = value_label), hjust = 0, size = 2.67,
      family = "sans", colour = "#25292D") +
    ggtext::geom_richtext(aes(x = stat_x, y = y, label = stat_label), hjust = 0, size = 2.58,
      family = "sans", colour = "#25292D", fill = NA, label.color = NA,
      label.padding = grid::unit(0, "pt")) +
    annotate("text", x = mean(c(forest_lo, forest_hi)), y = 2.69, label = "Effect",
      size = 2.82, family = "sans", fontface = "bold", colour = "#202428") +
    annotate("text", x = value_x, y = 2.69, label = "Difference [95% CI]", hjust = 0,
      size = 2.48, family = "sans", fontface = "bold", colour = "#202428") +
    annotate("richtext", x = stat_x, y = 2.69, label = "<i>P</i>; <i>q</i>", hjust = 0,
      size = 2.55, family = "sans", fontface = "bold", colour = "#202428",
      fill = NA, label.color = NA, label.padding = grid::unit(0, "pt")) +
    annotate("segment", x = forest_lo, xend = forest_hi, y = axis_y, yend = axis_y,
      linewidth = 0.40, colour = "#34383C") +
    geom_segment(data = data.table(x = axis_breaks), aes(x = x, xend = x, y = axis_y, yend = axis_y - 0.07),
      inherit.aes = FALSE, linewidth = 0.34, colour = "#34383C") +
    geom_text(data = data.table(x = axis_breaks), aes(x = x, y = axis_y - 0.14, label = fmt_num(x, 1)),
      inherit.aes = FALSE, vjust = 1, size = 2.52, family = "sans", colour = "#34383C") +
    annotate("text", x = mean(c(forest_lo, forest_hi)), y = 0.02, label = x_label,
      vjust = 1, size = 2.72, family = "sans", colour = "#24282C") +
    scale_y_continuous(breaks = d$y, labels = d$endpoint_label, limits = c(-0.32, 3.00), expand = expansion(mult = 0)) +
    scale_x_continuous(limits = c(forest_lo, full_hi), expand = expansion(mult = 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = "sans", base_size = 8.8) +
    theme(
      axis.text.y = element_text(size = 8.3, colour = "#24282C", hjust = 1,
        lineheight = 0.95, margin = margin(r = 7, unit = "pt")),
      plot.margin = margin(3, 7, 4, 7, unit = "pt")
    )
}

# Figure 3: shorten x-axis labels while preserving complete numerical columns.
p3a <- make_table_forest(f3, c("PGR_log2_FPKM", "ESR1_log2_FPKM"), "IUA - control (log2 scale)")
p3b <- make_table_forest(f3, c("frozen_WOI25_zscore", "stable_epithelial_core_zscore"), "IUA - control (z score)")

f3_manifest <- manifest_for(
  "Figure_3", c("primary", "primary"),
  c(point = "#3F6E8C", interval = "#858B91", null = "#747B82"),
  "Figure_3.yaml", "compact_table_forest_two_panel",
  "The untimed bulk result is not a clean replication of the primary epithelial contrast.",
  "Receptor and programme estimates retain distinct axes and units.",
  "Aligned columns provide each difference, 95% confidence interval, exact permutation P value, and four-endpoint BH q value.",
  use_forest_reference = TRUE
)

fig3 <- zhplot_compose_regular(
  list(p3a, p3b), ncol = 2, nrow = 1,
  panel_ids = c("receptor_transcripts", "programme_scores"),
  alignment = "independent", rel_widths = c(0.93, 1.07),
  tag_size = 10, tag_band_fraction = 0.095, font_family = "sans",
  gutter_x = 0.035, title_policy = "strip", figure_manifest = f3_manifest
)

f3_bundle <- zhplot_render_bundle(
  fig3, file.path(output_dir, "Figure_3"), width_mm = 183, height_mm = 68, dpi = 600,
  source_file = script_path,
  input_files = file.path(source_data_dir, "Figure_3_endpoint_estimates.tsv"),
  implementation_files = layout_file,
  include_tiff = FALSE
)

# Supplementary Figure S1A: threshold eligibility is an invariance matrix, not a trajectory.
threshold <- copy(s1_threshold)
ge_symbol <- intToUtf8(8805)
threshold[, threshold_label := factor(
  paste0(ge_symbol, threshold), levels = paste0(ge_symbol, c(100, 50, 20))
)]
threshold_long <- melt(
  threshold[, .(threshold_label, `PC donors` = as.character(PC), `NC donors` = as.character(NC),
    `All donors` = as.character(Total), `Identical set` = as.character(`Same donor IDs`))],
  id.vars = "threshold_label", variable.name = "column", value.name = "value"
)
threshold_long[, column := factor(column, levels = c("PC donors", "NC donors", "All donors", "Identical set"))]
threshold_long[, fill_group := fifelse(column == "Identical set", "identity", "count")]

pS1a <- ggplot(threshold_long, aes(column, threshold_label)) +
  geom_tile(aes(fill = fill_group), colour = "white", linewidth = 0.9, show.legend = FALSE) +
  geom_text(data = threshold_long[fill_group == "count"], aes(label = value),
    size = 3.0, family = "sans", fontface = "plain", colour = "#25313A") +
  geom_text(data = threshold_long[fill_group == "identity"], aes(label = value),
    size = 3.0, family = "sans", fontface = "bold", colour = "#25313A") +
  scale_fill_manual(values = c(count = "#EDF2F5", identity = "#D8EAF1")) +
  scale_x_discrete(position = "top", labels = c("PC\ndonors", "NC\ndonors", "All\ndonors", "Same donor\nset")) +
  labs(x = NULL, y = "Minimum nuclei per donor") +
  theme_void(base_family = "sans", base_size = 8.8) +
  theme(
    axis.text.x = element_text(size = 7.7, face = "bold", colour = "#25313A", lineheight = 0.92, margin = margin(b = 4)),
    axis.text.y = element_text(size = 8.2, colour = "#25313A", margin = margin(r = 5)),
    axis.title.y = element_text(size = 8.8, colour = "#25313A", angle = 90, margin = margin(r = 6)),
    plot.margin = margin(6, 5, 6, 5, unit = "pt")
  )

primary_comp <- f1[endpoint %in% c("PGR_logCPM", "frozen_candidate25_score"), .(
  endpoint, difference_PC_minus_NC, bootstrap_ci95_low, bootstrap_ci95_high,
  exact_BH_q_within_family, compartment = "primary_nonciliated"
)]
comp <- rbindlist(list(primary_comp, s1_compartment[, .(
  endpoint, difference_PC_minus_NC, bootstrap_ci95_low, bootstrap_ci95_high,
  exact_BH_q_within_family, compartment
)]), use.names = TRUE)

make_compartment_panel <- function(compartment_id, colour, right_margin_pt) {
  d <- copy(comp[compartment == compartment_id & endpoint %in% c("PGR_logCPM", "frozen_candidate25_score")])
  d[, endpoint_label := factor(endpoint,
    levels = c("PGR_logCPM", "frozen_candidate25_score"),
    labels = c("PGR", "Frozen WOI programme"))]
  d[, value_label := paste0(
    fmt_num(difference_PC_minus_NC), " [", fmt_num(bootstrap_ci95_low), ", ",
    fmt_num(bootstrap_ci95_high), "]<br>", fmt_stat_html("q", exact_BH_q_within_family)
  )]

  ggplot(d, aes(difference_PC_minus_NC, 1)) +
    geom_vline(xintercept = 0, linetype = "22", linewidth = 0.35, colour = "#7F858A") +
    geom_segment(aes(x = bootstrap_ci95_low, xend = bootstrap_ci95_high, yend = 1),
      linewidth = 0.86, lineend = "round", colour = colour) +
    geom_point(shape = 21, size = 2.75, fill = colour, colour = "white", stroke = 0.45) +
    geom_point(shape = 21, size = 2.75, fill = NA, colour = "#30363B", stroke = 0.35) +
    ggtext::geom_richtext(aes(x = Inf, label = value_label), hjust = -0.04,
      size = 2.42, family = "sans", colour = "#25292D", fill = NA,
      label.color = NA, label.padding = grid::unit(0, "pt"), lineheight = 1.30) +
    facet_wrap(~endpoint_label, ncol = 1, scales = "free_x", axes = "all_x", axis.labels = "all") +
    scale_y_continuous(limits = c(0.82, 1.18), breaks = NULL, expand = expansion(mult = 0)) +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.08))) +
    labs(x = "PC - NC difference", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_pub() +
    theme(
      strip.text = element_text(size = 7.9, face = "bold", hjust = 0),
      strip.placement = "outside",
      panel.spacing.y = grid::unit(2.8, "mm"),
      axis.title.x = element_text(size = 8.5, margin = margin(t = 4)),
      axis.text.x = element_text(size = 7.5),
      plot.margin = margin(3, right_margin_pt, 4, 5, unit = "pt")
    )
}

pS1b <- make_compartment_panel("glandular_epithelial", "#2C7FB8", right_margin_pt = 62)
pS1c <- make_compartment_panel("broad_epithelial_audit", "#E67E22", right_margin_pt = 68)

# Supplementary Figure S1D: retain all donor omissions and add a forest-style summary column.
lodo <- merge(
  copy(s1_lodo),
  f1[, .(endpoint, full_effect = difference_PC_minus_NC)],
  by = "endpoint", all.x = TRUE, sort = FALSE
)
stopifnot(all(is.finite(lodo$full_effect)), all(lodo$full_effect != 0))
lodo[, relative_effect := difference_PC_minus_NC / full_effect]
endpoint_order <- c("PGR_logCPM", "frozen_candidate25_score", "ESR1_logCPM", "stable_epithelial_core_score")
endpoint_labels <- c(
  PGR_logCPM = "PGR",
  frozen_candidate25_score = "Frozen WOI programme",
  ESR1_logCPM = "ESR1",
  stable_epithelial_core_score = "Epithelial-core programme"
)
lodo[, y := rev(seq_along(endpoint_order))[match(endpoint, endpoint_order)]]
lodo_summary <- lodo[, .(
  ratio_median = median(relative_effect),
  ratio_min = min(relative_effect),
  ratio_max = max(relative_effect),
  q_min = min(exact_BH_q_within_family),
  q_max = max(exact_BH_q_within_family)
), by = endpoint]
lodo_summary[, y := rev(seq_along(endpoint_order))[match(endpoint, endpoint_order)]]
lodo_summary[, endpoint_label := unname(endpoint_labels[endpoint])]
lodo_summary[, ratio_label := paste0(fmt_num(ratio_median), " [", fmt_num(ratio_min), ", ", fmt_num(ratio_max), "]")]
lodo_summary[, q_label := paste0(fmt_prob(q_min), "-", fmt_prob(q_max))]

forest_lo <- min(0.84, min(lodo_summary$ratio_min) - 0.01)
forest_hi <- max(1.16, max(lodo_summary$ratio_max) + 0.01)
span <- forest_hi - forest_lo
ratio_x <- forest_hi + 0.15 * span
full_hi <- forest_hi + 1.00 * span
axis_y <- 0.38
axis_breaks <- c(0.9, 1.0, 1.1)
endpoint_cols <- c(
  PGR_logCPM = "#3F6E8C", frozen_candidate25_score = "#2C7FB8",
  ESR1_logCPM = "#7A6FB5", stable_epithelial_core_score = "#168B73"
)

pS1d <- ggplot() +
  annotate("rect", xmin = forest_lo, xmax = full_hi, ymin = 4.48, ymax = 4.90, fill = "#EDF2F5", colour = NA) +
  annotate("segment", x = 1, xend = 1, y = 0.56, yend = 4.42, linetype = "22", linewidth = 0.40, colour = "#B2182B") +
  geom_segment(data = lodo_summary, aes(x = ratio_min, xend = ratio_max, y = y, yend = y, colour = endpoint),
    linewidth = 0.92, lineend = "round", show.legend = FALSE) +
  geom_point(data = lodo, aes(relative_effect, y, colour = endpoint),
    position = position_jitter(width = 0, height = 0.105, seed = 20260922),
    size = 1.25, alpha = 0.55, show.legend = FALSE) +
  geom_point(data = lodo_summary, aes(ratio_median, y, fill = endpoint),
    shape = 21, size = 2.65, colour = "#30363B", stroke = 0.38, show.legend = FALSE) +
  geom_text(data = lodo_summary, aes(x = ratio_x, y = y, label = ratio_label),
    hjust = 0, size = 2.48, family = "sans", colour = "#25292D") +
  annotate("text", x = mean(c(forest_lo, forest_hi)), y = 4.69, label = "Influence range",
    size = 2.68, family = "sans", fontface = "bold", colour = "#202428") +
  annotate("text", x = ratio_x, y = 4.69, label = "Median [range]", hjust = 0,
    size = 2.48, family = "sans", fontface = "bold", colour = "#202428") +
  annotate("segment", x = forest_lo, xend = forest_hi, y = axis_y, yend = axis_y,
    linewidth = 0.40, colour = "#34383C") +
  geom_segment(data = data.table(x = axis_breaks), aes(x = x, xend = x, y = axis_y, yend = axis_y - 0.07),
    inherit.aes = FALSE, linewidth = 0.34, colour = "#34383C") +
  geom_text(data = data.table(x = axis_breaks), aes(x = x, y = axis_y - 0.14, label = fmt_num(x, 1)),
    inherit.aes = FALSE, vjust = 1, size = 2.50, family = "sans", colour = "#34383C") +
  annotate("text", x = mean(c(forest_lo, forest_hi)), y = -0.13, label = "Relative effect (LODO/full)",
    vjust = 1, size = 2.64, family = "sans", colour = "#24282C") +
  scale_colour_manual(values = endpoint_cols) +
  scale_fill_manual(values = endpoint_cols) +
  scale_y_continuous(
    breaks = lodo_summary$y,
    labels = lodo_summary$endpoint_label,
    limits = c(-0.42, 4.92), expand = expansion(mult = 0)
  ) +
  scale_x_continuous(limits = c(forest_lo, full_hi), expand = expansion(mult = 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = "sans", base_size = 8.8) +
  theme(
    axis.text.y = element_text(size = 7.7, colour = "#24282C", hjust = 1, margin = margin(r = 6)),
    plot.margin = margin(3, 6, 4, 0, unit = "pt")
  )

s1_manifest <- manifest_for(
  "Supplementary_Figure_S1", c("supporting", "primary", "primary", "supporting"),
  c(matrix = "#D8EAF1", glandular = "#2C7FB8", broad = "#E67E22", reference = "#B2182B"),
  "Supplementary_Figure_S1.yaml", "eligibility_matrix_numeric_forests_influence_range",
  "The same 16 donors are retained at every eligibility threshold.",
  "Glandular and broad-epithelial sensitivity estimates are shown with exact intervals and q values.",
  "All donor omissions remain directionally stable, with every recalculation and the median range exposed.",
  use_forest_reference = TRUE
)

figS1 <- zhplot_compose_regular(
  list(pS1a, pS1b, pS1c, pS1d),
  ncol = 2, nrow = 2,
  panel_ids = c("threshold_eligibility", "glandular_sensitivity", "broad_epithelial_sensitivity", "leave_one_donor_out"),
  alignment = "independent", rel_widths = c(1, 1), rel_heights = c(1, 1),
  tag_size = 10, tag_band_fraction = 0.080, font_family = "sans",
  gutter_x = 0.035, gutter_y = 0.040, title_policy = "strip", figure_manifest = s1_manifest
)

s1_bundle <- zhplot_render_bundle(
  figS1, file.path(output_dir, "Supplementary_Figure_S1"),
  width_mm = 183, height_mm = 136, dpi = 600,
  source_file = script_path,
  input_files = c(
    file.path(source_data_dir, "Figure_1_endpoint_estimates.tsv"),
    file.path(source_data_dir, "Supplementary_Figure_S1_threshold_eligibility.tsv"),
    file.path(source_data_dir, "Supplementary_Figure_S1_compartment_sensitivity.tsv"),
    file.path(source_data_dir, "Supplementary_Figure_S1_leave_one_donor_out.tsv")
  ),
  implementation_files = layout_file,
  include_tiff = FALSE
)

qa <- list(
  schema_version = "1.0",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  scientific_inputs = list(
    figure3_rows = nrow(f3),
    threshold_rows = nrow(s1_threshold),
    threshold_donor_totals = as.list(s1_threshold$Total),
    identical_donor_sets = all(s1_threshold$`Same donor IDs` == "Yes"),
    compartment_rows = nrow(s1_compartment),
    lodo_rows = nrow(s1_lodo),
    lodo_donors = uniqueN(s1_lodo$left_out_donor),
    lodo_endpoints = uniqueN(s1_lodo$endpoint)
  ),
  rendering = list(
    backend = "R",
    R_version = R.version.string,
    ggplot2 = as.character(packageVersion("ggplot2")),
    data_table = as.character(packageVersion("data.table")),
    figure3 = f3_bundle,
    supplementary_figure_s1 = s1_bundle
  )
)
write_json(qa, file.path(output_dir, "figure3_s1_render_qa.json"), pretty = TRUE, auto_unbox = TRUE, null = "null", digits = 12)

message("Rendered Figure 3 and Supplementary Figure S1 from locked source tables.")

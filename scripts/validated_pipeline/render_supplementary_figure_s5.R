#!/usr/bin/env Rscript

# Visual repair of Supplementary Figure S5 from the locked 24-row table.
# This script preserves all patient values, endpoint-specific scales, pairing,
# and the descriptive-only contract. It does not fit a model or calculate P values.

options(stringsAsFactors = FALSE, warn = 1)

needed <- c("cowplot", "data.table", "digest", "ggplot2", "jsonlite", "ragg", "svglite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript render_supplementary_figure_s5.R VALUES_TSV ROUNDTRIP_TSV DESIGN_YAML OUTPUT_BASE", call. = FALSE)
}
# Keep caller-supplied ASCII junctions unresolved on Windows.
values_path <- gsub("\\\\", "/", args[[1]])
roundtrip_path <- gsub("\\\\", "/", args[[2]])
design_path <- gsub("\\\\", "/", args[[3]])
if (!file.exists(values_path)) stop("Values table is missing: ", values_path, call. = FALSE)
if (!file.exists(roundtrip_path)) stop("Round-trip QA table is missing: ", roundtrip_path, call. = FALSE)
if (!file.exists(design_path)) stop("Design contract is missing: ", design_path, call. = FALSE)
output_base <- gsub("\\\\", "/", args[[4]])
dir.create(dirname(output_base), recursive = TRUE, showWarnings = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
layout_file <- normalizePath(
  "C:/Users/zh/.codex/skills/zhplot/scripts/zhplot_layout.R",
  winslash = "/", mustWork = TRUE
)
source(layout_file)

s5 <- fread(values_path, sep = "\t", header = TRUE, showProgress = FALSE)
qa <- fread(roundtrip_path, sep = "\t", header = TRUE, showProgress = FALSE)

required_cols <- c("Patient", "endpoint", "value", "stage", "patient_alias")
stopifnot(
  all(required_cols %in% names(s5)),
  nrow(s5) == 24L,
  uniqueN(s5$Patient) == 3L,
  uniqueN(s5$endpoint) == 4L,
  all(s5[, .N, by = endpoint]$N == 6L),
  all(s5[, uniqueN(stage), by = .(endpoint, Patient)]$V1 == 2L),
  setequal(s5$stage, c("pre", "post")),
  setequal(s5$patient_alias, c("P1", "P2", "P3"))
)
qa_status <- setNames(qa$status, qa$qa_check)
stopifnot(
  qa_status[["eligible_pair_count"]] == "PASS",
  qa_status[["formal_inference_allowed"]] == "PASS",
  qa_status[["endpoint_long_rows"]] == "PASS",
  qa_status[["endpoint_matrix_roundtrip_exact"]] == "PASS"
)

endpoint_order <- c(
  "PGR", "portable_WOI22_zscore",
  "ESR1", "portable_epithelial_core9_zscore"
)
axis_labels <- c(
  PGR = "PGR (logCPM)",
  portable_WOI22_zscore = "WOI programme score",
  ESR1 = "ESR1 (logCPM)",
  portable_epithelial_core9_zscore = "Epithelial-core score"
)
patient_colours <- c(P1 = "#2F6FA3", P2 = "#D97829", P3 = "#2A9D78")
patient_shapes <- c(P1 = 21, P2 = 22, P3 = 24)

s5[, patient_alias := factor(patient_alias, levels = names(patient_colours))]
s5[, stage_label := factor(stage, levels = c("pre", "post"), labels = c("Pre", "Post"))]
s5[, endpoint := factor(endpoint, levels = endpoint_order)]
setorder(s5, endpoint, patient_alias, stage_label)

theme_s5 <- function() {
  theme_classic(base_size = 8.2, base_family = "sans") +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8.6, colour = "#202326", margin = margin(r = 6)),
      axis.text = element_text(size = 7.6, colour = "#202326"),
      axis.text.x = element_text(size = 7.9, margin = margin(t = 3)),
      axis.line = element_line(linewidth = 0.34, colour = "#383B3E"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#383B3E"),
      axis.ticks.length = grid::unit(1.5, "mm"),
      panel.grid.major.y = element_line(linewidth = 0.24, colour = "#E7EBEE"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 7.4, face = "bold", colour = "#202326"),
      legend.text = element_text(size = 7.3, colour = "#202326", margin = margin(l = 2.5)),
      legend.key.width = grid::unit(6.2, "mm"),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.spacing.x = grid::unit(1.8, "mm"),
      legend.box.spacing = grid::unit(0, "pt"),
      legend.margin = margin(t = 1, b = 0, unit = "pt"),
      plot.margin = margin(3, 5, 4, 5, unit = "pt")
    )
}

make_panel <- function(ep, show_legend = FALSE) {
  d <- s5[endpoint == ep]
  stopifnot(nrow(d) == 6L)
  p <- ggplot(
    d,
    aes(x = stage_label, y = value, group = patient_alias,
        colour = patient_alias, shape = patient_alias)
  ) +
    geom_line(linewidth = 0.72, alpha = 0.90, lineend = "round") +
    geom_point(size = 2.75, stroke = 0.78, fill = "white") +
    scale_colour_manual(values = patient_colours, name = "Patient") +
    scale_shape_manual(values = patient_shapes, name = "Patient") +
    scale_x_discrete(expand = expansion(add = 0.42)) +
    scale_y_continuous(n.breaks = 4, expand = expansion(mult = c(0.10, 0.12))) +
    labs(x = NULL, y = axis_labels[[ep]], title = NULL, subtitle = NULL) +
    theme_s5()

  if (show_legend) {
    p +
      theme(legend.position = "bottom", legend.direction = "horizontal") +
      guides(
        colour = guide_legend(
          nrow = 1, byrow = TRUE,
          override.aes = list(linewidth = 0.65, size = 2.05, fill = "white")
        ),
        shape = "none"
      )
  } else {
    p + theme(legend.position = "none")
  }
}

plots <- lapply(seq_along(endpoint_order), function(i) {
  make_panel(endpoint_order[[i]], show_legend = i == 4L)
})

figure_manifest <- list(
  figureId = "Supplementary_Figure_S5",
  journalTarget = "International Journal of Women's Health candidate",
  panelRoles = rep("primary", 4),
  semanticColors = as.list(patient_colours),
  sharedLegends = list(list(
    panels = paste0("paired_", endpoint_order),
    mapping = "patient_alias",
    owner = "paired_portable_epithelial_core9_zscore"
  )),
  fontSizes = list(
    panelLabelPt = 10, axisTitlePt = 8.6, axisTickPt = 7.6,
    legendPt = 7.3, annotationMinPt = 7.6
  ),
  gridColumns = 12L,
  outerMarginMm = 3.0,
  designProvenance = list(
    plotcase_status = "USED",
    plotcase_queries = list(
      "paired longitudinal slope plot three patients two timepoints four endpoints direct labels",
      "small-n paired before after connected dot plot patient trajectories",
      "multi-group paired line plot three patients before after four endpoints"
    ),
    plotcase_candidates = list("PlotCase multi-group paired line plot"),
    plotcase_reuse_level = "ADAPT_DESIGN_ONLY",
    plotcase_reason = paste(
      "The compatible case informed thin-line hierarchy and whitespace only;",
      "its normalization, parallel-coordinate axes, faceting, and statistical logic were rejected."
    ),
    plotcase_candidate_sufficiency = "LIMITED",
    plotcase_compatible_candidate_count = 1L,
    user_reference_status = "NOT_FOUND",
    user_reference_ids = character(),
    reference_features_used = c("paired-line hierarchy", "compact marks", "restrained whitespace"),
    scidraw_status = "NOT_REQUIRED_WITH_REASON",
    scidraw_reason = "The figure is a quantitative paired display and needs no schematic asset.",
    design_brief = design_path,
    selected_concept = "paired_2x2",
    visual_intent = list(
      first_read = "Only three patient trajectories are shown.",
      second_read = "Directions differ by endpoint and sometimes by patient.",
      third_read = "Patient colours and shapes are fixed across panels without formal inference."
    ),
    panel_visual_review = "MANUAL_PANEL_VISUAL_REVIEW_REQUIRED"
  )
)

fig <- zhplot_compose_regular(
  plots,
  ncol = 2,
  nrow = 2,
  panel_ids = paste0("paired_", endpoint_order),
  alignment = "independent",
  rel_widths = c(1, 1),
  rel_heights = c(1, 1),
  tag_size = 10,
  tag_band_fraction = 0.076,
  font_family = "sans",
  gutter_x = 0.035,
  gutter_y = 0.034,
  side_pad_pt = 4,
  title_policy = "strip",
  figure_manifest = figure_manifest
)

# The layout helper owns all devices: svglite::svglite for SVG,
# grDevices::cairo_pdf for PDF, and ragg::agg_png at 600 dpi for PNG.
bundle <- zhplot_render_bundle(
  fig,
  output_base,
  width_mm = 183,
  height_mm = 128,
  dpi = 600,
  source_file = script_path,
  input_files = c(values_path, roundtrip_path, design_path),
  implementation_files = layout_file,
  include_tiff = FALSE
)

for (required in c("pdf", "png", "svg", "manifest", "figure_manifest")) {
  if (!file.exists(bundle[[required]])) stop("Missing render output: ", required, call. = FALSE)
}

cat(jsonlite::toJSON(list(
  status = "PASS",
  rows = nrow(s5),
  patients = uniqueN(s5$Patient),
  endpoints = uniqueN(s5$endpoint),
  source_sha256 = digest::digest(file = values_path, algo = "sha256"),
  outputs = bundle
), auto_unbox = TRUE, pretty = TRUE), "\n")

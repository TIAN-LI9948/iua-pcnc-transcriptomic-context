#!/usr/bin/env Rscript

# Re-render Supplementary Figure S4 from locked submission source tables.
# This script changes only the visual encoding. It does not recalculate the
# embedding, annotation classes, lineage scores, cell-count gate, or eligibility.
# Vector/raster devices are owned by zhplot_render_bundle: svglite::svglite,
# grDevices::cairo_pdf, and ragg::agg_png at 600 dpi.

options(stringsAsFactors = FALSE, warn = 1)

needed <- c(
  "data.table", "digest", "ggplot2", "jsonlite", "patchwork",
  "ragg", "scales", "scattermore", "svglite"
)
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript render_supplementary_figure_s4.R SOURCE_DIR OUTPUT_BASE", call. = FALSE)
}
# Preserve an ASCII junction supplied by the caller. Resolving it on Windows can
# reintroduce non-ASCII path components that R's file APIs cannot translate.
source_dir <- gsub("\\\\", "/", args[[1]])
if (!dir.exists(source_dir)) stop("Source directory does not exist: ", source_dir, call. = FALSE)
output_base <- gsub("\\\\", "/", args[[2]])
dir.create(dirname(output_base), recursive = TRUE, showWarnings = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
layout_file <- normalizePath(
  "C:/Users/zh/.codex/skills/zhplot/scripts/zhplot_layout.R",
  winslash = "/", mustWork = TRUE
)
source(layout_file)

paths <- list(
  umap = file.path(source_dir, "Supplementary_Figure_S4_umap_classes.tsv.gz"),
  heatmap = file.path(source_dir, "Supplementary_Figure_S4_lineage_heatmap.tsv"),
  eligibility = file.path(source_dir, "Supplementary_Figure_S4_eligibility_table.tsv")
)
if (!all(file.exists(unlist(paths)))) {
  stop("Missing locked S4 source table(s): ",
       paste(names(paths)[!file.exists(unlist(paths))], collapse = ", "), call. = FALSE)
}

read_tsv <- function(path) {
  fread(path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
}
s4_umap <- read_tsv(paths$umap)
heat_long <- read_tsv(paths$heatmap)
elig_wide <- read_tsv(paths$eligibility)

class_levels <- c(
  "Nonciliated epithelium", "Ciliated epithelium", "Excluded epithelium",
  "Stromal / perivascular", "Endothelial", "Immune", "Erythrocyte"
)
class_colours <- c(
  "Nonciliated epithelium" = "#D55E00",
  "Ciliated epithelium" = "#E69F00",
  "Excluded epithelium" = "#A7A9AC",
  "Stromal / perivascular" = "#009E73",
  "Endothelial" = "#0072B2",
  "Immune" = "#6F58A8",
  "Erythrocyte" = "#B64E82"
)
stage_colours <- c("Pre" = "#2C7FB8", "Post" = "#D97706")
status_fills <- c(
  "Eligible" = "#E4F0F7",
  "Below gate" = "#F5F6F7",
  "Post only" = "#FBE9E8"
)

# Scientific and row-count locks.
stopifnot(nrow(s4_umap) == 51018L)
stopifnot(setequal(unique(s4_umap$annotation_class), class_levels))
stopifnot(nrow(heat_long) == 49L)
stopifnot(setequal(unique(heat_long$annotation_class), class_levels))
stopifnot(nrow(elig_wide) == 9L)
stopifnot(sum(elig_wide$complete_prepost) == 8L)
stopifnot(sum(elig_wide$eligible_min20_each_stage) == 3L)
stopifnot(identical(sort(elig_wide[eligible_min20_each_stage == TRUE, Patient]), c("01-13", "01-14", "01-15")))
stopifnot(identical(elig_wide[complete_prepost == FALSE, Patient], "01-05"))

theme_pub <- function(base_size = 8.8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_blank(), plot.subtitle = element_blank(),
      axis.title = element_text(size = 9.0, colour = "#1F1F1F"),
      axis.text = element_text(size = 7.8, colour = "#1F1F1F"),
      axis.line = element_line(linewidth = 0.34, colour = "#333333"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#333333"),
      legend.title = element_text(size = 7.7, colour = "#222222"),
      legend.text = element_text(size = 7.5, colour = "#222222"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5, unit = "pt")
    )
}

# Panel A: locked UMAP with category-matched label fills.
s4_umap[, annotation_class := factor(annotation_class, levels = class_levels)]
draw_order <- c(
  "Excluded epithelium", "Stromal / perivascular", "Endothelial",
  "Immune", "Erythrocyte", "Ciliated epithelium", "Nonciliated epithelium"
)
s4_umap[, draw_rank := match(as.character(annotation_class), draw_order)]
setorder(s4_umap, draw_rank, row_id)

cluster_sizes <- s4_umap[, .N, by = .(annotation_class, recovery_cluster)]
largest_cluster <- cluster_sizes[, .SD[which.max(N)], by = annotation_class]
label_pos <- merge(
  s4_umap,
  largest_cluster[, .(annotation_class, recovery_cluster)],
  by = c("annotation_class", "recovery_cluster"), all = FALSE
)[, .(umap_1 = median(umap_1), umap_2 = median(umap_2)), by = annotation_class]
label_pos[, label := c(
  "Nonciliated\nepithelium", "Ciliated\nepithelium", "Excluded\nepithelium",
  "Stromal /\nperivascular", "Endothelial", "Immune", "Erythrocyte"
)[match(as.character(annotation_class), class_levels)]]
# Small deterministic optical offsets keep the two adjacent epithelial labels
# distinct while retaining each label within its corresponding point cloud.
label_pos[annotation_class == "Excluded epithelium", `:=`(umap_1 = 6.9, umap_2 = 0.35)]
label_pos[annotation_class == "Nonciliated epithelium", `:=`(umap_1 = 12.75, umap_2 = 1.55)]
label_pos[annotation_class == "Ciliated epithelium", `:=`(umap_1 = 5.9, umap_2 = -2.15)]
label_pos[, text_tone := fifelse(
  as.character(annotation_class) %in% c("Ciliated epithelium", "Excluded epithelium"),
  "dark", "light"
)]

pA <- ggplot(s4_umap, aes(umap_1, umap_2, colour = annotation_class)) +
  scattermore::geom_scattermore(
    pointsize = 1.55, pixels = c(1900, 1700), alpha = 0.72,
    show.legend = FALSE
  ) +
  geom_label(
    data = label_pos[text_tone == "light"],
    aes(x = umap_1, y = umap_2, label = label, fill = annotation_class),
    inherit.aes = FALSE, colour = "white", family = "sans", fontface = "bold",
    size = 2.42, lineheight = 0.92, label.padding = grid::unit(0.85, "mm"),
    label.r = grid::unit(0.8, "mm"), linewidth = 0, show.legend = FALSE
  ) +
  geom_label(
    data = label_pos[text_tone == "dark"],
    aes(x = umap_1, y = umap_2, label = label, fill = annotation_class),
    inherit.aes = FALSE, colour = "#202020", family = "sans", fontface = "bold",
    size = 2.42, lineheight = 0.92, label.padding = grid::unit(0.85, "mm"),
    label.r = grid::unit(0.8, "mm"), linewidth = 0, show.legend = FALSE
  ) +
  scale_colour_manual(values = class_colours, drop = FALSE) +
  scale_fill_manual(values = class_colours, drop = FALSE) +
  coord_fixed(clip = "off") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_pub() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    legend.position = "none",
    plot.margin = margin(4, 7, 4, 4, unit = "pt")
  )

# Panel B: the locked 7 x 7 lineage matrix as a dot heatmap. Colour encodes
# signed within-class z score and area redundantly encodes its absolute value.
lineage_levels <- c(
  "Epithelium", "Ciliated epithelium", "Stromal", "Endothelial",
  "Perivascular", "Immune", "Erythrocyte"
)
class_labels_wrapped <- c(
  "Nonciliated epithelium" = "Nonciliated epithelium",
  "Ciliated epithelium" = "Ciliated epithelium",
  "Excluded epithelium" = "Excluded epithelium",
  "Stromal / perivascular" = "Stromal / perivascular",
  "Endothelial" = "Endothelial",
  "Immune" = "Immune",
  "Erythrocyte" = "Erythrocyte"
)
lineage_labels_wrapped <- c(
  "Epithelium" = "Epithel.",
  "Ciliated epithelium" = "Ciliated",
  "Stromal" = "Stromal",
  "Endothelial" = "Endoth.",
  "Perivascular" = "Perivasc.",
  "Immune" = "Immune",
  "Erythrocyte" = "Erythro."
)
heat_long[, class_label := factor(as.character(annotation_class), levels = rev(class_levels))]
heat_long[, lineage_label := factor(as.character(lineage_label), levels = lineage_levels)]
heat_long[, abs_z := abs(within_class_z)]
heat_long[, abs_z_display := pmin(abs_z, 2)]
stopifnot(all(is.finite(heat_long$within_class_z)))

pB <- ggplot(heat_long, aes(lineage_label, class_label)) +
  geom_tile(fill = "#FAFAFA", colour = "#E5E7E9", linewidth = 0.30) +
  geom_point(
    aes(size = abs_z_display, fill = within_class_z), shape = 21,
    colour = "#5A5A5A", stroke = 0.24
  ) +
  scale_fill_gradient2(
    low = "#2B6CB0", mid = "#F7F7F7", high = "#C33D3D",
    midpoint = 0, limits = c(-2, 2), oob = squish,
    name = "Within-class z score"
  ) +
  scale_size_continuous(
    range = c(1.3, 5.4), limits = c(0, 2),
    breaks = c(0.5, 1.0, 2.0), name = "Absolute z score"
  ) +
  scale_x_discrete(labels = lineage_labels_wrapped, drop = FALSE) +
  scale_y_discrete(labels = class_labels_wrapped, drop = FALSE) +
  labs(x = "Lineage panel", y = "Annotation class") +
  theme_pub() +
  theme(
    axis.text.x = element_text(
      size = 6.5, angle = 0, hjust = 0.5, vjust = 1,
      lineheight = 0.90
    ),
    axis.text.y = element_text(size = 7.3),
    axis.ticks = element_blank(), axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom", legend.box = "vertical",
    legend.box.spacing = grid::unit(0.8, "mm"),
    legend.margin = margin(t = 0, unit = "pt"),
    legend.key.height = grid::unit(2.8, "mm"),
    legend.key.width = grid::unit(10, "mm"),
    plot.margin = margin(4, 4, 4, 7, unit = "pt")
  ) +
  guides(
    fill = guide_colourbar(
      order = 1, title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(29, "mm"), barheight = grid::unit(2.5, "mm")
    ),
    size = guide_legend(
      order = 2, title.position = "top", title.hjust = 0.5,
      nrow = 1, override.aes = list(fill = "#D7DEE5")
    )
  )

# Panel C: eligibility audit. Pre and post counts are deliberately unconnected
# so that sample-recovery counts cannot be read as a treatment-response slope.
elig_wide[, status := factor(status, levels = c("Eligible", "Below gate", "Post only"))]
patient_order <- sort(elig_wide$Patient)
elig_wide[, patient_y := match(Patient, rev(patient_order))]
elig_long <- melt(
  elig_wide,
  id.vars = c("Patient", "patient_y", "status", "complete_prepost", "eligible_min20_each_stage"),
  measure.vars = c("pre", "post"), variable.name = "stage", value.name = "cells"
)
elig_long <- elig_long[!is.na(cells)]
elig_long[, stage := factor(stage, levels = c("pre", "post"), labels = c("Pre", "Post"))]
elig_long[, y_plot := patient_y + fifelse(stage == "Pre", 0.13, -0.13)]
elig_long[, label_x := fcase(
  cells >= 1000, cells * 1.12,
  cells >= 500, cells / 1.05,
  default = cells * 1.18
)]
elig_long[, label_hjust := fifelse(cells >= 500 & cells < 1000, 1, 0)]
stopifnot(nrow(elig_long) == 17L)
stopifnot(all(is.finite(elig_long$cells)), all(elig_long$cells > 0))

x_min <- 2.2
x_max <- 2400
row_bands <- elig_wide[, .(
  xmin = x_min, xmax = x_max,
  ymin = patient_y - 0.42, ymax = patient_y + 0.42,
  status
)]

pC <- ggplot() +
  geom_rect(
    data = row_bands,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = status),
    colour = NA
  ) +
  geom_vline(xintercept = 20, linetype = "22", linewidth = 0.55, colour = "#303030") +
  geom_hline(
    yintercept = seq(1.5, 8.5, by = 1), linewidth = 0.22,
    colour = "white"
  ) +
  geom_point(
    data = elig_long,
    aes(x = cells, y = y_plot, colour = stage, shape = stage),
    size = 2.55, stroke = 0.62
  ) +
  geom_text(
    data = elig_long,
    aes(x = label_x, y = y_plot, label = cells, colour = stage, hjust = label_hjust),
    size = 2.38, family = "sans", show.legend = FALSE
  ) +
  annotate(
    "text", x = 21.8, y = 9.48, label = "Eligibility threshold: 20 cells",
    hjust = 0, vjust = 1, size = 2.38, family = "sans",
    colour = "#303030"
  ) +
  scale_x_log10(
    limits = c(x_min, x_max), breaks = c(3, 10, 20, 100, 1000),
    labels = c("3", "10", "20", "100", "1,000"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = elig_wide$patient_y,
    labels = elig_wide$Patient,
    limits = c(0.5, 9.55), expand = expansion(mult = c(0, 0))
  ) +
  scale_colour_manual(values = stage_colours, name = "Library") +
  scale_shape_manual(values = c("Pre" = 16, "Post" = 17), name = "Library") +
  scale_fill_manual(values = status_fills, name = "Pair status", drop = FALSE) +
  labs(x = "Nonciliated epithelial cells (log scale)", y = "Patient") +
  theme_pub() +
  theme(
    panel.grid.major.x = element_line(linewidth = 0.25, colour = "#D8DADD"),
    panel.grid.minor = element_blank(),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.box = "horizontal", legend.box.just = "left",
    legend.margin = margin(t = 0, unit = "pt"),
    legend.key.height = grid::unit(3.2, "mm"), legend.key.width = grid::unit(4.2, "mm"),
    plot.margin = margin(4, 7, 4, 7, unit = "pt")
  ) +
  guides(
    colour = guide_legend(order = 1, title.position = "left"),
    shape = guide_legend(order = 1, title.position = "left"),
    fill = guide_legend(order = 2, title.position = "left", nrow = 1)
  )

figure_manifest <- list(
  figureId = "Supplementary_Figure_S4",
  journalTarget = "International Journal of Women's Health candidate",
  panelRoles = c("primary", "secondary", "primary"),
  semanticColors = as.list(c(class_colours, stage_colours)),
  sharedLegends = list(),
  fontSizes = list(
    panelLabelPt = 10, axisTitlePt = 8.5, axisTickPt = 7.3,
    legendPt = 7.0, annotationMinPt = 6.2
  ),
  gridColumns = 12L,
  outerMarginMm = 3.0,
  designProvenance = list(
    plotcase_status = "NO_COMPATIBLE_RESULT",
    plotcase_queries = list(
      "single-cell UMAP direct labels colored label backgrounds annotation classes publication",
      "dot heatmap matrix signed z score diverging color bubble publication",
      "eligibility threshold dot plot paired stages without connecting lines log scale sample counts"
    ),
    plotcase_candidates = list(),
    plotcase_reuse_level = "NOT_USED",
    plotcase_reason = "Returned cases were wrong-family heatmaps, circular heatmaps, boxplots, or barplots and would distort the locked annotation and eligibility contracts.",
    plotcase_candidate_sufficiency = "NONE",
    plotcase_compatible_candidate_count = 0L,
    user_reference_status = "NOT_FOUND",
    scidraw_status = "NOT_REQUIRED_WITH_REASON",
    scidraw_reason = "All panels are quantitative data displays.",
    selected_concept = "balanced_annotation_audit_with_threshold_plot",
    visual_intent = list(
      first_read = "The seven endpoint-blind annotation classes are directly identifiable in the locked UMAP.",
      second_read = "The lineage-score matrix supports those displayed annotation classes.",
      third_read = "Only three of eight complete pairs pass the 20-cell gate, while one library is post-only."
    ),
    panel_visual_review = "MANUAL_PANEL_VISUAL_REVIEW_REQUIRED"
  )
)

fig <- zhplot_compose_irregular(
  list(pA, pB, pC),
  areas = data.frame(
    panel = 1:3,
    t = c(1, 1, 2), l = c(1, 6, 1),
    b = c(1, 1, 2), r = c(5, 12, 12)
  ),
  panel_ids = c("locked_umap", "lineage_dot_heatmap", "eligibility_threshold_plot"),
  alignment = "independent",
  rel_widths = rep(1, 12), rel_heights = c(1, 1),
  tag_size = 10, tag_band_fraction = 0.055, font_family = "sans",
  gutter_x = 0.035, gutter_y = 0.035,
  title_policy = "strip", figure_manifest = figure_manifest
)

bundle <- zhplot_render_bundle(
  fig, output_base,
  width_mm = 183, height_mm = 166, dpi = 600,
  source_file = script_path,
  input_files = unlist(paths),
  implementation_files = layout_file,
  include_tiff = FALSE
)

receipt <- list(
  status = "PASS_RENDER",
  figure = "Supplementary_Figure_S4",
  scientific_contract = list(
    cells = nrow(s4_umap), annotation_classes = length(class_levels),
    complete_pairs = sum(elig_wide$complete_prepost),
    eligible_pairs = sum(elig_wide$eligible_min20_each_stage),
    post_only_libraries = sum(!elig_wide$complete_prepost),
    eligibility_threshold_cells_per_stage = 20L,
    connecting_lines_in_panel_c = FALSE
  ),
  source_sha256 = lapply(paths, function(x) digest::digest(file = x, algo = "sha256")),
  outputs = bundle
)
jsonlite::write_json(
  receipt, paste0(output_base, ".render_receipt.json"),
  pretty = TRUE, auto_unbox = TRUE, null = "null", digits = 12
)

cat("PASS_RENDER_SUPPLEMENTARY_FIGURE_S4\n")
cat("PDF=", bundle$pdf, "\n", sep = "")
cat("PNG=", bundle$png, "\n", sep = "")
cat("SVG=", bundle$svg, "\n", sep = "")

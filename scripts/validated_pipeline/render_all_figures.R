#!/usr/bin/env Rscript

# Redesign of the locked IUA / PC-NC manuscript figure set.
# This script only reads validated result tables and frozen single-cell display
# coordinates/annotations. It does not refit models, change thresholds, rerun
# enrichment, recluster cells, or revise biological annotations.

options(stringsAsFactors = FALSE, warn = 1)

needed <- c(
  "data.table", "digest", "ggplot2", "ggrepel", "ggtext", "gridtext",
  "jsonlite", "patchwork", "ragg", "scales", "scattermore", "svglite"
)
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript render_all_figures.R OUTPUT_ROOT", call. = FALSE)
out_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
project <- normalizePath("E:/MR/QTLMR/project/YL/IUA_ESR1_PGR_WOI_publicdb", winslash = "/", mustWork = TRUE)
cv <- file.path(project, "revision_20260919", "context_validation")
revision <- file.path(project, "revision_20260919")
submit_source <- file.path(out_root, "locked_inputs")
if (!dir.exists(submit_source)) stop("ASCII staging copy of submitted source-data is missing.", call. = FALSE)
layout_file <- normalizePath("C:/Users/zh/.codex/skills/zhplot/scripts/zhplot_layout.R", winslash = "/", mustWork = TRUE)
source(layout_file)

pdf_dir <- file.path(out_root, "PDF")
svg_dir <- file.path(out_root, "SVG")
png_dir <- file.path(out_root, "PNG")
source_dir <- file.path(out_root, "source_data")
qa_dir <- file.path(out_root, "QA")
render_dir <- file.path(qa_dir, "render_work")
manifest_dir <- file.path(qa_dir, "manifests")
design_dir <- file.path(out_root, "design")
for (d in c(pdf_dir, svg_dir, png_dir, source_dir, qa_dir, render_dir, manifest_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

figure_stems <- c("Figure_1", "Figure_2", "Figure_3", paste0("Supplementary_Figure_S", 1:5))
final_targets <- unlist(lapply(figure_stems, function(x) c(
  file.path(pdf_dir, paste0(x, ".pdf")),
  file.path(svg_dir, paste0(x, ".svg")),
  file.path(png_dir, paste0(x, ".png"))
)))
if (any(file.exists(final_targets))) stop("Refusing to overwrite an existing redesign export.", call. = FALSE)

paths <- list(
  f1_donor = file.path(cv, "reframed_figures_20260920", "revised_r2", "source_data", "Figure_1_donor_pseudobulks.tsv"),
  f1_est = file.path(cv, "reframed_figures_20260920", "revised_r2", "source_data", "Figure_1_locked_endpoint_estimates.tsv"),
  f2_vec = file.path(cv, "results", "cross_context_effect_vectors_v2.tsv"),
  f2_concord = file.path(cv, "results", "cross_context_concordance_v2.tsv"),
  cross_endpoint = file.path(cv, "results", "cross_context_endpoint_comparison_v2.tsv"),
  f3_est = file.path(cv, "results", "GSE224093_endpoint_tests_bulk_v1.tsv"),
  s1_threshold = file.path(cv, "results", "PC_NC_threshold_equivalence_v2.tsv"),
  s1_compartment = file.path(cv, "reframed_figures_20260920", "pairedA_r5", "source_data", "Supplementary_Figure_1_compartment_sensitivity.tsv"),
  s1_lodo = file.path(cv, "results", "PC_NC_leave_one_donor_out_v2.tsv"),
  s2_gate = file.path(submit_source, "GSE215968_identity_gated_epithelial_sensitivity_cell_gate_summary.tsv"),
  s2_score = file.path(submit_source, "GSE215968_identity_gated_epithelial_sensitivity_candidate_program_scores.tsv"),
  s2_test = file.path(submit_source, "GSE215968_identity_gated_epithelial_sensitivity_candidate_program_tests.tsv"),
  s2_receptor = file.path(submit_source, "GSE215968_identity_gated_epithelial_sensitivity_receptor_edgeR_QL.tsv"),
  s3_score = file.path(submit_source, "GSE215968_primary_epithelial_external_E2MPA_response_scores.tsv"),
  s3_test = file.path(submit_source, "GSE215968_primary_epithelial_external_E2MPA_response_tests.tsv"),
  s3_mapping = file.path(submit_source, "GSE200805_E2MPA_vs_E2_response_module_target_mapping_audit.tsv"),
  s4_cells = file.path(cv, "gse311899_paired_extension_20260920", "interim", "recovery_annotation_cells_v3_lognorm_meanrank.tsv.gz"),
  s4_umap = file.path(cv, "gse311899_paired_extension_20260920", "interim", "recovery_annotation_umap_v3_lognorm_meanrank.tsv.gz"),
  s4_labels = file.path(cv, "gse311899_paired_extension_20260920", "metadata", "recovery_annotation_cluster_label_map_v1.tsv"),
  s4_ledger = file.path(cv, "gse311899_paired_repair_20260920", "results", "paired_extension_sample_eligibility_repaired_v1.tsv"),
  s5_values = file.path(cv, "gse311899_paired_repair_20260920", "results", "paired_extension_per_patient_endpoint_values_repaired_v1.tsv"),
  s5_roundtrip = file.path(cv, "gse311899_paired_repair_20260920", "results", "paired_extension_endpoint_roundtrip_qa_repaired_v1.tsv")
)
if (!all(file.exists(unlist(paths)))) {
  stop("Missing locked input(s): ", paste(names(paths)[!file.exists(unlist(paths))], collapse = ", "), call. = FALSE)
}

read_tsv <- function(path) fread(path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
f1_donor <- read_tsv(paths$f1_donor)
f1_est <- read_tsv(paths$f1_est)
f2_vec <- read_tsv(paths$f2_vec)
f2_concord <- read_tsv(paths$f2_concord)
cross_endpoint <- read_tsv(paths$cross_endpoint)
f3_est <- read_tsv(paths$f3_est)
s1_threshold <- read_tsv(paths$s1_threshold)
s1_compartment <- read_tsv(paths$s1_compartment)
s1_lodo <- read_tsv(paths$s1_lodo)
s2_gate <- read_tsv(paths$s2_gate)
s2_score <- read_tsv(paths$s2_score)
s2_test <- read_tsv(paths$s2_test)
s2_receptor <- read_tsv(paths$s2_receptor)
s3_score <- read_tsv(paths$s3_score)
s3_test <- read_tsv(paths$s3_test)
s3_mapping <- read_tsv(paths$s3_mapping)
s4_cells <- read_tsv(paths$s4_cells)
s4_umap <- read_tsv(paths$s4_umap)
s4_labels <- read_tsv(paths$s4_labels)
s4_ledger <- read_tsv(paths$s4_ledger)
s5_values <- read_tsv(paths$s5_values)
s5_roundtrip <- read_tsv(paths$s5_roundtrip)

# Persist the seven small submission-source inputs under an ASCII output path.
# zhplot can then fingerprint them without asking Windows R to normalize the
# Chinese submission-folder target behind the junction.
submit_keys <- c("s2_gate", "s2_score", "s2_test", "s2_receptor", "s3_score", "s3_test", "s3_mapping")
render_inputs <- paths
for (nm in submit_keys) {
  target <- file.path(source_dir, paste0("locked_input_", basename(paths[[nm]])))
  if (file.exists(target)) {
    if (!identical(digest::digest(file = paths[[nm]], algo = "sha256"), digest::digest(file = target, algo = "sha256"))) {
      stop("Existing staged input hash mismatch: ", nm, call. = FALSE)
    }
  } else if (!file.copy(paths[[nm]], target, overwrite = FALSE)) {
    stop("Could not stage locked input: ", nm, call. = FALSE)
  }
  render_inputs[[nm]] <- target
}

# Scientific/data locks.
stopifnot(nrow(f1_donor) == 16L, sum(f1_donor$CycleType == "Controlled") == 7L, sum(f1_donor$CycleType == "Natural") == 9L)
stopifnot(nrow(f1_est) == 4L, setequal(f1_est$endpoint, c("PGR_logCPM", "frozen_candidate25_score", "ESR1_logCPM", "stable_epithelial_core_score")))
stopifnot(nrow(f2_vec) == 10465L, all(c("PGR", "ESR1") %in% f2_vec$gene))
f2_main <- f2_concord[analysis == "primary_full_common_expression_universe"]
stopifnot(nrow(f2_main) == 1L, f2_main$genes == 10465L)
gse215 <- cross_endpoint[context == "GSE215968_disjoint_label_branch"]
stopifnot(nrow(gse215) == 4L, setequal(gse215$endpoint, c("PGR_logCPM", "ESR1_logCPM", "frozen_candidate25_score", "stable_epithelial_core_score")))
stopifnot(nrow(f3_est) == 4L, unique(f3_est$severe_IUA_n) == 7L, unique(f3_est$normal_control_n) == 7L)
stopifnot(setequal(unique(s1_threshold$minimum_nuclei_threshold), c(20L, 50L, 100L)))
threshold_ids <- unique(s1_threshold[, .(minimum_nuclei_threshold, retained_donor_ids, retained_donors, Controlled_n, Natural_n)])
stopifnot(nrow(threshold_ids) == 3L, uniqueN(threshold_ids$retained_donor_ids) == 1L, all(threshold_ids$retained_donors == 16L))
stopifnot(nrow(s1_lodo) == 64L, all(s1_lodo[, .N, by = endpoint]$N == 16L))
stopifnot(nrow(s2_score) == 14L, sum(s2_score$condition == "AS") == 8L, sum(s2_score$condition != "AS") == 6L)
stopifnot(nrow(s2_receptor) == 2L, setequal(s2_receptor$gene, c("PGR", "ESR1")))
portable <- s3_mapping[portable_membership %in% c("E2MPA_up_donor_consistent", "E2MPA_down_donor_consistent")]
stopifnot(nrow(portable) == 26L, sum(portable$portable_membership == "E2MPA_up_donor_consistent") == 13L,
          sum(portable$portable_membership == "E2MPA_down_donor_consistent") == 13L)
stopifnot(nrow(s4_cells) == 51018L, nrow(s4_umap) == 51018L, nrow(s4_ledger) == 17L)
stopifnot(nrow(s5_values) == 24L, uniqueN(s5_values$Patient) == 3L, uniqueN(s5_values$endpoint) == 4L)
roundtrip_status <- setNames(s5_roundtrip$status, s5_roundtrip$qa_check)
stopifnot(roundtrip_status[["endpoint_matrix_roundtrip_exact"]] == "PASS", roundtrip_status[["endpoint_long_rows"]] == "PASS")

# Labels and source tables used for plotting.
endpoint_labels <- c(
  PGR_logCPM = "PGR TMM logCPM",
  frozen_candidate25_score = "Frozen WOI programme score\n(22 of 25 genes mapped)",
  ESR1_logCPM = "ESR1 TMM logCPM",
  stable_epithelial_core_score = "Stable epithelial-core score\n(9 of 11 genes mapped)"
)
endpoint_short <- c(
  PGR_logCPM = "PGR", frozen_candidate25_score = "Frozen WOI programme",
  ESR1_logCPM = "ESR1", stable_epithelial_core_score = "Epithelial-core programme"
)
condition_colours <- c("Natural cycle" = "#2C7FB8", "Programmed cycle" = "#E67E22")
condition_shapes <- c("Natural cycle" = 22, "Programmed cycle" = 21)
context_colours <- c("Local NC control" = "#2C7FB8", "AS under HRT" = "#E67E22")

fmt_num <- function(x, digits = 2L) formatC(x, format = "f", digits = digits)
fmt_prob <- function(x) {
  ifelse(x < 0.001, formatC(x, format = "g", digits = 2), formatC(x, format = "f", digits = 3))
}
math_stat <- function(est, lo, hi, symbol = "q", p = NA_real_) {
  paste0("Delta == ", fmt_num(est), "~\"[", fmt_num(lo), ", ", fmt_num(hi), "]\"*\";\"~~italic(", symbol, ")==", fmt_prob(p))
}

theme_pub <- function(base_size = 8.2) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_blank(), plot.subtitle = element_blank(),
      axis.title = element_text(size = 8.5, colour = "#1E1E1E"),
      axis.text = element_text(size = 7.5, colour = "#1E1E1E"),
      axis.line = element_line(linewidth = 0.34, colour = "#333333"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#333333"),
      legend.title = element_blank(), legend.text = element_text(size = 7.3),
      legend.key.height = grid::unit(3.4, "mm"), legend.key.width = grid::unit(4.5, "mm"),
      panel.grid = element_blank(), strip.background = element_blank(),
      strip.text = element_text(size = 7.8, face = "bold", colour = "#222222"),
      plot.margin = margin(4, 5, 4, 5, unit = "pt")
    )
}

base_provenance <- function(brief, concept, first, second, third, user_ref = FALSE, plotcase_forest = FALSE) {
  if (plotcase_forest) {
    plot_status <- "USED"
    plot_reason <- "A compatible forest-plot case supplied style-only guidance for pale row bands, interval hierarchy, and aligned numeric columns; the locked two-panel statistical structure was retained."
    queries <- list(
      "compact two-panel forest plot aligned endpoint labels confidence intervals numeric estimate columns zero reference",
      "forest plot confidence interval numeric estimate table zero line"
    )
  } else {
    plot_status <- "NOT_REQUIRED_WITH_REASON"
    plot_reason <- "The user supplied a figure-specific scientific and structural contract; no external template was needed."
    queries <- list()
  }
  list(
    plotcase_status = plot_status, plotcase_queries = queries,
    plotcase_candidates = if (plotcase_forest) list("PlotCase effect-size forest plot") else list(),
    plotcase_reuse_level = if (plotcase_forest) "STYLE_ONLY" else "NOT_USED",
    plotcase_reason = plot_reason,
    plotcase_candidate_sufficiency = if (plotcase_forest) "LIMITED" else "NONE",
    plotcase_compatible_candidate_count = if (plotcase_forest) 1L else 0L,
    user_reference_status = if (user_ref) "FOUND" else "NOT_FOUND",
    user_reference_ids = if (user_ref) "user-subgroup-hr-ci-forest-001" else character(),
    user_reference_id = if (user_ref) "user-subgroup-hr-ci-forest-001" else NULL,
    user_reference_transferability = if (user_ref) "VISUAL_REFERENCE_ONLY" else NULL,
    reference_features_used = if (user_ref) "lightweight zero-line and numeric-column hierarchy" else character(),
    reference_adaptation_notes = if (user_ref) "Only presentation hierarchy was transferred; no statistical method or input contract was reused." else NULL,
    scidraw_status = "NOT_REQUIRED_WITH_REASON", scidraw_reason = "All panels are quantitative data displays.",
    design_brief = normalizePath(file.path(design_dir, brief), winslash = "/", mustWork = TRUE),
    selected_concept = concept,
    visual_intent = list(first_read = first, second_read = second, third_read = third),
    panel_visual_review = "MANUAL_PANEL_VISUAL_REVIEW_REQUIRED"
  )
}

manifest <- function(id, roles, colours, brief, concept, first, second, third, user_ref = FALSE, plotcase_forest = FALSE) {
  list(
    figureId = id, journalTarget = "International Journal of Women's Health candidate",
    panelRoles = roles, semanticColors = as.list(colours), sharedLegends = list(),
    fontSizes = list(panelLabelPt = 10, axisTitlePt = 8.5, axisTickPt = 7.5, legendPt = 7.3, annotationMinPt = 7.5),
    gridColumns = 12L, outerMarginMm = 3.0,
    designProvenance = base_provenance(brief, concept, first, second, third, user_ref, plotcase_forest)
  )
}

finalize_bundle <- function(fig, stem, width_mm, height_mm, inputs, man) {
  base <- file.path(render_dir, stem)
  bundle <- zhplot_render_bundle(
    fig, base, width_mm = width_mm, height_mm = height_mm, dpi = 600,
    source_file = script_path, input_files = inputs,
    implementation_files = layout_file, include_tiff = FALSE
  )
  dest <- list(
    pdf = file.path(pdf_dir, paste0(stem, ".pdf")),
    png = file.path(png_dir, paste0(stem, ".png")),
    svg = file.path(svg_dir, paste0(stem, ".svg"))
  )
  for (ext in c("pdf", "png", "svg")) {
    if (!file.rename(bundle[[ext]], dest[[ext]])) stop("Failed to move ", stem, " ", ext, call. = FALSE)
  }
  layout_json <- jsonlite::read_json(bundle$manifest, simplifyVector = FALSE)
  layout_json$render$files <- dest
  final_layout <- file.path(manifest_dir, paste0(stem, ".layout.json"))
  jsonlite::write_json(layout_json, final_layout, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = 12)
  fig_json <- jsonlite::read_json(bundle$figure_manifest, simplifyVector = FALSE)
  fig_json$outputFiles <- dest
  fig_json$layoutManifest <- final_layout
  fig_json$figureManifest <- file.path(manifest_dir, paste0(stem, ".figure_manifest.json"))
  fig_json$master_artifact <- dest$pdf
  fig_json$master_sha256 <- digest::digest(file = dest$pdf, algo = "sha256")
  fig_json$qaStatus$artifactBinding <- "HASH_BOUND_PENDING_VISUAL_REVIEW"
  final_manifest <- fig_json$figureManifest
  jsonlite::write_json(fig_json, final_manifest, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = 12)
  unlink(c(bundle$manifest, bundle$figure_manifest))
  invisible(c(dest, list(layout = final_layout, figure_manifest = final_manifest)))
}

# Figure 1 -------------------------------------------------------------------
make_f1_panel <- function(endpoint) {
  ep_id <- endpoint
  d <- f1_donor[, .(
    donor_id, context,
    value = get(ep_id)
  )]
  d[, context := factor(context, levels = c("Natural cycle", "Programmed cycle"))]
  e <- f1_est[endpoint == ep_id]
  stopifnot(nrow(e) == 1L)
  annotation <- math_stat(
    e$difference_PC_minus_NC, e$bootstrap_ci95_low, e$bootstrap_ci95_high,
    symbol = "q", p = e$exact_BH_q_within_family
  )
  ggplot(d, aes(context, value)) +
    geom_boxplot(aes(fill = context), width = 0.48, outlier.shape = NA, alpha = 0.16,
      linewidth = 0.42, colour = "#3B3B3B", show.legend = FALSE) +
    stat_summary(fun = mean, geom = "point", shape = 95, size = 7.0,
      colour = "#111111", show.legend = FALSE) +
    geom_point(aes(fill = context, shape = context),
      position = position_jitter(width = 0.085, height = 0, seed = 20260922),
      size = 2.15, stroke = 0.48, colour = "#252525", show.legend = FALSE) +
    geom_text(
      data = data.frame(x = 1.5, y = Inf, label = annotation),
      aes(x = x, y = y, label = label), inherit.aes = FALSE,
      hjust = 0.5, vjust = 1.08, size = 2.55, family = "sans", parse = TRUE
    ) +
    scale_fill_manual(values = condition_colours) +
    scale_shape_manual(values = condition_shapes) +
    scale_x_discrete(labels = c("Natural cycle" = "NC\nn = 9", "Programmed cycle" = "PC\nn = 7")) +
    scale_y_continuous(expand = expansion(mult = c(0.07, 0.24))) +
    labs(x = NULL, y = endpoint_labels[[endpoint]]) +
    coord_cartesian(clip = "off") + theme_pub()
}

f1_order <- c("PGR_logCPM", "frozen_candidate25_score", "ESR1_logCPM", "stable_epithelial_core_score")
f1_plots <- lapply(f1_order, make_f1_panel)
f1_manifest <- manifest(
  "Figure_1", c("primary", "primary", "secondary", "secondary"), condition_colours,
  "Figure_1.yaml", "distribution_2x2",
  "Four prespecified PC-NC endpoint directions are visible from all donors.",
  "Bootstrap intervals and family-specific q values qualify each panel.",
  "Programme membership and mapped-gene coverage remain explicit."
)
fig1 <- zhplot_compose_regular(
  f1_plots, ncol = 2, nrow = 2, panel_ids = paste0("endpoint_", f1_order),
  alignment = "independent", rel_widths = c(1, 1), rel_heights = c(1, 1),
  tag_size = 10, tag_band_fraction = 0.072, font_family = "sans",
  gutter_x = 0.035, gutter_y = 0.035, title_policy = "strip", figure_manifest = f1_manifest
)
finalize_bundle(fig1, "Figure_1", 183, 135, c(paths$f1_donor, paths$f1_est, file.path(design_dir, "Figure_1.yaml")), f1_manifest)

# Figure 2 -------------------------------------------------------------------
highlight <- f2_vec[gene %in% c("PGR", "ESR1")]
concord_label <- paste0(
  "10,465 shared genes\nSpearman rho = ", fmt_num(f2_main$spearman, 3),
  " (95% CI ", fmt_num(f2_main$spearman_bootstrap_ci95_low, 3), "-", fmt_num(f2_main$spearman_bootstrap_ci95_high, 3), ")",
  "\nSign concordance = ", percent(f2_main$sign_concordance, accuracy = 0.1),
  " (95% CI ", percent(f2_main$sign_concordance_bootstrap_ci95_low, accuracy = 0.1), "-",
  percent(f2_main$sign_concordance_bootstrap_ci95_high, accuracy = 0.1), ")"
)
p2a <- ggplot(f2_vec, aes(PC_minus_NC_logCPM_mean_difference, AS_HRT_minus_local_control_logCPM_mean_difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.34, colour = "#8A8A8A") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.34, colour = "#8A8A8A") +
  geom_point(size = 0.48, alpha = 0.24, colour = "#4C78A8") +
  geom_point(data = highlight, aes(colour = gene), size = 2.2, alpha = 1) +
  ggrepel::geom_label_repel(
    data = highlight, aes(label = gene, colour = gene), seed = 20260922,
    size = 2.65, family = "sans", fill = alpha("white", 0.92),
    linewidth = 0.18, label.padding = grid::unit(1.0, "mm"),
    point.padding = grid::unit(0.7, "mm"), box.padding = grid::unit(1.1, "mm"),
    min.segment.length = 0, segment.size = 0.25, show.legend = FALSE
  ) +
  annotate("text", x = -Inf, y = Inf, label = concord_label, hjust = -0.03, vjust = 1.05,
    size = 2.35, family = "sans", colour = "#252525", lineheight = 1.08) +
  scale_colour_manual(values = c(PGR = "#D55E00", ESR1 = "#0072B2")) +
  labs(
    x = "PC - NC mean logCPM difference",
    y = "AS under HRT - local NC mean logCPM difference"
  ) + coord_cartesian(clip = "off") + theme_pub() + theme(legend.position = "none")

make_gse215_panel <- function(d, xlab, x_breaks = waiver()) {
  d <- copy(d)
  d[, y := rev(seq_len(.N))]
  d[, endpoint_label := unname(endpoint_short[endpoint])]
  d[, prob_label := paste0("<i>q</i> = ", fmt_prob(exact_q), "; <i>P</i> = ", fmt_prob(exact_p))]
  xr <- range(c(0, d$ci_low, d$ci_high), finite = TRUE)
  pad <- diff(xr) * 0.14
  contrast <- "AS under HRT (n = 8) - local NC control (n = 6)"
  ggplot(d, aes(effect, y)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.34, colour = "#8A8A8A") +
    geom_segment(aes(x = ci_low, xend = ci_high, yend = y), linewidth = 0.75, colour = "#8A8A8A") +
    geom_point(size = 2.55, shape = 21, fill = "#E67E22", colour = "#262626", stroke = 0.42) +
    ggtext::geom_richtext(aes(x = Inf, label = prob_label), hjust = -0.04,
      size = 2.25, family = "sans", fill = NA, label.color = NA,
      label.padding = grid::unit(0, "pt")) +
    annotate("text", x = mean(xr), y = max(d$y) + 0.55, label = contrast,
      size = 2.15, family = "sans", colour = "#454545") +
    scale_y_continuous(breaks = d$y, labels = d$endpoint_label,
      limits = c(0.55, max(d$y) + 0.75), expand = expansion(mult = c(0, 0))) +
    scale_x_continuous(limits = c(xr[1] - pad, xr[2] + pad), breaks = x_breaks) +
    labs(x = xlab, y = NULL) + coord_cartesian(clip = "off") +
    theme_pub() + theme(plot.margin = margin(4, 61, 4, 5, unit = "pt"))
}

d2b <- gse215[endpoint %in% c("PGR_logCPM", "ESR1_logCPM")]
d2b[, order__ := match(endpoint, c("PGR_logCPM", "ESR1_logCPM"))]
setorder(d2b, order__)
d2b[, order__ := NULL]
d2c <- gse215[endpoint %in% c("frozen_candidate25_score", "stable_epithelial_core_score")]
d2c[, order__ := match(endpoint, c("frozen_candidate25_score", "stable_epithelial_core_score"))]
setorder(d2c, order__)
d2c[, order__ := NULL]
p2b <- make_gse215_panel(d2b, "Mean logCPM difference")
p2c <- make_gse215_panel(d2c, "Mean programme-score difference", x_breaks = c(-6, -3, 0))

f2_manifest <- manifest(
  "Figure_2", c("primary", "secondary", "secondary"),
  c(background_genes = "#4C78A8", PGR = "#D55E00", ESR1 = "#0072B2", endpoint = "#E67E22"),
  "Figure_2.yaml", "scatter_plus_stacked_endpoints",
  "The gene-level cloud shows moderate cross-context similarity rather than replication.",
  "PGR and programme endpoints are supported in GSE215968, whereas ESR1 is not.",
  "AS status and HRT regimen remain inseparable in the composite contrast."
)
fig2 <- zhplot_compose_irregular(
  list(p2a, p2b, p2c),
  areas = data.frame(panel = 1:3, t = c(1, 1, 2), l = c(1, 4, 4), b = c(2, 1, 2), r = c(3, 5, 5)),
  panel_ids = c("gene_effect_scatter", "receptor_endpoints", "programme_endpoints"),
  alignment = "independent", rel_widths = c(1, 1, 1, 1, 1), rel_heights = c(1, 1),
  tag_size = 10, tag_band_fraction = 0.067, font_family = "sans",
  gutter_x = 0.028, gutter_y = 0.035, title_policy = "strip", figure_manifest = f2_manifest
)
finalize_bundle(fig2, "Figure_2", 183, 128, c(paths$f2_vec, paths$f2_concord, paths$cross_endpoint, file.path(design_dir, "Figure_2.yaml")), f2_manifest)

# Figure 3 -------------------------------------------------------------------
make_table_forest <- function(d, xlab) {
  d <- copy(d)
  d[, y := rev(seq_len(.N))]
  d[, endpoint_label := unname(endpoint_short[sub("_log2_FPKM$", "_logCPM", sub("frozen_WOI25_zscore", "frozen_candidate25_score", sub("stable_epithelial_core_zscore", "stable_epithelial_core_score", endpoint)))])]
  d[endpoint == "PGR_log2_FPKM", endpoint_label := "PGR"]
  d[endpoint == "ESR1_log2_FPKM", endpoint_label := "ESR1"]
  d[endpoint == "frozen_WOI25_zscore", endpoint_label := "Frozen WOI\nprogramme"]
  d[endpoint == "stable_epithelial_core_zscore", endpoint_label := "Epithelial-core\nprogramme"]
  d[, diff_label := paste0(fmt_num(difference_severe_IUA_minus_control), " [",
    fmt_num(bootstrap_ci95_low), ", ", fmt_num(bootstrap_ci95_high), "]")]
  d[, q_label := paste0("<i>q</i> = ", fmt_prob(exact_BH_q_within_external_family))]
  xr <- range(c(0, d$bootstrap_ci95_low, d$bootstrap_ci95_high))
  pad <- max(diff(xr) * 0.10, 0.05)
  forest_lo <- xr[1] - pad
  forest_hi <- xr[2] + pad
  forest_span <- forest_hi - forest_lo
  diff_x <- forest_hi + 0.20 * forest_span
  q_x <- forest_hi + 1.57 * forest_span
  full_hi <- forest_hi + 2.12 * forest_span
  axis_y <- 0.48
  axis_breaks <- pretty(c(forest_lo, forest_hi), n = 3)
  axis_breaks <- axis_breaks[axis_breaks >= forest_lo & axis_breaks <= forest_hi]
  ggplot(d) +
    annotate("rect", xmin = forest_lo, xmax = full_hi, ymin = 2.47, ymax = 2.92,
      fill = "#EEF2F5", colour = NA) +
    annotate("rect", xmin = forest_lo, xmax = full_hi, ymin = 1.60, ymax = 2.40,
      fill = "#F7F9FA", colour = NA) +
    annotate("rect", xmin = forest_lo, xmax = full_hi, ymin = 0.60, ymax = 1.40,
      fill = "#FFFFFF", colour = NA) +
    annotate("segment", x = 0, xend = 0, y = 0.57, yend = 2.42,
      linetype = "22", linewidth = 0.42, colour = "#747B82") +
    geom_segment(aes(x = bootstrap_ci95_low, xend = bootstrap_ci95_high, y = y, yend = y),
      linewidth = 1.02, lineend = "round", colour = "#858B91") +
    geom_point(aes(x = difference_severe_IUA_minus_control, y = y), shape = 21, size = 3.05,
      fill = "#3F6E8C", colour = "white", stroke = 0.60) +
    geom_point(aes(x = difference_severe_IUA_minus_control, y = y), shape = 21, size = 3.05,
      fill = NA, colour = "#2A2F33", stroke = 0.42) +
    geom_text(aes(x = diff_x, y = y, label = diff_label), hjust = 0,
      colour = "#25292D", size = 2.58, family = "sans") +
    ggtext::geom_richtext(aes(x = q_x, y = y, label = q_label), hjust = 0,
      colour = "#25292D", size = 2.50, family = "sans", fill = NA,
      label.color = NA, label.padding = grid::unit(0, "pt")) +
    annotate("text", x = mean(c(forest_lo, forest_hi)), y = 2.69,
      label = "Effect", size = 2.75, family = "sans", fontface = "bold", colour = "#202428") +
    annotate("text", x = diff_x, y = 2.69, label = "Difference [95% CI]",
      hjust = 0, size = 2.45, family = "sans", fontface = "bold", colour = "#202428") +
    annotate("richtext", x = q_x, y = 2.69, label = "<i>q</i>", hjust = 0,
      size = 2.45, family = "sans", fontface = "bold", colour = "#202428",
      fill = NA, label.color = NA, label.padding = grid::unit(0, "pt")) +
    annotate("segment", x = forest_lo, xend = forest_hi, y = axis_y, yend = axis_y,
      linewidth = 0.42, colour = "#34383C") +
    geom_segment(data = data.frame(x = axis_breaks),
      aes(x = x, xend = x, y = axis_y, yend = axis_y - 0.07), inherit.aes = FALSE,
      linewidth = 0.36, colour = "#34383C") +
    geom_text(data = data.frame(x = axis_breaks),
      aes(x = x, y = axis_y - 0.14, label = formatC(x, format = "f", digits = 1)),
      inherit.aes = FALSE, vjust = 1, size = 2.48, family = "sans", colour = "#34383C") +
    annotate("text", x = mean(c(forest_lo, forest_hi)), y = 0.02, label = xlab,
      vjust = 1, size = 2.66, lineheight = 0.95, family = "sans", colour = "#24282C") +
    scale_y_continuous(breaks = d$y, labels = d$endpoint_label,
      limits = c(-0.34, 3.00), expand = expansion(mult = c(0, 0))) +
    scale_x_continuous(limits = c(forest_lo, full_hi), expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = "sans", base_size = 8.2) +
    theme(
      axis.text.y = element_text(size = 8.0, colour = "#24282C", hjust = 1,
        lineheight = 0.95, margin = margin(r = 7, unit = "pt")),
      plot.margin = margin(3, 7, 5, 7, unit = "pt")
    )
}

d3a <- f3_est[endpoint %in% c("PGR_log2_FPKM", "ESR1_log2_FPKM")]
d3a[, order__ := match(endpoint, c("PGR_log2_FPKM", "ESR1_log2_FPKM"))]
setorder(d3a, order__)
d3a[, order__ := NULL]
d3b <- f3_est[endpoint %in% c("frozen_WOI25_zscore", "stable_epithelial_core_zscore")]
d3b[, order__ := match(endpoint, c("frozen_WOI25_zscore", "stable_epithelial_core_zscore"))]
setorder(d3b, order__)
d3b[, order__ := NULL]
p3a <- make_table_forest(d3a, "Severe IUA - control difference\nlog2(FPKM + 1)")
p3b <- make_table_forest(d3b, "Severe IUA - control difference\nprogramme z score")
f3_manifest <- manifest(
  "Figure_3", c("primary", "primary"), c(point = "#3F6E8C", interval = "#858B91", null = "#747B82"),
  "Figure_3.yaml", "compact_table_forest_two_panel",
  "All four untimed bulk estimates are visible without significance-first colouring.",
  "Separate receptor and programme axes preserve their distinct units.",
  "Aligned numerical columns report differences, intervals, and four-endpoint BH q values.",
  user_ref = TRUE, plotcase_forest = TRUE
)
fig3 <- zhplot_compose_regular(
  list(p3a, p3b), ncol = 2, nrow = 1,
  panel_ids = c("receptor_transcripts", "programme_scores"),
  alignment = "independent", rel_widths = c(0.93, 1.07), tag_size = 10,
  tag_band_fraction = 0.095, font_family = "sans", gutter_x = 0.035,
  title_policy = "strip", figure_manifest = f3_manifest
)
finalize_bundle(fig3, "Figure_3", 183, 70, c(paths$f3_est, file.path(design_dir, "Figure_3.yaml")), f3_manifest)

# Supplementary Figure S1 ----------------------------------------------------
threshold_table <- threshold_ids[, .(
  threshold = as.character(minimum_nuclei_threshold),
  PC = Controlled_n, NC = Natural_n, Total = retained_donors,
  `Same donor IDs` = ifelse(retained_donor_ids == retained_donor_ids[1], "Yes", "No")
)]
threshold_table[, threshold := factor(threshold, levels = c("100", "50", "20"))]
threshold_long <- melt(threshold_table, id.vars = "threshold", variable.name = "column", value.name = "value")
threshold_long[, column := factor(column, levels = c("PC", "NC", "Total", "Same donor IDs"))]
pS1a <- ggplot(threshold_long, aes(column, threshold)) +
  geom_tile(aes(fill = threshold), colour = "white", linewidth = 0.8, show.legend = FALSE) +
  geom_text(aes(label = value), size = 2.75, family = "sans", colour = "#252525") +
  scale_fill_manual(values = c("100" = "#EEF4F8", "50" = "#F7F9FB", "20" = "#EEF4F8")) +
  scale_y_discrete(labels = c("100" = expression(group >= 100), "50" = expression(group >= 50), "20" = expression(group >= 20))) +
  scale_x_discrete(position = "top", labels = c("PC", "NC", "Total", "Same\ndonor IDs")) +
  labs(x = NULL, y = "Minimum nuclei") +
  theme_void(base_family = "sans", base_size = 8.2) +
  theme(
    axis.text.x = element_text(size = 7.1, face = "bold", colour = "#252525", margin = margin(b = 3), lineheight = 0.92),
    axis.text.y = element_text(size = 7.5, colour = "#252525", margin = margin(r = 4)),
    axis.title.y = element_text(size = 8.5, colour = "#252525", angle = 90, margin = margin(r = 5)),
    plot.margin = margin(5, 5, 5, 5, unit = "pt")
  )

primary_comp <- f1_est[endpoint %in% c("PGR_logCPM", "frozen_candidate25_score"), .(
  endpoint, difference_PC_minus_NC, bootstrap_ci95_low, bootstrap_ci95_high,
  exact_BH_q_within_family, compartment = "primary_nonciliated"
)]
comp <- rbindlist(list(
  primary_comp,
  s1_compartment[, .(endpoint, difference_PC_minus_NC, bootstrap_ci95_low, bootstrap_ci95_high,
    exact_BH_q_within_family, compartment)]
), use.names = TRUE)
comp <- comp[endpoint %in% c("PGR_logCPM", "frozen_candidate25_score")]
comp[, compartment_label := factor(compartment,
  levels = c("broad_epithelial_audit", "glandular_epithelial", "primary_nonciliated"),
  labels = c("Broad epithelium", "Glandular only", "Primary nonciliated"))]
comp[, endpoint_label := factor(endpoint,
  levels = c("PGR_logCPM", "frozen_candidate25_score"), labels = c("PGR", "Frozen WOI programme"))]
comp[, q_label := paste0("<i>q</i> = ", fmt_prob(exact_BH_q_within_family))]
pS1b <- ggplot(comp, aes(difference_PC_minus_NC, compartment_label, colour = compartment_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.32, colour = "#8A8A8A") +
  geom_segment(aes(x = bootstrap_ci95_low, xend = bootstrap_ci95_high, yend = compartment_label), linewidth = 0.72) +
  geom_point(size = 2.25) +
  ggtext::geom_richtext(aes(x = Inf, label = q_label), hjust = -0.05, size = 2.25,
    family = "sans", fill = NA, label.color = NA, label.padding = grid::unit(0, "pt")) +
  facet_wrap(~endpoint_label, ncol = 1, scales = "free_x") +
  scale_colour_manual(values = c("Primary nonciliated" = "#3B3B3B", "Glandular only" = "#2C7FB8", "Broad epithelium" = "#E67E22"), guide = "none") +
  labs(x = "PC - NC mean difference", y = NULL) + coord_cartesian(clip = "off") +
  theme_pub() + theme(plot.margin = margin(3, 43, 3, 5, unit = "pt"), strip.text = element_text(size = 7.2))

lodo <- copy(s1_lodo)
lodo <- merge(
  lodo,
  f1_est[, .(endpoint, full_effect = difference_PC_minus_NC)],
  by = "endpoint",
  all.x = TRUE,
  sort = FALSE
)
if (any(!is.finite(lodo$full_effect)) || any(lodo$full_effect == 0)) {
  stop("Cannot compute leave-one-donor-out relative effects from the locked full-sample estimates.", call. = FALSE)
}
lodo[, relative_effect := difference_PC_minus_NC / full_effect]
lodo[, endpoint_label := factor(endpoint, levels = rev(f1_order), labels = rev(unname(endpoint_short[f1_order])))]
pS1c <- ggplot(lodo, aes(relative_effect, endpoint_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.42, colour = "#B2182B") +
  geom_point(position = position_jitter(width = 0, height = 0.105, seed = 20260922),
    size = 1.45, alpha = 0.80, colour = "#3D596D") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  labs(x = "Leave-one-donor-out effect / full-sample effect", y = NULL) + theme_pub()

s1_manifest <- manifest(
  "Supplementary_Figure_S1", c("supporting", "primary", "primary"),
  c(primary = "#3B3B3B", glandular = "#2C7FB8", broad = "#E67E22", reference = "#B2182B"),
  "Supplementary_Figure_S1.yaml", "audit_forest_lodo",
  "The same 16 donors are eligible at all three thresholds.",
  "Compartment sensitivity retains all directions and exposes broad-PGR attenuation.",
  "All 16 leave-one-donor-out recalculations per endpoint are displayed around a ratio reference of one."
)
figS1 <- zhplot_compose_irregular(
  list(pS1a, pS1b, pS1c),
  areas = data.frame(panel = 1:3, t = c(1, 1, 2), l = c(1, 6, 1), b = c(1, 1, 2), r = c(5, 12, 12)),
  panel_ids = c("threshold_eligibility", "compartment_sensitivity", "leave_one_donor_out"),
  alignment = "independent", rel_widths = rep(1, 12), rel_heights = c(0.86, 1.14),
  tag_size = 10, tag_band_fraction = 0.075, font_family = "sans",
  gutter_x = 0.035, gutter_y = 0.045, title_policy = "strip", figure_manifest = s1_manifest
)
finalize_bundle(figS1, "Supplementary_Figure_S1", 183, 142,
  c(paths$s1_threshold, paths$s1_compartment, paths$s1_lodo, paths$f1_est, file.path(design_dir, "Supplementary_Figure_S1.yaml")), s1_manifest)

# Supplementary Figure S2 ----------------------------------------------------
gate_sum <- s2_gate[, .(candidate_cells = sum(candidate_cells), passed_cells = sum(passed_cells)), by = .(donor_id, condition)]
gate_sum <- gate_sum[donor_id %in% s2_score$donor_id]
gate_sum[, pass_fraction := passed_cells / candidate_cells]
gate_sum[, display_group := factor(ifelse(condition == "AS", "AS under HRT", "Local NC control"), levels = c("AS under HRT", "Local NC control"))]
gate_sum[, donor_label := gsub("_", "-", donor_id)]
gate_order <- gate_sum[order(display_group, pass_fraction), donor_label]
gate_sum[, donor_label := factor(donor_label, levels = gate_order)]
stopifnot(nrow(gate_sum) == 14L, table(gate_sum$display_group)[["AS under HRT"]] == 8L, table(gate_sum$display_group)[["Local NC control"]] == 6L)
pS2a <- ggplot(gate_sum, aes(pass_fraction, donor_label, fill = display_group)) +
  geom_col(width = 0.70, colour = "#333333", linewidth = 0.24) +
  facet_grid(display_group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = context_colours, guide = "none") +
  scale_x_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cells passing marker gate", y = NULL) +
  theme_pub() + theme(strip.text.y = element_text(angle = 0, size = 7.2), panel.spacing.y = grid::unit(2.5, "mm"))

s2_score[, display_group := factor(ifelse(condition == "AS", "AS under HRT", "Local NC control"), levels = c("Local NC control", "AS under HRT"))]
s2_prog_test <- s2_test[score_name == "candidate25_available"]
stopifnot(nrow(s2_prog_test) == 1L)
s2_annot <- math_stat(s2_prog_test$mean_difference_AS_minus_control, s2_prog_test$bootstrap_ci95_low,
  s2_prog_test$bootstrap_ci95_high, symbol = "P", p = s2_prog_test$exact_permutation_p_two_sided)
pS2b <- ggplot(s2_score, aes(display_group, candidate25_available)) +
  geom_boxplot(aes(fill = display_group), width = 0.48, outlier.shape = NA, alpha = 0.16,
    linewidth = 0.4, colour = "#333333", show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 95, size = 7, colour = "#111111") +
  geom_point(aes(fill = display_group), shape = 21,
    position = position_jitter(width = 0.075, height = 0, seed = 20260922),
    size = 2.2, stroke = 0.42, colour = "#222222", show.legend = FALSE) +
  geom_text(data = data.frame(x = 1.5, y = Inf, label = s2_annot),
    aes(x, y, label = label), inherit.aes = FALSE, hjust = 0.5, vjust = 1.05,
    size = 2.35, family = "sans", parse = TRUE) +
  scale_fill_manual(values = context_colours) +
  scale_x_discrete(labels = c("Local NC control" = "Local NC\nn = 6", "AS under HRT" = "AS under HRT\nn = 8")) +
  scale_y_continuous(expand = expansion(mult = c(0.07, 0.25))) +
  labs(x = NULL, y = "Frozen WOI programme score") + theme_pub()

rec <- copy(s2_receptor)
rec[, y := match(gene, c("ESR1", "PGR"))]
rec[, q_label := paste0("<i>q</i> = ", fmt_prob(target_receptor_FDR_BH2))]
pS2c <- ggplot(rec, aes(logFC, y)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#8A8A8A") +
  geom_segment(aes(x = 0, xend = logFC, yend = y), linewidth = 0.72, colour = "#B8B8B8") +
  geom_point(aes(fill = gene), shape = 21, size = 3.0, colour = "#252525", stroke = 0.44) +
  ggtext::geom_richtext(aes(x = Inf, label = q_label), hjust = -0.06, size = 2.45,
    family = "sans", fill = NA, label.color = NA, label.padding = grid::unit(0, "pt")) +
  scale_fill_manual(values = c(ESR1 = "#BDBDBD", PGR = "#D55E00"), guide = "none") +
  scale_y_continuous(breaks = c(1, 2), labels = c("ESR1", "PGR"), limits = c(0.5, 2.5)) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.10))) +
  labs(x = "edgeR logFC (AS under HRT - local NC)", y = NULL) +
  coord_cartesian(clip = "off") + theme_pub() + theme(plot.margin = margin(4, 46, 4, 5, unit = "pt"))

s2_manifest <- manifest(
  "Supplementary_Figure_S2", c("primary", "supporting", "supporting"),
  c(AS_under_HRT = "#E67E22", local_NC = "#2C7FB8", PGR = "#D55E00", ESR1 = "#BDBDBD"),
  "Supplementary_Figure_S2.yaml", "horizontal_gate_distribution_lollipop",
  "Donor-specific gate retention is visible for all 14 analysed donors.",
  "The marker-gated frozen programme remains lower in AS/HRT pseudobulks.",
  "The receptor panel retains point estimates only and labels the receptor-BH2 q values.",
  user_ref = TRUE
)
figS2 <- zhplot_compose_irregular(
  list(pS2a, pS2b, pS2c),
  areas = data.frame(panel = 1:3, t = c(1, 1, 2), l = c(1, 7, 7), b = c(2, 1, 2), r = c(6, 12, 12)),
  panel_ids = c("gate_fraction", "gated_programme", "gated_receptors"),
  alignment = "independent", rel_widths = rep(1, 12), rel_heights = c(1, 1),
  tag_size = 10, tag_band_fraction = 0.070, font_family = "sans",
  gutter_x = 0.04, gutter_y = 0.04, title_policy = "strip", figure_manifest = s2_manifest
)
finalize_bundle(figS2, "Supplementary_Figure_S2", 183, 132,
  c(render_inputs$s2_gate, render_inputs$s2_score, render_inputs$s2_test, render_inputs$s2_receptor, file.path(design_dir, "Supplementary_Figure_S2.yaml")), s2_manifest)

# Supplementary Figure S3 ----------------------------------------------------
# Keep the target-cohort result first to match the manuscript legend. The
# selected-gene key is drawn inside the empty upper-right part of panel B so
# that an external legend cannot shorten or vertically displace its axes.
s3_score[, display_group := factor(
  ifelse(condition == "AS", "AS under HRT", "Local NC control"),
  levels = c("AS under HRT", "Local NC control")
)]
s3t <- s3_test[score_name == "E2MPA_response_score_PGR_excluded"]
stopifnot(nrow(s3t) == 1L)
s3_annot <- math_stat(s3t$mean_difference_AS_minus_control, s3t$bootstrap_ci95_low,
  s3t$bootstrap_ci95_high, symbol = "P", p = s3t$exact_permutation_p_two_sided)
pS3a <- ggplot(s3_score, aes(display_group, E2MPA_response_score_PGR_excluded)) +
  geom_boxplot(aes(fill = display_group), width = 0.48, outlier.shape = NA, alpha = 0.16,
    linewidth = 0.40, colour = "#333333", show.legend = FALSE) +
  geom_point(aes(fill = display_group), shape = 21,
    position = position_jitter(width = 0.075, height = 0, seed = 20260922),
    size = 2.2, stroke = 0.42, colour = "#222222", show.legend = FALSE) +
  geom_text(data = data.frame(x = 1.5, y = Inf, label = s3_annot),
    aes(x, y, label = label), inherit.aes = FALSE, hjust = 0.5, vjust = 1.05,
    size = 2.35, family = "sans", parse = TRUE) +
  scale_fill_manual(values = context_colours) +
  scale_x_discrete(labels = c("AS under HRT" = "AS/HRT (n=8)", "Local NC control" = "Local NC (n=6)")) +
  scale_y_continuous(limits = c(0.8, 6.3), expand = expansion(mult = 0)) +
  labs(x = "Study group", y = "E2+MPA response-only score") +
  theme_pub() +
  theme(plot.margin = margin(4, 5, 4, 5, unit = "pt"))

portable[, membership := factor(
  ifelse(grepl("up", portable_membership), "E2+MPA up", "E2+MPA down"),
  levels = c("E2+MPA down", "E2+MPA up")
)]
s3_key <- data.table(
  x = c(0.40, 0.40), x_text = c(0.52, 0.52), y = c(1.08, 0.90),
  membership = factor(c("E2+MPA down", "E2+MPA up"), levels = levels(portable$membership))
)
pS3b <- ggplot(portable, aes(delta_donor1, delta_donor2, colour = membership, shape = membership)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.32, colour = "#8A8A8A") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.32, colour = "#8A8A8A") +
  geom_point(size = 2.15, alpha = 0.92) +
  geom_point(data = s3_key, aes(x = x, y = y, colour = membership, shape = membership),
    inherit.aes = FALSE, size = 2.15, alpha = 0.92, show.legend = FALSE) +
  geom_text(data = s3_key, aes(x = x_text, y = y, label = membership),
    inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 2.25,
    family = "sans", colour = "#202428") +
  scale_colour_manual(values = c("E2+MPA down" = "#D55E00", "E2+MPA up" = "#009E73"), guide = "none") +
  scale_shape_manual(values = c("E2+MPA down" = 16, "E2+MPA up" = 17), guide = "none") +
  scale_x_continuous(limits = c(-1.75, 1.58), breaks = c(-1, 0, 1), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(-1.18, 1.20), breaks = c(-1, -0.5, 0, 0.5, 1), expand = expansion(mult = 0)) +
  labs(x = "Organoid donor 1 change in log2 RPKM", y = "Organoid donor 2 change in log2 RPKM") +
  theme_pub() +
  theme(legend.position = "none", plot.margin = margin(4, 5, 4, 5, unit = "pt"))

s3_manifest <- manifest(
  "Supplementary_Figure_S3", c("primary", "primary"),
  c(module_down = "#D55E00", module_up = "#009E73", AS_under_HRT = "#E67E22", local_NC = "#2C7FB8"),
  "Supplementary_Figure_S3.yaml", "distribution_then_selected_scatter_inset_key",
  "The target-cohort response-only score does not distinguish the two groups.",
  "The 26 selected portable genes have the prespecified signed membership in both external donors.",
  "The display does not report a correlation, regression, or validation claim."
)
figS3 <- zhplot_compose_regular(
  list(pS3a, pS3b), ncol = 2, nrow = 1, panel_ids = c("target_score", "external_selected_genes"),
  alignment = "plot-region",
  alignment_reason = "The user requested identical physical x- and y-axis lengths and a shared baseline; the two scales remain independently labelled and are not interpreted as commensurate.",
  rel_widths = c(1, 1), tag_size = 10,
  tag_band_fraction = 0.085, font_family = "sans", gutter_x = 0.045,
  title_policy = "strip", figure_manifest = s3_manifest
)
finalize_bundle(figS3, "Supplementary_Figure_S3", 183, 90,
  c(render_inputs$s3_mapping, render_inputs$s3_score, render_inputs$s3_test, file.path(design_dir, "Supplementary_Figure_S3.yaml")), s3_manifest)

# Supplementary Figure S4 ----------------------------------------------------
display_class <- function(x) {
  fcase(
    x == "Nonciliated_epithelial_primary", "Nonciliated epithelium",
    x == "Ciliated_epithelial", "Ciliated epithelium",
    grepl("Mixed_epithelial|sample_dominant_excluded|Nonciliated_epithelial_sample", x), "Excluded epithelium",
    grepl("Stromal", x) | x == "Perivascular", "Stromal / perivascular",
    x == "Endothelial", "Endothelial",
    grepl("Immune", x), "Immune",
    grepl("Erythrocyte", x), "Erythrocyte",
    default = "Other / ambiguous"
  )
}
class_levels <- c(
  "Nonciliated epithelium", "Ciliated epithelium", "Excluded epithelium",
  "Stromal / perivascular", "Endothelial", "Immune", "Erythrocyte"
)
class_colours <- c(
  "Nonciliated epithelium" = "#D55E00", "Ciliated epithelium" = "#E69F00",
  "Excluded epithelium" = "#BDBDBD", "Stromal / perivascular" = "#009E73",
  "Endothelial" = "#0072B2", "Immune" = "#7B61A8", "Erythrocyte" = "#CC79A7"
)
s4_umap[, row_id__ := .I]
s4_umap <- merge(s4_umap,
  s4_labels[, .(recovery_cluster, frozen_recovery_label)], by = "recovery_cluster", all.x = TRUE, sort = FALSE)
setorder(s4_umap, row_id__)
s4_cells[, row_id__ := .I]
s4_cells <- merge(s4_cells,
  s4_labels[, .(recovery_cluster, frozen_recovery_label)], by = "recovery_cluster", all.x = TRUE, sort = FALSE)
setorder(s4_cells, row_id__)
stopifnot(!anyNA(s4_umap$frozen_recovery_label), !anyNA(s4_cells$frozen_recovery_label))
s4_umap[, annotation_class := display_class(frozen_recovery_label)]
s4_cells[, annotation_class := display_class(frozen_recovery_label)]
stopifnot(setequal(unique(s4_umap$annotation_class), class_levels), setequal(unique(s4_cells$annotation_class), class_levels))
s4_umap[, annotation_class := factor(annotation_class, levels = class_levels)]
s4_cells[, annotation_class := factor(annotation_class, levels = class_levels)]

draw_order <- c("Excluded epithelium", "Stromal / perivascular", "Endothelial", "Immune", "Erythrocyte", "Ciliated epithelium", "Nonciliated epithelium")
s4_umap[, draw_rank := match(as.character(annotation_class), draw_order)]
setorder(s4_umap, draw_rank, row_id__)
cluster_sizes <- s4_umap[, .N, by = .(annotation_class, recovery_cluster)]
largest_cluster <- cluster_sizes[, .SD[which.max(N)], by = annotation_class]
label_pos <- merge(s4_umap, largest_cluster[, .(annotation_class, recovery_cluster)],
  by = c("annotation_class", "recovery_cluster"), all = FALSE)[,
  .(umap_1 = median(umap_1), umap_2 = median(umap_2)), by = annotation_class]
label_pos[, label := c(
  "Nonciliated\nepithelium", "Ciliated\nepithelium", "Excluded\nepithelium",
  "Stromal /\nperivascular", "Endothelial", "Immune", "Erythrocyte"
)[match(as.character(annotation_class), class_levels)]]

pS4a <- ggplot(s4_umap, aes(umap_1, umap_2, colour = annotation_class)) +
  scattermore::geom_scattermore(pointsize = 1.75, pixels = c(1900, 1700), alpha = 0.76, show.legend = FALSE) +
  ggrepel::geom_label_repel(
    data = label_pos, aes(x = umap_1, y = umap_2, label = label), inherit.aes = FALSE, seed = 20260922,
    size = 2.4, family = "sans", colour = "#202020", fill = alpha("white", 0.90),
    linewidth = 0.2, label.padding = grid::unit(1.0, "mm"),
    box.padding = grid::unit(1.2, "mm"), point.padding = grid::unit(0.7, "mm"),
    min.segment.length = 0, segment.colour = "#5A5A5A", segment.size = 0.24, max.overlaps = Inf
  ) +
  scale_colour_manual(values = class_colours) + coord_fixed(clip = "off") +
  labs(x = "UMAP 1", y = "UMAP 2") + theme_pub() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")

panel_cols <- c("Epithelium", "Ciliated_epithelium", "Stromal", "Endothelial", "Perivascular", "Immune", "Erythrocyte")
stopifnot(all(panel_cols %in% names(s4_cells)))
heat <- s4_cells[, lapply(.SD, mean), by = annotation_class, .SDcols = panel_cols]
heat_long <- melt(heat, id.vars = "annotation_class", variable.name = "lineage_panel", value.name = "mean_score")
heat_long[, within_class_z := as.numeric(scale(mean_score)), by = annotation_class]
stopifnot(all(is.finite(heat_long$within_class_z)))
class_z_check <- heat_long[, .(mean_z = mean(within_class_z), sd_z = sd(within_class_z)), by = annotation_class]
stopifnot(max(abs(class_z_check$mean_z)) < 1e-12, max(abs(class_z_check$sd_z - 1)) < 1e-12)
panel_labels <- c(
  Epithelium = "Epithelium", Ciliated_epithelium = "Ciliated epithelium",
  Stromal = "Stromal", Endothelial = "Endothelial", Perivascular = "Perivascular",
  Immune = "Immune", Erythrocyte = "Erythrocyte"
)
heat_long[, lineage_label := factor(unname(panel_labels[lineage_panel]), levels = unname(panel_labels[panel_cols]))]
heat_long[, class_label := factor(annotation_class, levels = class_levels)]
pS4b <- ggplot(heat_long, aes(class_label, lineage_label, fill = within_class_z)) +
  geom_tile(colour = "white", linewidth = 0.42) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-2, 2), oob = squish, name = "Within-class\nz score") +
  labs(x = "Annotation class", y = "Lineage panel") +
  theme_pub() + theme(
    axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1, size = 6.5),
    axis.text.y = element_text(size = 6.8), axis.ticks = element_blank(), axis.line = element_blank(),
    legend.position = "bottom", legend.key.width = grid::unit(14, "mm"),
    legend.key.height = grid::unit(2.5, "mm"), legend.margin = margin(t = 0, unit = "pt")
  ) + guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

complete_ids <- s4_ledger[, .N, by = Patient][N == 2L, Patient]
eligible_ids <- sort(unique(s4_ledger[eligible_min20_each_stage == TRUE, Patient]))
stopifnot(length(complete_ids) == 8L, identical(eligible_ids, c("01-13", "01-14", "01-15")))
post_only <- setdiff(unique(s4_ledger$Patient), complete_ids)
stopifnot(identical(post_only, "01-05"))
elig_wide <- dcast(s4_ledger, Patient ~ stage, value.var = "primary_nonciliated_epithelial_cells")
elig_status <- unique(s4_ledger[, .(Patient, eligible_min20_each_stage, complete_prepost)])
elig_wide <- merge(elig_wide, elig_status, by = "Patient", all.x = TRUE)
elig_wide[, status := fcase(!complete_prepost, "Post only", eligible_min20_each_stage, "Eligible", default = "Below gate")]
setorder(elig_wide, Patient)
elig_wide[, patient_label := factor(Patient, levels = rev(Patient))]
table_cells <- rbindlist(list(
  elig_wide[, .(patient_label, column = "Pre cells", value = ifelse(is.na(pre), "NA", as.character(pre)), status)],
  elig_wide[, .(patient_label, column = "Post cells", value = ifelse(is.na(post), "NA", as.character(post)), status)],
  elig_wide[, .(patient_label, column = "Pair status", value = status, status)]
))
table_cells[, column := factor(column, levels = c("Pre cells", "Post cells", "Pair status"))]
table_cells[, row_fill := ifelse(as.integer(patient_label) %% 2L == 0L, "even", "odd")]
pS4c <- ggplot(table_cells, aes(column, patient_label)) +
  geom_tile(aes(fill = row_fill), colour = "white", linewidth = 0.55) +
  geom_text(aes(label = value, colour = status), size = 2.35, family = "sans") +
  scale_fill_manual(values = c(even = "#F1F4F6", odd = "#FFFFFF"), guide = "none") +
  scale_colour_manual(values = c("Eligible" = "#0067A5", "Below gate" = "#666666", "Post only" = "#8E1B1B"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = "Patient", caption = "Minimum: 20 nonciliated epithelial cells at each stage.") +
  theme_void(base_family = "sans", base_size = 8.2) +
  theme(
    axis.text.x = element_text(size = 7.2, face = "bold", colour = "#252525", margin = margin(b = 3)),
    axis.text.y = element_text(size = 7.0, colour = "#252525", margin = margin(r = 3)),
    axis.title.y = element_text(size = 8.2, angle = 90, colour = "#252525", margin = margin(r = 4)),
    plot.caption = element_text(size = 6.8, colour = "#4A4A4A", hjust = 0),
    plot.margin = margin(5, 5, 5, 5, unit = "pt")
  )

s4_manifest <- manifest(
  "Supplementary_Figure_S4", c("primary", "secondary", "primary"), class_colours,
  "Supplementary_Figure_S4.yaml", "umap_heatmap_eligibility_table",
  "The locked UMAP displays all 51,018 cells using endpoint-blind annotation classes.",
  "Lineage evidence is shown with the audited within-class standardization definition.",
  "The table separates eight complete pairs, three eligible pairs, and one post-only library."
)
figS4 <- zhplot_compose_irregular(
  list(pS4a, pS4b, pS4c),
  areas = data.frame(panel = 1:3, t = c(1, 1, 2), l = c(1, 8, 8), b = c(2, 1, 2), r = c(7, 12, 12)),
  panel_ids = c("locked_umap", "lineage_heatmap", "eligibility_table"),
  alignment = "independent", rel_widths = rep(1, 12), rel_heights = c(1.02, 0.98),
  tag_size = 10, tag_band_fraction = 0.065, font_family = "sans",
  gutter_x = 0.032, gutter_y = 0.04, title_policy = "strip", figure_manifest = s4_manifest
)
finalize_bundle(figS4, "Supplementary_Figure_S4", 183, 148,
  c(paths$s4_cells, paths$s4_umap, paths$s4_labels, paths$s4_ledger, file.path(design_dir, "Supplementary_Figure_S4.yaml")), s4_manifest)

# Supplementary Figure S5 ----------------------------------------------------
alias_map <- data.table(Patient = sort(unique(s5_values$Patient)), patient_alias = paste0("P", seq_len(uniqueN(s5_values$Patient))))
s5 <- merge(copy(s5_values), alias_map, by = "Patient", all.x = TRUE, sort = FALSE)
s5[, stage_label := factor(stage, levels = c("pre", "post"), labels = c("Pre", "Post"))]
patient_colours <- setNames(c("#0072B2", "#D55E00", "#009E73"), alias_map$patient_alias)
patient_shapes <- setNames(c(21, 22, 24), alias_map$patient_alias)
s5_endpoint_order <- c("PGR", "portable_WOI22_zscore", "ESR1", "portable_epithelial_core9_zscore")
s5_ylabels <- c(
  PGR = "PGR TMM logCPM", portable_WOI22_zscore = "Portable 22-gene\nWOI score",
  ESR1 = "ESR1 TMM logCPM", portable_epithelial_core9_zscore = "Portable 9-gene\nepithelial-core score"
)
make_s5_panel <- function(ep, show_legend = FALSE) {
  d <- s5[endpoint == ep]
  stopifnot(nrow(d) == 6L)
  p <- ggplot(d, aes(stage_label, value, group = patient_alias, colour = patient_alias, shape = patient_alias)) +
    geom_line(linewidth = 0.68, alpha = 0.88) +
    geom_point(aes(fill = stage_label), size = 2.65, stroke = 0.68) +
    scale_colour_manual(values = patient_colours, name = "Patient alias") +
    scale_shape_manual(values = patient_shapes, name = "Patient alias") +
    scale_fill_manual(values = c(Pre = "white", Post = "#F4B183"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))) +
    labs(x = NULL, y = s5_ylabels[[ep]]) + theme_pub()
  if (show_legend) {
    p + theme(legend.position = "bottom", legend.direction = "horizontal", legend.margin = margin(t = -2, unit = "pt")) +
      guides(colour = guide_legend(nrow = 1, override.aes = list(fill = "white", size = 2.5)), shape = "none")
  } else {
    p + theme(legend.position = "none")
  }
}
s5_plots <- lapply(seq_along(s5_endpoint_order), function(i) make_s5_panel(s5_endpoint_order[[i]], show_legend = i == 4L))
s5_manifest <- manifest(
  "Supplementary_Figure_S5", rep("primary", 4), patient_colours,
  "Supplementary_Figure_S5.yaml", "paired_2x2",
  "Only three eligible patient trajectories are shown on their original endpoint scales.",
  "Fixed patient aliases, colours, and shapes support cross-panel tracking.",
  "No pooled estimate, confidence interval, P value, q value, or recovery claim is displayed."
)
figS5 <- zhplot_compose_regular(
  s5_plots, ncol = 2, nrow = 2, panel_ids = paste0("paired_", s5_endpoint_order),
  alignment = "independent", rel_widths = c(1, 1), rel_heights = c(1, 1),
  tag_size = 10, tag_band_fraction = 0.090, font_family = "sans",
  gutter_x = 0.04, gutter_y = 0.04, title_policy = "strip", figure_manifest = s5_manifest
)
finalize_bundle(figS5, "Supplementary_Figure_S5", 183, 122,
  c(paths$s5_values, paths$s5_roundtrip, file.path(design_dir, "Supplementary_Figure_S5.yaml")), s5_manifest)

# Source-data package --------------------------------------------------------
write_src <- function(x, name) {
  out <- file.path(source_dir, name)
  fwrite(x, out, sep = "\t", quote = FALSE, na = "NA")
  out
}
source_outputs <- list(
  Figure_1_donor_values = write_src(f1_donor, "Figure_1_donor_values.tsv"),
  Figure_1_endpoint_estimates = write_src(f1_est, "Figure_1_endpoint_estimates.tsv"),
  Figure_2_effect_vectors = write_src(f2_vec, "Figure_2_effect_vectors.tsv"),
  Figure_2_concordance = write_src(f2_main, "Figure_2_concordance_summary.tsv"),
  Figure_2_GSE215968_endpoints = write_src(gse215, "Figure_2_GSE215968_endpoint_estimates.tsv"),
  Figure_3_endpoints = write_src(f3_est, "Figure_3_endpoint_estimates.tsv"),
  S1_threshold = write_src(threshold_table, "Supplementary_Figure_S1_threshold_eligibility.tsv"),
  S1_compartment = write_src(comp, "Supplementary_Figure_S1_compartment_sensitivity.tsv"),
  S1_lodo = write_src(s1_lodo, "Supplementary_Figure_S1_leave_one_donor_out.tsv"),
  S2_gate = write_src(gate_sum, "Supplementary_Figure_S2_gate_fraction.tsv"),
  S2_scores = write_src(s2_score, "Supplementary_Figure_S2_programme_scores.tsv"),
  S2_test = write_src(s2_prog_test, "Supplementary_Figure_S2_programme_test.tsv"),
  S2_receptor = write_src(s2_receptor, "Supplementary_Figure_S2_receptor_estimates.tsv"),
  S3_mapping = write_src(portable, "Supplementary_Figure_S3_external_selected_genes.tsv"),
  S3_scores = write_src(s3_score, "Supplementary_Figure_S3_target_scores.tsv"),
  S3_test = write_src(s3t, "Supplementary_Figure_S3_target_test.tsv"),
  S4_umap = write_src(s4_umap[, .(row_id = row_id__, recovery_cluster, annotation_class, umap_1, umap_2)], "Supplementary_Figure_S4_umap_classes.tsv.gz"),
  S4_heatmap = write_src(heat_long, "Supplementary_Figure_S4_lineage_heatmap.tsv"),
  S4_eligibility = write_src(elig_wide, "Supplementary_Figure_S4_eligibility_table.tsv"),
  S5_values = write_src(s5, "Supplementary_Figure_S5_paired_values.tsv"),
  S5_alias = write_src(alias_map, "Supplementary_Figure_S5_patient_aliases.tsv"),
  S5_roundtrip = write_src(s5_roundtrip, "Supplementary_Figure_S5_roundtrip_QA.tsv")
)

input_frames <- list(
  f1_donor = f1_donor, f1_est = f1_est, f2_vec = f2_vec, f2_concord = f2_concord,
  cross_endpoint = cross_endpoint, f3_est = f3_est, s1_threshold = s1_threshold,
  s1_compartment = s1_compartment, s1_lodo = s1_lodo, s2_gate = s2_gate,
  s2_score = s2_score, s2_test = s2_test, s2_receptor = s2_receptor,
  s3_score = s3_score, s3_test = s3_test, s3_mapping = s3_mapping,
  s4_cells = s4_cells, s4_umap = s4_umap, s4_labels = s4_labels, s4_ledger = s4_ledger,
  s5_values = s5_values, s5_roundtrip = s5_roundtrip
)
input_manifest <- rbindlist(lapply(names(paths), function(nm) {
  d <- input_frames[[nm]]
  data.table(
    input_id = nm, path = gsub("\\\\", "/", paths[[nm]]),
    sha256 = digest::digest(file = paths[[nm]], algo = "sha256"),
    rows = nrow(d), columns = ncol(d),
    key_fields = paste(head(names(d), 8), collapse = ";")
  )
}))
fwrite(input_manifest, file.path(source_dir, "input_manifest.tsv"), sep = "\t", quote = FALSE)

output_manifest <- rbindlist(lapply(names(source_outputs), function(nm) {
  f <- source_outputs[[nm]]
  d <- read_tsv(f)
  data.table(source_table = basename(f), role = nm, rows = nrow(d), columns = ncol(d),
    sha256 = digest::digest(file = f, algo = "sha256"))
}))
fwrite(output_manifest, file.path(source_dir, "source_data_manifest.tsv"), sep = "\t", quote = FALSE)

figure_manifest <- rbindlist(lapply(figure_stems, function(stem) {
  files <- c(
    PDF = file.path(pdf_dir, paste0(stem, ".pdf")),
    SVG = file.path(svg_dir, paste0(stem, ".svg")),
    PNG = file.path(png_dir, paste0(stem, ".png"))
  )
  rbindlist(lapply(names(files), function(fmt) data.table(
    figure = stem, format = fmt, path = normalizePath(files[[fmt]], winslash = "/", mustWork = TRUE),
    bytes = file.info(files[[fmt]])$size, sha256 = digest::digest(file = files[[fmt]], algo = "sha256")
  )))
}))
fwrite(figure_manifest, file.path(qa_dir, "figure_manifest.tsv"), sep = "\t", quote = FALSE)

heat_audit <- list(
  status = "PASS_LABEL_CORRECTION_ONLY",
  rows = "lineage panels", columns = "displayed annotation classes",
  upstream_aggregation = "arithmetic mean of each lineage-panel score across cells within each displayed annotation class",
  standardization = "z score across the seven lineage-panel means separately within each displayed annotation class",
  standardization_axis_in_display = "columns / within class",
  standardization_reference = "all seven displayed lineage panels for that class; all seven displayed classes retained",
  transpose_or_subset_after_standardization = FALSE,
  clipping = "display colour scale squished to [-2, 2] only; source z scores retained",
  midpoint = 0,
  maximum_absolute_within_class_mean_z = max(abs(class_z_check$mean_z)),
  maximum_absolute_within_class_sd_minus_one = max(abs(class_z_check$sd_z - 1)),
  action = "Retain numerical matrix; replace Row z score with Within-class z score and correct legend/caption."
)
jsonlite::write_json(heat_audit, file.path(qa_dir, "S4C_standardization_audit.json"), pretty = TRUE, auto_unbox = TRUE, digits = 12)

writeLines(capture.output(sessionInfo()), file.path(qa_dir, "R_session_info.txt"), useBytes = TRUE)
writeLines(c(
  "PASS_RENDER_ALL_FIGURES",
  paste0("OUTPUT_ROOT=", out_root),
  paste0("FIGURES=", length(figure_stems)),
  paste0("S4_UMAP_CELLS=", nrow(s4_umap)),
  paste0("S4_COMPLETE_PAIRS=", length(complete_ids)),
  paste0("S4_ELIGIBLE_PAIRS=", length(eligible_ids)),
  paste0("S5_LONG_ROWS=", nrow(s5_values)),
  paste0("S4C_MAX_ABS_CLASS_MEAN_Z=", format(max(abs(class_z_check$mean_z)), scientific = TRUE))
), file.path(qa_dir, "render_receipt.txt"), useBytes = TRUE)

cat("PASS_RENDER_ALL_FIGURES\n")
cat("OUTPUT_ROOT=", out_root, "\n", sep = "")
cat("FIGURES=", length(figure_stems), "\n", sep = "")

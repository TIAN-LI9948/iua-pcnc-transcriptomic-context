#!/usr/bin/env Rscript

# Portable reproduction of the donor-level primary PC/NC endpoint analysis,
# Figure 1, and the numerical rows reported in Table 2. No target-dataset
# gene selection or threshold tuning is performed here.

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
output_dir <- if (length(args) >= 1L) args[[1]] else file.path(repo_root, "reproduced_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

primary_dir <- file.path(repo_root, "data", "primary")
derived_dir <- file.path(repo_root, "data", "derived")

metadata <- fread(file.path(primary_dir, "PC_NC_nonciliated_epithelial_pseudobulk_metadata.tsv"), data.table = FALSE)
count_table <- fread(file.path(primary_dir, "PC_NC_nonciliated_epithelial_pseudobulk_counts.tsv.gz"), data.table = FALSE)
frozen_table <- fread(file.path(primary_dir, "reference_WOI_program_epithelial.tsv"), data.table = FALSE)
core_table <- fread(file.path(primary_dir, "reference_epithelial_WOI_program_stability_core_cell_specificity.tsv"), data.table = FALSE)

ids <- metadata$pseudobulk_id
stopifnot(!anyDuplicated(ids), all(ids %in% names(count_table)), !anyDuplicated(count_table[[1]]))
counts <- as.matrix(count_table[, ids, drop = FALSE])
rownames(counts) <- count_table[[1]]
storage.mode(counts) <- "integer"
stopifnot(all(counts >= 0L), all(colnames(counts) == ids))

condition <- factor(metadata$CycleType, levels = c("Natural", "Controlled"))
stopifnot(identical(as.integer(table(condition)), c(9L, 7L)))

frozen_genes <- frozen_table$gene
core_flag <- as.character(core_table$epithelial_enriched_2fold_reference_rule)
core_genes <- core_table$gene[tolower(core_flag) %in% c("true", "t", "1", "yes")]

y <- DGEList(counts = counts)
y <- normLibSizes(y, method = "TMMwsp")
logcpm <- cpm(y, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)

mapped_frozen <- intersect(frozen_genes, rownames(logcpm))
mapped_core <- intersect(core_genes, rownames(logcpm))
stopifnot(length(frozen_genes) == 25L, length(mapped_frozen) == 22L)
stopifnot(length(core_genes) == 11L, length(mapped_core) == 9L)
stopifnot(all(c("ESR1", "PGR") %in% rownames(logcpm)))

scores <- data.frame(
  metadata,
  condition = condition,
  PGR_logCPM = as.numeric(logcpm["PGR", ]),
  frozen_candidate25_score = colMeans(logcpm[mapped_frozen, , drop = FALSE]),
  ESR1_logCPM = as.numeric(logcpm["ESR1", ]),
  stable_epithelial_core_score = colMeans(logcpm[mapped_core, , drop = FALSE]),
  check.names = FALSE
)

exact_stat <- function(x, group, seed, boot_reps = 20000L) {
  pc <- which(group == "Controlled")
  nc <- which(group == "Natural")
  observed <- mean(x[pc]) - mean(x[nc])
  allocations <- combn(seq_along(x), length(pc))
  null <- apply(allocations, 2L, function(z) mean(x[z]) - mean(x[-z]))
  set.seed(seed)
  bootstrap <- replicate(
    boot_reps,
    mean(sample(x[pc], length(pc), replace = TRUE)) -
      mean(sample(x[nc], length(nc), replace = TRUE))
  )
  data.frame(
    PC_n = length(pc), NC_n = length(nc),
    PC_mean = mean(x[pc]), NC_mean = mean(x[nc]),
    difference_PC_minus_NC = observed,
    bootstrap_ci95_low = as.numeric(quantile(bootstrap, 0.025, names = FALSE)),
    bootstrap_ci95_high = as.numeric(quantile(bootstrap, 0.975, names = FALSE)),
    exact_permutation_p_two_sided = mean(abs(null) >= abs(observed) - 1e-12),
    permutation_assignments = ncol(allocations),
    bootstrap_replicates = boot_reps,
    stringsAsFactors = FALSE
  )
}

endpoint_order <- c(
  "PGR_logCPM", "frozen_candidate25_score",
  "ESR1_logCPM", "stable_epithelial_core_score"
)
families <- c("Primary", "Primary", "Secondary", "Secondary")
tests <- rbindlist(lapply(seq_along(endpoint_order), function(i) {
  z <- exact_stat(scores[[endpoint_order[[i]]]], condition, 202609191L + i)
  z$endpoint <- endpoint_order[[i]]
  z$family <- families[[i]]
  z
}), fill = TRUE)
tests[, BH_q_within_family := ave(
  exact_permutation_p_two_sided, family,
  FUN = function(z) p.adjust(z, method = "BH")
)]
tests[, mapped_genes := c(NA_integer_, length(mapped_frozen), NA_integer_, length(mapped_core))]
tests[, frozen_membership := c(NA_integer_, length(frozen_genes), NA_integer_, length(core_genes))]

fwrite(scores, file.path(output_dir, "primary_endpoint_scores.tsv"), sep = "\t")
fwrite(tests, file.path(output_dir, "primary_endpoint_tests.tsv"), sep = "\t")

endpoint_labels <- c(
  PGR_logCPM = "PGR logCPM",
  frozen_candidate25_score = "Frozen WOI score, 22/25 genes",
  ESR1_logCPM = "ESR1 logCPM",
  stable_epithelial_core_score = "Stable epithelial-core score, 9/11 genes"
)
interpretation <- c(
  PGR_logCPM = "Higher PGR transcript abundance in PC",
  frozen_candidate25_score = "Lower frozen program score in PC",
  ESR1_logCPM = "Higher ESR1 transcript abundance in PC",
  stable_epithelial_core_score = "Lower stable-core score in PC"
)

primary_table_rows <- data.frame(
  Context = "Primary PC-NC",
  Endpoint = unname(endpoint_labels[tests$endpoint]),
  Family = tests$family,
  Groups = "7 PC; 9 NC",
  Difference = tests$difference_PC_minus_NC,
  CI_low = tests$bootstrap_ci95_low,
  CI_high = tests$bootstrap_ci95_high,
  Exact_P = tests$exact_permutation_p_two_sided,
  BH_q = tests$BH_q_within_family,
  Interpretation = unname(interpretation[tests$endpoint]),
  stringsAsFactors = FALSE
)

context <- fread(file.path(derived_dir, "cross_context_endpoint_comparison.tsv"), data.table = FALSE)
context <- context[context$context == "GSE215968_disjoint_label_branch", , drop = FALSE]
context <- context[match(endpoint_order, context$endpoint), , drop = FALSE]
stopifnot(nrow(context) == 4L, !anyNA(context$endpoint))
context_family <- c(
  PGR_logCPM = "Contextual receptor family",
  frozen_candidate25_score = "Contextual program family",
  ESR1_logCPM = "Contextual receptor family",
  stable_epithelial_core_score = "Contextual program family"
)
context_table_rows <- data.frame(
  Context = "GSE215968 AS/HRT-local NC",
  Endpoint = unname(c(
    PGR_logCPM = "PGR logCPM",
    frozen_candidate25_score = "Frozen WOI score",
    ESR1_logCPM = "ESR1 logCPM",
    stable_epithelial_core_score = "Stable epithelial-core score"
  )[context$endpoint]),
  Family = unname(context_family[context$endpoint]),
  Groups = "8 AS/HRT; 6 local NC",
  Difference = context$effect,
  CI_low = context$ci_low,
  CI_high = context$ci_high,
  Exact_P = context$exact_p,
  BH_q = context$exact_q,
  Interpretation = ifelse(
    context$endpoint == "ESR1_logCPM",
    "Not supported",
    "Composite disease-and-regimen contrast"
  ),
  stringsAsFactors = FALSE
)
fwrite(rbind(primary_table_rows, context_table_rows), file.path(output_dir, "Table_2_endpoint_estimates.tsv"), sep = "\t")

PC_ORANGE <- "#D55E00"
NC_BLUE <- "#0072B2"
NEUTRAL_DARK <- "#333333"

theme_submission <- function(base_size = 9.2) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks.length = grid::unit(1.35, "mm"),
      axis.text = element_text(size = base_size - 0.55, colour = "black"),
      axis.title = element_text(size = base_size, colour = "black"),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      legend.position = "none", panel.grid = element_blank(),
      plot.margin = margin(7, 6, 5, 7)
    )
}

format_q_expr <- function(x) {
  if (x < 0.001) return("italic(q) < 0.001")
  paste0("italic(q) == ", formatC(x, format = "f", digits = 3L))
}

scores$context <- factor(
  ifelse(scores$CycleType == "Controlled", "Programmed cycle", "Natural cycle"),
  levels = c("Programmed cycle", "Natural cycle")
)
q_lookup <- setNames(tests$BH_q_within_family, tests$endpoint)

small_n_panel <- function(value_column, y_label, seed) {
  d <- data.frame(context = scores$context, value = scores[[value_column]])
  means <- aggregate(value ~ context, data = d, FUN = mean)
  span <- diff(range(d$value, na.rm = TRUE)); if (!is.finite(span) || span <= 0) span <- 1
  y_annot <- max(d$value, na.rm = TRUE) + 0.15 * span
  ggplot(d, aes(context, value, fill = context, shape = context)) +
    geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.15,
                 colour = NEUTRAL_DARK, linewidth = 0.42, show.legend = FALSE) +
    geom_point(position = position_jitter(width = 0.075, height = 0, seed = seed),
               size = 2.45, colour = NEUTRAL_DARK, stroke = 0.45, show.legend = FALSE) +
    geom_crossbar(data = means, aes(x = context, y = value, ymin = value, ymax = value),
                  inherit.aes = FALSE, width = 0.28, linewidth = 0.72,
                  colour = "black", show.legend = FALSE) +
    annotate("text", x = 1.5, y = y_annot, label = format_q_expr(q_lookup[[value_column]]),
             parse = TRUE, size = 2.9, family = "Helvetica") +
    scale_fill_manual(values = c("Programmed cycle" = PC_ORANGE, "Natural cycle" = NC_BLUE)) +
    scale_shape_manual(values = c("Programmed cycle" = 21, "Natural cycle" = 22)) +
    scale_x_discrete(labels = c("Programmed cycle" = "Programmed\n(n = 7)",
                                "Natural cycle" = "Natural\n(n = 9)")) +
    scale_y_continuous(expand = expansion(mult = c(0.07, 0.21))) +
    labs(x = NULL, y = y_label) + theme_submission() +
    theme(plot.margin = margin(10, 7, 5, 7))
}

panels <- list(
  small_n_panel("PGR_logCPM", "PGR TMM logCPM", 202609211L),
  small_n_panel("frozen_candidate25_score", "Frozen 22-gene WOI score", 202609212L),
  small_n_panel("ESR1_logCPM", "ESR1 TMM logCPM", 202609213L),
  small_n_panel("stable_epithelial_core_score", "Stable 9-gene epithelial-core score", 202609214L)
)
figure_1 <- plot_grid(
  plotlist = panels, labels = c("A", "B", "C", "D"), ncol = 2,
  align = "none", label_x = 0.012, label_y = 0.995,
  hjust = 0, vjust = 1, label_size = 13, label_fontface = "bold",
  label_fontfamily = "Helvetica"
)

width_in <- 183 / 25.4; height_in <- 132 / 25.4
cairo_pdf(file.path(output_dir, "Figure_1.pdf"), width = width_in, height = height_in,
          family = "Helvetica", bg = "white")
print(figure_1)
dev.off()
ggsave(file.path(output_dir, "Figure_1.png"), figure_1, width = width_in,
       height = height_in, units = "in", dpi = 600, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

cat("Reproduced donor endpoint scores, exact tests, Table 2 rows, and Figure 1 in:\n", output_dir, "\n", sep = "")

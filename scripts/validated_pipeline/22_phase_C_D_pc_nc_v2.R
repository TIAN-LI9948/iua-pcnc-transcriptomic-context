#!/usr/bin/env Rscript
# Phase C-D v2: donor-level nonciliated epithelial analysis and fixed robustness checks.
suppressPackageStartupMessages({library(edgeR); library(data.table); library(jsonlite)})

root <- normalizePath("E:/MR/QTLMR/project/YL/IUA_ESR1_PGR_WOI_publicdb/revision_20260919/context_validation", winslash="/", mustWork=TRUE)
project <- dirname(dirname(root))
out <- file.path(root,"results"); meta_dir <- file.path(root,"metadata")
frozen_path <- file.path(project,"results","reference_WOI_program_epithelial.tsv")
core_path <- file.path(project,"results","reference_epithelial_WOI_program_stability_core_cell_specificity.tsv")
primary_branch <- "nonciliated_epithelial"
branches <- c("nonciliated_epithelial","glandular_epithelial","broad_epithelial_audit")
seed_base <- 202609191L
boot_n <- 20000L

required <- c("PC_NC_primary_endpoint_tests_v2.tsv","PC_NC_secondary_endpoint_tests_v2.tsv","PC_NC_nonciliated_epithelial_edgeR_QL_all_genes_v2.tsv","PC_NC_leave_one_donor_out_v2.tsv","PC_NC_threshold_equivalence_v2.tsv","PC_NC_compartment_sensitivity_v2.tsv","PC_NC_gene_contributions_v2.tsv","PC_NC_leave_one_gene_out_v2.tsv","PC_NC_v1_BH_family_correction_v2.tsv")
if(any(file.exists(file.path(out,required)))) stop("Refusing to overwrite v2 Phase C-D results.")

read_branch <- function(branch) {
  m <- fread(file.path(meta_dir,paste0("PC_NC_",branch,"_pseudobulk_metadata_v2.tsv")),data.table=FALSE)
  x <- fread(file.path(root,"interim",paste0("PC_NC_",branch,"_pseudobulk_counts_v2.tsv.gz")),data.table=FALSE)
  ids <- m$pseudobulk_id
  stopifnot(!anyDuplicated(ids), all(ids %in% names(x)), !anyDuplicated(x[[1]]))
  counts <- as.matrix(x[,ids,drop=FALSE]); rownames(counts) <- x[[1]]
  if(any(counts < 0) || any(abs(counts-round(counts)) > 1e-8)) stop("Pseudobulk is not non-negative integer counts.")
  storage.mode(counts) <- "integer"
  m$condition <- factor(m$CycleType, levels=c("Natural","Controlled"))
  stopifnot(identical(as.integer(table(m$condition)),c(9L,7L)), all(colnames(counts)==m$pseudobulk_id))
  list(meta=m,counts=counts)
}

frozen <- fread(frozen_path,data.table=FALSE)$gene
coretab <- fread(core_path,data.table=FALSE)
core <- coretab$gene[as.logical(coretab$epithelial_enriched_2fold_reference_rule)]

make_logcpm <- function(counts) {
  y <- normLibSizes(DGEList(counts=counts),method="TMMwsp")
  cpm(y,log=TRUE,prior.count=2,normalized.lib.sizes=TRUE)
}

exact_stat <- function(x, condition, seed, boot_reps=boot_n) {
  controlled <- which(condition=="Controlled"); natural <- which(condition=="Natural")
  if(length(controlled)<3 || length(natural)<3) stop("At least three donors per condition are required.")
  observed <- mean(x[controlled])-mean(x[natural])
  allocations <- combn(seq_along(x),length(controlled))
  null <- apply(allocations,2,function(z) mean(x[z])-mean(x[-z]))
  set.seed(seed)
  boot <- replicate(boot_reps, mean(sample(x[controlled],length(controlled),replace=TRUE))-mean(sample(x[natural],length(natural),replace=TRUE)))
  data.frame(
    Controlled_n=length(controlled),Natural_n=length(natural),Controlled_mean=mean(x[controlled]),Natural_mean=mean(x[natural]),
    difference_PC_minus_NC=observed,bootstrap_ci95_low=as.numeric(quantile(boot,.025,names=FALSE)),bootstrap_ci95_high=as.numeric(quantile(boot,.975,names=FALSE)),
    exact_permutation_p_two_sided=mean(abs(null)>=abs(observed)-1e-12),permutation_assignments=ncol(allocations),bootstrap_replicates=boot_reps,
    stringsAsFactors=FALSE)
}

make_scores <- function(obj) {
  logcpm <- make_logcpm(obj$counts)
  mapped_frozen <- intersect(frozen,rownames(logcpm)); mapped_core <- intersect(core,rownames(logcpm))
  if(length(mapped_frozen)<20L || length(mapped_core)<ceiling(.8*length(core))) stop("Frozen outcome coverage failed.")
  data.frame(obj$meta,
    ESR1_logCPM=logcpm["ESR1",], PGR_logCPM=logcpm["PGR",],
    frozen_candidate25_score=colMeans(logcpm[mapped_frozen,,drop=FALSE]),
    stable_epithelial_core_score=colMeans(logcpm[mapped_core,,drop=FALSE]),
    check.names=FALSE)
}

test_endpoints <- function(scores, endpoints=c("PGR_logCPM","frozen_candidate25_score","ESR1_logCPM","stable_epithelial_core_score"), seed_offset=0L) {
  fam <- ifelse(endpoints %in% c("PGR_logCPM","frozen_candidate25_score"),"primary","secondary")
  ans <- lapply(seq_along(endpoints),function(i) {
    z <- exact_stat(scores[[endpoints[i]]],scores$condition,seed_base+seed_offset+i)
    z$endpoint <- endpoints[i]; z$endpoint_family <- fam[i]; z
  })
  ans <- rbindlist(ans,fill=TRUE)
  ans$exact_BH_q_within_family <- ave(ans$exact_permutation_p_two_sided,ans$endpoint_family,FUN=function(z)p.adjust(z,"BH"))
  as.data.frame(ans)
}

# Main compartment and edgeR QL support.
main_obj <- read_branch(primary_branch)
main_scores <- make_scores(main_obj)
fwrite(main_scores,file.path(out,"PC_NC_nonciliated_epithelial_endpoint_scores_v2.tsv"),sep="\t")
main_tests <- test_endpoints(main_scores,seed_offset=0L)
main_tests$compartment <- primary_branch; main_tests$minimum_nuclei_threshold <- 20L
main_tests$frozen_genes <- ifelse(main_tests$endpoint=="frozen_candidate25_score",length(frozen),ifelse(main_tests$endpoint=="stable_epithelial_core_score",length(core),NA))
main_tests$mapped_genes <- ifelse(main_tests$endpoint=="frozen_candidate25_score",length(intersect(frozen,rownames(main_obj$counts))),ifelse(main_tests$endpoint=="stable_epithelial_core_score",length(intersect(core,rownames(main_obj$counts))),NA))

design <- model.matrix(~condition,data=main_obj$meta)
y <- DGEList(counts=main_obj$counts)
keep <- filterByExpr(y,design)
yd <- normLibSizes(y[keep,,keep.lib.sizes=FALSE],method="TMMwsp")
yd <- estimateDisp(yd,design,robust=TRUE)
fit <- glmQLFit(yd,design,robust=TRUE)
qlf <- glmQLFTest(fit,coef="conditionControlled")
de <- topTags(qlf,n=Inf,sort.by="none")$table; de$gene <- rownames(de)
de$edgeR_FDR_scope <- "all_filterByExpr_retained_genes"
fwrite(de,file.path(out,"PC_NC_nonciliated_epithelial_edgeR_QL_all_genes_v2.tsv"),sep="\t")
main_tests$gene <- sub("_logCPM$","",main_tests$endpoint)
main_tests$edgeR_logFC <- de$logFC[match(main_tests$gene,de$gene)]
main_tests$edgeR_PValue <- de$PValue[match(main_tests$gene,de$gene)]
main_tests$edgeR_gene_wide_FDR <- de$FDR[match(main_tests$gene,de$gene)]
rec_idx <- which(main_tests$endpoint %in% c("PGR_logCPM","ESR1_logCPM"))
main_tests$edgeR_BH_q_receptor_support <- NA_real_
main_tests$edgeR_BH_q_receptor_support[rec_idx] <- p.adjust(main_tests$edgeR_PValue[rec_idx],"BH")
primary <- subset(main_tests,endpoint_family=="primary")
secondary <- subset(main_tests,endpoint_family=="secondary")
fwrite(primary,file.path(out,"PC_NC_primary_endpoint_tests_v2.tsv"),sep="\t")
fwrite(secondary,file.path(out,"PC_NC_secondary_endpoint_tests_v2.tsv"),sep="\t")

# Fixed threshold equivalence and glandular/broad compartment sensitivities.
thresholds <- c(20L,50L,100L)
threshold_rows <- list(); sensitivity_rows <- list(); ii <- 0L; jj <- 0L
for(branch in branches) {
  obj <- read_branch(branch); scores <- make_scores(obj)
  for(thr in thresholds) {
    retain <- obj$meta$nuclei >= thr
    if(any(table(obj$meta$condition[retain]) < 3L)) stop(paste("Threshold makes a group ineligible:",branch,thr))
    ss <- scores[retain,,drop=FALSE]
    tests <- test_endpoints(ss,seed_offset=1000L+ii); ii <- ii+10L
    tests$compartment <- branch; tests$minimum_nuclei_threshold <- thr
    tests$retained_donor_ids <- paste(ss$donor_id,collapse=";")
    tests$retained_donors <- nrow(ss)
    tests$threshold_equivalent_to_20 <- identical(sort(ss$donor_id),sort(scores$donor_id))
    threshold_rows[[length(threshold_rows)+1]] <- tests
    if(thr==20L) { sensitivity_rows[[length(sensitivity_rows)+1]] <- tests }
  }
}
threshold_table <- rbindlist(threshold_rows,fill=TRUE)
fwrite(threshold_table,file.path(out,"PC_NC_threshold_equivalence_v2.tsv"),sep="\t")
fwrite(rbindlist(sensitivity_rows,fill=TRUE),file.path(out,"PC_NC_compartment_sensitivity_v2.tsv"),sep="\t")

# Leave-one-donor-out for the primary compartment, retaining fixed endpoints and fixed test families.
lodo <- list(); kk <- 0L
for(d in main_scores$donor_id) {
  s <- main_scores[main_scores$donor_id != d,,drop=FALSE]
  z <- test_endpoints(s,seed_offset=2000L+kk); kk <- kk+10L
  z$left_out_donor <- d; z$compartment <- primary_branch; z$minimum_nuclei_threshold <- 20L
  lodo[[length(lodo)+1]] <- z
}
lodo <- rbindlist(lodo,fill=TRUE)
fwrite(lodo,file.path(out,"PC_NC_leave_one_donor_out_v2.tsv"),sep="\t")

# Score contributions and leave-one-gene-out are diagnostics, never score reselection.
logcpm_main <- make_logcpm(main_obj$counts); mapped_frozen <- intersect(frozen,rownames(logcpm_main))
condition <- main_scores$condition; controlled <- condition=="Controlled"; natural <- condition=="Natural"
gene_effect <- rowMeans(logcpm_main[mapped_frozen,controlled,drop=FALSE])-rowMeans(logcpm_main[mapped_frozen,natural,drop=FALSE])
contrib <- data.frame(gene=mapped_frozen,individual_gene_difference_PC_minus_NC=as.numeric(gene_effect),score_contribution=as.numeric(gene_effect)/length(mapped_frozen),stringsAsFactors=FALSE)
contrib$score_difference_PC_minus_NC <- mean(gene_effect)
contrib$contribution_sum_matches_score_difference <- abs(sum(contrib$score_contribution)-contrib$score_difference_PC_minus_NC[1]) < 1e-12
fwrite(contrib,file.path(out,"PC_NC_gene_contributions_v2.tsv"),sep="\t")
logo <- lapply(seq_along(mapped_frozen),function(i) {
  keep_genes <- mapped_frozen[-i]
  score <- colMeans(logcpm_main[keep_genes,,drop=FALSE])
  z <- exact_stat(score,condition,seed_base+4000L+i)
  z$left_out_gene <- mapped_frozen[i]; z$genes_retained <- length(keep_genes); z$compartment <- primary_branch; z
})
logo <- rbindlist(logo,fill=TRUE); logo$exact_BH_q_diagnostic <- p.adjust(logo$exact_permutation_p_two_sided,"BH")
fwrite(logo,file.path(out,"PC_NC_leave_one_gene_out_v2.tsv"),sep="\t")

# Correct the historical broad-output q labels without changing its historical files.
old_rec <- fread(file.path(out,"PC_NC_broad_epithelial_receptor_effects.tsv"),data.table=FALSE)
old_prog <- fread(file.path(out,"PC_NC_broad_epithelial_programme_tests.tsv"),data.table=FALSE)
old <- rbind(
  data.frame(source_endpoint=old_rec$endpoint, endpoint=old_rec$gene, exact_p=old_rec$exact_permutation_p_two_sided, endpoint_family=ifelse(old_rec$gene=="PGR","primary","secondary")),
  data.frame(source_endpoint=old_prog$endpoint, endpoint=ifelse(old_prog$endpoint=="frozen_candidate25_score","frozen_candidate25_score","stable_epithelial_core_score"), exact_p=old_prog$exact_permutation_p_two_sided, endpoint_family=ifelse(old_prog$endpoint=="frozen_candidate25_score","primary","secondary"))
)
old$corrected_BH_q_within_v2_family <- ave(old$exact_p,old$endpoint_family,FUN=function(z)p.adjust(z,"BH"))
old$interpretation <- "historical_broad_compartment_audit_only; corrected q family does not make it the v2 primary result"
fwrite(old,file.path(out,"PC_NC_v1_BH_family_correction_v2.tsv"),sep="\t")

qa <- list(
  analysis_id="IUA_CONTEXT_REPAIR_V2", completed_at=format(Sys.time(),tz="Asia/Shanghai",usetz=TRUE),
  inference_unit="donor", primary_compartment=primary_branch, excluded_from_primary="ciliated epithelial cell; epithelial cell of cervix",
  group_counts=as.list(table(main_obj$meta$condition)), threshold_values=as.list(thresholds), all_primary_donors_eligible_by_observation=all(main_obj$meta$nuclei>=max(thresholds)),
  exact_primary_family=c("PGR_logCPM","frozen_candidate25_score"), exact_secondary_family=c("ESR1_logCPM","stable_epithelial_core_score"),
  primary_permutations=choose(16,7), bootstrap_replicates=boot_n, random_seed=seed_base,
  normalization="TMMwsp", edger_version=as.character(packageVersion("edgeR")), edger_ql="filterByExpr -> normLibSizes(TMMwsp) -> estimateDisp(robust=TRUE) -> glmQLFit(robust=TRUE) -> glmQLFTest",
  filtered_genes=sum(keep), marker_panel_overlap_checked_in_phase_B=TRUE)
write_json(qa,file.path(meta_dir,"PC_NC_phase_C_D_analysis_qa_v2.json"),auto_unbox=TRUE,pretty=TRUE)
cat(toJSON(list(status="C_D_complete",primary_permutations=choose(16,7),filtered_genes=sum(keep)),auto_unbox=TRUE),"\n")

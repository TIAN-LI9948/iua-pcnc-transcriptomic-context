#!/usr/bin/env Rscript
# Phase E v2: bounded cross-context comparison, not a disease-effect estimate.
suppressPackageStartupMessages({library(edgeR); library(data.table); library(jsonlite)})
root <- normalizePath("E:/MR/QTLMR/project/YL/IUA_ESR1_PGR_WOI_publicdb/revision_20260919/context_validation",winslash="/",mustWork=TRUE)
project <- dirname(dirname(root)); rev <- dirname(root); out <- file.path(root,"results"); mdir <- file.path(root,"metadata")
need <- c("cross_context_endpoint_comparison_v2.tsv","cross_context_comparability_v2.tsv","cross_context_gene_universe_v2.tsv","cross_context_effect_vectors_v2.tsv","cross_context_concordance_v2.tsv","cross_context_bootstrap_qa_v2.json")
if(any(file.exists(file.path(out,need))) || file.exists(file.path(mdir,"cross_context_bootstrap_qa_v2.json"))) stop("Refusing to overwrite Phase E v2 output.")
pc_m <- fread(file.path(mdir,"PC_NC_nonciliated_epithelial_pseudobulk_metadata_v2.tsv"),data.table=FALSE)
pc_x <- fread(file.path(root,"interim","PC_NC_nonciliated_epithelial_pseudobulk_counts_v2.tsv.gz"),data.table=FALSE)
iua_m <- fread(file.path(rev,"metadata","GSE215968_disjoint_gate_epithelium_label_pseudobulk_metadata.tsv"),data.table=FALSE)
iua_x <- fread(file.path(rev,"interim","GSE215968_disjoint_gate_epithelium_label_pseudobulk_counts.tsv.gz"),data.table=FALSE)
mk_obj <- function(m,x,cond_levels,label) {
 ids <- m$pseudobulk_id; stopifnot(!anyDuplicated(ids),all(ids%in%names(x)),!anyDuplicated(x[[1]]))
 counts <- as.matrix(x[,ids,drop=FALSE]);rownames(counts)<-x[[1]]
 if(any(counts<0)||any(abs(counts-round(counts))>1e-8))stop(paste(label,"is not integer pseudobulk"));storage.mode(counts)<-"integer"
 m$condition<-factor(m$condition %||% m$CycleType,levels=cond_levels); if(any(is.na(m$condition)))stop(paste(label,"condition coding failed")); if(any(table(m$condition)<3))stop(paste(label,"has <3 donors in a group"))
 list(meta=m,counts=counts)
}
`%||%` <- function(x,y) if(is.null(x)) y else x
pc <- mk_obj(transform(pc_m,condition=CycleType),pc_x,c("Natural","Controlled"),"PC_NC")
iua <- mk_obj(iua_m,iua_x,c("WOI_control_local","AS"),"GSE215968_label_epithelium")
if(!all(iua$meta$locked_cell_label=="Epithelium")) stop("IUA comparator is not the declared Epithelium label branch.")
normalised_logcpm <- function(obj,method) {
 design<-model.matrix(~condition,obj$meta); y<-DGEList(counts=obj$counts); keep<-filterByExpr(y,design)
 yn<-normLibSizes(y,method=method); list(logcpm=cpm(yn,log=TRUE,prior.count=2,normalized.lib.sizes=TRUE),keep=rownames(y)[keep],method=method)
}
pc_n <- normalised_logcpm(pc,"TMMwsp"); iua_n <- normalised_logcpm(iua,"TMM")
universe <- sort(intersect(intersect(rownames(pc$counts),rownames(iua$counts)),intersect(pc_n$keep,iua_n$keep)))
if(length(universe)<1000L)stop("Predeclared common expression universe is unexpectedly small.")
frozen <- fread(file.path(project,"results","reference_WOI_program_epithelial.tsv"),data.table=FALSE)$gene
coretab <- fread(file.path(project,"results","reference_epithelial_WOI_program_stability_core_cell_specificity.tsv"),data.table=FALSE); core<-coretab$gene[as.logical(coretab$epithelial_enriched_2fold_reference_rule)]
excluded_sensitivity <- unique(c("ESR1","PGR",frozen,core))
u <- data.frame(gene=universe,in_PC_filterByExpr=universe%in%pc_n$keep,in_IUA_filterByExpr=universe%in%iua_n$keep,in_frozen_programme=universe%in%frozen,in_stable_core=universe%in%core,excluded_in_predeclared_sensitivity=universe%in%excluded_sensitivity)
fwrite(u,file.path(out,"cross_context_gene_universe_v2.tsv"),sep="\t")
effect <- function(logcpm,condition,positive,negative,genes) rowMeans(logcpm[genes,condition==positive,drop=FALSE])-rowMeans(logcpm[genes,condition==negative,drop=FALSE])
e_pc <- effect(pc_n$logcpm,pc$meta$condition,"Controlled","Natural",universe)
e_iua <- effect(iua_n$logcpm,iua$meta$condition,"AS","WOI_control_local",universe)
vectors <- data.frame(gene=universe,PC_minus_NC_logCPM_mean_difference=as.numeric(e_pc),AS_HRT_minus_local_control_logCPM_mean_difference=as.numeric(e_iua),same_direction=sign(e_pc)==sign(e_iua),stringsAsFactors=FALSE)
fwrite(vectors,file.path(out,"cross_context_effect_vectors_v2.tsv"),sep="\t")
boot_metric <- function(genes,B=2000L,seed=202609195L) {
 set.seed(seed); P<-pc_n$logcpm[genes,,drop=FALSE]; I<-iua_n$logcpm[genes,,drop=FALSE]
 pc_pos<-which(pc$meta$condition=="Controlled");pc_neg<-which(pc$meta$condition=="Natural");i_pos<-which(iua$meta$condition=="AS");i_neg<-which(iua$meta$condition=="WOI_control_local")
 vals<-matrix(NA_real_,nrow=B,ncol=2,dimnames=list(NULL,c("spearman","sign_concordance")));fail<-0L
 for(b in seq_len(B)){
  wp<-tabulate(sample(pc_pos,length(pc_pos),TRUE),nbins=ncol(P))/length(pc_pos)-tabulate(sample(pc_neg,length(pc_neg),TRUE),nbins=ncol(P))/length(pc_neg)
  wi<-tabulate(sample(i_pos,length(i_pos),TRUE),nbins=ncol(I))/length(i_pos)-tabulate(sample(i_neg,length(i_neg),TRUE),nbins=ncol(I))/length(i_neg)
  a<-as.numeric(P%*%wp);d<-as.numeric(I%*%wi)
  if(any(!is.finite(a))||any(!is.finite(d))){fail<-fail+1L;next}
  vals[b,"spearman"]<-suppressWarnings(cor(a,d,method="spearman")); nz<-a!=0&d!=0; vals[b,"sign_concordance"]<-if(any(nz))mean(sign(a[nz])==sign(d[nz])) else NA_real_
 }
 list(estimate_spearman=suppressWarnings(cor(effect(pc_n$logcpm,pc$meta$condition,"Controlled","Natural",genes),effect(iua_n$logcpm,iua$meta$condition,"AS","WOI_control_local",genes),method="spearman")),estimate_sign=mean(sign(effect(pc_n$logcpm,pc$meta$condition,"Controlled","Natural",genes))==sign(effect(iua_n$logcpm,iua$meta$condition,"AS","WOI_control_local",genes))),ci_spearman=quantile(vals[,"spearman"],c(.025,.975),na.rm=TRUE,names=FALSE),ci_sign=quantile(vals[,"sign_concordance"],c(.025,.975),na.rm=TRUE,names=FALSE),valid=sum(is.finite(vals[,"spearman"])),failed=fail)
}
main_b <- boot_metric(universe); sens_genes <- setdiff(universe,excluded_sensitivity); sens_b <- boot_metric(sens_genes,seed=202609196L)
con <- data.frame(analysis=c("primary_full_common_expression_universe","predeclared_exclude_receptors_and_programme_core"),genes=c(length(universe),length(sens_genes)),spearman=c(main_b$estimate_spearman,sens_b$estimate_spearman),spearman_bootstrap_ci95_low=c(main_b$ci_spearman[1],sens_b$ci_spearman[1]),spearman_bootstrap_ci95_high=c(main_b$ci_spearman[2],sens_b$ci_spearman[2]),sign_concordance=c(main_b$estimate_sign,sens_b$estimate_sign),sign_concordance_bootstrap_ci95_low=c(main_b$ci_sign[1],sens_b$ci_sign[1]),sign_concordance_bootstrap_ci95_high=c(main_b$ci_sign[2],sens_b$ci_sign[2]),valid_bootstrap=c(main_b$valid,sens_b$valid),failed_bootstrap=c(main_b$failed,sens_b$failed),stringsAsFactors=FALSE)
fwrite(con,file.path(out,"cross_context_concordance_v2.tsv"),sep="\t")

pc_primary <- fread(file.path(out,"PC_NC_primary_endpoint_tests_v2.tsv"),data.table=FALSE);pc_secondary<-fread(file.path(out,"PC_NC_secondary_endpoint_tests_v2.tsv"),data.table=FALSE)
iua_r<-fread(file.path(rev,"results","GSE215968_disjoint_gate_epithelium_label_receptor_effects.tsv"),data.table=FALSE);iua_p<-fread(file.path(rev,"results","GSE215968_disjoint_gate_epithelium_label_programme_tests.tsv"),data.table=FALSE)
ep <- rbind(
 data.frame(context="PC_NC_v2_primary",endpoint=pc_primary$endpoint,direction="Controlled_PC_minus_Natural_NC",effect=pc_primary$difference_PC_minus_NC,ci_low=pc_primary$bootstrap_ci95_low,ci_high=pc_primary$bootstrap_ci95_high,exact_p=pc_primary$exact_permutation_p_two_sided,exact_q=pc_primary$exact_BH_q_within_family,donors="7 Controlled; 9 Natural",compartment="nonciliated epithelial",platform="snRNA-seq",normalization="TMMwsp",inferential_role="healthy hormonal-context reference"),
 data.frame(context="PC_NC_v2_secondary",endpoint=pc_secondary$endpoint,direction="Controlled_PC_minus_Natural_NC",effect=pc_secondary$difference_PC_minus_NC,ci_low=pc_secondary$bootstrap_ci95_low,ci_high=pc_secondary$bootstrap_ci95_high,exact_p=pc_secondary$exact_permutation_p_two_sided,exact_q=pc_secondary$exact_BH_q_within_family,donors="7 Controlled; 9 Natural",compartment="nonciliated epithelial",platform="snRNA-seq",normalization="TMMwsp",inferential_role="secondary healthy-context reference"),
 data.frame(context="GSE215968_disjoint_label_branch",endpoint=iua_r$endpoint,direction="AS_HRT_minus_local_control",effect=iua_r$mean_difference_AS_minus_control,ci_low=iua_r$bootstrap_ci95_low,ci_high=iua_r$bootstrap_ci95_high,exact_p=iua_r$exact_permutation_p_two_sided,exact_q=iua_r$exact_permutation_BH_q_two_receptors,donors="8 AS/HRT; 6 local controls",compartment="Epithelium label, ciliated excluded",platform="scRNA-seq",normalization="TMM",inferential_role="composite disease-and-regimen contrast; not disease independent of HRT"),
 data.frame(context="GSE215968_disjoint_label_branch",endpoint=iua_p$endpoint,direction="AS_HRT_minus_local_control",effect=iua_p$mean_difference_AS_minus_control,ci_low=iua_p$bootstrap_ci95_low,ci_high=iua_p$bootstrap_ci95_high,exact_p=iua_p$exact_permutation_p_two_sided,exact_q=iua_p$exact_permutation_BH_q_two_scores,donors="8 AS/HRT; 6 local controls",compartment="Epithelium label, ciliated excluded",platform="scRNA-seq",normalization="TMM",inferential_role="composite disease-and-regimen contrast; not disease independent of HRT")
)
fwrite(ep,file.path(out,"cross_context_endpoint_comparison_v2.tsv"),sep="\t")
comp <- data.frame(domain=c("independent_donors","disease_vs_hormone_identifiability","cell_compartment","assay","normalization","WOI_timing","clinical_function","allowed_interpretation","forbidden_interpretation"),PC_NC=c("16 healthy/reference donors; 7 Controlled and 9 Natural","Controlled versus Natural context; not a randomised receptor perturbation","Publisher glandular+luminal labels; ciliated excluded","snRNA-seq","TMMwsp","not recoverable from H5AD donor metadata","No implantation outcome","Reference direction for a hormonal-context comparison","Causal receptor activity or clinical receptivity"),GSE215968=c("14 local donors; 8 AS/HRT and 6 local controls","AS is perfectly confounded with HRT/regimen in available contrast","Epithelium label only; ciliated label excluded","scRNA-seq","TMM historical; fixed-normalisation vector comparison","local control is WOI-labelled; full protocol comparability unavailable","No direct outcome validation","Composite AS/HRT contrast and contextual direction only","Independent disease effect, replication of a PC/NC result, or HRT adjustment"),comparability_decision=c("different cohorts; no donor overlap evidenced","not identifiable as disease-only effect","approximate nonciliated correspondence","platform difference remains","within-study normalization only","partial/uncertain","not comparable","descriptive direction comparison permitted","prohibited"))
fwrite(comp,file.path(out,"cross_context_comparability_v2.tsv"),sep="\t")
qa <- list(analysis_id="IUA_CONTEXT_REPAIR_V2",status="complete_conditional_branch",conditional_rule="Both raw count matrices integer; gene symbols and nonciliated label correspondence audited; no disease-only interpretation.",PC_normalization="TMMwsp",IUA_normalization="TMM",PC_donors=16,IUA_donors=14,common_expression_genes=length(universe),bootstrap_replicates=2000,seed=202609195,bootstrap_mode="fixed within-study normalized donor logCPM; donor resampling within condition",primary_branch="GSE215968 Epithelium label, ciliated excluded")
write_json(qa,file.path(mdir,"cross_context_bootstrap_qa_v2.json"),auto_unbox=TRUE,pretty=TRUE)
cat(toJSON(list(status="E_complete_conditional",common_genes=length(universe)),auto_unbox=TRUE),"\n")

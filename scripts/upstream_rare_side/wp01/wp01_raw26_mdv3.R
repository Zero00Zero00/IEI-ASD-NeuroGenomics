#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(brglm2);library(yaml);library(parallel);library(mgcv)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) stop("usage: wp01_raw26_mdv3.R CH3_ROOT")
ROOT <- normalizePath(args[1],mustWork=TRUE)
MA <- file.path(ROOT,"12_molecular_autism_revision")
OUT <- file.path(MA,"01_anchor_ingest","WP01_raw26_primary_v1")
source(file.path(MA,"workflow","scripts","wp01_mdv3_legacy_bridge.R"))
wp01_require_solver_locked_release(ROOT,OUT)

replay_pass <- file.path(OUT,"WP01_SOLVER_FULL_REPLAY_PASS.txt")
if (!file.exists(replay_pass)) stop("Full legacy solver replay PASS required before raw26 MDV3")
rp <- wp01_kv(replay_pass)
if (!identical(unname(rp["status"]),"PASS") || !identical(unname(rp["replay_pass_n"]),"6671"))
  stop("Full legacy solver replay is not certified")
if (!identical(unname(rp["model_frame_sha256"]),wp01_sha256(file.path(OUT,"WP01_model_frame.tsv.gz"))))
  stop("Model frame changed after solver replay")

legacy <- wp01_load_legacy_solver(ROOT,OUT)
if (!identical(unname(rp["legacy_config_sha256"]),wp01_sha256(legacy$cfg_path)) ||
    !identical(unname(rp["legacy_solver_sha256"]),wp01_sha256(legacy$solver_path)))
  stop("Legacy solver/config changed after full replay")

design <- wp01_prepare_mdv3_design(ROOT,OUT)

# Structural audit for the 26-gene primary anchor.
smd <- function(x1,x0) {
  n1<-length(x1); n0<-length(x0); s1<-var(x1); s0<-var(x0)
  sp<-sqrt(((n1-1)*s1+(n0-1)*s0)/(n1+n0-2))
  (mean(x1)-mean(x0))/sp
}
base <- design$base

# Explicit row indices: avoid data.table single-symbol logical-i NSE.
idx_raw26 <- which(base$raw26_flag %in% TRUE)
idx_nonraw26 <- which(base$raw26_flag %in% FALSE)

if (length(idx_raw26) != 26L)
    stop("raw26 explicit index count != 26")

if (length(idx_nonraw26) != 19241L)
    stop("nonraw26 explicit index count != 19241")

raw26_dt <- base[idx_raw26]
nonraw26_dt <- base[idx_nonraw26]

if (nrow(raw26_dt) != 26L ||
    nrow(nonraw26_dt) != 19241L)
    stop("explicit raw26/nonraw26 subset failure")

audit <- data.table::rbindlist(list(
  data.table::data.table(covariate="log_gene_length",
    n_case=sum(base$raw26_flag),n_control=sum(!base$raw26_flag),
    case_mean=mean(raw26_dt$log_gene_length),control_mean=mean(nonraw26_dt$log_gene_length),
    SMD=smd(raw26_dt$log_gene_length,nonraw26_dt$log_gene_length),
    wilcoxon_P=stats::wilcox.test(raw26_dt$log_gene_length,nonraw26_dt$log_gene_length,exact=FALSE)$p.value),
  data.table::data.table(covariate="GC",
    n_case=sum(base$raw26_flag),n_control=sum(!base$raw26_flag),
    case_mean=mean(raw26_dt$GC),control_mean=mean(nonraw26_dt$GC),
    SMD=smd(raw26_dt$GC,nonraw26_dt$GC),
    wilcoxon_P=stats::wilcox.test(raw26_dt$GC,nonraw26_dt$GC,exact=FALSE)$p.value),
  data.table::data.table(covariate="log_annotation",
    n_case=sum(base$raw26_flag),n_control=sum(!base$raw26_flag),
    case_mean=mean(raw26_dt$log_annotation),control_mean=mean(nonraw26_dt$log_annotation),
    SMD=smd(raw26_dt$log_annotation,nonraw26_dt$log_annotation),
    wilcoxon_P=stats::wilcox.test(raw26_dt$log_annotation,nonraw26_dt$log_annotation,exact=FALSE)$p.value)
))
data.table::fwrite(audit,file.path(OUT,"WP01_raw26_bias_audit.tsv"),sep="\t")

gamfit <- mgcv::gam(as.integer(raw26_flag) ~ s(log_gene_length),data=base,family=stats::binomial())
gs <- summary(gamfit)
data.table::fwrite(data.table::data.table(edf=gs$s.table[1,"edf"],P=gs$s.table[1,ncol(gs$s.table)],
  model_basis="canonical_model_frame; exact legacy MDV3 solver/config"),file.path(OUT,"WP01_raw26_length_GAM.tsv"),sep="\t")

res <- wp01_fit_all(design,design$base$raw26_flag,legacy$cfg,threads=8L,label="raw26 primary MDV3")
failed <- res[fit_status!="PASS"|converged!=TRUE]
if (nrow(failed)) {
  data.table::fwrite(failed,file.path(OUT,"WP01_raw26_mdv3_noncertified.tsv"),sep="\t")
  stop("raw26 MDV3 numerical certification failed for ",nrow(failed)," pathway(s)")
}
res[,q_Firth:=stats::p.adjust(P_Firth,method=legacy$cfg$primary_model$fdr_method),by=source]
alpha<-as.numeric(legacy$cfg$primary_model$alpha)
res[,enriched_flag:=fit_status=="PASS"&converged==TRUE&is.finite(OR_Firth)&OR_Firth>1&q_Firth<alpha]
res[,depleted_flag:=fit_status=="PASS"&converged==TRUE&is.finite(OR_Firth)&OR_Firth<1&q_Firth<alpha]
data.table::setorder(res,source,q_Firth,P_Firth,pathway_id)

# Provide current WP01 naming while preserving exact solver metadata.
out <- res[,.(source,pathway_id,pathway_name,
  pathway_size=global_pathway_size,
  raw26_n=outcome_n,raw26_pathway_n=outcome_pathway_n,nonraw26_pathway_n=nonoutcome_pathway_n,
  beta_Firth=beta_pathway,SE_Firth=SE,logCI95L,logCI95U,OR_Firth,CI95L,CI95U,P_Firth,
  q_Firth,fit_status,converged,solver_attempt,solver_type,iterations,prob01_warning,warning,
  enriched_flag,depleted_flag)]
con <- gzfile(file.path(OUT,"WP01_raw26_mdv3_all.tsv.gz"),"wt")
utils::write.table(out,con,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE,na="NA")
close(con)
data.table::fwrite(out[enriched_flag==TRUE],file.path(OUT,"WP01_raw26_mdv3_hits.tsv"),sep="\t")

qc <- data.table::data.table(
  metric=c("n_pathways","n_fit_fail","n_firth_positive","GO_n","Reactome_n",
           "solver_replay_n","solver_replay_pass_n","solver_attempts_used","solver_types_used",
           "prob01_warning_n","model_frame_sha256","legacy_config_sha256","legacy_solver_sha256"),
  value=c(nrow(out),sum(out$fit_status!="PASS"),sum(out$enriched_flag),
          out[source=="GO_BP",.N],out[source=="Reactome",.N],
          6671,6671,paste(sort(unique(out$solver_attempt)),collapse=","),
          paste(sort(unique(out$solver_type)),collapse=","),
          sum(out$prob01_warning),wp01_sha256(design$model_file),
          wp01_sha256(legacy$cfg_path),wp01_sha256(legacy$solver_path))
)
data.table::fwrite(qc,file.path(OUT,"WP01_raw26_mdv3_qc.tsv"),sep="\t")

pass <- c(
  "WP01_RAW26_MDV3",
  "status=PASS",
  "pathways_n=6671",
  "fit_fail_n=0",
  paste0("results_sha256=",wp01_sha256(file.path(OUT,"WP01_raw26_mdv3_all.tsv.gz"))),
  paste0("replay_pass_sha256=",wp01_sha256(replay_pass))
)
writeLines(pass,file.path(OUT,"WP01_RAW26_MDV3_PASS.txt"))
cat("[PASS] raw26 primary MDV3 completed with exact legacy solver and 0 noncertified pathways.\n")

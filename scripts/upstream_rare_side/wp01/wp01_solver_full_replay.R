#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(brglm2);library(yaml);library(parallel)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) stop("usage: wp01_solver_full_replay.R CH3_ROOT")
ROOT <- normalizePath(args[1],mustWork=TRUE)
MA <- file.path(ROOT,"12_molecular_autism_revision")
OUT <- file.path(MA,"01_anchor_ingest","WP01_raw26_primary_v1")
source(file.path(MA,"workflow","scripts","wp01_mdv3_legacy_bridge.R"))
wp01_require_solver_locked_release(ROOT,OUT)
legacy <- wp01_load_legacy_solver(ROOT,OUT)
design <- wp01_prepare_mdv3_design(ROOT,OUT)

# Verify the original frozen pilot/QC facts before replay.
mqc <- data.table::fread(file.path(ROOT,"04_bias_empirical","MDV3_model_qc.tsv"))
if (sum(mqc$N)!=6671L || any(mqc$fit_status!="PASS") || any(mqc$converged!=TRUE) ||
    any(mqc$solver_attempt!="A_ASMEAN_RA05") || any(mqc$solver_type!="AS_mean"))
  stop("Frozen MDV3 model QC no longer matches the certified legacy solver state")
pilot <- data.table::fread(cmd=paste("gzip -dc",shQuote(file.path(ROOT,"04_bias_empirical","MDV3_numerical_pilot.tsv.gz"))))
if (nrow(pilot)!=279L || sum(as.logical(pilot$certified))!=279L ||
    sum(pilot$CoreSeed_pathway_n==0L)!=198L ||
    sum(as.logical(pilot$certified[pilot$CoreSeed_pathway_n==0L]))!=198L)
  stop("Frozen MDV3 numerical pilot facts do not match expected certified state")

# Replay ALL 6,671 pathways under CoreSeed25 using the exact legacy solver/config.
rep <- wp01_fit_all(design,design$base$CoreSeed_flag,legacy$cfg,threads=8L,label="CoreSeed25 full replay")
frozen <- data.table::fread(cmd=paste("gzip -dc",shQuote(file.path(ROOT,"04_bias_empirical","04_pathway_all.tsv.gz"))))
keep <- c("source","pathway_id","beta_pathway","SE","P_Firth","q_Firth","OR_Firth",
          "solver_attempt","solver_type","converged","prob01_warning","enriched_flag","depleted_flag")
if (length(setdiff(keep,names(frozen)))) stop("Frozen legacy MDV3 authority columns missing")
fr <- frozen[,..keep]
data.table::setnames(fr,setdiff(names(fr),c("source","pathway_id")),
                     paste0("legacy_",setdiff(names(fr),c("source","pathway_id"))))
z <- merge(rep,fr,by=c("source","pathway_id"),all=TRUE,sort=FALSE)
if (nrow(z)!=6671L || anyNA(z$legacy_beta_pathway) || anyNA(z$beta_pathway))
  stop("Full replay/frozen authority identity failure")

z[,replay_q_Firth:=stats::p.adjust(P_Firth,method="BH"),by=source]
z[,replay_enriched_flag:=fit_status=="PASS"&converged==TRUE&is.finite(OR_Firth)&OR_Firth>1&replay_q_Firth<0.05]
z[,replay_depleted_flag:=fit_status=="PASS"&converged==TRUE&is.finite(OR_Firth)&OR_Firth<1&replay_q_Firth<0.05]
z[,`:=`(
  abs_beta_diff=abs(beta_pathway-legacy_beta_pathway),
  abs_SE_diff=abs(SE-legacy_SE),
  abs_log10p_diff=abs(log10(pmax(P_Firth,1e-300))-log10(pmax(legacy_P_Firth,1e-300))),
  abs_q_diff=abs(replay_q_Firth-legacy_q_Firth),
  solver_attempt_match=solver_attempt==legacy_solver_attempt,
  solver_type_match=solver_type==legacy_solver_type,
  converged_match=as.logical(converged)==as.logical(legacy_converged),
  prob01_match=as.logical(prob01_warning)==as.logical(legacy_prob01_warning),
  enriched_match=as.logical(replay_enriched_flag)==as.logical(legacy_enriched_flag),
  depleted_match=as.logical(replay_depleted_flag)==as.logical(legacy_depleted_flag)
)]
z[,replay_PASS :=
    fit_status=="PASS" & converged==TRUE &
    is.finite(abs_beta_diff) & abs_beta_diff<=5e-5 &
    is.finite(abs_SE_diff) & abs_SE_diff<=5e-5 &
    is.finite(abs_log10p_diff) & abs_log10p_diff<=1e-3 &
    is.finite(abs_q_diff) & abs_q_diff<=1e-3 &
    solver_attempt_match & solver_type_match & converged_match & prob01_match &
    enriched_match & depleted_match]

data.table::setorder(z,source,pathway_id)
con <- gzfile(file.path(OUT,"WP01_solver_full_replay.tsv.gz"),"wt")
utils::write.table(z,con,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE,na="NA")
close(con)

summary <- data.table::data.table(
  replay_n=nrow(z),
  replay_pass_n=sum(z$replay_PASS),
  replay_fail_n=sum(!z$replay_PASS),
  max_abs_beta_diff=max(z$abs_beta_diff,na.rm=TRUE),
  max_abs_SE_diff=max(z$abs_SE_diff,na.rm=TRUE),
  max_abs_log10p_diff=max(z$abs_log10p_diff,na.rm=TRUE),
  max_abs_q_diff=max(z$abs_q_diff,na.rm=TRUE),
  solver_attempt_mismatch_n=sum(!z$solver_attempt_match),
  solver_type_mismatch_n=sum(!z$solver_type_match),
  enriched_mismatch_n=sum(!z$enriched_match),
  depleted_mismatch_n=sum(!z$depleted_match),
  noncertified_n=sum(z$fit_status!="PASS"|z$converged!=TRUE),
  prob01_warning_n=sum(as.logical(z$prob01_warning))
)
data.table::fwrite(summary,file.path(OUT,"WP01_solver_full_replay_summary.tsv"),sep="\t")

# Backward-compatible 30-pathway audit table, now derived from the full replay.
cand <- frozen[is.finite(P_Firth)&is.finite(beta_pathway)][order(P_Firth,source,pathway_id)]
top <- head(cand,10)
rest <- cand[order(source,pathway_id)]
if (nrow(rest)>20) rest <- rest[unique(round(seq(1,nrow(rest),length.out=20)))]
certset <- unique(rbind(top,rest),by=c("source","pathway_id"))
cert <- z[certset,on=.(source,pathway_id)]
cert[,cert_PASS:=replay_PASS]
data.table::fwrite(cert,file.path(OUT,"WP01_solver_certification.tsv"),sep="\t")

if (!all(z$replay_PASS)) {
  bad <- z[replay_PASS!=TRUE]
  data.table::fwrite(bad,file.path(OUT,"WP01_solver_full_replay_failures.tsv"),sep="\t")
  stop("Legacy full replay certification failed: ",sum(z$replay_PASS),"/6671 pass")
}

# Freeze the replay checkpoint with all algorithmic authorities.
pass <- c(
  "WP01_SOLVER_FULL_REPLAY",
  "status=PASS",
  "replay_n=6671",
  "replay_pass_n=6671",
  paste0("preanalysis_lock_sha256=",wp01_sha256(file.path(OUT,"WP01_preanalysis_lock.json"))),
  paste0("model_frame_sha256=",wp01_sha256(file.path(OUT,"WP01_model_frame.tsv.gz"))),
  paste0("legacy_config_sha256=",wp01_sha256(legacy$cfg_path)),
  paste0("legacy_solver_sha256=",wp01_sha256(legacy$solver_path)),
  paste0("legacy_pathway_all_sha256=",wp01_sha256(file.path(ROOT,"04_bias_empirical","04_pathway_all.tsv.gz"))),
  paste0("replay_results_sha256=",wp01_sha256(file.path(OUT,"WP01_solver_full_replay.tsv.gz"))),
  paste0("summary_sha256=",wp01_sha256(file.path(OUT,"WP01_solver_full_replay_summary.tsv")))
)
writeLines(pass,file.path(OUT,"WP01_SOLVER_FULL_REPLAY_PASS.txt"))
cat("[PASS] Legacy CoreSeed25 full replay certified 6,671/6,671 pathways.\n")

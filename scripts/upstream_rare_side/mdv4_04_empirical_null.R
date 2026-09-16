#!/usr/bin/env Rscript
source("workflow/scripts/mdv4_common.R")
a <- mdv4_parse_args(); if (is.null(a$config)) mdv4_stop("--config required")
cfg <- mdv4_read_cfg(a$config); setwd(cfg$project_root); outdir <- cfg$output_dir; mdv4_dir(outdir)
ns <- mdv4_kv(file.path(outdir,"MDV4_NULLSETS_COMPLETE.txt")); if (!identical(unname(ns["status"]),"COMPLETE")) mdv4_stop("Null sets not complete")
u <- mdv4_prepare_universe(fread(cfg$upstream$universe),cfg)
mem <- fread(cfg$upstream$membership)
paths <- fread(cfg$upstream$mdv3_pathways)
assign <- fread(file.path(outdir,"MDV4_pseudo_assignments.tsv.gz"))
if (nrow(paths)!=cfg$expected$total_pathways_n) mdv4_stop("Pathway family mismatch")
M <- mdv4_make_membership_matrix(u,mem,paths,cfg)
# Global sizes and observed CoreSeed counts must agree exactly with frozen MDV-3 results.
m <- as.integer(Matrix::colSums(M)); core_idx <- which(u[[cfg$columns$coreseed]]==1L); aobs <- as.integer(Matrix::colSums(M[core_idx,,drop=FALSE]))
if (!all(m==as.integer(paths$global_pathway_size))) mdv4_stop("Pathway size mismatch vs MDV3")
if (!all(aobs==as.integer(paths$CoreSeed_pathway_n))) mdv4_stop("CoreSeed pathway count mismatch vs MDV3")
res <- mdv4_empirical_from_assignments(assign,u,M,paths,cfg,B_expected=cfg$matching$primary_B)
res <- mdv4_add_qemp(res)
# Join the frozen adjusted MDV-3 evidence without allowing it to alter empirical computation.
md3 <- paths[,.(source,pathway_id,OR_Firth,q_Firth,enriched_flag,depleted_flag,fit_status,converged,mdv2_short_id,mdv2_primary_reporting_flag)]
res <- merge(res,md3,by=c("source","pathway_id"),all.x=TRUE,sort=FALSE)
res[, .ord__ := match(paste(source,pathway_id,sep="||"), paste(paths$source,paths$pathway_id,sep="||"))]; setorder(res,.ord__); res[, .ord__ := NULL]
mdv4_atomic_fwrite(res,file.path(outdir,"MDV4_empirical_all.tsv.gz"),compress="gzip")
mdv4_atomic_lines(c(
  "MDV4_EMPIRICAL_COMPLETE",
  "status=COMPLETE",
  paste0("version=",cfg$version),
  paste0("timestamp=",mdv4_now()),
  paste0("pathways_n=",nrow(res)),
  paste0("GO_BP_n=",res[source=="GO_BP",.N]),
  paste0("Reactome_n=",res[source=="Reactome",.N]),
  paste0("B=",cfg$matching$primary_B),
  paste0("min_possible_P_emp=",format(1/(cfg$matching$primary_B+1),scientific=TRUE,digits=8)),
  paste0("q_emp_lt_0.05_n=",res[q_emp<cfg$primary_rule$alpha,.N]),
  "observed_and_permuted_statistic=IDENTICAL_CONDITIONAL_HALDANE_LOG_OR",
  "BH_family=COMPLETE_FROZEN_GO_BP_AND_REACTOME_SEPARATELY",
  "next=MDV4_TIER_ASSIGNMENT"
),file.path(outdir,"MDV4_EMPIRICAL_COMPLETE.txt"))
message("[PASS] MDV-4 10,000-set empirical null completed.")

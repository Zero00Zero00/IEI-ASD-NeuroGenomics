#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(brglm2);library(yaml);library(parallel)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<4) stop("usage: mdv3_pseudo.R CH3_ROOT GROUPS TMP OUT")
ROOT<-normalizePath(args[1],mustWork=TRUE); GROUPS<-args[2]; TMP<-args[3]; OUTF<-args[4]
MA<-file.path(ROOT,"12_molecular_autism_revision"); FROZEN<-file.path(MA,"01_anchor_ingest","WP01_raw26_primary_v1")
dir.create(TMP,recursive=TRUE,showWarnings=FALSE)
file.copy(file.path(FROZEN,"WP01_model_frame.tsv.gz"),file.path(TMP,"WP01_model_frame.tsv.gz"),overwrite=TRUE)
source(file.path(MA,"workflow","scripts","wp01_mdv3_legacy_bridge.R"))
legacy<-wp01_load_legacy_solver(ROOT,TMP); design<-wp01_prepare_mdv3_design(ROOT,TMP)
g<-fread(GROUPS); overlap<-g[R0_group=="Overlap_raw",gene]
outcome<-design$base$gene %in% overlap
res<-wp01_fit_all(design,outcome,legacy$cfg,threads=as.integer(Sys.getenv("MAR07_PATHWAY_THREADS","8")),label="MAR07 null MDV3")
res[,q_Firth:=p.adjust(P_Firth,method=legacy$cfg$primary_model$fdr_method),by=source]
alpha<-as.numeric(legacy$cfg$primary_model$alpha)
res[,enriched_flag:=fit_status=="PASS"&converged==TRUE&is.finite(OR_Firth)&OR_Firth>1&q_Firth<alpha]
fwrite(res,OUTF,sep="\t")

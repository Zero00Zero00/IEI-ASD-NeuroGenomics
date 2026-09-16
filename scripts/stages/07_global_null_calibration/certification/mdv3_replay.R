#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(brglm2);library(yaml);library(parallel)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<3) stop("usage: mdv3_replay.R CH3_ROOT OUTDIR OUTFILE")
ROOT<-normalizePath(args[1],mustWork=TRUE); TMP<-args[2]; OUTF<-args[3]
MA<-file.path(ROOT,"12_molecular_autism_revision")
FROZEN<-file.path(MA,"01_anchor_ingest","WP01_raw26_primary_v1")
dir.create(TMP,recursive=TRUE,showWarnings=FALSE)
file.copy(file.path(FROZEN,"WP01_model_frame.tsv.gz"),file.path(TMP,"WP01_model_frame.tsv.gz"),overwrite=TRUE)
source(file.path(MA,"workflow","scripts","wp01_mdv3_legacy_bridge.R"))
legacy<-wp01_load_legacy_solver(ROOT,TMP)
design<-wp01_prepare_mdv3_design(ROOT,TMP)
res<-wp01_fit_all(design,design$base$raw26_flag,legacy$cfg,
                  threads=as.integer(Sys.getenv("MAR07_PATHWAY_THREADS","8")),label="MAR07 MDV3 replay")
res[,q_Firth:=p.adjust(P_Firth,method=legacy$cfg$primary_model$fdr_method),by=source]
alpha<-as.numeric(legacy$cfg$primary_model$alpha)
res[,enriched_flag:=fit_status=="PASS"&converged==TRUE&is.finite(OR_Firth)&OR_Firth>1&q_Firth<alpha]
fwrite(res,OUTF,sep="\t")
cat("[PASS] MDV3 exact legacy-solver replay completed\n")

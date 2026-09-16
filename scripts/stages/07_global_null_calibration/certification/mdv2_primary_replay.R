#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(yaml);library(data.table);library(brglm2);library(parallel);library(splines)})

args<-commandArgs(trailingOnly=TRUE)
if(length(args)<4) stop("usage: mdv2_primary_replay.R GENE_UNIVERSE PATHWAY_MEMBERSHIP CONFIG OUT")
u_path<-args[1]; pm_path<-args[2]; cfg_path<-args[3]; out_path<-args[4]
cfg<-yaml::read_yaml(cfg_path); cc<-cfg$columns
u<-data.table::fread(u_path,sep="\t",showProgress=FALSE)
pm<-data.table::fread(cmd=paste("gzip -dc",shQuote(pm_path)),showProgress=FALSE)

pmv<-pm[,.(source=get(cc$pathway_source),pathway_id=get(cc$pathway_id),
           pathway_name=get(cc$pathway_name),approved_symbol=get(cc$pathway_symbol))]
paths<-pmv[,.(pathway_name=unique(pathway_name)[1],
              global_pathway_size=data.table::uniqueN(approved_symbol)),by=.(source,pathway_id)]
data.table::setorder(paths,source,pathway_id)
member_sets<-split(pmv$approved_symbol,paste(pmv$source,pmv$pathway_id,sep="\r"))
pkeys<-paste(paths$source,paths$pathway_id,sep="\r")
alpha<-as.numeric(cfg$primary_model$alpha)
zcrit<-qnorm(0.975)

build_contrast<-function(contrast){
  group<-u[[cc$group]]
  keep<-if(contrast=="C1") group%in%c("Overlap_raw","IEI_only") else group%in%c("Overlap_raw","ASD_only")
  y<-if(contrast=="C1") as.integer(group[keep]=="Overlap_raw") else as.integer(group[keep]=="Overlap_raw")
  d<-u[keep]; d[,outcome:=y]; d
}
formula_for<-function(){
  as.formula(sprintf("outcome ~ pathway_member + splines::ns(%s, df=%d) + %s + scale(log1p(%s))",
    cc$log_gene_length,as.integer(cfg$primary_model$spline_df),cc$gc_z,cc$annotation_degree))
}
fit_one<-function(d,ms,f,p){
  d<-copy(d); d[,pathway_member:=as.integer(get(cc$universe_symbol)%in%ms)]
  ncase<-sum(d$outcome==1L); nctrl<-sum(d$outcome==0L)
  pcase<-sum(d$pathway_member[d$outcome==1L]); pctrl<-sum(d$pathway_member[d$outcome==0L])
  if(length(unique(d$pathway_member))<2L)
    return(list(beta=NA_real_,SE=NA_real_,OR=NA_real_,P=1.0,fit_status="NO_VARIATION",
                n_case=ncase,n_control=nctrl,pathway_case=pcase,pathway_control=pctrl))
  fit<-tryCatch(suppressWarnings(glm(f,data=d,family=binomial("logit"),method=brglm2::brglmFit,
                      type=cfg$primary_model$brglm2_type)),error=function(e)e)
  if(inherits(fit,"error"))
    return(list(beta=NA_real_,SE=NA_real_,OR=NA_real_,P=1.0,fit_status="FIT_ERROR",
                n_case=ncase,n_control=nctrl,pathway_case=pcase,pathway_control=pctrl))
  cf<-coef(fit); vc<-vcov(fit)
  beta<-as.numeric(cf[["pathway_member"]]); se<-sqrt(as.numeric(vc["pathway_member","pathway_member"]))
  pp<-if(is.finite(beta)&&is.finite(se)&&se>0) 2*pnorm(abs(beta/se),lower.tail=FALSE) else 1.0
  list(beta=beta,SE=se,OR=exp(beta),P=pp,fit_status="PASS",
       n_case=ncase,n_control=nctrl,pathway_case=pcase,pathway_control=pctrl)
}
d1<-build_contrast("C1"); d2<-build_contrast("C2"); f<-formula_for()
worker<-function(i){
  ms<-member_sets[[pkeys[i]]]; p<-paths[i]
  a<-fit_one(d1,ms,f,p); b<-fit_one(d2,ms,f,p)
  data.table(source=p$source,pathway_id=p$pathway_id,pathway_name=p$pathway_name,
             global_pathway_size=p$global_pathway_size,
             C1_n_case=a$n_case,C1_n_control=a$n_control,C1_pathway_case=a$pathway_case,C1_pathway_control=a$pathway_control,
             beta_C1=a$beta,SE_C1=a$SE,OR_C1=a$OR,P_C1=a$P,fit_status_C1=a$fit_status,
             C2_n_case=b$n_case,C2_n_control=b$n_control,C2_pathway_case=b$pathway_case,C2_pathway_control=b$pathway_control,
             beta_C2=b$beta,SE_C2=b$SE,OR_C2=b$OR,P_C2=b$P,fit_status_C2=b$fit_status)
}
threads<-as.integer(Sys.getenv("MAR07_PATHWAY_THREADS","8"))
rr<-if(.Platform$OS.type=="unix"&&threads>1) mclapply(seq_len(nrow(paths)),worker,mc.cores=threads,mc.preschedule=TRUE) else lapply(seq_len(nrow(paths)),worker)
out<-rbindlist(rr,fill=TRUE)
out[,P_conj:=pmax(P_C1,P_C2)]
out[,q_conj:=p.adjust(P_conj,method=cfg$primary_model$fdr_method),by=source]
out[,direction_pass:=is.finite(OR_C1)&is.finite(OR_C2)&OR_C1>1&OR_C2>1]
out[,conjunction_testable:=fit_status_C1%in%c("PASS","NO_VARIATION") & fit_status_C2%in%c("PASS","NO_VARIATION")]
out[,primary_pass:=conjunction_testable & direction_pass & q_conj<alpha]
data.table::fwrite(out,out_path,sep="\t")
cat("[PASS] MDV2 primary-only replay adapter completed\n")

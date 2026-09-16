# WP01 v1.3.0 shared bridge to the exact frozen MDV3 numerical solver.
# The original solver source and original YAML config remain READ ONLY and are the algorithmic authority.

wp01_sha256 <- function(p) {
  z <- system2("sha256sum", shQuote(p), stdout=TRUE)
  strsplit(z,"[[:space:]]+")[[1]][1]
}
wp01_kv <- function(p) {
  z <- readLines(p,warn=FALSE); z <- z[grepl("=",z,fixed=TRUE)]
  x <- strsplit(z,"=",fixed=TRUE)
  setNames(vapply(x,function(v) paste(v[-1],collapse="="),""),vapply(x,`[`,"",1))
}

wp01_require_solver_locked_release <- function(ROOT, OUT) {
  pre <- file.path(OUT,"WP01_preanalysis_lock.json")
  man <- file.path(OUT,"WP01_SOLVER_LOCKED_IMPLEMENTATION_v1.3.2.json")
  rel <- file.path(OUT,"WP01_SOLVER_LOCKED_IMPLEMENTATION_RELEASE_v1.3.2.txt")
  if (!all(file.exists(c(pre,man,rel)))) stop("WP01 v1.3.0 solver-locked manifest/release missing")
  kv <- wp01_kv(rel)
  if (!identical(unname(kv["decision"]),"PASS") ||
      !identical(unname(kv["implementation"]),"SOLVER_LOCKED_v1.3.2"))
    stop("WP01 v1.3.0 solver-locked release invalid")
  if (!identical(unname(kv["preanalysis_lock_sha256"]),wp01_sha256(pre)))
    stop("WP01 v1.3.0 release bound to a different preanalysis lock")
  if (!identical(unname(kv["implementation_manifest_sha256"]),wp01_sha256(man)))
    stop("WP01 v1.3.0 release bound to a different implementation manifest")
  invisible(TRUE)
}

wp01_load_legacy_solver <- function(ROOT, OUT) {
  for (p in c("yaml","data.table","brglm2","parallel")) {
    if (!requireNamespace(p,quietly=TRUE)) stop("Missing R package: ",p)
  }
  cfg_path <- file.path(ROOT,"config","mdv3_coreseed_bias_firth_gate_v1.yaml")
  solver_path <- file.path(ROOT,"workflow","scripts","mdv3_firth_solver_v1p3.R")
  if (!file.exists(cfg_path) || !file.exists(solver_path)) stop("Legacy MDV3 config/solver missing")

  cfg <- yaml::read_yaml(cfg_path)
  if (!identical(as.character(cfg$version),"mdv3_coreseed_bias_firth_gate_v1.0"))
    stop("Legacy MDV3 config version mismatch: ",cfg$version)
  if (!identical(as.character(cfg$implementation_patch),"numerical_certification_v1p3"))
    stop("Legacy MDV3 implementation patch mismatch: ",cfg$implementation_patch)

  cur_brglm <- as.character(utils::packageVersion("brglm2"))
  if (!identical(cur_brglm,"1.1.0"))
    stop("Exact legacy brglm2 version 1.1.0 required; installed=",cur_brglm)
  if (!startsWith(R.version.string,"R version 4.5.1"))
    stop("Exact legacy R 4.5.1 runtime required; current=",R.version.string)

  source(solver_path, local=.GlobalEnv)
  mdv3_check_brglm2_version(cfg)

  attempts <- cfg$numerical_certification$attempts
  if (!length(attempts)) stop("Legacy numerical-certification attempts are empty")
  ids <- vapply(attempts,function(a) as.character(a$id),character(1))
  types <- vapply(attempts,function(a) as.character(a$type),character(1))
  if (!("A_ASMEAN_RA05" %in% ids)) stop("Legacy A_ASMEAN_RA05 attempt not present in YAML config")
  if (!identical(types[match("A_ASMEAN_RA05",ids)],"AS_mean"))
    stop("A_ASMEAN_RA05 is not AS_mean in the legacy config")

  # Flatten the exact live YAML solver contract for audit; no guessed parameters.
  rows <- lapply(seq_along(attempts),function(i) {
    a <- attempts[[i]]
    data.table::data.table(
      order=i,
      id=as.character(a$id),
      type=as.character(a$type),
      epsilon=as.numeric(a$epsilon),
      maxit=as.integer(a$maxit),
      slowit=as.numeric(a$slowit),
      response_adjustment=as.numeric(a$response_adjustment),
      max_step_factor=as.integer(a$max_step_factor),
      a=as.numeric(a$a)
    )
  })
  contract <- data.table::rbindlist(rows,fill=TRUE)
  data.table::fwrite(contract,file.path(OUT,"WP01_legacy_solver_contract.tsv"),sep="\t")

  env <- data.table::data.table(
    item=c("R.version.string","brglm2","config_sha256","solver_sha256"),
    value=c(R.version.string,cur_brglm,wp01_sha256(cfg_path),wp01_sha256(solver_path))
  )
  data.table::fwrite(env,file.path(OUT,"WP01_solver_environment.tsv"),sep="\t")

  list(cfg=cfg,cfg_path=cfg_path,solver_path=solver_path,contract=contract)
}

wp01_prepare_mdv3_design <- function(ROOT, OUT) {
  model_file <- file.path(OUT,"WP01_model_frame.tsv.gz")
  base <- data.table::fread(cmd=paste("gzip -dc",shQuote(model_file)))
  req <- c("gene","CoreSeed_flag","raw26_flag",
           "length_ns1","length_ns2","length_ns3","GC_z_model","annotation_z_model")
  miss <- setdiff(req,names(base))
  if (length(miss)) stop("Canonical model-frame columns missing: ",paste(miss,collapse=","))
  if (nrow(base)!=19267L || anyDuplicated(base$gene)) stop("Canonical model-frame row/key failure")
  if (sum(as.integer(base$CoreSeed_flag))!=25L) stop("CoreSeed count !=25")
  if (sum(as.integer(base$raw26_flag))!=26L) stop("raw26 count !=26")
  if (anyNA(base[,..req])) stop("Required canonical model fields contain NA")
  base[,`:=`(CoreSeed_flag=as.logical(CoreSeed_flag),raw26_flag=as.logical(raw26_flag))]

  pm <- data.table::fread(cmd=paste("gzip -dc",shQuote(file.path(ROOT,"02_pathways","03_pathway_membership.tsv.gz"))))
  reqp <- c("source","pathway_id","pathway_name","approved_symbol")
  if (length(setdiff(reqp,names(pm)))) stop("Frozen pathway-membership schema mismatch")
  pmv <- unique(pm[,.(source=as.character(source),pathway_id=as.character(pathway_id),
                      pathway_name=as.character(pathway_name),gene=as.character(approved_symbol))])
  paths <- pmv[,.(pathway_name=pathway_name[1],global_pathway_size=data.table::uniqueN(gene)),by=.(source,pathway_id)]
  data.table::setorder(paths,source,pathway_id)
  if (nrow(paths)!=6671L || paths[source=="GO_BP",.N]!=5450L || paths[source=="Reactome",.N]!=1221L)
    stop("Frozen pathway family mismatch")

  member_sets <- split(pmv$gene,paste(pmv$source,pmv$pathway_id,sep="\r"))
  path_keys <- paste(paths$source,paths$pathway_id,sep="\r")
  fml <- stats::as.formula(
    "outcome ~ pathway_member + length_ns1 + length_ns2 + length_ns3 + GC_z_model + annotation_z_model"
  )
  list(base=base,paths=paths,member_sets=member_sets,path_keys=path_keys,fml=fml,model_file=model_file)
}

wp01_fit_all <- function(design, outcome, cfg, threads=8L, label="MDV3") {
  base <- design$base
  paths <- design$paths
  member_sets <- design$member_sets
  path_keys <- design$path_keys
  fml <- design$fml

  base_d <- data.table::data.table(
    outcome=as.integer(outcome),
    length_ns1=base$length_ns1,length_ns2=base$length_ns2,length_ns3=base$length_ns3,
    GC_z_model=base$GC_z_model,annotation_z_model=base$annotation_z_model,
    gene=base$gene
  )
  fit_one <- function(i) {
    p <- paths[i]
    ms <- member_sets[[path_keys[i]]]
    d <- data.table::copy(base_d)
    d[,pathway_member:=as.integer(gene %in% ms)]
    fit <- mdv3_certified_fit(fml,d,cfg,coef_name="pathway_member",ci_level=cfg$primary_model$ci_level)
    data.table::data.table(
      source=p$source,pathway_id=p$pathway_id,pathway_name=p$pathway_name,
      global_pathway_size=as.integer(p$global_pathway_size),
      outcome_n=sum(d$outcome==1L),
      outcome_pathway_n=sum(d$pathway_member[d$outcome==1L]),
      nonoutcome_pathway_n=sum(d$pathway_member[d$outcome==0L]),
      beta_pathway=fit$beta,SE=fit$se,logCI95L=fit$log_ci_low,logCI95U=fit$log_ci_high,
      OR_Firth=fit$or,CI95L=fit$ci_low,CI95U=fit$ci_high,P_Firth=fit$p,
      fit_status=if (fit$certified) "PASS" else "NUMERICAL_NONCONVERGENCE",
      converged=isTRUE(fit$converged),solver_attempt=fit$solver_attempt,solver_type=fit$solver_type,
      iterations=fit$iterations,prob01_warning=isTRUE(fit$prob01_warning),warning=fit$warning
    )
  }
  idx <- seq_len(nrow(paths))
  cat(sprintf("[INFO] %s: fitting %d pathways with exact legacy certified solver...\n",label,length(idx)))
  rr <- if (.Platform$OS.type=="unix" && threads>1L) {
    parallel::mclapply(idx,fit_one,mc.cores=threads,mc.preschedule=TRUE)
  } else lapply(idx,fit_one)
  out <- data.table::rbindlist(rr,fill=TRUE)
  if (nrow(out)!=6671L) stop(label,": result row count !=6671")
  out
}

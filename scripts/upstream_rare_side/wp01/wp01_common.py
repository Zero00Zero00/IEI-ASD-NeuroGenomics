#!/usr/bin/env python3
from __future__ import annotations
import csv, gzip, hashlib, json, math, os, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

PACKAGE_VERSION = "1.0.0"

def sha256_file(path: Path, block=1<<20):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        while True:
            b=f.read(block)
            if not b: break
            h.update(b)
    return h.hexdigest()

def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path,"rt",encoding="utf-8",errors="replace")
    return open(path,"rt",encoding="utf-8",errors="replace")

def truthy(v):
    if v is None: return False
    s=str(v).strip().lower()
    return s in {"1","true","t","yes","y","pass"}

def norm(s):
    return re.sub(r"[^a-z0-9]+","",str(s).lower())

def choose_col(cols, aliases):
    cols=list(cols)
    nmap={norm(c):c for c in cols}
    for a in aliases:
        if a in cols: return a
        na=norm(a)
        if na in nmap: return nmap[na]
    return None

SYMBOL_ALIASES=["approved_symbol","HGNC_symbol","hgnc_symbol","gene","symbol","Gene","SYMBOL"]
IUIS_ALIASES=["IUIS_flag","IUIS","iuis_gene","is_iuis"]
SFARI_ALIASES=["SFARI_R0_highconf_flag","SFARI_R0_high_conf_flag","SFARI_flag","SFARI_HC",
               "SFARI_highconf","SFARI_high_conf_flag","sfari_r0_flag","SFARI"]
CORE_ALIASES=["CoreSeed_flag","CoreSeed","coreseed","historical_coreseed"]
GROUP_ALIASES=["R0_group","group","gene_group","Group","four_group"]
SOURCE_ALIASES=["source","database","db","pathway_source"]
PID_ALIASES=["pathway_id","term_id","id","pathway","set_id"]
PNAME_ALIASES=["pathway_name","term_name","name","description","pathway_label"]

def load_groups(root: Path):
    p=root/"01_gene_universe/02_gene_groups.tsv"
    df=pd.read_csv(p,sep="\t",dtype=str)
    sym=choose_col(df.columns,SYMBOL_ALIASES)
    iuis=choose_col(df.columns,IUIS_ALIASES)
    sfari=choose_col(df.columns,SFARI_ALIASES)
    core=choose_col(df.columns,CORE_ALIASES)
    group=choose_col(df.columns,GROUP_ALIASES)
    missing=[x for x,v in [("symbol",sym),("IUIS",iuis),("SFARI_R0",sfari),("CoreSeed",core),("R0_group",group)] if not v]
    if missing: raise RuntimeError("Unresolved 02_gene_groups columns: "+",".join(missing)+" header="+",".join(df.columns))
    out=pd.DataFrame({
        "gene":df[sym].astype(str).str.strip(),
        "IUIS_flag":df[iuis].map(truthy),
        "SFARI_R0_highconf_flag":df[sfari].map(truthy),
        "CoreSeed_flag":df[core].map(truthy),
        "R0_group":df[group].astype(str).str.strip()
    })
    if out["gene"].duplicated().any():
        raise RuntimeError("Duplicate gene symbols in 02_gene_groups.tsv")
    return out

def _first_numeric(df, aliases):
    c=choose_col(df.columns,aliases)
    if c is None: return None,None
    return c,pd.to_numeric(df[c],errors="coerce")

def load_structural_frame(root: Path):
    groups=load_groups(root)
    up=root/"01_gene_universe/01_gene_universe.tsv"
    u=pd.read_csv(up,sep="\t",low_memory=False)
    usym=choose_col(u.columns,SYMBOL_ALIASES)
    if usym is None: raise RuntimeError("Universe gene symbol column unresolved: "+",".join(u.columns))
    u=u.rename(columns={usym:"gene"})
    u["gene"]=u["gene"].astype(str).str.strip()

    logc,logv=_first_numeric(u,["log_gene_length","ln_gene_length","log_length","log_gene_body_length"])
    if logv is None:
        lc,lv=_first_numeric(u,["gene_length","gene_body_length","length","length_bp"])
        if lv is None: raise RuntimeError("Gene length column unresolved")
        if (lv<=0).any(): raise RuntimeError("Non-positive gene length")
        logv=np.log(lv.astype(float))
        logc=f"log({lc})"

    gcc,gcv=_first_numeric(u,["GC","gc","GC_fraction","gc_fraction","gene_body_gc","gene_body_GC","GC_content"])
    if gcv is None: raise RuntimeError("GC column unresolved")
    gcv=gcv.astype(float)
    if np.nanmedian(gcv)>1.5:
        gcv=gcv/100.0

    annc,annv=_first_numeric(u,["annotation_degree","pathway_annotation_degree","anno_degree"])
    if annv is None:
        ap=root/"02_pathways/annotation_degree.tsv"
        a=pd.read_csv(ap,sep="\t",low_memory=False)
        asym=choose_col(a.columns,SYMBOL_ALIASES)
        adeg=choose_col(a.columns,["annotation_degree","pathway_annotation_degree","anno_degree","degree"])
        if asym is None or adeg is None:
            raise RuntimeError("annotation_degree unresolved in universe and annotation_degree.tsv")
        a=a[[asym,adeg]].rename(columns={asym:"gene",adeg:"annotation_degree"})
        a["gene"]=a["gene"].astype(str).str.strip()
        a["annotation_degree"]=pd.to_numeric(a["annotation_degree"],errors="coerce")
        u=u.merge(a,on="gene",how="left")
        annv=u["annotation_degree"]
        annc="annotation_degree.tsv"
    else:
        u["annotation_degree"]=annv

    frame=pd.DataFrame({
        "gene":u["gene"].astype(str),
        "log_gene_length":logv.astype(float),
        "GC":gcv.astype(float),
        "annotation_degree":pd.to_numeric(u["annotation_degree"],errors="coerce").astype(float),
    })
    for c in ["chr","chromosome","seqname","seqnames"]:
        if c in u.columns:
            frame["chr"]=u[c].astype(str); break
    for c in ["start","gene_start","start_position"]:
        if c in u.columns:
            frame["start"]=pd.to_numeric(u[c],errors="coerce"); break
    for c in ["end","gene_end","end_position"]:
        if c in u.columns:
            frame["end"]=pd.to_numeric(u[c],errors="coerce"); break

    frame=frame.merge(groups,on="gene",how="inner",validate="one_to_one")
    if len(frame)!=len(groups):
        raise RuntimeError(f"Universe/groups merge mismatch: {len(frame)} vs {len(groups)}")
    if frame[["log_gene_length","GC","annotation_degree"]].isna().any().any():
        bad=frame[frame[["log_gene_length","GC","annotation_degree"]].isna().any(axis=1)]
        raise RuntimeError(f"Missing primary structural covariates for {len(bad)} genes")
    frame["log_annotation"]=np.log1p(frame["annotation_degree"].astype(float))
    for c in ["log_gene_length","GC","log_annotation"]:
        sd=frame[c].std(ddof=1)
        if not np.isfinite(sd) or sd<=0: raise RuntimeError(f"Invalid SD for {c}")
        frame["z_"+c]=(frame[c]-frame[c].mean())/sd
    return frame, {"log_length_source":logc,"gc_source":gcc,"annotation_source":annc}

def load_pathway_manifest(root: Path):
    p=root/"02_pathways/pathway_manifest.tsv"
    df=pd.read_csv(p,sep="\t",low_memory=False)
    src=choose_col(df.columns,SOURCE_ALIASES); pid=choose_col(df.columns,PID_ALIASES)
    primary=choose_col(df.columns,["primary_10_500_flag","primary_flag","filtered_flag","retain_flag"])
    pname=choose_col(df.columns,PNAME_ALIASES)
    if not src or not pid or not primary:
        raise RuntimeError(f"Pathway manifest schema unresolved: source={src}, pid={pid}, primary={primary}")
    d=df[df[primary].map(truthy)].copy()
    out=pd.DataFrame({"source":d[src].astype(str).str.strip(),"pathway_id":d[pid].astype(str).str.strip()})
    out["pathway_name"]=d[pname].astype(str) if pname else out["pathway_id"]
    ncol=choose_col(d.columns,["n_genes","size","gene_count","pathway_size"])
    out["pathway_size_manifest"]=pd.to_numeric(d[ncol],errors="coerce") if ncol else np.nan
    return out.drop_duplicates(["source","pathway_id"])

def load_pathway_membership(root: Path):
    p=root/"02_pathways/03_pathway_membership.tsv.gz"
    df=pd.read_csv(p,sep="\t",compression="gzip",low_memory=False)
    src=choose_col(df.columns,SOURCE_ALIASES); pid=choose_col(df.columns,PID_ALIASES)
    gene=choose_col(df.columns,SYMBOL_ALIASES); pname=choose_col(df.columns,PNAME_ALIASES)
    if not src or not pid or not gene:
        raise RuntimeError(f"Pathway membership schema unresolved: source={src}, pid={pid}, gene={gene}")
    out=pd.DataFrame({
        "source":df[src].astype(str).str.strip(),
        "pathway_id":df[pid].astype(str).str.strip(),
        "gene":df[gene].astype(str).str.strip()
    })
    out["pathway_name"]=df[pname].astype(str) if pname else out["pathway_id"]
    out=out.drop_duplicates(["source","pathway_id","gene"])
    return out

def bh(pvals):
    p=np.asarray(pvals,dtype=float)
    n=len(p); out=np.full(n,np.nan)
    ok=np.isfinite(p)
    idx=np.where(ok)[0]
    if len(idx)==0: return out
    q=p[idx]
    order=np.argsort(q)
    ranked=q[order]*len(q)/(np.arange(len(q))+1)
    ranked=np.minimum.accumulate(ranked[::-1])[::-1]
    ranked=np.minimum(ranked,1.0)
    back=np.empty(len(q)); back[order]=ranked
    out[idx]=back
    return out

def smd_pooled(x1,x0):
    x1=np.asarray(x1,dtype=float); x0=np.asarray(x0,dtype=float)
    n1,n0=len(x1),len(x0)
    if n1<2 or n0<2: return np.nan
    s1=np.var(x1,ddof=1); s0=np.var(x0,ddof=1)
    sp=((n1-1)*s1+(n0-1)*s0)/(n1+n0-2)
    if sp<=0: return 0.0 if np.mean(x1)==np.mean(x0) else np.nan
    return (np.mean(x1)-np.mean(x0))/math.sqrt(sp)

def hypergeom_sf(N,K,n,x):
    lo=max(x, max(0,n-(N-K))); hi=min(K,n)
    logs=[]
    denom=math.lgamma(N+1)-math.lgamma(n+1)-math.lgamma(N-n+1)
    for k in range(lo,hi+1):
        lp=(math.lgamma(K+1)-math.lgamma(k+1)-math.lgamma(K-k+1)
            +math.lgamma(N-K+1)-math.lgamma(n-k+1)-math.lgamma(N-K-(n-k)+1)
            -denom)
        logs.append(lp)
    m=max(logs)
    return math.exp(m)*sum(math.exp(z-m) for z in logs)

def write_json(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True,default=str)+"\n",encoding="utf-8")

def read_contract(ma_root: Path):
    return json.loads((ma_root/"config/wp01_contract.json").read_text())

def parse_kv(path: Path):
    d={}
    for line in path.read_text(errors="replace").splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            k,v=line.split("=",1); d[k.strip()]=v.strip()
    return d

def require_release(path: Path, locksha: str, gate: str, extra_key=None, extra_value=None):
    if not path.exists(): raise RuntimeError(f"Release token missing: {path}")
    d=parse_kv(path)
    if d.get("decision")!="PASS" or d.get("gate")!=gate:
        raise RuntimeError(f"Invalid release token decision/gate: {path}")
    if d.get("preanalysis_lock_sha256")!=locksha:
        raise RuntimeError("Release token bound to different preanalysis lock")
    if extra_key and d.get(extra_key)!=extra_value:
        raise RuntimeError(f"Release token {extra_key} mismatch")

def exact_set_hash(genes):
    s="\n".join(sorted(map(str,genes)))+"\n"
    return hashlib.sha256(s.encode()).hexdigest()

def summarize_status(df, status_col="status"):
    return df[status_col].value_counts(dropna=False).to_dict()

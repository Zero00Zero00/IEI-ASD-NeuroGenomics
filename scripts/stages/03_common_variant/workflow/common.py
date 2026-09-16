
from pathlib import Path
import csv, gzip, json, hashlib, os, re, math, subprocess, shutil, glob
from datetime import datetime, timezone

CH3 = Path(os.environ.get("CH3_ROOT","/home/h3021/chapter3"))
MA = Path(os.environ.get("MA_ROOT","/home/h3021/chapter3/12_molecular_autism_revision"))
PKG = Path(__file__).resolve().parent.parent
OUT = MA / "03_common_retest" / "MAR03_common_variant_retest_v1.0R2"
MAR02 = MA / "02_raw26_domains" / "MAR02_raw26_compression_v1"
PGC = CH3 / "06_magma_pgc" / "MDV6_v1_1"
SPARK = CH3 / "07_magma_spark" / "MDV7_v1_1"
MDV8 = CH3 / "08_specificity_controls" / "MDV8_v1_1"
PYTHON_BIN = Path(os.environ.get("MAR03_PYTHON", str(MA/"envs"/"mar02_mdv5_legacy"/"bin"/"python")))

def utc():
    return datetime.now(timezone.utc).isoformat()

def ensure_dirs():
    for d in ["logs","work","magma","matched","sensitivity","source_data","state"]:
        (OUT/d).mkdir(parents=True,exist_ok=True)

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def read_tsv(p):
    opener = gzip.open if str(p).endswith(".gz") else open
    with opener(p,"rt",encoding="utf-8-sig",errors="replace",newline="") as f:
        r=csv.DictReader(f,delimiter="\t")
        return r.fieldnames or [], list(r)

def write_tsv(p, fields, rows):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,"w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader()
        for z in rows: w.writerow(z)

def pick(header, names):
    m={x.strip().lower():x for x in header}
    for n in names:
        if n.lower() in m: return m[n.lower()]
    return None

def bh(pvals):
    n=len(pvals)
    if n==0: return []
    order=sorted(range(n),key=lambda i:(float(pvals[i]),i))
    q=[1.0]*n; cur=1.0
    for rank_rev,i in enumerate(reversed(order),1):
        rank=n-rank_rev+1
        val=min(1.0,float(pvals[i])*n/rank)
        cur=min(cur,val); q[i]=cur
    return q

def verify_manifest_entry(manifest,target):
    target=Path(target)
    if not Path(manifest).is_file(): return ("NO_MANIFEST","")
    for line in Path(manifest).read_text(encoding="utf-8",errors="replace").splitlines():
        parts=line.strip().split()
        if len(parts)>=2 and Path(parts[-1].lstrip("*")).name==target.name:
            exp=parts[0]; act=sha256_file(target)
            return ("PASS" if exp==act else "FAIL",exp)
    return ("NOT_LISTED","")

def find_magma():
    candidates=[]
    env=os.environ.get("MAGMA_BIN","").strip()
    if env: candidates.append(Path(env))
    w=shutil.which("magma")
    if w: candidates.append(Path(w))
    # Known historical location from real interface probe.
    candidates.append(Path("/home/h3021/chapter2/2.0_IEI_ASD/2.3_GWAS/2.3.2_magma/code/magma"))
    candidates += [Path(x) for x in glob.glob(str(CH3/".snakemake"/"conda"/"*"/"bin"/"magma"))]
    uniq=[]; seen=set()
    for p in candidates:
        try:
            p=p.resolve()
            if p.is_file() and os.access(p,os.X_OK) and str(p) not in seen:
                seen.add(str(p)); uniq.append(p)
        except: pass
    if not uniq:
        raise RuntimeError("MAGMA executable not found")
    return uniq[0],uniq

def detect_map(path):
    h,rows=read_tsv(path)
    symbol=pick(h,["gene","symbol","gene_symbol","hgnc_symbol","approved_symbol","HGNC_symbol"])
    hgnc=pick(h,["hgnc_id","hgnc","hgnc_gene_id"])
    gid=pick(h,["gene_id","entrez_id","entrez","magma_gene_id","MAGMA_GENE","entrezgene","ncbi_gene_id"])
    if not gid:
        raise RuntimeError(f"cannot resolve gene ID column in {path}: {h}")
    by_symbol={}; by_hgnc={}
    for z in rows:
        g=(z.get(gid) or "").strip()
        if not g: continue
        if symbol and (z.get(symbol) or "").strip():
            by_symbol[(z.get(symbol) or "").strip().upper()]=g
        if hgnc and (z.get(hgnc) or "").strip():
            by_hgnc[(z.get(hgnc) or "").strip().upper()]=g
    return {"gene_id_col":gid,"symbol_col":symbol,"hgnc_col":hgnc,
            "by_symbol":by_symbol,"by_hgnc":by_hgnc}

def resolve_pathway_membership():
    p=CH3/"02_pathways"/"03_pathway_membership.tsv.gz"
    h,rows=read_tsv(p)
    source=pick(h,["source","database","db"])
    pid=pick(h,["pathway_id","term_id","id"])
    gene=pick(h,["gene","symbol","gene_symbol","hgnc_symbol","approved_symbol","HGNC_symbol"])
    hgnc=pick(h,["hgnc_id","hgnc"])
    key=pick(h,["pathway_key"])
    if not gene and not hgnc:
        raise RuntimeError(f"pathway membership lacks gene/HGNC columns: {h}")
    return p,h,rows,source,pid,gene,hgnc,key

def parse_gsa(path):
    raw=[]
    with open(path,encoding="utf-8",errors="replace") as f:
        for line in f:
            s=line.rstrip("\n")
            if not s.strip(): continue
            # MAGMA headers may be comment-prefixed.
            raw.append(s.lstrip("#").strip())
    hi=None; header=None
    for i,s in enumerate(raw):
        toks=re.split(r"\s+",s)
        up={x.upper() for x in toks}
        if "BETA" in up and "SE" in up and "P" in up and ("VARIABLE" in up or "FULL_NAME" in up):
            hi=i; header=toks; break
    if hi is None:
        raise RuntimeError(f"cannot find MAGMA GSA header in {path}")
    out=[]
    for s in raw[hi+1:]:
        vals=re.split(r"\s+",s)
        if len(vals)<len(header): continue
        z=dict(zip(header,vals[:len(header)]))
        # skip non-result summary rows
        if "BETA" not in z or "P" not in z: continue
        out.append(z)
    if not out:
        raise RuntimeError(f"no GSA result rows parsed from {path}")
    return header,out

def get_result_row(rows,target):
    for col in ["FULL_NAME","VARIABLE","NAME","SET"]:
        for z in rows:
            if str(z.get(col,""))==target: return z
    raise RuntimeError(f"target {target} absent from GSA output")

def result_values(z):
    def g(keys):
        for k in keys:
            if k in z and str(z[k])!="": return z[k]
        raise KeyError(keys)
    return int(float(g(["NGENES"]))), float(g(["BETA","B"])), float(g(["BETA_STD"])), float(g(["SE"])), float(g(["P","P_VALUE","PVALUE"]))

def run_cmd(cmd,log):
    log=Path(log); log.parent.mkdir(parents=True,exist_ok=True)
    cp=subprocess.run([str(x) for x in cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    log.write_text("COMMAND: "+" ".join(map(str,cmd))+"\n\n"+cp.stdout,encoding="utf-8")
    if cp.returncode!=0:
        raise RuntimeError(f"command failed rc={cp.returncode}; see {log}")
    return cp

def gsa_output(prefix):
    p=Path(str(prefix)+".gsa.out")
    if not p.is_file(): raise RuntimeError(f"missing MAGMA GSA output: {p}")
    return p

def read_gene_scores(path):
    h,rows=read_tsv(path)
    gid=pick(h,["MAGMA_GENE","GENE","gene_id","entrez_id"])
    zc=pick(h,["ZSTAT","zstat","Z","gene_z"])
    pc=pick(h,["P","p","p_value"])
    if not gid or (not zc and not pc):
        raise RuntimeError(f"gene results schema unresolved: {path} {h}")
    out={}
    from scipy.stats import norm
    for z in rows:
        g=(z.get(gid) or "").strip()
        if not g: continue
        try:
            if zc and (z.get(zc) or "")!="":
                val=float(z[zc])
            else:
                p=min(max(float(z[pc]),1e-300),1-1e-16)
                val=float(norm.isf(p))
            if math.isfinite(val): out[g]=val
        except: pass
    return out

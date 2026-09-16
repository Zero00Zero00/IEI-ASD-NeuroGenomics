#!/usr/bin/env python3
from __future__ import annotations
import gzip, hashlib, json, os, re, shutil, subprocess
from pathlib import Path
from typing import Iterable
import numpy as np
import pandas as pd
import yaml


def load_cfg(path: str | Path):
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def root_path(cfg, p):
    p=Path(p)
    return p if p.is_absolute() else Path(cfg["project_root"])/p


def out_dir(cfg):
    p=root_path(cfg,cfg["output_dir"]); p.mkdir(parents=True,exist_ok=True); return p


def prov_dir(cfg):
    p=root_path(cfg,cfg["provenance_dir"]); p.mkdir(parents=True,exist_ok=True); return p


def sha256_file(path, block=8*1024*1024):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(block),b""): h.update(b)
    return h.hexdigest()


def write_json(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with open(path,"w",encoding="utf-8") as f: json.dump(obj,f,indent=2,sort_keys=True)


def write_text(path,text):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with open(path,"w",encoding="utf-8") as f: f.write(text)


def require_file(path,label=None):
    p=Path(path)
    if not p.is_file() or p.stat().st_size==0: raise RuntimeError(f"Missing/empty {label or p}: {p}")
    return p


def verify_sha256_manifest(path, root):
    path=Path(path); root=Path(root); rows=[]
    if not path.exists(): return [{"entry":"MANIFEST","status":"FAIL","detail":f"missing {path}"}]
    for raw in path.read_text(encoding="utf-8",errors="replace").splitlines():
        raw=raw.strip()
        if not raw or raw.startswith("#"): continue
        m=re.match(r"^([0-9a-fA-F]{64})\s+\*?(.*)$",raw)
        if not m:
            rows.append({"entry":raw,"status":"FAIL","detail":"unparseable checksum line"}); continue
        exp,fn=m.group(1).lower(),m.group(2).strip()
        fp=Path(fn)
        if not fp.is_absolute(): fp=root/fp
        if not fp.exists(): rows.append({"entry":fn,"status":"FAIL","detail":"missing"}); continue
        got=sha256_file(fp)
        rows.append({"entry":fn,"status":"PASS" if got==exp else "FAIL","detail":got})
    return rows


def discover_1kg_prefix(extracted_dir):
    root=Path(extracted_dir)
    beds={p.with_suffix("") for p in root.rglob("*.bed")}
    bims={p.with_suffix("") for p in root.rglob("*.bim")}
    fams={p.with_suffix("") for p in root.rglob("*.fam")}
    common=sorted(beds & bims & fams)
    if not common: raise RuntimeError(f"No common BED/BIM/FAM prefix under {root}")
    preferred=[p for p in common if p.name=="g1000_eur"]
    if len(preferred)==1: return str(preferred[0])
    if len(common)==1: return str(common[0])
    raise RuntimeError("Multiple 1KG prefixes found; set ld_reference.prefix explicitly: "+", ".join(map(str,common[:20])))


def resolve_magma(binary):
    if binary and binary!="AUTO":
        p=Path(binary)
        if p.exists(): return str(p.resolve())
        q=shutil.which(binary)
        if q: return q
        raise RuntimeError(f"MAGMA binary not found: {binary}")
    for candidate in ["magma","magma_v1.10"]:
        q=shutil.which(candidate)
        if q: return q
    raise RuntimeError("MAGMA v1.10 not found in PATH. Install the official binary, then rerun preflight.")


def magma_version(binary):
    cp=subprocess.run([binary,"--version"],capture_output=True,text=True)
    txt=(cp.stdout+"\n"+cp.stderr).strip()
    if not txt:
        cp=subprocess.run([binary,"-v"],capture_output=True,text=True)
        txt=(cp.stdout+"\n"+cp.stderr).strip()
    return txt


def bh(pvals):
    x=np.asarray(pvals,dtype=float); n=len(x); out=np.full(n,np.nan)
    ok=np.isfinite(x)
    vals=x[ok]; idx=np.where(ok)[0]
    if not len(vals): return out
    order=np.argsort(vals); ranked=vals[order]
    q=ranked*len(ranked)/np.arange(1,len(ranked)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]
    q=np.minimum(q,1.0)
    oo=np.empty_like(q); oo[order]=q
    out[idx]=oo
    return out


def read_magma_table(path):
    # MAGMA text outputs can contain comments and variable whitespace.
    rows=[]
    with open(path,encoding="utf-8",errors="replace") as f:
        lines=[x.strip() for x in f if x.strip() and not x.startswith("#")]
    if not lines: raise RuntimeError(f"Empty MAGMA output: {path}")
    header=re.split(r"\s+",lines[0])
    for line in lines[1:]:
        vals=re.split(r"\s+",line)
        if len(vals)<len(header): continue
        if len(vals)>len(header):
            # FULL_NAME can contain no whitespace in our sets, so extras are unexpected.
            vals=vals[:len(header)]
        rows.append(vals)
    return pd.DataFrame(rows,columns=header)


def hgnc_numeric(x):
    m=re.fullmatch(r"HGNC:(\d+)",str(x).strip())
    if not m: raise ValueError(f"Invalid HGNC id: {x}")
    return m.group(1)


def norm_chr(x):
    s=str(x).strip().upper().replace("CHR","")
    if s=="X": return 23
    if s=="Y": return 24
    if s in {"M","MT"}: return 25
    try: return int(float(s))
    except: return None


def complement(a):
    tab=str.maketrans("ACGTacgt","TGCAtgca")
    return str(a).translate(tab)


def allele_compatible(a1,a2,b1,b2):
    A={str(a1).upper(),str(a2).upper()}; B={str(b1).upper(),str(b2).upper()}
    if A==B: return True,"EXACT_SET"
    if all(len(x)==1 and x in "ACGT" for x in A|B):
        C={complement(x).upper() for x in A}
        if C==B: return True,"COMPLEMENT_SET"
    return False,"MISMATCH"


def run(cmd, log_path, env=None):
    Path(log_path).parent.mkdir(parents=True,exist_ok=True)
    with open(log_path,"w",encoding="utf-8") as log:
        log.write("COMMAND: "+" ".join(map(str,cmd))+"\n\n")
        cp=subprocess.run(list(map(str,cmd)),stdout=log,stderr=subprocess.STDOUT,text=True,env=env)
    if cp.returncode!=0: raise RuntimeError(f"Command failed ({cp.returncode}); see {log_path}")

#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json
from pathlib import Path

ACTIVE_VERSION="SOLVER_LOCKED_v1.3.2"

def sha256_file(path: Path, block=1<<20):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        while True:
            b=f.read(block)
            if not b: break
            h.update(b)
    return h.hexdigest()

def parse_kv(path: Path):
    d={}
    for line in path.read_text(errors="replace").splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            k,v=line.split("=",1)
            d[k.strip()]=v.strip()
    return d

def verify_stable(ma: Path, out: Path, prelock_sha: str):
    manifest=out/"WP01_SOLVER_LOCKED_IMPLEMENTATION_v1.3.2.json"
    release=out/"WP01_SOLVER_LOCKED_IMPLEMENTATION_RELEASE_v1.3.2.txt"
    if not manifest.exists() or not release.exists():
        raise RuntimeError("WP01 solver-locked v1.3.0 manifest/release missing")

    rel=parse_kv(release)
    msh=sha256_file(manifest)
    if rel.get("decision")!="PASS" or rel.get("stage")!="WP01" or rel.get("implementation")!="SOLVER_LOCKED_v1.3.2":
        raise RuntimeError("Invalid WP01 solver-locked implementation release")
    if rel.get("preanalysis_lock_sha256")!=prelock_sha:
        raise RuntimeError("Solver-locked release bound to a different preanalysis lock")
    if rel.get("implementation_manifest_sha256")!=msh:
        raise RuntimeError("Solver-locked release bound to a different implementation manifest")

    obj=json.loads(manifest.read_text())
    if obj.get("status")!="READY_FOR_RELEASE" or obj.get("version")!="SOLVER_LOCKED_v1.3.2":
        raise RuntimeError("Solver-locked manifest status/version invalid")
    if obj.get("preanalysis_lock_sha256")!=prelock_sha:
        raise RuntimeError("Solver-locked manifest bound to a different preanalysis lock")

    changed=[]
    for relpath,expected in obj["active_code_sha256"].items():
        p=ma/relpath
        if not p.exists():
            changed.append(relpath+":MISSING")
        elif sha256_file(p)!=expected:
            changed.append(relpath+":"+sha256_file(p))
    if changed:
        raise RuntimeError("Active WP01 code drift after v1.3.0 release: "+";".join(changed))

    for relpath,expected in obj["external_solver_authority_sha256"].items():
        p=Path(obj["project_root"])/relpath
        if not p.exists() or sha256_file(p)!=expected:
            raise RuntimeError("Legacy solver authority drift: "+relpath)

    for relpath,expected in obj["historical_release_sha256"].items():
        p=out/relpath
        if not p.exists() or sha256_file(p)!=expected:
            raise RuntimeError("Historical WP01 release lineage drift: "+relpath)
    return obj

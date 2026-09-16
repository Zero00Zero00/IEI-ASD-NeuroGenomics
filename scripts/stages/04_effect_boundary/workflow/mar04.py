
from pathlib import Path
import csv, json, hashlib, os, math, tarfile
from statistics import NormalDist
from datetime import datetime, timezone

MA = Path(os.environ.get("MA_ROOT","/home/h3021/chapter3/12_molecular_autism_revision"))
PKG = Path(__file__).resolve().parent.parent
OUT = MA / "04_effect_boundary" / "MAR04_effect_boundary_v1"
MAR03_ENV = os.environ.get("MAR03_FREEZE_ROOT","").strip()

def utc():
    return datetime.now(timezone.utc).isoformat()

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def read_tsv(p):
    with open(p,encoding="utf-8-sig",errors="replace",newline="") as f:
        r=csv.DictReader(f,delimiter="\t")
        return r.fieldnames or [], list(r)

def write_tsv(p,fields,rows):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,"w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

def resolve_mar03():
    if MAR03_ENV:
        p=Path(MAR03_ENV)
        if p.is_dir(): return p.resolve()
    candidates=[
        MA/"MAR03_FREEZE_RELEASE_v1.0R2",
        MA/"03_common_retest"/"MAR03_FREEZE_RELEASE_v1.0R2",
    ]
    candidates += [p for p in MA.rglob("MAR03_FREEZE_RELEASE_v1.0R2") if p.is_dir()]
    uniq=[]
    seen=set()
    for p in candidates:
        try:
            rp=p.resolve()
            if rp.is_dir() and str(rp) not in seen:
                seen.add(str(rp)); uniq.append(rp)
        except: pass
    if len(uniq)!=1:
        raise SystemExit("HOLD: expected exactly one MAR03_FREEZE_RELEASE_v1.0R2; set MAR03_FREEZE_ROOT explicitly")
    return uniq[0]

def parse_pass(p):
    d={}
    for line in Path(p).read_text(encoding="utf-8").splitlines():
        if "=" in line:
            k,v=line.split("=",1); d[k]=v
    return d

import sys
ND=NormalDist()
Z975=ND.inv_cdf(0.975)
Z95=ND.inv_cdf(0.95)
Z80=ND.inv_cdf(0.80)

def fam_const(m):
    alpha=0.05/m
    zfam=ND.inv_cdf(1-alpha)
    return alpha,zfam,zfam+Z80

def preflight():
    OUT.mkdir(parents=True,exist_ok=True)
    mar03=resolve_mar03()
    checks=[]
    def ck(cid,ok,obs,exp):
        checks.append({"check_id":cid,"status":"PASS" if ok else "HOLD","observed":obs,"expected":exp})
    passf=mar03/"manifest"/"MAR03_FREEZE_PASS.txt"
    ck("MAR03_PASS_exists",passf.is_file(),str(passf),"exists")
    if passf.is_file():
        d=parse_pass(passf); ck("MAR03_status",d.get("status")=="PASS",d.get("status",""),"PASS")
    inputs={
      "primary":mar03/"source_data"/"MAR03_SOURCE_PRIMARY_HIGHER_ORDER.tsv",
      "standalone":mar03/"source_data"/"MAR03_SOURCE_STANDALONE_RESULTS.tsv",
      "conditional":mar03/"source_data"/"MAR03_SOURCE_CONDITIONAL_RESULTS.tsv",
    }
    for k,p in inputs.items(): ck(k+"_exists",p.is_file(),str(p),"exists")
    if all(p.is_file() for p in inputs.values()):
        _,a=read_tsv(inputs["primary"]); _,b=read_tsv(inputs["standalone"]); _,c=read_tsv(inputs["conditional"])
        ck("primary_rows",len(a)==14,len(a),14)
        ck("standalone_rows",len(b)==26,len(b),26)
        ck("conditional_rows",len(c)==14,len(c),14)
        ck("primary_dataset_counts",sorted([sum(z["dataset"]==d for z in a) for d in ["PGC","SPARK"]])==[7,7],
           f"PGC={sum(z['dataset']=='PGC' for z in a)},SPARK={sum(z['dataset']=='SPARK' for z in a)}","7/7")
        ck("standalone_dataset_counts",sorted([sum(z["dataset"]==d for z in b) for d in ["PGC","SPARK"]])==[13,13],
           f"PGC={sum(z['dataset']=='PGC' for z in b)},SPARK={sum(z['dataset']=='SPARK' for z in b)}","13/13")
    write_tsv(OUT/"MAR04_preflight_report.tsv",["check_id","status","observed","expected"],checks)
    bad=[z for z in checks if z["status"]!="PASS"]
    resolved={"MAR03_ROOT":str(mar03),"inputs":{k:str(v) for k,v in inputs.items()},
              "input_hashes":{k:sha256_file(v) for k,v in inputs.items() if v.is_file()}}
    (OUT/"MAR04_resolved_inputs.json").write_text(json.dumps(resolved,indent=2),encoding="utf-8")
    if bad: raise SystemExit(f"HOLD: MAR04 preflight {len(bad)} failures")
    print("MAR04_PREFLIGHT=PASS")

def freeze():
    if not (OUT/"MAR04_preflight_report.tsv").is_file(): preflight()
    formula=json.loads((PKG/"config"/"MAR04_FORMULA_LOCK.json").read_text())
    resolved=json.loads((OUT/"MAR04_resolved_inputs.json").read_text())
    content={
      "stage":"MAR04_EFFECT_BOUNDARY","version":"1.0R1",
      "formula_lock":formula,
      "source_input_hashes":resolved["input_hashes"],
      "source_mar03_root":resolved["MAR03_ROOT"],
      "SESOI_prespecified":"NO","equivalence_claim_allowed":"NO",
      "conditional_status":"SUPPORTIVE_PENDING_MAR05",
    }
    lock=OUT/"MAR04_preanalysis_lock.json"
    if lock.is_file():
        old=json.loads(lock.read_text())
        old.pop("timestamp_utc",None); new=dict(content)
        if old!=new: raise SystemExit("HOLD: existing MAR04 lock differs; do not overwrite")
        print("MAR04_FREEZE=PASS_REUSED")
        return
    content["timestamp_utc"]=utc()
    lock.write_text(json.dumps(content,indent=2),encoding="utf-8")
    print("MAR04_FREEZE=PASS")
    print("HOLD_FOR_MANUAL_GATE")

def release():
    lock=OUT/"MAR04_preanalysis_lock.json"
    if not lock.is_file(): raise SystemExit("HOLD: freeze first")
    token=OUT/"MAR04_STOPA_RELEASE.txt"
    token.write_text(
      "decision=PASS\ntrack=MOLECULAR_AUTISM_REVISION\nstage=MAR04\ngate=A\n"
      "reviewer=MANUAL_EXPERT_REVIEW\n"
      f"timestamp_utc={utc()}\npreanalysis_lock_sha256={sha256_file(lock)}\n"
      "notes=Prespecified deterministic effect-boundary formulas; no SESOI/equivalence margin.\n",
      encoding="utf-8"
    )
    print("MAR04_GATE_A=PASS")

def derive(rows,family,m,role,mar05):
    alpha=0.05/m; zfam=ND.inv_cdf(1-alpha); mult=zfam+Z80; out=[]
    for r in rows:
        beta=float(r["beta"]); se=float(r["SE"]); q=float(r["q_BH"])
        lo=beta-Z975*se; hi=beta+Z975*se; up=beta+Z95*se; mde=mult*se
        det="DETECTED_POSITIVE" if beta>0 and q<0.05 else "NOT_DETECTED"
        out.append({
          "dataset":r["dataset"],"family":family,"target_id":r["target_id"],
          "biological_label":r["biological_label"],"inference_role":role,
          "mar05_gate_required":mar05,"frozen_n":r["frozen_n"],"mapped_n":r["mapped_n"],
          "coverage":r["coverage"],"beta":beta,"SE":se,"P_one":r["P_one"],"q_BH":r["q_BH"],
          "CI95_low":lo,"CI95_high":hi,"upper95_one":up,
          "family_m":m,"alpha_family_one":alpha,"z_family_one":zfam,"z_power_80":Z80,
          "MDE80":mde,"SESOI_prespecified":"NO","equivalence_allowed":"NO",
          "detection_class":det
        })
    return out

def analysis():
    lock=OUT/"MAR04_preanalysis_lock.json"; tok=OUT/"MAR04_STOPA_RELEASE.txt"
    if not lock.is_file() or not tok.is_file(): raise SystemExit("HOLD: freeze/release required")
    td=parse_pass(tok)
    if td.get("preanalysis_lock_sha256")!=sha256_file(lock): raise SystemExit("HOLD: release-lock mismatch")
    resolved=json.loads((OUT/"MAR04_resolved_inputs.json").read_text())
    _,a=read_tsv(resolved["inputs"]["primary"]); _,b=read_tsv(resolved["inputs"]["standalone"]); _,c=read_tsv(resolved["inputs"]["conditional"])
    pa=derive(a,"A_PRIMARY_HIGHER_ORDER",7,"PRIMARY","NO")
    pb=derive(b,"B_SECONDARY_STANDALONE",13,"SECONDARY_COMPLETENESS","NO")
    pc=derive(c,"A_CONDITIONAL_SUPPORTIVE",7,"SUPPORTIVE_PENDING_MAR05","YES")
    allr=pa+pb+pc
    fields=list(allr[0].keys())
    write_tsv(OUT/"MAR04_effect_boundary.tsv",fields,allr)
    write_tsv(OUT/"MAR04_primary_effect_boundary.tsv",fields,pa)
    write_tsv(OUT/"MAR04_standalone_effect_boundary.tsv",fields,pb)
    write_tsv(OUT/"MAR04_conditional_supportive.tsv",fields,pc)
    forest=[]
    for r in pa+pb:
        forest.append({k:r[k] for k in ["dataset","family","target_id","biological_label","beta","CI95_low","CI95_high","upper95_one","MDE80","q_BH","detection_class"]})
    write_tsv(OUT/"MAR04_forest_source.tsv",list(forest[0].keys()),forest)
    power=[]
    for r in pa+pb:
        power.append({
          "dataset":r["dataset"],"family":r["family"],"target_id":r["target_id"],"biological_label":r["biological_label"],
          "beta":r["beta"],"SE":r["SE"],"CI95_low":r["CI95_low"],"CI95_high":r["CI95_high"],
          "upper95_one":r["upper95_one"],"MDE80":r["MDE80"],"q_BH":r["q_BH"],
          "equivalence_allowed":"NO",
          "allowed_language":"No reproducible positive enrichment detected at the tested resolution; nonsignificance is not evidence of equivalence."
        })
    write_tsv(OUT/"MAR04_power_interpretation.tsv",list(power[0].keys()),power)
    print("MAR04_ANALYSIS=PASS")

def postflight():
    req=["MAR04_effect_boundary.tsv","MAR04_primary_effect_boundary.tsv","MAR04_standalone_effect_boundary.tsv",
         "MAR04_conditional_supportive.tsv","MAR04_forest_source.tsv","MAR04_power_interpretation.tsv"]
    checks=[]
    def ck(cid,ok,obs,exp): checks.append({"check_id":cid,"status":"PASS" if ok else "HOLD","observed":obs,"expected":exp})
    for f in req: ck("exists_"+f,(OUT/f).is_file(),f,"exists")
    if all((OUT/f).is_file() for f in req):
        _,a=read_tsv(OUT/"MAR04_primary_effect_boundary.tsv"); _,b=read_tsv(OUT/"MAR04_standalone_effect_boundary.tsv"); _,c=read_tsv(OUT/"MAR04_conditional_supportive.tsv")
        ck("primary_n",len(a)==14,len(a),14); ck("standalone_n",len(b)==26,len(b),26); ck("conditional_n",len(c)==14,len(c),14)
        ck("primary_m7",all(str(z["family_m"])=="7" for z in a),7,7)
        ck("standalone_m13",all(str(z["family_m"])=="13" for z in b),13,13)
        ck("equivalence_no",all(z["equivalence_allowed"]=="NO" for z in a+b+c),"NO","NO")
        # formula recomputation
        tol=1e-9
        ok=True
        for z in a+b+c:
            beta=float(z["beta"]); se=float(z["SE"])
            if abs(float(z["CI95_low"])-(beta-Z975*se))>tol: ok=False
            if abs(float(z["CI95_high"])-(beta+Z975*se))>tol: ok=False
            if abs(float(z["upper95_one"])-(beta+Z95*se))>tol: ok=False
        ck("formula_recompute",ok,"","exact within 1e-9")
    write_tsv(OUT/"MAR04_postflight_checks.tsv",["check_id","status","observed","expected"],checks)
    bad=[z for z in checks if z["status"]!="PASS"]
    if bad: raise SystemExit(f"HOLD: MAR04 postflight {len(bad)} failures")
    # source_data
    sd=OUT/"source_data"; sd.mkdir(exist_ok=True)
    for src,dst in [
      ("MAR04_primary_effect_boundary.tsv","MAR04_SOURCE_PRIMARY_EFFECT_BOUNDARY.tsv"),
      ("MAR04_standalone_effect_boundary.tsv","MAR04_SOURCE_STANDALONE_EFFECT_BOUNDARY.tsv"),
      ("MAR04_conditional_supportive.tsv","MAR04_SOURCE_CONDITIONAL_SUPPORTIVE.tsv"),
      ("MAR04_forest_source.tsv","MAR04_SOURCE_FOREST_PLOT.tsv"),
      ("MAR04_power_interpretation.tsv","MAR04_SOURCE_POWER_INTERPRETATION.tsv"),
    ]:
        import shutil; shutil.copy2(OUT/src,sd/dst)
    files=[p for p in OUT.rglob("*") if p.is_file() and p.name not in {"MAR04_checksums.sha256","MAR04_PASS.txt"}]
    lock={"stage":"MAR04","status":"PASS","timestamp_utc":utc(),"preanalysis_lock_sha256":sha256_file(OUT/"MAR04_preanalysis_lock.json"),
          "output_hashes":{str(p.relative_to(OUT)):sha256_file(p) for p in files}}
    (OUT/"MAR04_analysis_lock.json").write_text(json.dumps(lock,indent=2),encoding="utf-8")
    lines=[f"{sha256_file(p)}  {p.relative_to(MA)}" for p in sorted([p for p in OUT.rglob("*") if p.is_file() and p.name!="MAR04_checksums.sha256"])]
    (OUT/"MAR04_checksums.sha256").write_text("\n".join(lines)+"\n",encoding="utf-8")
    (OUT/"MAR04_PASS.txt").write_text("status=PASS\nstage=MAR04_EFFECT_BOUNDARY\nprimary_family_m=7\nstandalone_family_m=13\nequivalence_allowed=NO\nconditional_status=SUPPORTIVE_PENDING_MAR05\nnext_stage=MAR05\n",encoding="utf-8")
    print("MAR04_POSTFLIGHT=PASS")

def pack_review():
    import tarfile
    if not (OUT/"MAR04_PASS.txt").is_file(): raise SystemExit("HOLD: MAR04_PASS missing")
    handoff=MA/"90_handoff"; handoff.mkdir(parents=True,exist_ok=True)
    ts=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    p=handoff/f"MAR04_REVIEW_{ts}.tar.gz"
    include=["MAR04_preflight_report.tsv","MAR04_resolved_inputs.json","MAR04_preanalysis_lock.json","MAR04_STOPA_RELEASE.txt",
             "MAR04_effect_boundary.tsv","MAR04_primary_effect_boundary.tsv","MAR04_standalone_effect_boundary.tsv",
             "MAR04_conditional_supportive.tsv","MAR04_forest_source.tsv","MAR04_power_interpretation.tsv",
             "MAR04_postflight_checks.tsv","MAR04_analysis_lock.json","MAR04_checksums.sha256","MAR04_PASS.txt"]
    with tarfile.open(p,"w:gz") as tf:
        for n in include:
            f=OUT/n
            if f.exists(): tf.add(f,arcname=n)
        for f in (OUT/"source_data").glob("*"): tf.add(f,arcname="source_data/"+f.name)
    sha=sha256_file(p); Path(str(p)+".sha256").write_text(f"{sha}  {p.name}\n",encoding="utf-8")
    print(f"REVIEW_PACKAGE={p}"); print(f"SHA256={sha}"); print("PACK_REVIEW=PASS")

cmd=sys.argv[1] if len(sys.argv)>1 else "status"
if cmd=="preflight": preflight()
elif cmd=="freeze": freeze()
elif cmd=="release": release()
elif cmd=="analysis": analysis()
elif cmd=="postflight": postflight()
elif cmd=="pack-review": pack_review()
elif cmd=="status":
    print("MAR04_OUT=",OUT)
    for f in ["MAR04_preflight_report.tsv","MAR04_preanalysis_lock.json","MAR04_STOPA_RELEASE.txt","MAR04_effect_boundary.tsv","MAR04_postflight_checks.tsv","MAR04_PASS.txt"]:
        print("OK" if (OUT/f).is_file() else "MISS",OUT/f)
else:
    raise SystemExit("Unknown command")

#!/usr/bin/env python3
import argparse, csv, gzip, hashlib, json, os, platform, re, shutil, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path

PKG_VERSION = "1.0.1"
DEFAULT_ROOT = "/home/h3021/chapter3"


def utcnow():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256_file(path, chunk=1024*1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else open(path, "rt", encoding="utf-8", errors="replace")


def tsv_info(path):
    rows = 0
    header = []
    with open_text(path) as f:
        first = f.readline().rstrip("\n\r")
        if first:
            header = first.split("\t")
        for _ in f:
            rows += 1
    return rows, header


def norm(s):
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())


def choose_col(header, candidates):
    m = {norm(h): h for h in header}
    for c in candidates:
        if norm(c) in m:
            return m[norm(c)]
    return None


def truthy(v):
    return str(v).strip().lower() in {"1", "true", "t", "yes", "y", "present", "in", "coreseed"}


def write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def read_manifest(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def command_version(cmd):
    exe = shutil.which(cmd)
    if not exe:
        return {"tool": cmd, "available": "NO", "path": "", "version": "NOT_FOUND"}
    probes = {
        "python3": [exe, "--version"], "Rscript": [exe, "--version"], "bash": [exe, "--version"],
        "git": [exe, "--version"], "tar": [exe, "--version"], "sha256sum": [exe, "--version"], "magma": [exe]
    }
    try:
        p = subprocess.run(probes.get(cmd, [exe, "--version"]), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=10)
        first = (p.stdout or "").strip().splitlines()[0] if (p.stdout or "").strip() else f"exit={p.returncode}"
    except Exception as e:
        first = f"ERROR:{e}"
    return {"tool": cmd, "available": "YES", "path": exe, "version": first[:500]}


def stage_paths(root):
    return {
      "MDV0": root/"00_provenance/MDV0_PASS.txt",
      "MDV1": root/"01_gene_universe/MDV1_PASS.txt",
      "PATHWAY": root/"02_pathways/PATHWAY_PREP_PASS.txt",
      "MDV2": root/"03_intersection/MDV2_PASS.txt",
      "MDV3": root/"04_bias_empirical/MDV3_PASS.txt",
      "MDV4": root/"04_bias_empirical/MDV4_empirical/MDV4_PASS.txt",
      "MDV5": root/"05_domains/MDV5_v1_1/MDV5_PASS.txt",
      "MDV6": root/"06_magma_pgc/MDV6_v1_1/MDV6_PASS.txt",
      "MDV7": root/"07_magma_spark/MDV7_v1_1/MDV7_PASS.txt",
      "MDV8": root/"08_specificity_controls/MDV8_v1_1/MDV8_PASS.txt",
    }


def ensure_dirs(ma):
    for d in ["00_contract", "01_anchor_ingest", "02_raw26_domains", "03_common_retest", "04_effect_boundary", "05_magma_collinearity", "06_matched_null_audit", "07_chain_calibration", "08_driver_influence", "09_figures_tables", "10_manuscript", "90_handoff", "91_reproducibility", "bin", "config", "workflow/scripts", "logs"]:
        (ma/d).mkdir(parents=True, exist_ok=True)


def preflight(root, ma, cfgdir):
    ensure_dirs(ma)
    out = ma/"00_contract"
    manifest = read_manifest(cfgdir/"ma_r00_authority_manifest.tsv")
    checks=[]
    root_markers=["00_provenance","01_gene_universe","02_pathways","03_intersection","04_bias_empirical","05_domains","06_magma_pgc","07_magma_spark","08_specificity_controls","09_figures_tables","10_reproducibility","config","workflow"]
    for d in root_markers:
        p=root/d
        checks.append({"check_id":f"DIR::{d}","severity":"CRITICAL","status":"PASS" if p.is_dir() else "FAIL","details":str(p)})
    for row in manifest:
        p=root/row["relpath"]
        req=row["required"]=="1"
        exists=p.exists() and (p.is_file() or p.is_dir())
        nonzero=(p.is_dir() or (p.is_file() and p.stat().st_size>0)) if exists else False
        status="PASS" if exists and nonzero else ("FAIL" if req else "CAUTION")
        checks.append({"check_id":f"AUTH::{row['role']}","severity":"CRITICAL" if req else "INFO","status":status,"details":f"{p}; required={req}"})
    for s,p in stage_paths(root).items():
        checks.append({"check_id":f"STAGEPASS::{s}","severity":"CRITICAL","status":"PASS" if p.is_file() and p.stat().st_size>0 else "FAIL","details":str(p)})
    for cmd in ["python3","bash","tar","sha256sum","git","Rscript"]:
        cv=command_version(cmd)
        sev="CRITICAL" if cmd in {"python3","bash","tar","sha256sum"} else "INFO"
        checks.append({"check_id":f"CMD::{cmd}","severity":sev,"status":"PASS" if cv["available"]=="YES" else ("FAIL" if sev=="CRITICAL" else "CAUTION"),"details":f"{cv['path']} | {cv['version']}"})
    # MA output containment guard
    try:
        ma.resolve().relative_to(root.resolve())
        contain=True
    except Exception:
        contain=False
    checks.append({"check_id":"WRITE_GUARD::MA_ROOT_UNDER_CH3","severity":"CRITICAL","status":"PASS" if contain else "FAIL","details":str(ma)})
    if ma.resolve() in [ (root/d).resolve() for d in root_markers if (root/d).exists() ]:
        checks.append({"check_id":"WRITE_GUARD::NOT_LEGACY_DIR","severity":"CRITICAL","status":"FAIL","details":"MA_ROOT collides with legacy directory"})
    else:
        checks.append({"check_id":"WRITE_GUARD::NOT_LEGACY_DIR","severity":"CRITICAL","status":"PASS","details":str(ma)})
    write_tsv(out/"MAR00_preflight_report.tsv",checks,["check_id","severity","status","details"])
    fails=[x for x in checks if x["severity"]=="CRITICAL" and x["status"]!="PASS"]
    summary={"phase":"preflight","timestamp_utc":utcnow(),"root":str(root),"ma_root":str(ma),"package_version":PKG_VERSION,"critical_failures":len(fails),"status":"PASS" if not fails else "FAIL"}
    with open(out/"MAR00_preflight_summary.json","w") as f: json.dump(summary,f,indent=2)
    marker=out/("MAR00_PREFLIGHT_PASS.txt" if not fails else "MAR00_PREFLIGHT_FAIL.txt")
    marker.write_text(f"{summary['status']}\n{summary['timestamp_utc']}\ncritical_failures={len(fails)}\n")
    return 0 if not fails else 2


def freeze(root, ma, cfgdir):
    out=ma/"00_contract"
    if not (out/"MAR00_PREFLIGHT_PASS.txt").exists():
        print("FATAL: preflight PASS missing",file=sys.stderr); return 2
    manifest=read_manifest(cfgdir/"ma_r00_authority_manifest.tsv")
    resolved=[]; baseline=[]
    for row in manifest:
        p=root/row["relpath"]
        rec=dict(row); rec["absolute_path"]=str(p.resolve()) if p.exists() else str(p)
        if p.is_file():
            st=p.stat(); rec.update({"exists":True,"size_bytes":st.st_size,"mtime_ns":st.st_mtime_ns,"sha256":sha256_file(p),"policy":"READ_ONLY"})
            baseline.append({"role":row["role"],"path":str(p),"size_bytes":st.st_size,"mtime_ns":st.st_mtime_ns,"sha256":rec["sha256"]})
        elif p.is_dir():
            st=p.stat(); rec.update({"exists":True,"size_bytes":0,"mtime_ns":st.st_mtime_ns,"sha256":"DIRECTORY","policy":"READ_ONLY"})
        else:
            rec.update({"exists":False,"size_bytes":0,"mtime_ns":None,"sha256":"MISSING","policy":"READ_ONLY"})
        resolved.append(rec)
    with open(out/"MAR00_resolved_inputs.json","w") as f: json.dump(resolved,f,indent=2)
    write_tsv(out/"MAR00_legacy_hash_baseline.tsv",baseline,["role","path","size_bytes","mtime_ns","sha256"])
    # snapshot relevant code/config provenance without copying data
    code_patterns=["config/mdv*.yaml","Snakefile.mdv*","run_mdv*_staged.sh","mdv*_release_gate.sh","workflow/scripts/mdv*","workflow/rules/mdv*","workflow/envs/mdv*.yaml"]
    code=[]
    for pat in code_patterns:
        for p in sorted(root.glob(pat)):
            if p.is_file():
                code.append({"path":str(p.relative_to(root)),"size_bytes":p.stat().st_size,"mtime_ns":p.stat().st_mtime_ns,"sha256":sha256_file(p)})
    write_tsv(out/"MAR00_legacy_code_manifest.tsv",code,["path","size_bytes","mtime_ns","sha256"])
    lock={
      "track":"MOLECULAR_AUTISM_REVISION","stage":"MA-R00","gate":"A","package_version":PKG_VERSION,
      "timestamp_utc":utcnow(),"ch3_root":str(root.resolve()),"ma_root":str(ma.resolve()),
      "resolved_inputs_sha256":sha256_file(out/"MAR00_resolved_inputs.json"),
      "legacy_baseline_sha256":sha256_file(out/"MAR00_legacy_hash_baseline.tsv"),
      "authority_contract_sha256":sha256_file(cfgdir/"ma_r00_authority_manifest.tsv"),
      "expected_metrics_sha256":sha256_file(cfgdir/"ma_r00_expected_metrics.json"),
      "data_firewall":{"legacy_policy":"READ_ONLY","allowed_write_prefix":str(ma.resolve()),"forbidden_write_prefixes":[str((root/d).resolve()) for d in ["00_provenance","01_gene_universe","02_pathways","03_intersection","04_bias_empirical","05_domains","06_magma_pgc","07_magma_spark","08_specificity_controls","09_figures_tables","10_reproducibility"]]},
      "postfreeze_rule":"analysis requires human release token bound to this lock SHA256"
    }
    with open(out/"MAR00_preanalysis_lock.json","w") as f: json.dump(lock,f,indent=2)
    locksha=sha256_file(out/"MAR00_preanalysis_lock.json")
    (out/"MAR00_STOPA_READY.txt").write_text(f"READY_FOR_MANUAL_REVIEW\n{utcnow()}\npreanalysis_lock_sha256={locksha}\n")
    print(f"STOP A READY: {locksha}")
    return 0


def checksum_audit(root):
    rows=[]
    files=[]
    for pattern in ["03_intersection/MDV2_checksums.sha256","04_bias_empirical/MDV3_checksums.sha256","04_bias_empirical/MDV4_empirical/MDV4_checksums.sha256","05_domains/MDV5_v1_1/MDV5_checksums.sha256","06_magma_pgc/MDV6_v1_1/MDV6_checksums.sha256","07_magma_spark/MDV7_v1_1/MDV7_checksums.sha256","08_specificity_controls/MDV8_v1_1/MDV8_checksums.sha256"]:
        p=root/pattern
        if p.exists(): files.append(p)
    for cf in files:
        with open(cf,encoding="utf-8",errors="replace") as f:
            for line in f:
                line=line.strip()
                if not line or line.startswith("#"): continue
                m=re.match(r"^([0-9a-fA-F]{64})\s+[* ]?(.*)$",line)
                if not m:
                    rows.append({"checksum_file":str(cf.relative_to(root)),"declared_sha256":"","declared_path":line,"resolved_path":"","status":"UNPARSED","details":""}); continue
                h,raw=m.group(1).lower(),m.group(2).strip()
                candidates=[]
                rp=Path(raw)
                if rp.is_absolute(): candidates.append(rp)
                else:
                    candidates += [cf.parent/rp, root/rp]
                target=next((x for x in candidates if x.is_file()),None)
                if not target:
                    rows.append({"checksum_file":str(cf.relative_to(root)),"declared_sha256":h,"declared_path":raw,"resolved_path":"","status":"UNRESOLVED","details":"historical path not found"}); continue
                actual=sha256_file(target)
                rows.append({"checksum_file":str(cf.relative_to(root)),"declared_sha256":h,"declared_path":raw,"resolved_path":str(target),"status":"PASS" if actual==h else "FAIL","details":actual})
    return rows


def generic_rows(path):
    try: return tsv_info(path)[0]
    except Exception: return None


def semantic_metrics(root, expected):
    out=[]
    def add(metric, observed, exp, source, blocker=True, details=""):
        if observed is None: status="CAUTION"
        else: status="PASS" if observed==exp else "FAIL"
        out.append({"metric":metric,"observed":observed if observed is not None else "UNRESOLVED","expected":exp,"status":status,"blocker":"YES" if blocker else "NO","source":source,"details":details})
    # simple row-count authorities
    simple=[
      ("gene_universe_n","01_gene_universe/01_gene_universe.tsv","gene_universe_n"),
      ("mdv2_primary_hits_n","03_intersection/03_intersection_primary_hits.tsv","mdv2_primary_hits_n"),
      ("mdv4_tier1_n","04_bias_empirical/MDV4_empirical/05_tier1_pathways.tsv","mdv4_tier1_n"),
      ("mdv5_higher_order_units_n","05_domains/MDV5_v1_1/MDV5_domain_summary.tsv","mdv5_higher_order_units_n"),
      ("mdv5_standalones_n","05_domains/MDV5_v1_1/MDV5_standalone_pathways.tsv","mdv5_standalones_n"),
      ("pgc_gene_n","06_magma_pgc/MDV6_v1_1/PGC2019_gene_results_w10.tsv","pgc_gene_n"),
      ("spark_gene_n","07_magma_spark/MDV7_v1_1/SPARK_EUR_gene_results_w10.tsv","spark_gene_n"),
      ("mdv8_common_scored_n","08_specificity_controls/MDV8_v1_1/MDV8_common_scored_gene_universe.tsv","mdv8_common_scored_n"),
      ("mdv8_units_n","08_specificity_controls/MDV8_v1_1/MDV8_cross_stage_evidence_matrix.tsv","mdv8_units_n")]
    for m,rel,k in simple:
        p=root/rel; add(m,generic_rows(p),expected[k],rel,True)
    # Gene groups semantic counts using the canonical MDV1 schema.
    # Canonical frozen columns were created by mdv1_final_freeze:
    # IUIS_flag, SFARI_R0_highconf_flag, SFARI_R1_highconf_flag,
    # CoreSeed_flag, R0_group.
    gp=root/"01_gene_universe/02_gene_groups.tsv"
    try:
        with open_text(gp) as f:
            dr=csv.DictReader(f,delimiter="\t"); header=dr.fieldnames or []
            iuis=choose_col(header,["IUIS_flag","IUIS","iuis_gene","is_iuis"])
            sfari=choose_col(header,[
                "SFARI_R0_highconf_flag","SFARI_R0_high_conf_flag",
                "SFARI_flag","SFARI_HC","SFARI_highconf","SFARI_high_conf_flag",
                "sfari_r0_flag","SFARI"
            ])
            core=choose_col(header,["CoreSeed_flag","CoreSeed","coreseed","historical_coreseed"])
            group=choose_col(header,["R0_group","group","gene_group","Group","four_group"])
            missing=[x for x,v in [("IUIS_flag",iuis),("SFARI_R0_highconf_flag",sfari),("CoreSeed_flag",core),("R0_group",group)] if not v]
            if missing:
                raise RuntimeError("canonical gene-group column(s) unresolved: "+",".join(missing)+"; header="+",".join(header))
            counts={"iuis_n":0,"sfari_r0_n":0,"raw_overlap_n":0,"coreseed_n":0,"iei_only_n":0,"asd_only_n":0,"background_n":0}
            n=0; group_mismatch=0; core_outside_overlap=0
            group_counts={"raw_overlap_n":0,"iei_only_n":0,"asd_only_n":0,"background_n":0}
            for r in dr:
                n+=1
                a,b,c=truthy(r.get(iuis,"")),truthy(r.get(sfari,"")),truthy(r.get(core,""))
                expected_group="Overlap_raw" if (a and b) else ("IEI_only" if a else ("ASD_only" if b else "Background"))
                observed_group=str(r.get(group,"")).strip()
                if observed_group != expected_group:
                    group_mismatch += 1
                if c and not (a and b):
                    core_outside_overlap += 1
                counts["iuis_n"]+=int(a); counts["sfari_r0_n"]+=int(b)
                counts["raw_overlap_n"]+=int(a and b); counts["iei_only_n"]+=int(a and not b)
                counts["asd_only_n"]+=int((not a) and b); counts["background_n"]+=int((not a) and (not b))
                counts["coreseed_n"]+=int(c)
                if observed_group=="Overlap_raw": group_counts["raw_overlap_n"]+=1
                elif observed_group=="IEI_only": group_counts["iei_only_n"]+=1
                elif observed_group=="ASD_only": group_counts["asd_only_n"]+=1
                elif observed_group=="Background": group_counts["background_n"]+=1
            details=f"columns: iuis={iuis}, sfari={sfari}, core={core}, group={group}; rows={n}; group_mismatch={group_mismatch}; core_outside_overlap={core_outside_overlap}"
            add("gene_groups_rows_n",n,expected["gene_universe_n"],str(gp.relative_to(root)),True,details)
            add("gene_group_internal_mismatch_n",group_mismatch,0,str(gp.relative_to(root)),True,details)
            add("coreseed_outside_raw_overlap_n",core_outside_overlap,0,str(gp.relative_to(root)),True,details)
            for k in ["iuis_n","sfari_r0_n","raw_overlap_n","coreseed_n","iei_only_n","asd_only_n","background_n"]:
                cross_ok=True
                if k in group_counts:
                    cross_ok=(counts[k]==group_counts[k])
                add(k,counts[k],expected[k],str(gp.relative_to(root)),True,details+f"; R0_group_crosscheck={cross_ok}")
    except Exception as e:
        add("gene_groups_rows_n",None,expected["gene_universe_n"],str(gp.relative_to(root)),True,str(e))
        add("gene_group_internal_mismatch_n",None,0,str(gp.relative_to(root)),True,str(e))
        add("coreseed_outside_raw_overlap_n",None,0,str(gp.relative_to(root)),True,str(e))
        for k in ["iuis_n","sfari_r0_n","raw_overlap_n","coreseed_n","iei_only_n","asd_only_n","background_n"]:
            add(k,None,expected[k],str(gp.relative_to(root)),True,str(e))

    # Pathway counts: pathway_manifest contains BOTH raw and primary-filtered
    # pathways. Only primary_10_500_flag==1 belongs to the frozen 6,671 family.
    pm=root/"02_pathways/pathway_manifest.tsv"
    mem=root/"02_pathways/03_pathway_membership.tsv.gz"
    try:
        with open_text(pm) as f:
            dr=csv.DictReader(f,delimiter="\t"); header=dr.fieldnames or []
            source=choose_col(header,["source","database","db","pathway_source"])
            pid=choose_col(header,["pathway_id","term_id","id","pathway","set_id"])
            primary=choose_col(header,["primary_10_500_flag","primary_flag","filtered_flag","retain_flag"])
            if not source or not pid or not primary:
                raise RuntimeError(f"canonical pathway-manifest column(s) unresolved: pid={pid}; source={source}; primary={primary}; header={header}")
            raw_ids=set(); primary_ids=set(); go=set(); reacs=set()
            for r in dr:
                ident=(r.get(source,"").strip(), r.get(pid,"").strip())
                raw_ids.add(ident)
                if truthy(r.get(primary,"")):
                    primary_ids.add(ident)
                    s=norm(r.get(source,""))
                    if s in {"gobp","go"} or "gobp" in s: go.add(ident)
                    if "reactome" in s: reacs.add(ident)

        # Independent cross-check against the canonical filtered membership file.
        mem_ids=set()
        with open_text(mem) as f:
            dr=csv.DictReader(f,delimiter="\t"); h=dr.fieldnames or []
            msource=choose_col(h,["source","database","db","pathway_source"])
            mpid=choose_col(h,["pathway_id","term_id","id","pathway","set_id"])
            if not msource or not mpid:
                raise RuntimeError(f"filtered membership schema unresolved: source={msource}; pid={mpid}; header={h}")
            for r in dr:
                mem_ids.add((r.get(msource,"").strip(), r.get(mpid,"").strip()))
        symdiff=len(primary_ids.symmetric_difference(mem_ids))
        details=f"pid={pid}; source={source}; primary={primary}; raw_terms={len(raw_ids)}; primary_terms={len(primary_ids)}; membership_terms={len(mem_ids)}; primary_vs_membership_symdiff={symdiff}"
        add("pathway_primary_set_symmetric_diff_n",symdiff,0,str(pm.relative_to(root)),True,details)
        add("pathways_total_n",len(primary_ids),expected["pathways_total_n"],str(pm.relative_to(root)),True,details)
        add("go_bp_n",len(go),expected["go_bp_n"],str(pm.relative_to(root)),True,details)
        add("reactome_n",len(reacs),expected["reactome_n"],str(pm.relative_to(root)),True,details)
    except Exception as e:
        add("pathway_primary_set_symmetric_diff_n",None,0,str(pm.relative_to(root)),True,str(e))
        add("pathways_total_n",None,expected["pathways_total_n"],str(pm.relative_to(root)),True,str(e))
        add("go_bp_n",None,expected["go_bp_n"],str(pm.relative_to(root)),True,str(e))
        add("reactome_n",None,expected["reactome_n"],str(pm.relative_to(root)),True,str(e))
    return out


def analysis(root, ma, cfgdir):
    out=ma/"00_contract"
    lockp=out/"MAR00_preanalysis_lock.json"; release=out/"MA_R00_STOPA_RELEASE.txt"
    if not lockp.exists() or not release.exists(): print("FATAL: STOP A release missing",file=sys.stderr); return 2
    locksha=sha256_file(lockp)
    txt=release.read_text(errors="replace")
    if "decision=PASS" not in txt or f"preanalysis_lock_sha256={locksha}" not in txt:
        print("FATAL: release token invalid or bound to different lock",file=sys.stderr); return 2
    # v1.0.1 is a schema-only semantic-QC amendment. It requires a second
    # human release token bound to BOTH the original preanalysis lock and the
    # patched audit-script SHA256.
    if PKG_VERSION=="1.0.1":
        aready=out/"MAR00_SCHEMA_AMENDMENT_READY.txt"
        arel=out/"MA_R00_SCHEMA_AMENDMENT_RELEASE.txt"
        if not aready.exists() or not arel.exists():
            print("FATAL: schema amendment release missing",file=sys.stderr); return 2
        current_script_sha=sha256_file(Path(__file__))
        atxt=arel.read_text(errors="replace")
        if ("decision=PASS" not in atxt or
            f"preanalysis_lock_sha256={locksha}" not in atxt or
            f"patched_audit_sha256={current_script_sha}" not in atxt):
            print("FATAL: schema amendment release invalid or bound to different lock/script",file=sys.stderr); return 2
    manifest=read_manifest(cfgdir/"ma_r00_authority_manifest.tsv")
    authority=[]
    for row in manifest:
        p=root/row["relpath"]
        if p.is_file():
            st=p.stat(); status="PASS" if st.st_size>0 else "FAIL"
            authority.append({"role":row["role"],"path":str(p),"required":row["required"],"policy":"READ_ONLY","filesystem_writable":"YES" if os.access(p,os.W_OK) else "NO","size_bytes":st.st_size,"mtime_ns":st.st_mtime_ns,"sha256":sha256_file(p),"status":status})
        else:
            authority.append({"role":row["role"],"path":str(p),"required":row["required"],"policy":"READ_ONLY","filesystem_writable":"NA","size_bytes":0,"mtime_ns":"","sha256":"MISSING","status":"FAIL" if row["required"]=="1" else "CAUTION"})
    write_tsv(out/"MAR00_legacy_authority_manifest.tsv",authority,["role","path","required","policy","filesystem_writable","size_bytes","mtime_ns","sha256","status"])
    runtimes=[command_version(c) for c in ["python3","Rscript","bash","git","tar","sha256sum","magma"]]
    runtimes.append({"tool":"platform","available":"YES","path":"","version":platform.platform()})
    write_tsv(out/"MAR00_runtime_versions.tsv",runtimes,["tool","available","path","version"])
    guard=[]
    legacy_dirs=["00_provenance","01_gene_universe","02_pathways","03_intersection","04_bias_empirical","05_domains","06_magma_pgc","07_magma_spark","08_specificity_controls","09_figures_tables","10_reproducibility"]
    for d in legacy_dirs:
        p=root/d; guard.append({"path_prefix":str(p),"policy":"READ_ONLY","filesystem_writable":"YES" if os.access(p,os.W_OK) else "NO","allowed_ma_write":"NO","audit_result":"PASS"})
    guard.append({"path_prefix":str(ma),"policy":"MA_OUTPUT","filesystem_writable":"YES" if os.access(ma,os.W_OK) else "NO","allowed_ma_write":"YES","audit_result":"PASS" if os.access(ma,os.W_OK) else "FAIL"})
    write_tsv(out/"MAR00_path_write_guard.tsv",guard,["path_prefix","policy","filesystem_writable","allowed_ma_write","audit_result"])
    cs=checksum_audit(root); write_tsv(out/"MAR00_checksum_audit.tsv",cs,["checksum_file","declared_sha256","declared_path","resolved_path","status","details"])
    expected=json.load(open(cfgdir/"ma_r00_expected_metrics.json"))
    metrics=semantic_metrics(root,expected); write_tsv(out/"MAR00_semantic_metrics.tsv",metrics,["metric","observed","expected","status","blocker","source","details"])
    stage=[]
    for s,p in stage_paths(root).items():
        stage.append({"stage":s,"pass_file":str(p),"exists":"YES" if p.exists() else "NO","content_head":p.read_text(errors="replace")[:300].replace("\n"," | ") if p.exists() else "","status":"PASS" if p.exists() else "FAIL"})
    write_tsv(out/"MAR00_legacy_stage_status.tsv",stage,["stage","pass_file","exists","content_head","status"])
    # Build project lock
    blockers=[m for m in metrics if m["blocker"]=="YES" and m["status"]!="PASS"]
    critical_auth=[a for a in authority if a["required"]=="1" and a["status"]!="PASS"]
    checksum_fail=[x for x in cs if x["status"]=="FAIL"]
    project_lock={"stage":"MA-R00","timestamp_utc":utcnow(),"package_version":PKG_VERSION,"legacy_authority_manifest_sha256":sha256_file(out/"MAR00_legacy_authority_manifest.tsv"),"semantic_metrics_sha256":sha256_file(out/"MAR00_semantic_metrics.tsv"),"preanalysis_lock_sha256":locksha,"critical_authority_failures":len(critical_auth),"semantic_blocker_mismatches":len(blockers),"checksum_failures":len(checksum_fail),"legacy_policy":"READ_ONLY","ma_write_root":str(ma.resolve()),"status":"PASS_CANDIDATE" if not critical_auth and not blockers and not checksum_fail else "HOLD"}
    with open(out/"MA_REVISION_PROJECT_LOCK.json","w") as f: json.dump(project_lock,f,indent=2)
    (out/"MAR00_ANALYSIS_COMPLETE.txt").write_text(f"{project_lock['status']}\n{utcnow()}\n")
    return 0 if project_lock["status"]=="PASS_CANDIDATE" else 3


def postflight(root, ma, cfgdir):
    out=ma/"00_contract"
    required=["MAR00_resolved_inputs.json","MAR00_preanalysis_lock.json","MA_R00_STOPA_RELEASE.txt","MAR00_legacy_authority_manifest.tsv","MAR00_runtime_versions.tsv","MAR00_path_write_guard.tsv","MAR00_checksum_audit.tsv","MAR00_semantic_metrics.tsv","MAR00_legacy_stage_status.tsv","MA_REVISION_PROJECT_LOCK.json","MAR00_ANALYSIS_COMPLETE.txt"]
    if PKG_VERSION=="1.0.1":
        required += ["MAR00_SCHEMA_AMENDMENT_v1.0.1.json","MAR00_SCHEMA_AMENDMENT_READY.txt","MA_R00_SCHEMA_AMENDMENT_RELEASE.txt"]
    checks=[]
    for name in required:
        p=out/name; checks.append({"check_id":f"OUTPUT::{name}","status":"PASS" if p.is_file() and p.stat().st_size>0 else "FAIL","details":str(p)})
    # verify legacy baseline unchanged
    base=list(csv.DictReader(open(out/"MAR00_legacy_hash_baseline.tsv"),delimiter="\t"))
    changed=0
    for r in base:
        p=Path(r["path"])
        if not p.is_file(): status="FAIL"; details="MISSING_AFTER_FREEZE"; changed+=1
        else:
            st=p.stat(); actual=sha256_file(p)
            same=(str(st.st_size)==str(r["size_bytes"]) and str(st.st_mtime_ns)==str(r["mtime_ns"]) and actual==r["sha256"])
            status="PASS" if same else "FAIL"; details=f"size={st.st_size};mtime_ns={st.st_mtime_ns};sha256={actual}"
            if not same: changed+=1
        checks.append({"check_id":f"LEGACY_UNCHANGED::{r['role']}","status":status,"details":details})
    # project lock status
    pl=json.load(open(out/"MA_REVISION_PROJECT_LOCK.json"))
    checks.append({"check_id":"PROJECT_LOCK_STATUS","status":"PASS" if pl.get("status")=="PASS_CANDIDATE" else "FAIL","details":pl.get("status","")})
    write_tsv(out/"MAR00_postflight_checks.tsv",checks,["check_id","status","details"])
    fails=[x for x in checks if x["status"]!="PASS"]
    # checksums for MA-R00 outputs, excluding checksum itself and bundles
    checksum_lines=[]
    for p in sorted(out.glob("MAR00_*")) + sorted(out.glob("MA_REVISION_PROJECT_LOCK.json")) + sorted(out.glob("MA_R00_STOPA_RELEASE.txt")) + sorted(out.glob("MA_R00_SCHEMA_AMENDMENT_RELEASE.txt")):
        if p.is_file() and p.name not in {"MAR00_checksums.sha256","MAR00_PASS.txt","MAR00_HOLD.txt","MAR00_validation_summary.txt"}:
            checksum_lines.append(f"{sha256_file(p)}  {p.name}")
    (out/"MAR00_checksums.sha256").write_text("\n".join(checksum_lines)+"\n")
    status="PASS" if not fails else "HOLD"
    summary=[f"MA-R00 validation summary",f"timestamp_utc={utcnow()}",f"status={status}",f"postflight_failures={len(fails)}",f"legacy_changed={changed}",f"project_lock_status={pl.get('status')}","", "NEXT STEP: If PASS, build/upload final R00 review bundle before starting raw26 WP01."]
    (out/"MAR00_validation_summary.txt").write_text("\n".join(summary)+"\n")
    marker=out/("MAR00_PASS.txt" if status=="PASS" else "MAR00_HOLD.txt")
    marker.write_text(f"{status}\n{utcnow()}\n")
    return 0 if status=="PASS" else 4


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("phase",choices=["preflight","freeze","analysis","postflight"])
    ap.add_argument("--root",default=DEFAULT_ROOT)
    ap.add_argument("--ma-root",default=None)
    ap.add_argument("--config-dir",required=True)
    a=ap.parse_args()
    root=Path(a.root).expanduser().resolve(); ma=Path(a.ma_root).expanduser().resolve() if a.ma_root else root/"12_molecular_autism_revision"; cfg=Path(a.config_dir).resolve()
    if a.phase=="preflight": rc=preflight(root,ma,cfg)
    elif a.phase=="freeze": rc=freeze(root,ma,cfg)
    elif a.phase=="analysis": rc=analysis(root,ma,cfg)
    else: rc=postflight(root,ma,cfg)
    sys.exit(rc)

if __name__=="__main__": main()

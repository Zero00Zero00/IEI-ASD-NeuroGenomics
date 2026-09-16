
from common import *
import tarfile

if not (OUT/"MAR03_PASS.txt").is_file():
    raise SystemExit("HOLD: MAR03_PASS missing")
ts=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
handoff=MA/"90_handoff"; handoff.mkdir(parents=True,exist_ok=True)
out=handoff/f"MAR03_REVIEW_{ts}.tar.gz"

include=[
 "MAR03_preflight_report.tsv","MAR03_resolved_inputs.json","MAR03_HISTORICAL_INTERFACE_SNAPSHOT.tsv",
 "MAR03_preanalysis_lock.json","MAR03_METHOD_IMPLEMENTATION_NOTE.txt","MAR03_STOPA_RELEASE.txt",
 "MAR03_frozen_hypotheses.tsv","MAR03_preanalysis_coverage.tsv","MAR03_pgc_higher_order.tsv",
 "MAR03_spark_higher_order.tsv","MAR03_standalone_tests.tsv","MAR03_gene_set_results_all.tsv",
 "MAR03_matched_meanZ.tsv","MAR03_matching_qc.tsv","MAR03_matching_summary.tsv",
 "MAR03_sensitivity.tsv","MAR03_primary_classification.tsv","MAR03_postflight_checks.tsv",
 "MAR03_analysis_lock.json","MAR03_validation_summary.txt","MAR03_checksums.sha256","MAR03_PASS.txt"
]
with tarfile.open(out,"w:gz") as tf:
    for name in include:
        p=OUT/name
        if p.exists(): tf.add(p,arcname=name)
    for p in (OUT/"source_data").glob("*"):
        tf.add(p,arcname="source_data/"+p.name)
sha=sha256_file(out)
Path(str(out)+".sha256").write_text(f"{sha}  {out.name}\n",encoding="utf-8")
print(f"REVIEW_PACKAGE={out}")
print(f"SHA256={sha}")
print("PACK_REVIEW=PASS")

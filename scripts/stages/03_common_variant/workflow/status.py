
from common import *
print(f"MAR03_OUT={OUT}")
for f in ["MAR03_preflight_report.tsv","MAR03_preanalysis_lock.json","MAR03_STOPA_RELEASE.txt",
          "MAR03_gene_set_results_all.tsv","MAR03_matched_meanZ.tsv","MAR03_sensitivity.tsv",
          "MAR03_postflight_checks.tsv","MAR03_PASS.txt"]:
    print(("OK  " if (OUT/f).is_file() else "MISS"),OUT/f)

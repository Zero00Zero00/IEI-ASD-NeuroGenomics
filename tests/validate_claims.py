from pathlib import Path
import json
import math
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SD = ROOT / "source_data"
EXPECTED = json.loads((ROOT / "tests" / "expected_values.json").read_text())

def read_tsv(path):
    path = ROOT / path
    if not path.exists():
        raise AssertionError(f"Missing required source table: {path}")
    return pd.read_csv(path, sep="\t")

def choose_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise AssertionError(
        f"None of {candidates} found. Available columns: {list(df.columns)}"
    )

def assert_eq(name, observed, expected):
    if observed != expected:
        raise AssertionError(f"{name}: observed={observed}, expected={expected}")
    print(f"PASS {name}: {observed}")

def assert_close(name, observed, expected, tol=1e-6):
    if not math.isclose(float(observed), float(expected), rel_tol=tol, abs_tol=tol):
        raise AssertionError(f"{name}: observed={observed}, expected={expected}")
    print(f"PASS {name}: {observed}")

# 26-gene catalog intersection
raw26 = read_tsv(
    "source_data/manuscript/WP01/WP01_raw26_gene_set.tsv"
)
gene_col = choose_col(
    raw26,
    ["gene", "symbol", "gene_symbol", "HGNC_symbol", "hgnc_symbol"]
)
assert_eq(
    "catalog_intersection_genes",
    int(raw26[gene_col].dropna().nunique()),
    EXPECTED["catalog_intersection_genes"]
)

# 18 dual-background pathway candidates
hits = read_tsv(
    "source_data/derived_inputs/03_intersection_primary_hits.tsv"
)
assert_eq(
    "dual_background_candidates",
    len(hits),
    EXPECTED["dual_background_candidates"]
)

# 39 BH-prioritized pathways
tier = read_tsv(
    "source_data/manuscript/WP01/WP01_raw26_tier1.tsv"
)
assert_eq(
    "bh_prioritized_pathways",
    len(tier),
    EXPECTED["bh_prioritized_pathways"]
)

# Redundancy-compressed hypothesis space
reps = read_tsv(
    "source_data/manuscript/MAR02/MAR02_SOURCE_REPRESENTATIVES.tsv"
)
units = read_tsv(
    "source_data/manuscript/MAR02/MAR02_SOURCE_UNIT_SUMMARY.tsv"
)
stand = read_tsv(
    "source_data/manuscript/MAR02/MAR02_SOURCE_STANDALONES.tsv"
)

assert_eq("representatives", len(reps), EXPECTED["representatives"])
assert_eq("higher_order_units", len(units), EXPECTED["higher_order_units"])
assert_eq("standalones", len(stand), EXPECTED["standalones"])

# Common-variant primary family
primary = read_tsv(
    "source_data/manuscript/MAR03/MAR03_SOURCE_PRIMARY_HIGHER_ORDER.tsv"
)
qcol = choose_col(primary, ["q_BH", "q", "BH_q"])
dataset_col = choose_col(primary, ["dataset", "Dataset"])

for dataset, key in [
    ("PGC", "pgc_primary_positive"),
    ("SPARK", "spark_primary_positive")
]:
    d = primary[primary[dataset_col].astype(str).str.upper() == dataset]
    assert_eq(
        key,
        int((pd.to_numeric(d[qcol], errors="coerce") < 0.05).sum()),
        EXPECTED[key]
    )

# Secondary standalone family
stand_res = read_tsv(
    "source_data/manuscript/MAR03/MAR03_SOURCE_STANDALONE_RESULTS.tsv"
)
sqcol = choose_col(stand_res, ["q_BH", "q", "BH_q"])
sdataset = choose_col(stand_res, ["dataset", "Dataset"])

for dataset, key in [
    ("PGC", "pgc_standalone_positive"),
    ("SPARK", "spark_standalone_positive")
]:
    d = stand_res[
        stand_res[sdataset].astype(str).str.upper() == dataset
    ]
    assert_eq(
        key,
        int((pd.to_numeric(d[sqcol], errors="coerce") < 0.05).sum()),
        EXPECTED[key]
    )

# Dependency-robust sensitivity: flexible schema handling
stability = read_tsv(
    "source_data/manuscript/MAR07/MAR07_SOURCE_TIER1_STABILITY.tsv"
)

by_count = None

for c in stability.columns:
    lc = c.lower()
    if "by" in lc and ("retained" in lc or "tier1" in lc or "priority" in lc):
        s = stability[c]
        if pd.api.types.is_numeric_dtype(s):
            by_count = int((pd.to_numeric(s, errors="coerce") > 0).sum())
            break
        vals = s.astype(str).str.upper()
        by_count = int(vals.isin(["YES", "TRUE", "TIER1", "PRIORITY", "RETAINED"]).sum())
        break

if by_count is None:
    by_qcols = [c for c in stability.columns if "by" in c.lower() and "q" in c.lower()]
    if by_qcols:
        by_count = int(
            (pd.to_numeric(stability[by_qcols[0]], errors="coerce") < 0.05).sum()
        )

if by_count is None:
    raise AssertionError(
        "Could not infer BY-retention column from "
        f"MAR07_SOURCE_TIER1_STABILITY.tsv. Columns={list(stability.columns)}"
    )

assert_eq(
    "by_retained_priorities",
    by_count,
    EXPECTED["by_retained_priorities"]
)

# Full-chain global-null calibration
metrics = read_tsv(
    "source_data/manuscript/MAR07/MAR07_SOURCE_CHAIN_CALIBRATION_METRICS.tsv"
)

rcol = choose_col(metrics, ["R", "n_replicates", "replicate_n"])
ratecol = choose_col(
    metrics,
    ["any_false_tier1_rate", "any_false_priority_rate", "p_any_false_tier1"]
)
meancol = choose_col(
    metrics,
    ["mean_false_tier1_n", "mean_false_priority_n", "mean_tier1_n"]
)

final = metrics.loc[pd.to_numeric(metrics[rcol], errors="coerce").idxmax()]

assert_eq(
    "global_null_replicates",
    int(final[rcol]),
    EXPECTED["global_null_replicates"]
)

assert_close(
    "global_null_any_false_rate",
    float(final[ratecol]),
    EXPECTED["global_null_any_false_rate"]
)

assert_close(
    "global_null_mean_false_count",
    float(final[meancol]),
    EXPECTED["global_null_mean_false_count"],
    tol=1e-4
)

print("\nALL MANUSCRIPT CLAIM CHECKS PASSED")

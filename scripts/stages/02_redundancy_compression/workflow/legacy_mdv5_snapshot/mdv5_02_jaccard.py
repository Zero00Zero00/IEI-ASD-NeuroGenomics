#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import sys

import numpy as np
import pandas as pd

from mdv5_common import (
    MDV5Error, as_bool, canonicalize_membership, canonicalize_pathway_columns, canonicalize_universe,
    ensure_numeric, load_config, outdir, output_path, read_table, resolve, root_from_cfg, sentinel_passes,
    write_pass, write_table,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    root = root_from_cfg(cfg)
    outdir(cfg, root)
    if not sentinel_passes(output_path(cfg, root, "preflight_pass")):
        raise MDV5Error("MDV5_PREFLIGHT_PASS missing/invalid")

    tier1 = canonicalize_pathway_columns(read_table(resolve(root, cfg["upstream"]["mdv4_tier1"])))
    tier1 = ensure_numeric(tier1, ["OR_Firth", "or_firth", "OR"], "OR_Firth")
    tier1 = ensure_numeric(tier1, ["q_Firth", "q_firth", "Firth_q"], "q_Firth")
    tier1 = ensure_numeric(tier1, ["P_emp", "p_emp", "empirical_p"], "P_emp", required=False)
    tier1 = ensure_numeric(tier1, ["q_emp", "qEmp", "empirical_q"], "q_emp")
    tier1 = tier1.sort_values(["source", "pathway_id"], kind="mergesort").reset_index(drop=True)

    membership = canonicalize_membership(read_table(resolve(root, cfg["upstream"]["pathway_membership"])))
    universe = canonicalize_universe(read_table(resolve(root, cfg["upstream"]["universe"])))
    core = universe.loc[as_bool(universe["CoreSeed_flag"]), ["HGNC_id", "HGNC_symbol"]].copy()
    core_ids = set(core["HGNC_id"])
    core_sym = dict(zip(core["HGNC_id"], core["HGNC_symbol"]))

    mem_by_key = {k: set(g["HGNC_id"].astype(str)) for k, g in membership.groupby("pathway_key", sort=False)}
    keys = tier1["pathway_key"].tolist()
    if len(keys) != cfg["expected"]["tier1_n"]:
        raise MDV5Error(f"Tier-1 drift: {len(keys)} != {cfg['expected']['tier1_n']}")

    # Evidence annotation is descriptive only; it never enters clustering weights or representative ranking.
    sizes = []
    driver_ns = []
    driver_strs = []
    for key in keys:
        genes = mem_by_key.get(key)
        if genes is None:
            raise MDV5Error(f"Tier-1 pathway absent from full membership: {key}")
        drv = sorted(genes & core_ids, key=lambda x: core_sym.get(x, x))
        sizes.append(len(genes))
        driver_ns.append(len(drv))
        driver_strs.append(";".join(core_sym.get(x, x) for x in drv))
    tier1["pathway_size"] = sizes
    tier1["CoreSeed_driver_n"] = driver_ns
    tier1["CoreSeed_drivers"] = driver_strs
    tier1["single_driver_flag"] = (tier1["CoreSeed_driver_n"] == 1).astype(int)
    tier1["multi_driver_flag"] = (tier1["CoreSeed_driver_n"] >= 2).astype(int)
    write_table(tier1, output_path(cfg, root, "tier1_annotated"))

    n = len(keys)
    mat = np.eye(n, dtype=float)
    pair_rows = []
    for i, j in itertools.combinations(range(n), 2):
        a, b = mem_by_key[keys[i]], mem_by_key[keys[j]]
        inter = len(a & b)
        union = len(a | b)
        jac = inter / union if union else 0.0
        mat[i, j] = mat[j, i] = jac
        pair_rows.append({
            "pathway_key_1": keys[i], "pathway_key_2": keys[j],
            "n1": len(a), "n2": len(b), "intersection_n": inter, "union_n": union,
            "jaccard": jac,
        })

    matrix_df = pd.DataFrame(mat, columns=keys)
    matrix_df.insert(0, "pathway_key", keys)
    write_table(matrix_df, output_path(cfg, root, "jaccard_matrix"))
    write_table(pd.DataFrame(pair_rows), output_path(cfg, root, "jaccard_pairs"))

    tol = float(cfg["jaccard"].get("pairwise_tolerance", 1e-12))
    if not np.allclose(mat, mat.T, atol=tol, rtol=0):
        raise MDV5Error("Jaccard matrix is not symmetric")
    if not np.allclose(np.diag(mat), 1.0, atol=tol, rtol=0):
        raise MDV5Error("Jaccard diagonal is not exactly 1 within tolerance")
    if np.nanmin(mat) < -tol or np.nanmax(mat) > 1 + tol:
        raise MDV5Error("Jaccard values outside [0,1]")

    write_pass(output_path(cfg, root, "jaccard_pass"), {
        "status": "PASS",
        "technical_status": "PASS",
        "stage": "MDV5_FULL_MEMBERSHIP_JACCARD",
        "tier1_n": n,
        "pair_n": len(pair_rows),
        "jaccard_definition": "full_frozen_membership_intersection_over_union",
        "gwas_inputs_used": "NO",
    })
    print(f"[PASS] Full-membership Jaccard: n={n}, pairs={len(pair_rows)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MDV5Error as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        raise SystemExit(2)

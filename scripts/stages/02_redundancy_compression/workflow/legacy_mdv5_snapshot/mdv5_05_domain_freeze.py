#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import sys
from collections import defaultdict

import pandas as pd

from mdv5_common import (
    MDV5Error, as_bool, atomic_text, canonicalize_membership, canonicalize_universe, load_config,
    outdir, output_path, read_table, resolve, root_from_cfg, sentinel_passes, write_table,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    root = root_from_cfg(cfg)
    outdir(cfg, root)
    if not sentinel_passes(output_path(cfg, root, "leiden_pass")):
        raise MDV5Error("MDV5_LEIDEN_PASS missing/invalid")

    rmap = read_table(output_path(cfg, root, "redundancy_map"))
    reps = read_table(output_path(cfg, root, "representatives"))
    partition = read_table(output_path(cfg, root, "primary_partition"))
    membership = canonicalize_membership(read_table(resolve(root, cfg["upstream"]["pathway_membership"])))
    universe = canonicalize_universe(read_table(resolve(root, cfg["upstream"]["universe"])))
    core_ids = set(universe.loc[as_bool(universe["CoreSeed_flag"]), "HGNC_id"].astype(str))
    universe_n = len(universe)
    sym_map = dict(zip(universe["HGNC_id"].astype(str), universe["HGNC_symbol"].astype(str)))

    part_map = partition.set_index("pathway_key").to_dict("index")
    rep_to_family = reps.set_index("pathway_key")["redundancy_family_id"].to_dict()
    family_to_rep = rmap.groupby("redundancy_family_id")["representative_pathway_key"].first().to_dict()
    original_to_rep = rmap.set_index("pathway_key")["representative_pathway_key"].to_dict()
    mem_by_key = {k: set(g["HGNC_id"].astype(str)) for k, g in membership.groupby("pathway_key", sort=False)}

    # Map every original Tier-1 pathway exactly once through its representative to a primary frozen unit or non-primary status.
    ptu_rows = []
    for _, row in rmap.sort_values(["redundancy_family_id", "representative_rank"], kind="mergesort").iterrows():
        orig = str(row["pathway_key"])
        rep = str(row["representative_pathway_key"])
        p = part_map[rep]
        ptu_rows.append({
            "pathway_key": orig,
            "source": row["source"],
            "pathway_id": row["pathway_id"],
            "pathway_name": row["pathway_name"],
            "redundancy_family_id": row["redundancy_family_id"],
            "representative_pathway_key": rep,
            "representative_flag": row["representative_flag"],
            "medoid_cluster": p["medoid_cluster"],
            "community_size": p["community_size"],
            "median_within_community_coclustering": p["median_within_community_coclustering"],
            "unit_type": p["unit_type"],
            "unit_id": p["unit_id"] if pd.notna(p["unit_id"]) else "",
        })
    ptu = pd.DataFrame(ptu_rows)
    write_table(ptu, output_path(cfg, root, "pathway_to_unit"))

    primary_units = partition.loc[partition["unit_type"].isin(["stable_domain", "small_module"])].copy()
    unit_ids = sorted([x for x in primary_units["unit_id"].dropna().astype(str).unique() if x])
    stamp = dt.datetime.now(dt.timezone.utc).isoformat()
    label_prefix = str(cfg["domain_gene_sets"]["placeholder_label_prefix"])
    broad_n_thr = int(cfg["domain_gene_sets"]["oversized_expanded_gene_n"])
    broad_frac_thr = float(cfg["domain_gene_sets"]["oversized_expanded_universe_fraction"])

    membership_rows = []
    summary_rows = []
    locus_rows = []
    for uid in unit_ids:
        unit_reps = sorted(primary_units.loc[primary_units["unit_id"] == uid, "pathway_key"].astype(str).tolist())
        unit_type = str(primary_units.loc[primary_units["unit_id"] == uid, "unit_type"].iloc[0])
        stability = float(primary_units.loc[primary_units["unit_id"] == uid, "median_within_community_coclustering"].iloc[0])
        original_keys = sorted(ptu.loc[ptu["unit_id"].astype(str) == uid, "pathway_key"].astype(str).tolist())
        rep_counts = defaultdict(int)
        all_counts = defaultdict(int)
        for k in unit_reps:
            for g in mem_by_key[k]:
                rep_counts[g] += 1
        for k in original_keys:
            for g in mem_by_key[k]:
                all_counts[g] += 1
        expanded = set(rep_counts)
        compact = {g for g, c in rep_counts.items() if c >= 2}
        all_tier1 = set(all_counts)
        domain_core = all_tier1 & core_ids
        row_genes = sorted(all_tier1, key=lambda g: sym_map.get(g, g))
        broad = int(len(expanded) > broad_n_thr or (len(expanded) / universe_n) > broad_frac_thr)
        label = f"{label_prefix}_{uid}"
        for g in row_genes:
            membership_rows.append({
                "domain_id": uid,
                "domain_label": label,
                "unit_type": unit_type,
                "representative_pathway_ids": ";".join(unit_reps),
                "original_tier1_pathway_ids": ";".join(original_keys),
                "HGNC_id": g,
                "gene": sym_map.get(g, g),
                "CoreSeed_flag": int(g in core_ids),
                "membership_count": int(rep_counts.get(g, 0)),
                "membership_count_all_tier1": int(all_counts.get(g, 0)),
                "DomainCore_flag": int(g in domain_core),
                "Expanded_flag": int(g in expanded),
                "Compact_flag": int(g in compact),
                "AllTier1_flag": int(g in all_tier1),
                "broad_domain_flag": broad,
                "freeze_timestamp": stamp,
            })
        # Locus-dominance audit uses CoreSeed driver provenance only; it never changes clustering or membership.
        unit_rmap = rmap.loc[rmap["pathway_key"].astype(str).isin(original_keys)].copy()
        driver_counts = defaultdict(int)
        driver_union = set()
        for x in unit_rmap["CoreSeed_drivers"].fillna("").astype(str):
            genes = [g.strip() for g in x.split(";") if g.strip()]
            driver_union.update(genes)
            for g in genes:
                driver_counts[g] += 1
        dominant_gene = ""
        dominant_n = 0
        if driver_counts:
            dominant_gene, dominant_n = sorted(driver_counts.items(), key=lambda kv: (-kv[1], kv[0]))[0]
        original_n = max(1, len(original_keys))
        dominant_fraction = dominant_n / original_n
        single_n = int((pd.to_numeric(unit_rmap["CoreSeed_driver_n"], errors="coerce") == 1).sum())
        multi_n = int((pd.to_numeric(unit_rmap["CoreSeed_driver_n"], errors="coerce") >= 2).sum())
        locus_flag = int(dominant_fraction >= float(cfg["domain_gene_sets"].get("locus_dominance_fraction_flag", 0.50)))
        locus_rows.append({
            "domain_id": uid, "unit_type": unit_type,
            "original_tier1_pathway_n": len(original_keys),
            "unique_CoreSeed_driver_n": len(driver_union),
            "CoreSeed_drivers": ";".join(sorted(driver_union)),
            "single_driver_pathway_n": single_n,
            "single_driver_pathway_fraction": single_n / original_n,
            "multi_driver_pathway_n": multi_n,
            "multi_driver_pathway_fraction": multi_n / original_n,
            "dominant_driver_gene": dominant_gene,
            "dominant_driver_pathway_n": dominant_n,
            "dominant_driver_fraction": dominant_fraction,
            "locus_dominated_flag": locus_flag,
            "locus_dominance_flag_threshold": float(cfg["domain_gene_sets"].get("locus_dominance_fraction_flag", 0.50)),
            "inference_role": "descriptive_only_not_used_for_clustering_or_unit_selection",
        })
        summary_rows.append({
            "domain_id": uid,
            "domain_label": label,
            "unit_type": unit_type,
            "representative_pathway_n": len(unit_reps),
            "original_tier1_pathway_n": len(original_keys),
            "median_coclustering": stability,
            "DomainCore_gene_n": len(domain_core),
            "Expanded_gene_n": len(expanded),
            "Compact_gene_n": len(compact),
            "AllTier1_gene_n": len(all_tier1),
            "Expanded_universe_fraction": len(expanded) / universe_n,
            "broad_domain_flag": broad,
            "unique_CoreSeed_driver_n": len(driver_union),
            "dominant_driver_gene": dominant_gene,
            "dominant_driver_fraction": dominant_fraction,
            "locus_dominated_flag": locus_flag,
            "freeze_timestamp": stamp,
        })

    mem_df = pd.DataFrame(membership_rows)
    if mem_df.empty:
        mem_df = pd.DataFrame(columns=[
            "domain_id","domain_label","unit_type","representative_pathway_ids","original_tier1_pathway_ids",
            "HGNC_id","gene","CoreSeed_flag","membership_count","membership_count_all_tier1","DomainCore_flag",
            "Expanded_flag","Compact_flag","AllTier1_flag","broad_domain_flag","freeze_timestamp"
        ])
    summary_df = pd.DataFrame(summary_rows)
    if summary_df.empty:
        summary_df = pd.DataFrame(columns=[
            "domain_id","domain_label","unit_type","representative_pathway_n","original_tier1_pathway_n",
            "median_coclustering","DomainCore_gene_n","Expanded_gene_n","Compact_gene_n","AllTier1_gene_n",
            "Expanded_universe_fraction","broad_domain_flag","unique_CoreSeed_driver_n","dominant_driver_gene",
            "dominant_driver_fraction","locus_dominated_flag","freeze_timestamp"
        ])
    locus_df = pd.DataFrame(locus_rows)
    if locus_df.empty:
        locus_df = pd.DataFrame(columns=[
            "domain_id","unit_type","original_tier1_pathway_n","unique_CoreSeed_driver_n","CoreSeed_drivers",
            "single_driver_pathway_n","single_driver_pathway_fraction","multi_driver_pathway_n",
            "multi_driver_pathway_fraction","dominant_driver_gene","dominant_driver_pathway_n",
            "dominant_driver_fraction","locus_dominated_flag","locus_dominance_flag_threshold","inference_role"
        ])
    write_table(mem_df, output_path(cfg, root, "domain_membership"))
    write_table(summary_df, output_path(cfg, root, "domain_summary"))
    write_table(locus_df, output_path(cfg, root, "locus_dominance_audit"))

    standalone = partition.loc[~partition["unit_type"].isin(["stable_domain", "small_module"])].merge(
        reps[["pathway_key", "source", "pathway_id", "pathway_name", "q_emp", "q_Firth", "pathway_size", "CoreSeed_driver_n", "CoreSeed_drivers"]],
        on="pathway_key", how="left"
    )
    write_table(standalone, output_path(cfg, root, "standalone_pathways"))

    stable_n = int(summary_df.loc[summary_df["unit_type"] == "stable_domain", "domain_id"].nunique()) if len(summary_df) else 0
    module_n = int(summary_df.loc[summary_df["unit_type"] == "small_module", "domain_id"].nunique()) if len(summary_df) else 0
    mode = "DOMAIN" if stable_n else ("MODULE" if module_n else "PATHWAY_MODULES")
    atomic_text(output_path(cfg, root, "domain_freeze_summary"), (
        "MDV-5 domain/module gene-set freeze summary\n"
        f"stable_domain_n={stable_n}\n"
        f"small_module_n={module_n}\n"
        f"primary_frozen_unit_n={len(unit_ids)}\n"
        f"standalone_or_unstable_representative_n={len(standalone)}\n"
        f"next_stage_mode={mode}\n"
        f"freeze_timestamp={stamp}\n"
        "association_results_read=NO\n"
        "biological_labels=PLACEHOLDER_ONLY_GWAS_BLINDED\n"
    ))
    atomic_text(output_path(cfg, root, "domain_freeze_complete"), (
        "status=PASS\n"
        "technical_status=PASS\n"
        "stage=MDV5_DOMAIN_FREEZE\n"
        f"primary_frozen_unit_n={len(unit_ids)}\n"
        f"next_stage_mode={mode}\n"
        "gwas_inputs_used=NO\n"
    ))
    print(f"[PASS] Domain freeze: units={len(unit_ids)}, mode={mode}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MDV5Error as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        raise SystemExit(2)

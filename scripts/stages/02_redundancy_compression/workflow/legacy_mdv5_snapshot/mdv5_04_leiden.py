#!/usr/bin/env python3
from __future__ import annotations

import argparse
import inspect
import itertools
import math
import random
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

from mdv5_common import (
    MDV5Error, atomic_text, load_config, outdir, output_path, read_table, root_from_cfg,
    sentinel_passes, write_pass, write_table,
)


def canonical_membership(keys, raw_membership):
    groups = defaultdict(list)
    for k, lab in zip(keys, raw_membership):
        groups[int(lab)].append(k)
    ordered = sorted(groups.values(), key=lambda xs: min(xs))
    mapping = {}
    for idx, members in enumerate(ordered, start=1):
        for k in members:
            mapping[k] = idx
    return [mapping[k] for k in keys]


def pair_coclustering(memberships: np.ndarray) -> np.ndarray:
    runs, n = memberships.shape
    out = np.eye(n, dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            p = float(np.mean(memberships[:, i] == memberships[:, j]))
            out[i, j] = out[j, i] = p
    return out


def medoid_run(memberships: np.ndarray, coclust: np.ndarray, seeds: list[int]) -> tuple[int, float]:
    n = memberships.shape[1]
    upper = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if not upper:
        return 0, 0.0
    losses = []
    for r in range(memberships.shape[0]):
        vals = [abs((1.0 if memberships[r, i] == memberships[r, j] else 0.0) - coclust[i, j]) for i, j in upper]
        losses.append(float(np.mean(vals)))
    best = min(range(len(losses)), key=lambda r: (losses[r], seeds[r]))
    return best, losses[best]


def comb2(n: int) -> int:
    return n * (n - 1) // 2


def adjusted_rand_index(a: list[int], b: list[int]) -> float:
    if len(a) != len(b):
        raise ValueError("ARI vectors differ in length")
    n = len(a)
    if n < 2:
        return 1.0
    ct = defaultdict(int)
    ca = defaultdict(int)
    cb = defaultdict(int)
    for x, y in zip(a, b):
        ct[(x, y)] += 1
        ca[x] += 1
        cb[y] += 1
    sum_nij = sum(comb2(v) for v in ct.values())
    sum_ai = sum(comb2(v) for v in ca.values())
    sum_bj = sum(comb2(v) for v in cb.values())
    total = comb2(n)
    expected = (sum_ai * sum_bj) / total if total else 0.0
    max_index = 0.5 * (sum_ai + sum_bj)
    denom = max_index - expected
    if denom == 0:
        return 1.0
    return (sum_nij - expected) / denom


def run_leiden_for_threshold(keys, full_mat, idx_map, threshold, cfg):
    try:
        import igraph as ig
    except Exception as exc:
        raise MDV5Error(f"python-igraph is required for Leiden: {exc}")

    edges = []
    weights = []
    edge_rows = []
    for i, j in itertools.combinations(range(len(keys)), 2):
        jac = float(full_mat[idx_map[keys[i]], idx_map[keys[j]]])
        if jac + 1e-15 >= threshold:
            edges.append((i, j))
            weights.append(jac)
            edge_rows.append({"edge_threshold": threshold, "pathway_key_1": keys[i], "pathway_key_2": keys[j], "jaccard": jac})
    g = ig.Graph(n=len(keys), edges=edges, directed=False)
    g.vs["name"] = keys
    if edges:
        g.es["weight"] = weights

    n_runs = int(cfg["clustering"]["n_runs"])
    seed_base = int(cfg["clustering"]["seed_base"])
    seeds = [seed_base + r for r in range(n_runs)]
    run_memberships = []
    run_quality = []

    if len(edges) == 0:
        base = list(range(1, len(keys) + 1))
        run_memberships = [base[:] for _ in seeds]
        run_quality = [float("nan") for _ in seeds]
    else:
        sig = inspect.signature(g.community_leiden)
        for seed in seeds:
            random.seed(seed)
            np.random.seed(seed % (2**32 - 1))
            try:
                ig.set_random_number_generator(random)
            except Exception:
                pass
            kwargs = dict(
                objective_function=str(cfg["clustering"]["objective_function"]),
                weights="weight",
                beta=float(cfg["clustering"]["beta"]),
                n_iterations=int(cfg["clustering"]["n_iterations"]),
            )
            if "resolution" in sig.parameters:
                kwargs["resolution"] = float(cfg["clustering"]["resolution"])
            else:
                kwargs["resolution_parameter"] = float(cfg["clustering"]["resolution"])
            cl = g.community_leiden(**kwargs)
            mem = canonical_membership(keys, cl.membership)
            run_memberships.append(mem)
            try:
                run_quality.append(float(cl.quality))
            except Exception:
                run_quality.append(float("nan"))

    arr = np.asarray(run_memberships, dtype=int)
    co = pair_coclustering(arr)
    med_idx, loss = medoid_run(arr, co, seeds)
    med = arr[med_idx, :].tolist()
    return {
        "threshold": threshold,
        "edges": edge_rows,
        "memberships": arr,
        "quality": run_quality,
        "seeds": seeds,
        "coclust": co,
        "medoid_index": med_idx,
        "medoid_seed": seeds[med_idx],
        "medoid_loss": loss,
        "medoid_membership": med,
    }


def classify_primary(keys, med, co, cfg):
    groups = defaultdict(list)
    for i, c in enumerate(med):
        groups[int(c)].append(i)
    stable_min_n = int(cfg["clustering"]["stable_domain"]["min_representative_pathways"])
    stable_min_p = float(cfg["clustering"]["stable_domain"]["median_within_domain_coclustering_min"])
    module_n = int(cfg["clustering"]["small_module"]["exact_representative_pathways"])
    module_min_p = float(cfg["clustering"]["small_module"]["pair_coclustering_min"])
    rows = []
    stable_counter = 0
    module_counter = 0
    for c in sorted(groups):
        idxs = groups[c]
        vals = [co[i, j] for pos, i in enumerate(idxs) for j in idxs[pos+1:]]
        medp = float(np.median(vals)) if vals else float("nan")
        if len(idxs) >= stable_min_n and medp >= stable_min_p:
            stable_counter += 1
            unit_type = "stable_domain"
            unit_id = f"D{stable_counter:02d}"
        elif len(idxs) == module_n and medp >= module_min_p:
            module_counter += 1
            unit_type = "small_module"
            unit_id = f"M{module_counter:02d}"
        elif len(idxs) == 1:
            unit_type = "standalone_representative"
            unit_id = ""
        else:
            unit_type = "unstable_cluster"
            unit_id = ""
        for i in idxs:
            rows.append({
                "pathway_key": keys[i],
                "medoid_cluster": c,
                "community_size": len(idxs),
                "median_within_community_coclustering": medp,
                "unit_type": unit_type,
                "unit_id": unit_id,
            })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    root = root_from_cfg(cfg)
    outdir(cfg, root)
    if not sentinel_passes(output_path(cfg, root, "redundancy_pass")):
        raise MDV5Error("MDV5_REDUNDANCY_PASS missing/invalid")

    reps = read_table(output_path(cfg, root, "representatives"))
    reps = reps.sort_values("pathway_key", kind="mergesort").reset_index(drop=True)
    keys = reps["pathway_key"].astype(str).tolist()
    jdf = read_table(output_path(cfg, root, "jaccard_matrix"))
    all_keys = jdf.iloc[:, 0].astype(str).tolist()
    full_mat = jdf.iloc[:, 1:].to_numpy(dtype=float)
    idx_map = {k: i for i, k in enumerate(all_keys)}
    missing = [k for k in keys if k not in idx_map]
    if missing:
        raise MDV5Error(f"Representatives missing from Jaccard matrix: {missing}")

    min_rep = int(cfg["clustering"]["proceed_to_leiden_min_representatives"])
    thresholds = [float(cfg["clustering"]["primary_edge_threshold"])] + [float(x) for x in cfg["clustering"].get("sensitivity_edge_thresholds", [])]
    thresholds = list(dict.fromkeys(thresholds))

    # For rep_n < 5 we do not force Leiden; deterministic singleton/pathway-module mode is written.
    results = {}
    if len(keys) >= min_rep:
        for t in thresholds:
            results[t] = run_leiden_for_threshold(keys, full_mat, idx_map, t, cfg)
    else:
        for t in thresholds:
            n = len(keys)
            med = list(range(1, n + 1))
            arr = np.tile(np.asarray(med, dtype=int), (int(cfg["clustering"]["n_runs"]), 1))
            co = np.eye(n, dtype=float)
            results[t] = {
                "threshold": t, "edges": [], "memberships": arr,
                "quality": [float("nan")] * arr.shape[0],
                "seeds": [int(cfg["clustering"]["seed_base"]) + r for r in range(arr.shape[0])],
                "coclust": co, "medoid_index": 0,
                "medoid_seed": int(cfg["clustering"]["seed_base"]), "medoid_loss": 0.0,
                "medoid_membership": med,
            }

    primary_t = float(cfg["clustering"]["primary_edge_threshold"])
    primary = results[primary_t]
    primary_partition = classify_primary(keys, primary["medoid_membership"], primary["coclust"], cfg)

    # All run assignments, all thresholds.
    run_rows = []
    all_edges = []
    medoid_rows = []
    sens_rows = []
    for t in thresholds:
        res = results[t]
        all_edges.extend(res["edges"])
        for r in range(res["memberships"].shape[0]):
            for i, key in enumerate(keys):
                run_rows.append({
                    "edge_threshold": t, "run_index": r + 1, "seed": res["seeds"][r],
                    "pathway_key": key, "cluster_id": int(res["memberships"][r, i]),
                    "quality": res["quality"][r],
                })
        medoid_rows.append({
            "edge_threshold": t, "medoid_run_index": res["medoid_index"] + 1,
            "medoid_seed": res["medoid_seed"], "medoid_loss": res["medoid_loss"],
            "edge_n": len(res["edges"]), "community_n": len(set(res["medoid_membership"])),
        })
        ari = adjusted_rand_index(primary["medoid_membership"], res["medoid_membership"])
        tmp_part = classify_primary(keys, res["medoid_membership"], res["coclust"], cfg)
        sens_rows.append({
            "edge_threshold": t, "is_primary": int(t == primary_t), "edge_n": len(res["edges"]),
            "community_n": len(set(res["medoid_membership"])),
            "ARI_vs_primary_medoid": ari,
            "stable_domain_n": int(tmp_part.loc[tmp_part["unit_type"] == "stable_domain", "unit_id"].nunique()) if len(tmp_part) else 0,
            "small_module_n": int(tmp_part.loc[tmp_part["unit_type"] == "small_module", "unit_id"].nunique()) if len(tmp_part) else 0,
            "medoid_seed": res["medoid_seed"], "medoid_loss": res["medoid_loss"],
        })

    write_table(pd.DataFrame(all_edges), output_path(cfg, root, "network_edges_sensitivity"))
    primary_edges = pd.DataFrame(primary["edges"])
    if primary_edges.empty:
        primary_edges = pd.DataFrame(columns=["edge_threshold", "pathway_key_1", "pathway_key_2", "jaccard"])
    write_table(primary_edges, output_path(cfg, root, "network_edges"))
    write_table(pd.DataFrame(run_rows), output_path(cfg, root, "leiden_runs"))
    write_table(pd.DataFrame(medoid_rows), output_path(cfg, root, "medoid_runs"))
    write_table(pd.DataFrame(sens_rows), output_path(cfg, root, "threshold_sensitivity"))

    co = primary["coclust"]
    codf = pd.DataFrame(co, columns=keys)
    codf.insert(0, "pathway_key", keys)
    write_table(codf, output_path(cfg, root, "coclustering_matrix"))
    cpairs = []
    for i, j in itertools.combinations(range(len(keys)), 2):
        cpairs.append({"pathway_key_1": keys[i], "pathway_key_2": keys[j], "coclustering_probability": co[i, j]})
    write_table(pd.DataFrame(cpairs), output_path(cfg, root, "coclustering_pairs"))
    write_table(primary_partition, output_path(cfg, root, "primary_partition"))

    stable_n = int(primary_partition.loc[primary_partition["unit_type"] == "stable_domain", "unit_id"].nunique())
    module_n = int(primary_partition.loc[primary_partition["unit_type"] == "small_module", "unit_id"].nunique())
    standalone_n = int((primary_partition["unit_type"] == "standalone_representative").sum())
    unstable_n = int(primary_partition.loc[primary_partition["unit_type"] == "unstable_cluster", "medoid_cluster"].nunique())
    mode = "DOMAIN" if stable_n > 0 else ("MODULE" if module_n > 0 else "PATHWAY_MODULES")
    summary = (
        "MDV-5 Leiden consensus summary\n"
        f"representative_n={len(keys)}\n"
        f"leiden_executed={'YES' if len(keys) >= min_rep else 'NO_GATE_REP_LT_5'}\n"
        f"primary_edge_threshold={primary_t}\n"
        f"primary_edge_n={len(primary['edges'])}\n"
        f"runs={cfg['clustering']['n_runs']}\n"
        f"medoid_seed={primary['medoid_seed']}\n"
        f"stable_domain_n={stable_n}\n"
        f"small_module_n={module_n}\n"
        f"standalone_representative_n={standalone_n}\n"
        f"unstable_cluster_n={unstable_n}\n"
        f"next_stage_mode={mode}\n"
        "gwas_inputs_used=NO\n"
    )
    atomic_text(output_path(cfg, root, "leiden_summary"), summary)
    write_pass(output_path(cfg, root, "leiden_pass"), {
        "status": "PASS",
        "technical_status": "PASS",
        "stage": "MDV5_LEIDEN_CONSENSUS",
        "representative_n": len(keys),
        "stable_domain_n": stable_n,
        "small_module_n": module_n,
        "next_stage_mode": mode,
        "gwas_inputs_used": "NO",
    })
    print(f"[PASS] Leiden/consensus: reps={len(keys)}, stable_domains={stable_n}, modules={module_n}, mode={mode}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MDV5Error as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        raise SystemExit(2)

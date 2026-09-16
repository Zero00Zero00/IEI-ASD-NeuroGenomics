#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from mdv5_common import (
    MDV5Error, load_config, outdir, output_path, read_table, root_from_cfg, sentinel_passes,
    write_pass, write_table, atomic_text,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    root = root_from_cfg(cfg)
    outdir(cfg, root)
    if not sentinel_passes(output_path(cfg, root, "jaccard_pass")):
        raise MDV5Error("MDV5_JACCARD_PASS missing/invalid")

    tier = read_table(output_path(cfg, root, "tier1_annotated"))
    jdf = read_table(output_path(cfg, root, "jaccard_matrix"))
    keys = jdf.iloc[:, 0].astype(str).tolist()
    mat = jdf.iloc[:, 1:].to_numpy(dtype=float)
    if list(jdf.columns[1:]) != keys:
        raise MDV5Error("Jaccard matrix row/column key order mismatch")
    tier = tier.set_index("pathway_key", drop=False).loc[keys].reset_index(drop=True)

    n = len(keys)
    thr = float(cfg["redundancy"]["high_redundancy_jaccard"])
    if n == 1:
        raw_labels = np.array([1])
    else:
        dist = 1.0 - mat
        np.fill_diagonal(dist, 0.0)
        condensed = squareform(dist, checks=True)
        z = linkage(condensed, method="complete", optimal_ordering=True)
        raw_labels = fcluster(z, t=1.0 - thr + 1e-12, criterion="distance")

    groups = {}
    for i, lab in enumerate(raw_labels):
        groups.setdefault(int(lab), []).append(i)

    src_pri = cfg["redundancy"].get("source_priority", {"GO_BP": 1, "Reactome": 2})
    tier["source_priority"] = tier["source"].map(lambda x: int(src_pri.get(str(x), 99)))

    family_records = []
    family_defs = []
    # Canonical family order by deterministic representative key after selection.
    tmp_families = []
    for raw_lab, idxs in groups.items():
        sub = tier.iloc[idxs].copy()
        sub = sub.sort_values(
            ["q_emp", "q_Firth", "pathway_size", "source_priority", "pathway_id"],
            ascending=[True, True, True, True, True], kind="mergesort"
        ).reset_index(drop=True)
        rep_key = str(sub.iloc[0]["pathway_key"])
        tmp_families.append((rep_key, raw_lab, idxs, sub))
    tmp_families.sort(key=lambda x: x[0])

    for fidx, (rep_key, raw_lab, idxs, sub) in enumerate(tmp_families, start=1):
        fam_id = f"RF{fidx:03d}"
        # Complete-linkage cut must imply every within-family pair J >= threshold.
        min_j = 1.0
        if len(idxs) > 1:
            vals = [mat[i, j] for pos, i in enumerate(idxs) for j in idxs[pos+1:]]
            min_j = float(min(vals))
        if min_j + 1e-12 < thr:
            raise MDV5Error(f"Complete-linkage family violates threshold: {fam_id} minJ={min_j}")
        for rank, (_, row) in enumerate(sub.iterrows(), start=1):
            family_records.append({
                **row.to_dict(),
                "redundancy_family_id": fam_id,
                "redundancy_family_n": len(idxs),
                "family_min_pairwise_jaccard": min_j,
                "representative_rank": rank,
                "representative_flag": int(str(row["pathway_key"]) == rep_key),
                "representative_pathway_key": rep_key,
                "representative_rule": "q_emp->q_Firth->pathway_size->source_priority->pathway_id",
            })
        family_defs.append({
            "redundancy_family_id": fam_id,
            "family_n": len(idxs),
            "representative_pathway_key": rep_key,
            "family_min_pairwise_jaccard": min_j,
        })

    rmap = pd.DataFrame(family_records).sort_values(["redundancy_family_id", "representative_rank"], kind="mergesort")
    reps = rmap.loc[rmap["representative_flag"] == 1].copy()
    write_table(rmap, output_path(cfg, root, "redundancy_map"))
    write_table(reps, output_path(cfg, root, "representatives"))

    rep_n = len(reps)
    mode = "LEIDEN_ELIGIBLE" if rep_n >= int(cfg["clustering"]["proceed_to_leiden_min_representatives"]) else ("PATHWAY_MODULE_MODE" if rep_n >= 2 else "LIMITED_EVIDENCE_MODE")
    summary = (
        f"MDV-5 redundancy compression summary\n"
        f"tier1_n={len(tier)}\n"
        f"redundancy_family_n={len(tmp_families)}\n"
        f"representative_n={rep_n}\n"
        f"high_redundancy_jaccard={thr}\n"
        f"family_algorithm=complete_linkage_on_1_minus_Jaccard\n"
        f"next_mode={mode}\n"
        f"gwas_inputs_used=NO\n"
    )
    atomic_text(output_path(cfg, root, "redundancy_summary"), summary)
    write_pass(output_path(cfg, root, "redundancy_pass"), {
        "status": "PASS",
        "technical_status": "PASS",
        "stage": "MDV5_REDUNDANCY_COMPRESSION",
        "tier1_n": len(tier),
        "representative_n": rep_n,
        "next_mode": mode,
        "gwas_inputs_used": "NO",
    })
    print(f"[PASS] Redundancy compression: 40 Tier-1 -> {rep_n} representatives; mode={mode}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MDV5Error as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        raise SystemExit(2)

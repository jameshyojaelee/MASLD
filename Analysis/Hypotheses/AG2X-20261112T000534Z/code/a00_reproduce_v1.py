#!/usr/bin/env python3
"""Stage A step 5: reproduce the archived AlphaGenome version 1 numbers before any version 2 statistic.

Reuses the earlier program's loaders and statistics by import (read-only):
  scripts/analysis/alphagenome_program/c1e2_05_recompute_v2.py  (statistics, block bootstrap, Atlas loaders)
  scripts/analysis/alphagenome_atlas/lib_atlas.py              (gate_set for A3a)
Checks, each to 3 decimals against the deposit:
  A1 tier A4: alphagenome_atac_liver 0.7361; controls -0.0008 and -0.0183 (c1e2_model_metrics.tsv)
  A2 495 sites: weighted 0.4507; controls 0.1079 and 0.0708
  A3a 201 leads: AUROC 0.5304, sign concordance 0.5224 (ag_atlas_gate_replication.json)
The v1 score columns are also rebuilt from the raw Atlas archives and compared with the matched tables.
Writes tables/v1_reproduction.tsv and tables/v1_reproduction.json, plus the v1 columns for the
hepatocyte-lineage A2 sites (new numbers, no archived value).
"""
from __future__ import annotations

import json
import math
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

sys.path.insert(0, str(L.REPO / "scripts/analysis/alphagenome_program"))
sys.path.insert(0, str(L.REPO / "scripts/analysis/alphagenome_atlas"))
import c1e2_05_recompute_v2 as C  # noqa: E402
import lib_atlas as la  # noqa: E402

PROG = L.REPO / "GWAS/finemapping/results/alphagenome_program"
C1 = PROG / "c1-endpoint2-v2-20260915T084559Z/tables"
D1 = PROG / "d1-lineage-allelic-20260914T194213Z/tables/lineage_site_summary.tsv"
RUN = L.REPO / "GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z"
OUT = L.PKG / "tables"


def metric(rows, ep, model, m):
    for r in rows:
        if r["endpoint"] == ep and r["model"] == model and r["metric"] == m:
            return r["value"]
    return float("nan")


def main():
    res = {"checks": []}
    arch = pd.read_csv(C1 / "c1e2_model_metrics.tsv", sep="\t")

    def archived(ep, model, m):
        s = arch[(arch.endpoint == ep) & (arch.model == model) & (arch.metric == m)]
        return float(s.value.iloc[0])

    # ---------------- A1 ----------------
    a1 = pd.read_csv(C1 / "endpoint1_matched_tierA4.tsv.gz", sep="\t")
    atlas = C.load_atlas_chunks()
    chk = a1[["key", "alphagenome_atac_liver"]].merge(atlas[["key", "alphagenome_atac_liver"]], on="key", suffixes=("", "_raw"))
    res["A1_rebuild_n"] = int(len(chk))
    res["A1_rebuild_max_abs_diff"] = float((chk.alphagenome_atac_liver - chk.alphagenome_atac_liver_raw).abs().max())
    rows, boot = [], {}
    models = ["alphagenome_atac_liver", "alphagenome_dnase_liver", "control_allele_identity", "control_position"]
    C.run_endpoint(a1, models, "beta_alt", None, "endpoint1_currin_tierA4", rows, boot)
    for m in models:
        for met in ("signed_spearman_unweighted", "signed_spearman_unweighted_boot_ci_lo", "signed_spearman_unweighted_boot_ci_hi"):
            ours, theirs = metric(rows, "endpoint1_currin_tierA4", m, met), archived("endpoint1_currin_tierA4", m, met)
            res["checks"].append({"endpoint": "A1_tierA4", "model": m, "metric": met, "reproduced": ours, "archived": theirs,
                                  "n": len(a1), "blocks": int(a1.block_1mb.nunique()),
                                  "match_3dp": round(ours, 3) == round(theirs, 3)})

    # ---------------- A2 ----------------
    a2 = pd.read_csv(C1 / "endpoint2_matched.tsv", sep="\t")
    t0 = C.load_track0_ase(a2.uid.tolist())
    chk2 = a2[["uid", "alphagenome_atac_liver"]].merge(t0.rename(columns={"variant_uid": "uid"})[["uid", "alphagenome_atac_liver"]], on="uid", suffixes=("", "_raw"))
    res["A2_rebuild_n"] = int(len(chk2))
    res["A2_rebuild_max_abs_diff"] = float((chk2.alphagenome_atac_liver - chk2.alphagenome_atac_liver_raw).abs().max())
    rows2, boot2 = [], {}
    C.run_endpoint(a2, models, "mean_log2_alt_over_ref", "precision_weight", "endpoint2_allelic_imbalance", rows2, boot2)
    for m in models:
        for met in ("signed_spearman_weighted", "signed_spearman_weighted_boot_ci_lo", "signed_spearman_weighted_boot_ci_hi",
                    "signed_spearman_unweighted"):
            ours, theirs = metric(rows2, "endpoint2_allelic_imbalance", m, met), archived("endpoint2_allelic_imbalance", m, met)
            res["checks"].append({"endpoint": "A2_495", "model": m, "metric": met, "reproduced": ours, "archived": theirs,
                                  "n": len(a2), "blocks": int(a2.block_1mb.nunique()),
                                  "match_3dp": round(ours, 3) == round(theirs, 3)})

    # ---------------- A2 hepatocyte lineage (new v1 number; no archived value) ----------------
    d1 = pd.read_csv(D1, sep="\t")
    hep = d1[(d1.scope == "POOLED_30_DONORS") & (d1.lineage == "Hepatocyte")].copy()
    t0h = C.load_track0_ase(hep.uid.tolist()).rename(columns={"variant_uid": "uid"})
    hep = hep.merge(t0h, on="uid", how="left")
    hep["block_1mb"] = hep["block"]
    hep["precision_weight"] = 1.0 / hep.se_log2
    hep["control_allele_identity"] = [C.NT_RANK.get(a, np.nan) - C.NT_RANK.get(r, np.nan) for r, a in zip(hep.ref, hep.alt)]
    hep["control_position"] = (hep.pos % 1000) / 1000.0 - 0.5
    res["A2_hep_sites"] = int(len(hep))
    res["A2_hep_v1_atac_served"] = int(hep.alphagenome_atac_liver.notna().sum())
    hep.to_csv(OUT / "a2_hepatocyte_lineage_v1.tsv", sep="\t", index=False)

    # ---------------- A3a ----------------
    import anndata
    import pyfaidx
    import atlas_archive as aa
    fa = pyfaidx.Fasta("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa",
                       as_raw=True, sequence_always_upper=True)
    truth = la.gate_set(la.read_tsv(la.PROJECT / "GWAS/finemapping/results/seqfunc/broadaway_benchmark_truth.tsv"),
                        la.read_tsv(la.PROJECT / "GWAS/finemapping/results/seqfunc/alphagenome_eqtl_scores.tsv"),
                        lambda c, s, e: fa[c][s:e])
    gate = {}
    for c in sorted((RUN / "raw/atlas_gate").glob("chunk_*")):
        a = anndata.read_h5ad(c / "RNA_SEQ.h5ad")
        mask = a.var["ontology_curie"].isin(["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]).values
        obs = a.obs.reset_index(drop=True)
        x = np.asarray(a.X)[:, mask]
        for i in range(len(obs)):
            gate[(aa.variant_uid_from_str(str(obs.at[i, "variant"])), str(obs.at[i, "gene_id"]).split(".")[0])] = float(np.nanmean(x[i]))
    recs = []
    for r in truth:
        v = gate.get((r["variant_uid"], r["ensembl"]))
        recs.append({**r, "v1_atlas_liver_rna_ref_to_alt": v})
    a3 = pd.DataFrame(recs)
    a3.to_csv(OUT / "a3a_gate_set_v1.tsv", sep="\t", index=False)
    sc = a3.dropna(subset=["v1_atlas_liver_rna_ref_to_alt"])
    s_ = (sc.orientation_sign * sc.v1_atlas_liver_rna_ref_to_alt).values
    l_ = sc.label.values.astype(int)
    from scipy.stats import rankdata
    r_ = rankdata(s_)
    n1, n0 = int(l_.sum()), int((l_ == 0).sum())
    auroc = float((r_[l_ == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
    conc = float(((s_ > 0) == (l_ == 1)).mean())
    gj = json.load(open(RUN / "tables/ag_atlas_gate_replication.json"))
    res["checks"].append({"endpoint": "A3a_gate", "model": "alphagenome_v1_atlas_liver_rna", "metric": "auroc", "reproduced": auroc,
                          "archived": gj["auroc"], "n": int(len(sc)), "blocks": None, "match_3dp": round(auroc, 3) == round(gj["auroc"], 3)})
    res["checks"].append({"endpoint": "A3a_gate", "model": "alphagenome_v1_atlas_liver_rna", "metric": "sign_concordance", "reproduced": conc,
                          "archived": gj["sign_concordance"], "n": int(len(sc)), "blocks": None,
                          "match_3dp": round(conc, 3) == round(gj["sign_concordance"], 3)})
    res["A3a_leads"] = int(len(a3))
    res["A3a_scored_v1"] = int(len(sc))
    res["all_match_3dp"] = bool(all(c["match_3dp"] for c in res["checks"]))
    pd.DataFrame(res["checks"]).to_csv(OUT / "v1_reproduction.tsv", sep="\t", index=False)
    json.dump(res, open(OUT / "v1_reproduction.json", "w"), indent=1, default=str)
    L.log(json.dumps({k: v for k, v in res.items() if k != "checks"}, default=str))
    print(pd.DataFrame(res["checks"]).to_string())


if __name__ == "__main__":
    main()

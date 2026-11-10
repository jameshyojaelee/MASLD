#!/usr/bin/env python
"""CGB17 (+B14, B15): Mendelian, genetic-panel and drug-target lists vs expression-matched random sets
on |cont_z|, |disease_t| and colocalized status; Gobeil pQTL and Zhu MPRA annotation. Follows PRESPEC_CGB17.md."""
import argparse
import json
import os

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
SEED, NRAND = 20261022, 1000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True)
    ap.add_argument("--weights", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=False)
    rng = np.random.default_rng(SEED)
    g = pd.read_csv(a.genes, sep="\t")
    g = g[g.cont_z.notna() & g.disease_t.notna() & g.AveExpr.notna()].copy()
    g["abs_cont"], g["abs_dis"] = g.cont_z.abs(), g.disease_t.abs()
    g["coloc"] = g.coloc_gen.fillna(False).astype(float)
    sym = g.dropna(subset=["gene_name"]).drop_duplicates("gene_name").set_index("gene_name").g
    dt = pd.read_csv(f"{ROOT}/data/external/drug_targets/drug_target_classification.tsv", sep="\t", low_memory=False)
    dt["g"] = dt.ensembl_id.str.split(".").str[0]
    mon = pd.read_csv(f"{ROOT}/data/external/okur2026_monogenic/okur2026_panel_genes_partial.tsv", sep="\t")
    man = pd.read_csv(f"{ROOT}/data/published_gene_panels/mancina_2025_panel.tsv", sep="\t", comment="#")
    lists = {
        "L1_okur_mendelian_partial": set(mon.gene.map(sym).dropna()),
        "L2_mancina_germline": set(man[man.origin == "germline"].gene.map(sym).dropna()),
        "L3_masld_clinical_targets": set(dt[dt.drug_dev_status.isin(["masld_approved", "masld_clinical", "masld_discontinued"])].g),
        "L4_drugged_other_indication": set(dt[dt.drug_dev_status == "drugged_other_indication"].g),
    }
    g["dec"] = pd.qcut(g.AveExpr, 10, labels=False)
    bins = {k: g.index[g.dec == k].to_numpy() for k in range(10)}
    rows = []
    for name, members in lists.items():
        m = g[g.g.isin(members)]
        counts = m.dec.value_counts()
        draws = [np.concatenate([rng.choice(bins[k], n, replace=False) for k, n in counts.items()]) for _ in range(NRAND)]
        for stat in ["abs_cont", "abs_dis", "coloc"]:
            obs = m[stat].mean()
            null = np.array([g.loc[ix, stat].mean() for ix in draws])
            rows.append(dict(list=name, n_listed=len(members), n_tested=len(m), statistic=stat, observed=obs,
                             null_median=np.median(null), null_lo=np.quantile(null, 0.025), null_hi=np.quantile(null, 0.975),
                             emp_p=(1 + (null >= obs).sum()) / (NRAND + 1)))
    t = pd.DataFrame(rows)
    t["bh"] = multipletests(t.emp_p, method="fdr_bh")[1]
    t.to_csv(f"{a.out}/list_tests.tsv", sep="\t", index=False)

    w = pd.read_csv(a.weights, sep="\t")
    sysmap = w.groupby("gene_id_base").community_id.apply(lambda s: ";".join(sorted(set(s))))
    gob = pd.read_csv(f"{ROOT}/data/external/gobeil2026_liver_pqtl/gobeil2026_liver_pqtl_causal.csv")
    zhu = pd.read_csv(f"{ROOT}/data/external/zhu2026_masld_mpra/zhu2026_masld_mpra_calls.csv")
    cols = ["g", "gene_name", "cont_z", "cont_q", "disease_t", "disease_padj", "coloc_gen"]
    ann = []
    for src, df, extra in [("gobeil2026_pqtl", gob, ["layer_relationship", "gobeil_pqtl_risk_direction"]),
                           ("zhu2026_mpra", zhu, ["zhu_finding_type", "zhu_context", "zhu_direction"])]:
        x = df[["gene"] + extra].rename(columns={"gene": "gene_name"}).merge(g[cols], on="gene_name", how="left")
        x["source"] = src
        x["systems"] = x.g.map(sysmap)
        ann.append(x)
    pd.concat(ann).to_csv(f"{a.out}/b14_b15_annotation.tsv", sep="\t", index=False)
    json.dump({"tests": t.to_dict(orient="records"), "seed": SEED}, open(f"{a.out}/results.json", "w"), indent=1, default=float)
    print(t.to_string())


if __name__ == "__main__":
    main()

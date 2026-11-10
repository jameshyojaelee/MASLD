#!/usr/bin/env python
"""CGA7: candidate Catalog r9d sidecar with Wave 1/2 follow-up fields. Follows PRESPEC_CGA7.md."""
import argparse
import glob
import json
import os

import pandas as pd

H = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses"
R9D = f"{H}/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.parquet"
CG0 = f"{H}/CG0-20261022T122300Z/out-20261022T122649Z"
DT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/drug_targets/drug_target_classification.tsv"
MAGMA_TRAITS = ["ukbb_alt", "ukbb_ast", "mvp_alt_eur", "mvp_ast_eur", "pazoki_pdff"]


def latest(pattern):
    hits = sorted(glob.glob(pattern))
    if not hits:
        raise SystemExit(f"missing {pattern}")
    return hits[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=False)
    cat = pd.read_parquet(R9D, columns=["ensembl_id", "symbol", "primary_evidence_class"])
    cat["g"] = cat.ensembl_id.str.split(".").str[0]
    sources = {"r9d": R9D}

    g = pd.read_csv(f"{CG0}/genes.tsv", sep="\t")
    f = g[["g", "cont_z", "cont_q", "continuum_associated", "loeuf", "pli"]].rename(
        columns={"cont_z": "cg_cont_z", "cont_q": "cg_cont_q", "continuum_associated": "cg_continuum_associated",
                 "loeuf": "cg_loeuf", "pli": "cg_pli"})
    sources["cg0_genes"] = f"{CG0}/genes.tsv"

    p = latest(f"{H}/CGB8-20261022T122538Z/eqtl-*/broadaway_gene_counts.tsv")
    sources["broadaway"] = p
    f = f.merge(pd.read_csv(p, sep="\t")[["g", "n_p1e5"]].rename(columns={"n_p1e5": "cg_broadaway_eqtl_n_p1e5"}),
                on="g", how="outer")

    p = latest(f"{H}/CGB11-20261022T122538Z/out-*/gene_directions.tsv")
    sources["b11"] = p
    d = pd.read_csv(p, sep="\t")
    d["cg_coloc_direction"] = d.gen_sign.map({1: "risk_haplotype_increases_expression", -1: "risk_haplotype_decreases_expression"})
    d["cg_disease_direction_agrees"] = (d.gen_sign == d.disease_t.apply(lambda x: 1 if x > 0 else -1)).where(d.disease_padj < 0.05)
    d["cg_continuum_direction_agrees"] = (d.gen_sign == d.cont_z.apply(lambda x: 1 if x > 0 else -1)).where(d.cont_q < 0.05)
    f = f.merge(d[["g", "cg_coloc_direction", "cg_disease_direction_agrees", "cg_continuum_direction_agrees"]], on="g", how="outer")

    mrun = latest(f"{H}/CG0-20261022T122300Z/magma/run-*")
    sources["magma_run"] = mrun
    for t in MAGMA_TRAITS:
        m = pd.read_csv(f"{mrun}/{t}.genes.out", sep=r"\s+")
        f = f.merge(m[["GENE", "ZSTAT"]].rename(columns={"GENE": "g", "ZSTAT": f"cg_magma_z_{t}"}), on="g", how="outer")

    p = latest(f"{H}/CGB2-20261022T180846Z/out-*/allelic_results.tsv")
    sources["b2"] = p
    b2 = pd.read_csv(p, sep="\t")
    b2 = b2.rename(columns={"gene": "symbol"})[["symbol", "variant", "n_het", "mean_alt_fraction", "p1_beta", "p1_se", "p1_testable"]]
    b2.columns = ["symbol"] + [f"cg_allelic_share_{c}" for c in b2.columns[1:]]

    p = latest(f"{H}/CGB17-20261022T180846Z/out-*/b14_b15_annotation.tsv")
    sources["b14_b15"] = p
    ann = pd.read_csv(p, sep="\t")
    gob = ann[ann.source == "gobeil2026_pqtl"].dropna(subset=["g"])[["g", "layer_relationship"]].rename(
        columns={"layer_relationship": "cg_gobeil_layer_relationship"})
    zhu = (ann[ann.source == "zhu2026_mpra"].dropna(subset=["g"]).groupby("g")
           .apply(lambda x: ";".join(sorted(set(x.zhu_context.dropna())))).rename("cg_zhu_mpra_context").reset_index())
    f = f.merge(gob, on="g", how="outer").merge(zhu, on="g", how="outer")

    dt = pd.read_csv(DT, sep="\t", usecols=["ensembl_id", "drug_dev_status"], low_memory=False).dropna()
    dt["g"] = dt.ensembl_id.str.split(".").str[0]
    f = f.merge(dt.drop_duplicates("g")[["g", "drug_dev_status"]].rename(columns={"drug_dev_status": "cg_drug_dev_status"}),
                on="g", how="outer")

    out = cat.merge(f, on="g", how="left").merge(b2, on="symbol", how="left").drop(columns="g")
    out["cg_composition_sensitivity"] = "not_computed"
    assert len(out) == len(cat), "row count changed"
    out.to_csv(f"{a.out}/catalog_r9d_sidecar_cg_followup.tsv.gz", sep="\t", index=False)
    conc = out[out.primary_evidence_class.astype(str).str.contains("concordant", case=False)]
    conc[["ensembl_id", "symbol", "cg_coloc_direction", "cg_disease_direction_agrees", "cg_continuum_direction_agrees",
          "cg_cont_z", "cg_loeuf"]].to_csv(f"{a.out}/concordant_class_directions.tsv", sep="\t", index=False)
    summary = {"rows": len(out), "r9d_rows": len(cat), "non_missing": {c: int(out[c].notna().sum()) for c in out.columns if c.startswith("cg_")},
               "concordant_class_n": len(conc),
               "concordant_with_direction": int(conc.cg_coloc_direction.notna().sum()),
               "concordant_disease_agree": conc.cg_disease_direction_agrees.value_counts(dropna=False).to_dict(),
               "class_counts": out.primary_evidence_class.value_counts().to_dict(), "sources": sources}
    json.dump(summary, open(f"{a.out}/summary.json", "w"), indent=1, default=str)
    print(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()

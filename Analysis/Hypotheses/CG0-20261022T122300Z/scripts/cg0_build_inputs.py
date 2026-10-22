#!/usr/bin/env python
"""CG0: shared inputs for the continuum-genetics plan (A1, B3, B4, B8, B11).

Builds three tables from frozen inputs; computes no outcome contrast.
  units.tsv   774 genetic-identity units: cohort, role, stage, sex, continuum
              (signature_pc1, z within cohort), 43 system z scores
              (collection_balanced), 16 BayesPrism fractions.
  genes.tsv   one row per Ensembl gene: expression, canonical disease-state
              effect, stage-adjusted continuum effect, system multiplicity,
              COLOC r2 state (H57 universe rule), gnomAD v4.1 LOEUF, length,
              gene density, ABC liver enhancer count, GTEx liver eQTL counts.
  gene_system_weights.tsv.gz  gene x system membership weight (long).
Plus input_manifest.tsv (sha256 of every input).
"""
import argparse
import gzip
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CONT = "RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-corrected774-20261009T171017Z"
LAY = "RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z"
SYS = "RNA-seq/results/histology_anchored_continuum/molecular_systems_inference/hac-system-inference-corrected774-20261009T171017Z"
ATLAS = "figures/candidates/histology-continuum-molecular-systems-atlas-20260824T195148Z"
INPUTS = {
    "continuum_scores": f"{CONT}/unsupervised/participant_scores.tsv",
    "system_scores": f"{SYS}/tables/participant_system_scores.tsv.gz",
    "continuum_gene_meta": f"{LAY}/bulk/continuum_gene_meta.tsv.gz",
    "continuum_membership": f"{LAY}/bulk/continuum_membership.tsv.gz",
    "canonical_deg": "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv",
    "coloc_gene_level": "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv",
    "gene_feature_membership": f"{ATLAS}/source_tables/outcome_blind_gene_membership.tsv.gz",
    "feature_system": f"{ATLAS}/geometry/latent_coordinates.tsv.gz",
    "gnomad": "data/external/gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv",
    "gene_annot_grch37": "tools/pops/example/data/utils/gene_annot_jun10.txt",
    "abc": "data/external/abc_liver/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz",
    "gtex_liver_eqtl": "GWAS/MR_Data/GTEx_v8_Liver_eQTL.tsv.gz",
}
COHORTS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"]
BP = "Analysis/Deconvolution/corrected_774_20261009T141331Z/results/{c}/{c}_bayesprism_proportions.tsv"
for c in COHORTS:
    INPUTS[f"bayesprism_{c}"] = BP.format(c=c)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def zscore(x):
    return (x - x.mean()) / x.std(ddof=1)


def build_units():
    cs = pd.read_csv(INPUTS["continuum_scores"], sep="\t")
    cs = cs[cs.axis_id == "signature_pc1"][["sample_id", "dataset", "axis_raw"]]
    cs["continuum_z"] = cs.groupby("dataset").axis_raw.transform(zscore)
    ss = pd.read_csv(INPUTS["system_scores"], sep="\t")
    ss = ss[ss.aggregation_method == "collection_balanced"]
    meta = ss.drop_duplicates("sample_id")[["sample_id", "dataset", "analysis_unit_id", "inferred_sex",
                                            "fibrosis_stage", "nas_score"]]
    wide = ss.pivot(index="sample_id", columns="community_id", values="system_score_z")
    wide.columns = [f"sys_{c.split('_')[1]}" for c in wide.columns]
    reg = pd.read_csv(f"{CONT}/../hac-continuum-corrected774-20261011T152106Z/projection/score_registry.tsv",
                      sep="\t")[["sample_id", "evaluation_role"]]
    u = meta.merge(cs, on=["sample_id", "dataset"], how="left").merge(reg, on="sample_id", how="left")
    u = u.merge(wide, left_on="sample_id", right_index=True, how="left")
    frames = []
    for c in COHORTS:
        b = pd.read_csv(INPUTS[f"bayesprism_{c}"], sep="\t", index_col=0)
        b.columns = ["bp_" + x.replace(" ", "_").replace("/", "_").replace("+", "_") for x in b.columns]
        frames.append(b)
    bp = pd.concat(frames)
    u = u.merge(bp, left_on="sample_id", right_index=True, how="left")
    return u


def build_genes():
    deg = pd.read_csv(INPUTS["canonical_deg"])
    deg["g"] = deg.gene.str.split(".").str[0]
    deg = deg.drop_duplicates("g")
    deg["deg"] = (deg.padj < 0.05) & (deg.logFC.abs() > 0.5)  # H57 rule
    g = deg[["g", "symbol", "AveExpr", "logFC", "t", "padj", "deg"]].rename(
        columns={"symbol": "gene_name", "logFC": "disease_logFC", "t": "disease_t", "padj": "disease_padj"})

    cm = pd.read_csv(INPUTS["continuum_gene_meta"], sep="\t")
    cm = cm[cm.axis_id == "signature_pc1"][["gene_id_base", "beta", "z_value", "bh_q_value"]]
    cm.columns = ["g", "cont_beta", "cont_z", "cont_q"]
    mem = pd.read_csv(INPUTS["continuum_membership"], sep="\t")[["gene_id_base", "continuum_associated",
                                                                  "axis_component"]]
    mem.columns = ["g", "continuum_associated", "signature_gene"]
    g = g.merge(cm, on="g", how="outer").merge(mem, on="g", how="left")

    co = pd.read_csv(INPUTS["coloc_gene_level"])
    co["g"] = co.ensembl.astype(str).str.split(".").str[0]
    co = co.sort_values(["coloc_best_susie_pp4", "coloc_best_abf_pp4"], ascending=False).drop_duplicates("g")
    co = co[["g", "coloc_best_susie_pp4", "susie_state_t12", "susie_driving_gwas"]]
    g = g.merge(co, on="g", how="left")
    g["coloc_gen"] = (pd.to_numeric(g.coloc_best_susie_pp4, errors="coerce") >= 0.5).fillna(False)
    g["h57_ueval"] = g.susie_state_t12.isin(["positive", "tested_le_0.5"]) & g.disease_padj.notna()

    gm = pd.read_csv(INPUTS["gene_feature_membership"], sep="\t", usecols=["gene_id_base", "feature_id",
                                                                           "mapped_weight"])
    fs = pd.read_csv(INPUTS["feature_system"], sep="\t", usecols=["feature_id", "community_id"])
    gs = gm.merge(fs, on="feature_id", how="inner")
    gsw = gs.groupby(["gene_id_base", "community_id"], as_index=False).mapped_weight.sum()
    nsys = gsw.groupby("gene_id_base").size().rename("n_systems")
    nfeat = gs.groupby("gene_id_base").size().rename("n_features")
    g = g.merge(nsys, left_on="g", right_index=True, how="left").merge(nfeat, left_on="g", right_index=True,
                                                                       how="left")
    g[["n_systems", "n_features"]] = g[["n_systems", "n_features"]].fillna(0).astype(int)

    gn = pd.read_csv(INPUTS["gnomad"], sep="\t", usecols=["gene", "gene_id", "transcript", "canonical",
                                                          "lof.oe_ci.upper", "lof.pLI"], dtype=str)
    gn = gn[(gn.canonical == "true") & gn.gene_id.str.startswith("ENSG")]
    gn = gn.drop_duplicates("gene_id")
    gn = gn.rename(columns={"gene_id": "g", "lof.oe_ci.upper": "loeuf", "lof.pLI": "pli"})[["g", "loeuf", "pli"]]
    gn[["loeuf", "pli"]] = gn[["loeuf", "pli"]].apply(pd.to_numeric, errors="coerce")
    g = g.merge(gn, on="g", how="left")

    an = pd.read_csv(INPUTS["gene_annot_grch37"], sep="\t", dtype={"CHR": str})
    an = an.drop_duplicates("ENSGID")
    an["gene_length"] = an.END - an.START + 1
    dens = []
    for chrom, d in an.groupby("CHR"):
        tss = np.sort(d.TSS.values)
        lo = np.searchsorted(tss, d.TSS.values - 500_000, side="left")
        hi = np.searchsorted(tss, d.TSS.values + 500_000, side="right")
        dens.append(pd.Series(hi - lo - 1, index=d.ENSGID.values))
    an["gene_density_1mb"] = an.ENSGID.map(pd.concat(dens))
    an = an.rename(columns={"ENSGID": "g", "CHR": "chr", "START": "start", "END": "end", "TSS": "tss"})
    g = g.merge(an[["g", "chr", "start", "end", "tss", "gene_length", "gene_density_1mb"]], on="g", how="left")

    abc = pd.read_csv(INPUTS["abc"], sep="\t", usecols=["name", "class", "TargetGene", "CellType"])
    abc = abc[(abc.CellType == "liver-ENCODE") & (abc["class"] != "promoter")]
    nenh = abc.groupby("TargetGene").name.nunique().rename("abc_liver_n_enh")
    g = g.merge(nenh, left_on="gene_name", right_index=True, how="left")
    g["abc_liver_n_enh"] = g.abc_liver_n_enh.fillna(0).astype(int)
    return g, gsw


def gtex_counts():
    """Per-gene count of liver cis variants at p<1e-5 and min p (streamed)."""
    nsig, minp = {}, {}
    with gzip.open(INPUTS["gtex_liver_eqtl"], "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        ip, ig = header.index("pvalue"), header.index("gene_id")
        for line in f:
            parts = line.split("\t", ig + 1)
            p = float(parts[ip]) if parts[ip] not in ("NA", "") else 1.0
            gid = parts[ig].split(".")[0]
            if p < minp.get(gid, 2.0):
                minp[gid] = p
            if p < 1e-5:
                nsig[gid] = nsig.get(gid, 0) + 1
    out = pd.DataFrame({"g": list(minp), "gtex_liver_min_p": list(minp.values())})
    out["gtex_liver_n_p1e5"] = out.g.map(nsig).fillna(0).astype(int)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.chdir(ROOT)
    os.makedirs(args.out, exist_ok=False)
    man = pd.DataFrame([(k, v, sha256(v)) for k, v in INPUTS.items()], columns=["input_id", "path", "sha256"])
    man.to_csv(f"{args.out}/input_manifest.tsv", sep="\t", index=False)

    u = build_units()
    u.to_csv(f"{args.out}/units.tsv", sep="\t", index=False)
    g, gsw = build_genes()
    g = g.merge(gtex_counts(), on="g", how="left")
    g.to_csv(f"{args.out}/genes.tsv", sep="\t", index=False)
    gsw.to_csv(f"{args.out}/gene_system_weights.tsv.gz", sep="\t", index=False)

    tested = g[g.cont_z.notna() & (g.n_systems > 0)]
    checks = {
        "units_n": int(len(u)),
        "units_by_cohort": u.dataset.value_counts().to_dict(),
        "units_missing_continuum": int(u.continuum_z.isna().sum()),
        "units_missing_bayesprism": int(u.filter(like="bp_").isna().any(axis=1).sum()),
        "system_columns": int(u.filter(like="sys_").shape[1]),
        "genes_rows": int(len(g)),
        "continuum_associated": int(g.continuum_associated.fillna(False).astype(bool).sum()),
        "canonical_deg": int(g.deg.fillna(False).astype(bool).sum()),
        "coloc_gen_all": int(g.coloc_gen.sum()),
        "h57_ueval_N_G_D_O": [int(g.h57_ueval.sum()), int((g.h57_ueval & g.coloc_gen).sum()),
                              int((g.h57_ueval & g.deg.fillna(False)).sum()),
                              int((g.h57_ueval & g.coloc_gen & g.deg.fillna(False)).sum())],
        "multiplicity_tested_genes": [int(len(tested)), int((tested.n_systems == 1).sum()),
                                      int((tested.n_systems >= 10).sum()), int(tested.n_systems.max())],
        "coverage": {c: float(g[c].notna().mean()) for c in ["loeuf", "gene_length", "gtex_liver_min_p"]},
    }
    with open(f"{args.out}/checks.json", "w") as f:
        json.dump(checks, f, indent=1, default=str)
    print(json.dumps(checks, indent=1, default=str))
    if checks["h57_ueval_N_G_D_O"] != [3769, 426, 220, 26]:
        print("WARNING: H57 universe not reproduced", file=sys.stderr)


if __name__ == "__main__":
    main()

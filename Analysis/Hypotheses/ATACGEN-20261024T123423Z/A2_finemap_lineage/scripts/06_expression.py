#!/usr/bin/env python3
"""A2 step 6: donor-level lineage expression for COLOC r2 genes (Q2 expression lineage).

Unit: biological donor. Runs are collapsed to donors with the donor_pairing.csv files used
by the program label-vs-content analysis; samples without a pairing record are their own
donor. Donor-lineage pseudobulks with fewer than 20 nuclei are dropped.
Lineages are mapped to the five ATAC lineages plus endothelial (no ATAC profile).
Writes results/q2_expression_lineage.tsv.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
A2 = ROOT / "Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage"
SC = ROOT / "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
PRIMARY = ROOT / ("Analysis/Multimodal_Program_Projection/candidates/"
                  "atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/"
                  "genetic_lineage_context_primary_pairs.tsv")
SEED, NBOOT, MIN_NUCLEI, STABLE = 20261024, 1000, 20, 0.90
LINEAGE_FILES = {
    "hepatocyte": ["Hepatocytes"],
    "stellate": ["Fibroblasts"],
    "macrophage": ["Macrophages", "Mono_mono_derived_cells"],
    "cholangiocyte": ["Cholangiocytes"],
    "t_nk": ["T_cells", "Circulating_NK_NKT", "Resident_NK"],
    "endothelial": ["Endothelial_cells"],
}
FILE_TO_LABEL = {"Hepatocytes": "Hepatocytes", "Fibroblasts": "Fibroblasts",
                 "Macrophages": "Macrophages", "Mono_mono_derived_cells": "Mono+mono derived cells",
                 "Cholangiocytes": "Cholangiocytes", "T_cells": "T cells",
                 "Circulating_NK_NKT": "Circulating NK/NKT", "Resident_NK": "Resident NK",
                 "Endothelial_cells": "Endothelial cells"}
LIN = list(LINEAGE_FILES)

meta = pd.read_csv(SC / "sample_celltype_metadata.tsv", sep="\t", dtype=str)
meta["n_cells"] = pd.to_numeric(meta.n_cells)
srr2donor = {}
for ds in ("GSE244832", "GSE185477", "GSE202379", "GSE136103"):
    dp = pd.read_csv(ROOT / f"data/{ds}/metadata/donor_pairing.csv", dtype=str)
    for _, r in dp.iterrows():
        for s in str(r["rna_srrs"]).split(";"):
            if s and s != "nan":
                srr2donor[s.strip()] = f"{ds}_{r['donor_id']}"
meta["donor"] = meta["sample"].map(srr2donor).fillna(meta.dataset + "_" + meta["sample"])

prim = pd.read_csv(PRIMARY, sep="\t", usecols=["gene", "ensembl"]).drop_duplicates()
genes = sorted(prim.gene.unique())

donor_means, donor_mats, ndon, raw_counts = {}, {}, {}, {}
for lin, files in LINEAGE_FILES.items():
    acc = None
    for f in files:
        m = pd.read_csv(SC / f"pseudobulk_{f}.tsv.gz", sep="\t", index_col=0)
        m = m.groupby(level=0).sum()
        acc = m if acc is None else acc.add(m, fill_value=0)
    cols = pd.Series(acc.columns)
    donors = cols.map(srr2donor)
    md = meta[meta.cell_type.isin([FILE_TO_LABEL[f] for f in files])]
    ds_of = md.drop_duplicates("sample").set_index("sample").dataset
    donors = donors.fillna(cols.map(ds_of) + "_" + cols)
    assert donors.notna().all(), lin
    acc = acc.T.groupby(donors.values).sum().T
    ncell = md.groupby("donor").n_cells.sum()
    keep = [d for d in acc.columns if ncell.get(d, 0) >= MIN_NUCLEI]
    acc = acc[keep]
    raw_counts[lin] = acc
    cpm = acc.div(acc.sum(axis=0), axis=1) * 1e6
    lg = np.log2(cpm + 1)
    donor_mats[lin] = lg.reindex(genes)
    donor_means[lin] = lg.mean(axis=1)
    ndon[lin] = len(keep)
    print(lin, "donors", len(keep), flush=True)

P = pd.DataFrame(donor_means)

# Amendment 1 sensitivity: median-of-ratios size factors across all donor-lineage libraries,
# to check that hepatocyte CPM compression by very abundant transcripts does not drive calls.
allc = pd.concat({l: raw_counts[l] for l in LIN}, axis=1)
pos = (allc > 0).all(axis=1)
lg_ref = np.log(allc[pos]).mean(axis=1)
sf = np.exp((np.log(allc[pos]).sub(lg_ref, axis=0)).median(axis=0))
libmed = np.median(allc.sum(axis=0) / sf)
mor = np.log2(allc.div(sf, axis=1) * 1e6 / libmed + 1)
P_mor = pd.DataFrame({l: mor[l].mean(axis=1) for l in LIN})
print("median-of-ratios reference genes", int(pos.sum()), flush=True)
rng = np.random.default_rng(SEED)
arg_counts = pd.DataFrame(0, index=genes, columns=LIN)
arg5_counts = pd.DataFrame(0, index=genes, columns=LIN[:5])
present = [g for g in genes if g in P.index]
boot = {}
for lin in LIN:
    X = donor_mats[lin].loc[present].to_numpy(np.float64).T  # donors x genes
    n = X.shape[0]
    idx = rng.integers(0, n, size=(NBOOT, n))
    boot[lin] = X[idx].mean(axis=1)  # NBOOT x genes
B = np.stack([boot[l] for l in LIN], axis=2)  # NBOOT x genes x lineages
am = B.argmax(axis=2)
am5 = B[:, :, :5].argmax(axis=2)
for j, lin in enumerate(LIN):
    arg_counts.loc[present, lin] = (am == j).sum(axis=0)
for j, lin in enumerate(LIN[:5]):
    arg5_counts.loc[present, lin] = (am5 == j).sum(axis=0)

out = P.reindex(genes).add_prefix("log2cpm_")
out.insert(0, "gene", genes)
out["in_matrix"] = out.gene.isin(P.index)
sub = P.reindex(genes)
out["expr_argmax6"] = sub.idxmax(axis=1)
out["expr_argmax5"] = sub[LIN[:5]].idxmax(axis=1)
out["stab6"] = [arg_counts.loc[g, a] / NBOOT if isinstance(a, str) else np.nan
                for g, a in zip(genes, out.expr_argmax6)]
out["stab5"] = [arg5_counts.loc[g, a] / NBOOT if isinstance(a, str) else np.nan
                for g, a in zip(genes, out.expr_argmax5)]
srt = np.sort(sub[LIN].to_numpy(), axis=1)
out["margin6_log2"] = srt[:, -1] - srt[:, -2]
lin_cpm = 2 ** sub - 1
out["hep_share6"] = lin_cpm["hepatocyte"] / lin_cpm.sum(axis=1)
out["expr_call6"] = np.where(~out.in_matrix, "untestable_absent",
                     np.where(out.stab6 >= STABLE, out.expr_argmax6, "broad"))
out["expr_call5"] = np.where(~out.in_matrix, "untestable_absent",
                     np.where(out.stab5 >= STABLE, out.expr_argmax5, "broad"))
out["expr_hep_class"] = np.select(
    [~out.in_matrix, out.expr_call6 == "hepatocyte",
     out.expr_call6.isin(["stellate", "macrophage", "cholangiocyte", "t_nk", "endothelial"])],
    ["untestable_absent", "hepatocyte_dominant", "non_hepatocyte_dominant"], "broad")
sm = P_mor.reindex(genes)
for lin in LIN:
    out[f"mor_log2_{lin}"] = sm[lin].to_numpy()
out["expr_argmax6_mor"] = sm.idxmax(axis=1).to_numpy()
for lin in LIN:
    out[f"n_donors_{lin}"] = ndon[lin]
(A2 / "results").mkdir(exist_ok=True)
out.to_csv(A2 / "results/q2_expression_lineage.tsv", sep="\t", index=False)
print(out.expr_hep_class.value_counts().to_dict())
print(out.expr_call6.value_counts().to_dict())
print("argmax6 CPM vs median-of-ratios agreement",
      float((out.expr_argmax6 == out.expr_argmax6_mor)[out.in_matrix].mean()))
print(pd.crosstab(out.expr_argmax6, out.expr_argmax6_mor))

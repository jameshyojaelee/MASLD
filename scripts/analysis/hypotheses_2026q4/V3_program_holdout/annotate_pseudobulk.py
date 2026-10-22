"""V3 step 1: outcome-blind lineage assignment and donor pseudobulk for GSE256398.

Follows PRESPEC_V3.md. Reads the 17 MASLD-spectrum CellBender h5 files, applies
nucleus QC, clusters, assigns clusters to lineages by prespecified marker
scores (no program gene is a marker), and writes donor x gene pseudobulk UMI
sums for hepatocytes, stellate/fibroblasts and cholangiocytes.

Usage: python -I annotate_pseudobulk.py <RESULT_DIR>
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

SEED = 0
MARKERS = {
    "hepatocyte": ["ALB", "APOA1", "APOC3", "CYP2E1", "CYP3A4", "SERPINA1", "TF", "HP", "FGB", "APOH"],
    "cholangiocyte": ["KRT7", "KRT19", "SOX9", "EPCAM", "CFTR", "ANXA4"],
    "stellate_fibroblast": ["COL1A1", "COL1A2", "COL3A1", "PDGFRB", "LRAT", "RELN", "RGS5"],
    "endothelial": ["PECAM1", "VWF", "STAB2", "CLEC4G", "FLT1"],
    "macrophage": ["CD68", "CD163", "MARCO", "VSIG4", "C1QA"],
    "t_nk": ["CD3E", "CD2", "NKG7", "PTPRC"],
    "b_plasma": ["MS4A1", "JCHAIN", "MZB1"],
}
PSEUDOBULK_LINEAGES = ["hepatocyte", "stellate_fibroblast", "cholangiocyte"]
TITLE_STAGE = {
    "Healthy Control": ("Healthy", 0),
    "MASLD F0": ("MASLD_F0", 1),
    "MASH Fibrosis F2-3": ("MASH_F2_F3", 2),
    "MASH Fibrosis F3": ("MASH_F2_F3", 2),
    "MASH Cirrhosis": ("MASH_cirrhosis", 3),
}


def parse_soft(path: Path) -> pd.DataFrame:
    rows, cur = [], None
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                cur = {"gsm": line.split("=", 1)[1].strip()}
                rows.append(cur)
            elif cur is not None and line.startswith("!Sample_title"):
                cur["title"] = line.split("=", 1)[1].strip()
            elif cur is not None and line.startswith("!Sample_characteristics_ch1"):
                k, _, v = line.split("=", 1)[1].strip().partition(":")
                cur[k.strip().lower()] = v.strip()
    return pd.DataFrame(rows)


def main() -> None:
    res = Path(sys.argv[1])
    dl = res / "downloads_untrusted"
    out = res / "pseudobulk"
    if out.exists():
        raise FileExistsError(out)
    if not (res / "PRESPEC_V3.md").exists():
        raise FileNotFoundError("PRESPEC_V3.md must exist before analysis")
    out.mkdir()

    soft = parse_soft(dl / "GSE256398_family.soft.gz")
    files = sorted(dl.glob("GSM*_CB_raw_feature_bc_matrix_filtered.h5"))
    if len(files) != 17:
        raise RuntimeError(f"expected 17 h5 files, found {len(files)}")

    meta_rows, adatas = [], []
    for f in files:
        gsm, sample = f.name.split("_")[:2]
        rec = soft.loc[soft["gsm"] == gsm]
        if len(rec) != 1:
            raise RuntimeError(f"{gsm} not matched to one GEO record")
        rec = rec.iloc[0]
        m = re.match(r"^(S\d+), Human, (.+)$", rec["title"])
        if m is None or m.group(1) != sample or m.group(2) not in TITLE_STAGE:
            raise RuntimeError(f"unexpected title for {gsm}: {rec['title']}")
        stage, ordinal = TITLE_STAGE[m.group(2)]
        age = float(re.match(r"(\d+)", rec["age"]).group(1))
        sex = rec["sex"].strip().lower()
        meta_rows.append(dict(donor=sample, gsm=gsm, title=rec["title"], stage=stage,
                              stage_ordinal=ordinal, disease=int(stage != "Healthy"),
                              sex=sex, age=age))
        a = sc.read_10x_h5(f, gex_only=True)
        a.var_names_make_unique()
        a.obs_names = [f"{sample}_{b}" for b in a.obs_names]
        a.obs["donor"] = sample
        a.var["mt"] = a.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
        n0 = a.n_obs
        a = a[(a.obs["n_genes_by_counts"] >= 200) & (a.obs["pct_counts_mt"] < 5)].copy()
        meta_rows[-1].update(n_nuclei_raw=n0, n_nuclei_qc=a.n_obs)
        adatas.append(a)

    meta = pd.DataFrame(meta_rows)
    meta.to_csv(out / "donor_metadata.tsv", sep="\t", index=False)

    gene_ids = adatas[0].var["gene_ids"].copy()
    for a in adatas[1:]:
        if not a.var_names.equals(adatas[0].var_names):
            raise RuntimeError("feature sets differ between samples")
    adata = ad.concat(adatas, join="inner", merge="same")
    adata.var["gene_ids"] = gene_ids.loc[adata.var_names].values
    del adatas
    adata.X = sp.csr_matrix(adata.X)
    # CellBender output is integer-valued; keep the raw layer for pseudobulk.
    adata.layers["counts"] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                                batch_key="donor", layer="counts")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, n_comps=30, mask_var="highly_variable", random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=15, random_state=SEED)
    sc.tl.leiden(adata, resolution=1.0, random_state=SEED, flavor="igraph",
                 n_iterations=2, directed=False)

    present = {k: [g for g in v if g in adata.var_names] for k, v in MARKERS.items()}
    for k, genes in present.items():
        sc.tl.score_genes(adata, genes, score_name=f"mk_{k}", random_state=SEED)
    cols = [f"mk_{k}" for k in MARKERS]
    cl = adata.obs.groupby("leiden", observed=True)[cols].mean()
    cl["lineage"] = cl[cols].idxmax(axis=1).str.replace("mk_", "", regex=False)
    cl["n_nuclei"] = adata.obs["leiden"].value_counts().reindex(cl.index).values
    cl.to_csv(out / "cluster_lineage_assignment.tsv", sep="\t")
    adata.obs["lineage"] = adata.obs["leiden"].map(cl["lineage"]).astype(str)
    pd.Series({k: ";".join(v) for k, v in present.items()}).to_csv(
        out / "markers_used.tsv", sep="\t", header=False)

    counts_tab = pd.crosstab(adata.obs["donor"], adata.obs["lineage"])
    counts_tab.to_csv(out / "nuclei_per_donor_lineage.tsv", sep="\t")

    X = adata.layers["counts"]
    donors = meta["donor"].tolist()
    for lin in PSEUDOBULK_LINEAGES:
        cols_out = {}
        for d in donors:
            mask = ((adata.obs["donor"] == d) & (adata.obs["lineage"] == lin)).values
            cols_out[d] = np.asarray(X[mask].sum(axis=0)).ravel()
        pb = pd.DataFrame(cols_out, index=adata.var["gene_ids"].values)
        pb.index.name = "gene_id"
        pb.insert(0, "symbol", adata.var_names.values)
        pb.to_csv(out / f"pseudobulk_{lin}.tsv.gz", sep="\t", compression="gzip")
    print(counts_tab)
    print(cl[["lineage", "n_nuclei"]])


if __name__ == "__main__":
    main()

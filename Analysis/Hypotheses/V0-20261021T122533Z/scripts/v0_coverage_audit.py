"""V0: rerun V3 labelling and pseudobulk, check reproduction, then compute
program-gene detection in hepatocyte (and secondary-lineage) nuclei.

Follows PRESPEC_V0.md. No program score or test is computed.

Usage: python -I v0_coverage_audit.py <V0_RUN_DIR> [OUT_SUBDIR, default out]
"""
from __future__ import annotations

import csv
import importlib.util
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
V3_SCRIPT = ROOT / "scripts/analysis/hypotheses_2026q4/V3_program_holdout/annotate_pseudobulk.py"
V3_PSEUDOBULK = ROOT / "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk/pseudobulk_hepatocyte.tsv.gz"
MEMBERSHIP = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
PROGRAMS = {
    "hotspot_hepatocytes_48f39dd4d817a10e": ("Ductular injury (BICC1)", ["hepatocyte", "cholangiocyte"]),
    "hotspot_hepatocytes_f05c535ae5bbc0b9": ("Stromal ECM (IGFBP7)", ["hepatocyte", "stellate_fibroblast"]),
}
DETECTION_MIN = 0.01
MIN_GENES, MIN_WEIGHT = 8, 0.80


def load_v3_module():
    spec = importlib.util.spec_from_file_location("v3_annotate", V3_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_labelled(v3, dl: Path):
    """V3 steps, verbatim order and parameters (annotate_pseudobulk.main)."""
    soft = v3.parse_soft(dl / "GSE256398_family.soft.gz")
    files = sorted(dl.glob("GSM*_CB_raw_feature_bc_matrix_filtered.h5"))
    if len(files) != 17:
        raise RuntimeError(f"expected 17 h5 files, found {len(files)}")
    donors, adatas = [], []
    for f in files:
        gsm, sample = f.name.split("_")[:2]
        rec = soft.loc[soft["gsm"] == gsm].iloc[0]
        m = re.match(r"^(S\d+), Human, (.+)$", rec["title"])
        if m is None or m.group(1) != sample or m.group(2) not in v3.TITLE_STAGE:
            raise RuntimeError(f"unexpected title for {gsm}: {rec['title']}")
        donors.append(sample)
        a = sc.read_10x_h5(f, gex_only=True)
        a.var_names_make_unique()
        a.obs_names = [f"{sample}_{b}" for b in a.obs_names]
        a.obs["donor"] = sample
        a.var["mt"] = a.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
        a = a[(a.obs["n_genes_by_counts"] >= 200) & (a.obs["pct_counts_mt"] < 5)].copy()
        adatas.append(a)
    gene_ids = adatas[0].var["gene_ids"].copy()
    adata = ad.concat(adatas, join="inner", merge="same")
    adata.var["gene_ids"] = gene_ids.loc[adata.var_names].values
    del adatas
    adata.X = sp.csr_matrix(adata.X)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                                batch_key="donor", layer="counts")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, n_comps=30, mask_var="highly_variable", random_state=v3.SEED)
    sc.pp.neighbors(adata, n_neighbors=15, random_state=v3.SEED)
    sc.tl.leiden(adata, resolution=1.0, random_state=v3.SEED, flavor="igraph",
                 n_iterations=2, directed=False)
    present = {k: [g for g in v if g in adata.var_names] for k, v in v3.MARKERS.items()}
    for k, genes in present.items():
        sc.tl.score_genes(adata, genes, score_name=f"mk_{k}", random_state=v3.SEED)
    cols = [f"mk_{k}" for k in v3.MARKERS]
    cl = adata.obs.groupby("leiden", observed=True)[cols].mean()
    cl["lineage"] = cl[cols].idxmax(axis=1).str.replace("mk_", "", regex=False)
    adata.obs["lineage"] = adata.obs["leiden"].map(cl["lineage"]).astype(str)
    return adata, donors


def main() -> None:
    run = Path(sys.argv[1])
    out = run / (sys.argv[2] if len(sys.argv) > 2 else "out")
    out.mkdir(exist_ok=True)
    if any(out.iterdir()):
        raise FileExistsError(f"{out} is not empty")
    if not (run / "PRESPEC_V0.md").exists():
        raise FileNotFoundError("PRESPEC_V0.md must exist before analysis")
    v3 = load_v3_module()
    adata, donors = build_labelled(v3, run / "downloads_untrusted")
    X = adata.layers["counts"]

    # Reproduction check against V3's hepatocyte pseudobulk.
    mask_hep = (adata.obs["lineage"] == "hepatocyte").values
    new = {}
    for d in donors:
        m = mask_hep & (adata.obs["donor"] == d).values
        new[d] = np.asarray(X[m].sum(axis=0)).ravel()
    new_pb = pd.DataFrame(new, index=adata.var["gene_ids"].values)
    old_pb = pd.read_csv(V3_PSEUDOBULK, sep="\t", index_col=0).drop(columns="symbol")
    same_shape = old_pb.shape == new_pb.shape and list(old_pb.columns) == list(new_pb.columns)
    n_diff = int((old_pb.values != new_pb.loc[old_pb.index].values).sum()) if same_shape else -1
    repro = pd.DataFrame([{"v3_shape": str(old_pb.shape), "v0_shape": str(new_pb.shape),
                           "same_shape_and_donors": same_shape, "n_cells_differing": n_diff,
                           "reproduced": same_shape and n_diff == 0}])
    repro.to_csv(out / "v0_reproduction.tsv", sep="\t", index=False)
    pd.crosstab(adata.obs["donor"], adata.obs["lineage"]).to_csv(out / "nuclei_per_donor_lineage.tsv", sep="\t")
    if not bool(repro["reproduced"].iloc[0]):
        print(repro.to_string())
        raise SystemExit("STOP: V3 hepatocyte pseudobulk not reproduced; no coverage reported")

    members = pd.read_csv(MEMBERSHIP, sep="\t")
    members = members[members["mapped_symbol_status"] == "gencode_v49_unique_symbol_confirmed"]
    gene_rows, cov_rows = [], []
    for uid, (label, lineages) in PROGRAMS.items():
        mm = members[members["program_uid"] == uid]
        for lin in lineages:
            lmask = (adata.obs["lineage"] == lin).values
            Xl = X[lmask]
            n_nuc = int(lmask.sum())
            det_pool = np.asarray((Xl > 0).mean(axis=0)).ravel()
            donor_l = adata.obs["donor"].values[lmask]
            present_w = h90_w = 0.0
            n_present = n_h90 = 0
            for _, r in mm.iterrows():
                g, w = r["mapped_symbol"], float(r["original_l1_weight"])
                if g in adata.var_names:
                    j = adata.var_names.get_loc(g)
                    dp = float(det_pool[j])
                    per_donor = [float((Xl[donor_l == d][:, j] > 0).mean()) for d in donors if (donor_l == d).any()]
                    med = float(np.median(per_donor)) if per_donor else float("nan")
                    present_w += w; n_present += 1
                    if dp >= DETECTION_MIN:
                        h90_w += w; n_h90 += 1
                else:
                    dp, med = float("nan"), float("nan")
                gene_rows.append({"program": label, "lineage": lin, "gene": g, "l1_weight": w,
                                  "on_axis": g in adata.var_names, "detection_pooled": dp,
                                  "detection_donor_median": med, "passes_h90_rule": dp >= DETECTION_MIN})
            cov_rows.append({"program": label, "lineage": lin, "n_nuclei": n_nuc,
                             "presence_genes": n_present, "presence_weight": round(present_w, 4),
                             "h90_rule_genes": n_h90, "h90_rule_weight": round(h90_w, 4),
                             "testable_presence": n_present >= MIN_GENES and present_w >= MIN_WEIGHT,
                             "testable_h90_rule": n_h90 >= MIN_GENES and h90_w >= MIN_WEIGHT})
    pd.DataFrame(gene_rows).to_csv(out / "v0_program_gene_detection.tsv", sep="\t", index=False)
    pd.DataFrame(cov_rows).to_csv(out / "v0_coverage.tsv", sep="\t", index=False)
    print(pd.DataFrame(cov_rows).to_string())


if __name__ == "__main__":
    main()

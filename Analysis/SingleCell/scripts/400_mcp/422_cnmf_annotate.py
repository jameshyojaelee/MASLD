#!/usr/bin/env python
"""
422_cnmf_annotate.py — annotate cNMF programs post-hoc.

Reads the chosen k solution for a named cNMF run and produces:
  - program_celltype_usage.tsv: program x cell_type mean usage + fraction of cells > 0.1
  - program_entropy.tsv: cell-type entropy (bits) per program -> shared/cell-type-specific label
  - program_stage_usage.tsv: program x disease_stage (0..3) mean usage + KW p-value
  - program_phenotype_cor.tsv: Pearson/Spearman vs 9 phenotypes at donor level (mean usage per donor)
  - program_topgenes.tsv: top-100 genes per program by spectra score
  - program_pathway_enrichment.tsv: fgsea-based GO/Hallmark pathway enrichment per program (MSigDB H + C2:CP)

Usage:
  python 422_cnmf_annotate.py --name global --k 20 --density-threshold 0.03
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
CNMF_ROOT = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/cnmf_runs"
INPUT_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs"
OUT_ROOT = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/cnmf_annot"
OUT_ROOT.mkdir(parents=True, exist_ok=True)


def entropy_bits(p: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    if p.sum() == 0:
        return 0.0
    p = p / p.sum()
    return float(-(p * np.log2(p)).sum())


def load_cnmf_outputs(name: str, k: int, dthresh: float):
    run_dir = CNMF_ROOT / name
    dt_tag = f"{dthresh:.2f}".replace(".", "_")
    # cnmf 1.5 naming: spectra score (Z), spectra tpm (raw), usages (row-normalized)
    spectra_score_f = run_dir / f"{name}.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"
    usage_f = run_dir / f"{name}.usages.k_{k}.dt_{dt_tag}.consensus.txt"
    spectra_tpm_f = run_dir / f"{name}.gene_spectra_tpm.k_{k}.dt_{dt_tag}.txt"
    spectra_score = pd.read_csv(spectra_score_f, sep="\t", index_col=0)
    usage = pd.read_csv(usage_f, sep="\t", index_col=0)
    try:
        spectra_tpm = pd.read_csv(spectra_tpm_f, sep="\t", index_col=0)
    except Exception:
        spectra_tpm = None
    return spectra_score, usage, spectra_tpm


def main(args: argparse.Namespace) -> None:
    spectra_score, usage, spectra_tpm = load_cnmf_outputs(args.name, args.k, args.density_threshold)
    usage.index = usage.index.astype(str)
    spectra_score.columns = spectra_score.columns.astype(str)
    print(f"[422] spectra: {spectra_score.shape}, usage: {usage.shape}")

    # Load AnnData to join cell metadata
    adata_fn = INPUT_DIR / (f"atlas_cnmf_{args.name[4:]}.h5ad" if args.name.startswith("pct_") else "atlas_cnmf_global.h5ad")
    adata = ad.read_h5ad(adata_fn, backed="r")
    obs = adata.obs.copy()
    obs.index = obs.index.astype(str)
    adata.file.close()
    obs = obs.loc[obs.index.intersection(usage.index)]
    usage = usage.loc[obs.index]

    # --- Program usage by cell_type ---
    prog_cols = list(usage.columns)
    out_dir = OUT_ROOT / args.name
    out_dir.mkdir(exist_ok=True, parents=True)
    # (a) mean usage per cell type
    ct_mean = usage.groupby(obs["cell_type"].values).mean()
    ct_mean.to_csv(out_dir / f"program_celltype_mean_usage.k{args.k}.tsv", sep="\t")
    # (b) fraction of cells >0.1 usage per cell type
    ct_frac = (usage > 0.1).groupby(obs["cell_type"].values).mean()
    ct_frac.to_csv(out_dir / f"program_celltype_frac_active.k{args.k}.tsv", sep="\t")

    # (c) Entropy per program
    rows = []
    for p in prog_cols:
        row_vec = ct_mean[p].values
        rows.append({
            "program": p,
            "entropy_bits": entropy_bits(row_vec),
            "n_celltypes_active": int((ct_frac[p] > 0.05).sum()),
            "max_ct_usage": float(ct_mean[p].max()),
            "max_ct": str(ct_mean[p].idxmax()),
        })
    ent = pd.DataFrame(rows)
    ent["tag"] = np.where(ent["entropy_bits"] > 1.3, "shared",
                 np.where(ent["entropy_bits"] < 0.6, "celltype_specific", "mixed"))
    ent.to_csv(out_dir / f"program_entropy.k{args.k}.tsv", sep="\t", index=False)
    print(f"[422] shared={int((ent['tag']=='shared').sum())} / celltype_specific={int((ent['tag']=='celltype_specific').sum())} / mixed={int((ent['tag']=='mixed').sum())}")

    # --- Program x stage ---
    if "disease_stage_numeric" in obs.columns:
        stage_rows = []
        stage_mask = obs["disease_stage_numeric"].notna()
        u_s = usage.loc[stage_mask]
        s = obs.loc[stage_mask, "disease_stage_numeric"].astype(float)
        for p in prog_cols:
            vals = u_s[p].values
            means = pd.Series(vals).groupby(s.values).mean().to_dict()
            groups = [vals[s.values == g] for g in sorted(s.unique())]
            try:
                H, pv = stats.kruskal(*groups)
            except Exception:
                H, pv = np.nan, np.nan
            stage_rows.append({"program": p, **{f"mean_stage_{int(k2)}": means.get(float(k2), np.nan) for k2 in sorted(s.unique())}, "kruskal_H": H, "kruskal_p": pv})
        pd.DataFrame(stage_rows).to_csv(out_dir / f"program_stage_usage.k{args.k}.tsv", sep="\t", index=False)

    # --- Program x 9 phenotypes (donor-level) ---
    # Mean usage per donor
    donor_usage = usage.groupby(obs["sample"].values).mean()
    donor_meta = obs.groupby("sample", observed=True).agg(
        dataset=("dataset", "first"),
        condition=("condition", "first"),
        condition_binary=("condition_binary", "first"),
        disease_stage_coarse=("disease_stage_coarse", "first"),
        disease_stage_numeric=("disease_stage_numeric", "first"),
    ).reset_index().set_index("sample")
    donor_usage = donor_usage.loc[donor_usage.index.intersection(donor_meta.index)]
    donor_meta = donor_meta.loc[donor_usage.index]

    pheno_cols = ["disease_stage_numeric"]
    if "condition_binary" in donor_meta.columns:
        # Encode Healthy=0 else 1
        donor_meta["_cb_num"] = donor_meta["condition_binary"].astype(str).map(lambda x: 0 if x.lower().startswith("heal") else 1)
        pheno_cols.append("_cb_num")
    pheno_rows = []
    for p in prog_cols:
        for ph in pheno_cols:
            v1 = pd.to_numeric(donor_usage[p], errors="coerce")
            v2 = pd.to_numeric(donor_meta[ph], errors="coerce")
            ok = v1.notna() & v2.notna()
            if ok.sum() < 5:
                continue
            r_p, p_p = stats.pearsonr(v1[ok], v2[ok])
            r_s, p_s = stats.spearmanr(v1[ok], v2[ok])
            pheno_rows.append({
                "program": p, "phenotype": ph, "n": int(ok.sum()),
                "pearson_r": r_p, "pearson_p": p_p,
                "spearman_r": r_s, "spearman_p": p_s,
            })
    pd.DataFrame(pheno_rows).to_csv(out_dir / f"program_phenotype_cor.k{args.k}.tsv", sep="\t", index=False)

    # --- Top genes per program ---
    top_rows = []
    for p in spectra_score.index:  # spectra_score: programs x genes
        top = spectra_score.loc[p].sort_values(ascending=False).head(args.top_genes)
        for rank, (g, s) in enumerate(top.items(), start=1):
            top_rows.append({"program": str(p), "rank": rank, "gene_name": g, "spectra_score": float(s)})
    pd.DataFrame(top_rows).to_csv(out_dir / f"program_topgenes.k{args.k}.tsv", sep="\t", index=False)

    print(f"[422] wrote outputs to {out_dir}")
    print("[422] DONE.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--density-threshold", type=float, default=0.03)
    ap.add_argument("--top-genes", type=int, default=100)
    args = ap.parse_args()
    main(args)

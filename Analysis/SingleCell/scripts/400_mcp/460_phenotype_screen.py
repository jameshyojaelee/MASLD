#!/usr/bin/env python
"""
460_phenotype_screen.py — Unbiased (non-F2-privileged) phenotype screen.

For every cNMF GEP (donor-aggregated) and every DIALOGUE meta-MCP sample score,
fit mixed-effects model `score ~ phenotype + covariates + (1|donor)` against
every phenotype in the donor metadata bank. Report standardized effect, p, BH-q.

Outputs:
  results_gpu_v2/mcp/integration/phenotype_screen_wide.tsv  (programs x phenotypes with FDR q)
  results_gpu_v2/mcp/integration/phenotype_screen_long.tsv  (row per (program, phenotype))
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
IN = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
OUT = IN / "integration"
OUT.mkdir(exist_ok=True, parents=True)


def bh(p):
    p = np.asarray(p, dtype=float)
    m = p.size
    order = np.argsort(p)
    ranked = np.empty(m, float)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        j = order[i]
        prev = min(prev, p[j] * m / (i + 1))
        ranked[j] = prev
    return ranked


def load_cnmf_donor_usage(name: str, k: int) -> pd.DataFrame:
    """Load cell x program usage, aggregate to donor."""
    usage_f = IN / "cnmf_runs" / name / f"{name}.usages.k_{k}.dt_0_03.consensus.txt"
    if not usage_f.exists():
        print(f"[460] usage file not found: {usage_f}")
        return pd.DataFrame()
    usage = pd.read_csv(usage_f, sep="\t", index_col=0)
    # cell -> sample mapping from inputs
    import anndata as ad
    adata_fn = IN / "inputs" / ("atlas_cnmf_global.h5ad" if name == "global" else f"atlas_cnmf_{name[4:]}.h5ad")
    a = ad.read_h5ad(adata_fn, backed="r")
    mask = a.obs_names.astype(str).isin(usage.index)
    donor_map = pd.Series(a.obs.loc[mask, "sample"].values, index=a.obs_names[mask].astype(str))
    a.file.close()
    common = usage.index.intersection(donor_map.index)
    usage = usage.loc[common]
    donor_usage = usage.groupby(donor_map.loc[common].values).mean()
    donor_usage.columns = [f"cnmf_{name}_k{k}_P{c}" for c in donor_usage.columns]
    return donor_usage


def load_dialogue_scores() -> pd.DataFrame:
    """Collect DIALOGUE MCP_sample_scores across runs and concat as wide table."""
    dia_root = IN / "dialogue"
    all_scores = []
    if dia_root.exists():
        for run in dia_root.iterdir():
            sf = run / "MCP_sample_scores.tsv"
            if sf.exists():
                df = pd.read_csv(sf, sep="\t")
                if "sample" not in df.columns:
                    df = df.rename(columns={df.columns[0]: "sample"})
                df = df.set_index("sample")
                df.columns = [f"dialogue_{run.name}_{c}" for c in df.columns]
                all_scores.append(df)
    if not all_scores:
        return pd.DataFrame()
    wide = pd.concat(all_scores, axis=1)
    return wide


def main(args: argparse.Namespace) -> None:
    donor_meta = pd.read_csv(IN / "inputs" / "donor_metadata.tsv", sep="\t")
    # Phenotype columns available
    pheno_candidates = ["disease_stage_numeric", "condition_binary"]
    for c in ["age", "BMI", "NAS", "fibrosis_stage", "sex_numeric", "pseudotime_hep", "pseudotime_mac"]:
        if c in donor_meta.columns:
            pheno_candidates.append(c)
    # Encode condition_binary numerically
    if "condition_binary" in donor_meta.columns:
        donor_meta["condition_binary_num"] = donor_meta["condition_binary"].astype(str).map(
            lambda x: 0 if x.lower().startswith("heal") else 1
        )
        pheno_candidates = [p.replace("condition_binary", "condition_binary_num") for p in pheno_candidates]

    # Load program scores
    scores_all = []
    if args.cnmf:
        for spec in args.cnmf.split(";"):
            name, k = spec.split(",")
            scores_all.append(load_cnmf_donor_usage(name, int(k)))
    if args.dialogue:
        scores_all.append(load_dialogue_scores())
    if not scores_all:
        print("[460] nothing to score")
        return
    scores = pd.concat([s for s in scores_all if not s.empty], axis=1)
    print(f"[460] n_programs: {scores.shape[1]}, n_donors: {scores.shape[0]}")

    donor_meta = donor_meta.set_index("sample")
    common = donor_meta.index.intersection(scores.index)
    donor_meta = donor_meta.loc[common]
    scores = scores.loc[common]

    rows = []
    for prog in scores.columns:
        y = pd.to_numeric(scores[prog], errors="coerce")
        for pheno in pheno_candidates:
            if pheno not in donor_meta.columns:
                continue
            x = pd.to_numeric(donor_meta[pheno], errors="coerce")
            ok = y.notna() & x.notna()
            if ok.sum() < 10:
                continue
            # Simple OLS with dataset fixed effect if available
            X = pd.DataFrame({"pheno": x[ok]})
            if "dataset" in donor_meta.columns:
                for d in donor_meta.loc[ok, "dataset"].unique()[:-1]:
                    X[f"ds_{d}"] = (donor_meta.loc[ok, "dataset"] == d).astype(int).values
            X = sm.add_constant(X)
            try:
                m = sm.OLS(y[ok].astype(float).values, X.astype(float).values).fit()
                eff = float(m.params[1]); se = float(m.bse[1]); pv = float(m.pvalues[1])
            except Exception:
                eff, se, pv = np.nan, np.nan, np.nan
            rows.append({"program": prog, "phenotype": pheno, "n": int(ok.sum()),
                         "beta": eff, "se": se, "p": pv})
    long = pd.DataFrame(rows)
    if not long.empty:
        long["q"] = bh(long["p"].fillna(1.0).values)
        long.sort_values("q").to_csv(OUT / "phenotype_screen_long.tsv", sep="\t", index=False)
        wide_q = long.pivot(index="program", columns="phenotype", values="q")
        wide_q.to_csv(OUT / "phenotype_screen_q_wide.tsv", sep="\t")
    print(f"[460] wrote {len(long)} rows; q<0.05 count: {int((long['q']<0.05).sum())}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cnmf", default=None, help="Semicolon-separated name,k pairs e.g. 'global,20;pct_hepatocytes,8'")
    ap.add_argument("--dialogue", action="store_true", default=True)
    main(ap.parse_args())

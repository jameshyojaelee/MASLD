#!/usr/bin/env python
"""
460_phenotype_screen_clean.py — SENSITIVITY VARIANT of 460_phenotype_screen.py.

Excludes donors with `exclude_stage_analysis == True` (protocol-contaminated
GSE136103 + dubious Healthy donors) before fitting program ~ phenotype models.

Inputs (added vs. 460):
  Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv
  (joined onto donor_metadata.tsv by `sample` to attach exclude_stage_analysis)

Outputs:
  results_gpu_v2/mcp/integration/phenotype_screen_long_clean.tsv
  results_gpu_v2/mcp/integration/phenotype_screen_q_wide_clean.tsv
  results_gpu_v2/mcp/phenotype_screen_clean.tsv             (mirror of long for top-level access)
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

EXT_META = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"


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
        print(f"[460c] usage file not found: {usage_f}")
        return pd.DataFrame()
    usage = pd.read_csv(usage_f, sep="\t", index_col=0)
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

    # CLEAN: join extended metadata + drop excluded donors
    ext = pd.read_csv(EXT_META, sep="\t")
    keep_cols = ["sample", "exclude_stage_analysis", "protocol_group"]
    have = [c for c in keep_cols if c in ext.columns]
    donor_meta = donor_meta.merge(ext[have], on="sample", how="left")
    n_before = len(donor_meta)
    excluded_mask = donor_meta["exclude_stage_analysis"].fillna(False).astype(bool)
    n_excluded = int(excluded_mask.sum())
    donor_meta = donor_meta.loc[~excluded_mask].copy()
    print(f"[460c] CLEAN sensitivity: excluded {n_excluded}/{n_before} donors with exclude_stage_analysis==True")
    print(f"[460c] remaining donors: {len(donor_meta)}")

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
        print("[460c] nothing to score")
        return
    scores = pd.concat([s for s in scores_all if not s.empty], axis=1)
    print(f"[460c] n_programs: {scores.shape[1]}, n_donors (pre-exclusion): {scores.shape[0]}")

    donor_meta = donor_meta.set_index("sample")
    common = donor_meta.index.intersection(scores.index)
    donor_meta = donor_meta.loc[common]
    scores = scores.loc[common]
    print(f"[460c] n_donors after intersection: {len(donor_meta)}")

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
        long.sort_values("q").to_csv(OUT / "phenotype_screen_long_clean.tsv", sep="\t", index=False)
        wide_q = long.pivot(index="program", columns="phenotype", values="q")
        wide_q.to_csv(OUT / "phenotype_screen_q_wide_clean.tsv", sep="\t")
        # Also mirror to top-level for protocol-sensitivity bundle
        long.sort_values("q").to_csv(IN / "phenotype_screen_clean.tsv", sep="\t", index=False)
    print(f"[460c] wrote {len(long)} rows; q<0.05 count: {int((long['q']<0.05).sum())}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cnmf", default=None, help="Semicolon-separated name,k pairs e.g. 'global,16'")
    ap.add_argument("--dialogue", action="store_true", default=True)
    main(ap.parse_args())

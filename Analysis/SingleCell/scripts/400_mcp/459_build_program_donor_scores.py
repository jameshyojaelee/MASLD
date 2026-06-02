#!/usr/bin/env python
"""
459_build_program_donor_scores.py — Build consolidated donor x program score matrix.

Combines:
  - cNMF global donor usage (mean per donor per program)
  - cNMF per-cell-type donor usage
  - DIALOGUE MCP sample scores (all runs concatenated)

Outputs:
  program_donor_scores_wide.tsv.gz   (donors x programs wide matrix)
  program_donor_scores.tsv.gz        (long: sample, program, score)

Used downstream by 460 phenotype screen, 461 changepoint, 463 composition controls.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
OUT = MCP / "integration"
OUT.mkdir(parents=True, exist_ok=True)


def load_cnmf_donor_usage(name: str, k: int, dthresh: float = 0.03) -> pd.DataFrame:
    dt = f"{dthresh:.2f}".replace(".", "_")
    uf = MCP / "cnmf_runs" / name / f"{name}.usages.k_{k}.dt_{dt}.consensus.txt"
    if not uf.exists():
        print(f"[459] cNMF usage missing for {name} k={k}")
        return pd.DataFrame()
    usage = pd.read_csv(uf, sep="\t", index_col=0)
    ct_suffix = "" if name == "global" else name[4:]
    adata_fn = MCP / "inputs" / ("atlas_cnmf_global.h5ad" if name == "global" else f"atlas_cnmf_{ct_suffix}.h5ad")
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
    dia_root = MCP / "dialogue"
    all_scores = []
    if not dia_root.exists():
        return pd.DataFrame()
    for run in dia_root.iterdir():
        if not run.is_dir():
            continue
        for sf in run.glob("MCP_sample_scores_*.tsv"):
            ct = sf.stem.replace("MCP_sample_scores_", "")
            df = pd.read_csv(sf, sep="\t")
            if "sample" not in df.columns:
                df = df.rename(columns={df.columns[0]: "sample"})
            df = df.set_index("sample")
            df.columns = [f"dialogue_{run.name}_{ct}_{c}" for c in df.columns]
            all_scores.append(df)
    return pd.concat(all_scores, axis=1) if all_scores else pd.DataFrame()


def main(args) -> None:
    frames = []
    for spec in (args.cnmf or "").split(";"):
        if not spec.strip():
            continue
        name, k = spec.split(",")
        u = load_cnmf_donor_usage(name, int(k))
        if not u.empty:
            frames.append(u)
    d = load_dialogue_scores()
    if not d.empty:
        frames.append(d)
    if not frames:
        print("[459] nothing to consolidate")
        return
    wide = pd.concat(frames, axis=1)
    wide.index.name = "sample"
    print(f"[459] wide shape: {wide.shape}")
    wide.to_csv(OUT / "program_donor_scores_wide.tsv.gz", sep="\t", compression="gzip")
    long = wide.reset_index().melt(id_vars="sample", var_name="program", value_name="score").dropna()
    long.to_csv(OUT / "program_donor_scores.tsv.gz", sep="\t", compression="gzip", index=False)
    print(f"[459] wrote {long.shape[0]} long rows.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cnmf", default="global,20",
                    help="Semicolon-separated name,k pairs, e.g. 'global,20;pct_hepatocytes,8'")
    main(ap.parse_args())

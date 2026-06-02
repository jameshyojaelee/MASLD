#!/usr/bin/env python
"""
figS_scdrs_magma_aggregate.py
=============================
Combine per-GWAS MAGMA-anchored scDRS group-analysis outputs into a single
table for figure-building. Joins per-trait CT z-scores with the audit
table (GWAS family, ancestry, N).

Outputs:
  Analysis/SingleCell/results_gpu_v2/disease_signatures/
    scdrs_magma_celltype_enrichment.csv  (long format: gwas, family, CT, z, p, fdr)
    scdrs_magma_summary.csv              (one row per GWAS with n_ct_sig, top CT)
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import pandas as pd

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SIG_DIR    = BASE / "Analysis/SingleCell/results_gpu_v2/disease_signatures"
OUT_DIR    = SIG_DIR / "scdrs_magma_outputs_25pct"
N_TABLE    = BASE / "tools/magma/work/gwas_inputs/_N_table.tsv"
GS_SUMMARY = SIG_DIR / "MAGMA_anchored_scdrs_summary.csv"
OUT_LONG   = SIG_DIR / "scdrs_magma_celltype_enrichment.csv"
OUT_SUM    = SIG_DIR / "scdrs_magma_summary.csv"


def log(msg: str) -> None:
    print(f"[magma_aggregate] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


def main() -> None:
    log(f"loading audit {N_TABLE}")
    audit = pd.read_csv(N_TABLE, sep="\t")
    audit_ok = audit[audit["status"] == "OK"].copy()
    log(f"  {len(audit_ok)} GWAS expected")

    log(f"loading gene-set summary {GS_SUMMARY}")
    if GS_SUMMARY.exists():
        gs_sum = pd.read_csv(GS_SUMMARY)
    else:
        gs_sum = pd.DataFrame()

    log("reading per-GWAS CT enrichment files")
    parts = []
    summary_rows = []
    for _, row in audit_ok.iterrows():
        gwas_id = row["gwas_id"]
        family  = row["family"]
        N       = row["N"]
        f = OUT_DIR / f"{gwas_id}_celltype_enrichment.csv"
        if not f.exists() or os.path.getsize(f) < 100:
            log(f"  MISSING/EMPTY {gwas_id}")
            summary_rows.append({"gwas_id": gwas_id, "family": family, "N": N,
                                 "status": "missing", "n_ct_fdr05": 0,
                                 "top_ct": None, "top_z": None, "top_p": None})
            continue
        df = pd.read_csv(f)
        if "assoc_mcz" not in df.columns or "cell_type" not in df.columns:
            log(f"  bad columns {gwas_id}: {list(df.columns)}")
            continue
        df = df.copy()
        df["gwas_id"] = gwas_id
        df["family"]  = family
        df["N"]       = N
        parts.append(df)

        n_sig = int((df["assoc_mcp"] < 0.05).sum())
        top   = df.loc[df["assoc_mcz"].idxmax()] if not df.empty else None
        summary_rows.append({
            "gwas_id": gwas_id, "family": family, "N": N, "status": "OK",
            "n_ct_fdr05": n_sig,
            "top_ct": top["cell_type"] if top is not None else None,
            "top_z":  float(top["assoc_mcz"]) if top is not None else None,
            "top_p":  float(top["assoc_mcp"]) if top is not None else None,
        })
        log(f"  {gwas_id}: {len(df)} CTs;  n_ct p<0.05 = {n_sig};  top {top['cell_type']} z={top['assoc_mcz']:.2f}")

    if parts:
        long_df = pd.concat(parts, ignore_index=True)
        long_df.to_csv(OUT_LONG, index=False)
        log(f"wrote {OUT_LONG}  shape={long_df.shape}")
    else:
        log("no per-GWAS CT enrichment files found")
        return

    pd.DataFrame(summary_rows).to_csv(OUT_SUM, index=False)
    log(f"wrote {OUT_SUM}")
    log("done.")


if __name__ == "__main__":
    main()

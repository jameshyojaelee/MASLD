#!/usr/bin/env python
"""
magma_gwas_prep.py
==================
For each GWAS sumstats in GWAS/finemapping/data/sumstats/ that has a documented
sample size, build a MAGMA-compatible input file with columns SNP, P, N by
joining chr:pos → rsID against the 1000G EUR reference. Skips GWAS without a
documented N.

Outputs (one per GWAS):
  tools/magma/work/gwas_inputs/{gwas_id}.magma_input.tsv
  tools/magma/work/gwas_inputs/_N_table.tsv  (audit: gwas_id, N_total, n_snps_in)
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import pandas as pd

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SUMSTATS_DIR = BASE / "GWAS/finemapping/data/sumstats"
LOOKUP_PATH  = BASE / "tools/magma/work/chr_pos_to_rsid.tsv"
OUT_DIR      = BASE / "tools/magma/work/gwas_inputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
N_TABLE_PATH = OUT_DIR / "_N_table.tsv"

# (gwas_id, sumstats filename stem (no _reformatted_hg19.tsv), N_total, family)
# Documented sample sizes from docs/manuscript/working/numbers_cheatsheet.md
# and Pan-UKBB README. BS/M/F PanUKBB strata are skipped (undocumented).
GWAS_TABLE = [
    # gwas_id,             sumstats_stem,                       N,          family
    ("ukbb_alt",           "UKBB_ALT",                          344000,     "Liver enzymes (EUR)"),
    ("ukbb_ast",           "UKBB_AST",                          344000,     "Liver enzymes (EUR)"),
    ("ukbb_ggt",           "UKBB_GGT",                          344000,     "Liver enzymes (EUR)"),
    ("bbj_alt",            "BBJ_ALT",                           261406,     "Liver enzymes (EAS)"),
    ("bbj_ast",            "BBJ_AST",                           261270,     "Liver enzymes (EAS)"),
    ("bbj_ggt",            "BBJ_GGT",                           164006,     "Liver enzymes (EAS)"),
    ("panukbb_afr_alt",    "PanUKBB_AFR_ALT",                   6636,       "Liver enzymes (AFR)"),
    ("panukbb_afr_ast",    "PanUKBB_AFR_AST",                   6636,       "Liver enzymes (AFR)"),
    ("panukbb_afr_ggt",    "PanUKBB_AFR_GGT",                   6636,       "Liver enzymes (AFR)"),
    ("panukbb_csa_alt",    "PanUKBB_CSA_ALT",                   8876,       "Liver enzymes (CSA)"),
    ("panukbb_csa_ast",    "PanUKBB_CSA_AST",                   8876,       "Liver enzymes (CSA)"),
    ("panukbb_csa_ggt",    "PanUKBB_CSA_GGT",                   8876,       "Liver enzymes (CSA)"),
    ("anstee_nafld",       "Anstee_NAFLD",                      19264,      "NAFLD"),
    ("finngen_nafld",      "FinnGen_NAFLD",                     438857,     "NAFLD"),
    ("finngen_nash",       "FinnGen_NASH",                      340000,     "NASH"),
    # finngen_hcc / ghouse_hcc / ghouse_cirrhosis dropped 2026-06-06
    # (cirrhosis/HCC GWAS removed from canonical 23-GWAS COLOC portfolio)
    ("pazoki_pdff",        "Pazoki_PDFF",                       32858,      "PDFF (imaging)"),
]


def log(msg: str) -> None:
    print(f"[magma_gwas_prep] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


def main() -> None:
    log(f"loading chr:pos → rsID lookup from {LOOKUP_PATH}")
    lookup = pd.read_csv(LOOKUP_PATH, sep="\t", dtype=str)
    log(f"  lookup rows: {len(lookup):,}")
    lookup_set = dict(zip(lookup["chr_pos"], lookup["rsid"]))

    audit_rows = []
    for gwas_id, stem, N, family in GWAS_TABLE:
        ss_path = SUMSTATS_DIR / f"{stem}_reformatted_hg19.tsv"
        if not ss_path.exists():
            log(f"  SKIP {gwas_id}: missing {ss_path}")
            audit_rows.append({"gwas_id": gwas_id, "stem": stem, "N": N,
                               "family": family, "ss_path": str(ss_path),
                               "n_snps_in": 0, "n_snps_matched": 0,
                               "status": "missing_sumstats"})
            continue
        log(f"reading sumstats: {gwas_id}  ({ss_path.name})")
        t0 = time.time()
        ss = pd.read_csv(ss_path, sep="\t",
                         usecols=["chromosome", "position", "pval"],
                         dtype={"chromosome": str, "position": str, "pval": float})
        n_in = len(ss)
        log(f"  rows in: {n_in:,}  ({time.time()-t0:.1f}s)")

        # Build chr:pos key
        ss["chr_pos"] = ss["chromosome"].str.removeprefix("chr") + ":" + ss["position"]
        ss["SNP"] = ss["chr_pos"].map(lookup_set)
        n_matched = ss["SNP"].notna().sum()
        log(f"  rows joined to rsID: {n_matched:,}  ({n_matched/n_in*100:.1f}%)")

        out = ss.dropna(subset=["SNP", "pval"]).copy()
        out = out[(out["pval"] > 0) & (out["pval"] <= 1)]
        out["N"] = N
        out_path = OUT_DIR / f"{gwas_id}.magma_input.tsv"
        out[["SNP", "pval", "N"]].rename(columns={"pval": "P"}).to_csv(
            out_path, sep="\t", index=False)
        log(f"  wrote {out_path}  ({len(out):,} rows)")
        audit_rows.append({"gwas_id": gwas_id, "stem": stem, "N": N,
                           "family": family, "ss_path": str(ss_path),
                           "n_snps_in": n_in, "n_snps_matched": int(n_matched),
                           "status": "OK"})

    pd.DataFrame(audit_rows).to_csv(N_TABLE_PATH, sep="\t", index=False)
    log(f"wrote audit table {N_TABLE_PATH}  ({len(audit_rows)} rows)")
    log("done.")


if __name__ == "__main__":
    main()

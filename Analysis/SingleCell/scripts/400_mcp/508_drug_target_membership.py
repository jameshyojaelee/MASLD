#!/usr/bin/env python
"""
508_drug_target_membership.py — Echo Task 20 / A18.

For each named MASLD-pipeline drug target, report which cNMF programs (k=16) carry it
in their top-100 spectra-score gene list, and the rank within the best program.

Output: validation/drug/drug_target_program_membership.tsv
"""
from __future__ import annotations
import os
from pathlib import Path

import pandas as pd

ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
TOPG_F = MCP / "cnmf_annot/global/program_topgenes.k16.tsv"
OUT_DIR = MCP / "validation/drug"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_F = OUT_DIR / "drug_target_program_membership.tsv"

TARGETS = [
    "THRB", "NR1H4", "GLP1R", "ACACA", "SCD1", "DGAT2", "FASN",
    "PNPLA3", "TM6SF2", "MBOAT7", "HSD17B13", "GCGR", "FGFR1",
    "FGF21", "SLC10A1", "FXR", "PPARA", "PPARG",
    # SCD1 = SCD; FXR = NR1H4; keep aliases for searching robustness:
    "SCD",
]


def main():
    topg = pd.read_csv(TOPG_F, sep="\t")
    # Restrict to top-100 per program by rank (table is already sorted by spectra_score
    # within program, so use rank if available else nlargest)
    if "rank" in topg.columns:
        topg = topg[topg["rank"] <= 100].copy()
    else:
        topg = topg.groupby("program", group_keys=False).apply(
            lambda g: g.nlargest(100, "spectra_score")
        )

    rows = []
    for tgt in TARGETS:
        hits = topg[topg["gene_name"] == tgt]
        if len(hits) == 0:
            rows.append({
                "target_gene": tgt,
                "program_ids_containing": "",
                "n_programs_containing": 0,
                "best_program": "",
                "top_rank_in_best_program": "",
                "best_spectra_score": "",
            })
            continue
        progs = sorted(hits["program"].astype(str).unique().tolist(), key=lambda s: int(s) if s.isdigit() else s)
        best_row = hits.sort_values("spectra_score", ascending=False).iloc[0]
        rows.append({
            "target_gene": tgt,
            "program_ids_containing": ",".join(progs),
            "n_programs_containing": int(len(progs)),
            "best_program": str(best_row["program"]),
            "top_rank_in_best_program": int(best_row["rank"]) if "rank" in best_row.index else "",
            "best_spectra_score": float(best_row["spectra_score"]),
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT_F, sep="\t", index=False)
    print(f"[508] wrote {OUT_F}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()

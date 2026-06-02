#!/usr/bin/env python
"""
508b_drug_target_membership_v2.py — Echo Task 20 / A18 redo.

v1 (508_drug_target_membership.py) used the pre-computed top-100 table, which:
  - missed targets that land at rank 101-500 (still meaningful);
  - searched exact symbols only, missing HGNC alias mismatches (SCD1 vs SCD);
  - did not expand to target FAMILIES (ACC = ACACA+ACACB, PPAR = A/G/D, ...).

v2 loads the full `gene_spectra_score` matrix (16 programs × 36,972 genes), ranks
each gene within each program, and queries expanded target families with aliases.

Output: validation/drug/drug_target_program_membership_v2.tsv
Columns: target_family, gene_symbol, in_atlas, program, rank_in_program,
         spectra_score, in_top_100, in_top_300, n_genes_total_in_program.

Summary columns (aggregated) are written as `*_summary.tsv`.
"""
from __future__ import annotations
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
SPECTRA = MCP / "cnmf_runs/global/global.gene_spectra_score.k_16.dt_0_03.txt"
OUT_DIR = MCP / "validation/drug"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_F = OUT_DIR / "drug_target_program_membership_v2.tsv"
SUMMARY_F = OUT_DIR / "drug_target_program_membership_v2_summary.tsv"

# Adversarial-reviewer-specified families (plus aliases). `SCD` is the canonical HGNC
# symbol; `SCD1` is a legacy alias. Include both so aliases fall through gracefully.
TARGET_FAMILIES: dict[str, list[str]] = {
    "THRB":        ["THRB", "THRA"],
    "FXR":         ["NR1H4", "NR1H3"],
    "GLP1R":       ["GLP1R", "GIPR", "GCGR"],           # tri-agonist landscape
    "ACC":         ["ACACA", "ACACB"],                   # acetyl-CoA carboxylase α/β
    "SCD":         ["SCD", "SCD5", "SCD1"],              # SCD canonical; SCD1 alias
    "DGAT":        ["DGAT1", "DGAT2"],
    "FASN":        ["FASN"],
    "LIPID_GWAS":  ["PNPLA3", "TM6SF2", "MBOAT7", "HSD17B13", "PNPLA2"],
    "FGFR":        ["FGFR1", "FGFR2", "FGFR3", "FGFR4"],
    "FGF":         ["FGF21", "FGF19", "KLB"],            # FGF21/19 + β-Klotho coreceptor
    "PPAR":        ["PPARA", "PPARG", "PPARD"],
    "BILE":        ["NR1H4", "SLC10A1", "SLCO1B1", "SLC10A2"],  # FXR + NTCP + OATP1B1 + ASBT
    "STAT":        ["STAT3"],                            # JAK/STAT inflammation
}


def main() -> None:
    print(f"[508b] loading spectra matrix: {SPECTRA}")
    # File layout: program IDs as rows (index_col=0 = first column "1".."16"),
    # gene symbols as columns (36,972 of them).
    spec = pd.read_csv(SPECTRA, sep="\t", index_col=0)
    print(f"[508b] spectra shape: {spec.shape}  (programs × genes)")
    print(f"[508b] program IDs: {list(spec.index)}")
    print(f"[508b] example columns: {list(spec.columns[:5])}")

    # Transpose so rows are genes, columns are programs — friendlier for rank-per-program.
    sT = spec.T  # (genes, programs)
    atlas_genes = set(sT.index)

    # Precompute per-program ranks (1 = highest spectra_score within that program).
    ranks = sT.rank(ascending=False, method="min")  # (genes, programs)
    n_genes_total = sT.shape[0]

    rows = []
    for family, symbols in TARGET_FAMILIES.items():
        for sym in symbols:
            in_atlas = sym in atlas_genes
            if not in_atlas:
                rows.append({
                    "target_family": family,
                    "gene_symbol": sym,
                    "in_atlas": False,
                    "program": "",
                    "rank_in_program": "",
                    "spectra_score": "",
                    "in_top_100": "",
                    "in_top_300": "",
                    "n_genes_total_in_program": n_genes_total,
                })
                continue
            # Emit one row per program for this gene
            for prog in sT.columns:
                r = int(ranks.loc[sym, prog])
                s = float(sT.loc[sym, prog])
                rows.append({
                    "target_family": family,
                    "gene_symbol": sym,
                    "in_atlas": True,
                    "program": str(prog),
                    "rank_in_program": r,
                    "spectra_score": s,
                    "in_top_100": bool(r <= 100),
                    "in_top_300": bool(r <= 300),
                    "n_genes_total_in_program": n_genes_total,
                })

    df = pd.DataFrame(rows)
    # Sort: in-atlas first, then by best rank per gene, then by family
    df["_best_rank"] = df.groupby("gene_symbol")["rank_in_program"].transform(
        lambda x: pd.to_numeric(x, errors="coerce").min()
    )
    df = df.sort_values(["_best_rank", "gene_symbol", "rank_in_program"]).drop(columns="_best_rank")
    df.to_csv(OUT_F, sep="\t", index=False)
    print(f"[508b] wrote full table: {OUT_F}  ({len(df)} rows)")

    # ===== Summary =====
    # Per-gene: best rank + best program
    hits = df[df["in_atlas"] == True].copy()
    hits["rank_in_program"] = pd.to_numeric(hits["rank_in_program"])

    per_gene = hits.sort_values("rank_in_program").groupby("gene_symbol", as_index=False).first()
    per_gene = per_gene[["gene_symbol", "program", "rank_in_program", "spectra_score"]].rename(
        columns={"program": "best_program", "rank_in_program": "best_rank",
                 "spectra_score": "best_score"}
    )
    # Also: any_top_100, any_top_300
    per_gene["in_top_100_anywhere"] = per_gene["best_rank"] <= 100
    per_gene["in_top_300_anywhere"] = per_gene["best_rank"] <= 300
    # Merge back with family
    fam_map = df[["gene_symbol", "target_family", "in_atlas"]].drop_duplicates(subset="gene_symbol")
    per_gene = fam_map.merge(per_gene, on="gene_symbol", how="right")
    per_gene.to_csv(SUMMARY_F, sep="\t", index=False)
    print(f"[508b] wrote per-gene summary: {SUMMARY_F}")

    # ===== Headline stats =====
    n_symbols_total = sum(len(v) for v in TARGET_FAMILIES.values())
    n_symbols_unique = len(set(sym for v in TARGET_FAMILIES.values() for sym in v))
    n_in_atlas = per_gene["gene_symbol"].nunique()
    n_in_top_100 = per_gene[per_gene["in_top_100_anywhere"]]["gene_symbol"].nunique()
    n_in_top_300 = per_gene[per_gene["in_top_300_anywhere"]]["gene_symbol"].nunique()

    family_has_top_100 = set()
    family_has_top_300 = set()
    for _, r in df.iterrows():
        if r["in_atlas"] is True and r["in_top_100"] is True:
            family_has_top_100.add(r["target_family"])
        if r["in_atlas"] is True and r["in_top_300"] is True:
            family_has_top_300.add(r["target_family"])

    print("\n========== HEADLINE ==========")
    print(f"Unique target symbols queried: {n_symbols_unique}")
    print(f"  -> found in atlas (spectra matrix): {n_in_atlas}")
    print(f"  -> in top-100 of any program: {n_in_top_100}")
    print(f"  -> in top-300 of any program: {n_in_top_300}")
    print(f"Families with any member in top-100: {len(family_has_top_100)}/{len(TARGET_FAMILIES)} "
          f"-> {sorted(family_has_top_100)}")
    print(f"Families with any member in top-300: {len(family_has_top_300)}/{len(TARGET_FAMILIES)} "
          f"-> {sorted(family_has_top_300)}")
    print("\nPer-gene best rank (atlas members only):")
    show = per_gene[per_gene["in_atlas"] == True].sort_values("best_rank")
    print(show.to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Saunders 181-gene panel × MASLD evidence overlap.

Feeds the wet-lab Perturb-seq pilot proposal (task #44). For each of Saunders'
181 1:1 human-orth panel genes, annotate MASLD-relevance across:
  - multi-evidence atlas (coloc_best_susie_pp4, is_conserved, drug-target cols)
  - convergence_evidence.csv (tier, rank, score)
  - progression_driver_genetics.csv (driver, genetically_validated)

Output:
  data/perturbation/datasets/saunders2025/processed/saunders_panel_vs_masld_overlap.csv
  data/perturbation/datasets/saunders2025/processed/saunders_panel_overlap_summary.md
"""
from __future__ import annotations

import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_CSV = PROJECT_ROOT / "data/perturbation/datasets/saunders2025/processed/saunders_panel_vs_masld_overlap.csv"
OUT_MD = PROJECT_ROOT / "data/perturbation/datasets/saunders2025/processed/saunders_panel_overlap_summary.md"

ATLAS = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
CONVERGENCE = PROJECT_ROOT / "RNA-seq/results/multi_evidence/convergence_evidence.csv"
PROGRESSION = PROJECT_ROOT / "RNA-seq/results/stratified_causal/progression_driver_genetics.csv"
ORTH_MAP = PROJECT_ROOT / "data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz"

# Canonical MASLD drug-target panel (literature curation; n=22 used in
# Phase L drug validation / backtest framework)
CANONICAL_DRUG_TARGETS = {
    "THRB",       # resmetirom (Madrigal)
    "PNPLA3",     # AZD2693 / AZX
    "HSD17B13",   # AZD9202
    "GLP1R",      # semaglutide
    "GIPR",       # tirzepatide co-target
    "SCD",        # aramchol
    "NR1H4",      # OCA, cilofexor (FXR)
    "FGFR1",      # efruxifermin (FGF21 analog)
    "FGF21",      # efruxifermin
    "ACACA",      # firsocostat (Acetyl-CoA carboxylase α)
    "ACACB",      # firsocostat (Acetyl-CoA carboxylase β)
    "DGAT2",      # vupanorsen / PF-06865571 / IONIS-DGAT2
    "GCKR",       # protective GCKR.P446L variant
    "MBOAT7",     # MASLD risk variant
    "TM6SF2",     # MASLD risk variant
    "MTARC1",     # MASLD risk variant
    "APOB",       # mipomersen
    "CCR2",       # cenicriviroc
    "CCR5",       # cenicriviroc
    "TGFB1",      # pirfenidone-style anti-fibrotic
    "MTOR",       # rapalogs (also in Saunders panel)
    "RPTOR",     # mTORC1 component
}


def load_saunders_panel(orth_map_path):
    """Saunders mouse target genes (sgrna table) → 1:1 human orthologs."""
    sgrna_tsv = PROJECT_ROOT / "data/perturbation/datasets/saunders2025/processed/per_cell_sgrna.tsv.gz"
    df = pd.read_csv(sgrna_tsv, sep="\t", compression="gzip", low_memory=False)
    mouse_targets = sorted({g for g in df.loc[df["tier"]=="confident", "target_gene"].unique() if g != "control"})
    print(f"Saunders confident-tier unique mouse targets: {len(mouse_targets)}")

    orth = pd.read_csv(orth_map_path, sep="\t").dropna(subset=["mouse_gene_symbol","human_gene_symbol"])
    one2one = orth[orth["ortholog_type"] == "ortholog_one2one"]
    m2h = dict(zip(one2one["mouse_gene_symbol"], one2one["human_gene_symbol"]))
    rows = []
    for g in mouse_targets:
        h = m2h.get(g)
        if h:
            rows.append({"gene_mouse": g, "gene_human": h})
    panel = pd.DataFrame(rows).drop_duplicates(subset=["gene_human"])
    print(f"  → {len(panel)} unique 1:1 human orthologs (full Saunders panel)")
    return panel


def annotate_with_atlas(panel, atlas_path):
    cols = ["human_symbol", "is_conserved", "coloc_best_susie_pp4_polyfun", "coloc_best_susie_gwas_polyfun",
            "coloc_pp4", "dream_logFC", "dream_padj",
            "dgidb_druggable", "opentargets_drug", "drug_target_stratification"]
    print(f"\nLoading atlas (selected cols)...")
    atlas = pd.read_csv(atlas_path, usecols=cols, low_memory=False)
    atlas = atlas.rename(columns={
        "coloc_best_susie_pp4_polyfun": "coloc_best_susie_pp4",
        "coloc_best_susie_gwas_polyfun": "coloc_best_susie_gwas",
    })
    print(f"  atlas rows: {len(atlas):,}")
    return panel.merge(atlas, left_on="gene_human", right_on="human_symbol", how="left").drop(columns=["human_symbol"])


def annotate_with_convergence(panel, conv_path):
    print(f"Loading convergence evidence...")
    conv = pd.read_csv(conv_path, low_memory=False)
    keep = ["human_symbol", "convergence_score", "convergence_rank", "tier", "concordance_state"]
    return panel.merge(conv[keep], left_on="gene_human", right_on="human_symbol", how="left").drop(columns=["human_symbol"])


def annotate_with_progression(panel, prog_path):
    print(f"Loading progression drivers...")
    prog = pd.read_csv(prog_path, low_memory=False)
    # collapse per-gene (multiple transitions/rows) — take MIN driver_rank, MAX driver_score
    agg = prog.groupby("gene_symbol").agg(
        prog_driver_score_max=("driver_score", "max"),
        prog_driver_rank_best=("driver_rank", "min"),
        prog_genetically_validated=("genetically_validated", lambda x: any(bool(v) for v in x)),
        prog_transitions=("transition", lambda x: "|".join(sorted(set(map(str, x))))),
    ).reset_index().rename(columns={"gene_symbol":"_g"})
    panel = panel.merge(agg, left_on="gene_human", right_on="_g", how="left").drop(columns=["_g"])
    panel["is_progression_driver"] = panel["prog_driver_rank_best"].notna()
    return panel


def add_drug_target_flag(panel):
    panel["is_drug_target"] = panel["gene_human"].isin(CANONICAL_DRUG_TARGETS)
    # Combine atlas drug-target hints as 'has_existing_drug_evidence'
    panel["has_drug_evidence"] = (
        panel["is_drug_target"]
        | panel["dgidb_druggable"].astype(str).str.lower().isin(["true","1","yes"])
        | panel["opentargets_drug"].notna() & (panel["opentargets_drug"].astype(str).str.strip() != "")
    )
    return panel


def add_recommendation(panel):
    # Tier strings can be e.g. "1_Genetic_validated", "2_Strong" etc.
    def is_top_tier(t):
        if not isinstance(t, str):
            return False
        return t.startswith("1_") or t.startswith("2_") or t.startswith("3_")
    panel["_top_tier"] = panel["tier"].apply(is_top_tier)
    panel["_strong_coloc"] = panel["coloc_best_susie_pp4"].fillna(0).astype(float) > 0.5

    def rec(row):
        if row["is_drug_target"]:
            return "keep_as_positive_control"
        if row["_top_tier"] or row["_strong_coloc"] or row["is_progression_driver"]:
            return "keep_in_pilot"
        return "drop_for_panel"

    panel["recommendation"] = panel.apply(rec, axis=1)
    return panel.drop(columns=["_top_tier","_strong_coloc"])


def write_summary(panel, out_md):
    n = len(panel)
    counts = panel["recommendation"].value_counts().to_dict()
    n_keep = counts.get("keep_in_pilot", 0)
    n_pos = counts.get("keep_as_positive_control", 0)
    n_drop = counts.get("drop_for_panel", 0)
    n_in_atlas = panel["coloc_best_susie_pp4"].notna().sum()
    n_pp4_05 = (panel["coloc_best_susie_pp4"].fillna(0).astype(float) > 0.5).sum()
    n_pp4_08 = (panel["coloc_best_susie_pp4"].fillna(0).astype(float) > 0.8).sum()
    n_progression = panel["is_progression_driver"].sum()
    n_conserved = panel["is_conserved"].astype(str).str.lower().isin(["true","1","yes"]).sum()
    n_drug_canonical = panel["is_drug_target"].sum()
    n_drug_any = panel["has_drug_evidence"].sum()

    # Top-10 by MASLD relevance: prioritize by (tier rank, coloc_pp4 desc, convergence_rank asc)
    def tier_int(t):
        if not isinstance(t, str): return 99
        try: return int(t.split("_",1)[0])
        except Exception: return 99
    panel = panel.copy()
    panel["_tier_int"] = panel["tier"].apply(tier_int)
    top10 = panel.sort_values(
        ["_tier_int","coloc_best_susie_pp4","convergence_rank"],
        ascending=[True, False, True],
        na_position="last"
    ).head(10)

    lines = []
    lines.append("# Saunders 2025 Perturb-Multi Panel × MASLD Evidence Overlap")
    lines.append("")
    lines.append(f"Generated from `09_saunders_panel_vs_masld.py`. Saunders panel = 181 1:1 mouse→human orthologs from the confident-tier sgRNA table (the wet-lab-realized panel from GSM8478327).")
    lines.append("")
    lines.append("## Headline counts")
    lines.append("")
    lines.append(f"- **Total Saunders panel (1:1 ortholog space)**: {n}")
    lines.append(f"- **keep_in_pilot**: {n_keep} (MASLD-relevant via convergence tier 1-3, OR coloc_best_susie_pp4 > 0.5, OR progression driver)")
    lines.append(f"- **keep_as_positive_control**: {n_pos} (canonical MASLD drug-target genes — wet-lab calibration positives)")
    lines.append(f"- **drop_for_panel**: {n_drop} (pure essentiality / housekeeping without MASLD evidence)")
    lines.append("")
    lines.append("## Evidence stratification (within Saunders panel)")
    lines.append("")
    lines.append(f"- In multi-evidence atlas: {n_in_atlas}/{n}")
    lines.append(f"- coloc_best_susie_pp4 > 0.5: {n_pp4_05}/{n}")
    lines.append(f"- coloc_best_susie_pp4 > 0.8: {n_pp4_08}/{n}")
    lines.append(f"- In progression-driver table: {n_progression}/{n}")
    lines.append(f"- Conserved_Core: {n_conserved}/{n}")
    lines.append(f"- Canonical MASLD drug targets: {n_drug_canonical}/{n}")
    lines.append(f"- Any drug evidence (DGIdb/OpenTargets/canonical): {n_drug_any}/{n}")
    lines.append("")
    lines.append("## Top-10 overlap genes (by MASLD relevance)")
    lines.append("")
    lines.append("| gene_human | gene_mouse | tier | coloc_pp4 | progression_driver | drug_target | recommendation |")
    lines.append("|---|---|---|---|---|---|---|")
    for _, r in top10.iterrows():
        lines.append(f"| {r['gene_human']} | {r['gene_mouse']} | {r['tier']} | "
                     f"{r['coloc_best_susie_pp4']:.3f}" if pd.notna(r['coloc_best_susie_pp4']) else f"| {r['gene_human']} | {r['gene_mouse']} | {r['tier']} | NA")
        # rewrite as single-line — easier:
    # rebuild table cleanly
    lines = lines[:-len(top10)]
    for _, r in top10.iterrows():
        pp4 = r["coloc_best_susie_pp4"]
        pp4_s = f"{pp4:.3f}" if pd.notna(pp4) and isinstance(pp4, (int, float)) else "NA"
        prog = "Y" if r["is_progression_driver"] else "—"
        drug = "Y" if r["is_drug_target"] else ("evidence" if r["has_drug_evidence"] else "—")
        lines.append(f"| {r['gene_human']} | {r['gene_mouse']} | {r['tier']} | {pp4_s} | {prog} | {drug} | {r['recommendation']} |")
    lines.append("")
    lines.append("## Strategic verdict")
    lines.append("")
    coverage_pct = 100.0 * (n_keep + n_pos) / n
    if (n_keep + n_pos) >= 60:
        lines.append(f"**Reusable.** {n_keep + n_pos}/{n} ({coverage_pct:.0f}%) of Saunders panel is MASLD-relevant or canonical-drug-target. A focused subset of Saunders' existing sgRNA library + chemistry can carry the pilot for the overlapping genes, saving ~$10-15K library cost + 4-week pool generation. For the remaining ~30-40 sgRNAs to reach the 100-sgRNA pilot target, design de novo against the MASLD-specific gene set not in Saunders (top COLOC PP4 > 0.8 + canonical drug targets missing from Saunders).")
    elif (n_keep + n_pos) >= 30:
        lines.append(f"**Partially reusable.** {n_keep + n_pos}/{n} ({coverage_pct:.0f}%) is MASLD-relevant — too thin to anchor the full pilot but worth piggybacking the overlap into a mixed library (Saunders subset + ~50 de novo MASLD-specific sgRNAs).")
    else:
        lines.append(f"**Insufficient overlap.** Only {n_keep + n_pos}/{n} ({coverage_pct:.0f}%) Saunders panel is MASLD-relevant. Design a fresh 100-sgRNA panel from scratch (drop Saunders chemistry reuse); use the panel only as a positive-control sanity reference.")
    lines.append("")
    lines.append("Recommended next step: cross-reference the **MASLD-side gap set** = atlas top COLOC PP4 > 0.5 (~364 genes) MINUS Saunders panel. Those are the candidates for de novo sgRNA design.")
    lines.append("")
    lines.append(f"Data: `{OUT_CSV.relative_to(PROJECT_ROOT)}` ({n} rows). Re-run with `09_saunders_panel_vs_masld.py`.")

    out_md.write_text("\n".join(lines))
    print(f"Wrote {out_md}")


def main():
    panel = load_saunders_panel(ORTH_MAP)
    panel = annotate_with_atlas(panel, ATLAS)
    panel = annotate_with_convergence(panel, CONVERGENCE)
    panel = annotate_with_progression(panel, PROGRESSION)
    panel = add_drug_target_flag(panel)
    panel = add_recommendation(panel)

    print("\nSaunders panel + MASLD evidence:")
    print(f"  total rows: {len(panel)}")
    print(f"  recommendation counts: {panel['recommendation'].value_counts().to_dict()}")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV}")

    write_summary(panel, OUT_MD)


if __name__ == "__main__":
    main()

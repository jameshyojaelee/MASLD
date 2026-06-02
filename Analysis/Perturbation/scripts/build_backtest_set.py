#!/usr/bin/env python3
"""Build canonical Phase 1 backtest validation set.

Integrates four orthogonal evidence sources:
  (a) Conserved_Core cross-species concordant DEGs (n=1,108)
  (b) COLOC PP4_susie > 0.8 causal genes (~314)
  (c) Curated MASLD drug targets (approved + clinical trial + genetic risk)
  (d) Cas13 positive controls (curated steatosis literature panel)

Outputs:
  - backtest_set.csv (deduped by gene; sources concatenated with ';')
  - backtest_set_summary.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR = ROOT / "Analysis/Perturbation/results/validation_sets"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# (a) Conserved_Core: cross-species concordant DEGs
# ---------------------------------------------------------------------------
concordance = pd.read_csv(
    ROOT / "Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv",
    low_memory=False,
)
conserved = concordance.loc[
    concordance["primary_category"] == "Conserved",
    ["human_symbol", "mean_h_lfc", "translatability_tier"],
].dropna(subset=["human_symbol"]).copy()
conserved["source"] = "conserved_core"
conserved["expected_direction"] = conserved["mean_h_lfc"].apply(
    lambda x: "up" if x > 0 else ("down" if x < 0 else "either")
)
conserved["expected_pathway"] = "cross-species_conserved_DEG"
conserved["cell_type_relevance"] = "general"
conserved["evidence_summary"] = conserved["translatability_tier"].apply(
    lambda t: f"Conserved cross-species DEG (translatability={t}); mean human LFC concordant with mouse"
)
conserved = conserved.rename(columns={"human_symbol": "gene"})[
    ["gene", "source", "expected_direction", "expected_pathway",
     "cell_type_relevance", "evidence_summary"]
]
print(f"[a] Conserved_Core: {len(conserved)} genes", file=sys.stderr)

# ---------------------------------------------------------------------------
# (b) COLOC PP4_susie > 0.8
# ---------------------------------------------------------------------------
coloc = pd.read_csv(
    ROOT / "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv",
    low_memory=False,
)
coloc_hi = coloc.loc[coloc["coloc_best_susie_pp4"] > 0.8].copy()
coloc_hi["source"] = "coloc_pp4_high"
coloc_hi["expected_direction"] = "either"  # GWAS direction not predictive of expr direction
coloc_hi["expected_pathway"] = "GWAS_causal_locus"
coloc_hi["cell_type_relevance"] = "general"
coloc_hi["evidence_summary"] = coloc_hi.apply(
    lambda r: f"SuSiE-COLOC PP4={r['coloc_best_susie_pp4']:.3f} in {r['coloc_best_susie_gwas']} (Broadaway liver eQTL × GWAS)",
    axis=1,
)
coloc_hi = coloc_hi.rename(columns={"gene": "gene"})[
    ["gene", "source", "expected_direction", "expected_pathway",
     "cell_type_relevance", "evidence_summary"]
]
print(f"[b] COLOC PP4>0.8: {len(coloc_hi)} genes", file=sys.stderr)

# ---------------------------------------------------------------------------
# (c) Curated known MASLD drug targets
# ---------------------------------------------------------------------------
# Columns: gene, expected_direction (in disease vs healthy), pathway, cell_type, evidence, source_subtype
drug_targets = [
    # --- Approved / Phase 3 ---
    ("THRB",   "down", "thyroid hormone signaling",  "hepatocyte", "Phase 3 approved 2024 (THRB agonist resmetirom; first FDA-approved MASH drug)", "drug_target_approved"),
    ("NR1H4",  "down", "FXR/bile acid signaling",    "hepatocyte", "Obeticholic acid (FXR agonist; FLINT/REGENERATE Phase 3)", "drug_target_trial"),
    ("FGF21",  "down", "FGF21/metabolic homeostasis", "hepatocyte", "Pegozafermin, efruxifermin (FGF21 analogs; Phase 2b/3 anti-fibrotic)", "drug_target_trial"),
    ("FGFR1",  "either", "FGF21 receptor signaling",  "hepatocyte", "Efruxifermin engages FGFR1/KLB complex (Phase 2b)", "drug_target_trial"),
    ("GLP1R",  "either", "incretin signaling",        "multi",      "Semaglutide (NN9931, ESSENCE Phase 3 for MASH)", "drug_target_approved"),
    ("GIPR",   "either", "incretin signaling",        "multi",      "Tirzepatide (dual GIP/GLP1; SYNERGY-NASH Phase 2)", "drug_target_trial"),
    # --- Lipogenesis ---
    ("SCD",    "up",   "stearoyl-CoA desaturation",  "hepatocyte", "Aramchol (SCD1 inhibitor; ARMOR Phase 3)", "drug_target_trial"),
    ("ACACA",  "up",   "fatty acid synthesis",       "hepatocyte", "Firsocostat (ACC1/2 inhibitor; ATLAS Phase 2)", "drug_target_trial"),
    ("ACACB",  "up",   "fatty acid oxidation",       "hepatocyte", "Firsocostat dual ACC1/ACC2 target", "drug_target_trial"),
    ("DGAT2",  "up",   "triglyceride synthesis",     "hepatocyte", "Ervogastat (DGAT2 inhibitor; Pfizer Phase 2)", "drug_target_trial"),
    ("SREBF1", "up",   "lipogenic transcription",    "hepatocyte", "SREBP1c master lipogenic TF; indirect target of many MASH drugs", "drug_target_trial"),
    # --- FXR pathway alt agonists ---
    ("NR0B2",  "either", "FXR-SHP axis",             "hepatocyte", "SHP downstream of FXR; tropifexor/cilofexor pathway", "drug_target_trial"),
    # --- Fibrosis / immune ---
    ("TGFB1",  "up",   "TGF-beta fibrogenesis",      "stellate",   "Anti-TGFβ antibodies (PF-06480605 etc.) in MASH trials", "drug_target_trial"),
    ("CCR2",   "up",   "monocyte recruitment",       "macrophage", "Cenicriviroc (CCR2/CCR5 inhibitor; AURORA Phase 3 — failed but mechanistic anchor)", "drug_target_trial"),
    ("CCR5",   "up",   "macrophage trafficking",     "macrophage", "Cenicriviroc dual CCR2/CCR5 target", "drug_target_trial"),
    # --- Genetic risk variants (well-established MASLD GWAS / WES anchors) ---
    ("PNPLA3", "up",   "lipid droplet remodeling",   "hepatocyte", "rs738409 I148M; single strongest MASLD common-variant risk (multi-ancestry)", "drug_target_genetic_risk"),
    ("HSD17B13","down","retinol metabolism / lipid droplet", "hepatocyte", "rs72613567 protective LoF (Abul-Husn 2018 NEJM)", "drug_target_genetic_risk"),
    ("MTARC1", "down", "mitochondrial amidoxime",    "hepatocyte", "rs2642438 A165T protective (Emdin 2020 PLoS Genet)", "drug_target_genetic_risk"),
    ("TM6SF2", "up",   "VLDL secretion",             "hepatocyte", "rs58542926 E167K risk variant (Kozlitina 2014)", "drug_target_genetic_risk"),
    ("GCKR",   "up",   "glucokinase regulation",     "hepatocyte", "rs1260326 P446L MASLD/triglyceride risk", "drug_target_genetic_risk"),
    ("APOB",   "either", "VLDL assembly",            "hepatocyte", "Rare loss-of-function causes hepatic steatosis (familial hypobetalipoproteinemia)", "drug_target_genetic_risk"),
    ("MBOAT7", "up",   "phosphatidylinositol remodeling", "hepatocyte", "rs641738 MASLD/fibrosis risk", "drug_target_genetic_risk"),
]
drug_df = pd.DataFrame(
    drug_targets,
    columns=["gene", "expected_direction", "expected_pathway",
             "cell_type_relevance", "evidence_summary", "source"],
)
drug_df = drug_df[
    ["gene", "source", "expected_direction", "expected_pathway",
     "cell_type_relevance", "evidence_summary"]
]
print(f"[c] Drug targets (curated): {len(drug_df)} genes", file=sys.stderr)

# ---------------------------------------------------------------------------
# (d) Cas13 positive controls (curated)
# ---------------------------------------------------------------------------
pos = pd.read_csv(ROOT / "results/library/positive_control.csv")
# Steatosis_Change_upon_KD: "Decrease" => KD reduces steatosis => gene is pro-disease
# (expected UP in disease). "Increase" => KD worsens steatosis => protective (expected DOWN).
def steatosis_to_dir(s):
    if not isinstance(s, str):
        return "either"
    s_low = s.strip().lower()
    if s_low.startswith("decrease"):
        return "up"
    if s_low.startswith("increase"):
        return "down"
    return "either"

pos_clean = pd.DataFrame({
    "gene": pos["Gene symbol"].astype(str).str.strip(),
    "source": "cas13_pos_control",
    "expected_direction": pos["Steatosis_Change_upon_KD"].apply(steatosis_to_dir),
    "expected_pathway": pos["Pathway / role"].fillna("steatosis_literature").astype(str).str[:80],
    "cell_type_relevance": "hepatocyte",
    "evidence_summary": pos.apply(
        lambda r: (
            f"Cas13 library positive control; KD effect on steatosis = "
            f"{r['Steatosis_Change_upon_KD']}"
            + (f"; clinical drug={r['Clinical_Drug']} ({r['Clinical_Phase']})"
               if isinstance(r["Clinical_Drug"], str) and r["Clinical_Drug"].strip()
               else "")
        ),
        axis=1,
    ),
})
pos_clean = pos_clean[pos_clean["gene"].str.len() > 0]
print(f"[d] Cas13 positive controls: {len(pos_clean)} genes", file=sys.stderr)

# ---------------------------------------------------------------------------
# Combine + dedupe by gene (concat sources with ';')
# ---------------------------------------------------------------------------
combined = pd.concat([conserved, coloc_hi, drug_df, pos_clean], ignore_index=True)
combined["gene"] = combined["gene"].astype(str).str.strip()
combined = combined[combined["gene"].str.len() > 0]

# Aggregate per gene
def agg_unique(series):
    seen, out = set(), []
    for v in series:
        v = str(v).strip()
        if v and v not in seen and v.lower() != "nan":
            seen.add(v)
            out.append(v)
    return ";".join(out)

def collapse_direction(series):
    vals = set(str(v) for v in series if pd.notna(v))
    if vals == {"up"}:
        return "up"
    if vals == {"down"}:
        return "down"
    if vals <= {"up", "either"}:
        return "up"
    if vals <= {"down", "either"}:
        return "down"
    if "up" in vals and "down" in vals:
        return "conflicting"
    return "either"

agg = combined.groupby("gene", as_index=False).agg({
    "source": agg_unique,
    "expected_direction": collapse_direction,
    "expected_pathway": agg_unique,
    "cell_type_relevance": agg_unique,
    "evidence_summary": agg_unique,
})
print(f"[combined] {len(agg)} unique genes", file=sys.stderr)

# ---------------------------------------------------------------------------
# Cross-reference with multi_evidence_atlas (atlas_missing flag)
# Read only the gene-symbol column to keep memory tiny.
# ---------------------------------------------------------------------------
atlas_path = ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
atlas_genes = pd.read_csv(atlas_path, usecols=["human_symbol"], low_memory=False)
atlas_set = set(atlas_genes["human_symbol"].dropna().astype(str).str.strip())
print(f"[atlas] {len(atlas_set)} unique symbols in multi_evidence_atlas", file=sys.stderr)

agg["atlas_missing"] = ~agg["gene"].isin(atlas_set)

# Sort: drug targets first, then COLOC, then conserved, then pos controls (gene order alpha within)
def source_rank(s: str) -> int:
    if "drug_target_approved" in s:
        return 0
    if "drug_target_trial" in s:
        return 1
    if "drug_target_genetic_risk" in s:
        return 2
    if "coloc_pp4_high" in s:
        return 3
    if "cas13_pos_control" in s:
        return 4
    return 5

agg["_rank"] = agg["source"].apply(source_rank)
agg = agg.sort_values(["_rank", "gene"]).drop(columns=["_rank"]).reset_index(drop=True)

out_csv = OUT_DIR / "backtest_set.csv"
agg.to_csv(out_csv, index=False)
print(f"[out] wrote {out_csv}", file=sys.stderr)

# ---------------------------------------------------------------------------
# Summary markdown
# ---------------------------------------------------------------------------
def count_with(substr):
    return int(agg["source"].str.contains(substr, regex=False).sum())

n_total = len(agg)
n_conserved = count_with("conserved_core")
n_coloc = count_with("coloc_pp4_high")
n_drug_appr = count_with("drug_target_approved")
n_drug_trial = count_with("drug_target_trial")
n_drug_gen = count_with("drug_target_genetic_risk")
n_pos = count_with("cas13_pos_control")

n_multi_source = int((agg["source"].str.contains(";")).sum())
n_atlas_missing = int(agg["atlas_missing"].sum())

dir_counts = agg["expected_direction"].value_counts().to_dict()

summary = f"""# Phase 1 Backtest Validation Set — Summary

Generated by `Analysis/Perturbation/scripts/build_backtest_set.py`.
Output: `Analysis/Perturbation/results/validation_sets/backtest_set.csv`

## Totals

| Metric | N |
|---|---|
| Unique genes (total) | **{n_total}** |
| Multi-source (≥2 evidence streams) | {n_multi_source} |
| Missing from multi_evidence_atlas | {n_atlas_missing} |
| Atlas coverage | {(n_total - n_atlas_missing)/n_total*100:.1f}% |

## Per-source counts (genes can appear in multiple sources)

| Source | N |
|---|---|
| `conserved_core` (cross-species concordant DEG) | {n_conserved} |
| `coloc_pp4_high` (SuSiE PP4 > 0.8) | {n_coloc} |
| `drug_target_approved` | {n_drug_appr} |
| `drug_target_trial` | {n_drug_trial} |
| `drug_target_genetic_risk` | {n_drug_gen} |
| `cas13_pos_control` (curated steatosis panel) | {n_pos} |

## Expected direction (disease vs healthy)

| Direction | N |
|---|---|
""" + "\n".join(f"| {k} | {v} |" for k, v in sorted(dir_counts.items(), key=lambda x: -x[1])) + f"""

## Schema (columns)

| Column | Description |
|---|---|
| `gene` | HGNC symbol (canonical) |
| `source` | One or more of: conserved_core, coloc_pp4_high, drug_target_approved, drug_target_trial, drug_target_genetic_risk, cas13_pos_control. Multiple sources concatenated with `;`. |
| `expected_direction` | up / down / either / conflicting (in disease vs control) |
| `expected_pathway` | Short biological annotation (`;`-joined across sources) |
| `cell_type_relevance` | hepatocyte / stellate / macrophage / multi / general |
| `evidence_summary` | One-line provenance per source (`;`-joined) |
| `atlas_missing` | TRUE if gene symbol is not found in `multi_evidence_atlas.csv` |

## Provenance

- (a) `Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv` — `primary_category == "Conserved"` (n=1,108 expected; observed {n_conserved} carried through dedup).
- (b) `GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv` — `coloc_best_susie_pp4 > 0.8`. PolyFun EUR LD reference (production as of 2026-05-06); BBJ EAS via 1KG EAS LD. Both ABF and SuSiE posteriors preserved upstream; this set uses SuSiE only.
- (c) Curated MASLD drug targets (n={len(drug_df)}): approved (resmetirom/THRB, semaglutide/GLP1R), Phase 2/3 trial (obeticholic acid/NR1H4, efruxifermin/FGF21+FGFR1, pegozafermin/FGF21, tirzepatide/GIPR, aramchol/SCD, firsocostat/ACACA-ACACB, ervogastat/DGAT2, cenicriviroc/CCR2-CCR5, anti-TGFβ Abs/TGFB1, FXR-SHP), and established genetic-risk loci (PNPLA3, HSD17B13, MTARC1, TM6SF2, GCKR, APOB, MBOAT7).
- (d) `results/library/positive_control.csv` — Cas13 library steatosis literature panel (curated by James Lee). `expected_direction` derived from `Steatosis_Change_upon_KD` (Decrease → gene is pro-disease, expected UP in disease; Increase → protective, expected DOWN).

## Use as canonical backtest set

Phase 1 zero-shot model evaluation: any perturbation-prediction model whose backtest recovery (recall at top-N) is at chance on this gene set is ejected from the consensus pool for that arm. The set is intentionally heterogeneous (transcriptomic + genetic + pharmacological + literature) so models cannot game any single axis.
"""

(OUT_DIR / "backtest_set_summary.md").write_text(summary)
print(f"[out] wrote {OUT_DIR/'backtest_set_summary.md'}", file=sys.stderr)
print(f"[done] total={n_total}, atlas_missing={n_atlas_missing}", file=sys.stderr)

#!/usr/bin/env python
"""
461b_compare_clean_vs_orig.py — Side-by-side comparison of original vs CLEAN
(GSE136103 + dubious-Healthy excluded) changepoint + phenotype-screen outputs.

Writes:
  results_gpu_v2/mcp/protocol_sensitivity_comparison.tsv
  results_gpu_v2/mcp/protocol_sensitivity_summary.md
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
INT = MCP / "integration"

# --- Changepoint ---
cp_orig  = pd.read_csv(INT / "program_changepoint.tsv", sep="\t")
cp_clean = pd.read_csv(INT / "program_changepoint_clean.tsv", sep="\t")

cp_orig  = cp_orig.rename(columns={c: f"{c}_orig"  for c in cp_orig.columns  if c != "program"})
cp_clean = cp_clean.rename(columns={c: f"{c}_clean" for c in cp_clean.columns if c != "program"})
cp = cp_orig.merge(cp_clean, on="program", how="outer")

def stage_bin(bp):
    if pd.isna(bp):
        return "none"
    if bp <= 0.5:
        return "<=0.5_pre"        # below Healthy/Steatosis transition
    if bp <= 1.0:
        return "0.5-1.0_H-St"     # Healthy→Steatosis quarter
    if bp <= 1.5:
        return "1.0-1.5_St-edge"  # approaching Steatosis-SH transition
    if bp <= 2.0:
        return "1.5-2.0_St-SH"    # crossed the F1->F2 / Steatosis→SH transition
    if bp <= 2.5:
        return "2.0-2.5_SH-Cir"   # past SH towards Cirrhosis
    return ">2.5_post"

cp["bp_bin_orig"]  = cp["breakpoint_stage_orig"].apply(stage_bin)
cp["bp_bin_clean"] = cp["breakpoint_stage_clean"].apply(stage_bin)
cp["bp_moved"] = (cp["best_model_orig"] != cp["best_model_clean"]) | \
                 (cp["bp_bin_orig"] != cp["bp_bin_clean"])

n_moved_bin = int(cp["bp_moved"].sum())
n_model_changed = int((cp["best_model_orig"] != cp["best_model_clean"]).sum())

# Distance moved among those that were segmented in BOTH
both_seg = cp[(cp["best_model_orig"] == "seg1") & (cp["best_model_clean"] == "seg1")]
moves_within_seg = (both_seg["breakpoint_stage_orig"] - both_seg["breakpoint_stage_clean"]).abs()

# --- Phenotype screen ---
ph_orig  = pd.read_csv(INT / "phenotype_screen_long.tsv", sep="\t")
ph_clean = pd.read_csv(INT / "phenotype_screen_long_clean.tsv", sep="\t")

m_ph = ph_orig.merge(
    ph_clean, on=["program", "phenotype"], how="outer",
    suffixes=("_orig", "_clean"),
)
m_ph["sig_orig"]  = (m_ph["q_orig"]  < 0.05).fillna(False)
m_ph["sig_clean"] = (m_ph["q_clean"] < 0.05).fillna(False)
m_ph["flip"] = m_ph["sig_orig"] != m_ph["sig_clean"]

n_flipped = int(m_ph["flip"].sum())
n_lost    = int((m_ph["sig_orig"] & ~m_ph["sig_clean"]).sum())
n_gained  = int((~m_ph["sig_orig"] & m_ph["sig_clean"]).sum())

# --- Write combined comparison table ---
cp_cols = ["program",
           "best_model_orig", "best_model_clean",
           "breakpoint_stage_orig", "breakpoint_stage_clean",
           "bp_bin_orig", "bp_bin_clean", "bp_moved",
           "n_orig", "n_clean",
           "delta_aic_orig", "delta_aic_clean"]
cp_cols = [c for c in cp_cols if c in cp.columns]
cp[cp_cols].to_csv(MCP / "protocol_sensitivity_comparison_changepoint.tsv", sep="\t", index=False)

ph_cols = ["program", "phenotype",
           "p_orig", "p_clean", "q_orig", "q_clean",
           "beta_orig", "beta_clean",
           "sig_orig", "sig_clean", "flip",
           "n_orig", "n_clean"]
ph_cols = [c for c in ph_cols if c in m_ph.columns]
m_ph[ph_cols].sort_values(
    by=["flip", "q_clean"], ascending=[False, True]
).to_csv(MCP / "protocol_sensitivity_comparison_phenotype.tsv", sep="\t", index=False)

# --- Console + summary numbers ---
print("=" * 70)
print("Protocol sensitivity comparison")
print("=" * 70)
print(f"Programs total: {len(cp)}")
print(f"  best_model changed (orig vs clean): {n_model_changed}")
print(f"  breakpoint-bin moved (incl. model change): {n_moved_bin}")
print(f"  in both seg1: {len(both_seg)}")
if len(both_seg) > 0:
    print(f"  abs(orig - clean) stage move among shared seg1 fits — "
          f"mean={moves_within_seg.mean():.3f}, max={moves_within_seg.max():.3f}")

print(f"\nPhenotype rows: {len(m_ph)}")
print(f"  significance flipped (q<0.05 boundary): {n_flipped}")
print(f"    lost (sig→ns)   : {n_lost}")
print(f"    gained (ns→sig) : {n_gained}")

# Show flips
flips = m_ph.loc[m_ph["flip"], ph_cols].copy()
if len(flips):
    print("\n--- Significance flips ---")
    print(flips.to_string(index=False, max_colwidth=22))

# Show movers
movers = cp.loc[cp["bp_moved"], cp_cols].copy()
if len(movers):
    print("\n--- Changepoint movers ---")
    print(movers.to_string(index=False, max_colwidth=22))

# --- Verdict ---
# Heuristic: major if any bp crosses the 1.5 (Steatosis→SH) boundary by >0.25
#   or any phenotype flip involves disease_stage_numeric / condition_binary
key_phenos = {"disease_stage_numeric", "condition_binary_num"}
key_flips = m_ph.loc[m_ph["flip"] & m_ph["phenotype"].isin(key_phenos)]
major_bp_move = False
for _, row in both_seg.iterrows():
    if pd.isna(row["breakpoint_stage_orig"]) or pd.isna(row["breakpoint_stage_clean"]):
        continue
    # crossed the 1.5 boundary
    if (row["breakpoint_stage_orig"] - 1.5) * (row["breakpoint_stage_clean"] - 1.5) < 0:
        major_bp_move = True

if major_bp_move or len(key_flips) > 0:
    verdict = "MAJOR"
elif n_moved_bin > 0 or n_flipped > 0:
    verdict = "MODERATE"
else:
    verdict = "MINOR"

print(f"\nVERDICT: {verdict}")

# --- Markdown summary ---
md = []
md.append("# cNMF Phase G protocol-contamination sensitivity\n")
md.append("**Question.** Did `460_phenotype_screen.py` / `461_changepoint.R` change conclusions when GSE136103 NPC-enriched donors and dubious-Healthy donors were removed?\n")
md.append("**Exclusion rule.** `exclude_stage_analysis == TRUE` in `donor_metadata_extended.tsv` (58 of 269 donors; 49 dubious-Healthy + 9 GSE136103 Cirrhosis).\n")
md.append(f"**Donors used.** Original: n={int(cp['n_orig'].max()) if 'n_orig' in cp else 'NA'} (per-program). Clean: n={int(cp['n_clean'].max()) if 'n_clean' in cp else 'NA'}.\n")
md.append("\n## Results\n")
md.append(f"- **Changepoint best-model changed**: {n_model_changed} / {len(cp)} programs.\n")
md.append(f"- **Breakpoint bin moved (incl. model change)**: {n_moved_bin} / {len(cp)} programs.\n")
if len(both_seg) > 0:
    md.append(f"- Among the {len(both_seg)} programs segmented in both, the mean |Δ breakpoint stage| = {moves_within_seg.mean():.3f} (max {moves_within_seg.max():.3f}).\n")
md.append(f"- **Phenotype significance flipped (q<0.05 boundary)**: {n_flipped} of {len(m_ph)} rows. Lost {n_lost} (sig→ns). Gained {n_gained} (ns→sig).\n")
md.append(f"- Disease-axis-relevant flips (disease_stage_numeric / condition_binary_num): {len(key_flips)}.\n")
md.append(f"- Any breakpoint crossed the Steatosis→SH (1.5) boundary direction: {major_bp_move}.\n")
md.append("\n## Verdict\n")
if verdict == "MAJOR":
    md.append("**MAJOR** — manuscript text and figures referencing cNMF program changepoints / phenotype enrichments need updating.\n")
elif verdict == "MODERATE":
    md.append("**MODERATE** — a footnote/methods sentence should disclose the protocol-contamination sensitivity. Headline statements unchanged.\n")
else:
    md.append("**MINOR** — protocol contamination does not alter cNMF Phase G conclusions; manuscript text unchanged.\n")
md.append("\n## Caveat\n")
md.append("After exclusion, the disease axis collapses to stages 0–2 (Healthy / Steatosis / Steatohepatitis); all donors with `disease_stage_numeric == 3` (Cirrhosis) come from GSE136103 and are dropped. Cirrhosis donors from GSE202379 carry NA `disease_stage_numeric` so they don't restore the F3 bin. The clean segmented breakpoints therefore live on a [0,2] range and are not directly numerically comparable to original breakpoints (which spanned [0,3]). The comparison interpretation focuses on (a) whether the inflection sits below or above the Steatosis→SH (1.5) midpoint, and (b) whether phenotype-screen significance is preserved.\n")

(MCP / "protocol_sensitivity_summary.md").write_text("".join(md))
print(f"\nWrote: {MCP / 'protocol_sensitivity_summary.md'}")
print(f"Wrote: {MCP / 'protocol_sensitivity_comparison_changepoint.tsv'}")
print(f"Wrote: {MCP / 'protocol_sensitivity_comparison_phenotype.tsv'}")

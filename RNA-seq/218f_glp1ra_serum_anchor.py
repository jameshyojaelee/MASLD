#!/usr/bin/env python
"""
Analysis (human serum clinical anchor) -- presented SEPARATELY, not fused with the
mouse liver target partition (gene-level cross-compartment concordance is chance; RevB/RevC).

Summarizes the human serum-proteomic reversion (Jara Nat Med 2025): semaglutide reverts
the MASH serum proteome toward healthy; decomposes weight-dependence; reports overlap of
the 72-protein signature with our convergence atlas (small, by construction -- serum
proteins are secreted, most liver targets are not).

Output: RNA-seq/results/glp1ra/glp1ra_serum_anchor_summary.txt
"""
import os
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
G = os.path.join(ROOT, "RNA-seq/results/glp1ra")
SIG = os.path.join(G, "human_semaglutide_signature/signature_genes.csv")
FULL = os.path.join(G, "human_semaglutide_signature/somascan_full_4979_proteins.csv")
REF = os.path.join(G, "disease_reference.csv")

sig = pd.read_csv(SIG)
full = pd.read_csv(FULL).rename(columns={"gene": "human_symbol"})
ref = pd.read_csv(REF)
conv_t1 = set(ref.loc[ref["conv_tier1"] == True, "human_symbol"])

# reversion of the 72-protein signature
sig["reverts"] = sig["reversal"].astype(str).str.lower().eq("reversal")
n = len(sig)
n_down = int((sig["direction"].astype(str).str.contains("down")).sum())
n_up = n - n_down
n_revert = int(sig["reverts"].sum())
n_wind = int(sig["weight_loss_independent"].astype(str).str.upper().isin(["Y", "TRUE", "1"]).sum())
n_indep_sig = int(sig["indep_signif"].astype(str).str.upper().isin(["Y", "TRUE", "1"]).sum())

# overlap of signature with our convergence Tier-1 atlas
sig_in_conv = sorted(set(sig["gene_symbol"]) & conv_t1)

# whole-proteome reversion rate (disease-sig serum proteins whose sema effect opposes disease)
full["s_dis"] = pd.to_numeric(full["indep_MASHvsHealthy_logFC"], errors="coerce")
full["s_trt"] = pd.to_numeric(full["sema_effect_est1"], errors="coerce")
full["s_dis_sig"] = full["indep_signif"].astype(str).str.upper().isin(["Y", "TRUE", "1"])
ds = full[full["s_dis_sig"] & full["s_dis"].notna() & full["s_trt"].notna()]
revert_rate = (np.sign(ds["s_trt"]) != np.sign(ds["s_dis"])).mean()

L = ["Human serum clinical anchor (Jara Nat Med 2025) -- SEPARATE line", "=" * 60,
     f"72-protein semaglutide treatment signature: {n_down} down / {n_up} up with semaglutide",
     f"  revert the MASH-vs-healthy serum direction: {n_revert}/{n}",
     f"  significant in independent MASH-vs-healthy cohort (CoCoMASLD): {n_indep_sig}/{n}",
     f"  weight-loss-INDEPENDENT: {n_wind}/{n}",
     f"whole disease-serum-proteome (n={len(ds)} disease-sig proteins): "
     f"{revert_rate*100:.0f}% reverted in direction by semaglutide",
     "",
     "Weight-dependence (mediation, from paper):",
     "  MASH resolution: 69.3% weight-loss-mediated | Steatosis 82.8% | Ballooning 71.6% | "
     "Fibrosis 25.1% (largely weight-INDEPENDENT)",
     "",
     f"Overlap of the 72-protein signature with our convergence Tier-1 atlas: "
     f"{len(sig_in_conv)} genes -> {sig_in_conv}",
     "  (small by construction: serum signature = secreted proteins; most liver convergence",
     "   targets are intracellular. This is WHY serum is a separate clinical anchor, not fused.)",
     "",
     "Hallmark GSEA (paper): 14 sets, ALL down with semaglutide, all METABOLIC",
     "  (fatty-acid metab, OXPHOS, xenobiotic, adipogenesis, mTORC1, bile-acid, glycolysis...).",
     "Named: SERPINC1/APOF up-in-MASH raised by sema; ADAMTSL2/ACY1/AKR1B10/TREM2 down."]
summary = "\n".join(L)
with open(os.path.join(G, "glp1ra_serum_anchor_summary.txt"), "w") as fh:
    fh.write(summary + "\n")
print(summary)

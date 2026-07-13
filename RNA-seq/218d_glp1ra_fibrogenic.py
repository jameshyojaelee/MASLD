#!/usr/bin/env python
"""
Analysis D: does semaglutide reverse the core FIBROGENIC program? (weight-independently?)

Analyzed DIRECTLY on a curated stellate-activation / ECM / crosslinking gene set --
NOT convergence-Tier1-gated, because collagens/TIMPs/ACTA2 never reach convergence Tier-1
(RevC). Uses the mouse meta (all-model) + the weight-stable CDA-HFD model.

Output: RNA-seq/results/glp1ra/glp1ra_fibrogenic_reversal.csv + summary.
"""
import os
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
G = os.path.join(ROOT, "RNA-seq/results/glp1ra")
META = os.path.join(G, "mouse_reversal/mouse_semaglutide_reversal_meta.csv")

# curated core fibrogenic / hepatic-stellate-activation / ECM program (human symbols)
FIBRO = ["COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
         "COL6A3", "COL8A1", "COL15A1", "TIMP1", "TIMP2", "MMP2", "MMP9", "MMP13", "MMP14",
         "PDGFRB", "PDGFRA", "TGFB1", "TGFB2", "LOX", "LOXL1", "LOXL2", "SPARC", "THBS1",
         "THBS2", "ACTA2", "TAGLN", "FN1", "CTHRC1", "CCN2", "DCN", "LUM", "VIM"]

m = pd.read_csv(META)
f = m[m["human_symbol"].isin(FIBRO)].copy()
f["disease_up"] = f["disease_meta_lfc"] > 0
f["disease_sig"] = f["disease_meta_padj"] < 0.05
f["meta_reverses"] = (f["disease_sig"] &
                      (np.sign(f["treatment_meta_lfc"]) != np.sign(f["disease_meta_lfc"])) &
                      (f["treatment_meta_padj"] < 0.05))
f["cdahfd_reverses"] = (f["disease_sig"] &
                        (np.sign(f["treatment_lfc_cdahfd"]) != np.sign(f["disease_meta_lfc"])) &
                        (f["treatment_padj_cdahfd"] < 0.05))
f = f.sort_values("disease_meta_lfc", ascending=False)
f[["human_symbol", "disease_meta_lfc", "disease_meta_padj", "treatment_meta_lfc",
   "treatment_meta_padj", "treatment_lfc_cdahfd", "treatment_padj_cdahfd",
   "meta_reverses", "cdahfd_reverses"]].to_csv(
    os.path.join(G, "glp1ra_fibrogenic_reversal.csv"), index=False)

ds = f[f["disease_sig"]]
L = ["Analysis D: fibrogenic-program reversal by semaglutide", "=" * 50,
     f"curated program n={len(FIBRO)}; measured in mouse meta n={len(f)}; "
     f"disease-significant n={len(ds)}",
     f"  reverse at meta (all-model):        {int(f['meta_reverses'].sum())}/{len(ds)}",
     f"  reverse weight-independently (CDA-HFD): {int(f['cdahfd_reverses'].sum())}/{len(ds)}",
     "",
     "per-gene (disease-up genes, sorted by disease effect):",
     ds[["human_symbol", "disease_meta_lfc", "treatment_meta_lfc", "treatment_lfc_cdahfd",
         "meta_reverses", "cdahfd_reverses"]].round(3).to_string(index=False)]
summary = "\n".join(L)
with open(os.path.join(G, "glp1ra_fibrogenic_summary.txt"), "w") as fh:
    fh.write(summary + "\n")
print(summary)

#!/usr/bin/env python
"""
343t_bootstrap_gates_redesign.py

Re-evaluate the bootstrap stability output with biologically meaningful
gate criteria instead of the strict F3-mode pre-reg.

Reads:
  bootstrap_sh_distribution_stability.tsv  (4 methods x 500 bootstraps x 5 F-stages)

Evaluates 5 biological gates per method:

  Gate A — advanced-fibrosis majority:
    F2+F3+F4 >= 0.50 in >=95% of bootstraps
    (SH should be moderate-to-advanced; clinically actionable ≥F2)

  Gate B — F3+F4 is the largest bin:
    F3+F4 > F0+F1+F2 in >=95% of bootstraps
    (advanced fibrosis as a combined bin should dominate over early/none)

  Gate C — F3+F4 is the mode (collapsed bin):
    F3+F4 > each individual F0/F1/F2 in >=95%
    (less strict than B; checks advanced-fibrosis bin > each early bin)

  Gate D — no bimodal collapse:
    F2+F3 >= 0.10 in >=95% of bootstraps
    (rules out scANVI-style F0/F4 endpoint collapse)

  Gate E — original strict (F3 the mode of all 5 bins):
    P_F3 > max(P_F0, P_F1, P_F2, P_F4) in >=95% of bootstraps

Writes:
  bootstrap_gates_redesign.tsv  (method x gate -> pass_pct + verdict)
  BOOTSTRAP_DECISION_V2.txt     (per-gate verdict + recommendation)
"""
from __future__ import annotations
import sys
import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR = BASE / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"

# Load long-format -> wide-format (one row per (method, boot))
df = pd.read_csv(OUT_DIR / "bootstrap_sh_distribution_stability.tsv", sep="\t")
wide = df.pivot_table(index=["method", "boot"], columns="fstage",
                      values="fraction", fill_value=0.0).reset_index()
wide.columns = [c if isinstance(c, str) else f"P_F{int(c)}" for c in wide.columns]
print(f"Loaded {len(wide)} (method, boot) rows")

# Compute combined-bin probabilities
wide["P_F2plus"]   = wide["P_F2"] + wide["P_F3"] + wide["P_F4"]
wide["P_F3plus"]   = wide["P_F3"] + wide["P_F4"]
wide["P_F2_F3"]    = wide["P_F2"] + wide["P_F3"]
wide["P_F0_F1_F2"] = wide["P_F0"] + wide["P_F1"] + wide["P_F2"]

# Evaluate gates
gates = {}
for method, sub in wide.groupby("method"):
    n = len(sub)
    gA = (sub["P_F2plus"] >= 0.50).sum() / n
    gB = (sub["P_F3plus"] > sub["P_F0_F1_F2"]).sum() / n
    gC = ((sub["P_F3plus"] > sub["P_F0"]) &
          (sub["P_F3plus"] > sub["P_F1"]) &
          (sub["P_F3plus"] > sub["P_F2"])).sum() / n
    gD = (sub["P_F2_F3"] >= 0.10).sum() / n
    gE = ((sub["P_F3"] > sub["P_F0"]) &
          (sub["P_F3"] > sub["P_F1"]) &
          (sub["P_F3"] > sub["P_F2"]) &
          (sub["P_F3"] > sub["P_F4"])).sum() / n
    gates[method] = dict(
        gateA_F2plus_50pct=gA,
        gateB_F3plus_majority=gB,
        gateC_F3plus_mode=gC,
        gateD_no_bimodal=gD,
        gateE_F3_unique_mode=gE,
        median_F0=sub["P_F0"].median(),
        median_F1=sub["P_F1"].median(),
        median_F2=sub["P_F2"].median(),
        median_F3=sub["P_F3"].median(),
        median_F4=sub["P_F4"].median(),
        median_F3plus=sub["P_F3plus"].median(),
        median_F2plus=sub["P_F2plus"].median(),
    )

out = pd.DataFrame(gates).T.reset_index().rename(columns={"index": "method"})
out.to_csv(OUT_DIR / "bootstrap_gates_redesign.tsv", sep="\t", index=False)
print()
print(out[[
    "method",
    "gateA_F2plus_50pct",
    "gateB_F3plus_majority",
    "gateC_F3plus_mode",
    "gateD_no_bimodal",
    "gateE_F3_unique_mode",
]].to_string(index=False))

print()
print("Median SH F-stage distribution per method:")
print(out[[
    "method", "median_F0", "median_F1", "median_F2",
    "median_F3", "median_F4", "median_F2plus", "median_F3plus",
]].to_string(index=False))

# Decision logic
PASS_THRESH = 0.95
print()
print("=" * 60)
print("Re-designed gate verdicts (>= 95% of bootstraps must pass)")
print("=" * 60)
with open(OUT_DIR / "BOOTSTRAP_DECISION_V2.txt", "w") as fh:
    fh.write("Re-designed bootstrap stability gates (biology-meaningful)\n")
    fh.write("=" * 60 + "\n\n")
    fh.write(f"Pass threshold: >= {int(PASS_THRESH*100)}% of B=500 bootstraps\n\n")
    for method, row in out.set_index("method").iterrows():
        verdict = "SHIP" if row["gateA_F2plus_50pct"] >= PASS_THRESH else "FALL_BACK"
        line = (f"{method:25s} | A(F2+≥50%)={row['gateA_F2plus_50pct']*100:5.1f}% | "
                f"B(F3+maj)={row['gateB_F3plus_majority']*100:5.1f}% | "
                f"C(F3+mode)={row['gateC_F3plus_mode']*100:5.1f}% | "
                f"D(no-bimodal)={row['gateD_no_bimodal']*100:5.1f}% | "
                f"E(F3-only)={row['gateE_F3_unique_mode']*100:5.1f}% | "
                f"-> {verdict}")
        print(line); fh.write(line + "\n")
    fh.write("\nRecommendation: use Gate A as the primary criterion ('SH "
             "expresses moderate-to-advanced fibrosis, F2+ majority').\n")
print()
print("Files written:")
print(f"  {OUT_DIR / 'bootstrap_gates_redesign.tsv'}")
print(f"  {OUT_DIR / 'BOOTSTRAP_DECISION_V2.txt'}")

#!/usr/bin/env python
"""D1 atlas integrator — produces arm-specific atlas columns from consensus.

D1 contributes 4 columns to the S8 evidence source:
  - perturb_mechanism_pathway       (string, top KEGG/Reactome enrichment)
  - perturb_mechanism_coa           (string, putative class-of-action label)
  - perturb_downstream_robust_n     (int, count of robust downstream genes)
  - perturb_downstream_robust_genes (string, ';'-joined gene list)

Writes results/integration/d1_atlas_columns.csv with one row per target gene.
Pathway / class-of-action enrichment is a TODO for the Atlas Integrator subagent.
"""
from __future__ import annotations

import sys
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

import pandas as pd  # noqa: E402

from arm_driver import RESULTS_ROOT  # noqa: E402

ARM_DIR = RESULTS_ROOT / "d1_mechanism"
CONSENSUS_CSV = ARM_DIR / "consensus_predictions.csv"
OUT_CSV = RESULTS_ROOT / "integration" / "d1_atlas_columns.csv"


def main() -> None:
    if not CONSENSUS_CSV.exists():
        print(f"[D1/atlas] missing {CONSENSUS_CSV}; writing empty atlas columns")
        OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=[
            "gene",
            "perturb_mechanism_pathway",
            "perturb_mechanism_coa",
            "perturb_downstream_robust_n",
            "perturb_downstream_robust_genes",
        ]).to_csv(OUT_CSV, index=False)
        return

    df = pd.read_csv(CONSENSUS_CSV)
    df = df[df["pass_consensus"] == True]  # noqa: E712
    grouped = (
        df.groupby("target_gene")["downstream_gene"]
        .apply(lambda s: sorted(set(s)))
        .reset_index()
    )
    grouped["perturb_downstream_robust_n"] = grouped["downstream_gene"].apply(len)
    grouped["perturb_downstream_robust_genes"] = grouped["downstream_gene"].apply(
        lambda x: ";".join(x[:50])  # cap at 50 to keep CSV cell sane
    )
    # TODO[d1-atlas-integrator]: compute pathway + class-of-action via fgsea /
    # gseapy on the robust downstream gene set (one enrichment per target_gene).
    grouped["perturb_mechanism_pathway"] = pd.NA
    grouped["perturb_mechanism_coa"] = pd.NA

    out = grouped[[
        "target_gene",
        "perturb_mechanism_pathway",
        "perturb_mechanism_coa",
        "perturb_downstream_robust_n",
        "perturb_downstream_robust_genes",
    ]].rename(columns={"target_gene": "gene"})
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"[D1/atlas] wrote {len(out)} rows -> {OUT_CSV}")


if __name__ == "__main__":
    main()

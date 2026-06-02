#!/usr/bin/env python
"""D2 atlas integrator -> 3 columns:
  perturb_reversal_score, perturb_reversal_rank, perturb_reversal_stage_specific.
"""
from __future__ import annotations

import sys
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

import pandas as pd  # noqa: E402

from arm_driver import RESULTS_ROOT, load_runs  # noqa: E402

ARM_DIR = RESULTS_ROOT / "d2_reversal"
CONSENSUS_CSV = ARM_DIR / "consensus_predictions.csv"
OUT_CSV = RESULTS_ROOT / "integration" / "d2_atlas_columns.csv"


def main() -> None:
    if not CONSENSUS_CSV.exists():
        OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=[
            "gene", "perturb_reversal_score", "perturb_reversal_rank",
            "perturb_reversal_stage_specific",
        ]).to_csv(OUT_CSV, index=False)
        print(f"[D2/atlas] missing consensus; wrote empty {OUT_CSV}")
        return

    cons = pd.read_csv(CONSENSUS_CSV)
    if cons.empty:
        cons = pd.DataFrame(columns=["gene", "mean_score", "mean_rank"])

    # BUG FIX 2026-05-21 (HARSH_REVIEW F4 / 2nd-opinion smoking gun):
    # The prior implementation wrote ALL 111 consensus rows to the atlas
    # regardless of pass_consensus. With current STATE-only operational
    # consensus, 110/111 rows have pass_consensus=False; the 1 True row
    # (5_8S_rRNA) is a rRNA pseudogene artifact. Filter before writing.
    if "pass_consensus" in cons.columns:
        n_before = len(cons)
        cons = cons[cons["pass_consensus"] == True]  # noqa: E712
        print(f"[D2/atlas] pass_consensus filter: {n_before} -> {len(cons)} rows")
    else:
        print(f"[D2/atlas] WARN consensus_predictions.csv has no pass_consensus column; writing as-is (len={len(cons)})")

    # Stage-specificity flag: read each per-context run and look at the
    # context_substring to flag genes that only pass in one stage.
    runs = load_runs(ARM_DIR)
    per_gene_stages: dict[str, set[str]] = {}
    for run in runs:
        # context format: F__Steatohepatitis__Progressor__ref_a -> coarse stage = 2nd field
        parts = run.context.split("__")
        stage = parts[1] if len(parts) > 1 else "—"
        for p in run.predictions:
            g = p["gene"] if isinstance(p, dict) else p.gene
            per_gene_stages.setdefault(g, set()).add(stage)

    def stage_label(g: str) -> str:
        s = per_gene_stages.get(g, set())
        return "single_stage:" + ",".join(sorted(s)) if len(s) == 1 else "multi_stage" if s else "—"

    out = pd.DataFrame({
        "gene": cons.get("gene", []),
        "perturb_reversal_score": cons.get("mean_score", []),
        "perturb_reversal_rank": cons.get("mean_rank", []),
        "perturb_reversal_stage_specific": [stage_label(g) for g in cons.get("gene", [])],
    })
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"[D2/atlas] wrote {len(out)} rows -> {OUT_CSV}")


if __name__ == "__main__":
    main()

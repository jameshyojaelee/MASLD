#!/usr/bin/env python
"""D3 atlas integrator -> 3 columns:
  perturb_synergy_partners, perturb_synergy_top_n, perturb_synergy_class.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

import pandas as pd  # noqa: E402

from arm_driver import RESULTS_ROOT  # noqa: E402

ARM_DIR = RESULTS_ROOT / "d3_synergy"
CONSENSUS_CSV = ARM_DIR / "consensus_predictions.csv"
OUT_CSV = RESULTS_ROOT / "integration" / "d3_atlas_columns.csv"


def main() -> None:
    if not CONSENSUS_CSV.exists():
        OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=[
            "gene", "perturb_synergy_partners", "perturb_synergy_top_n", "perturb_synergy_class",
        ]).to_csv(OUT_CSV, index=False)
        print(f"[D3/atlas] missing consensus; wrote empty {OUT_CSV}")
        return

    cons = pd.read_csv(CONSENSUS_CSV)
    cons = cons[cons.get("pass_consensus", False).astype(bool)] if "pass_consensus" in cons else cons

    partners: dict[str, list[str]] = defaultdict(list)
    classes: dict[str, list[str]] = defaultdict(list)
    for _, row in cons.iterrows():
        genes = str(row["genes"]).split("+")
        cls = row.get("synergy_class", "additive")
        for g in genes:
            others = [x for x in genes if x != g]
            partners[g].extend(others)
            classes[g].append(cls)

    rows = []
    for g, plist in partners.items():
        plist = sorted(set(plist))
        cls_counts = pd.Series(classes[g]).value_counts()
        rows.append(dict(
            gene=g,
            perturb_synergy_partners=";".join(plist[:50]),
            perturb_synergy_top_n=len(plist),
            perturb_synergy_class=cls_counts.idxmax() if len(cls_counts) else "—",
        ))

    out = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"[D3/atlas] wrote {len(out)} rows -> {OUT_CSV}")


if __name__ == "__main__":
    main()

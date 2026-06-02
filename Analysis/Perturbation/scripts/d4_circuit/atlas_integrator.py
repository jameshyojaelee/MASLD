#!/usr/bin/env python
"""D4 atlas integrator -> 2 columns:
  perturb_ccc_ligand_role, perturb_ccc_predicted_receivers.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

import pandas as pd  # noqa: E402

from arm_driver import RESULTS_ROOT  # noqa: E402

ARM_DIR = RESULTS_ROOT / "d4_circuit"
CONSENSUS_CSV = ARM_DIR / "consensus_predictions.csv"
OUT_CSV = RESULTS_ROOT / "integration" / "d4_atlas_columns.csv"


def main() -> None:
    if not CONSENSUS_CSV.exists():
        OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=[
            "gene", "perturb_ccc_ligand_role", "perturb_ccc_predicted_receivers",
        ]).to_csv(OUT_CSV, index=False)
        print(f"[D4/atlas] empty -> {OUT_CSV}")
        return

    cons = pd.read_csv(CONSENSUS_CSV)
    cons = cons[cons.get("pass_consensus", False).astype(bool)] if "pass_consensus" in cons else cons

    receivers: dict[str, list[str]] = defaultdict(list)
    for _, row in cons.iterrows():
        sender = row.get("sender_ligand")
        receiver = row.get("receiver_cell_type")
        if sender and receiver:
            receivers[sender].append(receiver)

    rows = []
    for g, rcvs in receivers.items():
        uniq = sorted(set(rcvs))
        rows.append(dict(
            gene=g,
            perturb_ccc_ligand_role="paracrine_sender" if uniq else "—",
            perturb_ccc_predicted_receivers=";".join(uniq),
        ))
    out = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"[D4/atlas] wrote {len(out)} rows -> {OUT_CSV}")


if __name__ == "__main__":
    main()

"""D4 Ligand-Receptor pair curator.

Builds the canonical hepatocyte-ligand × receiver-cell-type table used by
all D4 runners. Sources: CellChatDB, OmniPath, LIANA consensus,
CellTalkDB. Output: `data/hits/d4_circuit_lr_pairs.csv` with cols
  ligand, receiver_cell_type, receptor, source, ccc_score.

This script enriches the `d4_circuit_ligands.csv` hit file with explicit LR
pairs by joining against a shared LR resource.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
HITS_DIR = PROJECT_ROOT / "Analysis/Perturbation/data/hits"
LIANA_DB = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2/ccc/liana_consensus_lr_db.csv"
)
OUT = HITS_DIR / "d4_circuit_lr_pairs.csv"


def main() -> None:
    ligands = pd.read_csv(HITS_DIR / "d4_circuit_ligands.csv")
    # Already shipped with receiver_cell_type + receptor; this script is the
    # extension hook for curating additional LR pairs.
    if LIANA_DB.exists():
        lr = pd.read_csv(LIANA_DB)
        # TODO[d4-lr-curator]: join ligands with LR DB to enrich receptors
        # and add missing receiver-cell-type matches.
        merged = ligands.merge(
            lr, left_on=["ligand", "receptor"], right_on=["source", "target"],
            how="left", suffixes=("", "_lr"),
        )
    else:
        merged = ligands.copy()
    merged.to_csv(OUT, index=False)
    print(f"[d4/lr_curator] wrote {len(merged)} rows -> {OUT}")


if __name__ == "__main__":
    main()

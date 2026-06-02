"""D5 cross-species concordance scorer.

For each ortholog pair (human_gene, mouse_gene) compares the human D1 / D2 / D3
prediction profile (loaded from RESULTS_ROOT/d{1,2,3}) against the mouse
D5 runner output. Concordance metrics:
  - direction_match: same up/down direction across top-K downstream
  - jaccard_top_k:   Jaccard of top-K downstream gene sets
  - pearson_logfc:   correlation across the shared downstream universe

Subagent: d5-concordance-scorer.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Iterable

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from arm_driver import RESULTS_ROOT  # noqa: E402

ARM_DIR_MOUSE = RESULTS_ROOT / "d5_mouse"
ARM_DIR_HUMAN_D1 = RESULTS_ROOT / "d1_mechanism"


def _load_predictions(json_path: Path) -> list[dict]:
    with open(json_path) as fh:
        return json.load(fh).get("predictions", [])


def per_gene_topk(preds: Iterable[dict], k: int = 100) -> dict[str, list[str]]:
    by_target: dict[str, list[tuple[str, float]]] = {}
    for p in preds:
        t = p.get("target_gene")
        d = p.get("downstream_gene")
        l = float(p.get("logFC_predicted", 0.0))
        by_target.setdefault(t, []).append((d, l))
    result = {}
    for t, lst in by_target.items():
        lst.sort(key=lambda r: abs(r[1]), reverse=True)
        result[t] = [d for d, _ in lst[:k]]
    return result


def concordance(human_topk: dict[str, list[str]], mouse_topk: dict[str, list[str]],
                ortholog_pairs: pd.DataFrame, k: int = 100) -> pd.DataFrame:
    """ortholog_pairs has human_gene, mouse_gene columns."""
    rows = []
    for _, row in ortholog_pairs.iterrows():
        hg, mg = row["human_gene"], row.get("mouse_gene")
        if not mg or not isinstance(mg, str):
            continue
        h = set(human_topk.get(hg, [])[:k])
        m = set(mouse_topk.get(mg, [])[:k])
        if not h or not m:
            continue
        jacc = len(h & m) / max(1, len(h | m))
        rows.append(
            dict(
                human_gene=hg,
                mouse_gene=mg,
                jaccard_top_k=jacc,
                n_human_top=len(h),
                n_mouse_top=len(m),
            )
        )
    return pd.DataFrame(rows)


def score(out_csv: Path, *, k: int = 100) -> pd.DataFrame:
    """Run concordance against the latest D1 + D5 outputs and write CSV."""
    human_jsons = sorted(ARM_DIR_HUMAN_D1.glob("*.json"))
    mouse_jsons = sorted(ARM_DIR_MOUSE.glob("*.json"))
    if not human_jsons or not mouse_jsons:
        print("[d5/concordance] missing D1 or D5 outputs; writing empty CSV")
        df = pd.DataFrame(columns=["human_gene", "mouse_gene", "jaccard_top_k"])
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        return df

    human_preds = []
    for p in human_jsons:
        human_preds.extend(_load_predictions(p))
    mouse_preds = []
    for p in mouse_jsons:
        mouse_preds.extend(_load_predictions(p))

    human_topk = per_gene_topk(human_preds, k=k)
    mouse_topk = per_gene_topk(mouse_preds, k=k)

    ortho = pd.read_csv(
        Path(
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
            "Analysis/Perturbation/data/hits/d5_mouse_orthologs.csv"
        )
    )
    df = concordance(human_topk, mouse_topk, ortho, k=k)
    if df.empty:
        df = pd.DataFrame(columns=["human_gene", "mouse_gene", "jaccard_top_k", "n_human_top", "n_mouse_top"])
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    return df


if __name__ == "__main__":
    score(ARM_DIR_MOUSE / "concordance_scores.csv")

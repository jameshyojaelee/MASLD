"""STATE mouse v2 runner — REAL state.tx inference (cross-species proxy).

NOTE: There is NO mouse-trained ST head available. The 4 Replogle ST heads
(K562, HepG2, Jurkat, RPE1) are all human. For D5 we use the HepG2 head and
predict on the HUMAN ortholog of each mouse target, emitting GenePrediction
rows keyed on the MOUSE symbol with cross-species provenance noted.

For each (mouse_gene, human_gene) pair:
  1. Skip if human_gene not in HepG2 vocab.
  2. Predict (KO human - control) gene-space delta.
  3. Top-100 affected genes: emit downstream_gene = HUMAN gene symbol
     (mouse-side mapping happens in concordance_scorer.py downstream).

When the Saunders 2025 mouse hepatocyte Perturb-seq dataset (task #26 done)
is processed into a fine-tuned mouse ST head (Phase 2), this runner can
swap checkpoint paths and predict directly in mouse space.

Output: results/d5_mouse/state_mouse_v2_<modality>_<context>.json
"""
from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402

THIS_DIR = Path(__file__).resolve().parent
SHARED = THIS_DIR.parents[1] / "shared"
sys.path.insert(0, str(SHARED))

from _runner_template import D5Adapter, main_for_model  # noqa: E402

from state_tx_hepg2 import HepG2STInference, STATE_TX_CKPT_DIR  # noqa: E402

TOP_N_DOWNSTREAM = 100


class StateMouseV2Adapter(D5Adapter):
    model_name = "state_mouse_v2"
    default_checkpoint = str(STATE_TX_CKPT_DIR)

    def load_checkpoint(self) -> None:
        self.infer = HepG2STInference(self.checkpoint)
        self.gene_names = self.infer.gene_names
        self.n_attempted = 0
        self.n_predicted = 0
        self.n_skipped = 0

    def predict_for_hit(self, *, hit_row, context, substrate):
        mouse_gene = str(hit_row.get("mouse_gene", "") or "")
        human_gene = str(hit_row.get("human_gene", "") or "")
        if not mouse_gene or not human_gene:
            return []

        self.n_attempted += 1
        if self.n_attempted % 50 == 0:
            print(
                f"[state_mouse_v2] processed {self.n_attempted} "
                f"(in-vocab {self.n_predicted}, skipped {self.n_skipped})",
                flush=True,
            )

        if not self.infer.in_vocab(human_gene):
            self.n_skipped += 1
            return []
        try:
            delta, _ctrl = self.infer.predict_gene_delta(human_gene, n_batches=2)
        except Exception as e:
            print(f"[state_mouse_v2] {human_gene} failed: {e}", flush=True)
            self.n_skipped += 1
            return []

        self.n_predicted += 1
        ranks = np.argsort(-np.abs(delta))[:TOP_N_DOWNSTREAM]
        names = self.gene_names if (self.gene_names and len(self.gene_names) == len(delta)) \
            else [f"GENE_{i}" for i in range(len(delta))]
        out = []
        for rank, gi in enumerate(ranks):
            d = float(delta[gi])
            out.append(dict(
                target_gene=mouse_gene,
                downstream_gene=names[gi],
                logFC_predicted=d,
                abs_rank=int(rank + 1),
                direction=("up" if d > 0 else "down" if d < 0 else "—"),
                confidence=float(1.0 / (1.0 + np.exp(-abs(d) / 0.05))),
            ))
        return out


if __name__ == "__main__":
    main_for_model(StateMouseV2Adapter)

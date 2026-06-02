"""STATE v2 runner for D1 mechanism — REAL state.tx inference.

Replaces ``state_runner.py``'s SE-600M embedding-similarity proxy with a
true HepG2 ST forward pass on Replogle-trained gene-KO head.

For each hit gene `g`:
  1. Skip if `g` not in HepG2 pert vocab (record in skip stats).
  2. Sample 64 control cells × 2 batches. Predict expression under
     non-targeting and KO(g); take the per-gene delta (mean across cells).
  3. Rank the 6,546 output genes by |delta|; emit top-100 as
     GenePrediction rows with logFC_predicted = delta, direction = sign(delta).
  4. Confidence = sigmoid(|delta| / 0.05) so larger deltas → higher confidence.

LIMITATIONS:
  * Replogle vocab is ~2k essential genes; D1 hits show ~7% coverage.
  * `delta` is in the model's internal expression scale (~clipped 0..14),
    NOT a true logFC. Treat magnitudes as RELATIVE.
  * Phase 2 will fine-tune on hepatocyte Perturb-seq (Saunders 2025).

OUT: results/d1_mechanism/state_v2_<modality>_<context>.json
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

from _runner_template import D1Adapter, TOP_N_DOWNSTREAM, main_for_model  # noqa: E402

from state_tx_hepg2 import HepG2STInference, STATE_TX_CKPT_DIR  # noqa: E402


class StateV2Adapter(D1Adapter):
    model_name = "state_v2"
    default_checkpoint = str(STATE_TX_CKPT_DIR)

    def load_checkpoint(self) -> None:
        self.infer = HepG2STInference(self.checkpoint)
        self.gene_names = self.infer.gene_names
        # Stats
        self.n_attempted = 0
        self.n_skipped = 0
        self.n_predicted = 0
        self._missing: set[str] = set()

    def predict_for_hit(self, *, hit_row, context, substrate):
        target = str(hit_row.get("gene", ""))
        if not target:
            return []
        self.n_attempted += 1
        if self.n_attempted % 20 == 0:
            print(
                f"[state_v2] processed {self.n_attempted} hits "
                f"(in-vocab {self.n_predicted}, skipped {self.n_skipped})",
                flush=True,
            )

        if not self.infer.in_vocab(target):
            self.n_skipped += 1
            self._missing.add(target)
            return []
        try:
            delta, _ctrl = self.infer.predict_gene_delta(target, n_batches=2)
        except Exception as e:
            print(f"[state_v2] {target} predict failed: {e}", flush=True)
            self.n_skipped += 1
            return []

        self.n_predicted += 1

        # Rank by |delta|. The HepG2 head outputs 6,546 gene-space predictions.
        # delta is clipped to [0, 14] internally so deltas can be ±-bounded.
        ranks = np.argsort(-np.abs(delta))[:TOP_N_DOWNSTREAM]
        if not self.gene_names or len(self.gene_names) != len(delta):
            # Embedding-space fallback (rare): synthesize gene labels.
            names = [f"GENE_{i}" for i in range(len(delta))]
        else:
            names = self.gene_names

        out = []
        for rank, gi in enumerate(ranks):
            d = float(delta[gi])
            out.append(
                dict(
                    target_gene=target,
                    downstream_gene=names[gi],
                    logFC_predicted=d,
                    abs_rank=int(rank + 1),
                    direction=("up" if d > 0 else "down" if d < 0 else "—"),
                    confidence=float(1.0 / (1.0 + np.exp(-abs(d) / 0.05))),
                )
            )
        return out


if __name__ == "__main__":
    main_for_model(StateV2Adapter)

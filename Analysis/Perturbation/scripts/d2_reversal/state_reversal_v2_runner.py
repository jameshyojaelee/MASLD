"""STATE reversal v2 runner — REAL state.tx inference.

D2 task: given a Diseased→Healthy reference signature (LFC per gene), score
each candidate KO by how much its predicted delta REVERSES the disease LFC.

Reversal score = -cosine(delta_KO, ref_LFC). Genes whose KO induces an
expression shift OPPOSITE to the disease direction score highest.

For each hit:
  1. Skip if not in HepG2 ST vocab.
  2. Sample 64 control cells × 2 batches. Predict KO(gene) vs non-targeting
     in gene space (6,546-dim).
  3. Align the ST output gene names with the reference signature's gene
     symbols; intersect.
  4. reversal_score = -cosine(delta[intersect], ref_LFC[intersect]).
  5. Emit a ReversalPrediction row per (gene, reference). reversal_rank is
     filled after the full hit list runs (sorted desc).

Pre-condition: D2 hits CSV has a `gene` column.

Output: results/d2_reversal/state_reversal_v2_<modality>_<context>_<ref>.json
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

from _runner_template import D2Adapter, main_for_model  # noqa: E402

from state_tx_hepg2 import HepG2STInference, STATE_TX_CKPT_DIR  # noqa: E402


def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na < 1e-12 or nb < 1e-12:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


class StateReversalV2Adapter(D2Adapter):
    model_name = "state_reversal_v2"
    default_checkpoint = str(STATE_TX_CKPT_DIR)

    def load_checkpoint(self) -> None:
        self.infer = HepG2STInference(self.checkpoint)
        self.gene_names = self.infer.gene_names
        # Build ref-LFC array aligned to model gene order.
        ref = self._ref_df
        if ref.empty or "gene" not in ref.columns or "LFC" not in ref.columns:
            print(
                f"[state_reversal_v2] reference {self.reference_signature} empty / missing columns",
                flush=True,
            )
            self._ref_aligned = None
        else:
            lfc_by_sym = dict(zip(ref["gene"].astype(str), ref["LFC"].astype(float)))
            aligned = np.array([lfc_by_sym.get(g, 0.0) for g in self.gene_names], dtype=np.float32)
            self._ref_aligned = aligned
            # how many genes have non-zero ref LFC in the intersection
            self._n_ref_overlap = int(np.sum(np.array([g in lfc_by_sym for g in self.gene_names])))
            print(
                f"[state_reversal_v2] ref={self.reference_signature}: "
                f"{self._n_ref_overlap}/{len(self.gene_names)} ST genes have ref LFC",
                flush=True,
            )
        self.n_attempted = 0
        self.n_predicted = 0
        self.n_skipped = 0

    def predict_for_hit(self, *, hit_row, context, substrate):
        gene = str(hit_row.get("gene", ""))
        if not gene:
            return []
        self.n_attempted += 1
        if self.n_attempted % 50 == 0:
            print(
                f"[state_reversal_v2] processed {self.n_attempted} "
                f"(in-vocab {self.n_predicted}, skipped {self.n_skipped})",
                flush=True,
            )

        if not self.infer.in_vocab(gene):
            self.n_skipped += 1
            return []
        if self._ref_aligned is None:
            self.n_skipped += 1
            return []

        try:
            delta, _ctrl = self.infer.predict_gene_delta(gene, n_batches=2)
        except Exception as e:
            print(f"[state_reversal_v2] {gene} failed: {e}", flush=True)
            self.n_skipped += 1
            return []

        # Use only genes where ref has non-zero LFC (avoids dilution by zero-padding).
        nz = self._ref_aligned != 0
        cos = _cosine(delta[nz], self._ref_aligned[nz])
        reversal_score = float(-cos)

        self.n_predicted += 1
        return [dict(
            gene=gene,
            reversal_score=reversal_score,
            reversal_rank=0,  # filled by consensus / post-processing
            stage_specific="—",
            cell_type="hepatocyte_progressor",
            reference_signature=self.reference_signature,
        )]


if __name__ == "__main__":
    main_for_model(StateReversalV2Adapter)

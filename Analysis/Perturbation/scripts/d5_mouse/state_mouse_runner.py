#!/usr/bin/env python
"""STATE mouse runner for D5. Subagent: mouse-runner-state.

Phase 1 zero-shot implementation: STATE ships a single SE-600M checkpoint
trained on a broad ~22k gene vocabulary (Arc Institute). There is no
dedicated mouse fine-tune; the SE-600M vocab already covers most mouse
symbols (lower-case, e.g. ``Mrpl32``) plus uppercase human symbols.

For each mouse hit row we:
  1. Look up the SE-600M embedding for the ``mouse_gene`` symbol.
  2. Look up the SE-600M embedding for the human ortholog (``human_gene``).
  3. Compute cross-species cosine similarity → per-row concordance score.
  4. Emit a single ``GenePrediction`` dict keyed on ``target_gene=mouse_gene``
     with ``downstream_gene=human_gene`` (the cross-species pairing partner).
     ``logFC_predicted`` carries the cosine similarity (range [-1, 1]); the
     downstream concordance scorer reads this back as the cross-species score.

Tensor handling follows ``d3_synergy/state_pair_runner.py``:
  - ``emb.detach().cpu().numpy()`` before ``np.asarray``
  - mean-pool across isoforms when ``arr.ndim >= 2``
  - CPU-forced via ``CUDA_VISIBLE_DEVICES=""`` (torch 2.5.1+cu124 vs sm_60
    cluster GPU mismatch)
  - first 3 ``get_gene_embedding`` errors logged with truncated traceback
"""
from __future__ import annotations

import os
import sys
import traceback
from pathlib import Path

# Force CPU before any torch/state imports — STATE ships torch 2.5.1+cu124
# which can't JIT-compile on the cluster's older GPU SM levels. SE-600M is
# 2.86 GB → fits comfortably on CPU; embedding lookups are cheap.
os.environ["CUDA_VISIBLE_DEVICES"] = ""

import numpy as np  # noqa: E402
import torch  # noqa: E402

from _runner_template import D5Adapter, main_for_model  # noqa: E402

STATE_SE600M_DIR = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
    "data/perturbation/checkpoints/state/SE-600M"
)


class StateMouseAdapter(D5Adapter):
    model_name = "state_mouse"
    default_checkpoint = str(STATE_SE600M_DIR)

    # ------------------------------------------------------------------
    # Model load
    # ------------------------------------------------------------------
    def load_checkpoint(self) -> None:
        from state.emb import Inference
        from omegaconf import OmegaConf

        ckpt = Path(self.checkpoint)
        if not ckpt.exists():
            raise FileNotFoundError(f"STATE checkpoint not at {ckpt}")

        cfg_path = ckpt / "config.yaml"
        cfg = OmegaConf.load(cfg_path) if cfg_path.exists() else None

        pe_path = ckpt / "protein_embeddings.pt"
        protein_embeds = (
            torch.load(pe_path, weights_only=False) if pe_path.exists() else None
        )

        self.model = Inference(cfg=cfg, protein_embeds=protein_embeds)

        # Lightning .ckpt — prefer epoch16 (better-trained), fallback to epoch4.
        for candidate in [ckpt / "se600m_epoch16.ckpt", ckpt / "se600m_epoch4.ckpt"]:
            if candidate.exists():
                self.model.load_model(str(candidate))
                print(f"[state_mouse] loaded weights from {candidate.name}")
                break
        else:
            raise FileNotFoundError(f"No .ckpt weights under {ckpt}")

        print(f"[state_mouse] STATE SE-600M loaded from {ckpt.name}")

        # Skip / error tracking
        self._n_attempted = 0
        self._n_skipped = 0
        self._n_emb_errors = 0
        self._missing_mouse: set[str] = set()
        self._missing_human: set[str] = set()
        # Embedding cache — STATE lookups are cheap but ortholog list reuses
        # the same human symbols (e.g., MT-* family). 1180 rows × 2 lookups.
        self._emb_cache: dict[str, np.ndarray | None] = {}

    # ------------------------------------------------------------------
    # Embedding helper
    # ------------------------------------------------------------------
    def _get_embedding(self, gene: str) -> np.ndarray | None:
        """Return SE-600M embedding for ``gene`` (mean across isoforms) or None."""
        if not isinstance(gene, str) or not gene:
            return None
        if gene in self._emb_cache:
            return self._emb_cache[gene]
        try:
            emb = self.model.get_gene_embedding(gene)
        except Exception as e:
            if self._n_emb_errors < 3:
                tb = traceback.format_exc().splitlines()
                print(
                    f"[state_mouse] get_gene_embedding({gene!r}) raised "
                    f"{type(e).__name__}: {str(e)[:120]}"
                )
                if tb:
                    print(f"  -> {tb[-1][:200]}")
                self._n_emb_errors += 1
            self._emb_cache[gene] = None
            return None

        if emb is None:
            self._emb_cache[gene] = None
            return None
        if hasattr(emb, "detach"):
            emb = emb.detach().cpu().numpy()
        arr = np.asarray(emb)
        if arr.ndim >= 2:
            arr = arr.mean(axis=0)
        arr = arr.ravel().astype(np.float32, copy=False)
        self._emb_cache[gene] = arr
        return arr

    @staticmethod
    def _cosine(a: np.ndarray, b: np.ndarray) -> float:
        na = float(np.linalg.norm(a)) + 1e-12
        nb = float(np.linalg.norm(b)) + 1e-12
        return float(np.dot(a, b) / (na * nb))

    # ------------------------------------------------------------------
    # Prediction
    # ------------------------------------------------------------------
    def predict_for_hit(self, *, hit_row, context, substrate):
        mouse_gene = hit_row.get("mouse_gene")
        human_gene = hit_row.get("human_gene")
        if not isinstance(mouse_gene, str) or not mouse_gene:
            return []
        if not isinstance(human_gene, str) or not human_gene:
            return []

        self._n_attempted += 1
        if self._n_attempted % 100 == 0:
            print(
                f"[state_mouse] processed {self._n_attempted}, "
                f"skipped {self._n_skipped} "
                f"(missing_mouse={len(self._missing_mouse)}, "
                f"missing_human={len(self._missing_human)})"
            )

        mouse_emb = self._get_embedding(mouse_gene)
        human_emb = self._get_embedding(human_gene)

        if mouse_emb is None:
            self._missing_mouse.add(mouse_gene)
        if human_emb is None:
            self._missing_human.add(human_gene)

        if mouse_emb is None or human_emb is None:
            self._n_skipped += 1
            return []

        cos = self._cosine(mouse_emb, human_emb)
        # Bound cos for confidence (cos ∈ [-1, 1] → [0, 1]).
        confidence = max(0.0, min(1.0, (cos + 1.0) / 2.0))
        # Direction: STATE embeddings are sign-agnostic; concordance is
        # bidirectional. Schema requires up/down/either/—.
        direction = "either"

        return [
            dict(
                target_gene=str(mouse_gene),
                # Use the cross-species partner as the "downstream" gene —
                # this is what concordance_scorer.py reads back.
                downstream_gene=str(human_gene),
                # logFC_predicted carries the raw cosine similarity (signed).
                # Downstream scorers convert to per-row concordance.
                logFC_predicted=float(cos),
                abs_rank=1,
                direction=direction,
                confidence=float(confidence),
            )
        ]


if __name__ == "__main__":
    try:
        main_for_model(StateMouseAdapter)
    finally:
        # Defensive: surface skip-stats even if main raises during write-out.
        # Useful for debugging vocab coverage early in a fresh run.
        adapter_dbg = sys.modules.get("__main__")
        if adapter_dbg is not None:
            print("[state_mouse] run finished.")

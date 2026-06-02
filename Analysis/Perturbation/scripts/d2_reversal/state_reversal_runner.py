#!/usr/bin/env python
"""STATE (Arc Institute SE-600M) runner for D2 reversal — zero-shot embedding proxy.

D2 task: rank ALL atlas genes by predicted reversal of Diseased→Healthy
hepatocyte signature. Reversal score for gene g vs reference signature R:

    delta_g       = emb(g) - pop_avg
    disease_axis  = sum_{h in R}  ref_LFC[h] * emb(h)     # LFC-weighted sum
    reversal_score(g) = - cosine(delta_g, disease_axis)

Interpretation: when delta_g points *opposite* the embedding-space disease axis
(positive ref_LFC direction = upregulated in disease), cosine is negative and
``reversal_score`` is positive → gene g sits on the "healthy" side of the axis
and its perturbation is predicted to push the system back toward healthy.

This is the SAME limitation acknowledged in D1/D5: SE-600M is an *embedding*
model; without state.tx (real perturbation inference) we cannot simulate the
KO state directly. The LFC-weighted axis is the standard zero-shot proxy and
validates the D2 pipeline end-to-end. State.tx, when ready, will replace
``predict_for_hit`` with real predicted-Δ-expression.

Tensor handling follows ``d1_mechanism/state_runner.py`` +
``d5_mouse/state_mouse_runner.py``:
  - CPU-forced via ``CUDA_VISIBLE_DEVICES=""`` (torch 2.5.1+cu124 vs cluster
    SM mismatch)
  - ``emb.detach().cpu().numpy()`` before ``np.asarray``
  - mean-pool across isoforms when ``arr.ndim >= 2``
  - first 3 ``get_gene_embedding`` errors logged with truncated traceback

CLI:
    python state_reversal_runner.py --modality zero_shot --context all --reference ref_a
"""
from __future__ import annotations

import os
import sys
import traceback
from pathlib import Path

# Force CPU before any torch/state imports — SE-600M is 2.86 GB → fits on CPU.
os.environ["CUDA_VISIBLE_DEVICES"] = ""

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

from _runner_template import D2Adapter, main_for_model  # noqa: E402

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
STATE_SE600M_DIR = PROJECT_ROOT / "data/perturbation/checkpoints/state/SE-600M"


class StateReversalAdapter(D2Adapter):
    model_name = "state_reversal"
    default_checkpoint = str(STATE_SE600M_DIR)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Initialise stats + cache up-front so dry-run (which skips
        # load_checkpoint) doesn't AttributeError on first predict_for_hit.
        self._n_attempted = 0
        self._n_skipped = 0
        self._n_emb_errors = 0
        self._missing_genes: set[str] = set()
        self._emb_cache: dict[str, np.ndarray | None] = {}
        self._pop_avg = np.zeros(2048, dtype=np.float32)
        self._disease_axis = np.zeros(2048, dtype=np.float32)
        self._axis_norm = 0.0
        self._n_ref_genes_used = 0

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
                print(f"[state_reversal] loaded weights from {candidate.name}")
                break
        else:
            raise FileNotFoundError(f"No .ckpt weights under {ckpt}")

        print(f"[state_reversal] STATE SE-600M loaded from {ckpt.name}")

        # Build the embedding-space disease axis from the reference signature
        self._pop_avg, self._disease_axis, self._axis_norm = self._build_disease_axis()
        print(
            f"[state_reversal] disease axis built from ref={self.reference_signature} "
            f"({getattr(self, '_n_ref_genes_used', 0)} ref genes embedded; "
            f"axis_norm={self._axis_norm:.4f})"
        )

    # ------------------------------------------------------------------
    # Disease-axis construction
    # ------------------------------------------------------------------
    def _build_disease_axis(self) -> tuple[np.ndarray, np.ndarray, float]:
        """Construct (pop_avg, disease_axis, axis_norm) in SE-600M embedding space.

        Strategy:
          1. From the reference signature, restrict to genes with a finite LFC.
          2. (Optional) restrict to high-confidence axis genes by |LFC| top-K
             to avoid noisy near-zero coefficients diluting the axis.
          3. Look up embeddings for those genes; drop misses.
          4. pop_avg     = mean of retrieved embeddings.
          5. disease_axis= LFC-weighted sum of (emb - pop_avg) over kept genes.

        Returns
        -------
        (pop_avg, disease_axis, axis_norm)
            pop_avg, disease_axis: 1-D np.float32 arrays of dim 2048.
            axis_norm: float, used for cosine denominator + sanity-check log.
        """
        ref = self._ref_df.copy()
        if ref.empty or "LFC" not in ref.columns or "gene" not in ref.columns:
            print(
                f"[state_reversal] WARN: reference {self.reference_signature} empty / missing "
                f"LFC|gene columns → returning zero-axis"
            )
            self._n_ref_genes_used = 0
            return np.zeros(2048, dtype=np.float32), np.zeros(2048, dtype=np.float32), 0.0

        ref = ref.dropna(subset=["LFC", "gene"]).copy()
        ref["gene"] = ref["gene"].astype(str)
        # Drop Ensembl-only ref_b rows that don't map to a symbol (they remain
        # ENSG… in the gene column — STATE expects gene symbols and will miss).
        ref = ref[~ref["gene"].str.startswith("ENSG")]

        # Cap the axis to top-K by |LFC| to reduce noise from low-signal genes.
        # 2,500 is roughly the size of ref_a / ref_c's HVG set so this is a no-op
        # for them; for ref_b (34k) it filters to the strongest 2,500 by |LFC|.
        ref["abs_lfc"] = ref["LFC"].abs()
        ref = ref.sort_values("abs_lfc", ascending=False).head(2500)

        kept_emb: list[np.ndarray] = []
        kept_lfc: list[float] = []
        for g, lfc in zip(ref["gene"], ref["LFC"]):
            e = self._get_embedding(g)
            if e is None:
                continue
            kept_emb.append(e)
            kept_lfc.append(float(lfc))

        if not kept_emb:
            print(f"[state_reversal] WARN: 0 ref genes resolved to embeddings")
            self._n_ref_genes_used = 0
            return np.zeros(2048, dtype=np.float32), np.zeros(2048, dtype=np.float32), 0.0

        E = np.stack(kept_emb).astype(np.float32)  # (N, 2048)
        L = np.asarray(kept_lfc, dtype=np.float32)
        pop_avg = E.mean(axis=0)
        # Disease axis = LFC-weighted sum of mean-centered embeddings.
        # Centering matters: without it the bulk of the axis is just pop_avg.
        axis = (L[:, None] * (E - pop_avg)).sum(axis=0)
        norm = float(np.linalg.norm(axis))
        self._n_ref_genes_used = len(kept_emb)
        return pop_avg, axis, norm

    # ------------------------------------------------------------------
    # Embedding helper
    # ------------------------------------------------------------------
    def _get_embedding(self, gene: str) -> np.ndarray | None:
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
                    f"[state_reversal] get_gene_embedding({gene!r}) raised "
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
        gene = hit_row.get("gene")
        if not isinstance(gene, str) or not gene:
            return []

        self._n_attempted += 1
        if self._n_attempted % 100 == 0:
            print(
                f"[state_reversal] processed {self._n_attempted}, skipped "
                f"{self._n_skipped} (missing={len(self._missing_genes)})"
            )

        if self._axis_norm <= 0.0:
            # No usable axis → predictions are meaningless; emit nothing.
            self._n_skipped += 1
            return []

        emb = self._get_embedding(gene)
        if emb is None:
            self._missing_genes.add(gene)
            self._n_skipped += 1
            return []

        delta = emb - self._pop_avg
        cos = self._cosine(delta, self._disease_axis)
        # Reversal score: -cosine; high positive = pushes opposite to disease axis.
        reversal_score = float(-cos)

        return [
            dict(
                gene=str(gene),
                reversal_score=reversal_score,
                # rank assigned by template post-processor; placeholder here
                reversal_rank=0,
                stage_specific="—",
                cell_type="hepatocyte_progressor",
                reference_signature=self.reference_signature,
            )
        ]


if __name__ == "__main__":
    main_for_model(StateReversalAdapter)

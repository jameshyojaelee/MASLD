#!/usr/bin/env python
"""STATE (Arc Institute SE-600M) runner for D1 mechanism.

D1 task: for each hit gene, predict top-K downstream genes affected by its
perturbation. Phase 1 zero-shot proxy:

  1. Use STATE SE-600M to get a 2048-dim embedding for the hit gene
     (mean-pooled across isoforms).
  2. Pre-compute embeddings for a candidate downstream pool: the union of
     hit genes + curated MASLD-relevant genes from the multi-evidence atlas.
  3. For each hit, compute cosine similarity vs all candidates; the top-K
     nearest neighbours in embedding space are emitted as "downstream"
     predictions. This is a functional-similarity proxy (genes that move
     together when one is perturbed are co-functional / co-regulated).

Phase 2 will replace this with state.tx (real perturbation inference) once
the state-tx-implementer agent finishes.

Tensor handling follows ``d3_synergy/state_pair_runner.py`` +
``d5_mouse/state_mouse_runner.py``:
  - CPU-forced via ``CUDA_VISIBLE_DEVICES=""`` (torch 2.5.1+cu124 vs cluster
    SM mismatch)
  - ``emb.detach().cpu().numpy()`` before ``np.asarray``
  - mean-pool across isoforms when ``arr.ndim >= 2``
  - first 3 ``get_gene_embedding`` errors logged with truncated traceback

CLI:
    python state_runner.py --modality zero_shot --context all
    python state_runner.py --modality zero_shot --context all --dry-run
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

from _runner_template import D1Adapter, TOP_N_DOWNSTREAM, main_for_model  # noqa: E402

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
STATE_SE600M_DIR = PROJECT_ROOT / "data/perturbation/checkpoints/state/SE-600M"
HITS_PATH = PROJECT_ROOT / "Analysis/Perturbation/data/hits/d1_mechanism_hits.csv"
ATLAS_DREAM_PATH = PROJECT_ROOT / (
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/"
    "dream_results_ashr.csv"
)


class StateAdapter(D1Adapter):
    model_name = "state"
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
                print(f"[state] loaded weights from {candidate.name}")
                break
        else:
            raise FileNotFoundError(f"No .ckpt weights under {ckpt}")

        print(f"[state] STATE SE-600M loaded from {ckpt.name}")

        # Skip / error tracking
        self._n_attempted = 0
        self._n_skipped = 0
        self._n_emb_errors = 0
        self._missing_genes: set[str] = set()
        self._emb_cache: dict[str, np.ndarray | None] = {}

        # Build candidate downstream pool
        self._candidate_genes, self._candidate_embs = self._build_candidate_pool()
        print(
            f"[state] candidate downstream pool: {len(self._candidate_genes)} genes "
            f"(emb matrix {self._candidate_embs.shape})"
        )

    # ------------------------------------------------------------------
    # Candidate pool construction
    # ------------------------------------------------------------------
    def _build_candidate_pool(self) -> tuple[list[str], np.ndarray]:
        """Return (gene_list, normalized_embedding_matrix) for downstream candidates.

        Pool = union of (a) D1 hits and (b) MASLD-DEG top ~5K from dream atlas.
        We pre-compute & L2-normalize embeddings for fast cosine top-K.
        """
        candidates: set[str] = set()

        # (a) D1 hits themselves — guarantees every hit has potential downstream
        # candidates that aren't itself.
        if HITS_PATH.exists():
            try:
                hits = pd.read_csv(HITS_PATH)
                if "gene" in hits.columns:
                    candidates.update(
                        str(g) for g in hits["gene"].dropna().astype(str).tolist()
                    )
            except Exception as e:
                print(f"[state] WARN: failed reading hits ({e})")

        # (b) Dream-mega DEG top ~5K by |stat|. Filter to padj<0.1 first.
        if ATLAS_DREAM_PATH.exists():
            try:
                dr = pd.read_csv(ATLAS_DREAM_PATH, low_memory=False)
                # Prefer 'symbol' over 'gene' — dream atlas 'gene' col is Ensembl ID
                # (ENSG…) while STATE expects gene symbols.
                gene_col = next(
                    (c for c in dr.columns if c.lower() in ("symbol", "gene_symbol", "hgnc_symbol")),
                    None,
                )
                if gene_col is None:
                    # fallback to 'gene' if no symbol col present
                    gene_col = next(
                        (c for c in dr.columns if c.lower() in ("gene",)),
                        None,
                    )
                stat_col = next(
                    (c for c in dr.columns if c.lower() in ("t", "tstat", "stat", "logfc", "log2foldchange")),
                    None,
                )
                padj_col = next(
                    (c for c in dr.columns if c.lower() in ("padj", "adj.p.val", "fdr")),
                    None,
                )
                if gene_col and stat_col:
                    sub = dr[[gene_col, stat_col] + ([padj_col] if padj_col else [])].dropna()
                    if padj_col:
                        sub = sub[sub[padj_col] < 0.1]
                    sub = sub.copy()
                    sub["_abs"] = sub[stat_col].abs()
                    top = sub.sort_values("_abs", ascending=False).head(5000)
                    candidates.update(top[gene_col].astype(str).tolist())
            except Exception as e:
                print(f"[state] WARN: failed reading dream atlas ({e})")

        if not candidates:
            print("[state] WARN: empty candidate pool; falling back to hits only")

        # Get embedding for each candidate, drop misses
        gene_list: list[str] = []
        emb_list: list[np.ndarray] = []
        for g in sorted(candidates):
            e = self._get_embedding(g)
            if e is None:
                continue
            # L2-normalize so cosine = dot product
            n = float(np.linalg.norm(e)) + 1e-12
            emb_list.append(e / n)
            gene_list.append(g)

        if not emb_list:
            return [], np.zeros((0, 0), dtype=np.float32)

        return gene_list, np.stack(emb_list).astype(np.float32)

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
                    f"[state] get_gene_embedding({gene!r}) raised "
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

    # ------------------------------------------------------------------
    # Prediction
    # ------------------------------------------------------------------
    def predict_for_hit(self, *, hit_row, context, substrate):
        target = hit_row.get("gene")
        if not isinstance(target, str) or not target:
            return []

        self._n_attempted += 1
        if self._n_attempted % 50 == 0:
            print(
                f"[state] processed {self._n_attempted}, skipped {self._n_skipped} "
                f"(missing={len(self._missing_genes)})"
            )

        # Need a non-empty candidate pool
        if self._candidate_embs.shape[0] == 0:
            self._n_skipped += 1
            return []

        target_emb = self._get_embedding(target)
        if target_emb is None:
            self._missing_genes.add(target)
            self._n_skipped += 1
            return []

        # L2-normalize target; compute cosine against pre-normalized candidates
        t_norm = target_emb / (float(np.linalg.norm(target_emb)) + 1e-12)
        sims = self._candidate_embs @ t_norm  # (N_candidates,) cosines in [-1, 1]

        # Mask self (target may be in candidate pool)
        mask = np.array(
            [g != target for g in self._candidate_genes], dtype=bool
        )
        sims_masked = np.where(mask, sims, -np.inf)

        # Top-K by similarity magnitude (we use signed similarity here — high cos = closely
        # co-functional and likely co-up-regulated under perturbation; STATE proxy is
        # direction-agnostic so we tag direction as "either").
        k = min(TOP_N_DOWNSTREAM, int(mask.sum()))
        if k <= 0:
            self._n_skipped += 1
            return []
        top_idx = np.argpartition(-sims_masked, k - 1)[:k]
        # Sort the top-K block by descending sim
        top_idx = top_idx[np.argsort(-sims_masked[top_idx])]

        rows: list[dict] = []
        for rank, idx in enumerate(top_idx, start=1):
            cos = float(sims_masked[idx])
            # cos ∈ [-1, 1] → confidence in [0, 1]
            confidence = max(0.0, min(1.0, (cos + 1.0) / 2.0))
            rows.append(
                dict(
                    target_gene=str(target),
                    downstream_gene=str(self._candidate_genes[idx]),
                    # logFC_predicted carries the cosine similarity (signed).
                    # Phase 2 (state.tx) will replace with real predicted logFC.
                    logFC_predicted=float(cos),
                    abs_rank=int(rank),
                    direction="either",
                    confidence=float(confidence),
                )
            )
        return rows


if __name__ == "__main__":
    main_for_model(StateAdapter)

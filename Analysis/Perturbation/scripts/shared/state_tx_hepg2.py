"""Shared state.tx HepG2 ST head loader + inference helper.

Used by:
  * d3_synergy/state_pair_v2_runner.py (combinatorial KO synergy)
  * d1_mechanism/state_v2_runner.py (downstream effects)
  * d2_reversal/state_reversal_v2_runner.py (signature reversal)

The checkpoint is the Arc Institute ST-SE-Replogle HepG2-0.99 head
(2024 essential-gene KO vocab, 6546-gene output, X_state embedding input).

KEY LIMITATION: vocab = 2024 Replogle essential genes; coverage of
MASLD-relevant hits is small (D1 ~7%, D3 tier1 ~4%, D5 ~9%). Out-of-vocab
hits are skipped; counts recorded. Phase 2 fine-tunes on hepatocyte
Perturb-seq (Saunders 2025) to broaden vocab — see task #16.
"""
from __future__ import annotations

import os
import pickle
from pathlib import Path
from typing import Any

import numpy as np

# Force CPU before any cuda-touching import. Cluster GPUs predate sm_80 builds.
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
STATE_TX_CKPT_DIR = (
    PROJECT_ROOT
    / "data/perturbation/checkpoints/state/st-se-replogle-full/hepg2_0.99"
)


def patch_llama_validate() -> None:
    """Disable HF strict LlamaConfig validation (one-shot, idempotent).

    HepG2 ST head trained with hidden_size=328, num_attention_heads=12 (not a
    multiple). Newer transformers raises in validate(); the arch itself runs
    fine because head_dim=64 is set independently. Patch validate to no-op.
    """
    from transformers.models.llama.configuration_llama import LlamaConfig

    if getattr(LlamaConfig, "_state_tx_patched", False):
        return
    orig = LlamaConfig.validate

    def lenient(self):
        try:
            orig(self)
        except Exception:
            pass

    LlamaConfig.validate = lenient
    LlamaConfig._state_tx_patched = True


class HepG2STInference:
    """Wrap StateTransitionPerturbationModel for repeated inference calls."""

    CELL_SET_LEN = 64

    def __init__(
        self,
        ckpt_dir: str | Path = STATE_TX_CKPT_DIR,
        n_ctrl_cache: int | None = None,
    ):
        import torch
        import anndata as ad

        patch_llama_validate()

        from state.tx.models.state_transition import StateTransitionPerturbationModel

        ckpt_dir = Path(ckpt_dir)
        if not ckpt_dir.exists():
            raise FileNotFoundError(f"Missing ST checkpoint dir: {ckpt_dir}")

        print(f"[state_tx] loading from {ckpt_dir}", flush=True)
        self.model = StateTransitionPerturbationModel.load_from_checkpoint(
            ckpt_dir / "checkpoints" / "final.ckpt",
            map_location="cpu",
            strict=False,
        )
        self.model.eval()

        self.pert_map: dict[str, Any] = torch.load(
            ckpt_dir / "pert_onehot_map.pt", weights_only=False
        )
        with open(ckpt_dir / "batch_onehot_map.pkl", "rb") as fh:
            self.batch_map = pickle.load(fh)
        with open(ckpt_dir / "var_dims.pkl", "rb") as fh:
            self.var_dims = pickle.load(fh)
        self.gene_names: list[str] = list(self.var_dims.get("gene_names", []))

        # Control basal embeddings (X_state) from packaged eval anndata.
        eval_real = ckpt_dir / "eval_best.ckpt" / "adata_real.h5ad"
        a = ad.read_h5ad(eval_real)
        ctrl_mask = (a.obs["gene"] == "non-targeting").values
        self._ctrl_xs = a.obsm["X_state"][ctrl_mask].astype(np.float32)
        self._ctrl_gem = a.obs["gem_group"][ctrl_mask].astype(str).values
        if n_ctrl_cache is not None and n_ctrl_cache < len(self._ctrl_xs):
            rng_cache = np.random.RandomState(0)
            sel = rng_cache.choice(
                len(self._ctrl_xs), size=n_ctrl_cache, replace=False
            )
            self._ctrl_xs = self._ctrl_xs[sel]
            self._ctrl_gem = self._ctrl_gem[sel]

        self._rng = np.random.RandomState(42)
        self._torch = torch
        self.nontarget_oh = self.pert_map["non-targeting"].float()
        print(
            f"[state_tx] loaded: vocab={len(self.pert_map)} genes={len(self.gene_names)} "
            f"ctrl_cells={self._ctrl_xs.shape}",
            flush=True,
        )

    # ------------------------------------------------------------------
    def _batch_idx(self, gem: str) -> int:
        v = self.batch_map.get(gem)
        if v is None:
            return 0
        import torch
        if torch.is_tensor(v) and v.ndim == 1:
            return int(torch.argmax(v).item())
        if isinstance(v, (int, np.integer)):
            return int(v)
        return 0

    def in_vocab(self, gene: str) -> bool:
        return gene in self.pert_map

    def sample_basal(self, n: int | None = None, seed_offset: int = 0):
        """Return (basal[float32,N,2058], batch_idx[N], gem[str,N])."""
        torch = self._torch
        n = n or self.CELL_SET_LEN
        rng = (
            np.random.RandomState(42 + seed_offset)
            if seed_offset
            else self._rng
        )
        sel = rng.choice(len(self._ctrl_xs), size=n, replace=False)
        basal = self._ctrl_xs[sel]
        gems = self._ctrl_gem[sel]
        basal_t = torch.tensor(basal, dtype=torch.float32)
        batch_idx_t = torch.tensor(
            [self._batch_idx(g) for g in gems], dtype=torch.long
        )
        return basal_t, batch_idx_t

    def forward(self, pert_oh, basal_t, batch_idx_t) -> np.ndarray:
        """Return gene-space prediction [N, n_genes] for the supplied pert one-hot.

        Falls back to embedding-space [N, input_dim] if gene-space counts head
        isn't populated (rare; HepG2 head has output_space=all).
        """
        torch = self._torch
        with torch.no_grad():
            batch = {
                "ctrl_cell_emb": basal_t,
                "pert_emb": pert_oh.float().unsqueeze(0).repeat(basal_t.shape[0], 1),
                "pert_name": ["pert"] * basal_t.shape[0],
                "batch": batch_idx_t,
            }
            out = self.model.predict_step(batch, batch_idx=0, padded=False)
        preds = out["preds"].detach().cpu().numpy()
        counts = None
        if "pert_cell_counts_preds" in out and out["pert_cell_counts_preds"] is not None:
            counts = out["pert_cell_counts_preds"].detach().cpu().numpy()
        return counts if counts is not None else preds

    def predict_gene_delta(
        self, gene: str, n_batches: int = 2
    ) -> tuple[np.ndarray, np.ndarray]:
        """Predict mean (perturbed - control) gene-space response for `gene`.

        Returns (delta[n_genes], ctrl[n_genes]) — both averaged across batches × cells.
        Raises KeyError if gene not in pert_map.
        """
        if gene not in self.pert_map:
            raise KeyError(gene)
        pert_oh = self.pert_map[gene].float()

        delta = None
        ctrl_sum = None
        for b in range(n_batches):
            basal_t, batch_idx_t = self.sample_basal(seed_offset=b)
            ctrl_pred = self.forward(self.nontarget_oh, basal_t, batch_idx_t).mean(axis=0)
            ko_pred = self.forward(pert_oh, basal_t, batch_idx_t).mean(axis=0)
            d = ko_pred - ctrl_pred
            if delta is None:
                delta = d.copy()
                ctrl_sum = ctrl_pred.copy()
            else:
                delta += d
                ctrl_sum += ctrl_pred
        return delta / n_batches, ctrl_sum / n_batches

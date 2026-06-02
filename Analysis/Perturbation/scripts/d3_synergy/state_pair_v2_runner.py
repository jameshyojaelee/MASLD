"""STATE-pair v2 runner — REAL state.tx inference.

Replaces state_pair_runner.py's SE-600M embedding-distance heuristic with a
real ST (state-transition) forward pass on Replogle HepG2-trained head.

Approach:
  1. Load `StateTransitionPerturbationModel` (HepG2 0.99 head) + control
     X_state embeddings + pert one-hot map.
  2. For each pair (g1, g2):
      - Skip if either gene is NOT in the model's perturbation vocab.
      - Run 4 forward passes on 64 control cells:
          * pert = non-targeting       (baseline)
          * pert = one_hot(g1)         (single KO 1)
          * pert = one_hot(g2)         (single KO 2)
          * pert = one_hot(g1)+one_hot(g2)   (synthetic 2-hot double KO)
      - Aggregate per-pair (mean over the 64 cells in gene-space output).
      - delta_g1 = g1_pred - baseline
      - delta_g2 = g2_pred - baseline
      - delta_combo = double_pred - baseline
      - delta_additive = delta_g1 + delta_g2
      - synergy_magnitude = ||delta_combo - delta_additive|| / (||delta_additive|| + eps)
      - sigma_above_additive = 2 * tanh(synergy_magnitude)   # 0..2 scale
      - synergy_class = "synergistic" if synergy_magnitude > 0.3 else
                        "additive"     if synergy_magnitude < 0.1 else
                        "additive"     (mid-band defaults to additive)

LIMITATIONS:
  * Replogle HepG2 vocab is 2,024 essential-housekeeping genes.
    D3 tier1 (top-100 convergence) shares only 4 genes → only 6/4,950 pairs
    have both genes in vocab. Tier2 (top-770 convergence) shares 81 →
    3,240/296,065 pairs with full coverage.
  * The 2-hot double-KO encoding is OOD vs training (which saw only
    single-hot one-hots). Treat synergy magnitudes as RELATIVE rankings;
    DO NOT interpret absolute sigma values as calibrated probabilities.
  * For Phase 2 the right move is to fine-tune a ST head on hepatocyte
    Perturb-seq (Saunders 2025) — see task #16 in the campaign tracker.

OUT: Analysis/Perturbation/results/d3_synergy/state_pair_v2_zero_shot_<context>.json
"""
from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

# Force CPU before any CUDA-touching imports.
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
warnings.filterwarnings("ignore")

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from _runner_template import D3Adapter, main_for_model  # noqa: E402


STATE_TX_CKPT_DIR = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
    "data/perturbation/checkpoints/state/st-se-replogle-full/hepg2_0.99"
)


def _patch_llama_validate() -> None:
    """LlamaConfig (HF transformers) added strict validation in newer versions.

    The HepG2 ST head uses non-standard dims (hidden_size=328, n_heads=12;
    not a multiple) — the saved config fails strict validation. The model
    architecture itself works because head_dim is set independently. Patch
    the validator to log + skip instead of raising.
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


class StatePairV2Adapter(D3Adapter):
    """Real state.tx HepG2 forward pass for combinatorial KO synergy."""

    model_name = "state_pair_v2"
    default_checkpoint = str(STATE_TX_CKPT_DIR)

    # Number of control cells per forward pass. HepG2 head was trained with
    # cell_sentence_len=64, so this matches model expectations.
    CELL_SET_LEN = 64

    # Cap on number of distinct ctrl-cell mini-batches per pair. Each is one
    # forward pass; more = more averaging at higher compute cost.
    N_CTRL_BATCHES = 2

    def load_checkpoint(self) -> None:
        import torch
        import pickle
        import anndata as ad

        _patch_llama_validate()

        from state.tx.models.state_transition import StateTransitionPerturbationModel

        ckpt_dir = Path(self.checkpoint)
        if not ckpt_dir.exists():
            raise FileNotFoundError(f"Missing ST checkpoint dir: {ckpt_dir}")

        ckpt_file = ckpt_dir / "checkpoints" / "final.ckpt"
        print(f"[state_pair_v2] loading checkpoint: {ckpt_file}", flush=True)
        self._model = StateTransitionPerturbationModel.load_from_checkpoint(
            ckpt_file, map_location="cpu", strict=False
        )
        self._model.eval()
        print(
            f"[state_pair_v2] loaded; cell_sentence_len={self._model.cell_sentence_len}, "
            f"output_space={self._model.output_space}",
            flush=True,
        )

        self._pert_map = torch.load(
            ckpt_dir / "pert_onehot_map.pt", weights_only=False
        )
        with open(ckpt_dir / "batch_onehot_map.pkl", "rb") as fh:
            self._batch_map = pickle.load(fh)

        # Control basal embeddings — pull from packaged eval anndata.
        eval_real = ckpt_dir / "eval_best.ckpt" / "adata_real.h5ad"
        a = ad.read_h5ad(eval_real)
        ctrl_mask = (a.obs["gene"] == "non-targeting").values
        self._ctrl_xs = a.obsm["X_state"][ctrl_mask].astype(np.float32)
        self._ctrl_gem = a.obs["gem_group"][ctrl_mask].astype(str).values
        print(
            f"[state_pair_v2] vocab={len(self._pert_map)}, batch_vocab={len(self._batch_map)}, "
            f"ctrl_cells={self._ctrl_xs.shape}",
            flush=True,
        )

        self._torch = torch  # cached for predict_for_hit
        self._rng = np.random.RandomState(42)
        self._nontarget_oh = self._pert_map["non-targeting"].float()

        # Bookkeeping
        self.n_attempted = 0
        self.n_skipped = 0
        self.n_in_vocab = 0
        self._missing_genes: set[str] = set()

    # ---- internals -------------------------------------------------------
    def _batch_idx_for_gem(self, gem: str) -> int:
        v = self._batch_map.get(gem)
        if v is None:
            return 0
        import torch
        if torch.is_tensor(v) and v.ndim == 1:
            return int(torch.argmax(v).item())
        if isinstance(v, (int, np.integer)):
            return int(v)
        return 0

    def _forward(self, pert_oh, basal_t, batch_idx_t):
        """One forward pass; return (preds_emb, preds_gene_counts).

        preds_emb : ndarray [N, input_dim] — the ST decoder's residual output.
        preds_gene_counts : ndarray [N, gene_dim] or None — the gene-space
            decoder output (only present when output_space != embedding).
        """
        torch = self._torch
        with torch.no_grad():
            batch = {
                "ctrl_cell_emb": basal_t,
                "pert_emb": pert_oh.float().unsqueeze(0).repeat(basal_t.shape[0], 1),
                "pert_name": ["pert"] * basal_t.shape[0],
                "batch": batch_idx_t,
            }
            out = self._model.predict_step(batch, batch_idx=0, padded=False)
        preds = out["preds"].detach().cpu().numpy()
        counts = None
        if "pert_cell_counts_preds" in out and out["pert_cell_counts_preds"] is not None:
            counts = out["pert_cell_counts_preds"].detach().cpu().numpy()
        # Prefer counts (gene-space) for biology; fallback to emb-space.
        return counts if counts is not None else preds

    def _score_pair(self, genes: list[str]) -> dict:
        """Score a single combinatorial hit (2-4 genes).

        Returns a dict of metrics, or None if any gene is OOV.
        """
        import torch

        for g in genes:
            if g not in self._pert_map:
                self._missing_genes.add(g)
                return None

        ohs = [self._pert_map[g].float() for g in genes]

        # Aggregate across N_CTRL_BATCHES of 64 ctrl cells each (then mean over cells).
        delta_singles = [np.zeros(0, dtype=np.float32) for _ in genes]
        delta_combo = np.zeros(0, dtype=np.float32)
        delta_additive_norm_sq_sum = 0.0
        norm_combo_sum = 0.0

        for batch_i in range(self.N_CTRL_BATCHES):
            sel = self._rng.choice(
                len(self._ctrl_xs), size=self.CELL_SET_LEN, replace=False
            )
            basal = self._ctrl_xs[sel]
            gems = self._ctrl_gem[sel]
            basal_t = torch.tensor(basal, dtype=torch.float32)
            batch_idx_t = torch.tensor(
                [self._batch_idx_for_gem(g) for g in gems], dtype=torch.long
            )

            ctrl_pred = self._forward(self._nontarget_oh, basal_t, batch_idx_t).mean(axis=0)

            # Single-gene predictions
            single_preds = []
            for oh in ohs:
                p = self._forward(oh, basal_t, batch_idx_t).mean(axis=0)
                single_preds.append(p)

            # Combo (n-hot)
            oh_combo = ohs[0].clone()
            for o in ohs[1:]:
                oh_combo = oh_combo + o
            combo_pred = self._forward(oh_combo, basal_t, batch_idx_t).mean(axis=0)

            delta_combo_i = combo_pred - ctrl_pred
            delta_singles_i = [s - ctrl_pred for s in single_preds]
            delta_additive_i = np.sum(np.stack(delta_singles_i, axis=0), axis=0)

            # Accumulate
            if len(delta_combo) == 0:
                delta_combo = delta_combo_i.copy()
                delta_singles = [d.copy() for d in delta_singles_i]
            else:
                delta_combo += delta_combo_i
                for k, d in enumerate(delta_singles_i):
                    delta_singles[k] += d

            delta_additive_norm_sq_sum += float(np.linalg.norm(delta_additive_i) ** 2)
            norm_combo_sum += float(np.linalg.norm(delta_combo_i))

        # Mean over batches
        delta_combo /= self.N_CTRL_BATCHES
        delta_singles = [d / self.N_CTRL_BATCHES for d in delta_singles]
        delta_additive = np.sum(np.stack(delta_singles, axis=0), axis=0)

        norm_combo = float(np.linalg.norm(delta_combo))
        norm_additive = float(np.linalg.norm(delta_additive))
        residual = delta_combo - delta_additive
        synergy_magnitude = float(np.linalg.norm(residual) / (norm_additive + 1e-9))

        # Signed component: dot product of residual with the additive direction
        # (positive = same direction as additive = "super-additive";
        #  negative = opposing = "antagonistic").
        if norm_additive > 1e-9:
            cos = float(np.dot(residual, delta_additive) / (norm_additive * (np.linalg.norm(residual) + 1e-9)))
        else:
            cos = 0.0

        # Map magnitude to a [0, 2] sigma proxy.
        sigma = 2.0 * float(np.tanh(synergy_magnitude))

        if synergy_magnitude > 0.3:
            cls = "synergistic" if cos >= 0 else "antagonistic"
        else:
            cls = "additive"

        return dict(
            additive_baseline=norm_additive,
            observed_double_or_higher=norm_combo,
            synergy_magnitude=synergy_magnitude,
            synergy_class=cls,
            sigma_above_additive=sigma,
            # Replogle HepG2 vocab — record provenance in case downstream wants
            # to gate on it. NOT a K562-bias signal (the HepG2 head is hepatocyte
            # cancer line, not K562), so set to 0.5 = unknown / neutral.
            k562_bias_confidence=0.5,
        )

    # ---- public API ------------------------------------------------------
    def predict_for_hit(self, *, hit_row, context, substrate):
        genes = [hit_row.get(f"gene{i+1}") for i in range(4) if hit_row.get(f"gene{i+1}")]
        genes = [g for g in genes if isinstance(g, str) and g]
        if len(genes) < 2:
            return []

        self.n_attempted += 1
        if self.n_attempted % 20 == 0:
            print(
                f"[state_pair_v2] processed {self.n_attempted} pairs "
                f"(skipped {self.n_skipped}, in-vocab {self.n_in_vocab})",
                flush=True,
            )

        out = self._score_pair(genes)
        if out is None:
            self.n_skipped += 1
            return []
        self.n_in_vocab += 1

        # Pad gene1..gene4
        gene_dict = {f"gene{i+1}": (genes[i] if i < len(genes) else None) for i in range(4)}
        return [dict(**gene_dict, **out)]


if __name__ == "__main__":
    main_for_model(StatePairV2Adapter)

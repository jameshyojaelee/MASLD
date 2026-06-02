"""Shared state.tx Tahoe-100M ST head loader + inference helper.

Used by:
  * d1_mechanism/tahoe_v2_runner.py (drug-induced downstream effects)
  * d2_reversal/tahoe_reversal_v2_runner.py (drug signature reversal)

The checkpoint is the Arc Institute ST-SE-Tahoe head, trained on Tahoe-100M
(~100M drug-perturbation cells across 50 cancer cell lines, 379 unique drugs
at 3 concentrations = 1,138 (drug, conc) perturbation tokens). Outputs:

  * preds                      [N, 2058]   embedding-space delta (SE-600M)
  * pert_cell_counts_preds     [N, 2000]   HVG counts (HVG names NOT in
                                           checkpoint; treat indices as
                                           opaque ranking dims).

VOCAB / API DIFFERENCES FROM HepG2 ST head
-----------------------------------------
  HepG2 ST                        Tahoe ST
  --------                        --------
  control = "non-targeting"       control = "[('DMSO_TF', 0.0, 'uM')]"
  pert_map keys = "GENE_SYMBOL"   pert_map keys = "[('DRUG', conc, 'uM')]"
  ships eval_real.h5ad with        no eval anndata; we reuse the
  X_state basal for HepG2 cells    HepG2 ST head's basal X_state since
  (~7000 control cells)            both SE heads share 2058-dim embedding
                                   space (input compatible).
  output_dim 6546 = full gene      output_dim 2058 (emb) + decoder 2000-HVG
  space (gene_names list matches)  HVG NAMES are NOT in checkpoint;
                                   indices are opaque (see Gotcha below).
  cell_type token absent           cell_type token "HepG2/C3A" passed in
                                   batch (50 cell-line vocab).

KEY GOTCHA: HVG names
---------------------
The Tahoe ST gene_decoder outputs 2,000 HVG counts but the HVG selection
was done in scanpy on the source h5ad (saved as adata.var['highly_variable']);
that mask is not in the released checkpoint. So `pert_cell_counts_preds[i]`
indices are opaque — we expose them as `HVG_0`, `HVG_1`, ... but they cannot
be cross-referenced to gene symbols without the source h5ad.

  → For D1 (gene-level downstream prediction): not directly usable.
    Runner uses the existing Tahoe pseudobulk-DE lookup for gene-level
    output and uses this class's `predict_drug_emb_delta()` magnitude as
    an additional drug-level confidence score.
  → For D2 (reversal): runner uses embedding-delta and pseudobulk-DE cosine
    in tandem (state.tx provides drug-perturbation strength multiplier;
    pseudobulk DE provides directional reversal vs gene-space reference).

Phase 2 plan: re-extract HVG names from a Tahoe-100M subset (rerun
sc.pp.highly_variable_genes(n_top_genes=2000) with the same seed) so the
2000-dim HVG output can be mapped to gene symbols.

"""
from __future__ import annotations

import os
import pickle
import re
from pathlib import Path
from typing import Any

import numpy as np

# Force CPU before any cuda-touching import. Cluster GPUs predate sm_80 builds.
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

# Tahoe ST checkpoint subdirs — pick zeroshot (broader, no fine-tuning to
# held-out drugs) as default; users can switch via env var.
TAHOE_ST_ROOT = (
    PROJECT_ROOT / "data/perturbation/checkpoints/state/ST-SE-Tahoe"
)
TAHOE_ZEROSHOT_DIR = TAHOE_ST_ROOT / "zeroshot/state_generalization_zeroshot_X_state"
TAHOE_FEWSHOT_DIR = TAHOE_ST_ROOT / "fewshot/state_generalization_X_state"

# Reuse the HepG2 ST head's eval anndata to source basal X_state embeddings.
# Both heads share the SE-600M 2058-dim input space — HepG2 control cells
# are valid basal X_state vectors for any downstream ST head.
HEPG2_ST_CKPT_DIR = (
    PROJECT_ROOT
    / "data/perturbation/checkpoints/state/st-se-replogle-full/hepg2_0.99"
)
HEPG2_BASAL_H5AD = HEPG2_ST_CKPT_DIR / "eval_best.ckpt" / "adata_real.h5ad"

# Tahoe-100M drug metadata — list of 379 unique drugs.
TAHOE_DRUG_META = (
    PROJECT_ROOT / "data/perturbation/datasets/tahoe-100m/metadata/drug_metadata.parquet"
)

DEFAULT_CELL_TYPE = "HepG2/C3A"  # Tahoe's only liver-derived cell line
DEFAULT_CONC = 0.5  # μM — middle of the 3 Tahoe doses {0.05, 0.5, 5.0}
DEFAULT_PLATE = "plate1"
CONTROL_PERT_KEY = "[('DMSO_TF', 0.0, 'uM')]"


def _patch_llama_validate() -> None:
    """Disable HF strict LlamaConfig validation (idempotent).

    Reused from HepG2STInference — same architecture quirk.
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


def make_pert_key(drug: str, conc: float = DEFAULT_CONC) -> str:
    """Render a (drug, conc) tuple in the exact Tahoe pert_map key format.

    The pert_map stores keys produced by Python's `repr` of a single-element
    list of tuples — `[('DRUG', conc, 'uM')]`. Concentration is a Python
    float (no trailing zero strip), so `0.5` → `0.5`, `5.0` → `5.0`,
    `0.05` → `0.05`. We rebuild exactly to match.
    """
    return f"[('{drug}', {conc}, 'uM')]"


class TahoeSTInference:
    """Wrap StateTransitionPerturbationModel (Tahoe variant) for inference."""

    CELL_SET_LEN = 64

    def __init__(
        self,
        ckpt_dir: str | Path | None = None,
        n_ctrl_cache: int | None = None,
        basal_h5ad: str | Path = HEPG2_BASAL_H5AD,
        cell_type: str = DEFAULT_CELL_TYPE,
        plate: str = DEFAULT_PLATE,
    ):
        import torch
        import anndata as ad

        _patch_llama_validate()

        # Allow env-var override for zeroshot vs fewshot
        if ckpt_dir is None:
            mode = os.environ.get("TAHOE_ST_MODE", "zeroshot").lower()
            ckpt_dir = (
                TAHOE_FEWSHOT_DIR if mode.startswith("few") else TAHOE_ZEROSHOT_DIR
            )

        ckpt_dir = Path(ckpt_dir)
        if not ckpt_dir.exists():
            raise FileNotFoundError(f"Missing Tahoe ST checkpoint dir: {ckpt_dir}")
        self.ckpt_dir = ckpt_dir

        from state.tx.models.state_transition import StateTransitionPerturbationModel

        print(f"[tahoe_st] loading from {ckpt_dir}", flush=True)
        self.model = StateTransitionPerturbationModel.load_from_checkpoint(
            ckpt_dir / "checkpoints" / "final.ckpt",
            map_location="cpu",
            strict=False,
        )
        self.model.eval()
        self._torch = torch

        # ---- pert / batch / cell-type vocabs ---------------------------
        self.pert_map: dict[str, Any] = torch.load(
            ckpt_dir / "pert_onehot_map.pt", weights_only=False
        )
        with open(ckpt_dir / "batch_onehot_map.pkl", "rb") as fh:
            self.batch_map = pickle.load(fh)
        with open(ckpt_dir / "cell_type_onehot_map.pkl", "rb") as fh:
            self.cell_type_map = pickle.load(fh)
        with open(ckpt_dir / "var_dims.pkl", "rb") as fh:
            self.var_dims = pickle.load(fh)

        if CONTROL_PERT_KEY not in self.pert_map:
            raise RuntimeError(
                f"control pert {CONTROL_PERT_KEY!r} not in Tahoe vocab"
            )
        self.dmso_oh = self.pert_map[CONTROL_PERT_KEY].float()

        if cell_type not in self.cell_type_map:
            raise RuntimeError(
                f"cell_type {cell_type!r} not in Tahoe vocab; "
                f"available: {list(self.cell_type_map.keys())[:5]} ..."
            )
        self.cell_type = cell_type
        if plate not in self.batch_map:
            raise RuntimeError(
                f"plate {plate!r} not in Tahoe vocab; "
                f"available: {list(self.batch_map.keys())}"
            )
        self.plate = plate
        self._batch_idx_scalar = self._batch_idx(plate)

        # ---- Drug → in-vocab concentrations ---------------------------
        # Parse pert_map keys into {drug_name -> [conc, ...]} for fast lookup.
        self._drug_to_concs: dict[str, list[float]] = {}
        self._key_pattern = re.compile(
            r"^\[\('(?P<drug>.+)',\s*(?P<conc>[0-9.eE+-]+),\s*'(?P<unit>[^']+)'\)\]$"
        )
        for key in self.pert_map.keys():
            m = self._key_pattern.match(key)
            if not m:
                continue
            d = m.group("drug")
            c = float(m.group("conc"))
            self._drug_to_concs.setdefault(d, []).append(c)
        # Drop DMSO from the drug list (treated separately as control).
        self._drug_to_concs.pop("DMSO_TF", None)

        # ---- Basal X_state (reuse HepG2 SE basal embeddings) ----------
        # Both heads share 2058-dim SE input space; HepG2 non-targeting
        # control cells are valid basal X_state for any 2058-dim ST head.
        a = ad.read_h5ad(basal_h5ad)
        ctrl_mask = (a.obs["gene"] == "non-targeting").values
        self._ctrl_xs = a.obsm["X_state"][ctrl_mask].astype(np.float32)
        if n_ctrl_cache is not None and n_ctrl_cache < len(self._ctrl_xs):
            rng_cache = np.random.RandomState(0)
            sel = rng_cache.choice(
                len(self._ctrl_xs), size=n_ctrl_cache, replace=False
            )
            self._ctrl_xs = self._ctrl_xs[sel]
        self._rng = np.random.RandomState(42)

        # ---- Output naming ---------------------------------------------
        # HVG names are NOT in the checkpoint; emit placeholders.
        n_hvg = int(self.var_dims.get("hvg_dim", 2000))
        self.hvg_names: list[str] = [f"HVG_{i}" for i in range(n_hvg)]
        n_emb = int(self.var_dims.get("output_dim", 2058))
        self.emb_dim_names: list[str] = [f"STATE_DIM_{i}" for i in range(n_emb)]

        print(
            f"[tahoe_st] loaded: drugs={len(self._drug_to_concs)} "
            f"perts={len(self.pert_map)} cell_types={len(self.cell_type_map)} "
            f"plates={len(self.batch_map)} basal={self._ctrl_xs.shape} "
            f"hvg_dim={n_hvg} emb_dim={n_emb} cell_type={self.cell_type}",
            flush=True,
        )

    # ------------------------------------------------------------------
    def _batch_idx(self, plate: str) -> int:
        import torch
        v = self.batch_map.get(plate)
        if v is None:
            return 0
        if torch.is_tensor(v) and v.ndim == 1:
            return int(torch.argmax(v).item())
        if isinstance(v, (int, np.integer)):
            return int(v)
        return 0

    # ---- Vocab helpers ------------------------------------------------
    def available_drugs(self) -> set[str]:
        """Set of unique drug names (case-preserved) in the Tahoe vocab."""
        return set(self._drug_to_concs.keys())

    def in_vocab(self, drug: str, conc: float | None = None) -> bool:
        """True iff `drug` (case-sensitive) is in Tahoe vocab at this conc.

        If `conc is None`, returns True iff ANY concentration is available.
        """
        concs = self._drug_to_concs.get(drug)
        if concs is None:
            return False
        if conc is None:
            return True
        return any(abs(c - conc) < 1e-9 for c in concs)

    def best_conc(self, drug: str, prefer: float = DEFAULT_CONC) -> float | None:
        """Pick the available conc closest to `prefer` (default 0.5 uM)."""
        concs = self._drug_to_concs.get(drug)
        if not concs:
            return None
        return min(concs, key=lambda c: abs(c - prefer))

    def case_insensitive_lookup(self, drug: str) -> str | None:
        """Resolve a user-supplied drug name to the case-preserved vocab key."""
        if drug in self._drug_to_concs:
            return drug
        low = drug.lower()
        for k in self._drug_to_concs.keys():
            if k.lower() == low:
                return k
        return None

    # ---- Sampling / forward ------------------------------------------
    def sample_basal(self, n: int | None = None, seed_offset: int = 0):
        """Return (basal[float32,N,2058], batch_idx[N], gem[str,N])."""
        torch = self._torch
        n = n or self.CELL_SET_LEN
        if seed_offset:
            rng = np.random.RandomState(42 + seed_offset)
        else:
            rng = self._rng
        sel = rng.choice(len(self._ctrl_xs), size=n, replace=False)
        basal = self._ctrl_xs[sel]
        basal_t = torch.tensor(basal, dtype=torch.float32)
        batch_idx_t = torch.tensor(
            [self._batch_idx_scalar] * n, dtype=torch.long
        )
        return basal_t, batch_idx_t

    def forward(self, pert_oh, basal_t, batch_idx_t) -> dict:
        """Return raw model output dict {preds, pert_cell_counts_preds, ...}.

        preds                  → [N, 2058]  embedding-space residual
        pert_cell_counts_preds → [N, 2000]  HVG counts (post-ReLU)
        """
        torch = self._torch
        n = basal_t.shape[0]
        with torch.no_grad():
            batch = {
                "ctrl_cell_emb": basal_t,
                "pert_emb": pert_oh.float().unsqueeze(0).repeat(n, 1),
                "pert_name": ["pert"] * n,
                "batch": batch_idx_t,
                "cell_type": [self.cell_type] * n,
            }
            out = self.model.predict_step(batch, batch_idx=0, padded=False)
        result: dict = {}
        for k in ("preds", "pert_cell_counts_preds"):
            v = out.get(k, None)
            if v is None:
                result[k] = None
            else:
                result[k] = v.detach().cpu().numpy()
        return result

    # ---- Drug delta predictions --------------------------------------
    def predict_drug_delta(
        self,
        drug: str,
        conc: float | None = None,
        n_batches: int = 2,
    ) -> dict:
        """Predict mean (treated − DMSO) responses for a drug.

        Args
        ----
        drug : str        — case-sensitive vocab name (use ``case_insensitive_lookup``).
        conc : float|None — μM dose; None → best available (closest to 0.5 uM).
        n_batches : int   — number of 64-cell basal samples to average over.

        Returns
        -------
        dict with keys:
          drug, conc, pert_key, n_cells,
          emb_delta[2058]              embedding-space delta (residual)
          emb_ctrl[2058]               DMSO embedding (post model)
          hvg_delta[2000]              HVG counts delta (or None if no decoder)
          hvg_ctrl[2000]               DMSO HVG counts (or None)
          emb_l2                       ||emb_delta||
          hvg_l2                       ||hvg_delta|| (or None)

        Raises KeyError if drug not in vocab.
        """
        torch = self._torch
        canonical = self.case_insensitive_lookup(drug)
        if canonical is None:
            raise KeyError(f"drug not in Tahoe vocab: {drug!r}")
        if conc is None:
            conc = self.best_conc(canonical)
        if conc is None:
            raise KeyError(f"no concentration for {canonical!r}")
        if not self.in_vocab(canonical, conc):
            # snap to nearest available
            concs = self._drug_to_concs[canonical]
            conc = min(concs, key=lambda c: abs(c - conc))

        pert_key = make_pert_key(canonical, conc)
        if pert_key not in self.pert_map:
            raise KeyError(f"pert_key not in Tahoe vocab: {pert_key!r}")
        drug_oh = self.pert_map[pert_key].float()

        emb_delta_sum = None
        emb_ctrl_sum = None
        hvg_delta_sum = None
        hvg_ctrl_sum = None
        for b in range(n_batches):
            basal_t, batch_idx_t = self.sample_basal(seed_offset=b)
            out_ctrl = self.forward(self.dmso_oh, basal_t, batch_idx_t)
            out_drug = self.forward(drug_oh, basal_t, batch_idx_t)

            ctrl_emb = out_ctrl["preds"].mean(axis=0)
            drug_emb = out_drug["preds"].mean(axis=0)
            d_emb = drug_emb - ctrl_emb
            if emb_delta_sum is None:
                emb_delta_sum = d_emb.copy()
                emb_ctrl_sum = ctrl_emb.copy()
            else:
                emb_delta_sum += d_emb
                emb_ctrl_sum += ctrl_emb

            if out_ctrl["pert_cell_counts_preds"] is not None:
                ctrl_h = out_ctrl["pert_cell_counts_preds"].mean(axis=0)
                drug_h = out_drug["pert_cell_counts_preds"].mean(axis=0)
                d_hvg = drug_h - ctrl_h
                if hvg_delta_sum is None:
                    hvg_delta_sum = d_hvg.copy()
                    hvg_ctrl_sum = ctrl_h.copy()
                else:
                    hvg_delta_sum += d_hvg
                    hvg_ctrl_sum += ctrl_h

        emb_delta = emb_delta_sum / n_batches
        emb_ctrl = emb_ctrl_sum / n_batches
        hvg_delta = (hvg_delta_sum / n_batches) if hvg_delta_sum is not None else None
        hvg_ctrl = (hvg_ctrl_sum / n_batches) if hvg_ctrl_sum is not None else None

        return dict(
            drug=canonical,
            conc=float(conc),
            pert_key=pert_key,
            n_cells=int(self.CELL_SET_LEN * n_batches),
            emb_delta=emb_delta,
            emb_ctrl=emb_ctrl,
            hvg_delta=hvg_delta,
            hvg_ctrl=hvg_ctrl,
            emb_l2=float(np.linalg.norm(emb_delta)),
            hvg_l2=(float(np.linalg.norm(hvg_delta)) if hvg_delta is not None else None),
        )

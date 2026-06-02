"""Tahoe v2 runner for D2 reversal — REAL state.tx Tahoe ST inference.

Augments ``tahoe_reversal_runner.py`` (v1) with the Arc Institute
ST-SE-Tahoe head. v1 uses pseudobulk-DE cosine reversal vs a reference
LFC vector (mostly HepG2/C3A pseudobulk). v2 layers ST-head drug-strength
signal on top, so reversal predictions become:

    reversal_score_v2 = - cosine(drug_de_lfc, ref_lfc) * st_strength_mult

with ``st_strength_mult`` ∈ [1.0, 1.5] derived from the Tahoe ST head's
embedding-space delta L2 norm for the drug-of-interest at HepG2/C3A
basal. Drugs NOT in the ST vocab → st_strength_mult = 1.0 → identical to
v1 (graceful fallback).

Why hybrid?
-----------
Same gotcha as D1: Tahoe ST gene_decoder outputs 2,000 HVG dims but the
HVG SYMBOLS are not in the released checkpoint (see
``shared/state_tx_tahoe.py`` docstring "KEY GOTCHA"). So we can't run
cosine of the ST head's gene-space output against a gene-symbol-indexed
reference LFC vector. We CAN use the ST head's embedding-space delta
magnitude (basis-free L2 norm in SE-600M's 2058-dim space) as a drug-level
perturbation-strength signal — that captures "does this drug actually
shift the transcriptome at HepG2/C3A?" which is exactly what we want as
a multiplier on the cosine direction.

For each hit gene:
  1. (Identical to v1) Lookup drugs targeting `g`. Intersect with Tahoe
     pseudobulk-DE vocab and require ≥1 drug to have a cached DE profile.
  2. (Identical to v1) Aggregate per-drug DE → drug-level LFC vector.
     Cosine vs reference LFC → reversal_score_v1.
  3. (NEW) For each drug-in-ST-vocab, run Tahoe ST head forward and
     compute emb_l2. Mean across in-vocab drugs → mean_st.
  4. (NEW) st_strength_mult = 1.0 + 0.5 * tanh(mean_st / 1.0).
  5. reversal_score = reversal_score_v1 * st_strength_mult.

Output schema unchanged; ``st_emb_l2_mean`` + ``st_drugs_covered`` added
as auxiliary fields (best-effort, schema-tolerant).

Set ``TAHOE_ST_DRY=1`` to disable ST head (fall back to v1).
"""
from __future__ import annotations

import os
import sys
import time
import warnings
from pathlib import Path

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

THIS_DIR = Path(__file__).resolve().parent
SHARED = THIS_DIR.parents[1] / "shared"
sys.path.insert(0, str(SHARED))
sys.path.insert(0, str(THIS_DIR))

from tahoe_reversal_runner import TahoeReversalAdapter  # noqa: E402
from _runner_template import main_for_model  # noqa: E402

from state_tx_tahoe import TahoeSTInference  # noqa: E402


class TahoeReversalV2Adapter(TahoeReversalAdapter):
    """Hybrid: ST head strength × pseudobulk-DE cosine reversal."""

    model_name = "tahoe_reversal_v2"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # super().__init__ already initialises most v1 trackers + DE caches —
        # but NOT `_n_predicted` (only set by v1 base `_n_attempted` /
        # `_n_skipped_*`). Add it plus v2-specific ST trackers so neither
        # dry-run nor real-run AttributeErrors.
        self._n_predicted = 0
        self._st_dry = False
        self._st_cache: dict[str, float] = {}
        self.infer = None
        self._st_n_predicted = 0
        self._st_n_oov = 0
        self._st_n_failed = 0

    def load_checkpoint(self) -> None:
        super().load_checkpoint()
        self._st_dry = os.environ.get("TAHOE_ST_DRY", "0") in ("1", "true", "True")
        self._st_cache: dict[str, float] = {}
        if self._st_dry or self.dry_run:
            print(
                "[tahoe_reversal_v2] TAHOE_ST_DRY=1 (or --dry-run); ST disabled, "
                "pseudobulk-DE only.",
                flush=True,
            )
            self.infer = None
        else:
            print("[tahoe_reversal_v2] loading Tahoe ST head (~30-60s)...", flush=True)
            t = time.time()
            self.infer = TahoeSTInference()
            print(f"[tahoe_reversal_v2] ST head loaded in {time.time()-t:.1f}s", flush=True)
        self._st_n_predicted = 0
        self._st_n_oov = 0
        self._st_n_failed = 0

    def _drug_st_emb_l2(self, drug: str) -> float | None:
        """Mean Tahoe ST emb_l2 for `drug` (cached). None → not in ST vocab."""
        if self.infer is None:
            return None
        canonical = self.infer.case_insensitive_lookup(drug)
        if canonical is None:
            self._st_n_oov += 1
            return None
        if canonical in self._st_cache:
            return self._st_cache[canonical]
        try:
            result = self.infer.predict_drug_delta(canonical, n_batches=1)
        except Exception as e:
            print(f"[tahoe_reversal_v2] ST predict failed {drug!r}: {e}", flush=True)
            self._st_n_failed += 1
            return None
        self._st_n_predicted += 1
        self._st_cache[canonical] = float(result["emb_l2"])
        return self._st_cache[canonical]

    def predict_for_hit(self, *, hit_row, context, substrate):
        gene = hit_row.get("gene")
        if not isinstance(gene, str) or not gene:
            return []

        self._n_attempted += 1
        if self._n_attempted % 50 == 0:
            print(
                f"[tahoe_reversal_v2] processed {self._n_attempted} "
                f"(no_drug={self._n_skipped_no_drug}, no_de={self._n_skipped_no_de}, "
                f"no_ov={self._n_skipped_empty_overlap}, "
                f"st_pred={self._st_n_predicted}, st_oov={self._st_n_oov}, "
                f"st_fail={self._st_n_failed})",
                flush=True,
            )

        drugs = self._gene_to_drugs.get(gene, [])
        if not drugs:
            self._n_skipped_no_drug += 1
            return []
        drugs_tahoe = [d for d in drugs if d in self._tahoe_drug_index]
        if not drugs_tahoe:
            self._n_skipped_no_drug += 1
            return []
        drugs_with_de = [d for d in drugs_tahoe if d in self._drug_de]
        if not drugs_with_de:
            self._n_skipped_no_de += 1
            return []

        # ---- Reversal score (identical to v1) ----------------------------
        frames = [self._drug_de[d] for d in drugs_with_de]
        combo = pd.concat(frames, axis=1).mean(axis=1)
        common = combo.index.intersection(self._ref_series.index)
        if len(common) < 50:
            self._n_skipped_empty_overlap += 1
            return []
        drug_vec = combo.loc[common].astype(np.float32).values
        ref_vec = self._ref_series.loc[common].astype(np.float32).values
        cos = self._cosine(drug_vec, ref_vec)
        reversal_v1 = float(-cos)

        # ---- ST strength multiplier (NEW) --------------------------------
        st_l2s = []
        for d in drugs_with_de:
            v = self._drug_st_emb_l2(d)
            if v is not None:
                st_l2s.append(v)
        if st_l2s:
            mean_st = float(np.mean(st_l2s))
            st_mult = 1.0 + 0.5 * float(np.tanh(mean_st / 1.0))
        else:
            mean_st = 0.0
            st_mult = 1.0

        reversal_score = reversal_v1 * st_mult

        self._n_predicted += 1
        return [dict(
            gene=str(gene),
            reversal_score=reversal_score,
            reversal_rank=0,
            stage_specific="—",
            cell_type="hepatocyte_progressor",
            reference_signature=self.reference_signature,
            st_emb_l2_mean=mean_st,
            st_drugs_covered=len(st_l2s),
            st_mult=float(st_mult),
            reversal_score_v1=reversal_v1,
        )]


if __name__ == "__main__":
    main_for_model(TahoeReversalV2Adapter)

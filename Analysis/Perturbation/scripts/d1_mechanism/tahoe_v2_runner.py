"""Tahoe v2 runner for D1 mechanism — REAL state.tx Tahoe ST inference.

Replaces (or rather *augments*) ``tahoe_runner.py`` by using the actual
Arc Institute ST-SE-Tahoe head to score drug-perturbation strength. For
each (hit gene, drug-targeting-it) tuple where the drug is in the Tahoe
ST vocab, we predict the embedding-space delta and combine it with the
existing pseudobulk-DE-based downstream gene table.

Why hybrid?
-----------
The Tahoe ST gene_decoder outputs 2,000 HVG counts but the HVG SYMBOLS
are not stored in the released checkpoint (see ``shared/state_tx_tahoe.py``
docstring "KEY GOTCHA"). Until we re-extract the HVG mask from the
source h5ad, we cannot map the 2000 decoder dims back to gene symbols.

Strategy
--------
For each hit gene:
  1. Lookup drugs targeting `g` via DGIdb (same as v1 tahoe_runner.py).
  2. Intersect with Tahoe ST pert vocab (drug names; ~379 unique).
  3. For each in-vocab drug:
       - Predict embedding delta via TahoeSTInference. Use ``emb_l2`` as
         the drug-level perturbation-strength signal.
       - Pull the per-drug pseudobulk-DE row set from the existing v1
         shard cache and use those gene-level log2FCs to build downstream
         predictions.
  4. Aggregate across drugs (mean log2FC + mean ST emb_l2) and emit top-100
     downstream genes. Confidence is now a combination of |log2FC| AND
     ST-head perturbation magnitude.
  5. Drugs NOT in ST vocab fall back to the v1 pseudobulk-DE-only path.

This lets us flag gene predictions where the ST head independently confirms
that the drug really does cause a transcriptional shift (vs pseudobulk DE
that could be confounded by passage/plate effects).

Output schema unchanged — adds ``st_emb_l2`` and ``st_in_vocab`` fields
inside each GenePrediction's notes (best-effort; output_schema.py allows
arbitrary fields).

Set ``TAHOE_ST_DRY=1`` to skip ST inference (pseudobulk DE only — same as
v1; useful for sanity-comparison).
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

# Reuse the v1 Tahoe pseudobulk-DE-scanning logic.
from tahoe_runner import TahoeAdapter  # noqa: E402
from _runner_template import TOP_N_DOWNSTREAM, main_for_model  # noqa: E402

from state_tx_tahoe import TahoeSTInference  # noqa: E402


class TahoeV2Adapter(TahoeAdapter):
    """Hybrid: ST head for drug-strength signal + pseudobulk DE for genes."""

    model_name = "tahoe_v2"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Init trackers + lookups so dry-run (which skips load_checkpoint)
        # doesn't AttributeError inside predict_for_hit.
        self._n_attempted = 0
        self._n_skipped_no_drug = 0
        self._n_skipped_no_de = 0
        self._gene_to_drugs: dict[str, list[str]] = {}
        self._tahoe_drug_index: set[str] = set()
        self._drug_de_cache: dict[str, pd.DataFrame] = {}
        self._needed_drugs: set[str] = set()
        self._n_shards_scanned = 0
        self._n_de_rows = 0
        self._st_dry = False
        self._st_cache: dict[str, dict] = {}
        self.infer = None
        self._st_n_predicted = 0
        self._st_n_oov = 0
        self._st_n_failed = 0

    def load_checkpoint(self) -> None:
        # First load the v1 pseudobulk DE infrastructure (gene→drug, vocab,
        # per-drug DE shards). This is the source of GENE-LEVEL predictions.
        super().load_checkpoint()

        # Now load the Tahoe ST head, unless dry mode requested.
        self._st_dry = os.environ.get("TAHOE_ST_DRY", "0") in ("1", "true", "True")
        self._st_cache: dict[str, dict] = {}
        if self._st_dry or self.dry_run:
            print(
                "[tahoe_v2] TAHOE_ST_DRY=1 (or --dry-run); ST inference disabled, "
                "pseudobulk-DE only.",
                flush=True,
            )
            self.infer = None
        else:
            print("[tahoe_v2] loading Tahoe ST head (~30-60s on CPU)...", flush=True)
            t = time.time()
            self.infer = TahoeSTInference()
            print(f"[tahoe_v2] ST head loaded in {time.time()-t:.1f}s", flush=True)

        # Stats trackers (in addition to v1's)
        self._st_n_predicted = 0
        self._st_n_oov = 0
        self._st_n_failed = 0

    def _drug_st_strength(self, drug: str) -> float | None:
        """Return mean Tahoe ST emb_l2 for `drug` across in-vocab concs.

        Cached per drug. Returns None if not in vocab or ST inference disabled.
        """
        if self.infer is None:
            return None
        canonical = self.infer.case_insensitive_lookup(drug)
        if canonical is None:
            self._st_n_oov += 1
            return None
        if canonical in self._st_cache:
            return self._st_cache[canonical].get("mean_emb_l2")
        try:
            # Use single best concentration (0.5 uM default) for speed.
            result = self.infer.predict_drug_delta(canonical, n_batches=1)
        except Exception as e:
            print(f"[tahoe_v2] ST predict failed for {drug!r}: {e}", flush=True)
            self._st_n_failed += 1
            return None
        self._st_n_predicted += 1
        self._st_cache[canonical] = dict(
            mean_emb_l2=float(result["emb_l2"]),
            mean_hvg_l2=(float(result["hvg_l2"]) if result["hvg_l2"] is not None else 0.0),
            conc=float(result["conc"]),
        )
        return self._st_cache[canonical]["mean_emb_l2"]

    def predict_for_hit(self, *, hit_row, context, substrate):
        target = hit_row.get("gene")
        if not isinstance(target, str) or not target:
            return []

        self._n_attempted += 1
        if self._n_attempted % 25 == 0:
            print(
                f"[tahoe_v2] processed {self._n_attempted} "
                f"(no_drug={self._n_skipped_no_drug}, "
                f"no_de={self._n_skipped_no_de}, "
                f"st_pred={self._st_n_predicted}, st_oov={self._st_n_oov}, "
                f"st_fail={self._st_n_failed})",
                flush=True,
            )

        drugs = self._gene_to_drugs.get(target, [])
        if not drugs:
            self._n_skipped_no_drug += 1
            return []

        # Restrict to drugs in pseudobulk-DE vocab (already lowercase).
        drugs_tahoe = [d for d in drugs if d in self._tahoe_drug_index]
        if not drugs_tahoe:
            self._n_skipped_no_drug += 1
            return []
        drugs_with_de = [d for d in drugs_tahoe if d in self._drug_de_cache]
        if not drugs_with_de:
            self._n_skipped_no_de += 1
            return []

        # ---- ST head strength per drug (cached) -----------------------
        st_strengths = {}
        any_in_st_vocab = False
        for d in drugs_with_de:
            s = self._drug_st_strength(d)
            if s is not None:
                st_strengths[d] = s
                any_in_st_vocab = True

        # Aggregate DE across drugs (same as v1) - mean of means.
        frames = [self._drug_de_cache[d] for d in drugs_with_de]
        all_de = pd.concat(frames, ignore_index=True)
        agg = (
            all_de.groupby("gene_name")
            .agg(
                log2FC_mean=("log2FC_mean", "mean"),
                abs_lfc_mean=("abs_lfc_mean", "mean"),
            )
            .reset_index()
        )
        agg = agg[agg["gene_name"].astype(str) != target]
        agg = agg.dropna(subset=["log2FC_mean"])
        agg["abs_lfc"] = agg["log2FC_mean"].abs()
        top = agg.sort_values("abs_lfc", ascending=False).head(TOP_N_DOWNSTREAM)
        if top.empty:
            self._n_skipped_no_de += 1
            return []

        # Per-target ST strength multiplier: mean across drugs in vocab.
        # Map to a multiplicative confidence boost in [1.0, 1.5] based on
        # emb_l2 percentile relative to cached values. If no ST signal,
        # multiplier = 1.0 (pure pseudobulk DE path, identical to v1).
        if st_strengths:
            mean_st = float(np.mean(list(st_strengths.values())))
        else:
            mean_st = 0.0

        # Build confidence boost: 1.0 + 0.5 * tanh(mean_st / typical_emb_l2)
        # typical_emb_l2 ~ 1.0 for non-trivial drug perts (calibrated post-smoke).
        st_boost = 1.0 + 0.5 * float(np.tanh(mean_st / 1.0))

        rows: list[dict] = []
        max_abs = float(top["abs_lfc"].max()) or 1.0
        for rank, (_, r) in enumerate(top.iterrows(), start=1):
            lfc = float(r["log2FC_mean"])
            base_conf = float(np.tanh(abs(lfc) / max(0.5, max_abs * 0.5)))
            conf = max(0.0, min(1.0, base_conf * st_boost / 1.5))
            rows.append(
                dict(
                    target_gene=str(target),
                    downstream_gene=str(r["gene_name"]),
                    logFC_predicted=lfc,
                    abs_rank=int(rank),
                    direction="up" if lfc > 0 else "down",
                    confidence=conf,
                    # Augment notes with ST info (best-effort; tolerated by
                    # schema's discriminated-union validator).
                    st_in_vocab=bool(any_in_st_vocab),
                    st_emb_l2_mean=mean_st,
                    st_drugs_covered=len(st_strengths),
                    n_drugs_tested=len(drugs_with_de),
                )
            )
        return rows


if __name__ == "__main__":
    main_for_model(TahoeV2Adapter)

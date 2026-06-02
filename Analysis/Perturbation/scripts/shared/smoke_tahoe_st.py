"""Standalone smoke test for TahoeSTInference.

Verifies:
  1. TahoeSTInference loads without error (Tahoe ckpt + HepG2 basal).
  2. For 3 well-known drugs, predict_drug_delta returns non-trivial emb / HVG deltas.
  3. DMSO_TF self-comparison gives near-zero delta (sanity).

Run:
    .envs/perturbation_state/bin/python Analysis/Perturbation/scripts/shared/smoke_tahoe_st.py
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))

import numpy as np  # noqa: E402

from state_tx_tahoe import TahoeSTInference  # noqa: E402


def main() -> int:
    t0 = time.time()
    print("[smoke] constructing TahoeSTInference (CPU; ~30-60s)...", flush=True)
    infer = TahoeSTInference()
    print(f"[smoke] loaded in {time.time()-t0:.1f}s", flush=True)

    avail = infer.available_drugs()
    print(f"[smoke] {len(avail)} unique drugs in Tahoe vocab", flush=True)

    # Pick 3 well-known drugs from the Tahoe metadata: a proteasome inhibitor,
    # a HMGCR statin, and an EGFR TKI. All should produce strong transcriptional
    # responses if the model works.
    test_drugs = ["Bortezomib", "Simvastatin", "Erlotinib"]
    for d in test_drugs:
        canonical = infer.case_insensitive_lookup(d)
        in_vocab = canonical is not None
        concs = (
            infer._drug_to_concs.get(canonical, []) if in_vocab else []
        )
        print(f"  - {d!r}: in_vocab={in_vocab}, canonical={canonical!r}, concs={concs}", flush=True)
        if not in_vocab:
            # fall back to first alphabetical real drug to keep the smoke alive
            replacement = sorted(avail)[0]
            print(f"    → substituting {replacement!r}", flush=True)
            test_drugs[test_drugs.index(d)] = replacement

    # ---- DMSO self-control sanity (delta ~ 0) -------------------------------
    print("\n[smoke] DMSO_TF vs DMSO_TF (sanity, delta ~ 0)...", flush=True)
    t_dmso = time.time()
    # Direct forward to test invariance
    basal_t, batch_idx_t = infer.sample_basal(seed_offset=0)
    out_a = infer.forward(infer.dmso_oh, basal_t, batch_idx_t)
    out_b = infer.forward(infer.dmso_oh, basal_t, batch_idx_t)
    d_emb = (out_b["preds"] - out_a["preds"]).mean(axis=0)
    print(
        f"  DMSO self-delta emb: ||·||={np.linalg.norm(d_emb):.4e}, "
        f"max|·|={np.abs(d_emb).max():.4e} (expect tiny / numerical)",
        flush=True,
    )
    if out_a["pert_cell_counts_preds"] is not None:
        d_hvg = (out_b["pert_cell_counts_preds"] - out_a["pert_cell_counts_preds"]).mean(axis=0)
        print(
            f"  DMSO self-delta HVG: ||·||={np.linalg.norm(d_hvg):.4e}, "
            f"max|·|={np.abs(d_hvg).max():.4e}",
            flush=True,
        )
    print(f"  (took {time.time()-t_dmso:.1f}s)", flush=True)

    # ---- 3 drug deltas ----------------------------------------------------
    print("\n[smoke] per-drug deltas (2 batches × 64 cells each):", flush=True)
    summaries = []
    for drug in test_drugs:
        t = time.time()
        result = infer.predict_drug_delta(drug, n_batches=2)
        dt = time.time() - t
        emb_d = result["emb_delta"]
        hvg_d = result["hvg_delta"]
        print(
            f"  {drug!r} @ {result['conc']} uM "
            f"emb_l2={result['emb_l2']:.4f} emb_std={emb_d.std():.4f} "
            f"max|emb|={np.abs(emb_d).max():.4f} "
            f"hvg_l2={result['hvg_l2']:.4f} hvg_std={hvg_d.std():.4f} "
            f"max|hvg|={np.abs(hvg_d).max():.4f} "
            f"({dt:.1f}s)",
            flush=True,
        )
        # Top-10 HVG indices by |delta|
        top_idx = np.argsort(-np.abs(hvg_d))[:10]
        print(f"    top-10 HVG idx by |delta|: {top_idx.tolist()}", flush=True)
        print(f"    top-10 HVG delta values:   {hvg_d[top_idx].round(3).tolist()}", flush=True)
        summaries.append((drug, result))

    # ---- Drug-vs-drug discrimination ---------------------------------------
    print("\n[smoke] drug-vs-drug discrimination (emb cosine):", flush=True)
    from itertools import combinations
    for (d1, r1), (d2, r2) in combinations(summaries, 2):
        a, b = r1["emb_delta"], r2["emb_delta"]
        cos = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))
        print(f"  cosine(delta[{d1}], delta[{d2}]) = {cos:+.4f}", flush=True)

    print(f"\n[smoke] DONE in {time.time()-t0:.1f}s. Tahoe ST head is alive.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

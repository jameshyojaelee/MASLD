#!/usr/bin/env python
"""GEARS Phase 2 GI_predict — load Saunders fine-tuned checkpoint and score
D3 invocab tier1/tier2 pairs.

Output schema matches Phase 1 gears_zero_shot output (D3 SynergyPrediction):
  gene1, gene2, additive_baseline, observed_double_or_higher,
  synergy_magnitude, synergy_class, sigma_above_additive, k562_bias_confidence.

For Phase 2:
  - k562_bias_confidence = 0.2 (Saunders mouse hep -> lower K562 bias than 1.0 in Phase 1)
  - skipped pairs (≥1 gene not in perturbable vocab) are emitted as
    `predictions = []` rows for that hit (matching Phase 1 conventions).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint-dir", required=True,
                    help="Path to Saunders fine-tune output dir "
                         "(contains best.ckpt + dataset subdir from gears_finetune.py)")
    ap.add_argument("--dataset-name", default="saunders_hep",
                    help="PertData dataset_name used in fine-tune (e.g. saunders_hep / saunders_hep_smoke)")
    ap.add_argument("--hits", required=True,
                    help="Hits CSV with gene1, gene2 columns")
    ap.add_argument("--tier-label", required=True,
                    help="Tier label for output filename (e.g. 'tier1_pairs', 'tier2_pairs')")
    ap.add_argument("--out-json", required=True,
                    help="Output path for prediction JSON")
    ap.add_argument("--k562-bias", type=float, default=0.2,
                    help="k562_bias_confidence to stamp on each prediction (default: 0.2 for Saunders-tuned)")
    args = ap.parse_args()

    t0 = time.time()
    print(f"=== GEARS Phase 2 GI_predict ===")
    print(f"args: {vars(args)}")

    # ---- Load model ----
    print("\nImporting GEARS...")
    from gears import GEARS, PertData
    import torch

    ckpt_dir = Path(args.checkpoint_dir)
    pert_data = PertData(data_path=str(ckpt_dir))
    pert_data.load(data_name=args.dataset_name)
    try:
        pert_data.prepare_split(split="simulation", seed=42)
        pert_data.get_dataloader(batch_size=32, test_batch_size=64)
    except Exception as e:
        warnings.warn(f"prepare_split/get_dataloader: {e} (probably already split)")

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    model = GEARS(pert_data, device=device)
    model.model_initialize(hidden_size=64)
    best_ckpt = ckpt_dir / "best.ckpt"
    model.load_pretrained(str(best_ckpt))
    print(f"Loaded checkpoint from {best_ckpt}")

    # ---- Perturbable universe ----
    # GEARS' true predict() universe is essential_all_data_pert_genes.pkl (~9,975 genes)
    # — the GO-graph extrapolation set. Prefer this over pert_data.pert_names which
    # may be narrower (e.g. gene_names ~5K).
    pert_universe = set()
    import pickle
    ess_pkl = Path("Analysis/Perturbation/results/finetuned_checkpoints/essential_all_data_pert_genes.pkl")
    if ess_pkl.exists():
        with open(ess_pkl, "rb") as fh:
            ess = pickle.load(fh)
        ess_set = set(ess.tolist() if hasattr(ess, "tolist") else ess)
        ess_set.discard("ctrl")
        pert_universe = ess_set
        print(f"[essential_all_data_pert_genes.pkl] perturbable universe: {len(pert_universe)} genes")
    else:
        for attr in ("pert_names", "perturbable_genes", "gene_names"):
            if hasattr(pert_data, attr):
                v = getattr(pert_data, attr)
                if v is not None and len(v) > 0:
                    pert_universe = set(v)
                    print(f"[{attr}] fallback perturbable universe: {len(pert_universe)} genes")
                    break

    # ---- Load hits ----
    hits = pd.read_csv(args.hits)
    print(f"\nLoaded {len(hits)} hits from {args.hits}")
    hit_genes = sorted(set(hits["gene1"]) | set(hits["gene2"]))
    covered = sorted(set(hit_genes) & pert_universe)
    print(f"  unique genes in hits: {len(hit_genes)}; covered by perturb vocab: {len(covered)}")
    if covered:
        print(f"  covered genes: {covered}")

    # ---- Predict each pair ----
    predictions = []
    n_skipped = 0
    n_attempted = 0
    for _, row in hits.iterrows():
        g1, g2 = str(row["gene1"]), str(row["gene2"])
        if g1 not in pert_universe or g2 not in pert_universe:
            n_skipped += 1
            continue

        n_attempted += 1
        try:
            combo_pred = model.predict([[g1, g2]])
            single1 = model.predict([[g1]])
            single2 = model.predict([[g2]])
        except Exception as e:
            warnings.warn(f"GEARS predict failed for ({g1}, {g2}): {e}")
            continue

        combo_key = "+".join([g1, g2])
        combo_vec = combo_pred.get(combo_key)
        s1 = single1.get(g1)
        s2 = single2.get(g2)
        if combo_vec is None or s1 is None or s2 is None:
            continue

        additive = s1 + s2
        diff = combo_vec - additive
        mag = float(np.linalg.norm(diff))
        additive_mag = float(np.linalg.norm(additive)) + 1e-9
        sigma = mag / additive_mag

        if sigma > 0.3:
            synergy_class = "synergistic"
        elif sigma < -0.3:
            synergy_class = "antagonistic"
        else:
            synergy_class = "additive"

        predictions.append(dict(
            gene1=g1,
            gene2=g2,
            additive_baseline=additive_mag,
            observed_double_or_higher=float(np.linalg.norm(combo_vec)),
            synergy_magnitude=mag,
            synergy_class=synergy_class,
            sigma_above_additive=sigma,
            k562_bias_confidence=args.k562_bias,
        ))

    elapsed = time.time() - t0
    print(f"\nProcessed {len(hits)} hits: attempted={n_attempted}, "
          f"skipped (vocab miss)={n_skipped}, predicted={len(predictions)}")
    print(f"Elapsed: {elapsed:.1f}s")

    # ---- Emit envelope (matches Phase 1 schema) ----
    envelope = {
        "arm": "D3",
        "model": "gears",
        "modality": "fine_tuned",
        "context": f"{args.tier_label}__invocab",
        "timestamp": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "checkpoint_hash": _hash_dir(best_ckpt),
        "n_predictions": len(predictions),
        "predictions": predictions,
        "runtime_seconds": elapsed,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "notes": (
            f"D3 Phase 2 fine-tuned on Saunders 2025 Perturb-Multi mouse hepatocyte. "
            f"perturbable_universe={len(pert_universe)} genes; "
            f"hits_attempted={n_attempted} hits_skipped_vocab={n_skipped}."
        ),
        "model_ejected": False,
    }

    out_path = Path(args.out_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(envelope, fh, indent=2)
    print(f"\nWrote {out_path}")


def _hash_dir(p):
    import hashlib
    if not p.exists():
        return "missing"
    h = hashlib.sha256()
    for f in sorted(p.rglob("*")):
        if f.is_file():
            h.update(f.relative_to(p).as_posix().encode())
            h.update(b":")
            try:
                h.update(str(f.stat().st_size).encode())
            except OSError:
                h.update(b"err")
            h.update(b"\n")
    return h.hexdigest()[:16]


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Create fold-0-free labels and two-fold adapter fits for honest stacking.

For outer held fold f, each meta-training fold g has an adapter fitted without
either f or g. The saved three-fold adapter for f is its evaluation predictor.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
CAM = ROOT / "scripts/analysis/alphagenome_campaign"
SOURCE = ROOT / "GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz"
CONFIG = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model/configs/adapter_r16_last5_symmetric.json"
SEEDS = range(20260917, 20260922)
FOLDS = range(1, 5)


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("GPU compute-node allocation required")
    args.out.mkdir(parents=True, exist_ok=False)
    labels = pd.read_csv(SOURCE, sep="\t")
    if set(labels.heldout_fold.unique()) != set(range(5)):
        raise ValueError("Unexpected source folds")
    labels = labels.loc[labels.heldout_fold.ne(0)].copy()
    if 0 in labels.heldout_fold.unique():
        raise ValueError("Fold 0 survived exclusion")
    if labels.lead_variant_id.duplicated().any():
        raise ValueError("Repeated variant identity")
    filtered = args.out / "development_labels_folds1to4.tsv.gz"
    labels.to_csv(filtered, sep="\t", index=False, compression="gzip")
    del labels
    receipts = []
    for outer in FOLDS:
        for meta in FOLDS:
            if outer == meta:
                continue
            for seed in SEEDS:
                name = f"outer{outer}_meta{meta}_seed{seed}"
                target = args.out / name
                command = [sys.executable, str(CAM/"model_train.py"),
                           "--config", str(CONFIG), "--labels", str(filtered),
                           "--out", str(target), "--steps", "5000", "--microbatch", "4",
                           "--validation-limit", "0", "--held-fold", str(outer),
                           "--validation-fold", str(meta), "--seed", str(seed),
                           "--max-hours", "1.0"]
                tick = time.monotonic()
                result = subprocess.run(command, check=False)
                elapsed = time.monotonic() - tick
                if result.returncode:
                    raise RuntimeError(f"{name} failed with exit {result.returncode}")
                split = json.loads((target/"split.json").read_text())
                expected = sorted(set(FOLDS) - {outer, meta})
                if split["training_folds"] != expected or split["validation_fold"] != meta or split["held_fold"] != outer:
                    raise ValueError(f"{name}: two-fold training split mismatch")
                fit = json.loads((target/"feasibility.json").read_text())
                if fit["completed_steps"] != 5000:
                    raise ValueError(f"{name}: incomplete training")
                pred = pd.read_csv(target/"validation_predictions.tsv", sep="\t", usecols=["validation_fold"])
                if set(pred.validation_fold.unique()) != {meta}:
                    raise ValueError(f"{name}: validation file contains wrong fold")
                receipts.append({"name": name, "outer_fold": outer, "meta_fold": meta,
                                 "seed": seed, "training_folds": expected,
                                 "validation_rows": len(pred), "wall_seconds": elapsed,
                                 "fit_seconds": fit.get("training_wall_seconds_including_batching_and_progress_writes")})
                (args.out/"progress.json").write_text(json.dumps(receipts, indent=2) + "\n")
                print(json.dumps(receipts[-1]), flush=True)
    (args.out/"completion.json").write_text(json.dumps({"status": "complete",
        "fits": len(receipts), "expected_fits": 60, "fold0_labels_used": 0,
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "slurm_job_id": os.environ["SLURM_JOB_ID"]}, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())

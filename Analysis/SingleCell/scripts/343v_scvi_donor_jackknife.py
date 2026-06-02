#!/usr/bin/env python
"""
343v_scvi_donor_jackknife.py  (S2 deliverable #1)

Honest donor-level leave-one-out for the scVI F-stage projection.

Why this script exists
----------------------
The QWK = 0.76 number reported in `memory/MEMORY.md` for F-stage inference
was produced by kNN-5 leave-one-donor-out on a FIXED scVI latent space —
i.e. scVI was trained on the full atlas (all 269 donors), and only the
classifier head was re-fit per fold. The latent representation itself
already saw every held-out donor's cells, so per-donor LOOCV on the
classifier head is not a valid generalization estimate (S2 scANVI labels
audit, P0).

This script does the proper version: for ONE held-out donor it
  1. loads the hepatocyte atlas
  2. removes ALL cells from that donor
  3. trains scVI from scratch on the remaining cells
  4. projects the held-out donor's cells (latent-only inference, no
     parameter updates) using `model.get_latent_representation(query)`
  5. computes donor-level pseudobulk in latent space
  6. trains kNN-5 on the OTHER Andrews donors with documented F-stage,
     predicts the held-out donor, records (donor, true_F, pred_F).

Run 10 such jobs (one per held-out donor), then aggregate via
`343x_aggregate_loocv.R` to get the honest QWK distribution.

Usage
-----
    python 343v_scvi_donor_jackknife.py --hold-out-donor GSM6112210

CLI
---
--hold-out-donor   donor sample id (string, matches adata.obs['sample'])
--out-dir          override default output dir
--max-epochs       scVI max_epochs (default 100; matches 308 setup but capped
                   for wall-time on jackknife)
--n-latent         scVI n_latent (default 20)
--seed             RNG seed (default 42)

Outputs
-------
results_gpu_v2/scvi_validation/donor_jackknife/<donor>/
    scvi_model/                              # the retrained scVI checkpoint
    held_out_latent.npz                      # latent for held-out cells
    donor_pseudobulk.tsv                     # held-out donor pseudobulk latent
    knn_prediction.tsv                       # (sample, F_true, F_pred, P_F0..P_F4)
    summary.json                             # one-line metrics
"""
from __future__ import annotations
import argparse
import gc
import json
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [343v] %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
HEP_H5AD   = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad"
FSTAGE_DOC = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv"
WT_ROOT    = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation")
DEFAULT_OUT = WT_ROOT / "Analysis/SingleCell/results_gpu_v2/scvi_validation/donor_jackknife"


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hold-out-donor", required=True)
    ap.add_argument("--out-dir", default=str(DEFAULT_OUT))
    ap.add_argument("--max-epochs", type=int, default=100)
    ap.add_argument("--n-latent", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    return ap.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir) / args.hold_out_donor
    out_dir.mkdir(parents=True, exist_ok=True)

    import torch
    import anndata as ad
    import scvi
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.metrics import cohen_kappa_score

    log.info("hold-out donor: %s", args.hold_out_donor)
    log.info("output dir:     %s", out_dir)
    log.info("torch CUDA:     %s", torch.cuda.is_available())

    scvi.settings.seed = args.seed
    np.random.seed(args.seed)

    t0 = time.time()
    log.info("loading %s", HEP_H5AD.name)
    adata = ad.read_h5ad(HEP_H5AD)
    log.info("loaded %s in %.1fs", adata.shape, time.time() - t0)

    if "counts" not in adata.layers:
        sys.exit("FATAL: 'counts' layer missing — scVI needs raw counts")

    if args.hold_out_donor not in set(adata.obs["sample"].astype(str)):
        sys.exit(f"FATAL: donor {args.hold_out_donor} not in adata.obs['sample']")

    held_mask = (adata.obs["sample"].astype(str) == args.hold_out_donor).values
    log.info("held-out cells: %d   training cells: %d",
             held_mask.sum(), (~held_mask).sum())

    # Subset to highly variable genes (use the var flag baked in by 308)
    if "highly_variable" in adata.var.columns:
        adata = adata[:, adata.var["highly_variable"].values].copy()
        log.info("subsetted to HVGs: %s", adata.shape)
    else:
        log.warning("no highly_variable flag — using all genes (matches 308 if HVGs already subset upstream)")

    a_train = adata[~held_mask].copy()
    a_held  = adata[held_mask].copy()
    del adata
    gc.collect()

    log.info("setting up SCVI on training subset (%d cells)", a_train.n_obs)
    scvi.model.SCVI.setup_anndata(a_train, layer="counts", batch_key="dataset")
    model = scvi.model.SCVI(a_train, n_latent=args.n_latent, n_layers=2,
                            gene_likelihood="nb")
    log.info("training scVI: max_epochs=%d", args.max_epochs)
    t0 = time.time()
    model.train(max_epochs=args.max_epochs,
                early_stopping=True, early_stopping_patience=10,
                train_size=0.9, batch_size=256)
    log.info("scVI trained in %.1fs", time.time() - t0)
    model.save(str(out_dir / "scvi_model"), overwrite=True)

    log.info("projecting held-out donor cells into latent space")
    scvi.model.SCVI.prepare_query_anndata(a_held, model)
    z_held = model.get_latent_representation(a_held)
    np.savez_compressed(out_dir / "held_out_latent.npz", z=z_held,
                        cell_ids=a_held.obs.index.values.astype(str))

    z_train = model.get_latent_representation(a_train)

    # Donor-level pseudobulk
    cols = [f"z{i}" for i in range(z_train.shape[1])]
    pb_train = pd.DataFrame(z_train, columns=cols)
    pb_train["sample"]  = a_train.obs["sample"].astype(str).values
    pb_train["dataset"] = a_train.obs["dataset"].astype(str).values
    pb_train = pb_train.groupby(["sample", "dataset"], observed=True)[cols].mean().reset_index()

    pb_held = pd.DataFrame(z_held, columns=cols)
    pb_held["sample"]  = a_held.obs["sample"].astype(str).values
    pb_held["dataset"] = a_held.obs["dataset"].astype(str).values
    pb_held = pb_held.groupby(["sample", "dataset"], observed=True)[cols].mean().reset_index()
    pb_held.to_csv(out_dir / "donor_pseudobulk.tsv", sep="\t", index=False)

    fs = pd.read_csv(FSTAGE_DOC, sep="\t")
    fs["F_stage_documented"] = pd.to_numeric(fs["F_stage_documented"], errors="coerce")

    pb_train = pb_train.merge(fs[["sample", "F_stage_documented"]], on="sample", how="left")
    train_lab = pb_train.dropna(subset=["F_stage_documented"]).copy()
    log.info("labelled training donors (F-stage known): %d", len(train_lab))

    X_train = train_lab[cols].to_numpy()
    y_train = train_lab["F_stage_documented"].astype(int).to_numpy()

    # kNN-5 (matches 343c finding: kNN-5 in scVI latent = best per MEMORY.md)
    clf = KNeighborsClassifier(n_neighbors=5, weights="distance")
    clf.fit(X_train, y_train)

    X_held = pb_held[cols].to_numpy()
    pred = int(clf.predict(X_held)[0])
    proba_raw = clf.predict_proba(X_held)[0]
    proba_full = np.zeros(5)
    for j, c in enumerate(clf.classes_):
        proba_full[int(c)] = proba_raw[j]

    true_F = fs.loc[fs["sample"] == args.hold_out_donor, "F_stage_documented"]
    true_F = int(true_F.iloc[0]) if len(true_F) and not np.isnan(true_F.iloc[0]) else None

    row = dict(
        sample=args.hold_out_donor,
        F_true=true_F,
        F_pred=pred,
        abs_err=(abs(true_F - pred) if true_F is not None else None),
        **{f"P_F{k}": float(proba_full[k]) for k in range(5)},
        n_train_donors_labelled=len(train_lab),
        n_held_cells=int(held_mask.sum()),
        scvi_epochs=int(args.max_epochs),
    )
    pd.DataFrame([row]).to_csv(out_dir / "knn_prediction.tsv", sep="\t", index=False)

    with open(out_dir / "summary.json", "w") as f:
        json.dump(row, f, indent=2, default=str)

    log.info("DONE  sample=%s  true=%s  pred=%d  abs_err=%s",
             args.hold_out_donor, true_F, pred,
             abs(true_F - pred) if true_F is not None else "NA")


if __name__ == "__main__":
    main()

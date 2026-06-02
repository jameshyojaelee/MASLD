"""352b_fstage_on_latent.py — train F-stage classifier on a chosen latent.

Generates donor-level F-stage predictions from any chosen integration latent
(obsm key). Used for the three-way method-dependence test (critique #9):

  --latent X_scVI       (existing in v2 atlas; reproduces 04_fstage v2)
  --latent X_harmony    (existing in v2 atlas obsm)
  --latent X_scanorama  (from 352a_scanorama_integration.py)

For each, donor-mean of the latent across hepatocytes is computed; F-stage
classifier (sklearn LogisticRegression multinomial) is trained on documented
F-stage anchors (n=58 Andrews) plus clean-healthy F0 anchors plus cirrhosis
F4 anchors (matches the v2 augmented design); predictions written for all
269 donors.

Outputs (suffix = chosen latent):
  - results_gpu_v2_phase05/ccc/stage_trajectory_v3/donor_fstage_<suffix>.tsv
    Columns: sample, dataset, latent, F_stage_predicted_argmax, P_F0, ..., P_F4
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS_DEFAULT = ROOT / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2_annotated.h5ad"
SCAN_ATLAS = ROOT / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2_scanorama.h5ad"
META_V1 = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
FSTAGE_DOC = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv"
DUBIOUS = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/dubious_healthy_donors.tsv"
OUT_DIR = ROOT / "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3"


def donor_mean_of_latent(adata: ad.AnnData, latent_key: str) -> pd.DataFrame:
    """Compute donor-level mean of an obsm latent. Returns wide DataFrame:
    sample, dataset, n_cells, L0, L1, ..., L{d-1}."""
    if latent_key not in adata.obsm:
        raise KeyError(f"latent {latent_key} not in obsm; keys = {list(adata.obsm)}")
    Z = np.asarray(adata.obsm[latent_key])
    if np.isnan(Z).all():
        raise ValueError(f"latent {latent_key} is all-NaN")
    n_lat = Z.shape[1]
    cols = [f"L{i}" for i in range(n_lat)]
    df = pd.DataFrame(Z, columns=cols)
    df["sample"] = adata.obs["sample"].astype(str).to_numpy()
    df["dataset"] = adata.obs["dataset"].astype(str).to_numpy()
    # Drop rows that are all-NaN (e.g. cells without Scanorama alignment)
    df = df.dropna(subset=cols, how="all")
    log.info(f"  {len(df):,} cells with non-NaN latent")
    pb = df.groupby(["sample", "dataset"], observed=True)[cols].mean().reset_index()
    n_cells = (df.groupby(["sample", "dataset"], observed=True)
                  .size().rename("n_cells_for_latent").reset_index())
    pb = pb.merge(n_cells, on=["sample", "dataset"])
    return pb, cols


def build_anchors(pb: pd.DataFrame, latent_cols: list[str]) -> pd.DataFrame:
    """Build training-anchor table:
      - documented F-stage from FSTAGE_DOC (n=58 Andrews)
      - clean-healthy (disease_stage_coarse=Healthy, not in dubious list) -> F0
      - cirrhosis (disease_stage_coarse=Cirrhosis) -> F4 (only if not documented)
    """
    meta = pd.read_csv(META_V1, sep="\t")
    if "sample" not in meta.columns:
        raise ValueError(f"META_V1 missing 'sample' column")
    doc = pd.read_csv(FSTAGE_DOC, sep="\t")
    doc["F_stage_documented"] = pd.to_numeric(doc["F_stage_documented"],
                                              errors="coerce")
    if not DUBIOUS.exists():
        log.error(f"DUBIOUS file not found: {DUBIOUS}")
        log.error("F0 anchor pool would be contaminated; refusing to proceed")
        sys.exit(2)
    dub = pd.read_csv(DUBIOUS, sep="\t")
    dubious_set = set(dub["sample"].astype(str))
    log.info(f"loaded {len(dubious_set)} dubious-healthy donors from {DUBIOUS.name}")

    keep = ["sample", "dataset", "disease_stage_coarse"]
    df = meta[keep].copy()
    df = df.merge(pb, on=["sample", "dataset"], how="left")
    df = df.merge(doc[["sample", "dataset", "F_stage_documented"]],
                  on=["sample", "dataset"], how="left")
    df["dubious"] = df["sample"].astype(str).isin(dubious_set)

    df["y_train"] = np.nan
    df["origin"] = "none"
    m_doc = df["F_stage_documented"].notna()
    df.loc[m_doc, "y_train"] = df.loc[m_doc, "F_stage_documented"]
    df.loc[m_doc, "origin"] = "documented"
    m_ch = ((df["disease_stage_coarse"] == "Healthy") &
            (~df["dubious"]) & (~m_doc))
    df.loc[m_ch, "y_train"] = 0.0
    df.loc[m_ch, "origin"] = "clean_healthy_anchor"
    m_cirr = ((df["disease_stage_coarse"] == "Cirrhosis") & (~m_doc))
    df.loc[m_cirr, "y_train"] = 4.0
    df.loc[m_cirr, "origin"] = "cirrhosis_anchor"
    return df


def train_and_predict(df: pd.DataFrame, latent_cols: list[str],
                      latent_name: str) -> pd.DataFrame:
    # Training set
    train_mask = (df["y_train"].notna()
                  & df[latent_cols].notna().all(axis=1))
    train = df[train_mask].copy()
    n_doc = (train["origin"] == "documented").sum()
    n_ch = (train["origin"] == "clean_healthy_anchor").sum()
    n_cirr = (train["origin"] == "cirrhosis_anchor").sum()
    log.info(f"[{latent_name}] training anchors: {n_doc} documented + "
             f"{n_ch} clean-healthy + {n_cirr} cirrhosis = {len(train)} total")
    X_train = train[latent_cols].to_numpy()
    y_train = train["y_train"].to_numpy().astype(int)
    scaler = StandardScaler().fit(X_train)
    X_train_s = scaler.transform(X_train)
    clf = LogisticRegression(
        class_weight="balanced",
        max_iter=2000,
        solver="lbfgs",
        C=1.0,
    )
    clf.fit(X_train_s, y_train)

    pred_mask = df[latent_cols].notna().all(axis=1)
    pred = df[pred_mask].copy()
    X_pred = scaler.transform(pred[latent_cols].to_numpy())
    proba = clf.predict_proba(X_pred)
    classes = clf.classes_
    # Pad to F0..F4 in case some class missing in training
    proba_full = np.zeros((proba.shape[0], 5), dtype=np.float32)
    for j, c in enumerate(classes):
        proba_full[:, int(c)] = proba[:, j]
    argmax = np.argmax(proba_full, axis=1)
    out = pred[["sample", "dataset", "disease_stage_coarse", "origin"]].copy()
    out["latent"] = latent_name
    out["F_stage_predicted_argmax"] = argmax
    for k in range(5):
        out[f"P_F{k}"] = proba_full[:, k]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--latent", required=True,
                    choices=["X_scVI", "X_harmony", "X_scanorama"])
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    suffix = args.latent.replace("X_", "").lower()
    out_path = OUT_DIR / f"donor_fstage_{suffix}.tsv"

    # Pick atlas based on latent
    if args.latent == "X_scanorama":
        atlas_path = SCAN_ATLAS
    else:
        atlas_path = ATLAS_DEFAULT
    if not atlas_path.exists():
        log.error(f"atlas not found: {atlas_path}")
        sys.exit(1)

    log.info(f"loading atlas: {atlas_path}")
    adata = ad.read_h5ad(atlas_path)
    log.info(f"  shape: {adata.shape}; obsm keys: {list(adata.obsm)}")

    pb, latent_cols = donor_mean_of_latent(adata, args.latent)
    log.info(f"donor pseudobulk: {len(pb)} (sample, dataset) rows")

    df = build_anchors(pb, latent_cols)
    n_with = df[latent_cols].notna().all(axis=1).sum()
    log.info(f"donors with latent values: {n_with} / {len(df)}")

    out = train_and_predict(df, latent_cols, args.latent)
    out.to_csv(out_path, sep="\t", index=False)
    log.info(f"wrote -> {out_path} ({len(out)} rows)")

    # Quick sanity: distribution of predictions
    pred_dist = out.groupby("F_stage_predicted_argmax").size()
    log.info(f"pred F-stage distribution: {pred_dist.to_dict()}")
    log.info("DONE")


if __name__ == "__main__":
    main()

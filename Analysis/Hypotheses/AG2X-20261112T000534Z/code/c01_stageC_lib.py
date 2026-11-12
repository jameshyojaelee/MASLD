"""AG2X Stage C shared definitions (PRESPEC_AG2X_stageC.md). No network access in this module."""
from __future__ import annotations

import glob
import json
import pathlib

import numpy as np
import pandas as pd

PKG = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z")
C = PKG / "stageC"
RAW = C / "raw"
TAB = C / "tables"
CELLS = C / "inputs/cellsets"
LEDGER = RAW / "call_ledger.jsonl"
SEED = 20261111
STAGEC_CAP = 10_000
SETS = ["set1", "set2", "set3"]
STATES = ["ccl20pos", "otherhep", "stellate"]
SCORED_STATES = ["otherhep", "ccl20pos"]  # scoring order: otherhep first (control b), then ccl20pos
# SC_ATAC protocols observed in the reference hepatocyte condition (Stage B HEP_SCATAC)
PROTOCOLS = ["SHARE-seq-cellranger-nuclear", "scATAC-seq 10X-bowtie2-nuclear", "sciATAC-seq-bowtie2-nuclear"]
USABLE_MODES = {"observed", "imputed_supported"}
REF_HEP = PKG / "stageB/raw/t1_hepatocyte"
REF_STELLATE = PKG / "stageB/raw/t1_stellate"
# fixed thresholds
CTRL_A_MIN_BAL_ACC = 0.90
CTRL_B_MIN_RHO = 0.80
CALL_PCT = 0.99


def cond_name(s: str, state: str) -> str:
    return f"{s}_{state}"


def spent() -> int:
    if not LEDGER.exists():
        return 0
    return sum(json.loads(x)["n"] for x in open(LEDGER) if x.strip())


def protocol_of(track_name: str, suffix: str) -> str | None:
    for p in PROTOCOLS:
        if track_name == f"{p}_{suffix}":
            return p
    return None


def readout(tidy: pd.DataFrame, suffix: str, all_tracks: bool = False) -> pd.Series:
    """Per variant: mean raw_score over usable SC_ATAC tracks (primary: the three reference protocols)."""
    t = tidy[(tidy.output_type == "SC_ATAC") & tidy.imputation_mode.isin(USABLE_MODES)]
    if not all_tracks:
        t = t[t.track_name.map(lambda n: protocol_of(n, suffix) is not None)]
    return t.groupby(t.variant_id.str.replace(">", ":", regex=False)).raw_score.mean().rename_axis("uid")


def read_dir(d: pathlib.Path) -> pd.DataFrame:
    fs = sorted(glob.glob(str(d / "chunk_*.parquet")))
    return pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()


def reference(d: pathlib.Path) -> pd.Series:
    """Stage B reference readout: mean raw_score over observed SC_ATAC tracks (as Stage B HEP_SCATAC)."""
    t = read_dir(d)
    t = t[(t.output_type == "SC_ATAC") & (t.imputation_mode == "observed")]
    return t.groupby(t.variant_id.str.replace(">", ":", regex=False)).raw_score.mean().rename_axis("uid")


def loo_nearest_centroid_bal_acc(E: np.ndarray, y: np.ndarray) -> float:
    """Leave-one-out nearest-centroid (cosine) balanced accuracy for two labels. No fitted model beyond means."""
    En = E / np.linalg.norm(E, axis=1, keepdims=True)
    labs = np.unique(y)
    assert len(labs) == 2
    sums = {l: En[y == l].sum(0) for l in labs}
    cnt = {l: int((y == l).sum()) for l in labs}
    correct = {l: 0 for l in labs}
    for i in range(len(y)):
        best, bl = -np.inf, None
        for l in labs:
            s, n = sums[l], cnt[l]
            if y[i] == l:
                s, n = s - En[i], n - 1
            c = s / n
            sim = float(En[i] @ c / np.linalg.norm(c))
            if sim > best:
                best, bl = sim, l
        correct[y[i]] += int(bl == y[i])
    return float(np.mean([correct[l] / cnt[l] for l in labs]))


def cosine_dist(a: np.ndarray, b: np.ndarray) -> float:
    return float(1 - a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))

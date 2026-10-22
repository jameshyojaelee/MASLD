"""Shared paths, design table and checks for V5 (GSE298719). Follows PRESPEC_V5.md.

Copied helpers (not imported, sources unedited):
  log_normalize, hc3 <- Analysis/Hypotheses/H90-20261015T120337Z/scripts/
"""
from __future__ import annotations

import gzip
import hashlib
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = Path(__file__).resolve().parent.parent
OUT = RUN / "outputs"
DL = RUN / "downloads_untrusted"
MAT = DL / "matrices"
SUPP = DL / "MOESM4.xlsx"
DERIVED = RUN / "derived"

HOTSPOT = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
MEMBERSHIP = HOTSPOT / "program_membership_v2.tsv"
MEMBERSHIP_SHA = "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b"
CT_MODEL = Path("/gpfs/commons/home/jameslee/.celltypist/data/models/Healthy_Human_Liver.pkl")
CT_MODEL_SHA = "1b9999ba474e9423c583802a3482a900537c488d0e5973dbe20fa702446b8431"
V3_SCRIPT = ROOT / "scripts/analysis/hypotheses_2026q4/V3_program_holdout/annotate_pseudobulk.py"

SEED = 20260923
CONFIRMED = "gencode_v49_unique_symbol_confirmed"
PRIMARY = "hotspot_hepatocytes_48f39dd4d817a10e"    # Ductular injury (BICC1), module 20
SECONDARY = "hotspot_hepatocytes_f05c535ae5bbc0b9"  # Stromal ECM (IGFBP7), module 8
NC1 = "hotspot_hepatocytes_bece6cd6d8fc9c67"        # Sinusoidal endothelial (STAB2), module 2
LINEAGE_LABEL = {
    "hepatocytes": "Hepatocytes",
    "fibroblasts": "Fibroblasts",
    "cholangiocytes": "Cholangiocytes",
    "macrophages": "Macrophages",
    "tcells": "T cells",
}
MIN_CELLS = 30
DET_MIN = 0.01
MIN_GENES, MIN_WEIGHT = 8, 0.80
Y_PROBES = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y"]

# Supplementary Table 12, copied into PRESPEC_V5.md before any expression was read.
DESIGN = pd.DataFrame(
    [("01", "M", "Healthy", np.nan, "E30230224_2"), ("02", "M", "Healthy", np.nan, "E36240625_1"),
     ("03", "F", "Healthy", np.nan, "E35250325_1"), ("04", "M", "Healthy", np.nan, "E35250325_1"),
     ("05", "F", "Healthy", np.nan, "E36240625_1"), ("06", "M", "Healthy", np.nan, "E36240625_1"),
     ("07", "F", "MASLD", 3, "E35250325_1"), ("08", "M", "MASLD", 3, "E35250325_1"),
     ("09", "M", "MASLD", 2, "E35250325_1"), ("10", "F", "MASLD", 3, "E33260824_1"),
     ("11", "M", "MASLD", 3, "E35250325_1"), ("12", "F", "MASLD", 4, "E33260824_1"),
     ("13", "F", "MASLD", 4, "E36240625_1")],
    columns=["id", "sex", "group", "fibrosis", "plate"])


def sha256(path: Path, block: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while b := fh.read(block):
            h.update(b)
    return h.hexdigest()


def check_prespec() -> str:
    exp = (RUN / "PRESPEC_V5.md.sha256").read_text().split()[0]
    got = sha256(RUN / "PRESPEC_V5.md")
    if got != exp:
        sys.exit(f"STOP: PRESPEC_V5.md hash differs: {got}")
    if sha256(MEMBERSHIP) != MEMBERSHIP_SHA:
        sys.exit("STOP: membership hash differs")
    return got


def load_design() -> pd.DataFrame:
    """GSM records joined to Supplementary Table 12; stops on any disagreement."""
    rows = []
    for soft in sorted(DL.glob("GSM98599*.soft")):
        txt = soft.read_text()
        title = re.search(r"^!Sample_title = (.+)$", txt, re.M).group(1).strip()
        m = re.fullmatch(r"FRPBiopsy_(\d\d)_(Male|Female)_(Healthy|MASLD)", title)
        if m is None:
            sys.exit(f"STOP G0: unexpected title {title}")
        rows.append(dict(gsm=soft.stem, title=title, id=m.group(1), sex_geo=m.group(2)[0], group_geo=m.group(3)))
    geo = pd.DataFrame(rows)
    supp = pd.read_excel(SUPP, sheet_name="Suppl. Table 12", header=2)
    supp = supp[supp["ID"].astype(str).str.startswith("FRPBiopsy_")].copy()
    supp["id"] = supp["ID"].str[-2:]
    d = DESIGN.merge(geo, on="id", how="inner").merge(
        supp[["id", "Sex", "Diagnosis.Simplified", "Plate", "Age", "BMI"]], on="id", how="inner")
    ok = (len(d) == 13 and (d.sex == d.sex_geo).all() and (d.sex == d.Sex).all()
          and (d.group == d.group_geo).all() and (d.group == d["Diagnosis.Simplified"]).all()
          and (d.plate == d.Plate).all())
    if not ok:
        print(d.to_string())
        sys.exit("STOP G0: design table disagrees with GEO or Supplementary Table 12")
    d["sex"] = d.sex.map({"M": "male", "F": "female"})
    d["sex_f"] = (d.sex == "female").astype(int)
    d["disease"] = (d.group == "MASLD").astype(int)
    d = d.rename(columns={"Age": "age", "BMI": "bmi"})
    return d[["gsm", "title", "id", "sex", "sex_f", "group", "disease", "fibrosis", "plate", "age", "bmi"]] \
        .sort_values("gsm").reset_index(drop=True)


def membership() -> pd.DataFrame:
    mem = pd.read_csv(MEMBERSHIP, sep="\t")
    return mem[mem.mapped_symbol_status == CONFIRMED]


def weights(mem: pd.DataFrame, uid: str) -> pd.DataFrame:
    m = mem[mem.program_uid == uid]
    return m.groupby("mapped_symbol", as_index=False)["original_l1_weight"].sum()


def log_normalize(counts: sparse.csr_matrix, library: np.ndarray) -> sparse.csr_matrix:
    scale = np.divide(1e4, library, out=np.zeros_like(library, dtype=float), where=library > 0)
    norm = (sparse.diags(scale.astype(np.float32)) @ counts).tocsr()
    norm.data = np.log1p(norm.data)
    return norm


def hc3(y: np.ndarray, X: np.ndarray):
    XtX = np.linalg.inv(X.T @ X)
    b = XtX @ X.T @ y
    e = y - X @ b
    h = np.einsum("ij,jk,ik->i", X, XtX, X)
    meat = X.T @ (X * (e ** 2 / (1 - h) ** 2)[:, None])
    V = XtX @ meat @ XtX
    df = len(y) - X.shape[1]
    sigma = np.sqrt((e ** 2).sum() / df)
    return b, np.sqrt(np.diag(V)), sigma, df, h


def read_gz_lines(p: Path) -> list[str]:
    with gzip.open(p, "rt") as fh:
        return fh.read().splitlines()

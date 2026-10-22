"""Shared paths, hash checks and copied helper functions for H90.

Copied functions (not imported, sources unedited):
  build_strip_mask   <- Analysis/SingleCell/scripts/400_mcp/410_prep_atlas_cnmf.py
  strip rules        <- Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py (strip_confounders)
  build_sample_to_donor <- Analysis/SingleCell/scripts/hotspot_modules/511_pooled_cell_donor_scores.py
  log_normalize / standardize_columns <- Analysis/Multimodal_Program_Projection/scripts/04_spatial_projection.py
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
PRESPEC_SRC = ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H90.md"
PRESPEC_SHA = "ce61a0eb5c6820a5b1f0bb3a0b7cd6d3a63358d1cb8b82dc5b12815e8d62da27"
GSE_DIR = ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gse289173"
GATE = GSE_DIR / "gate"
WRITE = GSE_DIR / "write"
DATA = ROOT / "data/public_r4/GSE289173"
RAW = DATA / "raw"
DERIVED = DATA / f"derived_{RUN.name}"
TAR = RAW / "GSE289173_RAW.tar"
TAR_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289173/suppl/GSE289173_RAW.tar"
TAR_BYTES = 2_076_282_880
FEATURES_FL514 = DATA / "features/GSM8785686_FL514_features.tsv.gz"
FEATURES_SHA = "6af33116525d7a1f3b9d13ff6d637228092357ae8503710bf2e5bc95d6751e0b"

HOTSPOT = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot"
MEMBERSHIP = HOTSPOT / "program_membership_v2.tsv"
FIG4C_SCORES = HOTSPOT / "donor_program_scores_primary.tsv"
DISC_META = ROOT / "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
REGISTRY = ROOT / "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv"
PRIMARY_EFFECTS = ROOT / "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/primary_effects.tsv"
DISC_H5AD = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/atlas_cnmf_hepatocytes.h5ad"
CT_MODEL = Path("/gpfs/commons/home/jameslee/.celltypist/data/models/Healthy_Human_Liver.pkl")
GENE_META = ROOT / "data/gencode_v49_gene_metadata.tsv.gz"

SEED = 20260923
PRIMARY = ("hepatocytes", 20)   # Ductular injury (BICC1)
SECONDARY = ("hepatocytes", 8)  # Stromal ECM (IGFBP7)
NC1 = ("hepatocytes", 2)        # Sinusoidal endothelial (STAB2)
LINEAGE_LABEL = {
    "hepatocytes": "Hepatocytes",
    "fibroblasts": "Fibroblasts",
    "cholangiocytes": "Cholangiocytes",
    "macrophages": "Macrophages",
    "tcells": "T cells",
}
MIN_NUCLEI = 30  # MIN_CELLS_PER_DONOR_CT in 410_prep_atlas_cnmf.py


def sha256(path: Path, block: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def check_inputs(log=print) -> None:
    """Stop if the prespec copy or any listed input hash differs."""
    copy = RUN / "PRESPEC_H90.md"
    for p in (copy, PRESPEC_SRC):
        got = sha256(p)
        if got != PRESPEC_SHA:
            sys.exit(f"STOP: prespec hash differs for {p}: {got}")
    log(f"prespec sha256 OK ({PRESPEC_SHA}) for run copy and plans file")
    n = 0
    for line in (RUN / "inputs.sha256").read_text().splitlines():
        if not line.strip():
            continue
        exp, rel = line.split(None, 1)
        p = Path(rel.strip())
        p = p if p.is_absolute() else ROOT / p
        got = sha256(p)
        if got != exp:
            sys.exit(f"STOP: input hash differs: {p} expected {exp} got {got}")
        n += 1
    # inputs.sha256 lists 25 files: the 24 checked by the audit with sha256sum -c plus the
    # discovery hepatocyte h5ad, which the audit verified separately.
    if n != 25:
        sys.exit(f"STOP: expected 25 input hashes, checked {n}")
    log(f"{n} input hashes OK")


# ---- strip rules (copied) -------------------------------------------------
RIBO_RE = re.compile(r"^(RPS|RPL|MRPS|MRPL)[0-9]")
IG_RE = re.compile(r"^IG[HKL][VJCD]")
HB_LIST = {"HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"}
X_ESCAPEE = {"XIST", "TSIX", "KDM6A", "DDX3X", "EIF1AX", "UTX", "UBA1", "RLIM", "ZFX",
             "RPS4X", "EIF2S3", "TXLNG"}
MT_PREFIX = "MT-"
CONFOUNDER_PREFIXES = ("RPS", "RPL", "MT-", "MRPS", "MRPL")
CONFOUNDER_HB = ["HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBM", "HBQ1", "HBZ"]
CONFOUNDER_ESC = ["XIST", "TSIX", "KDM6A", "DDX3X", "DDX3Y", "EIF1AY", "USP9Y", "UTY", "RPS4Y1",
                  "ZFY", "KDM5D", "NLGN4Y"]


def build_strip_mask(symbols: pd.Index, gene_meta: pd.DataFrame) -> np.ndarray:
    """410 build_strip_mask on a symbol index. True = DROP."""
    var = pd.DataFrame(index=pd.Index(symbols, name="gene_name"))
    gm = gene_meta[["gene_name", "chromosome", "gene_biotype"]].drop_duplicates("gene_name")
    var = var.merge(gm, how="left", left_index=True, right_on="gene_name").set_index("gene_name")
    sym = pd.Index(var.index.astype(str))
    drop = np.zeros(len(sym), dtype=bool)
    drop |= var["chromosome"].eq("chrY").fillna(False).values
    drop |= sym.str.startswith(MT_PREFIX)
    drop |= var["chromosome"].eq("chrM").fillna(False).values
    drop |= var["gene_biotype"].isin(["rRNA", "Mt_rRNA", "Mt_tRNA"]).fillna(False).values
    drop |= sym.to_series().str.match(RIBO_RE).fillna(False).values
    drop |= sym.to_series().str.match(IG_RE).fillna(False).values
    drop |= sym.isin(HB_LIST)
    drop |= sym.isin(X_ESCAPEE)
    return drop


def strip_501_keep(symbols: pd.Index, gene_meta: pd.DataFrame) -> np.ndarray:
    """501 strip_confounders keep mask (protein_coding/lncRNA, not chrY/chrM, prefix/list strip)."""
    keep_gene = set(gene_meta.loc[gene_meta["gene_biotype"].isin({"protein_coding", "lncRNA"}), "gene_name"])
    keep_gene -= set(gene_meta.loc[gene_meta["chromosome"].isin(["chrY", "chrM"]), "gene_name"])
    v = np.asarray(symbols, dtype=str)
    return (np.isin(v, list(keep_gene))
            & ~np.array([n.startswith(CONFOUNDER_PREFIXES) for n in v])
            & ~np.isin(v, CONFOUNDER_HB) & ~np.isin(v, CONFOUNDER_ESC))


def library_gene_set(gse_symbols: pd.Index, disc_var: pd.Index) -> list[str]:
    """Genes used for library size on both sides: GSE289173 symbols after both
    discovery strips, intersected with the discovery hepatocyte var."""
    gm = pd.read_csv(GENE_META, sep="\t")
    sym = pd.Index(pd.unique(np.asarray(gse_symbols, dtype=str)))
    keep = (~build_strip_mask(sym, gm)) & strip_501_keep(sym, gm)
    kept = sym[keep]
    return sorted(set(kept).intersection(set(map(str, disc_var))))


# ---- run -> donor map (copied from 511) -------------------------------------
DONOR_PAIRING = {
    "GSE244832": "data/GSE244832/metadata/donor_pairing.csv",
    "GSE185477": "data/GSE185477/metadata/donor_pairing.csv",
    "GSE202379": "data/GSE202379/metadata/donor_pairing.csv",
    "GSE136103": "data/GSE136103/metadata/donor_pairing.csv",
}


def build_sample_to_donor() -> dict[str, str]:
    m: dict[str, str] = {}
    for ds, rel in DONOR_PAIRING.items():
        dp = pd.read_csv(ROOT / rel, dtype=str)
        for _, row in dp.iterrows():
            donor = f"{ds}_{row['donor_id']}"
            for srr in str(row["rna_srrs"]).split(";"):
                srr = srr.strip()
                if srr:
                    m[srr] = donor
    runs = pd.read_csv(ROOT / "Liver_Atlas/metadata/SraRunTable.csv", dtype=str, usecols=["Run", "shortfilename"])
    samples = pd.read_csv(ROOT / "Liver_Atlas/metadata/GSE192740_sampleInfo_scRNAseq.tsv", sep="\t", dtype=str,
                          usecols=["characteristics: shortFileName", "title"]).rename(
        columns={"characteristics: shortFileName": "shortfilename"})
    liver = runs.merge(samples, on="shortfilename", how="inner")
    liver["donor_short"] = liver["title"].map(
        lambda value: (mm.group(0) if (mm := re.search(r"H[0-9]+", value)) else None))
    liver = liver.dropna(subset=["donor_short"])
    if liver["Run"].nunique() != 48 or liver["donor_short"].nunique() != 19:
        raise RuntimeError("Liver Atlas source run-to-donor cardinality drift")
    liver_map = dict(zip(liver["Run"], "Liver_Atlas_" + liver["donor_short"]))
    if set(m).intersection(liver_map):
        raise RuntimeError("A sequencing library maps to more than one biological donor")
    m.update(liver_map)
    return m


# ---- projection helpers (copied) ---------------------------------------------
def log_normalize(counts: sparse.csr_matrix, library: np.ndarray) -> sparse.csr_matrix:
    scale = np.divide(1e4, library, out=np.zeros_like(library, dtype=float), where=library > 0)
    norm = sparse.diags(scale.astype(np.float32)) @ counts
    norm = norm.tocsr()
    norm.data = np.log1p(norm.data)
    return norm


def donor_mean_z(norm: sparse.csr_matrix, donor_codes: np.ndarray, n_donor: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Donor means of per-nucleus z (genes standardized over all rows, ddof=1).
    Exact by linearity: mean_d(z) = (mean_d(x) - mu) / sd. Returns (donor x gene, mu, sd)."""
    n = norm.shape[0]
    norm = norm.astype(np.float64)
    mu = np.asarray(norm.sum(axis=0)).ravel() / n
    sq = norm.copy()
    sq.data = sq.data.astype(np.float64) ** 2
    ex2 = np.asarray(sq.sum(axis=0)).ravel() / n
    var = (ex2 - mu ** 2) * n / (n - 1)
    sd = np.sqrt(np.clip(var, 0, None))
    sd[~np.isfinite(sd) | (sd <= 0)] = 1.0
    ind = sparse.csr_matrix((np.ones(n), (donor_codes, np.arange(n))), shape=(n_donor, n))
    cnt = np.asarray(ind.sum(axis=1)).ravel()
    dm = np.asarray((ind @ norm).todense(), dtype=np.float64) / cnt[:, None]
    return (dm - mu[None, :]) / sd[None, :], mu, sd


def program_table() -> pd.DataFrame:
    mem = pd.read_csv(MEMBERSHIP, sep="\t")
    return mem


def program_weights(mem: pd.DataFrame, ct: str, mod: int) -> pd.DataFrame:
    """Frozen L1 weights by mapped symbol (unmapped rows dropped), duplicates summed."""
    m = mem[(mem.cell_type == ct) & (mem.module == mod)]
    m = m[m.mapped_symbol.notna()]
    return m.groupby("mapped_symbol", as_index=False)["original_l1_weight"].sum()


def read_gz_lines(p: Path) -> list[str]:
    with gzip.open(p, "rt") as fh:
        return fh.read().splitlines()

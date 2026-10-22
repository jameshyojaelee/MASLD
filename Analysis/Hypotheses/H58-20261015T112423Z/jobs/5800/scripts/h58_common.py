"""Shared loaders and estimators for H58 (PRESPEC_H58.md).

The adopted producer (job 22218816 snapshot) is imported read-only; its
clustering, cluster table and context-mass functions are used unchanged.
"""
from __future__ import annotations

import hashlib
import sys
import zlib
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816"
PRODUCER_DIR = RUN / "snapshot/scripts/figures"
sys.path.insert(0, str(PRODUCER_DIR))
import fig2_multiplicity_correction as m  # noqa: E402

PRESPEC_SHA256 = "adbd01d8e77bef61dd9f588dd67c44d5dfc94494221e2ba64dc8d61d8ae56a58"
SEED = 20260923
N_BOOT = 10_000
N_SPLIT = 200
D, E = m.SCOPES  # direct, enzyme
STRATA = ["EUR", "AMR"]

CORRECTED_COORDS = RUN / "multiplicity/corrected_locus_coordinates.tsv"
FAMILY_TABLE = RUN / "multiplicity/fig2_family_multiplicity.tsv"
BASE = ROOT / "GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321"
ARCH = BASE / "credible_set_pip_architecture.tsv"
MEMBERS = BASE / "credible_set_members_annotated.tsv"
REGISTRY = RUN / "snapshot/GWAS/finemapping/config/gwas_registry.tsv"
LOOKUP = ROOT / "GWAS/finemapping/results/gwas_atac/variant_overlap_summary.csv"
TIER = ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv"
CREDSETS = ROOT / "GWAS/finemapping/results/credible_sets.csv"
COLOC = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"

Z_C = NormalDist().inv_cdf(1 - 0.05 / 6)   # 2.394
Z_POW = NormalDist().inv_cdf(0.80)          # 0.842
MDE_FACTOR = Z_C + Z_POW
LEVELS = {"95": (2.5, 97.5), "98.333": (100 * 0.05 / 6, 100 - 100 * 0.05 / 6),
          "99": (0.5, 99.5)}


def rng(fit_id: str) -> np.random.Generator:
    return np.random.default_rng([SEED, zlib.crc32(fit_id.encode())])


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_prespec(out: Path) -> None:
    got = sha256(out / "PRESPEC_H58.md")
    if got != PRESPEC_SHA256:
        raise SystemExit(f"PRESPEC_H58.md sha256 {got} != {PRESPEC_SHA256}; refusing to run")


# ---------------------------------------------------------------- loaders
def load_coding() -> pd.DataFrame:
    a = pd.read_csv(CORRECTED_COORDS, sep="\t", dtype={"locus": str})
    a["coding_mass"] = a[m.CODING_COLS].sum(axis=1)
    a["arch_total"] = a[m.ARCH_COLS].sum(axis=1)
    return a


def load_gated_members() -> pd.DataFrame:
    arch = pd.read_csv(ARCH, sep="\t", dtype={"locus": str})
    mem = pd.read_csv(MEMBERS, sep="\t", low_memory=False, dtype={"locus": str})
    arch = arch[arch.trait_scope.isin(m.SCOPES)].copy()
    mem = mem[mem.trait_scope.isin(m.SCOPES)].copy()
    arch, mem, _ = m.apply_gate(arch, mem)
    mem["chromosome"] = mem["chromosome"].astype(str)
    mem["position"] = mem["position"].astype(int)
    return mem


def context_per_set(mem: pd.DataFrame) -> pd.DataFrame:
    ps = m.parse_locus(m.build_context_masses(mem), mem)
    anc = mem.drop_duplicates("credible_set_uid").set_index("credible_set_uid").ancestry
    study = mem.drop_duplicates("credible_set_uid").set_index("credible_set_uid").study
    ps["ancestry"] = ps.credible_set_uid.map(anc)
    ps["study"] = ps.credible_set_uid.map(study)
    return ps


def repaired_members(mem: pd.DataFrame, sites: pd.DataFrame) -> pd.DataFrame:
    """Prespec 'Context repair' step 2."""
    rep = mem.merge(sites, on=["chromosome", "position"], how="left", validate="m:1")
    rep.index = mem.index  # left merge keeps row order; restore the member index
    if rep["n_map"].isna().any():
        raise AssertionError("member site missing from re-annotation table")
    tf = {True: "true", False: "false"}
    rep["abc_enhancer_overlap"] = rep["abc_new"].astype(bool).map(tf)
    one = rep["n_map"] == 1
    rep.loc[one, "promoter_proximal"] = rep.loc[one, "promoter_new"].astype(bool).map(tf)
    rep.loc[one, "lineage_accessible"] = rep.loc[one, "accessible_new"].astype(bool).map(tf)
    rep.loc[one, "promoter_context_status"] = "resolved"
    rep.loc[one, "lineage_accessibility_status"] = "resolved"
    return rep


def load_sites(out: Path) -> pd.DataFrame:
    s = pd.read_csv(out / "site_reannotated.tsv", sep="\t")
    s["chromosome"] = s["chromosome"].astype(str)
    s["position"] = s["position"].astype(int)
    for c in ["promoter_new", "accessible_new", "abc_new"]:
        s[c] = s[c].astype(str).str.upper().map({"TRUE": True, "FALSE": False, "NA": False,
                                                 "NAN": False})
    return s


# ---------------------------------------------------------------- stratified design
class Design:
    """Resampling design shared by a family.

    `frame`: instance table with chrom, pos, trait_scope, ancestry and the
    numerator / denominator columns. Pooled physical clusters are formed on
    ALL instances of `frame`; within-stratum clusters on each stratum's
    instances. Units = pooled clusters holding at least one row of `strata`.
    Arrays NUM[u, k], DEN[u, k] with k over (stratum, class) and one block per
    numerator column.
    """

    def __init__(self, frame: pd.DataFrame, num_cols: list[str], den_col: str,
                 strata: list[str], stratum_col: str = "ancestry",
                 pooled_frame: pd.DataFrame | None = None):
        self.num_cols, self.den_col, self.strata = num_cols, den_col, strata
        base = frame if pooled_frame is None else pooled_frame
        pooled = m.assign_clusters(base.copy())
        pid = pooled.set_index("credible_set_uid")["cluster_id"]
        rows = []
        self.within = {}
        for s in strata:
            sub = m.assign_clusters(frame[frame[stratum_col] == s].copy())
            sub["pooled_id"] = sub.credible_set_uid.map(pid)
            nest = sub.groupby("cluster_id").pooled_id.nunique()
            if nest.max() != 1:
                raise AssertionError(f"stratum {s}: within-stratum cluster spans pooled clusters")
            self.within[s] = sub
            g = sub.groupby(["trait_scope", "cluster_id"], as_index=False).agg(
                {**{c: "mean" for c in num_cols}, den_col: "mean", "pooled_id": "first"})
            g["stratum"] = s
            rows.append(g)
        tab = pd.concat(rows, ignore_index=True)
        self.tab = tab
        self.units = sorted(tab.pooled_id.unique())
        self.n_pooled_total = pooled.cluster_id.nunique()
        uidx = {u: i for i, u in enumerate(self.units)}
        self.keys = [(s, x) for s in strata for x in (D, E)]
        kidx = {k: i for i, k in enumerate(self.keys)}
        U, K = len(self.units), len(self.keys)
        self.NUM = {c: np.zeros((U, K)) for c in num_cols}
        self.DEN = np.zeros((U, K))
        for r in tab.itertuples(index=False):
            i, k = uidx[r.pooled_id], kidx[(r.stratum, r.trait_scope)]
            for c in num_cols:
                self.NUM[c][i, k] += getattr(r, c)
            self.DEN[i, k] += getattr(r, den_col)

    def n_clusters(self, s: str, x: str) -> int:
        t = self.tab
        return int(((t.stratum == s) & (t.trait_scope == x)).sum())

    def fractions(self, W: np.ndarray | None = None) -> dict:
        """W: (B, U) count matrix, or None for the observed data.
        Returns {col: (B or 1, K) fraction array} and den sums (B, K)."""
        if W is None:
            W = np.ones((1, len(self.units)))
        den = W @ self.DEN
        out = {}
        with np.errstate(divide="ignore", invalid="ignore"):
            for c in self.num_cols:
                out[c] = np.where(den > 0, (W @ self.NUM[c]) / den, np.nan)
        return out, den

    def k(self, s: str, x: str) -> int:
        return self.keys.index((s, x))


def count_matrix(g: np.random.Generator, n_units: int, n_boot: int = N_BOOT) -> np.ndarray:
    """Producer-style draw (rng.integers over units, n_boot x n_units), as counts."""
    idx = g.integers(0, n_units, size=(n_boot, n_units))
    off = idx + (np.arange(n_boot) * n_units)[:, None]
    return np.bincount(off.ravel(), minlength=n_boot * n_units).reshape(n_boot, n_units).astype(float)


def pct(x: np.ndarray, level: str) -> tuple[float, float]:
    lo, hi = LEVELS[level]
    return float(np.percentile(x, lo)), float(np.percentile(x, hi))


def covers0(x: np.ndarray, level: str) -> bool:
    lo, hi = pct(x, level)
    return lo <= 0 <= hi

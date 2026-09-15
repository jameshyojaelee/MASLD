#!/usr/bin/env python3
"""Verbatim copy of the P5 haplotype-additivity recipe (scripts/analysis/alphagenome_atlas/50_haplotype_additivity.py).

Copied, not imported, because the program boundary forbids touching the atlas tree, and copied verbatim
because F3's whole point is that observed and null must come from ONE loader whose statistic matches the
existing deposit. The only thing F2/F3 change is which pairs go in and how the null is drawn; every
function below is byte-identical in behaviour to the deposit's.

Functions copied: apply_variants, window_mask, fold_mask, summarise, effect_threshold, effect_stratum,
empirical_two_sided_p, bh_reject, stratum_is_testable. Constants copied: FASTA, GTF, LIVER_TERMS,
OUTPUTS, WINDOW, FLANK, EPS, BH_Q, MIN_MATCHED_NULL.
"""

from __future__ import annotations

import gzip
import math
import re

import numpy as np
from alphagenome.models import dna_client

FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
GTF = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/reference-build-21062075/gencode.v49.primary_analysis.annotation.gtf.gz"
LIVER_TERMS = ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]
OUTPUTS = ["RNA_SEQ", "ATAC", "DNASE", "CHIP_HISTONE"]
WINDOW = dna_client.SEQUENCE_LENGTH_1MB
FLANK = 1_000
EPS = 1e-6
BH_Q = 0.10
MIN_MATCHED_NULL = 20
CHANNELS = ("rna", "atac", "dnase", "h3k27ac")


class ContractError(RuntimeError):
    pass


def gene_spans() -> dict[str, tuple[str, int, int]]:
    out = {}
    with gzip.open(GTF, "rt") as h:
        for line in h:
            if line.startswith("#"):
                continue
            p = line.split("\t", 9)
            if p[2] != "gene":
                continue
            m = re.search(r'gene_id "([^"]+)"', p[8])
            if m:
                out[m.group(1).split(".")[0]] = (p[0], int(p[3]), int(p[4]))
    return out


def apply_variants(seq: str, start0: int, variants: list[tuple[int, str, str]]) -> str:
    """1-based positions; substitutions only (SNVs), so length is preserved."""
    s = list(seq)
    for pos, ref, alt in variants:
        i = pos - 1 - start0
        if s[i].upper() != ref.upper():
            raise ContractError(f"reference mismatch at {pos}: sequence has {s[i]}, expected {ref}")
        s[i] = alt
    return "".join(s)


def window_mask(length: int, start0: int, regions: list[tuple[int, int]]) -> np.ndarray:
    m = np.zeros(length, bool)
    for a, b in regions:
        lo, hi = max(0, a - start0), min(length, b - start0)
        if hi > lo:
            m[lo:hi] = True
    return m


def fold_mask(mask: np.ndarray, n_rows: int) -> np.ndarray:
    """Fold a base-resolution readout mask onto a track's own row count (bin kept when any base is selected)."""
    mask = np.asarray(mask, bool)
    if n_rows == mask.shape[0]:
        return mask
    if n_rows <= 0 or mask.shape[0] % n_rows:
        raise ContractError(f"track rows {n_rows} do not divide the base mask of {mask.shape[0]}")
    return mask.reshape(n_rows, mask.shape[0] // n_rows).any(axis=1)


def summarise(output, mask_rna: np.ndarray, mask_local: np.ndarray) -> dict[str, float]:
    vals = {}
    for name, td, mask in (("rna", output.rna_seq, mask_rna), ("atac", output.atac, mask_local),
                           ("dnase", output.dnase, mask_local), ("h3k27ac", output.chip_histone, mask_local)):
        if td is None:
            vals[name] = math.nan
            continue
        v = np.asarray(td.values, np.float64)
        keep = np.ones(v.shape[1], bool)
        if name == "h3k27ac" and td.metadata is not None and "name" in td.metadata:
            keep = td.metadata["name"].astype(str).str.contains("H3K27ac").to_numpy()
            if not keep.any():
                vals[name] = math.nan
                continue
        vals[name] = float(v[fold_mask(mask, v.shape[0])][:, keep].sum())
    return vals


def effect_threshold(null_v1: np.ndarray, null_v2: np.ndarray) -> float:
    """95th percentile of |single-variant effect| over the matched null, both arms pooled."""
    a = np.abs(np.concatenate([np.asarray(null_v1, float), np.asarray(null_v2, float)]))
    a = a[~np.isnan(a)]
    return float(np.quantile(a, 0.95)) if a.size else float("nan")


def effect_stratum(v1: np.ndarray, v2: np.ndarray, threshold: float) -> np.ndarray:
    """Pairs where BOTH variants carry a predicted effect above the null threshold."""
    a, b = np.abs(np.asarray(v1, float)), np.abs(np.asarray(v2, float))
    with np.errstate(invalid="ignore"):
        return (~np.isnan(a)) & (~np.isnan(b)) & (a > threshold) & (b > threshold)


def empirical_two_sided_p(residual: float, null: np.ndarray) -> float:
    """(count of |null| >= |residual|, plus one) / (draws plus one)."""
    nul = np.asarray(null, float)
    nul = nul[~np.isnan(nul)]
    if nul.size == 0 or residual != residual:
        return float("nan")
    return float((np.sum(np.abs(nul) >= abs(float(residual))) + 1) / (nul.size + 1))


def bh_reject(pvals: np.ndarray, q: float) -> np.ndarray:
    """Benjamini-Hochberg step-up over the whole family; returns the rejection mask."""
    p = np.asarray(pvals, float)
    ok = ~np.isnan(p)
    out = np.zeros(p.shape, bool)
    idx = np.flatnonzero(ok)
    if idx.size == 0:
        return out
    order = idx[np.argsort(p[idx], kind="mergesort")]
    m = order.size
    thresh = q * (np.arange(1, m + 1) / m)
    passing = np.flatnonzero(p[order] <= thresh)
    if passing.size:
        out[order[: passing[-1] + 1]] = True
    return out


def stratum_is_testable(n_observed: int, n_null_matched: int, min_null: int = MIN_MATCHED_NULL) -> dict:
    if n_null_matched < min_null:
        return {"testable": False, "n_null_matched": int(n_null_matched), "min_null": int(min_null),
                "reason": (f"only {n_null_matched} of the null pairs meet this stratum's selection "
                           f"(minimum {min_null}); there is no matched null, so no verdict is emitted.")}
    return {"testable": True, "n_null_matched": int(n_null_matched), "min_null": int(min_null), "reason": "matched null available"}

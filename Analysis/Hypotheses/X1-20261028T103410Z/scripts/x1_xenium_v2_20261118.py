#!/usr/bin/env python3
"""X1 (20261118 N09 remediation): Fig 5 genes in intact liver (SCP3359 Xenium, 9 Control vs 9 MASH-cirrhosis donors).

Remediates N09 defects relative to x1_xenium_v2.py:
  1. Removes acute-phase reactant HPX from hepatocyte spillover REF -> primary REF = ["KNG1", "FBP1"],
     with full sensitivity comparison against ["HPX", "KNG1", "FBP1"] written to ref_sensitivity_hpx.tsv.
  2. Eliminates carrier selection circularity by selecting carrier on the permutation-invariant 18-donor
     grand mean after enforcing made_all spillover + donor sufficiency BEFORE carrier selection.
  3. Excludes underpowered (<7 eligible donors, e.g. SPI1) and LABEL_DEFINING markers
     ({"KRT8", "KRT18", "MZB1", "LYVE1"}) from primary within_holm / within_bh multiplicity families
     while retaining all rows with explicit in_primary_within_holm and *_primary vs *_all_sens columns.
  4. Points FIG5B to WS4A's rebuilt proteomics_20261118 triage table and PB / snRNA lineage specificity
     to WS3A's donor-deduplicated single-cell outputs (NUCLEUS-FIX-20261107).
  5. Allows CATALOG override via X1_CATALOG_PATH (and optional FIG5B override via X1_FIG5B_PATH).

Unit = donor. Usage: x1_xenium_v2_20261118.py <out_dir>
"""
import importlib.util
import itertools
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io
from scipy.spatial import ConvexHull, cKDTree

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
SRC = f"{ROOT}/Analysis/Hypotheses/V7-20261021T235717Z/downloads_untrusted/scp3359/SCP3359"
EXPR = f"{SRC}/expression/68f4186f4e6da7a7e1b0edbe"
DEFAULT_CATALOG = (
    f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/"
    "catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.tsv.gz"
)
CATALOG = os.environ.get("X1_CATALOG_PATH", DEFAULT_CATALOG)
DEFAULT_FIG5B = (
    f"{ROOT}/Analysis/Multimodal_Program_Projection/results/proteomics_20261118/"
    "figure5/panels/data/fig5b_protein_triage_genes.tsv"
)
FIG5B = os.environ.get("X1_FIG5B_PATH", DEFAULT_FIG5B)

# Donor-deduplicated single-cell source (WS3A / NUCLEUS-FIX-20261107)
SC_CAND_DIR = f"{ROOT}/Analysis/Hypotheses/NUCLEUS-FIX-20261107/candidates/sc_20261118"
SC_FIG4C_SPEC = f"{SC_CAND_DIR}/fig4c_gene_specificity_16lineages.tsv.gz"
SC_UC_DIR = f"{ROOT}/Analysis/Hypotheses/NUCLEUS-FIX-20261107/unique_counts"
PB = f"{SC_UC_DIR}/pseudobulk/{{}}_pseudobulk.csv"
SC_DONOR_META = f"{SC_UC_DIR}/metadata/donor_metadata.tsv"
SC_LIB_META = f"{SC_UC_DIR}/metadata/library_metadata.tsv"

OUT = sys.argv[1]
Path(OUT).mkdir(parents=True, exist_ok=True)

SEED = 20261027
N_BOOT = 10_000
D5_DOWN = ["CYP3A4", "TAT", "GLS2", "MET", "LYVE1"]
D5_UP = ["FASN", "KRT8", "KRT18", "MZB1", "PLIN2", "PLIN4", "ACTA2", "ANXA5", "ERBB2"]
D5 = D5_DOWN + D5_UP
E5 = ["THRB", "SPI1", "KLRC1", "TNFSF10", "ADH4", "ADH1C", "ANPEP", "TM4SF4", "SCGN", "MEST", "MTOR", "CFHR1"]
REF = ["KNG1", "FBP1"]
REF_WITH_HPX = ["HPX", "KNG1", "FBP1"]
LABEL_DEFINING = {"KRT8", "KRT18", "MZB1", "LYVE1"}
HEP = "Hepatocyte"
TYPES = [
    "Hepatocyte", "Cholangiocyte", "LSEC", "Endothelial", "Stellate", "VSMC",
    "Kupffer", "Macrophage", "B", "T.CD4", "T.CD8", "T.Innate",
]
# v2: the source files DNASE1L3+FCN2+LYVE1-high, PDGFRB-low cells under broad "Stellate"; they are sinusoidal
# endothelium (label_check_refined_means.tsv). VSMC are split out of "Stellate" as well.
REFINED_TO_V2 = {"Stellate.DNASE1L3+CD36+": "LSEC", "VSMC.MYH11+": "VSMC"}
SN_LINEAGES = [
    "Hepatocytes", "Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Macrophages",
    "Mono+mono_derived_cells", "T_cells", "B_cells", "Plasma_cells", "Resident_NK",
    "Circulating_NK_NKT",
]
SN_TO_XEN = {
    "Hepatocytes": {"Hepatocyte"},
    "Cholangiocytes": {"Cholangiocyte"},
    "Endothelial_cells": {"Endothelial", "LSEC"},
    "Fibroblasts": {"Stellate", "VSMC"},
    "Macrophages": {"Kupffer", "Macrophage"},
    "Mono+mono_derived_cells": {"Macrophage", "Kupffer"},
    "T_cells": {"T.CD4", "T.CD8"},
    "Resident_NK": {"T.Innate"},
    "Circulating_NK_NKT": {"T.Innate"},
    "B_cells": {"B"},
    "Plasma_cells": {"B"},
}
OLD_CONTROLS = ["11-00C", "11-025"]
INTERIOR_NONHEP_UM = 15.0
MIN_INTERIOR_CELLS = 30
MIN_INTERIOR_HEP = 500
MIN_ARM_DONORS = 7
MIN_PRESENT_DONORS = 2 * MIN_ARM_DONORS  # 14 of 18 donors


def log(msg):
    print(msg, flush=True)


def write(df, name, index=False):
    df.to_csv(f"{OUT}/{name}", sep="\t", index=index)


def to_bool(series):
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


log(f"OUT={OUT}")
log(f"REF={REF} (sensitivity REF_WITH_HPX={REF_WITH_HPX})")
log(f"FIG5B={FIG5B}")
log(f"CATALOG={CATALOG}")
log(f"SC_UC_DIR={SC_UC_DIR}")
log(f"SC_FIG4C_SPEC={SC_FIG4C_SPEC}")

# ---------------------------------------------------------------- load
genes = [l.strip() for l in open(f"{EXPR}/counts_features.tsv") if l.strip()]
barcodes = [l.strip() for l in open(f"{EXPR}/counts_barcodes.tsv") if l.strip()]
gi = {g: i for i, g in enumerate(genes)}
for g in D5 + E5 + REF_WITH_HPX + ["KRT7"]:
    assert g in gi, g
log(f"genes {len(genes)}, cells {len(barcodes)}")

meta = pd.read_csv(
    f"{SRC}/cluster/coords.txt",
    sep="\t",
    skiprows=[1],
    usecols=[
        "NAME", "X", "Y", "broad_cell_type", "refined_cell_type", "condition",
        "age_at_sample_collection", "donor_id",
    ],
    dtype={"NAME": str, "broad_cell_type": str, "condition": str, "donor_id": str},
)
meta = meta.set_index("NAME")
assert meta.index.is_unique
meta = meta.loc[barcodes]
meta["broad_cell_type"] = [
    REFINED_TO_V2.get(r, b) for r, b in zip(meta.refined_cell_type, meta.broad_cell_type)
]
assert set(meta.broad_cell_type) == set(TYPES), set(meta.broad_cell_type)

X = io.mmread(f"{EXPR}/counts_matrix.mtx").T.tocsr().astype(np.float64)  # cells x genes
assert X.shape == (len(barcodes), len(genes))
log(f"matrix loaded, nnz {X.nnz}")

donor_tab = (
    meta.groupby("donor_id")
    .agg(condition=("condition", "first"), age=("age_at_sample_collection", "first"), n_cells=("X", "size"))
    .reset_index()
)
assert (meta.groupby("donor_id").condition.nunique() == 1).all()
donors = donor_tab.sort_values(["condition", "donor_id"]).donor_id.tolist()
is_mash = np.array([donor_tab.set_index("donor_id").condition[d] == "MASH" for d in donors])
log(f"donors {len(donors)}: {is_mash.sum()} MASH, {(~is_mash).sum()} Control")

# ---------------------------------------------------------------- interior cells
ctype = meta.broad_cell_type.values
xy_all = meta[["X", "Y"]].to_numpy()
interior = np.zeros(len(meta), bool)
area = {}
for d, idx in meta.groupby("donor_id").indices.items():
    xy = xy_all[idx]
    hep = ctype[idx] == HEP
    tree = cKDTree(xy)
    _, nn = tree.query(xy[hep], k=7)
    interior[idx[hep]] = hep[nn[:, 1:]].all(1)
    dist, _ = cKDTree(xy[hep]).query(xy[~hep], k=1)
    interior[idx[~hep]] = dist > INTERIOR_NONHEP_UM
    hull = ConvexHull(xy)
    occ = np.unique(np.floor(xy / 50.0).astype(np.int64), axis=0).shape[0] * 50.0 ** 2
    area[d] = (hull.volume / 1e6, occ / 1e6)  # mm2; 2-D hull volume is area

# ---------------------------------------------------------------- donor x type summaries
nD, nK, nG = len(donors), len(TYPES), len(genes)
n_all = np.zeros((nD, nK))
n_int = np.zeros((nD, nK))
m_all = np.zeros((nD, nK, nG))
m_int = np.full((nD, nK, nG), np.nan)
det_all = np.zeros((nD, nK, nG))
det_int = np.full((nD, nK, nG), np.nan)
tot_counts = np.zeros((nD, nG))
grp = meta.groupby(["donor_id", "broad_cell_type"]).indices
for a, d in enumerate(donors):
    didx = np.where(meta.donor_id.values == d)[0]
    tot_counts[a] = np.asarray(X[didx].sum(0)).ravel()
    for b, k in enumerate(TYPES):
        idx = grp.get((d, k), np.array([], int))
        n_all[a, b] = len(idx)
        if len(idx):
            sub = X[idx]
            m_all[a, b] = np.asarray(sub.mean(0)).ravel()
            det_all[a, b] = np.asarray((sub > 0).mean(0)).ravel()
        ii = idx[interior[idx]]
        n_int[a, b] = len(ii)
        if k != HEP and len(ii) < MIN_INTERIOR_CELLS:
            ii = idx  # per-donor fallback to all cells of this type (PRESPEC)
        if len(ii):
            sub = X[ii]
            m_int[a, b] = np.asarray(sub.mean(0)).ravel()
            det_int[a, b] = np.asarray((sub > 0).mean(0)).ravel()

N_d = n_all.sum(1)
share = n_all / N_d[:, None]
h = TYPES.index(HEP)

# Gate 1
g1 = pd.DataFrame({
    "donor_id": donors,
    "condition": np.where(is_mash, "MASH", "Control"),
    "n_cells": N_d,
    "n_hepatocytes": n_all[:, h],
    "n_interior_hepatocytes": n_int[:, h],
    "hepatocyte_share": share[:, h],
    "area_hull_mm2": [area[d][0] for d in donors],
    "area_occupancy50_mm2": [area[d][1] for d in donors],
})
g1["pass"] = g1.n_interior_hepatocytes >= MIN_INTERIOR_HEP
write(g1, "gate1_donors.tsv")
keep = g1["pass"].to_numpy()
log(f"Gate 1: {keep.sum()} of {nD} donors pass "
    f"({(keep & is_mash).sum()} MASH, {(keep & ~is_mash).sum()} Control)")
if (keep & is_mash).sum() < MIN_ARM_DONORS or (keep & ~is_mash).sum() < MIN_ARM_DONORS:
    log("GATE 1 FAILED: fewer than 7 donors in an arm. Stop.")
    sys.exit(0)
if not keep.all():
    log(f"Dropped donors: {[d for d, k in zip(donors, keep) if not k]}")
sel = np.where(keep)[0]
donors = [donors[i] for i in sel]
is_mash = is_mash[sel]
n_all, n_int, m_all, m_int = n_all[sel], n_int[sel], m_all[sel], m_int[sel]
det_all, det_int, tot_counts, share, N_d = det_all[sel], det_int[sel], tot_counts[sel], share[sel], N_d[sel]
nD = len(donors)
ctrl, mash = ~is_mash, is_mash

# Preserve donor aggregates before any disease-label-dependent cell-set choice.
# This supports exact replay of that choice and carrier selection under the null.
np.savez_compressed(
    f"{OUT}/donor_aggregates_selection.npz",
    genes=np.asarray(genes),
    donors=np.asarray(donors),
    types=np.asarray(TYPES),
    ref_genes=np.asarray(REF),
    ref_genes_with_hpx=np.asarray(REF_WITH_HPX),
    is_mash=is_mash,
    n_all=n_all,
    n_int=n_int,
    m_all=m_all,
    m_int=m_int,
    det_all=det_all,
    det_int=det_int,
    tot_counts=tot_counts,
)

# cell set per type
use_interior = {}
for b, k in enumerate(TYPES):
    if k == HEP:
        use_interior[k] = True
        continue
    ok = n_int[:, b] >= MIN_INTERIOR_CELLS
    use_interior[k] = bool((ok & mash).sum() >= MIN_ARM_DONORS and (ok & ctrl).sum() >= MIN_ARM_DONORS)
cellset = pd.DataFrame({
    "type": TYPES,
    "cell_set": ["interior" if use_interior[k] else "all" for k in TYPES],
    "control_median_interior_n": [np.median(n_int[ctrl, b]) for b in range(nK)],
    "mash_median_interior_n": [np.median(n_int[mash, b]) for b in range(nK)],
    "control_median_all_n": [np.median(n_all[ctrl, b]) for b in range(nK)],
    "mash_median_all_n": [np.median(n_all[mash, b]) for b in range(nK)],
})
write(cellset, "cell_set_by_type.tsv")
log(cellset.to_string())


def corrected_levels(m_set, hep_ref, ref_genes=REF):
    """Spillover-corrected per-cell levels; returns corr, expected, leak."""
    ridx = [gi[r] for r in ref_genes]
    leak = np.zeros((nD, nK))
    for b, k in enumerate(TYPES):
        if k == HEP:
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = m_set[:, b][:, ridx] / hep_ref[:, ridx]
        leak[:, b] = np.nanmedian(ratios, axis=1)
    expected = leak[:, :, None] * hep_ref[:, None, :]
    expected[:, h, :] = 0.0
    corr = np.clip(m_set - expected, 0.0, None)
    return corr, expected, leak


use_int_mask = np.array([use_interior[k] for k in TYPES], dtype=bool)
hep_ref = m_int[:, h, :]
m_use = np.where(use_int_mask[None, :, None], m_int, m_all)
m_use = np.nan_to_num(m_use, nan=0.0)
det_use = np.nan_to_num(np.where(use_int_mask[None, :, None], det_int, det_all), nan=0.0)

# Exact used cell count per donor x type (accounting for per-donor fallback when n_int < 30)
n_used_int = np.where(n_int >= MIN_INTERIOR_CELLS, n_int, n_all)
n_used_int[:, h] = n_int[:, h]
n_used_cells = np.where(use_int_mask[None, :], n_used_int, n_all)

corr, expected, leak = corrected_levels(m_use, hep_ref, ref_genes=REF)
m_all_h = m_all.copy()
m_all_h[:, h] = hep_ref
corr_allcells, expected_allcells, _ = corrected_levels(m_all_h, hep_ref, ref_genes=REF)
write(pd.DataFrame(leak, index=donors, columns=TYPES), "spillover_leak_by_donor_type.tsv", index=True)

# Document per-type spillover leak sensitivity to removing HPX from REF
_, _, leak_with_hpx = corrected_levels(m_use, hep_ref, ref_genes=REF_WITH_HPX)
leak_cmp = pd.DataFrame({
    "type": TYPES,
    "cell_set": ["interior" if use_interior[k] else "all" for k in TYPES],
    "leak_median_no_hpx": [float(np.median(leak[:, b])) for b in range(nK)],
    "leak_mean_no_hpx": [float(np.mean(leak[:, b])) for b in range(nK)],
    "leak_median_with_hpx": [float(np.median(leak_with_hpx[:, b])) for b in range(nK)],
    "leak_mean_with_hpx": [float(np.mean(leak_with_hpx[:, b])) for b in range(nK)],
    "delta_leak_median_no_minus_with_hpx": [
        float(np.median(leak[:, b]) - np.median(leak_with_hpx[:, b])) for b in range(nK)
    ],
    "delta_leak_mean_no_minus_with_hpx": [
        float(np.mean(leak[:, b]) - np.mean(leak_with_hpx[:, b])) for b in range(nK)
    ],
})
write(leak_cmp, "ref_leak_by_type_hpx_comparison.tsv")

# long table for later rendering
rows = []
for a, d in enumerate(donors):
    for b, k in enumerate(TYPES):
        for g in sorted(set(D5 + E5 + REF_WITH_HPX + ["KRT7"])):
            j = gi[g]
            rows.append((
                d, "MASH" if is_mash[a] else "Control", k, g, n_all[a, b], n_int[a, b],
                share[a, b], m_all[a, b, j], m_use[a, b, j], expected[a, b, j], corr[a, b, j],
                det_use[a, b, j],
            ))
write(
    pd.DataFrame(
        rows,
        columns=[
            "donor_id", "condition", "type", "gene", "n_cells", "n_interior",
            "share", "mean_all", "mean_used", "expected_spillover", "mean_corrected",
            "detect_frac_used",
        ],
    ),
    "donor_type_gene_levels.tsv.gz",
)

# ---------------------------------------------------------------- statistics helpers
A = np.array([np.isin(np.arange(nD), c) for c in itertools.combinations(range(nD), int(mash.sum()))])
log(f"permutation splits: {A.shape[0]}")
rng = np.random.default_rng(SEED)
mi, ci = np.where(mash)[0], np.where(ctrl)[0]
BM = rng.integers(0, len(mi), size=(N_BOOT, len(mi)))
BC = rng.integers(0, len(ci), size=(N_BOOT, len(ci)))


def lr(v_m, v_c):
    with np.errstate(divide="ignore", invalid="ignore"):
        m_val, c_val = np.mean(v_m), np.mean(v_c)
        if m_val == 0.0 and c_val == 0.0:
            return 0.0
        return np.log2(m_val / c_val)


def test(v):
    """log2 ratio of arm means, exact permutation p, bootstrap CI, leave-one-out range."""
    v = np.asarray(v, float)
    obs = lr(v[mash], v[ctrl])
    if not np.isfinite(obs):
        return dict(
            log2_ratio=obs, p_perm=np.nan, ci_lo=np.nan, ci_hi=np.nan,
            lodo_min=np.nan, lodo_max=np.nan, lodo_sign_flips=np.nan,
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        pm = (A @ v) / A.sum(1)
        pc = ((~A) @ v) / (~A).sum(1)
        tot = pm + pc
        bd = np.divide(pm - pc, tot, out=np.zeros_like(tot), where=tot > 0)
        obs_m, obs_c = np.mean(v[mash]), np.mean(v[ctrl])
        obs_bd = 0.0 if (obs_m + obs_c) == 0.0 else (obs_m - obs_c) / (obs_m + obs_c)
    p = float(np.mean(np.abs(bd) >= abs(obs_bd) - 1e-12))
    vm, vc = v[mi], v[ci]
    with np.errstate(divide="ignore", invalid="ignore"):
        bm_mean = vm[BM].mean(1)
        bc_mean = vc[BC].mean(1)
        boot = np.where((bm_mean == 0.0) & (bc_mean == 0.0), 0.0, np.log2(bm_mean / bc_mean))
    lo, hi = np.quantile(boot, [0.025, 0.975], method="inverted_cdf")
    lodo = []
    for drop in range(nD):
        k_mask = np.arange(nD) != drop
        lodo.append(lr(v[k_mask & mash], v[k_mask & ctrl]))
    lodo = np.array(lodo)
    finite_lodo = lodo[np.isfinite(lodo)]
    return dict(
        log2_ratio=float(obs),
        p_perm=p,
        ci_lo=float(lo),
        ci_hi=float(hi),
        lodo_min=float(np.min(finite_lodo)) if len(finite_lodo) else np.nan,
        lodo_max=float(np.max(finite_lodo)) if len(finite_lodo) else np.nan,
        lodo_sign_flips=int(np.sum(np.sign(finite_lodo) != np.sign(obs))) if len(finite_lodo) else np.nan,
    )


def holm(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    if len(q) == 0:
        return out
    order = np.argsort(q)
    adj = np.maximum.accumulate((len(q) - np.arange(len(q))) * q[order]).clip(max=1.0)
    res = np.empty_like(q)
    res[order] = adj
    out[ok] = res
    return out


def holm_subset(p, mask):
    p = np.asarray(p, float)
    mask = np.asarray(mask, bool)
    out = np.full_like(p, np.nan)
    if mask.any():
        out[mask] = holm(p[mask])
    return out


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    n = len(q)
    if n == 0:
        return out
    order = np.argsort(q)
    adj = np.minimum.accumulate((q[order] * n / np.arange(1, n + 1))[::-1])[::-1].clip(max=1.0)
    res = np.empty_like(q)
    res[order] = adj
    out[ok] = res
    return out


def bh_subset(p, mask):
    p = np.asarray(p, float)
    mask = np.asarray(mask, bool)
    out = np.full_like(p, np.nan)
    if mask.any():
        out[mask] = bh(p[mask])
    return out


def primary_exclusion_reasons(label_def, det_ctrl):
    reasons = []
    for ld, dc in zip(label_def, det_ctrl):
        if ld and dc < MIN_ARM_DONORS:
            reasons.append("label_defining_and_underpowered_lt7_donors")
        elif ld:
            reasons.append("label_defining")
        elif dc < MIN_ARM_DONORS:
            reasons.append("underpowered_lt7_donors")
        else:
            reasons.append("included")
    return reasons


old_ctrl = np.array([d in OLD_CONTROLS for d in donors])
panel_tot = m_use.sum(2)  # total panel counts per cell in the used cell set


# ---------------------------------------------------------------- snRNA lineage (donor-deduplicated WS3A)
def snrna_top(gs):
    """Donor-deduplicated snRNA lineage mean CPM (excl_GSE136103, prep_eligible, binary_eligible, lib >= 1e5).

    Uses the exact donor-deduplicated pseudobulks and WS3A primary donor eligibility cohort
    from NUCLEUS-FIX-20261107, collapsing all sequencing captures/lanes to biological donors
    before CPM calculation and enforcing the PRESPEC lib >= 1e5 filter. Also reads WS3A's
    fig4c_gene_specificity_16lineages.tsv.gz as an explicit cross-check.
    """
    dm = pd.read_csv(SC_DONOR_META, sep="\t")
    lm = pd.read_csv(SC_LIB_META, sep="\t").set_index("sample")
    elig_donors = set(
        dm[
            (~to_bool(dm["exclude"]))
            & to_bool(dm["binary_eligible"])
            & dm["condition_harmonized"].isin(["Healthy", "MASLD"])
        ]["donor"].astype(str)
    )
    prep_ok_caps = set(lm.index[to_bool(lm["prep_eligible"])])

    means = {}
    donor_counts = []
    for L in SN_LINEAGES:
        df = pd.read_csv(PB.format(L), index_col=0)
        if df.index.duplicated().any():
            df = df.groupby(level=0, sort=False).sum()
        keep_caps = [
            c for c in df.columns
            if (c in prep_ok_caps) and (str(lm.loc[c, "donor"]) in elig_donors)
        ]
        if not keep_caps:
            means[L] = pd.Series(np.nan, index=gs)
            donor_counts.append((L, 0, 0))
            continue
        d_map = lm.loc[keep_caps, "donor"].astype(str)
        # Collapse captures/lanes to 1 column per biological donor
        d_mat = df.loc[:, keep_caps].T.groupby(d_map).sum().T
        lib = d_mat.sum(0)
        k_donors = lib.index[lib >= 1e5]
        donor_counts.append((L, int(d_mat.shape[1]), int(len(k_donors))))
        if len(k_donors) == 0:
            means[L] = pd.Series(np.nan, index=gs)
            continue
        cpm = d_mat.loc[:, k_donors].div(lib[k_donors], axis=1) * 1e6
        means[L] = cpm.reindex(gs).mean(1)

    cpm_df = pd.DataFrame(means)
    top = cpm_df.idxmax(1)
    shr = cpm_df.max(1) / cpm_df.sum(1)
    dcount_df = pd.DataFrame(
        donor_counts,
        columns=["lineage", "n_elig_donors_any_umi", "n_elig_donors_lib_ge_1e5"],
    )

    spec = pd.read_csv(SC_FIG4C_SPEC, sep="\t").set_index("gene")
    spec_sub = spec.reindex(gs)[["top1_lineage", "top1_z", "tau_16lineages"]]
    return cpm_df, top, shr, dcount_df, spec_sub


# ---------------------------------------------------------------- per-gene results
def gene_rows(gs, family):
    out = []
    for g in gs:
        j = gi[g]
        contrib = share * corr[:, :, j]  # donors x types
        c_ctrl = contrib[ctrl].mean(0)
        frac_ctrl = c_ctrl / c_ctrl.sum() if c_ctrl.sum() > 0 else np.full(nK, np.nan)
        c_all = contrib.mean(0)
        frac_all = c_all / c_all.sum() if c_all.sum() > 0 else np.full(nK, np.nan)
        k_t = int(np.nanargmax(c_ctrl)) if c_ctrl.sum() > 0 else 0
        pc_ctrl = corr[ctrl, :, j].mean(0)
        pc_all = corr[:, :, j].mean(0)

        # Enforce permutation-invariant made_all (spillover >= 3x + donor sufficiency) BEFORE carrier selection
        made_all = np.zeros(nK, dtype=bool)
        made_all_donors = np.zeros(nK, dtype=int)
        for b_idx, k_name in enumerate(TYPES):
            if k_name == HEP:
                made_all[b_idx] = True
                made_all_donors[b_idx] = nD
            else:
                present_ok = int((n_used_cells[:, b_idx] > 0).sum()) >= MIN_PRESENT_DONORS
                mu_all = float(m_use[:, b_idx, j].mean())
                ex_all = float(expected[:, b_idx, j].mean())
                n_made_d = int(
                    ((m_use[:, b_idx, j] > 0.0) & (m_use[:, b_idx, j] >= 3.0 * expected[:, b_idx, j])).sum()
                )
                made_all_donors[b_idx] = n_made_d
                made_all[b_idx] = bool(
                    present_ok and (mu_all > 0.0) and (mu_all >= 3.0 * ex_all) and (n_made_d >= MIN_ARM_DONORS)
                )

        # Permutation-invariant 18-donor grand mean ranking with made_all enforced first
        max_gm = max(float(pc_all.max()), 0.0) + 1.0
        rank_score = np.where(made_all, pc_all, pc_all - max_gm)
        order = np.argsort(-rank_score, kind="stable")
        b = int(order[0])
        b2 = int(order[1])
        carrier = TYPES[b]
        assert made_all[b], f"carrier {carrier} for {g} failed made_all"

        # gate-3 detection rule (controls, plus 18-donor permutation-invariant count)
        det_ok_ctrl = det_use[ctrl, b, j] >= 0.05
        det_ok_all = det_use[:, b, j] >= 0.05
        if carrier != HEP:
            det_ok_ctrl &= m_use[ctrl, b, j] >= 3.0 * expected[ctrl, b, j]
            det_ok_all &= m_use[:, b, j] >= 3.0 * expected[:, b, j]

        w = test(corr[:, b, j])
        w_panel = test(corr[:, b, j] / panel_tot[:, b])
        w_all = test(corr_allcells[:, b, j])
        s = test(share[:, b])
        t_cell = test(tot_counts[:, j] / N_d)
        t_area = test(tot_counts[:, j] / np.array([area[d][0] for d in donors]))
        t_occ = test(tot_counts[:, j] / np.array([area[d][1] for d in donors]))
        age = lr(corr[mash, b, j], corr[old_ctrl, b, j])

        row = dict(
            family=family,
            gene=g,
            carrier=carrier,
            carrier_selection_rule="grand_mean_18donors_made_all_prefilter",
            carrier_frac_ctrl=float(frac_ctrl[b]),
            carrier_frac_all=float(frac_all[b]),
            top_contributor=TYPES[k_t],
            top_contributor_frac_ctrl=float(frac_ctrl[k_t]),
            made_by_carrier=bool(made_all[b]),
            made_all_donors_carrier=int(made_all_donors[b]),
            frac_ctrl_by_type=";".join(f"{TYPES[x]}:{frac_ctrl[x]:.3f}" for x in np.argsort(-frac_ctrl)[:4]),
            percell_ctrl_by_type=";".join(f"{TYPES[x]}:{pc_ctrl[x]:.3f}" for x in np.argsort(-pc_ctrl)[:4]),
            percell_all_by_type=";".join(f"{TYPES[x]}:{pc_all[x]:.3f}" for x in np.argsort(-pc_all)[:4]),
            carrier_cell_set="interior" if use_interior[carrier] else "all",
            all_mean_percell_carrier=float(corr[:, b, j].mean()),
            ctrl_mean_percell_carrier=float(corr[ctrl, b, j].mean()),
            mash_mean_percell_carrier=float(corr[mash, b, j].mean()),
            detect_ctrl_donors_pass=int(det_ok_ctrl.sum()),
            detect_all_donors_pass=int(det_ok_all.sum()),
            detect_frac_ctrl_median=float(np.median(det_use[ctrl, b, j])),
            detect_frac_all_median=float(np.median(det_use[:, b, j])),
            label_defining=g in LABEL_DEFINING,
        )
        row.update({f"within_{k}": v for k, v in w.items()})
        row.update({"within_panelshare_log2": w_panel["log2_ratio"], "within_allcells_log2": w_all["log2_ratio"]})
        row.update({f"share_{k}": v for k, v in s.items() if k in ("log2_ratio", "p_perm", "ci_lo", "ci_hi")})
        row.update({f"tissue_percell_{k}": v for k, v in t_cell.items() if k in ("log2_ratio", "p_perm", "ci_lo", "ci_hi")})
        row.update({
            "tissue_permm2_hull_log2": t_area["log2_ratio"],
            "tissue_permm2_occ_log2": t_occ["log2_ratio"],
            "within_vs_old_controls_log2": age,
        })
        row["within_units_agree"] = bool(np.sign(w["log2_ratio"]) == np.sign(w_panel["log2_ratio"]))
        row["tissue_units_agree"] = bool(np.sign(t_cell["log2_ratio"]) == np.sign(t_area["log2_ratio"]))

        w2 = test(corr[:, b2, j])
        row.update(dict(
            second_type=TYPES[b2],
            made_by_second_type=bool(made_all[b2]),
            second_type_frac_ctrl=float(frac_ctrl[b2]),
            second_within_log2=w2["log2_ratio"],
            second_within_p=w2["p_perm"],
            second_within_ci_lo=w2["ci_lo"],
            second_within_ci_hi=w2["ci_hi"],
        ))
        out.append(row)
    return pd.DataFrame(out)


b5 = pd.read_csv(FIG5B, sep="\t").set_index("gene")
res5d = gene_rows(D5, "5D")
res5d["bulk_logFC"] = b5.loc[res5d.gene, "bulk_logFC"].values
res5d["protein_logFC"] = b5.loc[res5d.gene, "protein_logFC"].values
assert (b5.loc[D5, "class"] == "confirmed").all()
assert (np.sign(res5d.bulk_logFC) == np.sign(res5d.protein_logFC)).all()
res5d["tissue_sign_agrees"] = np.sign(res5d.tissue_percell_log2_ratio) == np.sign(res5d.bulk_logFC)
res5d["in_primary_within_holm"] = (~res5d.label_defining) & (res5d.detect_ctrl_donors_pass >= MIN_ARM_DONORS)
res5d["primary_holm_exclusion_reason"] = primary_exclusion_reasons(res5d.label_defining, res5d.detect_ctrl_donors_pass)
res5d["within_holm_primary"] = holm_subset(res5d.within_p_perm, res5d.in_primary_within_holm)
res5d["within_holm_all_sens"] = holm(res5d.within_p_perm)
res5d["primary_holm_family_n"] = int(res5d.in_primary_within_holm.sum())
res5d["all_rows_family_n"] = len(res5d)

res5e = gene_rows(E5, "5E")
res5e["in_primary_within_holm"] = (~res5e.label_defining) & (res5e.detect_ctrl_donors_pass >= MIN_ARM_DONORS)
res5e["primary_holm_exclusion_reason"] = primary_exclusion_reasons(res5e.label_defining, res5e.detect_ctrl_donors_pass)
res5e["within_holm_primary"] = holm_subset(res5e.within_p_perm, res5e.in_primary_within_holm)
res5e["within_holm_all_sens"] = holm(res5e.within_p_perm)
res5e["primary_holm_family_n"] = int(res5e.in_primary_within_holm.sum())
res5e["all_rows_family_n"] = len(res5e)

cpm, top, shr, sn_dcounts, spec_sub = snrna_top(D5 + E5)
for r in (res5d, res5e):
    r["snrna_top_lineage"] = top.reindex(r.gene).values
    r["snrna_top_share"] = shr.reindex(r.gene).values
    r["snrna_matches_carrier"] = [
        ct in SN_TO_XEN.get(t, set()) if isinstance(t, str) else np.nan
        for ct, t in zip(r.carrier, r.snrna_top_lineage)
    ]
    r["snrna_fig4c_top1_lineage"] = spec_sub.reindex(r.gene)["top1_lineage"].values
    r["snrna_fig4c_top1_z"] = spec_sub.reindex(r.gene)["top1_z"].values
    r["snrna_fig4c_tau_16lineages"] = spec_sub.reindex(r.gene)["tau_16lineages"].values
    r["snrna_fig4c_matches_carrier"] = [
        ct in SN_TO_XEN.get(t, set()) if isinstance(t, str) else np.nan
        for ct, t in zip(r.carrier, r.snrna_fig4c_top1_lineage)
    ]

write(cpm.reset_index().rename(columns={"index": "gene"}), "snrna_lineage_mean_cpm.tsv")
write(sn_dcounts, "snrna_lineage_donor_counts.tsv")
write(res5d, "results_5d.tsv")
write(res5e, "results_5e.tsv")

# Gate 2
n_agree = int(res5d.tissue_sign_agrees.sum())
log(f"Gate 2: tissue sign agrees for {n_agree} of {len(res5d)} 5D genes (need >= 10)")

# Supplementary: on-panel canonical calls
cat = pd.read_csv(CATALOG, sep="\t", usecols=["symbol", "primary_evidence_class", "rna_log2fc", "rna_canonical_call"])
cat = cat[cat.rna_canonical_call.astype(str).str.lower() == "true"].drop_duplicates("symbol")
calls = sorted(set(cat.symbol) & set(genes))
log(f"on-panel canonical calls: {len(calls)}")
rc = gene_rows(calls, "calls")
rc["bulk_log2fc"] = cat.set_index("symbol").loc[rc.gene, "rna_log2fc"].values
rc["tissue_sign_agrees"] = np.sign(rc.tissue_percell_log2_ratio) == np.sign(rc.bulk_log2fc)
rc["in_primary_within_bh"] = (~rc.label_defining) & (rc.detect_ctrl_donors_pass >= MIN_ARM_DONORS)
rc["primary_bh_exclusion_reason"] = primary_exclusion_reasons(rc.label_defining, rc.detect_ctrl_donors_pass)
rc["within_bh_primary"] = bh_subset(rc.within_p_perm, rc.in_primary_within_bh)
rc["within_bh_all_sens"] = bh(rc.within_p_perm)
rc["primary_bh_family_n"] = int(rc.in_primary_within_bh.sum())
rc["all_rows_family_n"] = len(rc)
write(rc, "results_calls_supplementary.tsv")
log(f"canonical calls: tissue sign agrees {int(rc.tissue_sign_agrees.sum())} of {len(rc)}")

# Gate 3
n_det = int((res5e.detect_ctrl_donors_pass >= MIN_ARM_DONORS).sum())
log(f"Gate 3: {n_det} of 12 colocalized genes pass detection (need >= 3)")

gates = pd.DataFrame([
    ("gate1_donors_per_arm", f"{int(mash.sum())} MASH / {int(ctrl.sum())} Control", "pass"),
    ("gate2_sign_agree_5d", f"{n_agree}/14", "pass" if n_agree >= 10 else "FAIL"),
    ("gate2_sign_agree_calls", f"{int(rc.tissue_sign_agrees.sum())}/{len(rc)}", "report"),
    ("gate3_detect_5e", f"{n_det}/12", "pass" if n_det >= 3 else "FAIL: fold into 5D"),
    ("primary_holm_5d_n", f"{int(res5d.in_primary_within_holm.sum())}/14", "4 label_defining excluded"),
    ("primary_holm_5e_n", f"{int(res5e.in_primary_within_holm.sum())}/12", "SPI1 <7 donors excluded"),
], columns=["gate", "value", "status"])
write(gates, "gates.tsv")
log(gates.to_string())

# Positive checks
chk = []
for k in ("Cholangiocyte", HEP):
    b = TYPES.index(k)
    chk.append((
        f"share_{k}",
        share[ctrl, b].mean(),
        share[mash, b].mean(),
        "higher in MASH" if k == "Cholangiocyte" else "lower in MASH",
        (share[mash, b].mean() > share[ctrl, b].mean()) == (k == "Cholangiocyte"),
    ))
for g, want in (("LYVE1", "LSEC"), ("ACTA2", "VSMC"), ("MZB1", "B")):
    got = res5d.set_index("gene").loc[g, "carrier"]
    chk.append((f"carrier_{g}", got, want, "carrier type", got == want))
write(pd.DataFrame(chk, columns=["check", "observed", "expected_or_mash", "rule", "pass"]), "positive_checks.tsv")

# Cell shares by donor
write(
    pd.DataFrame(share, index=donors, columns=TYPES).assign(condition=np.where(mash, "MASH", "Control")),
    "cell_share_by_donor.tsv",
    index=True,
)

# 5A donor rule and cell-level export
up_gene = res5d[res5d.gene.isin(D5_UP)].sort_values("tissue_percell_log2_ratio", ascending=False).gene.iloc[0]
hep_depth = panel_tot[:, h]
pick = []
for arm in (False, True):
    idx = np.where(mash == arm)[0]
    lo, hi = np.percentile(hep_depth[idx], [25, 75])
    mid = idx[(hep_depth[idx] >= lo) & (hep_depth[idx] <= hi)]
    med = np.median(share[idx, h])
    best = mid[np.argmin(np.abs(share[mid, h] - med))]
    pick.append((donors[best], "MASH" if arm else "Control", share[best, h], med, hep_depth[best]))
pick = pd.DataFrame(
    pick,
    columns=["donor_id", "condition", "hepatocyte_share", "arm_median_share", "counts_per_interior_hepatocyte"],
)
pick["genes_shown"] = f"CYP3A4,{up_gene}"
write(pick, "panel_a_donors.tsv")
log(pick.to_string())
cells = []
for d in pick.donor_id:
    idx = np.where(meta.donor_id.values == d)[0]
    sub = X[idx][:, [gi[g] for g in D5 + E5 + ["KRT7"]]].toarray().astype(int)
    df = pd.DataFrame(sub, columns=D5 + E5 + ["KRT7"])
    df.insert(0, "interior", interior[idx])
    df.insert(0, "type", ctype[idx])
    df.insert(0, "Y", xy_all[idx, 1])
    df.insert(0, "X", xy_all[idx, 0])
    df.insert(0, "donor_id", d)
    cells.append(df)
write(pd.concat(cells), "panel_a_cells.tsv.gz")
log("x1_xenium_v2_20261118 base pass done; running selection-aware exact inference")

# Exact permutation + fixed-carrier conditional & unconditional bootstrap + HPX sensitivity
_selection_path = Path(__file__).with_name("x1_selection_aware_20261118.py")
_selection_spec = importlib.util.spec_from_file_location("x1_selection_aware_20261118", _selection_path)
_selection_module = importlib.util.module_from_spec(_selection_spec)
_selection_spec.loader.exec_module(_selection_module)
_selection_module.correct_results(Path(OUT), n_boot=N_BOOT, seed=SEED)

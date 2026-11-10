#!/usr/bin/env python3
"""X1: Fig 5 genes in intact liver (SCP3359 Xenium, 9 control vs 9 MASH-cirrhosis donors).

Rules are fixed in ../PRESPEC.md, with the two changes in ../DEVIATIONS.md (v2). Unit = donor. Usage: x1_xenium.py <out_dir>
"""
import itertools
import sys

import numpy as np
import pandas as pd
from scipy import io
from scipy.spatial import ConvexHull, cKDTree

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
SRC = f"{ROOT}/Analysis/Hypotheses/V7-20261021T235717Z/downloads_untrusted/scp3359/SCP3359"
EXPR = f"{SRC}/expression/68f4186f4e6da7a7e1b0edbe"
CATALOG = (f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/"
           "catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.tsv.gz")
FIG5B = f"{ROOT}/figures/main/fig5_molecular_context/panels/data/fig5b_protein_triage_genes.tsv"
PB = f"{ROOT}/Analysis/SingleCell/results_gpu_v2/pseudobulk/{{}}_pseudobulk.csv"
OUT = sys.argv[1]

SEED = 20261027
N_BOOT = 10_000
D5_DOWN = ["CYP3A4", "TAT", "GLS2", "MET", "LYVE1"]
D5_UP = ["FASN", "KRT8", "KRT18", "MZB1", "PLIN2", "PLIN4", "ACTA2", "ANXA5", "ERBB2"]
D5 = D5_DOWN + D5_UP
E5 = ["THRB", "SPI1", "KLRC1", "TNFSF10", "ADH4", "ADH1C", "ANPEP", "TM4SF4", "SCGN", "MEST", "MTOR", "CFHR1"]
REF = ["HPX", "KNG1", "FBP1"]
LABEL_DEFINING = {"KRT8", "KRT18", "MZB1", "LYVE1"}
HEP = "Hepatocyte"
TYPES = ["Hepatocyte", "Cholangiocyte", "LSEC", "Endothelial", "Stellate", "VSMC", "Kupffer", "Macrophage",
         "B", "T.CD4", "T.CD8", "T.Innate"]
# v2: the source files DNASE1L3+FCN2+LYVE1-high, PDGFRB-low cells under broad "Stellate"; they are sinusoidal
# endothelium (label_check_refined_means.tsv). VSMC are split out of "Stellate" as well.
REFINED_TO_V2 = {"Stellate.DNASE1L3+CD36+": "LSEC", "VSMC.MYH11+": "VSMC"}
SN_LINEAGES = ["Hepatocytes", "Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Macrophages",
               "Mono+mono_derived_cells", "T_cells", "B_cells", "Plasma_cells", "Resident_NK",
               "Circulating_NK_NKT"]
SN_TO_XEN = {"Hepatocytes": {"Hepatocyte"}, "Cholangiocytes": {"Cholangiocyte"},
             "Endothelial_cells": {"Endothelial", "LSEC"}, "Fibroblasts": {"Stellate", "VSMC"},
             "Macrophages": {"Kupffer", "Macrophage"}, "Mono+mono_derived_cells": {"Macrophage", "Kupffer"},
             "T_cells": {"T.CD4", "T.CD8"}, "Resident_NK": {"T.Innate"},
             "Circulating_NK_NKT": {"T.Innate"}, "B_cells": {"B"}, "Plasma_cells": {"B"}}
OLD_CONTROLS = ["11-00C", "11-025"]
INTERIOR_NONHEP_UM = 15.0
MIN_INTERIOR_CELLS = 30
MIN_INTERIOR_HEP = 500


def log(msg):
    print(msg, flush=True)


def write(df, name, index=False):
    df.to_csv(f"{OUT}/{name}", sep="\t", index=index)


# ---------------------------------------------------------------- load
genes = [l.strip() for l in open(f"{EXPR}/counts_features.tsv") if l.strip()]
barcodes = [l.strip() for l in open(f"{EXPR}/counts_barcodes.tsv") if l.strip()]
gi = {g: i for i, g in enumerate(genes)}
for g in D5 + E5 + REF:
    assert g in gi, g
log(f"genes {len(genes)}, cells {len(barcodes)}")

meta = pd.read_csv(f"{SRC}/cluster/coords.txt", sep="\t", skiprows=[1],
                   usecols=["NAME", "X", "Y", "broad_cell_type", "refined_cell_type", "condition",
                            "age_at_sample_collection", "donor_id"],
                   dtype={"NAME": str, "broad_cell_type": str, "condition": str, "donor_id": str})
meta = meta.set_index("NAME")
assert meta.index.is_unique
meta = meta.loc[barcodes]
meta["broad_cell_type"] = [REFINED_TO_V2.get(r, b) for r, b in zip(meta.refined_cell_type, meta.broad_cell_type)]
assert set(meta.broad_cell_type) == set(TYPES), set(meta.broad_cell_type)

X = io.mmread(f"{EXPR}/counts_matrix.mtx").T.tocsr().astype(np.float64)  # cells x genes
assert X.shape == (len(barcodes), len(genes))
log(f"matrix loaded, nnz {X.nnz}")

donor_tab = (meta.groupby("donor_id")
             .agg(condition=("condition", "first"), age=("age_at_sample_collection", "first"),
                  n_cells=("X", "size")).reset_index())
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
g1 = pd.DataFrame({"donor_id": donors, "condition": np.where(is_mash, "MASH", "Control"),
                   "n_cells": N_d, "n_hepatocytes": n_all[:, h], "n_interior_hepatocytes": n_int[:, h],
                   "hepatocyte_share": share[:, h],
                   "area_hull_mm2": [area[d][0] for d in donors],
                   "area_occupancy50_mm2": [area[d][1] for d in donors]})
g1["pass"] = g1.n_interior_hepatocytes >= MIN_INTERIOR_HEP
write(g1, "gate1_donors.tsv")
keep = g1["pass"].to_numpy()
log(f"Gate 1: {keep.sum()} of {nD} donors pass "
    f"({(keep & is_mash).sum()} MASH, {(keep & ~is_mash).sum()} Control)")
if (keep & is_mash).sum() < 7 or (keep & ~is_mash).sum() < 7:
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

# cell set per type
use_interior = {}
for b, k in enumerate(TYPES):
    if k == HEP:
        use_interior[k] = True
        continue
    ok = n_int[:, b] >= MIN_INTERIOR_CELLS
    use_interior[k] = bool((ok & mash).sum() >= 7 and (ok & ctrl).sum() >= 7)
cellset = pd.DataFrame({"type": TYPES, "cell_set": ["interior" if use_interior[k] else "all" for k in TYPES],
                        "control_median_interior_n": [np.median(n_int[ctrl, b]) for b in range(nK)],
                        "mash_median_interior_n": [np.median(n_int[mash, b]) for b in range(nK)],
                        "control_median_all_n": [np.median(n_all[ctrl, b]) for b in range(nK)],
                        "mash_median_all_n": [np.median(n_all[mash, b]) for b in range(nK)]})
write(cellset, "cell_set_by_type.tsv")
log(cellset.to_string())


def corrected_levels(m_set, hep_ref):
    """Spillover-corrected per-cell levels; returns corr, expected, leak."""
    ridx = [gi[r] for r in REF]
    leak = np.zeros((nD, nK))
    for b, k in enumerate(TYPES):
        if k == HEP:
            continue
        ratios = m_set[:, b][:, ridx] / hep_ref[:, ridx]
        leak[:, b] = np.nanmedian(ratios, axis=1)
    expected = leak[:, :, None] * hep_ref[:, None, :]
    expected[:, h, :] = 0.0
    corr = np.clip(m_set - expected, 0, None)
    return corr, expected, leak


hep_ref = m_int[:, h, :]
m_use = np.where(np.array([use_interior[k] for k in TYPES])[None, :, None], m_int, m_all)
m_use = np.nan_to_num(m_use, nan=0.0)
det_use = np.nan_to_num(np.where(np.array([use_interior[k] for k in TYPES])[None, :, None], det_int, det_all))
corr, expected, leak = corrected_levels(m_use, hep_ref)
m_all_h = m_all.copy()
m_all_h[:, h] = hep_ref
corr_allcells, expected_allcells, _ = corrected_levels(m_all_h, hep_ref)
write(pd.DataFrame(leak, index=donors, columns=TYPES), "spillover_leak_by_donor_type.tsv", index=True)

# long table for later rendering
rows = []
for a, d in enumerate(donors):
    for b, k in enumerate(TYPES):
        for g in sorted(set(D5 + E5 + REF + ["KRT7"])):
            j = gi[g]
            rows.append((d, "MASH" if is_mash[a] else "Control", k, g, n_all[a, b], n_int[a, b],
                         share[a, b], m_all[a, b, j], m_use[a, b, j], expected[a, b, j], corr[a, b, j],
                         det_use[a, b, j]))
write(pd.DataFrame(rows, columns=["donor_id", "condition", "type", "gene", "n_cells", "n_interior",
                                  "share", "mean_all", "mean_used", "expected_spillover", "mean_corrected",
                                  "detect_frac_used"]),
      "donor_type_gene_levels.tsv.gz")

# ---------------------------------------------------------------- statistics helpers
A = np.array([np.isin(np.arange(nD), c) for c in itertools.combinations(range(nD), int(mash.sum()))])
log(f"permutation splits: {A.shape[0]}")
rng = np.random.default_rng(SEED)
mi, ci = np.where(mash)[0], np.where(ctrl)[0]
BM = rng.integers(0, len(mi), size=(N_BOOT, len(mi)))
BC = rng.integers(0, len(ci), size=(N_BOOT, len(ci)))


def lr(v_m, v_c):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log2(np.mean(v_m) / np.mean(v_c))


def test(v):
    """log2 ratio of arm means, exact permutation p, bootstrap CI, leave-one-out range."""
    v = np.asarray(v, float)
    obs = lr(v[mash], v[ctrl])
    if not np.isfinite(obs):
        return dict(log2_ratio=obs, p_perm=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                    lodo_min=np.nan, lodo_max=np.nan, lodo_sign_flips=np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        pm = (A @ v) / A.sum(1)
        pc = ((~A) @ v) / (~A).sum(1)
        perm = np.log2(pm / pc)
    p = np.mean(np.abs(perm) >= abs(obs) - 1e-12)
    vm, vc = v[mi], v[ci]
    with np.errstate(divide="ignore", invalid="ignore"):
        boot = np.log2(vm[BM].mean(1) / vc[BC].mean(1))
    boot = boot[np.isfinite(boot)]
    lo, hi = np.percentile(boot, [2.5, 97.5]) if len(boot) else (np.nan, np.nan)
    lodo = []
    for drop in range(nD):
        k = np.arange(nD) != drop
        lodo.append(lr(v[k & mash], v[k & ctrl]))
    lodo = np.array(lodo)
    return dict(log2_ratio=obs, p_perm=p, ci_lo=lo, ci_hi=hi, lodo_min=np.nanmin(lodo),
                lodo_max=np.nanmax(lodo), lodo_sign_flips=int(np.sum(np.sign(lodo) != np.sign(obs))))


def holm(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    order = np.argsort(q)
    adj = np.maximum.accumulate((len(q) - np.arange(len(q))) * q[order]).clip(max=1)
    res = np.empty_like(q)
    res[order] = adj
    out[ok] = res
    return out


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    n = len(q)
    order = np.argsort(q)
    adj = np.minimum.accumulate((q[order] * n / np.arange(1, n + 1))[::-1])[::-1].clip(max=1)
    res = np.empty_like(q)
    res[order] = adj
    out[ok] = res
    return out


old_ctrl = np.array([d in OLD_CONTROLS for d in donors])
panel_tot = m_use.sum(2)  # total panel counts per cell in the used cell set


# ---------------------------------------------------------------- snRNA lineage
def snrna_top(gs):
    means = {}
    for L in SN_LINEAGES:
        m = pd.read_csv(PB.format(L), index_col=0)
        lib = m.sum(0)
        k = lib.index[lib >= 1e5]
        cpm = m.loc[:, k].div(lib[k], axis=1) * 1e6
        means[L] = cpm.reindex(gs).mean(1)
    cpm = pd.DataFrame(means)
    top = cpm.idxmax(1)
    shr = cpm.max(1) / cpm.sum(1)
    return cpm, top, shr


# ---------------------------------------------------------------- per-gene results
def gene_rows(gs, family):
    out = []
    for g in gs:
        j = gi[g]
        contrib = share * corr[:, :, j]  # donors x types
        c_ctrl = contrib[ctrl].mean(0)
        frac_ctrl = c_ctrl / c_ctrl.sum() if c_ctrl.sum() > 0 else np.full(nK, np.nan)
        k_t = int(np.nanargmax(c_ctrl)) if c_ctrl.sum() > 0 else 0
        pc_ctrl = corr[ctrl, :, j].mean(0)
        k_c = int(np.argmax(pc_ctrl))
        made = {}
        for b, k in enumerate(TYPES):
            made[k] = True if k == HEP else bool(m_use[ctrl, b, j].mean() >= 3 * expected[ctrl, b, j].mean())
        carrier = TYPES[k_c]  # v2: carrier = highest spillover-corrected level per cell in controls
        b = k_c
        # gate-3 detection rule (controls)
        det_ok = (det_use[ctrl, b, j] >= 0.05)
        if carrier != HEP:
            det_ok &= m_use[ctrl, b, j] >= 3 * expected[ctrl, b, j]
        w = test(corr[:, b, j])
        w_panel = test(corr[:, b, j] / panel_tot[:, b])
        w_all = test(corr_allcells[:, b, j])
        s = test(share[:, b])
        t_cell = test(tot_counts[:, j] / N_d)
        t_area = test(tot_counts[:, j] / np.array([area[d][0] for d in donors]))
        t_occ = test(tot_counts[:, j] / np.array([area[d][1] for d in donors]))
        age = lr(corr[mash, b, j], corr[old_ctrl, b, j])
        row = dict(family=family, gene=g, carrier=carrier, carrier_frac_ctrl=frac_ctrl[b],
                   top_contributor=TYPES[k_t], top_contributor_frac_ctrl=frac_ctrl[k_t],
                   made_by_carrier=made[carrier],
                   frac_ctrl_by_type=";".join(f"{TYPES[x]}:{frac_ctrl[x]:.3f}" for x in np.argsort(-frac_ctrl)[:4]),
                   percell_ctrl_by_type=";".join(f"{TYPES[x]}:{pc_ctrl[x]:.3f}" for x in np.argsort(-pc_ctrl)[:4]),
                   carrier_cell_set="interior" if use_interior[carrier] else "all",
                   ctrl_mean_percell_carrier=corr[ctrl, b, j].mean(), mash_mean_percell_carrier=corr[mash, b, j].mean(),
                   detect_ctrl_donors_pass=int(det_ok.sum()),
                   detect_frac_ctrl_median=float(np.median(det_use[ctrl, b, j])),
                   label_defining=g in LABEL_DEFINING)
        row.update({f"within_{k}": v for k, v in w.items()})
        row.update({"within_panelshare_log2": w_panel["log2_ratio"], "within_allcells_log2": w_all["log2_ratio"]})
        row.update({f"share_{k}": v for k, v in s.items() if k in ("log2_ratio", "p_perm", "ci_lo", "ci_hi")})
        row.update({f"tissue_percell_{k}": v for k, v in t_cell.items() if k in ("log2_ratio", "p_perm", "ci_lo", "ci_hi")})
        row.update({"tissue_permm2_hull_log2": t_area["log2_ratio"], "tissue_permm2_occ_log2": t_occ["log2_ratio"],
                    "within_vs_old_controls_log2": age})
        row["within_units_agree"] = bool(np.sign(w["log2_ratio"]) == np.sign(w_panel["log2_ratio"]))
        row["tissue_units_agree"] = bool(np.sign(t_cell["log2_ratio"]) == np.sign(t_area["log2_ratio"]))
        # second type: next-largest control contribution, for genes made by two types (e.g. PLIN2)
        order = [x for x in np.argsort(-pc_ctrl) if x != b]
        b2 = int(order[0])
        w2 = test(corr[:, b2, j])
        row.update(dict(second_type=TYPES[b2], second_type_frac_ctrl=frac_ctrl[b2],
                        second_within_log2=w2["log2_ratio"], second_within_p=w2["p_perm"],
                        second_within_ci_lo=w2["ci_lo"], second_within_ci_hi=w2["ci_hi"]))
        out.append(row)
    return pd.DataFrame(out)


b5 = pd.read_csv(FIG5B, sep="\t").set_index("gene")
res5d = gene_rows(D5, "5D")
res5d["bulk_logFC"] = b5.loc[res5d.gene, "bulk_logFC"].values
res5d["protein_logFC"] = b5.loc[res5d.gene, "protein_logFC"].values
assert (b5.loc[D5, "class"] == "confirmed").all()
assert (np.sign(res5d.bulk_logFC) == np.sign(res5d.protein_logFC)).all()
res5d["tissue_sign_agrees"] = np.sign(res5d.tissue_percell_log2_ratio) == np.sign(res5d.bulk_logFC)
res5d["within_holm14"] = holm(res5d.within_p_perm)
res5e = gene_rows(E5, "5E")
res5e["within_holm12"] = holm(res5e.within_p_perm)

cpm, top, shr = snrna_top(D5 + E5)
for r in (res5d, res5e):
    r["snrna_top_lineage"] = top.reindex(r.gene).values
    r["snrna_top_share"] = shr.reindex(r.gene).values
    r["snrna_matches_carrier"] = [ct in SN_TO_XEN.get(t, set()) if isinstance(t, str) else np.nan
                                  for ct, t in zip(r.carrier, r.snrna_top_lineage)]
write(cpm.reset_index(), "snrna_lineage_mean_cpm.tsv")
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
rc["within_bh"] = bh(rc.within_p_perm)
write(rc, "results_calls_supplementary.tsv")
log(f"canonical calls: tissue sign agrees {int(rc.tissue_sign_agrees.sum())} of {len(rc)}")

# Gate 3
n_det = int((res5e.detect_ctrl_donors_pass >= 7).sum())
log(f"Gate 3: {n_det} of 12 colocalized genes pass detection (need >= 3)")

gates = pd.DataFrame([
    ("gate1_donors_per_arm", f"{int(mash.sum())} MASH / {int(ctrl.sum())} Control", "pass"),
    ("gate2_sign_agree_5d", f"{n_agree}/14", "pass" if n_agree >= 10 else "FAIL"),
    ("gate2_sign_agree_calls", f"{int(rc.tissue_sign_agrees.sum())}/{len(rc)}", "report"),
    ("gate3_detect_5e", f"{n_det}/12", "pass" if n_det >= 3 else "FAIL: fold into 5D"),
], columns=["gate", "value", "status"])
write(gates, "gates.tsv")
log(gates.to_string())

# Positive checks
chk = []
for k in ("Cholangiocyte", HEP):
    b = TYPES.index(k)
    chk.append((f"share_{k}", share[ctrl, b].mean(), share[mash, b].mean(),
                "higher in MASH" if k == "Cholangiocyte" else "lower in MASH",
                (share[mash, b].mean() > share[ctrl, b].mean()) == (k == "Cholangiocyte")))
for g, want in (("LYVE1", "LSEC"), ("ACTA2", "VSMC"), ("MZB1", "B")):
    got = res5d.set_index("gene").loc[g, "carrier"]
    chk.append((f"carrier_{g}", got, want, "carrier type", got == want))
write(pd.DataFrame(chk, columns=["check", "observed", "expected_or_mash", "rule", "pass"]), "positive_checks.tsv")

# Cell shares by donor
write(pd.DataFrame(share, index=donors, columns=TYPES).assign(condition=np.where(mash, "MASH", "Control")),
      "cell_share_by_donor.tsv", index=True)

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
pick = pd.DataFrame(pick, columns=["donor_id", "condition", "hepatocyte_share", "arm_median_share",
                                   "counts_per_interior_hepatocyte"])
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
log("done")

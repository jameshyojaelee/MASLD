#!/usr/bin/env python3
"""Gate L01-04 pilot (job 3051). Feasibility only; not a result.

Within-hepatocyte, nucleus-level peak-gene association in GSE296875, SCENT-style:
Poisson GLM of gene count on binarized peak accessibility, covariates log nUMI,
fraction mitochondrial, donor fixed effects (donor is nested in well, so donor absorbs well).
Inference: Poisson model SE and robust sandwich (HC0) SE. The sandwich is the analytic
counterpart of SCENT's nonparametric nucleus bootstrap; no bootstrap is run here.
RNA: authors' decontX matrix rounded to integers (raw 10x counts were not re-read in this
8-minute job). Seed 20261015.
Outputs go to the job directory.
"""
import gzip, os, sys, time, json
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy import stats

T0 = time.time()
DEADLINE = 380  # seconds of wall time for the fitting stage cut-off
SEED = 20261015
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RECT = f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions/gse296875-crossmodal-recurrence-20260830T193000Z/rectangle"
ROSTER = f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions/gse296875-fragment-membership-21063829/cell_membership.tsv.gz"
SCR = f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-04"
OUT = os.getcwd()
MIN_DET = 0.05  # SCENT's own filter: peak and gene detected in >= 5% of analyzed cells
WINDOW = 500_000
SORT1_PEAK = "chr1-109274474-109275151"


def log(*a):
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)


nuclei = [l.strip() for l in open(f"{RECT}/nuclei.txt")]
n_nuc = len(nuclei)
ro = pd.read_csv(ROSTER, sep="\t", dtype=str).set_index("cell_id")
ro = ro.loc[nuclei]
hep = (ro["lineage_id"] == "hepatocyte").to_numpy()
donor = ro["donor_id"].to_numpy()
well = ro["well_id"].to_numpy()
hep_idx = np.flatnonzero(hep)
n_hep = hep.sum()
log("nuclei", n_nuc, "hepatocyte nuclei", n_hep, "donors", len(set(donor[hep])))

# ---------------- ATAC: binarized detection ----------------
peaks = [l.strip() for l in open(f"{RECT}/atac.features.txt")]
n_peak = len(peaks)
pk = pd.read_csv(f"{RECT}/atac_peaks.tsv", sep="\t").set_index("peak_id").loc[peaks].reset_index()
a_ptr = np.fromfile(f"{RECT}/atac.indptr.i32", dtype=np.int32).astype(np.int64)
a_ind = np.fromfile(f"{RECT}/atac.indices.i32", dtype=np.int32)
a_dat_head = np.fromfile(f"{RECT}/atac.data.i32", dtype=np.int32, count=20_000_000)
log("ATAC explicit zeros in first 2e7 entries:", int((a_dat_head <= 0).sum()))
del a_dat_head
col_hep_entry = np.repeat(hep, np.diff(a_ptr))
peak_det = np.bincount(a_ind[col_hep_entry], minlength=n_peak) / n_hep
del col_hep_entry
atac_nnz_per_nuc = np.diff(a_ptr)
pk["hep_det"] = peak_det
log("hepatocyte peak detection quantiles 5/25/50/75/95%:", np.round(np.quantile(peak_det, [.05, .25, .5, .75, .95]), 4),
    "peaks >= 5%:", int((peak_det >= MIN_DET).sum()))

# ---------------- RNA (decontX, rounded) ----------------
genes = pd.read_csv(f"{SCR}/job_gene_table.tsv", sep="\t")
n_gene = len(genes)
r_ptr = np.fromfile(f"{RECT}/expr_from_rds.indptr.i32", dtype=np.int32).astype(np.int64)
r_ind = np.fromfile(f"{RECT}/expr_from_rds.indices.i32", dtype=np.int32)
r_dat = np.round(np.fromfile(f"{RECT}/expr_from_rds.data.f64", dtype=np.float64)).astype(np.int32)
log("RNA entries", r_dat.size, "rounded-to-zero fraction", round(float((r_dat == 0).mean()), 4))
numi = np.add.reduceat(r_dat.astype(np.int64), r_ptr[:-1])
mt = genes["name_rect"].str.startswith("MT-").to_numpy()
mt_entry = mt[r_ind]
mt_counts = np.add.reduceat(np.where(mt_entry, r_dat, 0).astype(np.int64), r_ptr[:-1])
del mt_entry
pct_mt = mt_counts / np.maximum(numi, 1)
col_hep_entry = np.repeat(hep, np.diff(r_ptr))
pos = col_hep_entry & (r_dat > 0)
gene_det = np.bincount(r_ind[pos], minlength=n_gene) / n_hep
del pos, col_hep_entry
genes["hep_det"] = gene_det
log("hepatocyte nUMI (rounded decontX) median", int(np.median(numi[hep])), "genes >= 5%:", int((gene_det >= MIN_DET).sum()))


def extract_rows(ptr, ind, dat, rows, nuc_sel):
    """Return dense (len(rows) x len(nuc_sel)) int matrix for selected feature rows and nuclei."""
    rows = np.asarray(rows)
    lut = np.full(int(ind.max()) + 1 if ind.size else 1, -1, dtype=np.int32)
    lut[rows] = np.arange(rows.size, dtype=np.int32)
    loc = lut[ind]
    keep = np.flatnonzero(loc >= 0)
    col = np.searchsorted(ptr, keep, side="right") - 1
    colmap = np.full(ptr.size - 1, -1, dtype=np.int64)
    colmap[nuc_sel] = np.arange(nuc_sel.size)
    c = colmap[col]
    ok = c >= 0
    M = np.zeros((rows.size, nuc_sel.size), dtype=np.int32)
    vals = dat[keep[ok]] if dat is not None else 1
    M[loc[keep[ok]], c[ok]] = vals
    return M


# ---------------- choose pairs ----------------
rng = np.random.default_rng(SEED)
pk["center"] = (pk["start_1based"] + pk["end_1based"]) // 2
gok = genes[genes["chrom"].notna() & (genes["hep_det"] >= MIN_DET)].copy()
gok["tss"] = gok["tss"].astype(int)

# Locus pairs
loci = pd.read_csv(f"{SCR}/job_locus_inputs.tsv", sep="\t")
pk_idx = pd.Series(np.arange(n_peak), index=peaks)
loci["peak_row"] = loci["peak_id"].map(pk_idx)
loci["hep_det"] = peak_det[loci["peak_row"].to_numpy()]
loc_pairs = []
loc_rows = []
for _, L in loci.iterrows():
    p = pk.iloc[int(L["peak_row"])]
    allg = genes[(genes["chrom"] == p["chromosome"]) & genes["tss"].notna()]
    allg = allg[(allg["tss"] - p["center"]).abs() <= WINDOW]
    cset = set(str(L["coloc_ensembl"]).split(";"))
    expr = allg[allg["hep_det"] >= MIN_DET]
    loc_rows.append({"locus": L["locus"], "peak_id": L["peak_id"], "max_pp": L["max_pp"],
                     "peak_hep_det": L["hep_det"], "peak_passes_5pct": L["hep_det"] >= MIN_DET,
                     "window_genes": len(allg), "window_genes_expressed_5pct": len(expr),
                     "coloc_genes": len(cset), "coloc_genes_expressed_5pct": int(expr["ensembl"].isin(cset).sum()),
                     "noncoloc_expressed_5pct": int((~expr["ensembl"].isin(cset)).sum())})
    if L["hep_det"] >= MIN_DET:
        for _, g in expr.iterrows():
            loc_pairs.append(("locus", L["locus"], int(L["peak_row"]), int(g["row"]), g["ensembl"],
                              int(g["ensembl"] in cset), abs(int(g["tss"]) - int(p["center"]))))
loc_df = pd.DataFrame(loc_rows)
loc_df.to_csv(f"{OUT}/locus_testability.tsv", sep="\t", index=False)
log("loci:", len(loc_df), "peak >= 5% hepatocyte nuclei:", int(loc_df["peak_passes_5pct"].sum()),
    "and >=1 COLOC gene expressed >=5% and >=1 other gene:",
    int((loc_df["peak_passes_5pct"] & (loc_df["coloc_genes_expressed_5pct"] > 0) & (loc_df["noncoloc_expressed_5pct"] > 0)).sum()))

# SORT1 window
sp = int(pk_idx[SORT1_PEAK]); p = pk.iloc[sp]
sw = gok[(gok["chrom"] == p["chromosome"]) & ((gok["tss"] - p["center"]).abs() <= WINDOW)]
sort_pairs = [("sort1", "L013", sp, int(g["row"]), g["ensembl"], int(g["name_rect"] == "SORT1"),
               abs(int(g["tss"]) - int(p["center"]))) for _, g in sw.iterrows()]
log("SORT1 peak hepatocyte detection", round(float(peak_det[sp]), 4), "window genes >=5%:", len(sw))

# Own-promoter and trans
pc = gok[gok["gtype"] == "protein_coding"].copy()
prom = []
for chrom, block in pk[pk["hep_det"] >= MIN_DET].groupby("chromosome"):
    s = block["start_1based"].to_numpy(); e = block["end_1based"].to_numpy(); r = block.index.to_numpy()
    o = np.argsort(s); s, e, r = s[o], e[o], r[o]
    gg = pc[pc["chrom"] == chrom]
    for _, g in gg.iterrows():
        i = np.searchsorted(s, g["tss"] + 1000, side="right") - 1
        while i >= 0 and e[i] >= g["tss"] - 1000:
            if s[i] <= g["tss"] + 1000:
                prom.append((int(r[i]), int(g["row"]), g["ensembl"], chrom)); break
            i -= 1
prom = pd.DataFrame(prom, columns=["peak_row", "gene_row", "ensembl", "chrom"]).drop_duplicates("gene_row")
log("protein-coding genes >=5% with a >=5% promoter peak (TSS +-1 kb):", len(prom))
sel = prom.iloc[rng.choice(len(prom), size=min(150, len(prom)), replace=False)]
own_pairs = [("own", r.chrom, int(r.peak_row), int(r.gene_row), r.ensembl, 1, 0) for r in sel.itertuples()]
trans_pairs = []
for r in sel.itertuples():
    cand = pc[pc["chrom"] != r.chrom]
    g = cand.iloc[int(rng.integers(len(cand)))]
    trans_pairs.append(("trans", r.chrom, int(r.peak_row), int(g["row"]), g["ensembl"], 0, -1))

pairs = pd.DataFrame(sort_pairs + own_pairs + trans_pairs + loc_pairs,
                     columns=["set", "group", "peak_row", "gene_row", "ensembl", "flag", "dist"])
log("pairs to fit:", pairs["set"].value_counts().to_dict())

# Dense matrices for hepatocyte nuclei
P = extract_rows(a_ptr, a_ind, None, np.unique(pairs["peak_row"]), hep_idx)
P = (P > 0).astype(np.float64)
prow = {r: i for i, r in enumerate(np.unique(pairs["peak_row"]))}
del a_ind
G = extract_rows(r_ptr, r_ind, r_dat, np.unique(pairs["gene_row"]), hep_idx).astype(np.float64)
grow = {r: i for i, r in enumerate(np.unique(pairs["gene_row"]))}
del r_ind, r_dat
log("dense P", P.shape, "G", G.shape)

d_h = donor[hep_idx]; w_h = well[hep_idx]
dn = sorted(set(d_h))
D = (d_h[:, None] == np.array(dn)[None, :]).astype(np.float64)
cov = np.column_stack([np.log(numi[hep_idx].astype(float) + 1), pct_mt[hep_idx]])
cov = (cov - cov.mean(0)) / cov.std(0)
half_a = np.isin(w_h, ["well1", "well2", "well3", "well4"])


def fit(y, x, mask):
    y = y[mask]; Dm = D[mask]
    keep_d = Dm.T @ y > 0
    rows = (Dm[:, keep_d].sum(1) > 0)
    y = y[rows]
    X = np.column_stack([x[mask][rows], cov[mask][rows], Dm[rows][:, keep_d]])
    # Fix after job 22245251 crashed: a gene with zero counts in every donor of a half leaves no rows.
    if X.shape[0] == 0 or X[:, 0].min() == X[:, 0].max():
        return np.nan, np.nan, np.nan, np.nan
    b = np.zeros(X.shape[1])
    # donor intercepts start at the log mean count of each donor
    dm = X[:, 3:]
    b[3:] = np.log((dm.T @ y) / dm.sum(0))
    for _ in range(30):
        eta = np.clip(X @ b, -30, 30); mu = np.exp(eta)
        z = eta + (y - mu) / mu
        XtW = X.T * mu
        H = XtW @ X
        bn = np.linalg.solve(H, XtW @ z)
        if np.max(np.abs(bn - b)) < 1e-8:
            b = bn; break
        b = bn
    eta = np.clip(X @ b, -30, 30); mu = np.exp(eta)
    Hi = np.linalg.inv((X.T * mu) @ X)
    se_m = np.sqrt(Hi[0, 0])
    meat = (X.T * (y - mu) ** 2) @ X
    se_r = np.sqrt((Hi @ meat @ Hi)[0, 0])
    return b[0], se_m, se_r, int(len(y))


def work(i):
    if time.time() - T0 > DEADLINE:
        return i, None
    r = pairs.iloc[i]
    y = G[grow[r.gene_row]]; x = P[prow[r.peak_row]]
    out = {}
    out["full"] = fit(y, x, np.ones(y.size, bool))
    if r.set in ("own", "sort1", "locus", "trans"):
        out["A"] = fit(y, x, half_a)
        out["B"] = fit(y, x, ~half_a)
    return i, out


own_i = list(pairs.index[pairs["set"] == "own"]); trans_i = list(pairs.index[pairs["set"] == "trans"])
inter = [v for ab in zip(own_i, trans_i) for v in ab]
order = list(pairs.index[pairs["set"] == "sort1"]) + inter + list(pairs.index[pairs["set"] == "locus"])
res = {}
with Pool(4) as pool:
    for i, out in pool.imap_unordered(work, order, chunksize=4):
        if out is not None:
            res[i] = out
log("fitted pairs:", len(res), "of", len(pairs))
rows = []
for i, out in res.items():
    r = pairs.iloc[i].to_dict()
    for k in ("full", "A", "B"):
        if k in out:
            b, sm, sr, n = out[k]
            r[f"beta_{k}"] = b; r[f"se_model_{k}"] = sm; r[f"se_robust_{k}"] = sr; r[f"n_{k}"] = n
            r[f"p_robust_{k}"] = 2 * stats.norm.sf(abs(b / sr)) if np.isfinite(sr) and sr > 0 else np.nan
            r[f"p_model_{k}"] = 2 * stats.norm.sf(abs(b / sm)) if np.isfinite(sm) and sm > 0 else np.nan
    r["peak_id"] = peaks[r["peak_row"]]
    r["peak_hep_det"] = peak_det[r["peak_row"]]; r["gene_hep_det"] = gene_det[r["gene_row"]]
    r["gene"] = genes.loc[r["gene_row"], "name_rect"]
    rows.append(r)
R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/pilot_pairs.tsv", sep="\t", index=False)


def bh(p):
    p = np.asarray(p, float); n = np.isfinite(p).sum(); q = np.full(p.size, np.nan)
    idx = np.flatnonzero(np.isfinite(p)); o = idx[np.argsort(p[idx])]
    qq = p[o] * n / np.arange(1, n + 1); qq = np.minimum.accumulate(qq[::-1])[::-1]
    q[o] = np.minimum(qq, 1); return q


summ = {}
for kind in ("robust", "model"):
    R[f"q_{kind}"] = bh(R[f"p_{kind}_full"])
    for s in ("own", "trans", "locus", "sort1"):
        x = R[R["set"] == s]
        if len(x) == 0:
            continue
        rep = (x[f"q_{kind}"] < 0.1) & (np.sign(x["beta_A"]) == np.sign(x["beta_full"])) & (np.sign(x["beta_B"]) == np.sign(x["beta_full"]))
        summ[f"{s}_{kind}"] = {"n": int(len(x)), "rate_p05": float((x[f"p_{kind}_full"] < 0.05).mean()),
                               "rate_q10": float((x[f"q_{kind}"] < 0.1).mean()),
                               "rate_q10_pos": float(((x[f"q_{kind}"] < 0.1) & (x["beta_full"] > 0)).mean()),
                               "rate_q10_halves_agree": float(rep.mean())}
R.to_csv(f"{OUT}/pilot_pairs.tsv", sep="\t", index=False)
x = R[R["set"] == "sort1"].sort_values("p_robust_full")
summ["sort1_rows"] = x[["gene", "flag", "dist", "beta_full", "se_robust_full", "p_robust_full", "q_robust", "beta_A", "beta_B", "gene_hep_det", "peak_hep_det"]].to_dict("records")
L = R[R["set"] == "locus"].copy()
if len(L):
    L["link"] = (L["q_robust"] < 0.1) & (np.sign(L["beta_A"]) == np.sign(L["beta_full"])) & (np.sign(L["beta_B"]) == np.sign(L["beta_full"]))
    summ["locus_link_rate_coloc"] = float(L.loc[L["flag"] == 1, "link"].mean()) if (L["flag"] == 1).any() else None
    summ["locus_link_rate_other"] = float(L.loc[L["flag"] == 0, "link"].mean()) if (L["flag"] == 0).any() else None
    summ["locus_pairs_coloc"] = int((L["flag"] == 1).sum()); summ["locus_pairs_other"] = int((L["flag"] == 0).sum())
    summ["loci_fitted"] = int(L["group"].nunique())
summ["n_hep"] = int(n_hep); summ["elapsed_s"] = round(time.time() - T0, 1)
json.dump(summ, open(f"{OUT}/pilot_summary.json", "w"), indent=1, default=float)
print(json.dumps(summ, indent=1, default=float))
log("done")

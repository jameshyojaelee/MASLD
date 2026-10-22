"""F4abcd blind check of R04-09: donor bootstrap of the 117 content-lineage calls.
Content call implemented from the published plc-20260907T105813Z method (union of genes,
NaN-aware gene z across six lineages with ddof=1, L1-weighted mean z, >=5 genes, argmax).
Joint scheme: donors drawn with replacement from the union; each lineage uses the drawn donors
it has (with multiplicity). Within scheme: donors drawn independently within each lineage.
numpy default_rng(20260923), 1000 draws each.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

R = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1]); OUT.mkdir(exist_ok=True)
SCDIR = R / "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
HS = R / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
LIN = ["Hepatocytes", "Fibroblasts", "Cholangiocytes", "Macrophages", "T_cells", "Endothelial_cells"]
L2L = {"hepatocytes": "Hepatocytes", "fibroblasts": "Fibroblasts", "cholangiocytes": "Cholangiocytes",
       "macrophages": "Macrophages", "tcells": "T_cells"}
B = 1000

s2d = {}
for ds in ("GSE244832", "GSE185477", "GSE202379", "GSE136103"):
    dp = pd.read_csv(R / f"data/{ds}/metadata/donor_pairing.csv", dtype=str)
    for _, r in dp.iterrows():
        for s in str(r["rna_srrs"]).split(";"):
            if s and s != "nan":
                s2d[s.strip()] = f"{ds}_{r['donor_id']}"

mats = {}
for ln in LIN:
    m = pd.read_csv(SCDIR / f"pseudobulk_{ln}.tsv.gz", sep="\t", index_col=0)
    m = m.groupby(level=0).sum()
    cols = pd.Series(m.columns).map(s2d).fillna(pd.Series(m.columns)).values
    m = m.T.groupby(cols).sum().T
    cpm = m.div(m.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
    lc = np.log2(cpm + 1)
    mats[ln] = lc.loc[:, lc.notna().any(axis=0)]   # published mean(skipna) ignores zero-library donors
    print(ln, m.shape, flush=True)
genes = sorted(set().union(*[set(v.index) for v in mats.values()]))
gi = {g: i for i, g in enumerate(genes)}
donors = sorted(set().union(*[set(v.columns) for v in mats.values()]))
di = {d: i for i, d in enumerate(donors)}
A = {}
for ln in LIN:
    v = mats[ln].reindex(genes)  # NaN rows for genes absent from this lineage table
    A[ln] = (v.to_numpy(np.float64), np.array([di[c] for c in v.columns]))
print("genes", len(genes), "donors", len(donors), flush=True)

reg = pd.read_csv(HS / "program_registry_v2.tsv", sep="\t")
mem = pd.read_csv(HS / "program_membership_v2.tsv", sep="\t")


def profile(weights_per_lineage):
    P = np.empty((len(genes), len(LIN)))
    for j, ln in enumerate(LIN):
        X, idx = A[ln]
        w = weights_per_lineage[j]
        P[:, j] = (X @ w) / w.sum()
    return P


def zscore(P):
    mu = np.nanmean(P, axis=1, keepdims=True)
    sd = np.nanstd(P, axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = np.nan
    return (P - mu) / sd


full_w = [np.ones(len(A[ln][1])) for ln in LIN]
Pfull = profile(full_w)
keep = ~np.all(np.isnan(Pfull), axis=1)
Zf = zscore(Pfull)
# program weight vectors over genes with finite z in the full data (published: w in Z.index after dropna(how=all))
uids, Wrows, declared = [], [], []
for uid, g in mem.groupby("program_uid"):
    w = g.dropna(subset=["mapped_symbol"]).copy()
    w["wt"] = pd.to_numeric(w["original_l1_weight"], errors="coerce")
    w = w.dropna(subset=["wt"])
    w = w[w.mapped_symbol.map(lambda s: s in gi and keep[gi[s]])]
    if len(w) < 5:
        continue
    wt = w.groupby("mapped_symbol")["wt"].sum(); wt = wt / wt.sum()
    vec = np.zeros(len(genes)); vec[[gi[s] for s in wt.index]] = wt.values
    uids.append(uid); Wrows.append(vec)
    declared.append(L2L[reg.loc[reg.program_uid == uid, "cell_type"].iloc[0]])
W = np.vstack(Wrows); declared = np.array(declared)
print("programs scored", len(uids), flush=True)


def calls(Z):
    Zc = np.where(np.isnan(Z), 0.0, Z)   # pandas sum(skipna) treats NaN as 0
    S = W @ Zc
    return np.array(LIN)[S.argmax(axis=1)]


cf = calls(Zf)
pub = pd.read_csv(R / "Analysis/SingleCell/results_gpu_v2/program_label_vs_content/plc-20260907T105813Z/program_label_vs_content.tsv", sep="\t")
pubc = pub.set_index("program_uid").loc[uids, "content_lineage"].values
print("full-data calls equal published:", int((cf == pubc).sum()), "of", len(uids),
      "; match declared:", int((cf == declared).sum()), flush=True)

rng = np.random.default_rng(20260923)
res = {}
for scheme in ["joint", "within"]:
    C = np.empty((B, len(uids)), dtype=object)
    for b in range(B):
        if scheme == "joint":
            draw = rng.integers(0, len(donors), len(donors))
            cnt = np.bincount(draw, minlength=len(donors)).astype(float)
            ws = [cnt[A[ln][1]] for ln in LIN]
            ws = [w if w.sum() > 0 else np.ones_like(w) for w in ws]
        else:
            ws = []
            for ln in LIN:
                n = len(A[ln][1]); d = rng.integers(0, n, n)
                ws.append(np.bincount(d, minlength=n).astype(float))
        C[b] = calls(zscore(profile(ws)))
    same = (C == cf[None, :]).mean(axis=0)
    stable = same >= 0.95
    match_full = cf == declared
    nmatch = (C == declared[None, :]).sum(axis=1)
    lo, hi = np.percentile(nmatch, [2.5, 97.5])
    res[scheme] = dict(stable_match=int((stable & match_full).sum()), stable_mismatch=int((stable & ~match_full).sum()),
                       unstable=int((~stable).sum()), nmatch_lo=lo, nmatch_hi=hi, nmatch_median=np.median(nmatch))
    pd.DataFrame({"program_uid": uids, "declared": declared, "full_call": cf, "frac_same": same}).to_csv(
        OUT / f"content_bootstrap_{scheme}.tsv", sep="\t", index=False)
    print(scheme, res[scheme], flush=True)
    for u in ["hotspot_hepatocytes_f05c535ae5bbc0b9", "hotspot_hepatocytes_48f39dd4d817a10e"]:
        print("  ", u, cf[uids.index(u)], same[uids.index(u)])
pd.DataFrame(res).T.to_csv(OUT / "content_bootstrap_summary.tsv", sep="\t")

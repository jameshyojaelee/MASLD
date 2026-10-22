"""H90 job 6302 (a): donor-mean standardized genes per lineage, frozen-weight program
scores (117), NC3 matched random gene sets, E6 pseudobulk counts, E10 doublet-
sensitivity scores. Runs only after CHECKED_BEFORE_OUTCOME.txt says PROCEED: yes."""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).parent))
import h90_common as C

chk = C.RUN / "CHECKED_BEFORE_OUTCOME.txt"
if not chk.exists() or "PROCEED: yes" not in chk.read_text():
    sys.exit("STOP: CHECKED_BEFORE_OUTCOME.txt missing or not PROCEED: yes")
C.check_inputs()
out = C.OUT
A = ad.read_h5ad(C.DERIVED / "gse289173_qc_nuclei.h5ad")
X = A.X.tocsr() if sparse.issparse(A.X) else sparse.csr_matrix(A.X)
obs = A.obs.reset_index(drop=True)
genes = pd.Index(A.var_names.astype(str))
gi = {g: i for i, g in enumerate(genes)}
del A
dq = pd.read_csv(out / "donor_qc.tsv", sep="\t").sort_values("gsm").reset_index(drop=True)
pc1_excl = set(dq.gsm[~dq.pc1_concordant.astype(bool)])
mem = C.program_table()
cov = pd.read_csv(out / "D3_coverage.tsv", sep="\t")

score_rows, gene_rows, cov_rows = [], [], []
hep_state = {}
for ct, lab in C.LINEAGE_LABEL.items():
    lm = (obs.label.astype(str) == lab).to_numpy()
    det_all = np.asarray((X[lm] > 0).mean(axis=0)).ravel() if lm.any() else np.zeros(len(genes))
    lob = obs[lm].reset_index(drop=True)
    nper = lob.groupby("gsm").size()
    elig = sorted(g for g, n in nper.items() if n >= C.MIN_NUCLEI and g not in pc1_excl)
    for g in dq.gsm:
        cov_rows.append(dict(cell_type=ct, gsm=g, n_nuclei=int(nper.get(g, 0)), eligible=g in elig))
    if not elig:
        continue
    em = lob.gsm.isin(elig).to_numpy()
    Xl = X[np.flatnonzero(lm)[em]]
    eob = lob[em].reset_index(drop=True)
    codes = pd.Categorical(eob.gsm, categories=elig).codes
    progs = sorted(mem.loc[mem.cell_type == ct, "module"].unique())
    need = set()
    wtab = {}
    for mod in progs:
        w = C.program_weights(mem, ct, mod)
        w = w[w.mapped_symbol.isin(gi)]
        w = w[w.mapped_symbol.map(lambda g: det_all[gi[g]] >= 0.01)]
        wtab[mod] = w
        need |= set(w.mapped_symbol)
    if ct == "hepatocytes":
        need |= {g for g in genes if det_all[gi[g]] >= 0.01}
    need = sorted(need)
    cols = np.array([gi[g] for g in need])
    norm = C.log_normalize(Xl[:, cols].tocsr(), eob.lib.to_numpy(float))
    dmz, mu, sd = C.donor_mean_z(norm, codes, len(elig))
    ni = {g: i for i, g in enumerate(need)}
    for mod in progs:
        w = wtab[mod]
        retained = float(w.original_l1_weight.sum())
        testable = len(w) >= 8 and retained >= 0.20
        if testable:
            wv = (w.original_l1_weight / retained).to_numpy()
            sc = dmz[:, [ni[g] for g in w.mapped_symbol]] @ wv
            for g, s in zip(elig, sc):
                score_rows.append(dict(cell_type=ct, module=mod, gsm=g, score=s))
            for g, ww in zip(w.mapped_symbol, wv):
                gene_rows.append(dict(cell_type=ct, module=mod, gene=g, weight=ww))
    if ct == "hepatocytes":
        det_pool = np.asarray((norm > 0).mean(axis=0)).ravel()
        mean_pool = np.asarray(norm.mean(axis=0)).ravel()
        hep_state = dict(elig=elig, dmz=dmz, need=need, ni=ni, det_all=det_all, det_pool=det_pool,
                         mean_pool=mean_pool, Xl=Xl, eob=eob, codes=codes, wtab=wtab)

scores = pd.DataFrame(score_rows)
scores.to_csv(out / "donor_scores.tsv", sep="\t", index=False)
pd.DataFrame(gene_rows).to_csv(out / "scored_genes_weights.tsv", sep="\t", index=False)
lcov = pd.DataFrame(cov_rows)
lcov.to_csv(out / "lineage_donor_nuclei.tsv", sep="\t", index=False)

# ---- hepatocyte covariates ----------------------------------------------------------
H = hep_state
eob = H["eob"]
hc = eob.groupby("gsm").agg(n_hep=("umi", "size"), median_umi_hep=("umi", "median")).reindex(H["elig"]).reset_index()
hc = hc.merge(dq[["gsm", "title", "group", "stage", "sex", "series"]], on="gsm", how="left")
hc.to_csv(out / "donor_covariates_hep.tsv", sep="\t", index=False)

# ---- NC3 matched random gene sets (competitive null) -----------------------------
w20 = H["wtab"][C.PRIMARY[1]]
wv = (w20.original_l1_weight / w20.original_l1_weight.sum()).to_numpy()
m20 = mem[(mem.cell_type == C.PRIMARY[0]) & (mem.module == C.PRIMARY[1])]
excl = set(m20.source_gene.astype(str)) | set(m20.mapped_symbol.dropna().astype(str))
need = H["need"]
pool_df = pd.DataFrame({"gene": need, "mean": H["mean_pool"], "det": H["det_pool"],
                        "scored": [H["det_all"][gi[g]] >= 0.01 for g in need]})
pool_df = pool_df[pool_df.scored].reset_index(drop=True)
pool_df["mbin"] = pd.qcut(pool_df["mean"].rank(method="first"), 10, labels=False)
pool_df["dbin"] = pd.qcut(pool_df["det"].rank(method="first"), 10, labels=False)
pool_df["excluded"] = pool_df.gene.isin(excl)
cand = pool_df[~pool_df.excluded]
cells = {k: v.gene.to_numpy() for k, v in cand.groupby(["mbin", "dbin"])}
mrow = {k: v.gene.to_numpy() for k, v in cand.groupby("mbin")}
bins = pool_df.set_index("gene").loc[w20.mapped_symbol, ["mbin", "dbin"]]
rng = np.random.default_rng(C.SEED)
nc3 = np.zeros((1000, len(H["elig"])))
fallback = 0
sets = []
for k in range(1000):
    chosen = []
    for (mb, db) in bins.itertuples(index=False):
        c = [g for g in cells.get((mb, db), []) if g not in chosen]
        if not c:
            c = [g for g in mrow[mb] if g not in chosen]
            fallback += 1
        chosen.append(c[rng.integers(len(c))])
    sets.append(chosen)
    nc3[k] = H["dmz"][:, [H["ni"][g] for g in chosen]] @ wv
pd.DataFrame(nc3, columns=H["elig"]).to_csv(out / "NC3_set_scores.tsv.gz", sep="\t", index=False)
pd.DataFrame(sets).to_csv(out / "NC3_sets_genes.tsv.gz", sep="\t", index=False, header=False)
pool_df.to_csv(out / "NC3_pool_bins.tsv.gz", sep="\t", index=False)
pd.Series({"n_sets": 1000, "genes_per_set": len(wv), "pool_genes": len(cand), "fallback_draws": fallback}).to_csv(
    out / "NC3_meta.tsv", sep="\t", header=["value"])

# ---- E6 pseudobulk counts (primary scored genes) and library sizes ----------------
Xl, codes = H["Xl"], H["codes"]
ind = sparse.csr_matrix((np.ones(len(eob)), (codes, np.arange(len(eob)))), shape=(len(H["elig"]), len(eob)))
pb = np.asarray((ind @ Xl[:, [gi[g] for g in w20.mapped_symbol]]).todense())
e6 = pd.DataFrame(pb, columns=list(w20.mapped_symbol))
e6.insert(0, "lib_size", np.asarray(ind @ eob.lib.to_numpy(float)).ravel())
e6.insert(0, "gsm", H["elig"])
e6.to_csv(out / "E6_pseudobulk_counts.tsv", sep="\t", index=False)

# ---- E10 doublet sensitivity: drop hepatocyte nuclei with any KRT7 or KRT19 count --
krt = np.asarray((Xl[:, [gi["KRT7"], gi["KRT19"]]] > 0).sum(axis=1)).ravel() > 0
keep = ~krt
cols = np.array([gi[g] for g in w20.mapped_symbol])
norm10 = C.log_normalize(Xl[keep][:, cols].tocsr(), eob.lib.to_numpy(float)[keep])
codes10 = codes[keep]
present = sorted(set(codes10))
if len(present) != len(H["elig"]):
    print("E10: some donors lose all hepatocyte nuclei; they are dropped from E10")
remap = {c: i for i, c in enumerate(present)}
dmz10, _, _ = C.donor_mean_z(norm10, np.array([remap[c] for c in codes10]), len(present))
e10 = pd.DataFrame({"gsm": [H["elig"][c] for c in present], "score": dmz10 @ wv,
                    "n_hep_kept": np.bincount(codes10, minlength=len(H["elig"]))[present]})
e10.to_csv(out / "E10_scores.tsv", sep="\t", index=False)
print("04 done: programs scored", scores.groupby(["cell_type", "module"]).ngroups, "; NC3 fallback draws", fallback)

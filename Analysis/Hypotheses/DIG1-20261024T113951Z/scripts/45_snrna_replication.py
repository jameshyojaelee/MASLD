"""DIG1 job 45 (exploratory replication of jobs 34, 36b, 37b, 39). GSE256398 snRNA, 17 donors (6 healthy,
11 MASLD F0 to cirrhosis). Written before any GSE256398 cell-level coupling was computed.

Question: do hepatocyte nuclei with more KRT8/SQSTM1 RNA (MDB score, the GSE298719 genes) also carry more CCT, HSF1-set
and HSPB1 RNA and less PPARGC1A, as GSE298719 Flex hepatocytes do (§25-§28)? Nuclei hold mostly newly
made RNA, so a replication here would point to co-transcription rather than cytoplasmic mRNA
stability [Inference].
Labels: V3 route through V0's build_labelled, computed once by job 45L (out_snrna_labels/, which records
how far the labels are from V3's; exact reproduction needs an Intel v4 node).
Statistic (as job 39; re-derived on these nuclei, not ported): per donor, rho(MDB, test) with ranks
residualised on log UMI and log genes (depth) or depth + ranks of central and portal scores
(depth+zonation), minus the mean over 200 random sets matched gene by gene on hepatocyte-nucleus mean x
detection deciles (set, MDB and zonation genes excluded; seed 20261024). Exact sign-flip over donors.
Primary family (Holm over 4), depth+zonation, all 17 donors:
  R1 CCT, R2 HSF1 set (HSPA1A HSPA1B DNAJB1 HSPH1 BAG3 HSPA4L HSP90AA1), R3 HSPB1, R4 PPARGC1A.
  Replicated = same sign as GSE298719 (R1-R3 positive, R4 negative) and Holm < 0.05.
Descriptive: all job 39 sets, both adjustments (BH); excess by group; excess vs stage (Spearman over
donors); hepatocyte donor means, MASLD vs healthy + sex, exact relabeling (12,376 splits), and vs stage.
Usage: python -I 45_snrna_replication.py <new_out_dir>
"""
import sys
from itertools import combinations, product
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V3PB = ROOT / "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
rng = np.random.default_rng(20261024)
N_NULL = 200


LBL = RUN / "out_snrna_labels"
hep = ad.read_h5ad(LBL / "gse256398_hepatocyte_counts.h5ad")
donors = sorted(hep.obs["donor"].astype(str).unique())
assert len(donors) == 17

meta = pd.read_csv(V3PB / "donor_metadata.tsv", sep="\t").set_index("donor").loc[donors].reset_index()
X = sparse.csr_matrix(hep.X).astype(np.float32)
genes = pd.Index(hep.var_names); g = list(genes)
dH = hep.obs["donor"].astype(str).to_numpy()
del hep
umi = np.asarray(X.sum(axis=1)).ravel(); ngene = np.asarray((X > 0).sum(axis=1)).ravel()
L = X.multiply(1e4 / umi[:, None]).tocsr(); L.data = np.log1p(L.data)
ZON = {"central": ["CYP2E1", "GLUL", "CYP1A2", "OAT", "CYP3A4"], "portal": ["CPS1", "HAL", "ASS1", "SDS", "PCK1", "ALDOB"]}
MDB = ["KRT8", "SQSTM1"]  # same genes as GSE298719 (no KRT18 on Flex); job 46 adds KRT18 as a sensitivity
TESTS = {
    "CCT": ["TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8"],
    "HSF1": ["HSPA1A", "HSPA1B", "DNAJB1", "HSPH1", "BAG3", "HSPA4L", "HSP90AA1"],
    "chaperone_network": ["HSPA8", "HSP90AA1", "HSP90AB1", "STIP1", "AHSA1", "CHORDC1", "DNAJA1", "HSPA4", "HSPA4L", "HSPH1"],
    "tubulin": ["TUBA1A", "TUBA1B", "TUBA1C", "TUBA4A", "TUBB", "TUBB2A", "TUBB4B", "TUBB6"],
    "actin": ["ACTB", "ACTG1"],
    "CI": [x for x in g if x.startswith("NDUF") and x[4:5] in "ABSCV" and x[5:6].isdigit() and x != "NDUFA4"],
    "mito_chaperones": ["HSPD1", "HSPE1"],
    "HSPB1": ["HSPB1"], "HSPA8": ["HSPA8"], "HSP90AA1": ["HSP90AA1"], "HSPA1A": ["HSPA1A"], "PPARGC1A": ["PPARGC1A"], "CPT1A": ["CPT1A"],
}
PRIMARY = {"CCT": ("R1", 1), "HSF1": ("R2", 1), "HSPB1": ("R3", 1), "PPARGC1A": ("R4", -1)}
mean_all = np.asarray(L.mean(axis=0)).ravel(); det_all = np.asarray((X > 0).mean(axis=0)).ravel()
keep = det_all >= 0.01
setgenes = sorted(set(x for v in {**TESTS, **ZON, "MDB": MDB}.values() for x in v if x in genes))
pd.DataFrame(dict(gene=setgenes, hep_mean=[mean_all[genes.get_loc(x)] for x in setgenes], hep_detection=[det_all[genes.get_loc(x)] for x in setgenes],
                  binned=[bool(keep[genes.get_loc(x)]) for x in setgenes])).to_csv(OUT / "hep_gene_detection.tsv", sep="\t", index=False)
TESTS = {k: [x for x in v if x in genes] for k, v in TESTS.items()}
TESTS = {k: v for k, v in TESTS.items() if v}
ZON = {k: [x for x in v if x in genes] for k, v in ZON.items()}
MDB = [x for x in MDB if x in genes]
pd.DataFrame([dict(set=k, n=len(v), genes=",".join(v)) for k, v in {**TESTS, **ZON, "MDB": MDB}.items()]).to_csv(OUT / "sets_used.tsv", sep="\t", index=False)
excl = set(x for v in TESTS.values() for x in v) | set(x for v in ZON.values() for x in v) | set(MDB)
binm = pd.qcut(pd.Series(mean_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
bind = pd.qcut(pd.Series(det_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
bins = pd.Series(list(zip(binm, bind)), index=np.array(g)[keep])
pool = {b: [x for x in idx if x not in excl] for b, idx in bins.groupby(bins).groups.items()}
NULL = {t: [[rng.choice(pool[bins[x]]) for x in v if x in bins.index] for _ in range(N_NULL)] for t, v in TESTS.items()}
need = sorted(excl | set(x for r in NULL.values() for s in r for x in s))
Ld = L[:, [genes.get_loc(x) for x in need]].toarray(); gi = {x: i for i, x in enumerate(need)}
del L


def resid_rank(M, Z):
    R = np.apply_along_axis(stats.rankdata, 0, M)
    return R - Z @ np.linalg.lstsq(Z, R, rcond=None)[0]


def cor_cols(a, M):
    a = (a - a.mean()) / a.std(); sd = M.std(0); sd[sd == 0] = np.nan
    return ((M - M.mean(0)) / sd * a[:, None]).mean(0)


rows, dmeans = [], []
for d in donors:
    m = dH == d
    D = Ld[m].astype(np.float64); sd = D.std(0); sd[sd == 0] = np.nan
    Zg = (D - D.mean(0)) / sd
    sc = lambda gs: np.nanmean(Zg[:, [gi[x] for x in gs]], 1)
    mdb, cen, por = sc(MDB), sc(ZON["central"]), sc(ZON["portal"])
    dmeans.append(dict(donor=d, n_hep=int(m.sum()), **{k: float(D[:, [gi[x] for x in v]].mean()) for k, v in {**TESTS, "MDB": MDB}.items()}))
    base = np.column_stack([np.ones(m.sum()), stats.rankdata(np.log(umi[m])), stats.rankdata(np.log(ngene[m]))])
    zadj = np.column_stack([base, stats.rankdata(cen), stats.rankdata(por)])
    for adj, Zr in [("depth", base), ("depth+zonation", zadj)]:
        ri = resid_rank(mdb[:, None], Zr)[:, 0]
        for t, v in TESTS.items():
            s = sc(v)
            if not np.isfinite(s).all():
                continue
            ro = resid_rank(s[:, None], Zr)[:, 0]
            Rn = resid_rank(np.column_stack([np.nan_to_num(sc(r)) for r in NULL[t]]), Zr)  # null set undetected in a donor: constant column, dropped by cor_cols
            rows.append(dict(donor=d, adjust=adj, test=t, n_cells=int(m.sum()), rho=np.corrcoef(ri, ro)[0, 1], null_mean=np.nanmean(cor_cols(ri, Rn))))
    print(d, int(m.sum()), flush=True)
wd = pd.DataFrame(rows).merge(meta[["donor", "stage", "stage_ordinal", "disease"]], on="donor")
wd["excess"] = wd.rho - wd.null_mean
wd.to_csv(OUT / "within_donor_rho.tsv", sep="\t", index=False)
SIGNS = {}


def signflip(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]; n = len(v)
    if n not in SIGNS:
        SIGNS[n] = np.array(list(product([1, -1], repeat=n)), dtype=np.int8)
    return np.mean(np.abs((SIGNS[n] * v).mean(1)) >= abs(v.mean()) - 1e-12)


summ = []
for (adj, t), o in wd.groupby(["adjust", "test"]):
    o = o.set_index("donor").reindex(donors)
    e = o.excess
    summ.append(dict(adjust=adj, test=t, n_donors=int(e.notna().sum()), mean_rho=o.rho.mean(), mean_null=o.null_mean.mean(), mean_excess=e.mean(),
                     n_excess_positive=int((e > 0).sum()), p_excess=signflip(e),
                     excess_healthy=e[o.disease == 0].mean(), excess_masld=e[o.disease == 1].mean(),
                     rho_excess_stage=stats.spearmanr(o.stage_ordinal[e.notna()], e[e.notna()])[0]))
summ = pd.DataFrame(summ)
summ["bh"] = summ.groupby("adjust").p_excess.transform(stats.false_discovery_control)
summ.to_csv(OUT / "within_donor_summary.tsv", sep="\t", index=False)
prim = summ[(summ.adjust == "depth+zonation") & summ.test.isin(list(PRIMARY))].copy()
prim["id"] = prim.test.map(lambda t: PRIMARY[t][0]); prim["expected_sign"] = prim.test.map(lambda t: PRIMARY[t][1])
prim = prim.sort_values("id").reset_index(drop=True)
o = np.argsort(prim.p_excess.to_numpy()); ph = np.maximum.accumulate(prim.p_excess.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim["replicated"] = (np.sign(prim.mean_excess) == prim.expected_sign) & (prim.p_holm < 0.05)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)

dm = pd.DataFrame(dmeans).merge(meta[["donor", "stage_ordinal", "disease", "sex"]], on="donor")
dm.to_csv(OUT / "hep_donor_means.tsv", sep="\t", index=False)
n = len(dm); dis = dm.disease.to_numpy(float); sexv = (dm.sex == "male").to_numpy(float)
RELAB = []
for ix in combinations(range(n), int(dis.sum())):
    v = np.zeros(n); v[list(ix)] = 1; RELAB.append(v)
P = np.stack([np.linalg.pinv(np.column_stack([np.ones(n), v, sexv]))[1] for v in RELAB])
p_obs = np.linalg.pinv(np.column_stack([np.ones(n), dis, sexv]))[1]
sets = list(TESTS) + ["MDB"]
Y = dm[sets].to_numpy(float)
b_obs = p_obs @ Y; b_null = P @ Y
de = pd.DataFrame(dict(set=sets, beta_masld_minus_healthy=b_obs, p_exact=(np.abs(b_null) >= np.abs(b_obs) - 1e-12).mean(0),
                       rho_stage=[stats.spearmanr(dm.stage_ordinal, dm[s])[0] for s in sets], p_stage=[stats.spearmanr(dm.stage_ordinal, dm[s])[1] for s in sets]))
de["bh"] = stats.false_discovery_control(de.p_exact)
de.to_csv(OUT / "hep_disease_effect.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)
print(prim.round(4).to_string())
print(summ.round(4).to_string())
print(de.round(4).to_string())
print(pd.read_csv(OUT / "hep_gene_detection.tsv", sep="\t").round(4).to_string())

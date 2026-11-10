"""DIG1 job 24 (exploratory; follows jobs 22 and 22b). GSE298719 Flex hepatocytes, 13 donors.

In liver protein (PXD051911) Complex I, IV, V and mitoribosome protein fall with ballooning while
Complex II is spared, with no matching fall in bulk RNA. Question: within a donor, do the most
injured hepatocytes carry less OXPHOS / mitoribosome RNA? If not, the protein loss in injured
cells is post-transcriptional [Inference to be tested here, not assumed].

Injury scores (fixed before reading expression; no OXPHOS genes):
  primary   MDB = KRT8, KRT18, SQSTM1 (Mallory-Denk body proteins; protein axis of job 22)
  secondary stress = CDKN1A, GDF15, SOX4, LCN2 (job 7 stress set without KRT8)
Per-cell gene values: counts per 10k, log1p, z within donor; set score = mean z.
Within donor: Spearman of injury score with each set score after residualising ranks on log UMI
and log genes detected. Donor statistic = that rho.
Primary family (Holm over 2):
  P1  within-donor rho(MDB, CI - CII), exact sign-flip over 13 donors (8,192 patterns).
  P2  disease effect on donor-mean CI - CII, OLS ~ disease + sex, exact relabeling (1,716 splits).
Descriptive: every set vs both injury scores, with a per-donor null from 200 random gene sets
matched on hepatocyte mean and detection deciles (seed 20261024, mito genes excluded from the
pool); report the mean over donors of the observed rho minus the matched-null mean, BH over sets
x injury scores on sign-flip p. Disease effect on donor means of every set (exact relabeling).
Usage: python -I 24_hep_injury_oxphos_rna.py <new_out_dir>
"""
import sys
from itertools import combinations, product
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
V5 = ROOT / "Analysis/Hypotheses/V5-20261021T122723Z"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
rng = np.random.default_rng(20261024)
N_NULL = 200

A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str})
dq = dq[dq.included].sort_values("gsm").reset_index(drop=True)
hep = A[(A.obs.label == "Hepatocytes").to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
del A
genes = pd.Index(hep.var_names)
X = sparse.csr_matrix(hep.X).astype(np.float32)
gsm = hep.obs.gsm.astype(str).to_numpy()
umi = np.asarray(X.sum(axis=1)).ravel()
ngene = np.asarray((X > 0).sum(axis=1)).ravel()
L = X.multiply(1e4 / umi[:, None]).tocsr(); L.data = np.log1p(L.data)

g = list(genes)
CI = [x for x in g if x.startswith("NDUF") and x[4:5] in "ABSCV" and x[5:6].isdigit() and x != "NDUFA4"]
SETS = {
    "CI": CI, "CII": [x for x in ["SDHA", "SDHB", "SDHC", "SDHD"] if x in genes],
    "CIII": [x for x in g if x.startswith("UQCR") or x == "CYC1"],
    "CIV": [x for x in g if (x.startswith("COX") and x[3:4] in "45678")] + [x for x in ["NDUFA4"] if x in genes],
    "CV": [x for x in g if x.startswith("ATP5")],
    "mtDNA": [x for x in g if x.startswith("MT-")],
    "mt_ribosome": [x for x in g if (x.startswith("MRPL") or x.startswith("MRPS")) and x[4:5].isdigit()],
    "CI_assembly": [x for x in [f"NDUFAF{i}" for i in range(1, 9)] + ["TIMMDC1", "ACAD9", "ECSIT", "FOXRED1", "TMEM126B", "NUBPL", "DMAC1", "DMAC2"] if x in genes],
    "mt_nucleoid": [x for x in ["TFAM", "POLG", "SSBP1", "LRPPRC", "TWNK"] if x in genes],
    "mt_translation": [x for x in ["TUFM", "TSFM", "GFM1", "GFM2", "MTIF2", "MTIF3", "TACO1", "LRPPRC", "SLIRP"] if x in genes],
    "mito_mass": [x for x in ["VDAC1", "VDAC2", "VDAC3", "SLC25A4", "SLC25A5", "SLC25A6", "TOMM20", "TOMM40"] if x in genes],
    "hep_mito": [x for x in ["CPS1", "OTC", "HMGCS2", "GLS2"] if x in genes],
    "FAO": [x for x in ["CPT2", "ACADM", "ACADL", "ACADVL", "ACADS", "HADHA", "HADHB", "ECHS1", "HADH", "ETFA", "ETFB", "ETFDH"] if x in genes],
    "ATF4_targets": [x for x in ["PHGDH", "PSAT1", "PSPH", "SHMT2", "MTHFD2", "MTHFD1L", "ALDH18A1", "PYCR1", "ASNS", "DDIT3", "TRIB3", "CHAC1"] if x in genes],
    "mt_UPR": [x for x in ["HSPD1", "HSPE1", "HSPA9", "LONP1", "CLPP", "YME1L1"] if x in genes],
    "glycolysis": [x for x in ["PKM", "LDHA", "PFKP", "ALDOA", "TPI1", "GAPDH", "ENO1", "PGK1"] if x in genes],
    "LD_coat": [x for x in ["PLIN2", "PLIN4", "ABHD5"] if x in genes],
}
INJ = {"MDB": [x for x in ["KRT8", "KRT18", "SQSTM1"] if x in genes],
       "stress": [x for x in ["CDKN1A", "GDF15", "SOX4", "LCN2"] if x in genes]}
SETS = {k: v for k, v in SETS.items() if len(v) > 0}
pd.DataFrame([dict(set=k, n=len(v), genes=",".join(v)) for k, v in {**SETS, **INJ}.items()]).to_csv(OUT / "sets_used.tsv", sep="\t", index=False)

# matched-null pool: hepatocyte mean and detection deciles, all mito-related and injury genes excluded
mean_all = np.asarray(L.mean(axis=0)).ravel()
det_all = np.asarray((X > 0).mean(axis=0)).ravel()
keep = det_all >= 0.01
mito_like = set(x for v in SETS.values() for x in v) | set(x for v in INJ.values() for x in v)
mito_like |= set(x for x in g if x.startswith(("MRP", "NDUF", "COX", "ATP5", "UQCR", "SDH", "MT-", "TIMM", "TOMM", "SLC25")))
binm = pd.qcut(pd.Series(mean_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
bind = pd.qcut(pd.Series(det_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
kg = np.array(g)[keep]
bins = pd.Series(list(zip(binm, bind)), index=kg)
pool = {b: [x for x in idx if x not in mito_like] for b, idx in bins.groupby(bins).groups.items()}
NULL = {}
for s, gs in SETS.items():
    gs = [x for x in gs if x in bins.index]
    NULL[s] = [[rng.choice(pool[bins[x]]) for x in gs] for _ in range(N_NULL)] if gs else []
need = sorted(set(x for v in SETS.values() for x in v) | set(x for v in INJ.values() for x in v)
              | set(x for v in NULL.values() for r in v for x in r))
col = [genes.get_loc(x) for x in need]
Ld = L[:, col].toarray()
gi = {x: i for i, x in enumerate(need)}


def resid_rank(M, Z):
    R = np.apply_along_axis(stats.rankdata, 0, M)
    return R - Z @ np.linalg.lstsq(Z, R, rcond=None)[0]


donor_rows, donor_means = [], []
for d in dq.gsm:
    m = gsm == d
    D = Ld[m]
    sd = D.std(0); sd[sd == 0] = np.nan
    Zg = (D - D.mean(0)) / sd
    sc = lambda gs: np.nanmean(Zg[:, [gi[x] for x in gs]], 1)
    obs = {s: sc(v) for s, v in SETS.items()}
    obs["CI_minus_CII"] = obs["CI"] - obs["CII"]
    inj = {k: sc(v) for k, v in INJ.items()}
    Zc = np.column_stack([np.ones(m.sum()), stats.rankdata(np.log(umi[m])), stats.rankdata(np.log(ngene[m]))])
    names = list(obs)
    Ro = resid_rank(np.column_stack([obs[s] for s in names]), Zc)
    Rnull = {s: resid_rank(np.column_stack([sc(r) for r in NULL[s]]), Zc) for s in names if s in NULL and NULL[s]}
    for k, iv in inj.items():
        ri = resid_rank(iv[:, None], Zc)[:, 0]
        rho = {s: np.corrcoef(ri, Ro[:, j])[0, 1] for j, s in enumerate(names)}
        for s in names:
            nm = np.nan
            if s in Rnull:
                Rn = Rnull[s]
                nm = np.mean([np.corrcoef(ri, Rn[:, j])[0, 1] for j in range(Rn.shape[1])])
            donor_rows.append(dict(gsm=d, injury=k, set=s, n_cells=int(m.sum()), rho=rho[s], null_mean=nm))
    # donor-level mean expression (log1p CP10k averaged over cells), for the disease effect
    mu = D.mean(0)
    dm = {s: np.mean([mu[gi[x]] for x in v]) for s, v in SETS.items()}
    dm["CI_minus_CII"] = dm["CI"] - dm["CII"]
    donor_means.append(dict(gsm=d, **dm))
    print(d, int(m.sum()), flush=True)

dr = pd.DataFrame(donor_rows).merge(dq[["gsm", "id", "disease", "sex_f", "fibrosis"]], on="gsm")
dr["excess"] = dr.rho - dr.null_mean
dr.to_csv(OUT / "within_donor_rho.tsv", sep="\t", index=False)
dmn = pd.DataFrame(donor_means).merge(dq[["gsm", "id", "disease", "sex_f"]], on="gsm")
dmn.to_csv(OUT / "donor_means.tsv", sep="\t", index=False)

n = len(dq)
SIGNS = np.array(list(product([1, -1], repeat=n)))


def signflip(v):
    v = np.asarray(v, float)
    return np.mean(np.abs((SIGNS * v).mean(1)) >= abs(v.mean()) - 1e-12)


RELAB = []
for ix in combinations(range(n), int(dq.disease.sum())):
    v = np.zeros(n, int); v[list(ix)] = 1; RELAB.append(v)
RELAB = np.array(RELAB)


def dis(y):
    y = np.asarray(y, float)
    beta = lambda v: np.linalg.lstsq(np.column_stack([np.ones(n), v, dq.sex_f]).astype(float), y, rcond=None)[0][1]
    b = beta(dq.disease.to_numpy())
    return b, np.mean(np.abs([beta(v) for v in RELAB]) >= abs(b) - 1e-12)


p1v = dr[(dr.injury == "MDB") & (dr.set == "CI_minus_CII")].set_index("gsm").rho.reindex(dq.gsm)
b2, p2 = dis(dmn.set_index("gsm").CI_minus_CII.reindex(dq.gsm))
prim = pd.DataFrame([dict(id="P1", test="within-donor rho(MDB, CI-CII)", est=p1v.mean(), n_positive=int((p1v > 0).sum()), p=signflip(p1v)),
                     dict(id="P2", test="disease effect on donor-mean CI-CII (log1p CP10k)", est=b2, n_positive=np.nan, p=p2)])
o = np.argsort(prim.p.to_numpy()); ph = np.maximum.accumulate(prim.p.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)

summ = []
for (k, s), o in dr.groupby(["injury", "set"]):
    o = o.set_index("gsm").reindex(dq.gsm)
    summ.append(dict(injury=k, set=s, mean_rho=o.rho.mean(), n_rho_negative=int((o.rho < 0).sum()), p_rho=signflip(o.rho),
                     mean_null=o.null_mean.mean(), mean_excess=o.excess.mean(),
                     n_excess_negative=int((o.excess < 0).sum()) if o.excess.notna().all() else np.nan,
                     p_excess=signflip(o.excess) if o.excess.notna().all() else np.nan,
                     mean_rho_masld=o.rho[dq.set_index("gsm").disease == 1].mean(), mean_rho_healthy=o.rho[dq.set_index("gsm").disease == 0].mean()))
summ = pd.DataFrame(summ)
ok = summ.p_excess.notna()
summ.loc[ok, "bh_excess"] = stats.false_discovery_control(summ.loc[ok, "p_excess"])
summ.to_csv(OUT / "within_donor_summary.tsv", sep="\t", index=False)
de = pd.DataFrame([dict(set=s, beta=dis(dmn.set_index("gsm")[s].reindex(dq.gsm))[0], p_exact=dis(dmn.set_index("gsm")[s].reindex(dq.gsm))[1])
                   for s in list(SETS) + ["CI_minus_CII"]])
de["bh"] = stats.false_discovery_control(de.p_exact)
de.to_csv(OUT / "disease_effect_donor_means.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
print(prim.round(4).to_string()); print(summ.round(3).to_string()); print(de.round(3).to_string())

"""DIG1 job 48 (exploratory; follows job 47). Written before any value for the held-out genes below was looked at.
Job 47's top replicated genes after length control included five NF-kB targets (BIRC3, RELB, NFKB2, TNIP1,
SOD2) and the TWEAK receptor TNFRSF12A. Test whether NF-kB targets not in that list also rise with the
KRT8/SQSTM1 score inside hepatocytes.
Held-out NF-kB set (recalled standard targets [Unverified list]): NFKBIA TNFAIP3 NFKB1 ICAM1 IER3 BCL3 TRAF1
NFKBIE CCL2 CXCL2 TNFAIP2 REL.
Statistic and adjustment as job 47 A2: depth + zonation + long-gene share residualisation; 200 null sets
matched on mean x detection x span tertile (seed 20261024); exact sign-flip over donors; score KRT8 + SQSTM1.
Primary family (Holm over 2): N1 held-out set in GSE298719 Flex (13 donors); N2 in GSE256398 snRNA (17).
Holds = positive and Holm < 0.05.
Descriptive (BH within cohort): NF-kB seen set (BIRC3 RELB NFKB2 TNIP1 SOD2); TNFRSF12A; IL32; complement
(CFH CFHR3 CFHR4 C5 C3); fibrinogen (FGA FGB FGG); AP-1 (JUN JUNB FOS FOSB ATF3); each held-out gene alone;
hepatocyte donor means MASLD vs healthy (+ sex, exact relabeling) for these sets.
Usage: python -I 48_nfkb_heldout.py <new_out_dir>
"""
import gzip
import re
import sys
from itertools import combinations, product
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V5 = ROOT / "Analysis/Hypotheses/V5-20261021T122723Z"
GTF = ROOT / "Analysis/MASLD_Model_Benchmark/executions/gse267145-reference-crosswalk-21065154/Homo_sapiens.GRCh38.98.gtf.gz"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
N_NULL = 200
SCORE = ["KRT8", "SQSTM1"]
ZON = {"central": ["CYP2E1", "GLUL", "CYP1A2", "OAT", "CYP3A4"], "portal": ["CPS1", "HAL", "ASS1", "SDS", "PCK1", "ALDOB"]}
HELD = ["NFKBIA", "TNFAIP3", "NFKB1", "ICAM1", "IER3", "BCL3", "TRAF1", "NFKBIE", "CCL2", "CXCL2", "TNFAIP2", "REL"]
TESTS0 = {"NFkB_heldout": HELD, "NFkB_seen": ["BIRC3", "RELB", "NFKB2", "TNIP1", "SOD2"], "TNFRSF12A": ["TNFRSF12A"], "IL32": ["IL32"],
          "complement": ["CFH", "CFHR3", "CFHR4", "C5", "C3"], "fibrinogen": ["FGA", "FGB", "FGG"], "AP1": ["JUN", "JUNB", "FOS", "FOSB", "ATF3"],
          **{g: [g] for g in HELD}}

# gene lengths
span, exons, names = {}, {}, {}
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split("\t", 9)
        if p[2] not in ("gene", "exon"):
            continue
        gid = re.search(r'gene_id "([^"]+)"', p[8]).group(1)
        a, b = int(p[3]), int(p[4])
        if p[2] == "gene":
            span[gid] = b - a + 1
            m = re.search(r'gene_name "([^"]+)"', p[8]); names[gid] = m.group(1) if m else gid
        else:
            exons.setdefault(gid, []).append((a, b))


def union_len(iv):
    iv = sorted(iv); tot, (s, e) = 0, iv[0]
    for a, b in iv[1:]:
        if a > e:
            tot += e - s + 1; s, e = a, b
        else:
            e = max(e, b)
    return tot + e - s + 1


GL = pd.DataFrame(dict(gene_id=list(span), gene_name=[names[k] for k in span], span=[span[k] for k in span],
                       exon_len=[union_len(exons[k]) if k in exons else np.nan for k in span]))
by_id = GL.set_index("gene_id").span; by_sym = GL.sort_values("span", ascending=False).drop_duplicates("gene_name").set_index("gene_name").span
ex_id = GL.set_index("gene_id").exon_len; ex_sym = GL.sort_values("span", ascending=False).drop_duplicates("gene_name").set_index("gene_name").exon_len


def resid_rank(M, Z):
    R = stats.rankdata(M, axis=0)
    return R - Z @ np.linalg.lstsq(Z, R, rcond=None)[0]


def cor_cols(a, M):
    a = (a - a.mean()) / a.std(); sd = M.std(0); sd[sd == 0] = np.nan
    return ((M - M.mean(0)) / sd * a[:, None]).mean(0)


def signflip(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]; n = len(v)
    S = np.array(list(product([1, -1], repeat=n)), dtype=np.int8)
    return np.mean(np.abs((S * v).mean(1)) >= abs(v.mean()) - 1e-12)


def signflip_mat(E):
    n = E.shape[0]
    S = np.array(list(product([1.0, -1.0], repeat=n)), dtype=np.float32) / n
    obs = np.abs(E.mean(0)); p = np.empty(E.shape[1])
    for a in range(0, E.shape[1], 400):
        b = min(a + 400, E.shape[1])
        p[a:b] = (np.abs(S @ E[:, a:b].astype(np.float32)) >= obs[a:b] - 1e-6).mean(0)
    return p


def cohort(tag, X, gnames, glen, donor, seed):
    rng = np.random.default_rng(seed)
    genes = pd.Index(gnames); g = list(genes)
    umi = np.asarray(X.sum(axis=1)).ravel(); ngene = np.asarray((X > 0).sum(axis=1)).ravel()
    L = X.multiply(1e4 / umi[:, None]).tocsr(); L.data = np.log1p(L.data)
    mean_all = np.asarray(L.mean(axis=0)).ravel(); det_all = np.asarray((X > 0).mean(axis=0)).ravel()
    tests = {k: [x for x in v if x in genes] for k, v in TESTS0.items()}
    tests = {k: v for k, v in tests.items() if v}
    zon = {k: [x for x in v if x in genes] for k, v in ZON.items()}; sco = [x for x in SCORE if x in genes]
    excl = set(x for v in tests.values() for x in v) | set(x for v in zon.values() for x in v) | set(sco)
    keep = np.where((det_all >= 0.01) & np.isfinite(glen))[0]
    kg = np.array(g)[keep]
    binm = pd.qcut(pd.Series(mean_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
    bind = pd.qcut(pd.Series(det_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
    bint = pd.qcut(pd.Series(np.log10(glen[keep])).rank(method="first"), 3, labels=False).to_numpy()
    b0 = pd.Series(binm * 10 + bind, index=kg); b1 = pd.Series((binm * 10 + bind) * 3 + bint, index=kg)
    lq = pd.qcut(pd.Series(np.log10(glen[keep])).rank(method="first"), 5, labels=False).to_numpy()
    long_idx = keep[(lq == 4) & ~np.isin(kg, list(excl))]
    long_share = np.asarray(X[:, long_idx].sum(axis=1)).ravel() / umi
    mt_idx = [i for i, x in enumerate(g) if x.startswith("MT-")]
    mt_share = np.asarray(X[:, mt_idx].sum(axis=1)).ravel() / umi if mt_idx else np.full(len(umi), np.nan)
    det_tab = pd.DataFrame(dict(gene=[x for x in sorted(excl)], span=[glen[genes.get_loc(x)] for x in sorted(excl)],
                                hep_detection=[det_all[genes.get_loc(x)] for x in sorted(excl)]))
    det_tab.to_csv(OUT / f"{tag}_set_gene_detection.tsv", sep="\t", index=False)
    nulls = {}
    for nm, bins in [("A0", b0), ("A1", b1)]:
        pool = {b: [x for x in idx if x not in excl] for b, idx in bins.groupby(bins).groups.items()}
        nulls[nm] = {t: [[rng.choice(pool[bins[x]]) for x in v if x in bins.index] for _ in range(N_NULL)] for t, v in tests.items()}
    need = sorted(excl | set(x for nl in nulls.values() for r in nl.values() for s in r for x in s))
    Ld = L[:, [genes.get_loc(x) for x in need]].toarray(); gi = {x: i for i, x in enumerate(need)}
    rows, cellrho, dmeans = [], [], []
    for d in sorted(np.unique(donor)):
        m = donor == d
        D = Ld[m].astype(np.float64); sd = D.std(0); sd[sd == 0] = np.nan
        Zg = (D - D.mean(0)) / sd
        sc = lambda gs: np.nanmean(Zg[:, [gi[x] for x in gs]], 1)
        mdb = sc(sco)
        base = np.column_stack([np.ones(m.sum()), stats.rankdata(np.log(umi[m])), stats.rankdata(np.log(ngene[m])),
                                stats.rankdata(sc(zon["central"])), stats.rankdata(sc(zon["portal"]))])
        adjs = {"A2": (np.column_stack([base, stats.rankdata(long_share[m])]), "A1")}
        dmeans.append(dict(donor=d, **{t: float(np.nanmean(D[:, [gi[x] for x in v]])) for t, v in tests.items()}))
        cellrho.append(dict(donor=d, rho_score_long_share=stats.spearmanr(mdb, long_share[m])[0],
                            rho_score_mt_share=stats.spearmanr(mdb, mt_share[m])[0] if mt_idx else np.nan,
                            rho_score_umi=stats.spearmanr(mdb, umi[m])[0]))
        for adj, (Zr, nn) in adjs.items():
            ri = resid_rank(mdb[:, None], Zr)[:, 0]
            for t, v in tests.items():
                s = sc(v)
                if not np.isfinite(s).all():
                    continue
                ro = resid_rank(s[:, None], Zr)[:, 0]
                Rn = resid_rank(np.column_stack([np.nan_to_num(sc(r)) for r in nulls[nn][t]]), Zr)
                rows.append(dict(cohort=tag, donor=d, adjust=adj, test=t, rho=np.corrcoef(ri, ro)[0, 1], null_mean=np.nanmean(cor_cols(ri, Rn))))
        print(tag, d, int(m.sum()), flush=True)
    wd = pd.DataFrame(rows); wd["excess"] = wd.rho - wd.null_mean
    summ = []
    for (adj, t), o in wd.groupby(["adjust", "test"]):
        summ.append(dict(cohort=tag, adjust=adj, test=t, n_donors=int(o.excess.notna().sum()), mean_excess=o.excess.mean(),
                         n_excess_positive=int((o.excess > 0).sum()), p_excess=signflip(o.excess)))
    summ = pd.DataFrame(summ)

    return wd, summ, pd.DataFrame(dmeans).assign(cohort=tag)




def disease_effect(dm, meta):
    dm = dm.merge(meta, on="donor"); n = len(dm); dis = dm.disease.to_numpy(float); sx = dm.sex_f.to_numpy(float)
    rel = []
    for ix in combinations(range(n), int(dis.sum())):
        v = np.zeros(n); v[list(ix)] = 1; rel.append(np.linalg.pinv(np.column_stack([np.ones(n), v, sx]))[1])
    P = np.stack(rel); po = np.linalg.pinv(np.column_stack([np.ones(n), dis, sx]))[1]
    sets = [c for c in dm.columns if c in TESTS0]
    Y = dm[sets].to_numpy(float); ok = np.isfinite(Y).all(0)
    bo = np.full(len(sets), np.nan); pv = np.full(len(sets), np.nan)
    bo[ok] = po @ Y[:, ok]; pv[ok] = (np.abs(P @ Y[:, ok]) >= np.abs(bo[ok]) - 1e-12).mean(0)
    return pd.DataFrame(dict(set=sets, beta_masld_minus_healthy=bo, p_exact=pv))


A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str}); dq = dq[dq.included]
A = A[(A.obs.label == "Hepatocytes").to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
fg = list(A.var_names)
wf, sf, df_ = cohort("flex", sparse.csr_matrix(A.X).astype(np.float32), fg, by_sym.reindex(fg).to_numpy(float), A.obs.gsm.astype(str).to_numpy(), 20261024)
del A
H = ad.read_h5ad(RUN / "out_snrna_labels/gse256398_hepatocyte_counts.h5ad")
sg = list(H.var_names); sid = H.var.gene_ids.astype(str).to_numpy()
ws, ss, ds_ = cohort("snrna", sparse.csr_matrix(H.X).astype(np.float32), sg, by_id.reindex(sid).to_numpy(float), H.obs.donor.astype(str).to_numpy(), 20261024)
del H
pd.concat([wf, ws]).to_csv(OUT / "within_donor_rho.tsv.gz", sep="\t", index=False)
summ = pd.concat([sf, ss]); summ["bh"] = summ.groupby("cohort").p_excess.transform(stats.false_discovery_control)
summ.to_csv(OUT / "within_donor_summary.tsv", sep="\t", index=False)
prim = summ[summ.test == "NFkB_heldout"].copy(); prim["id"] = prim.cohort.map({"flex": "N1", "snrna": "N2"})
prim = prim.sort_values("id").reset_index(drop=True)
o = np.argsort(prim.p_excess.to_numpy()); ph = np.maximum.accumulate(prim.p_excess.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph); prim["holds"] = (prim.mean_excess > 0) & (prim.p_holm < 0.05)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)
mf = dq.rename(columns={"gsm": "donor"})[["donor", "disease", "sex_f"]]
v3m = pd.read_csv(ROOT / "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk/donor_metadata.tsv", sep="\t")
ms = v3m.assign(sex_f=(v3m.sex == "female").astype(int))[["donor", "disease", "sex_f"]]
de = pd.concat([disease_effect(df_.drop(columns="cohort"), mf).assign(cohort="flex"), disease_effect(ds_.drop(columns="cohort"), ms).assign(cohort="snrna")])
de.to_csv(OUT / "hep_disease_effect.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400)
print(prim.round(4).to_string())
print(summ.pivot_table(index="test", columns="cohort", values=["mean_excess", "n_excess_positive", "bh"]).round(4).to_string())
print(de.pivot_table(index="set", columns="cohort", values=["beta_masld_minus_healthy", "p_exact"]).round(4).to_string())

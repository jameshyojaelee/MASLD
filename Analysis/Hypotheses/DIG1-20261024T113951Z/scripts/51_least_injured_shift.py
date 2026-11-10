"""DIG1 job 51 (exploratory, post hoc to jobs 46-50). Do the least injured hepatocytes of a MASLD liver
already carry the MASLD hepatocyte shift? Job 50 found that the within-liver injury axis accounts for little
of the bulk disease transcriptome. If the shift is shared by all hepatocytes, it should be present in the
quarter of each liver's hepatocytes with the lowest KRT8+SQSTM1 score.
Cohorts: GSE298719 Flex (13 donors, 6 healthy / 7 MASLD), GSE256398 snRNA (17 donors, 6 / 11; job 45L labels).
Per cell: injury score = mean within-donor z of log1p CP10k KRT8 and SQSTM1, residualised on ranks of log UMI
and log genes within the donor. Strata: within-donor quartiles of that residual (Q1 least injured).
Per donor x stratum: pseudobulk counts, log2 CPM (prior 1). Genes: canonical bulk DEGs (774-unit release,
padj < 0.05, |logFC| > 0.5) detected in >= 1% of the cohort's hepatocytes; KRT8, KRT18, SQSTM1 excluded.
Flex genes matched by symbol, snRNA by Ensembl ID.
Disease effect per gene and stratum = mean(MASLD donors) - mean(healthy donors).
Primary family (Holm over 2), upper-tail exact relabeling of disease over donors (Flex 1,716; snRNA 12,376):
  T1 Flex Q1 sign agreement with the bulk direction.  T2 snRNA Q1 sign agreement.
Descriptive: the same for Q2-Q4 and all hepatocytes; retention = sum(eQ * eAll) / sum(eAll^2); Spearman
with bulk logFC; Q4 - Q1 difference per gene; injury-axis genes (job 47 A2 replicated) vs the rest.
Usage: python -I 51_least_injured_shift.py <new_out_dir>
"""
import sys
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V5 = ROOT / "Analysis/Hypotheses/V5-20261021T122723Z"
V3PB = ROOT / "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk"
BULK = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
SCORE = ["KRT8", "SQSTM1"]
DROP = {"KRT8", "KRT18", "SQSTM1"}
STRATA = ["Q1", "Q2", "Q3", "Q4", "all"]
SHOW = ["GNMT", "MAT1A", "AKR1B10", "INHBE", "ANGPTL8", "FMO1", "CA12", "CIDEC", "TRIM31", "IGF1", "CYP1A1", "CYP1A2", "SLCO1A2",
        "FASN", "SCD", "ME1", "IL32", "TNFRSF12A", "RELB", "SPP1", "PDGFA", "MRAS", "C9", "CFHR4", "HPR", "NUDT13"]

bulk = pd.read_csv(BULK)
bulk["ensembl"] = bulk.gene.str.replace(r"\..*$", "", regex=True)
canon = bulk[(bulk.padj < 0.05) & (bulk.logFC.abs() > 0.5)].copy()
ax = pd.read_csv(RUN / "out_length_confound/genome_A2_cross_cohort.tsv.gz", sep="\t")
axis = dict(zip(ax.gene, ax.replicated.fillna("")))


def resid_rank(M, Z):
    R = np.apply_along_axis(stats.rankdata, 0, M)
    return R - Z @ np.linalg.lstsq(Z, R, rcond=None)[0]


def cohort(tag, X, sym, ens, donor, dis_of):
    X = sparse.csr_matrix(X)
    umi = np.asarray(X.sum(axis=1)).ravel(); ngene = np.asarray((X > 0).sum(axis=1)).ravel()
    det = np.asarray((X > 0).mean(axis=0)).ravel()
    sidx = [list(sym).index(s) for s in SCORE]
    Ls = np.log1p(X[:, sidx].toarray() * (1e4 / umi[:, None]))
    dons = sorted(np.unique(donor)); stratum = np.empty(len(umi), dtype=object)
    for d in dons:
        m = np.where(donor == d)[0]
        S = Ls[m]; sd = S.std(0); sd[sd == 0] = np.nan
        sc = np.nanmean((S - S.mean(0)) / sd, 1)
        Z = np.column_stack([np.ones(len(m)), stats.rankdata(np.log(umi[m])), stats.rankdata(np.log(ngene[m]))])
        r = resid_rank(sc[:, None], Z)[:, 0]
        q = pd.qcut(pd.Series(r).rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"]).astype(str).to_numpy()
        stratum[m] = q
    # gene selection: canonical bulk DEGs measured here
    key = ens if ens is not None else sym
    ckey = canon.ensembl if ens is not None else canon.symbol
    cmap = dict(zip(ckey, canon.logFC))
    gsel = [i for i, k in enumerate(key) if k in cmap and det[i] >= 0.01 and sym[i] not in DROP]
    blfc = np.array([cmap[key[i]] for i in gsel]); gsym = np.array([sym[i] for i in gsel])
    # pseudobulk log2 CPM per donor x stratum
    E = {s: np.zeros((len(dons), len(gsel))) for s in STRATA}
    ncell = []
    for j, d in enumerate(dons):
        for s in STRATA:
            m = (donor == d) if s == "all" else ((donor == d) & (stratum == s))
            tot = np.asarray(X[m].sum(axis=0)).ravel()
            E[s][j] = np.log2(tot[gsel] / tot.sum() * 1e6 + 1)
            ncell.append(dict(cohort=tag, donor=d, stratum=s, n_cells=int(m.sum()), disease=dis_of[d]))
    dis = np.array([dis_of[d] for d in dons]).astype(int)
    n1 = dis.sum(); labs = [np.isin(np.arange(len(dons)), c) for c in combinations(range(len(dons)), n1)]
    L = np.array(labs, dtype=float)  # relabelings x donors (True = MASLD)
    obs_lab = dis.astype(bool)
    rows, gene_rows = [], []
    eff = {}
    for s in STRATA:
        M = E[s]
        # effects for every relabeling: (L/n1 - (1-L)/n0) @ M
        W = L / n1 - (1 - L) / (len(dons) - n1)
        EF = W @ M  # relabelings x genes
        agree = (np.sign(EF) == np.sign(blfc)).mean(1)
        e_obs = M[obs_lab].mean(0) - M[~obs_lab].mean(0); eff[s] = e_obs
        a_obs = (np.sign(e_obs) == np.sign(blfc)).mean()
        rows.append(dict(cohort=tag, stratum=s, n_genes=len(gsel), sign_agree=a_obs, null_mean=agree.mean(),
                         p_upper=(agree >= a_obs - 1e-12).mean(), n_relabelings=len(L),
                         spearman_bulk=stats.spearmanr(e_obs, blfc)[0]))
    R = pd.DataFrame(rows)
    ea = eff["all"]
    R["retention_vs_all"] = [np.sum(eff[s] * ea) / np.sum(ea ** 2) for s in STRATA]
    R["retention_vs_all_bulk_sign_genes"] = [np.sum((eff[s] * ea)[np.sign(ea) == np.sign(blfc)]) / np.sum((ea ** 2)[np.sign(ea) == np.sign(blfc)]) for s in STRATA]
    G = pd.DataFrame({"cohort": tag, "gene": gsym, "bulk_logFC": blfc, "axis": [axis.get(g, "") for g in gsym],
                      **{f"effect_{s}": eff[s] for s in STRATA}})
    G["q4_minus_q1"] = G.effect_Q4 - G.effect_Q1
    return R, G, pd.DataFrame(ncell)


# Flex
A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str}); dq = dq[dq.included]
A = A[(A.obs.label == "Hepatocytes").to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
rf, gf, nf = cohort("flex", A.X, np.array(A.var_names), None, A.obs.gsm.astype(str).to_numpy(), dict(zip(dq.gsm, dq.disease)))
del A
# snRNA
H = ad.read_h5ad(RUN / "out_snrna_labels/gse256398_hepatocyte_counts.h5ad")
md = pd.read_csv(V3PB / "donor_metadata.tsv", sep="\t")
rs, gs, ns = cohort("snrna", H.X, np.array(H.var_names), H.var.gene_ids.astype(str).to_numpy(), H.obs.donor.astype(str).to_numpy(),
                    dict(zip(md.donor.astype(str), md.disease)))
del H

R = pd.concat([rf, rs], ignore_index=True); G = pd.concat([gf, gs], ignore_index=True); N = pd.concat([nf, ns], ignore_index=True)
R.to_csv(OUT / "stratum_sign_agreement.tsv", sep="\t", index=False)
G.to_csv(OUT / "gene_effect_by_stratum.tsv.gz", sep="\t", index=False)
N.to_csv(OUT / "cells_per_donor_stratum.tsv", sep="\t", index=False)
prim = R[R.stratum == "Q1"][["cohort", "n_genes", "sign_agree", "null_mean", "p_upper"]].copy()
prim.insert(0, "id", ["T1", "T2"])
o = np.argsort(prim.p_upper.to_numpy()); ph = np.maximum.accumulate(prim.p_upper.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)
# axis vs non-axis genes
ax_rows = []
for (c, a), o in G.assign(axis_class=np.where(G.axis.isin(["up", "down"]), "axis", "other")).groupby(["cohort", "axis_class"]):
    ea = o.effect_all.to_numpy()
    ax_rows.append(dict(cohort=c, axis_class=a, n=len(o), **{f"retention_{s}": np.sum(o[f"effect_{s}"] * ea) / np.sum(ea ** 2) for s in ["Q1", "Q4"]},
                        median_abs_q4_minus_q1=o.q4_minus_q1.abs().median(),
                        frac_q1_bulk_sign=(np.sign(o.effect_Q1) == np.sign(o.bulk_logFC)).mean()))
pd.DataFrame(ax_rows).to_csv(OUT / "axis_vs_other_retention.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print(prim.round(4).to_string()); print(R.round(3).to_string()); print(pd.DataFrame(ax_rows).round(3).to_string())
print(N.groupby(["cohort", "stratum"]).n_cells.min().to_string())
print(G[G.gene.isin(SHOW)].sort_values(["gene", "cohort"]).round(3).to_string())
for c in ["flex", "snrna"]:
    o = G[G.cohort == c]
    print(c, "largest Q4-Q1 (disease shift concentrated in injured cells)")
    print(o.reindex(o.q4_minus_q1.abs().sort_values(ascending=False).index).head(25).round(3).to_string())

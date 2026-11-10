"""DIG1 job 47 (exploratory robustness check; follows jobs 45 and 46). Written after job 46 showed that the
genes rising with the KRT8/SQSTM1 score in both cohorts are mostly short, abundant cytoplasmic mRNAs
(apolipoproteins, FTL, mtDNA-encoded) and the falling genes are mostly long (CFH, DPYD, CPS1, PPARGC1A).
Before any length analysis was run.
Concern: a cell-level factor that shifts counts from long to short genes (cytoplasmic carry-over into
nuclei in snRNA; cell size or nascent-to-mature RNA balance in whole cells) would raise short chaperone,
CCT, tubulin and actin RNAs and lower long PPARGC1A in KRT8-high cells without any injury programme.
Gene length: Ensembl 98 GTF (GRCh38), genomic span and union exon length; Flex genes matched by symbol,
snRNA by Ensembl ID.
Statistic as job 45 (depth + zonation residualisation, 200 matched random sets, seed 20261024, exact
sign-flip over donors), in GSE298719 Flex (13 donors) and GSE256398 snRNA (17 donors), score KRT8 + SQSTM1.
Adjustments:
  A0 as jobs 39/45 (check that values reproduce).
  A1 null sets matched on hepatocyte mean x detection deciles x genomic-span tertile.
  A2 A1 + rank of each cell's long-gene share (counts in the longest-span quintile of detected genes,
     set/score/zonation genes excluded, over total counts) in the residualisation.
  A3 (snRNA only, descriptive) A2 + rank of the mtDNA-encoded count share.
Primary family (Holm over 8), A2: R1 CCT, R2 HSF1 set, R3 HSPB1, R4 PPARGC1A in each cohort. Holds = same
sign as GSE298719 under A0 and Holm < 0.05.
Descriptive: per-gene excess (job 46 outputs) vs log10 span and exon length, Spearman per cohort; per-cell
rho of the KRT8/SQSTM1 score with long-gene share and mt share; genome-wide scan rerun under A2
(bin-centring on mean x detection x span tertile) with cross-cohort replication counts.
Usage: python -I 47_length_confound.py <new_out_dir>
"""
import gzip
import re
import sys
from itertools import product
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
TESTS0 = {
    "CCT": ["TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8"],
    "HSF1": ["HSPA1A", "HSPA1B", "DNAJB1", "HSPH1", "BAG3", "HSPA4L", "HSP90AA1"],
    "chaperone_network": ["HSPA8", "HSP90AA1", "HSP90AB1", "STIP1", "AHSA1", "CHORDC1", "DNAJA1", "HSPA4", "HSPA4L", "HSPH1"],
    "tubulin": ["TUBA1A", "TUBA1B", "TUBA1C", "TUBA4A", "TUBB", "TUBB2A", "TUBB4B", "TUBB6"],
    "actin": ["ACTB", "ACTG1"],
    "mito_chaperones": ["HSPD1", "HSPE1"],
    "HSPB1": ["HSPB1"], "HSPA8": ["HSPA8"], "HSP90AA1": ["HSP90AA1"], "HSPA1A": ["HSPA1A"], "PPARGC1A": ["PPARGC1A"], "CPT1A": ["CPT1A"],
}
PRIMARY = {"CCT": ("R1", 1), "HSF1": ("R2", 1), "HSPB1": ("R3", 1), "PPARGC1A": ("R4", -1)}

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
GL.to_csv(OUT / "gene_lengths_ensembl98.tsv.gz", sep="\t", index=False)
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
    det_tab.to_csv(OUT / f"{tag}_set_gene_length.tsv", sep="\t", index=False)
    nulls = {}
    for nm, bins in [("A0", b0), ("A1", b1)]:
        pool = {b: [x for x in idx if x not in excl] for b, idx in bins.groupby(bins).groups.items()}
        nulls[nm] = {t: [[rng.choice(pool[bins[x]]) for x in v if x in bins.index] for _ in range(N_NULL)] for t, v in tests.items()}
    need = sorted(excl | set(x for nl in nulls.values() for r in nl.values() for s in r for x in s))
    Ld = L[:, [genes.get_loc(x) for x in need]].toarray(); gi = {x: i for i, x in enumerate(need)}
    rows, cellrho = [], []
    for d in sorted(np.unique(donor)):
        m = donor == d
        D = Ld[m].astype(np.float64); sd = D.std(0); sd[sd == 0] = np.nan
        Zg = (D - D.mean(0)) / sd
        sc = lambda gs: np.nanmean(Zg[:, [gi[x] for x in gs]], 1)
        mdb = sc(sco)
        base = np.column_stack([np.ones(m.sum()), stats.rankdata(np.log(umi[m])), stats.rankdata(np.log(ngene[m])),
                                stats.rankdata(sc(zon["central"])), stats.rankdata(sc(zon["portal"]))])
        adjs = {"A0": (base, "A0"), "A1": (base, "A1"), "A2": (np.column_stack([base, stats.rankdata(long_share[m])]), "A1")}
        if mt_idx:
            adjs["A3"] = (np.column_stack([base, stats.rankdata(long_share[m]), stats.rankdata(mt_share[m])]), "A1")
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

    # genome-wide under A2: per-gene rho, bin-centred on mean x detection x span tertile
    Lk = L[:, keep].tocsr(); is_score = np.isin(kg, sco)
    dons = sorted(np.unique(donor)); EX = np.full((len(dons), len(keep)), np.nan)
    sidx = [genes.get_loc(x) for x in sco]; Ls = L[:, sidx].toarray()
    for i, d in enumerate(dons):
        m = np.where(donor == d)[0]
        S = Ls[m]; sd = S.std(0); sd[sd == 0] = np.nan
        scv = np.nanmean((S - S.mean(0)) / sd, 1)
        Zc = Ld[m][:, [gi[x] for x in zon["central"]]]; Zp = Ld[m][:, [gi[x] for x in zon["portal"]]]
        zc_ = lambda Q: np.nanmean((Q - Q.mean(0)) / np.where(Q.std(0) == 0, np.nan, Q.std(0)), 1)
        Z = np.column_stack([np.ones(len(m)), stats.rankdata(np.log(umi[m])), stats.rankdata(np.log(ngene[m])),
                             stats.rankdata(zc_(Zc)), stats.rankdata(zc_(Zp)), stats.rankdata(long_share[m])])
        r_s = resid_rank(scv[:, None], Z)[:, 0]; r_s = (r_s - r_s.mean()) / r_s.std()
        rho = np.full(len(keep), np.nan); Lm = Lk[m].tocsc()
        for a in range(0, len(keep), 2000):
            b = min(a + 2000, len(keep))
            R = resid_rank(Lm[:, a:b].toarray().astype(np.float64), Z); sdr = R.std(0); sdr[sdr < 1e-9] = np.nan
            rho[a:b] = ((R - R.mean(0)) / sdr * r_s[:, None]).mean(0)
        bm = pd.Series(rho[~is_score]).groupby(b1.to_numpy()[~is_score]).mean()
        EX[i] = rho - bm.reindex(b1.to_numpy()).to_numpy()
    ok = np.isfinite(EX).all(0) & ~is_score
    p = np.full(len(keep), np.nan); p[ok] = signflip_mat(EX[:, ok])
    gw = pd.DataFrame(dict(gene=kg, span=glen[keep], mean_excess_A2=np.nanmean(EX, 0), n_pos=(EX > 0).sum(0), p=p))
    gw["bh"] = np.nan; gw.loc[ok, "bh"] = stats.false_discovery_control(p[ok])
    return wd, summ, pd.DataFrame(cellrho).assign(cohort=tag), gw


# Flex
A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str}); dq = dq[dq.included]
A = A[(A.obs.label == "Hepatocytes").to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
fg = list(A.var_names)
wf, sf, cf, gf = cohort("flex", sparse.csr_matrix(A.X).astype(np.float32), fg, by_sym.reindex(fg).to_numpy(float), A.obs.gsm.astype(str).to_numpy(), 20261024)
del A
# snRNA
H = ad.read_h5ad(RUN / "out_snrna_labels/gse256398_hepatocyte_counts.h5ad")
sg = list(H.var_names); sid = H.var.gene_ids.astype(str).to_numpy()
ws, ss, cs, gs = cohort("snrna", sparse.csr_matrix(H.X).astype(np.float32), sg, by_id.reindex(sid).to_numpy(float), H.obs.donor.astype(str).to_numpy(), 20261024)
del H

pd.concat([wf, ws]).to_csv(OUT / "within_donor_rho.tsv.gz", sep="\t", index=False)
summ = pd.concat([sf, ss]); summ["bh"] = summ.groupby(["cohort", "adjust"]).p_excess.transform(stats.false_discovery_control)
summ.to_csv(OUT / "within_donor_summary.tsv", sep="\t", index=False)
pd.concat([cf, cs]).to_csv(OUT / "score_vs_cell_composition.tsv", sep="\t", index=False)
prim = summ[(summ.adjust == "A2") & summ.test.isin(list(PRIMARY))].copy()
prim["id"] = prim.cohort + "_" + prim.test.map(lambda t: PRIMARY[t][0]); prim["expected_sign"] = prim.test.map(lambda t: PRIMARY[t][1])
prim = prim.sort_values("id").reset_index(drop=True)
o = np.argsort(prim.p_excess.to_numpy()); ph = np.maximum.accumulate(prim.p_excess.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim["holds"] = (np.sign(prim.mean_excess) == prim.expected_sign) & (prim.p_holm < 0.05)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)

# per-gene excess (job 46) vs length
J = pd.read_csv(RUN / "out_injury_coupling_genome/cross_cohort_gene_coupling.tsv.gz", sep="\t")
sn46 = pd.read_csv(RUN / "out_injury_coupling_genome/snrna_gene_coupling.tsv.gz", sep="\t")
fl46 = pd.read_csv(RUN / "out_injury_coupling_genome/flex_gene_coupling.tsv.gz", sep="\t")
sn46["span"] = sn46.gene.map(dict(zip(sg, by_id.reindex(sid).to_numpy(float)))); sn46["exon_len"] = sn46.gene.map(dict(zip(sg, ex_id.reindex(sid).to_numpy(float))))
fl46["span"] = by_sym.reindex(fl46.gene).to_numpy(float); fl46["exon_len"] = ex_sym.reindex(fl46.gene).to_numpy(float)
lr = []
for tag, T in [("flex", fl46), ("snrna", sn46)]:
    T = T[np.isfinite(T.span) & T.mean_excess.notna()]
    lr.append(dict(cohort=tag, n=len(T), rho_excess_log_span=stats.spearmanr(T.mean_excess, np.log10(T.span))[0],
                   rho_excess_log_exon_len=stats.spearmanr(T.mean_excess, np.log10(T.exon_len), nan_policy="omit")[0],
                   rho_excess_hep_mean=stats.spearmanr(T.mean_excess, T.hep_mean)[0]))
J["span"] = by_sym.reindex(J.gene).to_numpy(float)
for k in ["up", "down"]:
    s = J[J.replicated == k]
    lr.append(dict(cohort=f"job46 replicated {k}", n=len(s), median_span=float(np.nanmedian(s.span))))
lr.append(dict(cohort="job46 all shared", n=len(J), median_span=float(np.nanmedian(J.span))))
pd.DataFrame(lr).to_csv(OUT / "excess_vs_length.tsv", sep="\t", index=False)

# genome-wide under A2: cross-cohort
G2 = gf.merge(gs, on="gene", suffixes=("_flex", "_sn"))
both = (G2.bh_flex < 0.05) & (G2.bh_sn < 0.05)
G2["replicated"] = np.where(both & (G2.mean_excess_A2_flex > 0) & (G2.mean_excess_A2_sn > 0), "up",
                            np.where(both & (G2.mean_excess_A2_flex < 0) & (G2.mean_excess_A2_sn < 0), "down", np.where(both, "opposite", "")))
G2.to_csv(OUT / "genome_A2_cross_cohort.tsv.gz", sep="\t", index=False)
gsum = dict(n_shared=len(G2), spearman_excess=stats.spearmanr(G2.mean_excess_A2_flex, G2.mean_excess_A2_sn, nan_policy="omit")[0],
            n_up=int((G2.replicated == "up").sum()), n_down=int((G2.replicated == "down").sum()), n_opposite=int((G2.replicated == "opposite").sum()),
            rho_len_flex=stats.spearmanr(G2.mean_excess_A2_flex, np.log10(G2.span_flex), nan_policy="omit")[0],
            rho_len_sn=stats.spearmanr(G2.mean_excess_A2_sn, np.log10(G2.span_sn), nan_policy="omit")[0])
pd.DataFrame([gsum]).to_csv(OUT / "genome_A2_summary.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400)
print(prim.round(4).to_string())
print(summ.pivot_table(index="test", columns=["cohort", "adjust"], values="mean_excess").round(3).to_string())
print(summ.pivot_table(index="test", columns=["cohort", "adjust"], values="n_excess_positive").to_string())
print(pd.DataFrame(lr).round(3).to_string())
print(pd.concat([cf, cs]).groupby("cohort")[["rho_score_long_share", "rho_score_mt_share", "rho_score_umi"]].mean().round(3).to_string())
print(pd.DataFrame([gsum]).round(3).to_string())
cols = ["gene", "span_flex", "mean_excess_A2_flex", "bh_flex", "mean_excess_A2_sn", "bh_sn"]
print(G2[G2.replicated == "up"].assign(z=lambda x: x.mean_excess_A2_flex / G2.mean_excess_A2_flex.std() + x.mean_excess_A2_sn / G2.mean_excess_A2_sn.std()).sort_values("z", ascending=False)[cols].head(60).round(3).to_string())
print(G2[G2.replicated == "down"].assign(z=lambda x: x.mean_excess_A2_flex / G2.mean_excess_A2_flex.std() + x.mean_excess_A2_sn / G2.mean_excess_A2_sn.std()).sort_values("z")[cols].head(40).round(3).to_string())

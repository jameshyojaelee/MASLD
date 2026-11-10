"""DIG1 job 46 (exploratory discovery scan; follows jobs 34-39, 45). Written before any genome-wide
coupling was computed.
Question: across all genes, which RNAs rise or fall with the KRT8/SQSTM1 injury score inside the
hepatocytes of one liver, in two cohorts with different assays (GSE298719 Flex whole cells, 13 donors;
GSE256398 snRNA nuclei, 17 donors)? Which of them change in liver protein with ballooning?
Injury score: mean within-donor z of log1p CP10k of KRT8 and SQSTM1 (KRT18 is not on the Flex panel;
the snRNA score adds KRT18 as a sensitivity).
Per donor, per gene (pooled hepatocyte detection >= 1% in that cohort): rho = Pearson of ranks of the
gene and of the score, both residualised on ranks of log UMI and log genes. Excess = rho minus the mean
rho of the other genes in the same hepatocyte mean x detection decile bin (score genes excluded). This
is the per-gene form of the matched-set null of job 39, re-derived inside each cohort.
Across donors: mean excess, number of donors > 0, exact sign-flip p (8,192 and 131,072 patterns); BH
over genes within each cohort. Replicated = BH < 0.05 in both cohorts with the same sign.
GSE256398 labels: V3 route through V0's build_labelled, computed once by job 45L (out_snrna_labels/).
Descriptive, not tests (genes are correlated, so gene-level p values are anti-conservative):
  cross-cohort Spearman of excess; liver protein (PXD051911) and bulk RNA (GSE130970) ballooning rho and
  MASLD-vs-control protein logFC for replicated genes (jobs 19, 25); frozen 117-program member means
  against 1,000 size-matched random gene sets (seed 20261024).
Usage: python -I 46_injury_coupling_genome.py <new_out_dir>
"""
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
MEMB = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
rng = np.random.default_rng(20261024)
SCORE = ["KRT8", "SQSTM1"]


def rank_resid(M, Z):
    R = stats.rankdata(M, axis=0)
    return R - Z @ np.linalg.lstsq(Z, R, rcond=None)[0]


def signflip_p(E):
    """E: donors x genes, no NaN. Exact two-sided sign-flip p per gene, chunked over genes."""
    n = E.shape[0]
    S = np.array(list(product([1.0, -1.0], repeat=n)), dtype=np.float32) / n
    obs = np.abs(E.mean(0))
    p = np.empty(E.shape[1])
    for a in range(0, E.shape[1], 400):
        b = min(a + 400, E.shape[1])
        p[a:b] = (np.abs(S @ E[:, a:b].astype(np.float32)) >= obs[a:b] - 1e-6).mean(0)
    return p


def scan(X, gnames, donor, score_genes, tag):
    """X: hepatocytes x genes counts (csr). Returns per-gene summary and per-donor excess."""
    genes = pd.Index(gnames)
    umi = np.asarray(X.sum(axis=1)).ravel(); ngene = np.asarray((X > 0).sum(axis=1)).ravel()
    L = X.multiply(1e4 / umi[:, None]).tocsr(); L.data = np.log1p(L.data)
    mean_all = np.asarray(L.mean(axis=0)).ravel(); det_all = np.asarray((X > 0).mean(axis=0)).ravel()
    keep = np.where(det_all >= 0.01)[0]
    score_genes = [x for x in score_genes if x in genes]
    sidx = [genes.get_loc(x) for x in score_genes]
    binm = pd.qcut(pd.Series(mean_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
    bind = pd.qcut(pd.Series(det_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
    binid = binm * 10 + bind
    is_score = np.isin(keep, sidx)
    Lk = L[:, keep].tocsr(); Ls = L[:, sidx].toarray()
    donors = sorted(np.unique(donor))
    EX = np.full((len(donors), len(keep)), np.nan); RH = np.full_like(EX, np.nan)
    for i, d in enumerate(donors):
        m = np.where(donor == d)[0]
        S = Ls[m]; sd = S.std(0); sd[sd == 0] = np.nan
        sc = np.nanmean((S - S.mean(0)) / sd, 1)
        Z = np.column_stack([np.ones(len(m)), stats.rankdata(np.log(umi[m])), stats.rankdata(np.log(ngene[m]))])
        r_s = rank_resid(sc[:, None], Z)[:, 0]; r_s = (r_s - r_s.mean()) / r_s.std()
        rho = np.full(len(keep), np.nan); Lm = Lk[m].tocsc()
        for a in range(0, len(keep), 2000):
            b = min(a + 2000, len(keep))
            D = Lm[:, a:b].toarray().astype(np.float64)
            R = rank_resid(D, Z); sdr = R.std(0); sdr[sdr < 1e-9] = np.nan
            rho[a:b] = ((R - R.mean(0)) / sdr * r_s[:, None]).mean(0)
        bmean = pd.Series(rho[~is_score]).groupby(binid[~is_score]).mean()
        RH[i] = rho; EX[i] = rho - bmean.reindex(binid).to_numpy()
        print(tag, d, len(m), flush=True)
    ok = np.isfinite(EX).all(0) & ~is_score
    p = np.full(len(keep), np.nan); p[ok] = signflip_p(EX[:, ok])
    res = pd.DataFrame(dict(gene=genes[keep], hep_mean=mean_all[keep], hep_detection=det_all[keep], mean_rho=np.nanmean(RH, 0),
                            mean_excess=np.nanmean(EX, 0), n_pos=(EX > 0).sum(0), n_donors=np.isfinite(EX).sum(0), p_signflip=p))
    res["bh"] = np.nan; res.loc[ok, "bh"] = stats.false_discovery_control(p[ok])
    return res, pd.DataFrame(EX, index=donors, columns=genes[keep])


# GSE298719 Flex
A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str}); dq = dq[dq.included]
A = A[(A.obs.label == "Hepatocytes").to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
flex, flex_ex = scan(sparse.csr_matrix(A.X).astype(np.float32), list(A.var_names), A.obs.gsm.astype(str).to_numpy(), SCORE, "flex")
flex.to_csv(OUT / "flex_gene_coupling.tsv.gz", sep="\t", index=False)
flex_ex.to_csv(OUT / "flex_donor_excess.tsv.gz", sep="\t")
del A

# GSE256398 snRNA
hs = ad.read_h5ad(RUN / "out_snrna_labels/gse256398_hepatocyte_counts.h5ad")
Xs = sparse.csr_matrix(hs.X).astype(np.float32); gs = list(hs.var_names); dH = hs.obs["donor"].astype(str).to_numpy()
assert len(np.unique(dH)) == 17
del hs
sn, sn_ex = scan(Xs, gs, dH, SCORE, "snrna")
sn.to_csv(OUT / "snrna_gene_coupling.tsv.gz", sep="\t", index=False)
sn_ex.to_csv(OUT / "snrna_donor_excess.tsv.gz", sep="\t")
sn18, _ = scan(Xs, gs, dH, SCORE + ["KRT18"], "snrna_krt18")
sn18.to_csv(OUT / "snrna_gene_coupling_with_krt18.tsv.gz", sep="\t", index=False)
del Xs

# cross-cohort
T = flex.merge(sn, on="gene", suffixes=("_flex", "_sn"))
T["sens_krt18_excess"] = T.gene.map(sn18.set_index("gene").mean_excess)
both = (T.bh_flex < 0.05) & (T.bh_sn < 0.05)
T["replicated"] = np.where(both & (T.mean_excess_flex > 0) & (T.mean_excess_sn > 0), "up",
                           np.where(both & (T.mean_excess_flex < 0) & (T.mean_excess_sn < 0), "down",
                                    np.where(both, "opposite", "")))
T["z_avg"] = (T.mean_excess_flex / T.mean_excess_flex.std() + T.mean_excess_sn / T.mean_excess_sn.std()) / 2
bal = pd.read_csv(RUN / "out_bulk_rna_ballooning/per_gene_rna_vs_protein_ballooning.tsv.gz", sep="\t")[["gene", "rho_prot", "rho_rna"]]
dis = pd.read_csv(RUN / "out_prot_rna_discord/gene_discordance_disease.tsv.gz", sep="\t")[["gene", "logFC_prot", "logFC_rna", "d"]]
T = T.merge(bal, on="gene", how="left").merge(dis, on="gene", how="left")
T.to_csv(OUT / "cross_cohort_gene_coupling.tsv.gz", sep="\t", index=False)
cc = dict(n_shared=len(T), spearman_excess=stats.spearmanr(T.mean_excess_flex, T.mean_excess_sn)[0],
          n_bh_flex=int((T.bh_flex < 0.05).sum()), n_bh_sn=int((T.bh_sn < 0.05).sum()),
          n_up=int((T.replicated == "up").sum()), n_down=int((T.replicated == "down").sum()), n_opposite=int((T.replicated == "opposite").sum()),
          spearman_sn_vs_krt18_sensitivity=stats.spearmanr(T.mean_excess_sn, T.sens_krt18_excess, nan_policy="omit")[0])
pd.DataFrame([cc]).to_csv(OUT / "cross_cohort_summary.tsv", sep="\t", index=False)
mod = []
for k in ["up", "down", "opposite", ""]:
    s = T[T.replicated == k]
    mod.append(dict(group=k or "not replicated", n=len(s), n_protein=int(s.rho_prot.notna().sum()), mean_rho_prot_ballooning=s.rho_prot.mean(),
                    frac_prot_neg=(s.rho_prot < 0).mean(), mean_rho_bulk_rna_ballooning=s.rho_rna.mean(), frac_rna_pos=(s.rho_rna > 0).mean(),
                    mean_logFC_prot=s.logFC_prot.mean(), mean_logFC_rna=s.logFC_rna.mean(), mean_d=s.d.mean()))
pd.DataFrame(mod).to_csv(OUT / "replicated_vs_protein_bulk.tsv", sep="\t", index=False)

# frozen programs: member weighted mean of z_avg against size-matched random sets
mm = pd.read_csv(MEMB / "program_membership_v2.tsv", sep="\t")
mm = mm[mm.mapped_symbol_status == "gencode_v49_unique_symbol_confirmed"]
reg = pd.read_csv(MEMB / "program_registry_v2.tsv", sep="\t")[["program_uid", "cell_type", "module_name"]]
zs = T.set_index("gene").z_avg
prog = []
for uid, o in mm.groupby("program_uid"):
    o = o[o.mapped_symbol.isin(zs.index)]
    if len(o) < 5:
        continue
    w = o.original_l1_weight.to_numpy(); v = zs.loc[o.mapped_symbol].to_numpy()
    obs = float(np.sum(w * v) / w.sum())
    null = np.array([np.sum(w * zs.to_numpy()[rng.choice(len(zs), len(o), replace=False)]) / w.sum() for _ in range(1000)])
    prog.append(dict(program_uid=uid, n_genes_measured=len(o), weight_measured=float(w.sum()), weighted_mean_z=obs,
                     null_sd=null.std(), z_vs_null=(obs - null.mean()) / null.std(), p_two_sided=(1 + (np.abs(null - null.mean()) >= abs(obs - null.mean())).sum()) / 1001))
prog = pd.DataFrame(prog).merge(reg, on="program_uid", how="left")
prog["bh"] = stats.false_discovery_control(prog.p_two_sided)
prog.sort_values("z_vs_null").to_csv(OUT / "frozen_program_injury_coupling.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400)
print(pd.DataFrame([cc]).round(3).to_string()); print(pd.DataFrame(mod).round(3).to_string())
cols = ["gene", "mean_excess_flex", "n_pos_flex", "bh_flex", "mean_excess_sn", "n_pos_sn", "bh_sn", "rho_prot", "rho_rna", "logFC_prot"]
print(T[T.replicated == "up"].sort_values("z_avg", ascending=False)[cols].head(80).round(3).to_string())
print(T[T.replicated == "down"].sort_values("z_avg")[cols].head(80).round(3).to_string())
print(prog.sort_values("z_vs_null")[["cell_type", "module_name", "n_genes_measured", "weighted_mean_z", "z_vs_null", "bh"]].round(3).to_string())

"""DIG1 job 15 (exploratory). GSE298719 Flex, 13 donors. Two within-donor partner screens, same
method as job 9 Q2: per donor, partial Spearman of the anchor gene (CP-UMI) with every gene
detected in >= 10% of the lineage's cells, log UMI removed from ranks; Fisher z averaged over
donors with >= 50 cells; one-sample t; BH over genes.
  A. FGF7 in fibroblasts: which fibroblast state makes FGF7? (Takase et al. 2013: Thy1+ cells.)
  B. FGFR2 in hepatocytes: does FGFR2 travel with SOX9/EPCAM plasticity genes?
Also per-donor FGF7 detection and mean in fibroblasts (descriptive). Usage: python -I 15_fgf7_partners.py
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V5 = ROOT / "Analysis/Hypotheses/V5-20261021T122723Z"
OUT = RUN / "out_fgf7_partners"
OUT.mkdir(exist_ok=False)
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str})
dq = dq[dq.included].sort_values("gsm").reset_index(drop=True)
A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad", backed="r")
genes = pd.Index(A.var_names)


def partners(label, anchor, min_cells=50):
    m = (A.obs.label == label).to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()
    S = A[m].to_memory()
    X = S.X.tocsr().astype(np.float64)
    gsm = S.obs.gsm.astype(str).to_numpy()
    umi = np.asarray(X.sum(axis=1)).ravel()
    det = np.asarray((X > 0).mean(axis=0)).ravel()
    a = genes.get_loc(anchor)
    cand = np.where(det >= 0.10)[0]
    cand = cand[cand != a]
    donors = [g for g in dq.gsm if (gsm == g).sum() >= min_cells]
    Z = []
    for g in donors:
        ix = np.where(gsm == g)[0]
        C = np.column_stack([np.ones(len(ix)), np.log(umi[ix])])

        def rres(M):
            Rk = stats.rankdata(M, axis=0)
            return Rk - C @ np.linalg.lstsq(C, Rk, rcond=None)[0]
        ry = rres(np.asarray(X[ix, a].todense()).ravel() / umi[ix])
        rho = np.empty(len(cand))
        for s in range(0, len(cand), 500):
            cols = cand[s:s + 500]
            rx = rres(X[ix][:, cols].toarray() / umi[ix][:, None])
            rho[s:s + 500] = (rx * ry[:, None]).sum(0) / np.sqrt((rx ** 2).sum(0) * (ry ** 2).sum())
        Z.append(np.arctanh(np.clip(rho, -0.999, 0.999)))
    Z = np.array(Z)
    dis = dq.set_index("gsm").disease.reindex(donors).to_numpy()
    t = stats.ttest_1samp(Z, 0)
    res = pd.DataFrame(dict(gene=genes[cand], detection=det[cand], mean_z=Z.mean(0), n_pos=(Z > 0).sum(0), n_donors=len(donors),
                            mean_z_healthy=Z[dis == 0].mean(0), mean_z_masld=Z[dis == 1].mean(0), t=t.statistic, p=t.pvalue))
    res["bh_q"] = stats.false_discovery_control(res.p.fillna(1))
    res.sort_values("mean_z", ascending=False).to_csv(OUT / f"{anchor}_partners_in_{label}.tsv.gz", sep="\t", index=False)
    per = pd.DataFrame(dict(gsm=gsm, det=np.asarray(X[:, a].todense()).ravel() > 0,
                            cp10k=np.asarray(X[:, a].todense()).ravel() / umi * 1e4)).groupby("gsm").agg(
        n_cells=("det", "size"), frac_detected=("det", "mean"), mean_cp10k=("cp10k", "mean"))
    per = per.join(dq.set_index("gsm")[["disease", "fibrosis"]])
    per.to_csv(OUT / f"{anchor}_per_donor_in_{label}.tsv", sep="\t")
    return res, per, det[a]


SHOW = {"FGF7": ["THY1", "PDGFRA", "PDGFRB", "RGS5", "LRAT", "RELN", "HGF", "COL1A1", "COL3A1", "ACTA2", "LUM", "DCN",
                 "CCL19", "CXCL12", "NGFR", "ELN", "MGP", "PI16", "CD34", "SLIT2", "MSLN", "UPK3B", "FBLN1", "PTGDS"],
        "FGFR2": ["SOX9", "EPCAM", "KRT7", "KRT19", "KRT8", "KRT18", "CD24", "SPP1", "AKR1B10", "CDKN1A", "GNMT", "MAT1A",
                  "ITGA6", "TNFRSF12A", "ANXA4", "CFTR", "ALB", "CYP2E1"]}
for label, anchor in [("Fibroblasts", "FGF7"), ("Hepatocytes", "FGFR2")]:
    res, per, d = partners(label, anchor)
    print(f"== {anchor} in {label}: anchor detection {d:.3f}; {int((res.bh_q < 0.05).sum())}/{len(res)} genes BH < 0.05", flush=True)
    print(per.round(3).to_string())
    r = res.set_index("gene")
    print(r.reindex(SHOW[anchor])[["detection", "mean_z", "n_pos", "mean_z_healthy", "mean_z_masld", "bh_q"]].round(3).to_string())
    print("top +", ", ".join(res[res.bh_q < 0.05].head(30).gene))
    print("top -", ", ".join(res[res.bh_q < 0.05].sort_values("mean_z").head(30).gene), flush=True)

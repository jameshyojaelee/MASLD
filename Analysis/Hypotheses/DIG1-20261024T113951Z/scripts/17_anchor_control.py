"""DIG1 job 17 (exploratory). Control for jobs 9 and 15: do the within-donor partners of ANY
anchor gene predict the disease change, or only those of genes on the disease axis?
GSE298719 Flex; lineages Hepatocytes (anchor GNMT) and Fibroblasts (anchor FGF7). Anchors:
  null   = 40 genes with |logFC| < 0.15 and edgeR FDR > 0.5 in that lineage (job 1 v2 table),
           detection within 0.5x-2x of the real anchor's detection (seed 20261024);
  changed= 40 genes with FDR < 0.05, same detection window;
  plus the real anchor.
Partner z: same statistic as job 9 (per donor, ranks of CP-UMI residualised on log UMI, Pearson
of residual ranks, Fisher z, mean over donors), from all donors and from healthy donors only,
over genes detected in >= 10% of the lineage's cells. Read-out per anchor: Spearman of partner z
(anchor excluded) with GSE298719 logFC and with independent GSE256398 logFC; reported against the
anchor's own logFC. Descriptive. Usage: python -I 17_anchor_control.py
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V5 = ROOT / "Analysis/Hypotheses/V5-20261021T122723Z"
OUT = RUN / "out_anchor_control"
OUT.mkdir(exist_ok=False)
rng = np.random.default_rng(20261024)
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str})
dq = dq[dq.included].sort_values("gsm").reset_index(drop=True)
A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad", backed="r")
genes = pd.Index(A.var_names)
de = pd.read_csv(RUN / "out_pseudobulk_v2/gene_DE_edgeR_all_lineages.tsv.gz", sep="\t").drop_duplicates(["lkey", "symbol"])
rows = []
for label, key, anchor, vkey in [("Hepatocytes", "GSE298719:hepatocytes", "GNMT", "GSE256398:hepatocytes"),
                                 ("Fibroblasts", "GSE298719:fibroblasts", "FGF7", "GSE256398:fibroblasts")]:
    m = (A.obs.label == label).to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()
    S = A[m].to_memory()
    X = S.X.tocsr().astype(np.float64)
    gsm = S.obs.gsm.astype(str).to_numpy()
    umi = np.asarray(X.sum(axis=1)).ravel()
    det = pd.Series(np.asarray((X > 0).mean(axis=0)).ravel(), index=genes)
    cand = np.where(det.to_numpy() >= 0.10)[0]
    F = de[de.lkey == key].set_index("symbol")
    V = de[de.lkey == vkey].set_index("symbol")
    d0 = det[anchor]
    win = det[(det >= 0.5 * d0) & (det <= min(2 * d0, 1.0)) & det.index.isin(F.index)].index
    f = F.reindex(win)
    null = f[(f.logFC.abs() < 0.15) & (f.FDR > 0.5)].index.difference([anchor])
    chg = f[f.FDR < 0.05].index.difference([anchor])
    anchors = [anchor] + list(rng.choice(null, min(40, len(null)), replace=False)) + list(rng.choice(chg, min(40, len(chg)), replace=False))
    kind = ["real"] + ["null"] * min(40, len(null)) + ["changed"] * min(40, len(chg))
    a_idx = np.array([genes.get_loc(g) for g in anchors])
    cols = np.union1d(cand, a_idx)
    pos = {c: i for i, c in enumerate(cols)}
    Zall, dis = [], []
    for g in dq.gsm:
        ix = np.where(gsm == g)[0]
        if len(ix) < 50:
            continue
        C = np.column_stack([np.ones(len(ix)), np.log(umi[ix])])
        M = X[ix][:, cols].toarray() / umi[ix][:, None]
        Rk = stats.rankdata(M, axis=0)
        del M
        Rk -= C @ np.linalg.lstsq(C, Rk, rcond=None)[0]
        Rk /= np.sqrt((Rk ** 2).sum(0)) + 1e-12
        r = Rk.T @ Rk[:, [pos[c] for c in a_idx]]  # genes x anchors
        Zall.append(np.arctanh(np.clip(r, -0.999, 0.999)))
        dis.append(int(dq.set_index("gsm").disease[g]))
        del Rk
    Zall, dis = np.array(Zall), np.array(dis)
    zg = pd.Index(genes[cols])
    for k, (ga, kd) in enumerate(zip(anchors, kind)):
        for sub, msk in [("all", np.ones(len(dis), bool)), ("healthy", dis == 0)]:
            z = pd.Series(Zall[msk][:, :, k].mean(0), index=zg).drop(ga)
            z = z[z.index.isin(cand_names := genes[cand])]
            fz = F.logFC.reindex(z.index); vz = V.logFC.reindex(z.index)
            ok_f, ok_v = fz.notna(), vz.notna()
            rows.append(dict(lineage=label, anchor=ga, kind=kd, donors=sub, anchor_det=det[ga],
                             anchor_logFC_F=F.logFC.get(ga), anchor_FDR_F=F.FDR.get(ga), anchor_logFC_V=V.logFC.get(ga),
                             rho_vs_F=stats.spearmanr(z[ok_f], fz[ok_f])[0], rho_vs_V=stats.spearmanr(z[ok_v], vz[ok_v])[0],
                             n_F=int(ok_f.sum()), n_V=int(ok_v.sum())))
    print(label, "done", flush=True)
res = pd.DataFrame(rows)
res.to_csv(OUT / "anchor_control.tsv", sep="\t", index=False)
for (l, d), g in res.groupby(["lineage", "donors"]):
    print(f"== {l} / {d}")
    print(g.groupby("kind")[["rho_vs_F", "rho_vs_V"]].describe(percentiles=[0.05, 0.5, 0.95]).round(3).to_string())
    print(g[g.kind == "real"][["anchor", "anchor_logFC_F", "rho_vs_F", "rho_vs_V"]].round(3).to_string(index=False))
    gc = g[g.kind == "changed"]
    print("changed anchors: spearman(anchor logFC, rho_vs_V) =", round(stats.spearmanr(gc.anchor_logFC_F, gc.rho_vs_V)[0], 3))

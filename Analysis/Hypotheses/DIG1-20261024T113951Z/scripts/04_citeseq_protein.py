"""DIG1 job 4 (exploratory; EXPLORATORY_PLAN.md question F). GSE344087 CITE-seq.

Which surface proteins track fixed program scores inside a lineage, within donors? QC and
CellTypist labels repeat V8 01_gate.py. Program scores use the V8 testable programs (detection
>= 1%, >= 8 genes, >= 0.80 weight in pooled lineage cells). ADT: per-cell CLR over the 140
antibodies. Both program score and ADT are residualised on log10 RNA UMI and log10 ADT total
within donor before Spearman, so depth cannot drive the pairing. Donor is the unit: per-donor
rho, then a one-sample t on Fisher z across donors with >= 50 lineage cells; BH over all
program x ADT pairs within a lineage. Isotype controls are reported as the background floor.
Usage: python -I 04_citeseq_protein.py
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V8 = ROOT / "Analysis/Hypotheses/V8-20261021T235438Z"
DL = V8 / "downloads_untrusted"
HOT = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
CT_MODEL = Path("/gpfs/commons/home/jameslee/.celltypist/data/models/Healthy_Human_Liver.pkl")
LINEAGE_LABEL = {"hepatocytes": "Hepatocytes", "fibroblasts": "Fibroblasts", "cholangiocytes": "Cholangiocytes",
                 "macrophages": "Macrophages", "tcells": "T cells"}
SEED, MIN_CELLS = 20260923, 50
OUT = RUN / "out_citeseq"
OUT.mkdir(exist_ok=False)

adatas, adts = [], []
for f in sorted((DL / "raw").glob("GSM*_filtered_feature_bc_matrix.h5")):
    gsm = f.name.split("_")[0]
    a = sc.read_10x_h5(f, gex_only=False)
    a.var_names_make_unique()
    a.obs_names = [f"{gsm}_{b}" for b in a.obs_names]
    a.obs["gsm"] = gsm
    is_ab = (a.var.feature_types == "Antibody Capture").to_numpy()
    ab = a[:, is_ab].copy()
    a = a[:, ~is_ab].copy()
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
    keep = ((a.obs.n_genes_by_counts >= 200) & (a.obs.pct_counts_mt < 25)).to_numpy()
    adatas.append(a[keep].copy()); adts.append(ab[keep].copy())
A = ad.concat(adatas, join="inner")
B = ad.concat(adts, join="inner")
assert (A.obs_names == B.obs_names).all()
A.X = sparse.csr_matrix(A.X)
counts = A.X.copy()
rna_umi = np.asarray(counts.sum(axis=1)).ravel().astype(float)
sc.pp.normalize_total(A, target_sum=1e4)
sc.pp.log1p(A)
import celltypist  # noqa: E402

np.random.seed(SEED)
pred = celltypist.annotate(A, model=celltypist.models.Model.load(model=str(CT_MODEL)), majority_voting=True)
A.obs["label"] = pred.predicted_labels["majority_voting"].astype(str).values
v8 = pd.read_csv(V8 / "outputs/lineage_counts_per_donor.tsv", sep="\t", index_col=0)
now = pd.crosstab(A.obs.gsm, A.obs.label)
(now.reindex_like(v8).fillna(0) - v8).to_csv(OUT / "label_count_difference_vs_V8.tsv", sep="\t")

ab = np.asarray(sparse.csr_matrix(B.X).todense(), dtype=float)
ab_total = ab.sum(1)
lg = np.log1p(ab)
clr = lg - lg.mean(1, keepdims=True)
ab_names = pd.Index(B.var_names)
iso = [n for n in ab_names if n.startswith("Isotype")]

mem = pd.read_csv(HOT / "program_membership_v2.tsv", sep="\t")
mem = mem[mem.mapped_symbol_status == "gencode_v49_unique_symbol_confirmed"]
obsv = pd.read_csv(V8 / "outputs/observability_117_programs.tsv", sep="\t")
genes = pd.Index(A.var_names)
cov = np.column_stack([np.ones(A.n_obs), np.log10(rna_umi), np.log10(ab_total + 1)])


def resid_within(y, idx_by_donor):
    r = np.full_like(y, np.nan, dtype=float)
    for ix in idx_by_donor:
        C = cov[ix]
        r[ix] = y[ix] - C @ np.linalg.lstsq(C, y[ix], rcond=None)[0]
    return r


rows, pos = [], []
for ct, label in LINEAGE_LABEL.items():
    m = (A.obs.label == label).to_numpy()
    cells = np.where(m)[0]
    donors = A.obs.gsm.to_numpy()[cells]
    good = [d for d, n in pd.Series(donors).value_counts().items() if n >= MIN_CELLS]
    if len(good) < 5:
        continue
    keep_cells = cells[np.isin(donors, good)]
    by_donor = [keep_cells[A.obs.gsm.to_numpy()[keep_cells] == d] for d in good]
    det = np.asarray((counts[cells] > 0).mean(axis=0)).ravel()
    Xl = A.X[keep_cells]
    progs = obsv[(obsv.cell_type == ct) & (obsv.state == "testable")]
    S = {}
    for _, r in progs.iterrows():
        w = mem[mem.program_uid == r.program_uid].groupby("mapped_symbol").original_l1_weight.sum()
        w = w[[g in genes and det[genes.get_loc(g)] >= 0.01 for g in w.index]]
        M = Xl[:, [genes.get_loc(g) for g in w.index]].toarray()
        sd = M.std(0); sd[sd == 0] = 1
        S[r.module_name] = ((M - M.mean(0)) / sd) @ (w / w.sum()).to_numpy()
    RX = [resid_within(clr[:, j], by_donor) for j in range(len(ab_names))]
    for pname, s in S.items():
        y = np.full(A.n_obs, np.nan); y[keep_cells] = s
        ry = resid_within(y, by_donor)
        for j, an in enumerate(ab_names):
            rx = RX[j]
            rho = np.array([stats.spearmanr(ry[bd], rx[bd])[0] for bd in by_donor])
            z = np.arctanh(np.clip(rho, -0.999, 0.999))
            t = stats.ttest_1samp(z, 0)
            rows.append(dict(cell_type=ct, program=pname, adt=an, isotype=an in iso, n_donors=len(good),
                             median_rho=np.median(rho), n_pos=int((rho > 0).sum()), mean_z=z.mean(), t=t.statistic, p=t.pvalue))
    # RNA-protein positive controls in this lineage
    for gene, an in [("CD8A", "Hu.CD8"), ("CD4", "Hu.CD4_RPA.T4"), ("CD14", "Hu.CD14_M5E2"), ("CD163", "Hu.CD163"),
                     ("PTPRC", "Hu.CD45_HI30"), ("ITGAM", "Hu.CD11b"), ("CD36", "Hu.CD36"), ("NT5E", "Hu.CD73")]:
        if gene not in genes or an not in ab_names:
            continue
        y = np.asarray(A.X[:, genes.get_loc(gene)].todense()).ravel()
        ry, rx = resid_within(y, by_donor), resid_within(clr[:, ab_names.get_loc(an)], by_donor)
        rho = np.array([stats.spearmanr(ry[bd], rx[bd])[0] for bd in by_donor])
        pos.append(dict(cell_type=ct, gene=gene, adt=an, median_rho=np.nanmedian(rho), n_pos=int((rho > 0).sum()), n=len(rho)))
res = pd.DataFrame(rows)
res["bh_q_within_lineage"] = res.groupby("cell_type").p.transform(lambda p: stats.false_discovery_control(p.fillna(1)))
res.sort_values(["cell_type", "p"]).to_csv(OUT / "program_adt_within_donor.tsv", sep="\t", index=False)
pd.DataFrame(pos).to_csv(OUT / "rna_protein_positive_controls.tsv", sep="\t", index=False)
print(res[res.bh_q_within_lineage < 0.05].groupby("cell_type").size().to_string(), flush=True)

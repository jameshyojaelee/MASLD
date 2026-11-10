"""Verifier: compare two-sided exact rules (distance-from-E used by H57 exact_p vs min-likelihood) for T1, D3, PC-B and G3 power."""
import json, numpy as np, pandas as pd
from scipy.stats import hypergeom, nchypergeom_fisher
R = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
deg = pd.read_csv(R + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
gl = pd.read_csv(R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
deg["ens"] = deg.gene.str.split(".").str[0]; deg["DEG"] = (deg.padj < 0.05) & (deg.logFC.abs() > 0.5)
U = deg.merge(gl, left_on="ens", right_on="ensembl")
U = U[U.susie_state_t12.isin(["positive", "tested_le_0.5"])].sort_values("ens").reset_index(drop=True)
U["gen"] = U.susie_state_t12 == "positive"
ce = pd.read_csv(R + "Analysis/Hypotheses/H57-20261015T112421Z/run-20261015T113107Z/control_expression.tsv", sep="\t")
U = U.merge(ce[["gene", "ave_control"]], on="gene")
U = U.sort_values("ens").reset_index(drop=True)
def null(var):
    b = pd.qcut(U[var].rank(method="first"), 10, labels=False)
    pmf = np.array([1.0]); E = 0
    for k in range(10):
        s = U[b == k]; n, g, d = len(s), int(s.gen.sum()), int(s.DEG.sum())
        E += g * d / n; pmf = np.convolve(pmf, hypergeom.pmf(np.arange(min(g, d) + 1), n, g, d))
    return pmf, E
x = None
def pd_(pmf, E, O):
    xx = np.arange(len(pmf)); return float(pmf[np.abs(xx - E) >= abs(O - E) - 1e-9].sum())
def pml(pmf, O):
    return float(pmf[pmf <= pmf[O] * (1 + 1e-7)].sum())
out = {}
for var in ("AveExpr", "ave_control"):
    pmf, E = null(var)
    out[var] = dict(E=round(E, 3), p_distance=round(pd_(pmf, E, 26), 4), p_minlik=round(pml(pmf, 26), 4))
    if var == "AveExpr":
        rej = np.array([pd_(pmf, E, o) <= 0.05 for o in range(len(pmf))])
        pcb = pd.read_csv(R + "Analysis/Hypotheses/H57-20261015T112421Z/run-20261015T113107Z/pcb_replicates.tsv", sep="\t")
        out["PCB_distance_rule"] = dict(enrich=float(rej[pcb.total_enrich14].mean()), deplete=float(rej[pcb.total_deplete12].mean()),
                                         crit=[int(np.where(rej[:21])[0].max()), int(np.where(rej[21:])[0].min() + 21)])
# unstratified G3 under the distance rule
N, G, D = 3769, 426, 220
xs = np.arange(0, D + 1); p0 = hypergeom.pmf(xs, N, G, D); E0 = G * D / N
rej = np.array([p0[np.abs(xs - E0) >= abs(o - E0) - 1e-9].sum() <= 0.05 for o in xs])
out["G3_distance_crit"] = [int(xs[rej & (xs < 25)].max()), int(xs[rej & (xs > 25)].min()), round(float(p0[rej].sum()), 4)]
out["G3_distance_power_up"] = {str(o): round(float(nchypergeom_fisher.pmf(xs, N, G, D, o)[rej].sum()), 3) for o in (1.70, 1.71, 1.75, 1.76, 1.77)}
out["G3_distance_power_dn"] = {str(o): round(float(nchypergeom_fisher.pmf(xs, N, G, D, o)[rej].sum()), 3) for o in (0.45, 0.46, 0.47)}
print(json.dumps(out, indent=1))

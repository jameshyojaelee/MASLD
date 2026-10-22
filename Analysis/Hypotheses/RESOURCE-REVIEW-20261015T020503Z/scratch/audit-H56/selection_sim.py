"""Audit check (no expression data): does selecting calls on the unadjusted estimate b_u bias the
half-retained test? Real 774-unit design, rebuilt A16 coordinates, stored sample weights; simulated genes.
Compares the plan's test t = (b_a - 0.5 b_u)/sqrt(s2(u_a^2 - 0.75 u_u^2)) with the conditional test
t = (D + 0.5 b_u)/sqrt(s2(u_a^2 - u_u^2)), D = b_a - b_u. Also checks the two variance identities."""
import numpy as np, pandas as pd
from scipy.linalg import helmert
from scipy.stats import norm
rng = np.random.default_rng(20260923)
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
P = ROOT + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
dz = pd.read_csv(P + "model_design.tsv", sep="\t").set_index("sample_id")
w = pd.read_csv(P + "lvqw_sample_weights.tsv", sep="\t").set_index("sample_id").sample_weight.loc[dz.index].to_numpy()
new = pd.read_csv(ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz", sep="\t").set_index("sample_id")
ct = [c for c in new.columns if c != "dataset"]; F = new.loc[dz.index, ct].to_numpy(float)
Z = np.log(np.maximum(F, 1e-6)) @ (-helmert(16).T)
Z = (Z - Z.mean(0)) / Z.std(0)
Xu = dz.to_numpy(float)                      # intercept, 4 dataset, sex, group (last)
Xa = np.hstack([Xu[:, :-1], Z, Xu[:, -1:]])  # group last
sw = np.sqrt(w)
def coef_op(X):
    Xw = X * sw[:, None]; XtXi = np.linalg.inv(Xw.T @ Xw)
    return XtXi, (XtXi @ Xw.T) * sw[None, :]   # rows map y -> coefficients (WLS)
XtXu, Hu = coef_op(Xu); XtXa, Ha = coef_op(Xa)
uu2, ua2 = XtXu[-1, -1], XtXa[-1, -1]
hu, ha = Hu[-1], Ha[-1]
Winv = 1 / w
# exact identities under Var(y) = sigma^2 W^-1
cov_au = (ha * hu * Winv).sum(); var_d = ((ha - hu) ** 2 * Winv).sum(); var_k = ((ha - 0.5 * hu) ** 2 * Winv).sum()
print(f"u_u^2 {uu2:.6g}  u_a^2 {ua2:.6g}  SE inflation {np.sqrt(ua2/uu2):.3f}")
print(f"Cov(b_a,b_u)/u_u^2 = {cov_au/uu2:.10f}   Var(D)/(u_a^2-u_u^2) = {var_d/(ua2-uu2):.10f}   Var(b_a-0.5b_u)/(u_a^2-0.75u_u^2) = {var_k/(ua2-0.75*uu2):.10f}")
# simulation: true model y = X0 g + group*beta_g + Z theta + e ; theta chosen so beta_a/beta_u = 0.5 exactly (boundary of T_keep null)
sigma = 0.107 / np.sqrt(uu2)                 # SE(b_u) = 0.107, the median call SE
pi = Ha[:, :] @ Xu[:, -1]                    # not used; omitted-variable map computed below
# regression of Z on Xu (WLS): group coefficient of each Z column
Gz = (Hu @ Z)[-1]                            # how much each Z column shifts with group, given dataset and sex
d = Gz / np.linalg.norm(Gz)
n_g = 120000
rows = []
for beta_u in (0.45, 0.5, 0.55, 0.6, 0.8, 1.0):
    beta_g = 0.5 * beta_u                    # adjusted (true) effect = half of the unadjusted estimand
    theta = d * (beta_u - beta_g) / (Gz @ d) # so that Gz . theta = beta_u - beta_g
    mu = Xu[:, -1] * beta_g + Z @ theta
    hits = dict(n=0, sel=0, plan=0, cond=0)
    for _ in range(8):
        E = mu[None, :] + rng.standard_normal((n_g // 8, len(w))) * (sigma / sw)[None, :]
        bu = E @ hu; ba = E @ ha
        ra = E - E @ Ha.T @ Xa.T; s2 = (ra ** 2 * w).sum(1) / (len(w) - Xa.shape[1])
        se_u = np.sqrt(s2 * uu2)
        sel = (np.abs(bu / se_u) > 3.5) & (bu > 0.5)   # stand-in for padj<0.05 and logFC>0.5 (positive calls)
        t_plan = (ba - 0.5 * bu) / np.sqrt(s2 * (ua2 - 0.75 * uu2))
        t_cond = ((ba - bu) + 0.5 * bu) / np.sqrt(s2 * (ua2 - uu2))
        z = norm.ppf(0.95)
        hits["n"] += len(bu); hits["sel"] += sel.sum()
        hits["plan"] += (sel & (t_plan > z)).sum(); hits["cond"] += (sel & (t_cond > z)).sum()
    check_bu = bu.mean()
    rows.append([beta_u, hits["sel"] / hits["n"], hits["plan"] / max(hits["sel"], 1), hits["cond"] / max(hits["sel"], 1)])
print("true beta_u, true ratio 0.5 (null boundary of T_keep); one-sided alpha 0.05 in one cell (A16 rebuilt)")
print(pd.DataFrame(rows, columns=["beta_u", "P(selected)", "plan_test_reject_rate_given_selected", "conditional_test_reject_rate_given_selected"]).round(4).to_string(index=False))

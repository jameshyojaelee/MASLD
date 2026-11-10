"""Independent re-derivation of Figure 1B numbers from the adopted 774-unit bulk fit and COLOC r2."""
import json
import numpy as np, pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.contingency_tables import StratifiedTable
A = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
deg = pd.read_csv(A + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
col = pd.read_csv(A + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
out = {}
deg["ens"] = deg["gene"].str.replace(r"\.\d+$", "", regex=True)
deg["is_deg"] = (deg["padj"] < 0.05) & (deg["logFC"].abs() > 0.5)
out["n_tested"] = len(deg); out["n_deg_all"] = int(deg.is_deg.sum())
out["deg_ens_dup"] = int(deg.ens.duplicated().sum())
m = deg.merge(col, left_on="ens", right_on="ensembl", how="left")
out["deg_by_state"] = m[m.is_deg].susie_state_t12.fillna("NO_RECORD").value_counts().to_dict()
out["overlap_universe_any_record"] = int(m.ensembl.notna().sum())
u = m[m.susie_state_t12.isin(["positive", "tested_le_0.5"])].copy()
u["is_gen"] = u.susie_state_t12.eq("positive")
N, G, D = len(u), int(u.is_gen.sum()), int(u.is_deg.sum()); O = int((u.is_gen & u.is_deg).sum())
out.update(N=N, G=G, D=D, O=O, E_unstrat=G * D / N)
u["bin"] = pd.qcut(u.AveExpr, 10, labels=False)
tabs, pmf = [], np.array([1.0]); E = 0.0
for k, s in u.groupby("bin"):
    n, g, d, o = len(s), int(s.is_gen.sum()), int(s.is_deg.sum()), int((s.is_gen & s.is_deg).sum())
    E += g * d / n
    x = np.arange(0, min(g, d) + 1); pk = hypergeom.pmf(x, n, g, d)
    pmf = np.convolve(pmf, pk)
    tabs.append(np.array([[o, g - o], [d - o, n - g - d + o]]))
out["E_strat10"] = E
cdf = np.cumsum(pmf)
out["null_q025"] = int(np.searchsorted(cdf, 0.025)); out["null_q975"] = int(np.searchsorted(cdf, 0.975))
out["null_median"] = int(np.searchsorted(cdf, 0.5))
out["P_ge_obs"] = float(pmf[O:].sum()); out["P_le_obs"] = float(pmf[:O + 1].sum())
out["P_two_minlike"] = float(pmf[pmf <= pmf[O] * (1 + 1e-7)].sum())
out["P_two_doubled"] = float(min(1, 2 * min(out["P_ge_obs"], out["P_le_obs"])))
st = StratifiedTable(tabs)
lo, hi = st.oddsratio_pooled_confint()
out.update(MH_OR=st.oddsratio_pooled, MH_lo=lo, MH_hi=hi, MH_chi2_p=float(st.test_null_odds(correction=True).pvalue))
# outside the comparison
ddeg = m[m.is_deg]
out["deg_no_multisignal_result"] = int((~ddeg.susie_state_t12.isin(["positive", "tested_le_0.5"])).sum())
out["deg_no_record"] = int(ddeg.ensembl.isna().sum())
out["deg_untestable"] = int(ddeg.susie_state_t12.eq("untestable").sum())
gm = pd.read_csv(A + "data/gencode_v49_gene_metadata.tsv.gz", sep="\t")
print(gm.columns.tolist()[:12])
idc = [c for c in gm.columns if "gene_id" in c.lower()][0]; btc = [c for c in gm.columns if "type" in c.lower() or "biotype" in c.lower()][0]
gm["ens"] = gm[idc].str.replace(r"\.\d+$", "", regex=True)
bt = dict(zip(gm.ens, gm[btc]))
nor = ddeg[ddeg.ensembl.isna()]
out["deg_no_record_lncRNA"] = int(nor.ens.map(bt).eq("lncRNA").sum())
out["deg_no_record_biotypes"] = nor.ens.map(bt).value_counts().head(6).to_dict()
out["bins"] = [dict(bin=int(k), n=len(s), G=int(s.is_gen.sum()), D=int(s.is_deg.sum()), O=int((s.is_gen & s.is_deg).sum())) for k, s in u.groupby("bin")]
print(json.dumps(out, indent=1, default=float))
json.dump(out, open("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig1/rederive_1b.json", "w"), indent=1, default=float)

"""Skeptic re-derivation of Figure 1C, plus H57 D3/D10/PC-B and the 935/699 sub-claims (verifier, 2026-10-26)."""
import json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

R = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
C = R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/"
tier = pd.read_csv(R + "GWAS/finemapping/config/gwas_trait_tier.tsv", sep="\t")
tier = tier[tier.placement == "main"]
cls = dict(zip(tier.study_name, tier.tier.astype(int).map({1: "direct", 2: "enzyme"})))
res = {"n_studies": tier.tier.astype(int).value_counts().to_dict()}

a = pd.read_csv(C + "susie_coloc_all_gwas.csv",
                usecols=["gwas_name", "gene", "ensembl", "PP.H4.susie", "PP.H4.abf", "top_snp", "method", "eligibility_source"],
                dtype={"top_snp": str})
res["rows_all"] = len(a)
a = a[a.gwas_name.isin(cls)].copy()
res["rows_t12"] = len(a)
res["studies_present"] = int(a.gwas_name.nunique())
a["cls"] = a.gwas_name.map(cls)
a["ens"] = a.ensembl.astype(str).str.split(".").str[0]
t = a[np.isfinite(a["PP.H4.susie"])]
tested = {c: set(t.ens[t.cls == c]) for c in ("enzyme", "direct")}
supp = {c: set(t.ens[(t.cls == c) & (t["PP.H4.susie"] > 0.5)]) for c in ("enzyme", "direct")}
res["tested_enzyme"] = len(tested["enzyme"]); res["tested_direct"] = len(tested["direct"])
res["tested_any"] = len(tested["enzyme"] | tested["direct"])
res["supp_enzyme"] = len(supp["enzyme"]); res["supp_direct"] = len(supp["direct"])
both = supp["enzyme"] & supp["direct"]
res["supp_both"] = len(both); res["supp_union"] = len(supp["enzyme"] | supp["direct"])
for c, o in (("enzyme", "direct"), ("direct", "enzyme")):
    s = supp[c]
    res[f"{c}_bar"] = dict(also=len(s & supp[o]), tested_le=len((s & tested[o]) - supp[o]), not_tested=len(s - tested[o]))
tb = tested["enzyme"] & tested["direct"]
res["tested_both"] = len(tb)
res["tested_both_enz_supp"] = len(tb & supp["enzyme"]); res["tested_both_dir_supp"] = len(tb & supp["direct"])
res["tested_both_both"] = len(tb & both)
# proportions supported among tested, per class
res["frac_enzyme"] = round(res["supp_enzyme"] / res["tested_enzyme"], 4)
res["frac_direct"] = round(res["supp_direct"] / res["tested_direct"], 4)
res["frac_in_tb_enzyme"] = round(res["tested_both_enz_supp"] / len(tb), 4)
res["frac_in_tb_direct"] = round(res["tested_both_dir_supp"] / len(tb), 4)

# gene-level table agreement
gl = pd.read_csv(C + "gene_level_coloc_tier12.csv")
pos = set(gl.ensembl[gl.susie_state_t12 == "positive"])
res["gl_pos_equals_union"] = pos == (supp["enzyme"] | supp["direct"])
unt = gl[gl.susie_state_t12 == "untestable"]
res["untestable_genes"] = len(unt)
res["untestable_have_finite_le05"] = int(sum(e in (tested["enzyme"] | tested["direct"]) for e in unt.ensembl))

# DEG sub-claims (F1-05, 699 record check)
deg = pd.read_csv(R + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
deg["ens"] = deg.gene.str.split(".").str[0]
deg["DEG"] = (deg.padj < 0.05) & (deg.logFC.abs() > 0.5)
m = deg.merge(gl, left_on="ens", right_on="ensembl", how="left")
d935 = m[m.DEG & (m.susie_state_t12 == "no_susie_result")]
res["n935"] = len(d935)
res["n935_abf_finite"] = int(d935.coloc_best_abf_pp4.notna().sum())
res["n935_abf_gt05"] = int((d935.coloc_best_abf_pp4 > 0.5).sum())
# methods for the 935 in tier-1/2 rows
mm = a[a.ens.isin(set(d935.ens))].method.value_counts().to_dict()
res["n935_rows_by_method"] = mm
d699 = m[m.DEG & m.susie_state_t12.isna()]
res["n699"] = len(d699)
# any record in the all-tier gene table?
gall = pd.read_csv(C + "gene_level_coloc.csv", usecols=lambda c: c in ("ensembl", "ensembl_id", "gene"))
key = "ensembl" if "ensembl" in gall.columns else gall.columns[0]
gall_ens = set(gall[key].astype(str).str.split(".").str[0])
res["n699_with_any_tier_record"] = int(d699.ens.isin(gall_ens).sum())
# any row at all in susie_coloc_all_gwas (all studies)
res["gene_level_all_cols"] = list(gall.columns)

# D10: one gene per 500 kb cluster
U = m[m.susie_state_t12.isin(["positive", "tested_le_0.5"])].copy()
U["gen"] = U.susie_state_t12 == "positive"
gpos = U[U.gen][["ens", "coloc_best_susie_pp4", "susie_driving_gwas"]]
top = a[["ens", "gwas_name", "top_snp"]].drop_duplicates(["ens", "gwas_name"])
gpos = gpos.merge(top, left_on=["ens", "susie_driving_gwas"], right_on=["ens", "gwas_name"], how="left")
res["D10_pos_with_topsnp"] = int(gpos.top_snp.notna().sum())
sp = gpos.top_snp.str.split(":", expand=True)
gpos["c"] = sp[0]; gpos["p"] = pd.to_numeric(sp[1])
lab = {}
for ch, g in gpos.sort_values(["c", "p", "ens"]).groupby("c"):
    k = 0; prev = None
    for e, p in zip(g.ens, g.p):
        if prev is not None and p - prev > 500_000:
            k += 1
        lab[e] = f"{ch}_{k}"; prev = p
gpos["cl"] = gpos.ens.map(lab)
res["D10_clusters"] = int(gpos.cl.nunique())
keep = set(gpos.sort_values(["cl", "coloc_best_susie_pp4", "ens"], ascending=[True, False, True]).drop_duplicates("cl").ens)
U10 = U[~(U.gen & ~U.ens.isin(keep))].copy()


def strat(df, binvar, nb=10):
    df = df.sort_values("ens").reset_index(drop=True)
    df["b"] = pd.qcut(df[binvar].rank(method="first"), nb, labels=False)
    E = 0; pmf = np.array([1.0]); sR = sS = sP = sPQ = sQ = 0.0
    for _, s in df.groupby("b"):
        n = len(s); g = int(s.gen.sum()); d = int(s.DEG.sum()); aa = int((s.gen & s.DEG).sum())
        b = g - aa; c = d - aa; dd = n - g - d + aa
        E += g * d / n
        pmf = np.convolve(pmf, hypergeom.pmf(np.arange(0, min(g, d) + 1), n, g, d))
        Rk, Sk, Pk, Qk = aa * dd / n, b * c / n, (aa + dd) / n, (b + c) / n
        sR += Rk; sS += Sk; sP += Pk * Rk; sPQ += Pk * Sk + Qk * Rk; sQ += Qk * Sk
    O = int((df.gen & df.DEG).sum())
    orr = sR / sS
    se = np.sqrt(sP / (2 * sR ** 2) + sPQ / (2 * sR * sS) + sQ / (2 * sS ** 2))
    p = float(pmf[pmf <= pmf[O] * (1 + 1e-7)].sum())
    return dict(N=len(df), G=int(df.gen.sum()), D=int(df.DEG.sum()), O=O, E=round(E, 3), OR=round(orr, 3),
                lo=round(orr * np.exp(-1.96 * se), 3), hi=round(orr * np.exp(1.96 * se), 3), p=round(p, 4)), pmf


res["D10_test"], _ = strat(U10, "AveExpr")

# D3: control-only stratifier
ce = pd.read_csv(R + "Analysis/Hypotheses/H57-20261015T112421Z/run-20261015T113107Z/control_expression.tsv", sep="\t")
U3 = U.merge(ce[["gene", "ave_control"]], on="gene", how="left")
res["D3_missing"] = int(U3.ave_control.isna().sum())
res["D3_test"], _ = strat(U3, "ave_control")

# PC-B detection: bins and margins are unchanged by within-bin label swaps, so the p of a total depends only on the total
_, pmf10 = strat(U, "AveExpr")
pminl = np.array([pmf10[pmf10 <= v * (1 + 1e-7)].sum() for v in pmf10])
pcb = pd.read_csv(R + "Analysis/Hypotheses/H57-20261015T112421Z/run-20261015T113107Z/pcb_replicates.tsv", sep="\t")
res["PCB_reps"] = len(pcb)
res["PCB_detect_enrich"] = round(float((pminl[pcb.total_enrich14.to_numpy()] <= 0.05).mean()), 3)
res["PCB_detect_deplete"] = round(float((pminl[pcb.total_deplete12.to_numpy()] <= 0.05).mean()), 3)
res["pmf10_rejection_counts"] = [int(x) for x in np.where(pminl <= 0.05)[0] if x < 22][-1:] + [int(x) for x in np.where(pminl <= 0.05)[0] if x > 22][:1]

# NC-1 control count: recompute strat O/E expectation 13 / 0.6105
nc = pd.read_csv(R + "Analysis/Hypotheses/H57-20261015T112421Z/run-20261015T113107Z/nc1_lists.tsv", sep="\t")
ue = nc[nc.universe == "U_eval"]
res["NC1_Ueval_lists"] = len(ue); res["NC1_O_unique"] = sorted(ue.O.unique().tolist())
res["NC1_p_unique"] = sorted(ue.p_exact.round(4).unique().tolist())
res["NC1_E_implied"] = round(float((ue.O / ue.strat_OE).iloc[0]), 2)

json.dump(res, open(__file__.replace(".py", ".json"), "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))

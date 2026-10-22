# Re-derive the Myojin HLF palmitate class test (Plan 40 spec) with v1 classes
# (to reproduce the stored estimates) and with r8 Catalog classes.
# Unit: gene (CRISPR knockout screen); outcome PA_vs_vehicle positive-selection LFC.
# Model (spec): lfc ~ class + log1p_HLF_TPM + guide_count_min + HLF_Chronos, OLS, HC3.
import numpy as np, pandas as pd, statsmodels.api as sm, sys
MY = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf"
R8 = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/passport_gene_index.parquet"
OUT = sys.argv[1]
sr = pd.read_csv(f"{MY}/screen_results.tsv", sep="\t", low_memory=False)
gu = pd.read_csv(f"{MY}/screen_gene_universe.tsv", sep="\t", low_memory=False, usecols=["source_row","gencode_ensembl_base"])
sr = sr.merge(gu, on="source_row", how="left")
prim = sr[sr["primary_eligible"].astype(str).str.upper()=="TRUE"].copy()
print("primary eligible genes:", len(prim))
print("v1 classes:", prim["primary_evidence_class"].value_counts().to_dict())
y = "PA_vs_vehicle_positive_selection_LFC"
cov = ["log1p_HLF_TPM","guide_count_min","HLF_Chronos"]
def fit(df, cls, ref, levels, label):
    X = pd.DataFrame({f"c_{l}": (df[cls]==l).astype(float) for l in levels})
    X = pd.concat([X, df[cov].astype(float)], axis=1); X = sm.add_constant(X)
    m = sm.OLS(df[y].astype(float), X).fit(cov_type="HC3")
    m0 = sm.OLS(df[y].astype(float), X).fit()
    rsd = np.sqrt(m0.ssr/(len(df)-X.shape[1]))
    rows=[]
    for l in levels:
        b=m.params[f"c_{l}"]; se=m.bse[f"c_{l}"]
        rows.append(dict(analysis=label, contrast=f"{l}_vs_{ref}", n_total=len(df), n_class=int((df[cls]==l).sum()),
                         n_ref=int((df[cls]==ref).sum()), estimate=b, HC3_SE=se, ci_low=b-1.959964*se, ci_high=b+1.959964*se,
                         std_effect=b/rsd, wald_p_HC3=m.pvalues[f"c_{l}"], residual_SD=rsd))
    # omnibus HC3 Wald for the class indicators
    R = np.zeros((len(levels), X.shape[1]))
    for i,l in enumerate(levels): R[i, list(X.columns).index(f"c_{l}")]=1
    w = m.wald_test(R, scalar=True)
    rows.append(dict(analysis=label, contrast="omnibus_HC3_wald", n_total=len(df), estimate=float(w.statistic), wald_p_HC3=float(w.pvalue)))
    return rows
out=[]
v1 = prim[prim["primary_evidence_class"].isin(["neither","genetic_only","disease_state_only","convergent"])]
out += fit(v1, "primary_evidence_class", "neither", ["genetic_only","disease_state_only","convergent"], "v1_classes_reproduction")
r8 = pd.read_parquet(R8, columns=["ensembl_id","primary_evidence_class"]).rename(columns={"primary_evidence_class":"r8_class"})
r8["ens"] = r8["ensembl_id"].str.replace(r"\..*$","",regex=True)
prim["ens"] = prim["gencode_ensembl_base"].astype(str).str.replace(r"\..*$","",regex=True)
p8 = prim.merge(r8[["ens","r8_class"]], on="ens", how="left")
p8["r8_class"] = p8["r8_class"].fillna("not_in_r8")
print("r8 classes among primary-eligible screen genes:", p8["r8_class"].value_counts().to_dict())
print("crosstab v1 x r8:\n", pd.crosstab(p8["primary_evidence_class"], p8["r8_class"]))
lv = ["genetically_anchored","established_state_associated","concordant"]
a = p8[p8["r8_class"].isin(["unresolved"]+lv)]
out += fit(a, "r8_class", "unresolved", lv, "r8_classes_ref_unresolved_excl_untested")
b = p8[p8["r8_class"].isin(["unresolved","untested"]+lv)].copy(); b.loc[b["r8_class"]=="untested","r8_class"]="unresolved"
out += fit(b, "r8_class", "unresolved", lv, "r8_classes_ref_unresolved_plus_untested")
res = pd.DataFrame(out); res.to_csv(OUT, sep="\t", index=False)
pd.set_option("display.width",250); print(res.to_string())

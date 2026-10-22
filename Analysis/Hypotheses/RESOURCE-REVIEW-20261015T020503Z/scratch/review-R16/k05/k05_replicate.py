# K05: replicate the Figure 2B x-axis and recompute it with SuSiE-only PIP.
import pandas as pd, numpy as np, sys
S = sys.argv[1]
src = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc_source.csv"
cols = ["source","study","chromosome","position","susie_pip","recommended_pip","carma_pip","susiex_pip","mesusie_pip"]
fm = pd.concat([pd.read_csv(f"{S}/k05/fm_{t}.tsv", sep="\t", header=None, names=cols, dtype=str)
                for t in ["canonical_finemapping","coloc_targeted_finemapping"]])
for c in cols[4:]:
    fm[c] = pd.to_numeric(fm[c], errors="coerce")
fm["chromosome"] = fm["chromosome"].astype(int); fm["position"] = fm["position"].astype(int)
fm["causal_pip"] = fm[["recommended_pip","susie_pip"]].max(axis=1, skipna=True)
fm = fm[np.isfinite(fm["causal_pip"])]
fm["rank"] = np.where(fm["source"]=="canonical_finemapping",1,2)
fm = fm.sort_values(["study","chromosome","position","causal_pip","rank"], ascending=[True,True,True,False,True])
first = fm.drop_duplicates(["study","chromosome","position"])          # producer's fm_unique row
susie_any = fm.groupby(["study","chromosome","position"])["susie_pip"].max().rename("susie_key_max")
first = first.join(susie_any, on=["study","chromosome","position"])
co = pd.read_csv(f"{S}/k05/coloc_main_rows.tsv", sep="\t", header=None,
                 names=["gwas","gene","chromosome","position","pp4_susie","pp4_abf"], dtype={"gene":str})
co["chromosome"] = pd.to_numeric(co["chromosome"], errors="coerce"); co["position"] = pd.to_numeric(co["position"], errors="coerce")
co = co.dropna(subset=["chromosome","position"]); co[["chromosome","position"]] = co[["chromosome","position"]].astype(int)
co["pp4_susie"] = pd.to_numeric(co["pp4_susie"], errors="coerce")
m = co.merge(first[["study","chromosome","position","causal_pip","susie_pip","susie_key_max","source"]],
             left_on=["gwas","chromosome","position"], right_on=["study","chromosome","position"], how="inner")
m = m.sort_values(["gene","causal_pip","source"], ascending=[True,False,True])
gx = m.drop_duplicates("gene")[["gene","gwas","chromosome","position","causal_pip","susie_pip","source"]].rename(columns={"gwas":"x_gwas"})
gs = m.groupby("gene")["susie_key_max"].max().rename("x_susie_only")          # SuSiE-only, best over the same lead variants
gy = co.dropna(subset=["pp4_susie"]).sort_values("pp4_susie", ascending=False).drop_duplicates("gene")[["gene","gwas","pp4_susie"]].rename(columns={"gwas":"y_gwas"})
s = pd.read_csv(src); p = s[s["plotted"]==True].copy()
r = p.merge(gx, on="gene", how="left").merge(gs, on="gene", how="left").merge(gy, on="gene", how="left")
r["repro"] = (r["causal_pip_x"] - r["causal_pip_y"].round(4)).abs() <= 6e-5 if "causal_pip_y" in r else False
r = r.rename(columns={"causal_pip_x":"x_plotted","causal_pip_y":"x_replicated"})
r["repro"] = (r["x_plotted"] - r["x_replicated"].round(4)).abs() <= 6e-5
d = r["x_plotted"] - r["x_susie_only"].fillna(0)
out = []
out.append(("plotted_points", len(r)))
out.append(("x_reproduced_within_rounding", int(r["repro"].sum())))
out.append(("y_gwas_equals_source_pp4_gwas", int((r["y_gwas"]==r["pp4_gwas"]).sum())))
out.append(("x_gwas_differs_from_y_gwas", int((r["x_gwas"]!=r["pp4_gwas"]).sum())))
out.append(("x_susie_only_missing", int(r["x_susie_only"].isna().sum())))
for t in [0.01,0.1,0.25]:
    out.append((f"x_drops_by_more_than_{t}_under_susie_only", int((d>t).sum())))
out.append(("x_ge_0.5_plotted", int((r["x_plotted"]>=0.5).sum())))
out.append(("x_ge_0.5_susie_only", int((r["x_susie_only"]>=0.5).sum())))
out.append(("x_ge_0.99_plotted", int((r["x_plotted"]>=0.99).sum())))
out.append(("x_ge_0.99_susie_only", int((r["x_susie_only"]>=0.99).sum())))
out.append(("x_exactly_1_plotted", int((r["x_plotted"]>=0.99995).sum())))
for k,v in out: print(f"{k}\t{v}")
r[["gene","pp4_gwas","pp4","x_gwas","x_top_snp","x_plotted","x_replicated","x_susie_only","repro","is_label"]].to_csv(f"{S}/k05/fig2b_x_susie_only.tsv", sep="\t", index=False)
print(r[r["is_label"]==True][["gene","pp4_gwas","pp4","x_gwas","x_top_snp","x_plotted","x_replicated","x_susie_only"]].to_string())

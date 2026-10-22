"""H58 plan: count instances (credible sets) by trait class x ancestry x study x LD panel.
Reads metadata columns only; computes no coding fraction."""
import pandas as pd
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
D=ROOT+"figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/"
a=pd.read_csv(D+"multiplicity/corrected_locus_coordinates.tsv",sep="\t",dtype={"locus":str},
              usecols=["credible_set_uid","study","trait","trait_scope","ancestry","pip_sum_gate_passed","chrom","pos"])
print("rows",len(a),"gate passed",a.pip_sum_gate_passed.astype(str).str.lower().eq("true").sum())
reg=pd.read_csv(D+"snapshot/GWAS/finemapping/config/gwas_registry.tsv",sep="\t")
print("registry cols:",list(reg.columns))
a=a.merge(reg[["study_name","ld_panel"]].rename(columns={"study_name":"study"}),on="study",how="left")
print("missing ld_panel",a.ld_panel.isna().sum())
pd.set_option("display.width",200)
print(pd.crosstab(a.trait_scope,a.ancestry,margins=True))
print(pd.crosstab([a.ancestry],a.ld_panel))
print(a.groupby(["trait_scope","ancestry","study","ld_panel"]).size().to_string())

"""audit-H58: fine-mapped loci vs registry lead rows for EUR and AMR studies in the coding set. Counts only."""
import pandas as pd, os
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
reg=pd.read_csv(ROOT+"figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/snapshot/GWAS/finemapping/config/gwas_registry.tsv",sep="\t")
cs=pd.read_csv(ROOT+"GWAS/finemapping/results/credible_sets.csv",usecols=["study","locus"],dtype=str)
fm=cs.groupby("study").locus.nunique()
a=pd.read_csv(ROOT+"figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv",sep="\t",usecols=["study","trait_scope","ancestry"])
st=a[a.ancestry.isin(["EUR","AMR"])].drop_duplicates("study")
for _,r in st.sort_values(["ancestry","trait_scope","study"]).iterrows():
    rr=reg[reg.study_name==r.study].iloc[0]
    lp=ROOT+"GWAS/finemapping/"+str(rr.leadsnps_path)
    nl=None
    if os.path.exists(lp):
        l=pd.read_csv(lp,sep="\t"); nl=len(l)
    print(f"{r.ancestry} {r.trait_scope[:6]} {r.study:<32} lead_rows={nl} finemapped_loci={fm.get(r.study,0)} N={rr.N_tot} lead_path={rr.leadsnps_path}")

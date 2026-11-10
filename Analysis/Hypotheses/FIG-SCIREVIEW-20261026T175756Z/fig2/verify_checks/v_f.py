import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
x=pd.read_pickle("t12.pkl")
f=pd.read_csv(R+"figures/candidates/panel-deep-fixes-20261017T035633Z/fig2/source_tables/2F.tsv",sep="\t")
genes=f.gene.unique()
s=x[x.gene.isin(genes)]
anc=["EUR","AFR","AMR","EAS","SAS"]
rows=[]
for g in genes:
    for a in anc:
        c=s[(s.gene==g)&(s.reg_anc==a)]
        sus=c["PP.H4.susie"]; abf=c["PP.H4.abf"]
        n_fin=sus.notna().sum()
        mx_s=sus.max() if n_fin else np.nan
        mx_a=abf.max()
        unt=(c.method=="susie_untestable_insufficient_shared_posterior").sum()
        state="dark" if (n_fin and mx_s>0.5) else ("light" if mx_a>0.5 else "gray")
        rows.append(dict(gene=g,anc=a,n_rows=len(c),n_susie_finite=n_fin,max_susie=mx_s,max_abf=mx_a,n_untestable=unt,state=state))
d=pd.DataFrame(rows)
print(pd.crosstab(d.state, d.n_susie_finite>0))
# per gene: SuSiE groups supported / tested
gs=d.groupby("gene").apply(lambda z: pd.Series(dict(sup=(z.state=="dark").sum(), tested=(z.n_susie_finite>0).sum(), ge05=((z.max_susie>0.5)|(z.max_abf>0.5)).sum())))
print(gs.sort_values("sup",ascending=False).to_string())
# compare with 2F.tsv evidence_state
m=d.merge(f[["gene","ancestry","evidence_state","display_method","multi_pp4","single_pp4"]],left_on=["gene","anc"],right_on=["gene","ancestry"])
print(pd.crosstab(m.state,m.evidence_state))
print(m[m.gene=="CHEK2"].to_string())
print(m[m.gene=="EPHA2"].to_string())
# gray cells: any SuSiE tested?
print("gray cells w/ susie finite:", ((d.state=="gray")&(d.n_susie_finite>0)).sum(), " light w/ susie finite:",((d.state=="light")&(d.n_susie_finite>0)).sum())

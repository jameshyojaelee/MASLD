import pandas as pd, numpy as np
RUN='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/'
m=pd.read_pickle(RUN+'scratch/review-R02/r2_main.pkl')
p=m[m['PP.H4.susie']>0.5].copy()
print('SuSiE-positive rows',len(p),'genes',p.ensembl.nunique())
print(p[['gwas_name','top_snp','ld_panel']].head(3).to_string())
print(p.groupby('gwas_name').ensembl.nunique().sort_values(ascending=False).to_string())
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/'
c=pd.read_csv(D+'corrected_locus_coordinates.tsv',sep='\t')
c['fit']=c.study+'|'+c.chrom.astype(str)+':'+c.pos.astype(str)
c['k']=c.groupby('fit').credible_set_uid.transform('size')
cap=c[c.k==10]
# top_snp position
pos=p.top_snp.astype(str).str.extract(r'^(?:chr)?([0-9XY]+)[:_](\d+)')
p['tchr']=pos[0]; p['tpos']=pd.to_numeric(pos[1],errors='coerce')
print('top_snp parsed',p.tpos.notna().sum(),'of',len(p))
hit=[]
for i,r in p.iterrows():
    s=cap[(cap.study==r.gwas_name)&(cap.chrom.astype(str)==str(r.tchr))]
    hit.append(bool(((s.member_min_position-500000<=r.tpos)&(s.member_max_position+500000>=r.tpos)).any()) if pd.notna(r.tpos) else False)
p['in_capped_fit_window']=hit
print('positive rows near a capped GWAS fit (+-500kb):',p.in_capped_fit_window.sum())
g=p.groupby('ensembl').in_capped_fit_window.all()
print('genes whose every positive row is near a capped fit:',g.sum(),'of',len(g))

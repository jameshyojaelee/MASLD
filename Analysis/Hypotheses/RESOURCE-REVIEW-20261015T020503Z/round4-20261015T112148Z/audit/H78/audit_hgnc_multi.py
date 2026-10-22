"""H78 audit 3: Hart HGNC IDs with more than one primary-chromosome ENSG in GENCODE v49 (keys only)."""
import re, gzip, pandas as pd
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
rows=[]
with gzip.open('/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz','rt') as f:
    for line in f:
        if line.startswith('#'): continue
        p=line.split('\t',9)
        if p[2]!='gene': continue
        a=p[8]; hg=re.search(r'hgnc_id "([^"]+)"',a)
        rows.append((re.search(r'gene_id "([^"]+)"',a).group(1).split('.')[0], hg.group(1) if hg else None, p[0]))
gc=pd.DataFrame(rows,columns=['ens','hgnc','chr']).drop_duplicates('ens')
gp=gc[gc.chr.str.match(r'^chr([0-9]+|X|Y|M)$')].dropna(subset=['hgnc'])
mult=gp.groupby('hgnc').ens.nunique(); mult=set(mult[mult>1].index)
gn=pd.read_csv(R+'data/external/gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv',sep='\t',usecols=['gene_id','mane_select','lof.oe_ci.upper'],low_memory=False)
mane=set(gn[gn.mane_select.astype(str).str.lower().eq('true') & gn['lof.oe_ci.upper'].notna()].gene_id)
for nm in ['CEGv2','NEGv1']:
    h=pd.read_csv(R+f'data/public_r4/Hart2017_CEG2/{nm}.txt',sep='\t')
    x=h[h.HGNC_ID.isin(mult)]
    print(nm,'HGNC IDs with >1 primary ENSG',len(x), x.GENE.tolist()[:20])
    for hg in x.HGNC_ID:
        es=gp[gp.hgnc==hg].ens.tolist(); print('  ',hg,es,'MANE-LOEUF:',[e in mane for e in es])

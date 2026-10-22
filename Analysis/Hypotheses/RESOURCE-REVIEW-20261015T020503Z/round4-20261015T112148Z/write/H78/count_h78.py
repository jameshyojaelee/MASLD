"""H78 prespec counts: coverage and key mapping only. No route-level outcome is computed."""
import re, gzip, pandas as pd, numpy as np
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
g=pd.read_parquet(R+'figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet')
g['ens']=g.ensembl_id.astype(str).str.split('.').str[0]
print('catalog rows',len(g),'unique ens',g.ens.nunique(),'non-ENSG ids',(~g.ens.str.startswith('ENSG')).sum())
print('route x class'); print(pd.crosstab(g.next_experiment_rule_id,g.primary_evidence_class,margins=True).to_string())
gn=pd.read_csv(R+'data/external/gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv',sep='\t',
   usecols=['gene','gene_id','transcript','mane_select','canonical','lof.oe_ci.upper','lof.pLI','lof.oe_ci.upper_bin_decile','constraint_flags','cds_length'],low_memory=False)
gn['mane']=gn.mane_select.astype(str).str.lower().eq('true')
print('gnomAD rows',len(gn),'mane rows',gn.mane.sum())
m=gn[gn.mane]
ensm=m[m.gene_id.astype(str).str.startswith('ENSG')].copy(); ezm=m[m.gene_id.astype(str).str.fullmatch(r'\d+')].copy()
print('MANE ENSG rows',len(ensm),'unique ENSG',ensm.gene_id.nunique(),'| MANE Entrez rows',len(ezm),'unique',ezm.gene_id.nunique(),'| other',len(m)-len(ensm)-len(ezm))
print('MANE ENSG with LOEUF',ensm['lof.oe_ci.upper'].notna().sum(),'with pLI',ensm['lof.pLI'].notna().sum(),'with official decile',ensm['lof.oe_ci.upper_bin_decile'].notna().sum())
print('MANE ENSG constraint_flags non-empty',ensm.constraint_flags.astype(str).str.strip('[]').str.len().gt(0).sum())
L=ensm.dropna(subset=['lof.oe_ci.upper']).drop_duplicates('gene_id')
# ENSG->Entrez through MANE symbol pairs
s2e=ezm.drop_duplicates('gene').set_index('gene').gene_id.astype(int)
L['entrez']=L.gene.map(s2e)
print('ENSG MANE rows with an Entrez partner by symbol',L.entrez.notna().sum(),'symbols with >1 Entrez MANE row',int((ezm.gene.value_counts()>1).sum()))
dmcols=pd.read_csv(R+'Analysis/downstream_analysis/essentiality/CRISPRGeneEffect.csv',nrows=0).columns[1:]
dm=pd.DataFrame({'col':dmcols}); dm['sym']=dm.col.str.split(' \(').str[0]; dm['entrez']=dm.col.str.extract(r'\((\d+)\)')[0].astype(int)
print('DepMap genes',len(dm),'unique entrez',dm.entrez.nunique(),'unique symbol',dm.sym.nunique())
# GENCODE v49 ENSG->symbol, HGNC
rows=[]
with gzip.open('/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz','rt') as f:
    for line in f:
        if line.startswith('#'): continue
        p=line.split('\t',9)
        if p[2]!='gene': continue
        a=p[8]
        gid=re.search(r'gene_id "([^"]+)"',a).group(1).split('.')[0]
        nm=re.search(r'gene_name "([^"]+)"',a); hg=re.search(r'hgnc_id "([^"]+)"',a)
        rows.append((gid,nm.group(1) if nm else None,hg.group(1) if hg else None,p[0]))
gc=pd.DataFrame(rows,columns=['ens','sym','hgnc','chr'])
gcp=gc[~gc.ens.str.endswith('_PAR_Y')].drop_duplicates('ens')
print('GENCODE v49 gene rows',len(gc),'unique ENSG',gc.ens.nunique(),'with hgnc_id',gc.hgnc.notna().sum())
g=g.merge(gcp[['ens','sym','hgnc']],on='ens',how='left')
print('Catalog genes found in GENCODE v49',g.sym.notna().sum(),'with hgnc_id',g.hgnc.notna().sum())
# DepMap via Entrez (from gnomAD symbol pair) and via GENCODE symbol
g=g.merge(L[['gene_id','gene','entrez','lof.oe_ci.upper']].rename(columns={'gene_id':'ens','gene':'gn_sym'}),on='ens',how='left')
g['dm_by_entrez']=g.entrez.isin(dm.entrez); g['dm_by_sym']=g.sym.isin(dm.sym)
e2c=dm.set_index('entrez').col; s2c=dm.drop_duplicates('sym').set_index('sym').col
g['c_e']=g.entrez.map(e2c); g['c_s']=g.sym.map(s2c)
both=g.c_e.notna()&g.c_s.notna()
print('DepMap mapped: by entrez',g.c_e.notna().sum(),'by symbol',g.c_s.notna().sum(),'either',(g.c_e.notna()|g.c_s.notna()).sum(),'both and disagree',int((both&(g.c_e!=g.c_s)).sum()))
hpa=pd.read_csv(R+'data/external/hpa/proteinatlas.tsv',sep='\t',usecols=['Gene','Ensembl','Secretome location','Blood concentration - Conc. blood IM [pg/L]','Blood concentration - Conc. blood MS [pg/L]'])
print('HPA rows',len(hpa),'unique Ensembl',hpa.Ensembl.nunique())
hpa=hpa.drop_duplicates('Ensembl').rename(columns={'Ensembl':'ens'})
g=g.merge(hpa,on='ens',how='left')
def cov(x):
    return pd.Series({'n':len(x),'LOEUF':x['lof.oe_ci.upper'].notna().sum(),'DepMap_entrez':x.c_e.notna().sum(),'DepMap_either':(x.c_e.notna()|x.c_s.notna()).sum(),
      'HPA':x.Gene.notna().sum(),'secretome':x['Secretome location'].notna().sum(),'blood_MS':x['Blood concentration - Conc. blood MS [pg/L]'].notna().sum(),
      'blood_IM':x['Blood concentration - Conc. blood IM [pg/L]'].notna().sum()})
t=g.groupby('next_experiment_rule_id').apply(cov); t.loc['ALL']=cov(g)
print(t.astype(int).to_string())
for nm in ['CEGv2','NEGv1']:
    h=pd.read_csv(R+f'data/public_r4/Hart2017_CEG2/{nm}.txt',sep='\t'); h['ENTREZ_ID']=h.ENTREZ_ID.astype('Int64')
    hm=h.merge(gcp.dropna(subset=['hgnc']).drop_duplicates('hgnc'),left_on='HGNC_ID',right_on='hgnc',how='left')
    print(nm,'n',len(h),'missing entrez',h.ENTREZ_ID.isna().sum(),'HGNC->ENSG',hm.ens.notna().sum(),'with MANE LOEUF',hm.ens.isin(L.gene_id).sum(),
          'DepMap by entrez',h.ENTREZ_ID.isin(dm.entrez).sum(),'in Catalog',hm.ens.isin(g.ens).sum())

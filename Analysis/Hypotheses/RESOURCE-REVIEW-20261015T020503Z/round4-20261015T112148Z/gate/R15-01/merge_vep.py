import json,glob,os,sys,pandas as pd
D=sys.argv[1]
PA={'missense_variant','stop_gained','stop_lost','start_lost','frameshift_variant','inframe_insertion','inframe_deletion','protein_altering_variant','splice_donor_variant','splice_acceptor_variant'}
rows=[]
for f in glob.glob(f'{D}/tmp/vep/*.json'):
    c,p=os.path.basename(f)[:-5].split('_'); j=json.load(open(f))
    gs=set()
    for r in j:
        for t in r.get('transcript_consequences',[]):
            if t.get('biotype')=='protein_coding' and PA & set(t['consequence_terms']): gs.add((t['gene_id'],t.get('gene_symbol','')))
    for g in gs or [('NA','NA')]: rows.append(dict(site=f'{c}:{p}',ensembl_id=g[0],symbol=g[1]))
v=pd.DataFrame(rows); print('VEP sites',v.site.nunique(),'without protein-coding PA gene',v[v.ensembl_id=='NA'].site.nunique())
F='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv'
m=pd.read_csv(F,sep='\t',dtype={'chromosome':str}); m=m[m.consequence_category.isin(['protein_altering','canonical_splice'])]
m['site']=m.chromosome+':'+m.position.astype(str)
a=pd.read_csv(f'{D}/out/coding_member_gene_map.tsv',sep='\t')[['site','ensembl_id','symbol']].drop_duplicates()
mp=pd.concat([a,v]).drop_duplicates()
x=m.merge(mp,on='site',how='left')
print('sites total',x.site.nunique(),'sites with no gene',x[x.ensembl_id.isna()|(x.ensembl_id=='NA')].site.nunique())
x=x[x.ensembl_id.notna()&(x.ensembl_id!='NA')]
g=x.groupby(['ensembl_id','symbol']).agg(max_pip=('normalized_susie_pip','max'),n_studies=('study','nunique'),n_sites=('site','nunique')).reset_index()
cat=pd.read_parquet('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/catalog_local/passport_gene_index.parquet')
g=g.merge(cat[['ensembl_id','primary_evidence_class']],on='ensembl_id',how='left')
g.sort_values('max_pip',ascending=False).to_csv(f'{D}/out/coding_gene_summary_full.tsv',sep='\t',index=False)
h=x[x.normalized_susie_pip>=0.5]; print('PIP>=0.5 sites',h.site.nunique(),'mapped')
for t in [0,0.1,0.5]:
    s=g[g.max_pip>=t]; print(f'genes max PIP>={t}: {len(s)}; in Catalog {s.primary_evidence_class.notna().sum()};',s.primary_evidence_class.value_counts(dropna=False).to_dict())
print(g[g.max_pip>=0.5].sort_values('max_pip',ascending=False).to_string())
for k in ['PNPLA3','TM6SF2','MTARC1','GCKR','APOE','SERPINA1','HSD17B13','MBOAT7','TRIB1','PPP1R3B','CIDEB']:
    r=g[g.symbol==k]; print(k, 'absent' if r.empty else r[['max_pip','n_studies','primary_evidence_class']].values.tolist())

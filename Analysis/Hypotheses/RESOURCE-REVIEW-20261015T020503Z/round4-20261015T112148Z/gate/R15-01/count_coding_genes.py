# R15-01 gate: map protein-altering / canonical-splice credible-set members to GENCODE v49 genes; count Catalog coverage.
import csv, sys, collections, json
import pandas as pd
D=sys.argv[1]
F='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv'
LIFT='/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/gwas_atac/variant_overlap_summary.csv'
GTF='/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/gencode.v49.primary_chr.gtf'
CAT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/catalog_local/passport_gene_index.parquet'
m=pd.read_csv(F,sep='\t',dtype={'chromosome':str})
m=m[m.consequence_category.isin(['protein_altering','canonical_splice'])].copy()
lift=pd.read_csv(LIFT,dtype={'chromosome':str})
lift=lift[['chromosome','position','chr_hg38','pos_hg38']].drop_duplicates(['chromosome','position'])
m=m.merge(lift,on=['chromosome','position'],how='left')
print('coding member rows',len(m),'distinct sites',m[['chromosome','position']].drop_duplicates().shape[0],'sites without hg38',m[m.pos_hg38.isna()][['chromosome','position']].drop_duplicates().shape[0])
sites=m.dropna(subset=['pos_hg38'])[['chr_hg38','pos_hg38']].drop_duplicates()
by=collections.defaultdict(list)
for c,p in sites.itertuples(index=False): by[c].append(int(p))
hits=collections.defaultdict(set)
with open(GTF) as fh:
    for line in fh:
        if line[0]=='#': continue
        f=line.split('\t',9)
        if f[0] not in by or f[2] not in ('CDS','exon'): continue
        a=f'gene_type "protein_coding"'
        if a not in f[8]: continue
        s,e=int(f[3]),int(f[4])
        for p in by[f[0]]:
            ok = (f[2]=='CDS' and s<=p<=e) or (f[2]=='exon' and (s-2<=p<s or e<p<=e+2))
            if ok:
                g=f[8].split('gene_name "')[1].split('"')[0]; gid=f[8].split('gene_id "')[1].split('"')[0].split('.')[0]
                hits[(f[0],p)].add((gid,g))
rows=[]
for r in m.dropna(subset=['pos_hg38']).itertuples():
    for gid,g in hits.get((r.chr_hg38,int(r.pos_hg38)),[('NA','NA')]):
        rows.append(dict(ensembl_id=gid,symbol=g,site=f'{r.chromosome}:{r.position}',study=r.study,tier=r.tier,trait=r.trait,pip=r.normalized_susie_pip,cat=r.consequence_category,csq=r.consequence))
x=pd.DataFrame(rows); x.to_csv(f'{D}/out/coding_member_gene_map.tsv',sep='\t',index=False)
print('site-gene rows',len(x),'sites unmapped to a protein-coding CDS/splice',x[x.ensembl_id=='NA'].site.nunique())
g=x[x.ensembl_id!='NA'].groupby(['ensembl_id','symbol']).agg(max_pip=('pip','max'),n_studies=('study','nunique'),n_sites=('site','nunique')).reset_index()
g.to_csv(f'{D}/out/coding_gene_summary.tsv',sep='\t',index=False)
cat=pd.read_parquet(CAT)
g['in_catalog']=g.ensembl_id.isin(cat.ensembl_id)
gc=g.merge(cat[['ensembl_id','primary_evidence_class']],on='ensembl_id',how='left')
for t in [0,0.1,0.5]:
    s=gc[gc.max_pip>=t]; print(f'genes max PIP>={t}: {len(s)}; in Catalog {s.in_catalog.sum()}; class counts',s.primary_evidence_class.value_counts(dropna=False).to_dict())
known_pos=['PNPLA3','TM6SF2','MTARC1','GCKR','APOE','SERPINA1','HSD17B13','PNPLA2','TMC4','APOH','CIDEB','GPAM','PPP1R3B','ADH1B']
known_neg=['MBOAT7','TRIB1','FADS1','SLC30A10','TOMM40']
print(gc[gc.symbol.isin(known_pos+known_neg)].to_string())
print('known_pos absent:',[k for k in known_pos if k not in set(gc.symbol)])
print('known_neg present:',[k for k in known_neg if k in set(gc.symbol)])
print('catalog genes',cat.ensembl_id.nunique())

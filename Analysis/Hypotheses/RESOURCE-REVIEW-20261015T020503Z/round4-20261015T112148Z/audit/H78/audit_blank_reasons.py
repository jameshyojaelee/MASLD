"""H78 audit 2: why Catalog genes lack a gnomAD Ensembl MANE row (counts by GENCODE gene_type and gnomAD row type). No control or route outcome."""
import re, gzip, pandas as pd, numpy as np
np.random.seed(20260923)
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
g=pd.read_parquet(R+'figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'); g['ens']=g.ensembl_id.str.split('.').str[0]
gn=pd.read_csv(R+'data/external/gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv',sep='\t',usecols=['gene','gene_id','transcript','canonical','mane_select','lof.oe_ci.upper'],low_memory=False)
e=gn[gn.gene_id.astype(str).str.startswith('ENSG')]
mane=e[e.mane_select.astype(str).str.lower().eq('true')]
can=e[e.canonical.astype(str).str.lower().eq('true')]
rows=[]
with gzip.open('/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz','rt') as f:
    for line in f:
        if line.startswith('#'): continue
        p=line.split('\t',9)
        if p[2]!='gene': continue
        a=p[8]; rows.append((re.search(r'gene_id "([^"]+)"',a).group(1).split('.')[0], re.search(r'gene_type "([^"]+)"',a).group(1), p[0]))
gc=pd.DataFrame(rows,columns=['ens','gtype','chr']).drop_duplicates('ens')
g=g.merge(gc,on='ens',how='left')
g['in_mane']=g.ens.isin(mane.gene_id); g['mane_loeuf']=g.ens.isin(mane.dropna(subset=['lof.oe_ci.upper']).gene_id)
g['any_ensg_row']=g.ens.isin(e.gene_id); g['canon_loeuf']=g.ens.isin(can.dropna(subset=['lof.oe_ci.upper']).gene_id)
g['pc']=g.gtype.eq('protein_coding')
print('Catalog', len(g), '| in Ensembl MANE rows', g.in_mane.sum(), '| MANE with LOEUF', g.mane_loeuf.sum(), '| MANE row but LOEUF missing', int((g.in_mane & ~g.mane_loeuf).sum()))
nm=g[~g.in_mane]
print('Catalog genes with no Ensembl MANE row', len(nm), '| protein_coding', int(nm.pc.sum()), '| non-coding/other', int((~nm.pc).sum()))
print('  of no-MANE: any gnomAD ENSG row', int(nm.any_ensg_row.sum()), '| canonical row with LOEUF', int(nm.canon_loeuf.sum()), '| protein_coding with canonical LOEUF', int((nm.pc & nm.canon_loeuf).sum()))
print('  no-MANE protein_coding by chromosome class:', nm[nm.pc].chr.str.match(r'^chr([0-9]+|X|Y|M)$').value_counts().to_dict())
print('gene_type of no-MANE (top 10):'); print(nm.gtype.value_counts().head(10).to_string())
print('by route: no-MANE protein_coding'); print(nm[nm.pc].groupby('next_experiment_rule_id').size().to_string())
print('by route: no-MANE protein_coding with canonical LOEUF'); print(nm[nm.pc & nm.canon_loeuf].groupby('next_experiment_rule_id').size().to_string())
print('MANE-row-but-no-LOEUF by route'); print(g[g.in_mane & ~g.mane_loeuf].groupby('next_experiment_rule_id').size().to_string())
print('POLR2A rows (ENSG):'); print(e[e.gene=='POLR2A'].to_string())

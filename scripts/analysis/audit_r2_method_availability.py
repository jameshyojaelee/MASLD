#!/usr/bin/env python3
"""Finite r2 method availability from existing master, never inferred from labels."""
import collections,csv,hashlib,json,math,os,sys
from pathlib import Path
ROOT=Path(os.environ['MASLD_PROJECT_ROOT']);RELEASE=ROOT/'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z';OUT=Path(os.environ['R2_AVAILABILITY_ROOT']);OUT.mkdir(parents=True,exist_ok=False)
master=RELEASE/'susie_coloc_all_gwas.csv';tierpath=ROOT/'GWAS/finemapping/config/gwas_trait_tier.tsv';regpath=ROOT/'GWAS/finemapping/config/gwas_registry.tsv'
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as fh:
  for b in iter(lambda:fh.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
inputs={str(p):sha(p) for p in [master,tierpath,regpath]}
with tierpath.open() as fh:tier={x['study_name']:x for x in csv.DictReader(fh,delimiter='\t')}
with regpath.open() as fh:registry={x['study_name']:x for x in csv.DictReader(fh,delimiter='\t')}
expected_master=next(x.split()[0] for x in (RELEASE/'output_sha256.tsv').read_text().splitlines() if x.split()[-1].endswith('/susie_coloc_all_gwas.csv'));assert inputs[str(master)]==expected_master
scopes={k:dict(count=collections.Counter(),sets=collections.defaultdict(set),states=collections.Counter()) for k in ['full50','tier12_main','tier1_main','tier2_main','tier3','tier4']};seen=set();study_gene=set();duplicate_study_gene=0; duplicate_groups=collections.defaultdict(list)
num=lambda x:float(x) if x not in ['', 'NA','NaN','nan'] else math.nan
with master.open() as fh:
 for row in csv.DictReader(fh):
  study=row['gwas_name'];assert study in tier and study in registry
  key=(study,row['chr'],row['ensembl']);assert key not in seen;seen.add(key)
  sg=(study,row['ensembl']);duplicate_study_gene+=sg in study_gene;study_gene.add(sg);duplicate_groups[sg].append(dict(chr=row['chr'],gene=row['gene'],method=row['method'],PP4_susie=row['PP.H4.susie'],PP4_abf=row['PP.H4.abf']))
  t=int(tier[study]['tier']);main=tier[study]['placement']=='main';which=['full50']
  if main and t in [1,2]:which+=['tier12_main',f'tier{t}_main']
  if t in [3,4]:which+=[f'tier{t}']
  abf=num(row['PP.H4.abf']);susie=num(row['PP.H4.susie']);af=math.isfinite(abf);sf=math.isfinite(susie);gene=row['ensembl'];valid=gene.startswith('ENSG');symbol=row['gene'].strip();method=row['method'];source=row['eligibility_source'];rerun=source=='rerun_shared_posterior_check';untestable=method=='susie_untestable_insufficient_shared_posterior'
  if untestable:assert not sf
  if rerun and sf:assert int(row['n_cs_pairs_eligible'])>0
  for name in which:
   s=scopes[name];c=s['count'];g=s['sets'];c['rows']+=1;c['missing_or_invalid_ensembl_rows']+=int(not valid);c['blank_symbol_rows']+=int(not symbol);g['studies'].add(study);g['all_stable_genes'].add(gene) if valid else None
   for label,flag in [('abf_finite',af),('susie_finite',sf),('abf_positive',af and abf>.5),('susie_positive',sf and susie>.5),('union_positive',(af and abf>.5)or(sf and susie>.5)),('rerun_finite',rerun and sf),('untestable_pair',untestable)]:
    c[label+'_pairs']+=int(flag)
    if flag and valid:g[label+'_genes'].add(gene)
    if flag and symbol:g[label+'_symbols'].add(symbol)
   c['not_rerun_negative_finite']+=int(source=='not_rerun_negative' and sf);c['not_rerun_negative_nonfinite']+=int(source=='not_rerun_negative' and not sf)
   c['rerun_positive']+=int(rerun and sf and susie>.5);c['rerun_tested_le_0.5']+=int(rerun and sf and susie<=.5);c['abf_fallback_pairs']+=int(method=='abf_fallback');c['abf_only_pairs']+=int(method=='abf_only')
   s['states'][(source,method,'finite' if sf else 'nonfinite')]+=1
rows=[];states=[]
for name,s in scopes.items():
 c=s['count'];g=s['sets'];r=dict(scope=name,studies=len(g['studies']),**dict(c))
 for key in ['abf_finite','susie_finite','abf_positive','susie_positive','union_positive','rerun_finite','untestable_pair']:r[key+'_genes']=len(g[key+'_genes']);r[key+'_symbols']=len(g[key+'_symbols'])
 r['susie_finite_fraction']=c['susie_finite_pairs']/c['rows'];r['abf_finite_fraction']=c['abf_finite_pairs']/c['rows'];r['untestable_genes_with_other_finite_susie']=len(g['untestable_pair_genes']&g['susie_finite_genes']);rows.append(r)
 for (source,method,finite),n in sorted(s['states'].items()):states.append(dict(scope=name,eligibility_source=source,method=method,susie_state=finite,rows=n))
full,primary=rows[:2];assert full['studies']==50 and primary['studies']==35
assert full['susie_positive_genes']==708 and primary['susie_positive_genes']==436 and primary['abf_positive_genes']==822 and primary['union_positive_genes']==996
assert full['untestable_pair_pairs']==25 and primary['untestable_pair_pairs']==11
for filename,data in [('method_availability.tsv',rows),('method_state_detail.tsv',states)]:
 fields=list(dict.fromkeys(k for r in data for k in r))
 with (OUT/filename).open('w') as fh:w=csv.DictWriter(fh,fields,delimiter='\t');w.writeheader();w.writerows(data)
assert inputs=={str(p):sha(p) for p in [master,tierpath,regpath]}
duplicates=[dict(study=st,ensembl=gene,valid_stable_id=gene.startswith('ENSG'),rows=values) for (st,gene),values in duplicate_groups.items() if len(values)>1];(OUT/'collapsed_study_gene_duplicates.json').write_text(json.dumps(duplicates,indent=2));report=dict(status='PASS',duplicate_collapsed_identifiers=sorted({x['ensembl'] for x in duplicates}),duplicate_collapsed_study_groups=len(duplicates),rows=len(seen),duplicate_study_gene_rows=duplicate_study_gene,input_sha256=inputs,finite_definition='math.isfinite(PP.H4.method), never eligibility_source or method label alone',eligibility_boundary='Only rerun_shared_posterior_check finite rows were freshly shared-posterior-qualified; finite carried negative rows were not rerun, so complete fresh eligibility availability is unavailable from this bounded rerun',scopes=rows,python=sys.version)
(OUT/'availability_check.json').write_text(json.dumps(report,indent=2));(OUT/'CHECKED.txt').write_text('PASS: all50studies/35primary; original master checksum unchanged; numeric finite availability separately from labels; explicit25full/11primary untestable pairs; no new fits.\n');print(json.dumps(rows,indent=2))

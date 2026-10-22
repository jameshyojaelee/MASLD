import csv, gzip, math, os, sys, subprocess, statistics as st
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
G=R+'/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG11'
L=R+'/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L06'
out=open(G+'/counts.txt','w')
def p(*a):
    print(*a,file=out,flush=True)
anch={}
for r in csv.DictReader(open(G+'/anchored_genes.tsv'),delimiter='\t'):
    anch[r['symbol']]=r['primary_evidence_class']
p('anchored+concordant genes',len(anch))
# COLOC r2 rows for UKBB enzyme traits among anchored genes
best={}  # gene -> (pp, trait, snp)
anyrow={}
with open(R+'/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv') as f:
    rd=csv.DictReader(f)
    for r in rd:
        if r['gwas_name'] not in ('UKBB_ALT','UKBB_AST','UKBB_GGT'): continue
        if r['gene'] not in anch: continue
        pp=r['PP.H4.susie'] if r['PP.H4.susie'] not in ('','NA') else r['PP.H4.abf']
        try: pp=float(pp)
        except: continue
        anyrow.setdefault(r['gene'],[]).append((pp,r['gwas_name'],r['top_snp']))
        if r['gene'] not in best or pp>best[r['gene']][0]: best[r['gene']]=(pp,r['gwas_name'],r['top_snp'])
p('anchored genes with any UKBB enzyme COLOC row',len(anyrow))
for thr in (0.5,0.8):
    p(f'anchored genes with UKBB enzyme PP.H4>={thr}',sum(1 for g,v in best.items() if v[0]>=thr))
# clumped leads
leads={}
for t in ('alt','ast','ggt'):
    leads[t]=[(l.split('\t')[0],int(l.split('\t')[1])) for l in open(f'{L}/ukbb_{t}_clumped500kb.tsv')]
allleads=[x for t in leads for x in leads[t]]
near=0
for g,v in best.items():
    c,pos=v[2].split(':')[:2]; pos=int(pos)
    if any(lc==c and abs(lp-pos)<=500000 for lc,lp in allleads): near+=1
p('anchored genes whose best UKBB top_snp lies within 500kb of a clumped UKBB enzyme lead',near)
# lookup F/M at best top SNP (hg19 reformatted files)
need={}
for g,v in best.items():
    t=v[1].split('_')[1]; c,pos=v[2].split(':')[:2]
    need.setdefault(t,{})[(c,pos)]=g
for t in need:
    se={}
    for sx in 'FM':
        fn=f'{R}/GWAS/finemapping/data/sumstats/PanUKBB_{sx}_{t}_reformatted_hg19.tsv'
        with open(fn) as f:
            next(f)
            for l in f:
                a=l.split('\t')
                k=(a[0],a[1])
                if k in need[t]:
                    try: se.setdefault(k,{})[sx]=float(a[5])
                    except: pass
    both=[math.sqrt(v['F']**2+v['M']**2) for v in se.values() if 'F' in v and 'M' in v]
    p(f'{t}: anchored genes with best trait {t}: {len(need[t])}; with F and M estimate at top SNP: {len(both)}; median SE(diff) {st.median(both) if both else "NA"}')
# chrX counts in Neale sex files (no z_diff computed)
S=R+'/GWAS/MR_Data/PanUKBB/sex_stratified'
for pc in ('30620','30650','30730'):
    for sx in ('female','male'):
        cmd=f"zcat {S}/{pc}_irnt.gwas.imputed_v3.{sx}.varorder.tsv.bgz | awk -F'\\t' '$1 ~ /^X:/ {{n++; if($4==\"false\" && $11!=\"NaN\" && $11+0<5e-8) g++}} END{{print n+0, g+0}}'"
        p(pc,sx,'chrX rows, chrX high-confidence p<5e-8:',subprocess.run(cmd,shell=True,capture_output=True,text=True).stdout.strip())
out.close()

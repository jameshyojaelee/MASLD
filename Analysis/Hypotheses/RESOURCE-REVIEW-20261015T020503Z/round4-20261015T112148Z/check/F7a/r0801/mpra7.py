# Independent re-derivation of Fig 7B reporter-assay rows: all 7 displayed rows vs allele-identity ridge.
import sys, csv, gzip, math
from pathlib import Path
import numpy as np
from scipy.stats import rankdata, norm
ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark')
sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'src'))
from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes
C=("HepG2_control","HepG2_PAOA"); S=(1103,2909,4721,6673,8111)
CANDS={('available_simple_controls','allele_identity_ridge'):'605',
 ('hyenadna','delta_ridge'):'605',('nucleotide_transformer','delta_ridge'):'605',('caduceus','delta_ridge'):'605',('dnabert2','delta_ridge'):'605',
 ('enformer_crested_restricted_port','rc_ensemble_sad_ridge'):'604',('enformer_crested_restricted_port','rc_ensemble_sar_ridge'):'604',('sei','native_40class_ridge'):'604'}
FIT={'605':'executions/model-cpu-train-605-21099008/fit/oof_predictions.tsv.gz','604':'executions/model-cpu-train-604-21097189/fit/oof_predictions.tsv.gz'}
vals={}; blk={}
for run,p in FIT.items():
  with gzip.open(ROOT/p,'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
      k=(r['model_id'],r['head_id'])
      if CANDS.get(k)!=run or r['assay_context_id'] not in C: continue
      vals.setdefault((k,int(r['seed']),r['assay_context_id']),{})[r['element_id']]=float(r['prediction'])
      blk[r['element_id']]=r['long_range_block_id']
E=sorted(blk); print('elements',len(E),'blocks',len(set(blk.values())))
out=load_outcomes(ROOT/'executions/gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz',set(E))
obs={c:np.array([out[(e,c)]['mean'] for e in E]) for c in C}
ens={k:{c:np.mean([[vals[(k,s,c)][e] for e in E] for s in S],axis=0) for c in C} for k in CANDS}
def sp(a,b):
  l=rankdata(a);r=rankdata(b);l-=l.mean();r-=r.mean();return float((l*r).sum()/math.sqrt((l*l).sum()*(r*r).sum()))
def mac(v): return float(np.tanh(np.mean(np.arctanh(np.clip(v,-.999999,.999999)))))
pt={k:mac([sp(obs[c],ens[k][c]) for c in C]) for k in CANDS}
bl=np.array([blk[e] for e in E]); U=sorted(set(bl)); byb={b:np.flatnonzero(bl==b) for b in U}
rng=np.random.default_rng(20260825); B=10000; d={k:np.empty(B) for k in CANDS}
for i in range(B):
  s=rng.integers(0,len(U),size=len(U)); ii=np.concatenate([byb[U[j]] for j in s])
  for k in CANDS: d[k][i]=mac([sp(obs[c][ii],ens[k][c][ii]) for c in C])
ref=('available_simple_controls','allele_identity_ridge'); a=0.05/7; z=norm.ppf(1-a/2)
print('ref',pt[ref])
w=open(sys.argv[1],'w'); w.write('model\thead\tscore\tgain\tlo95\thi95\tbonf7_pct_lo\tbonf7_pct_hi\tbonf7_normal_lo\tfrac_le0\n')
for k in CANDS:
  if k==ref: continue
  g=d[k]-d[ref]; G=pt[k]-pt[ref]; q=np.quantile(g,[.025,.975,a/2,1-a/2])
  w.write(f"{k[0]}\t{k[1]}\t{pt[k]:.4f}\t{G:.4f}\t{q[0]:.4f}\t{q[1]:.4f}\t{q[2]:.4f}\t{q[3]:.4f}\t{G-z*g.std(ddof=1):.4f}\t{np.mean(g<=0):.4f}\n")
w.close()

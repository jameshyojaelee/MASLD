# Re-derive MPRA gain bootstrap (seed 20260825, 10000 long-range-block resamples) for
# hyenadna / nucleotide_transformer / caduceus / dnabert2 delta_ridge vs allele_identity_ridge,
# then report 95% and Bonferroni (alpha 0.05/7) percentile intervals. Read-only on repo files.
import sys, csv, gzip, math
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark')
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'src'))
from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes
CONTEXTS = ("HepG2_control", "HepG2_PAOA"); SEEDS=(1103,2909,4721,6673,8111)
rows = list(csv.DictReader(open(ROOT/'executions/model-check-220-21088696/contract/row_universe.tsv'), delimiter='\t'))
elements = sorted({r['element_id'] for r in rows}); idx={e:i for i,e in enumerate(elements)}
meta={}
for r in rows: meta.setdefault(r['element_id'], r['long_range_block_id'])
blocks=np.asarray([meta[e] for e in elements]); print('elements',len(elements),'blocks',len(set(blocks)))
out = load_outcomes(ROOT/'executions/gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz', set(elements))
obs={c:np.asarray([out[(e,c)]['mean'] for e in elements]) for c in CONTEXTS}
CANDS=[('available_simple_controls','allele_identity_ridge'),('hyenadna','delta_ridge'),('nucleotide_transformer','delta_ridge'),('caduceus','delta_ridge'),('dnabert2','delta_ridge')]
vals={(m,h,s,c):np.full(len(elements),np.nan) for m,h in CANDS for s in SEEDS for c in CONTEXTS}
with gzip.open(ROOT/'executions/model-cpu-train-605-21099008/fit/oof_predictions.tsv.gz','rt') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        k=(r['model_id'],r['head_id'],int(r['seed']),r['assay_context_id'])
        if k in vals: vals[k][idx[r['element_id']]]=float(r['prediction'])
assert all(np.isfinite(v).all() for v in vals.values())
ens={cand:{c:np.mean(np.vstack([vals[(*cand,s,c)] for s in SEEDS]),axis=0) for c in CONTEXTS} for cand in CANDS}
def sp(a,b):
    l=rankdata(a);r=rankdata(b);l-=l.mean();r-=r.mean();d=math.sqrt((l*l).sum()*(r*r).sum());return float((l*r).sum()/d)
def macro(v):
    v=np.clip(np.asarray(v),-0.999999,0.999999);return float(np.tanh(np.mean(np.arctanh(v))))
pts={cand:macro([sp(obs[c],ens[cand][c]) for c in CONTEXTS]) for cand in CANDS}
uniq=sorted(set(blocks.tolist())); byb={b:np.flatnonzero(blocks==b) for b in uniq}
rng=np.random.default_rng(20260825); B=10000
samp={cand:np.full(B,np.nan) for cand in CANDS}
for it in range(B):
    s=rng.integers(0,len(uniq),size=len(uniq)); ii=np.concatenate([byb[uniq[j]] for j in s])
    for cand in CANDS: samp[cand][it]=macro([sp(obs[c][ii],ens[cand][c][ii]) for c in CONTEXTS])
ref=CANDS[0]; a=0.05/7
print('reference point',pts[ref])
for cand in CANDS[1:]:
    g=samp[cand]-samp[ref]
    lo95,hi95=np.quantile(g,(0.025,0.975)); lob,hib=np.quantile(g,(a/2,1-a/2))
    print(f"{cand[0]}\tgain={pts[cand]-pts[ref]:.4f}\t95%=[{lo95:.4f},{hi95:.4f}]\tBonf7(99.29%)=[{lob:.4f},{hib:.4f}]\tfrac_draws<=0={np.mean(g<=0):.4f}\tsd={g.std(ddof=1):.4f}")

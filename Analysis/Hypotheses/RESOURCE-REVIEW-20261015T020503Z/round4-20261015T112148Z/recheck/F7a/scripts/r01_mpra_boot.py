"""Independent re-derivation of Fig 7B MPRA gains and Bonferroni lower bounds.
Own ensemble (pandas), own Spearman (scipy), own bootstrap loop; outcome definition from the
release loader. Runs the release bootstrap seed (20260825) and an independent seed (20260923)."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
BENCH = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark")
sys.path.insert(0, str(BENCH)); sys.path.insert(0, str(BENCH / "src"))
from scripts.evaluate_gse281364_dna_lm_common_lane import CONTEXTS, load_outcomes
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
B = 10000
cols = ["seed","assay_context_id","element_id","long_range_block_id","model_id","head_id","prediction"]
frames = []
for run in ["605-21099008", "604-21097189"]:
    d = pd.read_csv(BENCH/f"executions/model-cpu-train-{run}/fit/oof_predictions.tsv.gz", sep="\t", usecols=cols)
    d["run"] = run[:3]; frames.append(d)
d = pd.concat(frames)
print("contexts in OOF:", d.assay_context_id.value_counts().to_dict())
d = d[d.assay_context_id.isin(CONTEXTS)]
d = d[~((d.run=="604") & (d.model_id=="available_simple_controls"))]
print("seeds per head:", d.groupby(["model_id","head_id"]).seed.nunique().to_dict())
ens = d.groupby(["model_id","head_id","assay_context_id","element_id"]).prediction.mean().unstack("element_id")
elements = sorted(d.element_id.unique())
blkmap = d.drop_duplicates("element_id").set_index("element_id").long_range_block_id
blk = blkmap.loc[elements].values
print("elements", len(elements), "blocks", len(set(blk)))
oc = load_outcomes(BENCH/"executions/gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz", set(elements))
obs = {c: np.array([oc[(e,c)]["mean"] for e in elements]) for c in CONTEXTS}
heads = sorted(set(map(tuple, ens.index.droplevel(2).unique())))
heads = [h for h in heads if h[0] != "available_simple_controls" or h[1]=="allele_identity_ridge"]
REF = ("available_simple_controls","allele_identity_ridge")
pred = {h: {c: ens.loc[(h[0],h[1],c), elements].values.astype(float) for c in CONTEXTS} for h in heads}
def macro(rs): return float(np.tanh(np.mean(np.arctanh(np.clip(rs,-0.999999,0.999999)))))
point = {h: macro([spearmanr(obs[c], pred[h][c]).statistic for c in CONTEXTS]) for h in heads}
print("REF point", point[REF])
# ranks precomputed not possible under resampling with ties; compute each draw
ub = sorted(set(blk)); idx_by = {b: np.flatnonzero(blk==b) for b in ub}
def sp(a,b):
    ra = rankdata(a); rb = rankdata(b); ra-=ra.mean(); rb-=rb.mean()
    return float((ra*rb).sum()/np.sqrt((ra**2).sum()*(rb**2).sum()))
F7 = [("caduceus","delta_ridge"),("dnabert2","delta_ridge"),("hyenadna","delta_ridge"),("nucleotide_transformer","delta_ridge"),
      ("enformer_crested_restricted_port","rc_ensemble_sad_ridge"),("enformer_crested_restricted_port","rc_ensemble_sar_ridge"),("sei","native_40class_ridge")]
F11 = [h for h in heads if h != REF]
rows=[]
for seed in [20260825, 20260923]:
    rng = np.random.default_rng(seed)
    draws = {h: np.empty(B) for h in heads}
    for i in range(B):
        s = rng.integers(0, len(ub), size=len(ub))
        ix = np.concatenate([idx_by[ub[j]] for j in s])
        for h in heads:
            draws[h][i] = macro([sp(obs[c][ix], pred[h][c][ix]) for c in CONTEXTS])
    for fam, H in (("F7",F7),("F11",F11)):
        m=len(H)
        for h in H:
            g = draws[h]-draws[REF]
            rows.append(dict(boot_seed=seed, family=fam, m=m, model=h[0], head=h[1], gain=point[h]-point[REF],
                lo95=np.quantile(g,0.025), hi95=np.quantile(g,0.975), bonf_lo=np.quantile(g,0.025/m),
                frac_le0=float(np.mean(g<=0))))
pd.DataFrame(rows).to_csv(OUT/"mpra_recheck.tsv", sep="\t", index=False, float_format="%.6f")
print(pd.DataFrame(rows).to_string())

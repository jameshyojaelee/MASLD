"""Record the prespecified 1000-draw matched-pair label control.

Controls section2 relabels matched pairs. Swap genetic/state labels within each
realized pair, retaining its expression-matched members and outcome measurements.
This diagnostic does not recalibrate the invalid class relabel/rematch null.
"""
import hashlib,json,os
from pathlib import Path
import numpy as np
import pandas as pd
source=Path(os.environ["CLASS_COMPARISON_ROOT"])
out=Path(os.environ["CLASS_NEGATIVE_CONTROL_ROOT"])
if out.exists():raise SystemExit("Refusing existing output directory")
out.mkdir(parents=True)
SEED=20260907;N_NEG=1000
prespec=source/"prespec/PRESPEC.md"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(prespec)=="f5213ae6cf23db5dbf6ef9cd8f0c3077101dca22983fa7ba05b00bf8469c4ff5"
manifest={str(prespec):sha(prespec)};results=[]
for arm in ["singlecell_hep_specificity","protein_stage_response_residualised"]:
    path=source/"results"/f"pairs_{arm}.tsv";manifest[str(path)]=sha(path)
    p=pd.read_csv(path,sep="\t")
    assert p.g_genetic.is_unique and p.g_state.is_unique
    assert not set(p.g_genetic)&set(p.g_state)
    d=(p.y_genetic-p.y_state).to_numpy(float)
    assert np.isfinite(d).all()
    rng=np.random.default_rng(SEED)
    # +/-1 implements exactly the two possible label orders of a fixed pair.
    labels=rng.choice([-1.,1.],size=(N_NEG,len(d)))
    null=(labels*d).mean(axis=1)
    np.save(out/f"negative_control_draws_{arm}.npy",null)
    sd=float(null.std(ddof=1));lo,hi=np.quantile(null,[.025,.975])
    results.append(dict(arm=arm,n_pairs=len(d),seed=SEED,n_draws=N_NEG,
        observed_paired_mean=float(d.mean()),negative_control_mean=float(null.mean()),
        negative_control_sd=sd,monte_carlo_se_of_mean=sd/np.sqrt(N_NEG),
        analytical_expected_mean=0.,central95_lower=float(lo),central95_upper=float(hi),
        interpretation="centering_diagnostic_no_new_pass_threshold",
        primary_relabel_rematch_calibrated=False,
        gene_dependence_adjusted=False))
    assert sha(path)==manifest[str(path)]
pd.DataFrame(results).to_csv(out/"negative_control_summary.tsv",sep="\t",index=False)
(out/"inputs_sha256.json").write_text(json.dumps(manifest,indent=2))
(out/"method.txt").write_text("Literal matched-pair label relabeling control: 1000 draws per primary outcome, seed 20260907; realized pairs fixed, identical mean genetic-minus-state statistic. No model refits or new matching. Mean and Monte Carlo SE reported without inventing a pass tolerance. This does not validate class-expression exchangeability or adjust gene dependence.\n")
print(pd.DataFrame(results).to_string(index=False))
(out/"COMPLETE.txt").write_text("Two prespecified matched-pair negative-control summaries completed.\n")

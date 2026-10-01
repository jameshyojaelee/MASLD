#!/usr/bin/env python3
"""Independently reconstruct matched sweep performance from deposited predictions.

Only completed fixed-exposure runs on identical validation identities enter a
paired comparison. One split/seed remains development evidence. Native
AlphaGenome is admitted only after the complete population attempt accounting
exists, and its calibration uses the same permitted training chromosomes.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge

PROJ=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(PROJ/"scripts/analysis/alphagenome_program"))
import i1_common as C


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    runs=json.loads((args.sweep/"runs.json").read_text())
    successful=[r for r in runs if r.get("status") == "completed_fixed_exposure" and r.get("steps_requested") == 5000]
    disposition=[r for r in runs if r not in successful]
    if not successful:
        raise ValueError("No full-exposure models to compare")
    labels=pd.read_csv(C.C2_LABELS,sep="\t").set_index("lead_variant_id")
    predictions=None
    compute={}
    for run in successful:
        directory=Path(run["output"])
        table=pd.read_csv(directory/"validation_predictions.tsv",sep="\t").set_index("variant_id").sort_index()
        if table.index.duplicated().any():
            raise ValueError("Duplicate variant predictions")
        if predictions is None:
            predictions=pd.DataFrame(index=table.index)
            predictions["observed_beta"]=table.observed_beta
        if not predictions.index.equals(table.index):
            raise ValueError("Completed arms have unequal validation identities")
        np.testing.assert_allclose(predictions.observed_beta,table.observed_beta,atol=1e-7)
        predictions[run["name"]]=table.predicted_beta
        cost=json.loads((directory/"feasibility.json").read_text())
        compute[run["name"]]=cost
    joined=labels.loc[predictions.index]
    np.testing.assert_allclose(predictions.observed_beta,joined.beta_alt,atol=1e-7)
    if not (joined.heldout_fold == 1).all():
        raise ValueError("Unexpected validation fold")
    train=labels.loc[~labels.heldout_fold.isin([0,1])]
    bases={b:i for i,b in enumerate("ACGT")}
    def identity(frame):
        # Signed one-hot allele difference preserves same/swap identities.
        eye=np.eye(4)
        return np.stack([eye[bases[a]]-eye[bases[r]] for r,a in zip(frame.ref,frame.alt)])
    model=Ridge(alpha=1000,fit_intercept=False).fit(identity(train),train.beta_alt)
    predictions["allele_identity"]=model.predict(identity(joined))
    predictions["zero_effect"]=0.0
    native_status="not_requested"
    if args.native:
        if not (args.native/"completion.json").exists():
            native_status="full_population_attempt_accounting_incomplete"
        else:
            scores=pd.read_csv(args.native/"native_1048576.tsv",sep="\t").set_index("key")
            scores.index="chr"+scores.index
            common=predictions.index.intersection(scores.index)
            ntrain=train.index.intersection(scores.index)
            native_table=predictions.loc[common].copy()
            calibration=[]
            for output in ("atac","dnase"):
                x=scores.loc[ntrain,f"local_{output}_liver"].to_numpy()
                y=train.loc[ntrain,"beta_alt"].to_numpy()
                slope=float(np.dot(x,y)/max(np.dot(x,x),1e-20))
                native_table[f"native_{output}_calibrated"]=slope*scores.loc[common,f"local_{output}_liver"]
                calibration.append({"output":output,"training_n":len(ntrain),"slope":slope,"intercept":0})
            native_table.to_csv(args.out/"native_matched_predictions.tsv.gz",sep="\t",index_label="variant_id")
            pd.DataFrame(calibration).to_csv(args.out/"native_training_calibration.tsv",sep="\t",index=False)
            native_status=f"matched_subset_{len(common)}_of_{len(predictions)}; separate_table"
    metrics=[]
    for name in predictions.columns.drop("observed_beta"):
        y,p=predictions.observed_beta.to_numpy(),predictions[name].to_numpy()
        if not np.isfinite(p).all():
            raise ValueError("Nonfinite prediction entered matched comparison")
        norm=np.dot(p,p)
        cost=compute.get(name,{})
        metrics.append({"arm":name,"n":len(y),"locus_blocks":joined.block_1mb.nunique(),
            "RMSE_beta_units":float(np.sqrt(np.mean((p-y)**2))),"MAE_beta_units":float(np.mean(np.abs(p-y))),
            "signed_spearman":float(spearmanr(y,p).statistic) if np.std(p) else np.nan,
            "diagnostic_calibration_slope":float(np.dot(p,y)/norm) if norm else np.nan,
            "calibration_slope_use":"diagnostic_only_predictions_unchanged",
            "seconds":cost.get("elapsed_seconds_including_restore_checks_evaluation_checkpoint"),
            "training_steps":cost.get("completed_steps"),"trainable_parameters":cost.get("trainable_parameter_count"),
            "evidence":"single_split_seed_development"})
    predictions["block_1mb"]=joined.block_1mb
    predictions.to_csv(args.out/"matched_predictions.tsv.gz",sep="\t")
    pd.DataFrame(metrics).to_csv(args.out/"performance.tsv",sep="\t",index=False)
    (args.out/"disposition.json").write_text(json.dumps({"other_runs":disposition,"native_comparison":native_status,
        "protected_outcomes_read":False,"five_seed_nested_finalists":"not_run",
        "uncertainty":"paired_locus_uncertainty_not_run_here","native_function":"per_run_drift_not_measured_accuracy"},indent=2)+"\n")
    print(pd.DataFrame(metrics).to_string(index=False))


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sweep",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--native",type=Path)
    main(p.parse_args())

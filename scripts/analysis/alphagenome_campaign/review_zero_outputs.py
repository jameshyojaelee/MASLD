#!/usr/bin/env python3
"""Review zero-adapter output agreement without loading or fitting a model.

Scalar-output differences are independently reconstructed from deposited rows.
Pooled-feature differences and unchanged parameter/state assertions are checked
as GPU-produced summaries; their raw arrays are not deposited by this audit.
A completed GPU process is not a successful numerical-agreement result.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROUTES = ("fixed", "adapter_r4_last3", "adapter_r4_last5", "adapter_r16_last3", "adapter_r16_last5")
POOLINGS = ("variant", "symmetric")
SCALARS = ("production_effect", "diagnostic_effect", "ref_score", "alt_score")
OUTPUTS = (*SCALARS, "ref_pooled", "alt_pooled", "signed_pooled")


def boolean(value):
    if isinstance(value, (bool,np.bool_)):
        return bool(value)
    if str(value).lower() in ("true","false"):
        return str(value).lower() == "true"
    raise ValueError(f"Invalid Boolean annotation {value!r}")


def comparison(reference, observed):
    x,y = np.asarray(reference,dtype=float),np.asarray(observed,dtype=float)
    assert np.isfinite(x).all() and np.isfinite(y).all()
    delta = y-x
    return dict(max_absolute=float(np.max(np.abs(delta))), RMS_absolute=float(np.sqrt(np.mean(delta**2))),
        relative_L2=float(np.linalg.norm(delta)/max(np.linalg.norm(x),1e-20)), exact=np.array_equal(x,y))


def review(directory, out):
    design = json.loads((directory/"design.json").read_text())
    complete = json.loads((directory/"complete.json").read_text())
    examples = pd.read_csv(directory/"sequence_examples.tsv",sep="\t")
    effects = pd.read_csv(directory/"effect_outputs.tsv",sep="\t")
    agreement = pd.read_csv(directory/"route_agreement.tsv",sep="\t")
    repeat = json.loads((directory/"repeatability.json").read_text())
    assert design["variants"] == complete["variants"] == 32
    assert design["training_steps"] == complete["training_steps"] == 0
    assert complete["labels_read"] is False and complete["frozen_parameters_state_unchanged"] is True
    assert complete["routes"] == 5 and design["microbatch"] == 4 and design["length"] == 2048
    assert design["persistent_compilation_cache_enabled"] is False
    assert len(examples) == 32 and not examples.key.duplicated().any()
    assert examples.heldout_fold.isin([2,3,4]).all()
    assert not set(examples.columns).intersection({"beta_alt","beta_source","effect","pvalue","p_nominal"})
    assert len(effects) == 320 and not effects.duplicated(["route","pooling","key"]).any()
    assert set(effects.route) == set(ROUTES) and set(effects.pooling) == set(POOLINGS)
    assert np.isfinite(effects[list(SCALARS)].to_numpy()).all()
    expected = {(route,pool,output) for route in ROUTES for pool in POOLINGS for output in OUTPUTS}
    by_repeat = {(r["route"],r["pooling"],r["output"]):r for r in repeat}
    assert len(repeat) == len(by_repeat) == 70 and set(by_repeat) == expected
    for key,row in by_repeat.items():
        for repetition in ("same_compilation","fresh_compilation"):
            values = row[repetition]
            assert np.isfinite([values[c] for c in ("max_absolute","RMS_absolute","relative_L2")]).all()
            assert all(values[c] >= 0 for c in ("max_absolute","RMS_absolute","relative_L2"))
            if boolean(values["exact"]):
                assert values["max_absolute"] == values["RMS_absolute"] == values["relative_L2"] == 0
        calculated = max(row["same_compilation"]["max_absolute"],row["fresh_compilation"]["max_absolute"])
        np.testing.assert_allclose(row["repeat_floor_max_absolute"],calculated,atol=0,rtol=0)
    agreement = agreement.set_index(["route","pooling","output"])
    assert not agreement.index.duplicated().any()
    assert set(agreement.index) == {key for key in expected if key[0] != "fixed"}
    numeric = ["max_absolute","RMS_absolute","relative_L2","strict_float32_tolerance","repeat_floor","repeat_based_tolerance"]
    assert np.isfinite(agreement[numeric].to_numpy()).all() and (agreement[numeric] >= 0).all().all()
    reconstructed, instrumentation = [], []
    for pooling in POOLINGS:
        reference = effects.loc[(effects.route == "fixed")&(effects.pooling == pooling)].set_index("key").loc[examples.key]
        for route in ROUTES:
            observed = effects.loc[(effects.route == route)&(effects.pooling == pooling)].set_index("key")
            assert len(observed) == 32 and set(observed.index) == set(examples.key)
            observed = observed.loc[examples.key]
            # Diagnostic subtraction is explicitly float32; production and
            # diagnostic graphs may differ and are not forced into agreement.
            native_difference = observed.alt_score.to_numpy(dtype=np.float32)-observed.ref_score.to_numpy(dtype=np.float32)
            np.testing.assert_array_equal(observed.diagnostic_effect.to_numpy(dtype=np.float32),native_difference)
            instrument = comparison(observed.production_effect,observed.diagnostic_effect)
            instrumentation.append(dict(route=route,pooling=pooling,**instrument))
            if route == "fixed":
                continue
            for output in SCALARS:
                calculated = comparison(reference[output],observed[output])
                row = agreement.loc[(route,pooling,output)]
                for column in ("max_absolute","RMS_absolute","relative_L2"):
                    np.testing.assert_allclose(row[column],calculated[column],atol=1e-12,rtol=1e-9)
                assert boolean(row["exact"]) == bool(calculated["exact"])
                strict = 8*np.finfo(np.float32).eps*max(1,float(np.max(np.abs(reference[output]))))
                np.testing.assert_allclose(row.strict_float32_tolerance,strict,atol=1e-14,rtol=1e-12)
                reconstructed.append(dict(route=route,pooling=pooling,output=output,**calculated,strict_tolerance=strict))
    strict_pass, floor_pass = [], []
    for (route,pooling,output),row in agreement.iterrows():
        floor = max(by_repeat[("fixed",pooling,output)]["repeat_floor_max_absolute"],by_repeat[(route,pooling,output)]["repeat_floor_max_absolute"])
        np.testing.assert_allclose(row.repeat_floor,floor,atol=0,rtol=0)
        np.testing.assert_allclose(row.repeat_based_tolerance,floor+row.strict_float32_tolerance,atol=1e-14,rtol=1e-12)
        strict = row.max_absolute <= row.strict_float32_tolerance
        within_floor = row.max_absolute <= row.repeat_based_tolerance
        assert boolean(row.strict_output_agreement) == strict
        assert boolean(row.within_measured_repeat_floor) == within_floor
        strict_pass.append(bool(strict));floor_pass.append(bool(within_floor))
    assert complete["all_strict_output_agreement"] == all(strict_pass)
    assert complete["all_within_measured_repeat_floor"] == all(floor_pass)
    nonzero_floor = any(r["repeat_floor_max_absolute"] > 0 for r in repeat)
    assert complete["nonzero_repeat_floor"] == nonzero_floor
    pd.DataFrame(reconstructed).to_csv(out/"scalar_agreement_reconstructed.tsv",sep="\t",index=False)
    pd.DataFrame(instrumentation).to_csv(out/"production_vs_diagnostic.tsv",sep="\t",index=False)
    agreement.loc[~agreement.strict_output_agreement.map(boolean)].to_csv(out/"strict_agreement_failures.tsv",sep="\t")
    return dict(status="pass" if all(strict_pass) else "failed_output_agreement", artifact_consistency="pass",
        scalar_contrasts_independently_reconstructed=len(reconstructed), all_comparison_summaries=len(agreement),
        strict_failures=sum(not flag for flag in strict_pass), beyond_repeat_floor=sum(not flag for flag in floor_pass),
        nonzero_measured_repeat_floor=nonzero_floor, variants_per_route_pool=32, routes=5, poolings=2,
        training_steps=0, protected_outcomes_read=False,
        parameter_state_check="GPU_assertion_declared_unchanged; parameter_arrays_not_reloaded_by_CPU_review",
        feature_contrast_check="aggregate_consistency_only; raw_pooled_arrays_not_deposited",
        interpretation="a_nonzero_repeat_floor_or_successful_process_is_not_proof_that_small_adaptation_gains_are_robust",
        scope="these32sequences_common_head_and_compilation_routes_only")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run numerical review on a compute node")
    args.out.mkdir(parents=True,exist_ok=False)
    try:
        result=review(args.audit,args.out)
    except Exception as error:
        result=dict(status="unresolved_or_invalid_audit",reason=repr(error),protected_outcomes_read=False)
        (args.out/"checks.json").write_text(json.dumps(result,indent=2)+"\n")
        raise
    (args.out/"checks.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    if result["status"] != "pass":
        raise SystemExit(2)


if __name__ == "__main__":
    main()

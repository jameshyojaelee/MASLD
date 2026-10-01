#!/usr/bin/env python3
"""Scientific invariants for the bounded native-plus-sequence comparison.

Run on a compute node using rnaseq Python, which supplies scikit-learn for an
independent primal ridge implementation. No real molecular outcomes are read.
"""
import argparse
import json
from pathlib import Path
import tempfile
import time
from types import SimpleNamespace

import numpy as np
import sklearn
from sklearn.linear_model import Ridge

import model_native_residual as M


def same_models(a, b):
    if set(a) != set(b):
        raise AssertionError("Model identities changed")
    for name in a:
        for key in a[name]:
            if key == "metadata":
                if a[name][key] != b[name][key]:
                    raise AssertionError("Training-only metadata changed")
            elif not np.array_equal(a[name][key], b[name][key]):
                raise AssertionError(f"Training-only state changed: {name}/{key}")


def synthetic_nested(xs, z, y, folds):
    """Use the production selection and fitting functions on six small views."""
    native = M.native_cv(z, y, folds)
    candidates, sufficient = [], {}
    for (length, pool), x in xs.items():
        cc, _, ss = M.feature_cv(x, z, y, folds, length, pool)
        candidates.extend(cc)
        sufficient[length, pool] = ss
    recipes = M.choose_recipes(candidates)
    train = np.isin(folds, M.TRAIN)
    native_fits = M.native_candidates(z[train], y[train])
    models = {"native_selected": native_fits[native["selected"]],
              "native_ad_ridge": native_fits[native["ridge"]]}
    outer = {}
    for name, recipe in recipes.items():
        key = recipe["length"], recipe["pooling"]
        if key not in outer:
            outer[key] = M.feature_path(sufficient[key], {"length": key[0], "pooling": key[1]})
        models[name] = outer[key][recipe["family"], recipe["penalty"]]
    choices = {"recipes": recipes, "native": native["selected"], "ridge": native["ridge"],
               "native_MSE": native["scores"], "margin": 0.05*native["scores"][native["selected"]]}
    return choices, models


def run(out):
    tick = time.monotonic()
    rng = np.random.default_rng(M.SEED)
    folds = np.repeat([1, 2, 3, 4], 24)
    z = rng.normal(size=(len(folds), 2))
    x0 = rng.normal(size=(len(folds), 11)).astype(np.float32)
    xs = {(length, pool): x0 + np.float32(0.02*j)*rng.normal(size=x0.shape).astype(np.float32)
          for j, (length, pool) in enumerate((l, p) for l in M.LENGTHS for p in M.POOLS)}
    y = 0.6*z[:, 0] - 0.15*z[:, 1] + 0.35*x0[:, 3] + rng.normal(0, 0.15, len(folds))
    train, held = folds != 1, folds == 1
    checks = []

    selection, fitted = synthetic_nested(xs, z, y, folds)
    corrupted = y.copy()
    corrupted[held] = rng.normal(1e8, 1e7, held.sum())
    other_selection, other_fitted = synthetic_nested(xs, z, corrupted, folds)
    assert selection == other_selection
    same_models(fitted, other_fitted)
    checks.append({"name": "held_development_target_perturbation", "pass": True,
                   "changes": "All fold-1 target values replaced; all inner choices, preprocessing and coefficients unchanged"})

    other_x = {k: value.copy() for k, value in xs.items()}
    other_z = z.copy()
    other_z[held] += 999
    for value in other_x.values():
        value[held] *= 1000
    changed_selection, changed_fitted = synthetic_nested(other_x, other_z, y, folds)
    assert selection == changed_selection
    same_models(fitted, changed_fitted)
    checks.append({"name": "held_development_input_preprocessing_independence", "pass": True})

    x = x0.copy()
    x[train, 0] = 0
    x[held, 0] = 1000
    x[train, 1] *= np.float32(1e-12)
    path = M.feature_path(M.stats(x[train], z[train], y[train]))
    maximum = 0.
    for (family, penalty), model in path.items():
        assert not model["feature_admitted"][:2].any()
        assert np.count_nonzero(model["feature_coef_raw"][:2]) == 0
        xx = x[train].astype(float)/model["feature_scale"]
        zz = z[train]/model["native_scale"]
        xx[:, ~model["feature_admitted"]] = 0
        zz[:, ~model["native_admitted"]] = 0
        xt = x[held].astype(float)/model["feature_scale"]
        zt = z[held]/model["native_scale"]
        xt[:, ~model["feature_admitted"]] = 0
        zt[:, ~model["native_admitted"]] = 0
        if family == "joint":
            direct = Ridge(alpha=penalty*train.sum(), fit_intercept=False, solver="cholesky")
            direct.fit(np.column_stack([xx, zz]), y[train])
            expected = direct.predict(np.column_stack([xt, zt]))
        else:
            if family == "residual":
                a = np.linalg.lstsq(zz, y[train], rcond=1e-12)[0]
                target = y[train]-zz@a
                baseline = zt@a
            else:
                target, baseline = y[train], np.zeros(held.sum())
            direct = Ridge(alpha=penalty*train.sum(), fit_intercept=False, solver="cholesky").fit(xx, target)
            expected = baseline+direct.predict(xt)
        observed = M.predict(model, x[held], z[held])
        difference = float(np.max(abs(expected-observed)))
        maximum = max(maximum, difference)
        np.testing.assert_allclose(observed, expected, rtol=1e-8, atol=1e-9)
    checks.append({"name": "independent_sklearn_primal_ridge_all_18_paths", "pass": True,
                   "max_absolute_prediction_difference": maximum, "sklearn": sklearn.__version__})
    checks.append({"name": "zero_and_tiny_training_RMS_excluded_at_inference", "pass": True})

    zero_z = z.copy()
    zero_z[train, 1] = 0
    zero_native = M.native_candidates(zero_z[train], y[train])
    zero_path = M.feature_path(M.stats(x[train], zero_z[train], y[train]))
    assert all(m["native_coef_raw"][1] == 0 for m in [*zero_native.values(), *zero_path.values()])
    assert all(np.isfinite(M.predict(m, x[held], zero_z[held])).all() for m in zero_path.values())
    checks.append({"name": "zero_native_channel_safe_OLS_and_ridge", "pass": True})

    reload_max = 0.
    with tempfile.TemporaryDirectory(prefix="masld-native-residual-") as folder:
        for j, model in enumerate(path.values()):
            fp = Path(folder)/f"state{j}.npz"
            M.save_state(fp, model)
            restored = M.load_state(fp)
            observed = M.predict(model, x[held], z[held])
            other = M.predict(restored, x[held], z[held])
            assert np.array_equal(observed, other)
            reload_max = max(reload_max, float(np.max(abs(observed-other))))
            assert np.array_equal(-observed, M.predict(restored, -x[held], -z[held]))
            assert np.count_nonzero(M.predict(restored, np.zeros_like(x[held]), np.zeros_like(z[held]))) == 0
    checks.append({"name": "save_reload_zero_and_allele_swap", "pass": True, "reload_max_absolute": reload_max})

    with tempfile.TemporaryDirectory(prefix="masld-native-residual-interface-") as folder:
        directory = Path(folder)
        model = path["joint", M.LAMBDAS[0]]
        model["metadata"].update(length=2048, pooling="variant", native_length=1048576,
            native_readout_bp=501, native_score_order=M.NATIVE_COLUMNS, native_aggregation="DIFF_LOG2_SUM",
            checkpoint="synthetic_checkpoint", fasta="synthetic_reference")
        M.save_state(directory/"model.npz", model)
        payload = {"key": np.asarray([str(i) for i in range(held.sum())]),
            "native_atac_dnase": z[held], "native_score_order": np.asarray(M.NATIVE_COLUMNS),
            "native_length": 1048576, "native_readout_bp": 501, "native_aggregation": "DIFF_LOG2_SUM",
            "checkpoint": "synthetic_checkpoint", "fasta": "synthetic_reference",
            "length": 2048, "pooling": "variant", "allele_order": np.asarray(["REF", "ALT"]),
            "ref": np.zeros_like(x[held]), "alt": x[held]}
        np.savez(directory/"inputs.npz", **payload)
        arguments = SimpleNamespace(model=directory/"model.npz", inputs=directory/"inputs.npz", out=directory/"prediction.tsv")
        M.inference(arguments)
        observed = M.pd.read_csv(arguments.out, sep="\t").predicted_beta_alt.to_numpy()
        np.testing.assert_allclose(observed, M.predict(model, x[held], z[held]), atol=1e-12, rtol=1e-12)
        payload["native_score_order"] = np.asarray(M.NATIVE_COLUMNS[::-1])
        np.savez(directory/"wrong_order.npz", **payload)
        arguments.inputs, arguments.out = directory/"wrong_order.npz", directory/"rejected.tsv"
        try:
            M.inference(arguments)
        except ValueError as exc:
            assert "native channel order" in str(exc)
        else:
            raise AssertionError("Unknown/reversed native channel identity was accepted")
        assert not arguments.out.exists()
    checks.append({"name": "executable_inference_recipe_and_native_channel_identity", "pass": True})

    with tempfile.TemporaryDirectory(prefix="masld-native-residual-labels-") as folder:
        rows = []
        for fold in range(5):
            rows.append({"key": f"{fold+1}:100:A:C", "chr": f"chr{fold+1}", "pos_hg38": 100,
                         "ref": "A", "alt": "C", "peak_id": f"peak{fold}", "peak_start_hg38": 50,
                         "peak_stop_hg38": 150, "heldout_fold": fold, "block_1mb": f"chr{fold+1}:0",
                         "beta_alt": "forbidden_non_numeric_outcome" if fold < 2 else float(fold)})
        source = M.pd.DataFrame(rows)
        target = Path(folder)/"labels.tsv"
        source.to_csv(target, sep="\t", index=False)
        loaded = M.read_labels(target, source[M.IDENTITY], M.TRAIN)
        assert np.isnan(loaded[:2]).all()
        assert np.array_equal(loaded[2:], [2., 3., 4.])
    checks.append({"name": "unrequested_label_strings_never_converted_or_used", "pass": True})

    bad = M.stats(x[train, :2]+rng.normal(size=(train.sum(), 2)), z[train], y[train])
    bad["xx"] = train.sum()*np.array([[1., 2.], [2., 1.]])
    try:
        M.feature_path(bad)
    except ValueError as exc:
        assert "not PSD" in str(exc)
    else:
        raise AssertionError("Indefinite Gram was accepted")
    near = {"xx": np.array([[1., 1.+1e-12], [1.+1e-12, 1.]]),
            "xz": np.zeros((2, 2)), "xy": np.zeros(2), "zz": np.eye(2), "zy": np.zeros(2), "n": 2}
    roundoff = M.feature_path(near)
    assert all(m["metadata"]["eigen_negative_clipped"] == 1 for m in roundoff.values())
    checks.append({"name": "PSD_failure_and_explicit_roundoff_clipping", "pass": True})

    pairs = M.contrast_pairs()
    assert len(pairs) == len(set(pairs)) == 18
    blocks = np.repeat(["a", "b", "c", "d"], 6)
    pred = np.column_stack([y[held], y[held], np.zeros(held.sum())])
    mse, ranks = M.bootstrap(y[held], pred, blocks, 20, M.SEED)
    assert np.array_equal(mse[:, 0], mse[:, 1])
    assert np.array_equal(ranks[:, 0], ranks[:, 1])
    assert np.isnan(ranks[:, 2]).all()
    checks.append({"name": "paired_resampling_complete_family_and_undefined_zero_ranking", "pass": True})
    if out.exists():
        raise FileExistsError(out)
    M.write_json(out, {"pass": True, "checks": checks, "count": len(checks),
                      "seconds": time.monotonic()-tick, "seed": M.SEED,
                      "real_outcomes_read": False, "implementation_sha256": M.sha256(M.__file__)})
    print(json.dumps({"pass": True, "checks": len(checks), "receipt": str(out)}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    run(parser.parse_args().out)

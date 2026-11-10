"""Scientific identity checks for reuse of completed donor-held-out predictions."""
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse


def load(name):
    path = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


helper, core = load("network_external_readout"), load("network_compare")


def test_requested_reusable_methods_cannot_be_missing_or_incomplete():
    frame = pd.DataFrame([dict(outer_fold=fold, gene_id="g", endpoint=endpoint,
        donor_id=str(fold), method=method, squared_error=.1, standardized_squared_error=.2)
        for fold in [0, 1] for endpoint in ["rna", "promoter_accessibility"]
        for method in ["distance", "cicero", "TRIPOD"]])
    assert len(helper.reusable_predictions(frame, "g", ["cicero", "TRIPOD"], 2)) == 8
    with pytest.raises(ValueError, match="Missing or incomplete"):
        helper.reusable_predictions(frame, "g", ["absent"], 2)
    with pytest.raises(ValueError, match="Missing or incomplete"):
        helper.reusable_predictions(frame.iloc[:-1], "g", ["TRIPOD"], 2)


def test_reused_predictions_require_completed_noninvalidated_run(tmp_path):
    (tmp_path / "completion.json").write_text(json.dumps(dict(status="completed")))
    assert helper.completed(tmp_path)["status"] == "completed"
    (tmp_path / "INVALIDATED_input_units.json").write_text("{}")
    with pytest.raises(ValueError, match="Invalidated"):
        helper.completed(tmp_path)


def test_anchor_checks_detect_changed_losses_and_units():
    frame = pd.DataFrame([dict(outer_fold=0, gene_id="g", endpoint="rna" if method == "promoter_only" else "promoter_accessibility",
        method=method, donor_id="331", n_metacells=2, input_representation="binary",
        target_unit="log1p_counts_per_10000", observed=.1, predicted=.2, squared_error=.01,
        standardized_squared_error=.03, ridge_alpha=1.) for method in ["intercept", "promoter_only", "distance", "directnet_adapted"]])
    assert helper.anchor_predictions(frame, frame.copy(), "g")["maximum_absolute_difference"] == 0
    wrong = frame.copy()
    wrong.loc[0, "standardized_squared_error"] += .01
    with pytest.raises(ValueError, match="numerical"):
        helper.anchor_predictions(wrong, frame, "g")
    wrong = frame.copy()
    wrong.loc[0, "target_unit"] = "CPM_no_log"
    with pytest.raises(ValueError, match="units"):
        helper.anchor_predictions(wrong, frame, "g")


def test_cache_reuse_requires_exact_aggregate_counts_and_donor_order(tmp_path):
    matrix = sparse.csr_matrix([[1, 0, 2], [0, 0, 3], [0, 4, 2], [0, 5, 0]])
    metadata = pd.DataFrame(dict(donor_id=["331", "331", "334", "334"]))
    peaks = pd.DataFrame(dict(peak_id=["chr1:1-10", "chr1:20-30", "chr1:40-50"],
                              chrom=["chr1"] * 3, start0=[1, 20, 40], end=[10, 30, 50]))
    folds = pd.DataFrame(dict(outer_fold=[0, 0, 1, 1]))
    for fold in [0, 1]:
        training = np.flatnonzero(folds.outer_fold.to_numpy() != fold)
        nonzero = np.flatnonzero(np.asarray(matrix[training].sum(axis=0)).ravel() > 0)
        stage = tmp_path / f"cicero_fold{fold}_outer"
        stage.mkdir()
        (stage / "links.tsv.method.txt").write_text("covariance=donor_equal_reliability_weighted\ndistance_samples=10\n")
        peaks.iloc[nonzero].to_csv(stage / "peaks.tsv", sep="\t", index=False)
        with core.threadpool_limits(limits=1):
            core.mmwrite(stage / "counts.mtx", matrix[training][:, nonzero].T)
        donors = pd.DataFrame(dict(metacell_id=["metacell_1", "metacell_2"], donor_id=metadata.iloc[training].donor_id.to_numpy()))
        donors.to_csv(stage / "donors.tsv", sep="\t", index=False)
        donors["covariance_weight"] = .5
        donors.to_csv(stage / "links.tsv.donor_weights.tsv", sep="\t", index=False)
    checks = helper.verify_atac_caches(core, tmp_path, matrix, metadata, peaks, folds, 10)
    assert len(checks) == 2 and all(row["exact_count_match"] for row in checks)
    wrong = matrix.copy()
    wrong[0, 0] += 1
    with pytest.raises(ValueError, match="counts differ"):
        helper.verify_atac_caches(core, tmp_path, wrong, metadata, peaks, folds, 10)
    saved = pd.read_csv(tmp_path / "cicero_fold0_outer/donors.tsv", sep="\t", dtype=str)
    saved.loc[0, "donor_id"] = "wrong_donor"
    saved.to_csv(tmp_path / "cicero_fold0_outer/donors.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="donor order"):
        helper.verify_atac_caches(core, tmp_path, matrix, metadata, peaks, folds, 10)


def test_native_raw_input_fingerprint_rejects_changed_contents(tmp_path):
    raw = tmp_path / "counts"
    raw.write_bytes(b"original counts")
    recorded_time = raw.stat().st_mtime_ns
    baseline = tmp_path / "baseline"
    baseline.mkdir()
    (baseline / "run_settings.json").write_text("{}")
    native = tmp_path / "native"
    native.mkdir()
    (native / "input_sha256.json").write_text(json.dumps({str(raw.resolve()): helper.sha256(raw)}))
    assert helper.verify_native_inputs(baseline, [native], [raw])[0]["n_verified_files"] == 1
    raw.write_bytes(b"changed counts")
    os.utime(raw, ns=(recorded_time, recorded_time))
    with pytest.raises(ValueError, match="hash changed"):
        helper.verify_native_inputs(baseline, [native], [raw])

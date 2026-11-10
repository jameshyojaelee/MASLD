"""Reference recovery retains unlabeled backgrounds and separates evidence types."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

SPEC = importlib.util.spec_from_file_location("evaluate_link_reference", Path(__file__).resolve().parents[1] / "scripts/evaluate_link_reference.py")
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def inputs():
    context = dict(reference_type="cis_link", cohort="fixture", lineage="hepatocyte")
    reference = pd.DataFrame([context | dict(source_id="positive", target_id="gene1", reference_id="ref1",
        source_uri="https://example.org/qualified-fixture", qualification_status="qualified")])
    universe = pd.DataFrame([context | dict(source_id=source, target_id="gene1", distance_bin="near",
        observability_bin="high", observable="true") for source in ["positive", "background1", "background2"]] +
        [context | dict(source_id="wrong_distance", target_id="gene1", distance_bin="far", observability_bin="high", observable="true")])
    prediction = pd.DataFrame([context | dict(source_id="positive", target_id="gene1", method="fixture", outer_fold="0", score=1.0),
        context | dict(source_id="wrong_distance", target_id="gene1", method="fixture", outer_fold="0", score=0.5)])
    return reference, universe, prediction


def test_matching_is_reference_recovery_not_false_positive_accuracy():
    reference, universe, prediction = inputs()
    result, audit = module.evaluate_group(prediction, reference, universe, None, np.random.default_rng(42), 1, 20)
    assert result["status"] == "evaluated"
    assert result["n_matched_reference_edges"] == result["n_recovered_reference_edges"] == 1
    assert result["qualified_reference_recovery_fraction"] == 1
    assert result["matched_unlabeled_recovery_fraction"] == 0
    assert result["reference_recovery_enrichment"] == "undefined_zero_background_recovery"
    assert result["nominal_matched_background_tail_probability"] == 1 / 21
    assert audit[0]["prior_overlap_state"] == "unassessed"
    assert "auroc" not in result and "precision" not in result and "false_positive_rate" not in result
    # Wrong-distance selected background must not contaminate matching.
    result2, _ = module.evaluate_group(prediction, reference, universe, None, np.random.default_rng(42), 2, 20)
    assert result2["matched_unlabeled_recovery_fraction"] == 0


def test_prior_edge_and_reference_source_exposures_are_excluded():
    reference, universe, prediction = inputs()
    prior = prediction.iloc[[0]].copy()
    prior["prior_id"], prior["reference_id"] = "prior1", "ref1"
    result, audit = module.evaluate_group(prediction, reference, universe, prior, np.random.default_rng(42), 1, 20)
    assert result["reason"] == "all_reference_edges_exposed_to_prior"
    assert audit[0]["prior_overlap_state"] == "exact_edge_overlap"
    prior["source_id"] = "other_edge"
    _, audit = module.evaluate_group(prediction, reference, universe, prior, np.random.default_rng(42), 1, 20)
    assert audit[0]["prior_overlap_state"] == "reference_source_overlap"


def test_no_matched_background_is_untestable_not_negative():
    reference, universe, prediction = inputs()
    universe.loc[universe.source_id.str.startswith("background"), "observable"] = "false"
    result, _ = module.evaluate_group(prediction, reference, universe, None, np.random.default_rng(42), 1, 20)
    assert result["status"] == "not_run"
    assert result["reason"] == "no_distance_observability_matched_unlabeled_background"
    assert "qualified_reference_recovery_fraction" not in result


def test_native_adapter_is_cis_only_and_rejects_ambiguous_tf(tmp_path):
    path = tmp_path / "GSE244832_hepatocyte/links.tsv"
    path.parent.mkdir()
    pd.DataFrame([dict(peak_id="p", gene_id="g", method="m", outer_fold="0", score=1)]).to_csv(path, sep="\t", index=False)
    adapted = module.adapt_links(path)
    assert adapted.reference_type.eq("cis_link").all()
    adapted["reference_type"], adapted["motif_not_unique_tf"] = "TF_binding", "True"
    adapted.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="Ambiguous"):
        module.adapt_links(path)


def test_no_supplied_reference_writes_explicit_not_run(tmp_path):
    args = SimpleNamespace(output=tmp_path / "out", reference=[], seed=42)
    module.run(args)
    result = json.loads((args.output / "status.json").read_text())
    assert result["status"] == "not_run" and result["reason"] == "no_qualified_supplied_reference"
    assert result["auroc_auprc"] == "not_computed"
    with pytest.raises(FileExistsError):
        module.run(args)


def test_full_adapter_keeps_binding_without_tf_predictions_unassessed(tmp_path):
    reference, universe, prediction = inputs()
    binding = reference.copy()
    binding["reference_type"], binding["source_id"], binding["target_id"] = "TF_binding", "fixture_tf", "fixture_peak"
    reference = pd.concat([reference, binding], ignore_index=True)
    ref_path, uni_path, pred_path = [tmp_path / name for name in ["reference.tsv", "universe.tsv", "predictions.tsv"]]
    reference.to_csv(ref_path, sep="\t", index=False)
    universe.to_csv(uni_path, sep="\t", index=False)
    prediction.to_csv(pred_path, sep="\t", index=False)
    args = SimpleNamespace(output=tmp_path / "out", reference=[ref_path], universe=uni_path,
        prior=None, links=[pred_path], network=[], seed=42, top_k=1, background_repeats=20)
    module.run(args)
    result = module.table(args.output / "reference_recovery.tsv")
    assert result.loc[result.reference_type.eq("cis_link"), "status"].tolist() == ["evaluated"]
    assert result.loc[result.reference_type.eq("TF_binding"), "reason"].tolist() == ["no_compatible_supplied_predictions"]
    assert json.loads((args.output / "status.json").read_text())["n_evaluated_groups"] == 1

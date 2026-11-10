"""Small scientific-invariant fixtures, run inside a compute allocation."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

SPEC = importlib.util.spec_from_file_location("network_compare", Path(__file__).resolve().parents[1] / "scripts/network_compare.py")
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_metacells_preserve_counts_and_never_mix_donors():
    cells = pd.DataFrame({"cohort": ["x"] * 250, "lineage": ["hepatocyte"] * 250,
                          "donor_id": ["a"] * 120 + ["b"] * 130})
    counts = sparse.csr_matrix(np.arange(250 * 3).reshape(250, 3))
    aggregate, meta, groups, eligibility = module.build_metacells(counts, cells)
    assert aggregate.shape == (4, 3)
    assert len(set(np.concatenate(groups))) == 200
    for row, group in enumerate(groups):
        assert cells.iloc[group].donor_id.nunique() == 1
        np.testing.assert_array_equal(aggregate[row].toarray(), counts[group].sum(axis=0))
    again = module.build_metacells(counts, cells)
    for first, second in zip(groups, again[2]):
        np.testing.assert_array_equal(first, second)


def test_heldout_donors_are_disjoint():
    donors = np.repeat([f"d{x}" for x in range(15)], 2)
    folds = module.donor_folds(donors, 5, 42)
    for fold in range(5):
        assert not set(donors[folds == fold]) & set(donors[folds != fold])


def test_promoters_exclude_target_features_on_both_strands():
    peaks = pd.DataFrame({"peak_id": ["p", "l", "r", "far"], "chrom": ["chr1"] * 4,
                          "start0": [9900, 5000, 15000, 900000], "end": [10100, 5500, 15500, 900500]})
    for strand in ["+", "-"]:
        promoter, cis = module.gene_peak_indices(peaks, {"chrom": "chr1", "tss0": 10000, "strand": strand})
        assert promoter.tolist() == [0]
        assert cis.tolist() == [1, 2]
        assert not set(promoter) & set(cis)


def test_donor_equal_weighting_ignores_metacell_imbalance():
    donors = np.array(["a"] * 10 + ["b"])
    weights = module.donor_weights(donors)
    assert np.isclose(weights[:10].sum(), weights[-1])
    result = module.weighted_ridge(np.empty((11, 0)), np.r_[np.zeros(10), 2], donors, 1, np.empty((2, 0)))
    np.testing.assert_allclose(result, [1, 1])


def test_signflip_null_is_not_bootstrap_tail_fraction():
    result = module.paired_inference(np.ones(6), seed=42, bootstrap_reps=100)
    assert result["pvalue"] == 1 / 64
    assert result["ci_low"] == result["ci_high"] == 1
    null = module.paired_inference(np.zeros(6), seed=42, bootstrap_reps=100)
    assert null["pvalue"] == 1


def test_holm_preserves_family_and_monotonicity():
    np.testing.assert_allclose(module.holm([0.01, 0.04, 1]), [0.03, 0.08, 1])


def test_missing_methods_stay_in_prespecified_family():
    frame = pd.DataFrame([dict(outer_fold=0, gene_id="g", endpoint="promoter_accessibility", method="intercept", donor_id="d", metacell_id="m", standardized_squared_error=1)])
    comparison = module.compare_predictions(frame, ["distance", "directnet_adapted", "cicero"], 42, 20, ["promoter_accessibility", "rna"])
    assert len(comparison) == 12
    assert ((comparison.method == "directnet_adapted") & (comparison.baseline == "cicero")).sum() == 2
    assert comparison.pvalue.eq(1).all()


def test_source_dependent_lineage_propagates():
    cells = pd.DataFrame({"cohort": ["x"] * 100, "lineage": ["stellate"] * 100,
                          "donor_id": ["a"] * 100, "lineage_state": ["source_dependent"] * 100,
                          "source_lineage": ["mesenchymal"] * 100})
    _, meta, _, _ = module.build_metacells(sparse.csr_matrix(np.ones((100, 2))), cells)
    assert meta.lineage_state.eq("source_dependent").all()
    assert meta.source_lineage.eq("mesenchymal").all()


def test_small_nested_distance_run(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "inventory", lambda: {"xgboost": None})
    rng = np.random.default_rng(42)
    counts = sparse.csr_matrix(rng.poisson(5, size=(20, 4)))
    meta = pd.DataFrame({"donor_id": np.repeat([f"d{i}" for i in range(10)], 2),
                         "metacell_id": [f"m{i}" for i in range(20)]})
    peaks = pd.DataFrame({"peak_id": ["p", "d1", "d2", "far"], "chrom": ["chr1"] * 4,
                          "start0": [9900, 5000, 15000, 900000], "end": [10100, 5500, 15500, 900500]})
    genes = pd.DataFrame([dict(gene_id="g", chrom="chr1", tss0=10000, strand="+")])
    out = tmp_path / "result"
    comparisons = module.fit_group(counts, meta, peaks, genes, out,
                                   {"cicero_rscript": "unavailable", "bootstrap_reps": 20}, 42)
    pred = pd.read_csv(out / "predictions.tsv", sep="\t")
    assert set(pred.method) == {"intercept", "distance"}
    assert pred.groupby("donor_id").outer_fold.nunique().eq(1).all()
    assert comparisons.n_common_predictions.max() == 10
    assert len(comparisons) == 6
    fold_map = {f"d{i}": 0 if i < 5 else 4 for i in range(10)}
    module.fit_group(counts, meta, peaks, genes, tmp_path / "fixed_folds",
                     {"cicero_rscript": "unavailable", "bootstrap_reps": 20,
                      "outer_fold_map": fold_map}, 42)
    fixed = pd.read_csv(tmp_path / "fixed_folds" / "predictions.tsv", sep="\t")
    assert set(fixed.outer_fold) == {0, 4}
    assert all(row.outer_fold == fold_map[row.donor_id] for row in fixed.itertuples())


def test_donor_loss_is_mean_squared_error_not_squared_mean_error():
    rows = []
    meta = pd.DataFrame({"donor_id": ["a", "a", "b"], "metacell_id": ["m1", "m2", "m3"]})
    module.save_predictions(rows, 0, "g", "rna", "distance", np.array([0, 1, 2]),
                            meta, np.array([0., 2., 3.]), np.array([1., 1., 4.]), 2., 1.)
    first = rows[0]
    assert first["observed"] == first["predicted"] == 1
    assert first["squared_error"] == 1
    assert first["standardized_squared_error"] == .5
    assert first["n_metacells"] == 2


def test_donor_first_reduction_preserves_paired_inference():
    rng = np.random.default_rng(42)
    donors = np.repeat([f"d{i}" for i in range(6)], 3)
    meta = pd.DataFrame({"donor_id": donors})
    y, candidate, baseline = rng.normal(size=(3, 18))
    rows = []
    for method, pred in [("distance", candidate), ("intercept", baseline)]:
        module.save_predictions(rows, 0, "g", "promoter_accessibility", method,
                                np.arange(18), meta, y, pred, 2., 1.)
    reduced = pd.DataFrame(rows)
    result = module.compare_predictions(reduced, ["distance"], 42, 100).iloc[0]
    original = pd.DataFrame({"donor": donors, "difference": ((baseline-y)**2-(candidate-y)**2)/2}).groupby("donor").difference.mean()
    expected = module.paired_inference(original, module.stable_seed(42, ("promoter_accessibility", "distance", "intercept")), 100)
    for name in ["improvement", "ci_low", "ci_high", "pvalue"]:
        assert np.isclose(result[name], expected[name])


def test_explicit_prediction_disk_row_limit_is_respected(tmp_path):
    counts = sparse.csr_matrix(np.ones((20, 4)))
    meta = pd.DataFrame({"donor_id": np.repeat([f"d{i}" for i in range(10)], 2),
                         "metacell_id": [f"m{i}" for i in range(20)]})
    genes = pd.DataFrame([dict(gene_id="g", chrom="chr1", tss0=10000, strand="+")])
    with pytest.raises(MemoryError, match="max_prediction_rows"):
        module.fit_group(counts, meta, pd.DataFrame(), genes, tmp_path / "large",
                         {"max_genes": 0, "max_prediction_rows": 1}, 42)


def test_source_pilot_and_zero_gene_config_cannot_become_research(tmp_path, monkeypatch):
    import json
    import sys
    source = tmp_path / "input"
    source.mkdir()
    counts = sparse.csr_matrix(np.ones((1000, 2)))
    sparse.save_npz(source / "cell_counts.npz", counts)
    pd.DataFrame({"cohort": ["fixture"] * 1000, "lineage": ["hepatocyte"] * 1000,
                  "donor_id": np.repeat([f"d{i}" for i in range(10)], 100),
                  "cell_id": [f"c{i}" for i in range(1000)]}).to_csv(source / "cells.tsv", sep="\t", index=False)
    pd.DataFrame({"peak_id": ["chr1:9900-10100", "chr1:5000-5500"], "chrom": ["chr1", "chr1"],
                  "start0": [9900, 5000], "end": [10100, 5500]}).to_csv(source / "peaks.tsv", sep="\t", index=False)
    (source / "diagnostics.json").write_text(json.dumps({"pilot": True}))
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"network": {"max_genes": 0}}))
    genes = tmp_path / "genes.tsv"
    pd.DataFrame([dict(gene_id="g", chrom="chr1", tss0=10000, strand="+")]).to_csv(genes, sep="\t", index=False)
    captured = []
    def fake_fit(*args, **kwargs):
        captured.append(args[5])
        return pd.DataFrame()
    monkeypatch.setattr(module, "require_compute", lambda: None)
    monkeypatch.setattr(module, "fit_group", fake_fit)
    output = tmp_path / "output"
    monkeypatch.setattr(sys, "argv", ["network_compare", "--config", str(config), "--input-root", str(source), "--genes", str(genes), "--output", str(output)])
    module.main()
    assert captured[0]["max_genes"] == 25
    settings = json.loads((output / "run_settings.json").read_text())
    assert settings["pilot"] is True and settings["source_pilot"] is True
    assert json.loads((output / "completion.json").read_text())["interpretation"] == "technical_pilot_no_research_calls"


def test_streaming_comparisons_equal_full_donor_table():
    rng = np.random.default_rng(7)
    rows = []
    accumulator = module.ComparisonAccumulator(["promoter_accessibility", "rna"], ["distance", "directnet_adapted", "cicero"])
    for fold in range(2):
        for gene in ["g1", "g2", "g3"]:
            group = []
            for endpoint, baseline in [("rna", "promoter_only"), ("promoter_accessibility", "intercept")]:
                for donor in [f"d{2*fold}", f"d{2*fold+1}"]:
                    for method in [baseline, "distance", "directnet_adapted", "cicero"]:
                        if gene == "g2" and method == "cicero":
                            continue
                        group.append(dict(outer_fold=fold, gene_id=gene, endpoint=endpoint, method=method,
                                          donor_id=donor, n_metacells=2, standardized_squared_error=float(rng.uniform())))
            accumulator.append(group)
            rows.extend(group)
    actual, donors = accumulator.finalize(42, 100)
    expected = module.compare_predictions(pd.DataFrame(rows), ["distance", "directnet_adapted", "cicero"], 42, 100, ["promoter_accessibility", "rna"])
    pd.testing.assert_frame_equal(actual, expected)
    assert len(donors) == 48


def test_sparse_cicero_query_preserves_reverse_negative_and_zero_scores(tmp_path):
    peaks = pd.DataFrame(dict(peak_id=["p0", "p1", "p2", "p3"], start0=[0, 100, 200, 300], end=[50, 150, 250, 350]))
    edges = pd.DataFrame(dict(peak1=["p0", "p2", "p0", "p1"], peak2=["p1", "p0", "p2", "p3"], coaccess=[-.2, -.1, -.4, 0.]))
    path = tmp_path / "links.tsv"
    edges.to_csv(path, sep="\t", index=False)
    graph = module.CiceroGraph.from_tsv(path, peaks)
    assert graph.promoter_scores([0]) == {1: -.2, 2: -.1}
    assert graph.promoter_scores([1])[3] == 0
    x = sparse.csr_matrix([[1, 2, 3, 4], [2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]])
    arguments = ("cicero", np.array([1, 2, 3]), np.array([0]), peaks, {"tss0": 0}, x, np.arange(4), np.array(["a", "b", "c", "d"]), 42, 10)
    original = module.selected_indices(*arguments, cicero=edges)
    optimized = module.selected_indices(*arguments, cicero=graph)
    np.testing.assert_array_equal(original[0], optimized[0])
    np.testing.assert_array_equal(original[1], optimized[1])
    pd.concat([edges, edges.iloc[:1]]).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="Duplicate"):
        module.CiceroGraph.from_tsv(path, peaks)


def test_indexed_geometry_matches_exact_interval_masks():
    rng = np.random.default_rng(42)
    starts = rng.integers(0, 1000000, size=1000)
    peaks = pd.DataFrame(dict(chrom=rng.choice(["chr1", "chr2"], 1000), start0=starts,
                              end=starts + rng.integers(1, 10000, 1000)))
    geometry = module.PeakGeometry(peaks)
    for chrom in ["chr1", "chr2", "chrX"]:
        for strand in ["+", "-"]:
            for tss in [0, 500000, 1000000]:
                gene = dict(chrom=chrom, strand=strand, tss0=tss)
                actual = geometry.query(gene)
                expected = module.gene_peak_indices(peaks, gene)
                for a, b in zip(actual, expected):
                    np.testing.assert_array_equal(a, b)


def test_condition_roster_excludes_masl_without_filtering_paired_context():
    cells = pd.DataFrame(dict(cohort=["GSE244832"]*3+["GSE296875"], lineage=["hepatocyte"]*4,
                              donor_id=["a", "b", "c", "d"], condition=[" normal ", "MASL", "MASH", "not_qualified"]))
    settings = {"allowed_conditions_by_cohort": {"GSE244832": ["NORMAL", "MASH"]}}
    keep, excluded = module.condition_eligibility(cells, settings)
    np.testing.assert_array_equal(keep, [True, False, True, True])
    assert excluded.donor_id.tolist() == ["b"]
    assert excluded.n_excluded_cells.tolist() == [1]
    pilot, exclusions = module.condition_eligibility(cells, settings, technical_pilot=True)
    assert pilot.all() and exclusions.empty


def test_external_links_require_exact_inner_and_outer_training_donors(tmp_path, monkeypatch):
    import json
    monkeypatch.setattr(module, "inventory", lambda: {"xgboost": None})
    rng = np.random.default_rng(42)
    counts = sparse.csr_matrix(rng.poisson(5, size=(20, 4)))
    donors = np.repeat([f"d{i}" for i in range(10)], 2)
    meta = pd.DataFrame(dict(donor_id=donors))
    peaks = pd.DataFrame(dict(peak_id=["p", "d1", "d2", "far"], chrom=["chr1"]*4,
                              start0=[9900, 5000, 15000, 900000], end=[10100, 5500, 15500, 900500]))
    genes = pd.DataFrame([dict(gene_id="g", chrom="chr1", tss0=10000, strand="+")])
    mapping = {f"d{i}": 0 if i < 5 else 4 for i in range(10)}
    outer = np.array([mapping[x] for x in donors])
    paths = {}
    (tmp_path / "completion.json").write_text(json.dumps(dict(status="completed")))
    for fold in [0, 4]:
        train = np.flatnonzero(outer != fold)
        inner = module.donor_folds(donors[train], 3, module.stable_seed(42, ("inner", fold)))
        for tag, fit in [("outer", train)] + [(f"inner{k}", train[inner != k]) for k in sorted(set(inner))]:
            path = tmp_path / f"links_{fold}_{tag}.tsv"
            pd.DataFrame(dict(gene_id=["g", "g"], peak_id=["d1", "d2"], score=[.8, .4])).to_csv(path, sep="\t", index=False)
            path.with_suffix(".training_donors.json").write_text(json.dumps(dict(donor_ids=sorted(set(donors[fit])))))
            paths[f"{fold}:{tag}"] = str(path)
    # A wrong outer donor manifest must suppress this method for that entire fold.
    Path(paths["0:outer"]).with_suffix(".training_donors.json").write_text(json.dumps(dict(donor_ids=["held_out_donor"])))
    output = tmp_path / "external_run"
    module.fit_group(counts, meta, peaks, genes, output,
                     dict(cicero_rscript="unavailable", bootstrap_reps=20, outer_fold_map=mapping,
                          external_stage_links={"TRIPOD_fixture": paths},
                          external_stage_roots={"TRIPOD_fixture": str(tmp_path)}), 42)
    predictions = pd.read_csv(output / "predictions.tsv", sep="\t")
    assert set(predictions[predictions.method == "TRIPOD_fixture"].outer_fold) == {4}
    failure = pd.read_csv(output / "failures.tsv", sep="\t")
    assert failure.reason.str.contains("donor set").any()
    exploratory = pd.read_csv(output / "exploratory_comparisons.tsv", sep="\t")
    assert set(exploratory.method) == {"TRIPOD_fixture"}
    assert not exploratory.predictive_improvement.any()
    assert len(pd.read_csv(output / "comparisons.tsv", sep="\t")) == 6


def test_external_invalidated_or_uncompleted_stages_are_never_consumed(tmp_path, monkeypatch):
    import json
    monkeypatch.setattr(module, "inventory", lambda: {"xgboost": None})
    counts = sparse.csr_matrix(np.random.default_rng(1).poisson(5, (20, 3)))
    meta = pd.DataFrame(dict(donor_id=np.repeat([f"d{i}" for i in range(10)], 2)))
    peaks = pd.DataFrame(dict(peak_id=["p", "a", "b"], chrom=["chr1"]*3,
                              start0=[9900, 5000, 15000], end=[10100, 5500, 15500]))
    genes = pd.DataFrame([dict(gene_id="g", chrom="chr1", tss0=10000, strand="+")])
    source = tmp_path / "invalid_native"
    source.mkdir()
    (source / "completion.json").write_text(json.dumps(dict(status="completed")))
    (source / "INVALIDATED_TEST.json").write_text("{}")
    pd.DataFrame(dict(gene_id=["g"], peak_id=["a"], score=[.5])).to_csv(source / "links.tsv", sep="\t", index=False)
    settings = dict(cicero_rscript="unavailable", bootstrap_reps=20,
                    external_stage_roots={"native_fixture": str(source)},
                    external_stage_links={"native_fixture": {f"{fold}:{tag}": str(source / "links.tsv") for fold in range(5) for tag in ["outer", "inner0", "inner1", "inner2"]}})
    output = tmp_path / "rejected"
    module.fit_group(counts, meta, peaks, genes, output, settings, 42)
    failure = pd.read_csv(output / "failures.tsv", sep="\t")
    assert failure.reason.str.contains("invalidated").all()
    assert "native_fixture" not in set(pd.read_csv(output / "predictions.tsv", sep="\t").method)
    (source / "INVALIDATED_TEST.json").unlink()
    (source / "completion.json").write_text(json.dumps(dict(status="running")))
    output2 = tmp_path / "not_completed"
    module.fit_group(counts, meta, peaks, genes, output2, settings, 42)
    failure = pd.read_csv(output2 / "failures.tsv", sep="\t")
    assert failure.reason.str.contains("not completed").all()


def test_peak_identifier_must_match_coordinate_columns(tmp_path):
    sparse.save_npz(tmp_path / "cell_counts.npz", sparse.csr_matrix([[1, 2]]))
    pd.DataFrame(dict(cohort=["x"], lineage=["hepatocyte"], donor_id=["d"], cell_id=["c"])).to_csv(tmp_path / "cells.tsv", sep="\t", index=False)
    pd.DataFrame(dict(peak_id=["chr1:1-10", "chr1:20-30"], chrom=["chr1", "chr1"],
                      start0=[1, 2000], end=[10, 3000])).to_csv(tmp_path / "peaks.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="canonical"):
        module.load_input(tmp_path)


def test_explicit_native_granges_ids_map_one_to_one_to_bed(tmp_path):
    import json
    sparse.save_npz(tmp_path / "cell_counts.npz", sparse.csr_matrix([[1, 2]]))
    pd.DataFrame(dict(cohort=["GSE296875"], lineage=["hepatocyte"], donor_id=["d"], cell_id=["c"])).to_csv(tmp_path / "cells.tsv", sep="\t", index=False)
    peaks = pd.DataFrame(dict(peak_id=["chr1-2-10", "chr1-21-30"], chrom=["chr1", "chr1"], start0=[1, 20], end=[10, 30]))
    peaks.to_csv(tmp_path / "peaks.tsv", sep="\t", index=False)
    (tmp_path / "paired_export.json").write_text(json.dumps(dict(peak_id_coordinate_system="GRanges_1based_closed")))
    _, _, mapped = module.load_input(tmp_path)
    assert mapped.canonical_peak_id.tolist() == ["chr1:1-10", "chr1:20-30"]
    assert mapped.peak_id.tolist() == peaks.peak_id.tolist()
    peaks.loc[1, "start0"] = 21
    peaks.to_csv(tmp_path / "peaks.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="mapping disagrees"):
        module.load_input(tmp_path)


def test_fold_seed_is_identical_for_python_and_numpy_scalar_ids():
    assert module.stable_seed(42, (0, "outer")) == module.stable_seed(42, (np.int64(0), "outer"))
    assert module.stable_seed(42, (0, 2, "g")) == module.stable_seed(42, (np.int64(0), np.int64(2), "g"))
    assert module.stable_seed(42, ("inner", 4)) == module.stable_seed(42, ("inner", np.int64(4)))


def test_learned_links_never_pad_zero_scores_and_keep_negative_coupling():
    peaks = pd.DataFrame(dict(peak_id=["p", "a", "b", "c"], start0=[100, 500, 1000, 1500], end=[150, 550, 1050, 1550]))
    x = sparse.csr_matrix([[1, 2, 3, 4], [2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]])
    arguments = (np.array([1, 2, 3]), np.array([0]), peaks, dict(gene_id="g", tss0=125), x, np.arange(4), np.array(["a", "b", "c", "d"]), 42, 10)
    edges = pd.DataFrame(dict(peak1=["p"]*3, peak2=["a", "b", "c"], coaccess=[0., -.2, .1]))
    selected, scores = module.selected_indices("cicero", *arguments, cicero=edges)
    np.testing.assert_array_equal(selected, [3, 2])
    np.testing.assert_array_equal(scores, [.1, -.2])
    edges.coaccess = 0.
    assert len(module.selected_indices("cicero", *arguments, cicero=edges)[0]) == 0
    external = pd.DataFrame(dict(gene_id=["g"]*3, peak_id=["a", "b", "c"], score=[0., 0., 0.]))
    assert len(module.selected_indices("TRIPOD_fixture", *arguments, external=external)[0]) == 0
    assert len(module.selected_indices("distance", *arguments)[0]) == 3


def test_feature_observability_gives_small_donors_equal_weight():
    donors = np.array([f"d{i}" for i in range(12) for _ in range(20)] + ["small", "small"])
    counts = np.ones((len(donors), 2))
    counts[:, 1] = 0
    counts[-1, 1] = 1
    peaks = pd.DataFrame(dict(peak_id=["p", "rare"], start0=[100, 500], end=[150, 550]))
    chosen, _ = module.selected_indices("distance", np.array([1]), np.array([0]), peaks,
        dict(gene_id="g", tss0=125), sparse.csr_matrix(counts), np.arange(len(donors)), donors, 42, 10)
    # One of a small donor's two aggregates is observable. Extra aggregates in
    # other donors cannot erase its contribution to the 1% donor-mean filter.
    np.testing.assert_array_equal(chosen, [1])

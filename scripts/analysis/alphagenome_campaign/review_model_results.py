#!/usr/bin/env python3
"""Independently reconstruct deposited reporter and adaptation comparisons.

No producer metric functions are imported. This checks molecular labels, split
identities, experimental pairing, and the numerical comparison; it does not
turn a development comparison into protected confirmation.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[3]
ARMS = ("direct_interaction", "separate_conditions", "shared_condition_residual")


def review_reporter(directory, out):
    metadata = json.loads((directory/"analysis.json").read_text())
    source = pd.read_csv(metadata["source"], sep="\t")
    authority = pd.read_csv(metadata["split_authority"], sep="\t",
        usecols=["element_id", "outer_fold", "long_range_block_id"]).drop_duplicates()
    assert not authority.element_id.duplicated().any()
    authority = authority.set_index("element_id")
    metrics = pd.read_csv(directory/"performance.tsv", sep="\t").set_index(["cell_line", "arm"])
    intervals = pd.read_csv(directory/"paired_intervals.tsv", sep="\t").set_index(["cell_line", "arm"])
    selections = pd.read_csv(directory/"inner_selection.tsv", sep="\t")
    rows = []
    for cell in ("LX2", "HepG2"):
        table = pd.read_csv(directory/f"{cell}_oof.tsv", sep="\t").set_index("element_id")
        assert not table.index.duplicated().any()
        auth = authority.loc[table.index]
        folds = auth.outer_fold.astype(str).str.replace("fold-", "", regex=False).astype(int)
        np.testing.assert_array_equal(table.fold, folds)
        np.testing.assert_array_equal(table.long_range_block_id, auth.long_range_block_id)
        assert table.groupby("long_range_block_id").fold.nunique().max() == 1
        native = source.loc[source.cell_line.eq(cell) & source.eligible.astype(str).str.lower().eq("true")]
        assert not native.duplicated(["element_id", "experimental_replicate"]).any()
        control = native.pivot(index="element_id", columns="experimental_replicate", values="effect_control").loc[table.index]
        treated = native.pivot(index="element_id", columns="experimental_replicate", values="effect_treated").loc[table.index]
        deposited = native.pivot(index="element_id", columns="experimental_replicate", values="effect").loc[table.index]
        assert control.columns.equals(treated.columns) and control.columns.equals(deposited.columns)
        assert control.shape[1] == 4
        paired = treated.to_numpy()-control.to_numpy()
        assert np.isfinite(paired).all()
        np.testing.assert_allclose(paired, deposited, atol=1e-10, rtol=0)
        observed = paired.mean(axis=1)
        np.testing.assert_allclose(table.observed_interaction, observed, atol=1e-10, rtol=0)
        np.testing.assert_allclose(table.paired_mean_variance, paired.var(axis=1, ddof=1)/4, atol=1e-10, rtol=0)
        c, t = control.to_numpy(), treated.to_numpy()
        covariance = np.sum((c-c.mean(axis=1)[:, None])*(t-t.mean(axis=1)[:, None]), axis=1)/3
        np.testing.assert_allclose(table.control_treated_sample_covariance, covariance, atol=1e-10, rtol=0)
        np.testing.assert_array_equal(table.zero_interaction, 0)
        for selection in selections.loc[selections.cell_line.eq(cell)].itertuples():
            assert selection.inner_fold != selection.held_fold
            assert selection.test_n == sum(table.fold == selection.held_fold)
            assert selection.train_n == sum(table.fold != selection.held_fold)
        predicted = table[list(ARMS)].to_numpy()
        zero_mse = np.mean(observed**2)
        mse = np.mean((predicted-observed[:, None])**2, axis=0)
        np.testing.assert_allclose(np.sqrt(zero_mse), metrics.loc[(cell, "zero_interaction"), "RMSE_log2_activity_interaction"], atol=1e-12, rtol=0)
        # A shared draw of the four replicate indices preserves covariance
        # among every construct and condition. Resample original locus blocks.
        rng = np.random.default_rng(20260915)
        blocks = table.long_range_block_id.to_numpy()
        keys = np.asarray(sorted(set(blocks)))
        block_rows = {key: np.flatnonzero(blocks == key) for key in keys}
        resamples = int(intervals.loc[(cell, ARMS[0]), "resamples"])
        improvements = np.empty((resamples, len(ARMS)))
        for draw in range(resamples):
            index = np.concatenate([block_rows[key] for key in rng.choice(keys, len(keys), replace=True)])
            repindex = rng.integers(0, 4, 4)
            sampled = paired[index][:, repindex].mean(axis=1)
            # Algebraically independent expression for reduction from zero.
            improvements[draw] = np.mean(2*sampled[:, None]*predicted[index]-predicted[index]**2, axis=0)
        ci = np.quantile(improvements, [0.025, 0.975], axis=0)
        for j, arm in enumerate(ARMS):
            rank = float(spearmanr(observed, predicted[:, j]).statistic)
            np.testing.assert_allclose([np.sqrt(mse[j]), rank], metrics.loc[(cell, arm), ["RMSE_log2_activity_interaction", "signed_spearman"]].to_numpy(dtype=float), atol=1e-12, rtol=0)
            original = intervals.loc[(cell, arm), ["MSE_improvement_vs_zero_low95", "MSE_improvement_vs_zero_high95"]].to_numpy(dtype=float)
            np.testing.assert_allclose(ci[:, j], original, atol=1e-12, rtol=0)
            rows.append(dict(cell_line=cell, arm=arm, constructs=len(table), loci=len(keys),
                replicates=4, rmse=np.sqrt(mse[j]), zero_rmse=np.sqrt(zero_mse), signed_spearman=rank,
                mse_reduction=zero_mse-mse[j], ci95_low=ci[0, j], ci95_high=ci[1, j]))
    pd.DataFrame(rows).to_csv(out/"reporter_reconstructed.tsv", sep="\t", index=False)
    return dict(contrasts=len(rows), source_pairs_reconstructed=True, covariance_reconstructed=True,
        operative_folds_match=True, source_native_units=metadata["measurement_unit"],
        interval_reconstruction="joint_experimental_replicate_and_locus_2000_seed20260915",
        population_counts=[dict(cell_line=cell, constructs=next(r["constructs"] for r in rows if r["cell_line"] == cell),
            loci=next(r["loci"] for r in rows if r["cell_line"] == cell)) for cell in ("LX2", "HepG2")],
        population_rule=metadata.get("representation_population_rule", "historical_1033_subset"),
        limitation="fixed_OOF_fits; fitting_and_selection_uncertainty_not_included")


def review_sweep(directory, out, require_complete):
    records = json.loads((directory/"runs.json").read_text())
    labels = pd.read_csv(ROOT/"GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz", sep="\t").set_index("lead_variant_id")
    assert not labels.index.duplicated().any()
    completed = [r for r in records if r.get("status") == "completed_fixed_exposure" and r["name"] != "probe_microbatch4"]
    if require_complete:
        assert len(completed) == 12, "The twelve-arm comparison is incomplete"
    expected_keys, rows = None, []
    for record in completed:
        path = Path(record["output"])
        feasibility = json.loads((path/"feasibility.json").read_text())
        split = json.loads((path/"split.json").read_text())
        config = json.loads((path/"config.json").read_text())
        checks = json.loads((path/"pretraining_checks.json").read_text())
        drift = json.loads((path/"native_function_drift.json").read_text())
        assert record["exit_code"] == 0 and record["steps_completed"] == record["steps_requested"]
        assert feasibility["completed_steps"] == record["steps_requested"]
        assert all(feasibility[k] for k in ("frozen_parameters_unchanged", "running_statistics_unchanged", "save_reload_agreement"))
        assert all(checks[k] for k in ("identical_zero", "swapped_sign", "zero_adapter_or_initial_partial_exact"))
        assert not split["held_outcomes_evaluated"]
        table = pd.read_csv(path/"validation_predictions.tsv", sep="\t").set_index("variant_id").sort_index()
        assert not table.index.duplicated().any()
        assert len(table) == feasibility["validation_examples"]
        assert np.isfinite(table[["observed_beta", "predicted_beta"]]).all().all()
        known = labels.loc[table.index]
        assert (known.heldout_fold == split["validation_fold"]).all()
        assert (table.validation_fold == split["validation_fold"]).all()
        assert not (known.heldout_fold == split["held_fold"]).any()
        # Trainer explicitly stores float32 supervision.
        np.testing.assert_array_equal(table.observed_beta.to_numpy(dtype=np.float32), known.beta_alt.to_numpy(dtype=np.float32))
        if expected_keys is None:
            expected_keys = table.index
        else:
            assert expected_keys.equals(table.index), "Arms have different eligible validation rows"
        for assay in ("atac", "dnase"):
            assert assay in drift and drift[assay]["tracks"] > 0
            assert np.isfinite([drift[assay]["rms_log2_track_change"], drift[assay]["max_abs_log2_track_change"]]).all()
            if config["mode"] == "frozen":
                assert drift[assay]["max_abs_log2_track_change"] == 0
        y, p = table.observed_beta.to_numpy(), table.predicted_beta.to_numpy()
        rows.append(dict(arm=record["name"], n=len(table), steps=record["steps_completed"],
            rmse=float(np.sqrt(np.mean((p-y)**2))), spearman=float(spearmanr(y, p).statistic),
            signed_concordance=float(np.mean(np.sign(y) == np.sign(p))),
            seed=config["seed"], held_fold=split["held_fold"], validation_fold=split["validation_fold"],
            native_drift_atac=drift["atac"]["rms_log2_track_change"],
            native_drift_dnase=drift["dnase"]["rms_log2_track_change"]))
    pd.DataFrame(rows).to_csv(out/"sweep_reconstructed.tsv", sep="\t", index=False)
    return dict(completed_arms=len(completed), intended_arms=12, common_eligible_rows=len(expected_keys) if expected_keys is not None else 0,
        held_outcomes_read=False, limitation="single_split_single_seed_feasibility; no_protected_test_or_formal_finalist",
        unfinished=[r["name"] for r in records if r.get("status") != "completed_fixed_exposure"])


def review_figures(directory, reporter, donor):
    import hashlib
    import subprocess
    sources = json.loads((directory/"sources.json").read_text())
    for path, expected in sources.items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected
    plotted = pd.read_csv(directory/"plotted_values.tsv", sep="\t")
    performance = pd.read_csv(reporter/"performance.tsv", sep="\t").set_index(["cell_line", "arm"])
    population = pd.read_csv(reporter/"population_coverage.tsv", sep="\t").set_index("cell_line")
    reporter_ci = pd.read_csv(reporter/"paired_intervals.tsv", sep="\t").set_index(["cell_line", "arm"])
    donor_ci = pd.read_csv(donor/"paired_uncertainty.tsv", sep="\t")
    donor_ci = donor_ci.loc[donor_ci.endpoint == "donor_difference_mse_log2cpm"].set_index(["evaluation", "model"])
    reporter_names = dict(zip(("Direct interaction", "Separate condition heads", "Shared + condition residual"), ARMS))
    objective_names = {"Profile objective": "all_rna:shared_profile", "Training-mean residual": "all_rna:shared_training_mean_residual",
        "Donor-pair objective": "all_rna:shared_pairwise", "Additive sequence + RNA": "all_rna:additive", "Region-specific RNA": "all_rna:rna_only_region_specific"}
    exclusion_names = dict(zip(("All eligible RNA", "Overlapping genes excluded", "RNA within ±1 Mb excluded", "Target chromosome excluded",
        "Matched random: overlap", "Matched random: ±1 Mb", "Matched random: chromosome", "Technical inputs only", "RNA context shuffled"),
        ("all_rna", "overlapping_genes_excluded", "within_1mb_excluded", "chromosome_excluded", "random_overlap", "random_1mb", "random_chromosome", "technical_only", "context_shuffled")))
    for row in plotted.itertuples():
        if row.context in ("LX2", "HepG2"):
            arm = reporter_names[row.row]
            metric = performance.loc[(row.context, arm)]
            zero = performance.loc[(row.context, "zero_interaction"), "RMSE_log2_activity_interaction"]**2
            ci = reporter_ci.loc[(row.context, arm)]
            expected = [zero-metric.RMSE_log2_activity_interaction**2, ci.MSE_improvement_vs_zero_low95, ci.MSE_improvement_vs_zero_high95]
            assert row.n_constructs == metric.constructs and row.n_blocks == metric.locus_blocks and row.experimental_replicates == 4
        else:
            model = objective_names[row.row] if row.panel.endswith("objectives") else exclusion_names[row.row]+":shared_pairwise"
            ci = donor_ci.loc[(row.context, model)]
            expected = [-ci.mse_difference, -ci.ci95_high, -ci.ci95_low]
            assert row.n_donors == 99 and row.n_regions == (128 if row.context == "held_donors_trained_regions" else 64)
        np.testing.assert_allclose([row.point, row.low95, row.high95], expected, atol=1e-12, rtol=0)
    pdfs = sorted(directory.glob("*.pdf"))
    assert len(pdfs) == 4 and len(plotted) == 33
    pdf_text = {}
    for path in pdfs:
        content = subprocess.check_output(["pdftotext", "-layout", str(path), "-"], text=True)
        assert "MSE reduction" in content
        if "reporter" in path.stem:
            cell = "LX2" if "LX2" in path.stem else "HepG2"
            measured = performance.loc[(cell, "zero_interaction")]
            assert f"{int(measured.constructs):,} of {int(population.loc[cell, 'full_source_constructs']):,}" in content
            assert f"{int(measured.locus_blocks)} locus blocks" in content and "4 paired experimental replicates" in content
        else:
            assert "99 participants" in content
        pdf_text[path.name] = content
    return dict(pdfs=len(pdfs), plotted_rows=len(plotted), point_and_interval_match=True,
        source_hashes_match=True, pdf_population_and_unit_text_checked=True, pdf_text=pdf_text,
        visual_inspection="separate_rendered_panel_review_required")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lx2", type=Path)
    parser.add_argument("--sweep", type=Path)
    parser.add_argument("--require-sweep-complete", action="store_true")
    parser.add_argument("--figures", type=Path)
    parser.add_argument("--donor", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run scientific reconstruction on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    result = dict(status="pass", producer_metric_imports=False, protected_outcomes_read=False)
    if args.lx2:
        result["reporter"] = review_reporter(args.lx2, args.out)
    if args.sweep:
        result["sweep"] = review_sweep(args.sweep, args.out, args.require_sweep_complete)
    if args.figures:
        if not args.lx2 or not args.donor:
            raise ValueError("Figure checks require reporter and donor source directories")
        result["figures"] = review_figures(args.figures, args.lx2, args.donor)
    (args.out/"checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

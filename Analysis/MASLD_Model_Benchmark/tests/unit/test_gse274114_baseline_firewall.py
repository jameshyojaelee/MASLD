#!/usr/bin/env python3
"""Synthetic, outcome-firewalled tests for GSE274114 bulk baselines."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import tomllib
import unittest
from unittest.mock import patch
import warnings

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.contracts import PredictionBundle
from masld_bench.hashing import sha256_file
from scripts import gse274114_baseline_firewall as firewall


ROOT = Path(__file__).resolve().parents[2]
GROUP_COUNTS = {"CTRL": 9, "ENEG": 11, "NASH": 10, "ENEG_NASH": 9}


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def synthetic_authorities(root: Path) -> tuple[Path, Path, list[dict[str, object]]]:
    labels: list[dict[str, object]] = []
    folds: list[dict[str, object]] = []
    rows: list[dict[str, object]] = []
    for group, count in GROUP_COUNTS.items():
        for index in range(count):
            row_id = f"g274_synthetic_{group.lower()}_{index:02d}"
            fold = index % 5
            rows.append({"row_id": row_id, "source_group": group, "outer_fold": fold})
            labels.append(
                {
                    "row_id": row_id,
                    "source_group": group,
                    "source_group_semantics": f"synthetic_{group.lower()}",
                    "mash_state": "observed" if "NASH" in group else "not_observed",
                    "hbv_state": "observed" if "ENEG" in group else "not_observed",
                    "within_instrument_contrast": (
                        "hiseq_CTRL_vs_ENEG" if group in {"CTRL", "ENEG"} else "novaseq_NASH_vs_ENEG_NASH"
                    ),
                    "expected_instrument": (
                        "Illumina HiSeq 4000" if group in {"CTRL", "ENEG"} else "Illumina NovaSeq 6000"
                    ),
                }
            )
            folds.append(
                {
                    "row_id": row_id,
                    "outer_fold": fold,
                    "source_group": group,
                    "fold_assignment_used_labels": True,
                    "available_to_model_input": False,
                }
            )
    labels_path = root / "evaluator" / "source_labels.tsv"
    folds_path = root / "evaluator" / "participant_folds.tsv"
    write_tsv(
        labels_path,
        (
            "row_id",
            "source_group",
            "source_group_semantics",
            "mash_state",
            "hbv_state",
            "within_instrument_contrast",
            "expected_instrument",
        ),
        labels,
    )
    write_tsv(
        folds_path,
        (
            "row_id",
            "outer_fold",
            "source_group",
            "fold_assignment_used_labels",
            "available_to_model_input",
        ),
        folds,
    )
    return labels_path, folds_path, rows


def synthetic_feature_root(
    root: Path,
    rows: list[dict[str, object]],
    *,
    perturb_query_fold: int | None = None,
) -> Path:
    rng = np.random.default_rng(274_114)
    genes = 64
    counts = rng.poisson(8.0, size=(len(rows), genes)).astype(np.float64) + 1.0
    for index, row in enumerate(rows):
        group = row["source_group"]
        signal = {"CTRL": 0, "ENEG": 1, "NASH": 2, "ENEG_NASH": 3}[str(group)]
        counts[index, signal] += 120.0
        if perturb_query_fold is not None and row["outer_fold"] == perturb_query_fold:
            counts[index] *= 17.0
            counts[index, 10] += 10_000.0
    lengths = np.full_like(counts, 1_000.0)
    feature_root = root / (
        "features" if perturb_query_fold is None else f"features_perturbed_fold_{perturb_query_fold}"
    )
    feature_root.mkdir(parents=True)
    write_tsv(
        feature_root / "participants.tsv",
        ("row_id",),
        [{"row_id": row["row_id"]} for row in rows],
    )
    write_tsv(
        feature_root / "genes.tsv",
        ("gene_id",),
        [{"gene_id": f"ENSG{index + 1:011d}.1"} for index in range(genes)],
    )
    with (feature_root / "numreads.npy").open("xb") as handle:
        np.save(handle, counts, allow_pickle=False)
    with (feature_root / "effective_length.npy").open("xb") as handle:
        np.save(handle, lengths, allow_pickle=False)
    receipt = {
        "schema_version": "masld-bench-gse274114-unlabeled-features-v1",
        "participants": 39,
        "genes": genes,
        "labels_opened": False,
        "folds_opened": False,
        "participant_metadata_opened": False,
        "molecular_values_opened": True,
        "normalization_run": False,
        "model_fit_or_scoring_run": False,
        "effective_length_retained": True,
    }
    write_json_exclusive(feature_root / "receipt.json", receipt)
    freeze_tree(
        feature_root,
        {"artifact_class": "gse274114_unlabeled_feature_matrix", **receipt},
    )
    return feature_root


class GSE274114BaselineFirewallTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.labels, self.folds, self.rows = synthetic_authorities(self.root)
        self.features = synthetic_feature_root(self.root, self.rows)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def make_view(self, task_id: str, fold: int, *, features: Path | None = None) -> Path:
        name = f"view_{task_id}_{fold}_{'base' if features is None else features.name}"
        output = self.root / name
        firewall.prepare_fold_view(
            feature_root=self.features if features is None else features,
            labels_path=self.labels,
            folds_path=self.folds,
            task_id=task_id,
            outer_fold=fold,
            output=output,
        )
        return output

    def test_trainer_view_contains_training_labels_but_no_query_labels_or_covariates(self) -> None:
        view = self.make_view("gse274114_hiseq_healthy_vs_hbv", 0)
        manifest = verify_frozen_tree(view)
        query_header = (view / "query_participants.tsv").read_text(encoding="utf-8").splitlines()[0]
        self.assertEqual(query_header, "row_id")
        self.assertFalse(manifest["metadata"]["query_labels_included"])
        self.assertFalse(manifest["metadata"]["participant_covariates_included"])
        self.assertFalse(manifest["metadata"]["age_included"])
        self.assertFalse(manifest["metadata"]["sex_included"])
        self.assertFalse(manifest["metadata"]["group_aggregate_metadata_included"])
        self.assertNotIn("source_group", query_header)

    def test_query_values_cannot_change_training_preprocessing(self) -> None:
        task_id = "gse274114_hiseq_healthy_vs_hbv"
        base_view = self.make_view(task_id, 0)
        perturbed_features = synthetic_feature_root(
            self.root, self.rows, perturb_query_fold=0
        )
        perturbed_view = self.make_view(task_id, 0, features=perturbed_features)
        outputs = []
        for index, view in enumerate((base_view, perturbed_view)):
            output = self.root / f"nearest_centroid_{index}"
            bundle = firewall.fit_predict(
                view_root=view,
                expected_view_sha256=sha256_file(view / "ARTIFACTS.json"),
                model_kind="gene_rank_nearest_centroid",
                seed=0,
                output=output,
            )
            outputs.append(bundle)
        self.assertEqual(
            outputs[0]["metadata"]["training_preprocessing_sha256"],
            outputs[1]["metadata"]["training_preprocessing_sha256"],
        )
        self.assertFalse(outputs[0]["metadata"]["cohort_wide_normalization_run"])

    def test_prediction_bundle_is_frozen_and_contains_no_outcome(self) -> None:
        view = self.make_view("gse274114_novaseq_mash_vs_mash_hbv", 2)
        output = self.root / "elastic_net"
        firewall.fit_predict(
            view_root=view,
            expected_view_sha256=sha256_file(view / "ARTIFACTS.json"),
            model_kind="hvg_pca_elastic_net",
            seed=1103,
            output=output,
        )
        verify_frozen_tree(output)
        bundle = PredictionBundle.load_json(output / "prediction_bundle.json")
        bundle.validate_artifacts(output)
        fields, _ = firewall.read_tsv(output / "predictions.tsv")
        self.assertEqual(fields, ("row_hash", "unit_hash", "predicted"))
        self.assertFalse(bundle.metadata["query_labels_opened"])
        self.assertFalse(bundle.metadata["scoring_run"])
        self.assertEqual(bundle.metadata["hyperparameter_candidates_total"], 9)
        self.assertEqual(
            bundle.metadata["hyperparameter_candidates_valid"]
            + bundle.metadata["hyperparameter_candidates_invalid"],
            9,
        )
        self.assertEqual(
            bundle.metadata["no_valid_candidate_policy"],
            "rank_inner_valid_candidates_then_fail_outer_fold_if_none_converge",
        )
        self.assertGreaterEqual(bundle.metadata["selected_inner_candidate_rank"], 1)
        self.assertGreaterEqual(bundle.metadata["outer_training_candidates_attempted"], 1)

    def test_convergence_warning_invalidates_candidates_and_fails_closed(self) -> None:
        from sklearn.exceptions import ConvergenceWarning
        from sklearn.linear_model import LogisticRegression

        view = self.make_view("gse274114_hiseq_healthy_vs_hbv", 1)
        train_rows = firewall.read_tsv(view / "training_participants.tsv")[1]
        train_counts = np.load(view / "training_numreads.npy", allow_pickle=False)
        query_counts = np.load(view / "query_numreads.npy", allow_pickle=False)
        targets = np.asarray([int(row["target"]) for row in train_rows])

        def emit_warning(model: object, values: object, labels: object) -> object:
            warnings.warn("synthetic nonconvergence", ConvergenceWarning)
            return model

        with patch.object(LogisticRegression, "fit", new=emit_warning):
            with self.assertRaisesRegex(
                firewall.GSE274114BaselineFirewallError,
                "no nested hyperparameter candidate converged",
            ):
                firewall._predict_baseline(
                    train_counts,
                    query_counts,
                    targets,
                    model_kind="hvg_pca_elastic_net",
                    seed=1103,
                )

    def test_outer_training_nonconvergence_falls_back_by_inner_rank(self) -> None:
        view = self.make_view("gse274114_hiseq_healthy_vs_hbv", 1)
        train_rows = firewall.read_tsv(view / "training_participants.tsv")[1]
        train_counts = np.load(view / "training_numreads.npy", allow_pickle=False)
        query_counts = np.load(view / "query_numreads.npy", allow_pickle=False)
        targets = np.asarray([int(row["target"]) for row in train_rows])
        _, reference = firewall._predict_baseline(
            train_counts,
            query_counts,
            targets,
            model_kind="hvg_pca_elastic_net",
            seed=1103,
        )
        selected = (reference["selected_c"], reference["selected_l1_ratio"])
        original = firewall._fit_model_with_convergence_audit

        def reject_selected_outer_fit(model: object, values: object, labels: object) -> bool:
            if len(labels) == len(targets) and (
                getattr(model, "C", None), getattr(model, "l1_ratio", None)
            ) == selected:
                return False
            return original(model, values, labels)

        with patch.object(
            firewall,
            "_fit_model_with_convergence_audit",
            side_effect=reject_selected_outer_fit,
        ):
            _, receipt = firewall._predict_baseline(
                train_counts,
                query_counts,
                targets,
                model_kind="hvg_pca_elastic_net",
                seed=1103,
            )
        self.assertTrue(receipt["training_only_convergence_fallback_used"])
        self.assertEqual(receipt["selected_inner_candidate_rank"], 2)
        self.assertEqual(receipt["outer_training_candidates_attempted"], 2)
        self.assertEqual(
            receipt["outer_training_candidate_audit"][0]["outer_training_status"],
            "convergence_warning",
        )
        self.assertEqual(
            receipt["outer_training_candidate_audit"][1]["outer_training_status"],
            "converged_without_warning",
        )

    def test_outer_training_nonconvergence_fails_when_all_ranked_candidates_fail(self) -> None:
        view = self.make_view("gse274114_hiseq_healthy_vs_hbv", 1)
        train_rows = firewall.read_tsv(view / "training_participants.tsv")[1]
        train_counts = np.load(view / "training_numreads.npy", allow_pickle=False)
        query_counts = np.load(view / "query_numreads.npy", allow_pickle=False)
        targets = np.asarray([int(row["target"]) for row in train_rows])
        original = firewall._fit_model_with_convergence_audit

        def reject_all_outer_fits(model: object, values: object, labels: object) -> bool:
            if len(labels) == len(targets):
                return False
            return original(model, values, labels)

        with patch.object(
            firewall,
            "_fit_model_with_convergence_audit",
            side_effect=reject_all_outer_fits,
        ):
            with self.assertRaisesRegex(
                firewall.GSE274114BaselineFirewallError,
                "no inner-valid hyperparameter candidate converged on outer training partition",
            ):
                firewall._predict_baseline(
                    train_counts,
                    query_counts,
                    targets,
                    model_kind="hvg_pca_elastic_net",
                    seed=1103,
                )

    def test_incomplete_prediction_set_fails_before_missing_labels_are_opened(self) -> None:
        view = self.make_view("gse274114_hiseq_healthy_vs_hbv", 0)
        output = self.root / "prior_fold_0"
        firewall.fit_predict(
            view_root=view,
            expected_view_sha256=sha256_file(view / "ARTIFACTS.json"),
            model_kind="training_class_prior",
            seed=0,
            output=output,
        )
        index = self.root / "incomplete.tsv"
        write_tsv(
            index,
            ("outer_fold", "seed", "bundle_root", "artifacts_sha256"),
            [
                {
                    "outer_fold": 0,
                    "seed": 0,
                    "bundle_root": output,
                    "artifacts_sha256": sha256_file(output / "ARTIFACTS.json"),
                }
            ],
        )
        with self.assertRaisesRegex(
            firewall.GSE274114BaselineFirewallError, "prediction set is incomplete"
        ):
            firewall.evaluate_model(
                bundle_index=index,
                labels_path=self.root / "does_not_exist.tsv",
                folds_path=self.root / "also_missing.tsv",
                model_kind="training_class_prior",
                task_id="gse274114_hiseq_healthy_vs_hbv",
                bootstrap_replicates=10,
                output=self.root / "must_not_exist",
            )
        self.assertFalse((self.root / "must_not_exist").exists())

    def test_complete_deterministic_prediction_set_scores_only_after_freeze(self) -> None:
        task_id = "gse274114_hiseq_healthy_vs_hbv"
        index_rows: list[dict[str, object]] = []
        for fold in range(5):
            view = self.make_view(task_id, fold)
            output = self.root / f"prior_fold_{fold}"
            firewall.fit_predict(
                view_root=view,
                expected_view_sha256=sha256_file(view / "ARTIFACTS.json"),
                model_kind="training_class_prior",
                seed=0,
                output=output,
            )
            index_rows.append(
                {
                    "outer_fold": fold,
                    "seed": 0,
                    "bundle_root": output,
                    "artifacts_sha256": sha256_file(output / "ARTIFACTS.json"),
                }
            )
        index = self.root / "complete.tsv"
        write_tsv(
            index,
            ("outer_fold", "seed", "bundle_root", "artifacts_sha256"),
            index_rows,
        )
        output = self.root / "evaluation"
        receipt = firewall.evaluate_model(
            bundle_index=index,
            labels_path=self.labels,
            folds_path=self.folds,
            model_kind="training_class_prior",
            task_id=task_id,
            bootstrap_replicates=50,
            output=output,
        )
        verify_frozen_tree(output)
        self.assertTrue(receipt["prediction_bundles_verified_before_labels_opened"])
        self.assertEqual(receipt["prediction_bundle_count"], 5)
        self.assertFalse(receipt["global_four_class_scoring_run"])
        self.assertFalse(receipt["sealed_or_champion_claim_eligible"])

    def test_seed_and_campaign_contracts_are_fail_closed(self) -> None:
        with self.assertRaises(firewall.GSE274114BaselineFirewallError):
            firewall._seed_for("hvg_pca_elastic_net", 7)
        with self.assertRaises(firewall.GSE274114BaselineFirewallError):
            firewall._seed_for("training_class_prior", 1103)
        campaign = tomllib.loads(
            (ROOT / "config/campaigns/gse274114_within_platform_baselines_v2.toml").read_text(
                encoding="utf-8"
            )
        )
        self.assertFalse(campaign["firewall"]["cohort_wide_normalization_allowed"])
        self.assertFalse(campaign["covariates"]["participant_level_metadata_available"])
        self.assertFalse(campaign["global_four_class_task_allowed"])
        self.assertFalse(campaign["sealed_or_champion_claim_eligible"])
        self.assertFalse(
            campaign["firewall"]["outer_candidate_fallback_uses_query_outcomes"]
        )

    def test_feature_staging_path_allowlist_rejects_protected_members(self) -> None:
        valid = firewall._quant_relative_path(
            "model_inputs/quant_v49/g274_example.quant.genes.sf.gz"
        )
        self.assertEqual(valid.parts[:2], ("model_inputs", "quant_v49"))
        for forbidden in (
            "evaluator_only/source_labels.tsv",
            "evaluator_only/participant_folds.tsv",
            "descriptive_only/group_aggregate_metadata.tsv",
            "model_inputs/quant_manifest.tsv",
            "../evaluator_only/source_labels.tsv",
        ):
            with self.subTest(forbidden=forbidden):
                with self.assertRaises(firewall.GSE274114BaselineFirewallError):
                    firewall._quant_relative_path(forbidden)


if __name__ == "__main__":
    unittest.main()

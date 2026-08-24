"""Adversarial tests for deterministic selection and sealed promotion policy."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import masld_bench.tournament as tournament_module
from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.conditional_model import (
    STACKING_EVIDENCE_SCHEMA_VERSION,
    ComplementarityEvidence,
)
from masld_bench.contracts import ContractError, SelectionLock
from masld_bench.evaluators.stats import (
    benjamini_hochberg,
    holm_correction,
    paired_cluster_bootstrap,
    prospective_power_gate,
)
from masld_bench.evaluators.sealed_metrics import (
    SealedMetricError,
    recompute_sealed_metrics,
)
from masld_bench.evaluators.metrics import fisher_z_mean
from masld_bench.hashing import canonical_sha256
from masld_bench.firewall import authorize_one_time_outcome_join
from masld_bench.tournament import (
    verify_development_shortlist,
    freeze_development_shortlist,
    freeze_finalist_candidate_ledger,
    DEVELOPMENT_ENSEMBLE_POLICY_ID,
    UNIVERSE_KIND_SCREENING,
    UNIVERSE_KIND_FINALIST,
    LEDGER_AUTHORITY_SCREENING,
    LEDGER_AUTHORITY_FINALIST,
    freeze_screening_campaign_universe,
    verify_finalist_campaign_universe,
    freeze_finalist_campaign_universe,
    SHORTLIST_SOURCE_WAVE_CONTRACTS,
    family_model_slots,
    BLOCKED_MISSING_EXTERNAL_FAMILY,
    BLOCKED_MISSING_SEALED_SOURCE,
    BULK_TASK,
    CELL_TASK,
    CHAMPION_GATE_THRESHOLDS,
    CONDITIONAL_CONTEXT_TASK,
    CONDITIONAL_SEALED,
    DEVELOPMENT_ONLY,
    EXPLORATORY_ONLY,
    GRAPH_TASK,
    LOCKED_MULTIPLICITY_PLAN,
    PERTURBATION_TASK,
    RNA_ATAC_TASK,
    SEALED_CONFIRMATORY,
    SelectionLockError,
    TournamentError,
    UNIFIED_TASK,
    VARIANT_TASK,
    freeze_champion_gate_decision,
    freeze_confirmatory_multiplicity_bundle,
    freeze_scientific_run_receipt,
    freeze_sealed_metric_bundle,
    freeze_terminal_evaluation_authorization,
    one_standard_error_candidates,
    pareto_frontier,
    require_selection_lock,
    select_family_top_two,
    verify_champion_gate_decision,
    verify_confirmatory_multiplicity_bundle,
    verify_finalist_development_metric_bundle,
    verify_selection_candidate_ledger,
    verify_sealed_metric_bundle,
    verify_terminal_evaluation_authorization,
    verify_variant_secondary_run_receipt,
)
from tests.unit._sealed_fixtures import (
    BASELINE_MODEL as FIXTURE_BASELINE_MODEL,
    SECONDARY_COMPARATOR_SEED,
    SELECTED_MODEL as FIXTURE_SELECTED_MODEL,
    TASK_IDS as FIXTURE_TASK_IDS,
    commit_prediction_set,
    create_joint_outcome_bundle,
    create_scientific_selection,
    freeze_power_set,
)


TASK_IDS = tuple(
    sorted(
        (
            CELL_TASK,
            VARIANT_TASK,
            RNA_ATAC_TASK,
            BULK_TASK,
            GRAPH_TASK,
            PERTURBATION_TASK,
            UNIFIED_TASK,
        )
    )
)
SEEDS = (1103, 2909, 4721, 6673, 8111)
CELL_CLASSES = ("hepatocyte", "cholangiocyte", "stellate")
VARIANT_LINEAGES = ("hepatocyte", "stellate", "immune")
VARIANT_SECONDARY_MODELS = ("abc", "nearest_gene", "re2g")
RNA_ATAC_LINEAGES = (
    "hepatocyte",
    "stellate",
    "macrophage",
    "cholangiocyte",
    "t_nk",
)
BULK_ENDPOINTS = ("fibrosis", "nas")


def digest(label: str) -> str:
    return sha256(label.encode("utf-8")).hexdigest()


def variant_primary_capability() -> dict[str, object]:
    selected_model_id = f"{VARIANT_TASK}_selected"
    baseline_model_id = f"{VARIANT_TASK}_strongest_baseline"

    def record(
        model_id: str, *, role: str, mandatory: bool
    ) -> dict[str, object]:
        return {
            "model_id": model_id,
            "role": role,
            "native_outputs": ["gene_expression_delta"],
            "allowed_endpoints": [
                "allelic_direction",
                "eqtl_retrieval",
                "signed_cell_type_eqtl_effect",
            ],
            "primary_eligible": True,
            "requires_fitted_head": True,
            "requires_observed_target_context": False,
            "is_mandatory_baseline": mandatory,
        }

    selected = record(
        selected_model_id, role="primary_candidate", mandatory=False
    )
    baseline = record(baseline_model_id, role="baseline", mandatory=True)
    return {
        "schema_version": "masld-bench-locked-variant-primary-capability-v1",
        "task_id": VARIANT_TASK,
        "primary_endpoint_id": "signed_cell_type_eqtl_effect",
        "primary_evaluator_id": "variant_ld_block_fisher_z_spearman_gain_v1",
        "capability_registry_sha256": digest("variant-capability-registry"),
        "selected": {
            "model_id": selected_model_id,
            "capability": selected,
            "capability_sha256": canonical_sha256(selected),
        },
        "baseline": {
            "model_id": baseline_model_id,
            "capability": baseline,
            "capability_sha256": canonical_sha256(baseline),
        },
    }


def variant_secondary_evaluation() -> dict[str, object]:
    comparators: list[dict[str, object]] = []
    for model_id in VARIANT_SECONDARY_MODELS:
        capability = {
            "model_id": model_id,
            "role": "baseline" if model_id == "nearest_gene" else "link_only",
            "native_outputs": ["enhancer_gene_link_score"],
            "allowed_endpoints": ["enhancer_gene_link", "eqtl_retrieval"],
            "primary_eligible": False,
            "requires_fitted_head": False,
            "requires_observed_target_context": False,
            "is_mandatory_baseline": True,
        }
        comparators.append(
            {
                "model_id": model_id,
                "run_id": digest(f"{VARIANT_TASK}:secondary:{model_id}"),
                "capability": capability,
                "capability_sha256": canonical_sha256(capability),
                "score_transform_id": "identity_link_score_v1",
            }
        )
    return {
        "schema_version": (
            "masld-bench-locked-variant-secondary-evaluation-v1"
        ),
        "task_id": VARIANT_TASK,
        "endpoint_id": "eqtl_retrieval",
        "evaluator_id": "variant_ld_block_eqtl_retrieval_auprc_v1",
        "capability_registry_sha256": digest("variant-capability-registry"),
        "candidate_transform_id": "absolute_signed_effect_v1",
        "primary_baseline_transform_id": "absolute_signed_effect_v1",
        "mandatory_model_ids": list(VARIANT_SECONDARY_MODELS),
        "comparators": comparators,
    }


RUN_IDS = {
    task_id: {
        role: tuple(digest(f"{task_id}:{role}:{seed}") for seed in SEEDS)
        for role in ("selected", "baseline")
    }
    for task_id in TASK_IDS
}


def power_spec(task_id: str) -> dict[str, object]:
    policy = CHAMPION_GATE_THRESHOLDS[task_id]
    return {
        "metric_id": policy["power_metric"],
        "minimum_effect": float(policy["power_minimum_effect"]),
        "alpha": 0.05,
        "target_power": 0.80,
        "n_primary_claims": 5,
        "min_units": 10,
        "two_sided": True,
    }


def promotion_thresholds(task_id: str, mode: str | None = None) -> dict[str, object]:
    policy = CHAMPION_GATE_THRESHOLDS[task_id]
    if mode is None:
        mode = {
            CELL_TASK: SEALED_CONFIRMATORY,
            VARIANT_TASK: SEALED_CONFIRMATORY,
            RNA_ATAC_TASK: CONDITIONAL_SEALED,
            BULK_TASK: CONDITIONAL_SEALED,
            GRAPH_TASK: BLOCKED_MISSING_SEALED_SOURCE,
            PERTURBATION_TASK: EXPLORATORY_ONLY,
            UNIFIED_TASK: CONDITIONAL_SEALED,
        }[task_id]
    promotion: dict[str, object] = {"promotion_mode": mode}
    if task_id == CELL_TASK:
        promotion.update(
            primary_min=policy["primary_min"],
            delta_min=policy["delta_min"],
            delta_ci_low_strict_min=policy["delta_ci_low_strict_min"],
            stratum_min=policy["stratum_min"],
            brier_delta_max=policy["brier_delta_max"],
            required_class_ids=list(CELL_CLASSES),
        )
    elif task_id == VARIANT_TASK:
        promotion.update(
            delta_min=policy["delta_min"],
            delta_ci_low_strict_min=policy["delta_ci_low_strict_min"],
            stratum_min=policy["stratum_min"],
            required_major_lineage_ids=list(VARIANT_LINEAGES),
            eqtl_auprc_noninferiority_margin=0.01,
        )
    elif task_id == RNA_ATAC_TASK:
        promotion.update(
            delta_min=policy["delta_min"],
            delta_ci_low_strict_min=policy["delta_ci_low_strict_min"],
            stratum_min=policy["stratum_min"],
            min_improved_strata=policy["min_improved_strata"],
            required_major_lineage_ids=list(RNA_ATAC_LINEAGES),
        )
    elif task_id == BULK_TASK:
        promotion.update(
            delta_min=policy["delta_min"],
            delta_ci_low_strict_min=policy["delta_ci_low_strict_min"],
            required_endpoint_ids=list(BULK_ENDPOINTS),
        )
    elif task_id == GRAPH_TASK and mode == CONDITIONAL_SEALED:
        promotion.update(
            delta_min=policy["delta_min"],
            delta_ci_low_strict_min=policy["delta_ci_low_strict_min"],
            brier_delta_max=policy["brier_delta_max"],
            sealed_source_available=True,
        )
    elif task_id == UNIFIED_TASK and mode == CONDITIONAL_SEALED:
        promotion.update(
            min_superior_source_families=policy["min_superior_source_families"],
            min_external_superior_source_families=policy[
                "min_external_superior_source_families"
            ],
            min_eligible_external_source_families=policy[
                "min_eligible_external_source_families"
            ],
            eligible_specialist_ids=[CELL_TASK, VARIANT_TASK],
            specialist_noninferiority_margins={CELL_TASK: 0.01, VARIANT_TASK: 0.02},
            eligible_source_family_ids=[
                "gse289173_cell",
                "gse289173_regulatory",
                "gse267031",
                "fnih_86_liver",
            ],
            eligible_external_source_family_ids=["gse267031", "fnih_86_liver"],
            source_independence_group_by_family={
                "gse289173_cell": "gse289173_joint",
                "gse289173_regulatory": "gse289173_joint",
                "gse267031": "gse267031",
                "fnih_86_liver": "fnih_86_liver",
            },
        )
    if mode in {SEALED_CONFIRMATORY, CONDITIONAL_SEALED}:
        bundle = {
            CELL_TASK: "gse289173_joint_cell_and_regulatory_bundle",
            VARIANT_TASK: "gse289173_joint_cell_and_regulatory_bundle",
            RNA_ATAC_TASK: "fnih_86_liver_rna_atac_bundle",
            BULK_TASK: "gse267031_bulk_bundle",
            GRAPH_TASK: "future_graph_bundle",
            UNIFIED_TASK: "unified_three_source_bundle",
        }[task_id]
        promotion.update(
            confirmatory_family_id=f"{task_id}:primary",
            confirmatory_fwer=0.05,
            sealed_outcome_bundle_id=bundle,
            power=power_spec(task_id),
        )
    if mode == CONDITIONAL_SEALED:
        promotion["conditional_activation"] = {
            "activated": True,
            "activation_receipt_sha256": digest(f"{task_id}:activation"),
        }
    return promotion


def sealed_ids(task_id: str, promotion: dict[str, object]) -> tuple[str, ...]:
    mode = promotion["promotion_mode"]
    if mode not in {SEALED_CONFIRMATORY, CONDITIONAL_SEALED}:
        return ()
    if task_id == GRAPH_TASK:
        return ("future_graph_seal",)
    return tuple(CHAMPION_GATE_THRESHOLDS[task_id]["sealed_holdout_dataset_ids"])


def stable_task_decision(
    task_id: str,
    *,
    selected_ids: tuple[str, ...] | None = None,
    baseline_ids: tuple[str, ...] | None = None,
    promotion: dict[str, object] | None = None,
) -> dict[str, object]:
    selected_ids = RUN_IDS[task_id]["selected"] if selected_ids is None else selected_ids
    baseline_ids = RUN_IDS[task_id]["baseline"] if baseline_ids is None else baseline_ids
    promotion = promotion_thresholds(task_id) if promotion is None else promotion
    frozen_sealed = sealed_ids(task_id, promotion)
    candidate_datasets = tuple(sorted((f"development_{task_id}", *frozen_sealed)))
    return {
        "task_id": task_id,
        "selected_model_id": f"{task_id}_selected",
        "baseline_model_id": f"{task_id}_strongest_baseline",
        "selected_adaptation_regime": "five_seed_final",
        "baseline_adaptation_regime": "fixed_baseline",
        "checkpoint_sha256": digest(f"{task_id}:checkpoint"),
        "preprocessing_sha256": digest(f"{task_id}:preprocessing"),
        "task_head_sha256": digest(f"{task_id}:head"),
        "seeds": list(SEEDS),
        "calibration_sha256": digest(f"{task_id}:calibration"),
        "thresholds": {"promotion": promotion},
        "evaluator_sha256": digest(f"{task_id}:evaluator"),
        "promotion_gate": CHAMPION_GATE_THRESHOLDS[task_id]["promotion_gate_id"],
        "open_champion": promotion["promotion_mode"]
        in {SEALED_CONFIRMATORY, CONDITIONAL_SEALED},
        "selected_runs": [
            {"seed": seed, "run_id": run_id}
            for seed, run_id in zip(SEEDS, selected_ids, strict=True)
        ],
        "baseline_runs": [
            {"seed": seed, "run_id": run_id}
            for seed, run_id in zip(SEEDS, baseline_ids, strict=True)
        ],
        "candidate_dataset_ids": list(candidate_datasets),
        "sealed_dataset_ids": list(frozen_sealed),
        "dataset_registry_sha256s": {
            dataset_id: digest(f"dataset:{dataset_id}")
            for dataset_id in candidate_datasets
        },
        "variant_primary_capability": (
            variant_primary_capability() if task_id == VARIANT_TASK else None
        ),
        "variant_secondary_evaluation": (
            variant_secondary_evaluation() if task_id == VARIANT_TASK else None
        ),
    }


def build_lock(
    decisions: tuple[dict[str, object], ...] | None = None,
    *,
    extra_candidate_ids: tuple[str, ...] = (),
    outcomes_unlocked: bool = False,
) -> SelectionLock:
    if decisions is None:
        decisions = tuple(stable_task_decision(task_id) for task_id in TASK_IDS)
    decisions = tuple(sorted(decisions, key=lambda item: item["task_id"]))
    selected_ids = tuple(
        sorted(
            record["run_id"]
            for decision in decisions
            for record in decision["selected_runs"]
        )
    )
    candidate_ids = tuple(
        sorted(
            {
                *(record["run_id"] for decision in decisions for record in decision["selected_runs"]),
                *(record["run_id"] for decision in decisions for record in decision["baseline_runs"]),
                *(
                    comparator["run_id"]
                    for decision in decisions
                    if isinstance(
                        decision.get("variant_secondary_evaluation"), dict
                    )
                    for comparator in decision[
                        "variant_secondary_evaluation"
                    ]["comparators"]
                ),
                *extra_candidate_ids,
            }
        )
    )
    identity: dict[str, object] = {
        "schema_version": "masld-bench-selection-lock-v1",
        "campaign_id": "fixture-campaign",
        "plan_sha256": digest("plan"),
        "candidate_manifest_sha256": digest("candidate-manifest"),
        "registry_sha256": digest("registry"),
        "metrics_sha256": digest("development-metrics"),
        "locked": True,
        "outcomes_unlocked": outcomes_unlocked,
        "release_state": "candidate",
        "candidate_run_ids": list(candidate_ids),
        "selected_run_ids": list(selected_ids),
        "task_decisions": list(decisions),
        "conditional_model_decision": {
            "evidence_receipt_sha256": digest("conditional:evidence"),
            "comparison_task_ids": [RNA_ATAC_TASK, VARIANT_TASK],
            "created_before_outcome_unblind": True,
            "terminal_track": CONDITIONAL_CONTEXT_TASK,
        },
        "power_decisions": {},
        "multiplicity_plan": LOCKED_MULTIPLICITY_PLAN,
        "terminal_policy": "no_reselection_recalibration_threshold_change_or_repair",
        "metadata": {},
    }
    return SelectionLock.from_dict(
        {"lock_id": canonical_sha256(identity), **identity}
    )


def lock_payload(lock: SelectionLock) -> dict[str, object]:
    payload = lock.to_dict()
    payload.pop("lock_id")
    return payload


def rebuild_lock(payload: dict[str, object]) -> SelectionLock:
    payload = dict(payload)
    return SelectionLock.from_dict(
        {"lock_id": canonical_sha256(payload), **payload}
    )


def change_task(lock: SelectionLock, task_id: str, mutate) -> SelectionLock:
    payload = lock_payload(lock)
    for decision in payload["task_decisions"]:
        if decision["task_id"] == task_id:
            mutate(decision)
            break
    return rebuild_lock(payload)


def decision(lock: SelectionLock, task_id: str) -> dict[str, object]:
    return dict(next(item for item in lock.task_decisions if item["task_id"] == task_id))


def promotion(lock: SelectionLock, task_id: str) -> dict[str, object]:
    return dict(decision(lock, task_id)["thresholds"]["promotion"])


def frozen_run_ids(lock: SelectionLock, task_id: str, role: str) -> tuple[str, ...]:
    return tuple(record["run_id"] for record in decision(lock, task_id)[f"{role}_runs"])


class MultipleTestingAndPowerTests(unittest.TestCase):
    def test_holm_bh_and_cluster_bootstrap_are_deterministic(self) -> None:
        p_values = {"a": 0.01, "b": 0.04, "c": 0.03, "d": 0.002}
        self.assertEqual(
            holm_correction(p_values),
            {"a": 0.03, "b": 0.06, "c": 0.06, "d": 0.008},
        )
        self.assertEqual(
            benjamini_hochberg(p_values),
            {"a": 0.02, "b": 0.04, "c": 0.04, "d": 0.008},
        )
        first = paired_cluster_bootstrap(
            [2, 2, 0, 3, 3, 3],
            [1, 1, 1, 1, 1, 1],
            ["d1", "d1", "d2", "d3", "d3", "d3"],
            n_resamples=500,
            seed=4317,
        )
        second = paired_cluster_bootstrap(
            [2, 2, 0, 3, 3, 3],
            [1, 1, 1, 1, 1, 1],
            ["d1", "d1", "d2", "d3", "d3", "d3"],
            n_resamples=500,
            seed=4317,
        )
        self.assertEqual(first, second)
        self.assertAlmostEqual(first.estimate, 2 / 3)

    def test_power_gate_uses_holm_family_and_unit_floor(self) -> None:
        passing = prospective_power_gate(
            n_units=100,
            paired_sd=0.05,
            minimum_effect=0.05,
            n_primary_claims=5,
            min_units=10,
        )
        self.assertTrue(passing.passed)
        self.assertEqual(passing.effective_alpha, 0.01)
        underpowered = prospective_power_gate(
            n_units=12,
            paired_sd=0.2,
            minimum_effect=0.02,
            n_primary_claims=5,
            min_units=10,
        )
        self.assertFalse(underpowered.passed)


class SelectionTests(unittest.TestCase):
    def selection_fixture(self):
        selected = tuple(digest(f"selected-{index}") for index in range(5))
        baseline = tuple(digest(f"baseline-{index}") for index in range(5))
        candidates = []
        selected_specs = [
            (selected[0], "a", 0.90, 0.10, 10.0),
            (selected[1], "a", 0.88, 0.10, 1.0),
            (selected[2], "b", 0.80, 0.05, 2.0),
            (selected[3], "b", 0.79, 0.05, 1.0),
            (selected[4], "c", 0.75, 0.07, 1.0),
        ]
        for run_id, family, score, calibration, complexity in selected_specs:
            candidates.append(
                {
                    "candidate_id": run_id,
                    "family": family,
                    "metrics": {"score": score, "calibration": calibration},
                    "standard_errors": {"score": 0.03},
                    "complexity": complexity,
                    "eligible": True,
                }
            )
        # The first baseline is Pareto but outside one-SE; the rest are ineligible.
        candidates.append(
            {
                "candidate_id": baseline[0],
                "family": "a",
                "metrics": {"score": 0.70, "calibration": 0.01},
                "standard_errors": {"score": 0.01},
                "complexity": 0.1,
                "eligible": True,
            }
        )
        for index, run_id in enumerate(baseline[1:], start=1):
            candidates.append(
                {
                    "candidate_id": run_id,
                    "family": f"baseline-{index}",
                    "metrics": {"score": 0.99, "calibration": 0.01},
                    "standard_errors": {"score": 0.01},
                    "complexity": 0.1,
                    "eligible": False,
                }
            )
        return candidates, selected

    def test_shortlist_is_top_two_union_pareto_union_one_se(self) -> None:
        candidates, selected = self.selection_fixture()
        family_a = [item for item in candidates if item["family"] == "a"]
        self.assertEqual(
            set(
                pareto_frontier(
                    family_a, {"score": "max", "calibration": "min"}
                )
            ),
            {selected[0], digest("baseline-0")},
        )
        self.assertEqual(
            one_standard_error_candidates(family_a, primary_metric="score"),
            (selected[0], selected[1]),
        )
        observed = select_family_top_two(
            candidates,
            primary_metric="score",
            objectives={"score": "max", "calibration": "min"},
        )
        self.assertEqual(
            observed["a"],
            (selected[0], selected[1], digest("baseline-0")),
        )
        self.assertEqual(
            observed,
            select_family_top_two(
                reversed(candidates),
                primary_metric="score",
                objectives={"score": "max", "calibration": "min"},
            ),
        )

    def test_two_recipes_of_one_model_cannot_take_both_family_slots(self) -> None:
        """Family slots are chosen over models, then resolved to a best recipe."""

        candidates = [
            {
                "candidate_id": digest("m1-recipe-a"),
                "model_id": "model_one",
                "family": "f",
                "metrics": {"score": 0.99, "calibration": 0.01},
                "standard_errors": {"score": 0.001},
                "eligible": True,
            },
            {
                "candidate_id": digest("m1-recipe-b"),
                "model_id": "model_one",
                "family": "f",
                "metrics": {"score": 0.98, "calibration": 0.01},
                "standard_errors": {"score": 0.001},
                "eligible": True,
            },
            {
                "candidate_id": digest("m2-recipe-a"),
                "model_id": "model_two",
                "family": "f",
                "metrics": {"score": 0.50, "calibration": 0.01},
                "standard_errors": {"score": 0.001},
                "eligible": True,
            },
        ]
        slots = family_model_slots(candidates, primary_metric="score")
        self.assertEqual(
            slots["f"], (digest("m1-recipe-a"), digest("m2-recipe-a"))
        )
        self.assertNotIn(digest("m1-recipe-b"), slots["f"])

    def test_every_family_contributes_at_most_two_before_global_union(self) -> None:
        candidates, _ = self.selection_fixture()
        extra = (
            {
                "candidate_id": digest("x-primary"),
                "family": "x",
                "metrics": {"score": 0.90, "calibration": 0.20},
                "standard_errors": {"score": 0.01},
                "complexity": 1.0,
                "eligible": True,
            },
            {
                "candidate_id": digest("x-calibration"),
                "family": "x",
                "metrics": {"score": 0.50, "calibration": 0.01},
                "standard_errors": {"score": 0.01},
                "complexity": 1.0,
                "eligible": True,
            },
        )
        observed = select_family_top_two(
            [*candidates, *extra],
            primary_metric="score",
            objectives={"calibration": "min"},
        )
        self.assertEqual(
            observed["x"],
            (digest("x-primary"), digest("x-calibration")),
        )


class ConditionalTriggerContractTests(unittest.TestCase):
    @staticmethod
    def evidence(
        root: Path,
        *,
        residual_correlation: float = 0.70,
        task_keyed: bool = True,
    ) -> ComplementarityEvidence:
        relative_key = RNA_ATAC_TASK if task_keyed else "profile_deviance"
        absolute_key = VARIANT_TASK if task_keyed else "variant_correlation"
        source_bindings = []
        for role, count in (
            ("development_gain_bundle", 2),
            ("development_residual_bundle", 2),
            ("development_shortlist", 2),
        ):
            artifact_class, filename = {
                "development_gain_bundle": (
                    "development_stack_gain_bundle",
                    "development_stack_gain_bundle.json",
                ),
                "development_residual_bundle": (
                    "development_residual_bundle",
                    "development_residual_bundle.json",
                ),
                "development_shortlist": (
                    "development_shortlist",
                    "development_shortlist.json",
                ),
            }[role]
            for index in range(count):
                source = root / f"stacking-source-{role}-{index}"
                source.mkdir()
                write_json_exclusive(
                    source / filename,
                    {"fixture": True, "role": role, "index": index},
                )
                freeze_tree(source, {"artifact_class": artifact_class})
                source_bindings.append(
                    {
                        "role": role,
                        "artifact_class": artifact_class,
                        "path": source.resolve().as_posix(),
                        "manifest_sha256": tournament_module.sha256_file(
                            source / "ARTIFACTS.json"
                        ),
                        "document_filename": filename,
                        "document_sha256": tournament_module.sha256_file(
                            source / filename
                        ),
                    }
                )
        source_bindings.sort(key=lambda item: (item["role"], item["path"]))
        identity = {
            "schema_version": STACKING_EVIDENCE_SCHEMA_VERSION,
            "derived_evidence": {
                "residual_correlation": residual_correlation,
                "relative_deviance_gains": {relative_key: 0.06},
                "absolute_correlation_or_f1_gains": {absolute_key: 0.03},
                "qualifying_seed_count": 4,
                "evaluated_seed_count": 5,
                "study_count": 2,
                "cross_fitted": True,
                "nonnegative_stack": True,
                "best_open_models_compared": True,
            },
            "source_bindings": source_bindings,
            "created_before_outcome_unblind": True,
            "sealed_results_used": False,
        }
        evidence_id = canonical_sha256(identity)
        evidence_root = root / "stacking-evidence"
        evidence_root.mkdir()
        write_json_exclusive(
            evidence_root / "stacking_evidence.json",
            {"evidence_id": evidence_id, **identity},
        )
        freeze_tree(
            evidence_root,
            {"artifact_class": "stacking_evidence", "evidence_id": evidence_id},
        )
        return ComplementarityEvidence.from_stacking_evidence(evidence_root)

    @staticmethod
    def fixture(root: Path) -> dict[str, object]:
        paths = {
            name: root / name
            for name in (
                "variant-shortlist",
                "rna-atac-shortlist",
                "first-residual",
                "second-residual",
            )
        }
        for path in paths.values():
            path.mkdir()
        variant_shortlist_id = digest("conditional:variant-shortlist")
        rna_atac_shortlist_id = digest("conditional:rna-atac-shortlist")
        sequence_candidate = digest("conditional:sequence-candidate")
        context_candidate = digest("conditional:context-candidate")
        rna_atac_candidate = digest("conditional:rna-atac-candidate")
        variant_shortlist = {
            "task_id": VARIANT_TASK,
            "shortlist_id": variant_shortlist_id,
            "selected_candidate_ids": [sequence_candidate, context_candidate],
            "candidates": [
                {"candidate_id": sequence_candidate, "family": "sequence"},
                {"candidate_id": context_candidate, "family": "context"},
            ],
        }
        rna_atac_shortlist = {
            "task_id": RNA_ATAC_TASK,
            "shortlist_id": rna_atac_shortlist_id,
            "selected_candidate_ids": [rna_atac_candidate],
            "candidates": [
                {"candidate_id": rna_atac_candidate, "family": "profile"}
            ],
        }

        def shortlist_binding(path: Path, shortlist_id: str) -> dict[str, str]:
            return {
                "path": path.resolve().as_posix(),
                "manifest_sha256": digest(f"{shortlist_id}:manifest"),
                "document_sha256": digest(f"{shortlist_id}:document"),
                "shortlist_id": shortlist_id,
            }

        variant_binding = shortlist_binding(
            paths["variant-shortlist"], variant_shortlist_id
        )
        rna_atac_binding = shortlist_binding(
            paths["rna-atac-shortlist"], rna_atac_shortlist_id
        )
        common_residual = {
            "task_id": VARIANT_TASK,
            "development_shortlist_binding": variant_binding,
            "endpoint_evaluator_id": "variant_development_evaluator_v1",
            "endpoint_evaluator_sha256": digest("conditional:variant-evaluator"),
            "dataset_ids": ["variant-development"],
            "dataset_registry_sha256s": {
                "variant-development": digest("conditional:variant-dataset")
            },
            "fold": "held-donor-locus-fold",
            "row_id_field": "row_hash",
            "unit_id_namespace": "variant_to_regulation:held_development_residual_v1",
            "biological_unit": "donor+ld_block",
            "n_rows": 4,
            "row_set_sha256": digest("conditional:held-row-set"),
        }
        residuals = [
            {
                **common_residual,
                "residual_bundle_id": digest("conditional:first-residual"),
                "candidate_id": sequence_candidate,
                "source_family_id": "sequence",
            },
            {
                **common_residual,
                "residual_bundle_id": digest("conditional:second-residual"),
                "candidate_id": context_candidate,
                "source_family_id": "context",
            },
        ]
        return {
            "paths": paths,
            "shortlists": [variant_shortlist, rna_atac_shortlist],
            "residuals": residuals,
            "variant_binding": variant_binding,
            "rna_atac_binding": rna_atac_binding,
            "sequence_candidate": sequence_candidate,
            "context_candidate": context_candidate,
            "rna_atac_candidate": rna_atac_candidate,
        }

    def trigger(self, fixture: dict[str, object], evidence: ComplementarityEvidence):
        paths = fixture["paths"]
        with (
            patch.object(
                tournament_module,
                "verify_development_shortlist",
                side_effect=fixture["shortlists"],
            ),
            patch.object(
                tournament_module,
                "verify_development_residual_bundle",
                side_effect=fixture["residuals"],
            ),
            patch.object(
                tournament_module,
                "_pearson_residual_correlation",
                return_value=evidence.residual_correlation,
            ) as correlation,
        ):
            decision = tournament_module.conditional_new_model_trigger(
                [paths["variant-shortlist"], paths["rna-atac-shortlist"]],
                development_residual_bundle_paths=[
                    paths["first-residual"],
                    paths["second-residual"],
                ],
                complementarity_evidence=evidence,
                evidence_receipt_sha256=evidence.stacking_evidence_binding[
                    "manifest_sha256"
                ],
            )
        return decision, correlation

    def test_within_task_pair_records_candidate_ids_and_binds_both_task_gains(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self.fixture(root)
            decision, correlation = self.trigger(fixture, self.evidence(root))
        self.assertTrue(decision.triggered)
        self.assertEqual(
            decision.passing_tracks,
            tuple(sorted((VARIANT_TASK, RNA_ATAC_TASK))),
        )
        self.assertEqual(
            decision.complementary_pairs,
            (
                tuple(
                    sorted(
                        (
                            fixture["sequence_candidate"],
                            fixture["context_candidate"],
                        )
                    )
                ),
            ),
        )
        self.assertNotEqual(
            decision.complementary_pairs,
            (tuple(sorted((VARIANT_TASK, RNA_ATAC_TASK))),),
        )
        correlation.assert_called_once()

    def test_old_cross_task_residual_correlation_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self.fixture(root)
            cross_task = {
                **fixture["residuals"][1],
                "task_id": RNA_ATAC_TASK,
                "development_shortlist_binding": fixture["rna_atac_binding"],
                "candidate_id": fixture["rna_atac_candidate"],
                "source_family_id": "profile",
            }
            fixture["residuals"] = [fixture["residuals"][0], cross_task]
            with self.assertRaisesRegex(TournamentError, "cross-task"):
                self.trigger(fixture, self.evidence(root))

    def test_residual_pair_must_bind_one_same_comparison_shortlist(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self.fixture(root)
            second = {
                **fixture["residuals"][1],
                "development_shortlist_binding": fixture["rna_atac_binding"],
            }
            fixture["residuals"] = [fixture["residuals"][0], second]
            with self.assertRaisesRegex(TournamentError, "same frozen"):
                self.trigger(fixture, self.evidence(root))

    def test_residual_pair_requires_distinct_selected_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self.fixture(root)
            duplicate = {
                **fixture["residuals"][1],
                "candidate_id": fixture["sequence_candidate"],
                "source_family_id": "sequence",
            }
            fixture["residuals"] = [fixture["residuals"][0], duplicate]
            with self.assertRaisesRegex(TournamentError, "distinct selected"):
                self.trigger(fixture, self.evidence(root))

    def test_qualifying_gain_keys_must_name_both_comparison_tasks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self.fixture(root)
            decision, _ = self.trigger(
                fixture, self.evidence(root, task_keyed=False)
            )
        self.assertFalse(decision.triggered)
        self.assertEqual(decision.passing_tracks, ())
        self.assertEqual(decision.complementary_pairs, ())
        self.assertTrue(
            any("both comparison task IDs" in reason for reason in decision.reasons)
        )

    def test_strong_negative_residual_correlation_passes_signed_threshold(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self.fixture(root)
            decision, _ = self.trigger(
                fixture,
                self.evidence(root, residual_correlation=-0.95),
            )
        self.assertTrue(decision.triggered)
        self.assertEqual(len(decision.complementary_pairs), 1)

    def test_old_variant_and_cell_comparison_contract_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            evidence = self.evidence(Path(tmp))
            with self.assertRaisesRegex(
                TournamentError, "exactly variant_to_regulation"
            ):
                tournament_module.conditional_new_model_trigger(
                    [],
                    development_residual_bundle_paths=[],
                    complementarity_evidence=evidence,
                    evidence_receipt_sha256=evidence.stacking_evidence_binding[
                        "manifest_sha256"
                    ],
                    comparison_task_ids=(VARIANT_TASK, CELL_TASK),
                )
        self.assertEqual(
            tournament_module.CONDITIONAL_DECISION_SCHEMA_VERSION,
            "masld-bench-conditional-decision-v4",
        )


class StrictLockTests(unittest.TestCase):
    def test_selected_run_ids_must_be_candidate_run_ids(self) -> None:
        payload = lock_payload(build_lock())
        removed = payload["selected_run_ids"][0]
        payload["candidate_run_ids"].remove(removed)
        with self.assertRaisesRegex(ContractError, "subset"):
            SelectionLock.from_dict(
                {"lock_id": canonical_sha256(payload), **payload}
            )

    def test_exact_mapping_accepted_permissive_mapping_rejected(self) -> None:
        lock = build_lock()
        self.assertEqual(require_selection_lock(lock.to_dict()), lock.lock_id)
        with self.assertRaises(SelectionLockError):
            require_selection_lock(
                {"schema_version": "old", "locked": True, "lock_id": "lock"}
            )

    def test_variant_signed_endpoint_rejects_ineligible_capability(self) -> None:
        def make_ineligible(task_decision: dict[str, object]) -> None:
            binding = task_decision["variant_primary_capability"]
            selected = binding["selected"]
            capability = selected["capability"]
            capability["primary_eligible"] = False
            capability["allowed_endpoints"].remove(
                "signed_cell_type_eqtl_effect"
            )
            selected["capability_sha256"] = canonical_sha256(capability)

        lock = change_task(build_lock(), VARIANT_TASK, make_ineligible)
        with self.assertRaisesRegex(
            SelectionLockError, "cannot enter signed-effect scoring"
        ):
            require_selection_lock(lock)

    def test_nonvariant_task_rejects_variant_capability_binding(self) -> None:
        lock = change_task(
            build_lock(),
            CELL_TASK,
            lambda task_decision: task_decision.__setitem__(
                "variant_primary_capability", variant_primary_capability()
            ),
        )
        with self.assertRaisesRegex(
            SelectionLockError, "variant_primary_capability must be null"
        ):
            require_selection_lock(lock)

    def test_post_unblind_tampered_and_alias_task_locks_fail_closed(self) -> None:
        with self.assertRaises(SelectionLockError):
            require_selection_lock(build_lock(outcomes_unlocked=True))

        payload = lock_payload(build_lock())
        payload["metrics_sha256"] = digest("tampered")
        tampered = SelectionLock.from_dict(
            {"lock_id": digest("stale-lock-id"), **payload}
        )
        with self.assertRaisesRegex(SelectionLockError, "identity hash mismatch"):
            require_selection_lock(tampered)

        payload = lock_payload(build_lock())
        payload["task_decisions"][-1]["task_id"] = "invented_alias"
        payload["task_decisions"].sort(key=lambda item: item["task_id"])
        alias = rebuild_lock(payload)
        with self.assertRaisesRegex(SelectionLockError, "unsupported task_id"):
            require_selection_lock(alias)

    def test_stable_task_schema_and_five_seed_arrays_are_required(self) -> None:
        lock = build_lock()
        extra = change_task(lock, CELL_TASK, lambda item: item.update(extra_field=True))
        with self.assertRaisesRegex(SelectionLockError, "stable schema"):
            require_selection_lock(extra)

        def shorten(item):
            item["selected_runs"] = item["selected_runs"][:4]

        short = change_task(lock, CELL_TASK, shorten)
        with self.assertRaisesRegex(SelectionLockError, "all five"):
            require_selection_lock(short)


class SealedMetricCapabilityTests(unittest.TestCase):
    def test_graph_development_endpoint_excludes_sealed_retrieval_payload(self) -> None:
        self.assertEqual(
            tournament_module._JOINED_ENDPOINT_FIELDS[GRAPH_TASK],
            (
                "row_hash",
                "unit_hash",
                "block_hash",
                "observed_binary",
                "candidate",
                "baseline",
            ),
        )

    def test_variant_capability_is_checked_before_any_outcome_rows(self) -> None:
        with self.assertRaisesRegex(
            SealedMetricError, "locked primary-capability binding"
        ):
            recompute_sealed_metrics(
                task_id=VARIANT_TASK,
                selected_rows_by_seed={},
                baseline_rows_by_seed={},
                outcome_rows=[],
                secondary_outcome_rows=[],
                endpoint_parameters={"strata": []},
                n_resamples=10_000,
                bootstrap_seed=4317,
            )

    def test_variant_evaluator_rejects_observed_context_capability(self) -> None:
        binding = variant_primary_capability()
        selected = binding["selected"]
        capability = selected["capability"]
        capability["role"] = "observed_context_only"
        capability["requires_observed_target_context"] = True
        selected["capability_sha256"] = canonical_sha256(capability)
        with self.assertRaisesRegex(
            SealedMetricError, "invalid primary role|ineligible"
        ):
            recompute_sealed_metrics(
                task_id=VARIANT_TASK,
                selected_rows_by_seed={},
                baseline_rows_by_seed={},
                outcome_rows=[],
                secondary_outcome_rows=[],
                endpoint_parameters={"strata": []},
                n_resamples=10_000,
                bootstrap_seed=4317,
                variant_primary_capability=binding,
            )

    def test_variant_evaluator_requires_exact_secondary_binding_before_outcomes(
        self,
    ) -> None:
        with self.assertRaisesRegex(
            SealedMetricError, "locked secondary-evaluation binding"
        ):
            recompute_sealed_metrics(
                task_id=VARIANT_TASK,
                selected_rows_by_seed={},
                baseline_rows_by_seed={},
                outcome_rows=[],
                secondary_outcome_rows=[],
                endpoint_parameters={"strata": []},
                n_resamples=10_000,
                bootstrap_seed=4317,
                variant_primary_capability=variant_primary_capability(),
            )

        reordered = variant_secondary_evaluation()
        reordered["comparators"].reverse()
        with self.assertRaisesRegex(
            SealedMetricError, "not canonical"
        ):
            recompute_sealed_metrics(
                task_id=VARIANT_TASK,
                selected_rows_by_seed={},
                baseline_rows_by_seed={},
                outcome_rows=[],
                secondary_outcome_rows=[],
                endpoint_parameters={"strata": []},
                n_resamples=10_000,
                bootstrap_seed=4317,
                variant_primary_capability=variant_primary_capability(),
                variant_secondary_evaluation=reordered,
            )

    def test_variant_retrieval_uses_native_scores_and_stronger_link_baseline(
        self,
    ) -> None:
        strata = tuple(sorted(VARIANT_LINEAGES))
        selected_rows: list[dict[str, object]] = []
        baseline_rows: list[dict[str, object]] = []
        outcome_rows: list[dict[str, object]] = []
        secondary_outcomes: list[dict[str, object]] = []
        secondary_rows = {
            model_id: [] for model_id in VARIANT_SECONDARY_MODELS
        }
        for block_index in range(4):
            block_hash = digest(f"variant-metric-block:{block_index}")
            for stratum_index, stratum in enumerate(strata):
                direction = -1.0 if stratum_index == 0 else 1.0
                for position in range(2):
                    magnitude = float(block_index * 2 + position + 1)
                    observed = direction * magnitude
                    row_hash = digest(
                        f"variant-metric-row:{block_index}:{stratum}:{position}"
                    )
                    metadata = {
                        "row_hash": row_hash,
                        "unit_hash": digest(
                            f"variant-metric-unit:{block_index}:{stratum}:{position}"
                        ),
                        "block_hash": block_hash,
                        "stratum": stratum,
                    }
                    selected_rows.append({**metadata, "predicted": observed})
                    baseline_rows.append({**metadata, "predicted": -observed})
                    outcome_rows.append({**metadata, "observed": observed})
                    secondary_outcomes.append(
                        {
                            **metadata,
                            "observed_binary": (
                                "1" if block_index % 2 == 0 else "0"
                            ),
                            "summary_statistics_complete": "true",
                        }
                    )
                    for model_id in VARIANT_SECONDARY_MODELS:
                        secondary_rows[model_id].append(
                            {
                                **metadata,
                                "predicted": (
                                    (1.0 if block_index % 2 == 0 else 0.0)
                                    if model_id == "abc"
                                    else 1.0 / magnitude
                                ),
                            }
                        )
        selected_by_seed = {seed: selected_rows for seed in SEEDS}
        baseline_by_seed = {seed: baseline_rows for seed in SEEDS}
        result = recompute_sealed_metrics(
            task_id=VARIANT_TASK,
            selected_rows_by_seed=selected_by_seed,
            baseline_rows_by_seed=baseline_by_seed,
            outcome_rows=outcome_rows,
            secondary_outcome_rows=secondary_outcomes,
            endpoint_parameters={"strata": list(strata)},
            n_resamples=10_000,
            bootstrap_seed=4317,
            variant_primary_capability=variant_primary_capability(),
            variant_secondary_evaluation=variant_secondary_evaluation(),
            variant_secondary_rows_by_model=secondary_rows,
        )
        expected_signed_delta = fisher_z_mean([1.0] * len(strata)) - fisher_z_mean(
            [-1.0] * len(strata)
        )
        self.assertAlmostEqual(
            result.metrics["delta_mean_fisher_z_celltype_spearman"],
            expected_signed_delta,
        )
        selected_model = f"{VARIANT_TASK}_selected"
        baseline_model = f"{VARIANT_TASK}_strongest_baseline"
        self.assertEqual(
            result.metrics["retrieval_score_transform_by_model"],
            {
                selected_model: "absolute_signed_effect_v1",
                baseline_model: "absolute_signed_effect_v1",
                "abc": "identity_link_score_v1",
                "nearest_gene": "identity_link_score_v1",
                "re2g": "identity_link_score_v1",
            },
        )
        self.assertEqual(
            result.metrics["retrieval_auprc_by_model"][selected_model],
            result.metrics["retrieval_auprc_by_model"][baseline_model],
        )
        self.assertEqual(result.metrics["retrieval_auprc_by_model"]["abc"], 1.0)
        self.assertEqual(
            result.metrics["best_eqtl_retrieval_comparator_model_id"], "abc"
        )
        self.assertLess(
            result.metrics["delta_vs_best_retrieval_comparator_auprc"], 0.0
        )
        self.assertLess(
            result.metrics[
                "delta_vs_best_retrieval_comparator_auprc_ci_low"
            ],
            -0.01,
        )
        negative = next(
            row for row in result.ensemble_rows if float(row["candidate"]) < 0.0
        )
        negative_scores = json.loads(negative["retrieval_scores_json"])
        self.assertEqual(
            negative_scores["candidate"], abs(float(negative["candidate"]))
        )


class SealedArtifactEnvelopeTamperTests(unittest.TestCase):
    def test_multiplicity_document_tamper_is_rejected_before_semantic_parse(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "multiplicity"
            root.mkdir()
            document = root / "confirmatory_multiplicity_bundle.json"
            write_json_exclusive(document, {"method": "Holm"})
            freeze_tree(
                root,
                {
                    "artifact_class": "confirmatory_multiplicity_bundle",
                    "multiplicity_bundle_id": canonical_sha256({"fixture": True}),
                },
            )
            document.chmod(0o640)
            document.write_text(
                document.read_text(encoding="utf-8") + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                TournamentError, "invalid confirmatory multiplicity bundle"
            ):
                verify_confirmatory_multiplicity_bundle(
                    root, reverify_sources=True
                )


class DevelopmentEnsembleRankingTests(unittest.TestCase):
    """The ranked quantity must be the deployed five-seed ensemble."""

    @classmethod
    def setUpClass(cls) -> None:
        cls._temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls._temporary.name)
        cls.fixture = create_scientific_selection(cls.root)

    @classmethod
    def tearDownClass(cls) -> None:
        cls._temporary.cleanup()

    def test_finalist_bundle_reports_the_ensemble_not_a_mean_of_seeds(
        self,
    ) -> None:
        for task_id, bundle_dir in self.fixture.finalist_metric_bundle_dirs.items():
            payload = verify_finalist_development_metric_bundle(
                bundle_dir, reverify_sources=True
            )
            result = json.loads(
                (
                    bundle_dir / payload["endpoint_result"]["path"]
                ).read_text(encoding="utf-8")
            )
            self.assertEqual(
                result["schema_version"],
                "masld-bench-finalist-endpoint-distribution-v2",
                task_id,
            )
            self.assertEqual(
                result["aggregation_method"],
                "task_native_endpoint_on_mean_of_exactly_five_locked_seed_"
                "predictions",
                task_id,
            )
            self.assertEqual(
                result["development_ensemble_policy_id"],
                DEVELOPMENT_ENSEMBLE_POLICY_ID,
                task_id,
            )
            stability = result["seed_stability"]
            self.assertFalse(stability["used_for_ranking"], task_id)
            self.assertFalse(
                stability["seeds_are_biological_replicates"], task_id
            )
            self.assertEqual(stability["evaluated_seed_count"], 5, task_id)
            self.assertGreaterEqual(stability["positive_seed_count"], 4, task_id)

    def test_cell_ensemble_beats_the_mean_of_per_seed_macro_f1(self) -> None:
        """The fixture is built so the two estimands provably disagree.

        Each seed misclassifies a rotating sixth of the rows, so every per-seed
        macro-F1 is below one, while averaging the class probabilities first
        recovers the correct class on every row.  If ranking ever reverts to a
        mean of per-seed metrics this test fails.
        """

        task_id = "cell_state_mapping"
        payload = verify_finalist_development_metric_bundle(
            self.fixture.finalist_metric_bundle_dirs[task_id],
            reverify_sources=False,
        )
        result = json.loads(
            (
                self.fixture.finalist_metric_bundle_dirs[task_id]
                / payload["endpoint_result"]["path"]
            ).read_text(encoding="utf-8")
        )
        per_seed = result["seed_stability"]["per_seed_absolute_primary_metric"]
        self.assertEqual(len(per_seed), 5)
        mean_of_per_seed = sum(per_seed.values()) / len(per_seed)
        self.assertGreater(
            result["candidate_primary_metric"],
            mean_of_per_seed,
            "ensemble macro-F1 must exceed the mean of per-seed macro-F1",
        )

    def test_ledger_candidates_declare_the_ensemble_metric_basis(self) -> None:
        ledger = verify_selection_candidate_ledger(
            self.fixture.selection_candidate_ledger_dir
        )
        self.assertEqual(ledger["ledger_authority"], LEDGER_AUTHORITY_FINALIST)
        self.assertTrue(ledger["lock_construction_allowed"])
        self.assertEqual(ledger["campaign_universe_kind"], UNIVERSE_KIND_FINALIST)
        self.assertEqual(ledger["expected_seed_count"], 5)
        five_seed = [
            candidate
            for candidate in ledger["model_candidates"]
            if candidate["five_seed_universe_complete"]
        ]
        self.assertTrue(five_seed)
        for candidate in five_seed:
            self.assertEqual(
                candidate["metric_basis"],
                DEVELOPMENT_ENSEMBLE_POLICY_ID,
                candidate["candidate_id"],
            )
            artifacts = candidate["selection_artifact_ensemble_sha256s"]
            self.assertEqual(
                artifacts["checkpoint_sha256"], artifacts["task_head_sha256"]
            )
            self.assertEqual(
                artifacts["checkpoint_sha256"], artifacts["calibration_sha256"]
            )
            self.assertNotEqual(
                artifacts["checkpoint_sha256"], artifacts["preprocessing_sha256"]
            )

    def test_finalist_universe_was_frozen_before_scoring(self) -> None:
        universe = verify_finalist_campaign_universe(
            self.fixture.finalist_campaign_universe_dir
        )
        self.assertTrue(universe["created_before_development_scoring"])
        self.assertFalse(universe["development_outcomes_used"])
        self.assertFalse(universe["sealed_results_used"])
        self.assertEqual(universe["expected_seed_count"], 5)
        ledger = verify_selection_candidate_ledger(
            self.fixture.selection_candidate_ledger_dir
        )
        self.assertEqual(
            ledger["reviewed_universe_id"], universe["universe_id"]
        )

    def test_a_wrong_reviewed_universe_id_is_rejected(self) -> None:
        with self.assertRaises(TournamentError):
            freeze_finalist_candidate_ledger(
                finalist_campaign_universe_dir=(
                    self.fixture.finalist_campaign_universe_dir
                ),
                reviewed_universe_id=digest("not-the-reviewed-universe"),
                scientific_run_receipt_dirs=self.fixture.scientific_receipt_dirs,
                variant_secondary_run_receipt_dirs=(
                    self.fixture.variant_secondary_receipt_dirs
                ),
                output_root=self.root / "rejected_ledgers",
            )


class DevelopmentShortlistTests(unittest.TestCase):
    """The shortlist path had zero callers and could never succeed.

    freeze_development_shortlist wrote selection_rule v3 while its verifier
    demanded v2, and the verifier runs after publish_directory_noreplace, so
    every attempt stranded an immutable identity.  These tests execute the path
    for the first time.
    """

    @classmethod
    def setUpClass(cls) -> None:
        cls._temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls._temporary.name)
        cls.fixture = create_scientific_selection(cls.root)
        cls.objectives = {"absolute_primary_metric": "max"}

    @classmethod
    def tearDownClass(cls) -> None:
        cls._temporary.cleanup()

    def test_shortlist_freezes_and_reverifies_from_its_source_wave(self) -> None:
        frozen = freeze_development_shortlist(
            selection_candidate_ledger_dir=(
                self.fixture.selection_candidate_ledger_dir
            ),
            task_id=FIXTURE_TASK_IDS[0],
            source_wave="full_specialist_screen",
            primary_metric="absolute_primary_metric",
            objectives=self.objectives,
            output_root=self.root / "shortlists",
        )
        payload = verify_development_shortlist(frozen)
        self.assertEqual(payload["source_wave"], "full_specialist_screen")
        self.assertEqual(payload["authorized_downstream_wave"], "adaptation")
        self.assertEqual(
            list(payload["source_task_seeds"]), [1103, 2909, 4721, 6673, 8111]
        )
        self.assertTrue(payload["selected_candidate_ids"])

    def test_a_wave_that_may_not_shortlist_is_rejected(self) -> None:
        output = self.root / "rejected_wave_shortlists"
        with self.assertRaisesRegex(
            TournamentError, "may not create a development shortlist"
        ):
            freeze_development_shortlist(
                selection_candidate_ledger_dir=(
                    self.fixture.selection_candidate_ledger_dir
                ),
                task_id=FIXTURE_TASK_IDS[0],
                source_wave="adaptation",
                primary_metric="absolute_primary_metric",
                objectives=self.objectives,
                output_root=output,
            )
        # The guard precedes any mkdir, so nothing may have been published.
        self.assertFalse(output.exists() and any(output.iterdir()))


class LedgerCherryPickGuardTests(unittest.TestCase):
    """The most direct seed-cherry-picking attack had no test at all."""

    @classmethod
    def setUpClass(cls) -> None:
        cls._temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls._temporary.name)
        cls.fixture = create_scientific_selection(cls.root)
        cls.universe = verify_finalist_campaign_universe(
            cls.fixture.finalist_campaign_universe_dir
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls._temporary.cleanup()

    def test_dropping_a_scheduled_receipt_is_rejected(self) -> None:
        """Omitting one seed's receipt must not yield a smaller ledger."""

        tournament_module._clear_ledger_verification_cache_for_testing()
        with self.assertRaises(TournamentError) as caught:
            freeze_finalist_candidate_ledger(
                finalist_campaign_universe_dir=(
                    self.fixture.finalist_campaign_universe_dir
                ),
                reviewed_universe_id=self.universe["universe_id"],
                scientific_run_receipt_dirs=(
                    self.fixture.scientific_receipt_dirs[:-1]
                ),
                variant_secondary_run_receipt_dirs=(
                    self.fixture.variant_secondary_receipt_dirs
                ),
                output_root=self.root / "omitted_run_ledgers",
            )
        self.assertIn("omits scheduled primary runs", str(caught.exception))


class LedgerVerificationCacheTests(unittest.TestCase):
    """Memoizing verification must not weaken what verification proves."""

    @classmethod
    def setUpClass(cls) -> None:
        cls._temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls._temporary.name)
        cls.fixture = create_scientific_selection(cls.root)

    @classmethod
    def tearDownClass(cls) -> None:
        cls._temporary.cleanup()

    def test_a_different_tree_is_rederived_rather_than_served(self) -> None:
        tournament_module._clear_ledger_verification_cache_for_testing()
        source = self.fixture.selection_candidate_ledger_dir
        first = verify_selection_candidate_ledger(source)
        copy = self.root / "ledger_copy"
        shutil.rmtree(copy, ignore_errors=True)
        shutil.copytree(source, copy)
        second = verify_selection_candidate_ledger(copy)
        self.assertEqual(first["ledger_id"], second["ledger_id"])
        self.assertEqual(
            len(tournament_module._LEDGER_VERIFICATION_CACHE), 2
        )

    def test_in_place_tampering_is_never_served_from_cache(self) -> None:
        """The cache key must not be blind to a mutated member file.

        ARTIFACTS.json only LISTS member digests, so editing the ledger
        document in place leaves that manifest's own bytes untouched.  Keying
        on it alone would return a cached payload for a tampered tree.  Operate
        on a copy: mutating the shared fixture would corrupt every later test.
        """

        tournament_module._clear_ledger_verification_cache_for_testing()
        copy = self.root / "ledger_tamper"
        shutil.rmtree(copy, ignore_errors=True)
        shutil.copytree(self.fixture.selection_candidate_ledger_dir, copy)
        verify_selection_candidate_ledger(copy)
        document = copy / "selection_candidate_ledger.json"
        document.chmod(0o600)
        with document.open("a", encoding="utf-8") as handle:
            handle.write(" ")
        with self.assertRaises(TournamentError):
            verify_selection_candidate_ledger(copy)


class FrozenTournamentChainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls._temporary.name)
        cls.fixture = create_scientific_selection(cls.root)
        cls.state = cls.root / "evaluator_state"
        cls.commits = commit_prediction_set(
            cls.root / "sealed_predictions",
            selection_lock_dir=cls.fixture.selection_lock_dir,
            evaluator_state_dir=cls.state,
        )
        cls.powers = freeze_power_set(
            evaluator_state_dir=cls.state,
            commits=cls.commits,
        )
        cls.ready_gates: dict[str, Path] = {}
        for task_id in FIXTURE_TASK_IDS:
            cls.ready_gates[task_id] = freeze_champion_gate_decision(
                selection_lock_dir=cls.fixture.selection_lock_dir,
                task_id=task_id,
                prediction_commit_dirs=[
                    cls.state
                    / "prediction_commits"
                    / commit.prediction_commit_sha256
                    for commit in cls.commits[task_id]
                ],
                power_decision_dir=(
                    cls.state
                    / "power_decisions"
                    / cls.powers[task_id].power_decision_sha256
                ),
                output_root=cls.root / "ready_gates" / task_id,
            )
        cls.outcome = create_joint_outcome_bundle(cls.root / "joint_outcome")
        consumption_ledger = cls.root / "consumption_ledger"
        consumption_ledger.mkdir()
        cls.consumption = authorize_one_time_outcome_join(
            evaluator_state_dir=cls.state,
            power_decision_sha256s=(
                cls.powers[task_id].power_decision_sha256
                for task_id in FIXTURE_TASK_IDS
            ),
            sealed_outcome_bundle_dir=cls.outcome,
            consumption_ledger_dir=consumption_ledger,
        )
        cls.opened_gates: dict[str, Path] = {}
        for task_id in FIXTURE_TASK_IDS:
            cls.opened_gates[task_id] = freeze_champion_gate_decision(
                selection_lock_dir=cls.fixture.selection_lock_dir,
                task_id=task_id,
                prediction_commit_dirs=[
                    cls.state
                    / "prediction_commits"
                    / commit.prediction_commit_sha256
                    for commit in cls.commits[task_id]
                ],
                power_decision_dir=(
                    cls.state
                    / "power_decisions"
                    / cls.powers[task_id].power_decision_sha256
                ),
                sealed_outcome_bundle_dir=cls.outcome,
                outcome_consumption_path=cls.consumption.marker_path,
                output_root=cls.root / "opened_gates" / task_id,
            )
        cls.sealed_metrics: dict[str, Path] = {}
        for task_id in FIXTURE_TASK_IDS:
            cls.sealed_metrics[task_id] = freeze_sealed_metric_bundle(
                selection_lock_dir=cls.fixture.selection_lock_dir,
                task_id=task_id,
                prediction_commit_dirs=[
                    cls.state
                    / "prediction_commits"
                    / commit.prediction_commit_sha256
                    for commit in cls.commits[task_id]
                ],
                power_decision_dir=(
                    cls.state
                    / "power_decisions"
                    / cls.powers[task_id].power_decision_sha256
                ),
                sealed_outcome_bundle_dir=cls.outcome,
                outcome_consumption_path=cls.consumption.marker_path,
                output_root=cls.root / "sealed_metrics" / task_id,
            )
        cls.multiplicity = freeze_confirmatory_multiplicity_bundle(
            selection_lock_dir=cls.fixture.selection_lock_dir,
            sealed_metric_bundle_dirs=cls.sealed_metrics.values(),
            output_root=cls.root / "confirmatory_multiplicity",
        )
        cls.promoted_gates: dict[str, Path] = {}
        for task_id in FIXTURE_TASK_IDS:
            cls.promoted_gates[task_id] = freeze_champion_gate_decision(
                selection_lock_dir=cls.fixture.selection_lock_dir,
                task_id=task_id,
                prediction_commit_dirs=[
                    cls.state
                    / "prediction_commits"
                    / commit.prediction_commit_sha256
                    for commit in cls.commits[task_id]
                ],
                power_decision_dir=(
                    cls.state
                    / "power_decisions"
                    / cls.powers[task_id].power_decision_sha256
                ),
                sealed_outcome_bundle_dir=cls.outcome,
                outcome_consumption_path=cls.consumption.marker_path,
                sealed_metric_bundle_dir=cls.sealed_metrics[task_id],
                confirmatory_multiplicity_bundle_dir=cls.multiplicity,
                output_root=cls.root / "promoted_gates" / task_id,
            )

    @classmethod
    def tearDownClass(cls) -> None:
        cls._temporary.cleanup()

    def test_scientific_ledger_and_finalists_rederive(self) -> None:
        ledger = verify_selection_candidate_ledger(
            self.fixture.selection_candidate_ledger_dir
        )
        self.assertEqual(len(ledger["scientific_receipt_bindings"]), 10)
        self.assertEqual(len(ledger["variant_secondary_receipt_bindings"]), 3)
        self.assertEqual(len(ledger["variant_secondary_runs"]), 3)
        auxiliary_model_ids = {
            record["model_id"] for record in ledger["variant_secondary_runs"]
        }
        self.assertEqual(auxiliary_model_ids, set(VARIANT_SECONDARY_MODELS))
        self.assertTrue(
            all(
                record["outcomes_read"] is False
                and record["metrics_computed"] is False
                and record["auxiliary_status"] == "succeeded_non_scoring"
                for record in ledger["variant_secondary_runs"]
            )
        )
        self.assertTrue(
            auxiliary_model_ids.isdisjoint(
                record["model_id"] for record in ledger["scientific_runs"]
            )
        )
        self.assertTrue(
            auxiliary_model_ids.isdisjoint(
                record["model_id"] for record in ledger["model_candidates"]
            )
        )
        for receipt_dir in self.fixture.variant_secondary_receipt_dirs:
            receipt = verify_variant_secondary_run_receipt(receipt_dir)
            self.assertFalse(receipt["outcomes_read"])
            self.assertFalse(receipt["metrics_computed"])
            self.assertFalse(receipt["primary_endpoint_scoring_allowed"])
            self.assertNotIn("development_outcome_binding", receipt)
            self.assertNotIn("metrics", receipt)
        for task_id in FIXTURE_TASK_IDS:
            finalist = verify_finalist_development_metric_bundle(
                self.fixture.finalist_metric_bundle_dirs[task_id],
                reverify_sources=True,
            )
            self.assertEqual(len(finalist["selected_run_ids"]), 5)
            self.assertEqual(len(finalist["baseline_run_ids"]), 5)
            self.assertEqual(
                finalist["power_method_id"],
                "empirical_paired_endpoint_bootstrap_v1",
            )

    def test_variant_secondary_candidate_cannot_enter_primary_receipt(self) -> None:
        ledger = verify_selection_candidate_ledger(
            self.fixture.selection_candidate_ledger_dir
        )
        auxiliary = next(
            record
            for record in ledger["variant_secondary_runs"]
            if record["model_id"] == VARIANT_SECONDARY_MODELS[0]
        )
        baseline = next(
            run
            for run in ledger["runs"]
            if run["task_id"] == VARIANT_TASK
            and run["model_id"] == FIXTURE_BASELINE_MODEL
            and run["seed"] == SECONDARY_COMPARATOR_SEED
        )
        with self.assertRaisesRegex(
            TournamentError, "primary_endpoint_scoring_allowed"
        ):
            freeze_scientific_run_receipt(
                candidate_campaign_dir=self.fixture.candidate_dir,
                candidate_execution_attempt_dir=self.fixture.execution_attempt_dirs[
                    auxiliary["run_id"]
                ],
                candidate_prediction_bundle_path=(
                    self.fixture.development_prediction_bundle_paths[
                        auxiliary["run_id"]
                    ]
                ),
                baseline_campaign_dir=self.fixture.candidate_dir,
                baseline_execution_attempt_dir=self.fixture.execution_attempt_dirs[
                    baseline["run_id"]
                ],
                baseline_prediction_bundle_path=(
                    self.fixture.development_prediction_bundle_paths[
                        baseline["run_id"]
                    ]
                ),
                development_outcome_bundle_dir=(
                    self.fixture.development_outcome_bundle_dirs[VARIANT_TASK]
                ),
                output_root=self.root / "invalid_secondary_candidate_receipts",
            )

    def test_variant_secondary_baseline_cannot_enter_primary_receipt(self) -> None:
        ledger = verify_selection_candidate_ledger(
            self.fixture.selection_candidate_ledger_dir
        )
        candidate = next(
            run
            for run in ledger["runs"]
            if run["task_id"] == VARIANT_TASK
            and run["model_id"] == FIXTURE_SELECTED_MODEL
            and run["seed"] == SECONDARY_COMPARATOR_SEED
        )
        auxiliary = next(
            record
            for record in ledger["variant_secondary_runs"]
            if record["model_id"] == VARIANT_SECONDARY_MODELS[0]
        )
        with self.assertRaisesRegex(
            TournamentError, "primary_endpoint_scoring_allowed"
        ):
            freeze_scientific_run_receipt(
                candidate_campaign_dir=self.fixture.candidate_dir,
                candidate_execution_attempt_dir=self.fixture.execution_attempt_dirs[
                    candidate["run_id"]
                ],
                candidate_prediction_bundle_path=(
                    self.fixture.development_prediction_bundle_paths[
                        candidate["run_id"]
                    ]
                ),
                baseline_campaign_dir=self.fixture.candidate_dir,
                baseline_execution_attempt_dir=self.fixture.execution_attempt_dirs[
                    auxiliary["run_id"]
                ],
                baseline_prediction_bundle_path=(
                    self.fixture.development_prediction_bundle_paths[
                        auxiliary["run_id"]
                    ]
                ),
                development_outcome_bundle_dir=(
                    self.fixture.development_outcome_bundle_dirs[VARIANT_TASK]
                ),
                output_root=self.root / "invalid_secondary_baseline_receipts",
            )

    def test_ready_gate_is_nonterminal_and_cannot_authorize_release(self) -> None:
        for path in self.ready_gates.values():
            gate = verify_champion_gate_decision(path, reverify_sources=True)
            self.assertEqual(gate["disposition"], "ready_for_sealed_inference")
            self.assertFalse(gate["terminal"])
        with self.assertRaisesRegex(TournamentError, "nonterminal gate"):
            freeze_terminal_evaluation_authorization(
                selection_lock_dir=self.fixture.selection_lock_dir,
                champion_gate_decision_dirs=self.ready_gates.values(),
                output_root=self.root / "invalid_terminal_authorization",
            )

    def test_opened_outcome_is_terminal_but_never_promoted_without_metric_bundle(self) -> None:
        for path in self.opened_gates.values():
            gate = verify_champion_gate_decision(path, reverify_sources=True)
            self.assertEqual(
                gate["disposition"],
                "terminal_failure_missing_sealed_metric_bundle",
            )
            self.assertTrue(gate["terminal"])
            self.assertFalse(gate["passed"])
            self.assertEqual(gate["failures"], ["sealed_metric_bundle"])
        authorization = freeze_terminal_evaluation_authorization(
            selection_lock_dir=self.fixture.selection_lock_dir,
            champion_gate_decision_dirs=self.opened_gates.values(),
            output_root=self.root / "terminal_authorization",
        )
        verified = verify_terminal_evaluation_authorization(
            authorization, reverify_sources=True
        )
        self.assertTrue(all(item["terminal"] for item in verified["decisions"]))
        self.assertFalse(any(item["promoted"] for item in verified["decisions"]))
        self.assertIn(
            self.fixture.selection_candidate_ledger_dir.as_posix(),
            verified["_protected_source_roots"],
        )

    def test_rederived_metrics_holm_family_and_positive_promotion(self) -> None:
        multiplicity = verify_confirmatory_multiplicity_bundle(
            self.multiplicity, reverify_sources=True
        )
        self.assertEqual(multiplicity["task_ids"], list(FIXTURE_TASK_IDS))
        self.assertEqual(multiplicity["method"], "Holm")
        for task_id in FIXTURE_TASK_IDS:
            metric = verify_sealed_metric_bundle(
                self.sealed_metrics[task_id], reverify_sources=True
            )
            self.assertEqual(metric["n_resamples"], 10_000)
            self.assertEqual(metric["metrics"]["positive_seed_count"], 5)
            evaluator_sources = metric["evaluator_source_bundle"]
            self.assertEqual(
                evaluator_sources["schema_version"],
                "masld-bench-sealed-evaluator-source-bundle-v1",
            )
            self.assertEqual(
                {record["path"] for record in evaluator_sources["sources"]},
                {
                    "src/masld_bench/tournament.py",
                    "src/masld_bench/evaluators/sealed_metrics.py",
                    "src/masld_bench/evaluators/metrics.py",
                    "src/masld_bench/hashing.py",
                },
            )
            self.assertEqual(
                evaluator_sources["bundle_sha256"],
                canonical_sha256(evaluator_sources["sources"]),
            )
            self.assertEqual(
                metric["evaluator_source_sha256"],
                evaluator_sources["bundle_sha256"],
            )
            if task_id == VARIANT_TASK:
                self.assertTrue(
                    metric["variant_primary_capability"]["selected"][
                        "capability"
                    ]["primary_eligible"]
                )
                self.assertTrue(
                    metric["variant_primary_capability"]["baseline"][
                        "capability"
                    ]["primary_eligible"]
                )
                self.assertEqual(
                    metric["variant_secondary_evaluation"][
                        "mandatory_model_ids"
                    ],
                    list(VARIANT_SECONDARY_MODELS),
                )
                # Keyed by real model_id, not by role label: the evaluator
                # indexes this map with decision["selected_model_id"] and
                # decision["baseline_model_id"]
                # (sealed_metrics._variant_metrics).  This class builds its
                # chain from _sealed_fixtures, whose decision names
                # SELECTED_MODEL / BASELINE_MODEL -- not the
                # f"{task_id}_selected" ids the SelectionLock fixture in this
                # module uses.  This assertion had never executed, because
                # setUpClass failed before reaching it.
                self.assertEqual(
                    set(metric["metrics"]["retrieval_auprc_by_model"]),
                    {
                        FIXTURE_SELECTED_MODEL,
                        FIXTURE_BASELINE_MODEL,
                        *VARIANT_SECONDARY_MODELS,
                    },
                )
            else:
                self.assertIsNone(metric["variant_primary_capability"])
                self.assertIsNone(metric["variant_secondary_evaluation"])
            gate = verify_champion_gate_decision(
                self.promoted_gates[task_id], reverify_sources=True
            )
            self.assertEqual(gate["disposition"], "promoted")
            self.assertTrue(gate["passed"])
            self.assertTrue(gate["terminal"])
            self.assertEqual(gate["failures"], [])

        authorization = freeze_terminal_evaluation_authorization(
            selection_lock_dir=self.fixture.selection_lock_dir,
            champion_gate_decision_dirs=self.promoted_gates.values(),
            output_root=self.root / "promoted_terminal_authorization",
        )
        verified = verify_terminal_evaluation_authorization(
            authorization, reverify_sources=True
        )
        self.assertTrue(all(item["promoted"] for item in verified["decisions"]))

    def test_sealed_metric_and_multiplicity_tamper_fail_closed(self) -> None:
        metric_copy = self.root / "tampered_sealed_metric"
        shutil.copytree(self.sealed_metrics[CELL_TASK], metric_copy)
        table = metric_copy / "joined_ensemble.tsv"
        table.chmod(0o640)
        table.write_text(table.read_text(encoding="utf-8") + "tamper\n", encoding="utf-8")
        with self.assertRaises(TournamentError):
            verify_sealed_metric_bundle(metric_copy, reverify_sources=True)

        multiplicity_copy = self.root / "tampered_multiplicity"
        shutil.copytree(self.multiplicity, multiplicity_copy)
        document = multiplicity_copy / "confirmatory_multiplicity_bundle.json"
        document.chmod(0o640)
        original = document.read_text(encoding="utf-8")
        mutated = original.replace('"method":"Holm"', '"method":"BH"')
        self.assertNotEqual(mutated, original)
        document.write_text(
            mutated,
            encoding="utf-8",
        )
        with self.assertRaises(TournamentError):
            verify_confirmatory_multiplicity_bundle(
                multiplicity_copy, reverify_sources=True
            )

    def test_caller_metrics_are_not_a_gate_input(self) -> None:
        task_id = FIXTURE_TASK_IDS[0]
        with self.assertRaises(TypeError):
            freeze_champion_gate_decision(
                selection_lock_dir=self.fixture.selection_lock_dir,
                task_id=task_id,
                output_root=self.root / "forbidden_metrics",
                metrics={"invented": 1.0},
            )

    def test_verifiers_reject_symlink_aliases_and_publications_are_immutable(self) -> None:
        alias = self.root / "finalist_alias"
        alias.symlink_to(
            self.fixture.finalist_metric_bundle_dirs[FIXTURE_TASK_IDS[0]],
            target_is_directory=True,
        )
        with self.assertRaisesRegex(TournamentError, "symlink"):
            verify_finalist_development_metric_bundle(alias)

        task_id = FIXTURE_TASK_IDS[0]
        with self.assertRaisesRegex(TournamentError, "already exists"):
            freeze_champion_gate_decision(
                selection_lock_dir=self.fixture.selection_lock_dir,
                task_id=task_id,
                prediction_commit_dirs=[
                    self.state
                    / "prediction_commits"
                    / commit.prediction_commit_sha256
                    for commit in self.commits[task_id]
                ],
                power_decision_dir=(
                    self.state
                    / "power_decisions"
                    / self.powers[task_id].power_decision_sha256
                ),
                output_root=self.root / "ready_gates" / task_id,
            )


if __name__ == "__main__":
    unittest.main()

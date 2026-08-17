"""Generate immutable, unique run specifications without submitting jobs."""

from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .config import canonical_json_bytes, write_json_exclusive
from .firewall import load_selection_lock


class RunPlanError(RuntimeError):
    pass


def _method(ewc_lambda: float, replay: float) -> str:
    if ewc_lambda == 0 and replay == 0:
        return "fine_tune"
    if ewc_lambda == 0:
        return "replay_only"
    if replay == 0:
        return "ewc_only"
    return "continual_learning"


def screen_settings(config: dict[str, Any], best_lambda: float | None = None):
    values = []
    for ewc_lambda in config["screen"]["lambda_values"]:
        values.append((float(ewc_lambda), 0.20, "lambda_screen"))
    if best_lambda is not None:
        allowed = {float(x) for x in config["screen"]["lambda_values"] if float(x) > 0}
        if float(best_lambda) not in allowed:
            raise RunPlanError("best lambda must be a positive prespecified lambda")
        for replay in config["screen"]["replay_values"]:
            setting = (float(best_lambda), float(replay), "replay_screen")
            if setting[:2] != (float(best_lambda), 0.20):
                values.append(setting)
    pairs = [(x[0], x[1]) for x in values]
    if len(pairs) != len(set(pairs)):
        raise RunPlanError("screen plan contains a duplicate parameter pair")
    expected = 6 if best_lambda is None else 9
    if len(values) != expected:
        raise RunPlanError(f"screen plan must contain {expected} unique settings")
    return [{
        "setting_id": f"lambda_{lam:g}__replay_{replay:g}",
        "ewc_lambda": lam,
        "replay_fraction": replay,
        "method": _method(lam, replay),
        "phase": phase,
        "seed": config["screen"]["seed"],
    } for lam, replay, phase in values]


def create_run_plan(
    config: dict[str, Any], output: str | Path, *, best_lambda: float | None = None,
    selection_lock: str | Path | None = None,
) -> dict[str, Any]:
    settings = screen_settings(config, best_lambda)
    model_kinds = ["all_lineage", *config["lineages"]]
    screen_runs = [
        {**setting, "model_kind": model_kind,
         "resource_class": "large" if model_kind in {"all_lineage", "Hepatocytes"} else "small"}
        for setting in settings for model_kind in model_kinds
    ]
    plan: dict[str, Any] = {
        "schema_version": "masld-cl-run-plan-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_sha256": config["_config_sha256"],
        "screen_settings": settings,
        "screen_runs": screen_runs,
        "baseline_methods": [
            "incumbent_harmony", "de_novo", "architecture_surgery", "fine_tune",
            "replay_only", "ewc_only", "continual_learning",
        ],
        "selection_required_before": [
            "confirmation", "outcome_preservation", "held_study", "secondary_stress",
            "acquisition_order", "bi_guided_replay",
        ],
    }
    if selection_lock is not None:
        selection = load_selection_lock(selection_lock, config)
        selected = selection["selected"]
        runner = selection.get("pareto_runner_up")
        settings_by_role = {
            "selected": (selected["ewc_lambda"], selected["replay_fraction"]),
            "paper": (100.0, 0.20),
        }
        if runner is not None:
            settings_by_role["pareto_runner_up"] = (
                runner["ewc_lambda"], runner["replay_fraction"]
            )
        settings_by_pair: dict[tuple[float, float], list[str]] = {}
        for role, pair in settings_by_role.items():
            settings_by_pair.setdefault((float(pair[0]), float(pair[1])), []).append(role)
        confirmation = []
        for (ewc_lambda, replay), roles in settings_by_pair.items():
            for seed in config["screen"]["confirmation_seeds"]:
                for model_kind in model_kinds:
                    confirmation.append({
                        "roles": sorted(roles), "setting_id": f"lambda_{ewc_lambda:g}__replay_{replay:g}",
                        "ewc_lambda": ewc_lambda, "replay_fraction": replay,
                        "method": _method(ewc_lambda, replay), "seed": seed,
                        "model_kind": model_kind,
                        "resource_class": "large" if model_kind in {"all_lineage", "Hepatocytes"} else "small",
                    })
        plan.update({
            "selection_lock_sha256": selection["lock_sha256"],
            "confirmation_runs": confirmation,
            "held_study_tests": [
                {
                    "held_out": study,
                    "adaptation_datasets": [
                        value for value in config["acquisition_orders"]["accession"]
                        if value != study
                    ],
                    "training_labels": "query_cell_type_and_stage_hidden",
                    "held_out_datasets": [study],
                    "control_fisher_datasets": [
                        value for value in config["evaluation"]["powered_query_studies"]
                        if value != study
                    ],
                    "control_fisher": "powered_query_controls_excluding_held_study",
                }
                for study in config["evaluation"]["powered_query_studies"]
            ],
            "secondary_stress_tests": [
                {
                    "stress_dataset": study,
                    "adaptation_datasets": [
                        *config["evaluation"]["powered_query_studies"], study,
                    ],
                    "training_labels": "query_cell_type_and_stage_hidden",
                    "held_out_datasets": [],
                    "control_fisher_datasets": config["evaluation"]["powered_query_studies"],
                    "control_fisher": "powered_primary_query_controls_only",
                }
                for study in config["evaluation"]["secondary_stress_studies"]
            ],
            "mapping_only": config["evaluation"]["mapping_only_studies"],
            "acquisition_orders": config["acquisition_orders"],
            "bi_guided_replay": {
                "prerequisite": "random_replay_promotion_decision_passed",
                "required_unlock_schema": "masld-cl-bi-replay-unlock-v1",
                "sensitivity_only": True,
                "n_augmentations": 200,
                "mask_fraction": 0.50,
                "selection_modes": ["bottom", "top", "step"],
                "buffer_size": "exact_selected_random_replay_size",
            },
        })
    plan["plan_sha256"] = hashlib.sha256(canonical_json_bytes(plan)).hexdigest()
    write_json_exclusive(output, plan)
    return plan

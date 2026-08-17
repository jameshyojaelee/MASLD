"""True leave-query-study-out control selection and held-study evaluation for V9."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .benchmark import _bundle_centroids, _paired_improvement
from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .embedding import load_embedding
from .firewall import load_retargeted_bridge_selection_lock, validate_program_firewall
from .metrics import bh_adjust


def _grid(selection: dict[str, Any]) -> dict[tuple[float, str], dict[str, Any]]:
    result = {}
    for source in selection["control_results"]:
        with Path(source["path"]).open() as handle:
            row = json.load(handle)
        result[(float(row["global_reference_weight"]), row["model_kind"])] = row
    return result


def _study_result(config, row, study: str, cache, seed_offset: int):
    def bundle(path):
        path = str(Path(path).resolve())
        if path not in cache:
            cache[path] = load_embedding(path)
        return cache[path]

    candidate = bundle(row["sources"]["candidate_embedding"])
    harmony = bundle(row["sources"]["harmony_embedding"])
    architecture = bundle(row["sources"]["architecture_embedding"])
    lineage = None if row["model_kind"] == "all_lineage" else row["model_kind"]
    try:
        harmony_result = _paired_improvement(
            _bundle_centroids(candidate, lineage, control_study=study),
            _bundle_centroids(harmony, lineage, control_study=study),
            config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + seed_offset,
        )
        architecture_result = _paired_improvement(
            _bundle_centroids(candidate, lineage, control_study=study),
            _bundle_centroids(architecture, lineage, control_study=study),
            config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + seed_offset + 1,
        )
        protocol_harmony = _paired_improvement(
            _bundle_centroids(candidate, lineage, control_study=study, preparation="unsorted"),
            _bundle_centroids(harmony, lineage, control_study=study, preparation="unsorted"),
            config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + seed_offset + 2,
        )
        protocol_architecture = _paired_improvement(
            _bundle_centroids(candidate, lineage, control_study=study, preparation="unsorted"),
            _bundle_centroids(architecture, lineage, control_study=study, preparation="unsorted"),
            config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + seed_offset + 3,
        )
    except ContractError as error:
        if "at least three donors" not in str(error):
            raise
        return {
            "model_kind": row["model_kind"], "study": study, "testable": False,
            "reason": str(error),
            "prespecified_underpowered": f"{study}|{row['model_kind']}|control_alignment"
            in config["evaluation"]["prespecified_underpowered"],
        }
    return {
        "model_kind": row["model_kind"], "study": study, "testable": True,
        "harmony": harmony_result, "architecture_surgery": architecture_result,
        "protocol_harmony": protocol_harmony,
        "protocol_architecture_surgery": protocol_architecture,
    }


def _summarize(config, rows):
    all_lineage = next((row for row in rows if row["model_kind"] == "all_lineage"), None)
    if all_lineage is None or not all_lineage["testable"]:
        return {"eligible": False, "reason": "all_lineage_untestable"}
    lineage_rows = {row["model_kind"]: row for row in rows if row["model_kind"] in config["lineages"]}
    if set(lineage_rows) != set(config["lineages"]):
        raise ContractError("held-study summary lacks the five frozen lineages")
    p_harmony = np.asarray([
        lineage_rows[lineage]["harmony"]["p_one_sided"] if lineage_rows[lineage]["testable"] else 1.0
        for lineage in config["lineages"]
    ])
    p_architecture = np.asarray([
        lineage_rows[lineage]["architecture_surgery"]["p_one_sided"]
        if lineage_rows[lineage]["testable"] else 1.0
        for lineage in config["lineages"]
    ])
    q_harmony, q_architecture = bh_adjust(p_harmony), bh_adjust(p_architecture)
    threshold = config["gates"]["minimum_control_alignment_improvement"]
    improved = []
    for index, lineage in enumerate(config["lineages"]):
        row = lineage_rows[lineage]
        if not row["testable"]:
            continue
        if (
            row["harmony"]["improvement"] >= threshold
            and row["harmony"]["ci_low"] > 0 and q_harmony[index] <= 0.05
            and row["architecture_surgery"]["improvement"] >= threshold
            and row["architecture_surgery"]["ci_low"] > 0 and q_architecture[index] <= 0.05
        ):
            improved.append(lineage)
    protocol = all(
        not row["testable"] or (
            row["protocol_harmony"]["improvement"] >= -config["gates"]["maximum_stratum_worsening"]
            and row["protocol_architecture_surgery"]["improvement"] >= -config["gates"]["maximum_stratum_worsening"]
        ) for row in rows
    )
    all_pass = (
        all_lineage["harmony"]["improvement"] >= threshold
        and all_lineage["harmony"]["ci_low"] > 0
        and all_lineage["architecture_surgery"]["improvement"] >= threshold
        and all_lineage["architecture_surgery"]["ci_low"] > 0
    )
    mean_shift = float(np.mean([
        row["harmony"]["candidate_shift"] for row in rows if row["testable"]
    ]))
    return {
        "eligible": bool(all_pass and protocol and len(improved) >= config["gates"]["minimum_improved_lineages"]),
        "all_lineage_pass": bool(all_pass), "protocol_no_material_worsening": bool(protocol),
        "improved_lineages": improved, "improved_lineage_count": len(improved),
        "mean_candidate_shift": mean_shift,
        "bh_adjustment": {
            lineage: {"harmony": float(q_harmony[index]), "architecture_surgery": float(q_architecture[index])}
            for index, lineage in enumerate(config["lineages"])
        },
    }


def evaluate_leave_query_study_out(
    config: dict[str, Any], selection_lock: str | Path, output: str | Path,
) -> dict[str, Any]:
    selection_path = Path(selection_lock).resolve()
    selection = load_retargeted_bridge_selection_lock(selection_path, config)
    validate_program_firewall(config)
    grid = _grid(selection)
    weights = sorted({key[0] for key in grid})
    kinds = ("all_lineage", *config["lineages"])
    cache, held_results = {}, {}
    studies = list(config["evaluation"]["powered_query_studies"])
    if len(studies) != 2:
        raise ContractError("V9 leave-query-study-out requires exactly two powered query studies")
    for held in studies:
        training = next(study for study in studies if study != held)
        training_weights = {}
        training_rows_by_weight = {}
        for weight_index, weight in enumerate(weights):
            rows = [
                _study_result(
                    config, grid[(weight, kind)], training, cache,
                    8101 + 100 * weight_index + 10 * kind_index,
                ) for kind_index, kind in enumerate(kinds)
            ]
            training_rows_by_weight[str(weight)] = rows
            training_weights[str(weight)] = _summarize(config, rows)
        eligible = [weight for weight in weights if training_weights[str(weight)]["eligible"]]
        if not eligible:
            raise ContractError(f"no control-only weight survives when training on {training}")
        selected_weight = min(
            eligible, key=lambda weight: (training_weights[str(weight)]["mean_candidate_shift"], weight)
        )
        held_rows = [
            _study_result(
                config, grid[(selected_weight, kind)], held, cache,
                9101 + 10 * kind_index,
            ) for kind_index, kind in enumerate(kinds)
        ]
        held_summary = _summarize(config, held_rows)
        held_results[held] = {
            "training_study": training, "selected_global_reference_weight": selected_weight,
            "training_weight_summaries": training_weights,
            "training_control_rows": training_rows_by_weight,
            "held_control_rows": held_rows, "held_summary": held_summary,
            "query_cell_type_and_stage_labels_hidden": True,
            "held_control_indicator_used_only_for_control_calibration": True,
            "disease_geometry_pass_inherited_from_uniform_within_study_translation": True,
        }
    passed = bool(all(value["held_summary"]["eligible"] for value in held_results.values()))
    result = {
        "schema_version": "masld-cl-orthogonal-bridge-leave-query-study-out-v9",
        "config_sha256": config["_config_sha256"],
        "selection_lock": str(selection_path), "selection_lock_file_sha256": sha256_path(selection_path),
        "selection_lock_sha256": selection["lock_sha256"],
        "held_studies": held_results, "leave_query_study_out_pass": passed,
        "case_stage_program_hero_gene_and_cas13_not_used": True,
    }
    write_json_exclusive(output, result)
    return result

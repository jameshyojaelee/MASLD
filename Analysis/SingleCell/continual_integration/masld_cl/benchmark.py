"""Post-lock baseline comparison for identity and control-alignment gates."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .control_evaluation import _neighbor_jaccard_loss, _reference_f1_change_ci
from .embedding import load_embedding, matched_rows
from .execution import execution_owned_files, verify_execution_record
from .firewall import load_selection_lock, validate_program_firewall
from .metrics import (
    bh_adjust, donor_centroids, macro_f1, normalized_shift_control,
    positive_class_f1,
)


def _bundle_centroids(bundle, lineage=None, control_study=None, preparation=None):
    _, latent, cells = bundle
    mask = cells["analysis_eligible"].to_numpy(dtype=bool, copy=True)
    mask &= cells["strict_reference"].to_numpy(dtype=bool) | (
        cells["primary_query"].to_numpy(dtype=bool)
        & cells["query_control"].to_numpy(dtype=bool)
    )
    if lineage is not None:
        mask &= cells["audit_cell_type"].astype(str).to_numpy() == lineage
    if preparation is not None:
        mask &= cells["preparation"].astype(str).to_numpy() == preparation
    if control_study is not None:
        is_reference = cells["strict_reference"].to_numpy(dtype=bool)
        in_study = cells["dataset"].astype(str).to_numpy() == control_study
        mask &= is_reference | in_study
    selected = cells.loc[mask].reset_index(drop=True)
    donors, values = donor_centroids(np.asarray(latent[mask], float), selected["donor_id"].astype(str))
    metadata = selected.drop_duplicates("donor_id").set_index("donor_id").loc[donors]
    reference = metadata["strict_reference"].to_numpy(dtype=bool)
    controls = metadata["primary_query"].to_numpy(dtype=bool) & metadata["query_control"].to_numpy(dtype=bool)
    return donors.astype(str), values, metadata["dataset"].astype(str).to_numpy(), reference, controls


def _shift_case_supporting(bundle, lineage=None) -> float:
    _, latent, cells = bundle
    mask = cells["analysis_eligible"].to_numpy(dtype=bool)
    if lineage is not None:
        mask &= cells["audit_cell_type"].astype(str).to_numpy() == lineage
    reference = cells["strict_reference"].to_numpy(dtype=bool)
    cases = (
        cells["primary_query"].to_numpy(dtype=bool)
        & ~cells["query_control"].to_numpy(dtype=bool)
    )
    mask &= reference | cases
    selected = cells.loc[mask].reset_index(drop=True)
    donors, values = donor_centroids(
        np.asarray(latent[mask], dtype=float), selected["donor_id"].astype(str)
    )
    metadata = selected.drop_duplicates("donor_id").set_index("donor_id").loc[donors]
    is_reference = metadata["strict_reference"].to_numpy(dtype=bool)
    is_case = (
        metadata["primary_query"].to_numpy(dtype=bool)
        & ~metadata["query_control"].to_numpy(dtype=bool)
    )
    if is_reference.sum() < 3 or is_case.sum() < 3:
        raise ContractError("shift-case requires at least three donors per group")
    return normalized_shift_control(values[is_reference], values[is_case])


def _align_centroids(candidate, baseline):
    c_donor, c_value, c_study, c_ref, c_ctrl = candidate
    b_donor, b_value, b_study, b_ref, b_ctrl = baseline
    if set(c_donor) != set(b_donor):
        raise ContractError("candidate and baseline donor rosters differ")
    index = {donor: i for i, donor in enumerate(b_donor)}
    order = np.asarray([index[x] for x in c_donor], dtype=np.int64)
    if not (
        np.array_equal(c_study, b_study[order])
        and np.array_equal(c_ref, b_ref[order])
        and np.array_equal(c_ctrl, b_ctrl[order])
    ):
        raise ContractError("candidate and baseline donor roles differ")
    return c_value, b_value[order], c_study, c_ref, c_ctrl


def _stratified_draw(studies: np.ndarray, rng) -> np.ndarray:
    output = []
    for study in sorted(set(studies)):
        indices = np.flatnonzero(studies == study)
        output.extend(rng.choice(indices, len(indices), replace=True))
    return np.asarray(output, dtype=np.int64)


def _paired_improvement(candidate, baseline, replicates: int, seed: int):
    c, b, studies, reference, controls = _align_centroids(candidate, baseline)
    if reference.sum() < 3 or controls.sum() < 3:
        raise ContractError("control shift requires at least three donors per group")
    candidate_shift = normalized_shift_control(c[reference], c[controls])
    baseline_shift = normalized_shift_control(b[reference], b[controls])
    rng = np.random.default_rng(seed)
    draws = []
    r_study, c_study = studies[reference], studies[controls]
    c_ref, c_ctrl = c[reference], c[controls]
    b_ref, b_ctrl = b[reference], b[controls]
    for _ in range(replicates):
        r = _stratified_draw(r_study, rng)
        q = _stratified_draw(c_study, rng)
        try:
            cand = normalized_shift_control(c_ref[r], c_ctrl[q])
            base = normalized_shift_control(b_ref[r], b_ctrl[q])
        except ValueError:
            continue
        if base > 0:
            draws.append((base - cand) / base)
    if len(draws) < max(100, int(0.9 * replicates)):
        raise ContractError("too many degenerate paired control bootstraps")
    draws = np.asarray(draws)
    return {
        "candidate_shift": candidate_shift,
        "baseline_shift": baseline_shift,
        "improvement": (baseline_shift - candidate_shift) / baseline_shift,
        "ci_low": float(np.quantile(draws, 0.025)),
        "ci_high": float(np.quantile(draws, 0.975)),
        "p_one_sided": float((1 + np.sum(draws <= 0)) / (len(draws) + 1)),
        "valid_replicates": len(draws),
    }


def _mean_donor_f1(cells) -> float:
    values = [
        macro_f1(group["audit_cell_type"], group["predicted_cell_type"])
        for _, group in cells.groupby("donor_id", sort=True)
    ]
    return float(np.mean(values))


NON_CL_METHODS = {
    "de_novo", "architecture_surgery", "fine_tune", "replay_only", "ewc_only",
}


def _validate_harmony_info(info: dict[str, Any], config: dict[str, Any]) -> None:
    if (
        info.get("config_sha256") != config["_config_sha256"]
        or info.get("method") != "incumbent_scalesc_pca_harmony"
        or info.get("model_kind") != "all_lineage"
    ):
        raise ContractError("Harmony comparator does not match the current config and incumbent method")


def _improved_lineage_count(
    lineages: list[str], results: dict[str, dict[str, Any]],
    p_harmony: np.ndarray, p_architecture: np.ndarray, threshold: float,
) -> int:
    return sum(
        results[lineage]["harmony"]["improvement"] >= threshold
        and results[lineage]["harmony"]["ci_low"] > 0
        and results[lineage]["architecture_surgery"]["improvement"] >= threshold
        and results[lineage]["architecture_surgery"]["ci_low"] > 0
        and p_harmony[index] <= 0.05
        and p_architecture[index] <= 0.05
        for index, lineage in enumerate(lineages)
    )


def _label_scores(bundle, lineages: list[str]) -> tuple[float, dict[str, float]]:
    _, _, cells = bundle
    query = cells["analysis_eligible"].to_numpy(dtype=bool) & ~cells["strict_reference"].to_numpy(dtype=bool)
    query_cells = cells.loc[query].reset_index(drop=True)
    if query_cells.empty:
        raise ContractError("label concordance has no eligible query cells")
    overall = _mean_donor_f1(query_cells)
    lineage_scores: dict[str, float] = {}
    for lineage in lineages:
        values = []
        for _, donor_cells in query_cells.groupby("donor_id", sort=True):
            truth = donor_cells["audit_cell_type"].astype(str).to_numpy() == lineage
            if not truth.any():
                continue
            predicted = donor_cells["predicted_cell_type"].astype(str).to_numpy() == lineage
            values.append(positive_class_f1(truth, predicted))
        if not values:
            raise ContractError(f"major lineage is not evaluable for label concordance: {lineage}")
        lineage_scores[lineage] = float(np.mean(values))
    return overall, lineage_scores


def _label_concordance(candidate_bundle, baseline_bundle, lineages: list[str]):
    _, _, candidate = candidate_bundle
    _, _, baseline = baseline_bundle
    left, right = matched_rows(candidate, baseline)
    if len(left) != len(candidate):
        raise ContractError("non-CL label comparator lacks candidate cells")
    candidate = candidate.iloc[left].reset_index(drop=True)
    baseline = baseline.iloc[right].reset_index(drop=True)
    if not np.array_equal(
        candidate["audit_cell_type"].astype(str).to_numpy(),
        baseline["audit_cell_type"].astype(str).to_numpy(),
    ):
        raise ContractError("non-CL comparator has different frozen audit labels")
    candidate_scores = _label_scores((candidate_bundle[0], candidate_bundle[1], candidate), lineages)
    baseline_scores = _label_scores((baseline_bundle[0], baseline_bundle[1], baseline), lineages)
    changes = [
        candidate_scores[1][lineage] - baseline_scores[1][lineage]
        for lineage in lineages
    ]
    return candidate_scores[0] - baseline_scores[0], min(changes)


def evaluate_benchmarks(
    config: dict[str, Any], selection_lock: str | Path,
    reference_embedding: str | Path, harmony_embedding: str | Path,
    non_cl_embeddings: dict[str, str | Path],
    candidate_embeddings: dict[str, str | Path],
    architecture_embeddings: dict[str, str | Path],
    execution_records: list[str | Path], output: str | Path,
) -> list[dict[str, Any]]:
    selection = load_selection_lock(selection_lock, config)
    validate_program_firewall(config)
    harmony = load_embedding(harmony_embedding)
    reference = load_embedding(reference_embedding)
    _validate_harmony_info(harmony[0], config)
    expected = {"all_lineage", *config["lineages"]}
    if set(candidate_embeddings) != expected or set(architecture_embeddings) != expected:
        raise ContractError(f"benchmark requires candidate and architecture bundles for {sorted(expected)}")
    if set(non_cl_embeddings) != NON_CL_METHODS:
        raise ContractError(f"benchmark requires exact non-CL roles: {sorted(NON_CL_METHODS)}")
    pipeline_root = Path(config["_config_path"]).resolve().parent
    record_paths = [str(Path(value).resolve()) for value in execution_records]
    if not record_paths or len(record_paths) != len(set(record_paths)):
        raise ContractError("benchmark requires unique GPU execution records")
    verified_records = {
        path: verify_execution_record(path, pipeline_root, config["_config_sha256"])
        for path in record_paths
    }
    embedding_sources: dict[str, Any] = {}

    def owners_for(*paths: str | Path) -> list[str]:
        resolved = {str(Path(path).resolve()) for path in paths}
        owners = sorted(
            record_path for record_path, record in verified_records.items()
            if resolved.issubset(execution_owned_files(record))
        )
        if len(owners) != 1:
            raise ContractError(
                f"benchmark artifacts require exactly one GPU owner; paths={sorted(resolved)}"
            )
        return owners

    def register_source(
        key: str, embedding_path: Path, run_path: Path | None, *, gpu_owned: bool,
    ) -> None:
        source = {
            "embedding_manifest": str(embedding_path),
            "embedding_manifest_sha256": sha256_path(embedding_path),
        }
        if run_path is not None:
            source.update({
                "run_manifest": str(run_path),
                "run_manifest_sha256": sha256_path(run_path),
            })
        if gpu_owned:
            owners = owners_for(embedding_path, *([] if run_path is None else [run_path]))
            source["execution_records"] = [
                {"path": path, "sha256": sha256_path(path)} for path in owners
            ]
        embedding_sources[key] = source

    selected = selection["selected"]
    if float(selected["ewc_lambda"]) <= 0 or float(selected["replay_fraction"]) <= 0:
        raise ContractError("promotion benchmark requires a selected replay-plus-EWC setting")
    candidate_context = None
    for role, mapping in (
        ("candidate", candidate_embeddings), ("architecture", architecture_embeddings)
    ):
        for model_kind, value in sorted(mapping.items()):
            path = Path(value).resolve()
            run_name = "update_manifest.json"
            run_path = path.parent / run_name
            if not run_path.is_file():
                raise ContractError(f"benchmark embedding lacks sibling {run_name}: {path}")
            with run_path.open() as handle:
                run = json.load(handle)
            with path.open() as handle:
                embedding_info = json.load(handle)
            if (
                run.get("config_sha256") != config["_config_sha256"]
                or run.get("model_kind") != model_kind
                or run.get("selection_lock_sha256") != selection["lock_sha256"]
                or run.get("sensitivity_only") is True
                or run.get("embedding") != embedding_info
            ):
                raise ContractError(f"benchmark run provenance mismatch: {run_path}")
            if role == "candidate" and (
                run.get("method") != "continual_learning"
                or run.get("replay_mode") != "random"
                or float(run.get("ewc_lambda")) != float(selected["ewc_lambda"])
                or float(run.get("replay_fraction")) != float(selected["replay_fraction"])
            ):
                raise ContractError("benchmark candidate is not the selected random-replay CL fit")
            if role == "architecture" and (
                run.get("method") != "architecture_surgery"
                or float(run.get("ewc_lambda")) != 0
                or float(run.get("replay_fraction")) != 0
            ):
                raise ContractError("architecture comparator is not architecture surgery")
            context = {
                key: run.get(key) for key in (
                    "seed", "production", "query_datasets", "held_out_datasets",
                )
            }
            if role == "candidate":
                if candidate_context is None:
                    candidate_context = context
                elif context != candidate_context:
                    raise ContractError("six candidate models do not share one training context")
            elif context != candidate_context:
                raise ContractError("architecture comparator has a different training context")
            register_source(f"{role}|{model_kind}", path, run_path.resolve(), gpu_owned=True)
    candidate_context.update({
        "ewc_lambda": float(selected["ewc_lambda"]),
        "replay_fraction": float(selected["replay_fraction"]),
        "method": "continual_learning",
        "replay_mode": "random",
    })
    reference_path = Path(reference_embedding).resolve()
    reference_run_path = reference_path.parent / "reference_manifest.json"
    if not reference_run_path.is_file():
        raise ContractError("reference embedding lacks sibling reference manifest")
    with reference_run_path.open() as handle:
        reference_run = json.load(handle)
    if (
        reference_run.get("schema_version") != "masld-cl-reference-v1"
        or reference_run.get("model_kind") != "all_lineage"
        or reference_run.get("config_sha256") != config["_config_sha256"]
        or reference_run.get("embedding") != reference[0]
    ):
        raise ContractError("benchmark reference provenance mismatch")
    register_source("reference", reference_path, reference_run_path.resolve(), gpu_owned=True)
    harmony_path = Path(harmony_embedding).resolve()
    register_source("harmony", harmony_path, None, gpu_owned=False)

    non_cl_bundles: dict[str, Any] = {}
    non_cl_scores: dict[str, float] = {}
    for method, value in sorted(non_cl_embeddings.items()):
        path = Path(value).resolve()
        if method == "architecture_surgery" and path != Path(
            architecture_embeddings["all_lineage"]
        ).resolve():
            raise ContractError("non-CL architecture comparator must reuse the all-lineage architecture bundle")
        candidates = [path.parent / "update_manifest.json", path.parent / "denovo_manifest.json"]
        run_paths = [candidate for candidate in candidates if candidate.is_file()]
        if len(run_paths) != 1:
            raise ContractError(f"non-CL embedding lacks one unambiguous training manifest: {path}")
        run_path = run_paths[0].resolve()
        with run_path.open() as handle:
            run = json.load(handle)
        with path.open() as handle:
            embedding_info = json.load(handle)
        expected_schema = "masld-cl-denovo-v1" if method == "de_novo" else "masld-cl-update-v1"
        if (
            run.get("schema_version") != expected_schema
            or run.get("config_sha256") != config["_config_sha256"]
            or run.get("model_kind") != "all_lineage"
            or run.get("selection_lock_sha256") != selection["lock_sha256"]
            or run.get("sensitivity_only") is True
            or run.get("seed") != candidate_context["seed"]
            or run.get("production") != candidate_context["production"]
            or run.get("query_datasets") != candidate_context["query_datasets"]
            or run.get("embedding") != embedding_info
        ):
            raise ContractError(f"non-CL comparator has different provenance or context: {method}")
        if method != "de_novo" and (
            run.get("method") != method
            or run.get("held_out_datasets") != candidate_context["held_out_datasets"]
        ):
            raise ContractError(f"non-CL comparator method/context mismatch: {method}")
        bundle = load_embedding(path)
        score, _ = _label_scores(bundle, config["lineages"])
        non_cl_bundles[method] = bundle
        non_cl_scores[method] = score
        register_source(f"non_cl|{method}", path, run_path, gpu_owned=True)
    best_non_cl_method = sorted(non_cl_scores, key=lambda key: (-non_cl_scores[key], key))[0]
    best_non_cl = non_cl_bundles[best_non_cl_method]
    used_record_paths = {
        source["path"]
        for value in embedding_sources.values()
        for source in value.get("execution_records", [])
    }
    if used_record_paths != set(record_paths):
        raise ContractError("supplied benchmark GPU records differ from exact artifact owners")
    details: dict[str, Any] = {
        "alignment": {}, "alignment_strata": [], "untestable_strata": [],
    }
    improvement_by_lineage: dict[str, dict[str, Any]] = {}
    candidate_all = load_embedding(candidate_embeddings["all_lineage"])
    reference_cells = reference[2].loc[reference[2]["strict_reference"]].reset_index(drop=True)
    candidate_reference = candidate_all[2].loc[candidate_all[2]["strict_reference"]].reset_index(drop=True)
    left, right = matched_rows(reference_cells, candidate_reference)
    f1 = _reference_f1_change_ci(
        reference_cells.iloc[left].reset_index(drop=True),
        candidate_reference.iloc[right].reset_index(drop=True),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"],
    )
    candidate_reference_indices = np.flatnonzero(candidate_all[2]["strict_reference"].to_numpy())[right]
    jaccard_loss = _neighbor_jaccard_loss(
        np.asarray(reference[1])[left], np.asarray(candidate_all[1])[candidate_reference_indices],
        reference_cells.iloc[left].reset_index(drop=True),
        k=config["evaluation"]["reference_neighborhood_k"],
        maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
        seed=config["screen"]["seed"],
    )
    query_f1, major_f1 = _label_concordance(candidate_all, best_non_cl, config["lineages"])
    all_results = {}
    for model_kind in sorted(expected):
        candidate = candidate_all if model_kind == "all_lineage" else load_embedding(candidate_embeddings[model_kind])
        architecture = load_embedding(architecture_embeddings[model_kind])
        lineage = None if model_kind == "all_lineage" else model_kind
        candidate_centroids = _bundle_centroids(candidate, lineage)
        harmony_centroids = _bundle_centroids(harmony, lineage)
        architecture_centroids = _bundle_centroids(architecture, lineage)
        harmony_result = _paired_improvement(
            candidate_centroids, harmony_centroids,
            config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 101,
        )
        architecture_result = _paired_improvement(
            candidate_centroids, architecture_centroids,
            config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 201,
        )
        all_results[model_kind] = {"harmony": harmony_result, "architecture_surgery": architecture_result}
        case_shift = _shift_case_supporting(candidate, lineage)
        all_results[model_kind]["shift_case_supporting"] = case_shift
        all_results[model_kind]["shift_case_to_control_ratio_supporting"] = (
            case_shift / harmony_result["candidate_shift"]
        )
        details["alignment_strata"].append({
            "model_kind": model_kind, "stratum": "pooled_primary_query",
            "worsening_vs_harmony": -harmony_result["improvement"],
            "worsening_vs_architecture_surgery": -architecture_result["improvement"],
        })
        if model_kind != "all_lineage":
            improvement_by_lineage[model_kind] = all_results[model_kind]
        for preparation in sorted(set(candidate[2]["preparation"].astype(str))):
            try:
                c = _bundle_centroids(candidate, lineage, preparation=preparation)
                h = _bundle_centroids(harmony, lineage, preparation=preparation)
                a = _bundle_centroids(architecture, lineage, preparation=preparation)
                h_result = _paired_improvement(c, h, config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 301)
                a_result = _paired_improvement(c, a, config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 401)
            except ContractError as error:
                if "at least three donors" not in str(error):
                    raise
                details["untestable_strata"].append({
                    "model_kind": model_kind,
                    "preparation": preparation,
                    "reason": str(error),
                    "prespecified_minimum_rule": True,
                })
                continue
            details["alignment_strata"].append({
                "model_kind": model_kind, "stratum": f"preparation={preparation}",
                "worsening_vs_harmony": -h_result["improvement"],
                "worsening_vs_architecture_surgery": -a_result["improvement"],
            })
        for study in config["evaluation"]["powered_query_studies"]:
            try:
                c = _bundle_centroids(candidate, lineage, control_study=study)
                h = _bundle_centroids(harmony, lineage, control_study=study)
                a = _bundle_centroids(architecture, lineage, control_study=study)
                h_result = _paired_improvement(c, h, config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 501)
                a_result = _paired_improvement(c, a, config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 601)
                details["alignment_strata"].append({
                    "model_kind": model_kind, "stratum": f"control_study={study}",
                    "worsening_vs_harmony": -h_result["improvement"],
                    "worsening_vs_architecture_surgery": -a_result["improvement"],
                })
            except ContractError as error:
                key = f"{study}|{model_kind}|control_alignment"
                prespecified = key in config["evaluation"]["prespecified_underpowered"]
                details["untestable_strata"].append({
                    "model_kind": model_kind, "study": study,
                    "reason": str(error), "prespecified_underpowered": prespecified,
                })
                if not prespecified:
                    raise
    details["alignment"] = all_results
    threshold = config["gates"]["minimum_control_alignment_improvement"]
    p_harmony = bh_adjust([improvement_by_lineage[x]["harmony"]["p_one_sided"] for x in config["lineages"]])
    p_architecture = bh_adjust([improvement_by_lineage[x]["architecture_surgery"]["p_one_sided"] for x in config["lineages"]])
    details["bh_adjustment"] = {
        lineage: {"harmony": float(p_harmony[i]), "architecture_surgery": float(p_architecture[i])}
        for i, lineage in enumerate(config["lineages"])
    }
    improved = _improved_lineage_count(
        config["lineages"], improvement_by_lineage,
        p_harmony, p_architecture, threshold,
    )
    maximum_worsening = max(
        [0.0, *[
            max(x["worsening_vs_harmony"], x["worsening_vs_architecture_surgery"])
            for x in details["alignment_strata"]
        ]]
    )
    all_lineage = all_results["all_lineage"]
    rows = [
        {"metric": "reference_macro_f1_change_ci_low", "scope": "all_lineage", "value": f1["ci_low"]},
        {"metric": "reference_neighborhood_jaccard_loss", "scope": "all_lineage", "value": jaccard_loss},
        {"metric": "query_macro_f1_change", "scope": "all_lineage", "value": query_f1},
        {"metric": "major_lineage_f1_change", "scope": "worst_lineage", "value": major_f1},
        {"metric": "alignment_improvement_vs_harmony", "scope": "all_lineage", "value": all_lineage["harmony"]["improvement"]},
        {"metric": "alignment_improvement_vs_harmony_ci_low", "scope": "all_lineage", "value": all_lineage["harmony"]["ci_low"]},
        {"metric": "alignment_improvement_vs_architecture_surgery", "scope": "all_lineage", "value": all_lineage["architecture_surgery"]["improvement"]},
        {"metric": "alignment_improvement_vs_architecture_surgery_ci_low", "scope": "all_lineage", "value": all_lineage["architecture_surgery"]["ci_low"]},
        {"metric": "improved_lineage_count", "scope": "five_frozen_lineages", "value": improved},
        {"metric": "all_lineage_alignment_pass", "scope": "all_lineage", "value": int(
            all_lineage["harmony"]["improvement"] >= threshold
            and all_lineage["harmony"]["ci_low"] > 0
            and all_lineage["architecture_surgery"]["improvement"] >= threshold
            and all_lineage["architecture_surgery"]["ci_low"] > 0
        )},
        {"metric": "maximum_protocol_stratum_worsening", "scope": "all_evaluable", "value": maximum_worsening},
        {"metric": "shift_case_supporting", "scope": "all_lineage", "value": all_lineage["shift_case_supporting"]},
        {"metric": "shift_case_to_control_ratio_supporting", "scope": "all_lineage", "value": all_lineage["shift_case_to_control_ratio_supporting"]},
    ]
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "scope", "value"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_json_exclusive(output.with_suffix(output.suffix + ".details.json"), {
        "schema_version": "masld-cl-benchmark-v1",
        "config_sha256": config["_config_sha256"],
        "selection_lock_sha256": selection["lock_sha256"],
        "metrics_realpath": str(output.resolve()),
        "metrics_sha256": sha256_path(output),
        "candidate_context": candidate_context,
        "non_cl_label_scores": non_cl_scores,
        "best_non_cl_method": best_non_cl_method,
        "embedding_sources": embedding_sources,
        "execution_records": [
            {"path": path, "sha256": sha256_path(path)} for path in record_paths
        ],
        **details,
    })
    return rows

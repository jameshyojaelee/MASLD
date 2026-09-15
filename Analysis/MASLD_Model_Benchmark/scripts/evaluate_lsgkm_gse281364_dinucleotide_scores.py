#!/usr/bin/env python3
"""Evaluate fixed LS-GKM raw rank scores on exposed-development GSE281364."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.stats import rankdata, spearmanr

from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes


SEEDS = (1103, 2909, 4721, 6673, 8111)
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
CANDIDATES = (
    ("available_simple_controls", "allele_identity_ridge"),
    ("lsgkm_exact_dinucleotide_null", "direct_gkmsvm"),
    ("lsgkm_exact_dinucleotide_null", "deltasvm"),
)
REFERENCE = CANDIDATES[0]
FAMILIES = {
    CANDIDATES[0]: "partial_control",
    CANDIDATES[1]: "lsgkm_shared_fit",
    CANDIDATES[2]: "lsgkm_shared_fit",
}
ROW_FIELDS = (
    "seed", "row_hash", "unit_hash", "block_hash", "stratum", "outer_fold",
    "study_id", "assay_context_id", "element_id", "source_locus_group_id",
    "long_range_block_id",
)
CONTROL_FIELDS = (
    *ROW_FIELDS, "model_id", "head_id", "prediction", "experimental_replicates",
    "biological_donors", "outcome_role",
)
RAW_FIELDS = (
    *ROW_FIELDS, "model_id", "head_id", "prediction", "prediction_unit",
    "calibration_status", "context_specific_prediction", "experimental_replicates",
    "biological_donors", "outcome_role",
)


class LSGKMEvaluationError(RuntimeError):
    """Raised when an evaluator authority, identity, or metric requirement differs."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise LSGKMEvaluationError(f"JSON object differs: {path}")
    return value


def project_file(root: Path, relative_text: str, expected: str) -> Path:
    relative = Path(relative_text)
    if relative.is_absolute() or ".." in relative.parts:
        raise LSGKMEvaluationError("unsafe project-relative path")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    if path.is_symlink() or not path.is_file() or sha256(path) != expected:
        raise LSGKMEvaluationError(f"file authority differs: {relative_text}")
    return path


def authority_member(root: Path, authority: Mapping[str, Any]) -> Path:
    artifacts = project_file(root, authority["artifacts_path"], authority["artifacts_sha256"])
    manifest = load_json(artifacts)
    roster = {
        row["path"]: row
        for row in manifest.get("artifacts", [])
        if isinstance(row, dict) and isinstance(row.get("path"), str)
    }
    member = roster.get(authority["member"])
    if member is None or member.get("sha256") != authority["member_sha256"]:
        raise LSGKMEvaluationError("artifact member record differs")
    path = (artifacts.parent / authority["member"]).resolve(strict=True)
    path.relative_to(artifacts.parent)
    if path.is_symlink() or not path.is_file() or sha256(path) != authority["member_sha256"]:
        raise LSGKMEvaluationError("artifact member bytes differ")
    return path


def read_rows(path: Path, fields: Sequence[str], compressed: bool) -> list[dict[str, str]]:
    if compressed:
        context = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        context = path.open("r", encoding="utf-8", newline="")
    with context as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise LSGKMEvaluationError(f"TSV header differs: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise LSGKMEvaluationError(f"TSV empty: {path}")
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows or path.exists():
        raise LSGKMEvaluationError(f"cannot write table: {path}")
    fields = list(rows[0])
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def finite(text: str) -> float:
    try:
        value = float(text)
    except ValueError as error:
        raise LSGKMEvaluationError("prediction is not numeric") from error
    if not math.isfinite(value):
        raise LSGKMEvaluationError("prediction is not finite")
    return value


def fisher_macro(values: Sequence[float | None]) -> float | None:
    if any(value is None or not math.isfinite(value) for value in values):
        return None
    clipped = np.clip(np.asarray(values, dtype=np.float64), -0.999999, 0.999999)
    return float(np.tanh(np.mean(np.arctanh(clipped))))


def fast_spearman(observed: np.ndarray, predicted: np.ndarray) -> float | None:
    if observed.size < 3 or np.all(observed == observed[0]) or np.all(predicted == predicted[0]):
        return None
    left = rankdata(observed)
    right = rankdata(predicted)
    left -= left.mean()
    right -= right.mean()
    denominator = float(np.sqrt(np.sum(left**2) * np.sum(right**2)))
    return None if denominator == 0 else float(np.sum(left * right) / denominator)


def bootstrap_primary(
    observed: Mapping[str, np.ndarray],
    predictions: Mapping[tuple[str, str], Mapping[str, np.ndarray]],
    block_ids: np.ndarray,
    resamples: int,
    seed: int,
) -> dict[tuple[str, str], np.ndarray]:
    unique_blocks = sorted(set(block_ids.tolist()))
    by_block = {block: np.flatnonzero(block_ids == block) for block in unique_blocks}
    rng = np.random.default_rng(seed)
    output = {candidate: np.full(resamples, np.nan) for candidate in predictions}
    for bootstrap_index in range(resamples):
        sampled = rng.integers(0, len(unique_blocks), size=len(unique_blocks))
        indices = np.concatenate([by_block[unique_blocks[index]] for index in sampled])
        for candidate, by_context in predictions.items():
            value = fisher_macro(
                [fast_spearman(observed[context][indices], by_context[context][indices]) for context in CONTEXTS]
            )
            if value is not None:
                output[candidate][bootstrap_index] = value
    return output


def display(value: float | None) -> str:
    return "not_applicable" if value is None or not math.isfinite(value) else format(value, ".17g")


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("status") != "freeze_development_evaluator_after_raw_prediction_lock"
        or [(row.get("model_id"), row.get("head_id")) for row in config.get("candidates", [])]
        != list(CANDIDATES)
        or config.get("metric_contract", {}).get("reference_candidate")
        != "available_simple_controls/allele_identity_ridge"
        or config.get("uncertainty", {}).get("resamples") != 10000
        or config.get("uncertainty", {}).get("seed") != 20260825
        or config.get("reporting_limits", {}).get("development_only") is not True
        or config.get("reporting_limits", {}).get("external_evaluation") is not False
        or config.get("reporting_limits", {}).get("champion_claim") is not False
    ):
        raise LSGKMEvaluationError("evaluation config differs")


def evaluate(
    root: Path,
    config_path: Path,
    selection_lock_artifacts: Path,
    selection_lock_sha: str,
    output: Path,
) -> dict[str, Any]:
    config = load_json(config_path)
    validate_config(config)
    if output.exists():
        raise LSGKMEvaluationError("evaluation output exists")
    output.mkdir(parents=True, mode=0o750)
    lock_relative = selection_lock_artifacts.resolve(strict=True).relative_to(root).as_posix()
    lock_artifacts = project_file(root, lock_relative, selection_lock_sha)
    lock_manifest = load_json(lock_artifacts)
    lock_records = {
        row["path"]: row
        for row in lock_manifest.get("artifacts", [])
        if isinstance(row, dict) and isinstance(row.get("path"), str)
    }
    lock_member = lock_records.get("lock/selection_lock.json")
    if lock_member is None:
        raise LSGKMEvaluationError("selection lock member absent")
    lock_path = lock_artifacts.parent / "lock/selection_lock.json"
    if sha256(lock_path) != lock_member["sha256"]:
        raise LSGKMEvaluationError("selection lock member differs")
    lock = load_json(lock_path)
    if (
        lock.get("status") != "locked_before_evaluator_control_and_outcome_access"
        or lock.get("config_sha256") != sha256(config_path)
        or lock.get("raw_prediction_member_sha256")
        != config["raw_prediction_authority"]["member_sha256"]
        or lock.get("prediction_values_opened") is not False
        or lock.get("control_prediction_values_opened") is not False
        or lock.get("outcomes_opened") is not False
    ):
        raise LSGKMEvaluationError("selection lock differs")

    row_path = authority_member(root, config["row_authority"])
    raw_path = authority_member(root, config["raw_prediction_authority"])
    control_path = authority_member(root, config["shared_control_authority"])
    outcome_path = authority_member(root, config["outcome_authority"])
    row_rows = read_rows(row_path, ROW_FIELDS, False)
    if len(row_rows) != 10330:
        raise LSGKMEvaluationError("row denominator differs")
    row_lookup = {(int(row["seed"]), row["row_hash"]): row for row in row_rows}
    elements = sorted({row["element_id"] for row in row_rows})
    element_index = {element: index for index, element in enumerate(elements)}
    metadata: dict[str, tuple[str, str]] = {}
    for row in row_rows:
        value = (row["long_range_block_id"], row["outer_fold"])
        if metadata.setdefault(row["element_id"], value) != value:
            raise LSGKMEvaluationError("element metadata differs")
    if len(elements) != 1033 or len({value[0] for value in metadata.values()}) != 239:
        raise LSGKMEvaluationError("element or block census differs")

    outcomes = load_outcomes(outcome_path, set(elements))
    observed = {
        context: np.asarray([outcomes[(element, context)]["mean"] for element in elements])
        for context in CONTEXTS
    }
    block_ids = np.asarray([metadata[element][0] for element in elements])
    values = {
        (model, head, seed, context): np.full(1033, np.nan)
        for model, head in CANDIDATES
        for seed in SEEDS
        for context in CONTEXTS
    }

    raw_rows = read_rows(raw_path, RAW_FIELDS, True)
    raw_identities = set()
    for row in raw_rows:
        candidate = (row["model_id"], row["head_id"])
        seed = int(row["seed"])
        identity = (*candidate, seed, row["row_hash"])
        authority = row_lookup.get((seed, row["row_hash"]))
        if (
            candidate not in CANDIDATES[1:]
            or seed not in SEEDS
            or identity in raw_identities
            or authority is None
            or any(row[field] != authority[field] for field in ROW_FIELDS if field != "seed")
            or row["calibration_status"] != "uncalibrated_static_rank_score"
        ):
            raise LSGKMEvaluationError("raw prediction identity differs")
        raw_identities.add(identity)
        slot = values[(*candidate, seed, row["assay_context_id"])]
        index = element_index[row["element_id"]]
        if math.isfinite(slot[index]):
            raise LSGKMEvaluationError("duplicate raw prediction")
        slot[index] = finite(row["prediction"])
    if len(raw_identities) != 20660:
        raise LSGKMEvaluationError("raw prediction coverage differs")

    selected_bytes: list[bytes] = []
    control_identities = set()
    with gzip.open(control_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != CONTROL_FIELDS:
            raise LSGKMEvaluationError("control prediction header differs")
        for row in reader:
            if (row["model_id"], row["head_id"]) != REFERENCE:
                continue
            seed = int(row["seed"])
            identity = (*REFERENCE, seed, row["row_hash"])
            authority = row_lookup.get((seed, row["row_hash"]))
            if (
                identity in control_identities
                or authority is None
                or any(row[field] != authority[field] for field in ROW_FIELDS if field != "seed")
            ):
                raise LSGKMEvaluationError("control prediction identity differs")
            control_identities.add(identity)
            selected_bytes.append(("\t".join(row[field] for field in CONTROL_FIELDS) + "\n").encode("utf-8"))
            slot = values[(*REFERENCE, seed, row["assay_context_id"])]
            index = element_index[row["element_id"]]
            if math.isfinite(slot[index]):
                raise LSGKMEvaluationError("duplicate control prediction")
            slot[index] = finite(row["prediction"])
    control_digest = hashlib.sha256(b"".join(selected_bytes)).hexdigest()
    if (
        len(control_identities) != 10330
        or control_digest
        != config["shared_control_authority"]["canonical_selected_rows_sha256"]
        or any(not np.isfinite(array).all() for array in values.values())
    ):
        raise LSGKMEvaluationError("control identity or complete prediction matrix differs")

    per_seed_rows: list[dict[str, object]] = []
    per_seed_primary: dict[tuple[str, str, int], float | None] = {}
    ensemble_predictions: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    ensemble_rows: list[dict[str, object]] = []
    point_primary: dict[tuple[str, str], float | None] = {}
    for candidate in CANDIDATES:
        model, head = candidate
        for seed in SEEDS:
            correlations = []
            for context in CONTEXTS:
                array = values[(model, head, seed, context)]
                value = None if np.all(array == array[0]) else float(spearmanr(observed[context], array).statistic)
                correlations.append(value)
                per_seed_rows.append({
                    "family": FAMILIES[candidate], "model_id": model, "head_id": head,
                    "seed": seed, "assay_context_id": context, "elements": 1033,
                    "spearman": display(value),
                })
            primary = fisher_macro(correlations)
            per_seed_primary[(*candidate, seed)] = primary
            per_seed_rows.append({
                "family": FAMILIES[candidate], "model_id": model, "head_id": head,
                "seed": seed, "assay_context_id": "macro", "elements": 2066,
                "spearman": display(primary),
            })
        by_context = {
            context: np.mean(np.vstack([values[(model, head, seed, context)] for seed in SEEDS]), axis=0)
            for context in CONTEXTS
        }
        ensemble_predictions[candidate] = by_context
        correlations = []
        for context in CONTEXTS:
            array = by_context[context]
            value = None if np.all(array == array[0]) else float(spearmanr(observed[context], array).statistic)
            correlations.append(value)
            ensemble_rows.append({
                "family": FAMILIES[candidate], "model_id": model, "head_id": head,
                "assay_context_id": context, "elements": 1033, "spearman": display(value),
            })
        point_primary[candidate] = fisher_macro(correlations)
        ensemble_rows.append({
            "family": FAMILIES[candidate], "model_id": model, "head_id": head,
            "assay_context_id": "macro", "elements": 2066,
            "spearman": display(point_primary[candidate]),
        })

    bootstrap = bootstrap_primary(
        observed,
        ensemble_predictions,
        block_ids,
        int(config["uncertainty"]["resamples"]),
        int(config["uncertainty"]["seed"]),
    )
    reference_samples = bootstrap[REFERENCE]
    summary_rows: list[dict[str, object]] = []
    result_records: list[dict[str, object]] = []
    for candidate in CANDIDATES:
        samples = bootstrap[candidate]
        valid = np.isfinite(samples)
        if int(valid.sum()) < int(config["uncertainty"]["minimum_valid_resamples"]):
            raise LSGKMEvaluationError("insufficient valid bootstrap resamples")
        low, high = (float(value) for value in np.quantile(samples[valid], (0.025, 0.975)))
        standard_error = float(np.std(samples[valid], ddof=1))
        point = point_primary[candidate]
        if point is None:
            raise LSGKMEvaluationError("primary metric unavailable")
        gain = float(point - float(point_primary[REFERENCE]))
        gain_samples = samples - reference_samples
        gain_valid = np.isfinite(gain_samples)
        gain_low, gain_high = (
            float(value) for value in np.quantile(gain_samples[gain_valid], (0.025, 0.975))
        )
        positive_seeds = sum(
            per_seed_primary[(*candidate, seed)] is not None
            and per_seed_primary[(*candidate, seed)] > per_seed_primary[(*REFERENCE, seed)]
            for seed in SEEDS
        )
        record = {
            "family": FAMILIES[candidate],
            "model_id": candidate[0],
            "head_id": candidate[1],
            "primary_metric": config["metric_contract"]["primary"],
            "primary_score": display(point),
            "primary_ci_low": display(low),
            "primary_ci_high": display(high),
            "primary_bootstrap_se": display(standard_error),
            "reference_candidate": f"{REFERENCE[0]}/{REFERENCE[1]}",
            "gain_vs_reference": display(gain),
            "gain_ci_low": display(gain_low),
            "gain_ci_high": display(gain_high),
            "positive_gain_seeds": positive_seeds,
            "bootstrap_resamples": int(valid.sum()),
            "bootstrap_unit": "long_range_block_id",
            "development_only": "true",
            "shortlist": "false",
            "external_claim": "false",
            "champion_claim": "false",
        }
        summary_rows.append(record)
        result_records.append({
            "model_id": candidate[0], "head_id": candidate[1],
            "primary_score": point, "primary_ci_low": low, "primary_ci_high": high,
            "gain_vs_reference": gain, "gain_ci_low": gain_low, "gain_ci_high": gain_high,
            "positive_gain_seeds": positive_seeds,
        })

    write_tsv(output / "per_seed_metrics.tsv", per_seed_rows)
    write_tsv(output / "ensemble_metrics.tsv", ensemble_rows)
    write_tsv(output / "paired_block_bootstrap.tsv", summary_rows)
    receipt = {
        "schema_version": "masld-bench-lsgkm-gse281364-dinucleotide-development-evaluation-v1",
        "status": "pass_locked_exposed_development_evaluation",
        "selection_lock_artifacts_sha256": selection_lock_sha,
        "raw_prediction_artifacts_sha256": config["raw_prediction_authority"]["artifacts_sha256"],
        "candidate_count": 3,
        "prediction_rows": 30990,
        "elements": 1033,
        "long_range_blocks": 239,
        "fixed_seeds": list(SEEDS),
        "contexts": list(CONTEXTS),
        "biological_donors": 0,
        "experimental_replicates_per_context": 4,
        "experimental_replicates_used_as_independent_donors": False,
        "results": result_records,
        "reference_candidate": {"model_id": REFERENCE[0], "head_id": REFERENCE[1]},
        "shared_control_canonical_selected_rows_sha256": control_digest,
        "shared_control_project_exposed": True,
        "protocol_incident_id": config["protocol_incident"]["path"].rsplit("/", 1)[-1].removesuffix(".json"),
        "bootstrap_resamples": 10000,
        "bootstrap_unit": "long_range_block_id",
        "readouts_are_independent_model_families": False,
        "mandatory_baselines_complete": False,
        "candidate_selection_or_shortlist_created": False,
        "external_evaluation": False,
        "champion_claim": False,
        "confirmatory_p_values_reported": False,
        "outcomes_read_in_isolated_evaluator": True,
        "control_predictions_read_in_isolated_evaluator": True,
        "sealed_assets_read": False,
    }
    (output / "evaluation_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--selection-lock-artifacts", type=Path, required=True)
    parser.add_argument("--selection-lock-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    evaluate(
        args.project_root.resolve(strict=True),
        args.config.resolve(strict=True),
        args.selection_lock_artifacts,
        args.selection_lock_artifacts_sha256,
        args.output,
    )


if __name__ == "__main__":
    main()

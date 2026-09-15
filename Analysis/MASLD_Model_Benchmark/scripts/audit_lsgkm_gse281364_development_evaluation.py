#!/usr/bin/env python3
"""Independently rederive the frozen exposed-development LS-GKM evaluation."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.stats import rankdata


SEEDS = (1103, 2909, 4721, 6673, 8111)
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
CANDIDATES = (
    ("available_simple_controls", "allele_identity_ridge"),
    ("lsgkm_exact_dinucleotide_null", "direct_gkmsvm"),
    ("lsgkm_exact_dinucleotide_null", "deltasvm"),
)
REFERENCE = CANDIDATES[0]
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
OUTCOME_FIELDS = (
    "element_id", "allele", "context_id", "cell_line", "condition",
    "experimental_replicate", "sample_id", "DNA", "RNA", "n_barcodes",
    "assay_state", "missing_reason", "pairing", "biological_unit", "donor_id",
)
SUMMARY_FIELDS = (
    "family", "model_id", "head_id", "primary_metric", "primary_score",
    "primary_ci_low", "primary_ci_high", "primary_bootstrap_se",
    "reference_candidate", "gain_vs_reference", "gain_ci_low", "gain_ci_high",
    "positive_gain_seeds", "bootstrap_resamples", "bootstrap_unit",
    "development_only", "shortlist", "external_claim", "champion_claim",
)
ENSEMBLE_FIELDS = (
    "family", "model_id", "head_id", "assay_context_id", "elements", "spearman",
)
PER_SEED_FIELDS = (
    "family", "model_id", "head_id", "seed", "assay_context_id", "elements",
    "spearman",
)


class EvaluationAuditError(RuntimeError):
    """Raised when independent evaluation rederivation differs."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise EvaluationAuditError(f"JSON object differs: {path}")
    return value


def project_file(root: Path, relative_text: str, expected: str) -> Path:
    relative = Path(relative_text)
    if relative.is_absolute() or ".." in relative.parts:
        raise EvaluationAuditError("unsafe project-relative path")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    if path.is_symlink() or not path.is_file() or sha256(path) != expected:
        raise EvaluationAuditError(f"file authority differs: {relative_text}")
    return path


def artifact_member(
    root: Path, authority: Mapping[str, Any], member_key: str = "member"
) -> Path:
    artifacts = project_file(root, authority["artifacts_path"], authority["artifacts_sha256"])
    manifest = load_json(artifacts)
    roster = {
        row["path"]: row
        for row in manifest.get("artifacts", [])
        if isinstance(row, dict) and isinstance(row.get("path"), str)
    }
    relative = authority[member_key]
    record = roster.get(relative)
    expected = authority[f"{member_key}_sha256"]
    if record is None or record.get("sha256") != expected:
        raise EvaluationAuditError(f"artifact member record differs: {relative}")
    path = (artifacts.parent / relative).resolve(strict=True)
    path.relative_to(artifacts.parent)
    if path.is_symlink() or not path.is_file() or sha256(path) != expected:
        raise EvaluationAuditError(f"artifact member bytes differ: {relative}")
    return path


def frozen_member(artifacts: Path, artifacts_sha: str, relative: str) -> Path:
    if sha256(artifacts) != artifacts_sha:
        raise EvaluationAuditError("frozen evaluation authority differs")
    manifest = load_json(artifacts)
    roster = {
        row["path"]: row
        for row in manifest.get("artifacts", [])
        if isinstance(row, dict) and isinstance(row.get("path"), str)
    }
    record = roster.get(relative)
    if record is None:
        raise EvaluationAuditError(f"evaluation member absent: {relative}")
    path = artifacts.parent / relative
    if path.is_symlink() or not path.is_file() or sha256(path) != record.get("sha256"):
        raise EvaluationAuditError(f"evaluation member differs: {relative}")
    return path


def read_tsv(path: Path, fields: Sequence[str], compressed: bool = False) -> list[dict[str, str]]:
    context = (
        gzip.open(path, "rt", encoding="utf-8", newline="")
        if compressed
        else path.open("r", encoding="utf-8", newline="")
    )
    with context as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise EvaluationAuditError(f"TSV header differs: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise EvaluationAuditError(f"TSV empty: {path}")
    return rows


def finite(text: str) -> float:
    try:
        value = float(text)
    except ValueError as error:
        raise EvaluationAuditError("numeric value differs") from error
    if not math.isfinite(value):
        raise EvaluationAuditError("non-finite value")
    return value


def activity_delta(ref: tuple[int, int], alt: tuple[int, int]) -> float:
    ref_dna, ref_rna = ref
    alt_dna, alt_rna = alt
    return math.log2((alt_rna + 0.5) / (alt_dna + 0.5)) - math.log2(
        (ref_rna + 0.5) / (ref_dna + 0.5)
    )


def load_outcomes(path: Path, elements: set[str]) -> dict[str, np.ndarray]:
    pairs: dict[tuple[str, str, int], dict[str, tuple[int, int]]] = defaultdict(dict)
    for row in read_tsv(path, OUTCOME_FIELDS, compressed=True):
        if row["element_id"] not in elements or row["context_id"] not in CONTEXTS:
            continue
        replicate = int(row["experimental_replicate"])
        allele = row["allele"]
        if (
            replicate not in {1, 2, 3, 4}
            or allele not in {"ref", "alt"}
            or row["cell_line"] != "HepG2"
            or row["biological_unit"] != "experimental_replicate"
            or row["donor_id"] != "not_applicable"
            or row["pairing"] != "same_sample_different_aliquot"
            or row["assay_state"] not in {"observed", "below_qc"}
        ):
            raise EvaluationAuditError("outcome topology differs")
        if row["assay_state"] == "observed":
            key = (row["element_id"], row["context_id"], replicate)
            if allele in pairs[key]:
                raise EvaluationAuditError("duplicate observed outcome allele")
            pairs[key][allele] = (int(row["DNA"]), int(row["RNA"]))
    ordered = sorted(elements)
    output: dict[str, np.ndarray] = {}
    for context in CONTEXTS:
        means = []
        for element in ordered:
            replicates = []
            for replicate in range(1, 5):
                pair = pairs.get((element, context, replicate), {})
                if set(pair) != {"ref", "alt"}:
                    raise EvaluationAuditError("four complete replicate pairs required")
                replicates.append(activity_delta(pair["ref"], pair["alt"]))
            means.append(float(np.mean(replicates)))
        output[context] = np.asarray(means, dtype=np.float64)
    return output


def spearman(left: np.ndarray, right: np.ndarray) -> float:
    if left.size < 3 or np.all(left == left[0]) or np.all(right == right[0]):
        raise EvaluationAuditError("undefined Spearman correlation")
    left_rank = rankdata(left, method="average")
    right_rank = rankdata(right, method="average")
    left_centered = left_rank - left_rank.mean()
    right_centered = right_rank - right_rank.mean()
    denominator = np.sqrt(np.sum(left_centered**2) * np.sum(right_centered**2))
    if denominator == 0:
        raise EvaluationAuditError("undefined ranked covariance")
    return float(np.sum(left_centered * right_centered) / denominator)


def fisher_macro(values: Sequence[float]) -> float:
    clipped = np.clip(np.asarray(values, dtype=np.float64), -0.999999, 0.999999)
    return float(np.tanh(np.mean(np.arctanh(clipped))))


def close(left: float, right: float, tolerance: float = 2e-14) -> bool:
    return math.isfinite(left) and math.isfinite(right) and abs(left - right) <= tolerance


def audit(
    root: Path,
    config_path: Path,
    evaluation_artifacts: Path,
    evaluation_sha: str,
    output: Path,
) -> dict[str, Any]:
    config = load_json(config_path)
    if (
        config.get("status") != "freeze_development_evaluator_after_raw_prediction_lock"
        or config.get("uncertainty", {}).get("resamples") != 10000
        or config.get("uncertainty", {}).get("seed") != 20260825
        or config.get("reporting_limits", {}).get("development_only") is not True
        or config.get("reporting_limits", {}).get("external_evaluation") is not False
        or config.get("reporting_limits", {}).get("champion_claim") is not False
    ):
        raise EvaluationAuditError("evaluation config differs")
    if output.exists():
        raise EvaluationAuditError("audit output exists")
    output.mkdir(parents=True, mode=0o750)

    evaluation_artifacts = evaluation_artifacts.resolve(strict=True)
    evaluation_artifacts.relative_to(root)
    receipt_path = frozen_member(
        evaluation_artifacts, evaluation_sha, "evaluation/evaluation_receipt.json"
    )
    summary_path = frozen_member(
        evaluation_artifacts, evaluation_sha, "evaluation/paired_block_bootstrap.tsv"
    )
    ensemble_path = frozen_member(
        evaluation_artifacts, evaluation_sha, "evaluation/ensemble_metrics.tsv"
    )
    per_seed_path = frozen_member(
        evaluation_artifacts, evaluation_sha, "evaluation/per_seed_metrics.tsv"
    )
    evaluation_receipt = load_json(receipt_path)
    reported_summary = {
        (row["model_id"], row["head_id"]): row
        for row in read_tsv(summary_path, SUMMARY_FIELDS)
    }
    reported_ensemble = {
        (row["model_id"], row["head_id"], row["assay_context_id"]): row
        for row in read_tsv(ensemble_path, ENSEMBLE_FIELDS)
    }
    reported_per_seed = {
        (row["model_id"], row["head_id"], int(row["seed"]), row["assay_context_id"]): row
        for row in read_tsv(per_seed_path, PER_SEED_FIELDS)
    }
    if (
        evaluation_receipt.get("status") != "pass_locked_exposed_development_evaluation"
        or set(reported_summary) != set(CANDIDATES)
        or len(reported_ensemble) != 9
        or len(reported_per_seed) != 45
    ):
        raise EvaluationAuditError("reported evaluation census differs")

    row_path = artifact_member(root, config["row_authority"])
    raw_path = artifact_member(root, config["raw_prediction_authority"])
    control_path = artifact_member(root, config["shared_control_authority"])
    outcome_path = artifact_member(root, config["outcome_authority"])
    row_rows = read_tsv(row_path, ROW_FIELDS)
    row_lookup = {(int(row["seed"]), row["row_hash"]): row for row in row_rows}
    elements = sorted({row["element_id"] for row in row_rows})
    element_index = {element: index for index, element in enumerate(elements)}
    metadata: dict[str, str] = {}
    for row in row_rows:
        block = row["long_range_block_id"]
        if metadata.setdefault(row["element_id"], block) != block:
            raise EvaluationAuditError("element block differs")
    if len(row_rows) != 10330 or len(elements) != 1033 or len(set(metadata.values())) != 239:
        raise EvaluationAuditError("row universe census differs")
    outcomes = load_outcomes(outcome_path, set(elements))
    block_ids = np.asarray([metadata[element] for element in elements])

    values = {
        (model, head, seed, context): np.full(1033, np.nan, dtype=np.float64)
        for model, head in CANDIDATES
        for seed in SEEDS
        for context in CONTEXTS
    }
    raw_identities: set[tuple[str, str, int, str]] = set()
    for row in read_tsv(raw_path, RAW_FIELDS, compressed=True):
        candidate = (row["model_id"], row["head_id"])
        seed = int(row["seed"])
        authority = row_lookup.get((seed, row["row_hash"]))
        identity = (*candidate, seed, row["row_hash"])
        if (
            candidate not in CANDIDATES[1:]
            or identity in raw_identities
            or authority is None
            or any(row[field] != authority[field] for field in ROW_FIELDS if field != "seed")
        ):
            raise EvaluationAuditError("raw prediction identity differs")
        raw_identities.add(identity)
        values[(*candidate, seed, row["assay_context_id"])][element_index[row["element_id"]]] = finite(
            row["prediction"]
        )
    selected_control_bytes: list[bytes] = []
    control_identities: set[tuple[str, str, int, str]] = set()
    for row in read_tsv(control_path, CONTROL_FIELDS, compressed=True):
        if (row["model_id"], row["head_id"]) != REFERENCE:
            continue
        seed = int(row["seed"])
        authority = row_lookup.get((seed, row["row_hash"]))
        identity = (*REFERENCE, seed, row["row_hash"])
        if (
            identity in control_identities
            or authority is None
            or any(row[field] != authority[field] for field in ROW_FIELDS if field != "seed")
        ):
            raise EvaluationAuditError("control prediction identity differs")
        control_identities.add(identity)
        selected_control_bytes.append(
            ("\t".join(row[field] for field in CONTROL_FIELDS) + "\n").encode("utf-8")
        )
        values[(*REFERENCE, seed, row["assay_context_id"])][element_index[row["element_id"]]] = finite(
            row["prediction"]
        )
    control_digest = hashlib.sha256(b"".join(selected_control_bytes)).hexdigest()
    if (
        len(raw_identities) != 20660
        or len(control_identities) != 10330
        or control_digest != config["shared_control_authority"]["canonical_selected_rows_sha256"]
        or any(not np.isfinite(array).all() for array in values.values())
    ):
        raise EvaluationAuditError("prediction rectangle differs")

    seed_primary: dict[tuple[str, str, int], float] = {}
    ensemble_predictions: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    point_primary: dict[tuple[str, str], float] = {}
    for candidate in CANDIDATES:
        model, head = candidate
        for seed in SEEDS:
            correlations = []
            for context in CONTEXTS:
                value = spearman(outcomes[context], values[(model, head, seed, context)])
                correlations.append(value)
                reported = finite(reported_per_seed[(*candidate, seed, context)]["spearman"])
                if not close(value, reported):
                    raise EvaluationAuditError("per-seed context metric differs")
            primary = fisher_macro(correlations)
            seed_primary[(*candidate, seed)] = primary
            reported = finite(reported_per_seed[(*candidate, seed, "macro")]["spearman"])
            if not close(primary, reported):
                raise EvaluationAuditError("per-seed primary metric differs")
        by_context = {
            context: np.mean(
                np.vstack([values[(model, head, seed, context)] for seed in SEEDS]), axis=0
            )
            for context in CONTEXTS
        }
        ensemble_predictions[candidate] = by_context
        correlations = []
        for context in CONTEXTS:
            value = spearman(outcomes[context], by_context[context])
            correlations.append(value)
            reported = finite(reported_ensemble[(*candidate, context)]["spearman"])
            if not close(value, reported):
                raise EvaluationAuditError("ensemble context metric differs")
        point_primary[candidate] = fisher_macro(correlations)
        reported = finite(reported_ensemble[(*candidate, "macro")]["spearman"])
        if not close(point_primary[candidate], reported):
            raise EvaluationAuditError("ensemble primary metric differs")

    blocks = sorted(set(block_ids.tolist()))
    by_block = {block: np.flatnonzero(block_ids == block) for block in blocks}
    rng = np.random.default_rng(20260825)
    bootstrap = {candidate: np.empty(10000, dtype=np.float64) for candidate in CANDIDATES}
    for bootstrap_index in range(10000):
        selected = rng.integers(0, len(blocks), size=len(blocks))
        indices = np.concatenate([by_block[blocks[index]] for index in selected])
        for candidate in CANDIDATES:
            bootstrap[candidate][bootstrap_index] = fisher_macro(
                [
                    spearman(outcomes[context][indices], ensemble_predictions[candidate][context][indices])
                    for context in CONTEXTS
                ]
            )

    audited_results: list[dict[str, Any]] = []
    reference_samples = bootstrap[REFERENCE]
    reference_point = point_primary[REFERENCE]
    for candidate in CANDIDATES:
        samples = bootstrap[candidate]
        low, high = (float(value) for value in np.quantile(samples, (0.025, 0.975)))
        standard_error = float(np.std(samples, ddof=1))
        gain_samples = samples - reference_samples
        gain_low, gain_high = (float(value) for value in np.quantile(gain_samples, (0.025, 0.975)))
        gain = point_primary[candidate] - reference_point
        positive = sum(seed_primary[(*candidate, seed)] > seed_primary[(*REFERENCE, seed)] for seed in SEEDS)
        reported = reported_summary[candidate]
        comparisons = {
            "primary_score": point_primary[candidate],
            "primary_ci_low": low,
            "primary_ci_high": high,
            "primary_bootstrap_se": standard_error,
            "gain_vs_reference": gain,
            "gain_ci_low": gain_low,
            "gain_ci_high": gain_high,
        }
        if any(not close(value, finite(reported[field])) for field, value in comparisons.items()):
            raise EvaluationAuditError(f"bootstrap summary differs: {candidate}")
        if int(reported["positive_gain_seeds"]) != positive:
            raise EvaluationAuditError(f"positive-gain seed count differs: {candidate}")
        audited_results.append({
            "model_id": candidate[0],
            "head_id": candidate[1],
            "primary_score": point_primary[candidate],
            "primary_ci_low": low,
            "primary_ci_high": high,
            "gain_vs_reference": gain,
            "gain_ci_low": gain_low,
            "gain_ci_high": gain_high,
            "positive_gain_seeds": positive,
        })
    reported_receipt_results = {
        (row["model_id"], row["head_id"]): row for row in evaluation_receipt.get("results", [])
    }
    for row in audited_results:
        reported = reported_receipt_results.get((row["model_id"], row["head_id"]))
        if reported is None or any(
            not close(float(value), float(reported[field]))
            for field, value in row.items()
            if field not in {"model_id", "head_id", "positive_gain_seeds"}
        ) or row["positive_gain_seeds"] != reported["positive_gain_seeds"]:
            raise EvaluationAuditError("evaluation receipt result differs")

    receipt = {
        "schema_version": "masld-bench-lsgkm-gse281364-development-evaluation-audit-v1",
        "status": "pass_independent_exposed_development_evaluation_rederivation",
        "evaluation_artifacts_sha256": evaluation_sha,
        "config_sha256": sha256(config_path),
        "candidate_count": len(CANDIDATES),
        "prediction_rows": len(raw_identities) + len(control_identities),
        "elements": len(elements),
        "long_range_blocks": len(blocks),
        "bootstrap_resamples": 10000,
        "bootstrap_seed": 20260825,
        "bootstrap_unit": "long_range_block_id",
        "results": audited_results,
        "shared_control_canonical_selected_rows_sha256": control_digest,
        "biological_donors": 0,
        "experimental_replicates_per_context": 4,
        "experimental_replicates_used_as_independent_donors": False,
        "candidate_selection_or_shortlist_created": False,
        "external_evaluation": False,
        "champion_claim": False,
        "confirmatory_p_values_reported": False,
        "outcomes_opened_for_independent_audit": True,
        "control_predictions_opened_for_independent_audit": True,
        "sealed_assets_read": False,
    }
    (output / "evaluation_audit.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--evaluation-artifacts", type=Path, required=True)
    parser.add_argument("--evaluation-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(
        args.project_root.resolve(strict=True),
        args.config.resolve(strict=True),
        args.evaluation_artifacts,
        args.evaluation_artifacts_sha256,
        args.output,
    )


if __name__ == "__main__":
    main()

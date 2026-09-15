#!/usr/bin/env python3
"""Evaluate frozen restricted-Enformer allele scores on exposed MPRA outcomes."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
from hashlib import sha256
import io
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from scipy.stats import pearsonr, spearmanr


PREDICTION_FIELDS = (
    "fixture_id",
    "element_id",
    "outer_locus_sequence_group_id",
    "outer_fold",
    "hepg2_accessibility_sad",
    "hepg2_accessibility_sar",
    "liver_accessibility_sad",
    "liver_accessibility_sar",
    "all_accessibility_mean_sad",
    "all_accessibility_mean_sar",
)
OUTCOME_FIELDS = (
    "element_id",
    "allele",
    "context_id",
    "cell_line",
    "condition",
    "experimental_replicate",
    "sample_id",
    "DNA",
    "RNA",
    "n_barcodes",
    "assay_state",
    "missing_reason",
    "pairing",
    "biological_unit",
    "donor_id",
)
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
SCORES = {
    "hepg2_accessibility_sad": {
        "transform": "signed_ALT_minus_REF_sum",
        "track_scope": "three_prespecified_ENCODE_HepG2_DNase_tracks",
        "biological_match": "cell_line_and_accessibility_assay_family_matched",
    },
    "hepg2_accessibility_sar": {
        "transform": "signed_log2_ALT_plus1_minus_log2_REF_plus1",
        "track_scope": "three_prespecified_ENCODE_HepG2_DNase_tracks",
        "biological_match": "cell_line_and_accessibility_assay_family_matched",
    },
    "liver_accessibility_sad": {
        "transform": "signed_ALT_minus_REF_sum",
        "track_scope": "prespecified_hepatocyte_and_adult_liver_DNase_tracks",
        "biological_match": "hepatic_lineage_transfer_not_cell_line_matched",
    },
    "liver_accessibility_sar": {
        "transform": "signed_log2_ALT_plus1_minus_log2_REF_plus1",
        "track_scope": "prespecified_hepatocyte_and_adult_liver_DNase_tracks",
        "biological_match": "hepatic_lineage_transfer_not_cell_line_matched",
    },
}
EXCLUDED_SCORES = (
    "all_accessibility_mean_sad",
    "all_accessibility_mean_sar",
)
PSEUDOCOUNT = 0.5
BOOTSTRAP_REPLICATES = 1000
SEED = 20260824


class EnformerEvaluationError(RuntimeError):
    """Raised when an Enformer development-evaluation requirement is not met."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise EnformerEvaluationError(f"JSON object required: {path}")
    return value


def verify_frozen_tree(root: Path, expected_sha256: str) -> dict[str, Any]:
    if not root.is_dir() or root.is_symlink():
        raise EnformerEvaluationError(f"invalid frozen root: {root}")
    manifest_path, complete_path = root / "ARTIFACTS.json", root / "COMPLETE"
    if (
        not manifest_path.is_file()
        or not complete_path.is_file()
        or manifest_path.is_symlink()
        or complete_path.is_symlink()
        or sha256_file(manifest_path) != expected_sha256
    ):
        raise EnformerEvaluationError(f"frozen control files differ: {root}")
    manifest = _load_json(manifest_path)
    complete = _load_json(complete_path)
    artifacts = manifest.get("artifacts")
    if (
        set(manifest) != {"schema_version", "metadata", "artifacts"}
        or manifest.get("schema_version") != "masld-bench-artifacts-v1"
        or not isinstance(manifest.get("metadata"), dict)
        or not isinstance(artifacts, list)
        or complete
        != {
            "artifact_count": len(artifacts),
            "manifest_sha256": expected_sha256,
            "schema_version": "masld-bench-complete-v1",
        }
    ):
        raise EnformerEvaluationError(f"frozen control contract differs: {root}")
    expected: dict[str, tuple[str, int]] = {}
    for item in artifacts:
        if not isinstance(item, dict) or set(item) != {"path", "sha256", "size_bytes"}:
            raise EnformerEvaluationError("artifact member differs")
        relative = Path(str(item["path"]))
        key = relative.as_posix()
        if relative.is_absolute() or ".." in relative.parts or key in expected:
            raise EnformerEvaluationError("unsafe or duplicate artifact path")
        expected[key] = (str(item["sha256"]), int(item["size_bytes"]))
    observed: set[str] = set()
    for path in root.rglob("*"):
        if path.is_symlink():
            raise EnformerEvaluationError(f"symlink in frozen tree: {path}")
        relative = path.relative_to(root).as_posix()
        if path.is_file() and relative not in {"ARTIFACTS.json", "COMPLETE"}:
            observed.add(relative)
    if observed != set(expected):
        raise EnformerEvaluationError(f"frozen inventory differs: {root}")
    for relative, (digest, size) in expected.items():
        path = root / relative
        if path.stat().st_size != size or sha256_file(path) != digest:
            raise EnformerEvaluationError(f"frozen member differs: {relative}")
    return manifest


def activity_delta(
    ref_dna: int,
    ref_rna: int,
    alt_dna: int,
    alt_rna: int,
    pseudocount: float = PSEUDOCOUNT,
) -> float:
    if min(ref_dna, ref_rna, alt_dna, alt_rna) < 0 or pseudocount <= 0:
        raise EnformerEvaluationError("invalid reporter count")
    return math.log2((alt_rna + pseudocount) / (alt_dna + pseudocount)) - math.log2(
        (ref_rna + pseudocount) / (ref_dna + pseudocount)
    )


def read_predictions(root: Path) -> list[dict[str, Any]]:
    receipt = _load_json(root / "predictions/receipt.json")
    if (
        receipt.get("status") != "pass_outcome_blind_restricted_prediction"
        or receipt.get("model_id") != "enformer_crested_restricted_port"
        or receipt.get("registered_scientific_identity")
        != "restricted_conversion_not_native_sonnet"
        or receipt.get("elements") != 1033
        or receipt.get("outer_locus_sequence_groups") != 1033
        or receipt.get("outer_folds") != 5
        or receipt.get("tracks") != 5313
        or receipt.get("accessibility_tracks") != 684
        or receipt.get("allele_delta_sign") != "ALT_minus_REF"
        or receipt.get("deterministic_repeat_max_abs") != 0.0
        or receipt.get("outcomes_read")
        or receipt.get("reporter_counts_read")
        or receipt.get("sealed_outcomes_read")
        or receipt.get("model_fitted_or_adapted")
        or receipt.get("native_sonnet_parity_established")
        or receipt.get("champion_eligible")
        or receipt.get("terms")
        != "internal_academic_noncommercial_nontransferable_restricted_comparator"
    ):
        raise EnformerEvaluationError("prediction receipt differs")
    path = root / "predictions/predictions.tsv"
    if sha256_file(path) != receipt.get("predictions_sha256"):
        raise EnformerEvaluationError("prediction table checksum differs")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != PREDICTION_FIELDS:
            raise EnformerEvaluationError("prediction fields differ")
        raw = [dict(row) for row in reader]
    rows: list[dict[str, Any]] = []
    for row in raw:
        fold = int(row["outer_fold"])
        scores = {field: float(row[field]) for field in SCORES}
        excluded = {field: float(row[field]) for field in EXCLUDED_SCORES}
        if fold not in range(5) or not all(math.isfinite(value) for value in (*scores.values(), *excluded.values())):
            raise EnformerEvaluationError("prediction fold or scalar differs")
        rows.append({**row, "outer_fold": fold, **scores, **excluded})
    if (
        len(rows) != 1033
        or len({row["fixture_id"] for row in rows}) != 1033
        or len({row["element_id"] for row in rows}) != 1033
        or len({row["outer_locus_sequence_group_id"] for row in rows}) != 1033
        or {row["outer_fold"] for row in rows} != set(range(5))
    ):
        raise EnformerEvaluationError("prediction census differs")
    return rows


def verify_fixture_alignment(
    predictions: Sequence[Mapping[str, Any]], fixture_root: Path
) -> None:
    with (fixture_root / "fixture/manifest.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        fixture = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    observed = [
        (
            row["fixture_id"],
            row["element_id"],
            row["outer_locus_sequence_group_id"],
            int(row["outer_fold"]),
        )
        for row in predictions
    ]
    expected = [
        (
            row["fixture_id"],
            row["element_id"],
            row["outer_locus_sequence_group_id"],
            int(row["outer_fold"]),
        )
        for row in fixture
    ]
    if observed != expected:
        raise EnformerEvaluationError("prediction-to-fixture alignment differs")


def load_outcomes(
    path: Path, eligible_ids: set[str]
) -> dict[tuple[str, str], dict[str, Any]]:
    pairs: dict[tuple[str, str, int], dict[str, tuple[int, int]]] = defaultdict(dict)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != OUTCOME_FIELDS:
            raise EnformerEvaluationError("outcome fields differ")
        for row in reader:
            if row["element_id"] not in eligible_ids or row["context_id"] not in CONTEXTS:
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
                raise EnformerEvaluationError("outcome topology differs")
            if row["assay_state"] != "observed":
                continue
            key = (row["element_id"], row["context_id"], replicate)
            if allele in pairs[key]:
                raise EnformerEvaluationError("duplicate observed outcome allele")
            pairs[key][allele] = (int(row["DNA"]), int(row["RNA"]))
    output: dict[tuple[str, str], dict[str, Any]] = {}
    for element in sorted(eligible_ids):
        for context in CONTEXTS:
            values: list[float] = []
            for replicate in range(1, 5):
                pair = pairs.get((element, context, replicate), {})
                if set(pair) != {"ref", "alt"}:
                    raise EnformerEvaluationError("four complete experimental replicate pairs required")
                values.append(activity_delta(*pair["ref"], *pair["alt"]))
            output[(element, context)] = {
                "mean": float(np.mean(values)),
                "sd": float(np.std(values, ddof=1)),
                "replicates": tuple(values),
            }
    return output


def correlations(observed: Sequence[float], predicted: Sequence[float]) -> tuple[float, float]:
    truth = np.asarray(observed, dtype=np.float64)
    score = np.asarray(predicted, dtype=np.float64)
    if (
        truth.shape != score.shape
        or truth.ndim != 1
        or truth.size < 3
        or not np.isfinite(truth).all()
        or not np.isfinite(score).all()
        or np.all(truth == truth[0])
        or np.all(score == score[0])
    ):
        raise EnformerEvaluationError("correlation input differs or is constant")
    return (
        float(pearsonr(truth, score).statistic),
        float(spearmanr(truth, score).statistic),
    )


def group_bootstrap(
    observed: Sequence[float],
    predicted: Sequence[float],
    groups: Sequence[str],
    *,
    seed: int = SEED,
    replicates: int = BOOTSTRAP_REPLICATES,
) -> dict[str, float | int]:
    truth = np.asarray(observed, dtype=np.float64)
    score = np.asarray(predicted, dtype=np.float64)
    group = np.asarray(groups, dtype=str)
    unique = np.unique(group)
    if truth.shape != score.shape or truth.shape != group.shape or len(unique) < 3:
        raise EnformerEvaluationError("bootstrap input differs")
    locations = {value: np.flatnonzero(group == value) for value in unique}
    rng = np.random.default_rng(seed)
    pearsons: list[float] = []
    spearmans: list[float] = []
    for _ in range(replicates):
        sampled = rng.choice(unique, size=len(unique), replace=True)
        indices = np.concatenate([locations[value] for value in sampled])
        try:
            pearson, spearman = correlations(truth[indices], score[indices])
        except EnformerEvaluationError:
            continue
        pearsons.append(pearson)
        spearmans.append(spearman)
    if len(pearsons) < int(0.95 * replicates):
        raise EnformerEvaluationError("too many invalid bootstrap replicates")
    return {
        "pearson_ci_low": float(np.quantile(pearsons, 0.025)),
        "pearson_ci_high": float(np.quantile(pearsons, 0.975)),
        "spearman_ci_low": float(np.quantile(spearmans, 0.025)),
        "spearman_ci_high": float(np.quantile(spearmans, 0.975)),
        "valid_bootstrap_replicates": len(pearsons),
    }


def matched_rows(
    predictions: Sequence[Mapping[str, Any]],
    outcomes: Mapping[tuple[str, str], Mapping[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for prediction in predictions:
        for score_id, contract in SCORES.items():
            for context in CONTEXTS:
                outcome = outcomes[(str(prediction["element_id"]), context)]
                rows.append(
                    {
                        "model_id": "enformer_crested_restricted_port",
                        "score_id": score_id,
                        "score_transform": contract["transform"],
                        "track_scope": contract["track_scope"],
                        "biological_match": contract["biological_match"],
                        "context_id": context,
                        "context_position": (
                            "direct_control_context"
                            if context == "HepG2_control"
                            else "PAOA_condition_transport_not_condition_aware_prediction"
                        ),
                        "element_id": prediction["element_id"],
                        "outer_locus_sequence_group_id": prediction[
                            "outer_locus_sequence_group_id"
                        ],
                        "outer_fold": prediction["outer_fold"],
                        "prediction_score": format(float(prediction[score_id]), ".17g"),
                        "observed_mean_signed_log2_activity_delta": format(
                            float(outcome["mean"]), ".17g"
                        ),
                        "observed_replicate_sd": format(float(outcome["sd"]), ".17g"),
                        "experimental_replicates": 4,
                        "biological_donors": 0,
                        "assay_transfer": "accessibility_prediction_to_MPRA_reporter_activity",
                        "champion_eligible": "false",
                    }
                )
    return rows


def metric_rows(rows: Sequence[Mapping[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    aggregate: list[dict[str, Any]] = []
    per_fold: list[dict[str, Any]] = []
    for score_id in SCORES:
        for context in CONTEXTS:
            selected = [
                row for row in rows if row["score_id"] == score_id and row["context_id"] == context
            ]
            observed = [float(row["observed_mean_signed_log2_activity_delta"]) for row in selected]
            predicted = [float(row["prediction_score"]) for row in selected]
            groups = [str(row["outer_locus_sequence_group_id"]) for row in selected]
            pearson, spearman = correlations(observed, predicted)
            interval = group_bootstrap(
                observed,
                predicted,
                groups,
                seed=SEED + int(sha256(f"{score_id}|{context}".encode()).hexdigest()[:8], 16),
            )
            contract = SCORES[score_id]
            aggregate.append(
                {
                    "model_id": "enformer_crested_restricted_port",
                    "score_id": score_id,
                    "context_id": context,
                    "elements": len(selected),
                    "locus_sequence_groups": len(set(groups)),
                    "pearson": format(pearson, ".17g"),
                    "pearson_ci_low": format(float(interval["pearson_ci_low"]), ".17g"),
                    "pearson_ci_high": format(float(interval["pearson_ci_high"]), ".17g"),
                    "spearman": format(spearman, ".17g"),
                    "spearman_ci_low": format(float(interval["spearman_ci_low"]), ".17g"),
                    "spearman_ci_high": format(float(interval["spearman_ci_high"]), ".17g"),
                    "valid_group_bootstrap_replicates": interval[
                        "valid_bootstrap_replicates"
                    ],
                    "track_scope": contract["track_scope"],
                    "biological_match": contract["biological_match"],
                    "models_ranked": "false",
                    "multiple_testing": "none_descriptive_development_metrics",
                    "champion_eligible": "false",
                }
            )
            for fold in range(5):
                scoped = [row for row in selected if int(row["outer_fold"]) == fold]
                fold_pearson, fold_spearman = correlations(
                    [float(row["observed_mean_signed_log2_activity_delta"]) for row in scoped],
                    [float(row["prediction_score"]) for row in scoped],
                )
                per_fold.append(
                    {
                        "model_id": "enformer_crested_restricted_port",
                        "score_id": score_id,
                        "context_id": context,
                        "outer_fold": fold,
                        "elements": len(scoped),
                        "locus_sequence_groups": len(
                            {row["outer_locus_sequence_group_id"] for row in scoped}
                        ),
                        "pearson": format(fold_pearson, ".17g"),
                        "spearman": format(fold_spearman, ".17g"),
                        "models_ranked": "false",
                    }
                )
    return aggregate, per_fold


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise EnformerEvaluationError(f"cannot write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_gzip_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise EnformerEvaluationError(f"cannot write empty table: {path}")
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)


def evaluate(arguments: argparse.Namespace) -> dict[str, Any]:
    if arguments.output.exists():
        raise EnformerEvaluationError("evaluation output exists")
    prediction_manifest = verify_frozen_tree(
        arguments.predictions, arguments.predictions_sha256
    )
    outcome_manifest = verify_frozen_tree(arguments.outcomes, arguments.outcomes_sha256)
    fixture_manifest = verify_frozen_tree(arguments.fixture, arguments.fixture_sha256)
    if (
        prediction_manifest["metadata"].get("artifact_class")
        != "enformer_gse281364_outcome_blind_prediction"
        or prediction_manifest["metadata"].get("outcomes_read")
        or prediction_manifest["metadata"].get("model_fitted_or_adapted")
        or prediction_manifest["metadata"].get("native_sonnet_parity_established")
        or prediction_manifest["metadata"].get("champion_eligible")
        or outcome_manifest["metadata"].get("artifact_class")
        != "gse281364_replicate_safe_outcomes"
        or outcome_manifest["metadata"].get("donor_count") != 0
        or outcome_manifest["metadata"].get("outcome_role")
        != "exposed_development_MPRA_only"
        or outcome_manifest["metadata"].get("sealed_outcomes_loaded")
        or fixture_manifest["metadata"].get("artifact_class")
        != "gse281364_outcome_blind_enformer_fixture"
        or fixture_manifest["metadata"].get("outcomes_read")
    ):
        raise EnformerEvaluationError("input role or firewall differs")
    predictions = read_predictions(arguments.predictions)
    verify_fixture_alignment(predictions, arguments.fixture)
    outcomes = load_outcomes(
        arguments.outcomes / "outcomes/replicate_outcomes.tsv.gz",
        {str(row["element_id"]) for row in predictions},
    )
    matched = matched_rows(predictions, outcomes)
    aggregate, per_fold = metric_rows(matched)
    if len(matched) != 8264 or len(aggregate) != 8 or len(per_fold) != 40:
        raise EnformerEvaluationError("evaluation census differs")

    arguments.output.mkdir(parents=True, mode=0o750)
    write_gzip_tsv(arguments.output / "matched_endpoints.tsv.gz", matched)
    write_tsv(arguments.output / "aggregate_metrics.tsv", aggregate)
    write_tsv(arguments.output / "outer_fold_metrics.tsv", per_fold)
    caveats = [
        {
            "item": "model_identity",
            "disposition": "restricted_comparator_only",
            "reason": "CREsted Keras conversion; native Sonnet numerical parity not established",
        },
        {
            "item": "assay_transfer",
            "disposition": "descriptive_correlation_only",
            "reason": "Enformer accessibility predictions are not calibrated MPRA reporter activities",
        },
        {
            "item": "HepG2_PAOA",
            "disposition": "condition_transport_only",
            "reason": "prediction is sequence-only and has no PAOA treatment input",
        },
        {
            "item": "LX2_control;LX2_TGFb",
            "disposition": "excluded_biological_context_mismatch",
            "reason": "no prespecified LX2 or stellate accessibility track scalarization",
        },
        {
            "item": ";".join(EXCLUDED_SCORES),
            "disposition": "excluded_biological_specificity_mismatch",
            "reason": "mean across 684 accessibility tracks is not a HepG2 or liver-specific score",
        },
        {
            "item": "experimental_replicates",
            "disposition": "averaged_within_element_context",
            "reason": "four repeats are not independent biological donors; donor n=0",
        },
    ]
    write_tsv(arguments.output / "caveats.tsv", caveats)
    receipt = {
        "schema_version": "masld-bench-gse281364-enformer-development-evaluation-v1",
        "status": "pass_descriptive_development_assay_transfer",
        "dataset_id": "gse281364",
        "model_id": "enformer_crested_restricted_port",
        "registered_scientific_identity": "restricted_conversion_not_native_sonnet",
        "prediction_artifacts_sha256": arguments.predictions_sha256,
        "outcome_artifacts_sha256": arguments.outcomes_sha256,
        "fixture_artifacts_sha256": arguments.fixture_sha256,
        "elements": 1033,
        "outer_locus_sequence_groups": 1033,
        "outer_folds": 5,
        "evaluated_contexts": list(CONTEXTS),
        "excluded_contexts": ["LX2_control", "LX2_TGFb"],
        "evaluated_scores": list(SCORES),
        "excluded_scores": list(EXCLUDED_SCORES),
        "endpoint": "mean_signed_log2_ALT_minus_REF_RNA_over_DNA_activity",
        "experimental_replicates_per_context": 4,
        "experimental_replicates_used_as_independent_donors": False,
        "donor_count": 0,
        "bootstrap_unit": "outer_locus_sequence_group",
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "outcome_role": "exposed_development_MPRA_only",
        "assay_transfer": "accessibility_prediction_to_MPRA_reporter_activity",
        "models_ranked": False,
        "confirmatory_inference": False,
        "multiple_testing": "none_descriptive_development_metrics",
        "native_sonnet_parity_established": False,
        "restricted_comparator": True,
        "external_evaluation": False,
        "sealed_outcomes_read": False,
        "champion_eligible": False,
        "clinical_claim_supported": False,
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--predictions-sha256", required=True)
    parser.add_argument("--outcomes", type=Path, required=True)
    parser.add_argument("--outcomes-sha256", required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--fixture-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    evaluate(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

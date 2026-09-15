#!/usr/bin/env python3
"""Evaluate frozen Sei and MPRALegNet predictions on exposed MPRA outcomes."""

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

PSEUDOCOUNT = 0.5
BOOTSTRAP_REPLICATES = 1000
SEED = 20260824
CONTEXTS = ("HepG2_control", "HepG2_PAOA", "LX2_control", "LX2_TGFb")
SEI_FIELDS = (
    "fixture_id",
    "element_id",
    "outer_locus_sequence_group_id",
    "outer_fold",
    "max_abs_sequence_class_index0",
    "signed_max_abs_score",
    "max_abs_score",
    "mean_signed_score",
)
MPRA_FIELDS = (
    "element_id",
    "outer_locus_sequence_group_id",
    "outer_fold",
    "contig",
    "variant_pos0",
    "ref",
    "alt",
    "cell_context",
    "input_lane",
    "reference_forward_score",
    "alternative_forward_score",
    "reference_reverse_complement_score",
    "alternative_reverse_complement_score",
    "reference_orientation_mean_score",
    "alternative_orientation_mean_score",
    "alt_minus_ref_score",
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


class GSE281364EvaluationError(RuntimeError):
    """Raised when evaluation would violate a frozen data or endpoint requirement."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def verify_frozen_tree(root: Path) -> dict[str, Any]:
    if not root.is_dir() or root.is_symlink():
        raise GSE281364EvaluationError(f"invalid frozen-tree root: {root}")
    manifest_path = root / "ARTIFACTS.json"
    complete_path = root / "COMPLETE"
    if (
        not manifest_path.is_file()
        or manifest_path.is_symlink()
        or not complete_path.is_file()
        or complete_path.is_symlink()
    ):
        raise GSE281364EvaluationError(f"incomplete frozen tree: {root}")
    manifest = _load_json(manifest_path)
    complete = _load_json(complete_path)
    raw_artifacts = manifest.get("artifacts")
    if (
        set(manifest) != {"schema_version", "metadata", "artifacts"}
        or manifest.get("schema_version") != "masld-bench-artifacts-v1"
        or not isinstance(manifest.get("metadata"), dict)
        or not isinstance(raw_artifacts, list)
        or complete
        != {
            "artifact_count": len(raw_artifacts),
            "manifest_sha256": sha256_file(manifest_path),
            "schema_version": "masld-bench-complete-v1",
        }
    ):
        raise GSE281364EvaluationError(f"frozen-tree control file differs: {root}")
    expected: dict[str, tuple[str, int]] = {}
    for item in raw_artifacts:
        if not isinstance(item, dict) or set(item) != {"path", "sha256", "size_bytes"}:
            raise GSE281364EvaluationError("artifact manifest member differs")
        relative = Path(str(item["path"]))
        if relative.is_absolute() or ".." in relative.parts or relative.as_posix() in expected:
            raise GSE281364EvaluationError("unsafe or duplicate artifact path")
        expected[relative.as_posix()] = (str(item["sha256"]), int(item["size_bytes"]))
    observed: set[str] = set()
    for path in root.rglob("*"):
        if path.is_symlink():
            raise GSE281364EvaluationError(f"symlink in frozen tree: {path}")
        relative = path.relative_to(root).as_posix()
        if path.is_file() and relative not in {"ARTIFACTS.json", "COMPLETE"}:
            observed.add(relative)
    if observed != set(expected):
        raise GSE281364EvaluationError("frozen-tree file inventory differs")
    for relative, (expected_hash, expected_size) in expected.items():
        path = root / relative
        if path.stat().st_size != expected_size or sha256_file(path) != expected_hash:
            raise GSE281364EvaluationError(f"frozen artifact differs: {relative}")
    return manifest


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE281364EvaluationError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def _load_manifest(root: Path, expected_sha256: str) -> dict[str, Any]:
    manifest = verify_frozen_tree(root)
    if sha256_file(root / "ARTIFACTS.json") != expected_sha256:
        raise GSE281364EvaluationError(f"artifact manifest differs: {root}")
    return manifest


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise GSE281364EvaluationError(f"JSON object required: {path}")
    return value


def activity_delta(
    reference_dna: int,
    reference_rna: int,
    alternative_dna: int,
    alternative_rna: int,
    pseudocount: float = PSEUDOCOUNT,
) -> float:
    values = (reference_dna, reference_rna, alternative_dna, alternative_rna)
    if min(values) < 0 or pseudocount <= 0:
        raise GSE281364EvaluationError("invalid reporter count or pseudocount")
    reference = math.log2(
        (reference_rna + pseudocount) / (reference_dna + pseudocount)
    )
    alternative = math.log2(
        (alternative_rna + pseudocount) / (alternative_dna + pseudocount)
    )
    return alternative - reference


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
        raise GSE281364EvaluationError("correlation inputs differ or are constant")
    pearson = float(pearsonr(truth, score).statistic)
    spearman = float(spearmanr(truth, score).statistic)
    if not math.isfinite(pearson) or not math.isfinite(spearman):
        raise GSE281364EvaluationError("correlation is not finite")
    return pearson, spearman


def group_bootstrap_interval(
    observed: Sequence[float],
    predicted: Sequence[float],
    groups: Sequence[str],
    *,
    replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = SEED,
) -> dict[str, float]:
    truth = np.asarray(observed, dtype=np.float64)
    score = np.asarray(predicted, dtype=np.float64)
    group = np.asarray(groups, dtype=str)
    if truth.shape != score.shape or truth.shape != group.shape or replicates < 10:
        raise GSE281364EvaluationError("bootstrap inputs differ")
    unique = np.unique(group)
    if unique.size < 3:
        raise GSE281364EvaluationError("bootstrap requires at least three locus groups")
    indices = {value: np.flatnonzero(group == value) for value in unique}
    generator = np.random.default_rng(seed)
    pearsons: list[float] = []
    spearmans: list[float] = []
    for _ in range(replicates):
        sampled = generator.choice(unique, size=unique.size, replace=True)
        selected = np.concatenate([indices[value] for value in sampled])
        try:
            pearson, spearman = correlations(truth[selected], score[selected])
        except GSE281364EvaluationError:
            continue
        pearsons.append(pearson)
        spearmans.append(spearman)
    if len(pearsons) < int(0.95 * replicates):
        raise GSE281364EvaluationError("too many invalid group-bootstrap replicates")
    return {
        "pearson_ci_low": float(np.quantile(pearsons, 0.025)),
        "pearson_ci_high": float(np.quantile(pearsons, 0.975)),
        "spearman_ci_low": float(np.quantile(spearmans, 0.025)),
        "spearman_ci_high": float(np.quantile(spearmans, 0.975)),
        "valid_bootstrap_replicates": len(pearsons),
    }


def _prediction_rows(
    sei_root: Path, mpralegnet_root: Path
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    sei_receipt = _load_json(sei_root / "predictions/receipt.json")
    mpralegnet_receipt = _load_json(mpralegnet_root / "predictions/receipt.json")
    if (
        sei_receipt.get("status") != "pass_outcome_blind_native_prediction"
        or sei_receipt.get("elements") != 1033
        or sei_receipt.get("outcomes_read")
        or sei_receipt.get("sealed_outcomes_read")
        or sei_receipt.get("champion_eligible")
        or mpralegnet_receipt.get("status") != "pass_outcome_blind_prediction"
        or mpralegnet_receipt.get("elements") != 4359
        or mpralegnet_receipt.get("outcomes_read")
        or mpralegnet_receipt.get("sealed_outcomes_read")
        or mpralegnet_receipt.get("standalone_champion_eligible")
        or mpralegnet_receipt.get("input_length_bp") != 230
        or mpralegnet_receipt.get("checkpoint_provenance")
        != "third_party_Hugging_Face_repackaging_author_equivalence_unresolved"
    ):
        raise GSE281364EvaluationError("prediction receipt differs")
    sei_fields, sei_raw = read_tsv(sei_root / "predictions/predictions.tsv")
    mpralegnet_fields, mpralegnet_raw = read_tsv(
        mpralegnet_root / "predictions/allele_scores.tsv"
    )
    if sei_fields != SEI_FIELDS or mpralegnet_fields != MPRA_FIELDS:
        raise GSE281364EvaluationError("prediction table schema differs")
    if len(sei_raw) != 1033 or len(mpralegnet_raw) != 4359:
        raise GSE281364EvaluationError("prediction row census differs")

    def normalize(rows: list[dict[str, str]], score_field: str) -> list[dict[str, Any]]:
        output: list[dict[str, Any]] = []
        for row in rows:
            fold = int(row["outer_fold"])
            score = float(row[score_field])
            if fold not in range(5) or not math.isfinite(score):
                raise GSE281364EvaluationError("prediction fold or score differs")
            output.append({**row, "outer_fold": fold, "prediction_score": score})
        if len({row["element_id"] for row in output}) != len(output):
            raise GSE281364EvaluationError("prediction element is duplicated")
        return output

    sei = normalize(sei_raw, "max_abs_score")
    mpralegnet = normalize(mpralegnet_raw, "alt_minus_ref_score")
    if any(
        row["cell_context"] != "HepG2"
        or row["input_lane"] != "GRCh38p14_230bp_genomic_window_transfer"
        for row in mpralegnet
    ):
        raise GSE281364EvaluationError("MPRALegNet context or transfer lane differs")
    sei_keys = {
        row["element_id"]: (
            row["outer_locus_sequence_group_id"],
            row["outer_fold"],
        )
        for row in sei
    }
    mpralegnet_keys = {
        row["element_id"]: (
            row["outer_locus_sequence_group_id"],
            row["outer_fold"],
        )
        for row in mpralegnet
    }
    common = set(sei_keys).intersection(mpralegnet_keys)
    if len(common) != 1033 or any(sei_keys[key] != mpralegnet_keys[key] for key in common):
        raise GSE281364EvaluationError("shared prediction grouping differs")
    return sei, mpralegnet


def _outcome_deltas(
    path: Path, eligible_ids: set[str]
) -> tuple[dict[tuple[str, str], dict[str, Any]], dict[str, int]]:
    values: dict[tuple[str, str, int], dict[str, tuple[int, int]]] = defaultdict(dict)
    states: dict[tuple[str, str], dict[str, int]] = defaultdict(
        lambda: {"observed": 0, "below_qc": 0}
    )
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != OUTCOME_FIELDS:
            raise GSE281364EvaluationError("outcome table schema differs")
        for row in reader:
            element = row["element_id"]
            if element not in eligible_ids:
                continue
            context = row["context_id"]
            allele = row["allele"]
            replicate = int(row["experimental_replicate"])
            if (
                context not in CONTEXTS
                or allele not in {"ref", "alt"}
                or replicate not in {1, 2, 3, 4}
                or row["biological_unit"] != "experimental_replicate"
                or row["donor_id"] != "not_applicable"
            ):
                raise GSE281364EvaluationError("outcome topology differs")
            state = row["assay_state"]
            if state not in {"observed", "below_qc"}:
                raise GSE281364EvaluationError("outcome assay state differs")
            states[(element, context)][state] += 1
            if state != "observed":
                continue
            key = (element, context, replicate)
            if allele in values[key]:
                raise GSE281364EvaluationError("duplicate outcome allele row")
            values[key][allele] = (int(row["DNA"]), int(row["RNA"]))
    output: dict[tuple[str, str], dict[str, Any]] = {}
    coverage = {context: 0 for context in CONTEXTS}
    for element in sorted(eligible_ids):
        for context in CONTEXTS:
            deltas: list[float] = []
            for replicate in range(1, 5):
                pair = values.get((element, context, replicate), {})
                if set(pair) != {"ref", "alt"}:
                    break
                deltas.append(activity_delta(*pair["ref"], *pair["alt"]))
            if len(deltas) != 4:
                continue
            coverage[context] += 1
            output[(element, context)] = {
                "replicate_deltas": deltas,
                "mean_signed_delta": float(np.mean(deltas)),
                "sd_signed_delta": float(np.std(deltas, ddof=1)),
                "replicates": 4,
                "observed_rows": states[(element, context)]["observed"],
                "below_qc_rows": states[(element, context)]["below_qc"],
            }
    return output, coverage


def _matched_endpoint_rows(
    predictions: list[dict[str, Any]],
    outcomes: Mapping[tuple[str, str], Mapping[str, Any]],
    *,
    model_id: str,
    endpoint_id: str,
    contexts: Sequence[str],
    magnitude: bool,
    comparability: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for prediction in predictions:
        for context in contexts:
            outcome = outcomes.get((prediction["element_id"], context))
            if outcome is None:
                continue
            signed = float(outcome["mean_signed_delta"])
            rows.append(
                {
                    "model_id": model_id,
                    "endpoint_id": endpoint_id,
                    "context_id": context,
                    "element_id": prediction["element_id"],
                    "outer_locus_sequence_group_id": prediction[
                        "outer_locus_sequence_group_id"
                    ],
                    "outer_fold": prediction["outer_fold"],
                    "prediction_score": prediction["prediction_score"],
                    "observed_score": abs(signed) if magnitude else signed,
                    "observed_mean_signed_log2_activity_delta": signed,
                    "observed_replicate_sd": outcome["sd_signed_delta"],
                    "experimental_replicates": outcome["replicates"],
                    "outcome_unit": "element_mean_of_four_experimental_replicate_allele_deltas",
                    "donor_inference": "not_applicable",
                    "comparability": comparability,
                    "champion_eligible": "false",
                }
            )
    return rows


def _metric_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    endpoint_contexts = sorted(
        {(row["model_id"], row["endpoint_id"], row["context_id"]) for row in rows}
    )
    for model_id, endpoint_id, context in endpoint_contexts:
        selected = [
            row
            for row in rows
            if (row["model_id"], row["endpoint_id"], row["context_id"])
            == (model_id, endpoint_id, context)
        ]
        for fold in (None, 0, 1, 2, 3, 4):
            scoped = selected if fold is None else [row for row in selected if row["outer_fold"] == fold]
            pearson, spearman = correlations(
                [row["observed_score"] for row in scoped],
                [row["prediction_score"] for row in scoped],
            )
            interval: dict[str, float | int] = {
                "pearson_ci_low": math.nan,
                "pearson_ci_high": math.nan,
                "spearman_ci_low": math.nan,
                "spearman_ci_high": math.nan,
                "valid_bootstrap_replicates": 0,
            }
            if fold is None:
                interval = group_bootstrap_interval(
                    [row["observed_score"] for row in scoped],
                    [row["prediction_score"] for row in scoped],
                    [row["outer_locus_sequence_group_id"] for row in scoped],
                )
            output.append(
                {
                    "model_id": model_id,
                    "endpoint_id": endpoint_id,
                    "context_id": context,
                    "scope": "all_outer_folds" if fold is None else "outer_fold",
                    "outer_fold": "all" if fold is None else fold,
                    "elements": len(scoped),
                    "locus_sequence_groups": len(
                        {row["outer_locus_sequence_group_id"] for row in scoped}
                    ),
                    "pearson": format(pearson, ".17g"),
                    "pearson_ci_low": (
                        "not_applicable"
                        if fold is not None
                        else format(float(interval["pearson_ci_low"]), ".17g")
                    ),
                    "pearson_ci_high": (
                        "not_applicable"
                        if fold is not None
                        else format(float(interval["pearson_ci_high"]), ".17g")
                    ),
                    "spearman": format(spearman, ".17g"),
                    "spearman_ci_low": (
                        "not_applicable"
                        if fold is not None
                        else format(float(interval["spearman_ci_low"]), ".17g")
                    ),
                    "spearman_ci_high": (
                        "not_applicable"
                        if fold is not None
                        else format(float(interval["spearman_ci_high"]), ".17g")
                    ),
                    "valid_group_bootstrap_replicates": interval[
                        "valid_bootstrap_replicates"
                    ],
                    "multiple_testing": "none_descriptive_development_metrics",
                    "models_ranked": "false",
                }
            )
    return output


def evaluate(arguments: argparse.Namespace) -> dict[str, Any]:
    if arguments.output.exists():
        raise GSE281364EvaluationError("evaluation output exists")
    sei_manifest = _load_manifest(arguments.sei, arguments.sei_artifacts_sha256)
    mpralegnet_manifest = _load_manifest(
        arguments.mpralegnet, arguments.mpralegnet_artifacts_sha256
    )
    outcome_manifest = _load_manifest(
        arguments.outcomes, arguments.outcomes_artifacts_sha256
    )
    if (
        sei_manifest["metadata"].get("artifact_class")
        != "sei_gse281364_outcome_blind_prediction"
        or sei_manifest["metadata"].get("outcomes_read")
        or mpralegnet_manifest["metadata"].get("artifact_class")
        != "outcome_blind_MPRALegNet_HepG2_reporter_score_comparator"
        or mpralegnet_manifest["metadata"].get("outcomes_read")
        or mpralegnet_manifest["metadata"].get("standalone_champion_eligible")
        or outcome_manifest["metadata"].get("artifact_class")
        != "gse281364_replicate_safe_outcomes"
        or outcome_manifest["metadata"].get("donor_count") != 0
        or outcome_manifest["metadata"].get("outcome_role")
        != "exposed_development_MPRA_only"
        or outcome_manifest["metadata"].get("sealed_outcomes_loaded")
    ):
        raise GSE281364EvaluationError("artifact metadata or outcome firewall differs")
    sei, mpralegnet = _prediction_rows(arguments.sei, arguments.mpralegnet)
    eligible_ids = {row["element_id"] for row in sei + mpralegnet}
    outcomes, context_coverage = _outcome_deltas(
        arguments.outcomes / "outcomes/replicate_outcomes.tsv.gz", eligible_ids
    )
    if any(context_coverage[context] != len(eligible_ids) for context in CONTEXTS):
        raise GSE281364EvaluationError("matched prediction outcomes are incomplete")
    matched: list[dict[str, Any]] = []
    matched.extend(
        _matched_endpoint_rows(
            sei,
            outcomes,
            model_id="sei",
            endpoint_id="unsigned_regulatory_effect_magnitude",
            contexts=CONTEXTS,
            magnitude=True,
            comparability=(
                "sequence-class magnitude versus absolute reporter allele activity; "
                "not reporter calibration or signed cell-context prediction"
            ),
        )
    )
    matched.extend(
        _matched_endpoint_rows(
            mpralegnet,
            outcomes,
            model_id="mpralegnet_hepg2_test1_val2",
            endpoint_id="signed_hepg2_reporter_activity_delta",
            contexts=("HepG2_control", "HepG2_PAOA"),
            magnitude=False,
            comparability=(
                "signed HepG2 reporter transfer; PAOA is condition transport; 230bp genomic "
                "window differs from native 107bp assay insert and checkpoint equivalence is unresolved"
            ),
        )
    )
    metrics = _metric_rows(matched)
    coverage = [
        {
            "model_id": "sei",
            "prediction_rows": len(sei),
            "matched_elements_per_context": len(sei),
            "locus_sequence_groups": len(
                {row["outer_locus_sequence_group_id"] for row in sei}
            ),
            "eligible_contexts": 4,
            "excluded_contexts_for_native_mismatch": 0,
            "unmatched_prediction_rows": 0,
            "incomplete_replicate_pairs": 0,
        },
        {
            "model_id": "mpralegnet_hepg2_test1_val2",
            "prediction_rows": len(mpralegnet),
            "matched_elements_per_context": len(mpralegnet),
            "locus_sequence_groups": len(
                {row["outer_locus_sequence_group_id"] for row in mpralegnet}
            ),
            "eligible_contexts": 2,
            "excluded_contexts_for_native_mismatch": 2,
            "unmatched_prediction_rows": 0,
            "incomplete_replicate_pairs": 0,
        },
    ]
    comparability = [
        {
            "model_id": "sei",
            "native_prediction": "maximum_absolute_native_sequence_class_score",
            "outcome_endpoint": "absolute_mean_signed_log2_MPRA_allele_activity_delta",
            "contexts": ";".join(CONTEXTS),
            "position": "unsigned_regulatory_magnitude_transport",
            "cross_model_rankable": "false",
            "champion_eligible": "false",
            "limitation": "not an assay-calibrated reporter score and not cell-context signed",
        },
        {
            "model_id": "mpralegnet_hepg2_test1_val2",
            "native_prediction": "ALT_minus_REF_uncalibrated_HepG2_reporter_score",
            "outcome_endpoint": "mean_signed_log2_MPRA_allele_activity_delta",
            "contexts": "HepG2_control;HepG2_PAOA",
            "position": "secondary_nonchampion_reporter_transfer",
            "cross_model_rankable": "false",
            "champion_eligible": "false",
            "limitation": (
                "230bp genomic transfer differs from native 107bp assay insert; "
                "checkpoint author equivalence unresolved"
            ),
        },
    ]
    arguments.output.mkdir(parents=True, mode=0o750)
    write_tsv(arguments.output / "coverage.tsv", tuple(coverage[0]), coverage)
    write_tsv(
        arguments.output / "comparability.tsv", tuple(comparability[0]), comparability
    )
    metric_fields = tuple(metrics[0])
    write_tsv(arguments.output / "metrics.tsv", metric_fields, metrics)
    matched_fields = tuple(matched[0])
    matched_path = arguments.output / "matched_endpoints.tsv.gz"
    with matched_path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(
                    text,
                    fieldnames=list(matched_fields),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(matched)
    summary = {
        f"{row['model_id']}::{row['context_id']}": {
            "pearson": float(row["pearson"]),
            "spearman": float(row["spearman"]),
            "elements": row["elements"],
            "locus_sequence_groups": row["locus_sequence_groups"],
        }
        for row in metrics
        if row["scope"] == "all_outer_folds"
    }
    receipt: dict[str, Any] = {
        "schema_version": "masld-bench-gse281364-sei-mpralegnet-development-evaluation-v1",
        "status": "pass_descriptive_development_evaluation",
        "dataset_id": "gse281364",
        "prediction_artifacts_sha256": {
            "sei": arguments.sei_artifacts_sha256,
            "mpralegnet": arguments.mpralegnet_artifacts_sha256,
        },
        "outcome_artifacts_sha256": arguments.outcomes_artifacts_sha256,
        "outcome_role": "exposed_development_MPRA_only",
        "sealed_outcomes_read": False,
        "experimental_replicates_per_context": 4,
        "replicates_used_as_independent_donors": False,
        "donor_count": 0,
        "activity_transform": "ALT_minus_REF_log2((RNA+0.5)/(DNA+0.5))_per_replicate",
        "replicate_aggregation": "element_context_mean_of_four_experimental_replicate_deltas",
        "locus_sequence_grouping": "outer_locus_sequence_group_cluster_bootstrap",
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "seed": SEED,
        "sei_elements": len(sei),
        "mpralegnet_elements": len(mpralegnet),
        "locus_sequence_groups": 1033,
        "models_forced_to_common_endpoint": False,
        "models_ranked": False,
        "mpralegnet_secondary_nonchampion": True,
        "champion_eligible": False,
        "multiple_testing": "none_descriptive_development_metrics",
        "metric_summary": summary,
        "coverage_sha256": sha256_file(arguments.output / "coverage.tsv"),
        "comparability_sha256": sha256_file(arguments.output / "comparability.tsv"),
        "metrics_sha256": sha256_file(arguments.output / "metrics.tsv"),
        "matched_endpoints_sha256": sha256_file(matched_path),
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sei", type=Path, required=True)
    parser.add_argument("--sei-artifacts-sha256", required=True)
    parser.add_argument("--mpralegnet", type=Path, required=True)
    parser.add_argument("--mpralegnet-artifacts-sha256", required=True)
    parser.add_argument("--outcomes", type=Path, required=True)
    parser.add_argument("--outcomes-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    evaluate(parser.parse_args())


if __name__ == "__main__":
    main()

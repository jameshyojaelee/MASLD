#!/usr/bin/env python3
"""Fail-closed GSE274114 within-platform baseline trainer/evaluator separation."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
import warnings

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.contracts import ArtifactRef, PredictionBundle
from masld_bench.hashing import canonical_sha256, sha256_file


FIXTURE_ARTIFACTS_SHA256 = (
    "b334ff4fe77903be82bb228fc479a28d92d74f4cec0b06ea98c3d6a6100c9300"
)
DATASET_ID = "gse274114_mash_hbv"
FIXED_STOCHASTIC_SEEDS = (1103, 2207, 3301, 4409, 5501)
DETERMINISTIC_SEED = 0
MODEL_KINDS = (
    "training_class_prior",
    "gene_rank_nearest_centroid",
    "hvg_pca_elastic_net",
    "hvg_pca_linear_svm",
)
STOCHASTIC_MODELS = frozenset(
    {"hvg_pca_elastic_net", "hvg_pca_linear_svm"}
)
TASKS: dict[str, dict[str, Any]] = {
    "gse274114_hiseq_healthy_vs_hbv": {
        "groups": {"CTRL": 0, "ENEG": 1},
        "participants": 20,
        "group_counts": {"CTRL": 9, "ENEG": 11},
        "split_id": "gse274114_hiseq_ctrl_eneg_stratified_participant_outer_v1",
        "model_prefix": "hiseq_",
    },
    "gse274114_novaseq_mash_vs_mash_hbv": {
        "groups": {"NASH": 0, "ENEG_NASH": 1},
        "participants": 19,
        "group_counts": {"NASH": 10, "ENEG_NASH": 9},
        "split_id": "gse274114_novaseq_nash_eneg_nash_stratified_participant_outer_v1",
        "model_prefix": "novaseq_",
    },
}


class GSE274114BaselineFirewallError(RuntimeError):
    """Raised before a trainer/evaluator boundary can be crossed unsafely."""


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE274114BaselineFirewallError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(row))


def save_array(path: Path, value: Any) -> None:
    import numpy as np

    with path.open("xb") as handle:
        np.save(handle, np.asarray(value), allow_pickle=False)


def artifact_ref(path: Path, *, root: Path, role: str, media_type: str) -> ArtifactRef:
    return ArtifactRef.from_path(
        path, relative_to=root, role=role, media_type=media_type
    )


def _fixture_manifest_without_member_reads(fixture_root: Path) -> dict[str, Any]:
    manifest_path = fixture_root / "ARTIFACTS.json"
    if sha256_file(manifest_path) != FIXTURE_ARTIFACTS_SHA256:
        raise GSE274114BaselineFirewallError("fixture ARTIFACTS identity differs")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    completion = json.loads((fixture_root / "COMPLETE").read_text(encoding="utf-8"))
    if (
        manifest.get("schema_version") != "masld-bench-artifacts-v1"
        or completion.get("manifest_sha256") != FIXTURE_ARTIFACTS_SHA256
        or completion.get("artifact_count") != len(manifest.get("artifacts", ()))
    ):
        raise GSE274114BaselineFirewallError("fixture freeze controls differ")
    return manifest


def _quant_relative_path(value: str) -> Path:
    relative = Path(value)
    if (
        relative.is_absolute()
        or ".." in relative.parts
        or relative.parts[:2] != ("model_inputs", "quant_v49")
        or len(relative.parts) != 3
        or not relative.name.endswith(".quant.genes.sf.gz")
    ):
        raise GSE274114BaselineFirewallError(
            "quant path is outside the molecular-only allowlist"
        )
    return relative


def stage_unlabeled_features(*, fixture_root: Path, output: Path) -> dict[str, Any]:
    """Materialize molecular values without opening labels, folds, or metadata."""

    import numpy as np

    if output.exists():
        raise GSE274114BaselineFirewallError("refusing to overwrite feature stage")
    manifest = _fixture_manifest_without_member_reads(fixture_root)
    artifacts = {row["path"]: row for row in manifest["artifacts"]}
    quant_manifest_path = fixture_root / "model_inputs/quant_manifest.tsv"
    quant_record = artifacts.get("model_inputs/quant_manifest.tsv")
    if (
        not isinstance(quant_record, Mapping)
        or sha256_file(quant_manifest_path) != quant_record.get("sha256")
    ):
        raise GSE274114BaselineFirewallError("quant manifest identity differs")
    fields, rows = read_tsv(quant_manifest_path)
    expected_fields = (
        "row_id",
        "quant_path",
        "quant_sha256",
        "fields",
        "gene_axis",
        "source_annotation",
        "target_annotation",
        "par_y_x_aggregation_applied",
        "unmapped_source_rows_masked",
        "normalization_or_fit_applied",
        "labels_present",
    )
    if fields != expected_fields or len(rows) != 39:
        raise GSE274114BaselineFirewallError("unlabeled quant manifest differs")
    row_ids: list[str] = []
    gene_ids: list[str] | None = None
    counts: list[list[float]] = []
    effective_lengths: list[list[float]] = []
    for row in rows:
        if (
            row["labels_present"] != "False"
            or row["normalization_or_fit_applied"] != "False"
            or row["gene_axis"] != "GENCODE_v49_60324_unique_targets"
        ):
            raise GSE274114BaselineFirewallError("quant input semantics differ")
        relative = _quant_relative_path(row["quant_path"])
        path = fixture_root / relative
        record = artifacts.get(relative.as_posix())
        if (
            not isinstance(record, Mapping)
            or row["quant_sha256"] != record.get("sha256")
            or path.stat().st_size != record.get("size_bytes")
            or sha256_file(path) != record.get("sha256")
        ):
            raise GSE274114BaselineFirewallError("quant artifact identity differs")
        current_genes: list[str] = []
        current_counts: list[float] = []
        current_lengths: list[float] = []
        with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != (
                "Name",
                "Length",
                "EffectiveLength",
                "TPM",
                "NumReads",
            ):
                raise GSE274114BaselineFirewallError("Salmon header differs")
            for value in reader:
                length = float(value["EffectiveLength"])
                count = float(value["NumReads"])
                if not math.isfinite(length) or length <= 0 or not math.isfinite(count) or count < 0:
                    raise GSE274114BaselineFirewallError("Salmon value is invalid")
                current_genes.append(value["Name"])
                current_lengths.append(length)
                current_counts.append(count)
        if len(current_genes) != 60_324 or len(set(current_genes)) != 60_324:
            raise GSE274114BaselineFirewallError("GENCODE v49 gene axis differs")
        if gene_ids is None:
            gene_ids = current_genes
        elif current_genes != gene_ids:
            raise GSE274114BaselineFirewallError("participant gene order differs")
        row_ids.append(row["row_id"])
        counts.append(current_counts)
        effective_lengths.append(current_lengths)
    assert gene_ids is not None
    count_array = np.asarray(counts, dtype=np.float64)
    length_array = np.asarray(effective_lengths, dtype=np.float64)
    if count_array.shape != (39, 60_324) or length_array.shape != count_array.shape:
        raise GSE274114BaselineFirewallError("staged feature matrix shape differs")
    output.mkdir(parents=True)
    write_tsv(output / "participants.tsv", ("row_id",), ({"row_id": row} for row in row_ids))
    write_tsv(output / "genes.tsv", ("gene_id",), ({"gene_id": gene} for gene in gene_ids))
    save_array(output / "numreads.npy", count_array)
    save_array(output / "effective_length.npy", length_array)
    receipt = {
        "schema_version": "masld-bench-gse274114-unlabeled-features-v1",
        "participants": 39,
        "genes": 60_324,
        "labels_opened": False,
        "folds_opened": False,
        "participant_metadata_opened": False,
        "molecular_values_opened": True,
        "normalization_run": False,
        "model_fit_or_scoring_run": False,
        "effective_length_retained": True,
    }
    write_json_exclusive(output / "receipt.json", receipt)
    freeze_tree(output, {"artifact_class": "gse274114_unlabeled_feature_matrix", **receipt})
    return receipt


def prepare_fold_view(
    *,
    feature_root: Path,
    labels_path: Path,
    folds_path: Path,
    task_id: str,
    outer_fold: int,
    output: Path,
) -> dict[str, Any]:
    """Evaluator-side view builder; held labels are never written to the view."""

    import numpy as np

    if output.exists() or task_id not in TASKS or outer_fold not in range(5):
        raise GSE274114BaselineFirewallError("invalid or existing fold-view target")
    manifest = verify_frozen_tree(feature_root)
    if (
        manifest["metadata"].get("artifact_class")
        != "gse274114_unlabeled_feature_matrix"
        or manifest["metadata"].get("labels_opened") is not False
        or manifest["metadata"].get("normalization_run") is not False
    ):
        raise GSE274114BaselineFirewallError("feature-stage firewall differs")
    participant_fields, participants = read_tsv(feature_root / "participants.tsv")
    gene_fields, genes = read_tsv(feature_root / "genes.tsv")
    label_fields, labels = read_tsv(labels_path)
    fold_fields, folds = read_tsv(folds_path)
    if participant_fields != ("row_id",) or gene_fields != ("gene_id",):
        raise GSE274114BaselineFirewallError("feature-stage tables differ")
    if label_fields != (
        "row_id",
        "source_group",
        "source_group_semantics",
        "mash_state",
        "hbv_state",
        "within_instrument_contrast",
        "expected_instrument",
    ) or fold_fields != (
        "row_id",
        "outer_fold",
        "source_group",
        "fold_assignment_used_labels",
        "available_to_model_input",
    ):
        raise GSE274114BaselineFirewallError("evaluator table schema differs")
    row_ids = [row["row_id"] for row in participants]
    if len(row_ids) != len(set(row_ids)) or set(row_ids) != {row["row_id"] for row in labels}:
        raise GSE274114BaselineFirewallError("participant/label join differs")
    fold_by_row = {row["row_id"]: row for row in folds}
    label_by_row = {row["row_id"]: row for row in labels}
    if set(fold_by_row) != set(row_ids) or any(
        row["available_to_model_input"] != "False"
        or row["fold_assignment_used_labels"] != "True"
        or row["source_group"] != label_by_row[row_id]["source_group"]
        for row_id, row in fold_by_row.items()
    ):
        raise GSE274114BaselineFirewallError("evaluator fold contract differs")
    task = TASKS[task_id]
    eligible = [row for row in row_ids if label_by_row[row]["source_group"] in task["groups"]]
    observed_counts = Counter(label_by_row[row]["source_group"] for row in eligible)
    if len(eligible) != task["participants"] or observed_counts != Counter(task["group_counts"]):
        raise GSE274114BaselineFirewallError("task participant census differs")
    train_rows = [row for row in eligible if int(fold_by_row[row]["outer_fold"]) != outer_fold]
    query_rows = [row for row in eligible if int(fold_by_row[row]["outer_fold"]) == outer_fold]
    train_targets = [task["groups"][label_by_row[row]["source_group"]] for row in train_rows]
    query_targets = [task["groups"][label_by_row[row]["source_group"]] for row in query_rows]
    if not query_rows or set(train_targets) != {0, 1} or set(query_targets) != {0, 1}:
        raise GSE274114BaselineFirewallError("outer fold lacks a binary class")
    index = {row: position for position, row in enumerate(row_ids)}
    counts = np.load(feature_root / "numreads.npy", allow_pickle=False)
    lengths = np.load(feature_root / "effective_length.npy", allow_pickle=False)
    expected_shape = (len(row_ids), len(genes))
    if counts.shape != expected_shape or lengths.shape != expected_shape:
        raise GSE274114BaselineFirewallError("feature arrays differ")
    output.mkdir(parents=True)
    write_tsv(
        output / "training_participants.tsv",
        ("row_id", "target"),
        ({"row_id": row, "target": target} for row, target in zip(train_rows, train_targets)),
    )
    write_tsv(
        output / "query_participants.tsv",
        ("row_id",),
        ({"row_id": row} for row in query_rows),
    )
    write_tsv(output / "genes.tsv", ("gene_id",), genes)
    save_array(output / "training_numreads.npy", counts[[index[row] for row in train_rows]])
    save_array(output / "query_numreads.npy", counts[[index[row] for row in query_rows]])
    save_array(output / "training_effective_length.npy", lengths[[index[row] for row in train_rows]])
    save_array(output / "query_effective_length.npy", lengths[[index[row] for row in query_rows]])
    receipt = {
        "schema_version": "masld-bench-gse274114-trainer-view-v1",
        "task_id": task_id,
        "split_id": task["split_id"],
        "outer_fold": outer_fold,
        "training_participants": len(train_rows),
        "query_participants": len(query_rows),
        "training_labels_included": True,
        "query_labels_included": False,
        "source_group_strings_included": False,
        "participant_covariates_included": False,
        "age_included": False,
        "sex_included": False,
        "group_aggregate_metadata_included": False,
        "normalization_or_feature_selection_run": False,
        "model_fit_or_scoring_run": False,
        "biological_unit": "participant",
        "technical_runs_are_replicates": False,
        "evaluator_labels_sha256": sha256_file(labels_path),
        "evaluator_folds_sha256": sha256_file(folds_path),
    }
    write_json_exclusive(output / "receipt.json", receipt)
    freeze_tree(output, {"artifact_class": "gse274114_trainer_fold_view", **receipt})
    return receipt


def _log_cpm(counts: Any) -> Any:
    import numpy as np

    values = np.asarray(counts, dtype=np.float64)
    totals = values.sum(axis=1)
    if values.ndim != 2 or not np.all(np.isfinite(values)) or np.any(values < 0) or np.any(totals <= 0):
        raise GSE274114BaselineFirewallError("NumReads matrix is invalid")
    return np.log1p(values / totals[:, None] * 1_000_000.0)


def _fit_transform(train: Any, query: Any, *, maximum_features: int = 2_000) -> tuple[Any, Any, dict[str, Any]]:
    import numpy as np

    train_log = _log_cpm(train)
    query_log = _log_cpm(query)
    cpm = np.expm1(train_log)
    minimum_samples = max(2, math.ceil(len(train_log) * 0.20))
    eligible = np.flatnonzero(np.sum(cpm >= 1.0, axis=0) >= minimum_samples)
    if len(eligible) < 2:
        raise GSE274114BaselineFirewallError("too few training-filtered genes")
    variances = np.var(train_log[:, eligible], axis=0, ddof=1)
    order = np.lexsort((eligible, -variances))[: min(maximum_features, len(eligible))]
    selected = eligible[order]
    mean = np.mean(train_log[:, selected], axis=0)
    scale = np.std(train_log[:, selected], axis=0, ddof=1)
    scale = np.where(scale > 0, scale, 1.0)
    train_scaled = (train_log[:, selected] - mean) / scale
    query_scaled = (query_log[:, selected] - mean) / scale
    components = min(10, len(train_scaled) - 2, train_scaled.shape[1])
    if components < 2:
        raise GSE274114BaselineFirewallError("training fold cannot support PCA")
    _, _, right = np.linalg.svd(train_scaled, full_matrices=False)
    loadings = right[:components]
    signature = canonical_sha256(
        {
            "selected_indices": selected.tolist(),
            "training_mean": mean.tolist(),
            "training_scale": scale.tolist(),
            "pca_loadings": loadings.tolist(),
        }
    )
    return train_scaled @ loadings.T, query_scaled @ loadings.T, {
        "selected_features": int(len(selected)),
        "pca_components": int(components),
        "training_preprocessing_sha256": signature,
    }


def _fit_gene_rank_transform(
    train: Any, query: Any, *, maximum_features: int = 2_000
) -> tuple[Any, Any, dict[str, Any]]:
    """Select genes on training participants, then rank within each participant."""

    import numpy as np

    train_log = _log_cpm(train)
    query_log = _log_cpm(query)
    cpm = np.expm1(train_log)
    minimum_samples = max(2, math.ceil(len(train_log) * 0.20))
    eligible = np.flatnonzero(np.sum(cpm >= 1.0, axis=0) >= minimum_samples)
    if len(eligible) < 2:
        raise GSE274114BaselineFirewallError("too few training-filtered genes")
    variances = np.var(train_log[:, eligible], axis=0, ddof=1)
    order = np.lexsort((eligible, -variances))[: min(maximum_features, len(eligible))]
    selected = eligible[order]

    def average_ranks(values: Any) -> Any:
        matrix = np.asarray(values, dtype=np.float64)
        result = np.empty_like(matrix)
        for row_index, row in enumerate(matrix):
            ranked_order = np.argsort(row, kind="mergesort")
            sorted_values = row[ranked_order]
            start = 0
            while start < len(row):
                end = start + 1
                while end < len(row) and sorted_values[end] == sorted_values[start]:
                    end += 1
                result[row_index, ranked_order[start:end]] = (start + end - 1) / 2
                start = end
        return result / max(1, matrix.shape[1] - 1)

    train_rank = average_ranks(train_log[:, selected])
    query_rank = average_ranks(query_log[:, selected])
    return train_rank, query_rank, {
        "selected_features": int(len(selected)),
        "pca_components": 0,
        "training_preprocessing_sha256": canonical_sha256(
            {"selected_indices": selected.tolist(), "transform": "within_participant_average_gene_rank"}
        ),
    }


def _seed_for(model_kind: str, seed: int) -> None:
    expected = FIXED_STOCHASTIC_SEEDS if model_kind in STOCHASTIC_MODELS else (0,)
    if seed not in expected:
        raise GSE274114BaselineFirewallError(
            f"{model_kind} seed {seed} is outside the frozen seed policy"
        )


def _fit_model_with_convergence_audit(
    model: Any, values: Any, labels: Any
) -> bool:
    from sklearn.exceptions import ConvergenceWarning

    with warnings.catch_warnings(record=True) as caught:
        warnings.filterwarnings("ignore", category=FutureWarning)
        warnings.simplefilter("always", ConvergenceWarning)
        model.fit(values, labels)
    return not any(
        issubclass(record.category, ConvergenceWarning) for record in caught
    )


def _predict_baseline(
    train_counts: Any,
    query_counts: Any,
    targets: Any,
    *,
    model_kind: str,
    seed: int,
) -> tuple[Any, dict[str, Any]]:
    import numpy as np

    _seed_for(model_kind, seed)
    y = np.asarray(targets, dtype=np.int64)
    if set(y.tolist()) != {0, 1}:
        raise GSE274114BaselineFirewallError("training labels are not binary")
    if model_kind == "training_class_prior":
        probability = np.repeat(float(np.mean(y)), len(query_counts))
        return probability, {"training_preprocessing_sha256": canonical_sha256({"targets": y.tolist()}), "selected_features": 0, "pca_components": 0}

    if model_kind == "gene_rank_nearest_centroid":
        train_rank, query_rank, preprocessing = _fit_gene_rank_transform(
            train_counts, query_counts
        )
        centroids = np.vstack(
            [np.mean(train_rank[y == value], axis=0) for value in (0, 1)]
        )
        training_distances = np.stack(
            [np.sum((train_rank - centroid) ** 2, axis=1) for centroid in centroids],
            axis=1,
        )
        distances = np.stack(
            [np.sum((query_rank - centroid) ** 2, axis=1) for centroid in centroids],
            axis=1,
        )
        scale = max(
            float(
                np.median(
                    np.abs(training_distances[:, 1] - training_distances[:, 0])
                )
            ),
            1e-8,
        )
        logits = np.clip((distances[:, 0] - distances[:, 1]) / scale, -30.0, 30.0)
        return 1.0 / (1.0 + np.exp(-logits)), preprocessing

    train_pca, query_pca, preprocessing = _fit_transform(train_counts, query_counts)

    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.svm import SVC

    candidates: list[tuple[float, float | None]]
    if model_kind == "hvg_pca_elastic_net":
        candidates = [(c, ratio) for c in (0.1, 1.0, 10.0) for ratio in (0.0, 0.5, 1.0)]
    elif model_kind == "hvg_pca_linear_svm":
        candidates = [(c, None) for c in (0.1, 1.0, 10.0)]
    else:
        raise GSE274114BaselineFirewallError(f"unsupported model kind: {model_kind}")
    splitter = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    scored: list[tuple[float, float, float | None]] = []
    candidate_audit: list[dict[str, Any]] = []
    for c, ratio in candidates:
        losses: list[float] = []
        invalid_reason: str | None = None
        completed_inner_folds = 0
        for inner_train, inner_valid in splitter.split(train_counts, y):
            inner_x, valid_x, _ = _fit_transform(
                np.asarray(train_counts)[inner_train], np.asarray(train_counts)[inner_valid]
            )
            if model_kind == "hvg_pca_elastic_net":
                fit = LogisticRegression(
                    penalty="elasticnet",
                    solver="saga",
                    C=c,
                    l1_ratio=ratio,
                    class_weight="balanced",
                    max_iter=20_000,
                    random_state=seed,
                    tol=1e-6,
                )
            else:
                fit = SVC(
                    C=c,
                    kernel="linear",
                    probability=True,
                    class_weight="balanced",
                    random_state=seed,
                )
            try:
                converged = _fit_model_with_convergence_audit(
                    fit, inner_x, y[inner_train]
                )
            except (ArithmeticError, FloatingPointError, ValueError) as error:
                invalid_reason = f"fit_error:{type(error).__name__}"
                losses = []
                break
            if not converged:
                invalid_reason = "convergence_warning"
                losses = []
                break
            p = np.clip(fit.predict_proba(valid_x)[:, 1], 1e-8, 1 - 1e-8)
            truth = y[inner_valid]
            losses.append(float(np.mean(-(truth * np.log(p) + (1 - truth) * np.log(1 - p)))))
            completed_inner_folds += 1
        if losses:
            scored.append((float(np.mean(losses)), c, ratio))
            candidate_audit.append(
                {
                    "C": c,
                    "l1_ratio": ratio,
                    "status": "valid",
                    "reason": "all_inner_folds_converged",
                    "completed_inner_folds": len(losses),
                    "mean_inner_log_loss": float(np.mean(losses)),
                }
            )
        else:
            candidate_audit.append(
                {
                    "C": c,
                    "l1_ratio": ratio,
                    "status": "invalid",
                    "reason": invalid_reason or "fit_failed_before_metric",
                    "completed_inner_folds": completed_inner_folds,
                    "mean_inner_log_loss": None,
                }
            )
    if not scored:
        raise GSE274114BaselineFirewallError(
            "no nested hyperparameter candidate converged"
        )
    ranked = sorted(
        scored,
        key=lambda item: (
            item[0], item[1], -1 if item[2] is None else item[2]
        ),
    )
    final: Any | None = None
    best_c: float | None = None
    best_ratio: float | None = None
    selected_inner_rank: int | None = None
    outer_candidate_audit: list[dict[str, Any]] = []
    for rank, (mean_loss, candidate_c, candidate_ratio) in enumerate(
        ranked, start=1
    ):
        if model_kind == "hvg_pca_elastic_net":
            candidate_fit = LogisticRegression(
                penalty="elasticnet",
                solver="saga",
                C=candidate_c,
                l1_ratio=candidate_ratio,
                class_weight="balanced",
                max_iter=20_000,
                random_state=seed,
                tol=1e-6,
            )
        else:
            candidate_fit = SVC(
                C=candidate_c,
                kernel="linear",
                probability=True,
                class_weight="balanced",
                random_state=seed,
            )
        outer_status = "converged_without_warning"
        try:
            converged = _fit_model_with_convergence_audit(
                candidate_fit, train_pca, y
            )
        except (ArithmeticError, FloatingPointError, ValueError) as error:
            converged = False
            outer_status = f"fit_error:{type(error).__name__}"
        if not converged and outer_status == "converged_without_warning":
            outer_status = "convergence_warning"
        outer_candidate_audit.append(
            {
                "inner_rank": rank,
                "C": candidate_c,
                "l1_ratio": candidate_ratio,
                "mean_inner_log_loss": mean_loss,
                "outer_training_status": outer_status,
            }
        )
        if converged:
            final = candidate_fit
            best_c = candidate_c
            best_ratio = candidate_ratio
            selected_inner_rank = rank
            break
    if final is None or best_c is None or selected_inner_rank is None:
        raise GSE274114BaselineFirewallError(
            "no inner-valid hyperparameter candidate converged on outer training partition"
        )
    probability = final.predict_proba(query_pca)[:, 1]
    preprocessing.update(
        {
            "selected_c": best_c,
            "selected_l1_ratio": best_ratio,
            "hyperparameter_candidates_total": len(candidate_audit),
            "hyperparameter_candidates_valid": sum(
                item["status"] == "valid" for item in candidate_audit
            ),
            "hyperparameter_candidates_invalid": sum(
                item["status"] == "invalid" for item in candidate_audit
            ),
            "hyperparameter_candidate_audit": candidate_audit,
            "no_valid_candidate_policy": (
                "rank_inner_valid_candidates_then_fail_outer_fold_if_none_converge"
            ),
            "outer_training_candidate_audit": outer_candidate_audit,
            "outer_training_candidates_attempted": len(outer_candidate_audit),
            "selected_inner_candidate_rank": selected_inner_rank,
            "training_only_convergence_fallback_used": selected_inner_rank > 1,
            "selected_outer_fit_converged_without_warning": True,
        }
    )
    return probability, preprocessing


def fit_predict(
    *, view_root: Path, expected_view_sha256: str, model_kind: str, seed: int, output: Path
) -> dict[str, Any]:
    """Trainer entry point. It cannot accept an evaluator label path."""

    import numpy as np

    if output.exists() or model_kind not in MODEL_KINDS:
        raise GSE274114BaselineFirewallError("invalid or existing prediction target")
    if sha256_file(view_root / "ARTIFACTS.json") != expected_view_sha256:
        raise GSE274114BaselineFirewallError("trainer view identity differs")
    manifest = verify_frozen_tree(view_root)
    metadata = manifest["metadata"]
    if (
        metadata.get("artifact_class") != "gse274114_trainer_fold_view"
        or metadata.get("query_labels_included") is not False
        or metadata.get("source_group_strings_included") is not False
        or metadata.get("participant_covariates_included") is not False
        or metadata.get("normalization_or_feature_selection_run") is not False
        or metadata.get("model_fit_or_scoring_run") is not False
    ):
        raise GSE274114BaselineFirewallError("trainer view leaks protected inputs")
    train_fields, train_rows = read_tsv(view_root / "training_participants.tsv")
    query_fields, query_rows = read_tsv(view_root / "query_participants.tsv")
    if train_fields != ("row_id", "target") or query_fields != ("row_id",):
        raise GSE274114BaselineFirewallError("trainer participant schema differs")
    if {row["row_id"] for row in train_rows} & {row["row_id"] for row in query_rows}:
        raise GSE274114BaselineFirewallError("participant appears in train and query")
    train_counts = np.load(view_root / "training_numreads.npy", allow_pickle=False)
    query_counts = np.load(view_root / "query_numreads.npy", allow_pickle=False)
    targets = np.asarray([int(row["target"]) for row in train_rows], dtype=np.int64)
    if train_counts.shape[0] != len(train_rows) or query_counts.shape[0] != len(query_rows):
        raise GSE274114BaselineFirewallError("trainer matrix/row count differs")
    probability, fit_receipt = _predict_baseline(
        train_counts, query_counts, targets, model_kind=model_kind, seed=seed
    )
    if probability.shape != (len(query_rows),) or not np.all(np.isfinite(probability)) or np.any((probability < 0) | (probability > 1)):
        raise GSE274114BaselineFirewallError("predicted probability is invalid")
    task_id = str(metadata["task_id"])
    task = TASKS[task_id]
    model_id = task["model_prefix"] + model_kind
    run_id = canonical_sha256(
        {"view_sha256": expected_view_sha256, "task_id": task_id, "model_id": model_id, "seed": seed}
    )
    output.mkdir(parents=True)
    table_fields = ("row_hash", "unit_hash", "predicted")
    prediction_rows = [
        {"row_hash": row["row_id"], "unit_hash": row["row_id"], "predicted": format(float(value), ".17g")}
        for row, value in zip(query_rows, probability)
    ]
    write_tsv(output / "predictions.tsv", table_fields, prediction_rows)
    write_tsv(
        output / "prediction_rows.tsv",
        ("row_hash", "unit_hash"),
        ({"row_hash": row["row_hash"], "unit_hash": row["unit_hash"]} for row in prediction_rows),
    )
    table_ref = artifact_ref(
        output / "predictions.tsv",
        root=output,
        role=f"standardized_prediction_table:{task_id}",
        media_type="text/tab-separated-values",
    )
    row_ref = artifact_ref(
        output / "prediction_rows.tsv",
        root=output,
        role=f"prediction_row_ids:{task_id}",
        media_type="text/tab-separated-values",
    )
    join_sha = canonical_sha256(
        {
            "task_id": task_id,
            "dataset_ids": [DATASET_ID],
            "split_id": task["split_id"],
            "row_id_field": "row_hash",
            "unit_id_field": "unit_hash",
            "unit_id_namespace": f"{task_id}:participant:v1",
            "biological_unit": "participant",
        }
    )
    bundle = PredictionBundle(
        schema_version="masld-bench-prediction-bundle-v1",
        bundle_id=f"{task_id}.{model_id}.fold{metadata['outer_fold']}.seed{seed}",
        run_id=run_id,
        task_id=task_id,
        model_id=model_id,
        dataset_ids=(DATASET_ID,),
        split_id=task["split_id"],
        artifacts=(table_ref, row_ref),
        standardized_table=table_ref,
        row_ids=row_ref,
        n_predictions=len(query_rows),
        row_id_field="row_hash",
        unit_id_field="unit_hash",
        unit_id_namespace=f"{task_id}:participant:v1",
        biological_unit="participant",
        table_schema_sha256=canonical_sha256(list(table_fields)),
        source_join_key_sha256=join_sha,
        format_version="gse274114_binary_probability_v1",
        missing_state="observed",
        metadata={
            "outer_fold": metadata["outer_fold"],
            "seed": seed,
            "model_kind": model_kind,
            "query_labels_opened": False,
            "evaluator_files_opened": False,
            "scoring_run": False,
            "preprocessing_fit_on_training_participants_only": True,
            "cohort_wide_normalization_run": False,
            "effective_length_retained_in_view": True,
            "effective_length_used_by_this_log_cpm_baseline": False,
            **fit_receipt,
        },
    )
    write_json_exclusive(output / "prediction_bundle.json", bundle.to_dict())
    freeze_tree(
        output,
        {
            "artifact_class": "gse274114_frozen_prediction_bundle",
            "task_id": task_id,
            "model_id": model_id,
            "outer_fold": metadata["outer_fold"],
            "seed": seed,
            "query_labels_opened": False,
            "scoring_run": False,
        },
    )
    return bundle.to_dict()


def verify_complete_prediction_set(
    *, bundle_index: Path, model_kind: str, task_id: str
) -> tuple[
    list[dict[str, Any]], dict[str, dict[int, list[tuple[int, float]]]]
]:
    """Verify all bundles before an evaluator is permitted to open labels."""

    if task_id not in TASKS or model_kind not in MODEL_KINDS:
        raise GSE274114BaselineFirewallError("invalid evaluation identity")
    fields, rows = read_tsv(bundle_index)
    if fields != ("outer_fold", "seed", "bundle_root", "artifacts_sha256"):
        raise GSE274114BaselineFirewallError("bundle index schema differs")
    seeds = FIXED_STOCHASTIC_SEEDS if model_kind in STOCHASTIC_MODELS else (0,)
    expected = {(fold, seed) for fold in range(5) for seed in seeds}
    observed = {(int(row["outer_fold"]), int(row["seed"])) for row in rows}
    if observed != expected or len(rows) != len(expected):
        raise GSE274114BaselineFirewallError("prediction set is incomplete")
    task = TASKS[task_id]
    model_id = task["model_prefix"] + model_kind
    verified: list[dict[str, Any]] = []
    probabilities: dict[str, dict[int, list[tuple[int, float]]]] = {}
    for row in sorted(rows, key=lambda item: (int(item["outer_fold"]), int(item["seed"]))):
        root = Path(row["bundle_root"])
        if sha256_file(root / "ARTIFACTS.json") != row["artifacts_sha256"]:
            raise GSE274114BaselineFirewallError("prediction manifest identity differs")
        manifest = verify_frozen_tree(root)
        if manifest["metadata"].get("scoring_run") is not False:
            raise GSE274114BaselineFirewallError("prediction bundle contains scoring")
        bundle = PredictionBundle.load_json(root / "prediction_bundle.json")
        bundle.validate_artifacts(root)
        if bundle.task_id != task_id or bundle.model_id != model_id:
            raise GSE274114BaselineFirewallError("prediction bundle identity differs")
        prediction_fields, prediction_rows = read_tsv(root / bundle.standardized_table.path)
        if prediction_fields != ("row_hash", "unit_hash", "predicted") or any(
            "label" in field.lower() or "observed" in field.lower() for field in prediction_fields
        ):
            raise GSE274114BaselineFirewallError("prediction table exposes outcomes")
        seed = int(row["seed"])
        outer_fold = int(row["outer_fold"])
        for prediction in prediction_rows:
            if prediction["row_hash"] != prediction["unit_hash"]:
                raise GSE274114BaselineFirewallError("participant unit identity differs")
            probability = float(prediction["predicted"])
            if not math.isfinite(probability) or not 0 <= probability <= 1:
                raise GSE274114BaselineFirewallError("prediction is invalid")
            probabilities.setdefault(prediction["row_hash"], {}).setdefault(
                seed, []
            ).append((outer_fold, probability))
        verified.append(bundle.to_dict())
    if any(
        set(seed_records) != set(seeds)
        or any(len(values) != 1 for values in seed_records.values())
        for seed_records in probabilities.values()
    ):
        raise GSE274114BaselineFirewallError(
            "participant predictions are incomplete or duplicated across folds"
        )
    return verified, probabilities


def evaluate_model(
    *,
    bundle_index: Path,
    labels_path: Path,
    folds_path: Path,
    model_kind: str,
    task_id: str,
    bootstrap_replicates: int,
    output: Path,
) -> dict[str, Any]:
    """Evaluator entry point. Bundle completeness is checked before labels open."""

    import numpy as np
    from sklearn.metrics import (
        average_precision_score,
        brier_score_loss,
        f1_score,
        log_loss,
        roc_auc_score,
    )

    if output.exists() or bootstrap_replicates < 1:
        raise GSE274114BaselineFirewallError("invalid or existing evaluator target")
    bundles, by_row = verify_complete_prediction_set(
        bundle_index=bundle_index, model_kind=model_kind, task_id=task_id
    )
    label_fields, labels = read_tsv(labels_path)
    fold_fields, folds = read_tsv(folds_path)
    if label_fields[:2] != ("row_id", "source_group") or fold_fields[:3] != (
        "row_id", "outer_fold", "source_group"
    ):
        raise GSE274114BaselineFirewallError("evaluator label/fold schema differs")
    task = TASKS[task_id]
    target_by_row = {
        row["row_id"]: task["groups"][row["source_group"]]
        for row in labels
        if row["source_group"] in task["groups"]
    }
    if set(by_row) != set(target_by_row):
        raise GSE274114BaselineFirewallError("prediction/outcome row inventory differs")
    fold_by_row = {
        row["row_id"]: int(row["outer_fold"])
        for row in folds
        if row["source_group"] in task["groups"]
    }
    if set(fold_by_row) != set(target_by_row) or any(
        fold != fold_by_row[row]
        for row, seed_records in by_row.items()
        for fold, _ in (values[0] for values in seed_records.values())
    ):
        raise GSE274114BaselineFirewallError(
            "prediction rows do not follow evaluator-held participant folds"
        )
    seeds = FIXED_STOCHASTIC_SEEDS if model_kind in STOCHASTIC_MODELS else (0,)
    rows = sorted(target_by_row)
    truth = np.asarray([target_by_row[row] for row in rows], dtype=np.int64)
    probability = np.asarray(
        [np.mean([by_row[row][seed][0][1] for seed in seeds]) for row in rows],
        dtype=np.float64,
    )
    predicted = (probability >= 0.5).astype(np.int64)

    def metrics(indices: Any) -> dict[str, float]:
        y = truth[indices]
        p = probability[indices]
        hard = (p >= 0.5).astype(np.int64)
        return {
            "participant_macro_f1": float(f1_score(y, hard, average="macro")),
            "participant_auprc": float(average_precision_score(y, p)),
            "participant_auroc": float(roc_auc_score(y, p)),
            "participant_brier": float(brier_score_loss(y, p)),
            "participant_log_loss": float(log_loss(y, np.clip(p, 1e-8, 1 - 1e-8), labels=[0, 1])),
        }

    observed = metrics(np.arange(len(rows)))
    rng = np.random.default_rng(274_114)
    class_indices = [np.flatnonzero(truth == value) for value in (0, 1)]
    bootstrap = {key: [] for key in observed}
    for _ in range(bootstrap_replicates):
        indices = np.concatenate(
            [rng.choice(values, size=len(values), replace=True) for values in class_indices]
        )
        current = metrics(indices)
        for key, value in current.items():
            bootstrap[key].append(value)
    intervals = {
        key: {
            "lower_2_5": float(np.quantile(values, 0.025)),
            "upper_97_5": float(np.quantile(values, 0.975)),
        }
        for key, values in bootstrap.items()
    }
    receipt = {
        "schema_version": "masld-bench-gse274114-within-platform-evaluation-v1",
        "status": "development_only_scored_after_complete_prediction_freeze",
        "task_id": task_id,
        "model_id": task["model_prefix"] + model_kind,
        "participants": len(rows),
        "biological_unit": "participant",
        "prediction_bundles_verified_before_labels_opened": True,
        "prediction_bundle_count": len(bundles),
        "seeds_ensembled": list(seeds),
        "threshold": 0.5,
        "metrics": observed,
        "participant_stratified_bootstrap_replicates": bootstrap_replicates,
        "intervals": intervals,
        "global_four_class_scoring_run": False,
        "global_ood_scoring_run": False,
        "mash_vs_non_mash_scoring_run": False,
        "sealed_or_champion_claim_eligible": False,
        "diagnostic_or_clinical_claim_eligible": False,
    }
    output.mkdir(parents=True)
    write_tsv(
        output / "joined_evaluator_rows.tsv",
        ("row_hash", "observed", "predicted_probability", "predicted_class"),
        (
            {
                "row_hash": row,
                "observed": int(y),
                "predicted_probability": format(float(p), ".17g"),
                "predicted_class": int(hard),
            }
            for row, y, p, hard in zip(rows, truth, probability, predicted)
        ),
    )
    write_json_exclusive(output / "receipt.json", receipt)
    freeze_tree(output, {"artifact_class": "gse274114_development_evaluation", **receipt})
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="action", required=True)
    stage = subparsers.add_parser("stage-features")
    stage.add_argument("--fixture-root", type=Path, required=True)
    stage.add_argument("--output", type=Path, required=True)
    view = subparsers.add_parser("prepare-view")
    view.add_argument("--feature-root", type=Path, required=True)
    view.add_argument("--labels", type=Path, required=True)
    view.add_argument("--folds", type=Path, required=True)
    view.add_argument("--task-id", choices=sorted(TASKS), required=True)
    view.add_argument("--outer-fold", type=int, required=True)
    view.add_argument("--output", type=Path, required=True)
    fit = subparsers.add_parser("fit-predict")
    fit.add_argument("--view-root", type=Path, required=True)
    fit.add_argument("--view-sha256", required=True)
    fit.add_argument("--model-kind", choices=MODEL_KINDS, required=True)
    fit.add_argument("--seed", type=int, required=True)
    fit.add_argument("--output", type=Path, required=True)
    evaluate = subparsers.add_parser("evaluate")
    evaluate.add_argument("--bundle-index", type=Path, required=True)
    evaluate.add_argument("--labels", type=Path, required=True)
    evaluate.add_argument("--folds", type=Path, required=True)
    evaluate.add_argument("--model-kind", choices=MODEL_KINDS, required=True)
    evaluate.add_argument("--task-id", choices=sorted(TASKS), required=True)
    evaluate.add_argument("--bootstrap-replicates", type=int, default=10_000)
    evaluate.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.action == "stage-features":
        result = stage_unlabeled_features(fixture_root=args.fixture_root, output=args.output)
    elif args.action == "prepare-view":
        result = prepare_fold_view(
            feature_root=args.feature_root,
            labels_path=args.labels,
            folds_path=args.folds,
            task_id=args.task_id,
            outer_fold=args.outer_fold,
            output=args.output,
        )
    elif args.action == "fit-predict":
        result = fit_predict(
            view_root=args.view_root,
            expected_view_sha256=args.view_sha256,
            model_kind=args.model_kind,
            seed=args.seed,
            output=args.output,
        )
    else:
        result = evaluate_model(
            bundle_index=args.bundle_index,
            labels_path=args.labels,
            folds_path=args.folds,
            model_kind=args.model_kind,
            task_id=args.task_id,
            bootstrap_replicates=args.bootstrap_replicates,
            output=args.output,
        )
    print(json.dumps(result, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()

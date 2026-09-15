#!/usr/bin/env python3
"""Fit fold-safe common heads for GSE281364 DNA-language-model features."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Ridge


MODELS = ("dnabert2", "hyenadna", "nucleotide_transformer")
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
MODEL_ELIGIBILITY = {
    "dnabert2": ("Apache-2.0", False, True),
    "hyenadna": ("BSD-3-Clause", False, True),
    "nucleotide_transformer": ("CC-BY-NC-SA-4.0", True, False),
}
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
SEED = 20260824
PSEUDOCOUNT = 0.5
BOOTSTRAP_REPLICATES = 1000
RIDGE_ALPHAS = (
    0.1,
    1.0,
    10.0,
    100.0,
    1000.0,
    10000.0,
    100000.0,
    1000000.0,
    math.inf,
)
MLP_MAX_EPOCHS = 120
MLP_PATIENCE = 15
MLP_LEARNING_RATE = 1.0e-3
MLP_WEIGHT_DECAY = 1.0e-3


class CommonLaneError(RuntimeError):
    """Raised when a common-lane data or evaluation requirement is not met."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise CommonLaneError(f"JSON object required: {path}")
    return value


def verify_frozen_tree(root: Path, expected_sha256: str) -> dict[str, Any]:
    if not root.is_dir() or root.is_symlink():
        raise CommonLaneError(f"invalid frozen root: {root}")
    manifest_path = root / "ARTIFACTS.json"
    complete_path = root / "COMPLETE"
    if not manifest_path.is_file() or not complete_path.is_file():
        raise CommonLaneError(f"incomplete frozen root: {root}")
    if manifest_path.is_symlink() or complete_path.is_symlink():
        raise CommonLaneError(f"symlinked control file: {root}")
    if sha256_file(manifest_path) != expected_sha256:
        raise CommonLaneError(f"manifest checksum differs: {root}")
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
        raise CommonLaneError(f"frozen control contract differs: {root}")
    expected: dict[str, tuple[str, int]] = {}
    for item in artifacts:
        if not isinstance(item, dict) or set(item) != {"path", "sha256", "size_bytes"}:
            raise CommonLaneError("artifact member schema differs")
        relative = Path(str(item["path"]))
        key = relative.as_posix()
        if relative.is_absolute() or ".." in relative.parts or key in expected:
            raise CommonLaneError("unsafe or duplicate artifact path")
        expected[key] = (str(item["sha256"]), int(item["size_bytes"]))
    observed: set[str] = set()
    for path in root.rglob("*"):
        if path.is_symlink():
            raise CommonLaneError(f"symlink in frozen tree: {path}")
        relative = path.relative_to(root).as_posix()
        if path.is_file() and relative not in {"ARTIFACTS.json", "COMPLETE"}:
            observed.add(relative)
    if observed != set(expected):
        raise CommonLaneError(f"frozen inventory differs: {root}")
    for relative, (digest, size) in expected.items():
        path = root / relative
        if path.stat().st_size != size or sha256_file(path) != digest:
            raise CommonLaneError(f"frozen artifact differs: {relative}")
    return manifest


def array_digest(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value)
    value_hash = sha256()
    value_hash.update(str(array.dtype).encode())
    value_hash.update(str(array.shape).encode())
    value_hash.update(array.tobytes())
    return value_hash.hexdigest()


def activity_delta(
    ref_dna: int,
    ref_rna: int,
    alt_dna: int,
    alt_rna: int,
    pseudocount: float = PSEUDOCOUNT,
) -> float:
    if min(ref_dna, ref_rna, alt_dna, alt_rna) < 0 or pseudocount <= 0:
        raise CommonLaneError("invalid reporter count")
    return math.log2((alt_rna + pseudocount) / (alt_dna + pseudocount)) - math.log2(
        (ref_rna + pseudocount) / (ref_dna + pseudocount)
    )


def load_fixture(root: Path) -> tuple[list[dict[str, str]], np.ndarray]:
    path = root / "fixture/sequence_manifest.tsv"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    if (
        len(rows) != 1033
        or len({row["fixture_id"] for row in rows}) != 1033
        or len({row["element_id"] for row in rows}) != 1033
        or len({row["outer_locus_sequence_group_id"] for row in rows}) != 1033
        or {int(row["outer_fold"]) for row in rows} != set(range(5))
    ):
        raise CommonLaneError("common fixture census differs")
    bases = {base: index for index, base in enumerate("ACGT")}
    allele_features = np.zeros((len(rows), 16), dtype=np.float64)
    for index, row in enumerate(rows):
        ref, alt = row["ref"], row["alt"]
        if ref not in bases or alt not in bases or ref == alt:
            raise CommonLaneError("fixture allele identity differs")
        allele_features[index, 4 * bases[ref] + bases[alt]] = 1.0
    return rows, allele_features


def load_outcomes(
    path: Path, eligible_ids: set[str]
) -> dict[tuple[str, str], dict[str, Any]]:
    pairs: dict[tuple[str, str, int], dict[str, tuple[int, int]]] = defaultdict(dict)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != OUTCOME_FIELDS:
            raise CommonLaneError("outcome table schema differs")
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
                raise CommonLaneError("outcome topology differs")
            if row["assay_state"] != "observed":
                continue
            key = (row["element_id"], row["context_id"], replicate)
            if allele in pairs[key]:
                raise CommonLaneError("duplicate observed outcome allele")
            pairs[key][allele] = (int(row["DNA"]), int(row["RNA"]))
    output: dict[tuple[str, str], dict[str, Any]] = {}
    for element in sorted(eligible_ids):
        for context in CONTEXTS:
            replicate_values: list[float] = []
            for replicate in range(1, 5):
                pair = pairs.get((element, context, replicate), {})
                if set(pair) != {"ref", "alt"}:
                    raise CommonLaneError("four complete experimental replicate pairs required")
                replicate_values.append(activity_delta(*pair["ref"], *pair["alt"]))
            output[(element, context)] = {
                "mean": float(np.mean(replicate_values)),
                "sd": float(np.std(replicate_values, ddof=1)),
                "replicates": tuple(replicate_values),
            }
    return output


def _projected_features(raw: np.ndarray, parameters: Mapping[str, np.ndarray]) -> np.ndarray:
    mean = np.asarray(parameters["mean"], dtype=np.float64)
    std = np.asarray(parameters["standard_deviation"], dtype=np.float64)
    components = np.asarray(parameters["components"], dtype=np.float64)
    whitening = np.asarray(parameters["whitening_scale"], dtype=np.float64)
    projected = ((raw.astype(np.float64) - mean) / std) @ components.T
    projected = (projected / whitening).astype(np.float32)
    ref = (projected[:, 0] + projected[:, 2]) / 2.0
    alt = (projected[:, 1] + projected[:, 3]) / 2.0
    delta = alt - ref
    return np.concatenate((ref, alt, delta, np.abs(delta)), axis=1)


def audit_and_load_features(
    root: Path, fixture_rows: Sequence[Mapping[str, str]]
) -> tuple[dict[tuple[str, int], np.ndarray], dict[str, Any]]:
    receipt = _load_json(root / "projected/receipt.json")
    if (
        receipt.get("status") != "pass_outcome_blind_fold_safe_projection"
        or receipt.get("elements") != 1033
        or receipt.get("outer_locus_sequence_groups") != 1033
        or receipt.get("outer_folds") != 5
        or receipt.get("projection_width") != 256
        or receipt.get("head_input_width") != 1024
        or receipt.get("projection_fit_on_held_out_fold")
        or receipt.get("reporter_outcomes_read")
        or receipt.get("sealed_labels_read")
        or receipt.get("downstream_head_fit")
    ):
        raise CommonLaneError("projection receipt differs")
    expected_ids = np.asarray([row["fixture_id"] for row in fixture_rows])
    expected_groups = np.asarray(
        [row["outer_locus_sequence_group_id"] for row in fixture_rows]
    )
    expected_folds = np.asarray([int(row["outer_fold"]) for row in fixture_rows])
    features: dict[tuple[str, int], np.ndarray] = {}
    checks: list[dict[str, Any]] = []
    for model in MODELS:
        license_name, restricted, open_eligible = MODEL_ELIGIBILITY[model]
        raw_receipt = _load_json(root / f"raw/{model}/receipt.json")
        if (
            raw_receipt.get("status") != "pass_outcome_blind_embedding_extraction"
            or raw_receipt.get("model_id") != model
            or raw_receipt.get("elements") != 1033
            or raw_receipt.get("license") != license_name
            or bool(raw_receipt.get("restricted_comparator")) != restricted
            or bool(raw_receipt.get("open_champion_eligible_after_task_gates"))
            != open_eligible
            or raw_receipt.get("reporter_outcomes_read")
            or raw_receipt.get("sealed_labels_read")
            or raw_receipt.get("downstream_head_fit")
        ):
            raise CommonLaneError(f"{model} raw receipt differs")
        with np.load(root / f"raw/{model}/allele_embeddings.npz", allow_pickle=False) as data:
            ids = data["fixture_ids"]
            groups = data["outer_locus_sequence_group_ids"]
            folds = data["outer_folds"].astype(np.int64)
            raw = data["embeddings"]
            alleles = tuple(data["allele_order"].tolist())
        if (
            not np.array_equal(ids, expected_ids)
            or not np.array_equal(groups, expected_groups)
            or not np.array_equal(folds, expected_folds)
            or alleles != ("REF", "ALT", "REF_RC", "ALT_RC")
            or raw.shape[0:2] != (1033, 4)
            or not np.isfinite(raw).all()
        ):
            raise CommonLaneError(f"{model} raw alignment differs")
        model_receipt = receipt["models"][model]
        if (
            model_receipt.get("license") != license_name
            or bool(model_receipt.get("restricted_comparator")) != restricted
        ):
            raise CommonLaneError(f"{model} projected eligibility differs")
        for held_fold in range(5):
            fold_root = root / f"projected/{model}/heldout_fold{held_fold}"
            with np.load(fold_root / "projection_parameters.npz", allow_pickle=False) as data:
                parameters = {key: data[key] for key in data.files}
            with np.load(fold_root / "head_features.npz", allow_pickle=False) as data:
                fold_ids = data["fixture_ids"]
                fold_groups = data["outer_locus_sequence_group_ids"]
                fold_values = data["outer_folds"].astype(np.int64)
                block_order = tuple(data["feature_block_order"].tolist())
                saved_features = data["features"]
            train_raw = raw[folds != held_fold].reshape(-1, raw.shape[2]).astype(np.float64)
            expected_mean = train_raw.mean(axis=0)
            expected_std = np.maximum(train_raw.std(axis=0), 1.0e-6)
            components = np.asarray(parameters["components"], dtype=np.float64)
            if (
                int(parameters["held_out_fold"]) != held_fold
                or int(parameters["training_elements"]) != int(np.sum(folds != held_fold))
                or not np.array_equal(fold_ids, expected_ids)
                or not np.array_equal(fold_groups, expected_groups)
                or not np.array_equal(fold_values, expected_folds)
                or block_order != ("REF", "ALT", "ALT_minus_REF", "absolute_ALT_minus_REF")
                or saved_features.shape != (1033, 1024)
                or not np.allclose(parameters["mean"], expected_mean, rtol=0.0, atol=1e-12)
                or not np.allclose(parameters["standard_deviation"], expected_std, rtol=0.0, atol=1e-12)
                or components.shape != (256, raw.shape[2])
                or not np.allclose(components @ components.T, np.eye(256), atol=2e-5)
                or any(
                    row[int(np.argmax(np.abs(row)))] < 0
                    for row in components
                )
            ):
                raise CommonLaneError(f"{model} fold {held_fold} projection contract differs")
            reconstructed = _projected_features(raw, parameters)
            fold_receipt = model_receipt["folds"][held_fold]
            if (
                not np.allclose(reconstructed, saved_features, rtol=2e-6, atol=2e-6)
                or array_digest(saved_features) != fold_receipt["features_sha256"]
                or array_digest(components)
                != fold_receipt["projection_components_sha256"]
            ):
                raise CommonLaneError(f"{model} fold {held_fold} feature reconstruction differs")
            features[(model, held_fold)] = saved_features.astype(np.float64)
            checks.append(
                {
                    "model_id": model,
                    "held_out_fold": held_fold,
                    "training_elements": int(np.sum(folds != held_fold)),
                    "held_out_elements": int(np.sum(folds == held_fold)),
                    "training_only_mean_std_rederived": True,
                    "head_features_reconstructed": True,
                    "orthonormal_projection": True,
                }
            )
    return features, {"checks": checks, "status": "pass_independent_projection_consistency_audit"}


def _standardize_fit(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = values.mean(axis=0)
    scale = np.maximum(values.std(axis=0), 1.0e-6)
    return mean, scale


def _inner_validation_fold(held_fold: int) -> int:
    return (held_fold + 1) % 5


def fit_ridge_outer(
    features: np.ndarray,
    outcomes: np.ndarray,
    folds: np.ndarray,
    held_fold: int,
    *,
    alphas: Sequence[float] = RIDGE_ALPHAS,
) -> tuple[np.ndarray, dict[str, Any], dict[str, np.ndarray]]:
    outer_train = folds != held_fold
    outer_test = folds == held_fold
    inner_fold = _inner_validation_fold(held_fold)
    inner_train = outer_train & (folds != inner_fold)
    inner_valid = folds == inner_fold
    if np.any(inner_train & outer_test) or np.any(inner_valid & outer_test):
        raise CommonLaneError("held fold entered ridge fitting")
    mean, scale = _standardize_fit(features[inner_train])
    x_train = (features[inner_train] - mean) / scale
    x_valid = (features[inner_valid] - mean) / scale
    validation: list[tuple[float, float]] = []
    for alpha in alphas:
        if math.isinf(alpha):
            prediction = np.full(np.sum(inner_valid), outcomes[inner_train].mean())
        else:
            model = Ridge(alpha=float(alpha), fit_intercept=True, solver="lsqr", tol=1.0e-7)
            model.fit(x_train, outcomes[inner_train])
            prediction = model.predict(x_valid)
        rmse = float(np.sqrt(np.mean((outcomes[inner_valid] - prediction) ** 2)))
        validation.append((float(alpha), rmse))
    selected_alpha = min(validation, key=lambda value: (value[1], value[0]))[0]
    full_mean, full_scale = _standardize_fit(features[outer_train])
    if math.isinf(selected_alpha):
        coefficient = np.zeros(features.shape[1], dtype=np.float64)
        intercept = float(outcomes[outer_train].mean())
        prediction = np.full(np.sum(outer_test), intercept)
    else:
        final = Ridge(alpha=selected_alpha, fit_intercept=True, solver="lsqr", tol=1.0e-7)
        final.fit((features[outer_train] - full_mean) / full_scale, outcomes[outer_train])
        prediction = final.predict((features[outer_test] - full_mean) / full_scale)
        coefficient = np.asarray(final.coef_, dtype=np.float64)
        intercept = float(final.intercept_)
    receipt = {
        "held_out_fold": held_fold,
        "outer_training_folds": sorted(set(folds[outer_train].tolist())),
        "inner_validation_fold": inner_fold,
        "selected_alpha": selected_alpha,
        "inner_validation_rmse": {format(alpha, ".17g"): rmse for alpha, rmse in validation},
        "held_out_outcomes_read_during_fit": False,
    }
    state = {
        "feature_mean": full_mean,
        "feature_scale": full_scale,
        "coefficient": coefficient,
        "intercept": np.asarray(intercept, dtype=np.float64),
        "selected_alpha": np.asarray(selected_alpha, dtype=np.float64),
    }
    return np.asarray(prediction, dtype=np.float64), receipt, state


def fit_mlp_outer(
    features: np.ndarray,
    outcomes: np.ndarray,
    folds: np.ndarray,
    held_fold: int,
    *,
    seed: int,
    hidden_width: int = 256,
    max_epochs: int = MLP_MAX_EPOCHS,
    patience: int = MLP_PATIENCE,
) -> tuple[np.ndarray, dict[str, Any], dict[str, np.ndarray]]:
    import torch
    from torch import nn

    class Head(nn.Module):
        def __init__(self, width: int) -> None:
            super().__init__()
            self.layer1 = nn.Linear(width, hidden_width)
            self.activation = nn.GELU()
            self.layer2 = nn.Linear(hidden_width, 1)

        def forward(self, value: Any) -> Any:
            return self.layer2(self.activation(self.layer1(value))).squeeze(1)

    outer_train = folds != held_fold
    outer_test = folds == held_fold
    inner_fold = _inner_validation_fold(held_fold)
    inner_train = outer_train & (folds != inner_fold)
    inner_valid = folds == inner_fold
    if np.any(inner_train & outer_test) or np.any(inner_valid & outer_test):
        raise CommonLaneError("held fold entered MLP fitting")

    def selection_tensors(
        fit_mask: np.ndarray, validation_mask: np.ndarray
    ) -> tuple[Any, Any, Any, Any]:
        feature_mean, feature_scale = _standardize_fit(features[fit_mask])
        outcome_mean = float(outcomes[fit_mask].mean())
        outcome_scale = max(float(outcomes[fit_mask].std()), 1.0e-6)
        return (
            torch.from_numpy(((features[fit_mask] - feature_mean) / feature_scale).astype(np.float32)),
            torch.from_numpy(((outcomes[fit_mask] - outcome_mean) / outcome_scale).astype(np.float32)),
            torch.from_numpy(
                ((features[validation_mask] - feature_mean) / feature_scale).astype(np.float32)
            ),
            torch.from_numpy(
                ((outcomes[validation_mask] - outcome_mean) / outcome_scale).astype(np.float32)
            ),
        )

    torch.manual_seed(seed)
    x_train, y_train, x_valid, y_valid = selection_tensors(inner_train, inner_valid)
    selector = Head(features.shape[1])
    optimizer = torch.optim.AdamW(
        selector.parameters(), lr=MLP_LEARNING_RATE, weight_decay=MLP_WEIGHT_DECAY
    )
    best_epoch = 1
    best_loss = math.inf
    stale = 0
    for epoch in range(1, max_epochs + 1):
        selector.train()
        optimizer.zero_grad(set_to_none=True)
        loss = torch.mean((selector(x_train) - y_train) ** 2)
        loss.backward()
        optimizer.step()
        selector.eval()
        with torch.no_grad():
            valid_loss = float(torch.mean((selector(x_valid) - y_valid) ** 2).item())
        if valid_loss < best_loss - 1.0e-7:
            best_loss = valid_loss
            best_epoch = epoch
            stale = 0
        else:
            stale += 1
        if stale >= patience:
            break

    torch.manual_seed(seed)
    mean, scale = _standardize_fit(features[outer_train])
    outcome_mean = float(outcomes[outer_train].mean())
    outcome_scale = max(float(outcomes[outer_train].std()), 1.0e-6)
    x_full = torch.from_numpy(
        ((features[outer_train] - mean) / scale).astype(np.float32)
    )
    y_full = torch.from_numpy(
        ((outcomes[outer_train] - outcome_mean) / outcome_scale).astype(np.float32)
    )
    # Held-out features are transformed for inference. Held-out outcomes are not
    # indexed anywhere in model selection or final fitting.
    x_test = torch.from_numpy(
        ((features[outer_test] - mean) / scale).astype(np.float32)
    )
    final = Head(features.shape[1])
    optimizer = torch.optim.AdamW(
        final.parameters(), lr=MLP_LEARNING_RATE, weight_decay=MLP_WEIGHT_DECAY
    )
    for _ in range(best_epoch):
        final.train()
        optimizer.zero_grad(set_to_none=True)
        loss = torch.mean((final(x_full) - y_full) ** 2)
        loss.backward()
        optimizer.step()
    final.eval()
    with torch.no_grad():
        prediction = final(x_test).numpy().astype(np.float64) * outcome_scale + outcome_mean
    state: dict[str, np.ndarray] = {
        "feature_mean": mean,
        "feature_scale": scale,
        "outcome_mean": np.asarray(outcome_mean, dtype=np.float64),
        "outcome_scale": np.asarray(outcome_scale, dtype=np.float64),
        "selected_epochs": np.asarray(best_epoch, dtype=np.int64),
    }
    for name, value in final.state_dict().items():
        state[name.replace(".", "__")] = value.detach().cpu().numpy()
    receipt = {
        "held_out_fold": held_fold,
        "outer_training_folds": sorted(set(folds[outer_train].tolist())),
        "inner_validation_fold": inner_fold,
        "selected_epochs": best_epoch,
        "best_inner_standardized_mse": best_loss,
        "max_epochs": max_epochs,
        "patience": patience,
        "seed": seed,
        "held_out_outcomes_read_during_fit": False,
    }
    return prediction, receipt, state


def metric_values(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float | None]:
    if observed.shape != predicted.shape or observed.ndim != 1 or observed.size < 3:
        raise CommonLaneError("metric inputs differ")
    residual = observed - predicted
    result: dict[str, float | None] = {
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "mae": float(np.mean(np.abs(residual))),
        "r2": float(1.0 - np.sum(residual**2) / np.sum((observed - observed.mean()) ** 2)),
        "pearson": None,
        "spearman": None,
        "calibration_intercept": None,
        "calibration_slope": None,
    }
    if not np.all(predicted == predicted[0]):
        result["pearson"] = float(pearsonr(observed, predicted).statistic)
        result["spearman"] = float(spearmanr(observed, predicted).statistic)
        slope, intercept = np.polyfit(predicted, observed, 1)
        result["calibration_intercept"] = float(intercept)
        result["calibration_slope"] = float(slope)
    return result


def bootstrap_metrics(
    observed: np.ndarray,
    predicted: np.ndarray,
    groups: np.ndarray,
    *,
    seed: int,
    replicates: int = BOOTSTRAP_REPLICATES,
) -> dict[str, float | int | str]:
    unique = np.unique(groups)
    if len(unique) != len(groups):
        raise CommonLaneError("each smoke outcome must have one locus group")
    rng = np.random.default_rng(seed)
    values: dict[str, list[float]] = {"rmse": [], "pearson": [], "spearman": []}
    constant = np.all(predicted == predicted[0])
    for _ in range(replicates):
        indices = rng.integers(0, len(unique), size=len(unique))
        result = metric_values(observed[indices], predicted[indices])
        values["rmse"].append(float(result["rmse"]))
        if not constant and result["pearson"] is not None and result["spearman"] is not None:
            values["pearson"].append(float(result["pearson"]))
            values["spearman"].append(float(result["spearman"]))
    output: dict[str, float | int | str] = {"valid_bootstrap_replicates": replicates}
    for name in ("rmse", "pearson", "spearman"):
        if not values[name]:
            output[f"{name}_ci_low"] = "not_applicable_constant_prediction"
            output[f"{name}_ci_high"] = "not_applicable_constant_prediction"
        else:
            output[f"{name}_ci_low"] = float(np.quantile(values[name], 0.025))
            output[f"{name}_ci_high"] = float(np.quantile(values[name], 0.975))
    return output


def _write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise CommonLaneError(f"cannot write empty table: {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_gzip_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise CommonLaneError(f"cannot write empty table: {path}")
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            import io

            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)


def _value(value: float | None) -> str:
    return "not_applicable_constant_prediction" if value is None else format(value, ".17g")


def evaluate(arguments: argparse.Namespace) -> dict[str, Any]:
    if arguments.output.exists():
        raise CommonLaneError("evaluation output exists")
    embedding_manifest = verify_frozen_tree(arguments.embeddings, arguments.embeddings_sha256)
    outcome_manifest = verify_frozen_tree(arguments.outcomes, arguments.outcomes_sha256)
    fixture_manifest = verify_frozen_tree(arguments.fixture, arguments.fixture_sha256)
    if (
        embedding_manifest["metadata"].get("artifact_class")
        != "outcome_blind_common_representation_embeddings"
        or embedding_manifest["metadata"].get("reporter_outcomes_read")
        or embedding_manifest["metadata"].get("downstream_head_fit")
        or outcome_manifest["metadata"].get("artifact_class")
        != "gse281364_replicate_safe_outcomes"
        or outcome_manifest["metadata"].get("donor_count") != 0
        or outcome_manifest["metadata"].get("outcome_role")
        != "exposed_development_MPRA_only"
        or outcome_manifest["metadata"].get("sealed_outcomes_loaded")
        or fixture_manifest["metadata"].get("artifact_class")
        != "gse281364_outcome_blind_dna_lm_common_fixture"
        or fixture_manifest["metadata"].get("reporter_outcomes_read")
    ):
        raise CommonLaneError("input role or firewall differs")

    fixture_rows, allele_features = load_fixture(arguments.fixture)
    projected, projection_audit = audit_and_load_features(arguments.embeddings, fixture_rows)
    ids = np.asarray([row["element_id"] for row in fixture_rows])
    groups = np.asarray([row["outer_locus_sequence_group_id"] for row in fixture_rows])
    folds = np.asarray([int(row["outer_fold"]) for row in fixture_rows], dtype=np.int64)
    outcomes = load_outcomes(
        arguments.outcomes / "outcomes/replicate_outcomes.tsv.gz", set(ids.tolist())
    )
    target = {
        context: np.asarray([outcomes[(element, context)]["mean"] for element in ids])
        for context in CONTEXTS
    }
    target_sd = {
        context: np.asarray([outcomes[(element, context)]["sd"] for element in ids])
        for context in CONTEXTS
    }

    arguments.output.mkdir(parents=True, mode=0o750)
    (arguments.output / "heads").mkdir()
    prediction_rows: list[dict[str, Any]] = []
    selection_rows: list[dict[str, Any]] = []
    predictions: dict[tuple[str, str, str], np.ndarray] = {}

    for context_index, context in enumerate(CONTEXTS):
        y = target[context]
        for baseline_id in ("zero", "training_mean"):
            oof = np.empty(len(ids), dtype=np.float64)
            for held_fold in range(5):
                test = folds == held_fold
                train = ~test
                oof[test] = 0.0 if baseline_id == "zero" else float(y[train].mean())
            predictions[("task_native_baseline", baseline_id, context)] = oof
        allele_oof = np.empty(len(ids), dtype=np.float64)
        allele_root = arguments.output / f"heads/task_native_baseline/{context}"
        allele_root.mkdir(parents=True)
        for held_fold in range(5):
            test = folds == held_fold
            prediction, selection, state = fit_ridge_outer(
                allele_features, y, folds, held_fold
            )
            allele_oof[test] = prediction
            np.savez_compressed(allele_root / f"fold{held_fold}_allele_ridge.npz", **state)
            selection_rows.append(
                {
                    "model_id": "task_native_baseline",
                    "head_id": "allele_identity_ridge",
                    "context_id": context,
                    **selection,
                }
            )
        predictions[("task_native_baseline", "allele_identity_ridge", context)] = allele_oof

        for model_index, model_id in enumerate(MODELS):
            model_root = arguments.output / f"heads/{model_id}/{context}"
            model_root.mkdir(parents=True)
            for head_id in ("linear_ridge", "two_layer_gelu"):
                oof = np.empty(len(ids), dtype=np.float64)
                for held_fold in range(5):
                    test = folds == held_fold
                    feature_values = projected[(model_id, held_fold)]
                    if head_id == "linear_ridge":
                        prediction, selection, state = fit_ridge_outer(
                            feature_values, y, folds, held_fold
                        )
                    else:
                        prediction, selection, state = fit_mlp_outer(
                            feature_values,
                            y,
                            folds,
                            held_fold,
                            seed=SEED + context_index * 10 + held_fold,
                        )
                    oof[test] = prediction
                    np.savez_compressed(
                        model_root / f"fold{held_fold}_{head_id}.npz", **state
                    )
                    selection_rows.append(
                        {
                            "model_id": model_id,
                            "head_id": head_id,
                            "context_id": context,
                            **selection,
                        }
                    )
                predictions[(model_id, head_id, context)] = oof

    for (model_id, head_id, context), prediction in sorted(predictions.items()):
        for index, element in enumerate(ids):
            prediction_rows.append(
                {
                    "model_id": model_id,
                    "head_id": head_id,
                    "context_id": context,
                    "element_id": element,
                    "outer_locus_sequence_group_id": groups[index],
                    "outer_fold": int(folds[index]),
                    "prediction": format(float(prediction[index]), ".17g"),
                    "observed_mean_signed_log2_activity_delta": format(
                        float(target[context][index]), ".17g"
                    ),
                    "observed_replicate_sd": format(float(target_sd[context][index]), ".17g"),
                    "experimental_replicates": 4,
                    "biological_donors": 0,
                    "outcome_role": "exposed_development_MPRA_only",
                    "champion_eligible": "false",
                }
            )

    aggregate_rows: list[dict[str, Any]] = []
    fold_rows: list[dict[str, Any]] = []
    aggregate_lookup: dict[tuple[str, str, str], dict[str, float | None]] = {}
    for key, prediction in sorted(predictions.items()):
        model_id, head_id, context = key
        observed = target[context]
        metrics = metric_values(observed, prediction)
        aggregate_lookup[key] = metrics
        interval = bootstrap_metrics(
            observed,
            prediction,
            groups,
            seed=SEED + int(sha256("|".join(key).encode()).hexdigest()[:8], 16),
        )
        license_name, restricted, open_eligible = MODEL_ELIGIBILITY.get(
            model_id, ("not_applicable", False, False)
        )
        aggregate_rows.append(
            {
                "model_id": model_id,
                "head_id": head_id,
                "context_id": context,
                "elements": len(ids),
                "locus_sequence_groups": len(set(groups.tolist())),
                "pearson": _value(metrics["pearson"]),
                "pearson_ci_low": interval["pearson_ci_low"],
                "pearson_ci_high": interval["pearson_ci_high"],
                "spearman": _value(metrics["spearman"]),
                "spearman_ci_low": interval["spearman_ci_low"],
                "spearman_ci_high": interval["spearman_ci_high"],
                "rmse": _value(metrics["rmse"]),
                "rmse_ci_low": interval["rmse_ci_low"],
                "rmse_ci_high": interval["rmse_ci_high"],
                "mae": _value(metrics["mae"]),
                "r2": _value(metrics["r2"]),
                "calibration_intercept": _value(metrics["calibration_intercept"]),
                "calibration_slope": _value(metrics["calibration_slope"]),
                "license": license_name,
                "restricted_comparator": str(restricted).lower(),
                "open_champion_eligible_after_task_gates": str(open_eligible).lower(),
                "current_task_champion_eligible": "false",
                "multiple_testing": "none_descriptive_one_seed_smoke",
            }
        )
        for held_fold in range(5):
            selected = folds == held_fold
            fold_metric = metric_values(observed[selected], prediction[selected])
            fold_rows.append(
                {
                    "model_id": model_id,
                    "head_id": head_id,
                    "context_id": context,
                    "outer_fold": held_fold,
                    "elements": int(np.sum(selected)),
                    "locus_sequence_groups": int(np.sum(selected)),
                    **{name: _value(value) for name, value in fold_metric.items()},
                }
            )

    promotion_rows: list[dict[str, Any]] = []
    for model_id in MODELS:
        candidate_heads = []
        for head_id in ("linear_ridge", "two_layer_gelu"):
            context_metrics = [aggregate_lookup[(model_id, head_id, context)] for context in CONTEXTS]
            mean_spearman = float(np.mean([float(value["spearman"]) for value in context_metrics]))
            mean_rmse = float(np.mean([float(value["rmse"]) for value in context_metrics]))
            candidate_heads.append((mean_spearman, -mean_rmse, head_id, context_metrics))
        mean_spearman, _, best_head, context_metrics = max(candidate_heads)
        baseline_spearman = float(
            np.mean(
                [
                    float(aggregate_lookup[("task_native_baseline", "allele_identity_ridge", context)]["spearman"])
                    for context in CONTEXTS
                ]
            )
        )
        positive_both = all(float(value["spearman"]) > 0 for value in context_metrics)
        mean_baseline_better = all(
            float(context_metrics[index]["rmse"])
            < float(aggregate_lookup[("task_native_baseline", "training_mean", context)]["rmse"])
            for index, context in enumerate(CONTEXTS)
        )
        recommended = (
            mean_spearman >= 0.10
            and mean_spearman - baseline_spearman >= 0.02
            and positive_both
            and mean_baseline_better
        )
        license_name, restricted, open_eligible = MODEL_ELIGIBILITY[model_id]
        promotion_rows.append(
            {
                "model_id": model_id,
                "selected_head": best_head,
                "mean_two_context_spearman": format(mean_spearman, ".17g"),
                "allele_baseline_mean_spearman": format(baseline_spearman, ".17g"),
                "spearman_gain": format(mean_spearman - baseline_spearman, ".17g"),
                "positive_spearman_both_contexts": str(positive_both).lower(),
                "rmse_better_than_training_mean_both_contexts": str(mean_baseline_better).lower(),
                "three_seed_screening_recommended": str(recommended).lower(),
                "promotion_gate": (
                    "mean_spearman>=0.10;gain_over_allele_baseline>=0.02;"
                    "positive_both_contexts;rmse_better_than_training_mean_both_contexts"
                ),
                "license": license_name,
                "restricted_comparator": str(restricted).lower(),
                "open_champion_eligible_after_task_gates": str(open_eligible).lower(),
                "current_task_champion_eligible": "false",
            }
        )

    _write_gzip_tsv(arguments.output / "oof_predictions.tsv.gz", prediction_rows)
    _write_tsv(arguments.output / "aggregate_metrics.tsv", aggregate_rows)
    _write_tsv(arguments.output / "outer_fold_metrics.tsv", fold_rows)
    _write_tsv(arguments.output / "head_selection.tsv", selection_rows)
    _write_tsv(arguments.output / "promotion.tsv", promotion_rows)
    audit = {
        "schema_version": "masld-bench-gse281364-dna-lm-input-audit-v1",
        "status": "pass_independent_frozen_and_projection_audit",
        "embedding_artifacts_sha256": arguments.embeddings_sha256,
        "outcome_artifacts_sha256": arguments.outcomes_sha256,
        "fixture_artifacts_sha256": arguments.fixture_sha256,
        "recursive_member_checksums_verified": True,
        "complete_markers_verified": True,
        "training_only_projection_statistics_rederived": True,
        "head_features_reconstructed": True,
        **projection_audit,
    }
    (arguments.output / "input_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt = {
        "schema_version": "masld-bench-gse281364-dna-lm-common-lane-evaluation-v1",
        "status": "pass_descriptive_one_seed_development_smoke",
        "dataset_id": "gse281364",
        "models": list(MODELS),
        "contexts": list(CONTEXTS),
        "endpoint": "mean_signed_log2_ALT_minus_REF_RNA_over_DNA_activity",
        "elements": len(ids),
        "outer_locus_sequence_groups": len(set(groups.tolist())),
        "outer_folds": 5,
        "experimental_replicates_per_context": 4,
        "experimental_replicates_used_as_independent_donors": False,
        "donor_count": 0,
        "outcome_role": "exposed_development_MPRA_only",
        "heads": ["linear_ridge", "two_layer_gelu"],
        "baselines": ["zero", "training_mean", "allele_identity_ridge"],
        "projection_fit_scope": "outer_training_locus_sequence_groups_only",
        "head_fit_scope": "outer_training_locus_sequence_groups_only",
        "hyperparameter_selection_scope": "one_prespecified_inner_outer_fold_within_outer_training_only",
        "held_out_outcomes_read_during_fit": False,
        "bootstrap_unit": "outer_locus_sequence_group",
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "seed": SEED,
        "models_ranked_for_smoke_promotion_only": True,
        "confirmatory_inference": False,
        "multiple_testing": "none_descriptive_one_seed_smoke",
        "sealed_outcomes_read": False,
        "external_evaluation": False,
        "champion_eligible": False,
        "clinical_claim_supported": False,
        "three_seed_screening_recommended": [
            row["model_id"]
            for row in promotion_rows
            if row["three_seed_screening_recommended"] == "true"
        ],
        "restricted_comparator": ["nucleotide_transformer"],
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--embeddings-sha256", required=True)
    parser.add_argument("--outcomes", type=Path, required=True)
    parser.add_argument("--outcomes-sha256", required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--fixture-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    evaluate(arguments)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

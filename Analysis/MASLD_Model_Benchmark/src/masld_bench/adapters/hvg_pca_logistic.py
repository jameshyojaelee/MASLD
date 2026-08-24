#!/usr/bin/env python
"""Donor-held-out HVG/PCA classical cell-state baselines.

This file is deliberately executable as a script.  The scientific runtime is
Python 3.10, while the dependency-free control plane is Python 3.11.  No
benchmark metric is calculated here; labels are used only to fit the training
fold, and predictions contain no observed label.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any, Iterable, Mapping, Sequence


class ScientificAdapterError(RuntimeError):
    """Raised when the baseline's scientific or provenance contract fails."""


MODEL_IDS = (
    "hvg_pca_elastic_net",
    "hvg_pca_knn",
    "hvg_pca_linear_svm",
    "hvg_pca_logistic",
    "hvg_pca_nearest_centroid",
)
TASK_ID = "cell_state_mapping"
DATASET_ID = "resource_atlas_current"
RUNTIME_ID = "cpu_baseline_smoke"
ENVIRONMENT_SCHEMA = "masld-bench-python-environment-lock-v1"
RECEIPT_SCHEMA = "masld-bench-adapter-receipt-v1"
PREDICTION_SCHEMA = "masld-bench-prediction-bundle-v1"
_HEX = frozenset("0123456789abcdef")
_PARAMETER_FIELDS = frozenset(
    {
        "hvg_mean_bins",
        "join_namespace",
        "calibration_c",
        "calibration_max_iter",
        "elastic_net_l1_ratio",
        "knn_neighbors",
        "normalization_target_sum",
        "n_pca_components",
        "n_top_hvg",
        "outer_folds",
        "pca_svd_solver",
        "regularization_c",
        "svm_calibration_folds",
        "classifier_max_iter",
        "training_weight_policy",
    }
)


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _canonical_hash(value: Any) -> str:
    return sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _artifact_record(path: Path, *, relative_to: Path) -> dict[str, Any]:
    return {
        "path": path.relative_to(relative_to).as_posix(),
        "sha256": _sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def _write_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        handle.write(_canonical_json(value))
        handle.write("\n")


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
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


def _read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ScientificAdapterError(f"TSV has no header: {path}")
        fields = tuple(reader.fieldnames)
        rows = [dict(row) for row in reader]
    return fields, rows


def _require_sha256(value: Any, label: str) -> str:
    text = str(value)
    if len(text) != 64 or any(character not in _HEX for character in text):
        raise ScientificAdapterError(f"{label} must be a lowercase SHA-256")
    return text


def fold_index(unit_id: str, *, seed: int, outer_folds: int) -> int:
    """Assign one biological unit to one deterministic outer fold."""

    if not unit_id or outer_folds < 2 or seed < 0:
        raise ScientificAdapterError("invalid deterministic fold inputs")
    digest = sha256(f"{seed}\0{unit_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % outer_folds


def join_hash(namespace: str, kind: str, identifier: str) -> str:
    if not namespace or kind not in {"row", "unit"} or not identifier:
        raise ScientificAdapterError("invalid join-hash input")
    return sha256(f"{namespace}\0{kind}\0{identifier}".encode("utf-8")).hexdigest()


def donor_class_weights(donors: Sequence[str], labels: Sequence[str]) -> list[float]:
    """Give classes equal mass and donors within each class equal mass."""

    if not donors or len(donors) != len(labels):
        raise ScientificAdapterError("training donors and labels must align")
    counts: dict[tuple[str, str], int] = {}
    donors_by_class: dict[str, set[str]] = {}
    for donor, label in zip(donors, labels):
        counts[(donor, label)] = counts.get((donor, label), 0) + 1
        donors_by_class.setdefault(label, set()).add(donor)
    n_classes = len(donors_by_class)
    raw = [
        1.0
        / (
            n_classes
            * len(donors_by_class[label])
            * counts[(donor, label)]
        )
        for donor, label in zip(donors, labels)
    ]
    scale = len(raw) / sum(raw)
    return [weight * scale for weight in raw]


def _normalized_distribution_records() -> list[str]:
    records = {
        f"{str(distribution.metadata.get('Name')).strip().lower().replace('_', '-')}=={distribution.version}"
        for distribution in importlib.metadata.distributions()
        if distribution.metadata.get("Name")
    }
    return sorted(records)


def _validate_file_artifact(artifact: Mapping[str, Any], label: str) -> Path:
    path = Path(str(artifact.get("path", "")))
    if not path.is_absolute() or path.is_symlink() or not path.is_file():
        raise ScientificAdapterError(f"{label} must be an absolute regular file")
    if path.stat().st_size != artifact.get("size_bytes"):
        raise ScientificAdapterError(f"{label} size changed")
    if _sha256_file(path) != _require_sha256(artifact.get("sha256"), label):
        raise ScientificAdapterError(f"{label} SHA-256 changed")
    return path


def _input_by_role(request: Mapping[str, Any], role: str) -> Mapping[str, Any]:
    run_spec = request.get("run_spec")
    if not isinstance(run_spec, Mapping):
        raise ScientificAdapterError("request has no RunSpec")
    inputs = run_spec.get("inputs")
    if not isinstance(inputs, list):
        raise ScientificAdapterError("RunSpec inputs must be an array")
    matches = [item for item in inputs if isinstance(item, Mapping) and item.get("role") == role]
    if len(matches) != 1:
        raise ScientificAdapterError(f"request must contain exactly one {role} input")
    return matches[0]


def _validate_environment(request: Mapping[str, Any]) -> tuple[Mapping[str, Any], str]:
    artifact = _input_by_role(request, f"environment:{RUNTIME_ID}")
    lock_path = _validate_file_artifact(artifact, "environment lock")
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificAdapterError("environment lock is invalid JSON") from error
    if not isinstance(lock, Mapping) or lock.get("schema_version") != ENVIRONMENT_SCHEMA:
        raise ScientificAdapterError("environment lock schema differs")
    if lock.get("runtime_id") != RUNTIME_ID or lock.get("package_installation_performed") is not False:
        raise ScientificAdapterError("environment lock identity differs")
    if os.environ.get("PYTHONNOUSERSITE") != "1":
        raise ScientificAdapterError("PYTHONNOUSERSITE=1 is required")
    if os.path.realpath(sys.prefix) != lock.get("environment_prefix"):
        raise ScientificAdapterError("active Python prefix differs from the lock")
    if os.path.realpath(sys.executable) != lock.get("python_executable"):
        raise ScientificAdapterError("active Python executable differs from the lock")
    if platform.python_version() != lock.get("python_version"):
        raise ScientificAdapterError("active Python version differs from the lock")
    critical = lock.get("critical_versions")
    if not isinstance(critical, Mapping):
        raise ScientificAdapterError("environment lock has no critical versions")
    observed = {name: importlib.metadata.version(name) for name in critical}
    if observed != dict(critical):
        raise ScientificAdapterError("critical scientific package versions changed")
    distribution_text = "\n".join(_normalized_distribution_records()) + "\n"
    if sha256(distribution_text.encode("utf-8")).hexdigest() != lock.get(
        "python_distributions_sha256"
    ):
        raise ScientificAdapterError("Python distribution inventory changed")
    process = subprocess.run(
        [sys.executable, "-m", "pip", "freeze", "--all"],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env={**os.environ, "PYTHONNOUSERSITE": "1"},
    )
    pip_text = process.stdout.rstrip() + "\n"
    if process.returncode != 0 or sha256(pip_text.encode("utf-8")).hexdigest() != lock.get(
        "pip_freeze_sha256"
    ):
        raise ScientificAdapterError("pip freeze inventory changed")
    return lock, str(artifact["sha256"])


def _validate_request(
    request: Mapping[str, Any], action: str
) -> tuple[Mapping[str, Any], tuple[str, ...], str]:
    if request.get("schema_version") != "masld-bench-adapter-request-v1":
        raise ScientificAdapterError("unsupported adapter request schema")
    if request.get("action") != action:
        raise ScientificAdapterError("request action differs")
    _require_sha256(request.get("run_id"), "run_id")
    run_spec = request.get("run_spec")
    if not isinstance(run_spec, Mapping):
        raise ScientificAdapterError("request has no RunSpec")
    model_id = str(run_spec.get("model_id", ""))
    if model_id not in MODEL_IDS:
        raise ScientificAdapterError("RunSpec names an unsupported classical baseline")
    exact = {
        "task_id": TASK_ID,
        "split_id": "donor_outer",
        "stage": "smoke",
        "adaptation_regime": "native_lane",
        "runtime_id": RUNTIME_ID,
    }
    for field, expected in exact.items():
        if run_spec.get(field) != expected:
            raise ScientificAdapterError(f"RunSpec {field} differs from {expected}")
    if run_spec.get("dataset_ids") != [DATASET_ID]:
        raise ScientificAdapterError("baseline smoke requires only the Resource Atlas view")
    if request.get("fit_dataset_ids") != [DATASET_ID] or request.get("action_dataset_ids") != [DATASET_ID]:
        raise ScientificAdapterError("baseline smoke action scope differs")
    adapter_fields = {
        "input_matrix": "raw_unselected_counts",
        "row_id_source": "obs_names",
        "unit_id_source": "obs.donor_id",
        "label_source": "obs.broad_label",
        "dataset_view_id": "resource_atlas_geneformer_smoke_1000_v1",
    }
    for field, expected in adapter_fields.items():
        if request.get(field) != expected:
            raise ScientificAdapterError(f"adapter request {field} differs")
    parameters = run_spec.get("hyperparameters")
    if not isinstance(parameters, Mapping) or set(parameters) != _PARAMETER_FIELDS:
        raise ScientificAdapterError("baseline hyperparameter inventory differs")
    if parameters.get("pca_svd_solver") != "full" or parameters.get(
        "training_weight_policy"
    ) != "donor_class_balanced_rescaled_to_n_cells":
        raise ScientificAdapterError("baseline algorithm policy differs")
    metadata = run_spec.get("metadata")
    if not isinstance(metadata, Mapping):
        raise ScientificAdapterError("RunSpec metadata is missing")
    split = metadata.get("split_contract")
    if not isinstance(split, Mapping) or split.get("seed") != 20260821:
        raise ScientificAdapterError("frozen donor split seed differs")
    evaluator = metadata.get("evaluator_parameters")
    if not isinstance(evaluator, Mapping):
        raise ScientificAdapterError("frozen evaluator parameters are missing")
    roster = evaluator.get("class_roster")
    if not isinstance(roster, list) or roster != sorted(set(roster)) or len(roster) < 2:
        raise ScientificAdapterError("frozen class roster is invalid")
    numeric_positive = (
        "calibration_c",
        "calibration_max_iter",
        "knn_neighbors",
        "normalization_target_sum",
        "n_pca_components",
        "n_top_hvg",
        "outer_folds",
        "regularization_c",
        "svm_calibration_folds",
        "classifier_max_iter",
    )
    if any(
        isinstance(parameters.get(field), bool)
        or not isinstance(parameters.get(field), (int, float))
        or float(parameters[field]) <= 0
        for field in numeric_positive
    ):
        raise ScientificAdapterError("baseline numeric hyperparameters must be positive")
    l1_ratio = parameters.get("elastic_net_l1_ratio")
    if (
        isinstance(l1_ratio, bool)
        or not isinstance(l1_ratio, (int, float))
        or not 0.0 < float(l1_ratio) < 1.0
    ):
        raise ScientificAdapterError("elastic_net_l1_ratio must be strictly between zero and one")
    if int(parameters["svm_calibration_folds"]) < 2:
        raise ScientificAdapterError("linear-SVM calibration requires at least two donor folds")
    return parameters, tuple(str(item) for item in roster), model_id


def _dataset_path(request: Mapping[str, Any]) -> Path:
    role = "dataset_view_data:resource_atlas_geneformer_smoke_1000_v1"
    return _validate_file_artifact(_input_by_role(request, role), "Atlas smoke H5AD")


def _load_adata(path: Path) -> Any:
    import anndata
    import numpy as np
    from scipy import sparse

    adata = anndata.read_h5ad(path)
    required_obs = {"donor_id", "broad_label", "n_counts"}
    if not required_obs.issubset(adata.obs.columns):
        raise ScientificAdapterError("Atlas H5AD lacks required observation fields")
    if "ensembl_id" not in adata.var.columns:
        raise ScientificAdapterError("Atlas H5AD lacks the Ensembl feature axis")
    if adata.n_obs != 1000 or adata.n_vars != 37533:
        raise ScientificAdapterError("Atlas smoke dimensions changed")
    if len(set(map(str, adata.obs_names))) != adata.n_obs:
        raise ScientificAdapterError("Atlas global cell IDs are not unique")
    gene_ids = list(map(str, adata.var["ensembl_id"]))
    if len(set(gene_ids)) != adata.n_vars or any(not item.startswith("ENSG") for item in gene_ids):
        raise ScientificAdapterError("Atlas Ensembl axis is invalid")
    matrix = sparse.csr_matrix(adata.X, dtype=np.float64)
    if matrix.data.size and (
        not np.all(np.isfinite(matrix.data))
        or np.any(matrix.data < 0)
        or np.any(matrix.data != np.floor(matrix.data))
    ):
        raise ScientificAdapterError("Atlas matrix is not finite nonnegative raw counts")
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(totals <= 0) or np.any(totals != np.asarray(adata.obs["n_counts"], dtype=float)):
        raise ScientificAdapterError("Atlas raw-count totals differ from obs.n_counts")
    adata.X = matrix
    return adata


def _split_rows(adata: Any, parameters: Mapping[str, Any], run_spec: Mapping[str, Any]) -> list[dict[str, Any]]:
    outer_folds = int(parameters["outer_folds"])
    fold = int(run_spec["fold"])
    if outer_folds != 5 or not 0 <= fold < outer_folds:
        raise ScientificAdapterError("outer-fold contract differs")
    namespace = str(parameters["join_namespace"])
    rows: list[dict[str, Any]] = []
    for row_id, donor in zip(map(str, adata.obs_names), map(str, adata.obs["donor_id"])):
        assigned = fold_index(donor, seed=20260821, outer_folds=outer_folds)
        rows.append(
            {
                "row_hash": join_hash(namespace, "row", row_id),
                "unit_hash": join_hash(namespace, "unit", donor),
                "fold": assigned,
                "held_out": assigned == fold,
            }
        )
    rows.sort(key=lambda item: item["row_hash"])
    if any(len({item["fold"] for item in rows if item["unit_hash"] == unit}) != 1 for unit in {item["unit_hash"] for item in rows}):
        raise ScientificAdapterError("one donor appears in multiple folds")
    return rows


def _prior_output(request_path: Path, request: Mapping[str, Any], action: str) -> Path:
    prior = request.get("prior_action_outputs")
    if not isinstance(prior, list):
        raise ScientificAdapterError("prior action inventory is invalid")
    matches = [item for item in prior if isinstance(item, Mapping) and item.get("action") == action]
    if len(matches) != 1:
        raise ScientificAdapterError(f"exactly one prior {action} action is required")
    relative = Path(str(matches[0].get("output_path", "")))
    if relative.is_absolute() or ".." in relative.parts:
        raise ScientificAdapterError("prior output path escapes its attempt")
    root = request_path.parent.parent
    path = root / relative
    if path.is_symlink() or not path.is_dir():
        raise ScientificAdapterError(f"prior {action} output is missing")
    return path


def _receipt(
    *,
    action: str,
    request: Mapping[str, Any],
    output: Path,
    artifacts: Sequence[Path],
    environment_sha256: str,
    extra_metadata: Mapping[str, Any],
) -> None:
    run_spec = request["run_spec"]
    metadata = {
        "adapter": "hvg_pca_classical_v2",
        "environment_artifact_sha256": environment_sha256,
        "runtime_id": RUNTIME_ID,
        "fit_dataset_ids": list(request["fit_dataset_ids"]),
        **dict(extra_metadata),
    }
    receipt = {
        "schema_version": RECEIPT_SCHEMA,
        "action": action,
        "run_id": request["run_id"],
        "status": "complete",
        "artifacts": [_artifact_record(path, relative_to=output) for path in artifacts],
        "metadata": metadata,
    }
    if run_spec.get("task_id") != TASK_ID:
        raise ScientificAdapterError("receipt task binding changed")
    _write_json(output / "adapter_receipt.json", receipt)


def prepare(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    parameters, roster, model_id = _validate_request(request, "prepare")
    _, environment_sha256 = _validate_environment(request)
    adata = _load_adata(_dataset_path(request))
    observed = set(map(str, adata.obs["broad_label"]))
    if observed != set(roster):
        raise ScientificAdapterError("Atlas label set differs from the frozen roster")
    rows = _split_rows(adata, parameters, request["run_spec"])
    split_path = output / "split_rows.tsv"
    _write_tsv(
        split_path,
        ("row_hash", "unit_hash", "fold", "held_out"),
        (
            {
                **row,
                "held_out": "true" if row["held_out"] else "false",
            }
            for row in rows
        ),
    )
    held_rows = [row for row in rows if row["held_out"]]
    summary_path = output / "split_summary.json"
    _write_json(
        summary_path,
        {
            "schema_version": "masld-bench-donor-split-v1",
            "split_id": "donor_outer",
            "split_seed": 20260821,
            "outer_folds": int(parameters["outer_folds"]),
            "held_out_fold": int(request["run_spec"]["fold"]),
            "row_count": len(rows),
            "unit_count": len({row["unit_hash"] for row in rows}),
            "held_out_row_count": len(held_rows),
            "held_out_unit_count": len({row["unit_hash"] for row in held_rows}),
            "row_inventory_sha256": _canonical_hash(rows),
            "labels_exported": False,
        },
    )
    _receipt(
        action="prepare",
        request=request,
        output=output,
        artifacts=(split_path, summary_path),
        environment_sha256=environment_sha256,
        extra_metadata={
            "held_out_fold": int(request["run_spec"]["fold"]),
            "model_id": model_id,
        },
    )


def _log_normalize(matrix: Any, target_sum: float) -> Any:
    import numpy as np
    from scipy import sparse

    value = sparse.csr_matrix(matrix, dtype=np.float64, copy=True)
    totals = np.asarray(value.sum(axis=1)).ravel()
    if not math.isfinite(target_sum) or target_sum <= 0 or np.any(totals <= 0):
        raise ScientificAdapterError("invalid normalization inputs")
    value = sparse.diags(target_sum / totals).dot(value).tocsr()
    value.data = np.log1p(value.data)
    return value


def _select_hvgs(matrix: Any, gene_ids: Sequence[str], *, n_top: int, n_bins: int) -> tuple[list[int], list[dict[str, Any]]]:
    import numpy as np

    if n_top < 1 or n_bins < 2 or n_top > len(gene_ids):
        raise ScientificAdapterError("invalid HVG parameters")
    means = np.asarray(matrix.mean(axis=0)).ravel()
    mean_squares = np.asarray(matrix.power(2).mean(axis=0)).ravel()
    variances = np.maximum(mean_squares - means * means, 0.0)
    dispersions = np.log1p(variances / np.maximum(means, 1e-12))
    mean_order = sorted(range(len(gene_ids)), key=lambda index: (means[index], gene_ids[index]))
    bin_ids = np.empty(len(gene_ids), dtype=np.int64)
    for rank, index in enumerate(mean_order):
        bin_ids[index] = min(n_bins - 1, rank * n_bins // len(gene_ids))
    scores = np.zeros(len(gene_ids), dtype=np.float64)
    for bin_id in range(n_bins):
        members = np.flatnonzero(bin_ids == bin_id)
        values = dispersions[members]
        standard_deviation = float(values.std())
        if standard_deviation > 0:
            scores[members] = (values - values.mean()) / standard_deviation
    selected = sorted(
        range(len(gene_ids)),
        key=lambda index: (-float(scores[index]), gene_ids[index], index),
    )[:n_top]
    records = [
        {
            "feature_index": index,
            "ensembl_id": gene_ids[index],
            "training_mean": format(float(means[index]), ".17g"),
            "training_variance": format(float(variances[index]), ".17g"),
            "normalized_dispersion": format(float(scores[index]), ".17g"),
        }
        for index in selected
    ]
    return selected, records


def _save_array(path: Path, value: Any) -> None:
    import numpy as np

    with path.open("xb") as handle:
        np.save(handle, np.asarray(value), allow_pickle=False)


def _softmax(logits: Any) -> Any:
    import numpy as np

    values = np.asarray(logits, dtype=np.float64)
    if values.ndim != 2 or values.shape[1] < 2 or not np.all(np.isfinite(values)):
        raise ScientificAdapterError("classifier logits are invalid")
    values = values - values.max(axis=1, keepdims=True)
    probabilities = np.exp(values)
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    if not np.all(np.isfinite(probabilities)) or np.any(probabilities < 0):
        raise ScientificAdapterError("classifier produced invalid probabilities")
    return probabilities


def fit(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import numpy as np
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedGroupKFold
    from sklearn.svm import LinearSVC

    parameters, roster, model_id = _validate_request(request, "fit")
    lock, environment_sha256 = _validate_environment(request)
    adata = _load_adata(_dataset_path(request))
    rows = _split_rows(adata, parameters, request["run_spec"])
    prepare_root = _prior_output(request_path, request, "prepare")
    split_fields, prepared_rows = _read_tsv(prepare_root / "split_rows.tsv")
    if split_fields != ("row_hash", "unit_hash", "fold", "held_out"):
        raise ScientificAdapterError("prepared split schema changed")
    expected_prepared = [
        {
            "row_hash": row["row_hash"],
            "unit_hash": row["unit_hash"],
            "fold": str(row["fold"]),
            "held_out": "true" if row["held_out"] else "false",
        }
        for row in rows
    ]
    if prepared_rows != expected_prepared:
        raise ScientificAdapterError("prepared donor split differs from rederivation")
    row_lookup = {row["row_hash"]: row for row in rows}
    namespace = str(parameters["join_namespace"])
    row_hashes = [join_hash(namespace, "row", str(item)) for item in adata.obs_names]
    train_mask = np.asarray([not row_lookup[item]["held_out"] for item in row_hashes], dtype=bool)
    labels = np.asarray(list(map(str, adata.obs["broad_label"])), dtype=str)
    donors = np.asarray(list(map(str, adata.obs["donor_id"])), dtype=str)
    if set(labels[train_mask]) != set(roster):
        raise ScientificAdapterError("training fold does not cover every frozen class")
    normalized = _log_normalize(
        adata.X[train_mask], float(parameters["normalization_target_sum"])
    )
    gene_ids = list(map(str, adata.var["ensembl_id"]))
    selected, feature_records = _select_hvgs(
        normalized,
        gene_ids,
        n_top=int(parameters["n_top_hvg"]),
        n_bins=int(parameters["hvg_mean_bins"]),
    )
    dense = normalized[:, selected].toarray()
    n_components = min(
        int(parameters["n_pca_components"]),
        dense.shape[0] - 1,
        dense.shape[1],
    )
    if n_components < 2:
        raise ScientificAdapterError("training fold cannot support PCA")
    pca = PCA(n_components=n_components, svd_solver="full", whiten=False)
    transformed = pca.fit_transform(dense)
    weights = np.asarray(
        donor_class_weights(list(donors[train_mask]), list(labels[train_mask])),
        dtype=np.float64,
    )
    arrays = {
        "pca_components.npy": np.asarray(pca.components_, dtype=np.float64),
        "pca_mean.npy": np.asarray(pca.mean_, dtype=np.float64),
        "pca_explained_variance.npy": np.asarray(pca.explained_variance_, dtype=np.float64),
    }
    algorithm_parameters: dict[str, Any]
    seed = int(request["run_spec"]["seed"])
    y_train = labels[train_mask]
    donor_train = donors[train_mask]
    if model_id in {"hvg_pca_logistic", "hvg_pca_elastic_net"}:
        elastic = model_id == "hvg_pca_elastic_net"
        classifier = LogisticRegression(
            C=float(parameters["regularization_c"]),
            penalty="elasticnet" if elastic else "l2",
            solver="saga" if elastic else "lbfgs",
            l1_ratio=(float(parameters["elastic_net_l1_ratio"]) if elastic else None),
            max_iter=int(parameters["classifier_max_iter"]),
            fit_intercept=True,
            random_state=seed,
        )
        classifier.fit(transformed, y_train, sample_weight=weights)
        if tuple(map(str, classifier.classes_)) != roster:
            raise ScientificAdapterError("fitted classifier class order differs from roster")
        if np.any(classifier.n_iter_ >= int(parameters["classifier_max_iter"])):
            raise ScientificAdapterError("logistic classifier reached the iteration limit")
        arrays.update(
            {
                "classifier_coef.npy": np.asarray(classifier.coef_, dtype=np.float64),
                "classifier_intercept.npy": np.asarray(
                    classifier.intercept_, dtype=np.float64
                ),
            }
        )
        algorithm_parameters = {
            "classifier": "multinomial_elastic_net_logistic" if elastic else "multinomial_l2_logistic",
            "regularization_c": float(parameters["regularization_c"]),
            "elastic_net_l1_ratio": (
                float(parameters["elastic_net_l1_ratio"]) if elastic else None
            ),
            "solver": "saga" if elastic else "lbfgs",
            "max_iter": int(parameters["classifier_max_iter"]),
        }
    elif model_id == "hvg_pca_nearest_centroid":
        centroids = np.vstack(
            [
                np.average(transformed[y_train == class_id], axis=0, weights=weights[y_train == class_id])
                for class_id in roster
            ]
        )
        class_indices = np.asarray([roster.index(item) for item in y_train], dtype=np.int64)
        squared = np.sum((transformed - centroids[class_indices]) ** 2, axis=1)
        temperature = max(float(np.average(squared, weights=weights)), 1e-12)
        arrays["class_centroids.npy"] = np.asarray(centroids, dtype=np.float64)
        algorithm_parameters = {
            "classifier": "weighted_nearest_centroid",
            "distance": "squared_euclidean_in_training_pca_space",
            "probability_transform": "softmax_negative_squared_distance",
            "temperature": temperature,
            "temperature_fit": "weighted_training_within_class_mean_squared_distance",
        }
    elif model_id == "hvg_pca_knn":
        train_row_hashes = np.asarray(
            [row_hashes[index] for index in np.flatnonzero(train_mask)], dtype=str
        )
        order = np.argsort(train_row_hashes, kind="stable")
        arrays.update(
            {
                "training_embeddings.npy": np.asarray(transformed[order], dtype=np.float64),
                "training_class_indices.npy": np.asarray(
                    [roster.index(item) for item in y_train[order]], dtype=np.int64
                ),
                "training_weights.npy": np.asarray(weights[order], dtype=np.float64),
                "training_row_hashes.npy": np.asarray(train_row_hashes[order], dtype="<U64"),
            }
        )
        neighbors = min(int(parameters["knn_neighbors"]), transformed.shape[0])
        algorithm_parameters = {
            "classifier": "donor_class_weighted_knn",
            "distance": "euclidean_in_training_pca_space",
            "neighbors": neighbors,
            "tie_break": "canonical_training_row_hash_then_class_id",
        }
    elif model_id == "hvg_pca_linear_svm":
        calibration_folds = int(parameters["svm_calibration_folds"])
        splitter = StratifiedGroupKFold(
            n_splits=calibration_folds,
            shuffle=True,
            random_state=seed,
        )
        oof_scores = np.full((len(y_train), len(roster)), np.nan, dtype=np.float64)
        for inner_train, inner_validation in splitter.split(
            transformed, y_train, groups=donor_train
        ):
            if set(y_train[inner_train]) != set(roster):
                raise ScientificAdapterError(
                    "linear-SVM inner donor fold does not cover the class roster"
                )
            inner = LinearSVC(
                C=float(parameters["regularization_c"]),
                dual="auto",
                fit_intercept=True,
                max_iter=int(parameters["classifier_max_iter"]),
                random_state=seed,
            )
            inner.fit(
                transformed[inner_train],
                y_train[inner_train],
                sample_weight=weights[inner_train],
            )
            if tuple(map(str, inner.classes_)) != roster:
                raise ScientificAdapterError("inner linear-SVM class order differs")
            if int(inner.n_iter_) >= int(parameters["classifier_max_iter"]):
                raise ScientificAdapterError("inner linear-SVM reached the iteration limit")
            oof_scores[inner_validation] = inner.decision_function(
                transformed[inner_validation]
            )
        if not np.all(np.isfinite(oof_scores)):
            raise ScientificAdapterError("linear-SVM donor cross-fit calibration is incomplete")
        calibrator = LogisticRegression(
            C=float(parameters["calibration_c"]),
            penalty="l2",
            solver="lbfgs",
            max_iter=int(parameters["calibration_max_iter"]),
            fit_intercept=True,
            random_state=seed,
        )
        calibrator.fit(oof_scores, y_train, sample_weight=weights)
        if tuple(map(str, calibrator.classes_)) != roster:
            raise ScientificAdapterError("linear-SVM calibrator class order differs")
        if np.any(calibrator.n_iter_ >= int(parameters["calibration_max_iter"])):
            raise ScientificAdapterError("linear-SVM calibrator reached the iteration limit")
        classifier = LinearSVC(
            C=float(parameters["regularization_c"]),
            dual="auto",
            fit_intercept=True,
            max_iter=int(parameters["classifier_max_iter"]),
            random_state=seed,
        )
        classifier.fit(transformed, y_train, sample_weight=weights)
        if tuple(map(str, classifier.classes_)) != roster:
            raise ScientificAdapterError("fitted linear-SVM class order differs")
        if int(classifier.n_iter_) >= int(parameters["classifier_max_iter"]):
            raise ScientificAdapterError("linear-SVM reached the iteration limit")
        arrays.update(
            {
                "svm_coef.npy": np.asarray(classifier.coef_, dtype=np.float64),
                "svm_intercept.npy": np.asarray(classifier.intercept_, dtype=np.float64),
                "calibration_coef.npy": np.asarray(calibrator.coef_, dtype=np.float64),
                "calibration_intercept.npy": np.asarray(
                    calibrator.intercept_, dtype=np.float64
                ),
            }
        )
        algorithm_parameters = {
            "classifier": "linear_svm_with_cross_fitted_multinomial_calibration",
            "regularization_c": float(parameters["regularization_c"]),
            "calibration_c": float(parameters["calibration_c"]),
            "calibration_folds": calibration_folds,
            "calibration_group": "donor",
            "max_iter": int(parameters["classifier_max_iter"]),
            "calibration_max_iter": int(parameters["calibration_max_iter"]),
        }
    else:  # pragma: no cover - request validation closes this branch
        raise ScientificAdapterError("unsupported classical classifier")

    array_paths: list[Path] = []
    for name, value in arrays.items():
        path = output / name
        _save_array(path, value)
        array_paths.append(path)
    features_path = output / "selected_features.tsv"
    _write_tsv(
        features_path,
        (
            "feature_index",
            "ensembl_id",
            "training_mean",
            "training_variance",
            "normalized_dispersion",
        ),
        feature_records,
    )
    model_path = output / "fitted_model.json"
    model = {
        "schema_version": "masld-bench-hvg-pca-classical-model-v2",
        "run_id": request["run_id"],
        "model_id": model_id,
        "task_id": TASK_ID,
        "algorithm": algorithm_parameters,
        "class_roster": list(roster),
        "selected_feature_indices": selected,
        "selected_feature_ids_sha256": _canonical_hash(
            [gene_ids[index] for index in selected]
        ),
        "normalization_target_sum": float(parameters["normalization_target_sum"]),
        "n_pca_components": n_components,
        "pca_svd_solver": "full",
        "training_weight_policy": parameters["training_weight_policy"],
        "training_row_count": int(train_mask.sum()),
        "training_unit_count": len(set(donors[train_mask])),
        "held_out_fold": int(request["run_spec"]["fold"]),
        "parameter_sha256": _canonical_hash(dict(parameters)),
        "environment_lock_sha256": environment_sha256,
        "critical_versions": dict(lock["critical_versions"]),
        "arrays": {
            path.name: _artifact_record(path, relative_to=output) for path in array_paths
        },
        "features": _artifact_record(features_path, relative_to=output),
        "observed_labels_exported": False,
    }
    _write_json(model_path, model)
    _receipt(
        action="fit",
        request=request,
        output=output,
        artifacts=(*array_paths, features_path, model_path),
        environment_sha256=environment_sha256,
        extra_metadata={
            "held_out_fold": int(request["run_spec"]["fold"]),
            "training_row_count": int(train_mask.sum()),
            "training_unit_count": len(set(donors[train_mask])),
            "model_id": model_id,
        },
    )


def _load_array(root: Path, record: Mapping[str, Any]) -> Any:
    import numpy as np

    path = root / str(record.get("path", ""))
    if path.parent != root or path.is_symlink() or not path.is_file():
        raise ScientificAdapterError("fitted array path is invalid")
    if path.stat().st_size != record.get("size_bytes") or _sha256_file(path) != record.get("sha256"):
        raise ScientificAdapterError("fitted array artifact changed")
    return np.load(path, allow_pickle=False)


def predict(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import numpy as np

    parameters, roster, model_id = _validate_request(request, "predict")
    _, environment_sha256 = _validate_environment(request)
    adata = _load_adata(_dataset_path(request))
    rows = _split_rows(adata, parameters, request["run_spec"])
    fit_root = _prior_output(request_path, request, "fit")
    model_path = fit_root / "fitted_model.json"
    try:
        model = json.loads(model_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificAdapterError("fitted model manifest is invalid") from error
    if (
        not isinstance(model, Mapping)
        or model.get("run_id") != request["run_id"]
        or model.get("model_id") != model_id
        or model.get("class_roster") != list(roster)
        or model.get("parameter_sha256") != _canonical_hash(dict(parameters))
        or model.get("environment_lock_sha256") != environment_sha256
    ):
        raise ScientificAdapterError("fitted model binding differs")
    arrays = model.get("arrays")
    if not isinstance(arrays, Mapping):
        raise ScientificAdapterError("fitted model array inventory is missing")
    pca_components = _load_array(fit_root, arrays["pca_components.npy"])
    pca_mean = _load_array(fit_root, arrays["pca_mean.npy"])
    selected = [int(item) for item in model["selected_feature_indices"]]
    namespace = str(parameters["join_namespace"])
    row_hashes = [join_hash(namespace, "row", str(item)) for item in adata.obs_names]
    row_lookup = {row["row_hash"]: row for row in rows}
    held_indices = [index for index, item in enumerate(row_hashes) if row_lookup[item]["held_out"]]
    if not held_indices:
        raise ScientificAdapterError("held-out fold contains no rows")
    normalized = _log_normalize(
        adata.X[held_indices], float(model["normalization_target_sum"])
    )
    dense = normalized[:, selected].toarray()
    transformed = (dense - pca_mean) @ pca_components.T
    algorithm = model.get("algorithm")
    if not isinstance(algorithm, Mapping):
        raise ScientificAdapterError("fitted model algorithm contract is missing")
    if model_id in {"hvg_pca_logistic", "hvg_pca_elastic_net"}:
        coefficients = _load_array(fit_root, arrays["classifier_coef.npy"])
        intercept = _load_array(fit_root, arrays["classifier_intercept.npy"])
        probabilities = _softmax(transformed @ coefficients.T + intercept)
    elif model_id == "hvg_pca_nearest_centroid":
        centroids = _load_array(fit_root, arrays["class_centroids.npy"])
        squared_distances = np.sum(
            (transformed[:, np.newaxis, :] - centroids[np.newaxis, :, :]) ** 2,
            axis=2,
        )
        temperature = float(algorithm.get("temperature", 0.0))
        if not math.isfinite(temperature) or temperature <= 0:
            raise ScientificAdapterError("nearest-centroid temperature is invalid")
        probabilities = _softmax(-squared_distances / temperature)
    elif model_id == "hvg_pca_knn":
        training = _load_array(fit_root, arrays["training_embeddings.npy"])
        training_classes = _load_array(
            fit_root, arrays["training_class_indices.npy"]
        ).astype(np.int64, copy=False)
        training_weights = _load_array(fit_root, arrays["training_weights.npy"])
        training_row_hashes = _load_array(fit_root, arrays["training_row_hashes.npy"])
        if (
            training.ndim != 2
            or training.shape[1] != transformed.shape[1]
            or training.shape[0]
            != training_classes.shape[0]
            != training_weights.shape[0]
            != training_row_hashes.shape[0]
            or np.any(training_classes < 0)
            or np.any(training_classes >= len(roster))
            or np.any(training_weights <= 0)
            or list(map(str, training_row_hashes))
            != sorted(map(str, training_row_hashes))
        ):
            raise ScientificAdapterError("fitted kNN training inventory is invalid")
        neighbors = int(algorithm.get("neighbors", 0))
        if not 0 < neighbors <= training.shape[0]:
            raise ScientificAdapterError("fitted kNN neighbor count is invalid")
        probabilities = np.zeros((transformed.shape[0], len(roster)), dtype=np.float64)
        for row_index, embedding in enumerate(transformed):
            squared = np.sum((training - embedding) ** 2, axis=1)
            nearest = np.argsort(squared, kind="stable")[:neighbors]
            votes = np.bincount(
                training_classes[nearest],
                weights=training_weights[nearest],
                minlength=len(roster),
            )
            probabilities[row_index] = votes / votes.sum()
        if not np.all(np.isfinite(probabilities)) or np.any(probabilities < 0):
            raise ScientificAdapterError("kNN produced invalid probabilities")
    elif model_id == "hvg_pca_linear_svm":
        svm_coefficients = _load_array(fit_root, arrays["svm_coef.npy"])
        svm_intercept = _load_array(fit_root, arrays["svm_intercept.npy"])
        calibration_coefficients = _load_array(
            fit_root, arrays["calibration_coef.npy"]
        )
        calibration_intercept = _load_array(
            fit_root, arrays["calibration_intercept.npy"]
        )
        scores = transformed @ svm_coefficients.T + svm_intercept
        probabilities = _softmax(
            scores @ calibration_coefficients.T + calibration_intercept
        )
    else:  # pragma: no cover - request validation closes this branch
        raise ScientificAdapterError("unsupported classical classifier")
    records: list[dict[str, Any]] = []
    for position, source_index in enumerate(held_indices):
        row_hash = row_hashes[source_index]
        unit_hash = row_lookup[row_hash]["unit_hash"]
        values = probabilities[position]
        winner = min(
            range(len(roster)),
            key=lambda index: (-float(values[index]), roster[index]),
        )
        records.append(
            {
                "row_hash": row_hash,
                "unit_hash": unit_hash,
                "predicted_class": roster[winner],
                "probabilities": [float(value) for value in values],
            }
        )
    records.sort(key=lambda item: item["row_hash"])

    prediction_fields = ("row_hash", "unit_hash", "predicted_class")
    prediction_path = output / "predictions.tsv"
    _write_tsv(
        prediction_path,
        prediction_fields,
        ({field: row[field] for field in prediction_fields} for row in records),
    )
    probability_fields = (
        "row_hash",
        "unit_hash",
        *(f"probability::{class_id}" for class_id in roster),
    )
    probability_path = output / "class_probabilities.tsv"
    _write_tsv(
        probability_path,
        probability_fields,
        (
            {
                "row_hash": row["row_hash"],
                "unit_hash": row["unit_hash"],
                **{
                    f"probability::{class_id}": format(row["probabilities"][index], ".17g")
                    for index, class_id in enumerate(roster)
                },
            }
            for row in records
        ),
    )
    row_ids_path = output / "row_ids.tsv"
    _write_tsv(
        row_ids_path,
        ("row_hash", "unit_hash"),
        ({"row_hash": row["row_hash"], "unit_hash": row["unit_hash"]} for row in records),
    )
    run_spec = request["run_spec"]
    standardized = {
        **_artifact_record(prediction_path, relative_to=output),
        "media_type": "text/tab-separated-values",
        "role": f"standardized_prediction_table:{TASK_ID}",
    }
    row_artifact = {
        **_artifact_record(row_ids_path, relative_to=output),
        "media_type": "text/tab-separated-values",
        "role": f"prediction_row_ids:{TASK_ID}",
    }
    probability_artifact = {
        **_artifact_record(probability_path, relative_to=output),
        "media_type": "text/tab-separated-values",
        "role": f"class_probabilities:{TASK_ID}",
    }
    source_join = _canonical_hash(
        {
            "task_id": TASK_ID,
            "dataset_ids": list(run_spec["dataset_ids"]),
            "split_id": run_spec["split_id"],
            "row_id_field": "row_hash",
            "unit_id_field": "unit_hash",
            "unit_id_namespace": namespace,
            "biological_unit": "donor",
        }
    )
    bundle = {
        "schema_version": PREDICTION_SCHEMA,
        "bundle_id": f"{model_id}-{str(request['run_id'])[:16]}",
        "run_id": request["run_id"],
        "task_id": TASK_ID,
        "model_id": model_id,
        "dataset_ids": list(run_spec["dataset_ids"]),
        "split_id": run_spec["split_id"],
        "artifacts": [standardized, row_artifact, probability_artifact],
        "standardized_table": standardized,
        "row_ids": row_artifact,
        "n_predictions": len(records),
        "row_id_field": "row_hash",
        "unit_id_field": "unit_hash",
        "unit_id_namespace": namespace,
        "biological_unit": "donor",
        "table_schema_sha256": _canonical_hash(
            {"format": "tsv", "fields": list(prediction_fields)}
        ),
        "source_join_key_sha256": source_join,
        "format_version": "tsv-v1",
        "missing_state": "observed",
        "metadata": {
            "held_out_fold": int(run_spec["fold"]),
            "observed_labels_exported": False,
            "smoke_only": True,
            "algorithm": algorithm.get("classifier"),
        },
    }
    bundle_path = output / "prediction_bundle.json"
    _write_json(bundle_path, bundle)
    _receipt(
        action="predict",
        request=request,
        output=output,
        artifacts=(prediction_path, probability_path, row_ids_path, bundle_path),
        environment_sha256=environment_sha256,
        extra_metadata={
            "held_out_fold": int(run_spec["fold"]),
            "prediction_row_count": len(records),
            "prediction_unit_count": len({row["unit_hash"] for row in records}),
            "observed_labels_exported": False,
            "model_id": model_id,
        },
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--action", required=True, choices=("prepare", "fit", "predict"))
    parser.add_argument("--request", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args(argv)
    try:
        request = json.loads(arguments.request.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificAdapterError("adapter request is invalid JSON") from error
    if not isinstance(request, Mapping):
        raise ScientificAdapterError("adapter request must be an object")
    if arguments.output.exists():
        raise ScientificAdapterError(f"adapter output already exists: {arguments.output}")
    arguments.output.mkdir(parents=True, exist_ok=False)
    actions = {"prepare": prepare, "fit": fit, "predict": predict}
    actions[arguments.action](arguments.request.resolve(strict=True), request, arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

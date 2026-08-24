#!/usr/bin/env python
"""Donor-held-out classical RNA-to-ATAC profile baselines.

The adapter never calculates benchmark metrics. ``prepare`` is the sole action
allowed to receive the paired development HDF5. It emits training RNA/ATAC and
held-donor RNA, but never held-donor ATAC. ``fit`` and ``predict`` operate only
on those frozen prior-action artifacts.
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


class RNAATACAdapterError(RuntimeError):
    """Raised when an RNA-to-ATAC baseline contract is violated."""


MODEL_IDS = ("assay_native_pseudobulk", "mean_track", "nearest_context")
TASK_ID = "rna_conditioned_atac"
DATASET_ID = "gse296875"
VIEW_ID = "gse296875_rna_atac_smoke_1000_v1"
RUNTIME_ID = "cpu_baseline_smoke"
DATA_ROLE = f"dataset_view_data:{VIEW_ID}"
LINEAGES = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)
ENVIRONMENT_SCHEMA = "masld-bench-python-environment-lock-v1"
RECEIPT_SCHEMA = "masld-bench-adapter-receipt-v1"
PREDICTION_SCHEMA = "masld-bench-prediction-bundle-v1"
_HEX = frozenset("0123456789abcdef")
_PARAMETER_FIELDS = frozenset(
    {
        "join_namespace",
        "normalization_target_sum",
        "n_context_hvg",
        "n_smoke_peaks",
        "outer_folds",
        "profile_pseudocount",
        "split_seed",
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


def _require_sha256(value: Any, label: str) -> str:
    text = str(value)
    if len(text) != 64 or any(character not in _HEX for character in text):
        raise RNAATACAdapterError(f"{label} must be a lowercase SHA-256")
    return text


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


def _write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
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
        for row in rows:
            writer.writerow(dict(row))


def _read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RNAATACAdapterError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def fold_index(unit_id: str, *, seed: int, outer_folds: int) -> int:
    if not unit_id or seed < 0 or outer_folds < 2:
        raise RNAATACAdapterError("invalid donor fold inputs")
    digest = sha256(f"{seed}\0{unit_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % outer_folds


def join_hash(namespace: str, kind: str, identifier: str) -> str:
    if not namespace or kind not in {"row", "unit", "block"} or not identifier:
        raise RNAATACAdapterError("invalid join-hash input")
    return sha256(f"{namespace}\0{kind}\0{identifier}".encode("utf-8")).hexdigest()


def deterministic_peak_indices(peak_ids: Sequence[str], n_peaks: int) -> tuple[int, ...]:
    """Select an ATAC-blind fixed smoke universe using identifiers only."""

    if len(set(peak_ids)) != len(peak_ids) or not 1 <= n_peaks <= len(peak_ids):
        raise RNAATACAdapterError("invalid peak inventory or smoke-peak budget")
    ranked = sorted(
        range(len(peak_ids)),
        key=lambda index: (
            sha256(f"masld-rna-atac-smoke-v1\0{peak_ids[index]}".encode()).digest(),
            peak_ids[index],
            index,
        ),
    )
    return tuple(sorted(ranked[:n_peaks]))


def _normalized_distribution_records() -> list[str]:
    return sorted(
        {
            f"{str(distribution.metadata.get('Name')).strip().lower().replace('_', '-')}=={distribution.version}"
            for distribution in importlib.metadata.distributions()
            if distribution.metadata.get("Name")
        }
    )


def _input_by_role(request: Mapping[str, Any], role: str) -> Mapping[str, Any]:
    run_spec = request.get("run_spec")
    if not isinstance(run_spec, Mapping) or not isinstance(run_spec.get("inputs"), list):
        raise RNAATACAdapterError("RunSpec inputs are missing")
    matches = [
        item
        for item in run_spec["inputs"]
        if isinstance(item, Mapping) and item.get("role") == role
    ]
    if len(matches) != 1:
        raise RNAATACAdapterError(f"request must contain exactly one {role} input")
    return matches[0]


def _validate_file_artifact(artifact: Mapping[str, Any], label: str) -> Path:
    path = Path(str(artifact.get("path", "")))
    if not path.is_absolute() or path.is_symlink() or not path.is_file():
        raise RNAATACAdapterError(f"{label} must be an absolute regular file")
    if path.stat().st_size != artifact.get("size_bytes"):
        raise RNAATACAdapterError(f"{label} size changed")
    if _sha256_file(path) != _require_sha256(artifact.get("sha256"), label):
        raise RNAATACAdapterError(f"{label} SHA-256 changed")
    return path


def _validate_environment(request: Mapping[str, Any]) -> tuple[Mapping[str, Any], str]:
    artifact = _input_by_role(request, f"environment:{RUNTIME_ID}")
    lock_path = _validate_file_artifact(artifact, "environment lock")
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RNAATACAdapterError("environment lock is invalid JSON") from error
    if not isinstance(lock, Mapping) or lock.get("schema_version") != ENVIRONMENT_SCHEMA:
        raise RNAATACAdapterError("environment lock schema differs")
    if lock.get("runtime_id") != RUNTIME_ID or lock.get("package_installation_performed") is not False:
        raise RNAATACAdapterError("environment lock identity differs")
    if os.environ.get("PYTHONNOUSERSITE") != "1":
        raise RNAATACAdapterError("PYTHONNOUSERSITE=1 is required")
    if os.path.realpath(sys.prefix) != lock.get("environment_prefix"):
        raise RNAATACAdapterError("active Python prefix differs from the lock")
    if os.path.realpath(sys.executable) != lock.get("python_executable"):
        raise RNAATACAdapterError("active Python executable differs from the lock")
    if platform.python_version() != lock.get("python_version"):
        raise RNAATACAdapterError("active Python version differs from the lock")
    critical = lock.get("critical_versions")
    if not isinstance(critical, Mapping) or {
        name: importlib.metadata.version(name) for name in critical
    } != dict(critical):
        raise RNAATACAdapterError("critical scientific package versions changed")
    distribution_text = "\n".join(_normalized_distribution_records()) + "\n"
    if sha256(distribution_text.encode()).hexdigest() != lock.get(
        "python_distributions_sha256"
    ):
        raise RNAATACAdapterError("Python distribution inventory changed")
    process = subprocess.run(
        [sys.executable, "-m", "pip", "freeze", "--all"],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env={**os.environ, "PYTHONNOUSERSITE": "1"},
    )
    if (
        process.returncode != 0
        or sha256((process.stdout.rstrip() + "\n").encode()).hexdigest()
        != lock.get("pip_freeze_sha256")
    ):
        raise RNAATACAdapterError("pip freeze inventory changed")
    return lock, str(artifact["sha256"])


def _validate_request(
    request: Mapping[str, Any], action: str
) -> tuple[Mapping[str, Any], str]:
    if request.get("schema_version") != "masld-bench-adapter-request-v1":
        raise RNAATACAdapterError("unsupported adapter request schema")
    if request.get("action") != action:
        raise RNAATACAdapterError("request action differs")
    _require_sha256(request.get("run_id"), "run_id")
    run_spec = request.get("run_spec")
    if not isinstance(run_spec, Mapping):
        raise RNAATACAdapterError("request has no RunSpec")
    model_id = str(run_spec.get("model_id", ""))
    if model_id not in MODEL_IDS:
        raise RNAATACAdapterError("unsupported RNA-ATAC classical baseline")
    exact = {
        "task_id": TASK_ID,
        "split_id": "donor_outer",
        "stage": "smoke",
        "adaptation_regime": "native_lane",
        "runtime_id": RUNTIME_ID,
    }
    for field, expected in exact.items():
        if run_spec.get(field) != expected:
            raise RNAATACAdapterError(f"RunSpec {field} differs from {expected}")
    if run_spec.get("dataset_ids") != [DATASET_ID]:
        raise RNAATACAdapterError("RNA-ATAC smoke requires only GSE296875")
    if request.get("fit_dataset_ids") != [DATASET_ID] or request.get(
        "action_dataset_ids"
    ) != [DATASET_ID]:
        raise RNAATACAdapterError("RNA-ATAC action scope differs")
    if request.get("dataset_view_id") != VIEW_ID:
        raise RNAATACAdapterError("RNA-ATAC dataset view differs")
    expected_withheld = {"fit": [DATA_ROLE], "predict": [DATA_ROLE]}
    if request.get("withheld_input_roles_by_action") != expected_withheld:
        raise RNAATACAdapterError("held-ATAC input firewall differs")
    roles = {
        str(item.get("role"))
        for item in run_spec.get("inputs", [])
        if isinstance(item, Mapping)
    }
    if action == "prepare":
        if DATA_ROLE not in roles:
            raise RNAATACAdapterError("prepare requires the paired development HDF5")
    elif DATA_ROLE in roles:
        raise RNAATACAdapterError(f"{action} request exposes the paired ATAC HDF5")
    parameters = run_spec.get("hyperparameters")
    if not isinstance(parameters, Mapping) or set(parameters) != _PARAMETER_FIELDS:
        raise RNAATACAdapterError("RNA-ATAC hyperparameter inventory differs")
    integer_fields = ("n_context_hvg", "n_smoke_peaks", "outer_folds", "split_seed")
    if any(
        isinstance(parameters.get(field), bool)
        or not isinstance(parameters.get(field), int)
        or int(parameters[field]) < 1
        for field in integer_fields
    ):
        raise RNAATACAdapterError("integer hyperparameters must be positive")
    for field in ("normalization_target_sum", "profile_pseudocount"):
        value = parameters.get(field)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)) or float(value) <= 0:
            raise RNAATACAdapterError(f"{field} must be finite and positive")
    if int(parameters["outer_folds"]) != 5 or int(parameters["split_seed"]) != 20260821:
        raise RNAATACAdapterError("donor split contract differs")
    metadata = run_spec.get("metadata")
    if not isinstance(metadata, Mapping):
        raise RNAATACAdapterError("RunSpec metadata is missing")
    evaluator = metadata.get("evaluator_parameters")
    if not isinstance(evaluator, Mapping) or evaluator.get("strata") != list(LINEAGES):
        raise RNAATACAdapterError("RNA-ATAC lineage roster differs")
    fold = run_spec.get("fold")
    if isinstance(fold, bool) or not isinstance(fold, int) or not 0 <= fold < 5:
        raise RNAATACAdapterError("held donor fold differs")
    return parameters, model_id


def _decode(values: Any) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def _read_csr(group: Any) -> Any:
    import numpy as np
    from scipy import sparse

    shape = tuple(int(item) for item in group["shape"][:])
    matrix = sparse.csr_matrix(
        (
            np.asarray(group["data"][:]),
            np.asarray(group["indices"][:], dtype=np.int64),
            np.asarray(group["indptr"][:], dtype=np.int64),
        ),
        shape=shape,
    )
    if matrix.data.size and (
        not np.all(np.isfinite(matrix.data))
        or np.any(matrix.data < 0)
        or np.any(matrix.data != np.floor(matrix.data))
    ):
        raise RNAATACAdapterError("multimodal matrix is not raw nonnegative counts")
    return matrix


def _aggregate(matrix: Any, groups: Sequence[tuple[str, str]]) -> tuple[Any, list[dict[str, Any]]]:
    import numpy as np
    from scipy import sparse

    by_group: dict[tuple[str, str], list[int]] = {}
    for index, group in enumerate(groups):
        by_group.setdefault(group, []).append(index)
    ordered = sorted(by_group)
    rows = []
    metadata = []
    for donor, lineage in ordered:
        indices = by_group[(donor, lineage)]
        rows.append(sparse.csr_matrix(matrix[indices].sum(axis=0)))
        metadata.append(
            {"donor_id": donor, "lineage": lineage, "n_nuclei": len(indices)}
        )
    result = sparse.vstack(rows, format="csr")
    totals = np.asarray(result.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise RNAATACAdapterError("donor-lineage pseudobulk contains an empty profile")
    return result, metadata


def _save_sparse(path: Path, matrix: Any) -> None:
    from scipy import sparse

    with path.open("xb") as handle:
        sparse.save_npz(handle, sparse.csr_matrix(matrix), compressed=True)


def _load_sparse(path: Path) -> Any:
    from scipy import sparse

    if path.is_symlink() or not path.is_file():
        raise RNAATACAdapterError(f"sparse artifact is missing: {path}")
    return sparse.load_npz(path).tocsr()


def _prior_output(request_path: Path, request: Mapping[str, Any], action: str) -> Path:
    prior = request.get("prior_action_outputs")
    if not isinstance(prior, list):
        raise RNAATACAdapterError("prior action inventory is invalid")
    matches = [
        item
        for item in prior
        if isinstance(item, Mapping) and item.get("action") == action
    ]
    if len(matches) != 1:
        raise RNAATACAdapterError(f"exactly one prior {action} output is required")
    relative = Path(str(matches[0].get("output_path", "")))
    if relative.is_absolute() or ".." in relative.parts:
        raise RNAATACAdapterError("prior output path escapes its attempt")
    path = request_path.parent.parent / relative
    if path.is_symlink() or not path.is_dir():
        raise RNAATACAdapterError(f"prior {action} output is missing")
    return path


def _receipt(
    *,
    action: str,
    request: Mapping[str, Any],
    output: Path,
    artifacts: Sequence[Path],
    environment_sha256: str,
    metadata: Mapping[str, Any],
) -> None:
    _write_json(
        output / "adapter_receipt.json",
        {
            "schema_version": RECEIPT_SCHEMA,
            "action": action,
            "run_id": request["run_id"],
            "status": "complete",
            "artifacts": [
                _artifact_record(path, relative_to=output) for path in artifacts
            ],
            "metadata": {
                "adapter": "rna_atac_classical_v1",
                "runtime_id": RUNTIME_ID,
                "environment_artifact_sha256": environment_sha256,
                "fit_dataset_ids": list(request["fit_dataset_ids"]),
                **dict(metadata),
            },
        },
    )


def prepare(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import h5py
    import numpy as np

    parameters, model_id = _validate_request(request, "prepare")
    _, environment_sha256 = _validate_environment(request)
    source = _validate_file_artifact(_input_by_role(request, DATA_ROLE), "paired HDF5")
    with h5py.File(source, "r") as handle:
        if handle.attrs.get("schema_version") != "masld-bench-multimodal-h5-v1":
            raise RNAATACAdapterError("paired HDF5 schema differs")
        donors = _decode(handle["obs/donor_id"][:])
        lineages = _decode(handle["obs/broad_label"][:])
        if set(lineages) != set(LINEAGES) or len(donors) != 1000:
            raise RNAATACAdapterError("paired HDF5 rows or lineages differ")
        rna = _read_csr(handle["rna/counts_csr"])
        atac = _read_csr(handle["atac/counts_csr"])
        gene_ids = _decode(handle["rna/ensembl_id"][:])
        peak_ids = _decode(handle["atac/peak_id"][:])
        chromosomes = _decode(handle["atac/chromosome"][:])
        starts = np.asarray(handle["atac/bed_start_0based"][:], dtype=np.int64)
        ends = np.asarray(handle["atac/bed_end_half_open"][:], dtype=np.int64)
    if rna.shape != (1000, len(gene_ids)) or atac.shape != (1000, len(peak_ids)):
        raise RNAATACAdapterError("paired HDF5 matrix axes differ")
    fold = int(request["run_spec"]["fold"])
    assigned = [
        fold_index(
            donor,
            seed=int(parameters["split_seed"]),
            outer_folds=int(parameters["outer_folds"]),
        )
        for donor in donors
    ]
    train_indices = np.asarray([value != fold for value in assigned], dtype=bool)
    query_indices = ~train_indices
    if not train_indices.any() or not query_indices.any():
        raise RNAATACAdapterError("donor split produced an empty partition")
    for donor in set(donors):
        donor_folds = {assigned[index] for index, value in enumerate(donors) if value == donor}
        if len(donor_folds) != 1:
            raise RNAATACAdapterError("one donor appears in multiple folds")
    selected = deterministic_peak_indices(peak_ids, int(parameters["n_smoke_peaks"]))
    train_groups = [
        (donors[index], lineages[index]) for index in np.flatnonzero(train_indices)
    ]
    query_groups = [
        (donors[index], lineages[index]) for index in np.flatnonzero(query_indices)
    ]
    train_rna, train_rows = _aggregate(rna[train_indices], train_groups)
    train_atac, train_atac_rows = _aggregate(
        atac[train_indices][:, selected], train_groups
    )
    query_rna, query_rows = _aggregate(rna[query_indices], query_groups)
    if train_rows != train_atac_rows:
        raise RNAATACAdapterError("training RNA and ATAC pseudobulk rows differ")
    paths = {
        "training_rna": output / "training_rna.npz",
        "training_atac": output / "training_atac.npz",
        "query_rna": output / "query_rna.npz",
    }
    _save_sparse(paths["training_rna"], train_rna)
    _save_sparse(paths["training_atac"], train_atac)
    _save_sparse(paths["query_rna"], query_rna)
    row_fields = ("donor_id", "lineage", "n_nuclei")
    train_rows_path = output / "training_rows.tsv"
    query_rows_path = output / "query_rows.tsv"
    _write_tsv(train_rows_path, row_fields, train_rows)
    _write_tsv(query_rows_path, row_fields, query_rows)
    peaks_path = output / "selected_peaks.tsv"
    _write_tsv(
        peaks_path,
        ("selected_index", "source_index", "peak_id", "chromosome", "bed_start", "bed_end"),
        (
            {
                "selected_index": selected_index,
                "source_index": source_index,
                "peak_id": peak_ids[source_index],
                "chromosome": chromosomes[source_index],
                "bed_start": int(starts[source_index]),
                "bed_end": int(ends[source_index]),
            }
            for selected_index, source_index in enumerate(selected)
        ),
    )
    axes_path = output / "axes.json"
    axes = {
        "schema_version": "masld-bench-rna-atac-prepared-v1",
        "run_id": request["run_id"],
        "model_id": model_id,
        "held_out_fold": fold,
        "split_seed": int(parameters["split_seed"]),
        "gene_ids_sha256": _canonical_hash(gene_ids),
        "selected_peak_ids_sha256": _canonical_hash([peak_ids[index] for index in selected]),
        "training_rows_sha256": _canonical_hash(train_rows),
        "query_rows_sha256": _canonical_hash(query_rows),
        "training_donors": len({row["donor_id"] for row in train_rows}),
        "query_donors": len({row["donor_id"] for row in query_rows}),
        "selected_peak_count": len(selected),
        "held_atac_exported": False,
        "query_modalities": {"rna": "observed", "atac": "withheld_sealed"},
    }
    _write_json(axes_path, axes)
    _receipt(
        action="prepare",
        request=request,
        output=output,
        artifacts=(
            *paths.values(),
            train_rows_path,
            query_rows_path,
            peaks_path,
            axes_path,
        ),
        environment_sha256=environment_sha256,
        metadata={
            "model_id": model_id,
            "held_out_fold": fold,
            "held_atac_exported": False,
            "training_donor_count": axes["training_donors"],
            "query_donor_count": axes["query_donors"],
        },
    )


def _row_normalize(matrix: Any) -> Any:
    import numpy as np
    from scipy import sparse

    value = sparse.csr_matrix(matrix, dtype=np.float64, copy=True)
    totals = np.asarray(value.sum(axis=1)).ravel()
    if np.any(~np.isfinite(totals)) or np.any(totals <= 0):
        raise RNAATACAdapterError("profile normalization received an empty row")
    return sparse.diags(1.0 / totals).dot(value).tocsr()


def _log_normalize(matrix: Any, target_sum: float) -> Any:
    from scipy import sparse

    value = _row_normalize(matrix) * float(target_sum)
    value = sparse.csr_matrix(value)
    value.data = __import__("numpy").log1p(value.data)
    return value


def _select_context_hvgs(matrix: Any, n_top: int) -> tuple[int, ...]:
    import numpy as np

    means = np.asarray(matrix.mean(axis=0)).ravel()
    variances = np.maximum(
        np.asarray(matrix.power(2).mean(axis=0)).ravel() - means * means,
        0.0,
    )
    scores = variances / np.maximum(means, 1e-12)
    selected = sorted(
        range(matrix.shape[1]), key=lambda index: (-float(scores[index]), index)
    )[:n_top]
    return tuple(sorted(selected))


def fit(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import numpy as np
    from scipy import sparse

    parameters, model_id = _validate_request(request, "fit")
    _, environment_sha256 = _validate_environment(request)
    prepared = _prior_output(request_path, request, "prepare")
    train_rna = _load_sparse(prepared / "training_rna.npz")
    train_atac = _load_sparse(prepared / "training_atac.npz")
    fields, rows = _read_tsv(prepared / "training_rows.tsv")
    if fields != ("donor_id", "lineage", "n_nuclei") or len(rows) != train_rna.shape[0] or train_rna.shape[0] != train_atac.shape[0]:
        raise RNAATACAdapterError("prepared training rows or matrices differ")
    profile_rows = _row_normalize(train_atac)
    profile_labels: list[dict[str, str]]
    if model_id == "mean_track":
        profiles = sparse.csr_matrix(profile_rows.mean(axis=0))
        profile_labels = [{"profile_id": "global", "lineage": "all", "donor_id": "training_mean"}]
    elif model_id == "assay_native_pseudobulk":
        values = []
        profile_labels = []
        for lineage in LINEAGES:
            indices = [index for index, row in enumerate(rows) if row["lineage"] == lineage]
            if not indices:
                raise RNAATACAdapterError(f"training fold lacks lineage {lineage}")
            values.append(sparse.csr_matrix(profile_rows[indices].mean(axis=0)))
            profile_labels.append({"profile_id": lineage, "lineage": lineage, "donor_id": "training_mean"})
        profiles = sparse.vstack(values, format="csr")
    else:
        profiles = profile_rows
        profile_labels = [
            {
                "profile_id": f"{row['donor_id']}::{row['lineage']}",
                "lineage": row["lineage"],
                "donor_id": row["donor_id"],
            }
            for row in rows
        ]
    profiles_path = output / "profiles.npz"
    _save_sparse(profiles_path, profiles)
    labels_path = output / "profile_index.tsv"
    _write_tsv(labels_path, ("profile_id", "lineage", "donor_id"), profile_labels)
    artifacts: list[Path] = [profiles_path, labels_path]
    context_contract: dict[str, Any] = {"used": False}
    if model_id == "nearest_context":
        normalized = _log_normalize(
            train_rna, float(parameters["normalization_target_sum"])
        )
        n_hvg = min(int(parameters["n_context_hvg"]), normalized.shape[1])
        selected = _select_context_hvgs(normalized, n_hvg)
        context = normalized[:, selected].toarray()
        mean = context.mean(axis=0)
        scale = context.std(axis=0)
        scale[scale == 0] = 1.0
        context = (context - mean) / scale
        arrays = {
            "context_reference.npy": context,
            "context_mean.npy": mean,
            "context_scale.npy": scale,
            "selected_gene_indices.npy": np.asarray(selected, dtype=np.int64),
        }
        array_records = {}
        for name, value in arrays.items():
            path = output / name
            with path.open("xb") as handle:
                np.save(handle, value, allow_pickle=False)
            artifacts.append(path)
            array_records[name] = _artifact_record(path, relative_to=output)
        context_contract = {
            "used": True,
            "distance": "squared_euclidean_training_zscored_log1p_cpm_hvg",
            "tie_break": "profile_id_lexicographic",
            "selected_gene_count": len(selected),
            "selected_gene_indices_sha256": _canonical_hash(list(selected)),
            "arrays": array_records,
        }
    state_path = output / "fitted_model.json"
    _write_json(
        state_path,
        {
            "schema_version": "masld-bench-rna-atac-classical-model-v1",
            "run_id": request["run_id"],
            "model_id": model_id,
            "held_out_fold": int(request["run_spec"]["fold"]),
            "profile_count": profiles.shape[0],
            "peak_count": profiles.shape[1],
            "training_stratum_count": len(rows),
            "training_donor_count": len({row["donor_id"] for row in rows}),
            "profile_scale": "donor_lineage_equal_weight_depth_free_composition",
            "profile_pseudocount": float(parameters["profile_pseudocount"]),
            "context": context_contract,
            "parameter_sha256": _canonical_hash(dict(parameters)),
            "environment_lock_sha256": environment_sha256,
            "profiles": _artifact_record(profiles_path, relative_to=output),
            "profile_index": _artifact_record(labels_path, relative_to=output),
            "query_rna_used_for_fit": False,
            "held_atac_used_for_fit": False,
        },
    )
    artifacts.append(state_path)
    _receipt(
        action="fit",
        request=request,
        output=output,
        artifacts=artifacts,
        environment_sha256=environment_sha256,
        metadata={
            "model_id": model_id,
            "held_out_fold": int(request["run_spec"]["fold"]),
            "query_rna_used_for_fit": False,
            "held_atac_used_for_fit": False,
        },
    )


def _load_array(root: Path, record: Mapping[str, Any]) -> Any:
    import numpy as np

    path = root / str(record.get("path", ""))
    if path.parent != root or path.is_symlink() or not path.is_file():
        raise RNAATACAdapterError("fitted array path is invalid")
    if path.stat().st_size != record.get("size_bytes") or _sha256_file(path) != record.get("sha256"):
        raise RNAATACAdapterError("fitted array artifact changed")
    return np.load(path, allow_pickle=False)


def predict(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import numpy as np

    parameters, model_id = _validate_request(request, "predict")
    _, environment_sha256 = _validate_environment(request)
    prepared = _prior_output(request_path, request, "prepare")
    fitted = _prior_output(request_path, request, "fit")
    try:
        model = json.loads((fitted / "fitted_model.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RNAATACAdapterError("fitted model manifest is invalid") from error
    if (
        not isinstance(model, Mapping)
        or model.get("run_id") != request["run_id"]
        or model.get("model_id") != model_id
        or model.get("parameter_sha256") != _canonical_hash(dict(parameters))
        or model.get("environment_lock_sha256") != environment_sha256
        or model.get("held_atac_used_for_fit") is not False
    ):
        raise RNAATACAdapterError("fitted model binding differs")
    profiles = _load_sparse(fitted / "profiles.npz")
    fields, profile_rows = _read_tsv(fitted / "profile_index.tsv")
    query_fields, query_rows = _read_tsv(prepared / "query_rows.tsv")
    peak_fields, peak_rows = _read_tsv(prepared / "selected_peaks.tsv")
    if fields != ("profile_id", "lineage", "donor_id") or len(profile_rows) != profiles.shape[0]:
        raise RNAATACAdapterError("fitted profile index differs")
    if query_fields != ("donor_id", "lineage", "n_nuclei"):
        raise RNAATACAdapterError("query row schema differs")
    if peak_fields != ("selected_index", "source_index", "peak_id", "chromosome", "bed_start", "bed_end") or len(peak_rows) != profiles.shape[1]:
        raise RNAATACAdapterError("selected peak schema differs")
    selected_profiles: list[int] = []
    if model_id == "mean_track":
        selected_profiles = [0] * len(query_rows)
    elif model_id == "assay_native_pseudobulk":
        by_lineage = {row["lineage"]: index for index, row in enumerate(profile_rows)}
        if set(by_lineage) != set(LINEAGES):
            raise RNAATACAdapterError("lineage profile roster differs")
        selected_profiles = [by_lineage[row["lineage"]] for row in query_rows]
    else:
        query_rna = _load_sparse(prepared / "query_rna.npz")
        context = model.get("context")
        if not isinstance(context, Mapping) or context.get("used") is not True:
            raise RNAATACAdapterError("nearest-context state is missing")
        arrays = context.get("arrays")
        if not isinstance(arrays, Mapping):
            raise RNAATACAdapterError("nearest-context arrays are missing")
        reference = _load_array(fitted, arrays["context_reference.npy"])
        mean = _load_array(fitted, arrays["context_mean.npy"])
        scale = _load_array(fitted, arrays["context_scale.npy"])
        selected = _load_array(fitted, arrays["selected_gene_indices.npy"]).astype(np.int64)
        query_context = _log_normalize(
            query_rna, float(parameters["normalization_target_sum"])
        )[:, selected].toarray()
        query_context = (query_context - mean) / scale
        for query_index, row in enumerate(query_rows):
            eligible = [
                index for index, item in enumerate(profile_rows) if item["lineage"] == row["lineage"]
            ]
            if not eligible:
                raise RNAATACAdapterError("nearest context has no same-lineage reference")
            selected_profiles.append(
                min(
                    eligible,
                    key=lambda index: (
                        float(np.sum((query_context[query_index] - reference[index]) ** 2)),
                        profile_rows[index]["profile_id"],
                    ),
                )
            )
    namespace = str(parameters["join_namespace"])
    pseudocount = float(parameters["profile_pseudocount"])
    prediction_fields = (
        "row_hash",
        "donor_hash",
        "block_hash",
        "stratum",
        "predicted",
    )
    prediction_path = output / "predictions.tsv"
    row_ids_path = output / "row_ids.tsv"
    prediction_rows: list[dict[str, str]] = []
    row_id_rows: list[dict[str, str]] = []
    for query_index, query in enumerate(query_rows):
        profile = profiles.getrow(selected_profiles[query_index]).toarray().ravel()
        donor_hash = join_hash(namespace, "unit", query["donor_id"])
        for peak_index, peak in enumerate(peak_rows):
            row_hash = join_hash(
                namespace,
                "row",
                f"{query['donor_id']}\0{query['lineage']}\0{peak['peak_id']}",
            )
            prediction_rows.append(
                {
                    "row_hash": row_hash,
                    "donor_hash": donor_hash,
                    "block_hash": join_hash(namespace, "block", peak["chromosome"]),
                    "stratum": query["lineage"],
                    "predicted": format(float(profile[peak_index]) + pseudocount, ".17g"),
                }
            )
            row_id_rows.append({"row_hash": row_hash, "donor_hash": donor_hash})
    prediction_rows.sort(key=lambda row: row["row_hash"])
    row_id_rows.sort(key=lambda row: row["row_hash"])
    _write_tsv(prediction_path, prediction_fields, prediction_rows)
    _write_tsv(row_ids_path, ("row_hash", "donor_hash"), row_id_rows)
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
    source_join = _canonical_hash(
        {
            "task_id": TASK_ID,
            "dataset_ids": [DATASET_ID],
            "split_id": "donor_outer",
            "row_id_field": "row_hash",
            "unit_id_field": "donor_hash",
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
        "dataset_ids": [DATASET_ID],
        "split_id": "donor_outer",
        "artifacts": [standardized, row_artifact],
        "standardized_table": standardized,
        "row_ids": row_artifact,
        "n_predictions": len(prediction_rows),
        "row_id_field": "row_hash",
        "unit_id_field": "donor_hash",
        "unit_id_namespace": namespace,
        "biological_unit": "donor",
        "table_schema_sha256": _canonical_hash(
            {"format": "tsv", "fields": list(prediction_fields)}
        ),
        "source_join_key_sha256": source_join,
        "format_version": "tsv-v1",
        "missing_state": "observed",
        "metadata": {
            "held_out_fold": int(request["run_spec"]["fold"]),
            "query_stratum_count": len(query_rows),
            "selected_peak_count": len(peak_rows),
            "prediction_scale": "depth_free_profile_score_plus_uniform_training_frozen_pseudocount",
            "profile_pseudocount": pseudocount,
            "held_atac_input_exposed": False,
            "observed_atac_exported": False,
            "smoke_only": True,
        },
    }
    bundle_path = output / "prediction_bundle.json"
    _write_json(bundle_path, bundle)
    _receipt(
        action="predict",
        request=request,
        output=output,
        artifacts=(prediction_path, row_ids_path, bundle_path),
        environment_sha256=environment_sha256,
        metadata={
            "model_id": model_id,
            "held_out_fold": int(request["run_spec"]["fold"]),
            "prediction_row_count": len(prediction_rows),
            "prediction_donor_count": len({row["donor_hash"] for row in prediction_rows}),
            "held_atac_input_exposed": False,
            "observed_atac_exported": False,
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
        raise RNAATACAdapterError("adapter request is invalid JSON") from error
    if not isinstance(request, Mapping):
        raise RNAATACAdapterError("adapter request must be an object")
    if arguments.output.exists():
        raise RNAATACAdapterError(f"adapter output already exists: {arguments.output}")
    arguments.output.mkdir(parents=True, exist_ok=False)
    if arguments.action == "prepare":
        prepare(arguments.request, request, arguments.output)
    elif arguments.action == "fit":
        fit(arguments.request, request, arguments.output)
    else:
        predict(arguments.request, request, arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

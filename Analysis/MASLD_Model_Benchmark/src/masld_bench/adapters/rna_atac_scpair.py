#!/usr/bin/env python
"""Donor-held-out scPair v0.1 RNA-to-ATAC direct tensor adapter.

The pinned upstream high-level object expects both modalities in one combined
object and derives preprocessing and covariates before splitting. This adapter
reproduces the v0.1 RNA-to-ATAC feature predictor through its audited tensor
equations. Only ``prepare`` may read the paired development HDF5. ``fit`` sees
outer-training same-nucleus RNA/ATAC pairs, and ``predict`` sees held RNA only.

This 1,000-nucleus view uses a source-fixed consensus peak inventory that was
created by the data authors across the cohort. The run is therefore a smoke
fixture and cannot support a best-model or external-transfer claim.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import importlib.util
import json
import math
import os
from pathlib import Path
import random
import sys
from typing import Any, Iterable, Mapping, Sequence


class RNAATACScPairError(RuntimeError):
    """Raised when a scPair smoke requirement is not met."""


def _load_classical() -> Any:
    path = Path(__file__).with_name("rna_atac_classical.py")
    spec = importlib.util.spec_from_file_location(
        "masld_bench_standalone_rna_atac_classical_for_scpair", path
    )
    if spec is None or spec.loader is None:
        raise RNAATACScPairError(f"cannot load classical adapter helpers: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


base = _load_classical()

MODEL_ID = "scpair"
TASK_ID = base.TASK_ID
DATASET_ID = base.DATASET_ID
VIEW_ID = base.VIEW_ID
DATA_ROLE = base.DATA_ROLE
LINEAGES = base.LINEAGES
RUNTIME_ID = "gpu_rna_atac_torch_smoke"
RECEIPT_SCHEMA = base.RECEIPT_SCHEMA
PREDICTION_SCHEMA = base.PREDICTION_SCHEMA
UPSTREAM_REVISION = "52d1ca1e26c8127b2023f3d80cfd1c140a0ea80a"
base.RUNTIME_ID = RUNTIME_ID

_PARAMETER_FIELDS = frozenset(
    {
        "batch_size",
        "dropout_rate",
        "early_stopping_patience",
        "hidden_layers",
        "inference_batch_size",
        "join_namespace",
        "learning_rate",
        "max_epochs",
        "n_hvg",
        "n_smoke_peaks",
        "outer_folds",
        "profile_pseudocount",
        "split_seed",
        "validation_fraction",
        "weight_decay",
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


def _write_array(path: Path, value: Any) -> None:
    import numpy as np

    with path.open("xb") as handle:
        np.save(handle, value, allow_pickle=False)


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
                base._artifact_record(path, relative_to=output) for path in artifacts
            ],
            "metadata": {
                "adapter": "rna_atac_scpair_v1",
                "runtime_id": RUNTIME_ID,
                "environment_artifact_sha256": environment_sha256,
                "fit_dataset_ids": list(request["fit_dataset_ids"]),
                "upstream_revision": UPSTREAM_REVISION,
                **dict(metadata),
            },
        },
    )


def _validate_parameters(parameters: Any) -> Mapping[str, Any]:
    if not isinstance(parameters, Mapping) or set(parameters) != _PARAMETER_FIELDS:
        raise RNAATACScPairError("scPair hyperparameter inventory differs")
    integer_fields = (
        "batch_size",
        "early_stopping_patience",
        "inference_batch_size",
        "max_epochs",
        "n_hvg",
        "n_smoke_peaks",
        "outer_folds",
        "split_seed",
    )
    for field in integer_fields:
        value = parameters.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise RNAATACScPairError(f"{field} must be a positive integer")
    float_fields = (
        "dropout_rate",
        "learning_rate",
        "profile_pseudocount",
        "validation_fraction",
        "weight_decay",
    )
    for field in float_fields:
        value = parameters.get(field)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
        ):
            raise RNAATACScPairError(f"{field} must be finite")
    if parameters["hidden_layers"] != [900, 40]:
        raise RNAATACScPairError("scPair native hidden layers must be [900, 40]")
    if not 0.0 <= float(parameters["dropout_rate"]) < 1.0:
        raise RNAATACScPairError("dropout_rate is outside [0, 1)")
    if not 0.0 < float(parameters["validation_fraction"]) < 0.5:
        raise RNAATACScPairError("validation_fraction is outside (0, 0.5)")
    if (
        float(parameters["learning_rate"]) <= 0
        or float(parameters["profile_pseudocount"]) <= 0
        or float(parameters["weight_decay"]) < 0
    ):
        raise RNAATACScPairError("positive scale hyperparameters differ")
    if (
        int(parameters["outer_folds"]) != 5
        or int(parameters["split_seed"]) != 20260821
        or int(parameters["inference_batch_size"]) != 1
    ):
        raise RNAATACScPairError("donor split or query-invariance contract differs")
    return parameters


def _validate_request(
    request: Mapping[str, Any], action: str
) -> Mapping[str, Any]:
    if request.get("schema_version") != "masld-bench-adapter-request-v1":
        raise RNAATACScPairError("unsupported adapter request schema")
    if request.get("action") != action:
        raise RNAATACScPairError("request action differs")
    base._require_sha256(request.get("run_id"), "run_id")
    run_spec = request.get("run_spec")
    if not isinstance(run_spec, Mapping):
        raise RNAATACScPairError("request has no RunSpec")
    exact = {
        "model_id": MODEL_ID,
        "task_id": TASK_ID,
        "split_id": "donor_outer",
        "stage": "smoke",
        "adaptation_regime": "native_lane",
        "runtime_id": RUNTIME_ID,
    }
    for field, expected in exact.items():
        if run_spec.get(field) != expected:
            raise RNAATACScPairError(f"RunSpec {field} differs from {expected}")
    if run_spec.get("dataset_ids") != [DATASET_ID]:
        raise RNAATACScPairError("scPair smoke requires only GSE296875")
    if request.get("fit_dataset_ids") != [DATASET_ID] or request.get(
        "action_dataset_ids"
    ) != [DATASET_ID]:
        raise RNAATACScPairError("scPair action scope differs")
    if request.get("dataset_view_id") != VIEW_ID:
        raise RNAATACScPairError("scPair dataset view differs")
    expected_withheld = {"fit": [DATA_ROLE], "predict": [DATA_ROLE]}
    if request.get("withheld_input_roles_by_action") != expected_withheld:
        raise RNAATACScPairError("held-ATAC firewall differs")
    roles = {
        str(item.get("role"))
        for item in run_spec.get("inputs", [])
        if isinstance(item, Mapping)
    }
    if action == "prepare":
        if DATA_ROLE not in roles:
            raise RNAATACScPairError("prepare requires the paired development HDF5")
    elif DATA_ROLE in roles:
        raise RNAATACScPairError(f"{action} request exposes the paired ATAC HDF5")
    fold = run_spec.get("fold")
    seed = run_spec.get("seed")
    if (
        isinstance(fold, bool)
        or not isinstance(fold, int)
        or not 0 <= fold < 5
        or isinstance(seed, bool)
        or not isinstance(seed, int)
        or seed < 0
    ):
        raise RNAATACScPairError("fold or seed differs")
    metadata = run_spec.get("metadata")
    evaluator = metadata.get("evaluator_parameters") if isinstance(metadata, Mapping) else None
    if not isinstance(evaluator, Mapping) or evaluator.get("strata") != list(LINEAGES):
        raise RNAATACScPairError("scPair lineage roster differs")
    return _validate_parameters(run_spec.get("hyperparameters"))


def _validate_torch_runtime() -> Mapping[str, Any]:
    import torch

    if os.environ.get("CUBLAS_WORKSPACE_CONFIG") != ":4096:8":
        raise RNAATACScPairError("CUBLAS_WORKSPACE_CONFIG=:4096:8 is required")
    if torch.__version__ != "2.3.1+cu121":
        raise RNAATACScPairError("torch version differs from the frozen runtime")
    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise RNAATACScPairError("exactly one visible CUDA device is required")
    device_name = torch.cuda.get_device_name(0)
    capability = tuple(int(value) for value in torch.cuda.get_device_capability(0))
    if "L40S" not in device_name or capability != (8, 9):
        raise RNAATACScPairError("runtime is not the admitted L40S device")
    torch.use_deterministic_algorithms(True)
    torch.backends.cudnn.enabled = False
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    return {
        "torch_version": torch.__version__,
        "cuda_runtime": torch.version.cuda,
        "device_name": device_name,
        "device_capability": list(capability),
        "deterministic_algorithms": torch.are_deterministic_algorithms_enabled(),
    }


def _set_seed(seed: int) -> None:
    import numpy as np
    import torch

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def _make_model(
    *, input_dim: int, output_dim: int, hidden_layers: Sequence[int], dropout: float
) -> Any:
    import torch
    from torch import nn

    if list(hidden_layers) != [900, 40]:
        raise RNAATACScPairError("model constructor received non-native hidden layers")

    def fcnn(dimensions: Sequence[int]) -> nn.Sequential:
        layers: list[nn.Module] = []
        for source, target in zip(dimensions[:-1], dimensions[1:], strict=True):
            layers.extend(
                (
                    nn.Linear(int(source), int(target)),
                    nn.LayerNorm(int(target)),
                    nn.LeakyReLU(),
                    nn.Dropout(float(dropout)),
                )
            )
        return nn.Sequential(*layers)

    class Predictor(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            dimensions = [int(input_dim), *map(int, hidden_layers)]
            self.hidden_rep = fcnn(dimensions)
            self.extra_linear_rep = nn.Linear(hidden_layers[-1], hidden_layers[-1])
            self.hidden_lib = nn.Sequential(
                fcnn(dimensions), nn.Linear(hidden_layers[-1], 1)
            )
            self.output_layer = nn.Linear(hidden_layers[-1], int(output_dim))

        def forward(self, values: Any) -> tuple[Any, Any, Any]:
            latent_rep = self.extra_linear_rep(self.hidden_rep(values))
            latent_lib = self.hidden_lib(values)
            output_px = torch.sigmoid(self.output_layer(latent_rep))
            output_result = output_px * torch.sigmoid(latent_lib)
            return output_result, output_px, latent_rep

    model = Predictor()
    for module in model.modules():
        if isinstance(module, nn.Linear):
            nn.init.xavier_normal_(module.weight)
            nn.init.zeros_(module.bias)
    return model


def inner_validation_indices(
    rows: Sequence[Mapping[str, str]], *, seed: int, validation_fraction: float
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return a deterministic donor-disjoint inner train/validation split."""

    donors = sorted({row["donor_id"] for row in rows})
    if len(donors) < 5:
        raise RNAATACScPairError("inner validation requires at least five donors")
    ranked = sorted(
        donors,
        key=lambda donor: (
            sha256(f"scpair-inner-v1\0{seed}\0{donor}".encode()).digest(),
            donor,
        ),
    )
    n_validation = max(1, int(round(len(donors) * validation_fraction)))
    validation_donors = frozenset(ranked[:n_validation])
    training = tuple(
        index for index, row in enumerate(rows) if row["donor_id"] not in validation_donors
    )
    validation = tuple(
        index for index, row in enumerate(rows) if row["donor_id"] in validation_donors
    )
    if not training or not validation:
        raise RNAATACScPairError("inner donor split is empty")
    if {rows[index]["donor_id"] for index in training} & {
        rows[index]["donor_id"] for index in validation
    }:
        raise RNAATACScPairError("inner donor split overlaps")
    for lineage in LINEAGES:
        if not any(rows[index]["lineage"] == lineage for index in training):
            raise RNAATACScPairError(f"inner training lacks {lineage}")
        if not any(rows[index]["lineage"] == lineage for index in validation):
            raise RNAATACScPairError(f"inner validation lacks {lineage}")
    return training, validation


def _selected_hvgs(matrix: Any, n_hvg: int) -> tuple[int, ...]:
    normalized = base._log_normalize(matrix, 10000.0)
    return base._select_context_hvgs(normalized, min(n_hvg, matrix.shape[1]))


def prepare(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import h5py
    import numpy as np

    parameters = _validate_request(request, "prepare")
    _, environment_sha256 = base._validate_environment(request)
    runtime = _validate_torch_runtime()
    source = base._validate_file_artifact(
        base._input_by_role(request, DATA_ROLE), "paired HDF5"
    )
    with h5py.File(source, "r") as handle:
        if (
            handle.attrs.get("schema_version") != "masld-bench-multimodal-h5-v1"
            or handle.attrs.get("pairing_level") != "same_nucleus"
        ):
            raise RNAATACScPairError("paired HDF5 schema or topology differs")
        cells = base._decode(handle["obs/cell_id"][:])
        donors = base._decode(handle["obs/donor_id"][:])
        wells = base._decode(handle["obs/well_id"][:])
        lineages = base._decode(handle["obs/broad_label"][:])
        rna_status = base._decode(handle["obs/rna_status"][:])
        atac_status = base._decode(handle["obs/atac_status"][:])
        rna_mask = np.asarray(handle["obs/rna_observed_mask"][:])
        atac_mask = np.asarray(handle["obs/atac_observed_mask"][:])
        rna = base._read_csr(handle["rna/counts_csr"])
        atac = base._read_csr(handle["atac/counts_csr"])
        gene_ids = base._decode(handle["rna/ensembl_id"][:])
        peak_ids = base._decode(handle["atac/peak_id"][:])
        chromosomes = base._decode(handle["atac/chromosome"][:])
        starts = np.asarray(handle["atac/bed_start_0based"][:], dtype=np.int64)
        ends = np.asarray(handle["atac/bed_end_half_open"][:], dtype=np.int64)
    if (
        len(cells) != 1000
        or len(set(cells)) != 1000
        or set(lineages) != set(LINEAGES)
        or set(rna_status) != {"observed"}
        or set(atac_status) != {"observed"}
        or not np.all(rna_mask == 1)
        or not np.all(atac_mask == 1)
        or rna.shape[0] != 1000
        or atac.shape[0] != 1000
    ):
        raise RNAATACScPairError("same-nucleus row or missingness contract differs")
    if len(gene_ids) != rna.shape[1] or len(peak_ids) != atac.shape[1]:
        raise RNAATACScPairError("feature identities differ from matrix axes")

    fold = int(request["run_spec"]["fold"])
    assigned = [
        base.fold_index(
            donor,
            seed=int(parameters["split_seed"]),
            outer_folds=int(parameters["outer_folds"]),
        )
        for donor in donors
    ]
    training_mask = np.asarray([value != fold for value in assigned], dtype=bool)
    query_mask = ~training_mask
    for donor in set(donors):
        if len({assigned[index] for index, item in enumerate(donors) if item == donor}) != 1:
            raise RNAATACScPairError("one donor appears in multiple outer folds")
    selected_peaks = base.deterministic_peak_indices(
        peak_ids, int(parameters["n_smoke_peaks"])
    )
    selected_genes = _selected_hvgs(rna[training_mask], int(parameters["n_hvg"]))

    training_rna = rna[training_mask][:, selected_genes].tocsr()
    query_rna = rna[query_mask][:, selected_genes].tocsr()
    training_atac = atac[training_mask][:, selected_peaks].tocsr()
    training_atac.data = np.ones_like(training_atac.data, dtype=np.float32)
    training_atac.eliminate_zeros()
    if training_rna.shape[0] == 0 or query_rna.shape[0] == 0 or training_atac.nnz == 0:
        raise RNAATACScPairError("outer split produced an empty tensor")

    def metadata(mask: Any) -> list[dict[str, str]]:
        return [
            {
                "cell_id": cells[index],
                "donor_id": donors[index],
                "lineage": lineages[index],
                "well_id": wells[index],
            }
            for index in np.flatnonzero(mask)
        ]

    training_rows = metadata(training_mask)
    query_rows = metadata(query_mask)
    paths = {
        "training_rna": output / "training_rna.npz",
        "training_atac": output / "training_atac_binary.npz",
        "query_rna": output / "query_rna.npz",
    }
    base._save_sparse(paths["training_rna"], training_rna)
    base._save_sparse(paths["training_atac"], training_atac)
    base._save_sparse(paths["query_rna"], query_rna)
    row_fields = ("cell_id", "donor_id", "lineage", "well_id")
    training_rows_path = output / "training_rows.tsv"
    query_rows_path = output / "query_rows.tsv"
    _write_tsv(training_rows_path, row_fields, training_rows)
    _write_tsv(query_rows_path, row_fields, query_rows)
    genes_path = output / "selected_genes.tsv"
    _write_tsv(
        genes_path,
        ("selected_index", "source_index", "ensembl_id"),
        (
            {
                "selected_index": selected_index,
                "source_index": source_index,
                "ensembl_id": gene_ids[source_index],
            }
            for selected_index, source_index in enumerate(selected_genes)
        ),
    )
    peaks_path = output / "selected_peaks.tsv"
    _write_tsv(
        peaks_path,
        (
            "selected_index",
            "source_index",
            "peak_id",
            "chromosome",
            "bed_start",
            "bed_end",
        ),
        (
            {
                "selected_index": selected_index,
                "source_index": source_index,
                "peak_id": peak_ids[source_index],
                "chromosome": chromosomes[source_index],
                "bed_start": int(starts[source_index]),
                "bed_end": int(ends[source_index]),
            }
            for selected_index, source_index in enumerate(selected_peaks)
        ),
    )
    axes_path = output / "axes.json"
    axes = {
        "schema_version": "masld-bench-rna-atac-scpair-prepared-v1",
        "run_id": request["run_id"],
        "model_id": MODEL_ID,
        "held_out_fold": fold,
        "split_seed": int(parameters["split_seed"]),
        "training_cell_count": len(training_rows),
        "query_cell_count": len(query_rows),
        "training_donors": len({row["donor_id"] for row in training_rows}),
        "query_donors": len({row["donor_id"] for row in query_rows}),
        "selected_gene_count": len(selected_genes),
        "selected_peak_count": len(selected_peaks),
        "selected_gene_ids_sha256": _canonical_hash(
            [gene_ids[index] for index in selected_genes]
        ),
        "selected_peak_ids_sha256": _canonical_hash(
            [peak_ids[index] for index in selected_peaks]
        ),
        "training_rows_sha256": _canonical_hash(training_rows),
        "query_rows_sha256": _canonical_hash(query_rows),
        "pairing": "same_nucleus_exact_row_identity",
        "training_atac_scale": "binary_peak_accessibility",
        "query_modalities": {"rna": "observed", "atac": "structurally_missing"},
        "held_atac_exported": False,
        "query_atac_placeholder_created": False,
        "hvg_selection": "outer_training_only_log1p_cpm_dispersion",
        "peak_selection": "source_fixed_identifier_only_sha256_rank",
        "source_feature_inventory_includes_held_donor_peak_discovery": True,
        "smoke_only": True,
        "champion_claim_allowed": False,
        "runtime": runtime,
    }
    _write_json(axes_path, axes)
    _receipt(
        action="prepare",
        request=request,
        output=output,
        artifacts=(
            *paths.values(),
            training_rows_path,
            query_rows_path,
            genes_path,
            peaks_path,
            axes_path,
        ),
        environment_sha256=environment_sha256,
        metadata={
            "model_id": MODEL_ID,
            "held_out_fold": fold,
            "training_cell_count": len(training_rows),
            "query_cell_count": len(query_rows),
            "held_atac_exported": False,
            "query_atac_placeholder_created": False,
            "smoke_only": True,
        },
    )


def _state_records(state: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "key": key,
            "shape": list(value.shape),
            "dtype": str(value.dtype),
            "finite": bool(value.is_floating_point() is False or value.isfinite().all()),
        }
        for key, value in sorted(state.items())
    ]


def _optimizer(model: Any, *, learning_rate: float, weight_decay: float) -> Any:
    import torch

    decay = []
    no_decay = []
    for name, parameter in model.named_parameters():
        if not parameter.requires_grad:
            continue
        if parameter.ndim == 1 or name.startswith("output"):
            no_decay.append(parameter)
        else:
            decay.append(parameter)
    return torch.optim.AdamW(
        (
            {"params": decay, "weight_decay": weight_decay},
            {"params": no_decay, "weight_decay": 0.0},
        ),
        lr=learning_rate,
    )


def fit(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import numpy as np
    import torch
    from torch.utils.data import DataLoader, TensorDataset

    parameters = _validate_request(request, "fit")
    _, environment_sha256 = base._validate_environment(request)
    runtime = _validate_torch_runtime()
    prepared = base._prior_output(request_path, request, "prepare")
    training_rna = base._load_sparse(prepared / "training_rna.npz")
    training_atac = base._load_sparse(prepared / "training_atac_binary.npz")
    fields, rows = base._read_tsv(prepared / "training_rows.tsv")
    if (
        fields != ("cell_id", "donor_id", "lineage", "well_id")
        or len(rows) != training_rna.shape[0]
        or training_rna.shape[0] != training_atac.shape[0]
        or training_atac.shape[1] != int(parameters["n_smoke_peaks"])
        or (training_atac.data.size and not np.all(training_atac.data == 1))
    ):
        raise RNAATACScPairError("prepared training tensors or identities differ")
    inner_train, inner_validation = inner_validation_indices(
        rows,
        seed=int(request["run_spec"]["seed"]),
        validation_fraction=float(parameters["validation_fraction"]),
    )
    train_x = torch.from_numpy(
        training_rna[list(inner_train)].toarray().astype(np.float32, copy=False)
    )
    train_y = torch.from_numpy(
        training_atac[list(inner_train)].toarray().astype(np.float32, copy=False)
    )
    validation_x = torch.from_numpy(
        training_rna[list(inner_validation)].toarray().astype(np.float32, copy=False)
    ).cuda()
    validation_y = torch.from_numpy(
        training_atac[list(inner_validation)].toarray().astype(np.float32, copy=False)
    ).cuda()
    if not torch.isfinite(train_x).all() or not torch.isfinite(train_y).all():
        raise RNAATACScPairError("training tensors are non-finite")

    seed = int(request["run_spec"]["seed"])
    _set_seed(seed)
    model = _make_model(
        input_dim=training_rna.shape[1],
        output_dim=training_atac.shape[1],
        hidden_layers=parameters["hidden_layers"],
        dropout=float(parameters["dropout_rate"]),
    ).cuda()
    optimizer = _optimizer(
        model,
        learning_rate=float(parameters["learning_rate"]),
        weight_decay=float(parameters["weight_decay"]),
    )
    generator = torch.Generator(device="cpu")
    generator.manual_seed(seed)
    loader = DataLoader(
        TensorDataset(train_x, train_y),
        batch_size=int(parameters["batch_size"]),
        shuffle=True,
        num_workers=0,
        generator=generator,
    )
    criterion = torch.nn.BCELoss()
    best_loss = math.inf
    best_epoch = 0
    epochs_without_improvement = 0
    best_state: dict[str, Any] | None = None
    history: list[dict[str, Any]] = []
    for epoch in range(1, int(parameters["max_epochs"]) + 1):
        model.train()
        train_loss_sum = 0.0
        train_observations = 0
        for batch_x, batch_y in loader:
            batch_x = batch_x.cuda(non_blocking=False)
            batch_y = batch_y.cuda(non_blocking=False)
            optimizer.zero_grad(set_to_none=True)
            output_result, _, _ = model(batch_x)
            loss = criterion(output_result, batch_y)
            if not torch.isfinite(loss):
                raise RNAATACScPairError("training loss is non-finite")
            loss.backward()
            optimizer.step()
            train_loss_sum += float(loss.detach().cpu()) * batch_x.shape[0]
            train_observations += int(batch_x.shape[0])
        model.eval()
        with torch.no_grad():
            validation_result, _, _ = model(validation_x)
            validation_loss = float(criterion(validation_result, validation_y).cpu())
        history.append(
            {
                "epoch": epoch,
                "training_bce": train_loss_sum / train_observations,
                "validation_bce": validation_loss,
            }
        )
        if validation_loss < best_loss:
            best_loss = validation_loss
            best_epoch = epoch
            epochs_without_improvement = 0
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
        else:
            epochs_without_improvement += 1
        if epochs_without_improvement >= int(parameters["early_stopping_patience"]):
            break
    if best_state is None or best_epoch < 1 or not math.isfinite(best_loss):
        raise RNAATACScPairError("early stopping did not produce a valid checkpoint")

    state_path = output / "state_dict.pt"
    with state_path.open("xb") as handle:
        torch.save(best_state, handle)
    state_records = _state_records(best_state)
    state_records_path = output / "state_dict_manifest.json"
    _write_json(state_records_path, state_records)
    history_path = output / "training_history.tsv"
    _write_tsv(history_path, ("epoch", "training_bce", "validation_bce"), history)
    model_path = output / "fitted_model.json"
    model_manifest = {
        "schema_version": "masld-bench-rna-atac-scpair-model-v1",
        "run_id": request["run_id"],
        "model_id": MODEL_ID,
        "upstream_revision": UPSTREAM_REVISION,
        "held_out_fold": int(request["run_spec"]["fold"]),
        "seed": seed,
        "input_gene_count": training_rna.shape[1],
        "output_peak_count": training_atac.shape[1],
        "hidden_layers": list(parameters["hidden_layers"]),
        "dropout_rate": float(parameters["dropout_rate"]),
        "training_loss": "binary_cross_entropy_on_output_result",
        "prediction_output": "output_px_before_source_RNA_library_scaling",
        "inner_training_donor_count": len(
            {rows[index]["donor_id"] for index in inner_train}
        ),
        "inner_validation_donor_count": len(
            {rows[index]["donor_id"] for index in inner_validation}
        ),
        "inner_training_rows_sha256": _canonical_hash(
            [rows[index]["cell_id"] for index in inner_train]
        ),
        "inner_validation_rows_sha256": _canonical_hash(
            [rows[index]["cell_id"] for index in inner_validation]
        ),
        "best_epoch": best_epoch,
        "best_validation_bce": best_loss,
        "epochs_completed": len(history),
        "state_dict": base._artifact_record(state_path, relative_to=output),
        "state_dict_manifest": base._artifact_record(
            state_records_path, relative_to=output
        ),
        "state_key_shape_dtype_sha256": _canonical_hash(state_records),
        "parameter_sha256": _canonical_hash(dict(parameters)),
        "environment_lock_sha256": environment_sha256,
        "runtime": runtime,
        "query_rna_used_for_fit": False,
        "held_atac_used_for_fit": False,
        "whole_module_pickle_saved": False,
        "smoke_only": True,
        "champion_claim_allowed": False,
    }
    _write_json(model_path, model_manifest)
    _receipt(
        action="fit",
        request=request,
        output=output,
        artifacts=(state_path, state_records_path, history_path, model_path),
        environment_sha256=environment_sha256,
        metadata={
            "model_id": MODEL_ID,
            "held_out_fold": int(request["run_spec"]["fold"]),
            "best_epoch": best_epoch,
            "best_validation_bce": best_loss,
            "query_rna_used_for_fit": False,
            "held_atac_used_for_fit": False,
            "whole_module_pickle_saved": False,
            "smoke_only": True,
        },
    )


def _load_fitted_model(
    fitted: Path,
    *,
    request: Mapping[str, Any],
    parameters: Mapping[str, Any],
    environment_sha256: str,
) -> tuple[Any, Mapping[str, Any]]:
    import torch

    try:
        manifest = json.loads((fitted / "fitted_model.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RNAATACScPairError("fitted model manifest is invalid") from error
    if (
        not isinstance(manifest, Mapping)
        or manifest.get("run_id") != request["run_id"]
        or manifest.get("model_id") != MODEL_ID
        or manifest.get("upstream_revision") != UPSTREAM_REVISION
        or manifest.get("parameter_sha256") != _canonical_hash(dict(parameters))
        or manifest.get("environment_lock_sha256") != environment_sha256
        or manifest.get("held_atac_used_for_fit") is not False
        or manifest.get("whole_module_pickle_saved") is not False
    ):
        raise RNAATACScPairError("fitted model binding differs")
    state_record = manifest.get("state_dict")
    state_path = fitted / str(state_record.get("path", "")) if isinstance(state_record, Mapping) else fitted / "INVALID"
    if (
        state_path.parent != fitted
        or state_path.is_symlink()
        or not state_path.is_file()
        or state_path.stat().st_size != state_record.get("size_bytes")
        or base._sha256_file(state_path) != state_record.get("sha256")
    ):
        raise RNAATACScPairError("state-dict artifact changed")
    state = torch.load(state_path, map_location="cpu", weights_only=True)
    if not isinstance(state, Mapping) or not all(torch.is_tensor(value) for value in state.values()):
        raise RNAATACScPairError("state-dict is not a tensor mapping")
    records = _state_records(state)
    if _canonical_hash(records) != manifest.get("state_key_shape_dtype_sha256"):
        raise RNAATACScPairError("state-dict key/shape/dtype manifest differs")
    model = _make_model(
        input_dim=int(manifest["input_gene_count"]),
        output_dim=int(manifest["output_peak_count"]),
        hidden_layers=manifest["hidden_layers"],
        dropout=float(manifest["dropout_rate"]),
    )
    model.load_state_dict(state, strict=True)
    model.cuda().eval()
    return model, manifest


def profile_normalize(values: Any, pseudocount: float) -> Any:
    import numpy as np

    scores = np.asarray(values, dtype=np.float64)
    if scores.ndim != 2 or np.any(~np.isfinite(scores)) or np.any(scores < 0):
        raise RNAATACScPairError("profile scores are invalid")
    scores = scores + float(pseudocount)
    totals = scores.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise RNAATACScPairError("profile scores contain an empty row")
    result = scores / totals
    if not np.allclose(result.sum(axis=1), 1.0, rtol=0.0, atol=1e-12):
        raise RNAATACScPairError("profile normalization failed")
    return result


def predict(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import numpy as np
    import torch

    parameters = _validate_request(request, "predict")
    _, environment_sha256 = base._validate_environment(request)
    runtime = _validate_torch_runtime()
    prepared = base._prior_output(request_path, request, "prepare")
    fitted = base._prior_output(request_path, request, "fit")
    query_rna = base._load_sparse(prepared / "query_rna.npz")
    query_fields, query_rows = base._read_tsv(prepared / "query_rows.tsv")
    peak_fields, peak_rows = base._read_tsv(prepared / "selected_peaks.tsv")
    if (
        query_fields != ("cell_id", "donor_id", "lineage", "well_id")
        or len(query_rows) != query_rna.shape[0]
        or peak_fields
        != (
            "selected_index",
            "source_index",
            "peak_id",
            "chromosome",
            "bed_start",
            "bed_end",
        )
        or len(peak_rows) != int(parameters["n_smoke_peaks"])
    ):
        raise RNAATACScPairError("query RNA or feature identities differ")
    model, manifest = _load_fitted_model(
        fitted,
        request=request,
        parameters=parameters,
        environment_sha256=environment_sha256,
    )
    if query_rna.shape[1] != int(manifest["input_gene_count"]):
        raise RNAATACScPairError("query gene axis differs from fitted model")
    cell_scores = np.empty((query_rna.shape[0], len(peak_rows)), dtype=np.float32)
    with torch.no_grad():
        for index in range(query_rna.shape[0]):
            values = torch.from_numpy(
                query_rna[index].toarray().astype(np.float32, copy=False)
            ).cuda()
            _, output_px, _ = model(values)
            cell_scores[index] = output_px.detach().cpu().numpy()[0]
    if np.any(~np.isfinite(cell_scores)) or np.any(cell_scores < 0) or np.any(cell_scores > 1):
        raise RNAATACScPairError("scPair output_px is outside [0, 1]")

    groups: dict[tuple[str, str], list[int]] = {}
    for index, row in enumerate(query_rows):
        groups.setdefault((row["donor_id"], row["lineage"]), []).append(index)
    ordered_groups = sorted(groups)
    grouped_scores = np.vstack(
        [cell_scores[groups[group]].mean(axis=0) for group in ordered_groups]
    )
    profiles = profile_normalize(
        grouped_scores, float(parameters["profile_pseudocount"])
    )

    namespace = str(parameters["join_namespace"])
    prediction_fields = (
        "row_hash",
        "donor_hash",
        "block_hash",
        "stratum",
        "predicted",
    )
    prediction_rows: list[dict[str, str]] = []
    row_id_rows: list[dict[str, str]] = []
    for group_index, (donor, lineage) in enumerate(ordered_groups):
        donor_hash = base.join_hash(namespace, "unit", donor)
        for peak_index, peak in enumerate(peak_rows):
            row_hash = base.join_hash(
                namespace, "row", f"{donor}\0{lineage}\0{peak['peak_id']}"
            )
            prediction_rows.append(
                {
                    "row_hash": row_hash,
                    "donor_hash": donor_hash,
                    "block_hash": base.join_hash(
                        namespace, "block", peak["chromosome"]
                    ),
                    "stratum": lineage,
                    "predicted": format(float(profiles[group_index, peak_index]), ".17g"),
                }
            )
            row_id_rows.append({"row_hash": row_hash, "donor_hash": donor_hash})
    prediction_rows.sort(key=lambda row: row["row_hash"])
    row_id_rows.sort(key=lambda row: row["row_hash"])
    prediction_path = output / "predictions.tsv"
    row_ids_path = output / "row_ids.tsv"
    _write_tsv(prediction_path, prediction_fields, prediction_rows)
    _write_tsv(row_ids_path, ("row_hash", "donor_hash"), row_id_rows)
    standardized = {
        **base._artifact_record(prediction_path, relative_to=output),
        "media_type": "text/tab-separated-values",
        "role": f"standardized_prediction_table:{TASK_ID}",
    }
    row_artifact = {
        **base._artifact_record(row_ids_path, relative_to=output),
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
        "bundle_id": f"{MODEL_ID}-{str(request['run_id'])[:16]}",
        "run_id": request["run_id"],
        "task_id": TASK_ID,
        "model_id": MODEL_ID,
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
            "query_cell_count": len(query_rows),
            "query_stratum_count": len(ordered_groups),
            "selected_peak_count": len(peak_rows),
            "prediction_scale": "donor_lineage_mean_output_px_then_depth_free_profile",
            "training_loss_scale": "output_result_with_source_RNA_inferred_scalar",
            "output_px_distinct_from_output_result": True,
            "profile_pseudocount": float(parameters["profile_pseudocount"]),
            "inference_batch_size": int(parameters["inference_batch_size"]),
            "query_atac_tensor_created": False,
            "held_atac_input_exposed": False,
            "observed_atac_exported": False,
            "profile_fixture_passed": True,
            "query_composition_invariant_by_construction": True,
            "source_feature_inventory_includes_held_donor_peak_discovery": True,
            "runtime": runtime,
            "smoke_only": True,
            "champion_claim_allowed": False,
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
            "model_id": MODEL_ID,
            "held_out_fold": int(request["run_spec"]["fold"]),
            "prediction_row_count": len(prediction_rows),
            "prediction_donor_count": len(
                {row["donor_hash"] for row in prediction_rows}
            ),
            "query_atac_tensor_created": False,
            "held_atac_input_exposed": False,
            "observed_atac_exported": False,
            "profile_fixture_passed": True,
            "smoke_only": True,
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
        raise RNAATACScPairError("adapter request is invalid JSON") from error
    if not isinstance(request, Mapping):
        raise RNAATACScPairError("adapter request must be an object")
    if arguments.output.exists():
        raise RNAATACScPairError(f"adapter output already exists: {arguments.output}")
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

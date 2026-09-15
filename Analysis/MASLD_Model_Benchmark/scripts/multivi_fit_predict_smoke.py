#!/usr/bin/env python3
"""Fit five donor-held MultiVI folds and predict ATAC from held-donor RNA only."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


MODEL_ID = "multivi"
TASK_ID = "rna_conditioned_atac"
DATASET_ID = "gse296875"
JOIN_NAMESPACE = "gse296875:rna_conditioned_atac:donor_outer:smoke:v1"
PREDICTION_FIELDS = (
    "row_hash",
    "donor_hash",
    "block_hash",
    "stratum",
    "predicted",
)


class MultiVIFitPredictError(ValueError):
    """Raised when the donor-held MultiVI execution requirement differs."""


def canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def canonical_hash(value: Any) -> str:
    return sha256(canonical_json(value).encode()).hexdigest()


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def join_hash(namespace: str, kind: str, identifier: str) -> str:
    return sha256(f"{namespace}\0{kind}\0{identifier}".encode()).hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise MultiVIFitPredictError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def write_tsv(
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
        writer.writerows(rows)


def artifact_record(path: Path, relative_to: Path, role: str) -> dict[str, Any]:
    return {
        "path": path.relative_to(relative_to).as_posix(),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
        "media_type": "text/tab-separated-values",
        "role": role,
    }


def validate_modality_states(
    row_sums: Sequence[float], states: Sequence[str]
) -> tuple[bool, ...]:
    if len(row_sums) != len(states) or not row_sums:
        raise MultiVIFitPredictError("modality state and row-sum axes differ")
    result: list[bool] = []
    for value, state in zip(row_sums, states, strict=True):
        if state not in {"observed", "structurally_missing"}:
            raise MultiVIFitPredictError("unsupported modality state")
        if not math.isfinite(float(value)) or value < 0:
            raise MultiVIFitPredictError("modality row sum is invalid")
        if state == "observed" and value <= 0:
            raise MultiVIFitPredictError("observed modality row is empty")
        if state == "structurally_missing" and value != 0:
            raise MultiVIFitPredictError("missing modality row contains values")
        result.append(state == "observed")
    if result != [value > 0 for value in row_sums]:
        raise MultiVIFitPredictError("explicit state differs from runtime mask")
    return tuple(result)


def normalize_profiles(values: Any, pseudocount: float = 1.0e-8) -> Any:
    import numpy as np

    matrix = np.asarray(values, dtype=np.float64)
    if (
        matrix.ndim != 2
        or matrix.shape[0] == 0
        or matrix.shape[1] == 0
        or not math.isfinite(pseudocount)
        or pseudocount <= 0
        or np.any(~np.isfinite(matrix))
        or np.any(matrix < 0)
    ):
        raise MultiVIFitPredictError("profile matrix is invalid")
    matrix = matrix + pseudocount
    matrix /= matrix.sum(axis=1, keepdims=True)
    if np.any(matrix <= 0) or not np.allclose(
        matrix.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-12
    ):
        raise MultiVIFitPredictError("normalized profiles are invalid")
    return matrix


def make_model(torch: Any, n_obs: int, n_genes: int, n_peaks: int) -> Any:
    from scvi.module import MULTIVAE

    return MULTIVAE(
        n_input_regions=n_peaks,
        n_input_genes=n_genes,
        n_input_proteins=0,
        modality_weights="equal",
        modality_penalty="Jeffreys",
        n_batch=1,
        n_obs=n_obs,
        n_labels=1,
        gene_likelihood="poisson",
        gene_dispersion="gene",
        n_hidden=128,
        n_latent=16,
        n_layers_encoder=1,
        n_layers_decoder=1,
        dropout_rate=0.1,
        region_factors=True,
        use_batch_norm="none",
        use_layer_norm="both",
        latent_distribution="normal",
    ).cuda()


def tensors(torch: Any, rna: Any, atac: Any, indices: Any) -> dict[str, Any]:
    from scvi import REGISTRY_KEYS

    n_rows = rna.shape[0]
    return {
        REGISTRY_KEYS.X_KEY: rna,
        REGISTRY_KEYS.ATAC_X_KEY: atac,
        REGISTRY_KEYS.BATCH_KEY: torch.zeros(
            (n_rows, 1), dtype=torch.long, device=rna.device
        ),
        REGISTRY_KEYS.INDICES_KEY: indices.to(device=rna.device, dtype=torch.long),
        REGISTRY_KEYS.LABELS_KEY: torch.zeros(
            (n_rows, 1), dtype=torch.long, device=rna.device
        ),
    }


def dense_batch(torch: Any, matrix: Any, rows: Any, *, binary: bool) -> Any:
    import numpy as np

    values = matrix[rows].toarray().astype(np.float32, copy=False)
    if binary:
        values = (values > 0).astype(np.float32, copy=False)
    return torch.from_numpy(values).cuda(non_blocking=False)


def state_dict_hash(state: Mapping[str, Any]) -> str:
    digest = sha256()
    for name in sorted(state):
        tensor = state[name].detach().contiguous().cpu()
        digest.update(name.encode())
        digest.update(str(tensor.dtype).encode())
        digest.update(canonical_json(list(tensor.shape)).encode())
        digest.update(tensor.numpy().tobytes(order="C"))
    return digest.hexdigest()


def train_fold(
    torch: Any,
    model: Any,
    rna: Any,
    atac: Any,
    *,
    seed: int,
    epochs: int,
    batch_size: int,
) -> list[float]:
    import numpy as np

    optimizer = torch.optim.AdamW(model.parameters(), lr=1.0e-3, weight_decay=1.0e-6)
    losses: list[float] = []
    generator = np.random.default_rng(seed)
    for _epoch in range(epochs):
        permutation = generator.permutation(rna.shape[0])
        epoch_loss = 0.0
        epoch_rows = 0
        model.train()
        for start in range(0, len(permutation), batch_size):
            positions = permutation[start : start + batch_size]
            rna_batch = dense_batch(torch, rna, positions, binary=False)
            atac_batch = dense_batch(torch, atac, positions, binary=True)
            optimizer.zero_grad(set_to_none=True)
            _inference, _generative, output = model(
                tensors=tensors(
                    torch,
                    rna_batch,
                    atac_batch,
                    torch.from_numpy(positions.astype(np.int64, copy=False)),
                ),
                compute_loss=True,
            )
            if not torch.isfinite(output.loss):
                raise MultiVIFitPredictError("MultiVI training loss is non-finite")
            output.loss.backward()
            if not any(
                parameter.grad is not None and torch.isfinite(parameter.grad).all()
                for parameter in model.parameters()
            ):
                raise MultiVIFitPredictError("MultiVI finite gradients are absent")
            optimizer.step()
            count = len(positions)
            epoch_loss += float(output.loss.detach().cpu()) * count
            epoch_rows += count
        losses.append(epoch_loss / epoch_rows)
    if not losses or any(not math.isfinite(value) for value in losses):
        raise MultiVIFitPredictError("MultiVI loss history is invalid")
    return losses


def predict_query(torch: Any, model: Any, query_rna: Any, batch_size: int) -> Any:
    import numpy as np

    profiles: list[Any] = []
    model.eval()
    with torch.inference_mode():
        for start in range(0, query_rna.shape[0], batch_size):
            stop = min(start + batch_size, query_rna.shape[0])
            positions = np.arange(start, stop, dtype=np.int64)
            rna_batch = dense_batch(torch, query_rna, positions, binary=False)
            atac_batch = torch.zeros(
                (len(positions), model.n_input_regions),
                dtype=torch.float32,
                device="cuda",
            )
            validate_modality_states(
                atac_batch.sum(1).cpu().tolist(),
                ["structurally_missing"] * len(positions),
            )
            _inference, generative = model(
                tensors=tensors(
                    torch,
                    rna_batch,
                    atac_batch,
                    torch.from_numpy(positions),
                ),
                generative_kwargs={"use_z_mean": True},
                compute_loss=False,
            )
            profile = generative["p"]
            if (
                profile.shape != atac_batch.shape
                or not torch.isfinite(profile).all()
                or torch.any(profile < 0)
                or torch.any(profile > 1)
            ):
                raise MultiVIFitPredictError("MultiVI query profile is invalid")
            profiles.append(profile.cpu().numpy())
    return normalize_profiles(np.vstack(profiles))


def export_bundle(
    *,
    fold: int,
    run_id: str,
    query_rows: Sequence[Mapping[str, str]],
    peaks: Sequence[Mapping[str, str]],
    cell_profiles: Any,
    output: Path,
    seed: int,
    training_nuclei: int,
    state_hash: str,
) -> Path:
    import numpy as np

    groups: dict[tuple[str, str], list[int]] = {}
    for index, row in enumerate(query_rows):
        groups.setdefault((row["donor_id"], row["lineage"]), []).append(index)
    ordered = sorted(groups)
    profiles = normalize_profiles(
        np.vstack(
            [cell_profiles[groups[group]].mean(axis=0, dtype=np.float64) for group in ordered]
        )
    )
    prediction_rows: list[dict[str, str]] = []
    row_id_rows: list[dict[str, str]] = []
    for group_index, (donor, lineage) in enumerate(ordered):
        donor_hash = join_hash(JOIN_NAMESPACE, "unit", donor)
        for peak_index, peak in enumerate(peaks):
            row_hash = join_hash(
                JOIN_NAMESPACE,
                "row",
                f"{donor}\0{lineage}\0{peak['peak_id']}",
            )
            prediction_rows.append(
                {
                    "row_hash": row_hash,
                    "donor_hash": donor_hash,
                    "block_hash": join_hash(
                        JOIN_NAMESPACE, "block", peak["chromosome"]
                    ),
                    "stratum": lineage,
                    "predicted": format(float(profiles[group_index, peak_index]), ".17g"),
                }
            )
            row_id_rows.append({"row_hash": row_hash, "donor_hash": donor_hash})
    prediction_rows.sort(key=lambda row: row["row_hash"])
    row_id_rows.sort(key=lambda row: row["row_hash"])
    write_tsv(output / "predictions.tsv", PREDICTION_FIELDS, prediction_rows)
    write_tsv(output / "row_ids.tsv", ("row_hash", "donor_hash"), row_id_rows)
    prediction = artifact_record(
        output / "predictions.tsv",
        output,
        f"standardized_prediction_table:{TASK_ID}",
    )
    row_ids = artifact_record(
        output / "row_ids.tsv", output, f"prediction_row_ids:{TASK_ID}"
    )
    bundle = {
        "schema_version": "masld-bench-prediction-bundle-v1",
        "bundle_id": f"{MODEL_ID}-{run_id[:16]}",
        "run_id": run_id,
        "task_id": TASK_ID,
        "model_id": MODEL_ID,
        "dataset_ids": [DATASET_ID],
        "split_id": "donor_outer",
        "artifacts": [prediction, row_ids],
        "standardized_table": prediction,
        "row_ids": row_ids,
        "n_predictions": len(prediction_rows),
        "row_id_field": "row_hash",
        "unit_id_field": "donor_hash",
        "unit_id_namespace": JOIN_NAMESPACE,
        "biological_unit": "donor",
        "table_schema_sha256": canonical_hash(
            {"format": "tsv", "fields": list(PREDICTION_FIELDS)}
        ),
        "source_join_key_sha256": canonical_hash(
            {
                "task_id": TASK_ID,
                "dataset_ids": [DATASET_ID],
                "split_id": "donor_outer",
                "row_id_field": "row_hash",
                "unit_id_field": "donor_hash",
                "unit_id_namespace": JOIN_NAMESPACE,
                "biological_unit": "donor",
            }
        ),
        "format_version": "tsv-v1",
        "missing_state": "observed",
        "metadata": {
            "held_out_fold": fold,
            "seed": seed,
            "training_nuclei": training_nuclei,
            "query_nuclei": len(query_rows),
            "query_stratum_count": len(ordered),
            "selected_peak_count": len(peaks),
            "prediction_scale": "donor_lineage_mean_depth_free_multinomial_peak_composition",
            "query_rna_state": "observed",
            "query_atac_state": "structurally_missing",
            "query_atac_tensor_created": True,
            "query_atac_tensor_nonzero_values": 0,
            "held_atac_input_exposed": False,
            "observed_atac_exported": False,
            "source_hdf5_available_to_fit_predict": False,
            "same_nucleus_training": True,
            "state_dict_sha256": state_hash,
            "smoke_only": True,
            "champion_claim_allowed": False,
        },
    }
    path = output / "prediction_bundle.json"
    path.write_text(canonical_json(bundle) + "\n", encoding="utf-8")
    return path


def run(prepared: Path, output: Path, seed: int, epochs: int) -> dict[str, Any]:
    import numpy as np
    from scipy import sparse
    import scvi
    import torch

    if output.exists() or seed < 0 or epochs < 1:
        raise MultiVIFitPredictError("output, seed, or epoch contract differs")
    if scvi.__version__ != "1.5.0.post1" or torch.__version__ != "2.6.0+cu124":
        raise MultiVIFitPredictError("MultiVI runtime version differs")
    if torch.cuda.device_count() != 1 or "L40S" not in torch.cuda.get_device_name(0):
        raise MultiVIFitPredictError("exactly one L40S is required")
    receipt = json.loads((prepared / "prepare_receipt.json").read_text())
    if (
        receipt.get("status") != "pass"
        or receipt.get("pairing_topology") != "same_nucleus"
        or receipt.get("held_atac_exported_from_prepare") is not False
        or receipt.get("fit_predict_source_hdf5_required") is not False
        or receipt.get("outcomes_read") is not False
    ):
        raise MultiVIFitPredictError("prepared same-nucleus contract differs")
    output.mkdir(parents=True, mode=0o750)
    fold_receipts: dict[str, Any] = {}
    torch.use_deterministic_algorithms(True)
    for fold in range(5):
        fold_seed = seed + fold
        torch.manual_seed(fold_seed)
        torch.cuda.manual_seed_all(fold_seed)
        root = prepared / f"fold_{fold}"
        train_rna = sparse.load_npz(root / "training_rna.npz").tocsr()
        train_atac = sparse.load_npz(root / "training_atac.npz").tocsr()
        query_rna = sparse.load_npz(root / "query_rna.npz").tocsr()
        train_fields, train_rows = read_tsv(root / "training_rows.tsv")
        query_fields, query_rows = read_tsv(root / "query_rows.tsv")
        peak_fields, peaks = read_tsv(root / "selected_peaks.tsv")
        if (
            train_fields != ("cell_id", "donor_id", "lineage", "rna_state", "atac_state")
            or query_fields != train_fields
            or peak_fields
            != (
                "selected_index",
                "source_index",
                "peak_id",
                "chromosome",
                "bed_start",
                "bed_end",
            )
            or len(train_rows) != train_rna.shape[0]
            or len(query_rows) != query_rna.shape[0]
            or train_rna.shape[1] != 2000
            or train_atac.shape != (train_rna.shape[0], len(peaks))
            or query_rna.shape[1] != train_rna.shape[1]
            or len(peaks) != 10_000
            or set(row["rna_state"] for row in train_rows + query_rows) != {"observed"}
            or set(row["atac_state"] for row in train_rows) != {"observed"}
            or set(row["atac_state"] for row in query_rows) != {"structurally_missing"}
        ):
            raise MultiVIFitPredictError("prepared fold axes or states differ")
        validate_modality_states(
            np.asarray(train_rna.sum(axis=1)).ravel().tolist(),
            [row["rna_state"] for row in train_rows],
        )
        validate_modality_states(
            np.asarray(train_atac.sum(axis=1)).ravel().tolist(),
            [row["atac_state"] for row in train_rows],
        )
        validate_modality_states(
            np.asarray(query_rna.sum(axis=1)).ravel().tolist(),
            [row["rna_state"] for row in query_rows],
        )
        if {row["donor_id"] for row in train_rows} & {
            row["donor_id"] for row in query_rows
        }:
            raise MultiVIFitPredictError("donor crosses fit and query")
        model = make_model(
            torch, train_rna.shape[0], train_rna.shape[1], train_atac.shape[1]
        )
        losses = train_fold(
            torch,
            model,
            train_rna,
            train_atac,
            seed=fold_seed,
            epochs=epochs,
            batch_size=64,
        )
        state = {name: value.detach().cpu() for name, value in model.state_dict().items()}
        state_hash = state_dict_hash(state)
        fold_output = output / f"fold_{fold}"
        fold_output.mkdir(mode=0o750)
        state_path = fold_output / "state_dict.pt"
        torch.save(state, state_path)
        restored = make_model(
            torch, train_rna.shape[0], train_rna.shape[1], train_atac.shape[1]
        )
        restored.load_state_dict(
            torch.load(state_path, map_location="cpu", weights_only=True), strict=True
        )
        restored.cuda().eval()
        if state_dict_hash(restored.state_dict()) != state_hash:
            raise MultiVIFitPredictError("strict state-dict resume differs")
        cell_profiles = predict_query(torch, restored, query_rna, batch_size=128)
        repeat = predict_query(torch, restored, query_rna, batch_size=128)
        if not np.array_equal(cell_profiles, repeat):
            raise MultiVIFitPredictError("repeated query prediction differs")
        run_id = canonical_hash(
            {
                "model_id": MODEL_ID,
                "fold": fold,
                "seed": fold_seed,
                "epochs": epochs,
                "prepared_receipt_sha256": sha256_file(
                    prepared / "prepare_receipt.json"
                ),
                "state_dict_sha256": state_hash,
            }
        )
        bundle = export_bundle(
            fold=fold,
            run_id=run_id,
            query_rows=query_rows,
            peaks=peaks,
            cell_profiles=cell_profiles,
            output=fold_output,
            seed=fold_seed,
            training_nuclei=len(train_rows),
            state_hash=state_hash,
        )
        fold_receipts[str(fold)] = {
            "run_id": run_id,
            "training_nuclei": len(train_rows),
            "query_nuclei": len(query_rows),
            "training_donors": len({row["donor_id"] for row in train_rows}),
            "query_donors": len({row["donor_id"] for row in query_rows}),
            "epochs": epochs,
            "first_epoch_loss": losses[0],
            "last_epoch_loss": losses[-1],
            "state_dict_sha256": state_hash,
            "state_dict_file_sha256": sha256_file(state_path),
            "prediction_bundle": bundle.relative_to(output).as_posix(),
            "prediction_bundle_sha256": sha256_file(bundle),
            "repeat_prediction_bit_identical": True,
            "query_atac_state": "structurally_missing",
            "held_atac_read": False,
        }
        del model, restored, state, cell_profiles, repeat
        torch.cuda.empty_cache()
    final = {
        "schema_version": "masld-bench-multivi-smoke-fit-predict-v1",
        "status": "pass",
        "model_id": MODEL_ID,
        "task_id": TASK_ID,
        "dataset_id": DATASET_ID,
        "pairing_topology": "same_nucleus",
        "outer_unit": "donor",
        "folds": fold_receipts,
        "seed": seed,
        "epochs": epochs,
        "runtime": {
            "scvi_tools": scvi.__version__,
            "torch": torch.__version__,
            "cuda_build": torch.version.cuda,
            "gpu": torch.cuda.get_device_name(0),
        },
        "source_hdf5_available_to_fit_predict": False,
        "query_atac_read": False,
        "query_atac_state": "structurally_missing",
        "query_atac_tensor_nonzero_values": 0,
        "outcomes_read": False,
        "metrics_calculated_by_adapter": False,
        "smoke_only": True,
        "champion_claim_allowed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(final, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(canonical_json(final))
    return final


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=1103)
    parser.add_argument("--epochs", type=int, default=50)
    arguments = parser.parse_args()
    run(arguments.prepared, arguments.output, arguments.seed, arguments.epochs)


if __name__ == "__main__":
    main()

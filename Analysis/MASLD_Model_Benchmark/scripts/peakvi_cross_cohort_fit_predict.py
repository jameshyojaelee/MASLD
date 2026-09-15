#!/usr/bin/env python3
"""Fit source-only PeakVI models and export blind GSE281367 projections."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-peakvi-cross-cohort-fit-predict-v1"
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage")
ROTATIONS = {
    "valid_context_test_target": ("valid", "test"),
    "test_context_valid_target": ("test", "valid"),
}


class PeakVICrossCohortError(ValueError):
    """Raised when source fitting or blind prediction does not meet its requirements."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise PeakVICrossCohortError(f"JSON object required: {path}")
    return value


def read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise PeakVICrossCohortError(f"TSV fields differ: {path}")
        return list(reader)


def verify_bound_tree(
    root: Path, record: Mapping[str, Any], artifact_class: str
) -> Path:
    path = (root / str(record.get("path", ""))).resolve(strict=True)
    try:
        path.relative_to(root)
    except ValueError as error:
        raise PeakVICrossCohortError(f"{artifact_class} escapes benchmark root") from error
    manifest = verify_frozen_tree(path)
    if (
        digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256")
        or manifest.get("metadata", {}).get("artifact_class") != artifact_class
    ):
        raise PeakVICrossCohortError(f"{artifact_class} authority differs")
    return path


def feature_rows(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(
        path,
        ("feature_index", "window_id", "contig", "start", "end"),
    )
    if (
        len(rows) != 16_000
        or [int(row["feature_index"]) for row in rows] != list(range(16_000))
        or any(int(row["end"]) - int(row["start"]) != 1000 for row in rows)
        or len({row["window_id"] for row in rows}) != 16_000
    ):
        raise PeakVICrossCohortError("fixed-window feature axis differs")
    return rows


def cell_rows(path: Path, donors: int) -> list[dict[str, str]]:
    rows = read_tsv(
        path,
        (
            "cell_index",
            "cell_hash",
            "donor_id",
            "lineage_id",
            "outer_fold",
            "selection_rank",
            "available_nuclei",
        ),
    )
    if (
        [int(row["cell_index"]) for row in rows] != list(range(len(rows)))
        or len({row["cell_hash"] for row in rows}) != len(rows)
        or len({row["donor_id"] for row in rows}) != donors
        or {row["lineage_id"] for row in rows} != set(LINEAGES)
        or len({(row["donor_id"], row["lineage_id"]) for row in rows}) != donors * 4
    ):
        raise PeakVICrossCohortError("cell axis differs")
    return rows


def validate_config(root: Path, config_path: Path) -> tuple[dict[str, Any], dict[str, Path]]:
    config = read_json(config_path.resolve(strict=True))
    if (
        config.get("schema_version") != SCHEMA
        or config.get("source_dataset_id") != "gse296875"
        or config.get("target_dataset_id") != "gse281367"
        or config.get("seeds") != [1103, 2207, 3301]
        or config.get("model", {}).get("native_output_family") != "joint_representation"
        or config.get("model", {}).get("projected_output_family")
        != "masked_accessibility_count"
        or config.get("model", {}).get("target_query_adaptation_allowed") is not False
        or config.get("model", {}).get("dimensions")
        != {"hidden": 128, "latent": 16, "encoder_layers": 2, "decoder_layers": 2}
        or config.get("training", {}).get("peakvi_epochs") != 30
        or config.get("training", {}).get("decoder_epochs") != 20
        or config.get("training", {}).get("batch_size") != 128
        or config.get("training", {}).get("query_transform_scope")
        != "source_training_only"
    ):
        raise PeakVICrossCohortError("campaign identity, seeds, or model contract differs")
    for field in (
        "condition_values_read",
        "phenotype_values_read",
        "development_outcomes_read",
        "sealed_data_read",
        "metrics_calculated",
        "target_role_atac_consumed",
        "target_query_adaptation_performed",
        "target_query_calibration_performed",
        "scored_target_contigs_enter_query_transform",
        "missing_as_zero",
    ):
        if config.get("firewall", {}).get(field) is not False:
            raise PeakVICrossCohortError(f"campaign firewall opened: {field}")
    if config.get("firewall", {}).get("prediction_commit_before_scoring_required") is not True:
        raise PeakVICrossCohortError("blind prediction commit is not required")
    paths = {
        "terminal_gate": verify_bound_tree(
            root,
            config["terminal_gate"],
            "gse296875_observed_multiome_specialist_gate_readiness",
        ),
        "runtime_probe": verify_bound_tree(
            root,
            config["runtime_probe"],
            "multivi_peakvi_no_outcome_runtime_probe",
        ),
        "source": verify_bound_tree(
            root, config["source_input"], "peakvi_gse296875_source_cell_input"
        ),
        "target_valid": verify_bound_tree(
            root,
            config["target_inputs"]["valid"],
            "peakvi_gse281367_target_cell_context",
        ),
        "target_test": verify_bound_tree(
            root,
            config["target_inputs"]["test"],
            "peakvi_gse281367_target_cell_context",
        ),
        "axis": verify_bound_tree(
            root,
            config["exchange_axis"],
            "gse281367_label_free_atac_exchange_axis",
        ),
    }
    if read_json(paths["target_valid"] / "contract.json").get("context_role") != "valid":
        raise PeakVICrossCohortError("valid context artifact differs")
    if read_json(paths["target_test"] / "contract.json").get("context_role") != "test":
        raise PeakVICrossCohortError("test context artifact differs")
    model_path = Path(str(config["model_manifest"]["path"])).resolve(strict=True)
    if (
        model_path.is_symlink()
        or not model_path.is_file()
        or digest(model_path) != config["model_manifest"].get("sha256")
    ):
        raise PeakVICrossCohortError("PeakVI model-manifest binding differs")
    paths["model_manifest"] = model_path
    for name, record in config.get("implementation", {}).items():
        path = (root / str(record.get("path", ""))).resolve(strict=True)
        try:
            path.relative_to(root)
        except ValueError as error:
            raise PeakVICrossCohortError(f"implementation escapes root: {name}") from error
        if path.is_symlink() or not path.is_file() or digest(path) != record.get("sha256"):
            raise PeakVICrossCohortError(f"implementation binding differs: {name}")
    return config, paths


def make_peakvae(n_regions: int, dimensions: Mapping[str, int]) -> Any:
    from scvi.module import PEAKVAE

    return PEAKVAE(
        n_input_regions=n_regions,
        n_batch=1,
        n_hidden=int(dimensions["hidden"]),
        n_latent=int(dimensions["latent"]),
        n_layers_encoder=int(dimensions["encoder_layers"]),
        n_layers_decoder=int(dimensions["decoder_layers"]),
        dropout_rate=0.1,
        model_depth=True,
        region_factors=True,
        use_batch_norm="none",
        use_layer_norm="both",
        latent_distribution="normal",
    ).cuda()


def peak_tensors(torch: Any, values: Any) -> dict[str, Any]:
    from scvi import REGISTRY_KEYS

    return {
        REGISTRY_KEYS.X_KEY: values,
        REGISTRY_KEYS.BATCH_KEY: torch.zeros(
            (values.shape[0], 1), dtype=torch.long, device=values.device
        ),
    }


def dense_binary_batch(torch: Any, matrix: Any, positions: Any) -> Any:
    import numpy as np

    values = matrix[positions].toarray().astype(np.float32, copy=False)
    return torch.from_numpy(values > 0).to(device="cuda", dtype=torch.float32)


def fit_peakvae(
    torch: Any,
    model: Any,
    matrix: Any,
    *,
    seed: int,
    epochs: int,
    batch_size: int,
) -> list[float]:
    import numpy as np

    optimizer = torch.optim.AdamW(model.parameters(), lr=1.0e-3, weight_decay=1.0e-6)
    generator = np.random.default_rng(seed)
    losses: list[float] = []
    for _epoch in range(epochs):
        permutation = generator.permutation(matrix.shape[0])
        total = 0.0
        rows = 0
        model.train()
        for start in range(0, len(permutation), batch_size):
            positions = permutation[start : start + batch_size]
            values = dense_binary_batch(torch, matrix, positions)
            optimizer.zero_grad(set_to_none=True)
            _inference, _generative, output = model(
                tensors=peak_tensors(torch, values), compute_loss=True
            )
            if not torch.isfinite(output.loss):
                raise PeakVICrossCohortError("PeakVI loss is non-finite")
            output.loss.backward()
            if not all(
                parameter.grad is None or torch.isfinite(parameter.grad).all()
                for parameter in model.parameters()
            ):
                raise PeakVICrossCohortError("PeakVI gradient is non-finite")
            optimizer.step()
            total += float(output.loss.detach().cpu())
            rows += len(positions)
        losses.append(total / rows)
    if not losses or any(not math.isfinite(value) for value in losses):
        raise PeakVICrossCohortError("PeakVI loss history differs")
    return losses


def encode(
    torch: Any, model: Any, matrix: Any, *, batch_size: int, latent_dimension: int
) -> Any:
    import numpy as np

    results: list[Any] = []
    model.eval()
    with torch.inference_mode():
        for start in range(0, matrix.shape[0], batch_size):
            positions = np.arange(start, min(start + batch_size, matrix.shape[0]))
            values = dense_binary_batch(torch, matrix, positions)
            inference, _generative = model(
                tensors=peak_tensors(torch, values),
                generative_kwargs={"use_z_mean": True},
                compute_loss=False,
            )
            means = inference["qz"].loc
            if means.shape != (len(positions), latent_dimension) or not torch.isfinite(means).all():
                raise PeakVICrossCohortError("PeakVI posterior means differ")
            results.append(means.cpu().numpy())
    return np.vstack(results).astype(np.float32, copy=False)


def parameter_hash(state: Mapping[str, Any]) -> str:
    value = sha256()
    for name, tensor in sorted(state.items()):
        contiguous = tensor.detach().cpu().contiguous()
        value.update(name.encode())
        value.update(b"\0")
        value.update(str(contiguous.dtype).encode())
        value.update(b"\0")
        value.update(json.dumps(list(contiguous.shape)).encode())
        value.update(b"\0")
        value.update(contiguous.numpy().tobytes(order="C"))
    return value.hexdigest()


def fit_decoder(
    torch: Any,
    source_latent: Any,
    source_targets: Any,
    *,
    seed: int,
    epochs: int,
    batch_size: int,
) -> tuple[Any, list[float]]:
    import numpy as np

    torch.manual_seed(seed)
    decoder = torch.nn.Linear(source_latent.shape[1], source_targets.shape[1]).cuda()
    prevalence = np.asarray((source_targets > 0).mean(axis=0)).ravel()
    prevalence = np.clip(prevalence, 1.0e-4, 1.0 - 1.0e-4)
    with torch.no_grad():
        decoder.weight.zero_()
        decoder.bias.copy_(
            torch.from_numpy(np.log(prevalence / (1.0 - prevalence)).astype(np.float32)).cuda()
        )
    optimizer = torch.optim.AdamW(decoder.parameters(), lr=1.0e-3, weight_decay=1.0e-6)
    criterion = torch.nn.BCEWithLogitsLoss(reduction="mean")
    generator = np.random.default_rng(seed)
    losses: list[float] = []
    latent = torch.from_numpy(source_latent).cuda()
    for _epoch in range(epochs):
        permutation = generator.permutation(source_targets.shape[0])
        total = 0.0
        rows = 0
        decoder.train()
        for start in range(0, len(permutation), batch_size):
            positions = permutation[start : start + batch_size]
            targets = dense_binary_batch(torch, source_targets, positions)
            optimizer.zero_grad(set_to_none=True)
            loss = criterion(decoder(latent[positions]), targets)
            if not torch.isfinite(loss):
                raise PeakVICrossCohortError("source-only decoder loss is non-finite")
            loss.backward()
            if not all(
                parameter.grad is None or torch.isfinite(parameter.grad).all()
                for parameter in decoder.parameters()
            ):
                raise PeakVICrossCohortError("source-only decoder gradient is non-finite")
            optimizer.step()
            total += float(loss.detach().cpu()) * len(positions)
            rows += len(positions)
        losses.append(total / rows)
    if not losses or any(not math.isfinite(value) for value in losses):
        raise PeakVICrossCohortError("decoder loss history differs")
    return decoder, losses


def donor_lineage_predictions(
    torch: Any,
    decoder: Any,
    target_latent: Any,
    target_cells: Sequence[Mapping[str, str]],
    *,
    batch_size: int,
) -> Any:
    import numpy as np

    donors = sorted({row["donor_id"] for row in target_cells})
    if donors != [f"Z{index:02d}" for index in range(1, 13)]:
        raise PeakVICrossCohortError("target donor axis differs")
    donor_index = {value: index for index, value in enumerate(donors)}
    lineage_index = {value: index for index, value in enumerate(LINEAGES)}
    sums = np.zeros((12, 4, decoder.out_features), dtype=np.float64)
    counts = np.zeros((12, 4), dtype=np.int64)
    latent = torch.from_numpy(target_latent).cuda()
    decoder.eval()
    with torch.inference_mode():
        for start in range(0, len(target_cells), batch_size):
            stop = min(start + batch_size, len(target_cells))
            probabilities = torch.sigmoid(decoder(latent[start:stop])).cpu().numpy()
            for offset, row in enumerate(target_cells[start:stop]):
                group = (donor_index[row["donor_id"]], lineage_index[row["lineage_id"]])
                sums[group] += probabilities[offset]
                counts[group] += 1
    if np.any(counts == 0):
        raise PeakVICrossCohortError("target donor-lineage group is empty")
    result = sums / counts[..., None]
    if result.shape != (12, 4, 16_000) or not np.all(np.isfinite(result)):
        raise PeakVICrossCohortError("projected prediction geometry differs")
    if np.any(result < 0) or np.any(result > 1):
        raise PeakVICrossCohortError("projected accessibility probability differs")
    return result.astype(np.float32)


def array_record(root: Path, path: Path, dtype: str) -> dict[str, Any]:
    import numpy as np

    value = np.load(path, mmap_mode="r", allow_pickle=False)
    if str(value.dtype) != dtype:
        raise PeakVICrossCohortError(f"array dtype differs: {path}")
    return {
        "path": path.relative_to(root).as_posix(),
        "sha256": digest(path),
        "dtype": dtype,
        "shape": list(value.shape),
    }


def make_prediction_manifest(
    *,
    config: Mapping[str, Any],
    config_sha256: str,
    model_manifest: Path,
    rotation_id: str,
    predictions: Path,
    missing_state: Path,
    output: Path,
) -> dict[str, Any]:
    context_role, target_role = ROTATIONS[rotation_id]
    return {
        "schema_version": "masld-bench-atac-transport-prediction-bundle-v1",
        "bundle_id": f"peakvi-source-fit-{rotation_id}-{config_sha256[:16]}",
        "model_id": "peakvi_source_fit_fixed_window_decoder",
        "input_regime": "observed_atac",
        "output_family": "masked_accessibility_count",
        "prediction_layout": "donor_lineage",
        "rotation_id": rotation_id,
        "axis_artifacts_sha256": config["exchange_axis"]["artifacts_sha256"],
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse281367",
        "model_manifest": {
            "path": str(model_manifest),
            "sha256": digest(model_manifest),
        },
        "arrays": {
            "predictions": array_record(output, predictions, "float32"),
            "missing_state": array_record(output, missing_state, "uint8"),
        },
        "native_family_contract": {
            "native_output_family": "joint_representation",
            "native_target_latents_preserved": True,
            "projected_output_family": "masked_accessibility_count",
            "projection": "source_only_logistic_fixed_window_decoder",
            "direct_native_to_profile_ranking_allowed": False,
        },
        "prediction_scale": "mean_predicted_single_nucleus_binary_accessibility_probability",
        "seeds": list(config["seeds"]),
        "evidence_contract": {
            "condition_labels_consumed": False,
            "phenotype_values_consumed": False,
            "development_outcomes_consumed": False,
            "sealed_data_consumed": False,
            "target_role_atac_consumed": False,
            "target_role_contigs_entered_query_transform": False,
            "missing_as_zero": False,
            "scored_target_role": target_role,
            "observed_atac_context_role": context_role,
            "query_atac_consumed": True,
            "receptive_field_bp": 1000,
            "query_transform_scope": "source_training_only",
        },
        "metrics_calculated": False,
        "champion_claim_allowed": False,
    }


def fit_rotation(
    *,
    torch: Any,
    config: Mapping[str, Any],
    config_sha256: str,
    paths: Mapping[str, Path],
    rotation_id: str,
    output: Path,
) -> dict[str, Any]:
    import numpy as np
    from scipy import sparse

    context_role, target_role = ROTATIONS[rotation_id]
    target_input = paths[f"target_{context_role}"]
    source_cells = cell_rows(paths["source"] / "cell_axis.tsv", 39)
    target_cells = cell_rows(target_input / "cell_axis.tsv", 12)
    source_context_features = feature_rows(paths["source"] / f"{context_role}.windows.tsv")
    target_context_features = feature_rows(target_input / "context.windows.tsv")
    source_target_features = feature_rows(paths["source"] / f"{target_role}.windows.tsv")
    if source_context_features != target_context_features:
        raise PeakVICrossCohortError("source and target context feature order differs")
    if {row["contig"] for row in source_context_features} & {
        row["contig"] for row in source_target_features
    }:
        raise PeakVICrossCohortError("context and target contigs overlap")
    source_context = sparse.load_npz(
        paths["source"] / f"{context_role}.counts.npz"
    ).tocsr()
    source_targets = sparse.load_npz(
        paths["source"] / f"{target_role}.counts.npz"
    ).tocsr()
    target_context = sparse.load_npz(target_input / "context.counts.npz").tocsr()
    if (
        source_context.shape != (len(source_cells), 16_000)
        or source_targets.shape != source_context.shape
        or target_context.shape != (len(target_cells), 16_000)
        or any((matrix.data < 0).any() for matrix in (source_context, source_targets, target_context))
        or any((matrix.getnnz(axis=1) == 0).any() for matrix in (source_context, source_targets, target_context))
    ):
        raise PeakVICrossCohortError("source or target sparse matrix differs")
    output.mkdir(parents=True, mode=0o750)
    seed_predictions: list[Any] = []
    seed_receipts: list[dict[str, Any]] = []
    dimensions = config["model"]["dimensions"]
    batch_size = int(config["training"]["batch_size"])
    for seed in config["seeds"]:
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        model = make_peakvae(16_000, dimensions)
        losses = fit_peakvae(
            torch,
            model,
            source_context,
            seed=seed,
            epochs=int(config["training"]["peakvi_epochs"]),
            batch_size=batch_size,
        )
        model_state = {name: value.detach().cpu() for name, value in model.state_dict().items()}
        model_hash = parameter_hash(model_state)
        seed_root = output / f"seed-{seed}"
        seed_root.mkdir(mode=0o750)
        model_state_path = seed_root / "peakvi.state_dict.pt"
        torch.save(model_state, model_state_path)
        restored_model = make_peakvae(16_000, dimensions)
        restored_model.load_state_dict(
            torch.load(model_state_path, map_location="cpu", weights_only=True),
            strict=True,
        )
        if parameter_hash(restored_model.state_dict()) != model_hash:
            raise PeakVICrossCohortError("strict PeakVI state restore differs")
        del model
        source_latent = encode(
            torch,
            restored_model,
            source_context,
            batch_size=batch_size,
            latent_dimension=int(dimensions["latent"]),
        )
        target_latent = encode(
            torch,
            restored_model,
            target_context,
            batch_size=batch_size,
            latent_dimension=int(dimensions["latent"]),
        )
        repeated = encode(
            torch,
            restored_model,
            target_context,
            batch_size=batch_size,
            latent_dimension=int(dimensions["latent"]),
        )
        if not np.array_equal(target_latent, repeated):
            raise PeakVICrossCohortError("target posterior-mean encoding is not repeatable")
        decoder, decoder_losses = fit_decoder(
            torch,
            source_latent,
            source_targets,
            seed=seed + 10_000,
            epochs=int(config["training"]["decoder_epochs"]),
            batch_size=batch_size,
        )
        decoder_state = {
            name: value.detach().cpu() for name, value in decoder.state_dict().items()
        }
        decoder_state_path = seed_root / "decoder.state_dict.pt"
        torch.save(decoder_state, decoder_state_path)
        restored_decoder = torch.nn.Linear(
            source_latent.shape[1], source_targets.shape[1]
        ).cuda()
        restored_decoder.load_state_dict(
            torch.load(decoder_state_path, map_location="cpu", weights_only=True),
            strict=True,
        )
        if parameter_hash(restored_decoder.state_dict()) != parameter_hash(decoder_state):
            raise PeakVICrossCohortError("strict decoder state restore differs")
        del decoder
        prediction = donor_lineage_predictions(
            torch,
            restored_decoder,
            target_latent,
            target_cells,
            batch_size=batch_size,
        )
        np.save(seed_root / "source_latent.float32.npy", source_latent, allow_pickle=False)
        np.save(seed_root / "target_latent.float32.npy", target_latent, allow_pickle=False)
        np.save(seed_root / "prediction.float32.npy", prediction, allow_pickle=False)
        seed_predictions.append(prediction)
        seed_receipts.append(
            {
                "seed": seed,
                "peakvi_parameter_sha256": model_hash,
                "decoder_parameter_sha256": parameter_hash(decoder_state),
                "peakvi_state_file_sha256": digest(model_state_path),
                "decoder_state_file_sha256": digest(decoder_state_path),
                "first_peakvi_loss_per_source_cell": losses[0],
                "last_peakvi_loss_per_source_cell": losses[-1],
                "first_decoder_bce": decoder_losses[0],
                "last_decoder_bce": decoder_losses[-1],
                "target_encoding_repeat_bit_identical": True,
                "strict_state_restore_passed": True,
            }
        )
        del restored_model, restored_decoder, model_state, decoder_state
        del source_latent, target_latent, repeated
        torch.cuda.empty_cache()
    ensemble = np.mean(np.stack(seed_predictions, axis=0), axis=0, dtype=np.float64).astype(np.float32)
    missing = np.zeros((12, 4, 16_000), dtype=np.uint8)
    predictions_path = output / "predictions.float32.npy"
    missing_path = output / "missing_state.uint8.npy"
    np.save(predictions_path, ensemble, allow_pickle=False)
    np.save(missing_path, missing, allow_pickle=False)
    write_json_exclusive(
        output / "native_representation_manifest.json",
        {
            "schema_version": "masld-bench-peakvi-native-representation-v1",
            "model_id": "peakvi",
            "output_family": "joint_representation",
            "rotation_id": rotation_id,
            "context_role": context_role,
            "source_cell_axis_sha256": digest(paths["source"] / "cell_axis.tsv"),
            "target_cell_axis_sha256": digest(target_input / "cell_axis.tsv"),
            "source_context_features_sha256": digest(
                paths["source"] / f"{context_role}.windows.tsv"
            ),
            "target_context_features_sha256": digest(target_input / "context.windows.tsv"),
            "seed_latent_paths": [
                {
                    "seed": seed,
                    "source": f"seed-{seed}/source_latent.float32.npy",
                    "target": f"seed-{seed}/target_latent.float32.npy",
                }
                for seed in config["seeds"]
            ],
            "direct_profile_ranking_allowed": False,
            "target_query_adaptation_performed": False,
            "target_role_atac_consumed": False,
        },
    )
    manifest = make_prediction_manifest(
        config=config,
        config_sha256=config_sha256,
        model_manifest=paths["model_manifest"],
        rotation_id=rotation_id,
        predictions=predictions_path,
        missing_state=missing_path,
        output=output,
    )
    write_json_exclusive(output / "prediction_bundle.json", manifest)
    write_json_exclusive(
        output / "fit_receipt.json",
        {
            "schema_version": "masld-bench-peakvi-cross-cohort-rotation-receipt-v1",
            "status": "passed",
            "rotation_id": rotation_id,
            "source_cells": len(source_cells),
            "target_cells": len(target_cells),
            "source_donors": 39,
            "target_donors": 12,
            "source_donor_lineage_units": 156,
            "target_donor_lineage_units": 48,
            "context_role": context_role,
            "scored_target_role": target_role,
            "seeds": seed_receipts,
            "native_output_family": "joint_representation",
            "projected_output_family": "masked_accessibility_count",
            "target_query_adaptation_performed": False,
            "condition_labels_read": False,
            "phenotype_values_read": False,
            "development_outcomes_read": False,
            "sealed_data_read": False,
            "metrics_calculated": False,
            "champion_claim_allowed": False,
        },
    )
    artifact_sha = freeze_tree(
        output,
        {
            "artifact_class": "peakvi_gse281367_blind_prediction",
            "model_id": "peakvi_source_fit_fixed_window_decoder",
            "rotation_id": rotation_id,
            "native_output_family": "joint_representation",
            "projected_output_family": "masked_accessibility_count",
            "development_outcomes_read": False,
            "metrics_calculated": False,
            "status": "passed",
        },
    )
    return {
        "rotation_id": rotation_id,
        "prediction_artifacts_sha256": artifact_sha,
        "seed_count": len(seed_receipts),
        "metrics_calculated": False,
    }


def run(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    import scvi
    import torch

    root = root.resolve(strict=True)
    if output.exists():
        raise PeakVICrossCohortError("refusing to overwrite fit-predict output")
    config, paths = validate_config(root, config_path)
    if scvi.__version__ != "1.5.0.post1" or torch.__version__ != "2.6.0+cu124":
        raise PeakVICrossCohortError("PeakVI runtime version differs")
    if torch.cuda.device_count() != 1 or "L40S" not in torch.cuda.get_device_name(0):
        raise PeakVICrossCohortError("exactly one L40S is required")
    torch.use_deterministic_algorithms(True)
    output.mkdir(parents=True, mode=0o750)
    config_sha256 = digest(config_path.resolve(strict=True))
    rotations = []
    for rotation_id in ROTATIONS:
        rotations.append(
            fit_rotation(
                torch=torch,
                config=config,
                config_sha256=config_sha256,
                paths=paths,
                rotation_id=rotation_id,
                output=output / rotation_id,
            )
        )
    receipt = {
        "schema_version": "masld-bench-peakvi-cross-cohort-fit-predict-receipt-v1",
        "status": "passed",
        "campaign_sha256": config_sha256,
        "rotations": rotations,
        "runtime": {
            "scvi_tools": scvi.__version__,
            "torch": torch.__version__,
            "cuda_build": torch.version.cuda,
            "gpu": torch.cuda.get_device_name(0),
        },
        "native_output_family": "joint_representation",
        "projected_output_family": "masked_accessibility_count",
        "condition_labels_read": False,
        "phenotype_values_read": False,
        "development_outcomes_read": False,
        "sealed_data_read": False,
        "metrics_calculated": False,
        "champion_claim_allowed": False,
    }
    write_json_exclusive(output / "receipt.json", receipt)
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    run(arguments.root, arguments.config, arguments.output)


if __name__ == "__main__":
    main()

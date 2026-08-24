#!/usr/bin/env python3
"""Fit PeakVI on training ATAC and decode training-lineage latent prototypes."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import multivi_fit_predict_smoke as common


MODEL_ID = "peakvi_training_context"
LINEAGES = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)


class PeakVIContextTransferError(ValueError):
    """Raised when leakage-safe PeakVI context transfer differs."""


def make_model(n_peaks: int) -> Any:
    from scvi.module import PEAKVAE

    return PEAKVAE(
        n_input_regions=n_peaks,
        n_batch=1,
        n_hidden=128,
        n_latent=16,
        n_layers_encoder=1,
        n_layers_decoder=1,
        dropout_rate=0.1,
        model_depth=True,
        region_factors=True,
        use_batch_norm="none",
        use_layer_norm="both",
        latent_distribution="normal",
    ).cuda()


def peak_tensors(torch: Any, atac: Any) -> dict[str, Any]:
    from scvi import REGISTRY_KEYS

    return {
        REGISTRY_KEYS.X_KEY: atac,
        REGISTRY_KEYS.BATCH_KEY: torch.zeros(
            (atac.shape[0], 1), dtype=torch.long, device=atac.device
        ),
    }


def train_model(
    torch: Any,
    model: Any,
    atac: Any,
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
        permutation = generator.permutation(atac.shape[0])
        epoch_loss = 0.0
        epoch_rows = 0
        model.train()
        for start in range(0, len(permutation), batch_size):
            positions = permutation[start : start + batch_size]
            values = common.dense_batch(torch, atac, positions, binary=True)
            optimizer.zero_grad(set_to_none=True)
            _inference, _generative, output = model(
                tensors=peak_tensors(torch, values), compute_loss=True
            )
            if not torch.isfinite(output.loss):
                raise PeakVIContextTransferError("PeakVI loss is non-finite")
            output.loss.backward()
            if not any(
                parameter.grad is not None and torch.isfinite(parameter.grad).all()
                for parameter in model.parameters()
            ):
                raise PeakVIContextTransferError("PeakVI finite gradients are absent")
            optimizer.step()
            count = len(positions)
            epoch_loss += float(output.loss.detach().cpu())
            epoch_rows += count
        losses.append(epoch_loss / epoch_rows)
    if not losses or any(not math.isfinite(value) for value in losses):
        raise PeakVIContextTransferError("PeakVI loss history differs")
    return losses


def training_latents(torch: Any, model: Any, atac: Any, batch_size: int) -> Any:
    import numpy as np

    results: list[Any] = []
    model.eval()
    with torch.inference_mode():
        for start in range(0, atac.shape[0], batch_size):
            positions = np.arange(start, min(start + batch_size, atac.shape[0]))
            values = common.dense_batch(torch, atac, positions, binary=True)
            inference, _generative = model(
                tensors=peak_tensors(torch, values),
                generative_kwargs={"use_z_mean": True},
                compute_loss=False,
            )
            means = inference["qz"].loc
            if means.shape != (len(positions), 16) or not torch.isfinite(means).all():
                raise PeakVIContextTransferError("PeakVI latent means differ")
            results.append(means.cpu().numpy())
    return np.vstack(results)


def decode_lineage_profiles(
    torch: Any,
    model: Any,
    latents: Any,
    training_rows: Sequence[Mapping[str, str]],
) -> Any:
    import numpy as np

    groups: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(training_rows):
        groups[row["lineage"]].append(index)
    if set(groups) != set(LINEAGES):
        raise PeakVIContextTransferError("training lineage roster differs")
    prototypes = np.vstack(
        [latents[groups[lineage]].mean(axis=0, dtype=np.float64) for lineage in LINEAGES]
    ).astype(np.float32)
    value = torch.from_numpy(prototypes).cuda()
    batch = torch.zeros((len(LINEAGES), 1), dtype=torch.long, device="cuda")
    model.eval()
    with torch.inference_mode():
        generated = model.generative(
            z=value,
            qz_m=value,
            batch_index=batch,
            use_z_mean=True,
        )["px"]
        if model.region_factors is not None:
            generated = generated * torch.sigmoid(model.region_factors)
    profiles = generated.cpu().numpy()
    if (
        profiles.shape != (len(LINEAGES), model.n_input_regions)
        or not np.all(np.isfinite(profiles))
        or np.any(profiles < 0)
    ):
        raise PeakVIContextTransferError("decoded PeakVI lineage profiles differ")
    return common.normalize_profiles(profiles)


def export_bundle(
    *,
    fold: int,
    run_id: str,
    query_rows: Sequence[Mapping[str, str]],
    peaks: Sequence[Mapping[str, str]],
    profiles: Any,
    output: Path,
    seed: int,
    training_nuclei: int,
    state_hash: str,
) -> Path:
    groups = sorted({(row["donor_id"], row["lineage"]) for row in query_rows})
    lineage_index = {lineage: index for index, lineage in enumerate(LINEAGES)}
    prediction_rows: list[dict[str, str]] = []
    row_id_rows: list[dict[str, str]] = []
    for donor, lineage in groups:
        donor_hash = common.join_hash(common.JOIN_NAMESPACE, "unit", donor)
        for peak_index, peak in enumerate(peaks):
            row_hash = common.join_hash(
                common.JOIN_NAMESPACE,
                "row",
                f"{donor}\0{lineage}\0{peak['peak_id']}",
            )
            prediction_rows.append(
                {
                    "row_hash": row_hash,
                    "donor_hash": donor_hash,
                    "block_hash": common.join_hash(
                        common.JOIN_NAMESPACE, "block", peak["chromosome"]
                    ),
                    "stratum": lineage,
                    "predicted": format(
                        float(profiles[lineage_index[lineage], peak_index]), ".17g"
                    ),
                }
            )
            row_id_rows.append({"row_hash": row_hash, "donor_hash": donor_hash})
    prediction_rows.sort(key=lambda row: row["row_hash"])
    row_id_rows.sort(key=lambda row: row["row_hash"])
    common.write_tsv(output / "predictions.tsv", common.PREDICTION_FIELDS, prediction_rows)
    common.write_tsv(output / "row_ids.tsv", ("row_hash", "donor_hash"), row_id_rows)
    prediction = common.artifact_record(
        output / "predictions.tsv",
        output,
        f"standardized_prediction_table:{common.TASK_ID}",
    )
    row_ids = common.artifact_record(
        output / "row_ids.tsv",
        output,
        f"prediction_row_ids:{common.TASK_ID}",
    )
    bundle = {
        "schema_version": "masld-bench-prediction-bundle-v1",
        "bundle_id": f"{MODEL_ID}-{run_id[:16]}",
        "run_id": run_id,
        "task_id": common.TASK_ID,
        "model_id": MODEL_ID,
        "dataset_ids": [common.DATASET_ID],
        "split_id": "donor_outer",
        "artifacts": [prediction, row_ids],
        "standardized_table": prediction,
        "row_ids": row_ids,
        "n_predictions": len(prediction_rows),
        "row_id_field": "row_hash",
        "unit_id_field": "donor_hash",
        "unit_id_namespace": common.JOIN_NAMESPACE,
        "biological_unit": "donor",
        "table_schema_sha256": common.canonical_hash(
            {"format": "tsv", "fields": list(common.PREDICTION_FIELDS)}
        ),
        "source_join_key_sha256": common.canonical_hash(
            {
                "task_id": common.TASK_ID,
                "dataset_ids": [common.DATASET_ID],
                "split_id": "donor_outer",
                "row_id_field": "row_hash",
                "unit_id_field": "donor_hash",
                "unit_id_namespace": common.JOIN_NAMESPACE,
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
            "query_stratum_count": len(groups),
            "selected_peak_count": len(peaks),
            "prediction_scale": "training_lineage_peakvi_latent_prototype_depth_free_multinomial",
            "donor_context": "none_training_lineage_prototype_repeated_across_held_donors",
            "query_rna_used": False,
            "query_atac_input_exposed": False,
            "held_atac_input_exposed": False,
            "observed_atac_exported": False,
            "held_atac_outcome_available_to_adapter": False,
            "state_dict_sha256": state_hash,
            "rna_conditioned_atac_eligible": False,
            "sealed_inference_eligible": False,
            "smoke_only": True,
            "champion_claim_allowed": False,
        },
    }
    path = output / "prediction_bundle.json"
    path.write_text(common.canonical_json(bundle) + "\n", encoding="utf-8")
    return path


def run(prepared: Path, output: Path, seed: int, epochs: int) -> dict[str, Any]:
    import numpy as np
    from scipy import sparse
    import scvi
    import torch

    if output.exists() or seed < 0 or epochs < 1:
        raise PeakVIContextTransferError("output, seed, or epochs differ")
    if scvi.__version__ != "1.5.0.post1" or torch.__version__ != "2.6.0+cu124":
        raise PeakVIContextTransferError("runtime version differs")
    if torch.cuda.device_count() != 1 or "L40S" not in torch.cuda.get_device_name(0):
        raise PeakVIContextTransferError("exactly one L40S is required")
    prepare_receipt = json.loads((prepared / "prepare_receipt.json").read_text())
    if (
        prepare_receipt.get("status") != "pass"
        or prepare_receipt.get("pairing_topology") != "same_nucleus"
        or prepare_receipt.get("held_atac_exported_from_prepare") is not False
    ):
        raise PeakVIContextTransferError("prepared contract differs")
    output.mkdir(parents=True, mode=0o750)
    torch.use_deterministic_algorithms(True)
    receipts: dict[str, Any] = {}
    for fold in range(5):
        fold_seed = seed + fold
        torch.manual_seed(fold_seed)
        torch.cuda.manual_seed_all(fold_seed)
        root = prepared / f"fold_{fold}"
        train_atac = sparse.load_npz(root / "training_atac.npz").tocsr()
        train_fields, train_rows = common.read_tsv(root / "training_rows.tsv")
        query_fields, query_rows = common.read_tsv(root / "query_rows.tsv")
        peak_fields, peaks = common.read_tsv(root / "selected_peaks.tsv")
        expected_rows = ("cell_id", "donor_id", "lineage", "rna_state", "atac_state")
        if (
            train_fields != expected_rows
            or query_fields != expected_rows
            or train_atac.shape != (len(train_rows), 10_000)
            or len(peaks) != 10_000
            or any(row["atac_state"] != "observed" for row in train_rows)
            or any(row["atac_state"] != "structurally_missing" for row in query_rows)
            or {row["donor_id"] for row in train_rows}
            & {row["donor_id"] for row in query_rows}
        ):
            raise PeakVIContextTransferError("fold input or donor firewall differs")
        common.validate_modality_states(
            np.asarray(train_atac.sum(axis=1)).ravel().tolist(),
            ["observed"] * len(train_rows),
        )
        model = make_model(train_atac.shape[1])
        losses = train_model(
            torch,
            model,
            train_atac,
            seed=fold_seed,
            epochs=epochs,
            batch_size=64,
        )
        state = {name: value.detach().cpu() for name, value in model.state_dict().items()}
        state_hash = common.state_dict_hash(state)
        fold_output = output / f"fold_{fold}"
        fold_output.mkdir(mode=0o750)
        state_path = fold_output / "state_dict.pt"
        torch.save(state, state_path)
        restored = make_model(train_atac.shape[1])
        restored.load_state_dict(
            torch.load(state_path, map_location="cpu", weights_only=True), strict=True
        )
        restored.cuda().eval()
        if common.state_dict_hash(restored.state_dict()) != state_hash:
            raise PeakVIContextTransferError("strict state resume differs")
        latents = training_latents(torch, restored, train_atac, batch_size=128)
        profiles = decode_lineage_profiles(torch, restored, latents, train_rows)
        repeat = decode_lineage_profiles(torch, restored, latents, train_rows)
        if not np.array_equal(profiles, repeat):
            raise PeakVIContextTransferError("repeated prototype decode differs")
        run_id = common.canonical_hash(
            {
                "model_id": MODEL_ID,
                "fold": fold,
                "seed": fold_seed,
                "epochs": epochs,
                "prepared_receipt_sha256": common.sha256_file(
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
            profiles=profiles,
            output=fold_output,
            seed=fold_seed,
            training_nuclei=len(train_rows),
            state_hash=state_hash,
        )
        receipts[str(fold)] = {
            "run_id": run_id,
            "training_nuclei": len(train_rows),
            "query_nuclei": len(query_rows),
            "first_epoch_loss_per_training_nucleus": losses[0],
            "last_epoch_loss_per_training_nucleus": losses[-1],
            "state_dict_sha256": state_hash,
            "state_dict_file_sha256": common.sha256_file(state_path),
            "prediction_bundle": bundle.relative_to(output).as_posix(),
            "prediction_bundle_sha256": common.sha256_file(bundle),
            "repeat_prototype_decode_bit_identical": True,
            "held_atac_read": False,
        }
        del model, restored, state, latents, profiles, repeat
        torch.cuda.empty_cache()
    receipt = {
        "schema_version": "masld-bench-peakvi-training-context-smoke-v1",
        "status": "pass",
        "model_id": MODEL_ID,
        "task_id": common.TASK_ID,
        "dataset_id": common.DATASET_ID,
        "pairing_topology": "same_nucleus_training",
        "outer_unit": "donor",
        "folds": receipts,
        "seed": seed,
        "epochs": epochs,
        "runtime": {
            "scvi_tools": scvi.__version__,
            "torch": torch.__version__,
            "cuda_build": torch.version.cuda,
            "gpu": torch.cuda.get_device_name(0),
        },
        "query_rna_used": False,
        "held_atac_available_to_adapter": False,
        "outcomes_read": False,
        "metrics_calculated_by_adapter": False,
        "rna_conditioned_atac_eligible": False,
        "sealed_inference_eligible": False,
        "smoke_only": True,
        "champion_claim_allowed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(common.canonical_json(receipt))
    return receipt


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

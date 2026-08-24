#!/usr/bin/env python3
"""Train a leakage-safe scButterfly-B core and predict held ATAC from RNA only."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence


MODEL_ID = "scbutterfly_b_safe"


class SCButterflyError(ValueError):
    """Raised when the scButterfly smoke execution differs."""


def load_component(path: Path) -> Any:
    spec = importlib.util.spec_from_file_location("scbutterfly_component_runtime", path)
    if spec is None or spec.loader is None:
        raise SCButterflyError("scButterfly component spec differs")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    required = {
        "NetBlock",
        "Split_Chrom_Encoder_block",
        "Split_Chrom_Decoder_block",
        "Translator",
    }
    if not required.issubset(vars(module)):
        raise SCButterflyError("scButterfly component class census differs")
    return module


def log_cpm(matrix: Any) -> Any:
    import numpy as np

    depth = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(depth <= 0):
        raise SCButterflyError("RNA depth differs")
    value = matrix.multiply((10_000.0 / depth)[:, None]).tocsr().astype(np.float32)
    value.data = np.log1p(value.data)
    return value


def training_tfidf(training: Any) -> tuple[Any, Any, float]:
    import numpy as np

    binary = training.copy().tocsr().astype(np.float32)
    binary.data[:] = 1.0
    depth = np.asarray(binary.sum(axis=1)).ravel()
    detection = np.asarray((binary > 0).sum(axis=0)).ravel()
    if np.any(depth <= 0):
        raise SCButterflyError("ATAC depth differs")
    inverse = np.log1p(binary.shape[0] / np.maximum(detection, 1)).astype(np.float32)
    value = binary.multiply((1.0 / depth)[:, None]).multiply(inverse).tocsr()
    scale = float(value.data.max())
    if not math.isfinite(scale) or scale <= 0:
        raise SCButterflyError("ATAC TF-IDF scale differs")
    value.data /= scale
    return value, inverse, scale


def chromosome_contract(peaks: Sequence[Mapping[str, str]]) -> list[int]:
    chromosomes = [row["chromosome"] for row in peaks]
    groups: list[tuple[str, int]] = []
    for chromosome in chromosomes:
        if not groups or groups[-1][0] != chromosome:
            groups.append((chromosome, 1))
        else:
            groups[-1] = (chromosome, groups[-1][1] + 1)
    if len({chromosome for chromosome, _count in groups}) != len(groups):
        raise SCButterflyError("peak chromosome blocks are not contiguous")
    counts = [count for _chromosome, count in groups]
    if len(counts) < 2 or sum(counts) != len(peaks):
        raise SCButterflyError("peak chromosome census differs")
    return counts


def make_model(component: Any, torch: Any, genes: int, peaks: int, chroms: list[int]) -> Any:
    nn = torch.nn

    class SafeButterfly(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            split_width = 32 * len(chroms)
            self.rna_encoder = component.NetBlock(
                2, [genes, 256, 128], [nn.LeakyReLU(), nn.LeakyReLU()], 0.1, 0.5
            )
            self.atac_encoder = component.Split_Chrom_Encoder_block(
                2,
                [peaks, split_width, 128],
                [nn.LeakyReLU(), nn.LeakyReLU()],
                chroms,
                0.1,
                0.3,
            )
            self.rna_decoder = component.NetBlock(
                2, [128, 256, genes], [nn.LeakyReLU(), nn.LeakyReLU()], 0.1, 0.0
            )
            self.atac_decoder = component.Split_Chrom_Decoder_block(
                2,
                [128, split_width, peaks],
                [nn.LeakyReLU(), nn.Sigmoid()],
                chroms,
                0.1,
                0.0,
            )
            self.translator = component.Translator(
                128,
                128,
                64,
                [nn.LeakyReLU(), nn.LeakyReLU(), nn.LeakyReLU()],
            )

        def paired(self, rna: Any, atac: Any, mode: str) -> tuple[Any, ...]:
            rna_encoded = self.rna_encoder(rna)
            atac_encoded = self.atac_encoder(atac)
            if mode == "train":
                rna_to_rna, rna_to_atac, rna_mu, rna_logvar = (
                    self.translator.train_model(rna_encoded, "RNA")
                )
                atac_to_rna, atac_to_atac, atac_mu, atac_logvar = (
                    self.translator.train_model(atac_encoded, "ATAC")
                )
            else:
                rna_to_rna, rna_to_atac, rna_mu, rna_logvar = (
                    self.translator.test_model(rna_encoded, "RNA")
                )
                atac_to_rna, atac_to_atac, atac_mu, atac_logvar = (
                    self.translator.test_model(atac_encoded, "ATAC")
                )
            return (
                self.rna_decoder(rna_to_rna),
                self.atac_decoder(rna_to_atac),
                self.rna_decoder(atac_to_rna),
                self.atac_decoder(atac_to_atac),
                rna_mu,
                rna_logvar,
                atac_mu,
                atac_logvar,
            )

        def rna_to_atac(self, rna: Any) -> Any:
            encoded = self.rna_encoder(rna)
            _rna, atac, _mu, _logvar = self.translator.test_model(encoded, "RNA")
            return self.atac_decoder(atac)

    return SafeButterfly()


def dense_rows(torch: Any, matrix: Any, indices: Any, device: str) -> Any:
    import numpy as np

    value = matrix[indices].toarray().astype(np.float32, copy=False)
    return torch.from_numpy(value).to(device)


def train_model(
    torch: Any,
    model: Any,
    training_rna: Any,
    training_atac: Any,
    *,
    seed: int,
    epochs: int,
    batch_size: int,
) -> list[float]:
    import numpy as np

    model.cuda().train()
    optimizer = torch.optim.Adam(model.parameters(), lr=1.0e-3)
    generator = torch.Generator(device="cpu")
    generator.manual_seed(seed)
    losses: list[float] = []
    for _epoch in range(epochs):
        order = torch.randperm(training_rna.shape[0], generator=generator).numpy()
        total = 0.0
        rows = 0
        for start in range(0, len(order), batch_size):
            positions = order[start : start + batch_size]
            if len(positions) < 2:
                continue
            rna = dense_rows(torch, training_rna, positions, "cuda")
            atac = dense_rows(torch, training_atac, positions, "cuda")
            optimizer.zero_grad(set_to_none=True)
            outputs = model.paired(rna, atac, "train")
            rna_rna, rna_atac, atac_rna, atac_atac = outputs[:4]
            rna_mu, rna_logvar, atac_mu, atac_logvar = outputs[4:]
            reconstruction = (
                torch.nn.functional.mse_loss(rna_rna, rna)
                + torch.nn.functional.mse_loss(atac_rna, rna)
                + torch.nn.functional.binary_cross_entropy(rna_atac, atac)
                + torch.nn.functional.binary_cross_entropy(atac_atac, atac)
            )
            kl = -0.5 * (
                torch.mean(1 + rna_logvar - rna_mu.square() - rna_logvar.exp())
                + torch.mean(1 + atac_logvar - atac_mu.square() - atac_logvar.exp())
            )
            loss = reconstruction + 1.0e-3 * kl
            if not torch.isfinite(loss):
                raise SCButterflyError("scButterfly training loss is non-finite")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            if not any(
                parameter.grad is not None and torch.isfinite(parameter.grad).all()
                for parameter in model.parameters()
            ):
                raise SCButterflyError("scButterfly finite gradients are absent")
            optimizer.step()
            total += float(loss.detach().cpu()) * len(positions)
            rows += len(positions)
        losses.append(total / rows)
    if len(losses) != epochs or any(not math.isfinite(value) for value in losses):
        raise SCButterflyError("scButterfly loss history differs")
    return losses


def predict(torch: Any, model: Any, query_rna: Any, batch_size: int) -> Any:
    import numpy as np

    model.eval()
    outputs: list[Any] = []
    with torch.inference_mode():
        for start in range(0, query_rna.shape[0], batch_size):
            positions = np.arange(start, min(start + batch_size, query_rna.shape[0]))
            rna = dense_rows(torch, query_rna, positions, "cuda")
            value = model.rna_to_atac(rna)
            if (
                value.shape[1] != 10_000
                or not torch.isfinite(value).all()
                or torch.any(value < 0)
                or torch.any(value > 1)
            ):
                raise SCButterflyError("scButterfly held profile differs")
            outputs.append(value.cpu().numpy())
    result = np.vstack(outputs)
    result += np.float32(1.0e-8)
    result /= result.sum(axis=1, keepdims=True)
    return result


def run(
    prepared: Path,
    component_path: Path,
    helper_path: Path,
    output: Path,
    seed: int,
    epochs: int,
    folds: Sequence[int] = tuple(range(5)),
) -> dict[str, Any]:
    import numpy as np
    from scipy import sparse
    import torch

    if (
        output.exists()
        or epochs < 1
        or not folds
        or len(set(folds)) != len(folds)
        or any(fold not in range(5) for fold in folds)
        or torch.cuda.device_count() != 1
    ):
        raise SCButterflyError("output, epoch, or GPU contract differs")
    if "L40S" not in torch.cuda.get_device_name(0):
        raise SCButterflyError("one L40S is required")
    component = load_component(component_path)
    spec = importlib.util.spec_from_file_location("multivi_export_helper", helper_path)
    if spec is None or spec.loader is None:
        raise SCButterflyError("prediction export helper spec differs")
    helper = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = helper
    spec.loader.exec_module(helper)
    helper.MODEL_ID = MODEL_ID
    prepare_receipt = json.loads((prepared / "prepare_receipt.json").read_text())
    if (
        prepare_receipt.get("status") != "pass"
        or prepare_receipt.get("pairing_topology") != "same_nucleus"
        or prepare_receipt.get("held_atac_exported_from_prepare") is not False
        or prepare_receipt.get("outcomes_read") is not False
    ):
        raise SCButterflyError("prepared contract differs")
    output.mkdir(parents=True, mode=0o750)
    torch.use_deterministic_algorithms(True)
    receipts: dict[str, Any] = {}
    for fold in sorted(folds):
        fold_seed = seed + fold
        torch.manual_seed(fold_seed)
        torch.cuda.manual_seed_all(fold_seed)
        root = prepared / f"fold_{fold}"
        train_rna_raw = sparse.load_npz(root / "training_rna.npz").tocsr()
        train_atac_raw = sparse.load_npz(root / "training_atac.npz").tocsr()
        query_rna_raw = sparse.load_npz(root / "query_rna.npz").tocsr()
        train_fields, train_rows = helper.read_tsv(root / "training_rows.tsv")
        query_fields, query_rows = helper.read_tsv(root / "query_rows.tsv")
        peak_fields, peaks = helper.read_tsv(root / "selected_peaks.tsv")
        if (
            train_fields
            != ("cell_id", "donor_id", "lineage", "rna_state", "atac_state")
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
            or train_atac_raw.shape != (train_rna_raw.shape[0], 10_000)
            or query_rna_raw.shape[1] != train_rna_raw.shape[1]
            or train_rna_raw.shape[1] != 2000
            or {row["donor_id"] for row in train_rows}
            & {row["donor_id"] for row in query_rows}
            or set(row["atac_state"] for row in query_rows) != {"structurally_missing"}
        ):
            raise SCButterflyError("fold axes, state, or donor firewall differs")
        chroms = chromosome_contract(peaks)
        train_rna = log_cpm(train_rna_raw)
        query_rna = log_cpm(query_rna_raw)
        train_atac, _inverse, _scale = training_tfidf(train_atac_raw)
        model = make_model(
            component, torch, train_rna.shape[1], train_atac.shape[1], chroms
        )
        losses = train_model(
            torch,
            model,
            train_rna,
            train_atac,
            seed=fold_seed,
            epochs=epochs,
            batch_size=64,
        )
        state = {
            name: value.detach().cpu().clone() for name, value in model.state_dict().items()
        }
        state_hash = helper.state_dict_hash(state)
        fold_output = output / f"fold_{fold}"
        fold_output.mkdir(mode=0o750)
        state_path = fold_output / "state_dict.pt"
        torch.save(state, state_path)
        restored = make_model(
            component, torch, train_rna.shape[1], train_atac.shape[1], chroms
        )
        restored.load_state_dict(
            torch.load(state_path, map_location="cpu", weights_only=True), strict=True
        )
        restored.cuda().eval()
        if helper.state_dict_hash(restored.state_dict()) != state_hash:
            raise SCButterflyError("strict scButterfly state resume differs")
        profiles = predict(torch, restored, query_rna, 128)
        repeat = predict(torch, restored, query_rna, 128)
        if not np.array_equal(profiles, repeat):
            raise SCButterflyError("repeated scButterfly query prediction differs")
        run_id = helper.canonical_hash(
            {
                "model_id": MODEL_ID,
                "fold": fold,
                "seed": fold_seed,
                "epochs": epochs,
                "state_dict_sha256": state_hash,
            }
        )
        bundle_path = helper.export_bundle(
            fold=fold,
            run_id=run_id,
            query_rows=query_rows,
            peaks=peaks,
            cell_profiles=profiles,
            output=fold_output,
            seed=fold_seed,
            training_nuclei=train_rna.shape[0],
            state_hash=state_hash,
        )
        bundle = json.loads(bundle_path.read_text())
        bundle["metadata"].update(
            {
                "adapter": "direct_low_level_scbutterfly_b_safe",
                "epochs": epochs,
                "loss_first": losses[0],
                "loss_last": losses[-1],
                "query_atac_tensor_created": False,
                "query_atac_tensor_nonzero_values": 0,
                "inverse_tfidf_executed": False,
                "high_level_combined_preprocessing_executed": False,
                "held_query_donor_or_lineage_used_by_model": False,
                "strict_state_resume": True,
                "repeated_query_prediction_bit_identical": True,
            }
        )
        bundle_path.write_text(helper.canonical_json(bundle) + "\n", encoding="utf-8")
        receipts[str(fold)] = {
            "training_nuclei": train_rna.shape[0],
            "query_nuclei": query_rna.shape[0],
            "genes": train_rna.shape[1],
            "peaks": train_atac.shape[1],
            "chromosome_blocks": len(chroms),
            "epochs": epochs,
            "loss_first": losses[0],
            "loss_last": losses[-1],
            "state_dict_sha256": state_hash,
            "strict_state_resume": True,
            "repeated_query_prediction_bit_identical": True,
            "held_atac_input_exposed": False,
            "inverse_tfidf_executed": False,
        }
    receipt = {
        "schema_version": "masld-bench-scbutterfly-b-safe-fit-v1",
        "status": "pass",
        "model_id": MODEL_ID,
        "exact_distribution": "scButterfly 0.0.9",
        "dataset_id": "gse296875",
        "pairing_topology": "same_nucleus",
        "outer_unit": "donor",
        "folds": receipts,
        "folds_requested": sorted(folds),
        "query_atac_state": "structurally_missing",
        "held_atac_input_exposed": False,
        "high_level_combined_preprocessing_executed": False,
        "inverse_tfidf_executed": False,
        "test_outcomes_read": False,
        "champion_claim_allowed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def self_test(component_path: Path) -> None:
    import numpy as np
    from scipy.sparse import csr_matrix

    component = load_component(component_path)
    value = log_cpm(csr_matrix(np.asarray([[1, 0], [0, 2]], dtype=np.int32)))
    if value.shape != (2, 2) or component.Translator.__name__ != "Translator":
        raise SCButterflyError("scButterfly self-test differs")
    print("status\tpass")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path)
    parser.add_argument("--component-path", type=Path, required=True)
    parser.add_argument("--helper-path", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--seed", type=int, default=1907)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--fold", type=int, action="append")
    parser.add_argument("--self-test", action="store_true")
    arguments = parser.parse_args()
    if arguments.self_test:
        self_test(arguments.component_path)
        return
    if arguments.prepared is None or arguments.helper_path is None or arguments.output is None:
        parser.error("--prepared, --helper-path, and --output are required")
    run(
        arguments.prepared,
        arguments.component_path,
        arguments.helper_path,
        arguments.output,
        arguments.seed,
        arguments.epochs,
        tuple(arguments.fold) if arguments.fold is not None else tuple(range(5)),
    )


if __name__ == "__main__":
    main()

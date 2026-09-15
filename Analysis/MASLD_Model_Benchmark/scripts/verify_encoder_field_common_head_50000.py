#!/usr/bin/env python3
"""Prove the encoder-field guards can fail before trusting anything they pass.

Every check below is run twice: once on the real inputs, and once on inputs
tampered on purpose.  A guard is reported ``demonstrably_live`` only when the
tamper actually changed the bytes it inspects AND the guard then rejected them.
A tamper that leaves the bytes identical, or that perturbs a quantity the guard
never reads, is reported as ``tamper_was_a_noop`` and the guard is NOT credited.

Two of these are exact-reproduction controls rather than digest comparisons:

  head_reproduction
      The head loop in this experiment must reproduce an already-frozen
      TranscriptFormer shard to within float noise.  If it cannot, "the same
      head on every encoder" is false and no delta in the panel is comparable.

  scvi_latent_placement
      ScviLatentBlock places reference_latent on the outer-training rows and
      query_latent on the held-out rows.  Nothing in the frozen bundle records
      that ordering, so it is proved by rebuilding the weighted-kNN baseline
      from those latents and reproducing the frozen scvi_baseline table.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any, Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for entry in (str(ROOT), str(ROOT / "src")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from masld_bench.hashing import sha256_file  # noqa: E402
from scripts.fit_predict_encoder_field_common_head_50000 import (  # noqa: E402
    EncoderFieldError,
    FrozenMatrixBlock,
    ScviLatentBlock,
    build_blocks,
    fit_shard,
    load_split,
)
from scripts.fit_predict_common_cell_heads_study_50000 import ROSTER  # noqa: E402


def _read_probabilities(path: Path) -> tuple[list[str], np.ndarray]:
    import csv

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    row_ids = [row["row_id"] for row in rows]
    values = np.asarray(
        [[float(row[f"probability::{label}"]) for label in ROSTER] for row in rows],
        dtype=np.float64,
    )
    return row_ids, values


def _tamper_record(
    *,
    name: str,
    guard: str,
    before_digest: str,
    after_digest: str,
    rejected: bool,
    detail: str,
) -> dict[str, Any]:
    bytes_changed = before_digest != after_digest
    if not bytes_changed:
        status = "tamper_was_a_noop"
    elif rejected:
        status = "demonstrably_live"
    else:
        status = "GUARD_DID_NOT_FIRE"
    return {
        "check": name,
        "guard": guard,
        "tamper_changed_bytes": bytes_changed,
        "digest_before": before_digest,
        "digest_after": after_digest,
        "guard_rejected_tampered_input": rejected,
        "status": status,
        "detail": detail,
    }


def _expect_raise(callable_: Callable[[], Any]) -> tuple[bool, str]:
    try:
        callable_()
    except EncoderFieldError as error:
        return True, str(error)
    except Exception as error:  # noqa: BLE001
        return True, f"{type(error).__name__}: {error}"
    return False, "no exception raised"


def check_head_reproduction(
    *, registry: dict[str, Any], context: dict[str, Any], source: Path
) -> list[dict[str, Any]]:
    spec = registry["reproduction_control"]
    frozen_path = ROOT / spec["frozen_shard_relative"]
    frozen_row_ids, frozen = _read_probabilities(frozen_path)
    blocks = build_blocks(
        registry=registry, context=context, fixture=source, only=[spec["block_id"]]
    )
    block = blocks[0]
    result = fit_shard(
        block=block,
        head_id=spec["head_id"],
        screen_seed=int(spec["screen_seed"]),
        held_outer=int(spec["outer_fold"]),
        context=context,
    )
    produced_row_ids = [context["row_ids"][index] for index in result["test_indices"]]
    produced = result["probabilities"]
    aligned = produced_row_ids == frozen_row_ids
    max_abs = float(np.max(np.abs(produced - frozen))) if aligned else float("inf")
    records = [
        {
            "check": "head_reproduction",
            "guard": "the experiment's head loop reproduces the frozen common-head shard",
            "frozen_shard": spec["frozen_shard_relative"],
            "frozen_shard_sha256": sha256_file(frozen_path),
            "rows": len(frozen_row_ids),
            "row_order_matches": aligned,
            "max_absolute_probability_difference": max_abs,
            "status": "verified_positively"
            if aligned and max_abs < 1.0e-9
            else "REPRODUCTION_FAILED",
        }
    ]

    # Tamper: zero one feature column the head demonstrably uses.  Zeroing a
    # column the model ignores would be a no-op, so the digest of the matrix is
    # compared before and after and the difference is asserted downstream.
    matrix = block.features(screen_seed=int(spec["screen_seed"]), held_outer=0)
    before = float(np.abs(matrix).sum())
    tampered = matrix.copy()
    variances = tampered.var(axis=0)
    loud_column = int(np.argmax(variances))
    tampered[:, loud_column] = 0.0
    after = float(np.abs(tampered).sum())

    class _TamperedBlock:
        block_id = block.block_id
        family = block.family
        width = block.width

        def identity(self) -> dict[str, Any]:
            return block.identity()

        def features(self, *, screen_seed: int, held_outer: int) -> np.ndarray:
            return tampered

    tampered_result = fit_shard(
        block=_TamperedBlock(),
        head_id=spec["head_id"],
        screen_seed=int(spec["screen_seed"]),
        held_outer=int(spec["outer_fold"]),
        context=context,
    )
    tampered_max_abs = float(
        np.max(np.abs(tampered_result["probabilities"] - frozen))
    )
    records.append(
        _tamper_record(
            name="head_reproduction_tamper",
            guard="zeroing the highest-variance feature column breaks the reproduction",
            before_digest=f"abs_sum={before:.6f}",
            after_digest=f"abs_sum={after:.6f}",
            rejected=tampered_max_abs >= 1.0e-9,
            detail=(
                f"zeroed column {loud_column} (variance {float(variances[loud_column]):.6g}); "
                f"max abs probability difference vs frozen shard {tampered_max_abs:.6g}"
            ),
        )
    )
    return records


def check_scvi_latent_placement(
    *, registry: dict[str, Any], context: dict[str, Any]
) -> list[dict[str, Any]]:
    from masld_bench.adapters.scvi_scanvi_inductive import (
        donor_class_weights as scvi_donor_class_weights,
        weighted_knn_probabilities,
    )

    spec = registry["scvi_latent_placement_control"]
    frozen_path = ROOT / spec["frozen_prediction_relative"]
    frozen_row_ids, frozen = _read_probabilities(frozen_path)
    outer = context["outer"]
    donors = context["donors"]
    targets = context["targets"]
    row_ids = context["row_ids"]
    bundle = ROOT / registry["scvi_bundle_relative"]
    seed = int(spec["seed"])
    neighbours = int(spec["knn_neighbors"])
    jobs = max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))

    # The frozen adapter is label-string native; targets are roster indices.
    labels = np.asarray([ROSTER[int(index)] for index in targets], dtype=str)
    produced = np.full((len(outer), len(ROSTER)), np.nan, dtype=np.float64)
    for fold in range(5):
        fold_root = bundle / f"seed-{seed}" / "folds" / f"fold{fold}"
        reference = np.asarray(np.load(fold_root / "reference_latent.npy"), np.float32)
        query = np.asarray(np.load(fold_root / "query_latent.npy"), np.float32)
        train = np.flatnonzero(outer != fold)
        test = np.flatnonzero(outer == fold)
        weights = scvi_donor_class_weights(
            donors[train].tolist(), labels[train].tolist()
        )
        produced[test] = weighted_knn_probabilities(
            reference,
            labels[train].tolist(),
            weights,
            query,
            neighbors=neighbours,
            n_jobs=jobs,
        )
    order = {row_id: index for index, row_id in enumerate(row_ids.tolist())}
    frozen_positions = [order[row_id] for row_id in frozen_row_ids]
    reordered = produced[frozen_positions]
    max_abs = float(np.max(np.abs(reordered - frozen)))
    records = [
        {
            "check": "scvi_latent_placement",
            "guard": (
                "reference_latent occupies the outer-training rows and query_latent "
                "the held-out rows, proved by rebuilding the frozen kNN baseline"
            ),
            "frozen_prediction": spec["frozen_prediction_relative"],
            "frozen_prediction_sha256": sha256_file(frozen_path),
            "seed": seed,
            "max_absolute_probability_difference": max_abs,
            "status": "verified_positively" if max_abs < 1.0e-9 else "PLACEMENT_FAILED",
        }
    ]

    # Tamper: reverse the reference rows within one fold.  The rows are the same
    # rows and the matrix norm is unchanged, so only the placement changes.
    fold = 0
    fold_root = bundle / f"seed-{seed}" / "folds" / f"fold{fold}"
    reference = np.asarray(np.load(fold_root / "reference_latent.npy"), np.float32)
    query = np.asarray(np.load(fold_root / "query_latent.npy"), np.float32)
    train = np.flatnonzero(outer != fold)
    test = np.flatnonzero(outer == fold)
    weights = scvi_donor_class_weights(donors[train].tolist(), labels[train].tolist())
    shuffled = reference[::-1].copy()
    tampered = weighted_knn_probabilities(
        shuffled,
        labels[train].tolist(),
        weights,
        query,
        neighbors=neighbours,
        n_jobs=jobs,
    )
    frozen_fold = {row_id: index for index, row_id in enumerate(frozen_row_ids)}
    fold_rows = [frozen_fold[row_id] for row_id in row_ids[test].tolist()]
    tampered_max_abs = float(np.max(np.abs(tampered - frozen[fold_rows])))
    records.append(
        _tamper_record(
            name="scvi_latent_placement_tamper",
            guard="reversing the reference-row order breaks the kNN reproduction",
            before_digest=f"first_row_sum={float(reference[0].sum()):.9f}",
            after_digest=f"first_row_sum={float(shuffled[0].sum()):.9f}",
            rejected=tampered_max_abs >= 1.0e-9,
            detail=(
                "reference rows reversed; identical row multiset, different "
                f"placement; max abs probability difference {tampered_max_abs:.6g}"
            ),
        )
    )
    return records


def check_row_order_guard(
    *, registry: dict[str, Any], context: dict[str, Any]
) -> list[dict[str, Any]]:
    entry = registry["frozen_matrix_blocks"][0]
    bundle_path = ROOT / entry["bundle_relative"]
    receipt_path = ROOT / entry["receipt_relative"]
    before = sha256_file(bundle_path)
    with tempfile.TemporaryDirectory() as scratch:
        target = Path(scratch) / "common_embeddings.npz"
        with np.load(bundle_path, allow_pickle=False) as source_bundle:
            arrays = {key: source_bundle[key] for key in source_bundle.files}
        swapped = arrays["row_ids"].copy()
        swapped[[0, 1]] = swapped[[1, 0]]
        arrays["row_ids"] = swapped
        with target.open("xb") as handle:
            np.savez(handle, **arrays)
        after = sha256_file(target)
        rejected, detail = _expect_raise(
            lambda: FrozenMatrixBlock(
                block_id=entry["block_id"],
                bundle_path=target,
                receipt_path=receipt_path,
                row_ids=context["row_ids"],
                exposure_status=entry["exposure_status"],
            )
        )
    return [
        _tamper_record(
            name="embedding_row_order_tamper",
            guard="an embedding whose row order differs from the fixture is refused",
            before_digest=before,
            after_digest=after,
            rejected=rejected,
            detail=f"swapped row_ids[0] and row_ids[1]; {detail}",
        )
    ]


def check_scvi_exposure_firewall(
    *, registry: dict[str, Any], context: dict[str, Any]
) -> list[dict[str, Any]]:
    """A held-out study appearing in the liver encoder's training roster must abort."""

    bundle = ROOT / registry["scvi_bundle_relative"]
    receipt_path = bundle / "seed-20260824" / "folds" / "fold0" / "fold_receipt.json"
    before = sha256_file(receipt_path)
    with tempfile.TemporaryDirectory() as scratch:
        shadow = Path(scratch) / "bundle"
        for seed in (20260824, 20260825, 20260826):
            for fold in range(5):
                source_fold = bundle / f"seed-{seed}" / "folds" / f"fold{fold}"
                target_fold = shadow / f"seed-{seed}" / "folds" / f"fold{fold}"
                target_fold.mkdir(parents=True)
                for name in (
                    "fold_receipt.json",
                    "reference_latent.npy",
                    "query_latent.npy",
                ):
                    os.symlink(source_fold / name, target_fold / name)
        tampered_receipt = shadow / "seed-20260824" / "folds" / "fold0" / "fold_receipt.json"
        payload = json.loads(tampered_receipt.read_text(encoding="utf-8"))
        tampered_receipt.unlink()
        # Put the held-out study into the training roster: the exact leak the
        # firewall exists to catch.
        payload["training_studies"] = sorted(
            set(payload["training_studies"]) | set(payload["query_studies"])
        )
        tampered_receipt.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        after = sha256_file(tampered_receipt)
        rejected, detail = _expect_raise(
            lambda: ScviLatentBlock(
                bundle_root=shadow,
                row_ids=context["row_ids"],
                outer=context["outer"],
                datasets=context["datasets"],
            )
        )
    return [
        _tamper_record(
            name="scvi_training_study_firewall_tamper",
            guard="the liver encoder must not list a held-out study among its training studies",
            before_digest=before,
            after_digest=after,
            rejected=rejected,
            detail=f"held-out study added to training_studies; {detail}",
        )
    ]


def check_prediction_tables_label_free(
    *, registry: dict[str, Any]
) -> list[dict[str, Any]]:
    import csv

    offenders: list[str] = []
    inspected = 0
    for entry in registry["prior_frozen_prediction_sets"]:
        path = ROOT / entry["predictions_relative"]
        with path.open("r", encoding="utf-8", newline="") as handle:
            fields = next(csv.reader(handle, delimiter="\t"))
        inspected += 1
        for field in fields:
            lowered = field.lower()
            if ("label" in lowered and not lowered.startswith("probability::")) or (
                "histolog" in lowered
            ):
                offenders.append(f"{entry['model_id']}:{field}")
    return [
        {
            "check": "prediction_tables_label_free",
            "guard": "no reused prediction table carries an observed label column",
            "tables_inspected": inspected,
            "offending_fields": offenders,
            "status": "verified_positively" if not offenders else "LABEL_LEAK",
        }
    ]


def main() -> int:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--source", required=True, type=Path)
    value.add_argument("--split", required=True, type=Path)
    value.add_argument("--registry", required=True, type=Path)
    value.add_argument("--output", required=True, type=Path)
    value.add_argument("--skip-gpu-checks", action="store_true")
    arguments = value.parse_args()
    registry = json.loads(arguments.registry.read_text(encoding="utf-8"))
    context = load_split(arguments.source, arguments.split)
    records: list[dict[str, Any]] = []
    records.extend(check_prediction_tables_label_free(registry=registry))
    records.extend(check_row_order_guard(registry=registry, context=context))
    records.extend(check_scvi_exposure_firewall(registry=registry, context=context))
    records.extend(check_scvi_latent_placement(registry=registry, context=context))
    if not arguments.skip_gpu_checks:
        records.extend(
            check_head_reproduction(
                registry=registry, context=context, source=arguments.source
            )
        )
    summary = {
        "schema_version": "masld-bench-encoder-field-guard-report-v1",
        "checks": records,
        "counts": {
            status: sum(1 for record in records if record["status"] == status)
            for status in sorted({record["status"] for record in records})
        },
        "all_guards_live_or_verified": all(
            record["status"] in {"demonstrably_live", "verified_positively"}
            for record in records
        ),
    }
    arguments.output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["all_guards_live_or_verified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

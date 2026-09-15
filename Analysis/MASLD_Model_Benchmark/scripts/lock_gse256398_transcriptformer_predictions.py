#!/usr/bin/env python3
"""Freeze GSE256398 prediction identities before donor metadata are joined."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


EXPECTED_ROWS = 172_997
EXPECTED_SOURCE_EXACT_ROWS = 165_372
EXPECTED_LEAVE_S35_OUT_ROWS = 163_801
POLICIES = ("native", "common")
ARMS = (
    "prospective_capped",
    "source_exact",
    "prospective_capped_leave_s35_out",
    "source_exact_leave_s35_out",
)


class GSE256398PredictionLockError(RuntimeError):
    """Raised when an embedding stream cannot be frozen without metadata."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def arm_mask(
    *, arm: str, donor: str, source_mask: np.ndarray, rows: int
) -> np.ndarray:
    if arm == "prospective_capped":
        return np.ones(rows, dtype=bool)
    if arm == "source_exact":
        return source_mask
    if arm == "prospective_capped_leave_s35_out":
        return np.full(rows, donor != "S35", dtype=bool)
    if arm == "source_exact_leave_s35_out":
        return source_mask & (donor != "S35")
    raise GSE256398PredictionLockError(f"unknown arm: {arm}")


def leave_s35_masks_identical(
    *, donor: str, source_mask: np.ndarray, rows: int
) -> bool:
    prospective = arm_mask(
        arm="prospective_capped_leave_s35_out",
        donor=donor,
        source_mask=source_mask,
        rows=rows,
    )
    source_exact = arm_mask(
        arm="source_exact_leave_s35_out",
        donor=donor,
        source_mask=source_mask,
        rows=rows,
    )
    return bool(np.array_equal(prospective, source_exact))


def update_prediction_digest(
    digest: Any, row_ids: Iterable[str], embeddings: np.ndarray
) -> int:
    count = 0
    for row_id, embedding in zip(row_ids, embeddings, strict=True):
        encoded = row_id.encode("utf-8")
        digest.update(len(encoded).to_bytes(4, "little"))
        digest.update(encoded)
        digest.update(np.asarray(embedding, dtype="<f4").tobytes(order="C"))
        count += 1
    return count


def build(
    *, embeddings_root: Path, embeddings_sha256: str, output: Path
) -> dict[str, Any]:
    if output.exists():
        raise GSE256398PredictionLockError("prediction lock output already exists")
    verify_frozen_tree(embeddings_root)
    if sha256_file(embeddings_root / "ARTIFACTS.json") != embeddings_sha256:
        raise GSE256398PredictionLockError("frozen embedding artifact differs")
    receipt = json.loads(
        (embeddings_root / "embeddings/embedding_receipt.json").read_text(
            encoding="utf-8"
        )
    )
    records = receipt.get("prediction_shards", [])
    if (
        receipt.get("status") != "pass_outcome_blind_cross_cohort_embeddings"
        or receipt.get("rows") != EXPECTED_ROWS
        or receipt.get("phenotype_metadata_read") is not False
        or receipt.get("sealed_outcomes_read") is not False
        or len(records) != 52
    ):
        raise GSE256398PredictionLockError("embedding firewall differs")

    digests = {
        (policy, arm): sha256(
            f"gse256398-transcriptformer-prediction-v1\0{policy}\0{arm}\0".encode(
                "utf-8"
            )
        )
        for policy in POLICIES
        for arm in ARMS
    }
    row_counts = {(policy, arm): 0 for policy in POLICIES for arm in ARMS}
    references: list[dict[str, Any]] = []
    for policy in POLICIES:
        policy_records = [record for record in records if record["policy"] == policy]
        if len(policy_records) != 26:
            raise GSE256398PredictionLockError("policy shard census differs")
        for record in policy_records:
            path = embeddings_root / "embeddings" / record["path"]
            if sha256_file(path) != record["sha256"]:
                raise GSE256398PredictionLockError("prediction shard hash differs")
            with np.load(path, allow_pickle=False) as values:
                if set(values.files) != {
                    "embeddings",
                    "in_source_exact",
                    "outer_folds",
                    "row_ids",
                }:
                    raise GSE256398PredictionLockError("prediction shard fields differ")
                embeddings = values["embeddings"]
                source_mask = values["in_source_exact"]
                row_ids = values["row_ids"].astype(str)
            if (
                embeddings.shape != (len(row_ids), 2048)
                or source_mask.shape != (len(row_ids),)
                or not np.all(np.isfinite(embeddings))
            ):
                raise GSE256398PredictionLockError("prediction shard values differ")
            if not leave_s35_masks_identical(
                donor=record["source_sample_id"],
                source_mask=source_mask,
                rows=len(row_ids),
            ):
                raise GSE256398PredictionLockError(
                    "leave-S35-out membership streams differ"
                )
            for arm in ARMS:
                mask = arm_mask(
                    arm=arm,
                    donor=record["source_sample_id"],
                    source_mask=source_mask,
                    rows=len(row_ids),
                )
                row_counts[(policy, arm)] += update_prediction_digest(
                    digests[(policy, arm)], row_ids[mask], embeddings[mask]
                )
            references.append(
                {
                    "policy": policy,
                    "source_sample_id": record["source_sample_id"],
                    "path": record["path"],
                    "sha256": record["sha256"],
                    "rows": record["rows"],
                }
            )

    expected_counts = {
        "prospective_capped": EXPECTED_ROWS,
        "source_exact": EXPECTED_SOURCE_EXACT_ROWS,
        "prospective_capped_leave_s35_out": EXPECTED_LEAVE_S35_OUT_ROWS,
        "source_exact_leave_s35_out": EXPECTED_LEAVE_S35_OUT_ROWS,
    }
    for policy in POLICIES:
        for arm, expected in expected_counts.items():
            if row_counts[(policy, arm)] != expected:
                raise GSE256398PredictionLockError("prediction arm row count differs")
    prediction_streams = {
        policy: {
            arm: {
                "rows": row_counts[(policy, arm)],
                "sha256": digests[(policy, arm)].hexdigest(),
            }
            for arm in ARMS
        }
        for policy in POLICIES
    }
    value = {
        "schema_version": "masld-bench-gse256398-transcriptformer-prediction-lock-v1",
        "status": "pass_predictions_locked_before_metadata",
        "dataset_id": "gse256398",
        "embedding_artifacts_sha256": embeddings_sha256,
        "prediction_stream_encoding": "namespace_then_repeated_uint32le_row_id_length_utf8_row_id_float32le_embedding",
        "prediction_streams": prediction_streams,
        "prediction_shard_references": references,
        "phenotype_metadata_read": False,
        "barcode_cell_labels_read": False,
        "sealed_outcomes_read": False,
        "supervised_accuracy_claim_allowed": False,
    }
    output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(output / "prediction_lock.json", value)
    freeze_tree(
        output,
        {
            "artifact_class": "gse256398_transcriptformer_prediction_lock",
            "dataset_id": "gse256398",
            "phenotype_metadata_read": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )
    verify_frozen_tree(output)
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings-root", type=Path, required=True)
    parser.add_argument("--embeddings-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    value = build(
        embeddings_root=args.embeddings_root,
        embeddings_sha256=args.embeddings_sha256,
        output=args.output,
    )
    print(json.dumps(value, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

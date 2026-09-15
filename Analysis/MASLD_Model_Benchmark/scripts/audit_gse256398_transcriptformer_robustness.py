#!/usr/bin/env python3
"""Describe GSE256398 TranscriptFormer QC robustness after fixing the predictions."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np

from masld_bench.artifacts import verify_frozen_tree, write_json_exclusive


EXPECTED_METADATA_SHA256 = (
    "e64fb9374ac96c6255b702084525c0b79b07cb4d449933e1f08870e960b652ba"
)
EXPECTED_ROWS = 172_997
EXPECTED_DONORS = {"masld_relevant": 17, "etiology_ood": 9}
POLICIES = ("native", "common")


class GSE256398RobustnessError(RuntimeError):
    """Raised when a frozen input or descriptive-only boundary differs."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def cosine_distance(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator == 0.0:
        raise GSE256398RobustnessError("zero-norm donor representation")
    return float(1.0 - np.dot(left, right) / denominator)


def summarize(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise GSE256398RobustnessError("summary values differ")
    return {
        "minimum": float(array.min()),
        "median": float(np.median(array)),
        "maximum": float(array.max()),
    }


def read_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 26:
        raise GSE256398RobustnessError("donor metadata census differs")
    allowed_roles = set(EXPECTED_DONORS)
    if set(row["development_role"] for row in rows) != allowed_roles:
        raise GSE256398RobustnessError("donor development roles differ")
    by_sample = {row["source_sample_id"]: row for row in rows}
    if len(by_sample) != 26:
        raise GSE256398RobustnessError("donor metadata key is not unique")
    for role, expected in EXPECTED_DONORS.items():
        if sum(row["development_role"] == role for row in rows) != expected:
            raise GSE256398RobustnessError("donor role cardinality differs")
    return by_sample


def donor_means(
    embeddings_root: Path, records: list[dict[str, Any]]
) -> dict[str, dict[str, dict[str, Any]]]:
    observed: dict[str, dict[str, dict[str, Any]]] = {policy: {} for policy in POLICIES}
    for record in records:
        policy = record["policy"]
        donor = record["source_sample_id"]
        if policy not in observed or donor in observed[policy]:
            raise GSE256398RobustnessError("prediction shard census differs")
        path = embeddings_root / "embeddings" / record["path"]
        if sha256_file(path) != record["sha256"]:
            raise GSE256398RobustnessError("prediction shard hash differs")
        with np.load(path, allow_pickle=False) as values:
            if set(values.files) != {
                "embeddings",
                "in_source_exact",
                "outer_folds",
                "row_ids",
            }:
                raise GSE256398RobustnessError("prediction shard fields differ")
            embeddings = values["embeddings"].astype(np.float64)
            source_mask = values["in_source_exact"].astype(bool)
        if (
            embeddings.ndim != 2
            or embeddings.shape[1] != 2048
            or source_mask.shape != (embeddings.shape[0],)
            or not source_mask.any()
            or not np.all(np.isfinite(embeddings))
        ):
            raise GSE256398RobustnessError("prediction shard values differ")
        observed[policy][donor] = {
            "prospective_rows": int(embeddings.shape[0]),
            "source_exact_rows": int(source_mask.sum()),
            "prospective": embeddings.mean(axis=0),
            "source_exact": embeddings[source_mask].mean(axis=0),
        }
    if any(len(observed[policy]) != 26 for policy in POLICIES):
        raise GSE256398RobustnessError("prediction donor census differs")
    return observed


def build(
    *,
    embeddings_root: Path,
    embeddings_sha256: str,
    prediction_lock_root: Path,
    prediction_lock_sha256: str,
    metadata_root: Path,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise GSE256398RobustnessError("audit output already exists")
    for root, expected in (
        (embeddings_root, embeddings_sha256),
        (prediction_lock_root, prediction_lock_sha256),
        (metadata_root, EXPECTED_METADATA_SHA256),
    ):
        verify_frozen_tree(root)
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise GSE256398RobustnessError("frozen input artifact differs")

    lock = json.loads(
        (prediction_lock_root / "prediction_lock.json").read_text(encoding="utf-8")
    )
    if (
        lock.get("status") != "pass_predictions_locked_before_metadata"
        or lock.get("embedding_artifacts_sha256") != embeddings_sha256
        or lock.get("phenotype_metadata_read") is not False
        or lock.get("sealed_outcomes_read") is not False
        or lock.get("supervised_accuracy_claim_allowed") is not False
    ):
        raise GSE256398RobustnessError("prediction lock boundary differs")
    receipt = json.loads(
        (embeddings_root / "embeddings/embedding_receipt.json").read_text(
            encoding="utf-8"
        )
    )
    records = receipt.get("prediction_shards", [])
    if (
        receipt.get("status") != "pass_outcome_blind_cross_cohort_embeddings"
        or receipt.get("rows") != EXPECTED_ROWS
        or len(records) != 52
        or receipt.get("phenotype_metadata_read") is not False
        or receipt.get("sealed_outcomes_read") is not False
    ):
        raise GSE256398RobustnessError("embedding firewall differs")

    metadata = read_metadata(metadata_root / "donor_metadata.tsv")
    representations = donor_means(embeddings_root, records)
    if set(representations[POLICIES[0]]) != set(metadata):
        raise GSE256398RobustnessError("prediction-to-metadata donor join differs")

    donor_rows: list[dict[str, Any]] = []
    distance_rows: list[dict[str, Any]] = []
    summaries: dict[str, Any] = {}
    for policy in POLICIES:
        policy_values = representations[policy]
        role_summaries: dict[str, Any] = {}
        for donor, values in sorted(policy_values.items()):
            source = values["source_exact"]
            prospective = values["prospective"]
            drift = float(np.linalg.norm(prospective - source) / np.linalg.norm(source))
            donor_rows.append(
                {
                    "policy": policy,
                    "source_sample_id": donor,
                    "development_role": metadata[donor]["development_role"],
                    "source_exact_rows": values["source_exact_rows"],
                    "prospective_rows": values["prospective_rows"],
                    "added_rows": values["prospective_rows"] - values["source_exact_rows"],
                    "source_vs_prospective_cosine_distance": cosine_distance(
                        source, prospective
                    ),
                    "source_vs_prospective_relative_l2_drift": drift,
                }
            )
        for role, expected in EXPECTED_DONORS.items():
            selected = [row for row in donor_rows if row["policy"] == policy and row["development_role"] == role]
            if len(selected) != expected:
                raise GSE256398RobustnessError("role-specific audit census differs")
            role_summaries[role] = {
                "donors": expected,
                "donors_with_membership_change": sum(row["added_rows"] > 0 for row in selected),
                "cosine_distance": summarize([row["source_vs_prospective_cosine_distance"] for row in selected]),
                "relative_l2_drift": summarize([row["source_vs_prospective_relative_l2_drift"] for row in selected]),
            }

        masld_donors = sorted(
            donor for donor, row in metadata.items() if row["development_role"] == "masld_relevant"
        )
        alcohol_donors = sorted(
            donor for donor, row in metadata.items() if row["development_role"] == "etiology_ood"
        )
        distance_summaries = {}
        for membership in ("source_exact", "prospective"):
            masld_vectors = np.stack(
                [policy_values[donor][membership] for donor in masld_donors]
            )
            masld_centroid = masld_vectors.mean(axis=0)
            for donor in alcohol_donors:
                vector = policy_values[donor][membership]
                distances = [
                    cosine_distance(vector, candidate) for candidate in masld_vectors
                ]
                distance_rows.append(
                    {
                        "policy": policy,
                        "membership": membership,
                        "source_sample_id": donor,
                        "development_role": "etiology_ood",
                        "distance_to_masld_centroid": cosine_distance(
                            vector, masld_centroid
                        ),
                        "nearest_masld_donor_distance": min(distances),
                    }
                )
            for index, donor in enumerate(masld_donors):
                vector = masld_vectors[index]
                other = np.delete(masld_vectors, index, axis=0)
                distance_rows.append(
                    {
                        "policy": policy,
                        "membership": membership,
                        "source_sample_id": donor,
                        "development_role": "masld_relevant_leave_one_out_context",
                        "distance_to_masld_centroid": cosine_distance(
                            vector, other.mean(axis=0)
                        ),
                        "nearest_masld_donor_distance": min(
                            cosine_distance(vector, candidate) for candidate in other
                        ),
                    }
                )
            alcohol_selected = [
                row
                for row in distance_rows
                if row["policy"] == policy
                and row["membership"] == membership
                and row["development_role"] == "etiology_ood"
            ]
            masld_selected = [
                row
                for row in distance_rows
                if row["policy"] == policy
                and row["membership"] == membership
                and row["development_role"]
                == "masld_relevant_leave_one_out_context"
            ]
            distance_summaries[membership] = {
                "alcohol_associated_donors": len(alcohol_selected),
                "masld_leave_one_out_donors": len(masld_selected),
                "alcohol_distance_to_masld_centroid": summarize(
                    [row["distance_to_masld_centroid"] for row in alcohol_selected]
                ),
                "masld_leave_one_out_distance_to_centroid": summarize(
                    [row["distance_to_masld_centroid"] for row in masld_selected]
                ),
                "alcohol_nearest_masld_distance": summarize(
                    [row["nearest_masld_donor_distance"] for row in alcohol_selected]
                ),
                "masld_leave_one_out_nearest_distance": summarize(
                    [row["nearest_masld_donor_distance"] for row in masld_selected]
                ),
            }
        summaries[policy] = {
            "membership_robustness": role_summaries,
            "descriptive_distance_context": distance_summaries,
        }

    output.mkdir(parents=True, exist_ok=False)
    for name, rows in (("donor_qc_robustness.tsv", donor_rows), ("etiology_ood_distances.tsv", distance_rows)):
        with (output / name).open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    value = {
        "schema_version": "masld-bench-gse256398-transcriptformer-descriptive-audit-v1",
        "status": "pass_descriptive_only_no_executable_task_spec",
        "dataset_id": "gse256398",
        "model_id": "transcriptformer_tf_sapiens",
        "unit_of_replication": "donor",
        "masld_relevant_donors": 17,
        "alcohol_associated_ood_donors": 9,
        "summaries": summaries,
        "prediction_lock_artifacts_sha256": prediction_lock_sha256,
        "embeddings_artifacts_sha256": embeddings_sha256,
        "metadata_artifacts_sha256": EXPECTED_METADATA_SHA256,
        "prediction_locked_before_metadata_join": True,
        "sealed_outcomes_read": False,
        "barcode_cell_labels_read": False,
        "supervised_accuracy_claim_allowed": False,
        "model_selection_allowed": False,
        "alcohol_donors_used_as_masld_negative_controls": False,
        "confirmatory_transfer_claim_allowed": False,
        "confirmatory_blocker": "No executable preregistered biological endpoint TaskSpec exists, and no authoritative barcode-level cell labels were released.",
        "interpretation_limit": "Distances are descriptive whole-donor embedding diagnostics confounded by cell composition and source disease context; they are not accuracy, diagnosis, or a biological replicate beyond donor.",
    }
    write_json_exclusive(output / "receipt.json", value)
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings-root", type=Path, required=True)
    parser.add_argument("--embeddings-sha256", required=True)
    parser.add_argument("--prediction-lock-root", type=Path, required=True)
    parser.add_argument("--prediction-lock-sha256", required=True)
    parser.add_argument("--metadata-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    value = build(
        embeddings_root=args.embeddings_root,
        embeddings_sha256=args.embeddings_sha256,
        prediction_lock_root=args.prediction_lock_root,
        prediction_lock_sha256=args.prediction_lock_sha256,
        metadata_root=args.metadata_root,
        output=args.output,
    )
    print(json.dumps(value, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build a training-donor-only locus-specific ATAC mean baseline."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from scripts.build_chrombpnet_development_outcomes import CCRE_FIELDS, SPLIT_FIELDS


MANIFEST_FIELDS = (
    "donor_test_fold",
    "donor_valid_fold",
    "donor_train_folds",
    "lineage_id",
    "donors",
    "nuclei",
    "unique_fragments",
    "read_support",
    "tn5_insertions",
    "nonzero_positions",
    "max_pending_positions",
    "fragment_path",
    "fragment_size_bytes",
    "fragment_sha256",
    "bigwig_path",
    "bigwig_size_bytes",
    "bigwig_sha256",
)


class TrainingMeanBaselineError(RuntimeError):
    """Raised when the training-only mean baseline contract is violated."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise TrainingMeanBaselineError(f"{path.name} fields differ")
        return [dict(row) for row in reader]


def select_training_row(
    rows: Sequence[Mapping[str, str]], *, split: Mapping[str, str], lineage: str
) -> dict[str, str]:
    matches = [
        dict(row)
        for row in rows
        if row["donor_test_fold"] == split["donor_test_fold"]
        and row["donor_valid_fold"] == split["donor_valid_fold"]
        and row["donor_train_folds"] == split["donor_train_folds"]
        and row["lineage_id"] == lineage
    ]
    if len(matches) != 1 or int(matches[0]["donors"]) < 2:
        raise TrainingMeanBaselineError("training pseudobulk row differs")
    return matches[0]


def role_windows(
    rows: Sequence[Mapping[str, str]], *, split: Mapping[str, str], role: str
) -> list[dict[str, str]]:
    fold_field = {"valid": "genomic_valid_fold", "test": "genomic_test_fold"}[role]
    fold = split[fold_field]
    selected = [dict(row) for row in rows if row["genomic_fold"] == fold]
    if len(selected) != 16_000 or any(
        int(row["output_end"]) - int(row["output_start"]) != 1_000
        for row in selected
    ):
        raise TrainingMeanBaselineError(f"{role} fixed-window contract differs")
    return selected


def _extract_role(arguments: Mapping[str, Any]) -> tuple[str, Any]:
    import numpy as np
    import pyBigWig

    path = Path(arguments["path"])
    rows = list(arguments["rows"])
    values = np.empty((len(rows), 1_000), dtype=np.uint32)
    bigwig = pyBigWig.open(str(path))
    try:
        for index, row in enumerate(rows):
            raw = np.asarray(
                bigwig.values(
                    row["contig"],
                    int(row["output_start"]),
                    int(row["output_end"]),
                    numpy=True,
                ),
                dtype=np.float64,
            )
            raw = np.nan_to_num(raw, nan=0.0, posinf=np.inf, neginf=-np.inf)
            rounded = np.rint(raw)
            if (
                raw.shape != (1_000,)
                or not np.isfinite(raw).all()
                or (raw < 0).any()
                or not np.allclose(raw, rounded, rtol=0, atol=1e-6)
                or (rounded > np.iinfo(np.uint32).max).any()
            ):
                raise TrainingMeanBaselineError("training track is not uint32-valued")
            values[index] = rounded.astype(np.uint32)
    finally:
        bigwig.close()
    return str(arguments["role"]), values


def build(
    *,
    pseudobulk: Path,
    pseudobulk_artifacts_sha256: str,
    split_contract: Path,
    split_artifacts_sha256: str,
    split_id: str,
    lineage: str,
    output: Path,
) -> dict[str, Any]:
    import h5py
    import numpy as np

    if output.exists():
        raise TrainingMeanBaselineError(f"output exists: {output}")
    pseudobulk = pseudobulk.resolve(strict=True)
    split_contract = split_contract.resolve(strict=True)
    if (
        sha256_file(pseudobulk / "ARTIFACTS.json")
        != pseudobulk_artifacts_sha256
        or sha256_file(split_contract / "ARTIFACTS.json")
        != split_artifacts_sha256
    ):
        raise TrainingMeanBaselineError("input ARTIFACTS SHA-256 differs")
    split_matches = [
        row
        for row in read_tsv(
            split_contract / "crossed_outer_splits.tsv", SPLIT_FIELDS
        )
        if row["split_id"] == split_id
    ]
    if len(split_matches) != 1:
        raise TrainingMeanBaselineError("split identifier differs")
    split = split_matches[0]
    training = select_training_row(
        read_tsv(pseudobulk / "training_pseudobulk_manifest.tsv", MANIFEST_FIELDS),
        split=split,
        lineage=lineage,
    )
    bigwig_path = (pseudobulk / training["bigwig_path"]).resolve(strict=True)
    if (
        bigwig_path.is_symlink()
        or bigwig_path.stat().st_size != int(training["bigwig_size_bytes"])
        or sha256_file(bigwig_path) != training["bigwig_sha256"]
    ):
        raise TrainingMeanBaselineError("training pseudobulk bigWig differs")
    windows = read_tsv(
        split_contract / "ccre_evaluation_windows.tsv", CCRE_FIELDS
    )
    by_role = {
        role: role_windows(windows, split=split, role=role)
        for role in ("valid", "test")
    }
    tasks = [
        {"role": role, "rows": rows, "path": bigwig_path}
        for role, rows in by_role.items()
    ]
    with ProcessPoolExecutor(max_workers=2) as executor:
        extracted = dict(executor.map(_extract_role, tasks))
    output.mkdir(mode=0o750)
    string_type = h5py.string_dtype(encoding="utf-8")
    with h5py.File(output / "training_mean_baseline.h5", "x") as handle:
        handle.attrs["schema_version"] = (
            "masld-bench-sequence-training-mean-baseline-v1"
        )
        handle.attrs["model_id"] = "training_pseudobulk_mean"
        handle.attrs["training_donors"] = int(training["donors"])
        handle.attrs["held_donor_outcomes_exposed"] = False
        for role, rows in by_role.items():
            counts = extracted[role]
            group = handle.create_group(role)
            group.create_dataset(
                "pseudobulk_counts",
                data=counts,
                dtype="u4",
                chunks=(min(32, len(rows)), 1_000),
                compression="gzip",
                compression_opts=1,
            )
            group.create_dataset(
                "mean_expected_counts",
                data=counts.sum(axis=1, dtype=np.float64) / int(training["donors"]),
                dtype="f8",
            )
            for field in ("window_id", "selection_hash", "contig"):
                group.create_dataset(
                    field,
                    data=np.asarray([row[field] for row in rows], dtype=object),
                    dtype=string_type,
                )
            group.attrs["genomic_fold"] = int(
                split[{"valid": "genomic_valid_fold", "test": "genomic_test_fold"}[role]]
            )
    summary = {
        "schema_version": "masld-bench-sequence-training-mean-baseline-v1",
        "status": "pass",
        "model_id": "training_pseudobulk_mean",
        "dataset_id": "gse296875",
        "split_id": split_id,
        "lineage_id": lineage,
        "training_donor_folds": split["donor_train_folds"],
        "training_donors": int(training["donors"]),
        "roles": {role: len(rows) for role, rows in by_role.items()},
        "source_signal": "training_donor_only_deduplicated_Tn5_pseudobulk",
        "held_donor_outcomes_exposed": False,
        "benchmark_metrics_calculated": False,
        "profile_smoothing_owned_by_evaluator": True,
        "pseudobulk_artifacts_sha256": pseudobulk_artifacts_sha256,
        "split_artifacts_sha256": split_artifacts_sha256,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pseudobulk", type=Path, required=True)
    parser.add_argument("--pseudobulk-artifacts-sha256", required=True)
    parser.add_argument("--split-contract", type=Path, required=True)
    parser.add_argument("--split-artifacts-sha256", required=True)
    parser.add_argument("--split-id", required=True)
    parser.add_argument("--lineage", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build(
        pseudobulk=args.pseudobulk,
        pseudobulk_artifacts_sha256=args.pseudobulk_artifacts_sha256,
        split_contract=args.split_contract,
        split_artifacts_sha256=args.split_artifacts_sha256,
        split_id=args.split_id,
        lineage=args.lineage,
        output=args.output,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

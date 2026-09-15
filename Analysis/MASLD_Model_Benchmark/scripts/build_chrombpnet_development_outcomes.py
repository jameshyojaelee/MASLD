#!/usr/bin/env python3
"""Build evaluator-only donor ATAC outcomes on fixed cCRE windows."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


BIGWIG_FIELDS = (
    "donor_id",
    "outer_fold",
    "lineage_id",
    "analysis_role",
    "nuclei",
    "unique_fragments",
    "read_support",
    "tn5_insertions",
    "nonzero_positions",
    "max_pending_positions",
    "path",
    "size_bytes",
    "sha256",
)
CCRE_FIELDS = (
    "contig",
    "output_start",
    "output_end",
    "window_id",
    "genomic_fold",
    "window_class",
    "ccre_class",
    "ccre_id",
    "ccre_start",
    "ccre_end",
    "input_start",
    "input_end",
    "selection_hash",
)
SPLIT_FIELDS = (
    "split_id",
    "donor_train_folds",
    "donor_valid_fold",
    "donor_test_fold",
    "genomic_train_folds",
    "genomic_valid_fold",
    "genomic_test_fold",
)


class DevelopmentOutcomeError(ValueError):
    """Raised when evaluator-only ATAC outcome extraction does not meet its requirements."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise DevelopmentOutcomeError(f"{path.name} fields differ")
        return [dict(row) for row in reader]


def _extract_track(arguments: Mapping[str, Any]) -> dict[str, object]:
    import numpy as np
    import pyBigWig

    path = Path(arguments["path"])
    if (
        path.is_symlink()
        or path.stat().st_size != int(arguments["size_bytes"])
        or sha256_file(path) != arguments["sha256"]
    ):
        raise DevelopmentOutcomeError("donor bigWig input differs")
    values = []
    bigwig = pyBigWig.open(str(path))
    try:
        for contig, start, end in arguments["windows"]:
            value = bigwig.stats(
                contig, int(start), int(end), type="sum", exact=True
            )[0]
            values.append(0.0 if value is None else float(value))
    finally:
        bigwig.close()
    counts = np.asarray(values, dtype=np.float64)
    rounded = np.rint(counts)
    if (
        not np.isfinite(counts).all()
        or (counts < 0).any()
        or not np.allclose(counts, rounded, rtol=0, atol=1e-6)
        or (rounded > np.iinfo(np.uint32).max).any()
    ):
        raise DevelopmentOutcomeError("donor cCRE counts are not uint32-valued")
    return {
        "donor_id": str(arguments["donor_id"]),
        "outer_fold": int(arguments["outer_fold"]),
        "nuclei": int(arguments["nuclei"]),
        "unique_fragments": int(arguments["unique_fragments"]),
        "values": rounded.astype(np.uint32),
    }


def build(
    *,
    donor_bigwigs: Path,
    donor_bigwigs_artifacts_sha256: str,
    split_contract: Path,
    split_artifacts_sha256: str,
    split_id: str,
    lineage: str,
    output: Path,
    workers: int,
    expected_windows_per_fold: int = 16_000,
) -> dict[str, object]:
    import h5py
    import numpy as np

    if output.exists() or not 1 <= workers <= 16 or expected_windows_per_fold < 1:
        raise DevelopmentOutcomeError("invalid output or worker contract")
    donor_bigwigs = donor_bigwigs.resolve(strict=True)
    split_contract = split_contract.resolve(strict=True)
    if (
        sha256_file(donor_bigwigs / "ARTIFACTS.json")
        != donor_bigwigs_artifacts_sha256
        or sha256_file(split_contract / "ARTIFACTS.json")
        != split_artifacts_sha256
    ):
        raise DevelopmentOutcomeError("input ARTIFACTS SHA-256 differs")
    split_matches = [
        row
        for row in read_tsv(
            split_contract / "crossed_outer_splits.tsv", SPLIT_FIELDS
        )
        if row["split_id"] == split_id
    ]
    if len(split_matches) != 1:
        raise DevelopmentOutcomeError("crossed split identifier differs")
    split = split_matches[0]
    bigwig_rows = [
        row
        for row in read_tsv(donor_bigwigs / "bigwig_manifest.tsv", BIGWIG_FIELDS)
        if row["lineage_id"] == lineage
    ]
    if not bigwig_rows:
        raise DevelopmentOutcomeError("lineage has no donor bigWigs")
    windows = read_tsv(
        split_contract / "ccre_evaluation_windows.tsv", CCRE_FIELDS
    )
    role_contract = {
        "valid": (
            int(split["donor_valid_fold"]),
            int(split["genomic_valid_fold"]),
        ),
        "test": (
            int(split["donor_test_fold"]),
            int(split["genomic_test_fold"]),
        ),
    }
    output.mkdir(mode=0o750)
    string_type = h5py.string_dtype(encoding="utf-8")
    role_summaries = {}
    with h5py.File(output / "outcomes.h5", "x") as handle:
        handle.attrs["schema_version"] = (
            "masld-bench-chrombpnet-development-outcomes-v1"
        )
        handle.attrs["dataset_id"] = "gse296875"
        handle.attrs["split_id"] = split_id
        handle.attrs["lineage_id"] = lineage
        handle.attrs["biological_unit"] = "donor"
        handle.attrs["model_inference_input_eligible"] = False
        for role, (donor_fold, genomic_fold) in role_contract.items():
            role_rows = sorted(
                [row for row in bigwig_rows if int(row["outer_fold"]) == donor_fold],
                key=lambda row: int(row["donor_id"]),
            )
            role_windows = [
                row for row in windows if int(row["genomic_fold"]) == genomic_fold
            ]
            if len(role_windows) != expected_windows_per_fold or not role_rows:
                raise DevelopmentOutcomeError("held donor or window census differs")
            extraction = [
                {
                    **row,
                    "path": str(donor_bigwigs / row["path"]),
                    "windows": [
                        (
                            window["contig"],
                            int(window["output_start"]),
                            int(window["output_end"]),
                        )
                        for window in role_windows
                    ],
                }
                for row in role_rows
            ]
            with ProcessPoolExecutor(max_workers=workers) as executor:
                results = list(executor.map(_extract_track, extraction))
            counts = np.vstack([row["values"] for row in results])
            group = handle.create_group(role)
            group.create_dataset(
                "counts",
                data=counts,
                dtype="u4",
                chunks=(1, min(4096, counts.shape[1])),
                compression="gzip",
                compression_opts=1,
            )
            group.create_dataset(
                "donor_id",
                data=np.asarray([row["donor_id"] for row in results], dtype=object),
                dtype=string_type,
            )
            group.create_dataset(
                "outer_fold",
                data=np.asarray([row["outer_fold"] for row in results], dtype=np.int8),
            )
            group.create_dataset(
                "nuclei",
                data=np.asarray([row["nuclei"] for row in results], dtype=np.int32),
            )
            for field in ("window_id", "selection_hash", "ccre_class"):
                group.create_dataset(
                    field,
                    data=np.asarray([row[field] for row in role_windows], dtype=object),
                    dtype=string_type,
                )
            group.attrs["donor_fold"] = donor_fold
            group.attrs["genomic_fold"] = genomic_fold
            group.attrs["model_inference_input_eligible"] = False
            role_summaries[role] = {
                "donor_fold": donor_fold,
                "genomic_fold": genomic_fold,
                "donors": len(results),
                "windows": len(role_windows),
                "total_tn5_insertions_in_windows": int(counts.sum()),
                "donors_with_positive_window_mass": int(
                    np.count_nonzero(counts.sum(axis=1) > 0)
                ),
            }
    if any(
        row["donors_with_positive_window_mass"] != row["donors"]
        for row in role_summaries.values()
    ):
        raise DevelopmentOutcomeError("held donor has no cCRE outcome mass")
    summary = {
        "schema_version": "masld-bench-chrombpnet-development-outcomes-v1",
        "status": "pass",
        "dataset_id": "gse296875",
        "split_id": split_id,
        "lineage_id": lineage,
        "biological_unit": "donor",
        "roles": role_summaries,
        "source_signal": "deduplicated_Cell_Ranger_ARC_Tn5_insertions",
        "observed_atac_outcomes": True,
        "model_inference_input_eligible": False,
        "benchmark_metrics_calculated": False,
        "donor_bigwigs_artifacts_sha256": donor_bigwigs_artifacts_sha256,
        "split_artifacts_sha256": split_artifacts_sha256,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--donor-bigwigs", type=Path, required=True)
    parser.add_argument("--donor-bigwigs-artifacts-sha256", required=True)
    parser.add_argument("--split-contract", type=Path, required=True)
    parser.add_argument("--split-artifacts-sha256", required=True)
    parser.add_argument("--split-id", required=True)
    parser.add_argument("--lineage", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    arguments = parser.parse_args()
    result = build(
        donor_bigwigs=arguments.donor_bigwigs,
        donor_bigwigs_artifacts_sha256=arguments.donor_bigwigs_artifacts_sha256,
        split_contract=arguments.split_contract,
        split_artifacts_sha256=arguments.split_artifacts_sha256,
        split_id=arguments.split_id,
        lineage=arguments.lineage,
        output=arguments.output,
        workers=arguments.workers,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

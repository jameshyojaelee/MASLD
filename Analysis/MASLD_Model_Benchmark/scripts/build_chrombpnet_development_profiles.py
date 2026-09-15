#!/usr/bin/env python3
"""Build evaluator-only base-resolution donor ATAC profiles on fixed cCREs."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from scripts.build_chrombpnet_development_outcomes import (
    BIGWIG_FIELDS,
    CCRE_FIELDS,
    SPLIT_FIELDS,
)


class DevelopmentProfileError(ValueError):
    """Raised when evaluator-only base-profile extraction does not meet its requirements."""


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
            raise DevelopmentProfileError(f"{path.name} fields differ")
        return [dict(row) for row in reader]


def _extract_track_profile(arguments: Mapping[str, Any]) -> dict[str, object]:
    import numpy as np
    import pyBigWig

    path = Path(arguments["path"])
    if (
        path.is_symlink()
        or path.stat().st_size != int(arguments["size_bytes"])
        or sha256_file(path) != arguments["sha256"]
    ):
        raise DevelopmentProfileError("donor bigWig input differs")
    windows = list(arguments["windows"])
    width = int(arguments["width"])
    values = np.empty((len(windows), width), dtype=np.uint32)
    bigwig = pyBigWig.open(str(path))
    try:
        for index, (contig, start, end) in enumerate(windows):
            raw = np.asarray(
                bigwig.values(contig, int(start), int(end), numpy=True),
                dtype=np.float64,
            )
            if raw.shape != (width,):
                raise DevelopmentProfileError("donor base-profile geometry differs")
            raw = np.nan_to_num(raw, nan=0.0, posinf=np.inf, neginf=-np.inf)
            rounded = np.rint(raw)
            if (
                not np.isfinite(raw).all()
                or (raw < 0).any()
                or not np.allclose(raw, rounded, rtol=0, atol=1e-6)
                or (rounded > np.iinfo(np.uint32).max).any()
            ):
                raise DevelopmentProfileError(
                    "donor base-profile values are not uint32-valued"
                )
            values[index] = rounded.astype(np.uint32)
    finally:
        bigwig.close()
    return {
        "donor_id": str(arguments["donor_id"]),
        "outer_fold": int(arguments["outer_fold"]),
        "nuclei": int(arguments["nuclei"]),
        "values": values,
    }


def _decode(values: Any) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def build(
    *,
    donor_bigwigs: Path,
    donor_bigwigs_artifacts_sha256: str,
    split_contract: Path,
    split_artifacts_sha256: str,
    regional_outcomes: Path,
    regional_outcomes_artifacts_sha256: str,
    split_id: str,
    lineage: str,
    output: Path,
    workers: int,
    expected_windows_per_fold: int = 16_000,
    expected_width: int = 1_000,
) -> dict[str, object]:
    import h5py
    import numpy as np

    if (
        output.exists()
        or not 1 <= workers <= 16
        or expected_windows_per_fold < 1
        or expected_width < 1
    ):
        raise DevelopmentProfileError("invalid output or worker contract")
    donor_bigwigs = donor_bigwigs.resolve(strict=True)
    split_contract = split_contract.resolve(strict=True)
    regional_outcomes = regional_outcomes.resolve(strict=True)
    expected_artifacts = (
        (donor_bigwigs, donor_bigwigs_artifacts_sha256),
        (split_contract, split_artifacts_sha256),
        (regional_outcomes, regional_outcomes_artifacts_sha256),
    )
    if any(
        sha256_file(root / "ARTIFACTS.json") != expected
        for root, expected in expected_artifacts
    ):
        raise DevelopmentProfileError("input ARTIFACTS SHA-256 differs")
    split_matches = [
        row
        for row in read_tsv(
            split_contract / "crossed_outer_splits.tsv", SPLIT_FIELDS
        )
        if row["split_id"] == split_id
    ]
    if len(split_matches) != 1:
        raise DevelopmentProfileError("crossed split identifier differs")
    split = split_matches[0]
    bigwig_rows = [
        row
        for row in read_tsv(donor_bigwigs / "bigwig_manifest.tsv", BIGWIG_FIELDS)
        if row["lineage_id"] == lineage
    ]
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
    with h5py.File(
        regional_outcomes / "outcomes/outcomes.h5", "r"
    ) as regional, h5py.File(output / "profiles.h5", "x") as handle:
        if (
            regional.attrs.get("split_id") != split_id
            or regional.attrs.get("lineage_id") != lineage
            or bool(regional.attrs.get("model_inference_input_eligible"))
        ):
            raise DevelopmentProfileError("regional outcome binding differs")
        handle.attrs["schema_version"] = (
            "masld-bench-chrombpnet-development-base-profiles-v1"
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
                raise DevelopmentProfileError("held donor or window census differs")
            if any(
                int(row["output_end"]) - int(row["output_start"]) != expected_width
                for row in role_windows
            ):
                raise DevelopmentProfileError("fixed cCRE output width differs")
            extraction = [
                {
                    **row,
                    "path": str(donor_bigwigs / row["path"]),
                    "width": expected_width,
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
                results = list(executor.map(_extract_track_profile, extraction))
            profiles = np.stack([row["values"] for row in results], axis=0)
            donor_ids = [row["donor_id"] for row in results]
            window_ids = [row["window_id"] for row in role_windows]
            regional_group = regional[role]
            if (
                donor_ids != _decode(regional_group["donor_id"][:])
                or window_ids != _decode(regional_group["window_id"][:])
                or not np.array_equal(
                    profiles.sum(axis=2, dtype=np.uint64),
                    np.asarray(regional_group["counts"][:], dtype=np.uint64),
                )
            ):
                raise DevelopmentProfileError(
                    "base profiles do not reproduce frozen regional outcomes"
                )
            group = handle.create_group(role)
            group.create_dataset(
                "counts",
                data=profiles,
                dtype="u4",
                chunks=(1, min(32, profiles.shape[1]), expected_width),
                compression="gzip",
                compression_opts=1,
            )
            group.create_dataset(
                "donor_id",
                data=np.asarray(donor_ids, dtype=object),
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
                "positions_per_window": expected_width,
                "total_tn5_insertions_in_profiles": int(
                    profiles.sum(dtype=np.uint64)
                ),
                "maximum_position_count": int(profiles.max()),
                "regional_outcomes_reproduced": True,
            }
            del profiles
    summary = {
        "schema_version": "masld-bench-chrombpnet-development-base-profiles-v1",
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
        "regional_outcomes_reproduced": True,
        "donor_bigwigs_artifacts_sha256": donor_bigwigs_artifacts_sha256,
        "split_artifacts_sha256": split_artifacts_sha256,
        "regional_outcomes_artifacts_sha256": regional_outcomes_artifacts_sha256,
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
    parser.add_argument("--regional-outcomes", type=Path, required=True)
    parser.add_argument("--regional-outcomes-artifacts-sha256", required=True)
    parser.add_argument("--split-id", required=True)
    parser.add_argument("--lineage", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    arguments = parser.parse_args()
    result = build(
        donor_bigwigs=arguments.donor_bigwigs,
        donor_bigwigs_artifacts_sha256=arguments.donor_bigwigs_artifacts_sha256,
        split_contract=arguments.split_contract,
        split_artifacts_sha256=arguments.split_artifacts_sha256,
        regional_outcomes=arguments.regional_outcomes,
        regional_outcomes_artifacts_sha256=(
            arguments.regional_outcomes_artifacts_sha256
        ),
        split_id=arguments.split_id,
        lineage=arguments.lineage,
        output=arguments.output,
        workers=arguments.workers,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

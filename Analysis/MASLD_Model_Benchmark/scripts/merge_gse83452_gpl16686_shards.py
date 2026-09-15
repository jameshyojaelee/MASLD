#!/usr/bin/env python3
"""Merge GSE83452 shard matrices, verifying cross-job reproducibility first.

Sharding is only legitimate if an array's summary does not depend on which job
computed it.  The sentinel arrays appear in every shard, so this module can
check that directly: independent jobs, on different nodes, must produce
bit-identical digests for them.  If they do not, the shards are discarded rather
than merged, because the whole single-array guarantee would be void.

Each array's column is also re-digested from the shard file it was read out of
and required to match the digest the shard recorded, so a merge cannot silently
transpose or drop a column.

Every owned array must appear in exactly one shard, and the merged axis must
match the frozen crosswalk exactly.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import struct
from typing import Any, Mapping

SERIES = "GSE83452"
PLATFORM_ID = "GPL16686"
EXPECTED_ARRAYS = 231
EXPECTED_PLATFORM_FEATURES = 53_617
EXPECTED_PACKAGE_VERSIONS = {
    "oligo": "1.74.0",
    "affyio": "1.80.0",
    "affxparser": "1.82.0",
    "pd.hugene.2.0.st": "3.14.1",
}


class ShardMergeError(RuntimeError):
    """Raised when the shards cannot be merged into a trustworthy matrix."""


def _scalar(value: Any) -> Any:
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


def digest_column(values: list[float]) -> str:
    """Digest a column exactly as the R stage does: raw little-endian doubles."""
    payload = struct.pack(f"<{len(values)}d", *values)
    return sha256(payload).hexdigest()


def read_matrix(path: Path, id_column: str) -> tuple[list[str], list[str], list[list[float]]]:
    """Return (row ids, column ids, columns) with columns as float lists."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[0] != id_column:
            raise ShardMergeError(f"{path.name} first column is not {id_column}")
        columns = header[1:]
        rows: list[str] = []
        data: list[list[float]] = [[] for _ in columns]
        for record in reader:
            if len(record) != len(header):
                raise ShardMergeError(f"{path.name} is ragged at row {len(rows) + 1}")
            rows.append(record[0])
            for index, cell in enumerate(record[1:]):
                data[index].append(float(cell))
    if len(set(rows)) != len(rows):
        raise ShardMergeError(f"{path.name} repeats a feature identifier")
    if len(set(columns)) != len(columns):
        raise ShardMergeError(f"{path.name} repeats a sample identifier")
    return rows, columns, data


def verify_shard(receipt: Mapping[str, Any]) -> None:
    if _scalar(receipt.get("status")) != "pass_label_blind_single_array_shard":
        raise ShardMergeError("a shard receipt is not a pass")
    if _scalar(receipt.get("series")) != SERIES or _scalar(receipt.get("platform_id")) != PLATFORM_ID:
        raise ShardMergeError("a shard receipt names a different source")
    if _scalar(receipt.get("arrays_read_per_summarization_call")) != 1:
        raise ShardMergeError("a shard read more than one array per call")
    if _scalar(receipt.get("matrices_round_trip_bitwise")) is not True:
        raise ShardMergeError("a shard did not assert a bitwise round trip")
    if _scalar(receipt.get("label_blind")) is not True:
        raise ShardMergeError("a shard is not marked label blind")
    if _scalar(receipt.get("samples_excluded")) != 0:
        raise ShardMergeError("a shard excluded a sample")
    for flag in (
        "all_sample_RMA_run",
        "across_array_quantile_normalization_run",
        "GEO_series_matrix_read",
        "labels_read",
        "model_training_activated",
        "sealed_outcomes_read",
    ):
        if _scalar(receipt.get(flag)) is not False:
            raise ShardMergeError(f"a shard did not clear {flag}")
    versions = {k: _scalar(v) for k, v in (receipt.get("package_versions") or {}).items()}
    if versions != EXPECTED_PACKAGE_VERSIONS:
        raise ShardMergeError(f"shard package versions differ: {versions!r}")


def merge(
    *, shard_dirs: list[Path], plan: Mapping[str, Any], output: Path
) -> dict[str, Any]:
    sentinels = set(plan["sentinels"])
    if not sentinels:
        raise ShardMergeError("shard plan carries no sentinel arrays")

    columns_by_accession: dict[str, list[float]] = {}
    owner_of: dict[str, int] = {}
    sentinel_digests: dict[str, dict[int, str]] = {name: {} for name in sentinels}
    platform_axis: list[str] | None = None
    gene_axis: list[str] | None = None
    gene_by_accession: dict[str, list[float]] = {}
    qc_rows: list[dict[str, str]] = []

    for shard_dir in sorted(shard_dirs):
        receipts = sorted(shard_dir.glob("shard_*_receipt.json"))
        if len(receipts) != 1:
            raise ShardMergeError(f"{shard_dir} does not carry exactly one shard receipt")
        receipt = json.loads(receipts[0].read_text(encoding="utf-8"))
        verify_shard(receipt)
        shard_index = int(_scalar(receipt["shard"]))

        rows, cols, data = read_matrix(
            shard_dir / f"shard_{shard_index:02d}_platform_feature_matrix.tsv",
            "platform_feature_id",
        )
        if platform_axis is None:
            platform_axis = rows
        elif rows != platform_axis:
            raise ShardMergeError("shards disagree on the platform feature axis")

        grows, gcols, gdata = read_matrix(
            shard_dir / f"shard_{shard_index:02d}_gene_matrix.tsv", "ensembl_gene_id"
        )
        if gene_axis is None:
            gene_axis = grows
        elif grows != gene_axis:
            raise ShardMergeError("shards disagree on the gene axis")
        if gcols != cols:
            raise ShardMergeError("a shard's gene and platform matrices disagree on samples")

        recorded = {
            _scalar(record["sample_accession"]): _scalar(record["array_digest"])
            for record in receipt["per_array"]
        }
        for position, accession in enumerate(cols):
            column = data[position]
            observed = digest_column(column)
            if observed != recorded.get(accession):
                raise ShardMergeError(
                    f"{accession} column digest does not match the shard receipt"
                )
            if accession in sentinels:
                sentinel_digests[accession][shard_index] = observed
                continue
            if accession in owner_of:
                raise ShardMergeError(f"{accession} appears in more than one shard as an owner")
            owner_of[accession] = shard_index
            columns_by_accession[accession] = column
            gene_by_accession[accession] = gdata[position]

        for record in receipt["per_array"]:
            qc_rows.append({k: str(_scalar(v)) for k, v in record.items()})

    # Cross-job reproducibility: independent jobs must agree bitwise.
    for accession, digests in sentinel_digests.items():
        if len(digests) < 2:
            raise ShardMergeError(f"sentinel {accession} was not summarized in at least two shards")
        if len(set(digests.values())) != 1:
            raise ShardMergeError(
                f"sentinel {accession} differs between shards; single-array isolation does not "
                f"hold across jobs: {digests!r}"
            )

    # A sentinel is also an owned array in exactly one shard; recover it there.
    for accession in sentinels:
        if accession in columns_by_accession:
            continue
        for shard_dir in sorted(shard_dirs):
            receipts = sorted(shard_dir.glob("shard_*_receipt.json"))
            receipt = json.loads(receipts[0].read_text(encoding="utf-8"))
            shard_index = int(_scalar(receipt["shard"]))
            rows, cols, data = read_matrix(
                shard_dir / f"shard_{shard_index:02d}_platform_feature_matrix.tsv",
                "platform_feature_id",
            )
            _, gcols, gdata = read_matrix(
                shard_dir / f"shard_{shard_index:02d}_gene_matrix.tsv", "ensembl_gene_id"
            )
            if accession in cols:
                position = cols.index(accession)
                columns_by_accession[accession] = data[position]
                gene_by_accession[accession] = gdata[gcols.index(accession)]
                owner_of[accession] = shard_index
                break

    accessions = sorted(columns_by_accession)
    if len(accessions) != EXPECTED_ARRAYS:
        raise ShardMergeError(f"merged cohort is {len(accessions)} arrays, not {EXPECTED_ARRAYS}")
    if platform_axis is None or len(platform_axis) != EXPECTED_PLATFORM_FEATURES:
        raise ShardMergeError("merged platform axis differs from the frozen crosswalk")
    if gene_axis is None or len(gene_axis) < 1:
        raise ShardMergeError("merged gene axis is empty")

    output.mkdir(parents=True, exist_ok=False)

    def write(path: Path, id_column: str, axis: list[str], source: dict[str, list[float]]) -> str:
        with path.open("w", encoding="utf-8", newline="") as handle:
            handle.write("\t".join([id_column, *accessions]) + "\n")
            for index, feature in enumerate(axis):
                cells = ["%.17g" % source[a][index] for a in accessions]
                handle.write("\t".join([feature, *cells]) + "\n")
        # Re-read and require bitwise equality, same guard as the R stage.
        rows, cols, data = read_matrix(path, id_column)
        if rows != axis or cols != accessions:
            raise ShardMergeError(f"{path.name} axis changed on write")
        for position, accession in enumerate(accessions):
            if data[position] != source[accession]:
                raise ShardMergeError(f"{path.name} does not round-trip bitwise for {accession}")
        digest = sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    platform_sha = write(
        output / "gse83452_gpl16686_platform_feature_matrix.tsv",
        "platform_feature_id", platform_axis, columns_by_accession,
    )
    gene_sha = write(
        output / "gse83452_gpl16686_gene_matrix.tsv",
        "ensembl_gene_id", gene_axis, gene_by_accession,
    )
    (output / "summarized_columns.txt").write_text("\n".join(accessions) + "\n", encoding="utf-8")

    if qc_rows:
        fieldnames = list(qc_rows[0])
        with (output / "gse83452_qc_metrics.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(qc_rows)

    receipt = {
        "schema_version": "masld-bench-gse83452-gpl16686-merged-summary-v1",
        "status": "pass_label_blind_single_array_summarization",
        "series": SERIES,
        "platform_id": PLATFORM_ID,
        "cohort_family_id": "antwerp_inserm_shared",
        "accession_aliases": ["GSE106737", "GSE83452"],
        "independent_evaluation_between_accession_aliases": False,
        "method": "single_array_core_summary_normalize_false",
        "arrays_summarized": len(accessions),
        "arrays_read_per_summarization_call": 1,
        "platform_features": len(platform_axis),
        "gene_features": len(gene_axis),
        "shards": len(shard_dirs),
        "sentinels": sorted(sentinels),
        "sentinel_digests_agree_across_shards": True,
        "cross_job_single_array_reproducibility_verified": True,
        "every_column_digest_matched_its_shard_receipt": True,
        "platform_feature_matrix_sha256": platform_sha,
        "gene_matrix_sha256": gene_sha,
        "matrices_round_trip_bitwise": True,
        "samples_excluded": 0,
        "label_blind": True,
        "all_sample_RMA_run": False,
        "across_array_quantile_normalization_run": False,
        "GEO_series_matrix_read": False,
        "labels_read": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
    }
    (output / "merged_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard-dir", action="append", dest="shard_dirs", required=True, type=Path)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    receipt = merge(
        shard_dirs=arguments.shard_dirs,
        plan=json.loads(arguments.plan.read_text(encoding="utf-8")),
        output=arguments.output,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

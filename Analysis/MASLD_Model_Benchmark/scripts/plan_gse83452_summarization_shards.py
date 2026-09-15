#!/usr/bin/env python3
"""Split the GSE83452 arrays into independent, label-blind summarization shards.

The 231 arrays are independent by construction once the separation has shown that
each array's summary does not depend on what else the run processed, so they can
be summarized in parallel jobs.  The split is purely positional over sorted GSM
accessions: no label, phenotype, intensity, timepoint, or participant grouping
takes part, so the shard an array lands in carries no information.

Every shard also carries the same sentinel arrays.  Those are summarized
repeatedly across independent jobs on different nodes, and the merge requires
all shards to agree on them bitwise.  That extends the separation's guarantee from
"within one run" to "across independent jobs", which is the property the sharded
design actually needs and which a single serial pass would never test.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

SERIES = "GSE83452"
EXPECTED_RECORDS = 231


class ShardPlanError(RuntimeError):
    """Raised when the shard plan would not cover the cohort exactly once."""


def read_accessions(manifest: Path) -> list[str]:
    with manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "sample_accession" not in reader.fieldnames:
            raise ShardPlanError("manifest lacks sample_accession")
        rows = [row["sample_accession"] for row in reader if row["series"] == SERIES]
    if len(rows) != EXPECTED_RECORDS:
        raise ShardPlanError(f"{SERIES} manifest is not {EXPECTED_RECORDS} records")
    if len(set(rows)) != len(rows):
        raise ShardPlanError("manifest repeats a sample accession")
    return sorted(rows)


def plan_shards(
    accessions: list[str], *, shards: int, sentinels: int
) -> tuple[list[list[str]], list[str]]:
    if shards < 1:
        raise ShardPlanError("shard count must be positive")
    if not 0 < sentinels < len(accessions):
        raise ShardPlanError("sentinel count is out of range")

    # Sentinels are spread across the sorted axis so they are not all from one
    # end of the accession range.
    step = len(accessions) // sentinels
    sentinel_positions = sorted({min(index * step, len(accessions) - 1) for index in range(sentinels)})
    sentinel_list = [accessions[position] for position in sentinel_positions]

    # Round-robin assignment keeps shards balanced without regard to content.
    buckets: list[list[str]] = [[] for _ in range(shards)]
    for index, accession in enumerate(accessions):
        buckets[index % shards].append(accession)

    covered = [accession for bucket in buckets for accession in bucket]
    if sorted(covered) != accessions:
        raise ShardPlanError("shard plan does not cover the cohort exactly once")

    # Sentinels are appended to every shard that does not already own them, so
    # each sentinel is summarized in every shard.
    plan = []
    for bucket in buckets:
        owned = set(bucket)
        extra = [accession for accession in sentinel_list if accession not in owned]
        plan.append(sorted(bucket + extra))
    return plan, sentinel_list


def write_plan(
    *, manifest: Path, output: Path, shards: int, sentinels: int
) -> dict[str, object]:
    accessions = read_accessions(manifest)
    plan, sentinel_list = plan_shards(accessions, shards=shards, sentinels=sentinels)
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for index, bucket in enumerate(plan):
        path = output / f"shard_{index:02d}.txt"
        path.write_text("\n".join(bucket) + "\n", encoding="utf-8")
        owned = [a for a in bucket if a not in sentinel_list]
        records.append(
            {
                "shard": index,
                "path": str(path),
                "arrays": len(bucket),
                "owned_arrays": len(owned),
                "sentinel_arrays": len([a for a in bucket if a in sentinel_list]),
            }
        )
    receipt = {
        "schema_version": "masld-bench-gse83452-shard-plan-v1",
        "status": "pass_shard_plan_covers_cohort_exactly_once",
        "series": SERIES,
        "platform_id": "GPL16686",
        "cohort_family_id": "antwerp_inserm_shared",
        "total_arrays": len(accessions),
        "shards": len(plan),
        "sentinels": sentinel_list,
        "sentinel_count": len(sentinel_list),
        "assignment_rule": "round_robin_over_sorted_GSM_accession",
        "assignment_used_labels": False,
        "assignment_used_intensity": False,
        "assignment_used_timepoint_or_participant": False,
        "sentinels_present_in_every_shard": True,
        "cross_job_reproducibility_checked_by": "sentinel digest agreement across shards",
        "shard_records": records,
        "labels_read": False,
        "model_training_activated": False,
    }
    (output / "shard_plan.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shards", type=int, default=8)
    parser.add_argument("--sentinels", type=int, default=3)
    arguments = parser.parse_args()
    receipt = write_plan(
        manifest=arguments.manifest,
        output=arguments.output,
        shards=arguments.shards,
        sentinels=arguments.sentinels,
    )
    print(json.dumps({k: receipt[k] for k in receipt if k != "shard_records"}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

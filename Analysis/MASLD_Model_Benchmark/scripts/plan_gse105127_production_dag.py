#!/usr/bin/env python3
"""Plan participant-contained GSE105127 RRBS and RNA production bundles."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


EXPECTED_ACTIVATION_ARTIFACTS_SHA256 = (
    "fc529784b5f6df53af7ba82036ff617507ca9d551d64d79aac6cb9061faf62bb"
)
EXPECTED_RRBS_ARTIFACTS_SHA256 = (
    "5e3f698a5172d612b38bfe0a78b57f7199c10031ae9740190ead1835ef90e216"
)
EXPECTED_ROWS = 57
EXPECTED_PARTICIPANTS = 19
EXPECTED_RNA_BYTES = 79_576_393_854
EXPECTED_RNA_READS = 1_140_779_295
EXPECTED_RNA_BASES = 86_699_226_420
EXPECTED_RRBS_BYTES = 8_192_644_943
FORBIDDEN = {
    "phenotype", "phenotype_name", "label", "outcome", "disease",
    "fibrosis", "nas", "sex", "age", "bmi", "outer_fold",
}


class GSE105127PlanError(RuntimeError):
    """Raised when the label-free production plan differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE105127PlanError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise GSE105127PlanError(f"refusing to write an empty plan: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = tuple(rows[0])
    if FORBIDDEN & set(fields):
        raise GSE105127PlanError("outcome-bearing field entered a production plan")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def assign_participant_bundles(
    rows: Sequence[Mapping[str, Any]], bundle_count: int
) -> tuple[dict[str, int], list[dict[str, Any]]]:
    if bundle_count < 1:
        raise GSE105127PlanError("bundle count must be positive")
    weights: defaultdict[str, int] = defaultdict(int)
    row_counts: defaultdict[str, int] = defaultdict(int)
    for row in rows:
        participant = str(row["participant_group_id"])
        weights[participant] += int(row["combined_source_bytes"])
        row_counts[participant] += 1
    if len(weights) != EXPECTED_PARTICIPANTS or set(row_counts.values()) != {3}:
        raise GSE105127PlanError("participant-zone grouping differs")
    if bundle_count > len(weights):
        raise GSE105127PlanError("bundle count exceeds participant count")
    bundles = [
        {"bundle_id": index, "source_bytes": 0, "participants": []}
        for index in range(bundle_count)
    ]
    for participant in sorted(weights, key=lambda value: (-weights[value], value)):
        target = min(
            bundles,
            key=lambda value: (
                int(value["source_bytes"]), len(value["participants"]), int(value["bundle_id"])
            ),
        )
        target["participants"].append(participant)
        target["source_bytes"] += weights[participant]
    assignment = {
        participant: int(bundle["bundle_id"])
        for bundle in bundles
        for participant in bundle["participants"]
    }
    summary = [
        {
            "bundle_id": int(bundle["bundle_id"]),
            "participants": len(bundle["participants"]),
            "participant_zone_rows": sum(row_counts[value] for value in bundle["participants"]),
            "combined_source_bytes": int(bundle["source_bytes"]),
        }
        for bundle in bundles
    ]
    return assignment, summary


def build_plan(
    *, activation_root: Path,
    rrbs_root: Path,
    output: Path,
    bundle_count: int,
) -> dict[str, Any]:
    if output.exists():
        raise GSE105127PlanError(f"refusing to overwrite plan: {output}")
    verify_frozen_tree(activation_root)
    if sha256_file(activation_root / "ARTIFACTS.json") != EXPECTED_ACTIVATION_ARTIFACTS_SHA256:
        raise GSE105127PlanError(f"frozen activation identity differs: {activation_root}")
    # The RRBS execution also contains a historical QC table with phenotypes.
    # Production must never open that table. Bind its read-only ARTIFACTS file,
    # then resolve raw members only from the checksum inventory and label-free
    # activation manifest.
    rrbs_artifacts_path = rrbs_root / "ARTIFACTS.json"
    if sha256_file(rrbs_artifacts_path) != EXPECTED_RRBS_ARTIFACTS_SHA256:
        raise GSE105127PlanError(f"frozen RRBS inventory identity differs: {rrbs_root}")
    completion = json.loads((rrbs_root / "COMPLETE").read_text(encoding="utf-8"))
    if (
        completion.get("schema_version") != "masld-bench-complete-v1"
        or completion.get("manifest_sha256") != EXPECTED_RRBS_ARTIFACTS_SHA256
    ):
        raise GSE105127PlanError("RRBS completion marker differs")
    rrbs_artifacts = json.loads(rrbs_artifacts_path.read_text(encoding="utf-8"))
    raw_inventory = {
        Path(str(item["path"])).name.split("__", 1)[0]: item
        for item in rrbs_artifacts.get("artifacts", [])
        if str(item.get("path", "")).startswith("raw/")
        and str(item.get("path", "")).endswith(".CG.bed.gz")
    }
    if len(raw_inventory) != EXPECTED_ROWS:
        raise GSE105127PlanError("RRBS raw checksum inventory differs")
    fixture = activation_root / "fixture"
    axis_fields, axis = read_tsv(fixture / "model_inputs/row_axis.tsv")
    rna_fields, rna = read_tsv(fixture / "authorities/rna_run_manifest.tsv")
    rrbs_fields, rrbs = read_tsv(fixture / "authorities/rrbs_source_manifest.tsv")
    for fields in (axis_fields, rna_fields, rrbs_fields):
        if FORBIDDEN & set(fields):
            raise GSE105127PlanError("activation fixture violates the label firewall")
    if len(axis) != EXPECTED_ROWS or len(rna) != EXPECTED_ROWS or len(rrbs) != EXPECTED_ROWS:
        raise GSE105127PlanError("activation row census differs")
    by_row = {row["row_id"]: row for row in axis}
    rna_by_row = {row["row_id"]: row for row in rna}
    rrbs_by_row = {row["row_id"]: row for row in rrbs}
    if set(by_row) != set(rna_by_row) or set(by_row) != set(rrbs_by_row):
        raise GSE105127PlanError("RNA/RRBS row axes differ")
    combined = []
    for row_key in sorted(by_row):
        axis_row = by_row[row_key]
        rna_row = rna_by_row[row_key]
        rrbs_row = rrbs_by_row[row_key]
        member = raw_inventory.get(rrbs_row["sample_accession"])
        if (
            member is None
            or axis_row["pairing_topology"] != "adjacent_section"
            or axis_row["same_section"] != "False"
            or rna_row["model_input_state"] != "derivable_not_processed"
            or rrbs_row["coordinate_state"] != "native_hg19_remap_pending"
        ):
            raise GSE105127PlanError("source state differs from no-fit activation")
        local = rrbs_root / str(member["path"])
        if (
            str(member["sha256"]) != rrbs_row["processed_bed_sha256"]
            or int(member["size_bytes"]) != int(rrbs_row["processed_bed_bytes"])
            or not local.is_file()
            or local.stat().st_size != int(member["size_bytes"])
        ):
            raise GSE105127PlanError("frozen local RRBS member is absent")
        combined.append(
            {
                "row_id": row_key,
                "participant_group_id": axis_row["participant_group_id"],
                "zone": axis_row["zone"],
                "pairing_topology": "adjacent_section",
                "rna_fastq_bytes": int(rna_row["fastq_bytes"]),
                "rrbs_bed_bytes": int(rrbs_row["processed_bed_bytes"]),
                "combined_source_bytes": int(rna_row["fastq_bytes"])
                + int(rrbs_row["processed_bed_bytes"]),
            }
        )
    if (
        sum(int(row["rna_fastq_bytes"]) for row in combined) != EXPECTED_RNA_BYTES
        or sum(int(row["rrbs_bed_bytes"]) for row in combined) != EXPECTED_RRBS_BYTES
    ):
        raise GSE105127PlanError("source byte census differs")
    assignment, bundle_summary = assign_participant_bundles(combined, bundle_count)
    rna_plan = []
    rrbs_plan = []
    for row_key in sorted(by_row):
        axis_row = by_row[row_key]
        rna_row = rna_by_row[row_key]
        rrbs_row = rrbs_by_row[row_key]
        member = raw_inventory[rrbs_row["sample_accession"]]
        common = {
            "bundle_id": assignment[axis_row["participant_group_id"]],
            "row_id": row_key,
            "participant_group_id": axis_row["participant_group_id"],
            "zone": axis_row["zone"],
            "pairing_topology": "adjacent_section",
        }
        rna_plan.append(
            {
                **common,
                "sample_accession": rna_row["sample_accession"],
                "run_accession": rna_row["run_accession"],
                "library_layout": "SINGLE",
                "read_length": 76,
                "read_count": int(rna_row["read_count"]),
                "base_count": int(rna_row["base_count"]),
                "fastq_url": rna_row["fastq_url"],
                "fastq_bytes": int(rna_row["fastq_bytes"]),
                "fastq_md5": rna_row["fastq_md5"],
                "relative_fastq_path": f"participants/{row_key}/reads.fastq.gz",
            }
        )
        rrbs_plan.append(
            {
                **common,
                "sample_accession": rrbs_row["sample_accession"],
                "source_bed_path": str((rrbs_root / str(member["path"])).resolve()),
                "source_bed_bytes": int(rrbs_row["processed_bed_bytes"]),
                "source_bed_sha256": rrbs_row["processed_bed_sha256"],
                "source_cytosine_rows": int(rrbs_row["cytosine_rows"]),
                "relative_collapsed_path": f"participants/{row_key}/cpg_counts.hg19.tsv.gz",
            }
        )
    rna_plan.sort(key=lambda row: (int(row["bundle_id"]), row["participant_group_id"], row["zone"]))
    rrbs_plan.sort(key=lambda row: (int(row["bundle_id"]), row["participant_group_id"], row["zone"]))
    if (
        sum(int(row["read_count"]) for row in rna_plan) != EXPECTED_RNA_READS
        or sum(int(row["base_count"]) for row in rna_plan) != EXPECTED_RNA_BASES
        or any(int(row["base_count"]) != 76 * int(row["read_count"]) for row in rna_plan)
    ):
        raise GSE105127PlanError("RNA read/base census differs")
    output.mkdir(parents=True)
    write_tsv(output / "bundles.tsv", bundle_summary)
    write_tsv(output / "rna_rows.tsv", rna_plan)
    write_tsv(output / "rrbs_rows.tsv", rrbs_plan)
    references = {
        "schema_version": "masld-bench-gse105127-reference-plan-v1",
        "source_fasta_url": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz",
        "source_fasta_compressed_bytes": 892_331_003,
        "source_fasta_uncompressed_md5": "0ce84c872fc0072a885926823dcd0338",
        "forward_chain_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        "forward_chain_md5": "35887f73fe5e2231656504d1f6430900",
        "reverse_chain_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
        "reverse_chain_md5": "ff3031d93792f4cbb86af44055efd903",
        "target_fasta_path": "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz",
        "target_fasta_sha256": "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215",
        "target_gtf_path": "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz",
        "target_gtf_sha256": "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9",
    }
    (output / "references.json").write_text(
        json.dumps(references, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt = {
        "schema_version": "masld-bench-gse105127-production-plan-v1",
        "status": "planned_no_fit",
        "participants": EXPECTED_PARTICIPANTS,
        "participant_zone_rows": EXPECTED_ROWS,
        "bundles": bundle_count,
        "participant_split_across_bundles": False,
        "adjacent_section_only": True,
        "labels_accessed": False,
        "fit_or_score_performed": False,
        "rna_fastq_bytes": EXPECTED_RNA_BYTES,
        "rna_reads": EXPECTED_RNA_READS,
        "rna_bases": EXPECTED_RNA_BASES,
        "rrbs_bed_bytes": EXPECTED_RRBS_BYTES,
        "activation_artifacts_sha256": EXPECTED_ACTIVATION_ARTIFACTS_SHA256,
        "rrbs_artifacts_sha256": EXPECTED_RRBS_ARTIFACTS_SHA256,
        "historical_outcome_qc_table_opened": False,
    }
    (output / "plan.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(output, {"artifact_class": "gse105127_production_plan", "status": "passed"})
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activation-root", type=Path, required=True)
    parser.add_argument("--rrbs-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bundles", type=int, default=8)
    args = parser.parse_args()
    print(json.dumps(build_plan(
        activation_root=args.activation_root,
        rrbs_root=args.rrbs_root,
        output=args.output,
        bundle_count=args.bundles,
    ), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Create deterministic participant-contained download bundles for GSE268273."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import verify_frozen_tree


EXPECTED_PREPROCESSING_ARTIFACTS_SHA256 = (
    "2892cb792064719a95f1ac7edde917bff7283b97d8347a378081f84a49257374"
)
EXPECTED_MODEL_INPUT_ARTIFACTS_SHA256 = (
    "75047dff8fd1e1bc994bf1f915f80c81c87038af5a78c65931fe2faa6c6d2294"
)
EXPECTED_ACTIVATION_ARTIFACTS_SHA256 = (
    "d9a9e8b1a99739908587f0b854bb3adc472a3543b769eeeb0813544e78d57c55"
)
EXPECTED_PARTICIPANTS = 109
EXPECTED_RUNS = 824
EXPECTED_FASTQ_FILES = 968
EXPECTED_FASTQ_BYTES = 512_881_099_727
EXPECTED_FILE_COUNT_DISTRIBUTION = {1: 680, 2: 144}
EXPECTED_EFFECTIVE_LAYOUT_RUNS = {"single_end": 680, "paired_end": 144}
EXPECTED_EFFECTIVE_LAYOUT_PARTICIPANTS = {"single_end": 73, "paired_end": 36}
GENCODE_V36_GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/"
    "gencode.v36.primary_assembly.annotation.gtf.gz"
)
GENCODE_V36_GTF_MD5 = "5038acd7158686c69a1f88a5e26a5934"
GENCODE_V36_GTF_BYTES = 44_517_553
GENOME_PATH = Path(
    "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/"
    "genome/GRCh38.p14.genome.fa.gz"
)
GENOME_SHA256 = "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215"
GENCODE_V49_GTF_PATH = Path(
    "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/"
    "gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
)
GENCODE_V49_GTF_SHA256 = (
    "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
)
PROHIBITED_REFERENCE_PATHS = (
    "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/GRCh38.p14.genome.fa.gz",
    "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.transcripts.fa.gz",
)


class GSE268273RawPlanError(RuntimeError):
    """Raised when a raw campaign cannot be planned without ambiguity."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE268273RawPlanError(f"TSV lacks a header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def expand_run_manifest(
    rows: Sequence[Mapping[str, str]],
    ena_read_rows: Sequence[Mapping[str, str]],
    *,
    production_contract: bool = True,
) -> list[dict[str, Any]]:
    required = {
        "row_id",
        "run_accession",
        "experiment_accession",
        "biosample_accession",
        "library_layout",
        "fastq_ftp",
        "fastq_bytes",
        "fastq_md5",
        "replicate_role",
    }
    if not rows or not required <= set(rows[0]):
        raise GSE268273RawPlanError("raw run-manifest schema differs")
    forbidden = {
        "fibrosis",
        "label",
        "outcome",
        "disease",
        "source_group",
        "sex",
        "nas",
        "bmi",
    }
    if forbidden & set(rows[0]):
        raise GSE268273RawPlanError("raw run manifest violates the label firewall")
    ena_required = {
        "run_accession",
        "experiment_accession",
        "library_layout",
        "read_count",
        "base_count",
        "fastq_ftp",
        "fastq_bytes",
        "fastq_md5",
    }
    if not ena_read_rows or not ena_required <= set(ena_read_rows[0]):
        raise GSE268273RawPlanError("ENA read-metadata schema differs")
    if forbidden & set(ena_read_rows[0]):
        raise GSE268273RawPlanError("ENA read metadata violates the label firewall")
    runs = [row["run_accession"] for row in rows]
    if len(runs) != len(set(runs)):
        raise GSE268273RawPlanError("raw run accession is duplicated")
    ena_by_run = {row["run_accession"]: row for row in ena_read_rows}
    if len(ena_by_run) != len(ena_read_rows) or set(ena_by_run) != set(runs):
        raise GSE268273RawPlanError("ENA read-metadata run axis differs")
    expanded: list[dict[str, Any]] = []
    file_count_distribution: defaultdict[int, int] = defaultdict(int)
    effective_layout_runs: defaultdict[str, int] = defaultdict(int)
    experiment_by_participant: dict[str, str] = {}
    biosample_by_participant: dict[str, str] = {}
    effective_layout_by_participant: dict[str, str] = {}
    for row in rows:
        if (
            row["library_layout"] != "SINGLE"
            or row["replicate_role"] != "technical_run_partition"
        ):
            raise GSE268273RawPlanError("ENA layout or technical-run role differs")
        row_id = row["row_id"]
        experiment = row["experiment_accession"]
        biosample = row["biosample_accession"]
        if row_id in experiment_by_participant and (
            experiment_by_participant[row_id] != experiment
            or biosample_by_participant[row_id] != biosample
        ):
            raise GSE268273RawPlanError(
                "one participant maps to multiple experiment or BioSample identifiers"
            )
        experiment_by_participant[row_id] = experiment
        biosample_by_participant[row_id] = biosample
        ena = ena_by_run[row["run_accession"]]
        if (
            any(
                ena[field] != row[field]
                for field in (
                    "experiment_accession",
                    "library_layout",
                    "fastq_ftp",
                    "fastq_bytes",
                )
            )
            or ena["fastq_md5"].lower() != row["fastq_md5"].lower()
        ):
            raise GSE268273RawPlanError("frozen and current ENA run evidence differs")
        paths = row["fastq_ftp"].split(";")
        sizes_text = row["fastq_bytes"].split(";")
        md5s = row["fastq_md5"].lower().split(";")
        if not paths or len(paths) != len(sizes_text) or len(paths) != len(md5s):
            raise GSE268273RawPlanError("parallel ENA FASTQ fields differ")
        try:
            read_count = int(ena["read_count"])
            base_count = int(ena["base_count"])
        except ValueError as error:
            raise GSE268273RawPlanError(
                "ENA read_count or base_count is not an integer"
            ) from error
        names = [PurePosixPath(path).name for path in paths]
        run = row["run_accession"]
        if (
            len(paths) == 1
            and names == [f"{run}.fastq.gz"]
            and read_count > 0
            and read_count * 101 - 100 <= base_count <= read_count * 101
        ):
            effective_layout = "single_end"
            read_roles = ["single_end_read"]
            layout_evidence = "one_fastq_and_approximately_101_bases_per_ENA_read_count"
        elif (
            len(paths) == 2
            and names == [f"{run}_1.fastq.gz", f"{run}_2.fastq.gz"]
            and read_count > 0
            and base_count == read_count * 202
        ):
            effective_layout = "paired_end"
            read_roles = ["paired_end_mate_1", "paired_end_mate_2"]
            layout_evidence = (
                "two_mate_fastqs_and_202_bases_per_ENA_read_count_"
                "despite_declared_SINGLE"
            )
        else:
            raise GSE268273RawPlanError(
                "FASTQ topology cannot resolve the effective read layout"
            )
        if row_id in effective_layout_by_participant and (
            effective_layout_by_participant[row_id] != effective_layout
        ):
            raise GSE268273RawPlanError(
                "one participant contains mixed effective read layouts"
            )
        effective_layout_by_participant[row_id] = effective_layout
        effective_layout_runs[effective_layout] += 1
        file_count_distribution[len(paths)] += 1
        for file_index, (ftp_path, size_text, md5, read_role) in enumerate(
            zip(paths, sizes_text, md5s, read_roles, strict=True)
        ):
            try:
                size = int(size_text)
            except ValueError as error:
                raise GSE268273RawPlanError("ENA FASTQ size is not an integer") from error
            if (
                not ftp_path.startswith("ftp.sra.ebi.ac.uk/")
                or PurePosixPath(ftp_path).name in {"", ".", ".."}
                or size <= 0
                or len(md5) != 32
                or any(character not in "0123456789abcdef" for character in md5)
            ):
                raise GSE268273RawPlanError("ENA FASTQ member identity differs")
            expanded.append(
                {
                    "row_id": row_id,
                    "run_accession": row["run_accession"],
                    "experiment_accession": experiment,
                    "biosample_accession": biosample,
                    "declared_ena_library_layout": "SINGLE",
                    "effective_library_layout": effective_layout,
                    "layout_evidence": layout_evidence,
                    "ena_read_count": read_count,
                    "ena_base_count": base_count,
                    "file_part_index": file_index,
                    "file_part_count": len(paths),
                    "read_role": read_role,
                    "https_url": "https://" + ftp_path,
                    "relative_path": (
                        f"fastq/{row['run_accession']}/{PurePosixPath(ftp_path).name}"
                    ),
                    "fastq_bytes": size,
                    "fastq_md5": md5,
                }
            )
    relative_paths = [row["relative_path"] for row in expanded]
    urls = [row["https_url"] for row in expanded]
    if len(relative_paths) != len(set(relative_paths)) or len(urls) != len(set(urls)):
        raise GSE268273RawPlanError("expanded FASTQ member is duplicated")
    if production_contract and (
        len(rows) != EXPECTED_RUNS
        or len(experiment_by_participant) != EXPECTED_PARTICIPANTS
        or len(set(experiment_by_participant.values())) != EXPECTED_PARTICIPANTS
        or len(set(biosample_by_participant.values())) != EXPECTED_PARTICIPANTS
        or len(expanded) != EXPECTED_FASTQ_FILES
        or sum(int(row["fastq_bytes"]) for row in expanded) != EXPECTED_FASTQ_BYTES
        or dict(file_count_distribution) != EXPECTED_FILE_COUNT_DISTRIBUTION
        or dict(effective_layout_runs) != EXPECTED_EFFECTIVE_LAYOUT_RUNS
        or dict(Counter(effective_layout_by_participant.values()))
        != EXPECTED_EFFECTIVE_LAYOUT_PARTICIPANTS
    ):
        raise GSE268273RawPlanError("production ENA topology or byte census differs")
    return expanded


def assign_participant_bundles(
    files: Sequence[Mapping[str, Any]], bundle_count: int
) -> tuple[
    list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]
]:
    if bundle_count < 1:
        raise GSE268273RawPlanError("bundle count must be positive")
    by_participant: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in files:
        by_participant[str(row["row_id"])].append(row)
    if bundle_count > len(by_participant):
        raise GSE268273RawPlanError("bundle count exceeds participant count")
    totals = {
        row_id: sum(int(row["fastq_bytes"]) for row in members)
        for row_id, members in by_participant.items()
    }
    bundles: list[dict[str, Any]] = [
        {"bundle_id": index, "bytes": 0, "participants": []}
        for index in range(bundle_count)
    ]
    for row_id in sorted(totals, key=lambda value: (-totals[value], value)):
        target = min(
            bundles,
            key=lambda value: (
                int(value["bytes"]),
                len(value["participants"]),
                int(value["bundle_id"]),
            ),
        )
        target["participants"].append(row_id)
        target["bytes"] += totals[row_id]
    assignment = {
        row_id: int(bundle["bundle_id"])
        for bundle in bundles
        for row_id in bundle["participants"]
    }
    file_rows = [
        {"bundle_id": assignment[str(row["row_id"])], **row}
        for row in files
    ]
    file_rows.sort(
        key=lambda row: (
            int(row["bundle_id"]),
            str(row["row_id"]),
            str(row["run_accession"]),
            int(row["file_part_index"]),
        )
    )
    participant_rows = []
    for row_id in sorted(by_participant, key=lambda value: (assignment[value], value)):
        members = by_participant[row_id]
        participant_rows.append(
            {
                "bundle_id": assignment[row_id],
                "row_id": row_id,
                "experiment_accession": next(
                    iter({str(row["experiment_accession"]) for row in members})
                ),
                "biosample_accession": next(
                    iter({str(row["biosample_accession"]) for row in members})
                ),
                "technical_runs": len({str(row["run_accession"]) for row in members}),
                "fastq_files": len(members),
                "fastq_bytes": sum(int(row["fastq_bytes"]) for row in members),
                "effective_library_layout": next(
                    iter({str(row["effective_library_layout"]) for row in members})
                ),
                "aggregation_rule": (
                    "matewise_comma_ordered_technical_runs_then_one_participant_RSEM_library"
                    if next(
                        iter({str(row["effective_library_layout"]) for row in members})
                    )
                    == "paired_end"
                    else "comma_ordered_single_end_technical_runs_then_one_participant_RSEM_library"
                ),
            }
        )
    summary_rows = []
    for bundle in bundles:
        bundle_id = int(bundle["bundle_id"])
        members = [row for row in file_rows if int(row["bundle_id"]) == bundle_id]
        summary_rows.append(
            {
                "bundle_id": bundle_id,
                "participants": len(bundle["participants"]),
                "technical_runs": len({row["run_accession"] for row in members}),
                "fastq_files": len(members),
                "fastq_bytes": sum(int(row["fastq_bytes"]) for row in members),
            }
        )
    return file_rows, participant_rows, summary_rows


def build_plan(
    *,
    preprocessing_root: Path,
    model_input_root: Path,
    activation_root: Path,
    ena_read_metadata: Path,
    output: Path,
    bundle_count: int,
) -> dict[str, Any]:
    if output.exists():
        raise GSE268273RawPlanError(f"refusing to overwrite campaign plan: {output}")
    for path, expected in (
        (preprocessing_root, EXPECTED_PREPROCESSING_ARTIFACTS_SHA256),
        (model_input_root, EXPECTED_MODEL_INPUT_ARTIFACTS_SHA256),
        (activation_root, EXPECTED_ACTIVATION_ARTIFACTS_SHA256),
    ):
        verify_frozen_tree(path)
        if sha256_file(path / "ARTIFACTS.json") != expected:
            raise GSE268273RawPlanError(f"frozen artifact identity differs: {path}")
    _, run_rows = read_tsv(preprocessing_root / "raw_run_manifest.tsv")
    _, ena_rows = read_tsv(ena_read_metadata)
    files = expand_run_manifest(run_rows, ena_rows)
    file_rows, participant_rows, bundle_rows = assign_participant_bundles(
        files, bundle_count
    )
    if (
        len(participant_rows) != EXPECTED_PARTICIPANTS
        or len(bundle_rows) != bundle_count
        or sum(int(row["fastq_bytes"]) for row in bundle_rows)
        != EXPECTED_FASTQ_BYTES
    ):
        raise GSE268273RawPlanError("bundle accounting differs")
    output.mkdir(parents=True)
    write_tsv(output / "fastq_files.tsv", tuple(file_rows[0]), file_rows)
    write_tsv(output / "participants.tsv", tuple(participant_rows[0]), participant_rows)
    write_tsv(output / "bundles.tsv", tuple(bundle_rows[0]), bundle_rows)
    source_evidence = output / "source_evidence"
    source_evidence.mkdir()
    shutil.copyfile(ena_read_metadata, source_evidence / "ena_read_metadata.tsv")
    reference = {
        "schema_version": "masld-bench-gse268273-rsem-reference-v1",
        "alignment_genome": {
            "assembly": "GRCh38.p14",
            "path": str(GENOME_PATH),
            "sha256": GENOME_SHA256,
            "required_checks": [
                "gzip_test",
                "compressed_sha256",
                "fasta_contig_inventory",
            ],
        },
        "source_quantification_annotation": {
            "release": "GENCODE_v36_GRCh38p13_primary_assembly",
            "url": GENCODE_V36_GTF_URL,
            "bytes": GENCODE_V36_GTF_BYTES,
            "md5": GENCODE_V36_GTF_MD5,
            "purpose": (
                "source_native_gene_ids_before_v49_crosswalk_including_the_"
                "seven_frozen_primary_scaffold_features"
            ),
        },
        "target_annotation": {
            "release": "GENCODE_v49_GRCh38p14",
            "path": str(GENCODE_V49_GTF_PATH),
            "sha256": GENCODE_V49_GTF_SHA256,
            "required_checks": [
                "gzip_test",
                "compressed_sha256",
                "source_to_target_crosswalk_membership",
            ],
        },
        "prohibited_reference_paths": list(PROHIBITED_REFERENCE_PATHS),
        "quantifier": "RSEM_1.3.1_observed_binary_expected_count_with_STAR_2.7.10b",
        "runtime_module_label": "RSEM_1.3.3_foss_2022a",
        "runtime_module_label_differs_from_binary_reported_version": True,
        "source_paper_quantifier": "RSEM_1.3.0_with_STAR_2.5.3a",
        "version_difference_is_transport_limitation": True,
        "alignment_parameters": "RSEM_integrated_STAR_defaults",
        "source_paper_alignment_description": "STAR_with_ENCODE_parameters",
        "source_exact_alignment_invocation_available": False,
        "single_end_fragment_distribution_parameters": "RSEM_defaults_source_did_not_deposit_values",
        "strandedness": "unstranded_source_not_deposited",
    }
    (output / "reference_contract.json").write_text(
        json.dumps(reference, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    canonical = json.dumps(
        {"algorithm": "participant_lpt_bytes_v1", "files": file_rows},
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    campaign_hash = hashlib.sha256(canonical).hexdigest()
    receipt = {
        "schema_version": "masld-bench-gse268273-raw-campaign-plan-v1",
        "status": "planned_not_submitted",
        "campaign_hash": campaign_hash,
        "bundle_algorithm": "largest_participant_first_to_lowest_byte_total_then_bundle_id",
        "bundles": bundle_count,
        "participants": len(participant_rows),
        "technical_runs": len(run_rows),
        "fastq_files": len(file_rows),
        "single_file_runs": EXPECTED_FILE_COUNT_DISTRIBUTION[1],
        "declared_single_two_fastq_runs": EXPECTED_FILE_COUNT_DISTRIBUTION[2],
        "effective_single_end_runs": EXPECTED_EFFECTIVE_LAYOUT_RUNS["single_end"],
        "effective_paired_end_runs": EXPECTED_EFFECTIVE_LAYOUT_RUNS["paired_end"],
        "effective_single_end_participants": EXPECTED_EFFECTIVE_LAYOUT_PARTICIPANTS[
            "single_end"
        ],
        "effective_paired_end_participants": EXPECTED_EFFECTIVE_LAYOUT_PARTICIPANTS[
            "paired_end"
        ],
        "declared_layout_discrepancy_admitted": True,
        "fastq_bytes": EXPECTED_FASTQ_BYTES,
        "fastq_gib": EXPECTED_FASTQ_BYTES / (1024**3),
        "minimum_bundle_bytes": min(int(row["fastq_bytes"]) for row in bundle_rows),
        "maximum_bundle_bytes": max(int(row["fastq_bytes"]) for row in bundle_rows),
        "participant_contained_in_one_bundle": True,
        "labels_included": False,
        "outcomes_accessed": False,
        "fit_or_score_performed": False,
        "source_preprocessing_artifacts_sha256": EXPECTED_PREPROCESSING_ARTIFACTS_SHA256,
        "source_model_input_artifacts_sha256": EXPECTED_MODEL_INPUT_ARTIFACTS_SHA256,
        "source_activation_artifacts_sha256": EXPECTED_ACTIVATION_ARTIFACTS_SHA256,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preprocessing-root", type=Path, required=True)
    parser.add_argument("--model-input-root", type=Path, required=True)
    parser.add_argument("--activation-root", type=Path, required=True)
    parser.add_argument("--ena-read-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bundles", type=int, default=8)
    args = parser.parse_args()
    receipt = build_plan(
        preprocessing_root=args.preprocessing_root,
        model_input_root=args.model_input_root,
        activation_root=args.activation_root,
        ena_read_metadata=args.ena_read_metadata,
        output=args.output,
        bundle_count=args.bundles,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

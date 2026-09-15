#!/usr/bin/env python3
"""Quantify one GSE268273 raw bundle as participant-level RSEM expected counts."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


class GSE268273QuantificationError(RuntimeError):
    """Raised when participant-level raw quantification is not eligible."""


def tool_environment() -> dict[str, str]:
    path = os.environ.get("GSE268273_TOOL_PATH", "")
    if not path:
        raise GSE268273QuantificationError("isolated RSEM/STAR tool PATH is absent")
    environment = dict(os.environ)
    environment["PATH"] = path
    library_path = os.environ.get("GSE268273_TOOL_LD_LIBRARY_PATH", "")
    if library_path:
        environment["LD_LIBRARY_PATH"] = library_path
    else:
        environment.pop("LD_LIBRARY_PATH", None)
    environment.pop("PYTHONPATH", None)
    return environment


def tool_executable(name: str, environment: Mapping[str, str]) -> str:
    executable = shutil.which(name, path=environment["PATH"])
    if executable is None:
        raise GSE268273QuantificationError(
            f"isolated tool executable is absent: {name}"
        )
    return executable


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
            raise GSE268273QuantificationError(f"TSV lacks a header: {path}")
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


def parse_rsem_gene_results(path: Path) -> dict[str, float]:
    fields, rows = read_tsv(path)
    required = {
        "gene_id",
        "transcript_id(s)",
        "length",
        "effective_length",
        "expected_count",
        "TPM",
        "FPKM",
    }
    if not required <= set(fields) or not rows:
        raise GSE268273QuantificationError("RSEM gene-results schema differs")
    output: dict[str, float] = {}
    for row in rows:
        gene_id = row["gene_id"]
        if gene_id in output:
            raise GSE268273QuantificationError("RSEM gene_id is duplicated")
        try:
            value = float(row["expected_count"])
        except ValueError as error:
            raise GSE268273QuantificationError(
                "RSEM expected_count is not numeric"
            ) from error
        if not np.isfinite(value) or value < 0:
            raise GSE268273QuantificationError("RSEM expected_count is invalid")
        output[gene_id] = value
    return output


def aggregate_v49_expected_counts(
    rsem_counts: Mapping[str, float], candidate_rows: Sequence[Mapping[str, str]]
) -> np.ndarray:
    values = np.empty(len(candidate_rows), dtype=np.float64)
    for index, row in enumerate(candidate_rows):
        members = row["source_feature_ids"].split(";")
        if len(members) != int(row["source_feature_count"]):
            raise GSE268273QuantificationError("candidate source-member count differs")
        missing = [member for member in members if member not in rsem_counts]
        if missing:
            raise GSE268273QuantificationError(
                "RSEM lacks a source feature required by the v49 crosswalk"
            )
        if len(members) == 1:
            if row["raw_count_aggregation"] != "identity":
                raise GSE268273QuantificationError("one-to-one aggregation rule differs")
        elif row["raw_count_aggregation"] != "deterministic_sum_before_normalization":
            raise GSE268273QuantificationError("duplicate aggregation rule differs")
        values[index] = sum(float(rsem_counts[member]) for member in members)
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise GSE268273QuantificationError("aggregated v49 count vector is invalid")
    return values


def arrange_participant_fastqs(
    members: Sequence[Mapping[str, str]],
) -> tuple[str, list[Mapping[str, str]], list[Mapping[str, str]]]:
    """Return RSEM upstream/downstream inputs without confusing mates with lanes."""

    layouts = {row["effective_library_layout"] for row in members}
    if len(layouts) != 1:
        raise GSE268273QuantificationError(
            "one participant contains mixed effective read layouts"
        )
    layout = next(iter(layouts))
    by_run: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in members:
        by_run[row["run_accession"]].append(row)
    upstream: list[Mapping[str, str]] = []
    downstream: list[Mapping[str, str]] = []
    for run in sorted(by_run):
        run_members = sorted(
            by_run[run], key=lambda row: int(row["file_part_index"])
        )
        roles = [row["read_role"] for row in run_members]
        indices = [int(row["file_part_index"]) for row in run_members]
        counts = {int(row["file_part_count"]) for row in run_members}
        if layout == "single_end":
            if (
                len(run_members) != 1
                or roles != ["single_end_read"]
                or indices != [0]
                or counts != {1}
            ):
                raise GSE268273QuantificationError(
                    "single-end technical-run topology differs"
                )
            upstream.append(run_members[0])
        elif layout == "paired_end":
            if (
                len(run_members) != 2
                or roles != ["paired_end_mate_1", "paired_end_mate_2"]
                or indices != [0, 1]
                or counts != {2}
            ):
                raise GSE268273QuantificationError(
                    "paired-end technical-run topology differs"
                )
            upstream.append(run_members[0])
            downstream.append(run_members[1])
        else:
            raise GSE268273QuantificationError("effective read layout differs")
    if layout == "paired_end" and len(upstream) != len(downstream):
        raise GSE268273QuantificationError("paired-end mate cardinality differs")
    return layout, upstream, downstream


def _load_plan_bundle(plan_root: Path, bundle_id: int) -> list[dict[str, str]]:
    verify_frozen_tree(plan_root)
    fields, rows = read_tsv(plan_root / "fastq_files.tsv")
    required = {
        "bundle_id",
        "row_id",
        "run_accession",
        "declared_ena_library_layout",
        "effective_library_layout",
        "file_part_index",
        "file_part_count",
        "read_role",
        "relative_path",
        "fastq_bytes",
    }
    forbidden = {"fibrosis", "outcome", "label", "disease", "sex", "nas", "bmi"}
    if not required <= set(fields) or forbidden & set(fields):
        raise GSE268273QuantificationError("quantification plan violates label firewall")
    selected = [row for row in rows if int(row["bundle_id"]) == bundle_id]
    if not selected:
        raise GSE268273QuantificationError("quantification bundle is absent")
    return selected


def quantify_bundle(
    *,
    plan_root: Path,
    raw_state_root: Path,
    reference_root: Path,
    model_input_root: Path,
    bundle_id: int,
    output: Path,
    attempt_root: Path,
    threads: int,
    seed: int,
) -> dict[str, Any]:
    if (output / "COMPLETE").exists():
        verify_frozen_tree(output)
        return json.loads((output / "bundle_receipt.json").read_text(encoding="utf-8"))
    if output.exists() and (output / "ARTIFACTS.json").exists():
        raise GSE268273QuantificationError("incomplete bundle is already frozen")
    rows = _load_plan_bundle(plan_root, bundle_id)
    planned_layout_sets: defaultdict[str, set[str]] = defaultdict(set)
    for row in rows:
        planned_layout_sets[row["row_id"]].add(row["effective_library_layout"])
    if any(len(layouts) != 1 for layouts in planned_layout_sets.values()):
        raise GSE268273QuantificationError("planned participant layout is mixed")
    planned_layout_participants = Counter(
        next(iter(layouts)) for layouts in planned_layout_sets.values()
    )
    raw_bundle = raw_state_root / f"bundle_{bundle_id:02d}"
    if not (raw_bundle / "COMPLETE").is_file():
        raise GSE268273QuantificationError("verified raw bundle is incomplete")
    verify_frozen_tree(raw_bundle)
    verify_frozen_tree(reference_root)
    raw_receipt = json.loads(
        (raw_bundle / "bundle_receipt.json").read_text(encoding="utf-8")
    )
    if (
        raw_receipt.get("status") != "complete_verified"
        or raw_receipt.get("bundle_id") != bundle_id
        or raw_receipt.get("fastq_files") != len(rows)
        or raw_receipt.get("fastq_bytes")
        != sum(int(row["fastq_bytes"]) for row in rows)
        or raw_receipt.get("effective_layout_participants")
        != dict(sorted(planned_layout_participants.items()))
        or raw_receipt.get("ena_md5_verified_for_every_file") is not True
        or raw_receipt.get("gzip_integrity_verified_for_every_file") is not True
        or raw_receipt.get("fastq_first_record_layout_verified_for_every_file")
        is not True
        or raw_receipt.get("paired_first_record_identity_verified_for_every_run")
        is not True
        or raw_receipt.get("labels_accessed") is not False
    ):
        raise GSE268273QuantificationError("verified raw-bundle receipt differs")
    reference_receipt = json.loads(
        (reference_root / "receipt.json").read_text(encoding="utf-8")
    )
    if (
        not (reference_root / "COMPLETE").is_file()
        or reference_receipt.get("schema_version")
        != "masld-bench-gse268273-rsem-reference-build-v2"
        or reference_receipt.get("status") != "passed"
        or reference_receipt.get("source_annotation")
        != "GENCODE_v36_primary_assembly"
        or reference_receipt.get("reference_contract_revision")
        != "gse268273_reference_primary_assembly_revision_v2"
        or reference_receipt.get("source_annotation_gene_ids") != 60_719
        or reference_receipt.get("source_gtf_contigs") != 47
        or reference_receipt.get("source_features_admitted") != 14_149
        or reference_receipt.get("source_primary_scaffold_features_admitted")
        != 7
        or reference_receipt.get("target_stable_genes_admitted") != 14_078
        or reference_receipt.get("target_duplicate_sum_groups") != 11
        or "RSEM v1.3.1" not in str(reference_receipt.get("rsem_version", ""))
        or "2.7.10b" not in str(reference_receipt.get("star_version", ""))
    ):
        raise GSE268273QuantificationError("frozen RSEM reference contract differs")
    verify_frozen_tree(model_input_root)
    candidate_fields, candidates = read_tsv(
        model_input_root / "candidate_gene_axis.tsv"
    )
    if (
        len(candidates) != 14_078
        or "stable_gene_id" not in candidate_fields
        or sum(int(row["source_feature_count"]) > 1 for row in candidates) != 11
    ):
        raise GSE268273QuantificationError("candidate v49 axis differs")
    tools = tool_environment()
    rsem_calculate = tool_executable("rsem-calculate-expression", tools)
    star_executable = tool_executable("STAR", tools)
    output.mkdir(parents=True, exist_ok=True)
    participants_root = output / "participants"
    participants_root.mkdir(exist_ok=True)
    attempt_root.mkdir(parents=True)
    by_participant: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_participant[row["row_id"]].append(row)
    completed = []
    completed_layouts: Counter[str] = Counter()
    reference_prefix = reference_root / "reference" / "gse268273_gencode_v36"
    star_path = str(Path(star_executable).parent)
    for row_id in sorted(by_participant):
        final = participants_root / row_id
        if (final / "COMPLETE").exists():
            verify_frozen_tree(final)
            existing_receipt = json.loads(
                (final / "receipt.json").read_text(encoding="utf-8")
            )
            completed_layouts[str(existing_receipt["effective_library_layout"])] += 1
            completed.append(row_id)
            continue
        stage = attempt_root / f"{row_id}.staging"
        if stage.exists():
            raise GSE268273QuantificationError(
                f"attempt-local participant stage already exists: {stage}"
            )
        stage.mkdir()
        members = by_participant[row_id]
        effective_layout, upstream_members, downstream_members = (
            arrange_participant_fastqs(members)
        )
        paths_by_identity: dict[tuple[str, int], str] = {}
        for row in sorted(
            members,
            key=lambda value: (
                value["run_accession"],
                int(value["file_part_index"]),
            ),
        ):
            relative = PurePosixPath(row["relative_path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise GSE268273QuantificationError("raw FASTQ relative path is unsafe")
            path = raw_bundle.joinpath(*relative.parts)
            if not path.is_file() or path.stat().st_size != int(row["fastq_bytes"]):
                raise GSE268273QuantificationError("verified FASTQ member is absent")
            paths_by_identity[(row["run_accession"], int(row["file_part_index"]))] = str(
                path
            )
        upstream_paths = [
            paths_by_identity[(row["run_accession"], int(row["file_part_index"]))]
            for row in upstream_members
        ]
        downstream_paths = [
            paths_by_identity[(row["run_accession"], int(row["file_part_index"]))]
            for row in downstream_members
        ]
        prefix = stage / "rsem"
        command = [
            rsem_calculate,
            "--star",
            "--star-path",
            star_path,
            "--star-gzipped-read-file",
            "--num-threads",
            str(threads),
            "--no-bam-output",
            "--strandedness",
            "none",
            "--seed",
            str(seed),
        ]
        if effective_layout == "paired_end":
            command.append("--paired-end")
        command.append(",".join(upstream_paths))
        if effective_layout == "paired_end":
            command.append(",".join(downstream_paths))
        command.extend([str(reference_prefix), str(prefix)])
        with (stage / "command.stdout.txt").open("x", encoding="utf-8") as stdout, (
            stage / "command.stderr.txt"
        ).open("x", encoding="utf-8") as stderr:
            result = subprocess.run(
                command,
                check=False,
                stdout=stdout,
                stderr=stderr,
                text=True,
                env=tools,
            )
        if result.returncode != 0:
            raise GSE268273QuantificationError(f"RSEM failed for {row_id}")
        result_path = prefix.with_suffix(".genes.results")
        rsem_counts = parse_rsem_gene_results(result_path)
        values = aggregate_v49_expected_counts(rsem_counts, candidates)
        np.save(stage / "v49_expected_counts.npy", values, allow_pickle=False)
        retained = stage / "rsem.genes.results"
        if result_path != retained:
            result_path.rename(retained)
        receipt = {
            "schema_version": "masld-bench-gse268273-participant-rsem-v1",
            "status": "passed_raw_count_scale",
            "row_id": row_id,
            "technical_runs": len({row["run_accession"] for row in members}),
            "fastq_files": len(members),
            "fastq_bytes": sum(int(row["fastq_bytes"]) for row in members),
            "declared_ena_library_layout": "SINGLE",
            "effective_library_layout": effective_layout,
            "input_layout": (
                "paired_end_matewise_comma_ordered_technical_runs"
                if effective_layout == "paired_end"
                else "single_end_comma_ordered_technical_runs"
            ),
            "declared_layout_discrepancy": effective_layout == "paired_end",
            "biological_library": "one_participant",
            "rsem_source_genes": len(rsem_counts),
            "v49_target_genes": len(values),
            "v49_duplicate_sum_groups": 11,
            "measurement": "RSEM_expected_count_raw_count_scale",
            "strandedness": "none_source_default_not_deposited",
            "alignment_parameters": "RSEM_integrated_STAR_defaults",
            "single_end_fragment_distribution_parameters": (
                "RSEM_defaults_source_did_not_deposit_values"
                if effective_layout == "single_end"
                else "learned_from_paired_reads"
            ),
            "normalization_applied": False,
            "total_expected_count": float(values.sum()),
            "seed": seed,
            "labels_accessed": False,
            "fit_or_score_performed": False,
            "plan_artifacts_sha256": sha256_file(plan_root / "ARTIFACTS.json"),
            "raw_bundle_artifacts_sha256": sha256_file(
                raw_bundle / "ARTIFACTS.json"
            ),
            "reference_artifacts_sha256": sha256_file(
                reference_root / "ARTIFACTS.json"
            ),
        }
        (stage / "receipt.json").write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        freeze_tree(
            stage,
            {
                "artifact_class": "gse268273_participant_rsem_expected_counts",
                "row_id": row_id,
                "target_genes": 14_078,
                "status": "passed",
            },
        )
        verify_frozen_tree(stage)
        if final.exists():
            raise GSE268273QuantificationError(
                f"refusing to overwrite completed participant: {final}"
            )
        stage.rename(final)
        completed_layouts[effective_layout] += 1
        completed.append(row_id)
    expected_layouts = Counter(
        next(iter({row["effective_library_layout"] for row in members}))
        for members in by_participant.values()
    )
    if set(completed) != set(by_participant) or completed_layouts != expected_layouts:
        raise GSE268273QuantificationError(
            "completed participant or effective-layout census differs"
        )
    receipt = {
        "schema_version": "masld-bench-gse268273-rsem-bundle-v1",
        "status": "complete",
        "bundle_id": bundle_id,
        "participants": len(completed),
        "effective_layout_participants": dict(sorted(completed_layouts.items())),
        "target_genes": 14_078,
        "target_duplicate_sum_groups": 11,
        "measurement": "RSEM_expected_count_raw_count_scale",
        "strandedness": "none_source_default_not_deposited",
        "alignment_parameters": "RSEM_integrated_STAR_defaults",
        "normalization_applied": False,
        "labels_accessed": False,
        "fit_or_score_performed": False,
        "participant_level_checkpointing": True,
        "mid_participant_resume": False,
    }
    (output / "bundle_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(
        output,
        {
            "artifact_class": "gse268273_rsem_expected_count_bundle",
            "bundle_id": bundle_id,
            "participants": len(completed),
            "status": "passed",
        },
    )
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--raw-state-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--model-input-root", type=Path, required=True)
    parser.add_argument("--bundle-id", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--attempt-root", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--seed", type=int, default=268273)
    args = parser.parse_args()
    receipt = quantify_bundle(
        plan_root=args.plan_root,
        raw_state_root=args.raw_state_root,
        reference_root=args.reference_root,
        model_input_root=args.model_input_root,
        bundle_id=args.bundle_id,
        output=args.output,
        attempt_root=args.attempt_root,
        threads=args.threads,
        seed=args.seed,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

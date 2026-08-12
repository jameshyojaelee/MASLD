#!/usr/bin/env python3
"""Verify and snapshot the accepted F_five and corrected nine-cohort substrate."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

from accepted_run_provenance import (
    EXPECTED_GRAPH_BYTES,
    EXPECTED_GRAPH_FILES,
    Selection,
    collect_accepted_run_selections,
    verify_accepted_graph,
)
from resource_contract import (
    BG_ROOT,
    CANONICAL_REFERENCE,
    CANONICAL_REFERENCE_SHA256,
    CANDIDATE_ID,
    CANDIDATE_ROOT,
    EXPECTED_HASHES,
    EXPECTED_POOLED,
    EXPECTED_RESOURCE,
    F_FIVE_ROOT,
    F_LEGACY_ROOT,
    GENCODE_GTF,
    GENCODE_GTF_BYTES,
    GENCODE_GTF_SHA256,
    GSE193066_METADATA,
    GSE193066_METADATA_SHA256,
    PROJECT_ROOT,
    REVIEWED_CODE_FILES,
    RNASEQ_PYTHON,
    RNASEQ_PYTHON_SHA256,
    UNIFIED_METADATA,
    UNIFIED_METADATA_SHA256,
    ContractError,
    read_tsv,
    require,
    require_regular_file,
    sha256,
    verify_manifest,
)
from snapshot_io import (
    SnapshotIOError,
    copy_exclusive,
    parse_sha256_manifest,
    publish_directory_noreplace,
)


SCRIPT_DIR = Path(__file__).resolve().parent
REVIEWED_MANIFEST = SCRIPT_DIR / "reviewed_code.sha256"


def write_exclusive(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o440)
    with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def tsv_payload(rows: list[dict[str, object]], fields: list[str]) -> str:
    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def verify_expected_hashes() -> None:
    for relative, expected in EXPECTED_HASHES.items():
        source = BG_ROOT / relative
        require_regular_file(source)
        require(sha256(source) == expected, f"expected hash mismatch: {relative}")


def verify_deg_counts() -> None:
    path = F_FIVE_ROOT / "results/integration/deg_results.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"gene", "logFC", "treat_lfc", "treat_fdr"}
        require(required.issubset(reader.fieldnames or []), "F_five DEG schema is incomplete")
        genes: set[str] = set()
        n_treat = n_up = n_down = 0
        for row in reader:
            require(None not in row, "malformed F_five DEG row width")
            gene = row["gene"]
            require(gene and gene not in genes, f"blank or duplicate DEG gene: {gene!r}")
            genes.add(gene)
            require(float(row["treat_lfc"]) == 0.25, "TREAT lfc contract drift")
            if float(row["treat_fdr"]) < 0.05:
                n_treat += 1
                effect = float(row["logFC"])
                require(effect != 0, "significant gene has zero effect")
                n_up += int(effect > 0)
                n_down += int(effect < 0)
    observed = {
        "n_genes": len(genes),
        "n_treat": n_treat,
        "n_up": n_up,
        "n_down": n_down,
    }
    for key, value in observed.items():
        require(value == EXPECTED_POOLED[key], f"F_five {key} drift: {value}")


def verify_sample_censuses() -> None:
    model = read_tsv(F_FIVE_ROOT / "results/integration/model_design.tsv")
    require(len(model) == EXPECTED_POOLED["n_samples"], "pooled model sample-count drift")
    sample_ids = [row["sample_id"] for row in model]
    require(len(set(sample_ids)) == len(sample_ids), "duplicate pooled model sample")
    disease = sum(float(row["group_binaryDisease"]) == 1.0 for row in model)
    require(disease == EXPECTED_POOLED["n_disease"], "pooled disease count drift")
    require(
        len(model) - disease == EXPECTED_POOLED["n_control"],
        "pooled control count drift",
    )

    with (F_LEGACY_ROOT / "qc/sample_qc_report.csv").open(
        newline="",
        encoding="utf-8",
    ) as handle:
        qc = list(csv.DictReader(handle))
    require(len(qc) == EXPECTED_RESOURCE["n_matched"], "nine-cohort QC row-count drift")
    require(len({row["sample_id"] for row in qc}) == len(qc), "duplicate QC sample")
    n_pass = sum(row["pass_technical"].strip().upper() == "TRUE" for row in qc)
    require(n_pass == EXPECTED_RESOURCE["n_pass_technical"], "pass_technical drift")

    pre = read_tsv(F_LEGACY_ROOT / "provenance/preprocessing_manifest.tsv")
    require(len(pre) == 1, "F_legacy preprocessing manifest cardinality drift")
    require(
        int(pre[0]["n_samples"]) == EXPECTED_RESOURCE["n_pass_technical"],
        "F_legacy DGE sample drift",
    )
    require(
        int(pre[0]["n_genes"]) == EXPECTED_RESOURCE["n_f_legacy_genes"],
        "F_legacy gene drift",
    )


def rederive_bulk_object_census() -> dict[str, int]:
    require(
        os.environ.get("SLURM_JOB_PARTITION") == "io"
        and os.environ.get("SLURM_JOB_ID", "").isdigit(),
        "bulk RDS census must run in the ffreeze I/O job",
    )
    runtime_path = BG_ROOT / "contract/analysis_runtime_contract.json"
    require_regular_file(runtime_path)
    runtime = json.loads(runtime_path.read_text(encoding="utf-8"))
    require(
        isinstance(runtime, dict) and isinstance(runtime.get("Rscript"), dict),
        "runtime contract drift",
    )
    rscript = Path(str(runtime["Rscript"]["path"]))
    require_regular_file(rscript)
    require(sha256(rscript) == runtime["Rscript"]["sha256"], "Rscript identity drift")
    script = SCRIPT_DIR / "preflight_bulk_census.R"
    require_regular_file(script)
    result = subprocess.run(
        [
            str(rscript),
            "--vanilla",
            str(script),
            str(F_LEGACY_ROOT / "results/integration/merged_counts_raw.rds"),
            str(F_LEGACY_ROOT / "results/integration/merged_dge.rds"),
            str(F_FIVE_ROOT / "results/integration/merged_dge.rds"),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=3600,
    )
    require(
        result.returncode == 0,
        f"bulk RDS census subprocess failed: {result.stderr.strip()}",
    )
    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    require(reader.fieldnames == ["metric", "value"], "bulk RDS census schema drift")
    rows = list(reader)
    require(len(rows) == 8, "bulk RDS census cardinality drift")
    require(len({row["metric"] for row in rows}) == len(rows), "duplicate census metric")
    observed = {row["metric"]: int(row["value"]) for row in rows}
    expected = {
        "raw_n_genes": EXPECTED_RESOURCE["n_raw_genes"],
        "raw_n_samples": EXPECTED_RESOURCE["n_matched"],
        "legacy_n_genes": EXPECTED_RESOURCE["n_f_legacy_genes"],
        "legacy_n_samples": EXPECTED_RESOURCE["n_pass_technical"],
        "legacy_n_cohorts": EXPECTED_RESOURCE["n_cohorts"],
        "five_n_genes": EXPECTED_POOLED["n_genes"],
        "five_n_samples": EXPECTED_POOLED["n_samples"],
        "five_n_cohorts": EXPECTED_POOLED["n_cohorts"],
    }
    require(observed == expected, f"bulk object census drift: {observed}")
    return observed


def verify_gse193066_metadata() -> None:
    require_regular_file(GSE193066_METADATA)
    require(
        sha256(GSE193066_METADATA) == GSE193066_METADATA_SHA256,
        "GSE193066 metadata identity drift",
    )
    rows = read_tsv(GSE193066_METADATA)
    require(len(rows) == 164, "GSE193066 biopsy-row count drift")
    required = {"sample_id", "!Sample_title", "biopsy", "fibrosis stage"}
    require(required.issubset(rows[0]), "GSE193066 metadata schema drift")
    require(len({row["sample_id"] for row in rows}) == 164, "duplicate GSE193066 sample ID")
    require(len({row["!Sample_title"] for row in rows}) == 164, "duplicate GSE193066 title")
    with (F_LEGACY_ROOT / "qc/sample_qc_report.csv").open(
        newline="",
        encoding="utf-8",
    ) as handle:
        qc_rows = list(csv.DictReader(handle))
    qc_by_sample = {row["sample_id"]: row for row in qc_rows}
    require(
        all(row["sample_id"] in qc_by_sample for row in rows),
        "GSE193066 sample missing from corrected QC",
    )
    require(
        sum(
            qc_by_sample[row["sample_id"]]["pass_technical"].upper() == "TRUE"
            for row in rows
        )
        == 160,
        "GSE193066 corrected-QC pass count drift",
    )
    by_participant: dict[str, list[str]] = {}
    for row in rows:
        title = row["!Sample_title"]
        key = title[:-2] if title.endswith(("_1", "_2")) else title
        by_participant.setdefault(key, []).append(row["biopsy"])
        if row["biopsy"] == "2nd biopsy":
            require(title.endswith("_2"), "second-biopsy title suffix drift")
        else:
            require(row["biopsy"] == "1st biopsy", "unexpected GSE193066 biopsy timing")
            require(
                title.endswith("_1") or not title.endswith(("_1", "_2")),
                "first-biopsy title suffix drift",
            )
        require(
            row["fibrosis stage"] in {"0", "1", "2", "3", "4"},
            "GSE193066 fibrosis-stage schema drift",
        )
    require(len(by_participant) == 106, "GSE193066 participant-key count drift")
    paired = [
        key
        for key, biopsies in by_participant.items()
        if sorted(biopsies) == ["1st biopsy", "2nd biopsy"]
    ]
    require(len(paired) == 58, "GSE193066 paired-participant count drift")
    require(
        all(
            biopsies == ["1st biopsy"]
            or sorted(biopsies) == ["1st biopsy", "2nd biopsy"]
            for biopsies in by_participant.values()
        ),
        "GSE193066 participant has duplicate or invalid biopsy timing",
    )


def verify_reviewed_code_manifest() -> tuple[str, dict[str, str]]:
    require_regular_file(REVIEWED_MANIFEST)
    expected_manifest_hash = os.environ.get(
        "MASLD_REVIEWED_CODE_MANIFEST_SHA256",
        "",
    )
    require(
        len(expected_manifest_hash) == 64
        and all(character in "0123456789abcdef" for character in expected_manifest_hash),
        "missing or malformed reviewed-code manifest hash",
    )
    require(
        sha256(REVIEWED_MANIFEST) == expected_manifest_hash,
        "reviewed-code manifest identity drift",
    )
    rows = parse_sha256_manifest(
        REVIEWED_MANIFEST,
        separator="  ",
        expected_rows=len(REVIEWED_CODE_FILES),
    )
    expected_paths = {
        (SCRIPT_DIR / name).relative_to(PROJECT_ROOT).as_posix()
        for name in REVIEWED_CODE_FILES
    }
    observed_paths = {raw for _, raw in rows}
    require(observed_paths == expected_paths, "reviewed-code file-set drift")
    by_path: dict[str, str] = {}
    for digest, raw in rows:
        relative = PurePosixPath(raw)
        require(
            not relative.is_absolute() and ".." not in relative.parts,
            f"unsafe reviewed-code path: {raw}",
        )
        source = PROJECT_ROOT.joinpath(*relative.parts)
        require_regular_file(source)
        require(sha256(source) == digest, f"reviewed-code hash drift: {raw}")
        by_path[raw] = digest
    require(
        Path(sys.executable).resolve(strict=True) == RNASEQ_PYTHON,
        "ffreeze Python executable path drift",
    )
    require(sha256(RNASEQ_PYTHON) == RNASEQ_PYTHON_SHA256, "ffreeze Python hash drift")
    return expected_manifest_hash, by_path


def add_selection(
    selections: dict[str, Selection],
    source: Path,
    destination: str,
    role: str,
    *,
    expected_size: int | None = None,
    expected_hash: str | None = None,
) -> None:
    relative = PurePosixPath(destination)
    require(
        not relative.is_absolute() and ".." not in relative.parts,
        f"unsafe snapshot destination: {destination}",
    )
    require_regular_file(source)
    size = source.stat().st_size
    digest = sha256(source)
    if expected_size is not None:
        require(size == expected_size, f"source size drift: {source}")
    if expected_hash is not None:
        require(digest == expected_hash, f"source hash drift: {source}")
    if destination in selections:
        old_source, old_role, old_size, old_hash = selections[destination]
        require(
            old_source == source
            and old_size == size
            and old_hash == digest,
            f"snapshot destination collision: {destination}",
        )
        roles = sorted(set(old_role.split(";")) | {role})
        selections[destination] = (source, ";".join(roles), size, digest)
        return
    selections[destination] = (source, role, size, digest)


def collect_selections(
    reviewed_manifest_hash: str,
    reviewed_rows: dict[str, str],
) -> dict[str, Selection]:
    selections = collect_accepted_run_selections()

    runtime_source = (
        BG_ROOT
        / "source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
    )
    for name in ("analysis_runtime_contract.py", "safe_io.py"):
        add_selection(
            selections,
            runtime_source / name,
            f"code/frozen_bg001_runtime/{name}",
            "frozen_bg001_runtime_verifier",
        )

    audit_root = PROJECT_ROOT / "results/code_review/2026-08-04_full_codebase"
    for relative in (
        "BG-001_ESCALATION_RESPONSE_2026-08-09.md",
        "BG-001_FOUR_ARM_RESULTS_2026-08-09.md",
        "BG-001_HEADLINE_IMPACT.md",
        "BG-001_analysis_stage_deviations_2026-08-09.md",
        "bg001_escalation_refits/escalation_refit_comparisons.tsv",
        "bg001_escalation_refits/escalation_refit_summary.tsv",
    ):
        add_selection(
            selections,
            audit_root / relative,
            f"inputs/BG001-DECISION/external_review/{relative}",
            "bg001_external_acceptance_support",
        )

    add_selection(
        selections,
        CANONICAL_REFERENCE,
        "inputs/CANONICAL-REFERENCE/canonical_deg_results.csv",
        "unchanged_superseded_reference",
        expected_hash=CANONICAL_REFERENCE_SHA256,
    )
    add_selection(
        selections,
        GSE193066_METADATA,
        "inputs/BULK-NINE-COHORT/source_metadata/GSE193066_metadata.tsv",
        "authoritative_first_second_biopsy_key",
        expected_hash=GSE193066_METADATA_SHA256,
    )

    for source, destination, role, digest in (
        (
            UNIFIED_METADATA,
            "inputs/BG001-DECISION/source_snapshot/RNA-seq/Human/Patient_Cohorts/"
            "analysis/integration/metadata/unified_metadata.csv",
            "corrected_stage_metadata",
            UNIFIED_METADATA_SHA256,
        ),
        (
            BG_ROOT / "source_snapshot/config/human_datasets.yaml",
            "inputs/BULK-F-FIVE/frozen_model_inputs/human_datasets.yaml",
            "pooled_config",
            None,
        ),
        (
            BG_ROOT / "source_snapshot/data/gencode_v49_gene_metadata.tsv.gz",
            "inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz",
            "gene_annotation",
            None,
        ),
        (
            PROJECT_ROOT
            / "docs/archive/documentation_consolidation_2026-08-11/originals/docs/manuscript/working/F_FIVE_RESOURCE_ADOPTION.md",
            "inputs/BG001-DECISION/F_FIVE_RESOURCE_ADOPTION.md",
            "owner_adoption",
            None,
        ),
    ):
        add_selection(
            selections,
            source,
            destination,
            role,
            expected_hash=digest,
        )

    add_selection(
        selections,
        GENCODE_GTF,
        "inputs/BG001-DECISION/external_reference/"
        "gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz",
        "external_reference_gtf",
        expected_size=GENCODE_GTF_BYTES,
        expected_hash=GENCODE_GTF_SHA256,
    )

    for name in REVIEWED_CODE_FILES:
        source = SCRIPT_DIR / name
        raw = source.relative_to(PROJECT_ROOT).as_posix()
        add_selection(
            selections,
            source,
            f"code/{name}",
            "resource_candidate_code",
            expected_hash=reviewed_rows[raw],
        )
    add_selection(
        selections,
        REVIEWED_MANIFEST,
        "code/reviewed_code.sha256",
        "reviewed_code_manifest",
        expected_hash=reviewed_manifest_hash,
    )
    return selections


def source_provenance(source: Path) -> tuple[str, str]:
    resolved = source.resolve(strict=True)
    try:
        return "project_relative", resolved.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return "external_absolute", resolved.as_posix()


def preflight() -> tuple[dict[str, Selection], dict[str, int], str]:
    require(PROJECT_ROOT.is_dir(), "project root is missing")
    require(BG_ROOT.is_dir(), "accepted BG-001 run is missing")
    require(
        not CANDIDATE_ROOT.exists() and not CANDIDATE_ROOT.is_symlink(),
        "candidate root already exists",
    )
    reviewed_manifest_hash, reviewed_rows = verify_reviewed_code_manifest()
    require_regular_file(CANONICAL_REFERENCE)
    require(
        sha256(CANONICAL_REFERENCE) == CANONICAL_REFERENCE_SHA256,
        "live canonical reference drift",
    )
    require_regular_file(UNIFIED_METADATA)
    require(
        sha256(UNIFIED_METADATA) == UNIFIED_METADATA_SHA256,
        "accepted unified metadata drift",
    )
    verify_expected_hashes()
    verify_manifest(F_FIVE_ROOT, F_FIVE_ROOT / "provenance/artifact_manifest.tsv")
    verify_manifest(F_LEGACY_ROOT, F_LEGACY_ROOT / "provenance/artifact_manifest.tsv")
    verify_manifest(BG_ROOT / "comparisons", BG_ROOT / "comparisons/artifact_manifest.tsv")
    verify_deg_counts()
    verify_sample_censuses()
    verify_gse193066_metadata()
    selections = collect_selections(reviewed_manifest_hash, reviewed_rows)
    require(selections, "source selection is empty")
    return selections, rederive_bulk_object_census(), reviewed_manifest_hash


def prepare() -> None:
    selections, object_census, reviewed_manifest_hash = preflight()
    parent = CANDIDATE_ROOT.parent
    parent.mkdir(parents=True, exist_ok=True)
    temporary = parent / f".{CANDIDATE_ID}.tmp.{os.getpid()}.{uuid.uuid4().hex}"
    require(
        not temporary.exists() and not temporary.is_symlink(),
        "temporary root collision",
    )
    require(temporary.mkdir(mode=0o750) is None, "failed to create temporary root")

    try:
        manifest_rows: list[dict[str, object]] = []
        for destination, (source, role, expected_size, expected_hash) in sorted(
            selections.items()
        ):
            target = temporary / destination
            copied_size, copied_hash = copy_exclusive(source, target)
            require(
                copied_size == expected_size and copied_hash == expected_hash,
                f"descriptor copy differs from sealed identity: {source}",
            )
            require(
                target.stat().st_size == expected_size
                and sha256(target) == expected_hash,
                f"candidate copy differs from sealed identity: {destination}",
            )
            scope, source_path = source_provenance(source)
            manifest_rows.append(
                {
                    "source_scope": scope,
                    "source_path": source_path,
                    "candidate_path": destination,
                    "role": role,
                    "size_bytes": expected_size,
                    "sha256": expected_hash,
                }
            )

        mirrored_graph = verify_accepted_graph(
            temporary / "inputs/BG001-DECISION"
        )
        require(len(mirrored_graph) == EXPECTED_GRAPH_FILES, "temporary graph count drift")
        require(
            sum(size for _, size, _ in mirrored_graph.values())
            == EXPECTED_GRAPH_BYTES,
            "temporary graph byte drift",
        )

        source_manifest = temporary / "manifests/source_selection.tsv"
        write_exclusive(
            source_manifest,
            tsv_payload(
                manifest_rows,
                [
                    "source_scope",
                    "source_path",
                    "candidate_path",
                    "role",
                    "size_bytes",
                    "sha256",
                ],
            ),
        )
        object_census_rows = [
            {"metric": metric, "value": value}
            for metric, value in object_census.items()
        ]
        object_census_path = temporary / "manifests/bulk_object_census.tsv"
        write_exclusive(
            object_census_path,
            tsv_payload(object_census_rows, ["metric", "value"]),
        )
        census_rows = [
            {
                "layer": "matched_raw_nine_cohort",
                "n_samples": object_census["raw_n_samples"],
                "n_genes": object_census["raw_n_genes"],
                "role": "resource_source",
            },
            {
                "layer": "pass_technical_nine_cohort",
                "n_samples": object_census["legacy_n_samples"],
                "n_genes": object_census["legacy_n_genes"],
                "role": "stage_context_dge",
            },
            {
                "layer": "pooled_five_cohort",
                "n_samples": object_census["five_n_samples"],
                "n_genes": object_census["five_n_genes"],
                "role": "disease_control_model",
            },
        ]
        write_exclusive(
            temporary / "manifests/sample_census.tsv",
            tsv_payload(census_rows, ["layer", "n_samples", "n_genes", "role"]),
        )
        contract = {
            "schema": "masld-resource-f-five-candidate-v1",
            "candidate_id": CANDIDATE_ID,
            "bg001_run_id": BG_ROOT.name,
            "bg001_graph_file_count": EXPECTED_GRAPH_FILES,
            "bg001_graph_size_bytes": EXPECTED_GRAPH_BYTES,
            "approved_pooled_layer": "human_bulk_disease_control_fivecohort_candidate",
            "resource_source_layer": "human_bulk_corrected_ninecohort_candidate",
            "canonical_promotion_authorized": False,
            "figure_main_write_authorized": False,
            "portal_deployment_authorized": False,
            "corrected_coloc_status": "pending",
            "pooled_reproduction_required": True,
            "source_manifest_sha256": sha256(source_manifest),
            "bulk_object_census_sha256": sha256(object_census_path),
            "reviewed_code_manifest_sha256": reviewed_manifest_hash,
            "gse193066_metadata_sha256": GSE193066_METADATA_SHA256,
            "external_gtf_sha256": GENCODE_GTF_SHA256,
            "external_gtf_size_bytes": GENCODE_GTF_BYTES,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        write_exclusive(
            temporary / "manifests/adoption_contract.json",
            json.dumps(contract, indent=2, sort_keys=True) + "\n",
        )
        for directory in (
            temporary / "workstreams",
            temporary / "logs",
            temporary / "manifests/jobs",
        ):
            directory.mkdir(parents=True, exist_ok=True, mode=0o750)
        write_exclusive(temporary / ".resource_candidate_root", CANDIDATE_ID + "\n")
        sentinel = temporary / ".resource_candidate_root"
        base_rows = []
        for path in sorted(item for item in temporary.rglob("*") if item.is_file()):
            relative = path.relative_to(temporary).as_posix()
            if relative in {
                "BASE_SNAPSHOT_COMPLETE.json",
                "manifests/base_artifact_manifest.tsv",
            }:
                continue
            base_rows.append(
                {
                    "candidate_path": relative,
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
        base_manifest = temporary / "manifests/base_artifact_manifest.tsv"
        write_exclusive(
            base_manifest,
            tsv_payload(base_rows, ["candidate_path", "size_bytes", "sha256"]),
        )
        completion = {
            "status": "BASE_SNAPSHOT_COMPLETE",
            "candidate_id": CANDIDATE_ID,
            "source_file_count": len(manifest_rows),
            "source_manifest_sha256": sha256(source_manifest),
            "base_file_count": len(base_rows),
            "base_artifact_manifest_sha256": sha256(base_manifest),
            "sample_census_sha256": sha256(
                temporary / "manifests/sample_census.tsv"
            ),
            "bulk_object_census_sha256": sha256(object_census_path),
            "reviewed_code_manifest_sha256": reviewed_manifest_hash,
            "adoption_contract_sha256": sha256(
                temporary / "manifests/adoption_contract.json"
            ),
            "sentinel_sha256": sha256(sentinel),
            "canonical_promotion_authorized": False,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        write_exclusive(
            temporary / "BASE_SNAPSHOT_COMPLETE.json",
            json.dumps(completion, indent=2, sort_keys=True) + "\n",
        )
        require(
            not CANDIDATE_ROOT.exists() and not CANDIDATE_ROOT.is_symlink(),
            "candidate root appeared during preparation",
        )
        publication = publish_directory_noreplace(
            temporary,
            CANDIDATE_ROOT,
            marker="BASE_SNAPSHOT_COMPLETE.json",
        )
        require(
            publication.get("source_retained") is True,
            "candidate publisher did not report retained source",
        )
        require(
            temporary.is_dir() and not temporary.is_symlink(),
            "candidate publication source was not retained safely",
        )
        require(
            CANDIDATE_ROOT.is_dir() and not CANDIDATE_ROOT.is_symlink(),
            "published candidate root is invalid",
        )
        source_marker = temporary / "BASE_SNAPSHOT_COMPLETE.json"
        official_marker = CANDIDATE_ROOT / "BASE_SNAPSHOT_COMPLETE.json"
        require_regular_file(source_marker)
        require_regular_file(official_marker)
        require(
            sha256(source_marker) == sha256(official_marker),
            "published base completion-marker hash drift",
        )
    except Exception:
        print(
            f"FAILED temporary candidate retained for audit: {temporary}",
            file=sys.stderr,
        )
        raise

    print(json.dumps(completion, indent=2, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight-only", action="store_true")
    arguments = parser.parse_args()
    try:
        if arguments.preflight_only:
            selected, object_census, manifest_hash = preflight()
            print(
                json.dumps(
                    {
                        "status": "PREFLIGHT_PASS",
                        "candidate_id": CANDIDATE_ID,
                        "source_file_count": len(selected),
                        "bulk_object_census": object_census,
                        "reviewed_code_manifest_sha256": manifest_hash,
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
        else:
            prepare()
    except (ContractError, SnapshotIOError) as error:
        raise SystemExit(f"CONTRACT ERROR: {error}") from error

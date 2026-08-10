#!/usr/bin/env python3
"""Adjudicate Plan 44 without reopening its frozen hypotheses.

The two required outcome branches have already run.  This script copies their
immutable source tables into the documented candidate contract, records the
mechanical no-promotion verdict, and fails closed on any upstream byte drift.
It does not search for a replacement cue, lineage, cohort, program, or gene.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import platform
import shutil
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT_ROOT = Path(__file__).resolve().parent
CANDIDATE = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates/"
    "multicellular-assembly-response-2026-08-09"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    temporary.replace(path)


def truth(value: str) -> bool:
    return str(value).strip().lower() == "true"


def atomic_json(path: Path, value: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def copy_contract(source: Path, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    shutil.copyfile(source, temporary)
    temporary.replace(destination)


def main() -> None:
    preflight = json.loads((CANDIDATE / "PREFLIGHT_VALIDATED.json").read_text())
    if preflight["status"] != "PREFLIGHT_VALIDATED_OUTCOMES_UNOPENED":
        raise RuntimeError("validated outcome-blind preflight is required")
    if sha256(CANDIDATE / "release_manifest.tsv") != preflight["release_manifest_sha256"]:
        raise RuntimeError("preflight release manifest drifted after outcome access")

    for row in read_tsv(CANDIDATE / "upstream_immutability_baseline.tsv"):
        source = ROOT / row["path"]
        if (
            not source.is_file()
            or source.stat().st_size != int(row["bytes"])
            or sha256(source) != row["sha256"]
        ):
            raise RuntimeError(f"immutable predecessor drift: {row['path']}")

    expression_manifest = read_tsv(
        CANDIDATE / "source_expression/expression_source_manifest.tsv"
    )
    for row in expression_manifest:
        source = CANDIDATE / row["relative_path"]
        if (
            not source.is_file()
            or source.stat().st_size != int(row["bytes"])
            or sha256(source) != row["sha256"]
        ):
            raise RuntimeError(f"expression source drift: {row['relative_path']}")

    mps_gate = read_tsv(CANDIDATE / "mps/mps_gate_status.tsv")[0]
    human_gate = read_tsv(CANDIDATE / "human_response/human_gate_status.tsv")[0]
    if mps_gate["status"] != "complete_nonconfirmatory":
        raise RuntimeError(f"unexpected MPS gate: {mps_gate}")
    if human_gate["status"] != "complete_nonconfirmatory":
        raise RuntimeError(f"unexpected human gate: {human_gate}")

    # Materialize the names promised by the workstream contract without
    # altering the branch-native, fully detailed tables.
    copy_contract(
        CANDIDATE / "mps/mps_sample_manifest.tsv",
        CANDIDATE / "mps_sample_manifest.tsv",
    )
    primary_condition_rows = [
        row
        for row in read_tsv(CANDIDATE / "mps/mps_condition_scores.tsv")
        if row["mode"] == "primary"
    ]
    write_tsv(
        CANDIDATE / "mps_condition_design.tsv",
        primary_condition_rows,
        list(primary_condition_rows[0]),
    )
    copy_contract(
        CANDIDATE / "mps/mps_lineage_effects.tsv",
        CANDIDATE / "mps_lineage_effects.tsv",
    )
    copy_contract(
        CANDIDATE / "human_response/human_pair_manifest.tsv",
        CANDIDATE / "human_pair_manifest.tsv",
    )
    human_topology_rows = [
        row
        for row in read_tsv(
            CANDIDATE / "human_response/human_participant_changes.tsv"
        )
        if row["mode"] == "primary" and row["normalization"] == "TMM_logCPM"
    ]
    write_tsv(
        CANDIDATE / "human_response_topology.tsv",
        human_topology_rows,
        list(human_topology_rows[0]),
    )

    mps_effects = read_tsv(CANDIDATE / "mps/mps_lineage_effects.tsv")
    human_sensitivity = read_tsv(CANDIDATE / "human_response/human_sensitivity.tsv")
    sensitivity_rows: list[dict[str, object]] = []
    for row in mps_effects:
        sensitivity_rows.append(
            {
                "branch": "GSE168285_MPS",
                "sensitivity": f"{row['contrast']}__{row['mode']}",
                "estimate": row["oriented_estimate"],
                "statistic": row["oriented_statistic"],
                "p": row["maxT_adjusted_p_one_sided"]
                if row["mode"] == "primary"
                else row["parametric_p"],
                "expected_direction_met": row["expected_direction_met"],
                "status": "complete"
                if row["mode"] != "shared_gene_removed"
                else "not_testable_prespecified_coverage_gate",
                "detail": "oriented effect; positive is the frozen expected direction",
            }
        )
    for row in human_sensitivity:
        sensitivity_rows.append(
            {
                "branch": "GSE175448_paired_human",
                "sensitivity": row["sensitivity"],
                "estimate": row["estimate"],
                "statistic": row["statistic"],
                "p": row["parametric_p"],
                "expected_direction_met": row["expected_direction_met"],
                "status": "complete",
                "detail": "negative is the frozen expected direction",
            }
        )
    for item in (
        "sample_level_voom_clustered_bootstrap",
        "leave_one_condition_family_out",
        "source_batch_adjustment",
        "family_overlap_reduction",
    ):
        sensitivity_rows.append(
            {
                "branch": "GSE168285_MPS",
                "sensitivity": item,
                "estimate": "",
                "statistic": "",
                "p": "",
                "expected_direction_met": "",
                "status": "not_run_after_primary_direction_falsification",
                "detail": (
                    "nondecisive sensitivity stopped after the frozen primary TGF-beta "
                    "contrast was strongly opposite and the conjunction p was 1"
                ),
            }
        )
    write_tsv(
        CANDIDATE / "sensitivity.tsv",
        sensitivity_rows,
        [
            "branch",
            "sensitivity",
            "estimate",
            "statistic",
            "p",
            "expected_direction_met",
            "status",
            "detail",
        ],
    )

    perturbseq = [
        {
            "dataset": "GSE238219",
            "analysis": "complete_11_target_hepatocyte_vs_remodeling_observability",
            "status": "not_run_nonblocking_after_primary_failure",
            "n_source_targets": 11,
            "n_corrected_genetic_only_targets": 4,
            "inference_authorized": "false",
            "estimate": "",
            "p": "",
            "q": "",
            "detail": (
                "supplementary-only arm could not rescue either required branch; "
                "stopped under the sealed public-data termination rule"
            ),
        }
    ]
    write_tsv(
        CANDIDATE / "perturbseq_observability.tsv",
        perturbseq,
        list(perturbseq[0]),
    )

    primary_mps = {
        row["contrast"]: row
        for row in mps_effects
        if row["mode"] == "primary"
    }
    verdict_rows = [
        {
            "component_id": "controlled_multicellular_assembly",
            "status": "fail_frozen_direction_and_multiplicity",
            "gate_pass": "false",
            "evidence": (
                f"fat_oriented_beta={primary_mps['fat_hepatocyte_specificity']['oriented_estimate']};"
                f"fat_maxT_p={primary_mps['fat_hepatocyte_specificity']['maxT_adjusted_p_one_sided']};"
                f"npc_oriented_beta={primary_mps['npc_remodeling_specificity']['oriented_estimate']};"
                f"npc_maxT_p={primary_mps['npc_remodeling_specificity']['maxT_adjusted_p_one_sided']};"
                f"tgfb_oriented_beta={primary_mps['tgfb_remodeling_specificity']['oriented_estimate']};"
                f"tgfb_maxT_p={primary_mps['tgfb_remodeling_specificity']['maxT_adjusted_p_one_sided']};"
                f"conjunction_p={mps_gate['conjunction_p']}"
            ),
            "permitted_claim": (
                "model-boundary result: the frozen source-lineage axes do not encode "
                "the proposed fat-versus-profibrotic assembly topology"
            ),
        },
        {
            "component_id": "paired_human_remodeling_disassembly",
            "status": "fail_opposite_direction",
            "gate_pass": "false",
            "evidence": (
                f"improved_minus_not_beta={human_gate['primary_estimate']};"
                f"exact_p={human_gate['primary_exact_p']};"
                "expected_beta_negative=true"
            ),
            "permitted_claim": (
                "no lineage-selective disassembly support in this outcome-unseen "
                "19-participant cohort"
            ),
        },
        {
            "component_id": "supplementary_perturbseq_observability",
            "status": "not_run_nonblocking_after_primary_failure",
            "gate_pass": "false",
            "evidence": "supplementary-only analysis had no path to paper-level promotion",
            "permitted_claim": "none",
        },
        {
            "component_id": "complete_lineage_assembly_disassembly",
            "status": "complete_no_promotion",
            "gate_pass": "false",
            "evidence": "controlled_MPS=false; paired_human=false; conjunction=false",
            "permitted_claim": "none beyond the pre-existing Resource thesis",
        },
        {
            "component_id": "public_data_escalation_program",
            "status": "terminal_stop",
            "gate_pass": "false",
            "evidence": "Plans 41-44 exhausted their frozen public-data promotion gates",
            "permitted_claim": (
                "a top-journal risk-to-state mechanism now requires a new, "
                "source-independent multicellular perturbation"
            ),
        },
    ]
    write_tsv(
        CANDIDATE / "promotion_verdict.tsv",
        verdict_rows,
        ["component_id", "status", "gate_pass", "evidence", "permitted_claim"],
    )

    terminal_source_rows = []
    for row in read_tsv(CANDIDATE / "source_gate_status.tsv"):
        output = dict(row)
        output["outcomes_accessed"] = (
            "true" if row["dataset"] in {"GSE168285", "GSE175448"} else "false"
        )
        output["terminal_status"] = (
            "analyzed_nonconfirmatory"
            if row["dataset"] in {"GSE168285", "GSE175448"}
            else "not_run_nonblocking_after_primary_failure"
        )
        terminal_source_rows.append(output)
    write_tsv(
        CANDIDATE / "terminal_source_gate_status.tsv",
        terminal_source_rows,
        list(terminal_source_rows[0]),
    )

    combined_sources = []
    for row in read_tsv(CANDIDATE / "source_manifest.tsv"):
        combined_sources.append(
            {
                "relative_path": row["relative_path"],
                "bytes": row["bytes"],
                "sha256": row["sha256"],
                "content_class": row["content_class"],
                "source_url": row["source_url"],
                "retrieval_date": row["retrieval_date"],
            }
        )
    for row in expression_manifest:
        combined_sources.append(
            {
                "relative_path": row["relative_path"],
                "bytes": row["bytes"],
                "sha256": row["sha256"],
                "content_class": "processed_expression_outcome",
                "source_url": row["source_url"],
                "retrieval_date": row["retrieval_date"],
            }
        )
    write_tsv(
        CANDIDATE / "terminal_source_manifest.tsv",
        combined_sources,
        [
            "relative_path",
            "bytes",
            "sha256",
            "content_class",
            "source_url",
            "retrieval_date",
        ],
    )

    script_rows = []
    for path in sorted(SCRIPT_ROOT.iterdir()):
        if path.is_file() and path.suffix in {".py", ".R", ".sh", ".sbatch"}:
            script_rows.append(
                {
                    "artifact": str(path.relative_to(ROOT)),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    write_tsv(
        CANDIDATE / "execution_manifest.tsv",
        script_rows,
        ["artifact", "bytes", "sha256"],
    )
    write_tsv(
        CANDIDATE / "environment_manifest.tsv",
        [
            {
                "component": "python",
                "version": platform.python_version(),
                "detail": sys.executable,
            },
            {
                "component": "platform",
                "version": platform.platform(),
                "detail": "",
            },
            {
                "component": "preflight_manifest",
                "version": preflight["release_manifest_sha256"],
                "detail": "outcome-blind specification",
            },
        ],
        ["component", "version", "detail"],
    )

    adjudicated = {
        "release_id": "multicellular-assembly-response-2026-08-09",
        "status": "complete_no_promotion",
        "adjudicated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "mps_gate_pass": False,
        "paired_human_gate_pass": False,
        "figure5_promotion": False,
        "canonical_promotion_authorized": False,
        "public_data_escalation_ended": True,
        "journal_posture": "Cell_Genomics_Resource",
        "next_decisive_evidence": (
            "source-independent risk-oriented perturbation in a multicellular human "
            "liver model with temporal and orthogonal remodeling readouts"
        ),
    }
    atomic_json(CANDIDATE / "ADJUDICATED.json", adjudicated)
    print(json.dumps(adjudicated, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

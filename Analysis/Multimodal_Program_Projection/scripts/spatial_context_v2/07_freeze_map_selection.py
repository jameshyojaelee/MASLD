#!/usr/bin/env python3
"""Freeze Figure 4 map identities from the sealed Figure 2 source only.

This producer is deliberately outcome-blind with respect to Plan 11 and all
external spatial adapters.  It selects every sealed Plan 20 program that is
both robust for display and eligible for external testing.  In the current v2
registry that is the complete two-program hepatocyte family, so retaining both
is less discretionary than choosing a single visually favourable program.
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

from contract_lib import (
    ContractError,
    RELEASE_ID,
    parse_bool,
    parse_float,
    program_universe,
    read_tsv,
    sha256_file,
    validate_hotspot_seal,
    write_tsv,
)


MAP_COLUMNS = (
    "release_id",
    "selection_freeze_id",
    "selection_rank",
    "program_uid",
    "membership_sha256",
    "legacy_program_id",
    "program_label",
    "cell_type",
    "fig2_beta",
    "fig2_qvalue",
    "fig2_hc3_qvalue",
    "fig2_stability_score",
    "fig2_source_stability",
    "selection_rule",
    "selection_reason",
    "fig2_source_sha256",
    "registry_sha256",
    "ready_sha256",
    "external_outcomes_read",
    "yakubovsky_program_outcomes_present_at_freeze",
    "frozen_at_utc",
    "producer",
    "producer_sha256",
)


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _yakubovsky_program_outcomes_present(project_root: Path) -> list[Path]:
    root = (
        project_root
        / "Analysis/Spatial/candidates"
        / RELEASE_ID
        / "yakubovsky2026"
    )
    protected = (
        "per_sample_program_scores.tsv",
        "donor_program_effects.tsv",
        "program_effects.tsv",
        "model_stage_seal.tsv",
        "terminal_stage_seal.tsv",
        "READY",
    )
    return [root / filename for filename in protected if (root / filename).exists()]


def _manifest_hash_for(relative_path: str, hotspot_root: Path) -> str:
    _, rows = read_tsv(
        hotspot_root / "release_manifest.tsv",
        ("relative_path", "sha256", "release_id"),
    )
    matches = [
        row for row in rows
        if row["relative_path"] == relative_path and row["release_id"] == RELEASE_ID
    ]
    if len(matches) != 1:
        raise ContractError(
            f"Hotspot release manifest must contain one row for {relative_path}, "
            f"found {len(matches)}"
        )
    return matches[0]["sha256"]


def expected_selection(
    project_root: Path,
    hotspot_root: Path,
    frozen_at_utc: str,
) -> list[dict[str, object]]:
    seal = validate_hotspot_seal(hotspot_root)
    fig2_path = hotspot_root / "fig2_program_source.tsv"
    fig2_relative = fig2_path.relative_to(project_root).as_posix()
    fig2_sha256 = sha256_file(fig2_path)
    if _manifest_hash_for(fig2_relative, hotspot_root) != fig2_sha256:
        raise ContractError("Figure 2 source hash does not match the sealed release manifest")

    _, fig2_rows = read_tsv(
        fig2_path,
        (
            "release_id",
            "program_uid",
            "membership_sha256",
            "cell_type",
            "module",
            "module_name",
            "beta",
            "qvalue",
            "hc3_qvalue",
            "stability_score",
            "robust_display",
            "source_stability",
        ),
    )
    fig2_by_uid = {row["program_uid"]: row for row in fig2_rows}
    if len(fig2_by_uid) != len(fig2_rows):
        raise ContractError("Figure 2 source contains duplicate program_uid values")

    eligible = program_universe(seal, "external_test_eligible")
    selected = []
    for uid, registry_row in eligible.items():
        if uid not in fig2_by_uid:
            raise ContractError(f"external-test program absent from Figure 2 source: {uid}")
        row = fig2_by_uid[uid]
        if not parse_bool(row["robust_display"], f"fig2.robust_display[{uid}]"):
            raise ContractError(f"external-test program is not robust_display in Figure 2: {uid}")
        if row["membership_sha256"] != registry_row["membership_sha256"]:
            raise ContractError(f"Figure 2 membership hash drift for {uid}")
        selected.append(row)

    robust_fig2 = {
        row["program_uid"]
        for row in fig2_rows
        if parse_bool(row["robust_display"], f"fig2.robust_display[{row['program_uid']}]")
    }
    if robust_fig2 != set(eligible):
        raise ContractError(
            "Map freeze requires the complete robust-display/external-test family; "
            f"Figure2={sorted(robust_fig2)}, external={sorted(eligible)}"
        )
    if len(selected) != 2:
        raise ContractError(
            f"Current prespecified map rule expects the sealed two-program family, found {len(selected)}"
        )

    selected.sort(
        key=lambda row: (
            -abs(parse_float(row["beta"], f"fig2.beta[{row['program_uid']}]") or 0.0),
            -(parse_float(row["stability_score"], f"fig2.stability[{row['program_uid']}]") or 0.0),
            row["program_uid"],
        )
    )
    producer = Path(__file__).resolve()
    producer_relative = producer.relative_to(project_root).as_posix()
    producer_sha256 = sha256_file(producer)
    selection_rule = (
        "complete_sealed_external_test_eligible_family;"
        "order_abs_fig2_beta_desc_then_stability_desc_then_program_uid"
    )
    rows = []
    for rank, row in enumerate(selected, start=1):
        rows.append(
            {
                "release_id": RELEASE_ID,
                "selection_freeze_id": "SP-INT-06-map-freeze-v1",
                "selection_rank": rank,
                "program_uid": row["program_uid"],
                "membership_sha256": row["membership_sha256"],
                "legacy_program_id": f"{row['cell_type']}:{row['module']}",
                "program_label": row["module_name"],
                "cell_type": row["cell_type"],
                "fig2_beta": row["beta"],
                "fig2_qvalue": row["qvalue"],
                "fig2_hc3_qvalue": row["hc3_qvalue"],
                "fig2_stability_score": row["stability_score"],
                "fig2_source_stability": row["source_stability"],
                "selection_rule": selection_rule,
                "selection_reason": (
                    "retain_both_only_sealed_robust_hepatocyte_programs;"
                    "no_external_spatial_outcome_used"
                ),
                "fig2_source_sha256": fig2_sha256,
                "registry_sha256": seal.registry_sha256,
                "ready_sha256": seal.ready_sha256,
                "external_outcomes_read": False,
                "yakubovsky_program_outcomes_present_at_freeze": False,
                "frozen_at_utc": frozen_at_utc,
                "producer": producer_relative,
                "producer_sha256": producer_sha256,
            }
        )
    return rows


def check_existing(project_root: Path, candidate_root: Path, hotspot_root: Path) -> None:
    path = candidate_root / "map_selection_manifest.tsv"
    _, observed = read_tsv(path, MAP_COLUMNS)
    if not observed:
        raise ContractError("map_selection_manifest.tsv is empty")
    frozen_values = {row["frozen_at_utc"] for row in observed}
    if len(frozen_values) != 1:
        raise ContractError("map selection rows do not share one frozen_at_utc")
    expected = expected_selection(project_root, hotspot_root, next(iter(frozen_values)))
    normalized_expected = [{column: str(row[column]) for column in MAP_COLUMNS} for row in expected]
    for row in normalized_expected:
        row["external_outcomes_read"] = "FALSE"
        row["yakubovsky_program_outcomes_present_at_freeze"] = "FALSE"
    if observed != normalized_expected:
        raise ContractError("map_selection_manifest.tsv does not rederive from the sealed Figure 2 source")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("write", "check"))
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--candidate-root", type=Path, default=None)
    parser.add_argument("--hotspot-root", type=Path, default=None)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    candidate_root = args.candidate_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "spatial_context"
    )
    hotspot_root = args.hotspot_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "hotspot"
    )
    try:
        if args.mode == "check":
            check_existing(project_root, candidate_root, hotspot_root)
            print("PASS: outcome-blind map selection rederived")
            return 0

        output = candidate_root / "map_selection_manifest.tsv"
        if output.exists():
            raise ContractError(f"refusing to overwrite frozen map selection: {output}")
        unexpected = _yakubovsky_program_outcomes_present(project_root)
        if unexpected:
            raise ContractError(
                "refusing outcome-blind map freeze after Yakubovsky program outcomes exist: "
                + ", ".join(str(path) for path in unexpected)
            )
        candidate_root.mkdir(parents=True, exist_ok=True)
        rows = expected_selection(project_root, hotspot_root, utc_now())
        staging = candidate_root / f".map_selection_manifest.incomplete.{os.getpid()}.tsv"
        write_tsv(staging, MAP_COLUMNS, rows)
        os.replace(staging, output)
        check_existing(project_root, candidate_root, hotspot_root)
        print(f"PASS: froze {len(rows)} map identities before external program outcomes")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

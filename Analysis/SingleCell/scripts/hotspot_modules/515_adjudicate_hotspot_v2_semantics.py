#!/usr/bin/env python3
"""Add a non-mutating semantic adjudication layer to the frozen Hotspot v2 registry.

The original registry is an outcome-blind release artifact and is never rewritten.
This companion table deprecates only the historical interpretation of
``tested_negative`` as evidence of absence; it does not change selection,
external-test eligibility, membership, weights, or any fitted statistic.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import tempfile
from pathlib import Path
from typing import Iterable, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]
RELEASE_ID = "program-context-v2-candidate-2026-08-07"
CANDIDATE_REL = Path(
    "Analysis/Multimodal_Program_Projection/candidates"
) / RELEASE_ID / "hotspot"
EXPECTED_REGISTRY_SHA256 = (
    "b134b5f8e46a716e857af22271aa54e1294432b9f9ab735559d7f603cdd29fe7"
)
EXPECTED_REGISTRY_BYTES = 188_469

OUTPUT_NAME = "program_registry_v2_semantic_adjudication.tsv"
MANIFEST_NAME = "semantic_adjudication_manifest.tsv"


class ContractError(RuntimeError):
    """Raised when the frozen candidate no longer matches this adjudication."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ContractError(f"missing TSV header: {path}")
        return [dict(row) for row in reader]


def atomic_write_tsv(
    path: Path,
    rows: Iterable[Mapping[str, object]],
    fieldnames: Sequence[str],
) -> None:
    materialized = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(materialized)
    os.replace(temporary, path)


def parse_bool(value: str, field: str) -> bool:
    if value == "TRUE":
        return True
    if value == "FALSE":
        return False
    raise ContractError(f"{field} is not canonical TRUE/FALSE: {value!r}")


def parse_float(value: str, field: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise ContractError(f"{field} is not numeric: {value!r}") from error
    if not math.isfinite(result):
        raise ContractError(f"{field} is not finite: {value!r}")
    return result


def assert_candidate_root(project_root: Path, candidate_root: Path) -> Path:
    expected = (project_root / CANDIDATE_REL).resolve()
    observed = candidate_root.resolve()
    if observed != expected:
        raise ContractError(
            f"candidate root must be exactly {expected}; received {observed}"
        )
    return observed


def adjudicate(row: Mapping[str, str]) -> tuple[str, str, str, str]:
    estimable = parse_bool(row["primary_estimable"], "primary_estimable")
    if not estimable:
        return (
            "untestable",
            "not_selected_untestable",
            "false",
            "the donor-level stage coefficient was not estimable",
        )

    primary_q = parse_float(row["primary_qvalue"], "primary_qvalue")
    primary_selected = parse_bool(row["primary_selected"], "primary_selected")
    selected_unstable = parse_bool(row["selected_unstable"], "selected_unstable")
    selected_hc3_fragile = parse_bool(
        row["selected_hc3_fragile"], "selected_hc3_fragile"
    )
    robust_display = parse_bool(row["robust_display"], "robust_display")
    external_eligible = parse_bool(
        row["external_test_eligible"], "external_test_eligible"
    )

    if primary_q >= 0.05:
        if primary_selected or selected_unstable or robust_display or external_eligible:
            raise ContractError(
                f"q-failing row has incompatible selection flags: {row['program_uid']}"
            )
        return (
            "indeterminate_nonconfirmatory",
            "not_selected_q_nonsignificant",
            "false",
            "finite q>=0.05 is nonconfirmatory; no equivalence or informative-negative gate was specified",
        )

    if selected_unstable:
        if primary_selected or robust_display or external_eligible:
            raise ContractError(
                f"stability-failed row has incompatible selection flags: {row['program_uid']}"
            )
        return (
            "indeterminate_stability_failed",
            "q_significant_stability_failed",
            "false",
            "primary q<0.05 but the prespecified Hotspot LOO median stability threshold failed",
        )

    if not primary_selected:
        raise ContractError(
            f"q-significant stable row is not primary-selected: {row['program_uid']}"
        )

    if robust_display:
        if not external_eligible or row["cell_type"] != "hepatocytes":
            raise ContractError(
                f"robust row is not the sealed hepatocyte external subset: {row['program_uid']}"
            )
        return (
            "supported_internal_stage_association",
            "primary_selected_internal_robust",
            "false",
            "selected donor-level association passed scoring-direction and HC3 gates; external outcomes are not adjudicated here",
        )

    if selected_hc3_fragile:
        if external_eligible:
            raise ContractError(
                f"HC3-fragile row entered the external subset: {row['program_uid']}"
            )
        return (
            "indeterminate_hc3_fragile",
            "primary_selected_hc3_fragile",
            "false",
            "primary-selected association did not pass the prespecified HC3 robustness gate",
        )

    raise ContractError(
        f"primary-selected row has no recognized robustness state: {row['program_uid']}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument(
        "--candidate-root", type=Path, default=PROJECT_ROOT / CANDIDATE_REL
    )
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    if project_root != PROJECT_ROOT.resolve():
        raise ContractError(f"project-root drift: {project_root}")
    candidate = assert_candidate_root(project_root, args.candidate_root)
    registry_path = candidate / "program_registry_v2.tsv"
    ready_path = candidate / "READY"
    if not registry_path.is_file() or not ready_path.is_file():
        raise ContractError("frozen registry or READY seal is missing")
    if registry_path.stat().st_size != EXPECTED_REGISTRY_BYTES:
        raise ContractError("frozen registry byte-size drift")
    registry_sha = sha256_file(registry_path)
    if registry_sha != EXPECTED_REGISTRY_SHA256:
        raise ContractError(
            f"frozen registry hash drift: {registry_sha} != {EXPECTED_REGISTRY_SHA256}"
        )

    ready_rows = read_tsv(ready_path)
    if len(ready_rows) != 1:
        raise ContractError("READY must contain exactly one data row")
    ready = ready_rows[0]
    if ready.get("registry_sha256") != registry_sha:
        raise ContractError("READY does not anchor the frozen registry")
    if ready.get("status") != "ready_for_external_testing":
        raise ContractError("frozen registry did not pass its original validation")

    rows = read_tsv(registry_path)
    required = {
        "release_id",
        "program_uid",
        "cell_type",
        "module",
        "module_name",
        "membership_sha256",
        "primary_estimable",
        "primary_qvalue",
        "primary_hc3_qvalue",
        "stability_median",
        "primary_selected",
        "selected_unstable",
        "selected_hc3_fragile",
        "robust_display",
        "external_test_eligible",
        "tested_negative",
        "registry_state",
    }
    missing = required - set(rows[0]) if rows else required
    if missing:
        raise ContractError(f"registry missing fields: {sorted(missing)}")
    if len(rows) != 117 or len({row["program_uid"] for row in rows}) != 117:
        raise ContractError("registry is not 117 unique program_uid rows")

    output_rows: list[dict[str, object]] = []
    for row in rows:
        state, tier, negative_authorized, reason = adjudicate(row)
        legacy_negative = parse_bool(row["tested_negative"], "tested_negative")
        output_rows.append(
            {
                "release_id": row["release_id"],
                "program_uid": row["program_uid"],
                "cell_type": row["cell_type"],
                "module": row["module"],
                "module_name": row["module_name"],
                "membership_sha256": row["membership_sha256"],
                "registry_sha256": registry_sha,
                "primary_estimable": row["primary_estimable"],
                "primary_qvalue": row["primary_qvalue"],
                "primary_hc3_qvalue": row["primary_hc3_qvalue"],
                "stability_median": row["stability_median"],
                "primary_selected": row["primary_selected"],
                "selected_unstable": row["selected_unstable"],
                "selected_hc3_fragile": row["selected_hc3_fragile"],
                "robust_display": row["robust_display"],
                "external_test_eligible": row["external_test_eligible"],
                "legacy_tested_negative": row["tested_negative"],
                "legacy_registry_state": row["registry_state"],
                "legacy_negative_semantics_deprecated": (
                    "true" if legacy_negative else "false"
                ),
                "adjudicated_state": state,
                "selection_tier": tier,
                "tested_negative_authorized": negative_authorized,
                "evidence_scope": "internal_donor_level_stage_association_only",
                "external_outcome_state": "not_adjudicated_here",
                "adjudication_reason": reason,
                "semantic_contract": "hotspot_v2_postfreeze_semantics_v1",
            }
        )

    output_path = candidate / OUTPUT_NAME
    atomic_write_tsv(output_path, output_rows, list(output_rows[0]))

    code_paths = [
        SCRIPT_DIR / "515_adjudicate_hotspot_v2_semantics.py",
        SCRIPT_DIR / "516_validate_hotspot_v2_semantics.py",
        SCRIPT_DIR / "run_hotspot_v2_semantic_adjudication.sbatch",
    ]
    for path in code_paths:
        if not path.is_file():
            raise ContractError(f"required reviewed code is missing: {path}")
    manifest_rows = []
    for path, role, mutable in [
        (registry_path, "frozen_registry_input", "false"),
        (ready_path, "frozen_registry_ready_seal", "false"),
        (output_path, "derived_semantic_adjudication", "true"),
        (code_paths[0], "semantic_adjudication_producer", "true"),
        (code_paths[1], "independent_semantic_validator", "true"),
        (code_paths[2], "review_only_slurm_wrapper", "true"),
    ]:
        manifest_rows.append(
            {
                "relative_path": str(path.resolve().relative_to(project_root)),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "role": role,
                "frozen_registry_mutated": mutable if role.startswith("frozen_") else "not_applicable",
                "canonical_promotion_status": "not_promoted",
            }
        )
    manifest_path = candidate / MANIFEST_NAME
    atomic_write_tsv(manifest_path, manifest_rows, list(manifest_rows[0]))

    if sha256_file(registry_path) != registry_sha:
        raise ContractError("frozen registry changed during semantic adjudication")
    print(
        f"PASS: wrote {len(output_rows)} derived semantic rows; "
        f"frozen registry remains {registry_sha}"
    )


if __name__ == "__main__":
    main()

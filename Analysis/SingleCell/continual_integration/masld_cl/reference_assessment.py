"""Fail-closed healthy-reference roster assessment.

This module evaluates reference composition only. It does not read case, stage,
program, hero-gene, or Cas13 outcomes and it never changes the atlas contract.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


class ReferenceAssessmentError(ValueError):
    """Raised when a proposed reference roster violates the frozen policy."""


CURRENT_REFERENCE_DONORS = {
    "GSE185477_D02",
    "GSE185477_D03",
    "GSE185477_D04",
    "Liver_Atlas_H06",
    "Liver_Atlas_H07",
    "Liver_Atlas_H10",
    "Liver_Atlas_H22",
}


def sha256_path(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: str | Path) -> dict[str, Any]:
    with Path(path).open() as handle:
        return json.load(handle)


def _read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_policy(path: str | Path) -> dict[str, Any]:
    policy = _read_json(path)
    if policy.get("schema_version") != "masld-reference-assessment-v1":
        raise ReferenceAssessmentError("invalid reference-policy schema_version")
    if policy.get("primary_reference_unchanged") is not True:
        raise ReferenceAssessmentError("the current seven-donor primary reference must remain unchanged")
    firewall = policy.get("selection_firewall", {})
    if firewall.get("lock_before_lambda_replay_screen") is not True:
        raise ReferenceAssessmentError("reference selection must lock before lambda/replay tuning")
    forbidden = set(firewall.get("forbidden_inputs", []))
    required_forbidden = {
        "case distance", "disease stage", "program scores", "hero genes",
        "Cas13 results", "UMAP appearance",
    }
    if forbidden != required_forbidden:
        raise ReferenceAssessmentError("reference selection firewall is incomplete or changed")
    external = policy.get("external_arm", {})
    if external.get("source") != "HLiCA_v1_CELLxGENE":
        raise ReferenceAssessmentError("external reference source must be the pinned HLiCA v1 asset")
    expected = (
        int(external["current_reference_donors"])
        + int(external["new_unique_donors"])
    )
    if expected != int(external["expected_reference_donors"]):
        raise ReferenceAssessmentError("external reference donor arithmetic is inconsistent")
    audit = external.get("compatibility_audit", {})
    if (
        external.get("status")
        != "phenotype_qualified_but_ineligible_for_identical_current_4000_hvg_comparison"
        or audit.get("eligible_for_identical_current_4000_hvgs") is not False
        or int(audit.get("mapped_current_hvgs", -1))
        + int(audit.get("missing_current_hvgs", -1)) != 4000
        or int(audit.get("ambiguous_current_hvgs", -1)) != 0
    ):
        raise ReferenceAssessmentError("external feature-compatibility gate changed")
    return policy


def validate_registry(path: str | Path) -> list[dict[str, str]]:
    rows = _read_tsv(path)
    required = {
        "candidate_id", "study", "reference_use", "phenotype_evidence",
        "source_url", "data_url",
    }
    if not rows or required - rows[0].keys():
        raise ReferenceAssessmentError("external reference registry is empty or missing columns")
    ids = [row["candidate_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ReferenceAssessmentError("external reference registry has duplicate candidate_id values")
    hlica = next((row for row in rows if row["candidate_id"] == "hlica_all"), None)
    if hlica is None or hlica["reference_use"] != "exclude_wholesale":
        raise ReferenceAssessmentError("whole-HLiCA ingestion must be explicitly forbidden")
    permitted = {"andrews_2022", "andrews_2024"}
    observed = {
        row["candidate_id"] for row in rows
        if row["reference_use"] == "include_external_clean"
    }
    if observed != permitted:
        raise ReferenceAssessmentError("the external clean subset changed without a policy revision")
    return rows


def verify_frozen_inputs(policy: dict[str, Any], repository_root: str | Path) -> dict[str, str]:
    root = Path(repository_root)
    observed: dict[str, str] = {}
    for name, spec in policy["frozen_inputs"].items():
        path = root / spec["path"]
        if not path.is_file():
            raise ReferenceAssessmentError(f"missing frozen input: {name}: {path}")
        digest = sha256_path(path)
        if digest != spec["sha256"]:
            raise ReferenceAssessmentError(
                f"frozen input hash mismatch for {name}: {digest} != {spec['sha256']}"
            )
        observed[name] = digest
    return observed


def _bool(value: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise ReferenceAssessmentError(f"invalid manifest boolean: {value!r}")


def _manifest_by_donor(rows: Iterable[dict[str, str]]) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        donor = row["donor_id"]
        if donor in result:
            raise ReferenceAssessmentError(f"duplicate donor manifest row: {donor}")
        result[donor] = row
    return result


def validate_local_arms(
    policy: dict[str, Any], donor_manifest: str | Path
) -> tuple[dict[str, Any], dict[str, set[str]]]:
    rows = _read_tsv(donor_manifest)
    donors = _manifest_by_donor(rows)
    observed_current = {
        donor for donor, row in donors.items() if _bool(row["strict_reference"])
    }
    if observed_current != CURRENT_REFERENCE_DONORS:
        raise ReferenceAssessmentError("current strict-reference donor roster changed")
    query_controls = {
        donor for donor, row in donors.items() if _bool(row["query_control"])
    }
    summaries: dict[str, Any] = {}
    arm_donors: dict[str, set[str]] = {}
    for arm_name, arm in policy["local_arms"].items():
        selected = observed_current | set(arm["added_donors"])
        missing = sorted(selected - donors.keys())
        if missing:
            raise ReferenceAssessmentError(f"{arm_name} has missing donors: {missing}")
        overlap = sorted(selected & query_controls)
        if overlap:
            raise ReferenceAssessmentError(f"{arm_name} leaks query controls: {overlap}")
        for donor in selected:
            row = donors[donor]
            if not _bool(row["analysis_eligible"]):
                raise ReferenceAssessmentError(f"{arm_name} includes ineligible donor {donor}")
            if row["harmonized_stage"] != "Healthy":
                raise ReferenceAssessmentError(f"{arm_name} includes non-healthy donor {donor}")
        cells = sum(int(donors[donor]["n_cells_analyzed"]) for donor in selected)
        if len(selected) != int(arm["expected_donors"]):
            raise ReferenceAssessmentError(f"{arm_name} donor total changed")
        if cells != int(arm["expected_cells"]):
            raise ReferenceAssessmentError(f"{arm_name} cell total changed")
        by_dataset: dict[str, dict[str, int]] = defaultdict(lambda: {"donors": 0, "cells": 0})
        for donor in selected:
            row = donors[donor]
            dataset = row["dataset"]
            by_dataset[dataset]["donors"] += 1
            by_dataset[dataset]["cells"] += int(row["n_cells_analyzed"])
        donor_cells = sorted(
            ((donor, int(donors[donor]["n_cells_analyzed"])) for donor in selected),
            key=lambda value: (-value[1], value[0]),
        )
        summaries[arm_name] = {
            "status": arm["status"],
            "donors": len(selected),
            "cells": cells,
            "datasets": dict(sorted(by_dataset.items())),
            "largest_donor": donor_cells[0][0],
            "largest_donor_cells": donor_cells[0][1],
            "largest_donor_fraction": donor_cells[0][1] / cells,
            "donor_ids": sorted(selected),
        }
        arm_donors[arm_name] = selected
    return summaries, arm_donors


def summarize_lineages(
    cell_manifest: str | Path, arm_donors: dict[str, set[str]]
) -> dict[str, list[dict[str, Any]]]:
    path = Path(cell_manifest)
    opener = gzip.open if path.suffix == ".gz" else open
    counts: dict[tuple[str, str, str], int] = Counter()
    donor_sets: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    donor_to_arms: dict[str, list[str]] = defaultdict(list)
    for arm, donors in arm_donors.items():
        for donor in donors:
            donor_to_arms[donor].append(arm)
    with opener(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            donor = row["donor_id"]
            for arm in donor_to_arms.get(donor, []):
                key = (arm, row["dataset"], row["cell_type"])
                counts[key] += 1
                donor_sets[key].add(donor)
    result: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for (arm, dataset, cell_type), cells in sorted(counts.items()):
        result[arm].append({
            "dataset": dataset,
            "cell_type": cell_type,
            "cells": cells,
            "donors": len(donor_sets[(arm, dataset, cell_type)]),
        })
    return dict(result)


def build_assessment(
    policy_path: str | Path,
    registry_path: str | Path,
    repository_root: str | Path,
) -> dict[str, Any]:
    policy = load_policy(policy_path)
    registry = validate_registry(registry_path)
    hashes = verify_frozen_inputs(policy, repository_root)
    root = Path(repository_root)
    donor_manifest = root / policy["frozen_inputs"]["donor_manifest"]["path"]
    cell_manifest = root / policy["frozen_inputs"]["cell_manifest"]["path"]
    local, arm_donors = validate_local_arms(policy, donor_manifest)
    lineages = summarize_lineages(cell_manifest, arm_donors)
    uses = Counter(row["reference_use"] for row in registry)
    return {
        "schema_version": policy["schema_version"],
        "primary_reference_unchanged": True,
        "frozen_input_sha256": hashes,
        "local_arms": local,
        "lineage_composition": lineages,
        "external_registry_counts": dict(sorted(uses.items())),
        "external_arm": policy["external_arm"],
        "feature_comparison": policy["feature_comparison"],
        "selection_firewall": policy["selection_firewall"],
        "recommendation": (
            "Finish the current seven-donor pilot, then lock the reference roster using "
            "control-only metrics across the identical-4000-HVG local arms. Use the "
            "nine-donor GSE174748 expansion as the primary challenger and keep GSE136103 "
            "as a non-hepatocyte protocol sensitivity. The deduplicated "
            "Andrews_2022/Andrews_2024 arm is post-lock sensitivity only and requires a "
            "separate common-universe 4,000-HVG lock plus identical refits of every arm."
        ),
    }

#!/usr/bin/env python3
"""Independently validate the Plan 45 lineage/accessibility routing release."""

from __future__ import annotations

import json
import math

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def main() -> None:
    seal_path = CANDIDATE_ROOT / "LINEAGE_ACCESSIBILITY_ROUTING_SEALED.json"
    if not seal_path.is_file():
        raise RuntimeError("Missing routing seal")
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "lineage_accessibility_routing_complete_targets_not_selected":
        raise RuntimeError("Invalid routing status")
    if seal.get("experimental_targets_frozen") is not False:
        raise RuntimeError("Routing release froze targets")

    output_names = {
        "shared_signal_accessibility": "shared_signal_accessibility.tsv",
        "lineage_accessibility_routing": "lineage_accessibility_routing.tsv",
        "routing_input_manifest": "routing_input_manifest.tsv",
        "routing_gate_status": "routing_gate_status.tsv",
    }
    for key, filename in output_names.items():
        path = CANDIDATE_ROOT / filename
        if not path.is_file() or sha256_file(path) != seal["output_sha256"][key]:
            raise RuntimeError(f"Routing output hash mismatch: {filename}")

    variants = read_tsv(CANDIDATE_ROOT / "shared_signal_accessibility.tsv")
    routes = read_tsv(CANDIDATE_ROOT / "lineage_accessibility_routing.tsv")
    manifests = read_tsv(CANDIDATE_ROOT / "routing_input_manifest.tsv")
    gates = read_tsv(CANDIDATE_ROOT / "routing_gate_status.tsv")
    if len(routes) != int(seal["n_prespecified_pairs"]) or len(gates) != 1:
        raise RuntimeError("Routing cardinality mismatch")
    if len({row["orientation_uid"] for row in routes}) != len(routes):
        raise RuntimeError("Duplicate routed loci")
    if any(row["target_freeze_status"] != "prohibited_until_guideability_and_exact_edit_gate" for row in routes):
        raise RuntimeError("A routed locus opened target selection")

    route_by_uid = {row["orientation_uid"]: row for row in routes}
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in variants:
        grouped.setdefault(row["orientation_uid"], []).append(row)
    for uid, route in route_by_uid.items():
        rows = grouped.get(uid, [])
        pass_variants = [row for row in rows if yes(row["accessibility_variant_pass"])]
        if len(pass_variants) != int(route["n_accessibility_variant_pass"]):
            raise RuntimeError(f"Accessible-variant count mismatch: {uid}")
        if any(float(row["shared_posterior"]) < 0.05 for row in pass_variants):
            raise RuntimeError(f"Posterior threshold violation: {uid}")
        if any(not yes(row["primary_lineage_peak_overlap"]) for row in pass_variants):
            raise RuntimeError(f"Primary peak overlap violation: {uid}")
        if route["dominant_celltype"] == "Hepatocytes" and any(
            not yes(row["independent_hepatocyte_peak_overlap"]) for row in pass_variants
        ):
            raise RuntimeError(f"Hepatocyte replication violation: {uid}")
        expected_pass = (
            yes(route["orientation_gate_pass"])
            and yes(route["pair_family_orientation_concordant"])
            and yes(route["lineage_expression_gate"])
            and yes(route["lineage_specificity_gate"])
            and route["mapping_status"] == "exact"
            and len(pass_variants) > 0
        )
        if yes(route["routing_gate_pass"]) != expected_pass:
            raise RuntimeError(f"Routing gate arithmetic mismatch: {uid}")

    for row in manifests:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Routing source drift: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Routing source hash drift: {path}")

    n_orientation = sum(yes(row["orientation_gate_pass"]) for row in routes)
    n_routing = sum(yes(row["routing_gate_pass"]) for row in routes)
    if n_orientation != int(seal["n_orientation_pair_gate_pass"]):
        raise RuntimeError("Orientation gate count mismatch")
    if n_routing != int(seal["n_routing_pair_gate_pass"]):
        raise RuntimeError("Routing gate count mismatch")
    if int(gates[0]["n_routing_pair_pass"]) != n_routing:
        raise RuntimeError("Routing gate table count mismatch")
    if len({row["pair_family_uid"] for row in routes}) != int(seal["n_pair_families"]):
        raise RuntimeError("Routing pair-family count mismatch")
    if len({row["coarse_locus_uid"] for row in routes}) != int(seal["n_physical_loci"]):
        raise RuntimeError("Routing physical-locus count mismatch")
    if yes(gates[0]["experimental_targets_frozen"]):
        raise RuntimeError("Routing gate table froze targets")

    print(
        "LINEAGE_ACCESSIBILITY_ROUTING_VALIDATION_PASS "
        f"pairs={len(routes)} pair_families={len({row['pair_family_uid'] for row in routes})} "
        f"orientation_pass={n_orientation} "
        f"routing_pass={n_routing} targets_frozen=false"
    )


if __name__ == "__main__":
    main()

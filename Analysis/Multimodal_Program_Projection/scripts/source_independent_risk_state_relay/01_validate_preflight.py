#!/usr/bin/env python3
"""Independent structural validation for the Plan 45 preliminary preflight."""

from __future__ import annotations

import json

from relay_common import CANDIDATE_ROOT, read_tsv, sha256_file


def main() -> None:
    audit = CANDIDATE_ROOT / "provisional_target_audit.tsv"
    gates = CANDIDATE_ROOT / "source_gate_status.tsv"
    manifest = CANDIDATE_ROOT / "preflight_input_manifest.tsv"
    preliminary = CANDIDATE_ROOT / "PRELIMINARY.json"
    for path in [audit, gates, manifest, preliminary]:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty Plan 45 preflight artifact: {path}")

    rows = read_tsv(audit)
    if len(rows) != 9 or len({row["gene_symbol"] for row in rows}) != 9:
        raise RuntimeError("Plan 45 legacy universe must contain exactly nine unique genes")
    if any(row["target_freeze_status"] != "prohibited_until_corrected_coloc_registry" for row in rows):
        raise RuntimeError("A target was improperly frozen before the corrected coloc gate")
    if any(row["phenotype_stratum"] not in {"direct_masld_mash_diagnosis", "mri_pdff_or_histologic_steatosis"} for row in rows):
        raise RuntimeError("Proxy-trait candidate entered the direct-disease preflight")

    gate_rows = {row["gate_id"]: row for row in read_tsv(gates)}
    if gate_rows["target_selection"]["status"] != "prohibited":
        raise RuntimeError("Target-selection firewall is not closed")
    observed = int(gate_rows["corrected_coloc_portfolio"]["observed"])
    required = int(gate_rows["corrected_coloc_portfolio"]["required"])
    if not 0 <= observed <= required == 1100:
        raise RuntimeError("Invalid corrected-coloc progress count")

    payload = json.loads(preliminary.read_text(encoding="utf-8"))
    expected = {
        "provisional_target_audit_sha256": sha256_file(audit),
        "source_gate_status_sha256": sha256_file(gates),
        "preflight_input_manifest_sha256": sha256_file(manifest),
    }
    for field, value in expected.items():
        if payload.get(field) != value:
            raise RuntimeError(f"Preflight hash mismatch: {field}")
    if payload.get("target_list_frozen") is not False:
        raise RuntimeError("PRELIMINARY.json incorrectly claims a frozen target list")
    if payload.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("PRELIMINARY.json indicates premature outcome access")
    print(f"Plan 45 preflight validation passed: 9 candidates; coloc {observed}/1100; firewall closed")


if __name__ == "__main__":
    main()

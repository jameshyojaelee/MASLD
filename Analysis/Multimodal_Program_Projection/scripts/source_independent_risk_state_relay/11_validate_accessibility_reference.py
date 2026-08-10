#!/usr/bin/env python3
"""Independently validate the sealed Plan 45 accessibility reference."""

from __future__ import annotations

import json

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def main() -> None:
    seal_path = CANDIDATE_ROOT / "ACCESSIBILITY_REFERENCE_SEALED.json"
    if not seal_path.is_file():
        raise RuntimeError(f"Missing accessibility seal: {seal_path}")
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "accessibility_routing_reference_targets_not_selected":
        raise RuntimeError("Invalid accessibility release status")
    if seal.get("experimental_targets_frozen") is not False:
        raise RuntimeError("Accessibility reference opened target freeze")

    for name, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Output hash mismatch: {name}")

    manifest = read_tsv(CANDIDATE_ROOT / "accessibility_source_manifest.tsv")
    mappings = read_tsv(CANDIDATE_ROOT / "lineage_accessibility_map.tsv")
    provenance = read_tsv(CANDIDATE_ROOT / "accessibility_provenance_manifest.tsv")
    gate = read_tsv(CANDIDATE_ROOT / "accessibility_gate_status.tsv")
    if len(manifest) != 10 or len(mappings) != 16 or len(gate) != 1:
        raise RuntimeError("Accessibility table cardinality drift")
    if sum(row["mapping_status"] == "exact" for row in mappings) != 9:
        raise RuntimeError("Exact lineage-map count drift")
    if gate[0]["status"] != "pass" or gate[0]["experimental_targets_frozen"] != "false":
        raise RuntimeError("Invalid accessibility gate row")
    if {row["reference_build"] for row in manifest} != {"GRCh38"}:
        raise RuntimeError("Mixed accessibility reference builds")
    if sum(row["cohort"] == "GSE244832" for row in manifest) != 9:
        raise RuntimeError("GSE244832 source count drift")
    if sum(row["cohort"] == "GSE281367" for row in manifest) != 1:
        raise RuntimeError("GSE281367 source count drift")

    for row in manifest:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file():
            raise RuntimeError(f"Missing peak source: {path}")
        if path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Peak size drift: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Peak hash drift: {path}")
        n_rows = int(row["n_peak_rows"])
        n_unique = int(row["n_unique_peak_coordinates"])
        n_duplicate = int(row["n_duplicate_peak_rows"])
        if n_rows <= 0 or n_unique <= 0 or n_rows - n_unique != n_duplicate:
            raise RuntimeError(f"Peak arithmetic invalid: {row['source_id']}")

    for row in provenance:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Provenance source drift: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Provenance hash drift: {path}")

    hepatocyte = [row for row in mappings if row["expression_lineage"] == "Hepatocytes"]
    if len(hepatocyte) != 1 or hepatocyte[0]["mapping_status"] != "exact":
        raise RuntimeError("Hepatocyte mapping is not exact")
    if hepatocyte[0]["hepatocyte_independent_replication_required"] != "true":
        raise RuntimeError("Independent hepatocyte replication is not required")

    print(
        "ACCESSIBILITY_REFERENCE_VALIDATION_PASS "
        f"peak_sources={len(manifest)} expression_lineages={len(mappings)} "
        f"exact_maps={sum(row['mapping_status'] == 'exact' for row in mappings)} "
        "targets_frozen=false"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Validate the frozen Plan 45 RSR-02 orientation worklist."""

from __future__ import annotations

import json
import os
from pathlib import Path

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def main() -> None:
    seal_path = CANDIDATE_ROOT / "ORIENTATION_WORKLIST_SEALED.json"
    worklist_path = CANDIDATE_ROOT / "orientation_worklist.tsv"
    for path in [seal_path, worklist_path]:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing orientation artifact: {path}")
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("Orientation worklist indicates premature outcome access")
    if seal.get("experimental_targets_frozen") is not False:
        raise RuntimeError("Orientation worklist improperly freezes experimental targets")
    if sha256_file(worklist_path) != seal.get("orientation_worklist_sha256"):
        raise RuntimeError("Orientation worklist hash mismatch")

    rows = read_tsv(worklist_path)
    if len(rows) != int(seal["eligible_pair_count"]):
        raise RuntimeError("Orientation worklist row count differs from its seal")
    if len({row["orientation_uid"] for row in rows}) != len(rows):
        raise RuntimeError("Orientation worklist contains duplicate UIDs")
    pair_keys = {
        (row["coarse_locus_uid"], row["ensembl_id"], row["gwas_name"])
        for row in rows
    }
    if len(pair_keys) != len(rows):
        raise RuntimeError("Orientation worklist contains duplicate locus-gene-GWAS pairs")
    locus_count = len({row["coarse_locus_uid"] for row in rows})
    family_count = len({row["pair_family_uid"] for row in rows})
    if locus_count != int(seal["eligible_locus_count"]):
        raise RuntimeError("Orientation worklist locus count differs from its seal")
    if family_count != int(seal["eligible_pair_family_count"]):
        raise RuntimeError("Orientation worklist pair-family count differs from its seal")
    observed_family_counts: dict[str, int] = {}
    for row in rows:
        observed_family_counts[row["pair_family_uid"]] = (
            observed_family_counts.get(row["pair_family_uid"], 0) + 1
        )
    for row in rows:
        if int(row["n_eligible_gwas_pairs_for_locus_gene"]) != observed_family_counts[
            row["pair_family_uid"]
        ]:
            raise RuntimeError(
                f"Pair-family support count mismatch: {row['pair_family_uid']}"
            )
    for row in rows:
        if row["pair_export_status"] != "required_not_run":
            raise RuntimeError(f"Pair export unexpectedly marked complete: {row['orientation_uid']}")
        if row["target_freeze_status"] != "prohibited_until_orientation_and_lineage_gates":
            raise RuntimeError(f"Target firewall open: {row['orientation_uid']}")
        for field in [
            "required_pair_summary",
            "required_variant_posterior",
            "required_allele_audit",
        ]:
            if (CANDIDATE_ROOT / row[field]).exists():
                raise RuntimeError(
                    f"Outcome-like pair export existed before worklist freeze: {row[field]}"
                )

    raw_registry = os.environ.get("PLAN45_REGISTRY_ROOT", "").strip()
    if raw_registry:
        registry = Path(raw_registry)
        registry = registry if registry.is_absolute() else PROJECT_ROOT / registry
        registry = registry.resolve()
        if sha256_file(registry / "REGISTRY_V3_SEALED.json") != seal[
            "source_registry_seal_sha256"
        ]:
            raise RuntimeError("Live registry source differs from the worklist seal")

    print(
        "Plan 45 orientation-worklist validation passed: "
        f"eligible_pairs={len(rows)}; pair_families={family_count}; "
        f"loci={locus_count}; pair exports pending; target firewall closed"
    )


if __name__ == "__main__":
    main()

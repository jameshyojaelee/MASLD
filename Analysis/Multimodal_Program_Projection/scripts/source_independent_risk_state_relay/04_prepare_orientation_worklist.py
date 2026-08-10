#!/usr/bin/env python3
"""Freeze the complete, outcome-blind RSR-02 pair-export worklist.

Registry v3 only identifies eligible gene--GWAS pairs.  It cannot orient them.
This script creates a complete worklist for targeted reruns that must preserve
the selected SuSiE signal pair and variant-level SNP.PP.H4 distribution.
"""

from __future__ import annotations

import datetime as dt
import json
import math
import os
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_json,
    read_tsv,
    sha256_file,
    write_tsv,
)


DIRECT_STRATA = {
    "direct_masld_mash_diagnosis",
    "mri_pdff_or_histologic_steatosis",
}
FIELDS = [
    "orientation_uid",
    "pair_family_uid",
    "coarse_locus_uid",
    "gene_symbol",
    "ensembl_id",
    "gwas_name",
    "trait",
    "tier",
    "phenotype_stratum",
    "chromosome",
    "coarse_top_snp",
    "susie_pp4",
    "abf_pp4",
    "ancestry",
    "regulatory_ancestry_status",
    "source_dependence",
    "gwas_evidence_family",
    "gwas_independence_class",
    "counts_for_replication_breadth",
    "n_eligible_gwas_pairs_for_locus_gene",
    "pair_export_status",
    "required_pair_summary",
    "required_variant_posterior",
    "required_allele_audit",
    "target_freeze_status",
]


def finite(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def source_registry_root() -> Path:
    raw = os.environ.get("PLAN45_REGISTRY_ROOT", "").strip()
    if not raw:
        raise RuntimeError("PLAN45_REGISTRY_ROOT must point to a sealed registry-v3 candidate")
    path = Path(raw)
    path = path if path.is_absolute() else PROJECT_ROOT / path
    path = path.resolve()
    allowed = (
        PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates"
    ).resolve()
    if allowed not in path.parents:
        raise RuntimeError(f"Registry source escapes candidate root: {path}")
    return path


def eligible(row: dict[str, str]) -> bool:
    pp4 = finite(row.get("susie_pp4", ""))
    return (
        row.get("pair_primary_genetic") == "true"
        and row.get("is_representative_gene_in_locus") == "true"
        and row.get("static_class") == "genetic_only"
        and row.get("tier") == "1"
        and row.get("phenotype_stratum") in DIRECT_STRATA
        and row.get("regulatory_ancestry_status") == "ancestry_matched_eur"
        and row.get("method") == "susie"
        and pp4 is not None
        and pp4 > 0.5
        and bool(row.get("chromosome"))
    )


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite orientation candidate: {CANDIDATE_ROOT}")
    registry_root = source_registry_root()
    source_seal = registry_root / "REGISTRY_V3_SEALED.json"
    source_pairs = registry_root / "corrected_gene_gwas_pair_registry_v3.tsv"
    for path in [source_seal, source_pairs]:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing registry-v3 source: {path}")
    seal = json.loads(source_seal.read_text(encoding="utf-8"))
    if seal.get("status") != "registry_v3_frozen_orientation_pending":
        raise RuntimeError("Source registry is not a sealed orientation-pending registry v3")
    expected_hash = seal.get("output_sha256", {}).get(source_pairs.name)
    if not expected_hash or sha256_file(source_pairs) != expected_hash:
        raise RuntimeError("Registry-v3 gene-GWAS pair source does not match its seal")

    all_pairs = read_tsv(source_pairs)
    selected = sorted(
        (row for row in all_pairs if eligible(row)),
        key=lambda row: (row["coarse_locus_uid"], row["gene_symbol"], row["gwas_name"]),
    )
    family_counts: dict[tuple[str, str], int] = {}
    for row in selected:
        key = (row["coarse_locus_uid"], row["ensembl_id"])
        family_counts[key] = family_counts.get(key, 0) + 1
    worklist: list[dict[str, object]] = []
    for index, row in enumerate(selected, start=1):
        uid = f"orientation-{index:04d}-{row['gene_symbol']}-{row['gwas_name']}"
        pair_family_uid = f"{row['coarse_locus_uid']}::{row['ensembl_id']}"
        base = f"pair_exports/{uid}"
        worklist.append(
            {
                "orientation_uid": uid,
                "pair_family_uid": pair_family_uid,
                "coarse_locus_uid": row["coarse_locus_uid"],
                "gene_symbol": row["gene_symbol"],
                "ensembl_id": row["ensembl_id"],
                "gwas_name": row["gwas_name"],
                "trait": row["trait"],
                "tier": row["tier"],
                "phenotype_stratum": row["phenotype_stratum"],
                "chromosome": row["chromosome"],
                "coarse_top_snp": row["top_snp"],
                "susie_pp4": row["susie_pp4"],
                "abf_pp4": row["abf_pp4"],
                "ancestry": row["ancestry"],
                "regulatory_ancestry_status": row["regulatory_ancestry_status"],
                "source_dependence": row["source_dependence"],
                "gwas_evidence_family": row["gwas_evidence_family"],
                "gwas_independence_class": row["gwas_independence_class"],
                "counts_for_replication_breadth": row[
                    "counts_for_replication_breadth"
                ],
                "n_eligible_gwas_pairs_for_locus_gene": family_counts[
                    (row["coarse_locus_uid"], row["ensembl_id"])
                ],
                "pair_export_status": "required_not_run",
                "required_pair_summary": f"{base}/selected_pair_summary.tsv",
                "required_variant_posterior": f"{base}/variant_shared_posterior.tsv",
                "required_allele_audit": f"{base}/allele_harmonization.tsv",
                "target_freeze_status": "prohibited_until_orientation_and_lineage_gates",
            }
        )

    CANDIDATE_ROOT.mkdir(parents=True)
    worklist_path = CANDIDATE_ROOT / "orientation_worklist.tsv"
    write_tsv(worklist_path, worklist, FIELDS)
    payload = {
        "status": "orientation_worklist_frozen_pair_exports_pending",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_registry_root": str(registry_root.relative_to(PROJECT_ROOT)),
        "source_registry_seal_sha256": sha256_file(source_seal),
        "source_pair_registry_sha256": sha256_file(source_pairs),
        "eligible_pair_count": len(worklist),
        "eligible_locus_count": len({row["coarse_locus_uid"] for row in worklist}),
        "eligible_pair_family_count": len({row["pair_family_uid"] for row in worklist}),
        "selection_rule": (
            "all qualifying GWAS pairs for each representative genetic_only locus gene; "
            "Tier 1 direct diagnosis or PDFF/histologic steatosis; "
            "EUR ancestry-matched regulatory evidence; corrected SuSiE PP.H4>0.5; "
            "SuSiE method; resolved chromosome"
        ),
        "pair_export_required": True,
        "scientific_outcomes_inspected": False,
        "experimental_targets_frozen": False,
        "orientation_worklist_sha256": sha256_file(worklist_path),
    }
    atomic_write_json(CANDIDATE_ROOT / "ORIENTATION_WORKLIST_SEALED.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

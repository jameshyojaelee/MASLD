#!/usr/bin/env python3
"""Prepare a non-overwriting Plan 42 CLCC1 error-correction workspace."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

from bridge_common import CANDIDATE_ROOT, PROJECT_ROOT, require_validated_seal, sha256_file, write_tsv, atomic_write_text


OLD = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/risk-state-double-dissociation-2026-08-09"
ROOT = CANDIDATE_ROOT / "clcc1_integrity_correction"


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    require_validated_seal()
    ROOT.mkdir(parents=True, exist_ok=True)
    for name in ("SEALED.json", "SEAL_VALIDATED.json"):
        shutil.copy2(CANDIDATE_ROOT / name, ROOT / name)
    copies = {
        OLD / "source/CLCC1_ST1/41586_2025_10064_MOESM3_ESM.xlsx": ROOT / "source/CLCC1_ST1/41586_2025_10064_MOESM3_ESM.xlsx",
        OLD / "source_gate_status.tsv": ROOT / "source_gate_status.tsv",
        OLD / "frozen_inputs/phenotype_registry__phenotype_registry.tsv": ROOT / "frozen_inputs/phenotype_registry__phenotype_registry.tsv",
        OLD / "preprocessed/huh7_observability_covariates.tsv": ROOT / "preprocessed/huh7_observability_covariates.tsv",
    }
    for source, target in copies.items():
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, target)

    classes = read(CANDIDATE_ROOT / "integrity_correction_v2/corrected_evidence_classes_v2.tsv")
    for row in classes:
        if row.get("static_class") in {"genetic_only", "convergent"}:
            row["driving_gwas"] = row.get("tier12_driving_gwas", "")
            row["driving_trait"] = row.get("tier12_driving_trait", "")
    write_tsv(ROOT / "frozen_inputs/evidence_classes__frozen_evidence_classes.tsv", classes, list(classes[0]))

    source_loci = read(CANDIDATE_ROOT / "corrected_genetic_locus_registry.tsv")
    fields = ["locus_uid", "gene_symbol", "ensembl_id", "driving_gwas", "driving_trait", "top_snp",
              "chromosome", "position", "source_susie_pp4", "source_abf_pp4", "locus_resolution",
              "n_genetic_genes_in_locus", "is_representative", "representative_gene"]
    loci = [{
        "locus_uid": row["tier12_locus_uid"], "gene_symbol": row["gene_symbol"], "ensembl_id": row["ensembl_id"],
        "driving_gwas": row["tier12_driving_gwas"], "driving_trait": row["tier12_driving_trait"], "top_snp": row["tier12_top_snp"],
        "chromosome": row["chromosome"], "position": row["position"], "source_susie_pp4": row["tier12_susie_pp4"],
        "source_abf_pp4": row["tier12_abf_pp4"], "locus_resolution": row["locus_resolution"],
        "n_genetic_genes_in_locus": row["n_genetic_genes_in_locus"], "is_representative": row["is_representative"],
        "representative_gene": row["representative_gene"],
    } for row in source_loci]
    if len({row["locus_uid"] for row in loci}) != 326 or sum(row["is_representative"] == "true" for row in loci) != 326:
        raise RuntimeError("Corrected CLCC1 locus registry must contain 326 representatives")
    write_tsv(ROOT / "genetic_locus_registry.tsv", loci, fields)
    manifest = [{"path": str(path.relative_to(ROOT)), "sha256": sha256_file(path)} for path in sorted(ROOT.rglob("*")) if path.is_file()]
    write_tsv(ROOT / "correction_input_manifest.tsv", manifest, ["path", "sha256"])
    atomic_write_text(ROOT / "CORRECTION_SCOPE.json", json.dumps({
        "status": "prepared", "n_corrected_loci": 326,
        "interpretation": "retrospective_error_correction_supplementary_nonconfirmatory_regardless_of_outcome",
        "previous_plan42_release_immutable": True,
    }, indent=2, sort_keys=True) + "\n")
    print("CLCC1_CORRECTION_PREPARED")


if __name__ == "__main__": main()

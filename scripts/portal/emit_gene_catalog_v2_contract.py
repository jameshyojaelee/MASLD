#!/usr/bin/env python3
"""Emit the candidate-only MASLD Gene Catalog v2 public contract.

This writes schema documentation and deterministic experiment-routing rules.
It does not populate gene records or read genetics outcomes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd

import gene_catalog_v2 as catalog


FIELD_CONTRACT: tuple[dict[str, Any], ...] = (
    {
        "field": "catalog_schema_version",
        "type": "string",
        "required_v2": True,
        "allowed_values": catalog.CATALOG_SCHEMA_VERSION,
        "description": "Version of the public MASLD Gene Catalog contract.",
    },
    {
        "field": "ensembl_id",
        "type": "string",
        "required_v2": True,
        "allowed_values": "versioned GENCODE v49 Ensembl gene ID when the object is a transcript",
        "description": "Stable primary identifier; symbols are display aliases only.",
    },
    {
        "field": "gene_biotype",
        "type": "enum",
        "required_v2": True,
        "allowed_values": ";".join(sorted(catalog.GENE_BIOTYPES)),
        "description": "GENCODE transcript biotype or not-applicable state for a DNA element.",
    },
    {
        "field": "candidate_object",
        "type": "enum",
        "required_v2": True,
        "allowed_values": ";".join(sorted(catalog.CANDIDATE_OBJECTS - {"unspecified"})),
        "description": "Biological object nominated by the Resource; determines follow-up logic.",
    },
    {
        "field": "lncrna_genomic_class",
        "type": "enum",
        "required_v2": True,
        "allowed_values": ";".join(sorted(catalog.LNCRNA_GENOMIC_CLASSES)),
        "description": "Deterministic GENCODE-coordinate class for an lncRNA transcript.",
    },
    {
        "field": "annotation_release",
        "type": "string",
        "required_v2": True,
        "allowed_values": "GENCODE v49 for this candidate release",
        "description": "Annotation authority used for identity and genomic classification.",
    },
    {
        "field": "strand_audit_status",
        "type": "enum",
        "required_v2": True,
        "allowed_values": ";".join(sorted(catalog.STRAND_AUDIT_STATES)),
        "description": "Whether the RNA measurement passed the frozen strandedness gate.",
    },
    {
        "field": "mapping_status",
        "type": "enum",
        "required_v2": True,
        "allowed_values": ";".join(sorted(catalog.MAPPING_STATES)),
        "description": "Stable-ID and locus-mapping eligibility for claim-bearing use.",
    },
    {
        "field": "noncoding_dna_context",
        "type": "semicolon-delimited enum",
        "required_v2": True,
        "allowed_values": ";".join(sorted(catalog.NONCODING_DNA_CONTEXTS)),
        "description": "Nonexclusive source-qualified regulatory contexts for the DNA locus.",
    },
    {
        "field": "credible_set_id",
        "type": "string",
        "required_v2": True,
        "allowed_values": "empty or reliable study/trait/locus/ancestry/credible-set identifier",
        "description": "Fine-mapped DNA object; it does not itself establish a target transcript.",
    },
    {
        "field": "credible_set_coding_pip_mass",
        "type": "number or null",
        "required_v2": True,
        "allowed_values": "0..1; all three PIP masses must sum to 1",
        "description": "Credible-set PIP mass assigned to protein-altering coding consequence.",
    },
    {
        "field": "credible_set_noncoding_pip_mass",
        "type": "number or null",
        "required_v2": True,
        "allowed_values": "0..1; all three PIP masses must sum to 1",
        "description": "Credible-set PIP mass assigned to annotated noncoding consequence.",
    },
    {
        "field": "credible_set_unresolved_pip_mass",
        "type": "number or null",
        "required_v2": True,
        "allowed_values": "0..1; all three PIP masses must sum to 1",
        "description": "Credible-set PIP mass lacking a resolved coding/noncoding annotation.",
    },
    {
        "field": "credible_set_to_gene_link_status",
        "type": "enum",
        "required_v2": True,
        "allowed_values": ";".join(sorted(catalog.LINK_STATES)),
        "description": "Strength and type of the link between a fine-mapped signal and transcript.",
    },
    {
        "field": "target_link_basis",
        "type": "enum",
        "required_v2": True,
        "allowed_values": ";".join(sorted(catalog.TARGET_LINK_BASES)),
        "description": "Assay or mapping basis for the proposed target link.",
    },
    {
        "field": "assay_applicability",
        "type": "canonical JSON object",
        "required_v2": True,
        "allowed_values": ";".join(sorted(catalog.ASSAY_APPLICABILITY_STATES)),
        "description": "Per-assay applicability; proteomics is not applicable to lncRNA molecules.",
    },
    {
        "field": "open_mechanistic_question",
        "type": "enum",
        "required_v2": True,
        "allowed_values": ";".join(sorted(catalog.OPEN_QUESTIONS - {"unspecified"})),
        "description": "Unresolved biological question that the next experiment must distinguish.",
    },
    {
        "field": "recommended_experiment_rule_id",
        "type": "derived enum",
        "required_v2": True,
        "allowed_values": ";".join(catalog.EXPERIMENT_RULES_V2),
        "description": "Deterministic outcome-blind route derived from the candidate object and question.",
    },
    {
        "field": "next_experiment_rule_id",
        "type": "derived string",
        "required_v2": True,
        "allowed_values": "same value as recommended_experiment_rule_id during compatibility period",
        "description": "Legacy reader alias; not an independent scientific field.",
    },
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _schema() -> dict[str, Any]:
    properties: dict[str, Any] = {}
    for row in FIELD_CONTRACT:
        properties[row["field"]] = {
            "description": row["description"],
            "contract_type": row["type"],
            "allowed_values": row["allowed_values"],
        }
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "MASLD Gene Catalog v2 extension contract",
        "description": (
            "Candidate-only molecular-object and experiment-routing extension. "
            "Genetics-dependent biological entries are intentionally absent until "
            "the corrected COLOC release is promoted."
        ),
        "catalog_schema_version": catalog.CATALOG_SCHEMA_VERSION,
        "type": "object",
        "required": [row["field"] for row in FIELD_CONTRACT if row["required_v2"]],
        "properties": properties,
        "additionalProperties": True,
    }


def emit(output_dir: Path) -> None:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    output_dir.mkdir(parents=True)

    dictionary_path = output_dir / "masld_gene_catalog_v2_data_dictionary.tsv"
    schema_path = output_dir / "masld_gene_catalog_v2_schema.json"
    rules_path = output_dir / "masld_gene_catalog_v2_experiment_rules.tsv"
    readme_path = output_dir / "README.md"

    pd.DataFrame(FIELD_CONTRACT).to_csv(dictionary_path, sep="\t", index=False)
    with schema_path.open("w", encoding="utf-8") as handle:
        json.dump(_schema(), handle, indent=2, sort_keys=True)
        handle.write("\n")
    pd.DataFrame(
        [
            {"experiment_rule_id": rule_id, **rule}
            for rule_id, rule in catalog.EXPERIMENT_RULES_V2.items()
        ]
    ).to_csv(rules_path, sep="\t", index=False)
    readme_path.write_text(
        "# MASLD Gene Catalog v2 candidate contract\n\n"
        "This candidate contains the public schema, field dictionary, and "
        "deterministic experiment-routing rules. It contains no gene outcomes. "
        "Genetics-dependent entries remain blocked until the corrected COLOC "
        "release is complete and promoted. Regulatory DNA, transcription through "
        "a lncRNA locus, and the mature lncRNA molecule remain separate objects.\n",
        encoding="utf-8",
    )

    artifact_paths = [dictionary_path, schema_path, rules_path, readme_path]
    manifest = {
        "catalog_schema_version": catalog.CATALOG_SCHEMA_VERSION,
        "status": "candidate_contract_only",
        "contains_gene_outcomes": False,
        "genetics_dependent_entries": "blocked_pending_corrected_coloc",
        "artifacts": [
            {
                "path": path.name,
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for path in artifact_paths
        ],
        "producer": str(Path(__file__).resolve()),
        "python": sys.version,
    }
    manifest_path = output_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    emit(args.output_dir.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

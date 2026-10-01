#!/usr/bin/env python3
"""Versioned noncoding-aware contract for the MASLD Gene Catalog.

This module is an adapter beside the sealed v1 ``passport_*`` file contract.
It can read a v1 gene-index record without changing its route, and it provides
the explicit molecular-object routing required by Catalog v2.  It does not
read, select, or score experimental follow-up screen results.
"""

from __future__ import annotations

import json
import math
import re
from collections import OrderedDict
from typing import Any, Mapping


CATALOG_SCHEMA_VERSION = "masld_gene_catalog_v2"
LEGACY_SCHEMA_VERSION = "masld_gene_catalog_v1_compat"

V2_FIELDS = (
    "gene_biotype",
    "candidate_object",
    "lncrna_genomic_class",
    "annotation_release",
    "strand_audit_status",
    "mapping_status",
    "noncoding_dna_context",
    "credible_set_id",
    "credible_set_coding_pip_mass",
    "credible_set_noncoding_pip_mass",
    "credible_set_unresolved_pip_mass",
    "credible_set_to_gene_link_status",
    "target_link_basis",
    "assay_applicability",
    "open_mechanistic_question",
    "recommended_experiment_rule_id",
)

CANDIDATE_OBJECTS = {
    "protein_coding_gene",
    "lncrna_transcript",
    "noncoding_regulatory_element",
    "ambiguous_lncrna_locus",
    "untestable_nomination",
    "unspecified",
}

GENE_BIOTYPES = {"protein_coding", "lncRNA", "other", "not_applicable", "unknown"}

LNCRNA_GENOMIC_CLASSES = {
    "antisense_exonic_overlap",
    "same_strand_exonic_overlap",
    "antisense_gene_body_overlap",
    "same_strand_gene_body_overlap",
    "opposite_strand_promoter_proximal",
    "intergenic",
    "complex",
    "not_applicable",
    "unknown",
}

STRAND_AUDIT_STATES = {
    "passed_reverse_stranded",
    "failed",
    "not_applicable",
    "not_audited",
}

MAPPING_STATES = {
    "unique_stable_ensembl",
    "mapping_ambiguous",
    "not_applicable",
    "unknown",
}

ASSAY_APPLICABILITY_STATES = {
    "applicable",
    "untestable",
    "not_applicable",
    "source_dependent",
    "indeterminate",
}

OPEN_QUESTIONS = {
    "regulatory_element_to_target",
    "lncrna_rna_product",
    "lncrna_dna_transcription_or_rna",
    "protein_coding_state_function",
    "missing_decisive_assay",
    "unspecified",
}

NONCODING_DNA_CONTEXTS = {
    "promoter_proximal",
    "liver_abc_enhancer",
    "liver_lineage_accessible_chromatin",
    "no_deposited_regulatory_context",
    "unresolved",
    "not_applicable",
}

LINK_STATES = {
    "colocalized_egene",
    "source_qualified_context",
    "no_authoritative_link",
    "not_applicable",
    "unknown",
}

TARGET_LINK_BASES = {
    "coloc_susie",
    "coloc_abf",
    "promoter_overlap",
    "abc_enhancer_gene",
    "lineage_accessible_chromatin",
    "multiple_source_qualified",
    "none",
    "unknown",
}


EXPERIMENT_RULES_V2: OrderedDict[str, dict[str, str]] = OrderedDict(
    [
        (
            "EXP_REGULATORY_DNA_V2",
            {
                "biological_model": "human liver-lineage model in which the element is accessible",
                "context": "basal and source-relevant metabolic challenge",
                "perturbation": "allele-aware reporter or base editing of the regulatory element",
                "primary_readout": "proposed target expression plus an assay-native cellular phenotype",
                "falsifying_outcome": "the allele or element does not alter the proposed target in the observable lineage",
            },
        ),
        (
            "EXP_LNCRNA_RNA_PRODUCT_V2",
            {
                "biological_model": "human liver-lineage model expressing the lncRNA",
                "context": "the disease state in which the transcript was observed",
                "perturbation": "RNA depletion using an RNA-targeting nuclease or antisense oligonucleotide, with an independent rescue when feasible",
                "primary_readout": "lncRNA depletion, neighboring-gene expression, and the prespecified cellular phenotype",
                "falsifying_outcome": "adequate RNA depletion changes neither the phenotype nor the proposed downstream readout",
            },
        ),
        (
            "EXP_LNCRNA_LOCUS_DISAMBIGUATION_V2",
            {
                "biological_model": "human liver-lineage model expressing the lncRNA locus",
                "context": "basal and source-relevant disease challenge",
                "perturbation": "paired RNA depletion and locus-level CRISPRi or allele perturbation",
                "primary_readout": "lncRNA abundance, neighboring-gene expression, and the prespecified phenotype",
                "falsifying_outcome": "RNA and locus perturbations do not reproduce the proposed regulatory or phenotypic effect",
            },
        ),
        (
            "EXP_PROTEIN_STATE_CONTEXT_V2",
            {
                "biological_model": "human liver-lineage model matched to the observed disease context",
                "context": "source-relevant metabolic, inflammatory, or fibrotic challenge",
                "perturbation": "gene perturbation with an orthogonal rescue",
                "primary_readout": "the frozen disease-state program and an assay-native cellular phenotype",
                "falsifying_outcome": "perturbation changes neither the disease-state program nor the phenotype",
            },
        ),
        (
            "EXP_MEASURE_MISSING_ASSAY_V2",
            {
                "biological_model": "source-appropriate human tissue or liver-lineage model",
                "context": "the unresolved biological context",
                "perturbation": "acquire the missing decisive measurement before functional interpretation",
                "primary_readout": "prespecified assay applicability, detectability, and native effect estimate",
                "falsifying_outcome": "adequate direct measurement does not support the proposed evidence layer",
            },
        ),
    ]
)


LEGACY_EXPERIMENT_ROUTES = {
    "genetically_anchored": "EXP_ALLELE_AWARE_V1",
    "context_supported": "EXP_CONTEXT_PERTURB_V1",
    "established_state_associated": "EXP_STRESS_PERTURB_V1",
    "concordant": "EXP_MULTI_LAYER_V1",
    "discordant": "EXP_MATCHED_RNA_PROTEIN_V1",
    "untested": "EXP_DIRECT_ASSAY_V1",
    "unresolved": "EXP_UNRESOLVED_V1",
}


_FORBIDDEN_SCREEN_KEYS = {
    "cas13_screen_target",
    "cas13_screen_targets",
    "cas13_guide",
    "cas13_guides",
    "guide_sequence",
    "guide_sequences",
    "cas13_screen_phenotype",
    "cas13_screen_hit",
    "cas13_screen_result",
    "screen_effect",
    "screen_rank",
    "screen_hit",
    "pool_membership",
}

_FORBIDDEN_SCREEN_TEXT = re.compile(
    r"\b(?:cas13\s+screen|screen\s+(?:target|hit|rank|effect)|"
    r"guide\s+sequence|pool\s+membership)\b",
    flags=re.IGNORECASE,
)


class GeneCatalogV2Error(ValueError):
    """Fail-closed v2 contract violation with a stable error code."""

    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


def _clean(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _enum(value: Any, allowed: set[str], field: str) -> str:
    cleaned = _clean(value)
    if cleaned not in allowed:
        raise GeneCatalogV2Error("VOCABULARY", f"{field}={cleaned!r}")
    return cleaned


def _optional_mass(value: Any, field: str) -> float | None:
    if value is None or _clean(value) == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise GeneCatalogV2Error("PIP_MASS", f"{field} is not numeric") from exc
    if not math.isfinite(number) or not 0.0 <= number <= 1.0:
        raise GeneCatalogV2Error("PIP_MASS", f"{field} must be within [0, 1]")
    return number


def canonical_assay_applicability(value: Any) -> str:
    """Return a stable JSON assay-to-applicability map."""

    if value is None or _clean(value) == "":
        return "{}"
    if isinstance(value, str):
        try:
            parsed = json.loads(value)
        except json.JSONDecodeError as exc:
            raise GeneCatalogV2Error("ASSAY_APPLICABILITY_JSON", str(exc)) from exc
    else:
        parsed = value
    if not isinstance(parsed, Mapping):
        raise GeneCatalogV2Error("ASSAY_APPLICABILITY_TYPE", "expected a JSON object")
    normalized: dict[str, str] = {}
    for assay, state in parsed.items():
        assay_name = _clean(assay)
        if not assay_name:
            raise GeneCatalogV2Error("ASSAY_APPLICABILITY_ASSAY", "assay name is empty")
        normalized[assay_name] = _enum(
            state, ASSAY_APPLICABILITY_STATES, f"assay_applicability[{assay_name}]"
        )
    return json.dumps(normalized, sort_keys=True, separators=(",", ":"))


def assert_resource_firewall(record: Mapping[str, Any]) -> None:
    """Reject future experimental screen contents while allowing generic RNA routing."""

    def inspect(value: Any, path: str) -> None:
        if isinstance(value, Mapping):
            forbidden_keys = sorted(
                str(key)
                for key in value
                if str(key).strip().lower() in _FORBIDDEN_SCREEN_KEYS
            )
            if forbidden_keys:
                raise GeneCatalogV2Error(
                    "CAS13_SCREEN_FIREWALL",
                    f"prohibited screen fields at {path}: {forbidden_keys}",
                )
            for key, nested in value.items():
                inspect(nested, f"{path}.{key}")
        elif isinstance(value, (list, tuple, set)):
            for index, nested in enumerate(value):
                inspect(nested, f"{path}[{index}]")
        elif isinstance(value, str) and _FORBIDDEN_SCREEN_TEXT.search(value):
            raise GeneCatalogV2Error(
                "CAS13_SCREEN_FIREWALL", f"prohibited screen content at {path}"
            )

    inspect(record, "record")


def _noncoding_context(value: Any) -> str:
    cleaned = _clean(value) or "not_applicable"
    contexts = [token.strip() for token in cleaned.split(";") if token.strip()]
    unknown = sorted(set(contexts) - NONCODING_DNA_CONTEXTS)
    if unknown:
        raise GeneCatalogV2Error("NONCODING_CONTEXT", f"unknown contexts: {unknown}")
    if "not_applicable" in contexts and len(contexts) > 1:
        raise GeneCatalogV2Error(
            "NONCODING_CONTEXT",
            "not_applicable cannot be combined with another context",
        )
    return ";".join(sorted(set(contexts)))


def route_experiment(record: Mapping[str, Any]) -> str:
    """Route one Catalog entry without using outcomes or screen information."""

    assert_resource_firewall(record)
    candidate_object = _enum(
        record.get("candidate_object", "unspecified"),
        CANDIDATE_OBJECTS,
        "candidate_object",
    )
    question = _enum(
        record.get("open_mechanistic_question", "unspecified"),
        OPEN_QUESTIONS,
        "open_mechanistic_question",
    )
    mapping_status = _enum(
        record.get("mapping_status", "unknown"), MAPPING_STATES, "mapping_status"
    )
    link_status = _enum(
        record.get("credible_set_to_gene_link_status", "unknown"),
        LINK_STATES,
        "credible_set_to_gene_link_status",
    )
    evidence_class = _clean(record.get("primary_evidence_class"))

    if (
        candidate_object == "untestable_nomination"
        or question == "missing_decisive_assay"
        or mapping_status == "mapping_ambiguous"
        and candidate_object not in {"ambiguous_lncrna_locus", "lncrna_transcript"}
    ):
        return "EXP_MEASURE_MISSING_ASSAY_V2"
    if candidate_object == "noncoding_regulatory_element":
        return "EXP_REGULATORY_DNA_V2"
    if candidate_object == "ambiguous_lncrna_locus":
        return "EXP_LNCRNA_LOCUS_DISAMBIGUATION_V2"
    if candidate_object == "lncrna_transcript":
        requires_locus_disambiguation = (
            question == "lncrna_dna_transcription_or_rna"
            or link_status == "colocalized_egene"
            or evidence_class
            in {"genetically_anchored", "genetic_only", "concordant", "convergent"}
        )
        return (
            "EXP_LNCRNA_LOCUS_DISAMBIGUATION_V2"
            if requires_locus_disambiguation
            else "EXP_LNCRNA_RNA_PRODUCT_V2"
        )
    if candidate_object == "protein_coding_gene":
        return "EXP_PROTEIN_STATE_CONTEXT_V2"
    if candidate_object == "unspecified":
        existing = _clean(record.get("next_experiment_rule_id"))
        if existing:
            return existing
        evidence_class = _clean(record.get("primary_evidence_class")) or "unresolved"
        try:
            return LEGACY_EXPERIMENT_ROUTES[evidence_class]
        except KeyError as exc:
            raise GeneCatalogV2Error("LEGACY_ROUTE", evidence_class) from exc
    raise GeneCatalogV2Error("ROUTE_UNDEFINED", candidate_object)


def adapt_record(record: Mapping[str, Any]) -> dict[str, Any]:
    """Read a v1 or v2 gene record and emit a complete compatibility record."""

    assert_resource_firewall(record)
    output = dict(record)
    is_v2 = all(field in record for field in V2_FIELDS[:-1])
    if not is_v2:
        defaults = {
            "gene_biotype": "unknown",
            "candidate_object": "unspecified",
            "lncrna_genomic_class": "unknown",
            "annotation_release": "",
            "strand_audit_status": "not_audited",
            "mapping_status": "unknown",
            "noncoding_dna_context": "not_applicable",
            "credible_set_id": "",
            "credible_set_coding_pip_mass": None,
            "credible_set_noncoding_pip_mass": None,
            "credible_set_unresolved_pip_mass": None,
            "credible_set_to_gene_link_status": "unknown",
            "target_link_basis": "unknown",
            "assay_applicability": "{}",
            "open_mechanistic_question": "unspecified",
        }
        for field, value in defaults.items():
            output.setdefault(field, value)

    output["catalog_schema_version"] = (
        CATALOG_SCHEMA_VERSION if is_v2 else LEGACY_SCHEMA_VERSION
    )
    output["gene_biotype"] = _enum(
        output["gene_biotype"], GENE_BIOTYPES, "gene_biotype"
    )
    output["candidate_object"] = _enum(
        output["candidate_object"], CANDIDATE_OBJECTS, "candidate_object"
    )
    output["lncrna_genomic_class"] = _enum(
        output["lncrna_genomic_class"], LNCRNA_GENOMIC_CLASSES, "lncrna_genomic_class"
    )
    output["strand_audit_status"] = _enum(
        output["strand_audit_status"], STRAND_AUDIT_STATES, "strand_audit_status"
    )
    output["mapping_status"] = _enum(
        output["mapping_status"], MAPPING_STATES, "mapping_status"
    )
    output["noncoding_dna_context"] = _noncoding_context(
        output["noncoding_dna_context"]
    )
    output["credible_set_to_gene_link_status"] = _enum(
        output["credible_set_to_gene_link_status"],
        LINK_STATES,
        "credible_set_to_gene_link_status",
    )
    output["target_link_basis"] = _enum(
        output["target_link_basis"], TARGET_LINK_BASES, "target_link_basis"
    )
    output["assay_applicability"] = canonical_assay_applicability(
        output["assay_applicability"]
    )
    output["open_mechanistic_question"] = _enum(
        output["open_mechanistic_question"], OPEN_QUESTIONS, "open_mechanistic_question"
    )
    if is_v2:
        if not _clean(output["annotation_release"]):
            raise GeneCatalogV2Error(
                "ANNOTATION_RELEASE", "v2 records require annotation_release"
            )
        expected_questions = {
            "protein_coding_gene": {"protein_coding_state_function"},
            "lncrna_transcript": {
                "lncrna_rna_product",
                "lncrna_dna_transcription_or_rna",
            },
            "noncoding_regulatory_element": {"regulatory_element_to_target"},
            "ambiguous_lncrna_locus": {"lncrna_dna_transcription_or_rna"},
            "untestable_nomination": {"missing_decisive_assay"},
        }
        if output["candidate_object"] == "unspecified":
            raise GeneCatalogV2Error("CANDIDATE_OBJECT", "unspecified is legacy-only")
        if (
            output["open_mechanistic_question"]
            not in expected_questions[output["candidate_object"]]
        ):
            raise GeneCatalogV2Error(
                "OBJECT_QUESTION",
                f"{output['candidate_object']} is incompatible with "
                f"{output['open_mechanistic_question']}",
            )
        if output["candidate_object"] == "lncrna_transcript":
            genetically_anchored = output[
                "credible_set_to_gene_link_status"
            ] == "colocalized_egene" or _clean(
                output.get("primary_evidence_class")
            ) in {"genetically_anchored", "genetic_only", "concordant", "convergent"}
            if (
                genetically_anchored
                and output["open_mechanistic_question"]
                != "lncrna_dna_transcription_or_rna"
            ):
                raise GeneCatalogV2Error(
                    "LNCRNA_MECHANISM",
                    "a genetically anchored lncRNA must distinguish DNA, transcription, and RNA-product effects",
                )
    for field in (
        "credible_set_coding_pip_mass",
        "credible_set_noncoding_pip_mass",
        "credible_set_unresolved_pip_mass",
    ):
        output[field] = _optional_mass(output[field], field)

    masses = [
        output["credible_set_coding_pip_mass"],
        output["credible_set_noncoding_pip_mass"],
        output["credible_set_unresolved_pip_mass"],
    ]
    if any(value is not None for value in masses):
        if any(value is None for value in masses):
            raise GeneCatalogV2Error(
                "PIP_MASS_COMPLETENESS", "all three PIP masses are required"
            )
        if not math.isclose(sum(masses), 1.0, abs_tol=1e-6):
            raise GeneCatalogV2Error("PIP_MASS_SUM", f"observed sum={sum(masses):.9g}")
        if not _clean(output["credible_set_id"]):
            raise GeneCatalogV2Error(
                "CREDIBLE_SET_ID", "PIP masses require credible_set_id"
            )

    biotype = output["gene_biotype"]
    candidate_object = output["candidate_object"]
    if (
        candidate_object in {"lncrna_transcript", "ambiguous_lncrna_locus"}
        and biotype != "lncRNA"
    ):
        raise GeneCatalogV2Error(
            "OBJECT_BIOTYPE", f"{candidate_object} requires lncRNA"
        )
    if candidate_object == "protein_coding_gene" and biotype != "protein_coding":
        raise GeneCatalogV2Error(
            "OBJECT_BIOTYPE", "protein_coding_gene requires protein_coding"
        )
    if biotype == "lncRNA":
        if output["lncrna_genomic_class"] in {"not_applicable", "unknown"}:
            raise GeneCatalogV2Error(
                "LNCRNA_GENOMIC_CLASS", "v2 lncRNA records require a genomic class"
            )
        applicability = json.loads(output["assay_applicability"])
        if applicability.get("proteomics") not in {None, "not_applicable"}:
            raise GeneCatalogV2Error(
                "LNCRNA_PROTEOMICS",
                "proteomics must be not_applicable for an lncRNA molecule",
            )

    routed = route_experiment(output)
    supplied = _clean(record.get("recommended_experiment_rule_id"))
    if is_v2 and supplied and supplied != routed:
        raise GeneCatalogV2Error(
            "ROUTE_MUTATION", f"supplied {supplied}; deterministic route is {routed}"
        )
    output["recommended_experiment_rule_id"] = routed
    # Keep the old reader field populated during the compatibility period.
    output["next_experiment_rule_id"] = routed
    assert_resource_firewall(output)
    return output


def experiment_for_record(record: Mapping[str, Any]) -> dict[str, str]:
    """Return the v2 experiment card for a validated v2 record."""

    adapted = adapt_record(record)
    rule_id = adapted["recommended_experiment_rule_id"]
    if rule_id not in EXPERIMENT_RULES_V2:
        raise GeneCatalogV2Error(
            "LEGACY_EXPERIMENT_TEXT",
            "v1 routes remain readable but their experiment text belongs to the sealed v1 rulebook",
        )
    return {"experiment_rule_id": rule_id, **EXPERIMENT_RULES_V2[rule_id]}

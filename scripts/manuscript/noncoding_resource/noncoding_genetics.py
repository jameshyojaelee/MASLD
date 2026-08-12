#!/usr/bin/env python3
"""Fail-closed noncoding genetics adapter for a promoted Resource release."""

from __future__ import annotations

import bisect
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
from typing import Iterable, Mapping, Sequence


PP4_THRESHOLD = 0.5
PIP_SUM_MIN = 0.90
PIP_SUM_MAX = 1.10
MIN_CS_GATE_FRACTION = 0.80
MIN_ANNOTATED_PIP_MASS = 0.95
PROMOTER_WINDOW_BP = 2_000
EXPECTED_COLOC_TASKS = 1_100

INPUT_ROLES = (
    "credible_sets",
    "consequence_annotation",
    "trait_registry",
    "gene_identity",
    "coloc",
    "evidence_classes",
    "variant_liftover",
    "abc_context",
    "atac_context",
)

CS_REQUIRED = (
    "chromosome",
    "position",
    "allele1",
    "allele2",
    "trait",
    "locus",
    "study",
    "ancestry",
    "susie_pip",
    "susie_cs",
    "susie_converged",
    "susie_reliable",
    "variant_id",
)
CONSEQUENCE_REQUIRED = ("variant_id", "Consequence", "class")
TRAIT_REQUIRED = ("study_name", "trait", "tier", "tier_label", "placement")
IDENTITY_REQUIRED = (
    "annotation_release",
    "gene_id_versioned",
    "gene_id_base",
    "gene_name",
    "gene_type",
    "chromosome",
    "start_1based",
    "end_1based",
    "strand",
    "is_canonical_chromosome",
    "base_id_mapping_status",
)
COLOC_REQUIRED = (
    "gwas_name",
    "gene",
    "ensembl",
    "PP.H4.abf",
    "PP.H4.susie",
    "ancestry",
)
EVIDENCE_REQUIRED = (
    "gene_id_versioned",
    "primary_evidence_class",
    "joint_testable",
)
LIFTOVER_REQUIRED = (
    "chromosome",
    "position",
    "allele1",
    "allele2",
    "chr_hg38",
    "pos_hg38",
    "overlaps_any_peak",
    "cell_types_overlapping",
)
ABC_REQUIRED = ("variant_id", "abc_is_self_promoter")
ATAC_REQUIRED = ("variant_id", "cell_type")
PROMOTION_REQUIRED = (
    "release_id",
    "input_role",
    "source_path",
    "source_sha256",
    "promotion_status",
    "n_expected_tasks",
    "n_terminal_tasks",
)

MEMBER_FIELDS = (
    "release_id",
    "credible_set_uid",
    "study",
    "trait",
    "trait_scope",
    "tier",
    "ancestry",
    "locus",
    "susie_cs",
    "variant_id",
    "chromosome",
    "position",
    "allele1",
    "allele2",
    "susie_pip",
    "raw_cs_pip_sum",
    "pip_sum_gate_passed",
    "normalized_susie_pip",
    "consequence",
    "consequence_category",
    "consequence_annotation_status",
    "promoter_proximal",
    "promoter_context_status",
    "promoter_gene_ids",
    "abc_enhancer_overlap",
    "lineage_accessible",
    "lineage_accessibility_status",
    "lineage_cell_types",
)

ARCHITECTURE_FIELDS = (
    "release_id",
    "credible_set_uid",
    "study",
    "trait",
    "trait_scope",
    "tier",
    "ancestry",
    "locus",
    "susie_cs",
    "n_members",
    "raw_cs_pip_sum",
    "pip_sum_gate_passed",
    "protein_altering_pip_mass",
    "canonical_splice_pip_mass",
    "synonymous_or_utr_pip_mass",
    "other_noncoding_pip_mass",
    "unresolved_pip_mass",
    "annotated_pip_mass",
    "mass_sum_check",
)

REGULATORY_FIELDS = (
    "release_id",
    "credible_set_uid",
    "study",
    "trait",
    "trait_scope",
    "tier",
    "ancestry",
    "locus",
    "susie_cs",
    "pip_sum_gate_passed",
    "promoter_proximal_pip_mass",
    "abc_enhancer_pip_mass",
    "lineage_accessible_pip_mass",
    "no_deposited_context_pip_mass",
    "unresolved_context_pip_mass",
    "lineage_cell_types",
)

COLOC_BIOTYPE_FIELDS = (
    "release_id",
    "method",
    "trait_scope",
    "tier",
    "gwas_name",
    "ancestry",
    "gene_id_versioned",
    "gene_id_base",
    "gencode_gene_name",
    "source_gene_symbol",
    "gene_biotype",
    "pp_h4",
    "threshold",
)

EVIDENCE_BIOTYPE_FIELDS = (
    "release_id",
    "primary_evidence_class",
    "joint_testable",
    "gene_biotype",
    "n_genes",
    "fraction_within_evidence_class",
)

LINK_STATUS_FIELDS = (
    "release_id",
    "method",
    "trait_scope",
    "gwas_name",
    "gene_id_versioned",
    "gene_biotype",
    "credible_set_id",
    "credible_set_to_gene_link_status",
    "target_link_basis",
    "interpretation",
)

VERDICT_FIELDS = (
    "release_id",
    "record_type",
    "trait_scope",
    "method",
    "metric",
    "value",
    "threshold",
    "passed",
    "detail",
)


class NoncodingGeneticsError(RuntimeError):
    """Raised when a promoted input or scientific invariant is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise NoncodingGeneticsError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def delimiter_for(path: Path) -> str:
    return "\t" if path.suffix in {".tsv", ".txt"} else ","


def read_table(path: Path, required: Sequence[str]) -> list[dict[str, str]]:
    require(path.is_file(), f"missing input: {path}")
    require(not path.is_symlink(), f"symlinked input is prohibited: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter_for(path))
        fields = tuple(reader.fieldnames or ())
        missing = [field for field in required if field not in fields]
        require(not missing, f"{path} is missing required columns: {missing}")
        return list(reader)


def write_tsv(
    path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            require(set(row) == set(fields), f"schema drift while writing {path.name}")
            writer.writerow(row)


def parse_bool(value: str, label: str) -> bool:
    normalized = value.strip().lower()
    require(
        normalized in {"true", "false", "1", "0"},
        f"invalid boolean for {label}: {value!r}",
    )
    return normalized in {"true", "1"}


def parse_float(value: str, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise NoncodingGeneticsError(
            f"invalid numeric value for {label}: {value!r}"
        ) from error
    require(math.isfinite(number), f"non-finite numeric value for {label}: {value!r}")
    return number


def parse_optional_float(value: str, label: str) -> float | None:
    if value.strip() in {"", "NA", "NaN", "nan"}:
        return None
    return parse_float(value, label)


def parse_int(value: str, label: str) -> int:
    number = parse_float(value, label)
    require(number.is_integer(), f"non-integral value for {label}: {value!r}")
    return int(number)


def variant_id(row: Mapping[str, str]) -> str:
    return ":".join(
        (row["chromosome"], row["position"], row["allele1"], row["allele2"])
    )


def stable_uid(parts: Iterable[str]) -> str:
    return "cs_" + hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:20]


def trait_scope(row: Mapping[str, str]) -> str:
    tier = parse_int(row["tier"], "trait tier")
    label = row["tier_label"]
    if tier == 1 and label == "direct_MASLD":
        return "tier1_direct_masld_pdff"
    if tier == 2 and label == "liver_enzyme":
        return "tier2_liver_enzyme"
    raise NoncodingGeneticsError(
        f"main Tier-1/2 trait has unexpected tier/label pair: {tier}/{label!r}"
    )


def validate_promotion_manifest(
    manifest_path: Path,
    inputs: Mapping[str, Path],
) -> tuple[str, list[dict[str, str]]]:
    rows = read_table(manifest_path, PROMOTION_REQUIRED)
    by_role: dict[str, dict[str, str]] = {}
    for row in rows:
        role = row["input_role"]
        require(role not in by_role, f"duplicate promotion-manifest role: {role}")
        by_role[role] = row
    require(
        set(by_role) == set(INPUT_ROLES),
        "promotion manifest input-role family is incomplete or expanded",
    )
    release_ids = {row["release_id"] for row in rows}
    require(
        len(release_ids) == 1 and "" not in release_ids,
        "promotion manifest must have one release_id",
    )
    for role in INPUT_ROLES:
        row = by_role[role]
        source = inputs[role].resolve()
        require(row["promotion_status"] == "promoted", f"{role} is not promoted")
        require(
            Path(row["source_path"]).resolve() == source,
            f"{role} path differs from promotion manifest",
        )
        require(
            row["source_sha256"] == sha256_file(source),
            f"{role} SHA256 differs from promotion manifest",
        )
    coloc_gate = by_role["coloc"]
    require(
        parse_int(coloc_gate["n_expected_tasks"], "coloc n_expected_tasks")
        == EXPECTED_COLOC_TASKS,
        f"corrected COLOC expected-task count must be {EXPECTED_COLOC_TASKS}",
    )
    require(
        parse_int(coloc_gate["n_terminal_tasks"], "coloc n_terminal_tasks")
        == EXPECTED_COLOC_TASKS,
        "all corrected COLOC tasks must be terminal before this adapter runs",
    )
    return next(iter(release_ids)), rows


def load_trait_registry(path: Path) -> dict[str, dict[str, str]]:
    rows = read_table(path, TRAIT_REQUIRED)
    registry: dict[str, dict[str, str]] = {}
    for row in rows:
        study = row["study_name"]
        require(study not in registry, f"duplicate trait-registry study: {study}")
        registry[study] = row
    return registry


def consequence_category(source_class: str, consequence: str) -> str:
    tokens = {token.strip() for token in consequence.split(",") if token.strip()}
    if source_class == "coding_protein_altering":
        return "protein_altering"
    if tokens & {"splice_acceptor_variant", "splice_donor_variant"}:
        return "canonical_splice"
    if source_class == "coding_synonymous" or tokens & {
        "synonymous_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
    }:
        return "synonymous_or_utr"
    if source_class or tokens:
        return "other_noncoding"
    return "unresolved"


def load_unique_by_variant(
    path: Path,
    required: Sequence[str],
) -> dict[str, dict[str, str]]:
    output: dict[str, dict[str, str]] = {}
    for row in read_table(path, required):
        key = row.get("variant_id", "") or variant_id(row)
        require(key not in output, f"duplicate variant in {path.name}: {key}")
        output[key] = row
    return output


def build_identity_maps(
    rows: Sequence[Mapping[str, str]],
) -> tuple[
    dict[str, Mapping[str, str]],
    dict[str, Mapping[str, str]],
    dict[str, list[tuple[int, str]]],
]:
    versioned: dict[str, Mapping[str, str]] = {}
    base: dict[str, Mapping[str, str]] = {}
    tss_by_chromosome: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for row in rows:
        require(
            row["annotation_release"] == "GENCODE v49",
            "gene identity input must use GENCODE v49",
        )
        gene_id = row["gene_id_versioned"]
        require(gene_id not in versioned, f"duplicate versioned GENCODE ID: {gene_id}")
        versioned[gene_id] = row
        if row["base_id_mapping_status"] == "unique":
            base_id = row["gene_id_base"]
            require(
                base_id not in base, f"duplicate supposedly unique base ID: {base_id}"
            )
            base[base_id] = row
        if parse_bool(row["is_canonical_chromosome"], f"canonical flag {gene_id}"):
            start = parse_int(row["start_1based"], f"start {gene_id}")
            end = parse_int(row["end_1based"], f"end {gene_id}")
            require(
                row["strand"] in {"+", "-"}, f"invalid GENCODE strand for {gene_id}"
            )
            tss_by_chromosome[row["chromosome"]].append(
                (start if row["strand"] == "+" else end, gene_id)
            )
    for values in tss_by_chromosome.values():
        values.sort()
    return versioned, base, tss_by_chromosome


def promoter_gene_ids(
    chromosome: str,
    position: int,
    tss_by_chromosome: Mapping[str, Sequence[tuple[int, str]]],
) -> list[str]:
    indexed = tss_by_chromosome.get(chromosome, ())
    left = bisect.bisect_left(indexed, (position - PROMOTER_WINDOW_BP, ""))
    right = bisect.bisect_right(indexed, (position + PROMOTER_WINDOW_BP, "\U0010ffff"))
    return sorted(gene_id for _, gene_id in indexed[left:right])


def format_number(value: float | None) -> str:
    return "" if value is None else format(value, ".12g")


def normalized_biotype(gene_type: str) -> str:
    if gene_type in {"protein_coding", "lncRNA"}:
        return gene_type
    return "other"


def build_credible_set_products(
    release_id: str,
    cs_rows: Sequence[Mapping[str, str]],
    trait_registry: Mapping[str, Mapping[str, str]],
    consequence_by_variant: Mapping[str, Mapping[str, str]],
    liftover_by_variant: Mapping[str, Mapping[str, str]],
    abc_enhancer_variants: set[str],
    atac_cell_types: Mapping[str, set[str]],
    tss_by_chromosome: Mapping[str, Sequence[tuple[int, str]]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    groups: dict[tuple[str, str, str, str, int], list[Mapping[str, str]]] = defaultdict(
        list
    )
    seen_members: set[tuple[str, str, str, str, int, str]] = set()
    for row in cs_rows:
        study = row["study"]
        registry = trait_registry.get(study)
        require(
            registry is not None,
            f"credible-set study absent from trait registry: {study}",
        )
        if registry["placement"] != "main" or parse_int(
            registry["tier"], "tier"
        ) not in {1, 2}:
            continue
        if not parse_bool(row["susie_converged"], f"susie_converged {study}"):
            continue
        if not parse_bool(row["susie_reliable"], f"susie_reliable {study}"):
            continue
        cs_id = parse_int(row["susie_cs"], f"susie_cs {study}")
        if cs_id <= 0:
            continue
        require(row["trait"] == registry["trait"], f"trait mismatch for {study}")
        pip = parse_float(row["susie_pip"], f"susie_pip {study}")
        require(0 <= pip <= 1, f"SuSiE PIP outside [0,1] for {study}")
        require(
            row["variant_id"] == variant_id(row),
            f"variant_id fields disagree for {study}",
        )
        key = (study, row["trait"], row["ancestry"], row["locus"], cs_id)
        member_key = key + (row["variant_id"],)
        require(
            member_key not in seen_members,
            f"duplicate credible-set member: {member_key}",
        )
        seen_members.add(member_key)
        groups[key].append(row)
    require(groups, "no reliable main Tier-1/2 SuSiE credible sets")

    member_output: list[dict[str, object]] = []
    architecture_output: list[dict[str, object]] = []
    context_output: list[dict[str, object]] = []
    for key in sorted(groups):
        study, trait, ancestry, locus, cs_id = key
        registry = trait_registry[study]
        scope = trait_scope(registry)
        uid = stable_uid((study, trait, ancestry, locus, str(cs_id)))
        members = sorted(groups[key], key=lambda row: row["variant_id"])
        raw_sum = sum(parse_float(row["susie_pip"], "susie_pip") for row in members)
        gate = PIP_SUM_MIN <= raw_sum <= PIP_SUM_MAX
        category_mass = Counter()
        context_mass = Counter()
        all_cell_types: set[str] = set()
        for row in members:
            vid = row["variant_id"]
            pip = parse_float(row["susie_pip"], f"susie_pip {vid}")
            normalized = pip / raw_sum if gate else None
            annotation = consequence_by_variant.get(vid)
            if annotation is None:
                consequence = ""
                category = "unresolved"
                annotation_status = "unresolved_missing_variant_annotation"
            else:
                consequence = annotation["Consequence"]
                category = consequence_category(annotation["class"], consequence)
                annotation_status = "annotated"
            if normalized is not None:
                category_mass[category] += normalized
                if vid in abc_enhancer_variants:
                    context_mass["abc"] += normalized

            lift = liftover_by_variant.get(vid)
            if lift is None or not lift["chr_hg38"] or not lift["pos_hg38"]:
                promoter_hits: list[str] = []
                promoter_value = ""
                promoter_status = "unresolved_liftover"
                lineage_value = ""
                lineage_status = "unresolved_liftover"
                cell_types: set[str] = set()
                if normalized is not None:
                    context_mass["unresolved"] += normalized
            else:
                promoter_hits = promoter_gene_ids(
                    lift["chr_hg38"],
                    parse_int(lift["pos_hg38"], f"pos_hg38 {vid}"),
                    tss_by_chromosome,
                )
                promoter_value = str(bool(promoter_hits)).lower()
                promoter_status = "resolved"
                summary_accessible = parse_bool(
                    lift["overlaps_any_peak"], f"overlaps_any_peak {vid}"
                )
                cell_types = set(atac_cell_types.get(vid, set()))
                require(
                    summary_accessible == bool(cell_types),
                    f"ATAC summary/detail disagreement for {vid}",
                )
                lineage_value = str(summary_accessible).lower()
                lineage_status = "resolved"
                all_cell_types.update(cell_types)
                if normalized is not None:
                    if promoter_hits:
                        context_mass["promoter"] += normalized
                    if summary_accessible:
                        context_mass["accessible"] += normalized
                    if (
                        not promoter_hits
                        and vid not in abc_enhancer_variants
                        and not summary_accessible
                    ):
                        context_mass["none"] += normalized

            member_output.append(
                {
                    "release_id": release_id,
                    "credible_set_uid": uid,
                    "study": study,
                    "trait": trait,
                    "trait_scope": scope,
                    "tier": registry["tier"],
                    "ancestry": ancestry,
                    "locus": locus,
                    "susie_cs": str(cs_id),
                    "variant_id": vid,
                    "chromosome": row["chromosome"],
                    "position": row["position"],
                    "allele1": row["allele1"],
                    "allele2": row["allele2"],
                    "susie_pip": format_number(pip),
                    "raw_cs_pip_sum": format_number(raw_sum),
                    "pip_sum_gate_passed": str(gate).lower(),
                    "normalized_susie_pip": format_number(normalized),
                    "consequence": consequence,
                    "consequence_category": category,
                    "consequence_annotation_status": annotation_status,
                    "promoter_proximal": promoter_value,
                    "promoter_context_status": promoter_status,
                    "promoter_gene_ids": ";".join(promoter_hits),
                    "abc_enhancer_overlap": str(vid in abc_enhancer_variants).lower(),
                    "lineage_accessible": lineage_value,
                    "lineage_accessibility_status": lineage_status,
                    "lineage_cell_types": ";".join(sorted(cell_types)),
                }
            )
        masses = {
            name: category_mass.get(name, 0.0) if gate else None
            for name in (
                "protein_altering",
                "canonical_splice",
                "synonymous_or_utr",
                "other_noncoding",
                "unresolved",
            )
        }
        mass_sum = (
            sum(value for value in masses.values() if value is not None)
            if gate
            else None
        )
        if gate:
            require(
                abs((mass_sum or 0) - 1) <= 1e-9,
                f"consequence masses do not sum to one for {uid}",
            )
        architecture_output.append(
            {
                "release_id": release_id,
                "credible_set_uid": uid,
                "study": study,
                "trait": trait,
                "trait_scope": scope,
                "tier": registry["tier"],
                "ancestry": ancestry,
                "locus": locus,
                "susie_cs": str(cs_id),
                "n_members": str(len(members)),
                "raw_cs_pip_sum": format_number(raw_sum),
                "pip_sum_gate_passed": str(gate).lower(),
                "protein_altering_pip_mass": format_number(masses["protein_altering"]),
                "canonical_splice_pip_mass": format_number(masses["canonical_splice"]),
                "synonymous_or_utr_pip_mass": format_number(
                    masses["synonymous_or_utr"]
                ),
                "other_noncoding_pip_mass": format_number(masses["other_noncoding"]),
                "unresolved_pip_mass": format_number(masses["unresolved"]),
                "annotated_pip_mass": format_number(
                    None if masses["unresolved"] is None else 1 - masses["unresolved"]
                ),
                "mass_sum_check": format_number(mass_sum),
            }
        )
        context_output.append(
            {
                "release_id": release_id,
                "credible_set_uid": uid,
                "study": study,
                "trait": trait,
                "trait_scope": scope,
                "tier": registry["tier"],
                "ancestry": ancestry,
                "locus": locus,
                "susie_cs": str(cs_id),
                "pip_sum_gate_passed": str(gate).lower(),
                "promoter_proximal_pip_mass": format_number(
                    context_mass.get("promoter", 0.0) if gate else None
                ),
                "abc_enhancer_pip_mass": format_number(
                    context_mass.get("abc", 0.0) if gate else None
                ),
                "lineage_accessible_pip_mass": format_number(
                    context_mass.get("accessible", 0.0) if gate else None
                ),
                "no_deposited_context_pip_mass": format_number(
                    context_mass.get("none", 0.0) if gate else None
                ),
                "unresolved_context_pip_mass": format_number(
                    context_mass.get("unresolved", 0.0) if gate else None
                ),
                "lineage_cell_types": ";".join(sorted(all_cell_types)),
            }
        )
    return member_output, architecture_output, context_output


def build_coloc_products(
    release_id: str,
    coloc_rows: Sequence[Mapping[str, str]],
    trait_registry: Mapping[str, Mapping[str, str]],
    versioned_identity: Mapping[str, Mapping[str, str]],
    base_identity: Mapping[str, Mapping[str, str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    targets: list[dict[str, object]] = []
    links: list[dict[str, object]] = []
    summary_counts: dict[tuple[str, str], Counter] = defaultdict(Counter)
    positive_ids: dict[tuple[str, str], set[str]] = defaultdict(set)
    seen: set[tuple[str, str, str]] = set()
    for row in coloc_rows:
        study = row["gwas_name"]
        registry = trait_registry.get(study)
        if (
            registry is None
            or registry["placement"] != "main"
            or parse_int(registry["tier"], "tier") not in {1, 2}
        ):
            continue
        scope = trait_scope(registry)
        raw_gene_id = row["ensembl"].strip()
        key = (study, row["ancestry"], raw_gene_id)
        require(key not in seen, f"duplicate gene-study COLOC row: {key}")
        seen.add(key)
        identity = versioned_identity.get(raw_gene_id) or base_identity.get(raw_gene_id)
        for method, field in (("abf", "PP.H4.abf"), ("susie", "PP.H4.susie")):
            value = parse_optional_float(row[field], f"{field} {study}/{raw_gene_id}")
            counts = summary_counts[(scope, method)]
            if value is None:
                continue
            require(0 <= value <= 1, f"{field} outside [0,1] for {study}/{raw_gene_id}")
            counts["n_tested"] += 1
            if value <= PP4_THRESHOLD:
                continue
            counts["n_positive"] += 1
            require(
                identity is not None,
                f"positive COLOC target has unmapped GENCODE ID: {raw_gene_id}",
            )
            biotype = normalized_biotype(identity["gene_type"])
            counts[f"n_{biotype}"] += 1
            positive_ids[(scope, method)].add(identity["gene_id_versioned"])
            target = {
                "release_id": release_id,
                "method": method,
                "trait_scope": scope,
                "tier": registry["tier"],
                "gwas_name": study,
                "ancestry": row["ancestry"],
                "gene_id_versioned": identity["gene_id_versioned"],
                "gene_id_base": identity["gene_id_base"],
                "gencode_gene_name": identity["gene_name"],
                "source_gene_symbol": row["gene"],
                "gene_biotype": biotype,
                "pp_h4": format_number(value),
                "threshold": format_number(PP4_THRESHOLD),
            }
            targets.append(target)
            links.append(
                {
                    "release_id": release_id,
                    "method": method,
                    "trait_scope": scope,
                    "gwas_name": study,
                    "gene_id_versioned": identity["gene_id_versioned"],
                    "gene_biotype": biotype,
                    "credible_set_id": "",
                    "credible_set_to_gene_link_status": "not_established_no_explicit_signal_pair",
                    "target_link_basis": "gene_level_colocalization_only",
                    "interpretation": "noncoding_variant_location_does_not_establish_lncrna_or_neighboring_gene_mediation",
                }
            )
    verdict_rows: list[dict[str, object]] = []
    for (scope, method), counts in sorted(summary_counts.items()):
        metrics = {
            metric: counts.get(metric, 0)
            for metric in (
                "n_tested",
                "n_positive",
                "n_protein_coding",
                "n_lncRNA",
                "n_other",
            )
        }
        metrics["n_unique_positive_targets"] = len(positive_ids[(scope, method)])
        for metric, count in metrics.items():
            verdict_rows.append(
                {
                    "release_id": release_id,
                    "record_type": "coloc_method_census",
                    "trait_scope": scope,
                    "method": method,
                    "metric": metric,
                    "value": str(count),
                    "threshold": "",
                    "passed": "",
                    "detail": "PP.H4 methods retain separate tested and positive denominators",
                }
            )
    targets.sort(
        key=lambda row: (
            row["trait_scope"],
            row["method"],
            row["gwas_name"],
            row["gene_id_versioned"],
        )
    )
    links.sort(
        key=lambda row: (
            row["trait_scope"],
            row["method"],
            row["gwas_name"],
            row["gene_id_versioned"],
        )
    )
    return targets, links, verdict_rows


def build_evidence_biotype(
    release_id: str,
    evidence_rows: Sequence[Mapping[str, str]],
    versioned_identity: Mapping[str, Mapping[str, str]],
) -> list[dict[str, object]]:
    counts: Counter[tuple[str, bool, str]] = Counter()
    seen: set[str] = set()
    for row in evidence_rows:
        gene_id = row["gene_id_versioned"]
        require(gene_id not in seen, f"duplicate evidence-class gene: {gene_id}")
        seen.add(gene_id)
        identity = versioned_identity.get(gene_id)
        require(
            identity is not None,
            f"evidence-class gene absent from GENCODE identity: {gene_id}",
        )
        joint = parse_bool(row["joint_testable"], f"joint_testable {gene_id}")
        counts[
            (
                row["primary_evidence_class"],
                joint,
                normalized_biotype(identity["gene_type"]),
            )
        ] += 1
    totals = Counter()
    for (evidence_class, joint, _), count in counts.items():
        totals[(evidence_class, joint)] += count
    output = []
    for (evidence_class, joint, biotype), count in sorted(counts.items()):
        output.append(
            {
                "release_id": release_id,
                "primary_evidence_class": evidence_class,
                "joint_testable": str(joint).lower(),
                "gene_biotype": biotype,
                "n_genes": str(count),
                "fraction_within_evidence_class": format_number(
                    count / totals[(evidence_class, joint)]
                ),
            }
        )
    return output


def gate_verdicts(
    release_id: str,
    architecture: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    n_sets = len(architecture)
    passed_sets = [row for row in architecture if row["pip_sum_gate_passed"] == "true"]
    fraction = len(passed_sets) / n_sets
    total_unresolved = sum(
        parse_float(str(row["unresolved_pip_mass"]), "unresolved mass")
        for row in passed_sets
    )
    annotated_mass = 1 - total_unresolved / len(passed_sets) if passed_sets else 0.0
    coverage_passed = fraction >= MIN_CS_GATE_FRACTION
    annotation_passed = annotated_mass >= MIN_ANNOTATED_PIP_MASS
    gates = (
        (
            "reliable_credible_set_pip_sum_coverage",
            fraction,
            MIN_CS_GATE_FRACTION,
            coverage_passed,
            f"{len(passed_sets)}/{n_sets} reliable credible sets have raw SuSiE PIP sums in [0.90,1.10]",
        ),
        (
            "aggregate_consequence_annotation_coverage",
            annotated_mass,
            MIN_ANNOTATED_PIP_MASS,
            annotation_passed,
            "aggregate normalized PIP mass with a mutually exclusive consequence annotation",
        ),
        (
            "top_snp_credible_set_to_gene_join",
            0.0,
            0.0,
            True,
            "no credible-set-to-gene link is inferred from COLOC top_snp",
        ),
        (
            "main_figure_architecture_eligible",
            float(coverage_passed and annotation_passed),
            1.0,
            coverage_passed and annotation_passed,
            "requires both credible-set PIP-sum coverage and consequence-annotation gates",
        ),
    )
    return [
        {
            "release_id": release_id,
            "record_type": "architecture_gate",
            "trait_scope": "all_tier1_tier2",
            "method": "susie",
            "metric": metric,
            "value": format_number(value),
            "threshold": format_number(threshold),
            "passed": str(passed).lower(),
            "detail": detail,
        }
        for metric, value, threshold, passed, detail in gates
    ]


def build_release(
    *,
    output_dir: Path,
    promotion_manifest: Path,
    credible_sets: Path,
    consequence_annotation: Path,
    trait_registry: Path,
    gene_identity: Path,
    coloc: Path,
    evidence_classes: Path,
    variant_liftover: Path,
    abc_context: Path,
    atac_context: Path,
) -> dict[str, object]:
    require(not output_dir.exists(), f"output directory already exists: {output_dir}")
    inputs = {
        "credible_sets": credible_sets,
        "consequence_annotation": consequence_annotation,
        "trait_registry": trait_registry,
        "gene_identity": gene_identity,
        "coloc": coloc,
        "evidence_classes": evidence_classes,
        "variant_liftover": variant_liftover,
        "abc_context": abc_context,
        "atac_context": atac_context,
    }
    release_id, promotion_rows = validate_promotion_manifest(promotion_manifest, inputs)
    registry = load_trait_registry(trait_registry)
    identity_rows = read_table(gene_identity, IDENTITY_REQUIRED)
    versioned_identity, base_identity, tss_by_chromosome = build_identity_maps(
        identity_rows
    )
    consequences = load_unique_by_variant(consequence_annotation, CONSEQUENCE_REQUIRED)
    liftover = load_unique_by_variant(variant_liftover, LIFTOVER_REQUIRED)
    abc_enhancer_variants = {
        row["variant_id"]
        for row in read_table(abc_context, ABC_REQUIRED)
        if not parse_bool(
            row["abc_is_self_promoter"], f"ABC promoter {row['variant_id']}"
        )
    }
    atac_cell_types: dict[str, set[str]] = defaultdict(set)
    for row in read_table(atac_context, ATAC_REQUIRED):
        require(
            row["cell_type"].strip() != "",
            f"empty ATAC cell type for {row['variant_id']}",
        )
        atac_cell_types[row["variant_id"]].add(row["cell_type"])

    members, architecture, regulatory = build_credible_set_products(
        release_id,
        read_table(credible_sets, CS_REQUIRED),
        registry,
        consequences,
        liftover,
        abc_enhancer_variants,
        atac_cell_types,
        tss_by_chromosome,
    )
    targets, link_status, coloc_verdicts = build_coloc_products(
        release_id,
        read_table(coloc, COLOC_REQUIRED),
        registry,
        versioned_identity,
        base_identity,
    )
    evidence_biotype = build_evidence_biotype(
        release_id,
        read_table(evidence_classes, EVIDENCE_REQUIRED),
        versioned_identity,
    )
    verdicts = gate_verdicts(release_id, architecture) + coloc_verdicts

    output_dir.mkdir(parents=True, exist_ok=False)
    products = (
        ("credible_set_members_annotated.tsv", members, MEMBER_FIELDS),
        ("credible_set_pip_architecture.tsv", architecture, ARCHITECTURE_FIELDS),
        ("credible_set_regulatory_context.tsv", regulatory, REGULATORY_FIELDS),
        ("coloc_target_biotype.tsv", targets, COLOC_BIOTYPE_FIELDS),
        ("evidence_class_by_biotype.tsv", evidence_biotype, EVIDENCE_BIOTYPE_FIELDS),
        ("noncoding_dna_lncrna_link_status.tsv", link_status, LINK_STATUS_FIELDS),
        ("genetics_noncoding_verdict.tsv", verdicts, VERDICT_FIELDS),
    )
    for filename, rows, fields in products:
        write_tsv(output_dir / filename, rows, fields)

    source_rows = []
    for row in sorted(promotion_rows, key=lambda item: item["input_role"]):
        source_rows.append(
            {
                "release_id": release_id,
                "input_role": row["input_role"],
                "source_path": str(Path(row["source_path"]).resolve()),
                "size_bytes": str(Path(row["source_path"]).stat().st_size),
                "sha256": row["source_sha256"],
                "promotion_status": row["promotion_status"],
            }
        )
    source_rows.append(
        {
            "release_id": release_id,
            "input_role": "promotion_manifest",
            "source_path": str(promotion_manifest.resolve()),
            "size_bytes": str(promotion_manifest.stat().st_size),
            "sha256": sha256_file(promotion_manifest),
            "promotion_status": "control_manifest",
        }
    )
    source_rows.sort(key=lambda row: row["input_role"])
    write_tsv(
        output_dir / "source_manifest.tsv",
        source_rows,
        (
            "release_id",
            "input_role",
            "source_path",
            "size_bytes",
            "sha256",
            "promotion_status",
        ),
    )
    execution = {
        "release_id": release_id,
        "producer": str(Path(__file__).resolve()),
        "producer_sha256": sha256_file(Path(__file__).resolve()),
        "python_version": platform.python_version(),
        "n_reliable_credible_sets": len(architecture),
        "n_credible_set_members": len(members),
        "n_colocalized_target_rows": len(targets),
        "pp4_threshold": PP4_THRESHOLD,
        "pip_sum_gate": [PIP_SUM_MIN, PIP_SUM_MAX],
        "top_snp_join_used": False,
    }
    with (output_dir / "execution_manifest.json").open("x", encoding="utf-8") as handle:
        json.dump(execution, handle, indent=2, sort_keys=True)
        handle.write("\n")
    checksum_rows = []
    for path in sorted(output_dir.iterdir(), key=lambda item: item.name):
        if path.name == "checksum_manifest.tsv":
            continue
        checksum_rows.append(
            {
                "relative_path": path.name,
                "size_bytes": str(path.stat().st_size),
                "sha256": sha256_file(path),
            }
        )
    write_tsv(
        output_dir / "checksum_manifest.tsv",
        checksum_rows,
        ("relative_path", "size_bytes", "sha256"),
    )
    return {
        "release_id": release_id,
        "n_reliable_credible_sets": len(architecture),
        "n_credible_set_members": len(members),
        "n_colocalized_target_rows": len(targets),
        "output_dir": str(output_dir.resolve()),
    }

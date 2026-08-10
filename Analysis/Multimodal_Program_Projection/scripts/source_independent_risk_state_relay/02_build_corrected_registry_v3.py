#!/usr/bin/env python3
"""Build Plan 45 registry v3 from the complete corrected-coloc portfolio.

This is deliberately fail-closed.  It writes nothing until all 50 studies by
22 autosomes have one exact, nonempty, schema-valid output.  It rebuilds the
genetic/state classes instead of patching Plan 43's July-derived registry.
The resulting locus registry is suitable for candidate routing, but it does
not orient an experimental perturbation; RSR-02 must re-export the selected
SuSiE credible-set pair and variant-level shared posterior separately.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

from relay_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_json,
    read_tsv,
    sha256_file,
    write_tsv,
)


RERUN_ROOT = PROJECT_ROOT / "GWAS/finemapping/results/susie_coloc_rerun"
GWAS_REGISTRY = PROJECT_ROOT / "GWAS/finemapping/config/gwas_registry.tsv"
TRAIT_TIER = PROJECT_ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv"
PHENOTYPE_REGISTRY = (
    PROJECT_ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / "chronic-state-risk-bridge-2026-08-09/frozen_inputs"
    / "phenotype_registry__phenotype_registry.tsv"
)
CANONICAL_BULK = (
    PROJECT_ROOT
    / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration"
    / "canonical_deg_results.csv"
)
PLAN43_CLASSES = (
    PROJECT_ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / "chronic-state-risk-bridge-2026-08-09/integrity_correction_v2"
    / "corrected_evidence_classes_v2.tsv"
)

EXPECTED_STUDIES = 50
EXPECTED_CHROMOSOMES = 22
EXPECTED_PRODUCTS = EXPECTED_STUDIES * EXPECTED_CHROMOSOMES
LOCUS_WINDOW_BP = 1_000_000
REQUIRED_COLOC_COLUMNS = {
    "gene",
    "ensembl",
    "chr",
    "gwas_name",
    "PP.H4.abf",
    "PP.H4.susie",
    "method",
    "top_snp",
    "top_snp_PP",
}


class RegistryGateError(RuntimeError):
    """Raised when a source or scientific registry gate is not satisfied."""


def clean(value: object) -> str:
    text = "" if value is None else str(value).strip()
    return "" if text in {"NA", "NaN", "nan", "None", "NULL", "."} else text


def finite_float(value: object) -> float | None:
    try:
        result = float(clean(value))
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def ensembl_base(value: object) -> str:
    return clean(value).split(".", 1)[0]


def parse_top_snp(value: object) -> tuple[str, int] | None:
    parts = clean(value).replace("_", ":").split(":")
    if len(parts) < 2:
        return None
    try:
        position = int(parts[1])
    except ValueError:
        return None
    return parts[0].removeprefix("chr"), position


def rank_value(value: object) -> float:
    parsed = finite_float(value)
    return -1.0 if parsed is None else parsed


def gwas_evidence_family(
    study: str, metadata: dict[str, str]
) -> tuple[str, str, bool]:
    """Collapse known participant-overlap families without inventing independence."""

    note = metadata.get("source_note", "").lower()
    stratum = metadata.get("phenotype_stratum", "")
    if (
        stratum == "mri_pdff_or_histologic_steatosis"
        and "uk biobank imaging cohort" in note
    ):
        return (
            "UKBB_PDFF_imaging_overlap_family",
            "known_overlapping_ukbb_imaging",
            True,
        )
    if "meta-analysis includes participant pools reused" in note:
        return (
            f"overlap_unresolved_meta::{study}",
            "participant_overlap_unresolved_meta_analysis",
            False,
        )
    return (
        f"portfolio_study::{study}",
        "distinct_portfolio_study_not_proven_independent",
        True,
    )


def expected_product_paths(studies: Iterable[str]) -> set[Path]:
    return {
        RERUN_ROOT / study / f"susie_coloc_chr{chrom}.csv"
        for study in studies
        for chrom in range(1, EXPECTED_CHROMOSOMES + 1)
    }


def validate_complete_portfolio(
    studies: set[str],
) -> tuple[list[Path], list[dict[str, object]]]:
    expected = expected_product_paths(studies)
    observed = set(RERUN_ROOT.glob("*/susie_coloc_chr*.csv"))
    missing = expected - observed
    unexpected = observed - expected
    empty = {path for path in observed if path.stat().st_size == 0}
    if missing or unexpected or empty or len(observed) != EXPECTED_PRODUCTS:
        raise RegistryGateError(
            "Corrected-coloc portfolio incomplete: "
            f"valid={len((observed & expected) - empty)}/{EXPECTED_PRODUCTS}, "
            f"missing={len(missing)}, unexpected={len(unexpected)}, empty={len(empty)}"
        )

    manifest: list[dict[str, object]] = []
    paths = sorted(expected)
    for path in paths:
        study = path.parent.name
        chrom_text = path.stem.removeprefix("susie_coloc_chr")
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            header = set(reader.fieldnames or [])
            absent = REQUIRED_COLOC_COLUMNS - header
            if absent:
                raise RegistryGateError(f"{path} lacks columns: {sorted(absent)}")
            n_rows = 0
            for row in reader:
                n_rows += 1
                if clean(row.get("gwas_name")) != study:
                    raise RegistryGateError(f"Study mismatch in {path} row {n_rows}")
                if clean(row.get("chr")) != chrom_text:
                    raise RegistryGateError(f"Chromosome mismatch in {path} row {n_rows}")
            if n_rows == 0:
                raise RegistryGateError(f"Header-only corrected-coloc product: {path}")
        manifest.append(
            {
                "study_name": study,
                "chromosome": chrom_text,
                "source_path": str(path.relative_to(PROJECT_ROOT)),
                "size_bytes": path.stat().st_size,
                "n_rows": n_rows,
                "sha256": sha256_file(path),
            }
        )
    return paths, manifest


def collapse_bulk() -> dict[str, dict[str, object]]:
    by_symbol: dict[str, list[dict[str, str]]] = defaultdict(list)
    with CANONICAL_BULK.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            symbol = clean(row.get("symbol"))
            if symbol:
                by_symbol[symbol].append(row)

    collapsed: dict[str, dict[str, object]] = {}
    for symbol, candidates in by_symbol.items():
        selected = sorted(
            candidates,
            key=lambda row: (
                finite_float(row.get("treat_fdr"))
                if finite_float(row.get("treat_fdr")) is not None
                else float("inf"),
                -abs(finite_float(row.get("t")) or 0.0),
                clean(row.get("gene")),
            ),
        )[0]
        collapsed[symbol] = {
            "ensembl_bulk": ";".join(
                sorted(
                    {
                        ensembl_base(row.get("gene"))
                        for row in candidates
                        if ensembl_base(row.get("gene"))
                    }
                )
            ),
            "bulk_row_count": len(candidates),
            "bulk_logFC": clean(selected.get("logFC")),
            "bulk_t": clean(selected.get("t")),
            "bulk_AveExpr": clean(selected.get("AveExpr")),
            "bulk_treat_fdr": clean(selected.get("treat_fdr")),
            "established_state_associated": any(
                finite_float(row.get("treat_fdr")) is not None
                and float(row["treat_fdr"]) < 0.05
                for row in candidates
            ),
        }
    return collapsed


def load_corrected_genetics(
    paths: list[Path], tier12_studies: set[str], phenotype: dict[str, dict[str, str]]
) -> dict[str, dict[str, object]]:
    by_symbol: dict[str, list[dict[str, str]]] = defaultdict(list)
    for path in paths:
        study = path.parent.name
        if study not in tier12_studies:
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                symbol = clean(row.get("gene"))
                ensembl = ensembl_base(row.get("ensembl"))
                if not symbol or not ensembl:
                    continue
                enriched = dict(row)
                enriched["ensembl_base"] = ensembl
                enriched["study_name"] = study
                by_symbol[symbol].append(enriched)

    genetics: dict[str, dict[str, object]] = {}
    for symbol, rows in by_symbol.items():
        selected = sorted(
            rows,
            key=lambda row: (
                -rank_value(row.get("PP.H4.susie")),
                -rank_value(row.get("PP.H4.abf")),
                row["study_name"],
                row["ensembl_base"],
            ),
        )[0]
        max_susie = max((rank_value(row.get("PP.H4.susie")) for row in rows), default=-1.0)
        max_abf = max((rank_value(row.get("PP.H4.abf")) for row in rows), default=-1.0)
        study = selected["study_name"]
        metadata = phenotype[study]
        genetics[symbol] = {
            "ensembl_genetic": ";".join(sorted({row["ensembl_base"] for row in rows})),
            "driving_ensembl": selected["ensembl_base"],
            "genetic_row_count": len(rows),
            "coloc_best_susie_pp4": "" if max_susie < 0 else max_susie,
            "coloc_best_abf_pp4": "" if max_abf < 0 else max_abf,
            "driving_gwas": study,
            "driving_trait": metadata["trait"],
            "driving_tier": metadata["tier"],
            "driving_phenotype_stratum": metadata["phenotype_stratum"],
            "driving_ancestry": metadata["ancestry"],
            "driving_regulatory_ancestry_status": metadata[
                "regulatory_ancestry_status"
            ],
            "driving_top_snp": clean(selected.get("top_snp")),
            "driving_top_snp_pp": clean(selected.get("top_snp_PP")),
            "driving_method": clean(selected.get("method")),
            "primary_genetic": max_susie > 0.5,
        }
    return genetics


def load_corrected_gene_gwas_pairs(
    paths: list[Path], tier12_studies: set[str], phenotype: dict[str, dict[str, str]]
) -> list[dict[str, object]]:
    """Retain every gene-by-GWAS record needed for cross-study orientation.

    The gene-level evidence registry intentionally uses the strongest Tier-1/2
    result to classify a gene.  That collapse must not be reused for target
    orientation: it would discard every corroborating or contradictory GWAS.
    """

    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for path in paths:
        study = path.parent.name
        if study not in tier12_studies:
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                symbol = clean(row.get("gene"))
                ensembl = ensembl_base(row.get("ensembl"))
                if symbol and ensembl:
                    grouped[(symbol, ensembl, study)].append(dict(row))

    pairs: list[dict[str, object]] = []
    for (symbol, ensembl, study), rows in sorted(grouped.items()):
        selected = sorted(
            rows,
            key=lambda row: (
                -rank_value(row.get("PP.H4.susie")),
                -rank_value(row.get("PP.H4.abf")),
                clean(row.get("top_snp")),
            ),
        )[0]
        metadata = phenotype[study]
        evidence_family, independence_class, breadth_eligible = gwas_evidence_family(
            study, metadata
        )
        parsed = parse_top_snp(selected.get("top_snp"))
        chrom, position = parsed if parsed is not None else ("", "")
        susie = finite_float(selected.get("PP.H4.susie"))
        pairs.append(
            {
                "gene_symbol": symbol,
                "ensembl_id": ensembl,
                "gwas_name": study,
                "trait": metadata["trait"],
                "tier": metadata["tier"],
                "phenotype_stratum": metadata["phenotype_stratum"],
                "ancestry": metadata["ancestry"],
                "regulatory_ancestry_status": metadata[
                    "regulatory_ancestry_status"
                ],
                "source_dependence": metadata.get("source_dependence", ""),
                "gwas_evidence_family": evidence_family,
                "gwas_independence_class": independence_class,
                "counts_for_replication_breadth": bool_text(breadth_eligible),
                "susie_pp4": clean(selected.get("PP.H4.susie")),
                "abf_pp4": clean(selected.get("PP.H4.abf")),
                "method": clean(selected.get("method")),
                "top_snp": clean(selected.get("top_snp")),
                "top_snp_pp": clean(selected.get("top_snp_PP")),
                "chromosome": chrom,
                "position": position,
                "pair_primary_genetic": bool_text(susie is not None and susie > 0.5),
                "pair_row_count": len(rows),
                "coarse_locus_uid": "",
                "static_class": "",
                "is_representative_gene_in_locus": "false",
                "representative_gene": "",
                "orientation_status": (
                    "requires_targeted_credible_set_pair_export"
                    if susie is not None and susie > 0.5
                    else "not_applicable"
                ),
            }
        )
    return pairs


def build_classes(
    bulk: dict[str, dict[str, object]], genetics: dict[str, dict[str, object]]
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for symbol in sorted(set(bulk) | set(genetics)):
        b = bulk.get(symbol, {})
        g = genetics.get(symbol, {})
        bulk_tested = symbol in bulk
        genetic_tested = symbol in genetics
        joint = bulk_tested and genetic_tested
        primary = bool(g.get("primary_genetic", False))
        state = bool(b.get("established_state_associated", False))
        if not joint:
            evidence_class = "indeterminate_not_jointly_testable"
        elif primary and state:
            evidence_class = "convergent"
        elif primary:
            evidence_class = "genetic_only"
        elif state:
            evidence_class = "disease_state_only"
        else:
            evidence_class = "neither"
        rows.append(
            {
                "gene_symbol": symbol,
                "ensembl_bulk": b.get("ensembl_bulk", ""),
                "ensembl_genetic": g.get("ensembl_genetic", ""),
                "driving_ensembl": g.get("driving_ensembl", ""),
                "bulk_tested": bool_text(bulk_tested),
                "primary_genetic_map_tested": bool_text(genetic_tested),
                "joint_testable": bool_text(joint),
                "primary_genetic": bool_text(primary),
                "established_state_associated": bool_text(state),
                "static_class": evidence_class,
                "bulk_row_count": b.get("bulk_row_count", 0),
                "bulk_logFC": b.get("bulk_logFC", ""),
                "bulk_t": b.get("bulk_t", ""),
                "bulk_AveExpr": b.get("bulk_AveExpr", ""),
                "bulk_treat_fdr": b.get("bulk_treat_fdr", ""),
                "genetic_row_count": g.get("genetic_row_count", 0),
                "coloc_best_susie_pp4": g.get("coloc_best_susie_pp4", ""),
                "coloc_best_abf_pp4": g.get("coloc_best_abf_pp4", ""),
                "driving_gwas": g.get("driving_gwas", ""),
                "driving_trait": g.get("driving_trait", ""),
                "driving_tier": g.get("driving_tier", ""),
                "driving_phenotype_stratum": g.get("driving_phenotype_stratum", ""),
                "driving_ancestry": g.get("driving_ancestry", ""),
                "driving_regulatory_ancestry_status": g.get(
                    "driving_regulatory_ancestry_status", ""
                ),
                "driving_top_snp": g.get("driving_top_snp", ""),
                "driving_top_snp_pp": g.get("driving_top_snp_pp", ""),
                "driving_method": g.get("driving_method", ""),
                "orientation_status": (
                    "requires_targeted_credible_set_pair_export" if primary else "not_applicable"
                ),
            }
        )
    return rows


def build_locus_registry(classes: list[dict[str, object]]) -> list[dict[str, object]]:
    provisional: list[dict[str, object]] = []
    unresolved_counter = 0
    for row in classes:
        if row["primary_genetic"] != "true":
            continue
        parsed = parse_top_snp(row["driving_top_snp"])
        if parsed is None:
            unresolved_counter += 1
            chrom, position = "", ""
            coarse_uid = f"unresolved:{row['ensembl_genetic']}:{unresolved_counter:04d}"
        else:
            chrom, position = parsed
            coarse_uid = ""
        provisional.append(
            {
                "gene_symbol": row["gene_symbol"],
                "ensembl_id": row["driving_ensembl"],
                "static_class": row["static_class"],
                "driving_gwas": row["driving_gwas"],
                "driving_trait": row["driving_trait"],
                "phenotype_stratum": row["driving_phenotype_stratum"],
                "ancestry": row["driving_ancestry"],
                "regulatory_ancestry_status": row[
                    "driving_regulatory_ancestry_status"
                ],
                "coarse_top_snp": row["driving_top_snp"],
                "chromosome": chrom,
                "position": position,
                "susie_pp4": row["coloc_best_susie_pp4"],
                "abf_pp4": row["coloc_best_abf_pp4"],
                "method": row["driving_method"],
                "coarse_locus_uid": coarse_uid,
                "locus_resolution": (
                    "cross_study_top_snp_1mb_component"
                    if parsed is not None
                    else "unresolved_singleton"
                ),
                "orientation_status": "requires_targeted_credible_set_pair_export",
            }
        )

    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in provisional:
        if row["position"] != "":
            grouped[str(row["chromosome"])].append(row)
    for chrom, rows in grouped.items():
        rows.sort(key=lambda item: (int(item["position"]), str(item["ensembl_id"])))
        component = 1
        previous: int | None = None
        for row in rows:
            position = int(row["position"])
            if previous is not None and position - previous > LOCUS_WINDOW_BP:
                component += 1
            row["coarse_locus_uid"] = f"chr{chrom}:component{component:04d}"
            previous = position

    by_locus: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in provisional:
        by_locus[str(row["coarse_locus_uid"])].append(row)
    for locus_rows in by_locus.values():
        selected = sorted(
            locus_rows,
            key=lambda item: (
                -rank_value(item["susie_pp4"]),
                -rank_value(item["abf_pp4"]),
                str(item["ensembl_id"]),
            ),
        )[0]
        for item in locus_rows:
            item["n_primary_genetic_genes_in_locus"] = len(locus_rows)
            item["is_representative"] = bool_text(item is selected)
            item["representative_gene"] = selected["gene_symbol"]
    return sorted(
        provisional,
        key=lambda row: (str(row["coarse_locus_uid"]), str(row["ensembl_id"])),
    )


def build_pair_level_locus_registries(
    classes: list[dict[str, object]], pairs: list[dict[str, object]]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Assign all primary gene-GWAS pairs to physical signal components.

    Unlike the legacy helper above, this component construction happens before
    collapsing to one driving GWAS. It therefore preserves same-locus studies
    for the required cross-GWAS orientation check.
    """

    class_by_gene = {str(row["gene_symbol"]): row for row in classes}
    primary_indices: list[int] = []
    by_chromosome: dict[str, list[int]] = defaultdict(list)
    unresolved_counter = 0
    for index, row in enumerate(pairs):
        cls = class_by_gene.get(str(row["gene_symbol"]), {})
        row["static_class"] = cls.get("static_class", "indeterminate_not_jointly_testable")
        if row["pair_primary_genetic"] != "true":
            continue
        primary_indices.append(index)
        if row["chromosome"] == "" or row["position"] == "":
            unresolved_counter += 1
            row["coarse_locus_uid"] = (
                f"unresolved:{row['ensembl_id']}:{row['gwas_name']}:{unresolved_counter:04d}"
            )
        else:
            by_chromosome[str(row["chromosome"])].append(index)

    for chromosome, indices in by_chromosome.items():
        indices.sort(
            key=lambda idx: (
                int(pairs[idx]["position"]),
                str(pairs[idx]["ensembl_id"]),
                str(pairs[idx]["gwas_name"]),
            )
        )
        component = 1
        previous: int | None = None
        for index in indices:
            position = int(pairs[index]["position"])
            if previous is not None and position - previous > LOCUS_WINDOW_BP:
                component += 1
            pairs[index]["coarse_locus_uid"] = f"chr{chromosome}:component{component:04d}"
            previous = position

    by_locus_gene: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for index in primary_indices:
        row = pairs[index]
        by_locus_gene[(str(row["coarse_locus_uid"]), str(row["ensembl_id"]))].append(row)

    locus_rows: list[dict[str, object]] = []
    by_locus: dict[str, list[dict[str, object]]] = defaultdict(list)
    for (locus_uid, ensembl), gene_pairs in sorted(by_locus_gene.items()):
        selected = sorted(
            gene_pairs,
            key=lambda row: (
                -rank_value(row["susie_pp4"]),
                -rank_value(row["abf_pp4"]),
                str(row["gwas_name"]),
            ),
        )[0]
        directions_pending = len({str(row["gwas_name"]) for row in gene_pairs})
        summary = {
            "gene_symbol": selected["gene_symbol"],
            "ensembl_id": ensembl,
            "static_class": selected["static_class"],
            "driving_gwas": selected["gwas_name"],
            "driving_trait": selected["trait"],
            "phenotype_stratum": selected["phenotype_stratum"],
            "ancestry": selected["ancestry"],
            "regulatory_ancestry_status": selected["regulatory_ancestry_status"],
            "coarse_top_snp": selected["top_snp"],
            "chromosome": selected["chromosome"],
            "position": selected["position"],
            "susie_pp4": selected["susie_pp4"],
            "abf_pp4": selected["abf_pp4"],
            "method": selected["method"],
            "coarse_locus_uid": locus_uid,
            "locus_resolution": (
                "cross_study_pair_top_snp_1mb_component"
                if not locus_uid.startswith("unresolved:")
                else "unresolved_pair_singleton"
            ),
            "n_primary_gene_gwas_pairs": len(gene_pairs),
            "n_primary_gwas_for_gene_in_locus": directions_pending,
            "orientation_status": "requires_all_eligible_pair_exports",
        }
        locus_rows.append(summary)
        by_locus[locus_uid].append(summary)

    for locus_uid, rows in by_locus.items():
        selected = sorted(
            rows,
            key=lambda row: (
                -rank_value(row["susie_pp4"]),
                -rank_value(row["abf_pp4"]),
                str(row["ensembl_id"]),
            ),
        )[0]
        for row in rows:
            row["n_primary_genetic_genes_in_locus"] = len(rows)
            row["is_representative"] = bool_text(row is selected)
            row["representative_gene"] = selected["gene_symbol"]
        for pair in pairs:
            if pair["coarse_locus_uid"] == locus_uid:
                pair["representative_gene"] = selected["gene_symbol"]
                pair["is_representative_gene_in_locus"] = bool_text(
                    pair["ensembl_id"] == selected["ensembl_id"]
                )

    return (
        sorted(
            locus_rows,
            key=lambda row: (str(row["coarse_locus_uid"]), str(row["ensembl_id"])),
        ),
        sorted(
            pairs,
            key=lambda row: (
                str(row["coarse_locus_uid"]),
                str(row["ensembl_id"]),
                str(row["gwas_name"]),
            ),
        ),
    )


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RegistryGateError(f"Refusing to overwrite candidate root: {CANDIDATE_ROOT}")

    required = [
        GWAS_REGISTRY,
        TRAIT_TIER,
        PHENOTYPE_REGISTRY,
        CANONICAL_BULK,
        PLAN43_CLASSES,
    ]
    absent = [str(path) for path in required if not path.is_file()]
    if absent:
        raise RegistryGateError("Missing registry-v3 inputs: " + "; ".join(absent))

    gwas_rows = read_tsv(GWAS_REGISTRY)
    studies = {row["study_name"] for row in gwas_rows}
    if len(gwas_rows) != EXPECTED_STUDIES or len(studies) != EXPECTED_STUDIES:
        raise RegistryGateError(
            f"Expected {EXPECTED_STUDIES} unique GWAS studies; observed {len(studies)}"
        )
    tier_rows = read_tsv(TRAIT_TIER)
    tier12 = {
        row["study_name"]
        for row in tier_rows
        if row["tier"] in {"1", "2"} and row["placement"] == "main"
    }
    phenotype = {row["study_name"]: row for row in read_tsv(PHENOTYPE_REGISTRY)}
    if tier12 != set(phenotype):
        raise RegistryGateError(
            "Tier-1/2 phenotype registry mismatch: "
            f"missing={sorted(tier12 - set(phenotype))}, "
            f"unexpected={sorted(set(phenotype) - tier12)}"
        )

    # The complete portfolio is checked before candidate-root creation or any write.
    coloc_paths, coloc_manifest = validate_complete_portfolio(studies)
    bulk = collapse_bulk()
    genetics = load_corrected_genetics(coloc_paths, tier12, phenotype)
    classes = build_classes(bulk, genetics)
    pairs = load_corrected_gene_gwas_pairs(coloc_paths, tier12, phenotype)
    loci, pairs = build_pair_level_locus_registries(classes, pairs)
    if not loci:
        raise RegistryGateError("Corrected registry produced no primary genetic loci")

    CANDIDATE_ROOT.mkdir(parents=True)
    source_manifest = [
        {
            "source_id": path.stem,
            "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in required
    ]
    write_tsv(
        CANDIDATE_ROOT / "corrected_coloc_product_manifest.tsv",
        coloc_manifest,
        ["study_name", "chromosome", "source_path", "size_bytes", "n_rows", "sha256"],
    )
    write_tsv(
        CANDIDATE_ROOT / "registry_v3_source_manifest.tsv",
        source_manifest,
        ["source_id", "source_path", "size_bytes", "sha256"],
    )
    write_tsv(
        CANDIDATE_ROOT / "corrected_evidence_registry_v3.tsv",
        classes,
        list(classes[0]),
    )
    write_tsv(
        CANDIDATE_ROOT / "corrected_primary_genetic_locus_registry_v3.tsv",
        loci,
        list(loci[0]),
    )
    write_tsv(
        CANDIDATE_ROOT / "corrected_gene_gwas_pair_registry_v3.tsv",
        pairs,
        list(pairs[0]),
    )

    previous_counts = Counter(
        row["static_class"] for row in read_tsv(PLAN43_CLASSES)
    )
    current_counts = Counter(str(row["static_class"]) for row in classes)
    audit_rows: list[dict[str, object]] = []
    for evidence_class in sorted(set(previous_counts) | set(current_counts)):
        audit_rows.append(
            {
                "metric": f"class_count::{evidence_class}",
                "plan43_value": previous_counts.get(evidence_class, 0),
                "registry_v3_value": current_counts.get(evidence_class, 0),
                "delta": current_counts.get(evidence_class, 0)
                - previous_counts.get(evidence_class, 0),
                "interpretation": "reported_drift_not_acceptance_threshold",
            }
        )
    audit_rows.extend(
        [
            {
                "metric": "primary_genetic_gene_rows",
                "plan43_value": sum(
                    row["primary_genetic"] == "true" for row in read_tsv(PLAN43_CLASSES)
                ),
                "registry_v3_value": sum(
                    row["primary_genetic"] == "true" for row in classes
                ),
                "delta": "",
                "interpretation": "corrected portfolio result",
            },
            {
                "metric": "coarse_physical_loci",
                "plan43_value": 326,
                "registry_v3_value": len({row["coarse_locus_uid"] for row in loci}),
                "delta": "",
                "interpretation": "not expected to equal Plan43 study-specific loci",
            },
        ]
    )
    write_tsv(
        CANDIDATE_ROOT / "registry_v3_audit.tsv",
        audit_rows,
        ["metric", "plan43_value", "registry_v3_value", "delta", "interpretation"],
    )

    outputs = [
        CANDIDATE_ROOT / "corrected_coloc_product_manifest.tsv",
        CANDIDATE_ROOT / "registry_v3_source_manifest.tsv",
        CANDIDATE_ROOT / "corrected_evidence_registry_v3.tsv",
        CANDIDATE_ROOT / "corrected_primary_genetic_locus_registry_v3.tsv",
        CANDIDATE_ROOT / "corrected_gene_gwas_pair_registry_v3.tsv",
        CANDIDATE_ROOT / "registry_v3_audit.tsv",
    ]
    payload = {
        "status": "registry_v3_frozen_orientation_pending",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "candidate_id": CANDIDATE_ROOT.name,
        "git_head": subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=PROJECT_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip(),
        "n_corrected_coloc_products": len(coloc_paths),
        "n_tier12_studies": len(tier12),
        "n_evidence_rows": len(classes),
        "class_counts": dict(sorted(current_counts.items())),
        "n_primary_genetic_rows": len(loci),
        "n_coarse_physical_loci": len({row["coarse_locus_uid"] for row in loci}),
        "n_gene_gwas_pair_rows": len(pairs),
        "n_primary_gene_gwas_pairs": sum(
            row["pair_primary_genetic"] == "true" for row in pairs
        ),
        "orientation_complete": False,
        "experimental_target_list_frozen": False,
        "orientation_gate": "targeted selected-pair export with variant-level SNP.PP.H4 required",
        "output_sha256": {
            path.name: sha256_file(path) for path in outputs
        },
    }
    atomic_write_json(CANDIDATE_ROOT / "REGISTRY_V3_SEALED.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

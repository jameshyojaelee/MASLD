#!/usr/bin/env python3
"""Reproduce source-positive eGenes and assemble power/observability metadata."""

from __future__ import annotations

import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path

from genetics_common import (
    ContractError,
    DEFAULT_CANDIDATE_ROOT,
    PROJECT_ROOT,
    assert_candidate_root,
    atomic_write_tsv,
    bool_text,
    clean,
    ensembl_base,
    iter_table,
    parse_float,
    read_table,
    write_stage_seal,
)


SOURCE_EXPECTED = {
    "deposited_signals": 9013,
    "source_defined_egenes": 6564,
    "primary_signal_rows": 6564,
    "signals_joint_p_le_1e-5": 9013,
    "primary_egenes_fdr_le_0_05": 6556,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--candidate-root", type=Path, default=DEFAULT_CANDIDATE_ROOT)
    return parser.parse_args()


def manifest_inputs(candidate: Path) -> dict[str, Path]:
    rows = read_table(candidate / "input_manifest.tsv", "tsv")
    result = {}
    for row in rows:
        result[row["input_id"]] = (
            candidate / row["snapshot_path"]
            if row["snapshot_path"]
            else Path(row["source_path"])
        )
    return result


def parse_attributes(text: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]+)"', text))


def read_gene_spans(gtf_path: Path) -> dict[str, dict[str, object]]:
    spans: dict[str, dict[str, object]] = {}
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = parse_attributes(fields[8])
            ens = ensembl_base(attrs.get("gene_id"))
            if not ens:
                continue
            start, end = int(fields[3]), int(fields[4])
            span = end - start + 1
            current = spans.get(ens)
            if current is None or span > int(current["gene_span_bp"]):
                spans[ens] = {
                    "gene_span_bp": span,
                    "gtf_chromosome": fields[0],
                    "gtf_gene_name": clean(attrs.get("gene_name")),
                }
    return spans


def reproduce_source_leads(rows: list[dict[str, str]]):
    by_ensembl: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        ens = ensembl_base(row.get("Ensembl"))
        if not ens:
            raise ContractError("Broadaway lead row has no Ensembl identifier")
        by_ensembl[ens].append(row)

    primary_rows = [row for row in rows if clean(row.get("Signal")).startswith("1:")]
    observed = {
        "deposited_signals": len(rows),
        "source_defined_egenes": len(by_ensembl),
        "primary_signal_rows": len(primary_rows),
        "signals_joint_p_le_1e-5": sum(
            parse_float(row.get("Joint_Pvalue")) is not None
            and float(row["Joint_Pvalue"]) <= 1e-5
            for row in rows
        ),
        "primary_egenes_fdr_le_0_05": sum(
            parse_float(row.get("FDR")) is not None
            and float(row["FDR"]) <= 0.05
            for row in primary_rows
        ),
    }
    failures = [
        f"{metric}: expected {expected}, observed {observed[metric]}"
        for metric, expected in SOURCE_EXPECTED.items()
        if observed[metric] != expected
    ]
    if failures:
        raise ContractError("Broadaway source reproduction failed: " + "; ".join(failures))

    summary: dict[str, dict[str, object]] = {}
    symbol_to_ensembl: dict[str, set[str]] = defaultdict(set)
    for ens, signals in by_ensembl.items():
        primaries = [row for row in signals if clean(row.get("Signal")).startswith("1:")]
        if len(primaries) != 1:
            raise ContractError(f"expected one primary Broadaway signal for {ens}; found {len(primaries)}")
        primary = primaries[0]
        symbols = {clean(row.get("Gene")) for row in signals if clean(row.get("Gene"))}
        for symbol in symbols:
            symbol_to_ensembl[symbol].add(ens)
        joint_values = [parse_float(row.get("Joint_Pvalue")) for row in signals]
        stepwise_values = [parse_float(row.get("Stepwise_Pvalue")) for row in signals]
        n_values = [parse_float(row.get("N")) for row in signals]
        summary[ens] = {
            "source_symbol": clean(primary.get("Gene")),
            "source_biotype": clean(primary.get("Gene_BioType")),
            "source_signal_count": len(signals),
            "source_primary_fdr": clean(primary.get("FDR")),
            "source_primary_eigenmt_p": clean(primary.get("eigenMT_adjusted_Pvalue")),
            "source_min_joint_p": min(value for value in joint_values if value is not None),
            "source_min_stepwise_p": min(value for value in stepwise_values if value is not None),
            "source_eqtl_n_min": int(min(value for value in n_values if value is not None)),
            "source_eqtl_n_max": int(max(value for value in n_values if value is not None)),
            "source_defined_egene": True,
            "source_primary_fdr_05": (
                parse_float(primary.get("FDR")) is not None
                and float(primary["FDR"]) <= 0.05
            ),
        }
    return summary, symbol_to_ensembl, observed


def coverage_from_coloc(
    all_gwas_path: Path, target_ensembl: set[str]
) -> dict[str, dict[str, object]]:
    coverage: dict[str, dict[str, object]] = defaultdict(
        lambda: {
            "coloc_rows": 0,
            "gwas": set(),
            "n_snps": [],
            "susie_rows": 0,
            "ld_status": set(),
        }
    )
    for row in iter_table(all_gwas_path, "csv"):
        ens = ensembl_base(row.get("ensembl"))
        if ens not in target_ensembl:
            continue
        record = coverage[ens]
        record["coloc_rows"] += 1
        if clean(row.get("gwas_name")):
            record["gwas"].add(row["gwas_name"])
        n_snps = parse_float(row.get("n_snps"))
        if n_snps is not None:
            record["n_snps"].append(int(n_snps))
        if parse_float(row.get("PP.H4.susie")) is not None:
            record["susie_rows"] += 1
        if clean(row.get("ld_reliability")):
            record["ld_status"].add(row["ld_reliability"])
    return coverage


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    candidate = assert_candidate_root(project, args.candidate_root)
    stage2 = candidate / "work" / "stage_seals" / "02_build_phenotype_registry.json"
    if not stage2.is_file():
        raise ContractError("phenotype-registry seal is missing")
    inputs = manifest_inputs(candidate)
    frozen = read_table(candidate / "frozen_evidence_classes.tsv", "tsv")
    leads = read_table(inputs["broadaway_leads"], "tsv")
    source, source_symbol_map, source_observed = reproduce_source_leads(leads)

    gencode_rows = read_table(inputs["gencode_metadata"], "tsv_gz")
    gencode = {
        ensembl_base(row["ensembl_base"]): {
            "gencode_symbol": clean(row["gene_name"]),
            "gencode_biotype": clean(row["gene_biotype"]),
            "gencode_chromosome": clean(row["chromosome"]),
        }
        for row in gencode_rows
        if ensembl_base(row["ensembl_base"])
    }
    spans = read_gene_spans(inputs["gencode_gtf"])

    target_ensembl: set[str] = set()
    for row in frozen:
        for column in ("ensembl_bulk", "ensembl_genetic"):
            target_ensembl.update(x for x in row[column].split(";") if x)
    coverage = coverage_from_coloc(inputs["all_gwas_coloc"], target_ensembl)

    observability_rows: list[dict[str, object]] = []
    mapping_counts: defaultdict[str, int] = defaultdict(int)
    for row in frozen:
        symbol = row["gene_symbol"]
        gene_ensembl = sorted(
            {
                ens
                for column in ("ensembl_bulk", "ensembl_genetic")
                for ens in row[column].split(";")
                if ens
            }
        )
        source_hits = [ens for ens in gene_ensembl if ens in source]
        mapping_basis = "ensembl"
        if not source_hits:
            symbol_hits = sorted(source_symbol_map.get(symbol, set()))
            if len(symbol_hits) == 1:
                source_hits = symbol_hits
                mapping_basis = "unique_symbol_fallback"
            elif len(symbol_hits) > 1:
                mapping_basis = "ambiguous_symbol_no_call"
            else:
                mapping_basis = "not_in_significant_lead_deposit"
        mapping_counts[mapping_basis] += 1

        source_records = [source[ens] for ens in source_hits]
        annotation_ens = gene_ensembl[0] if gene_ensembl else (source_hits[0] if source_hits else "")
        annotation = gencode.get(annotation_ens, {})
        span = spans.get(annotation_ens, {})
        cov_records = [coverage[ens] for ens in gene_ensembl if ens in coverage]
        n_snps = [value for record in cov_records for value in record["n_snps"]]
        gwas_names = sorted({value for record in cov_records for value in record["gwas"]})
        ld_status = sorted({value for record in cov_records for value in record["ld_status"]})
        observability_rows.append(
            {
                "gene_symbol": symbol,
                "ensembl_ids": ";".join(gene_ensembl),
                "static_class": row["static_class"],
                "joint_testable": row["joint_testable"],
                "primary_genetic": row["primary_genetic"],
                "established_state_associated": row["established_state_associated"],
                "bulk_AveExpr": row["bulk_AveExpr"],
                "bulk_t": row["bulk_t"],
                "bulk_treat_fdr": row["bulk_treat_fdr"],
                "coloc_best_susie_pp4": row["coloc_best_susie_pp4"],
                "coloc_best_abf_pp4": row["coloc_best_abf_pp4"],
                "source_eqtl_mapping_basis": mapping_basis,
                "source_eqtl_ensembl": ";".join(source_hits),
                "source_defined_egene": bool_text(bool(source_records)),
                "source_egene_absence_interpretation": (
                    "source_positive"
                    if source_records
                    else "indeterminate_complete_tested_universe_not_deposited"
                ),
                "source_signal_count": sum(int(x["source_signal_count"]) for x in source_records),
                "source_primary_fdr_min": min(
                    (parse_float(x["source_primary_fdr"]) for x in source_records),
                    default=None,
                    key=lambda value: float("inf") if value is None else value,
                )
                if source_records
                else "",
                "source_min_joint_p": min(
                    (float(x["source_min_joint_p"]) for x in source_records), default=""
                ),
                "source_min_stepwise_p": min(
                    (float(x["source_min_stepwise_p"]) for x in source_records), default=""
                ),
                "source_eqtl_n_min": min(
                    (int(x["source_eqtl_n_min"]) for x in source_records), default=""
                ),
                "source_eqtl_n_max": max(
                    (int(x["source_eqtl_n_max"]) for x in source_records), default=""
                ),
                "source_expression_value": "",
                "source_expression_status": "not_deposited_in_significant_lead_table",
                "gencode_biotype": annotation.get("gencode_biotype", ""),
                "gene_span_bp": span.get("gene_span_bp", ""),
                "gene_length_build": "GRCh38_GENCODE_v49",
                "local_tested_variant_density": "",
                "local_locus_density_status": "not_reconstructable_from_deposited_lead_table",
                "coloc_coverage_rows": sum(int(x["coloc_rows"]) for x in cov_records),
                "coloc_gwas_n": len(gwas_names),
                "coloc_gwas_names": ";".join(gwas_names),
                "coloc_n_snps_min_coverage_only": min(n_snps) if n_snps else "",
                "coloc_n_snps_max_coverage_only": max(n_snps) if n_snps else "",
                "coloc_susie_rows": sum(int(x["susie_rows"]) for x in cov_records),
                "coloc_ld_status": ";".join(ld_status),
                "adequate_static_eqtl_negative_authorized": "false",
                "adequate_detection_reason": (
                    "source_positive_egene_only"
                    if source_records
                    else "complete_source_tested_gene_universe_unavailable"
                ),
            }
        )

    observability_path = candidate / "gene_observability.tsv"
    observability_fields = list(observability_rows[0])
    atomic_write_tsv(observability_path, observability_rows, observability_fields)

    reproduction_rows = [
        {
            "metric": metric,
            "expected": expected,
            "observed": source_observed[metric],
            "status": "pass",
            "source_definition": (
                "all deposited conditionally distinct signals; joint p <= 1e-5"
                if metric in {"deposited_signals", "signals_joint_p_le_1e-5"}
                else "unique source eGenes, one primary Signal beginning 1: per Ensembl gene"
                if metric in {"source_defined_egenes", "primary_signal_rows"}
                else "primary-signal FDR <= 0.05; eight source eGenes remain source-defined but exceed this secondary FDR rule"
            ),
        }
        for metric, expected in SOURCE_EXPECTED.items()
    ]
    reproduction_path = candidate / "source_eqtl_reproduction.tsv"
    atomic_write_tsv(
        reproduction_path,
        reproduction_rows,
        ["metric", "expected", "observed", "status", "source_definition"],
    )

    universe_path = candidate / "source_tested_gene_universe_status.tsv"
    atomic_write_tsv(
        universe_path,
        [
            {
                "source": "Broadaway_liver_eQTL_meta",
                "reported_expression_genes": 18322,
                "deposited_source_significant_egenes": 6564,
                "complete_tested_gene_universe_available": "false",
                "complete_non_egene_background_available": "false",
                "source_expression_values_available": "false",
                "local_tested_variant_density_available": "false",
                "status": "coverage_limited",
                "authorized_use": "source_positive_eGene_annotation_only",
                "prohibited_use": "non_eGene_calls; adequate-negative claims; source-tested enrichment; power matching",
                "reason": "The lead deposit contains significant signals only. The paper reports 18,322 expression genes, but the exact complete tested-gene rows and source expression values are not present in this deposit; the local 19,072-gene OTTERS annotation and 18,975 downstream COLOC rows are derived analysis universes, not substitutes.",
            }
        ],
        [
            "source",
            "reported_expression_genes",
            "deposited_source_significant_egenes",
            "complete_tested_gene_universe_available",
            "complete_non_egene_background_available",
            "source_expression_values_available",
            "local_tested_variant_density_available",
            "status",
            "authorized_use",
            "prohibited_use",
            "reason",
        ],
    )

    readiness_rows = [
        {
            "requirement": "source_significant_egene_status",
            "available_genes": sum(row["source_defined_egene"] == "true" for row in observability_rows),
            "total_genes": len(observability_rows),
            "status": "partial_positive_only",
            "impact": "positive-source annotations are allowed; absence is indeterminate",
        },
        {
            "requirement": "complete_source_tested_gene_background",
            "available_genes": 0,
            "total_genes": len(observability_rows),
            "status": "missing",
            "impact": "nested source-tested universe and powered negatives are skipped",
        },
        {
            "requirement": "source_liver_expression",
            "available_genes": 0,
            "total_genes": len(observability_rows),
            "status": "missing",
            "impact": "source-expression matching is not authorized",
        },
        {
            "requirement": "local_tested_variant_density",
            "available_genes": 0,
            "total_genes": len(observability_rows),
            "status": "missing",
            "impact": "balanced power matching is not authorized",
        },
        {
            "requirement": "gene_span_and_biotype",
            "available_genes": sum(
                bool(row["gene_span_bp"]) and bool(row["gencode_biotype"])
                for row in observability_rows
            ),
            "total_genes": len(observability_rows),
            "status": "available_with_mapping_missingness",
            "impact": "descriptive annotation only until all matching covariates are available",
        },
        {
            "requirement": "downstream_n_snps",
            "available_genes": sum(bool(row["coloc_n_snps_max_coverage_only"]) for row in observability_rows),
            "total_genes": len(observability_rows),
            "status": "coverage_only",
            "impact": "must never be interpreted as eGene strength",
        },
    ]
    readiness_path = candidate / "observability_readiness.tsv"
    atomic_write_tsv(
        readiness_path,
        readiness_rows,
        ["requirement", "available_genes", "total_genes", "status", "impact"],
    )

    mapping_path = candidate / "gene_mapping_audit.tsv"
    atomic_write_tsv(
        mapping_path,
        [
            {
                "mapping_basis": basis,
                "gene_count": count,
                "interpretation": (
                    "source-positive eGene mapping"
                    if basis in {"ensembl", "unique_symbol_fallback"}
                    else "no source-negative inference"
                ),
            }
            for basis, count in sorted(mapping_counts.items())
        ],
        ["mapping_basis", "gene_count", "interpretation"],
    )
    write_stage_seal(
        candidate,
        "03_build_eqtl_observability",
        [
            observability_path,
            reproduction_path,
            universe_path,
            readiness_path,
            mapping_path,
        ],
        [stage2],
    )
    print(
        "PASS: reproduced 9,013 source signals / 6,564 eGenes; "
        "complete source-tested background remains coverage_limited"
    )


if __name__ == "__main__":
    main()

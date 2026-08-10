#!/usr/bin/env python3
"""Independently validate and seal the genetics source/power preflight."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from genetics_common import (
    ContractError,
    DEFAULT_CANDIDATE_ROOT,
    PROJECT_ROOT,
    SCRIPT_DIR,
    assert_candidate_root,
    atomic_write_tsv,
    clean,
    read_table,
    sha256_file,
    write_stage_seal,
)


EXPECTED_CONTRACT = {
    "treat_rows": 1918,
    "treat_symbols": 1915,
    "primary_susie_named": 473,
    "primary_susie_bulk_testable": 447,
    "primary_overlap": 34,
}
EXPECTED_SOURCE = {
    "deposited_signals": 9013,
    "source_defined_egenes": 6564,
    "primary_signal_rows": 6564,
    "signals_joint_p_le_1e-5": 9013,
    "primary_egenes_fdr_le_0_05": 6556,
}


def f(value: str) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result


def base(value: str) -> str:
    return clean(value).split(".", 1)[0]


def inputs(candidate: Path) -> dict[str, Path]:
    result = {}
    for row in read_table(candidate / "input_manifest.tsv", "tsv"):
        path = candidate / row["snapshot_path"] if row["snapshot_path"] else Path(row["source_path"])
        if sha256_file(path) != row["sha256"]:
            raise ContractError(f"input drift after freeze: {row['input_id']}")
        result[row["input_id"]] = path
    return result


def independent_contract(source_paths: dict[str, Path]) -> dict[str, int]:
    with source_paths["canonical_treat"].open(newline="") as handle:
        bulk = list(csv.DictReader(handle))
    treat_rows = [row for row in bulk if f(row["treat_fdr"]) is not None and f(row["treat_fdr"]) < 0.05]
    treat_symbols = {clean(row["symbol"]) for row in treat_rows if clean(row["symbol"])}
    bulk_symbols = {clean(row["symbol"]) for row in bulk if clean(row["symbol"])}

    with source_paths["full_genelevel"].open(newline="") as handle:
        full = list(csv.DictReader(handle))
    ens_to_gene: dict[str, set[str]] = defaultdict(set)
    for row in full:
        if base(row["ensembl"]) and clean(row["gene"]):
            ens_to_gene[base(row["ensembl"])].add(clean(row["gene"]))

    with source_paths["tier12_genelevel"].open(newline="") as handle:
        tier = csv.DictReader(handle)
        primary = set()
        for row in tier:
            pp4 = f(row["coloc_best_susie_pp4"])
            ens = base(row["ensembl"])
            if pp4 is None or pp4 <= 0.5 or not ens:
                continue
            symbols = ens_to_gene.get(ens, set())
            if len(symbols) != 1:
                raise ContractError(f"validator cannot uniquely name primary Ensembl gene {ens}")
            primary.update(symbols)

    primary_joint = primary & bulk_symbols
    return {
        "treat_rows": len(treat_rows),
        "treat_symbols": len(treat_symbols),
        "primary_susie_named": len(primary),
        "primary_susie_bulk_testable": len(primary_joint),
        "primary_overlap": len(primary_joint & treat_symbols),
    }


def independent_source(source_path: Path) -> dict[str, int]:
    with source_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    primary = [row for row in rows if row["Signal"].startswith("1:")]
    return {
        "deposited_signals": len(rows),
        "source_defined_egenes": len({base(row["Ensembl"]) for row in rows}),
        "primary_signal_rows": len(primary),
        "signals_joint_p_le_1e-5": sum(f(row["Joint_Pvalue"]) is not None and f(row["Joint_Pvalue"]) <= 1e-5 for row in rows),
        "primary_egenes_fdr_le_0_05": sum(f(row["FDR"]) is not None and f(row["FDR"]) <= 0.05 for row in primary),
    }


def require_equal(observed: dict[str, int], expected: dict[str, int], label: str) -> None:
    failures = [
        f"{metric}: expected {value}, observed {observed.get(metric)}"
        for metric, value in expected.items()
        if observed.get(metric) != value
    ]
    if failures:
        raise ContractError(f"{label} validation failed: " + "; ".join(failures))


def check_phenotypes(candidate: Path) -> list[dict[str, object]]:
    rows = read_table(candidate / "phenotype_registry.tsv", "tsv")
    if len(rows) != 35 or len({row["study_name"] for row in rows}) != 35:
        raise ContractError("phenotype registry is not one row per 35 primary GWAS")
    counts = Counter(row["phenotype_stratum"] for row in rows)
    expected = Counter(
        {
            "direct_masld_mash_diagnosis": 12,
            "mri_pdff_or_histologic_steatosis": 3,
            "alt_ast_or_ggt": 20,
        }
    )
    if counts != expected:
        raise ContractError(f"phenotype stratum drift: {dict(counts)}")
    non_eur = [row for row in rows if row["ancestry"] != "EUR"]
    if any(row["regulatory_ancestry_status"] != "cross_ancestry_eqtl_limited" for row in non_eur):
        raise ContractError("non-EUR GWAS missing European-eQTL limitation")
    if any(not row["source_citation"] or not row["source_url"] for row in rows):
        raise ContractError("phenotype row lacks source citation or URL")
    bbj = [row for row in rows if row["study_name"].startswith("BBJ_")]
    if len(bbj) != 3 or any(row["sample_size_audit"] != "source_registry_mismatch" for row in bbj):
        raise ContractError("BBJ sample-size conflicts were not preserved")
    intermountain = next(
        row for row in rows if row["study_name"] == "2023_36280732_NAFLD_Intermountain_EUR"
    )
    if intermountain["sample_size_audit"] != "rounded_source_registry_mismatch":
        raise ContractError("Intermountain rounded source/registry size conflict was not preserved")
    return [
        {
            "check_id": "phenotype_registry",
            "status": "pass",
            "observed": len(rows),
            "expected": 35,
            "note": "12 direct diagnosis, 3 PDFF, 20 enzymes; all non-EUR rows are cross_ancestry_eqtl_limited",
        }
    ]


def check_fail_closed(candidate: Path) -> list[dict[str, object]]:
    universe = read_table(candidate / "source_tested_gene_universe_status.tsv", "tsv")
    if len(universe) != 1 or universe[0]["status"] != "coverage_limited":
        raise ContractError("source tested-gene universe is not coverage_limited")
    if universe[0]["complete_tested_gene_universe_available"] != "false":
        raise ContractError("complete source tested-gene universe was incorrectly asserted")
    interface = read_table(candidate / "power_stratified_interface.tsv", "tsv")
    by_name = {row["universe"]: row for row in interface}
    try:
        all_joint_overlap = int(by_name["all_joint_testable"]["n_overlap"].strip())
        all_joint_denominator = int(
            by_name["all_joint_testable"]["overlap_denominator"].strip()
        )
    except (KeyError, ValueError) as error:
        raise ContractError("all-joint interface has malformed integer fields") from error
    if all_joint_overlap != 34 or all_joint_denominator != 447:
        raise ContractError("all-joint interface does not preserve 34/447")
    for name in (
        "bulk_expressed_and_source_eqtl_tested",
        "source_significant_egene_complete_covariates",
    ):
        if not by_name[name]["status"].startswith("skipped_"):
            raise ContractError(f"incomplete source universe was not skipped: {name}")
        if by_name[name]["source_negative_authorized"] != "false":
            raise ContractError(f"source negative incorrectly authorized: {name}")
    gates = {row["gate_id"]: row for row in read_table(candidate / "gate_status.tsv", "tsv")}
    expected_gates = {
        "GEN02_COMPLETE_SOURCE_TESTED_UNIVERSE": "coverage_limited",
        "GEN02_MATCHED_SENSITIVITY": "skipped_failed_entry_gate",
        "GEN03_HONG_SOURCE_AUDIT": "not_started_outside_preflight_scope",
        "NEGATIVE_CLAIM_AUTHORIZATION": "prohibited",
        "PREFLIGHT_OVERALL": "preflight_complete_claim_integration_blocked",
    }
    for gate, expected in expected_gates.items():
        if gates.get(gate, {}).get("status") != expected:
            raise ContractError(f"gate drift for {gate}")
    return [
        {
            "check_id": "fail_closed_source_power",
            "status": "pass",
            "observed": "coverage_limited_and_skipped",
            "expected": "coverage_limited_and_skipped",
            "note": "No source-tested negative, matched contrast, Hong rescue, or rescued fraction was emitted.",
        }
    ]


def check_prohibited(candidate: Path) -> list[dict[str, object]]:
    manifest_text = (candidate / "input_manifest.tsv").read_text(encoding="utf-8")
    quarantined = (
        "RNA-seq/results/causal_inference/sceqtl/ieqtl_disease_genes.csv",
        "RNA-seq/results/causal_inference/sceqtl/ieqtl_summary.csv",
    )
    if any(item in manifest_text for item in quarantined):
        raise ContractError("quarantined 9,007/4,703 local ieQTL product entered manifest")
    banned_phrases = (
        "causal_vs_reactive",
        "absence_of_genetic_biology",
        "ancestry_matched_validation",
        "powered_genetic_null",
    )
    offenders = []
    for path in candidate.glob("*.tsv"):
        content = path.read_text(encoding="utf-8").lower()
        if any(phrase in content for phrase in banned_phrases):
            offenders.append(path.name)
    if offenders:
        raise ContractError("prohibited semantic claims in: " + ", ".join(offenders))
    return [
        {
            "check_id": "prohibited_inputs_and_semantics",
            "status": "pass",
            "observed": "none",
            "expected": "none",
            "note": "Quarantined local ieQTL outputs and prohibited claim strings are absent.",
        }
    ]


def build_manifest(candidate: Path, manifest_path: Path) -> None:
    rows = []
    for path in sorted(candidate.rglob("*")):
        if not path.is_file() or path == manifest_path:
            continue
        relative = path.relative_to(candidate)
        if str(relative).startswith("work/input_snapshots/"):
            role = "immutable_input_snapshot"
        elif str(relative).startswith("work/stage_seals/"):
            role = "execution_seal"
        elif relative.parts[0] == "work":
            role = "work_product"
        else:
            role = "preflight_result"
        rows.append(
            {
                "relative_path": str(relative),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "role": role,
                "canonical_promotion_status": "not_promoted",
            }
        )
    atomic_write_tsv(
        manifest_path,
        rows,
        ["relative_path", "sha256", "bytes", "role", "canonical_promotion_status"],
    )


def build_producer_manifest(path: Path) -> None:
    allowed_suffixes = {".py", ".R", ".sbatch", ".sh", ".tsv", ".csv", ".md"}
    rows = []
    for source in sorted(SCRIPT_DIR.rglob("*")):
        if not source.is_file() or source.suffix not in allowed_suffixes:
            continue
        rows.append(
            {
                "relative_path": str(source.relative_to(PROJECT_ROOT)),
                "sha256": sha256_file(source),
                "bytes": source.stat().st_size,
                "role": (
                    "test_or_fixture"
                    if "tests" in source.relative_to(SCRIPT_DIR).parts
                    else "producer_or_contract"
                ),
            }
        )
    atomic_write_tsv(
        path,
        rows,
        ["relative_path", "sha256", "bytes", "role"],
    )


def main() -> None:
    candidate = assert_candidate_root(PROJECT_ROOT, DEFAULT_CANDIDATE_ROOT)
    stage4 = candidate / "work" / "stage_seals" / "04_power_interface.json"
    if not stage4.is_file():
        raise ContractError("power-interface stage seal is missing")
    source_paths = inputs(candidate)
    contract = independent_contract(source_paths)
    source = independent_source(source_paths["broadaway_leads"])
    require_equal(contract, EXPECTED_CONTRACT, "frozen primary contract")
    require_equal(source, EXPECTED_SOURCE, "Broadaway source-positive deposit")

    report_rows = [
        *[
            {
                "check_id": f"frozen_contract:{metric}",
                "status": "pass",
                "observed": contract[metric],
                "expected": expected,
                "note": "independently re-derived from pinned canonical input",
            }
            for metric, expected in EXPECTED_CONTRACT.items()
        ],
        *[
            {
                "check_id": f"source_reproduction:{metric}",
                "status": "pass",
                "observed": source[metric],
                "expected": expected,
                "note": "independently re-derived from pinned Broadaway significant-lead deposit",
            }
            for metric, expected in EXPECTED_SOURCE.items()
        ],
        *check_phenotypes(candidate),
        *check_fail_closed(candidate),
        *check_prohibited(candidate),
    ]
    report_path = candidate / "validation_report.tsv"
    atomic_write_tsv(
        report_path,
        report_rows,
        ["check_id", "status", "observed", "expected", "note"],
    )
    producer_manifest = candidate / "producer_manifest.tsv"
    build_producer_manifest(producer_manifest)
    manifest_path = candidate / "preflight_manifest.tsv"
    build_manifest(candidate, manifest_path)
    write_stage_seal(
        candidate,
        "05_validate_preflight",
        [report_path, producer_manifest, manifest_path],
        [stage4],
    )
    print(
        "PASS: genetics source/power preflight validated; "
        "claim-bearing source-tested/context branches remain blocked"
    )


if __name__ == "__main__":
    main()

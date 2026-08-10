#!/usr/bin/env python3
"""Independent SP-INT-05 validation for protein and ATAC candidate lanes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import math
import statistics
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from contract_lib import (
    ADAPTER_REGISTRY_COLUMNS,
    Check,
    ContractError,
    RELEASE_ID,
    SCHEMAS,
    parse_bool,
    parse_float,
    parse_int,
    program_universe,
    read_tsv,
    sha256_file,
    validate_candidate_contract,
    validate_hotspot_seal,
    write_tsv,
)


REPORT_COLUMNS = ("scope", "check_id", "status", "detail")
BUNDLE_COLUMNS = ("path_scope", "relative_path", "bytes", "sha256", "role")
STATUS_COLUMNS = (
    "release_id",
    "status",
    "protein_program_status",
    "fixed_25_status",
    "gse281367_status",
    "gse244832_status",
    "static_atac_status",
    "adapter_registry_sha256",
    "native_input_manifest_sha256",
    "bundle_manifest_sha256",
    "validation_report_sha256",
    "v1_preservation_unchanged",
    "proteomics_imported",
    "atac_imported",
    "plan13_complete",
    "validated_at_utc",
)
READY_COLUMNS = (
    "release_id",
    "status",
    "adapter_registry_sha256",
    "native_input_manifest_sha256",
    "bundle_manifest_sha256",
    "validation_report_sha256",
    "status_sha256",
    "builder_sha256",
    "native_producer_sha256",
    "validator_sha256",
    "n_adapters",
    "n_dynamic_program_effects",
    "n_fixed_protein_rows",
    "n_static_context_rows",
    "v1_preservation_unchanged",
    "proteomics_imported",
    "atac_imported",
    "plan13_complete",
    "validated_at_utc",
)


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ContractError(f"missing CSV header: {path}")
        return [dict(row) for row in reader]


def exact_two_group(score_by_id: Mapping[str, float], case_ids: set[str]) -> tuple[float, float, float]:
    ids = sorted(score_by_id)
    values = [float(score_by_id[sample]) for sample in ids]
    n = len(ids)
    n_case = len(case_ids)
    if n < 6 or n_case == 0 or n_case == n or not case_ids <= set(ids):
        raise ContractError("invalid exact two-group design")
    case = [score_by_id[sample] for sample in ids if sample in case_ids]
    control = [score_by_id[sample] for sample in ids if sample not in case_ids]
    effect = statistics.mean(case) - statistics.mean(control)
    total = sum(values)
    exceed = 0
    n_perm = 0
    for combo in itertools.combinations(range(n), n_case):
        case_sum = sum(values[index] for index in combo)
        perm_effect = case_sum / n_case - (total - case_sum) / (n - n_case)
        exceed += abs(perm_effect) >= abs(effect) - 1e-12
        n_perm += 1
    pvalue = exceed / n_perm
    std_error = math.sqrt(statistics.variance(case) / len(case) + statistics.variance(control) / len(control))
    return effect, pvalue, std_error


def bh_adjust(pvalues: Mapping[str, float]) -> dict[str, float]:
    if not pvalues:
        return {}
    ordered = sorted(pvalues.items(), key=lambda item: (item[1], item[0]))
    m = len(ordered)
    adjusted: dict[str, float] = {}
    running = 1.0
    for reverse_index in range(m - 1, -1, -1):
        key, value = ordered[reverse_index]
        rank = reverse_index + 1
        running = min(running, value * m / rank, 1.0)
        adjusted[key] = running
    return adjusted


def assert_complete_family(observed: Iterable[str], expected: Iterable[str], label: str) -> None:
    observed_values = list(observed)
    expected_values = list(expected)
    observed_set = set(observed_values)
    expected_set = set(expected_values)
    if (
        observed_set != expected_set
        or len(observed_values) != len(observed_set)
        or len(expected_values) != len(expected_set)
    ):
        raise ContractError(
            f"{label} is not the complete unique family: observed={sorted(observed_set)}, expected={sorted(expected_set)}"
        )


def assert_byte_identical(observed: Path, expected: Path, label: str) -> None:
    if observed.stat().st_size != expected.stat().st_size or sha256_file(observed) != sha256_file(expected):
        raise ContractError(f"{label} is not byte-identical")


def validate_absolute_input_manifest(path: Path, status_path: Path) -> list[Check]:
    _, rows = read_tsv(path, ("absolute_path", "bytes", "sha256", "input_role"))
    if not rows:
        raise ContractError("native ATAC input manifest is empty")
    for row in rows:
        source = Path(row["absolute_path"])
        if not source.is_absolute() or not source.is_file():
            raise ContractError(f"native input manifest points to missing/nonabsolute path: {source}")
        if parse_int(row["bytes"], f"input.bytes[{source}]") != source.stat().st_size:
            raise ContractError(f"native input byte drift: {source}")
        if row["sha256"] != sha256_file(source):
            raise ContractError(f"native input hash drift: {source}")
    _, status_rows = read_tsv(
        status_path,
        ("status", "input_manifest_sha256", "producer_sha256", "n_programs", "n_cohorts", "all_v1_regression_pass"),
    )
    if len(status_rows) != 1:
        raise ContractError("native ATAC producer status must contain one row")
    status = status_rows[0]
    if status["status"] != "passed_v1_selected_numeric_regression_and_v2_projection":
        raise ContractError("native ATAC producer did not pass its v1 regression")
    if status["input_manifest_sha256"] != sha256_file(path):
        raise ContractError("native ATAC input-manifest seal drift")
    producer = Path(__file__).resolve().parent / "10_run_atac_v2_projection.R"
    if status["producer_sha256"] != sha256_file(producer):
        raise ContractError("native ATAC producer hash drift")
    if status["n_programs"] != "2" or status["n_cohorts"] != "2" or not parse_bool(
        status["all_v1_regression_pass"], "native.all_v1_regression_pass"
    ):
        raise ContractError("native ATAC producer family/status drift")
    return [Check("native_input_manifest", "pass", f"{len(rows)} exact inputs rehashed, including external GTF")]


def _validate_native_regression(native_root: Path) -> list[Check]:
    _, rows = read_tsv(
        native_root / "v1_selected_regression.tsv",
        (
            "cohort", "legacy_program_id", "effect_abs_diff", "pvalue_abs_diff",
            "equal_effect_abs_diff", "leave_top_effect_abs_diff", "max_score_abs_diff",
            "regression_pass",
        ),
    )
    expected = {
        (cohort, program)
        for cohort in ("GSE281367", "GSE244832")
        for program in ("hepatocytes::8", "hepatocytes::20")
    }
    if {(row["cohort"], row["legacy_program_id"]) for row in rows} != expected or len(rows) != 4:
        raise ContractError("native v1 selected regression does not cover two programs x two cohorts")
    fields = (
        "effect_abs_diff", "pvalue_abs_diff", "equal_effect_abs_diff",
        "leave_top_effect_abs_diff", "max_score_abs_diff",
    )
    for row in rows:
        if not parse_bool(row["regression_pass"], "v1.regression_pass"):
            raise ContractError("native v1 selected regression contains a failure")
        if any((parse_float(row[field], f"v1.{field}") or 0.0) > 1e-10 for field in fields):
            raise ContractError("native v1 selected regression exceeds numeric tolerance")
    return [Check("atac_v1_numeric_regression", "pass", "4/4 selected cohort-program rows and donor scores reproduce")]


def _mapped_membership(hotspot_root: Path, uids: set[str]) -> dict[str, dict[str, float]]:
    _, rows = read_tsv(
        hotspot_root / "program_membership_v2.tsv",
        ("program_uid", "mapped_symbol", "mapped_symbol_status", "original_l1_weight"),
    )
    result = {uid: {} for uid in uids}
    for row in rows:
        uid = row["program_uid"]
        if uid not in result:
            continue
        symbol = row["mapped_symbol"].strip()
        if not symbol or symbol.upper() == "NA" or row["mapped_symbol_status"] == "unmapped_or_ambiguous":
            continue
        result[uid][symbol] = result[uid].get(symbol, 0.0) + float(row["original_l1_weight"])
    return result


def _validate_protein(project_root: Path, candidate_root: Path, hotspot_root: Path) -> list[Check]:
    root = candidate_root / "protein_atac/adapters/pxd051911_diams_program"
    fixed_source = project_root / "Analysis/Multimodal_Program_Projection/results/proteomics/panel4c_mrna_protein.tsv"
    fixed_candidate = root / "fixed_25_protein_context.tsv"
    assert_byte_identical(fixed_candidate, fixed_source, "fixed 25-protein descriptive table")
    _, fixed_rows = read_tsv(
        fixed_candidate,
        ("gene", "selection_conditioned", "interpretation"),
    )
    if len(fixed_rows) != 25 or len({row["gene"] for row in fixed_rows}) != 25:
        raise ContractError("fixed 25-protein table does not contain 25 unique rows")
    if any(
        not parse_bool(row["selection_conditioned"], f"fixed25.selection_conditioned[{row['gene']}]")
        or row["interpretation"] != "descriptive_same_cohort_reestimate"
        for row in fixed_rows
    ):
        raise ContractError("fixed 25-protein rows lost selection-conditioned descriptive semantics")

    seal = validate_hotspot_seal(hotspot_root)
    universe = program_universe(seal, "external_test_eligible")
    membership = _mapped_membership(hotspot_root, set(universe))
    _, de_rows = read_tsv(
        project_root / "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv",
        ("gene",),
    )
    measured = {row["gene"] for row in de_rows}
    _, mapping_rows = read_tsv(root / "gene_mapping_audit.tsv", SCHEMAS["gene_mapping_audit.tsv"])
    mapping = {row["program_uid"]: row for row in mapping_rows}
    _, effects = read_tsv(root / "program_effects.tsv", SCHEMAS["program_effects.tsv"])
    _, scores = read_tsv(root / "per_sample_program_scores.tsv", SCHEMAS["per_sample_program_scores.tsv"])
    assert_complete_family([row["program_uid"] for row in effects], universe, "protein effect family")
    if scores:
        raise ContractError("untestable protein adapter contains program scores")
    for uid in universe:
        expected_genes = set(membership[uid]) & measured
        expected_weight = sum(membership[uid][gene] for gene in expected_genes)
        row = mapping[uid]
        if parse_int(row["n_genes_measured"], f"protein.n_genes[{uid}]") != len(expected_genes):
            raise ContractError(f"protein measured-gene count does not rederive for {uid}")
        observed_weight = parse_float(row["retained_l1_weight"], f"protein.weight[{uid}]")
        if observed_weight is None or abs(observed_weight - expected_weight) > 1e-12:
            raise ContractError(f"protein retained weight does not rederive for {uid}")
    for row in effects:
        if row["evidence_state"] != "untestable" or row["testable"] != "FALSE":
            raise ContractError("protein program row is not terminal untestable")
        if any(row[field].strip() for field in ("estimate", "std_error", "pvalue", "padj")):
            raise ContractError("untestable protein program carries a numeric outcome")
    return [
        Check("protein_program_gate", "pass", "2/2 v2 programs untestable by prespecified DIA coverage; no outcome test"),
        Check("fixed_25_protein_context", "pass", "25 byte-identical, unique, selection-conditioned descriptive rows"),
    ]


def _rows_by_uid(path: Path) -> dict[str, dict[str, str]]:
    _, rows = read_tsv(path, ("program_uid",))
    result = {row["program_uid"]: row for row in rows}
    if len(result) != len(rows):
        raise ContractError(f"duplicate program_uid in {path}")
    return result


def _validate_atac_cohort(candidate_root: Path, hotspot_root: Path, cohort: str) -> list[Check]:
    root = candidate_root / f"protein_atac/adapters/{cohort.lower()}_snatac_dynamic"
    seal = validate_hotspot_seal(hotspot_root)
    universe = program_universe(seal, "external_test_eligible")
    effects = _rows_by_uid(root / "program_effects.tsv")
    assert_complete_family(effects, universe, f"{cohort} ATAC effect family")
    _, native_scores = read_tsv(
        root / "scoring_sensitivity_scores.tsv",
        (
            "program_uid", "sample_id", "condition", "score", "equal_score",
            "leave_top_score", "cohort", "membership_sha256",
        ),
    )
    score_by_uid: dict[str, list[dict[str, str]]] = {uid: [] for uid in universe}
    for row in native_scores:
        if row["cohort"] != cohort or row["program_uid"] not in score_by_uid:
            raise ContractError(f"{cohort} sensitivity score references wrong cohort/program")
        score_by_uid[row["program_uid"]].append(row)
    pvalues = {}
    rederived = {}
    for uid in sorted(universe):
        rows = score_by_uid[uid]
        expected_n = 12 if cohort == "GSE281367" else 14
        if len(rows) != expected_n or len({row["sample_id"] for row in rows}) != expected_n:
            raise ContractError(f"{cohort} incomplete donor sensitivity scores for {uid}")
        cases = {row["sample_id"] for row in rows if row["condition"] == "MASH"}
        primary = {row["sample_id"]: float(row["score"]) for row in rows}
        equal = {row["sample_id"]: float(row["equal_score"]) for row in rows}
        leave = {row["sample_id"]: float(row["leave_top_score"]) for row in rows}
        effect, pvalue, std_error = exact_two_group(primary, cases)
        equal_effect, _, _ = exact_two_group(equal, cases)
        leave_effect, _, _ = exact_two_group(leave, cases)
        pvalues[uid] = pvalue
        rederived[uid] = (effect, std_error, equal_effect, leave_effect)
    qvalues = bh_adjust(pvalues)
    n_robust = 0
    for uid in sorted(universe):
        row = effects[uid]
        effect, std_error, equal_effect, leave_effect = rederived[uid]
        observed = {
            "estimate": effect,
            "std_error": std_error,
            "pvalue": pvalues[uid],
            "padj": qvalues[uid],
        }
        for field, expected in observed.items():
            value = parse_float(row[field], f"{cohort}.{field}[{uid}]")
            if value is None or abs(value - expected) > 1e-12:
                raise ContractError(f"{cohort} {field} does not rederive for {uid}: {value} != {expected}")
        direction = "positive" if effect > 0 else "negative" if effect < 0 else "zero"
        sign_agree = all(
            ("positive" if value > 0 else "negative" if value < 0 else "zero") == direction
            for value in (equal_effect, leave_effect)
        )
        robust = qvalues[uid] < 0.05 and sign_agree
        expected_state = "robust" if robust else "indeterminate"
        if row["direction_observed"] != direction:
            raise ContractError(f"{cohort} direction drift for {uid}")
        if row["descriptive_effect_direction"] != direction or row["inferential_test_direction"] != direction:
            raise ContractError(f"{cohort} descriptive/inferential direction drift for {uid}")
        if not parse_bool(row["direction_agreement"], f"{cohort}.direction_agreement[{uid}]"):
            raise ContractError(f"{cohort} direction agreement unexpectedly false for {uid}")
        if row["negative_call_rule_id"].strip():
            raise ContractError(f"{cohort} indeterminate row carries a negative-call rule for {uid}")
        if parse_bool(row["sensitivity_sign_agree"], f"{cohort}.sign_agree[{uid}]") != sign_agree:
            raise ContractError(f"{cohort} sensitivity sign drift for {uid}")
        if parse_bool(row["robustness_pass"], f"{cohort}.robust[{uid}]") != robust:
            raise ContractError(f"{cohort} robust flag drift for {uid}")
        if row["evidence_state"] != expected_state:
            raise ContractError(f"{cohort} evidence state drift for {uid}")
        n_robust += int(robust)
        _, sensitivities = read_tsv(root / "sensitivity.tsv", SCHEMAS["sensitivity.tsv"])
        this_sensitivity = {r["sensitivity_id"]: r for r in sensitivities if r["program_uid"] == uid}
        for sensitivity_id, expected in (
            ("equal_weight", equal_effect),
            ("leave_highest_weight_gene_out", leave_effect),
        ):
            if sensitivity_id not in this_sensitivity:
                raise ContractError(f"{cohort} missing {sensitivity_id} for {uid}")
            value = parse_float(this_sensitivity[sensitivity_id]["estimate"], f"{cohort}.{sensitivity_id}[{uid}]")
            if value is None or abs(value - expected) > 1e-12:
                raise ContractError(f"{cohort} {sensitivity_id} does not rederive for {uid}")
    return [
        Check(
            f"{cohort.lower()}_dynamic_atac",
            "pass",
            f"2-program BH family; donor exact permutations rederived; {n_robust} robust (positivity not required)",
        )
    ]


def _validate_static(project_root: Path, candidate_root: Path, hotspot_root: Path) -> list[Check]:
    root = candidate_root / "protein_atac/source_native/atac_static_context"
    _, contexts = read_tsv(
        root / "static_atac_context.tsv",
        (
            "program_uid", "legacy_program_id", "trait_scope", "n_source_loci",
            "n_study_loci", "n_variants", "n_genes", "n_studies", "n_ancestries",
            "effect_unit", "source_dependence", "interpretation",
        ),
    )
    _, hits = read_tsv(
        root / "static_atac_hits.tsv",
        (
            "program_uid", "legacy_program_id", "gene_symbol", "study", "ancestry",
            "trait_scope", "source_locus", "variant_id", "source_dependence", "interpretation",
        ),
    )
    seal = validate_hotspot_seal(hotspot_root)
    universe = program_universe(seal, "external_test_eligible")
    expected_keys = {(uid, scope) for uid in universe for scope in ("direct_disease_PDFF", "liver_enzyme")}
    if {(row["program_uid"], row["trait_scope"]) for row in contexts} != expected_keys or len(contexts) != 4:
        raise ContractError("static ATAC context is not the complete two-program x two-scope grid")
    for row in contexts:
        uid, scope = row["program_uid"], row["trait_scope"]
        subset = [hit for hit in hits if hit["program_uid"] == uid and hit["trait_scope"] == scope]
        expected = {
            "n_source_loci": len({hit["source_locus"] for hit in subset}),
            "n_study_loci": len({(hit["study"], hit["source_locus"]) for hit in subset}),
            "n_variants": len({hit["variant_id"] for hit in subset}),
            "n_genes": len({hit["gene_symbol"] for hit in subset}),
            "n_studies": len({hit["study"] for hit in subset}),
            "n_ancestries": len({hit["ancestry"] for hit in subset}),
        }
        for field, value in expected.items():
            if parse_int(row[field], f"static.{field}[{uid}/{scope}]") != value:
                raise ContractError(f"static ATAC {field} does not rederive for {uid}/{scope}")
        if row["source_dependence"] != "descriptive" or "not_dynamic_accessibility" not in row["interpretation"]:
            raise ContractError("static ATAC context lost descriptive/non-dynamic semantics")
    if any(row["source_dependence"] != "descriptive" for row in hits):
        raise ContractError("static ATAC hit row lost descriptive source label")
    header, _ = read_tsv(root / "static_atac_context.tsv")
    if any(column in header for column in ("pvalue", "padj", "robust", "testable")):
        raise ContractError("static descriptive context invented inferential columns")
    assert_byte_identical(
        root / "static_study_coverage.tsv",
        project_root / "Analysis/Multimodal_Program_Projection/results/atac/static_study_coverage.tsv",
        "static 35-study coverage",
    )
    return [Check("static_atac_context", "pass", f"4 descriptive scope rows and {len(hits)} source-native gene-link rows rederived")]


def _validate_v1(candidate_root: Path) -> list[Check]:
    _, rows = read_tsv(
        candidate_root / "v1_preservation_summary.tsv",
        ("n_unchanged", "n_changed", "n_missing", "n_added", "aggregate_sha256"),
    )
    if len(rows) != 1:
        raise ContractError("v1 preservation summary must contain one row")
    row = rows[0]
    if row["n_unchanged"] != "591" or any(row[field] != "0" for field in ("n_changed", "n_missing", "n_added")):
        raise ContractError("v1 preservation failed during SP-INT-05")
    return [Check("v1_preservation", "pass", f"591/591 unchanged; aggregate={row['aggregate_sha256']}")]


def run_checks(project_root: Path, candidate_root: Path, hotspot_root: Path) -> list[Check]:
    native_root = candidate_root / "protein_atac_native/atac_v2"
    checks = validate_absolute_input_manifest(native_root / "input_manifest.tsv", native_root / "producer_status.tsv")
    checks.extend(_validate_native_regression(native_root))
    checks.extend(
        validate_candidate_contract(
            candidate_root,
            candidate_root / "protein_atac/adapter_registry.tsv",
            hotspot_root,
            project_root,
        )
    )
    _, registry = read_tsv(
        candidate_root / "protein_atac/adapter_registry.tsv",
        ADAPTER_REGISTRY_COLUMNS,
    )
    expected_adapters = {
        "pxd051911_diams_program",
        "gse281367_snatac_dynamic",
        "gse244832_snatac_dynamic",
    }
    if {row["adapter_id"] for row in registry} != expected_adapters or len(registry) != 3:
        raise ContractError("SP-INT-05 adapter registry is not the bounded three-adapter family")
    checks.append(Check("bounded_sp_int_05_registry", "pass", "exact protein + two dynamic ATAC adapters"))
    checks.extend(_validate_protein(project_root, candidate_root, hotspot_root))
    checks.extend(_validate_atac_cohort(candidate_root, hotspot_root, "GSE281367"))
    checks.extend(_validate_atac_cohort(candidate_root, hotspot_root, "GSE244832"))
    checks.extend(_validate_static(project_root, candidate_root, hotspot_root))
    checks.extend(_validate_v1(candidate_root))
    return checks


def _bundle_manifest(candidate_root: Path) -> list[dict[str, object]]:
    roots = (
        (candidate_root / "protein_atac", "candidate_relative", "sp_int_05_adapter_bundle"),
        (candidate_root / "protein_atac_native/atac_v2", "candidate_relative", "sp_int_05_native_atac"),
    )
    rows = []
    for root, scope, role in roots:
        for path in sorted(root.rglob("*")):
            if not path.is_file() or path.name in {"bundle_manifest.tsv", "validation.tsv", "status.tsv", "READY"}:
                continue
            rows.append(
                {
                    "path_scope": scope,
                    "relative_path": path.relative_to(candidate_root).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "role": role,
                }
            )
    if len(rows) < 40:
        raise ContractError(f"SP-INT-05 bundle unexpectedly small: {len(rows)} files")
    return rows


def _write_outputs(project_root: Path, candidate_root: Path, checks: Sequence[Check]) -> None:
    root = candidate_root / "protein_atac"
    outputs = {
        "bundle": root / "bundle_manifest.tsv",
        "report": root / "validation.tsv",
        "status": root / "status.tsv",
        "ready": root / "READY",
    }
    existing = [str(path) for path in outputs.values() if path.exists()]
    if existing:
        raise ContractError("refusing to overwrite SP-INT-05 validation outputs: " + ", ".join(existing))
    write_tsv(outputs["bundle"], BUNDLE_COLUMNS, _bundle_manifest(candidate_root))
    write_tsv(
        outputs["report"],
        REPORT_COLUMNS,
        [
            {"scope": "SP-INT-05", "check_id": row.check_id, "status": row.status, "detail": row.detail}
            for row in checks
        ],
    )
    native_manifest = candidate_root / "protein_atac_native/atac_v2/input_manifest.tsv"
    validated = utc_now()
    status_row = {
        "release_id": RELEASE_ID,
        "status": "protein_atac_ready_candidate_only_plan13_incomplete",
        "protein_program_status": "untestable_two_of_two",
        "fixed_25_status": "descriptive_selection_conditioned_byte_identical",
        "gse281367_status": "donor_dynamic_test_complete",
        "gse244832_status": "source_dependent_donor_dynamic_test_complete",
        "static_atac_status": "descriptive_gene_link_context_complete",
        "adapter_registry_sha256": sha256_file(root / "adapter_registry.tsv"),
        "native_input_manifest_sha256": sha256_file(native_manifest),
        "bundle_manifest_sha256": sha256_file(outputs["bundle"]),
        "validation_report_sha256": sha256_file(outputs["report"]),
        "v1_preservation_unchanged": True,
        "proteomics_imported": True,
        "atac_imported": True,
        "plan13_complete": False,
        "validated_at_utc": validated,
    }
    write_tsv(outputs["status"], STATUS_COLUMNS, [status_row])
    _, dynamic_281 = read_tsv(root / "adapters/gse281367_snatac_dynamic/program_effects.tsv", SCHEMAS["program_effects.tsv"])
    _, dynamic_244 = read_tsv(root / "adapters/gse244832_snatac_dynamic/program_effects.tsv", SCHEMAS["program_effects.tsv"])
    _, fixed_rows = read_tsv(root / "adapters/pxd051911_diams_program/fixed_25_protein_context.tsv")
    _, static_rows = read_tsv(root / "source_native/atac_static_context/static_atac_context.tsv")
    ready_row = {
        "release_id": RELEASE_ID,
        "status": "ready_sp_int_05_candidate_only_not_plan13_complete",
        "adapter_registry_sha256": status_row["adapter_registry_sha256"],
        "native_input_manifest_sha256": status_row["native_input_manifest_sha256"],
        "bundle_manifest_sha256": status_row["bundle_manifest_sha256"],
        "validation_report_sha256": status_row["validation_report_sha256"],
        "status_sha256": sha256_file(outputs["status"]),
        "builder_sha256": sha256_file(Path(__file__).resolve().parent / "11_build_protein_atac_adapters.py"),
        "native_producer_sha256": sha256_file(Path(__file__).resolve().parent / "10_run_atac_v2_projection.R"),
        "validator_sha256": sha256_file(Path(__file__).resolve()),
        "n_adapters": 3,
        "n_dynamic_program_effects": len(dynamic_281) + len(dynamic_244),
        "n_fixed_protein_rows": len(fixed_rows),
        "n_static_context_rows": len(static_rows),
        "v1_preservation_unchanged": True,
        "proteomics_imported": True,
        "atac_imported": True,
        "plan13_complete": False,
        "validated_at_utc": validated,
    }
    write_tsv(outputs["ready"], READY_COLUMNS, [ready_row])


def _check_ready(candidate_root: Path) -> None:
    root = candidate_root / "protein_atac"
    _, rows = read_tsv(root / "READY", READY_COLUMNS)
    if len(rows) != 1:
        raise ContractError("SP-INT-05 READY must contain one row")
    row = rows[0]
    expected = {
        "adapter_registry_sha256": sha256_file(root / "adapter_registry.tsv"),
        "native_input_manifest_sha256": sha256_file(candidate_root / "protein_atac_native/atac_v2/input_manifest.tsv"),
        "bundle_manifest_sha256": sha256_file(root / "bundle_manifest.tsv"),
        "validation_report_sha256": sha256_file(root / "validation.tsv"),
        "status_sha256": sha256_file(root / "status.tsv"),
        "builder_sha256": sha256_file(Path(__file__).resolve().parent / "11_build_protein_atac_adapters.py"),
        "native_producer_sha256": sha256_file(Path(__file__).resolve().parent / "10_run_atac_v2_projection.R"),
        "validator_sha256": sha256_file(Path(__file__).resolve()),
    }
    if row["release_id"] != RELEASE_ID or row["status"] != "ready_sp_int_05_candidate_only_not_plan13_complete":
        raise ContractError("SP-INT-05 READY status drift")
    for field, value in expected.items():
        if row[field] != value:
            raise ContractError(f"SP-INT-05 READY hash mismatch: {field}")
    if row["n_adapters"] != "3" or row["n_dynamic_program_effects"] != "4" or row["n_fixed_protein_rows"] != "25" or row["n_static_context_rows"] != "4":
        raise ContractError("SP-INT-05 READY count drift")
    if not parse_bool(row["v1_preservation_unchanged"], "READY.v1"):
        raise ContractError("SP-INT-05 READY does not preserve v1")
    if not parse_bool(row["proteomics_imported"], "READY.protein") or not parse_bool(row["atac_imported"], "READY.atac"):
        raise ContractError("SP-INT-05 READY omits completed lanes")
    if parse_bool(row["plan13_complete"], "READY.plan13"):
        raise ContractError("SP-INT-05 READY overstates Plan 13 completion")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--candidate-root", type=Path, default=None)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    candidate_root = args.candidate_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "spatial_context"
    )
    hotspot_root = project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "hotspot"
    try:
        checks = run_checks(project_root, candidate_root, hotspot_root)
        if any(row.status != "pass" for row in checks):
            raise ContractError("one or more SP-INT-05 checks did not pass")
        if args.check_only:
            _check_ready(candidate_root)
            print(f"PASS: rederived {len(checks)} SP-INT-05 checks and READY hashes")
            return 0
        _write_outputs(project_root, candidate_root, checks)
        _check_ready(candidate_root)
        print(f"PASS: {len(checks)} SP-INT-05 checks; protein/ATAC candidate READY")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

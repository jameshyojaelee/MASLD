#!/usr/bin/env python3
"""Independently validate and seal the candidate true-Kleiner LVQW bundle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


CANDIDATE_ID = "program-context-v2-candidate-2026-08-07"
ANALYSIS_ID = "fibrosis-adjacent-true-kleiner-lvqw-v1"
ALLOWLIST = (
    "GSE130970",
    "GSE135251",
    "GSE162694",
    "GSE174478",
    "GSE193066",
    "GSE240729",
)
EXCLUDED = {"GSE213621", "PRJNA512027"}
EXPECTED = {
    "F1_vs_F0": ("F0_to_F1", 0, 1, 126, 187, 313, 6),
    "F2_vs_F1": ("F1_to_F2", 1, 2, 187, 174, 361, 6),
    "F3_vs_F2": ("F2_to_F3", 2, 3, 174, 133, 307, 6),
    "F4_vs_F3": ("F3_to_F4", 3, 4, 108, 44, 152, 5),
}
EXPECTED_COHORTS = {
    "F0_to_F1": set(ALLOWLIST),
    "F1_to_F2": set(ALLOWLIST),
    "F2_to_F3": set(ALLOWLIST),
    "F3_to_F4": set(ALLOWLIST) - {"GSE193066"},
}
EXPECTED_N_GENES = 27638
METHOD = "limma_voom_qw_C2_true_kleiner_first_biopsy"
MODEL = "~ dataset + inferred_sex + fib_group"
REQUIRED_PRODUCER_FILES = {
    "results/fibrosis_consecutive_lvqw_true_kleiner.csv",
    "results/transition_summary.tsv",
    "audits/null_shuffle_results.csv",
    "audits/null_shuffle_summary.tsv",
    "audits/sample_manifest.tsv",
    "audits/cohort_arm_audit.tsv",
    "audits/design_audit.tsv",
    "audits/model_design.tsv",
    "audits/shuffle_audit.tsv",
    "audits/gse193066_donor_biopsy_crosswalk.tsv",
    "manifests/input_manifest.tsv",
    "manifests/code_manifest.tsv",
    "manifests/environment_manifest.tsv",
    "manifests/sessionInfo.txt",
    "manifests/run_contract.tsv",
    "manifests/producer_status.tsv",
}


class ValidationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_rows(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def read_key_values(path: Path) -> dict[str, str]:
    rows = read_rows(path)
    require(all(set(row) == {"key", "value"} for row in rows), "invalid key/value table")
    values = {row["key"]: row["value"] for row in rows}
    require(len(values) == len(rows), "duplicate run-contract key")
    return values


def parse_int(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise ValidationError(f"{label} is not an integer: {value!r}") from error
    return parsed


def parse_float(value: str, label: str) -> float:
    try:
        parsed = float(value)
    except ValueError as error:
        raise ValidationError(f"{label} is not numeric: {value!r}") from error
    require(math.isfinite(parsed), f"{label} is not finite")
    return parsed


def parse_bool(value: str, label: str) -> bool:
    normalized = value.strip().lower()
    require(normalized in {"true", "false"}, f"{label} is not TRUE/FALSE: {value!r}")
    return normalized == "true"


def gse193066_donor_id(source_title: str) -> str:
    require(
        re.fullmatch(r"HUnafld[0-9]{3}(?:_[12])?", source_title) is not None,
        f"invalid GSE193066 source title: {source_title!r}",
    )
    return "GSE193066::" + re.sub(r"_[12]$", "", source_title)


def bh_adjust(p_values: list[float]) -> list[float]:
    n_values = len(p_values)
    order = sorted(range(n_values), key=lambda index: p_values[index])
    adjusted = [1.0] * n_values
    running = 1.0
    for reverse_index in range(n_values - 1, -1, -1):
        original_index = order[reverse_index]
        rank = reverse_index + 1
        running = min(running, p_values[original_index] * n_values / rank, 1.0)
        adjusted[original_index] = running
    return adjusted


def write_exclusive(path: Path, payload: str) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o640)
    with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
        handle.write(payload)


def tsv_payload(rows: list[dict[str, object]], fields: list[str]) -> str:
    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def expected_root(project: Path) -> Path:
    return (
        project
        / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/candidates"
        / CANDIDATE_ID
        / ANALYSIS_ID
    )


def validate(project: Path, root: Path) -> dict[str, object]:
    checks: list[dict[str, object]] = []

    def passed(check_id: str, detail: str) -> None:
        checks.append({"check_id": check_id, "status": "PASS", "detail": detail})

    require(root.exists() and not root.is_symlink(), "candidate root is missing or symlinked")
    project = project.resolve(strict=True)
    root = root.resolve(strict=True)
    require(root == expected_root(project).resolve(strict=True), "candidate root is not the fixed isolated path")
    require(not root.is_symlink(), "candidate root is a symlink")
    for path in root.rglob("*"):
        require(not path.is_symlink(), f"symlink is prohibited in candidate bundle: {path}")
    passed("V01_FIXED_ROOT", "candidate root and every component are non-symlinked")

    for forbidden in ("READY", "manifests/validation_checks.tsv", "manifests/validation_report.json", "manifests/validated_artifact_manifest.tsv"):
        require(not (root / forbidden).exists(), f"validator refuses existing output: {forbidden}")

    producer_manifest_path = root / "manifests/producer_artifact_manifest.tsv"
    require(producer_manifest_path.is_file(), "producer artifact manifest is missing")
    producer_rows = read_rows(producer_manifest_path)
    require(producer_rows, "producer artifact manifest is empty")
    manifest_paths = {row["artifact_path"] for row in producer_rows}
    require(len(manifest_paths) == len(producer_rows), "duplicate producer artifact path")
    require(manifest_paths == REQUIRED_PRODUCER_FILES, "producer artifact set is not exact")
    observed_prevalidation = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path != producer_manifest_path
    }
    require(observed_prevalidation == manifest_paths, "unmanifested or missing producer artifact")
    for row in producer_rows:
        path = root / row["artifact_path"]
        require(path.is_file(), f"missing producer artifact: {path}")
        require(path.stat().st_size == parse_int(row["bytes"], f"{path} bytes"), f"byte drift: {path}")
        require(sha256(path) == row["sha256"], f"SHA256 drift: {path}")
    passed("V02_PRODUCER_MANIFEST", f"{len(producer_rows)} exact producer artifacts verified")

    source_paths: dict[str, Path] = {}
    for manifest_name, expected_ids in (
        (
            "input_manifest.tsv",
            {"merged_dge", "meta_matched", "sample_qc", "gse193066_source_metadata"},
        ),
        ("code_manifest.tsv", {"producer", "validator", "slurm_wrapper", "lvqw_engine"}),
    ):
        rows = read_rows(root / "manifests" / manifest_name)
        require({row["artifact_id"] for row in rows} == expected_ids, f"{manifest_name} IDs drift")
        for row in rows:
            source = (project / row["path"]).resolve(strict=True)
            require(source.is_file() and not source.is_symlink(), f"unsafe manifest source: {source}")
            require(project == source or project in source.parents, f"manifest source escapes project: {source}")
            require(source.stat().st_size == parse_int(row["bytes"], f"{source} bytes"), f"source byte drift: {source}")
            require(sha256(source) == row["sha256"], f"source SHA256 drift: {source}")
            require(row["artifact_id"] not in source_paths, "duplicate source artifact ID")
            source_paths[row["artifact_id"]] = source
    passed("V03_SOURCE_CODE_HASHES", "all input and code hashes match live frozen sources")

    contract = read_key_values(root / "manifests/run_contract.tsv")
    expected_contract = {
        "candidate_id": CANDIDATE_ID,
        "analysis_id": ANALYSIS_ID,
        "model": MODEL,
        "coefficient": "fib_grouphigh",
        "method": METHOD,
        "tested_universe": str(EXPECTED_N_GENES),
        "significance_rule": "BH padj < 0.05 per complete contrast family",
        "lfc_gate": "none",
        "allowlist": ";".join(ALLOWLIST),
        "hard_exclusions": "GSE213621;PRJNA512027",
        "gse193066_timepoint": "1st biopsy only",
        "gse193066_donor_key": "deposited !Sample_title with terminal _1/_2 removed",
        "shuffle_seed": "20260808",
        "replication_unit": "one first-biopsy biological donor/sample",
        "canonical_write": "false",
    }
    require(contract == expected_contract, "run contract differs from the prespecified contract")
    passed("V04_RUN_CONTRACT", "fixed-effect LVQW, complete BH family, and no LFC gate verified")

    environment = read_rows(root / "manifests/environment_manifest.tsv")
    require({row["component"] for row in environment} == {"R", "data.table", "edgeR", "limma", "ashr"}, "environment components drift")
    require(all(row["slurm_job_id"].isdigit() for row in environment), "producer did not run under SLURM")
    require(all(row["slurm_partition"] == "cpu" for row in environment), "unexpected SLURM partition")
    passed("V05_ENVIRONMENT", f"SLURM job {environment[0]['slurm_job_id']} on cpu with pinned package versions")

    # Independently re-derive the GSE193066 donor/timepoint key from the
    # deposited local metadata rather than trusting the producer crosswalk.
    gse_source_rows = read_rows(source_paths["gse193066_source_metadata"])
    required_source_fields = {
        "Run",
        "sample_id",
        "!Sample_title",
        "biopsy",
        "fibrosis stage",
    }
    require(len(gse_source_rows) == 164, "GSE193066 source row count drift")
    require(
        all(required_source_fields.issubset(row) for row in gse_source_rows),
        "GSE193066 source schema drift",
    )
    gse_source_by_sample: dict[str, dict[str, object]] = {}
    gse_source_by_donor: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in gse_source_rows:
        sample_id = row["Run"]
        require(sample_id == row["sample_id"] and sample_id, "GSE193066 Run/sample_id mismatch")
        require(sample_id not in gse_source_by_sample, "duplicate GSE193066 run")
        source_title = row["!Sample_title"]
        donor_id = gse193066_donor_id(source_title)
        biopsy = row["biopsy"]
        require(biopsy in {"1st biopsy", "2nd biopsy"}, "invalid GSE193066 biopsy label")
        require(not (biopsy == "1st biopsy" and source_title.endswith("_2")), "first biopsy has _2 title")
        require(biopsy != "2nd biopsy" or source_title.endswith("_2"), "second biopsy lacks _2 title")
        stage = parse_int(row["fibrosis stage"], "GSE193066 source fibrosis stage")
        require(stage in {0, 1, 2, 3, 4}, "invalid GSE193066 source fibrosis stage")
        parsed = {
            "sample_id": sample_id,
            "source_title": source_title,
            "donor_id": donor_id,
            "biopsy": biopsy,
            "stage": stage,
        }
        gse_source_by_sample[sample_id] = parsed
        gse_source_by_donor[donor_id].append(parsed)
    require(len(gse_source_by_donor) == 106, "GSE193066 donor count drift")
    require(
        sum(row["biopsy"] == "1st biopsy" for row in gse_source_by_sample.values()) == 106
        and sum(row["biopsy"] == "2nd biopsy" for row in gse_source_by_sample.values()) == 58,
        "GSE193066 first/second-biopsy census drift",
    )
    paired_donors = 0
    for donor_id, rows in gse_source_by_donor.items():
        first = [row for row in rows if row["biopsy"] == "1st biopsy"]
        second = [row for row in rows if row["biopsy"] == "2nd biopsy"]
        require(len(first) == 1 and len(second) <= 1, f"invalid biopsy multiplicity for {donor_id}")
        if second:
            paired_donors += 1
            require(first[0]["source_title"].endswith("_1"), f"paired first biopsy lacks _1 for {donor_id}")
    require(paired_donors == 58, "GSE193066 paired-donor count drift")

    qc_rows = read_rows(source_paths["sample_qc"], delimiter=",")
    qc_pass: dict[str, bool] = {}
    for row in qc_rows:
        sample_id = row["sample_id"]
        require(sample_id not in qc_pass, "duplicate sample in source QC table")
        qc_pass[sample_id] = parse_bool(row["pass_technical"], "source pass_technical")
    require(set(gse_source_by_sample).issubset(qc_pass), "GSE193066 source sample absent from QC table")

    crosswalk_rows = read_rows(root / "audits/gse193066_donor_biopsy_crosswalk.tsv")
    require(len(crosswalk_rows) == 164, "GSE193066 crosswalk row count drift")
    crosswalk_by_sample: dict[str, dict[str, str]] = {}
    included_gse_samples: set[str] = set()
    for row in crosswalk_rows:
        sample_id = row["sample_id"]
        require(sample_id not in crosswalk_by_sample, "duplicate GSE193066 crosswalk sample")
        require(sample_id in gse_source_by_sample, "crosswalk sample absent from deposited metadata")
        source = gse_source_by_sample[sample_id]
        require(row["source_title"] == source["source_title"], "crosswalk source-title drift")
        require(row["donor_id"] == source["donor_id"], "crosswalk donor-key drift")
        require(row["biopsy"] == source["biopsy"], "crosswalk biopsy drift")
        require(
            parse_int(row["source_fibrosis_stage"], "crosswalk source stage") == source["stage"]
            and parse_int(row["harmonized_fibrosis_stage"], "crosswalk harmonized stage") == source["stage"],
            "crosswalk fibrosis-stage drift",
        )
        is_first = parse_bool(row["is_first_biopsy"], "crosswalk is_first_biopsy")
        pass_technical = parse_bool(row["pass_technical"], "crosswalk pass_technical")
        dge_present = parse_bool(row["dge_present"], "crosswalk dge_present")
        included = parse_bool(row["included_primary"], "crosswalk included_primary")
        require(is_first == (source["biopsy"] == "1st biopsy"), "crosswalk first-biopsy flag drift")
        require(pass_technical == qc_pass[sample_id], "crosswalk QC flag drift")
        require(dge_present == pass_technical, "GSE193066 DGE/pass_technical universe drift")
        require(included == (is_first and pass_technical and dge_present), "crosswalk inclusion-rule drift")
        expected_reason = (
            "included_first_biopsy"
            if included
            else "excluded_second_biopsy"
            if not is_first
            else "excluded_failed_technical_qc"
        )
        require(row["exclusion_reason"] == expected_reason, "crosswalk exclusion-reason drift")
        if included:
            included_gse_samples.add(sample_id)
        crosswalk_by_sample[sample_id] = row
    require(set(crosswalk_by_sample) == set(gse_source_by_sample), "GSE193066 crosswalk/source universe drift")
    require(len(included_gse_samples) == 105, "GSE193066 first-biopsy pass-QC census drift")
    require(
        len({crosswalk_by_sample[sample]["donor_id"] for sample in included_gse_samples}) == 105,
        "GSE193066 primary donor IDs are not unique",
    )
    passed(
        "V06_GSE193066_DONORS",
        "164 source rows resolve to 106 donors; only 105 pass-QC first biopsies are eligible",
    )

    expected_by_transition = {values[0]: (contrast, *values[1:]) for contrast, values in EXPECTED.items()}
    sample_rows = read_rows(root / "audits/sample_manifest.tsv")
    samples: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    donors: dict[str, set[str]] = defaultdict(set)
    all_primary_donors: set[str] = set()
    observed_gse_samples: set[str] = set()
    for row in sample_rows:
        transition = row["transition"]
        require(transition in expected_by_transition, f"unknown sample transition: {transition}")
        require(row["contrast"] == expected_by_transition[transition][0], "sample contrast mapping drift")
        require(row["dataset"] in ALLOWLIST and row["dataset"] not in EXCLUDED, "invalid sample cohort")
        require(row["pass_technical"] == "TRUE", "non-pass_technical sample entered")
        require(row["sample_id"] not in samples[transition], "duplicate sample within transition")
        require(row["donor_id"] not in donors[transition], "duplicate biological donor within transition")
        if row["dataset"] == "GSE193066":
            require(row["sample_id"] in included_gse_samples, "ineligible GSE193066 biopsy entered model")
            source = gse_source_by_sample[row["sample_id"]]
            require(row["donor_id"] == source["donor_id"], "GSE193066 model donor ID drift")
            require(row["source_title"] == source["source_title"], "GSE193066 model source title drift")
            require(row["biopsy"] == "1st biopsy", "GSE193066 second biopsy entered model")
            require(row["donor_id_source"] == "deposited_title_root", "GSE193066 donor provenance drift")
            observed_gse_samples.add(row["sample_id"])
        else:
            require(
                row["donor_id"] == f"{row['dataset']}::{row['sample_id']}",
                "single-biopsy cohort donor ID drift",
            )
            require(row["donor_id_source"] == "single_biopsy_sample_id", "donor provenance drift")
        samples[transition][row["sample_id"]] = row
        donors[transition].add(row["donor_id"])
        all_primary_donors.add(row["donor_id"])
    for transition, (contrast, low_stage, high_stage, n_low, n_high, n_samples, n_cohorts) in expected_by_transition.items():
        rows = list(samples[transition].values())
        require(len(rows) == n_samples, f"{transition} sample count drift")
        require(len(donors[transition]) == n_samples, f"{transition} donor count drift")
        require({row["dataset"] for row in rows} == EXPECTED_COHORTS[transition], f"{transition} cohort set drift")
        require(sum(row["arm"] == "low" for row in rows) == n_low, f"{transition} low-arm drift")
        require(sum(row["arm"] == "high" for row in rows) == n_high, f"{transition} high-arm drift")
        for row in rows:
            stage = parse_int(row["fibrosis_stage"], f"{transition} fibrosis_stage")
            require((row["arm"] == "low" and stage == low_stage) or (row["arm"] == "high" and stage == high_stage), f"{transition} stage/arm mismatch")
    require(observed_gse_samples == included_gse_samples, "eligible GSE193066 first-biopsy universe is incomplete")
    require(len(all_primary_donors) == 664, "primary exact-stage donor universe drift")
    passed("V07_SAMPLE_CENSUS", f"{len(sample_rows)} donor-by-contrast rows match the transition-specific census")

    cohort_rows = read_rows(root / "audits/cohort_arm_audit.tsv")
    require(len(cohort_rows) == 23, "cohort-arm audit must contain 6+6+6+5 rows")
    observed_cohort_keys = set()
    for row in cohort_rows:
        key = (row["transition"], row["dataset"])
        require(key not in observed_cohort_keys, "duplicate cohort-arm row")
        observed_cohort_keys.add(key)
        require(row["dataset"] in ALLOWLIST, "cohort-arm row outside allowlist")
        require(row["both_arms_present"] == "TRUE", "single-arm cohort entered model")
        low = parse_int(row["n_low"], "cohort n_low")
        high = parse_int(row["n_high"], "cohort n_high")
        total = parse_int(row["n_total"], "cohort n_total")
        require(low > 0 and high > 0 and low + high == total, "cohort-arm counts are invalid")
    require(
        observed_cohort_keys
        == {
            (transition, cohort)
            for transition, cohorts in EXPECTED_COHORTS.items()
            for cohort in cohorts
        },
        "cohort-arm transition membership drift",
    )
    passed("V08_COHORT_ARM_GATE", "every included transition/cohort contains both arms")

    design_rows = read_rows(root / "audits/design_audit.tsv")
    require(len(design_rows) == 4, "design audit must contain four rows")
    for row in design_rows:
        transition = row["transition"]
        require(transition in expected_by_transition, "unknown design transition")
        require(row["formula"] == MODEL and row["coefficient"] == "fib_grouphigh", "model formula drift")
        require(row["dropped_terms"] == "", "a model term was dropped")
        require(parse_int(row["design_rank"], "design rank") == parse_int(row["n_design_columns"], "design columns"), "rank-deficient design")
        require(parse_int(row["n_genes_tested"], "design tested genes") == EXPECTED_N_GENES, "tested family drift")
        expected_cohorts = EXPECTED_COHORTS[transition]
        require(parse_int(row["n_cohorts"], "design cohorts") == len(expected_cohorts), "design cohort count drift")
        require(set(row["cohort_set"].split(";")) == expected_cohorts, "design cohort membership drift")
        require(set(row["expected_cohort_set"].split(";")) == expected_cohorts, "expected design cohort set drift")
        require(parse_int(row["n_unique_donors"], "design donors") == len(samples[transition]), "design donor count drift")
    passed("V09_DESIGN", "all four full-rank fixed-effect designs retain dataset, sex, and target coefficient")

    model_rows = read_rows(root / "audits/model_design.tsv")
    model_seen: dict[str, set[str]] = defaultdict(set)
    for row in model_rows:
        transition = row["transition"]
        sample_id = row["sample_id"]
        require(sample_id in samples[transition], "model-design sample absent from manifest")
        require(sample_id not in model_seen[transition], "duplicate model-design sample")
        model_seen[transition].add(sample_id)
        indicator = parse_float(row["fib_grouphigh"], "fib_grouphigh")
        require(indicator in {0.0, 1.0}, "fib_grouphigh is not binary")
        require((row["arm"] == "high") == (indicator == 1.0), "coefficient coding is reversed")
    require(all(model_seen[t] == set(samples[t]) for t in expected_by_transition), "model-design census drift")
    passed("V10_COEFFICIENT_CODING", "fib_grouphigh=1 exactly encodes the higher stage")

    summary_rows = {row["contrast"]: row for row in read_rows(root / "results/transition_summary.tsv")}
    require(set(summary_rows) == set(EXPECTED), "transition summary contrast set drift")
    result_path = root / "results/fibrosis_consecutive_lvqw_true_kleiner.csv"
    result_counts: dict[str, dict[str, int]] = defaultdict(lambda: {"tested": 0, "sig": 0, "up": 0, "down": 0})
    result_genes: dict[str, set[str]] = defaultdict(set)
    p_values: dict[str, list[float]] = defaultdict(list)
    padj_values: dict[str, list[float]] = defaultdict(list)
    with result_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required_fields = {"gene", "logFC", "SE", "t", "P.Value", "padj", "contrast", "transition", "method", "n_cohorts", "n_samples", "n_low", "n_high"}
        require(required_fields.issubset(reader.fieldnames or []), "main result schema is incomplete")
        for row in reader:
            contrast = row["contrast"]
            require(contrast in EXPECTED, f"unknown result contrast: {contrast}")
            transition, low_stage, high_stage, n_low, n_high, n_samples, n_cohorts = EXPECTED[contrast]
            require(row["transition"] == transition and row["method"] == METHOD, "result provenance drift")
            require(parse_int(row["n_samples"], "result n_samples") == n_samples, "result sample count drift")
            require(parse_int(row["n_cohorts"], "result n_cohorts") == n_cohorts, "result cohort count drift")
            require(parse_int(row["n_low"], "result n_low") == n_low and parse_int(row["n_high"], "result n_high") == n_high, "result arm count drift")
            gene = row["gene"]
            require(gene and gene not in result_genes[contrast], "blank or duplicate gene within contrast")
            result_genes[contrast].add(gene)
            logfc = parse_float(row["logFC"], f"{contrast} logFC")
            se = parse_float(row["SE"], f"{contrast} SE")
            t_value = parse_float(row["t"], f"{contrast} t")
            p_value = parse_float(row["P.Value"], f"{contrast} P.Value")
            padj = parse_float(row["padj"], f"{contrast} padj")
            require(se > 0 and 0 <= p_value <= 1 and 0 <= padj <= 1, "invalid coefficient statistics")
            require(abs(logfc - t_value * se) <= 5e-10 * (1 + abs(logfc)), "logFC/SE/t identity failed")
            p_values[contrast].append(p_value)
            padj_values[contrast].append(padj)
            count = result_counts[contrast]
            count["tested"] += 1
            if padj < 0.05:
                count["sig"] += 1
                require(logfc != 0, "significant coefficient has exactly zero logFC")
                count["up" if logfc > 0 else "down"] += 1
    require(len(result_genes) == 4, "main result lacks a contrast")
    first_genes = result_genes[next(iter(EXPECTED))]
    for contrast in EXPECTED:
        require(result_counts[contrast]["tested"] == EXPECTED_N_GENES, f"{contrast} tested family is incomplete")
        require(result_genes[contrast] == first_genes, f"{contrast} gene universe differs")
        rederived = bh_adjust(p_values[contrast])
        max_delta = max(abs(left - right) for left, right in zip(rederived, padj_values[contrast]))
        require(max_delta <= 5e-12, f"{contrast} BH rederivation failed: {max_delta}")
        summary = summary_rows[contrast]
        counts = result_counts[contrast]
        require(parse_int(summary["n_genes_tested"], "summary tested") == counts["tested"], "summary tested count drift")
        require(parse_int(summary["n_deg_05"], "summary significant") == counts["sig"], "summary DEG count drift")
        require(parse_int(summary["n_up"], "summary up") == counts["up"], "summary up count drift")
        require(parse_int(summary["n_down"], "summary down") == counts["down"], "summary down count drift")
        require(summary["significance_rule"] == "BH padj < 0.05" and summary["lfc_gate"] == "none", "post-hoc LFC gate entered summary")
    passed("V11_RESULTS_AND_BH", "110,552 coefficients and four complete BH families independently reproduce")

    shuffle_rows = read_rows(root / "audits/shuffle_audit.tsv")
    require(len(shuffle_rows) == 23, "shuffle audit must contain 6+6+6+5 rows")
    observed_shuffle_keys: set[tuple[str, str]] = set()
    changed_by_transition: dict[str, int] = defaultdict(int)
    for row in shuffle_rows:
        require(row["transition"] in EXPECTED_COHORTS, "unknown shuffle transition")
        key = (row["transition"], row["dataset"])
        require(key not in observed_shuffle_keys, "duplicate shuffle transition/cohort")
        observed_shuffle_keys.add(key)
        require(row["dataset"] in EXPECTED_COHORTS[row["transition"]], "shuffle cohort outside transition universe")
        require(parse_int(row["original_low"], "shuffle original low") == parse_int(row["shuffled_low"], "shuffle low"), "shuffle changed low count")
        require(parse_int(row["original_high"], "shuffle original high") == parse_int(row["shuffled_high"], "shuffle high"), "shuffle changed high count")
        changed_by_transition[row["transition"]] += parse_int(row["n_changed"], "shuffle changed")
    require(all(changed_by_transition[t] > 0 for t in expected_by_transition), "shuffle changed no labels in a transition")
    require(
        observed_shuffle_keys == observed_cohort_keys,
        "shuffle transition/cohort universe differs from the fitted universe",
    )

    null_summary = {row["contrast"]: row for row in read_rows(root / "audits/null_shuffle_summary.tsv")}
    require(set(null_summary) == set(EXPECTED), "null summary contrast set drift")
    null_counts: dict[str, int] = defaultdict(int)
    null_sig: dict[str, int] = defaultdict(int)
    null_p: dict[str, list[float]] = defaultdict(list)
    null_q: dict[str, list[float]] = defaultdict(list)
    with (root / "audits/null_shuffle_results.csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            contrast = row["contrast"]
            require(contrast in EXPECTED and row["transition"] == EXPECTED[contrast][0], "null contrast mapping drift")
            p_value = parse_float(row["P.Value"], "null P.Value")
            padj = parse_float(row["padj"], "null padj")
            require(0 <= p_value <= 1 and 0 <= padj <= 1, "invalid null probability")
            null_p[contrast].append(p_value)
            null_q[contrast].append(padj)
            null_counts[contrast] += 1
            null_sig[contrast] += padj < 0.05
    for contrast in EXPECTED:
        require(null_counts[contrast] == EXPECTED_N_GENES, f"{contrast} null family is incomplete")
        max_delta = max(abs(left - right) for left, right in zip(bh_adjust(null_p[contrast]), null_q[contrast]))
        require(max_delta <= 5e-12, f"{contrast} null BH rederivation failed")
        require(parse_int(null_summary[contrast]["n_deg_05"], "null summary n_deg") == null_sig[contrast], "null summary count drift")
        require(null_summary[contrast]["interpretation"].startswith("diagnostic shuffle only"), "null interpretation drift")
    passed("V12_SHUFFLED_NULL", "within-cohort donor shuffles preserve arm counts and four null BH families reproduce")

    status_rows = read_rows(root / "manifests/producer_status.tsv")
    require(len(status_rows) == 1 and status_rows[0]["status"] == "PRODUCER_COMPLETE_PENDING_VALIDATION", "producer status drift")
    require(status_rows[0]["canonical_write"] == "FALSE", "producer claims canonical write")
    require(parse_int(status_rows[0]["n_unique_primary_donors"], "producer donor count") == 664, "producer donor count drift")
    passed("V13_NONCANONICAL", "producer explicitly records candidate-only, noncanonical status")

    return {
        "checks": checks,
        "job_id": environment[0]["slurm_job_id"],
        "summary": {
            contrast: {
                "n_genes_tested": result_counts[contrast]["tested"],
                "n_deg_05": result_counts[contrast]["sig"],
                "n_up": result_counts[contrast]["up"],
                "n_down": result_counts[contrast]["down"],
                "n_samples": EXPECTED[contrast][5],
                "n_cohorts": EXPECTED[contrast][6],
                "null_n_deg_05": null_sig[contrast],
            }
            for contrast in EXPECTED
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    arguments = parser.parse_args()

    result = validate(arguments.project_root, arguments.output_root)
    root = arguments.output_root.resolve(strict=True)
    manifest_dir = root / "manifests"
    checks_path = manifest_dir / "validation_checks.tsv"
    report_path = manifest_dir / "validation_report.json"
    validated_manifest_path = manifest_dir / "validated_artifact_manifest.tsv"
    ready_path = root / "READY"

    checks_payload = tsv_payload(result["checks"], ["check_id", "status", "detail"])
    write_exclusive(checks_path, checks_payload)
    report = {
        "status": "READY",
        "candidate_id": CANDIDATE_ID,
        "analysis_id": ANALYSIS_ID,
        "canonical_promotion_status": "not_promoted",
        "job_id": result["job_id"],
        "n_checks": len(result["checks"]),
        "summary": result["summary"],
        "validated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_exclusive(report_path, json.dumps(report, indent=2, sort_keys=True) + "\n")

    manifest_rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file() and p not in {validated_manifest_path, ready_path}):
        manifest_rows.append(
            {
                "artifact_path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    manifest_payload = tsv_payload(manifest_rows, ["artifact_path", "bytes", "sha256"])
    write_exclusive(validated_manifest_path, manifest_payload)
    ready = {
        "status": "READY",
        "candidate_id": CANDIDATE_ID,
        "analysis_id": ANALYSIS_ID,
        "canonical_promotion_status": "not_promoted",
        "validation_report_sha256": sha256(report_path),
        "validation_checks_sha256": sha256(checks_path),
        "validated_artifact_manifest_sha256": sha256(validated_manifest_path),
    }
    write_exclusive(ready_path, json.dumps(ready, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

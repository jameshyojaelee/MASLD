#!/usr/bin/env python3
"""Verify the immutable BG-001 analysis/scientific contract and baseline seal."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
from pathlib import Path


RUN_RE = re.compile(r"^bg001-fragment-v211-gencode49-[0-9]{8}T[0-9]{6}Z$")
EXPECTED_THRESHOLDS = {
    "qc_symmetric_difference_max": 4,
    "gene_jaccard_min": 0.99,
    "sample_weight_spearman_min": 0.99,
    "offset_median_abs_delta_max": 0.01,
    "offset_p95_abs_delta_max": 0.03,
    "offset_max_abs_delta_max": 0.10,
    "offset_group_shift_max": 0.02,
    "logfc_spearman_min": 0.995,
    "expressed_logfc_median_abs_delta_max": 0.01,
    "expressed_logfc_p95_abs_delta_max": 0.05,
    "treat_jaccard_min": 0.95,
    "treat_count_relative_change_max": 0.05,
}
EXPECTED_ARMS = ("R0", "F_locked", "F_legacy", "F_five")
EXPECTED_ARM_MODES = {
    "R0": {"counts": "read", "qc": "recompute", "filter_scope": "legacy_all", "gene_mode": "native"},
    "F_locked": {"counts": "fragment", "qc": "locked", "filter_scope": "legacy_all", "gene_mode": "locked"},
    "F_legacy": {"counts": "fragment", "qc": "recompute", "filter_scope": "legacy_all", "gene_mode": "native"},
    "F_five": {"counts": "fragment", "qc": "reuse_F_legacy", "filter_scope": "canonical_five", "gene_mode": "native"},
}
EXPECTED_ANALYSIS_MODEL = {
    "formula": "~ dataset + inferred_sex + group_binary",
    "engine": "voomWithQualityWeights -> lmFit -> treat",
    "disease_coefficient": "group_binaryDisease",
    "treat_lfc": 0.25,
    "qc_seed": 42,
}
EXPECTED_ACTIVE_COHORTS = (
    "GSE130970", "GSE135251", "GSE174478", "GSE213621", "GSE240729",
)
EXPECTED_CANONICAL_COHORTS = ("GSE130970", "GSE135251", "GSE213621")
EXPECTED_FEATURECOUNTS = {
    "version": "2.1.1",
    "arguments": ["-p", "--countReadPairs", "-B", "-s", "2"],
    "gtf": "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz",
    "gtf_sha256": "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9",
}
BAM_GUARD_INPUTS = (
    "manifests/bam_manifest_finalize.preflight.json",
    "manifests/bam_manifest_finalize.postflight.json",
)
EXPECTED_BASELINE_INPUTS = (
    "manifests/read_count_sources.tsv",
    "manifests/read_count_overrides.tsv",
    "manifests/protected_source_manifest.tsv",
    "contract/run_contract.json",
    "contract/source_manifest.tsv",
    "contract/analysis_runtime_contract.json",
    "contract/analysis_environment.explicit.txt",
    "contract/analysis_environment.conda.json",
    "contract/analysis_r_dependency_files.tsv",
    "contract/analysis_native_dependencies.tsv",
    *BAM_GUARD_INPUTS,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_analysis_contract(contract: dict, root: Path) -> None:
    run_id = root.name
    expected_arms = {arm: str(root / "arms" / arm) for arm in EXPECTED_ARMS}
    failures: list[str] = []
    if contract.get("schema_version") != "1.2":
        failures.append("schema_version")
    if contract.get("run_id") != run_id:
        failures.append("run_id")
    if contract.get("candidate_root") != str(root):
        failures.append("candidate_root")
    if contract.get("reference_root") != str(root / "source_snapshot"):
        failures.append("reference_root")
    if contract.get("thresholds") != EXPECTED_THRESHOLDS:
        failures.append("thresholds")
    if contract.get("arms") != expected_arms:
        failures.append("arms")
    if contract.get("arm_modes") != EXPECTED_ARM_MODES:
        failures.append("arm_modes")
    if contract.get("analysis_model") != EXPECTED_ANALYSIS_MODEL:
        failures.append("analysis_model")
    if contract.get("active_affected_cohorts") != list(EXPECTED_ACTIVE_COHORTS):
        failures.append("active_affected_cohorts")
    if contract.get("canonical_affected_cohorts") != list(EXPECTED_CANONICAL_COHORTS):
        failures.append("canonical_affected_cohorts")
    if contract.get("featurecounts") != EXPECTED_FEATURECOUNTS:
        failures.append("featurecounts")
    if contract.get("bam_manifest_finalization_required") is not True:
        failures.append("bam_manifest_finalization_required")
    if failures:
        raise SystemExit("Analysis contract differs from frozen specification: " + ",".join(failures))


def parse_sha_manifest(path: Path, root: Path) -> dict[str, str]:
    if b"\r" in path.read_bytes():
        raise SystemExit(f"SHA manifest is not LF-only: {path}")
    rows: dict[str, str] = {}
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        parts = line.split(maxsplit=1)
        if len(parts) != 2 or not re.fullmatch(r"[0-9a-f]{64}", parts[0]):
            raise SystemExit(f"Malformed SHA manifest row {path}:{number}")
        candidate = Path(parts[1])
        if not candidate.is_absolute():
            candidate = root / candidate
        try:
            resolved_parent = candidate.parent.resolve(strict=True)
            resolved_parent.relative_to(root)
        except (FileNotFoundError, ValueError) as exc:
            raise SystemExit(f"SHA manifest path escapes run root: {candidate}") from exc
        relative = str(candidate.relative_to(root))
        if relative in rows:
            raise SystemExit(f"Duplicate SHA manifest path: {relative}")
        try:
            mode = os.lstat(candidate).st_mode
        except FileNotFoundError as exc:
            raise SystemExit(f"SHA manifest target is missing: {candidate}") from exc
        if stat.S_ISLNK(mode) or not stat.S_ISREG(mode):
            raise SystemExit(f"SHA manifest target is symlinked or non-regular: {candidate}")
        rows[relative] = parts[0]
    return rows


def validate_baseline_seal(root: Path, contract: dict) -> None:
    # Reconstruct the guarded BAM-manifest transaction before trusting hashes
    # embedded in the baseline seal. This proves semantic validity as well as
    # exact byte binding of both pre/postflight JSON records.
    from bam_manifest_finalization_guard import GuardError, verify_artifacts

    try:
        verify_artifacts(root)
    except GuardError as exc:
        raise SystemExit(f"BAM-manifest guard verification failed: {exc}") from exc
    files_manifest = root / "manifests/baseline_frozen_files.sha256"
    inputs_manifest = root / "manifests/baseline_manifests.sha256"
    marker_path = root / "BASELINE_FROZEN.json"
    compatibility_path = root / "BASELINE_FROZEN"
    for path in (files_manifest, inputs_manifest, marker_path, compatibility_path):
        try:
            mode = os.lstat(path).st_mode
        except FileNotFoundError as exc:
            raise SystemExit(f"Baseline transaction is incomplete: {path}") from exc
        if stat.S_ISLNK(mode) or not stat.S_ISREG(mode):
            raise SystemExit(f"Baseline transaction path is symlinked or non-regular: {path}")

    frozen_rows = parse_sha_manifest(files_manifest, root)
    observed_frozen: set[str] = set()
    for path in (root / "frozen_sets").rglob("*"):
        mode = os.lstat(path).st_mode
        if stat.S_ISLNK(mode) or (not stat.S_ISDIR(mode) and not stat.S_ISREG(mode)):
            raise SystemExit(f"Frozen baseline contains symlink/special path: {path}")
        if stat.S_ISREG(mode):
            observed_frozen.add(str(path.relative_to(root)))
    if set(frozen_rows) != observed_frozen:
        raise SystemExit("Baseline frozen-file manifest does not cover the exact file set")
    for relative, expected in frozen_rows.items():
        if sha256(root / relative) != expected:
            raise SystemExit(f"Frozen baseline hash drift: {relative}")

    input_rows = parse_sha_manifest(inputs_manifest, root)
    if set(input_rows) != set(EXPECTED_BASELINE_INPUTS):
        raise SystemExit("Baseline input manifest does not bind the exact required inputs")
    for relative, expected in input_rows.items():
        if sha256(root / relative) != expected:
            raise SystemExit(f"Baseline input hash drift: {relative}")

    marker = json.loads(marker_path.read_text())
    expected_marker = {
        "schema": "bg001-baseline-seal-v1",
        "run_id": root.name,
        "baseline_frozen_files_manifest": str(files_manifest.relative_to(root)),
        "baseline_frozen_files_manifest_sha256": sha256(files_manifest),
        "baseline_manifests_manifest": str(inputs_manifest.relative_to(root)),
        "baseline_manifests_manifest_sha256": sha256(inputs_manifest),
        "frozen_file_count": len(frozen_rows),
        "bound_input_count": len(input_rows),
        "run_contract_sha256": sha256(root / "contract/run_contract.json"),
        "source_manifest_sha256": sha256(root / "contract/source_manifest.tsv"),
        "compatibility_marker": str(compatibility_path.relative_to(root)),
        "compatibility_marker_contract": "lowercase_sha256(BASELINE_FROZEN.json)+LF",
    }
    if marker != expected_marker:
        raise SystemExit("BASELINE_FROZEN transaction marker differs from its manifests")
    if marker["source_manifest_sha256"] != contract["source_manifest_sha256"]:
        raise SystemExit("Baseline marker/source contract hash mismatch")
    expected_compatibility = sha256(marker_path) + "\n"
    try:
        observed_compatibility = compatibility_path.read_text(encoding="ascii")
    except UnicodeError as exc:
        raise SystemExit("BASELINE_FROZEN compatibility marker is not ASCII") from exc
    if observed_compatibility != expected_compatibility:
        raise SystemExit(
            "BASELINE_FROZEN compatibility marker does not bind BASELINE_FROZEN.json"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--require-baseline-frozen", action="store_true")
    args = parser.parse_args()
    root = args.run_root.resolve(strict=True)
    if not RUN_RE.fullmatch(root.name):
        raise SystemExit("Invalid analysis run root name")
    contract_path = root / "contract/run_contract.json"
    contract = json.loads(contract_path.read_text())
    validate_analysis_contract(contract, root)
    if sha256(root / "contract/source_manifest.tsv") != contract.get("source_manifest_sha256"):
        raise SystemExit("Analysis contract source-manifest hash drift")
    if args.require_baseline_frozen or (root / "BASELINE_FROZEN.json").exists():
        validate_baseline_seal(root, contract)
    print("PASS immutable BG-001 analysis contract", root.name)


if __name__ == "__main__":
    main()

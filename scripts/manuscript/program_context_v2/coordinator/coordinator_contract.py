#!/usr/bin/env python3
"""Fail-closed contracts and presentation adapters for Plan 60.

This module is deliberately upstream-read-only.  Workstream owners publish
their own signed ``plan60_terminal_artifacts.tsv`` handoffs.  The coordinator
validates those handoffs, derives presentation-only rows, and prepares inputs
for REL00--05 without changing a scientific call.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import importlib
import json
import math
import platform
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Iterable, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from release_common import (  # noqa: E402, F401 -- intentional coordinator facade
    ARTIFACT_FIELDS,
    CANDIDATE_ID,
    CLOSURE_FIELDS,
    PLAN50_FROZEN_PRODUCER_PATHS,
    PLAN50_MANIFEST_PRODUCER_PATHS,
    PROTECTED_BASELINE_FIELDS,
    PROTECTED_SCOPE_FIELDS,
    ReleaseContractError,
    allowed_workstream_roots,
    atomic_write_json,
    atomic_write_tsv,
    parse_nonnegative_int,
    project_relative,
    read_tsv_exact,
    require_sha256,
    require_within_any,
    sha256_file,
    validate_artifact_manifest,
)
from release_products import (  # noqa: E402, F401 -- intentional coordinator facade
    BASE_SELECTION_FIELDS,
    EXPECTED_FIGURE_TITLES,
    SOURCE_ROW_FIELDS,
    validate_source_dependency_contract,
)
from fibrosis_candidate_contract import (  # noqa: E402
    FIBROSIS_PRIMARY_ARTIFACT_ID,
    FIBROSIS_PRIMARY_RELATIVE,
    fibrosis_snapshot_relpath,
)
from publish_plan60_workstream_handoff import (  # noqa: E402
    PLAN50_PAYLOAD_MANIFEST_FIELDS,
    workstream_contract,
)


class CoordinatorContractError(ReleaseContractError):
    """Raised before any coordinator preparation output is written."""


SOURCE_EVIDENCE_FIELDS = SOURCE_ROW_FIELDS[:-2]
HANDOFF_SIGNATURE_VERSION = "plan60_workstream_handoff_v1"
COORDINATOR_SIGNATURE_VERSION = "plan60_coordinator_attestation_v1"
PREPARATION_CONTRACT_VERSION = 1
FIXED_HERO_GENES = ("THRB", "HKDC1", "GLP1R", "MTARC1")
PLAN50_ANALYSIS_RELEASE_ID = f"{CANDIDATE_ID}-passports-v1"
EXPECTED_COHORT_COUNTS = {
    "qc_samples": 1260,
    "qc_cohorts": 9,
    "pooled_samples": 846,
    "pooled_cohorts": 5,
}
EXPECTED_COMPOSITION_COUNTS = {
    "n_pass_technical": 1260,
    "n_deconvolved_qc_pass": 1221,
    "n_assay_missing": 39,
    "n_qc_failures_excluded": 21,
    "n_declared_celltypes": 22,
    "n_testable_celltypes": 16,
    "n_untestable_celltypes": 6,
}
COMPOSITION_MODEL = "asin(sqrt(proportion)) ~ group_binary + dataset + inferred_sex"
COMPOSITION_MULTIPLICITY_FAMILY = (
    "BH across all 16 testable deposited cell-type columns"
)
COMPOSITION_RESULT_FIELDS = (
    "celltype",
    "testability",
    "n_total_qc",
    "n_analyzed",
    "n_control",
    "n_disease",
    "effect_arcsin_sqrt",
    "t",
    "pvalue",
    "padj",
    "model_formula",
    "biological_unit",
    "multiplicity_family",
    "source_assay",
    "allowed_wording",
    "prohibited_wording",
)
COMPOSITION_TESTABILITY_FIELDS = (
    "celltype",
    "n_nonmissing_qc",
    "testability",
    "reason",
)
COMPOSITION_READY_FIELDS = (
    "candidate_id",
    "status",
    "n_pass_technical",
    "n_deconvolved_qc_pass",
    "n_assay_missing",
    "n_qc_failures_excluded",
    "n_declared_celltypes",
    "n_testable_celltypes",
    "n_untestable_celltypes",
    "model",
    "biological_unit",
    "multiplicity_family",
    "results_sha256",
    "testability_sha256",
    "audit_sha256",
    "validation_sha256",
    "validator_sha256",
    "canonical_promotion_authorized",
)
CORE_RUNTIME_FIELDS = (
    "environment_id",
    "record_type",
    "name",
    "version",
    "executable_path",
    "executable_sha256",
    "source_scope",
)
CORE_PACKAGE_FIELDS = (
    "environment_id",
    "distribution_id",
    "package",
    "version",
    "installation_location",
    "metadata_path",
    "resolution_status",
)
CORE_SYSPATH_FIELDS = (
    "environment_id",
    "path_order",
    "path",
    "exists",
    "is_directory",
)
CORE_IMPORT_FIELDS = (
    "environment_id",
    "module",
    "imported_version",
    "imported_path",
    "imported_file_sha256",
    "metadata_version",
)
CORE_PRODUCER_FIELDS = (
    "producer_id",
    "repository_path",
    "sha256",
    "bytes",
)
NMF_CONTINUOUS_FIELDS = (
    "sample_id",
    "k",
    "program_code",
    "program_label",
    "continuous_loading",
    "source_sha256",
)
NMF_SOURCE_MANIFEST_FIELDS = (
    "source_id",
    "source_role",
    "origin_project_path",
    "origin_sha256",
    "origin_bytes",
    "bundle_relative_path",
    "bundle_sha256",
    "bundle_bytes",
)
NMF_VALIDATION_FIELDS = (
    "check_id",
    "status",
    "observed",
    "criterion",
)
NMF_READY_FIELDS = (
    "candidate_id",
    "status",
    "source_manifest_sha256",
    "loadings_sha256",
    "evidence_sha256",
    "validation_sha256",
    "producer_manifest_sha256",
    "n_samples",
    "n_loading_rows",
    "n_evidence_rows",
    "n_s2a_marks",
    "n_s2b_marks",
    "canonical_promotion_authorized",
    "frozen_on_date",
)
EXPECTED_PLAN13_TERMINAL_SHA256 = (
    "f97a3a349e5c2ffdb484e06679a2ed12a061c6f5f121f6f9a144787b004cbc8e"
)
FORBIDDEN_PRESENTATION_TOKENS = (
    "universal_score",
    "combined_score",
    "overall_score",
    "overall_rank",
    "cross_assay_rank",
    "modality_count",
    "evidence_count",
    "priority_rank",
)
WORKSTREAM_ORDER = ("PLAN13", "PLAN20", "PLAN30", "PLAN40", "PLAN50")
WORKSTREAM_LANES = {
    "PLAN13": "SP-INT",
    "PLAN20": "HS-V2",
    "PLAN30": "GEN",
    "PLAN40": "HLF",
    "PLAN50": "PASS",
}


@dataclass(frozen=True)
class RealPaths:
    project_root: Path
    multimodal_root: Path
    plan13_root: Path
    plan20_root: Path
    plan30_root: Path
    plan40_root: Path
    plan50_root: Path

    @classmethod
    def build(cls, project_root: Path) -> "RealPaths":
        project = project_root.resolve()
        candidate = (
            project / "Analysis/Multimodal_Program_Projection/candidates" / CANDIDATE_ID
        )
        return cls(
            project_root=project,
            multimodal_root=candidate,
            plan13_root=candidate / "spatial_context_semantic_v2_2026-08-08",
            plan20_root=candidate / "hotspot",
            plan30_root=candidate / "genetics_context",
            plan40_root=candidate / "myojin_hlf",
            plan50_root=(
                project / "RNA-seq/results/evidence_passports/candidates" / CANDIDATE_ID
            ),
        )

    def workstream_root(self, workstream_id: str) -> Path:
        return {
            "PLAN13": self.plan13_root,
            "PLAN20": self.plan20_root,
            "PLAN30": self.plan30_root,
            "PLAN40": self.plan40_root,
            "PLAN50": self.plan50_root,
        }[workstream_id]

    def owner_manifest(self, workstream_id: str) -> Path:
        return self.workstream_root(workstream_id) / "plan60_terminal_artifacts.tsv"

    def owner_signature(self, workstream_id: str) -> Path:
        return (
            self.workstream_root(workstream_id)
            / "plan60_terminal_artifacts.signature.json"
        )


@dataclass(frozen=True)
class Handoff:
    workstream_id: str
    root: Path
    manifest: Path
    signature: Path
    manifest_sha256: str
    artifacts: tuple[object, ...]


def one_row_tsv(path: Path, fields: Sequence[str]) -> dict[str, str]:
    rows = read_tsv_exact(path, fields)
    if len(rows) != 1:
        raise CoordinatorContractError(
            f"expected exactly one row in {path}; found {len(rows)}"
        )
    return rows[0]


def read_tsv_flexible(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    if not path.is_file() or path.is_symlink():
        raise CoordinatorContractError(
            f"required regular table is missing or symlinked: {path}"
        )
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = tuple(reader.fieldnames or ())
        if not fields:
            raise CoordinatorContractError(f"table has no header: {path}")
        return fields, [dict(row) for row in reader]


def parse_bool(value: object, context: str) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value)
    if text in {"TRUE", "true"}:
        return True
    if text in {"FALSE", "false"}:
        return False
    raise CoordinatorContractError(f"{context} is not an exact Boolean: {value!r}")


def parse_float(value: object, context: str, allow_blank: bool = False) -> float | None:
    text = str(value).strip()
    if allow_blank and not text:
        return None
    try:
        result = float(text)
    except ValueError as error:
        raise CoordinatorContractError(
            f"{context} is not numeric: {value!r}"
        ) from error
    if not math.isfinite(result):
        raise CoordinatorContractError(f"{context} is non-finite: {value!r}")
    return result


def bh_adjust(pvalues: Sequence[float]) -> list[float]:
    """Re-derive Benjamini-Hochberg values without a package dependency."""

    if not pvalues or any(
        not math.isfinite(value) or value < 0 or value > 1 for value in pvalues
    ):
        raise CoordinatorContractError(
            "BH family contains no values or an invalid p value"
        )
    order = sorted(range(len(pvalues)), key=lambda index: pvalues[index])
    adjusted = [1.0] * len(pvalues)
    running = 1.0
    total = len(pvalues)
    for rank_index in range(total - 1, -1, -1):
        index = order[rank_index]
        running = min(running, pvalues[index] * total / (rank_index + 1), 1.0)
        adjusted[index] = running
    return adjusted


def require_numeric_equal(left: object, right: object, context: str) -> None:
    observed = parse_float(left, f"{context} observed")
    expected = parse_float(right, f"{context} expected")
    assert observed is not None and expected is not None
    if not math.isclose(observed, expected, rel_tol=1e-12, abs_tol=1e-12):
        raise CoordinatorContractError(
            f"{context} drift: observed={observed!r}, expected={expected!r}"
        )


def validate_hotspot_registry_projection(
    registry_rows: Sequence[Mapping[str, str]],
    figure_rows: Sequence[Mapping[str, str]],
) -> None:
    """Prove the complete Figure-2 source is a lossless registry projection."""

    registry = {row["program_uid"]: row for row in registry_rows}
    figure = {row["program_uid"]: row for row in figure_rows}
    if (
        len(registry_rows) != 117
        or len(figure_rows) != 117
        or len(registry) != 117
        or len(figure) != 117
        or set(registry) != set(figure)
    ):
        raise CoordinatorContractError(
            "Hotspot registry/Figure2 projection is not one common 117-program family"
        )
    exact_map = {
        "release_id": "release_id",
        "membership_sha256": "membership_sha256",
        "cell_type": "cell_type",
        "module": "module",
        "module_name": "module_name",
        "module_hallmark": "module_hallmark",
        "module_top_pathway": "module_top_pathway",
        "direction": "primary_direction",
        "primary_selected": "primary_selected",
        "hc3_supported": "hc3_supported",
        "scoring_direction_agree": "scoring_direction_agree",
        "selected_hc3_fragile": "selected_hc3_fragile",
        "robust_display": "robust_display",
        "source_stability": "source_stability",
        "source_stability_reason": "source_stability_reason",
        "n_donors": "primary_n_donors",
        "n_datasets": "primary_n_datasets",
        "n_healthy": "primary_n_stage0",
        "n_steatosis": "primary_n_stage1",
        "n_steatohepatitis": "primary_n_stage2",
    }
    numeric_map = {
        "beta": "primary_beta",
        "se": "primary_se",
        "statistic": "primary_statistic",
        "pvalue": "primary_pvalue",
        "qvalue": "primary_qvalue",
        "hc3_se": "primary_hc3_se",
        "hc3_statistic": "primary_hc3_statistic",
        "hc3_pvalue": "primary_hc3_pvalue",
        "hc3_qvalue": "primary_hc3_qvalue",
        "max_leverage": "primary_max_leverage",
        "stability_score": "stability_score",
        "stability_mean": "stability_mean",
        "stability_median": "stability_median",
    }
    for uid, observed in figure.items():
        source = registry[uid]
        for figure_field, registry_field in exact_map.items():
            if observed[figure_field] != source[registry_field]:
                raise CoordinatorContractError(
                    f"Hotspot Figure2/registry mapped-field drift for "
                    f"{uid}:{figure_field}->{registry_field}"
                )
        for figure_field, registry_field in numeric_map.items():
            require_numeric_equal(
                observed[figure_field],
                source[registry_field],
                f"Hotspot {uid}:{figure_field}->{registry_field}",
            )


def validate_composition_contract(
    result_rows: Sequence[Mapping[str, str]],
    testability_rows: Sequence[Mapping[str, str]],
    ready: Mapping[str, str],
    audit_rows: Sequence[Mapping[str, str]] | None = None,
    validation_rows: Sequence[Mapping[str, str]] | None = None,
) -> None:
    """Validate the complete 22-column, QC-filtered sample composition release."""

    expected_ready = {
        "candidate_id": CANDIDATE_ID,
        "status": "QC_FILTERED_SAMPLE_LEVEL_COMPOSITION_READY",
        "model": COMPOSITION_MODEL,
        "biological_unit": "human liver sample",
        "multiplicity_family": COMPOSITION_MULTIPLICITY_FAMILY,
        "canonical_promotion_authorized": "false",
    }
    for key, value in expected_ready.items():
        if ready.get(key) != value:
            raise CoordinatorContractError(
                f"Plan20 composition READY {key} drift: {ready.get(key)!r}"
            )
    for key, expected in EXPECTED_COMPOSITION_COUNTS.items():
        if int(str(ready.get(key, "-1"))) != expected:
            raise CoordinatorContractError(
                f"Plan20 composition READY {key} drift: {ready.get(key)!r}"
            )

    by_result = {row["celltype"]: row for row in result_rows}
    by_testability = {row["celltype"]: row for row in testability_rows}
    if len(result_rows) != 16 or len(by_result) != 16:
        raise CoordinatorContractError(
            "Plan20 composition result is not the complete 16-test family"
        )
    if len(testability_rows) != 22 or len(by_testability) != 22:
        raise CoordinatorContractError(
            "Plan20 composition testability does not enumerate 22 columns"
        )
    testable = {
        celltype
        for celltype, row in by_testability.items()
        if row["testability"] == "testable"
    }
    untestable = {
        celltype
        for celltype, row in by_testability.items()
        if row["testability"] == "untestable_structurally_unavailable"
    }
    if (
        len(testable) != 16
        or len(untestable) != 6
        or testable | untestable != set(by_testability)
    ):
        raise CoordinatorContractError(
            "Plan20 composition 16-testable/6-untestable partition drift"
        )
    if set(by_result) != testable:
        raise CoordinatorContractError(
            "Plan20 composition results differ from its testable family"
        )
    for celltype in sorted(testable):
        row = by_result[celltype]
        checks = {
            "testability": "testable",
            "n_total_qc": "1260",
            "n_analyzed": "1221",
            "n_control": "152",
            "n_disease": "1069",
            "model_formula": COMPOSITION_MODEL,
            "biological_unit": "QC-passing human liver sample with MuSiC deconvolution",
            "multiplicity_family": COMPOSITION_MULTIPLICITY_FAMILY,
            "source_assay": "bulk RNA-seq MuSiC deconvolution",
            "allowed_wording": "QC-filtered sample-level cross-sectional composition association",
            "prohibited_wording": (
                "donor-level, longitudinal change, lineage transition, or causal composition effect"
            ),
        }
        for field, expected in checks.items():
            if row[field] != expected:
                raise CoordinatorContractError(
                    f"Plan20 composition {celltype} {field} drift: {row[field]!r}"
                )
        if int(by_testability[celltype]["n_nonmissing_qc"]) != 1221:
            raise CoordinatorContractError(
                f"Plan20 composition {celltype} nonmissing n drift"
            )
        parse_float(row["effect_arcsin_sqrt"], f"composition {celltype} effect")
        parse_float(row["t"], f"composition {celltype} t")
    for celltype in untestable:
        row = by_testability[celltype]
        if int(row["n_nonmissing_qc"]) != 0 or row["reason"] != (
            "deposited merged column is structurally all-NA"
        ):
            raise CoordinatorContractError(
                f"Plan20 composition untestable contract drift: {celltype}"
            )

    ordered = [by_result[celltype] for celltype in sorted(by_result)]
    adjusted = bh_adjust([float(str(row["pvalue"])) for row in ordered])
    for row, expected in zip(ordered, adjusted, strict=True):
        require_numeric_equal(
            row["padj"], expected, f"Plan20 composition BH {row['celltype']}"
        )

    if audit_rows is not None:
        audit = {row["audit_id"]: row for row in audit_rows}
        expected_audit = {
            "pass_technical_census": "1260",
            "qc_join": "1221",
            "deconvolution_missingness": "39",
            "excluded_qc_failures": "21",
            "declared_celltypes": "22",
            "testable_family": "16",
            "untestable_columns": "6",
            "bh_family": "16",
            "direction_sensitivity": "15",
            "fdr_sensitivity": "16",
            "biological_unit": "human liver sample",
            "canonical_write": "false",
        }
        if set(audit) != set(expected_audit):
            raise CoordinatorContractError("Plan20 composition audit universe drift")
        for audit_id, expected in expected_audit.items():
            if (
                audit[audit_id]["status"] != "pass"
                or audit[audit_id]["value"] != expected
            ):
                raise CoordinatorContractError(
                    f"Plan20 composition audit failed: {audit_id}"
                )
    if validation_rows is not None:
        required_checks = {
            "qc_census",
            "assay_join",
            "declared_universe",
            "ols_effects",
            "bh_rederivation",
            "sample_level_wording",
            "candidate_only",
        }
        if {row["check_id"] for row in validation_rows} != required_checks or any(
            row["status"] != "pass" for row in validation_rows
        ):
            raise CoordinatorContractError(
                "Plan20 composition validation report is incomplete"
            )


def canonical_json_sha256(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def core_runtime_environment_rows() -> list[dict[str, object]]:
    executable = Path(sys.executable).resolve()
    if not executable.is_file() or executable.is_symlink():
        raise CoordinatorContractError(
            f"core/coordinator Python executable is missing or unsafe: {executable}"
        )
    return [
        {
            "environment_id": "plan60_core_coordinator_python",
            "record_type": "runtime_executable",
            "name": "python",
            "version": platform.python_version(),
            "executable_path": str(executable),
            "executable_sha256": sha256_file(executable),
            "source_scope": "REL00-REL05 core plus coordinator adapters/validators",
        }
    ]


def core_package_export_rows() -> list[dict[str, object]]:
    distributions: dict[str, list[tuple[str, str, str]]] = {}
    for distribution in importlib.metadata.distributions():
        name = str(distribution.metadata.get("Name") or "").strip()
        version = str(distribution.version).strip()
        if not name or not version:
            continue
        normalized = re.sub(r"[-_.]+", "-", name).lower()
        installation = str(Path(distribution.locate_file("")).resolve())
        raw_metadata_path = getattr(distribution, "_path", None)
        metadata_path = (
            str(Path(raw_metadata_path).resolve())
            if raw_metadata_path
            else "unavailable"
        )
        distributions.setdefault(normalized, []).append(
            (version, installation, metadata_path)
        )
    if not distributions:
        raise CoordinatorContractError(
            "core/coordinator installed-package export is empty"
        )
    rows: list[dict[str, object]] = []
    for package in sorted(distributions):
        visible = sorted(set(distributions[package]))
        versions = {item[0] for item in visible}
        status = (
            "resolved_unique"
            if len(visible) == 1
            else "duplicate_same_version"
            if len(versions) == 1
            else "ambiguous_multiple_versions"
        )
        for version, installation, metadata_path in visible:
            identity = canonical_json_sha256(
                [package, version, installation, metadata_path]
            )[:16]
            rows.append(
                {
                    "environment_id": "plan60_core_coordinator_python",
                    "distribution_id": identity,
                    "package": package,
                    "version": version,
                    "installation_location": installation,
                    "metadata_path": metadata_path,
                    "resolution_status": status,
                }
            )
    return rows


def core_sys_path_rows() -> list[dict[str, object]]:
    rows = []
    for order, raw_path in enumerate(sys.path):
        resolved = Path(raw_path or ".").resolve()
        rows.append(
            {
                "environment_id": "plan60_core_coordinator_python",
                "path_order": order,
                "path": str(resolved),
                "exists": str(resolved.exists()).lower(),
                "is_directory": str(resolved.is_dir()).lower(),
            }
        )
    if not rows:
        raise CoordinatorContractError("core/coordinator sys.path export is empty")
    return rows


def core_import_resolution_rows() -> list[dict[str, object]]:
    rows = []
    for module_name in ("pyarrow",):
        try:
            module = importlib.import_module(module_name)
            metadata_version = importlib.metadata.version(module_name)
        except (ImportError, importlib.metadata.PackageNotFoundError) as error:
            raise CoordinatorContractError(
                f"required coordinator runtime dependency is unavailable: {module_name}"
            ) from error
        imported = Path(str(getattr(module, "__file__", ""))).resolve()
        if not imported.is_file() or imported.is_symlink():
            raise CoordinatorContractError(
                f"imported coordinator dependency has no safe file: {module_name}"
            )
        imported_version = str(getattr(module, "__version__", "")).strip()
        if not imported_version:
            raise CoordinatorContractError(
                f"imported coordinator dependency lacks __version__: {module_name}"
            )
        rows.append(
            {
                "environment_id": "plan60_core_coordinator_python",
                "module": module_name,
                "imported_version": imported_version,
                "imported_path": str(imported),
                "imported_file_sha256": sha256_file(imported),
                "metadata_version": metadata_version,
            }
        )
    return rows


def recursive_release_producer_rows(project_root: Path) -> list[dict[str, object]]:
    project = project_root.resolve()
    program_root = project / "scripts/manuscript/program_context_v2"
    if not program_root.is_dir() or program_root.is_symlink():
        raise CoordinatorContractError(
            f"Plan60 release-program root is missing or unsafe: {program_root}"
        )
    permitted_suffixes = {".py", ".R", ".sh", ".sbatch"}
    paths = [
        path
        for path in sorted(program_root.rglob("*"))
        if path.is_file()
        and not path.is_symlink()
        and path.suffix in permitted_suffixes
        and "__pycache__" not in path.parts
        and "prepared" not in path.parts
    ]
    passport_producers = [
        project / relative for relative in PLAN50_FROZEN_PRODUCER_PATHS
    ]
    unsafe_passport_producers = [
        str(path)
        for path in passport_producers
        if not path.is_file() or path.is_symlink()
    ]
    if unsafe_passport_producers:
        raise CoordinatorContractError(
            "Plan50 frozen producer universe is missing or unsafe: "
            f"{unsafe_passport_producers}"
        )
    nmf_freezer = (
        project
        / "Analysis/SingleCell/scripts/hotspot_modules"
        / "519_freeze_nmf_continuous_supplement.py"
    )
    if not nmf_freezer.is_file() or nmf_freezer.is_symlink():
        raise CoordinatorContractError(
            f"current NMF rederivation producer is missing or unsafe: {nmf_freezer}"
        )
    paths = sorted(set(paths + passport_producers + [nmf_freezer]))
    if not paths:
        raise CoordinatorContractError("recursive Plan60 producer inventory is empty")
    rows = []
    for path in paths:
        relative = project_relative(project, path)
        rows.append(
            {
                "producer_id": relative.replace("/", "__"),
                "repository_path": relative,
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
        )
    required = {
        "scripts/manuscript/program_context_v2/publish_plan60_workstream_handoff.py",
        "scripts/manuscript/program_context_v2/release_common.py",
        "scripts/manuscript/program_context_v2/coordinator/coordinator_contract.py",
        "scripts/manuscript/program_context_v2/coordinator/prepare_real_release.py",
        "scripts/manuscript/program_context_v2/coordinator/build_scientific_validation_registry.py",
        "Analysis/SingleCell/scripts/hotspot_modules/519_freeze_nmf_continuous_supplement.py",
        *PLAN50_FROZEN_PRODUCER_PATHS,
    }
    observed_paths = {str(row["repository_path"]) for row in rows}
    if not required.issubset(observed_paths):
        raise CoordinatorContractError(
            f"recursive producer inventory lacks: {sorted(required - observed_paths)}"
        )
    return rows


def load_json_object(path: Path, context: str) -> dict[str, object]:
    if not path.is_file() or path.is_symlink():
        raise CoordinatorContractError(f"{context} is missing or symlinked: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise CoordinatorContractError(f"{context} is invalid JSON: {path}") from error
    if not isinstance(payload, dict):
        raise CoordinatorContractError(f"{context} must be a JSON object: {path}")
    return payload


def validate_owner_signature(
    signature_path: Path,
    workstream_id: str,
    manifest_sha256: str,
) -> dict[str, object]:
    payload = load_json_object(signature_path, f"{workstream_id} handoff signature")
    required = {
        "attestation_version",
        "candidate_id",
        "workstream_id",
        "manifest_sha256",
        "signed_by",
        "signed_at_utc",
        "canonical_promotion_authorized",
    }
    if set(payload) != required:
        raise CoordinatorContractError(
            f"{workstream_id} handoff signature schema mismatch: "
            f"missing={sorted(required - set(payload))}, extra={sorted(set(payload) - required)}"
        )
    expected = {
        "attestation_version": HANDOFF_SIGNATURE_VERSION,
        "candidate_id": CANDIDATE_ID,
        "workstream_id": workstream_id,
        "manifest_sha256": manifest_sha256,
        "canonical_promotion_authorized": False,
    }
    for key, value in expected.items():
        if payload[key] != value:
            raise CoordinatorContractError(
                f"{workstream_id} handoff signature {key} mismatch: {payload[key]!r}"
            )
    if type(payload["canonical_promotion_authorized"]) is not bool:
        raise CoordinatorContractError(
            f"{workstream_id} canonical_promotion_authorized must be a JSON Boolean"
        )
    if (
        not str(payload["signed_by"]).strip()
        or not str(payload["signed_at_utc"]).strip()
    ):
        raise CoordinatorContractError(
            f"{workstream_id} handoff signature lacks signer/time"
        )
    return payload


def require_manifest_artifact(
    handoff: Handoff,
    source_path: Path,
    role: str,
) -> None:
    source = source_path.resolve()
    matches = [
        artifact
        for artifact in handoff.artifacts
        if artifact.source_path.resolve() == source and artifact.artifact_role == role
    ]
    if len(matches) != 1:
        raise CoordinatorContractError(
            f"{handoff.workstream_id} manifest requires exactly one {role} row for "
            f"{source_path}; observed {len(matches)}"
        )


def require_exact_manifest_artifacts(
    handoff: Handoff,
    expected: Sequence[tuple[Path, str]],
) -> None:
    """Require one exact, workstream-owned artifact set and no hidden extras."""

    observed = Counter(
        (artifact.source_path.resolve(), artifact.artifact_role)
        for artifact in handoff.artifacts
    )
    required = Counter((source.resolve(), role) for source, role in expected)
    if observed != required:
        missing = sorted(
            f"{path}:{role}"
            for (path, role), count in (required - observed).items()
            for _ in range(count)
        )
        extra = sorted(
            f"{path}:{role}"
            for (path, role), count in (observed - required).items()
            for _ in range(count)
        )
        raise CoordinatorContractError(
            f"{handoff.workstream_id} Plan60 artifact set drift: "
            f"missing={missing}, extra={extra}"
        )
    for source, role in expected:
        require_manifest_artifact(handoff, source, role)


def published_artifact_pairs(
    paths: RealPaths, workstream_id: str
) -> tuple[tuple[Path, str], ...]:
    root, specifications = workstream_contract(paths.project_root)[workstream_id]
    if root.resolve() != paths.workstream_root(workstream_id).resolve():
        raise CoordinatorContractError(
            f"shared publisher/coordinator root drift for {workstream_id}"
        )
    pairs = tuple((root / relative, role) for _, relative, role, _ in specifications)
    if len(pairs) != len(set((path.resolve(), role) for path, role in pairs)):
        raise CoordinatorContractError(
            f"shared publisher contract duplicates a path/role for {workstream_id}"
        )
    return pairs


def validate_handoff(paths: RealPaths, workstream_id: str) -> Handoff:
    manifest = paths.owner_manifest(workstream_id)
    signature = paths.owner_signature(workstream_id)
    roots = allowed_workstream_roots(paths.project_root, workstream_id)
    require_within_any(manifest.resolve(), roots, f"{workstream_id} owner manifest")
    if not manifest.is_file() or manifest.is_symlink():
        raise CoordinatorContractError(
            f"{workstream_id} workstream-owned Plan60 manifest is missing: {manifest}"
        )
    digest = sha256_file(manifest)
    validate_owner_signature(signature, workstream_id, digest)
    artifacts = validate_artifact_manifest(paths.project_root, workstream_id, manifest)
    return Handoff(
        workstream_id=workstream_id,
        root=paths.workstream_root(workstream_id),
        manifest=manifest,
        signature=signature,
        manifest_sha256=digest,
        artifacts=artifacts,
    )


def validate_nmf_continuous_bundle(
    paths: RealPaths,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Rebuild Figure S2 from the frozen Plan20 source copies and compare bytes."""

    root = paths.plan20_root / "nmf_continuous_supplement"
    ready = one_row_tsv(root / "NMF_CONTINUOUS_SUPPLEMENT_READY", NMF_READY_FIELDS)
    expected_ready = {
        "candidate_id": CANDIDATE_ID,
        "status": "continuous_nmf_supplement_frozen",
        "n_samples": "1104",
        "n_loading_rows": "11040",
        "n_evidence_rows": "18",
        "n_s2a_marks": "10",
        "n_s2b_marks": "8",
        "canonical_promotion_authorized": "false",
    }
    for field, expected in expected_ready.items():
        if ready[field] != expected:
            raise CoordinatorContractError(
                f"Plan20 continuous NMF READY {field} drift: {ready[field]!r}"
            )

    output_paths = {
        "source_manifest_sha256": root / "source_manifest.tsv",
        "loadings_sha256": root / "nmf_continuous_loadings.tsv",
        "evidence_sha256": root / "nmf_continuous_supplement.tsv",
        "validation_sha256": root / "validation_report.tsv",
        "producer_manifest_sha256": root / "producer_manifest.tsv",
    }
    for hash_field, path in output_paths.items():
        if (
            not path.is_file()
            or path.is_symlink()
            or sha256_file(path) != ready[hash_field]
        ):
            raise CoordinatorContractError(
                f"Plan20 continuous NMF {hash_field} drifted after READY"
            )

    sources = read_tsv_exact(root / "source_manifest.tsv", NMF_SOURCE_MANIFEST_FIELDS)
    expected_sources = {
        "k4_loadings": (
            "RNA-seq/results/subtypes/nmf_assignments_k4_pre_k6restore.csv",
            "source_inputs/nmf_assignments_k4_pre_k6restore.csv",
        ),
        "k6_loadings": (
            "RNA-seq/results/subtypes/nmf_assignments.csv",
            "source_inputs/nmf_assignments_k6.csv",
        ),
        "k4_labels": (
            "RNA-seq/results/subtypes/program_labels_k4_pre_k6restore.csv",
            "source_inputs/program_labels_k4_pre_k6restore.csv",
        ),
        "k6_labels": (
            "RNA-seq/results/subtypes/program_labels.csv",
            "source_inputs/program_labels_k6.csv",
        ),
        "three_seed_metrics": (
            "RNA-seq/results/subtypes/nmf_3seed_canonical/per_seed_metric_summary.csv",
            "source_inputs/per_seed_metric_summary.csv",
        ),
        "three_seed_stability": (
            "RNA-seq/results/subtypes/nmf_3seed_canonical/cross_seed_stability_summary.csv",
            "source_inputs/cross_seed_stability_summary.csv",
        ),
    }
    if {row["source_id"] for row in sources} != set(expected_sources) or len(
        sources
    ) != 6:
        raise CoordinatorContractError("Plan20 continuous NMF source universe drift")
    source_by_id = {row["source_id"]: row for row in sources}
    for source_id, (origin_relative, bundle_relative) in expected_sources.items():
        row = source_by_id[source_id]
        if (
            row["origin_project_path"] != origin_relative
            or row["bundle_relative_path"] != bundle_relative
            or row["origin_sha256"] != row["bundle_sha256"]
            or row["origin_bytes"] != row["bundle_bytes"]
        ):
            raise CoordinatorContractError(
                f"Plan20 continuous NMF source contract drift: {source_id}"
            )
        origin = paths.project_root / origin_relative
        bundle = root / bundle_relative
        for source_path, digest, byte_text, label in (
            (origin, row["origin_sha256"], row["origin_bytes"], "origin"),
            (bundle, row["bundle_sha256"], row["bundle_bytes"], "bundle"),
        ):
            if (
                not source_path.is_file()
                or source_path.is_symlink()
                or source_path.stat().st_size != int(byte_text)
                or sha256_file(source_path) != digest
            ):
                raise CoordinatorContractError(
                    f"Plan20 continuous NMF {label} bytes drift: {source_id}"
                )

    validation = read_tsv_exact(root / "validation_report.tsv", NMF_VALIDATION_FIELDS)
    expected_checks = {
        "sample_universe",
        "continuous_loading_rows",
        "supplement_evidence_rows",
        "figure_s2a_marks",
        "figure_s2b_marks",
        "supplement_only",
        "continuous_only_semantics",
    }
    if {row["check_id"] for row in validation} != expected_checks or any(
        row["status"] != "pass" for row in validation
    ):
        raise CoordinatorContractError("Plan20 continuous NMF validation is incomplete")

    # This is immutable *historical* provenance for the code that created the
    # Plan20 freeze.  Its exact bytes are already sealed by READY above and by
    # the Plan20 owner handoff in normal release preparation.  The historical
    # hashes must remain well-formed and attached to the original IDs, but are
    # deliberately not compared with later live coordinator bytes.
    historical_producers = read_tsv_exact(
        root / "producer_manifest.tsv", CORE_PRODUCER_FIELDS
    )
    expected_historical_producers = {
        "Analysis/SingleCell/scripts/hotspot_modules/519_freeze_nmf_continuous_supplement.py": "plan20_nmf_freezer",
        "scripts/manuscript/program_context_v2/coordinator/coordinator_contract.py": "plan60_nmf_contract",
    }
    if (
        len(historical_producers) != len(expected_historical_producers)
        or {row["repository_path"]: row["producer_id"] for row in historical_producers}
        != expected_historical_producers
    ):
        raise CoordinatorContractError(
            "Plan20 continuous NMF historical producer universe drift"
        )
    for row in historical_producers:
        try:
            producer_bytes = int(row["bytes"])
        except ValueError as error:
            raise CoordinatorContractError(
                "Plan20 continuous NMF historical producer byte count is invalid: "
                f"{row['producer_id']}"
            ) from error
        if re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is None or producer_bytes <= 0:
            raise CoordinatorContractError(
                "Plan20 continuous NMF historical producer record is invalid: "
                f"{row['producer_id']}"
            )

    # Current rederivation code has a separate temporal identity.  Require it
    # in the recursive Plan60 producer universe and independently verify its
    # live bytes.  prepare_real_release later binds these current hashes, the
    # frozen validation-report inputs, and the rederived adapter output in one
    # exact adapter-provenance row.
    current_producers = recursive_release_producer_rows(paths.project_root)
    current_by_path = {str(row["repository_path"]): row for row in current_producers}
    required_current_producers = set(expected_historical_producers)
    if not required_current_producers.issubset(current_by_path):
        raise CoordinatorContractError(
            "Plan20 continuous NMF current producer universe incomplete: "
            f"{sorted(required_current_producers - set(current_by_path))}"
        )
    for relative in sorted(required_current_producers):
        row = current_by_path[relative]
        producer = paths.project_root / relative
        if (
            not producer.is_file()
            or producer.is_symlink()
            or producer.stat().st_size != int(row["bytes"])
            or sha256_file(producer) != row["sha256"]
        ):
            raise CoordinatorContractError(
                f"Plan20 continuous NMF current producer drift: {relative}"
            )

    def source_rows(source_id: str) -> list[dict[str, str]]:
        _, rows = read_tsv_flexible(
            root / source_by_id[source_id]["bundle_relative_path"]
        )
        return rows

    rebuilt_loadings, rebuilt_evidence = build_nmf_continuous_supplement(
        source_rows("k4_loadings"),
        source_rows("k6_loadings"),
        source_rows("k4_labels"),
        source_rows("k6_labels"),
        source_rows("three_seed_metrics"),
        source_rows("three_seed_stability"),
        source_by_id["k4_loadings"]["bundle_sha256"],
        source_by_id["k6_loadings"]["bundle_sha256"],
    )
    frozen_loadings = read_tsv_exact(
        root / "nmf_continuous_loadings.tsv", NMF_CONTINUOUS_FIELDS
    )
    frozen_evidence = read_tsv_exact(
        root / "nmf_continuous_supplement.tsv", SOURCE_EVIDENCE_FIELDS
    )
    normalized_loadings = [
        {field: str(row[field]) for field in NMF_CONTINUOUS_FIELDS}
        for row in rebuilt_loadings
    ]
    normalized_evidence = [
        {field: str(row[field]) for field in SOURCE_EVIDENCE_FIELDS}
        for row in rebuilt_evidence
    ]
    if frozen_loadings != normalized_loadings or frozen_evidence != normalized_evidence:
        raise CoordinatorContractError(
            "Plan20 continuous NMF frozen products do not rederive from source copies"
        )
    validate_source_rows(frozen_evidence, "continuous_nmf_handoff")
    if (
        len(frozen_evidence) != 18
        or sum(
            row["panel_id"] == "S2A" and row["plot_role"] == "mark"
            for row in frozen_evidence
        )
        != 10
        or sum(
            row["panel_id"] == "S2B" and row["plot_role"] == "mark"
            for row in frozen_evidence
        )
        != 8
        or any(row["manuscript_included"] != "false" for row in frozen_evidence)
    ):
        raise CoordinatorContractError("Plan20 continuous NMF FigureS2 contract drift")
    return frozen_loadings, frozen_evidence


def validate_plan20(paths: RealPaths, handoff: Handoff) -> None:
    root = paths.plan20_root
    ready_fields = (
        "release_id",
        "status",
        "validation_status_sha256",
        "release_manifest_sha256",
        "registry_sha256",
        "membership_table_sha256",
        "validator_sha256",
        "external_outcomes_read",
        "sealed_at_utc",
    )
    semantic_fields = (
        "release_id",
        "status",
        "registry_sha256",
        "adjudication_sha256",
        "manifest_sha256",
        "validation_sha256",
        "producer_sha256",
        "validator_sha256",
        "n_programs",
        "n_q_significant",
        "n_q_significant_stability_failed",
        "n_primary_selected",
        "n_robust_display",
        "n_indeterminate_nonconfirmatory",
        "n_tested_negative_authorized",
        "validated_at_utc",
    )
    ready = one_row_tsv(root / "READY", ready_fields)
    semantic = one_row_tsv(root / "SEMANTIC_ADJUDICATION_READY", semantic_fields)
    composition_ready = one_row_tsv(
        root / "COMPOSITION_QC_READY", COMPOSITION_READY_FIELDS
    )
    if (
        ready["release_id"] != CANDIDATE_ID
        or ready["status"] != "ready_for_external_testing"
    ):
        raise CoordinatorContractError("Plan20 READY release/status drift")
    if parse_bool(ready["external_outcomes_read"], "Plan20 external_outcomes_read"):
        raise CoordinatorContractError("Plan20 READY is not outcome blind")
    if semantic["status"] != "semantic_adjudication_validated":
        raise CoordinatorContractError("Plan20 semantic adjudication is not terminal")
    if semantic["registry_sha256"] != ready["registry_sha256"]:
        raise CoordinatorContractError("Plan20 READY/semantic registry hashes disagree")
    if int(semantic["n_programs"]) != 117 or int(semantic["n_robust_display"]) != 2:
        raise CoordinatorContractError(
            "Plan20 reviewed 117-program/two-robust contract drift"
        )
    if int(semantic["n_tested_negative_authorized"]) != 0:
        raise CoordinatorContractError("Plan20 authorizes a tested-negative call")
    if sha256_file(root / "program_registry_v2.tsv") != ready["registry_sha256"]:
        raise CoordinatorContractError("Plan20 registry bytes drifted after READY")
    if (
        sha256_file(root / "program_membership_v2.tsv")
        != ready["membership_table_sha256"]
    ):
        raise CoordinatorContractError("Plan20 membership bytes drifted after READY")
    if (
        sha256_file(root / "program_registry_v2_semantic_adjudication.tsv")
        != semantic["adjudication_sha256"]
    ):
        raise CoordinatorContractError(
            "Plan20 semantic table drifted after adjudication"
        )
    _, registry = read_tsv_flexible(root / "program_registry_v2.tsv")
    _, figure_source = read_tsv_flexible(root / "fig2_program_source.tsv")
    validate_hotspot_registry_projection(registry, figure_source)

    composition_results = read_tsv_exact(
        root / "composition_sample_qc_v2.tsv", COMPOSITION_RESULT_FIELDS
    )
    composition_testability = read_tsv_exact(
        root / "composition_sample_qc_testability.tsv", COMPOSITION_TESTABILITY_FIELDS
    )
    composition_audit = read_tsv_exact(
        root / "composition_sample_qc_audit.tsv",
        ("audit_id", "status", "value", "detail"),
    )
    composition_validation = read_tsv_exact(
        root / "composition_sample_qc_validation.tsv",
        ("check_id", "status", "observed", "criterion"),
    )
    composition_paths = {
        "results_sha256": root / "composition_sample_qc_v2.tsv",
        "testability_sha256": root / "composition_sample_qc_testability.tsv",
        "audit_sha256": root / "composition_sample_qc_audit.tsv",
        "validation_sha256": root / "composition_sample_qc_validation.tsv",
    }
    for hash_field, path in composition_paths.items():
        if sha256_file(path) != composition_ready[hash_field]:
            raise CoordinatorContractError(
                f"Plan20 composition {hash_field} drifted after READY"
            )
    validate_composition_contract(
        composition_results,
        composition_testability,
        composition_ready,
        composition_audit,
        composition_validation,
    )
    validate_nmf_continuous_bundle(paths)
    require_exact_manifest_artifacts(handoff, published_artifact_pairs(paths, "PLAN20"))


def validate_plan30(paths: RealPaths, handoff: Handoff) -> None:
    root = paths.plan30_root
    fields = (
        "release_id",
        "status",
        "closure_sha256",
        "terminal_gate_sha256",
        "terminal_validation_sha256",
        "historical_preflight_sha256",
        "closure_manifest_sha256",
        "closure_validation_sha256",
        "validator_sha256",
        "source_wide_unique_source_defined_egenes",
        "mapped_symbol_annotation_rows",
        "represented_unique_source_ensgs",
        "context_rescue_authorized",
        "negative_claim_authorized",
        "canonical_promotion_status",
        "validated_at_utc",
    )
    ready = one_row_tsv(root / "GEN_TERMINAL_CLOSURE_READY", fields)
    expected = {
        "release_id": CANDIDATE_ID,
        "status": "coverage_limited_terminal_validated",
        "source_wide_unique_source_defined_egenes": "6564",
        "mapped_symbol_annotation_rows": "6583",
        "represented_unique_source_ensgs": "6562",
        "context_rescue_authorized": "false",
        "negative_claim_authorized": "false",
        "canonical_promotion_status": "not_promoted",
    }
    for key, value in expected.items():
        if ready[key] != value:
            raise CoordinatorContractError(
                f"Plan30 terminal {key} drift: {ready[key]!r}"
            )
    if sha256_file(root / "terminal_closure.tsv") != ready["closure_sha256"]:
        raise CoordinatorContractError("Plan30 terminal closure drifted after READY")
    require_exact_manifest_artifacts(handoff, published_artifact_pairs(paths, "PLAN30"))


def validate_plan40(paths: RealPaths, handoff: Handoff) -> None:
    root = paths.plan40_root
    terminal_fields = (
        "release_id",
        "status",
        "specification_sha256",
        "phase_c_code_bundle_sha256",
        "validation_sha256",
        "external_outcomes_read",
        "biological_positivity_required",
    )
    verdict_fields = (
        "main_figure_eligible",
        "source_gate_pass",
        "class_omnibus_pass",
        "known_hit_class_omnibus_pass",
        "class_pairwise_branch_pass",
        "program_branch_pass",
        "known_hit_survival_pass",
        "class_minimum_standardized_effect",
        "program_minimum_signed_effect",
        "observed_qualifying_effect",
        "qualifying_class_contrasts",
        "qualifying_programs",
        "verdict_reason_codes",
        "specification_sha256",
        "phase_c_code_bundle_sha256",
        "biological_positivity_required_for_software_acceptance",
        "release_id",
    )
    terminal = one_row_tsv(root / "PHASE_C_VALIDATED", terminal_fields)
    verdict = one_row_tsv(root / "fig5_verdict.tsv", verdict_fields)
    if (
        terminal["status"] != "phase_c_release_validated"
        or terminal["release_id"] != CANDIDATE_ID
    ):
        raise CoordinatorContractError("Plan40 PHASE_C_VALIDATED status drift")
    if not parse_bool(terminal["external_outcomes_read"], "Plan40 outcomes read"):
        raise CoordinatorContractError(
            "Plan40 terminal does not record held-out outcome evaluation"
        )
    if parse_bool(
        terminal["biological_positivity_required"], "Plan40 positivity acceptance"
    ):
        raise CoordinatorContractError(
            "Plan40 software acceptance requires biological positivity"
        )
    if parse_bool(verdict["main_figure_eligible"], "Plan40 main figure eligibility"):
        raise CoordinatorContractError("Plan40 unexpectedly enters a main figure")
    if parse_bool(
        verdict["biological_positivity_required_for_software_acceptance"],
        "Plan40 software gate",
    ):
        raise CoordinatorContractError(
            "Plan40 fig5 verdict uses positivity as software acceptance"
        )
    if verdict["specification_sha256"] != terminal["specification_sha256"]:
        raise CoordinatorContractError("Plan40 terminal/verdict specification mismatch")
    require_exact_manifest_artifacts(handoff, published_artifact_pairs(paths, "PLAN40"))


def validate_plan13(paths: RealPaths, handoff: Handoff) -> None:
    root = paths.plan13_root
    fields = (
        "release_id",
        "semantic_contract_id",
        "status",
        "manifest_sha256",
        "final_ready_sha256",
        "yak_adapter_ready_sha256",
        "n_files",
        "n_robust",
        "n_indeterminate",
        "n_tested_negative",
        "historical_bundles_unchanged",
        "canonical_promotion_authorized",
    )
    ready = one_row_tsv(root / "SEMANTIC_V2_READY", fields)
    if sha256_file(root / "SEMANTIC_V2_READY") != EXPECTED_PLAN13_TERMINAL_SHA256:
        raise CoordinatorContractError("Plan13 semantic-v2 terminal seal hash drift")
    if ready["release_id"] != CANDIDATE_ID:
        raise CoordinatorContractError("Plan13 semantic-v2 release ID drift")
    if ready["status"] != "sealed_semantic_v2_candidate_no_canonical_promotion":
        raise CoordinatorContractError("Plan13 semantic-v2 seal is not terminal")
    if int(ready["n_tested_negative"]) != 0:
        raise CoordinatorContractError(
            "Plan13 semantic-v2 contains tested-negative rows"
        )
    if int(ready["n_robust"]) != 2 or int(ready["n_indeterminate"]) != 8:
        raise CoordinatorContractError(
            "Plan13 reviewed 2-robust/8-indeterminate contract drift"
        )
    if not parse_bool(
        ready["historical_bundles_unchanged"], "Plan13 historical bundles"
    ):
        raise CoordinatorContractError(
            "Plan13 historical bundles are not byte-identical"
        )
    if parse_bool(ready["canonical_promotion_authorized"], "Plan13 promotion"):
        raise CoordinatorContractError("Plan13 terminal seal authorizes promotion")
    manifest = root / "semantic_v2_release_manifest.tsv"
    if sha256_file(manifest) != ready["manifest_sha256"]:
        raise CoordinatorContractError("Plan13 semantic release manifest drift")
    final_root = root / "final_integration"
    if sha256_file(final_root / "READY") != ready["final_ready_sha256"]:
        raise CoordinatorContractError("Plan13 final READY drift")
    _, matrix = read_tsv_flexible(final_root / "figure4_program_matrix.tsv")
    if len(matrix) != 18:
        raise CoordinatorContractError(
            f"Plan13 Figure4 matrix must contain 18 rows; found {len(matrix)}"
        )
    if any(row.get("evidence_state") == "tested_negative" for row in matrix):
        raise CoordinatorContractError(
            "Plan13 Figure4 matrix reintroduces tested-negative semantics"
        )
    require_exact_manifest_artifacts(handoff, published_artifact_pairs(paths, "PLAN13"))


PASS50_FIELDS = (
    "analysis_release_id",
    "status",
    "selection_sha256",
    "manifest_sha256",
    "validation_report_sha256",
    "n_genes",
    "n_evidence_rows",
    "n_programs",
    "n_program_context_rows",
    "automated_validation",
    "manual_acceptance",
    "handoff_allowed",
    "canonical_promotion_authorized",
    "scientific_call_recomputed",
    "validated_at_utc",
)
PLAN50_MANUAL_ACCEPTANCE_FIELDS = (
    "review_id",
    "analysis_release_id",
    "selection_sha256",
    "reviewer",
    "reviewed_at_utc",
    "ui_index_sha256",
    "ui_contract_sha256",
    "ui_source_manifest_sha256",
    "ui_review_manifest_sha256",
    "call_states_reviewed",
    "provenance_states_reviewed",
    "hero_genes_reviewed",
    "visible_boundary_pass",
    "visible_testability_pass",
    "visible_source_dependence_pass",
    "visible_falsifier_pass",
    "decision",
)
PLAN50_TERMINAL_PROVENANCE_FIELDS = (
    "analysis_release_id",
    "bundle_uri",
    "artifact_role",
    "relative_path",
    "sha256",
    "bytes",
    "producer",
    "producer_sha256",
)
PLAN50_INNER_HANDOFF_FIELDS = (
    "analysis_release_id",
    "bundle_path",
    "manifest_sha256",
    "selection_sha256",
    "automated_validation",
    "manual_acceptance",
    "promotion_allowed",
    "build_command",
    "environment",
    "scientific_call_recomputed",
    "omitted_input_inventory",
)
PLAN50_VALIDATION_REPORT_FIELDS = ("check_id", "status", "detail")


def _plan50_frozen_producer_index(project_root: Path) -> dict[str, dict[str, object]]:
    rows = recursive_release_producer_rows(project_root)
    indexed = {str(row["repository_path"]): row for row in rows}
    if len(indexed) != len(rows):
        raise CoordinatorContractError("Plan50 frozen producer paths are duplicated")
    missing = set(PLAN50_FROZEN_PRODUCER_PATHS) - set(indexed)
    if missing:
        raise CoordinatorContractError(
            f"Plan50 frozen producer inventory lacks: {sorted(missing)}"
        )
    return indexed


def _validate_plan50_producer_binding(
    row: Mapping[str, str],
    frozen: Mapping[str, Mapping[str, object]],
    context: str,
) -> None:
    producer = row["producer"]
    clean = PurePosixPath(producer)
    if (
        clean.is_absolute()
        or ".." in clean.parts
        or producer not in PLAN50_MANIFEST_PRODUCER_PATHS
    ):
        raise CoordinatorContractError(
            f"{context} names a producer outside the Plan50 allowlist: {producer}"
        )
    current = frozen.get(producer)
    if current is None:
        raise CoordinatorContractError(
            f"{context} producer is absent from the frozen REL01 universe: {producer}"
        )
    expected_hash = require_sha256(row["producer_sha256"], f"{context} producer_sha256")
    if expected_hash != str(current["sha256"]):
        raise CoordinatorContractError(
            f"{context} producer hash differs from the frozen REL01 producer: {producer}"
        )


def _validate_plan50_acceptance_chain(
    root: Path,
    ready: Mapping[str, str],
    manifest_rows: Sequence[Mapping[str, str]],
    frozen_producers: Mapping[str, Mapping[str, object]],
) -> None:
    manifest = root / "passport_release_manifest.tsv"
    manifest_index = {row["relative_path"]: row for row in manifest_rows}
    required_payload = {
        "passport_manual_acceptance.tsv",
        "passport_manual_acceptance.signature.json",
        "passport_domain_coverage.tsv",
        "passport_assay_status.tsv",
        "portal_candidate/index.html",
        "portal_candidate/ui_contract.json",
        "portal_candidate/ui_source_manifest.tsv",
        "portal_candidate/review/review_manifest.tsv",
    }
    missing = required_payload - set(manifest_index)
    if missing:
        raise CoordinatorContractError(
            f"Plan50 manifest lacks acceptance/UI payload files: {sorted(missing)}"
        )

    manual_path = root / "passport_manual_acceptance.tsv"
    manual_rows = read_tsv_exact(manual_path, PLAN50_MANUAL_ACCEPTANCE_FIELDS)
    if len(manual_rows) != 1:
        raise CoordinatorContractError("Plan50 manual acceptance must have one row")
    manual = manual_rows[0]
    if (
        manual["analysis_release_id"] != PLAN50_ANALYSIS_RELEASE_ID
        or manual["selection_sha256"] != ready["selection_sha256"]
        or manual["decision"] != "accepted"
        or set(manual["hero_genes_reviewed"].split(";")) != set(FIXED_HERO_GENES)
    ):
        raise CoordinatorContractError("Plan50 manual acceptance semantics drift")
    for field in (
        "call_states_reviewed",
        "provenance_states_reviewed",
        "visible_boundary_pass",
        "visible_testability_pass",
        "visible_source_dependence_pass",
        "visible_falsifier_pass",
    ):
        if manual[field] != "true":
            raise CoordinatorContractError(
                f"Plan50 manual acceptance did not pass {field}"
            )
    ui_bindings = {
        "ui_index_sha256": "portal_candidate/index.html",
        "ui_contract_sha256": "portal_candidate/ui_contract.json",
        "ui_source_manifest_sha256": "portal_candidate/ui_source_manifest.tsv",
        "ui_review_manifest_sha256": ("portal_candidate/review/review_manifest.tsv"),
    }
    for field, relative in ui_bindings.items():
        expected = require_sha256(manual[field], f"Plan50 manual {field}")
        if expected != sha256_file(root / relative):
            raise CoordinatorContractError(
                f"Plan50 manual acceptance does not bind {relative}"
            )
    signature_path = root / "passport_manual_acceptance.signature.json"
    signature = load_json_object(signature_path, "Plan50 manual acceptance signature")
    required_signature = {
        "attestation_version",
        "acceptance_sha256",
        *ui_bindings,
        "signed_by",
        "signed_at_utc",
        "decision_register_id",
        "fixture_only",
    }
    if set(signature) != required_signature:
        raise CoordinatorContractError(
            "Plan50 manual acceptance signature schema drift"
        )
    if (
        signature["attestation_version"] != "passport_manual_acceptance_v2"
        or signature["acceptance_sha256"] != sha256_file(manual_path)
        or signature["fixture_only"] is not False
        or signature["signed_by"] != manual["reviewer"]
        or signature["signed_at_utc"] != manual["reviewed_at_utc"]
        or not str(signature["decision_register_id"]).strip()
    ):
        raise CoordinatorContractError("Plan50 manual acceptance signature drift")
    for field in ui_bindings:
        if signature[field] != manual[field]:
            raise CoordinatorContractError(
                f"Plan50 manual signature does not bind {field}"
            )

    inner_path = root / "passport_plan60_handoff.tsv"
    inner_rows = read_tsv_exact(inner_path, PLAN50_INNER_HANDOFF_FIELDS)
    if len(inner_rows) != 1:
        raise CoordinatorContractError("Plan50 inner handoff must have one row")
    inner = inner_rows[0]
    expected_inner = {
        "analysis_release_id": PLAN50_ANALYSIS_RELEASE_ID,
        "bundle_path": f"release://{PLAN50_ANALYSIS_RELEASE_ID}/",
        "manifest_sha256": ready["manifest_sha256"],
        "selection_sha256": ready["selection_sha256"],
        "automated_validation": "pass",
        "manual_acceptance": "passed_signed_manual_acceptance",
        "promotion_allowed": "false",
        "scientific_call_recomputed": "false",
    }
    for field, expected in expected_inner.items():
        if inner[field] != expected:
            raise CoordinatorContractError(
                f"Plan50 inner handoff {field} drift: {inner[field]!r}"
            )

    provenance_path = root / "passport_terminal_provenance.tsv"
    provenance_rows = read_tsv_exact(provenance_path, PLAN50_TERMINAL_PROVENANCE_FIELDS)
    provenance = {row["artifact_role"]: row for row in provenance_rows}
    if len(provenance_rows) != 2 or set(provenance) != {
        "payload_manifest",
        "plan60_handoff",
    }:
        raise CoordinatorContractError("Plan50 terminal provenance roles drift")
    expected_terminal = {
        "payload_manifest": ("passport_release_manifest.tsv", manifest),
        "plan60_handoff": ("passport_plan60_handoff.tsv", inner_path),
    }
    bundle_uri = f"release://{PLAN50_ANALYSIS_RELEASE_ID}/"
    for role, (relative, artifact) in expected_terminal.items():
        row = provenance[role]
        if (
            row["analysis_release_id"] != PLAN50_ANALYSIS_RELEASE_ID
            or row["bundle_uri"] != bundle_uri
            or row["relative_path"] != relative
            or row["sha256"] != sha256_file(artifact)
            or parse_nonnegative_int(row["bytes"], f"Plan50 {role} bytes")
            != artifact.stat().st_size
        ):
            raise CoordinatorContractError(f"Plan50 terminal provenance drift: {role}")
        _validate_plan50_producer_binding(
            row, frozen_producers, f"Plan50 terminal provenance {role}"
        )

    validation_path = root / "passport_validation_report.tsv"
    validation_rows = read_tsv_exact(validation_path, PLAN50_VALIDATION_REPORT_FIELDS)
    validation = {row["check_id"]: row for row in validation_rows}
    if len(validation) != len(validation_rows):
        raise CoordinatorContractError("Plan50 validation report has duplicate IDs")
    terminal_check = validation.get("PASS06_TERMINAL_PROVENANCE")
    if terminal_check is None or terminal_check["status"] != "pass":
        raise CoordinatorContractError(
            "Plan50 validation lacks a passing terminal-provenance check"
        )
    try:
        detail = json.loads(terminal_check["detail"])
    except json.JSONDecodeError as error:
        raise CoordinatorContractError(
            "Plan50 terminal-provenance validation detail is invalid JSON"
        ) from error
    if (
        detail.get("contract_version") != "passport_terminal_provenance_v1"
        or detail.get("terminal_provenance_path") != "passport_terminal_provenance.tsv"
        or detail.get("terminal_provenance_sha256") != sha256_file(provenance_path)
    ):
        raise CoordinatorContractError(
            "Plan50 validation report does not bind terminal provenance"
        )


def validate_plan50(paths: RealPaths, handoff: Handoff) -> None:
    root = paths.plan50_root
    ready = one_row_tsv(root / "PASS06_VALIDATED", PASS50_FIELDS)
    expected = {
        "analysis_release_id": PLAN50_ANALYSIS_RELEASE_ID,
        "status": "passport_release_validated",
        "automated_validation": "pass",
        "manual_acceptance": "passed_signed_manual_acceptance",
        "handoff_allowed": "true",
        "canonical_promotion_authorized": "false",
        "scientific_call_recomputed": "false",
    }
    for key, value in expected.items():
        if ready[key] != value:
            raise CoordinatorContractError(
                f"Plan50 terminal {key} drift: {ready[key]!r}"
            )
    manifest = root / "passport_release_manifest.tsv"
    if sha256_file(manifest) != ready["manifest_sha256"]:
        raise CoordinatorContractError("Plan50 release manifest drift")
    rows = read_tsv_exact(manifest, PLAN50_PAYLOAD_MANIFEST_FIELDS)
    indexed = {row["relative_path"]: row for row in rows}
    if len(indexed) != len(rows):
        raise CoordinatorContractError(
            "Plan50 release manifest contains duplicate paths"
        )
    required_files = {
        "passport_gene_index.parquet",
        "passport_evidence_long.parquet",
        "passport_coverage_long.parquet",
        "passport_program_context.parquet",
        "passport_next_experiment.tsv",
        "portal_candidate/index.html",
        "portal_candidate/ui_contract.json",
        "portal_candidate/ui_source_manifest.tsv",
        "portal_candidate/review/overview.png",
        "portal_candidate/review/THRB.png",
        "portal_candidate/review/HKDC1.png",
        "portal_candidate/review/GLP1R.png",
        "portal_candidate/review/MTARC1.png",
        "portal_candidate/review/review_manifest.tsv",
    }
    missing = required_files - set(indexed)
    if missing:
        raise CoordinatorContractError(
            f"Plan50 release manifest lacks: {sorted(missing)}"
        )
    if any(row["release_status"] != "validated_candidate_handoff" for row in rows):
        raise CoordinatorContractError(
            "Plan50 manifest includes a nonterminal release status"
        )
    frozen_producers = _plan50_frozen_producer_index(paths.project_root)
    observed_producers = {row["producer"] for row in rows}
    if observed_producers != set(PLAN50_MANIFEST_PRODUCER_PATHS):
        raise CoordinatorContractError(
            "Plan50 payload producer universe drift: "
            f"observed={sorted(observed_producers)}"
        )
    for relative, row in indexed.items():
        clean = PurePosixPath(relative)
        if clean.is_absolute() or ".." in clean.parts:
            raise CoordinatorContractError(
                f"Plan50 manifest path is unsafe: {relative}"
            )
        source = root / relative
        if (
            not source.is_file()
            or source.is_symlink()
            or source.stat().st_size != int(row["bytes"])
            or sha256_file(source) != row["sha256"]
        ):
            raise CoordinatorContractError(
                f"Plan50 manifested artifact drift: {source}"
            )
        _validate_plan50_producer_binding(
            row,
            frozen_producers,
            f"Plan50 payload manifest {relative}",
        )
    selection = root / "passport_input_selection.tsv"
    validation_report = root / "passport_validation_report.tsv"
    if (
        not selection.is_file()
        or selection.is_symlink()
        or sha256_file(selection) != ready["selection_sha256"]
        or "passport_input_selection.tsv" not in indexed
    ):
        raise CoordinatorContractError("Plan50 PASS06 input-selection hash drift")
    if (
        not validation_report.is_file()
        or validation_report.is_symlink()
        or sha256_file(validation_report) != ready["validation_report_sha256"]
    ):
        raise CoordinatorContractError("Plan50 PASS06 validation-report hash drift")
    _validate_plan50_acceptance_chain(root, ready, rows, frozen_producers)
    require_exact_manifest_artifacts(handoff, published_artifact_pairs(paths, "PLAN50"))


def validate_all_handoffs(paths: RealPaths) -> dict[str, Handoff]:
    handoffs = {
        workstream: validate_handoff(paths, workstream)
        for workstream in WORKSTREAM_ORDER
    }
    validate_plan13(paths, handoffs["PLAN13"])
    validate_plan20(paths, handoffs["PLAN20"])
    validate_plan30(paths, handoffs["PLAN30"])
    validate_plan40(paths, handoffs["PLAN40"])
    validate_plan50(paths, handoffs["PLAN50"])
    return handoffs


def source_row(**updates: object) -> dict[str, object]:
    row: dict[str, object] = {field: "" for field in SOURCE_EVIDENCE_FIELDS}
    row.update(updates)
    return row


def source_dependency(
    dependence_class: str,
    discovery_sources: Iterable[str],
    evaluation_sources: Iterable[str],
    reuse_detail: str,
    independence_boundary: str,
) -> dict[str, str]:
    """Build and immediately validate an explicit discovery/evaluation edge."""

    discovery = ";".join(sorted(set(discovery_sources)))
    evaluation = ";".join(sorted(set(evaluation_sources)))
    validate_source_dependency_contract(
        dependence_class,
        discovery,
        evaluation,
        reuse_detail,
        independence_boundary,
        "coordinator source-dependency adapter",
    )
    return {
        "source_dependence": dependence_class,
        "discovery_sources": discovery,
        "evaluation_sources": evaluation,
        "reuse_detail": reuse_detail,
        "independence_boundary": independence_boundary,
    }


def validate_source_rows(rows: Sequence[Mapping[str, object]], context: str) -> None:
    if not rows:
        raise CoordinatorContractError(f"{context} evidence-row adapter is empty")
    ids: set[str] = set()
    claim_dependencies: dict[str, tuple[str, str, str, str, str]] = {}
    for row in rows:
        if tuple(row) != SOURCE_EVIDENCE_FIELDS:
            raise CoordinatorContractError(f"{context} source-row schema drift")
        record_id = str(row["record_id"])
        if not record_id or record_id in ids:
            raise CoordinatorContractError(
                f"{context} duplicate/blank record ID: {record_id!r}"
            )
        ids.add(record_id)
        parse_float(row["value"], f"{context}:{record_id}:value")
        for field in (
            "group_id",
            "artifact_id",
            "claim_id",
            "section_id",
            "figure_id",
            "panel_id",
            "panel_title",
            "panel_role",
            "label",
            "number_role",
            "plot_role",
            "display_value",
            "unit",
            "biological_unit",
            "model_contrast",
            "effect_unit",
            "evidence_status",
            "tested_universe",
            "source_dependence",
            "discovery_sources",
            "evaluation_sources",
            "reuse_detail",
            "independence_boundary",
            "allowed_wording",
            "prohibited_wording",
            "claim_text",
        ):
            if not str(row[field]).strip():
                raise CoordinatorContractError(f"{context}:{record_id} blank {field}")
        validate_source_dependency_contract(
            str(row["source_dependence"]),
            str(row["discovery_sources"]),
            str(row["evaluation_sources"]),
            str(row["reuse_detail"]),
            str(row["independence_boundary"]),
            f"{context}:{record_id}",
        )
        dependency = (
            str(row["source_dependence"]),
            str(row["discovery_sources"]),
            str(row["evaluation_sources"]),
            str(row["reuse_detail"]),
            str(row["independence_boundary"]),
        )
        claim_id = str(row["claim_id"])
        if (
            claim_id in claim_dependencies
            and claim_dependencies[claim_id] != dependency
        ):
            raise CoordinatorContractError(
                f"{context}:{claim_id} has inconsistent source-dependency edges"
            )
        claim_dependencies[claim_id] = dependency
        status = str(row["evidence_status"])
        rule = str(row["negative_adequacy_criterion"]).strip()
        if status == "tested_negative" and not rule:
            raise CoordinatorContractError(
                f"{context}:{record_id} tested_negative lacks an adequate-negative rule"
            )
        if status != "tested_negative" and rule:
            raise CoordinatorContractError(
                f"{context}:{record_id} has a negative rule without tested_negative status"
            )
        text = " ".join(
            str(row[field]).lower()
            for field in (
                "label",
                "category",
                "number_role",
                "model_contrast",
                "effect_unit",
                "allowed_wording",
                "claim_text",
            )
        )
        if any(token in text for token in FORBIDDEN_PRESENTATION_TOKENS):
            raise CoordinatorContractError(
                f"{context}:{record_id} contains a forbidden score/rank token"
            )
        if "broad validation" in text:
            raise CoordinatorContractError(
                f"{context}:{record_id} overstates external evidence"
            )


def build_cohort_overview_rows(
    metadata_rows: Sequence[Mapping[str, str]],
    qc_rows: Sequence[Mapping[str, str]],
    pooled_datasets: set[str],
) -> list[dict[str, object]]:
    metadata = {row["sample_id"]: row for row in metadata_rows}
    if len(metadata) != len(metadata_rows):
        raise CoordinatorContractError("cohort metadata has duplicate sample_id values")
    qc = {row["sample_id"]: row for row in qc_rows}
    if len(qc) != len(qc_rows):
        raise CoordinatorContractError("QC table has duplicate sample_id values")
    eligible = [
        metadata[sample]
        for sample, row in qc.items()
        if parse_bool(row["pass_technical"], f"QC {sample} pass_technical")
        and sample in metadata
        and metadata[sample]["dataset"] != "PRJNA512027"
    ]
    observed = {
        "qc_samples": len(eligible),
        "qc_cohorts": len({row["dataset"] for row in eligible}),
        "pooled_samples": sum(row["dataset"] in pooled_datasets for row in eligible),
        "pooled_cohorts": len(
            {row["dataset"] for row in eligible if row["dataset"] in pooled_datasets}
        ),
    }
    if observed != EXPECTED_COHORT_COUNTS:
        raise CoordinatorContractError(
            "cohort source rederivation disagrees with the canonical numerical authority: "
            f"observed={observed}, expected={EXPECTED_COHORT_COUNTS}"
        )
    rows: list[dict[str, object]] = []
    census_dependency = source_dependency(
        "reused_source",
        ("canonical_human_metadata_v2026_08_08", "canonical_sample_qc_v2026_08_08"),
        ("canonical_human_metadata_v2026_08_08", "canonical_sample_qc_v2026_08_08"),
        "The resource census and displayed denominator are rederived from the same frozen metadata and QC tables.",
        "This is a descriptive inventory, not an independent validation analysis.",
    )
    observability_dependency = source_dependency(
        "mixed",
        ("paper_directional_plan_2026_08_07", "program_context_v1"),
        ("PLAN13", "PLAN20", "PLAN30", "PLAN40", "PLAN50"),
        "The interface reuses the prespecified thesis while the five workstreams supply the displayed role boundaries.",
        "No workstream result is treated as an independent evaluation of the directional plan itself.",
    )
    inventory = (
        ("qc_samples", "QC-passing human samples", "human biological sample"),
        ("qc_cohorts", "Human cohorts", "human cohort"),
        ("pooled_samples", "Pooled disease-control samples", "human biological sample"),
        ("pooled_cohorts", "Control-bearing pooled cohorts", "human cohort"),
    )
    for order, (key, label, unit) in enumerate(inventory, 1):
        rows.append(
            source_row(
                record_id=f"resource_{key}",
                group_id=f"resource_{key}",
                artifact_id="cohort_overview",
                claim_id=f"claim_resource_{key}",
                section_id="resource_interface",
                figure_id="Figure1",
                panel_id="1A",
                panel_title="Cohort and analysis-unit overview",
                panel_role="interface",
                label=label,
                category="resource inventory",
                number_role=key,
                plot_role="mark",
                value=str(observed[key]),
                display_value=f"{observed[key]:,}",
                numerator=str(observed[key]),
                denominator="",
                unit=unit,
                biological_unit=unit,
                model_contrast="descriptive frozen cohort census",
                effect_unit="count",
                evidence_status="descriptive",
                tested_universe="frozen public human cohort registry",
                **census_dependency,
                allowed_wording="descriptive resource scale with explicit analysis denominator",
                prohibited_wording="largest, broadest, most modalities, or prospective cohort claim",
                claim_text=f"The frozen resource census contains {observed[key]:,} {label.lower()}.",
                next_experiment="not applicable to a descriptive census",
                manuscript_included="true",
                plot_order=str(order),
                is_control="false",
            )
        )
    roles = (
        "Inherited regulatory evidence",
        "Established-state transcriptomics",
        "Physical context",
        "Prespecified functional challenge",
        "MASLD Gene Catalog",
    )
    for order, label in enumerate(roles, 1):
        rows.append(
            source_row(
                record_id=f"observability_role_{order}",
                group_id=f"observability_role_{order}",
                artifact_id="cohort_overview",
                claim_id=f"claim_observability_role_{order}",
                section_id="resource_interface",
                figure_id="Figure1",
                panel_id="1B",
                panel_title="Evidence observability and source firewall",
                panel_role="interface",
                label=label,
                category="distinct evidence role",
                number_role="categorical_role_indicator",
                plot_role="mark",
                value="1",
                display_value="preserved as a distinct role",
                numerator="",
                denominator="",
                unit="display-only categorical indicator, not a score",
                biological_unit="evidence role",
                model_contrast="observability interface",
                effect_unit="not an effect",
                evidence_status="descriptive",
                tested_universe="five predeclared evidence roles",
                **observability_dependency,
                allowed_wording="coverage-dependent complementary evidence role",
                prohibited_wording="modality total, convergence score, probability, or rank",
                claim_text=f"{label} remains a distinct evidence role under the source firewall.",
                next_experiment="measure the decisive missing assay for a named hypothesis",
                manuscript_included="true",
                plot_order=str(order),
                is_control="false",
            )
        )
    validate_source_rows(rows, "cohort_overview")
    return rows


def build_hotspot_rows(
    registry_rows: Sequence[Mapping[str, str]],
    figure_rows: Sequence[Mapping[str, str]],
    semantic_rows: Sequence[Mapping[str, str]],
    design_rows: Sequence[Mapping[str, str]],
    cohort_lodo_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, object]]:
    validate_hotspot_registry_projection(registry_rows, figure_rows)
    by_uid = {row["program_uid"]: row for row in figure_rows}
    semantic_by_uid = {row["program_uid"]: row for row in semantic_rows}
    if (
        len(by_uid) != 117
        or len(semantic_by_uid) != 117
        or set(by_uid) != set(semantic_by_uid)
    ):
        raise CoordinatorContractError(
            "Hotspot Figure2/semantic registries are not the same 117-program family"
        )
    selected = [
        row
        for row in semantic_rows
        if parse_bool(
            row["robust_display"], f"Hotspot {row['program_uid']} robust_display"
        )
    ]
    if len(selected) != 2:
        raise CoordinatorContractError(
            f"expected two robust Hotspot programs; found {len(selected)}"
        )
    output = []
    dependency = source_dependency(
        "reused_source",
        ("hotspot_module_registry_v1", "integrated_scrna_donors_7datasets"),
        ("hotspot_module_registry_v1", "integrated_scrna_donors_7datasets"),
        "The frozen module definitions, donor-stage evaluation, and LODO audit reuse the integrated single-cell source and pre-existing Hotspot registry.",
        "This is an internal donor-level refit and sensitivity analysis, not an external validation cohort.",
    )
    for order, semantic in enumerate(
        sorted(selected, key=lambda row: int(row["module"])), 1
    ):
        if semantic["adjudicated_state"] != "supported_internal_stage_association":
            raise CoordinatorContractError(
                "robust Hotspot row lacks supported semantic adjudication"
            )
        if parse_bool(
            semantic["tested_negative_authorized"],
            "Hotspot tested-negative authorization",
        ):
            raise CoordinatorContractError(
                "Hotspot robust presentation row authorizes tested-negative semantics"
            )
        row = by_uid[semantic["program_uid"]]
        shared = {
            "group_id": semantic["program_uid"],
            "artifact_id": "hotspot_programs",
            "claim_id": f"claim_{semantic['program_uid']}",
            "section_id": "established_state_transcriptomics",
            "figure_id": "Figure2",
            "panel_id": "2C",
            "panel_title": "Robust donor-level established-state programs",
            "panel_role": "transcriptomic_hero",
            "label": semantic["module_name"],
            "category": f"{semantic['cell_type']} Hotspot program",
            "biological_unit": "human biological donor",
            "model_contrast": (
                "program_score ~ stage_ordinal + dataset; Healthy=0, Steatosis=1, "
                "Steatohepatitis=2; cirrhosis excluded"
            ),
            "evidence_status": "state_associated",
            "tested_universe": "117 prespecified Hotspot modules; one BH family",
            **dependency,
            "allowed_wording": "robust cross-sectional donor-level stage association",
            "prohibited_wording": (
                "longitudinal progression, cell-intrinsic mechanism, subtype, or tested negative"
            ),
            "claim_text": (
                f"{semantic['module_name']} is robustly associated with cross-sectional "
                "stage order in the donor-level model."
            ),
            "next_experiment": "prospective donor-resolved perturbation or longitudinal validation",
            "plot_order": str(order),
            "is_control": "false",
        }
        output.append(
            source_row(
                **shared,
                record_id=f"hotspot_{semantic['program_uid']}_beta",
                number_role="stage_ordinal_beta",
                plot_role="mark",
                value=row["beta"],
                display_value=f"{float(row['beta']):.3f}",
                numerator="",
                denominator=row["n_donors"],
                unit="donor-level score beta per stage unit",
                effect_unit="donor-score beta per cross-sectional stage unit",
                p_value=row["pvalue"],
                q_value=row["qvalue"],
                manuscript_included="true",
            )
        )
        output.append(
            source_row(
                **shared,
                record_id=f"hotspot_{semantic['program_uid']}_se",
                number_role="stage_ordinal_se",
                plot_role="annotation",
                value=row["se"],
                display_value=f"{float(row['se']):.3f}",
                numerator="",
                denominator=row["n_donors"],
                unit="donor-level model standard error",
                effect_unit="sampling standard error of stage-ordinal beta",
                p_value="",
                q_value="",
                manuscript_included="false",
            )
        )
        output.append(
            source_row(
                **shared,
                record_id=f"hotspot_{semantic['program_uid']}_robust_display",
                number_role="robust_display_flag",
                plot_role="annotation",
                value="1",
                display_value="TRUE",
                numerator="1",
                denominator="117",
                unit="frozen semantic-adjudication selection flag",
                effect_unit="binary selection flag; not an effect",
                p_value="",
                q_value="",
                manuscript_included="false",
            )
        )
    selected_uids = {row["program_uid"] for row in selected}
    full_stage_datasets: dict[str, set[str]] = {}
    for uid in selected_uids:
        per_dataset: dict[str, set[str]] = {}
        program_design = [row for row in design_rows if row["program_uid"] == uid]
        if not program_design:
            raise CoordinatorContractError(f"Hotspot design audit omits {uid}")
        for row in program_design:
            if int(row["stage_dataset_n_donors"]) > 0 and row["stage_ordinal"] in {
                "0",
                "1",
                "2",
            }:
                per_dataset.setdefault(row["dataset"], set()).add(row["stage_ordinal"])
        full_stage_datasets[uid] = {
            dataset
            for dataset, stages in per_dataset.items()
            if stages == {"0", "1", "2"}
        }
    if any(datasets != {"GSE244832"} for datasets in full_stage_datasets.values()):
        raise CoordinatorContractError(
            "Hotspot stage-by-dataset design boundary changed: GSE244832 is no longer "
            "the sole cohort spanning Healthy/Steatosis/Steatohepatitis"
        )
    figure_direction = {row["program_uid"]: row["direction"] for row in figure_rows}
    gse244832_lodo = [
        row
        for row in cohort_lodo_rows
        if row["program_uid"] in selected_uids
        and row["analysis_type"] == "lodo"
        and row["held_out_dataset"] == "GSE244832"
    ]
    if len(gse244832_lodo) != len(selected_uids) or any(
        not parse_bool(row["estimable"], f"Hotspot {row['program_uid']} LODO estimable")
        for row in gse244832_lodo
    ):
        raise CoordinatorContractError(
            "Hotspot leave-GSE244832-out sensitivity is incomplete"
        )
    n_direction_agree = sum(
        row["direction"] == figure_direction[row["program_uid"]]
        for row in gse244832_lodo
    )
    if n_direction_agree != len(selected_uids):
        raise CoordinatorContractError(
            "a displayed Hotspot program changes direction when GSE244832 is held out"
        )
    output.append(
        source_row(
            record_id="hotspot_stage_dataset_lodo_boundary",
            group_id="hotspot_stage_dataset_lodo_boundary",
            artifact_id="hotspot_programs",
            claim_id="claim_hotspot_stage_dataset_lodo_boundary",
            section_id="established_state_transcriptomics",
            figure_id="Figure2",
            panel_id="2C",
            panel_title="Robust donor-level established-state programs",
            panel_role="transcriptomic_hero",
            label="Stage-by-dataset separability and leave-cohort-out boundary",
            category="design and sensitivity audit",
            number_role="leave_gse244832_direction_agreement_count",
            plot_role="annotation",
            value=str(n_direction_agree),
            display_value=(
                f"{n_direction_agree}/{len(selected_uids)} retain direction without "
                "GSE244832; only GSE244832 spans Healthy/Steatosis/Steatohepatitis"
            ),
            numerator=str(n_direction_agree),
            denominator=str(len(selected_uids)),
            unit="displayed programs retaining direction in leave-GSE244832-out donor models",
            biological_unit="human biological donor",
            model_contrast=(
                "leave GSE244832 out; program_score ~ stage_ordinal + dataset; "
                "Healthy=0, Steatosis=1, Steatohepatitis=2"
            ),
            effect_unit="direction-agreement count; not an independent effect estimate",
            p_value="",
            q_value="",
            evidence_status="descriptive",
            tested_universe="two prespecified robust displayed programs and the complete stage-by-dataset design audit",
            **dependency,
            allowed_wording="weak stage-cohort separability, possible cohort influence, and direction-preserving leave-one-dataset-out sensitivity",
            prohibited_wording="stage-cohort confounding eliminated, cohort-independent effect, external replication, or longitudinal progression",
            claim_text=(
                "Only GSE244832 spans all three primary stage levels, so possible cohort "
                "influence remains despite both displayed programs retaining direction "
                "when GSE244832 is held out."
            ),
            next_experiment="balanced donor recruitment with every stage represented within multiple independent cohorts",
            manuscript_included="true",
            plot_order="3",
            is_control="false",
        )
    )
    validate_source_rows(output, "hotspot")
    return output


def composition_slug(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    if not slug:
        raise CoordinatorContractError(
            f"composition cell type has no stable slug: {value!r}"
        )
    return slug


def build_composition_rows(
    result_rows: Sequence[Mapping[str, str]],
    testability_rows: Sequence[Mapping[str, str]],
    ready: Mapping[str, str],
) -> list[dict[str, object]]:
    """Adapt the frozen Plan20 sample-level refit without recomputing a call."""

    validate_composition_contract(result_rows, testability_rows, ready)
    by_result = {row["celltype"]: row for row in result_rows}
    output: list[dict[str, object]] = []
    tested_universe = (
        "22 declared deposited cell-type columns; BH across all 16 testable columns; "
        "6 structurally unavailable"
    )
    composition_dependency = source_dependency(
        "reused_source",
        ("canonical_bulk_human_cohorts_9", "music_deconvolution_persample_v1"),
        ("canonical_bulk_human_cohorts_9", "music_deconvolution_persample_v1"),
        "The composition hypotheses and Disease-minus-Control estimates reuse the frozen bulk cohort and deposited MuSiC deconvolution table.",
        "The refit corrects the biological unit and QC census but is not an independent composition replication.",
    )
    for order, testability in enumerate(
        sorted(testability_rows, key=lambda row: row["celltype"]), 1
    ):
        celltype = testability["celltype"]
        slug = composition_slug(celltype)
        shared = {
            "record_id": f"composition_{slug}",
            "group_id": f"composition_{slug}",
            "artifact_id": "sample_composition",
            "claim_id": f"claim_composition_{slug}",
            "section_id": "established_state_transcriptomics",
            "figure_id": "Figure2",
            "panel_id": "2B",
            "panel_title": "QC-passing sample-level cell composition",
            "panel_role": "transcriptomic_hero",
            "label": celltype,
            "plot_role": "mark",
            "tested_universe": tested_universe,
            **composition_dependency,
            "next_experiment": "independent sample-resolved tissue imaging validation",
            "manuscript_included": "false",
            "plot_order": str(order),
            "is_control": "false",
        }
        if celltype in by_result:
            result = by_result[celltype]
            effect = float(result["effect_arcsin_sqrt"])
            qvalue = float(result["padj"])
            output.append(
                source_row(
                    **shared,
                    category="sample-level cell composition",
                    number_role="composition_effect_arcsin_sqrt",
                    value=result["effect_arcsin_sqrt"],
                    display_value=f"{effect:.3f}",
                    numerator="",
                    denominator=result["n_analyzed"],
                    unit="arcsine-square-root proportion difference",
                    biological_unit=result["biological_unit"],
                    model_contrast=(
                        "Disease minus Control coefficient; " + result["model_formula"]
                    ),
                    effect_unit="arcsine-square-root proportion difference",
                    p_value=result["pvalue"],
                    q_value=result["padj"],
                    evidence_status=(
                        "state_associated" if qvalue < 0.05 else "indeterminate"
                    ),
                    allowed_wording=result["allowed_wording"],
                    prohibited_wording=result["prohibited_wording"],
                    claim_text=(
                        f"{celltype} has a QC-filtered sample-level cross-sectional "
                        f"composition association in the prespecified model."
                        if qvalue < 0.05
                        else f"{celltype} sample-level composition evidence is "
                        "indeterminate in the prespecified model."
                    ),
                )
            )
        else:
            output.append(
                source_row(
                    **shared,
                    category="sample-level cell composition testability",
                    number_role="composition_testability_indicator",
                    value="1",
                    display_value="untestable: structurally unavailable",
                    numerator="",
                    denominator="",
                    unit=(
                        "display-only categorical testability indicator; value 1 is not "
                        "an effect and not a zero-effect estimate"
                    ),
                    biological_unit=(
                        "QC-passing human liver sample; deposited assay column "
                        "structurally unavailable"
                    ),
                    model_contrast="not applicable: no nonmissing deposited values",
                    effect_unit="not applicable",
                    p_value="",
                    q_value="",
                    evidence_status="untestable",
                    allowed_wording="structurally unavailable deposited cell-type column",
                    prohibited_wording=(
                        "donor-level, zero effect, tested negative, longitudinal change, "
                        "lineage transition, or causal composition effect"
                    ),
                    claim_text=(
                        f"{celltype} is untestable because its deposited merged "
                        "composition column is structurally unavailable."
                    ),
                )
            )
    n_q_significant = sum(float(row["padj"]) < 0.05 for row in result_rows)
    output.append(
        source_row(
            record_id="composition_family_summary",
            group_id="composition_family_summary",
            artifact_id="sample_composition",
            claim_id="claim_composition_family_summary",
            section_id="established_state_transcriptomics",
            figure_id="Figure2",
            panel_id="2B",
            panel_title="QC-passing sample-level cell composition",
            panel_role="transcriptomic_hero",
            label="Complete QC-filtered composition family",
            category="sample-level cell composition family summary",
            number_role="bh_significant_testable_celltype_count",
            plot_role="annotation",
            value=str(n_q_significant),
            display_value=(f"{n_q_significant}/16 q<0.05; 6/22 unavailable; n=1,221"),
            numerator=str(n_q_significant),
            denominator="16",
            unit=(
                "BH-significant testable cell types; 22 declared columns and "
                "6 structurally unavailable"
            ),
            biological_unit="QC-passing human liver sample with MuSiC deconvolution",
            model_contrast="Disease minus Control coefficient; " + COMPOSITION_MODEL,
            effect_unit="count within the complete prespecified composition family",
            p_value="",
            q_value="",
            evidence_status="descriptive",
            tested_universe=tested_universe,
            **composition_dependency,
            allowed_wording="QC-filtered sample-level cross-sectional composition family",
            prohibited_wording=(
                "donor-level, longitudinal change, lineage transition, causal composition "
                "effect, or largest composition atlas"
            ),
            claim_text=(
                "The complete QC-filtered sample-level composition family retains all "
                "testable and structurally unavailable deposited cell-type columns."
            ),
            next_experiment="independent sample-resolved tissue imaging validation",
            manuscript_included="true",
            plot_order="23",
            is_control="false",
        )
    )
    validate_source_rows(output, "sample_composition")
    if len(output) != 23 or sum(row["plot_role"] == "mark" for row in output) != 22:
        raise CoordinatorContractError(
            "composition adapter does not preserve 22 columns plus one ledger-backed summary"
        )
    return output


def quantile_linear(values: Sequence[float], probability: float) -> float:
    if not values or probability < 0 or probability > 1:
        raise CoordinatorContractError("invalid continuous-loading quantile request")
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def build_nmf_continuous_supplement(
    k4_rows: Sequence[Mapping[str, str]],
    k6_rows: Sequence[Mapping[str, str]],
    k4_label_rows: Sequence[Mapping[str, str]],
    k6_label_rows: Sequence[Mapping[str, str]],
    metric_rows: Sequence[Mapping[str, str]],
    stability_rows: Sequence[Mapping[str, str]],
    k4_source_sha256: str,
    k6_source_sha256: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Extract continuous axes only and build supplement-only descriptive rows."""

    loading_sources = {4: k4_rows, 6: k6_rows}
    label_sources = {4: k4_label_rows, 6: k6_label_rows}
    source_hashes = {4: k4_source_sha256, 6: k6_source_sha256}
    sample_sets: dict[int, set[str]] = {}
    loading_output: list[dict[str, object]] = []
    evidence_output: list[dict[str, object]] = []
    axis_dependency = source_dependency(
        "reused_source",
        ("bulk_nmf_continuous_k4_seed42", "bulk_nmf_continuous_k6_seed42"),
        ("bulk_nmf_continuous_k4_seed42", "bulk_nmf_continuous_k6_seed42"),
        "The FigureS2 continuous-axis display is a descriptive extraction from the same frozen sample-by-program loading files.",
        "These are reused internal bulk data, not independent validation and not patient partitions.",
    )
    for k in (4, 6):
        rows = loading_sources[k]
        programs = tuple(f"P{index}" for index in range(1, k + 1))
        labels = {
            row["program_code"]: row["biological_label"] for row in label_sources[k]
        }
        if set(labels) != set(programs):
            raise CoordinatorContractError(f"k={k} NMF program-label universe drift")
        sample_ids = [row["sample_id"] for row in rows]
        if len(rows) != 1104 or len(set(sample_ids)) != len(sample_ids):
            raise CoordinatorContractError(
                f"k={k} continuous NMF table must contain 1,104 unique samples"
            )
        sample_sets[k] = set(sample_ids)
        for program_order, program in enumerate(programs, 1):
            values: list[float] = []
            for row in sorted(rows, key=lambda item: item["sample_id"]):
                value = parse_float(row[program], f"k={k} {program} continuous loading")
                if value < 0:
                    raise CoordinatorContractError(
                        f"k={k} {program} contains a negative NMF loading"
                    )
                values.append(value)
                loading_output.append(
                    {
                        "sample_id": row["sample_id"],
                        "k": str(k),
                        "program_code": program,
                        "program_label": labels[program],
                        "continuous_loading": f"{value:.17g}",
                        "source_sha256": source_hashes[k],
                    }
                )
            median = quantile_linear(values, 0.5)
            lower = quantile_linear(values, 0.25)
            upper = quantile_linear(values, 0.75)
            evidence_output.append(
                source_row(
                    record_id=f"nmf_k{k}_{program.lower()}_continuous",
                    group_id=f"nmf_k{k}_{program.lower()}",
                    artifact_id="nmf_continuous_supplement",
                    claim_id=f"claim_nmf_k{k}_{program.lower()}_continuous",
                    section_id="supplementary_nmf",
                    figure_id="FigureS2",
                    panel_id="S2A",
                    panel_title="Continuous k4/k6 NMF loading distributions",
                    panel_role="supplementary_continuous_programs",
                    label=f"k={k} {program}: {labels[program]}",
                    category=f"continuous k={k} NMF axis",
                    number_role="median_continuous_loading",
                    plot_role="mark",
                    value=f"{median:.17g}",
                    display_value=(
                        f"median={median:.3f}; IQR {lower:.3f}-{upper:.3f}; n=1,104"
                    ),
                    numerator="",
                    denominator="1104",
                    unit="nonnegative continuous NMF usage/loading",
                    biological_unit="QC-filtered human bulk RNA-seq sample in the frozen NMF fit",
                    model_contrast="descriptive sample-level continuous loading distribution",
                    effect_unit="median continuous loading; not an inferential effect or class assignment",
                    p_value="",
                    q_value="",
                    evidence_status="descriptive",
                    tested_universe=f"all {k} continuous axes across the same 1,104 samples",
                    **axis_dependency,
                    allowed_wording="continuous interpretive NMF axis with descriptive loading distribution",
                    prohibited_wording="subtype, hard cluster, patient class, reproducible stratification, transition, or biomarker",
                    claim_text=(
                        f"k={k} {program} is retained only as a continuous interpretive "
                        "axis across the frozen sample set."
                    ),
                    next_experiment="prospective outcome-linked validation before any patient-stratification use",
                    manuscript_included="false",
                    plot_order=str((k * 10) + program_order),
                    is_control="false",
                )
            )
    if sample_sets[4] != sample_sets[6]:
        raise CoordinatorContractError(
            "k=4 and k=6 continuous NMF sample universes differ"
        )
    if len(loading_output) != 11040:
        raise CoordinatorContractError("continuous NMF long table row count drift")

    metric_by_k = {row["k"]: row for row in metric_rows if row["k"] in {"4", "6"}}
    stability_by_k = {row["k"]: row for row in stability_rows if row["k"] in {"4", "6"}}
    if set(metric_by_k) != {"4", "6"} or set(stability_by_k) != {"4", "6"}:
        raise CoordinatorContractError(
            "canonical three-seed NMF stability summaries omit k=4 or k=6"
        )
    stability_dependency = source_dependency(
        "reused_source",
        ("bulk_nmf_three_seed_metrics_v2", "bulk_nmf_three_seed_stability_v2"),
        ("bulk_nmf_three_seed_metrics_v2", "bulk_nmf_three_seed_stability_v2"),
        "The displayed stability metrics are read directly from the canonical three-seed internal audit.",
        "Factor stability does not establish stable sample partitions or reproducible patient classes.",
    )
    stability_metrics = (
        ("mean_cophenetic", "mean_cophenetic", "mean cophenetic across three seeds"),
        ("mean_silhouette", "mean_silhouette", "mean silhouette across three seeds"),
        ("mean_cosine", "mean_matched_cosine", "mean Hungarian-matched factor cosine"),
        (
            "min_top50_jaccard",
            "min_top50_jaccard",
            "minimum matched top-50-gene Jaccard",
        ),
    )
    for k in (4, 6):
        combined = {**metric_by_k[str(k)], **stability_by_k[str(k)]}
        for metric_order, (source_field, number_role, label) in enumerate(
            stability_metrics, 1
        ):
            value = parse_float(combined[source_field], f"k={k} {source_field}")
            evidence_output.append(
                source_row(
                    record_id=f"nmf_k{k}_{number_role}",
                    group_id=f"nmf_k{k}_stability",
                    artifact_id="nmf_continuous_supplement",
                    claim_id=f"claim_nmf_k{k}_stability",
                    section_id="supplementary_nmf",
                    figure_id="FigureS2",
                    panel_id="S2B",
                    panel_title="Three-seed factor stability and hard-partition boundary",
                    panel_role="supplementary_continuous_programs",
                    label=f"k={k}: {label}",
                    category=f"k={k} three-seed factor audit",
                    number_role=number_role,
                    plot_role="mark",
                    value=f"{value:.17g}",
                    display_value=f"{value:.3f}",
                    numerator="",
                    denominator="3" if source_field.startswith("mean_") else "",
                    unit=label,
                    biological_unit="NMF factorization seed/factor comparison",
                    model_contrast="three-seed internal stability audit",
                    effect_unit="factor-stability statistic; not a sample-partition reproducibility estimate",
                    p_value="",
                    q_value="",
                    evidence_status="descriptive",
                    tested_universe="prespecified k=4 and k=6 three-seed NMF audits",
                    **stability_dependency,
                    allowed_wording="factor-axis stability with continuous interpretation only",
                    prohibited_wording="stable subtype, reproducible hard partition, patient class, clinical stratifier, or biomarker",
                    claim_text=(
                        f"The k={k} factor audit is shown with an explicit boundary that "
                        "factor stability does not imply stable patient partitions."
                    ),
                    next_experiment="prospective outcome-linked validation before any patient-stratification use",
                    manuscript_included="false",
                    plot_order=str((k * 10) + metric_order),
                    is_control="false",
                )
            )
    validate_source_rows(evidence_output, "continuous_nmf_supplement")
    return loading_output, evidence_output


def build_genetics_rows(
    phenotype_rows: Sequence[Mapping[str, str]],
    power_rows: Sequence[Mapping[str, str]],
    terminal_row: Mapping[str, str],
) -> list[dict[str, object]]:
    strata = Counter(row["phenotype_stratum"] for row in phenotype_rows)
    expected_strata = {
        "direct_masld_mash_diagnosis": 12,
        "mri_pdff_or_histologic_steatosis": 3,
        "alt_ast_or_ggt": 20,
    }
    if dict(strata) != expected_strata:
        raise CoordinatorContractError(
            f"genetic phenotype strata drift: {dict(strata)}"
        )
    if len(phenotype_rows) != 35:
        raise CoordinatorContractError(
            "genetic phenotype registry no longer contains 35 strata"
        )
    all_joint = [row for row in power_rows if row["universe"] == "all_joint_testable"]
    if len(all_joint) != 1 or all_joint[0]["status"] != "pass":
        raise CoordinatorContractError(
            "all-joint genetic/state interface is unavailable"
        )
    joint = all_joint[0]
    if joint["n_primary_genetic"] != "447" or joint["n_overlap"] != "34":
        raise CoordinatorContractError("joint genetic/state count drift")
    if not math.isclose(
        float(joint["fisher_odds_ratio"]), 0.889461813691, rel_tol=0, abs_tol=1e-12
    ) or not math.isclose(
        float(joint["fisher_p"]), 0.603934955413, rel_tol=0, abs_tol=1e-12
    ):
        raise CoordinatorContractError("joint genetic/state enrichment audit drift")
    if (
        terminal_row["context_rescue_authorized"] != "false"
        or terminal_row["negative_claim_authorized"] != "false"
    ):
        raise CoordinatorContractError(
            "coverage-limited genetics terminal authorizes a prohibited claim"
        )
    output: list[dict[str, object]] = []
    labels = {
        "direct_masld_mash_diagnosis": "Direct MASLD/MASH diagnosis",
        "mri_pdff_or_histologic_steatosis": "PDFF or histologic steatosis",
        "alt_ast_or_ggt": "ALT/AST/GGT proxy traits",
    }
    phenotype_dependency = source_dependency(
        "partially_dependent",
        ("genetic_evidence_source_gwas_set",),
        ("genetic_evidence_source_gwas_set", "phenotype_provenance_audit_v2"),
        "The provenance registry classifies the same GWAS sources used by the genetic analyses and adds a separate audit layer.",
        "The labels are audit-derived, but the displayed counts are not independent of the underlying GWAS source set.",
    )
    interface_dependency = source_dependency(
        "source_dependent",
        ("gen_joint_testable_universe_v2", "liver_eqtl_source_positive_universe"),
        ("gen_joint_testable_universe_v2", "liver_eqtl_source_positive_universe"),
        "The observability counts and interpretation are computed directly from the frozen joint and source-positive liver-eQTL universes.",
        "The complete source-tested non-eGene background is absent, so this panel cannot independently test genetic nulls or context rescue.",
    )
    bbj_rows = [row for row in phenotype_rows if row["study_name"].startswith("BBJ_")]
    if (
        len(bbj_rows) != 3
        or {row["ancestry"] for row in bbj_rows} != {"EAS"}
        or {row["eqtl_panel"] for row in bbj_rows}
        != {"Broadaway_liver_meta_N1183_EUR_hg19"}
        or {row["regulatory_ancestry_status"] for row in bbj_rows}
        != {"cross_ancestry_eqtl_limited"}
    ):
        raise CoordinatorContractError("BBJ ancestry/eQTL observability boundary drift")
    for order, key in enumerate(expected_strata, 1):
        output.append(
            source_row(
                record_id=f"phenotype_{key}",
                group_id=f"phenotype_{key}",
                artifact_id="genetics_evidence",
                claim_id=f"claim_phenotype_{key}",
                section_id="genetics_context",
                figure_id="Figure3",
                panel_id="3A",
                panel_title="Phenotype provenance",
                panel_role="genetics_boundary",
                label=labels[key],
                category="GWAS phenotype provenance",
                number_role="gwas_stratum_count",
                plot_role="mark",
                value=str(expected_strata[key]),
                display_value=str(expected_strata[key]),
                numerator=str(expected_strata[key]),
                denominator="35",
                unit="prespecified GWAS strata",
                biological_unit="GWAS stratum",
                model_contrast="descriptive phenotype registry",
                effect_unit="count",
                evidence_status="descriptive",
                tested_universe="35 prespecified GWAS strata",
                **phenotype_dependency,
                allowed_wording="phenotype-provenance-separated inherited regulatory evidence",
                prohibited_wording="all traits are direct MASLD susceptibility or independent GWAS",
                claim_text=f"The genetic portfolio includes {expected_strata[key]} {labels[key].lower()} strata.",
                next_experiment="ancestry-matched replication with phenotype-specific ascertainment",
                manuscript_included="true",
                plot_order=str(order),
                is_control="false",
            )
        )
    output.append(
        source_row(
            record_id="phenotype_bbj_eas_eur_eqtl_boundary",
            group_id="phenotype_bbj_eas_eur_eqtl_boundary",
            artifact_id="genetics_evidence",
            claim_id="claim_phenotype_bbj_eas_eur_eqtl_boundary",
            section_id="genetics_context",
            figure_id="Figure3",
            panel_id="3A",
            panel_title="Phenotype provenance",
            panel_role="genetics_boundary",
            label="BBJ regulatory-ancestry boundary",
            category="cross-ancestry eQTL observability",
            number_role="bbj_eas_gwas_with_eur_liver_eqtl_count",
            plot_role="annotation",
            value="3",
            display_value="3 BBJ EAS liver-enzyme GWAS paired to EUR liver eQTL",
            numerator="3",
            denominator="3",
            unit="BBJ GWAS strata with cross-ancestry regulatory evidence",
            biological_unit="GWAS stratum",
            model_contrast="BBJ EAS GWAS joined to Broadaway EUR bulk-liver eQTL",
            effect_unit="descriptive ancestry-panel pairing count; not a causal effect",
            p_value="",
            q_value="",
            evidence_status="descriptive",
            tested_universe="three prespecified BBJ ALT/AST/GGT proxy-trait strata",
            **phenotype_dependency,
            allowed_wording="EAS GWAS with ancestry-matched EAS LD but EUR bulk-liver eQTL due to unavailable ancestry-matched EAS liver eQTL",
            prohibited_wording="ancestry-matched regulatory replication, cross-ancestry causality, direct MASLD phenotype, or source-negative conclusion",
            claim_text=(
                "The three BBJ EAS liver-enzyme GWAS use ancestry-matched EAS LD but "
                "EUR Broadaway bulk-liver eQTL because no ancestry-matched EAS liver "
                "eQTL panel is available."
            ),
            next_experiment="repeat regulatory colocalization with an adequately powered EAS liver eQTL panel",
            manuscript_included="true",
            plot_order="4",
            is_control="false",
        )
    )
    interface_rows = (
        ("joint_genetic", "Jointly testable primary genetic genes", 447, "gene"),
        (
            "joint_overlap",
            "Primary genetic genes also established-state associated",
            34,
            "gene",
        ),
        (
            "source_egenes",
            "Source-defined significant liver eGenes",
            6564,
            "source-defined unique Ensembl gene",
        ),
        (
            "tested_background",
            "Complete source-tested non-eGene background",
            0,
            "display-only missingness placeholder",
        ),
    )
    for order, (key, label, value, unit) in enumerate(interface_rows, 1):
        missing = key == "tested_background"
        output.append(
            source_row(
                record_id=f"genetic_boundary_{key}",
                group_id=f"genetic_boundary_{key}",
                artifact_id="genetics_evidence",
                claim_id=f"claim_genetic_boundary_{key}",
                section_id="genetics_context",
                figure_id="Figure3",
                panel_id="3B",
                panel_title="Static interface and eQTL observability boundary",
                panel_role="genetics_boundary",
                label=label,
                category="genetic observability",
                number_role="coverage_gate" if missing else "descriptive_count",
                plot_role="mark",
                value=str(value),
                display_value="not deposited" if missing else f"{value:,}",
                numerator="" if missing else str(value),
                denominator="447" if key == "joint_overlap" else "",
                unit=unit,
                biological_unit=unit,
                model_contrast="coverage-limited terminal genetic audit",
                effect_unit="not an observed zero" if missing else "count",
                evidence_status="untestable" if missing else "descriptive",
                tested_universe="terminal GEN audit; source-positive eGenes only",
                **interface_dependency,
                allowed_wording="static genetically anchored interface with explicit eQTL observability limit",
                prohibited_wording="context rescue, powered genetic null, source negative, or interchange of identifier units",
                claim_text=(
                    "The complete source-tested non-eGene background was not deposited, so context rescue and source-negative claims are not authorized."
                    if missing
                    else f"The terminal genetics audit records {value:,} {label.lower()}."
                ),
                next_experiment="obtain the complete source-tested eQTL universe with expression and local variant-density covariates",
                manuscript_included="true",
                plot_order=str(order),
                is_control="false",
            )
        )
    output.append(
        source_row(
            record_id="genetic_boundary_overlap_enrichment",
            group_id="genetic_boundary_overlap_enrichment",
            artifact_id="genetics_evidence",
            claim_id="claim_genetic_boundary_overlap_enrichment",
            section_id="genetics_context",
            figure_id="Figure3",
            panel_id="3B",
            panel_title="Static interface and eQTL observability boundary",
            panel_role="genetics_boundary",
            label="Genetic/state overlap enrichment audit",
            category="descriptive genetic-state interface",
            number_role="fisher_overlap_odds_ratio",
            plot_role="annotation",
            value=joint["fisher_odds_ratio"],
            display_value=(
                f"34/447 overlap; OR={float(joint['fisher_odds_ratio']):.2f}; "
                f"Fisher p={float(joint['fisher_p']):.3f}"
            ),
            numerator=joint["n_overlap"],
            denominator=joint["n_primary_genetic"],
            unit="Fisher overlap odds ratio in the 14,931-gene joint-testable universe",
            biological_unit="joint-testable Ensembl-mapped gene",
            model_contrast="primary genetic evidence versus established-state association",
            effect_unit="Fisher exact overlap odds ratio",
            p_value=joint["fisher_p"],
            q_value="",
            evidence_status="indeterminate",
            tested_universe="14,931 all-joint-testable genes",
            **interface_dependency,
            allowed_wording="descriptive non-enriched overlap in the frozen joint-testable universe",
            prohibited_wording="replicated cross-ancestry causality, causal-versus-reactive separation, powered genetic null, or biological equivalence",
            claim_text=(
                "The 34-of-447 genetic/state overlap is descriptive and non-enriched in "
                "the frozen joint universe (OR 0.89; Fisher p 0.604), not evidence of "
                "replicated cross-ancestry causality."
            ),
            next_experiment="ancestry-matched regulatory mapping and mechanistic perturbation for named loci",
            manuscript_included="true",
            plot_order="5",
            is_control="false",
        )
    )
    validate_source_rows(output, "genetics")
    return output


def direction_value(value: str, context: str) -> int:
    mapping = {"positive": 1, "negative": -1, "zero": 0, "none": 0, "": 0}
    if value not in mapping:
        raise CoordinatorContractError(f"{context} has unknown direction: {value!r}")
    return mapping[value]


def spatial_source_dependency(row: Mapping[str, str]) -> dict[str, str]:
    dataset_source = "dataset_" + re.sub(
        r"[^A-Za-z0-9_.:-]+", "_", row["dataset"]
    ).strip("_")
    dependence = row["source_dependence"]
    if dependence == "independent":
        return source_dependency(
            "independent",
            ("hotspot_registry_v2",),
            (dataset_source,),
            "The program identity is frozen upstream; the assay effect is estimated only in the named external dataset.",
            "The evaluation dataset did not contribute to program discovery or registry selection; independence is asserted under the declared accession/source IDs because cross-alias equivalence is not independently audited.",
        )
    if dependence == "source_dependent":
        return source_dependency(
            "source_dependent",
            ("hotspot_registry_v2", dataset_source),
            (dataset_source,),
            "The named source dataset contributed to an upstream context or program-support role and is reused for this assay-native display.",
            "The row is retained as source-dependent context evidence and is not called independent replication.",
        )
    raise CoordinatorContractError(
        f"unsupported Plan13 source-dependence class: {dependence!r}"
    )


def build_spatial_rows(
    matrix_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, object]]:
    if len(matrix_rows) != 18:
        raise CoordinatorContractError(
            "spatial adapter requires the complete 18-row Plan13 matrix"
        )
    robust_rows = [row for row in matrix_rows if row["evidence_state"] == "robust"]
    if len(robust_rows) != 2 or len({row["program_uid"] for row in robust_rows}) != 1:
        raise CoordinatorContractError(
            "Plan13 must retain two robust dataset-by-program rows for one unique program; "
            "the row count is not a count of distinct externally robust programs"
        )
    panels = {
        "native": ("4A", "Assay-native physical context"),
        "yak": ("4B", "Zonation-adjusted lipid challenge"),
        "atac": ("4A", "Assay-native physical context"),
        "gate": ("4B", "Zonation-adjusted lipid challenge"),
    }
    output: list[dict[str, object]] = []
    for order, row in enumerate(matrix_rows, 1):
        state = row["evidence_state"]
        if state == "tested_negative":
            raise CoordinatorContractError(
                "semantic-v2 spatial matrix contains tested_negative"
            )
        if row["dataset"] in {"GSE192741", "Vu_et_al_2025"}:
            lane = "native"
        elif row["dataset"] == "Yakubovsky_2026":
            lane = "yak"
        elif row["dataset"] in {"GSE281367", "GSE244832"}:
            lane = "atac"
        else:
            lane = "gate"
        panel_id, title = panels[lane]
        estimate = str(row["estimate"]).strip()
        has_effect = bool(estimate)
        value = estimate if has_effect else "1"
        display = f"{float(estimate):.3f}" if has_effect else state
        if lane == "native" and str(row["std_error"]).strip():
            raise CoordinatorContractError(
                "native Moran matched-null dispersion is labeled as sampling SE"
            )
        if (
            row["biological_unit_resolution"] == "unresolved"
            and str(row["n_biological"]).strip()
        ):
            raise CoordinatorContractError(
                "unresolved spatial unit reports a biological n"
            )
        status = {
            "robust": "context_supported",
            "indeterminate": "indeterminate",
            "untestable": "untestable",
            "not_applicable": "not_applicable",
            "skipped": "skipped",
        }.get(state)
        if status is None:
            raise CoordinatorContractError(
                f"unsupported Plan13 evidence state: {state}"
            )
        unit = (
            row["effect_unit"]
            if has_effect
            else "display-only gate indicator, not an inferential effect"
        )
        bio = row["biological_unit"]
        if row["biological_unit_resolution"] == "unresolved":
            bio = "donor unknown; technical/reporting unit only"
        allowed_wording = (
            "assay-native physical context or explicit observability boundary"
        )
        prohibited_wording = (
            "cross-assay score, broad validation, cell-intrinsic mechanism, or donor n "
            "inferred from technical units"
        )
        claim_text = (
            f"{row['dataset']} places the prespecified program in assay-native physical "
            f"context with state {state}."
        )
        if row["dataset"] == "Vu_et_al_2025":
            if (
                row["source_dependence"] != "source_dependent"
                or row["biological_unit_resolution"] != "unresolved"
                or row["n_technical"] != "10"
                or str(row["n_biological"]).strip()
            ):
                raise CoordinatorContractError(
                    "Vu evidence must remain source-dependent with ten technical sections "
                    "and unresolved donor identity"
                )
            allowed_wording = (
                "source-dependent reproducible section-level spatial organization with "
                "unresolved donor identity"
            )
            prohibited_wording = (
                "independent spatial replication, donor-level external replication, "
                "population generalizability, donor count, or broad validation"
            )
            claim_text = (
                "Vu supports source-dependent section-level spatial organization across "
                "ten technical sections; unresolved donor identity prevents donor-level "
                "external-replication or population-generalizability claims."
            )
        elif row["dataset"] == "GSE192741":
            if (
                row["source_dependence"] != "independent"
                or row["biological_unit_resolution"] != "resolved"
                or row["n_biological"] != "4"
            ):
                raise CoordinatorContractError(
                    "GSE192741 must remain the independent four-donor spatial dataset"
                )
            allowed_wording = "independent donor-resolved spatial-context replication within GSE192741"
            prohibited_wording = (
                "population generalizability, universal spatial validation, cross-assay "
                "score, or cell-intrinsic mechanism"
            )
        mark_group = f"spatial_{row['dataset']}_{row['program_uid']}"
        dependency = spatial_source_dependency(row)
        output.append(
            source_row(
                record_id=f"spatial_{order}",
                group_id=mark_group,
                artifact_id="spatial_evidence",
                claim_id=f"claim_spatial_{order}",
                section_id="physical_context",
                figure_id="Figure4",
                panel_id=panel_id,
                panel_title=title,
                panel_role="prespecified_external_challenge",
                label=f"{row['dataset']}: {row['program_label']}",
                category=f"{row['assay']} ({state})",
                number_role={
                    "native": "spatial_native_effect",
                    "yak": "lipid_context_effect",
                    "atac": "chromatin_context_effect",
                    "gate": "assay_gate_status",
                }[lane],
                plot_role="mark",
                value=value,
                display_value=display,
                numerator="",
                denominator=row["n_biological"]
                if row["biological_unit_resolution"] == "resolved"
                else "",
                unit=unit,
                biological_unit=bio,
                model_contrast=row["contrast_or_exposure"],
                effect_unit=unit,
                p_value=row["pvalue"] if has_effect else "",
                q_value=row["padj"] if has_effect else "",
                evidence_status=status,
                tested_universe="two prespecified hepatocyte programs within each assay/dataset family",
                **dependency,
                allowed_wording=allowed_wording,
                prohibited_wording=prohibited_wording,
                claim_text=claim_text,
                next_experiment="matched donor-resolved tissue and perturbation in the source-relevant context",
                manuscript_included="true",
                plot_order=str(order),
                is_control="false",
            )
        )
        output.append(
            source_row(
                record_id=f"spatial_{order}_technical_units",
                group_id=mark_group,
                artifact_id="spatial_evidence",
                claim_id=f"claim_spatial_{order}",
                section_id="physical_context",
                figure_id="Figure4",
                panel_id=panel_id,
                panel_title=title,
                panel_role="prespecified_external_challenge",
                label=f"{row['dataset']}: {row['program_label']}",
                category=f"{row['assay']} technical-unit ledger",
                number_role="technical_unit_count",
                plot_role="annotation",
                value=row["n_technical"],
                display_value=f"{row['n_technical']} {row['technical_unit']}",
                numerator=row["n_technical"],
                denominator="",
                unit=f"technical/reporting unit count ({row['technical_unit']}); not a biological denominator",
                biological_unit="technical-unit ledger entry; biological n is recorded separately",
                model_contrast="descriptive assay-unit census",
                effect_unit="count; not an inferential effect",
                p_value="",
                q_value="",
                evidence_status="descriptive",
                tested_universe="complete technical/reporting units declared by the Plan13 assay adapter",
                **dependency,
                allowed_wording="technical-unit count displayed separately from biological replication",
                prohibited_wording="donor count, biological denominator, independent replicate count, or effect estimate",
                claim_text=f"{row['dataset']} records {row['n_technical']} {row['technical_unit']} technical/reporting units separately from biological n.",
                next_experiment="resolve donor identity for any assay whose biological unit remains unresolved",
                manuscript_included="false",
                plot_order=str(order),
                is_control="false",
            )
        )
        uncertainty_rows = (
            (
                "assay_native_standard_error",
                row["std_error"],
                "sampling standard error in the source assay's native effect unit",
            ),
            (
                "matched_null_standard_deviation",
                row["matched_null_sd"],
                "matched-gene null standard deviation; not sampling SE",
            ),
            (
                "assay_native_interval_lower",
                row["interval_low"],
                f"{row['interval_type']} lower bound in the source assay's native unit",
            ),
            (
                "assay_native_interval_upper",
                row["interval_high"],
                f"{row['interval_type']} upper bound in the source assay's native unit",
            ),
        )
        for uncertainty_order, (
            number_role,
            uncertainty,
            uncertainty_unit,
        ) in enumerate(uncertainty_rows, 1):
            if not str(uncertainty).strip():
                continue
            output.append(
                source_row(
                    record_id=f"spatial_{order}_{number_role}",
                    group_id=mark_group,
                    artifact_id="spatial_evidence",
                    claim_id=f"claim_spatial_{order}",
                    section_id="physical_context",
                    figure_id="Figure4",
                    panel_id=panel_id,
                    panel_title=title,
                    panel_role="prespecified_external_challenge",
                    label=f"{row['dataset']}: {row['program_label']}",
                    category=f"{row['assay']} ({state})",
                    number_role=number_role,
                    plot_role="annotation",
                    value=str(uncertainty),
                    display_value=f"{float(uncertainty):.3f}",
                    numerator="",
                    denominator=row["n_biological"]
                    if row["biological_unit_resolution"] == "resolved"
                    else "",
                    unit=uncertainty_unit,
                    biological_unit=bio,
                    model_contrast=row["contrast_or_exposure"],
                    effect_unit=uncertainty_unit,
                    p_value="",
                    q_value="",
                    evidence_status=status,
                    tested_universe="two prespecified hepatocyte programs within each assay/dataset family",
                    **dependency,
                    allowed_wording=(
                        allowed_wording
                        + "; assay-native uncertainty shown separately from the effect estimate"
                    ),
                    prohibited_wording=(
                        prohibited_wording
                        + "; sampling error substituted for matched-null dispersion"
                    ),
                    claim_text=claim_text,
                    next_experiment="matched donor-resolved tissue and perturbation in the source-relevant context",
                    manuscript_included="false",
                    plot_order=str(order),
                    is_control="false",
                )
            )
        if row["dataset"] == "Yakubovsky_2026" and str(
            row["legacy_program_id"]
        ).endswith(":8"):
            directions = (
                (
                    "descriptive_median_slope_direction",
                    row["descriptive_effect_direction"],
                    "descriptive donor-median slope direction",
                ),
                (
                    "inferential_signed_stouffer_direction",
                    row["inferential_test_direction"],
                    "signed-Stouffer inferential direction",
                ),
            )
            for offset, (number_role, direction, label) in enumerate(directions, 1):
                output.append(
                    source_row(
                        record_id=f"yak_module8_{number_role}",
                        group_id="yak_module8",
                        artifact_id="spatial_evidence",
                        claim_id=f"claim_yak_module8_{number_role}",
                        section_id="physical_context",
                        figure_id="Figure4",
                        panel_id="4B",
                        panel_title=panels["yak"][1],
                        panel_role="prespecified_external_challenge",
                        label=f"Stromal ECM: {label}",
                        category="direction/heterogeneity audit",
                        number_role=number_role,
                        plot_role="annotation",
                        value=str(direction_value(direction, number_role)),
                        display_value=direction,
                        numerator="",
                        denominator=row["n_biological"],
                        unit="signed direction indicator",
                        biological_unit="public lipid donor",
                        model_contrast=row["contrast_or_exposure"],
                        effect_unit="direction only; not an effect magnitude",
                        p_value=row["pvalue"]
                        if number_role.startswith("inferential")
                        else "",
                        q_value=row["padj"]
                        if number_role.startswith("inferential")
                        else "",
                        evidence_status="indeterminate",
                        tested_universe="two prespecified hepatocyte programs",
                        **dependency,
                        allowed_wording="descriptive and inferential directions shown separately with donor heterogeneity",
                        prohibited_wording="direction agreement, robust lipid gradient, hepatocyte-specific effect, or continuous lipid exposure",
                        claim_text="Yakubovsky module 8 has discordant descriptive-median and signed-Stouffer directions with donor heterogeneity.",
                        next_experiment="donor-resolved continuous lipid measurement with cell-type-resolved spatial expression",
                        manuscript_included="true",
                        plot_order=str(order * 10 + offset),
                        is_control="false",
                    )
                )
    validate_source_rows(output, "spatial")
    return output


def build_myojin_rows(
    class_rows: Sequence[Mapping[str, str]],
    program_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    myojin_dependency = source_dependency(
        "independent",
        ("genetics_context_v2", "hotspot_registry_v2"),
        ("myojin_hlf_palmitate_screen_v1",),
        "The target/program hypotheses were frozen from genetics and Hotspot work before the held-out HLF screen outcomes were inspected.",
        "The Myojin HLF screen did not contribute to target discovery, class assignment, expected direction, or program selection; independence is asserted under the declared source IDs because cross-alias equivalence is not independently audited.",
    )
    for order, row in enumerate(class_rows, 1):
        is_omnibus = row["statistic_type"] == "partial_F"
        is_pairwise = row["statistic_type"] == "absolute_OLS_coefficient"
        if is_omnibus == is_pairwise:
            raise CoordinatorContractError(
                f"Myojin class row has unsupported statistic type: {row['statistic_type']}"
            )
        estimate = row["statistic"] if is_omnibus else row["standardized_effect"]
        if not estimate or (is_omnibus and row["standardized_effect"]):
            raise CoordinatorContractError(
                f"Myojin class row has no displayable statistic: {row['contrast']}"
            )
        if is_pairwise and (not row["standardized_effect"] or not row["BH_q"]):
            raise CoordinatorContractError(
                f"Myojin pairwise row lacks standardized effect/BH q: {row['contrast']}"
            )
        output.append(
            source_row(
                record_id=f"myojin_class_{order}",
                group_id=f"myojin_class_{order}",
                artifact_id="myojin_supplement",
                claim_id=f"claim_myojin_class_{order}",
                section_id="supplementary_myojin",
                figure_id="FigureS1",
                panel_id="S1A",
                panel_title="Complete Myojin evidence-class stress test",
                panel_role="supplementary_functional",
                label=f"{row['analysis_variant']}: {row['contrast']}",
                category="HLF palmitate-specific screen",
                number_role=(
                    "class_omnibus_partial_f"
                    if is_omnibus
                    else "class_pairwise_standardized_effect"
                ),
                plot_role="annotation" if is_omnibus else "mark",
                value=estimate,
                display_value=f"{float(estimate):.3f}",
                numerator=row["n_class"],
                denominator=row["n_total"],
                unit=(
                    "partial F omnibus statistic; annotation only and not comparable to "
                    "standardized pairwise effects"
                    if is_omnibus
                    else "standardized evidence-class contrast effect"
                ),
                biological_unit="screened gene in one HLF cell line",
                model_contrast=row["contrast"],
                effect_unit=(
                    "partial F omnibus statistic"
                    if is_omnibus
                    else "standardized evidence-class contrast effect"
                ),
                p_value=row["permutation_p"],
                q_value="" if is_omnibus else row["BH_q"],
                evidence_status="indeterminate",
                tested_universe="complete frozen Myojin evidence-class analysis table",
                **myojin_dependency,
                allowed_wording="assay-specific non-support or nonconfirmatory functional challenge",
                prohibited_wording="precise null, tested negative, broad functional validation, or in-vivo mechanism",
                claim_text="The complete prespecified Myojin evidence-class branch was nonconfirmatory in HLF cells.",
                next_experiment="repeat in independent primary-like human hepatocyte models with orthogonal phenotypes",
                manuscript_included="false",
                plot_order=str(order),
                is_control="false",
            )
        )
    for order, row in enumerate(program_rows, 1):
        testable = parse_bool(row["testable"], f"Myojin {row['program_uid']} testable")
        estimate = row["signed_effect"] if testable else "1"
        output.append(
            source_row(
                record_id=f"myojin_program_{order}",
                group_id=f"myojin_program_{order}",
                artifact_id="myojin_supplement",
                claim_id=f"claim_myojin_program_{order}",
                section_id="supplementary_myojin",
                figure_id="FigureS1",
                panel_id="S1B",
                panel_title="Complete Myojin program stress test",
                panel_role="supplementary_functional",
                label=f"{row['analysis_variant']} / {row['weight_mode']} / {row['program_uid']}",
                category="HLF palmitate-specific program challenge",
                number_role="program_signed_effect"
                if testable
                else "program_testability_indicator",
                plot_role="mark",
                value=estimate,
                display_value=f"{float(estimate):.3f}" if testable else "untestable",
                numerator=row["n_eligible_genes"],
                denominator="",
                unit=(
                    "matched-null signed effect"
                    if testable
                    else "display-only categorical testability indicator; value 1 is not an effect and not a zero-effect estimate"
                ),
                biological_unit="screened gene in one HLF cell line",
                model_contrast=f"{row['analysis_variant']} {row['weight_mode']}",
                effect_unit="matched-null signed effect"
                if testable
                else "not applicable",
                p_value=row["empirical_p"] if testable else "",
                q_value=row["BH_q"] if testable else "",
                evidence_status="indeterminate" if testable else "untestable",
                tested_universe="complete frozen Myojin program analysis table",
                **myojin_dependency,
                allowed_wording="assay-specific non-support or nonconfirmatory program challenge",
                prohibited_wording="expression reversal, tested negative, broad validation, or gene-level expansion",
                claim_text="The complete prespecified Myojin program branch was nonconfirmatory or untestable in HLF cells.",
                next_experiment="repeat the program perturbation in independent primary-like hepatocyte models",
                manuscript_included="false",
                plot_order=str(order),
                is_control="false",
            )
        )
    validate_source_rows(output, "myojin")
    class_output = [row for row in output if row["panel_id"] == "S1A"]
    if (
        len(class_output) != 28
        or sum(row["number_role"] == "class_omnibus_partial_f" for row in class_output)
        != 7
        or sum(
            row["number_role"] == "class_pairwise_standardized_effect"
            for row in class_output
        )
        != 21
        or any(
            row["plot_role"] != "annotation"
            for row in class_output
            if row["number_role"] == "class_omnibus_partial_f"
        )
        or any(
            row["plot_role"] != "mark"
            for row in class_output
            if row["number_role"] == "class_pairwise_standardized_effect"
        )
    ):
        raise CoordinatorContractError(
            "Myojin class presentation must separate 7 omnibus annotations from "
            "21 standardized pairwise-effect marks"
        )
    return output


def validate_passport_negative_semantics(
    evidence_rows: Sequence[Mapping[str, object]],
) -> None:
    for row in evidence_rows:
        state = str(row.get("call_state", ""))
        rule = str(row.get("negative_call_rule_id", "") or "").strip()
        passed = row.get("negative_call_passed")
        boundary = row.get("negative_decision_boundary")
        margin = row.get("negative_margin")
        if state == "tested_negative":
            if str(row.get("testability_state", "")) != "testable":
                raise CoordinatorContractError(
                    f"Gene Catalog tested_negative row is not testable: {row.get('evidence_result_id')}"
                )
            boundary_value = parse_float(
                boundary, f"Gene Catalog {row.get('evidence_result_id')} negative boundary"
            )
            margin_value = parse_float(
                margin, f"Gene Catalog {row.get('evidence_result_id')} negative margin"
            )
            if (
                not rule
                or passed is not True
                or boundary_value is None
                or margin_value is None
            ):
                raise CoordinatorContractError(
                    f"Gene Catalog tested_negative row lacks a complete passing rule: {row.get('evidence_result_id')}"
                )
            if margin_value < 0:
                raise CoordinatorContractError(
                    f"Gene Catalog tested_negative row has a negative adequacy margin: {row.get('evidence_result_id')}"
                )
        elif (
            rule
            or boundary is not None
            or margin is not None
            or passed not in {None, False}
        ):
            raise CoordinatorContractError(
                f"Gene Catalog non-negative row carries reserved negative metadata: {row.get('evidence_result_id')}"
            )


def nullable_number_text(value: object) -> str:
    if value is None:
        return ""
    parsed = parse_float(value, "nullable presentation number")
    return f"{parsed:.12g}"


def stable_record_token(value: object) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    if not token:
        raise CoordinatorContractError(
            f"cannot form a stable record token from {value!r}"
        )
    return token


def passport_graph_indices(
    source_nodes: Sequence[Mapping[str, object]],
    source_edges: Sequence[Mapping[str, object]],
) -> tuple[dict[str, Mapping[str, object]], dict[str, Mapping[str, object]]]:
    nodes = {str(row["source_node_id"]): row for row in source_nodes}
    edges = {str(row["source_edge_id"]): row for row in source_edges}
    if (
        len(nodes) != len(source_nodes)
        or len(edges) != len(source_edges)
        or not nodes
        or not edges
    ):
        raise CoordinatorContractError(
            "Plan50 source graph has blank/duplicate/empty IDs"
        )
    for edge_id, edge in edges.items():
        if (
            str(edge["from_node_id"]) not in nodes
            or str(edge["to_node_id"]) not in nodes
        ):
            raise CoordinatorContractError(
                f"Plan50 source edge {edge_id} references an unknown source node"
            )
    return nodes, edges


def passport_evidence_graph_refs(
    result: Mapping[str, object],
    nodes: Mapping[str, Mapping[str, object]],
    edges: Mapping[str, Mapping[str, object]],
) -> set[str]:
    result_id = str(result["evidence_result_id"])
    source_node = str(result.get("source_node_id") or "")
    if source_node not in nodes:
        raise CoordinatorContractError(
            f"Gene Catalog entry {result_id} lacks a valid source_node_id"
        )
    matches = [
        (edge_id, edge)
        for edge_id, edge in edges.items()
        if str(edge.get("edge_type")) == "tested_by"
        and str(edge.get("evidence_result_id")) == result_id
    ]
    if not matches:
        raise CoordinatorContractError(
            f"Gene Catalog entry {result_id} has no exact tested_by source-graph edge"
        )
    refs = {source_node}
    for edge_id, edge in matches:
        refs.update((edge_id, str(edge["from_node_id"]), str(edge["to_node_id"])))
    return refs


def passport_coverage_graph_refs(
    coverage: Mapping[str, object],
    nodes: Mapping[str, Mapping[str, object]],
    edges: Mapping[str, Mapping[str, object]],
) -> set[str]:
    source_release = str(coverage.get("source_release_id") or "")
    if not source_release:
        raise CoordinatorContractError(
            f"Gene Catalog coverage {coverage.get('coverage_result_id')} lacks source_release_id"
        )
    matched_nodes = {
        node_id
        for node_id, node in nodes.items()
        if str(node.get("source_release_id")) == source_release
    }
    if not matched_nodes:
        raise CoordinatorContractError(
            f"Gene Catalog coverage source release is absent from the source graph: {source_release}"
        )
    refs = set(matched_nodes)
    for edge_id, edge in edges.items():
        if (
            str(edge["from_node_id"]) in matched_nodes
            or str(edge["to_node_id"]) in matched_nodes
        ):
            refs.update((edge_id, str(edge["from_node_id"]), str(edge["to_node_id"])))
    return refs


def passport_vignette_dependency(
    graph_refs: Iterable[str],
    source_provenance_state: str,
) -> dict[str, str]:
    evaluation = tuple(sorted(set(graph_refs)))
    if not evaluation:
        raise CoordinatorContractError(
            "Gene Catalog example has no exact source-graph references"
        )
    return source_dependency(
        "source_dependent",
        ("figure5_vignette_selection_posthoc", *evaluation),
        evaluation,
        (
            "The illustrative Figure5 gene was fixed after project results were known; "
            f"the underlying source row retains provenance_state={source_provenance_state}, "
            "but the vignette reuses the exact named Plan50 source nodes/edges."
        ),
        "Figure 5 contains illustrative Gene Catalog entries, not held-out target discovery or independent validation.",
    )


def build_passport_rows(
    gene_rows: Sequence[Mapping[str, object]],
    evidence_rows: Sequence[Mapping[str, object]],
    coverage_rows: Sequence[Mapping[str, object]],
    experiment_rows: Sequence[Mapping[str, object]],
    program_context_rows: Sequence[Mapping[str, object]],
    source_node_rows: Sequence[Mapping[str, object]],
    source_edge_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    for table_name, table_rows in (
        ("gene index", gene_rows),
        ("gene evidence", evidence_rows),
        ("gene coverage", coverage_rows),
        ("next experiment", experiment_rows),
        ("program context", program_context_rows),
    ):
        if any(
            str(row.get("analysis_release_id", "")) != PLAN50_ANALYSIS_RELEASE_ID
            for row in table_rows
        ):
            raise CoordinatorContractError(
                f"Plan50 {table_name} contains an unexpected nested analysis release ID"
            )
    validate_passport_negative_semantics(evidence_rows)
    validate_passport_negative_semantics(program_context_rows)
    source_nodes, source_edges = passport_graph_indices(
        source_node_rows, source_edge_rows
    )
    if any(
        row.get("gene_call_expansion_authorized") is not False
        for row in program_context_rows
    ):
        raise CoordinatorContractError(
            "Plan50 program context authorizes program-to-gene call expansion"
        )
    by_symbol: dict[str, list[Mapping[str, object]]] = {}
    for row in gene_rows:
        by_symbol.setdefault(str(row["symbol"]), []).append(row)
    evidence_by_passport: dict[str, list[Mapping[str, object]]] = {}
    for row in evidence_rows:
        evidence_by_passport.setdefault(str(row["passport_id"]), []).append(row)
    coverage_by_passport: dict[str, list[Mapping[str, object]]] = {}
    for row in coverage_rows:
        coverage_by_passport.setdefault(str(row["passport_id"]), []).append(row)
    experiment_by_passport: dict[str, Mapping[str, object]] = {}
    for row in experiment_rows:
        passport = str(row["passport_id"])
        if passport in experiment_by_passport:
            raise CoordinatorContractError(
                f"duplicate Gene Catalog next experiment: {passport}"
            )
        experiment_by_passport[passport] = row
    output: list[dict[str, object]] = []
    hero_ensembl_ids: set[str] = set()
    for order, symbol in enumerate(FIXED_HERO_GENES, 1):
        matches = by_symbol.get(symbol, [])
        if len(matches) != 1:
            raise CoordinatorContractError(
                f"Gene Catalog example symbol {symbol} maps to {len(matches)} rows"
            )
        gene = matches[0]
        passport = str(gene["passport_id"])
        ensembl = str(gene["ensembl_id"])
        if not re.fullmatch(r"ENSG\d+(?:\.\d+)?", ensembl):
            raise CoordinatorContractError(
                f"Gene Catalog example {symbol} lacks Ensembl-first identity"
            )
        if ensembl in hero_ensembl_ids:
            raise CoordinatorContractError(
                f"Gene Catalog example Ensembl identity is duplicated: {ensembl}"
            )
        hero_ensembl_ids.add(ensembl)
        if gene.get("symbol_collision") is not False:
            raise CoordinatorContractError(
                f"Gene Catalog example {symbol} has an unresolved symbol collision"
            )
        evidence = evidence_by_passport.get(passport, [])
        if not evidence:
            raise CoordinatorContractError(
                f"Gene Catalog example {symbol} has no gene-level evidence rows"
            )
        if any(str(row["ensembl_id"]) != ensembl for row in evidence):
            raise CoordinatorContractError(
                f"Gene Catalog example {symbol} evidence identity drift"
            )
        coverage = coverage_by_passport.get(passport, [])
        if not coverage or any(str(row["ensembl_id"]) != ensembl for row in coverage):
            raise CoordinatorContractError(
                f"Gene Catalog example {symbol} coverage identity drift"
            )
        experiment = experiment_by_passport.get(passport)
        if experiment is None or str(experiment["ensembl_id"]) != ensembl:
            raise CoordinatorContractError(
                f"Gene Catalog example {symbol} lacks an Ensembl-linked next experiment"
            )
        provenance = sorted({str(row["provenance_state"]) for row in evidence})
        provenance_text = provenance[0] if len(provenance) == 1 else "mixed"
        hero_graph_refs: set[str] = set()
        for result in evidence:
            hero_graph_refs.update(
                passport_evidence_graph_refs(result, source_nodes, source_edges)
            )
        shared_dependency = passport_vignette_dependency(
            hero_graph_refs,
            provenance_text,
        )
        shared = {
            "artifact_id": "passport_evidence",
            "section_id": "evidence_passports",
            "figure_id": "Figure5",
            "panel_role": "passport_boundary",
            "category": str(gene["primary_evidence_class"]),
            "plot_role": "mark",
            "value": "1",
            "numerator": "",
            "denominator": "",
            "unit": "fixed display-only Gene Catalog example indicator; not a score or rank",
            "biological_unit": "Ensembl-identified gene",
            "effect_unit": "not an effect",
            "p_value": "",
            "q_value": "",
            "evidence_status": "passport_example",
            "tested_universe": "four predeclared manual-acceptance hero genes",
            **shared_dependency,
            "prohibited_wording": "best target, leaderboard, score, rank, probability, treatment recommendation, or program-to-gene call expansion",
            "manuscript_included": "true",
            "plot_order": str(order),
            "is_control": "false",
        }
        output.append(
            source_row(
                **shared,
                record_id=f"passport_{symbol}_state",
                group_id=f"passport_{symbol}",
                claim_id=f"claim_passport_{symbol}_state",
                panel_id="5A",
                panel_title="Evidence role, testability, and source dependence",
                label=f"{symbol} ({ensembl})",
                number_role="categorical_passport_example",
                display_value=str(gene["primary_evidence_class"]),
                model_contrast="source-preserving MASLD Gene Catalog display",
                allowed_wording="representative source-preserving Gene Catalog entry with visible testability and provenance",
                claim_text=f"{symbol} is shown as a fixed representative Gene Catalog entry without a combined score or rank.",
                next_experiment=str(experiment["primary_readout"]),
            )
        )
        for evidence_order, result in enumerate(
            sorted(evidence, key=lambda item: str(item["evidence_result_id"])), 1
        ):
            call_state = str(result["call_state"])
            testability = str(result["testability_state"])
            result_id = str(result["evidence_result_id"])
            estimate = nullable_number_text(result.get("estimate"))
            value = estimate or "1"
            negative_rule = str(result.get("negative_call_rule_id") or "")
            negative_criterion = ""
            allowed = (
                "source-preserving assay call with visible testability and provenance"
            )
            if call_state == "tested_negative":
                negative_criterion = (
                    f"{negative_rule}; boundary={nullable_number_text(result.get('negative_decision_boundary'))}; "
                    f"margin={nullable_number_text(result.get('negative_margin'))}; passed=true"
                )
                allowed = "tested negative only under the frozen adequate-negative rule and margin"
            effect_unit = str(result.get("effect_unit") or "").strip()
            if not estimate:
                effect_unit = (
                    "display-only call-state indicator; not an inferential effect"
                )
            result_dependency = passport_vignette_dependency(
                passport_evidence_graph_refs(result, source_nodes, source_edges),
                str(result["provenance_state"]),
            )
            output.append(
                source_row(
                    record_id=f"passport_{symbol}_evidence_{stable_record_token(result_id)}",
                    group_id=f"passport_{symbol}",
                    artifact_id="passport_evidence",
                    claim_id=f"claim_passport_{symbol}_evidence_{evidence_order}",
                    section_id="evidence_passports",
                    figure_id="Figure5",
                    panel_id="5A",
                    panel_title="Evidence role, testability, and source dependence",
                    panel_role="passport_boundary",
                    label=f"{symbol}: {result['evidence_domain']} / {result['dataset_id']}",
                    category=f"{result['assay']} assay call",
                    number_role=f"assay_call_{evidence_order}",
                    plot_role="annotation",
                    value=value,
                    display_value=f"{call_state} / {testability}",
                    numerator="",
                    denominator=""
                    if result.get("n_biological_units") is None
                    else str(result["n_biological_units"]),
                    unit=effect_unit,
                    biological_unit=str(result["biological_unit"]),
                    model_contrast=str(result["contrast_or_exposure"]),
                    effect_unit=effect_unit,
                    p_value=nullable_number_text(result.get("p_value")),
                    q_value=nullable_number_text(result.get("q_value")),
                    evidence_status=call_state,
                    negative_adequacy_criterion=negative_criterion,
                    tested_universe=f"{result['dataset_id']} / {result['assay']} accepted source result",
                    **result_dependency,
                    allowed_wording=allowed,
                    prohibited_wording="cross-assay score, hidden rank, treatment recommendation, or stronger call than the source",
                    claim_text=f"{symbol} retains the source assay call {call_state} with {testability} testability.",
                    next_experiment=str(experiment["primary_readout"]),
                    manuscript_included="false",
                    plot_order=str(order * 100 + evidence_order),
                    is_control="false",
                )
            )
        for coverage_order, coverage_row in enumerate(
            sorted(coverage, key=lambda item: str(item["coverage_result_id"])), 1
        ):
            coverage_id = str(coverage_row["coverage_result_id"])
            call_state = str(coverage_row["call_state"])
            testability = str(coverage_row["testability_state"])
            coverage_dependency = passport_vignette_dependency(
                passport_coverage_graph_refs(coverage_row, source_nodes, source_edges),
                "source_dependent",
            )
            output.append(
                source_row(
                    record_id=f"passport_{symbol}_coverage_{stable_record_token(coverage_id)}",
                    group_id=f"passport_{symbol}",
                    artifact_id="passport_evidence",
                    claim_id=f"claim_passport_{symbol}_coverage_{coverage_order}",
                    section_id="evidence_passports",
                    figure_id="Figure5",
                    panel_id="5A",
                    panel_title="Evidence role, testability, and source dependence",
                    panel_role="passport_boundary",
                    label=f"{symbol}: {coverage_row['evidence_domain']} coverage",
                    category=f"{coverage_row['assay']} coverage",
                    number_role=f"coverage_state_{coverage_order}",
                    plot_role="annotation",
                    value="1",
                    display_value=f"{call_state} / {testability}",
                    numerator="",
                    denominator=(
                        ""
                        if coverage_row.get("n_biological_units") is None
                        else str(coverage_row["n_biological_units"])
                    ),
                    unit="display-only assay coverage state; not a score",
                    biological_unit="assay coverage record",
                    model_contrast=f"{coverage_row['dataset_id']} assay coverage",
                    effect_unit="not an inferential effect",
                    p_value="",
                    q_value="",
                    evidence_status="descriptive",
                    tested_universe=str(coverage_row["coverage_denominator"]),
                    **coverage_dependency,
                    allowed_wording="visible assay coverage and testability boundary",
                    prohibited_wording="missing coverage as negative evidence, cross-assay score, or rank",
                    claim_text=f"{symbol} retains the deposited {testability} coverage state for {coverage_row['assay']}.",
                    next_experiment=str(experiment["primary_readout"]),
                    manuscript_included="false",
                    plot_order=str(order * 1000 + coverage_order),
                    is_control="false",
                )
            )
        output.append(
            source_row(
                **shared,
                record_id=f"passport_{symbol}_experiment",
                group_id=f"passport_{symbol}_experiment",
                claim_id=f"claim_passport_{symbol}_experiment",
                panel_id="5B",
                panel_title="Next discriminating experiments",
                label=f"{symbol}: {experiment['biological_model']}",
                number_role="categorical_next_experiment",
                display_value=str(experiment["primary_readout"]),
                model_contrast="falsifiable next-experiment rulebook",
                allowed_wording="next discriminating experiment with explicit falsifying outcome",
                claim_text=f"The {symbol} Gene Catalog entry names a discriminating experiment and a falsifying outcome.",
                next_experiment=(
                    f"{experiment['perturbation']}; readout: {experiment['primary_readout']}; "
                    f"falsifier: {experiment['falsifying_outcome']}"
                ),
            )
        )
    validate_source_rows(output, "passport")
    return output


def locked_release_blueprint() -> dict[str, object]:
    figures = [
        (
            "Figure1",
            (
                ("1A", "Cohort and analysis-unit overview", "interface"),
                ("1B", "Evidence observability and source firewall", "interface"),
            ),
        ),
        (
            "Figure2",
            (
                (
                    "2A",
                    "Cross-sectional stage-associated DEG counts",
                    "transcriptomic_hero",
                ),
                (
                    "2B",
                    "QC-passing sample-level cell composition",
                    "transcriptomic_hero",
                ),
                (
                    "2C",
                    "Robust donor-level established-state programs",
                    "transcriptomic_hero",
                ),
            ),
        ),
        (
            "Figure3",
            (
                ("3A", "Phenotype provenance", "genetics_boundary"),
                (
                    "3B",
                    "Static interface and eQTL observability boundary",
                    "genetics_boundary",
                ),
            ),
        ),
        (
            "Figure4",
            (
                (
                    "4A",
                    "Assay-native physical context",
                    "prespecified_external_challenge",
                ),
                (
                    "4B",
                    "Zonation-adjusted lipid challenge",
                    "prespecified_external_challenge",
                ),
            ),
        ),
        (
            "Figure5",
            (
                (
                    "5A",
                    "Evidence role, testability, and source dependence",
                    "passport_boundary",
                ),
                ("5B", "Next discriminating experiments", "passport_boundary"),
            ),
        ),
    ]
    return {
        "contract_version": 1,
        "candidate_id": CANDIDATE_ID,
        "journal_branch": "Cell Genomics Resource",
        "figure_count": 5,
        "myojin_role": "supplement_only",
        "results_sections": [
            "resource_interface",
            "established_state_transcriptomics",
            "genetics_context",
            "physical_context",
            "evidence_passports",
        ],
        "figures": [
            {
                "figure_id": figure_id,
                "title": EXPECTED_FIGURE_TITLES[figure_id],
                "panels": [
                    {"panel_id": panel, "title": title, "role": role}
                    for panel, title, role in panels
                ],
            }
            for figure_id, panels in figures
        ],
        "supplementary_figures": [
            {
                "figure_id": "FigureS1",
                "title": "Prespecified Myojin HLF functional stress test",
                "panels": [
                    {
                        "panel_id": "S1A",
                        "title": "Complete Myojin evidence-class stress test",
                        "role": "supplementary_functional",
                    },
                    {
                        "panel_id": "S1B",
                        "title": "Complete Myojin program stress test",
                        "role": "supplementary_functional",
                    },
                ],
            },
            {
                "figure_id": "FigureS2",
                "title": "Continuous NMF axes and factor stability",
                "panels": [
                    {
                        "panel_id": "S2A",
                        "title": "Continuous k4/k6 NMF loading distributions",
                        "role": "supplementary_continuous_programs",
                    },
                    {
                        "panel_id": "S2B",
                        "title": "Three-seed factor stability and hard-partition boundary",
                        "role": "supplementary_continuous_programs",
                    },
                ],
            },
        ],
        "sources": [
            {
                "source_key": "BASE:cohort_overview",
                "snapshot_path": "inputs/BASE/cohort_overview.tsv",
                "kind": "evidence_rows",
            },
            {
                "source_key": f"BASE:{FIBROSIS_PRIMARY_ARTIFACT_ID}",
                "snapshot_path": fibrosis_snapshot_relpath(FIBROSIS_PRIMARY_RELATIVE),
                "kind": "fibrosis_transition_raw",
            },
            {
                "source_key": "BASE:sample_composition",
                "snapshot_path": "inputs/BASE/sample_composition.tsv",
                "kind": "evidence_rows",
            },
            {
                "source_key": "BASE:hotspot_programs",
                "snapshot_path": "inputs/BASE/hotspot_programs.tsv",
                "kind": "evidence_rows",
            },
            {
                "source_key": "BASE:genetics_evidence",
                "snapshot_path": "inputs/BASE/genetics_evidence.tsv",
                "kind": "evidence_rows",
            },
            {
                "source_key": "BASE:spatial_evidence",
                "snapshot_path": "inputs/BASE/spatial_evidence.tsv",
                "kind": "evidence_rows",
            },
            {
                "source_key": "BASE:passport_evidence",
                "snapshot_path": "inputs/BASE/passport_evidence.tsv",
                "kind": "evidence_rows",
            },
            {
                "source_key": "BASE:myojin_supplement",
                "snapshot_path": "inputs/BASE/myojin_supplement.tsv",
                "kind": "evidence_rows",
            },
            {
                "source_key": "PLAN20:nmf_continuous_supplement",
                "snapshot_path": "inputs/HS-V2/nmf_continuous_supplement.tsv",
                "kind": "evidence_rows",
            },
        ],
        "directional_decisions": {
            "figure4_title": "Physical context, assay observability, and prespecified external challenges",
            "cross_sectional_not_longitudinal": True,
            "myojin_supplement_only": True,
            "nmf_continuous_axes_supplement_only": True,
            "nmf_hard_partition_claims_prohibited": True,
            "no_cross_assay_score": True,
            "no_program_to_gene_call_expansion": True,
            "no_tested_negative_without_frozen_rule": True,
            "canonical_promotion_authorized": False,
        },
    }


def closure_rows(
    project_root: Path,
    handoffs: Mapping[str, Handoff],
    coordinator: str,
    signing_date: str,
) -> list[dict[str, object]]:
    decisions = {
        "PLAN13": (
            "accepted_main",
            "SP_INT_SEMANTIC_V2",
            "SEMANTIC_V2_READY",
            "true",
            "true",
            "figure_4_physical_context",
            "assay-native physical context and explicit observability boundaries",
            "broad validation; cross-assay score; cell-intrinsic validation",
            "",
        ),
        "PLAN20": (
            "accepted_main",
            "HS_V2_READY_AND_SEMANTIC",
            "READY",
            "true",
            "true",
            "figure_2_programs",
            "robust cross-sectional donor-level established-state programs under the post-freeze semantic adjudication",
            "longitudinal progression; legacy tested_negative boolean; subtype",
            "",
        ),
        "PLAN30": (
            "accepted_main",
            "GEN_TERMINAL",
            "coverage_limited_terminal",
            "true",
            "true",
            "figure_3_genetics",
            "phenotype-separated static genetic/state interface with positive-only eGene annotation",
            "context rescue; powered genetic null; source negative",
            "",
        ),
        "PLAN40": (
            "accepted_supplement",
            "HLF_FINAL",
            "supplement_only_nonconfirmatory",
            "false",
            "true",
            "supplementary_functional",
            "assay-specific non-support/nonconfirmatory HLF challenge",
            "tested negative; precise null; broad functional validation",
            "main-figure gate false",
        ),
        "PLAN50": (
            "accepted_main",
            "PASS06",
            "passport_release_validated",
            "true",
            "true",
            "figure_5_passports",
            "source-preserving Gene Catalog entries and falsifiable next experiments",
            "score; rank; treatment recommendation; program-to-gene expansion",
            "",
        ),
    }
    rows = []
    for workstream in WORKSTREAM_ORDER:
        state, gate_id, verdict, main, supp, role, allowed, prohibited, reason = (
            decisions[workstream]
        )
        handoff = handoffs[workstream]
        rows.append(
            {
                "workstream_id": workstream,
                "terminal_state": state,
                "gate_id": gate_id,
                "gate_verdict": verdict,
                "artifact_manifest": project_relative(project_root, handoff.manifest),
                "manifest_sha256": handoff.manifest_sha256,
                "include_main": main,
                "include_supplement": supp,
                "figure_role": role,
                "allowed_wording": allowed,
                "prohibited_wording": prohibited,
                "exclusion_reason": reason,
                "adjudicator": coordinator,
                "date": signing_date,
            }
        )
    return rows


def coordinator_signature(
    artifact_path: Path,
    artifact_sha256: str,
    coordinator: str,
    signed_at_utc: str,
    decision_register_id: str,
    artifact_role: str,
) -> dict[str, object]:
    return {
        "attestation_version": COORDINATOR_SIGNATURE_VERSION,
        "candidate_id": CANDIDATE_ID,
        "artifact_role": artifact_role,
        "artifact_path": artifact_path.as_posix(),
        "artifact_sha256": artifact_sha256,
        "signed_by": coordinator,
        "signed_at_utc": signed_at_utc,
        "decision_register_id": decision_register_id,
        "authority_document": "docs/ROADMAP.md",
        "canonical_promotion_authorized": False,
    }


def protected_scope_rows(project_root: Path) -> list[dict[str, str]]:
    project = project_root.resolve()
    release_parent = project / "RNA-seq/results/manuscript_release"
    rows = []
    for root in sorted(
        path
        for path in release_parent.iterdir()
        if path.is_dir() and not path.is_symlink() and path.name != "candidates"
    ):
        rows.append(
            {
                "scope_id": f"manuscript_{root.name}",
                "protection_class": "named_manuscript_release",
                "root_path": project_relative(project, root),
            }
        )
    for scope_id, relative in (
        ("program_v1_results", "Analysis/Multimodal_Program_Projection/results"),
        ("program_v1_config", "Analysis/Multimodal_Program_Projection/config"),
        ("program_v1_figure5", "figures/main/fig5_molecular_context"),
    ):
        root = project / relative
        if not root.is_dir() or root.is_symlink():
            raise CoordinatorContractError(
                f"protected v1 root is missing or unsafe: {root}"
            )
        rows.append(
            {
                "scope_id": scope_id,
                "protection_class": "program_v1",
                "root_path": relative,
            }
        )
    if not any(row["protection_class"] == "named_manuscript_release" for row in rows):
        raise CoordinatorContractError(
            "no pre-existing named manuscript release was found"
        )
    return rows


def protected_baseline_rows(
    project_root: Path, scopes: Sequence[Mapping[str, str]]
) -> list[dict[str, object]]:
    project = project_root.resolve()
    rows: list[dict[str, object]] = []
    for scope in scopes:
        root = project / scope["root_path"]
        files = []
        for path in sorted(root.rglob("*")):
            if path.is_symlink():
                raise CoordinatorContractError(
                    f"protected scope contains symlink: {path}"
                )
            if path.is_file():
                files.append(path)
        if not files:
            raise CoordinatorContractError(f"protected scope is empty: {root}")
        for path in files:
            rows.append(
                {
                    "scope_id": scope["scope_id"],
                    "relative_path": path.relative_to(root).as_posix(),
                    "sha256": sha256_file(path),
                    "bytes": path.stat().st_size,
                }
            )
    return rows


def selection_row(
    project_root: Path,
    source: Path,
    artifact_id: str,
    snapshot_relpath: str,
    artifact_role: str,
    coordinator: str,
    signing_date: str,
    decision_register_id: str,
    allowed: str,
    prohibited: str,
) -> dict[str, object]:
    return {
        "selection_id": decision_register_id,
        "coordinator": coordinator,
        "date": signing_date,
        "artifact_id": artifact_id,
        "source_path": project_relative(project_root, source),
        "source_sha256": sha256_file(source),
        "source_bytes": source.stat().st_size,
        "snapshot_relpath": snapshot_relpath,
        "artifact_role": artifact_role,
        "allowed_wording": allowed,
        "prohibited_wording": prohibited,
    }

#!/usr/bin/env python3
"""Shared, outcome-blind contracts for the Yakubovsky Plan 11 pipeline.

This module contains no import-time filesystem reads or writes.  In particular,
importing it cannot open the Yakubovsky source objects or their lipid labels.
The production entry points require the Plan 10 source gate and the independently
sealed Plan 20 registry before any candidate product is written.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import struct
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import h5py
import numpy as np
import pandas as pd
from patsy import dmatrix
from scipy import sparse, stats
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


RELEASE_ID = "program-context-v2-candidate-2026-08-07"
DATASET = "Yakubovsky_2026"
ALLOWED_SOURCE_GATE = "pass_ordinal_lipid"
ALLOWED_LIPID_CLASSES = ("lipid_zone", "non_lipid_zone")
GENE_AXIS_HASH_ALGORITHM = "sha256_len64le_utf8_after_matlab_zero_drop_v1"
GENE_AXIS_ZERO_POLICY = "drop_all_zero_code_units_matching_decode_matlab_char"

PRODUCTION_CODE_FILENAMES = (
    "00_freeze_plan11.py",
    "01_build_adapter.py",
    "02_score_programs.py",
    "03_fit_binary_lipid.py",
    "04_validate_bundle.py",
    "05_build_plan13_adapter.py",
    "06_check_stage.py",
    "07_verify_gene_axis.py",
    "yakubovsky_common.py",
    "run_adapter.sbatch",
    "run_models.sbatch",
    "run_validate.sbatch",
    "submit_plan11.sh",
)

SPATIAL_BLOCK_AUDIT_COLUMNS = (
    "donor",
    "section",
    "island_id",
    "block_id",
    "block_uid",
    "tile_x",
    "tile_y",
    "n_spots",
    "median_nearest_distance",
    "graph_max_edge_distance",
    "block_width",
    "n_small_block_merges_in_island",
    "small_block_merge_trace",
    "excluded_small_island",
)

ADAPTER_STAGE_ARTIFACTS = (
    "yakubovsky_human_adapter.h5ad",
    "sample_manifest.tsv",
    "gene_mapping_audit.tsv",
    "design_audit.tsv",
    "source_manifest.tsv",
    "adapter_execution_manifest.tsv",
)
SCORING_STAGE_ARTIFACTS = (
    "program_testability.tsv",
    "per_sample_program_scores.tsv",
    "per_donor_program_scores.tsv",
    "zonation_reference.tsv",
    "scoring_audit.tsv",
    "scoring_execution_manifest.tsv",
)
MODEL_STAGE_ARTIFACTS = (
    "spatial_block_audit.tsv",
    "model_sample_manifest.tsv",
    "donor_program_effects.tsv",
    "spatial_bootstrap_summary.tsv",
    "donor_combination_audit.tsv",
    "program_effects.tsv",
    "sensitivity.tsv",
    "multiplicity_manifest.tsv",
    "model_design_reference.tsv",
    "bootstrap_reference_draws.tsv",
    "model_design_audit.tsv",
    "model_execution_manifest.tsv",
)
TERMINAL_STAGE_ARTIFACTS = (
    "validation_report.tsv",
    "execution_manifest.tsv",
    "candidate_release_manifest.tsv",
    "gate_status.tsv",
    "READY",
)
FREEZE_STAGE_ARTIFACTS = (
    "gene_axis_verification.tsv",
    "analysis_freeze_manifest.tsv",
    "analysis_specification.tsv",
    "registry_manifest.tsv",
    "freeze_execution_manifest.tsv",
)

STAGE_ARTIFACTS = {
    "freeze": FREEZE_STAGE_ARTIFACTS,
    "adapter": ADAPTER_STAGE_ARTIFACTS,
    "scoring": SCORING_STAGE_ARTIFACTS,
    "model": MODEL_STAGE_ARTIFACTS,
    "terminal": TERMINAL_STAGE_ARTIFACTS,
}
STAGE_ORDER = tuple(STAGE_ARTIFACTS)

EXPECTED_REGISTRY_SHA256 = (
    "b134b5f8e46a716e857af22271aa54e1294432b9f9ab735559d7f603cdd29fe7"
)
EXPECTED_MEMBERSHIP_SHA256 = (
    "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b"
)
EXPECTED_VALIDATION_SHA256 = (
    "e291f7060af704e3accbcaa3e1c27551a73c41bfb3b818c623d4a4888b5aecf5"
)
EXPECTED_RELEASE_MANIFEST_SHA256 = (
    "6e3701ba48341fba493a7cf737745ec0640f1625e00b6ab290f4f3a201b92101"
)
EXPECTED_EXTERNAL_MODULES = (8, 20)

MIN_PROGRAM_GENES = 8
MIN_RETAINED_L1 = 0.20
MIN_DONOR_SPOTS = 100
MIN_FINAL_BLOCKS = 8
MIN_BLOCK_SPOTS = 5
GRAPH_DISTANCE_MULTIPLIER = 2.5
BLOCK_WIDTH_MULTIPLIER = 5.0
PRIMARY_SPLINE_DF = 4
SPLINE_SENSITIVITY_DF = (3, 5)
N_BOOTSTRAP = 4_999
MAX_FAILED_BOOTSTRAP_FRACTION = 0.05
BASE_SEED = 20260807

PRIMARY_MODEL = (
    "program_score ~ lipid_zone + ns(zonation,df=4) + "
    "center(log1p(background_corrected_spot_sum)) + center(detected_genes)"
)
ROBUST_RULE = (
    "complete_assay_testable_BH_family; primary_BH_q<0.05; signed_Stouffer_"
    "expected_direction; median_slope_same_direction; ceil(0.75*n_source_gate_"
    "donors)_observed_slope_direction; equal_weight_leave_top_df3_df5_signed_"
    "Stouffer_expected_direction; all_leave_one_donor_out_signed_Stouffer_"
    "expected_direction"
)


class ContractError(RuntimeError):
    """A frozen input, source gate, or output contract was violated."""


def project_root() -> Path:
    return Path(
        os.environ.get(
            "MASLD_PROJECT_ROOT",
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
        )
    ).resolve()


def default_paths(base: Path | None = None) -> dict[str, Path]:
    base = (base or project_root()).resolve()
    candidate = (
        base
        / "Analysis/Spatial/candidates"
        / RELEASE_ID
        / "yakubovsky2026"
    )
    acquisition = (
        base
        / "Analysis/Spatial/candidates"
        / RELEASE_ID
        / "acquisition/yakubovsky2026"
    )
    data = base / "Analysis/Spatial/data/public_expansion_2026/yakubovsky2026"
    hotspot = (
        base
        / "Analysis/Multimodal_Program_Projection/candidates"
        / RELEASE_ID
        / "hotspot"
    )
    return {
        "base": base,
        "candidate": candidate,
        "gene_axis_verification": candidate / "gene_axis_verification.tsv",
        "acquisition": acquisition,
        "data": data,
        "source_repo": data / "Human-liver",
        "source_importer": data
        / "Human-liver/Matlab_functions/import_visium_data_funcion_for_github.m",
        "source_background": data
        / "Human-liver/Matlab_functions/compute_background_counts_for_github.m",
        "hotspot": hotspot,
        "registry": hotspot / "program_registry_v2.tsv",
        "membership": hotspot / "program_membership_v2.tsv",
        "external": hotspot / "external_test_programs.tsv",
        "tested_universe": hotspot / "tested_universe.tsv",
        "hotspot_validation": hotspot / "validation_status.tsv",
        "hotspot_release_manifest": hotspot / "release_manifest.tsv",
        "hotspot_ready": hotspot / "READY",
        "v1_preservation": hotspot / "v1_preservation.tsv",
        "source_gate": acquisition / "yakubovsky_source_gate.tsv",
        "source_join": acquisition / "yakubovsky_join_audit.tsv",
        "source_files": acquisition / "source_files.tsv",
        "source_git": acquisition / "github_source_manifest.tsv",
        "v_mat": data / "v.mat",
        "gencode": base / "data/gencode_v49_gene_metadata.tsv.gz",
        "gene_axis_manifest": base
        / "Analysis/Spatial/scripts/yakubovsky2026/fixtures/"
        / "yakubovsky_gene_axis_equivalence.tsv",
    }


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(block_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def deterministic_seed(*parts: object) -> int:
    token = "\x1f".join(str(part) for part in parts)
    return int(hashlib.sha256(token.encode()).hexdigest()[:8], 16)


def atomic_write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        delete=False,
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_write_frame(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        frame.to_csv(handle, sep="\t", index=False, na_rep="")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def validate_analysis_freeze(output: Path) -> pd.DataFrame:
    """Rederive every pre-outcome file hash, including all producer code."""
    path = output / "analysis_freeze_manifest.tsv"
    if not path.is_file():
        raise ContractError(f"Missing pre-outcome analysis freeze: {path}")
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(
        frame,
        ["record_type", "name", "value", "sha256", "lipid_outcome_read"],
        "Plan 11 analysis freeze",
    )
    files = frame[frame["record_type"] == "file_input"].copy()
    if files["name"].duplicated().any():
        raise ContractError("Plan 11 analysis freeze contains duplicate file roles")
    if any(bool_value(value) for value in files["lipid_outcome_read"]):
        raise ContractError("A pre-outcome frozen file row claims lipid outcomes were read")
    frozen_by_name = files.set_index("name")
    script_dir = Path(__file__).resolve().parent
    expected_producers = {
        f"producer_code::{filename}": (script_dir / filename).resolve()
        for filename in PRODUCTION_CODE_FILENAMES
    }
    missing_roles = sorted(set(expected_producers).difference(frozen_by_name.index))
    if missing_roles:
        raise ContractError(
            "Pre-outcome seal lacks producer(s): " + ", ".join(missing_roles)
        )
    for role, expected_path in expected_producers.items():
        observed_path = Path(frozen_by_name.loc[role, "value"]).resolve()
        if observed_path != expected_path:
            raise ContractError(
                f"Frozen producer path drift ({role}): {observed_path} != {expected_path}"
            )
    for row in files.itertuples(index=False):
        frozen_path = Path(row.value)
        if not frozen_path.is_file():
            raise ContractError(f"Frozen file disappeared ({row.name}): {frozen_path}")
        observed = sha256_file(frozen_path)
        if observed != row.sha256:
            raise ContractError(
                f"Pre-outcome file drift ({row.name}): {observed} != {row.sha256}"
            )
    return frame


def write_stage_seal(
    output: Path,
    stage: str,
    artifact_names: Sequence[str],
) -> Path:
    """Write a last-step hash seal for an immutable pipeline stage."""
    seal_path = output / f"{stage}_stage_seal.tsv"
    if seal_path.exists():
        raise ContractError(f"Refusing to overwrite stage seal: {seal_path}")
    rows = []
    for name in artifact_names:
        path = output / name
        if not path.is_file():
            raise ContractError(f"Cannot seal missing {stage} artifact: {path}")
        rows.append(
            {
                "release_id": RELEASE_ID,
                "stage": stage,
                "relative_path": name,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    atomic_write_frame(seal_path, pd.DataFrame(rows))
    return seal_path


def validate_stage_seal(
    output: Path,
    stage: str,
    artifact_names: Sequence[str],
) -> str:
    """Require a complete stage seal and rederive every sealed artifact."""
    validate_analysis_freeze(output)
    seal_path = output / f"{stage}_stage_seal.tsv"
    if not seal_path.is_file():
        raise ContractError(f"Missing {stage} stage seal: {seal_path}")
    frame = pd.read_csv(seal_path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(
        frame,
        ["release_id", "stage", "relative_path", "bytes", "sha256"],
        f"{stage} stage seal",
    )
    if frame["relative_path"].duplicated().any():
        raise ContractError(f"{stage} stage seal contains duplicate artifacts")
    if set(frame["relative_path"]) != set(artifact_names):
        raise ContractError(f"{stage} stage artifact universe drift")
    if set(frame["release_id"]) != {RELEASE_ID} or set(frame["stage"]) != {stage}:
        raise ContractError(f"{stage} stage seal provenance drift")
    for row in frame.itertuples(index=False):
        artifact = output / row.relative_path
        if not artifact.is_file():
            raise ContractError(f"Sealed {stage} artifact disappeared: {artifact}")
        if artifact.stat().st_size != int(row.bytes) or sha256_file(artifact) != row.sha256:
            raise ContractError(f"Sealed {stage} artifact drift: {artifact}")
    return sha256_file(seal_path)


def validate_stage_chain(output: Path, through_stage: str) -> dict[str, str]:
    """Rederive every immutable stage seal through ``through_stage``.

    A downstream seal is not treated as a transitive content hash: it can
    refer to a manifest that names upstream products without rehashing those
    products.  Every entry point therefore validates the complete chain it
    consumes rather than trusting only the most recent seal.
    """
    if through_stage not in STAGE_ARTIFACTS:
        raise ContractError(f"Unknown Plan 11 stage: {through_stage!r}")
    observed: dict[str, str] = {}
    for stage in STAGE_ORDER[: STAGE_ORDER.index(through_stage) + 1]:
        observed[stage] = validate_stage_seal(
            output, stage, STAGE_ARTIFACTS[stage]
        )
    return observed


def read_one_row(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if len(frame) != 1:
        raise ContractError(f"Expected exactly one row in {path}, found {len(frame)}")
    return frame.iloc[0].to_dict()


def bool_value(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes"}:
        return True
    if text in {"false", "f", "0", "no", "", "nan", "na"}:
        return False
    raise ContractError(f"Cannot parse Boolean value: {value!r}")


def require_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = sorted(set(columns).difference(frame.columns))
    if missing:
        raise ContractError(f"{label} lacks required columns: {', '.join(missing)}")


def check_hash(path: Path, expected: str, label: str) -> str:
    if not path.is_file():
        raise ContractError(f"Missing {label}: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise ContractError(
            f"{label} hash drift: observed {observed}, expected {expected}: {path}"
        )
    return observed


@dataclass(frozen=True)
class RegistryContract:
    family: pd.DataFrame
    membership: pd.DataFrame
    registry_sha256: str
    membership_sha256: str
    validation_sha256: str
    release_manifest_sha256: str
    tested_universe_sha256: str
    ready_sha256: str


def validate_registry_contract(paths: dict[str, Path]) -> RegistryContract:
    registry_hash = check_hash(
        paths["registry"], EXPECTED_REGISTRY_SHA256, "Plan 20 registry"
    )
    membership_hash = check_hash(
        paths["membership"], EXPECTED_MEMBERSHIP_SHA256, "Plan 20 membership"
    )
    validation_hash = check_hash(
        paths["hotspot_validation"],
        EXPECTED_VALIDATION_SHA256,
        "Plan 20 validation status",
    )
    release_hash = check_hash(
        paths["hotspot_release_manifest"],
        EXPECTED_RELEASE_MANIFEST_SHA256,
        "Plan 20 release manifest",
    )
    ready = read_one_row(paths["hotspot_ready"])
    required_ready = {
        "release_id": RELEASE_ID,
        "status": "ready_for_external_testing",
        "registry_sha256": registry_hash,
        "membership_table_sha256": membership_hash,
        "validation_status_sha256": validation_hash,
        "release_manifest_sha256": release_hash,
        "external_outcomes_read": "FALSE",
    }
    for field, expected in required_ready.items():
        if ready.get(field) != expected:
            raise ContractError(
                f"Plan 20 READY mismatch for {field}: {ready.get(field)!r} != {expected!r}"
            )

    family = pd.read_csv(paths["external"], sep="\t", dtype=str, keep_default_na=False)
    require_columns(
        family,
        [
            "release_id",
            "registry_sha256",
            "program_uid",
            "membership_sha256",
            "cell_type",
            "module",
            "expected_direction",
        ],
        "external-test family",
    )
    modules = tuple(sorted(pd.to_numeric(family["module"]).astype(int).tolist()))
    if modules != EXPECTED_EXTERNAL_MODULES:
        raise ContractError(
            f"External family drift: modules={modules}, expected={EXPECTED_EXTERNAL_MODULES}"
        )
    if len(family) != len(EXPECTED_EXTERNAL_MODULES) or family["program_uid"].duplicated().any():
        raise ContractError("External-test family is not the expected unique two-program set")
    if set(family["cell_type"]) != {"hepatocytes"}:
        raise ContractError("External-test family contains a non-hepatocyte program")
    if set(family["release_id"]) != {RELEASE_ID}:
        raise ContractError("External-test family release ID drift")
    if set(family["registry_sha256"]) != {registry_hash}:
        raise ContractError("External-test family registry hash drift")
    if not family["expected_direction"].isin(["positive", "negative"]).all():
        raise ContractError("External-test family contains an invalid expected direction")

    membership = pd.read_csv(
        paths["membership"], sep="\t", dtype=str, keep_default_na=False
    )
    require_columns(
        membership,
        [
            "program_uid",
            "mapped_symbol",
            "mapped_symbol_status",
            "original_l1_weight",
            "membership_sha256",
        ],
        "Plan 20 membership",
    )
    membership = membership[membership["program_uid"].isin(family["program_uid"])].copy()
    if set(membership["program_uid"]) != set(family["program_uid"]):
        raise ContractError("At least one frozen external program has no membership rows")
    membership["original_l1_weight"] = pd.to_numeric(
        membership["original_l1_weight"], errors="raise"
    )
    if (
        ~np.isfinite(membership["original_l1_weight"])
        | (membership["original_l1_weight"] <= 0)
    ).any():
        raise ContractError("Frozen external program has a non-positive/non-finite weight")
    family_hashes = family.set_index("program_uid")["membership_sha256"].to_dict()
    for uid, part in membership.groupby("program_uid", sort=False):
        if set(part["membership_sha256"]) != {family_hashes[uid]}:
            raise ContractError(f"Membership hash mismatch within {uid}")

    if not paths["tested_universe"].is_file():
        raise ContractError(f"Missing Plan 20 tested universe: {paths['tested_universe']}")
    return RegistryContract(
        family=family.sort_values(["cell_type", "module"]).reset_index(drop=True),
        membership=membership.reset_index(drop=True),
        registry_sha256=registry_hash,
        membership_sha256=membership_hash,
        validation_sha256=validation_hash,
        release_manifest_sha256=release_hash,
        tested_universe_sha256=sha256_file(paths["tested_universe"]),
        ready_sha256=sha256_file(paths["hotspot_ready"]),
    )


@dataclass(frozen=True)
class SourceGateContract:
    row: dict[str, str]
    passing_donors: tuple[str, ...]
    source_gate_sha256: str


def validate_source_gate_summary(path: Path) -> SourceGateContract:
    row = read_one_row(path)
    status = row.get("terminal_verdict", row.get("status", ""))
    required = {
        "release_id": RELEASE_ID,
        "dataset": "yakubovsky2026",
        "ordinal_effect_interpretation": "binary_lipid_zone_vs_non_lipid_zone",
    }
    for field, expected in required.items():
        if row.get(field) != expected:
            raise ContractError(
                f"Plan 10 source gate mismatch for {field}: {row.get(field)!r} != {expected!r}"
            )
    if status != ALLOWED_SOURCE_GATE:
        raise ContractError(
            f"Binary Plan 11 pipeline requires {ALLOWED_SOURCE_GATE}, observed {status!r}"
        )
    if row.get("status", "") != status:
        raise ContractError("Plan 10 status and terminal_verdict disagree")
    required_true = (
        "source_gate_pass",
        "integrity_gate_pass",
        "zonation_gate_pass",
        "ordinal_lipid_gate_pass",
        "authoritative_barcode_join_verified",
        "authoritative_ordinal_lipid_join_verified",
    )
    for field in required_true:
        if not bool_value(row.get(field, "")):
            raise ContractError(f"Plan 10 terminal gate did not certify {field}")
    if bool_value(row.get("v2_registry_read", "")):
        raise ContractError(
            "Plan 10 source gate is not outcome/registry independent (v2_registry_read=TRUE)"
        )
    if bool_value(row.get("continuous_lipid_gate_pass", "")):
        raise ContractError("Binary pipeline refuses a continuous-lipid source gate")
    donors = tuple(
        donor for donor in row.get("passing_ordinal_donors", "").split(";") if donor
    )
    if len(donors) < 3 or len(set(donors)) != len(donors):
        raise ContractError(f"Source gate has invalid passing-donor set: {donors}")
    if int(row.get("n_passing_ordinal_donors", "0")) != len(donors):
        raise ContractError("Passing binary-donor count disagrees with donor list")
    donor_counts: dict[str, tuple[int, int, int]] = {}
    count_text = row.get("ordinal_donor_spot_counts_lipid_nonlipid_total", "")
    for token in filter(None, count_text.split(";")):
        try:
            donor, encoded = token.split(":", 1)
            n_lipid, n_non_lipid, n_total = map(int, encoded.split("/"))
        except (TypeError, ValueError) as error:
            raise ContractError(
                f"Malformed ordinal donor spot-count token: {token!r}"
            ) from error
        if donor in donor_counts:
            raise ContractError(f"Duplicate ordinal donor spot-count row: {donor}")
        donor_counts[donor] = (n_lipid, n_non_lipid, n_total)
    if set(donor_counts) != set(donors):
        raise ContractError(
            "Ordinal donor spot-count universe disagrees with passing donor list"
        )
    for donor, (n_lipid, n_non_lipid, n_total) in donor_counts.items():
        if (
            n_lipid <= 0
            or n_non_lipid <= 0
            or n_total < MIN_DONOR_SPOTS
            or n_lipid + n_non_lipid != n_total
        ):
            raise ContractError(
                f"Ordinal donor fails the binary nondegeneracy/spot gate: {donor} "
                f"({n_lipid}/{n_non_lipid}/{n_total})"
            )
    return SourceGateContract(
        row=row,
        passing_donors=donors,
        source_gate_sha256=sha256_file(path),
    )


def bh_adjust(values: Sequence[float]) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    p = values[valid]
    if not len(p):
        return result
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result[valid] = restored
    return result


def signed_stouffer(slopes: Sequence[float], pvalues: Sequence[float]) -> tuple[float, float]:
    slopes = np.asarray(slopes, dtype=float)
    pvalues = np.asarray(pvalues, dtype=float)
    valid = np.isfinite(slopes) & np.isfinite(pvalues) & (pvalues > 0) & (pvalues <= 1)
    if not valid.any():
        return np.nan, np.nan
    z = np.sign(slopes[valid]) * stats.norm.isf(np.clip(pvalues[valid], 1e-300, 1) / 2)
    combined = float(np.sum(z) / math.sqrt(len(z)))
    return combined, float(2 * stats.norm.sf(abs(combined)))


def decode_matlab_char(dataset: h5py.Dataset) -> str:
    values = np.asarray(dataset[...]).reshape(-1, order="F")
    return "".join(chr(int(value)) for value in values if int(value) != 0)


@dataclass(frozen=True)
class BulkCellstrResult:
    """Decoded MATLAB cell-string axes plus resolver provenance."""

    values: dict[str, list[str]]
    n_axes: int
    n_references: int
    n_unique_objects: int
    n_objects_visited: int


def matlab_object_reference_addresses(dataset: h5py.Dataset) -> np.ndarray:
    """Read object references as their exact HDF5 object-header addresses.

    High-level h5py object-reference arrays are normally dereferenced one by
    one.  MATLAB cell arrays with hundreds of thousands of strings make that
    random traversal prohibitively slow on network storage.  HDF5 object
    references are eight-byte object-header addresses; reading them through
    the declared reference memory type lets ``bulk_decode_matlab_cellstr``
    resolve the unique target set in one sequential ``visititems`` pass.
    """

    if h5py.check_dtype(ref=dataset.dtype) is None:
        raise ContractError(f"Expected object-reference dataset at {dataset.name}")
    addresses = np.empty(dataset.shape, dtype=np.uint64)
    dataset.id.read(
        h5py.h5s.ALL,
        h5py.h5s.ALL,
        addresses,
        mtype=h5py.h5t.STD_REF_OBJ,
    )
    return addresses.reshape(-1, order="F")


def bulk_decode_matlab_cellstr(
    datasets: Mapping[str, h5py.Dataset],
    handle: h5py.File,
    *,
    allow_null: bool = False,
) -> BulkCellstrResult:
    """Decode multiple MATLAB cell-string axes with one object-tree walk.

    Only referenced character datasets are read.  Other objects encountered
    by ``visititems`` contribute metadata addresses but their values are never
    opened, which keeps this resolver suitable for outcome-blind schema work.
    """

    if not datasets:
        return BulkCellstrResult({}, 0, 0, 0, 0)
    address_axes: dict[str, np.ndarray] = {}
    for name, dataset in datasets.items():
        if name in address_axes:
            raise ContractError(f"Duplicate MATLAB cell-string axis name: {name}")
        address_axes[name] = matlab_object_reference_addresses(dataset)
    n_references = sum(len(values) for values in address_axes.values())
    null_count = sum(int(np.sum(values == 0)) for values in address_axes.values())
    if null_count and not allow_null:
        raise ContractError(
            f"MATLAB cell-string axes contain {null_count} null object reference(s)"
        )
    targets = {
        int(address)
        for values in address_axes.values()
        for address in values
        if int(address) != 0
    }
    unresolved = set(targets)
    paths: dict[int, str] = {}
    n_objects_visited = 0

    def visitor(_name: str, obj: h5py.Group | h5py.Dataset) -> bool | None:
        nonlocal n_objects_visited
        n_objects_visited += 1
        address = int(h5py.h5o.get_info(obj.id).addr)
        if address not in unresolved:
            return None
        if not isinstance(obj, h5py.Dataset):
            raise ContractError(
                f"MATLAB cell-string reference targets non-dataset object {obj.name}"
            )
        paths[address] = obj.name
        unresolved.remove(address)
        return True if not unresolved else None

    if unresolved:
        handle.visititems(visitor)
    if unresolved:
        preview = ",".join(str(value) for value in sorted(unresolved)[:5])
        raise ContractError(
            f"Could not resolve {len(unresolved)} MATLAB object reference(s): {preview}"
        )
    decoded_by_address: dict[int, str] = {}
    for address, path in paths.items():
        target = handle[path]
        if not isinstance(target, h5py.Dataset):
            raise ContractError(f"Resolved MATLAB character target is not a dataset: {path}")
        if target.dtype.kind not in {"u", "i"}:
            raise ContractError(
                f"MATLAB character target has non-integer dtype {target.dtype}: {path}"
            )
        decoded_by_address[address] = decode_matlab_char(target)
    values = {
        name: [
            "" if int(address) == 0 else decoded_by_address[int(address)]
            for address in addresses
        ]
        for name, addresses in address_axes.items()
    }
    return BulkCellstrResult(
        values=values,
        n_axes=len(address_axes),
        n_references=n_references,
        n_unique_objects=len(targets),
        n_objects_visited=n_objects_visited,
    )


def decode_matlab_cellstr(dataset: h5py.Dataset, handle: h5py.File) -> list[str]:
    strings: list[str] = []
    for reference in np.asarray(dataset[...]).reshape(-1, order="F"):
        if not reference:
            strings.append("")
            continue
        target = handle[reference]
        if not isinstance(target, h5py.Dataset):
            raise ContractError(f"Expected MATLAB character dataset at {target.name}")
        strings.append(decode_matlab_char(target))
    return strings


def gene_axis_sha256(values: Sequence[str]) -> str:
    """Hash an ordered decoded gene axis with an unambiguous serialization."""

    digest = hashlib.sha256()
    for value in values:
        encoded = str(value).encode("utf-8")
        digest.update(struct.pack("<Q", len(encoded)))
        digest.update(encoded)
    return digest.hexdigest()


@dataclass(frozen=True)
class GeneAxisContract:
    rows: pd.DataFrame
    master_donor: str
    master_group_path: str
    n_genes: int
    common_sha256: str
    v_mat_sha256: str
    manifest_sha256: str


def validate_gene_axis_manifest(
    path: Path,
    observed_v_mat_sha256: str,
) -> GeneAxisContract:
    """Validate the exact, source-SHA-bound all-donor gene-axis certificate."""

    if not path.is_file():
        raise ContractError(f"Missing Yakubovsky gene-axis manifest: {path}")
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(
        frame,
        [
            "axis_order",
            "donor",
            "group_path",
            "n_genes",
            "gene_axis_sha256",
            "common_gene_axis_sha256",
            "v_mat_sha256",
            "hash_algorithm",
            "zero_policy",
            "master_axis",
            "verification_status",
            "verified_outcome_fields_read",
        ],
        "Yakubovsky gene-axis manifest",
    )
    if frame.empty or frame["donor"].duplicated().any() or frame["group_path"].duplicated().any():
        raise ContractError("Yakubovsky gene-axis manifest has empty or duplicate axes")
    order = pd.to_numeric(frame["axis_order"], errors="raise").astype(int).tolist()
    if order != list(range(1, len(frame) + 1)):
        raise ContractError("Yakubovsky gene-axis order is not contiguous and deterministic")
    n_genes = pd.to_numeric(frame["n_genes"], errors="raise").astype(int)
    if (n_genes <= 0).any() or n_genes.nunique() != 1:
        raise ContractError("Yakubovsky donor gene-axis lengths are not identical and positive")
    common_hashes = set(frame["common_gene_axis_sha256"])
    if len(common_hashes) != 1 or set(frame["gene_axis_sha256"]) != common_hashes:
        raise ContractError("Yakubovsky donor gene-axis content hashes are not identical")
    if set(frame["v_mat_sha256"]) != {observed_v_mat_sha256}:
        raise ContractError("Yakubovsky gene-axis certificate targets a different v.mat")
    if set(frame["hash_algorithm"]) != {GENE_AXIS_HASH_ALGORITHM}:
        raise ContractError("Yakubovsky gene-axis hash algorithm drift")
    if set(frame["zero_policy"]) != {GENE_AXIS_ZERO_POLICY}:
        raise ContractError("Yakubovsky gene-axis zero-normalization policy drift")
    if set(frame["verification_status"]) != {"exact_decoded_content_equal"}:
        raise ContractError("Yakubovsky gene-axis certificate lacks exact-content status")
    if any(bool_value(value) for value in frame["verified_outcome_fields_read"]):
        raise ContractError("Yakubovsky gene-axis certificate is not outcome-blind")
    masters = frame[frame["master_axis"].str.lower().isin(["true", "1"])]
    if len(masters) != 1 or int(masters.iloc[0]["axis_order"]) != 1:
        raise ContractError("Yakubovsky gene-axis manifest must identify its first axis as master")
    return GeneAxisContract(
        rows=frame.copy(),
        master_donor=masters.iloc[0]["donor"],
        master_group_path=masters.iloc[0]["group_path"],
        n_genes=int(n_genes.iloc[0]),
        common_sha256=next(iter(common_hashes)),
        v_mat_sha256=observed_v_mat_sha256,
        manifest_sha256=sha256_file(path),
    )


def numeric_vector(dataset: h5py.Dataset) -> np.ndarray:
    return np.asarray(dataset[...], dtype=float).reshape(-1, order="F")


def parse_h5ad_sparse(group: h5py.Group) -> sparse.csr_matrix:
    required = {"data", "indices", "indptr"}
    if not required.issubset(group.keys()):
        raise ContractError(f"H5AD sparse group {group.name} lacks {sorted(required)}")
    shape = tuple(int(x) for x in np.asarray(group.attrs["shape"]).tolist())
    encoding = group.attrs.get("encoding-type", "")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    cls = sparse.csc_matrix if encoding == "csc_matrix" else sparse.csr_matrix
    matrix = cls(
        (
            np.asarray(group["data"]),
            np.asarray(group["indices"]),
            np.asarray(group["indptr"]),
        ),
        shape=shape,
    )
    return matrix.tocsr()


def _decode_string_array(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values)
    decoded: list[object] = []
    for item in array.reshape(-1):
        if isinstance(item, (bytes, np.bytes_)):
            decoded.append(bytes(item).decode("utf-8"))
        elif item is None:
            decoded.append(None)
        else:
            decoded.append(str(item))
    return np.asarray(decoded, dtype=object).reshape(array.shape)


def read_h5ad_dataframe_column(group: h5py.Group, column: str) -> np.ndarray:
    obj = group[column]
    if isinstance(obj, h5py.Dataset):
        values = np.asarray(obj)
        if values.dtype.kind in {"S", "U", "O"}:
            return _decode_string_array(values)
        return values
    if isinstance(obj, h5py.Group) and {"codes", "categories"}.issubset(obj.keys()):
        codes = np.asarray(obj["codes"], dtype=int)
        categories = _decode_string_array(np.asarray(obj["categories"]))
        return np.asarray(
            [categories[code] if code >= 0 else None for code in codes], dtype=object
        )
    raise ContractError(f"Unsupported H5AD dataframe encoding: {obj.name}")


def read_h5ad_dataframe_index(group: h5py.Group) -> np.ndarray:
    """Read an AnnData dataframe index using its encoded index-field name.

    AnnData stores the physical index dataset under the pandas index name when
    that name is non-null (for example ``gene_symbol``), and records that name
    in the dataframe group's ``_index`` attribute.  Hard-coding a dataset named
    ``_index`` therefore fails for a valid named index.
    """
    encoded_name: Any = group.attrs.get("_index", "_index")
    if isinstance(encoded_name, np.ndarray) and encoded_name.shape == ():
        encoded_name = encoded_name.item()
    if isinstance(encoded_name, bytes):
        encoded_name = encoded_name.decode()
    index_name = str(encoded_name)
    if index_name not in group:
        raise ContractError(
            f"H5AD dataframe index field {index_name!r} is absent from {group.name}"
        )
    return read_h5ad_dataframe_column(group, index_name)


def standardize_columns(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    mean = np.nanmean(values, axis=0, keepdims=True)
    sd = np.nanstd(values, axis=0, ddof=1, keepdims=True)
    valid = np.isfinite(sd) & (sd > 0)
    result = np.zeros_like(values, dtype=float)
    np.divide(values - mean, sd, out=result, where=valid)
    result[~np.isfinite(result)] = 0.0
    return result


def score_frozen_program(
    standardized_gene_values: np.ndarray,
    original_l1_weights: Sequence[float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, int]:
    """Return frozen-weight, equal-weight, and leave-top-gene scores.

    The input genes must already be restricted to the measured frozen members
    and ordered identically to ``original_l1_weights``.  No outcome is accepted
    by this function.
    """
    values = np.asarray(standardized_gene_values, dtype=float)
    weights = np.asarray(original_l1_weights, dtype=float)
    if values.ndim != 2 or values.shape[1] != len(weights):
        raise ContractError("Program-score matrix and weight dimensions disagree")
    if not len(weights) or not np.isfinite(weights).all() or (weights <= 0).any():
        raise ContractError("Program-score weights must be positive and finite")
    normalized = weights / weights.sum()
    primary = values @ normalized
    equal = np.mean(values, axis=1)
    top = int(np.argmax(weights))
    leave: np.ndarray | None
    if len(weights) - 1 < MIN_PROGRAM_GENES:
        leave = None
    else:
        keep = np.arange(len(weights)) != top
        leave_weights = weights[keep] / weights[keep].sum()
        leave = values[:, keep] @ leave_weights
    return primary, equal, leave, top


def natural_spline_design(zonation: np.ndarray, df: int) -> np.ndarray:
    zonation = np.asarray(zonation, dtype=float)
    if np.unique(zonation[np.isfinite(zonation)]).size < max(4, df):
        raise ContractError("Zonation is too degenerate for the prespecified spline")
    basis = np.asarray(
        dmatrix(
            f"cr(z, df={int(df)}, constraints='center')",
            {"z": zonation},
            return_type="dataframe",
        ),
        dtype=float,
    )
    if basis.shape[1] != df + 1:
        raise ContractError(
            f"Unexpected centered natural-spline design width {basis.shape[1]} for df={df}"
        )
    return basis


def build_model_design(
    lipid_zone: np.ndarray,
    zonation: np.ndarray,
    background_corrected_spot_sum: np.ndarray,
    detected_genes: np.ndarray,
    spline_df: int,
    zonation_landmark_score: np.ndarray | None = None,
) -> tuple[np.ndarray, list[str]]:
    lipid_zone = np.asarray(lipid_zone, dtype=float)
    zonation = np.asarray(zonation, dtype=float)
    background_corrected_spot_sum = np.asarray(
        background_corrected_spot_sum, dtype=float
    )
    detected_genes = np.asarray(detected_genes, dtype=float)
    spline = natural_spline_design(zonation, spline_df)
    columns = ["intercept", "lipid_zone"] + [
        f"zonation_ns_{index + 1}" for index in range(spline_df)
    ] + ["centered_log1p_background_corrected_spot_sum", "centered_detected_genes"]
    covariates = [
        np.log1p(background_corrected_spot_sum)
        - np.nanmean(np.log1p(background_corrected_spot_sum)),
        detected_genes - np.nanmean(detected_genes),
    ]
    if zonation_landmark_score is not None:
        zonation_landmark_score = np.asarray(zonation_landmark_score, dtype=float)
        covariates.append(
            zonation_landmark_score - np.nanmean(zonation_landmark_score)
        )
        columns.append("centered_source_zonation_landmark_expression")
    # spline[:, 0] is the intercept emitted by patsy; retain it once.
    design = np.column_stack([spline[:, 0], lipid_zone, spline[:, 1:], *covariates])
    return design, columns


@dataclass(frozen=True)
class OLSFit:
    beta: float
    se: float
    statistic: float
    pvalue: float
    residual_df: int
    rank: int
    n_columns: int
    condition_number: float


def ols_lipid_fit(design: np.ndarray, outcome: np.ndarray) -> OLSFit:
    design = np.asarray(design, dtype=float)
    outcome = np.asarray(outcome, dtype=float)
    if len(outcome) != len(design) or not np.isfinite(design).all() or not np.isfinite(outcome).all():
        raise ContractError("Non-finite or mismatched model arrays")
    rank = int(np.linalg.matrix_rank(design))
    if rank != design.shape[1]:
        raise ContractError(f"Model is rank deficient: rank={rank}, columns={design.shape[1]}")
    coef, _, _, singular = np.linalg.lstsq(design, outcome, rcond=None)
    residual = outcome - design @ coef
    residual_df = len(outcome) - design.shape[1]
    if residual_df <= 0:
        raise ContractError("Model has no residual degrees of freedom")
    sigma2 = float(np.sum(residual * residual) / residual_df)
    covariance = sigma2 * np.linalg.inv(design.T @ design)
    se = float(np.sqrt(max(covariance[1, 1], 0)))
    statistic = float(coef[1] / se) if se > 0 else np.nan
    pvalue = float(2 * stats.t.sf(abs(statistic), df=residual_df)) if np.isfinite(statistic) else np.nan
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else np.inf
    return OLSFit(
        beta=float(coef[1]),
        se=se,
        statistic=statistic,
        pvalue=pvalue,
        residual_df=residual_df,
        rank=rank,
        n_columns=design.shape[1],
        condition_number=condition,
    )


@dataclass(frozen=True)
class SpatialBlocks:
    labels: np.ndarray
    audit: pd.DataFrame
    n_blocks: int


def _connected_islands(coords: np.ndarray) -> tuple[np.ndarray, float, float]:
    if len(coords) < MIN_BLOCK_SPOTS:
        raise ContractError("Section has too few spots for a tissue graph")
    tree = cKDTree(coords)
    k = min(7, len(coords))
    distances, neighbors = tree.query(coords, k=k)
    nearest = distances[:, 1]
    nearest = nearest[np.isfinite(nearest) & (nearest > 0)]
    if not len(nearest):
        raise ContractError("Section has no positive nearest-neighbor distance")
    median_nearest = float(np.median(nearest))
    max_distance = GRAPH_DISTANCE_MULTIPLIER * median_nearest
    edges: set[tuple[int, int]] = set()
    for left, (distance_row, neighbor_row) in enumerate(zip(distances, neighbors, strict=True)):
        for distance, right in zip(distance_row[1:], neighbor_row[1:], strict=True):
            if np.isfinite(distance) and distance <= max_distance:
                a, b = sorted((int(left), int(right)))
                if a != b:
                    edges.add((a, b))
    if not edges:
        raise ContractError("No within-tissue graph edges remain")
    edge = np.asarray(sorted(edges), dtype=int)
    graph = sparse.coo_matrix(
        (
            np.ones(2 * len(edge), dtype=np.int8),
            (
                np.concatenate([edge[:, 0], edge[:, 1]]),
                np.concatenate([edge[:, 1], edge[:, 0]]),
            ),
        ),
        shape=(len(coords), len(coords)),
    ).tocsr()
    _, labels = connected_components(graph, directed=False, return_labels=True)
    return labels.astype(int), median_nearest, max_distance


def build_spatial_blocks(
    donor: str,
    section: Sequence[str],
    x: Sequence[float],
    y: Sequence[float],
) -> SpatialBlocks:
    section = np.asarray(section, dtype=object)
    coords = np.column_stack([np.asarray(x, dtype=float), np.asarray(y, dtype=float)])
    if not np.isfinite(coords).all():
        raise ContractError(f"Non-finite coordinates for donor {donor}")
    labels = np.full(len(coords), -1, dtype=int)
    audit_rows: list[dict[str, Any]] = []
    next_block = 0
    for section_id in pd.unique(section):
        section_idx = np.flatnonzero(section == section_id)
        islands, median_nearest, max_distance = _connected_islands(coords[section_idx])
        for island_id in np.unique(islands):
            local = np.flatnonzero(islands == island_id)
            global_idx = section_idx[local]
            island_coords = coords[global_idx]
            if len(island_coords) >= 2:
                island_nearest = cKDTree(island_coords).query(island_coords, k=2)[0][:, 1]
                island_nearest = island_nearest[
                    np.isfinite(island_nearest) & (island_nearest > 0)
                ]
                island_median_nearest = (
                    float(np.median(island_nearest))
                    if len(island_nearest)
                    else median_nearest
                )
            else:
                island_median_nearest = median_nearest
            if len(global_idx) < MIN_BLOCK_SPOTS:
                audit_rows.append(
                    {
                        "donor": donor,
                        "section": str(section_id),
                        "island_id": int(island_id),
                        "block_id": -1,
                        "block_uid": f"{donor}:{section_id}:excluded_island_{int(island_id)}",
                        "tile_x": "",
                        "tile_y": "",
                        "n_spots": len(global_idx),
                        "median_nearest_distance": island_median_nearest,
                        "graph_max_edge_distance": max_distance,
                        "block_width": BLOCK_WIDTH_MULTIPLIER * island_median_nearest,
                        "n_small_block_merges_in_island": 0,
                        "small_block_merge_trace": "",
                        "excluded_small_island": True,
                    }
                )
                continue
            width = BLOCK_WIDTH_MULTIPLIER * island_median_nearest
            origin = np.min(island_coords, axis=0)
            tiles = np.floor((island_coords - origin) / width).astype(int)
            keys = [tuple(row) for row in tiles]
            unique_keys = sorted(set(keys))
            block_members = {
                key: np.asarray([i for i, value in enumerate(keys) if value == key], dtype=int)
                for key in unique_keys
            }
            merge_rows: list[tuple[tuple[int, int], tuple[int, int]]] = []
            small = [key for key, members in block_members.items() if len(members) < MIN_BLOCK_SPOTS]
            for key in sorted(small, key=lambda item: (len(block_members[item]), item)):
                if key not in block_members or len(block_members) == 1:
                    continue
                # A small tile can already have absorbed another small tile
                # earlier in this deterministic pass.  Once it reaches the
                # prespecified minimum it must remain a block; merging it again
                # would make the result depend on the stale initial ``small``
                # list and needlessly erase spatial resolution.
                if len(block_members[key]) >= MIN_BLOCK_SPOTS:
                    continue
                source_centroid = np.mean(island_coords[block_members[key]], axis=0)
                candidates = [candidate for candidate in block_members if candidate != key]
                target = min(
                    candidates,
                    key=lambda candidate: (
                        np.linalg.norm(
                            source_centroid
                            - np.mean(island_coords[block_members[candidate]], axis=0)
                        ),
                        candidate,
                    ),
                )
                block_members[target] = np.concatenate([block_members[target], block_members[key]])
                del block_members[key]
                merge_rows.append((key, target))
            for key in sorted(block_members):
                members = global_idx[block_members[key]]
                labels[members] = next_block
                audit_rows.append(
                    {
                        "donor": donor,
                        "section": str(section_id),
                        "island_id": int(island_id),
                        "block_id": int(next_block),
                        "block_uid": f"{donor}:{section_id}:block_{int(next_block)}",
                        "tile_x": key[0],
                        "tile_y": key[1],
                        "n_spots": len(members),
                        "median_nearest_distance": island_median_nearest,
                        "graph_max_edge_distance": max_distance,
                        "block_width": width,
                        "n_small_block_merges_in_island": len(merge_rows),
                        "small_block_merge_trace": ";".join(
                            f"{source[0]},{source[1]}->{target[0]},{target[1]}"
                            for source, target in merge_rows
                        ),
                        "excluded_small_island": False,
                    }
                )
                next_block += 1
    audit = pd.DataFrame(audit_rows, columns=SPATIAL_BLOCK_AUDIT_COLUMNS)
    if audit["n_spots"].sum() != len(labels):
        raise ContractError(f"Block audit does not reconcile spot count for donor {donor}")
    if not np.array_equal(np.flatnonzero(labels < 0), np.flatnonzero(labels == -1)):
        raise ContractError(f"Invalid excluded-block labels for donor {donor}")
    return SpatialBlocks(labels=labels, audit=audit, n_blocks=next_block)


@dataclass(frozen=True)
class BootstrapResult:
    estimates: np.ndarray
    ci_low: float
    ci_high: float
    pvalue: float
    n_failed: int


def block_bootstrap_lipid_slope(
    design: np.ndarray,
    outcome: np.ndarray,
    blocks: np.ndarray,
    sections: Sequence[str],
    n_bootstrap: int,
    seed: int,
) -> BootstrapResult:
    design = np.asarray(design, dtype=float)
    outcome = np.asarray(outcome, dtype=float)
    blocks = np.asarray(blocks, dtype=int)
    sections = np.asarray(sections, dtype=object)
    if not (len(design) == len(outcome) == len(blocks) == len(sections)):
        raise ContractError("Bootstrap arrays have mismatched lengths")
    rng = np.random.default_rng(seed)
    sufficient_by_section: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for section in pd.unique(sections):
        section_name = str(section)
        section_blocks = np.unique(blocks[sections == section])
        xtx = []
        xty = []
        for block_id in section_blocks:
            member = (sections == section) & (blocks == block_id)
            block_x = design[member]
            block_y = outcome[member]
            xtx.append(block_x.T @ block_x)
            xty.append(block_x.T @ block_y)
        sufficient_by_section[section_name] = (
            section_blocks,
            np.stack(xtx),
            np.stack(xty),
        )
    estimates = np.full(n_bootstrap, np.nan, dtype=float)
    for draw in range(n_bootstrap):
        weighted_xtx = np.zeros((design.shape[1], design.shape[1]), dtype=float)
        weighted_xty = np.zeros(design.shape[1], dtype=float)
        for _section, (section_blocks, block_xtx, block_xty) in sufficient_by_section.items():
            counts = rng.multinomial(
                len(section_blocks), np.repeat(1 / len(section_blocks), len(section_blocks))
            )
            weighted_xtx += np.tensordot(counts, block_xtx, axes=(0, 0))
            weighted_xty += np.tensordot(counts, block_xty, axes=(0, 0))
        if np.linalg.matrix_rank(weighted_xtx) != design.shape[1]:
            continue
        try:
            estimates[draw] = np.linalg.solve(weighted_xtx, weighted_xty)[1]
        except np.linalg.LinAlgError:
            continue
    valid = estimates[np.isfinite(estimates)]
    n_failed = n_bootstrap - len(valid)
    if not len(valid):
        return BootstrapResult(estimates, np.nan, np.nan, np.nan, n_failed)
    lower = float(np.quantile(valid, 0.025))
    upper = float(np.quantile(valid, 0.975))
    non_positive = int(np.sum(valid <= 0))
    non_negative = int(np.sum(valid >= 0))
    pvalue = min(
        1.0,
        # Singular resamples are explicitly recorded and may be tolerated up
        # to the frozen failure-rate gate.  They carry no tail information,
        # so the finite-sample correction is based on the valid resamples,
        # not the requested count.  Using ``n_bootstrap`` here would make the
        # P value anti-conservative whenever even one draw failed.
        2 * (1 + min(non_positive, non_negative)) / (len(valid) + 1),
    )
    return BootstrapResult(estimates, lower, upper, float(pvalue), n_failed)


def expected_sign(direction: str) -> int:
    if direction == "positive":
        return 1
    if direction == "negative":
        return -1
    raise ContractError(f"Unknown expected direction: {direction!r}")


def direction_label(value: float) -> str:
    if not np.isfinite(value) or value == 0:
        return "none"
    return "positive" if value > 0 else "negative"


def source_v1_manifest_state(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(frame, ["path", "sha256"], "v1 baseline manifest")
    return frame

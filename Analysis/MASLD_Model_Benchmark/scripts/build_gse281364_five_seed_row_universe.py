#!/usr/bin/env python3
"""Build the outcome-blind common row universe for GSE281364 stacking."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import tomllib
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import (
    ArtifactError,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
    write_text_exclusive,
)


SCHEMA = "masld-bench-gse281364-five-seed-row-universe-v1"
STATUS = "row_universe_defined_stacking_blocked"
SEEDS = (1103, 2909, 4721, 6673, 8111)
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
SOURCE_IDS = (
    "dna_lm_common_heads",
    "nucleotide_transformer_common_head",
    "mpralegnet_native",
    "restricted_native_sequence",
    "borzoi_native",
    "corgi_regular_context",
    "blocked_context_candidates",
)
ROW_FIELDS = (
    "seed",
    "row_hash",
    "unit_hash",
    "block_hash",
    "stratum",
    "outer_fold",
    "study_id",
    "assay_context_id",
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
)
COMPONENT_FIELDS = (
    "seed",
    "row_hash",
    "unit_hash",
    "block_hash",
    "stratum",
    "outer_fold",
    "study_id",
    "observed",
    "sequence_prediction",
    "context_prediction",
)
SOURCE_FIELDS = (
    "source_id",
    "source_family_role",
    "model_ids",
    "current_row_count_per_context",
    "current_seed_ids",
    "current_fold_map",
    "current_scale",
    "current_open_trigger_eligible",
    "eligible_after_gates",
    "diagnostic_only",
    "blockers",
)
GROUP_FIELDS = (
    "outer_locus_sequence_group_id",
    "outer_fold",
    "contig",
    "elements",
    "minimum_variant_pos0",
    "maximum_variant_pos0",
)
SEI_FIELDS = (
    "fixture_id",
    "element_id",
    "outer_locus_sequence_group_id",
    "outer_fold",
    "contig",
    "variant_pos0",
    "variant_pos1",
    "input_start0",
    "input_end0",
    "variant_index0",
    "ref",
    "alt",
    "reference_sequence_sha256",
    "alternative_sequence_sha256",
    "reference_reverse_complement_sha256",
    "alternative_reverse_complement_sha256",
)
BORZOI_MAP_FIELDS = (
    "source_locus_group_id",
    "borzoi_long_range_group_id",
    "outer_fold",
    "contig",
    "minimum_variant_pos0",
    "maximum_variant_pos0",
    "source_elements",
)
BORZOI_MANIFEST_FIELDS = (
    "fixture_id",
    "element_id",
    "source_locus_group_id",
    "borzoi_long_range_group_id",
    "outer_fold",
    "contig",
    "variant_pos0",
    "variant_pos1",
    "input_start0",
    "input_end0",
    "variant_index0",
    "ref",
    "alt",
    "reference_sequence_sha256",
    "alternative_sequence_sha256",
    "reference_reverse_complement_sha256",
    "alternative_reverse_complement_sha256",
    "allele_effect_sign",
)
BORZOI_COMPONENT_FIELDS = (
    "borzoi_long_range_group_id",
    "outer_fold",
    "contig",
    "minimum_variant_pos0",
    "maximum_variant_pos0",
    "source_locus_groups",
    "source_elements",
)
FILE_AUTHORITIES = {
    "variant_task",
    "dna_language_registry",
    "regulatory_sequence_registry",
    "regulatory_local_registry",
    "stacking_implementation",
    "conditional_preflight",
}
TREE_AUTHORITIES = {
    "outcome_blind_split",
    "common_element_selection",
    "long_range_fold_map",
    "dna_lm_features",
    "caduceus_features",
    "mpralegnet_predictions",
    "sei_predictions",
    "enformer_predictions",
    "dna_lm_one_seed_evaluation",
    "caduceus_one_seed_evaluation",
    "borzoi_execution_preflight",
    "complementarity_readiness_audit",
}
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class RowUniverseError(RuntimeError):
    """Raised when a row, fold, source, or separation authority differs."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hash_id(*parts: str) -> str:
    return hashlib.sha256("\0".join(parts).encode("utf-8")).hexdigest()


def canonical_set_hash(label: str, values: Iterable[str]) -> str:
    payload = label + "\n" + "\n".join(sorted(values)) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise RowUniverseError(f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RowUniverseError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise RowUniverseError(f"{label} must be a JSON object")
    return value


def load_toml(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise RowUniverseError(f"{label} is not a regular file")
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise RowUniverseError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise RowUniverseError(f"{label} must be a TOML table")
    return value


def safe_relative_path(root: Path, configured: object, *, label: str) -> Path:
    if not isinstance(configured, str) or not configured:
        raise RowUniverseError(f"{label} path is missing")
    relative = Path(configured)
    if relative.is_absolute() or ".." in relative.parts:
        raise RowUniverseError(f"unsafe {label} path: {configured}")
    try:
        path = reject_symlink_components(root / relative, label=label).resolve(strict=True)
        path.relative_to(root)
    except (ArtifactError, OSError, ValueError) as error:
        raise RowUniverseError(f"{label} is missing or escapes project root") from error
    return path


def read_tsv(path: Path, fields: Sequence[str], *, label: str) -> list[dict[str, str]]:
    if path.is_symlink() or not path.is_file():
        raise RowUniverseError(f"{label} is not a regular file")
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != tuple(fields):
                raise RowUniverseError(f"{label} header differs")
            rows = [dict(row) for row in reader]
    except (OSError, UnicodeDecodeError, csv.Error) as error:
        raise RowUniverseError(f"cannot read {label}: {error}") from error
    if not rows:
        raise RowUniverseError(f"{label} is empty")
    return rows


def tsv_text(fields: Sequence[str], rows: Sequence[Mapping[str, object]]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row[field] for field in fields})
    return output.getvalue()


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != STATUS
        or config.get("dataset_id") != "gse281364"
        or config.get("task_id") != "variant_to_regulation"
        or config.get("endpoint_id") != "mpra_allelic_direction"
        or config.get("evaluator_id") != "locus_heldout_allelic_spearman_v1"
        or config.get("row_universe_build_authorized") is not True
        or config.get("prediction_generation_authorized") is not False
        or config.get("outcome_access_authorized") is not False
        or config.get("stack_fit_authorized") is not False
        or config.get("residual_correlation_authorized") is not False
        or config.get("conditional_model_build_authorized") is not False
        or config.get("sealed_asset_access_authorized") is not False
    ):
        raise RowUniverseError("top-level row-universe boundary differs")

    row = config.get("row_contract", {})
    if (
        row.get("study_id") != "gse281364"
        or row.get("stratum") != "all"
        or row.get("selected_elements_per_source_locus_group") != 1
        or row.get("source_locus_groups") != 1033
        or tuple(row.get("assay_context_ids", ())) != CONTEXTS
        or row.get("base_rows_per_seed") != 2066
        or tuple(row.get("fixed_seeds", ())) != SEEDS
        or row.get("seeded_rows") != 10330
        or row.get("biological_donors") != 0
        or row.get("experimental_replicates_per_assay_context") != 4
        or tuple(row.get("component_fields", ())) != COMPONENT_FIELDS
        or tuple(row.get("row_universe_fields", ())) != ROW_FIELDS
    ):
        raise RowUniverseError("row identity or denominator contract differs")

    fold = config.get("fold_contract", {})
    if (
        fold.get("outer_folds") != 5
        or fold.get("largest_receptive_field_buffer_bp") != 524_288
        or fold.get("source_locus_group_count") != 1033
        or fold.get("long_range_block_count") != 239
        or fold.get("expected_source_groups_reassigned_from_original_6kb_fold_map")
        != 830
        or fold.get("original_6kb_fold_predictions_directly_stackable") is not False
    ):
        raise RowUniverseError("long-range fold contract differs")
    expected_folds = {f"fold-{index}" for index in range(5)}
    if set(fold.get("expected_fold_source_locus_groups", {})) != expected_folds:
        raise RowUniverseError("long-range fold roster differs")

    outcome = config.get("outcome_contract", {})
    prediction = config.get("prediction_contract", {})
    gates = config.get("admission_gates", {})
    if (
        outcome.get("assay_native_unit")
        != "mean_signed_log2_ALT_minus_REF_RNA_over_DNA_activity"
        or outcome.get("allele_effect_sign") != "ALT_minus_REF"
        or outcome.get("all_4359_constructs_are_rows") is not False
        or outcome.get("within_locus_construct_aggregation_allowed") is not False
        or outcome.get("experimental_replicates_are_independent_donors") is not False
        or prediction.get("required_unit")
        != "out_of_fold_predicted_mean_signed_log2_ALT_minus_REF_RNA_over_DNA_activity"
        or prediction.get("required_sign") != "ALT_minus_REF"
        or tuple(prediction.get("required_contexts", ())) != CONTEXTS
        or prediction.get("static_native_score_directly_stackable") is not False
        or prediction.get("same_rows_all_models") is not True
        or prediction.get("same_rows_all_seeds") is not True
        or prediction.get("same_outer_fold_all_models") is not True
        or prediction.get("missing_prediction_allowed") is not False
        or prediction.get("prediction_as_zero_allowed") is not False
        or prediction.get("prediction_values_read_by_this_builder") is not False
        or tuple(gates.get("allowed_checkpoint_exposure_states", ()))
        != ("clean_declared", "target_label_unexposed")
        or gates.get("restricted_or_unknown_terms_may_enter_open_trigger") is not False
        or gates.get("five_fixed_seeds_required") is not True
        or gates.get("minimum_qualifying_seed_count") != 4
        or gates.get("conditional_model_may_supply_context_component") is not False
        or gates.get("sealed_results_may_select_components") is not False
    ):
        raise RowUniverseError("scale, missingness, contamination, or license gate differs")

    sources = config.get("source_status", ())
    if tuple(source.get("source_id") for source in sources) != SOURCE_IDS:
        raise RowUniverseError("source identity roster differs")
    if any(source.get("current_open_trigger_eligible") is not False for source in sources):
        raise RowUniverseError("no current source may be trigger eligible")
    context_sources = [source for source in sources if source.get("source_family_role") == "context"]
    if not context_sources or any(source.get("current_row_count_per_context") != 0 for source in context_sources):
        raise RowUniverseError("a current context row source cannot be implied")
    blocked = next(source for source in sources if source["source_id"] == "blocked_context_candidates")
    if "context_borzoi" not in blocked.get("model_ids", ()) or blocked.get("eligible_after_gates") is not False:
        raise RowUniverseError("conditional architecture circularity gate differs")

    ask = config.get("bundled_cpu_validation_ask", {})
    if ask != {
        "job_name": "model-check-220",
        "partition": "cpu",
        "qos": "nslab",
        "cpus_per_task": 2,
        "memory_gb": 16,
        "walltime": "02:00:00",
        "array": False,
        "job_count": 1,
        "total_cpu_hours": 4,
        "total_gpu_hours": 0,
        "submitted": False,
    }:
        raise RowUniverseError("CPU validation ask differs")


def validate_authorities(root: Path, config: Mapping[str, Any]) -> dict[str, dict[str, str]]:
    files = config.get("frozen_file_authorities", {})
    trees = config.get("frozen_tree_authorities", {})
    if set(files) != FILE_AUTHORITIES or set(trees) != TREE_AUTHORITIES:
        raise RowUniverseError("frozen authority census differs")
    result: dict[str, dict[str, str]] = {}
    for name, binding in files.items():
        path = safe_relative_path(root, binding.get("path"), label=f"file authority {name}")
        expected = str(binding.get("sha256", ""))
        if not SHA256.fullmatch(expected) or sha256_file(path) != expected:
            raise RowUniverseError(f"frozen file authority changed: {name}")
        result[name] = {"kind": "file", "path": path.relative_to(root).as_posix(), "sha256": expected}
    for name, binding in trees.items():
        path = safe_relative_path(root, binding.get("path"), label=f"tree authority {name}")
        expected = str(binding.get("artifacts_sha256", ""))
        if not SHA256.fullmatch(expected):
            raise RowUniverseError(f"invalid frozen tree hash: {name}")
        try:
            verify_frozen_tree(path)
        except ArtifactError as error:
            raise RowUniverseError(f"frozen tree verification failed: {name}: {error}") from error
        if sha256_file(path / "ARTIFACTS.json") != expected:
            raise RowUniverseError(f"frozen tree authority changed: {name}")
        result[name] = {"kind": "tree", "path": path.relative_to(root).as_posix(), "artifacts_sha256": expected}
    return result


def model_records(path: Path) -> dict[str, dict[str, Any]]:
    document = load_toml(path, label=f"model registry {path.name}")
    rows = document.get("models", ())
    if not isinstance(rows, list):
        raise RowUniverseError(f"model registry {path.name} has no model roster")
    return {str(row.get("model_id")): dict(row) for row in rows}


def validate_registry_states(root: Path, authorities: Mapping[str, Mapping[str, str]]) -> dict[str, str]:
    records: dict[str, dict[str, Any]] = {}
    for authority in (
        "dna_language_registry",
        "regulatory_sequence_registry",
        "regulatory_local_registry",
    ):
        records.update(model_records(root / authorities[authority]["path"]))
    expected = {
        "dnabert2": ("candidate", "target_label_unexposed"),
        "hyenadna": ("candidate", "target_label_unexposed"),
        "caduceus": ("candidate", "target_label_unexposed"),
        "nucleotide_transformer": ("restricted_comparator", "target_label_unexposed"),
        "mpralegnet": ("candidate", "target_label_unexposed"),
        "sei": ("restricted_comparator", "target_label_unexposed"),
        "enformer": ("restricted_comparator", "target_label_unexposed"),
        "borzoi_ensemble": ("blocked_terms", "target_label_unexposed"),
        "corgi_regular": ("candidate", "target_label_unexposed"),
        "corgi_plus": ("deferred", "target_label_unexposed"),
        "context_borzoi": ("deferred", "unknown"),
    }
    for model_id, identity in expected.items():
        record = records.get(model_id)
        if record is None or (record.get("status"), record.get("exposure_status")) != identity:
            raise RowUniverseError(f"registered status or exposure differs: {model_id}")
    if "conflict" not in str(records["corgi_regular"].get("license_status", "")):
        raise RowUniverseError("Corgi license conflict is no longer represented")
    if "UNDECLARED" not in str(records["borzoi_ensemble"].get("license_status", "")):
        raise RowUniverseError("Borzoi weight terms are no longer represented")
    if "NC" not in str(records["nucleotide_transformer"].get("license_status", "")):
        raise RowUniverseError("Nucleotide Transformer restriction differs")
    return {model_id: f"{status}|{exposure}" for model_id, (status, exposure) in expected.items()}


def validate_prediction_identities(root: Path, authorities: Mapping[str, Mapping[str, str]]) -> dict[str, Any]:
    def tree(name: str) -> Path:
        return root / authorities[name]["path"]

    dna = load_json(tree("dna_lm_features") / "projected/receipt.json", label="DNA LM feature receipt")
    caduceus = load_json(tree("caduceus_features") / "projected/receipt.json", label="Caduceus feature receipt")
    mpr = load_json(tree("mpralegnet_predictions") / "predictions/receipt.json", label="MPRALegNet prediction receipt")
    sei = load_json(tree("sei_predictions") / "predictions/receipt.json", label="Sei prediction receipt")
    enformer = load_json(tree("enformer_predictions") / "predictions/receipt.json", label="Enformer prediction receipt")
    dna_eval = load_json(tree("dna_lm_one_seed_evaluation") / "evaluation/receipt.json", label="DNA LM evaluation receipt")
    caduceus_eval = load_json(tree("caduceus_one_seed_evaluation") / "evaluation/receipt.json", label="Caduceus evaluation receipt")
    borzoi = load_json(tree("borzoi_execution_preflight") / "receipt.json", label="Borzoi execution receipt")
    complementarity = load_json(tree("complementarity_readiness_audit") / "audit/receipt.json", label="complementarity receipt")

    if (
        dna.get("elements") != 1033
        or set(dna.get("models", {})) != {"dnabert2", "hyenadna", "nucleotide_transformer"}
        or dna.get("downstream_head_fit") is not False
        or dna.get("reporter_outcomes_read") is not False
        or caduceus.get("elements") != 1033
        or set(caduceus.get("models", {})) != {"caduceus"}
        or caduceus.get("downstream_head_fit") is not False
        or caduceus.get("reporter_outcomes_read") is not False
    ):
        raise RowUniverseError("outcome-blind sequence feature identities differ")
    if (
        mpr.get("elements") != 4359
        or mpr.get("outer_locus_sequence_groups") != 1033
        or mpr.get("native_output") != "uncalibrated_lentiMPRA_reporter_expression_score"
        or mpr.get("model_fitted_or_adapted") is not False
        or mpr.get("outcomes_read") is not False
        or sei.get("elements") != 1033
        or sei.get("model_fitted_or_adapted") is not False
        or sei.get("outcomes_read") is not False
        or enformer.get("elements") != 1033
        or enformer.get("model_fitted_or_adapted") is not False
        or enformer.get("outcomes_read") is not False
    ):
        raise RowUniverseError("native prediction identity or denominator differs")
    for receipt in (dna_eval, caduceus_eval):
        if (
            receipt.get("seed") != 20260824
            or receipt.get("contexts") != list(CONTEXTS)
            or receipt.get("elements") != 1033
            or receipt.get("outer_folds") != 5
            or receipt.get("status") != "pass_descriptive_one_seed_development_smoke"
            or receipt.get("sealed_outcomes_read") is not False
        ):
            raise RowUniverseError("one-seed common-lane evaluation identity differs")
    if (
        borzoi.get("executable") is not False
        or borzoi.get("model_forward_executed") is not False
        or borzoi.get("metrics_calculated") is not False
        or borzoi.get("sealed_assets_read") is not False
        or complementarity.get("bound_stacking_source_rows") != 0
        or complementarity.get("trigger_ready_pair_rows") != 0
        or complementarity.get("prediction_values_read") is not False
        or complementarity.get("raw_or_summary_outcomes_read") is not False
        or complementarity.get("conditional_model_build_authorized") is not False
    ):
        raise RowUniverseError("current complementarity or Borzoi boundary differs")
    return {
        "common_head_seed_ids": [20260824],
        "common_head_fold_map": "original_6kb_source_locus_group_map",
        "dna_lm_models": sorted(dna["models"]),
        "caduceus_models": sorted(caduceus["models"]),
        "mpralegnet_construct_rows": mpr["elements"],
        "native_selected_element_rows": sei["elements"],
        "borzoi_prediction_available": False,
        "context_prediction_authority_bound": False,
    }


def build_rows(root: Path, authorities: Mapping[str, Mapping[str, str]], config: Mapping[str, Any]) -> tuple[list[dict[str, object]], dict[str, Any]]:
    split = root / authorities["outcome_blind_split"]["path"]
    selection = root / authorities["common_element_selection"]["path"]
    long_range = root / authorities["long_range_fold_map"]["path"]
    groups = read_tsv(split / "split/groups.tsv", GROUP_FIELDS, label="source group table")
    selected = read_tsv(selection / "fixture/manifest.tsv", SEI_FIELDS, label="selected element table")
    block_map = read_tsv(long_range / "fixture/source_group_map.tsv", BORZOI_MAP_FIELDS, label="long-range source map")
    borzoi_manifest = read_tsv(long_range / "fixture/manifest.tsv", BORZOI_MANIFEST_FIELDS, label="Borzoi selected element table")
    components = read_tsv(long_range / "fixture/long_range_components.tsv", BORZOI_COMPONENT_FIELDS, label="long-range component table")

    if len(groups) != 1033 or len(selected) != 1033 or len(block_map) != 1033 or len(borzoi_manifest) != 1033 or len(components) != 239:
        raise RowUniverseError("source locus, selection, or long-range block count differs")
    original = {row["outer_locus_sequence_group_id"]: row for row in groups}
    chosen = {row["outer_locus_sequence_group_id"]: row for row in selected}
    mapping = {row["source_locus_group_id"]: row for row in block_map}
    borzoi_chosen = {row["source_locus_group_id"]: row for row in borzoi_manifest}
    if len(original) != 1033 or set(original) != set(chosen) or set(original) != set(mapping) or set(original) != set(borzoi_chosen):
        raise RowUniverseError("source locus group sets differ")
    if any(chosen[group]["element_id"] != borzoi_chosen[group]["element_id"] for group in original):
        raise RowUniverseError("long-range and common fixtures select different elements")
    if any(mapping[group]["contig"] != original[group]["contig"] for group in original):
        raise RowUniverseError("long-range source contig differs")
    if any(mapping[group]["borzoi_long_range_group_id"] != borzoi_chosen[group]["borzoi_long_range_group_id"] for group in original):
        raise RowUniverseError("long-range block identity differs between fixtures")
    if any(mapping[group]["outer_fold"] != borzoi_chosen[group]["outer_fold"] for group in original):
        raise RowUniverseError("long-range fold differs between fixtures")

    component_by_id = {row["borzoi_long_range_group_id"]: row for row in components}
    if len(component_by_id) != 239 or set(component_by_id) != {row["borzoi_long_range_group_id"] for row in block_map}:
        raise RowUniverseError("long-range block roster differs")
    by_contig: dict[str, list[dict[str, str]]] = {}
    for component in components:
        by_contig.setdefault(component["contig"], []).append(component)
    buffer_bp = int(config["fold_contract"]["largest_receptive_field_buffer_bp"])
    for contig_rows in by_contig.values():
        ordered = sorted(contig_rows, key=lambda row: int(row["minimum_variant_pos0"]))
        for left, right in zip(ordered, ordered[1:]):
            distance = int(right["minimum_variant_pos0"]) - int(left["maximum_variant_pos0"])
            if distance < buffer_bp:
                raise RowUniverseError("long-range blocks violate the 524288-bp boundary")

    fold_counts = Counter(int(row["outer_fold"]) for row in block_map)
    expected_fold_counts = {
        int(name.removeprefix("fold-")): int(value)
        for name, value in config["fold_contract"]["expected_fold_source_locus_groups"].items()
    }
    if dict(sorted(fold_counts.items())) != expected_fold_counts:
        raise RowUniverseError("long-range fold source-locus counts differ")
    reassigned = sum(int(original[group]["outer_fold"]) != int(mapping[group]["outer_fold"]) for group in original)
    if reassigned != config["fold_contract"]["expected_source_groups_reassigned_from_original_6kb_fold_map"]:
        raise RowUniverseError("reassignment count from the original fold map differs")

    rows: list[dict[str, object]] = []
    for seed in SEEDS:
        for group in sorted(original):
            element_id = chosen[group]["element_id"]
            block_id = mapping[group]["borzoi_long_range_group_id"]
            outer_fold = f"fold-{mapping[group]['outer_fold']}"
            for context in CONTEXTS:
                rows.append(
                    {
                        "seed": seed,
                        "row_hash": hash_id("gse281364", "mpra_allelic_direction", element_id, group, context),
                        "unit_hash": hash_id("gse281364", "outer_locus_sequence_group", group),
                        "block_hash": hash_id("gse281364", "long_range_524288bp_component", block_id),
                        "stratum": "all",
                        "outer_fold": outer_fold,
                        "study_id": "gse281364",
                        "assay_context_id": context,
                        "element_id": element_id,
                        "source_locus_group_id": group,
                        "long_range_block_id": block_id,
                    }
                )
    if len(rows) != 10330 or len({(row["seed"], row["row_hash"]) for row in rows}) != 10330:
        raise RowUniverseError("seeded row count or identity uniqueness differs")
    row_sets = {
        seed: {str(row["row_hash"]) for row in rows if row["seed"] == seed}
        for seed in SEEDS
    }
    if any(values != row_sets[SEEDS[0]] for values in row_sets.values()):
        raise RowUniverseError("row sets differ across seeds")
    if len(row_sets[SEEDS[0]]) != 2066:
        raise RowUniverseError("base row set count differs")
    metadata = {
        "source_locus_groups": 1033,
        "selected_elements": 1033,
        "assay_contexts": list(CONTEXTS),
        "base_rows_per_seed": 2066,
        "fixed_seeds": list(SEEDS),
        "seeded_rows": 10330,
        "long_range_blocks": 239,
        "fold_source_locus_groups": {f"fold-{fold}": count for fold, count in sorted(fold_counts.items())},
        "source_groups_reassigned_from_original_6kb_fold_map": reassigned,
        "base_row_set_sha256": canonical_set_hash("gse281364_base_row_hashes", row_sets[SEEDS[0]]),
        "unit_set_sha256": canonical_set_hash("gse281364_unit_hashes", {str(row["unit_hash"]) for row in rows}),
        "block_set_sha256": canonical_set_hash("gse281364_block_hashes", {str(row["block_hash"]) for row in rows}),
        "selected_element_set_sha256": canonical_set_hash("gse281364_element_ids", {str(row["element_id"]) for row in rows}),
    }
    return rows, metadata


def source_rows(config: Mapping[str, Any]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for source in config["source_status"]:
        rows.append(
            {
                "source_id": source["source_id"],
                "source_family_role": source["source_family_role"],
                "model_ids": ",".join(source["model_ids"]),
                "current_row_count_per_context": source["current_row_count_per_context"],
                "current_seed_ids": ",".join(str(seed) for seed in source["current_seed_ids"]),
                "current_fold_map": source["current_fold_map"],
                "current_scale": source["current_scale"],
                "current_open_trigger_eligible": str(source["current_open_trigger_eligible"]).lower(),
                "eligible_after_gates": str(source["eligible_after_gates"]).lower(),
                "diagnostic_only": str(source["diagnostic_only"]).lower(),
                "blockers": "|".join(source["blockers"]),
            }
        )
    return rows


def preflight(*, root: Path, config_path: Path, output: Path | None = None) -> dict[str, Any]:
    root = root.resolve(strict=True)
    config_path = config_path.resolve(strict=True)
    try:
        config_path.relative_to(root)
    except ValueError as error:
        raise RowUniverseError("config must be inside project root") from error
    config = load_toml(config_path, label="row-universe config")
    validate_config(config)
    authorities = validate_authorities(root, config)
    registry_states = validate_registry_states(root, authorities)
    prediction_identities = validate_prediction_identities(root, authorities)
    rows, metadata = build_rows(root, authorities, config)
    sources = source_rows(config)
    row_text = tsv_text(ROW_FIELDS, rows)
    source_text = tsv_text(SOURCE_FIELDS, sources)
    receipt: dict[str, Any] = {
        "schema_version": SCHEMA,
        "status": STATUS,
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "endpoint_id": "mpra_allelic_direction",
        "evaluator_id": "locus_heldout_allelic_spearman_v1",
        **metadata,
        "row_universe_tsv_sha256": hashlib.sha256(row_text.encode("utf-8")).hexdigest(),
        "source_status_tsv_sha256": hashlib.sha256(source_text.encode("utf-8")).hexdigest(),
        "source_status_rows": len(sources),
        "current_open_sequence_source_rows": 0,
        "current_open_context_source_rows": 0,
        "current_trigger_compatible_pair_rows": 0,
        "prediction_identities": prediction_identities,
        "registry_states": registry_states,
        "authority_count": len(authorities),
        "authorities": authorities,
        "outcomes_read": False,
        "prediction_values_read": False,
        "reporter_counts_read": False,
        "sealed_assets_read": False,
        "predictions_generated": False,
        "stack_fit": False,
        "residual_correlation_calculated": False,
        "conditional_model_built_or_fit": False,
        "stacking_blockers": [
            "no_admitted_context_prediction_source_on_the_exact_row_universe",
            "no_sequence_or_context_source_has_all_five_fixed_seeds",
            "current_cross_fitted_heads_use_the_superseded_original_6kb_fold_map",
            "static_native_scores_are_not_calibrated_to_the_required_assay_native_outcome_unit",
            "open_sequence_and_context_development_shortlists_are_unbound",
        ],
    }
    if output is not None:
        if output.exists():
            raise RowUniverseError(f"refusing to overwrite output: {output}")
        output.mkdir(parents=True)
        write_text_exclusive(output / "row_universe.tsv", row_text)
        write_text_exclusive(output / "source_status.tsv", source_text)
        write_json_exclusive(output / "receipt.json", receipt)
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = preflight(root=args.root, config_path=args.config, output=args.output)
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

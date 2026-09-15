#!/usr/bin/env python3
"""Evaluate the exact five-fold by five-seed scBasset development rectangle."""

from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from statistics import fmean
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import verify_frozen_tree
from scripts import evaluate_scbasset_valid as base
from scripts.evaluate_sequence_task_native_five_seed_rectangle import (
    FiveSeedEvaluationError,
    stratified_two_way_bootstrap,
)


FOLDS = tuple(range(5))
SEEDS = (20260824, 20260825, 20260826, 20260827, 20260828)
LINEAGES = tuple(base.LINEAGES)
COMPLETION_FIELDS = (
    "outer_fold",
    "split_id",
    "seed",
    "disposition",
    "model_root",
    "model_artifacts_sha256",
    "prediction_root",
    "prediction_artifacts_sha256",
)
UNIT_FIELDS = (
    "outer_fold",
    "model_id",
    "seed",
    "donor_hash",
    "block_hash",
    "lineage_id",
    "observed_insertions",
    "regions",
    "deviance_per_insertion",
)
SECONDARY_FIELDS = (
    "outer_fold",
    "model_id",
    "seed",
    "donor_hash",
    "lineage_id",
    "peak_auprc",
    "profile_spearman",
)


class ScBassetFiveSeedEvaluationError(RuntimeError):
    """Raised when exact coverage, source integrity, or donor safety differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _json(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ScBassetFiveSeedEvaluationError(f"missing or linked {label}: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScBassetFiveSeedEvaluationError(f"invalid {label}: {path}") from error
    if not isinstance(value, dict):
        raise ScBassetFiveSeedEvaluationError(f"{label} must be an object")
    return value


def _read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    if path.is_symlink() or not path.is_file():
        raise ScBassetFiveSeedEvaluationError(f"missing or linked TSV: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise ScBassetFiveSeedEvaluationError(f"TSV fields differ: {path}")
        return [dict(row) for row in reader]


def _write_tsv(
    path: Path,
    fields: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def _resolve(root: Path, value: str) -> Path:
    configured = Path(value)
    path = configured if configured.is_absolute() else root / configured
    path = path.resolve(strict=True)
    if root not in path.parents:
        raise ScBassetFiveSeedEvaluationError(f"bound path escapes benchmark root: {path}")
    return path


def _verify_root(root: Path, expected: str, label: str) -> dict[str, Any]:
    if digest(root / "ARTIFACTS.json") != expected:
        raise ScBassetFiveSeedEvaluationError(f"{label} ARTIFACTS differs: {root}")
    verify_frozen_tree(root)
    return _json(root / "ARTIFACTS.json", f"{label} ARTIFACTS").get("metadata", {})


def load_evaluator_contract(path: Path, root: Path) -> dict[str, Any]:
    path = path.resolve(strict=True)
    if root not in path.parents:
        raise ScBassetFiveSeedEvaluationError("evaluator contract escapes benchmark root")
    contract = _json(path, "scBasset evaluator contract")
    rectangle = contract.get("rectangle", {})
    firewall = contract.get("firewall", {})
    claims = contract.get("claims", {})
    if (
        contract.get("schema_version")
        != "masld-bench-scbasset-five-seed-evaluator-contract-v1"
        or contract.get("dataset_id") != "gse296875"
        or contract.get("evaluation_role") != "valid"
        or rectangle.get("outer_folds") != list(FOLDS)
        or rectangle.get("fixed_seeds") != list(SEEDS)
        or rectangle.get("expected_fits") != 25
        or rectangle.get("expected_prediction_views") != 25
        or rectangle.get("model_id") != "scbasset"
        or firewall.get("all_25_prediction_trees_verify_before_outcome_authorities_resolve")
        is not True
        or firewall.get("partial_rectangle_ranking_allowed") is not False
        or firewall.get("test_role_authorized") is not False
        or firewall.get("sealed_data_authorized") is not False
        or claims.get("rna_conditioned_claim_allowed") is not False
        or claims.get("cross_task_family_ranking_allowed") is not False
        or claims.get("champion_claim_allowed") is not False
    ):
        raise ScBassetFiveSeedEvaluationError("scBasset evaluator contract differs")
    inputs = contract.get("input_artifacts", [])
    if (
        len(inputs) != 5
        or {int(row["outer_fold"]) for row in inputs} != set(FOLDS)
        or any(row["split_id"] != f"donor{row['outer_fold']}_genomic{row['outer_fold']}" for row in inputs)
    ):
        raise ScBassetFiveSeedEvaluationError("scBasset input authority roster differs")
    for section in contract.get("outcome_authorities", {}).values():
        if not isinstance(section, dict) or not section.get("root") or not section.get(
            "artifacts_sha256"
        ):
            raise ScBassetFiveSeedEvaluationError("outcome authority binding differs")
    return contract


def _prediction_bundle_record(
    prediction_root: Path,
    *,
    fold: int,
    seed: int,
    input_artifacts_sha256: str,
    model_artifacts_sha256: str,
) -> dict[str, Any]:
    bundle_path = prediction_root / "predictions/bundle/prediction_bundle.json"
    bundle = _json(bundle_path, "prediction bundle")
    metadata = bundle.get("metadata", {})
    donor_count = int(metadata.get("held_donor_count", -1))
    if (
        bundle.get("task_id") != "rna_conditioned_atac"
        or bundle.get("model_id") != "scbasset"
        or donor_count <= 0
        or bundle.get("n_predictions") != donor_count * len(LINEAGES) * 16_000
        or metadata.get("donor_test_fold") != fold
        or metadata.get("donor_valid_fold") != (fold + 1) % 5
        or metadata.get("evaluation_role") != "valid"
        or metadata.get("region_count") != 16_000
        or metadata.get("held_atac_input_exposed") is not False
        or metadata.get("observed_atac_exported") is not False
        or metadata.get("rna_input_exposed") is not False
        or metadata.get("held_donor_profiles_identical_within_lineage") is not True
        or metadata.get("model_artifact_sha256") != model_artifacts_sha256
    ):
        raise ScBassetFiveSeedEvaluationError("prediction bundle contract differs")
    record = bundle.get("standardized_table", {})
    table = bundle_path.parent / str(record.get("path", ""))
    artifacts = _json(prediction_root / "ARTIFACTS.json", "prediction ARTIFACTS")
    artifact_by_path = {
        row.get("path"): row for row in artifacts.get("artifacts", []) if isinstance(row, dict)
    }
    frozen_table = artifact_by_path.get("predictions/bundle/predictions.tsv", {})
    if (
        table.is_symlink()
        or not table.is_file()
        or table.stat().st_size != record.get("size_bytes")
        or frozen_table.get("size_bytes") != record.get("size_bytes")
        or frozen_table.get("sha256") != record.get("sha256")
    ):
        raise ScBassetFiveSeedEvaluationError("prediction table binding differs")
    return {
        "fold": fold,
        "seed": seed,
        "table": table,
        "table_sha256": record["sha256"],
        "expected_rows": int(bundle["n_predictions"]),
        "donor_count": donor_count,
        "namespace": str(bundle.get("unit_id_namespace", "")),
        "input_artifacts_sha256": input_artifacts_sha256,
        "model_artifacts_sha256": model_artifacts_sha256,
        "prediction_root": prediction_root,
    }


def load_locked_sources(
    *,
    root: Path,
    evaluator_contract: Path,
    completion_root: Path,
    completion_artifacts_sha256: str,
) -> tuple[dict[str, Any], dict[int, dict[str, Any]], dict[int, list[dict[str, Any]]]]:
    """Verify all model-side output files without resolving an outcome path."""

    contract = load_evaluator_contract(evaluator_contract, root)
    completion_root = completion_root.resolve(strict=True)
    if root not in completion_root.parents:
        raise ScBassetFiveSeedEvaluationError("completion lock escapes benchmark root")
    _verify_root(completion_root, completion_artifacts_sha256, "completion lock")
    completion = _json(
        completion_root / "completion/completion_lock.json", "completion receipt"
    )
    if (
        completion.get("status") != "ready_both_rectangles_locked"
        or completion.get("both_full_rectangles") is not True
        or completion.get("evaluation_launch_authorized") is not True
        or completion.get("partial_ranking_authorized") is not False
        or completion.get("prediction_values_read") is not False
        or completion.get("observed_outcomes_read") is not False
    ):
        raise ScBassetFiveSeedEvaluationError("completion lock does not authorize evaluation")
    completion_contract_path = (
        root
        / "config/evaluation/sequence_regulatory_full_rectangle_completion_20260825.json"
    ).resolve(strict=True)
    completion_contract = _json(completion_contract_path, "completion contract")
    binding = completion_contract.get("task_families", {}).get(
        "scbasset_sequence_only", {}
    )
    evaluator_contract = evaluator_contract.resolve(strict=True)
    if (
        completion.get("contract_sha256") != digest(completion_contract_path)
        or _resolve(root, str(binding.get("evaluator_contract", "")))
        != evaluator_contract
        or binding.get("evaluator_contract_sha256") != digest(evaluator_contract)
    ):
        raise ScBassetFiveSeedEvaluationError("completion-to-evaluator binding differs")
    rows = _read_tsv(
        completion_root / "completion/scbasset_fit_dispositions.tsv", COMPLETION_FIELDS
    )
    expected = {(fold, seed) for fold in FOLDS for seed in SEEDS}
    observed = {(int(row["outer_fold"]), int(row["seed"])) for row in rows}
    if (
        len(rows) != 25
        or observed != expected
        or any(row["disposition"] != "complete_compatible" for row in rows)
    ):
        raise ScBassetFiveSeedEvaluationError("scBasset completion rectangle is not exact")

    input_sources: dict[int, dict[str, Any]] = {}
    for binding in contract["input_artifacts"]:
        fold = int(binding["outer_fold"])
        input_root = _resolve(root, binding["root"])
        metadata = _verify_root(
            input_root, binding["artifacts_sha256"], f"fold {fold} input"
        )
        summary = _json(input_root / "inputs/summary.json", "input summary")
        if (
            summary.get("status") != "pass"
            or summary.get("dataset_id") != "gse296875"
            or summary.get("split_id") != binding["split_id"]
            or summary.get("donor_test_fold") != fold
            or summary.get("donor_valid_fold") != (fold + 1) % 5
            or summary.get("held_donor_atac_exported") is not False
            or summary.get("genomic_test_atac_exported") is not False
        ):
            raise ScBassetFiveSeedEvaluationError(f"fold {fold} input contract differs")
        input_sources[fold] = {
            "root": input_root,
            "artifacts_sha256": binding["artifacts_sha256"],
            "metadata": metadata,
            "summary": summary,
        }

    predictions: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in sorted(rows, key=lambda item: (int(item["outer_fold"]), int(item["seed"]))):
        fold = int(row["outer_fold"])
        seed = int(row["seed"])
        split_id = f"donor{fold}_genomic{fold}"
        if row["split_id"] != split_id:
            raise ScBassetFiveSeedEvaluationError("scBasset split identity differs")
        model_root = _resolve(root, row["model_root"])
        prediction_root = _resolve(root, row["prediction_root"])
        model_metadata = _verify_root(
            model_root, row["model_artifacts_sha256"], f"fold {fold} seed {seed} model"
        )
        prediction_metadata = _verify_root(
            prediction_root,
            row["prediction_artifacts_sha256"],
            f"fold {fold} seed {seed} predictions",
        )
        input_sha = input_sources[fold]["artifacts_sha256"]
        if (
            model_metadata.get("artifact_class") != "scbasset_fold_model"
            or prediction_metadata.get("artifact_class") != "scbasset_fold_predictions"
            or model_metadata.get("dataset_id") != "gse296875"
            or prediction_metadata.get("dataset_id") != "gse296875"
            or model_metadata.get("model_id") != "scbasset"
            or prediction_metadata.get("model_id") != "scbasset"
            or model_metadata.get("split_id") != split_id
            or prediction_metadata.get("split_id") != split_id
            or model_metadata.get("seed") != seed
            or prediction_metadata.get("seed") != seed
            or model_metadata.get("input_artifacts_sha256") != input_sha
            or prediction_metadata.get("input_artifacts_sha256") != input_sha
            or prediction_metadata.get("model_artifacts_sha256")
            != row["model_artifacts_sha256"]
            or prediction_metadata.get("evaluation_role") != "valid"
            or prediction_metadata.get("held_donor_atac_used") is not False
            or prediction_metadata.get("held_cell_embedding_available") is not False
        ):
            raise ScBassetFiveSeedEvaluationError("scBasset source metadata differs")
        predictions[fold].append(
            _prediction_bundle_record(
                prediction_root,
                fold=fold,
                seed=seed,
                input_artifacts_sha256=input_sha,
                model_artifacts_sha256=row["model_artifacts_sha256"],
            )
        )
    for fold in FOLDS:
        predictions[fold].sort(key=lambda item: item["seed"])
        if tuple(item["seed"] for item in predictions[fold]) != SEEDS:
            raise ScBassetFiveSeedEvaluationError("scBasset fold seed roster differs")
        if len({item["namespace"] for item in predictions[fold]}) != 1:
            raise ScBassetFiveSeedEvaluationError("scBasset seed namespaces differ")
        if any(
            item["donor_count"] != predictions[fold][0]["donor_count"]
            or item["expected_rows"] != predictions[fold][0]["expected_rows"]
            for item in predictions[fold]
        ):
            raise ScBassetFiveSeedEvaluationError("scBasset seed axes differ")
    return contract, input_sources, dict(predictions)


def _verify_outcome_authorities(
    root: Path, contract: Mapping[str, Any]
) -> tuple[Path, Path]:
    authorities = contract["outcome_authorities"]
    bigwigs = _resolve(root, authorities["donor_bigwigs"]["root"])
    split = _resolve(root, authorities["sequence_split"]["root"])
    _verify_root(
        bigwigs,
        authorities["donor_bigwigs"]["artifacts_sha256"],
        "donor bigWigs",
    )
    _verify_root(
        split,
        authorities["sequence_split"]["artifacts_sha256"],
        "sequence split",
    )
    return bigwigs, split


def _prediction_index(
    *, namespace: str, donor_ids: Sequence[str], regions: Sequence[Mapping[str, str]]
) -> dict[str, int]:
    expected: dict[str, int] = {}
    for donor_offset, donor in enumerate(donor_ids):
        for lineage_offset, lineage in enumerate(LINEAGES):
            for region_offset, region in enumerate(regions):
                row_hash = base.join_hash(
                    namespace,
                    "row",
                    f"{donor}\0{lineage}\0{region['region_id']}",
                )
                flat = (
                    (donor_offset * len(LINEAGES) + lineage_offset) * len(regions)
                    + region_offset
                )
                if row_hash in expected:
                    raise ScBassetFiveSeedEvaluationError("prediction row hash collides")
                expected[row_hash] = flat
    return expected


def load_prediction_array(
    source: Mapping[str, Any],
    *,
    donor_ids: Sequence[str],
    regions: Sequence[Mapping[str, str]],
    expected: Mapping[str, int],
) -> Any:
    """Stream one standardized prediction table into the common raw-axis tensor."""

    import numpy as np

    shape = (len(donor_ids), len(LINEAGES), len(regions))
    values = np.empty(math.prod(shape), dtype=np.float64)
    seen = np.zeros(values.size, dtype=np.bool_)
    table = Path(source["table"])
    count = 0
    with table.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(base.PREDICTION_FIELDS):
            raise ScBassetFiveSeedEvaluationError("prediction table fields differ")
        for row in reader:
            offset = expected.get(row["row_hash"])
            if offset is None or seen[offset]:
                raise ScBassetFiveSeedEvaluationError("prediction row identity differs")
            donor_offset, remainder = divmod(offset, len(LINEAGES) * len(regions))
            lineage_offset, region_offset = divmod(remainder, len(regions))
            donor = donor_ids[donor_offset]
            lineage = LINEAGES[lineage_offset]
            region = regions[region_offset]
            namespace = source["namespace"]
            if (
                row["donor_hash"] != base.join_hash(namespace, "unit", donor)
                or row["block_hash"]
                != base.join_hash(namespace, "block", region["block_id"])
                or row["stratum"] != lineage
            ):
                raise ScBassetFiveSeedEvaluationError("prediction join fields differ")
            value = float(row["predicted"])
            if not math.isfinite(value) or value <= 0:
                raise ScBassetFiveSeedEvaluationError("prediction value differs")
            values[offset] = value
            seen[offset] = True
            count += 1
    if count != source["expected_rows"] or not bool(seen.all()):
        raise ScBassetFiveSeedEvaluationError("prediction universe is incomplete")
    array = values.reshape(shape)
    if not np.allclose(array, array[:1], rtol=0, atol=1.0e-12):
        raise ScBassetFiveSeedEvaluationError("sequence-only held-donor profiles differ")
    return array


def _metric_rows(
    *,
    fold: int,
    model_id: str,
    seed: str,
    vector: Any,
    observed: Any,
    donor_ids: Sequence[str],
    regions: Sequence[Mapping[str, str]],
    namespace: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    import numpy as np

    block_positions = {
        block: np.asarray(
            [index for index, row in enumerate(regions) if row["block_id"] == block],
            dtype=np.int64,
        )
        for block in sorted({row["block_id"] for row in regions})
    }
    unit_rows: list[dict[str, object]] = []
    secondary_rows: list[dict[str, object]] = []
    for donor_offset, donor in enumerate(donor_ids):
        for lineage_offset, lineage in enumerate(LINEAGES):
            truth = observed[donor_offset, lineage_offset].astype(np.float64)
            candidate = vector[donor_offset, lineage_offset]
            donor_hash = base.join_hash(namespace, "donor", donor)
            secondary_rows.append(
                {
                    "outer_fold": fold,
                    "model_id": model_id,
                    "seed": seed,
                    "donor_hash": donor_hash,
                    "lineage_id": lineage,
                    "peak_auprc": format(base.average_precision(truth, candidate), ".17g"),
                    "profile_spearman": format(
                        base.spearman_or_zero(truth, candidate), ".17g"
                    ),
                }
            )
            for block, positions in block_positions.items():
                block_truth = truth[positions]
                if block_truth.sum() <= 0:
                    continue
                unit_rows.append(
                    {
                        "outer_fold": fold,
                        "model_id": model_id,
                        "seed": seed,
                        "donor_hash": donor_hash,
                        "block_hash": base.join_hash(
                            namespace, "block", f"fold{fold}:{block}"
                        ),
                        "lineage_id": lineage,
                        "observed_insertions": int(block_truth.sum()),
                        "regions": len(positions),
                        "deviance_per_insertion": format(
                            base.deviance_per_insertion(block_truth, candidate[positions]),
                            ".17g",
                        ),
                    }
                )
    return unit_rows, secondary_rows


def evaluate_fold(
    *,
    fold: int,
    input_source: Mapping[str, Any],
    prediction_sources: Sequence[Mapping[str, Any]],
    donor_bigwigs: Path,
    split_root: Path,
    workers: int,
    namespace: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]], set[str]]:
    import numpy as np
    from scipy import sparse

    input_root = Path(input_source["root"])
    inputs = input_root / "inputs"
    summary = input_source["summary"]
    donors = _read_tsv(inputs / "held_valid_donors.tsv", base.DONOR_FIELDS)
    regions = _read_tsv(inputs / "valid_regions.tsv", base.REGION_FIELDS)
    cells = _read_tsv(inputs / "training_cells.tsv", base.CELL_FIELDS)
    donor_ids = [row["donor_id"] for row in donors]
    valid_fold = (fold + 1) % 5
    if (
        len(donor_ids) != prediction_sources[0]["donor_count"]
        or len(set(donor_ids)) != len(donor_ids)
        or len(regions) != 16_000
        or len(cells) != int(summary["training_cells"])
        or any(int(row["outer_fold"]) != valid_fold for row in donors)
        or any(row["evaluation_role"] != "valid" for row in donors)
        or any(row["role"] != "valid" for row in regions)
    ):
        raise ScBassetFiveSeedEvaluationError(f"fold {fold} validation axes differ")

    split_windows = {
        row["window_id"]: row
        for row in _read_tsv(
            split_root / "ccre_evaluation_windows.tsv", base.CCRE_FIELDS
        )
        if int(row["genomic_fold"]) == valid_fold
    }
    if set(split_windows) != {row["region_id"] for row in regions}:
        raise ScBassetFiveSeedEvaluationError(f"fold {fold} cCRE axes differ")
    for row in regions:
        window = split_windows[row["region_id"]]
        midpoint = (int(window["output_start"]) + int(window["output_end"])) // 2
        if (
            row["block_id"] != window["contig"]
            or int(row["sequence_start"]) != midpoint - 672
            or int(row["sequence_end"]) != midpoint + 672
        ):
            raise ScBassetFiveSeedEvaluationError(f"fold {fold} sequence join differs")

    donor_index = {value: index for index, value in enumerate(donor_ids)}
    lineage_index = {value: index for index, value in enumerate(LINEAGES)}
    manifest = _read_tsv(donor_bigwigs / "bigwig_manifest.tsv", base.BIGWIG_FIELDS)
    manifest_index = {
        (row["donor_id"], row["lineage_id"]): row
        for row in manifest
        if row["donor_id"] in donor_index and row["lineage_id"] in lineage_index
    }
    if len(manifest_index) != len(donor_ids) * len(LINEAGES):
        raise ScBassetFiveSeedEvaluationError(f"fold {fold} bigWig census differs")
    windows = [
        (
            split_windows[row["region_id"]]["contig"],
            int(split_windows[row["region_id"]]["output_start"]),
            int(split_windows[row["region_id"]]["output_end"]),
        )
        for row in regions
    ]
    extraction = [
        {
            **row,
            "lineage_id": lineage,
            "path": str(donor_bigwigs / row["path"]),
            "windows": windows,
        }
        for (_, lineage), row in sorted(manifest_index.items())
    ]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        extracted = list(executor.map(base._extract_bigwig, extraction))
    observed = np.zeros((len(donor_ids), len(LINEAGES), len(regions)), dtype=np.uint32)
    for item in extracted:
        observed[
            donor_index[item["donor_id"]], lineage_index[item["lineage"]]
        ] = item["counts"]
    if np.any(observed.sum(axis=2) <= 0):
        raise ScBassetFiveSeedEvaluationError(f"fold {fold} held outcome is empty")

    training = sparse.load_npz(inputs / "m_valid.npz").tocsr()
    if training.shape != (len(regions), len(cells)):
        raise ScBassetFiveSeedEvaluationError(f"fold {fold} training matrix differs")
    training_labels = np.asarray([row["lineage"] for row in cells])
    lineage_mean = np.empty((len(LINEAGES), len(regions)), dtype=np.float64)
    for lineage, offset in lineage_index.items():
        positions = np.flatnonzero(training_labels == lineage)
        if positions.size == 0:
            raise ScBassetFiveSeedEvaluationError(f"fold {fold} lineage is empty")
        lineage_mean[offset] = np.asarray(training[:, positions].mean(axis=1)).ravel()
    global_mean = np.asarray(training.mean(axis=1)).ravel()
    lineage_mean += 1.0e-8
    global_mean += 1.0e-8
    lineage_mean /= lineage_mean.sum(axis=1, keepdims=True)
    global_mean /= global_mean.sum()
    baseline_lineage = np.broadcast_to(
        lineage_mean[None, :, :], observed.shape
    )
    baseline_global = np.broadcast_to(global_mean[None, None, :], observed.shape)

    units: list[dict[str, object]] = []
    secondary: list[dict[str, object]] = []
    for model_id, vector in (
        ("training_lineage_mean", baseline_lineage),
        ("training_global_mean", baseline_global),
    ):
        these_units, these_secondary = _metric_rows(
            fold=fold,
            model_id=model_id,
            seed="training_only",
            vector=vector,
            observed=observed,
            donor_ids=donor_ids,
            regions=regions,
            namespace=namespace,
        )
        units.extend(these_units)
        secondary.extend(these_secondary)

    source_namespace = prediction_sources[0]["namespace"]
    expected = _prediction_index(
        namespace=source_namespace, donor_ids=donor_ids, regions=regions
    )
    ensemble = np.zeros(observed.shape, dtype=np.float64)
    for source in prediction_sources:
        candidate = load_prediction_array(
            source, donor_ids=donor_ids, regions=regions, expected=expected
        )
        ensemble += candidate / len(SEEDS)
        these_units, these_secondary = _metric_rows(
            fold=fold,
            model_id="scbasset",
            seed=str(source["seed"]),
            vector=candidate,
            observed=observed,
            donor_ids=donor_ids,
            regions=regions,
            namespace=namespace,
        )
        units.extend(these_units)
        secondary.extend(these_secondary)
    if np.any(~np.isfinite(ensemble)) or np.any(ensemble <= 0):
        raise ScBassetFiveSeedEvaluationError("five-seed ensemble differs")
    these_units, these_secondary = _metric_rows(
        fold=fold,
        model_id="scbasset",
        seed="ensemble",
        vector=ensemble,
        observed=observed,
        donor_ids=donor_ids,
        regions=regions,
        namespace=namespace,
    )
    units.extend(these_units)
    secondary.extend(these_secondary)
    return units, secondary, set(donor_ids)


def _macro_mean(rows: Sequence[Mapping[str, Any]], field: str) -> float:
    strata: dict[tuple[int, str], list[float]] = defaultdict(list)
    for row in rows:
        strata[(int(row["outer_fold"]), str(row["lineage_id"]))].append(
            float(row[field])
        )
    if set(strata) != {(fold, lineage) for fold in FOLDS for lineage in LINEAGES}:
        raise ScBassetFiveSeedEvaluationError("macro stratum roster differs")
    return fmean(fmean(values) for values in strata.values())


def _paired_skill_rows(
    candidate: Sequence[Mapping[str, Any]],
    baseline: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    key_fields = ("outer_fold", "donor_hash", "block_hash", "lineage_id")
    baseline_index = {
        tuple(row[field] for field in key_fields): float(row["deviance_per_insertion"])
        for row in baseline
    }
    if len(baseline_index) != len(baseline):
        raise ScBassetFiveSeedEvaluationError("baseline unit identity is duplicated")
    result = []
    for row in candidate:
        key = tuple(row[field] for field in key_fields)
        value = float(row["deviance_per_insertion"])
        reference = baseline_index.pop(key, None)
        if reference is None or reference <= 0 or not math.isfinite(value):
            raise ScBassetFiveSeedEvaluationError("paired deviance unit differs")
        result.append({**row, "skill": (reference - value) / reference})
    if baseline_index:
        raise ScBassetFiveSeedEvaluationError("paired deviance universe differs")
    return result


def summarize(
    unit_rows: Sequence[Mapping[str, Any]],
    secondary_rows: Sequence[Mapping[str, Any]],
    *,
    n_resamples: int,
    bootstrap_seed: int,
) -> dict[str, Any]:
    baselines = {}
    for model_id in ("training_lineage_mean", "training_global_mean"):
        rows = [row for row in unit_rows if row["model_id"] == model_id]
        baselines[model_id] = {
            "rows": rows,
            "macro_deviance_per_insertion": _macro_mean(rows, "deviance_per_insertion"),
        }
    strongest = min(
        baselines,
        key=lambda model: (baselines[model]["macro_deviance_per_insertion"], model),
    )
    baseline_rows = baselines[strongest]["rows"]
    ensemble = [
        row
        for row in unit_rows
        if row["model_id"] == "scbasset" and row["seed"] == "ensemble"
    ]
    skills = _paired_skill_rows(ensemble, baseline_rows)
    try:
        uncertainty = stratified_two_way_bootstrap(
            skills, n_resamples=n_resamples, seed=bootstrap_seed
        )
    except FiveSeedEvaluationError as error:
        raise ScBassetFiveSeedEvaluationError(str(error)) from error
    seed_scores = {}
    for seed in SEEDS:
        rows = [
            row
            for row in unit_rows
            if row["model_id"] == "scbasset" and row["seed"] == str(seed)
        ]
        seed_scores[str(seed)] = _macro_mean(
            _paired_skill_rows(rows, baseline_rows), "skill"
        )
    lineage_scores = {
        lineage: fmean(
            fmean(
                float(row["skill"])
                for row in skills
                if int(row["outer_fold"]) == fold and row["lineage_id"] == lineage
            )
            for fold in FOLDS
        )
        for lineage in LINEAGES
    }
    secondary = {}
    for metric in ("peak_auprc", "profile_spearman"):
        secondary[metric] = {}
        for model_id, seed in (
            ("scbasset", "ensemble"),
            ("training_lineage_mean", "training_only"),
            ("training_global_mean", "training_only"),
        ):
            rows = [
                row
                for row in secondary_rows
                if row["model_id"] == model_id and row["seed"] == seed
            ]
            secondary[metric][model_id] = _macro_mean(rows, metric)
    positive_seeds = sum(value > 0 for value in seed_scores.values())
    ensemble_skill = _macro_mean(skills, "skill")
    return {
        "strongest_training_only_baseline": strongest,
        "baseline_macro_deviance_per_insertion": {
            model: values["macro_deviance_per_insertion"]
            for model, values in baselines.items()
        },
        "five_seed_ensemble_relative_deviance_reduction": ensemble_skill,
        "paired_two_way_bootstrap": uncertainty,
        "seed_relative_deviance_reduction": seed_scores,
        "positive_gain_seeds": positive_seeds,
        "seed_stability_requirement_met": positive_seeds >= 4,
        "lineage_relative_deviance_reduction": lineage_scores,
        "secondary_macro_metrics": secondary,
        "development_gate": {
            "minimum_relative_deviance_reduction": 0.05,
            "minimum_improved_lineages": 4,
            "maximum_allowed_lineage_worsening": -0.02,
            "minimum_positive_gain_seeds": 4,
            "overall_threshold_passed": ensemble_skill >= 0.05,
            "improved_lineages": sum(value > 0 for value in lineage_scores.values()),
            "no_lineage_worse_than_threshold": min(lineage_scores.values()) >= -0.02,
            "seed_stability_passed": positive_seeds >= 4,
        },
    }


def evaluate(
    *,
    root: Path,
    evaluator_contract: Path,
    completion_root: Path,
    completion_artifacts_sha256: str,
    output: Path,
    workers: int,
    n_resamples: int,
    bootstrap_seed: int,
) -> dict[str, Any]:
    if output.exists() or not 1 <= workers <= 16:
        raise ScBassetFiveSeedEvaluationError("output or worker contract differs")
    root = root.resolve(strict=True)
    contract, inputs, predictions = load_locked_sources(
        root=root,
        evaluator_contract=evaluator_contract,
        completion_root=completion_root,
        completion_artifacts_sha256=completion_artifacts_sha256,
    )

    # No outcome path is resolved until every one of the 25 prediction trees,
    # five input trees, and the full-rectangle completion selection record have verified.
    donor_bigwigs, split_root = _verify_outcome_authorities(root, contract)
    namespace = sha256(
        f"{completion_artifacts_sha256}\0scbasset-five-seed-evaluation".encode()
    ).hexdigest()
    unit_rows: list[dict[str, object]] = []
    secondary_rows: list[dict[str, object]] = []
    donor_sets: list[set[str]] = []
    for fold in FOLDS:
        fold_units, fold_secondary, donors = evaluate_fold(
            fold=fold,
            input_source=inputs[fold],
            prediction_sources=predictions[fold],
            donor_bigwigs=donor_bigwigs,
            split_root=split_root,
            workers=workers,
            namespace=namespace,
        )
        unit_rows.extend(fold_units)
        secondary_rows.extend(fold_secondary)
        donor_sets.append(donors)
    if (
        [len(values) for values in donor_sets] != [9, 7, 4, 8, 11]
        or len(set().union(*donor_sets)) != 39
        or sum(len(values) for values in donor_sets) != 39
    ):
        raise ScBassetFiveSeedEvaluationError("held donors overlap or differ")
    result_summary = summarize(
        unit_rows,
        secondary_rows,
        n_resamples=n_resamples,
        bootstrap_seed=bootstrap_seed,
    )
    output.mkdir(mode=0o750)
    _write_tsv(output / "donor_lineage_block_metrics.tsv", UNIT_FIELDS, unit_rows)
    _write_tsv(
        output / "donor_lineage_secondary_metrics.tsv",
        SECONDARY_FIELDS,
        secondary_rows,
    )
    result = {
        "schema_version": "masld-bench-scbasset-five-seed-evaluation-v1",
        "status": "pass",
        "dataset_id": "gse296875",
        "task_id": "sequence_only_outer_training_atac_profile_development",
        "evaluation_role": "valid",
        "outer_folds": list(FOLDS),
        "seeds": list(SEEDS),
        "fits": 25,
        "prediction_views": 25,
        "n_donors": 39,
        "n_lineages": len(LINEAGES),
        "n_regions_per_fold": 16_000,
        "seed_ensemble": "arithmetic_mean_of_exactly_five_aligned_positive_depth_free_peak_profiles_before_scoring",
        "streaming_memory_policy": "one_seed_prediction_table_in_memory_per_fold",
        "all_prediction_trees_verified_before_outcome_authorities_resolved": True,
        "held_valid_atac_exposed_to_model": False,
        "test_outcomes_read": False,
        "sealed_data_read": False,
        "partial_rectangle_used": False,
        "seeds_are_biological_replicates": False,
        "cells_are_biological_replicates": False,
        "cross_task_family_ranking_performed": False,
        "rna_conditioned_claim_allowed": False,
        "external_evaluation": False,
        "champion_claim_allowed": False,
        "universal_claim_allowed": False,
        "bootstrap_replicates": n_resamples,
        "bootstrap_seed": bootstrap_seed,
        "completion_artifacts_sha256": completion_artifacts_sha256,
        "evaluator_contract_sha256": digest(evaluator_contract.resolve(strict=True)),
        "summary": result_summary,
    }
    (output / "evaluation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--evaluator-contract", type=Path, required=True)
    parser.add_argument("--completion-root", type=Path, required=True)
    parser.add_argument("--completion-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--bootstrap-resamples", type=int, default=10_000)
    parser.add_argument("--bootstrap-seed", type=int, default=17)
    arguments = parser.parse_args()
    print(
        json.dumps(
            evaluate(
                root=arguments.root,
                evaluator_contract=arguments.evaluator_contract,
                completion_root=arguments.completion_root,
                completion_artifacts_sha256=arguments.completion_artifacts_sha256,
                output=arguments.output,
                workers=arguments.workers,
                n_resamples=arguments.bootstrap_resamples,
                bootstrap_seed=arguments.bootstrap_seed,
            ),
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

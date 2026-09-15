#!/usr/bin/env python3
"""Build an outcome-blind crossed-fold fixture for Corgi FiLM-head adaptation."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 in the frozen Corgi runtime.
    import tomli as tomllib


CONTRACT_SHA256 = "615f3b37d094a6d6031c7efd1fd68d8f1f37e0446727aab2bf5abd1fa4ff24b9"
CONTRACT_PREFLIGHT_SHA256 = "378a21b2ccb7b6e9cae865cd11aef61309f8a3f6e5f49bedbe640b5aca4a82a7"
TILE_ROSTER_SHA256 = "2617f3a85a9ced8b0f7cb9efe3819b1e57177629221a2013e39fc7c96e187f5b"
MAPPER_HASHES = {
    0: "0c2f6577c9981039a39409d5fb3593973e089fd9d2d3d665f72947c55d1f7c7a",
    1: "f964222d36d11755e44e77c07a1769f0ffebd3771dbd877909c3c2316911cc5a",
    2: "ceae9d4a53ed6a8bcfaf08aeffe6a396aad12ee1839ebe2da818a75332709507",
    3: "18ddf520a41db2b968a766ced416c85c02bef355becef1c12f5c9b808052b367",
    4: "a2ef2ae490c3cfe429c0c29755d33285085adbe8bf36418e938345c31c3e6053",
}
CONTEXT_ARMS = (
    "actual_released_rank_masked",
    "actual_length_adjusted_tpm_rank_masked",
    "training_lineage_mean_released_rank",
    "nearest_training_released_rank",
    "shuffled_valid_released_rank",
)
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
LINEAGE = "hepatocyte"
SEED_TEXT = "corgi-film-head-20260825-v1"


class CorgiFilmFixtureError(ValueError):
    """Raised when the outcome-blind adaptation fixture does not meet its requirements."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
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


def stable_key(*values: object) -> str:
    return sha256("|".join([SEED_TEXT, *(str(value) for value in values)]).encode()).hexdigest()


def stable_tiles(
    tile_rows: Sequence[dict[str, str]],
    *,
    limit: int,
    phase: str,
    evaluation_fold: int,
    epoch: int,
    donor_id: str,
) -> list[dict[str, str]]:
    if not 0 < limit <= len(tile_rows):
        raise CorgiFilmFixtureError("tile sample limit differs")
    return sorted(
        tile_rows,
        key=lambda row: stable_key(
            phase, evaluation_fold, epoch, donor_id, row["tile_id"]
        ),
    )[:limit]


def _load_mapper(
    mapper: Path, valid_fold: int
) -> tuple[list[dict[str, object]], np.ndarray, np.ndarray, dict[str, Any]]:
    outer_fold = (valid_fold - 1) % 5
    receipt = json.loads((mapper / "receipt.json").read_text(encoding="utf-8"))
    if (
        receipt.get("schema_version") != "masld-bench-corgi-masked-context-mapper-v1"
        or receipt.get("status") != "pass_outcome_free_mapper"
        or receipt.get("outer_fold") != outer_fold
        or receipt.get("valid_fold") != valid_fold
        or receipt.get("held_ATAC_or_other_outcomes_used") is not False
        or receipt.get("test_or_sealed_outcomes_read") is not False
        or receipt.get("structurally_missing_policy")
        != "explicit_mask_then_reference_median_neutral_imputation"
    ):
        raise CorgiFilmFixtureError(f"mapper receipt differs for valid fold {valid_fold}")
    units = read_tsv(mapper / "units_with_roles.tsv")
    with np.load(mapper / "mapped_contexts.npz", allow_pickle=False) as archive:
        required = {
            "released_rank_masked_neutral",
            "length_adjusted_tpm_rank_masked_neutral",
            "training_lineage_mean_released_rank",
            "nearest_training_unit",
            "shuffled_held_unit",
            "length_observed_gene_mask",
        }
        if not required.issubset(archive.files):
            raise CorgiFilmFixtureError("mapper context archive differs")
        mapped = {key: np.asarray(archive[key]) for key in required}
    valid_indices = [
        index
        for index, row in enumerate(units)
        if row["outer_role"] == "valid"
        and row["lineage_id"] == LINEAGE
        and int(row["outer_fold"]) == valid_fold
    ]
    if not valid_indices:
        raise CorgiFilmFixtureError(f"valid hepatocyte contexts absent for fold {valid_fold}")
    lineage_index = LINEAGES.index(LINEAGE)
    records: list[dict[str, object]] = []
    contexts: list[np.ndarray] = []
    masks: list[np.ndarray] = []
    length_mask = np.asarray(mapped["length_observed_gene_mask"], dtype=np.uint8)
    if length_mask.shape != (2_891,) or not np.isin(length_mask, [0, 1]).all():
        raise CorgiFilmFixtureError("length annotation mask differs")
    for unit_index in valid_indices:
        unit = units[unit_index]
        nearest = int(mapped["nearest_training_unit"][unit_index])
        shuffled = int(mapped["shuffled_held_unit"][unit_index])
        arms = (
            (CONTEXT_ARMS[0], mapped["released_rank_masked_neutral"][unit_index], unit_index),
            (
                CONTEXT_ARMS[1],
                mapped["length_adjusted_tpm_rank_masked_neutral"][unit_index],
                unit_index,
            ),
            (
                CONTEXT_ARMS[2],
                mapped["training_lineage_mean_released_rank"][lineage_index],
                -1,
            ),
            (CONTEXT_ARMS[3], mapped["released_rank_masked_neutral"][nearest], nearest),
            (CONTEXT_ARMS[4], mapped["released_rank_masked_neutral"][shuffled], shuffled),
        )
        for arm, vector, source_unit in arms:
            vector = np.asarray(vector, dtype=np.float32)
            if vector.shape != (2_891,) or not np.isfinite(vector).all():
                raise CorgiFilmFixtureError("mapped context vector differs")
            records.append(
                {
                    "context_index": -1,
                    "valid_fold": valid_fold,
                    "mapper_outer_fold": outer_fold,
                    "unit_index": unit_index,
                    "context_source_unit": source_unit,
                    "donor_id": unit["donor_id"],
                    "donor_hash": sha256(
                        f"gse296875-corgi-film-head-v1\0{unit['donor_id']}".encode()
                    ).hexdigest(),
                    "lineage_id": LINEAGE,
                    "context_arm": arm,
                    "rna_role": "valid_oof_mapper",
                    "structural_missing_policy": receipt["structurally_missing_policy"],
                }
            )
            contexts.append(vector)
            masks.append(length_mask)
    return records, np.vstack(contexts), np.vstack(masks), receipt


def _append_schedule(
    rows: list[dict[str, object]],
    *,
    evaluation_fold: int,
    phase: str,
    epoch: int,
    context: Mapping[str, object],
    tile: Mapping[str, str],
) -> None:
    rows.append(
        {
            "schedule_index": len(rows),
            "evaluation_fold": evaluation_fold,
            "phase": phase,
            "epoch": epoch,
            "valid_fold": context["valid_fold"],
            "genomic_fold": int(tile["genomic_fold"]),
            "donor_id": context["donor_id"],
            "donor_hash": context["donor_hash"],
            "context_index": context["context_index"],
            "context_arm": context["context_arm"],
            "tile_id": tile["tile_id"],
            "atac_access": (
                "training_development_atac"
                if phase in {"inner_train", "inner_valid", "refit"}
                else "none_prediction_only"
            ),
        }
    )


def build_fixture(
    *,
    contract_path: Path,
    contract_preflight: Path,
    mapper_roots: Mapping[int, Path],
    tile_roster: Path,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise CorgiFilmFixtureError("refusing to overwrite adaptation fixture")
    if digest(contract_path) != CONTRACT_SHA256:
        raise CorgiFilmFixtureError("FiLM-head contract hash differs")
    if digest(contract_preflight / "ARTIFACTS.json") != CONTRACT_PREFLIGHT_SHA256:
        raise CorgiFilmFixtureError("FiLM-head contract preflight hash differs")
    if digest(tile_roster / "ARTIFACTS.json") != TILE_ROSTER_SHA256:
        raise CorgiFilmFixtureError("tile roster hash differs")
    if set(mapper_roots) != set(range(5)):
        raise CorgiFilmFixtureError("five mapper roots are required")
    for outer_fold, mapper in mapper_roots.items():
        if digest(mapper / "ARTIFACTS.json") != MAPPER_HASHES[outer_fold]:
            raise CorgiFilmFixtureError(f"mapper hash differs for outer fold {outer_fold}")

    contract = tomllib.loads(contract_path.read_text(encoding="utf-8"))
    if (
        contract.get("status") != "frozen_prospective_preflight_execution_not_yet_authorized"
        or contract.get("training_execution_authorized") is not False
        or contract.get("prediction_execution_authorized") is not False
        or contract.get("evaluation_execution_authorized") is not False
        or contract.get("sampling", {}).get("sampling_uses_atac_signal") is not False
    ):
        raise CorgiFilmFixtureError("prospective contract gate differs")

    context_records: list[dict[str, object]] = []
    matrices: list[np.ndarray] = []
    masks: list[np.ndarray] = []
    mapper_receipts: dict[str, dict[str, Any]] = {}
    for valid_fold in range(5):
        outer_fold = (valid_fold - 1) % 5
        records, matrix, mask, receipt = _load_mapper(
            mapper_roots[outer_fold], valid_fold
        )
        for record in records:
            record["context_index"] = len(context_records)
            context_records.append(record)
        matrices.append(matrix)
        masks.append(mask)
        mapper_receipts[str(outer_fold)] = receipt
    context_matrix = np.vstack(matrices)
    length_masks = np.vstack(masks)
    if context_matrix.shape != (39 * len(CONTEXT_ARMS), 2_891):
        raise CorgiFilmFixtureError("39-donor context fixture shape differs")
    if length_masks.shape != context_matrix.shape:
        raise CorgiFilmFixtureError("context missingness-mask shape differs")

    tile_rows = read_tsv(tile_roster / "subset/tiles.tsv")
    by_fold_tiles = {
        fold: [row for row in tile_rows if int(row["genomic_fold"]) == fold]
        for fold in range(5)
    }
    if any(len(rows) != 128 for rows in by_fold_tiles.values()):
        raise CorgiFilmFixtureError("tile roster is not 128 per genomic fold")
    by_fold_context = {
        fold: [row for row in context_records if int(row["valid_fold"]) == fold]
        for fold in range(5)
    }
    by_fold_actual = {
        fold: [row for row in rows if row["context_arm"] == CONTEXT_ARMS[0]]
        for fold, rows in by_fold_context.items()
    }
    folds = contract.get("crossfit_folds", [])
    schedule: list[dict[str, object]] = []
    fold_receipts: list[dict[str, Any]] = []
    for fold in folds:
        evaluation_fold = int(fold["evaluation_fold"])
        inner_train = [int(value) for value in fold["inner_training_folds"]]
        inner_valid = int(fold["inner_validation_fold"])
        fit_folds = [int(value) for value in fold["fit_valid_folds"]]
        start = len(schedule)
        for epoch in range(1, 4):
            for valid_fold in inner_train:
                for context in by_fold_actual[valid_fold]:
                    for tile in stable_tiles(
                        by_fold_tiles[valid_fold],
                        limit=4,
                        phase="inner_train",
                        evaluation_fold=evaluation_fold,
                        epoch=epoch,
                        donor_id=str(context["donor_id"]),
                    ):
                        _append_schedule(
                            schedule,
                            evaluation_fold=evaluation_fold,
                            phase="inner_train",
                            epoch=epoch,
                            context=context,
                            tile=tile,
                        )
            for context in by_fold_actual[inner_valid]:
                for tile in stable_tiles(
                    by_fold_tiles[inner_valid],
                    limit=16,
                    phase="inner_valid",
                    evaluation_fold=evaluation_fold,
                    epoch=0,
                    donor_id=str(context["donor_id"]),
                ):
                    _append_schedule(
                        schedule,
                        evaluation_fold=evaluation_fold,
                        phase="inner_valid",
                        epoch=epoch,
                        context=context,
                        tile=tile,
                    )
            for valid_fold in fit_folds:
                for context in by_fold_actual[valid_fold]:
                    for tile in stable_tiles(
                        by_fold_tiles[valid_fold],
                        limit=4,
                        phase="refit",
                        evaluation_fold=evaluation_fold,
                        epoch=epoch,
                        donor_id=str(context["donor_id"]),
                    ):
                        _append_schedule(
                            schedule,
                            evaluation_fold=evaluation_fold,
                            phase="refit",
                            epoch=epoch,
                            context=context,
                            tile=tile,
                        )
        for context in by_fold_context[evaluation_fold]:
            for tile in by_fold_tiles[evaluation_fold]:
                _append_schedule(
                    schedule,
                    evaluation_fold=evaluation_fold,
                    phase="held_predict",
                    epoch=0,
                    context=context,
                    tile=tile,
                )
        fold_rows = schedule[start:]
        fit_rows = [row for row in fold_rows if row["phase"] != "held_predict"]
        held_rows = [row for row in fold_rows if row["phase"] == "held_predict"]
        if (
            any(
                int(row["valid_fold"]) == evaluation_fold
                or int(row["genomic_fold"]) == evaluation_fold
                for row in fit_rows
            )
            or any(row["context_arm"] != CONTEXT_ARMS[0] for row in fit_rows)
            or any(row["atac_access"] != "training_development_atac" for row in fit_rows)
            or any(row["atac_access"] != "none_prediction_only" for row in held_rows)
            or {row["context_arm"] for row in held_rows} != set(CONTEXT_ARMS)
        ):
            raise CorgiFilmFixtureError(f"fold firewall differs for {evaluation_fold}")
        fit_donors = {str(row["donor_id"]) for row in fit_rows}
        held_donors = {str(row["donor_id"]) for row in held_rows}
        if fit_donors & held_donors:
            raise CorgiFilmFixtureError("held donor entered a fit schedule")
        fold_receipts.append(
            {
                "evaluation_fold": evaluation_fold,
                "inner_training_folds": inner_train,
                "inner_validation_fold": inner_valid,
                "fit_valid_folds": fit_folds,
                "fit_donors": len(fit_donors),
                "held_donors": len(held_donors),
                "fit_held_donor_overlap": 0,
                "fit_held_genomic_fold_overlap": 0,
                "inner_train_rows": sum(row["phase"] == "inner_train" for row in fold_rows),
                "inner_valid_rows": sum(row["phase"] == "inner_valid" for row in fold_rows),
                "refit_rows_max_three_epochs": sum(row["phase"] == "refit" for row in fold_rows),
                "held_prediction_rows": len(held_rows),
            }
        )

    output.mkdir(mode=0o750)
    np.savez_compressed(
        output / "contexts.npz",
        context=context_matrix.astype(np.float32, copy=False),
        length_annotation_observed_mask=length_masks.astype(np.uint8, copy=False),
    )
    write_tsv(
        output / "context_records.tsv",
        tuple(context_records[0]),
        context_records,
    )
    write_tsv(output / "schedule.tsv", tuple(schedule[0]), schedule)
    result: dict[str, Any] = {
        "schema_version": "masld-bench-corgi-film-head-outcome-blind-fixture-v1",
        "status": "pass_outcome_blind_crossed_fold_fixture",
        "model_id": "corgi_regular",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "lineage_id": LINEAGE,
        "adaptation_rung": "film_plus_head",
        "contract_sha256": CONTRACT_SHA256,
        "contract_preflight_artifacts_sha256": CONTRACT_PREFLIGHT_SHA256,
        "tile_roster_artifacts_sha256": TILE_ROSTER_SHA256,
        "mapper_artifacts_sha256_by_outer_fold": {
            str(key): value for key, value in sorted(MAPPER_HASHES.items())
        },
        "folds": fold_receipts,
        "contexts": len(context_records),
        "context_width": int(context_matrix.shape[1]),
        "context_arms": list(CONTEXT_ARMS),
        "donors": len({str(row["donor_id"]) for row in context_records}),
        "schedule_rows": len(schedule),
        "contexts_sha256": digest(output / "contexts.npz"),
        "context_records_sha256": digest(output / "context_records.tsv"),
        "schedule_sha256": digest(output / "schedule.tsv"),
        "held_rna_contexts_materialized": True,
        "atac_bigwig_manifest_or_signal_read": False,
        "development_atac_outcomes_read": False,
        "histology_or_disease_labels_read": False,
        "test_or_sealed_features_or_outcomes_read": False,
        "model_checkpoint_loaded": False,
        "model_fit_performed": False,
        "model_prediction_performed": False,
        "benchmark_metrics_computed": False,
        "training_execution_authorized": False,
        "prediction_execution_authorized": False,
        "evaluation_execution_authorized": False,
        "global_census_modified": False,
        "global_promotion_gate_modified": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--contract",
        type=Path,
        default=root / "config/evaluation/corgi_regular_film_plus_head_smoke.toml",
    )
    parser.add_argument("--contract-preflight", type=Path, required=True)
    parser.add_argument("--tile-roster", type=Path, required=True)
    parser.add_argument("--mapper", action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    mapper_roots: dict[int, Path] = {}
    for value in arguments.mapper:
        fold_text, separator, path_text = value.partition("=")
        if not separator or not fold_text.isdigit():
            raise CorgiFilmFixtureError("mapper arguments must be OUTER_FOLD=PATH")
        mapper_roots[int(fold_text)] = Path(path_text).resolve(strict=True)
    result = build_fixture(
        contract_path=arguments.contract.resolve(strict=True),
        contract_preflight=arguments.contract_preflight.resolve(strict=True),
        mapper_roots=mapper_roots,
        tile_roster=arguments.tile_roster.resolve(strict=True),
        output=arguments.output.resolve(strict=False),
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

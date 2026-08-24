#!/usr/bin/env python3
"""Freeze one donor-fold Corgi mapper and outcome-free context ablations."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata


LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
SEED = 20260824


class CorgiMapperError(RuntimeError):
    """Raised when the split, mask, or context representation differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def masked_quantile_map(values: np.ndarray, observed: np.ndarray, reference: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    observed = np.asarray(observed, dtype=np.bool_)
    reference = np.asarray(reference, dtype=np.float64)
    if values.shape != observed.shape or values.ndim != 1 or reference.shape != values.shape:
        raise CorgiMapperError("masked quantile input shape differs")
    if observed.sum() < 2 or not np.isfinite(values[observed]).all() or np.any(values[observed] < 0):
        raise CorgiMapperError("masked quantile observed values differ")
    sorted_reference = np.sort(reference)
    percentiles = (rankdata(values[observed], method="average") - 1.0) / (observed.sum() - 1.0)
    mapped = np.full(values.shape, np.median(sorted_reference), dtype=np.float32)
    mapped[observed] = np.interp(
        percentiles,
        np.linspace(0.0, 1.0, len(sorted_reference)),
        sorted_reference,
    ).astype(np.float32)
    if not np.isfinite(mapped).all():
        raise CorgiMapperError("masked quantile output is non-finite")
    return mapped


def derangement(indices: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    indices = np.asarray(indices, dtype=np.int64)
    if len(indices) < 2:
        raise CorgiMapperError("context shuffle requires at least two units")
    for _ in range(1_000):
        candidate = rng.permutation(indices)
        if np.all(candidate != indices):
            return candidate
    raise CorgiMapperError("could not construct deterministic derangement")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contexts", type=Path, required=True)
    parser.add_argument("--outer-fold", type=int, choices=range(5), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with np.load(args.contexts / "raw_context_counts.npz", allow_pickle=False) as source:
        counts = source["counts"].astype(np.float64)
        observed = source["observed_gene_mask"].astype(np.bool_)
        reference = source["tf_reference"].astype(np.float64)
    units = read_tsv(args.contexts / "units.tsv")
    genes = read_tsv(args.contexts / "corgi_gene_crosswalk.tsv")
    if counts.shape != (195, 2_891) or observed.shape != (2_891,) or int(observed.sum()) != 2_857:
        raise CorgiMapperError("raw context shape or mask differs")
    if len(units) != 195 or len(genes) != 2_891:
        raise CorgiMapperError("context manifest dimensions differ")
    if any(int(row["unit_index"]) != index for index, row in enumerate(units)):
        raise CorgiMapperError("unit order differs")
    length_values = np.asarray(
        [float(row["gencode_v49_exon_union_length"]) if row["gencode_v49_exon_union_length"] else np.nan for row in genes],
        dtype=np.float64,
    )
    length_observed = observed & np.isfinite(length_values) & (length_values > 0)
    if length_observed.sum() < 2_800:
        raise CorgiMapperError("too few observed genes have GENCODE v49 exon lengths")
    count_contexts = np.vstack(
        [masked_quantile_map(row, observed, reference) for row in counts]
    ).astype(np.float32)
    rates = np.zeros_like(counts, dtype=np.float64)
    rates[:, length_observed] = counts[:, length_observed] / (length_values[length_observed] / 1_000.0)
    tpm_contexts = np.vstack(
        [masked_quantile_map(row, length_observed, reference) for row in rates]
    ).astype(np.float32)
    folds = np.asarray([int(row["outer_fold"]) for row in units], dtype=np.int8)
    valid_fold = (args.outer_fold + 1) % 5
    roles = np.full(len(units), "train", dtype="U5")
    roles[folds == args.outer_fold] = "test"
    roles[folds == valid_fold] = "valid"
    if set(roles.tolist()) != {"train", "valid", "test"}:
        raise CorgiMapperError("outer role assignment differs")
    lineage_values = np.asarray([row["lineage_id"] for row in units])
    means_count = np.zeros((len(LINEAGES), 2_891), dtype=np.float32)
    means_tpm = np.zeros_like(means_count)
    nearest = np.full(len(units), -1, dtype=np.int32)
    shuffled = np.full(len(units), -1, dtype=np.int32)
    observed_counts = counts[:, observed]
    library = observed_counts.sum(axis=1, keepdims=True)
    if np.any(library <= 0):
        raise CorgiMapperError("unit RNA library is empty")
    log_cpm = np.log1p(observed_counts * (1_000_000.0 / library))
    rng = np.random.default_rng(SEED + args.outer_fold)
    for lineage_index, lineage in enumerate(LINEAGES):
        train = np.flatnonzero((roles == "train") & (lineage_values == lineage))
        means_count[lineage_index] = count_contexts[train].mean(axis=0, dtype=np.float64)
        means_tpm[lineage_index] = tpm_contexts[train].mean(axis=0, dtype=np.float64)
        for role in ("valid", "test"):
            held = np.flatnonzero((roles == role) & (lineage_values == lineage))
            shuffled[held] = derangement(held, rng)
            for unit in held:
                distances = np.square(log_cpm[train] - log_cpm[unit]).mean(axis=1)
                nearest[unit] = int(train[int(np.argmin(distances))])
    held_mask = roles != "train"
    if np.any(nearest[held_mask] < 0) or np.any(shuffled[held_mask] < 0):
        raise CorgiMapperError("held context ablation mapping is incomplete")
    args.output.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(
        args.output / "mapped_contexts.npz",
        released_rank_masked_neutral=count_contexts,
        length_adjusted_tpm_rank_masked_neutral=tpm_contexts,
        training_lineage_mean_released_rank=means_count,
        training_lineage_mean_length_adjusted_tpm_rank=means_tpm,
        nearest_training_unit=nearest,
        shuffled_held_unit=shuffled,
        outer_role=roles,
        length_observed_gene_mask=length_observed,
    )
    role_rows = [
        {
            **row,
            "outer_role": roles[index],
            "nearest_training_unit": int(nearest[index]),
            "shuffled_held_unit": int(shuffled[index]),
        }
        for index, row in enumerate(units)
    ]
    with (args.output / "units_with_roles.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(role_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(role_rows)
    receipt = {
        "schema_version": "masld-bench-corgi-masked-context-mapper-v1",
        "status": "pass_outcome_free_mapper",
        "outer_fold": args.outer_fold,
        "valid_fold": valid_fold,
        "seed": SEED,
        "units": 195,
        "training_units": int((roles == "train").sum()),
        "valid_units": int((roles == "valid").sum()),
        "test_units": int((roles == "test").sum()),
        "observed_count_genes": int(observed.sum()),
        "length_observed_genes": int(length_observed.sum()),
        "structurally_missing_policy": "explicit_mask_then_reference_median_neutral_imputation",
        "native_missing_as_zero_mapper_run": False,
        "tie_policy": "average_rank",
        "training_only_objects": ["lineage_mean_contexts", "nearest_context_reference_bank"],
        "held_RNA_used_for_context": True,
        "held_ATAC_or_other_outcomes_used": False,
        "test_or_sealed_outcomes_read": False,
        "input_artifacts_sha256": digest(args.contexts / "ARTIFACTS.json"),
        "eligible_next_action": "valid_role_outcome_free_Corgi_prediction",
    }
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

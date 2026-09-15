#!/usr/bin/env python3
"""Validate the outcome-blind SGD elastic-net rescue on one outer training fold."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from scripts.fit_predict_cell_baselines_study_50000 import (
    N_PCA,
    N_TOP_HVG,
    ROSTER,
    TARGET_SUM,
    _log_normalize,
    _select_hvgs,
    canonical_sha256,
    donor_class_weights,
    fit_sgd_elastic_net,
    prepare_inner_projections,
    select_sgd_elastic_net,
    standardize_projection,
)


DIAGNOSTIC_504_SHA256 = "bbccba9f03aaa9736a9ac9146e25d00eb7c9fa9b6b2a06c11ce62e4c6b26306e"
DIAGNOSTIC_505_SHA256 = "8f00036292a60a1712e4da0f5dd8f9e9c4e8a47541e4af8a15778d145d2124c0"


class SGDElasticRescueError(RuntimeError):
    """Raised when the frozen diagnosis or training-only rescue differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def audit_frozen_saga_diagnostics(diagnostic_504: Path, diagnostic_505: Path) -> dict[str, object]:
    if (
        sha256_file(diagnostic_504 / "ARTIFACTS.json") != DIAGNOSTIC_504_SHA256
        or sha256_file(diagnostic_505 / "ARTIFACTS.json") != DIAGNOSTIC_505_SHA256
    ):
        raise SGDElasticRescueError("frozen SAGA diagnostic identity differs")
    initial = json.loads(
        (diagnostic_504 / "diagnostic/diagnostic.json").read_text(encoding="utf-8")
    )
    if (
        initial.get("outer_fold") != 0
        or initial.get("metrics_calculated") is not False
        or initial.get("test_partition_predicted") is not False
        or initial.get("sealed_outcomes_read") is not False
    ):
        raise SGDElasticRescueError("initial SAGA diagnostic firewall differs")
    failed_50k: list[int] = []
    failed_loose: list[int] = []
    max_weights: list[float] = []
    for fold in range(5):
        value = json.loads(
            (diagnostic_505 / f"fold-{fold}/diagnostic/diagnostic.json").read_text(
                encoding="utf-8"
            )
        )
        fits = {item["config_id"]: item for item in value.get("fits", [])}
        if (
            value.get("outer_fold") != fold
            or value.get("metrics_calculated") is not False
            or value.get("test_partition_predicted") is not False
            or value.get("sealed_outcomes_read") is not False
            or set(fits)
            != {"float64_tol1e-3_c1_max50000", "float64_tol1e-2_c1"}
        ):
            raise SGDElasticRescueError("all-fold SAGA diagnostic firewall differs")
        if not fits["float64_tol1e-3_c1_max50000"]["converged_before_ceiling"]:
            failed_50k.append(fold)
        if not fits["float64_tol1e-2_c1"]["converged_before_ceiling"]:
            failed_loose.append(fold)
        max_weights.append(float(value["training_weight_quantiles"]["max"]))
    if failed_50k != [1, 2, 3, 4] or failed_loose != [2, 3]:
        raise SGDElasticRescueError("frozen SAGA failure pattern differs")
    return {
        "saga_50000_failed_outer_folds": failed_50k,
        "saga_tol1e2_failed_outer_folds": failed_loose,
        "maximum_training_cell_weight": max(max_weights),
        "diagnostic_504_artifacts_sha256": DIAGNOSTIC_504_SHA256,
        "diagnostic_505_artifacts_sha256": DIAGNOSTIC_505_SHA256,
    }


def run(arguments: argparse.Namespace) -> dict[str, object]:
    import anndata
    from sklearn.decomposition import PCA

    if arguments.output.exists():
        raise SGDElasticRescueError("refusing to overwrite rescue diagnostic")
    history = audit_frozen_saga_diagnostics(arguments.diagnostic_504, arguments.diagnostic_505)
    if (
        sha256_file(arguments.source / "ARTIFACTS.json")
        != arguments.source_artifacts_sha256
        or sha256_file(arguments.split / "ARTIFACTS.json")
        != arguments.split_artifacts_sha256
    ):
        raise SGDElasticRescueError("source or split identity differs")
    split_rows = read_tsv(arguments.split / "row_outer_folds.tsv")
    inner_rows = read_tsv(arguments.split / "inner_donor_folds.tsv")
    adata = anndata.read_h5ad(
        arguments.source / "resource_atlas_frozen_screen_50000.h5ad"
    )
    row_ids = np.asarray(list(map(str, adata.obs["row_id"])), dtype=str)
    donors = np.asarray([row["donor_id"] for row in split_rows], dtype=str)
    outer = np.asarray([int(row["outer_fold"]) for row in split_rows], dtype=np.int8)
    labels = np.asarray(list(map(str, adata.obs["broad_label"])), dtype=str)
    if (
        len(split_rows) != 50_000
        or adata.n_obs != 50_000
        or row_ids.tolist() != [row["row_id"] for row in split_rows]
        or set(outer) != set(range(5))
        or any(len(set(outer[donors == donor])) != 1 for donor in set(donors))
    ):
        raise SGDElasticRescueError("source or donor-safe split differs")
    train = np.flatnonzero(outer != arguments.outer_fold)
    train_labels = labels[train]
    label_to_index = {label: index for index, label in enumerate(ROSTER)}
    if set(train_labels) != set(ROSTER):
        raise SGDElasticRescueError("outer training class roster differs")
    train_y = np.asarray([label_to_index[label] for label in train_labels], dtype=np.int64)
    train_donors = donors[train]
    inner_assignment = {
        row["donor_id"]: int(row["inner_fold"])
        for row in inner_rows
        if int(row["held_outer_fold"]) == arguments.outer_fold
    }
    if (
        set(inner_assignment) != set(train_donors)
        or set(inner_assignment.values()) != set(range(5))
    ):
        raise SGDElasticRescueError("inner donor assignment differs")
    raw_train = adata.X[train]
    gene_ids = list(map(str, adata.var["ensembl_id"]))
    prepared = prepare_inner_projections(
        raw_train=raw_train,
        gene_ids=gene_ids,
        train_y=train_y,
        train_donors=train_donors,
        inner_assignment=inner_assignment,
        fold=arguments.outer_fold,
        seed=arguments.seed,
    )
    selected, candidates = select_sgd_elastic_net(
        prepared, train_y, fold=arguments.outer_fold, seed=arguments.seed
    )
    normalized = _log_normalize(raw_train, TARGET_SUM)
    selected_features, _ = _select_hvgs(
        normalized, gene_ids, n_top=N_TOP_HVG, n_bins=20
    )
    dense = normalized[:, selected_features].toarray().astype(np.float32, copy=False)
    pca = PCA(
        n_components=N_PCA,
        svd_solver="randomized",
        whiten=False,
        random_state=arguments.seed + arguments.outer_fold,
    )
    training_x = pca.fit_transform(dense).astype(np.float32, copy=False)
    training_x, _, _, _ = standardize_projection(training_x, training_x[:1])
    weights = np.asarray(
        donor_class_weights(
            list(train_donors), [ROSTER[index] for index in train_y]
        ),
        dtype=np.float64,
    )
    final, final_state = fit_sgd_elastic_net(
        training_x,
        train_y,
        weights,
        alpha=selected["alpha"],
        l1_ratio=selected["l1_ratio"],
        seed=arguments.seed + arguments.outer_fold,
    )
    if not final_state["converged_before_ceiling"]:
        raise SGDElasticRescueError("selected outer-training SGD fit did not converge")
    training_probabilities = final.predict_proba(training_x)
    if (
        training_probabilities.shape != (len(train), len(ROSTER))
        or not np.all(np.isfinite(training_probabilities))
    ):
        raise SGDElasticRescueError("training-only probabilities differ")
    result = {
        "schema_version": "masld-bench-sgd-elastic-net-rescue-diagnostic-v1",
        "status": "passed_training_only_rescue",
        "outer_fold": arguments.outer_fold,
        "seed": arguments.seed,
        "training_rows": len(train),
        "training_donors": len(set(train_donors)),
        "training_classes": list(ROSTER),
        "selected": selected,
        "candidate_audit": candidates,
        "final_fit": final_state,
        "selected_feature_ids_sha256": canonical_sha256(
            [gene_ids[index] for index in selected_features]
        ),
        "training_probability_sha256": sha256(
            np.ascontiguousarray(training_probabilities).tobytes()
        ).hexdigest(),
        "historical_saga_diagnosis": history,
        "selection_partition": "donor_grouped_inner_training_only",
        "selection_metric": "donor_class_balanced_multiclass_log_loss",
        "outer_test_features_transformed": False,
        "outer_test_partition_predicted": False,
        "outer_test_metrics_calculated": False,
        "sealed_outcomes_read": False,
    }
    arguments.output.mkdir(parents=True, exist_ok=False)
    (arguments.output / "diagnostic.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--source", required=True, type=Path)
    value.add_argument("--source-artifacts-sha256", required=True)
    value.add_argument("--split", required=True, type=Path)
    value.add_argument("--split-artifacts-sha256", required=True)
    value.add_argument("--diagnostic-504", required=True, type=Path)
    value.add_argument("--diagnostic-505", required=True, type=Path)
    value.add_argument("--outer-fold", required=True, type=int, choices=range(5))
    value.add_argument("--seed", type=int, default=20260824)
    value.add_argument("--output", required=True, type=Path)
    return value


def main() -> int:
    result = run(parser().parse_args())
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Training-only convergence audit for frozen GSE274114 trainer views."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from typing import Any
import warnings

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import canonical_sha256, sha256_file
from scripts.gse274114_baseline_firewall import (
    FIXED_STOCHASTIC_SEEDS,
    GSE274114BaselineFirewallError,
    STOCHASTIC_MODELS,
    _fit_transform,
    read_tsv,
    write_tsv,
)


class GSE274114ConvergenceAuditError(RuntimeError):
    """Raised when a diagnostic would cross its training-only boundary."""


def _candidate_grid(model_kind: str) -> list[tuple[float, float | None]]:
    if model_kind == "hvg_pca_elastic_net":
        return [
            (c, ratio)
            for c in (0.1, 1.0, 10.0)
            for ratio in (0.0, 0.5, 1.0)
        ]
    if model_kind == "hvg_pca_linear_svm":
        return [(c, None) for c in (0.1, 1.0, 10.0)]
    raise GSE274114ConvergenceAuditError("diagnostic model is not stochastic")


def _make_model(model_kind: str, *, c: float, ratio: float | None, seed: int) -> Any:
    if model_kind == "hvg_pca_elastic_net":
        from sklearn.linear_model import LogisticRegression

        return LogisticRegression(
            penalty="elasticnet",
            solver="saga",
            C=c,
            l1_ratio=ratio,
            class_weight="balanced",
            max_iter=20_000,
            random_state=seed,
            tol=1e-6,
        )
    if model_kind == "hvg_pca_linear_svm":
        from sklearn.svm import SVC

        return SVC(
            C=c,
            kernel="linear",
            probability=True,
            class_weight="balanced",
            random_state=seed,
        )
    raise GSE274114ConvergenceAuditError("diagnostic model is not stochastic")


def _fit_trace(model: Any, values: Any, labels: Any) -> dict[str, Any]:
    from sklearn.exceptions import ConvergenceWarning

    caught: list[warnings.WarningMessage]
    error: str | None = None
    with warnings.catch_warnings(record=True) as caught:
        warnings.filterwarnings("ignore", category=FutureWarning)
        warnings.simplefilter("always", ConvergenceWarning)
        try:
            model.fit(values, labels)
        except (ArithmeticError, FloatingPointError, ValueError) as exception:
            error = f"{type(exception).__name__}:{exception}"
    convergence = [
        str(record.message)
        for record in caught
        if issubclass(record.category, ConvergenceWarning)
    ]
    n_iter = getattr(model, "n_iter_", None)
    if n_iter is not None:
        n_iter = [int(value) for value in __import__("numpy").ravel(n_iter)]
    return {
        "fit_error": error,
        "convergence_warning_count": len(convergence),
        "convergence_warning_text": convergence,
        "n_iter": n_iter,
        "converged_without_warning": error is None and not convergence,
    }


def _audit_one(arguments: dict[str, Any]) -> dict[str, Any]:
    import numpy as np
    from sklearn.model_selection import StratifiedKFold

    view_root = Path(arguments["view_root"])
    expected_sha256 = str(arguments["view_sha256"])
    model_kind = str(arguments["model_kind"])
    seed = int(arguments["seed"])
    outer_fold = int(arguments["outer_fold"])
    if model_kind not in STOCHASTIC_MODELS or seed not in FIXED_STOCHASTIC_SEEDS:
        raise GSE274114ConvergenceAuditError("diagnostic run identity differs")
    if sha256_file(view_root / "ARTIFACTS.json") != expected_sha256:
        raise GSE274114ConvergenceAuditError("trainer view identity differs")
    manifest = verify_frozen_tree(view_root)
    metadata = manifest["metadata"]
    if (
        metadata.get("artifact_class") != "gse274114_trainer_fold_view"
        or metadata.get("outer_fold") != outer_fold
        or metadata.get("training_labels_included") is not True
        or metadata.get("query_labels_included") is not False
        or metadata.get("participant_covariates_included") is not False
        or metadata.get("model_fit_or_scoring_run") is not False
    ):
        raise GSE274114ConvergenceAuditError("trainer view firewall differs")
    fields, rows = read_tsv(view_root / "training_participants.tsv")
    if fields != ("row_id", "target"):
        raise GSE274114ConvergenceAuditError("training participant schema differs")
    train_counts = np.load(view_root / "training_numreads.npy", allow_pickle=False)
    y = np.asarray([int(row["target"]) for row in rows], dtype=np.int64)
    if train_counts.shape[0] != len(y) or set(y.tolist()) != {0, 1}:
        raise GSE274114ConvergenceAuditError("training matrix or labels differ")

    outer_x, _, preprocessing = _fit_transform(train_counts, train_counts[:1])
    splitter = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    candidates: list[dict[str, Any]] = []
    for c, ratio in _candidate_grid(model_kind):
        inner_traces: list[dict[str, Any]] = []
        losses: list[float] = []
        for inner_fold, (inner_train, inner_valid) in enumerate(
            splitter.split(train_counts, y)
        ):
            inner_x, valid_x, inner_preprocessing = _fit_transform(
                train_counts[inner_train], train_counts[inner_valid]
            )
            model = _make_model(model_kind, c=c, ratio=ratio, seed=seed)
            trace = _fit_trace(model, inner_x, y[inner_train])
            loss: float | None = None
            if trace["converged_without_warning"]:
                probability = np.clip(
                    model.predict_proba(valid_x)[:, 1], 1e-8, 1 - 1e-8
                )
                truth = y[inner_valid]
                loss = float(
                    np.mean(
                        -(
                            truth * np.log(probability)
                            + (1 - truth) * np.log(1 - probability)
                        )
                    )
                )
                losses.append(loss)
            inner_traces.append(
                {
                    "inner_fold": inner_fold,
                    "training_participants": len(inner_train),
                    "validation_participants": len(inner_valid),
                    "preprocessing_sha256": inner_preprocessing[
                        "training_preprocessing_sha256"
                    ],
                    "log_loss": loss,
                    **trace,
                }
            )
        inner_valid = len(losses) == 3
        mean_loss = float(np.mean(losses)) if inner_valid else None
        outer_model = _make_model(model_kind, c=c, ratio=ratio, seed=seed)
        outer_trace = _fit_trace(outer_model, outer_x, y)
        candidates.append(
            {
                "C": c,
                "l1_ratio": ratio,
                "inner_candidate_valid": inner_valid,
                "mean_inner_log_loss": mean_loss,
                "inner_traces": inner_traces,
                "outer_trace": outer_trace,
            }
        )
    ranked = sorted(
        (row for row in candidates if row["inner_candidate_valid"]),
        key=lambda row: (
            row["mean_inner_log_loss"],
            row["C"],
            -1 if row["l1_ratio"] is None else row["l1_ratio"],
        ),
    )
    for rank, row in enumerate(ranked, start=1):
        row["inner_valid_rank"] = rank
    selected = ranked[0] if ranked else None
    return {
        "outer_fold": outer_fold,
        "model_kind": model_kind,
        "seed": seed,
        "training_participants": len(y),
        "training_class_counts": {
            "0": int(np.sum(y == 0)),
            "1": int(np.sum(y == 1)),
        },
        "outer_preprocessing_sha256": preprocessing[
            "training_preprocessing_sha256"
        ],
        "candidates": candidates,
        "selected_inner_candidate": (
            None
            if selected is None
            else {
                "C": selected["C"],
                "l1_ratio": selected["l1_ratio"],
                "mean_inner_log_loss": selected["mean_inner_log_loss"],
                "outer_converged_without_warning": selected["outer_trace"][
                    "converged_without_warning"
                ],
            }
        ),
    }


def run_audit(
    *,
    task_id: str,
    view_roots: list[Path],
    view_sha256s: list[str],
    workers: int,
    output: Path,
) -> dict[str, Any]:
    if (
        output.exists()
        or len(view_roots) != 5
        or len(view_sha256s) != 5
        or not 1 <= workers <= 16
    ):
        raise GSE274114ConvergenceAuditError("audit inputs differ")
    runs = [
        {
            "task_id": task_id,
            "outer_fold": fold,
            "view_root": view_roots[fold].as_posix(),
            "view_sha256": view_sha256s[fold],
            "model_kind": model_kind,
            "seed": seed,
        }
        for fold in range(5)
        for model_kind in sorted(STOCHASTIC_MODELS)
        for seed in FIXED_STOCHASTIC_SEEDS
    ]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(_audit_one, runs, chunksize=1))
    results.sort(key=lambda row: (row["model_kind"], row["outer_fold"], row["seed"]))
    if len(results) != 50:
        raise GSE274114ConvergenceAuditError("diagnostic run census differs")

    output.mkdir(parents=True)
    write_json_exclusive(
        output / "convergence_trace.json",
        {
            "schema_version": "masld-bench-gse274114-training-convergence-trace-v1",
            "task_id": task_id,
            "runs": results,
        },
    )
    summary_rows: list[dict[str, Any]] = []
    for result in results:
        selected = result["selected_inner_candidate"]
        summary_rows.append(
            {
                "outer_fold": result["outer_fold"],
                "model_kind": result["model_kind"],
                "seed": result["seed"],
                "inner_valid_candidates": sum(
                    row["inner_candidate_valid"] for row in result["candidates"]
                ),
                "selected_c": "" if selected is None else selected["C"],
                "selected_l1_ratio": (
                    "" if selected is None or selected["l1_ratio"] is None else selected["l1_ratio"]
                ),
                "selected_mean_inner_log_loss": (
                    "" if selected is None else selected["mean_inner_log_loss"]
                ),
                "selected_outer_converged_without_warning": (
                    False if selected is None else selected["outer_converged_without_warning"]
                ),
                "outer_convergent_inner_valid_candidates": sum(
                    row["inner_candidate_valid"]
                    and row["outer_trace"]["converged_without_warning"]
                    for row in result["candidates"]
                ),
            }
        )
    write_tsv(
        output / "selected_candidate_summary.tsv",
        (
            "outer_fold",
            "model_kind",
            "seed",
            "inner_valid_candidates",
            "selected_c",
            "selected_l1_ratio",
            "selected_mean_inner_log_loss",
            "selected_outer_converged_without_warning",
            "outer_convergent_inner_valid_candidates",
        ),
        summary_rows,
    )
    failures = [
        row for row in summary_rows if not row["selected_outer_converged_without_warning"]
    ]
    receipt = {
        "schema_version": "masld-bench-gse274114-training-convergence-audit-v1",
        "status": "passed_training_only_convergence_audit",
        "task_id": task_id,
        "diagnostic_runs": len(results),
        "selected_candidate_outer_nonconvergence_count": len(failures),
        "selected_candidate_outer_nonconvergence": failures,
        "training_views_opened": True,
        "training_labels_opened": True,
        "query_molecular_values_opened": False,
        "query_labels_opened": False,
        "participant_covariates_opened": False,
        "predictions_generated": False,
        "scoring_or_evaluation_run": False,
        "diagnostic_identity_sha256": canonical_sha256(
            {
                "task_id": task_id,
                "view_sha256s": view_sha256s,
                "model_kinds": sorted(STOCHASTIC_MODELS),
                "seeds": FIXED_STOCHASTIC_SEEDS,
            }
        ),
    }
    write_json_exclusive(output / "receipt.json", receipt)
    freeze_tree(output, {"artifact_class": "gse274114_training_only_convergence_audit", **receipt})
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--view-root", action="append", type=Path, required=True)
    parser.add_argument("--view-sha256", action="append", required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    receipt = run_audit(
        task_id=arguments.task_id,
        view_roots=arguments.view_root,
        view_sha256s=arguments.view_sha256,
        workers=arguments.workers,
        output=arguments.output,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

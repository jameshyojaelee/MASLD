#!/usr/bin/env python3
"""Refit the frozen GSE281364 HyenaDNA head under shortened input windows.

Every fitting choice is reused from the published campaign
(`scripts/fit_gse281364_dna_language_seeded_heads.py`) by importing its helpers: the
239-block outer folds, the (outer+1) mod 5 inner validation fold, the five seeds and
their bootstrap salts, the whitened 256-wide projection refit per fold, the alpha and
top-k grids, and the 10,000-draw paired block bootstrap.  Only the raw embedding matrix
changes between arms.

The allele-identity ridge reference is copied verbatim from the frozen Enformer/Sei
out-of-fold predictions, so all arms are compared against one identical control and the
bootstrap draws are shared.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.stats import rankdata, spearmanr

from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes
from scripts.fit_gse281364_dna_language_seeded_heads import fit_seeded_dual_projection_fold
from scripts.fit_gse281364_open_sequence_heads import fit_inner_projection
from scripts.gse281364_dna_lm_native_contract import apply_projection, head_features


CONTEXTS = ("HepG2_control", "HepG2_PAOA")
HEADS = ("delta_ridge", "full_ridge")
SEEDS = (1103, 2909, 4721, 6673, 8111)
ALPHAS = (0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)
DELTA_TOP_K = (32, 128, 256)
FULL_TOP_K = (32, 128, 512, 1024)
BOOTSTRAP_RESAMPLES = 10_000
BOOTSTRAP_SEED = 20_260_825
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
PUBLISHED_SPEARMAN = 0.21622332457761556
PUBLISHED_GAIN = 0.082007807505673547
REPRODUCTION_TOLERANCE = 1.0e-9
EXTRACTION_TOLERANCE = 0.010


class AblationError(RuntimeError):
    """Raised when a reproduction, alignment, or census requirement is not met."""


def _read_tsv(path: Path, fields: Sequence[str], *, compressed: bool = False) -> list[dict[str, str]]:
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise AblationError(f"TSV schema differs: {path}")
        return [dict(row) for row in reader]


def _write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _macro(values: Sequence[float | None]) -> float | None:
    if any(value is None or not math.isfinite(value) for value in values):
        return None
    clipped = np.clip(np.asarray(values), -0.999999, 0.999999)
    return float(np.tanh(np.mean(np.arctanh(clipped))))


def _fast_spearman(observed: np.ndarray, predicted: np.ndarray) -> float | None:
    if np.all(observed == observed[0]) or np.all(predicted == predicted[0]):
        return None
    left = rankdata(observed)
    right = rankdata(predicted)
    left -= left.mean()
    right -= right.mean()
    denominator = float(np.sqrt(np.sum(left**2) * np.sum(right**2)))
    return None if denominator == 0 else float(np.sum(left * right) / denominator)


def _load_arm(path: Path, fixture_to_element: Mapping[str, str], elements: Sequence[str]):
    """Return embeddings in npz order plus the permutation into fit element order."""
    with np.load(path, allow_pickle=False) as data:
        fixture_ids = data["fixture_ids"].astype(str)
        allele_order = tuple(data["allele_order"].astype(str).tolist())
        embeddings = data["embeddings"]
    if (
        allele_order != ("REF", "ALT", "REF_RC", "ALT_RC")
        or embeddings.shape != (1033, 4, 256)
        or not np.isfinite(embeddings).all()
    ):
        raise AblationError(f"arm embedding geometry differs: {path}")
    position = {fixture_to_element[fixture]: index for index, fixture in enumerate(fixture_ids)}
    if set(position) != set(elements):
        raise AblationError(f"arm element coverage differs: {path}")
    return embeddings, np.asarray([position[element] for element in elements])


def _arm_predictions(
    embeddings: np.ndarray,
    reorder: np.ndarray,
    folds: np.ndarray,
    blocks: np.ndarray,
    targets: Mapping[str, np.ndarray],
) -> dict[tuple[str, int, str], np.ndarray]:
    """Fit every head/seed/context for one arm under the frozen protocol."""
    # The published outer projection was fitted in fixture order and reordered
    # afterwards; the inner projection was fitted in element order.  Both are
    # reproduced here so the frozen arm matches bit for bit.
    element_raw = embeddings[reorder]
    npz_folds = np.empty_like(folds)
    npz_folds[reorder] = folds
    inner_features: dict[int, np.ndarray] = {}
    outer_features: dict[int, np.ndarray] = {}
    for held_fold in range(5):
        inner_fold = (held_fold + 1) % 5
        inner_mask = (folds != held_fold) & (folds != inner_fold)
        inner_features[held_fold] = head_features(
            apply_projection(element_raw, fit_inner_projection(element_raw, inner_mask))
        ).astype(np.float64)
        outer = head_features(
            apply_projection(embeddings, fit_inner_projection(embeddings, npz_folds != held_fold))
        )
        outer_features[held_fold] = outer[reorder].astype(np.float64)
    predictions: dict[tuple[str, int, str], np.ndarray] = {}
    for context in CONTEXTS:
        y = targets[context]
        for head in HEADS:
            feature_slice = slice(512, 768) if head == "delta_ridge" else slice(0, 1024)
            top_k = DELTA_TOP_K if head == "delta_ridge" else FULL_TOP_K
            for seed in SEEDS:
                oof = np.empty(1033)
                for held_fold in range(5):
                    prediction, _, _ = fit_seeded_dual_projection_fold(
                        inner_features=inner_features[held_fold][:, feature_slice],
                        outer_features=outer_features[held_fold][:, feature_slice],
                        outcomes=y,
                        folds=folds,
                        block_ids=blocks,
                        held_fold=held_fold,
                        seed=seed,
                        top_k_grid=top_k,
                        alphas=ALPHAS,
                    )
                    oof[folds == held_fold] = prediction
                if not np.isfinite(oof).all():
                    raise AblationError("OOF prediction coverage differs")
                predictions[(head, seed, context)] = oof
    return predictions


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.root.resolve(strict=True)
    output = arguments.output
    if output.exists():
        raise AblationError("ablation output exists")

    row_rows = _read_tsv(root / arguments.row_universe, ROW_FIELDS)
    if len(row_rows) != 10330:
        raise AblationError("row denominator differs")
    metadata: dict[str, tuple[str, str]] = {}
    for row in row_rows:
        value = (row["long_range_block_id"], row["outer_fold"])
        if metadata.setdefault(row["element_id"], value) != value:
            raise AblationError("element metadata differs")
    # Preserve the frozen row-universe insertion order used by the published fit.
    fit_elements = list(metadata)
    folds = np.asarray(
        [int(metadata[element][1].removeprefix("fold-")) for element in fit_elements],
        dtype=np.int64,
    )
    blocks = np.asarray([metadata[element][0] for element in fit_elements])
    if len(fit_elements) != 1033 or len(set(blocks.tolist())) != 239:
        raise AblationError("element or block census differs")

    outcome_data = load_outcomes(root / arguments.outcomes, set(fit_elements))
    targets = {
        context: np.asarray([outcome_data[(element, context)]["mean"] for element in fit_elements])
        for context in CONTEXTS
    }

    with (root / arguments.fixture_manifest).open(encoding="utf-8", newline="") as handle:
        fixture_to_element = {
            row["fixture_id"]: row["element_id"]
            for row in csv.DictReader(handle, delimiter="\t")
        }
    if len(fixture_to_element) != 1033 or set(fixture_to_element.values()) != set(fit_elements):
        raise AblationError("fixture to element mapping differs")

    control: dict[tuple[int, str], float] = {}
    with gzip.open(root / arguments.shared_control, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (
                row["model_id"] != "available_simple_controls"
                or row["head_id"] != "allele_identity_ridge"
            ):
                continue
            key = (int(row["seed"]), row["element_id"], row["assay_context_id"])
            if key in control:
                raise AblationError("duplicate shared-control row")
            control[key] = float(row["prediction"])
    if len(control) != 10330:
        raise AblationError("shared-control coverage differs")

    # Evaluation indexes elements in sorted order, as the frozen evaluator does.
    eval_elements = sorted(fit_elements)
    eval_index = {element: position for position, element in enumerate(eval_elements)}
    eval_order = np.asarray([eval_index[element] for element in fit_elements])
    eval_blocks = np.asarray([metadata[element][0] for element in eval_elements])
    observed = {
        context: np.asarray([outcome_data[(element, context)]["mean"] for element in eval_elements])
        for context in CONTEXTS
    }

    candidates: dict[tuple[str, str], dict[tuple[int, str], np.ndarray]] = {}
    control_by_seed: dict[tuple[int, str], np.ndarray] = {}
    for seed in SEEDS:
        for context in CONTEXTS:
            control_by_seed[(seed, context)] = np.asarray(
                [control[(seed, element, context)] for element in eval_elements]
            )
    candidates[("allele_identity_ridge", "control")] = control_by_seed

    arms = [("frozen_published_w4096", root / arguments.frozen_embeddings)]
    ablation_root = root / arguments.ablation_embeddings
    for arm_root in sorted(path for path in ablation_root.iterdir() if path.is_dir()):
        arms.append((arm_root.name, arm_root / "allele_embeddings.npz"))

    arm_receipts = []
    for arm_id, path in arms:
        embeddings, reorder = _load_arm(path, fixture_to_element, fit_elements)
        predictions = _arm_predictions(embeddings, reorder, folds, blocks, targets)
        for head in HEADS:
            by_seed = {}
            for seed in SEEDS:
                for context in CONTEXTS:
                    values = np.empty(1033)
                    values[eval_order] = predictions[(head, seed, context)]
                    by_seed[(seed, context)] = values
            candidates[(arm_id, head)] = by_seed
        arm_receipts.append({"arm_id": arm_id, "embeddings": str(path.relative_to(root))})
        print(f"fitted arm {arm_id}", flush=True)

    # Point estimates on the five-seed mean prediction, macro Fisher-z across contexts.
    ensemble = {
        candidate: {
            context: np.mean(np.vstack([by_seed[(seed, context)] for seed in SEEDS]), axis=0)
            for context in CONTEXTS
        }
        for candidate, by_seed in candidates.items()
    }
    points = {
        candidate: _macro(
            [
                float(spearmanr(observed[context], contexts[context]).statistic)
                for context in CONTEXTS
            ]
        )
        for candidate, contexts in ensemble.items()
    }
    constant = sorted(candidate for candidate, value in points.items() if value is None)
    if constant:
        raise AblationError(f"constant ensemble prediction, no Spearman defined: {constant}")
    per_seed = {
        (candidate, seed): _macro(
            [
                float(spearmanr(observed[context], by_seed[(seed, context)]).statistic)
                for context in CONTEXTS
            ]
        )
        for candidate, by_seed in candidates.items()
        for seed in SEEDS
    }
    constant_seeds = sorted(key for key, value in per_seed.items() if value is None)
    if constant_seeds:
        raise AblationError(f"constant per-seed prediction, no Spearman defined: {constant_seeds}")

    reference = ("allele_identity_ridge", "control")
    published = ("frozen_published_w4096", "delta_ridge")
    reproduction = {
        "frozen_arm_macro_spearman": points[published],
        "published_macro_spearman": PUBLISHED_SPEARMAN,
        "frozen_arm_gain": points[published] - points[reference],
        "published_gain": PUBLISHED_GAIN,
        "control_macro_spearman": points[reference],
        "spearman_abs_diff": abs(points[published] - PUBLISHED_SPEARMAN),
        "gain_abs_diff": abs((points[published] - points[reference]) - PUBLISHED_GAIN),
        "tolerance": REPRODUCTION_TOLERANCE,
    }
    reproduction["reproduced"] = bool(
        reproduction["spearman_abs_diff"] <= REPRODUCTION_TOLERANCE
        and reproduction["gain_abs_diff"] <= REPRODUCTION_TOLERANCE
    )
    extraction = None
    if ("w4096_proportional", "delta_ridge") in points:
        value = points[("w4096_proportional", "delta_ridge")]
        extraction = {
            "reextracted_macro_spearman": value,
            "published_macro_spearman": PUBLISHED_SPEARMAN,
            "abs_diff": abs(value - PUBLISHED_SPEARMAN),
            "tolerance": EXTRACTION_TOLERANCE,
            "reproduced": bool(abs(value - PUBLISHED_SPEARMAN) <= EXTRACTION_TOLERANCE),
        }

    output.mkdir(parents=True, mode=0o750)
    if not reproduction["reproduced"] and not arguments.continue_on_failed_reproduction:
        (output / "receipt.json").write_text(
            json.dumps(
                {
                    "status": "fail_frozen_reproduction",
                    "reproduction": reproduction,
                    "extraction_reproduction": extraction,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        raise AblationError(
            "D1 failed: refit of the frozen features did not reproduce the published score"
        )

    # Paired block bootstrap: identical draws for every candidate.
    order = list(candidates)
    samples = {candidate: np.full(BOOTSTRAP_RESAMPLES, np.nan) for candidate in order}
    unique = sorted(set(eval_blocks.tolist()))
    by_block = {block: np.flatnonzero(eval_blocks == block) for block in unique}
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    for iteration in range(BOOTSTRAP_RESAMPLES):
        sampled = rng.integers(0, len(unique), size=len(unique))
        indices = np.concatenate([by_block[unique[index]] for index in sampled])
        for candidate in order:
            contexts = ensemble[candidate]
            value = _macro(
                [
                    _fast_spearman(observed[context][indices], contexts[context][indices])
                    for context in CONTEXTS
                ]
            )
            if value is not None:
                samples[candidate][iteration] = value

    reference_samples = samples[reference]
    published_gain_samples = samples[published] - reference_samples
    rows = []
    for candidate in order:
        arm_id, head = candidate
        current = samples[candidate]
        valid = np.isfinite(current)
        low, high = (float(value) for value in np.quantile(current[valid], (0.025, 0.975)))
        gain_samples = current - reference_samples
        gain_low, gain_high = (
            float(value) for value in np.quantile(gain_samples[np.isfinite(gain_samples)], (0.025, 0.975))
        )
        # Paired contrast against the published 4,096 bp arm, on shared draws.
        versus = published_gain_samples - gain_samples
        versus_valid = np.isfinite(versus)
        versus_low, versus_high = (
            float(value) for value in np.quantile(versus[versus_valid], (0.025, 0.975))
        )
        rows.append(
            {
                "arm_id": arm_id,
                "head_id": head,
                "macro_spearman": format(points[candidate], ".17g"),
                "ci_low": format(low, ".17g"),
                "ci_high": format(high, ".17g"),
                "gain_vs_allele_identity_ridge": format(points[candidate] - points[reference], ".17g"),
                "gain_ci_low": format(gain_low, ".17g"),
                "gain_ci_high": format(gain_high, ".17g"),
                "gain_ci_excludes_zero": str(gain_low > 0.0 or gain_high < 0.0).lower(),
                "positive_gain_seeds": sum(
                    per_seed[(candidate, seed)] > per_seed[(reference, seed)] for seed in SEEDS
                ),
                "published_arm_gain_minus_this_gain": format(
                    (points[published] - points[reference]) - (points[candidate] - points[reference]),
                    ".17g",
                ),
                "published_minus_this_ci_low": format(versus_low, ".17g"),
                "published_minus_this_ci_high": format(versus_high, ".17g"),
                "published_minus_this_excludes_zero": str(
                    versus_low > 0.0 or versus_high < 0.0
                ).lower(),
                "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
                "bootstrap_unit": "long_range_block_id",
            }
        )
    _write_tsv(output / "arm_metrics.tsv", rows)
    _write_tsv(
        output / "per_seed_metrics.tsv",
        [
            {
                "arm_id": candidate[0],
                "head_id": candidate[1],
                "seed": seed,
                "macro_spearman": format(per_seed[(candidate, seed)], ".17g"),
            }
            for candidate in order
            for seed in SEEDS
        ],
    )
    np.savez_compressed(
        output / "bootstrap_samples.npz",
        **{f"{arm}__{head}": samples[(arm, head)] for arm, head in order},
    )

    receipt = {
        "schema_version": "masld-bench-gse281364-hyenadna-window-ablation-v1",
        "status": "pass_hyenadna_window_ablation"
        if reproduction["reproduced"]
        else "pass_with_failed_frozen_reproduction",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "model_id": "hyenadna",
        "elements": 1033,
        "long_range_blocks": 239,
        "outer_folds": 5,
        "seeds": list(SEEDS),
        "contexts": list(CONTEXTS),
        "heads": list(HEADS),
        "arms": arm_receipts,
        "frozen_reproduction": reproduction,
        "extraction_reproduction": extraction,
        "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_unit": "long_range_block_id",
        "shared_control_refit_performed": False,
        "shared_control_prediction_policy": "reuse_bit_exact_frozen_enformer_sei_oof",
        "biological_donors": 0,
        "experimental_replicates_per_context": 4,
        "mandatory_task_native_baselines_complete": False,
        "champion_claim": False,
        "cross_family_ranking": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--row-universe",
        default="executions/model-check-220-21088696/contract/row_universe.tsv",
    )
    parser.add_argument(
        "--outcomes",
        default="executions/gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz",
    )
    parser.add_argument(
        "--fixture-manifest",
        default="executions/gse281364-dna-lm-common-fixture-21069076/fixture/sequence_manifest.tsv",
    )
    parser.add_argument(
        "--shared-control",
        default="executions/model-cpu-train-604-21097189/fit/oof_predictions.tsv.gz",
    )
    parser.add_argument(
        "--frozen-embeddings",
        default="executions/model-work-001-21070132/raw/hyenadna/allele_embeddings.npz",
    )
    parser.add_argument("--ablation-embeddings", required=True)
    parser.add_argument("--continue-on-failed-reproduction", action="store_true")
    run(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Fit an enriched local-sequence allele baseline ladder for GSE281364 MPRA.

Adversarial re-test of the frozen ``variant_to_regulation`` claim. The frozen
control (``available_simple_controls/allele_identity_ridge``) encodes only a
16-column one-hot of (REF base, ALT base). This lane fits a nested ladder of
progressively richer *local sequence* controls and re-scores the frozen HyenaDNA
delta head against each one.

Fitting and scoring reuse the production functions rather than reimplementing
them, so folds, seeds, bootstrap salts, ridge selection, the Fisher-z macro
endpoint and the 10,000 paired-block bootstrap are the same code that produced
``executions/model-cpu-train-605-21099008``. HyenaDNA is never refitted; its
frozen out-of-fold rows are read as-is.

Behaviour is fixed by ``executions/mpra-enriched-baseline-prespec-*/PRESPEC.md``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import itertools
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes
from scripts.evaluate_gse281364_dna_language_seeded_heads import (
    _bootstrap,
    _fast_spearman,
    _macro,
    _metrics,
)
from scripts.fit_gse281364_enformer_sei_heads import fit_seeded_outer_fold


SCHEMA = "masld-bench-gse281364-enriched-allele-baseline-ladder-v1"
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
SEEDS = (1103, 2909, 4721, 6673, 8111)
ALPHAS = (0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)
VARIANT_INDEX0 = 62
ELEMENT_LENGTH_BP = 107
BOOTSTRAP_RESAMPLES = 10000
BOOTSTRAP_SEED = 20260825

FROZEN = {
    "control_primary_score": 0.13421551707194201,
    "hyenadna_primary_score": 0.21622332457761556,
    "hyenadna_gain": 0.082007807505673547,
    "hyenadna_gain_ci_low": 0.027427479204226245,
    "hyenadna_gain_ci_high": 0.13891768494934109,
}
GUARD_TOLERANCE = 1.0e-12

STATIC_MEMBER = "executions/model-check-274-21097044/binding/base_static_scores.tsv"
ELEMENTS_MEMBER = (
    "executions/gse281364-validated-reconstruction-21066470/validated/elements.tsv"
)
ROW_UNIVERSE_MEMBER = "executions/model-check-220-21088696/contract/row_universe.tsv"
OUTCOME_MEMBER = (
    "executions/gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz"
)
FROZEN_FIT_MEMBER = (
    "executions/model-cpu-train-605-21099008/fit/oof_predictions.tsv.gz"
)
OUTCOME_SHA256 = "783caf7e7b659d95d34f5227e3aa149f636faf92ea5d685c4a25b41cb1f0b8f1"

BASES = "ACGT"
DINUCLEOTIDES = tuple("".join(pair) for pair in itertools.product(BASES, repeat=2))
# The variant sits at index 62 of a 107 bp element, so it has 44 bases to its
# right: a centred window cannot exceed 2 * 44 + 1 = 89 bp. The widest
# composition feature is therefore the whole element rather than a centred
# window, which is also what the reporter construct actually carries.
# See PRESPEC_AMENDMENT_01.md.
GC_WINDOWS = (11, 21, 51)
FLANK_WINDOW = 21


class LadderError(RuntimeError):
    """Raised when a frozen identity, guard or data contract is not met."""


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_tsv(path: Path, *, compressed: bool = False) -> list[dict[str, str]]:
    if compressed:
        handle = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        handle = path.open("r", encoding="utf-8", newline="")
    with handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def _write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _display(value: float | None) -> str:
    return "not_applicable" if value is None else format(value, ".17g")


# ---------------------------------------------------------------------------
# Feature ladder
# ---------------------------------------------------------------------------


def _gc_fraction(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def _cpg_count(sequence: str) -> int:
    return sum(1 for i in range(len(sequence) - 1) if sequence[i : i + 2] == "CG")


def _has_cpg_at_variant(sequence: str) -> bool:
    p = VARIANT_INDEX0
    return sequence[p - 1 : p + 1] == "CG" or sequence[p : p + 2] == "CG"


def build_feature_blocks(
    sequences: Mapping[str, tuple[str, str]], elements: Sequence[str]
) -> tuple[dict[str, np.ndarray], dict[str, list[str]]]:
    """Return the four nested rungs as (matrix, column-name) pairs.

    B0 reproduces the frozen 16-column REF/ALT one-hot exactly, including the
    four constant-zero REF==ALT columns that the production standardizer drops.
    """
    p = VARIANT_INDEX0
    sub_names = [f"sub_{ref}{alt}" for ref in BASES for alt in BASES]
    cpg_names = ["cpg_created", "cpg_destroyed"]
    composition_names: list[str] = []
    for window in GC_WINDOWS:
        composition_names.append(f"gc_ref_w{window}")
        composition_names.append(f"cpg_ref_w{window}")
    composition_names.append("gc_ref_full107")
    composition_names.append("cpg_ref_full107")
    flank_names = [f"flank5_{base}" for base in BASES] + [
        f"flank3_{base}" for base in BASES
    ]
    di_names = [f"di21_{pair}" for pair in DINUCLEOTIDES]

    base_index = {base: position for position, base in enumerate(BASES)}
    sub_rows: list[np.ndarray] = []
    cpg_rows: list[np.ndarray] = []
    composition_rows: list[np.ndarray] = []
    flank_rows: list[np.ndarray] = []
    di_rows: list[np.ndarray] = []

    for element in elements:
        reference, alternate = sequences[element]
        if len(reference) != ELEMENT_LENGTH_BP or len(alternate) != ELEMENT_LENGTH_BP:
            raise LadderError("element sequence length differs")
        differing = [i for i in range(ELEMENT_LENGTH_BP) if reference[i] != alternate[i]]
        if differing != [p]:
            raise LadderError("element is not a single substitution at the frozen index")
        ref_base, alt_base = reference[p], alternate[p]
        if ref_base not in base_index or alt_base not in base_index or ref_base == alt_base:
            raise LadderError("allele identity differs")

        substitution = np.zeros(16, dtype=np.float64)
        substitution[4 * base_index[ref_base] + base_index[alt_base]] = 1.0
        sub_rows.append(substitution)

        reference_cpg = _has_cpg_at_variant(reference)
        alternate_cpg = _has_cpg_at_variant(alternate)
        cpg_rows.append(
            np.asarray(
                [
                    float(alternate_cpg and not reference_cpg),
                    float(reference_cpg and not alternate_cpg),
                ],
                dtype=np.float64,
            )
        )

        composition: list[float] = []
        for window in GC_WINDOWS:
            half = window // 2
            span = reference[p - half : p + half + 1]
            if len(span) != window:
                raise LadderError("composition window truncated")
            composition.append(_gc_fraction(span))
            composition.append(float(_cpg_count(span)))
        composition.append(_gc_fraction(reference))
        composition.append(float(_cpg_count(reference)))
        composition_rows.append(np.asarray(composition, dtype=np.float64))

        flank = [float(reference[p - 1] == base) for base in BASES]
        flank += [float(reference[p + 1] == base) for base in BASES]
        flank_rows.append(np.asarray(flank, dtype=np.float64))

        half = FLANK_WINDOW // 2
        span = reference[p - half : p + half + 1]
        pairs = [span[i : i + 2] for i in range(len(span) - 1)]
        denominator = float(len(pairs))
        di_rows.append(
            np.asarray(
                [pairs.count(pair) / denominator for pair in DINUCLEOTIDES],
                dtype=np.float64,
            )
        )

    substitution_matrix = np.vstack(sub_rows)
    cpg_matrix = np.vstack(cpg_rows)
    composition_matrix = np.vstack(composition_rows)
    flank_matrix = np.vstack(flank_rows)
    di_matrix = np.vstack(di_rows)

    b0 = substitution_matrix
    b1 = np.hstack([b0, cpg_matrix])
    b2 = np.hstack([b1, composition_matrix])
    b3 = np.hstack([b2, flank_matrix, di_matrix])

    features = {
        "allele_identity_ridge_b0": b0,
        "enriched_b1_cpg": b1,
        "enriched_b2_composition": b2,
        "enriched_b3_flank_dinucleotide": b3,
    }
    names = {
        "allele_identity_ridge_b0": list(sub_names),
        "enriched_b1_cpg": sub_names + cpg_names,
        "enriched_b2_composition": sub_names + cpg_names + composition_names,
        "enriched_b3_flank_dinucleotide": (
            sub_names + cpg_names + composition_names + flank_names + di_names
        ),
    }
    for rung, matrix in features.items():
        if matrix.shape != (len(elements), len(names[rung])):
            raise LadderError(f"feature shape differs for {rung}")
        if not np.isfinite(matrix).all():
            raise LadderError(f"non-finite feature for {rung}")
    return features, names


def _top_k_grid(width: int) -> tuple[int, ...]:
    """Prespecified grid [12, 24, 48, W], truncated and de-duplicated."""
    return tuple(sorted({value for value in (12, 24, 48, width) if 0 < value <= width}))


# ---------------------------------------------------------------------------
# Lane
# ---------------------------------------------------------------------------


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.root.resolve(strict=True)
    output = arguments.output
    if output.exists():
        raise LadderError("output already exists; refusing to overwrite")

    prespec = arguments.prespec.resolve(strict=True)
    prespec_sha256 = file_sha256(prespec)
    if prespec_sha256 != arguments.prespec_sha256:
        raise LadderError("prespec identity differs from the sealed digest")

    outcome_path = root / OUTCOME_MEMBER
    if file_sha256(outcome_path) != OUTCOME_SHA256:
        raise LadderError("outcome authority identity differs")

    # ---- fit-side universe: the frozen fit ordering, not the evaluator ordering
    base_rows = _read_tsv(root / STATIC_MEMBER)
    base_rows.sort(key=lambda row: int(row["enformer_feature_index0"]))
    fit_elements = [row["element_id"] for row in base_rows]
    if len(fit_elements) != 1033 or len(set(fit_elements)) != 1033:
        raise LadderError("element denominator differs")
    folds = np.asarray(
        [int(row["outer_fold"].removeprefix("fold-")) for row in base_rows], dtype=np.int64
    )
    fit_blocks = np.asarray([row["long_range_block_id"] for row in base_rows])
    if set(folds.tolist()) != set(range(5)) or len(set(fit_blocks.tolist())) != 239:
        raise LadderError("fold or block census differs")
    for block in set(fit_blocks.tolist()):
        if len(set(folds[fit_blocks == block].tolist())) != 1:
            raise LadderError("long-range block crosses outer folds")

    wanted = set(fit_elements)
    sequences: dict[str, tuple[str, str]] = {}
    for row in _read_tsv(root / ELEMENTS_MEMBER):
        if row["element_id"] not in wanted:
            continue
        if row["element_id"] in sequences:
            raise LadderError("duplicate element sequence")
        if row["pair_state"] != "paired_snv":
            raise LadderError("element is not a paired SNV")
        sequences[row["element_id"]] = (
            row["ref_sequence_107bp"],
            row["alt_sequence_107bp"],
        )
    if set(sequences) != wanted:
        raise LadderError("sequence coverage differs")

    features, feature_names = build_feature_blocks(sequences, fit_elements)
    outcome_data = load_outcomes(outcome_path, wanted)
    targets = {
        context: np.asarray(
            [outcome_data[(element, context)]["mean"] for element in fit_elements],
            dtype=np.float64,
        )
        for context in CONTEXTS
    }

    output.mkdir(parents=True, mode=0o750)

    # ---- fit every rung on the frozen protocol
    predictions: dict[tuple[str, int, str], np.ndarray] = {}
    selection_rows: list[dict[str, Any]] = []
    for rung, matrix in features.items():
        grid = _top_k_grid(matrix.shape[1])
        for context in CONTEXTS:
            y = targets[context]
            for seed in SEEDS:
                oof = np.empty(1033, dtype=np.float64)
                for held_fold in range(5):
                    prediction, selection, _state = fit_seeded_outer_fold(
                        features=matrix,
                        outcomes=y,
                        folds=folds,
                        block_ids=fit_blocks,
                        held_fold=held_fold,
                        seed=seed,
                        top_k_grid=grid,
                        alphas=ALPHAS,
                    )
                    oof[folds == held_fold] = prediction
                    selection_rows.append(
                        {
                            "rung": rung,
                            "feature_columns": matrix.shape[1],
                            "top_k_grid": ",".join(str(value) for value in grid),
                            "assay_context_id": context,
                            **selection,
                        }
                    )
                if not np.isfinite(oof).all():
                    raise LadderError(f"OOF coverage differs for {rung}")
                predictions[(rung, seed, context)] = oof

    fit_index = {element: position for position, element in enumerate(fit_elements)}

    # ---- G1: bit-exact reproduction of the frozen allele-identity control
    frozen_rows = _read_tsv(root / FROZEN_FIT_MEMBER, compressed=True)
    frozen_control: dict[tuple[int, str, str], float] = {}
    frozen_hyenadna: dict[tuple[int, str, str], float] = {}
    for row in frozen_rows:
        key = (int(row["seed"]), row["assay_context_id"], row["element_id"])
        if row["model_id"] == "available_simple_controls" and row["head_id"] == "allele_identity_ridge":
            frozen_control[key] = float(row["prediction"])
        elif row["model_id"] == "hyenadna" and row["head_id"] == "delta_ridge":
            frozen_hyenadna[key] = float(row["prediction"])
    if len(frozen_control) != 10330 or len(frozen_hyenadna) != 10330:
        raise LadderError("frozen prediction denominator differs")

    control_deltas = []
    for seed in SEEDS:
        for context in CONTEXTS:
            refit = predictions[("allele_identity_ridge_b0", seed, context)]
            for element in fit_elements:
                control_deltas.append(
                    abs(refit[fit_index[element]] - frozen_control[(seed, context, element)])
                )
    g1_max_absolute_difference = float(max(control_deltas))
    g1_pass = g1_max_absolute_difference <= GUARD_TOLERANCE

    # ---- evaluation-side universe: the frozen evaluator ordering
    row_rows = _read_tsv(root / ROW_UNIVERSE_MEMBER)
    if len(row_rows) != 10330:
        raise LadderError("row universe denominator differs")
    eval_elements = sorted({row["element_id"] for row in row_rows})
    block_of: dict[str, str] = {}
    for row in row_rows:
        if block_of.setdefault(row["element_id"], row["long_range_block_id"]) != row["long_range_block_id"]:
            raise LadderError("block metadata differs")
    eval_blocks = np.asarray([block_of[element] for element in eval_elements])
    if len(eval_elements) != 1033 or len(set(eval_blocks.tolist())) != 239:
        raise LadderError("evaluation element/block denominator differs")
    observed = {
        context: np.asarray(
            [outcome_data[(element, context)]["mean"] for element in eval_elements],
            dtype=np.float64,
        )
        for context in CONTEXTS
    }

    def _to_eval(vector: np.ndarray) -> np.ndarray:
        return np.asarray(
            [vector[fit_index[element]] for element in eval_elements], dtype=np.float64
        )

    per_seed: dict[tuple[str, int], np.ndarray] = {}
    for rung in features:
        for seed in SEEDS:
            for context in CONTEXTS:
                per_seed[(rung, seed, context)] = _to_eval(
                    predictions[(rung, seed, context)]
                )
    for seed in SEEDS:
        for context in CONTEXTS:
            per_seed[("hyenadna_delta_ridge", seed, context)] = np.asarray(
                [frozen_hyenadna[(seed, context, element)] for element in eval_elements],
                dtype=np.float64,
            )

    candidates = tuple(features) + ("hyenadna_delta_ridge",)
    ensemble: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    points: dict[tuple[str, str], float | None] = {}
    per_seed_primary: dict[tuple[str, int], float | None] = {}
    per_seed_rows: list[dict[str, Any]] = []
    ensemble_rows: list[dict[str, Any]] = []

    for candidate in candidates:
        key = (candidate, "ensemble")
        for seed in SEEDS:
            context_metrics = []
            for context in CONTEXTS:
                metrics = _metrics(observed[context], per_seed[(candidate, seed, context)])
                context_metrics.append(metrics)
                per_seed_rows.append(
                    {
                        "candidate": candidate,
                        "seed": seed,
                        "assay_context_id": context,
                        "elements": 1033,
                        **{name: _display(value) for name, value in metrics.items()},
                    }
                )
            primary = _macro([item["spearman"] for item in context_metrics])
            per_seed_primary[(candidate, seed)] = primary
            per_seed_rows.append(
                {
                    "candidate": candidate,
                    "seed": seed,
                    "assay_context_id": "macro",
                    "elements": 2066,
                    "spearman": _display(primary),
                }
            )
        by_context = {
            context: np.mean(
                np.vstack([per_seed[(candidate, seed, context)] for seed in SEEDS]), axis=0
            )
            for context in CONTEXTS
        }
        ensemble[key] = by_context
        context_metrics = []
        for context in CONTEXTS:
            metrics = _metrics(observed[context], by_context[context])
            context_metrics.append(metrics)
            ensemble_rows.append(
                {
                    "candidate": candidate,
                    "assay_context_id": context,
                    "elements": 1033,
                    **{name: _display(value) for name, value in metrics.items()},
                }
            )
        points[key] = _macro([item["spearman"] for item in context_metrics])
        ensemble_rows.append(
            {
                "candidate": candidate,
                "assay_context_id": "macro",
                "elements": 2066,
                "spearman": _display(points[key]),
            }
        )

    samples = _bootstrap(
        observed,
        ensemble,
        eval_blocks,
        resamples=BOOTSTRAP_RESAMPLES,
        seed=BOOTSTRAP_SEED,
    )

    hyenadna_key = ("hyenadna_delta_ridge", "ensemble")
    hyenadna_samples = samples[hyenadna_key]
    gain_rows: list[dict[str, Any]] = []
    for rung in features:
        key = (rung, "ensemble")
        reference_samples = samples[key]
        gain_samples = hyenadna_samples - reference_samples
        valid = np.isfinite(gain_samples)
        if int(valid.sum()) < 9500:
            gain_low = gain_high = None
        else:
            gain_low, gain_high = (
                float(value) for value in np.quantile(gain_samples[valid], (0.025, 0.975))
            )
        reference_valid = np.isfinite(reference_samples)
        if int(reference_valid.sum()) < 9500:
            reference_low = reference_high = reference_se = None
        else:
            reference_low, reference_high = (
                float(value)
                for value in np.quantile(reference_samples[reference_valid], (0.025, 0.975))
            )
            reference_se = float(np.std(reference_samples[reference_valid], ddof=1))
        positive_seeds = sum(
            per_seed_primary[("hyenadna_delta_ridge", seed)] is not None
            and per_seed_primary[(rung, seed)] is not None
            and per_seed_primary[("hyenadna_delta_ridge", seed)] > per_seed_primary[(rung, seed)]
            for seed in SEEDS
        )
        gain_rows.append(
            {
                "reference_rung": rung,
                "reference_feature_columns": features[rung].shape[1],
                "reference_primary_score": _display(points[key]),
                "reference_ci_low": _display(reference_low),
                "reference_ci_high": _display(reference_high),
                "reference_bootstrap_se": _display(reference_se),
                "hyenadna_primary_score": _display(points[hyenadna_key]),
                "gain_vs_reference": _display(
                    None
                    if points[key] is None or points[hyenadna_key] is None
                    else float(points[hyenadna_key] - points[key])
                ),
                "gain_ci_low": _display(gain_low),
                "gain_ci_high": _display(gain_high),
                "gain_ci_excludes_zero": "not_applicable"
                if gain_low is None
                else str(bool(gain_low > 0.0 or gain_high < 0.0)).lower(),
                "positive_gain_seeds": positive_seeds,
                "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
                "bootstrap_unit": "long_range_block_id",
            }
        )

    # ---- D3: is the enriched control measuring the same thing as HyenaDNA?
    agreement_rows = []
    for rung in features:
        for context in CONTEXTS:
            agreement_rows.append(
                {
                    "rung": rung,
                    "assay_context_id": context,
                    "spearman_vs_hyenadna_ensemble": _display(
                        _fast_spearman(
                            ensemble[(rung, "ensemble")][context],
                            ensemble[hyenadna_key][context],
                        )
                    ),
                }
            )

    # ---- G2/G3/G4
    b0_key = ("allele_identity_ridge_b0", "ensemble")
    b0_gain_row = next(row for row in gain_rows if row["reference_rung"] == "allele_identity_ridge_b0")
    guards = {
        "G1_control_oof_bit_exact": {
            "max_absolute_difference": g1_max_absolute_difference,
            "tolerance": GUARD_TOLERANCE,
            "pass": bool(g1_pass),
        },
        "G2_hyenadna_primary_score": {
            "observed": points[hyenadna_key],
            "frozen": FROZEN["hyenadna_primary_score"],
            "absolute_difference": abs(float(points[hyenadna_key]) - FROZEN["hyenadna_primary_score"]),
            "pass": bool(
                abs(float(points[hyenadna_key]) - FROZEN["hyenadna_primary_score"]) <= GUARD_TOLERANCE
            ),
        },
        "G3_control_primary_score": {
            "observed": points[b0_key],
            "frozen": FROZEN["control_primary_score"],
            "absolute_difference": abs(float(points[b0_key]) - FROZEN["control_primary_score"]),
            "pass": bool(
                abs(float(points[b0_key]) - FROZEN["control_primary_score"]) <= GUARD_TOLERANCE
            ),
        },
        "G4_control_gain_interval": {
            "observed_low": float(b0_gain_row["gain_ci_low"]),
            "observed_high": float(b0_gain_row["gain_ci_high"]),
            "frozen_low": FROZEN["hyenadna_gain_ci_low"],
            "frozen_high": FROZEN["hyenadna_gain_ci_high"],
            "absolute_difference_low": abs(
                float(b0_gain_row["gain_ci_low"]) - FROZEN["hyenadna_gain_ci_low"]
            ),
            "absolute_difference_high": abs(
                float(b0_gain_row["gain_ci_high"]) - FROZEN["hyenadna_gain_ci_high"]
            ),
            "pass": bool(
                abs(float(b0_gain_row["gain_ci_low"]) - FROZEN["hyenadna_gain_ci_low"]) <= GUARD_TOLERANCE
                and abs(float(b0_gain_row["gain_ci_high"]) - FROZEN["hyenadna_gain_ci_high"])
                <= GUARD_TOLERANCE
            ),
        },
        "G5_finite_oof_coverage": {
            "rows_per_rung": 10330,
            "pass": True,
        },
    }

    _write_tsv(
        output / "head_selection.tsv",
        list(selection_rows[0].keys()),
        selection_rows,
    )
    _write_tsv(
        output / "per_seed_metrics.tsv",
        [
            "candidate",
            "seed",
            "assay_context_id",
            "elements",
            "spearman",
            "pearson",
            "rmse",
            "mae",
            "r2",
            "calibration_intercept",
            "calibration_slope",
        ],
        per_seed_rows,
    )
    _write_tsv(
        output / "ensemble_metrics.tsv",
        [
            "candidate",
            "assay_context_id",
            "elements",
            "spearman",
            "pearson",
            "rmse",
            "mae",
            "r2",
            "calibration_intercept",
            "calibration_slope",
        ],
        ensemble_rows,
    )
    _write_tsv(output / "ladder_gains.tsv", list(gain_rows[0].keys()), gain_rows)
    _write_tsv(
        output / "hyenadna_agreement.tsv",
        list(agreement_rows[0].keys()),
        agreement_rows,
    )

    prediction_fields = ("rung", "seed", "assay_context_id", "element_id", "prediction")
    prediction_rows = [
        {
            "rung": rung,
            "seed": seed,
            "assay_context_id": context,
            "element_id": element,
            "prediction": format(
                float(predictions[(rung, seed, context)][fit_index[element]]), ".17g"
            ),
        }
        for rung in features
        for seed in SEEDS
        for context in CONTEXTS
        for element in fit_elements
    ]
    with gzip.open(output / "oof_predictions.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(prediction_fields), delimiter="\t")
        writer.writeheader()
        for row in prediction_rows:
            writer.writerow(row)

    (output / "feature_columns.json").write_text(
        json.dumps(feature_names, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_enriched_allele_baseline_ladder"
        if all(guard["pass"] for guard in guards.values())
        else "fail_guard",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "prespec_path": str(prespec.relative_to(root)),
        "prespec_sha256": prespec_sha256,
        "elements": 1033,
        "long_range_blocks": 239,
        "outer_folds": 5,
        "fixed_seeds": list(SEEDS),
        "contexts": list(CONTEXTS),
        "variant_index0": VARIANT_INDEX0,
        "rungs": {rung: features[rung].shape[1] for rung in features},
        "top_k_grids": {rung: list(_top_k_grid(features[rung].shape[1])) for rung in features},
        "alpha_grid": list(ALPHAS),
        "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_unit": "long_range_block_id",
        "primary_metric": "tanh_mean_Fisher_z_Spearman_across_contexts",
        "hyenadna_refitted": False,
        "hyenadna_source": FROZEN_FIT_MEMBER,
        "production_functions_reused": [
            "scripts.fit_gse281364_enformer_sei_heads.fit_seeded_outer_fold",
            "scripts.evaluate_gse281364_dna_language_seeded_heads._bootstrap",
            "scripts.evaluate_gse281364_dna_language_seeded_heads._macro",
            "scripts.evaluate_gse281364_dna_language_seeded_heads._metrics",
            "scripts.evaluate_gse281364_dna_lm_common_lane.load_outcomes",
        ],
        "guards": guards,
        "biological_donors": 0,
        "experimental_replicates_per_context": 4,
        "experimental_replicates_used_as_independent_donors": False,
        "outcome_role": "exposed_development_MPRA_only",
        "champion_claim": False,
        "shortlist_created": False,
        "external_evaluation": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--prespec", type=Path, required=True)
    parser.add_argument("--prespec-sha256", type=str, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    receipt = run(arguments)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if receipt["status"] != "pass_enriched_allele_baseline_ladder":
        raise SystemExit("guard failed")


if __name__ == "__main__":
    main()

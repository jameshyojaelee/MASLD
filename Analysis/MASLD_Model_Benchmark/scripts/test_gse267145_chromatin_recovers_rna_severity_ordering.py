#!/usr/bin/env python3
"""Does H3K27ac place a GSE267145 sample on the same severity ordering bulk RNA places it on?

Development-grade, single-cohort, participant-held-out. GSE267145 is project-exposed
downstream_demo substrate (champion_eligible=false); every number produced here is
development-grade and carries no external-transfer or clinical claim.

Target definition (declared before any chromatin fit):
  RNA severity position = the frozen five-seed out-of-fold prediction of
  nash_crn_component_sum by rna_hvg_pca_elastic_net (model-scoring-078). It is the
  supervised RNA model's own prediction, computed on this cohort inside the frozen
  participant-held-out folds, and it never saw H3K27ac.
  Fibrosis is excluded (97.1% binary here). stage3 is excluded because it is nearly
  the binary disease call that floor (i) must separate from.

H3K27ac regions are used only as opaque feature keys. No sequence, liftover, overlap
or motif work. RNA values are fractional and are never rounded; no integer likelihood.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import platform
import sys
from hashlib import sha256
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from masld_bench.artifacts import verify_frozen_tree  # noqa: E402
from masld_bench.hashing import canonical_sha256  # noqa: E402

BASELINES_PATH = ROOT / "scripts" / "fit_gse267145_histology_baselines.py"
_spec = importlib.util.spec_from_file_location("gse267145_baselines", BASELINES_PATH)
BASE = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(BASE)

# ---- declared, pre-fit configuration -------------------------------------------------
RNA_FEATURE_REQUEST = 5000        # from the frozen rna lane candidate_counts
H3_FEATURE_REQUEST = 20000        # from the frozen h3 lane candidate_counts
PCA_COMPONENTS = 20               # from the frozen pca_components grid
ALPHAS = [0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
L1_RATIOS = [0.0, 0.5, 1.0]
POSITION_CLIP = (0.0, 8.0)        # nash_crn_component_sum source range
BINARY_CLIP = (0.0, 1.0)
RESIDUAL_CLIP = (-8.0, 8.0)
RANDOM_REGION_DRAWS = 50
RANDOM_REGION_LADDER = [100, 500, 2000]
RANDOM_LADDER_DRAWS = 20
RANDOM_SCRAMBLED_DRAWS = 10
PERMUTATIONS = 200
FIT_SEED = 20260828
BOOTSTRAP_SEED = 20260824         # the GSE267145 frozen scoring seed, not 20260828
BOOTSTRAP_REPLICATES = 10000
BOOTSTRAP_INDICES_SHA256 = "fb7ea4b28ef37ee759a70bb3664868f6affa5accc8dadd96edf202238536b2c5"


def rank(values: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata

    return np.asarray(rankdata(values, method="average"), dtype=np.float64)


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra, rb = rank(a), rank(b)
    if ra.std() == 0.0 or rb.std() == 0.0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def main() -> None:
    global PERMUTATIONS, RANDOM_REGION_DRAWS, BOOTSTRAP_REPLICATES
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture", required=True)
    parser.add_argument("--activation", required=True)
    parser.add_argument("--fit-views", required=True)
    parser.add_argument("--scoring", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--permutations", type=int, default=PERMUTATIONS)
    parser.add_argument("--random-draws", type=int, default=RANDOM_REGION_DRAWS)
    parser.add_argument("--bootstrap-replicates", type=int, default=BOOTSTRAP_REPLICATES)
    args = parser.parse_args()

    PERMUTATIONS = int(args.permutations)
    RANDOM_REGION_DRAWS = int(args.random_draws)
    BOOTSTRAP_REPLICATES = int(args.bootstrap_replicates)

    fixture = Path(args.fixture).resolve(strict=True)
    activation = Path(args.activation).resolve(strict=True)
    fit_views = Path(args.fit_views).resolve(strict=True)
    scoring = Path(args.scoring).resolve(strict=True)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)

    provenance = {}
    for name, path in (
        ("fixture", fixture),
        ("fixture_molecular", fixture / "molecular"),
        ("fixture_folds", fixture / "folds"),
        ("activation", activation),
        ("scoring", scoring),
    ):
        verify_frozen_tree(path)
        provenance[name] = {
            "path": str(path),
            "artifacts_json_sha256": sha256((path / "ARTIFACTS.json").read_bytes()).hexdigest(),
        }

    # ---- axes -----------------------------------------------------------------------
    molecular = fixture / "molecular"
    axis_rows = read_tsv(molecular / "participant_axis.tsv")
    fold_rows = read_tsv(fixture / "folds" / "participant_outer_folds.tsv")
    participants = [row["participant_id"] for row in axis_rows]
    if participants != [row["participant_id"] for row in fold_rows]:
        raise SystemExit("participant order differs between molecular axis and folds")
    n = len(participants)
    if n != 99:
        raise SystemExit("participant count differs")
    index_of = {pid: i for i, pid in enumerate(participants)}
    outer_fold = np.asarray([int(row["outer_fold"]) for row in fold_rows], dtype=int)

    rna_raw = np.load(molecular / "rna_values.npy", allow_pickle=False)
    h3_raw = np.load(molecular / "h3k27ac_counts.npy", allow_pickle=False)
    if rna_raw.shape != (99, 42163) or h3_raw.shape != (99, 96460):
        raise SystemExit("molecular shapes differ")
    if not np.any(rna_raw != np.floor(rna_raw)):
        raise SystemExit("RNA values are not fractional; fixture identity differs")
    rna_ids = [r["stable_gene_id"] for r in read_tsv(molecular / "rna_feature_axis.tsv")]
    h3_ids = [r["opaque_source_feature_key"] for r in read_tsv(molecular / "h3k27ac_feature_axis.tsv")]

    h3_library = np.asarray(h3_raw, dtype=np.float64).sum(axis=1)
    rna = BASE._log_library_scale(rna_raw)
    h3 = BASE._log_library_scale(h3_raw)

    # ---- outcomes (evaluator side only) ---------------------------------------------
    endpoint_rows = {r["participant_id"]: r for r in read_tsv(activation / "participant_endpoints.tsv")}
    nas_true = np.asarray([float(endpoint_rows[p]["nash_crn_component_sum"]) for p in participants])
    stage3 = [endpoint_rows[p]["stage3"] for p in participants]
    binary = np.asarray([0.0 if s == "NOR" else 1.0 for s in stage3])
    fibrosis = np.asarray([float(endpoint_rows[p]["fibrosis"]) for p in participants])

    # ---- frozen predictions ---------------------------------------------------------
    bundles = scoring / "scores" / "prediction_bundles"

    def frozen(model_id: str, endpoint: str, column: str) -> np.ndarray:
        rows = read_tsv(bundles / f"{model_id}--endpoint-{endpoint}" / "predictions.tsv")
        if [r["participant_id"] for r in rows] != participants:
            raise SystemExit(f"frozen bundle order differs: {model_id}/{endpoint}")
        return np.asarray([float(r[column]) for r in rows], dtype=np.float64)

    rna_position = frozen("rna_hvg_pca_elastic_net", "nash_crn_component_sum",
                          "predicted_nash_crn_component_sum")
    h3_position_frozen = frozen("h3_variance_pca_elastic_net", "nash_crn_component_sum",
                                "predicted_nash_crn_component_sum")
    h3_stage_nor = frozen("h3_variance_pca_elastic_net", "stage3", "probability_NOR")
    h3_disease_score_frozen = 1.0 - h3_stage_nor
    rna_stage_nor = frozen("rna_hvg_pca_elastic_net", "stage3", "probability_NOR")

    # ---- bootstrap indices (frozen, GSE267145 scoring seed) --------------------------
    boot = np.random.default_rng(BOOTSTRAP_SEED).integers(
        0, n, size=(BOOTSTRAP_REPLICATES, n), endpoint=False
    )
    boot_digest = canonical_sha256(boot.tolist())
    if BOOTSTRAP_REPLICATES == 10000 and boot_digest != BOOTSTRAP_INDICES_SHA256:
        raise SystemExit(f"frozen bootstrap index digest differs: {boot_digest}")
    distinct = np.asarray([len(np.unique(row)) for row in boot], dtype=float)

    # ---- fold-honest refits ----------------------------------------------------------
    rna_position_refit = np.full(n, np.nan)
    h3_to_position = np.full(n, np.nan)
    h3_to_binary = np.full(n, np.nan)
    h3_to_residual = np.full(n, np.nan)
    h3_libsum_fit = np.full(n, np.nan)
    residual_target_eval = np.full(n, np.nan)
    random_region_pred = np.full((RANDOM_REGION_DRAWS, n), np.nan)
    random_ladder_pred = {
        size: np.full((RANDOM_LADDER_DRAWS, n), np.nan) for size in RANDOM_REGION_LADDER
    }
    random_scrambled_pred = np.full((RANDOM_SCRAMBLED_DRAWS, n), np.nan)
    permuted_pred = np.full((PERMUTATIONS, n), np.nan)
    fold_receipts = []

    for k in range(5):
        train_rows = read_tsv(fit_views / f"outer_{k}" / "training_endpoints.tsv")
        query_rows = read_tsv(fit_views / f"outer_{k}" / "query_participants.tsv")
        train = np.asarray([index_of[r["participant_id"]] for r in train_rows], dtype=int)
        test = np.asarray([index_of[r["participant_id"]] for r in query_rows], dtype=int)
        if np.any(outer_fold[train] == k) or np.any(outer_fold[test] != k):
            raise SystemExit(f"fit view crosses outer fold {k}")
        inner = np.asarray([int(r["inner_validation_fold"]) for r in train_rows], dtype=int)
        nas_train = np.asarray([float(r["nash_crn_component_sum"]) for r in train_rows])
        binary_train = binary[train]

        base_seed = FIT_SEED + k * 100_000
        # RepresentationCache adds a 32-bit sha256 offset to its seed before handing it
        # to PCA(random_state=...), so the cache seed itself must stay small or sklearn
        # rejects the resulting value. Small deterministic cache seeds only.
        rna_cache = BASE.RepresentationCache(
            rna, rna_ids, [RNA_FEATURE_REQUEST], [PCA_COMPONENTS], 1000 + k
        )
        h3_cache = BASE.RepresentationCache(
            h3, h3_ids, [H3_FEATURE_REQUEST], [PCA_COMPONENTS], 2000 + k
        )
        rna_rep = {"feature_request": RNA_FEATURE_REQUEST, "pca_components": PCA_COMPONENTS}
        h3_rep = {"feature_request": H3_FEATURE_REQUEST, "pca_components": PCA_COMPONENTS}

        def tune(cache, rep, target, clip, seed, maximize=True):
            return BASE.tune_regression(
                endpoint_id="ordering",
                cache=cache,
                representation=rep,
                outer_training_indices=train,
                outer_test_indices=test,
                inner_assignment=inner,
                target=np.asarray(target, dtype=np.float64),
                alphas=ALPHAS,
                l1_ratios=L1_RATIOS,
                maximize=maximize,
                score=BASE.spearman,
                clip=clip,
                seed=seed,
            )

        # Step 1: leak-free fold-specific RNA severity position.
        # tune_regression returns inner out-of-fold predictions on the training rows and
        # outer-test predictions; both use only fold-k training labels.
        rna_fit = tune(rna_cache, rna_rep, nas_train, POSITION_CLIP, base_seed + 1)
        position_train = np.asarray(rna_fit["oof"], dtype=np.float64)
        rna_position_refit[test] = rna_fit["test"]

        # Step 2: chromatin -> RNA severity position.
        h3_fit = tune(h3_cache, h3_rep, position_train, POSITION_CLIP, base_seed + 2)
        h3_to_position[test] = h3_fit["test"]

        # Floor (i): chromatin -> binary disease call, same folds, same pipeline.
        bin_fit = tune(h3_cache, h3_rep, binary_train, BINARY_CLIP, base_seed + 3)
        h3_to_binary[test] = bin_fit["test"]

        # Floor (ii): position residualised on the binary call (training group means).
        means_train = {}
        for value in (0.0, 1.0):
            sel = binary_train == value
            means_train[value] = float(position_train[sel].mean()) if sel.any() else 0.0
        residual_train = position_train - np.asarray([means_train[v] for v in binary_train])
        res_fit = tune(h3_cache, h3_rep, residual_train, RESIDUAL_CLIP, base_seed + 4)
        h3_to_residual[test] = res_fit["test"]
        eval_means = {}
        for value in (0.0, 1.0):
            sel = binary[train] == value
            eval_means[value] = float(rna_position[train][sel].mean()) if sel.any() else 0.0
        residual_target_eval[test] = rna_position[test] - np.asarray(
            [eval_means[v] for v in binary[test]]
        )

        # Floor (iii): H3K27ac library sum alone, sign fixed on training rows.
        lib = np.log(h3_library)
        design = np.column_stack([np.ones(len(train)), lib[train]])
        coef, *_ = np.linalg.lstsq(design, position_train, rcond=None)
        h3_libsum_fit[test] = np.clip(
            coef[0] + coef[1] * lib[test], POSITION_CLIP[0], POSITION_CLIP[1]
        )

        # Floor (iv): size-matched random-region sets at the same representation.
        train_matrix = h3[train]
        minimum_positive = max(3, int(np.ceil(0.10 * len(train))))
        eligible = np.flatnonzero(np.sum(train_matrix > 0, axis=0) >= minimum_positive)
        real_rep = h3_cache.get("outer_final", train, test, H3_FEATURE_REQUEST)
        matched_size = int(real_rep["actual_feature_count"])
        rng = np.random.default_rng(FIT_SEED + 7_000 + k)
        scramble_rng = np.random.default_rng(FIT_SEED + 6_000 + k)
        sub_rep = {"feature_request": "all_eligible", "pca_components": PCA_COMPONENTS}

        def random_region_fit(size, draw, cache_seed, target, seed):
            cols = np.sort(rng.choice(eligible, size=min(size, len(eligible)), replace=False))
            sub_cache = BASE.RepresentationCache(
                h3[:, cols], [h3_ids[c] for c in cols], ["all_eligible"],
                [PCA_COMPONENTS], cache_seed,
            )
            return tune(sub_cache, sub_rep, target, POSITION_CLIP, seed)["test"]

        for draw in range(RANDOM_REGION_DRAWS):
            random_region_pred[draw, test] = random_region_fit(
                matched_size, draw, 3000 + k * 100 + draw, position_train,
                base_seed + 900_000 + draw,
            )
        for size in RANDOM_REGION_LADDER:
            for draw in range(RANDOM_LADDER_DRAWS):
                random_ladder_pred[size][draw, test] = random_region_fit(
                    size, draw, 10_000 + k * 1000 + RANDOM_REGION_LADDER.index(size) * 100 + draw,
                    position_train, base_seed + 700_000 + size + draw,
                )
        # The random-region arm must itself be able to fail: same size-matched random
        # regions, training target scrambled inside the training rows.
        for draw in range(RANDOM_SCRAMBLED_DRAWS):
            random_scrambled_pred[draw, test] = random_region_fit(
                matched_size, draw, 20_000 + k * 100 + draw,
                position_train[scramble_rng.permutation(len(train))],
                base_seed + 600_000 + draw,
            )

        # Floor (v): permutation null. The training target is scrambled within the
        # training rows only; the evaluation target is untouched. A pipeline that
        # leaks would still score above zero here.
        perm_rng = np.random.default_rng(FIT_SEED + 8_000 + k)
        for rep_index in range(PERMUTATIONS):
            scrambled = position_train[perm_rng.permutation(len(train))]
            fit = tune(h3_cache, h3_rep, scrambled, POSITION_CLIP,
                       base_seed + 800_000 + rep_index)
            permuted_pred[rep_index, test] = fit["test"]

        fold_receipts.append(
            {
                "outer_fold": k,
                "training_participants": int(len(train)),
                "test_participants": int(len(test)),
                "h3_selected_feature_count": matched_size,
                "rna_position_selected": rna_fit["selected"]["parameters"],
                "h3_to_position_selected": h3_fit["selected"]["parameters"],
                "h3_to_binary_selected": bin_fit["selected"]["parameters"],
                "h3_to_residual_selected": res_fit["selected"]["parameters"],
                "h3_library_sum_slope": float(coef[1]),
            }
        )

    for name, vector in (
        ("rna_position_refit", rna_position_refit),
        ("h3_to_position", h3_to_position),
        ("h3_to_binary", h3_to_binary),
        ("h3_to_residual", h3_to_residual),
        ("h3_libsum_fit", h3_libsum_fit),
        ("residual_target_eval", residual_target_eval),
    ):
        if not np.all(np.isfinite(vector)):
            raise SystemExit(f"{name} has non-finite entries")

    # ---- arms ------------------------------------------------------------------------
    arms = {
        "A1_frozen_h3_position_vs_rna_position": (h3_position_frozen, rna_position),
        "A2_refit_h3_to_position_vs_rna_position": (h3_to_position, rna_position),
        "F1a_frozen_h3_disease_score_vs_rna_position": (h3_disease_score_frozen, rna_position),
        "F1b_refit_h3_binary_vs_rna_position": (h3_to_binary, rna_position),
        "F2_refit_h3_residual_vs_residual_position": (h3_to_residual, residual_target_eval),
        "F3_h3_library_sum_vs_rna_position": (h3_libsum_fit, rna_position),
        "REF_binary_label_vs_rna_position": (binary, rna_position),
        "REF_rna_position_vs_true_nas": (rna_position, nas_true),
        "REF_frozen_h3_position_vs_true_nas": (h3_position_frozen, nas_true),
        "REF_refit_h3_position_vs_true_nas": (h3_to_position, nas_true),
        "REF_rna_position_refit_vs_frozen": (rna_position_refit, rna_position),
        "REF_frozen_h3_position_vs_binary": (h3_position_frozen, binary),
        "REF_rna_position_vs_binary": (rna_position, binary),
        "REF_rna_position_vs_fibrosis": (rna_position, fibrosis),
        "REF_frozen_h3_position_vs_fibrosis": (h3_position_frozen, fibrosis),
    }
    point = {name: spearman(a, b) for name, (a, b) in arms.items()}

    replicate = {}
    for name, (a, b) in arms.items():
        values = np.empty(BOOTSTRAP_REPLICATES)
        for i in range(BOOTSTRAP_REPLICATES):
            idx = boot[i]
            values[i] = spearman(a[idx], b[idx])
        replicate[name] = values

    def interval(values: np.ndarray) -> dict:
        valid = values[np.isfinite(values)]
        return {
            "valid_replicates": int(valid.size),
            "percentile_2p5": float(np.percentile(valid, 2.5)),
            "percentile_97p5": float(np.percentile(valid, 97.5)),
        }

    contrasts = {
        "position_minus_frozen_disease_score": (
            "A1_frozen_h3_position_vs_rna_position",
            "F1a_frozen_h3_disease_score_vs_rna_position",
        ),
        "refit_position_minus_refit_binary": (
            "A2_refit_h3_to_position_vs_rna_position",
            "F1b_refit_h3_binary_vs_rna_position",
        ),
        "position_minus_library_sum": (
            "A2_refit_h3_to_position_vs_rna_position",
            "F3_h3_library_sum_vs_rna_position",
        ),
        "position_minus_binary_label": (
            "A1_frozen_h3_position_vs_rna_position",
            "REF_binary_label_vs_rna_position",
        ),
    }
    contrast_report = {}
    for name, (left, right) in contrasts.items():
        diff = replicate[left] - replicate[right]
        valid = diff[np.isfinite(diff)]
        contrast_report[name] = {
            "arm_a": left,
            "arm_b": right,
            "point_estimate": point[left] - point[right],
            "percentile_2p5": float(np.percentile(valid, 2.5)),
            "percentile_97p5": float(np.percentile(valid, 97.5)),
            "probability_arm_a_greater": float(np.mean(valid > 0)),
            "valid_replicates": int(valid.size),
        }

    # ---- nulls -----------------------------------------------------------------------
    random_rho = np.asarray(
        [spearman(random_region_pred[d], rna_position) for d in range(RANDOM_REGION_DRAWS)]
    )
    ladder_rho = {
        size: np.asarray(
            [spearman(random_ladder_pred[size][d], rna_position) for d in range(RANDOM_LADDER_DRAWS)]
        )
        for size in RANDOM_REGION_LADDER
    }
    random_scrambled_rho = np.asarray(
        [spearman(random_scrambled_pred[d], rna_position) for d in range(RANDOM_SCRAMBLED_DRAWS)]
    )
    perm_rho = np.asarray(
        [spearman(permuted_pred[p], rna_position) for p in range(PERMUTATIONS)]
    )
    observed = point["A2_refit_h3_to_position_vs_rna_position"]
    nulls = {
        "size_matched_random_regions": {
            "draws": RANDOM_REGION_DRAWS,
            "mean": float(np.mean(random_rho)),
            "median": float(np.median(random_rho)),
            "percentile_5": float(np.percentile(random_rho, 5)),
            "percentile_95": float(np.percentile(random_rho, 95)),
            "maximum": float(np.max(random_rho)),
            "observed_refit_rho": observed,
            "draws_at_or_above_observed": int(np.sum(random_rho >= observed)),
        },
        "size_matched_random_regions_scrambled_target": {
            "draws": RANDOM_SCRAMBLED_DRAWS,
            "mean": float(np.mean(random_scrambled_rho)),
            "median": float(np.median(random_scrambled_rho)),
            "minimum": float(np.min(random_scrambled_rho)),
            "maximum": float(np.max(random_scrambled_rho)),
            "note": (
                "Same size-matched random-region representation, training target scrambled. "
                "This is the evidence that the random-region arm is not a pipeline artifact."
            ),
        },
        "random_region_size_ladder": {
            str(size): {
                "draws": RANDOM_LADDER_DRAWS,
                "mean": float(np.mean(ladder_rho[size])),
                "median": float(np.median(ladder_rho[size])),
                "percentile_5": float(np.percentile(ladder_rho[size], 5)),
                "percentile_95": float(np.percentile(ladder_rho[size], 95)),
            }
            for size in RANDOM_REGION_LADDER
        },
        "scrambled_training_target": {
            "replicates": PERMUTATIONS,
            "mean": float(np.mean(perm_rho)),
            "median": float(np.median(perm_rho)),
            "percentile_2p5": float(np.percentile(perm_rho, 2.5)),
            "percentile_97p5": float(np.percentile(perm_rho, 97.5)),
            "maximum": float(np.max(perm_rho)),
            "replicates_at_or_above_observed": int(np.sum(perm_rho >= observed)),
            "null_validity_note": (
                "This null re-runs the identical fold machinery with the training target "
                "scrambled. A leaking pipeline would still score above zero here, so a "
                "distribution centred on zero is the evidence that the machinery can "
                "reject a target it must fail."
            ),
        },
    }

    per_fold = []
    for k in range(5):
        sel = outer_fold == k
        per_fold.append(
            {
                "outer_fold": k,
                "n": int(sel.sum()),
                "A1_frozen": spearman(h3_position_frozen[sel], rna_position[sel]),
                "A2_refit": spearman(h3_to_position[sel], rna_position[sel]),
                "F2_residual": spearman(h3_to_residual[sel], residual_target_eval[sel]),
            }
        )

    # ---- deposits --------------------------------------------------------------------
    with (output / "per_sample_vectors.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "participant_id", "outer_fold", "stage3", "binary_disease_call",
                "nash_crn_component_sum", "fibrosis", "rna_position_frozen",
                "rna_position_refit", "h3_position_frozen", "h3_to_position_refit",
                "h3_to_binary_refit", "h3_to_residual_refit", "h3_libsum_fit",
                "residual_target_eval", "h3_library_sum", "h3_disease_score_frozen",
                "rna_probability_NOR_frozen",
            ]
        )
        for i, pid in enumerate(participants):
            writer.writerow(
                [
                    pid, int(outer_fold[i]), stage3[i], int(binary[i]),
                    int(nas_true[i]), int(fibrosis[i]),
                    format(rna_position[i], ".17g"), format(rna_position_refit[i], ".17g"),
                    format(h3_position_frozen[i], ".17g"), format(h3_to_position[i], ".17g"),
                    format(h3_to_binary[i], ".17g"), format(h3_to_residual[i], ".17g"),
                    format(h3_libsum_fit[i], ".17g"), format(residual_target_eval[i], ".17g"),
                    format(h3_library[i], ".17g"), format(h3_disease_score_frozen[i], ".17g"),
                    format(rna_stage_nor[i], ".17g"),
                ]
            )

    np.savez_compressed(
        output / "bootstrap_replicates.npz",
        bootstrap_indices_sha256=np.asarray([boot_digest]),
        **{f"rho__{name}": values for name, values in replicate.items()},
        random_region_rho=random_rho,
        random_scrambled_rho=random_scrambled_rho,
        permutation_rho=perm_rho,
        **{f"ladder_rho__{size}": ladder_rho[size] for size in RANDOM_REGION_LADDER},
        random_region_predictions=random_region_pred,
        permutation_predictions=permuted_pred,
    )

    report = {
        "schema_version": "masld-bench-gse267145-chromatin-vs-rna-severity-ordering-v1",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "unit": "sample",
        "samples": n,
        "donor_key_exists": False,
        "claim_grade": "development_only_project_exposed_downstream_demo",
        "champion_eligible": False,
        "external_claim_eligible": False,
        "target_definition": {
            "rna_severity_position": (
                "frozen five-seed out-of-fold prediction of nash_crn_component_sum by "
                "rna_hvg_pca_elastic_net in model-scoring-078, in the frozen "
                "participant-held-out folds"
            ),
            "why_not_fibrosis": "97.1% binary in this cohort; excluded by task constraint",
            "why_not_stage3": "nearly the binary disease call that floor (i) must separate from",
            "chromatin_never_entered_any_rna_derived_axis": True,
            "training_target_for_the_chromatin_refit": (
                "fold-specific inner out-of-fold RNA position built from fold-k training "
                "labels only, so no held-out label reaches the chromatin training target"
            ),
        },
        "configuration": {
            "rna_feature_request": RNA_FEATURE_REQUEST,
            "h3_feature_request": H3_FEATURE_REQUEST,
            "pca_components": PCA_COMPONENTS,
            "alphas": ALPHAS,
            "l1_ratios": L1_RATIOS,
            "inner_selection": "four-fold inner CV with the frozen one-standard-error rule",
            "folds_recomputed": False,
            "fit_seed": FIT_SEED,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "bootstrap_replicates": BOOTSTRAP_REPLICATES,
            "bootstrap_indices_sha256": boot_digest,
            "bootstrap_distinct_units_mean": float(distinct.mean()),
            "bootstrap_distinct_units_fraction": float(distinct.mean() / n),
        },
        "provenance": provenance,
        "point_estimates": point,
        "intervals": {name: interval(values) for name, values in replicate.items()},
        "contrasts": contrast_report,
        "nulls": nulls,
        "per_fold": per_fold,
        "fold_receipts": fold_receipts,
        "reporting_rules": [
            "No p-values, no ranking, no champion or best-model language.",
            "An interval wholly on one side of zero is reported as such, never as significant.",
            "Cross-sectional histology-associated remodeling only.",
        ],
    }
    (output / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (output / "session_info.json").write_text(
        json.dumps(
            {
                "python": sys.version,
                "platform": platform.platform(),
                "numpy": np.__version__,
                "node": platform.node(),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(json.dumps({"output": str(output), "point_estimates": point}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

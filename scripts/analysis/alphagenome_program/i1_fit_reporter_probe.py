#!/usr/bin/env python3
"""Adaptation arm 2, stage 2: the published reporter difference head on frozen AlphaGenome embeddings.

Prespecified in
`GWAS/finemapping/results/alphagenome_program/ad-arm2-mpra-20260915T114407Z/PRESPECIFICATION.md`
(sha256 f44a912f73a1a30cddcfb60aacd5ddd0edb2e49838a0b8c6011c28be3ad267fc) plus addendum 01
(sha256 260914ce33d1aed30dd68fb856df64bf051979f84d34d1e07004824a2b5284bb).

Nested increment: the head, the folds, the blocks, the seeds, the projection width, the grids, the
selection rule and the metric are the published DNA-language arm's. The ONLY component that changes is the
backbone that produced the raw allele embeddings.

Every helper that defines the recipe is imported from the benchmark tree rather than retyped, so a
divergence is impossible:
  scripts.fit_gse281364_enformer_sei_heads._bootstrap_indices, ._training_transform
  scripts.fit_gse281364_open_sequence_heads.fit_inner_projection
  scripts.gse281364_dna_lm_native_contract.fit_projection, .apply_projection, .head_features
  scripts.evaluate_gse281364_dna_lm_common_lane.load_outcomes
Nothing in the benchmark tree is written.

Hard gate before any AlphaGenome number is produced: this script must reproduce the published
hyenadna/delta_ridge score 0.21622332457761556 and the allele-identity control 0.13421551707194201 from the
frozen campaign predictions using its own loader and its own metric. If it cannot, it exits non-zero.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import pathlib
import sys
import time

import numpy as np

PROJ = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BENCH = PROJ / "Analysis/MASLD_Model_Benchmark"
ROW_UNIVERSE = BENCH / "executions/model-check-220-21088696/contract/row_universe.tsv"
OUTCOMES = BENCH / "executions/gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz"
PUBLISHED_FIT = BENCH / "executions/model-cpu-train-605-21099008/fit/oof_predictions.tsv.gz"
PUBLISHED_EVAL = BENCH / "executions/model-cpu-train-605-21099008/evaluation/model_head_summary.tsv"

CONTEXTS = ("HepG2_control", "HepG2_PAOA")
SEEDS = (1103, 2909, 4721, 6673, 8111)
INNER_SALT = 10_000
OUTER_SALT = 20_000
PROJECTION_WIDTH = 256
DELTA_SLICE = slice(512, 768)
DELTA_TOP_K = (32, 128, 256)
ALPHA_GRID = (1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4)
BOOTSTRAP_RESAMPLES = 10_000
BOOTSTRAP_SEED_PRIMARY = 20260914   # the adaptation arm's own bootstrap seed
BOOTSTRAP_SEED_PUBLISHED = 20260825  # the published campaign's, reported beside it
MIN_USABLE_DRAWS = 9_500

COMPARATOR = ("hyenadna", "delta_ridge")
CONTROL = ("available_simple_controls", "allele_identity_ridge")
PUBLISHED_COMPARATOR_SCORE = 0.21622332457761556
PUBLISHED_CONTROL_SCORE = 0.13421551707194201

PRESPEC_SHA = "f44a912f73a1a30cddcfb60aacd5ddd0edb2e49838a0b8c6011c28be3ad267fc"
ADDENDUM_SHA = "260914ce33d1aed30dd68fb856df64bf051979f84d34d1e07004824a2b5284bb"

CELLS = (
    # (cell_id, input_length_bp, pooling_array, pooling_label, confirmatory)
    ("len2048_poolregion", 2048, "pooled_region", "pool_region_mean_12x128bp", True),
    ("len2048_centrebin", 2048, "pooled_centre_bin", "single_variant_bin_128bp", False),
    ("len16384_poolregion", 16384, "pooled_region", "pool_region_mean_12x128bp", False),
    ("len16384_centrebin", 16384, "pooled_centre_bin", "single_variant_bin_128bp", False),
)


class ProbeFitError(RuntimeError):
    pass


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}] {msg}", flush=True)


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------------------------------
# metric, verbatim from evaluate_gse281364_dna_language_seeded_heads.py
# --------------------------------------------------------------------------------------------------
def macro(values) -> float | None:
    if any(v is None or not math.isfinite(v) for v in values):
        return None
    clipped = np.clip(np.asarray(values), -0.999999, 0.999999)
    return float(np.tanh(np.mean(np.arctanh(clipped))))


def fast_spearman(observed: np.ndarray, predicted: np.ndarray) -> float | None:
    from scipy.stats import rankdata

    if np.all(observed == observed[0]) or np.all(predicted == predicted[0]):
        return None
    left = rankdata(observed)
    right = rankdata(predicted)
    left -= left.mean()
    right -= right.mean()
    denominator = float(np.sqrt(np.sum(left**2) * np.sum(right**2)))
    return None if denominator == 0 else float(np.sum(left * right) / denominator)


def point_macro(observed, predictions, keep=None) -> float | None:
    from scipy.stats import spearmanr

    per = []
    for context in CONTEXTS:
        o = observed[context] if keep is None else observed[context][keep]
        p = predictions[context] if keep is None else predictions[context][keep]
        per.append(float(spearmanr(o, p).statistic))
    return macro(per), per


# --------------------------------------------------------------------------------------------------
# substrate
# --------------------------------------------------------------------------------------------------
def load_split():
    """Row-universe insertion order, exactly as the published fit preserves it."""
    metadata: dict[str, tuple[int, str]] = {}
    with ROW_UNIVERSE.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        value = (int(r["outer_fold"].removeprefix("fold-")), r["long_range_block_id"])
        if metadata.setdefault(r["element_id"], value) != value:
            raise ProbeFitError("row universe fold/block is not unique per element")
    elements = list(metadata)
    if len(rows) != 10_330 or len(elements) != 1033:
        raise ProbeFitError("row universe census differs")
    folds = np.asarray([metadata[e][0] for e in elements], dtype=np.int64)
    blocks = np.asarray([metadata[e][1] for e in elements])
    if set(folds.tolist()) != set(range(5)) or len(set(blocks.tolist())) != 239:
        raise ProbeFitError("fold or block census differs")
    return elements, folds, blocks


def load_frozen_oof(elements, wanted):
    index = {e: i for i, e in enumerate(elements)}
    out = {(m, h, s, c): np.full(len(elements), np.nan) for m, h in wanted for s in SEEDS for c in CONTEXTS}
    with gzip.open(PUBLISHED_FIT, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["model_id"], row["head_id"])
            if key not in wanted:
                continue
            slot = out[(*key, int(row["seed"]), row["assay_context_id"])]
            pos = index[row["element_id"]]
            if math.isfinite(slot[pos]):
                raise ProbeFitError("duplicate frozen prediction")
            slot[pos] = float(row["prediction"])
    for k, v in out.items():
        if not np.isfinite(v).all():
            raise ProbeFitError(f"frozen prediction coverage differs for {k}")
    return out


def ensemble(oof, key):
    return {c: np.mean(np.vstack([oof[(*key, s, c)] for s in SEEDS]), axis=0) for c in CONTEXTS}


# --------------------------------------------------------------------------------------------------
# the head
# --------------------------------------------------------------------------------------------------
def fit_cell(raw, outcomes_by_context, folds, blocks, helpers):
    """Five seeds x two contexts x five outer folds of the published dual-projection delta ridge."""
    fit_inner_projection = helpers["fit_inner_projection"]
    fit_projection = helpers["fit_projection"]
    apply_projection = helpers["apply_projection"]
    head_features = helpers["head_features"]
    bootstrap_indices = helpers["_bootstrap_indices"]
    training_transform = helpers["_training_transform"]
    from sklearn.linear_model import Ridge

    n = raw.shape[0]
    inner_feat: dict[int, np.ndarray] = {}
    outer_feat: dict[int, np.ndarray] = {}
    projection_rows = []
    for held in range(5):
        inner_fold = (held + 1) % 5
        inner_mask = (folds != held) & (folds != inner_fold)
        p_inner = fit_inner_projection(raw, inner_mask)
        inner_feat[held] = head_features(apply_projection(raw, p_inner)).astype(np.float64)
        p_outer = fit_projection(raw, folds, held, width=PROJECTION_WIDTH)
        outer_feat[held] = head_features(apply_projection(raw, p_outer)).astype(np.float64)
        if inner_feat[held].shape != (n, 4 * PROJECTION_WIDTH):
            raise ProbeFitError("inner head-feature width differs")
        projection_rows.append(
            {
                "held_out_fold": held,
                "inner_validation_fold": inner_fold,
                "inner_projection_training_elements": int(inner_mask.sum()),
                "outer_projection_training_elements": int((folds != held).sum()),
                "inner_projection_fit_on_inner_training_only": "true",
                "outer_projection_fit_on_outer_training_only": "true",
                "held_outcomes_used": "false",
            }
        )

    predictions: dict[tuple[int, str], np.ndarray] = {}
    selections = []
    for context in CONTEXTS:
        y = outcomes_by_context[context]
        for seed in SEEDS:
            oof = np.full(n, np.nan)
            for held in range(5):
                inner_fold = (held + 1) % 5
                outer_test = folds == held
                outer_train = ~outer_test
                inner_valid = folds == inner_fold
                inner_train = outer_train & ~inner_valid
                if set(blocks[outer_test]) & set(blocks[outer_train]):
                    raise ProbeFitError("held long-range block overlaps training")
                xi = inner_feat[held][:, DELTA_SLICE]
                xo = outer_feat[held][:, DELTA_SLICE]

                inner_boot = bootstrap_indices(inner_train, blocks, seed=seed, salt=INNER_SALT + held)
                mean, scale, rank = training_transform(xi, y, inner_boot)
                choices = []
                for k in sorted({min(int(v), int(rank.size)) for v in DELTA_TOP_K if int(v) > 0}):
                    if k == 0:
                        continue
                    sel = rank[:k]
                    xt = (xi[inner_boot][:, sel] - mean[sel]) / scale[sel]
                    xv = (xi[inner_valid][:, sel] - mean[sel]) / scale[sel]
                    for alpha in ALPHA_GRID:
                        m = Ridge(alpha=float(alpha), fit_intercept=True, solver="lsqr", tol=1.0e-7)
                        m.fit(xt, y[inner_boot])
                        pred = m.predict(xv)
                        choices.append(
                            (float(np.sqrt(np.mean((y[inner_valid] - pred) ** 2))), k, float(alpha))
                        )
                constant_rmse = float(np.sqrt(np.mean((y[inner_valid] - y[inner_boot].mean()) ** 2)))
                best_rmse, best_k, best_alpha = min(
                    [*choices, (constant_rmse, 0, math.inf)], key=lambda v: (v[0], v[1], v[2])
                )
                outer_boot = bootstrap_indices(outer_train, blocks, seed=seed, salt=OUTER_SALT + held)
                o_mean, o_scale, o_rank = training_transform(xo, y, outer_boot)
                if best_k == 0 or o_rank.size == 0:
                    sel = np.empty(0, dtype=np.int64)
                    intercept = float(y[outer_boot].mean())
                    pred = np.full(int(outer_test.sum()), intercept)
                    selected_alpha = "constant_training_mean"
                else:
                    sel = o_rank[: min(best_k, int(o_rank.size))]
                    xt = (xo[outer_boot][:, sel] - o_mean[sel]) / o_scale[sel]
                    xte = (xo[outer_test][:, sel] - o_mean[sel]) / o_scale[sel]
                    m = Ridge(alpha=best_alpha, fit_intercept=True, solver="lsqr", tol=1.0e-7)
                    m.fit(xt, y[outer_boot])
                    pred = np.asarray(m.predict(xte), dtype=np.float64)
                    selected_alpha = best_alpha
                oof[outer_test] = pred
                selections.append(
                    {
                        "assay_context_id": context,
                        "seed": seed,
                        "held_out_fold": held,
                        "inner_validation_fold": inner_fold,
                        "inner_training_elements": int(inner_train.sum()),
                        "inner_validation_elements": int(inner_valid.sum()),
                        "outer_training_elements": int(outer_train.sum()),
                        "outer_test_elements": int(outer_test.sum()),
                        "outer_test_blocks": len(set(blocks[outer_test].tolist())),
                        "inner_bootstrap_rows": int(inner_boot.size),
                        "outer_bootstrap_rows": int(outer_boot.size),
                        "selected_feature_count": int(sel.size),
                        "selected_alpha": selected_alpha,
                        "selected_inner_rmse": best_rmse,
                        "constant_inner_rmse": constant_rmse,
                        "held_outcomes_used_for_fit_or_tuning": False,
                        "held_features_used_for_preprocessing_selection_or_tuning": False,
                        "held_blocks_overlap_training": False,
                        "posthoc_calibration": False,
                    }
                )
            if not np.isfinite(oof).all():
                raise ProbeFitError("OOF prediction coverage differs")
            predictions[(seed, context)] = oof
    return predictions, selections, projection_rows


# --------------------------------------------------------------------------------------------------
# bootstrap
# --------------------------------------------------------------------------------------------------
def paired_bootstrap(observed, candidates, blocks, *, resamples, seed):
    """One block draw per iteration, shared by every candidate, so differences are paired."""
    unique = sorted(set(blocks.tolist()))
    by_block = {b: np.flatnonzero(blocks == b) for b in unique}
    out = {name: np.full(resamples, np.nan) for name in candidates}
    rng = np.random.default_rng(seed)
    for it in range(resamples):
        sampled = rng.integers(0, len(unique), size=len(unique))
        idx = np.concatenate([by_block[unique[i]] for i in sampled])
        for name, pred in candidates.items():
            v = macro([fast_spearman(observed[c][idx], pred[c][idx]) for c in CONTEXTS])
            if v is not None:
                out[name][it] = v
    return out


def interval(samples):
    valid = np.isfinite(samples)
    if int(valid.sum()) < MIN_USABLE_DRAWS:
        return None, None, None, int(valid.sum())
    lo, hi = (float(v) for v in np.quantile(samples[valid], (0.025, 0.975)))
    return lo, hi, float(np.std(samples[valid], ddof=1)), int(valid.sum())


def two_sided_p(diff_samples):
    valid = np.isfinite(diff_samples)
    if int(valid.sum()) < MIN_USABLE_DRAWS:
        return None
    d = diff_samples[valid]
    n = d.size
    p_le = (np.sum(d <= 0) + 1) / (n + 1)
    p_ge = (np.sum(d >= 0) + 1) / (n + 1)
    return float(min(1.0, 2.0 * min(p_le, p_ge)))


def bh(pvalues):
    items = [(i, p) for i, p in enumerate(pvalues) if p is not None]
    out = [None] * len(pvalues)
    if not items:
        return out
    items.sort(key=lambda v: v[1])
    m = len(items)
    prev = 1.0
    for rank in range(m, 0, -1):
        i, p = items[rank - 1]
        prev = min(prev, p * m / rank)
        out[i] = float(prev)
    return out


def write_tsv(path, rows):
    fields: list[str] = []
    for r in rows:
        for f in r:
            if f not in fields:
                fields.append(f)
    with pathlib.Path(path).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({f: r.get(f, "") for f in fields})
    log(f"wrote {path}")


def fmt(v):
    return "" if v is None else format(v, ".17g")


# --------------------------------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--embeddings-2048", required=True)
    ap.add_argument("--embeddings-16384", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--card", required=True)
    args = ap.parse_args()

    sys.path.insert(0, str(BENCH))
    sys.path.insert(0, str(BENCH / "src"))
    from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes
    from scripts.fit_gse281364_enformer_sei_heads import _bootstrap_indices, _training_transform
    from scripts.fit_gse281364_open_sequence_heads import fit_inner_projection
    from scripts.gse281364_dna_lm_native_contract import apply_projection, fit_projection, head_features

    helpers = {
        "fit_inner_projection": fit_inner_projection,
        "fit_projection": fit_projection,
        "apply_projection": apply_projection,
        "head_features": head_features,
        "_bootstrap_indices": _bootstrap_indices,
        "_training_transform": _training_transform,
    }

    out = pathlib.Path(args.out)
    (out / "tables").mkdir(parents=True, exist_ok=True)

    elements, folds, blocks = load_split()
    log(f"split: {len(elements)} elements, {len(set(blocks.tolist()))} blocks, folds {np.bincount(folds)}")
    oc = load_outcomes(OUTCOMES, set(elements))
    observed = {c: np.asarray([oc[(e, c)]["mean"] for e in elements]) for c in CONTEXTS}

    # ---------------------------------------------------------------- hard reproduction gate
    frozen = load_frozen_oof(elements, {COMPARATOR, CONTROL})
    comp_ens = ensemble(frozen, COMPARATOR)
    ctrl_ens = ensemble(frozen, CONTROL)
    comp_point, comp_per = point_macro(observed, comp_ens)
    ctrl_point, ctrl_per = point_macro(observed, ctrl_ens)
    gate = {
        "comparator": "/".join(COMPARATOR),
        "comparator_reproduced": comp_point,
        "comparator_published": PUBLISHED_COMPARATOR_SCORE,
        "comparator_per_context": comp_per,
        "comparator_exact": comp_point == PUBLISHED_COMPARATOR_SCORE,
        "control": "/".join(CONTROL),
        "control_reproduced": ctrl_point,
        "control_published": PUBLISHED_CONTROL_SCORE,
        "control_per_context": ctrl_per,
        "control_exact": ctrl_point == PUBLISHED_CONTROL_SCORE,
        "published_summary_sha256": sha256_file(PUBLISHED_EVAL),
        "published_oof_sha256": sha256_file(PUBLISHED_FIT),
    }
    (out / "reproduction_gate.json").write_text(json.dumps(gate, indent=1) + "\n")
    if not (gate["comparator_exact"] and gate["control_exact"]):
        log(f"GATE FAILED: {gate}")
        return 2
    log(f"gate passed: comparator {comp_point!r}, control {ctrl_point!r} both digit-for-digit")

    # ---------------------------------------------------------------- embeddings
    npz = {}
    for length, path in ((2048, args.embeddings_2048), (16384, args.embeddings_16384)):
        d = np.load(path, allow_pickle=False)
        if int(d["input_length_bp"]) != length:
            raise ProbeFitError(f"{path} is not the {length} bp extraction")
        card = str(d["gpu_card"])
        if card != args.card or card != "l40s":
            raise ProbeFitError(f"{path} was produced on card {card!r}, expected l40s")
        ids = d["element_ids"].astype(str).tolist()
        if ids != elements:
            # reorder to the row-universe insertion order
            pos = {e: i for i, e in enumerate(ids)}
            if set(ids) != set(elements):
                raise ProbeFitError(f"{path} covers a different element set")
            order = np.asarray([pos[e] for e in elements])
        else:
            order = np.arange(len(elements))
        npz[length] = {
            "pooled_region": d["pooled_region"][order],
            "pooled_centre_bin": d["pooled_centre_bin"][order],
            "extracted": d["extracted"][order],
            "pool_bin_lo": int(d["pool_bin_lo"]),
            "pool_bin_hi": int(d["pool_bin_hi"]),
            "card": card,
            "path": path,
            "sha256": sha256_file(path),
        }
        log(
            f"{length} bp: {int(npz[length]['extracted'].sum())} of {len(elements)} elements extracted, "
            f"pool bins [{npz[length]['pool_bin_lo']},{npz[length]['pool_bin_hi']}), card {card}"
        )

    rejected_elements = sorted(
        {elements[i] for length in npz for i in np.flatnonzero(~npz[length]["extracted"])}
    )
    matched_keep = np.ones(len(elements), bool)
    for length in npz:
        matched_keep &= npz[length]["extracted"]
    log(f"rejected at some length: {rejected_elements}; matched set is {int(matched_keep.sum())} elements")

    # ---------------------------------------------------------------- fit every cell
    cell_pred: dict[str, dict[str, np.ndarray]] = {}
    cell_seed_pred: dict[str, dict[tuple[int, str], np.ndarray]] = {}
    cell_keep: dict[str, np.ndarray] = {}
    per_seed_rows = []
    selection_rows = []
    projection_rows = []
    oof_rows = []
    for cell_id, length, array_name, pooling_label, confirmatory in CELLS:
        keep = npz[length]["extracted"]
        raw = npz[length][array_name][keep].astype(np.float64)
        f = folds[keep]
        b = blocks[keep]
        y = {c: observed[c][keep] for c in CONTEXTS}
        log(
            f"=== cell {cell_id}: n={raw.shape[0]} elements, {len(set(b.tolist()))} blocks, "
            f"raw width {raw.shape[2]}"
        )
        t0 = time.time()
        preds, sels, projs = fit_cell(raw, y, f, b, helpers)
        log(f"    fitted in {time.time() - t0:.1f}s")
        cell_keep[cell_id] = keep
        cell_seed_pred[cell_id] = preds
        cell_pred[cell_id] = {
            c: np.mean(np.vstack([preds[(s, c)] for s in SEEDS]), axis=0) for c in CONTEXTS
        }
        for s in sels:
            selection_rows.append({"cell_id": cell_id, "input_length_bp": length, "pooling": pooling_label, "gpu_card": args.card, **s})
        for p in projs:
            projection_rows.append({"cell_id": cell_id, "input_length_bp": length, **p})
        for seed in SEEDS:
            per = [float(fast_spearman(y[c], preds[(seed, c)])) for c in CONTEXTS]
            per_seed_rows.append(
                {
                    "cell_id": cell_id,
                    "input_length_bp": length,
                    "pooling": pooling_label,
                    "seed": seed,
                    "gpu_card": args.card,
                    "n_elements": raw.shape[0],
                    "spearman_HepG2_control": fmt(per[0]),
                    "spearman_HepG2_PAOA": fmt(per[1]),
                    "macro": fmt(macro(per)),
                }
            )
        kept_elements = [e for e, k in zip(elements, keep) if k]
        for seed in SEEDS:
            for c in CONTEXTS:
                for e, v in zip(kept_elements, preds[(seed, c)]):
                    oof_rows.append(
                        {
                            "cell_id": cell_id,
                            "model_id": "alphagenome_local_frozen_trunk",
                            "head_id": "delta_ridge",
                            "input_length_bp": length,
                            "pooling": pooling_label,
                            "gpu_card": args.card,
                            "seed": seed,
                            "assay_context_id": c,
                            "element_id": e,
                            "prediction": format(float(v), ".17g"),
                            "experimental_replicates": 4,
                            "biological_donors": 0,
                            "outcome_role": "exposed_development_MPRA_only",
                        }
                    )

    # ---------------------------------------------------------------- score every cell and contrast it
    score_rows = []
    contrast_rows = []
    scopes = [("own", None)]
    if int(matched_keep.sum()) < len(elements):
        scopes.append((f"matched_{int(matched_keep.sum())}", matched_keep))
    for bseed in (BOOTSTRAP_SEED_PRIMARY, BOOTSTRAP_SEED_PUBLISHED):
        for scope, keep_mask in scopes:
            for cell_id, length, _a, pooling_label, confirmatory in CELLS:
                keep = cell_keep[cell_id] if keep_mask is None else (cell_keep[cell_id] & keep_mask)
                sub_cell = {c: cell_pred[cell_id][c] for c in CONTEXTS}
                # cell predictions live on the cell's own kept rows; expand back to full length
                full = {c: np.full(len(elements), np.nan) for c in CONTEXTS}
                for c in CONTEXTS:
                    full[c][cell_keep[cell_id]] = sub_cell[c]
                idx = np.flatnonzero(keep)
                cands = {
                    "alphagenome": {c: full[c][idx] for c in CONTEXTS},
                    "hyenadna": {c: comp_ens[c][idx] for c in CONTEXTS},
                    "control": {c: ctrl_ens[c][idx] for c in CONTEXTS},
                }
                obs = {c: observed[c][idx] for c in CONTEXTS}
                bl = blocks[idx]
                points = {}
                per_ctx = {}
                for name, pred in cands.items():
                    p, per = point_macro(obs, pred)
                    points[name] = p
                    per_ctx[name] = per
                samples = paired_bootstrap(
                    obs, cands, bl, resamples=BOOTSTRAP_RESAMPLES, seed=bseed
                )
                lo, hi, se, usable = interval(samples["alphagenome"])
                score_rows.append(
                    {
                        "cell_id": cell_id,
                        "scope": scope,
                        "bootstrap_seed": bseed,
                        "input_length_bp": length,
                        "pooling": pooling_label,
                        "confirmatory": str(bool(confirmatory and scope == "own" and bseed == BOOTSTRAP_SEED_PRIMARY)).lower(),
                        "gpu_card": args.card,
                        "n_elements": int(keep.sum()),
                        "n_blocks": len(set(bl.tolist())),
                        "primary_score": fmt(points["alphagenome"]),
                        "spearman_HepG2_control": fmt(per_ctx["alphagenome"][0]),
                        "spearman_HepG2_PAOA": fmt(per_ctx["alphagenome"][1]),
                        "ci_low": fmt(lo),
                        "ci_high": fmt(hi),
                        "bootstrap_se": fmt(se),
                        "usable_draws": usable,
                        "hyenadna_score_this_scope": fmt(points["hyenadna"]),
                        "control_score_this_scope": fmt(points["control"]),
                        "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
                        "bootstrap_unit": "long_range_block_id",
                    }
                )
                for against in ("hyenadna", "control"):
                    d = samples["alphagenome"] - samples[against]
                    dlo, dhi, dse, dusable = interval(d)
                    diff = (
                        None
                        if points["alphagenome"] is None or points[against] is None
                        else points["alphagenome"] - points[against]
                    )
                    p = two_sided_p(d)
                    contrast_rows.append(
                        {
                            "cell_id": cell_id,
                            "scope": scope,
                            "bootstrap_seed": bseed,
                            "against": against,
                            "confirmatory": str(
                                bool(
                                    confirmatory
                                    and scope == "own"
                                    and bseed == BOOTSTRAP_SEED_PRIMARY
                                    and against == "hyenadna"
                                )
                            ).lower(),
                            "input_length_bp": length,
                            "pooling": pooling_label,
                            "gpu_card": args.card,
                            "n_elements": int(keep.sum()),
                            "n_blocks": len(set(bl.tolist())),
                            "alphagenome_score": fmt(points["alphagenome"]),
                            "reference_score": fmt(points[against]),
                            "paired_difference": fmt(diff),
                            "paired_ci_low": fmt(dlo),
                            "paired_ci_high": fmt(dhi),
                            "paired_se": fmt(dse),
                            "interval_includes_zero": (
                                "" if dlo is None else str(bool(dlo <= 0.0 <= dhi)).lower()
                            ),
                            "two_sided_bootstrap_p": fmt(p),
                            "usable_draws": dusable,
                        }
                    )

    # BH inside the named secondary family: the non-confirmatory configuration cells, primary seed,
    # own scope, against the comparator.
    family = [
        r
        for r in contrast_rows
        if r["confirmatory"] == "false"
        and r["scope"] == "own"
        and int(r["bootstrap_seed"]) == BOOTSTRAP_SEED_PRIMARY
        and r["against"] == "hyenadna"
    ]
    ps = [None if r["two_sided_bootstrap_p"] == "" else float(r["two_sided_bootstrap_p"]) for r in family]
    for r, q in zip(family, bh(ps)):
        r["bh_q_within_ad_arm2_reporter_probe_configuration"] = fmt(q)
    for r in contrast_rows:
        r.setdefault("bh_q_within_ad_arm2_reporter_probe_configuration", "")

    write_tsv(out / "tables/cell_scores.tsv", score_rows)
    write_tsv(out / "tables/paired_contrasts.tsv", contrast_rows)
    write_tsv(out / "tables/per_seed_macro.tsv", per_seed_rows)
    write_tsv(out / "tables/head_selection.tsv", selection_rows)
    write_tsv(out / "tables/projection_audit.tsv", projection_rows)
    with gzip.open(out / "oof_predictions.tsv.gz", "wt", newline="") as fh:
        fields = list(oof_rows[0])
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(oof_rows)
    log(f"wrote {out / 'oof_predictions.tsv.gz'} ({len(oof_rows)} rows)")

    # ---------------------------------------------------------------- AD.2 as written
    conf_score = next(r for r in score_rows if r["confirmatory"] == "true")
    conf_contrast = next(r for r in contrast_rows if r["confirmatory"] == "true")
    diff = float(conf_contrast["paired_difference"])
    includes_zero = conf_contrast["interval_includes_zero"] == "true"
    ad2 = {
        "prediction": (
            "AD.2: on the reporter endpoint the frozen probe does not beat HyenaDNA's 0.2162 by more than "
            "0.05, and its paired interval includes 0 difference."
        ),
        "confirmatory_cell": conf_score["cell_id"],
        "alphagenome_score": float(conf_score["primary_score"]),
        "alphagenome_ci": [float(conf_score["ci_low"]), float(conf_score["ci_high"])],
        "n_elements": int(conf_score["n_elements"]),
        "n_blocks": int(conf_score["n_blocks"]),
        "hyenadna_score": float(conf_contrast["reference_score"]),
        "paired_difference": diff,
        "paired_ci": [float(conf_contrast["paired_ci_low"]), float(conf_contrast["paired_ci_high"])],
        "two_sided_bootstrap_p": float(conf_contrast["two_sided_bootstrap_p"]),
        "clause_1_difference_at_most_plus_0_05": bool(diff <= 0.05),
        "clause_2_paired_interval_includes_zero": bool(includes_zero),
        "met": bool(diff <= 0.05 and includes_zero),
    }
    ctrl_contrast = next(
        r
        for r in contrast_rows
        if r["cell_id"] == conf_score["cell_id"]
        and r["scope"] == "own"
        and int(r["bootstrap_seed"]) == BOOTSTRAP_SEED_PRIMARY
        and r["against"] == "control"
    )
    receipt = {
        "stage": "ad_arm2_reporter_frozen_probe_fit_and_evaluation",
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "prespecification_sha256": PRESPEC_SHA,
        "addendum_01_sha256": ADDENDUM_SHA,
        "prespecification_written_before_any_fit": True,
        "dataset_id": "gse281364",
        "endpoint": "reporter_allele_effect",
        "elements": len(elements),
        "long_range_blocks": 239,
        "outer_folds": 5,
        "contexts": list(CONTEXTS),
        "biological_donors": 0,
        "experimental_replicates_per_context": 4,
        "unit_of_replication": "experimental_replicate_not_donor",
        "seeds": list(SEEDS),
        "seeds_are_optimisation_variability_not_biological_replicates": True,
        "split_source": str(ROW_UNIVERSE),
        "split_source_note": (
            "folds and blocks come from the row authority, not the sequence fixture, whose outer_fold "
            "column is the obsolete 6-kb assignment and disagrees for 830 of 1033 elements"
        ),
        "head": "published dual-projection delta_ridge, imported helpers, unchanged",
        "component_that_changed": "backbone only",
        "gpu_card": args.card,
        "reproduction_gate": gate,
        "cells": [
            {"cell_id": c[0], "input_length_bp": c[1], "pooling": c[3], "confirmatory": c[4]}
            for c in CELLS
        ],
        "embeddings": {str(k): {kk: v[kk] for kk in ("path", "sha256", "card", "pool_bin_lo", "pool_bin_hi")} for k, v in npz.items()},
        "rejected_elements": rejected_elements,
        "matched_element_count": int(matched_keep.sum()),
        "bootstrap_seed_primary": BOOTSTRAP_SEED_PRIMARY,
        "bootstrap_seed_published_campaign": BOOTSTRAP_SEED_PUBLISHED,
        "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
        "bootstrap_unit": "long_range_block_id",
        "confirmatory_contrast": {
            "cell_id": conf_score["cell_id"],
            "against": "hyenadna/delta_ridge",
            "holm_family": "three allele-effect endpoints at family-wise 0.05; this is one member and the "
            "family does not close until the caQTL and allelic-imbalance endpoints are run",
        },
        "secondary_family": "ad_arm2_reporter_probe_configuration (BH)",
        "AD_2": ad2,
        "contrast_against_allele_identity_control": {
            "paired_difference": float(ctrl_contrast["paired_difference"]),
            "paired_ci": [
                float(ctrl_contrast["paired_ci_low"]),
                float(ctrl_contrast["paired_ci_high"]),
            ],
            "interval_includes_zero": ctrl_contrast["interval_includes_zero"] == "true",
            "two_sided_bootstrap_p": float(ctrl_contrast["two_sided_bootstrap_p"]),
        },
        "terms": "noncommercial_local_alphagenome_weights_derivatives_inherit",
        "hosted_atlas_or_api_outputs_used": False,
        "posthoc_calibration": False,
        "held_outcomes_used_for_fit_or_tuning": False,
        "champion_claim": False,
        "cas13_outcomes_used": False,
    }
    (out / "receipt.json").write_text(json.dumps(receipt, indent=1, sort_keys=True, default=str) + "\n")
    log(f"wrote {out / 'receipt.json'}")
    print(json.dumps({"AD_2": ad2}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

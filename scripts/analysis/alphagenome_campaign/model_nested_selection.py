"""Read the five-seed, four-fold adaptation fits and perform nested selection.

Every adapter number deposited before this came from one seed on one validation
fold, so recipe ordering was confounded with seed noise and with which
chromosomes happened to validate. This reads the 240 fits and reports three
separate quantities that a single-seed table cannot distinguish:

  1. what one run of a recipe gives (mean over seeds of the per-seed metric),
     with the across-seed spread, which is the quantity the earlier single-seed
     numbers actually estimated;
  2. what averaging five runs gives (the seed-ensemble prediction), which is a
     different and better-performing object and must not be quoted as if it were
     the first;
  3. what a recipe gives when it was chosen without seeing the fold it is scored
     on: leave-one-fold-out selection over the four validation folds.

Fold 0 is never read. Selection for outer fold f uses only the other three
validation folds, so no fold contributes to both choosing and scoring a recipe.
The native comparator is calibrated on each fit's own three training folds with
a zero intercept, the existing convention, so the adapter and the native score
are compared on identical rows under identical calibration exposure.

Nothing here is a finalist. The single confirmation on fold 0 is a separate step
that this script deliberately does not perform.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model"
NATIVE_2048 = BASE / "native_lengths_21775189/native_2048.tsv"
NATIVE_1MB = BASE / "native/native_1048576.tsv"
CLOSED_FOLD = 0
VALIDATION_FOLDS = (1, 2, 3, 4)
N_DRAWS = 2000
SEED = 20260917


def spearman(y, p) -> float:
    if len(y) < 3 or np.ptp(y) == 0 or np.ptp(p) == 0:
        return float("nan")
    return float(np.corrcoef(rankdata(y), rankdata(p))[0, 1])


def fisher_z_mean(values) -> float:
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if not v.size:
        return float("nan")
    return float(np.tanh(np.mean(np.arctanh(np.clip(v, -0.999999, 0.999999)))))


def adjust_bh(p) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    order = sorted(range(len(p)), key=lambda i: (p[i], i))
    out = np.ones(len(p))
    running = 1.0
    for rank in range(len(order), 0, -1):
        i = order[rank - 1]
        running = min(running, p[i] * len(p) / rank)
        out[i] = running
    return out


def zero_intercept_slope(x, y) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    denom = float(x @ x)
    if not np.isfinite(denom) or denom <= 1e-20:
        raise ValueError("Degenerate native calibration design")
    return float((x @ y) / denom)


def load_native(path: Path, column: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", usecols=["key", "heldout_fold", "block_1mb", "chr",
                                                 "beta_alt", column])
    frame = frame.loc[frame.heldout_fold.isin(VALIDATION_FOLDS)].copy()
    if (frame.heldout_fold == CLOSED_FOLD).any():
        raise ValueError("Native comparator table leaked the closed fold")
    if frame.key.duplicated().any():
        raise ValueError(f"Repeated variant in {path.name}")
    return frame.set_index("key")


def read_fits(sweeps: list[Path]) -> pd.DataFrame:
    rows = []
    for folder in sweeps:
        runs = json.loads((folder / "runs.json").read_text())
        for run in runs:
            if run.get("status") != "completed_fixed_exposure":
                continue
            path = Path(run["output"]) / "validation_predictions.tsv"
            split = json.loads((Path(run["output"]) / "split.json").read_text())
            if split["held_fold"] != CLOSED_FOLD or CLOSED_FOLD in split["training_folds"]:
                raise ValueError(f"{run['name']} exposed the closed fold")
            if not split["held_fold_never_trained_or_validated"]:
                raise ValueError(f"{run['name']} did not keep the closed fold closed")
            frame = pd.read_csv(path, sep="\t")
            if frame.validation_fold.nunique() != 1:
                raise ValueError(f"{run['name']} mixed validation folds")
            fold = int(frame.validation_fold.iloc[0])
            if fold != split["validation_fold"] or fold == CLOSED_FOLD:
                raise ValueError(f"{run['name']} fold disagreement")
            frame["key"] = frame.variant_id.str.removeprefix("chr")
            rows.append({"recipe": run["recipe"], "seed": int(run["seed"]), "fold": fold,
                         "training_folds": tuple(split["training_folds"]),
                         "predictions": frame.set_index("key")})
    if not rows:
        raise ValueError("No completed fits found")
    return pd.DataFrame(rows)


def block_bootstrap_pairs(diff, blocks, n_draws, rng):
    """Paired mean difference resampled over whole blocks; both arms move together."""
    order = pd.Series(np.arange(len(diff)), index=np.asarray(blocks))
    groups = [g.to_numpy() for _, g in order.groupby(level=0)]
    d = np.asarray(diff, dtype=float)
    draws = np.empty(n_draws)
    n = len(groups)
    for i in range(n_draws):
        pick = rng.integers(0, n, n)
        idx = np.concatenate([groups[j] for j in pick])
        draws[i] = float(np.mean(d[idx]))
    return draws


def contrast(diff, blocks, label, rng):
    observed = float(np.mean(diff))
    draws = block_bootstrap_pairs(diff, blocks, N_DRAWS, rng)
    centred = draws - np.mean(draws)
    p = float((np.sum(np.abs(centred) >= abs(observed)) + 1) / (N_DRAWS + 1))
    return {"contrast": label, "n": int(len(diff)), "blocks": int(len(set(blocks))),
            "mean_squared_error_improvement": observed,
            "ci_low95": float(np.quantile(draws, 0.025)),
            "ci_high95": float(np.quantile(draws, 0.975)),
            "p_nominal_centred_block_bootstrap": p}


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    fits = read_fits(sorted(args.sweeps))
    native2k = load_native(NATIVE_2048, "local_atac_liver")
    native1m = load_native(NATIVE_1MB, "local_atac_liver")

    per_cell, comparator_rows = [], []
    ensembles, natives = {}, {}
    for (fold, recipe), part in fits.groupby(["fold", "recipe"], sort=True):
        keys = None
        for pred in part.predictions:
            keys = set(pred.index) if keys is None else keys & set(pred.index)
        keys &= set(native2k.index) & set(native1m.index)
        shared = sorted(keys)
        if len(shared) < 100:
            raise ValueError(f"fold {fold} recipe {recipe}: only {len(shared)} shared rows")
        y = native2k.loc[shared, "beta_alt"].to_numpy(dtype=float)
        blocks = native2k.loc[shared, "block_1mb"].to_numpy()
        training = part.training_folds.iloc[0]
        if fold in training or CLOSED_FOLD in training:
            raise ValueError("Training folds overlap the validation fold or the closed fold")

        if fold not in natives:
            entry = {}
            for name, table, column in (("native_atac_2048", native2k, "local_atac_liver"),
                                        ("native_atac_1048576", native1m, "local_atac_liver")):
                train_mask = table.heldout_fold.isin(training)
                slope = zero_intercept_slope(table.loc[train_mask, column],
                                             table.loc[train_mask, "beta_alt"])
                entry[name] = {"slope": slope, "training_n": int(train_mask.sum()),
                               "prediction": table.loc[shared, column].to_numpy(dtype=float) * slope}
                comparator_rows.append({"fold": fold, "comparator": name, "training_folds": str(training),
                                        "training_n": int(train_mask.sum()), "zero_intercept_slope": slope})
            natives[fold] = entry

        stack = []
        for row in part.itertuples(index=False):
            pred = row.predictions.loc[shared, "predicted_beta"].to_numpy(dtype=float)
            observed = row.predictions.loc[shared, "observed_beta"].to_numpy(dtype=float)
            if not np.allclose(observed, y, atol=1.5e-6, rtol=1.5e-6):
                raise ValueError(f"{recipe} seed {row.seed} fold {fold}: label disagreement with source")
            stack.append(pred)
            per_cell.append({"fold": fold, "recipe": recipe, "seed": row.seed, "n": len(shared),
                             "mean_squared_error": float(np.mean((pred - y) ** 2)),
                             "root_mean_squared_error": float(np.sqrt(np.mean((pred - y) ** 2))),
                             "signed_spearman": spearman(y, pred)})
        ensembles[(fold, recipe)] = {"prediction": np.mean(np.vstack(stack), axis=0),
                                     "y": y, "blocks": blocks, "keys": shared,
                                     "n_seeds": len(stack)}

    cells = pd.DataFrame(per_cell).sort_values(["fold", "recipe", "seed"])
    cells.to_csv(args.out / "per_seed_cells.tsv", sep="\t", index=False)

    # (1) one run, and (2) the five-run ensemble, kept as separate quantities.
    recipe_rows = []
    for (fold, recipe), part in cells.groupby(["fold", "recipe"], sort=True):
        ens = ensembles[(fold, recipe)]
        recipe_rows.append({
            "fold": fold, "recipe": recipe, "n": int(part.n.iloc[0]), "seeds": len(part),
            "one_run_mean_MSE": float(part.mean_squared_error.mean()),
            "one_run_sd_MSE": float(part.mean_squared_error.std(ddof=1)),
            "one_run_min_MSE": float(part.mean_squared_error.min()),
            "one_run_max_MSE": float(part.mean_squared_error.max()),
            "one_run_mean_signed_spearman": float(part.signed_spearman.mean()),
            "one_run_sd_signed_spearman": float(part.signed_spearman.std(ddof=1)),
            "seed_ensemble_MSE": float(np.mean((ens["prediction"] - ens["y"]) ** 2)),
            "seed_ensemble_signed_spearman": spearman(ens["y"], ens["prediction"]),
            "native_atac_2048_MSE": float(np.mean((natives[fold]["native_atac_2048"]["prediction"] - ens["y"]) ** 2)),
            "native_atac_1048576_MSE": float(np.mean((natives[fold]["native_atac_1048576"]["prediction"] - ens["y"]) ** 2)),
        })
    recipes = pd.DataFrame(recipe_rows)
    recipes.to_csv(args.out / "per_fold_recipes.tsv", sep="\t", index=False)
    pd.DataFrame(comparator_rows).drop_duplicates().to_csv(
        args.out / "native_training_calibration.tsv", sep="\t", index=False)

    macro = recipes.groupby("recipe").agg(
        folds=("fold", "nunique"),
        one_run_mean_MSE=("one_run_mean_MSE", "mean"),
        one_run_mean_sd_across_seeds=("one_run_sd_MSE", "mean"),
        seed_ensemble_MSE=("seed_ensemble_MSE", "mean")).reset_index()
    macro["one_run_macro_signed_spearman"] = [
        fisher_z_mean(recipes.loc[recipes.recipe.eq(r), "one_run_mean_signed_spearman"])
        for r in macro.recipe]
    macro["seed_ensemble_macro_signed_spearman"] = [
        fisher_z_mean(recipes.loc[recipes.recipe.eq(r), "seed_ensemble_signed_spearman"])
        for r in macro.recipe]
    macro = macro.sort_values("one_run_mean_MSE")
    macro.to_csv(args.out / "macro_recipes.tsv", sep="\t", index=False)

    # (3) leave-one-fold-out selection: fold f never informs the choice scored on f.
    rng = np.random.default_rng(SEED)
    selection, contrasts = [], []
    for fold in sorted(recipes.fold.unique()):
        other = recipes.loc[recipes.fold.ne(fold)]
        ranked = other.groupby("recipe").one_run_mean_MSE.mean().sort_values()
        chosen = str(ranked.index[0])
        scored = recipes.loc[recipes.fold.eq(fold) & recipes.recipe.eq(chosen)].iloc[0]
        ens = ensembles[(fold, chosen)]
        row = {"outer_fold": fold, "selected_on_folds": str(sorted(set(recipes.fold) - {fold})),
               "selected_recipe": chosen,
               "selection_margin_MSE": float(ranked.iloc[1] - ranked.iloc[0]),
               "runner_up": str(ranked.index[1]),
               "outer_one_run_mean_MSE": float(scored.one_run_mean_MSE),
               "outer_seed_ensemble_MSE": float(scored.seed_ensemble_MSE),
               "outer_native_2048_MSE": float(scored.native_atac_2048_MSE),
               "outer_native_1048576_MSE": float(scored.native_atac_1048576_MSE)}
        selection.append(row)
        for name in ("native_atac_2048", "native_atac_1048576"):
            native_pred = natives[fold][name]["prediction"]
            diff = (native_pred - ens["y"]) ** 2 - (ens["prediction"] - ens["y"]) ** 2
            c = contrast(diff, ens["blocks"], f"{chosen}__minus__{name}", rng)
            c.update({"outer_fold": fold, "arm": chosen, "comparator": name,
                      "arm_form": "seed_ensemble_of_5", "selection": "leave_one_fold_out"})
            contrasts.append(c)
    pd.DataFrame(selection).to_csv(args.out / "nested_selection.tsv", sep="\t", index=False)

    # Complete exploratory family: every recipe against both comparators, every fold.
    for fold in sorted(recipes.fold.unique()):
        for recipe in sorted(recipes.recipe.unique()):
            ens = ensembles[(fold, recipe)]
            for name in ("native_atac_2048", "native_atac_1048576"):
                native_pred = natives[fold][name]["prediction"]
                diff = (native_pred - ens["y"]) ** 2 - (ens["prediction"] - ens["y"]) ** 2
                c = contrast(diff, ens["blocks"], f"{recipe}__minus__{name}", rng)
                c.update({"outer_fold": fold, "arm": recipe, "comparator": name,
                          "arm_form": "seed_ensemble_of_5", "selection": "none_complete_family"})
                contrasts.append(c)
    frame = pd.DataFrame(contrasts)
    family = frame.selection.eq("none_complete_family")
    frame.loc[family, "BH_q_complete_exploratory_family"] = adjust_bh(
        frame.loc[family, "p_nominal_centred_block_bootstrap"])
    frame.to_csv(args.out / "contrasts.tsv", sep="\t", index=False)

    best = macro.iloc[0]
    winners = {int(f): str(recipes.loc[recipes.fold.eq(f)].sort_values("one_run_mean_MSE").recipe.iloc[0])
               for f in sorted(recipes.fold.unique())}
    summary = {
        "status": "nested_selection_complete_fold0_closed",
        "closed_fold": CLOSED_FOLD,
        "closed_fold_rows_read": 0,
        "validation_folds": sorted(int(f) for f in recipes.fold.unique()),
        "seeds": sorted(int(s) for s in cells.seed.unique()),
        "recipes": int(cells.recipe.nunique()),
        "completed_fits": int(len(cells)),
        "planned_fits": len(VALIDATION_FOLDS) * 5 * 12,
        "best_recipe_by_one_run_macro_MSE": str(best.recipe),
        "best_one_run_mean_MSE": float(best.one_run_mean_MSE),
        "best_mean_across_seed_sd_MSE": float(best.one_run_mean_sd_across_seeds),
        "best_seed_ensemble_MSE": float(best.seed_ensemble_MSE),
        "per_fold_best_recipe": winners,
        "same_recipe_best_in_every_fold": len(set(winners.values())) == 1,
        "leave_one_fold_out_selected": [r["selected_recipe"] for r in selection],
        "bootstrap_draws": N_DRAWS,
        "bootstrap_seed": SEED,
        "resampling_unit": "historical_1Mb_bins_within_fold_not_verified_LD_units",
        "one_run_versus_ensemble": "the one-run mean is the quantity the earlier single-seed tables estimated; the seed ensemble is a different, stronger object and is labelled separately everywhere",
        "native_calibration": "zero intercept, fit on each cell's own three training folds, never on the scored fold",
        "finalist_declared": False,
        "fold0_confirmation_performed": False,
        "adopted": False,
        "sources": {str(p.relative_to(ROOT)): hashlib.sha256((p / "runs.json").read_bytes()).hexdigest()
                    for p in sorted(args.sweeps)},
    }
    (args.out / "selection_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweeps", type=Path, nargs="+", required=True)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())

"""Independent fixed selected analytical fits; no neural training or reselection."""
import time
import numpy as np

GRID = np.logspace(-5, 1, 13)
GRID_KEYS = {f"c{value:g}": float(value) for value in GRID}


def kernel(raw, train, keep=None):
    values = np.asarray(raw[:, keep] if keep is not None else raw, float)
    sd = values[train].std(axis=0)
    selected = sd > 1e-10
    standardized = (values[:, selected] - values[train][:, selected].mean(axis=0)) / sd[selected]
    return standardized @ standardized[train].T / max(1, selected.sum())


def solve_ridge(cross, train, target, coefficient):
    mean = target.mean(axis=0)
    matrix = cross[train]
    penalty = coefficient * max(float(np.linalg.eigvalsh(matrix).max()), 1e-12)
    return cross @ np.linalg.solve(matrix + penalty * np.eye(len(train)), target - mean) + mean


def compare(predictions, name, expected, donors, regions, records, tolerance=1e-6):
    observed = predictions[name][np.ix_(donors, regions)]
    wanted = expected[np.ix_(donors, regions)]
    delta = float(np.max(np.abs(observed - wanted)))
    if not np.isfinite(delta) or delta > tolerance:
        raise ValueError(f"Selected analytical baseline differs: {name}; max={delta}")
    records.append({"model": name, "donors": len(donors), "regions": len(regions), "max_absolute_difference": delta})


def reconstruct(x, sequence, y, fold, reg, predictions, selection, fits, gene_chrom, seconds_cap=3000):
    start = time.monotonic()
    records = []
    for mode in ("trained_regions", "held_regions"):
        regions = np.arange(len(reg)) if mode == "trained_regions" else np.flatnonzero(reg.region_role.eq("train"))
        evaluation = np.arange(len(reg)) if mode == "trained_regions" else np.flatnonzero(reg.region_role.eq("held"))
        seq_kernel = kernel(sequence, regions)
        for held in range(5):
            train = np.flatnonzero(fold != held); test = np.flatnonzero(fold == held)
            target = y[np.ix_(train, regions)]
            relevant_fits = [row for row in fits if row["mode"] == mode and row["held_fold"] == held]
            coefficients = {row["sequence_mean_c"] for row in relevant_fits}
            if mode == "trained_regions":
                assert coefficients == {None}
                baseline = target.mean(axis=0)
            else:
                assert len(coefficients) == 1
                coefficient = float(next(iter(coefficients)))
                assert any(abs(coefficient - c) < 1e-14 for c in GRID)
                baseline = solve_ridge(seq_kernel, regions, target.mean(axis=0)[:, None], coefficient)[:, 0]
            sequence_only = solve_ridge(seq_kernel, regions, target.mean(axis=0)[:, None], .01)[:, 0]
            for arm in ("all_rna", "annotated_all_rna", "chromosome", "random_chromosome"):
                for name, values in (("training_baseline", baseline), ("sequence_only", sequence_only)):
                    compare(predictions[mode], arm + "|" + name, np.broadcast_to(values, y.shape), test, evaluation, records)
            for arm in ("all_rna", "annotated_all_rna"):
                if time.monotonic() - start > seconds_cap:
                    return {"status": "partial_fixed_analytical_check_time_bound", "checked": records,
                            "selection_or_NN_refit_performed": False, "seconds": time.monotonic() - start}
                keep = np.ones(x.shape[1], bool) if arm == "all_rna" else gene_chrom != "unknown"
                rna_kernel = kernel(x, train, keep)
                chosen = selection.loc[selection["mode"].eq(mode) & selection.held_fold.eq(held)
                                       & selection.arm.eq(arm) & selection.group.eq("all")]
                for kind in ("rna_global", "additive_linear", "product_kernel"):
                    row = chosen.loc[chosen.model.eq(kind)]
                    assert len(row) == 1
                    coefficient = float(row.c.iloc[0])
                    assert any(abs(coefficient - c) < 1e-14 for c in GRID)
                    residual = target - target.mean(axis=0)
                    if kind in ("rna_global", "additive_linear"):
                        offset = np.full(len(reg), target.mean()) if kind == "rna_global" else baseline
                        donor = solve_ridge(rna_kernel, train, residual.mean(axis=1, keepdims=True), coefficient)[:, 0]
                        estimate = offset[None] + donor[:, None]
                    else:
                        regional = seq_kernel + 1
                        lv, lu = np.linalg.eigh(rna_kernel[train]); sv, su = np.linalg.eigh(regional[regions])
                        lv = np.maximum(lv, 0); sv = np.maximum(sv, 0)
                        penalty = coefficient * max(float(lv.max() * sv.max()), 1e-12)
                        fitted = lu @ ((lu.T @ residual @ su) / (lv[:, None] * sv[None] + penalty)) @ su.T
                        estimate = baseline[None] + rna_kernel @ fitted @ regional.T
                    compare(predictions[mode], arm + "|" + kind, estimate, test, evaluation, records)
                if mode != "trained_regions":
                    continue
                values, vectors = np.linalg.eigh(rna_kernel[train])
                order = np.argsort(values)[::-1]; values = values[order]; vectors = vectors[:, order]
                admitted = values > max(values[0] * 1e-10, 1e-12)
                scores = rna_kernel @ vectors[:, admitted] / np.sqrt(values[admitted])[None]
                standardized_target = (target - target.mean(0)) / np.maximum(target.std(0), 1e-8)
                directions, _, _ = np.linalg.svd(scores[train].T @ standardized_target, full_matrices=False)
                kernels = {"full_rna_ridge": rna_kernel}
                estimates = {}
                for family in ("full_rna_ridge", "pcr", "pls_svd_ridge", "selected_rna_linear"):
                    rows = chosen.loc[chosen.model.eq(family)].sort_values("region_index_in_panel")
                    assert len(rows) == 128
                    np.testing.assert_array_equal(rows.region_index_in_panel.to_numpy(int), np.arange(128))
                    recipes = rows.recipe.tolist()
                    for recipe in dict.fromkeys(recipes):
                        if recipe in estimates:
                            continue
                        terms = recipe.split(":"); kind = terms[0]
                        if kind == "full_rna_ridge":
                            cross = rna_kernel; coefficient = GRID_KEYS[terms[1]]
                        else:
                            requested = int(terms[1][1:]); key = ":".join(terms[:2])
                            if kind == "pcr":
                                assert requested in (4, 8, 16, 32, 64)
                                latent = scores[:, :min(requested, scores.shape[1])]
                                coefficient = GRID_KEYS[terms[2]]
                            else:
                                assert kind == "pls_svd_ridge" and requested in (2, 4, 8, 10, 16)
                                latent = scores @ directions[:, :min(requested, directions.shape[1])]
                                coefficient = float(terms[2][1:]); assert coefficient in (.001, .01, .1, 1.)
                            if key not in kernels:
                                kernels[key] = latent @ latent[train].T
                            cross = kernels[key]
                        estimates[recipe] = solve_ridge(cross, train, target, coefficient)
                    estimate = np.column_stack([estimates[recipe][:, j] for j, recipe in enumerate(recipes)])
                    compare(predictions[mode], arm + "|" + family, estimate, test, evaluation, records)
    return {"status": "fixed_selected_analytical_baselines_reconstructed", "checked": records,
            "scope": "Every training-baseline/sequence-only arm; unrestricted and annotated RNA analytical baselines only. Selection was not repeated; chromosome/random RNA analytical fits and all NN fits were not refit.",
            "selection_or_NN_refit_performed": False, "seconds": time.monotonic() - start}

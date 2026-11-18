#!/usr/bin/env python3
"""Fixed-fold hepatocyte refits (20261118 candidate, N19).

Updates over final_prediction_refits.py:
1. Adds `shuffle_cross_donor_{43,44,45}` scenarios that permute RNA rows
   strictly across training metacells/donors inside each outer training fold
   (leaving held-out test donors untouched and reusing exact ATAC Cicero caches).
2. Evaluates donor-mean-residualized within-donor MSE (`within_donor_squared_error`
   and `within_donor_standardized_squared_error`) alongside total OOF MSE for every
   held-out donor/gene/fold.
3. Supports single-command per-fold execution across `--array=0-4` into
   `Analysis/ATAC/Regulatory_Methods/results/20261118_prediction_refits/`.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits


A = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))
from promoter_lookup import PromoterLookup
from ridge_path import tune_path

WIDE = A / "results/20261024/ridge_grid_sensitivity_r1"
PREPARED = WIDE / "hepatocyte/prepared"
INPUTS_V2 = A / "results/20261025/final_masld_r1/prediction_inputs_v2"
OUT_ROOT = A / "results/20261118_prediction_refits"
GRID = [0.1, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]

SCENARIOS_20261118 = [
    {"name": "baseline", "source": "baseline", "cross_donor_seed": None},
    {"name": "pairing_shuffle_43", "source": "pairing_shuffle_43", "cross_donor_seed": None},
    {"name": "pairing_shuffle_44", "source": "pairing_shuffle_44", "cross_donor_seed": None},
    {"name": "pairing_shuffle_45", "source": "pairing_shuffle_45", "cross_donor_seed": None},
    {"name": "shuffle_cross_donor_43", "source": "baseline", "cross_donor_seed": 43},
    {"name": "shuffle_cross_donor_44", "source": "baseline", "cross_donor_seed": 44},
    {"name": "shuffle_cross_donor_45", "source": "baseline", "cross_donor_seed": 45},
]

CORE = None
CTX = None
LIMIT = None


def dump(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def load_core():
    global CORE
    source = PREPARED / "source_snapshot/network_compare.py"
    spec = importlib.util.spec_from_file_location("final_frozen_core", source)
    CORE = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(CORE)
    CORE.require_compute()
    receipt = json.loads((PREPARED / "preparation.json").read_text())
    expected = [v for k, v in receipt["source_inventory"].items() if Path(k).name == source.name]
    assert expected == [hashlib.sha256(source.read_bytes()).hexdigest()]
    return receipt


def init_worker() -> None:
    global LIMIT
    LIMIT = threadpool_limits(limits=1)


def fixed_donor_splits(donors, outer_map, inner_map, fold, omitted=""):
    outer = np.array([outer_map[d] for d in donors], int)
    train = np.flatnonzero((outer != fold) & (donors != omitted))
    test = np.flatnonzero(outer == fold)
    inner = np.array([inner_map[str(fold)][d] for d in donors[train]], int)
    assert len(set(inner)) == 3 and not set(donors[train]) & set(donors[test])
    return outer, train, test, inner


def save_predictions_20261118(
    rows, fold, gene, endpoint, method, indices, metadata, y, prediction, variance, within_train_var, alpha
):
    held_donors = metadata.iloc[indices].donor_id.to_numpy()
    observed = y[indices]
    for donor in sorted(set(held_donors)):
        selected = held_donors == donor
        obs_d = observed[selected]
        pred_d = prediction[selected]
        mse = float(np.mean((pred_d - obs_d) ** 2))
        obs_res = obs_d - float(np.mean(obs_d))
        pred_res = pred_d - float(np.mean(pred_d))
        within_mse = float(np.mean((pred_res - obs_res) ** 2))
        within_test_var = float(np.mean(obs_res ** 2))
        within_std = float(within_mse / within_train_var) if within_train_var > 1e-12 else np.nan
        rows.append(
            dict(
                outer_fold=fold,
                gene_id=gene,
                endpoint=endpoint,
                method=method,
                donor_id=donor,
                n_metacells=int(selected.sum()),
                observed=float(obs_d.mean()),
                predicted=float(pred_d.mean()),
                squared_error=mse,
                standardized_squared_error=mse / variance,
                within_donor_squared_error=within_mse,
                within_donor_standardized_squared_error=within_std,
                within_donor_test_variance=within_test_var,
                within_donor_train_variance=within_train_var,
                ridge_alpha=alpha,
            )
        )


def fit_gene(job):
    number, gene = job
    c = CTX
    gid = gene["gene_id"]
    rows, links, tuning = [], [], []
    if gid not in c["rna_lookup"]:
        return number, gid, rows, links, tuning
    prom, candidates = c["geometry"].query(gene, 250000, 2000)
    if not len(prom) or not len(candidates):
        return number, gid, rows, links, tuning
    promoter = c["promoter"].promoter(prom)
    train, test, inner, donors = [c[k] for k in ["train", "test", "inner", "donors"]]
    if np.var(promoter[train]) <= 1e-12:
        return number, gid, rows, links, tuning
    y = CORE.log_columns(c["rna"], [c["rna_lookup"][gid]]).ravel()
    weights = CORE.donor_weights(donors[train])
    variance = float(np.average((y[train] - np.average(y[train], weights=weights)) ** 2, weights=weights))
    if variance <= 1e-12:
        return number, gid, rows, links, tuning

    y_tr = y[train]
    d_tr = donors[train]
    y_tr_res = np.empty_like(y_tr)
    for d_id in set(d_tr):
        m_d = d_tr == d_id
        y_tr_res[m_d] = y_tr[m_d] - float(np.mean(y_tr[m_d]))
    within_train_var = float(np.average(y_tr_res ** 2, weights=weights))

    for method in ["promoter_only", "distance", "cicero"]:
        if method == "promoter_only":
            design = promoter[:, None]
            designs = {k: design for k in sorted(set(inner))}
        else:
            chosen, scores = CORE.selected_indices(
                method,
                candidates,
                prom,
                c["peaks"],
                gene,
                c["training_x"],
                promoter[train],
                donors[train],
                CORE.stable_seed(42, (c["fold"], gid)),
                10,
                c["cicero"]["outer"],
            )
            if not len(chosen):
                continue
            design = CORE.design_matrix(c["feature_x"], chosen, promoter, "rna")
            links += [
                dict(
                    outer_fold=c["fold"],
                    gene_id=gid,
                    method=method,
                    rank=rank,
                    peak_id=c["peaks"].iloc[i].peak_id,
                    score=float(score),
                )
                for rank, (i, score) in enumerate(zip(chosen, scores), 1)
            ]
            designs = {}
            for k in sorted(set(inner)):
                it = train[inner != k]
                selected, _ = CORE.selected_indices(
                    method,
                    candidates,
                    prom,
                    c["peaks"],
                    gene,
                    c["inner_x"][k],
                    promoter[it],
                    donors[it],
                    CORE.stable_seed(42, (c["fold"], k, gid)),
                    10,
                    c["cicero"][f"inner{k}"],
                )
                designs[k] = CORE.design_matrix(c["feature_x"], selected, promoter, "rna")
        losses = tune_path(CORE, designs, y, donors, train, inner, GRID)
        alpha = min(GRID, key=lambda a: (losses[a], a))
        prediction = CORE.weighted_ridge(design[train], y[train], donors[train], alpha, design[test])
        save_predictions_20261118(
            rows, c["fold"], gid, "rna", method, test, c["metadata"], y, prediction, variance, within_train_var, alpha
        )
        tuning.append(
            dict(
                outer_fold=c["fold"],
                gene_id=gid,
                method=method,
                ridge_alpha=alpha,
                inner_mse=losses[alpha],
                upper_boundary_mse=losses[GRID[-1]],
                n_features=design.shape[1],
            )
        )
    return number, gid, rows, links, tuning


def run_fold_scenario(inputs_root: Path, sc_cfg: dict, fold: int, output: Path, workers: int, max_genes: int) -> None:
    global CTX
    comp_file = output / "completion.json"
    if comp_file.exists():
        print(f"[skip] {sc_cfg['name']} fold={fold} already completed", flush=True)
        return

    receipt = load_core()
    started = time.monotonic()
    spec = json.loads((inputs_root / "specification.json").read_text())
    assert spec["ridge_alphas"] == GRID
    src_name = sc_cfg["source"]
    sc_meta = json.loads((inputs_root / src_name / "scenario.json").read_text())
    output.mkdir(parents=True, exist_ok=True)
    source = inputs_root / src_name

    aggregate = sparse.load_npz(source / "counts.npz")
    rna = sparse.load_npz(source / "rna.npz")
    meta = CORE.read_table(source / "metadata.tsv")
    assert CORE.matrix_digest(aggregate) == sc_meta["counts_sha256"]
    assert CORE.matrix_digest(rna) == sc_meta["rna_sha256"]

    donors = meta.donor_id.to_numpy()
    outer, train, test, inner = fixed_donor_splits(donors, spec["outer_map"], spec["inner_map"], fold, "")

    cross_seed = sc_cfg.get("cross_donor_seed")
    if cross_seed is not None:
        # Permute RNA rows ONLY across outer training metacells/donors; leave test untouched!
        rng = np.random.default_rng(int(cross_seed) + 1000 * fold)
        perm_all = np.arange(len(donors))
        perm_all[train] = rng.permutation(train)
        assert np.array_equal(perm_all[test], test)
        rna = rna[perm_all].tocsr()

    peaks = CORE.read_table(PREPARED / "peaks.tsv")
    peaks[["start0", "end"]] = peaks[["start0", "end"]].astype(np.int64)
    genes = CORE.read_table(PREPARED / "genes.tsv")
    genes.tss0 = genes.tss0.astype(np.int64)
    rgenes = CORE.read_table(PREPARED / "rna_genes.tsv")

    sys.path.insert(0, str(A / "scripts"))
    import ridge_grid_sensitivity as wider
    wider.CORE = CORE
    cache = {}
    stages = [("outer", train)] + [(f"inner{k}", train[inner != k]) for k in sorted(set(inner))]
    for tag, _indices in stages:
        old = WIDE / f"hepatocyte/fold{fold}/{tag}/source_verification.json"
        old_root = Path(json.loads(old.read_text())["source"])
        cache[tag], _checked = wider.cached_cicero(old_root, output / tag, aggregate, peaks, meta, _indices, 100)

    x = CORE.normalize_sparse(aggregate)
    CTX = dict(
        peaks=peaks,
        geometry=CORE.PeakGeometry(peaks),
        metadata=meta,
        donors=donors,
        train=train,
        test=test,
        inner=inner,
        fold=fold,
        promoter=PromoterLookup(x),
        feature_x=x.tocsc(),
        training_x=x[train].tocsc(),
        inner_x={k: x[train[inner != k]].tocsc() for k in sorted(set(inner))},
        rna=CORE.normalize_sparse(rna).tocsc(),
        rna_lookup={g: i for i, g in enumerate(rgenes.gene_id)},
        cicero=cache,
    )

    jobs = list(enumerate(genes.sort_values("gene_id").to_dict("records"), 1))
    if max_genes:
        jobs = jobs[:max_genes]

    columns = [
        "outer_fold", "gene_id", "endpoint", "method", "donor_id", "n_metacells",
        "observed", "predicted", "squared_error", "standardized_squared_error",
        "within_donor_squared_error", "within_donor_standardized_squared_error",
        "within_donor_test_variance", "within_donor_train_variance", "ridge_alpha",
    ]
    predictions = CORE.TableStream(output / "predictions.tsv", columns)
    links = CORE.TableStream(output / "links.tsv", ["outer_fold", "gene_id", "method", "rank", "peak_id", "score"])
    tuning = CORE.TableStream(
        output / "tuning.tsv",
        ["outer_fold", "gene_id", "method", "ridge_alpha", "inner_mse", "upper_boundary_mse", "n_features"],
    )
    pool = mp.get_context("fork").Pool(workers, initializer=init_worker)
    try:
        for number, gid, rows, edges, tune in pool.imap(fit_gene, jobs, chunksize=4):
            predictions.append(rows)
            links.append(edges)
            tuning.append(tune)
            if number % 500 == 0 or max_genes:
                print(f"[{sc_cfg['name']}] fold={fold} gene={number} rows={predictions.rows}", flush=True)
        pool.close()
        pool.join()
    except BaseException:
        pool.terminate()
        pool.join()
        raise
    finally:
        predictions.close()
        links.close()
        tuning.close()

    if sc_cfg["name"] == "baseline" and not max_genes:
        keys = ["outer_fold", "gene_id", "endpoint", "method", "donor_id"]
        actual = pd.read_csv(output / "predictions.tsv", sep="\t", dtype={"donor_id": str})
        expected = []
        gids = set(actual.gene_id)
        for chunk in pd.read_csv(WIDE / f"hepatocyte/fold{fold}/predictions.tsv", sep="\t", chunksize=100000, dtype={"donor_id": str}):
            expected.append(chunk[chunk.grid.eq("wide") & chunk.method.isin(["promoter_only", "distance", "cicero"]) & chunk.gene_id.isin(gids)])
        expected = pd.concat(expected).sort_values(keys).reset_index(drop=True)
        actual = actual.sort_values(keys).reset_index(drop=True)
        assert actual[keys].equals(expected[keys])
        numeric = ["observed", "predicted", "squared_error", "standardized_squared_error", "ridge_alpha"]
        assert np.allclose(actual[numeric], expected[numeric], rtol=1e-10, atol=1e-11)
        dump(
            output / "baseline_reproduction.json",
            dict(status="passed", rows=len(actual), max_abs_diff=float(np.max(np.abs(actual[numeric].to_numpy() - expected[numeric].to_numpy())))),
        )

    dump(
        comp_file,
        dict(
            status="completed",
            scenario=sc_cfg["name"],
            outer_fold=fold,
            cross_donor_seed=cross_seed,
            n_held_donors=len(set(donors[test])),
            n_training_donors=len(set(donors[train])),
            elapsed_seconds=time.monotonic() - started,
            predictions=predictions.rows,
            technical_pilot=bool(max_genes),
            ridge_alphas=GRID,
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=INPUTS_V2)
    parser.add_argument("--outdir", type=Path, default=OUT_ROOT)
    parser.add_argument("--fold", type=int, default=None, choices=range(5))
    parser.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    parser.add_argument("--max-genes", type=int, default=0)
    args = parser.parse_args()

    fold = args.fold if args.fold is not None else int(os.environ["SLURM_ARRAY_TASK_ID"])
    for sc_cfg in SCENARIOS_20261118:
        out_dir = args.outdir / sc_cfg["name"] / f"fold{fold}"
        run_fold_scenario(args.inputs, sc_cfg, fold, out_dir, args.workers, args.max_genes)


if __name__ == "__main__":
    main()

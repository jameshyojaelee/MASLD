#!/usr/bin/env python3
"""Retrospective reveal-and-refit test on an annotation-fixed chromatin panel.

Outer test participants are never candidates for acquisition. The panel is a
hash-selected 1,000-region sample of regions with an annotated nearby gene;
the hash and cis annotation are read before any H3K27ac outcome.
"""
import argparse
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT/"Analysis/MASLD_Model_Benchmark"
FIX = BENCH/"executions/model-data-064-21079902/fixture"
WEIGHTS = BENCH/"release/masld-liver-chromatin-state-v1.2/weights/chromatin_state_v1_2.npz"
JOIN = BENCH/"executions/gse267145-authoritative-join-21064930/participant_join.tsv"
SEED = 20260922
RANK = 12
PENALTY = 1e-7
LOCAL_PENALTY = 0.1
N_REGIONS = 1000
BUDGETS = (0.1, 0.2, 0.3)


def logcpm(raw):
    return np.log2(raw/raw.sum(1)[:, None]*1e6+1)


def concentration(raw):
    nreg = raw.shape[1]
    lib = raw.sum(1)
    frac = raw/lib[:, None]
    entropy = -np.sum(np.where(frac > 0, frac*np.log(np.maximum(frac, 1e-300)), 0), axis=1)
    sorted_raw = np.sort(raw, axis=1)
    gini = ((2*np.arange(1, nreg+1)-nreg-1)*sorted_raw).sum(1)/(nreg*lib)
    iqr = np.array([np.subtract(*np.percentile(np.log2(row[row > 0]), [75, 25])) for row in raw])
    top5 = sorted_raw[:, -int(round(.05*nreg)):].sum(1)/lib
    return np.column_stack([gini, entropy, iqr, top5])


def panel(weights):
    key = weights["prof_region_key"].astype(str)
    cis = np.asarray(weights["prof_Ocis_cis_idx"], int)
    eligible = np.flatnonzero((cis >= 0).any(1))
    if len(key) != 96460 or len(eligible) < N_REGIONS:
        raise ValueError("Annotation or region axis differs")
    # No prediction, target, reliable-region flag or previous skill enters.
    score = [hashlib.sha256((str(SEED)+"|"+key[j]).encode()).digest() for j in eligible]
    selected = eligible[np.argsort(score)[:N_REGIONS]]
    return selected, key[selected], cis[selected]


def rrr_fit_predict(x, y, tr, query, selected):
    """Dual-space reduced-rank regression; only requested regions materialized."""
    mean = x[tr].mean(0)
    sd = x[tr].std(0)
    sd[sd < 1e-8] = 1.0
    train = (x[tr]-mean)/sd
    other = (x[query]-mean)/sd
    gram = train@train.T/train.shape[1]
    eig, u = np.linalg.eigh(gram)
    eig = np.maximum(eig, 0)
    lam = PENALTY*max(float(eig.max()), 1e-12)
    ym = y[tr].mean(0)
    yc = y[tr]-ym
    uy = u.T@yc
    fitted = u@((eig/(eig+lam))[:, None]*uy)
    fa = fitted@fitted.T
    sv2, left = np.linalg.eigh(fa)
    order = np.argsort(sv2)[::-1][:min(RANK, len(tr)-1)]
    valid = sv2[order] > 1e-10
    order = order[valid]
    if len(order) < 4:
        raise ValueError("Insufficient response rank")
    v = fitted.T@left[:, order]/np.sqrt(sv2[order])
    alpha = u@((1/(eig+lam))[:, None]*uy)
    compressed = alpha@v
    kt = train@train.T/train.shape[1]
    kq = other@train.T/train.shape[1]
    loading = v[selected].T
    in_sample = kt@compressed@loading+ym[selected]
    predicted = kq@compressed@loading+ym[selected]
    return in_sample, predicted


def local_fit_predict(x, residual, tr, query, cis):
    pred = np.zeros((len(query), len(cis)), np.float64)
    for j, genes in enumerate(cis):
        genes = genes[genes >= 0]
        if len(genes) == 0:
            continue
        a, b = x[tr][:, genes], x[query][:, genes]
        mu, sd = a.mean(0), a.std(0)
        sd[sd < 1e-8] = 1.0
        za, zb = (a-mu)/sd, (b-mu)/sd
        target = residual[:, j]
        ym = target.mean()
        gram = za.T@za
        lam = LOCAL_PENALTY*max(float(np.linalg.eigvalsh(gram).max()), 1e-12)
        coef = np.linalg.solve(gram+lam*np.eye(len(genes)), za.T@(target-ym))
        pred[:, j] = zb@coef+ym
    return pred


def fit_profile(x, y, tr, query, cis, selected):
    # Cross-fitted global predictions for the local residual: each training
    # participant's own chromatin is absent from its offset prediction.
    train_idx = np.asarray(tr, int)
    rng = np.random.default_rng(SEED + len(train_idx))
    parts = np.array_split(rng.permutation(train_idx), 3)
    oof = np.zeros((len(train_idx), len(selected)), float)
    locations = {int(person): pos for pos, person in enumerate(train_idx)}
    for part in parts:
        other = np.setdiff1d(train_idx, part)
        _, p = rrr_fit_predict(x, y, other, part, selected)
        oof[[locations[int(person)] for person in part]] = p
    _, global_query = rrr_fit_predict(x, y, train_idx, query, selected)
    residual = y[train_idx][:, selected]-oof
    local_query = local_fit_predict(x, residual, train_idx, query, cis)
    return global_query, global_query+local_query


def diversity_order(x, train, pool):
    mu, sd = x[train].mean(0), x[train].std(0)
    sd[sd < 1e-8] = 1.0
    ztr, zp = (x[train]-mu)/sd, (x[pool]-mu)/sd
    gram = ztr@ztr.T
    eig, u = np.linalg.eigh(gram)
    order = np.argsort(eig)[::-1][:10]
    score = zp@ztr.T@u[:, order]/np.sqrt(np.maximum(eig[order], 1e-8))
    chosen = []
    distance = np.sum((score-score.mean(0))**2, axis=1)
    while len(chosen) < len(pool):
        distance[chosen] = -1
        j = int(np.argmax(distance))
        chosen.append(j)
        distance = np.minimum(distance, np.sum((score-score[j])**2, axis=1))
    return np.array(chosen)


def histology_order(hist, pool, ids):
    labels = [tuple(row) for row in hist[pool]]
    groups = sorted(set(labels))
    selected = []
    while len(selected) < len(pool):
        left = [j for j in range(len(pool)) if j not in selected]
        counts = {g: sum(labels[j] == g for j in selected) for g in groups}
        selected.append(min(left, key=lambda j: (counts[labels[j]], labels[j], ids[pool[j]])))
    return np.asarray(selected)


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.mkdir(parents=True)
    weights = np.load(WEIGHTS, allow_pickle=True)
    selected, keys, cis = panel(weights)
    pd.DataFrame({"region_index": selected, "region_key": keys,
                  "annotation": "at least one TSS within 100 kb; SHA256 outcome-blind sample"}).to_csv(
                      args.out/"fixed_region_panel.tsv", sep="\t", index=False)
    axis = pd.read_csv(FIX/"molecular/participant_axis.tsv", sep="\t")
    ids = axis.participant_id.astype(str).to_numpy()
    if len(ids) != 99 or len(set(ids)) != 99:
        raise ValueError("Expected 99 distinct participants")
    fold = pd.read_csv(FIX/"folds/participant_outer_folds.tsv", sep="\t").set_index("participant_id").loc[ids, "outer_fold"].to_numpy(int)
    join = pd.read_csv(JOIN, sep="\t").set_index("participant_id").loc[ids]
    hist = join[["steatosis", "ballooning", "lobular_inflammation", "fibrosis"]].to_numpy(float)
    rna = np.asarray(np.load(FIX/"molecular/rna_values.npy", mmap_mode="r"), np.float64)
    h3 = np.asarray(np.load(FIX/"molecular/h3k27ac_counts.npy", mmap_mode="r"), np.float64)
    x, yraw = logcpm(rna), logcpm(h3)
    conc = concentration(h3)
    rows, donors = [], []
    for outer in sorted(set(fold)):
        evaluation = np.flatnonzero(fold == outer)
        available = np.flatnonzero(fold != outer)
        # Base and pool assignment depends on donor ID hash only.
        available = np.array(sorted(available, key=lambda j: hashlib.sha256(
            (str(SEED)+"|acquisition|"+ids[j]).encode()).digest()))
        base, pool = available[:45], available[45:]
        if len(pool) < 29:
            raise ValueError("Candidate pool too small for 30% budget")
        Zbase = np.column_stack([np.ones(len(base)), conc[base]])
        coef = np.linalg.lstsq(Zbase, yraw[base], rcond=None)[0]
        y = yraw-np.column_stack([np.ones(len(yraw)), conc])@coef
        all_query = np.concatenate([pool, evaluation])
        base_global_all, base_all = fit_profile(x, y, base, all_query, cis, selected)
        base_global_pool, base_pool = base_global_all[:len(pool)], base_all[:len(pool)]
        base_eval = base_all[len(pool):]
        error_pool = np.sum((y[pool][:, selected]-base_pool)**2, axis=1)
        eval_truth = y[evaluation][:, selected]
        base_loss = np.sum((eval_truth-base_eval)**2, axis=1)
        disagreement = np.mean((base_pool-base_global_pool)**2, axis=1)
        orders = {"histology": histology_order(hist, pool, ids),
                  "rna_diversity": diversity_order(x, base, pool),
                  "global_local_disagreement": np.argsort(-disagreement)}
        rng = np.random.default_rng(SEED+outer)
        for repeat in range(3):
            orders[f"random_{repeat}"] = rng.permutation(len(pool))
        for fraction in BUDGETS:
            k = max(1, int(round(len(ids)*fraction)))
            if k > len(pool):
                raise ValueError("Requested acquisition exceeds pool")
            for name, order in orders.items():
                chosen = pool[order[:k]]
                training = np.concatenate([base, chosen])
                _, refit = fit_profile(x, y, training, evaluation, cis, selected)
                refit_loss = np.sum((eval_truth-refit)**2, axis=1)
                policy = "random" if name.startswith("random_") else name
                repeat = int(name.split("_")[-1]) if policy == "random" else 0
                rows.append({"outer_fold": int(outer), "policy": policy, "random_repeat": repeat,
                             "budget_fraction": fraction, "n_acquired": k, "n_base": len(base),
                             "n_pool": len(pool), "n_evaluation": len(evaluation),
                             "captured_pool_error_fraction": float(error_pool[order[:k]].sum()/error_pool.sum()),
                             "evaluation_base_MSE": float(base_loss.mean()/len(selected)),
                             "evaluation_refit_MSE": float(refit_loss.mean()/len(selected)),
                             "evaluation_MSE_reduction": float((base_loss.mean()-refit_loss.mean())/len(selected))})
                for j, donor in enumerate(evaluation):
                    donors.append({"outer_fold": int(outer), "participant_id": ids[donor],
                                   "policy": policy, "random_repeat": repeat,
                                   "budget_fraction": fraction, "base_SSE": float(base_loss[j]),
                                   "refit_SSE": float(refit_loss[j])})
            print(json.dumps({"completed_fold": int(outer), "budget": fraction,
                              "candidate_policies": len(orders)}), flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(args.out/"acquisition_per_fold.tsv", sep="\t", index=False)
    pd.DataFrame(donors).to_csv(args.out/"acquisition_per_participant.tsv.gz", sep="\t", index=False)
    summary = table.groupby(["budget_fraction", "policy"], sort=True).agg(
        n_folds=("outer_fold", "nunique"), n_fit_comparisons=("outer_fold", "size"),
        mean_error_capture=("captured_pool_error_fraction", "mean"),
        mean_refit_MSE_reduction=("evaluation_MSE_reduction", "mean"),
        mean_base_MSE=("evaluation_base_MSE", "mean")).reset_index()
    summary.to_csv(args.out/"acquisition_summary.tsv", sep="\t", index=False)
    report = {"status": "retrospective single-cohort development", "participants": len(ids),
              "regions": len(selected), "region_selection": "nearby-gene annotation then SHA256, no H3 value or predictive skill",
              "target_unit": "base-training-residual H3K27ac log2 CPM",
              "model": "rank-12 RRR at fixed shrinkage 1e-7 of leading RNA-Gram eigenvalue plus 10-nearby-gene ridge on cross-fitted global residual, local penalty fraction 0.1",
              "model_relation": "same form as released RRR+local head with fixed complexity and 1,000-region evaluation; this is a bounded refit experiment, not the full released 96,460-region recipe",
              "budgets": BUDGETS, "random_repeats": 3, "seed": SEED,
              "environment": {"python": sys.version, "numpy": np.__version__,
                              "pandas": pd.__version__, "platform": platform.platform(),
                              "slurm_job_id": os.environ["SLURM_JOB_ID"]}}
    (args.out/"results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())

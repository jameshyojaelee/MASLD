#!/usr/bin/env python3
"""Re-evaluate saved held-participant profiles beyond recorded histology.

This is development evidence from the already inspected GSE267145 folds.  It
does not refit the RNA model or claim a new external confirmation.
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
BENCH = ROOT / "Analysis/MASLD_Model_Benchmark"
FIX = BENCH / "executions/model-data-064-21079902/fixture"
OOF = BENCH / "executions/chromatin-stable-rrr-forms-20260908T192024Z/out/stable_oof.npz"
TECH = BENCH / "executions/chromatin-stable-rrr-forms-20260908T192024Z/out/technical_features.tsv"
JOIN = BENCH / "executions/gse267145-authoritative-join-21064930/participant_join.tsv"
DECONV = BENCH / "executions/chromatin-pc1-identity-20260906T005837Z/deconv_final_mRNA_fraction.tsv"
DESIGN = Path(__file__).with_name("design.json")
SEED = 20260922
MODELS = {"global_RRR": "pred_rrr_wide_stable", "RRR_local": "pred_rrr_offset_cis_wide"}


def fit_nuisance(train, test, y):
    """Mean/scale of covariates and all region coefficients use training people."""
    mu, sd = train.mean(0), train.std(0)
    sd[sd < 1e-9] = 1.0
    a = np.column_stack([np.ones(len(train)), (train-mu)/sd])
    b = np.column_stack([np.ones(len(test)), (test-mu)/sd])
    coef = np.linalg.lstsq(a, y, rcond=None)[0]
    return b @ coef


def fit_nuisance_ridge(train, test, y, seed):
    """Select a whole-profile penalty inside the outer training participants."""
    grid = (0.0, 0.01, 0.1, 1.0, 10.0, 100.0, float("inf"))
    rng = np.random.default_rng(seed)
    split = np.array_split(rng.permutation(len(train)), 4)
    errors = np.zeros(len(grid))

    def prepare(a, b, target):
        mean, sd = a.mean(0), a.std(0)
        sd[sd < 1e-9] = 1.0
        za, zb = (a-mean)/sd, (b-mean)/sd
        ym = target.mean(0)
        return za, zb, target-ym, ym

    def predict(a, b, target, fraction):
        za, zb, yc, ym = prepare(a, b, target)
        if np.isinf(fraction):
            return np.broadcast_to(ym, (len(b), len(ym)))
        gram = za.T@za
        lam = fraction*max(float(np.linalg.eigvalsh(gram).max()), 1e-12)
        coef = np.linalg.lstsq(gram+lam*np.eye(gram.shape[0]), za.T@yc, rcond=None)[0]
        return zb@coef+ym

    for valid in split:
        internal = np.setdiff1d(np.arange(len(train)), valid)
        for j, fraction in enumerate(grid):
            residual = y[valid]-predict(train[internal], train[valid], y[internal], fraction)
            errors[j] += float(np.sum(residual**2))
    chosen = grid[int(np.argmin(errors))]
    return predict(train, test, y, chosen), chosen


def bootstrap_skill(loss, reference, ids, draws=2000):
    rng = np.random.default_rng(SEED)
    n = len(ids)
    bi = rng.integers(0, n, (draws, n))
    obs = 1.0 - loss.sum()/reference.sum()
    bs = 1.0 - loss[bi].sum(1)/reference[bi].sum(1)
    return {"skill": float(obs), "ci95": np.quantile(bs, [0.025, 0.975]).tolist(),
            "biological_n": n, "bootstrap_unit": "distinct participant",
            "bootstrap_seed": SEED, "bootstrap_draws": draws}


def paired_delta(loss_a, loss_b, reference, draws=2000):
    rng = np.random.default_rng(SEED + 7)
    n = len(loss_a)
    bi = rng.integers(0, n, (draws, n))
    delta = (loss_b.sum() - loss_a.sum())/reference.sum()
    bs = (loss_b[bi].sum(1) - loss_a[bi].sum(1))/reference[bi].sum(1)
    centred = bs - np.mean(bs)
    p = (1 + np.sum(np.abs(centred) >= abs(delta)))/(draws+1)
    return {"delta_skill": float(delta), "ci95": np.quantile(bs, [0.025, 0.975]).tolist(),
            "p_nominal_centred_participant_bootstrap": float(p)}


def holm(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ans = [1.0]*len(values)
    previous = 0.0
    for j, i in enumerate(order):
        previous = max(previous, (len(values)-j)*values[i])
        ans[i] = min(1.0, previous)
    return ans


def diversity_order(x, candidates, k):
    """Farthest first in an outcome-blind RNA PCA space fit on this fold's training people."""
    if len(candidates) <= k:
        return candidates
    candidate = np.asarray(candidates, int)
    centre = x[candidate].mean(0)
    first = int(np.argmax(np.sum((x[candidate]-centre)**2, axis=1)))
    selected = [first]
    distance = np.sum((x[candidate]-x[candidate[first]])**2, axis=1)
    for _ in range(1, k):
        distance[selected] = -1.0
        j = int(np.argmax(distance))
        selected.append(j)
        distance = np.minimum(distance, np.sum((x[candidate]-x[candidate[j]])**2, axis=1))
    return candidate[selected]


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Compute-node allocation required")
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.mkdir(parents=True)
    design = json.loads(DESIGN.read_text())
    if design["seed"] != SEED:
        raise ValueError("Design seed disagrees with evaluation")
    axis = pd.read_csv(FIX / "molecular/participant_axis.tsv", sep="\t")
    ids = axis.participant_id.astype(str).to_numpy()
    if len(ids) != 99 or len(set(ids)) != 99:
        raise ValueError("Participant identity or n changed")
    join = pd.read_csv(JOIN, sep="\t").set_index("participant_id").loc[ids]
    fold = pd.read_csv(FIX / "folds/participant_outer_folds.tsv", sep="\t").set_index("participant_id").loc[ids, "outer_fold"].to_numpy(int)
    tech = pd.read_csv(TECH, sep="\t")
    if len(tech) != len(ids):
        raise ValueError("Technical descriptor rows differ")
    core = join[["steatosis", "ballooning", "lobular_inflammation", "fibrosis"]].to_numpy(float)
    core = np.column_stack([core, (join.sex.to_numpy() == "F").astype(float)])
    if not np.isfinite(core).all():
        raise ValueError("Incomplete histology or sex")
    rna_tech = tech[[c for c in tech if not c.startswith("h3k27ac_")]].to_numpy(float)
    hep = pd.read_csv(DECONV, sep="\t", index_col=0).loc[ids, "Hepatocytes"].to_numpy(float)
    if not np.isfinite(hep).all():
        raise ValueError("Incomplete composition adjustment")
    covars = {"histology_sex": core,
              "histology_sex_RNA_technical": np.column_stack([core, rna_tech]),
              "histology_sex_RNA_technical_hepatocytes": np.column_stack([core, rna_tech, hep])}
    stored = np.load(OOF, allow_pickle=False)
    if not np.array_equal(stored["fold"], fold):
        raise ValueError("Saved predictions and participant fold order differ")
    truth = stored["Cres_oof"].astype(np.float64)
    pred = {name: stored[column].astype(np.float64) for name, column in MODELS.items()}
    if any(p.shape != truth.shape for p in pred.values()):
        raise ValueError("Prediction shape differs from target")
    h3 = np.load(FIX / "molecular/h3k27ac_counts.npy", mmap_mode="r")
    if h3.shape != truth.shape:
        raise ValueError("Measured H3K27ac matrix shape differs")
    logcpm = np.log2(np.asarray(h3, np.float64)/np.asarray(h3.sum(1), np.float64)[:, None]*1e6 + 1)
    baseline = np.zeros_like(truth)
    conditional = {name: np.zeros_like(truth) for name in covars}
    # Four sample-level concentration descriptors. Their regression coefficients
    # below are refitted separately on every outer training fold.
    raw = np.asarray(h3, np.float64)
    lib = raw.sum(1)
    frac = raw/lib[:, None]
    entropy = -np.nansum(np.where(frac > 0, frac*np.log(np.maximum(frac, 1e-300)), 0), axis=1)
    sorted_raw = np.sort(raw, axis=1)
    nreg = raw.shape[1]
    gini = ((2*np.arange(1, nreg+1)-nreg-1)*sorted_raw).sum(1)/(nreg*lib)
    iqr = np.array([np.subtract(*np.percentile(np.log2(row[row > 0]), [75, 25])) for row in raw])
    top5 = sorted_raw[:, -int(round(.05*nreg)):].sum(1)/lib
    conc = np.column_stack([gini, entropy, iqr, top5])
    nuisance_choices = []
    for f in sorted(set(fold)):
        tr, te = fold != f, fold == f
        # Target residualization has the same training-fold convention as the saved lane.
        Zi = np.column_stack([np.ones(tr.sum()), conc[tr]])
        Zte = np.column_stack([np.ones(te.sum()), conc[te]])
        coef = np.linalg.lstsq(Zi, logcpm[tr], rcond=None)[0]
        ytr = logcpm[tr] - Zi @ coef
        cres = logcpm[te] - Zte @ coef
        if not np.allclose(cres, truth[te], atol=4e-5, rtol=0):
            raise ValueError(f"Fold {f}: measured residual does not reproduce saved target")
        baseline[te] = ytr.mean(0)
        for name, c in covars.items():
            conditional[name][te], chosen = fit_nuisance_ridge(c[tr], c[te], ytr, SEED+f)
            nuisance_choices.append({"outer_fold": int(f), "covariates": name,
                                     "inner_selected_penalty_fraction": chosen})
    ref = np.sum((truth-baseline)**2, axis=1)
    losses = {name: np.sum((truth-value)**2, axis=1) for name, value in pred.items()}
    for name, value in conditional.items():
        losses[name] = np.sum((truth-value)**2, axis=1)
    result = {"evidence_state": "inspected single-cohort development", "n_participants": len(ids),
              "n_regions": truth.shape[1], "target_unit": "residual H3K27ac log2 CPM",
              "nuisance_choices": nuisance_choices,
              "skill": {name: bootstrap_skill(loss, ref, ids) for name, loss in losses.items()},
              "contrasts": {}, "source_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (OOF, TECH, DESIGN)},
              "environment": {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
                              "platform": platform.platform(), "slurm_job_id": os.environ["SLURM_JOB_ID"]}}
    pairs = [("global_RRR", "histology_sex_RNA_technical"),
             ("RRR_local", "histology_sex_RNA_technical"),
             ("RRR_local", "global_RRR")]
    for a, b in pairs:
        result["contrasts"][f"{a}_vs_{b}"] = paired_delta(losses[a], losses[b], ref)
    corrected = holm([result["contrasts"][f"{a}_vs_{b}"]["p_nominal_centred_participant_bootstrap"] for a,b in pairs])
    for (a, b), q in zip(pairs, corrected):
        result["contrasts"][f"{a}_vs_{b}"]["p_holm_primary_family"] = q
    # Capture is an acquisition proxy. A separate refit is necessary for measured value.
    xraw = np.load(FIX / "molecular/rna_values.npy", mmap_mode="r")
    x = np.log2(np.asarray(xraw, np.float32)/np.asarray(xraw.sum(1), np.float32)[:, None]*1e6 + 1)
    triage = []
    rng = np.random.default_rng(SEED)
    for f in sorted(set(fold)):
        tr, te = np.flatnonzero(fold != f), np.flatnonzero(fold == f)
        # PCA via training-person Gram; no chromatin value enters the feature map.
        mu, sd = x[tr].mean(0), x[tr].std(0)
        sd[sd < 1e-7] = 1.0
        ztr = (x[tr]-mu)/sd
        zte = (x[te]-mu)/sd
        gram = ztr@ztr.T
        vals, vecs = np.linalg.eigh(gram)
        order = np.argsort(vals)[::-1][:min(10, len(tr)-1)]
        score = zte@ztr.T@vecs[:, order]/np.sqrt(np.maximum(vals[order], 1e-8))
        error = losses["RRR_local"][te]
        disagreement = np.mean((pred["RRR_local"][te]-pred["global_RRR"][te])**2, axis=1)
        for budget in design["chromatin"]["measurement_budgets"]:
            k = max(1, int(round(len(te)*budget)))
            picks = {"rna_diversity": diversity_order(score, np.arange(len(te)), k),
                     "global_local_disagreement": np.argsort(-disagreement)[:k]}
            # Proportional allocation across exact histology strata with deterministic ties.
            labels = [tuple(row) for row in core[te, :4]]
            groups = sorted(set(labels))
            hist_pick = []
            while len(hist_pick) < k:
                eligible = [j for j in range(len(te)) if j not in hist_pick]
                counts = {g: sum(labels[j] == g for j in hist_pick) for g in groups}
                hist_pick.append(min(eligible, key=lambda j: (counts[labels[j]], labels[j], ids[te[j]])))
            picks["histology"] = np.asarray(hist_pick)
            random_capture = []
            for _ in range(1000):
                random_capture.append(float(error[rng.choice(len(te), k, replace=False)].sum()/error.sum()))
            triage.append({"fold": int(f), "budget_fraction": budget, "selected": k,
                           "policy": "random", "capture_fraction": float(np.mean(random_capture)),
                           "random_sd": float(np.std(random_capture, ddof=1))})
            for policy, selected in picks.items():
                triage.append({"fold": int(f), "budget_fraction": budget, "selected": k,
                               "policy": policy, "capture_fraction": float(error[selected].sum()/error.sum()),
                               "random_sd": None})
    pd.DataFrame(triage).to_csv(args.out/"measurement_capture.tsv", sep="\t", index=False)
    (args.out/"results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"n_participants": len(ids), "n_regions": truth.shape[1],
                      "skill": {k: v["skill"] for k,v in result["skill"].items()},
                      "contrasts": result["contrasts"]}, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())

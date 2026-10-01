#!/usr/bin/env python3
"""Held-participant chromatin differences within exact histology/sex strata."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from chromatin_acquisition import concentration, logcpm
from chromatin_development import FIX, OOF, TECH, JOIN, DECONV, holm

ROOT = Path(__file__).resolve().parents[3]
WEIGHTS = ROOT/"Analysis/MASLD_Model_Benchmark/release/masld-liver-chromatin-state-v1.2/weights/chromatin_state_v1_2.npz"
SEED = 20260923


def pairs_from_metadata(join):
    groups = ["outer_fold", "steatosis", "ballooning", "lobular_inflammation", "fibrosis", "sex"]
    pairs = []
    for key, part in join.groupby(groups, sort=True, dropna=False):
        order = sorted(part.participant_id.astype(str), key=lambda p: hashlib.sha256(f"{SEED}|{p}".encode()).digest())
        for offset in range(0, len(order)-1, 2):
            pairs.append({"pair_id": len(pairs), "participant_a": order[offset],
                          "participant_b": order[offset+1], **dict(zip(groups, key))})
    result = pd.DataFrame(pairs)
    if len(result) != 22 or len(set(result.participant_a)|set(result.participant_b)) != 44:
        raise ValueError("Disjoint metadata pair census differs")
    return result


def nuisance_predictions(x, raw_y, descriptors, fold, choices, name):
    out = np.zeros_like(raw_y)
    for f in sorted(set(fold)):
        tr, te = np.flatnonzero(fold != f), np.flatnonzero(fold == f)
        basis = np.column_stack([np.ones(len(tr)), descriptors[tr]])
        coefficient = np.linalg.lstsq(basis, raw_y[tr], rcond=None)[0]
        ytr = raw_y[tr]-basis@coefficient
        mu, sd = x[tr].mean(0), x[tr].std(0)
        sd[sd < 1e-9] = 1
        a, b = (x[tr]-mu)/sd, (x[te]-mu)/sd
        ym = ytr.mean(0)
        gram = a.T@a
        penalty = choices[int(f), name]
        if np.isinf(penalty):
            out[te] = ym
        else:
            lam = penalty*max(float(np.linalg.eigvalsh(gram).max()), 1e-12)
            beta = np.linalg.lstsq(gram+lam*np.eye(len(gram)), a.T@(ytr-ym), rcond=None)[0]
            out[te] = b@beta+ym
    return out


def summary(loss, reference, samples, chrom_samples=None):
    # Matrices have disjoint donor-pair rows and chromosome columns.
    observed = float(1-loss.sum()/reference.sum())
    if chrom_samples is None:
        draws = 1-loss.sum(1)[samples].sum(1)/reference.sum(1)[samples].sum(1)
    else:
        draws = np.array([1-loss[np.ix_(p, c)].sum()/reference[np.ix_(p, c)].sum()
                          for p, c in zip(samples, chrom_samples)])
    return {"skill": observed, "ci95": np.quantile(draws, [.025, .975]).tolist()}


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.mkdir(parents=True)
    axis = pd.read_csv(FIX/"molecular/participant_axis.tsv", sep="\t").participant_id.astype(str)
    join = pd.read_csv(JOIN, sep="\t").set_index("participant_id").loc[axis].reset_index()
    pairs = pairs_from_metadata(join)
    pairs.to_csv(args.out/"pairs_fixed_from_metadata.tsv", sep="\t", index=False)
    positions = {p: i for i, p in enumerate(axis)}
    a = pairs.participant_a.map(positions).to_numpy(int)
    b = pairs.participant_b.map(positions).to_numpy(int)
    stored = np.load(OOF, allow_pickle=False)
    fold = stored["fold"].astype(int)
    if not np.array_equal(fold, join.outer_fold.to_numpy(int)) or not np.array_equal(fold[a], fold[b]):
        raise ValueError("Matched donors do not share target transformation/training fold")
    y = stored["Cres_oof"].astype(float)
    predictions = {"global_RRR": stored["pred_rrr_wide_stable"].astype(float),
                   "RRR_local": stored["pred_rrr_offset_cis_wide"].astype(float)}
    original = json.loads((args.conditional/"results.json").read_text())
    choices = {(x["outer_fold"], x["covariates"]): x["inner_selected_penalty_fraction"]
               for x in original["nuisance_choices"]}
    tech = pd.read_csv(TECH, sep="\t")
    z = np.column_stack([join[["steatosis", "ballooning", "lobular_inflammation", "fibrosis"]].to_numpy(float),
                         join.sex.eq("F").to_numpy(float),
                         tech[[c for c in tech if not c.startswith("h3k27ac_")]].to_numpy(float)])
    h3 = np.asarray(np.load(FIX/"molecular/h3k27ac_counts.npy", mmap_mode="r"), float)
    raw_y, descriptors = logcpm(h3), concentration(h3)
    name = "histology_sex_RNA_technical"
    predictions["technical_nuisance"] = nuisance_predictions(z, raw_y, descriptors, fold, choices, name)
    hep = pd.read_csv(DECONV, sep="\t", index_col=0).loc[axis, "Hepatocytes"].to_numpy(float)
    predictions["technical_plus_hepatocytes"] = nuisance_predictions(
        np.column_stack([z, hep]), raw_y, descriptors, fold, choices, name+"_hepatocytes")
    with np.load(WEIGHTS, allow_pickle=True) as w:
        chromosomes = np.array([x.split(":")[0] for x in w["prof_region_key"].astype(str)])
    groups = sorted(set(chromosomes))
    truth = y[a]-y[b]
    ref = np.column_stack([np.sum(truth[:, chromosomes == c]**2, axis=1) for c in groups])
    losses, rows = {}, []
    for model, pred in predictions.items():
        delta = pred[a]-pred[b]
        loss = np.column_stack([np.sum((truth[:, chromosomes == c]-delta[:, chromosomes == c])**2, axis=1)
                                for c in groups])
        losses[model] = loss
        for i in range(len(pairs)):
            for j, chromosome in enumerate(groups):
                rows.append({"pair_id": i, "chromosome": chromosome, "model": model,
                             "reference_SSE": ref[i, j], "model_SSE": loss[i, j]})
    pd.DataFrame(rows).to_csv(args.out/"pair_chromosome_errors.tsv", sep="\t", index=False)
    rng = np.random.default_rng(SEED)
    sampled = rng.integers(0, len(pairs), size=(10000, len(pairs)))
    sampled_chr = rng.integers(0, len(groups), size=(10000, len(groups)))
    report = {"status": "within-histology inspected development sensitivity",
              "n_pairs": len(pairs), "n_participants": 44, "regions": truth.shape[1],
              "matching": "identical four histology grades, sex and held-participant fold",
              "skill": {m: summary(v, ref, sampled) for m, v in losses.items()},
              "crossed_pair_chromosome_sensitivity": {m: summary(v, ref, sampled, sampled_chr) for m, v in losses.items()},
              "contrasts": {}, "seed": SEED, "bootstrap_draws": len(sampled),
              "limits": "Disjoint donor pairs; intervals conditional on fitted models and selected pairs. Shared training folds omit retraining uncertainty. Chromosomes are dependence sensitivity, not unseen-region evaluation."}
    for comparator in ("technical_nuisance", "global_RRR"):
        diff = losses[comparator].sum(1)-losses["RRR_local"].sum(1)
        denom = ref.sum(1)
        observed = float(diff.sum()/denom.sum())
        bootstrap = diff[sampled].sum(1)/denom[sampled].sum(1)
        p = float((1+np.sum(np.abs(bootstrap-bootstrap.mean()) >= abs(observed)))/(len(bootstrap)+1))
        report["contrasts"]["RRR_local_vs_"+comparator] = {"delta_skill": observed,
            "ci95": np.quantile(bootstrap, [.025, .975]).tolist(), "p_nominal_centered_bootstrap": p}
    corrected = holm([v["p_nominal_centered_bootstrap"] for v in report["contrasts"].values()])
    for entry, value in zip(report["contrasts"].values(), corrected):
        entry["p_holm_two_contrasts"] = value
    (args.out/"results.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--conditional", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    main(p.parse_args())

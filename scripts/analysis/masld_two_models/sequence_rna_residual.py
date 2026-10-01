#!/usr/bin/env python3
"""One bounded low-rank sequence-by-RNA residual experiment on paired liver data.

The task is held participants at 1,024 regions fixed before this run. Unseen
region transport is a separate endpoint and is not claimed from this fit.
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
from sklearn.utils.extmath import randomized_svd

from chromatin_acquisition import fit_profile, logcpm, concentration

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT/"Analysis/MASLD_Model_Benchmark"
FIX = BENCH/"executions/model-data-064-21079902/fixture"
PAIRED = ROOT/"GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model/donor_curve_fixtures_21784771/size_8192"
WEIGHTS = BENCH/"release/masld-liver-chromatin-state-v1.2/weights/chromatin_state_v1_2.npz"
SEED = 20260922
RANKS = (4, 8, 16)
PENALTIES = (.01, .1, 1.)


def context(rna, sequence, donors, regions, query_donors, query_regions):
    """PCA transforms see only fit donors and fit regions, never test H3."""
    mu, sd = rna[donors].mean(0), rna[donors].std(0)
    sd[sd < 1e-8] = 1
    ztr = (rna[donors]-mu)/sd
    _, _, vt = randomized_svd(ztr, n_components=16, n_iter=4, random_state=SEED)
    xd, xq = ztr@vt.T, ((rna[query_donors]-mu)/sd)@vt.T
    pc_sd = xd.std(0)
    pc_sd[pc_sd < 1e-8] = 1
    xd, xq = xd/pc_sd, xq/pc_sd
    smu, ssd = sequence[regions].mean(0), sequence[regions].std(0)
    ssd[ssd < 1e-8] = 1
    sz = (sequence[regions]-smu)/ssd
    _, _, svt = randomized_svd(sz, n_components=16, n_iter=4, random_state=SEED+1)
    sr, sq = sz@svt.T, ((sequence[query_regions]-smu)/ssd)@svt.T
    sr_sd = sr.std(0)
    sr_sd[sr_sd < 1e-8] = 1
    return xd, xq, sr/sr_sd, sq/sr_sd


def fit_predict(rna, sequence, residual, donors, query_donors, regions, query_regions, penalty):
    xd, xq, sr, sq = context(rna, sequence, donors, regions, query_donors, query_regions)
    r = residual[np.ix_(donors, regions)]
    region_mean = r.mean(0)
    # RNA-only is an independent per-region linear head, which is a demanding
    # comparator for the compact sequence-conditioned interaction.
    kr = xd.T@xd
    lam_rna = .1*max(float(np.linalg.eigvalsh(kr).max()), 1e-12)
    rna_coef = np.linalg.solve(kr+lam_rna*np.eye(len(kr)), xd.T@(r-region_mean))
    rna_pred = region_mean+(xq@rna_coef) if np.array_equal(regions, query_regions) else None
    donor_mean = r.mean(1)-r.mean()
    donor_coef = np.linalg.solve(kr+lam_rna*np.eye(len(kr)), xd.T@donor_mean)
    ks = sr.T@sr
    lam_seq = .1*max(float(np.linalg.eigvalsh(ks).max()), 1e-12)
    seq_coef = np.linalg.solve(ks+lam_seq*np.eye(len(ks)), sr.T@(region_mean-r.mean()))
    sequence_train = r.mean()+sr@seq_coef
    sequence_mean = r.mean()+sq@seq_coef
    additive_train = sequence_train[None, :]+(xd@donor_coef)[:, None]
    additive_query = sequence_mean[None, :]+(xq@donor_coef)[:, None]
    interaction_target = r-additive_train
    dx, ux = np.linalg.eigh(kr)
    ds, us = np.linalg.eigh(sr.T@sr)
    scale = max(float(dx.max()*ds.max()), 1e-12)
    cross = ux.T@(xd.T@interaction_target@sr)@us
    coefficient = ux@(cross/(dx[:, None]*ds[None, :]+penalty*scale))@us.T
    u, singular, vt = np.linalg.svd(coefficient, full_matrices=False)
    result = {"sequence_only": np.broadcast_to(sequence_mean, (len(query_donors), len(query_regions))).copy(),
              "additive": additive_query,
              "rna_only": rna_pred}
    for rank in RANKS:
        w = (u[:, :rank]*singular[:rank])@vt[:rank]
        result[f"interaction_r{rank}"] = additive_query+xq@w@sq.T
    return result


def choose(rna, sequence, residual, donors, regions):
    rng = np.random.default_rng(SEED+len(donors))
    parts = np.array_split(rng.permutation(donors), 3)
    errors = {(rank, lam): 0.0 for rank in RANKS for lam in PENALTIES}
    for held in parts:
        trained = np.setdiff1d(donors, held)
        target = residual[np.ix_(held, regions)]
        for lam in PENALTIES:
            candidate = fit_predict(rna, sequence, residual, trained, held, regions, regions, lam)
            for rank in RANKS:
                errors[rank, lam] += float(np.sum((target-candidate[f"interaction_r{rank}"])**2))
    selected = min(errors, key=lambda pair: (errors[pair], pair[0], pair[1]))
    return selected, {f"r{rank}_lambda{lam}": error for (rank, lam), error in errors.items()}


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.mkdir(parents=True)
    reg = pd.read_csv(PAIRED/"regions.tsv", sep="\t")
    subset = reg.loc[reg.training_1024 & reg.region_role.eq("train")]
    if len(subset) != 1024 or subset.region_index.duplicated().any():
        raise ValueError("Training-region annotation or genomic components changed")
    with np.load(PAIRED/"paired.npz", allow_pickle=True) as f:
        ids = f["participant_ids"].astype(str)
        donor_fold = f["donor_fold"].astype(int)
    with np.load(PAIRED/"features.npz", allow_pickle=False) as f:
        if not np.array_equal(f["region_key"].astype(str), reg.region_key.astype(str).to_numpy()):
            raise ValueError("Frozen sequence feature identity differs")
        sequence = np.asarray(f["mean_strands"], np.float32)
    if len(ids) != 99 or len(set(ids)) != 99:
        raise ValueError("Expected 99 participants")
    selected = subset.region_index.to_numpy(int)
    sequence_selected = sequence[subset.index.to_numpy()]
    with np.load(WEIGHTS, allow_pickle=True) as f:
        cis = np.asarray(f["prof_Ocis_cis_idx"], int)[selected]
        region_axis = f["prof_region_key"].astype(str)
    if not np.array_equal(region_axis[selected], subset.region_key.astype(str).to_numpy()):
        raise ValueError("Subset region keys differ from released annotation")
    rna = np.asarray(np.load(FIX/"molecular/rna_values.npy", mmap_mode="r"), float)
    h3 = np.asarray(np.load(FIX/"molecular/h3k27ac_counts.npy", mmap_mode="r"), float)
    x, yraw, conc = logcpm(rna), logcpm(h3), concentration(h3)
    rows, choices = [], []
    for outer in sorted(set(donor_fold)):
        tr, te = np.flatnonzero(donor_fold != outer), np.flatnonzero(donor_fold == outer)
        design = np.column_stack([np.ones(len(tr)), conc[tr]])
        coeff = np.linalg.lstsq(design, yraw[tr], rcond=None)[0]
        y = yraw-np.column_stack([np.ones(len(ids)), conc])@coeff
        # The strong RNA baseline's training residuals are generated without
        # each recipient donor's chromatin, inside the outer training fold.
        residual = np.zeros((len(ids), len(selected)), np.float32)
        parts = np.array_split(np.random.default_rng(SEED+outer).permutation(tr), 3)
        for part in parts:
            internal = np.setdiff1d(tr, part)
            _, estimated = fit_profile(x, y, internal, part, cis, selected)
            residual[part] = y[part][:, selected]-estimated
        _, base_test = fit_profile(x, y, tr, te, cis, selected)
        residual[te] = y[te][:, selected]-base_test
        (rank, penalty), inner = choose(x, sequence_selected, residual, tr, np.arange(len(selected)))
        choices.append({"outer_fold": int(outer), "selected_rank": rank,
                        "selected_penalty": penalty, "inner_errors": inner})
        fitted = fit_predict(x, sequence_selected, residual, tr, te,
                             np.arange(len(selected)), np.arange(len(selected)), penalty)
        # Correct-context and shuffled-context use the same selected rank and
        # penalty. Shuffle within recorded fibrosis stage among training donors.
        join = pd.read_csv(BENCH/"executions/gse267145-authoritative-join-21064930/participant_join.tsv",
                           sep="\t").set_index("participant_id").loc[ids]
        fibrosis = join.fibrosis.to_numpy()
        shuffled = x.copy()
        rng = np.random.default_rng(SEED+outer+100)
        for stage in np.unique(fibrosis[tr]):
            group = tr[fibrosis[tr] == stage]
            if len(group) > 1:
                shuffled[group] = x[rng.permutation(group)]
        null = fit_predict(shuffled, sequence_selected, residual, tr, te,
                           np.arange(len(selected)), np.arange(len(selected)), penalty)
        for j, donor in enumerate(te):
            base_sse = float(np.sum(residual[donor]**2))
            record = {"participant_id": ids[donor], "fold": int(outer),
                      "baseline_SSE": base_sse,
                      "rna_only_SSE": float(np.sum((residual[donor]-fitted["rna_only"][j])**2)),
                      "sequence_only_SSE": float(np.sum((residual[donor]-fitted["sequence_only"][j])**2)),
                      "additive_SSE": float(np.sum((residual[donor]-fitted["additive"][j])**2)),
                      "interaction_SSE": float(np.sum((residual[donor]-fitted[f"interaction_r{rank}"][j])**2)),
                      "shuffled_interaction_SSE": float(np.sum((residual[donor]-null[f"interaction_r{rank}"][j])**2))}
            rows.append(record)
        print(json.dumps({"completed_fold": int(outer), "selected_rank": rank,
                          "selected_penalty": penalty}), flush=True)
    table = pd.DataFrame(rows)
    if table.participant_id.duplicated().any() or len(table) != 99:
        raise ValueError("Participant evaluation census differs")
    table.to_csv(args.out/"held_participant_errors.tsv", sep="\t", index=False)
    denom = table.baseline_SSE.sum()
    skills = {name: float(1-table[name+"_SSE"].sum()/denom)
              for name in ("rna_only", "sequence_only", "additive", "interaction", "shuffled_interaction")}
    best_simple = max(skills[name] for name in ("rna_only", "sequence_only", "additive"))
    report = {"status": "single-cohort inspected development; one added model experiment",
              "biological_n": 99, "train_regions": len(selected), "unseen_region_claim": False,
              "target": "H3K27ac residual after cross-fitted rank-12 RRR plus local gene correction",
              "selected_inside_outer_training": choices, "held_participant_residual_skill": skills,
              "interaction_increment_over_best_simple": skills["interaction"]-best_simple,
              "development_complexity_margin_met": skills["interaction"]-best_simple >= .01,
              "shuffled_context_control": "training RNA permuted within fibrosis stage; evaluated with real held RNA",
              "uncertainty": "paired participant resampling pending; no superiority claim from point estimates alone",
              "environment": {"python": sys.version, "numpy": np.__version__,
                              "pandas": pd.__version__, "platform": platform.platform(),
                              "slurm_job_id": os.environ["SLURM_JOB_ID"]}}
    (args.out/"results.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"skills": skills, "increment": report["interaction_increment_over_best_simple"]}, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())

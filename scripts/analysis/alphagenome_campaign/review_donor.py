#!/usr/bin/env python3
"""Bounded GSE267145 donor-difference development comparison.

All fitting occurs in the training donors and, for the crossed evaluation,
training tiles. The target is log2(1+CPM) H3K27ac, not an integer RNA count.
RNA inputs remain continuous. This is a simple CPU feasibility comparison,
not a re-fit of Corgi or evidence of generalization to another cohort.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT / "Analysis/MASLD_Model_Benchmark"
FIX = BENCH / "executions/model-data-064-21079902/fixture"
ZS = BENCH / "executions/corgi-gse267145-h3k27ac-comparator-20260907T185333Z/out"
GTF = Path("/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz")
SEED = 20260915
ARMS = ("all_rna", "overlapping_genes_excluded", "within_1mb_excluded",
        "chromosome_excluded", "random_overlap", "random_1mb", "random_chromosome",
        "technical_only", "context_shuffled")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def stable_int(value):
    return int(hashlib.sha256(value.encode()).hexdigest()[:15], 16)


def logcpm(values):
    values = np.asarray(values, float)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("nonnegative finite continuous input required")
    total = values.sum(1, keepdims=True)
    if np.any(total <= 0):
        raise ValueError("empty library after exclusions")
    return np.log2(1 + values / total * 1e6)


def standardize(values, train):
    mean = values[train].mean(0)
    sd = values[train].std(0)
    return (values - mean) / np.where(sd > 1e-10, sd, 1)


def partners(n, epoch, seed=SEED):
    """Two distinct partners, exactly two appearances in each role per donor."""
    if n < 3:
        raise ValueError("three donors required")
    rng = np.random.default_rng(seed + epoch)
    ring = rng.permutation(n)
    shifts = rng.choice(np.arange(1, n), 2, replace=False)
    return np.tile(ring, 2), np.concatenate([np.roll(ring, int(k)) for k in shifts])


def shuffle_by_stratum(ids, strata, seed):
    """Shuffle only within this donor partition and observed source/batch strata."""
    ids = np.asarray(ids)
    output = ids.copy()
    rng = np.random.default_rng(seed)
    for group in sorted(set(strata[ids])):
        pos = np.flatnonzero(strata[ids] == group)
        if len(pos) > 1:
            ring = rng.permutation(pos)
            output[ring] = ids[np.roll(ring, 1)]
    return output


def ridge(x, y, alpha=10):
    penalty = np.eye(x.shape[1]) * alpha
    penalty[0, 0] = 0
    return np.linalg.solve(x.T @ x + penalty, x.T @ y)


def sequence_features(seq):
    a = np.frombuffer(seq.upper().encode(), dtype="S1")
    if not len(a) or np.mean(np.isin(a, list(map(str.encode, "ACGT")))) < 0.99:
        raise ValueError("ambiguous counted interval sequence")
    # Reverse-complement invariant dinucleotide composition is an explicit,
    # inexpensive sequence-only comparator, not an AlphaGenome embedding.
    rc = seq.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]
    kmers = [a + b for a in "ACGT" for b in "ACGT"]
    f = [(sum(s[i:i+2] == k for s in (seq.upper(), rc)
              for i in range(len(s)-1))) / (2*max(len(seq)-1, 1)) for k in kmers]
    return f + [np.log1p(len(seq))]


def load_genes(gene_ids):
    wanted = set(gene_ids)
    found = {}
    with gzip.open(GTF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if f[2] != "gene":
                continue
            match = re.search(r'gene_id "([^".]+)', f[8])
            if match and match.group(1) in wanted:
                # RNA is mapped to GENCODE49 primary assembly. Retain the
                # primary mapping; duplicated haplotype annotations do not
                # override a primary-contig gene.
                key = match.group(1)
                if f[0] in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}:
                    found[key] = (f[0], int(f[3])-1, int(f[4]))
    return [found.get(g, ("unknown", -1, -1)) for g in gene_ids]


def load_substrate(regions_per_tile):
    parts = pd.read_csv(FIX / "molecular/participant_axis.tsv", sep="\t")
    ids = parts.participant_id.astype(str).to_numpy()
    if len(ids) != 99 or len(set(ids)) != 99:
        raise ValueError("expected 99 independent participants")
    foldtab = pd.read_csv(FIX / "folds/participant_outer_folds.tsv", sep="\t")
    fold = foldtab.set_index("participant_id").loc[ids, "outer_fold"].to_numpy(int)
    genes = pd.read_csv(FIX / "molecular/rna_feature_axis.tsv", sep="\t")
    rna = np.load(FIX / "molecular/rna_values.npy")
    h3 = np.load(FIX / "molecular/h3k27ac_counts.npy")
    if rna.shape != (99, 42163) or h3.shape != (99, 96460):
        raise ValueError("fixture dimensions changed")
    if not np.issubdtype(h3.dtype, np.integer):
        raise ValueError("H3 labels must retain source integer counts")
    regions = pd.read_csv(ZS / "tile_regions.tsv", sep="\t")
    tiles = pd.read_csv(ZS / "tiles.tsv", sep="\t")
    sequences, names = [], []
    with gzip.open(ZS / "tiles.fa.gz", "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                names.append(line[1:].strip())
                sequences.append("")
            else:
                sequences[-1] += line.strip()
    if len(sequences) != 64 or len(tiles) != 64:
        raise ValueError("established 64-tile substrate changed")
    hashes = [int(hashlib.sha256(f"corgi-film-tile|{n}".encode()).hexdigest()[:16], 16) for n in names]
    train_tiles = set(np.argsort(np.asarray(hashes, dtype=np.uint64))[:32])
    # Sampling uses only interval identity, never activity or predictability.
    regions["sample_hash"] = regions.region_key.map(lambda k: stable_int("week1-b2|" + k))
    reg = regions.sort_values("sample_hash").groupby("tile_index", sort=True).head(regions_per_tile)
    reg = reg.sort_values(["tile_index", "region_index"]).reset_index(drop=True)
    sf = []
    for r in reg.itertuples():
        start = int(tiles.iloc[r.tile_index].window_start)
        subseq = sequences[r.tile_index][int(r.start0)-start:int(r.end0)-start]
        if len(subseq) != int(r.end0-r.start0):
            raise ValueError("counted interval does not map into tile")
        sf.append(sequence_features(subseq))
    reg["region_role"] = ["train" if t in train_tiles else "held" for t in reg.tile_index]
    # Effective sequence inputs are counted intervals. Existing 524kb Corgi
    # windows may overlap, but are not inputs to this CPU model.
    overlaps = 0
    for a in reg[reg.region_role == "train"].itertuples():
        for b in reg[reg.region_role == "held"].itertuples():
            overlaps += int(a.chrom == b.chrom and a.start0 < b.end0 and b.start0 < a.end0)
    if overlaps:
        raise ValueError("train/held counted intervals overlap")
    coords = load_genes(genes.stable_gene_id.astype(str).tolist())
    gc, gs, ge = map(np.asarray, zip(*coords))
    return ids, fold, rna, logcpm(h3)[:, reg.region_index.to_numpy()], reg, np.asarray(sf), (gc, gs, ge)


def exclusion_masks(reg, gene_coords):
    chrom, start, end = gene_coords
    unknown = chrom == "unknown"
    masks = {}
    rows = []
    for j, r in enumerate(reg.itertuples()):
        overlap = (chrom == r.chrom) & (start < r.end0) & (end > r.start0)
        nearby = (chrom == r.chrom) & (start < r.end0+1_000_000) & (end > r.start0-1_000_000)
        samechr = chrom == r.chrom
        per = {"all_rna": np.zeros(len(chrom), bool), "overlapping_genes_excluded": overlap,
               "within_1mb_excluded": nearby, "chromosome_excluded": samechr,
               "technical_only": np.zeros(len(chrom), bool), "context_shuffled": np.zeros(len(chrom), bool)}
        for control, mask in (("random_overlap", overlap), ("random_1mb", nearby), ("random_chromosome", samechr)):
            rng = np.random.default_rng(SEED + stable_int(r.region_key + control))
            random = np.zeros(len(chrom), bool)
            random[rng.choice(np.flatnonzero(~unknown), int(mask.sum()), replace=False)] = True
            per[control] = random
        for arm in ARMS:
            # Unknown-coordinate genes excluded in every arm, including
            # matched controls, so no unknown local gene escapes a mask.
            masks[j, arm] = ~(per[arm] | unknown)
            rows.append(dict(region_key=r.region_key, arm=arm, excluded=int(per[arm].sum()),
                             missing_coordinates=int(unknown.sum()), target_link_status="unknown",
                             target_exclusion_status="untestable_no_validated_target_links",
                             coordinate_semantics="counted_interval", coordinate_residual_bp=1))
    return masks, rows


def contexts(rna, mask, train, arm, seed, rank=8):
    if arm == "technical_only":
        lib = rna.sum(1)
        x = np.column_stack([np.log1p(lib), np.mean(rna > 0, axis=1),
                             np.sort(rna, axis=1)[:, -50:].sum(1)/lib])
        return standardize(x, train)
    raw = rna[:, mask]
    x = logcpm(raw)  # excluded values cannot enter normalization
    var = x[train].var(0)
    chosen = np.lexsort((np.arange(len(var)), -var))[:min(256, len(var))]
    chosen = chosen[var[chosen] > 1e-10]
    z = standardize(x[:, chosen], train)
    # A fixed gene-identity projection keeps the axes aligned between region
    # masks. Its coefficients and rank use no chromatin outcome.
    indices = np.flatnonzero(mask)[chosen]
    proj = np.stack([np.random.default_rng(SEED + int(i)).choice([-1., 1.], rank) for i in indices])
    projected = z @ proj / np.sqrt(len(chosen))
    return standardize(projected, train)


def design(x, s, interactions):
    n, nr, _ = x.shape
    tiled = np.broadcast_to(s[None], (n, nr, s.shape[1]))
    pieces = [np.ones((n, nr, 1)), tiled, x]
    if interactions:
        pieces.append((x[..., :, None] * tiled[..., None, :]).reshape(n, nr, -1))
    return np.concatenate(pieces, axis=2)


def fit_shared(x, s, y, donors, regions, objective, epochs):
    d = design(x, s, True)
    t = d[np.ix_(donors, regions)]
    yt = y[np.ix_(donors, regions)]
    if objective == "profile":
        beta = ridge(t.reshape(-1, t.shape[-1]), yt.reshape(-1))
        return np.einsum("drp,p->dr", d, beta)
    seq = np.column_stack([np.ones(len(s)), s])
    mu = yt.mean(0)
    base = mu if len(regions) == len(s) else seq @ ridge(seq[regions], mu)
    # Training-region centering never reads measured held-region means.
    residual = yt - mu
    tx = t - t.mean(0, keepdims=True)
    if objective == "training_mean_residual":
        xx = tx.reshape(-1, t.shape[-1])
        yy = residual.reshape(-1)
        penalty = np.eye(xx.shape[1]) * 10
        beta = np.linalg.solve(xx.T @ xx + penalty, xx.T @ yy)
    elif objective == "pairwise":
        gram = np.eye(t.shape[-1]) * 10
        rhs = np.zeros(t.shape[-1])
        for epoch in range(epochs):
            a, b = partners(len(donors), epoch)
            dx = (t[a] - t[b]).reshape(-1, t.shape[-1])
            dy = (yt[a] - yt[b]).reshape(-1)
            # Average over epochs preserves the fixed penalty scale.
            gram += dx.T @ dx / epochs
            rhs += dx.T @ dy / epochs
        beta = np.linalg.solve(gram, rhs)
    else:
        raise ValueError(objective)
    dc = d - d[donors].mean(0, keepdims=True)
    return base[None, :] + np.einsum("drp,p->dr", dc, beta)


def donor_difference_loss(errors, folds):
    """Mean squared difference across all distinct held pairs, within fold.

    The identity E[(e_i-e_j)^2] = 2 sample_variance(e) avoids enumerating
    dependent pairs. Fold weighting gives every index donor equal weight.
    """
    result = np.zeros(errors.shape[1])
    for fold in sorted(set(folds)):
        donors = np.flatnonzero(folds == fold)
        result += len(donors)/len(folds) * 2*errors[donors].var(0, ddof=1)
    return result


def metrics(pred, target, baseline, folds):
    residual = target-baseline
    error = (pred-target)**2
    denom = np.mean(residual**2)
    rhos = [spearmanr(pred[:, j], target[:, j]).statistic
            for j in range(target.shape[1]) if np.ptp(pred[:, j]) > 1e-12 and np.ptp(target[:, j]) > 1e-12]
    derror = donor_difference_loss(pred-target, folds).mean()
    return dict(mse_log2cpm=float(error.mean()), skill_vs_training_baseline=float(1-error.mean()/denom),
                donor_difference_mse_log2cpm=float(derror),
                donor_difference_skill_vs_zero=float(1-derror/donor_difference_loss(target, folds).mean()),
                mean_region_donor_spearman=float(np.mean(rhos)) if rhos else None,
                rank_defined_regions=len(rhos), n_donors=target.shape[0], n_regions=target.shape[1])


def uncertainty(diff, region_blocks, draws=1000):
    """Paired two-way donor-by-tile resampling; seeds are not replication."""
    rng = np.random.default_rng(SEED)
    blocks = np.unique(region_blocks)
    by = {b: np.flatnonzero(region_blocks == b) for b in blocks}
    values = []
    for _ in range(draws):
        di = rng.integers(0, len(diff), len(diff))
        ri = np.concatenate([by[b] for b in rng.choice(blocks, len(blocks), replace=True)])
        values.append(diff[np.ix_(di, ri)].mean())
    return np.quantile(values, [0.025, 0.975]).tolist()


def difference_uncertainty(pred, target, folds, region_blocks, draws=1000):
    """Resample original donors, never treat their paired contrasts as n."""
    rng = np.random.default_rng(SEED)
    blocks = np.unique(region_blocks)
    by = {b: np.flatnonzero(region_blocks == b) for b in blocks}
    fold_donors = [np.flatnonzero(folds == f) for f in sorted(set(folds))]
    values = []
    for _ in range(draws):
        di = np.concatenate([rng.choice(d, len(d), replace=True) for d in fold_donors])
        ri = np.concatenate([by[b] for b in rng.choice(blocks, len(blocks), replace=True)])
        e = (pred-target)[np.ix_(di, ri)]
        observed = target[np.ix_(di, ri)]
        delta = donor_difference_loss(e, folds[di])-donor_difference_loss(observed, folds[di])
        values.append(delta.mean())
    return np.quantile(values, [0.025, 0.975]).tolist()


def self_test():
    for n in (3, 8, 21):
        a, b = partners(n, 3)
        assert np.all(a != b)
        assert np.all(np.bincount(a) == 2) and np.all(np.bincount(b) == 2)
    x = np.arange(24, dtype=float).reshape(8, 3)
    train = np.arange(5)
    z = standardize(x, train)
    q = x.copy(); q[5:] += 1e6
    assert np.array_equal(z[train], standardize(q, train)[train])
    strata = np.array(["a"]*4 + ["b"]*4)
    sh = shuffle_by_stratum(np.arange(8), strata, 7)
    assert np.all(sh != np.arange(8)) and np.all(strata == strata[sh])
    rng = np.random.default_rng(SEED)
    xx = rng.normal(size=(10, 5, 2)); ss = rng.normal(size=(5, 2)); yy = rng.normal(size=(10, 5))
    tr = np.arange(7); regions = np.arange(3)
    for objective in ("profile", "training_mean_residual", "pairwise"):
        p = fit_shared(xx, ss, yy, tr, regions, objective, 2)
        changed = yy.copy(); changed[7:] += 1000; changed[:, 3:] -= 2000
        p2 = fit_shared(xx, ss, changed, tr, regions, objective, 2)
        assert np.allclose(p, p2, atol=1e-10), objective
    print(json.dumps({"tests": "pass", "checks": ["two_balanced_partners", "training_only_scaling",
          "within_stratum_shuffle", "held_donor_and_region_outcome_invariance"]}), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--regions-per-tile", type=int, default=2)
    ap.add_argument("--pair-epochs", type=int, default=20)
    ap.add_argument("--self-test-only", action="store_true")
    args = ap.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("run executable scientific checks on a compute node")
    self_test()
    if args.self_test_only:
        return
    args.output.mkdir(parents=True, exist_ok=False)
    t0 = time.monotonic()
    ids, fold, rna, y, reg, seq, coords = load_substrate(args.regions_per_tile)
    masks, coverage = exclusion_masks(reg, coords)
    pd.DataFrame(coverage).to_csv(args.output / "exclusion_coverage.tsv", sep="\t", index=False)
    reg.to_csv(args.output / "regions.tsv", sep="\t", index=False)
    metrics_rows, contrasts, all_predictions = [], [], {}
    # GEO supplies one source and no authoritative RNA batch labels in the
    # paired fixture. Shuffle within source and donor partition; batch-specific
    # shuffling cannot be claimed. This limit is recorded explicitly.
    strata = np.array(["GSE267145"] * len(ids))
    for mode in ("held_donors_trained_regions", "held_donors_held_regions"):
        region_train = np.arange(len(reg)) if mode.endswith("trained_regions") else np.flatnonzero(reg.region_role == "train")
        region_test = np.arange(len(reg)) if mode.endswith("trained_regions") else np.flatnonzero(reg.region_role == "held")
        mode_pred = {}
        baseline = np.zeros_like(y)
        for f in sorted(set(fold)):
            tr = np.flatnonzero(fold != f); te = np.flatnonzero(fold == f)
            s = standardize(seq, region_train)
            sx = np.column_stack([np.ones(len(s)), s])
            mu = y[tr].mean(0)
            seqpred = sx @ ridge(sx[region_train], mu[region_train])
            base = mu if mode.endswith("trained_regions") else seqpred
            baseline[te] = base
            for arm in ARMS:
                x = np.stack([contexts(rna, masks[j, arm], tr, arm, SEED+f)
                              for j in range(len(reg))], axis=1)
                if arm == "context_shuffled":
                    order = np.arange(len(ids))
                    order[tr] = shuffle_by_stratum(tr, strata, SEED+f)
                    order[te] = shuffle_by_stratum(te, strata, SEED+100+f)
                    x = x[order]
                for objective in ("profile", "training_mean_residual", "pairwise"):
                    pred = fit_shared(x, s, y, tr, region_train, objective, args.pair_epochs)
                    key = arm + ":shared_" + objective
                    mode_pred.setdefault(key, np.zeros_like(y))[te] = pred[te]
                d = design(x, s, False)
                beta = ridge(d[np.ix_(tr, region_train)].reshape(-1, d.shape[-1]), y[np.ix_(tr, region_train)].reshape(-1))
                mode_pred.setdefault(arm+":additive", np.zeros_like(y))[te] = np.einsum("drp,p->dr", d[te], beta)
                if mode.endswith("trained_regions"):
                    pred = np.zeros_like(y)
                    for j in range(len(reg)):
                        xr = np.column_stack([np.ones(len(ids)), x[:, j]])
                        pred[:, j] = xr @ ridge(xr[tr], y[tr, j])
                    mode_pred.setdefault(arm+":rna_only_region_specific", np.zeros_like(y))[te] = pred[te]
                print(f"{mode} fold={f} arm={arm} elapsed={time.monotonic()-t0:.1f}s", flush=True)
        mode_pred["sequence_only"] = np.zeros_like(y)
        # For trained regions, sequence-only has no participant information;
        # the stronger training-region mean is also retained as comparator.
        for f in sorted(set(fold)):
            tr = np.flatnonzero(fold != f); te = np.flatnonzero(fold == f)
            s = standardize(seq, region_train); sx = np.column_stack([np.ones(len(s)), s])
            mode_pred["sequence_only"][te] = sx @ ridge(sx[region_train], y[np.ix_(tr, region_train)].mean(0))
        mode_pred["training_baseline"] = baseline
        for key, pred in mode_pred.items():
            m = metrics(pred[:, region_test], y[:, region_test], baseline[:, region_test], fold)
            metrics_rows.append(dict(evaluation=mode, model=key, **m))
            if key != "training_baseline":
                # Negative difference favors the candidate. No threshold or
                # recipe selection is performed using these development CIs.
                delta = (pred[:, region_test]-y[:, region_test])**2 - (baseline[:, region_test]-y[:, region_test])**2
                lo, hi = uncertainty(delta, reg.tile_index.to_numpy()[region_test])
                contrasts.append(dict(evaluation=mode, model=key, comparator="training_baseline",
                                      endpoint="profile_mse_log2cpm",
                                      mse_difference=float(delta.mean()), ci95_low=lo, ci95_high=hi,
                                      family="descriptive_pilot_no_significance_claim"))
                dd = donor_difference_loss(pred[:, region_test]-y[:, region_test], fold) - donor_difference_loss(y[:, region_test], fold)
                dlo, dhi = difference_uncertainty(pred[:, region_test], y[:, region_test], fold, reg.tile_index.to_numpy()[region_test])
                contrasts.append(dict(evaluation=mode, model=key, comparator="zero_donor_difference",
                                      endpoint="donor_difference_mse_log2cpm",
                                      mse_difference=float(dd.mean()), ci95_low=dlo, ci95_high=dhi,
                                      family="descriptive_pilot_no_significance_claim"))
            all_predictions[mode+"|"+key] = pred[:, region_test].astype(np.float32)
    pd.DataFrame(metrics_rows).to_csv(args.output/"metrics.tsv", sep="\t", index=False)
    pd.DataFrame(contrasts).to_csv(args.output/"paired_uncertainty.tsv", sep="\t", index=False)
    # Predictions remain internal research files; no source measurements are
    # copied into this output or the public catalog.
    np.savez_compressed(args.output/"internal_predictions.npz", **all_predictions)
    inputs = [FIX/"molecular/rna_values.npy", FIX/"molecular/h3k27ac_counts.npy",
              FIX/"folds/participant_outer_folds.tsv", ZS/"tile_regions.tsv", ZS/"tiles.tsv", GTF]
    receipt = dict(status="development_pilot_not_independent_confirmation", seed=SEED,
                   biological_n=99, regions=len(reg), tile_n=int(reg.tile_index.nunique()),
                   donor_fold_counts={str(f): int(sum(fold == f)) for f in set(fold)},
                   region_split="established_sha256_corgi_film_32_train_32_held_tiles",
                   sequence_inputs="reverse_complement_invariant_counted_interval_dinucleotides_and_width",
                   counted_interval_uncertainty_bp=1, actual_input_train_held_overlap=0,
                   source_strata="one_GSE267145_source", batch_stratified_shuffle="untestable_no_batch_authority",
                   target_exclusion="untestable_no_validated_target_links; overlapping_gene_sensitivity_only",
                   objective_units="log2(1+H3K27ac_CPM)", rna_units="continuous_fractional_estimates",
                   calibration="direct_native_scale_training_fit; no_held_outcome_fitted_scaling",
                   region_specific_rna_comparator_crossed="ineligible_not_run",
                   covariates="none; separate_RNA_technical_only_control; no_histology_conditioning",
                   inference_unit="donor_and_tile; paired_two_way_bootstrap_1000_draws",
                   correction="none_descriptive_pilot; no_confirmatory_tests_or_claims",
                   protected_outcomes_read=False, seconds=time.monotonic()-t0,
                   software={"python":sys.version,"numpy":np.__version__,"pandas":pd.__version__},
                   input_sha256={str(p):sha(p) for p in inputs})
    (args.output/"receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()

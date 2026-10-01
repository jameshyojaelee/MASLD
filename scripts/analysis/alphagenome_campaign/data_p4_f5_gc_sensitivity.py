#!/usr/bin/env python3
"""F5g: repeat the F5 program TF profiles against a GC-matched random gene-set null.

F5's prespecified null matches on gene-set size alone, while F1-F4 match on
GC x width x signal-mean x signal-sd x promoter strata. The unstratified F5 hits
depleted promoter-generic factors and enriched liver-identity factors, which is
what a promoter CpG/GC composition difference produces, so the unstratified
family cannot separate predicted TF sensitivity from promoter base composition.

This script reuses F5's own region universe, promoter links, TF matrix, draw
count and statistic unchanged, and changes exactly one thing: each draw now
samples genes within promoter-GC deciles so the drawn set carries the program's
GC composition. Predictions are in data_p4_f5_gc_prespec.json, written before
any stratified score was computed. Nothing here re-opens the region universe and
no model score enters the stratification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control

import data_p4_families as fam
import data_p4_f5_f6_families as f5

PRESPEC = Path(__file__).with_name("data_p4_f5_gc_prespec.json")
N_DECILES = 10


def region_gc(regions: list[str], frame: pd.DataFrame) -> np.ndarray:
    """GC per region from the reference FASTA, on the keyed interval as written."""
    import pysam
    coords = frame.drop_duplicates("region_key").set_index("region_key")
    genome = pysam.FastaFile(f5.FASTA)
    try:
        out = np.empty(len(regions))
        for i, key in enumerate(regions):
            row = coords.loc[key]
            seq = genome.fetch(str(row.chrom), int(row.start0), int(row.end)).upper()
            out[i] = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
    finally:
        genome.close()
    return out


def profile_and_gc(gene_regions: sparse.csr_matrix, region_tf: sparse.csr_matrix,
                   gc: np.ndarray, rows: np.ndarray):
    """TF share and mean GC over the union of a gene set's promoter regions."""
    mask = np.asarray(gene_regions[rows].max(axis=0).todense()).ravel() > 0
    n = int(mask.sum())
    if not n:
        return None, 0, float("nan")
    share = np.asarray(mask[None, :].astype(float) @ region_tf).ravel() / n
    return share, n, float(gc[mask].mean())


def stratified_draw(rng, decile_members: list[np.ndarray], counts: np.ndarray) -> np.ndarray:
    picks = []
    for d, want in enumerate(counts):
        if want:
            picks.append(rng.choice(decile_members[d], size=want, replace=False))
    return np.concatenate(picks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reduced", type=Path, nargs="+", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)

    prespec = json.loads(PRESPEC.read_text())
    if prespec["family"] != "F5g_gc_matched_sensitivity":
        raise ValueError("F5g prespecification identity changed")
    seed, draws_n = prespec["seed"], prespec["draws"]
    if draws_n != f5.F5_DRAWS:
        raise ValueError("The sensitivity must use F5's own draw count to be comparable")

    feats, tfs = f5.read_features(list(args.reduced))
    da = f5.eligible_da()
    links = f5.promoter_links(da)
    available = set(feats.region_key)
    background_genes = sorted(links)
    regions = sorted(set().union(*links.values()) & available)
    region_index = {k: i for i, k in enumerate(regions)}
    gene_index = {g: i for i, g in enumerate(background_genes)}

    gr = sparse.lil_matrix((len(background_genes), len(regions)), dtype=bool)
    for gene, rs in links.items():
        for r in rs & available:
            gr[gene_index[gene], region_index[r]] = True
    gene_regions = gr.tocsr()

    relevant = tfs.loc[tfs.region_key.isin(region_index) & tfs.group.isin(f5.TF_GROUPS)].copy()
    relevant["track_key"] = relevant.group + "|" + relevant.track_name
    tracks = sorted(relevant.track_key.unique())
    track_index = {t: i for i, t in enumerate(tracks)}
    factor = dict(zip(relevant.track_key, relevant.transcription_factor))
    rt = sparse.lil_matrix((len(regions), len(tracks)), dtype=np.float64)
    for row in relevant.itertuples():
        rt[region_index[row.region_key], track_index[row.track_key]] = 1.0
    region_tf = rt.tocsr()

    gc = region_gc(regions, feats)
    # A gene's GC is the mean over its own promoter regions; genes are binned once,
    # over the whole eligible background universe, before any program is read.
    gene_gc = np.array([gc[gene_regions[i].indices].mean() if gene_regions[i].nnz else np.nan
                        for i in range(len(background_genes))])
    usable = ~np.isnan(gene_gc)
    edges = np.quantile(gene_gc[usable], np.linspace(0, 1, N_DECILES + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    decile = np.digitize(gene_gc, edges[1:-1], right=False)
    decile[~usable] = -1
    decile_members = [np.flatnonzero((decile == d) & usable) for d in range(N_DECILES)]

    programs = pd.read_csv(fam.PROGRAMS, sep="\t")
    rng = np.random.default_rng(seed)
    rows, matching = [], []
    for uid, part in programs.groupby("program_uid", sort=True):
        genes = sorted({g for g in part.canonical_gene.dropna() if g in gene_index})
        chosen = np.array([gene_index[g] for g in genes], dtype=int)
        chosen = chosen[usable[chosen]]
        program_regions = set().union(*(links.get(g, set()) for g in genes)) if genes else set()
        if len(program_regions & available) < fam.MIN_GROUP or len(chosen) < 2:
            continue
        counts = np.bincount(decile[chosen], minlength=N_DECILES)
        short = [(d, int(counts[d]), len(decile_members[d])) for d in range(N_DECILES)
                 if counts[d] > len(decile_members[d])]
        if short:
            matching.append({"program_uid": uid, "status": "indeterminate_decile_exhausted",
                             "detail": json.dumps(short)})
            continue
        observed, n_observed, gc_observed = profile_and_gc(gene_regions, region_tf, gc, chosen)
        strat = np.empty((draws_n, region_tf.shape[1]))
        strat_gc = np.empty(draws_n)
        unstrat_gc = np.empty(draws_n)
        rng_u = np.random.default_rng(seed + 1)
        for i in range(draws_n):
            pick = stratified_draw(rng, decile_members, counts)
            p, _, g = profile_and_gc(gene_regions, region_tf, gc, pick)
            strat[i] = p if p is not None else 0.0
            strat_gc[i] = g
            pick_u = rng_u.choice(np.flatnonzero(usable), size=len(chosen), replace=False)
            _, _, gu = profile_and_gc(gene_regions, region_tf, gc, pick_u)
            unstrat_gc[i] = gu
        # P5g_a and the matching-bites check, both on mean promoter GC.
        extreme_gc = int(np.sum(np.abs(unstrat_gc - np.median(unstrat_gc))
                                >= abs(gc_observed - np.median(unstrat_gc))))
        matching.append({
            "program_uid": uid, "status": "computed", "n_genes": int(len(chosen)),
            "n_profile_regions": int(n_observed),
            "program_mean_gc": gc_observed,
            "unstratified_null_mean_gc_median": float(np.median(unstrat_gc)),
            "unstratified_null_gc_lo": float(np.quantile(unstrat_gc, 0.025)),
            "unstratified_null_gc_hi": float(np.quantile(unstrat_gc, 0.975)),
            "stratified_null_mean_gc_median": float(np.median(strat_gc)),
            "gc_gap_unstratified": abs(gc_observed - float(np.median(unstrat_gc))),
            "gc_gap_stratified": abs(gc_observed - float(np.median(strat_gc))),
            "p_nominal_gc_vs_unstratified": float((extreme_gc + 1) / (draws_n + 1)),
        })
        for track, j in track_index.items():
            column = strat[:, j]
            centre = float(np.median(column))
            extreme = int(np.sum(np.abs(column - centre) >= abs(observed[j] - centre)))
            rows.append({"program_uid": uid, "group": track.split("|", 1)[0],
                         "track": track.split("|", 1)[1],
                         "transcription_factor": factor.get(track, ""),
                         "n_genes": int(len(chosen)), "n_profile_regions": int(n_observed),
                         "share_program": float(observed[j]),
                         "share_gc_matched_null_median": centre,
                         "null_lo": float(np.quantile(column, 0.025)),
                         "null_hi": float(np.quantile(column, 0.975)),
                         "share_minus_null": float(observed[j] - centre),
                         "n_draws": draws_n,
                         "p_nominal": float((extreme + 1) / (draws_n + 1))})

    frame = pd.DataFrame(rows)
    if len(frame):
        frame["BH_q_within_program"] = np.concatenate([
            false_discovery_control(part.p_nominal.to_numpy(), method="bh")
            for _, part in frame.groupby("program_uid", sort=False)])
    frame.to_csv(args.out / "F5g_program_tf_profiles_gc_matched.tsv.gz", sep="\t", index=False)
    match = pd.DataFrame(matching)
    match.to_csv(args.out / "F5g_gc_matching.tsv", sep="\t", index=False)

    computed = match.loc[match.status.eq("computed")] if len(match) else match
    bites = bool(len(computed)) and bool(
        (computed.gc_gap_stratified <= computed.gc_gap_unstratified).all())
    summary = {
        "status": "f5g_gc_matched_sensitivity_complete",
        "prespecification_sha256": hashlib.sha256(PRESPEC.read_bytes()).hexdigest(),
        "seed": seed, "draws": draws_n,
        "programs_computed": int(len(computed)),
        "contrasts": int(len(frame)),
        "contrasts_bh_q_under_0p05": int((frame.BH_q_within_program < 0.05).sum()) if len(frame) else 0,
        "p_floor": 1.0 / (draws_n + 1),
        "matching_bites": bites,
        "matching_bites_meaning": ("every computable program's mean promoter GC is at least as close to the "
                                   "stratified null as to the unstratified one; False voids the sensitivity"),
        "replaces_F5": False,
        "nominates_nothing_for_followup": True,
        "adopted": False,
    }
    (args.out / "f5g_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()

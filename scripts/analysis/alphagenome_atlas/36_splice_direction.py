#!/usr/bin/env python3
"""Step 36 (P3a, splice arm third form): a SIGNED junction effect, through the model API.

Prespecification 36_splice_direction_prespec.json. The Atlas serves no signed splice channel
(scorer_metadata is_signed=False for all three splice scorers, and every deposited junction quantile is
positive), so the deposited sQTL direction numbers could not have worked. The model API returns predicted
junction USAGE for reference and alternate sequence separately, which gives a direction.

The primary statistic is LeafCutter's own estimand, because that is what a GTEx sQTL slope is defined on:
the target junction's share of its cluster's usage, in the alternate arm against the reference arm. The
cluster is the junctions sharing the target's donor or acceptor. Arms are joined on junction coordinates,
never by row index: the model proposes a different junction set for each sequence.

usage: 36_splice_direction.py measure <p3a_run> [n_variants]   # coordinate convention only, no statistics
       36_splice_direction.py run <p3a_run> [n_variants]
Outputs (tables/): splice_direction_pairs.tsv, splice_direction.json, splice_direction_offsets.json
"""

from __future__ import annotations

import collections
import importlib.util
import json
import math
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import lib_atlas as la

HERE = pathlib.Path(__file__).resolve().parent
PRESPEC = HERE / "36_splice_direction_prespec.json"
FASTA = ("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
         "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")
LIVER = ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]
# GTEx closed intron (s, e) <-> model junction (s + START_OFFSET, e + END_OFFSET). The expected value comes
# from lib_atlas.leafcutter_intron_key, measured on two LeafCutter substrates; `measure` re-measures it here
# on the model API's own junction coordinates and the run refuses to proceed if it disagrees.
EXPECTED_OFFSETS = (0, -1)
NBOOT = 2000
SEED = 123
MIN_BLOCKS = 20
EPS = 1e-9


def load_step35():
    """Reuse step 35's archive loader and rank machinery; the coordinate conversion must not be duplicated."""
    spec = importlib.util.spec_from_file_location("step35", HERE / "35_splice_detection.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def junction_coords(j) -> tuple:
    """(start, end) of a model-API junction, window-relative. Objects first, string form as the fallback."""
    s, e = getattr(j, "start", None), getattr(j, "end", None)
    if s is not None and e is not None:
        return int(s), int(e)
    body = str(j).split(":")
    if len(body) < 2 or "-" not in body[1]:
        raise la.ContractError(f"cannot parse junction: {j}")
    a, b = body[1].split("-")
    return int(a), int(b)


def junction_usage(output, tag: str) -> tuple[dict, int]:
    """(start, end) -> per-liver-track usage vector, for one arm. Window-relative coordinates."""
    d = getattr(output, "splice_junctions", None)
    if d is None or d.values is None or not np.size(d.values):
        raise la.ContractError(f"{tag}: no splice_junctions in the response")
    values = np.asarray(d.values, dtype=float)
    out = {}
    for i, j in enumerate(d.junctions):
        out[junction_coords(j)] = values[i]
    return out, values.shape[1]


def site_usage_at(output, sites: list[int], pad: int = 1) -> float:
    """Mean over liver tracks of the per-base splice-site usage at each site and its immediate neighbours.

    The +/- 1 pad is not smoothing for its own sake: the junction coordinate and the per-base index agree to
    within the donor/acceptor convention, and a one-base frame difference would otherwise read the wrong
    base and turn a real effect into noise. Both ends of the intron are weighted equally.
    """
    d = getattr(output, "splice_site_usage", None)
    if d is None or d.values is None or not np.size(d.values):
        return math.nan
    v = np.asarray(d.values, dtype=float)
    keep = [o for s in sites for o in range(s - pad, s + pad + 1) if 0 <= o < v.shape[0]]
    if not keep:
        return math.nan
    return float(np.nanmean(v[keep, :]))


def cluster_index(keys) -> tuple[dict, dict]:
    """donor -> junctions and acceptor -> junctions, built once per variant.

    A 1-Mb window returns 13,000 to 265,000 junctions, so scanning the whole set for every junction whose
    cluster is needed is the difference between seconds and hours over the substrate.
    """
    by_start, by_end = collections.defaultdict(list), collections.defaultdict(list)
    for k in keys:
        by_start[k[0]].append(k)
        by_end[k[1]].append(k)
    return by_start, by_end


def cluster_of(target: tuple, index: tuple[dict, dict]) -> list:
    """Junctions sharing the target's donor (start) or acceptor (end). LeafCutter's definition."""
    by_start, by_end = index
    # a junction cannot appear in both lists unless it shares the target's start AND end, which makes it the
    # target itself, so `k != target` is the only exclusion needed and no de-duplication is possible here
    out = [k for k in by_start.get(target[0], ()) if k != target]
    out += [k for k in by_end.get(target[1], ()) if k != target]
    return out


def excision_ratio_log2(target: tuple, cluster: list, ref: dict, alt: dict) -> float:
    """Per track, log2 of the target's share of cluster usage in alt over the same share in ref.

    The share, not the level: an allele that lifts every junction in a cluster equally changes no excision
    ratio, which is exactly what a LeafCutter phenotype measures. Averaged over liver tracks afterwards so
    one track with a near-zero denominator cannot dominate.
    """
    members = [target] + list(cluster)
    if len(members) < 2:
        return math.nan
    rs = np.sum([ref[k] for k in members], axis=0)
    as_ = np.sum([alt[k] for k in members], axis=0)
    r = (ref[target] + EPS) / (rs + EPS)
    a = (alt[target] + EPS) / (as_ + EPS)
    with np.errstate(divide="ignore", invalid="ignore"):
        lr = np.log2((a + EPS) / (r + EPS))
    return float(np.nanmean(lr))


def usage_log2(target: tuple, ref: dict, alt: dict) -> float:
    """Mean over liver tracks of log2(alt usage / ref usage) for the target junction alone."""
    with np.errstate(divide="ignore", invalid="ignore"):
        lr = np.log2((alt[target] + EPS) / (ref[target] + EPS))
    return float(np.nanmean(lr))


def measure_offsets(pairs: pd.DataFrame, ref_keys_by_variant: dict, start0_by_variant: dict) -> dict:
    """Match counts for every (start, end) offset in -2..2, over the GTEx targets of every variant."""
    counts = {}
    for a in range(-2, 3):
        for b in range(-2, 3):
            hit = 0
            for _, r in pairs.iterrows():
                keys = ref_keys_by_variant.get(r["variant_uid"])
                if keys is None:
                    continue
                s0 = start0_by_variant[r["variant_uid"]]
                if (int(r["intron_start"]) + a - s0, int(r["intron_end"]) + b - s0) in keys:
                    hit += 1
            counts[f"{a},{b}"] = hit
    best = max(counts, key=lambda k: counts[k])
    ties = [k for k, v in counts.items() if v == counts[best]]
    return {"counts": counts, "best": best, "n_tied_at_best": len(ties),
            "expected": f"{EXPECTED_OFFSETS[0]},{EXPECTED_OFFSETS[1]}",
            "agrees_with_expected": best == f"{EXPECTED_OFFSETS[0]},{EXPECTED_OFFSETS[1]}" and len(ties) == 1}


def block_concordance(pred: np.ndarray, meas: np.ndarray, blocks: np.ndarray,
                      nboot: int = NBOOT, seed: int = SEED) -> dict:
    """Sign concordance with its marginal expectation, a per-block binomial, and a block bootstrap."""
    from scipy import stats
    ok = np.isfinite(pred) & np.isfinite(meas) & (pred != 0)
    p, m, b = pred[ok], meas[ok], blocks[ok]
    if not len(p):
        return {"n": 0, "concordance": None}
    agree = ((p > 0) == (m > 0)).astype(float)
    labels = pd.unique(b)
    per_block = pd.Series(agree).groupby(pd.Series(b)).mean()
    above = int((per_block > 0.5).sum())
    decided = int((per_block != 0.5).sum())
    rng = np.random.default_rng(seed)
    index = {lab: np.where(b == lab)[0] for lab in labels}
    draws = np.empty(nboot)
    for d in range(nboot):
        idx = np.concatenate([index[lab] for lab in rng.choice(labels, size=len(labels), replace=True)])
        draws[d] = float(np.mean(agree[idx]))
    return {"n": int(len(p)), "n_blocks": int(len(labels)), "concordance": float(np.mean(agree)),
            "marginal_expectation": float(la.marginal_expected_concordance(p, m)),
            "above_marginal": float(np.mean(agree) - la.marginal_expected_concordance(p, m)),
            "share_pred_positive": float(np.mean(p > 0)), "share_meas_positive": float(np.mean(m > 0)),
            "blocks_above_half": above, "blocks_decided": decided,
            "block_binomial_p": float(stats.binomtest(above, decided, 0.5).pvalue) if decided else None,
            "block_lo": float(np.percentile(draws, 2.5)), "block_hi": float(np.percentile(draws, 97.5))}


def column(df: pd.DataFrame, name: str) -> np.ndarray:
    """The column as float, or all-missing of the right length when no row produced it."""
    if name in df:
        return df[name].to_numpy(float)
    return np.full(len(df), np.nan)


def verdicts(out: dict) -> dict:
    """The written predictions, read off the numbers just computed."""
    ratio = out["direction"]["excision_ratio_log2"]
    raw = out["direction"]["target_usage_log2"]
    site = out["direction"]["site_usage_delta"]
    det = out["detection"]["rank_abs_excision_ratio"]
    atlas = out["detection"]["rank_abs_atlas_quantile"]
    c = ratio.get("concordance")
    p1 = bool(c is not None and c > 0.65 and ratio["above_marginal"] > 0.05
              and ratio["block_binomial_p"] is not None and ratio["block_binomial_p"] < 0.01)
    p2 = bool(c is not None and raw.get("concordance") is not None and raw["concordance"] <= c - 0.03)
    p3 = bool(det.get("mean") is not None and det["mean"] > 0.70 and det["mean"] > 0.6053)
    p4 = bool(c is not None and site.get("concordance") is not None and abs(site["concordance"] - c) < 0.05)
    p5 = bool(det.get("mean") is not None and atlas.get("mean") is not None
              and det["mean"] - atlas["mean"] >= 0.15)
    if ratio.get("n_blocks", 0) < MIN_BLOCKS:
        rule = "D3: fewer than 20 blocks; indeterminate, not negative"
    elif not p1:
        rule = ("D1: the splice channel shows no direction on measured liver splicing at the training-data "
                "ceiling. P6c is not built as a prioritisation, nothing is ranked by a predicted splice "
                "effect, and the HSD17B13 splice magnitude must be reported as an unvalidated channel.")
    else:
        rule = "D2: direction holds; P6c may be built on this model-API statistic with the proximity baseline printed"
    return {"P1_direction_recovered": p1, "P2_cluster_normalisation_helps": p2,
            "P3_detection_above_0.70_and_beats_proximity": p3,
            "P4_site_usage_agrees_with_ratio": p4, "P5_model_api_beats_atlas_detection": p5,
            "decision": rule}


def main() -> None:
    import hashlib

    import pysam
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    import atlas_query as aq

    if len(sys.argv) < 3 or sys.argv[1] not in ("measure", "run"):
        raise la.ContractError("usage: 36_splice_direction.py measure|run <p3a_run> [n_variants]")
    mode, p3a = sys.argv[1], pathlib.Path(sys.argv[2]).resolve()
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    tables = la.out_root() / "tables"
    s35 = load_step35()

    atlas_junc = s35.load_junctions(p3a / "raw" / "atlas_sqtl")
    pairs = s35.sqtl_pairs()
    rows35, counts35 = s35.build_rows(pairs, atlas_junc)
    # keyed on the phenotype, not on the intron: two GTEx phenotypes can name the same intron under
    # different gene ids, and step 35 scored one and dropped the other. Keying on the intron alone
    # re-admitted 15 pairs step 35 never scored, which would unpair the comparison P5 asks for.
    atlas_rank = {(r["variant_uid"], r["phenotype_id"]): r for r in rows35}
    pairs = pairs[[(r.variant_uid, r.phenotype_id) in atlas_rank
                   for r in pairs.itertuples()]].reset_index(drop=True)
    la.log(f"step 36: {len(pairs)} pairs / {pairs['variant_uid'].nunique()} variants carried over from step 35")

    variants = list(dict.fromkeys(pairs["variant_uid"]))
    if limit:
        variants = variants[:limit]
        pairs = pairs[pairs["variant_uid"].isin(set(variants))].reset_index(drop=True)

    key = la.load_api_key()
    la.log(f"step 36: api key length {len(key[0])} from {key[1]}")
    model = dna_client.create(key[0])
    window = dna_client.SEQUENCE_LENGTH_1MB
    outs = [dna_client.OutputType.SPLICE_JUNCTIONS, dna_client.OutputType.SPLICE_SITE_USAGE]
    fasta = pysam.FastaFile(FASTA)

    by_variant = collections.defaultdict(list)
    for _, r in pairs.iterrows():
        by_variant[r["variant_uid"]].append(r)
    # the candidate pool step 35 ranked against, indexed once: scanning 122,687 archive keys per row is
    # 80 million comparisons for no reason
    atlas_pool = collections.defaultdict(list)
    for (u, s, e), genes in atlas_junc.items():
        for g in genes:
            atlas_pool[(u, g)].append((s, e))

    ref_keys, start0s, rows, n_tracks_seen = {}, {}, [], set()
    for n, uid in enumerate(variants, 1):
        chrom, pos, ref, alt = uid.split(":")
        pos1 = int(pos)
        start0 = max(0, pos1 - 1 - window // 2)
        start0s[uid] = start0
        seq = fasta.fetch(chrom, start0, start0 + window).upper()
        if len(seq) != window:
            rows.append({"variant_uid": uid, "state": "window_off_chromosome"})
            continue
        at = pos1 - 1 - start0
        if seq[at] != ref:
            rows.append({"variant_uid": uid, "state": f"reference_mismatch: fasta {seq[at]} vs gtex {ref}"})
            continue
        alt_seq = seq[:at] + alt + seq[at + 1:]
        interval = genome.Interval(chromosome=chrom, start=start0, end=start0 + window)
        got = {}
        arms = (("ref", seq),) if mode == "measure" else (("ref", seq), ("alt", alt_seq))
        for arm, s in arms:
            o = aq.call_with_quota_retry(
                lambda q=s: model.predict_sequence(sequence=q, requested_outputs=outs,
                                                   ontology_terms=LIVER, interval=interval),
                label=f"{uid}:{arm}")
            got[arm] = o
        ju_ref, nt = junction_usage(got["ref"], f"{uid}:ref")
        ref_keys[uid] = set(ju_ref)
        if mode == "measure":
            n_tracks_seen.add((nt, nt))
            la.log(f"step 36 measure: {n}/{len(variants)} {uid} {len(ju_ref)} junctions")
            continue
        ju_alt, nt2 = junction_usage(got["alt"], f"{uid}:alt")
        n_tracks_seen.add((nt, nt2))
        shared = set(ju_ref) & set(ju_alt)
        index = cluster_index(shared)
        for r in by_variant[uid]:
            gs, ge = int(r["intron_start"]), int(r["intron_end"])
            target = (gs + EXPECTED_OFFSETS[0] - start0, ge + EXPECTED_OFFSETS[1] - start0)
            row = {"variant_uid": uid, "phenotype_id": r["phenotype_id"], "gene": r["gene"],
                   "intron_start": gs, "intron_end": ge, "chrom": r["chrom"], "position": pos1,
                   "slope": float(r["slope"]), "pval_nominal": float(r["pval_nominal"]),
                   "is_sgene_lead": bool(r["is_sgene_lead"]),
                   "n_junctions_ref": len(ju_ref), "n_junctions_alt": len(ju_alt), "n_shared": len(shared)}
            if target not in shared:
                row["state"] = "target_junction_not_predicted_in_both_arms"
                rows.append(row)
                continue
            cluster = cluster_of(target, index)
            row["n_cluster"] = len(cluster)
            row["excision_ratio_log2"] = excision_ratio_log2(target, cluster, ju_ref, ju_alt)
            row["target_usage_log2"] = usage_log2(target, ju_ref, ju_alt)
            row["target_usage_ref"] = float(np.nanmean(ju_ref[target]))
            # donor and acceptor bases of the target intron, window-relative
            sites = [target[0], target[1]]
            row["site_usage_delta"] = site_usage_at(got["alt"], sites) - site_usage_at(got["ref"], sites)
            # detection on the SAME pool step 35 used: the Atlas-served junctions of this variant and gene
            a = atlas_rank[(uid, r["phenotype_id"])]
            pool = [(s + EXPECTED_OFFSETS[0] - start0, e + EXPECTED_OFFSETS[1] - start0)
                    for (s, e) in atlas_pool.get((uid, r["gene"]), []) if (s, e) != (gs, ge)]
            pool = [k for k in pool if k in shared]
            row["n_pool_model"] = len(pool)
            row["rank_abs_atlas_quantile"] = a["rank_abs_quantile"]
            row["rank_proximity"] = a["rank_proximity"]
            if len(pool) >= 2:
                stat = [abs(excision_ratio_log2(k, cluster_of(k, index), ju_ref, ju_alt)) for k in pool]
                row["rank_abs_excision_ratio"] = s35.percentile_rank(abs(row["excision_ratio_log2"]), stat)
                row["rank_abs_usage_log2"] = s35.percentile_rank(
                    abs(row["target_usage_log2"]), [abs(usage_log2(k, ju_ref, ju_alt)) for k in pool])
            row["state"] = "scored"
            rows.append(row)
        if n % 25 == 0 or n == len(variants):
            la.log(f"step 36: {n}/{len(variants)} variants, {sum(1 for x in rows if x.get('state') == 'scored')} scored")

    offsets = measure_offsets(pairs, ref_keys, start0s)
    json.dump({"offsets": offsets, "n_variants_probed": len(ref_keys), "n_tracks": sorted(n_tracks_seen)},
              (tables / "splice_direction_offsets.json").open("w"), indent=1, default=float)
    la.log(f"step 36: offsets best {offsets['best']} (expected {offsets['expected']}), "
           f"agrees {offsets['agrees_with_expected']}")
    if mode == "measure":
        return
    if not offsets["agrees_with_expected"]:
        raise la.ContractError(f"junction coordinate convention disagrees with the measured one: {offsets}")

    cols = sorted({k for r in rows for k in r})
    la.write_tsv_once(tables / "splice_direction_pairs.tsv", rows, cols)
    df = pd.DataFrame([r for r in rows if r.get("state") == "scored"])
    if df.empty:
        raise la.ContractError("no pair scored")
    blocks = la.coarse_blocks([(f"{i}", r.chrom, int(r.position)) for i, r in enumerate(df.itertuples())])
    df["block"] = [blocks[f"{i}"] for i in range(len(df))]
    b = df["block"].to_numpy()
    meas = df["slope"].to_numpy(float)
    out = {"prespec_sha256": hashlib.sha256(PRESPEC.read_bytes()).hexdigest(),
           "p3a_run": str(p3a), "step35_counts": counts35, "n_pairs_carried": int(len(pairs)),
           "n_variants_predicted": len(variants), "n_scored": int(len(df)),
           "states": pd.Series([r.get("state", "") for r in rows]).value_counts().to_dict(),
           "median_cluster_size": float(df["n_cluster"].median()),
           "offsets": offsets,
           "direction": {c: block_concordance(column(df, c), meas, b)
                         for c in ("excision_ratio_log2", "target_usage_log2", "site_usage_delta")},
           "detection": {c: s35.block_mean(column(df, c), b)
                         for c in ("rank_abs_excision_ratio", "rank_abs_usage_log2",
                                   "rank_abs_atlas_quantile", "rank_proximity")}}
    out["detection_model_minus_atlas"] = s35.block_paired_difference(
        column(df, "rank_abs_excision_ratio"),
        column(df, "rank_abs_atlas_quantile"), b)
    out["detection_model_minus_proximity"] = s35.block_paired_difference(
        column(df, "rank_abs_excision_ratio"),
        column(df, "rank_proximity"), b)
    leads = df[df["is_sgene_lead"]]
    out["sgene_leads_only"] = {"direction": block_concordance(leads["excision_ratio_log2"].to_numpy(float),
                                                             leads["slope"].to_numpy(float),
                                                             leads["block"].to_numpy())}
    out["magnitude"] = {"spearman_abs_slope_vs_abs_ratio": float(
        pd.Series(np.abs(meas)).corr(pd.Series(np.abs(df["excision_ratio_log2"].to_numpy(float))),
                                     method="spearman"))}
    out["exposure"] = ("documented_overlap: GTEx v8 is an AlphaGenome training source; a ceiling, not "
                       "evidence of generalisation")
    out["route"] = ("model API predict_sequence, predicted track LEVELS; never pooled with an Atlas "
                    "calibrated quantile")
    out["verdicts"] = verdicts(out)
    json.dump(out, (tables / "splice_direction.json").open("w"), indent=1, default=float)
    la.log(f"step 36: {len(df)} pairs, concordance "
           f"{out['direction']['excision_ratio_log2']['concordance']:.4f} "
           f"(marginal {out['direction']['excision_ratio_log2']['marginal_expectation']:.4f}); "
           f"{out['verdicts']['decision'][:40]}")


if __name__ == "__main__":
    main()

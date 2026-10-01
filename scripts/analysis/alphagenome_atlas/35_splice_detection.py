#!/usr/bin/env python3
"""Step 35 (P3a, splice arm second form): does the splice channel find the RIGHT junction?

Prespecification 35_splice_detection_prespec.json, written before any rank was computed. The deposited
sQTL arm tested direction and failed (sign concordance 0.538 on the matched within-cluster contrast). P6c
needs only junction LOCATION, so location is measured here before any P6c code exists.

For each GTEx liver sQTL pair whose target intron the Atlas scored, the target is ranked against the other
introns the SAME variant was scored for in the SAME gene. Holding the variant fixed removes the variant's
own magnitude, the sequence context and the model from the comparison; what is left is whether the channel
points at the intron that actually moved. Three ranking statistics are compared on the same pairs: the
absolute junction quantile, the estimand-matched within-cluster contrast, and proximity to the variant.

usage: 35_splice_detection.py <p3a_run_root>
Outputs (tables/): splice_detection_pairs.tsv, splice_detection.json
"""

from __future__ import annotations

import collections
import json
import math
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import atlas_archive as aa
import lib_atlas as la

PRESPEC = pathlib.Path(__file__).resolve().parent / "35_splice_detection_prespec.json"
SQTL_DIR = la.PROJECT / "data/external/gtex_v8_liver_sqtl/GTEx_Analysis_v8_sQTL"
LIVER = ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]
NBOOT = 2000
SEED = 123
MIN_POOL = 2
MIN_BLOCKS = 20


def liver_track_mask(var: pd.DataFrame) -> np.ndarray:
    """Adult and embryonic liver-lineage tracks, the same filter the deposited sQTL arm used."""
    curie = var["ontology_curie"].astype(str) if "ontology_curie" in var else pd.Series([""] * len(var))
    return curie.isin(LIVER).to_numpy()


def gtex_intron_key(junction_start: int, junction_end: int) -> tuple:
    """Atlas (half-open) junction -> GTEx LeafCutter (closed) intron, measured offsets, no tolerance."""
    return la.leafcutter_intron_key(junction_start, junction_end)


def load_junctions(raw_dir: pathlib.Path) -> dict:
    """(variant_uid, intron_start, intron_end) -> {gene, quantile} over the liver tracks.

    The archive holds one row per variant x junction; a junction can be served under more than one gene id
    (overlapping annotations), so the gene is kept per row and the pool is built per gene, not per variant.
    """
    out: dict = {}
    chunks = sorted(raw_dir.glob("chunk_*"))
    if not chunks:
        raise la.ContractError(f"no archive chunks under {raw_dir}")
    for chunk in chunks:
        path = chunk / "SPLICE_JUNCTIONS.h5ad"
        if not path.exists():
            continue
        import anndata
        a = anndata.read_h5ad(path)
        if not a.shape[0] or "variant" not in a.obs:
            continue
        obs = a.obs.reset_index(drop=True)
        mask = liver_track_mask(a.var.reset_index(drop=True))
        if not mask.any():
            raise la.ContractError("no liver track survived the mask")
        q = np.asarray(a.layers["quantiles"], np.float32)[:, mask]
        for i in range(len(obs)):
            uid = aa.variant_uid_from_str(str(obs.at[i, "variant"]))
            key = (uid, *gtex_intron_key(int(float(obs.at[i, "junction_Start"])),
                                         int(float(obs.at[i, "junction_End"]))))
            gene = str(obs.at[i, "gene_id"]).split(".")[0]
            value = float(np.nanmedian(q[i]))
            out.setdefault(key, {})[gene] = value
    return out


def sqtl_pairs() -> pd.DataFrame:
    """Per LeafCutter cluster the 10 most significant liver sQTL pairs, SNVs only (step 33's selection)."""
    sg = pd.read_csv(SQTL_DIR / "Liver.v8.sgenes.txt.gz", sep="\t", dtype=str)
    pairs = pd.read_csv(SQTL_DIR / "Liver.v8.sqtl_signifpairs.txt.gz", sep="\t", dtype=str)
    pairs["cluster"] = pairs["phenotype_id"].str.split(":").str[3]
    pairs["pval_nominal"] = pd.to_numeric(pairs["pval_nominal"], errors="coerce")
    pairs = pairs.sort_values("pval_nominal").groupby("cluster").head(10)

    def uid(value: str):
        chrom, pos, ref, alt, _ = value.split("_")
        return f"{chrom}:{pos}:{ref}:{alt}" if len(ref) == 1 and len(alt) == 1 else None

    pairs["variant_uid"] = pairs["variant_id"].map(uid)
    pairs = pairs.dropna(subset=["variant_uid"]).copy()
    fields = pairs["phenotype_id"].str.split(":", expand=True)
    pairs["intron_start"] = fields[1].astype(int)
    pairs["intron_end"] = fields[2].astype(int)
    pairs["gene"] = fields[4].str.split(".").str[0]
    pairs["chrom"] = pairs["variant_uid"].str.split(":").str[0]
    pairs["position"] = pairs["variant_uid"].str.split(":").str[1].astype(int)
    pairs["slope"] = pd.to_numeric(pairs["slope"], errors="coerce")
    pairs["is_sgene_lead"] = pairs["variant_id"].isin(set(sg["variant_id"]))
    return pairs


def percentile_rank(target: float, pool: list) -> float:
    """(#pool strictly below + 0.5 x #ties) / #pool. Chance is 0.5 whatever the pool size."""
    values = [v for v in pool if v == v]
    if not values or target != target:
        return math.nan
    below = sum(1 for v in values if v < target)
    ties = sum(1 for v in values if v == target)
    return (below + 0.5 * ties) / len(values)


def intron_distance(position: int, start: int, end: int) -> int:
    """bp from the variant to the nearer end of the intron; 0 when the variant is inside it."""
    if start <= position <= end:
        return 0
    return int(min(abs(position - start), abs(position - end)))


def cluster_contrast(target_key: tuple, gene: str, by_uid: dict, junc: dict) -> float:
    """|target - mean(siblings)|; siblings share the target's donor or acceptor within the same gene."""
    uid, ts, te = target_key
    t = junc.get(target_key, {}).get(gene)
    if t is None or t != t:
        return math.nan
    sibs = [junc[(uid, s, e)][gene] for (s, e) in by_uid.get((uid, gene), [])
            if (s, e) != (ts, te) and (s == ts or e == te)
            and junc[(uid, s, e)].get(gene) == junc[(uid, s, e)].get(gene)]
    if not sibs:
        return math.nan
    return float(t - np.mean(sibs))


def build_rows(pairs: pd.DataFrame, junc: dict) -> tuple[list, dict]:
    """One row per (variant, target intron) pair that the Atlas scored, with all three ranks."""
    by_uid: dict = collections.defaultdict(list)
    for (uid, s, e), genes in junc.items():
        for gene in genes:
            by_uid[(uid, gene)].append((s, e))
    counts = {"pairs_considered": int(len(pairs)), "target_not_scored": 0, "pool_too_small": 0}
    rows = []
    for _, r in pairs.iterrows():
        uid, gene = r["variant_uid"], r["gene"]
        target_key = (uid, int(r["intron_start"]), int(r["intron_end"]))
        target_q = junc.get(target_key, {}).get(gene)
        if target_q is None or target_q != target_q:
            counts["target_not_scored"] += 1
            continue
        pool_keys = [(s, e) for (s, e) in by_uid.get((uid, gene), []) if (s, e) != target_key[1:]]
        if len(pool_keys) < MIN_POOL:
            counts["pool_too_small"] += 1
            continue
        pool_q = [junc[(uid, s, e)][gene] for (s, e) in pool_keys]
        target_c = cluster_contrast(target_key, gene, by_uid, junc)
        pool_c = [cluster_contrast((uid, s, e), gene, by_uid, junc) for (s, e) in pool_keys]
        pos = int(r["position"])
        target_d = intron_distance(pos, int(r["intron_start"]), int(r["intron_end"]))
        pool_d = [intron_distance(pos, s, e) for (s, e) in pool_keys]
        rows.append({
            "variant_uid": uid, "phenotype_id": r["phenotype_id"], "gene": gene,
            "intron_start": int(r["intron_start"]), "intron_end": int(r["intron_end"]),
            "chrom": r["chrom"], "position": pos, "slope": float(r["slope"]),
            "pval_nominal": float(r["pval_nominal"]), "is_sgene_lead": bool(r["is_sgene_lead"]),
            "n_pool": len(pool_keys), "target_junction_quantile": float(target_q),
            "target_leafcutter_contrast": target_c, "target_distance_bp": target_d,
            "rank_abs_quantile": percentile_rank(abs(target_q), [abs(v) for v in pool_q]),
            "rank_abs_contrast": percentile_rank(abs(target_c), [abs(v) for v in pool_c]) if target_c == target_c else math.nan,
            "rank_proximity": percentile_rank(-target_d, [-d for d in pool_d]),
            "n_pool_exact_zero_quantile": int(sum(1 for v in pool_q if v == 0.0)),
        })
    return rows, counts


def block_mean(values: np.ndarray, blocks: np.ndarray, nboot: int = NBOOT, seed: int = SEED) -> dict:
    """Mean and a block-bootstrap interval, resampling 1-Mb blocks (the caQTL arm's estimator)."""
    ok = np.isfinite(values)
    v, b = values[ok], blocks[ok]
    if not len(v):
        return {"n": 0, "n_blocks": 0, "mean": None}
    labels = pd.unique(b)
    rng = np.random.default_rng(seed)
    index = {lab: np.where(b == lab)[0] for lab in labels}
    draws = np.empty(nboot)
    for d in range(nboot):
        pick = rng.choice(labels, size=len(labels), replace=True)
        draws[d] = float(np.mean(v[np.concatenate([index[lab] for lab in pick])]))
    return {"n": int(len(v)), "n_blocks": int(len(labels)), "mean": float(np.mean(v)),
            "block_lo": float(np.percentile(draws, 2.5)), "block_hi": float(np.percentile(draws, 97.5)),
            "block_median": float(np.median(draws)),
            "share_first": float(np.mean(v == 1.0)), "share_above_half": float(np.mean(v > 0.5))}


def block_paired_difference(a: np.ndarray, b: np.ndarray, blocks: np.ndarray,
                            nboot: int = NBOOT, seed: int = SEED) -> dict:
    """Mean(a) - mean(b) on the pairs where BOTH are defined, with a block-bootstrap interval."""
    ok = np.isfinite(a) & np.isfinite(b)
    av, bv, bl = a[ok], b[ok], blocks[ok]
    if not len(av):
        return {"n": 0, "difference": None}
    labels = pd.unique(bl)
    rng = np.random.default_rng(seed)
    index = {lab: np.where(bl == lab)[0] for lab in labels}
    draws = np.empty(nboot)
    for d in range(nboot):
        idx = np.concatenate([index[lab] for lab in rng.choice(labels, size=len(labels), replace=True)])
        draws[d] = float(np.mean(av[idx]) - np.mean(bv[idx]))
    return {"n": int(len(av)), "n_blocks": int(len(labels)),
            "difference": float(np.mean(av) - np.mean(bv)),
            "block_lo": float(np.percentile(draws, 2.5)), "block_hi": float(np.percentile(draws, 97.5))}


def verdicts(out: dict) -> dict:
    """The prespecified predictions, each read off the numbers just computed."""
    primary = out["detection"]["rank_abs_quantile"]
    prox = out["detection"]["rank_proximity"]
    contrast = out["detection"]["rank_abs_contrast"]
    advantage = out["model_minus_proximity"]["rank_abs_quantile"]
    first = out["sign_among_model_first"]
    p1 = bool(primary["mean"] is not None and primary["mean"] > 0.60 and primary["block_lo"] > 0.50)
    p2 = bool(prox["mean"] is not None and prox["mean"] > 0.55)
    p3 = bool(advantage["difference"] is not None and 0 < advantage["difference"] < 0.10)
    # P4 asked whether ranking an intron first also buys its DIRECTION. It cannot be evaluated here: the
    # predictor is SPLICE_JUNCTIONS, which the server declares unsigned and whose every value is positive,
    # so this concordance is identically the share of positive measured slopes. The equality with the
    # marginal expectation is kept as the diagnostic that says so, and the prediction is withdrawn rather
    # than scored, because a boolean here would read as a finding.
    p4 = "withdrawn: the predictor is unsigned, so no direction test is possible (see sign_among_model_first)"
    p5 = bool(contrast["mean"] is not None and primary["mean"] is not None
              and abs(contrast["mean"] - primary["mean"]) < 0.05)
    indeterminate = bool(primary["n_blocks"] < MIN_BLOCKS)
    if indeterminate:
        rule = "D3: fewer than 20 blocks; the arm is indeterminate, not negative"
    elif not p1:
        rule = ("D1: detection FAILED at the training-data ceiling. P6c is not built as a prioritisation; "
                "report coverage only and rank no gene by a predicted splice effect.")
    elif advantage["difference"] is None or advantage["block_lo"] <= 0:
        # D2 turns on whether the advantage is distinguishable from zero, which is a different question
        # from P3's written bound on its size; a +0.02 point estimate whose interval covers 0 is no advantage.
        rule = "D2: detection holds but shows no advantage over proximity; print the proximity baseline beside every P6c number"
    else:
        rule = "D2 not triggered: detection holds and beats proximity; P6c may rank junctions, with the baseline printed"
    return {"P1_detection_above_chance": p1, "P2_proximity_above_chance": p2,
            "P3_model_beats_proximity_by_under_0.10": p3, "P4_direction_when_ranked_first": p4,
            "P5_contrast_ranks_like_absolute": p5, "indeterminate_by_D3": indeterminate,
            "decision": rule}


def main() -> None:
    import hashlib
    if len(sys.argv) < 2:
        raise la.ContractError("usage: 35_splice_detection.py <p3a_run_root>")
    p3a = pathlib.Path(sys.argv[1]).resolve()
    tables = la.out_root() / "tables"
    junc = load_junctions(p3a / "raw" / "atlas_sqtl")
    la.log(f"step 35: {len(junc)} archived variant x junction keys")
    pairs = sqtl_pairs()
    rows, counts = build_rows(pairs, junc)
    if not rows:
        raise la.ContractError("no scored sQTL pair survived the pool rule")
    blocks = la.coarse_blocks([(f"{i}", r["chrom"], r["position"]) for i, r in enumerate(rows)])
    for i, r in enumerate(rows):
        r["block"] = blocks[f"{i}"]
    df = pd.DataFrame(rows)
    b = df["block"].to_numpy()
    out = {"prespec_sha256": hashlib.sha256(PRESPEC.read_bytes()).hexdigest(),
           "p3a_run": str(p3a), "counts": counts, "n_pairs_scored": int(len(df)),
           "n_variants": int(df["variant_uid"].nunique()), "n_genes": int(df["gene"].nunique()),
           "median_pool_size": float(df["n_pool"].median()),
           "detection": {c: block_mean(df[c].to_numpy(float), b)
                         for c in ("rank_abs_quantile", "rank_abs_contrast", "rank_proximity")},
           "model_minus_proximity": {c: block_paired_difference(df[c].to_numpy(float),
                                                                df["rank_proximity"].to_numpy(float), b)
                                     for c in ("rank_abs_quantile", "rank_abs_contrast")}}
    leads = df[df["is_sgene_lead"]]
    out["sgene_leads_only"] = {c: block_mean(leads[c].to_numpy(float), leads["block"].to_numpy())
                               for c in ("rank_abs_quantile", "rank_abs_contrast", "rank_proximity")}
    top = df[df["rank_abs_quantile"] == 1.0]
    out["sign_among_model_first"] = {
        "not_a_direction_test": True,
        "why": ("SPLICE_JUNCTIONS is server-declared unsigned (is_signed=False) and every quantile in the "
                "archive is positive, so this concordance is identically the share of positive measured "
                "slopes. It is reported ONLY because its equality with the marginal expectation is the "
                "evidence for that statement."),
        "n": int(len(top)), "n_blocks": int(top["block"].nunique()),
        "concordance": float(((top["slope"] > 0) == (top["target_junction_quantile"] > 0)).mean()) if len(top) else None,
        "marginal_expectation": float(la.marginal_expected_concordance(
            top["target_junction_quantile"], top["slope"])) if len(top) else None}
    out["distance_diagnostic"] = {
        "median_target_distance_bp": float(df["target_distance_bp"].median()),
        "share_variant_inside_target_intron": float((df["target_distance_bp"] == 0).mean()),
        "share_targets_with_exact_zero_quantile": float((df["target_junction_quantile"] == 0).mean()),
        "share_pool_entries_exact_zero_quantile": float(df["n_pool_exact_zero_quantile"].sum() / df["n_pool"].sum())}
    out["exposure"] = ("documented_overlap: GTEx v8 is an AlphaGenome training source; this is a ceiling on "
                       "the channel, not evidence of generalisation")
    out["verdicts"] = verdicts(out)
    la.write_tsv_once(tables / "splice_detection_pairs.tsv", rows, list(rows[0].keys()))
    json.dump(out, (tables / "splice_detection.json").open("w"), indent=1, default=float)
    la.log(f"step 35: {len(df)} pairs, {out['detection']['rank_abs_quantile']['n_blocks']} blocks, "
           f"rank {out['detection']['rank_abs_quantile']['mean']:.4f} vs proximity "
           f"{out['detection']['rank_proximity']['mean']:.4f}; {out['verdicts']['decision'][:40]}")


if __name__ == "__main__":
    main()

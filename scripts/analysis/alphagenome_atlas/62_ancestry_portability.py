#!/usr/bin/env python3
"""Step 62 (P6b): is the Atlas annotation ancestry-portable, and if not, which step loses the mass?

Every variant receives the same static prediction whatever the ancestry of the study that fine-mapped it,
so nothing here can discover a locus or move a credible set. What can differ by ancestry is how much of a
credible set's posterior mass the annotation ever reaches. That quantity decomposes exactly:

    total posterior mass = queried_mass + excluded_mass_below_floor + unmapped_mass

`unmapped_mass` is hg19->hg38 liftover failure (a property of the source study's coordinates),
`excluded_mass_below_floor` is mass sitting on variants under the query weight floor (a property of
fine-mapping RESOLUTION - a diffuse credible set spreads its mass thin), and `queried_mass` is what the
Atlas was asked about. Within the queried mass, each scorer either served the variant or did not.

Non-EUR credible sets have no matched liver eQTL panel, so the static model is the only functional
prioritisation available for them; whether it covers them equally is the translational-equity question.

Prespecification: 62_ancestry_portability_prespec.json (written before any ancestry contrast).
Outputs (tables/): ancestry_portability.tsv, ancestry_portability_summary.json, ancestry_variant_dossier.tsv.gz
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict

import numpy as np

import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
TRACK0 = la.track0_root() / "tables"
PRESPEC = la.SCRIPT_DIR / "62_ancestry_portability_prespec.json"

MASS_TOLERANCE = 1e-6
MIN_SIGNALS = 10
MIN_BLOCKS = 10
BOOTSTRAP_DRAWS = 2000
SEED = 20260913
CORE_CHANNELS = ["ATAC", "DNASE", "CHIP_HISTONE", "RNA_SEQ", "CHIP_TF", "CAGE", "CONTACT_MAPS", "AVI_SCORE"]
CONDITIONAL_CHANNELS = ["SPLICE_SITE_USAGE", "POLYADENYLATION"]
PRIMARY_STRATUM = ("C_enzyme", "coloc_snp_pp_h4_conditional_on_H4")


def stable_seed(name: str) -> int:
    """Seeds must survive a restart: Python's hash() is salted per process."""
    return int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)


def excluded_mass_by_reason(pairs) -> dict:
    """Posterior mass lost per exclusion reason.

    `signal_query_coverage.unmapped_mass` is total minus MAPPED mass, and `mapped` in step 04 excludes
    indels as well as liftover failures. Over the whole portfolio the split is 207,254 pairs / 49.5 mass
    units `indel_deferred` against 158 pairs / 0.000 units `liftover_failed`: the term is the Atlas point
    query API refusing indels (UNIMPLEMENTED), not a coordinate problem. Calling it liftover would be wrong.
    """
    out = {}
    for r in pairs:
        if r["mapping_status"] == "excluded":
            reason = r["exclusion_reason"] or "unspecified"
            out[reason] = out.get(reason, 0.0) + float(r["weight"])
    return out


def effective_from_sums(sum_w: float, sum_w2: float) -> float:
    if sum_w <= 0 or sum_w2 <= 0:
        raise la.ContractError("effective set size is undefined for zero posterior mass")
    return float(sum_w * sum_w / sum_w2)


def effective_set_size(weights) -> float:
    """Participation ratio 1/sum(p^2): how many variants the posterior is effectively spread over.

    `n_cs_variants` is 0 for every coloc signal (a credible set is a SuSiE object), so it cannot serve as a
    fine-mapping resolution key: matching on it made every EUR signal a comparator for every target and the
    matched contrast reproduced the unmatched one exactly. This is defined for both posterior definitions.
    """
    vals = [float(w) for w in (weights.values() if hasattr(weights, "values") else weights)]
    return effective_from_sums(sum(vals), sum(v * v for v in vals))


def tail_share(values, floor: float) -> float:
    """Share of signals BELOW the coverage the query rule was built to reach.

    The median cannot see this: the query set targets a cumulative mass of 0.999, so an arm with no
    exclusions sits at the ceiling and the informative quantity is how often a signal falls off it.
    """
    v = [float(x) for x in values]
    return float(sum(x < floor for x in v)) / len(v) if v else float("nan")


COVERED_FLOOR = 0.999
UNCOVERED_CEILING = 1e-6


def coverage_state(queried_share: float, indel_share: float) -> str:
    """What a downstream consumer may say about this signal from the Atlas.

    A signal whose posterior sits on indels the point-query API refuses still receives a profile row, and
    that row's top-variant fields are unweighted maxima over the near-zero-weight SNVs that were queried
    instead. The posterior-weighted columns read ~1e-23 and are honest; the top-variant columns are not
    protected by anything. This state is the gate.
    """
    q, i = float(queried_share), float(indel_share)
    if q >= COVERED_FLOOR:
        return "covered"
    if q <= UNCOVERED_CEILING:
        return "uncovered_indel" if i > 0.5 else "uncovered_other"
    return "partial_indel" if i > 0.5 * (1.0 - q) else "partial_other"


def mass_decomposition(row: dict) -> dict:
    """Split a signal's posterior mass into the three places it can go, refusing a row that does not reconstruct."""
    total = float(row["total_mass"])
    queried = float(row["queried_mass"])
    floor = float(row["excluded_mass_below_floor"])
    unmapped = float(row["unmapped_mass"])
    if total <= 0:
        raise la.ContractError(f"signal carries no posterior mass: {row.get('signal_uid', '?')}")
    if abs(queried + floor + unmapped - total) > MASS_TOLERANCE * max(1.0, total):
        raise la.ContractError(
            f"mass identity broken for {row.get('signal_uid', '?')}: "
            f"{queried} + {floor} + {unmapped} != {total}")
    return {"total_mass": total, "queried_share": queried / total,
            "floor_share": floor / total, "excluded_share": unmapped / total}


def channel_available_share(weights: dict, served: dict) -> float:
    """Share of a signal's QUERIED mass that a scorer actually served. Mass-weighted, not variant-counted."""
    tot = float(sum(weights.values()))
    if tot <= 0:
        return float("nan")
    return float(sum(w for v, w in weights.items() if served.get(v))) / tot


def block_draw(rows: list, seed: int) -> list:
    """One row per 1-Mb block: signals at the same locus are not independent units."""
    rng = np.random.default_rng(seed)
    by_block = defaultdict(list)
    for r in rows:
        by_block[r["analysis_block"]].append(r)
    return [bl[int(rng.integers(len(bl)))] for _, bl in sorted(by_block.items())]


def block_bootstrap_median_difference(a_rows: list, b_rows: list, draws: int = BOOTSTRAP_DRAWS,
                                      seed: int = SEED) -> dict:
    """Difference in median value, a minus b, with blocks resampled with replacement in both arms."""
    def collapse(rows):
        by = defaultdict(list)
        for r in rows:
            by[r["analysis_block"]].append(float(r["value"]))
        return [np.asarray(v) for _, v in sorted(by.items())]

    A, B = collapse(a_rows), collapse(b_rows)
    if not A or not B:
        return {"observed": float("nan"), "lo": float("nan"), "hi": float("nan"), "draws": 0}
    rng = np.random.default_rng(seed)
    obs = float(np.median(np.concatenate(A)) - np.median(np.concatenate(B)))
    vals = []
    for _ in range(draws):
        ai = rng.integers(len(A), size=len(A))
        bi = rng.integers(len(B), size=len(B))
        av = [float(A[i][rng.integers(A[i].size)]) for i in ai]
        bv = [float(B[i][rng.integers(B[i].size)]) for i in bi]
        vals.append(float(np.median(av) - np.median(bv)))
    return {"observed": obs, "lo": float(np.quantile(vals, 0.025)), "hi": float(np.quantile(vals, 0.975)),
            "draws": int(draws), "n_blocks_a": len(A), "n_blocks_b": len(B)}


def decile_edges(values) -> np.ndarray:
    v = np.asarray([float(x) for x in values], dtype=float)
    return np.unique(np.quantile(v, np.arange(0.1, 1.0, 0.1)))


def matched_comparators(target: dict, pool: list, edges: np.ndarray, key: str = "effective_set_size") -> list:
    """Pool members in the target's size decile. A target outside the pool's observed range has none."""
    t = float(target[key])
    pv = [float(r[key]) for r in pool]
    if not pv or t < min(pv) or t > max(pv):
        return []
    tb = int(np.searchsorted(edges, t, side="right"))
    return [r for r in pool if int(np.searchsorted(edges, float(r[key]), side="right")) == tb]


def contrast_verdict(n_signals: int, n_blocks: int, min_signals: int = MIN_SIGNALS,
                     min_blocks: int = MIN_BLOCKS) -> dict:
    reasons = []
    if n_blocks < min_blocks:
        reasons.append(f"{n_blocks} independent 1-Mb blocks is below the prespecified minimum of {min_blocks}")
    if n_signals < min_signals:
        reasons.append(f"{n_signals} signals is below the prespecified minimum of {min_signals}")
    return {"conclusive": not reasons, "reason": "; ".join(reasons) or "adequate",
            "n_signals": int(n_signals), "n_blocks": int(n_blocks)}


def paired_matched_difference(targets: list, pool: list, edges: np.ndarray, value_key: str,
                              draws: int = BOOTSTRAP_DRAWS, seed: int = SEED) -> dict:
    """Median of PER-SIGNAL differences (target minus its matched comparators' median), not a difference of medians."""
    diffs, unmatched = [], 0
    for t in targets:
        comp = matched_comparators(t, pool, edges)
        if len(comp) < 3:
            unmatched += 1
            continue
        diffs.append({"analysis_block": t["analysis_block"],
                      "value": float(t[value_key]) - float(np.median([float(c[value_key]) for c in comp])),
                      "n_comparators": len(comp)})
    if not diffs:
        return {"observed": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "n_paired": 0, "n_unmatched": unmatched}
    rng = np.random.default_rng(seed)
    by = defaultdict(list)
    for d in diffs:
        by[d["analysis_block"]].append(d["value"])
    blocks = [np.asarray(v) for _, v in sorted(by.items())]
    vals = []
    for _ in range(draws):
        idx = rng.integers(len(blocks), size=len(blocks))
        vals.append(float(np.median([float(blocks[i][rng.integers(blocks[i].size)]) for i in idx])))
    return {"observed": float(np.median([d["value"] for d in diffs])),
            "lo": float(np.quantile(vals, 0.025)), "hi": float(np.quantile(vals, 0.975)),
            "n_paired": len(diffs), "n_blocks": len(blocks), "n_unmatched": unmatched,
            "median_comparators": float(np.median([d["n_comparators"] for d in diffs]))}


def main() -> None:
    prespec = json.load(PRESPEC.open())
    la.log(f"P6b prespec sha256 {la.sha256_file(PRESPEC)}")

    signals = {s["signal_uid"]: s for s in la.read_tsv(TRACK0 / "eligible_signals.tsv")}
    coverage = {c["signal_uid"]: c for c in la.read_tsv(TRACK0 / "signal_query_coverage.tsv")}

    served = {}
    with la.open_text(TRACK0 / "atlas_availability.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            served[r["variant_uid"]] = r

    qweights = defaultdict(dict)
    sums = defaultdict(lambda: [0.0, 0.0])          # signal -> [sum w, sum w^2] over ALL its variants
    reasons = defaultdict(lambda: defaultdict(float))
    with la.open_text(TRACK0 / "signal_variant_weights.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            w = float(r["weight"])
            acc = sums[r["signal_uid"]]
            acc[0] += w
            acc[1] += w * w
            if r["mapping_status"] == "excluded":
                reasons[r["signal_uid"]][r["exclusion_reason"] or "unspecified"] += w
            if r["in_query_set"] == "True":
                qweights[r["signal_uid"]][r["variant_uid"]] = w

    rows = []
    for uid, s in signals.items():
        cov = coverage.get(uid)
        if cov is None:
            raise la.ContractError(f"signal has no coverage row: {uid}")
        d = mass_decomposition({**cov, "signal_uid": uid})
        w = qweights.get(uid, {})
        row = {"signal_uid": uid, "universe": s["universe"], "ancestry": s["ancestry"],
               "posterior_definition": s["posterior_definition"], "analysis_block": s["analysis_block"],
               "trait": s["trait"], "trait_class": s["trait_class"], "gwas_name": s["gwas_name"],
               "gene": s["gene"], "n_cs_variants": int(s["n_cs_variants"]) if s["n_cs_variants"] else 0,
               "cs_coverage": float(s["cs_coverage"]) if s["cs_coverage"] else float("nan"),
               "n_queried": int(cov["n_queried"]),
               "effective_set_size": effective_from_sums(*sums[uid]), **d}
        row["atlas_coverage_state"] = coverage_state(d["queried_share"],
                                                     reasons.get(uid, {}).get("indel_deferred", 0.0) / d["total_mass"])
        rs = reasons.get(uid, {})
        for reason in ("indel_deferred", "liftover_failed", "multimapped", "unspecified"):
            row[f"excluded_share_{reason}"] = rs.get(reason, 0.0) / d["total_mass"]
        for ch in CORE_CHANNELS + CONDITIONAL_CHANNELS:
            row[f"avail_{ch}"] = channel_available_share(w, {v: served.get(v, {}).get(ch) == "1" for v in w})
        rows.append(row)
    rows.sort(key=lambda r: (r["universe"], r["ancestry"], r["signal_uid"]))
    la.write_tsv_once(TABLES / "ancestry_portability.tsv", rows, list(rows[0].keys()))

    stratum = [r for r in rows if (r["universe"], r["posterior_definition"]) == PRIMARY_STRATUM]
    eur = [r for r in stratum if r["ancestry"] == "EUR"]
    edges = decile_edges([r["effective_set_size"] for r in eur])
    contrasts = {}
    for anc in sorted({r["ancestry"] for r in stratum} - {"EUR"}):
        tgt = [r for r in stratum if r["ancestry"] == anc]
        verdict = contrast_verdict(len(tgt), len({r["analysis_block"] for r in tgt}))
        entry = {"verdict": verdict, "n_signals": len(tgt), "n_blocks": verdict["n_blocks"]}
        if verdict["conclusive"]:
            for key in ("queried_share", "floor_share", "excluded_share", "excluded_share_indel_deferred"):
                entry[f"unmatched_{key}"] = block_bootstrap_median_difference(
                    [{"analysis_block": r["analysis_block"], "value": r[key]} for r in tgt],
                    [{"analysis_block": r["analysis_block"], "value": r[key]} for r in eur],
                    seed=SEED + stable_seed(anc) % 1000)
                entry[f"matched_{key}"] = paired_matched_difference(tgt, eur, edges, key, seed=SEED)
        entry["median_queried_share"] = float(np.median([r["queried_share"] for r in tgt]))
        entry["median_effective_set_size"] = float(np.median([r["effective_set_size"] for r in tgt]))
        entry["POST_HOC_tail_share_below_0.999"] = tail_share([r["queried_share"] for r in tgt], 0.999)
        entry["POST_HOC_median_excluded_share_indel_deferred"] = float(
            np.median([r["excluded_share_indel_deferred"] for r in tgt]))
        entry["POST_HOC_n_signals_losing_over_half_their_mass_to_indels"] = int(
            sum(r["excluded_share_indel_deferred"] > 0.5 for r in tgt))
        for ch in CORE_CHANNELS + CONDITIONAL_CHANNELS:
            entry[f"median_avail_{ch}"] = float(np.nanmedian([r[f"avail_{ch}"] for r in tgt]))
        contrasts[anc] = entry

    eur_entry = {"n_signals": len(eur), "n_blocks": len({r["analysis_block"] for r in eur}),
                 "median_queried_share": float(np.median([r["queried_share"] for r in eur])),
                 "median_effective_set_size": float(np.median([r["effective_set_size"] for r in eur])),
                 "POST_HOC_tail_share_below_0.999": tail_share([r["queried_share"] for r in eur], 0.999),
                 "POST_HOC_median_excluded_share_indel_deferred": float(
                     np.median([r["excluded_share_indel_deferred"] for r in eur])),
                 "POST_HOC_n_signals_losing_over_half_their_mass_to_indels": int(
                     sum(r["excluded_share_indel_deferred"] > 0.5 for r in eur))}
    for ch in CORE_CHANNELS + CONDITIONAL_CHANNELS:
        eur_entry[f"median_avail_{ch}"] = float(np.nanmedian([r[f"avail_{ch}"] for r in eur]))

    listing = defaultdict(lambda: {"n_signals": 0, "blocks": set()})
    for r in rows:
        if r["ancestry"] != "EUR" and (r["universe"], r["posterior_definition"]) != PRIMARY_STRATUM:
            k = f'{r["universe"]}|{r["posterior_definition"]}|{r["ancestry"]}'
            listing[k]["n_signals"] += 1
            listing[k]["blocks"].add(r["analysis_block"])

    states = defaultdict(int)
    for r in rows:
        states[r["atlas_coverage_state"]] += 1
    flagged = [{"signal_uid": r["signal_uid"], "universe": r["universe"], "gwas_name": r["gwas_name"],
                "trait": r["trait"], "gene": r["gene"], "ancestry": r["ancestry"],
                "atlas_coverage_state": r["atlas_coverage_state"], "queried_share": r["queried_share"],
                "excluded_share_indel_deferred": r["excluded_share_indel_deferred"]}
               for r in rows if r["atlas_coverage_state"] != "covered"]
    flagged.sort(key=lambda r: r["queried_share"])
    la.write_tsv_once(TABLES / "atlas_coverage_flags.tsv", flagged, list(flagged[0].keys()))

    dossier = [r for r in rows if r["ancestry"] != "EUR"]
    la.write_tsv_once(TABLES / "ancestry_variant_dossier.tsv", dossier, list(rows[0].keys()))

    summary = {
        "prespec_sha256": la.sha256_file(PRESPEC),
        "predictions_before_results": prespec["predictions_before_results"],
        "primary_stratum": {"universe": PRIMARY_STRATUM[0], "posterior_definition": PRIMARY_STRATUM[1]},
        "EUR_reference": eur_entry,
        "contrasts": contrasts,
        "listing_only_strata": {k: {"n_signals": v["n_signals"], "n_blocks": len(v["blocks"]),
                                    "reason": "below the prespecified block or signal minimum; no contrast computed"}
                                for k, v in sorted(listing.items())},
        "coverage_states": dict(sorted(states.items())),
        "coverage_state_note": ("a signal that is not `covered` still has a profile row; its posterior-weighted "
                                "columns are ~0 and correct, but its TOP-VARIANT columns are unweighted maxima over "
                                "queried variants that carry no posterior mass and must not be read as mechanism"),
        "mass_identity": "total = queried + excluded_below_floor + excluded_unmappable, asserted per signal",
        "naming_correction": ("signal_query_coverage names the third term `unmapped_mass`, which reads as liftover "
                              "failure. It is total minus MAPPED mass, and step 04 marks indels unmapped because the "
                              "Atlas point-query API refuses them (UNIMPLEMENTED). Portfolio-wide the split is "
                              "207,254 pairs / 49.5 mass units indel_deferred against 158 pairs / 0.000 units "
                              "liftover_failed. Reported here as excluded_share with the reason attached."),
        "post_hoc_statistics": ("the prespecified contrast is a difference in MEDIAN coverage. Coverage turned out to "
                                "sit at a construction ceiling (the query rule targets cumulative mass 0.999), so the "
                                "median cannot see the effect and the informative quantity is the lower tail. Tail and "
                                "indel statistics are labelled POST_HOC and are not prespecified tests."),
        "claim_boundary": prespec["claim_boundary"],
    }
    json.dump(summary, (TABLES / "ancestry_portability_summary.json").open("w"), indent=1, default=float)
    la.log(f"P6b: {len(rows)} signals; contrasts {sorted(contrasts)}; dossier {len(dossier)} non-EUR signals")


if __name__ == "__main__":
    main()

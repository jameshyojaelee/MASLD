#!/usr/bin/env python3
"""Step 43 (P4): per-region sequence-sensitivity features from the saturation archive.

For every region of the step-41 universe whose archive is complete, the signed calibrated quantiles of every
possible substitution are reduced to a per-position profile and summarised. Prespecification:
43_saturation_prespec.json (written before any feature was computed).

Reduction, in this order, done twice: on the raw predicted effect divided by the track's own mean level (the
primary scale) and on the calibrated quantile (kept, but it saturates inside measured peaks):
  1. |effect| per (substitution, track);
  2. per (position, track), the MAXIMUM over the up-to-three substitutions at that position, so a profile does
     not depend on which bases a position allows;
  3. per position, the MEAN over the tracks of a group (tracks of one assay in one tissue are replicates);
  4. per region: mean, max and sum over positions, the share of the sum inside the best 50-bp window, and the
     fraction of positions above 0.10 (relative scale) or 0.25 (quantile scale).
Per region the five TF tracks with the largest relative step-2 mean are written, adult liver and HepG2 apart.

Track groups use lib_atlas.track_class, the rule every other step uses. CL:0000182 hepatocyte tracks are
in vitro differentiated in this Atlas build and are never pooled with adult liver.

usage: 43_saturation_features.py <saturation_run_dir> [n_workers] [max_regions]
Outputs (tables/): saturation_region_features.tsv.gz, saturation_region_tf_top5.tsv.gz,
saturation_features_summary.json
"""

from __future__ import annotations

import json
import pathlib
import re
import sys
from collections import Counter
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import lib_atlas as la

SCORERS = ("ATAC", "DNASE", "CHIP_HISTONE", "CHIP_TF", "CAGE", "AVI_SCORE")
VARIANT_RE = re.compile(r"^(chr[^:]+):(\d+):([ACGTN])>([ACGTN])$")
WINDOW_BP = 50
HIGH = {"rel": 0.10, "q": 0.25}
SCALES = ("rel", "q")
TOP_TF = 5
GROUPS = ("liver_atac", "liver_dnase", "liver_h3k27ac", "hepg2_atac", "cage_liver", "avi",
          "chip_tf_liver", "chip_tf_hepg2")
STATS = ("mean", "max", "sum", "share_best50", "frac_above_high")


def region_dir(raw_root: pathlib.Path, universe: str, chrom: str, start0: int, end: int) -> pathlib.Path:
    return raw_root / universe / f"{chrom}_{start0}-{end}"


def group_columns(scorer: str, var: dict) -> dict:
    """Column indices of each track group in one scorer's var table (a dict of equal-length lists)."""
    n = len(next(iter(var.values()))) if var else 0
    blank = [""] * n
    cls = [la.track_class(str(o), str(b), str(t)) for o, b, t in
           zip(var.get("ontology_curie", blank), var.get("biosample_name", blank), var.get("biosample_type", blank))]
    idx = np.arange(n)

    def pick(mask):
        return idx[np.asarray(mask, dtype=bool)] if n else idx

    if scorer == "ATAC":
        return {"liver_atac": pick([c == "primary_liver" for c in cls]),
                "hepg2_atac": pick([c == "HepG2" for c in cls])}
    if scorer == "DNASE":
        return {"liver_dnase": pick([c == "primary_liver" for c in cls])}
    if scorer == "CHIP_HISTONE":
        marks = [str(m).upper() for m in var.get("histone_mark", blank)]
        return {"liver_h3k27ac": pick([c == "primary_liver" and m == "H3K27AC" for c, m in zip(cls, marks)])}
    if scorer == "CAGE":
        return {"cage_liver": pick([c == "primary_liver" for c in cls])}
    if scorer == "AVI_SCORE":
        return {"avi": idx}
    if scorer == "CHIP_TF":
        return {"chip_tf_liver": pick([c == "primary_liver" for c in cls]),
                "chip_tf_hepg2": pick([c == "HepG2" for c in cls])}
    raise la.ContractError(f"unknown scorer {scorer!r}")


def position_index(variants, start0: int, end: int) -> tuple:
    """Row -> position index, and the sorted positions; refuses a variant outside the region or not an SNV."""
    pos = np.empty(len(variants), dtype=np.int64)
    for i, v in enumerate(variants):
        m = VARIANT_RE.match(str(v))
        if m is None:
            raise la.ContractError(f"not a substitution string: {v!r}")
        pos[i] = int(m.group(2))
    if len(pos) and (pos.min() <= start0 or pos.max() > end):
        raise la.ContractError(f"variant outside ({start0}, {end}]: {pos.min()}..{pos.max()}")
    positions, inverse, counts = np.unique(pos, return_inverse=True, return_counts=True)
    if len(counts) and counts.max() > 3:
        raise la.ContractError(f"more than three substitutions at one position ({counts.max()})")
    return inverse, positions


def per_position_max(abs_q: np.ndarray, inverse: np.ndarray, n_positions: int) -> np.ndarray:
    """(n_substitutions, n_tracks) -> (n_positions, n_tracks): max over the substitutions at each position."""
    out = np.full((n_positions, abs_q.shape[1]), -np.inf)
    np.fmax.at(out, inverse, abs_q)          # fmax: a missing quantile never masks a present one
    out[np.isneginf(out)] = np.nan
    return out


def relative_effects(x: np.ndarray, var: dict) -> np.ndarray:
    """|raw effect| divided by each track's own mean level, so tracks of different scale can be averaged.

    The calibrated quantile saturates inside a measured peak (median |quantile| 0.48-0.93 in sampled U1
    regions), so the relative raw effect is the primary scale; a track with no usable mean becomes NaN.
    """
    n = x.shape[1]
    mean = np.asarray([float(m) if str(m) not in ("", "nan", "None") else np.nan
                       for m in var.get("nonzero_mean", [np.nan] * n)], dtype=float)
    mean = np.where(mean > 0, mean, np.nan)
    return np.abs(x) / mean[None, :]


def profile_stats(profile: np.ndarray, positions: np.ndarray, high: float) -> dict:
    """Summaries of one per-position sensitivity profile. Positions are genomic, so gaps are respected."""
    ok = ~np.isnan(profile)
    if not ok.any():
        return {s: float("nan") for s in STATS}
    p, x = positions[ok], profile[ok]
    total = float(x.sum())
    # best 50-bp window by genomic coordinate: for each start position, the sum over positions in [start, start+50)
    csum = np.concatenate([[0.0], np.cumsum(x)])
    right = np.searchsorted(p, p + WINDOW_BP, side="left")
    best = float(np.max(csum[right] - csum[np.arange(len(p))]))
    return {"mean": float(x.mean()), "max": float(x.max()), "sum": total,
            "share_best50": (best / total) if total > 0 else float("nan"),
            "frac_above_high": float((x > high).mean())}


def summarise_scorer(scorer: str, variants, raw: np.ndarray, quantiles: np.ndarray, var: dict,
                     start0: int, end: int) -> tuple:
    """Features on both scales for one scorer of one region; TF top-5 ranked on the relative scale."""
    inverse, positions = position_index(variants, start0, end)
    raw = np.asarray(raw, dtype=float)
    # AVI_SCORE has no track mean level: its raw score is already on one scale, so it is used as is.
    rel = np.abs(raw) if scorer == "AVI_SCORE" else relative_effects(raw, var)
    pmax = {"rel": per_position_max(rel, inverse, len(positions)),
            "q": per_position_max(np.abs(np.asarray(quantiles, dtype=float)), inverse, len(positions))}
    feats, top = {}, []
    for group, cols in group_columns(scorer, var).items():
        feats[f"{group}_n_tracks"] = int(len(cols))
        for scale in SCALES:
            if len(cols) == 0:
                feats.update({f"{group}_{scale}_{s}": float("nan") for s in STATS})
                continue
            sub = pmax[scale][:, cols]
            with np.errstate(invalid="ignore"):
                profile = np.nanmean(sub, axis=1) if np.isfinite(sub).any() else np.full(len(positions), np.nan)
            stats = profile_stats(profile, positions, HIGH[scale])
            feats.update({f"{group}_{scale}_{k}": v for k, v in stats.items()})
        if group.startswith("chip_tf") and len(cols):
            with np.errstate(invalid="ignore"):
                rel_means = np.nanmean(pmax["rel"][:, cols], axis=0)
                q_means = np.nanmean(pmax["q"][:, cols], axis=0)
            ranked = np.argsort(-np.nan_to_num(rel_means, nan=-1.0))
            order = [i for i in ranked if rel_means[i] == rel_means[i]]
            n = len(next(iter(var.values())))
            names = var.get("name", [""] * n)
            factors = var.get("transcription_factor", names)
            for rank, i in enumerate(order[:TOP_TF], start=1):
                c = cols[i]
                top.append({"group": group, "rank": rank, "track_name": str(names[c]),
                            "transcription_factor": str(factors[c]),
                            "mean_rel_effect": float(rel_means[i]),
                            "mean_max_abs_quantile": float(q_means[i])})
    return feats, top

def read_scorer(path: pathlib.Path) -> tuple:
    """(variants, raw effects, quantiles, var-as-dict) from one archived h5ad; refuses one without quantiles."""
    import anndata

    ad = anndata.read_h5ad(path)
    if "quantiles" not in ad.layers:
        raise la.ContractError(f"no quantiles layer in {path}")
    var = {c: ad.var[c].astype(str).tolist() for c in ad.var.columns}
    var["name"] = ad.var["name"].astype(str).tolist() if "name" in ad.var else ad.var_names.astype(str).tolist()
    q = ad.layers["quantiles"]
    q = q.toarray() if hasattr(q, "toarray") else np.asarray(q)
    x = ad.X.toarray() if hasattr(ad.X, "toarray") else np.asarray(ad.X)
    return ad.obs["variant"].astype(str).tolist(), x, q, var


def process_region(task: tuple) -> tuple:
    """One region -> (feature row, TF rows, state). Never raises: a failure is returned as a state."""
    raw_root, key, universe, chrom, start0, end = task
    d = region_dir(pathlib.Path(raw_root), universe, chrom, start0, end)
    row = {"region_key": key, "universe": universe, "chrom": chrom, "start0": start0, "end": end}
    if not (d / "request.json").exists() or not all((d / f"{s}.h5ad").exists() for s in SCORERS):
        return row, [], "archive_incomplete"
    tf_rows = []
    try:
        for scorer in SCORERS:
            variants, x, q, var = read_scorer(d / f"{scorer}.h5ad")
            feats, top = summarise_scorer(scorer, variants, x, q, var, start0, end)
            if scorer == "ATAC":
                row["n_positions"] = len({v.split(":")[1] for v in variants})
            row.update(feats)
            tf_rows += [{"region_key": key, **t} for t in top]
    except Exception as exc:          # recorded per region, never silently dropped
        return row, [], f"failed: {type(exc).__name__}: {exc}"
    return row, tf_rows, "ok"


def main() -> None:
    sat = pathlib.Path(sys.argv[1])
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    tables = la.out_root() / "tables"
    universe = la.read_tsv(sat / "tables" / "saturation_universe.tsv.gz")
    seen, tasks = set(), []
    for u in universe:
        k = (u["universe"], u["chrom"], int(u["start0"]), int(u["end"]))
        if k in seen:
            continue
        seen.add(k)
        tasks.append((str(sat / "raw" / "saturation"), u["region_key"], *k))
    if limit:
        tasks = tasks[:limit]
    la.log(f"step 43: {len(tasks)} regions from {sat}")
    rows, tf_rows, states = [], [], Counter()
    with Pool(workers) as pool:
        for i, (row, top, state) in enumerate(pool.imap(process_region, tasks, chunksize=16), start=1):
            row["archive_state"] = state
            states[state.split(":")[0]] += 1
            rows.append(row)
            tf_rows += top
            if i % 5000 == 0:
                la.log(f"step 43: {i}/{len(tasks)} {dict(states)}")
    cols = ["region_key", "universe", "chrom", "start0", "end", "n_positions", "archive_state"]
    cols += [f"{g}_n_tracks" for g in GROUPS]
    cols += [f"{g}_{sc}_{s}" for g in GROUPS for sc in SCALES for s in STATS]
    la.write_tsv_once(tables / "saturation_region_features.tsv.gz", rows, cols)
    la.write_tsv_once(tables / "saturation_region_tf_top5.tsv.gz", tf_rows,
                      ["region_key", "group", "rank", "track_name", "transcription_factor", "mean_rel_effect",
                       "mean_max_abs_quantile"])
    failures = Counter(r["archive_state"] for r in rows if r["archive_state"].startswith("failed"))
    summary = {"saturation_run": str(sat), "n_regions": len(rows), "archive_states": dict(states),
               "first_failures": dict(failures.most_common(5)), "reduction": "see the module docstring; primary scale |raw effect| / track nonzero_mean, quantile kept",
               "prespec": "43_saturation_prespec.json", "limit": limit}
    json.dump(summary, (tables / "saturation_features_summary.json").open("w"), indent=1)
    la.log(f"step 43: done {dict(states)}")


if __name__ == "__main__":
    main()

"""R08-01 / R09-03 repair: Bonferroni lower bounds for the Figure 7B reporter-assay rows.

Uses the release evaluation code (scripts/evaluate_gse281364_dna_language_seeded_heads.py:
_bootstrap, _macro, _metrics; scripts/evaluate_gse281364_dna_lm_common_lane.py: load_outcomes)
on the release out-of-fold predictions of runs 605 and 604, with the release bootstrap seed
(config uncertainty.seed) and 10,000 resamples of long-range blocks. The 95% intervals are
first checked against the released 7B table; only then are Bonferroni percentile bounds read
from the same draws.

Families:
  F7  the seven displayed 7B reporter-assay rows
  F11 every fitted pretrained head in runs 605 and 604 (displayed rows plus the undisplayed heads)
Bonferroni two-sided: lower bound = quantile alpha/(2m) of the paired gain draws, alpha 0.05.
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import norm

BENCH = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark")
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(BENCH / "src"))
import math  # noqa: E402
from typing import Mapping, Sequence  # noqa: E402
from scipy.stats import pearsonr, rankdata, spearmanr  # noqa: E402

# --- verbatim copy of scripts/evaluate_gse281364_dna_language_seeded_heads.py lines 91-149
# (that module imports masld_bench, which needs Python >= 3.11 tomllib; the rnaseq env is older)
def _macro(values: Sequence[float | None]) -> float | None:
    if any(value is None or not math.isfinite(value) for value in values):
        return None
    clipped = np.clip(np.asarray(values), -0.999999, 0.999999)
    return float(np.tanh(np.mean(np.arctanh(clipped))))


def _metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float | None]:
    residual = observed - predicted
    denominator = float(np.sum((observed - observed.mean()) ** 2))
    values: dict[str, float | None] = {
        "spearman": None,
        "pearson": None,
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "mae": float(np.mean(np.abs(residual))),
        "r2": None if denominator == 0 else float(1.0 - np.sum(residual**2) / denominator),
        "calibration_intercept": None,
        "calibration_slope": None,
    }
    if not np.all(predicted == predicted[0]) and not np.all(observed == observed[0]):
        values["spearman"] = float(spearmanr(observed, predicted).statistic)
        values["pearson"] = float(pearsonr(observed, predicted).statistic)
        slope, intercept = np.polyfit(predicted, observed, 1)
        values["calibration_intercept"] = float(intercept)
        values["calibration_slope"] = float(slope)
    return values


def _fast_spearman(observed: np.ndarray, predicted: np.ndarray) -> float | None:
    if np.all(observed == observed[0]) or np.all(predicted == predicted[0]):
        return None
    left = rankdata(observed)
    right = rankdata(predicted)
    left -= left.mean()
    right -= right.mean()
    denominator = float(np.sqrt(np.sum(left**2) * np.sum(right**2)))
    return None if denominator == 0 else float(np.sum(left * right) / denominator)


def _bootstrap(
    observed: Mapping[str, np.ndarray],
    predictions: Mapping[tuple[str, str], Mapping[str, np.ndarray]],
    blocks: np.ndarray,
    *,
    resamples: int,
    seed: int,
) -> dict[tuple[str, str], np.ndarray]:
    unique = sorted(set(blocks.tolist()))
    by_block = {block: np.flatnonzero(blocks == block) for block in unique}
    output = {candidate: np.full(resamples, np.nan) for candidate in predictions}
    rng = np.random.default_rng(seed)
    for iteration in range(resamples):
        sampled = rng.integers(0, len(unique), size=len(unique))
        indices = np.concatenate([by_block[unique[index]] for index in sampled])
        for candidate, contexts in predictions.items():
            value = _macro([_fast_spearman(observed[context][indices], contexts[context][indices]) for context in CONTEXTS])
            if value is not None:
                output[candidate][iteration] = value
    return output
# --- end verbatim copy

from scripts.evaluate_gse281364_dna_lm_common_lane import CONTEXTS, load_outcomes  # noqa: E402


OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
cfg = json.load(open(BENCH / "config/gse281364_dna_language_seeded_head_campaign.json"))
SEED = int(cfg["uncertainty"]["seed"])
B = int(cfg["uncertainty"]["resamples"])
SEEDS = tuple(int(s) for s in cfg["task"]["seeds"])
FITS = {"605": BENCH / "executions/model-cpu-train-605-21099008/fit/oof_predictions.tsv.gz",
        "604": BENCH / "executions/model-cpu-train-604-21097189/fit/oof_predictions.tsv.gz"}
REF = ("available_simple_controls", "allele_identity_ridge")
CONTROLS = {"zero", "outer_training_mean", "allele_identity_ridge"}

vals, blocks, heads_by_run, ref604 = {}, {}, {}, {}
for run, path in FITS.items():
    heads = Counter()
    with gzip.open(path, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["assay_context_id"] not in CONTEXTS:
                continue
            k = (r["model_id"], r["head_id"])
            if run == "604" and k[0] == "available_simple_controls":
                if k == REF:
                    ref604.setdefault((int(r["seed"]), r["assay_context_id"]), {})[r["element_id"]] = float(r["prediction"])
                continue  # controls are taken from run 605; the 604 reference is checked equal below
            heads[k] += 1
            vals.setdefault((k, int(r["seed"]), r["assay_context_id"]), {})[r["element_id"]] = float(r["prediction"])
            prev = blocks.setdefault(r["element_id"], r["long_range_block_id"])
            assert prev == r["long_range_block_id"], "block map differs between runs"
    heads_by_run[run] = sorted(heads)
elements = sorted(blocks)
ref_maxdiff = max(abs(ref604[(s, c)][e] - vals[(REF, s, c)][e]) for (s, c) in ref604 for e in ref604[(s, c)])
print("run 604 vs 605 allele-identity ridge predictions, max |difference|:", ref_maxdiff, "rows", sum(len(v) for v in ref604.values()))
blk = np.asarray([blocks[e] for e in elements])
print("elements", len(elements), "blocks", len(set(blk.tolist())), "seeds", SEEDS, "boot seed", SEED, "B", B)
assert len(elements) == 1033 and len(set(blk.tolist())) == 239

outcomes = load_outcomes(BENCH / "executions/gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz",
                         set(elements))
observed = {c: np.asarray([outcomes[(e, c)]["mean"] for e in elements]) for c in CONTEXTS}

cands = [k for run in heads_by_run for k in heads_by_run[run]
         if k[0] != "available_simple_controls" or k == REF]
ens = {k: {c: np.mean(np.vstack([[vals[(k, s, c)][e] for e in elements] for s in SEEDS]), axis=0)
           for c in CONTEXTS} for k in cands}
point = {k: _macro([_metrics(observed[c], ens[k][c])["spearman"] for c in CONTEXTS]) for k in cands}
print("reference point", point[REF], "(receipt 0.13421551707194201)")

draws = _bootstrap(observed, ens, blk, resamples=B, seed=SEED)

# Released 7B rows for the reproduction check.
LABEL = {("caduceus", "delta_ridge"): "caduceus delta ridge",
         ("dnabert2", "delta_ridge"): "dnabert2 delta ridge",
         ("hyenadna", "delta_ridge"): "hyenadna delta ridge",
         ("nucleotide_transformer", "delta_ridge"): "nucleotide_transformer delta ridge [restricted comparator]",
         ("enformer_crested_restricted_port", "rc_ensemble_sad_ridge"): "enformer (restricted port) SAD ridge",
         ("enformer_crested_restricted_port", "rc_ensemble_sar_ridge"): "enformer (restricted port) SAR ridge",
         ("sei", "native_40class_ridge"): "sei native 40-class ridge"}
rel = {r["model"]: r for r in csv.DictReader(open(BENCH / "executions/fig7-source-exposure-correction-20261014/tables/panel_7b_gains.tsv"), delimiter="\t")
       if r["task"].startswith("MPRA")}
F7 = list(LABEL)
F11 = [k for k in cands if k != REF]
print("displayed rows", len(F7), "all pretrained heads", len(F11), F11)

rows = []
for fam_name, fam in (("F7_displayed_rows", F7), (f"F{len(F11)}_all_fitted_heads", F11)):
    m = len(fam)
    a = 0.05 / m
    z = norm.ppf(1 - a / 2)
    for k in fam:
        g = draws[k] - draws[REF]
        g = g[np.isfinite(g)]
        G = point[k] - point[REF]
        lo95, hi95, blo, bhi = np.quantile(g, [0.025, 0.975, a / 2, 1 - a / 2])
        r = {"family": fam_name, "m": m, "model": k[0], "head": k[1], "gain": G,
             "lo95": lo95, "hi95": hi95, "bonf_pct_lo": blo, "bonf_pct_hi": bhi,
             "bonf_normal_lo": G - z * np.std(g, ddof=1), "frac_draws_le0": float(np.mean(g <= 0)),
             "n_valid": len(g)}
        if k in LABEL:
            rr = rel[LABEL[k]]
            r["release_gain_diff"] = G - float(rr["gain"])
            r["release_lo_diff"] = lo95 - float(rr["lo"])
            r["release_hi_diff"] = hi95 - float(rr["hi"])
        rows.append(r)

fields = ["family", "m", "model", "head", "gain", "lo95", "hi95", "bonf_pct_lo", "bonf_pct_hi",
          "bonf_normal_lo", "frac_draws_le0", "n_valid", "release_gain_diff", "release_lo_diff", "release_hi_diff"]
with open(OUT / "mpra_bonferroni.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})

repro = max(abs(r.get(x, 0.0)) for r in rows for x in ("release_gain_diff", "release_lo_diff", "release_hi_diff"))
print("max |difference| from released 7B gain/lo/hi:", repro)

# LX-2 contexts: present in the outcome table or not.
ctx = Counter()
with gzip.open(BENCH / "executions/gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz", "rt") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        ctx[(r["context_id"], r["cell_line"])] += 1
print("outcome rows by context:", dict(ctx))
json.dump({"max_abs_diff_from_release": repro, "ref604_vs_605_maxdiff": ref_maxdiff, "outcome_rows_by_context": {f"{a}|{b}": n for (a, b), n in ctx.items()},
           "heads_by_run": {k: [list(x) for x in v] for k, v in heads_by_run.items()}, "boot_seed": SEED, "B": B},
          open(OUT / "mpra_bonferroni_meta.json", "w"), indent=1)

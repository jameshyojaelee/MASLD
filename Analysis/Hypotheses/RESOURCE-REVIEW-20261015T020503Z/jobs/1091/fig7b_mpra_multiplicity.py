#!/usr/bin/env python3
"""Review R09: re-derive the Fig 7B MPRA (GSE281364) gains over the allele-identity ridge and
add family-wise inference.

Reproduces the evaluator in
Analysis/MASLD_Model_Benchmark/scripts/evaluate_gse281364_dna_language_seeded_heads.py
(same OOF predictions, same seed ensemble, same 239 long-range blocks, same bootstrap seed 20260825
and 10,000 resamples), then reports for every DNA-language candidate (4 models x 2 prespecified heads):
  - point gain and percentile 95% interval (must match model_head_summary.tsv),
  - one-sided percentile-bootstrap p for gain <= 0, (k+1)/(N+1),
  - Holm adjustment over the 8 DNA-language candidates and over the 4 displayed delta heads,
  - max-|t| simultaneous 95% intervals over the 8 candidates.
Read-only on all inputs. Writes only to the current directory.
"""
import csv
import gzip
import math
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import rankdata

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/"
OOF = ROOT + "executions/model-cpu-train-605-21099008/fit/oof_predictions.tsv.gz"
OUTC = ROOT + "executions/gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz"
SUMMARY = ROOT + "executions/model-cpu-train-605-21099008/evaluation/model_head_summary.tsv"
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
SEEDS = (1103, 2909, 4721, 6673, 8111)
PSEUDO = 0.5
RESAMPLES = 10000
BOOT_SEED = 20260825

# Predictions
pred = defaultdict(dict)  # (model, head, seed, context) -> {element: value}
block_of = {}
with gzip.open(OOF, "rt") as h:
    for r in csv.DictReader(h, delimiter="\t"):
        key = (r["model_id"], r["head_id"], int(r["seed"]), r["assay_context_id"])
        pred[key][r["element_id"]] = float(r["prediction"])
        b = block_of.setdefault(r["element_id"], r["long_range_block_id"])
        assert b == r["long_range_block_id"]
elements = sorted(block_of)
blocks = np.asarray([block_of[e] for e in elements])
print("elements", len(elements), "blocks", len(set(blocks.tolist())))
candidates = sorted({(k[0], k[1]) for k in pred})
print("candidates", candidates)

# Outcomes: mean over 4 replicates of log2 alt activity minus log2 ref activity
pairs = defaultdict(dict)
eset = set(elements)
with gzip.open(OUTC, "rt") as h:
    for r in csv.DictReader(h, delimiter="\t"):
        if r["element_id"] not in eset or r["context_id"] not in CONTEXTS or r["assay_state"] != "observed":
            continue
        pairs[(r["element_id"], r["context_id"], int(r["experimental_replicate"]))][r["allele"]] = (int(r["DNA"]), int(r["RNA"]))
observed = {}
for c in CONTEXTS:
    vals = []
    for e in elements:
        reps = []
        for rep in range(1, 5):
            p = pairs[(e, c, rep)]
            (rd, rr), (ad, ar) = p["ref"], p["alt"]
            reps.append(math.log2((ar + PSEUDO) / (ad + PSEUDO)) - math.log2((rr + PSEUDO) / (rd + PSEUDO)))
        vals.append(np.mean(reps))
    observed[c] = np.asarray(vals)

ens = {}
for m, hd in candidates:
    ens[(m, hd)] = {c: np.mean(np.vstack([np.asarray([pred[(m, hd, s, c)][e] for e in elements]) for s in SEEDS]), axis=0) for c in CONTEXTS}


def spear(o, p):
    if np.all(o == o[0]) or np.all(p == p[0]):
        return None
    a = rankdata(o); b = rankdata(p)
    a -= a.mean(); b -= b.mean()
    d = math.sqrt(float(np.sum(a * a) * np.sum(b * b)))
    return None if d == 0 else float(np.sum(a * b) / d)


def macro(vals):
    if any(v is None or not math.isfinite(v) for v in vals):
        return None
    x = np.clip(np.asarray(vals), -0.999999, 0.999999)
    return float(np.tanh(np.mean(np.arctanh(x))))


points = {k: macro([spear(observed[c], ens[k][c]) for c in CONTEXTS]) for k in candidates}
ref = ("available_simple_controls", "allele_identity_ridge")
print("reference point", points[ref])

unique = sorted(set(blocks.tolist()))
by_block = {b: np.flatnonzero(blocks == b) for b in unique}
rng = np.random.default_rng(BOOT_SEED)
samples = {k: np.full(RESAMPLES, np.nan) for k in candidates}
for it in range(RESAMPLES):
    sampled = rng.integers(0, len(unique), size=len(unique))
    idx = np.concatenate([by_block[unique[i]] for i in sampled])
    for k in candidates:
        v = macro([spear(observed[c][idx], ens[k][c][idx]) for c in CONTEXTS])
        if v is not None:
            samples[k][it] = v

dna = [k for k in candidates if k[0] != "available_simple_controls"]
gain = {k: points[k] - points[ref] for k in dna}
gs = {k: samples[k] - samples[ref] for k in dna}

# Stored values for reproduction check
stored = {}
with open(SUMMARY) as h:
    for r in csv.DictReader(h, delimiter="\t"):
        stored[(r["model_id"], r["head_id"])] = r

rows = []
for k in dna:
    g = gs[k][np.isfinite(gs[k])]
    lo, hi = np.quantile(g, (0.025, 0.975))
    p_one = (int(np.sum(g <= 0)) + 1) / (len(g) + 1)
    sd = float(np.std(g, ddof=1))
    rows.append(dict(model=k[0], head=k[1], gain=gain[k], lo=lo, hi=hi, sd=sd, p_one_sided=p_one,
                     stored_gain=float(stored[k]["gain_vs_available_control"]),
                     stored_lo=float(stored[k]["gain_ci_low"]), stored_hi=float(stored[k]["gain_ci_high"])))


def holm(ps):
    order = np.argsort(ps)
    m = len(ps)
    adj = np.empty(m)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (m - rank) * ps[i]))
        adj[i] = run
    return adj


p_all = np.asarray([r["p_one_sided"] for r in rows])
for r, a in zip(rows, holm(p_all)):
    r["holm_8_dna_lm"] = a
delta_idx = [i for i, r in enumerate(rows) if r["head"] == "delta_ridge"]
adj_d = holm(p_all[delta_idx])
for i, a in zip(delta_idx, adj_d):
    rows[i]["holm_4_delta"] = a

# max-|t| simultaneous band over the 8 DNA-LM candidates
mat = np.vstack([gs[k] for k in dna])
ok = np.all(np.isfinite(mat), axis=0)
mat = mat[:, ok]
centre = np.asarray([gain[k] for k in dna])[:, None]
sds = np.asarray([r["sd"] for r in rows])[:, None]
tmax = np.max(np.abs((mat - centre) / sds), axis=0)
crit = float(np.quantile(tmax, 0.95))
for r in rows:
    r["simul_lo"] = r["gain"] - crit * r["sd"]
    r["simul_hi"] = r["gain"] + crit * r["sd"]
print("max-|t| critical value over 8 candidates:", round(crit, 4), "valid draws", int(ok.sum()))

fields = ["model", "head", "gain", "lo", "hi", "stored_gain", "stored_lo", "stored_hi", "sd", "p_one_sided",
          "holm_8_dna_lm", "holm_4_delta", "simul_lo", "simul_hi"]
with open("fig7b_mpra_multiplicity.tsv", "w") as h:
    w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow({f: (r.get(f, "") if isinstance(r.get(f, ""), str) else format(r.get(f, float("nan")), ".6g")) for f in fields})
for r in rows:
    print("\t".join(str(r.get(f, "")) if isinstance(r.get(f, ""), str) else format(r.get(f, float("nan")), ".5g") for f in fields))
print("max abs diff gain vs stored:", max(abs(r["gain"] - r["stored_gain"]) for r in rows))
print("max abs diff lo/hi vs stored:", max(max(abs(r["lo"] - r["stored_lo"]), abs(r["hi"] - r["stored_hi"])) for r in rows))

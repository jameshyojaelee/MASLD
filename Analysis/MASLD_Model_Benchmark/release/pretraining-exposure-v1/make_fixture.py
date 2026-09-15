#!/usr/bin/env python3
"""Build a synthetic corpus and query with KNOWN ground truth, so the test can fail.

A release test that only re-runs the tool on the data the tool was written against proves
nothing. This builds a corpus where we planted the answer, including the cases the method
is supposed to REFUSE, and `test_release.sh` checks the tool against what was planted.

THE GENERATIVE MODEL. Every sample is

    x = mu0  +  study_latent . S  +  sample_latent . F  +  eps

with the variance budget chosen so the synthetic corpus has the property that makes this
problem hard in real data: two different samples of the SAME study already correlate at
about 0.88, and two samples of different studies at about 0.80. A method that decides
membership by "is the correlation high" cannot work on that, which is the point.

    var(mu0) = 0.80   shared across every sample in the corpus
    var(study) = 0.08 shared within a study
    var(sample) = 0.10
    var(eps) = 0.02

THE PIPELINE PAIR. The corpus is written through "pipeline B" -- a per-gene systematic
offset plus per-sample white noise -- while the query stays in "pipeline A". That is the
real situation: the corpus was quantified by the model's authors and the query by whoever
is auditing. The `crosspipeline` scenario tunes it to the magnitude measured in the
reference application (same-sample band r in [0.9835, 0.9934]).

ONE ATTENUATION LEVEL PER SCENARIO. The band is a property of the pipeline PAIR, so every
cohort in one scan must share it. Mixing a verbatim cohort with an attenuated one in a
single scan is not a configuration that occurs, and the band would wrongly reject the
attenuated one. Cohorts quantified differently must be scanned separately.

THE SIX COHORTS, and which leg each one exercises:

  planted_contiguous        verbatim copy, deposit order = query order, one run.
                            Both legs pass. Calibrates the band.
  planted_shuffled_order    verbatim copy, one run, deposit order scrambled.
                            Both legs pass; the offset run does NOT, which is why leg 1
                            is order-free.
  planted_scattered_rows    verbatim copy at non-contiguous rows across the corpus.
                            Injectivity passes, the WINDOW fails -> the two halves of
                            leg 1 are not redundant.
  paired_sibling_scrambled  different samples, 1:1 paired to the query by a shared
                            sample-level latent, deposited scrambled. Leg 1 passes on
                            both halves; only MAGNITUDE rejects it -> leg 2 is not
                            redundant with leg 1.
  similar_only              different samples of the same study, unpaired. The plain
                            false-positive guard: nothing of this cohort is in the corpus.
  absent                    nothing of this cohort, and no sibling, is in the corpus.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

# Variance budget. Same-study r ~= (V0+VS)/1.0, different-study r ~= V0/1.0.
V0, VS, VF, VE = 0.80, 0.08, 0.10, 0.02
K_STUDY, K_SAMPLE = 6, 12
#: Correlation between a query sample's latent and its 1:1 sibling's latent. Sets the
#: paired sibling to ~0.97 -- above every other neighbour, below any same-sample band.
SIBLING_RHO = 0.90

#: Per-sample spread of the attenuation. Real cross-pipeline loss is not the same for
#: every library, so the band has width; this reproduces the width measured in the
#: reference application rather than giving every sample the identical noise level.
NOISE_LOGSD = 0.13

SCENARIOS = {
    # name           -> (per-gene offset sd, per-sample white noise sd)
    # `crosspipeline` is tuned so the same-sample band lands on the magnitude measured in
    # the reference application: r in [0.9835, 0.9934], median 0.9895.
    "verbatim":       (0.0, 0.0),
    "crosspipeline":  (0.030, 0.1430),
    "nothing_planted": (0.0, 0.0),
}

COHORTS = [
    "planted_contiguous", "planted_shuffled_order", "planted_scattered_rows",
    "paired_sibling_scrambled", "similar_only", "absent",
]


def _basis(rng, n_genes):
    return (rng.normal(0.0, 1.0, n_genes) * np.sqrt(V0),
            rng.normal(0.0, 1.0, (K_STUDY, n_genes)) * np.sqrt(VS / K_STUDY),
            rng.normal(0.0, 1.0, (K_SAMPLE, n_genes)) * np.sqrt(VF / K_SAMPLE))


def _draw(rng, mu0, S, F, study_latent, sample_latents):
    X = mu0[None, :] + (study_latent @ S)[None, :] + sample_latents @ F
    return X + rng.normal(0.0, np.sqrt(VE), X.shape)


def build(scenario: str, out_dir: Path, seed: int = 1103, n_genes: int = 1200,
          n_per_cohort: int = 60, n_corpus: int = 6000) -> dict:
    if scenario not in SCENARIOS:
        raise SystemExit(f"unknown scenario {scenario}; choose from {sorted(SCENARIOS)}")
    offset_sd, noise_sd = SCENARIOS[scenario]
    planted = scenario != "nothing_planted"
    rng = np.random.default_rng(seed)
    out_dir.mkdir(parents=True, exist_ok=True)
    genes = [f"G{i:06d}" for i in range(n_genes)]
    mu0, S, F = _basis(rng, n_genes)

    # ---- the query: six cohorts, each its own study ---------------------------------
    q_parts, q_lat, q_study = {}, {}, {}
    for c in COHORTS:
        sl = rng.normal(0.0, 1.0, K_STUDY)
        fl = rng.normal(0.0, 1.0, (n_per_cohort, K_SAMPLE))
        q_study[c], q_lat[c] = sl, fl
        q_parts[c] = _draw(rng, mu0, S, F, sl, fl)
    query = np.concatenate([q_parts[c] for c in COHORTS]).astype(np.float32)
    cohort_of = [c for c in COHORTS for _ in range(n_per_cohort)]
    sample_ids = [f"{c}_s{i:03d}" for c in COHORTS for i in range(n_per_cohort)]

    # ---- the corpus: filler studies, then the planted blocks written over them -------
    corpus = np.empty((n_corpus, n_genes), dtype=np.float64)
    at = 0
    while at < n_corpus:
        n = int(min(rng.integers(30, 121), n_corpus - at))
        corpus[at:at + n] = _draw(rng, mu0, S, F, rng.normal(0, 1, K_STUDY),
                                  rng.normal(0, 1, (n, K_SAMPLE)))
        at += n

    truth = {c: {"planted_rows": None, "deposit_order_matches_query": False}
             for c in COHORTS}
    starts = {"planted_contiguous": 1200, "planted_shuffled_order": 3000,
              "paired_sibling_scrambled": 4200, "similar_only": 5000}

    if planted:
        # verbatim copies, deposit order preserved
        s = starts["planted_contiguous"]
        corpus[s:s + n_per_cohort] = q_parts["planted_contiguous"]
        truth["planted_contiguous"].update(
            planted_rows=[s, s + n_per_cohort - 1], deposit_order_matches_query=True,
            constant_offset=s)

        # verbatim copy, one run, deposit order scrambled
        s = starts["planted_shuffled_order"]
        perm = rng.permutation(n_per_cohort)
        corpus[s:s + n_per_cohort] = q_parts["planted_shuffled_order"][perm]
        truth["planted_shuffled_order"].update(
            planted_rows=[s, s + n_per_cohort - 1], deposit_permutation=perm.tolist())

        # verbatim copy at scattered, non-contiguous rows
        reserved = set()
        for k, st in starts.items():
            reserved |= set(range(st, st + n_per_cohort))
        free = np.array([i for i in range(n_corpus) if i not in reserved])
        scat = np.sort(rng.choice(free, size=n_per_cohort, replace=False))
        corpus[scat] = q_parts["planted_scattered_rows"]
        truth["planted_scattered_rows"].update(planted_rows=scat.tolist())

        # 1:1-paired different samples, deposited scrambled
        s = starts["paired_sibling_scrambled"]
        base = q_lat["paired_sibling_scrambled"]
        sib = SIBLING_RHO * base + np.sqrt(1 - SIBLING_RHO ** 2) * \
            rng.normal(0, 1, base.shape)
        perm = rng.permutation(n_per_cohort)
        corpus[s:s + n_per_cohort] = _draw(
            rng, mu0, S, F, q_study["paired_sibling_scrambled"], sib)[perm]
        truth["paired_sibling_scrambled"].update(
            planted_rows=[s, s + n_per_cohort - 1], sibling_rho=SIBLING_RHO)

        # unpaired different samples of the same study
        s = starts["similar_only"]
        corpus[s:s + n_per_cohort] = _draw(
            rng, mu0, S, F, q_study["similar_only"],
            rng.normal(0, 1, (n_per_cohort, K_SAMPLE)))
        truth["similar_only"].update(sibling_block=[s, s + n_per_cohort - 1])

    # ---- pipeline B: the corpus's quantification is not the query's -----------------
    if offset_sd > 0:
        corpus += rng.normal(0.0, offset_sd, n_genes)[None, :]
    if noise_sd > 0:
        scale = np.exp(rng.normal(0.0, NOISE_LOGSD, n_corpus))[:, None]
        corpus += rng.normal(0.0, noise_sd, corpus.shape) * scale

    np.save(out_dir / "corpus.npy", corpus.astype(np.float32))
    (out_dir / "corpus_genes.txt").write_text("\n".join(genes) + "\n")
    with open(out_dir / "query.tsv", "w") as fh:      # rows = genes, cols = samples
        fh.write("gene\t" + "\t".join(sample_ids) + "\n")
        for g in range(n_genes):
            fh.write(genes[g] + "\t" + "\t".join(f"{v:.6g}" for v in query[:, g]) + "\n")
    with open(out_dir / "cohorts.tsv", "w") as fh:
        fh.write("sample_id\tcohort\n")
        for sid, c in zip(sample_ids, cohort_of):
            fh.write(f"{sid}\t{c}\n")

    expect = _expected_states(scenario)
    meta = {
        "scenario": scenario, "seed": seed, "n_genes": n_genes,
        "n_per_cohort": n_per_cohort, "n_corpus": n_corpus,
        "pipeline_b": {"per_gene_offset_sd": offset_sd, "white_noise_sd": noise_sd,
                       "per_sample_noise_logsd": NOISE_LOGSD if noise_sd else 0.0},
        "variance_budget": {"mu0": V0, "study": VS, "sample": VF, "eps": VE},
        "implied_same_study_r": round((V0 + VS) / (V0 + VS + VF + VE), 4),
        "implied_diff_study_r": round(V0 / (V0 + VS + VF + VE), 4),
        "cohorts": {c: truth[c] | expect[c] for c in COHORTS},
    }
    (out_dir / "truth.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def _expected_states(scenario: str) -> dict:
    if scenario == "nothing_planted":
        return {c: {"expected_state": "unknown", "expected_leg1": False,
                    "expected_leg2": False,
                    "why": "nothing of this cohort is in the corpus"} for c in COHORTS}
    return {
        "planted_contiguous": {
            "expected_state": "encoder_seen", "expected_leg1": True, "expected_leg2": True,
            "why": "verbatim contiguous copy in deposit order: both legs, exact offset"},
        "planted_shuffled_order": {
            "expected_state": "encoder_seen", "expected_leg1": True, "expected_leg2": True,
            "why": "verbatim copy in one run, order scrambled: leg 1 is order-free"},
        "planted_scattered_rows": {
            "expected_state": "unknown", "expected_leg1": False, "expected_leg2": True,
            "why": "injectivity passes, the window fails: leg 1's two halves differ"},
        "paired_sibling_scrambled": {
            "expected_state": "unknown", "expected_leg1": True, "expected_leg2": False,
            "why": "structure passes, magnitude rejects: leg 2 is not redundant"},
        "similar_only": {
            "expected_state": "unknown", "expected_leg1": False, "expected_leg2": False,
            "why": "merely similar: the false-positive guard"},
        "absent": {
            "expected_state": "unknown", "expected_leg1": False, "expected_leg2": False,
            "why": "not in the corpus at all"},
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--scenario", required=True, choices=sorted(SCENARIOS))
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=1103)
    p.add_argument("--n-genes", type=int, default=1200)
    p.add_argument("--n-per-cohort", type=int, default=60)
    p.add_argument("--n-corpus", type=int, default=6000)
    a = p.parse_args()
    m = build(a.scenario, Path(a.out), a.seed, a.n_genes, a.n_per_cohort, a.n_corpus)
    print(f"  fixture '{a.scenario}' -> {a.out}  "
          f"(same-study r ~ {m['implied_same_study_r']}, "
          f"diff-study r ~ {m['implied_diff_study_r']})")

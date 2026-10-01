"""Outcome-blind Gate A simulator smoke check against the v0 observation model.

This samples entirely synthetic genotypes and reads. Its cohort depth profiles
come from the public aggregate T4 depth summary. It makes no Gate A decision;
the final design simulation needs Gate 0's cohort × gene inputs and v1 rules.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEPTH = ROOT / "Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-tags-v2-20260929T134419Z/tag_depth_summary_by_cohort.tsv"
AFC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/public_afc_summary.json"
QUANTILES = ((0.0, "q0.1"), (0.10, "q0.1"), (0.25, "q0.25"), (0.50, "q0.5"),
             (0.75, "q0.75"), (0.90, "q0.9"), (0.95, "q0.95"), (0.99, "q0.99"),
             (0.999, "q0.999"), (1.0, "max_n"))


def read_profiles():
    with DEPTH.open() as handle:
        profiles = list(csv.DictReader(handle, delimiter="\t"))
    if len(profiles) != 10 or len({p["cohort"] for p in profiles}) != 10:
        raise ValueError("expected ten unique source cohorts")
    return profiles


def sample_depth(profile, rng, size):
    probabilities = np.array([x[0] for x in QUANTILES])
    depths = np.array([float(profile[x[1]]) for x in QUANTILES])
    if np.any(np.diff(depths) < 0):
        raise ValueError(f"depth quantiles are not monotone: {profile['cohort']}")
    # Interpolate on log depth, preserving the long upper tail without
    # treating the cohort maximum as a typical donor.
    draws = np.exp(np.interp(rng.random(size), probabilities, np.log(depths)))
    return np.maximum(1, np.rint(draws).astype(np.int32))


def beta_binomial(rng, n, mean, rho):
    total = (1 - rho) / rho
    probability = rng.beta(mean * total, (1 - mean) * total, size=n.shape)
    return rng.binomial(n, probability)


def simulate(profile, rng, genes, donors, alpha, kappa, tau, min_depth=8,
             min_minor=2, min_fraction=0.05):
    # Synthetic q, e, rho, omega, and balanced stage are provisional; Gate 0
    # will replace them with their prespecified aggregate distributions.
    n = sample_depth(profile, rng, (genes, donors))
    stage = np.broadcast_to(np.arange(donors) % 3, (genes, donors)) - 1
    q = 0.2
    genotype = rng.choice(3, size=n.shape, p=(2*q*(1-q), (1-q)**2, q**2))
    effect = np.full((genes, 1), alpha) * rng.choice((-1, 1), size=(genes, 1))
    deviation = rng.normal(0, tau, size=(genes, 1))
    eta = effect * (1 + kappa * stage) + deviation * stage - 0.02
    mean_het = 1 / (1 + np.exp(-eta))
    mean = np.where(genotype == 0, mean_het, np.where(genotype == 1, 0.01, 0.99))
    rho = np.where(genotype == 0, 0.02, 0.01)
    a = beta_binomial(rng, n, mean, rho)
    lower = np.maximum(min_minor, np.ceil(min_fraction * n)).astype(int)
    included = (n >= min_depth) & (a >= lower) & (n - a >= lower)
    return {
        "candidate_rows": int(n.size),
        "depth_pass": int(np.count_nonzero(n >= min_depth)),
        "called_het_rows": int(np.count_nonzero(included)),
        "true_het_among_called": int(np.count_nonzero(included & (genotype == 0))),
        "stage_called_het": [int(np.count_nonzero(included & (stage == s))) for s in (-1, 0, 1)],
        "genes_with_called_het": int(np.count_nonzero(included.any(axis=1))),
        "synthetic_median_depth": float(np.median(n)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", type=int, default=64)
    parser.add_argument("--donors-per-cohort", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20260929)
    args = parser.parse_args()
    if args.genes < 1 or args.donors_per_cohort < 3:
        raise ValueError("need positive genes and at least three donors")
    afc = json.loads(AFC.read_text())
    profiles = read_profiles()
    result = {"status": "synthetic plumbing check only; no Gate A decision",
              "seed": args.seed, "genes_per_cohort": args.genes,
              "donors_per_cohort": args.donors_per_cohort, "scenarios": []}
    for item in afc["thresholds"]:
        theta = item["theta"]
        alpha = item["median_abs_ln_afc"]
        tau_max = 0.125 * alpha
        for kappa, tau in ((0, 0), (-0.125, 0), (0.125, 0), (0, tau_max),
                           (-0.125, tau_max), (0.125, tau_max)):
            rng = np.random.default_rng(args.seed + int(theta * 1000) * 10
                                        + int((kappa + 0.125) * 1000)
                                        + int(tau > 0) * 10000)
            cohorts = {p["cohort"]: simulate(p, rng, args.genes,
                       args.donors_per_cohort, alpha, kappa, tau) for p in profiles}
            result["scenarios"].append({"theta": theta, "eligible_genes": item["eligible_genes"],
                     "alpha_median": alpha, "kappa": kappa, "tau": tau,
                     "cohorts": cohorts})
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

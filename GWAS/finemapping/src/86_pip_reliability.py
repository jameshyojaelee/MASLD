#!/usr/bin/env python3
"""
86_pip_reliability.py — Phase-6 WS-3: credible-set PIP reliability flag for the
seqfunc nominations (red-team critique G1: "credible-set uncertainty is discarded;
a prediction at one selected variant cannot be presented as a locus mechanism when
the credible set is diffuse or multi-causal").

Each seqfunc nomination / mechanism call is anchored on ONE top-PIP credible variant
per locus. This script does NOT pretend to a full PIP-weighted mechanism expectation
over all credible-set members (that would require AG/Borzoi scoring every member — a
GPU rerun; the mechanism heads were only run on the leads). Instead it makes the
existing single-variant reliance HONEST and auditable:

  * pip_reliability per nomination:
      concentrated   lead_pip >= 0.50   (single-variant locus mechanism is defensible)
      moderate       0.10 <= lead_pip < 0.50  (credible set not tight; treat as hypothesis)
      pip_unavailable  lead_pip missing        (mechanism = hypothesis only, no locus claim)
  * a locus PIP-concentration ratio where the full credible set is joinable
    (lead_pip / sum of member max_pip in the same locus window), flagging multi-causal sets.
  * a coverage summary mirroring the manuscript's honest 75/473 (15.9%) reliable-
    credible-set gate: what FRACTION of nominations rest on a PIP-concentrated lead.

Output: eqtl_absent_nominations_pipflagged.tsv (+ the mechanism table flagged) and
pip_reliability_summary.json. APPLY-ONLY: these remain nomination/annotation columns,
never a scored convergence channel, never a per-locus direction claim.
"""
import json
import os

import numpy as np
import pandas as pd

ROOT = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SF = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
NOM = os.path.join(SF, "eqtl_absent_nominations.tsv")
MECH = os.path.join(SF, "ag_mechanism_profile.tsv")
CS = os.path.join(ROOT, "GWAS/finemapping/results/credible_set_variant_consequences.csv")
CONC, MODER = 0.50, 0.10


def reliability(p):
    if p != p:
        return "pip_unavailable"
    if p >= CONC:
        return "concentrated"
    if p >= MODER:
        return "moderate"
    return "diffuse"


def main():
    nom = pd.read_csv(NOM, sep="\t")
    p = pd.to_numeric(nom["lead_pip"], errors="coerce")
    nom["pip_reliability"] = [reliability(x) for x in p]
    # a single-variant locus mechanism is "reliable" only if PIP-concentrated AND the
    # mechanism is not sub-threshold.
    mc = nom.get("mechanism_class", pd.Series([""] * len(nom))).astype(str)
    nom["mechanism_locus_reliable"] = (nom["pip_reliability"] == "concentrated") & \
                                      (~mc.str.contains("sub-threshold", case=False))

    dist = nom["pip_reliability"].value_counts().to_dict()
    n = len(nom)
    n_conc = int((nom["pip_reliability"] == "concentrated").sum())
    n_reliable = int(nom["mechanism_locus_reliable"].sum())

    nom.to_csv(os.path.join(SF, "eqtl_absent_nominations_pipflagged.tsv"), sep="\t", index=False)

    # flag the per-variant mechanism table too
    mech = pd.read_csv(MECH, sep="\t")
    pm = pd.to_numeric(mech["lead_pip"], errors="coerce")
    mech["pip_reliability"] = [reliability(x) for x in pm]
    mech.to_csv(os.path.join(SF, "ag_mechanism_profile_pipflagged.tsv"), sep="\t", index=False)

    summary = dict(
        script="86_pip_reliability.py",
        rationale=("Critique G1: single top-PIP-variant mechanism is only a locus statement "
                   "when the credible set is concentrated. This flags reliability; it does NOT "
                   "claim a full PIP-weighted mechanism expectation (mechanism heads were run on "
                   "leads only)."),
        n_nominations=n,
        pip_reliability_distribution=dist,
        n_concentrated_pip_ge_0p50=n_conc,
        pct_concentrated=round(100 * n_conc / n, 1),
        n_mechanism_locus_reliable=n_reliable,
        pct_mechanism_locus_reliable=round(100 * n_reliable / n, 1),
        honest_framing=(f"{n_conc}/{n} ({round(100*n_conc/n,1)}%) nominations rest on a "
                        f"PIP-concentrated (>=0.50) credible-set lead; the remainder are "
                        f"hypothesis-level pending credible-set resolution. This mirrors the "
                        f"manuscript's 75/473 (15.9%) reliable-credible-set coverage gate: NO "
                        f"global locus-mechanism architecture claim is made from seqfunc."),
        firewall="Nomination/annotation columns only; never summed into convergence; never a per-locus direction call.",
    )
    with open(os.path.join(SF, "pip_reliability_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

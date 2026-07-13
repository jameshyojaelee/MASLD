#!/usr/bin/env python3
"""
91_pip_weighted_fine_mechanism.py — Phase-6 WS-3 Tier-2 integration: within-credible-set
PIP-weighted FINE regulatory sub-mechanism (completes critique G1 at the fine layer).

Tier-1 (src/89) showed the COARSE mechanism class (coding/splice/regulatory) is robust to
credible-set uncertainty. The diffuse credible sets are diffuse in WHICH VARIANT but mostly
100%-regulatory in CLASS. This step asks the finer question: within a diffuse credible set, is
the regulatory SUB-mechanism (enhancer / promoter / TF-footprint / splice) CONSISTENT across
members, or does it vary? If consistent -> the sub-mechanism is a defensible locus statement
despite variant uncertainty; if it varies -> even the sub-mechanism is a hypothesis.

INPUTS:
  ag_member_finemechanism.tsv  (src/90; AlphaGenome zero-shot fine_mech for ~505 diffuse-CS members)
  ag_mechanism_profile.tsv     (src/72; mechanism_class for the ~64 already-scored leads)
  nominated_cs_members.tsv     (within-CS susie_pip)
  pip_weighted_mechanism_per_cs.tsv (src/89; coarse result + concentration)

Per credible set: PIP-weighted fine-mechanism posterior mass, scored-coverage (fraction of the
CS posterior on scored members), dominant fine mechanism + mass, and whether it is defensible
(dominant fine mass >= 0.5 AND coverage >= 0.5).

FRAMING: the fine taxonomy is a hand-built, unvalidated classifier (critique M6) → HYPOTHESIS-
LEVEL. AlphaGenome zero-shot ANNOTATION only (never trained on; non-commercial). APPLY-ONLY.
"""
import json
import os

import numpy as np
import pandas as pd

ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SF = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
MEMBERS = os.path.join(SF, "nominated_cs_members.tsv")
AGMEM = os.path.join(SF, "ag_member_finemechanism.tsv")
AGLEAD = os.path.join(SF, "ag_mechanism_profile.tsv")
COARSE = os.path.join(SF, "pip_weighted_mechanism_per_cs.tsv")
OUT_TSV = os.path.join(SF, "pip_weighted_fine_mechanism_per_cs.tsv")
OUT_JSON = os.path.join(SF, "pip_weighted_fine_mechanism_summary.json")

FINE = ["enhancer-disrupting", "promoter/TSS", "splice-altering",
        "TF-footprint-breaking", "3D-contact-rewiring", "sub-threshold"]


def main():
    mem = pd.read_csv(MEMBERS, sep="\t")
    mem["susie_pip"] = pd.to_numeric(mem["susie_pip"], errors="coerce").fillna(0.0)

    # fine mechanism per member: prefer the src/90 member score, fall back to the lead score
    fm = pd.read_csv(AGMEM, sep="\t")[["variant_id", "fine_mech"]]
    fm = fm[fm["fine_mech"] != "ERROR"].drop_duplicates("variant_id")
    fine = dict(zip(fm["variant_id"].astype(str), fm["fine_mech"].astype(str)))
    lead = pd.read_csv(AGLEAD, sep="\t")[["variant_id_hg19", "mechanism_class"]].dropna()
    for _, r in lead.iterrows():
        fine.setdefault(str(r["variant_id_hg19"]), str(r["mechanism_class"]))

    mem["fine_mech"] = mem["variant_id"].astype(str).map(fine).fillna("unscored")

    coarse = pd.read_csv(COARSE, sep="\t")[["cs_key", "concentration", "dominant_mech", "locus_id"]]
    conc = dict(zip(coarse["cs_key"], coarse["concentration"]))
    cdom = dict(zip(coarse["cs_key"], coarse["dominant_mech"]))
    cloc = dict(zip(coarse["cs_key"], coarse["locus_id"]))

    rows = []
    for csk, g in mem.groupby("cs_key"):
        w = g["susie_pip"].to_numpy(float)
        s = w.sum()
        if s <= 0:
            continue
        w = w / s
        scored = (g["fine_mech"] != "unscored").to_numpy()
        coverage = float(w[scored].sum())
        mass = {c: float(np.sum(w[(g["fine_mech"] == c).to_numpy()])) for c in FINE}
        mass_unscored = float(np.sum(w[~scored]))
        # renormalise the fine posterior over the SCORED mass (report coverage separately)
        if coverage > 0:
            fine_mass = {c: mass[c] / coverage for c in FINE}
        else:
            fine_mass = {c: np.nan for c in FINE}
        dom = max(fine_mass, key=lambda k: (fine_mass[k] if fine_mass[k] == fine_mass[k] else -1))
        dom_mass = fine_mass[dom]
        rows.append(dict(
            cs_key=csk, locus_id=cloc.get(csk, ""), concentration=conc.get(csk, ""),
            coarse_dominant=cdom.get(csk, ""), n_members=len(g),
            scored_coverage=round(coverage, 4), mass_unscored=round(mass_unscored, 4),
            dominant_fine_mech=dom, dominant_fine_mass=round(dom_mass, 4) if dom_mass == dom_mass else None,
            fine_defensible=bool(dom_mass == dom_mass and dom_mass >= 0.5 and coverage >= 0.5),
            **{f"fmass_{c.split('-')[0].split('/')[0]}": round(fine_mass[c], 4)
               if fine_mass[c] == fine_mass[c] else None for c in FINE}))
    res = pd.DataFrame(rows)
    res.to_csv(OUT_TSV, sep="\t", index=False)

    # summary: restrict to CS with adequate scored coverage
    ok = res[res["scored_coverage"] >= 0.5]
    diffuse = ok[ok["concentration"] == "diffuse"]
    n = len(ok)
    summary = dict(
        script="91_pip_weighted_fine_mechanism.py",
        model="AlphaGenome zero-shot fingerprint (src/90) + src/72 lead scores; ANNOTATION ONLY, never trained on",
        taxonomy_caveat="fine regulatory sub-mechanism is a hand-built UNVALIDATED classifier (critique M6) -> HYPOTHESIS-LEVEL",
        n_cs_total=len(res),
        n_cs_scored_coverage_ge_0p5=n,
        fine_defensible_dommass_ge_0p5=int(ok["fine_defensible"].sum()),
        pct_fine_defensible=round(100 * ok["fine_defensible"].sum() / n, 1) if n else None,
        diffuse_cs_scored=len(diffuse),
        diffuse_fine_defensible=int(diffuse["fine_defensible"].sum()),
        diffuse_fine_mech_mix=diffuse["dominant_fine_mech"].value_counts().to_dict(),
        overall_dominant_fine_mix=ok["dominant_fine_mech"].value_counts().to_dict(),
        honest_framing=(f"Of {n} credible sets with >=50% scored posterior coverage, "
                        f"{int(ok['fine_defensible'].sum())} have a PIP-weighted dominant fine "
                        f"sub-mechanism carrying >=50% mass (defensible despite variant uncertainty); "
                        f"the rest have a fine sub-mechanism that varies across the credible set "
                        f"(hypothesis only). Among the {len(diffuse)} diffuse credible sets, "
                        f"{int(diffuse['fine_defensible'].sum())} still resolve to a consistent fine "
                        f"sub-mechanism. Fine taxonomy is unvalidated (M6) -> hypothesis-level; the "
                        f"load-bearing result stays the Tier-1 coarse posterior (src/89)."),
        firewall="Annotation only; never summed into convergence; never a per-locus direction call.",
    )
    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
    print(f"\nwrote {OUT_TSV} ({len(res)} CS)")


if __name__ == "__main__":
    main()

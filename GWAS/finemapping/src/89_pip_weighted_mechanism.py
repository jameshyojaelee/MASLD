#!/usr/bin/env python3
"""
89_pip_weighted_mechanism.py — Phase-6 WS-3 FULL: within-credible-set PIP-weighted
mechanism posterior for the seqfunc nominations (red-team critique G1, the real fix).

WHY THIS SUPERSEDES src/86 (the interim reliability flag):
  86 flagged whether a nomination's single lead is PIP-concentrated. It did NOT do the
  thing G1 actually asks: "compute a PIP-weighted expectation over ALL variants in the
  credible set and carry the posterior mass that supports the relevant mechanism class ...
  Preserve GWAS, ancestry, trait, signal, fine-mapping method, credible-set id. Do NOT
  merge signals merely because they fall within 500 kb."

UNIT (correct): a credible set = (study, trait, locus, susie_cs). Within it, susie_pip
  sums to ~1 (verified 1.000). max_pip in the substrate is a MAX across all 55 GWAS x CS
  and is NOT a within-CS posterior — so 500-kb clumps were the WRONG unit for weighting.
  The 81 nominated leads sit in 253 credible sets (a lead recurs across GWAS/trait signals).

WHAT THIS COMPUTES, per credible set:
  * PIP-weighted coarse mechanism posterior mass P(coding), P(splice), P(regulatory),
    P(unannotated) = sum_i w_i [class_i == c], w_i = susie_pip_i (already sums to 1).
  * effective #variants = 1 / sum(w_i^2)  (concentration; ~1 = single-variant-resolved).
  * PIP-weighted Borzoi magnitude E[|logSED|] over members that carry a Borzoi score.
  * dominant mechanism + its posterior mass; flag dominant_mass < 0.5 = multi-mechanism CS.
  * the single lead-variant mechanism call, and whether PIP-weighting CHANGES it
    (dominant_mech != lead_mech) -- the loci where the old one-variant approach misled.

Coarse mechanism = substrate var_class (coding/splice/regulatory) primary, VEP consequence
  fallback, else 'unannotated' (NEVER silently assumed regulatory).

APPLY-ONLY: annotation columns only; never summed into convergence; never a per-locus
  direction claim. AlphaGenome is NOT used here (fine sub-mechanism rescore = src/90).
"""
import json
import os

import numpy as np
import pandas as pd

ROOT = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FIN = os.path.join(ROOT, "GWAS/finemapping")
SF = os.path.join(FIN, "results/seqfunc")
FM = os.path.join(FIN, "results/combined_finemapping.csv")
MEMBERS = os.path.join(SF, "nominated_cs_members.tsv")   # cached extraction
LEADS = os.path.join(SF, "borzoi_magnitude_leads.tsv")
SUB = os.path.join(SF, "variant_substrate_hg38.tsv")
VEP = os.path.join(FIN, "results/credible_set_variant_consequences.csv")
CAND = os.path.join(SF, "borzoi_magnitude_candidates.tsv")
OUT_TSV = os.path.join(SF, "pip_weighted_mechanism_per_cs.tsv")
OUT_JSON = os.path.join(SF, "pip_weighted_mechanism_summary.json")


def extract_cs_members():
    """Pull the (study,trait,locus,susie_cs) credible sets that contain a nominated lead,
    with all members + within-CS susie_pip. Cached to MEMBERS."""
    if os.path.exists(MEMBERS):
        return pd.read_csv(MEMBERS, sep="\t")
    leads = pd.read_csv(LEADS, sep="\t")
    lead_ids = set(leads["variant_id_hg19"].astype(str))
    usecols = ["chromosome", "position", "allele1", "allele2", "trait", "locus",
               "study", "ancestry", "susie_pip", "susie_cs", "variant_id", "recommended_pip"]
    keep = []
    for chunk in pd.read_csv(FM, usecols=usecols, chunksize=2_000_000, low_memory=False):
        cs = pd.to_numeric(chunk["susie_cs"], errors="coerce").fillna(-1)
        keep.append(chunk[cs > 0])
    fm = pd.concat(keep, ignore_index=True)
    fm["cs_key"] = (fm["study"].astype(str) + "|" + fm["trait"].astype(str) + "|" +
                    fm["locus"].astype(str) + "|cs" + fm["susie_cs"].astype(str))
    lead_cs = fm[fm["variant_id"].astype(str).isin(lead_ids)]["cs_key"].unique()
    mem = fm[fm["cs_key"].isin(lead_cs)].copy()
    mem.to_csv(MEMBERS, sep="\t", index=False)
    return mem


def coarse_mechanism(mem):
    """coding / splice / regulatory / unannotated per member. substrate var_class primary,
    VEP consequence fallback."""
    sub = pd.read_csv(SUB, sep="\t")[["variant_id_hg19", "var_class", "consequence"]]
    sub = sub.drop_duplicates("variant_id_hg19")
    m = mem.merge(sub, left_on="variant_id", right_on="variant_id_hg19", how="left")
    vep = pd.read_csv(VEP)[["variant_id", "class", "Consequence"]].drop_duplicates("variant_id")
    m = m.merge(vep, on="variant_id", how="left")

    def classify(r):
        vc = str(r.get("var_class", "")).lower()
        if vc in ("coding", "splice", "regulatory"):
            return vc
        cons = str(r.get("consequence", "")) + " " + str(r.get("Consequence", ""))
        cl = str(r.get("class", "")).lower()
        if "splice" in cons.lower():
            return "splice"
        if cl == "coding_protein_altering" or "missense" in cons.lower() or "stop" in cons.lower():
            return "coding"
        if cl in ("noncoding", "coding_synonymous") or "regulatory" in cons.lower():
            return "regulatory"
        return "unannotated"

    m["mech"] = m.apply(classify, axis=1)
    return m


def main():
    mem = extract_cs_members()
    mem["susie_pip"] = pd.to_numeric(mem["susie_pip"], errors="coerce").fillna(0.0)
    m = coarse_mechanism(mem)

    # Borzoi magnitude per member (best-gene row)
    cand = pd.read_csv(CAND, sep="\t")
    cand = cand.sort_values("borzoi_abs_logsed", ascending=False).drop_duplicates("variant_id_hg19")
    m = m.merge(cand[["variant_id_hg19", "borzoi_abs_logsed"]].rename(
        columns={"variant_id_hg19": "vid_b"}), left_on="variant_id", right_on="vid_b", how="left")

    # nominated lead -> locus_id
    leads = pd.read_csv(LEADS, sep="\t")[["variant_id_hg19", "locus_id"]]
    lead_locus = dict(zip(leads["variant_id_hg19"].astype(str), leads["locus_id"].astype(str)))

    CLASSES = ["coding", "splice", "regulatory", "unannotated"]
    rows = []
    for csk, g in m.groupby("cs_key"):
        w = g["susie_pip"].to_numpy(float)
        s = w.sum()
        if s <= 0:
            continue
        w = w / s
        eff_n = 1.0 / np.sum(w ** 2)
        mass = {c: float(np.sum(w[(g["mech"] == c).to_numpy()])) for c in CLASSES}
        dom = max(mass, key=mass.get)
        dom_mass = mass[dom]
        # lead member of this CS
        li = int(np.argmax(g["susie_pip"].to_numpy()))
        lead_var = g.iloc[li]["variant_id"]
        lead_mech = g.iloc[li]["mech"]
        lead_pip = float(g.iloc[li]["susie_pip"])
        # is this CS's lead one of the nominated leads?
        nom_leads_here = [v for v in g["variant_id"].astype(str) if v in lead_locus]
        locus_id = lead_locus.get(str(lead_var)) or (lead_locus.get(nom_leads_here[0]) if nom_leads_here else "")
        # pip-weighted magnitude
        mag = pd.to_numeric(g["borzoi_abs_logsed"], errors="coerce").to_numpy()
        ok = ~np.isnan(mag)
        pw_mag = float(np.sum(w[ok] * mag[ok]) / w[ok].sum()) if ok.any() and w[ok].sum() > 0 else np.nan
        lead_mag = float(g.iloc[li]["borzoi_abs_logsed"]) if pd.notna(g.iloc[li]["borzoi_abs_logsed"]) else np.nan
        st, tr, lo, _ = csk.split("|")
        rows.append(dict(
            cs_key=csk, study=st, trait=tr, gwas_locus=lo, locus_id=locus_id,
            n_members=len(g), eff_n_variants=round(eff_n, 3),
            concentration=("concentrated" if eff_n <= 1.5 else ("moderate" if eff_n < 3 else "diffuse")),
            lead_variant=lead_var, lead_pip=round(lead_pip, 4), lead_mech=lead_mech,
            mass_coding=round(mass["coding"], 4), mass_splice=round(mass["splice"], 4),
            mass_regulatory=round(mass["regulatory"], 4), mass_unannotated=round(mass["unannotated"], 4),
            dominant_mech=dom, dominant_mass=round(dom_mass, 4),
            mechanism_defensible=bool(dom_mass >= 0.5),
            pip_weighting_changes_call=bool(dom != lead_mech),
            pip_weighted_borzoi_mag=(round(pw_mag, 4) if pw_mag == pw_mag else None),
            lead_borzoi_mag=(round(lead_mag, 4) if lead_mag == lead_mag else None),
        ))
    res = pd.DataFrame(rows)
    res.to_csv(OUT_TSV, sep="\t", index=False)

    n = len(res)
    conc = (res["concentration"] == "concentrated").sum()
    diff = (res["concentration"] == "diffuse").sum()
    defensible = int(res["mechanism_defensible"].sum())
    changed = int(res["pip_weighting_changes_call"].sum())
    # coarse architecture: PIP-weighted vs lead-only mechanism mix (per CS)
    lead_mix = res["lead_mech"].value_counts().to_dict()
    dom_mix = res["dominant_mech"].value_counts().to_dict()
    summary = dict(
        script="89_pip_weighted_mechanism.py",
        supersedes="src/86 reliability flag (which did not do the PIP-weighted expectation)",
        unit="credible set = (study, trait, locus, susie_cs); susie_pip sums to 1 within-CS",
        n_credible_sets=n, n_unique_members=int(m["variant_id"].nunique()),
        n_nominated_leads=81,
        concentration={"concentrated_eff<=1.5": int(conc),
                       "moderate": int(n - conc - diff), "diffuse_eff>=3": int(diff),
                       "median_eff_n": round(float(res["eff_n_variants"].median()), 3)},
        mechanism_defensible_dommass_ge_0p5=defensible,
        pct_mechanism_defensible=round(100 * defensible / n, 1),
        n_pip_weighting_changes_dominant_call=changed,
        pct_call_changes=round(100 * changed / n, 1),
        lead_only_mechanism_mix=lead_mix,
        pip_weighted_dominant_mechanism_mix=dom_mix,
        honest_framing=(f"Of {n} nominated credible sets, {conc} are PIP-concentrated "
                        f"(effective #variants <=1.5) where the single-variant mechanism IS the "
                        f"credible-set mechanism; {diff} are diffuse (eff>=3) where it is not. "
                        f"{defensible}/{n} ({round(100*defensible/n,1)}%) have a PIP-weighted "
                        f"dominant mechanism carrying >=50% posterior mass. PIP-weighting changes "
                        f"the dominant mechanism call vs the single lead in {changed} CS. Mechanism "
                        f"is a defensible LOCUS statement only for the concentrated/high-mass CS; "
                        f"elsewhere it is a hypothesis. No global coding-vs-regulatory architecture "
                        f"claim is made (consistent with the manuscript 15.9% coverage gate)."),
        firewall="Annotation only; never summed into convergence; never a per-locus direction call.",
    )
    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
    print(f"\nwrote {OUT_TSV} ({n} credible sets)")


if __name__ == "__main__":
    main()

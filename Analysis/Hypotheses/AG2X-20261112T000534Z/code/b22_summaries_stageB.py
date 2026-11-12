#!/usr/bin/env python3
"""Stage B descriptive summaries for RESULT_AG2X_stageB.md (requested by the orchestrator after the analysis job).

No test, threshold or variant set is changed here. Outputs: per-locus table for tier 1 (and tiers 2-3),
member versus background pass shares, SORT1 values in every readout, and a caQTL orientation diagnostic
(lead-file EA orientation and proxy-file beta against the Stage A A1 hg38 ALT beta). No hosted call.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

B = L.PKG / "stageB"
IN, TAB = B / "inputs", B / "tables"
REPO = L.REPO
CURRIN = REPO / "GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1"
A1 = REPO / "GWAS/finemapping/results/alphagenome_program/c1-endpoint2-v2-20260915T084559Z/tables/endpoint1_matched_tierA4.tsv.gz"
ENZ = REPO / "Analysis/Hypotheses/ENZFAT1-20261110T204401Z/out_S2/locus_effects.tsv"
BENCH = {"BL_ATAC": "bulk_liver|ATAC", "BL_DNASE": "bulk_liver|DNASE", "HEP_SCATAC": "hepatocyte|SC_ATAC"}
LABEL = {"DIG1_15": "MTARC1 (MARC1)", "DIG1_34": "TNFSF10", "DIG1_97": "PANX1", "DIG1_68": "GPT", "DIG1_86": "GPAM",
         "named_HSD17B13": "HSD17B13", "named_PNPLA3": "PNPLA3", "DIG1_148": "TM6SF2", "DIG1_6": "SORT1 (control)"}
out = {"candidate_not_adopted": True}


def flip(u):
    c, p, r, a = u.split(":")
    return f"{c}:{p}:{a}:{r}"


pv = pd.read_csv(TAB / "stageB_variant_fields.tsv.gz", sep="\t", low_memory=False)
imap = dict(pd.read_csv(IN / "indel_uid_map.tsv", sep="\t").values)
mv = pd.read_csv(IN / "member_variants.tsv", sep="\t")
mv["uid"] = mv.uid.map(lambda u: imap.get(u, u))
maxpip = dict(zip(mv.uid, mv.max_pip))
pctcols = [c for c in pv.columns if c.startswith("pct:")]
is_mem = pv.roles.str.contains("(?:^|;)member(?:$|;)", regex=True)
is_bg = pv.roles.str.contains("background")

# ---------------- pass shares, members vs background
rows = []
for c in pctcols:
    r = c[4:]
    for grp, m in (("background", is_bg), ("member", is_mem)):
        x = pv.loc[m, c].dropna()
        rows.append({"readout": r, "benchmarked": r in BENCH.values(), "group": grp, "n": int(len(x)),
                     "n_pass": int((x >= 0.99).sum()), "share_pass": float((x >= 0.99).mean())})
any3 = lambda d: (d[["pct:" + r for r in BENCH.values()]] >= 0.99).any(axis=1)  # noqa: E731
for grp, m in (("background", is_bg), ("member", is_mem)):
    d = pv[m].dropna(subset=["pct:" + r for r in BENCH.values()])
    rows.append({"readout": "any of BL_ATAC, BL_DNASE, HEP_SCATAC", "benchmarked": True, "group": grp, "n": int(len(d)),
                 "n_pass": int(any3(d).sum()), "share_pass": float(any3(d).mean())})
sh = pd.DataFrame(rows)
bgsh = sh[sh.group == "background"].set_index("readout").share_pass
sh["expected_n_pass_at_background_share"] = sh.n * sh.readout.map(bgsh)
sh["candidate_not_adopted"] = True
sh.to_csv(TAB / "stageB_pass_shares.tsv", sep="\t", index=False)
# by tier, combined rule, members only
bt = []
for t, g in pv[is_mem].groupby("tier"):
    d = g.dropna(subset=["pct:" + r for r in BENCH.values()])
    bt.append({"tier": int(t), "members": int(len(d)), "sequence_supported": int(any3(d).sum()),
               "expected_at_background_share": float(len(d) * bgsh["any of BL_ATAC, BL_DNASE, HEP_SCATAC"])})
out["members_by_tier_combined_rule"] = bt

# ---------------- per-locus table
loci = pd.read_csv(IN / "loci.tsv", sep="\t")
enz = pd.read_csv(ENZ, sep="\t")
enz = enz[~enz.locus.astype(str).str.startswith("ctrl_")]
enz["locus_id"] = "DIG1_" + enz.locus.astype(int).astype(str)
ecols = ["postmean_pdff", "lfsr_pdff", "postmean_cirr_ghouse", "lfsr_cirr_ghouse", "postmean_nafld_finngen",
         "lfsr_nafld_finngen"]
loci = loci.drop(columns=[c for c in ecols if c in loci.columns]).merge(enz[["locus_id"] + ecols], on="locus_id", how="left")
pvi = pv.set_index("uid")


def describe(u):
    s = pvi.loc[u]
    bench = [f"{k} {s['pct:' + r]:.4f}" for k, r in BENCH.items() if s["pct:" + r] >= 0.99]
    other = sorted({c[4:].split("|")[0] for c in pctcols if c[4:] not in BENCH.values() and s[c] >= 0.99})
    sign = s.accessibility_sign_alt_minus_ref
    return (f"{u} (max PIP {maxpip.get(u, float('nan')):.3f}; benchmarked: {', '.join(bench)}; "
            f"not benchmarked: {', '.join(other) if other else 'none'}; accessibility {sign}; "
            f"coding {'yes' if bool(s.coding_pip05) else 'no'})")


res = []
for _, lc in loci.sort_values("locus_id").iterrows():
    lid = lc.locus_id
    here = pv[is_mem & pv.loci.fillna("").str.split(";").map(lambda x: lid in x)]
    sets = int(lc.n_sets) if pd.notna(lc.n_sets) else 0
    eff = "; ".join(f"{t} postmean {lc['postmean_' + t]:.3g}, lfsr {lc['lfsr_' + t]:.3g}"
                    for t in ("pdff", "cirr_ghouse", "nafld_finngen") if pd.notna(lc.get("lfsr_" + t)))
    row = {"locus_id": lid, "tier": int(lc.tier), "label": LABEL.get(lid, ""), "genes": lc.genes if pd.notna(lc.genes) else "",
           "disease_effect": bool(lc.disease_effect), "enzfat1_effects": eff, "scored_as": lc.scored_as,
           "n_credible_sets": sets, "n_members_scored": int(len(here))}
    if lc.scored_as != "credible_set_members":
        row["sequence_supported_variant"] = "no member scored (lead and European r2 > 0.8 proxies only)"
        prox = pv[pv.loci.fillna("").str.split(";").map(lambda x: lid in x) & ~is_bg]
        row["lead_proxy_variants_scored"] = int(len(prox))
        row["lead_proxy_any_benchmarked_pass"] = int(any3(prox.dropna(subset=["pct:" + r for r in BENCH.values()])).sum())
    else:
        ss = here[here.sequence_supported]
        row["sequence_supported_variant"] = ("none (not evidence against the locus)" if len(ss) == 0 else
                                             "one" if len(ss) == 1 else f"several ({len(ss)})")
        row["expected_supported_at_background_share"] = round(len(here) * bgsh["any of BL_ATAC, BL_DNASE, HEP_SCATAC"], 2)
        row["supported_variants"] = " | ".join(describe(u) for u in ss.uid)
        row["coding_members_pip05"] = ";".join(here.uid[here.coding_pip05.astype(bool)])
    row["candidate_not_adopted"] = True
    res.append(row)
lt = pd.DataFrame(res)
lt[(lt.tier == 1) | (lt.locus_id == "DIG1_6")].to_csv(TAB / "stageB_tier1_locus_table.tsv", sep="\t", index=False)
lt[(lt.tier > 1) & (lt.locus_id != "DIG1_6")].to_csv(TAB / "stageB_tier2_3_locus_table.tsv", sep="\t", index=False)

# ---------------- SORT1 in every readout
s1 = pvi.loc["chr1:109274968:G:T"]
out["SORT1_all_readouts"] = {c[4:]: {"raw": float(s1["raw:" + c[4:]]), "pct": float(s1[c])} for c in pctcols}

# ---------------- caQTL orientation diagnostic against Stage A A1 (hg38 ALT beta)
a1 = pd.read_csv(A1, sep="\t", usecols=["chr_x", "pos_hg38", "ref", "alt", "beta_alt"], dtype={"pos_hg38": str})
a1["chr_x"] = a1.chr_x.astype(str)
a1["uid"] = a1.chr_x.where(a1.chr_x.str.startswith("chr"), "chr" + a1.chr_x) + ":" + a1.pos_hg38 + ":" + a1.ref + ":" + a1.alt
a1b = dict(zip(a1.uid, a1.beta_alt))
lead = pd.read_csv(CURRIN / "liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz",
                   sep="\t", usecols=["lead_variant_ID", "peak_ID", "beta", "EA", "NEA", "allele_flip"])
lc_ = lead.lead_variant_ID.str.split(":", expand=True)
lead["beta_alt"] = np.where(lead.EA == lc_[3], lead.beta, np.where(lead.EA == lc_[2], -lead.beta, np.nan))
prox = pd.read_csv(CURRIN / "caQTL_variants_overlappingPeaks_LD-r2-0.8_withLead.bed.gz", sep="\t",
                   usecols=["proxy_ID", "lead_ID", "caPeak", "proxy_caQTL_beta"])


def agree(ids, vals):
    s = []
    for u, b in zip(ids, vals):
        if u in a1b:
            s.append(np.sign(b) == np.sign(a1b[u]))
        elif flip(u) in a1b:
            s.append(np.sign(b) == -np.sign(a1b[flip(u)]))
    return {"n": len(s), "sign_agreement_with_A1_beta_alt": float(np.mean(s)) if s else None}


diag = {"lead_file_EA_oriented": agree(lead.lead_variant_ID, lead.beta_alt),
        "lead_file_raw_beta": agree(lead.lead_variant_ID, lead.beta),
        "proxy_file_lead_rows": agree(prox.loc[prox.proxy_ID == prox.lead_ID, "proxy_ID"],
                                      prox.loc[prox.proxy_ID == prox.lead_ID, "proxy_caQTL_beta"]),
        "lead_EA_equals_ID_alt_share": float((lead.EA == lc_[3]).mean()),
        "allele_flip_counts": lead.allele_flip.value_counts().to_dict()}
out["caqtl_orientation_diagnostic"] = diag
json.dump(out, open(TAB / "stageB_summaries.json", "w"), indent=1, default=str)
L.log(json.dumps({k: v for k, v in out.items() if k != "SORT1_all_readouts"}, default=str))

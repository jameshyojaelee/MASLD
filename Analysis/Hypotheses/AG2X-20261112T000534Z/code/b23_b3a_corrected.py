#!/usr/bin/env python3
"""POST HOC, written after B3a was seen. Corrects the B3a caQTL orientation only.

PRESPEC_AG2X_stageB.md section 5 oriented the Currin lead beta to the hg38 ALT by negating it when EA is the
hg38 REF. The diagnostic in stageB_summaries.json shows that rule is wrong: the raw lead `beta` agrees in sign
with the Stage A A1 hg38-ALT beta (PRESPEC_AG2X.md section 3: "beta_alt is already the hg38 ALT effect") at
4,177 of 4,177 shared leads, and the EA-oriented beta at 48%. Here the lead beta is used as the ALT effect,
the prespecified 99% proxy rule is re-applied with that orientation, and B3a is rebuilt with the same member
sets, supported rule (BL_ATAC) and outcome definition. B3b and B3c are unchanged. Original outputs are kept.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

TAB = L.PKG / "stageB/tables"
CURRIN = L.REPO / "GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1"


def flip(u):
    c, p, r, a = u.split(":")
    return f"{c}:{p}:{a}:{r}"


lead = pd.read_csv(CURRIN / "liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz",
                   sep="\t", usecols=["lead_variant_ID", "peak_ID", "beta"])
prox = pd.read_csv(CURRIN / "caQTL_variants_overlappingPeaks_LD-r2-0.8_withLead.bed.gz", sep="\t",
                   usecols=["proxy_ID", "lead_ID", "caPeak", "proxy_caQTL_beta"])
chk = prox[prox.proxy_ID == prox.lead_ID].merge(lead, left_on=["lead_ID", "caPeak"], right_on=["lead_variant_ID", "peak_ID"])
agree = float((np.sign(chk.proxy_caQTL_beta) == np.sign(chk.beta)).mean())
use_proxies = bool(agree >= 0.99)
eff = lead.rename(columns={"lead_variant_ID": "id", "beta": "beta_alt"})[["id", "beta_alt"]]
if use_proxies:
    eff = pd.concat([eff, prox.rename(columns={"proxy_ID": "id", "proxy_caQTL_beta": "beta_alt"})[["id", "beta_alt"]]])
eff = eff.assign(a=eff.beta_alt.abs()).sort_values("a", ascending=False).drop_duplicates("id")
m = dict(zip(eff.id, eff.beta_alt))

pv = pd.read_csv(TAB / "stageB_variant_fields.tsv.gz", sep="\t", usecols=["uid", "raw:bulk_liver|ATAC"]).set_index("uid")
b3 = pd.read_csv(TAB / "stageB_B3_input.tsv", sep="\t")
d = b3[b3.arm == "B3a_caqtl_matching_sign"].copy()  # same members, same supported flags, same set filter


def outcome(u):
    b = m.get(u, -m[flip(u)] if flip(u) in m else None)
    if b is None:
        return 0.0
    return float(np.sign(b) == np.sign(pv.at[u, "raw:bulk_liver|ATAC"]))


d["outcome"] = d.uid.map(outcome)
d["arm"] = "B3a_caqtl_matching_sign_corrected_posthoc"
d.to_csv(TAB / "stageB_B3a_corrected_input.tsv", sep="\t", index=False)
note = {"posthoc": True, "proxy_check_lead_rows": int(len(chk)), "proxy_sign_agreement_with_lead_beta": agree,
        "proxies_used": use_proxies, "members_in_arm": int(len(d)),
        "members_with_caqtl_entry": int(sum((u in m) or (flip(u) in m) for u in d.uid)),
        "members_with_matching_sign": int(d.outcome.sum()), "candidate_not_adopted": True}
json.dump(note, open(TAB / "stageB_B3a_corrected_notes.json", "w"), indent=1)
L.log(json.dumps(note))

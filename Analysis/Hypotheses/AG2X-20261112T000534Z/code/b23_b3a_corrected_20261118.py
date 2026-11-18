#!/usr/bin/env python3
"""Post-review 2026-11-18 Stage B3a Currin-corrected input on out_20261118/stageB (zero API calls).

Rebuilds the orientation-corrected Currin B3a input across the same credible sets and
supported_BL_ATAC rule, returning np.nan (NOT 0.0) when uid is unmeasured in Currin
so unmeasured variants are dropped instead of treated as discordant (N06).
"""
from __future__ import annotations

import gzip
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

IN = L.PKG / "stageB/inputs"
TAB = L.PKG / "out_20261118/stageB"
CURRIN = L.REPO / "GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1"


def flip(u: str) -> str:
    c, p, r, a = u.split(":")
    return f"{c}:{p}:{a}:{r}"


def main():
    lead = pd.read_csv(
        CURRIN / "liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz",
        sep="\t",
        usecols=["lead_variant_ID", "peak_ID", "beta"],
    )
    prox = pd.read_csv(
        CURRIN / "caQTL_variants_overlappingPeaks_LD-r2-0.8_withLead.bed.gz",
        sep="\t",
        usecols=["proxy_ID", "lead_ID", "caPeak", "proxy_caQTL_beta"],
    )
    chk = prox[prox.proxy_ID == prox.lead_ID].merge(
        lead, left_on=["lead_ID", "caPeak"], right_on=["lead_variant_ID", "peak_ID"]
    )
    agree = float((np.sign(chk.proxy_caQTL_beta) == np.sign(chk.beta)).mean())
    use_proxies = bool(agree >= 0.99)
    eff = lead.rename(columns={"lead_variant_ID": "id", "beta": "beta_alt"})[["id", "beta_alt"]]
    if use_proxies:
        eff = pd.concat([eff, prox.rename(columns={"proxy_ID": "id", "proxy_caQTL_beta": "beta_alt"})[["id", "beta_alt"]]])
    eff = eff.assign(a=eff.beta_alt.abs()).sort_values("a", ascending=False).drop_duplicates("id")
    m = dict(zip(eff.id, eff.beta_alt))

    pv = pd.read_csv(
        TAB / "stageB_variant_fields.tsv.gz",
        sep="\t",
        usecols=["uid", "raw:bulk_liver|ATAC", "supported_BL_ATAC"],
    ).set_index("uid")

    imap_path = IN / "indel_uid_map.tsv"
    imap = dict(pd.read_csv(imap_path, sep="\t").values) if imap_path.exists() else {}
    with gzip.open(IN / "credible_set_member_rows.tsv.gz", "rt") as fh:
        cs = pd.read_csv(fh, sep="\t", usecols=["locus_id", "set_id", "uid"])
    cs["uid"] = cs.uid.map(lambda u: imap.get(u, u) if isinstance(u, str) else u)
    cs = cs[cs.uid.isin(pv.index)].drop_duplicates()

    d = cs.assign(supported=cs.uid.map(pv["supported_BL_ATAC"])).dropna(subset=["supported"]).copy()
    d["supported"] = d["supported"].map(lambda x: 1 if str(x).strip().lower() in ("true", "1", "1.0") else 0)
    ok = d.groupby("set_id").supported.agg(lambda x: x.min() == 0 and x.max() == 1)
    d = d[d.set_id.isin(ok[ok].index)].drop_duplicates(["locus_id", "uid"]).copy()

    def outcome(u):
        b = m.get(u, -m[flip(u)] if flip(u) in m else None)
        if b is None:
            return np.nan
        raw_val = pv.at[u, "raw:bulk_liver|ATAC"] if u in pv.index else np.nan
        if not np.isfinite(raw_val):
            return np.nan
        return float(np.sign(b) == np.sign(raw_val))

    n_before_drop = int(len(d))
    d["outcome"] = d.uid.map(outcome)
    d = d.dropna(subset=["outcome"]).copy()
    d["arm"] = "B3a_caqtl_matching_sign_corrected_posthoc_measured_only"
    d = d[["arm", "locus_id", "uid", "supported", "outcome"]]
    d.to_csv(TAB / "stageB_B3a_corrected_input.tsv", sep="\t", index=False)

    note = {
        "posthoc": True,
        "proxy_check_lead_rows": int(len(chk)),
        "proxy_sign_agreement_with_lead_beta": agree,
        "proxies_used": use_proxies,
        "members_in_contrasting_sets_before_currin_filter": n_before_drop,
        "members_measured_in_currin": int(len(d)),
        "members_unmeasured_dropped_as_nan": n_before_drop - int(len(d)),
        "supported_measured": int((d.supported == 1).sum()),
        "supported_matching_sign": int(((d.supported == 1) & (d.outcome == 1)).sum()),
        "other_measured": int((d.supported == 0).sum()),
        "other_matching_sign": int(((d.supported == 0) & (d.outcome == 1)).sum()),
        "candidate_not_adopted": True,
    }
    (TAB / "stageB_B3a_corrected_notes.json").write_text(json.dumps(note, indent=1))
    L.log(json.dumps(note))


if __name__ == "__main__":
    main()

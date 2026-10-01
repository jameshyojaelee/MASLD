#!/usr/bin/env python3
"""Model A-final tag table (spec v1 section 1.2): one row per T1 v2 tag row, exclusions as flags.

Inputs (all outcome-free):
  --t3prime     T3' v2 tag_exclusions_v2.tsv: steps 1-4 of spec 1.2 (non-autosomal, MHC/IG/TR span,
                alt_copy, pseudogene/IG/TR biotype, NUMT, imprinted, multi-gene exon, REDIportal, T3'
                fail or untestable); column kept_steps_1_4.
  --freq, --ld  tag_2pq_by_superpop.tsv and tag_ld_by_superpop.tsv for the full v2 list (step 5: a tag
                needs both rows).
  --ffpe-drop   comma list of cohorts whose transition tags are excluded (step 6: QC v2 null-site
                G1-FFPE verdict not pass).
Columns written: gene_id, gene_name, lead_variant_id, chrom, pos, ref, alt, orientation, s_t
(+1 when the lead-ALT allele is the tag REF, i.e. orientation -1), transition, kept_primary
(steps 1-5), drop_in_<cohort> for step 6, and sensitivity flags: multiallelic_1kgp, low_r2_larger_ref
(L16), het_excess_1kgp (obs/expected REF/ALT het > 1.2 with n_hap >= 400 in any superpopulation),
r2_ge_0.8_<POP> and sign_flip_<POP> (ancestry-matched validity), t3prime_alt_background (T3'
sensitivity flag). The table and its sha256 are written to --out.
"""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
TRANSITIONS = {frozenset("CT"), frozenset("AG")}
KEY = ["gene_id", "chrom", "pos", "ref", "alt"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t3prime", required=True)
    ap.add_argument("--freq", required=True)
    ap.add_argument("--ld", required=True)
    ap.add_argument("--ffpe-drop", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    t = pd.read_csv(a.t3prime, sep="\t", dtype={"pos": int})
    t["kept_steps_1_4"] = t["kept_steps_1_4"].astype(str) == "True"
    fr = pd.read_csv(a.freq, sep="\t", dtype={"pos": int})
    ld = pd.read_csv(a.ld, sep="\t", dtype={"pos": int})

    cols = ["gene_id", "gene_name", "lead_variant_id", "chrom", "pos", "ref", "alt", "orientation",
            "kept_steps_1_4", "sens_alt_background_any"]
    tab = t[cols].rename(columns={"sens_alt_background_any": "t3prime_alt_background"})
    tab["orientation"] = pd.to_numeric(tab["orientation"]).astype(int)
    tab["s_t"] = (tab["orientation"] == -1).map({True: 1, False: -1})
    tab["transition"] = [frozenset((r, x)) in TRANSITIONS for r, x in zip(tab["ref"], tab["alt"])]
    tab["t3prime_alt_background"] = tab["t3prime_alt_background"].astype(str) == "True"

    ratio = []
    for p in POPS:
        ok = fr[f"n_hap_{p}"] >= 400
        ratio.append(ok & (fr[f"obs_het_ref_alt_{p}"] > 1.2 * fr[f"het_ref_alt_{p}"]))
    fr["het_excess_1kgp"] = pd.concat(ratio, axis=1).any(axis=1)
    fr["multiallelic_1kgp"] = fr["multiallelic_in_1kgp"].astype(str) == "True"
    tab = tab.merge(fr[["chrom", "pos", "ref", "alt", "het_excess_1kgp", "multiallelic_1kgp"]].assign(has_freq=True),
                    on=["chrom", "pos", "ref", "alt"], how="left")
    ldc = ld[KEY + ["low_r2_eur_eas_larger_ref"] + [f"r2_{p}" for p in POPS] + [f"sign_flip_{p}" for p in POPS]]
    tab = tab.merge(ldc.assign(has_ld=True), on=KEY, how="left")
    for c in ["has_freq", "has_ld", "het_excess_1kgp", "multiallelic_1kgp"]:
        tab[c] = tab[c].astype("boolean").fillna(False).astype(bool)
    tab["low_r2_larger_ref"] = tab["low_r2_eur_eas_larger_ref"].astype(str) == "True"
    for p in POPS:
        tab[f"r2_ge_0.8_{p}"] = tab[f"r2_{p}"] >= 0.8
        tab[f"sign_flip_{p}"] = tab[f"sign_flip_{p}"].astype(str) == "True"
    tab = tab.drop(columns=["low_r2_eur_eas_larger_ref"] + [f"r2_{p}" for p in POPS])

    tab["kept_primary"] = tab["kept_steps_1_4"] & tab["has_freq"] & tab["has_ld"]
    ffpe = [c for c in a.ffpe_drop.split(",") if c]
    for c in ffpe:
        tab[f"drop_in_{c}"] = tab["transition"]

    tab = tab.sort_values(["chrom", "pos", "gene_id"]).reset_index(drop=True)
    p = out / "tag_table.tsv"
    tab.to_csv(p, sep="\t", index=False)
    sha = hashlib.sha256(p.read_bytes()).hexdigest()
    (out / "tag_table.tsv.sha256").write_text(f"{sha}  tag_table.tsv\n")
    k = tab[tab["kept_primary"]]
    counts = {"rows": len(tab), "kept_primary_rows": len(k), "kept_primary_snvs": int(k[["chrom", "pos"]].drop_duplicates().shape[0]),
              "kept_primary_genes": int(k["gene_id"].nunique()),
              "removed_step5_no_freq_or_ld": int((tab["kept_steps_1_4"] & ~tab["kept_primary"]).sum()),
              "kept_transition_rows": int(k["transition"].sum()),
              "kept_genes_without_transversion_tag": int(k.groupby("gene_id")["transition"].all().sum()),
              "ffpe_drop_cohorts": ffpe, "sha256": sha,
              "sensitivity_flags_among_kept": {c: int(k[c].sum()) for c in
                                               ["multiallelic_1kgp", "low_r2_larger_ref", "het_excess_1kgp", "t3prime_alt_background"]
                                               + [f"sign_flip_{q}" for q in POPS]},
              "kept_rows_r2_ge_0.8": {q: int(k[f"r2_ge_0.8_{q}"].sum()) for q in POPS}}
    (out / "tag_table_counts.json").write_text(json.dumps(counts, indent=1))
    print(json.dumps(counts, indent=1))


if __name__ == "__main__":
    main()

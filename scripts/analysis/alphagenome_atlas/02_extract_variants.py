#!/usr/bin/env python3
"""Step 02: stream per-variant posteriors for the eligible signals (memory-bounded)."""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict

import lib_atlas as la

P = la.prespec()
OUT = la.out_root() / "tables"
AUDIT = la.PROJECT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/genetics/context/variant_liftover_audit.tsv"
U35 = la.PROJECT / "GWAS/finemapping/runs/uniform35_v4_2026-08-04/aggregate/susie_variant_results.tsv.gz"

A_COLS = ["signal_uid", "hg19_variant_id", "allele1", "allele2", "snp_pp_h4", "mapping_status",
          "hg38_chrom", "hg38_position_1based", "hg38_ref", "hg38_alt", "orientation"]
B_COLS = ["signal_uid", "locus_id", "study_name", "chromosome", "position", "rsid", "effect_allele",
          "other_allele", "beta", "se", "pip", "cs_id", "palindromic", "strand_resolution",
          "strand_source", "primary_eligible"]


def main() -> None:
    signals = la.read_tsv(OUT / "eligible_signals_provisional.tsv")
    key_a = {(s["gwas_name"], s["ensembl"], s["signal_pair_index"]): s["signal_uid"]
             for s in signals if s["universe"] in ("A_direct", "C_enzyme")}
    key_b = {(s["locus_id"], s["cs_id"]): s["signal_uid"] for s in signals if s["universe"] == "B_direct"}

    # Universe A/C: single pass over the 22 M-row audit.
    mass = defaultdict(float)
    dup = set()
    n_keep = 0
    csv.field_size_limit(1 << 30)
    with open(AUDIT, newline="") as src, la.open_text(OUT / "variant_posteriors_A.tsv.gz", "wt") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(dst, fieldnames=A_COLS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for r in reader:
            uid = key_a.get((r["gwas_name"], r["ensembl"], r["signal_pair_index"]))
            if uid is None:
                continue
            pair = (uid, r["hg19_variant_id"], r["allele1"], r["allele2"])
            if pair in dup:
                raise la.ContractError(f"duplicated posterior row: {pair}")
            dup.add(pair)
            w = float(r["snp_pp_h4"])
            if not (0 <= w <= 1):
                raise la.ContractError(f"invalid probability {w} for {pair}")
            mass[uid] += w
            writer.writerow({c: r[c] if c != "signal_uid" else uid for c in A_COLS})
            n_keep += 1
    bad = {k: v for k, v in mass.items() if abs(v - 1) > 1e-6}
    if bad:
        raise la.ContractError(f"SNP.PP.H4 does not sum to one for {len(bad)} signals, e.g. {list(bad.items())[:3]}")
    missing = set(key_a.values()) - set(mass)
    if missing:
        raise la.ContractError(f"{len(missing)} A/C signals without posterior rows, e.g. {sorted(missing)[:3]}")
    la.log(f"A/C posterior rows kept: {n_keep} across {len(mass)} signals")

    # Universe B: credible-set members from the uniform35 aggregate.
    n_b = 0
    pip_mass = defaultdict(float)
    with gzip.open(U35, "rt") as src, open(OUT / "variant_posteriors_B_hg19.tsv", "w") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(dst, fieldnames=B_COLS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for r in reader:
            if r["tier"] != "1" or r["cs_id"] in ("0", "", "NA"):
                continue
            uid = key_b.get((r["locus_id"], r["cs_id"]))
            if uid is None:
                raise la.ContractError(f"credible-set member without signal: {r['locus_id']} cs {r['cs_id']}")
            pip_mass[uid] += float(r["pip"])
            writer.writerow({"signal_uid": uid, **{c: r[c] for c in B_COLS if c != "signal_uid"}})
            n_b += 1
    missing_b = set(key_b.values()) - set(pip_mass)
    if missing_b:
        raise la.ContractError(f"{len(missing_b)} B signals without members, e.g. {sorted(missing_b)[:3]}")
    la.log(f"B members: {n_b} across {len(pip_mass)} credible sets; PIP mass range "
           f"{min(pip_mass.values()):.3f}-{max(pip_mass.values()):.3f}")


if __name__ == "__main__":
    main()

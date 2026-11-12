#!/usr/bin/env python3
"""ENZFAT1 S0: allele-orientation checks from the raw rows written by 00_build_input.py.

Hard checks stop the job (exit 1). Soft checks are asserted only when |z| >= 2 in that file.
Usage: 01_orientation.py <out_dir>
"""
import sys

import pandas as pd

ALLELE = {"rs738409": "G", "rs58542926": "T", "rs1260326": "T", "rs641738": "T"}
FAT_ALT = ["pdff", "alt_ukb", "alt_mvp"]
HARD = [(s, t) for s in ("rs738409", "rs58542926") for t in FAT_ALT]
HARD += [("rs738409", t) for t in ("cirr_ghouse", "cirr_mvp", "nafld_mvp", "nafld_finngen")]
HARD += [("rs1260326", "pdff")]
SOFT = [("rs641738", t) for t in FAT_ALT] + [("rs1260326", t) for t in ("alt_ukb", "alt_mvp")]
TRAITS = ["pdff", "alt_ukb", "alt_mvp", "cirr_ghouse", "cirr_mvp", "nafld_mvp", "nafld_finngen", "ast_ukb",
          "ggt_ukb"]


def main(out_dir):
    raw = pd.concat([pd.read_csv(f"{out_dir}/checks/{t}.tsv", sep="\t") for t in TRAITS], ignore_index=True)
    rows, failed = [], False
    for snp, allele in ALLELE.items():
        for trait in TRAITS:
            r = raw[(raw["SNP"] == snp) & (raw["trait"] == trait)]
            kind = "hard" if (snp, trait) in HARD else "soft" if (snp, trait) in SOFT else "report"
            if len(r) != 1:
                status = "FAIL_missing" if kind == "hard" else "missing"
                rows.append(dict(snp=snp, allele=allele, trait=trait, kind=kind, z=None, beta=None, status=status))
                failed |= kind == "hard"
                continue
            r = r.iloc[0]
            if r["A1"] == allele:
                sign = 1.0
            elif r["A2"] == allele:
                sign = -1.0
            else:
                rows.append(dict(snp=snp, allele=allele, trait=trait, kind=kind, z=None, beta=None,
                                 status="FAIL_allele_absent" if kind != "report" else "allele_absent"))
                failed |= kind != "report"
                continue
            z, beta = sign * r["Z"], sign * r["BETA"]
            if kind == "hard":
                status = "pass" if z > 0 else "FAIL"
            elif kind == "soft":
                status = "not_asserted_absz_lt_2" if abs(z) < 2 else ("pass" if z > 0 else "FAIL")
            else:
                status = "reported"
            failed |= status == "FAIL"
            rows.append(dict(snp=snp, allele=allele, trait=trait, kind=kind, z=round(z, 3), beta=beta,
                             status=status))
    out = pd.DataFrame(rows)
    out.to_csv(f"{out_dir}/orientation_checks.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    if failed:
        sys.exit("orientation check failed; inputs are not munged")


if __name__ == "__main__":
    main(sys.argv[1])

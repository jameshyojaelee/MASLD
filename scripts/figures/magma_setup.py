#!/usr/bin/env python
"""
magma_setup.py
==============
One-time setup for the MAGMA-anchored scDRS pipeline.

Steps:
  1. Read the 1000G EUR .bim (~9.4M rsID-keyed SNPs).
  2. Write a SNP location file for MAGMA: cols = SNP, CHR, BP.
  3. Build a chr:pos → rsID lookup table for joining our GWAS sumstats
     (which use chr:pos, no rsID) to the 1000G reference.

Outputs:
  tools/magma/work/g1000_eur.snploc   (MAGMA --snp-loc input)
  tools/magma/work/chr_pos_to_rsid.tsv (chr_pos<TAB>rsid; for joining sumstats)
"""
from __future__ import annotations

import os
import time
from pathlib import Path

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
BIM_PATH       = BASE / "tools/magma/ref/g1000_eur.bim"
SNPLOC_PATH    = BASE / "tools/magma/work/g1000_eur.snploc"
LOOKUP_PATH    = BASE / "tools/magma/work/chr_pos_to_rsid.tsv"
SNPLOC_PATH.parent.mkdir(parents=True, exist_ok=True)


def log(msg: str) -> None:
    print(f"[magma_setup] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


def main() -> None:
    log(f"reading {BIM_PATH}")
    t0 = time.time()
    n = 0
    with open(BIM_PATH) as f_in, \
         open(SNPLOC_PATH, "w") as f_snp, \
         open(LOOKUP_PATH, "w") as f_lookup:
        f_snp.write("SNP\tCHR\tBP\n")
        f_lookup.write("chr_pos\trsid\n")
        for line in f_in:
            # .bim cols: chr, rsid, cm, bp, a1, a2
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, rsid, _, bp = parts[0], parts[1], parts[2], parts[3]
            f_snp.write(f"{rsid}\t{chrom}\t{bp}\n")
            f_lookup.write(f"{chrom}:{bp}\t{rsid}\n")
            n += 1
            if n % 1_000_000 == 0:
                log(f"  {n:,} SNPs processed...")
    log(f"wrote {SNPLOC_PATH}  ({n:,} SNPs;  {time.time()-t0:.1f}s)")
    log(f"wrote {LOOKUP_PATH}")
    log("done.")


if __name__ == "__main__":
    main()

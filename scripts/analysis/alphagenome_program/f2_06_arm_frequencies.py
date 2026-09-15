#!/usr/bin/env python3
"""F2 step 6: which of the four scored arms is a haplotype that EXISTS in Europeans.

The four-arm additivity decomposition (REF, V1, V2, V1+V2) treats the two single-variant arms as if they
were observable sequences. At r2 = 1 they are not: only two of the four two-locus haplotypes segregate.
This reads the PLINK 1.9 --ld haplotype tables written by f2_05_gnmt_ld.sbatch and labels each scored arm
with its measured frequency in the project's 1000 Genomes EUR panel (379 founders, 758 chromosomes).

Output: tables/gnmt_arm_haplotype_frequencies.tsv
"""

from __future__ import annotations

import os
import pathlib
import re

OUT = pathlib.Path(os.environ["F2_OUT_ROOT"])
W = OUT / "raw" / "gnmt_ld"

# hg19 position -> (rsid, REF on GRCh38, ALT on GRCh38)
V = {42920549: ("rs9471976", "G", "T"),
     42928017: ("rs11752813", "C", "G"),
     42928758: ("rs2296805", "T", "G"),
     42931261: ("rs2296804", "C", "G")}


def parse(path: pathlib.Path) -> tuple[str, str, float, float, dict[str, float]] | None:
    text = path.read_text()
    m = re.search(r"--ld (\S+) (\S+):", text)
    if not m:
        return None
    a, b = m.group(1), m.group(2)
    st = re.search(r"R-sq = (\S+)\s+D' = (\S+)", text)
    haps = {}
    for line in text.splitlines():
        mm = re.match(r"\s{5,}([ACGT]{2})\s+(-?[\d.eE+-]+)\s+[\d.eE+-]+\s*$", line)
        if mm:
            haps[mm.group(1)] = max(0.0, float(mm.group(2)))
    return a, b, float(st.group(1)), float(st.group(2)), haps


def main() -> None:
    rows = []
    seen = set()
    for path in sorted(W.glob("ld_*.log")):
        got = parse(path)
        if not got:
            continue
        a, b, r2, dprime, haps = got
        pa, pb = int(a.split(":")[1]), int(b.split(":")[1])
        if pa > pb or (pa, pb) in seen or pa not in V or pb not in V:
            continue
        seen.add((pa, pb))
        (rs1, ref1, alt1), (rs2, ref2, alt2) = V[pa], V[pb]
        for arm, al1, al2 in (("ref", ref1, ref2), ("v1", alt1, ref2), ("v2", ref1, alt2), ("joint", alt1, alt2)):
            f = haps.get(al1 + al2, haps.get(al2 + al1))
            rows.append({"pair": f"{rs1}|{rs2}", "rsid1": rs1, "rsid2": rs2, "r2_eur": r2, "dprime_eur": dprime,
                         "arm": arm, "allele1": al1, "allele2": al2,
                         "eur_haplotype_frequency": "" if f is None else f"{f:.6f}",
                         "exists_in_eur_panel": "" if f is None else ("yes" if f > 0 else "no")})
    cols = ["pair", "rsid1", "rsid2", "r2_eur", "dprime_eur", "arm", "allele1", "allele2",
            "eur_haplotype_frequency", "exists_in_eur_panel"]
    path = OUT / "tables" / "gnmt_arm_haplotype_frequencies.tsv"
    with open(path, "w") as h:
        h.write("\t".join(cols) + "\n")
        for r in rows:
            h.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"{len(rows)} arm rows over {len(seen)} pairs -> {path}")
    for r in rows:
        print(r["pair"], r["arm"], r["allele1"] + r["allele2"], r["eur_haplotype_frequency"], r["exists_in_eur_panel"])


if __name__ == "__main__":
    main()

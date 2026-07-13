#!/usr/bin/env python
"""build_all_guide_coords.py — genomic-coordinate table for EVERY targeting guide in
the final Cas13 library (data/guides/library.csv), excluding non-targeting controls.

Reuses the trusted sequence-search coordinate logic from spotcheck_100_guides.analyse()
(target_seq located in a real annotated transcript, then mapped transcript->genome), so
the coordinates match the manual-QC spot-check figures. Output feeds the plotgardener
whole-library coverage figure (spotcheck_guides_plotgardener.R).

Excludes control_class == "non_targeting" (500 guides with no genomic target); keeps
target / positive_control / essential guides.

Output:
  Cas13_Library_Design/data/guide_qc_all_library.tsv   (per-guide coordinates + region + scores)
"""
from __future__ import annotations
import csv, os, sys
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spotcheck_100_guides as sc   # reuse analyse() + coordinate helpers

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
LIB  = ROOT / "Cas13_Library_Design/data/guides/library.csv"
OUT  = Path(os.environ.get("OUT_TSV", str(ROOT / "Cas13_Library_Design/data/guide_qc_all_library.tsv")))


def main():
    rows = list(csv.DictReader(open(LIB)))
    tgt  = [r for r in rows if r["control_class"] != "non_targeting"]
    per  = defaultdict(list)
    for r in tgt:
        per[r["gene_symbol_mouse"]].append(r)
    genes = list(per.keys())          # library (file) order
    print(f"[all] {len(tgt)} targeting guides / {len(genes)} genes "
          f"({len(rows) - len(tgt)} non_targeting excluded)")

    res = sc.analyse(per, genes)

    with open(OUT, "w") as out:
        out.write("gene\tbiotype\ttier\tref_transcript\tguide\tregion_pool\tregion_recomputed\t"
                  "targets_gene\trc_matches_target\ton_reference_tx\tgenomic_pos\tmrna_pos\t"
                  "tiger\tcas13design\tcontrol_class\tgene_source\n")
        for g in genes:
            d = res[g]; cc = per[g][0]["control_class"]
            for x in d["guides"]:
                out.write(f"{g}\t{d['biotype']}\t{x['tier']}\t{d['ref']}\t{x['id']}\t{x['region_data']}\t"
                          f"{x['region_calc']}\t{x['targets']}\t{x['rc_ok']}\t{x['on_ref']}\t"
                          f"{x['gpos']}\t{x['ppos']}\t{x['tiger']}\t{x['cas13']}\t{cc}\tall_library\n")

    allg = [x for g in genes for x in res[g]["guides"]]
    n = len(allg)
    targ = sum(x["targets"] for x in allg)
    rcok = sum(x["rc_ok"] for x in allg)
    mapped = sum(1 for x in allg if x["gpos"] is not None)
    print(f"[all] guides={n}  target_found={targ}/{n}  rc_ok={rcok}/{n}  genomic_mapped={mapped}/{n}")
    print(f"[all] TSV -> {OUT}")


if __name__ == "__main__":
    main()

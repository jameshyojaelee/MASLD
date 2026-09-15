#!/usr/bin/env python3
"""F3 side check: does the existing P5 deposit's `rna_readout` column say which readout was actually used?

In 50_haplotype_additivity.py the readout CHOICE requires the gene to overlap the 1-Mb window:

    rna_regions = [(span[1], span[2])] if span and span[0] == chrom and span[2] > start0 and span[1] < end else [(mid-10_000, mid+10_000)]

but the LABEL written to the table only checks the chromosome:

    "rna_readout": "gene_span" if span and span[0] == chrom else "central_20kb"

For a null draw, which is a uniform position on the sampled pair's chromosome carrying that pair's gene, the
two conditions come apart. This recomputes the actual choice for every row of both deposited tables and
reports how many rows the column misdescribes.

Output: tables/p5_deposit_rna_readout_check.json
"""

from __future__ import annotations

import csv
import json
import os
import pathlib

import f2_recipe as R

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT = pathlib.Path(os.environ["F2_OUT_ROOT"])
P5 = PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/tables"


def check(path: pathlib.Path, spans: dict) -> dict:
    n = miss = 0
    actual = {"gene_span": 0, "central_20kb": 0}
    labelled = {"gene_span": 0, "central_20kb": 0}
    with open(path) as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            n += 1
            p1, p2 = int(r["pos1"]), int(r["pos2"])
            chrom = r["chrom"]
            mid = (p1 + p2) // 2
            start0 = max(0, mid - R.WINDOW // 2)
            end = start0 + R.WINDOW
            span = spans.get(r["ensembl"])
            used = "gene_span" if (span and span[0] == chrom and span[2] > start0 and span[1] < end) else "central_20kb"
            actual[used] += 1
            labelled[r["rna_readout"]] = labelled.get(r["rna_readout"], 0) + 1
            if used != r["rna_readout"]:
                miss += 1
    return {"rows": n, "rows_where_the_column_misdescribes_the_readout": miss,
            "readout_actually_used": actual, "readout_as_labelled_in_the_deposit": labelled}


def main() -> None:
    spans = R.gene_spans()
    out = {"observed_table": check(P5 / "haplotype_effects.tsv", spans),
           "null_table": check(P5 / "haplotype_null_effects.tsv", spans),
           "consequence": ("the RNA channel's observed arm and null arm sum predicted signal over different windows. "
                           "The deposit's own column hides that, because it reports the gene span whenever the gene is "
                           "on the same chromosome rather than when the gene is inside the scored window.")}
    json.dump(out, (OUT / "tables" / "p5_deposit_rna_readout_check.json").open("w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

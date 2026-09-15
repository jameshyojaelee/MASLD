#!/usr/bin/env python3
"""f2-haplotype-v2 amendment 01: match the null host gene's TYPE as well as its span width.

Written before any predict_sequence call in this package. Timestamps come from the machine clock at write
time and the digest is rewritten in the same second, so written_utc can never be later than the file mtime.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import pathlib

OUT = pathlib.Path(os.environ["F2V2_OUT_ROOT"])
HERE = pathlib.Path(__file__).resolve().parent

AMD = {
    "amendment": 1,
    "package": "f2-haplotype-v2",
    "written_utc": None,
    "state_at_writing": (
        "no predict_sequence call has been made in this package; raw/ does not exist and no residual, p-value, "
        "median or verdict produced by this package exists on disk. The draw list built under the original rule "
        "was inspected for its gene-type composition only."
    ),
    "changes": "defect_1_repair.draw_rule.host_gene, and the draw count per gene",
    "what_changed": {
        "was": ("host gene = any GENCODE v49 chr6 gene with span width in [2472, 4791] bp; 100 distinct host "
                "genes sampled without replacement."),
        "now": ("host gene = a chr6 gene with gene_type protein_coding and span width in [2472, 6388] bp "
                "(6388 = 2 x GNMT's 3,194 bp). 70 such genes host a qualifying pair, so the 100 draws are 70 "
                "genes contributing one pair each plus 30 of those genes, drawn without replacement, "
                "contributing a second pair whose two positions are disjoint from the first pair's."),
    },
    "why": (
        "the draw list under the original rule is 70 lncRNA, 6 pseudogene and 2 TEC of 100 hosts and only 22 "
        "protein_coding. GNMT is protein_coding and carries real liver RNA signal (its REF RNA sum over the "
        "3,194 bp span is 29,288 in v1's table). A 3-4 kb lncRNA or pseudogene span with almost no predicted RNA "
        "signal gives a log2 ratio of a near-zero sum, which is the same noise problem the repair exists to "
        "remove: it would replace an unmatched window WIDTH with an unmatched window SIGNAL LEVEL. Matching the "
        "gene type costs a factor of 1.33 in the width band (median span width rises from about 3.9 kb to about "
        "4-5 kb against GNMT's 3.2 kb) and buys a readout that carries transcribed signal."
    ),
    "measured_feasibility_before_writing_this": {
        "protein_coding_chr6_genes_width_2472_4791": {"genes": 73, "hosting_a_qualifying_pair": 39, "qualifying_pairs": 1247},
        "protein_coding_chr6_genes_width_2472_6388": {"genes": 107, "hosting_a_qualifying_pair": 70, "qualifying_pairs": 3251},
        "all_types_width_2472_4791": {"genes": 363, "hosting_a_qualifying_pair": 186, "qualifying_pairs": 6356},
        "why_not_the_tight_band": ("39 host genes cannot give 100 draws at one pair per gene, and 2-3 pairs inside "
                                   "one 3 kb span share their readout window and most of their sequence, so the "
                                   "draws would be strongly dependent."),
    },
    "why_100_draws_is_kept": (
        "the p floor is 1/(draws+1) = 0.0099 at 100 draws, identical to v1's gnmtnull set, so the repaired null "
        "and the superseded one have the same resolution and the verdicts are comparable on that axis."
    ),
    "dependence_this_introduces_and_how_it_is_reported": (
        "30 of the 100 draws share a host gene with another draw. The empirical p treats draws as exchangeable, so "
        "RESULTS.md reports the distinct-host-gene count (70) beside the draw count (100) and repeats the F2.1 "
        "verdict on the 70-draw one-pair-per-gene sub-null as a sensitivity check."
    ),
    "added_diagnostic_labelled_post_hoc": (
        "the REF RNA sum of each null draw over its own gene span, against GNMT's 29,288.250, as the direct measure "
        "of whether the repaired null matches the observed readout's signal level; and, if at least 20 draws fall "
        "within a factor of 2 of GNMT's REF RNA sum, the F2.1 RNA p recomputed inside that signal-matched subset. "
        "Both are labelled post hoc and neither replaces the primary verdict."
    ),
    "unchanged": (
        "the chromosome, the separation bin [2471.2, 3916.3], the requirement that both variants be common EUR SNVs "
        "inside the host span with REF taken from GRCh38, the exclusion of GNMT and of anything overlapping its span, "
        "the seed 20260915, the recipe, the +/-1 kb local readout, and every prediction P1 to P5."
    ),
    "predictions_unchanged": True,
}


def main() -> None:
    AMD["written_utc"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    src = HERE / "f2v2_00_amendment_01.json"
    src.write_text(json.dumps(AMD, indent=1) + "\n")
    dest = OUT / "prespec" / "f2v2_00_amendment_01.json"
    dest.write_text(src.read_text())
    lines = []
    for name in sorted(p.name for p in (OUT / "prespec").glob("*.json")):
        h = hashlib.sha256((OUT / "prespec" / name).read_bytes()).hexdigest()
        lines.append(f"{h}  {name}")
    (OUT / "prespec" / "f2v2_00_prespec.sha256").write_text("\n".join(lines) + "\n")
    stamp = {"written_utc": AMD["written_utc"],
             "file_mtime_utc": dt.datetime.fromtimestamp(dest.stat().st_mtime, dt.timezone.utc)
                                .strftime("%Y-%m-%dT%H:%M:%SZ")}
    print(json.dumps(stamp, indent=1))
    print((OUT / "prespec" / "f2v2_00_prespec.sha256").read_text())


if __name__ == "__main__":
    main()

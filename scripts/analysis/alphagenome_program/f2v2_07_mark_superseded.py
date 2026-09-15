#!/usr/bin/env python3
"""f2-haplotype-v2 step 7: write SUPERSEDED.md into the v1 directory. Adds one new file; changes nothing else."""

from __future__ import annotations

import datetime as dt
import json
import os
import pathlib

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
V1 = PROJECT / "GWAS/finemapping/results/alphagenome_program/f2-haplotype-20260914T230433Z"
OUT = pathlib.Path(os.environ["F2V2_OUT_ROOT"])


def main() -> None:
    S = json.loads((OUT / "tables" / "f2v2_summary.json").read_text())
    d1, d2 = S["defect_1_repair"], S["defect_2_repair"]
    rep, old = d1["REPAIRED_primary"], d1["V1_SUPERSEDED"]
    ro = d1["repaired_null_readout"]
    ch_names = {"rna": "RNA_SEQ", "atac": "ATAC", "dnase": "DNASE", "h3k27ac": "H3K27ac"}
    moved = d1["moved"]
    geo = S["readout_geometry"]
    dup = geo["duplicate_rows_in_the_199_pair_family"]
    gs = geo["gene_span_readout_bp"]
    chrom_ps = ", ".join("{} {:.4f} to {:.4f}".format(ch_names[c], moved[c]["v1_p"], moved[c]["v2_p"])
                         for c in ("atac", "dnase", "h3k27ac"))
    lines = [
        "# SUPERSEDED",
        "",
        f"This package is superseded by `{OUT.name}`, written {dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}.",
        "",
        "`RESULTS.md` and every table here are left exactly as published, so the record of what was reported stays",
        "readable. The numbers listed below are wrong here; take them from the superseding package.",
        "",
        "## Why",
        "",
        "**Defect 1, the primary GNMT verdict was not a matched comparison in the RNA channel.** All six",
        f"GNMT-region pairs read RNA out over the {ro['gnmt_rna_readout_bp']} bp GENCODE v49 GNMT span while all",
        f"{old['n_null_draws']} of their matched chr6 null draws read RNA out over a central",
        f"{ro['v1_null_rna_readout_bp']} bp window, a {ro['width_ratio_v1_null_over_observed']}-fold width mismatch,",
        "so the RNA arm of F2.1 was not a matched comparison. Rescored against a gene-span-matched chr6 null of",
        f"{rep['n_null_draws']} draws, the RNA verdict is "
        f"{'inside' if rep['F2.1_verdict_per_channel']['rna'] else 'OUTSIDE'} the central 95% with p "
        f"{rep['F2.1_p_per_channel']['rna']:.4f}, against "
        f"{'inside' if old['F2.1_verdict_per_channel']['rna'] else 'OUTSIDE'} with p "
        f"{old['F2.1_p_per_channel']['rna']:.4f} here: **the verdict and the p do not change.** The defect was in",
        "the comparison, not in the answer it produced, and F2.1 stands. What does change is the three chromatin",
        "channels' p-values, because the repaired draws sit in protein-coding gene spans rather than at arbitrary",
        f"chr6 positions: {chrom_ps}. Those verdicts also stand.",
        "",
        "**Defect 2, the both-variants-active stratum p-values in `tables/f3_summary_n2500.json` were computed",
        "against the full pooled 2,500-draw null instead of the stratum's own matched draws.** Recomputed with the",
        "null inside the stratum, BH q<0.10 rejections are:",
        "",
        "| channel | as published here | with the null inside the stratum |",
        "|---|---|---|",
    ]
    for ch in ("rna", "atac", "dnase", "h3k27ac"):
        a = d2["n2500"]["per_channel"][ch]
        v1b = a["V1_AS_PUBLISHED_p_against_the_pooled_null"]["n_bh_rejected"]
        v2b = a["REPAIRED_p_against_the_null_inside_the_stratum"]
        shown = (f"untestable, {v2b['n_bh_rejected_if_the_20_draw_gate_is_ignored']} if the 20-draw gate is ignored"
                 if v2b["n_bh_rejected_gated"] is None else str(v2b["n_bh_rejected_gated"]))
        lines.append(f"| {ch_names[ch]} | {v1b} of {a['n_observed_in_stratum']} | {shown} |")
    lines += [
        "",
        "## Numbers that are wrong here and are corrected in the superseding package",
        "",
        "- the both-variants-active stratum rejection counts above;",
        "- the cancelled job's call and checkpoint accounting, 45 calls and 11 pairs, not 86 and 21;",
        "- the claim that the prespecification hash was stamped before the first API call;",
        "- the amendment `written_utc` fields, all later than their own file mtimes;",
        "- the both-active headline, which called the stratum untestable above a table marking three of four",
        "  channels testable and mixed two stages' counts;",
        "- the 0.25% by-construction stratum rate, which assumes the two arms are independent;",
        "- the headline sentence attributing the three RNA rejections to overlapping chromatin readout windows.",
        "",
        "Everything else in this package, including F2.1's verdict, F2.2, F3.1, F3.2, F3.3, the EUR haplotype",
        "frequencies and every diagnostic, is unchanged and is restated in the superseding package's RESULTS.md,",
        f"which also records three further defects found while repairing these two: the {dup['rows']}-row family",
        f"holds only {dup['distinct_variant_pairs']} distinct variant pairs, the RNA readout contains neither",
        f"variant in {geo['n_pairs_whose_rna_readout_contains_neither_variant']} of the {dup['rows']} rows, and the",
        f"RNA readout width runs from {gs['min']:.0f} to {gs['max']:.0f} bp against a fixed",
        f"{geo['v1_null_rna_readout_bp_for_2470_of_2500_draws']} bp in almost every null draw.",
        "",
    ]
    (V1 / "SUPERSEDED.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {V1 / 'SUPERSEDED.md'}")


if __name__ == "__main__":
    main()

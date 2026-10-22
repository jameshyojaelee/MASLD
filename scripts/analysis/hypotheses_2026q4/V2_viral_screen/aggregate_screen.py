#!/usr/bin/env python3
"""V2 target screen, final step: per-cohort and per-person rates.

Reads every `libraries/<cohort>/<library>/screen_result.tsv` under the result
directory, applies the batch-surrogate index-hopping flag from PRESPEC_V2.md
section 6, and writes:

  screen_per_library.tsv   one row per library and target, with the flag
  screen_per_cohort.tsv    per-cohort counts and rates with the assessability
                           label from PRESPEC_V2.md section 3
  screen_per_person.tsv    person-level rate for the pooled five cohorts,
                           using the one-per-person roster

Assessability is declared in advance by chemistry, not by outcome. A rate is
not reported for a target and cohort combination declared not assessable.

Usage:
  aggregate_screen.py --result-dir <dir> [--roster <primary_bulk_metadata.tsv>]
"""

from __future__ import annotations

import argparse
import collections
import csv
import glob
import os
import sys

# Enrichment chemistry per cohort, from GEO !Sample_extract_protocol_ch1.
# See PRESPEC_V2.md section 3. "unstated" means GEO does not say whether the
# library was poly-A selected or ribosomally depleted.
ENRICHMENT = {
    "GSE126848": "polyA",
    "GSE130970": "polyA",
    "GSE135251": "unstated",
    "GSE162694": "unstated",
    "GSE213621": "unstated",
    "GSE167523": "polyA",
    "GSE174478": "rRNA_depletion",
    "GSE193066": "exome_capture",
    "GSE240729": "unstated",
    "PRJNA512027": "unstated",
}

POOLED_FIVE = ["GSE126848", "GSE130970", "GSE135251", "GSE162694",
               "GSE213621"]

HOP_FRACTION = 1.0 / 1000.0


def assessability(cohort: str, virus: str) -> str:
    """Prespecified label; never depends on an outcome."""
    chem = ENRICHMENT.get(cohort, "unstated")
    if chem == "exome_capture":
        return "not_assessable_human_exome_capture"
    if virus.startswith("CONTROL"):
        return "control"
    if virus == "HBV":
        # target class B transcripts are polyadenylated.
        return "assessable"
    # target class C and target class E genomes are not polyadenylated.
    if chem == "rRNA_depletion":
        return "assessable"
    if chem == "polyA":
        return "not_assessable_polyA_selection"
    return "sensitivity_unknown_enrichment_unstated"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--result-dir", required=True)
    ap.add_argument("--roster", default=(
        "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/"
        "resource-corrections-20261005-v1/identity_resource_qualified/"
        "bulk_execution_inputs/primary_bulk_metadata.tsv"))
    args = ap.parse_args()

    rows = []
    pattern = os.path.join(args.result_dir, "libraries", "*", "*",
                           "screen_result.tsv")
    for path in sorted(glob.glob(pattern)):
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                r["fragments"] = int(r["fragments"])
                r["criteria_met"] = int(r["criteria_met"])
                r["rpm"] = float(r["rpm"])
                rows.append(r)
    if not rows:
        sys.stderr.write("no screen_result.tsv files under %s\n"
                         % args.result_dir)
        return 2

    # Index-hopping flag: cohort as the batch surrogate, because the
    # SRA-normalized FASTQ headers carry no flowcell or lane.
    max_frag: dict[tuple[str, str], int] = {}
    for r in rows:
        k = (r["cohort"], r["virus"])
        max_frag[k] = max(max_frag.get(k, 0), r["fragments"])
    for r in rows:
        k = (r["cohort"], r["virus"])
        r["assessability"] = assessability(r["cohort"], r["virus"])
        hop = (r["call"] == "positive"
               and max_frag[k] > 0
               and r["fragments"] < HOP_FRACTION * max_frag[k])
        r["index_hop_flag"] = int(hop)
        r["final_call"] = ("flagged_possible_index_hopping_or_carryover"
                           if hop else r["call"])

    out_dir = args.result_dir
    cols = ["cohort", "library", "virus", "assessability", "best_contig",
            "fragments", "unique_sequences", "covered_bases",
            "windows_200bp", "usable_length", "rpm", "criteria_met", "call",
            "index_hop_flag", "final_call", "total_reads_or_pairs"]
    with open(os.path.join(out_dir, "screen_per_library.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: (x["cohort"], x["library"],
                                             x["virus"])):
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # Per-cohort rates.
    agg: dict[tuple[str, str], collections.Counter] = collections.defaultdict(
        collections.Counter)
    for r in rows:
        agg[(r["cohort"], r["virus"])][r["final_call"]] += 1
        agg[(r["cohort"], r["virus"])]["n_libraries"] += 1

    ccols = ["cohort", "virus", "enrichment", "assessability", "n_libraries",
             "positive", "flagged", "indeterminate", "negative",
             "positive_rate_libraries"]
    with open(os.path.join(out_dir, "screen_per_cohort.tsv"), "w") as fh:
        fh.write("\t".join(ccols) + "\n")
        for (cohort, virus), c in sorted(agg.items()):
            lab = assessability(cohort, virus)
            n = c["n_libraries"]
            pos = c["positive"]
            rate = ("%.4f" % (pos / n)) if (n and lab in ("assessable",
                                                          "control")) else "NA"
            fh.write("\t".join(str(x) for x in [
                cohort, virus, ENRICHMENT.get(cohort, "unstated"), lab, n, pos,
                c["flagged_possible_index_hopping_or_carryover"],
                c["indeterminate"], c["negative"], rate]) + "\n")

    # Person-level rate for the pooled five cohorts.
    roster: dict[str, str] = {}
    if os.path.exists(args.roster):
        with open(args.roster) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sid = r.get("sample_id", "").strip()
                unit = r.get("analysis_unit_id", "").strip()
                if sid and unit:
                    roster[sid] = unit

    pcols = ["cohort", "virus", "assessability", "n_persons_screened",
             "positive_persons", "positive_rate_persons"]
    with open(os.path.join(out_dir, "screen_per_person.tsv"), "w") as fh:
        fh.write("\t".join(pcols) + "\n")
        if not roster:
            fh.write("# roster not found at %s; person-level rate not "
                     "computed\n" % args.roster)
        else:
            per: dict[tuple[str, str], dict[str, str]] = \
                collections.defaultdict(dict)
            for r in rows:
                if r["cohort"] not in POOLED_FIVE:
                    continue
                unit = roster.get(r["library"])
                if not unit:
                    continue
                prev = per[(r["cohort"], r["virus"])].get(unit)
                # one library per person in this roster; keep the stronger call
                rank = {"positive": 3,
                        "flagged_possible_index_hopping_or_carryover": 2,
                        "indeterminate": 1, "negative": 0}
                if prev is None or rank[r["final_call"]] > rank[prev]:
                    per[(r["cohort"], r["virus"])][unit] = r["final_call"]
            for (cohort, virus), units in sorted(per.items()):
                lab = assessability(cohort, virus)
                n = len(units)
                pos = sum(1 for v in units.values() if v == "positive")
                rate = ("%.4f" % (pos / n)) if (n and lab in ("assessable",
                                                              "control")) else "NA"
                fh.write("\t".join(str(x) for x in [
                    cohort, virus, lab, n, pos, rate]) + "\n")

    sys.stderr.write("aggregated %d library-target rows from %d libraries\n"
                     % (len(rows), len({(r["cohort"], r["library"])
                                        for r in rows})))
    return 0


if __name__ == "__main__":
    sys.exit(main())

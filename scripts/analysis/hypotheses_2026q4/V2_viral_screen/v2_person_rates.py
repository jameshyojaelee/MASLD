#!/usr/bin/env python3
"""V2 target screen: person-level and library-level rates with intervals.

Reads screen_per_library.tsv written by aggregate_screen.py (positivity rule,
index-hopping flag and assessability labels come from that script and from
score_library.py, unchanged) and writes, into the same analysis directory:

  rates_by_cohort.tsv         counts, rate, Wilson 95% interval, rule-of-three
                              upper bound when no positives, per basis
  person_calls.tsv            one row per person unit and target
  unscreened_libraries.tsv    manifest libraries without a screen result
  controls.tsv                phiX, HPV16 and synthetic-spike summaries
  borderline_positives.tsv    positive libraries near the thresholds (descriptive)

Bases (never added together):
  library                 one row per library
  person_screen_component person = filtered genetic-identity cluster
                          (`screen_component`, all_screened_roster.tsv); a
                          library absent from that roster is its own unit
  person_roster774        pooled five cohorts only, person = analysis_unit_id
                          of the 774-unit primary roster (PRESPEC_V2 section 8)

A person is positive when any screened library of that person is positive
(rank positive > flagged > indeterminate > negative, as in aggregate_screen.py).
A rate is reported only for target x cohort combinations labelled assessable or
control in aggregate_screen.py; counts are reported for all.

Usage:
  v2_person_rates.py --analysis-dir D --manifest M --roster-all R --roster-774 P
                     --aggregator aggregate_screen.py
"""

from __future__ import annotations

import argparse
import collections
import csv
import importlib.util
import math
import os
import sys

Z95 = 1.959963984540054
RANK = {"positive": 3, "flagged_possible_index_hopping_or_carryover": 2,
        "indeterminate": 1, "negative": 0}
TARGETS = ["HBV", "HCV", "HEV"]
CONTROLS = ["CONTROL_phiX", "CONTROL_HPV16"]
SCREENED_CLASSES = TARGETS + CONTROLS
SYNTH = "SYNTHETIC_CONTROL"


def wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1.0 + Z95 * Z95 / n
    centre = (p + Z95 * Z95 / (2.0 * n)) / den
    half = Z95 * math.sqrt(p * (1.0 - p) / n + Z95 * Z95 / (4.0 * n * n)) / den
    return (max(0.0, centre - half), min(1.0, centre + half))


def fmt(x: float) -> str:
    return "NA" if x != x else "%.4f" % x


def read_tsv(path: str) -> list[dict[str, str]]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_aggregator(path: str):
    spec = importlib.util.spec_from_file_location("aggregate_screen", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis-dir", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--roster-all", required=True)
    ap.add_argument("--roster-774", required=True)
    ap.add_argument("--aggregator", required=True)
    a = ap.parse_args()

    agg = load_aggregator(a.aggregator)
    out = a.analysis_dir

    lib_rows = read_tsv(os.path.join(out, "screen_per_library.tsv"))
    manifest = read_tsv(a.manifest)
    roster_all = {r["sample_id"]: r for r in read_tsv(a.roster_all)}
    roster_774 = {r["sample_id"]: r["analysis_unit_id"]
                  for r in read_tsv(a.roster_774)}

    # call per (cohort, library, target)
    call: dict[tuple[str, str, str], str] = {}
    rows_by_key: dict[tuple[str, str, str], dict[str, str]] = {}
    for r in lib_rows:
        key = (r["cohort"], r["library"], r["virus"])
        call[key] = r["final_call"]
        rows_by_key[key] = r
    screened = {(c, l) for (c, l, _v) in call if c != SYNTH}
    manifest_set = {(m["cohort"], m["library"]) for m in manifest}

    # ---------------------------------------------------------- unscreened
    with open(os.path.join(out, "unscreened_libraries.tsv"), "w") as fh:
        fh.write("cohort\tlibrary\n")
        for c, l in sorted(manifest_set - screened):
            fh.write("%s\t%s\n" % (c, l))
    unscreened = manifest_set - screened
    extra = screened - manifest_set
    if extra:
        sys.stderr.write("screened libraries absent from manifest: %d\n"
                         % len(extra))

    # --------------------------------------------------------------- units
    # library basis
    units: dict[str, dict[tuple[str, str], list[str]]] = {
        "library": collections.defaultdict(list),
        "person_screen_component": collections.defaultdict(list),
        "person_roster774": collections.defaultdict(list),
    }
    comp_cohorts: dict[str, set[str]] = collections.defaultdict(set)
    for c, l in sorted(screened):
        units["library"][(c, "LIB::" + l)].append(l)
        ra = roster_all.get(l)
        if c == "GSE193066":
            # all 164 libraries are excluded from the identity roster
            # (repeat biopsies, unqualified crosswalk); no person unit exists
            continue
        if ra and ra.get("screen_component"):
            pid = ra["screen_component"]
            comp_cohorts[pid].add(c)
        else:
            pid = "LIB::" + l
        units["person_screen_component"][(c, pid)].append(l)
        if c in agg.POOLED_FIVE and l in roster_774:
            units["person_roster774"][(c, roster_774[l])].append(l)

    spanning = {k: v for k, v in comp_cohorts.items() if len(v) > 1}

    def unit_call(libs: list[str], cohort: str, virus: str) -> str:
        best = "negative"
        for l in libs:
            c = call.get((cohort, l, virus))
            if c is not None and RANK[c] > RANK[best]:
                best = c
        return best

    # person_calls.tsv
    with open(os.path.join(out, "person_calls.tsv"), "w") as fh:
        fh.write("basis\tcohort\tunit\tvirus\tcall\tn_libraries\tlibraries\n")
        for basis in ("person_screen_component", "person_roster774"):
            for (c, u), libs in sorted(units[basis].items()):
                for v in SCREENED_CLASSES:
                    fh.write("\t".join([basis, c, u, v, unit_call(libs, c, v),
                                        str(len(libs)), ",".join(libs)]) + "\n")

    # ---------------------------------------------------------------- rates
    cols = ["basis", "cohort", "virus", "enrichment", "assessability",
            "n_units", "n_positive", "n_flagged", "n_indeterminate",
            "n_negative", "rate", "wilson95_low", "wilson95_high",
            "rule_of_three_upper"]
    out_rows: list[list[str]] = []

    def add_row(basis: str, cohort: str, virus: str, calls: list[str],
                enrichment: str, lab: str) -> None:
        n = len(calls)
        cnt = collections.Counter(calls)
        pos = cnt["positive"]
        reportable = lab in ("assessable", "control") and n > 0
        if reportable:
            lo, hi = wilson(pos, n)
            rate = pos / n
            r3 = (3.0 / n) if pos == 0 else float("nan")
        else:
            lo = hi = rate = r3 = float("nan")
        out_rows.append([basis, cohort, virus, enrichment, lab, str(n),
                         str(pos),
                         str(cnt["flagged_possible_index_hopping_or_carryover"]),
                         str(cnt["indeterminate"]), str(cnt["negative"]),
                         fmt(rate), fmt(lo), fmt(hi), fmt(r3)])

    cohorts = sorted({c for (c, _u) in units["library"]})
    for basis in ("library", "person_screen_component", "person_roster774"):
        for cohort in cohorts:
            ulist = [(u, libs) for (c, u), libs in units[basis].items()
                     if c == cohort]
            if not ulist:
                continue
            for v in SCREENED_CLASSES:
                calls = [unit_call(libs, cohort, v) for _u, libs in ulist]
                add_row(basis, cohort, v, calls,
                        agg.ENRICHMENT.get(cohort, "unstated"),
                        agg.assessability(cohort, v))
    # pooled five, roster-774 persons, target class B only is assessable everywhere
    pooled = [(u, libs, c) for (c, u), libs in units["person_roster774"].items()]
    for v in SCREENED_CLASSES:
        calls = [unit_call(libs, c, v) for _u, libs, c in pooled]
        labs = {agg.assessability(c, v) for _u, _l, c in pooled}
        lab = "assessable" if labs == {"assessable"} else (
            "control" if labs == {"control"} else
            "mixed_" + "+".join(sorted(labs)))
        # pooled interval only for target class B and controls (assessable in all five)
        enr = "pooled_five_mixed_chemistry"
        if lab not in ("assessable", "control"):
            lab = "not_pooled_" + lab
        add_row("person_roster774", "POOLED_FIVE", v, calls, enr, lab)

    with open(os.path.join(out, "rates_by_cohort.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out_rows:
            fh.write("\t".join(r) + "\n")

    # ------------------------------------------------------------- controls
    ctrl = []
    # HPV16 negative control across all screened libraries
    hpv = collections.Counter(call[(c, l, "CONTROL_HPV16")]
                              for (c, l) in screened
                              if (c, l, "CONTROL_HPV16") in call)
    ctrl.append(["HPV16_negative_control", "all_screened_libraries",
                 str(sum(hpv.values())),
                 ";".join("%s=%d" % kv for kv in sorted(hpv.items()))])
    phix: dict[str, collections.Counter] = collections.defaultdict(
        collections.Counter)
    for (c, l) in screened:
        k = (c, l, "CONTROL_phiX")
        if k in call:
            phix[c][call[k]] += 1
    for c in sorted(phix):
        ctrl.append(["phiX_process_control", c, str(sum(phix[c].values())),
                     ";".join("%s=%d" % kv for kv in sorted(phix[c].items()))])
    # synthetic spike
    syn = [(l, v, call[(c, l, v)]) for (c, l, v) in call
           if c == SYNTH and v in SCREENED_CLASSES]
    for l, v, cl in sorted(syn):
        r = rows_by_key[(SYNTH, l, v)]
        ctrl.append(["synthetic_HBV_spike", l + ":" + v, r["fragments"],
                     "call=%s;covered_bases=%s;windows=%s;rpm=%s" % (
                         cl, r["covered_bases"], r["windows_200bp"], r["rpm"])])
    with open(os.path.join(out, "controls.tsv"), "w") as fh:
        fh.write("control\tscope\tn\tresult\n")
        for r in ctrl:
            fh.write("\t".join(r) + "\n")

    # ---------------------------------------------------- borderline positives
    max_frag: dict[tuple[str, str], int] = {}
    for (c, l, v), r in rows_by_key.items():
        if c == SYNTH:
            continue
        k = (c, v)
        max_frag[k] = max(max_frag.get(k, 0), int(r["fragments"]))
    with open(os.path.join(out, "borderline_positives.tsv"), "w") as fh:
        fh.write("cohort\tlibrary\tvirus\tfragments\tcohort_max_fragments\t"
                 "fraction_of_cohort_max\trpm\tfinal_call\treason\n")
        for (c, l, v), r in sorted(rows_by_key.items()):
            if c == SYNTH or v not in TARGETS:
                continue
            if r["call"] != "positive":
                continue
            frac = int(r["fragments"]) / max_frag[(c, v)] if max_frag[(c, v)] else float("nan")
            reasons = []
            if float(r["rpm"]) < 1.0:
                reasons.append("rpm_below_1")
            if frac < 0.01:
                reasons.append("below_1pct_of_cohort_max")
            if reasons:
                fh.write("\t".join([c, l, v, r["fragments"], str(max_frag[(c, v)]),
                                    fmt(frac), r["rpm"], r["final_call"],
                                    ";".join(reasons)]) + "\n")

    sys.stderr.write(
        "screened libraries %d; unscreened %d; person units (screen_component) %d; "
        "roster774 units %d; components spanning cohorts %d\n" % (
            len(screened), len(unscreened),
            len(units["person_screen_component"]),
            len(units["person_roster774"]), len(spanning)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

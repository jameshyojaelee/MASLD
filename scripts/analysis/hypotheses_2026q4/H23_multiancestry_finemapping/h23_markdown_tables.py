#!/usr/bin/env python3
"""Print the H23 result tables as Markdown from the analysis output directory (no recomputation).

Usage: python3 h23_markdown_tables.py <analysis_dir>
"""
import csv
import json
import sys
from pathlib import Path


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x, nd=3):
    if x in ("", "NA", None):
        return "NA"
    v = float(x)
    if v != 0 and abs(v) < 0.001:
        return f"{v:.1e}"
    return f"{v:.{nd}f}"


def ci(lo, hi, nd=2):
    if lo in ("", "NA") or hi in ("", "NA"):
        return "NA"
    return f"{float(lo):.{nd}f} to {float(hi):.{nd}f}"


def wilson(rate, n, z=1.959964):
    """Wilson 95% interval for a proportion given as rate over n units; closed form, not a bootstrap."""
    if rate in ("", "NA", None) or n in ("", "NA", None) or int(float(n)) == 0:
        return "NA"
    n = int(float(n)); p = float(rate)
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / den
    return f"{p:.2f} ({max(0.0, centre - half):.2f} to {min(1.0, centre + half):.2f})"


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


d = Path(sys.argv[1])
inv = read(d / "inventory.tsv")
tests = read(d / "tests.tsv")
rates = read(d / "rates_by_set_trait.tsv")

print("### Inventory (unit = locus window)\n")
print(table(
    ["Set", "Trait", "RDS loci", "With EUR arm", "After window collapse", "MA fit ok", "EUR-only ok", "Comparable MA", "Comparable baseline", "Paired",
     "MA CS, EUR-only 0 CS", "MA CS, EUR-only unusable", "Note"],
    [[r["set"], r["trait"], r["n_rds"], r["n_with_eur_arm"], r["n_after_window_collapse"], r["n_ma_ok"], r["n_A_ok"], r["n_comparable_MA"],
      r["n_comparable_B"], r["n_paired"], r.get("n_ma_cs_but_no_eur_cs", ""), r.get("n_ma_cs_eur_fit_unusable", ""), r["note"]] for r in inv]))

for fam, title in (("primary", "Primary family (paired by locus; BH over 12 tests)"),
                   ("secondary", "Secondary family (pooled ALT+AST baseline, unpaired; BH over 9 tests)")):
    rows = [r for r in tests if r["family"] == fam]
    print(f"\n### {title}\n")
    hdr = ["Set", "Trait", "Endpoint", "n MA", "n base", "n paired", "MA", "Baseline", "MA minus baseline (95% CI)", "Discordant MA-only/base-only",
           "p one-sided", "q BH", "Status"]
    out = []
    for r in rows:
        if r["status"] != "tested":
            out.append([r["set"], r["trait"], r["endpoint"], r["n_MA"], r["n_B"], r["n_paired"], "", "", "", "", "", "", r["status"]])
            continue
        disc = "" if r["disc_MA_only"] in ("", "NA") else f'{r["disc_MA_only"]}/{r["disc_B_only"]}'
        out.append([r["set"], r["trait"], r["endpoint"], r["n_MA"], r["n_B"], r["n_paired"], num(r["rate_MA"]), num(r["rate_B"]),
                    f'{num(r["effect"], 3)} ({ci(r["ci_lo"], r["ci_hi"])})', disc, num(r["p"]), num(r["q_BH"]), "tested"])
    print(table(hdr, out))

for fam in sorted({r["family"] for r in tests if r["family"].startswith(("sensitivity", "exploratory"))}):
    rows = [r for r in tests if r["family"] == fam]
    print(f"\n### {fam} (nominal p, outside BH)\n")
    hdr = ["Set", "Trait", "Endpoint", "n MA", "n base", "n paired", "MA", "Baseline", "MA minus baseline (95% CI)", "Discordant", "p one-sided", "p unpaired Fisher", "Status"]
    out = []
    for r in rows:
        if r["status"] != "tested":
            out.append([r["set"], r["trait"], r["endpoint"], r["n_MA"], r["n_B"], r["n_paired"], "", "", "", "", "", "", r["status"]])
            continue
        disc = "" if r["disc_MA_only"] in ("", "NA") else f'{r["disc_MA_only"]}/{r["disc_B_only"]}'
        out.append([r["set"], r["trait"], r["endpoint"], r["n_MA"], r["n_B"], r["n_paired"], num(r["rate_MA"]), num(r["rate_B"]),
                    f'{num(r["effect"], 3)} ({ci(r["ci_lo"], r["ci_hi"])})', disc, num(r["p"]), num(r["p_unpaired"]), "tested"])
    print(table(hdr, out))

print("\n### Descriptive rates per set and trait (comparable units)\n")
print(table(["Set", "Trait", "n MA", "n base", "Tie at max PIP, MA", "Tie at max PIP, base", "Lead replaced MA (Wilson 95%)", "Lead replaced base (Wilson 95%)",
             "Strict-identity change MA", "Strict-identity change base", "Median log2 CS ratio MA", "Median log2 CS ratio base"],
            [[r["set"], r["trait"], r["n_cmp_MA"], r["n_cmp_B"], num(r["tie_frac_MA"], 2), num(r["tie_frac_B"], 2),
              wilson(r["replaced_rate_MA"], r["n_cmp_MA"]), wilson(r["replaced_rate_B"], r["n_cmp_B"]),
              num(r["strict_rate_MA"], 2), num(r["strict_rate_B"], 2), num(r["median_log2_cs_ratio_MA"], 2), num(r["median_log2_cs_ratio_B"], 2)] for r in rates]))

# lead replaced by the number of credible sets in the EUR-only reference fit (descriptive, no test)
units = read(d / "units.tsv")
arms = [("MA", "rep_MA", "cmp_MA"), ("independent EUR", "rep_B", "cmp_B"), ("nested EUR, matched (Mm)", "rep_Mm", "cmp_Mm"), ("nested EUR, all of B (Mf)", "rep_Mf", "cmp_Mf")]
if units and "rep_Mm" in units[0]:
    print("\n### Lead replaced by number of credible sets in the EUR-only reference fit (units where all four comparisons are evaluable)\n")
    out = []
    for st in ("MVP", "4WAY"):
        for tr in ("ALT", "AST"):
            common = [u for u in units if u["set"] == st and u["trait"] == tr and all(u[c] == "TRUE" for _, _, c in arms)]
            for label, lo, hi in (("1 CS", 1, 1), ("2 or more CS", 2, 99)):
                g = [u for u in common if lo <= int(float(u["A_n_cs"])) <= hi]
                cells = [f'{sum(u[r] == "TRUE" for u in g)}/{len(g)}' for _, r, _ in arms]
                out.append([st, tr, label, len(g)] + cells)
    print(table(["Set", "Trait", "EUR-only reference CS count", "n units"] + [a for a, _, _ in arms], out))

print("\n### Controls\n")
ctl = read(d / "controls.tsv")
for r in ctl:
    print(f'- **{r["control"]}**: `{r["json"]}`')

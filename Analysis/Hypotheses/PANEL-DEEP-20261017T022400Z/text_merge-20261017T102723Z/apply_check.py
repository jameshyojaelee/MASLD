"""Ran against the pre-apply live files (all passed); the 137 edits were applied by team-lead on
2026-10-17 06:56, so it no longer runs against live text. Its 7E checks moved to followup_check.py.

Copy the live text files to mirror/ under this directory, apply edits.json in
order, and check the result. Live files are only read.

Checks: each `old` occurs exactly once in the mirror at the moment it is
applied; the edit count per file matches edits.json; no edit adds an em dash or
the words frozen, progression or cascade; key numbers printed by the new text
are re-read from their producing source tables. Writes mirror_diff.txt.
"""
import csv
import difflib
import json
import os
import re
import shutil
from collections import Counter

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/text_merge-20261017T102723Z")
CAND = os.path.join(ROOT, "figures/candidates/panel-deep-fixes-20261017T035633Z")
MIRROR = os.path.join(OUT, "mirror")
BANNED = re.compile(r"—|\bfrozen\b|\bprogression\b|\bcascade\b", re.IGNORECASE)

edits = json.load(open(os.path.join(OUT, "edits.json"), encoding="utf-8"))
files = sorted({e["file"] for e in edits})

if os.path.isdir(MIRROR):
    shutil.rmtree(MIRROR)  # mirror is a scratch copy owned by this script
original, text = {}, {}
for f in files:
    dst = os.path.join(MIRROR, f)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(os.path.join(ROOT, f), dst)
    original[f] = text[f] = open(dst, encoding="utf-8").read()

failures = []
applied = Counter()
for e in edits:
    n = text[e["file"]].count(e["old"])
    if n != 1:
        failures.append(f"{e['id']}: old occurs {n} times in {e['file']} when applied")
        continue
    added = len(BANNED.findall(e["new"])) - len(BANNED.findall(e["old"]))
    if added > 0:
        failures.append(f"{e['id']}: adds {added} banned term(s) or em dash: {BANNED.findall(e['new'])}")
    text[e["file"]] = text[e["file"]].replace(e["old"], e["new"], 1)
    applied[e["file"]] += 1

expected = Counter(e["file"] for e in edits)
assert applied == expected, (applied, expected)
assert sum(applied.values()) == len(edits) == 137, sum(applied.values())

diff_lines = []
for f in files:
    with open(os.path.join(MIRROR, f), "w", encoding="utf-8") as fh:
        fh.write(text[f])
    diff_lines += difflib.unified_diff(original[f].splitlines(keepends=True), text[f].splitlines(keepends=True),
                                       fromfile=f"live/{f}", tofile=f"mirror/{f}", n=0)
open(os.path.join(OUT, "mirror_diff.txt"), "w", encoding="utf-8").writelines(diff_lines)

# Banned terms remaining anywhere in the edited files (pre-existing text is reported, not failed).
residual = {f: len(BANNED.findall(text[f])) for f in files}
residual_before = {f: len(BANNED.findall(original[f])) for f in files}


# --------------------------------------------------------------------------
# Re-read key printed numbers from their producing tables.
# --------------------------------------------------------------------------
def rows(rel):
    with open(os.path.join(CAND, rel), encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


checks = []


def check(label, got, want):
    checks.append((label, got, want, got == want))


c1 = rows("fig1/source_tables/fig1c.tsv")
enz = {r["state"]: int(r["n"]) for r in c1 if r["trait_class"].startswith("liver_enzyme")}
dir_ = {r["state"]: int(r["n"]) for r in c1 if r["trait_class"].startswith("direct")}
check("1C enzyme supported/tested/none", (enz["supported"], enz["tested"], enz["none"]), (10, 46, 368))
check("1C direct supported/tested/none", (dir_["supported"], dir_["tested"], dir_["none"]), (10, 11, 1))
check("1C tested denominators", sorted({int(r["class_tested_genes"]) for r in c1}), [392, 3888])
cd = {r["quantity"]: int(r["n"]) for r in rows("fig1/source_tables/fig1c_common_denominator.tsv")}
check("1C common denominator", (cd["genes SuSiE-tested in both classes"], cd["enzyme-supported among them"],
                                cd["direct-supported among them"]), (362, 56, 21))

c2 = {r["class"]: int(r["n_genes"]) for r in rows("fig2/source_tables/2C.tsv")}
check("2C counts", [c2[k] for k in c2], [14, 21, 70, 104, 524, 263])
check("2C distinct leads", {r["n_distinct_lead_variants"] for r in rows("fig2/source_tables/2C.tsv")}, {"933"})

e2 = rows("fig2/source_tables/2E.tsv")
check("2E EUR-only split", [(int(r["n_genes"]), r["of_which_nonEUR_ABF_gt_0.5"]) for r in e2[:2]],
      [(127, "12.0"), (213, "22.0")])

h = rows("fig4/source_tables/4H.tsv")
check("4H donor-first rho, same direction", (round(float(h[0]["rho_donor_first"]), 3),
                                              int(h[0]["same_direction_donor_first"])), (0.670, 34))
check("4H feature-first rho, same direction", (round(float(h[0]["rho_feature_first_installed"]), 3),
                                                int(h[0]["same_direction_feature_first_installed"])), (0.631, 32))
check("4H bootstrap interval", h[0]["rho_donor_first_boot_95"], "0.382 to 0.842")
st = Counter(r["state"] for r in h)
check("4H states", (st["Both positive"], st["Both negative"], st["Not resolved"], len(h)), (19, 2, 21, 43))

e5 = rows("fig5/source_tables/5E.tsv")
open_ = [r for r in e5 if int(r["n_open_genes"] or 0) > 0]
robust = [r for r in e5 if r["robust"] == "TRUE"]
testable = [r for r in e5 if r["testable"] == "TRUE"]
check("5E lineage-open programs", len(open_), 31)
check("5E direct-open programs", sum(int(r["n_open_direct_genes"] or 0) > 0 for r in e5), 3)
check("5E testable programs", len(testable), 95)
check("5E robust (neg, pos)", (len(robust), sum(float(r["rho"]) < 0 for r in robust),
                               sum(float(r["rho"]) > 0 for r in robust)), (22, 14, 8))
check("5E both", sum(r in robust for r in open_), 8)
chol5 = [r for r in robust if "ETFB" in r["program_name"]]
check("5E ETFB boundary P, q", [(r["pvalue"], round(float(r["padj"]), 5)) for r in chol5], [("0.0121", 0.04998)])
efm = [r for r in e5 if r["program_name"].startswith("Elastic-fiber ECM")]
check("5E Elastic-fiber ECM measured", [r["n_measured"] for r in efm], ["7"])

g = rows("fig3/source_tables/fig3g.tsv")
uni = {r["universe"]: r for r in g if r["universe"]}
u69 = [uni["testable_3769"]]
u76 = [uni["all_finite_3776"]]
check("3G 3,769 classes", [tuple(int(float(r[k])) for k in ("neither", "tx_only", "gen_only", "both")) for r in u69],
      [(3149, 194, 400, 26)])
check("3G 3,769 rho, P, OR, P, expected",
      [(round(float(r["spearman_rho"]), 3), round(float(r["spearman_p"]), 2), round(float(r["fisher_or"]), 2),
        round(float(r["fisher_p"]), 2), round(float(r["expected_both"]), 1)) for r in u69],
      [(0.001, 0.95, 1.06, 0.83, 24.9)])
check("3G 3,776 rho, P, OR, P",
      [(round(float(r["spearman_rho"]), 4), round(float(r["spearman_p"]), 2), round(float(r["fisher_or"]), 2),
        round(float(r["fisher_p"]), 2)) for r in u76], [(0.0014, 0.93, 1.06, 0.74)])

f5 = {(r["program"], r["dataset"]): r for r in rows("fig5/source_tables/5F.tsv") if r["mark"] == "pooled_point"}
s5 = f5[("Stromal ECM", "GSE192741")]
check("5F Stromal ECM donors positive, sign-flip P", (s5["n_units_positive"], round(float(s5["signflip_p"]), 2)),
      ("3.0", 0.19))
check("5F Stromal ECM Vu arrays positive", f5[("Stromal ECM", "Vu_et_al_2025")]["n_units_positive"], "10.0")

# 7B and 7D producer values (team-lead 2026-10-17: legend uses producer values, not the redraw recomputes).
iv = json.load(open(os.path.join(ROOT, "Analysis/MASLD_Model_Benchmark/executions/dnalm_heads-20260831T183647Z/logs/03_intervals.json")))
pdz = iv["paired_deltas"]["abs_z"]
for model, want in (("hyenadna", (-0.022, -0.033, -0.010)), ("nucleotide_transformer", (-0.029, -0.042, -0.016))):
    r = pdz[f"{model}::m1_dist_emb  MINUS  {model}::b1_dist"]
    check(f"7B eQTL {model} point, ci", (round(r["point"], 3), round(r["ci"][0], 3), round(r["ci"][1], 3)), want)
check("7B eQTL seed, resamples", (iv["seed"], iv["bootstrap_resamples"]), (20260831, 2000))
with open(os.path.join(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F7a/exec/mpra/mpra_bonferroni.tsv")) as fh:
    mp = {r["model"]: r for r in csv.DictReader(fh, delimiter="\t") if r["family"] == "F7_displayed_rows"}
check("7B MPRA HyenaDNA gain, Bonferroni lo, hi", tuple(round(float(mp["hyenadna"][k]), 3) for k in ("gain", "bonf_pct_lo", "bonf_pct_hi")),
      (0.082, 0.009, 0.158))
check("7B MPRA NT gain, Bonferroni lo, hi", tuple(round(float(mp["nucleotide_transformer"][k]), 3) for k in ("gain", "bonf_pct_lo", "bonf_pct_hi")),
      (0.054, -0.018, 0.130))
with open(os.path.join(ROOT, "figures/candidates/review-r3-fixes-20261015T112148Z/F7b/repair/fig7_transfer_primary_units.tsv")) as fh:
    sev = [r for r in csv.DictReader(fh, delimiter="\t") if r["cohort"] == "GSE268273" and r["primary"] == "True"][0]
check("7D GSE268273 interval, B, seed", (round(float(sev["ci_lo"]), 4), round(float(sev["ci_hi"]), 4), sev["ci_B"], sev["ci_seed"]),
      (0.5685, 0.7483, "10000", "20260923"))

# Signed eQTL (7B, 7D, 7E) from eqtl_sign_reliability 03_strata_intervals.json (seed 20260831).
sr = json.load(open(os.path.join(ROOT, "Analysis/MASLD_Model_Benchmark/executions/eqtl_sign_reliability-20260831T195702Z/logs/03_strata_intervals.json")))
sa = sr["strata"]["all"]
check("7B signed eQTL AUC, lo, hi, chance (3 dp as printed)",
      tuple(round(x, 3) for x in (sa["best_auc"], *sa["best_auc_ci"], sa["chance_measured"]["median"])), (0.502, 0.490, 0.513, 0.500))
check("7B allele ridge sign AUC", round(sr["headline_provenance"]["allele_feature_ridge_sign_auc_recomputed"], 3), 0.503)
check("7B eQTL seed, resamples", (sr["seed"], sr["bootstrap_resamples"]), (20260831, 2000))
bound = sa["sequence_attributable_bound"]["C_sequence_attributable_upper_bound"]
check("7E bound above chance (AUC)", round(bound - sa["chance_measured"]["median"], 3), 0.006)
check("7E bound / threshold", round((bound - sa["chance_measured"]["median"]) / (2 * sa["resolution_half_width"]), 2), 0.24)
with open(os.path.join(CAND, "fig7/source_tables/fig7E.producer.tsv"), encoding="utf-8") as fh:
    e7 = {r["mark"]: r for r in csv.DictReader(fh, delimiter="\t") if r["row_label"].startswith("signed eQTL")}
est = (sa["best_auc"] - sa["chance_measured"]["median"]) / (2 * sa["resolution_half_width"])
check("7E signed estimate: producer table vs recomputed from 03_strata_intervals",
      (round(float(e7["estimate"]["x"]), 4), round(float(e7["estimate"]["lo"]), 4), round(float(e7["estimate"]["hi"]), 4)),
      (round(est, 4), -0.4168, 0.5832))
check("7E signed bound: producer table", round(float(e7["bound"]["x"]), 4), round((bound - sa["chance_measured"]["median"]) / (2 * sa["resolution_half_width"]), 4))
check("7E signed block and source verdict", (e7["estimate"]["block"], sa["verdict"]), ("untestable", "tested_negative"))

want_strata = {"all": (0.4901, 0.5133, 0.4998), "z_ge_1": (0.4830, 0.5199, 0.5003), "z_ge_2": (0.4822, 0.5143, 0.4999),
               "z_ge_3": (0.4581, 0.5344, 0.5003), "z_ge_5": (0.4688, 0.5121, 0.5019)}
for k, w in want_strata.items():
    s_ = sr["strata"][k]
    check(f"7D eQTL stratum {k} interval, chance", tuple(round(x, 4) for x in (*s_["best_auc_ci"], s_["chance_measured"]["median"])), w)
with open(os.path.join(ROOT, "figures/candidates/review-r3-fixes-20261015T112148Z/F7b/repair/fig7_transfer_primary_units.tsv")) as fh:
    prim = [r for r in csv.DictReader(fh, delimiter="\t") if r["primary"] == "True"]
check("7D severity primary intervals", [(r["cohort"], round(float(r["ci_lo"]), 4), round(float(r["ci_hi"]), 4)) for r in prim],
      [("GSE268273", 0.5685, 0.7483), ("GSE213621", 0.602, 0.7267), ("GSE276114", 0.6985, 0.8483),
       ("GSE276114", 0.3752, 0.6878)])

# Markers that scripts/manuscript/validate_resource_scope.py requires must survive the edits.
contract = json.load(open(os.path.join(ROOT, "config/resource_paper_scope.config.json")))
norm = lambda s: " ".join(s.split())
markers = {"docs/PAPER.md": ["## Ultimate goal", "## The central scientific contribution", "## Seven-figure story",
                             "## Claims we can make", "## Claims we cannot make", "## Two-paper firewall",
                             "## Release acceptance", contract["ultimate_goal"], contract["paper_title"],
                             "The paper records **liver-trait genetic evidence** and **disease-state remodeling** in human MASLD side by side for each gene."]
                            + [f"{i}. " + s for i, s in enumerate(["A multimodal human MASLD Resource", "Inherited genetic variation",
                               "Multicohort disease remodeling", "Multicellular programs", "Molecular and physical tissue context",
                               "MASLD Gene Catalog", "Model framework and transfer limits"], 1)],
           "docs/manuscript/draft/fig5_discussion.md": ["deterministic next-experiment rule"],
           "docs/manuscript/05_figure_legends.md": ["machine-readable evidence state"],
           "figures/main/INDEX.md": ["Shared Figure 1–7 decision rules"]}
for f in ("fig1", "fig2", "fig3", "fig4", "fig5"):
    markers[f"docs/manuscript/draft/{f}.md"] = ["next discriminating experiment"]
for f, ms in markers.items():
    if f not in text:
        continue
    for m in ms:
        before, after = norm(m) in norm(original[f]), norm(m) in norm(text[f])
        if before and not after:
            failures.append(f"validator marker removed from {f}: {m[:60]}")
    checks.append((f"validator markers kept in {f}", len(ms), len(ms), True))

with open(os.path.join(OUT, "value_checks.tsv"), "w", encoding="utf-8") as fh:
    fh.write("check\tread_from_table\tprinted_in_text\tagree\n")
    for c in checks:
        fh.write("\t".join(map(str, c)) + "\n")

print(f"applied {sum(applied.values())} edits")
for f in files:
    print(f"  {applied[f]:3d}  {f}  (banned terms or em dashes in file: {residual_before[f]} before, {residual[f]} after)")
print(f"value checks: {sum(c[3] for c in checks)}/{len(checks)} agree")
for c in checks:
    if not c[3]:
        print("  DISAGREE", c)
print("failures:", failures if failures else "none")
assert not failures

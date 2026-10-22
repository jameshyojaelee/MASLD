"""Follow-up to the 137 applied edits (team-lead, 2026-10-17): the 7E signed-eQTL row.

Writes edits_followup.json with exact old/new pairs against the CURRENT live files
(after the 137 edits were applied), applies them to followup_mirror/, checks that
each old occurs once, that no edit adds an em dash or a retired term, that the
validator markers survive, and re-reads the 7E values from their producer files.
Writes followup_diff.txt and followup_value_checks.tsv. Live files are only read.
"""
import csv
import difflib
import json
import os
import re
import shutil

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/text_merge-20261017T102723Z")
CAND = os.path.join(ROOT, "figures/candidates/panel-deep-fixes-20261017T035633Z")
MIRROR = os.path.join(OUT, "followup_mirror")
BANNED = re.compile(r"—|\bfrozen\b|\bprogression\b|\bcascade\b", re.IGNORECASE)
LEG = "docs/manuscript/05_figure_legends.md"
STRATA = "Analysis/MASLD_Model_Benchmark/executions/eqtl_sign_reliability-20260831T195702Z/logs/03_strata_intervals.json"

edits = [dict(
    id="F-7E-verdicts", file=LEG,
    old="Verdicts are those of the producing analyses.",
    new="Verdicts are those of the producing analyses, except signed eQTL direction (below).",
    reason="The signed-eQTL row is drawn untestable although its source records tested negative; the general "
           "sentence must not say every verdict is the producer's.",
    decision="3", source_ids=["M-7E", "team-lead 2026-10-17 7E producer"],
    producing_file=f"{STRATA} (strata.all.verdict)",
), dict(
    id="F-7E-verdict", file=LEG,
    old="Signed eQTL direction is drawn as untestable.",
    new="The source analysis records signed eQTL direction as tested negative, but its sequence-attributable bound "
        "is below the threshold, so under this panel's rule it is drawn as untestable.",
    reason="7E signed-eQTL row redrawn from fig7E.producer.tsv: estimate 0.0779 [-0.4168, 0.5832], bound 0.2436 "
           "(earlier recompute 0.0851 [-0.4099, 0.5901], 0.2596). The legend prints the bound as 0.24 (already "
           "applied) and does not print the estimate; this edit names the source verdict and the panel rule.",
    decision="3", source_ids=["M-7E", "team-lead 2026-10-17 7E producer"],
    producing_file=f"figures/candidates/panel-deep-fixes-20261017T035633Z/fig7/source_tables/fig7E.producer.tsv; {STRATA} (strata.all.verdict)",
)]
json.dump(edits, open(os.path.join(OUT, "edits_followup.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)

if os.path.isdir(MIRROR):
    shutil.rmtree(MIRROR)  # scratch copy owned by this script
files = sorted({e["file"] for e in edits})
original, text = {}, {}
for f in files:
    dst = os.path.join(MIRROR, f)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(os.path.join(ROOT, f), dst)
    original[f] = text[f] = open(dst, encoding="utf-8").read()

failures = []
for e in edits:
    n = text[e["file"]].count(e["old"])
    if n != 1:
        failures.append(f"{e['id']}: old occurs {n} times")
        continue
    if len(BANNED.findall(e["new"])) > len(BANNED.findall(e["old"])):
        failures.append(f"{e['id']}: adds a banned term or em dash")
    text[e["file"]] = text[e["file"]].replace(e["old"], e["new"], 1)

norm = lambda s: " ".join(s.split())
for f in files:
    if f == LEG and norm("machine-readable evidence state") not in norm(text[f]):
        failures.append("validator marker removed from legends")
    with open(os.path.join(MIRROR, f), "w", encoding="utf-8") as fh:
        fh.write(text[f])
diff = []
for f in files:
    diff += difflib.unified_diff(original[f].splitlines(True), text[f].splitlines(True),
                                 fromfile=f"live/{f}", tofile=f"followup_mirror/{f}", n=0)
open(os.path.join(OUT, "followup_diff.txt"), "w", encoding="utf-8").writelines(diff)

# Value checks: the 7E producer row, recomputed from the strata producer.
checks = []


def check(label, got, want):
    checks.append((label, got, want, got == want))


sa = json.load(open(os.path.join(ROOT, STRATA)))["strata"]["all"]
chance, res = sa["chance_measured"]["median"], sa["resolution_half_width"]
bound = sa["sequence_attributable_bound"]["C_sequence_attributable_upper_bound"]
scale = lambda auc: (auc - chance) / (2 * res)
with open(os.path.join(CAND, "fig7/source_tables/fig7E.producer.tsv"), encoding="utf-8") as fh:
    row = {r["mark"]: r for r in csv.DictReader(fh, delimiter="\t") if r["row_label"].startswith("signed eQTL")}
check("7E estimate, lo, hi (table)", tuple(round(float(row["estimate"][k]), 4) for k in ("x", "lo", "hi")),
      (0.0779, -0.4168, 0.5832))
check("7E estimate recomputed from strata", round(scale(sa["best_auc"]), 4), round(float(row["estimate"]["x"]), 4))
check("7E interval recomputed from strata", tuple(round(scale(x), 4) for x in sa["best_auc_ci"]),
      (round(float(row["estimate"]["lo"]), 4), round(float(row["estimate"]["hi"]), 4)))
check("7E bound (table) and recomputed", (round(float(row["bound"]["x"]), 4), round(scale(bound), 4)), (0.2436, 0.2436))
check("7E bound as printed in legend", f"{scale(bound):.2f}", "0.24")
check("7E bound below threshold 1", scale(bound) < 1, True)
check("7E block; source verdict", (row["estimate"]["block"], sa["verdict"]), ("untestable", "tested_negative"))
check("legend prints 0.24 of the threshold", "is 0.24 of the threshold" in text[LEG], True)
check("legend does not print the old 0.26 or 0.085", ("0.26 of the threshold" in text[LEG], "0.085" in text[LEG]),
      (False, False))

with open(os.path.join(OUT, "followup_value_checks.tsv"), "w", encoding="utf-8") as fh:
    fh.write("check\tread_from_file\texpected\tagree\n")
    for c in checks:
        fh.write("\t".join(map(str, c)) + "\n")
print(f"follow-up edits: {len(edits)}")
print(f"value checks: {sum(c[3] for c in checks)}/{len(checks)} agree")
for c in checks:
    if not c[3]:
        print("  DISAGREE", c)
print("failures:", failures if failures else "none")
assert not failures

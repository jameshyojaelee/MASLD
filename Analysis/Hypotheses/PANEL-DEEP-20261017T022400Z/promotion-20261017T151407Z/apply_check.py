"""Promotion text for release 2026-10-17-r2 (user approval relayed by team-lead, 2026-10-17).

Defines the edits (item 8 placement of the Figure 2D companions, ROADMAP final
release checklist, STATUS release paragraph and retired-checker row) and the
file moves, writes edits.json and file_ops.tsv, copies the current live files to
mirror/, applies the edits in order and checks them: each old occurs once, no
overlaps, no em dash or retired term added, validator marker strings kept, and
the move sources exist and targets do not. Writes mirror_diff.txt. Live files
are only read; no file is moved.
"""
import csv
import difflib
import hashlib
import json
import os
import re
import shutil
from collections import Counter, defaultdict

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/promotion-20261017T151407Z")
MIRROR = os.path.join(OUT, "mirror")
BANNED = re.compile(r"—|\bfrozen\b|\bprogression\b|\bcascade\b|\bpassports?\b", re.IGNORECASE)

LEG = "docs/manuscript/05_figure_legends.md"
RES = "docs/RESULTS.md"
STA = "docs/STATUS.md"
ROAD = "docs/ROADMAP.md"
FIG2MD = "docs/manuscript/draft/fig2.md"
IDX = "figures/main/INDEX.md"
F2R = "figures/main/fig2_genetics/README.md"
P2 = "figures/main/fig2_genetics/panels"
MAN = f"{P2}/Figure2_reindex_manifest.tsv"
REL = "RNA-seq/results/manuscript_release/2026-10-17-r2/release_manifest.tsv"
DEEP = "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z"

_live = {}


def live(path):
    if path not in _live:
        with open(os.path.join(ROOT, path), encoding="utf-8", newline="") as fh:
            _live[path] = fh.read()
    return _live[path]


def locate(path, text):
    pat = r"\s+".join(re.escape(p) for p in text.split())
    hits = [m.group(0) for m in re.finditer(pat, live(path))]
    if len(hits) != 1:
        raise SystemExit(f"locate: {len(hits)} matches in {path} for {text[:80]!r}")
    return hits[0]


EDITS = []


def add(id_, path, old, new, reason, decision, producing_file, exact=False):
    EDITS.append(dict(id=id_, file=path, old=old if exact else locate(path, old), new=new, reason=reason,
                      decision=decision, source_ids=["team-lead 2026-10-17 promotion request"],
                      producing_file=producing_file))


# --------------------------------------------------------------------------
# File moves (item 8): the direct-subclass companion goes to the Figure S2 set,
# which lives in the same panels/ directory under the FigS2<letter>_ prefix
# (FigS2A to FigS2R exist; S is the next letter).
# --------------------------------------------------------------------------
FILE_OPS = [
    (f"{P2}/Fig2D_direct_subclass_split.pdf", f"{P2}/FigS2S_direct_subclass_split.pdf"),
    (f"{P2}/data/panel_deep_20261017/2D-subclass.tsv", f"{P2}/FigS2S_direct_subclass_split_source.tsv"),
]

sha_anc = hashlib.sha256(open(os.path.join(ROOT, P2, "Fig2D_coding_by_ancestry.pdf"), "rb").read()).hexdigest()

# --------------------------------------------------------------------------
# Item 8
# --------------------------------------------------------------------------
add("P8-manifest", MAN,
    "Phenotype and ancestry provenance\r\n",
    "Phenotype and ancestry provenance\r\n"
    f"D_ancestry\t\tFig2D_coding_by_ancestry.pdf\t{sha_anc}\t"
    "Ancestry-stratified coding/splice posterior fractions; companion to D, adopted 2026-10-17\r\n",
    "Item 8: the ancestry companion stays in main Figure 2D and joins the Figure 2 panel roster.", "item 8",
    f"{P2}/Fig2D_coding_by_ancestry.pdf (sha256 computed here)", exact=True)

subclass_par = locate(LEG,
    "**(D, direct-subclass panel)** Three EUR PDFF GWAS (37 sets, 12 clusters) 40.14% and NAFLD plus NASH "
    "(170 sets, 13 clusters) 12.44%, each against all enzyme traits (3.87%): differences 36.3 (18.3 to 54.6) and "
    "8.6 (−1.6 to 21.0) percentage points, 95% intervals from 10,000 physical-cluster draws, seed 20260923. PDFF "
    "minus NAFLD plus NASH is 27.7 (7.2 to 47.7; subclasses resampled independently, 20,000 draws, seed "
    "20261017). Because all PDFF sets are EUR, PDFF is also compared with EUR enzyme sets only (82 sets, 14 "
    "clusters, 6.92%): 33.2 (14.4 to 52.2; 10,000 draws, seed 20261017). NAFLD plus NASH sets are 69 EUR, 95 "
    "AMR, 5 AFR and 1 EAS; enzyme sets are 82 EUR, 715 AMR, 1,664 AFR, 3,694 EAS and 53 SAS. ")
_end = live(LEG).find(subclass_par) + len(subclass_par)
if live(LEG)[_end] == " ":
    subclass_par += " "  # remove the separating space too, so the next sentence keeps its line start
moved_body = subclass_par.replace("**(D, direct-subclass panel)** ", "", 1).rstrip()
EDITS.append(dict(
    id="P8-legend-remove", file=LEG, old=subclass_par, new="",
    reason="Item 8: the direct-subclass paragraph moves to the Figure S2 legend (P8-legend-S2); text unchanged.",
    decision="item 8", source_ids=["team-lead 2026-10-17 promotion request"],
    producing_file=f"{P2}/FigS2S_direct_subclass_split_source.tsv (moved from data/panel_deep_20261017/2D-subclass.tsv)"))

add("P8-legend-S2", LEG,
    "## Figure 3. Bulk transcriptomic remodeling and genetic overlap",
    "## Supplementary Figure S2. Genetics detail\n\n"
    "**(S2S) Direct-subclass split of the Figure 2D coding/splice posterior fraction.** "
    f"{moved_body} The ancestry stratification is the Figure 2D ancestry companion.\n\n"
    "## Figure 3. Bulk transcriptomic remodeling and genetic overlap",
    "Item 8: new S2 legend section holding the moved direct-subclass paragraph (no S2 legend existed in this "
    "file; panels FigS2A to FigS2R have no legend text here).", "item 8",
    f"{P2}/FigS2S_direct_subclass_split_source.tsv")

add("P8-fig2md", FIG2MD,
    "Within the direct class, the three European PDFF GWAS",
    "Within the direct class (Figure S2S), the three European PDFF GWAS",
    "Item 8: the subclass values are cited to their new supplementary panel.", "item 8",
    f"{P2}/FigS2S_direct_subclass_split_source.tsv")

add("P8-results", RES,
    "**Ancestry and subclass stratification (H58 and review round 3, 2026-10-15; not part of the adopted table "
    "above).**",
    "**Ancestry and subclass stratification (H58 and review round 3, 2026-10-15).** The ancestry stratification "
    "(H58 family C) was adopted on 2026-10-17 as the Figure 2D ancestry companion "
    "(`figures/main/fig2_genetics/panels/Fig2D_coding_by_ancestry.pdf`); the direct-subclass split is "
    "supplementary (Figure S2S, `figures/main/fig2_genetics/panels/FigS2S_direct_subclass_split.pdf`). Neither "
    "is part of the adopted ratio table above, and no number changes.",
    "Item 8: release state of H58 family C set to adopted 2026-10-17; subclass split set to supplementary.",
    "item 8", "Analysis/Hypotheses/H58-20261015T112423Z/jobs/5800/H58-20261015T113519Z_22252009/family_C.tsv")

add("P8-index", IDX,
    "verified r2 A–I panels installed; D retains adopted coordinate correction;",
    "verified r2 A–I panels installed; D retains adopted coordinate correction; D ancestry companion adopted "
    "2026-10-17; direct-subclass split moved to Figure S2S;",
    "Item 8.", "item 8", MAN)

add("P8-fig2readme", F2R,
    "The A-I roster is checksummed in `panels/Figure2_reindex_manifest.tsv`.",
    "The A-I roster and the D ancestry companion (`panels/Fig2D_coding_by_ancestry.pdf`, adopted 2026-10-17) "
    "are checksummed in `panels/Figure2_reindex_manifest.tsv`. The direct-subclass split is supplementary, at "
    "`panels/FigS2S_direct_subclass_split.pdf` with its source table "
    "`panels/FigS2S_direct_subclass_split_source.tsv`.",
    "Item 8.", "item 8", MAN)

# --------------------------------------------------------------------------
# ROADMAP final release checklist
# --------------------------------------------------------------------------
add("R-release-id", ROAD,
    "- [ ] One release ID names every paper output file.",
    "- [x] One release ID names every paper output file: release 2026-10-17-r2, with manifest\n"
    f"      `{REL}`.\n"
    "      Release 2026-10-17-r1 was superseded the same day by r2 after text edits.",
    "Checklist: release ID.", "promotion", REL)
add("R-read-only", ROAD,
    "- [ ] Source files, environments, scripts, and checksums are read-only.",
    "- [x] Source files, cited producer files and the release manifest are read-only.\n"
    "      Environments are recorded by `pip freeze` and conda export in the release\n"
    "      directory rather than locked.",
    "Checklist: read-only (user approved making sources read-only).", "promotion", REL)
add("R-source-tables", ROAD,
    "- [x] Figure source tables reproduce every plotted value exactly.",
    "- [x] Main-figure source tables reproduce every plotted value exactly (46 panels,\n"
    "      panel-deep review 2026-10-17); supplementary panels were not rechecked.",
    "Checklist wording.", "promotion", f"{DEEP}/REPORT.md")
add("R-tested-negative", ROAD,
    "testability, ATAC state counts, and zero strict Catalog `tested_negative` calls are presented",
    "testability, ATAC state counts, and Catalog `tested_negative` calls (not assessed) are presented",
    "Checklist wording; RESULTS records Catalog tested_negative calls as not assessed.", "promotion",
    "docs/RESULTS.md ('| Catalog `tested_negative` calls | not assessed |')")
add("R-numbers", ROAD,
    "- [ ] Current manuscript numbers and figure legends agree with",
    "- [x] Current manuscript numbers and figure legends agree with",
    "Checklist: text merge and open-item edits applied 2026-10-17.", "promotion",
    f"{DEEP}/text_merge-20261017T102723Z/README.md; {DEEP}/open_items-20261017T105826Z/README.md")
add("R-checks", ROAD,
    "- [ ] The scope validator and release-consistency checks pass.",
    "- [x] The scope validator and the release manifest check\n"
    "      (`scripts/manuscript/check_release_manifest.py`: checksums plus main-panel\n"
    "      source tables) pass; the July-only `check_release_consistency.py` is retired\n"
    "      for this release.",
    "Checklist: replaces the July release check (user approved).", "promotion",
    "scripts/manuscript/check_release_manifest.py")
add("R-approval", ROAD,
    "- [ ] The user explicitly approves promotion.",
    "- [x] The user explicitly approves promotion: approved by the user 2026-10-17.",
    "Checklist: user approval relayed by team-lead 2026-10-17.", "promotion", "team-lead message 2026-10-17")

# --------------------------------------------------------------------------
# STATUS
# --------------------------------------------------------------------------
add("S-top", STA,
    "# Current paper and codebase status",
    "# Current paper and codebase status\n\n"
    "**Release 2026-10-17-r2, promoted 2026-10-17 with user approval.** Its manifest is\n"
    f"`{REL}`; it supersedes 2026-10-17-r1 from the same day. Panel-deep review changes\n"
    f"(`{DEEP}/REPORT.md`):\n\n"
    "- Decision 1: Figures 1C and 2E draw genes never SuSiE-tested in the other trait class or in a non-EUR study as their own untested state (368 of the 414 enzyme-only genes; 213 of the 340 EUR-only genes).\n"
    "- Decision 2: Figure 7B draws graded fibrosis against plain expression with matched degrees of freedom (−0.128 macro Spearman), and the Nucleotide Transformer reporter-assay gain is drawn as including zero.\n"
    "- Decision 3: Figure 7E draws all seven severity arms, each labelled precommitted or post hoc, and draws signed eQTL direction as untestable.\n"
    "- Decision 4: Figures 6B to 6D label the genetic axis eQTL colocalization and mark the 31 genes with a fine-mapped coding variant at PIP ≥ 0.5; r8 classes are unchanged, and the 7 one-study-untestable genes stay untestable with a footnote.\n"
    "- Decision 5: Figure 5F prints the donor-level result beside the gene-level null (Stromal ECM, 3 of 4 donors positive, sign-flip P 0.19) and draws Vu et al. as source-dependent.\n"
    "- Decision 6: Figure 5C computes protein-histology associations within the 46 MASLD participants (60 of 125 tests pass BH).\n"
    "- Decision 7: Figure 4H plots the donor-first system estimates (rho 0.67, 34/43 same direction).\n"
    "- Decision 8: Figure 3C shows the mixed-sign share of multi-cohort genes with and without GSE126848 (39.7% and 13.5%).\n"
    "- Decision 9: Figure 5E counts protein members once per primary-assembly gene, as its gene column does (95 testable programs, 22 robust).\n"
    "- Figure 4D stability counts use the joint donor bootstrap, which keeps each donor's lineages together (76 stable matches, 34 stable mismatches, 7 unstable).\n"
    "- The ancestry companion of Figure 2D is adopted in the main figure, and the direct-subclass split moves to Figure S2S.",
    "Dated release paragraph at the top of STATUS.", "promotion",
    f"{REL}; {DEEP}/REPORT.md")
add("S-consistency-row", STA,
    "| The checker fails on the release it is meant to certify. The release is probably correct — drug "
    "calibration and the coding-architecture fraction were both deliberately demoted to supplement — so the "
    "required-file list is stale and should follow the release rules. |",
    "| Retired for release 2026-10-17-r2. Its required-file list was written for the July 10 release and fails "
    "on later releases. Release 2026-10-17-r2 is checked by `scripts/manuscript/check_release_manifest.py` "
    "(checksums plus main-panel source tables) and the scope validator. |",
    "Release-consistency row: retired for this release (user approved replacing the July check).", "promotion",
    "scripts/manuscript/check_release_manifest.py")

# --------------------------------------------------------------------------
# Order, check, write
# --------------------------------------------------------------------------
failures = []
for e in EDITS:
    e["_pos"] = live(e["file"]).find(e["old"])
    if live(e["file"]).count(e["old"]) != 1:
        failures.append(f"{e['id']}: old occurs {live(e['file']).count(e['old'])} times")
EDITS.sort(key=lambda e: (e["file"], e["_pos"]))
spans = defaultdict(list)
for e in EDITS:
    spans[e["file"]].append((e["_pos"], e["_pos"] + len(e["old"]), e["id"]))
for f, s in spans.items():
    s.sort()
    for (a0, a1, ai), (b0, b1, bi) in zip(s, s[1:]):
        if b0 < a1:
            failures.append(f"overlap in {f}: {ai} and {bi}")
for e in EDITS:
    del e["_pos"]
json.dump(EDITS, open(os.path.join(OUT, "edits.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)

with open(os.path.join(OUT, "file_ops.tsv"), "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["op", "source", "target", "sha256", "note"])
    for src, dst in FILE_OPS:
        sp, dp = os.path.join(ROOT, src), os.path.join(ROOT, dst)
        if not os.path.isfile(sp):
            failures.append(f"move source missing: {src}")
            continue
        if os.path.exists(dp):
            failures.append(f"move target exists: {dst}")
        w.writerow(["move", src, dst, hashlib.sha256(open(sp, "rb").read()).hexdigest(),
                    "item 8: direct-subclass companion to Figure S2S; byte-identical move"])

text = {}
if os.path.isdir(MIRROR):
    shutil.rmtree(MIRROR)  # scratch copy owned by this script
for f in {e["file"] for e in EDITS}:
    dst = os.path.join(MIRROR, f)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(os.path.join(ROOT, f), dst)
    text[f] = live(f)
applied = Counter()
for e in EDITS:
    if text[e["file"]].count(e["old"]) != 1:
        failures.append(f"{e['id']}: old not unique when applied")
        continue
    if len(BANNED.findall(e["new"])) > len(BANNED.findall(e["old"])):
        failures.append(f"{e['id']}: adds banned term or em dash {BANNED.findall(e['new'])}")
    text[e["file"]] = text[e["file"]].replace(e["old"], e["new"], 1)
    applied[e["file"]] += 1

norm = lambda s: " ".join(s.split())
markers = {
    ROAD: ["## Critical path", "## Final release checklist", "Plans 40–44", "Plans 45–46C",
           "No follow-up screen result", "user explicitly approves promotion",
           "every main figure must expose the claim"],
    STA: ["Validated update selected for the paper", "Corrected rerun PROMOTED to canonical 2026-08-17",
          "Five-cohort synchronized candidate", "117 program memberships", "Separate future paper"],
    RES: ["| Catalog `tested_negative` calls | not assessed |", "COLOC rerun is COMPLETE and PROMOTED to canonical"],
    LEG: ["machine-readable evidence state"],
    FIG2MD: ["next discriminating experiment"],
    IDX: ["Shared Figure 1–7 decision rules"],
}
for f, ms in markers.items():
    for m in ms:
        if norm(m) in norm(live(f)) and norm(m) not in norm(text[f]):
            failures.append(f"validator marker removed: {f}: {m}")

# The moved paragraph must appear exactly once in the result, under the S2 heading.
leg = text[LEG]
if leg.count(moved_body) != 1 or leg.find(moved_body) < leg.find("## Supplementary Figure S2."):
    failures.append("moved subclass paragraph not found exactly once under the S2 heading")
if "direct-subclass panel)**" in leg:
    failures.append("old subclass paragraph still in the Figure 2 legend")

diff = []
for f in sorted(text):
    with open(os.path.join(MIRROR, f), "w", encoding="utf-8", newline="") as fh:
        fh.write(text[f])
    diff += difflib.unified_diff(live(f).splitlines(True), text[f].splitlines(True),
                                 fromfile=f"live/{f}", tofile=f"mirror/{f}", n=0)
open(os.path.join(OUT, "mirror_diff.txt"), "w", encoding="utf-8").writelines(diff)

notes = []
for p in (REL, "scripts/manuscript/check_release_manifest.py"):
    if not os.path.exists(os.path.join(ROOT, p)):
        notes.append(f"cited but not on disk yet: {p}")
print(f"applied {sum(applied.values())} edits")
for f in sorted(applied):
    print(f"  {applied[f]:3d}  {f}")
print(f"file moves: {len(FILE_OPS)}")
for n in notes:
    print("NOTE", n)
print("failures:", failures if failures else "none")
assert not failures

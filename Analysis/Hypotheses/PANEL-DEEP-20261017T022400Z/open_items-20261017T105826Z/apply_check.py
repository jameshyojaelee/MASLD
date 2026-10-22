"""Merge part_A/B/C.json into edits.json and resolution tsvs into resolution.tsv,
then copy the current live files to mirror/, apply edits.json in order, and check:
each old occurs once when applied, no two edits overlap in the live text, no edit
adds an em dash or a retired term, and the marker strings required by
scripts/manuscript/validate_resource_scope.py survive. Writes mirror_diff.txt.
Live files are only read.
"""
import csv
import difflib
import json
import os
import re
import shutil
from collections import Counter, defaultdict

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/open_items-20261017T105826Z")
MIRROR = os.path.join(OUT, "mirror")
BANNED = re.compile(r"—|\bfrozen\b|\bprogression\b|\bcascade\b", re.IGNORECASE)

edits = []
for part in "ABC":
    edits += json.load(open(os.path.join(OUT, f"part_{part}.json"), encoding="utf-8"))

# Wording adjustments made at merge time.
for e in edits:
    if e["id"] == "A-03":
        e["new"] = e["new"].replace("through two adopted consumers that read the same files",
                                    "through two adopted analyses that read the same files")
    if e["id"] == "B-6":
        e["decision"] = "item 6 (USER; applies the recommended joint scheme; drop this edit if declined)"

live = {}
for e in edits:
    live.setdefault(e["file"], open(os.path.join(ROOT, e["file"]), encoding="utf-8").read())
    e["_pos"] = live[e["file"]].find(e["old"])
failures = []
for e in edits:
    if live[e["file"]].count(e["old"]) != 1:
        failures.append(f"{e['id']}: old occurs {live[e['file']].count(e['old'])} times in live {e['file']}")
edits.sort(key=lambda e: (e["file"], e["_pos"]))
spans = defaultdict(list)
for e in edits:
    spans[e["file"]].append((e["_pos"], e["_pos"] + len(e["old"]), e["id"]))
for f, s in spans.items():
    s.sort()
    for (a0, a1, ai), (b0, b1, bi) in zip(s, s[1:]):
        if b0 < a1:
            failures.append(f"overlap in {f}: {ai} and {bi}")
for e in edits:
    del e["_pos"]
json.dump(edits, open(os.path.join(OUT, "edits.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)

# resolution.tsv
rows = []
for part in "ABC":
    with open(os.path.join(OUT, f"part_{part}_resolution.tsv"), encoding="utf-8") as fh:
        rows += list(csv.DictReader(fh, delimiter="\t"))
rows.sort(key=lambda r: int(r["item"]))
with open(os.path.join(OUT, "resolution.tsv"), "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["item", "status", "action", "producer_path", "note"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)
expected_items = {2, 3, 4, 6, 8, 9, 11, 12, 13, 18, 20, 21, 22, 23, 24, 25, 26, 28, 29, 30, 31, 32, 33, 34}
got_items = [int(r["item"]) for r in rows]
if set(got_items) != expected_items or len(got_items) != len(expected_items):
    failures.append(f"resolution items differ: missing {sorted(expected_items - set(got_items))}, "
                    f"extra {sorted(set(got_items) - expected_items)}")

# Apply to a mirror.
if os.path.isdir(MIRROR):
    shutil.rmtree(MIRROR)  # scratch copy owned by this script
text = {}
for f in live:
    dst = os.path.join(MIRROR, f)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(os.path.join(ROOT, f), dst)
    text[f] = live[f]
applied = Counter()
for e in edits:
    if text[e["file"]].count(e["old"]) != 1:
        failures.append(f"{e['id']}: old not unique when applied")
        continue
    if len(BANNED.findall(e["new"])) > len(BANNED.findall(e["old"])):
        failures.append(f"{e['id']}: adds banned term or em dash {BANNED.findall(e['new'])}")
    text[e["file"]] = text[e["file"]].replace(e["old"], e["new"], 1)
    applied[e["file"]] += 1
assert sum(applied.values()) == len(edits), (applied, len(edits))

norm = lambda s: " ".join(s.split())
contract = json.load(open(os.path.join(ROOT, "config/resource_paper_scope.config.json")))
markers = {
    "docs/manuscript/05_figure_legends.md": ["machine-readable evidence state"],
    "docs/manuscript/draft/fig2.md": ["next discriminating experiment"],
    "docs/manuscript/draft/fig4.md": ["next discriminating experiment"],
    "docs/manuscript/draft/fig5_discussion.md": ["deterministic next-experiment rule"],
}
for f, ms in markers.items():
    if f in text:
        for m in ms:
            if norm(m) in norm(live[f]) and norm(m) not in norm(text[f]):
                failures.append(f"validator marker removed: {f}: {m}")

diff = []
for f in sorted(text):
    with open(os.path.join(MIRROR, f), "w", encoding="utf-8") as fh:
        fh.write(text[f])
    diff += difflib.unified_diff(live[f].splitlines(True), text[f].splitlines(True),
                                 fromfile=f"live/{f}", tofile=f"mirror/{f}", n=0)
open(os.path.join(OUT, "mirror_diff.txt"), "w", encoding="utf-8").writelines(diff)

print(f"applied {sum(applied.values())} edits")
for f in sorted(applied):
    print(f"  {applied[f]:3d}  {f}")
print("status:", dict(Counter(r["status"] for r in rows)))
print("failures:", failures if failures else "none")
assert not failures

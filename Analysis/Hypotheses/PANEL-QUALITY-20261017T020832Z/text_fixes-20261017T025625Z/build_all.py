"""Merge the staged edit fragments in work/ into edits.json and check.tsv.

Checks, against the live files at run time:
  - each `old` occurs exactly once in its live file;
  - applying the edits of one file in order (in memory only), each `old` is
    still unique at its turn, so the list can be applied by a later script;
  - no `new` string adds an em dash or the words frozen/progression/cascade.
No live file is written."""
import csv
import json
import os
import re
from collections import defaultdict

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-QUALITY-20261017T020832Z/text_fixes-20261017T025625Z")
WORK = os.path.join(OUT, "work")
os.chdir(ROOT)

# N-02 cites a recomputation, not a producing table; it is routed to needs_user.md.
DROP = {"N-02"}
FRAGMENTS = ["lead_edits.json", "legends_edits.json", "index_edits.json", "nosource_edits.json"]
REREADS = ["lead_reread.tsv", "legends_reread.tsv", "index_reread.tsv"]
BANNED = re.compile(r"—|\bfrozen\b|\bprogression\b|\bcascade\b", re.IGNORECASE)
NUM = re.compile(r"(?<![\w.])[-−+]?\d[\d,]*(?:\.\d+)?%?")

edits = []
for frag in FRAGMENTS:
    for e in json.load(open(os.path.join(WORK, frag))):
        if e["id"] not in DROP:
            edits.append(e)

reread = defaultdict(list)
for fn in REREADS:
    with open(os.path.join(WORK, fn)) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            reread[row["edit_id"]].append(f"{row['column_or_row']}={row['value_reread']}")

live = {}
for e in edits:
    if e["file"] not in live:
        live[e["file"]] = open(e["file"], encoding="utf-8").read()

problems = []
rows = []
state = dict(live)
for e in edits:
    n_live = live[e["file"]].count(e["old"])
    n_turn = state[e["file"]].count(e["old"])
    if n_live != 1:
        problems.append(f"{e['id']}: old occurs {n_live} times in live {e['file']}")
    if n_turn != 1:
        problems.append(f"{e['id']}: old occurs {n_turn} times after earlier edits in {e['file']}")
    else:
        state[e["file"]] = state[e["file"]].replace(e["old"], e["new"], 1)
    added = BANNED.findall(e["new"])
    if len(added) > len(BANNED.findall(e["old"])):
        problems.append(f"{e['id']}: new adds banned term(s) {added}")
    old_nums = NUM.findall(e["old"])
    new_nums = [x for x in NUM.findall(e["new"]) if x not in old_nums]
    vr = reread.get(e["id"]) or [str(e.get("value_reread", ""))]
    if e["id"].startswith("T-09"):
        vr = reread.get("T-09", vr)
    if e["id"].startswith("T-28"):
        vr = reread.get("T-28", vr)
    rows.append([e["id"], e["file"], n_live, "; ".join(new_nums) if new_nums else "none",
                 " | ".join(vr), e["source"]])

with open(os.path.join(OUT, "edits.json"), "w", encoding="utf-8") as f:
    json.dump([{k: e[k] for k in ["id", "file", "old", "new", "reason", "source"]} for e in edits],
              f, indent=1, ensure_ascii=False)
with open(os.path.join(OUT, "check.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["id", "file", "old_count_in_live", "numbers_in_new", "value_reread", "source_path"])
    w.writerows(rows)

by_file = defaultdict(int)
for e in edits:
    by_file[e["file"]] += 1
print(f"{len(edits)} edits")
for k, v in sorted(by_file.items()):
    print(f"  {v:3d}  {k}")
print("problems:", problems if problems else "none")

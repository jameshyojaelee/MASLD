import json, re, sys
from pathlib import Path
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
S = json.load(open(ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F2d/out/structured_output.json"))
TAGS = re.compile(r"\[(Inference|Unverified|Speculation|Estimate|Contains Unverified Content|Observed pattern[^\]]*)\]")
BAD = ["—", "progression", "causal gene", "significant", "Mendelian", "mega"]
out = []
for e in S["entries"]:
    for i, t in enumerate(e["text_edits"]):
        f = ROOT / t["file"]
        txt = f.read_text()
        n_old = txt.count(t["old_text"])
        n_new = txt.count(t["new_text"])
        # line of old_text
        line = txt[: txt.find(t["old_text"])].count("\n") + 1 if n_old else None
        flags = [b for b in BAD if b in t["new_text"]]
        tags = TAGS.findall(t["new_text"])
        out.append(dict(fid=e["fid"], i=i, file=t["file"], line_given=t["line"], line_found=line,
                        n_old=n_old, n_new_present=n_new, bad=flags, tags=tags))
for o in out:
    print(o)

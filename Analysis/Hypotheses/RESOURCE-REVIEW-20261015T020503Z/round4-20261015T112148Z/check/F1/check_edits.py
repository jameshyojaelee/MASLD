"""Check each proposed old_text exists verbatim once, at the stated line range; scan new_text
for em dashes and simple style flags. Read-only on docs/."""
import re
import sys
sys.path.insert(0, "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F1")
from edits_input import EDITS

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
for i, (f, line, old, new) in enumerate(EDITS):
    txt = open(ROOT + f, encoding="utf-8").read()
    n = txt.count(old)
    if n:
        pos = txt.find(old)
        l0 = txt[:pos].count("\n") + 1
        l1 = l0 + old.rstrip("\n").count("\n")
        # what follows the old_text on the same line (join point)
        after = txt[pos + len(old): pos + len(old) + 60].split("\n")[0]
        before = txt[max(0, pos - 60):pos].split("\n")[-1]
    else:
        l0 = l1 = None; after = before = ""
    em = [c for c in new if c in "—"]
    print(f"[{i}] {f}:{line} found={n} actual_lines={l0}-{l1} em_dash_in_new={len(em)} "
          f"prefix_on_line={before!r} suffix_on_line={after!r}")
    if n == 0:
        # find longest prefix that matches to locate the mismatch
        lo, hi = 0, len(old)
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if old[:mid] in txt:
                lo = mid
            else:
                hi = mid - 1
        print(f"     longest matching prefix {lo}/{len(old)} chars; mismatch near: {old[max(0,lo-30):lo+30]!r}")

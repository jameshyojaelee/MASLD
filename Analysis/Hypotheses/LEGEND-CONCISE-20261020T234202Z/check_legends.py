#!/usr/bin/env python3
"""For each staged legend section: words before/after, numbers new to the staged text (must be none), banned words, panel letters lost."""
import re, sys
from pathlib import Path
RUN = Path(__file__).resolve().parent
NUM = re.compile(r"(?<![A-Za-z0-9])[−-]?\d[\d,]*(?:\.\d+)?(?:e-?\d+)?")
BANNED = re.compile(r"—|\bnotably\b|\bimportantly\b|\bremarkably\b|\bprogression\b|causal gene", re.I)
PANEL = re.compile(r"\*\*\(([A-Z]{1,2})\)\*\*")
def nums(t): return {n.replace(",", "").replace("−", "-").lstrip("-") for n in NUM.findall(t)}
bad = False
for base in sorted((RUN / "baseline").glob("Figure_*.md")):
    st = RUN / "staged" / base.name
    if not st.exists(): print(f"{base.stem}\tNOT_STAGED"); bad = True; continue
    b, s = base.read_text(), st.read_text()
    new = sorted(nums(s) - nums(b))
    lost_panels = sorted(set(PANEL.findall(b)) - set(PANEL.findall(s)))
    banned = sorted({m.group(0) for m in BANNED.finditer(s)})
    print(f"{base.stem}\twords {len(b.split())}->{len(s.split())}\tnew_numbers {new}\tlost_panels {lost_panels}\tbanned {banned}")
    if new or lost_panels or banned: bad = True
sys.exit(1 if bad else 0)

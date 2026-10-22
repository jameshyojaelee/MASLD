"""Figure 5A check: box text content vs the legend/analysis facts, text-in-box containment, arrow endpoints, text colours.
Usage: python geom_5a.py <pdf> <tag>
"""
import sys
import pdfplumber

pdf, tag = sys.argv[1], sys.argv[2]
p = pdfplumber.open(pdf).pages[0]

def hexc(c):
    if c is None: return ""
    if isinstance(c, (int, float)): c = (c,)
    if len(c) == 1: c = (c[0],) * 3
    return "#%02X%02X%02X" % tuple(int(round(v * 255)) for v in c)

boxes = [r for r in p.rects if r["width"] > 40]
words = p.extract_words(extra_attrs=["size", "non_stroking_color"])
print(tag, "boxes", len(boxes), "words", len(words))
issues = []
for b in boxes:
    half = b["linewidth"] / 2
    inner = (b["x0"] + half, b["x1"] - half, b["top"] + half, b["bottom"] - half)
    inside = [w for w in words if (w["top"] + w["bottom"]) / 2 > b["top"] - 2 and (w["top"] + w["bottom"]) / 2 < b["bottom"] + 2
              and w["x1"] > b["x0"] - 6 and w["x0"] < b["x1"] + 6 and not (w["x1"] < b["x0"] - 6)]
    for w in inside:
        ov = {"left": inner[0] - w["x0"], "right": w["x1"] - inner[1], "top": inner[2] - w["top"], "bottom": w["bottom"] - inner[3]}
        bad = {k: round(v, 2) for k, v in ov.items() if v > 0}
        if bad:
            issues.append((round(b["x0"], 1), round(b["top"], 1), w["text"], bad))
print(tag, "words crossing a box's inner edge (pt beyond edge):")
for i in issues:
    print("  ", i)
cols = {}
for w in words:
    cols.setdefault(hexc(w["non_stroking_color"]), []).append(w["text"])
print(tag, "text colours", {k: len(v) for k, v in cols.items()}, {k: " ".join(v)[:80] for k, v in cols.items()})
print(tag, "font sizes", sorted(set(round(w["size"], 2) for w in words)), "fonts", sorted(set(c["fontname"] for c in p.chars)))

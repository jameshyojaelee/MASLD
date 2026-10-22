"""True font size per char from the text matrix (handles rotated text), font names, and text fill colours, per panel PDF.
Usage: python fonts_check.py <pdf> [<pdf> ...]"""
import sys, pdfplumber, numpy as np, collections
from pdfmarks import col
for p in sys.argv[1:]:
    pg = pdfplumber.open(p).pages[0]
    sizes, fonts, fills = collections.Counter(), collections.Counter(), collections.Counter()
    for c in pg.chars:
        a, b, cc, d, e, f = c["matrix"]
        # glyph vertical scale = |(c,d)|; pdfminer matrix already includes font size
        sizes[round(float(np.hypot(cc, d)), 2)] += 1
        fonts[c["fontname"]] += 1
        fills[col(c.get("non_stroking_color"))] += 1
    print(p.split("/")[-1], "| n chars", len(pg.chars), "| sizes", dict(sizes), "| fonts", dict(fonts), "| text fills", dict(fills))

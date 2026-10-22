"""Dump every vector mark in a panel PDF (chars, rects, lines, curves) with colours.

Usage: python extract_geometry.py <pdf> <out.json>
"""
import json
import sys

import pdfplumber


def col(c):
    if c is None:
        return None
    if isinstance(c, (list, tuple)):
        vals = [float(v) for v in c]
        if len(vals) == 1:
            vals = vals * 3
        if len(vals) == 4:  # CMYK -> RGB approx
            C, M, Y, K = vals
            vals = [(1 - C) * (1 - K), (1 - M) * (1 - K), (1 - Y) * (1 - K)]
        return "#%02X%02X%02X" % tuple(int(round(v * 255)) for v in vals[:3])
    return str(c)


def main(pdf, out):
    res = {}
    with pdfplumber.open(pdf) as doc:
        page = doc.pages[0]
        res["width"] = float(page.width)
        res["height"] = float(page.height)
        res["chars"] = [
            {"text": c["text"], "x0": c["x0"], "x1": c["x1"], "top": c["top"], "bottom": c["bottom"],
             "size": c["size"], "font": c["fontname"], "fill": col(c.get("non_stroking_color"))}
            for c in page.chars
        ]
        res["words"] = [
            {"text": w["text"], "x0": w["x0"], "x1": w["x1"], "top": w["top"], "bottom": w["bottom"]}
            for w in page.extract_words(keep_blank_chars=False, use_text_flow=False, x_tolerance=1.5)
        ]
        for kind in ("rects", "lines", "curves"):
            items = []
            for o in getattr(page, kind):
                d = {"x0": o["x0"], "x1": o["x1"], "top": o["top"], "bottom": o["bottom"],
                     "lw": o.get("linewidth"), "stroke": col(o.get("stroking_color")),
                     "fill": col(o.get("non_stroking_color")), "filled": o.get("fill"),
                     "stroked": o.get("stroke")}
                if kind != "rects" and "pts" in o:
                    d["pts"] = [(float(p[0]), float(p[1])) for p in o["pts"]]
                dash = o.get("dash")
                if dash:
                    d["dash"] = str(dash)
                items.append(d)
            res[kind] = items
    with open(out, "w") as fh:
        json.dump(res, fh, indent=0, default=float)
    print(pdf, {k: len(v) for k, v in res.items() if isinstance(v, list)})


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

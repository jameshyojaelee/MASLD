"""Dump vector marks from a one-page PDF panel with pdfplumber.

Usage: python dump_geometry.py <pdf> <out_prefix>
Writes <out_prefix>_chars.tsv, _rects.tsv, _lines.tsv, _curves.tsv, _words.tsv.
Coordinates are pdfplumber top-left (x0, top); y_pdf = page height - top.
"""
import sys
import pandas as pd
import pdfplumber


def col(c):
    if c is None:
        return ""
    if isinstance(c, (list, tuple)):
        if len(c) == 3:
            return "#%02X%02X%02X" % tuple(int(round(v * 255)) for v in c)
        if len(c) == 1:
            v = int(round(c[0] * 255))
            return "#%02X%02X%02X" % (v, v, v)
        if len(c) == 4:  # CMYK
            C, M, Y, K = c
            return "#%02X%02X%02X" % tuple(int(round(255 * (1 - a) * (1 - K))) for a in (C, M, Y))
    return str(c)


def main(pdf, prefix):
    with pdfplumber.open(pdf) as doc:
        p = doc.pages[0]
        H = p.height
        chars = pd.DataFrame([{
            "text": c["text"], "x0": c["x0"], "x1": c["x1"], "top": c["top"], "bottom": c["bottom"],
            "size": round(c["size"], 2), "font": c["fontname"], "color": col(c.get("non_stroking_color")),
            "upright": c.get("upright"),
        } for c in p.chars])
        chars.to_csv(prefix + "_chars.tsv", sep="\t", index=False)
        words = p.extract_words(keep_blank_chars=False, use_text_flow=False, extra_attrs=["size"])
        pd.DataFrame(words).to_csv(prefix + "_words.tsv", sep="\t", index=False)
        rows = []
        for r in p.rects:
            rows.append({"x0": r["x0"], "x1": r["x1"], "top": r["top"], "bottom": r["bottom"],
                         "w": r["width"], "h": r["height"], "fill": col(r.get("non_stroking_color")) if r.get("fill") else "",
                         "stroke": col(r.get("stroking_color")) if r.get("stroke") else "", "lw": r.get("linewidth")})
        pd.DataFrame(rows).to_csv(prefix + "_rects.tsv", sep="\t", index=False)
        rows = []
        for l in p.lines:
            rows.append({"x0": l["x0"], "x1": l["x1"], "top": l["top"], "bottom": l["bottom"],
                         "stroke": col(l.get("stroking_color")), "lw": l.get("linewidth"),
                         "pts": ";".join("%.3f,%.3f" % (x, y) for x, y in l["pts"])})
        pd.DataFrame(rows).to_csv(prefix + "_lines.tsv", sep="\t", index=False)
        rows = []
        for c in p.curves:
            pts = c["pts"]
            rows.append({"x0": c["x0"], "x1": c["x1"], "top": c["top"], "bottom": c["bottom"],
                         "cx": (c["x0"] + c["x1"]) / 2, "cy": (c["top"] + c["bottom"]) / 2,
                         "w": c["x1"] - c["x0"], "h": c["bottom"] - c["top"],
                         "fill": col(c.get("non_stroking_color")) if c.get("fill") else "",
                         "stroke": col(c.get("stroking_color")) if c.get("stroke") else "", "lw": c.get("linewidth"),
                         "npts": len(pts)})
        pd.DataFrame(rows).to_csv(prefix + "_curves.tsv", sep="\t", index=False)
        print(pdf, "W", p.width, "H", H, "chars", len(p.chars), "rects", len(p.rects), "lines", len(p.lines),
              "curves", len(p.curves), "images", len(p.images))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

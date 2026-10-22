"""Dump vector marks (chars, words, rects, lines, curves) of a one-page PDF to TSVs with colours.

usage: pdf_dump.py <pdf> <outprefix>
Coordinates are pdfplumber top-left based (x0, top, x1, bottom) in pt; also y (bottom-left) = page_height - top.
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
            v = int(round(c[0] * 255)); return "#%02X%02X%02X" % (v, v, v)
        if len(c) == 4:  # cmyk
            C, M, Y, K = c
            return "#%02X%02X%02X" % tuple(int(round(255 * (1 - x) * (1 - K))) for x in (C, M, Y))
    return str(c)


def dump(pdf, prefix):
    with pdfplumber.open(pdf) as p:
        pg = p.pages[0]
        H = pg.height
        meta = {"width": pg.width, "height": pg.height}
        ch = pd.DataFrame([{"text": c["text"], "x0": c["x0"], "x1": c["x1"], "top": c["top"], "bottom": c["bottom"],
                            "size": round(c["size"], 2), "font": c["fontname"], "fill": col(c.get("non_stroking_color")),
                            "upright": c.get("upright", True)} for c in pg.chars])
        words = pg.extract_words(extra_attrs=["size", "fontname"], keep_blank_chars=False, use_text_flow=False)
        wd = pd.DataFrame(words)
        out = {"chars": ch, "words": wd}
        for kind in ["rects", "lines", "curves"]:
            objs = getattr(pg, kind)
            rows = []
            for o in objs:
                pts = o.get("pts") or []
                rows.append({"x0": o["x0"], "x1": o["x1"], "top": o["top"], "bottom": o["bottom"],
                             "w": o["x1"] - o["x0"], "h": o["bottom"] - o["top"],
                             "fill": col(o.get("non_stroking_color")) if o.get("fill") else "",
                             "stroke": col(o.get("stroking_color")) if o.get("stroke") else "",
                             "lw": o.get("linewidth"), "npts": len(pts),
                             "pts": ";".join(f"{x:.3f},{y:.3f}" for x, y in pts[:12])})
            out[kind] = pd.DataFrame(rows)
        for k, v in out.items():
            v.to_csv(f"{prefix}.{k}.tsv", sep="\t", index=False)
        print(pdf, meta, {k: len(v) for k, v in out.items()})
        return meta, out


if __name__ == "__main__":
    dump(sys.argv[1], sys.argv[2])

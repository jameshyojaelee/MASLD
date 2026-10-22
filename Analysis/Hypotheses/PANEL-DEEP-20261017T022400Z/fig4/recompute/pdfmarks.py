"""Shared helpers: extract vector marks from a one-page PDF with pdfplumber."""
import pdfplumber, numpy as np, pandas as pd


def col(c):
    if c is None:
        return None
    if isinstance(c, (list, tuple)):
        if len(c) == 1:
            v = c[0]
            return "#%02X%02X%02X" % (round(v * 255), round(v * 255), round(v * 255))
        if len(c) == 3:
            return "#%02X%02X%02X" % tuple(round(v * 255) for v in c)
        if len(c) == 4:  # cmyk
            C, M, Y, K = c
            return "#%02X%02X%02X" % tuple(round(255 * (1 - x) * (1 - K)) for x in (C, M, Y))
    return str(c)


def load(path):
    pdf = pdfplumber.open(path)
    pg = pdf.pages[0]
    H = pg.height
    words = pg.extract_words(keep_blank_chars=False, use_text_flow=False, extra_attrs=["size", "fontname"])
    chars = pd.DataFrame([{"text": c["text"], "x0": c["x0"], "x1": c["x1"], "top": c["top"], "bottom": c["bottom"],
                           "size": c["size"], "font": c["fontname"], "fill": col(c.get("non_stroking_color")),
                           "upright": c.get("upright", True)} for c in pg.chars])
    w = pd.DataFrame([{"text": x["text"], "x0": x["x0"], "x1": x["x1"], "top": x["top"], "bottom": x["bottom"],
                       "xc": (x["x0"] + x["x1"]) / 2, "yc": (x["top"] + x["bottom"]) / 2, "size": x["size"]}
                      for x in words])

    def obj(o, kind):
        d = {"kind": kind, "x0": o["x0"], "x1": o["x1"], "top": o["top"], "bottom": o["bottom"],
             "w": o["x1"] - o["x0"], "h": o["bottom"] - o["top"],
             "stroke": col(o.get("stroking_color")), "fill": col(o.get("non_stroking_color")),
             "lw": o.get("linewidth"), "filled": o.get("fill"), "stroked": o.get("stroke")}
        d["xc"] = (d["x0"] + d["x1"]) / 2
        d["yc"] = (d["top"] + d["bottom"]) / 2
        if "pts" in o and o["pts"]:
            d["pts"] = [(p[0], p[1]) for p in o["pts"]]
        return d
    rects = pd.DataFrame([obj(o, "rect") for o in pg.rects])
    lines = pd.DataFrame([obj(o, "line") for o in pg.lines])
    curves = pd.DataFrame([obj(o, "curve") for o in pg.curves])
    return pdf, pg, chars, w, rects, lines, curves


def axis_fit(pix, vals):
    """least-squares linear map pixel->data from >=2 tick (pixel, value) pairs; returns (a,b,maxresid)"""
    pix = np.asarray(pix, float); vals = np.asarray(vals, float)
    A = np.vstack([pix, np.ones_like(pix)]).T
    (a, b), *_ = np.linalg.lstsq(A, vals, rcond=None)
    resid = np.abs(a * pix + b - vals).max() if len(pix) > 2 else 0.0
    return a, b, resid


def hex2rgb(h):
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], float) / 255.0

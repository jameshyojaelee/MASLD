"""Geometry and printed-label check for Figure 1B and 1C PDFs against recomputed values.

Usage: python geometry_check.py <out_prefix> <pdf> [<pdf> ...]
Each PDF is classified as 1B or 1C from its text, axis ticks are read from tick marks
and tick labels, a linear transform is fitted, every data rect and dashed line is
inverted to a data value, and each is compared to the expected value recomputed in
recompute_1b.json / recompute_1c.json.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pdfplumber

HERE = Path(__file__).resolve().parent
R1B = json.load(open(HERE / "recompute_1b.json"))
R1C = json.load(open(HERE / "recompute_1c.json"))

WHITE = (1.0, 1.0, 1.0)


def hexcol(c):
    if c is None:
        return "none"
    if isinstance(c, (int, float)):
        c = (c, c, c)
    if len(c) == 1:
        c = (c[0],) * 3
    if len(c) == 4:  # cmyk
        C, M, Y, K = c
        c = ((1 - C) * (1 - K), (1 - M) * (1 - K), (1 - Y) * (1 - K))
    return "#%02X%02X%02X" % tuple(int(round(255 * v)) for v in c[:3])


def words(p):
    return p.extract_words(keep_blank_chars=True, extra_attrs=["fontname", "size", "non_stroking_color"])


def fit_axis(p, ws, kind):
    """Return (a, b, ticks) with x_pt = a + b * value, from tick marks + labels."""
    if kind == "1B":
        labs = [w for w in ws if re.fullmatch(r"\d+(\.\d+)?%", w["text"].strip())]
        vals = [float(w["text"].strip().rstrip("%")) for w in labs]
    else:
        labs = [w for w in ws if re.fullmatch(r"\d+", w["text"].strip())]
        # tick labels: the lowest row of numeric labels
        if not labs:
            return None
        ymax = max(w["top"] for w in labs)
        labs = [w for w in labs if abs(w["top"] - ymax) < 1.0]
        vals = [float(w["text"]) for w in labs]
    # vertical tick lines: short vertical segments
    ticks = [l for l in p.lines if abs(l["x0"] - l["x1"]) < 0.01 and 0.5 < (l["bottom"] - l["top"]) < 4]
    pairs = []
    for w, v in zip(labs, vals):
        xc = (w["x0"] + w["x1"]) / 2
        if ticks:
            t = min(ticks, key=lambda l: abs(l["x0"] - xc))
            if abs(t["x0"] - xc) < 4 and t["top"] < w["top"]:
                pairs.append((v, t["x0"], xc))
                continue
        pairs.append((v, None, xc))
    use = [(v, x) for v, x, _ in pairs if x is not None]
    if len(use) < 2:
        use = [(v, xc) for v, _, xc in pairs]
        src = "label_centres"
    else:
        src = "tick_marks"
    v = np.array([u[0] for u in use]); x = np.array([u[1] for u in use])
    b, a = np.polyfit(v, x, 1)
    resid = x - (a + b * v)
    return dict(a=a, b=b, src=src, n_ticks=len(use), max_resid_pt=float(np.abs(resid).max()),
                tick_values=v.tolist(), tick_x=x.tolist(),
                label_offsets=[float(xc - x0) for _, x0, xc in pairs if x0 is not None])


def data_rects(p, plot_top=None):
    out = []
    for r in p.rects:
        col = r.get("non_stroking_color")
        if not r.get("fill") or col is None:
            continue
        if tuple(round(c, 3) for c in (col if len(col) == 3 else (col[0],) * 3)) == WHITE:
            continue
        h = r["bottom"] - r["top"]
        out.append(dict(x0=r["x0"], x1=r["x1"], top=r["top"], bottom=r["bottom"], h=h, col=hexcol(col)))
    return out


def ylab_for(ws, yc, xmax):
    cands = [w for w in ws if w["x1"] <= xmax + 0.5 and abs((w["top"] + w["bottom"]) / 2 - yc) < 8]
    return " ".join(w["text"] for w in sorted(cands, key=lambda w: w["top"]))


def check_1b(path, p, ws, rows):
    ax = fit_axis(p, ws, "1B")
    inv = lambda x: (x - ax["a"]) / ax["b"]
    rng = 100.0
    tol_rng = 0.005 * rng
    rects = [r for r in data_rects(p) if r["h"] > 12]
    legend = [r for r in data_rects(p) if r["h"] <= 12]
    G, D, O = R1B["G"], R1B["D"], R1B["O"]
    E = R1B["strat_10bins"]["E"]
    exp = {
        "coloc": dict(both=100 * O / G, total=100.0, dashed=100 * E / G),
        "deg": dict(both=100 * O / D, total=100.0, dashed=100 * E / D),
    }
    # group bars by row
    rowsy = {}
    for r in rects:
        key = round((r["top"] + r["bottom"]) / 2, 1)
        rowsy.setdefault(key, []).append(r)
    axis_x0 = ax["a"]
    ycentres = sorted(rowsy)
    out_rows = []
    for yc in ycentres:
        segs = sorted(rowsy[yc], key=lambda r: r["x0"])
        # de-duplicate stroke-only overlays already filtered; merge identical extents
        uniq = []
        for s in segs:
            if not any(abs(s["x0"] - u["x0"]) < 0.01 and abs(s["x1"] - u["x1"]) < 0.01 for u in uniq):
                uniq.append(s)
        lab = ylab_for(ws, yc, axis_x0)
        which = "coloc" if re.search(r"olocali", lab) else ("deg" if re.search(r"DEG|disease", lab, re.I) else "unknown")
        out_rows.append((yc, which, lab, uniq))
    dashed = [l for l in p.lines if l.get("dash") and l["dash"][0]]
    for yc, which, lab, segs in out_rows:
        e = exp.get(which)
        for i, s in enumerate(segs):
            v0, v1 = inv(s["x0"]), inv(s["x1"])
            if i == 0:
                ev0, ev1, name = 0.0, (e["both"] if e else np.nan), "both_segment"
            else:
                ev0, ev1, name = (e["both"] if e else np.nan), 100.0, "only_segment"
            for nm, dv, evv in ((name + "_start", v0, ev0), (name + "_end", v1, ev1)):
                d = dv - evv
                rows.append(dict(file=path, panel="1B", row=which, ylabel=lab, mark=nm, colour=s["col"], drawn=round(dv, 4),
                                 expected=round(evv, 4), diff=round(d, 4), tol=tol_rng, off=bool(abs(d) > tol_rng)))
        # dashed line overlapping this row
        dl = [l for l in dashed if l["top"] <= yc <= l["bottom"]]
        if not dl:
            rows.append(dict(file=path, panel="1B", row=which, ylabel=lab, mark="expected_dashed", colour="", drawn=None,
                             expected=round(e["dashed"], 4) if e else None, diff=None, tol=tol_rng, off=True))
        for l in dl:
            dv = inv(l["x0"]); evv = e["dashed"] if e else np.nan
            rows.append(dict(file=path, panel="1B", row=which, ylabel=lab, mark="expected_dashed", colour=hexcol(l.get("stroking_color")),
                             drawn=round(dv, 4), expected=round(evv, 4), diff=round(dv - evv, 4), tol=tol_rng,
                             off=bool(abs(dv - evv) > tol_rng)))
    return ax, out_rows, legend


def check_1c(path, p, ws, rows):
    ax = fit_axis(p, ws, "1C")
    inv = lambda x: (x - ax["a"]) / ax["b"]
    plot_right = max(l["x1"] for l in p.lines)
    rng = inv(plot_right) - inv(ax["a"])
    tol_rng = 0.005 * rng
    cls = R1C["class_counts_positive436"]
    expmap = {"enzyme": cls["enzyme_only"], "both": cls["both"], "direct": cls["direct_only"]}
    bars = [r for r in data_rects(p) if r["h"] > 8]
    out_rows = []
    for r in bars:
        yc = (r["top"] + r["bottom"]) / 2
        lab = ylab_for(ws, yc, ax["a"])
        if re.search(r"enzyme", lab, re.I) and not re.search(r"both", lab, re.I):
            k = "enzyme"
        elif re.search(r"^\s*both", lab, re.I):
            k = "both"
        elif re.search(r"direct|PDFF|NAFLD", lab, re.I):
            k = "direct"
        else:
            k = "unknown"
        dv0, dv1 = inv(r["x0"]), inv(r["x1"])
        evv = expmap.get(k, np.nan)
        tol = max(0.5, tol_rng)
        rows.append(dict(file=path, panel="1C", row=k, ylabel=lab, mark="bar_start", colour=r["col"], drawn=round(dv0, 3), expected=0.0,
                         diff=round(dv0, 3), tol=round(tol, 3), off=bool(abs(dv0) > tol)))
        rows.append(dict(file=path, panel="1C", row=k, ylabel=lab, mark="bar_end", colour=r["col"], drawn=round(dv1, 3), expected=evv,
                         diff=round(dv1 - evv, 3), tol=round(tol, 3), off=bool(abs(dv1 - evv) > tol),
                         strict_off=bool(abs(dv1 - evv) > 0.5)))
        out_rows.append((yc, k, lab, r))
    return ax, out_rows, rng


def label_checks_1b(path, ws, labrows):
    G, D, O = R1B["G"], R1B["D"], R1B["O"]
    E = R1B["strat_10bins"]["E"]
    txt = "\n".join(w["text"] for w in ws)
    expected = {
        "gen_only_count": (str(G - O), r"\b%d\b" % (G - O)),
        "gen_only_pct": ("%.1f%%" % (100 * (G - O) / G), re.escape("%.1f%%" % (100 * (G - O) / G))),
        "deg_only_count": (str(D - O), r"\b%d\b" % (D - O)),
        "deg_only_pct": ("%.1f%%" % (100 * (D - O) / D), re.escape("%.1f%%" % (100 * (D - O) / D))),
        "overlap_pct_of_G": ("%d observed (%.1f%%)" % (O, 100 * O / G), re.escape("%d observed (%.1f%%)" % (O, 100 * O / G))),
        "overlap_pct_of_D": ("%d observed (%.1f%%)" % (O, 100 * O / D), re.escape("%d observed (%.1f%%)" % (O, 100 * O / D))),
        "expected_pct_of_G": ("%.1f%%" % (100 * E / G), r"expected[^\n]*" + re.escape("%.1f%%" % (100 * E / G))),
        "expected_pct_of_D": ("%.1f%%" % (100 * E / D), r"expected[^\n]*" + re.escape("%.1f%%" % (100 * E / D))),
        "no_susie_degs": (str(R1B["deg_no_susie_result_in_U_full"]), r"\b%d\b[^\n]*COLOC" % R1B["deg_no_susie_result_in_U_full"]),
        "deg_rule": ("BH q<0.05, |log2FC|>0.50", r"q\s*<\s*0\.05[^\n]*0\.5"),
        "n_coloc_list": (str(G), r"\b%d\b" % G),
        "n_deg_list": (str(D), r"\b%d\b" % D),
        "n_universe": ("3,769", r"3,?769"),
        "axis_title": ("x-axis title", r"[Ff]raction|[Pp]ercent|%\s*of"),
        "MH_OR_or_CI": ("OR/CI or P", r"OR|CI|P\s*=|P="),
    }
    out = []
    for k, (want, pat) in expected.items():
        found = re.search(pat, txt)
        out.append(dict(file=path, panel="1B", label=k, expected=want, found=bool(found), text=found.group(0) if found else ""))
    # any printed numbers that are not in the expected set
    nums = re.findall(r"\d+(?:\.\d+)?%?", txt)
    allowed = {str(G - O), "%.1f%%" % (100 * (G - O) / G), str(D - O), "%.1f%%" % (100 * (D - O) / D), str(O),
               "%.1f%%" % (100 * O / G), "%.1f%%" % (100 * O / D), "%.1f%%" % (100 * E / G), "%.1f%%" % (100 * E / D),
               str(R1B["deg_no_susie_result_in_U_full"]), "0.05", "0.50", "0.5", "2", "0%", "25%", "50%", "75%", "100%",
               str(G), str(D), "3", "769", "10"}
    stray = [n for n in nums if n not in allowed]
    out.append(dict(file=path, panel="1B", label="unexpected_numbers", expected="none", found=not stray, text=",".join(stray)))
    return out


def label_checks_1c(path, ws, bar_rows):
    cls = R1C["class_counts_positive436"]
    txt = "\n".join(w["text"] for w in ws)
    out = []
    # value labels next to bars: a numeric word on the same row as the bar, right of the axis
    for yc, k, lab, r in bar_rows:
        near = [w for w in ws if abs((w["top"] + w["bottom"]) / 2 - yc) < 4 and w["x0"] >= r["x0"] - 1 and re.fullmatch(r"\d+", w["text"].strip())]
        want = {"enzyme": cls["enzyme_only"], "both": cls["both"], "direct": cls["direct_only"]}.get(k)
        got = [int(w["text"]) for w in near]
        out.append(dict(file=path, panel="1C", label=f"value_label_{k}", expected=str(want), found=(want in got), text=",".join(map(str, got))))
        # label must sit near the bar end
        for w in near:
            out.append(dict(file=path, panel="1C", label=f"value_label_{k}_gap_pt", expected="<6", found=bool(abs(w["x0"] - r["x1"]) < 6),
                            text="%.2f" % (w["x0"] - r["x1"])))
    for k, pat in (("total_436", r"\b436\b"), ("axis_title", r"genes|Genes"), ("call_rule_PPH4", r"PP\.?H4|0\.5"),
                   ("trait_class_enzyme_label", r"[Ee]nzyme|ALT"), ("trait_class_direct_label", r"[Dd]irect|NAFLD|PDFF"),
                   ("both_label", r"[Bb]oth"), ("tested_denominator", r"tested|of \d")):
        m = re.search(pat, txt)
        out.append(dict(file=path, panel="1C", label=k, expected=pat, found=bool(m), text=m.group(0) if m else ""))
    return out


def text_style(path, panel, ws):
    sizes = sorted({round(w["size"], 2) for w in ws})
    fonts = sorted({w["fontname"] for w in ws})
    cols = sorted({hexcol(w["non_stroking_color"]) for w in ws})
    bold = [w["text"] for w in ws if re.search(r"Bold|Black|Heavy", w["fontname"])]
    return dict(file=path, panel=panel, sizes=sizes, fonts=fonts, text_colours=cols, bold_words=bold)


def main():
    prefix = sys.argv[1]
    marks, labels, axes, styles = [], [], [], []
    for path in sys.argv[2:]:
        pdf = pdfplumber.open(path)
        p = pdf.pages[0]
        ws = words(p)
        txt = " ".join(w["text"] for w in ws)
        kind = "1B" if re.search(r"DEG|observed|Fraction", txt) else "1C"
        if kind == "1B":
            ax, rws, legend = check_1b(path, p, ws, marks)
            labels += label_checks_1b(path, ws, rws)
            axes.append(dict(file=path, panel=kind, **{k: v for k, v in ax.items()}, legend_keys=[r["col"] for r in legend],
                             page=[p.width, p.height]))
        else:
            ax, rws, rng = check_1c(path, p, ws, marks)
            labels += label_checks_1c(path, ws, rws)
            axes.append(dict(file=path, panel=kind, **{k: v for k, v in ax.items()}, axis_range=rng, page=[p.width, p.height]))
        st = text_style(path, kind, ws)
        st["words"] = [(round(w["x0"], 1), round(w["top"], 1), w["text"], hexcol(w["non_stroking_color"])) for w in ws]
        styles.append(st)
    import pandas as pd
    pd.DataFrame(marks).to_csv(f"{prefix}_marks.tsv", sep="\t", index=False)
    pd.DataFrame(labels).to_csv(f"{prefix}_labels.tsv", sep="\t", index=False)
    json.dump(dict(axes=axes, styles=styles), open(f"{prefix}_axes_styles.json", "w"), indent=1, default=str)
    print(pd.DataFrame(marks).to_string())
    print(pd.DataFrame(labels).to_string())
    for a in axes:
        print({k: v for k, v in a.items() if k not in ("tick_x",)})
    for s in styles:
        print({k: v for k, v in s.items() if k != "words"})


if __name__ == "__main__":
    main()

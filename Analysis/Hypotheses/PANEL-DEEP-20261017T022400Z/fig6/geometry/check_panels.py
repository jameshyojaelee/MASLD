#!/usr/bin/env python3
"""Geometry and printed-label checks for Figure 6 panels against recomputed values.

Usage: check_panels.py <panel_dir> <out_prefix>
Reads only the panel PDFs in <panel_dir> and the recomputed values in ../recompute/.
Writes <out_prefix>_marks.tsv (one row per checked mark) and <out_prefix>_summary.json.
"""
import sys, os, re, json, math
import numpy as np
import pandas as pd
import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
REC = os.path.join(HERE, "..", "recompute")

# ---------------- recomputed expectations ----------------
mat = pd.read_csv(os.path.join(REC, "fig6b_matrix_recomputed_crosswalk.tsv"), sep="\t", index_col=0)
cls = pd.read_csv(os.path.join(REC, "fig6d_class_counts_recomputed_crosswalk.tsv"), sep="\t", index_col=0).iloc[:, 0]
uds = pd.read_csv(os.path.join(REC, "fig6d_untested_split_recomputed_crosswalk.tsv"), sep="\t", index_col=0).iloc[:, 0]
we = pd.read_csv(os.path.join(REC, "fig6c_worked_entries_recomputed.tsv"), sep="\t").set_index("symbol")
calls = pd.read_csv(os.path.join(REC, "catalog_calls_recomputed.tsv.gz"), sep="\t",
                    usecols=["susie_state", "gen_call"])
N_TOTAL = int(mat.loc["All", "All"])
N_NO_T12 = int((calls.susie_state.isna()).sum())
N_SUSIE_UNTESTABLE = int((calls.susie_state == "untestable").sum())
CELL_CLASS = {("supported", "supported"): "both supported", ("supported", "indeterminate"): "genetically anchored",
              ("indeterminate", "supported"): "established state associated", ("indeterminate", "indeterminate"): "unresolved"}
CLASS_KEYS = {"untested": "untested", "unresolved": "unresolved", "established": "established_state_associated",
              "genetically": "genetically_anchored", "both": "concordant", "concordant": "concordant"}


def hexcol(c):
    if c is None:
        return "none"
    if isinstance(c, (int, float)):
        c = (c, c, c)
    if len(c) == 1:
        c = (c[0], c[0], c[0])
    if len(c) == 4:  # CMYK
        C, M, Y, K = c
        c = ((1 - C) * (1 - K), (1 - M) * (1 - K), (1 - Y) * (1 - K))
    return "#%02X%02X%02X" % tuple(int(round(255 * v)) for v in c[:3])


def text_lines(chars, tol=1.0, gap=3.5):
    """Group chars into lines by baseline 'top'; returns list of (top, x0, x1, text)."""
    chars = [c for c in chars if c["text"].strip() or c["text"] == " "]
    chars = sorted(chars, key=lambda c: (round(c["top"], 1), c["x0"]))
    lines = []
    for c in chars:
        for L in lines:
            if abs(L["top"] - c["top"]) <= tol and c["upright"] == L["upright"]:
                L["chars"].append(c)
                break
        else:
            lines.append({"top": c["top"], "upright": c["upright"], "chars": [c]})
    out = []
    for L in lines:
        cs = sorted(L["chars"], key=lambda c: c["x0"])
        segs, cur = [], [cs[0]]
        for c in cs[1:]:
            if c["x0"] - cur[-1]["x1"] > gap:
                segs.append(cur); cur = [c]
            else:
                cur.append(c)
        segs.append(cur)
        for sg in segs:
            s, prev = "", None
            for c in sg:
                if prev is not None and c["x0"] - prev["x1"] > 0.6:
                    s += " "
                s += c["text"]
                prev = c
            out.append((min(c["top"] for c in sg), sg[0]["x0"], sg[-1]["x1"], re.sub(r"\s+", " ", s).strip(), L["upright"]))
    return sorted(out, key=lambda t: (t[0], t[1]))


def in_box(c, b, pad=0.0):
    cx = (c["x0"] + c["x1"]) / 2
    cy = (c["top"] + c["bottom"]) / 2
    return b[0] - pad <= cx <= b[2] + pad and b[1] - pad <= cy <= b[3] + pad


def num(s):
    return float(s.replace(",", ""))


def rotated_text(page):
    cs = [c for c in page.chars if not c["upright"]]
    if not cs:
        return ""
    cs = sorted(cs, key=lambda c: -c["top"])  # bottom-to-top for 90-degree rotated text
    return "".join(c["text"] for c in cs)


def fonts(page):
    out = {}
    for c in page.chars:
        # pdfplumber size = glyph box height; for 90-degree rotated text the em size is the box width
        sz = c["size"] if c["upright"] else (c["x1"] - c["x0"])
        key = (c["fontname"].split("+")[-1], round(sz, 2), hexcol(c.get("non_stroking_color")))
        out[str(key)] = out.get(str(key), 0) + 1
    return out


rows = []


def rec(panel, mark, kind, drawn, expected, tol, note=""):
    if isinstance(drawn, (int, float)) and isinstance(expected, (int, float)) and drawn is not None and expected is not None:
        off = abs(drawn - expected) > tol
    else:
        off = str(drawn) != str(expected)
    rows.append(dict(panel=panel, mark=mark, kind=kind, drawn=drawn, expected=expected, tol=tol, off=bool(off), note=note))


# ---------------- 6B ----------------
def check_6b(path):
    P = "6B"
    pg = pdfplumber.open(path).pages[0]
    W, H = pg.width, pg.height
    cells = [r for r in pg.rects if not (r["x0"] < 1 and r["top"] < 1 and r["width"] > W - 2)]
    xs = sorted(set(round(r["x0"], 1) for r in cells))
    ys = sorted(set(round(r["top"], 1) for r in cells))
    lines = text_lines(pg.chars)
    col_lab, row_lab = {}, {}
    for xi in xs:
        x1 = max(r["x1"] for r in cells if round(r["x0"], 1) == xi)
        hdr = [t for t in lines if t[0] < min(ys) and t[1] >= xi - 5 and t[2] <= x1 + 5 and t[4]]
        col_lab[xi] = " ".join(t[3] for t in hdr)
    for yi in ys:
        y1 = max(r["bottom"] for r in cells if round(r["top"], 1) == yi)
        lab = [t for t in lines if yi - 2 <= t[0] <= y1 and t[2] <= min(xs) + 1 and t[4]]
        row_lab[yi] = " ".join(t[3] for t in lab)

    def call_of(s):
        s = s.lower()
        for k in ("indeterminate", "untestable", "supported"):
            if s.startswith(k) or (k in s and not any(s.startswith(j) for j in ("indeterminate", "untestable", "supported"))):
                return k
        return None
    summary = {"col_labels": col_lab, "row_labels": row_lab}
    for r in cells:
        b = (r["x0"], r["top"], r["x1"], r["bottom"])
        txt = [t for t in text_lines([c for c in pg.chars if in_box(c, b)])]
        nums = [t[3] for t in txt if re.fullmatch(r"[\d,]+", t[3])]
        label = " ".join(t[3] for t in txt if not re.fullmatch(r"[\d,]+", t[3]))
        gc = call_of(row_lab[round(r["top"], 1)])
        rc = call_of(col_lab[round(r["x0"], 1)])
        exp = int(mat.loc[gc, rc]) if gc and rc else None
        drawn = num(nums[0]) if nums else None
        rec(P, f"cell gen={gc} rna={rc}", "printed count", drawn, exp, 0)
        exp_cls = CELL_CLASS.get((gc, rc), "untested")
        rec(P, f"cell gen={gc} rna={rc}", "printed class label", label.split(" (")[0].strip(), exp_cls, 0,
            note="label text: " + label)
        fill = hexcol(r.get("non_stroking_color"))
        rows.append(dict(panel=P, mark=f"cell gen={gc} rna={rc}", kind="fill colour", drawn=fill,
                         expected="gray family if untested else tinted", tol=0,
                         off=bool((exp_cls == "untested") == (fill not in ("#EEEEEE", "#E0E0E0", "#9E9E9E", "#F0F0F0"))),
                         note=f"stroke {hexcol(r.get('stroking_color'))}"))
    # margins / row-label numbers
    rl = " ".join(row_lab.values())
    m = re.search(r"COLOC\s*record:\s*([\d,]+)", rl)
    rec(P, "row label: no Tier-1/2 COLOC record", "printed count", num(m.group(1)) if m else None, N_NO_T12, 0)
    m = re.search(r"one study:\s*([\d,]+)", rl)
    rec(P, "row label: SuSiE untestable", "printed count", num(m.group(1)) if m else None, N_SUSIE_UNTESTABLE, 0)
    alltxt = " ".join(t[3] for t in lines)
    summary["n_printed"] = bool(re.search(r"28,233|n\s*=", alltxt))
    summary["margins_printed"] = any(s in alltxt for s in ("436", "18,478", "9,319", "1,853", "22,148", "4,232"))
    summary["axis_titles"] = {"x(top)": " ".join(t[3] for t in lines if t[0] < 8 and t[4]), "y(rotated)": rotated_text(pg)}
    summary["fonts"] = fonts(pg)
    summary["fills"] = sorted(set(hexcol(r.get("non_stroking_color")) for r in cells))
    return summary


# ---------------- 6D ----------------
def check_6d(path):
    P = "6D"
    pg = pdfplumber.open(path).pages[0]
    W, H = pg.width, pg.height
    hl = [l for l in pg.lines if abs(l["top"] - l["bottom"]) < 0.01]
    vl = [l for l in pg.lines if abs(l["x0"] - l["x1"]) < 0.01]
    axis = max(hl, key=lambda l: l["x1"] - l["x0"])
    ticks = sorted([l for l in vl if l["top"] >= axis["top"] - 0.01 and (l["bottom"] - l["top"]) < 6], key=lambda l: l["x0"])
    lines = text_lines(pg.chars)
    words = pg.extract_words()
    tl = [w for w in words if w["top"] > axis["top"] and w["top"] < axis["top"] + 12 and re.fullmatch(r"[\d,]+", w["text"])]
    pairs = []
    for t in ticks:
        w = min(tl, key=lambda w: abs((w["x0"] + w["x1"]) / 2 - t["x0"])) if tl else None
        if w is not None and abs((w["x0"] + w["x1"]) / 2 - t["x0"]) < 6:
            pairs.append((num(w["text"]), t["x0"]))
    v = np.array([p[0] for p in pairs]); x = np.array([p[1] for p in pairs])
    b, a = np.polyfit(v, x, 1)
    resid = x - (a + b * v)
    summary = {"ticks": pairs, "transform": {"x0_pt": a, "pt_per_gene": b, "max_tick_resid_pt": float(np.abs(resid).max())}}
    axis_range = (v.max() - v.min())
    tol = max(0.005 * axis_range, 1.0)
    bars = [r for r in pg.rects if r["height"] < 40 and r["width"] < W - 2 and hexcol(r.get("non_stroking_color")) not in ("#FFFFFF",)]
    # category labels: text lines left of axis start
    labs = [t for t in lines if t[2] < a - 1 and t[4]]
    bands = {}
    for r in bars:
        yc = (r["top"] + r["bottom"]) / 2
        near = sorted(labs, key=lambda t: abs(t[0] + 2.5 - yc))[:2]
        near = sorted(near, key=lambda t: t[0])
        key = near[0][3]
        bands.setdefault(key, []).append(r)
    for key, rs in bands.items():
        rs = sorted(rs, key=lambda r: r["x0"])
        k = next((CLASS_KEYS[w] for w in CLASS_KEYS if key.lower().startswith(w)), None)
        end_val = (rs[-1]["x1"] - a) / b
        rec(P, f"bar end: {key}", "bar length (inverted axis)", round(end_val, 1), int(cls.get(k, -1)), tol,
            note="fills " + ",".join(hexcol(r.get("non_stroking_color")) for r in rs))
        if len(rs) > 1:
            segtxt = [t for t in lines if any(in_box({"x0": t[1], "x1": t[2], "top": t[0], "bottom": t[0] + 6}, (r["x0"], r["top"], r["x1"], r["bottom"])) for r in rs)]
            for i, r in enumerate(rs):
                seg_val = (r["x1"] - r["x0"]) / b
                inside = [t[3] for t in lines if r["x0"] <= t[1] and t[2] <= r["x1"] + 0.5 and r["top"] - 1 <= t[0] <= r["bottom"]]
                name = " ".join(s for s in inside if not re.fullmatch(r"[\d,]+", s))
                dom = "no_rna_call_not_detected" if "rna" in name.lower() else ("no_genetic_call" if "genetic" in name.lower() else None)
                exp = int(uds.get(dom, -1)) if dom else None
                rec(P, f"segment {i}: {name}", "segment length (inverted axis)", round(seg_val, 1), exp, tol,
                    note=f"fill {hexcol(r.get('non_stroking_color'))}")
                pn = [s for s in inside if re.fullmatch(r"[\d,]+", s)]
                rec(P, f"segment {i}: {name}", "printed count", num(pn[0]) if pn else None, exp, 0)
        route = {"untested": "acquire the missing decisive assay", "unresolved": "measure the single most decisive missing layer",
                 "established_state_associated": "gene knockdown and rescue", "genetically_anchored": "allele-aware regulatory editing",
                 "concordant": "allele-aware edit plus gene perturbation"}.get(k)
        yc0 = (rs[0]["top"] + rs[0]["bottom"]) / 2
        lab2 = [t[3] for t in labs if abs(t[0] + 2.5 - yc0) < 9]
        rec(P, f"route label: {key}", "printed route vs rule", route if route and any(route == l for l in lab2) else " / ".join(lab2), route, 0)
        # printed value right of bar
        yc = (rs[0]["top"] + rs[0]["bottom"]) / 2
        right = [w for w in words if w["x0"] >= rs[-1]["x1"] - 0.5 and abs((w["top"] + w["bottom"]) / 2 - yc) < 4 and re.fullmatch(r"[\d,]+", w["text"])]
        rec(P, f"bar label: {key}", "printed count", num(right[0]["text"]) if right else None, int(cls.get(k, -1)), 0)
    alltxt = " ".join(t[3] for t in lines)
    summary["axis_title"] = " ".join(t[3] for t in lines if t[0] > axis["top"] + 8)
    summary["n_printed"] = bool(re.search(r"28,233|n\s*=", alltxt))
    summary["fonts"] = fonts(pg)
    summary["fills"] = sorted(set(hexcol(r.get("non_stroking_color")) for r in bars))
    summary["categories"] = list(bands.keys())
    summary["tol_genes"] = tol
    return summary


# ---------------- 6C ----------------
def fmt_pp(v):
    if v > 0.999:
        return ">0.999"
    return f"{v:.3f}"


def fmt_q(q):
    return f"{q:.1e}".replace("e-0", "e-")


def check_6c(path):
    P = "6C"
    pg = pdfplumber.open(path).pages[0]
    lines = text_lines(pg.chars)
    gw = [w for w in pg.extract_words(extra_attrs=["fontname"]) if "Oblique" in w["fontname"] or "Italic" in w["fontname"]]
    gw = sorted(gw, key=lambda w: w["top"])
    hdr = [t for t in lines if t[0] < gw[0]["top"] - 8]
    # column starts from header segments on the first header line
    starts = []
    for x0 in sorted(t[1] for t in hdr):
        if not starts or x0 - starts[-1] > 3:
            starts.append(x0)
    names = []
    for x0 in starts:
        names.append(" ".join(t[3] for t in sorted(hdr, key=lambda t: t[0]) if abs(t[1] - x0) < 3))
    summary = {"genes": [w["text"] for w in gw], "columns": names}

    def col_of(x):
        j = max([i for i, s0 in enumerate(starts) if x >= s0 - 3] or [0])
        return j
    for i, w in enumerate(gw):
        g = w["text"]
        lo = (gw[i - 1]["top"] + w["top"]) / 2 if i else w["top"] - 12
        hi = (gw[i + 1]["top"] + w["top"]) / 2 if i + 1 < len(gw) else pg.height
        cols = {}
        for t in lines:
            if lo <= t[0] < hi:
                cols.setdefault(col_of(t[1]), []).append(t[3])
        cols = {k: " ".join(v) for k, v in cols.items()}
        gen_txt, rna_txt, cls_txt = cols.get(1, ""), cols.get(2, ""), cols.get(3, "")
        route_txt = cols.get(4, "")
        summary[g] = cols
        if g not in we.index:
            rec(P, g, "gene present in recompute", g, "one of GLP1R/HKDC1/MTARC1/THRB", 0)
            continue
        e = we.loc[g]
        rec(P, g, "genetics call", gen_txt.split()[0] if gen_txt else None, e.gen_call, 0, note=gen_txt)
        if e.gen_call == "supported":
            exp = fmt_pp(e.susie_best) + " SuSiE, " + e.susie_study
            m = re.search(r"(>?\d\.\d+)\s*SuSiE,\s*(\S+)", gen_txt)
            rec(P, g, "best SuSiE PP.H4 + study", (m.group(1) + " SuSiE, " + m.group(2)) if m else None, exp, 0)
        else:
            exp = f"ABF {e.abf_best:.3f}, {e.abf_study}"
            m = re.search(r"ABF\s*(\d\.\d+),\s*(\S+)", gen_txt)
            rec(P, g, "best ABF PP.H4 + study", (f"ABF {m.group(1)}, {m.group(2)}") if m else None, exp, 0)
        rec(P, g, "RNA call", rna_txt.split()[0] if rna_txt else None, e.rna_call, 0, note=rna_txt)
        if e.rna_call != "untestable":
            m = re.search(r"([+-]?\d\.\d+)\s*\(q\s*([\d.e-]+)\)", rna_txt)
            exp = f"{e.rna_log2fc:+.2f} (q {fmt_q(e.rna_q)})"
            rec(P, g, "RNA log2FC (BH q)", (m.group(1) + " (q " + m.group(2) + ")") if m else None, exp, 0)
        disp = {"untested": "Untested", "concordant": "Both supported", "unresolved": "Unresolved",
                "genetically_anchored": "Genetically anchored", "established_state_associated": "Established state associated"}[e.class_recomputed]
        rec(P, g, "evidence class", cls_txt, disp, 0)
        gene_route_kw = {"untested": "missing decisive assay", "concordant": "gene perturbation", "unresolved": "decisive missing layer",
                         "genetically_anchored": "Allele-aware", "established_state_associated": "knockdown"}[e.class_recomputed]
        rec(P, g, "gene route keyword", gene_route_kw if gene_route_kw.lower() in route_txt.lower() else "missing: " + route_txt,
            gene_route_kw, 0, note=route_txt)
        obj_kw = "missing assay" if e.class_recomputed == "untested" else "perturbation"
        rec(P, g, "object route keyword", obj_kw if obj_kw in route_txt.split("\n")[-1].lower() else "missing: " + route_txt, obj_kw, 0,
            note="sidecar: GLP1R EXP_MEASURE_MISSING_ASSAY_V2; others EXP_PROTEIN_STATE_CONTEXT_V2 (gene perturbation with orthogonal rescue)")
    alltxt = " ".join(t[3] for t in lines)
    summary["n_printed"] = bool(re.search(r"774", alltxt))
    summary["fonts"] = fonts(pg)
    return summary


# ---------------- 6A / 6E schematics ----------------
def check_schematic(path, P):
    pg = pdfplumber.open(path).pages[0]
    lines = text_lines(pg.chars)
    boxes = [r for r in pg.rects if r["width"] < pg.width - 2]
    out = {"text": [t[3] for t in lines], "box_strokes": [hexcol(r.get("stroking_color")) for r in boxes],
           "n_curves(arrows)": len(pg.curves), "fonts": fonts(pg),
           "numbers_printed": re.findall(r"\b\d[\d,]*\b", " ".join(t[3] for t in lines))}
    if P == "6A":
        full = " ".join(t[3] for t in lines)
        rows.append(dict(panel=P, mark="object types listed", kind="schematic claim",
                         drawn="DNA element listed" if "DNA element" in full else "no DNA element",
                         expected="0 DNA-element records (protein_coding_gene 14,176; untestable_nomination 13,576; lncrna_transcript 423; ambiguous_lncrna_locus 58)",
                         tol=0, off="DNA element" in full, note="molecular_object_sidecar.jsonl"))
    return out


def main():
    d, prefix = sys.argv[1], sys.argv[2]
    f = {k: os.path.join(d, v) for k, v in {
        "6A": "fig6a_catalog_structure.pdf", "6B": "fig6b_coverage_observability.pdf", "6C": "fig6c_worked_entries.pdf",
        "6D": "fig6d_next_experiment.pdf", "6E": "fig6e_portal_source_access.pdf"}.items()}
    S = {}
    S["6A"] = check_schematic(f["6A"], "6A")
    S["6B"] = check_6b(f["6B"])
    S["6C"] = check_6c(f["6C"])
    S["6D"] = check_6d(f["6D"])
    S["6E"] = check_schematic(f["6E"], "6E")
    df = pd.DataFrame(rows)
    df.to_csv(prefix + "_marks.tsv", sep="\t", index=False)
    S["counts"] = df.groupby("panel").off.agg(["count", "sum"]).rename(columns={"count": "checked", "sum": "off"}).to_dict("index")
    json.dump(S, open(prefix + "_summary.json", "w"), indent=1, default=str)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
    print(df[["panel", "mark", "kind", "drawn", "expected", "off"]].to_string())
    print(json.dumps(S["counts"]))
    for k in ("6B", "6D", "6C"):
        print(k, {kk: vv for kk, vv in S[k].items() if kk not in ("fonts",)})
        print(k, "fonts", S[k]["fonts"])
    for k in ("6A", "6E"):
        print(k, S[k])


if __name__ == "__main__":
    main()

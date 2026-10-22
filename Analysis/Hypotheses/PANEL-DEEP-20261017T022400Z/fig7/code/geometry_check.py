"""Geometry check for the Figure 7 panels: build each x-axis transform from its tick marks and tick labels,
invert every data mark (point centre, interval ends, bound and reference glyphs, heatmap cells, bar widths),
and compare with this review's recomputed values.

Usage: python geometry_check.py <out_dir> <prefix: main|calib>
Writes <out_dir>/geometry_<prefix>.tsv (one row per mark) and geometry_<prefix>_summary.tsv.
A mark is 'off' when |drawn - recomputed| exceeds 0.5 pt on the page, converted to data units with the axis scale.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(sys.argv[1])
PREFIX = sys.argv[2]
NUM = re.compile(r"^[−\-–]?\d+(\.\d+)?$")
TOL_PT = 0.5
MC_TOL = {}  # (panel,row) -> extra tolerance in data units for bootstrap endpoints computed with a different seed


def val(t):
    return float(t.replace("−", "-").replace("–", "-"))


def load(panel):
    return json.load(open(OUT / f"geom_{PREFIX}_{panel}.json"))


def axes(g):
    """Return list of axes: dict(y, x0, x1, a, b, resid_pt, ticks) with value = a + b * x."""
    out = []
    lines = g["lines"]
    axl = [l for l in lines if l["stroke"] == "#000000" and abs(l["top"] - l["bottom"]) < 0.01 and l["x1"] - l["x0"] > 30
           and "2.56" not in l.get("dash", "")]
    for ax in axl:
        y = ax["top"]
        ticks = [l for l in lines if l["stroke"] == "#000000" and abs(l["x0"] - l["x1"]) < 0.01 and abs(l["top"] - y) < 0.05
                 and 0.5 < l["bottom"] - l["top"] < 3]
        labs = [w for w in g["words"] if NUM.match(w["text"]) and 0 < w["top"] - y < 6 and ax["x0"] - 8 < (w["x0"] + w["x1"]) / 2 < ax["x1"] + 8]
        pairs = []
        for w in labs:
            xc = (w["x0"] + w["x1"]) / 2
            t = min(ticks, key=lambda l: abs(l["x0"] - xc)) if ticks else None
            tx = t["x0"] if t is not None and abs(t["x0"] - xc) < 1.5 else None
            pairs.append((val(w["text"]), tx if tx is not None else xc, tx is not None, w["text"]))
        if len(pairs) < 2:
            continue
        v = np.array([p[0] for p in pairs]); x = np.array([p[1] for p in pairs])
        b, a = np.polyfit(x, v, 1)
        resid = np.abs(a + b * x - v) / abs(b)
        out.append({"y": y, "x0": ax["x0"], "x1": ax["x1"], "a": a, "b": b, "resid_pt": float(resid.max()),
                    "labels": ";".join(p[3] for p in pairs), "ticks_matched": int(sum(p[2] for p in pairs))})
    return sorted(out, key=lambda d: (d["y"], d["x0"]))


def axis_for(axs, x, y):
    cand = [a for a in axs if a["y"] >= y - 0.5 and a["x0"] - 3 <= x <= a["x1"] + 3]
    return min(cand, key=lambda a: a["y"]) if cand else None


def marks(g, ymax=None):
    """Points (curves), intervals (coloured horizontal lines), squares (small filled rects)."""
    res = []
    for c in g["curves"]:
        w = c["x1"] - c["x0"]
        if 2 < w < 4 and len(c.get("pts", [])) == 5:
            res.append({"kind": "point", "x": (c["x0"] + c["x1"]) / 2, "y": (c["top"] + c["bottom"]) / 2,
                        "stroke": c["stroke"], "fill": c["fill"] if c["filled"] else "none"})
    for l in g["lines"]:
        if abs(l["top"] - l["bottom"]) < 0.01 and l["stroke"] not in ("#000000",) and l["x1"] - l["x0"] > 0.01 and l["lw"] > 0.7:
            res.append({"kind": "interval", "x": l["x0"], "x_hi": l["x1"], "y": l["top"], "stroke": l["stroke"], "fill": ""})
    for r in g["rects"]:
        if r["filled"] and 2 < r["x1"] - r["x0"] < 3.5 and 2 < r["bottom"] - r["top"] < 3.5:
            res.append({"kind": "square", "x": (r["x0"] + r["x1"]) / 2, "y": (r["top"] + r["bottom"]) / 2, "stroke": r["stroke"], "fill": r["fill"]})
    if ymax is not None:
        res = [m for m in res if m["y"] < ymax]
    return res


def row_label(g, y, xmax):
    ws = [w for w in g["words"] if abs((w["top"] + w["bottom"]) / 2 - y) < 3.2 and w["x1"] < xmax]
    return " ".join(w["text"] for w in sorted(ws, key=lambda w: w["x0"]))


# ---------------- expected values (this review's recomputation)
mp = pd.read_csv(OUT / "mpra_recomputed.tsv", sep="\t").set_index("candidate")
eq = pd.read_csv(OUT / "eqtl_7b_rows_recomputed.tsv", sep="\t").set_index("row")
es = pd.read_csv(OUT / "eqtl_strata_recomputed.tsv", sep="\t").set_index("stratum")
fr = json.load(open(OUT / "fibrosis_routing_rrr_recomputed.json"))
em = json.load(open(OUT / "eqtl_misc_recomputed.json"))
sv = pd.read_csv(OUT / "severity_transfer_recomputed.tsv", sep="\t")
cc = pd.read_csv(OUT / "cellclass_scores_recomputed.tsv", sep="\t").set_index("model")
cg = pd.read_csv(OUT / "cellclass_gains_recomputed.tsv", sep="\t").set_index("model")
e7 = pd.read_csv(OUT / "panel7e_recomputed.tsv", sep="\t").set_index("row")
CHANCE = float(es.loc["all", "chance_median"])


def mprow(k):
    r = mp.loc[k]
    return (r.gain_vs_allele_ridge, r.gain_ci95_low, r.gain_ci95_high)


h3 = fr["rna_to_h3k27ac"]; at = fr["rna_to_atac_transfer"]; c1 = fr["graded_fibrosis"]["C1_c_minus_a"]
nt_m1 = eq.loc["nucleotide_transformer signAUC m1_dist_emb - b1_dist"]
EXP = {
    "fig7B_gain_over_task_native_baseline": [
        ("RRR rank 4", (h3["skill"], h3["fold_min"], h3["fold_max"])),
        ("RRR rank 4, transferred", (at["skill"], at["recorded_ci"][0], at["recorded_ci"][1])),
        ("hyenadna", mprow("hyenadna::delta_ridge")), ("sei", mprow("sei::native_40class_ridge")),
        ("nucleotide transformer", mprow("nucleotide_transformer::delta_ridge")), ("caduceus", mprow("caduceus::delta_ridge")),
        ("dnabert2", mprow("dnabert2::delta_ridge")),
        ("Enformer SAR", mprow("enformer_crested_restricted_port::rc_ensemble_sar_ridge")),
        ("Enformer SAD", mprow("enformer_crested_restricted_port::rc_ensemble_sad_ridge")),
        ("hyenadna", (eq.loc["hyenadna |z| m1_dist_emb - b1_dist", "point"],) + tuple(eq.loc["hyenadna |z| from producer draws", ["lo", "hi"]])),
        ("nucleotide transformer", (eq.loc["nucleotide_transformer |z| m1_dist_emb - b1_dist", "point"],) + tuple(eq.loc["nucleotide_transformer |z| from producer draws", ["lo", "hi"]])),
        ("nucleotide transformer", (nt_m1.m1_auc - CHANCE, nt_m1.m1_auc_lo - CHANCE, nt_m1.m1_auc_hi - CHANCE)),
        ("composition + expression", (c1["point"], c1["ci95"][0], c1["ci95"][1])),
        ("routed system", (fr["routing"]["routed_minus_ref"], 0.0, 0.0)),
    ],
    "fig7C_score_beside_training_exposure": [
        ("H3K27ac, fit here", (h3["skill"], h3["fold_min"], h3["fold_max"])),
        ("ATAC, transferred", (at["skill"], at["recorded_ci"][0], at["recorded_ci"][1])),
        ("hyenadna", mprow("hyenadna::delta_ridge")), ("sei", mprow("sei::native_40class_ridge")),
        ("nucleotide transformer", mprow("nucleotide_transformer::delta_ridge")), ("caduceus", mprow("caduceus::delta_ridge")),
        ("dnabert2", mprow("dnabert2::delta_ridge")),
        ("Enformer SAR", mprow("enformer_crested_restricted_port::rc_ensemble_sar_ridge")),
        ("Enformer SAD", mprow("enformer_crested_restricted_port::rc_ensemble_sad_ridge")),
        ("TranscriptFormer", (cc.loc["transcriptformer_tf_sapiens::two_layer_mlp", "all_studies"], None, None)),
        ("HVG+PCA kNN", (cc.loc["hvg_pca_knn", "all_studies"], None, None)),
        ("HVG+PCA SVM", (cc.loc["hvg_pca_linear_svm", "all_studies"], None, None)),
        ("HVG+PCA logistic", (cc.loc["hvg_pca_logistic", "all_studies"], None, None)),
        ("HVG+PCA elastic net", (cc.loc["hvg_pca_elastic_net", "all_studies"], None, None)),
        ("HVG+PCA centroid", (cc.loc["hvg_pca_nearest_centroid", "all_studies"], None, None)),
    ],
}
svk = sv.set_index("unit")
d_rows = []
for s in ["all", "z_ge_1", "z_ge_2", "z_ge_3", "z_ge_5"]:
    r = es.loc[s]
    d_rows.append((s, {"achieved": (r.best_auc, r.auc_lo, r.auc_hi), "bound": r.seq_bound, "reference": r.chance_median}))
for u in ["participant record", "genetic identity, primary library, discordant removed", "participant, MASLD", "participant, CVH"]:
    r = svk.loc[u]
    d_rows.append((u, {"achieved": (r.spearman, r.ci_strat_lo, r.ci_strat_hi), "bound": r.tie_ceiling}))
for s_ in ["all", "z_ge_1", "z_ge_2", "z_ge_3", "z_ge_5"]:
    # seed-to-seed Monte Carlo spread of a 2.5% bootstrap quantile is about 0.04 x half-width per run; allow 3 SE of a difference
    MC_TOL[("fig7D_measured_bounds_and_join_structure", s_)] = 0.13 * float(es.loc[s_, "resolution_half_width"])
ch = json.load(open(OUT / "eqtl_chance_model_score.json"))  # model-score null: reported as a display alternative only
for s_ in ch:
    # the drawn reference is a fixed-random-score within-block permutation median (the producer's method, reproduced
    # here with a different random score); its value moves with the realised random score, so allow 0.4 x permutation SD
    sd_ = (ch[s_]["hi"] - ch[s_]["lo"]) / 3.92
    MC_TOL[("fig7D_measured_bounds_and_join_structure", s_)] = max(MC_TOL[("fig7D_measured_bounds_and_join_structure", s_)], 0.4 * sd_)
for u in ["participant record", "genetic identity, primary library, discordant removed", "participant, MASLD", "participant, CVH"]:
    r = svk.loc[u]
    MC_TOL[("fig7D_measured_bounds_and_join_structure", u)] = 0.1 * (r.ci_strat_hi - r.ci_strat_lo) / 2
d_rows.append(("own block left in", {"bound": em["ld_join_left_in_pearson"]}))
d_rows.append(("own block held out", {"achieved": (em["ld_join_held_out_pearson"], None, None)}))
d_rows.append(("pre-declared floor", {"reference": fr["routing"]["floor"]}))
d_rows.append(("oracle router", {"bound": fr["routing"]["oracle_gain_over_ref"]}))
d_rows.append(("rule as fitted", {"achieved": (fr["routing"]["routed_minus_ref"], None, None)}))
E_ORDER = ["normal vs diseased map", "zone from expression", "NAFLD endpoint", "severity, female only", "ATAC map vs fibrosis",
           "NAFLD meta-analysis", "severity within disease", "cirrhosis endpoint", "ATAC map vs steatosis", "signed eQTL direction",
           "routing adds skill", "program h2, AST", "program h2, ALT", "treatment response"]
F_ORDER = ["geneformer_v2_316m", "transcriptformer_tf_sapiens", "scimilarity_v1_1", "uce_4l", "uce_33l", "scvi_liver_latent"]
F_EXPO = {"geneformer_v2_316m": {"GSE136103": "confounded", "GSE185477": "confounded", "Liver_Atlas": "confounded",
                                 "GSE174748": "unresolved", "GSE189600": "unresolved", "GSE202379": "unresolved", "GSE244832": "unresolved"},
          "transcriptformer_tf_sapiens": {"GSE185477": "confounded", "Liver_Atlas": "confounded"},
          "scimilarity_v1_1": {"GSE136103": "confounded", "GSE185477": "confounded"},
          "uce_4l": {"GSE185477": "confounded", "Liver_Atlas": "confounded"},
          "uce_33l": {"GSE185477": "confounded", "Liver_Atlas": "confounded"}, "hvg_pca": {}}
COL = {"#1565C0": "clean", "#C9265E": "confounded", "#9E9E9E": "unresolved"}
census = pd.read_csv(OUT / "census_recount.tsv", sep="\t")

rows = []


def rec(panel, row, label, kind, drawn, exp, scale, colour="", note=""):
    tol = TOL_PT * abs(scale)
    if kind.startswith("interval") or kind == "reference":
        tol = max(tol, MC_TOL.get((panel, row), 0.0))
    off = None if exp is None or drawn is None else abs(drawn - exp) > tol
    rows.append({"panel": panel, "row": row, "label_on_panel": label, "mark": kind, "drawn": drawn, "recomputed": exp,
                 "diff": None if exp is None or drawn is None else drawn - exp,
                 "diff_pt": None if exp is None or drawn is None else (drawn - exp) / abs(scale),
                 "off": off, "colour": colour, "note": note})


def forest(panel, expected, legend_cut):
    g = load(panel); axs = axes(g); ms = marks(g, ymax=legend_cut)
    ys = sorted({round(m["y"], 1) for m in ms})
    ys = [y for i, y in enumerate(ys) if i == 0 or y - ys[i - 1] > 1.0]
    for ax in axs:
        rows.append({"panel": panel, "row": "AXIS", "label_on_panel": ax["labels"], "mark": "axis", "drawn": ax["b"],
                     "recomputed": None, "diff": None, "diff_pt": ax["resid_pt"], "off": ax["resid_pt"] > 0.25,
                     "colour": "", "note": f"y={ax['y']:.2f}; value = {ax['a']:.5f} + {ax['b']:.6f}*x; tick residual max {ax['resid_pt']:.3f} pt; ticks matched {ax['ticks_matched']}"})
    if len(ys) != len(expected):
        rows.append({"panel": panel, "row": "ROWCOUNT", "label_on_panel": "", "mark": "rows", "drawn": len(ys), "recomputed": len(expected),
                     "diff": None, "diff_pt": None, "off": True, "colour": "", "note": "row count differs"})
    for (lab, (p, lo, hi)), y in zip(expected, ys):
        rm = [m for m in ms if abs(m["y"] - y) < 1.0]
        x_any = rm[0]["x"] if rm else 0
        ax = axis_for(axs, x_any, y)
        if ax is None:
            continue
        drawn_label = row_label(g, y, ax["x0"])
        pts = [m for m in rm if m["kind"] == "point"]
        ivs = [m for m in rm if m["kind"] == "interval"]
        for m in pts:
            rec(panel, lab, drawn_label, "point", ax["a"] + ax["b"] * m["x"], p, ax["b"], f"{m['stroke']}/{m['fill']}")
        if lo is not None:
            if ivs:
                m = max(ivs, key=lambda m: m["x_hi"] - m["x"])
                rec(panel, lab, drawn_label, "interval_lo", ax["a"] + ax["b"] * m["x"], lo, ax["b"], m["stroke"])
                rec(panel, lab, drawn_label, "interval_hi", ax["a"] + ax["b"] * m["x_hi"], hi, ax["b"], m["stroke"])
            elif lo == hi == p:
                rec(panel, lab, drawn_label, "interval (zero width, none drawn)", ax["a"] + ax["b"] * pts[0]["x"] if pts else None, lo, ax["b"], "")
            else:
                rec(panel, lab, drawn_label, "interval MISSING", None, lo, ax["b"], "", "interval expected but not drawn")
        elif ivs:
            rec(panel, lab, drawn_label, "interval UNEXPECTED", ax["a"] + ax["b"] * ivs[0]["x"], None, ax["b"], ivs[0]["stroke"])
        else:
            rows.append({"panel": panel, "row": lab, "label_on_panel": drawn_label, "mark": "no interval drawn", "drawn": None,
                         "recomputed": None, "diff": None, "diff_pt": None, "off": None, "colour": "", "note": "no uncertainty shown"})


forest("fig7B_gain_over_task_native_baseline", EXP["fig7B_gain_over_task_native_baseline"], legend_cut=180)
forest("fig7C_score_beside_training_exposure", EXP["fig7C_score_beside_training_exposure"], legend_cut=192)

# ---- 7D: achieved (filled magenta circle + interval), bound (open gray circle), reference (gray square)
panel = "fig7D_measured_bounds_and_join_structure"
g = load(panel); axs = axes(g); ms = marks(g, ymax=222)
ys = sorted({round(m["y"], 1) for m in ms}); ys = [y for i, y in enumerate(ys) if i == 0 or y - ys[i - 1] > 1.0]
for ax in axs:
    rows.append({"panel": panel, "row": "AXIS", "label_on_panel": ax["labels"], "mark": "axis", "drawn": ax["b"], "recomputed": None,
                 "diff": None, "diff_pt": ax["resid_pt"], "off": ax["resid_pt"] > 0.25, "colour": "",
                 "note": f"y={ax['y']:.2f}; tick residual max {ax['resid_pt']:.3f} pt; ticks matched {ax['ticks_matched']}"})
if len(ys) != len(d_rows):
    rows.append({"panel": panel, "row": "ROWCOUNT", "drawn": len(ys), "recomputed": len(d_rows), "off": True, "mark": "rows"})
for (lab, ex), y in zip(d_rows, ys):
    rm = [m for m in ms if abs(m["y"] - y) < 1.0]
    ax = axis_for(axs, rm[0]["x"], y); dl = row_label(g, y, ax["x0"])
    for m in rm:
        v = ax["a"] + ax["b"] * m["x"]
        if m["kind"] == "point" and m["fill"] == "#C9265E":
            rec(panel, lab, dl, "achieved", v, ex.get("achieved", (None,))[0], ax["b"], "#C9265E")
        elif m["kind"] == "point" and m["fill"] == "#FFFFFF":
            rec(panel, lab, dl, "bound", v, ex.get("bound"), ax["b"], "open #9E9E9E")
        elif m["kind"] == "square":
            rec(panel, lab, dl, "reference", v, ex.get("reference"), ax["b"], "#9E9E9E square")
        elif m["kind"] == "interval":
            a3 = ex.get("achieved", (None, None, None))
            rec(panel, lab, dl, "interval_lo", v, a3[1], ax["b"], m["stroke"])
            rec(panel, lab, dl, "interval_hi", ax["a"] + ax["b"] * m["x_hi"], a3[2], ax["b"], m["stroke"])
    kinds = {("achieved" if m["kind"] == "point" and m["fill"] == "#C9265E" else "bound" if m["kind"] == "point" else m["kind"]) for m in rm}
    for k in ex:
        if k == "achieved" and "achieved" not in kinds:
            rec(panel, lab, dl, "achieved MISSING", None, ex[k][0], ax["b"])
        if k == "bound" and "bound" not in kinds:
            rec(panel, lab, dl, "bound MISSING", None, ex[k], ax["b"])
        if k == "reference" and "square" not in kinds:
            rec(panel, lab, dl, "reference MISSING", None, ex[k], ax["b"])

# ---- 7E: one point per row
panel = "fig7E_verdict_vocabulary_and_refusal_thresholds"
g = load(panel); axs = axes(g); ms = [m for m in marks(g) if m["kind"] == "point"]
ms = sorted(ms, key=lambda m: m["y"])
for ax in axs:
    rows.append({"panel": panel, "row": "AXIS", "label_on_panel": ax["labels"], "mark": "axis", "drawn": ax["b"], "recomputed": None,
                 "diff": None, "diff_pt": ax["resid_pt"], "off": ax["resid_pt"] > 0.25, "colour": "", "note": f"tick residual max {ax['resid_pt']:.3f} pt"})
ax = axs[0]
for lab, m in zip(E_ORDER, ms):
    rec(panel, lab, row_label(g, m["y"], ax["x0"]), "point", ax["a"] + ax["b"] * m["x"], e7.loc[lab, "ratio"], ax["b"], m["fill"],
        f"verdict drawn: {e7.loc[lab, 'verdict_drawn']}")

# ---- 7F: heatmap cells and gain forest
panel = "fig7F_liver_encoder_against_the_field"
g = load(panel)
cells = [r for r in g["rects"] if r["filled"] and 25 < r["x1"] - r["x0"] < 35 and 8 < r["bottom"] - r["top"] < 12]
colx = sorted({round(r["x0"], 1) for r in cells}); rowy = sorted({round(r["top"], 1) for r in cells})
heads = [w for w in g["words"] if w["text"].startswith(("GSE", "Liver"))]
hx = sorted(heads, key=lambda w: w["x0"])
col_names = [w["text"] for w in hx]
model_rows = F_ORDER[:5] + ["hvg_pca"]
for r in cells:
    ci = colx.index(round(r["x0"], 1)); ri = rowy.index(round(r["top"], 1))
    st = col_names[ci]; mdl = model_rows[ri]
    exp = F_EXPO[mdl].get(st, "clean"); drawn = COL.get(r["fill"], r["fill"])
    rows.append({"panel": panel, "row": f"{mdl} x {st}", "label_on_panel": row_label(g, (r["top"] + r["bottom"]) / 2, 100) + " | " + st,
                 "mark": "heatmap cell", "drawn": drawn, "recomputed": exp, "diff": None, "diff_pt": None, "off": drawn != exp,
                 "colour": r["fill"], "note": "column header order read left to right; headers rotated"})
axs = axes(g); ms = marks(g, ymax=205)
for ax in axs:
    rows.append({"panel": panel, "row": "AXIS", "label_on_panel": ax["labels"], "mark": "axis", "drawn": ax["b"], "recomputed": None,
                 "diff": None, "diff_pt": ax["resid_pt"], "off": ax["resid_pt"] > 0.25, "colour": "", "note": f"x {ax['x0']:.1f}-{ax['x1']:.1f}; tick residual {ax['resid_pt']:.3f} pt"})
ys = sorted({round(m["y"], 1) for m in ms if m["y"] > 110}); ys = [y for i, y in enumerate(ys) if i == 0 or y - ys[i - 1] > 1.0]
for mdl, y in zip(F_ORDER, ys):
    for ax, h in zip(sorted(axs, key=lambda a: a["x0"]), ("linear", "two_layer_mlp")):
        rm = [m for m in ms if abs(m["y"] - y) < 1.0 and ax["x0"] - 3 <= m["x"] <= ax["x1"] + 3]
        key = f"{mdl}::{h}"; e = cg.loc[key]
        for m in rm:
            if m["kind"] == "point":
                rec(panel, key, row_label(g, y, 108 if h == "linear" else 0), "point", ax["a"] + ax["b"] * m["x"], e.gain_all_point, ax["b"], f"{m['stroke']}/{m['fill']}")
            elif m["kind"] == "interval":
                rec(panel, key, "", "interval_lo", ax["a"] + ax["b"] * m["x"], e.gain_all_lo, ax["b"], m["stroke"])
                rec(panel, key, "", "interval_hi", ax["a"] + ax["b"] * m["x_hi"], e.gain_all_hi, ax["b"], m["stroke"])

# ---- 7A: bar segment widths proportional to the printed counts and to the recount
panel = "fig7A_framework_schematic"
g = load(panel)
seg = sorted([r for r in g["rects"] if r["filled"] and r["bottom"] - r["top"] > 15 and r["fill"] != "#FFFFFF"], key=lambda r: r["x0"])
tot_w = sum(r["x1"] - r["x0"] for r in seg)
st = census.exposure_state.fillna("NONE")
unexp = int(st.str.startswith("target_label_unexposed").sum())
seen = int(st.isin(["encoder_seen", "continual_seen"]).sum())
unk = int((st == "unknown").sum())
nof = int((st == "NONE").sum())
other = int(len(st) - unexp - seen - unk - nof - (st == "NO_FINDING").sum() - (st == "NO_AUDIT_FILE").sum())
recount = [unexp, nof, unk, other, seen]
labels_in_bar = [w for w in g["words"] if NUM.match(w["text"]) and 80 < w["top"] < 95]
printed = [int(w["text"]) for w in sorted(labels_in_bar, key=lambda w: w["x0"])]
tot_printed = sum(printed)
for r, n_print, n_re, name in zip(seg, printed, recount, ["unexposed", "no finding", "unknown", "other", "encoder seen"]):
    w = r["x1"] - r["x0"]
    drawn_n = w / tot_w * tot_printed
    rows.append({"panel": panel, "row": name, "label_on_panel": str(n_print), "mark": "bar segment width (entries)", "drawn": drawn_n,
                 "recomputed": n_re, "diff": drawn_n - n_print, "diff_pt": (drawn_n - n_print) * tot_w / tot_printed,
                 "off": abs((drawn_n - n_print) * tot_w / tot_printed) > TOL_PT, "colour": r["fill"],
                 "note": f"width {w:.2f} pt; printed label {n_print}; recount from current audit files {n_re} (panel states a Sept 1 snapshot)"})
box_nums = [w["text"] for w in g["words"] if NUM.match(w["text"]) and w["top"] < 45]
rows.append({"panel": panel, "row": "boxes", "label_on_panel": " ".join(box_nums), "mark": "printed counts", "drawn": None, "recomputed": None,
             "diff": None, "diff_pt": None, "off": None, "colour": "", "note": "8 task files, 106 model dirs, 104 audit files, 24 = 17 zonation + 7 DNA-LM guards (two lanes), 4 states"})

t = pd.DataFrame(rows)
t.to_csv(OUT / f"geometry_{PREFIX}.tsv", sep="\t", index=False, float_format="%.5f")
data = t[~t["mark"].isin(["axis", "printed counts", "no interval drawn"]) & t["row"].ne("ROWCOUNT")]
summ = data.groupby("panel").agg(marks_checked=("off", lambda s: int(s.notna().sum())), marks_off=("off", lambda s: int((s == True).sum())))
summ.to_csv(OUT / f"geometry_{PREFIX}_summary.tsv", sep="\t")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(summ)
print(t[t.off == True][["panel", "row", "label_on_panel", "mark", "drawn", "recomputed", "diff", "diff_pt", "colour"]].to_string())

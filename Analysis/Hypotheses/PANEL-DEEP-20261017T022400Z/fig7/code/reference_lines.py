"""Reference lines (zero lines, thresholds) and legend glyphs: invert each vertical reference line through the
panel's axis transform and list legend keys with and without a glyph. Usage: python reference_lines.py <out_dir> <main|calib>"""
import sys, json, importlib.util
import pandas as pd
OUT, PRE = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location("gc", __file__.replace("reference_lines.py", "geometry_check.py"))
src = open(spec.origin).read().split("# ---------------- expected values")[0]
ns = {"__name__": "x"}; sys_argv = sys.argv; exec(src.replace("OUT = Path(sys.argv[1])", f"OUT = Path('{OUT}')").replace("PREFIX = sys.argv[2]", f"PREFIX = '{PRE}'"), ns)
rows = []
for p, expect in {"fig7B_gain_over_task_native_baseline": [0.0], "fig7E_verdict_vocabulary_and_refusal_thresholds": [0.0, 1.0],
                  "fig7F_liver_encoder_against_the_field": [0.0, 0.0]}.items():
    g = ns["load"](p); axs = ns["axes"](g)
    vl = [l for l in g["lines"] if abs(l["x0"] - l["x1"]) < 0.01 and l["bottom"] - l["top"] > 10]
    for l in vl:
        ax = min(axs, key=lambda a: 0 if a["x0"] - 3 <= l["x0"] <= a["x1"] + 3 else 1e9)
        v = ax["a"] + ax["b"] * l["x0"]
        rows.append({"panel": p, "x": l["x0"], "value": round(v, 5), "stroke": l["stroke"], "dash": l.get("dash", ""), "lw": l["lw"],
                     "nearest_expected": min(expect, key=lambda e: abs(e - v)), "off_pt": abs(v - min(expect, key=lambda e: abs(e - v))) / abs(ax["b"])})
t = pd.DataFrame(rows); t.to_csv(f"{OUT}/reference_lines_{PRE}.tsv", sep="\t", index=False)
print(t.to_string())
# legend glyphs: words in the legend line, and glyphs (curves/rects/lines) on the same line
for p, ylo in {"fig7B_gain_over_task_native_baseline": 180, "fig7C_score_beside_training_exposure": 192, "fig7D_measured_bounds_and_join_structure": 222,
               "fig7F_liver_encoder_against_the_field": 205, "fig7A_framework_schematic": 100}.items():
    g = ns["load"](p)
    words = [w for w in g["words"] if w["top"] > ylo]
    glyphs = [("curve", (c["x0"] + c["x1"]) / 2, c["stroke"], c["fill"] if c["filled"] else "none") for c in g["curves"] if c["top"] > ylo] + \
             [("rect", (r["x0"] + r["x1"]) / 2, r["stroke"], r["fill"]) for r in g["rects"] if r["top"] > ylo and r["x1"] - r["x0"] < 10]
    print(p, "legend words:", " ".join(w["text"] for w in sorted(words, key=lambda w: w["x0"])))
    for w in sorted(words, key=lambda w: w["x0"]):
        pass
    print("   glyphs:", [(k, round(x, 1), s, f) for k, x, s, f in sorted(glyphs, key=lambda z: z[1])])

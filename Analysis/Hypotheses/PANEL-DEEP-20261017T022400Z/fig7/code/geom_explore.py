import json, sys, re
g = json.load(open(sys.argv[1]))
num = re.compile(r"^[−\-–]?\d+(\.\d+)?$")
print("size", g["width"], g["height"])
print("NUMERIC WORDS")
for w in g["words"]:
    if num.match(w["text"]):
        print(f"  {w['text']:>6} xc={(w['x0']+w['x1'])/2:.2f} top={w['top']:.2f}")
print("LINES (non-degenerate)")
for l in g["lines"]:
    L = ((l["x1"]-l["x0"])**2 + (l["bottom"]-l["top"])**2) ** .5
    print(f"  x {l['x0']:.2f}-{l['x1']:.2f} y {l['top']:.2f}-{l['bottom']:.2f} len {L:.2f} lw {l['lw']} {l['stroke']} {l.get('dash','')}")
print("CURVES")
for c in g["curves"]:
    print(f"  xc {(c['x0']+c['x1'])/2:.3f} yc {(c['top']+c['bottom'])/2:.2f} w {c['x1']-c['x0']:.2f} stroke {c['stroke']} fill {c['fill']} filled {c['filled']} npts {len(c.get('pts',[]))}")
print("RECTS")
for r in g["rects"]:
    print(f"  x {r['x0']:.2f}-{r['x1']:.2f} y {r['top']:.2f}-{r['bottom']:.2f} fill {r['fill']} filled {r['filled']} stroke {r['stroke']} stroked {r['stroked']}")

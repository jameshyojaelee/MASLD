"""SSX-B (PRESPEC_SSX_B.md): intronic features for the 24,115 release genes.
Gene span minus the union of exons of any gene on the same strand; genes whose remaining intronic region overlaps
another gene's span on the same strand are removed; intervals < 50 bp dropped.
Usage: python3 -I ssx_b_intron_gtf.py <gencode gtf.gz> <release gene list> <out dir>"""
import gzip, sys, os, re, bisect, json, collections
gtf, genelist, out = sys.argv[1:4]
assert not os.path.exists(out); os.makedirs(out)
release = set(l.strip() for l in open(genelist) if l.strip())
genes = {}; exons = collections.defaultdict(list)
gid_re = re.compile(r'gene_id "([^"]+)"')
with gzip.open(gtf, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.split("\t", 8)
        if c[2] not in ("gene", "exon"): continue
        gid = gid_re.search(c[8]).group(1); key = (c[0], c[6]); s, e = int(c[3]), int(c[4])
        if c[2] == "gene": genes[gid] = (c[0], c[6], s, e)
        else: exons[key].append((s, e))
def merge(iv):
    iv.sort(); m = []
    for s, e in iv:
        if m and s <= m[-1][1] + 1: m[-1][1] = max(m[-1][1], e)
        else: m.append([s, e])
    return m
mex = {k: merge(v) for k, v in exons.items()}
mex_starts = {k: [x[0] for x in v] for k, v in mex.items()}
spans = collections.defaultdict(list)
for gid, (ch, st, s, e) in genes.items(): spans[(ch, st)].append((s, e, gid))
for k in spans: spans[k].sort()
span_starts = {k: [x[0] for x in v] for k, v in spans.items()}
missing = sorted(release - set(genes)); n = collections.Counter(); rows = []
for gid in sorted(release & set(genes)):
    ch, st, s, e = genes[gid]; key = (ch, st)
    ex = mex.get(key, []); i = max(bisect.bisect_left(mex_starts.get(key, []), s) - 1, 0)
    pos = s; iv = []
    while i < len(ex) and ex[i][0] <= e:
        a, b = ex[i]
        if b >= pos and a > pos: iv.append((pos, min(a - 1, e)))
        pos = max(pos, b + 1); i += 1
    if pos <= e: iv.append((pos, e))
    iv = [(a, b) for a, b in iv if b - a + 1 >= 50]
    if not iv: n["no_intron_ge50"] += 1; continue
    sp = spans[key]; j = bisect.bisect_right(span_starts[key], e); clash = False
    for a, b in iv:
        for k2 in range(j - 1, -1, -1):
            gs, ge, g2 = sp[k2]
            if g2 == gid: continue
            if gs <= b and ge >= a: clash = True; break
            if gs < a - 3_000_000: break
        if clash: break
    if clash: n["other_gene_overlap"] += 1; continue
    n["kept"] += 1
    for a, b in iv: rows.append(f'{ch}\tSSX\tintron\t{a}\t{b}\t.\t{st}\t.\tgene_id "{gid}";\n')
with open(os.path.join(out, "intronic_features.gtf"), "w") as f: f.writelines(rows)
n["release"] = len(release); n["release_missing_in_gtf"] = len(missing); n["intervals"] = len(rows)
json.dump(n, open(os.path.join(out, "summary.json"), "w"), indent=1); print(dict(n))

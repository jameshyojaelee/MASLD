#!/usr/bin/env python3
"""
Pre-registered before/after comparison: COLOC rerun vs canonical.

WRITTEN BEFORE THE RERUN FINISHED, deliberately. The predictions live in
tmp_diag/PREREGISTRATION_coloc_rerun.md and were fixed before any rerun output
existed; this script fixes the ANALYSIS before the full data exists, so neither
the thresholds nor the subgroups can be chosen after seeing results.

Anti-rationalisation rules carried over from the pre-registration:
  1. Compute every check BEFORE writing any narrative.
  2. Report every failed prediction explicitly.
  3. Do not introduce a new subgroup/threshold/exclusion after seeing numbers.
  4. Gene-level, restricted to genes present in BOTH runs; genes present in only
     one are reported separately as coverage, never silently dropped.
  5. If the rerun is incomplete, report partial coverage rather than extrapolate.

P1  ABF PP.H4.abf reproduces: median |d| < 0.001 AND r > 0.99
P2  ABF crossings of 0.5 < 1% of the ABF gene set
P3  SuSiE PP.H4.susie r > 0.98 AND crossings < 5%
P4  every crossing gene is within 0.05 of 0.5 in at least one arm
P5  chr6 SuSiE coverage drops ~62 genes; no OTHER chromosome loses > 5
P6  tier-1/2 headline moves by < 10 genes net      [needs 07; run separately]
P7  ALL of the above, split by whether the gene's top SNP is PALINDROMIC.
    P7 exists because pooling the ~13% of genes whose evidence actually changed
    with the ~87% that could not change is what produced my earlier misleading
    "median |delta| is tiny" summary. The affected subgroup is reported alone.

PERFORMANCE NOTE: a first version recovered top-SNP alleles by scanning the 50
per-study sumstats (48 GB, ~500M rows in Python) and was killed at 900 s. Whether
a variant is palindromic is a property of the ALLELE PAIR, not of the study, so
the 1000G panel AF sidecars (55 MB total, one small file per chromosome) answer
it identically for ~900x less I/O. Row storage is likewise reduced to the three
fields the checks use, instead of 25-field dicts for ~2.5M rows.
"""
import csv, glob, gzip, os, statistics, sys
from collections import defaultdict

FM = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
CANON = os.path.join(FM, "results/susie_coloc")
RERUN = os.path.join(FM, "results/susie_coloc_rerun")
PANEL_AF = os.path.join(FM, "data/ld_ref/panel_af")
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def f(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def load_side(root):
    """(study, chr, ensembl) -> (pp_abf, pp_susie, top_snp). Only what P1-P7 use."""
    out = {}
    for p in glob.glob(os.path.join(root, "*", "susie_coloc_chr*.csv")):
        study = os.path.basename(os.path.dirname(p))
        chrom = os.path.basename(p).replace("susie_coloc_chr", "").replace(".csv", "")
        try:
            with open(p) as fh:
                for r in csv.DictReader(fh):
                    g = r.get("ensembl")
                    if g:
                        out[(study, chrom, g)] = (
                            f(r.get("PP.H4.abf")),
                            f(r.get("PP.H4.susie")),
                            (r.get("top_snp") or "").strip(),
                        )
        except Exception as e:
            print(f"  WARN unreadable {p}: {e}", file=sys.stderr)
    return out


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / (sxx * syy) ** 0.5


def palindromic_top_snps(keys, rerun):
    """Classify each gene's top SNP as palindromic / not, from panel sidecars."""
    need = defaultdict(set)                    # chr -> {pos}
    top = {}
    for k in keys:
        ts = rerun[k][2]
        if ":" not in ts:
            continue
        c, _, pos = ts.partition(":")
        if pos.isdigit():
            need[c].add(pos)
            top[k] = (c, pos)

    alleles = {}
    for c, positions in need.items():
        path = os.path.join(PANEL_AF, "eur", f"chr{c}.af.tsv.gz")
        if not os.path.exists(path):
            continue
        with gzip.open(path, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            try:
                ip, i1, i2 = hdr.index("position"), hdr.index("bim_a1"), hdr.index("bim_a2")
            except ValueError:
                continue
            for line in fh:
                p = line.split("\t", ip + 1)[ip] if ip == 0 else None
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(ip, i1, i2):
                    continue
                p = parts[ip]
                if p in positions:
                    alleles[(c, p)] = (parts[i1].upper(), parts[i2].upper())

    pal, nonpal, unknown = set(), set(), set()
    for k, cp in top.items():
        a = alleles.get(cp)
        if not a or a[0] not in COMPLEMENT or a[1] not in COMPLEMENT:
            unknown.add(k)
        elif COMPLEMENT[a[0]] == a[1]:
            pal.add(k)
        else:
            nonpal.add(k)
    return pal, nonpal, unknown


IDX = {"PP.H4.abf": 0, "PP.H4.susie": 1}


def arm_stats(keys, canon, rerun, col, label):
    i = IDX[col]
    xs, ys, cross, near = [], [], [], 0
    for k in keys:
        a, b = canon[k][i], rerun[k][i]
        if a is None or b is None:
            continue
        xs.append(a); ys.append(b)
        if (a > 0.5) != (b > 0.5):
            cross.append((k, a, b))
            if abs(a - 0.5) <= 0.05 or abs(b - 0.5) <= 0.05:
                near += 1
    if not xs:
        print(f"  {label}: no comparable genes")
        return None
    d = [abs(a - b) for a, b in zip(xs, ys)]
    r = pearson(xs, ys)
    print(f"  {label}: n={len(xs):6d}  median|d|={statistics.median(d):.3e}  "
          f"max|d|={max(d):.3e}  r={'NA' if r is None else f'{r:.6f}'}  "
          f"crossings={len(cross)} ({100*len(cross)/len(xs):.2f}%)  "
          f"within 0.05 of cut: {near}")
    return dict(n=len(xs), med=statistics.median(d), mx=max(d), r=r,
                cross=cross, near=near)


def main():
    n_rerun_files = len(glob.glob(os.path.join(RERUN, "*", "susie_coloc_chr*.csv")))
    print(f"rerun output files: {n_rerun_files} / 1100")
    if n_rerun_files < 1100:
        print("  *** RERUN INCOMPLETE -- partial coverage per rule 5, NOT extrapolated ***")

    canon, rerun = load_side(CANON), load_side(RERUN)
    ck, rk = set(canon), set(rerun)
    both = ck & rk
    print(f"\ncoverage: canonical {len(ck):,} | rerun {len(rk):,} | both {len(both):,}"
          f" | canonical-only {len(ck-rk):,} | rerun-only {len(rk-ck):,}")

    print("\n=== P1/P2  ABF arm (uses no LD; provably sign-invariant) ===")
    abf = arm_stats(both, canon, rerun, "PP.H4.abf", "ABF all")
    print("\n=== P3/P4  SuSiE arm ===")
    sus = arm_stats(both, canon, rerun, "PP.H4.susie", "SuSiE all")

    print("\n=== P5  per-chromosome SuSiE coverage change ===")
    def susie_n(side, keys):
        c = defaultdict(int)
        for k in keys:
            if side[k][1] is not None:
                c[k[1]] += 1
        return c
    cs, rs = susie_n(canon, ck), susie_n(rerun, rk)
    p5_fail = []
    for c in map(str, range(1, 23)):
        delta = rs.get(c, 0) - cs.get(c, 0)
        flag = ""
        if c == "6":
            flag = "  <- MHC excluded by design"
        elif delta < -5:
            flag = "  *** P5 FAIL: loses > 5"; p5_fail.append(c)
        print(f"  chr{c:<3} canonical {cs.get(c,0):5d}  rerun {rs.get(c,0):5d}  delta {delta:+5d}{flag}")

    print("\n=== P7  split by whether the gene's TOP SNP is palindromic ===")
    pal, nonpal, unk = palindromic_top_snps(both, rerun)
    print(f"  top-SNP class: palindromic {len(pal):,} | non-palindromic {len(nonpal):,} | unresolved {len(unk):,}")
    for name, subset in (("PALINDROMIC top SNP", pal), ("non-palindromic top SNP", nonpal)):
        if subset:
            print(f"  -- {name} --")
            arm_stats(subset, canon, rerun, "PP.H4.abf", "    ABF  ")
            arm_stats(subset, canon, rerun, "PP.H4.susie", "    SuSiE")

    print("\n=== VERDICTS (thresholds fixed in advance) ===")
    def verdict(tag, ok, detail):
        print(f"  {tag}: {'PASS' if ok else '*** FAIL'}  {detail}")
    if abf:
        verdict("P1", abf["med"] < 0.001 and (abf["r"] or 0) > 0.99,
                f"median|d|={abf['med']:.3e} (<0.001), r={abf['r']:.6f} (>0.99)")
        verdict("P2", len(abf["cross"]) < 0.01 * abf["n"],
                f"{len(abf['cross'])} crossings = {100*len(abf['cross'])/abf['n']:.2f}% (<1%)")
    if sus:
        verdict("P3", (sus["r"] or 0) > 0.98 and len(sus["cross"]) < 0.05 * sus["n"],
                f"r={sus['r']:.6f} (>0.98), {100*len(sus['cross'])/sus['n']:.2f}% crossings (<5%)")
        allc = (abf["cross"] if abf else []) + sus["cross"]
        far = [c for c in allc if abs(c[1]-0.5) > 0.05 and abs(c[2]-0.5) > 0.05]
        verdict("P4", not far, f"{len(far)} crossing genes >0.05 from the cut in BOTH arms")
    verdict("P5", not p5_fail, f"non-chr6 chromosomes losing >5: {p5_fail or 'none'}")
    print("  P6: requires 07_combine_susie_coloc.R on the rerun. NOTE: that script")
    print("      hardcodes COLOC_DIR=results/susie_coloc and writes gene_level_coloc.csv")
    print("      -- run a redirected COPY only; never the original against the rerun.")


if __name__ == "__main__":
    main()

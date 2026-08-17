#!/usr/bin/env python3
"""
Why did the rerun lose ~3% of SuSiE coverage and 24 headline genes?

Established already, and NOT the answer:
  * the eQTL SuSiE substrate      (+12 fits net; chr11 lost 211 SuSiE genes
                                   while its eQTL fits moved by -1)
  * the rebuilt non-EUR LD panels (EUR lost the most in absolute terms, -502,
                                   and EUR runs on PolyFun, which was untouched)

06_susie_coloc.R records its own reason in `method`:
    susie         the SuSiE arm ran and colocalised
    abf_fallback  an eQTL .rds EXISTED but the SuSiE arm did not complete --
                  LD load failed, too few matched/common variants, or a
                  non-converged fit on either side
    abf_only      no eQTL .rds for this gene at all

So every gene that had a SuSiE result in canonical and lost it must appear in the
rerun as abf_fallback or abf_only, and which one localises the cause:
    -> abf_only     means the eQTL fit went missing (substrate)
    -> abf_fallback means the fit was there but the arm bailed (LD / n / convergence)

For the abf_fallback group this also compares n_snps and n_match against the
canonical row for the SAME gene, because MIN_TRIPLE_SNPS = 50 gates on matched
variant counts and a drop there would be a concrete, checkable mechanism.
"""
import csv, glob, os, statistics
from collections import defaultdict, Counter

FM = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"


def load(root):
    d = {}
    for p in glob.glob(os.path.join(root, "*", "susie_coloc_chr*.csv")):
        s = os.path.basename(os.path.dirname(p))
        c = os.path.basename(p).replace("susie_coloc_chr", "").replace(".csv", "")
        with open(p) as fh:
            for r in csv.DictReader(fh):
                g = r.get("ensembl")
                if g:
                    d[(s, c, g)] = r
    return d


def fnum(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def has_susie(r):
    return fnum(r.get("PP.H4.susie")) is not None


canon = load(os.path.join(FM, "results/susie_coloc"))
rerun = load(os.path.join(FM, "results/susie_coloc_rerun"))
both = set(canon) & set(rerun)

lost = [k for k in both if has_susie(canon[k]) and not has_susie(rerun[k])]
gained = [k for k in both if not has_susie(canon[k]) and has_susie(rerun[k])]
kept = [k for k in both if has_susie(canon[k]) and has_susie(rerun[k])]

print(f"gene-GWAS pairs in both runs : {len(both):,}")
print(f"  SuSiE kept                 : {len(kept):,}")
print(f"  SuSiE LOST in rerun        : {len(lost):,}")
print(f"  SuSiE GAINED in rerun      : {len(gained):,}")
print(f"  net                        : {len(gained)-len(lost):+,}")

print("\n=== what does the rerun say the LOST ones did? (`method`) ===")
for m, n in Counter(rerun[k].get("method", "?") for k in lost).most_common():
    print(f"  {m:<16} {n:6d}  ({100*n/max(len(lost),1):5.1f}%)")

print("\n=== and what did canonical call them? ===")
for m, n in Counter(canon[k].get("method", "?") for k in lost).most_common():
    print(f"  {m:<16} {n:6d}")

fb = [k for k in lost if rerun[k].get("method") == "abf_fallback"]
ao = [k for k in lost if rerun[k].get("method") == "abf_only"]
print(f"\n  -> abf_only  = {len(ao):,}  (eQTL fit absent: substrate)")
print(f"  -> abf_fallback = {len(fb):,}  (fit present, arm bailed: LD / n / convergence)")

if fb:
    print("\n=== for the abf_fallback losses: did matched-variant counts move? ===")
    print("    MIN_TRIPLE_SNPS = 50 gates on these.")
    for col in ("n_snps", "n_match"):
        dc, dr, dd = [], [], []
        for k in fb:
            a, b = fnum(canon[k].get(col)), fnum(rerun[k].get(col))
            if a is not None and b is not None:
                dc.append(a); dr.append(b); dd.append(b - a)
        if dd:
            below = sum(1 for v in dr if v < 50)
            print(f"  {col:<10} canonical median {statistics.median(dc):8.0f} | "
                  f"rerun median {statistics.median(dr):8.0f} | "
                  f"median delta {statistics.median(dd):+8.0f} | "
                  f"rerun below 50: {below}")

    print("\n=== are the losses concentrated in particular studies? ===")
    per = Counter(k[0] for k in fb)
    for s, n in per.most_common(10):
        tot = sum(1 for k in both if k[0] == s and has_susie(canon[k]))
        print(f"  {s:<34} {n:5d} lost of {tot:5d} canonical susie ({100*n/max(tot,1):4.1f}%)")

print("\n=== headline impact: lost genes that were ABOVE PP.H4.susie 0.5 in canonical ===")
hi = [k for k in lost if (fnum(canon[k].get("PP.H4.susie")) or 0) > 0.5]
print(f"  {len(hi)} gene-GWAS pairs crossed out of the >0.5 set by losing the arm entirely")
if hi:
    print("  top 10 by canonical PP.H4.susie:")
    for k in sorted(hi, key=lambda k: -(fnum(canon[k].get("PP.H4.susie")) or 0))[:10]:
        print(f"    {canon[k].get('gene','?'):<14} {k[0]:<30} chr{k[1]:<3} "
              f"canon susie={fnum(canon[k].get('PP.H4.susie')):.3f} "
              f"rerun method={rerun[k].get('method')} "
              f"n_match {canon[k].get('n_match')} -> {rerun[k].get('n_match')}")

#!/usr/bin/env python3
"""
Structural verification of the COLOC rerun outputs.

Everything here is checkable NOW, on whatever has landed -- only the final
"1100/1100" assertion has to wait for the last tasks. Separating the two matters:
if 1,092 files share a defect, that should surface now, not after the run ends.

Checks, per output file:
  1. header present and exactly the 25 expected columns, in the expected order
  2. every data row has the same field count as the header (no truncation)
  3. at least one data row
  4. ensembl IDs unique within the file (a duplicate implies a bad resume/merge --
     the checkpoint merge path dedups on ensembl, so a dup means that failed)
  5. PP.H0..H4.abf present and summing to ~1 per row (coloc posteriors are a
     partition; a row that does not sum to 1 is corrupt, not merely unusual)
  6. no negative or >1 posteriors
  7. lambda_s_locus, where present, below the LAMBDA_S_HIGH = 0.20 gate is
     reported as a count, not asserted -- values above it are legitimate, this
     is the QC flag itself

And across files:
  8. (study, chr) coverage matrix vs the expected 50 x 22
  9. gene counts vs the SAME study+chr in canonical (apples-to-apples; a rerun
     count BELOW canonical is the direction that would indicate lost work, since
     palindrome retention can only add testable genes)
"""
import csv, glob, os, sys
from collections import defaultdict

FM = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
RERUN = os.path.join(FM, "results/susie_coloc_rerun")
CANON = os.path.join(FM, "results/susie_coloc")

EXPECTED = ["gene", "ensembl", "chr", "gwas_name",
            "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf",
            "PP.H3.susie", "PP.H4.susie", "n_cs_pairs", "n_snps", "n_match",
            "n_flip", "n_unresolved", "n_ambiguous", "n_palindromic_kept",
            "n_palindromic_af_tested", "nonpal_pairmiss",
            "n_palindromic_af_dropped", "method", "lambda_s_locus",
            "top_snp", "top_snp_PP"]
PPCOLS = ["PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf"]


def num(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def main():
    files = sorted(glob.glob(os.path.join(RERUN, "*", "susie_coloc_chr*.csv")))
    print(f"outputs present: {len(files)} / 1100")
    if len(files) < 1100:
        print(f"  (verifying the {len(files)} that exist; the 1100/1100 assertion waits)\n")

    bad_header, bad_width, empty, dup, bad_sum, out_of_range = [], [], [], [], [], []
    total_rows = 0
    lam_above = 0
    lam_n = 0
    counts = {}

    for p in files:
        study = os.path.basename(os.path.dirname(p))
        chrom = os.path.basename(p).replace("susie_coloc_chr", "").replace(".csv", "")
        with open(p) as fh:
            rd = csv.reader(fh)
            try:
                hdr = next(rd)
            except StopIteration:
                empty.append(p); continue
            if hdr != EXPECTED:
                bad_header.append((p, len(hdr)))
                continue
            idx = {c: i for i, c in enumerate(hdr)}
            seen, n = set(), 0
            for row in rd:
                n += 1
                if len(row) != len(hdr):
                    bad_width.append((p, n, len(row))); continue
                g = row[idx["ensembl"]]
                if g in seen:
                    dup.append((p, g))
                seen.add(g)
                pps = [num(row[idx[c]]) for c in PPCOLS]
                if all(v is not None for v in pps):
                    s = sum(pps)
                    if abs(s - 1.0) > 1e-3:
                        bad_sum.append((p, g, s))
                    if any(v < -1e-9 or v > 1 + 1e-9 for v in pps):
                        out_of_range.append((p, g))
                lv = num(row[idx["lambda_s_locus"]])
                if lv is not None:
                    lam_n += 1
                    if lv > 0.20:
                        lam_above += 1
            total_rows += n
            counts[(study, chrom)] = n
            if n == 0:
                empty.append(p)

    print("=== per-file structural checks ===")
    def rep(name, lst, ok_msg):
        if lst:
            print(f"  *** {name}: {len(lst)}")
            for x in lst[:5]:
                print(f"        {x}")
        else:
            print(f"  {name}: {ok_msg}")
    rep("wrong header", bad_header, "0 — all files have the 25 expected columns in order")
    rep("truncated rows", bad_width, "0 — every row matches its header width")
    rep("empty files", empty, "0")
    rep("duplicate ensembl within a file", dup, "0 — checkpoint/merge dedup held")
    rep("PP.H0..H4 not summing to 1", bad_sum, "0 — posteriors form a valid partition")
    rep("posteriors outside [0,1]", out_of_range, "0")
    print(f"  total gene rows: {total_rows:,}")
    print(f"  lambda_s values: {lam_n:,}   above the 0.20 gate: {lam_above}")

    print("\n=== coverage matrix ===")
    studies = sorted({s for s, _ in counts})
    print(f"  studies with >=1 output: {len(studies)} / 50")
    per_chr = defaultdict(int)
    for (_, c) in counts:
        per_chr[c] += 1
    short = [c for c in map(str, range(1, 23)) if per_chr.get(c, 0) < 50]
    print(f"  chromosomes at 50/50: {22 - len(short)} / 22"
          + (f"   short: {', '.join('chr'+c+f'({per_chr.get(c,0)})' for c in short)}" if short else ""))

    print("\n=== gene counts vs the SAME study+chr in canonical ===")
    lower, equal, higher, nocanon = 0, 0, 0, 0
    worst = []
    for (study, chrom), n in counts.items():
        cp = os.path.join(CANON, study, f"susie_coloc_chr{chrom}.csv")
        if not os.path.exists(cp):
            nocanon += 1; continue
        with open(cp) as fh:
            cn = sum(1 for _ in fh) - 1
        if n < cn:
            lower += 1; worst.append((study, chrom, n, cn, n - cn))
        elif n == cn:
            equal += 1
        else:
            higher += 1
    print(f"  rerun == canonical : {equal}")
    print(f"  rerun  > canonical : {higher}   (expected direction: palindrome retention adds testable genes)")
    print(f"  rerun  < canonical : {lower}    (the direction that would indicate lost work)")
    print(f"  no canonical counterpart: {nocanon}")
    if worst:
        worst.sort(key=lambda r: r[4])
        print("  largest shortfalls:")
        for s, c, n, cn, d in worst[:10]:
            print(f"        {s} chr{c}: rerun {n} vs canonical {cn}  ({d:+d})")

    fatal = bad_header or bad_width or empty or dup or bad_sum or out_of_range
    print("\n=== VERDICT ===")
    print("  STRUCTURAL: " + ("*** PROBLEMS FOUND (see above)" if fatal
                              else "clean — every landed output is well-formed"))
    print(f"  COMPLETENESS: {len(files)}/1100 "
          + ("— COMPLETE" if len(files) >= 1100 else "— still filling, assertion deferred"))


if __name__ == "__main__":
    main()

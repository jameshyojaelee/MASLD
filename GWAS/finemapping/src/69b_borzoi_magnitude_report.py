#!/usr/bin/env python3
"""69b: honest summary of the Axis-4 Borzoi magnitude nominations. Read-only over
borzoi_magnitude_nominations.tsv + borzoi_magnitude_candidates.tsv + the eQTL-present
reference. Emits the headline counts AND the best-of-K / lncRNA-volatility caveats and
stricter strata, so the nomination is not over-sold."""
import csv, os, numpy as np
FM = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results"
S = os.path.join(FM, "seqfunc")
nom = list(csv.DictReader(open(os.path.join(S, "borzoi_magnitude_nominations.tsv")), delimiter="\t"))
ref = [abs(float(r["borzoi_logsed_liver"])) for r in
       csv.DictReader(open(os.path.join(S, "borzoi_eqtl_logsed_scores.tsv")), delimiter="\t")
       if r.get("is_anchor") != "TRUE" and r.get("borzoi_logsed_liver") not in ("", "nan", "NaN")]
ref = np.array(sorted(ref)); p80, p90, p95 = np.percentile(ref, [80, 90, 95])
SORT1 = 0.051150  # anchor magnitude from 64b (a bona-fide positive-control effect)

def fnum(r, k):
    v = r.get(k, "")
    try: return float(v)
    except (ValueError, TypeError): return float("nan")

scored = [r for r in nom if r.get("candidate_genes_scored") not in ("", "0", None)]
best = np.array([fnum(r, "borzoi_abs_logsed") for r in scored])
best = best[np.isfinite(best)]
n_p80 = sum(fnum(r, "borzoi_abs_logsed") > p80 for r in scored)
n_p90 = sum(fnum(r, "borzoi_abs_logsed") > p90 for r in scored)
n_sort1 = sum(fnum(r, "borzoi_abs_logsed") >= SORT1 for r in scored)
bt = {}
for r in scored: bt[r.get("best_gene_biotype", "?")] = bt.get(r.get("best_gene_biotype", "?"), 0) + 1
agree = {}
for r in scored: agree[r.get("fold_agreement", "?")] = agree.get(r.get("fold_agreement", "?"), 0) + 1
# apples-to-apples: substrate (pre-specified single) gene percentile, best-of-K-free
subp = [fnum(r, "substrate_gene_percentile") for r in scored]
subp = [x for x in subp if np.isfinite(x)]
n_sub_p80 = sum(x > 80 for x in subp)

print("================ AXIS-4 BORZOI MAGNITUDE NOMINATION SUMMARY ================")
print(f"eQTL-absent signal leads (input)        : {len(nom)}")
print(f"loci with >=1 scorable candidate gene    : {len(scored)}")
print(f"reference (eQTL-present) n / p80/p90/p95  : {len(ref)} / {p80:.4f} / {p90:.4f} / {p95:.4f}")
print(f"reference max |logSED| / SORT1 anchor     : {ref.max():.4f} / {SORT1:.4f}")
print("-- HEADLINE (as specified: best-gene |logSED| vs eQTL-present p80) --")
print(f"loci > 80th magnitude pct (best-of-K)     : {n_p80}/{len(scored)} ({100*n_p80/len(scored):.0f}%)")
print("-- STRICTER STRATA (guard the best-of-K + lncRNA-volatility inflation) --")
print(f"loci > 90th pct                           : {n_p90}/{len(scored)}")
print(f"loci >= SORT1 positive-control magnitude   : {n_sort1}/{len(scored)}")
print(f"best |logSED| distribution: median {np.median(best):.4f}  p75 {np.percentile(best,75):.4f}  max {best.max():.4f}")
print(f"best-gene biotype                         : {bt}")
print(f"best-gene fold_sign_agreement (0-4)       : {agree}")
print("-- APPLES-TO-APPLES (pre-specified substrate gene, no best-of-K selection) --")
print(f"substrate-gene percentile available        : {len(subp)} loci")
print(f"substrate-gene > 80th pct                  : {n_sub_p80}/{len(subp)}")
print("-- TOP 12 NOMINATIONS by |logSED| magnitude --")
top = sorted(scored, key=lambda r: fnum(r, "borzoi_abs_logsed"), reverse=True)[:12]
print(f"{'locus':16} {'variant':22} {'best_gene':16} {'biotype':14} {'|logSED|':>9} {'pct':>5} {'agree':>5} {'pip':>6}")
for r in top:
    print(f"{r['locus_id']:16} {r['variant_id_hg19']:22} {r['best_gene'][:15]:16} "
          f"{r['best_gene_biotype']:14} {fnum(r,'borzoi_abs_logsed'):9.4f} "
          f"{fnum(r,'magnitude_percentile'):5.0f} {r['fold_agreement']:>5} {fnum(r,'lead_pip'):6.2f}")
print("==========================================================================")
# machine-readable TSV for the parent
with open(os.path.join(S, "borzoi_magnitude_report.txt"), "w") as fh:
    fh.write(f"n_leads\t{len(nom)}\nn_scored\t{len(scored)}\nref_n\t{len(ref)}\n"
             f"ref_p80\t{p80:.6f}\nref_p90\t{p90:.6f}\nn_exceed_p80\t{n_p80}\n"
             f"n_exceed_p90\t{n_p90}\nn_ge_sort1\t{n_sort1}\n"
             f"biotype\t{bt}\nfold_agreement\t{agree}\n"
             f"n_substrate_gene_percentiles\t{len(subp)}\nn_substrate_gt_p80\t{n_sub_p80}\n")

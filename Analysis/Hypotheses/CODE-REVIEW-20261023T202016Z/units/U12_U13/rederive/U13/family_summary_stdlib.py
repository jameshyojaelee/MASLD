"""Summarize 114-program family runs vs AN4 with the stdlib only (pandas segfaults reading one run file on the login node)."""
import csv, sys
from pathlib import Path
U = Path(__file__).resolve().parent
AN4 = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN4/20261021T124414Z")
STROMAL, DUCT = "hotspot_hepatocytes_f05c535ae5bbc0b9", "hotspot_hepatocytes_48f39dd4d817a10e"
def rows(p): return list(csv.DictReader(open(p), delimiter="\t"))
an = {}
for r in rows(AN4 / "combined/family_all_programs.tsv"):
    if r["testable"] == "True": an[(r["dataset"], r["program_id"])] = r
summ, flips = [], []
runs = ["gse_s42_reg", "gse_s42_rev", "gse_s20261023", "gse_s31337", "vu_s1042_reg", "vu_s1042_rev", "vu_s20261023", "vu_s31337"]
calls_by_run = {}
for t in runs:
    rr = [r for r in rows(U / "out" / t / "family_results.tsv") if r["testable"] == "True"]
    ds = rr[0]["dataset"]
    called = {r["program_id"] for r in rr if float(r["q_bh_family"]) < 0.05}
    an_called = {pid for (d, pid), a in an.items() if d == ds and float(a["residual_padj_family"]) < 0.05}
    calls_by_run[t] = called
    ps = sorted(float(r["p"]) for r in rr); n = len(ps); c = len(called)
    rank = {r["program_id"]: int(float(r["rank_excess"])) for r in rr}
    summ.append(dict(run=t, dataset=ds, n=n, calls=c, an4_calls=len(an_called), gained=len(called - an_called),
                     lost=len(an_called - called), sha_match=sum(r["sha_match"] == "True" for r in rr),
                     max_abs_dp_vs_an4=round(max(abs(float(r["p"]) - float(an[(ds, r["program_id"])]["residual_pvalue"])) for r in rr), 4),
                     p_last_called=ps[c - 1], cut_last=round(0.05 * c / n, 6), p_next=ps[c], cut_next=round(0.05 * (c + 1) / n, 6),
                     stromal_rank=rank[STROMAL], ductular_rank=rank[DUCT]))
    for pid in sorted(called ^ an_called):
        r = next(x for x in rr if x["program_id"] == pid)
        flips.append(dict(run=t, dataset=ds, program_id=pid, program_name=r["program_name"], p=r["p"], q=r["q_bh_family"],
                          p_an4=an[(ds, pid)]["residual_pvalue"], q_an4=an[(ds, pid)]["residual_padj_family"],
                          called_here=pid in called))
for name, data in [("family_runs_summary.tsv", summ), ("family_call_flips_vs_an4.tsv", flips)]:
    with open(U / name, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(data[0].keys()) if data else ["run"], delimiter="\t"); w.writeheader(); w.writerows(data)
for d in summ: print(d)
for f in flips: print(f)
for ds in ("gse", "vu"):
    sets = [calls_by_run[t] for t in runs if t.startswith(ds)]
    print(ds, "called in all 4 runs:", len(set.intersection(*sets)), "in any:", len(set.union(*sets)))

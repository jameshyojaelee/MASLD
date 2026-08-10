#!/usr/bin/env python3
"""Source table for the DEG-threshold sensitivity supplement.

The paper's headline — that colocalized genes are largely not differentially
expressed — depends on how "differentially expressed" is defined. A more
permissive gate yields more DEGs, more overlap, and a lower non-overlap
percentage. Rather than defend one threshold, we publish the whole ladder.

Emits, for a sweep of |logFC| floors paired with padj<0.05, plus the canonical
TREAT interval-null gate:
  n_DEG, overlap with the jointly-testable SuSiE set, % non-DE,
  Fisher OR and p against the jointly-testable background.

Read-only. Writes one TSV under figures/supplementary/.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from scipy import stats

PROJ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEG = PROJ / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"
SUSIE = (PROJ / "results/remediation/bg001/bg001-fragment-v211-gencode49-20260807T195243Z"
         / "frozen_sets/finite_susie_447.tsv")
OUT = PROJ / "figures/supplementary/figS_deg_threshold_sensitivity"

FLOORS = [0.00, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.75, 1.00]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(DEG)
    d["g"] = d["gene"].astype(str).str.split(".").str[0]
    su = pd.read_csv(SUSIE, sep="\t")
    su["g"] = su["gene"].astype(str).str.split(".").str[0]

    uni = set(d["g"])
    jt = set(su["g"]) & uni
    assert len(jt) == 447, f"expected 447 jointly testable SuSiE genes, got {len(jt)}"

    def score(label, gate, mask):
        deg = set(d.loc[mask, "g"])
        a = len(jt & deg); b = len(jt) - a
        c = len(deg - jt); e = len(uni) - len(jt) - c
        orr, p = stats.fisher_exact([[a, b], [c, e]])
        return dict(label=label, gate=gate, lfc_floor=None, n_deg=len(deg),
                    overlap=a, pct_non_de=100 * (1 - a / len(jt)),
                    fisher_or=orr, fisher_p=p)

    rows = []
    for t in FLOORS:
        m = (d["padj"] < 0.05) & (d["logFC"].abs() > t) if t > 0 else (d["padj"] < 0.05)
        r = score(f"padj<0.05, |logFC|>{t:.2f}" if t > 0 else "padj<0.05 (no floor)",
                  "conventional", m)
        r["lfc_floor"] = t
        rows.append(r)
    r = score("TREAT FDR<0.05 @ lfc=0.25", "treat", d["treat_fdr"] < 0.05)
    r["lfc_floor"] = 0.25
    rows.append(r)

    out = pd.DataFrame(rows)
    out["n_jointly_testable"] = len(jt)
    out.to_csv(OUT / "figS_deg_threshold_sensitivity_source.tsv", sep="\t", index=False)

    tre = out.loc[out["gate"] == "treat"].iloc[0]
    conv = out[out["gate"] == "conventional"].copy()
    conv["gap"] = (conv["n_deg"] - tre["n_deg"]).abs()
    match = conv.nsmallest(1, "gap").iloc[0]
    print(out.to_string(index=False))
    print(f"\nTREAT n={tre.n_deg}  pct_non_de={tre.pct_non_de:.2f}  OR={tre.fisher_or:.3f}")
    print(f"closest conventional: {match.label}  n={match.n_deg}  "
          f"pct_non_de={match.pct_non_de:.2f}  OR={match.fisher_or:.3f}")
    print(f"\nWROTE {OUT}")


if __name__ == "__main__":
    main()

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

The Fisher background is the 14,931-gene JOINTLY TESTABLE universe -- genes
tested in both the bulk contrast and colocalization -- not all 27,638 bulk-tested
genes. A gene never tested for colocalization cannot be counted as evidence
against overlap; including such genes inflates the double-negative cell and
therefore inflates every odds ratio. The earlier version of this script used
`uni = set(d["g"])` (all 27,638), which is an asymmetric denominator: it
reported the TREAT gate at OR 1.106 (enrichment) where the sealed release
artifact `manuscript_release/2026-07-15-r2/orthogonality_audit.tsv` reports
0.8895 (depletion) for the same genes. The assertion below pins this script to
that artifact so the two can no longer disagree.

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
# Supplies `joint_testable`, which the DEG table does not carry. Conversely the
# release table carries only TREAT statistics, not conventional padj, so the
# floor sweep needs both files joined on version-stripped Ensembl ID.
EVIDENCE = PROJ / "RNA-seq/results/manuscript_release/2026-07-15-r2/evidence_class_table.tsv"
OUT = PROJ / "figures/supplementary/figS_deg_threshold_sensitivity"

FLOORS = [0.00, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.75, 1.00]

# Ground truth from orthogonality_audit.tsv, row (all, susie, PP4>=0.5).
# The integer cells are the claim and must match exactly. The odds ratio is a
# derived conditional-MLE whose last digits depend on the solver -- scipy and R
# disagree around 5e-6 on an identical table -- so it gets a tolerance.
AUDIT_TREAT_CELLS = dict(overlap=34, n_deg_joint=1261, n_genetic=447,
                         n_background=14931)
AUDIT_TREAT_STATS = dict(fisher_or=0.889461813691099, fisher_p=0.603934955412512)
STAT_TOL = 1e-4


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(DEG)
    d["g"] = d["gene"].astype(str).str.split(".").str[0]
    su = pd.read_csv(SUSIE, sep="\t")
    su["g"] = su["gene"].astype(str).str.split(".").str[0]

    ev = pd.read_csv(EVIDENCE, sep="\t", usecols=["ensembl_bulk", "joint_testable"])
    ev["g"] = ev["ensembl_bulk"].astype(str).str.split(".").str[0]
    # pandas coerces the TSV's TRUE/FALSE to bool, so str(True) == "True", not
    # "TRUE". Normalise case rather than matching the on-disk spelling.
    ev["jt"] = ev["joint_testable"].astype(str).str.upper() == "TRUE"

    # The Fisher background: genes tested in BOTH the bulk contrast and coloc.
    uni = set(ev.loc[ev["jt"], "g"]) & set(d["g"])
    assert len(uni) == 14931, f"expected 14931 jointly testable genes, got {len(uni)}"

    jt = set(su["g"]) & uni
    assert len(jt) == 447, f"expected 447 jointly testable SuSiE genes, got {len(jt)}"

    def score(label, gate, mask):
        deg_all = set(d.loc[mask, "g"])
        deg = deg_all & uni                       # the test runs inside the universe
        a = len(jt & deg); b = len(jt) - a
        c = len(deg - jt); e = len(uni) - len(jt) - c
        orr, p = stats.fisher_exact([[a, b], [c, e]])
        return dict(label=label, gate=gate, lfc_floor=None,
                    n_deg_all=len(deg_all), n_deg_joint=len(deg),
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
    out["n_genetic"] = len(jt)
    out["n_background"] = len(uni)

    # Pin to the sealed release artifact. If this fails, the two disagree again.
    tre = out.loc[out["gate"] == "treat"].iloc[0]
    for k, want in AUDIT_TREAT_CELLS.items():
        assert int(tre[k]) == want, (
            f"TREAT row {k}={int(tre[k])} does not reproduce "
            f"orthogonality_audit.tsv ({want}). The Fisher background has "
            f"drifted from the release.")
    for k, want in AUDIT_TREAT_STATS.items():
        assert abs(float(tre[k]) - want) <= STAT_TOL, (
            f"TREAT row {k}={float(tre[k]):.9f} differs from "
            f"orthogonality_audit.tsv ({want:.9f}) by more than {STAT_TOL}.")

    out.to_csv(OUT / "figS_deg_threshold_sensitivity_source.tsv", sep="\t", index=False)

    conv = out[out["gate"] == "conventional"].copy()
    conv["gap"] = (conv["n_deg_joint"] - tre["n_deg_joint"]).abs()
    match = conv.nsmallest(1, "gap").iloc[0]
    print(out.to_string(index=False))
    print(f"\nbackground = {len(uni)} jointly testable; genetic set = {len(jt)}")
    print(f"TREAT n_joint={tre.n_deg_joint}  pct_non_de={tre.pct_non_de:.2f}  "
          f"OR={tre.fisher_or:.4f}  [reproduces orthogonality_audit.tsv]")
    print(f"closest conventional: {match.label}  n_joint={match.n_deg_joint}  "
          f"pct_non_de={match.pct_non_de:.2f}  OR={match.fisher_or:.4f}")
    print(f"\nWROTE {OUT}")


if __name__ == "__main__":
    main()

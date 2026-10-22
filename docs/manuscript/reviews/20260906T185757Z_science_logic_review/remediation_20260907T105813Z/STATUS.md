# Remediation status, 2026-09-07

## Done

| review takeaway | status |
|---|---|
| **1. Headline is the independence null** | Done. The chance expectation (92.9%) now sits beside the observed 95.4% in `RESULTS.md`, `01_intro.md`, `fig5_discussion.md` and the rebuilt Figure 1B. "Organizing observation" and "biologically informative" are gone. A failable replacement test was run and is reported. |
| **2. Retired July numbers** | Done. OR 2.1/9.1, gsMap, drug benchmark and "75 of 473" deleted with an audit comment recording what was removed and why. Figure 1B rebuilt on the adopted release as a candidate. Figure 1C has the same defect and is **blocked** (needs stratum-scope classification). |
| **3. F0 mixes controls and disease** | **Done, and larger than expected.** 31 control / 73 disease, nested in cohort. F0→F1 goes **1,471 → 7** on disease-only F0, and **1,471 → 258 at unchanged n** with a control-status term, so it is the confound and not the sample size. Retires "fibrosis remodeling expands early" and the GNMT/MAT1A F0→F1 effects. Arm A retained as primary per the PI, quoted with arm B beside it. Also fixed the 257,070 NA confidence intervals. |
| **4. 117 programs ≈ random on fibrosis** | Investigated and partly fixed. 37 of 117 programs are mislabelled by lineage; the cause is measured ambient RNA (0.432 vs 0.061 median burden, rank-biserial 0.998). Panels move to content labels. The 113/27 vs 114/26 conflict is reconciled by recomputation. The null itself stands and is now reported with its construction, band and axes. |
| **5. Abstract has no number** | Done. Coding-mass effect size, its robustness, the overlap-versus-chance sentence and the class-test null are in. Placeholders removed. Lead position left open pending the PI's headline decision. |
| **6. Worked examples are known biology** | Partly. GNMT retained; the ECM/BICC1 pair is now presented as an instance of a measured registry-wide property rather than a curiosity. |
| **8. Title's third clause unsupported** | Done. `draft/fig7.md` and a Figure 7 legend section written from verified producing files; seven-figure consistency fixed across STATUS, README and legends. |
| **9. Hedge density** | Partly. Rewritten passages state one claim with its number. A full pass is pending the headline decision. |

## Analyses run

1. **Class-difference test** — colocalized vs disease-state genes are **not distinguishable** on cell-type expression or protein response at this n. Prespecified, sealed, hash-verified. Found and corrected an anti-conservative null (1.59x too narrow) and a circular protein comparison that would have produced a confident wrong claim.
2. **Coding-mass robustness** — all three owed checks plus LD-block clustering **pass**; 62 refits, ratio 4.99-9.58, every p <= 8e-4. Dropping the study that carries 95 of 208 direct credible sets **raises** the ratio.
3. **All-117 label-vs-content audit** — 37 mislabelled, ambient-RNA mechanism measured, gold-standard validated.
4. **7F donor weighting** — 2 of 12 model blocks change sign; a 2-donor, 149-cell study carries a seventh of the metric.
5. **F0 three-arm rerun** — F0→F1 is 99.5% control-versus-disease; two published claims retired; the stage-extension confidence-interval defect fixed.

## Pending

- **GSE193066**, the sealed holdout. The PI approved opening it **once** for two prespecified tests (within-person instrument change; coordination-community replication). Both prespecifications must be written and hashed before the cohort is read. Not started.
- **Figure 1C** rebuild — blocked on the direct-disease vs liver-enzyme stratum-scope classification, which `HIGH_IMPACT_PATTERNS.md` row 12 records as not derivable from the adopted registry alone.
- **Headline, title and abstract lead** — deferred to the PI.
- Hedge-density pass over the full Discussion; SAF aspect test; continuum panels 4G/4H demotion.

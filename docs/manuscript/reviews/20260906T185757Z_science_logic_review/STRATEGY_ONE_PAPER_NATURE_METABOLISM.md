# One-paper strategy for Nature Metabolism (2026-09-06)

Supersedes the two-paper recommendation in REVIEW.md. The PI decided: one paper, target Nature Metabolism or higher.

## Thesis (one sentence)

Inherited MASLD risk and established-disease remodeling are separable layers of human liver biology: direct-disease risk concentrates in protein-altering variants of metabolic genes and is largely invisible to liver expression, whereas established disease is one diffuse, linear, fibro-inflammatory severity axis that saturates the reproducibility ceiling of histology, and the two layers meet in a hepatocyte metabolic-loss state.

## Why this thesis and not the current one

- It is positive, quantitative, and every leg already has a producing file.
- It is metabolic, which is what the journal is for.
- It explains the null that the current draft uses as its headline (chance overlap) instead of presenting it.
- It absorbs Figure 7: pretrained models add nothing because the tissue signal is diffuse and linear, and the exposure audit is why that conclusion can be trusted. Figure 7 stops being a benchmark and becomes the measurement-limit leg of the biology.
- It positions against Kamzolas (same journal, 2026) as complementary: they ordered participants along a continuum; this paper shows what that continuum is made of (one diffuse axis, no system-specific programs after competitive nulls), what bounds it (repeat-biopsy sigma_e 0.47 to 0.64 stages), and what inherited risk does instead.

## Three legs, each with the number that carries it

| leg | number | producing file | status |
|---|---|---|---|
| direct-disease risk is coding-heavy | 22.4% vs 3.4% protein-altering PIP mass; ratio 6.55; BH q 6e-4; survives effect-size matching | RESULTS.md:922-995 | done; owed: leave-one-study-out, LD-block clustering, window sensitivity, metabolic-gene annotation of the coding mass |
| established disease is one diffuse linear axis | 14,651 of 23,231 transcripts continuum-associated after stage; no system continuum-specific under the connected null; instrument transfers to 585 external participants at 0.68 to 0.79 of ceiling; nonlinear and pretrained models add 0.000 ± 0.067 | RESULTS.md:516-802; MODEL_RESULTS_SUMMARY.md:59-62 | done; owed: disease-only F0 rerun; GSE193066 within-person change |
| the layers meet in hepatocyte metabolic loss | GNMT, MAT1A, CYP2C19 lower in RNA, stage, and protein; hepatocyte-metabolism NMF axis falls -0.37 to -1.70 SD F1 to F4; the two gene sets intersect at chance (20 vs 31.2) | RESULTS.md:74-99; HIGH_IMPACT row 24 | done; owed: a class-difference test that can fail (colocalized vs disease-state genes by cell-type expression and protein response), so "separable" is tested, not asserted |

## Figure plan (seven, within Nature Portfolio limits [Unverified: recheck the current author guide])

1. Cohorts and the two layers. Drop the observability counts from the main panel.
2. Inherited risk: coding mass at direct-disease loci versus enzyme loci; consequence and ancestry context; PNPLA3/TM6SF2/MTARC1 as the anchors, not GNMT.
3. Established disease: five-cohort DEGs, corrected adjacent stage contrasts, the hepatocyte metabolic-loss state with protein confirmation.
4. One diffuse axis: 4F focal programs, the competitive-null result as the panel, the vocabulary null (117 programs do not beat random on fibrosis) in the supplement; 4G/4H to supplement.
5. Measurement limits: severity instrument transfer to 585 participants, repeat-biopsy ceiling, exposure audit, models add nothing. This is the former Figure 7, rebuilt around the thesis; fix 7F by donor weighting or drop it.
6. Where the layers meet and where they do not: chance overlap, the failable class-difference test, tissue context for the classes (the best of current Figure 5).
7. Gene Catalog and portal: one worked entry per class; evidence states as a legend, not a panel.

Cut or supplement: the seven-state accounting as a headline; GeoMx/CosMx/GSE287826 untestability counts; ATAC beyond the four stellate peaks; lncRNA package; the continuum-translation therapeutic-hypothesis extension of 2026-09-06 (candidate, no headline, adds a fourth vocabulary); ligand-receptor and TF activity.

## Analyses to run before writing, in order

1. Rebuild the Discussion and Figure 1B on the adopted release; delete the July gsMap, drug-benchmark, 75-of-473 and OR 2.1/9.1 claims. Days.
2. Rerun adjacent fibrosis with disease-only F0 and NAS with per-cohort effects. Days.
3. Class-difference test with a null that can fail: do colocalized genes and disease-state genes differ in hepatocyte versus non-parenchymal expression (donor-level atlas), in protein response (PXD051911), and in coding architecture. Prespecify. One week.
4. Coding-mass robustness and annotation: leave-one-study-out, LD-block clustering, window sensitivity; fraction of coding mass in lipid-metabolism and hepatocyte-expressed genes. One to two weeks.
5. GSE193066 paired biopsies: within-person change in the severity instrument against within-person change in Kleiner stage, 58 pairs, prespecified. The only longitudinal substrate in the project. One week.
6. Optional upside: frozen-community holdout test for the fibroblast-cholangiocyte coordination result. If it reproduces, it becomes Figure 4's biology; if not, delete HIGH_IMPACT rows 2 and 3.

## Writing rules for this submission

- Abstract: two numbers minimum (22.4% vs 3.4%; 585 participants at 0.68 to 0.79 of ceiling), no placeholders, at most one hedge.
- Title states the result. Candidates: "Inherited risk and tissue remodeling are separable layers of human MASLD"; "Protein-altering variation and a diffuse severity axis partition human MASLD".
- Discussion: at most 0.25 hedges per sentence; one Limitations paragraph.
- Update title.md, STATUS.md, README.md, legends to seven figures in one edit.

## Submission path

1. Presubmission inquiry to Nature Metabolism (27.5) citing Kamzolas as the precedent and stating the complementary claim in three sentences.
2. In parallel, presubmission inquiry to Cell Metabolism (37.0) and Journal of Hepatology (40.1); both outrank Nature Metabolism on 2025 IF and both publish human multi-cohort MASLD work. J Hepatol needs the patient-level framing of leg 2 up front.
3. Fallbacks: Nature Communications (18.1), then Cell Genomics (9.0).

## Risks

- Leg 3's class-difference test may come back null; then "separable" becomes "complementary at chance", and the paper leans on legs 1 and 2 only.
- Nature Metabolism may read the diffuse-axis result as contradicting Kamzolas's systems; state it as a bound on within-stage resolution, with their axis reproduced.
- A public-data-only paper with no new cohort remains a desk-reject risk at every journal above 20; the exposure audit and the 585-participant external transfer are the counter.

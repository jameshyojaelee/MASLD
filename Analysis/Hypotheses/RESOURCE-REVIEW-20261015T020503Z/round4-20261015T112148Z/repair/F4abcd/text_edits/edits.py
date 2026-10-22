"""F4abcd text edits, returned as data only. Each old_text must occur exactly once
in its file; this script checks that and records the real line span."""
import json
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
ABS = "docs/manuscript/draft/abstract.md"
F4 = "docs/manuscript/draft/fig4.md"
F1 = "docs/manuscript/draft/fig1.md"
F5D = "docs/manuscript/draft/fig5_discussion.md"
INTRO = "docs/manuscript/draft/01_intro.md"
TITLE = "docs/manuscript/draft/title.md"
LEG = "docs/manuscript/05_figure_legends.md"
METH = "docs/manuscript/METHODS.md"
MSR = "docs/manuscript/README.md"
PAPER = "docs/PAPER.md"
RES = "docs/RESULTS.md"
STAT = "docs/STATUS.md"
READ = "docs/README.md"
ROAD = "docs/ROADMAP.md"
FIG4README = "figures/main/fig4_singlecell_programs/README.md"

E = []


def add(fid, file, old, new, claim=False, tag=""):
    E.append(dict(fid=fid, file=file, old_text=old, new_text=new, claim_level=claim, tag=tag))


# ---------------- R04-01: H59 null propagated (also carries R04-02, R04-03, G1-01 numbers)
add("R04-01", ABS,
    "Donor-aware single-cell analysis and a prespecified 117-program registry\n"
    "resolved remodeling into multicellular programs while ambient-RNA sensitivity\n"
    "limited lineage assignment. After source-diagnosis correction, Ductular injury\n"
    "alone met the primary 117-program rule; no program passed HC3 family correction.",
    "Donor-aware single-cell analysis scored a prespecified 117-program registry\n"
    "in five cell lineages, while ambient-RNA sensitivity limited lineage\n"
    "assignment. With dataset-specific variances modelled, no program passed the\n"
    "117-program family for a cross-sectional stage association.",
    claim=True)

add("R04-01", F4,
    "Figure 4C shows the three positive and three negative program effects with the\n"
    "smallest complete-family BH q per lineage, with all 117 programs retained in\n"
    "the supplementary view. ECM/IGFBP7 and ductular-injury/BICC1, both discovered in\n"
    "hepatocyte profiles, are among the 37 whose content contradicts that label.\n"
    "After correcting two GSE174748 explants to source-documented Cirrhosis and\n"
    "applying the unchanged cirrhosis exclusion, ordinal-stage models use at most\n"
    "62 biological donors across five datasets with dataset as a covariate.\n"
    "Only BICC1 passes the primary 117-program BH/stability rule (beta=0.4021,\n"
    "q=0.03763); ECM does not (beta=0.3515, q=0.08546). Neither passes HC3 family\n"
    "correction (BICC1 q=0.1558; ECM q=0.1652). The previous 64-donor two-program\n"
    "selection is historical.",
    "Figure 4C shows the three positive and three negative program effects with the\n"
    "smallest 117-program permutation BH q per lineage, with all 117 programs\n"
    "retained in the supplementary view. ECM/IGFBP7 and ductular-injury/BICC1, both\n"
    "discovered in hepatocyte profiles, are among the 37 whose content contradicts\n"
    "that label. After correcting two GSE174748 explants to source-documented\n"
    "Cirrhosis and applying the unchanged cirrhosis exclusion, the stage model uses\n"
    "the 60 donors whose datasets could be ambient-corrected. Only GSE202379 and\n"
    "GSE244832 contain more than one stage, so with dataset as a covariate they\n"
    "alone inform the slope: 55 hepatocyte, 38 cholangiocyte, 43 fibroblast, 25\n"
    "macrophage and 53 T-cell donors. Residual variance differs between these two\n"
    "datasets, so the model fits a separate variance per dataset (generalized least\n"
    "squares) and compares each program with 19,999 within-dataset stage\n"
    "permutations. No program passes the 117-program BH/stability rule (smallest\n"
    "q=0.323). Ductular injury has beta=0.137 (95% CI 0.012 to 0.261; permutation\n"
    "p=0.037, q=0.323); detection probability at the earlier estimate of 0.402 is\n"
    "0.999. Stromal ECM has beta=0.084 (−0.067 to 0.235; q=0.578); detection\n"
    "probability at its earlier 0.352 is 0.79, so this analysis cannot exclude an\n"
    "effect of that size. Stromal ECM is no longer among the rows shown in Figure\n"
    "4C. The two dataset slopes differ for both programs (Ductular injury 0.107 and\n"
    "0.742, Cochran Q p=0.016; Stromal ECM 0.037 and 0.756, Q p=0.0068), so each\n"
    "pooled beta is a precision-weighted average, not one stage effect. The\n"
    "62-donor equal-variance fit, which adds two GSE189600 donors without ambient\n"
    "correction, gave Ductular injury q=0.03763 (HC3 q=0.1558) and is a\n"
    "supplementary sensitivity; it has no call with the within-dataset permutation\n"
    "null (q=0.123), without those two donors (q=0.0985) or with a 50-cell lineage\n"
    "gate (q=0.0765). This null bounds only a\n"
    "pooled ordinal stage association in these two datasets; it does not show that\n"
    "the programs are absent. The previous 64-donor two-program selection is\n"
    "historical.")

add("R04-01", PAPER,
    "On at most 62 source-reported donors,\n"
    "   only Ductular injury passes the primary 117-program BH/stability rule\n"
    "   (β=0.402086, q=0.037633); Stromal ECM does not (β=0.351525, q=0.085458).\n"
    "   Neither passes HC3 family correction.",
    "With a separate residual variance per\n"
    "   dataset and a within-dataset permutation null on the 60 ambient-correctable\n"
    "   donors (stage-informative donors only in GSE202379 and GSE244832), no\n"
    "   program passes the 117-program BH/stability rule: Ductular injury β=0.1365,\n"
    "   q=0.323375; Stromal ECM β=0.0837, q=0.577711, which cannot exclude its\n"
    "   earlier effect size. The 62-donor equal-variance fit (Ductular injury\n"
    "   q=0.037633, HC3 q=0.155826) is a supplementary sensitivity.",
    claim=True)

add("R04-01", RES,
    "The source-diagnosis correction completed October 14 uses at most 62\n"
    "source-reported donors across five datasets, tests all 117 programs in one BH\n"
    "family, and adjusts for dataset. The two GSE174748 disease explants are\n"
    "cirrhotic in the original [Filliol study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9949942/),\n"
    "so the unchanged three-stage rule excludes them. Stromal ECM has disease\n"
    "β=0.351525, primary q=0.085458 and HC3 q=0.165240; it no longer passes primary\n"
    "selection. Ductular injury has β=0.402086, primary q=0.037633 and HC3 q=0.155826;\n"
    "it is the sole primary-selected program. Neither passes HC3 family correction.",
    "The source-diagnosis correction completed October 14 excludes the two\n"
    "GSE174748 disease explants, which are cirrhotic in the original [Filliol study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9949942/),\n"
    "under the unchanged three-stage rule. The primary stage model\n"
    "([H59](../Analysis/Hypotheses/H59-20261015T112432Z/RESULT_H59.md))\n"
    "uses the 60 donors whose datasets could be ambient-corrected. Only GSE202379\n"
    "and GSE244832 contain more than one stage, so they alone inform the slope\n"
    "(fitted donors: hepatocytes 55, cholangiocytes 38, fibroblasts 43,\n"
    "macrophages 25, T cells 53). It fits generalized least squares with a\n"
    "separate residual variance per dataset and tests each program against 19,999\n"
    "within-dataset stage permutations (seed 20260923), with BH over all 117\n"
    "programs and the unchanged stability rule. No program is selected; the\n"
    "smallest permutation q is 0.323375. Ductular injury has β=0.1365 (95% CI\n"
    "0.0121 to 0.2609), permutation p=0.0367, q=0.323375; detection probability at\n"
    "its earlier β=0.402086 is 0.999. Stromal ECM has β=0.0837 (−0.0672 to\n"
    "0.2346), q=0.577711; detection probability at its earlier β=0.351525 is\n"
    "0.794, so this analysis cannot exclude an effect of that size. Dataset slopes\n"
    "differ for both (Ductular injury 0.1074 and 0.7423, Cochran Q p=0.0164;\n"
    "Stromal ECM 0.0373 and 0.7556, Q p=0.00676). The 62-donor equal-variance fit\n"
    "(Stromal ECM β=0.351525, q=0.085458, HC3 q=0.165240; Ductular injury\n"
    "β=0.402086, q=0.037633, HC3 q=0.155826) is a sensitivity that depends on\n"
    "equal variances, on its parametric null and on two GSE189600 donors without\n"
    "ambient correction. It has zero calls with the same within-dataset\n"
    "permutation null (Ductular injury permutation q=0.12285), without the two\n"
    "GSE189600 donors (q=0.098523) or with a 50-cell lineage gate (q=0.076511).")

add("R04-01", RES,
    "| Programs meeting primary association/stability criteria | 1 | Diagnosis-qualified 62-donor maximum on the ambient-uncorrected atlas; BICC1 only; BH over all 117 programs |\n"
    "| Programs passing HC3 family correction | 0 | Corrected ECM HC3 q=0.165240 and BICC1 q=0.155826 |",
    "| Programs meeting primary association/stability criteria | 0 | H59: 60 ambient-correctable donors, per-dataset residual variance, within-dataset permutation BH over all 117 programs; smallest q=0.323375 |\n"
    "| Programs meeting the 62-donor equal-variance sensitivity | 1 | BICC1 only (q=0.037633); zero without the two GSE189600 donors (q=0.098523) |\n"
    "| Programs passing HC3 family correction (62-donor sensitivity) | 0 | Corrected ECM HC3 q=0.165240 and BICC1 q=0.155826 |")

add("R04-01", RES,
    "| Atlas label | Program | Primary q | HC3 q | Source-diagnosis-refit state |\n"
    "|---|---|---:|---:|---|\n"
    "| Hepatocyte | Stromal ECM | 0.085458 | 0.165240 | Not primary-selected |\n"
    "| Hepatocyte | Ductular injury | 0.037633 | 0.155826 | Primary-selected; fails HC3 family correction |\n"
    "\n"
    "These two rows are the October 14 diagnosis-qualified, 62-donor maximum\n"
    "single-cell refit.",
    "| Atlas label | Program | H59 permutation q | 62-donor OLS q | 62-donor HC3 q | State |\n"
    "|---|---|---:|---:|---:|---|\n"
    "| Hepatocyte | Stromal ECM | 0.577711 | 0.085458 | 0.165240 | Not selected; cannot exclude the earlier effect (detection 0.794) |\n"
    "| Hepatocyte | Ductular injury | 0.323375 | 0.037633 | 0.155826 | Not selected; the 62-donor call is a sensitivity only |\n"
    "\n"
    "The H59 column is the 60-donor primary fit. The other two columns are the\n"
    "October 14 diagnosis-qualified, 62-donor maximum single-cell refit, now a\n"
    "sensitivity.")

add("R04-01", RES,
    "lncRNA-independent. Both source-diagnosis-refit primary-selected programs passed this\n"
    "lncRNA-removal sensitivity:",
    "lncRNA-independent. Both focal programs passed this\n"
    "lncRNA-removal sensitivity:")

add("R04-01", RES,
    "117-program stage context has one primary-selected program, BICC1/Ductular\n"
    "injury, zero HC3-selected programs, 116 ordinary nonsignificant programs and\n"
    "zero strict negative calls.",
    "117-program stage context has no selected program in the H59 primary fit\n"
    "(smallest permutation q=0.323375); the 62-donor equal-variance sensitivity has\n"
    "one (BICC1/Ductular injury) and zero HC3-selected programs. There are zero\n"
    "strict negative calls.")

add("R04-01", STAT,
    "disease explants to Cirrhosis. The unchanged three-stage rule excludes them:\n"
    "maximum primary n=62, BICC1 alone primary-selected, zero HC3-selected;\n"
    "maximum ambient common n=60. Figure 4C and the full 117-program supplement\n"
    "use the corrected effects.",
    "disease explants to Cirrhosis. The unchanged three-stage rule excludes them.\n"
    "H59 (2026-10-15) refits all 117 programs on the 60-donor ambient-correctable\n"
    "universe with a per-dataset residual variance and a within-dataset\n"
    "permutation null: zero programs selected (smallest q=0.323375). The 62-donor\n"
    "equal-variance fit (BICC1 alone, q=0.037633; zero HC3-selected) becomes a\n"
    "sensitivity. Figure 4C and S4C show the H59 effects; their source panels are\n"
    "`figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/repair/` and\n"
    "replace the 62-donor panels together with this text.")

add("R04-01", STAT,
    "primary maximum n=62, BICC1 alone selected (q=0.037633), ECM not selected (q=0.085458), zero HC3-selected.",
    "H59 primary fit (60 donors, per-dataset variance, permutation BH over 117) selects zero programs (Ductular injury q=0.323375, ECM q=0.577711); the 62-donor equal-variance fit (BICC1 q=0.037633, ECM q=0.085458, zero HC3-selected) is a sensitivity.")

add("R04-01", READ,
    "Original-source diagnosis correction also refits the fixed 117-program\n"
    "stage models: maximum n=62, BICC1 alone primary-selected, zero HC3-selected.",
    "Original-source diagnosis correction also refits the fixed 117-program\n"
    "stage models. With per-dataset variances on the 60-donor ambient-correctable\n"
    "universe (H59), zero programs are selected; the 62-donor BICC1 call is a\n"
    "sensitivity only.")

add("R04-01", ROAD,
    "The fixed-score refit has maximum primary n=62,\n"
    "BICC1 alone selected and zero HC3-selected; its ambient common universe has\n"
    "maximum n=60.",
    "The fixed-score refit has maximum n=62\n"
    "(BICC1 alone selected, zero HC3-selected). H59 refits the 60-donor ambient\n"
    "common universe with per-dataset variances and selects zero programs, so the\n"
    "62-donor call becomes a sensitivity.")

add("R04-01", MSR,
    "Diagnosis-qualified 62/60-donor primary/ambient fits; Ductular injury alone primary-selected, zero HC3-selected; fixed weights and unresolved lineage origin",
    "H59 60-donor per-dataset-variance fit selects zero of 117 programs; the 62-donor Ductular injury call is a sensitivity; fixed weights and unresolved lineage origin")

add("R04-01", F5D,
    "The diagnosis-qualified refit retains primary support only for Ductular injury;\n"
    "Stromal ECM no longer passes the 117-program BH/stability rule, and neither\n"
    "passes HC3 family correction.",
    "With dataset-specific variances modelled, neither program passes the\n"
    "117-program BH/stability rule for a pooled stage association; the analysis\n"
    "cannot exclude the earlier Stromal ECM effect size.")

add("R04-01", INTRO,
    "Source-diagnosis correction retains only Ductular injury\n"
    "under the primary 117-program rule; neither focal program passes HC3 family\n"
    "correction.",
    "Neither focal program passes the 117-program rule for a\n"
    "cross-sectional stage association once dataset-specific variances are\n"
    "modelled.")

# ---------------- R04-02: model description (variance model)
add("R04-02", METH,
    "The corrected primary ordinal models therefore use at most 62 donors across\n"
    "five datasets, adjusting for dataset; the former 64-donor fit is historical.",
    "The corrected primary ordinal models use the 60 donors whose datasets could\n"
    "be ambient-corrected (GSE189600 excluded). Each program score is modelled as\n"
    "score ~ stage + dataset by generalized least squares with a separate residual\n"
    "variance per dataset (REML; nlme gls with varIdent), because the residual\n"
    "standard deviation differed between GSE244832 and GSE202379: the ratio is\n"
    "3.6 for Ductular injury and 2.8 for Stromal ECM (variance ratios about 13 and\n"
    "8), with a median of 1.6 and a range of 0.07 to 15.3 across the 117 programs.\n"
    "Donors in the two single-stage datasets (GSE174748 and GSE185477, Healthy\n"
    "only) carry no stage information and are not fitted. The p-value is\n"
    "(1 + the number of 19,999 within-dataset stage permutations, seed 20260923,\n"
    "with |t| at least the observed value) / 20,000; Benjamini-Hochberg correction\n"
    "is applied over all 117 programs, and selection also requires the unchanged\n"
    "stability median of at least 0.5. In simulation, the model-based 95% interval\n"
    "covered the true effect in 94.2% (true beta 0) and 95.5% (true beta 0.30) of\n"
    "draws. The 62-donor equal-variance OLS fit, with BH over 117 programs and an\n"
    "HC3 sensitivity, is reported as a supplementary sensitivity; the former\n"
    "64-donor fit is historical.")

# ---------------- R04-03: informative donors in the 62-donor design
add("R04-03", RES,
    "The corrected disease-state refit contains at most 62\n"
    "biological donors across GSE185477, GSE202379, GSE244832, GSE174748, and the\n"
    "two analyzed GSE189600 donors, with no Liver Atlas rows.",
    "The corrected disease-state refit contains at most 62\n"
    "biological donors across GSE185477, GSE202379, GSE244832, GSE174748, and the\n"
    "two analyzed GSE189600 donors, with no Liver Atlas rows. With dataset as a\n"
    "covariate, GSE174748 (2 Healthy) and GSE185477 (3 Healthy) contain one stage\n"
    "and do not inform the stage slope. In that 62-donor fit the slope is estimated\n"
    "from 57 hepatocyte donors in three datasets and from 38 cholangiocyte, 43 fibroblast,\n"
    "25 macrophage and 53 T-cell donors in two datasets (GSE202379, GSE244832),\n"
    "with four macrophage Steatosis donors in total. The H59 primary fit also\n"
    "drops GSE189600, leaving 55 hepatocyte donors.")

# ---------------- R04-04: ambient attenuation, nominal p beside q
add("R04-04", F4,
    "ECM raw beta=0.3509, corrected beta=0.1279, and\n"
    "paired delta=-0.2230 (delta BH q=0.4781); BICC1 raw beta=0.3846,\n"
    "corrected beta=0.4010, and delta=0.01642 (delta BH q=1). The former\n"
    "62-donor ambient values are historical. BICC1 remains\n"
    "`attenuation_indeterminate`; ECM is\n"
    "`not_applicable_not_frozen_disease_association` after losing primary selection.",
    "ECM raw beta=0.3509, corrected beta=0.1279, and\n"
    "paired delta=-0.2230, a 64% reduction (nominal p=0.0064, the smallest\n"
    "ordinary p of 106 computable paired changes; HC3 p=0.024, eighth by HC3 p;\n"
    "117-program delta BH q=0.4781);\n"
    "BICC1 raw beta=0.3846, corrected beta=0.4010, and delta=0.01642 (nominal\n"
    "p=0.67; delta BH q=1). Under the primary per-dataset-variance model the\n"
    "ambient-corrected betas are 0.052 for ECM and 0.143 for BICC1 (permutation\n"
    "q=0.559 and 0.270). The former 62-donor ambient values are historical.\n"
    "Neither program is selected in the primary fit, so the selection-conditioned\n"
    "attenuation classification applies to neither\n"
    "(`not_applicable_not_frozen_disease_association`).")

add("R04-04", PAPER,
    "Neither paired change passes the 117-program\n"
    "   BH family; ECM's selection-conditioned attenuation classification is no\n"
    "   longer applicable.",
    "Neither paired change passes the 117-program\n"
    "   BH family, although ECM's 64% reduction has nominal p=0.0064, the smallest\n"
    "   ordinary p of 106 computable paired changes (eighth by HC3 p). The\n"
    "   selection-conditioned attenuation classification applies to neither\n"
    "   program.")

add("R04-04", RES,
    "Its paired correction-induced change is β=-0.223006\n"
    "(BH q=0.478071; HC3 q=0.337646), and its",
    "Its paired correction-induced change is β=-0.223006,\n"
    "a 63.6% reduction (nominal p=0.006417, the smallest ordinary p of 106\n"
    "computable paired changes; HC3 p=0.024386, eighth by HC3 p, where\n"
    "cholangiocyte Inflammatory / oxidative (SOD2) is first at 0.00307; BH\n"
    "q=0.478071; HC3 q=0.337646), and its")

add("R04-04", RES,
    "BICC1 remains\n"
    "`attenuation_indeterminate`; ECM's selection-conditioned attenuation state is\n"
    "not applicable after losing primary selection.",
    "Neither program is\n"
    "selected in the H59 primary fit, so the selection-conditioned attenuation\n"
    "state is not applicable to either. Under the H59 model the ambient-corrected\n"
    "betas are 0.0518 for Stromal ECM (permutation q=0.55887) and 0.1428 for\n"
    "Ductular injury (q=0.26975).")

# ---------------- R04-05: Fig 4C legend (panel CI meaning, rows, heatmap columns); also carries G1-03 clause
add("R04-05", LEG,
    "**(C)** The three positive and three\n"
    "negative Hotspot effects with the smallest complete-family BH q in each of the\n"
    "four non-T lineages, shown across Healthy, Steatosis, and Steatohepatitis and\n"
    "paired with signed disease-stage β estimates and 95% confidence intervals. T\n"
    "cell programs are omitted from this compact main panel; the complete 117-program\n"
    "registry view, including T cell programs, is shown in Figure S4. The corrected\n"
    "62-donor primary fit excludes the two GSE174748 cirrhotic explants under the\n"
    "original stage rule. Only Ductular injury passes the 117-program BH/stability\n"
    "rule (beta 0.4021, q=0.03763); Stromal ECM does not (beta 0.3515, q=0.08546).\n"
    "Neither passes HC3 family correction (q=0.1558 and 0.1652, respectively).\n"
    "The prior 64-donor two-program selection is historical.",
    "**(C)** The three positive and three\n"
    "negative Hotspot effects with the smallest 117-program permutation BH q in\n"
    "each of the four non-T lineages. Heatmap columns are median dataset-centered\n"
    "donor scores in Healthy, Steatosis and Steatohepatitis, using only GSE202379\n"
    "and GSE244832, the two datasets with more than one stage. Points are ordinal\n"
    "stage β from generalized least squares with a separate residual variance per\n"
    "dataset, on the 60 donors whose datasets could be ambient-corrected; fitted\n"
    "donors are hepatocytes 55, cholangiocytes 38, fibroblasts 43 and macrophages\n"
    "25. Black bars are nominal 95% confidence intervals. Gray bars are β ± the\n"
    "program's within-dataset permutation critical |t| at 0.05/117 times its\n"
    "standard error. No program passes the 117-program permutation BH/stability\n"
    "rule (smallest q=0.323), and no gray bar excludes zero. T cell programs are\n"
    "omitted from this compact main panel; the complete 117-program registry view,\n"
    "including T cell programs, is shown in Figure S4. In the T-cell fits, 17 of\n"
    "53 donors have fewer than 20 T cells (minimum 1). The two GSE174748 cirrhotic\n"
    "explants are excluded under the original stage rule. The 62-donor\n"
    "equal-variance fit, in which Ductular injury had q=0.03763 (HC3 q=0.1558),\n"
    "is a supplementary sensitivity. The prior 64-donor two-program selection is\n"
    "historical.")

# ---------------- R04-06: S4 native forest legend
add("R04-06", LEG,
    "The native forest\n"
    "shows cohort and cohort-omission disease-grade effects in hepatocytes; its\n"
    "intervals are nominal sensitivity intervals, not new complete-family calls.",
    "The native forest\n"
    "shows hepatocyte disease-grade effects in the corrected 62-donor\n"
    "equal-variance sensitivity: the pooled fit (n=62), the two datasets with more\n"
    "than one stage (GSE202379, n=37; GSE244832, n=18) and five dataset omissions\n"
    "(n=60, 59, 60, 25 and 44 when omitting GSE174748, GSE185477, GSE189600,\n"
    "GSE202379 and GSE244832). GSE174748 and GSE185477 contain Healthy donors only\n"
    "and GSE189600 one donor per stage, so they have no dataset row. Intervals are\n"
    "nominal HC3 sensitivity intervals, not new complete-family calls.")

# ---------------- R04-07: V3 independent-donor test reported
add("R04-07", F4,
    "Taken with the content audit, this prevents any\n"
    "hepatocyte-intrinsic claim from the original atlas.",
    "Taken with the content audit, this prevents any\n"
    "hepatocyte-intrinsic claim from the original atlas. In independent donors\n"
    "(GSE256398 single-nucleus hepatocyte pseudobulk; score ~ stage + sex + age or\n"
    "disease + sex + age; HC3; Holm over four tests), no test passes Holm\n"
    "correction. Ductular injury has stage beta=0.694 (95% CI 0.062 to 1.326; 13\n"
    "non-cirrhotic donors; Holm p=0.139) and disease-versus-healthy beta=1.080\n"
    "(−0.091 to 2.250; Holm p=0.203) in 17 donors, whose 11 disease donors include\n"
    "4 with MASH cirrhosis; Stromal ECM has 0.300 (−0.709 to 1.309) and 0.176\n"
    "(−0.840 to 1.192), Holm p=1. With 13 and 17 donors these tests have low\n"
    "power, and no positive control was run.")

add("R04-07", RES,
    "These results support same-atlas descriptive localization, not transcript\n"
    "origin or cell autonomy.",
    "These results support same-atlas descriptive localization, not transcript\n"
    "origin or cell autonomy.\n"
    "\n"
    "The only completed test of the two fixed programs in new donors (V3,\n"
    "GSE256398 single-nucleus hepatocyte pseudobulk, at least 50 nuclei per donor;\n"
    "HC3; Holm over four tests) passes no test. Ductular injury: stage beta 0.694\n"
    "(95% CI 0.062 to 1.326; n=13 non-cirrhotic donors; p=0.0347, Holm 0.139)\n"
    "and disease versus healthy 1.080 (−0.091 to 2.250; n=17; p=0.0677, Holm\n"
    "0.203). The disease group of that second contrast has 11 donors: 3 MASLD F0,\n"
    "4 MASH F2-F3 and 4 MASH cirrhosis, so unlike the three-stage atlas models it\n"
    "includes cirrhosis. Stromal ECM: 0.300 (−0.709 to 1.309) and 0.176 (−0.840\n"
    "to 1.192), Holm 1. Power is low at these donor counts and no positive control\n"
    "was named ([V3 table](../Analysis/Hypotheses/V3-20261008T015239Z/results/v3_primary_family.tsv)).\n"
    "GSE289173 (48 donors, single-nucleus RNA-seq) was opened to test whether the\n"
    "Ductular injury score is higher in later-stage donors, but the test could not\n"
    "be run because only 14 of 21 program genes, carrying 0.64 of the frozen\n"
    "program weight, were detected in at least 1% of hepatocyte nuclei\n"
    "(prespecified minimum 0.80; [H90](../Analysis/Hypotheses/H90-20261015T120337Z/RESULT_H90.md)).")

add("R04-07", METH,
    "reported. This analysis describes cellular context within the same atlas and\n"
    "does not test transcript origin, secretion, or cell autonomy.",
    "reported. This analysis describes cellular context within the same atlas and\n"
    "does not test transcript origin, secretion, or cell autonomy.\n"
    "\n"
    "The two fixed programs were also scored in an independent single-nucleus\n"
    "dataset (GSE256398) as hepatocyte pseudobulk from donors with at least 50\n"
    "hepatocyte nuclei. Two contrasts were prespecified: ordinal stage in\n"
    "non-cirrhotic donors (score ~ stage + sex + age; n=13) and disease versus\n"
    "healthy (score ~ disease + sex + age; n=17, with 6 healthy and 11 disease\n"
    "donors, 4 of them with MASH cirrhosis). Inference used HC3 standard errors\n"
    "with Holm correction over the four program-by-contrast tests.")

# ---------------- R04-08: 4D text matches the installed panel
add("R04-08", F4,
    "Figure 4D asks the separate question of cellular context. In donor-paired\n"
    "same-atlas contrasts, Stromal ECM/IGFBP7 is highest in fibroblasts and Ductular\n"
    "injury/BICC1 is highest in cholangiocytes. The ten prespecified contrasts pass\n"
    "both ordinary and HC3 family criteria for those focal lineages. This is\n"
    "descriptive localization, not transcript origin or cell autonomy.",
    "Figure 4D asks the separate question of cellular context. It places each of\n"
    "the 117 programs by discovery lineage and content-associated lineage: 80 match\n"
    "and 37 do not. In separate donor-paired same-atlas contrasts (Figure S4),\n"
    "Stromal ECM/IGFBP7 is highest in fibroblasts (beta=2.271 over hepatocytes; 53\n"
    "donors) and Ductular injury/BICC1 in cholangiocytes (beta=1.574; 51 donors);\n"
    "both pass the ten-contrast BH and HC3 families. This is descriptive\n"
    "localization, not transcript origin or cell autonomy.")

# ---------------- R04-09: bootstrap uncertainty on 80/37
add("R04-09", F4,
    "Scored against donor-level per-lineage expression, **80 of 117\n"
    "programs match their declared lineage and 37 do not**, and the mismatch is\n"
    "systematic:",
    "Scored against donor-level per-lineage expression, **80 of 117\n"
    "programs match their declared lineage and 37 do not**. In 1,000 donor\n"
    "bootstrap draws, 76 matches and 34 mismatches keep the same call in at least\n"
    "95% of draws and 7 programs are unstable (10 when donors are resampled within\n"
    "each lineage); the per-draw match count has a central 95% range of 77 to 80.\n"
    "The mismatch is\n"
    "systematic:")

add("R04-09", LEG,
    "Eighty match their discovery lineage and\n"
    "37 do not.",
    "Eighty match their discovery lineage and\n"
    "37 do not. In 1,000 donor bootstrap draws (seed 20260923), 76 matches and 34\n"
    "mismatches keep the same call in at least 95% of draws and 7 are unstable;\n"
    "Stromal ECM (fibroblast content) and Ductular injury (cholangiocyte content)\n"
    "keep their call in every draw.")

# ---------------- R04-10: ECM 'disease-associated program'
add("R04-10", INTRO,
    "ECM/IGFBP7 and ductular-injury/BICC1 test whether a disease-associated program,\n"
    "its cell-of-origin assignment, and its intact-tissue organization support the\n"
    "same conclusion.",
    "ECM/IGFBP7 and ductular-injury/BICC1 test whether a program's disease\n"
    "association, its cell-of-origin assignment, and its intact-tissue organization\n"
    "support the same conclusion. The bulk fibrosis association of ECM/IGFBP7\n"
    "reaches BH q=0.045 (224 program-by-axis tests) only in the 457-record\n"
    "four-cohort arm without composition adjustment. It is q=0.214 in the\n"
    "364-unit three-cohort arm, and q=0.23 to 0.25 in the four-cohort arm after\n"
    "composition adjustment, depending on the cell-fraction set.")

# ---------------- R10-01: vocabulary null on the 457-record Fig 4E arm
add("R10-01", F4,
    "Likewise, the existing Hotspot matched-background comparison\n"
    "does not establish exceptional disease specificity (fibrosis 26/114 versus\n"
    "null mean 18.445, nominal P=0.134; NAS 4 versus 4.225, P=0.512).",
    "Likewise, the Hotspot matched-background comparison does not\n"
    "establish disease specificity. With unweighted scores on the 457-record\n"
    "four-cohort arm, 20 of 114 programs reach fibrosis BH q<0.05 (228-test\n"
    "family) against a mean of 26.9 for size- and expression-matched random gene\n"
    "sets (empirical P=0.871, 200 draws), and 2 reach NAS against 2.57 (P=0.687).\n"
    "On the 364-unit three-cohort arm, 1 program reaches fibrosis q<0.05 against\n"
    "1.24 (P=0.473).")

add("R10-01", RES,
    "the vocabulary null (26 of 114\n"
    "fibrosis-supported against 18.4 expected, p = 0.134) now has a cause,",
    "the vocabulary null (469-record roster: 26 of 114\n"
    "fibrosis-supported against 18.4 expected, p = 0.134; 457-record four-cohort\n"
    "arm: 20 against 26.9, p = 0.871) now has a cause,")

add("R10-01", RES,
    "The 27/113 belongs to the 469-record working roster;\n"
    "the 774-unit arms are in the Figure 4E paragraph below.",
    "The 27/113 belongs to the 469-record working roster.\n"
    "The 26/114 vocabulary count also comes from that roster. On the 457-record\n"
    "four-cohort arm the unweighted vocabulary count is 20/114 against a\n"
    "matched-random mean of 26.91 (p = 0.871), and NAS 2 against 2.565\n"
    "(p = 0.687). On the 364-unit three-cohort arm it is 1 against 1.24\n"
    "(p = 0.473). The 774-unit weighted arms are in the Figure 4E paragraph below.")

# ---------------- R13-04: 'resolved remodeling into multicellular programs'
add("R13-04", INTRO,
    "Donor-aware\n"
    "single-cell analyses and prespecified gene programs resolve multicellular\n"
    "remodeling, while protein,",
    "Donor-aware\n"
    "single-cell analyses score prespecified gene programs in each cell lineage,\n"
    "while protein,")

add("R13-04", F1,
    "Figure 4\n"
    "resolves donor-level multicellular programs and regulatory activity.",
    "Figure 4\n"
    "reports donor-level program scores in five cell lineages and their bulk-tissue\n"
    "context; transcription-factor activity is supplementary.")

add("R13-04", PAPER,
    "2. Multicohort bulk and donor-aware single-cell analyses resolve established\n"
    "   disease as multicellular, cross-sectional remodeling.",
    "2. Multicohort bulk analyses describe established disease as cross-sectional\n"
    "   remodeling, and donor-aware single-cell analyses place fixed programs in\n"
    "   cell lineages; no single-cell program has a family-level stage association.",
    claim=True)

# Title rows removed in repair: clause 2 is removed under TITLE_CLAUSE_DECISION.md (row H52 untestable /
# H54 supports); the replacement title is drafted by Claims R13-01 across title.md, PAPER, README and STATUS.

# ---------------- Conflict 12: Fig 4 README 4C row (stale single-call statement)
add("R04-01", FIG4README,
    "| 4C | Twenty-four-program Hotspot stage heatmap across four non-T lineages; all 117 programs, including T cell programs, in Figure S4 | source-diagnosis-corrected model, maximum 62 donors; Ductular injury is the only selected primary call |",
    "| 4C | Twenty-four-program Hotspot stage heatmap across four non-T lineages; all 117 programs, including T cell programs, in Figure S4 | H59 primary fit (60 donors, per-dataset variance, permutation BH over 117): no program selected; the 62-donor Ductular injury call is a sensitivity only |")

# ---------------- G1-01: gate scope
add("G1-01", METH,
    "The 50-cell gate is primary; 20 cells and\n"
    "annotation confidence ≥0.90 are sensitivities.",
    "The 50-cell gate is primary; 20 cells and\n"
    "annotation confidence ≥0.90 are sensitivities. The ten donor-paired\n"
    "lineage-localization contrasts use the same 50-cell gate. The 117-program\n"
    "stage models do not: they score the cells retained for program discovery,\n"
    "which keeps a library's hepatocytes, cholangiocytes, fibroblasts or\n"
    "macrophages only when that library has at least 30 cells of the lineage and\n"
    "keeps T cells without a minimum. In the 62-donor equal-variance sensitivity,\n"
    "a 50-cell lineage gate removes the single Ductular injury call (q 0.0376 to\n"
    "0.0765).")

# ---------------- G1-02: GSE189600 provenance
add("G1-02", METH,
    "after two prespecified low-cell\n"
    "GSE189600 libraries are excluded.",
    "after two prespecified low-cell\n"
    "GSE189600 libraries are excluded. GSE189600 deposits two sequencing runs for\n"
    "each of three healthy and three NASH replicate samples; the atlas ingests one\n"
    "run from four of them, and two libraries (one Healthy, one Steatohepatitis;\n"
    "142 nuclei in the scored lineages) are analyzed. The deposit does not state\n"
    "whether the replicate labels are distinct people. Because ambient correction\n"
    "failed for this dataset, it is outside the primary stage-model universe and\n"
    "enters only the 62-donor sensitivity.")

# ---------------- G1-03: cells per donor and lineage
add("G1-03", RES,
    "| Cell types represented | 5 | Hepatocyte, macrophage, fibroblast, cholangiocyte, T cell |",
    "| Cell types represented | 5 | Hepatocyte, macrophage, fibroblast, cholangiocyte, T cell |\n"
    "| Fitted donors with fewer than 50 cells of the tested lineage (H59 fit) | T cells 28 of 53 (17 below 20; minimum 1); cholangiocytes 6 of 38; fibroblasts 10 of 43; macrophages 5 of 25; hepatocytes 0 of 55 | Donor means enter the stage models unweighted; no T-cell program is near selection (smallest T-cell permutation p 0.1947, q 0.520) |")

# ---------------- G1-04: per-cell QC disclosure
add("G1-04", METH,
    "Dataset identity was retained during\n"
    "integration.",
    "Dataset identity was retained during\n"
    "integration. Program scoring used the integrated atlas with its original\n"
    "per-cell filters: a 200-count floor, no mitochondrial-fraction filter, no\n"
    "doublet removal and no ambient-RNA correction. A later harmonized QC\n"
    "(mitochondrial fraction <20%, at least 200 genes and 500 counts, Scrublet\n"
    "score <0.25) retains 1,030,882 cells. The share of scored-lineage cells that\n"
    "pass it is 49.4% in GSE185477, 79.3% in GSE244832, 89.1% in GSE174748,\n"
    "97.0% in GSE202379 and 7.7% in GSE189600. Rescoring on passing cells keeps\n"
    "Ductular injury donor scores at Spearman 0.995 with the published scores.")


def main():
    out = []
    bad = 0
    for e in E:
        txt = (ROOT / e["file"]).read_text()
        n = txt.count(e["old_text"])
        if n != 1 and not e["tag"]:
            print("COUNT", n, e["fid"], e["file"], repr(e["old_text"][:80]))
            bad += 1
            continue
        if n != 1 and e["tag"]:
            print("COUNT", n, e["fid"], e["file"], "conditional")
            bad += 1
            continue
        i = txt.index(e["old_text"])
        l0 = txt[:i].count("\n") + 1
        l1 = l0 + e["old_text"].count("\n")
        line = f"{l0}-{l1}" if l1 > l0 else f"{l0}"
        if e["tag"]:
            line += " [" + e["tag"] + "]"
        if "—" in e["new_text"]:
            print("EMDASH", e["fid"])
            bad += 1
        out.append(dict(fid=e["fid"], file=e["file"], line=line, old_text=e["old_text"],
                        new_text=e["new_text"], claim_level=e["claim_level"]))
    # Sequential application per file: every old_text must still occur exactly once after the
    # edits listed before it are applied (shared-line edits are plain string replacements).
    files = {}
    for e in E:
        files.setdefault(e["file"], []).append(e)
    for f, es in files.items():
        txt = (ROOT / f).read_text()
        for e in es:
            n = txt.count(e["old_text"])
            if n != 1:
                print("SEQ", n, e["fid"], f, repr(e["old_text"][:60]))
                bad += 1
                continue
            txt = txt.replace(e["old_text"], e["new_text"], 1)
    for e in E:
        for bad_str in ("[Inference", "[Unverified", "[Estimate", "[Speculation", "rests on", "\u2014"):
            if bad_str in e["new_text"]:
                print("BADSTR", bad_str, e["fid"], e["file"])
                bad += 1
    Path(__file__).with_name("edits.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print("edits", len(out), "bad", bad)
    for o in out:
        print(o["fid"], o["file"], o["line"][:60], o["claim_level"])


if __name__ == "__main__":
    main()

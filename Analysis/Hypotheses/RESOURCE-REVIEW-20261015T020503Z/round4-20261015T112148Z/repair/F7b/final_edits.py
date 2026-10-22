"""F7b repair: final text edits. Checks every old_text occurs exactly once in the current
file, that no new_text has an em dash or an evidence tag, and writes final_edits.json."""
import json
import re

DOCS = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/docs/"
OUT = ("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/"
       "RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F7b/final_edits.json")
LEG = "docs/manuscript/05_figure_legends.md"
F7 = "docs/manuscript/draft/fig7.md"
ME = "docs/manuscript/METHODS.md"
RE_ = "docs/RESULTS.md"
ST = "docs/STATUS.md"
PA = "docs/PAPER.md"
DI = "docs/manuscript/draft/fig5_discussion.md"
TAB = "figures/candidates/review-r3-fixes-20261015T112148Z/F7b/repair/fig7_transfer_primary_units.tsv"

E = {}

E["R08-02"] = [
    dict(file=LEG, line="508-511", claim_level=True,
         old_text=("The severity instrument these panels bound transfers out of cohort to 585\n"
                   "participants in three cohorts it never trained on, reaching 0.68 to 0.79 of the\n"
                   "tie ceiling each cohort's own fibrosis coding allows, and in each cohort it\n"
                   "matches rather than beats a ridge fitted inside that cohort. Historical"),
         new_text=("The severity instrument these panels bound transfers out of cohort to three\n"
                   "cohorts whose records it was not fitted on (panel D): 109 GSE268273 participant\n"
                   "records, 263 GSE213621 genetic identities (one record each, after removing 11\n"
                   "identities whose two records carry different stage labels) and 81 GSE276114\n"
                   "participants with MASLD. It reaches 0.70 to 0.87 of the tie ceiling each\n"
                   "cohort's own fibrosis coding allows. GSE268273 was also the cohort on which the\n"
                   "shipped expression transform was chosen among four eligible arms (Spearman\n"
                   "0.6506 to 0.6631). In the 82 GSE276114 participants with chronic viral\n"
                   "hepatitis, whose fibrosis the source study staged with METAVIR rather than the\n"
                   "Kleiner system, it reaches 0.68 of the ceiling. The source study's participants\n"
                   "with hepatocellular carcinoma, sampled from peritumoral liver, are folded into\n"
                   "these etiology groups without a per-sample label: 9 or 10 of the 81 with MASLD\n"
                   "and 24 of the 82 with viral hepatitis. Paired differences from a ridge fitted\n"
                   "inside each cohort, computed on the published record sets, are +0.054, -0.040\n"
                   "and -0.011 Spearman, with 95% intervals that include zero and minimum detectable\n"
                   "effects of 0.08 to 0.12. No difference was detected, and equivalence is not\n"
                   "established. Historical")),
    dict(file=RE_, line="3175", claim_level=False,
         old_text=("| Same instrument, two additional external cohorts and three contrasts, Spearman / tie ceiling "
                   "| **0.677 / 0.774 / 0.794** | GSE213621 299, GSE213621 4-level 367, GSE276114 177 participant "
                   "records; `executions/external-validation-expansion-20260901T190000Z/results.json`. No demonstrated "
                   "advantage over the recorded transductive in-cohort ridge; intervals including zero do not "
                   "establish equivalence |"),
         new_text=("| Same instrument, two additional external cohorts, Spearman / tie ceiling "
                   "| **0.709** (GSE213621, one record per genetic identity; 0.789 control-inclusive) / **0.866** "
                   "(GSE276114 MASLD); 0.680 in GSE276114 chronic viral hepatitis "
                   "| GSE213621: 263 genetic identities, Spearman 0.6680 (95% CI 0.6020 to 0.7267), keeping each "
                   "identity's primary-selected record and removing 11 identities with discordant stage labels "
                   "(control-inclusive 331). GSE276114: 81 MASLD participants, Spearman 0.7860 (0.6985 to 0.8483); "
                   "82 chronic viral hepatitis, 0.5443 (0.3752 to 0.6878), staged by METAVIR rather than Kleiner in "
                   "the source study. Intervals resample units within stage (10,000 draws, seed 20260923). The "
                   "source folds 34 hepatocellular carcinoma patients (peritumoral liver; 10 MASLD and 24 viral "
                   "hepatitis among its 178 samples) into the etiology labels without a per-sample flag. "
                   "Published record-set sensitivity: 0.677 / 0.774 / 0.794 for GSE213621 299, GSE213621 4-level "
                   "367 and GSE276114 177 records "
                   "(`executions/external-validation-expansion-20260901T190000Z/results.json`). "
                   f"`{TAB}`. No demonstrated advantage over the recorded transductive in-cohort ridge "
                   "(computed on the published record sets); intervals including zero do not establish "
                   "equivalence |")),
    dict(file=ST, line="831", claim_level=False,
         old_text=("two additional external cohorts, three contrasts, Spearman/ceiling 0.677–0.794 "
                   "(no demonstrated advantage"),
         new_text=("two additional external cohorts, Spearman/ceiling 0.709 for GSE213621 at one record per "
                   "genetic identity (263) and 0.866 for the 81 GSE276114 participants with MASLD (0.680 in its "
                   "82 with chronic viral hepatitis, staged by METAVIR; published record sets 0.677–0.794) "
                   "(paired comparisons on published record sets: no demonstrated advantage")),
    dict(file=ME, line="835-840", claim_level=False,
         old_text=("Its GSE268273 evaluation\n"
                   "uses 109 participant records. The subsequent GSE213621 disease-only and\n"
                   "control-inclusive contrasts use 299 and 367 records from the same cohort;\n"
                   "GSE276114 contributes 177. These are two additional cohorts and three\n"
                   "contrasts, not three independent cohorts. Source-declared participant keys and\n"
                   "count-to-phenotype joins were checked; genetic identity was not verified."),
         new_text=("Its GSE268273 evaluation\n"
                   "uses 109 participant records. The subsequent GSE213621 disease-only and\n"
                   "control-inclusive contrasts use 299 and 367 records from the same cohort;\n"
                   "GSE276114 contributes 177. These are two additional cohorts and three\n"
                   "contrasts, not three independent cohorts. Source-declared participant keys and\n"
                   "count-to-phenotype joins were checked. GSE213621 libraries were genotype-screened\n"
                   "with the Resource identity screen: the 299 staged records are 274 genetic\n"
                   "identities and the 367 records are 342, with no identity shared with another\n"
                   "screened cohort. The primary GSE213621 estimate keeps each identity's\n"
                   "primary-selected record and removes the 11 identities whose two records carry\n"
                   "different stage labels (263 and 331 records). GSE268273 and GSE276114 were not\n"
                   "genotype-screened. GSE276114 records carry an etiology label (81 MASLD, 82\n"
                   "chronic viral hepatitis, 14 alcohol-related liver disease). The primary\n"
                   "GSE276114 estimate uses the 81 MASLD records. The viral hepatitis records are\n"
                   "reported separately: they are not MASLD, and the source study staged them with\n"
                   "METAVIR rather than the Kleiner system. The 14 alcohol-related records are all\n"
                   "F4 and carry no within-group stage ordering. The source study folds 34\n"
                   "participants with hepatocellular carcinoma, sampled from peritumoral liver,\n"
                   "into these etiology labels without a per-sample label (10 MASLD and 24 viral\n"
                   "hepatitis among its 178 samples; 177 are analysed here). Intervals for these\n"
                   "primary estimates resample participants or genetic identities within stage\n"
                   "(10,000 draws, seed 20260923), the scheme of the published record-level\n"
                   "intervals.")),
]

E["R08-03"] = [
    dict(file=F7, line="37-40", claim_level=False,
         old_text=("Its three primary evaluations comprise 109 GSE268273,\n"
                   "299 GSE213621 and 177 GSE276114 participant records, with source-declared\n"
                   "identity rather than independent genetic confirmation. Raw Spearman\n"
                   "correlations are 0.6631, 0.6382 and 0.6739."),
         new_text=("Its three published evaluations use 109 GSE268273 records, 299 GSE213621\n"
                   "records and 177 GSE276114 records. GSE268273 and GSE276114 identities are\n"
                   "source-declared and were not genotype-checked. In GSE213621 the Resource's\n"
                   "genotype screen places the 299 staged records in 274 genetic identities, 25 of\n"
                   "them with two libraries and 11 of those with two different stage labels.\n"
                   "GSE276114 contains 81 participants with MASLD, 82 with chronic viral hepatitis\n"
                   "and 14 with alcohol-related liver disease, all 14 at F4. The source study\n"
                   "staged viral hepatitis with METAVIR and MASLD with the Kleiner system, and it\n"
                   "folds 34 participants with hepatocellular carcinoma, sampled from peritumoral\n"
                   "liver, into these groups without a per-sample label. The primary estimates\n"
                   "therefore use 109 GSE268273 records, one record for each of 263 GSE213621\n"
                   "identities after removing the 11 discordant ones, and the 81 GSE276114\n"
                   "participants with MASLD. Their Spearman correlations are 0.6631 (95% interval\n"
                   "0.568 to 0.748), 0.668 (0.602 to 0.727) and 0.786 (0.698 to 0.848), which are\n"
                   "0.696, 0.709 and 0.866 of each cohort's tie ceiling. In the 82 participants with\n"
                   "chronic viral hepatitis it is 0.544 (0.680 of the ceiling). On the published\n"
                   "record sets the raw Spearman correlations are 0.6631, 0.6382 and 0.6739.\n"
                   "GSE268273 was also used to choose the shipped expression transform among four\n"
                   "eligible arms, so it is not an unused test cohort for that choice.")),
    dict(file=F7, line="42-43", claim_level=False,
         old_text=("cohort. In the two additional cohorts, recorded paired differences from native\n"
                   "ridge are -0.040"),
         new_text=("cohort. In the two additional cohorts, recorded paired differences from native\n"
                   "ridge on the published record sets, which count repeat GSE213621 libraries and\n"
                   "non-MASLD GSE276114 participants, are -0.040")),
]

E["R08-07"] = [
    dict(file=F7, line="46-47", claim_level=False,
         old_text=("The reference-cohort point difference is +0.054, without a newly reproduced\n"
                   "paired interval in this assessment."),
         new_text=("The reference-cohort paired difference is +0.054 (recorded 95% interval -0.025\n"
                   "to +0.137; minimum detectable effect 0.115).")),
]

E["R08-08"] = [
    dict(file=RE_, line="3383-3387", claim_level=False,
         old_text=("TranscriptFormer minus the liver encoder macro-F1 difference is -0.004170\n"
                   "(95% conditional interval -0.008732 to -0.000436) for the linear head and\n"
                   "-0.003490 (-0.008311 to -0.000181) for the two-layer head. These are descriptive\n"
                   "fixed-model comparisons, conditional on the five studies and evaluable cell\n"
                   "classes, not population-of-studies generalization or new significance claims."),
         new_text=("TranscriptFormer minus the liver encoder macro-F1 difference is -0.004170\n"
                   "(95% conditional interval -0.008732 to -0.000436) for the linear head and\n"
                   "-0.003490 (-0.008311 to -0.000181) for the two-layer head as equal-study means.\n"
                   "One of the five studies, GSE189600, has 2 donors; its own linear-head difference\n"
                   "is -0.0653. From the same saved draws, the donor-weighted linear-head\n"
                   "difference is +0.013781 (+0.008469 to +0.019264) and the four-study mean without\n"
                   "GSE189600 is +0.011123 (+0.005672 to +0.015973), so the linear-head sign is not\n"
                   "resolved. The two-layer head is -0.001274 (-0.004236 to +0.001679) donor-weighted\n"
                   "and -0.002254 (-0.007141 to +0.002323) without GSE189600. These are descriptive\n"
                   "fixed-model comparisons, conditional on the five studies and evaluable cell\n"
                   "classes, not population-of-studies generalization or new significance claims.\n"
                   "[Source table](../figures/candidates/review-r3-fixes-20261015T112148Z/F7b/"
                   "fig7_transcriptformer_clean_subset_weightings.tsv).")),
    dict(file=ME, line="828", claim_level=False,
         old_text="population of studies. No models or scores were refitted for this correction.",
         new_text=("population of studies. The same saved draws give donor-weighted means and means\n"
                   "without the 2-donor study GSE189600. The sign of the linear-head difference\n"
                   "depends on that choice and is not reported as resolved. No models or scores\n"
                   "were refitted for this correction.")),
    dict(file=F7, line="59-61", claim_level=False,
         old_text=("Two of twelve model-by-head blocks change the sign of their delta when studies\n"
                   "are weighted by donors instead, so we report both averages and make no sign\n"
                   "claim that rests on the unweighted one."),
         new_text=("Two of twelve model-by-head blocks change the sign of their delta when studies\n"
                   "are weighted by donors instead, so we report both averages and make no sign\n"
                   "claim that rests on the unweighted one. The same holds for the five-study clean\n"
                   "subset for TranscriptFormer. Its linear-head difference from the liver encoder\n"
                   "is -0.0042 (95% conditional interval -0.0087 to -0.0004) with equal study\n"
                   "weights, +0.0138 (+0.0085 to +0.0193) with donor weights and +0.0111 (+0.0057 to\n"
                   "+0.0160) without the 2-donor study, so its sign is not resolved. The two-layer\n"
                   "head gives -0.0035, -0.0013 and -0.0023, and only the equal-weight interval\n"
                   "excludes zero.")),
]

E["R08-09"] = [
    dict(file=LEG, line="495-496", claim_level=False,
         old_text=("each held-out study. Study-level deltas are reported both unweighted and\n"
                   "weighted by donors."),
         new_text=("each held-out study. Per-study deltas against the liver encoder and their\n"
                   "unweighted and donor-weighted means are not drawn in 7F; their ranges are\n"
                   "given below.")),
    dict(file=LEG, line="501-502", claim_level=False,
         old_text=("an unweighted mean over seven studies, and it contributes -0.45 to -0.49 where\n"
                   "every other study contributes -0.065 to -0.002."),
         new_text=("an unweighted mean over seven studies. Its delta is -0.45 to -0.49 in three of\n"
                   "twelve model-by-head blocks and -0.065 to -0.002 in the other nine, while the\n"
                   "other six studies range from -0.043 to +0.097, with 43 of their 72 deltas positive.")),
]

E["R08-10"] = [
    dict(file=ME, line="845-846", claim_level=False,
         old_text=("Cohorts used to inspect refusal thresholds or other\n"
                   "acceptance criteria are development-exposed for those criteria."),
         new_text=("Cohorts used to inspect refusal thresholds or other\n"
                   "acceptance criteria are development-exposed for those criteria. GSE268273 is\n"
                   "development-exposed for the choice of expression transform. Among the four\n"
                   "eligible arms, which have a coherent expression transform and can score one\n"
                   "sample at a time (out-of-cohort Spearman 0.6506 to 0.6631 on these 109\n"
                   "records), the release ships the one with the highest Spearman. That rule comes\n"
                   "from the original prespecification. Amendment 03 made coherence and\n"
                   "reliability preconditions after the first three transform families had been\n"
                   "read and before the quantile-normalization family, which the release ships,\n"
                   "had been run.")),
    dict(file=RE_, line="3174", claim_level=False,
         old_text=("| 109 participants, GSE268273; `release/masld-severity-v1.2/expected_ood.json` "
                   "(bitwise = v1.1) |"),
         new_text=("| 109 participants, GSE268273; `release/masld-severity-v1.2/expected_ood.json` "
                   "(bitwise = v1.1). GSE268273 also chose the shipped transform: the highest Spearman of four "
                   "coherent, deployable arms (0.6506 to 0.6631), under the prespec 2c rule that "
                   "`release-instrument-20260901T211017Z/prespec/AMENDMENT_03.md:53-54` kept; the amendment "
                   "was written before the quantile arms were run |")),
    dict(file=RE_, line="3220", claim_level=False,
         old_text="stage in the three cohorts never trained on (109 / 299 / 177 participants; Spearman",
         new_text=("stage in the three cohorts whose records it was not fitted on (109 / 299 / 177 records; "
                   "GSE268273 also chose the shipped transform; Spearman")),
    dict(file=F7, line="34", claim_level=False,
         old_text="For tissue-level severity, the released linear instrument was fitted to 521",
         new_text=("For tissue-level severity, the released linear instrument, whose expression "
                   "transform was chosen on GSE268273, was fitted to 521")),
    dict(file=DI, line="268", claim_level=False,
         old_text=("baseline: the severity instrument at 0.6631 against 0.6094 for the in-cohort ridge, "
                   "and the"),
         new_text=("baseline: the severity instrument at 0.6631 against 0.6094 for the in-cohort ridge\n"
                   "(GSE268273 was also used to choose the expression transform; the paired difference,\n"
                   "+0.054, has a 95% interval of -0.025 to +0.137, which includes zero), and the")),
]

E["R11-01"] = [
    dict(file=PA, line="181-190", claim_level=True,
         old_text=("**Evidence state: separability `tested_negative` in one cohort; attribution\n"
                   "`indeterminate`.** The separability of the three NASH-CRN components has now\n"
                   "been tested in GSE267145, in 99 participants carrying steatosis, ballooning and\n"
                   "lobular inflammation as separately recorded columns, against a prespecification\n"
                   "recorded before the analysis. The three components span a participation ratio of\n"
                   "1.42 of a possible 3.0 against a 2.0 floor, so they behave as one axis rather\n"
                   "than three, and an adequately observed, applicable, properly joined test returns\n"
                   "a negative. Attribution remains `indeterminate`: nothing can be assigned to\n"
                   "steatosis, ballooning, or lobular inflammation individually, now for a measured\n"
                   "reason rather than for want of a test."),
         new_text=("**Evidence state: separability `indeterminate` in one cohort; attribution\n"
                   "`indeterminate`.** The separability of the three NASH-CRN components was\n"
                   "examined in GSE267145, in 99 participants carrying steatosis, ballooning and\n"
                   "lobular inflammation as separately recorded columns, against an RNA\n"
                   "prespecification recorded before the analysis. Two of its three criteria are\n"
                   "properties of the histology grades alone. The grades correlate at Spearman\n"
                   "0.68, 0.73 and 0.82 and span a participation ratio of 1.42 of a possible 3.0,\n"
                   "below the 2.0 floor, so the rule could not return a positive result whatever\n"
                   "the RNA showed. The RNA criterion also failed: the per-gene association vectors\n"
                   "of the three grades correlate at 0.83, 0.87 and 0.87 against a 0.80 threshold.\n"
                   "In simulations on these 99 grades with independent genes, one shared program\n"
                   "and three programs that each follow one recorded grade give the same values\n"
                   "(0.82 to 0.83 when a fifth of genes carry signal). Programs that follow only the\n"
                   "part of each grade not shared with the other two give values below 0.66, and an\n"
                   "equal mixture of the two forms gives values below 0.75; the observed values\n"
                   "match neither. The separate H3K27ac analysis agrees: its per-region association\n"
                   "vectors give a participation ratio of 1.15, against a null centred at 1.46,\n"
                   "inside the simulated range for the first two forms (1.09 to 1.43) and below\n"
                   "that for grade-specific programs (1.48 to 2.08). The test therefore cannot tell\n"
                   "one shared program from separate programs that follow the correlated grades.\n"
                   "Attribution remains `indeterminate`: nothing can be assigned to steatosis,\n"
                   "ballooning, or lobular inflammation individually.")),
    dict(file=PA, line="192", claim_level=False,
         old_text="This settles one cohort and one assay. It does not settle the SAF partition,",
         new_text="This concerns one cohort. It does not address the SAF partition,"),
    dict(file=PA, line="195-196", claim_level=False,
         old_text="A negative on the NASH-CRN partition in one\nbulk cohort is not a negative on SAF.",
         new_text="An unresolved result on the NASH-CRN\npartition in one bulk cohort says nothing about SAF."),
    dict(file=ST, line="1115-1122", claim_level=False,
         old_text=("Separability of the\n"
                   "  three NASH-CRN components is now **`tested_negative` in GSE267145**: tested in\n"
                   "  99 participants against a prespecification recorded before the analysis, the\n"
                   "  components span a participation ratio of 1.42 of a possible 3.0 against a 2.0\n"
                   "  floor and behave as one axis rather than three. Aspect **attribution stays\n"
                   "  `indeterminate`**, now for a measured reason rather than for want of a test;\n"
                   "  nothing can be assigned to a single aspect on this evidence. The result\n"
                   "  settles one cohort and one assay and does not settle the SAF partition of"),
         new_text=("Separability of the\n"
                   "  three NASH-CRN components is **`indeterminate` in GSE267145**. In 99\n"
                   "  participants the histology grades alone fail two of the three prespecified\n"
                   "  criteria (Spearman 0.68/0.73/0.82; label participation ratio 1.42 against a 2.0\n"
                   "  floor), so the rule could not return a positive; the RNA criterion also fails\n"
                   "  (association-vector correlations 0.83/0.87/0.87 against 0.80). In simulation\n"
                   "  on these grades with independent genes, one shared program and three programs\n"
                   "  that each follow one recorded grade give the same values; programs tied to the\n"
                   "  part of each grade not shared with the other two give values the data do not\n"
                   "  show (RNA criterion below 0.66; H3K27ac participation ratio 1.48 to 2.08\n"
                   "  against the observed 1.15). Aspect **attribution stays `indeterminate`**;\n"
                   "  nothing can be assigned to a single aspect on this evidence. The result\n"
                   "  concerns one cohort and does not address the SAF partition of")),
    dict(file=RE_, line="3252-3253", claim_level=False,
         old_text=("does not revise the recorded `tested_negative` on NAS-component separability, which concerns\n"
                   "per-region association patterns rather than a donor-level factor conditional on the grades."),
         new_text=("does not change the NAS-component separability state, now `indeterminate`: its per-region\n"
                   "association statistic takes the same values under one shared program and under three programs\n"
                   "that each follow their recorded grade. That statistic concerns per-region association patterns\n"
                   "rather than a donor-level factor conditional on the grades.")),
]

TAGS = re.compile(r"\[(Unverified|Inference|Speculation|Estimate|Contains Unverified Content)\]")
ok = True
for fid, eds in E.items():
    for t in eds:
        path = DOCS + t["file"][len("docs/"):]
        txt = open(path, encoding="utf-8").read()
        n = txt.count(t["old_text"])
        start = txt.find(t["old_text"])
        l0 = txt[:start].count("\n") + 1 if start >= 0 else -1
        l1 = l0 + t["old_text"].count("\n")
        bad = []
        if n != 1:
            bad.append(f"old_text count {n}")
        if "—" in t["new_text"]:
            bad.append("em dash")
        if TAGS.search(t["new_text"]):
            bad.append("evidence tag")
        if t["old_text"] == t["new_text"]:
            bad.append("no change")
        print(f"{fid:7s} {t['file']:45s} declared {t['line']:10s} found {l0}-{l1} "
              f"{'OK' if not bad else 'FAIL ' + ', '.join(bad)}")
        ok &= not bad
json.dump(E, open(OUT, "w"), indent=1, ensure_ascii=False)
print("ALL OK" if ok else "PROBLEMS")

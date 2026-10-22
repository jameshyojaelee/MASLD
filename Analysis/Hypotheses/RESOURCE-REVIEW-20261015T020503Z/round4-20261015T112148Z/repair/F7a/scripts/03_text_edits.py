"""Repair-round text edits for F7a. Checks every old_text against the CURRENT file:
exactly one occurrence, and no two edits in the same file overlap. Writes text_edits.json."""
import json
import os
import sys

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEG = "docs/manuscript/05_figure_legends.md"
F7 = "docs/manuscript/draft/fig7.md"
PAP = "docs/PAPER.md"
RDM = "figures/main/fig7_model_framework/README.md"
MET = "docs/manuscript/METHODS.md"
M = "−"

E = {}

E["R08-01"] = [
    dict(file=LEG, claim_level=True,
         old_text=("**(B)** Gain over each task's own baseline. The claim is\n"
                   "task-conditional: pretrained sequence models improve on the task-native\n"
                   "baseline for variant-to-regulatory-activity on a synthetic reporter assay,\n"
                   "where the signal is local and sequence-encoded, and do not improve on it for\n"
                   "tissue-level severity, where the signal is diffuse and tissue-level. Neither\n"
                   "half alone is the result. **(C)**"),
         new_text=("**(B)** Gain over each task's own baseline. On a synthetic reporter assay\n"
                   "(GSE281364, 1,033 elements, HepG2 control and palmitate-oleate contexts only),\n"
                   "one of seven pretrained sequence-model rows (six models), HyenaDNA, exceeds a\n"
                   "ridge on allele identity alone after Bonferroni correction across the seven\n"
                   "rows (gain +0.082, 95% block-bootstrap interval +0.027 to +0.139; corrected\n"
                   "lower bound +0.009 from the same bootstrap draws). Sequence-based reporter\n"
                   "baselines (deltaSVM, gkm-SVM, a sequence CNN or transformer, and the source\n"
                   "study's MPRALegNet) were not run, so this is a gain over a sequence-free\n"
                   "control. On the endogenous variant-to-gene |z| task both pretrained embeddings\n"
                   f"score below distance alone (HyenaDNA {M}0.022, Nucleotide Transformer {M}0.029;\n"
                   "intervals exclude zero). No pretrained sequence model was scored on\n"
                   "tissue-level severity: the severity rows are a composition-plus-expression\n"
                   "model and a routed multi-modal system. **(C)**")),
    dict(file=F7, claim_level=False,
         old_text=("The answer is task-conditional, and both halves are required. Where the signal\n"
                   "is local and sequence-encoded, pretrained sequence models improve on the\n"
                   "task-native baseline: on a synthetic reporter assay they exceed a baseline that\n"
                   "uses allele identity alone, by a small margin. Where the signal is diffuse and\n"
                   "tissue-level, they do not. On the endogenous variant-to-gene effect the\n"
                   "embeddings are inert, in that permuting the correspondence between a variant\n"
                   "and its embedding does not degrade the score; that arm also carries its own\n"
                   "scope limit, because its labels are marginal, un-fine-mapped associations, so\n"
                   "it cannot separate a model failure from a target that is not sequence-determined\n"
                   "at the variants tested."),
         new_text=("The sequence-model evidence is narrow. On a synthetic reporter assay in two\n"
                   "HepG2 contexts, one of seven pretrained sequence-model rows (six models),\n"
                   "HyenaDNA, exceeds a ridge on allele identity alone after Bonferroni correction\n"
                   "across the seven rows, by 0.082 Spearman (95% interval 0.027 to 0.139;\n"
                   "corrected lower bound 0.009). Sequence-based reporter baselines were not run,\n"
                   "so the comparison is with a sequence-free control. On the endogenous\n"
                   "variant-to-gene effect both pretrained embeddings score below distance alone\n"
                   "(HyenaDNA -0.022, Nucleotide Transformer -0.029; intervals exclude zero), and\n"
                   "permuting the correspondence between a variant and its embedding does not\n"
                   "degrade the score. That arm carries its own scope limit: its labels are\n"
                   "marginal, un-fine-mapped associations, so it cannot separate a model failure\n"
                   "from a target that is not sequence-determined at the variants tested. No\n"
                   "pretrained sequence model was scored on tissue-level severity. The only\n"
                   "biologically pretrained model scored on severity, BulkFormer, is an expression\n"
                   "model, and all 109 GSE268273 evaluation samples are in its pretraining set.")),
    dict(file=PAP, claim_level=True,
         old_text=("- Whether a pretrained model improves on a task-native baseline depends on the\n"
                   "  task. On a synthetic reporter assay it does, by a small margin over a baseline\n"
                   "  that uses allele identity alone. On the endogenous variant-to-gene effect it\n"
                   "  does not: the embeddings are inert, in that permuting the correspondence\n"
                   "  between a variant and its embedding does not degrade the score. It also does\n"
                   "  not improve tissue-level severity transfer."),
         new_text=("- On a synthetic reporter assay in two HepG2 contexts, one of seven pretrained\n"
                   "  sequence-model rows (six models), HyenaDNA, exceeds a ridge on allele identity\n"
                   "  alone after Bonferroni correction across the seven rows; sequence-based\n"
                   "  reporter baselines were not run. On the endogenous variant-to-gene effect both\n"
                   "  pretrained embeddings score below distance alone, and permuting the\n"
                   "  correspondence between a variant and its embedding does not degrade the score.\n"
                   "  No pretrained sequence model was scored on tissue-level severity.")),
    dict(file=PAP, claim_level=True,
         old_text=("| Pretrained models help where the signal is local and sequence-encoded and do "
                   "not help where it is diffuse and tissue-level; the difference is measurable, and "
                   "neither result is readable without exposure control |"),
         new_text=("| Of six pretrained sequence models scored in seven rows, one (HyenaDNA) beats a "
                   "sequence-free control on a HepG2 reporter assay after correction; pretrained "
                   "embeddings score below distance on endogenous eQTL effects; no sequence model was "
                   "scored on severity; no score is readable without exposure control |")),
]

E["R08-04"] = [
    dict(file=LEG, claim_level=False,
         old_text=("**(E)** Verdict vocabulary and refusal thresholds,\n"
                   "including the batch sizes below which the framework declines to answer. This\n"
                   "characterizes the framework's own limits and is not a performance claim for any\n"
                   "clinical or diagnostic use."),
         new_text=("**(E)** Verdict vocabulary and refusal thresholds,\n"
                   "including the batch sizes below which the framework declines to answer. Panel E\n"
                   "characterizes the framework's own limits and is not a performance claim for any\n"
                   "clinical or diagnostic use. Each statistic is oriented so that positive is the\n"
                   "prespecified direction. The treatment-response row uses the GSE83452 diet arm\n"
                   "(6 resolved vs 16 persisted), the only arm with more than one persister. There\n"
                   "the point estimate lies opposite to the direction prespecified for the pooled\n"
                   f"contrast, at {M}0.43 of the diet arm's own minimum detectable effect; its 95%\n"
                   "interval spans zero (permutation p 0.23), and the source analysis records the\n"
                   "result as negative and underpowered.")),
]

E["R08-05"] = [
    dict(file=LEG, claim_level=False,
         old_text="Sign statements from the unweighted mean alone are not\nreported.",
         new_text=("Sign statements from the unweighted mean alone are not\nreported. Geneformer V2, "
                   "with no clean held-out study, is drawn hollow gray in\nthe gain forest: its gains "
                   "(+0.049 linear, +0.017 two-layer) are uninterpretable\nwhatever their intervals.")),
]

E["R08-06"] = [
    dict(file=LEG, claim_level=False,
         old_text="**(F)** A liver-specific encoder against the field",
         new_text=("Signed eQTL direction keeps its prespecified verdict, tested negative: in every\n"
                   "|z| stratum the headroom between the label-reliability ceiling and measured\n"
                   "chance exceeds twice the test's resolution. If each gene and LD block holds at\n"
                   "most one causal variant, a sequence model could gain at most 0.006 AUC above chance on\n"
                   "all rows, below the 0.012 resolution, so the test cannot separate a model\n"
                   "failure from a target that is not sequence-determined at these variants.\n"
                   "**(F)** A liver-specific encoder against the field")),
]

E["R09-01"] = [
    dict(file=PAP, claim_level=True,
         old_text=("8. What a predictive model can learn from these maps depends on the task, and\n"
                   "   the difference is measurable. Pretrained models improve on a task-native\n"
                   "   baseline where the signal is local and sequence-encoded, and do not where it\n"
                   "   is diffuse and tissue-level. Neither result is interpretable without knowing"),
         new_text=("8. What a predictive model can learn from these maps depends on the task. Of six\n"
                   "   pretrained sequence models scored in seven rows, one (HyenaDNA) beats a\n"
                   "   sequence-free allele-identity ridge on a HepG2 reporter assay after\n"
                   "   correction; sequence-based reporter baselines were not run, pretrained\n"
                   "   embeddings score below distance alone on endogenous eQTL effects, and no\n"
                   "   sequence model was scored on tissue-level severity. No score is\n"
                   "   interpretable without knowing")),
    dict(file=RDM, claim_level=False,
         old_text=("The claim is task-conditional and it is measured, not asserted. Pretrained\n"
                   "models improve on a task-native baseline for variant-to-regulatory-activity,\n"
                   "where the signal is local and sequence-encoded. They do not improve on it for\n"
                   "tissue-level severity transfer. Both halves are required; neither alone is the\n"
                   "result."),
         new_text=("On a HepG2 reporter assay, one of six pretrained sequence models scored in\n"
                   "seven rows (HyenaDNA) beats a ridge on allele identity alone after Bonferroni\n"
                   "correction across the seven rows. Sequence-based reporter baselines (deltaSVM,\n"
                   "gkm-SVM, sequence CNN/transformer, MPRALegNet) were not run. Pretrained\n"
                   "embeddings score below distance alone on endogenous eQTL |z|. No pretrained\n"
                   "sequence model was scored on tissue-level severity.")),
]

E["R09-03"] = [
    dict(file=MET, claim_level=False,
         old_text=("spaces and are not a model-superiority comparison. Development-only extensions\n"
                   "remain separate from the working Figure 7 claims."),
         new_text=("spaces and are not a model-superiority comparison. The Figure 7B reporter-assay\n"
                   "rows score 1,033 GSE281364 elements in the HepG2 control and palmitate-oleate\n"
                   "contexts only (tanh of the mean Fisher-z Spearman over the two contexts; 10,000\n"
                   "resamples of 239 long-range blocks). LX-2 control and TGF-beta outcomes are in\n"
                   "the same outcome table but were not scored. The comparator is a ridge on allele\n"
                   "identity; deltaSVM, gkm-SVM, sequence CNN or transformer baselines and\n"
                   "MPRALegNet were not run. The Bonferroni family is the seven displayed\n"
                   "reporter-assay rows, and the corrected lower bound is the 0.025/7 quantile of\n"
                   "the same paired bootstrap draws. HyenaDNA's corrected lower bound stays above\n"
                   "zero (+0.006) when the family is widened to all eleven fitted pretrained heads.\n"
                   "Development-only extensions\nremain separate from the working Figure 7 claims.")),
]

E["R08-11"] = [
    dict(file=LEG, claim_level=False,
         old_text=("A nonlinear model with no biological\n"
                   f"pretraining, evaluated on the same folds, differs from ridge by {M}0.046 to\n"
                   "+0.001 Spearman across six comparisons against a minimum detectable effect of\n"
                   "0.067; none exceeds it, and an increment below the minimum detectable effect is\n"
                   "not evidence of no effect."),
         new_text=("A nonlinear model with no biological\n"
                   "pretraining, evaluated on the same folds and the same 30 or 100 principal\n"
                   f"components, differs from ridge by {M}0.014 to +0.001 Spearman across four\n"
                   "comparisons in 109 GSE268273 participants, against a minimum detectable effect\n"
                   "of 0.067; none exceeds it, and an increment below the minimum detectable effect\n"
                   "is not evidence of no effect.")),
    dict(file=F7, claim_level=False,
         old_text=("nonlinear model carrying no biological pretraining, evaluated on the same\n"
                   "folds and the same features, differs from ridge by -0.046 to +0.001 across six\n"
                   "comparisons against a minimum detectable effect of 0.067. None exceeds it."),
         new_text=("nonlinear model carrying no biological pretraining, evaluated on the same\n"
                   "folds and the same 30 or 100 principal components, differs from ridge by -0.014\n"
                   "to +0.001 across four comparisons in 109 GSE268273 participants, against a\n"
                   "minimum detectable effect of 0.067. None exceeds it.")),
]

E["R08-12"] = [
    dict(file=LEG, claim_level=False,
         old_text=("**(D)** Measured bounds and join structure. A bound below the\n"
                   "improvement declared worthwhile in advance is shown with both the bound and the\n"
                   "threshold visible, so an excluded family is legible as excluded on these data\n"
                   "rather than as a failed attempt."),
         new_text=("**(D)** Measured bounds and join structure. A bound below the\n"
                   "improvement declared worthwhile in advance is shown with both the bound and the\n"
                   "threshold visible, so an excluded family is legible as excluded on these data\n"
                   "rather than as a failed attempt. The severity-transfer block shows the released\n"
                   "instrument's Spearman, with its 95% interval, beside each cohort's own\n"
                   "tie-structure ceiling in the published participant records of its three\n"
                   "transfer cohorts (GSE268273 109, GSE213621 299, GSE276114 177).")),
    dict(file=LEG, claim_level=False,
         old_text=("The severity instrument these panels bound transfers out of cohort to 585\n"
                   "participants in three cohorts it never trained on, reaching 0.68 to 0.79 of the\n"
                   "tie ceiling each cohort's own fibrosis coding allows,"),
         new_text=("The severity instrument transfers out of cohort to 585 participant records in\n"
                   "three cohorts it never trained on (panel D), reaching 0.68 to 0.79 of the tie\n"
                   "ceiling each cohort's own stage coding allows,")),
]

# ---- checks against the current files
cache, spans, bad = {}, {}, []
for fid, edits in E.items():
    for e in edits:
        s = cache.setdefault(e["file"], open(os.path.join(ROOT, e["file"])).read())
        n = s.count(e["old_text"])
        if n != 1:
            bad.append((fid, e["file"], n))
            continue
        a = s.find(e["old_text"])
        b = a + len(e["old_text"])
        l0, l1 = s[:a].count("\n") + 1, s[:b].count("\n") + 1
        e["line"] = f"{l0}-{l1}" if l1 > l0 else f"{l0}"
        spans.setdefault(e["file"], []).append((a, b, fid))
        for t in ("—", "[Inference]", "[Unverified]", "[Speculation]", "[Estimate]"):
            if t in e["new_text"]:
                bad.append((fid, e["file"], f"forbidden token {t!r}"))
for f, sp in spans.items():
    sp.sort()
    for (a1, b1, f1), (a2, b2, f2) in zip(sp, sp[1:]):
        if a2 < b1:
            bad.append(("overlap", f, f1, f2))
print("problems:", bad)
for fid, edits in E.items():
    for e in edits:
        print(fid, e["file"], e.get("line"))
json.dump(E, open(os.path.join(RD, "exec", "text_edits.json"), "w"), indent=1, ensure_ascii=False)
sys.exit(1 if bad else 0)

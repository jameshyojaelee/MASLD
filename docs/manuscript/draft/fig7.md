## Model scores are reported against task-native baselines and each model's training exposure

For a Resource of this kind, a practical question is how much of the mapped
biology a predictive model can recover, and how much of an apparent recovery
reflects data the model saw during its own training. To separate the two, we
scored each model family against its task's own baseline and reported every
score beside a per-cohort record of that model's training exposure (Figure 7A).
We resolved exposure for each study from the model's own audit and development
crosswalk, not from a publication. Where the record was unresolved, we reported
the score as uninterpretable, not as evidence of skill. Of the six
encoders in Figure 7C, Geneformer V2 [ref: Geneformer V2 model paper] had no qualified
clean held-out study, so it has no clean-only comparison.

We first examined broad cell-class mapping in seven held-out single-cell
studies. For the exact TranscriptFormer TF-Sapiens checkpoint [ref:
TranscriptFormer model paper], the complete 644-dataset source inventory
identified five clean studies and two training-seen studies. With the two-layer
head, its macro F1 across all seven studies was 0.9389, a value that includes
both training-seen studies (Figure 7C). On its five clean studies, the same head
scored 0.963 (Figure 7C). Figures 7C and 7F draw each encoder's score and gain
twice, once on all seven held-out studies and once on that encoder's clean
studies. This task maps broad cell classes, not disease stages.

The study-level deltas are means over the held-out studies, so small studies
carry a large share of each mean. One held-out study contributes 2 donors and
149 cells, 0.3% of the evaluation set, but receives one seventh of the weight.
When we weighted studies by donors instead, two of twelve model-by-head blocks
changed the sign of their delta (Figure 7F). We therefore report both means and
make no sign claim that rests on the unweighted mean.

We next tested pretrained sequence models on two variant-level tasks (Figure
7B). The first was a synthetic reporter assay in two HepG2 contexts [ref:
GSE281364 reporter assay source study]. Of seven pretrained sequence-model rows
(six models), one, HyenaDNA [ref: HyenaDNA model paper], exceeded a ridge on
allele identity alone after Bonferroni correction across the seven rows. Its
gain was 0.082 Spearman (95% interval 0.027 to 0.139; Bonferroni-corrected
lower bound 0.009) (Figure 7B). We did not run sequence-based reporter
baselines, so this gain is over a sequence-free control. The second task was
the strength (|z|) of endogenous variant-to-gene associations. Adding either
pretrained embedding to distance lowered the held-out Spearman correlation below
that of distance alone (HyenaDNA -0.022, Nucleotide Transformer [ref:
Nucleotide Transformer model paper] -0.029; intervals exclude zero). Permuting
which embedding went with which variant did not lower the score (Figure 7B).
The labels in this task are marginal associations that were not fine-mapped.
The task therefore cannot separate a model failure from a target that is not
sequence-determined at the variants tested. Taken together, the evidence that
pretrained sequence models add skill on these tasks rests on one reporter-assay
row, measured against a sequence-free control.

No pretrained sequence model was scored on tissue-level severity. We did score
BulkFormer [ref: BulkFormer model paper], a biologically pretrained expression
model, on severity, but all 109 GSE268273 evaluation samples are in its
pretraining set. To evaluate severity, we used the released linear instrument.
It was fitted to 521 sample records from five cohorts (GSE130970, GSE135251,
GSE162694, GSE174478 and GSE240729), and its expression
transform was chosen on GSE268273 [ref: GSE268273 source study]. This fit is
separate from the 774-unit differential-expression fit used for the Resource
bulk analysis (Supplementary Note 14). Only the first three of these cohorts are
in the Resource bulk fit. The model's training pool records no donor key, so
its unit is the sample record, not the genetic-identity cluster of the Resource
bulk fit. In the inferential units drawn in Figure
7D, the Spearman correlation with fibrosis stage, as each cohort codes it, was
0.6631 in 109 GSE268273 participant records (source-declared identity), 0.6680 in 263
GSE213621 genetic identities [ref: GSE213621 source study] and 0.7860 in 81
GSE276114 participants with MASLD [ref: GSE276114 source study]. GSE268273 is
scored on Kleiner stage. GSE213621 and GSE276114 are scored on their own
three-level fibrosis codings. The GSE213621
count excludes 11 identities with conflicting stage labels. In 82 GSE276114
participants with chronic viral hepatitis, the correlation was 0.5443 (Figure
7D). On the published participant records, 299 in GSE213621 and 177 in
GSE276114, the raw Spearman correlations were 0.6382 and 0.6739. The 367-record
GSE213621 analysis that includes controls is a sensitivity analysis in the same
cohort, not another external cohort.

We then compared the instrument with a ridge fitted inside each cohort (native
ridge). In the two additional cohorts, GSE213621 (299 records) and GSE276114
(177 records), the recorded paired differences, instrument minus native ridge,
were -0.040 (95% CI -0.099 to +0.017) and -0.011 (-0.079 to +0.054). These
native fits use expression scaling from the evaluation cohort, so they are
transductive. Their conditional intervals do not establish equivalence. In the
reference cohort, GSE268273, the paired difference was +0.054 (recorded 95%
interval -0.025 to +0.137; minimum detectable effect 0.115). These comparisons
do not establish that any model class is generally superior. We also tested a
nonlinear model with no biological pretraining on the same folds and the same
30 or 100 principal components. Across four comparisons in 109 GSE268273
participants, it differed from ridge by -0.014 to +0.001, and none of the four
differences exceeded the minimum detectable effect of 0.067. An increment below
a minimum detectable effect is not evidence of no effect, so we report it as a
bound, not as a negative result. Across Figure 7, each statistic is placed
against its own null and threshold, and an underpowered arm is labelled
untestable, not negative (Figure 7E).

[Figure 7 model results to be added]

We next asked how much the recorded histology limits interpretation. The
current data do not calibrate its measurement error. Within-participant stage
SDs from source-reported repeat biopsies were 0.4655 in PXD051911 [ref:
PXD051911 source study] and 0.6433 in the 58 source-declared GSE193066 repeat-biopsy pairs (every participant with a second biopsy) [ref:
GSE193066 source study]. In GSE193066, 28 of the 58 pairs agreed exactly.
GSE193066 identities remain unqualified, so the 0.6433 estimate and its score
ICC cannot establish primary within-person reliability or label-error
calibration. The 58-pair score ICC also requires scoring outside the
instrument's released acceptance bounds, because only three pairs have both
biopsies inside those bounds. We therefore treat these values as inputs to
sensitivity scenarios, not as an empirically qualified error range or a
calibrated ceiling. Rank agreement is not calibration to absolute stage. Serial
biopsies can also reflect true change, sampling and shared technical effects.
Interpretation also depends on what each model already saw.

The next discriminating experiment is a within-person one. Cross-sectional
transfer cannot distinguish an instrument that tracks severity from one that
tracks between-participant differences correlated with severity. Paired
biopsies from the same participant can make this distinction. We specify that
test in advance against the measurement bound above, rather than reading the
result off after the fact.

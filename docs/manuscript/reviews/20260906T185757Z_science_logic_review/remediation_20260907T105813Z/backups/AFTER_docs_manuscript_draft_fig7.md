## What predictive models learn from these maps, and what bounds them

> Working Figure 7 prose. Panels are candidates until their numbers are
> re-derived by someone who did not produce them. The lead framing of the paper
> is a pending decision; this section states only what is measured.

A resource of this kind invites a question it cannot answer by inspection:
how much of the biology it maps is recoverable by a model, and how much of an
apparent recovery is the model recognizing data it has already seen. We
therefore evaluated model families against each task's own baseline and
reported every score beside a per-cohort record of that model's training
exposure. Exposure is resolved per study from each model's audit and
development crosswalk rather than transcribed from a publication. Where the
record is unresolved, the score is reported as uninterpretable rather than as
evidence of skill; two pretrained single-cell models have no clean held-out
study at all and therefore report no comparison, which is an outcome of the
exposure rule rather than a missing measurement.

The answer is task-conditional, and both halves are required. Where the signal
is local and sequence-encoded, pretrained sequence models improve on the
task-native baseline: on a synthetic reporter assay they exceed a baseline that
uses allele identity alone, by a small margin. Where the signal is diffuse and
tissue-level, they do not. On the endogenous variant-to-gene effect the
embeddings are inert, in that permuting the correspondence between a variant
and its embedding does not degrade the score; that arm also carries its own
scope limit, because its labels are marginal, un-fine-mapped associations, so
it cannot separate a model failure from a target that is not sequence-determined
at the variants tested.

For tissue-level severity the limit is not the model family. In every cohort,
the instrument matches a ridge fitted inside that cohort rather than beating
one. In the two cohorts where the comparison is paired and
interval-estimated the differences are -0.040 (95% CI -0.099 to +0.017) and
-0.011 (-0.079 to +0.054); in the reference cohort the instrument sits +0.054
above the in-cohort ridge. The registered primary arm is the weakest of the
three, and the correct statement is that the instrument matches an in-cohort
baseline, not that it beats one. A
nonlinear model carrying no biological pretraining, evaluated on the same
folds and the same features, differs from ridge by -0.046 to +0.001 across six
comparisons against a minimum detectable effect of 0.067. None exceeds it. An
increment below a minimum detectable effect is not evidence of no effect, and
we report it as a bound rather than as a negative.

One panel-level caution belongs with these comparisons rather than in a
methods note. Study-level deltas across held-out single-cell studies are
averaged over studies, and one held-out study contributes 2 donors and 149
cells, 0.3% of the evaluation set, while receiving a seventh of that average.
Two of twelve model-by-head blocks change the sign of their delta when studies
are weighted by donors instead, so we report both averages and make no sign
claim that rests on the unweighted one.

The bound that organizes these results is the histology label. Repeat biopsies
in two independent deposits give a within-participant Kleiner standard
deviation of 0.4655 to 0.6433 stages, with 28 of 58 pairs agreeing exactly in
the first deposit; the same error variance yields reliabilities that differ
threefold between deposits because the two populations differ in stage
variance. That is a range, not a transportable constant, and it means the
ceiling on any severity score is a property of the cohort being scored. Read
together, these panels say that a diffuse, linear tissue signal measured
against a noisy label leaves little for a larger model to add, and that the
statement is only interpretable because every score is paired with what its
model already saw.

The next discriminating experiment is a within-person one. Cross-sectional
transfer cannot distinguish an instrument that tracks severity from one that
tracks between-participant differences correlated with it; paired biopsies from
the same participant can, and that test is specified in advance against the
measurement bound above rather than read off after the fact.

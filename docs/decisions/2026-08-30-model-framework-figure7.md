# Seventh figure: model framework and transfer limits

Date: 2026-08-30

Decision: the Resource manuscript adds a seventh main figure reporting what
predictive models learn from the Resource maps. The paper stays one paper. There
is no separate methods or benchmark manuscript.

What changed, as one coordinated edit:

- `config/resource_paper_scope.config.json` gains a seventh figure role,
  `model_framework_and_transfer_limits`; the goal paragraph now names five
  linked dimensions instead of four; the title now names the model measurement.
- `scripts/manuscript/validate_resource_scope.py` expects seven figure roles
  instead of six, requires the `## Seven-figure story` heading, requires
  `7. Model framework and transfer limits` in the ordered spine, and requires
  the two `figures/main/fig7_model_framework/` files.
- `docs/PAPER.md` mirrors the goal and title word for word, adds the seventh
  row to the figure table, and adds five lines to "Claims we can make".
- `figures/main/fig7_model_framework/` is created with the same two files every
  other figure directory carries.

What did **not** change, deliberately:

- `central_scientific_contribution` stays
  `inherited_genetic_variation_and_multicellular_disease_remodeling`. Figure 7
  measures what models learn from the maps; it does not replace what the paper
  is about.
- The two lines in "Claims we cannot make" that rule out a universal target rank
  and prospective biomarker performance stay exactly as they are. Figure 7 is
  written to sit inside them. It reports measured transfer and measured limits,
  never clinical or diagnostic performance.

Why this figure exists. The Resource assembles maps in Figures 1 to 6, and the
obvious question is whether a model can learn anything from them that the maps
do not already state. The answer is conditional and the condition is measurable.
On variant-to-regulatory-activity, pretrained sequence models beat a trivial
allele-identity baseline. On tissue-level severity transfer they do not, and an
upper bound over the whole family of combination strategies falls well below the
improvement declared worthwhile beforehand, so the family is excluded on these
data rather than merely untried. The structural reason is visible in the
Resource itself: no donor in this project has three molecular layers measured on
them, so there was never much to combine.

The part with no published equivalent we could find is the exposure control.
Every model score is reported beside a per-cohort record of whether that cohort
was in the model's own training data, checked against the published training
tables with their hashes. That makes it possible to hold a model fixed and vary
only what it had already seen.

Revoking this: remove the seventh role from the config, restore six in the
validator, restore the `## Six-figure story` heading and drop the seventh table
row and the five claim lines from `docs/PAPER.md`, and delete
`figures/main/fig7_model_framework/`. The validator fails closed on any partial
revert, so it will tell you if one of the five is missed.

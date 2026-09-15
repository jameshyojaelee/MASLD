> # ⛔ SUPERSEDED -- DO NOT CITE SECTION 4d
>
> **Section 4d of this card is WRONG, and so is part of section 4b. USE
> `masld-liver-protein-severity-v1.3` INSTEAD** -- that is the current card, and it is the
> only one you need. Do not quote any number from sections 4b or 4d below.
>
> Two things v1.1 got wrong: section 4d proposed a mechanism for an inversion that **the
> deployed scorer does not exhibit at all**, and section 4b's replacement null was
> **over-inverted**, making its bar easier to clear rather than harder -- the opposite of
> what this card claims a few paragraphs down. Section 4c-bis's liveness measurements stand,
> but its reading of the spread across constructions as a real dependence does not: v1.2
> shows it is UNDERPOWERED.
>
> **What is wrong.** `score.py:250` computes `severity = intercept + z @ w` where
> `intercept` is a SINGLE CONSTANT, exactly the training mean 1.0689655172413792. A
> constant added to every sample cannot change a Spearman correlation at all, so the
> deployed ordering is `z @ w` and nothing else. The coverage curves in sections 4d and 5
> instead used a PER-SAMPLE leave-one-out intercept `(n*ybar - y_i)/(n-1)`, whose Spearman
> against Kleiner is exactly -1.0000. At 2% coverage that term's spread (0.01389) exceeds
> the restricted dot product's (0.01167), so the curve measures the intercept.
>
> **A per-fold intercept is CORRECT for the leave-one-out headline** -- it is the honest
> out-of-fold predictor, and at full coverage it contributes about a thirty-sixth of the
> dot term's spread, so **0.7291 is untouched**. It is **WRONG for a deployment curve**,
> because the shipped scorer adds a fixed constant. The two uses were conflated. Since a
> constant cannot move a Spearman, the deployed curve needs no intercept term at all.
>
> **v1 blamed the weight vector. v1.1 blamed the leave-one-out restriction geometry. Both
> were explaining a number that does not describe `score.py`.** The deployed scorer does
> not appear to invert.
>
> **Section 5's 30% threshold is NOT withdrawn.** The spurious intercept DEPRESSES the
> deployed column at every coverage and most at low coverage, so removing it raises the
> curve everywhere and section 5's conditions would hold below 30%. The threshold is
> therefore CONSERVATIVE rather than permissive and no user of this release is exposed by
> it. [That inference was later confirmed by measurement: the deployed curve is positive at
> every coverage tested and the rule's margin is 0.25 -- see v1.3.] What needs restating is the
> threshold's JUSTIFICATION and section 5's "INVERTS" language, not the number.
>
> **Two further corrections in v1.2.** Section 4b's null is described as
> "geometry-preserving"; the accurate description is **marginal-preserving and
> covariance-destroying** -- permuting each protein across participants keeps per-protein
> marginals and destroys the between-protein covariance, collapsing the top Gram
> eigenvalue 76,013 -> 6,471 and raising effective df 27.4 -> 35.2. The null therefore
> shrinks LESS than a covariance-matched one and is less negative than it should be, which
> makes section 4b's bar CONSERVATIVE; the section's conclusion is unaffected. Section
> 4c-bis's construction-dependence table omits a fifth construction, the per-fold-lambda
> plain arm at +0.3316 [+0.0870, +0.5948], which also excludes zero; its omission
> UNDERSTATED that section's own finding.
>
> **What still stands:** the headline 0.7291, every weight, `score.py`, `expected.json`,
> the 30% threshold, section 4b's conclusion, and every section 4c-bis point estimate and
> liveness value. v1.1's measured numbers are left exactly as they were computed; this
> notice is added on top and nothing below it has been altered.

# masld-liver-protein-severity-v1.1

**v1.1 corrects three sections of v1's model card. The model is unchanged** -- same weights,
same `score.py`, same 30% coverage threshold, same leave-one-out Spearman
0.7291. Scores from v1 and v1.1 agree to
`0.0e+00`.

v1's sections 4b and 4d measured their nulls by fitting on a permuted label and scoring
against the unpermuted one, which makes the statistic centred at zero by construction at any
lambda. Section 4c-bis described a model that was never fitted. Sections 4b, 4c-bis and 4d
are corrected, **and section 4d's conclusion is reversed**: the inversion at 2% coverage is
the leave-one-out restriction geometry, not the weight vector. `MODEL_CARD.md` opens with the
supersession table.

`release/masld-liver-protein-severity-v1/` is untouched and remains the record of what
shipped.


A fibrosis severity score for bulk liver PROTEOMICS. Give it intensities, get one number
per sample.

**Before anything else: this instrument has NO EXTERNAL VALIDATION.** One cohort, 58
participants, zero held-out cohorts. `MODEL_CARD.md` section 0 leads with it.

## Install

Nothing to build. Python 3.9+, `numpy` and `scipy` (`pandas` only for the smoke test).

    pip install numpy scipy

## Run

    python score.py --abundance my_intensities.tsv --out scores.tsv

`my_intensities.tsv` is proteins-by-samples RAW INTENSITIES: first column protein ids,
header row sample ids, `NA` for missing. NPZ also works (`--abundance x.npz` with arrays
`abundance`, `ids`, optionally `samples`).

Output columns: `sample_id`, `severity`, `axis_coverage`. There is deliberately **no
interval** -- see `MODEL_CARD.md` section 7.

If your identifiers are gene symbols rather than UniProt accessions:

    python score.py --abundance x.tsv --id-type symbol --out scores.tsv

Add `--json report.json` for join diagnostics: how many of the 5,377 axis proteins
matched, how many duplicate rows were summed, how many values were imputed.

## Check it works

    bash test_release.sh

Scores the 58 training proteomes through the released code path only, asserts the scores
reproduce to 1e-8, and asserts the scorer *refuses* five inputs it should refuse.

## What this does, honestly

It **orders** liver proteomes by fibrosis severity, at one point in time, in research.
Leave-one-participant-out Spearman **0.7291** against a tie-structure ceiling of
0.9107 on 58 participants (80% of it), against a single-best-protein baseline of
0.5805. All of that is INTERNAL: there is no second cohort.

It does **not** give you a Kleiner stage, says nothing about one person over time, nothing
about cell fractions, and nothing clinical. It is not a diagnostic.

## The cross-modal result

The RNA instrument's released weights, mapped by gene symbol onto protein abundance and
applied WITHOUT REFITTING, reach Spearman **+0.6128** against Kleiner -- 84% of the
protein-native score, from 19.6% of the RNA axis. A permuted-weight control sits at
-0.0043.

**Read section 6 before quoting that number.** Intensity complexity alone reaches
+0.4772 here, and once it is in the model the cross-modal score's
nested increment is +0.1927, CI95 [-0.0146, +0.4121] -- it includes zero. The protein-native score's increment
does not. At n=58 this cohort cannot separate the transfer from intensity complexity.

Note that 19.6% is below `masld-severity-v1`'s own 30% coverage threshold, so its
`score.py` refuses this matrix; the transfer was run with that threshold overridden.

The nested-gate liveness of that increment was never measured in v1; it is
0.6272 against a >0.99 dead-gate
threshold, so the comparison is a genuine nested one. See section 4c-bis.

## Files

    score.py          the scorer
    weights/          coefficients, frozen constants, the protein axis
    MODEL_CARD.md     read this before using any number
    test_release.sh   smoke test through the released path
    expected.json     the values the smoke test asserts

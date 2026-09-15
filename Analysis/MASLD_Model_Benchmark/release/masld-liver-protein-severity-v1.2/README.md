> # ⛔ SUPERSEDED BY v1.3 -- SECTION 4b QUOTES THE WRONG STATISTIC
>
> **Use `masld-liver-protein-severity-v1.3`.** Every measured number below is correct and
> unchanged in v1.3, and no verdict changes. What is wrong is which statistic section 4b
> puts the weight on.
>
> Section 4b says *"the bar is the row-permutation null: +0.4506 at its maximum over 200
> draws."* **A maximum over N draws is an increasing, unbounded function of N.** Take more
> draws and it rises on its own, with no change to the model, the data or the null -- so the
> claim as written depends on the analyst's draw count rather than on the instrument. Anyone
> rerunning this null at 1,000 or 10,000 draws would get a higher "bar" and could read it as
> a refutation when nothing had changed.
>
> v1.3 restates the same result with statistics that are stable in N: **0 of 200
> row-permutation draws reach the observed 0.7291, so p <= 1/201 = 0.00498**, and the
> standing bar is the null's **97.5th percentile +0.3322**, which the observed exceeds by
> +0.3969. The maximum is kept as a descriptive statistic with its draw count attached and is
> never called the bar. More draws can only tighten that p.
>
> Sections 4c-bis and 4d are carried into v1.3 unchanged. This card is left exactly as it
> shipped; the notice is added on top and nothing below it has been altered.

# masld-liver-protein-severity-v1.2

**v1.2 corrects three sections of v1's model card. The model is unchanged** -- same weights,
same `score.py`, same 30% coverage threshold, same leave-one-out Spearman
0.7291. Scores from v1 and v1.2 agree to
`0.0e+00`.

v1's sections 4b and 4d measured their nulls by fitting on a permuted label and scoring
against the unpermuted one, which makes the statistic centred at zero by construction at any
lambda. Section 4c-bis described a model that was never fitted.

**v1.1 corrected those and got two of them wrong in new ways**, and is superseded. v1.2
replaces section 4b's null with a covariance-preserving one, reads section 4c-bis's spread as
UNDERPOWERED rather than construction-dependent, and establishes that **the deployed scorer
does not invert at low coverage at all** -- `score.py` adds a constant intercept, which cannot
change a rank correlation. `MODEL_CARD.md` opens with the full supersession table.

`release/masld-liver-protein-severity-v1/` and `-v1.1/` are untouched and remain the record of
what shipped.


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

# masld-liver-protein-severity-v1

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
nested increment is +0.1927, CI95 [-0.0176,
+0.4137] -- it includes zero. The protein-native score's increment
does not. At n=58 this cohort cannot separate the transfer from intensity complexity.

Note that 19.6% is below `masld-severity-v1`'s own 30% coverage threshold, so its
`score.py` refuses this matrix; the transfer was run with that threshold overridden.

## Files

    score.py          the scorer
    weights/          coefficients, frozen constants, the protein axis
    MODEL_CARD.md     read this before using any number
    test_release.sh   smoke test through the released path
    expected.json     the values the smoke test asserts

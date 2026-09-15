# masld-severity-v1

A latent fibrosis severity score for bulk liver RNA-seq. Give it counts, get one
number per sample.

## Install

Nothing to build. You need Python 3.9+, `numpy` and `scipy` (`pandas` only for the
smoke test).

    pip install numpy scipy

## Run

    python score.py --counts my_counts.tsv --out scores.tsv

`my_counts.tsv` is genes-by-samples **counts**: first column gene ids, header row
sample ids. Ensembl stable gene ids, with or without a `.7` version suffix. NPZ
also works (`--counts x.npz` with arrays `counts`, `genes`, optionally `samples`).

Output columns: `sample_id`, `latent_severity`, `severity_vs_batch_median`,
`axis_coverage`, plus `vs_batch_median_lo80/hi80/lo90/hi90` **when you supply at
least 20 samples from the same batch**. Below that, the scores come back
with no interval and `score.py` says why on stderr — an interval at that batch size
would not be calibrated, so none is given.

There is deliberately **no interval on the absolute severity level**. It was
measured uncalibrated on both held-out cohorts, and the score does not recover
where a cohort sits: across the five comparable cohorts the slope of predicted on
true cohort level is +0.0120, against
60% of the within-cohort spread retained. So no absolute
interval can be calibrated by any method. `MODEL_CARD.md` section 7 has the
measurements, including the two artifacts that had to be removed to get that
number right.

Add `--json report.json` for the join diagnostics: how many of the 26,629 axis
genes matched, how many version suffixes were stripped, how many duplicate stable
ids were summed.

## Check it works

    bash test_release.sh

This scores GSE268273 — 109 participants the model never trained on — from raw
files through the released code path only, asserts the out-of-cohort Spearman
0.6631, and asserts that the scorer *refuses*
four inputs it should refuse.

## The measurement-error calculator

    python kleiner.py --demo
    python kleiner.py --stage-sd 1.2 --r 0.62

Independent of the model. Corrects correlations and slopes against the Kleiner
stage for the stage's own read error, sigma_e = 0.6433 stages. `--demo` prints the
worked example in which a reported r = 0.50 implies an impossible r_latent = 1.15.

## What this does, honestly

It **orders** bulk liver RNA-seq samples by fibrosis severity, at one point in
time, in research. Out of cohort it reaches Spearman
0.6631 against a tie-structure ceiling of
0.9523 on 109 participants, which is
70% of the most any score could reach
against that outcome, and it is not worse than a ridge trained inside that cohort
(paired interval [-0.0253, +0.1373], which
includes zero). It does **not** give you a Kleiner stage: the absolute level
carries an uncorrected per-cohort offset, so only the ordering is meaningful. It
says nothing about one person over time — the intervention arm reverses sign 5 of
5 — nothing about cell fractions, and nothing clinical. It is not a diagnostic.
Every one of those limits is a measurement, not caution; `MODEL_CARD.md` gives the
numbers behind each.

## Files

    score.py          the scorer
    kleiner.py        the disattenuation calculator, standalone
    weights/          coefficients, the frozen standardization constants, the gene axis
    MODEL_CARD.md     read this before using any number
    test_release.sh   smoke test through the released path
    expected_ood.json the value the smoke test asserts

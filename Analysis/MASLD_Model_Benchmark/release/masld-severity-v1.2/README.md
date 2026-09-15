# masld-severity-v1.2

A latent fibrosis severity score for bulk liver RNA-seq, plus (v1.2) an activity ordering and
a label-free shift diagnostic. Give it counts, get two numbers per sample and a diagnostic line.

## What changed from v1.1

- **`latent_activity`**, appended as the LAST column so every v1.1 column keeps its position.
  It orders samples by NASH-CRN activity (steatosis + ballooning + lobular inflammation),
  fitted leave-one-cohort-out on 464 training samples from four cohorts and shipped as primal
  weights. Its sealed gate passed at the minimum 2 of 4 cohorts; it carries **no two-axis and
  no subtype claim**. `MODEL_CARD.md` 4.4. Disable with `--activity-weights ''`.
- **A shift diagnostic on stderr and in the JSON**, before any score: Hotelling T² and the
  off-basis residual Q in the training pool's 20-PC space, cohort displacement per PC, and
  ribosomal / mitochondrial fraction z. **It is not a gate.** A refusal bound on the same
  statistics was prespecified, measured, and not shipped: it admitted GSE281797 (0 of 94)
  and refused GSE225740 (73 of 93) for coverage loss. `MODEL_CARD.md` 6d has the table.
  Disable with `--diagnostic-weights ''`.
- `latent_severity`, every refusal, the input contract and the interval are unchanged;
  `test_release.sh` asserts bitwise identity with v1.1 on GSE268273.

## What changed from v1 to v1.1

## What changed from v1, and why v1 is superseded

**v1 could not refuse input it had never been measured on.** The score is a pure
function of the within-sample rank vector over the matched genes, so it is exactly
invariant to sequencing depth and it responds only to the ordering and the tie
structure. Hand v1 a matrix with a different sparsity structure — pseudobulked
single-nucleus counts, a shallow library, a thresholded matrix — and every guard it
had passed while the number it returned was outside anything that had been measured.

v1.1 adds a measured **input contract** on two per-sample statistics: the zero
fraction over the matched axis (bound **25.44%**), and the divergence of the sample's
realized post-quantile marginal from the frozen reference (bound **0.036006**).
Out-of-contract input **raises**, like every other refusal here.

**The bound is tight and will refuse some legitimate input.** It sits 1.10x above the
largest value seen in 630 in-contract samples, and on real single-nucleus donor
pseudobulk it refuses 4 of 46 and 2 of 39 donors. **A refusal means your input is
outside the measured envelope, not that your data are bad.**

    python score.py --counts my_counts.tsv --out scores.tsv --acknowledge-out-of-contract

scores it anyway, records `out_of_contract: true` in the JSON with both observed
statistics and both bounds, and warns on stderr. Those scores have no measured ordering
and no measured error. `MODEL_CARD.md` sections 2, 6b and 6c have the degradation curve
the bounds were measured on, the looser alternative that was computed and rejected, and
where pseudobulk lands.

Three further defects are fixed:

- **`_PAR_Y` identifiers were dropped, not summed.** `VERSION_SUFFIX` is anchored at end,
  so a GENCODE id like `ENSG00000002586.21_PAR_Y` never matched the axis and its counts
  were discarded rather than summed into the X row — the exact mechanism the model card
  advertises, silently not working for the case it names. 23 axis genes are affected in
  principle; **zero** `_PAR_Y` rows exist in any cohort this instrument has scored, so no
  published number moves. `MODEL_CARD.md` 6c has the measurement and both branches of the
  fix.
- the denylist read only the **first 2000** sample ids, so a denylisted accession in
  column 2001 or beyond passed. It now reads every id.
- `MODEL_CARD.md` section 7 announced output columns `latent_lo95`, `latent_hi95` and
  `next_read_*_sigma_e_*` that **do not exist in the code**, in v1 or here. Section 7c
  was and is correct; the opening paragraph was not, and is rewritten.

**Scores on in-contract input are bitwise identical to v1**, and `test_release.sh`
reproduces the same out-of-cohort Spearman 0.6631190836675596. v1 remains on disk,
unmodified, with its own manifest. Nothing in it was overwritten.

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
ids were summed, and where your input sits against the two input-contract bounds
(`zero_fraction_observed_max`, `qn_w1_observed_max` and the bounds they are checked
against, reported on every run whether or not the contract fired).

## Check it works

    bash test_release.sh

This scores GSE268273 — 109 participants the model never trained on — from raw
files through the released code path only, asserts the out-of-cohort Spearman
0.6631, and asserts that the scorer *refuses*
six inputs it should refuse: a symbol-keyed matrix, a matrix carrying 5% of the
axis, a CPM matrix, a denylisted accession, a denylisted sample id at column 2401,
and a 70%-masked matrix that trips the input contract while every other guard
passes. It then asserts that `--acknowledge-out-of-contract` lets that last one
through *and records the acknowledgement in the output*, and that a pseudoautosomal
gene split into an X row and a `_PAR_Y` row is summed back to a bitwise-identical score.

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

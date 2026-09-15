# encoder-benchmark-v1

An exposure-resolved benchmark for single-cell encoders. Give it predictions and
an exposure declaration, get paired deltas against a native baseline — and, for
any model whose pretraining corpus overlaps the held-out data, no number at all.

## Install

Nothing to build. Python 3.10+, `numpy` and `scipy`. `torch` only if you run the
baselines.

    pip install numpy scipy

## Run

    python evaluate.py \
      --predictions reference/predictions.npz \
      --row-contract reference/row_contract.tsv \
      --exposure reference/exposure_declaration.json \
      --reference scvi_liver_latent::two_layer_mlp \
      --out results/

You get `leaderboard.tsv`, `not_comparable.tsv`, `skipped_for_head.tsv` and
`results.json`.

**Input contract.** `row_contract.tsv` needs `row_id`, `donor_id`, `study_id`,
`outer_fold`, `true_class`. `predictions.npz` needs a `row_id` array and one
`n_rows x n_classes` array per model. `exposure_declaration.json` maps each model
to `head_id` and a `per_study_exposure` state for **every** study in the roster.

## What it reports, and what it refuses

The reportable quantity is never a bare score. It is the **paired delta against a
reference you name**, restricted to that model's **clean** held-out studies, with
a 95% interval from a **donor** cluster bootstrap using multiplicities shared
across models so the deltas are paired on the same draws.

Exposure is resolved **per study, not per model** — the two disagree in real
data. Geneformer declares its checkpoint `target_label_unexposed` while being
`encoder_seen` on three of seven held-out studies.

**A model with no clean held-out study gets no number.** Not a number with an
asterisk, not an empty cell. `exposure.py` enforces this by type: the comparable
result cannot be constructed without a clean study, the non-comparable one has no
numeric attribute and raises if you reach for one, and the leaderboard refuses to
rank it. In the shipped reference that withholds **Geneformer and
TranscriptFormer, both heads** — which are two of the three highest raw scorers.

By default a model is compared only against a reference carrying the **same
head**. Comparing a linear-head model to an MLP-head reference measures the head
and calls it the encoder, which is the exact error this benchmark exists to
document. `--head-policy any` lifts it if you mean to.

## Check it works

    bash test_release.sh

Runs the evaluator over the shipped reference through the released path only,
asserts it reproduces the published deltas, asserts the four exposure-withheld
encoders are absent from the leaderboard, and asserts that asking them for a
number **raises**.

## The baselines you have to beat

    python baseline.py --counts my_counts.npz --row-contract my_rows.tsv --out preds/

`hvg_pca` is a 50-dimensional PCA of highly variable genes. It is not a strawman:
it is one of only two blocks with all seven held-out studies clean, and with a
two-layer head it is a **dead tie** with a liver-specific encoder trained on 102
donors (+0.000015, interval [−0.0050, +0.0060]).

`library_shape` is ten library-complexity descriptors and no gene identity at
all. It exists so the benchmark can state its own floor. Here it reaches 0.2456
against a random-scorer floor of 0.1808 and a PCA at 0.9254 — a floor, not a
competitor. On this task. `MODEL_CARD.md` says why that does not generalise.

## What this measures, honestly

Where a clean comparison is possible the margins are small and **the head matters
more than the encoder**. UCE-4L beats the liver-specific encoder by **+0.0063**
with a linear head and **loses by −0.0927** with an MLP — a swing far larger than
any gap between encoders. Every number here is one cohort family, 102 donors,
five broad cell classes. Nothing here is replicated, and none of it is a
statement about clinical use.

## Files

    evaluate.py             the evaluator
    exposure.py             the data model that enforces the refusal
    metrics.py              metric and bootstrap machinery, vendored verbatim
    baseline.py             the baselines, runnable on your own data
    reference/              21 models x 50,000 rows, the row contract, the exposure declaration
    expected_benchmark.json the values test_release.sh asserts
    MODEL_CARD.md           read this before quoting any number
    test_release.sh         smoke test through the released path

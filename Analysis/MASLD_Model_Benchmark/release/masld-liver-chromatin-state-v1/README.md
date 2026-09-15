# masld-liver-chromatin-state-v1

Ten axes of liver H3K27ac chromatin state, predicted from bulk liver RNA counts. Trained with paired
RNA and H3K27ac from the same biopsies; **needs only RNA at use**.

```
python score.py --counts my_counts.tsv --out states.tsv --json report.json
bash test_release.sh          # re-derives the shipped fit and its self-consistency guard
```

`my_counts.tsv` is genes-by-samples **counts** with Ensembl stable gene ids in the first column
(versioned ids and `_PAR_Y` suffixes are handled). You get ten `chromatin_state_k1..k10` columns.

## The one-paragraph version

Bulk liver RNA carries at least ten separable axes of H3K27ac that are **not** explained by
histology. Each was verified out of fold, on held-out participants, after conditioning on all four
NASH-CRN grades, sex and 14 RNA library-complexity descriptors. The strongest two are hepatocyte
synthetic function (out-of-fold partial ρ **+0.770**) and hepatocyte identity/zonation versus an
immediate-early stress program (**+0.666**). All ten have a lower bound above zero.

## What the scores are, and are not

**Are:** a within-batch ordering of chromatin state, in units of the training cohort's component
standard deviation, with a per-component out-of-fold residual sd shipped as the interval basis.

**Are not:** a stage, a subtype, a diagnosis, a measurement of any individual enhancer, or anything
longitudinal. The absolute level does not transport between batches. `MODEL_CARD.md` section 6 lists
every limit in force; they are not boilerplate.

## Why it is a linear model

Four alternatives were tested against this exact target under sealed prespecifications, and all four
lost: a paired teacher–student (−0.057), MOFA+ (−1.74, worse than the training mean), BulkFormer
embeddings (pretrained loses to its own random-init twin, 15/15 comparisons), and Corgi — the closest
published RNA-conditioned epigenome model — both zero-shot (+0.0003) and **fine-tuned on this paired
cohort** (+0.0001), against this model's +0.119 on the same regions. Tuned kernel-ridge, gradient
boosting and MLP heads never beat ridge either. The choice is measured, not assumed;
`MODEL_CARD.md` section 5 has the table.

## Files

| file | what |
|---|---|
| `score.py` | the only entry point; every refusal raises rather than returning a number |
| `weights/chromatin_state_v1.npz` | gene axis, frozen PCA basis (81 dims), ten ridge heads, PLS chromatin loadings for low-rank profile reconstruction, per-component out-of-fold residual sd |
| `weights/fit_report.json` | out-of-fold accuracy per component, the fit's hyperparameters, and the representation choice with its reason |
| `MODEL_CARD.md` | accuracy, contract, what lost to it, and the scope limits |
| `test_release.sh` | end-to-end run plus the self-consistency guard |

A dense raw-gene → 96,460-region profile head is deliberately **not** shipped: it is ~32 GB, and the
same reconstruction is available in low-rank form through `pls_V`.

## Provenance and terms

Trained on GSE267145 (99 participants). Weight release is permitted under
`docs/decisions/2026-09-07-derivative-weight-release.md`. **Source data are not redistributed** —
this release ships weights, code and aggregate metrics only.

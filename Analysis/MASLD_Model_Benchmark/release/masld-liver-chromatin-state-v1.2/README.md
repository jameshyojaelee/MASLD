# masld-liver-chromatin-state-v1.2

Predicts, from a bulk liver RNA counts matrix, ten directions of liver H3K27ac covariance learned in one
paired cohort (GSE267145, 99 participants) plus one **fixed steatosis-associated chromatin axis**.
Trained with paired RNA and H3K27ac; **needs only RNA at use**.

```
python score.py --counts my_counts.tsv --out states.tsv --json report.json
python score.py --counts my_counts.tsv --out states.tsv --profile profile.npz      # v1.2: per-region H3K27ac profile
bash test_release.sh
```

`my_counts.tsv`: genes-by-samples **counts**, Ensembl stable ids in the first column (versions and
`_PAR_Y` handled). Output: `chromatin_state_k1..k10`, `steatosis_chromatin_axis`, `axis_coverage`.

## New in v1.2

Primary profile form: `rrr_offset_cis` (see card section 3b for how it was chosen). `--profile profile.npz` writes, per sample, the predicted residualised H3K27ac over the 96,460 training regions
(float64) with the 4,847-region reliably-predictable mask. In-cohort out-of-fold skill **+0.1429** of the residual variance
(+0.024 over the best global model, every fold); the gene block is **local** (+0.036 over ten random genes); it **transfers**
to 39 multiome donors' ATAC (mean per-region ρ +0.205 on the reliable set vs a within-well null of +0.054, p 0.000999;
shuffled-fit control at its null). Needs a batch of ≥ 10 samples under the default marginal transport. Card section 3b.

## Read this before using a number

This is a **development-grade** model and its card says so in every section. What changed from v1:

- The ten heads are bitwise identical to v1. Their accuracies are now the out-of-fold numbers of the
  **exact shipped recipe**: k1 **+0.596**, k2 **+0.612** (v1 printed 0.666 / 0.770, copied from a
  different estimator). Whole-matrix reconstruction of the shipped form is **+0.045**, not +0.119.
- k2 is an **activity axis** (ρ +0.55 with steatosis, +0.50 lobular inflammation), not "orthogonal to
  histology"; k1, k3, k5 carry sex in the training cohort. Only k1–k3 keep a stable identity across folds.
- **Outside GSE267145 no axis is a reproducible donor property**: on 58 repeat-biopsy pairs k1 has ICC
  +0.008 and no axis beats a random direction in the same RNA space; k3, k5, k10 are sex axes in three
  external cohorts; k2, k7, k8 track fibrosis grade externally; k7, k8 track age.
- New: the **steatosis-associated head** (out-of-fold partial **+0.636** [+0.47, +0.76] given all four
  grades, sex and 14 RNA descriptors), the one component with a motivated biological target.

What a user can legitimately do with it: order samples **within one batch** by RNA-predicted position on
each direction, as a development tool, knowing that ten unsupervised RNA principal components do nearly
as well for k1–k7. What they cannot do is call the outputs liver chromatin states, compare across
batches, or read them as anything about an individual regulatory element.

## Files

| file | what |
|---|---|
| `score.py` | the only entry point; every refusal raises rather than returning a number |
| `weights/chromatin_state_v1_2.npz` | every v1.1 array (bitwise) + `prof_*` profile head (three forms, transforms, reliable mask, per-fold skills) |
| `weights/predictable_regions.bed` | the 4,847 reliably predictable regions with per-fold out-of-fold skills |
| `weights/fit_report.json` | every number on the card, each traceable to a lane JSON key and sha256 |
| `MODEL_CARD.md` | accuracy, what the axes are, external evaluation, limits, and the list of v1 numbers this release kills |
| `test_release.sh` | end-to-end run; asserts the ten columns are bitwise v1 and the reported numbers are the lane's, not the in-sample fit |

Source data are not redistributed. Weight release is permitted under
`docs/decisions/2026-09-07-derivative-weight-release.md`.

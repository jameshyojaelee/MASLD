# masld-liver-chromatin-state-v1.1

Predicts, from a bulk liver RNA counts matrix, ten directions of liver H3K27ac covariance learned in one
paired cohort (GSE267145, 99 participants) plus one **fixed steatosis-associated chromatin axis**.
Trained with paired RNA and H3K27ac; **needs only RNA at use**.

```
python score.py --counts my_counts.tsv --out states.tsv --json report.json
bash test_release.sh
```

`my_counts.tsv`: genes-by-samples **counts**, Ensembl stable ids in the first column (versions and
`_PAR_Y` handled). Output: `chromatin_state_k1..k10`, `steatosis_chromatin_axis`, `axis_coverage`.

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
| `weights/chromatin_state_v1_1.npz` | every v1 array (bitwise) + `ste_*` steatosis head |
| `weights/fit_report.json` | every number on the card, each traceable to a lane JSON key and sha256 |
| `MODEL_CARD.md` | accuracy, what the axes are, external evaluation, limits, and the list of v1 numbers this release kills |
| `test_release.sh` | end-to-end run; asserts the ten columns are bitwise v1 and the reported numbers are the lane's, not the in-sample fit |

Source data are not redistributed. Weight release is permitted under
`docs/decisions/2026-09-07-derivative-weight-release.md`.

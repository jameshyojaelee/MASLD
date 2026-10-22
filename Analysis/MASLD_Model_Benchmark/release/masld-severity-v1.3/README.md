# masld-severity-v1.3 (masld-pairwise-conformal-v1): candidate, not adopted

**Status: candidate, not adopted.** Installed 2026-10-17 by user decision ("yes to all") from the B1 candidate package
`/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-program-20261015T223347Z/trackB/B1/package-20261016T013445Z/`. Result: `/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-program-20261015T223347Z/trackB/B1/RESULT.md`.
It does not supersede masld-severity-v1.2 and is not a paper release. It has a MODELS.md row under
"Candidates 2026-10 (not adopted)".

Read `MODEL_CARD.md` before using any output.

```
python score.py --counts counts.tsv --out scores.tsv --pairs pairs.tsv [--batches batches.tsv] [--json report.json]
bash test_release.sh
```

| File | Content |
|---|---|
| `score.py` | counts -> tie-averaged plane -> refit score -> within-cohort pair calls at the sealed 0.90 q |
| `weights/b1_pairwise_refit_weights.npz` | frozen axis, mu, sd, qn_reference; per plane: axis positions, w, intercept, lambda, q(0.90) |
| `weights/plane_<name>_genes.tsv` | the exact gene set of each plane |
| `thresholds.json` | the three sealed 0.90 q values, the rule, and the no-guarantee statement |
| `expected.json` | GSE213621 and GSE281797 coverage, call rate, accuracy and pair counts from B1 `results.json` |
| `test_fixtures/` | evaluation rosters (GSE213621 genetic individuals, GSE281797 participants) and reference per-sample scores |
| `test_release.sh` | reproduces `expected.json` from deposited counts through `score.py` only |
| `MANIFEST.sha256.json` | sha256 of every file above |

The GSE213621 roster fixture lists the 322 representatives kept after collapsing same-person libraries (genotype
identity file), dropping 11 same-person groups with discordant stage and the 9 near-duplicate flags of the sealed run.
The test removes Control from it by the deposited `fibrotic_stage`, giving 257 units. Stage labels are read by the test
from the deposited metadata, not shipped.

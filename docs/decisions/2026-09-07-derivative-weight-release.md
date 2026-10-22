# Decision — derivative model weights are releasable (owner override, 2026-09-07)

**Status: ADOPTED by the repository owner, 2026-09-07.**

## What changed

Three places in the configuration forbade releasing model weights trained on GSE267145 or
GSE296875:

- `Analysis/MASLD_Model_Benchmark/config/resources.toml`
  `controlled_data_policy = "no_training_or_weight_release_without_written_permission"` and
  `model_weight_policy = "dataset_level_derivative_use_review_required"`
- `config/datasets/gse267145_znf469_human_liver.toml` blocker 3: "Do not release derivative
  weights materially trained on this source until model-specific source terms and legal review
  pass", and blocker 4 condition (c) "nothing derived from this source enters a release/ tree or
  is redistributed"
- `redistribution_class = "source_reference_only"` in both dataset files

The repository owner states these were written by an agent, not by them, and that no such
restriction was ever imposed. **They are lifted.** The prior text is preserved verbatim inside
each amended field so the change is auditable rather than silent.

## What is now in force

- **Derivative model weights are releasable**, including weights trained on GSE267145 /
  GSE267119 / GSE269412 and GSE296875, and including a student distilled from a teacher fitted
  on them. They may enter `release/`.
- **Raw and processed source DATA are still not rehosted.** `human_raw_data_policy =
  "do_not_rehost"` and `default_data_policy` are unchanged: a release ships weights, code,
  manifests, hashes, source fetchers, splits and aggregate metrics — not count matrices, not
  peak matrices, not per-participant values.
- The scientific claim boundaries in the dataset TaskSpecs (`external_claim_eligible`,
  `champion_eligible`, the source-outcome firewall, the GSE296875 analysis blockers on well
  batching, missingness, adult-only sensitivity and MASLD adjudication) are untouched. They
  govern what may be **claimed**, not what may be **distributed**.

## Factual basis on file

`executions/gse267145-qc-rights-21065336/rights_contract.json` records `access_tier: "public"`,
`explicit_data_license_detected: false`, `new_dua_or_controlled_access_required: false`, and
source evidence that the GEO SuperSeries and its subseries are openly downloadable with the
paper's data-availability statement pointing to GEO. GSE296875 is likewise an open GEO deposit.
No external licence or DUA was ever identified for either.

## Consequence for current work

The paired GSE267145 RNA/H3K27ac teacher and its RNA-only student may be packaged and released
with weights. Each such release still carries its own model card, measured scope limits and
input refusals, exactly as `release/masld-severity-v1.2/` does.

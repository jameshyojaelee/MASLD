# Blind-phase output and unseal contract

The only outputs authorized before unsealing are source/provenance metadata,
neutral screen coverage fields, missingness indicators, frozen prediction
objects, generic HLF covariates, specifications, hashes, and gate states.

`masked_screen_schema.tsv` has an exact allow-list enforced independently by
the validator. It cannot contain any MAGeCK result value. The raw publisher
workbook stays under `source/raw/` with read-only permissions and is never an
input to the specification freezer.

`detectability_covariates.tsv` is header-only while the source gate is blocked.
It becomes eligible for use only when `01_extract_depmap_hlf.py` writes it from
the complete version-matched DepMap 24Q4 triplet and writes `DEPMAP_READY`.

`SEALED` and `FIREWALL_READY` jointly authorize only the separately hashed
one-pass executor recorded by `PHASE_C_CODE_READY` and
`phase_c_code_manifest.tsv`. `BLOCKED` explicitly forbids outcome access.
Before unsealing, that executor independently reauthenticates the exact sealed
specification, all twelve source gates, official DepMap files, HLF identity,
raw workbook, code/environment bundle, and estimability of every prespecified
design.

Every claim-bearing blind artifact is recorded in `release_manifest.tsv`, is marked
`canonical=FALSE`, and records `external_outcomes_read=FALSE`. Canonical atlas,
figure, and manuscript paths are forbidden.

The pre-unseal bundle owns the following scientific contracts:

| Contract | Required pre-unseal content |
|---|---|
| Source | authenticated workbook, masked exact allow-list, source schema, and public-access status |
| Universe | one row per source gene with mapping, class, neutral guide counts, HLF covariates, eligibility, and exclusion reasons |
| Prediction | four mutually exclusive classes; sealed hepatocyte membership, weights, stage direction, and registry hashes |
| Statistics | complete formulas, estimands, confidence intervals, permutation/matched-null construction, multiplicity, seeds, and accepted draw counts |
| Audit | producer/input hashes, source gates, specification hash, deviations, unseal record, and noncanonical release manifest |

Cross-file DepMap Entrez conflicts must remain explicit in both
`depmap_gene_mapping_audit.tsv` and the screen universe. Conflicted values are
blank, `primary_eligible=FALSE`, and no symbol-only fallback is permitted.

`covariate_coverage_audit.tsv` must separately report finite TPM+Chronos source
coverage and TPM>=1 eligibility attrition for every frozen evidence class. Only
the former is used for the 80% source gate; the latter remains mandatory for
the primary gene and program universes.

The following Plan-40 result files are deliberately forbidden before the frozen
one-pass executor populates `unseal_record.tsv`: `class_effects.tsv`,
`class_protective_hits.tsv`, `program_effects.tsv`, `matched_null_audit.tsv`,
`sensitivity.tsv`, and `fig5_verdict.tsv`. Their eventual schemas must preserve
source row and mapping state at gene level; estimate, effect unit, confidence
interval, attempted/accepted null draws, p/q family, robustness and known-hit
sensitivity at class/program level; and every Boolean component of the
mechanical Figure-5 gate. Creating empty files with these result names is not
allowed because it could be mistaken for a valid null analysis.

`release_manifest.tsv` is the prevalidation manifest and therefore does not
self-hash or include validator-generated `validation_status.tsv`. Source
binaries are authenticated in `source_manifest.tsv`; all claim-bearing blind
tables, status markers, and producer/input manifests are included directly.

After unsealing, the noncanonical phase-C contract additionally requires:

- all 18,343 source rows in `screen_results.tsv`, with source/mapping state,
  frozen eligibility, direct D21 and PA-versus-vehicle positive-selection
  fields, and mechanical protective-hit annotations;
- seven continuous class families, six binary families, six program
  variant/weight families, 100,000 accepted class draws, and 10,000 matched
  program draws where testable;
- explicit permutation and matched-null audits, all mandatory sensitivities,
  and a one-row mechanical `fig5_verdict.tsv`;
- `phase_c_release_manifest.tsv`, `PHASE_C_COMPLETE`, independent
  `phase_c_validation_status.tsv`, and `PHASE_C_VALIDATED`; and
- a valid complete null when no class or program passes. Biological positivity
  is never a release-validation criterion.

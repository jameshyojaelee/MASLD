# Decision: scoped exemption from the ordering-method guard

**Date:** 2026-08-14
**Granted by:** james.hyojae.lee@gmail.com (PI)
**Exemption id:** `ORDERING-EXEMPTION-2026-08-14`
**Workstream:** `CONTINUOUS-AXIS-BENCHMARK-v1`
**Status:** active

## What was decided

Ordering methods (`transcriptomic ordering`, `pseudotime`, `sliding windows`) and their
vocabulary are permitted **inside one named workstream only**, whose object of study is
those methods. They remain forbidden everywhere else.

## Why

Kamzolas et al. (*Nat Metab* 2026) order patients along a data-driven molecular axis,
cut it into overlapping windows, and conclude that MASLD programmes activate in an
ordered sequence that discrete histology staging misses. Our own result from
`RNA-seq/results/stage_axis_geometry/20260813T142216Z/` says the opposite: after
correcting for measurement noise and for the algebraic bias from contrasts sharing donor
groups, stage-associated remodeling is scalar amplification of a single axis.

These are competing, testable models. Testing theirs requires implementing their method.
The guard was written to stop ordering claims from leaking into the Resource, not to stop
the project from evaluating an external method — but as written it forbids both, so an
exemption is needed to run the benchmark at all.

## Why an exemption rather than deleting the guard

The comment at `scripts/manuscript/program_observability_map/systems_contract.R:97-98`
records that these terms "have each leaked back into this project after being retired".
Deleting the guard would remove the only thing standing between a benchmark result and a
Resource claim. The exemption is therefore scoped by path, revocable, and paired with a
claim firewall that actively blocks the leak the guard exists to prevent.

## What changed

| File | Change |
|---|---|
| `config/method_exemptions.json` | New. Single source of truth for exemptions. Missing file or empty array restores unconditional enforcement. |
| `scripts/manuscript/program_observability_map/systems_contract.R` | `assert_language()` gains optional `exemption_id` and `caller_path`. `FORBIDDEN_TERMS` unchanged. New `assert_no_exemption()`. |
| `scripts/manuscript/program_observability_map/run_tests.R` | Seven new cases covering the exemption logic (52 -> 59 checks). |
| `scripts/manuscript/validate_resource_scope.py` | New `check_claim_firewall()`: fails if any of the five human authorities contains an exempt term in prose or cites an exemption scope path. |
| `docs/decisions/2026-08-14-ordering-exemption.md` | This record. |

## What deliberately did NOT change

`scripts/manuscript/fibrosis_nas_map/freeze_contract.py` and its eight frozen
`analysis_contract.json` outputs are **untouched**. Two reasons:

1. Its `forbidden_methods` list at lines 76-79 is read by no code. It is a declarative
   record of what BULK-PROGRAM-MAP-v8 committed to, not an active guard. The v8
   workstream genuinely did forbid those methods, and that record should stand.
2. The script refuses to overwrite a frozen contract that differs from what it would
   generate. Editing the literal would have made `freeze_contract.py` raise on its next
   run against an already-frozen, mode-440 contract.

Rewriting a frozen provenance record to suit a later decision is the failure mode this
repository is built to prevent. The exemption is recorded here and in the config instead.

## Boundaries

- `progression` and `molecular clock` are listed in `never_exemptible_terms` and cannot be
  licensed by any exemption. `progression` is a clinical-course word, not a method word.
- An exemption licenses terms only when the caller path matches its `scope_paths`. An id
  with no path, an out-of-scope path, or an unknown id licenses nothing.
- No output produced under this exemption may appear in Resource claim text, figure
  captions, release tables, or provenance chains. Enforced by `check_claim_firewall()`.
- This exemption does not authorise any Resource-level ordering claim. If the benchmark
  finds that a data-driven ordering adds information beyond recorded histology, that
  result requires its own scope decision before it can enter the paper.

## Version-control status (checked, not assumed)

`.gitignore:8` is a blanket `*` and `.gitignore:114` ignores all of `/docs/`, so:

| File | Status |
|---|---|
| `config/method_exemptions.json` | **gitignored** |
| `docs/decisions/2026-08-14-ordering-exemption.md` | **gitignored** |
| `scripts/manuscript/program_observability_map/*` | untracked, not ignored |

So this exemption cannot be revoked by `git revert`, and any mechanism in this workstream
that depends on "committed to git" does not work. Provenance is carried by content hashes
in read-only mode-440 files instead — the same idiom `freeze_contract.py` already uses.
Every benchmark output directory writes `EXEMPTION.json` recording the exemption id and
the sha256 of `config/method_exemptions.json` at run time, and that hash is included in
`output_checksums.tsv`. A changed exemption config is therefore detectable after the fact
even though it is not versioned.

## Revocation

Set `status` to `"revoked"` in `config/method_exemptions.json`, or delete the entry, or
delete the file. Then re-run:

```bash
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
"${MASLD_PROJECT_ROOT}"/../../../../home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla \
  scripts/manuscript/program_observability_map/run_tests.R
python3 scripts/manuscript/validate_resource_scope.py
```

Guards default to zero exemptions when the config is absent, so a partial revert also
fails closed. Because the config is gitignored, deleting the file is the complete
revocation — there is no tracked copy to restore it.

## Verification performed at grant time

- `run_tests.R` 59/59, `t2_unit_tests.R` 42/42, `t3_run_tests.R` 12/12, `t6_run_tests.R` pass.
- All 14 pre-existing `assert_language()` call sites pass neither new argument, so their
  behaviour is unchanged; confirmed by the suites above passing without edits.
- Claim firewall confirmed to fire on an exempt term in prose and on a cited scope path,
  and confirmed not to fire on a substring inside a longer word.

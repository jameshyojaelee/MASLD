# Program-context v2 release scaffold

The public product is the **MASLD Gene Catalog**. Internal `passport_*` paths,
columns, and artifact IDs below belong to the sealed v1 machine contract and
are retained only so existing candidates remain verifiable.

This directory implements the candidate-only Plan 60 `REL-00` through `REL-05`
release boundary. Real scientific assembly remains blocked until all upstream
workstreams close. The included REL-02--04 implementation is exercised only on
synthetic fixtures in this revision; it never promotes a release or writes a
canonical manuscript or figure path.

The real candidate root is fixed to:

```text
RNA-seq/results/manuscript_release/candidates/
  program-context-v2-candidate-2026-08-07/
```

No script accepts a substitute candidate root. A real snapshot must not be
created until Plans 13, 20, 30, 40, and 50 each have exactly one signed
terminal row and exactly one `terminal_verdict` artifact.

## Components

- `rel00_close_workstreams.py` validates the signed closure and performs no
  writes.
- `rel01_snapshot_candidate.py` copies closed inputs into the exact candidate
  root. It rejects symlink components, verifies source SHA256 and byte size
  before and after copying, uses no-overwrite materialization, records safe
  runtime identity, and checks protected releases before and after the copy.
  The Plan 50 snapshot includes every `relative_path` in the sealed
  `passport_release_manifest.tsv` plus the five-file terminal chain:
  `PASS06_VALIDATED`, the payload manifest, the validation report, terminal
  provenance, and the inner Plan 60 handoff. Its artifact count is therefore
  derived at runtime as the payload-manifest row count plus five. The sealed
  payload includes the self-contained local review UI and its review artifacts;
  these are candidate evidence artifacts, not a hosted deployment.
- `rel02_stage_base_inputs.py` extends the immutable snapshot with a
  coordinator-selected `inputs/BASE/` manifest. It requires exact source paths,
  SHA256, byte sizes, a named coordinator, and a signed-file hash. The fibrosis
  input is the complete 21-file true-Kleiner candidate bundle beneath
  `inputs/BASE/fibrosis_candidate/`, with exactly one primary coefficient-table
  role and 20 validation/provenance roles. A stand-alone CSV is never accepted.
- `fibrosis_candidate_contract.py` revalidates READY identity, accepted-release
  SHA256 pins, the exact producer/validated manifests, all artifact hashes,
  27,638 genes in each of four BH families, the 6/6/6/5-cohort census, 664
  biological donors, and GSE193066 first-biopsy-only handling before preparation
  and again from the frozen BASE snapshot.
- `rel02_build_candidate_tables.py` rederives the cross-sectional fibrosis
  counts, freezes the number/claim/source-dependency/gate ledgers, and writes one
  candidate source table per panel.
- `rel03_render_candidate_panels.py` reads only candidate `figure_sources/` and
  emits deterministic individual PDF panels under the isolated
  `figures/misc/candidates/<candidate>/` branch. Fixture PDFs use normal
  Helvetica at 6 pt; PNG and SVG substitutes are prohibited.
- `rel04_build_candidate_manuscript.py` creates the fixed five-section,
  five-figure Cell Genomics Resource recut below
  `docs/manuscript/candidates/<candidate>/` and seals a deterministic
  post-snapshot product transition. Myojin remains supplementary.
- `rel05_validate_candidate.py` validates candidate-only downstream reads,
  exact snapshot membership, protected-release integrity, producer hashes,
  plotted/manuscript number linkage, scientific check registration, and warning
  adjudications. Its strongest real status remains infrastructure-level; it
  deliberately reports `full_release_pass=false` and never authorizes canonical
  promotion.
- `release_products.py` defines candidate-only path constants and REL-02--04
  schemas. These mirror `FIG_MISC` from `load_figure_data.R` without editing any
  existing figure constant or script.
- `release_common.py` defines the exact TSV schemas and shared guards.
- `tests/` contains synthetic-only fixtures and adversarial tests. Test payloads
  carry no biological result.

## Terminal-state contract

The closure distinguishes these states:

- `accepted_main`
- `accepted_supplement`
- `skipped_source_gate`
- `rejected_qc`

A complete nonconfirmatory analysis is retained as `accepted_main` or
`accepted_supplement` according to its prespecified paper role and records
`indeterminate` evidence calls; it is never relabeled a tested negative.
Coverage limitation is recorded in `gate_verdict`, allowed/prohibited wording,
and the workstream artifacts; it is not a fifth terminal state. A source-gate skip cannot be promoted to main
inference, and a QC rejection cannot enter either figure tier. Plan 20 inclusion
also requires a `READY` gate and exactly one `ready_seal` artifact.

Each terminal artifact manifest has the exact columns:

```text
artifact_id,source_path,source_sha256,source_bytes,snapshot_relpath,
artifact_role,downstream_read_path
```

`snapshot_relpath` and `downstream_read_path` must be identical and must be
beneath the workstream's candidate `inputs/<lane>/` directory. Original source
paths are retained only as checksummed provenance; consumers read the copy.

## Scientific validation boundary

`REL-05` requires one registered report for each of the following future
scientific checks:

```text
hotspot_registry_consistency
spatial_assay_native_consistency
numbers_claim_ledger
figure_panel_manifest
passport_semantics
rebuild_determinism
null_compatibility
manuscript_retired_claims
biological_unit_audit
multiplicity_universe_audit
```

These reports are inputs to this validator; this scaffold does not fabricate or
perform their scientific analyses. A warning is accepted only with a written,
checksummed adjudication. Missing checks, null-only outcomes, and absent panels
must fail or route according to their signed scientific report rather than be
silently converted into a positive release.

REL-02--04 additionally enforce these scientific language boundaries:

- stage panels are cross-sectional true-Kleiner contrasts and are rederived at
  `padj < 0.05` as 798/2,190/3,490/783 DEGs (442/1,510/2,131/504 up and
  356/680/1,359/279 down), with 313/361/307/152 biological donors and
  6/6/6/5 cohorts;
- non-significance maps to `indeterminate`, not `tested_negative`;
- `tested_negative` requires an explicit adequate-negative or equivalence
  criterion in the frozen source row;
- Yakubovsky module 8 keeps descriptive median-slope and inferential
  signed-Stouffer directions as separate rows;
- unresolved Vu/CosMx arrays remain technical/reporting units with donor count
  unknown;
- matched-null Moran dispersion is never labeled sampling SE or CI;
- Myojin is described as assay-specific non-support/nonconfirmatory and remains
  supplementary; and
- Figure 4 is “Physical context, assay observability, and prespecified external
  challenges,” not broad validation.

## Synthetic verification

Run the lightweight tests from any directory with:

```bash
bash scripts/manuscript/program_context_v2/tests/run_tests.sh
```

The tests create isolated temporary project roots. They never touch the fixed
real candidate root or any existing named release.

## Real execution remains gated

Once the coordinating agent has signed a real closure, protected-release
baseline, and base-input selection, the CLIs run serially: closure validation,
snapshot, base-input extension, REL-02 tables, REL-03 panels, REL-04 manuscript,
then REL-05 validation. Do not infer placeholder or “latest” paths from this
README. The coordinator must pass the exact signed files and hashes explicitly.
Promotion remains a separate, explicitly user-approved Plan 60 operation that
is not implemented here.

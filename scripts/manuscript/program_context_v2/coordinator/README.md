# Plan 60 real-release coordinator

This directory prepares and audits candidate-only inputs for the Plan 60
`REL-00`--`REL-05` pipeline. It does not overwrite an upstream workstream,
write a canonical result, authorize promotion, or manufacture a workstream
owner's signature.

The coordinator is upstream-read-only. Each workstream owner publishes one
exact artifact manifest and detached non-promotion signature in that
workstream's candidate root. Only after all five handoffs pass does the
coordinator materialize one no-overwrite preparation bundle under:

```text
scripts/manuscript/program_context_v2/coordinator/prepared/<decision-id>/
```

The fixed release-candidate path is:

```text
RNA-seq/results/manuscript_release/candidates/
  program-context-v2-candidate-2026-08-07/
```

## Locked scientific branch

The locked branch is a five-main-figure Cell Genomics Resource posture with
11 main panels and four supplementary panels (15 panels total):

| Figure | Panels | Locked role |
|---|---:|---|
| Figure 1 | 1A--1B | Dataset interface and evidence observability |
| Figure 2 | 2A--2C | Cross-sectional established-state transcriptomics |
| Figure 3 | 3A--3B | Genetics, phenotype provenance, and eQTL observability |
| Figure 4 | 4A--4B | Assay-native physical context and zonation-adjusted lipid challenge |
| Figure 5 | 5A--5B | Evidence passports and next discriminating experiments |
| Figure S1 | S1A--S1B | Complete nonconfirmatory Myojin HLF stress test |
| Figure S2 | S2A--S2B | Continuous k4/k6 NMF axes and factor-stability audit |

Scientific boundaries enforced by the adapters and validator:

- Figure 2A is cross-sectional stage-associated remodeling, not observed
  longitudinal progression.
- Figure 2B uses the deposited sample-level model
  `asin(sqrt(proportion)) ~ group_binary + dataset + inferred_sex`. Its complete
  universe is 22 cell-type columns: 16 testable in one BH family and six
  structurally unavailable. The audited denominator is 1,221 analyzed of 1,260
  `pass_technical` samples; 39 lack deconvolution and 21 deconvolved QC failures
  are excluded.
- Figure 2C reports two internally robust donor-level stage-associated Hotspot
  programs from the full 117-module registry. These are not external
  validations. The GSE244832 stage/cohort limitation remains explicit.
- Figure S2 preserves k4/k6 as continuous axes only. The frozen supplement has
  1,104 samples, 11,040 loading rows, 10 S2A marks, and eight S2B marks. S2B
  reports factorization/factor-stability metrics and an explicit non-inference
  boundary; it does not report or imply a measured stable patient partition or
  a hard-partition ARI.
- Figure 4 retains all 18 prespecified assay-native dataset-program tests and
  never creates a cross-assay score. Exactly one unique program is robust in
  the current frozen integration. GSE192741 has four independent donors; Vu
  has ten technical sections and unresolved donors. Yakubovsky's descriptive
  and inferential module-8 directions remain separate.
- Myojin is complete but nonconfirmatory and supplementary only.
- Hero genes are post-results illustrations, not a discovery or validation
  universe.
- `tested_negative` is forbidden without a frozen, adequate negative-decision
  rule. Non-significance otherwise remains `indeterminate`.
- Evidence dependence is explicit. Stable accessions and the Plan 50 source
  graph are resolved where available; all other independence statements are
  limited to declared source IDs because cross-alias equivalence has not been
  independently audited.

The immutable 22-program v1 release is never modified by this workflow.

## Workstream-owned handoff contract

Every workstream root must contain:

```text
plan60_terminal_artifacts.tsv
plan60_terminal_artifacts.signature.json
```

The exact path/role/snapshot contract is machine-authoritative in
[`publish_plan60_workstream_handoff.py`](../publish_plan60_workstream_handoff.py)
and is imported by the coordinator validator; there is no second handwritten
artifact list. The current contracts are:

| Workstream | Terminal artifact | Exact artifact rows |
|---|---|---:|
| Plan 13 | `SEMANTIC_V2_READY` | 13 |
| Plan 20 | `SEMANTIC_ADJUDICATION_READY` plus `READY` | 37 |
| Plan 30 | `GEN_TERMINAL_CLOSURE_READY` | 15 |
| Plan 40 | `PHASE_C_VALIDATED` | 21 |
| Plan 50 | `PASS06_VALIDATED` | `N_payload + 5` (derived at runtime) |

Plan 20's 37-row contract includes the six NMF source copies, source manifest,
continuous loading table, figure-source table, validation, producer manifest,
and `NMF_CONTINUOUS_SUPPLEMENT_READY`. The publisher rejects missing rows,
extras, duplicate artifact IDs, duplicate snapshot paths, source drift, unsafe
paths, or an existing handoff. It prevalidates in a private directory and
publishes the manifest/signature pair with no-overwrite hard links; a failed
post-link validation removes only links created by that invocation.

Plan 50 has no static row-count subset. Its publisher contract contains
every `relative_path` in the sealed `passport_release_manifest.tsv` plus
`PASS06_VALIDATED`, `passport_release_manifest.tsv`,
`passport_validation_report.tsv`, `passport_terminal_provenance.tsv`, and
`passport_plan60_handoff.tsv`. The exact count is the payload-manifest row
count plus five. This full bundle includes the self-contained review interface
and its sealed review artifacts; they remain candidate-only and are never
hosted or deployed by this workflow.

The detached signature contains exactly:

```text
attestation_version = plan60_workstream_handoff_v1
candidate_id
workstream_id
manifest_sha256
signed_by
signed_at_utc
canonical_promotion_authorized = false  # JSON Boolean
```

The workstream owner, not the coordinator, runs:

```bash
micromamba run -n spatial python \
  scripts/manuscript/program_context_v2/publish_plan60_workstream_handoff.py \
  --project-root /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design \
  --workstream PLAN13 \
  --signed-by '<owner>' \
  --signed-at-utc YYYY-MM-DDTHH:MM:SSZ
```

Repeat with the correct owner for each workstream. Never copy one owner's
signature to another workstream or use the coordinator as the owner.

## Coordinator preparation

Run the coordinator contract suite in the existing `spatial` environment:

```bash
micromamba run -n spatial bash \
  scripts/manuscript/program_context_v2/coordinator/tests/run_tests.sh
```

After all five real owner handoffs exist, create one unique decision bundle;
the output directory must not already exist:

```bash
micromamba run -n spatial python \
  scripts/manuscript/program_context_v2/coordinator/prepare_real_release.py \
  --project-root /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design \
  --output-root /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/manuscript/program_context_v2/coordinator/prepared/<decision-id> \
  --coordinator '<name>' \
  --signing-date YYYY-MM-DD \
  --signed-at-utc YYYY-MM-DDTHH:MM:SSZ \
  --decision-register-id '<decision-id>'
```

The preparation command validates all scientific and source products in
memory before its first write. It rederives the NMF supplement from Plan
20-owned frozen raw copies and requires exact agreement with the frozen
figure-source and loading tables. The 18 rederived NMF evidence rows are
included in the `hotspot_programs` presentation adapter so they participate in
the same source-dependency and provenance audits.

A successful preparation contains:

- signed workstream closure, protected scopes, and protected baseline;
- the complete accepted true-Kleiner fibrosis bundle: 21 selected files under
  one BASE prefix, exactly one primary result role, 20 validation/provenance
  roles, and fixed READY/validated-manifest/primary-result SHA256 pins;
- locked 11-main/4-supplement blueprint;
- seven presentation adapters and `adapter_provenance.tsv`;
- `source_alias_audit.tsv` with the declared-ID-only limitation where needed;
- recursive producer, runtime, package, `sys.path`, and import-resolution
  manifests;
- base-input selection and detached coordinator attestations; and
- a checksummed preparation manifest plus `COORDINATOR_PREPARATION_READY`.

The accepted fibrosis pins are READY
`f6a6d274124749b5c8bbdbf204629c847a33b45bc9d2a570ba98352d895fbcc9`,
validated manifest
`9e8dc350e90bdeac45633007f079bcb5cf5e08f21789078adc50212bcbb64994`,
and primary coefficient table
`01a942d33fb321631db0a3cb0579b2b89b96373c399d61f87030a49ed9788f18`.
A different but internally self-consistent bundle at the same path fails before
any preparation write and fails again if substituted in frozen BASE.

Every adapter provenance row binds the adapter output hash to exact input
hashes and to the required current producer hashes in the frozen recursive
producer manifest. All adapters bind `coordinator_contract.py` and
`prepare_real_release.py`; the Hotspot/NMF adapter additionally binds
`519_freeze_nmf_continuous_supplement.py`. The seal must state:

```text
rel00_05_executed = false
candidate_root_created = false
canonical_promotion_authorized = false
```

## Serial candidate build and retained clean rebuilds

The following is an operator checklist, not promotion authorization. Every
path and hash must come from one accepted preparation bundle; never substitute
`latest`, infer a path, or read a mutable live result after snapshotting.

1. Run `rel00_close_workstreams.py` read-only.
2. Run `rel01_snapshot_candidate.py` to create the immutable authority
   candidate snapshot.
3. Run `rel02_stage_base_inputs.py` with the exact prepared selection and hash.
4. Run `rel02_build_candidate_tables.py`.
5. Run `rel03_render_candidate_panels.py`.
6. Run `rel04_build_candidate_manuscript.py`.
7. Materialize two retained source projects with
   `materialize_clean_rebuild_project.py`, then independently repeat
   REL-01--REL-04 in each root using its rehydrated, hash-checked closure and
   BASE selection. The materializer never copies the authority candidate.
8. Compare the authority candidate with both retained builds.
9. Build the scientific validation registry.
10. Run REL-05, then rerun the registry with `--check-only` and finish with
    `rel05_validate_candidate.py --no-write-report`.

Retained-build topology is fixed. For comparison ID `<comparison-id>`, the
filesystem slots are literally `build_a` and `build_b`; comparator build IDs
remain explicit report labels:

```text
RNA-seq/results/manuscript_release/candidates/
  program-context-v2-candidate-2026-08-07/
    retained_clean_rebuilds/<comparison-id>/
      build_a/project/
        RNA-seq/results/manuscript_release/candidates/
          program-context-v2-candidate-2026-08-07/
      build_b/project/
        RNA-seq/results/manuscript_release/candidates/
          program-context-v2-candidate-2026-08-07/
```

The namespace, comparison, build, `project`, and candidate components must be
real directories, not symlinks. Both isolated builds must be retained after
comparison because REL-05 revalidates their live product hashes.

After the authority REL-01--04 candidate is complete, materialize each source
project (each target must be absent):

```bash
micromamba run -n spatial python \
  scripts/manuscript/program_context_v2/coordinator/materialize_clean_rebuild_project.py \
  --authority-project-root /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design \
  --preparation-root /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/manuscript/program_context_v2/coordinator/prepared/<decision-id> \
  --comparison-id '<comparison-id>' \
  --build-slot build_a

# Repeat with --build-slot build_b.
```

The helper creates a detached sparse Git worktree at the frozen authority
commit and deliberately does not replay the authority workspace's global
tracked diff. It restores exact producer bytes from the frozen producer
manifest, rehydrates every closure/BASE source from immutable snapshots, and
reconstructs protected files against the frozen baseline. The deterministic
`rebuild_source_spec` binds the commit, exact closure/artifact universe, BASE
selection, protected-source contracts, and frozen current producer hashes; it
excludes the authority tracked-diff hash, runtime inventory, and execution
context. The isolated worktree may contain only those exact source overlays
and may never contain a dirty gitlink. The helper writes
`MATERIALIZATION_READY.json` outside the isolated project only after repository,
closure, BASE-source, producer, and protected-baseline checks pass. A failed
build remains unsealed and inspectable; the helper never deletes it.

The authority working tree's tracked binary diff SHA256 (audited 2026-08-08)
is retained as provenance:
`125f35ce19d48d2258fe9e8a73e9d3212cb6d0a7292095dc12054e1e7d125078`
and it includes the unrelated dirty gitlink
`GWAS/finemapping/src/chrombpnet_variant_scorer` at gitlink
`0e1e34199e63112aa618748bb79a206fc491300a-dirty`. That unrelated state is not
copied into or used to identify the scientific A/B rebuilds. A dirty gitlink
inside either isolated rebuild still fails closed.

After authority and retained A/B REL-01--REL-04 builds are complete, run:

```bash
micromamba run -n spatial python \
  scripts/manuscript/program_context_v2/coordinator/compare_clean_rebuilds.py \
  --authority-project-root /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design \
  --build-a-id '<build-a-id>' \
  --build-a-process-id '<derived-execution-identity-a>' \
  --build-b-id '<build-b-id>' \
  --build-b-process-id '<derived-execution-identity-b>' \
  --comparison-id '<comparison-id>'
```

The report is written only in the authority candidate at:

```text
manifests/two_clean_rebuild_comparison.tsv
```

The process-ID arguments are checks, not operator labels. Each must exactly
equal the identity independently derived from that build's immutable REL-01
`execution_context.json`:
`slurm:<SLURM_JOB_ID>:array:<SLURM_ARRAY_TASK_ID-or-none>:step:<SLURM_STEP_ID-or-none>`
when present, otherwise `local:<host>:<actual-pid>:<UTC-recorded-time>`.
Authority, A, and B must have three distinct context hashes and identities.
Each REL-02--04 transition must
bind its product manifest to that candidate's context-bound candidate state;
a copied candidate/product tree plus invented CLI labels is rejected.

The comparator requires identical `rebuild_source_spec` and closure identity,
not identical global snapshot provenance. Authority/A/B snapshot-spec hashes,
execution contexts, repository tracked-diff provenance, and context-bound
transition files may differ and must each validate internally. Scientific
tables, panel-source tables, PDFs, manuscript products, and their exact product
phase, producer, role, path, SHA256, and byte size must agree across all three
builds. The context and transition paths and hashes are frozen into every
report row. `revalidate_clean_rebuild_report()`
recomputes all three live roots; copying a prior report without retaining its
builds cannot pass.

## Scientific registry and REL-05

After the three-way comparison passes:

```bash
micromamba run -n spatial python \
  scripts/manuscript/program_context_v2/coordinator/build_scientific_validation_registry.py \
  --project-root /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
```

The registry rederives the ten required scientific reports rather than
accepting hand-authored pass rows. It also freezes the candidate snapshot
specification hash and exact producer hashes for the registry implementation.
`--check-only` rederives those predicates and verifies the existing frozen
reports without writing.

Key fail-closed checks include:

- all 117 Hotspot modules, the full BH family, the exact Figure 2 projection,
  and the GSE244832 estimability boundary;
- all 22 composition columns and an independent rederivation of the 16-test BH
  family;
- the 1,104-sample/11,040-loading NMF supplement and all 18 S2 marks;
- all 18 assay-native spatial rows, one unique robust program, independent
  donor/technical-unit annotations, and separate Yakubovsky directions;
- the BBJ EAS-GWAS/EUR-eQTL limitation and the nonenriched genetic/state
  interface result;
- Plan 50 graph-resolved passport semantics, exact Figure 5 graph references,
  and post-results hero-gene selection;
- exact figure-source-to-PDF/manuscript number linkage;
- three-way clean-rebuild identity and producer bindings; and
- null/nonconfirmatory outcomes that leave the locked figure branch intact.

The final read-only check is:

```bash
micromamba run -n spatial python \
  scripts/manuscript/program_context_v2/coordinator/build_scientific_validation_registry.py \
  --project-root /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design \
  --check-only
```

For exact REL-00--REL-05 flags, use each script's `--help` and the parent
[`README.md`](../README.md). Do not use `--fixture-mode` for a real candidate.

## Promotion boundary

This directory contains no promotion command. A fully validated candidate is
still not canonical. Promotion remains a separate Plan 60 decision requiring
explicit user authorization, protected-release revalidation, and a dedicated
implementation that is intentionally absent here.

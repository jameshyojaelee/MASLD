# Spatial-context v2 contract tooling

This directory implements Plan 13 tasks `SP-INT-01` through `SP-INT-07`.
Yakubovsky outcomes are imported only from the passed, frozen Plan 11 adapter;
they are never refit here.

- `01_v1_preservation.py` freezes and later compares the byte-level v1 artifact
  census. Candidate-v2 tooling and generated Python bytecode are outside the
  v1 scope; the scope itself is written to `v1_preservation_scope.tsv`.
- `02_build_synthetic_contract.py` writes the required adapter schemas and
  synthetic-only pass/skip fixtures using the sealed Hotspot registry.
- `03_validate_contract.py` independently rederives Hotspot `READY` linkage,
  validates adapter artifacts and biological-unit counts, and proves that
  malformed registry hashes, counts, evidence states, uncertainty intervals,
  and universal-score columns fail closed. Interval bounds are generic and
  must carry an explicit type; `donor_effect_range` is descriptive donor
  heterogeneity and never a confidence interval.
- `04_freeze_visium_rerun.py` freezes the v1/v2 registry comparison, tested
  families, source paths, exact engine parameters, and the fail-closed reuse
  verdict. The frozen v2 family is the two Plan 20 `robust_display` programs;
  the 22-program v1 matched-control pool and BH family therefore cannot be
  reused.
- `visium_rerun_lib.py` translates the sealed v2 registry to the exact input
  schema expected by the byte-pinned v1 engine. Candidate code calls the
  original engine functions directly and never calls its canonical-writing
  `main()` function.
- `05_run_visium_rerun.py` writes atomically to `native_spatial/v1_regression`
  or `native_spatial/v2_candidate`. A v2 run is impossible until the
  independent validator has sealed a passing v1 regression.
- `06_validate_visium_rerun.py` independently rederives BH, robust flags,
  tested families, biological-unit counts, null and matched-set linkage, and
  v1 semantic equivalence. It writes `READY` only after all checks pass.
- `07_freeze_map_selection.py` freezes both members of the complete sealed
  Figure 2 external-test family before Yakubovsky program outcomes exist. It
  reads no external spatial outcome.
- `08_build_real_adapters.py` closes GSE287826 as `skipped_no_donor_key`,
  preserves the repository's explicit Govaere GeoMx exclusion, and imports the
  accepted Govaere CosMx IL32 proximity result at its native physical-array
  unit. It performs no whole-program refit.
- `09_validate_real_adapters.py` validates those three adapters against the
  ten-artifact contract, independently rederives the four CosMx array rows and
  panel coverage, proves no program score/effect was invented, and emits only
  a partial-seam `REAL_ADAPTERS_READY` marker. It explicitly records that Plan
  13 is incomplete.
- `10_run_atac_v2_projection.R` reuses the validated donor-pseudobulk promoter
  accessibility logic in an isolated candidate root. Before testing the
  frozen two-program v2 family, it must reproduce the corresponding accepted
  v1 effects, exact permutation P values, sensitivities, coverage, and donor
  scores to numerical tolerance.
- `11_build_protein_atac_adapters.py` emits the common adapter contract for
  PXD051911 and the two snATAC cohorts. The two v2 programs remain terminal
  `untestable` in DIA-MS under the frozen 8-gene/20%-L1 gate; the protected
  25-protein table is copied byte-for-byte as a selection-conditioned,
  descriptive source-native display. Static GWAS/open-promoter links are also
  descriptive and never receive inferential program statistics.
- `12_validate_protein_atac_adapters.py` independently re-enumerates donor
  label permutations, rederives Welch standard errors and complete two-program
  BH families, rechecks protein observability and static link counts, and
  verifies the 591-file v1 census.
- `final_integration_lib.py`, `13_build_final_integration.py`, and
  `14_validate_final_integration.py` implement the final assay-native long
  contract, Figure 4 source matrix, dataset-level verdict, independent
  rederivation, and candidate READY seal. They prohibit cross-assay scoring and
  retain both frozen hepatocyte programs in robust, null, untestable, and
  skipped branches.
- `15_freeze_final_integration.py` freezes the final producer, validator,
  semantic tests, and wrappers before real native-spatial outcomes are
  available. The real wrapper hard-requires passing v1-regression and v2-native
  READY seals and explicitly imports Plan 11's validation-only hotfix manifest.
- `tests/test_visium_rerun.py` hand-checks graph islands, Moran's I,
  donor-first section collapse, complete-union matched-control exclusion,
  deterministic registry translation, and score sensitivities.
- `tests/test_real_adapters.py` checks source-native row isolation, unique-gene
  coverage accounting, and fail-closed MASH-array direction semantics.
- `tests/test_protein_atac_adapters.py` hand-checks exact permutation and BH
  arithmetic, complete-family uniqueness, and byte-identity failure behavior.

The scripts use only the Python standard library and are lightweight enough for
login-node contract testing. Run the wrapper once on a new candidate root:

```bash
bash Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_sp_int_01_02.sh
```

`01_v1_preservation.py baseline` intentionally refuses to overwrite an existing
baseline. Use the individual build/validate/check commands for later rechecks.

Freeze `SP-INT-03` once, before any heavy compute:

```bash
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python \
  Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/04_freeze_visium_rerun.py write
```

The production sequence is intentionally serial at the scientific gate. Both
jobs use the `cpu` partition, single-word job name `spatial`, 16 CPUs, 128 GB,
and 72 hours. They have been created but are not submitted by this workstream.

```bash
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_visium_v1_regression.sbatch
# Submit only after the first job writes native_spatial/v1_regression/READY:
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_visium_v2.sbatch
```

The v2 candidate remains isolated. No file under
`Analysis/Multimodal_Program_Projection/results/`, `Analysis/Spatial/results/`,
or `figures/main/` is a permitted output.

The bounded real-adapter seam is run separately because it does not complete
cross-assay integration or render Figure 4:

```bash
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_real_adapters.sbatch
```

Its real registry is intentionally `real_adapter_registry.tsv`; the original
`adapter_registry.tsv` remains the synthetic contract fixture until all real
Plan 13 lanes are ready for final assembly.

The bounded protein/ATAC lane runs independently of the compatible-Visium jobs
and never invokes a canonical v1 producer:

```bash
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_protein_atac_v2.sbatch
```

Its `protein_atac/READY` marker is candidate-only and deliberately records
`plan13_complete=FALSE`; final integration and promotion remain Plan 60 work.

The outcome-free final-integration fixture is run separately:

```bash
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_final_integration_fixture.sbatch
```

Freeze the implementation only after that fixture and its independent
validator pass:

```bash
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python \
  Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/15_freeze_final_integration.py write
```

Run the real assembly only after both native spatial seals exist and pass:

```bash
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_final_integration.sbatch
```

This final Plan 13 lane writes source tables and a panel verdict only. It
renders no PDF and never promotes into canonical results; those actions belong
to Plan 60.

## Spatial Resource candidate layer

Scripts `20`–`28` extend the sealed semantic-v2 work into the manuscript-facing
spatial transportability and observability layer. They write only to
`Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11`.
They never mutate the frozen 117-program registry, the semantic-v2 source
bundle, canonical portal data, or `figures/main/`.

- `spatial_resource_datasets.tsv` is the manifest-driven source and claim
  contract. It records biological and technical units, join resolution,
  source dependence, permitted estimands, prohibited claims, license, and
  processing provenance. Donor, clinical, histology, and metabolite requests
  fail closed when their corresponding unit or join is unresolved.
- `20_build_resource_registry.py` seals that registry and checksums every
  registered assay gene-axis source.
- `21_build_resource_coverage.py` builds the complete outcome-free
  117-program-by-assay coverage product. `observable` means at least eight
  genes detected in at least 1% of assay units and at least 20% of the frozen
  positive L1 weight; `partial` and `untestable` are retained explicitly. Gene
  axis presence and adequate detection are separate fields. This table contains
  no P values, validation labels, or winner selection.
- `22_build_hmsma_label_blind.py` processes all 35 HMSMA arrays and 130,097
  retained spots as technical coverage. It computes within-array descriptive
  organization for the two predeclared programs on exact Visium hex-lattice
  neighbors, but emits no population P value and no biological sample count.
- `23_build_resource_effects.py` reuses the sealed 9,999-draw matched-null
  results for the two predeclared programs, qualifies Vu as source-dependent,
  and appends HMSMA only as metadata-pending descriptive context.
- `24_assemble_resource_candidate.py` produces the four release-linked portal
  products plus a categorical display contract. Legacy `s5_spatial` is
  compatibility-only and explicitly excluded from rankings, active-layer
  counts, strength bars, and default decisions.
- `25_validate_resource_candidate.py` independently checks the 117-by-assay
  and two-by-source Cartesian families, registry and membership hashes,
  donor/unit gates, the asymmetric IGFBP7-supported/BICC1-indeterminate result,
  portal semantics, and v1 numerical preservation.
- `26_render_resource_figure4.py` renders candidate-only 4A, 4E, and 4F panels
  from sealed tables. Accepted panels 4B–4D are regenerated under a candidate
  figure-output root. `27_validate_resource_figure4.py` requires six one-page
  PDFs and six exact source tables before writing a candidate READY marker.
- `28_validate_hmsma_metadata_join.py` is the dormant metadata-branch gate. It
  enumerates the 35 registered array IDs directly from the H5AD and rejects
  duplicated/unmatched arrays, provisional donor identities, ambiguous repeated
  donor relationships, incomplete clinical covariates, and unregistered H&E or
  MALDI coordinates. It performs no outcome analysis and writes no release.
- `29_render_resource_downstream_figures.py` consumes the sealed registry,
  coverage, effects, and HMSMA label-blind tables to render four candidate-only
  downstream panels: the Figure 1 assay-capability matrix, the Figure 5 spatial
  Gene Catalog examples, the all-117-program coverage heatmap, and the HMSMA
  descriptive organization distribution. It writes a new immutable integration
  candidate and never writes under canonical `figures/main/`.
- `30_validate_resource_downstream_figures.py` independently rederives the
  plotted values and complete families, rejects donor-count leakage and HMSMA
  inferential P values, checks one-page PDFs and non-Type-3 fonts, and seals exact
  panel/source-table manifests with canonical promotion disabled.

Run the real candidate build on a compute node:

```bash
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_spatial_resource_candidate.sbatch
```

The resulting READY markers explicitly state that canonical promotion is not
authorized. Corrected-COLOC-dependent gene-class tables, gsMap refreshes, and
final Gene Catalog entries remain blocked on their upstream synchronized release and are
not synthesized by this wrapper.

Render the release-linked Figure 1/Figure 5 candidate insets and supplement
panels only after the spatial candidate and Figure 4 are sealed:

```bash
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_spatial_downstream_figures.sbatch
```

This writes only to
`Analysis/Multimodal_Program_Projection/candidates/spatial-figure-integration-candidate-2026-08-11-r2`.
Final assembly into Figures 1 and 5 remains synchronized-release work requiring
explicit adjudication; genetics-dependent Gene Catalog content is not synthesized.

Candidate portal consumers are implemented in the source-capability audit,
gene spatial tab, and program spatial-observability tab. They read the four
release-linked products directly and fall closed when those files have not yet
been staged. The public data directory is intentionally untouched. Figure 1's
capability matrix and Figure 5's spatial Gene Catalog fields must consume the same
sealed registry/effect tables during synchronized assembly; their final panels
are not rendered here because their non-spatial inputs remain release-blocked.

## Additional spatial falsification and diagnostic lanes

Scripts `31`–`38` are bounded, candidate-only analyses that test or explain the
two-program result without expanding the confirmatory family. They do not
change Figure 4 calls, portal rankings, the 117-program coverage product, or
canonical outputs.

- `31_run_full_composition_sensitivity.py` replaces the primary
  hepatocyte-abundance adjustment with all 16 shared cell2location q05 factors
  plus the same two QC covariates. It reuses the frozen program weights,
  graphs, matched-control hashes, seeds, null counts, and donor-collapse rules.
  `32_validate_full_composition_sensitivity.py` independently rederives the
  nulls and BH families and permits bounded floating tolerance only for two
  recomputed matching-correlation summaries; keys, bins, counts, graphs, and
  matched-set hashes must remain unchanged.
- `33_render_unit_null_diagnostics.py` renders all physical-unit estimates,
  matched-gene null calibration, and frozen/equal/leave-top weight
  sensitivities. `34_validate_unit_null_diagnostics.py` enforces the five
  GSE192741 sections to four-donor collapse, keeps all ten Vu arrays
  source-dependent, and forbids section- or array-level P values.
- `35_build_cell_context_attribution.py` reports the complete two-program by
  16-factor descriptive association family within each section or array, then
  collapses GSE192741 donor-first. `36_validate_cell_context_attribution.py`
  rederives all 480 correlations from stored spot values. Both the source
  tables and panel state that program scores and abundance estimates share the
  same Visium RNA matrix, so this is context covariation rather than lineage
  origin or independent validation.
- `37_build_hmsma_mask_sensitivity.py` evaluates the sealed 500-UMI proxy mask
  against stricter 1,000- and 2,000-UMI masks for all 35 arrays and both
  programs. It cannot evaluate a looser mask because the available H5ADs were
  already filtered at 500 UMIs. `38_validate_hmsma_mask_sensitivity.py`
  rederives the 210-row family, reproduces the sealed baseline exactly, hashes
  all 35 input H5ADs, and prohibits population inference.
- `39`–`44` preserve the full-composition figure's visual-review history. The
  provenance-complete paper-eligible panel is rendered by
  `43_render_full_composition_figure_r4.py` and checked by
  `44_validate_full_composition_figure_r4.py`; earlier r2/r3 panels remain
  immutable visual provenance.
- `45_build_spatial_diagnostic_bundle.py` and
  `46_validate_spatial_diagnostic_bundle.py` index the five accepted diagnostic
  components, six PDFs, all candidate artifacts, and their executable code
  dependencies under one candidate release. The bundle does not copy or
  promote canonical outputs.
- `47_render_spatial_publication_figures.py` rebuilds the spatial figure family
  around one decision per panel. It simplifies 11 panels, retains the already
  compact Figure 4E maps and 117-program coverage heatmap byte-for-byte, and
  moves long claim boundaries into exact source tables. The focused
  cell-context panel shows only the prespecified IGFBP7–fibroblast and
  BICC1–cholangiocyte comparisons; the complete 64-row factor table remains
  sealed upstream. `48_validate_spatial_publication_figures.py` checks all 13
  one-page PDFs, exact source rederivation, compact canvas bounds, non-Type-3
  fonts, retained-panel identity, canonical-write prohibition, and the measured
  reduction in plot text. The visually accepted immutable output is
  `spatial-publication-figures-candidate-2026-08-11-r4`.
- `49_build_spatial_impact_figures.py` turns that compact family into a
  result-led assembly without changing inference. Figure 4A now exposes the
  calibration logic, Figure 4F plots observed residual Moran statistics against
  matched-gene null quantiles with biological/physical-unit sign counts, and
  the former complete source-state matrix moves to the supplement.
  `50_validate_spatial_impact_figures.py` independently rederives the exact two
  frozen program hashes, all four 9,999-draw null results, 4/4, 10/10, 4/4, and
  9/10 positive-unit counts, composition sensitivity, panel identity, and claim
  boundaries. The visually reviewed immutable output is
  `spatial-impact-figures-candidate-2026-08-11`; it contains 14 one-page PDFs,
  adds no statistical family, and never writes canonical figures.

Run the lanes independently on `cpu`:

```bash
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_full_composition_sensitivity.sbatch
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_unit_null_diagnostics.sbatch
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_cell_context_attribution.sbatch
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_hmsma_mask_sensitivity.sbatch
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_full_composition_figure_r4.sbatch
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_spatial_diagnostic_bundle.sbatch
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_spatial_publication_figures.sbatch
sbatch Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/run_spatial_impact_figures.sbatch
```

Only the latest candidate carrying `READY` and a completed visual review is
paper-eligible. Failed validators and visually superseded candidates remain
provenance, not results.

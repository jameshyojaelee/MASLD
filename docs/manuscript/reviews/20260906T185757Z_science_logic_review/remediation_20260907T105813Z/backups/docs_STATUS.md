# Current paper and codebase status

**Verified:** 2026-08-17 EDT  
**Framing updated:** 2026-08-14 EDT  
**Continuum integration verified:** 2026-08-21 EDT  
**Scope:** standalone MASLD Resource paper  
**Important:** live SLURM state must be queried directly; the snapshot below is dated.

## Executive assessment

The binding title and framing are **A multimodal human MASLD Resource maps
inherited genetic variation and multicellular disease remodeling**. The two
primary maps are connected to donor-aware cell programs and assay-native
protein, chromatin, and spatial tissue context, then delivered through the
MASLD Gene Catalog. Observability, provenance, and evidence states remain the
methodological guardrails, not the headline. The immediate task is
**reconciliation, not expansion**: propagate the corrected COLOC release adopted
on 2026-08-17 through its ~30 downstream consumers, propagate the chosen
five-cohort fragment-count bulk substrate through dependent analyses, regenerate
figures that have changed, and settle one synchronized release. The paper does not need a
new public-data mechanism or Cas13 result to be complete. The recurring
biological examples are now fixed by question. A provisional direct-NAFLD GNMT
locus bridges inherited variation to the severity-linked
GNMT–MAT1A–CYP2C19 RNA–stage–protein state. The Stromal ECM versus
Ductular injury program contrast separates disease association, tissue organization, and
lineage attribution. Current official COLOC does not support a primary
SuSiE-colocalized mechanism. The aligned ambient-RNA and donor-paired
cross-lineage candidates are independently validated but not adopted. A
Kamzolas-adapted molecular continuum complements recorded fibrosis stage with
stage-adjusted donor-level resolution. Its additive Figure 4F–4H and complete
S3/S4 panel package were adopted into the working figure authorities on
2026-08-24; the synchronized Resource release has not been adopted.

## What is official, updated, or open

**Continuum translation extension, 2026-09-06: implemented candidate.** The
complete Gene Catalog sidecar, integrated 117-program figure, signed genetic
records, composition sensitivity and matched liver/plasma protein analysis are
available. Main inherited memberships remain 462 SuSiE and 822 ABF genes.
Only 83 of 925 supported main-trait signal pairs resolve direction under the
declared conservative marginal-effect rules; multiple signals and ambiguous
alleles remain unresolved. Of 117 programs, 93 meet the existing stage/sex
two-axis support rule, 68 meet that rule with composition adjustment, and 26
also pass BH separately in each evaluation cohort and axis after composition
adjustment. These counts are different criteria, not independent replications.
The public protein data link 41 participants; 49 of 345 shared proteins pass
the adjusted matched-association BH test. Genetic cross-study direction
conflicts remain explicit. Details and completed checks are in the
[technical review](technical/continuum_translation_20260906T135538Z/IMPLEMENTATION.md).
The candidate is not an adopted figure or a treatment-response claim.

| Area | Current state | What may be used now | Next dependency |
|---|---|---|---|
| Paper scope | **Settled** | Resource/framework framing and six-figure story | Keep all drafts inside [PAPER.md](PAPER.md) |
| Bulk disease-vs-control | **Validated update selected for the paper; main Figure 3C alluvial and Figure 3D endpoint candidate complete; synchronized release pending** | Canonical inference remains 844 participants (157 controls, 687 disease) across five cohorts. Main Figure 3C shows per-cohort DEG membership; the four-cohort 418-participant endpoint analysis is sensitivity evidence only, recovering 1,307/1,347 canonical DEGs with all-gene log2FC rho 0.938 | Rebuild remaining dependent fields and adopt them only as one synchronized release |
| Noncoding RNA source check | **Passed under a documented amended cohort-level rule** | All five cohorts retain every library; the stricter all-sample rule remains preserved as a failed sensitivity | Use only the passed, hashed source-check candidate for lncRNA outcome work |
| lncRNA identity and program content | **Validated candidates; not adopted** | Versioned GENCODE v49 identities, seven genomic classes, mapping exclusions, unchanged 117-program content, exact pooled biotype subdivisions, 271 high-confidence disease-associated lncRNAs, and 48/49 passing native leave-all-lncRNA-out sensitivities | Integrate these candidates with composition and the synchronized release |
| Stage analysis | **Five-cohort synchronized candidate validated** | Clean job `19966029` produced checksum-verified adjacent fibrosis, adjacent NAS, F-versus-F0, and per-cohort models with complete 23,370-gene BH families from the validated 844-participant DGE. Candidate v1 remains failed execution provenance after job `19957749` was edited during execution | Render and visually adjudicate candidate v2; hold off on adoption |
| Histology-anchored molecular continuum | **Working Figure 4F–4H and S3/S4 panels adopted; synchronized release pending** | The two co-primary scores are modeled inferentially in GSE162694 and GSE213621. Descriptive windows show all five Resource cohorts, while stage panels show the four cohorts with recorded fibrosis. The signature-excluded atlas retains 5,329 pathways and 116 Hotspot programs; 5,439 mapped nodes form 43 size-ranked systems and six nodes remain explicit microcomponents. The separate donor-level 43-system analysis adjusts for fibrosis stage and sex within each evaluation cohort: 33 systems pass both co-primary fixed-effect meta-analysis BH/direction criteria, 29 replicate with feature-balanced scoring, 22 pass BH-adjusted within-stage/sex permutations in all four tests, and 18 pass maxT family-wise P<0.05 in all four tests. Figure 4H shows all 43 fibrosis-stage versus stage-adjusted continuum effects, labels the 15 largest systems fixed before outcome overlay, and outlines the 18 strict maxT-supported systems; rho=0.658 and 36/43 same-direction summaries remain descriptive. Active S4U preserves all 43 rows, equal-cohort window values, and the display-only pooled late-minus-early row order; the contrast tile is omitted and one neutral `*` marks each of the 18 strict maxT-supported systems, while the 33 primary BH calls remain table-only. Figure 4G/S4Q–S4U contain no embedded explanatory footnotes; those explanations remain in the external caption authority. Window tiles remain descriptive; source-overlap cohorts remain outside primary inference, and the systems are not a biological ontology. | Retain Figure 3; carry Figure 4G–4H and S4Q–S4U into the synchronized release audit with their source tables and claim boundaries |
| Noncoding RNA figure package | **Nine-panel candidate independently validated and visually reviewed** | Revised Figure 3 pooled/stage panels, six lncRNA/program supplementary panels, and the Figure 6 molecular-object routing panel under `noncoding-resource-2026-08-12/plot-upgrade-v3`; final repository validation job `19779385` passed | Assemble only after synchronized adoption |
| Noncoding DNA figure package | **Figure 2D synchronized and adopted with main Figure 2 A-I; remaining package independently validated** | The title-free upright locus-balanced consequence panel is canonical in the 2026-08-17 round-2 Figure 2 release. Promoter/ABC/open-chromatin context, six-lineage accessibility, the 30-variant heatmap, and the full noncoding-mass distribution remain validated candidates under `noncoding-resource-2026-08-12/noncoding-dna-v5` | Build target-transcript biotype from the adopted COLOC release and retain the broader DNA package as supplementary candidate material until synchronized selection |
| Cell composition | **Synchronized candidate rebuilt; render validation active** | Donor-aware hepatocyte and non-parenchymal summaries use the accepted composition source and the same stage sample manifests as Figure S3 | Validate candidate v3 and keep canonical output unchanged |
| Hotspot programs | **Aligned ambient and cross-lineage candidates independently validated; working main art selected, not synchronized** | All 117 program memberships and weights remain fixed. The common-universe ambient model corrects five datasets and exposes GSE189600 and Liver Atlas as failed decontX datasets; 106 programs pass integer-count transport, while 11 T-cell programs are untestable. Same-atlas donor-paired localization meets the predefined family criteria for the Stromal ECM and Ductular injury programs. The shared-axis 4D geometry passed adoption validation in job 20878079; the gradual 4E bulk line plot now uses uniform points and no program-significance size encoding, validated in job 20878757 | Adjudicate with the synchronized manuscript release; do not describe localization as transcript origin or cell autonomy |
| Continual single-cell integration | **Ran to a terminal decision; primary adoption REJECTED on 2026-08-15; permitted only as supplementary sensitivity** | The corrected scvi-tools 1.3.3 replay-plus-EWC port ran end to end on real GPU fits (V100/L40S, seed 17, one all-lineage SCANVI plus five lineage models) against the read-only 275/104 descriptive and 273/102 analyzed census, the Harmony export built to the current requirements, and the fixed 1,232,318-cell x 4,000-HVG count object; the final regression suite passed 136 tests with 0 failures and 1 skip. `candidates/final-adoption-decision-v42-2026-08-15T215000Z/promotion_decision.json` is terminal: 13/14 criteria met pass and `independent_external_routing` FAILS on the untouched GSE212837 cohort (12 donors, 128,851 cells; donor-balanced macro-F1 0.2756 against the predefined 0.70 floor). The method identity is a **replay-plus-EWC-guided geometry-preserving adapter**, not the unmodified preprint method. `candidates/router-v2-final-decision-v46-2026-08-17T233700Z/decision.json` then rejects the pretrained-CellTypist Router V2: on the untouched GSE296875 cohort (39 donors, 68,398 cells) overall donor-balanced macro-F1 0.7071 clears the 0.70 floor and paired improvement over V40 is 0.2033 (95% CI 0.1847-0.2219), but B cells 0.2888 and NK-T 0.4229 fail the predefined 0.50 per-class floor. Harmony remains the primary representation, and the strict seven-donor reference is retained after both the local nine-donor arm and the common-universe 26-donor HLiCA Andrews_2022/Andrews_2024 expansion failed the reference-preservation criteria; wholesale HLiCA ingestion remains forbidden. The canonical H5AD, manuscript, paper authorities, and the byte-identical 117-program separation are unchanged. ⚠ the documented full-scale sbatch path never executed. B6K jobs `20208726` and `20215473` both sat pending from 2026-08-13 and were cancelled at 00:00:00 elapsed, so every fit ran instead inside `gpu_interactive` V100/L40S allocations and the persistent CPU allocation `19664223`. ⚠ `secondary_technical_stress` passed on technical stress alone with `independent_confirmation_pass: false` because GSE136103 is not independent and GSE174748 is underpowered. ⚠ reference fits use the prespecified 2,000-cell per donor-by-label cap, so 61,960 of 216,957 strict-reference cells entered training while embeddings still cover all 216,957 | Nothing on the critical path. Cite the candidate only as supplementary sensitivity and never as the Figure 4 representation. GSE212837 and GSE296875 are diagnostic-only and cannot tune any router or threshold. Any continuation is a new hierarchical immune-aware router, settled in advance, validated on a third untouched cohort, which is a new method study rather than a rescue of Router V2 |
| Recurring biological examples | **Non-genetic program sequence validated; genetics remains provisional** | GNMT is the provisional locus-to-state bridge; MAT1A and CYP2C19 define the severity-linked metabolic state. The Stromal ECM/Ductular injury contrast now has aligned donor association, ambient recalibration, descriptive cross-lineage localization, bulk transport, and matched-null spatial states; IGFBP7 and BICC1 remain member-gene provenance only for these program analyses | Recheck all five loci against the adopted 2026-08-17 COLOC release, then adjudicate the worked examples in one synchronized release |
| Single-cell ATAC context | **Underlying v3 package remains candidate-only; its variant-consistent primary-pair context now feeds canonical Figure 5D** | The isolated v3 candidate contains a 676,658-peak condition-blind atlas, exact donor fragment-record counts, matched DA, all 117 program projections in both cohorts, official chromVAR, and `GENETIC_READY`/`FULL_READY` completion markers. The replay reproduced all **816 gene-study pairs over 462 unique genes**, exporting 5,362 signal pairs and 22,356,697 variant posterior rows; liftover and aggregation over five shared lineages passed job `20948527`. Canonical 5D uses only the recorded primary-pair table and keeps direct MASLD/liver-fat (n=29) separate from liver-enzyme pairs (n=787) | Retain the candidate package and its amendment record without treating static peak overlap as target assignment or functional validation. The per-variant posterior vector has no adopted per-variant reference and inherits the documented replay-tolerance limitation |
| Inherited genetic variation | **Corrected rerun PROMOTED to canonical 2026-08-17; main Figure 2 A-I canonical** | The adopted release in `GWAS/finemapping/results/susie_coloc/`: tier-1/2 multi-signal 462 (was 473), single-signal 822, union 1,013, full portfolio 745. Canonical Figure 2 A-I was produced by `fig2_round2_promotion_20260817T211348Z` and job `20879541`: exact page boxes; a 1.488-by-1.970-inch portrait plot area in B; 36/73/107/455/342 consequence counts; a 352/70/26/14 ancestry partition; and all 13 genes supported in exactly four ancestry panels in F. Structural, source, font, label, raster-equivalence, visual, and scope checks passed | Rebuild remaining downstream consumers and target-transcript biotype. Supplementary E, K, and N remain explicitly withheld rather than carrying unsupported or unlinked historical results |
| Evidence classes | **Waiting on the rebuild, no longer on compute** | Nothing yet; the July classes are stale and the adopted classes do not exist | Rederive from the adopted COLOC and the bulk release |
| Spatial context | **Validated candidate awaiting synchronized release** | Complete 117-program observability, the two-program candidate family, label-blind HMSMA coverage, and categorical portal requirements | Keep HMSMA clinical inference withheld; rebuild genetics-dependent spatial classes from the adopted 2026-08-17 COLOC release |
| Figure 5 context | **All six panels and the vector assembly adopted as working canonical outputs on 2026-08-24** | 5A uses current display-set wording; 5B is the prioritized-gene protein triage; 5C is explicitly selection-conditioned; 5D is variant-consistent across five shared ATAC lineages; 5E separates assay-native program tracks without a composite score; and 5F pairs all-cell-adjusted matched-gene calibration with an illustrative outcome-blind Guilliams et al. H&E overlay. Candidate v6 passed 80 checks in job `21068950`; post-copy validation passed job `21068987` | Carry the hash-pinned six-panel set into the synchronized Resource release without upgrading peak overlap to target assignment or the spatial overlay to named histologic localization |
| Figure 3G genetics–state interface | **Adopted working main panel** | `fig3g-promoted-coloc-2026-08-17-v17-equal-points` keeps the 2.45-inch-wide threshold-free scatter in main 3G at 6 pt text, without an in-panel correlation annotation. Every point is 0.65; coral-red Both points are drawn above violet Transcriptomics sig. and blue Genetics sig. points. Across 3,617 finite SuSiE pairs, rho=0.014 (P=0.40). The expression-matched TREAT sensitivity is separate Figure S3G: 30 Genetics sig.-only, 21 matched-control-only, 377 neither, and 0 both calls across 428 pairs (paired OR=1.43; P=0.26). Jobs `20871856` and `20877891` passed recompute, PDF, font, cardinality, checksum, and scope checks | Retain panel-level status until the synchronized Resource release is adopted |
| Figures 1/6 and spatial supplement | **Impact-layout candidate validated; not adopted** | Compact spatial assay-capability and Gene Catalog insets, the retained 117-program coverage heatmap, simplified HMSMA panel, and supplementary source-state matrix under `spatial-impact-figures-candidate-2026-08-11` | Do not add an adjusted coding-versus-lncRNA coverage panel until complete versioned assay eligibility exists; assemble the remaining panels after inputs are settled |
| Spatial robustness diagnostics | **Validated candidates; supplementary only** | Full source tables plus compact composition, unit/matched-null/weight, focused cell-context, and HMSMA mask panels; one bundle preserves the complete diagnostic provenance | Use the impact candidate for assembly and `spatial-diagnostic-bundle-candidate-2026-08-11` for complete provenance; do not change the two-program family or main Figure 5 calls |
| MASLD Gene Catalog | **Noncoding-aware v2 requirements independently validated; genetics-dependent entries stale** | The read-only candidate documents 19 fields and five deterministic follow-up routes without gene outcomes; the spatial requirements, state vocabulary, and portal consumers are also ready | Rebuild entries after COLOC and bulk adoption; stage the same read-only tables for Figure 6 and the portal |
| Public perturbation plans | **Complete supplementary boundaries** | Valid results and nulls in supplement/provenance | No replacement-dataset search |
| Manuscript | **Working draft** | Six-figure Resource story and bounded claims | Reconcile prose and legends to the synchronized release |
| Cas13 screen | **Separate future paper** | Resource-to-screen prospective rationale only | No screen result may feed back into this paper |

## Active compute

**The corrected COLOC array is finished and adopted.** Array `19648525` held
1,100 study-by-chromosome tasks and every one has reached a terminal state
(`sacct -j 19648525 -X`, read 2026-08-17):

`1,100/1,100  ██████████  1,057 COMPLETED · 42 FAILED · 1 TIMEOUT`

The 42 scheduler-failed tasks have structurally complete result files in the
read-only ATAC handoff audit (40 of them from a 35-minute window on 2026-08-11
when `06_susie_coloc.R` was transiently corrupted): 17 have no exceptional rows
and 23 contain only one anonymous, low-support `abf_only` row with no SuSiE
posterior, so none can enter the primary SuSiE replay. The read-only audit is
under `atac-context-v3-candidate-2026-08-11-r1/upstream_monitor/coloc-failure-audit-2026-08-12-r2/`.

⚠ **Task `19648525_176` hit TIMEOUT after 3d18h on 2026-08-14 and is not covered
by that audit.** Its disposition is unrecorded, so the ROADMAP box "all 1,100
tasks have a terminal, audited disposition" cannot yet be ticked.

Running as of 2026-08-17 10:20 EDT (`squeue -u $USER`):

| Job | Name | Resources | What it is |
|---|---|---|---|
| `20871418` | `cabA2bnull` | 8 CPU / 140 GB / cpu | Permutation null, ComBat arm; perm 5/20 at last read |
| `20871423` | `ambient` | 16 CPU / 200 GB / cpu | Stromal ECM program ambient sensitivity rescore (`dec_depth`); registry and membership hashes verified against READY, stored BH family reproduced exactly |

```bash
squeue -u "$USER" -o '%.10i %.9P %.24j %.2t %.10M %R'
sacct -j 19648525 -X --format=JobID%20,State,Elapsed,End%20
```

The adopted cache optimization changes only the `lambda_s_locus` diagnostic at
approximately floating-point precision; all tested PP values and reported
fields were bit-identical. The pre-cache products do not need to be rerun, but
the final release must record the two execution code versions. See the dated
technical note in the archived progress index.

## Figure readiness

| Figure | Status | Remaining work |
|---|---|---|
| **1. Multimodal MASLD Resource** | Resource architecture, two-map interface, seven-state vocabulary, spatial capability matrix, and experiment-routing grammar validated; coding-versus-lncRNA adjusted matching excluded by the source check | Lead with inherited variation, multicellular remodeling, and their cell/tissue context; retain coverage and evidence states as interpretive guardrails |
| **2. Inherited genetic variation** | Main A-I is canonical. Early and late supplementary sets were independently rederived and validated as candidates; supplementary E, K, and N are explicitly withheld rather than rendered from unsupported adopted inputs | Add the target-transcript-biotype panel and adjudicate supplementary adoption separately; the canonical RORA/FABP1 examples are EUR and support method-resolution evidence, not ancestry-specific biology |
| **3. Bulk remodeling and genetic overlap** | Five-cohort synchronized analysis passed manifest, design, and complete-family checks; GNMT–MAT1A–CYP2C19 is the fixed RNA–stage–protein example; render job `19972754` failed on 2026-08-13 and is **closed, not pending**: `figures/candidates/pi-figure-redesign-2026-08-12-v2/FAILED_RENDER_JOB_19972754.txt` records it as `failed_render_provenance_only` ("data.table rejected a scalar axis label in the supplementary composition grouping expression") with the disposition **"Do not review or promote. The validated analysis was copied byte-for-byte to the 2026-08-13 v3 candidate before rerendering."** That lineage ran v3→v8; job `20161528` passed PDF, checksum, cardinality, composite-proof and scope validation | **Do not resubmit `19972754`.** Visually adjudicate the rendered, `VALIDATED` candidate `figures/candidates/aging-inspired-resource-rebuild-candidate-2026-08-15-v8/` (15 PDFs, validation_report.tsv, manifests, sessionInfo); 3F is no longer waiting now that the release is adopted, but must be built from the adopted release |
| **4. Single-cell programs and multicellular context** | The user-approved working selection retains the atlas UMAP and Hotspot heatmap, donor-paired descriptive lineage localization, fixed-member bulk tissue-state transport, focal continuum trajectories, the outcome-blind molecular-systems continuum map in Figure 4G, and the complete system-level stage-versus-continuum comparison in Figure 4H. Complete transcript, pathway, labeled Hotspot-window, NMF, and windowed systems comparisons are adopted working supplements. The ambient audit stays in source-table provenance | Retain the 11 T-cell ambient results as untestable and reconcile the adopted working panels in the synchronized release audit |
| **5. Molecular/physical context** | Six canonical panels and `fig5_molecular_context.pdf` are adopted from validated candidate v6; 5D is variant-consistent and 5E/5F retain assay-native program-level inference | Retain the panel set unchanged for the synchronized Resource release audit |
| **6. MASLD Gene Catalog** | Candidate schema, machine-readable evidence states, human-readable data dictionary, molecular-object panel, deterministic next-experiment routes, and spatial Catalog inset validated | Recompute genetics-dependent entries, then assemble worked examples from the same release without changing accepted calls |

## What is still missing

1. **Corrected COLOC propagation.** The rerun is complete and adopted
   (2026-08-17), so this is no longer a compute dependency. What is missing now is that
   some downstream consumers may still carry July values. **This list is
   measured, not assumed — several of its original entries proved already
   current.** Verified on 2026-08-20/21:
   - **Atlas: was already current.** 0 of 27,187 genes differed from the
     adopted release before any rebuild. The apparent "706 vs master 745" gap
     is a gene-universe difference (38 adopted SuSiE-positive symbols are
     absent from the atlas universe), not staleness. The chain
     `27a→45a→75→217→46d→46e→27b` was rerun anyway on 2026-08-20 (jobs
     `20951364`-`20951370`), refreshing 12 `multi_evidence/` products.
   - **Number ledger: already current.** Every 473/736 in
     [RESULTS.md](RESULTS.md) sits in an explicitly superseded row; 462/745 are
     the stated values. There is no `docs/manuscript/NUMBERS.md` at the path the
     adoption note names.
   - **`relayreg`/`orientsub`: inputs verified equivalent.** The master
     aggregate and 8/8 sampled shards are byte-identical between
     `susie_coloc_rerun/` and canonical, and the relay reads shards only, never
     `gene_level_coloc.csv` (which does differ: `b1c890c4…` vs `fcb593df…`).
     Their outputs are stale only because they *ran* on 2026-08-15, before the
     08-16 eQTL merge. ⚠ Regeneration is **not scheduled**: steps 00-07 head a
     64-step pipeline whose steps 14 and 22-64 are Cas13 Stage A execution,
     guide design, and collaborator/clinical outreach — the separate-paper
     boundary. Needs an explicit scope decision before rerunning.
   - Still genuinely open: `fig4a_overview_candidates.py`, evidence classes,
     Gene Catalog genetics fields, and the worked-example posteriors. The latter
     two additionally wait on the bulk adoption, not on COLOC.
   - The "~30 consumers" figure understates *references*: 234 live code files
     mention `gene_level_coloc` or `results/susie_coloc` (GWAS 80, scripts 75,
     RNA-seq 48, Analysis 20, other 11), excluding archived snapshot copies. A
     reference is not an active consumer, and a token/mtime sweep cannot tell
     them apart — only a per-consumer value-level join against the adopted
     release can. Do not trust a blanket staleness list.
   Rebuilt already: the portal ancestry export, the Figure 1 GWAS panels,
   canonical Figure 2 A-I, and independently validated early/late supplementary
   Figure 2 candidates.
2. **One synchronized bulk adoption.** The updated five-cohort result is the
   intended paper substrate, but stage, composition, figures, and ledgers must
   be regenerated from the same read-only inputs.
3. **HMSMA metadata.** Clinical, histologic, donor, and MALDI claims remain
   unavailable until an authoritative array-to-donor and registration table arrives;
   this does not hold up label-blind coverage or the rest of Figure 5.
4. **Final release assembly.** The July release is historical and incomplete for
   the present Resource requirements.
5. **Noncoding genetics synchronization.** The lncRNA strand rule and a complete
   DNA-only candidate are resolved. Target-transcript biotype, final examples,
   and genetics-dependent Catalog fields are no longer waiting on compute; the
   DNA tables must now be rebuilt from the adopted 2026-08-17 release.

## Canonical DEG threshold change (2026-08-12)

The canonical bulk DEG definition changed from the TREAT interval-null threshold
(`treat_fdr < 0.05` at `lfc = 0.25`) to the conventional `padj < 0.05` with
`|log2FC| > 0.50`, for comparability with the published RNA-seq literature.
On the current named release this moves the count from 1,918 to **1,853**
(1,438 up / 415 down). The genetic non-overlap that accompanies it depends on
which COLOC release is in the denominator, so quote it with its release:

| | July COLOC (447 / 14,931) | adopted 2026-08-17 (437 / 14,920) |
|---|---|---|
| non-overlap, canonical threshold | 94.85% | **95.4%** |
| overlap odds ratio | 0.690, p = 0.094 | **0.616, nominal p = 0.037** |
| same test, fragment substrate | — | 0.661, p = 0.101 |

The earlier figure of 0.752 (p = 0.214) in this section predated the 2026-08-13
background correction and is retired. **Do not read the nominal p = 0.037 as a
finding**: it is one test in a family, it is not significant on the intended
fragment substrate, and the expression-matched ladder is flat, so colocalized
genes must still not be described as buffered or constrained. See the DEG-threshold
sensitivity section of [RESULTS.md](RESULTS.md).

- The single definition lives in `load_figure_data.R::is_canonical_deg()`.
  `is_dream_deg()` is now an alias for it, so scripts calling the shared helper
  follow automatically. `is_treat_deg()` exposes the retired threshold for the
  sensitivity arm. Do not re-derive a threshold inline.
- The `treat_*` columns stay in `canonical_deg_results.csv`. They back the
  calibrated comparator and the recorded F_five and lncRNA candidate
  disagreement assertions; deleting them breaks those checks.
- **Current candidate state:** the fragment-substrate primary-threshold count is
  validated. Figure redesign v2 stopped during supplementary composition
  rendering; v3 render job `20161509` failed closed after identifying that the
  UMAP renderer had not applied the two prespecified low-cell GSE189600 library
  exclusions. Candidate v4 validated that census and then failed closed on a
  program-renderer syntax error. Candidate v5 job `20161525` rendered every
  panel that was not withheld but failed closed in Ghostscript log decoding. Candidate v6
  completed PDF, checksum, cardinality, composite-proof, and scope validation
  in job `20161528`; manual proof review identified three label collisions.
  Candidate v7 resolved those collisions; candidate v8 narrows the final 4D
  legend to prevent edge clipping. Canonical PDFs remain
  unchanged until review and adoption.
- Deliberately left on other thresholds, by design, not oversight: the
  master-regulator and drug-target panels (`lfsr < 0.05`, so small-effect
  targets such as THRB survive) and the Cas13 library (its own analytic TREAT
  at `lfc = 0.2` for library sizing).
- RORA was convergent under TREAT (logFC −0.381) and fails the 0.50 floor. It
  must no longer be described as a convergent hit.

## Open defects

| Defect | Location | Effect |
|---|---|---|
| Release-consistency checker requires two output files the current official release does not contain | `scripts/manuscript/check_release_consistency.py:78` requires `strict_drug_verdict.json` and `credible_set_coding_summary.tsv`; both exist only in the July 10 release, not the July 15 one | The checker fails on the release it is meant to certify. The release is probably correct — drug calibration and the coding-architecture fraction were both deliberately demoted to supplement — so the required-file list is stale and should follow the release rules. |
| Two rendered stage panels carry superseded numbers | `figures/main/fig3_bulk_transcriptomics/_legacy/2026-08-11_pre_fragment_stage/` | Moved out of the main panel directory on 2026-08-11; see the stage section of RESULTS. |
| Source comment quotes a retired diagnostic margin | `GWAS/finemapping/src/06_susie_coloc.R:43-45` still justifies the LD cache with a 2.5e12× `lambda_s` margin; the true value is 1.66e9× | Overstates the safety margin by ~1,500×. Corrected in `docs/archive/progress/2026-08-11_coloc_ld_cache.md`; fix the comment when that file is next touched. |
| One ledger row has no recorded producer | `docs/RESULTS.md`, the TREAT-comparator overlap row (33 / 92.5%) | The 2026-08-17 interface recompute emits the canonical threshold only, so this number is not reproduced by any run on the adopted release. Recompute or drop it before it enters prose. |
| The finalized ATAC replay submitter assumes a node-local scratch parent that does not exist | `Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/slurm/08_genetics_batch.sbatch` calls `mktemp -d /scratch/$USER/...`, but `/scratch/jameslee` is absent on the compute nodes and `mktemp` does not create parents | Every replay batch dies on its first command. The file is hash-pinned in `genetics/implementation_validation/source_manifest.tsv`, so it was not edited; batches were run on 2026-08-19 through `candidates/atac-v3-genetics-upstream-manifest-20260819T144901Z/01_replay_batch_wrapper.sbatch`, which re-checks the recorded sha256, creates the directory, and `exec bash`es the recorded script unmodified. Fix the recorded script and record its new hash, or keep the wrapper as the documented submission path. |
| ATAC replay reproduction tolerance was amended after it was finalized | `Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/22_validate_replay_batch.py:130` now asserts `rel_tol=1e-4` (was `1e-7`); `GENETICS_IMPLEMENTATION_READY` and `genetics/implementation_validation/` still carry the pre-amendment hash `31837b0c…` rather than `d5ad5db2…` | Non-EUR strata fall back to rebuilt 1000G Gram LD panels that are not positive semidefinite, so `susie_rss` zeroes negative eigenvalues and the adopted run is not bit-reproducible: across the complete set of all 816 validated pairs the median difference is 5.4e-12, the max relative 6.6e-5 and max absolute 4.5e-5, with 42/816 over 1e-7 and 0 over 1e-4 (an earlier 116-pair mid-run sample reported max 1.9e-5 and is superseded — it understated the tail 3.5×). The tolerance is bounded by the selection margin — the adopted pair closest to the 0.5 cut is 0.503716, **82×** the largest absolute difference — so no gene-level membership call can move. ⚠ Two residual risks: the observed max is only 1.5× below the 1e-4 assertion, and the assertion covers pair-level `PP.H4` only — the per-variant `SNP.PP.H4` vectors that Figure 5D consumes have no adopted reference to check against at any tolerance. Recorded in `REPLAY_TOLERANCE_AMENDMENT.md`, `replay_tolerance_amendment.tsv`, and `replay_tolerance_amendment_hashes.tsv` in the candidate root, following the `PEAK_CALL_AMENDMENT.md` convention. `20_validate_genetics_implementation.py:140` refuses to overwrite an existing hash record, so re-recording requires deleting the 2026-08-12 fixture record and explicit approval. |
| ATAC liftover unpacking was vectorized after the same finalization | `Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/08_liftover_variants.R:126-149` replaced a per-element `lapply` over the `liftOver` `GRangesList` with `elementNROWS`/`unlist`/`cumsum`; sha256 `90ac8e71…` → `cc74fea8…`, and the recorded hash is still the pre-amendment one | The loop cost 5.62 ms per unique variant (962 s for 171,301 real variants) and retained 1,845,061 single-row `data.table`s before `rbindlist`. Over 22,356,697 posterior rows that killed job `20918763` OUT_OF_MEMORY at 30.7 GB of 32 GB after 15:03:40 with no output written. The replacement is `identical()`-verified against the recorded implementation on synthetic edge cases and on 171,301 real variants, and is 50,635× faster; its one uncovered path is multi-mapping variants, absent from the real-data test and covered synthetically only. Recorded in `LIFTOVER_VECTORIZATION_AMENDMENT.md` and `liftover_vectorization_amendment.tsv`. Rerun is job `20948247` at 128 GB on `--qos=interactive`. |
| Genetic figure page-size assertion contradicted its own renderer | `11_validate_genetics.py:275-276` asserted 388.8 x 201.6 pts within `abs_tol=0.2`, but `13_render_genetic.R:110-111` renders via `cairo_pdf`, which truncates the page box to whole points and wrote 388 x 201 | The check could never pass on a cairo host; no fixture rendered a real PDF, so it surfaced only on first real-data execution (job `20948356`). Amended to `abs_tol=1.0` (1/72 inch), recorded in `FIGURE_PAGESIZE_AMENDMENT.md` and `figure_pagesize_amendment.tsv`, sha256 `78376e1c…` → `503e98be…`. The renderer was deliberately left on `cairo_pdf`: base `pdf()` emits exact fractional points but risks the Type 3 fonts the next assertion rejects, and the recorded non-genetic validator asserts no numeric page dimensions at all (`10_validate_non_genetic.py:375`). |
| Retracted probe run wrote real readiness files | A probe copy of `11_validate_genetics.py`, differing only in the page-size tolerance, was run to enumerate remaining failures; it passed and wrote `GENETIC_READY`, `FULL_READY`, and `genetics/validation/` | Those readiness files were produced by an unmanifested ad-hoc file and did not record the code that produced them. They were moved to `candidates/atac-v3-genetics-upstream-manifest-20260819T144901Z/superseded_probe_run_20260820T192500Z/` and superseded by the recorded rerun `20948527`. Comparing the two runs afterwards showed every data file byte-identical, so the retraction cost nothing scientifically, but a validator that writes readiness files must not be used as a dry run. |

## Work that can proceed now

- Rewrite introduction, discussion, and non-numerical Results around the
  Resource contribution and claim boundaries.
- Finalize Methods sections whose inputs and models are already settled.
- Prepare figure layouts, legends, and data dictionaries with placeholders for
  values that are explicitly not yet available.
- Rebuild bulk-dependent stage and composition analyses in isolated candidate
  paths.
- Retain the adopted program-level Figure 5E/5F panels; positive Moran
  direction alone is not the main spatial result.
- Draft supplementary legends for the validated composition, unit/null,
  cell-context, and HMSMA proxy-mask diagnostics without upgrading their claim
  boundaries.
- Audit portal/Gene Catalog schemas without publishing genetics-dependent entries;
  the spatial candidate products are ready for synchronized integration.
- Retain the finalized ATAC v3 replay and canonical Figure 5D against the adopted
  2026-08-17 COLOC manifest; do not convert peak overlap into target assignment.
- Use the validated GENCODE identity, pooled biotype subdivisions,
  fragment-native stage counts, and prespecified-program lncRNA annotations in
  candidate figures and Methods; do not adopt them before release
  synchronization.
- Reconcile references and competitor positioning.
- Integrate the validated histology-anchored continuum qualitatively as a
  complementary molecular-resolution layer; keep exact candidate values tied
  to the checksummed result ledger until synchronized adoption.

Exact genetic counts and overlap percentages are now authorized **from the
adopted 2026-08-17 release only**, as ledgered in [RESULTS.md](RESULTS.md); any
July value is superseded. Do not settle new gene examples or Gene Catalog classes
until their own rebuilds against that release have run and been checked.

## Known interpretation boundaries

- The updated bulk analysis contains **844** pooled disease-vs-control samples;
  the complete corrected fragment source has **1,257** QC-passing samples.
  These are different denominators for different purposes.
- The prior six-cohort fragment-native stage candidate passed its own checks but
  is provenance only. Figure 3 requires the active five-cohort synchronized
  rebuild and one shared release ID.
- The single-cell analysis must use biological donors, not sequencing runs or sort
  fractions. Cross-cohort inferred fibrosis stage is not authorized.
- NMF k4/k6 factors are continuous interpretive axes, not reproducible discrete
  patient subtypes.
- The three disease aspects of metabolism, inflammation, and fibrosis are **not
  resolved as separate outcome axes**. NAS enters as a graded sum and as an
  adjustment set; fibrosis stage is the severity anchor. Separability of the
  three NASH-CRN components is now **`tested_negative` in GSE267145**: tested in
  99 participants against a prespecification recorded before the analysis, the
  components span a participation ratio of 1.42 of a possible 3.0 against a 2.0
  floor and behave as one axis rather than three. Aspect **attribution stays
  `indeterminate`**, now for a measured reason rather than for want of a test;
  nothing can be assigned to a single aspect on this evidence. The result
  settles one cohort and one assay and does not settle the SAF partition of
  steatosis, activity and fibrosis (SAF grades in GSE202379 and PXD051911 remain
  unused; GSE135251 carries only the NAS sum).
- **Variant combinations are untested** beyond a single characterised local
  haplotype at GNMT, which is reported as unresolved. The colocalization design
  is single-variant and does not test combinations or interaction.
- **Chromatin supplies context, not effect.** The large majority of jointly
  testable peaks resolve to `indeterminate` and a minority to `source_dependent`.
  No epigenetic effect or causal claim is available.
- Yakubovsky is a normal liver zonation/lipid-context reference, not a MASLD
  disease-validation cohort.
- Vu arrays remain unresolved physical units; directional replication is
  `source_dependent`, not donor-level inference.
- HMSMA contributes label-blind coverage and descriptive per-array organization
  only until the authoritative metadata and image-registration joins arrive.
- GSE287826 cannot support donor-level inference without a public donor key.
- No currently usable spatial cohort supplies every desired MASLD control,
  donor, histology, and metabolite join.
- Proteomics is not applicable to lncRNA molecules. Chromatin observes a locus
  or regulatory context, not the RNA molecule itself.
- ATAC v3 separates static observability from MASH-associated accessibility.
  Its four strictly replicated DA peaks are stellate-lineage chromatin findings,
  not causal-gene or causal-cell-type proof. The GSE281367 T/NK label remains a
  combined, source-dependent observability context.

## Source of truth

- Scientific scope and claims: [PAPER.md](PAPER.md)
- Current numeric ledger: [RESULTS.md](RESULTS.md)
- Execution sequence: [ROADMAP.md](ROADMAP.md)
- Detailed release history: [technical/RELEASE_PROVENANCE.md](technical/RELEASE_PROVENANCE.md)
- Full superseded handoffs and dated logs: [archive/INDEX.md](archive/INDEX.md)

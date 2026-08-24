# Review: MASLD Multimodal Model Benchmark plan

**Date:** 2026-08-21 · **Reviewed against:** repo state, GEO/PRIDE records, live SLURM, local reference files.

## Verdict

The machinery is right and the arithmetic is not. Sealing, leakage control, mandatory baselines,
and the preregistered null are better than most benchmark papers ever get. But six of the seven
tracks are disqualified from producing a champion by the plan's own rules, and the one that
survives reuses pass/fail floors this project has already failed twice. Expected outcome after
~120 model adapters is zero champions.

## Ranked issues

**1. One eligible track, with a failing in-project prior.** Variant-to-regulation: "no champion"
without sumstats. RNA-conditioned ATAC: retrospective without FNIH. Graph: "no champion without a
sealed source." Perturbation: exploratory. Unified: "no unified champion is possible." Spatial:
deferred. Bulk: dead on power (#4). That leaves cell-state mapping, whose gate (macro-F1 ≥0.70,
every class F1 ≥0.50) is numerically identical to the floors that terminated ROADMAP WP6a: V42
scored 0.2756 on untouched GSE212837; Router V2 hit 0.7071 overall on GSE296875 but failed
per-class at B cells 0.2888 and NK-T 0.4229. The plan never cites this.
→ Cite WP6a as the pilot with those numbers. The binding constraint is rare immune-class F1, not
the overall floor, so either argue why 70 models change that specific behavior or swap the
conjunctive pass/fail gate for a preregistered per-class estimation endpoint that is interpretable
either way. Drop "champion" as the deliverable.

**2. GSE289173 cannot be split into two seals, and its outcomes are not on GEO.** The GEO deposit is
one `GSE289173_RAW.tar` of MTX/TSV matrices, no eQTL/ieQTL summary statistics, no genotypes; the
outcomes live in a single 6.6 GB Zenodo archive (10.5281/zenodo.14586466). Pulling task-1 cell-state
labels puts task-2 sumstats in hand at the same moment, which contradicts "separately sealed" and
"one-time sealed inference."
→ Inventory the Zenodo archive before locking, confirm the sumstats are cell-type stratified (if not,
"no major lineage worse by >0.03" is uncomputable), and state that GSE289173 supports one
unblinding event.

**3. The seal has no enforcement and no external timestamp.** `id` returns uid=91825 only; no second
account or group; `setfacl`/`getfacl` absent; GPFS `mmputacl` is owner-revocable; no sudo.
"Independent evaluator" means an independent code module in the same uid's directory. And the one
primitive that would make it tamper-evident does not land: `git check-ignore -v` confirms
`.gitignore:8:*` swallows `locks/selection_lock.json`, `locks/predictions.sha256`, and
`locks/manifest.tsv`, so no lock ever reaches origin.
→ Add `!Analysis/MASLD_Model_Benchmark/locks/*` re-includes and make an external stamp (pushed
commit SHA, or OSF/Zenodo deposit of the lock hashes) a hard precondition of wave 9. If no second
party holds the archive, call it a prospectively-locked retrospective evaluation, not a seal.

**4. The bulk gate is unreachable by 5-9×.** At n=18 with 5 ordinal stages, simulated SD of the
paired Spearman difference is 0.08-0.14; 80%-power MDE is 0.22-0.44 rank-correlation units, and
power at the registered 0.05 gain is 0.05-0.13, i.e. at alpha. Under the exact null the 0.05
effect-size half of the gate is crossed 27-37% of the time. The plan's own "activate only after
prospective 80% power" clause therefore blocks the track by construction.
→ Delete the 0.05. Register the detectable effect (~0.25-0.35) or demote GSE267031 to a
CI-reporting descriptive consistency check.

**5. GSE296875 has donor-linked steatosis and fibrosis histology, but not an adjudicated MASLD
label.** The GSM pages expose only tissue and well, but the paper's Data S1 Table S1 contains BMI,
categorical and numeric macrovesicular steatosis, and categorical and free-text fibrosis for the
de-identified donors. The steatosis field is nonmissing for 38 of 39 analyzed donors and the
fibrosis field for 37 of 39; these counts do not mean all 38 or 37 donors are phenotype-positive.
The public fields do not include NAS, ballooning, lobular inflammation, alcohol exposure, or an
etiology adjudication, so steatosis must not be relabeled as confirmed MASLD or MASH. Donors were
pooled five per 10x inlet, so donor-grouped outer folds alone do not test well-level batch
transport. STATUS.md also records prior project use as a diagnostic-only router cohort, which
makes it project-exposed rather than a new external holdout for this paper.
→ Use the histology as development-only cross-sectional phenotype evidence, retain donor as the
inferential unit and well as a technical batch, freeze the public barcode-to-donor assignment, and
prespecify leave-one-well-out sensitivity. Do not claim disease etiology from these fields.

**6. The reference contract pins a file nothing can read.** Both SHA-256s verify exactly, and the
truncation call on `gencode_v49/` is correct (3.7 MB genome, 71 KB transcripts, both fail `gzip -t`).
But `refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz` is plain gzip (`1f8b 0808`, no BGZF
`BC` field) with no `.fai`/`.gzi`, so no sequence model can window-extract from it; the repo's 17
existing sequence scripts all read an indexed `fasta/genome.fa` that the contract does not name.
It is also the patch/ALT assembly (706 contigs; 503 of 528 GTF seqnames non-primary; 2,945 symbols
on more than one contig), so a locus duplicated onto a PATCH contig escapes the chromosome/LD-block
holdout and breaks the PredictionBundle row-ID guarantee. `contig_policy` in `reference.py:207` is
inert free text.
→ Pin the primary-only analysis set plus `gencode.v49.primary_assembly.annotation.gtf`, with
fasta + .fai + chrom.sizes + blacklist + per-model bias files each hashed, and implement an actual
contig-membership check. Add `native_annotation_release` to DatasetManifest: the project runs three
live annotations (v49 bulk, v47 in refdata/salmon, v44 for the multiome-ATAC path).

**7. Compute is budgeted; engineering is not.** Campaign approval requires GPU-hours, CPU-hours,
memory, and a hash, and says nothing about person-effort. The census is ~92 checkpoints plus ~28
baselines, each needing an isolated env and a probe/prepare/fit/predict/export adapter. You have 26
micromamba envs, 4 of which touch this census, already spanning six mutually exclusive framework
pins. GPU reality: 18 L40S + 22 B6K institute-wide, 27 pending vs 15 running right now, and your own
B6K sbatch jobs 20208726/20215473 pended from 2026-08-13 and were cancelled at 0:00 elapsed.
[Estimate] adapter/env engineering, not GPU-hours, is the 10× oversubscription.
→ Add a per-model engineering budget line and cap the census by that budget, not by scientific
completeness. State at the top what Resource work this defers, given STATUS.md's standing
"reconciliation, not expansion" and ~30 unpropagated COLOC consumers.

## Cheap fixes before census freeze

- Corgi's released path rank-maps the trans vector to a packaged 2,891-gene reference
  (`predict.py` L61-78), so raw UMI, CPM, and log-CPM give byte-identical inputs. Two of your four
  mapper arms are the same experiment. Keep length-adjusted TPM and the learned encoder. The
  "log-transformed, quantile-normalized, batch-corrected" description of Corgi's training input is
  [Unverified] from open sources; source it or cut it.
- 16 RNA-ATAC integration models are admitted with no task row, no endpoint, and no gate, and the
  latent-factor members cannot emit the only chromatin endpoint on offer. Add a retrospective
  same-nucleus retrieval row on GSE296875 or move the island to "registered but deferred."
- "scPRINT-2 small/medium": no medium-v2 checkpoint exists (small-v1, medium-v1, medium-v1.5,
  large-v1, small-v2 only). Change to small-v2 with SHA at freeze.
- GSE313774 has no resmetirom arm: all 33 GSMs carry only Condition 1/2 ± Ins across 3 batches, no
  donor or sex field. The CXCL1/IL8 exclusion guards a leak that cannot occur, and
  `expected_biological_units = 33` violates your own pseudoreplication rule (n = 3 batches).
- Device eligibility is per-environment, not global: cu13.0 builds do not run on L40S under driver
  550.107.02, cu12.4 cannot target sm_120, TF 2.15.1 has no sm_120 kernels. A cu12.8 torch build
  spans both pools. Standardize on cu12.8 and add `device_eligibility` to ModelManifest.
- The prospective-power rule fixes the variance source but not the effect estimate. State that the
  numerator uses a shrunk or cross-fitted gain (one-SE-below-max, or bootstrap-debiased), not the
  development maximum.
- Immutability has no retention tier. A 32-bp genome-wide track is ~97M bins (~19 GB per
  model-fold-seed at 100 tracks fp16; ~1.5 TB for a full 7,611-track Borzoi pass), /gpfs is 89%
  full, and `mmlsquota -j sanjana_lab` returns no fileset. Tier retention (metrics, predictions, and
  hashes forever; embeddings and tracks with a declared expiry plus a hashed regeneration recipe)
  and get the real quota before campaign 1.

## Do not let anyone talk you out of these

- Preregistering the null and publishing "no foundation model beat the baseline" as an accepted outcome.
- Mandatory strong baselines retained through final selection, with a baseline allowed to win.
- Explicit pairing topologies and missingness states, with masks instead of missing-as-zero.
- Checkpoint-exposure classification gating champion eligibility, and the frozen-census no-latest-tags rule.
- Scientific language is already compliant: donors as the unit, "cross-sectional stage-associated
  remodeling," "colocalized"/"genetically anchored candidate," no MR, no "TREAT."

## If you change one thing

Cut the census to the one or two tracks that have a real external endpoint and a reachable gate, and
budget adapter and environment engineering explicitly. A 120-model sweep whose own rules permit at
most one champion is the failure mode, not the ambition.

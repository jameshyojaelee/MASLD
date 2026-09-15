> ⛔ **SUPERSEDED 2026-08-30 — see `AGENT_HANDOFF_2026-08-30.md`.** This file lists job 21083341 as
> RUNNING; it COMPLETED 2026-08-25. Reading it as current state caused four tasks to be dispatched
> on 2026-08-30 for work that was already finished. Treat it as history, not state.

# MASLD model benchmark takeover handoff

**Takeover snapshot: 2026-08-25 09:20 EDT**

This document hands the entire MASLD model workflow to the next primary agent.
Continue only inside `Analysis/MASLD_Model_Benchmark/`. Preserve every unrelated
worktree change. Do not cancel, replace, detach, or resubmit any running or
pending SLURM job during takeover. Inspect live state first because scheduler
state and work-unit counts will change after this snapshot.

## Read order

Before changing anything, read:

1. Repository `AGENTS.md` and the five Resource authorities under `docs/`:
   `README.md`, `PAPER.md`, `STATUS.md`, `RESULTS.md`, and `ROADMAP.md`.
2. This handoff.
3. `OVERALL_PLAN.md`, the binding model campaign plan and scientific requirements.
4. `README.md`, the checkpoint-, topology-, and implementation-level census.
5. `MODEL_RESULTS_SUMMARY.md`, the living results and progress dashboard.
6. Frozen evaluator output files named in the dashboard before quoting a metric.

Do not reconstruct direction from older dated handoffs. They are history only.

## Ultimate objective

Build the strongest useful, publishable, and externally transferable
multimodal MASLD research model suite for mapping cross-sectional
liver-environment remodeling across fibrosis, NAS components, and a descriptive
molecular continuum. Systematically test and adapt leading cell, sequence,
chromatin, multimodal, perturbation, and regulatory families across independent
cohorts and assay modalities against strong task-native baselines. Outputs
should be uncertainty-aware, lineage-specific predictions of altered cell
states, regulatory activity, genes, pathways, and variant effects that rank
mechanistic hypotheses for wet-lab validation. Use donor-safe and study-held-out
validation, explicit modality topology and missingness, contamination audits,
and project-reserved external evaluation. Build one context-conditioned
sequence-plus-cell-state architecture only if complementary specialist evidence
justifies it. Release reproducible per-task best models without universal, causal,
diagnostic, prognostic, treatment-response, or other clinical claims beyond
what external data directly validate.

The target is a donor- and cell-context-aware map of remodeling, not a disease
classifier. A collection of broad cell labels, MPRA scores, or assay-specific
embeddings alone does not complete the objective.

## Non-negotiable user instructions

- Work only under `Analysis/MASLD_Model_Benchmark/`.
- Leave unrelated edits and output files untouched. The shared worktree is very
  dirty and contains multiple agents' uncommitted work.
- Leave every currently running and pending SLURM job in place. The next agent
  inherits all of them.
- Do not ask the user for routine campaign approval. Continue autonomously
  within the frozen scientific and compute requirements.
- GPU production uses `gpu`, `--qos=nslab`, and generic/masked job names. Do not
  use `innovation` for future GPU jobs unless the user explicitly changes the
  policy again.
- There is no current GPU pending-job limit. CPU and I/O jobs may also remain
  broadly parallel. The last user ceiling was 300 total jobs.
- Never manufacture RNA-ATAC pairing. GSE296875 is same nucleus; GSE244832 is
  different aliquots with an unresolved RNA-ATAC participant crosswalk.
- Cells, nuclei, spots, guides, arrays, and libraries are not biological
  replicates. Use donors or participants for inference.
- Missing assays use explicit masks/states and are never encoded as zero.
- Do not train on or optimize the frozen 117 Hotspot programs. They are read-only
  post hoc projections.
- Do not change Resource evidence classes, claims, figures, Gene Catalog, portal,
  or the separate Cas13-screen paper from model results.
- Keep held-back outcomes evaluator-only. Freeze prediction hashes, selection,
  calibration, and thresholds before one-time outcome joins.
- Restricted or exposure-unknown models can be comparators but cannot become an
  best openly licensed model.

## Scientific phase and headline state

The project is in overlapping phase 3 frozen screens and phase 4 full
specialist/native-lane screens. Phase 5 complementarity is incomplete. Phase 6
conditional architecture is not justified. Phase 7 finalist selection record and held-back
evaluation have not begun.

There is no task-held-back MASLD best-model and no universal model. The primary
GSE289173 cell-state and cell-type eQTL/ieQTL outcomes have not been opened. A
prior metadata incident exposed disease-group and sex fields for two samples,
so disease/control AUROC is ineligible as a held-back claim; the primary cell-state
and variant-effect outcomes remain isolated.

GSE296875 is not literally phenotype-free. It has donor-linked BMI, age, sex,
macrovesicular steatosis for 38/39 donors, and fibrosis fields for 37/39. It does
not establish adjudicated MASLD/MASH, NAS, standard fibrosis stage, alcohol
exclusion, etiology, or longitudinal progression. Use it for development,
histology-stratified audits, and topology-matched tasks, not a final MASLD
phenotype or clinical claim.

## Current quantitative results

The complete tables, dataset roles, progress bars, uncertainty, and metric
authority paths are in `MODEL_RESULTS_SUMMARY.md`. Update that document in place
whenever a material result or exact work-unit count changes.

| Task | Strongest completed evidence | Interpretation |
|---|---|---|
| 1,000-cell broad state | Geneformer V2 316M macro-F1 0.9522 versus tuned HVG/PCA 0.9341; gain +0.0181, 95% donor-bootstrap CI +0.0044 to +0.0320 | Positive but misses the +0.02 check; broad identity only |
| 50,000-cell broad state | TF-Sapiens MLP macro-F1 0.938875, Brier 0.091985; linear 0.926011, Brier 0.108977 | Matching baseline and other models incomplete; no winner yet |
| Cell-model complementarity | Geneformer plus baseline stack gains only +0.0009 | Does not justify a context architecture for broad identity |
| MPRA variant effect | HyenaDNA delta 0.216223 versus allele ridge 0.134216; gain +0.082008, 95% block CI +0.027427 to +0.138918 | Strongest positive regulatory development result; not eQTL transfer |
| MPRA secondary | Sei 0.193865, gain +0.059649, CI -0.003487 to +0.119520; restricted NT delta 0.187921, positive gain CI | Promising/restricted; mandatory raw controls incomplete |
| Observed multiome | Observed-ATAC GLM deviance skill 0.162667, two-way CI 0.151999-0.173259 | Strongest completed accessibility control; it uses observed ATAC |
| RNA-only ATAC | MultiVI improves deviance 0.794%, 3/5 lineages; RNA-only factorized control 0.017609 with CI crossing zero | No strong RNA-conditioned gain yet |
| Regular Corgi head | Relative deviance skill 2.53e-7, effectively zero | Do not scale this rung; stronger FiLM rung is queued |
| Histology source | GSE267145 H3 SVM stage macro-F1 0.691; RNA elastic-net fibrosis F1 0.586; RNA centroid NASH-CRN component-sum rho 0.776 | Participant-level development, not lineage-specific |
| Histology transfer | Frozen RNA SVM macro-F1 0.419 on 16 GSE260666 participants, 95% CI 0.222-0.652 | Weak external development; no adoption |

Other important negatives are retained: Enformer SAD/SAR underperform the MPRA
allele baseline; DNABERT-2, gkm-SVM, and deltaSVM do not adopt; PeakVI's first
training-context smoke is negative; scBasset's preliminary screen is negative;
scGLUE and StabMap do not beat linear CCA retrieval; multimodal late fusion does
not lead the GSE267145 registered endpoints.

## Exact live SLURM instruction

**Leave all jobs below untouched.** Do not cancel them during takeover. Do not
submit replacement copies merely to shorten a dependency wait. First inspect:

```bash
squeue -u "$USER" -o '%.18i %.26j %.10T %.10P %.10q %.12M %.12l %R'
sacct -j <jobid> --format=JobID,JobName%28,State,Elapsed,MaxRSS,ExitCode -P
scontrol show job -o <jobid>
```

Snapshot at 2026-08-25 09:20 EDT: 42 visible GPU jobs, one running and 41
pending. The running GPU job was Cobolt `21100294`. Three GPU roots were pending
on priority; the other 38 were held until its dependencies are ready. Queue counts are volatile.

### Running model/data jobs

| Job | Purpose | Snapshot state | Output/notes |
|---|---|---|---|
| `21100294` `model-work-212` | Cobolt 25-fit GSE296875 model bundle | RUNNING gpu/nslab; 6h limit | `executions/model-work-212-21100294.{out,err}`; 150/150 CPU controls already complete; 0/25 verified model receipts at last refresh |
| `21097161` `model-cpu-train-503` | Coherent 50k classical cell-baseline rerun | RUNNING cpu/nslab; 10/15 fold receipts | `executions/model-cpu-train-503-21097161.staging`; scorer `21097344` waits afterok |
| `21100574` `model-dispatcher` | Central GPU bundle dispatcher | RUNNING cpu/nslab | Leave running; GPU bundles are already submitted |
| `21083341_0..7` | GSE268273 RSEM quantification array | RUNNING cpu/nslab | `slurm/quantify_gse268273_rsem_bundles_cpu.sbatch` |
| `21091639_0,2,4,6` | GSE105127 RNA quantification | RUNNING cpu/nslab | Other tasks may already be terminal; consolidate job waits for full array |
| `20917657`, `20917658` | Long-lived user session jobs | RUNNING | Not model jobs; still leave untouched |

### Independent pending GPU roots

| Job | Purpose | State |
|---|---|---|
| `21083052` | UCE 4L/33L 50k frozen screen | PENDING Priority |
| `21100293` | scVI/scANVI 50k study-held, 15 logical tasks | PENDING Priority |
| `21100295` | 50k common cell-head bundle using frozen embeddings | PENDING Priority |

### Serialized GPU chain after Cobolt

The current queue intentionally serializes several large bundles. It uses
`afterany`, so every wrapper must verify its required predecessor output files and
fail closed; dependency satisfaction alone is not evidence of predecessor
success.

| Dependency | Job | Purpose |
|---|---|---|
| `21100294` -> | `21100296` | 50 raw-sequence CNN/transformer MPRA controls |
| `21100296` -> | `21100297` | Regular Corgi FiLM-plus-head rung |
| `21100297` -> | `21100298` | Borzoi converted full-window probe |
| `21100298` -> | `21100299` | scooby full-backbone synthetic probe |
| `21100299` -> | `21100300` | scBasset three-seed rectangle bundle |
| `21100300` -> | `21100301` | scBasset five-seed extension |
| thereafter | `21100302` through `21100331` | Sequence task-native five-seed bundles for BPNet, full ChromBPNet, CNN, and transformer across lane-serialized chains |

The frozen sequence completion check `21101029` waits afterany on
`21100328`-`21100331`. It must emit a terminal non-ranking disposition if any
rectangle cell is missing, failed, zero-filled, or dropped. Do not open partial
metrics.

### PeakVI and GSE296875 transfer chains

| Job/dependency | Purpose |
|---|---|
| `21100328` -> `21101210` | PeakVI GSE296875 source fit and GSE281367 cross-cohort revision-2 bundle; currently pending dependency |
| `21101210` -> `21101384` | TF-Sapiens GSE296875 query embedding extraction; dependency is lane serialization, not scientific coupling |
| `21097161` -> `21101362` | Classical Atlas-fit predictions on GSE296875 |
| `21097161` -> `21101446` | Final GSE296875 transfer input inclusion |
| `21101384` -> `21101402` | TF linear/MLP Atlas-fit predictions on GSE296875 |
| `21101362`,`21101402`,`21101446` -> `21101405` | Prediction hash record |
| `21101405` -> `21101410` | Evaluator-side donor-balanced GSE296875 scoring |

The TF transfer path corrected a pre-execution defect: activation-fold metadata
differed from the authoritative study-held split for 42,584/50,000 cells. The
predictor now ignores that stale field, binds split SHA `10927e16...`, and
replays both heads on 21,572 fold-0 held cells within maximum probability error
1.129e-6 linear and 5.413e-7 MLP. Verification output file:
`executions/model-check-337-21101451`.

### CPU acquisition dependency chains

| Chain | Purpose |
|---|---|
| `21083341_*` -> `21083342` -> `21088737` -> `21089146` -> `21096845` | Consolidate, admit, independently audit, and freeze GSE268273 development readiness |
| `21091639_*` -> `21091641` -> `21096478` -> `21096844` | Consolidate and freeze GSE105127 development readiness |
| `21097161` -> `21097344` | Score the complete 50k classical baseline only after all 15 folds succeed |

## Current rectangle counts

These are outcome-blind counts from the last frozen census. Recompute from the
completion check before changing them.

| Lane | Complete | Remaining |
|---|---:|---:|
| BPNet, ChromBPNet, CNN, transformer task-native fits | 62/500 | 438 |
| Task-native evaluated views | 81/625 | 544 |
| BPNet | 20/125 | 105 |
| Full ChromBPNet | 19/125 | 106; 38/250 views complete |
| scBasset fixed rectangle | 4/25 | 21 |
| Corgi FiLM fold-seed fits | 0/25 | 25 |
| Cobolt model fits | 0/25 verified at snapshot | GPU bundle running |
| 50k classical cell folds | 10/15 | 5 |
| GSE296875 transfer candidates | 0/7 scored | Entire prediction/lock/evaluator chain pending |
| Held-back tracks | 0 | No SelectionLock or co-unblinding |

## GSE244832 PeakVI: unsafe/incomplete handoff state

**Do not queue or run the production GSE244832 PeakVI bundle yet. Do not score
GSE244832 outcomes.**

Scientifically valid target materialization is complete:

- Inclusion job `21101456`; I/O job `21101459`.
- 3,049 eligible ATAC-only cells.
- 48 observed donor-lineage units; four structural and 20 below-QC units remain
  explicit missing/NaN.
- Valid cell SHA `1a5d40...54be4`; test SHA `91cb4a...d3df9`; work SHA
  `493f60...e2b89`.
- No RNA, outcome, or metric was read.

The inference-only runner itself was compute-tested before the takeover edits:

- `21101584` complete; output file SHA `c18675...19845`.
- Runner SHA `daf9fee44e6efcdf7cf6124b3c90a23414d15f90923a418748c68fe760ddb7f8`.
- It contains no optimizer or refit path.
- Synthetic generic commit integration `21101585` completed, but it predates the
  new PeakVI-specific committer and does not validate that new layer.

Latest untracked in-flight files:

| File | SHA-256/status |
|---|---|
| `scripts/peakvi_gse244832_inference_only.py` | `daf9fee...db7f8`; runner compute-tested before specialized-commit changes |
| `scripts/commit_peakvi_gse244832_atac_transport_predictions.py` | `66a7e6c...1c983`; new; AST only, not compute-tested |
| `tests/unit/test_commit_peakvi_gse244832_atac_transport_predictions.py` | `dbf7fb7...85a1`; new; not compute-tested |
| `slurm/run_peakvi_gse244832_inference_only_bundle.sbatch` | `fc2c99e...9d8a6`; calls specialized committer; shell syntax only |
| `slurm/test_peakvi_gse244832_inference_only_cpu.sbatch` | `3328fcb...b549`; updated for eight tests; unsubmitted |
| `slurm/test_peakvi_gse244832_commit_integration_cpu.sbatch` | `78ce2af...3498`; builds fake frozen states and specialized commit; unsubmitted |
| `tests/unit/test_peakvi_gse244832_inference_only.py` | Existing untracked prior-agent tests |
| `scripts/commit_gse244832_atac_transport_predictions.py` | Generic committer left unchanged because active task-native jobs depend on it |

The new specialized committer is intended to verify exact source outer output file,
source fit receipts, all six state hashes and parameter hashes, strict restore,
and `fit_or_refit_performed=false` before generic bundle validation. It has not
yet been compute-verified.

Remaining missing pieces and next safe actions:

1. The final campaign config
   `config/campaigns/peakvi_gse296875_to_gse244832_inference_only_20260825.json`
   is intentionally absent. No GPU queue item exists.
2. The exact production source output file and six states do not exist until
   `21101210` completes.
3. Submit the two CPU-only test wrappers above through SLURM and verify their new
   read-only outputs. Peer-audit the specialized committer.
4. Resolve output semantics prospectively. PeakVI emits 0-1 accessibility
   probabilities; do not label them `masked_accessibility_count` or claim a
   count-only endpoint without a scientifically coherent TaskSpec revision.
5. Verify every source input payload before reading window TSVs. Ensure the
   committer itself, not only the producer, rejects target refit and wrong source
   states.
6. After `21101210` completes, verify
   `executions/peakvi-cross-cohort-fit-predict-21101210`, its exact outer
   `ARTIFACTS.json`, both fit receipts, all six state hashes, and parameter
   hashes.
7. Only then create a read-only config and CPU inclusion, add a generic
   gpu/nslab queue item, and commit predictions. Never open or score outcomes in
   the producer environment.

The task-native GSE244832 baseline lane is separate. Attempts `21101288`,
`21101327`, and `21101353` failed before source values, predictions, commits,
scoring, or outcomes were read. Latest failure was CLI `config` versus
`config_path`. Revision-4 config SHA
`4b82d55552f4feaf2f1902f3617092f3695f44898ba410733302600848e9761f`
is staged but unsubmitted. Progress remains 0/6 prediction bundles, commits,
and verifications.

## Phenotype-linked remodeling lane

This is the largest scientific gap relative to the revised goal. Existing
histology models are bulk participant-level. Existing cell models predict broad
identity. Existing GSE296875 regulatory tasks do not currently use donor
histology as a training or selection endpoint.

The next safe design should use included GSE296875 donor metadata only after
authoritative endpoint and missingness verification:

- Donor pseudobulk or donor-state summaries are the inferential units.
- Outcomes are cross-sectional macrovesicular steatosis and the source fibrosis
  field, with explicit masks; do not relabel either as NAS or standard F stage.
- Prespecify an adult-only sensitivity because five of 39 donors are under 18.
- Include age, recorded sex, BMI, and technical/study variables only where the
  TaskSpec permits them.
- Compare molecular-only, metadata-only, and molecular-plus-metadata paths.
- Keep phenotype labels evaluator-only for model paths and fit every transform
  within donor-training folds.
- Do not use the frozen 117 programs as targets, selected features, or tuning
  endpoints. They may be projected only after selection.
- Treat this as project-exposed development and an error/complementarity audit,
  not external confirmation.
- The stronger external phenotype path remains GSE267145 source development,
  GSE260666 small transfer, and future GSE268273/GSE49541/GSE83452 activation.

No verified phenotype-linked GSE296875 campaign was present at this handoff.
The brief phenotype audit subtask was interrupted for takeover before producing
a completion report; inspect `git status` and current files rather than assuming
it made no changes.

## Living dashboard requirement

`MODEL_RESULTS_SUMMARY.md` is the single human-readable dashboard requested by
the user. It is versionable through `.gitignore`. Update it in place whenever:

- a fit/prediction/evaluation denominator changes;
- a model gains a frozen metric or terminal disposition;
- a cross-cohort transfer completes;
- a missing piece changes;
- scheduler state is materially refreshed;
- a SelectionLock or held-back event occurs.

Progress bars must carry their exact numerator and denominator. Do not use them
as scientific confidence. Frozen evaluator output files remain authoritative if a
dashboard row becomes stale.

## Worktree state

Snapshot:

- Branch: `main`
- HEAD: `7ba421fea5352f1df18c7031665b0738d4eeecfd`
- Under `Analysis/MASLD_Model_Benchmark/`: 57 modified/staged paths, 1,226
  untracked paths, zero deletions at 09:20 EDT.
- No commit, branch, push, reset, revert, or cleanup was performed for this
  handoff.
- The new living dashboard and this handoff are untracked/versionable files.
- `executions/` and many generated output files are intentionally ignored. Their
  absence from `git status` does not mean they do not exist.

Do not run `git reset`, `git checkout --`, broad cleanup, or remove untracked
files. Many are active agents' source, configuration, manifests, and read-only
receipts. Use path-specific diffs and hashes before editing.

Current handoff-related changes:

- New `MODEL_RESULTS_SUMMARY.md`.
- New `AGENT_HANDOFF_2026-08-25.md`.
- `.gitignore` has whitelist entries for those two documents.
- GSE244832 PeakVI untracked files listed above are incomplete in-flight work.
- `OVERALL_PLAN.md` was updated by prior agents through the latest transfer and
  PeakVI audit state but may be ignored by Git; treat file contents and output file
  citations as current, not tracking status.

## Immediate takeover sequence

1. Read the authorities in the stated order and run live `squeue`/`sacct`.
2. Leave every inherited running and pending job untouched. Check Cobolt
   `21100294` and baseline `21097161` for terminal state and read-only receipts.
3. If Cobolt completes, independently verify all 25 model fits before evaluating
   them against the frozen 150-control matrix. If it fails, preserve the failed
   attempt and diagnose before any new revision.
4. If the 50k baseline completes, allow `21097344` to score. Update the dashboard
   only from the frozen evaluator output.
5. Compute-test and peer-audit the new PeakVI specialized committer and wrappers.
   Do not queue production and do not score GSE244832.
6. Monitor GPU roots and the sequence chains. Do not open partial rectangle
   metrics. Let `21101029` enforce completion.
7. Finish an outcome-blind, donor-safe phenotype-linked GSE296875 TaskSpec and
   fixture, then verify it on CPU before any model campaign.
8. Monitor GSE268273 and GSE105127 acquisition chains and activate only after
   their own read-only readiness checks pass.
9. Once specialist rectangles and phenotype-linked audits are complete, run the
   prespecified complementarity audit. Activate a context-conditioned model only
   if residual correlation and cross-fitted stacking checks pass.
10. Freeze finalists and SelectionLocks before one-time held-back evaluation. A
    failed hold back remains a failure; do not repair or reselect on the same holdout.

## Useful takeover commands

```bash
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

squeue -u "$USER" -o '%.18i %.26j %.10T %.10P %.10q %.12M %.12l %R'

find Analysis/MASLD_Model_Benchmark/executions/model-cpu-train-503-21097161.staging/predictions/folds \
  -name fold_receipt.json | wc -l

sacct -j 21100294,21097161,21100574 \
  --format=JobID,JobName%28,State,Elapsed,MaxRSS,ExitCode -P

git status --short -- Analysis/MASLD_Model_Benchmark

sha256sum \
  Analysis/MASLD_Model_Benchmark/MODEL_RESULTS_SUMMARY.md \
  Analysis/MASLD_Model_Benchmark/AGENT_HANDOFF_2026-08-25.md
```

## Metric and output file authorities

- Cell-state benchmark: `RESULTS_2026-08-23_cell_state_mapping.md`
- Transferability audit: `TRANSFERABILITY_2026-08-23.md`
- Contamination audit: `CONTAMINATION_AUDIT_2026-08-23.md`
- TF 50k evaluator:
  `executions/model-cpu-score-504-21083585/evaluation/head_metrics.tsv`
- DNA-language models:
  `executions/model-cpu-train-605-21099008/evaluation/model_head_summary.tsv`
- Enformer/Sei:
  `executions/model-cpu-train-604-21097189/evaluation/family_native_summary.tsv`
- LS-GKM:
  `executions/model-eval-626-21100570/evaluation/ensemble_metrics.tsv`
- Observed-multiome complete rectangle:
  `executions/model-cpu-score-328-21100332/receipt.json`
- Histology screen:
  `executions/model-scoring-078-21092130/scores/standardized_metrics.tsv`
- GSE260666 transfer:
  `executions/model-eval-084-21096365`
- Sequence completion snapshot:
  `executions/sequence-regulatory-completion-snapshot-21101005`
- GSE296875 transfer fixture:
  `executions/model-data-332-21100842`
- Corrected TF replay:
  `executions/model-check-337-21101451`

## Final caution

The campaign has excellent breadth, but breadth is not the paper's central
scientific result. The next agent should keep the family-native tournament
running while directing new work toward the missing connection: donor phenotype
to lineage-specific expression and regulatory remodeling across independent
cohorts. A high score on broad cell identity, observed-ATAC reconstruction, or
MPRA alone cannot support the requested MASLD remodeling claim.

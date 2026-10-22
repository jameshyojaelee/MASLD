# MASLD model program started 2026-10-15: every model, its status and what is waiting

Written 2026-10-16 for the PI. Nothing in this document is adopted. Every model from this program is a
candidate, not adopted, until you decide otherwise.

Path shorthand used below:
- P = `Analysis/MASLD_Model_Benchmark/executions/model-program-20261015T223347Z`
- Every number is copied from the file named beside it. Where two files disagree, both are given.

## 1. All models on one screen

Class words: improves / comparable / worse / inconclusive / no call. "Exposure" means whether the test cohort could
have been seen before (by a pretrained model, or by earlier analyses in this project).

| Model | Plain question (inputs → output) | Data, n | Track | Status | Headline result (class; exposure) | Decision pending from you |
|---|---|---|---|---|---|---|
| **Model A-final**, in-biopsy allelic | Tag-SNV allele counts in each person's liver RNA → does the GTEx cis-eQTL effect change in proportion with fibrosis stage (κ), and does that change vary by gene (τ)? | 558 staged development participants, 5 cohorts; sealed replication set of 237 | A | Running. Gate A passed. Gate B (simulation check of the test sizes) not yet started at full scale | No result yet (no call). Gate A: κ detectable down to 0.0265 at 80% power, 392 genes; P2 (persistence) not testable | Gate B cells beyond the judged tier; opening the sealed set; any override of the wave-1 lead rulings before the first Gate B dataset |
| **B1**, pairwise conformal severity calls | Bulk RNA severity score for two people → call which one has more fibrosis, or abstain | Confirmatory GSE213621, 257 participants (21,986 discordant pairs); development pool 449 | B | Result in; candidate package built | Comparable to released v1.2 under the 0.02 margin (coverage difference −0.0106, 95% CI −0.0197 to −0.0021, entirely below 0). No pretrained model | Whether to put the package in `release/` and add a MODELS.md row |
| **B2**, Geneformer V2-104M on bulk RNA | Bulk rank-value expression → fibrosis order; pretrained vs random-weight twins and vs linear ridge | 454 genetic individuals, 4 cohorts; transfer GSE281797 n = 94 | B | Result in, reviewed | Pretrained vs twins: no call (Δ −0.08). Pretrained vs ridge: worse (−0.44). Exposure unknown in every cohort, so descriptive | MODELS.md candidate row; whether to keep the lead's amended masked-token check (if reversed, B2 closes with no Stage 2 number) |
| **B3**, foundation-model cell embeddings | Donor-mean single-cell embeddings → cross-sectional stage 0/1/2, vs a linear pseudobulk baseline | T1: 55 evaluation donors (GSE202379, GSE244832); T2: GSE289173, 46 donors (S computed on the 26 FL-series donors) | B | Result in, reviewed | T1: inconclusive for TranscriptFormer and scGPT (exposure clean_declared). T2: descriptive only, second use of GSE289173 labels | MODELS.md candidate row |
| **B4**, RNA → H3K27ac transport | Frozen rank-4 RNA-to-H3K27ac map vs PCA ridge → per-donor H3K27ac profile in a new cohort | Trained on GSE267145 (99); tested on GSE128072 (17 donors) | B | Result in, reviewed | Worse (Δ −0.068, 95% CI −0.113 to −0.014). Test donors were seen before: development evidence | MODELS.md candidate row (draft exists) |
| **B5**, missing-component program prediction | 714 CosMx panel genes → 8 held-out program scores; reduced-rank vs full ridge | GSE268273, 109 participants | B | Result in, reviewed | Comparable but degenerate: reduced-rank picked full rank, so it is the same model as ridge (Δ 0.000). Spent development cohort | MODELS.md candidate row (draft exists) |
| **C1-01**, within-person allelic RNA at credible-set variants | Allele counts at transcribed fine-mapped liver-trait variants in heterozygous people's MASLD biopsy RNA → does the variant-carrying haplotype change its gene's RNA by more than 0.2 log2? Pooled over stage; stage dependence is not estimated | 936 development units (660 discovery, 276 holdout) | C | Stopped by kill line 8 | No outcome computed (no call). Projected 164.7 GPU-h against a 110 GPU-h kill line | Whether to raise its GPU limit to about 180 GPU-h and restart |
| **Track C ideas not built** | 58 other ideas from 40 researchers | (various) | C | Idea not built | Wave 1: 40 ideas, 1 viable. Wave 2: 19 ideas, 0 viable. Six wave-2 ideas sit just below the bar | Whether to build any of the six near-bar ideas |
| **Model D**, within-stage fibrosis decomposition | Fallback if Model A failed Gate A | n/a | B (closure list) | Dormant | Not opened, because Model A passed Gate A | None |

Sources: `P/decisions_pending.md`; `P/program_state.json`; `P/trackA/A4b_gateA/GATE_A_DECISION.md`;
`P/trackB/B1/RESULT.md`; `P/trackB/B2/RESULT.md`; `P/trackB/B3/RESULT.md`; `P/trackB/MODELS_ROWS_DRAFT.md`;
`P/trackC/ideas/C1-01/state.json` (`kill8_record_2026-10-16`); `P/WAVE_LOG.md` line 199 (C1-01 split counts);
`P/trackB/CLOSURES.md` §6 (Model D); `P/trackC/ideas/C1-01/PRESPEC_FROZEN.md` §1 (C1-01 question).

## 2. Track by track

### 2.1 Model A-final (Track A)

**What it is.** A beta-binomial model of each participant's tag-SNV allele counts in bulk liver RNA. Genotypes come
from a prior, not from DNA. κ is the proportional change in the GTEx cis effect per fibrosis step (F0-1, F2, F3-4).
τ is how much that change varies between genes. Three tests, Holm-corrected as one family: P1 (κ = 0), P2
(persistence, |κ/λ| < 0.125 plus a bound on τ), P3 (τ = 0). Data: 558 development participants in five cohorts
(208, 160, 83, 54, 53) and a sealed replication set of 237. Source: `P/trackA/CONTEXT_modelA.md`.

**Gate A: passed (2026-10-16).** Source: `P/trackA/A4b_gateA/GATE_A_DECISION.md`.
- Smallest detectable κ at 80% power: 0.0265 (0.0306 at α 0.05/3), at θ = 0, information I = 11,158.02.
- 392 κ-eligible genes.
- Branch: P1 + P3 + per-gene fields. P2 stays in the Holm family with p = 1, "not testable", because the GSE193066
  label-reliability arm was undecidable.
- The primary and blind recomputations agree to below 1e-17 on every decision quantity.
- Model D does not become the lead.

**Settings completed before Gate B (none uses the single allowed amendment).**
- v1.8 and addendum v1.8a: settings the spec left open. Sources: `P/decisions_pending.md` ("Model A spec v1.8",
  "Model A addendum v1.8a").
- v1.8.1: two exact speed-ups adopted (bitwise equal, 120/120). A looser section-3 tolerance was rejected because it
  raised the fit failure rate (fixture 0 to 2 of 121). Source: `P/trackA/LEAD_RULINGS_W4_2026-10-16.md` H-9.
- v1.8.2: you adopted the "asc" convergence rule. Under the old rule, observed section-3 fits failed in 5 of 80 and
  3 of 80 simulated null datasets. Source: same file, H-10.
- v1.8.3: the observed fit uses two starts and reports the converged one with the higher likelihood. Estimated cost
  +812 CPU-h for the judged tier [Estimate]. Review verdict: pass, no major items. On the 60-person fixture, P1 is
  identical to v1.8.2 (p 0.17647, 50 of 50 draws equal, 0 failures). Sources: H-11 in the W4 rulings;
  `P/trackA/A6_gateB/v183-REVIEW.md`.
- v1.8.3a (pending): the lead ruled three review points (H-12). An unconverged start with higher likelihood makes
  the fit count as failed. The second-start seeds are fixed. It is applied after the GPU version returns, so the CPU
  and GPU copies get the same change. Source: W4 rulings H-12.

**Gate B plan.** Gate B simulates data where the answer is known and checks that each test rejects at the right rate.
- Pass rule: size in [0.04, 0.06] with at least 4,000 datasets per null cell; 95% κ interval coverage in
  [0.93, 0.97]; 200 full-CRT replicates. Source: `P/trackA/CONTEXT_modelA.md` "Gate rules".
- Your decision (2026-10-16): run the judged tier first (null cells 1, 2, 3, 40, 27 and the full CRT). Cap about
  100,000 CPU-h, pause at 80,000. Remaining cells are decided after real P3 costs are measured.
  Source: `P/budget.json` `user_decision_2026-10-16_gateB`.
- The 600-dataset pilot was held and then cancelled. Smoke costs projected it at 2,700 CPU-h or more against a
  400 CPU-h estimate; one dataset in each of cells 1 to 3 was still running past about 27,960 CPU-s. It is to be
  rerun on the new code. Source: `P/trackA/A6_gateB/state.json` `pilot_hold_2026-10-16`.

**Gate B cost.**
- Cost profile before the speed-ups (`P/trackA/A6_gateB/cost-profile-20261016T104945Z/options/options_cost.json`):
  - Baseline: full Gate B 273,587 to 524,608 CPU-h; judged tier 179,545 to 304,860 CPU-h.
  - Cheapest option at the required 4,000 datasets per cell (labelled `c_exact_tol34`): full 82,897 to
    160,723 CPU-h; judged tier 54,369 to 93,226 CPU-h. Options with fewer datasets cost less (down to 59,247 CPU-h
    full at 1,000). `c_exact_tol34` adds looser optimizer stopping in sections 3 and 4 to the exact speed-ups
    [Inference: from the option name and `scripts/prototype_lbfgs.py`]; the section-3 tolerance change was later
    rejected (H-9), so this cost is not available as costed.
- Optimizer test, 20 datasets (`P/trackA/A6_gateB/newton-spike-20261016T165340Z/collected-20261016T171020Z/summary.json`):
  median section-3 fit 264.6 CPU-s with the current optimizer; 147.7 CPU-s with the trust-region method ("tc"),
  which matched the current optimum within 1e-8 in 19 of 20; 19.5 CPU-s with the "pn" method, which reached a lower
  optimum in 10 of 20. Throwaway test only; nothing adopted.

**GPU budget.** Your decision (2026-10-16): Track A GPU cap 2,500 GPU-h, pause at 2,000, up to 14 GPUs. Plan: a
float64 GPU version of v1.8.3, checked against the CPU version, with the judged tier run on mixed GPU and CPU.
Source: `P/budget.json` `user_decision_2026-10-16_gateB_gpu`. The GPU version checks now running are estimated at
15 to 25 GPU-h [Estimate] (`P/trackA/A6_gateB/gpu-backend-20261016T170127Z/state.json`). Their results were not in
the files read for this document.

### 2.2 Track B (B1 to B5)

**B1. Pairwise conformal severity calls.** Comparable to the released v1.2 scorer. On GSE213621 (257 participants),
the 0.90 rule covers the true order in 0.9016 of discordant pairs (95% interval 0.8736 to 0.9278) and makes a call on
0.7863 of them, with call accuracy 0.8749. The interval reaches below the 0.875 bar, so the pass does not show
coverage of at least 0.875 with confidence. Against v1.2 at the same call rate, the coverage difference is −0.0106
(95% CI −0.0197 to −0.0021): comparable under the 0.02 margin, with the interval entirely below 0. The refit score adds nothing over v1.2; B1 adds a calibrated rule for when to abstain.
A candidate package is built and its test reproduces the numbers. It is not in `release/` and has no MODELS.md row.
Sources: `P/trackB/B1/RESULT.md`; `P/trackB/B1/state.json`.

**B2. Geneformer V2-104M on bulk RNA.** 454 genetic individuals in 4 cohorts. Pretrained vs random-weight twins:
Δ −0.08, no call, because shifting the ridge penalty one step moves Δ by more than 0.01. Pretrained vs linear ridge:
worse, −0.44 (97.5% CI −0.59 to −0.17). Transfer to GSE281797 (n = 94): Δ +0.18, no call; pretrained minus ridge
−0.15, inconclusive. Every number is for bulk rank-value input to a single-cell-trained encoder. Exposure is unknown
for every cohort, so every number is descriptive. Two independent recomputes reproduced every call and label. The
embedding-arm values move by up to 0.0075 between float32 and float64 scoring, so they are given to 2 decimals. Sources: `P/trackB/B2/RESULT.md`; `P/decisions_pending.md` ("Track B, B2 complete").

**B3. Foundation-model cell embeddings vs a linear pseudobulk baseline.** T1 used 55 evaluation donors. Stage score S:
TranscriptFormer 0.0859, scGPT 0.0791, baseline 0.5367. Every delta for both arms is inconclusive, so there is no
headline. Inconclusive here is not "no difference": the delta against the baseline is about −0.45, and its
conditional interval lies below 0 while the refit interval crosses 0. Exposure is clean_declared on both evaluation cohorts. UCE was dropped before the readout because its
decoder failed a check on this liver data; that is not a result on the MASLD task. T2 on GSE289173 (46 donors) is
descriptive only and is the second use of those labels. S is computed on the 26 FL-series donors: baseline 0.60,
TranscriptFormer 0.17 (exposure unknown), scGPT 0.02 (descriptive class worse against the baseline), Geneformer 0.13. Sources: `P/trackB/B3/RESULT.md`; `P/decisions_pending.md`
("Track B, B3 complete"); `P/WAVE_LOG.md` line 200 (46 donors).

**B4. RNA to H3K27ac transport.** Rank-4 map vs PCA ridge, trained on GSE267145 (99), tested on 17 GSE128072 donors.
Worse: Δ −0.068 (95% CI −0.113 to −0.014), sign-flip p 0.020. Both maps scored below the training reference in 17 of
17 donors. Dropping chrX and chrY regions gives Δ −0.054 (−0.101 to +0.004), which includes 0. The test donors had
been seen by earlier jobs, so this is development evidence. Sources: `P/trackB/B4/RESULT.md`;
`P/trackB/MODELS_ROWS_DRAFT.md`.

**B5. Missing-component program prediction (GSE268273, 109 participants).** Comparable, but degenerate: the
reduced-rank model chose full rank in 5 of 5 folds, so it is ridge. Skill 0.777 (95% 0.705 to 0.836); Δ 0.000.
Technical covariates alone give 0.564; panel genes add 0.216. Not shown to be program-specific: p 0.178 against
matched gene sets. The acceptable-error line was untestable (power 0.764 < 0.8). Sources: `P/trackB/B5/RESULT.md`;
`P/trackB/MODELS_ROWS_DRAFT.md`.

**Closed in Track B without a new run** (`P/trackB/CLOSURES.md`): caQTL target-conditioned head (+3.15% against a
5% margin); NAS activity head (needs an unused cohort with NAS and n ≥ 163 [Estimate]); within-histology
RNA-to-H3K27ac contrast (22 pairs, almost all F0 female); GSE223843 RNA-to-ATAC transport (n = 6, both source models
null); GSE244832 RNA-to-ATAC (no person-level link between snRNA and snATAC); Model D (dormant).

### 2.3 Track C (new ideas)

- Two research waves. Wave 1: 24 researchers, 40 distinct ideas, 10 passed triage, 1 viable (C1-01). Wave 2: 16
  researchers, 19 ideas, 12 passed triage, 0 viable. The prespecified stop rule (two waves in a row with at most one
  new viable idea) was met; 40 of 120 researchers were used. Source: `P/program_state.json` `trackC.research_loop`.
- Why the other 58 were not built (counted from the `reason` field of `P/trackC/idea_registry.jsonl`): 37 failed
  triage, most often because the smallest detectable effect at the available n was larger than the planned margin
  or there were too few donors; a few on controlled-access labels (C1-29, C1-34, C1-35) or because the analysis
  already exists (C1-21, run as B4; C1-28; C1-32). 18 passed triage but scored below the prespecified bar on the
  judges' medians (6 in wave 1, 12 in wave 2). 3 passed triage and were vetoed (C1-05 and C1-40 on statistics,
  C1-24 on novelty). Per-idea reasons are in
  `P/trackC/wave1/WAVE1_SUMMARY.md` and `P/trackC/wave2/WAVE2_SUMMARY.md`.
- Near the bar in wave 2 (not built because the bar was prespecified): C2-05 image fat and fibrosis vs programs in
  GTEx; C2-08 consequence-partitioned COLOC Genetics field; C2-10 MAGMA and burden on disease-state genes; C2-16
  protein severity v1.3 external test; C2-17 zonal protein loss; C2-18 MESC. Source: `P/decisions_pending.md`.
- **C1-01 stopped by kill line 8.** The simulation plan projected 164.7 GPU-h after all four prespecified fallback
  steps, above the 110 GPU-h kill line and the 150 GPU-h per-model cap. Spent 182.0 CPU-h and 0.12 GPU-h at the stop
  (`kill8_record_2026-10-16`). The same state file's running totals are 215.92 CPU-h and 0.85 GPU-h. The GPU total
  equals 0.73 GPU-h before the operator stages plus 0.12 at the stop. The CPU fields do not add up the same way
  (37.05 before the operator + 3.2 other jobs + 182.0 = 222.25, not 215.92), and the file does not explain the gap. No stage-2 outcome was computed. Source: `P/trackC/ideas/C1-01/state.json`.

### 2.4 Models from before this program

Source for this table: `P/CONTEXT_models.md`.

| Model | What it does | State | Number |
|---|---|---|---|
| masld-severity-v1.2 | Orders bulk RNA by fibrosis and NAS activity | Released | |
| masld-panel-v1 | Counts or Cq panel | Released; qPCR use not established | |
| masld-liver-protein-severity-v1.3 | Orders liver proteomes by fibrosis | Released | |
| masld-liver-chromatin-state-v1.2 | RNA → H3K27ac state | Released | In-cohort skill +0.1413 [+0.1125, +0.1686]; transfer to 39 donors ρ +0.256 vs null +0.071 |
| Model B, noise-aware aspects | | Closed | K = 2 preferred (p 0.005); sampler did not converge |
| Model C, H&E → programs | | Done | 16 of 79 programs readable; TCGA endpoint not met (1 of 16) |
| Model F, mouse diet vs human stage | | Done; mouse, outside this paper | No specificity |
| Ordinal cumulative-logit on RNA | | Rejected by kill line | −0.0076 |
| Anchor regression | | Rejected | Δ 0.0000 |
| Pairwise L2 | | Rejected | Equal to ridge |
| RankNet | | Rejected | −0.0131 |
| Absolute conformal | | Rejected | Coverage failure |
| DEBIAS-M | | Rejected | No change |
| Teacher-student | | Rejected | −0.057 |
| MOFA+ | | Rejected | −1.74 |
| BulkFormer | | Rejected | Loses to its random-init twin |
| Corgi zero-shot / fine-tuned | | Rejected | +0.00028 / +0.00016 |
| GPU MLP | | Rejected | −0.0014 |
| Plasma → liver protein | | Stopped | |
| GSE296875 context increment | | Rejected | Collapses under leave-one-well-out |

Names only, with no prespec or code, listed in `P/CONTEXT_models.md` as seeds for new ideas: A-T7 GTEx transport;
PRESPEC_O cross-layer coupling; conditional SAF remodeling (GSE202379, PXD051911); zonation distortion (GSE105127);
haplotype or context at a few loci.

Power limit for any program-level test: smallest detectable difference 0.0667 at K = 14 on 585 participants;
reaching 0.05 needs about 1,041 participants (`P/CONTEXT_models.md`).

## 3. Decisions waiting for you

1. **B1:** put the candidate package into `release/` and add a MODELS.md row, or leave it as a candidate, not adopted.
2. **B2, B3, B4, B5:** add MODELS.md "Candidates (not adopted)" rows. Drafts for B4 and B5 are in
   `P/trackB/MODELS_ROWS_DRAFT.md`; B2 and B3 have none yet.
3. **B2 masked-token check:** keep the lead's amendment, or reverse it, which closes B2 with no Stage 2 number.
4. **C1-01:** raise its GPU limit to about 180 GPU-h and restart, or leave it stopped.
5. **Track C near-bar ideas:** build any of C2-05, C2-08, C2-10, C2-16, C2-17, C2-18, or none.
6. **Model A Gate B, after the judged tier:** which remaining cells to run, once real P3 costs are measured
   (`P/budget.json`).
7. **Model A lead rulings:** you may override any wave-1 lead ruling (D-SPEC, D-QC, D-LABEL items) before the first
   Gate B dataset is generated (`P/decisions_pending.md`). The wave-2 to wave-4 rulings (V, Q, H items) are listed
   there and in `P/trackA/LEAD_RULINGS_W2..W4`, but none of those files states the same override window for them.
8. **Model A sealed set:** opening the 237-person sealed replication set needs your explicit authorization
   (`P/trackA/CONTEXT_modelA.md`).

Decisions 1-5 above were answered "yes to all" on 2026-10-17 (`P/decisions_pending.md`). They were then overtaken by the
deep pass below.

## 4. Deep pass, 2026-10-17: what each model could change, its ceiling, raw-data checks, cut

D = `P/deep-20261017T022651Z`. Per-model assessment and critic: `D/models/<model>/ASSESSMENT.md`, `CRITIC.md`.
Label rebuild from GEO and source supplements: `D/labels/<cohort>/`. Reliability table: `D/reliability/`.

| Model | Claim it could change | Ceiling and share | Raw-data check | Outside evidence | Call |
|---|---|---|---|---|---|
| Model A | New allele-by-stage cis result (PAPER.md:599-603); Catalog field `cis_allelic_effect_in_MASLD_liver` | Detectable true κ 0.042-0.053 at point λ (0.092 at the λ 2.5th percentile), below Δ 0.125. Per-gene κ: median MDE 0.83, 0 of 392 genes reach Δ | 83 of 2,343 kept tag SNVs never piled up (80 genes; 28 with no rows). Depth CCC 0.9992 at piled sites. REF/ALT: R3 passes; R5 fails (WASP raises pooled ALT share 0.0085). Labels match source. Beta-binomial GLMM agrees (κ̂ r 0.88-0.96); mixQTL fails the α field | Sealed set (needs user); spec §9 controls; G2 and PNPLA3 at Gate C | Keep: P1, P3, α field |
| C2-05 | New field `gtex_postmortem_fibrosis_state` | Note arm 0.42-0.51 of ceiling; image arm 0.54-0.71 [Estimate]; stage-2 power 0.21 with note error | Labels match source; scorer stopped at its near-duplicate rule, no outcome | TCGA-LIHC holdout (MDE 0.48); outside pathologist read | Keep, after fixes |
| C2-10 | One ALT supplementary sentence | T1 MDE 0.119 vs margin 0.15 | Inputs reproduce from Pan-UKB | BBJ holdout (power 0.26) | Stopped by its kill lines 1a and 3b; no outcome |
| B1 | None | v1.2 0.948 of ceiling, B1 0.937 | No label disagreement | None unused | Stop; release stays a candidate |
| B2 | None | Ridge at 0.80 of a 0.80 ceiling | Composition-only ridge 0.40 beats every Geneformer arm | None | Stop |
| B3 | Candidate Figure 7 row only | Baseline 72%, TF 12%, scGPT 11%; room 0.21 < MDE 0.50 | 3 donor labels differ from source; scGPT +0.02 exclusion fails under source labels | None suitable | Stop |
| B4 | None | Read-sampling ceiling 0.94 | Composition removes 64% / 52% of development skill | Sex-chromosome check, not run | Stop |
| B5 | None | Δ MDE 0.119 > margin 0.05 | Numbers reproduce | Spent cohort | Stop |
| C1-01 | Weak field beside H87 | Kill line 2 unreachable | n/a | GTEx eGene panel | Stopped |
| C2-08 | Fig 2C class count | P1 power <= 0.054 | n/a | n/a | Stop |
| C2-16, C2-17, C2-18 | None or one conditional sentence | C2-17 MDE 0.78 > 0.70; C2-18 power < 7% | Labels match source | n/a | Stop |

User decisions on 2026-10-17 (second "yes to all", `P/decisions_pending.md`): re-pile the never-piled tag SNVs and
rebuild Gate A before Gate B; drop the per-gene κ fields; report the beta-binomial GLMM beside Model A; release the
REF/ALT recount aggregates (`D/raw/modelA_recount/20261017T035409Z/RELEASE_NOTE.md`); C2-05 fixes, amendment and
outside read; MODELS.md B3 label note. Track A ruling: `P/trackA/LEAD_RULINGS_W4_2026-10-16.md` H-16.

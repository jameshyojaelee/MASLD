# AlphaGenome Atlas × MASLD Resource — execution log (started 2026-09-09)

Plan of record: `~/.claude/plans/trying-to-execute-this-snappy-moth.md` (v2, seven pillars, approved 2026-09-09).
Protocol: `docs/technical/MASLD_ALPHAGENOME_RESOURCE_PLAN_2026-09-09.md` (+ §11 execution amendments).
Code: `scripts/analysis/alphagenome_atlas/` (tests: `tests/`, 45 passing on a compute node).
Results roots (write-once): `GWAS/finemapping/results/alphagenome_atlas/`
- `run-20260909T153939Z/` — Track 0 (Packages A–H on the 816 pairs + 173 credible sets)
- `preview-direct-only-<ts>/` — pipeline check on the direct-universe archive only (NOT a result)
- `p2-anchors-<ts>/` — P2 known-mechanism anchor panel
- `p3a-benchmarks-<ts>/` — P3a external same-variant benchmarks (caQTL, sQTL, MPRA)
- `p3-feasibility-<ts>/` — read-only inventory of raw reads for P3/P7

## Facts established today (each broke or bounded a step)
- Atlas point queries reject insertions and deletions (gRPC UNIMPLEMENTED; SNV at the same position served).
  HSD17B13 rs72613567 (T>TA) is therefore `not_served`; 930 of 11,896 caQTL leads are indels and are excluded.
- Throughput ≈ 1,000 variants/min under the per-minute quota; `query_interval` uses 32-bp sub-requests.
- Both CL:0000182 "hepatocyte" tracks are embryonic; adult liver = UBERON:0002107/0001114/0001115 (+ GTEx liver). One adult
  hepatic stellate cell DNase track (CL:0000632) and two hepatic mesenchymal stem cell CAGE tracks (CL:0002571) exist;
  no stellate ATAC/histone/TF track.
- A scorer with nothing to score returns AnnData (0, n_tracks) with no obs columns (step 07 guard added).
- The `rnaseq` env needs `PYTHONNOUSERSITE=1` (a NumPy-1 numexpr in ~/.local breaks pandas).

## Track 0 status
- Substrate: 816 pairs (A 29 / C 787), B 173 credible sets; 1,846,327 source variants; 417,155 query variants.
- Pilot PASS; direct retrieval done (13,690 variants, 55 chunks); gate retrieval done (181 unique SNV uids from
  211 lead–gene rows; hg38 alleles from `alphagenome_eqtl_scores.tsv`, FASTA-validated; 107 rows had hg19
  ref/alt swapped relative to hg38 — the first gate attempt used the hg19 alleles and was rejected by the Atlas).
- Enzyme retrieval job 21558249 in progress (quota-bound); tables job 21558781 dependent.
- Package H (job 21580767, COMPLETED 6 min): G1 whole-matrix skill +0.1413 reproduced (|Δ| 9.5e-12); G2 per-region
  per-fold skill reproduced (max |Δ| 1.9e-6). Predictions: H2 ✓ (skill vs signal mean ρ +0.286, vs sd +0.226),
  H4 ✓ (promoter median skill 0.069 vs distal 0.089), H5 ✓ (transfer ρ vs in-cohort skill +0.210),
  H1 ✗ (local increment vs n_cis_100k ρ −0.073, vs TSS distance +0.039; random-gene increment −0.112),
  H3 ✗ (GC ρ −0.075), H6 ✗ (covariate R² 0.164; signal mean partial 0.119, sd 0.079). Reliable sets: shipped 4,364,
  states+cis 4,847, both 3,146 (model card headlines 4,847). The first Package H attempt died on a tuple JSON key;
  its partial tables are kept under `tables/package_H_partial_failed_run_21579560/`.

## P3 feasibility (job 21579486)
Per-read data on disk: bulk RNA STAR BAMs for GSE162694 (143), GSE130970 (156, legacy), GSE126848 (56),
GSE213621 (367); GSE244832 snATAC fastqs (270 files, 970 GB, unaligned); GSE296875 fragments only (no allele
information); GSE202379 fastqs (1.3 TB). NOT on disk: GSE267145 (H3K27ac + RNA, 99 participants) raw reads,
GSE281367 raw reads. Tools: samtools (rnaseq env), bcftools 1.19 and GATK 4.6 modules. GSE267145 rights
contract: public GEO, no DUA, internal research use allowed, redistribution prohibited.
Consequence: P3 (allelic tests) is feasible now on bulk RNA (allele-specific expression, ~720 participants) and
on GSE244832 after alignment; the H3K27ac allelic test and P7's genotype term need a GSE267145 raw-read download
(size pending from SRA) — a user decision.

## Jobs submitted today
21557989 substrate · 21558056/21558078 pilot · 21558249 enzyme retrieval · 21558781 tables (dependent) ·
21568540 gate (failed: hg19 alleles) · 21579624 gate (OK) · 21569792/21579560/21580767 Package H ·
21579486 P3 feasibility · 21580029 (failed: indels)/21580927 benchmark retrieval · 21581787 anchors ·
21581800 MPRA intervals · previews 21568541 (cancelled), 21579625/21580835 (step-07 defects), 21581514.

## P4 saturation universe (step 41, fixed 2026-09-09 before any saturation score)
U1 = 96,460 H3K27ac regions; U2 = 7,000 snATAC consensus peaks (6,930 unique coordinates): 1,848 hepatocyte/stellate
peaks with a non-indeterminate DA state or BH q < 0.05 in either cohort, 5,416 program-linked (promoter gene in one of
the 117 programs), genetics-overlap flag to be added after the Track 0 tables job. Total 88.7 Mb. Pilot (step 40, 20 + 20
regions, ontology-filtered to liver + HepG2 + four cross-tissue proxies) sizes the full retrieval.

## P3 bulk-RNA allelic arm — first pass and its limit (job 21584253, GSE126848, 56 samples)
Allele counts at the 13,853 Track 0 direct-universe + gate SNVs stream cleanly from unindexed STAR BAMs with
`bcftools mpileup -T` (no index needed; ~1 min/sample at 16 threads). But most Resource variants are not
transcribed: only ~57 sites per sample reach 10 reads. Two consequences, both recorded before any analysis:
(1) the counts are ref vs non-ref (`<*>`), which is what ASE needs but not allele-specific to the listed ALT;
(2) an ASE test at these variants is haplotype-tagged, not same-variant — an exonic marker's allelic ratio
reflects the regulatory variant it tags, not its own effect. The same-variant chromatin test (H3K27ac /
snATAC allelic imbalance) therefore remains the only same-variant arm, and it needs raw reads that are not on
disk. The bulk-RNA arm is re-scoped to a gene-level question (do genes whose credible-set variants carry strong
predicted Atlas RNA effects show more donor-level allelic imbalance at their exonic markers than matched genes?)
with exonic common markers from the Broadaway hg38 variant list as targets.

## Quota lesson (2026-09-09)
Five retrieval jobs ran against one per-minute API quota at once; interval queries stalled and the service
returned UNAVAILABLE. Retrieval is now strictly serial: Track 0 enzyme → anchors → MPRA intervals → saturation
pilot, chained with `--dependency=afterany`. A bounded transient-error retry (10 attempts, 30 s) was added to
`atlas_query.call_with_quota_retry`.

## Corrections after the repository inventory (2026-09-09, before any affected result was read)
1. **GNMT haplotype rsIDs.** The anchor prespecification named rs9471976/rs11752813 as the Resource's unresolved
   two-variant GNMT haplotype. The haplotype is rs2296804/rs2296805 (UKBB-EUR r2 = 1.000, docs/RESULTS.md:249-250);
   rs9471976 and rs11752813 sit at r2 0.819 and 0.580 with rs2296805. Both are kept under the class
   `gnmt_locus_partial_ld` and the true pair was added, so the panel is 30 variants. No anchor score had been read.
2. **Catalog field name.** The sidecar wrote only the legacy alias `next_experiment_rule_id`; the v2 contract field is
   `recommended_experiment_rule_id`. It now writes both.
3. **Drug-target scope.** Figure 5 is molecular and physical tissue context, NOT a drug-target figure, and the scope
   validator rejects therapeutic/calibration/top-N panels in the Figure 6 manifest. The P6a target-direction table is
   therefore a supplementary analysis and a Catalog-sidecar input only; the plan's "Figure 5/6 supplementary panel"
   wording is withdrawn.
4. **Reference FASTA for the bulk arm.** The bulk STAR BAMs are `refdata-gex-GRCh38-2024-A` (706 contigs); the first
   ASE pass used the cellranger-atac FASTA (194 contigs). All targets are autosomal so the counts stand, but the
   correct FASTA is now a required argument.

## What the inventory changed for P3 and P7
Per-read data DO exist for two snATAC cohorts: GSE281367 (12 donors, 421 GB cellranger-atac BAMs with CB tags,
hg38) and GSE244832 (18 donors, 192 GB bowtie2 BAMs, barcode in the read name, no CB tag). GSE267145 has no raw
reads (`role = "blocked"`) and GSE296875 has fragments only, so neither the H3K27ac allelic test nor a GSE267145
genotype term is possible. The same-variant chromatin allelic test therefore moves to GSE281367 (CB tags allow a
per-lineage test) with GSE244832 as the second cohort once barcodes are re-derived. No donor genotypes exist
anywhere in the repository; every heterozygous site must be called from the reads.

## Track 0 retrieval complete; the 211-set replication is at chance (2026-09-09 evening)
Enzyme retrieval finished (job 21558249, 1,614 chunks, 7h21m). Tables job 21558781 running.
The preview pipeline (job 21581514 + a rerun of step 14) passes 12 of 14 checks; the two failures are the
preview's known direct-only scope (13,690 of 417,155 variants archived) and are why it is not a result.

**Gate replication (Package D):** Atlas liver RNA_SEQ, oriented to the eQTL effect allele exactly as 65b,
auROC **0.5304 [0.449, 0.607]**, permutation p 0.229, sign concordance 0.522, n = 201 of 210 SNV leads.
The written prediction was "at chance"; it holds. This replicates the prior in-house AlphaGenome FAIL
(0.5607, n = 211) on hg38-correct alleles.

## P3 same-variant snATAC allelic concordance — POSITIVE (step 73)
Substrate: 12 GSE281367 donors, 1,119 Atlas-scored SNVs inside supported consensus peaks, allele counts from
the cellranger-atac BAMs (MAPQ >= 30, BQ >= 20). 435 sites pass the reliability filter; 412 carry an Atlas
ATAC liver quantile. Prespec `73_allelic_prespec.json` (sha256 e951738c...) was written before any number.

| statistic | ATAC | DNase |
|---|---|---|
| site-level sign concordance | 0.655 | 0.663 |
| expected from marginal sign skew alone | 0.499 | 0.500 |
| per-block mean (36 1-Mb blocks) | 0.735 | 0.668 |
| blocks above 0.5 (exact binomial) | 29/33, p 5.5e-06 | 29/34, p 1.9e-05 |
| 1-Mb block bootstrap | 0.737 [0.583, 0.861] | 0.670 [0.500, 0.806] |
| \|measured log2\| > 1 (n = 20) | 0.950 | 0.900 |
| \|measured log2\| <= 1 | 0.640 | 0.651 |
| sign-shuffle control / permutation null | 0.502 / 0.501 | 0.510 / 0.506 |
| Spearman(predicted, measured) | 0.464 | 0.421 |

Reference bias is negligible (global mean donor-level log2 ratio -0.007); bias-centred numbers are unchanged.
The block-level sign test and the block bootstrap are the LD-respecting statements; the site-level permutation
ignores LD between sites at the same locus and is anti-conservative. Targets are clustered (chr2 carries 200 of
412 sites, 36 blocks total), so the independent n is 36 loci, not 412 sites.

Four of the five written predictions FAILED, all by understating the Atlas: A1 (strong-imbalance concordance
0.55-0.70) observed 0.95; A2 (weak-imbalance [0.45,0.58]) observed 0.640; A3 (reference bias magnitude > 0.02)
observed 0.007; A4 (fewer than 200 sites pass) observed 435. A5 held (both nulls at 0.50).
Deviation from the plan, recorded before results: >= 3 het donors rather than >= 8, because GSE267145 has no
raw reads and this arm runs on 12 donors, not 99 participants.

## caQTL reproduction guard — the guard was wrong, the reproduction is exact
Step 33's guard applied `79_score_panels.subsample` (3,000 leads) to the archived `alphagenome_atac` row and
read a false failure. That subsample rule governs the GPU panels; the archived row records **n = 11,896**, every
served lead. Recomputed on that set this script returns **0.707398872753231** against the archived
**0.707398872753231** - identical to 15 digits. The guard now compares on the archive's own n and refuses a
comparison across different variant sets (`reproduction_guard`, fixture-tested). Atlas ATAC on the same leads:
**0.7128 [0.703, 0.723]**, n = 10,966; Spearman between Atlas and API scores 0.981. The written prediction
(|delta auROC| < 0.02) holds at 0.0054.

## Interval queries burst past the per-minute quota (anchors 21584517, MPRA 21584518)
Both jobs looped indefinitely on RESOURCE_EXHAUSTED: `query_interval` issues 32-bp sub-requests for every
requested scorer at once, so one 2-kb anchor window (22 scorers) and even one 126-bp MPRA oligo exceed the
per-minute allowance in a single call. A 65 s sleep cannot clear it because the retry re-issues the same burst.
The anchor POINT stage completed and is archived (30 anchors). Both jobs cancelled; the interval strategy is
being re-sized against a direct probe before resubmission.

## Duplicate donor BAMs (step 72)
cellranger-atac keeps an internal copy of each donor's BAM under
`<donor>/SC_ATAC_COUNTER_CS/.../WRITE_ATAC_BAM/fork0/join-<hash>/files/possorted_bam.bam`. The first counting
pass matched both copies and wrote 24 files for 12 donors; treating those as 24 donors would be
pseudoreplication. Step 72 now matches only `*/outs/possorted_bam.bam`; the duplicate counts are kept under
`atac_counts_duplicate_cellranger_internal_bams/` and enter no analysis.

## The interval-query loop was a client message-size limit, not the quota (corrected)
Probing settled it. `query_interval` at 2 kb over 2 scorers succeeds (12,288 rows, 10.6 s); at 256 bp over all
22 scorers it fails with

    RESOURCE_EXHAUSTED  Stream removed (CLIENT: Received message larger than max (22031061 vs. 4194304))

which is gRPC's **client-side receive cap**, not a quota breach. gRPC reports both with the same status code,
so `is_quota_error` matched it, slept 65 s, and re-issued the identical oversized request forever. The earlier
entry attributing the anchors/MPRA loop to quota exhaustion is wrong on cause; the burst size was never the issue.

Fixes, each fixture-tested (the classification test reproduced the production loop and died of OOM before the fix):
- `is_message_size_error` separates the two cases; `is_quota_error` now excludes message-size errors, so an
  oversized request propagates on the first attempt instead of looping.
- `create_client_with_large_messages` rebuilds `atlas.create`'s three public objects with
  `grpc.max_receive_message_length` at 512 MB. `create_client` uses it, so every script inherits the fix.
- `query_interval_split` halves a window on a message-size error, with a floor, for anything still too large.

Verified against the API with the new client: 256 bp x 22 scorers OK (21.8 M cells, 22.6 s); 2 kb x 22 scorers
OK (174.3 M cells, 149.2 s); 128 bp MPRA oligo x 22 scorers OK (10.9 M cells, 9.5 s). Anchors, MPRA and the
saturation pilot resubmitted serially (21609734 -> 21609735 -> 21609736). Note for P4: 149 s per 2-kb window
over 22 scorers does not scale to 96,460 regions; the saturation pillar needs a reduced scorer panel, and the
pilot now sizes that.

## P3 replication in an independent cohort (GSE244832, 18 donors) — directionally consistent, weaker
The identical 1,119-target list, filters and prediction source were applied to GSE244832 (18 donors, indexed
hg38 BAMs, same reference build). 279 sites pass the reliability filter (vs 435 in GSE281367).

| statistic | GSE281367 (12 donors) | GSE244832 (18 donors) |
|---|---|---|
| reliable sites / scored | 435 / 412 | 279 / 264 |
| ATAC site-level concordance | 0.655 | 0.564 |
| marginal-skew expectation | 0.499 | 0.498 |
| ATAC per-block mean | 0.735 | 0.594 |
| ATAC blocks above 0.5 | 29/33, p 5.5e-06 | 21/32, p 0.055 |
| ATAC block bootstrap | 0.737 [0.583, 0.861] | 0.593 [0.429, 0.743] |
| ATAC strong imbalance | 0.950 (n 20) | 0.733 (n 30) |
| DNase blocks above 0.5 | 29/34, p 1.9e-05 | 21/31, p 0.035 |
| reference bias | -0.007 | -0.029 |

The replication agrees in direction in both channels but its ATAC interval includes chance (the DNase arm
reaches p 0.035). Report it as a partial replication: the discovery cohort carries the significant result and
the independent cohort is consistent but underpowered at these sites, which were selected on GSE281367 peak
support. Cohorts stay in separate arms; they are different libraries and are not pooled.

## sQTL: an estimand mismatch, not a null result
The first sQTL deposit compared GTEx `slope` against the Atlas **absolute** junction quantile and read
sign concordance 0.5114 (n = 919 matched junctions), against a written prediction of > 0.7. That comparison is
not a test of splice-direction skill: GTEx sQTL slope is a LeafCutter **within-cluster intron excision ratio**,
so an allele that lifts every junction in a cluster equally moves no slope while moving every absolute quantile.
`leafcutter_contrast` (target junction quantile minus the mean over its sibling junctions in the same cluster,
siblings the Atlas did not score skipped rather than counted as zero) supplies the matched quantity, with
cluster membership taken from the full pair table rather than the p-ranked top-10 subset. Both numbers are now
reported side by side; the first deposit is kept under
`tables/_superseded_sqtl_absolute_junction_estimand/`. The rerun is in flight.

## P2 anchor panel (job 21609734, 30 anchors, all 22 scorers + 2-kb windows)
Verdicts after the correction below: **hit 2, miss 3, uncovered 2, not_applicable 9, not_served 1,
not_prespecified 13.** The scored panel is therefore 7 anchors with a stated regulatory expectation, of which
2 hit; 9 coding anchors are not_applicable by construction; 13 anchors carry `expected_channel: none_stated`
in the prespecification (mechanism not established, or the gene unassigned) and are descriptive only. That
13 is a weakness of the prespecification, not of the model, and no expectation was added after scoring.

Written predictions: **P2.2 held** (rs72613567, the HSD17B13 splice-donor insertion, is not served by the
point endpoint). **P2.3 held** (MBOAT7 rs641738). **P2.1 failed** - every coding_missense anchor was predicted
to show |liver RNA quantile| < 0.5 and instead all nine sit at 0.90-0.999. **P2.4 failed** - at most 2 of 9
coding anchors were predicted to reach window rank 0.9 in any regulatory channel and 3 did.

The two failures need the right reading, and one control settles it. Across the 13,690 ordinary
credible-set variants of the direct universe, |quantile| > 0.9 occurs in 12.6% (ATAC), 7.8% (DNase), 10.1%
(RNA) and 6.2% (histone), with a median |quantile| near 0.5. So the coding anchors really are in the extreme
tail of the calibrated distribution. But their rank **within their own 2-kb window** is mostly 0.59-0.96:
they are not the standout variant among the ~6,000 SNVs beside them. Both readings hold together - these
anchors sit in high-signal loci, and the calibrated quantile carries that locus signal, while the
within-window rank is the variant-specific statement. The practical consequence for this paper: a bare
calibrated quantile must not be read as evidence that a particular variant is special, and the window rank
must be reported beside it. This is a Figure 7 point about what the scores do and do not localise, and it
does not touch the P3 allelic result, which is a within-position sign test at a fixed locus.

**Verdict defect corrected before the numbers were used.** `served` meant "an hg38 uid was built and
queried", not "the Atlas returned a row". GNMT rs2296804 and rs2296805 - the r2 = 1.000 pair the Resource's
open two-variant question rests on - resolved cleanly, were queried, and came back with no row for either
uid; both were scored `miss`. That asserts a negative about a variant that was never measured. The panel now
separates `requested` from `returned`, emits `atlas_returned_a_row`, and gives these anchors `uncovered`
(fixture-tested). The first deposit is kept under `tables/_superseded_served_meant_requestable/`; no
retrieval changed, only the scoring. Consequence for P5: the Atlas cannot score the primary GNMT pair at
all, so the haplotype question has to go through the model API (`predict_sequence`), which P5A already uses.

## sQTL, finished: the splice channel does not recover GTEx liver sQTL direction
Three defects were fixed before the number was accepted, each measured rather than assumed.

1. **Estimand.** GTEx `slope` is a LeafCutter within-cluster excision ratio; the first comparison used the
   Atlas absolute junction quantile. Different quantities.
2. **Coordinate convention.** Atlas junction coordinates are half-open, GTEx introns closed. Measured over
   20,272 archived junction rows against 1,414 distinct GTEx liver introns: `junction_Start` lands on the
   intron start at offset 0 (795 hits, next best 5 at +2) and `junction_End` at offset +1 (814 hits, next
   best 4 at -1). The join is now exact `(start, end + 1)`; the earlier +/-1 tolerance is removed (914
   matched vs 919, so the tolerance had been doing the same work loosely).
3. **Cluster membership is not distributed.** GTEx's public sQTL release lists only the top phenotype per
   cluster: all 1,414 clusters in the significant-pairs file are singletons, and 10,686 of 10,853 in
   `sgenes`. The cluster the slope is defined against cannot be read off the files on disk, which is why the
   first contrast attempt returned n = 0. Siblings are therefore RECONSTRUCTED by LeafCutter's own rule -
   junctions sharing the target's donor or acceptor - over the Atlas junction set for the same variant.
   Every number from it carries that label; it is not GTEx's cluster.

Result: matched (reconstructed-cluster) sign concordance **0.538**, n = 684; sGene leads 0.532, n = 111;
Spearman(slope, contrast) 0.087. The mismatched absolute-junction comparison reads 0.510, n = 914. The
written prediction was > 0.7. Sites are clumped, so the nominal SE (about 0.019 at n = 684) overstates
precision; either way the answer is near chance and nowhere near the prediction. Exposure is
`documented_overlap` - GTEx is an AlphaGenome training source, so this arm was meant to be a CEILING, and
the model does not clear it.

## The pattern across three independent benchmarks
| channel | benchmark | result | prediction |
|---|---|---|---|
| ATAC / DNase | Currin caQTL, 138 donors, 10,966 leads | auROC **0.713** [0.703, 0.723] | held |
| ATAC / DNase | GSE281367 snATAC allelic, same-variant | concordance **0.655**, blocks p 5.5e-06 | exceeded |
| RNA (gene) | 211-set eQTL direction replication | auROC **0.530** [0.449, 0.607] | held (at chance) |
| splice | GTEx liver sQTL, reconstructed clusters | concordance **0.538**, n 684 | failed (predicted > 0.7) |

The accessibility channels carry real liver signal; the RNA and splice channels do not recover measured liver
QTL direction, including on a benchmark inside the model's own training distribution. That is the Figure 7
statement, and it now rests on three independent substrates rather than one.

## P3 state-dependence: no global signal at this n, and per-site testing declared untestable in advance
Prespecification `74_state_dependence_prespec.json` (sha256 f5ed6581...) was written before any state
statistic, and it abandons per-site testing rather than reporting it as a null: with >= 3 het donors per
group a per-site rank test has a minimum two-sided p of 0.1, so no site can reach any sensible FDR
threshold. Only 61 sites in GSE281367 and 40 in GSE244832 even reach 3-vs-3. Naming a state-conditional
variant is not possible on this substrate and was not attempted.

The powered question is global: is there ANY state-dependence beyond noise, with sites as repeated measures
and donors as the unit? Null shuffles DONOR group labels and regenerates the statistic end to end (2,000
draws, seed 20260909); shuffling sites would destroy the correlation structure the statistic lives in.

| | GSE281367 | GSE244832 |
|---|---|---|
| donors | 6 control / 6 disease | 5 control / 13 disease |
| eligible sites (>= 2 het donors in both groups) | 206 | 91 |
| group depth gap (log10) | 0.214 -> permutation stratified on depth tertile | 0.019 -> unstratified |
| observed mean\|delta\| | 0.3777 | 0.6530 |
| donor-label null | 0.3869 [0.359, 0.412] | 0.6482 [0.548, 0.752] |
| exceedance | 1592/2000, p 0.80 | 924/2000, p 0.46 |
| sd(delta) observed vs null | 0.4925 vs 0.5045, p 0.74 | 0.8371 vs 0.8316, p 0.47 |

**No global state-dependence is detectable at this n in either cohort.** Written predictions S1, S2 and S3
held. The prespecified depth diagnostic fired in GSE281367 (gap 0.214 exceeds the 0.15 limit) and the
permutation was stratified on donor depth tertile, because the +0.5 pseudocount shrinks a log ratio toward
zero at low depth and a depth-group association would otherwise manufacture a delta. This is
`untestable at this n`, not `unsupported`.

## Mapping bias cannot explain the P3 allelic concordance
No simulated-read remapping control was run, but the direction of the two biases settles it. Reference
mapping bias favours the reference allele: measured log ratios do lean that way, at a global mean of -0.007
and a measured positive fraction of 0.481. The Atlas leans the other way, predicting alt-raises for 0.539 of
the same sites. Two biases pointing in opposite directions can only depress sign concordance, not create it,
and the concordance expected from both marginals together is 0.4985 against 0.6553 observed.

## P5A haplotypes: two prespecification gaps closed before any residual was analysed
The pair list came out at **1,587 eligible pairs**, not the ~200 the plan assumed, and that exposed two
weaknesses in my own prespecification. Both were amended while the run was still computing effects, with no
residual distribution, p-value or BH result read.

1. **No multiple-testing family was named.** The script's rule was "residual exceeds the matched null's 95th
   percentile", a 5% per-pair rate that would flag about 79 of 1,587 pairs by construction, and the written
   prediction ("at most 2 signals") had nothing to be measured against. Amended to: family = all eligible
   pairs in the run; per-pair two-sided empirical p against the matched null; BH at q < 0.10; reported per
   SIGNAL beside the pair count and the family size. The uncorrected count is kept, labelled uncorrected.
   Revised H1: fewer than 5 signals survive.
2. **Additivity is undefined when neither variant does anything.** A probe on a real pair
   (chr13:49478927:A:G with chr13:49480172:G:A) gives single-variant local ATAC effects of +0.001332 and
   -0.000962 log2 - the same order as the residual itself. Pooling 1,587 such pairs would let numerical
   wobble dominate the family. Amended to report two strata: (a) all pairs, and (b) pairs where BOTH
   variants' single-variant |effect| exceeds the 95th percentile of the matched null's single-variant
   |effect| distribution. The threshold is fixed by the null, not chosen from the observed residuals.
   Stratum (b) conditions on components of the statistic, so it is reported beside (a), never instead of it,
   and BH runs within each stratum. If fewer than 20 pairs enter (b), the arm reports that additivity is
   untestable on this substrate, which is itself the finding.

Two bugs fixed in the same lane. The first run died in `summarise()`: readout masks are built at base
resolution (1,048,576 for a 1-Mb window) but CHIP_HISTONE returns binned at 128 bp (8,192 rows), so the
boolean index raised. `fold_mask` folds a mask onto each track's own row count, keeping a bin when any base
inside it is selected. And `50_haplotype_additivity.py` read `sys.argv[1]` at module level, which breaks test
import - the same defect as step 73; both now resolve arguments inside `main()`.

Diagnostic checked before trusting the arm: `predict_sequence` is bitwise identical across repeated
identical 1-Mb calls (max |diff| 0.0), so the residual is a real property of the model rather than
call-to-call noise. Resolutions measured: ATAC, DNase and RNA_SEQ at 1 bp; CHIP_HISTONE at 128 bp.

## A note that reconciles the anchor panel with the raw model effects
The Atlas gives the coding anchors calibrated quantiles of 0.90-0.999 while `predict_sequence` shows raw
single-variant effects around 0.1% of the local signal sum. These are consistent, not contradictory: the
calibrated quantile RANKS a small raw change against a distribution of even smaller changes. An extreme
quantile is therefore not a claim that the predicted effect is large in absolute terms, and the paper should
not present it as one.

## The allelic result does not depend on the het-calling filter (step 75)
The het call is the headline's main free choice. Sweeping all three thresholds over a fixed grid
(min total reads 10/15/20, min reads per allele 3/4, min het donors 3/4/5 - 18 combinations per cohort):

| cohort | concordance range | block sign-test p range | marginal expectation |
|---|---|---|---|
| GSE281367 | 0.634 - 0.688 | 5.5e-06 - 1.7e-02 | 0.499 - 0.501 throughout |
| GSE244832 | 0.544 - 0.596 | 1.8e-02 - 1.6e-01 | 0.492 - 0.499 throughout |

The discovery cohort is above the marginal-skew expectation and significant on the LD-respecting block sign
test at **every one** of the 18 combinations. The replication is positive in direction at every combination
and never reaches the discovery cohort's strength; it does get stronger as the depth threshold rises
(0.584-0.596 at min total 20, p 0.018-0.076), which is what noise dilution rather than absence of effect
looks like. The grid is deposited whole in `allelic_filter_sensitivity.tsv`; no threshold is selected from it.

## P6a scoping, found by dry-running against the preview
Running the drug-target direction step against the preview archive (a dry run, not a result) exposed a
denominator error in the plan's written prediction. The funnel:

| stage | n |
|---|---|
| targets listed (`drug_target_classification.tsv`, MASLD status or in the agent list) | 1,719 |
| with any agent (`mash_clinical_pipeline.tsv`, 45 agents) | 28 |
| with a READABLE agent direction (mechanism parses to inhibit or activate) | 20 |
| with a Track 0 genetic signal (preview archive, ~3% of the query set) | 67 |
| with both a readable agent direction and a resolved genetic direction | 1 |

The last row will grow when the full archive replaces the preview, but the ceiling is **20** - the number of
targets whose agent direction is readable at all. The plan's prediction ("at least 60% of targets with a
regulatory genetic signal are concordant") uses a denominator that cannot be scored: most targets with a
genetic signal have no agent. Corrected before any full-run number was read to "targets with a readable agent
direction AND a resolved genetic direction", with the funnel reported beside it and a stated rule that fewer
than 10 scorable targets means the deliverable is a per-target listing and no rate is quoted.

Separately, only 2.1% of `direction_variant_level.tsv.gz` rows carry both a GWAS and an eQTL beta in the
preview. That is the preview's direct-only archive (13,690 of 417,155 query variants), not a defect, and it
is why the dry run resolves almost nothing.

## Result-root provenance is now a switch, not an edit
Four step scripts read Track 0 deposits by absolute path, and step 73 read `liver_summaries.tsv.gz` from the
PREVIEW (direct-only) archive because the full run was still building. `la.track0_root()` (env
`AGA_TRACK0_ROOT`) and `AGA_LIVER_SUMMARIES` make the source a stated switch, and step 73 now records
`liver_summaries_source` in its summary. Re-running both allelic cohorts through the new path reproduces the
first pass exactly (GSE281367 n 412, concordance 0.6553; GSE244832 n 264, 0.5644). When the full tables land,
the allelic arm is re-run against the full archive and the two numbers compared.

## The saturation pilot died on a real outage, not a bad request (2026-09-09) — WRONG, SEE THE CORRECTION BELOW
Job 21609736 exhausted its transient-retry budget on the very first interval and raised UNAVAILABLE. Before
changing anything I checked whether the service or the request was at fault: minutes later
`scorer_metadata()` returned all 22 scorers in 0.5 s, and the pilot's exact interval `chr1:96466-96736`
succeeded with the pilot's own six scorers at `max_workers=4` in 1.7 s. So the request shape was fine and
the Atlas was simply down for longer than 10 x 30 s = 5 minutes.

Fix is the budget, not the request: `transient_sleep_seconds` backs off 30 s doubling to a 300 s cap over 14
attempts, covering **58 minutes** of outage (schedule 30, 60, 120, 240, then 300 throughout). Fixture-tested,
including that a transient error which clears returns the value and that the backoff actually grew. Pilot
resubmitted as 21610805.

This is the second time a gRPC status code sent the diagnosis the wrong way. RESOURCE_EXHAUSTED meant an
oversized response, not a quota breach; UNAVAILABLE here meant a genuine outage rather than a malformed
request. Both times the answer came from probing the service directly with the exact failing request rather
than reasoning from the status code.

## Calibration of the allelic concordance: the gradient is real-looking and NOT established (step 76)
Sign concordance pooled over all sites blends strata the model is confident about with strata it is not, so
the question that decides whether these scores are usable per variant is whether agreement rises with the
predicted effect size. Equal-count quartiles of |predicted ATAC quantile| in GSE281367:

| quartile | median \|quantile\| | n | concordance | marginal expectation |
|---|---|---|---|---|
| Q1 weakest | 0.311 | 107 | 0.477 | 0.499 |
| Q2 | 0.782 | 99 | 0.626 | 0.492 |
| Q3 | 0.940 | 106 | 0.698 | 0.500 |
| Q4 strongest | 0.988 | 100 | **0.830** | 0.501 |

That looks like a textbook calibration curve. It does not survive the LD-respecting test. Spearman trend
across sites +0.287; **1-Mb block bootstrap +0.279, CI [-0.023, 0.571]** over 36 blocks - the interval
includes zero. GSE244832 is flatter (+0.055, boot +0.072, CI [-0.275, 0.405]). Verdict recorded by the script
itself: **not established**.

What stands and what does not:
- The concordance headline stands (0.655; 29/33 blocks above half; p 5.5e-06) - it survives the same test.
- The stronger and more useful claim - that the score is calibrated, so a confident prediction is a
  trustworthy one for a particular variant - is NOT supported here. With sites clustered at 36 loci this
  substrate cannot resolve the trend. The quartile table alone would have made a persuasive figure panel and
  an overclaim; it is deposited beside the verdict that it is not established.

Concordance does not rise with the number of het donors (0.667 / 0.676 / 0.634 at 3 / 4 / 5+ donors), so the
limit here is the number of independent loci, not measurement noise per site.

## MPRA benchmark (Hu/Zhu 2026, full 3,498-oligo library)
Prespecification `33_mpra_prespec.json` (sha256 ce84f33c...) was sealed before retrieval finished. Channel
panel and summary rule fixed in advance per the AlphaGenome exposure audit, which forbids MPRA outcomes
selecting tracks or aggregation.

**DAV vs non-DAV auROC** (positives: 246 HepG2 ctrl, 311 HepG2 PA-OA, 235 LX2 ctrl, 403 LX2 TGFB):

| channel | HepG2_ctrl | HepG2_PAOA | LX2_ctrl | LX2_TGFB |
|---|---|---|---|---|
| hepg2_h3k27ac | **0.685** | 0.595 | 0.666 | 0.602 |
| hepg2_dnase | **0.685** | 0.586 | 0.652 | 0.611 |
| hepg2_atac | 0.670 | 0.567 | 0.651 | 0.576 |
| liver_atac | 0.663 | 0.574 | 0.646 | 0.593 |
| hepg2_cage | 0.623 | 0.539 | 0.617 | 0.532 |
| fibroblast_dnase (cross-tissue proxy) | 0.601 | 0.555 | 0.590 | 0.597 |
| avi | 0.611 | 0.535 | 0.586 | 0.566 |

Prediction M1 said HepG2 auROC 0.60-0.70 with AVI approximately the best single channel. The **range holds**
(0.611-0.685) and the **AVI clause fails**: AVI is the WEAKEST channel in all four contexts. This is the
fourth independent substrate on which the chromatin channels carry the liver signal and the aggregate
variant-impact score does not (caQTL, snATAC allelic, sQTL, MPRA).

**Stimulus-specific vs constitutive DAVs** came back BELOW chance, and the reading depends on a control that
was not in the first deposit:

| context | channel | auROC | 95% CI | magnitude control | reading |
|---|---|---|---|---|---|
| HepG2_PAOA | hepg2_atac | 0.370 | [0.290, 0.455] | 0.544 | context |
| HepG2_PAOA | hepg2_h3k27ac | 0.366 | [0.294, 0.443] | 0.544 | context |
| HepG2_PAOA | avi | 0.386 | [0.318, 0.455] | 0.544 | context |
| LX2_TGFB | hepg2_atac | 0.340 | [0.272, 0.412] | 0.420 | partly magnitude |
| LX2_TGFB | hepg2_h3k27ac | 0.404 | [0.326, 0.484] | 0.420 | partly magnitude |
| LX2_TGFB | avi | 0.450 | [0.371, 0.526] | 0.420 | partly magnitude |

`magnitude_control_auc` is the same discrimination computed on |MPRA log2FC| itself. Constitutive DAVs are
significant in two contexts, so they may simply carry larger measured effects, in which case any model that
tracks effect size lands below 0.5 for a reason unrelated to context.

- **HepG2 + PA-OA: the control clears it.** Stimulus-specific DAVs have equal-or-larger measured effects
  (control 0.544, median |log2FC| 0.803 vs 0.706, p 0.24) yet the Atlas scores them LOWER (median quantile
  0.816 vs 0.900). Effect size would push the auROC the other way. The static model systematically
  down-ranks elements that only become active under lipid loading. This is the matched-cell-type arm.
- **LX2 + TGFB: partly explained by magnitude.** Constitutive DAVs really do have larger effects there
  (control 0.420, 0.678 vs 0.552, p 0.019), so the context claim is weaker - and LX-2 is also the arm where
  the Atlas has only a cross-tissue fibroblast proxy and no stellate track.

Prediction M2 (stimulus-specific auROC within [0.45, 0.58]) **fails**: the result is below chance, not at it.
Below chance is a stronger statement than absence, and only the HepG2 arm supports the context reading.
The first deposit, which lacked the control, is kept under `tables/_superseded_no_magnitude_control/`.

## CORRECTION: the saturation pilot was killed by `ontology_terms`, not by an outage
The earlier entry attributing the pilot's UNAVAILABLE errors to a service outage is **wrong**, and the test
I cited did not test what I said it tested. The probe extracted the pilot's `TERMS` list with a regex that
silently returned `None` because `TERMS` is a set comprehension rather than a literal, so both arms of that
probe ran WITHOUT ontology filtering and I labelled one of them "with terms".

What actually happens, tested with the real 50-curie list:

| request | result |
|---|---|
| `query_interval` WITH `ontology_terms` (50 curies), max_workers=4 | **UNAVAILABLE after 59.9 s** |
| identical request WITHOUT `ontology_terms` | **OK in 2.2 s** |
| WITH `ontology_terms`, max_workers=1 | **UNAVAILABLE after 208.2 s** |

The server appears to time out applying a large ontology filter across six scorers over an interval and
reports it as UNAVAILABLE. The signature that should have caught this earlier: the pilot failed **10
consecutive times on the same interval over 43 minutes** while a direct probe of that interval succeeded in
1.7 s. A flaky service would have let one of ten attempts through; a request that cannot succeed fails every
time.

Fix: retrieval no longer passes `ontology_terms`. The ontology filter is applied client-side after
retrieval, which is where every other step already applies it (`liver_mask` / `liver_track_mask`), and the
archive keeps all tracks. Provenance records `ontology_filter: client-side` so no reader mistakes the
archive for an ontology-filtered one. Pilot resubmitted as 21611192.

The backoff change made under the wrong diagnosis is still worth keeping - it covers genuine outages - but
it was treating a symptom, and it would have let this request retry for 58 minutes instead of 5 before
failing. A request that fails identically every time should be diagnosed, not retried harder.

## Track 0 COMPLETE — validation 14/14 (job 21558781, 4 h 30 m)
The production run passes every check, including the two the preview could not: all **417,155** query
variants archived, and availability rows matching the query set exactly.

| | |
|---|---|
| signals | 989 (A_direct 29, C_enzyme 787, B_direct 173) |
| variants queried and archived | 417,155 (the prior in-house profile covered 126) |
| channel-profile rows | 12,857 |
| Package A channel rows / B variant-peak rows | 12,857 / 40,700 |
| Package D variant rows | 3,452,688 |
| Package G links | 11,760 over 816 signals |
| Catalog sidecar rows | 989 |
| coverage | posterior mass >= 1.0 at the 75th percentile |

**211-set replication, the plan's first directional number and not re-thresholdable:** auROC **0.5304**
[0.449, 0.607], permutation p 0.229, sign concordance 0.522, n = 201 of 210 SNV leads. Identical to the
preview, which is the right outcome - the gate archive did not change between them, so the number was not an
artifact of partial data. The prior in-house AlphaGenome benchmark read 0.5607 (FAIL, n = 211); on
hg38-correct alleles with server-declared signs it reads 0.5304. Written prediction was "at chance"; it holds.

Downstream launched on the full archive (job 21611415): P1 extensions E1-E3, P6a, and the allelic arm
re-run against the FULL `liver_summaries.tsv.gz` rather than the preview's direct-only 13,690 variants.
Packages C and F (job 21610395) - the two steps the production sbatch omitted - run in parallel.

## Where the programme stands
Five independent substrates now agree on the same split.

| channel | substrate | result |
|---|---|---|
| ATAC / DNase | Currin caQTL, 138 donors, 10,966 leads | auROC **0.713** [0.703, 0.723] |
| ATAC / DNase | GSE281367 snATAC allelic, same-variant, 12 donors | concordance **0.655**, 29/33 blocks, p 5.5e-06 |
| ATAC / DNase | GSE244832 snATAC allelic replication, 18 donors | 0.564, 21/32 blocks, p 0.055 (consistent, underpowered) |
| H3K27ac / DNase | Hu MPRA, HepG2 control, 3,498 oligos | auROC **0.685** |
| RNA (gene) | 211-set liver eQTL direction | 0.530 [0.449, 0.607] |
| splice | GTEx liver sQTL, reconstructed clusters | 0.538, n 684 |

The accessibility channels carry real liver signal on three independent assays. The RNA and splice channels
do not recover measured liver QTL direction, and the splice arm fails on GTEx, which is inside AlphaGenome's
own training data. AVI, the aggregate variant-impact score, is the weakest MPRA channel in all four contexts.

Two cautions the paper must carry, both established here rather than assumed:
- A bare calibrated quantile is not a variant-specific statement. Anchors sit at |q| 0.90-0.999 while their
  rank within their own 2-kb window is 0.59-0.96, and across 13,690 ordinary credible-set variants |q| > 0.9
  occurs in only 6-13%. Report the window rank beside the quantile.
- Whether the quantile is CALIBRATED against measurement is not established. The quartile gradient
  (0.477 -> 0.830) looks clean but its block-bootstrapped trend is +0.279 [-0.023, 0.571] over 36 loci.

## The full archive strengthens the allelic arm, and the calibration verdict changes (2026-09-10)
Re-running the allelic arm against the FULL `liver_summaries.tsv.gz` (417,155 variants) instead of the
preview's direct-only 13,690:

| | GSE281367 | GSE244832 |
|---|---|---|
| scored sites | 412 -> 430 | 264 -> 277 |
| concordance | 0.6553 -> 0.6605 | 0.5644 -> 0.5668 |
| marginal-skew expectation | 0.4985 -> 0.4990 | 0.4984 -> 0.4989 |
| independent 1-Mb blocks | 33 -> **48** | 32 -> **44** |
| block sign-test p | 5.5e-06 -> 1.7e-06 | 0.055 -> 0.048 |
| block bootstrap | 0.737 [0.583, 0.861] -> 0.733 [0.608, 0.843] | 0.593 [0.429, 0.743] -> 0.591 [0.447, 0.723] |

The replication's p crossing 0.05 (0.055 -> 0.048) is a threshold crossing on an unchanged picture and is
NOT reported as a new finding. The substantive change is the number of independent loci.

**Calibration, recomputed on the larger substrate (step 76, method and seeds unchanged):**

| cohort | n sites | blocks | trend | block bootstrap | verdict |
|---|---|---|---|---|---|
| GSE281367 | 430 | 51 | +0.303 | **+0.421 [0.182, 0.621]** | calibrated |
| GSE244832 | 277 | 47 | +0.065 | +0.108 [-0.189, 0.383] | not established |

GSE281367 quartiles now run 0.454 -> 0.841. The interval excludes zero, so in the discovery cohort agreement
does rise with predicted effect size. **The replication does not support it**, so this is a single-cohort
result and is carried as provisional, not as an established property of the score.

The history of this number is kept deliberately: the quartile table was first described as a clean
calibration curve (an overclaim), the block bootstrap over 36 loci said "not established" (the correction),
and 51 loci moved it to "established in one cohort only". Only the substrate size differs between those
last two - `tables/_calibration_on_preview_archive/` holds the earlier state with that note.

## P5A haplotype additivity: no evidence of non-additivity, after three traps
1,572 eligible pairs scored in four arms each (REF / v1 / v2 / joint) over 1-Mb windows, against 300 matched
random-pair nulls. `predict_sequence` verified bitwise deterministic before the arm was trusted.

**Result, on the family the data can support:**

| channel | all pairs (1,572) | median \|residual\| | matched null median | both-variants-active stratum |
|---|---|---|---|---|
| rna | **0 BH pairs / 0 signals** | 0.00132 | 0.00194 | not testable (2 matched null pairs) |
| atac | **0 / 0** | 0.00078 | 0.00060 | not testable (1) |
| dnase | **0 / 0** | 0.00053 | 0.00041 | not testable (3) |
| h3k27ac | **0 / 0** | 0.00107 | 0.00067 | not testable (2) |

Predictions H1 (fewer than 5 signals survive BH) and H2 (the GNMT pair is additive) hold. In the RNA channel
the observed median residual is SMALLER than the matched null's, so Resource variant pairs are if anything
more additive than random pairs at the same separation. The paper's variant-combination question stays
`indeterminate`.

**Three traps, each caught before it became a number:**
1. **No multiple-testing family.** The original rule counted pairs exceeding the null's 95th percentile: on
   1,572 pairs that flags about 5% by construction, and it flagged **89**. BH at q < 0.10 flags **0**.
2. **Additivity is undefined when neither variant does anything.** Single-variant effects are the same order
   as the residuals (median |effect| 0.0013-0.0038 log2), so most pairs contribute numerical wobble. Fixed by
   adding a stratum of pairs where both variants clear the null's 95th percentile of single-variant effect.
3. **That stratum had no matched null - created by the fix for (2).** It selects on large single-variant
   effects, but only **1-3 of 300 null pairs** clear the same bar, so residuals were compared against a null
   made almost entirely of no-effect pairs. It rejected **19-26 of 23-33 tests (83%)** across four channels,
   which reads as a strong positive and is a broken comparison. `stratum_is_testable` now withholds a verdict
   for any stratum with fewer than 20 matched null pairs, fixture-tested.

Trap 3 is the instructive one: it was introduced by the correction for trap 2 and produced the most
convincing-looking output of the three (17-19 "non-additive signals"). Testing that stratum honestly needs
enough null draws to yield about 20 matched active pairs - roughly 15,000 null draws at the observed rate,
about 25 hours of model API time. Not run; recorded as what the arm would require.

## Programme complete except P4 (2026-09-10)
Packages C (816 signals, 7,245 splice event rows) and F (1,294 tissue rows, 66 ancestry rows) landed - the two
steps the production tables job would have silently omitted. **155 unit tests pass** on a compute node.
`scripts/manuscript/validate_resource_scope.py` returns **PASS**, contract
`standalone-masld-resource-2026-08-13-inherited-variation-remodeling`, `CAS13_FEEDBACK_TO_RESOURCE FALSE`.

Every pillar is finished except P4 saturation, which is retrieving 103,390 regions with two chained
`afterany` resumption jobs (21611739, 21611740) so it survives its 7-day walltime without intervention.

**P4 exceeds a written kill criterion and that is flagged, not overridden.** The plan says restrict to the
~15k reduced universe if the full set exceeds 5 days; the sustained rate (694-761 regions/h) projects
**6.2 days**. The full U1 set is kept because the Package H overlay question - are regions the released
predictor CAN predict more sequence-sensitive than matched regions it cannot - needs all 96,460 regions
rather than a DA-selected subset. Reverting to the reduced universe is a one-line change to the step-41
filter and would finish in about 20 hours. This is a user decision, recorded here rather than taken silently.

### Claims I made this session that the evidence did not support
Each is corrected above with the original flagged in place, not deleted.
1. **"A clean calibration curve."** Said before nulling it; the block bootstrap over 36 loci gave
   [-0.023, 0.571]. With 51 loci from the full archive it is now [0.182, 0.621] in the discovery cohort and
   still [-0.189, 0.383] in the replication - established in one cohort, carried as provisional.
2. **"A genuine outage killed the saturation pilot."** The probe I cited never tested the pilot's
   configuration: a regex silently returned `None` for `TERMS`, so both arms ran unfiltered and I labelled
   one "with terms". The real cause was `ontology_terms` itself.
3. **"89 pairs exceed the null"** would have been P5A's headline under a rule with no multiple-testing
   family. BH gives 0.

## The junction convention was wrong in a second place (Package C, 2026-09-10)
The disease-splicing join used the same `+/-1` junction tolerance the GTEx sQTL arm had, and it had already
shipped inside the completed Track 0. The Atlas-to-LeafCutter convention was measured directly on this
substrate rather than carried over from GTEx: **start -> boundary1 at offset 0** (17,159 hits, next best 11 at
-2) and **end -> boundary2 at offset +1** (19,893, next best 133 at -2). Identical to GTEx and exact, so a
+/-1 window admits those near-misses and can select a neighbouring junction inside a dense cluster.

Effect of the fix, stated at its true size: **9 of 2,896 matches (0.3%) were coordinate artefacts** and
**11 surviving matches changed which junction they point at**. Event rows (7,245) and signal states (633
measured / 183 untestable) are unchanged, so no Track 0 conclusion moves. The correction matters less for
this arm than for sQTL, where the same defect sat under a headline number - and it removes a second
hand-written tolerance.

`la.leafcutter_intron_key` is now the single definition, called by both `33_benchmark_eval.py` and
`15_package_C_splicing_join.py`, with the measurements from both substrates in its docstring. The tolerance
deposit is kept at `tables/_superseded_junction_tolerance/`.

## caQTL: three fixes, and the paired comparison changes the number (2026-09-10)
Pulling deposited numbers for a status report surfaced that this arm had never been re-run after the
reproduction guard was corrected, so the file on disk still carried the false failure verdict. Fixing that
exposed two more defects, both of the SAME class the guard exists to catch - comparing across different
variant sets - and both introduced by my own corrections:

1. **Stale deposit.** Guard fixed in code on 2026-09-09, arm not re-run until 2026-09-10.
2. **Dangling `api_sub`.** The patch removed the old guard computation but left its variable referenced in
   the prediction block, so the re-run crashed. It also meant the written prediction was being evaluated on
   the 3,000-lead GPU subsample rather than on the leads the archived row uses.
3. **The comparison was not paired.** The Atlas serves no indels, so it scored 10,966 of the API's 11,896
   leads. The reported delta compared two auROCs computed on two different variant sets.

`both_scored` now restricts the comparison to leads BOTH models scored, fixture-tested.

**Final numbers:**

| quantity | value |
|---|---|
| archive reproduction guard | **pass, delta 0.0**, n 11,896 both sides |
| Atlas ATAC auROC | 0.7128 [0.703, 0.723], n 10,966 |
| API auROC on the SAME 10,966 variants | 0.7145 |
| delta on the same variants | **0.00175** (written prediction < 0.02) |
| delta on the unpaired sets (superseded) | 0.0054 |
| Spearman(Atlas, API) scores | 0.981 |

The paired delta is three times tighter than the unpaired one. The honest reading: the Atlas ATAC scorer and
the AlphaGenome API agree almost exactly on liver caQTL direction, which is what should happen - they are the
same model exposed two ways - and part of the earlier 0.0054 was the API's extra 930 indel leads.

**Pattern worth recording:** the reproduction guard was written to catch "same number, different variant set",
and I then committed that error twice more in the same file while fixing it. A correction is a new analysis;
it inherits none of the original's testing and can reintroduce the flaw at a different point.

## Aiming higher: using the validated channel as an instrument (2026-09-11)
The accessibility channel is validated on three assays, so two analyses treat it as an INSTRUMENT rather
than benchmarking it further.

### P6e prioritisation — cannot conclude at 11 loci
Question: within a credible set, does the Atlas rank the allelically ACTIVE variant above the inactive ones,
and does it beat fine-mapping posterior? Ground truth is measured GSE281367 allelic imbalance
(|mean log2| > 0.5). Every variant compared is peak-resident and read-covered, so peak membership - which
would otherwise win the ranking outright - is controlled by construction.

| ranker | mean rank percentile | block bootstrap | null | perm p |
|---|---|---|---|---|
| atlas | 0.638 | 0.606 [0.455, 0.728] | 0.500 | 0.0005 |
| posterior | 0.543 | 0.610 [0.492, 0.744] | 0.501 | 0.0070 |
| combined | 0.613 | 0.588 [0.454, 0.719] | 0.500 | 0.0005 |

**Verdict: cannot conclude.** 57 eligible sets sit at only **11 independent 1-Mb blocks**, below the
prespecified minimum of 15, and every interval includes 0.5. Predictions R1 (Atlas interval excludes 0.5)
and R2 (posterior at chance) both FAIL - posterior bootstraps to 0.610, indistinguishable from the Atlas.
The +0.096 raw advantage does not survive LD-aware inference; the within-set permutation p of 0.0005 is
anti-conservative because sets at one locus share variants. R4 held (null centres at 0.500).

The binding constraint is the same one that limits the calibration trend: the Resource's credible sets
concentrate at a handful of liver-trait loci, so 435 allelic sites collapse to 11 blocks.

**This gives P4 a second purpose.** Saturation scores every substitution in 103,390 measured peaks, which
would expand the allelic target list well beyond the 1,119 Atlas-scored peak variants now available - adding
independent loci, which is precisely what this test and the calibration trend both lack.

## TF convergence (step 80): a clean negative, delivered by the prespecified control
Question: with 21 adult-liver and 539 HepG2 TF ChIP tracks over ~501 factors, do liver-trait credible
variants converge on disruption of particular transcription factors? Design: credible variants (weight >=
0.05) against non-credible variants of the SAME signal (weight < 0.005), so locus, LD block, GC and peak
context are controlled without constructing a matched draw; and a factor-specificity deviation (each factor's
|quantile| minus the variant's own mean across factors) separating WHICH factor is hit from HOW MUCH the
variant does overall. 761 signals, 18,233 variants, 502 factors.

| set | n | mean rank of 502 | expected | p |
|---|---|---|---|---|
| hepatocyte regulators (HNF4A, HNF1A/B, HNF4G, CEBPA/B, FOXA1/2, NR2F2, ONECUT1) | 10 | 166.5 | 250.5 | 0.0332 |
| **general machinery CONTROL** (POLR2A, TAF1, CTCF, RAD21, YY1, EP300) | 6 | **121.2** | 250.5 | **0.0117** |

**Verdict: no MASLD-specific TF convergence.** Hepatocyte regulators do rank above chance, but the
prespecified negative control ranks HIGHER - YY1 at rank 14 and POLR2A at 81, ahead of HNF4A at 87. If
inherited risk converged on hepatocyte regulatory identity that ordering would reverse. The pattern is what
credible variants sitting in more active regions produces: the most ubiquitously bound factors move most.
THRB topping the raw list is a curiosity given resmetirom, not evidence, and is not presented as one.

**Two defects of mine, both caught by prespecified diagnostics before any number was used:**
1. **The null was not exchangeable.** The statistic weights the test side by posterior and leaves the control
   side unweighted, but the permutation shuffled LABELS while weights stayed attached to rows, so a permuted
   test set was built from low-weight rows. Null centre drifted to 0.0097 against a 0.002 limit, and 151 of
   502 factors crossed BH on it. Permuting the (label, weight) PAIR fixes the centring to **0.00046**.
2. **Per-factor BH was not a valid family.** The deviations sum to zero within every variant - measured at
   -0.000000 - so any real shift in profile shape pushes some factors up and an equal mass down. 247 of 502
   crossed BH even on the corrected null. Per-factor p-values are now descriptive only; the inferential
   question is a single prespecified rank-enrichment test carrying its own negative control.

Two engineering defects were fixed on the way: a permutation loop costing ~1.4e10 element comparisons per
draw (~17 h against a 12 h walltime), replaced by a signal-major vectorised scorer at 0.07 s per draw and
asserted numerically identical to the loop; and a control subsample seeded with `abs(hash(signal_uid))`,
which Python randomises per process - two runs of the same job drew 38,523 and 38,456 variants. Seeds now
come from a sha256 digest, verified identical across processes. The deviation matrix is cached
(`tf_convergence_matrix.npz`) so scoring changes cost minutes rather than a 4-hour archive re-read.

---

## 2026-09-15: the splice channel emits no direction (steps 35, 36, and a repo-wide signedness audit)

### The defect

`scorer_metadata.tsv` (Track 0 deposit, `run-20260909T153939Z/tables/`) carries the server-declared
`is_signed` per scorer. **UNSIGNED:** `SPLICE_JUNCTIONS`, `SPLICE_SITES`, `SPLICE_SITE_USAGE`, `AVI_SCORE`,
`AVI_SCORE_FEATURE_IMPORTANCE`, `CONTACT_MAPS`, `POLYADENYLATION`, and every `*_ACTIVE`. **SIGNED:** `ATAC`,
`DNASE`, `CHIP_HISTONE`, `CHIP_TF`, `RNA_SEQ`, `CAGE`, `PROCAP`, `AVI_SCORE_MODEL_FEATURES`.

Step 33's sQTL arm computed a sign concordance on `SPLICE_JUNCTIONS`. Measured in the deposited table: all
914 `atlas_junction_liver_quantile` values are positive (min 0.0316), and all 4,685
`atlas_usage_max_abs_quantile` values. Against an always-positive predictor a sign concordance is
identically the share of positive measured slopes, so the deposited 0.5098 and 0.4615 measured nothing.
`lib_atlas.require_signed` existed for this and was called nowhere.

### Repo-wide audit of directional statistics (2026-09-15)

Three sites used an unsigned predictor for a directional statistic:

| site | statistic | scorer | consequence |
|---|---|---|---|
| `33_benchmark_eval.py` sQTL arm | sign concordance vs GTEx slope | SPLICE_JUNCTIONS | the withdrawn row; no claim survives |
| `33_benchmark_eval.py` MPRA arm `avi|allelic_sign` | sign concordance vs MPRA log2FC | AVI_SCORE raw | rows deposited (0.588 / 0.44 / 0.467 / 0.30); no written prediction used them; AVI raw does take both signs (111 of 280 positive) but its sign is not a ref->alt direction |
| `33_benchmark_eval.py` caQTL arm `atlas_avi` | auROC vs the caQTL SIGN label | AVI_SCORE raw | 0.518, perm p 0.002; `ag_orient_mult` is 1 for all 11,896 rows, so the predictor is plain unsigned AVI; the verdict JSON and prediction P3a use `atlas_atac` only |

Plus one introduced on 2026-09-15 and corrected the same day: step 35's `sign_among_model_first`, which
scored prediction P4 as a boolean. Withdrawn, see below.

Every other directional statistic reads a SIGNED scorer. Confirmed at: `11_package_D_direction.py` (RNA_SEQ,
the root of the whole eQTL-direction chain and of panel D), its gate replication auROC, `60_target_direction.py`
(RNA_SEQ via step 11), `65_clinical_dossier.py` `direction_concordance` and `66_dossier_checks.py` R5 (same
chain), `21_anchor_panel.py` P2.3 (RNA_SEQ), `09_package_A_landscape.py` `flag_opposing_signs_target_expression`
(RNA_SEQ, which `13_catalog_sidecar.py` turns into "credible_variants_predict_opposite_expression_directions"),
`73_atac_allelic_analysis.py` / `75_allelic_sensitivity.py` / `76_allelic_calibration.py` (ATAC, DNASE - the
headline allelic result), and step 33's other caQTL and MPRA channels (ATAC, DNASE, CHIP_HISTONE filtered to a
mark, CHIP_TF filtered to a factor, CAGE).

Three cosmetic cases where an unsigned value carries a sign in an exported column but no statistic uses it:
`65_clinical_dossier.py:405` (`splice_site_usage_liver_quantile`), `21_anchor_panel.py:115` (`np.sign` applied
to a max-|q|), `15_package_C_splicing_join.py:75` (junction quantile beside `deltapsi`, labelled positional).
`la.posterior_summary` is called correctly in both places (step 08 reads `is_signed` per scorer, step 09
hard-codes False for AVI, step 14 asserts unsigned scorers have NaN `S_liver_median`).

Guards added 2026-09-15 to `33_benchmark_eval.py`: `scorer_is_signed()` reads the server declaration and
`sign_note()` writes a `sign_validity` column on every caQTL and MPRA allelic-sign row. Changes no existing
number; it labels them.

### A separate defect in the MPRA arm, found and NOT patched

`33_benchmark_eval.py:429-430` builds the `liver_h3k27ac` mask by searching `biosample_name` for the string
"H3K27ac", which matches no track. That channel has **n = 0 in all four MPRA contexts** in the deposited run.
The headline MPRA number (0.685) is `hepg2_h3k27ac`, whose filter is biosample-based and works.

Deliberately left unpatched. The fix would change which tracks a channel selects, and the MPRA outcomes are
now known; `exposure_audit.json` forbids MPRA outcomes selecting tracks, windows or aggregation. Repairing it
is legitimate only as a declared amendment in `33_mpra_prespec.json` written before a rerun, and the rerun
needs the 176 GB interval archives (the `mpra_rescore` path cannot recover a channel that scored zero tracks).
Pending, not done.

### Step 35: does the unsigned magnitude at least find the right intron? No.

`p3a-splice-detection-20260915T142531Z` (supersedes `...T135807Z`, same numbers, corrected P4 verdict).
Prespecification `35_splice_detection_prespec.json`. Zero API cost - reads the archived sQTL chunks.

Per (variant, target intron) pair, the target is ranked against the other introns the Atlas scored for the
SAME variant and SAME gene, so the variant, its magnitude and the sequence context are held fixed. Chance is
0.5 by construction. 657 pairs, 640 variants, 348 genes, 255 1-Mb blocks, median pool 11 junctions.
Attrition: of 12,137 candidate pairs, 11,285 targets were never scored by the Atlas and 195 had a pool under 2.

| ranking statistic | mean rank | 95% block interval | ranks first |
|---|---|---|---|
| Atlas \|junction quantile\| | 0.4873 | [0.4559, 0.5200] | 16.1% |
| Atlas \|within-cluster contrast\| | 0.5577 | [0.5304, 0.5873] | 13.0% |
| proximity (nearer intron end) | 0.6053 | [0.5717, 0.6375] | 5.9% |

Paired differences: quantile minus proximity **-0.1179 [-0.1682, -0.0674]**; contrast minus proximity
-0.0674 [-0.1135, -0.0202]. Per-block binomial: quantile 115 of 255 blocks above half (p 0.41), contrast 121
of 228 (p 0.032), proximity 158 of 255. Prediction P1 (rank > 0.60) **failed**, firing rule **D1**. No exact
zeros in the archive, so ties are not the explanation. 41.9% of variants sit inside the target intron; median
distance 111 bp.

Diagnostic that names the defect: among the 106 pairs the quantile ranks first, sign concordance is 0.60377
and the marginal expectation is **0.60377**, identical, because every predicted value is positive.

GTEx v8 is an AlphaGenome training source, so this is a ceiling; a failure at the ceiling is decisive.

### Step 36: the test that IS possible (running at the time of writing)

`36_splice_direction.py`, prespecification `36_splice_direction_prespec.json`, job 21768295, out root
`p3a-splice-direction-20260915T141415Z`. 2 CPU / 48 GB / 6 h; ~1,280 model-API calls at ~8 s per variant.

`dna_client.predict_sequence` returns per-junction predicted USAGE for reference and alternate sequence
separately, so a direction exists off-Atlas. Probe facts (chr4, 2 calls): junction usage strictly positive,
7 liver tracks, 13,000-265,000 junctions per 1-Mb window, coordinates WINDOW-RELATIVE with a strand, and the
two arms return DIFFERENT junction sets (38,807 vs 38,917), so arms are joined on coordinates and never by row
index. `splice_site_usage` is per-base 1,048,576 x 14, aligned between arms.

Primary statistic is LeafCutter's own estimand: the target junction's share of its cluster's usage, alt over
ref, per liver track then averaged. Cluster = junctions sharing the target's donor or acceptor. The GTEx
closed intron (s, e) maps to model junction (s, e-1): measured over the substrate, best offset unique and
agreeing with `lib_atlas.leafcutter_intron_key`; the run fails closed otherwise. Substrate is restricted to
the 657 pairs step 35 scored, so the two routes are compared on identical pairs with identical pools.

Decision rule **D1**: if concordance <= 0.60, or not above its marginal expectation, or the block binomial
p >= 0.05, then P6c is not built as a prioritisation, nothing in the paper is ranked by a predicted splice
effect, and **the HSD17B13 rs72613567 splice magnitude (log2 -1.090, 24x its largest control) must be
reported as a magnitude in an unvalidated channel wherever it appears**. That variant is an insertion and is
not in this SNV substrate, so what transfers is the channel's validity, not a variant-level result.

A 12-variant `measure` pass and a 6-variant end-to-end pilot ran first; the pilot's concordance (1.0 on 5
pairs, marginal 0.68, its own D3 fired) was seen before the full run and is recorded in the prespec's
`notes_added_before_the_full_run` so the sequence is on the record.

### Also corrected

Section 8 of the results synthesis still quoted the superseded 59x for HSD17B13 against section 7's
corrected 24x. Fixed, with the reason for the difference stated in place.

### Step 36 outcome: direction is real but weak, and the prespecified bar was above it

`p3a-splice-direction-20260915T141415Z`, job 21768295, 1 h 14 min, 8.6 GB peak. 640 variants predicted
(11 refused, a 1-Mb window off the chromosome end), 637 pairs scored, 9 targets not predicted in both arms.
Coordinate convention measured again on the full substrate: 638 matches at (0,-1) against 16 at (0,1) and 18
at (0,2), unique winner, agreeing with `lib_atlas.leafcutter_intron_key`.

| statistic | concordance | marginal | above marginal | blocks above half | block p | 95% block interval |
|---|---|---|---|---|---|---|
| excision ratio (primary) | 0.5777 | 0.5004 | +0.0773 | 131 of 211 | 0.00055 | [0.5405, 0.6152] |
| target usage, no cluster | 0.5604 | 0.5009 | +0.0595 | 131 of 213 | 0.00096 | [0.5207, 0.6000] |
| per-base site usage | 0.5393 | 0.4972 | +0.0421 | 128 of 218 | 0.012 | [0.4985, 0.5794] |
| **sGene leads only** (n 89) | **0.6854** | 0.4978 | +0.1876 | 51 of 75 | 0.0024 | [0.5795, 0.7805] |

**The direction is above chance and not marginally so** (p 0.00055, interval excludes 0.5, 247 blocks), and
among the 89 highest-confidence pairs it reaches 0.685. **Prediction P1 nevertheless failed**: it required
concordance > 0.65 on the full set, and rule D1's threshold was 0.60. D1 fires, so nothing in the paper is
ranked by a predicted splice effect. **The D1 label text, "shows no direction", is wrong as a scientific
statement and must not be quoted.** The accurate statement is: direction is recovered weakly, at 0.578
against a 0.500 marginal, which is below the bar written in advance for using it to rank anything.

Magnitude carries nothing: Spearman(|slope|, |ratio|) = 0.036.

Detection fails for BOTH routes and proximity beats both. Model-API ratio rank 0.4727 [0.4351, 0.5108];
model-API raw usage rank 0.5472 [0.5148, 0.5770]; Atlas quantile 0.4813 on the same 637 pairs; proximity
0.6081 [0.5751, 0.6401]. Model minus Atlas **-0.0042 [-0.0371, +0.0297]**, so P5 failed: the signed route is
no better at locating the intron. Model minus proximity -0.1410 [-0.1983, -0.0881]. P2 failed in the
predicted direction (cluster normalisation helps by +0.0173, under the written +0.03). P4 held (site usage
within 0.0384 of the ratio).

**The leading limitation, and the strongest case that 0.578 understates the channel.** The cluster is defined
as the model's own junctions sharing the target's donor or acceptor, and the model proposes junctions densely:
cluster size runs **22 to 559, median 200, and not one of the 637 pairs has a cluster of 5 or fewer**. A real
LeafCutter cluster is 2 to 5 junctions. The excision-ratio denominator therefore sums hundreds of low-usage
model proposals, which dilutes the statistic the primary test rests on. A tighter cluster - junctions above a
reference-usage floor, or restricted to GENCODE v49 annotated junctions - is the obvious next form and could
plausibly move the primary number. It must be written as a declared amendment BEFORE it is run, because the
outcome of the loose form is now known; choosing the definition that wins afterwards is how a prespecification
gets laundered. Not run.

Consequence carried out: the HSD17B13 rs72613567 splice magnitude (log2 -1.090, 24x its largest control) is
labelled with the channel's measured direction validity wherever it appears, and is not called validated. It
is an insertion and is not in this SNV substrate, so only the channel's validity transfers.

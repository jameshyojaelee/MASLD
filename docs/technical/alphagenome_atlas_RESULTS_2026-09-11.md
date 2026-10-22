# AlphaGenome Atlas × MASLD Resource — results synthesis

Every number below is read from a deposit under `GWAS/finemapping/results/alphagenome_atlas/`, with the
control that licenses it and the written prediction it was tested against. Execution history, defects and
corrections are in `alphagenome_atlas_execution_2026-09-09.md`.

## 1. What the Atlas can and cannot do in liver

Five independent substrates, each with a prespecified prediction and a matched null.

| channel | substrate | n | result | prediction |
|---|---|---|---|---|
| ATAC/DNase | Currin bulk-liver caQTL, 138 donors | 10,966 leads | auROC **0.7128** [0.703, 0.723] | held |
| ATAC | GSE281367 snATAC allelic, **same variant** | 430 sites / 48 loci | concordance **0.6605**, p 1.7e-06 | exceeded |
| ATAC | GSE244832 allelic replication | 277 sites / 44 loci | 0.5668, p 0.048 | consistent, underpowered |
| H3K27ac/DNase | Hu MPRA HepG2 control | 3,498 oligos | auROC **0.685** | held (0.60-0.70) |
| RNA gene-level | 211-set liver eQTL direction | 201 leads | 0.5304 [0.449, 0.607] | held ("at chance") |
| splice (Atlas) | GTEx liver sQTL | 914 pairs | **no direction test was possible**: the scorer is unsigned | withdrawn |
| splice (model API, signed) | GTEx liver sQTL, LeafCutter estimand | 637 pairs / 247 blocks | concordance **0.5777** vs marginal 0.5004, p 5.5e-04 | **failed** (predicted > 0.65) |

**Accessibility carries real liver signal on three assays. Expression does not recover measured liver eQTL
direction.** Splicing is a separate case: the Atlas emits no signed splice score at all, so that row was
never a test; the model API does, and it recovers direction weakly (0.578, below the 0.65 bar written in
advance). Both are worked through in the subsection below. AVI, the aggregate variant-impact score, is the WEAKEST MPRA channel in all four
contexts (0.535-0.611) — and it is also unsigned, so it is used only as a magnitude.

### The caQTL arm reproduces the archive exactly
Re-implementation vs archived AlphaGenome-API row: **delta 0.0**, n 11,896 both sides. On the 10,966 variants
both models scored, Atlas 0.71278 vs API 0.71454, **delta 0.00175** (predicted < 0.02). Spearman between the
two score sets 0.981. The Atlas scorer and the API are the same model exposed two ways, and they agree.

### The splice arm tested a direction the Atlas never emitted (corrected 2026-09-15)

**The Atlas serves no signed splice score.** `scorer_metadata.tsv` records `is_signed=False` for
`SPLICE_JUNCTIONS`, `SPLICE_SITES` and `SPLICE_SITE_USAGE`, and the deposited table confirms it: all 914
`atlas_junction_liver_quantile` values are positive (min 0.0316), as are all 4,685 splice-site-usage values.
A sign concordance against an always-positive predictor equals the share of positive measured slopes, so the
deposited 0.5098 (and 0.4615 among sGene leads) measured nothing about direction. The one splice statistic
with a real sign, `atlas_leafcutter_contrast` (48.0% negative, 0.5380), signs *the target junction's
magnitude moved more than its siblings'*, which is not the direction of an effect either.
`lib_atlas.require_signed` exists for this check and was called nowhere.

**What replaces the withdrawn claim.** Not "the model predicts splicing badly" but "this API exposes no
splice direction". Two tests follow, both prespecified after the defect was found and before any statistic
was read.

**Step 35: does the unsigned magnitude at least find the right intron?** It does not.
(`p3a-splice-detection-20260915T142531Z`, which supersedes `...T135807Z` with identical numbers; prespecification `35_splice_detection_prespec.json`.) For each of
657 sQTL pairs the target intron was ranked, by |junction quantile|, against the other introns the Atlas
scored for the **same variant and same gene** — the variant, its magnitude and the sequence context all held
fixed. Chance is 0.5 by construction.

| ranking statistic | mean percentile rank of the measured intron | 95% block interval | ranks it first |
|---|---|---|---|
| Atlas \|junction quantile\| | **0.4873** | [0.456, 0.520] | 16.1% |
| Atlas \|within-cluster contrast\| | 0.5577 | [0.530, 0.587] | 13.0% |
| proximity to the variant | **0.6053** | [0.572, 0.638] | 5.9% |

657 pairs, 640 variants, 348 genes, 255 independent 1-Mb blocks, median pool 11 junctions. The quantile is
at chance and **below the proximity baseline: paired difference -0.118 [-0.168, -0.067]**, so "which intron
is the variant nearest" beats the model. The contrast form carries a little location signal (121 of 228
blocks above half, p 0.032) and is still below proximity, -0.067 [-0.114, -0.020]. Prediction P1 (rank > 0.60) **failed**, which fired
decision rule D1. GTEx v8 is an AlphaGenome training source, so this is a ceiling, and a failure at the
ceiling is decisive.

One diagnostic worth keeping: among the 106 pairs the quantile does rank first, sign concordance is 0.6038
and the marginal expectation is **0.6038** — identical to four decimals, because every predicted value is
positive. That equality is the signature of a vacuous direction test.

**Step 36: a signed splice effect exists off-Atlas, and it recovers direction weakly.**
(`p3a-splice-direction-20260915T141415Z`, prespecification `36_splice_direction_prespec.json`.) The model
API returns predicted junction usage for reference and alternate sequence separately, so a direction exists
there even though the Atlas exposes none. On the same 637 pairs, scoring LeafCutter's own estimand — the
target junction's share of its cluster's usage, alt over ref:

| statistic | concordance | marginal | blocks above half | block p |
|---|---|---|---|---|
| excision ratio (primary) | **0.5777** [0.5405, 0.6152] | 0.5004 | 131 of 211 | **5.5e-04** |
| target usage, no cluster normalisation | 0.5604 | 0.5009 | 131 of 213 | 9.6e-04 |
| per-base splice-site usage | 0.5393 | 0.4972 | 128 of 218 | 0.012 |
| **sGene leads only** (n 89) | **0.6854** [0.5795, 0.7805] | 0.4978 | 51 of 75 | 0.0024 |

**Direction is above chance and not marginally so**, and among the 89 highest-confidence pairs it reaches
0.685. Prediction P1 still **failed**: it required 0.65 on the full set, and decision rule D1's threshold was
0.60. **D1 therefore fires and nothing in the paper is ranked by a predicted splice effect** — but D1's own
label text, "shows no direction", is wrong as a scientific statement and is not quoted anywhere. The accurate
statement is that direction is recovered weakly, below the bar written in advance for acting on it. Magnitude
carries nothing (Spearman |slope| vs |ratio| = 0.036).

Detection fails for **both** routes, and proximity beats both: model-API ratio rank 0.4727 [0.4351, 0.5108],
Atlas quantile 0.4813 on the same pairs, proximity 0.6081. Model minus Atlas **-0.0042 [-0.0371, +0.0297]**,
so the signed route is no better at locating the intron (P5 failed).

**The limitation that most likely understates this.** The cluster is the model's own junctions sharing the
target's splice site, and the model proposes them densely: cluster size 22 to 559, median 200, with **no pair
below 6**, where a real LeafCutter cluster is 2 to 5 junctions. The ratio's denominator is therefore diluted
by hundreds of low-usage proposals. A usage floor or a GENCODE-annotated cluster is the obvious next form and
could move the primary number; it must be declared as an amendment before being run, since the loose form's
outcome is now known. Not run.

## 2. The new evidence: same-variant allelic concordance

The only test in this paper that puts an Atlas liver prediction against a measured allelic effect **at the
identical variant**, within donor, in adult MASLD liver.

- concordance **0.6605**; marginal-skew expectation **0.499**; 40 of 48 1-Mb blocks above half, exact
  binomial **p 1.7e-06**; block bootstrap 0.733 [0.608, 0.843]
- strong imbalance (|log2| > 1, n 22): **0.955**
- robust to every one of **18** het-calling filter combinations (0.634-0.688, block p 5.5e-06 to 1.7e-02)
- reference bias -0.007, and it points the WRONG WAY to explain the result: the aligner leans reference
  (measured positives 0.481) while the Atlas leans alt (0.539), so bias can only depress concordance
- independent cohort GSE244832 agrees in direction (0.5668) and strengthens with read depth

## 3. What is NOT established

| claim | status | why |
|---|---|---|
| the quantile is calibrated against measurement | **one cohort only** | GSE281367 trend +0.421 [0.182, 0.621] over 51 loci; GSE244832 +0.108 [-0.189, 0.383] |
| the channel prioritises the functional variant | **cannot conclude** | 57 sets sit at 11 loci; all intervals include 0.5; posterior is not distinguishable from Atlas |
| disease state modifies allelic effects | **untestable at this n** | per-site testing bottoms out at p 0.1 with 6v6 donors; global null clean (p 0.80 / 0.46) |
| variant combinations are non-additive | **no evidence** | 1,572 pairs, 0 BH signals in four channels; observed residuals SMALLER than the matched null |
| drug-target direction concordance | **3 targets, no rate** | THRB, RORA, ACACB all concordant; ceiling is 20 targets, rate withheld below n=10 |

## 4. Two constraints on presenting these scores

**A calibrated quantile is a locus-level statement, not a variant-level one.** The clinically validated
anchors sit at |q| 0.90-0.999 while ranking only 0.59-0.96 among the ~6,000 SNVs in their own 2-kb window,
and |q| > 0.9 occurs in just 6-13% of ordinary credible-set variants. The within-window rank must be
reported beside the quantile wherever a variant-specific claim is made.

**An extreme quantile is not a large predicted effect.** The model API shows raw single-variant effects of
~0.1% of local signal. The quantile RANKS a small change against a distribution of even smaller ones.

## 5. Coverage delivered

417,155 variants scored across 989 signals (the prior in-house profile covered 126), 12,857 channel-profile
rows, 3,452,688 Package D variant rows, 40,700 variant-peak rows, 7,245 splice event rows, 1,294 tissue and
66 ancestry rows, and a 989-row Gene Catalog sidecar. Track 0 validation **14/14**.

## 6. The coverage gap is a variant class, not an ancestry (P6b, 2026-09-13)

`GWAS/finemapping/results/alphagenome_atlas/p6b-ancestry-20260913T153316Z/`

P6b asked whether the annotation covers non-EUR credible sets as well as EUR ones, since non-EUR sets have
no matched liver eQTL panel and a static sequence model is the only functional prioritisation available for
them. Posterior mass decomposes exactly into what was queried, what fell below the query weight floor, and
what could not be mapped, and each of the eight core channels (ATAC, DNase, H3K27ac and the other histone
marks, RNA, TF, CAGE, contact maps, AVI) served **100% of queried mass in every ancestry**. At the variant
level the model is ancestry-neutral, which is what the prediction said.

The prediction that non-EUR sets would be less covered is **wrong, and wrong in direction**: matched on
fine-mapping resolution, EAS coverage runs +0.0044 [+0.0001, +0.0045] above EUR (98 signals, 56 blocks),
AMR +0.0094 [-0.0000, +0.0098] (21 signals, 15 blocks) and AFR +0.0000 [-0.0000, +0.0045] (19 signals, 13
blocks). SAS (4 signals, 4 blocks) and every non-EUR Universe B stratum fall below the prespecified block
minimum and are listed without a contrast.

**What the median hid.** Coverage sits at a construction ceiling, because the query rule targets cumulative
mass 0.999. The informative quantity is the lower tail, and it is entirely EUR: 51% of EUR liver-enzyme
signals fall below 0.999 while **no** non-EUR signal in any arm does. The mass is lost to one thing.

**The Atlas point-query API refuses indels (UNIMPLEMENTED), and that is where the posterior goes.**
Across the portfolio, 207,254 variant-signal pairs carrying 49.5 units of posterior mass are deferred
indels, against 158 pairs carrying 0.000 units that fail liftover. Per signal:

| Atlas coverage state | signals | what it means |
|---|---|---|
| covered | 623 | ≥ 0.999 of posterior mass queried |
| partial_indel | 353 | mass lost, majority of the loss to unscorable indels |
| uncovered_indel | 13 | essentially all posterior mass sits on indels the API will not serve |

The 13 uncovered signals are **not** a random tail. Ten are direct MASLD traits (FinnGen NAFLD/NASH ×7,
MVP NAFLD ×3); 58 signals lose more than 20% of their mass and 217 lose more than 5%. The heaviest deferred
indels sit at the best-known MASLD loci: hg19 chr22:44385219 C>CT (PNPLA3/SAMM50 region, 2.00 mass units
over 2 signals), chr19:19.7-20.2 Mb (the TM6SF2/SUGP1/GATAD2A region, seven separate indels at 1.00 each),
chr15:40.33 Mb (1.30 and 1.01), and **HSD17B13 rs72613567 (chr4:88231392 T>TA), the best-established
protective splice variant in MASLD, which carries 0.032-0.049 of the posterior in several signals and
cannot be scored at all.**

**Consequence for reading Track 0.** An uncovered signal still has a profile row. Its posterior-weighted
columns correctly read ~1e-23 and are honest. Its TOP-VARIANT columns (`avi_top_quantile` and the rest) are
unweighted maxima over the 2,000 near-zero-weight SNVs that were queried instead, and nothing in the row
protects a reader who takes one as the signal's mechanism. `atlas_coverage_flags.tsv` (366 rows) names
every signal that is not `covered`; gate on it before quoting a per-signal channel number.

This belongs in Figure 7's "what models can and cannot learn" arm as a coverage statement, not an accuracy
one: the limit is not that the model predicts indels badly, it is that this API will not accept them.

Corrections made in place before these numbers were read: the third mass term was first reported as
liftover failure, inheriting the misleading name `unmapped_mass` from `signal_query_coverage`; and the
matched contrast was vacuous, having matched on `n_cs_variants`, which is 0 for every coloc signal because
a credible set is a SuSiE object. It was replaced with the participation ratio sum(w)^2/sum(w^2). Both are
recorded as amendments in `62_ancestry_portability_prespec.json`; the written predictions were not touched.
Tail and indel statistics are labelled POST_HOC in the summary and are descriptions, not prespecified tests.

## 7. The indels are recoverable through the model API (steps 63, 64, 64b, 68, 68b; 2026-09-14)

The Atlas cannot reach indels by **either** route. Point queries return gRPC UNIMPLEMENTED, and interval
saturation returns only the three substitutions per position: 0 multi-base alleles in 104,352 variants
across 40 archived regions, checked 2026-09-13. That closes a question left open on 2026-09-09.

The model API is a different service and accepts an arbitrary sequence, so the indel can be built into the
1-Mb window directly, length-compensated at the far end (~500 kb from the readout). Probe on the variant
that matters most:

**HSD17B13 rs72613567 (hg38 chr4:87310240 T>TA), predicted splice-site usage log2 −1.090** — a 2.1-fold
reduction — against **−0.0044, +0.0048, +0.0071 and +0.0457** for its four matched SNV controls at its
signals scored the same way, i.e. 24x the largest control in magnitude (59x in the per-id run, whose
controls came from one of the variant's three source ids only). Its orientation is also the
test that caught a wrong fix: a rule that orients an indel by its longer allele calls this variant a
deletion, and it is a known T>TA insertion.
The T insertion disrupts the exon-6 splice donor, and a halving of predicted donor usage is the direction
its established mechanism predicts. The variant carries 0.120 of the posterior across three direct-trait
signals (FinnGen NAFLD, FinnGen NASH, MVP NAFLD) and the Atlas scores it `not_served`.

**CORRECTED 2026-09-13.** The first run (`p6f-indel-rescue-20260913T154340Z`, marked SUPERSEDED)
scored 140 of 152 targets and recorded 12 as reference mismatches. That was a defect in this step, not a
property of the data: it tested only the SOURCE id's `ref` field against hg38, and a GWAS source id carries
effect/other alleles in no particular order. **All 12 resolve with the alleles the other way round**; the
heaviest two rank 20th and 21st of the 152 by posterior mass (0.748 and 0.710) and the twelve together carry
2.79 of the list's 39.97 mass units. The earlier statement that the refusals were genuine liftover mismatches
was wrong — the check made was for coordinate
shifts, never for allele order. Both the target alleles and the SNV controls' alleles now come from the
reference-validated orientation (`lib_atlas.resolve_indel_alleles`: the source's allele order first, the
other order as fallback; mapped SNVs read their `hg38_ref`/`hg38_alt` straight from the crosswalk). The
numbers below are from the corrected run.

**A limitation the correction exposed, and dbSNP then closed it.** A VCF indel's short allele is a prefix of
its long one, so whenever the long allele is at the reference BOTH readings are representable: `C/CGT` at a
position reading `CGT` is equally a deletion `CGT>C` and an insertion `C>CGT`, giving different alternate
sequences. The reference cannot orient 100,468 of the 114,389 deferred indel ids (87.8%). These ids name
96,143 distinct variants: no hg19 position carries two different allele pairs, and 18,223 positions carry the
same indel under two or three ids written in opposite allele orders or with and without `chr`.

**Step 68 orients 112,863 of the 114,389 ids against dbSNP** (`p0-dbsnp-orientation-20260914T210732Z`).
The arbiter is build 157 on GRCh38.p14, with dbSNP records and our alleles both normalised against this
repository's single FASTA (the normaliser matches `bcftools norm`). dbSNP stores a repeat as one record whose
REF spans the run, so string matching had called most of these "absent"; after normalisation 13 are. In a
homopolymer dbSNP often holds BOTH reciprocal events (the run one base longer and one base shorter); there
the orientation whose alternate allele reaches a dbSNP FREQ of at least 0.01 is taken when it is the only one.
Build 156 on GRCh37 confirms by exact record, and two arms that resolve to different alleles are refused.

| step 68 outcome | ids |
|---|---|
| oriented by a single dbSNP record | 45,489 |
| oriented by the frequency tie-break | 67,374 |
| unresolved: dbSNP holds both events and both or neither are common, or has no record | 1,523 |
| refused: the two assemblies give different alleles | 3 |

**Step 68b checked this against the allele order of the GWAS files themselves** (`p0-orientation-check-20260914T230807Z`,
prespecification `68b_orientation_check_prespec.json`, all five predictions met, no decision rule fired).
Some GWAS files store alleles in a fixed order, which can be measured on their SNVs: UKBB_ALT/AST/GGT list
the alternate allele first (0 of 1,129,990 SNVs have allele 1 = reference), MVP and BBJ list the reference
first (99.3-100%), and FinnGen and most single-cohort files are mixed. Applied to the same files' indels, that
order matches the genome on 13,845 of 13,845 indels the reference can orient. On the reference-ambiguous
indels it gives:

| rule, reference-ambiguous ids with a file-order truth | n | correct |
|---|---|---|
| single dbSNP record | 31,439 | 99.86% |
| frequency tie-break | 66,621 | 99.85% |
| the source file's allele order alone (the previous fallback) | 99,548 | 47.36% |

Split by consortium, all four cells are at least 99.79%, and MVP and UKBB use opposite conventions, so a
truth set that one pipeline had oriented against dbSNP cannot explain the agreement. A seeded sample of 943
rows re-derived with `bcftools query` and `bcftools norm` and a separately written matcher reproduces step 68
on every row. Of the 3 refused rows, the study truth sides with the hg19 arm on 1 and the frequency tie-break
on 2.

**CORRECTED 2026-09-14: four earlier statements in this section were wrong.** (1) "dbSNP flips 38,438 of
the 62,960 it resolves, 61.05%" and the 50.06% "coin flip" came from the table before normalisation; on the
final table dbSNP flips the step-67 call on 53,014 of the 98,947 ambiguous ids it resolves (53.6%) and on 0
of the 13,916 the reference had already decided. (2) The explanation that "the effect allele is listed first
and is usually the non-reference one" was an inference I did not check; the measured cause is the mix of file
conventions above, and the 71.7% alternate-first rate describes that mix, not a rule. (3) "The assemblies
agree on 87,857 and disagree on 334" counted differing STATES; in the run before the tie-break, 0 rows had
two resolved arms with different alleles, and the final run has 3. (4) gnomAD was set aside because
frequency "would be a guess"; dbSNP's own FREQ field does the job, at 99.85%.

**One target per variant.** The per-id run scored 152 source ids that name 144 variants; at 3 sites the
two ids were scored as both reciprocal variants. Step 63 now merges ids naming one variant (same hg38 site,
same unordered allele pair), sums their mass and pools their signals, and step 64 refuses two rows for one
variant. With dbSNP applied, 38 of the 152 ids had been scored on the other orientation.

**The excluded variant class carries larger predicted effects than the included one at the same signals.**
`p6f-indel-rescue-20260914T231617Z`: **156 of 156 variants scored**, 526 SNV controls. The paired test
rests on the **142** with at least one SNV control at their own signals, in **97 independent 1-Mb blocks**;
the other 14 have no control because the indel IS the posterior there, which biases the comparison toward
the weaker cases. Absolute effects, block as the resampling unit, family = the five channels, BH within
stratum.

| channel | median paired Δ log2 | 95% block CI | blocks + | BH q |
|---|---|---|---|---|
| RNA | +0.0220 | [0.0138, 0.0276] | 80/97 | 1.5e-10 |
| splice | +0.0092 | [0.0051, 0.0167] | 80/97 | 1.5e-10 |
| H3K27ac | +0.0069 | [0.0037, 0.0129] | 74/97 | 3.4e-07 |
| ATAC | +0.0057 | [0.0034, 0.0107] | 73/97 | 8.0e-07 |
| DNase | +0.0036 | [0.0014, 0.0056] | 66/97 | 4.9e-04 |

**It survives the orientation strata, which are nested.** dbSNP orients 155 of the 156 (141 paired, 96
blocks): RNA +0.0219 q 2.5e-10, splice +0.0091 q 2.5e-10, ATAC +0.0060 q 3.9e-07, H3K27ac +0.0069 q
3.9e-07, DNase +0.0036 q 6.7e-04. Without the frequency tie-break, 60 remain (55 paired, 50 blocks) and all
five still clear q = 0.05 with larger medians: RNA +0.0305 q 1.2e-04, splice +0.0153 q 7.6e-04, H3K27ac
+0.0089 q 2.6e-03, DNase +0.0056 q 2.6e-03, ATAC +0.0038 q 1.6e-03. In the 19 the reference alone orients
(16 paired, 16 blocks), four clear q = 0.05 (RNA +0.0428, H3K27ac +0.0192, DNase +0.0073, ATAC +0.0072,
all q 0.027) and splice does not (+0.0328, 12 of 16 blocks positive, q 0.077).

Step 64b re-derived all four strata from `indel_rescue_effects.tsv` with pandas and scipy, blocks recomputed
from coordinates, and reproduces every count, median, sign-test p and BH q exactly; its own bootstrap
intervals move at most 0.29 of their width (Monte Carlo, different seed). The model API is deterministic:
the 375 controls and 108 variants scored on the same alleles in both runs are identical to the last digit.

Three limits travel with these numbers. The 156 targets are the variants carrying >= 0.05 posterior mass, not a
random sample of deferred indels. The median indel is short, so the splice result is how a sequence model
reads a frame shift at a donor motif, not a claim that indels matter more biologically. And these are
predicted track ratios from `dna_client.predict_sequence`, a different service on a different scale from the
Atlas quantiles; they are never pooled with them.

## 8. A per-variant clinical dossier, gated on coverage (P6d, steps 65-66, 2026-09-13)

`GWAS/finemapping/results/alphagenome_atlas/p6d-dossier-20260915T000218Z/` (rerun 2026-09-15 on the final orientations)
Prespecification `scripts/analysis/alphagenome_atlas/65_clinical_dossier_prespec.json`
(sha256 e720767483ac..., predictions written before any row was read; four amendments, all before results).

**1,422 rows over 202 direct MASLD signals and 927 distinct posterior entries.** One row per
(signal, variant) for every variant carrying >= 0.01 of its signal's posterior, plus each signal's
highest-weight variant whatever its weight, so a signal whose mass sits entirely on deferred indels still
appears. Each row joins the Atlas adult-liver channels, the chromatin measured in donors, the liver eQTL
direction, the ancestries whose credible sets contain the variant, the P2 anchor verdict where one exists,
the model-API rescue where the Atlas refused the variant, and the Gene Catalog next-experiment rule.

Adult liver only: 3 ATAC, 3 DNase, 2 H3K27ac, 5 RNA and 5 splice-site-usage tracks from the
`primary_liver` group. The CL:0000182 hepatocyte tracks are embryonic in this Atlas build, so these
numbers will not equal `liver_summaries.tsv.gz`, which pools them.

**951 of 1,422 rows are readable.** A row is readable only when the variant itself was served AND its
parent signal's posterior mass is covered. 461 rows sit at `partial_indel` signals and 10 at
`uncovered_indel` signals; 49 rows are deferred indels the Atlas does not serve, 48 of them oriented by
dbSNP and 32 carrying a model-API effect from step 63. The final orientations moved 4 of them from deletion
to insertion, where the earlier unnormalised dbSNP table had flipped them.

**All 12 recovered indels sit at liver-enzyme signals, not direct-MASLD ones**, so the correction above
strengthens the indel-vs-SNV comparison without moving any number in this dossier.

**A defect this lane found, now repaired.** A variant the Atlas refused has no hg38 `variant_uid`: step 04 leaves the
field empty and records the exclusion reason. Keying a posterior on `variant_uid` alone collapses every
deferred indel of a signal into one entry and sums their weights. No previously written table is affected
(step 08 already falls back to the source id; steps 09, 34 and 62 filter to mapped variants first), but the
dossier now keys on `variant_uid` when present and the normalised hg19 source id otherwise, and a guard
refuses two entries sharing an identity.

`lib_atlas.crosswalk_variant` now writes the uid, so any future run is correct, and **step 67 repairs the
run already on disk: all 114,389 deferred indels recover an hg38 identity** (141 rows keep none, having no
coordinate at all). All 49 unserved dossier rows now carry one. Recovering an identity does not make the
point-query API accept an indel: `atlas_served` stays False and every Atlas quantile stays empty.

### The five written predictions

| | prediction | observed | verdict |
|---|---|---|---|
| R1 | all 202 direct signals appear | 202 | met |
| R2 | 20-60% of rows inside a measured liver peak | **1.9%** | **not met, low by an order of magnitude** |
| R3 | every non-covered signal shows an unserved row | **34 of 37** | **not met** |
| R4 | >= 10 signals whose top variant the Atlas cannot serve | 14 | met |
| R5 | direction concordance in [0.45, 0.65] | 0.588 (n 204) | met |

**R2 failed because the interval was anchored on nothing.** The measured snATAC consensus peaks cover on
the order of 1% of the genome, so a fine-mapped variant landing in one is the exception. Marginally,
high-weight variants do overlap more often than low-weight ones at the same signals (20/886 = 2.26% vs
103/12,804 = 0.80%), but **that difference does not survive taking the 1-Mb locus as the unit: 5 of 14
informative blocks positive, median block effect 0.000, two-sided p 0.42.** The marginal ratio is
between-locus structure, not a within-locus preference. Post hoc either way; reported as an observation.

**R3 failed in an informative direction.** Three signals lose mass to deferred indels
(0.15%, 0.62%, 5.2%) without any deferred indel clearing the 0.01 floor or ranking first, so the row set
alone does not show the reader anything is missing. What gates a reader is the per-row
`signal_atlas_coverage_state` column, not the presence of an unserved row.

**R5 is inside the interval and is not a directional claim.** Observed concordance 0.588 against a
marginal-skew baseline of 0.506, consistent with the prior in-house AlphaGenome directional benchmark
(auROC 0.5607, n 211). The 204 rows sit in 19 independent 1-Mb blocks and one locus contributes 45 of them;
the block sign test is 13 of 18 positive, **two-sided p 0.096**, which does not reach significance.

**HSD17B13 rs72613567 was the anchor panel's only `not_served` verdict**, and step 63 now scores it
through the model API in its prespecified expected channel. The scale is different (a log2 track ratio, not
a calibrated quantile), so the P2 prediction of |quantile| > 0.8 is not tested by it; what can be said is
that the expected channel is the channel that moves, by **24x** the largest of its four matched SNV controls
(the 59x figure was the per-id run, whose controls came from one of the variant's three source ids only;
section 7 carries the final number).

## 9. Open

P4 saturation (103,390 regions, 71,000 retrieved at 2026-09-14 23:00 UTC, 614/h, ~53 h remaining) is the
only arm that can still change the picture: it asks whether regions the released chromatin predictor CAN
predict are sequence-sensitive or state-driven, and it would expand the allelic target list beyond 1,119
variants, adding the independent loci that both the calibration trend and the prioritisation test lack.
The user extended its walltime to 14 days on 2026-09-13, so it finishes inside the one job with about a
week of slack and the two chained resumption jobs are insurance rather than necessity. That extension also
settles the kill criterion the plan had written at 5 days: the full 96,460-region universe is kept.

**The P4 analysis is prespecified, and its scale was fixed before any test** (`43_saturation_prespec.json`).
The Atlas quantile is calibrated genome-wide, so inside a measured liver peak it saturates: in 6 sampled
regions the median |quantile| of a single substitution in liver ATAC was 0.48-0.93, while raw effects
spanned a 17-fold range. Step 43 therefore uses the raw effect divided by each track's mean level as the
primary scale and keeps the quantile beside it; its reduction reproduces by a separate code path to 1e-16.

Complete: steps 63-64 indel rescue with 64b re-derivation and 68-68b orientation with its check (section 7),
steps 65-66 clinical variant dossier (section 8).

**Pending work, as of 2026-09-15 15:00 UTC** (full reasoning in
`alphagenome_atlas_execution_2026-09-09.md`, section 2026-09-15):

| item | state | what it needs |
|---|---|---|
| step 43 P4 features + step 44 tests | blocked | U1 retrieval, the first 96,460 of 103,390 regions; at 80,000 and 614/h, about 27 h. Then `sbatch --cpus-per-task=16 --mem=96G --time=12:00:00 run_saturation_features.sbatch <new ts root> 16` |
| step 36 signed splice direction | **done** | 0.5777 vs marginal 0.5004, p 5.5e-04; D1 fired on the 0.60 bar |
| **P6c splicing vulnerability** | **not built as a prioritisation; D1 fired** | it becomes a coverage paragraph only: which Catalog genes have measured disease splicing and whether the channel covers their junctions. Nothing is ranked by a predicted splice effect. The HSD17B13 log2 -1.090 carries the channel's measured validity (0.578 direction, no detection) wherever it appears |
| tighter LeafCutter cluster for step 36 | optional, **not run** | the loose cluster (median 200 junctions, none under 6) dilutes the primary statistic; a usage floor or GENCODE-annotated cluster needs a declared amendment written BEFORE the run, since the loose form's outcome is known |
| P4 F4 (DA grammar) and F5 (program profiles) | blocked | U2 retrieval (the last 7,000 regions). Note U2c contributed 0 peaks, so the genetic-overlap arm is empty by construction |
| MPRA `liver_h3k27ac` channel is empty (n=0, filter searches the wrong field) | found, NOT patched | a declared amendment in `33_mpra_prespec.json` written BEFORE a rerun, plus a full re-read of the 176 GB interval archives; the rescore path cannot recover it |
| figure panels, exposure audits, authority updates | not started | after the arms above settle |

The deposited `avi|allelic_sign` MPRA rows and the `atlas_avi` caQTL row are not direction tests (unsigned
scorer) and now carry a `sign_validity` column explaining so on any rerun. No written prediction used them.

The thinnest number anywhere in this document is the count of independent 1-Mb blocks. The dossier's
direction check has 19; the indel-vs-SNV comparison has 91 (15 in its orientation-resolved stratum);
the post-hoc peak contrast has 14 informative ones and returns nothing. Every one of those would gain from P4's expanded target list, and none of them
should be strengthened by adding rows at loci already represented.

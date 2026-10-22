# Scientific review corrections: outcomes

The [scientific review of 2026-10-12](#scientific-review-2026-10-12) below records
findings and proposed documentation improvements. It does not implement the
recommendations or adopt new scientific results.
Its [2026-10-13 recheck](#review-of-this-review-2026-10-13) corrects the review
after independent adjudication of the 13 criticisms.
The subsequently authorized [coordinate repair and reruns](#coordinate-parser-repair-and-reruns-2026-10-13)
completed all four affected Figure 2 analyses on October 13. The user then
adopted those corrected results and requested panel regeneration. The historical
review's statements that these analyses were not run, and the initial rerun
record's candidate state, describe their earlier inspection states.

**Started:** 2026-09-23  
**Plan:** [SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23.md](SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23.md)  
**State:** dated candidate analyses and user decisions. The 774-unit bulk refit
was adopted on 2026-10-11, and r8 was selected as the paper's Gene Catalog that
day. The corrected Figure 2 coordinate-dependent results and regenerated panels
were adopted on 2026-10-13. Verified COLOC r2 is selected and the remaining
authorized corrections are recorded in the October 14 section below.
The complete publication bundle still requires the ROADMAP checklist.
Historical sections retain their dated input states; current states come from
[STATUS](../STATUS.md) and [RESULTS](../RESULTS.md).

Items 9–11 (cross-assay module candidate) are deferred by user decision. That
candidate stays out of the paper until it is corrected.

## Decisions taken during the work

| Decision | Choice | Why |
|---|---|---|
| GSE130970 healthy controls | 6 biopsies reconstructed from Hoang et al. 2019 Supplementary Table 1 (SRR9036375, 376, 377, 378, 379, 384); strict NAS-0 set of 4 as sensitivity | The supplement has no per-sample flag. Exactly one 6-subset of the 8 steatosis-0 biopsies matches Table 1 on sex (2 male), age mean (49.8) and SD (21); its complement reproduces the Table 1 NAFLD steatosis, lobular and fibrosis distributions exactly. Two of the six have lobular grade 1 (NAS 1) in GEO, while Table 1 lists all controls at NAS 0. One of the two reported males, SRR9036376, is female by expression (the model uses expression-inferred sex), so its RNA may not come from that control; the arm without the three sex-discordant GSE130970 biopsies bounds the effect (1,770 against 1,725 DEGs). |
| Figure 2 coding-mass inference | Intervals only, no p-value; corrected-coordinate results adopted by the user October 13 | The September 24-direct/19-shared census used truncated locus positions. Member-supported restoration gives 19 direct and 14 shared clusters. The 79% denominator statement and old intervals must not be reused; unequal GWAS ascertainment still prevents an automatic trait-label exchange. |
| Figure 7 items 12–14 | Relabel only; no retraining | Fold 1–4 selection results are exploratory; confirmation rests on the spent fold 0. |
| Catalog disease-state rule | Canonical padj < 0.05 and abs(log2FC) > 0.5; interval-null test kept as a sensitivity | The 2026-08-12 threshold change lists no Catalog exemption. |
| SuSiE-untestable genes | Separate state `genetic_untestable_shared_posterior`; kept in all denominators as not positive; exclusion reported only as a labelled sensitivity | Only adopted positives were rechecked for eligibility, so excluding untestable genes would clean one side of the comparison only. |
| Figure 2A display of untestable genes | Drawn as single-signal only, with a source-table state column and a caption sentence | ABF evidence stands; SuSiE gave no support and could not test them. |

## Items 3 and 4: Figure 2 coding posterior probability (coordinate recheck, 2026-10-13)

Historical outputs: `GWAS/finemapping/results/fig2_multiplicity/20260923T215748Z/` and
`GWAS/finemapping/results/fig2_robustness/20260923T215800Z/`. Seed 20260817,
10,000 bootstrap draws. An identical earlier run (`20260923T211611Z`) differs only
in column names. These outputs used truncated positions and their confidence
intervals and robustness claims are not valid evidence for the restored
coordinate analysis. At the initial coordinate recheck, new bootstrap and
robustness runs were **not run**; the subsequently authorized October 13
reruns and adoption below supersede that inspection state.

The locus field stores chromosome and position as decimal text. Missing trailing
zeros cause the parser to interpret GCKR locus 2.2773094 as chromosome 2,
position 2,773,094 instead of 27,730,940. Independent inspection of the gated
architecture and member tables identifies 715 truncated positions among 6,415
sets. For each set, exactly one trailing-zero restoration lies within 1 Mb of
every native member coordinate. This is member-supported restoration of the
sentinel encoding, not replacement by a member mean; 5,907 sentinels are not
literal members of their credible sets.

| Quantity | Historical September calculation, truncated coordinates | Member-supported position restoration, 2026-10-13 |
|---|---:|---:|
| Credible sets entering | 207 direct, 6,208 enzyme; 52 gate failures excluded | Same 6,415 gated sets |
| Protein-altering plus canonical-splice fraction, direct | 22.68% | 23.34% |
| Protein-altering plus canonical-splice fraction, enzyme | 3.39% | 3.87% |
| Ratio | 6.690613; old interval 3.27–12.99 must not be used | 6.026671; interval not recomputed |
| Difference | 0.192906; old interval 0.084–0.315 must not be used | 0.194695; interval not recomputed |
| Physical 1-Mb clusters | 264 total; 24 direct, 259 enzyme, 19 shared | 233 total; 19 direct, 228 enzyme, 14 shared |
| 62 robustness refits | Historical 5.13–9.76 range uses the faulty position parser | Not rerun |

Inputs are the [architecture table](../../GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_pip_architecture.tsv)
and [native member coordinates](../../GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv).
The direction of the descriptive coding contrast persists. Both parsing
implementations, derived cluster expectations, bootstrap intervals, noncoding
context families and 62 robustness refits require correction before the
quantitative Figure 2 inference can be used.

The following noncoding context estimates are historical and also require
coordinate-corrected refits. Their September values (difference, 95% interval): promoter
0.107 (0.006–0.240), ABC enhancer −0.009 (−0.015 to −0.004), lineage-accessible
0.046 (−0.015 to 0.118), resolved with no context 0.197 (0.050–0.356), unresolved
context −0.345 (−0.552 to −0.137). The historical ABC interval excludes zero, but direct
ABC mass is 0.1%, close to the boundary, so the percentile interval there is not
reliable. The ABC row is not described as a depletion.

## Item 6: coding enrichment as an explanation for sparse overlap (complete)

The claim that protein-altering variation cannot produce an expression QTL, and
that this explains the chance-level overlap, is removed from `docs/PAPER.md`,
`docs/RESULTS.md`, the abstract, introduction, Figure 2 and Figure 5/6 discussion
drafts. The text now states the descriptive coding-mass difference and gives any
role in the overlap as an untested hypothesis (Rivas et al. 2015, PMID 25954003).

## Items 7 and 8: Figure 5C (historical 844-record candidate)

The historical marker updates below do not apply to the adopted 774-unit fit.
Its primary HLA-DQA1 annotation, ENSG00000196735.14, has log2FC 0.379801 and
BH padj 0.024969, failing the canonical absolute log2FC >0.5 criterion. The
September filled-marker call came from ENSG00000206305.12 (log2FC 0.540029,
padj 0.026826). The current working 5C PDF actually has a hollow HLA-DQA1 RNA
marker; the prose legend's claim of a filled marker describes the historical
candidate. No 774-unit 5B/5C counterpart was found in the inspected candidate
outputs. The quantitative marker source must be refreshed before claiming
alignment to the current bulk fit.

- HLA-DQA1 mRNA marker used interval-null FDR 0.426 under the canonical rule. With
  conventional padj it is significant. On the retired 846-sample fit the marker is
  filled at log2FC 0.651 (padj 0.016); on the pre-correction F_five fit it is
  0.536 (padj 0.017). The final value comes from the corrected fit.
- Partial Spearman now uses one complete-case participant set per pair, both
  residualizations on the same design, df = n − 6. BH over all 125 pairs: 113 → 111
  significant. ABHD5–ballooning q 0.043 → 0.052; CD74–steatosis q 0.048 → 0.057.
  Independent reproduction of all 125 rho and p to 1e-15.
- The panel key now reads "padj < 0.05 and, for mRNA, abs(log2FC) > 0.5".

## Items 12–14: Figure 7 wording (complete)

- Item 12: fold 1–4 selected-recipe results relabelled exploratory; confirmation
  rests on fold 0 only (RESULTS, Figure 7 candidate captions, and the 7D caption
  producer).
- Item 13: "LD-tagging bound" renamed "restricted causal-lead oracle
  sensitivity". AlphaGenome 0.7441 is compared with 0.8774 on the same 3,627
  leads; 0.807 is reported separately for 32,322 leads.
- Item 14: "almost no headroom" removed. Both repeat-biopsy error scenarios are
  shown as sensitivity scenarios (58 and 30 pairs); conservative ratios 0.757,
  0.731 and 0.885.

## Item 1: GSE130970 controls and the bulk refit (fit complete; reproduction and downstream running)

Source: `RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20260923T231939Z/`
(`arm_deg_counts.tsv`, `design_rank.tsv`, per-arm `cohort_group_counts.tsv` and
`old_vs_new_*`). Model unchanged: `filterByExpr(~ 0 + group_binary)`, RLE,
voomWithQualityWeights, `~ dataset + inferred_sex + group_binary`, eBayes, interval-null testing
lfc 0.25; canonical DEG = padj < 0.05 and abs(log2FC) > 0.5. One biopsy per
participant. All three designs are full rank (7/7).

Before any corrected arm was fitted, a reproduction guard refitted the old labels
and matched the adopted fit exactly (23,370 genes; 1,347 canonical; 1,616 interval-null calls;
logFC, SE, padj and interval-null FDR within 1e-8; 26/26 checks).

| Arm | Participants | Control / disease | GSE130970 control / disease | Genes | Canonical DEGs (up / down) | Interval-null DEGs |
|---|---:|---:|---:|---:|---:|---:|
| Adopted 2026-08-10 | 844 | 157 / 687 | 23 / 53 | 23,370 | 1,347 (1,003 / 344) | 1,616 |
| Primary: 6 source controls | 844 | 140 / 704 | 6 / 70 | 23,702 | 1,725 (1,261 / 464) | 2,036 |
| Sensitivity: 4 controls with steatosis 0 and NAS 0 | 844 | 138 / 706 | 4 / 72 | 23,733 | 1,775 (1,305 / 470) | 2,054 |
| Sensitivity: 3 sex-discordant GSE130970 biopsies removed | 841 | 139 / 702 | 5 / 68 | 23,699 | 1,770 (1,282 / 488) | 2,068 |

Primary arm against the adopted fit: of the canonical DEGs, 1,326 are retained,
366 gained among the 23,370 shared genes, 33 gained among the 332 genes newly
passing the expression filter, and 21 lost. Across the 23,370 shared genes, log2FC
Spearman 0.983; median absolute change 0.024; maximum 0.498. Relabelling 17
NAFLD biopsies out of the control group increases the number of DEGs; the three
arms agree within 50 DEGs, and their DEG sets overlap closely: Jaccard 0.936 between
the primary and strict arms (33 and 83 genes unique to each) and 0.958 between the
primary and sex-discordance arms (15 and 60). The choice between the 4- and 6-control
definitions and the sex-discordant biopsies therefore changes few calls.

Unit: gene counts in this section are Ensembl gene rows, as in the adopted release.
The bulk universe includes genes on alternate contigs and patches (426 of 23,702
rows). On primary chromosomes the canonical DEG count is 1,700 (adopted: 1,325 of
1,347).

Worked-example genes (log2FC, padj; adopted → primary): GNMT −0.813 (6.7e-22) →
−0.896 (1.1e-22); MAT1A −0.433 → −0.463; CYP2C19 −1.438 → −1.694; HLA-DQA1 0.536
(0.017) → 0.540 (0.027); IGFBP7 0.390 → 0.364; BICC1 0.586 → 0.586. Directions are
unchanged.

The first attempt (`...-20260923T213011Z`) is failed provenance: same guard and
primary-arm values, then a summary-table bug stopped it.

## Items 2 and 5: COLOC eligibility and trait drivers (release built; figures running)

Source: `GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/` (`r2_summary.json`,
`qc/`). The patched `06_susie_coloc.R` calls `coloc.susie` on the full SuSiE fits,
so coloc 5.2.3's 50% shared-posterior check applies, and writes one row per signal
pair. Only the 1,307 adopted SuSiE-positive (study, gene) rows were rerun (279
study-by-chromosome tasks, 5 batch jobs). Dropping pairs cannot raise a gene's best
PP4, so adopted negatives cannot become positive. ABF values reproduce the adopted
master exactly.

Checks: an offline recomputation from the saved eQTL fits and the ATAC replay
shared-SNP sets agrees with the rerun on every eligibility call for all 816 tier-1/2
gene-studies and 5,362 signal pairs (PP4 within 4.5e-5). The only differences are
two top-SNP labels for ADH1C in the AMR study (4:100276771 vs 4:100277558) with
identical eligibility and PP4; non-EUR fits are not bit-reproducible, and these are
recorded in `qc/offline_hit_label_differences.tsv`. Of 74 changed calls, the 49
tier-1/2 calls reproduce the adopted PP4 from the replay before the drop; 25 outside
tier 1/2 have no replay table. An independent recount from the master files gives
the same membership.

| Quantity | Adopted 2026-08-17 | Corrected |
|---|---:|---:|
| Tier-1/2 multi-signal (SuSiE PP4 > 0.5) genes | 462 | 436 (26 lost, 0 gained) |
| Tier-1/2 single-signal only | 551 | 560 |
| Tier-1/2 union | 1,013 | 996 |
| Tier-1/2 SuSiE untestable (not negative) | – | 7 (ADAM15, ASAP3, CKB, GSTM2, LILRB3, PSMC5, RITA1) |
| Full portfolio, all 50 studies | 745 | 708 |
| SuSiE split, enzyme / direct (SuSiE driver) | 448 / 14 | 423 / 13 |
| ABF split, enzyme / direct (ABF driver; item 5) | 758 / 65 of 823 (SuSiE driver, blank ID included) | 705 / 117 of 822 |
| Tier-1/2 gene-study pairs (Fig 5D) | 816 | 767 (direct 28, enzyme 739) |

The 26 genes that lose tier-1/2 SuSiE support: ADAM15, ARL17B, ASAP3, CDK2AP1,
CERS5, CKB, DRAM2, FLCN, GGT7, GSTM1, GSTM2, GSTT1, GTF2H5, JAZF1, LILRB3, NMRAL1,
NUDT2, PDHB, PPIL3, PSMC5, RITA1, RPP14, S100A13, SLC17A1, TMEM116, VKORC1L1.
RORA and FABP1 (the Figure 2 examples) keep PP4 0.983 and 0.998. GNMT, MAT1A,
CYP2C19, IGFBP7 and BICC1 have no SuSiE-positive row in any study, so item 2 does
not change them. For four retained gene-studies (PRELID1, LNPEP, EVI2A, OMG) the
best pair is now ineligible, so Figure 5D re-selects their primary pair among
eligible pairs.

Under coloc's rule the last SNP column carries the null prior, so two eQTL signals
with negligible uniform-prior shared mass stay eligible (IDI2, ZNF584). Their pair
PP4 is at most 0.005, so no call changes; both rules are reported in
`qc/eligibility_rule_comparison.tsv`.

### Genetic and disease-state overlap on the corrected inputs

Source: `RNA-seq/results/manuscript_release/candidates/evidence-classes-corrected-20260925T175215Z/`
(supersedes `...-20260924T053244Z`). Corrected bulk primary arm (canonical DEGs) and
corrected COLOC. Duplicate bulk rows per symbol are resolved to the primary-assembly
row, then by padj, and COLOC symbols are mapped through Ensembl IDs, the same rule
as the Catalog inputs; the two now agree on all 1,709 DEG symbols. Untestable genes
stay in the denominator as not genetically positive.

| Scope (PP4 > 0.5) | Jointly testable | Genetic | DEG | Overlap | Genetic not DEG | Fisher OR (nominal p) |
|---|---:|---:|---:|---:|---:|---:|
| All traits, SuSiE | 14,243 | 409 | 1,082 | 24 | 94.1% | 0.75 (0.22) |
| All traits, SuSiE or ABF | 14,243 | 930 | 1,082 | 56 | 94.0% | 0.77 (0.06) |
| Direct disease, SuSiE or ABF | 14,243 | 150 | 1,082 | 10 | 93.3% | 0.87 (0.88) |

The orthogonality gate passes (at least 85% genetic-not-DE in all three scopes).
The class-validation gate still fails, as before.

## Rebuilt figures and Catalog inputs on corrected bulk and COLOC

**Figure 2 A–I** (`figures/candidates/fig2-coloc-r2-20260924T053244Z/`; validation passes).
Values read from the rendered PDFs.

| Panel | Adopted (round 2, 2026-08-17) | Corrected |
|---|---|---|
| 2A cascade: multi-signal / single-signal only / union | 462 / 551 / 1,013 | 436 / 560 / 996; SuSiE PP4 bands >0.9, 0.8–0.9, 0.5–0.8 = 244 / 89 / 103 |
| 2A caption | – | 3 genes drawn as single-signal only could not be tested by multi-signal COLOC; not negatives |
| 2C consequence counts | 36 / 73 / 107 / 455 / 342 (1,013) | 35 / 70 / 104 / 450 / 337 (996) |
| 2E ancestry partition: EUR-only / shared / non-EUR sub-threshold / non-EUR GWS-unique | 352 / 70 / 26 / 14 | 340 / 64 / 23 / 9 (GWS ancestry 8 EAS, 1 AFR) |
| 2F cross-ancestry genes | 13 | 13 (same genes) |
| 2G/2H RORA, FABP1 | PP4 0.983, 0.998 | unchanged |

**Figure 3G** (`figures/candidates/fig3g-coloc-r2-20260924T053244Z/`; interface tables in
`RNA-seq/results/genetics_state_interface/20260924T053244Z_coloc_r2_corrected_controls/`).
The fragment arm is the corrected bulk primary arm.

| Quantity | Adopted | Corrected |
|---|---:|---:|
| Threshold-free SuSiE pairs; Spearman on raw log2FC (nominal p) | 3,617; 0.014 (0.40) | 3,622; 0.021 (0.20) |
| Expression-matched pairs; paired OR (p) | 428; 1.43 (0.26) | 406; 1.32 (0.36) |
| Interface overlap, SuSiE PP4 > 0.5: genetic / DEG / observed vs expected | – | 406 / 1,077 / 23 vs 30.8; OR 0.72 (nominal p 0.15) |

**Figure 5 v7, historical 844-record candidate** (`figures/candidates/fig5-completion-2026-09-24-v7/`; 81/81 checks pass).
5B, 5C and 5D are rebuilt; 5A, 5E and 5F are copied unchanged (the dependency sweep found no
COLOC or control-label input to them).
- 5B, mRNA effects from the corrected fit joined by Ensembl ID (gene universe locked
  2026-07-10, not reselected): 210 protein-significant (211), 183 concordant 87.1%
  (185, 87.7%), 27 discordant 12.9% (26, 12.3%), Spearman 0.44 (0.44), 65% same
  direction (64%). 37 measured prioritized genes have no row in the corrected gene set.
- 5C: three mRNA markers change from hollow to filled under the canonical rule: LAMA5 (log2FC 0.569), CD151 (0.557) and HLA-DQA1 (0.540, padj 0.027). Filled mRNA markers go from 6 to 9.
- 5B and 5D use different COLOC releases within one figure: the 5B gene universe was locked on 2026-07-10 from the July COLOC release, while 5D uses the corrected release. The 5B source table records this; the caption must say it.
- 5D: 767 gene-study pairs (28 direct, 739 enzyme) from 816 (29, 787).
- Cosmetic: 5D prints "n =  28" with a double space.

**Gene Catalog inputs** (`Analysis/Multimodal_Program_Projection/candidates/gen-corrected-20260925T175215Z/`;
unsigned selection at `RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-20260925T175215Z/passport_input_selection.tsv`;
supersedes the `20260924T053244Z` roots). Canonical disease-state rule. Moving from
the canonical to the interval-null rule changes the class of 464 genes: 365 neither →
disease-state-only, 12 genetic-only → convergent, 85 disease-state-only → neither,
and 2 convergent → genetic-only. The Catalog build needs
the user's signature on the selection. Open items: the GEN_CLOSURE text
("34 of 447") and the Myojin class test (v1 classes) are marked stale.

## Gene Catalog corrected release (2026-09-25; historical, superseded by r8)

The section below records the September decision and build. The user selected
r8 for the paper on 2026-10-11, superseding this 844-record Catalog. r8 contains
28,233 genes: 400 genetically anchored, 26 concordant, 1,129 established-state
associated, 13,127 unresolved and 13,551 untested. Its genetic classes depend
on unadopted COLOC r2. See [the current Catalog numbers](../RESULTS.md#masld-gene-catalog)
and the October completion record below.

On 2026-09-25 the user instructed adoption of the corrected Gene Catalog release. The
input selection was signed that day (decision register
`2026-09-25-gene-catalog-corrected-release-adoption`; selection sha256
`0bb1d627…`), and the production build is
`RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-20260925T175215Z/`.
Its inputs are the corrected bulk primary arm, the corrected COLOC release
(`susie_coloc_r2_20260924T053244Z`) and the canonical disease-state rule.

The first build failed the validator rule that an untestable record may not carry a
statistic: six SuSiE-untestable genetics records showed the best PP4 of another
study or the ABF value. The builder now leaves those records without an estimate
or study driver. The failed build is kept at
`...gene-catalog-corrected-20260925T175215Z_failed-build-21904900/`.

The rebuilt Catalog passes PASS-00 to PASS-05: 28,017 genes, 56,034 evidence rows,
117 programs unchanged (6,633 memberships), 574 audit-only identity quarantine rows
and no fatal rows, zero tested-negative rows, no combined score or rank.

| Gene Catalog primary class | Previous Catalog (2026-08-07 bundle) | Corrected |
|---|---:|---:|
| Genetically anchored | 398 | 383 |
| Concordant | 33 | 23 |
| Established-state associated | 1,214 | 1,054 |
| Unresolved | 12,998 | 12,700 |
| Untested | – | 13,857 |

The previous bundle predated the adopted 2026-08-17 COLOC release and used the interval-null
disease-state rule, so the change combines items 1, 2 and 5 with the canonical rule.
The six SuSiE-untestable records show "untestable / insufficient shared posterior"
and no statistic.

Remaining step: PASS-06, the recorded manual review of the five portal review
images and the terminal seal. The images were inspected (overview, MTARC1, HKDC1,
GLP1R, THRB): call states, testability, the interpretation and grain boundaries,
source dependence and falsifying outcomes are all visible. Figures 6B and 6D render
from this bundle once it is sealed.

## Class-difference test on corrected inputs (2026-09-25)

Source: `RNA-seq/results/class_difference_test/cdt-corrected-20260925T203546Z/`
(`VERDICT.md`). Same sealed prespecification (sha256 `f5213ae6…`); substrate is the
corrected COLOC release and the corrected-control F_five fit; class counts checked
against the corrected Figure 3G interface (14,183 / 406 / 1,077 / 23).

| Hepatocyte-specificity arm | 2026-09-07 (adopted inputs) | Corrected inputs |
|---|---:|---:|
| Matched pairs (share of genetic-only genes) | 242 (58%) | 272 (71%) |
| Mean AveExpr, genetic-only vs state-only | 4.05 vs −0.29 | 4.16 vs 0.47 |
| Balance after matching (SMD) | 0.005 | 0.003 |
| Mean paired difference | −0.203 | −0.435 |
| Paired sign-flip p | 0.064 | 2.0e-5 |
| Wilcoxon p | 0.088 | 3.7e-5 |
| Prespecified permutation p (BH q) | 0.0036 (0.0072) | 1.0e-4 (2.0e-4), at the floor |
| Prespecified verdict | not separable (direction) | not separable (direction) |

The protein-level arm stays null (71 pairs, permutation p 0.14, sign-flip p 0.18).
The prespecified prediction (colocalized genes more hepatocyte-specific) is not
supported. The classes do differ, in the opposite direction: colocalized genes are
less hepatocyte-specific than expression-matched disease-state genes. This direction
was not prespecified and is reported as a post hoc result. The paired sign-flip null
is formalized as `04_signflip_diagnostic.py` and reproduces the 2026-09-07
hand-recorded values on the old pairs.

The disease-versus-control outputs of `20_pathway_continuum.R` (old labels) are not
read by any script or quoted in any paper document, so no reported number depends
on them.

## Figure 3 bulk panels, F0 arms, Figure 4E and worked loci on the corrected fit

Sources: `figures/candidates/bulk-release-corrected-controls-20260925T175215Z/`
(stage extensions, NAS grid, Fig 3A–E, S3), `figures/candidates/f0-reference-arms-corrected-controls-20260923T231939Z/`,
`figures/candidates/fig4e-corrected-controls-20260925T175215Z/` and
`figures/candidates/worked-loci-corrected-20260925T175215Z/` (validated). The first
bulk-release attempt (`...-20260923T231939Z`) stopped on a leftover 1,347-DEG check
in `render_fig3c_robustness.R` and is failed provenance.

| Quantity | Adopted | Corrected |
|---|---:|---:|
| Fig 3C genes unique to one study | 9,777 | 9,721 |
| Fig 3D endpoint sensitivity: participants (strict control / advanced) | 418 (115 / 303) | 416 (113 / 303); GSE130970 keeps 6 strict controls |
| Endpoint recovery of canonical DEGs | 1,307 / 1,347 (97.0%) | 1,630 / 1,725 (94.5%) |
| Endpoint all-gene log2FC Spearman | 0.938 | 0.940 |

F0-referenced stage contrasts (limma, `~ dataset + inferred_sex [+ control_status]`,
BH over the complete gene family; 3 cohorts with recorded F0; F0 reference = 104
biopsies, 14 of them source controls):

| F0 → F1 arm | Adopted | Corrected |
|---|---:|---:|
| A: all F0 | 1,471 | 1,426 |
| B: disease-only F0 (source controls removed) | 7 | 929 |
| C: all F0 with a control-status term | 258 | 975 |

The adopted explanation, that the F0-to-F1 signal is mostly control confounding
because it collapsed to 7 genes without controls, depended on labelling 23
GSE130970 F0 biopsies as controls. With source labels, about two thirds of the
F0-to-F1 signal remains among disease-only F0 references, so that explanation is
withdrawn. Later stage steps (F1→F2 544, F2→F3 518, F3→F4 234) are unchanged
across arms because their references contain at most one control.

Figure 4E: all eight BICC1/IGFBP7 program stage effects change by less than 0.001
(F4: 1.099 → 1.099; 0.761 → 0.762). Worked loci: GNMT is shown as ABF-supported
(PP.H4 0.706, 2021_34841290_NAFLD_EUR; SuSiE 0.128, primary rule not met).

## Continuum bulk transcript layer (Figure 4F–H inputs)

Source: `RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected-controls-20260925T192142Z/bulk/`
against the adopted `hac-molecular-layers-20260818T173348Z/bulk/`. The frozen axes
still use the same 139 signature genes. One further signature gene, DPEP1, now
passes the expression filter; it is not an axis component and is tested as an
ordinary gene. The first attempt stopped on an exact-universe check for the
signature and was changed to require only that the 139 frozen genes are present.

| Quantity | Adopted | Corrected |
|---|---:|---:|
| Genes tested | 23,370 | 23,702 |
| Continuum-associated genes | 14,651 | 14,818 |
| Canonical DEGs that are continuum-associated | 923 / 1,347 (68.5%) | 1,173 / 1,725 (68.0%) |
| Canonical DEG direction agreement with the signature axis | 0.832 | 0.808 |
| Held-out check: disease signature discovered in source-overlap cohorts (343 participants) | 1,146 genes | 3,669 genes |
| Held-out direction agreement in the evaluation cohorts (501 participants) | 0.873 | 0.739 |
| Held-out BH-supported discovery genes | 767 | 1,897 |
| Held-out all-gene Spearman | 0.426 | 0.323 |
| Stage vs continuum Spearman, F0→F1 / F1→F2 / F2→F3 / F3→F4 (signature axis) | 0.352 / 0.106 / 0.378 / 0.369 | 0.348 / 0.105 / 0.379 / 0.366 |

The discovery cohorts (GSE126848, GSE130970, GSE135251) include GSE130970, whose
control group changed from 23 to 6 biopsies. The corrected discovery signature is
three times larger and agrees less often in direction in the held-out cohorts,
although more discovery genes are supported there. Any text that quotes the adopted
held-out agreement (0.87) must be updated.

The 43-system participant-level analysis behind Figure 4G–H was checked on
2026-09-25 by tracing every input (reviewer trace, scripts 02, 03, 20, 31, 62, 66,
67, 70, 72, 73). No step that produces a Figure 4G/4H number reads a control label
or a disease-versus-control result: stage models use fibrosis stage and inferred sex
(`62:173-196`), inference uses dataset, sex, fibrosis stage and NAS
(`70:109-117`), and inference is restricted to GSE162694 and GSE213621 (501
participants; GSE130970 enters only the label-free descriptive five-cohort window
summaries). Its only link to the correction is the gene universe: it is built on the
adopted F_five 23,370-gene expression-filtered set, while the corrected fit keeps
23,702 genes (332 added, none lost). Rebuilding on 23,702 would re-derive the system
geometry (DPEP1 enters the signature, 721 of 5,323 pathway nodes and 40 of 43 systems
gain genes) rather than reproduce it, so the 33 / 29 / 22 / 18 maxT, rho 0.658 and
36/43 values stand, and the Figure 4 caption must state the 23,370-gene universe.
The disease-versus-control output of `20_pathway_continuum.R` (lines 160-176) used
the old GSE130970 labels; no Figure 4G/4H script reads it, but any other consumer of
it is out of date.

## Figure 3 lncRNA inset (rebuilt 2026-09-25)

The two lncRNA producers reran on the corrected bulk fit
(`RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-corrected-controls-20260923T231956Z/`);
the inset is `figures/candidates/fig3-lncrna-corrected-controls-20260925T195928Z/` (validator passes).

The inset's 1:1 propensity match of disease-associated lncRNAs to protein-coding
genes (0.2-SD caliper, without replacement) no longer met the 0.10 balance limit on
the corrected fit: 528 of 533 lncRNAs matched, with standardized differences of
−0.157 for gene length and 0.129 for expression variability. At the user's request
(2026-09-25) the match now tries calipers of 0.20, 0.15, 0.10 and 0.05 SD in order
and keeps the first that meets the unchanged limit. The choice uses the covariates
only; the direction outcome is computed afterwards. On the adopted inputs the search
keeps 0.20 SD and reproduces the adopted result exactly (368 pairs, OR 0.203); on the
corrected inputs it selects 0.10 SD (499 pairs; largest difference 0.099, gene
length). A Mahalanobis match within the 0.2-SD caliper was also tried and failed the
limit (0.169).

| Quantity | Adopted | Corrected |
|---|---:|---:|
| lncRNA share of tested genes / of canonical DEGs | – | 27.0% / 30.9% |
| High-confidence disease-associated lncRNAs | 279 | 358 |
| Unadjusted odds of being a DEG, lncRNA vs coding | – | 1.14 (1.02–1.27) |
| Adjusted for expression and detection | 0.439 | 0.507 (0.447–0.575) |
| Upregulated share: coding / lncRNA DEGs | 744/886 / 220/402 | 895/1,090 (82.1%) / 300/533 (56.3%) |
| Matched direction OR, lncRNA vs coding | 0.203 (0.136–0.302), 368 pairs | 0.210 (0.150–0.294), 499 pairs; discordant 41 vs 195 |

## Open items from the second review (2026-09-24)

- Resolved 2026-09-25: both builders now resolve duplicate bulk rows to the
  primary-assembly row first, and the evidence classes map COLOC symbols through
  Ensembl IDs. HLA-DRA is `neither` and SORD2P `disease_state_only` in both.
- Two overlap estimates exist for the same corrected inputs. The evidence classes
  use the bulk and atlas gene universe (14,243 jointly testable); the Figure 3G
  interface uses bulk genes present in the tier-1/2 COLOC table (14,183). They are
  different denominators for different panels, and the text must say which one it
  quotes.
- Resolved 2026-09-24: Figure 3G colours genes by the canonical rule
  (`figures/candidates/fig3g-coloc-r2-20260924T061535Z/`: 23 both, 180
  transcriptomics only, 383 genetics only); the S3G tiles use the historical
  interval-null significance label.
- `driving_trait` in the corrected `gene_level_coloc_tier12.csv` is the legacy
  column: the SuSiE driver when one exists. For an ABF-only gene (GNMT reads
  UKBB_ALT) the ABF driver is in `abf_driving_trait`. The Catalog builder uses the
  ABF column; any other consumer that reads `driving_trait` for ABF-only genes gets
  the wrong trait class. The atlas chain (`46d_convergence_evidence.R`) reads it and
  has not been rerun.

## Other findings

- Three GSE130970 biopsies (SRR9036338, 344, 376) and three GSE162694 biopsies
  (SRR13199861, 865, 917) have a GEO sex that disagrees with expression-inferred
  sex. The bulk model uses inferred sex. A sensitivity arm drops the three
  GSE130970 biopsies.
- The 2026-09-09 bulk release was never packaged (job 21610123 failed in
  `build_bulk_release.py`).
- Figure 5B plotted bulk logFC from the retired 846-sample fit; it is rebuilt from
  the corrected fit in the Figure 5 v7 candidate.

## Resource corrections implementation, 2026-10-05

The user requested implementation of the reviewed items 1–5 and a Markdown log
of all changes. This section is the implementation record; the five human
authorities remain the scientific and current-result entry points. Outputs are
candidate-only under `figures/candidates/resource-corrections-20261005-v1/`.
No synchronized release or named historical result is overwritten.

### Decisions and execution boundaries

- Exclude unresolved identity or relevant histology conflicts from primary
  inference; retain explicitly conditional original-record sensitivities.
  Neither the other-allele fraction cut nor expression-sex markers alone prove
  mixed samples. Source-linked donors and conservative genetic identity clusters
  must be distinguished.
- Require complete ALT/scaffold annotation correspondence before the primary
  bulk refit. Equivalent representations, distinct loci, alternative-only
  genes and duplications are separate dispositions. Irreducibly cross-gene
  fragments need not be assigned. The 1,833-DEG name-fold result is a sensitivity,
  not an accepted count correction or a final DEG target.
- Prepare both six-panel Figure 7 candidates for review. Preserve weights,
  model recipes and protected Currin fold 0; no additional selection or fitting.
- Keep the 117 programs, their weights, 43 system definitions, fixed external
  continuum loadings and k4/k6 bases unchanged. Significance and counts must be
  rederived when their inputs change.
- Work uses compute allocation `22037179` on `ne1dc4-043` (16 CPUs, 200 GB).
  No new SLURM job has been submitted for this correction work.

### Modified code and completed artifacts

| Change | Files or outputs | Scientific purpose and state |
|---|---|---|
| Qualified bulk inputs | `scripts/manuscript/resource_f_five/qualified_bulk_inputs.R`; `refit_corrected_controls.R`; `test_qualified_bulk_inputs.R` | Require complete annotation and counting evidence, bind the count matrix to its SHA-256, reject repeated analysis units and noninteger fragment counts, and carry the unit definition into model outputs. Preserve the unchanged historical reproduction guard. The corrected scientific fit has not run. |
| Corrected Figure 1C candidate | `scripts/figures/fig1c_corrected_genetics.R`; candidate `figure1/` | Reads the finished September genetics and any-eligible-trait support: 414 enzyme-only, 12 direct-disease/liver-fat-only and 10 both, totaling 436. No bulk effect enters this panel. The existing main PDF is unchanged. |
| Annotation correspondence | `scripts/analysis/alt_contig_recount/resolve_annotation_correspondence.py`; candidate `annotation_resolution/` | Enumerates the exact GTF and uses official GRCh38.p14 scaffold alignments. Initial analysis found unresolved correspondence; full refinement is in progress. No completed-count certificate has been issued. |
| Resource identity | Candidate `identity_resource/` and `scripts/analysis/fig34_corrections/` | Filtered fingerprints, source-resolution conflicts and analysis-specific rosters are in progress. Model A identities and protected partitions remain unchanged. |
| Figure 7 exposure and two review sets | `Analysis/MASLD_Model_Benchmark/config/artifacts/models/` crosswalks and exposure references; candidate `figure7_review/` | Dataset identity corrections and saved-prediction clean-study comparisons are in progress. Both Figure 7 families remain candidates; final numerical and visual checks are pending. |
| Catalog molecular-object records and missing Figure 6 panels | Existing Catalog v2 adapter with candidate-specific producers | Implementation is in progress. Final bulk-dependent classes wait for complete annotation/counting and identity results; no reviewer attestation has been fabricated. |

### Commands and checks executed

Commands below ran inside the existing compute allocation, using the existing
`/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/` environment unless stated.

| Command or inspection | Result |
|---|---|
| `Rscript --vanilla scripts/manuscript/resource_f_five/test_qualified_bulk_inputs.R` | PASS: nine supported-unit/count-certificate invariants, including rejection of incomplete annotation, incomplete counting, altered counts, repeated units, fractional contributions and samples outside the original model. No disease model fitted. |
| R parser checks for `refit_corrected_controls.R` and `qualified_bulk_inputs.R` | PASS. |
| `Rscript --vanilla scripts/figures/fig1c_corrected_genetics.R <September evidence classes> <September r2_summary.json> <candidate figure1 PDF>` | PASS: source partition 414 / 10 / 12, sum 436. Input hashes and sessionInfo saved beside the panel. |
| `pdftotext` and independent pypdf checks for candidate Figure 1C | PASS: one page, expected counts, no Type 3 fonts. |
| Figure 1C temporary PNG rendering with `pdftoppm` and visual inspection | PASS: count labels, trait labels and axis are legible and unclipped. Temporary raster proof is not a publication output. |
| `pdffonts` | Not run: command unavailable. The font subtype check above used installed pypdf instead. |

### Remaining dependencies

The primary bulk refit, changed stage/composition/continuum estimates, revised
maxT markings and final bulk-dependent Catalog classes have not run. They must
wait for a complete annotation/count certificate and the analysis-specific
primary unit roster. An incomplete correspondence table must not be replaced
with a source-dependent primary fallback under the user's selected requirement.
Final candidate checks, authority reconciliation and reviewer attestation are
recorded here when completed.

### Resumed work after the usage cutoff, 2026-10-05

The original local conversation was recovered from the other account's session
record, `01a109b4-d723-7bf1-920f-d4e7d9a33f5c`: the user had requested review of
items 1–5, selected the unresolved-record exclusion and full annotation
correspondence policies, then instructed implementation and a Markdown change
log. The three child records and saved artifacts showed that several analyses
had finished after the last parent update. The checks below supersede the
earlier “in progress” descriptions; they do not adopt a scientific release.
The user's subsequent direction is to present findings outside these five
items for their decision before fixing them.

Work continued inside existing allocation `21981199` on `ne1dc4-049`, without
submitting another job. The earlier allocation `22037179` remains running;
completed identity and Figure 7 computations were not restarted.

| Item | Recovered evidence and resumed change | What remains |
|---|---|---|
| Identity and inferential units | Independently reproduced all recorded source hashes, seven focused identity invariants, representative selections and the corrected source-control join. The primary screen selects 774 bulk/continuum genetic clusters, 392 endpoint records and 418 stage records; strict-depth counts are 738, 385 and 407. Receipts are under `identity_resource_qualified/resume_verification/`. | These are source-qualified genetic clusters, not verified independent people. Missing source labels cannot establish concordance. GSE193066 corroborates only 24/58 source repeat pairs in the primary screen, with four genetic pairs crossing source-person labels; its replication, ICC and derived ceilings remain excluded. |
| Qualified bulk inputs | `qualified_bulk_inputs.R` now checks both the actual source annotation summary and the selected roster against hashes recorded with the counts. Altered annotation evidence, a changed roster and an internally inconsistent completion claim are rejected. | Fragment counting and the scientific fit are not complete; no final DEG count, changed stage/composition/continuum estimate or maxT marking is available. |
| Counting annotation | The completed r2 enumeration contained 1,239 unresolved relationships. Resumed review found 16 supposedly ALT-only annotations whose mapped exons overlap same-strand primary gene bodies, without sufficient identity evidence; these remain unresolved, giving 1,255. Eleven records on four rDNA/NOR patches have no matching alignment files in the inspected official directory. | Annotation overlap alone does not establish equivalence or a distinct gene. Missing coordinate evidence and ambiguous gene identity must be resolved before primary counting/refitting. Shared symbols or primary-chromosome realignment cannot substitute for those relationships. |
| Catalog and Figure 6 | Candidate molecular-object records and local review views were recovered. The assembler's unit wording is corrected to source-linked participants or conservative genetic clusters; unsupported/repeated/empty units are rejected. Static gene-level COLOC support does not establish a signal-specific target link. | Existing class totals still depend on the historical September bulk fit; no final attestation, final class count or publication is supported. |
| Figure 1C | Independently rechecked input hashes and source partition 414 enzyme-only / 10 shared / 12 direct-disease/liver-fat-only = 436; read the PDF text and inspected a fresh raster proof. | The main artwork and synchronized result remain unchanged. |
| Both Figure 7 families | All 12 single-page PDFs were recovered and rechecked. `resumed_review_v2/` corrects two stale Figure 7F source-table verdicts, the clean-study exposure label and caption color wording. The rendered values already used the correct interval limits. Full study states and source paths are retained. | Comparisons remain conditional on saved models and eligible study identities. No new model selection, weights, protected Currin fold-0 outcomes or accuracy claims were introduced. |

Figure 7 verification independently reproduces the Guilliams identity join
(42 libraries, 19 donor labels and 315,007 cells), exposure counts of
54 clean / 22 confounded / eight unresolved encoder–head–study cells, and
24 saved clean-study comparisons, 20 with numerical support. The 10,000-draw
study-stratified donor-cluster bootstrap reproduces to a maximum absolute
difference of 3.78e-16. Two-model checks cover all six Figure 7C contrasts,
20 Figure 7D pooled MSE rows, 15 participant-pooled Figure 7E reductions,
the six source-table hashes, Figure 7A's 99-by-24 display and Figure 7B's
1,500 development rows restricted to folds 1–4. These checks establish
reproduction of the saved development evidence, not independent confirmation
of model performance or every source participant's identity.

The corrected Figure 7F table records UCE-33L MLP's interval
[-0.0884023, -0.0809825] as excluding zero and TranscriptFormer's MLP interval
[-0.00079267, 0.00830383] as crossing zero. The clean-study subset is explicit;
the original global exposure label is not reused for this comparison.

The revised annotation enumeration is in `annotation_resolution_r3/`, with
86,369 annotations: 78,691 primary genes, 5,980 equivalent representations,
152 duplicated representations, 266 distinct loci, 25 alternative-only
annotations and 1,255 unresolved relationships. `review_checks.json` records
the 16 changed dispositions and confirms that no grouped counting GTF or
primary-fit permission was emitted. Among the equivalent representations,
438 have partial exon-coordinate coverage but a unique HGNC/Havana identity
at the mapped primary locus. Gene-object correspondence does not imply that
every ALT-specific exon base has a verified coordinate projection. The producer
must restrict the partial-equivalent ALT features used for fragment rescue to
mapped assembly blocks; omitted insertion/tail sequence cannot support rescued
counts. This is the original conservative counting requirement, not a new
permission to count unverified exon sequence.

The clipping rule is implemented in
`scripts/analysis/alt_contig_recount/resolve_annotation_correspondence.py` and
tested by an insertion fixture in `test_annotation_correspondence.py`.
`annotation_resolution_r3/policy_receipt.json` and
`count_checks/exon_rescue_policy_checks.json` distinguish the implemented
restriction from actual fragment counting, which has not run.

The completed provisional Catalog review is `catalog_review_r4/`, preserving
all 28,017 source gene memberships and evidence classes in the molecular-object
sidecar. Every record retains `no_authoritative_link` with target-link basis
`none`; gene-level COLOC support is not converted into an inferred target.
Figure 6A/C/E are saved under its `panels/` directory and describe the Catalog
structure, existing worked entries and the actual local access interface.
The Figure 6E label says “View worked gene entries”; this candidate does not
claim an implemented search capability. The producers are
`scripts/portal/build_resource_catalog_sidecar.py` and
`scripts/portal/render_resource_catalog_review.py`. Source classes remain
provisional until the corrected bulk results exist.

| Executed check | Result |
|---|---|
| `Rscript --vanilla scripts/manuscript/resource_f_five/test_qualified_bulk_inputs.R` using the existing `rnaseq` environment | PASS, 12 invariants after adding the annotation/roster checks; no scientific fit. |
| R parser checks for `refit_corrected_controls.R` and `qualified_bulk_inputs.R` | PASS. |
| Identity focused tests and source/roster reconstruction | PASS, seven tests and exact recorded-source hashes; receipt under `identity_resource_qualified/resume_verification/`. |
| `verify_figure7_review.py` on `figure7_review/resumed_review_v2/` | PASS; source reconstruction and all 12 PDFs. |
| Figure 7 crosswalk/config tests | PASS, six focused tests using installed `tomli` compatibility for Python 3.10; no dependency added. |
| Annotation/fragment-counting Python invariants | PASS, 14 tests including exclusion of unaligned inserted exon sequence; no real recount. |
| `Rscript --vanilla scripts/analysis/alt_contig_recount/test_assemble_verified_counts.R` | PASS, valid mixed-unit assembly and rejection of unsupported, repeated and empty units. |
| Catalog membership/class check against its source bundle | PASS, all 28,017 records; no unsupported target-link attribution. |
| Figure 6A/C/E PDF text, font and visual checks | PASS, three single-page PDFs, no Type 3 font, legible and unclipped. |
| `pdftotext`, raster proofs and visual inspection | PASS for Figure 1C and both Figure 7 families; each Figure 7 PDF has one page and no Type 3 font. |
| `python3 scripts/manuscript/validate_resource_scope.py` after authority reconciliation | PASS. |
| `git diff --check` on the changed authority and bulk-input files | PASS. |

Figure 7 details and commands are saved in
`figure7_review/resumed_review_v2/RESUMED_REVIEW.md`, `verification.json` and
`logs/resumed_checks.json`. Named results, main figures and model weights were
preserved. The bulk fit, dependent Catalog classes and final scientific
adoption remain unfinished under the original full-correspondence requirement.

### Autonomous bulk-refit execution, 2026-10-05

The user approved proceeding with the detailed plan autonomously. New scientific
choices outside that plan still go to the user. Work continues within allocation
`21981199`; no new jobs, dependencies, named-release changes or publication were
introduced.

The exact execution rosters are in
`figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/bulk_execution_inputs/`.
Primary bulk/continuum uses 774 source-qualified genetic clusters, including
138 controls and 636 disease records. Endpoint and exact-stage rosters contain
392 and 418 records. Strict-depth counts are 738, 385 and 407 respectively.
The primary and strict bulk union requires 776 BAMs, including two strict-only
representatives. Missing source diagnosis or sex information remains qualified;
neither inferred sex nor a genetic cluster proves independent participant
identity. The close-relative screen is separately recorded and outcome-free.

`bulk_count_inputs/library_manifest_776.tsv` checks all required BAM paths,
layout/strand flags and source count headers. The 776 BAMs total
4,941,087,264,669 bytes (4.941 TB decimal). `pilot_libraries.tsv` selects the lower
median genetic-coverage library within each primary cohort, with SRR tie-breaking.
Selection does not use disease effects. The preparation and independent rank
verification are recorded alongside those tables. Full historical count-file
hashes were retained from their ledger, not rechecked during this header-only
preparation. Real diagnostic pilots use a deterministic 1/1,000 QNAME/RG sample
from a complete source-BAM scan, preserving every selected mate and alignment.
Their native unresolved GTF outputs cannot certify primary corrected counts.

The official Ensembl release 115 metadata identifies GENCODE 49 and
GRCh38.p14. Exact gene versions, strands and assembly aliases were checked
against the counting GTF before using the alternative-allele relationships.
Membership also requires independently mapped locus support, consistent HGNC
identity and copy checks; a group alone does not justify equivalence.
Authoritative patch and primary AGP components show that all 11 NOR gene bodies
lie in added sequence whose placements replace only primary chromosome `N`
short-arm gaps. These are alternative-only annotations, with family copies
retained separately.

The reviewed result is `annotation_resolution_r5/`, not the intermediate r4
directories. Its 86,369 annotations comprise 78,691 primary genes, 6,346
equivalent representations, 152 duplicated representations, 266 distinct loci,
36 alternative-only annotations and 878 unresolved records. Relative to r3,
367 unresolved records gain supported equivalence and 11 gain added-patch-only
dispositions. `ENSG00000282091.1` loses an earlier equivalence because its exact
exon correspondence and official allele-group target disagree. The 878 holds
include 870 unresolved coordinate/identity relationships, that conflict, and
seven annotations on KI270825.1. All 468 used official source files were
independently rehashed; the full annotation ID set is unchanged. No grouped
counting GTF or primary-fit permission was emitted.

The official GRC exclusions BED labels the entire KI270825.1 scaffold
contamination, despite the recovered exact official alignment and seven
strand-consistent allele correspondences. The user has been asked whether to
exclude those seven annotations as non-applicable for primary counting while
leaving fragments ambiguous between valid and excluded alignments unassigned.
This decision has not been applied. Other GRC false-duplication regions were
not automatically reclassified. Source URLs, hashes and gene-specific evidence
are in `bulk_execution_sources/`.

An independent exact-transcript metadata pass found no HGNC identifiers missing
from the unresolved gene rows. Shared RefSeq/Havana accessions are not used to
force identity: 21 shared-accession links have conflicting HGNC IDs and 19
mapped transcript links disagree with nonempty official allele-group targets.
These tables are supporting evidence, not additional automatic resolutions.

The fragment reader now rejects chimeric primary segments marked by `SA` or
STAR's `ch` tag, noncanonical HI sets, cross-contig pairs and nonreciprocal mate
coordinates. It does not add an improper-pair flag exclusion to the original
featureCounts `-B` rule. The installed featureCounts single/paired fixture
still accepts a complete within-gene multihit fragment once and rejects
incomplete or cross-gene alignments. The first real single-end pilot,
SRR8601578, reproduces all 35,458 native unique-read counts gene by gene. Its
native `-M` table contains 51,281 alignment assignments, whereas complete
one-fragment accounting accepts 35,764 fragments, including 306 same-gene
multihit fragments, and rejects 5,315 cross-gene groups. This is diagnostic
evidence on native annotations, not a corrected bulk estimate. Remaining pilots
are running with at most three workers (six CPUs/24 GB) in the same allocation.

Qualified stage and system inputs preserve source analysis units, accept dynamic
roster sizes and require corrected-expression lineage in projected scores.
Missing comparators or controls remain explicitly unestimable. Historical modes
retain their reproduction counts and frozen definitions. Corrected fixed-W NMF
scores and corrected program/continuum/system inputs do not yet exist, so real
stage/system inference has not run. An independent identity-only comparison
uses historical unique-only native counts and the 774/738 rosters; it remains
explicitly count-uncorrected and cannot substitute for the primary fit.

Executed checks: 37 focused Python tests pass; the installed single/paired
featureCounts fixture passes; annotation transition, source-hash, ID-universe
and held-contamination checks pass. Stage/system synthetic input checks and
launcher syntax checks pass. Primary fragment counting, the corrected bulk fit,
dependent figure estimates and final Catalog classes have not run. Remaining
pilot and identity-only results will be recorded here when available.

Further source review found 203 of the 240 incomplete/mixed exon-map records
with complete gene bodies inside an internal official GFF `Gap I` interval,
with matched sequence on both sides and no overlapping `M` segment elsewhere
in the scaffold alignment. Fifty retain existing identity/copy conflicts.
Of the other 153, 28 have declared alternative-allele groups without primary
members and 125 lack recorded group relationships. This proves local inserted
sequence, not global gene uniqueness or absence of other copies. Their
dispositions remain held; the unresolved total is still 878. Independent
coordinate consumption, anchors and 72 source hashes are recorded under
`bulk_execution_sources/internal_alt_insertions_r5_review/`.

NMF introduces a separate scientific decision. The original disease-only
preprocessing matrix was reconstructed from frozen inputs to a maximum absolute
difference of 2.66e-15; the cached nonnegative transform agrees exactly.
Recovered batch coefficients, gene minima and historical TMM factors are saved
under `bulk_execution_sources/fixed_nmf_preprocessing_review/`. Stage P1–P6
use the full fixed-W basis; the signature-excluded projection in the continuum
producer is a different prespecified use and must not replace those stage axes.

All 5,000 original feature versions occur in the GTF, but 233 ALT feature objects
merge into primary counting targets under r5. Ninety-seven duplicate target
groups contain 229 basis rows; 27 rows retain unresolved identities. Repeating
one corrected count across those old rows or collapsing the basis would change
the projection definition. The original preprocessing also normalized only
1,104 disease libraries; it contains no original TMM factors for controls.
No corrected projection was emitted, and no normalization or missing-feature
policy was invented. The user has been asked to choose between keeping NMF
held/historical and preparing a revised projection protocol for review. The
bases remain unchanged. Exact feature tables are in
`bulk_execution_sources/fixed_nmf_feature_units_r5_review/`.

All five native-library diagnostics are now complete. They scanned 546,274,850
source alignment records and selected 547,533 records with the same QNAME/RG
hash rule. Every native unique gene count matches the fragment reader exactly
in each cohort. Together they accept 174,469 unique reads/pairs plus 1,840
complete same-gene multihit groups and reject 17,828 cross-gene groups. The
native `-M` alignment tables sum to 242,274, compared with 176,309 accepted
read/fragment groups. Sampled NH/HI, supplementary/chimeric and mate defects
were not observed; this does not establish their absence in all libraries.
Earlier loaded helper versions were retained and their small assignment BAMs
rescored with the strengthened helper, without rescanning source BAMs.
`count_checks/native_pilots/five_cohort_diagnostic_summary.tsv` and its JSON
record this diagnostic-only evidence. No primary count certificate was emitted.

The primary 774-unit identity-only fit on unchanged historical unique-only
native counts is complete: 23,535 tested annotation rows; 1,756 canonical calls
(1,243 up, 513 down); 2,014 interval-null calls (1,364 up, 650 down). Canonical BH and
interval-null BH each cover their full expression-filtered family separately. Compared
with the corrected-control 844-library fit's 1,725 canonical calls, it retains
1,572, loses 153 and gains 184. LogFC Pearson correlation is 0.9881, median
absolute delta 0.0277 and maximum absolute delta 0.3919. This sensitivity jointly
changes roster/cohort/case composition, filtering, RLE and quality weights;
it cannot attribute the difference solely to repeat-library removal.

These are native annotation calls, not a corrected molecular-gene total.
Four targets, COL3A1, PRSS51, IGHA2 and FCGBP, each have two canonical source
rows under r5's equivalence map, and three canonical rows have unresolved
identities. Do not deduplicate that result table into a supposed corrected DEG
count: new fragment counts, filtering and model fitting are still required.
The actual inputs, exact original extracted statistical functions, full-rank
seven-column design, finite positive quality weights and independently
rederived point-null and interval-null BH values are verified under
`identity_only_comparison_r2/primary_774/`. The strict 738-unit fit subsequently
completed; the final checks and results are recorded below.

The same official GRC BED identifies 12 false-duplication intervals intersecting
57 gene bodies. Fifty-six entire gene bodies and exon unions lie inside those
intervals; TRPM3 has 1.96% intronic body overlap and no exon overlap, so its
expression features are not an exclusion case. The 56 complete overlaps
include CCNYL1B and hsa-mir-8069-1. They are primary-assembly annotations,
not ALT loci. No Figure 1C 436-gene members, worked Catalog examples or the five
paper hero genes intersect. Fifteen overlapping annotations are tested in the
774-unit historical-count diagnostic, but none has a canonical call. Excluding
reference artifacts would still change the eligible gene family and may alter
BH or normalization after recounting. The user has been asked whether to apply
the conservative reference exclusion to these 56 annotations; no disposition,
count, score or figure was changed. Gene/exon coordinates, source aliases,
hashes and existing-evidence joins are in
`bulk_execution_sources/false_duplication_r5_review/`.

Both identity-only fits are now complete, with no statistical model rerun for
verification. The strict 738-unit fit tests 23,940 native annotation rows and
gives 1,807 canonical calls (1,254 up, 553 down), plus 2,009 interval-null calls
(1,337 up, 672 down). Against the historical 844-library fit, it retains 1,539
canonical calls, loses 186 and gains 268, with no sign flip among retained
calls. LogFC Pearson correlation is 0.97565; median absolute delta is 0.03911
and maximum absolute delta is 0.67198. Among the strict canonical calls,
1,786 are primary annotations, 18 are mapped equivalents and three remain
unresolved. Five counting targets each have two called annotation rows.
Neither the 1,807 nor the 1,756 total is a corrected physical-gene count.

The root agent independently checked all input file hashes, exact source
roster/unit/group/sex joins, saved model-matrix values and rank seven, gene and
sample alignment, and every saved filtered count against the historical raw
matrix. Quality weights are finite and positive: primary range
0.07868329–2.85561784; strict range 0.07882073–2.76520049. A separate sorted-p
BH calculation reproduces both point-null and interval-null adjusted values to a maximum
absolute error below 4e-15. Called-row counts, directions, historical
membership changes, filtering changes and effect-delta summaries also agree.
Commands used the existing 21981199 compute allocation with one CPU and 8 GB;
no new allocation, model refit or dependency installation was needed. The
check script and final evidence are
`identity_only_comparison_r2/root_verify_results.R` and
`identity_only_comparison_r2/root_independent_verification.json`.

The corrected primary recount and fit have not run. All 878 unresolved
annotation dispositions still need scientific resolution; approval of the
seven contamination exclusions alone would not resolve the other 871.
Decisions on those exclusions, the additional 56 false-duplication exclusions,
and the fixed NMF projection remain unanswered. No user choice has been
inferred from a preselected answer or elapsed time. Existing named results,
fixed bases, primary counts and final Catalog classes were not replaced.

### Authorized reference exclusions and revised projection, 2026-10-05

The user's subsequent instruction, "fix all those", authorizes the two
documented reference exclusions and preparation of the revised NMF protocol.
The earlier unanswered-decision state above is historical. No new permission
was requested for those actions.

`annotation_resolution_r6/` applies the official GRC BED by verified assembly
aliases and zero-based half-open coordinates. Seven contamination annotations
and 56 complete false-duplication gene-body/exon-union overlaps become
`excluded_reference`, with empty counting targets and retained records. TRPM3
remains unchanged because its invalid-region overlap is intronic. Excluded
features are omitted from the eventual grouped counting GTF, not deleted from
the correspondence table. Count reading and matrix assembly omit empty targets.
Installed featureCounts single- and paired-end fixtures confirm that reads/pairs
with valid and excluded alignments are unassigned; there is no fractional
redistribution or alignment rescue through an excluded region.

A further exact-transcript pass produces `annotation_resolution_r7/`.
Shared RefSeq metadata alone does not establish gene identity. Classification
requires a versioned curated NM/NR RNA accession associated with only one
primary-chromosome gene, matching strand and biotype, no conflicting HGNC or
declared primary gene group, and a complete independently projected isoform
exon set identical to the primary isoform. Missing alignment bases, conflicting
groups and invalid reference flags cannot be rescued by a shared accession.
Thirteen unresolved records become equivalents and three become explicit
same-haplotype copies. The latter expose three earlier equivalents that also
must remain separate copies. The independent coordinate pass verifies every
new transcript proof, all 469 source records are rehashed, and the 86,369-ID
universe is unchanged. Current totals are 78,653 primary, 6,356 equivalent,
248 distinct, 36 alternative-only, 158 duplicated, 63 excluded-reference and
855 unresolved. No grouped counting annotation was emitted. Detailed source
checks are in `bulk_execution_sources/exact_shared_transcript_r6_review/` and
`annotation_resolution_r7/review_checks.json`.

`catalog_review_r6_reference_exclusions/` retains all 28,017 object records and
their historical source classes. Fourteen represented reference artifacts now
have non-applicable bulk RNA calls/applicability, retaining the original RNA
call in a separate source field. Gene IDs join by exact versions or the source
Catalog's explicit versionless stable IDs when the same GENCODE release has
only one exact version; symbols are never used. Other assays, target links,
gene routes and the four worked examples are unchanged. The public Catalog
has not been replaced. Reference records absent from that source Catalog remain
in the complete annotation table.

The revised NMF protocol is concrete but is a sensitivity definition for review,
not the original 5,000-feature projection. It retains 4,517 unchanged one-to-one
counting objects for full stage bases and 4,438 after the existing continuum
signature exclusion. Any target receiving an equivalent ALT annotation is
masked even if that ALT was outside the original basis. Unresolved, excluded,
duplicated and merged source objects cannot be filled with zeros or replicated
across rows. All retained k4/k6 weights are bit-identical to their original
values; no basis fitting or outcome modeling was performed.

Normalization uses one serialized historical disease-library reference,
SRR21623103, selected by highest `full_called` in the qualified original
training roster with lexical sample-ID ties. Pairwise edgeR TMM ratios on
26,613 unchanged objects are anchored to that reference's original factor.
The ratio cancels edgeR's pairwise geometric centering, so adding or removing
other samples cannot change a sample's normalization. The installed voom-style
0.5-count/one-library-count offsets, historical gene-by-cohort subtraction,
row minima and clipping at zero are retained. Unknown cohorts are rejected.

Masking alone passes all ten historical reconstruction checks, but the complete
normalization diagnostic must govern the protocol. Across the 1,260 native
libraries, including 1,104 original disease-training libraries and 156 other
libraries, all scores are finite and nonnegative. Eight full-stage axes and
all ten signature-excluded continuum axes meet the prespecified criteria:
global Spearman >=0.99, every historical cohort >=0.90, retained basis mass
>=0.80 and full retained basis rank. Stage k4-P2 (rho 0.9897164) and k6-P6
fall below the global threshold and remain unavailable; the criteria were
not relaxed. These are native-count diagnostics, not corrected scores or
control-group inference. Source tables and the final serialized protocol are
under `bulk_execution_sources/revised_fixed_nmf_projection_r7_review/` and
`bulk_execution_sources/revised_fixed_nmf_projection_r7_transport_checks/`.
Use `revised_fixed_nmf_projection_protocol_with_transport_checks_r2.rds`
and `final_protocol_review_r2.json` in the review directory. The earlier
metadata serialization rendered an empty unavailable-axis vector as `k-`;
the corrected serialization checks vector lengths and reports no unavailable
continuum axes. Projection values and eligibility decisions are unchanged.
The partially unavailable stage scores must not be passed to the historical
finite-score stage input contract as if all six programs were supported.

Executed verification uses the existing 21981199 allocation, with no new job,
array or dependency installation: all 50 focused Python tests pass; the
installed single/paired exclusion fixture passes; the synthetic matrix assembly
retains supported units and omits the excluded row; and the NMF helper checks
confirm masking, fixed-reference sample independence, finite transforms and
unknown-cohort rejection. The final NMF protocol independently checks unchanged
weights and the historical reconstruction decisions. The Resource scope check
passes. The updated implementation files are
`resolve_annotation_correspondence.py`, `count_verified_fragments.py`,
`assemble_verified_counts.R`, their focused checks, the source-review scripts,
`scripts/portal/apply_catalog_reference_exclusions.py` and the fixed NMF
transport/preparation/check scripts under the molecular-layer analysis directory.

**Historical status at the end of the 2026-10-05 work:** the following
unavailability statement precedes the October policy decision and recount.

The remaining 855 identities are a scientific source-evidence limitation.
Excluded reference artifacts do not license exclusion or equivalence of the
other unknown genes. The required zero-unresolved primary annotation,
corrected 776-library union counts, primary bulk fit and dependent corrected
stage/continuum/Catalog estimates remain unavailable. No claim of complete
annotation resolution or synchronized scientific adoption is made.

## Counting and refit completion, 2026-10-07 to 2026-10-11

On 2026-10-07 the user decided to hold the remaining 855 annotations as
untestable instead of requiring individual resolution before the primary fit.
The r8 annotation retains their counts but excludes them from testing; the other
85,514 annotation dispositions and counting targets are unchanged. The primary
rule excludes reads shared between a tested gene and a held annotation. The
single-candidate credit arm is a sensitivity, covering 403 of the held records.

The completed 776-library recount is in
`figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/`
and contains 79,095 counting rows. The corrected-count fit dated 2026-10-08 is
`RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/`.
It tests 24,115 genes in 774 genetic-identity clusters, one library each:
138 controls and 636 disease units; 1,854 canonical calls, 1,333 up and 521 down.
The 738-unit strict-depth sensitivity gives 1,906 calls. The corrected-count
fixed-NMF projection was scored on 2026-10-09. These completions supersede the
earlier statements that recount and refit were unavailable; the 855 annotations
remain untestable, not resolved.

On 2026-10-11 the user adopted the 774-unit bulk refit and selected Gene Catalog
r8 for the paper. The r8 classes are 400 genetically anchored, 26 concordant,
1,129 established-state associated, 13,127 unresolved and 13,551 untested over
28,233 genes. Genetic memberships still depend on unadopted COLOC r2. The
decisions and source references are recorded in [STATUS](../STATUS.md#what-is-official-updated-or-open)
and [RESULTS](../RESULTS.md#human-bulk-rna-seq); they do not constitute adoption
of the synchronized Resource release.

## Scientific review, 2026-10-12

**Amended 2026-10-13 after review of the review.** The adjudication below records
which criticisms were confirmed, which require qualification, and which checks
remain unperformed. Numerical reproduction of an erroneous coordinate parser
did not validate the Figure 2 scientific units.

This is a review of the current Resource, with five parallel assessments of
genetics, bulk transcriptomics, single-cell programs/continuum, spatial/Catalog
evidence, and predictive models. The review inspected the five human authorities,
manuscript drafts, relevant analysis code, small result tables, and rendered
main panels. Selected new hypothesis results were checked where they affect
the main conclusions. It is not an exhaustive reanalysis of every candidate.
No scientific code, figure, result, or evidence class was changed by these
review sessions. The repository edit is this document, including corrections
to its historical labels and completion history.

The review uses the 774-unit bulk refit adopted on 2026-10-11. Those units are
genetic-identity clusters, not independently verified people. September COLOC
r2 remains unadopted. The r8 Catalog is the user's selected paper Catalog, but
its genetic classifications depend on that COLOC candidate. Earlier sections
of this dated outcomes document retain their original historical status.

**Overall assessment.** The Resource has substantial support for cross-sectional
multicellular tissue remodeling, a phenotype-qualified inherited-trait map,
and calibrated assay context. Its strongest contribution is making those
different measurements interpretable together. The reviewed evidence does not
establish cellular regulation from tissue covariance, superior longitudinal
disease measurement, or biological equivalence from nonsignificant comparisons.
Several manuscript sentences still exceed boundaries already recognized in
PAPER and RESULTS. These require scientific wording corrections, not a broad
new analysis campaign.

### Findings that hold up

| Main result | Evidence checked | Defensible interpretation and limit |
|---|---|---|
| Bulk disease remodeling | Independently recounted 24,115 tested genes, 1,854 canonical calls, 1,333 up and 521 down; 774 records, 138 controls and 636 disease units. [Bulk authority](../RESULTS.md#human-bulk-rna-seq). | Cohort- and inferred-sex-adjusted disease association. The canonical list uses BH plus a post-test effect-size filter; it is not a calibrated interval-null list. |
| Broad directional robustness | [Leave-one-cohort-out table](../../RNA-seq/results/manuscript_release/candidates/bulk-loco-summary-20261009T143029Z/loco_summary.tsv) confirms 66.7–96.8% same-sign call retention. Of 510 still-tested calls lost on removing GSE213621, 505 retain direction; 1,026 calls survive all five removals. | A broad directional signal recurs. GSE213621 supplies 326/774 units, 68/138 controls and 45.1% of quality weight; significance is materially cohort-dependent. |
| Different fine-mapped architectures across trait classes; amended coordinate check | The initial 6.6906 ratio repeated the faulty decimal-position parser. Unique restoration against [native member coordinates](../../GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv) gives 23.3428% versus 3.8732%, ratio 6.026671; 233 pooled clusters, 19 direct and 228 enzyme, with 14 shared; the same 207/6,208 credible-set instances. | The descriptive direction persists. The old interval and 62-refit range are not validated for the corrected coordinates; reruns are required. This is a portfolio-conditional coding/splice fraction, not a coding share of MASLD heritability. |
| Predominantly enzyme-trait regulatory map | Independently counted 436 SuSiE-supported genes: strongest trait AST 158, ALT 137, GGT 128, versus NAFLD 8, PDFF 4 and NASH 1. | 423/436 strongest signals concern enzymes, not direct disease or liver fat. The phenotype qualification belongs beside every inherited-risk headline. GNMT remains a provisional ABF-supported direct-NAFLD example, not a primary SuSiE-supported mechanism. |
| Spatially organized Stromal ECM program | [All-composition comparison](../../Analysis/Multimodal_Program_Projection/candidates/spatial-full-composition-sensitivity-candidate-2026-08-11-r2/comparison_to_primary.tsv) confirms residual Moran excess 0.04050 and 0.05034, with BH q=0.001 and 0.0006. | Organization exceeds matched-gene background after substantial composition attenuation, 70.86% and 65.78%. Biological replication is four donors in GSE192741; Vu has ten arrays without a donor join. Matched-gene precision is not uncertainty across people. |
| Additional donor-resolved spatial candidate | Independently recounted 32/32 positive centered effects in V4, mean 0.098272. The one-sided exact sign-flip nominal p is 2^-32=2.3283e-10; Holm-adjusted p=4.6566e-10 across the prespecified HMSMA/Vu two-test family. [Exact summary](../../Analysis/Hypotheses/V4-20261008T015029Z/summary_signflip_exact.tsv). | Useful additional matched-background tissue evidence. [Inference] The claimed 28 MASLD/MASH plus four-control roster depends on the inferred exclusion of three controls, without an authoritative person-to-clinical-label join. Adjustment uses a hepatocyte-marker proxy rather than the canonical 16-factor model. No disease/control or stage comparison is established. |
| Task-specific model value and limits | [Two-model results](../RESULTS.md#two-model-molecular-value-development-2026-09-22-not-adopted) retain global-plus-local RNA and native 1-Mb sequence. The combination's 3.15% MSE reduction misses its fixed 5% complexity margin. | Same-source molecular prediction and cross-sectional severity ranking are useful candidate findings. Superiority to all task-native methods, external H3K27ac transfer, and longitudinal tracking remain unestablished. |

The fixed 117-program definitions, complete testing families, donor collapse,
explicit unavailable assays, matched spatial nulls, and retention of failed
model comparisons are strengths. Source-owned results and the separate
experimental follow-up must remain distinct from Resource discoveries.

### Scientific findings requiring correction or clearer limits

**1. Repeat-biopsy variability is still used to support a measurement ceiling
despite unresolved identities.** PAPER excludes GSE193066 from primary
replication, ICC and derived ceiling claims, but
[Figure 7 prose](../manuscript/draft/fig7.md) still organizes performance around
0.4655–0.6433 stages of repeat-biopsy variation and 28/58 exact agreements.
Only three of those 58 pairs satisfy the released score's input bounds.
Different-visit variation includes true change, tissue sampling, reading error,
and technical effects; it does not identify reading-error variance alone.
Independent primary evidence also distinguishes biopsy sampling from observer
variability ([Ratziu et al., 2005](https://pubmed.ncbi.nlm.nih.gov/15940625/)).

**Suggested documentation:** remove GSE193066-derived ceiling and primary
reliability conclusions from manuscript claims while identity remains unresolved.
Describe other repeat-biopsy estimates as conditional variability scenarios,
with their assumptions, rather than a universal measurement bound. Keep
mathematical ceilings caused by ordinal ties separate from error-based ceilings.

**2. The new N4 random-pair comparison uses different quality eligibility from
the real pairs.** The
[executed script](../../Analysis/Hypotheses/N4-20261007T210855Z/n4_within_person_agreement.as_run.R)
requires real GSE193066 pairs to pass N1 and the existing same-person call
(lines 33–50), but its random-pair pool omits the N1 restriction (lines 152–160).
The independently reconstructed pool contains 48 libraries, including 15 called
mixed. [Inference] This mismatch can change the claim that 87/108 programs
exceed the null; its direction is not established. N1's non-call is not proof
of purity or identity: the
[simulation table](../../Analysis/Hypotheses/N1-20261007T210855Z/positive_control_summary.tsv)
has 38% sensitivity at a 2% mixture in the mid-depth stratum, despite 88% at 10%.

The 22-pair N4 result comprises 19 GSE193066 and three GSE130970 different-visit
pairs. Median program ICC 0.706 versus stage ICC 0.074 describes selected-pair
agreement; a persistent person-specific signal could also give a higher ICC.
It does not demonstrate better histology measurement or tracking of disease
change. No direct selection on program agreement was found.

**Suggested documentation:** mark the null comparison provisional in STATUS,
name the asymmetric eligibility, and use “greater agreement across selected
different-visit pairs” for the descriptive result. A later correction should
apply the same library-quality rules to observed and random pairs before the
null claim is used. No corrected null was run in this review.

**3. Continuum association is shared RNA covariance, and does not establish
system-specific biology.** The
[Hotspot continuum producer](../../scripts/analysis/histology_anchored_continuum/molecular_layers/31_hotspot_continuum.R)
excludes signature genes (line 175) and explicitly records the same expression
substrate (line 225). This removes direct component overlap; it does not make
the remaining outcomes independent measurements. The two continuum scores are
correlated at approximately 0.999. The [corrected 774-unit stress summary](../../RNA-seq/results/histology_anchored_continuum/molecular_systems_stress/hac-system-stress-corrected774-20261009T171017Z/tables/system_stress_summary.tsv)
records 33/43 systems meeting joint maxT on both axes, 0/43 exceeding the
connected-region competitive null on both axes, and 1/43 exceeding the
unconstrained size-matched null on both axes. Those tests answer different
questions; the stage-adjusted roster is 435 units.

Composition sensitivities further qualify the interpretation. The standard
16-type [774-unit composition candidate](../../figures/candidates/continuum-composition-774-20261012T012137Z/basis_774/continuum/continuum_meta_effects.tsv)
gives ECM 0.294 with rebuilt fractions and 0.479 with earlier fractions, versus
0.760 unadjusted. The low end 0.203 belongs to the separate
[rare-type pooled sensitivity](../../figures/candidates/continuum-rare-type-arms-20261012T020132Z/arms/summary_arm_by_fraction_set.tsv)
(arm B, 13 parts/12 ILR coordinates). Across those distinct composition
codings Ductular injury remains about 0.60–0.73. Both candidates use 367
disease-only units and HC3-based two-cohort models. The organizational
interpretation already in PAPER is appropriate.
The inherited continuum method itself was published by
[Kamzolas et al., 2026](https://www.nature.com/articles/s42255-026-01543-7);
the Resource contribution is its qualified cross-resource use.

**Suggested documentation:** bring these limits into the abstract, Figure 4
Results and title rationale. Use “membership-defined systems summarize shared
cross-sectional transcriptional variation.” Do not use the number of associated
systems as the number of independent biological discoveries or treat RNA
covariance as independent validation of severity or cell-to-cell regulation.

**4. A failed dimensionality criterion is overinterpreted as biological
non-separability.** [PAPER](../PAPER.md), lines 170–179, treats a participation
ratio of 1.42 below the prespecified 2.0 floor as grounds that nothing can be
attributed to individual histology components. The
[saved result](../../Analysis/MASLD_Model_Benchmark/executions/model-check-1044-21181334-aspect-separability/separability/aspect_separability.json)
contains eigenvalues 2.4858, 0.3376 and 0.1766. Independent arithmetic reproduces
the ratio 1.423053 and correlation-matrix determinant 0.148187: the measured
label matrix has full rank, with 17.14% of variance outside the leading
component. This does not establish useful distinct biological axes, but neither
does it make conditional component effects mathematically undefined.

**Suggested documentation:** preserve the prespecified stopping decision as
“the proposed three-axis model did not meet the separability criteria in this
99-participant cohort.” Restrict any tested-negative label to that exact
criterion. Component-specific conditional associations are a separate question,
with their own uncertainty and validation requirements. The later H8 candidate
must be evaluated on those requirements, not accepted or rejected from this
global participation ratio.

**5. A male-only cohort is incorrectly described as lacking an estimable disease
contrast.** The per-cohort routine in
[build_pi_stage_extensions.R](../../scripts/figures/build_pi_stage_extensions.R)
(line 288) rejects a cohort if inferred sex has only one level. The source
[screened roster](../../figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv)
has 25 male controls and 31 disease records, including ten women. The identity
screen excludes 24 disease records (14 men and all ten women) and two controls,
all for insufficient filtered identity sites. The selected roster therefore has
23 controls and seven disease units, all male; the source cohort was not
male-only. [Inference] This unequal eligibility may affect representativeness.
Within the selected roster, an intercept-plus-disease coefficient is estimable
when the constant sex term is omitted. That modeling correction does not resolve
the upstream selection concern.

**Suggested documentation:** replace “does not support a per-cohort estimate”
with “the identity-selected roster is male-only and omitted under the fixed
sex-adjusted modeling rule; its disease contrast remains estimable.” State the
source sex composition and exclusion counts alongside it. A later authorized
analysis can omit the constant nuisance term; separate sensitivity is needed
to assess the selected roster's representativeness.

**6. The claim vocabulary sometimes changes the measured object or null.**
[Figure 5 prose](../manuscript/draft/fig5.md) calls lower GNMT–MAT1A–CYP2C19
RNA/protein abundance loss of metabolic function, without an activity or flux
measurement. [Catalog prose](../manuscript/draft/fig5_discussion.md) assigns
spatial indeterminacy to BICC1, although the spatial test concerns the weighted
Ductular injury program; BICC1 contributes 5.44% of its weight. The
[DEG-threshold sensitivity discussion](../RESULTS.md#deg-threshold-sensitivity)
says colocalized genes respond “exactly as much” after nonsignificant matched
comparisons, without an equivalence margin. It also says “complementary evidence
maps without evidence of biological antagonism.” The original review wrongly
called that biological complementarity. The phrase should be qualified as
operational organization, because the class-difference language rule prohibits
presenting complementarity of the two gene classes as a biological finding.

**Suggested documentation:** use “lower metabolic gene/protein abundance,”
“the Ductular injury program remains spatially indeterminate,” and “no difference
was detected under the specified comparison.” Keep functional loss, a named
gene's spatial behavior, and biological equivalence as separate hypotheses.
Likewise, canonical 1,854-gene membership must always be described as BH
point-null significance plus an observed effect filter. The existing interval-null
sensitivity tests a different threshold, 0.25 rather than 0.5, and is not a
calibration of that exact canonical list.

**7. Focal program disease support is less secure than “robust programs” implies.**
The [ambient table](../../Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/results/ambient_program_effects.tsv)
confirms original 64-donor/five-dataset focal HC3 BH q about 0.105. In the
62-donor/four-dataset common universe, the raw HC3 BH q is 0.123 for both
programs. In that paired comparison ECM changes from 0.401 to 0.130 and
Ductular injury from 0.396 to 0.389;
paired-change BH q is 0.271 and 1.00. Neither corrected native-lineage disease
association clears the complete family. The lack of a significant change does
not establish ambient resistance. Targeted fibroblast follow-up addresses a
different lineage and selected two-program family, not independent replication
of the original hepatocyte interpretation.

**Suggested documentation:** retain the attenuation-indeterminate state and
distinguish discovery lineage, content-associated lineage, and disease-tested
lineage in each claim. Correct the focal denominator in the
[spatial-program discussion](../RESULTS.md#prespecified-programs-with-different-spatial-organization-stromal-ecm-and-ductular-injury):
the 62-donor models have four contributing datasets, although five datasets
underwent correction.

**8. Working figures and manuscript legends presently state different numbers.**
These are actual reader-visible discrepancies, several already acknowledged
in RESULTS. Text was extracted from the PDFs and the relevant panels visually
inspected; this review did not install replacements.

| Panel | Current rendered content or source-table values | Comparison, with its analysis basis |
|---|---|---|
| 3C | Five cohorts including GSE126848; 9,777 study-unique genes | Current four-fit candidate has 8,065 study-unique genes; its omission of GSE126848 also needs finding 5 addressed. |
| 3E | Bar source counts F1–F4 versus F0: 1,338 / 4,266 / 6,190 / 7,043; five displayed k6 axes. The exact counts are not printed inside the PDF. | Each stage count uses BH q<0.05 with no fold-change filter. These differ from the canonical disease/control rule. Caption still says 23,370-gene families and six axes; current family is 24,115, and k6-P6 is unavailable. |
| 3F | 61 genes, 244 adjacent-stage effects | 774-unit candidate: 54 genes, 216 effects; GNMT F1→F2 q=0.0562. |
| 3G | Neither 3,036; RNA only 180; genetics only 383; both 23 | 774-unit candidate: 3,156 / 194 / 400 / 26. |
| 5B | 211 protein-significant; 185 concordant, 26 discordant; 64% direction-concordant | Legend describes September's 844-record candidate: 210 / 183 / 27 and 65%. No 774-unit counterpart was located in inspected candidate outputs. |
| 5C | Current working HLA-DQA1 RNA marker is hollow | Legend claims nine filled RNA markers including HLA-DQA1 from September's 844-record candidate. The 774-unit primary annotation has log2FC 0.379801, padj 0.024969 and fails the canonical effect-size criterion. |
| 5D | 29 direct-trait and 787 enzyme pairs | Legend describes the September COLOC r2 candidate's 28 / 739; this is a genetics correction, not a 774-unit bulk rerun. |
| 6B/6D | 30,754 genes; 33 concordant | r8 uses the 774-unit fit and candidate COLOC r2: 28,233 genes; 26 concordant. Legend already claims corrected inputs. |

The original review also missed [Figure 3 manuscript prose](../manuscript/draft/fig3.md),
which still describes 844 participants, 140 controls, 704 disease records and
23,702 genes. Those are historical 844-record inputs, not the adopted 774-unit
fit; neither roster should be called verified participants without identity
qualification.

**Suggested documentation:** legends must first describe the PDF actually
shown and identify pending replacements. After an authorized replacement,
update the caption and assembly together. In the Figure 3 stage interpretation,
state that F0 contains 14 source controls plus 89 F0 disease biopsies: the F1
contrast yields 1,338 calls with all F0, 783 with disease-only F0, and 814 with
a control-status term. This qualifies the early-stage magnitude without
turning cross-sectional comparisons into temporal onset.

**9. Distinct identity and evidence-state definitions need explicit names.**
The categorical overlap really does contain 14,720 jointly testable genes in
its [source audit](../../RNA-seq/results/manuscript_release/candidates/evidence-classes-corrected-r8-20261009T142027Z/orthogonality_audit.tsv).
The [scatter interface](../../RNA-seq/results/genetics_state_interface/20261009T142027Z_coloc_r2_corrected_controls_r8/overlap_by_arm.tsv)
really does contain 14,704. The former joins by symbol with primary-assembly
preference; the latter joins stable Ensembl IDs. Independent comparison found
24 categorical-only IDs and eight interface-only IDs. This is not an arithmetic
error. [Unverified] The biological equivalence of those annotation mappings
was not established. Both give 26 observed versus about 33.4 expected overlaps,
nominal p about 0.20, so the chance-compatible interpretation is unchanged.

Separately, r8 calls 708 genes “untested” despite support in one domain, because
the other domain is untestable. It calls 11,301 BH-significant smaller RNA
effects indeterminate because of the effect-size rule; 555 indeterminate
genetics rows carry ABF PP.H4 above 0.5 without a supported SuSiE call. These
are documented in [RESULTS](../RESULTS.md#masld-gene-catalog), but a reader could
misread the class names as absence of any evidence.

**Suggested documentation:** name both overlap universes and reconcile the
32 differing identity records before describing one shared universe. Do not
automatically replace 14,720 with 14,704. Describe Catalog “untested” as incomplete
joint characterization; show assay-specific evidence and distinguish ABF from
SuSiE statistics beside each gene.

**10. Remaining genetic and model qualifications should accompany their headline
numbers.** The blanket PDFF coordinate uncertainty in the original review is
withdrawn for SuSiE and COLOC r2. Actual September r2 logs for
[chromosome 7](../../GWAS/finemapping/logs/coloc_eligibility/coloc_elig_21862335.out)
and [chromosome 2](../../GWAS/finemapping/logs/coloc_eligibility/coloc_elig_21862339.out)
read the AF-augmented copy identified in the SuSiE study manifest. The
[derived-copy checks](../../GWAS/orientation_checks/pazoki_pdff_derived_copies-20261010T145917Z/shell_validation.tsv)
show that correction changed allele labels only; the non-allele columns were
byte-identical. The double-lifted coordinates belong to two separate hg19 copies
not read by these pipelines. Their unsigned posteriors are invariant to the
uniform sign reversal. Only C2orf16 and RABGEF1 among the 436 primary genes have
Pazoki SuSiE support, and both also pass the primary threshold in another study.
[Unverified] Direction-sensitive and model outputs beyond those inspected
pipelines remain untraced.

The review missed a separate input-state mix in [PAPER](../PAPER.md). Its current
overlap paragraph identifies the 774-unit fit, then quotes the historical
class-difference test's 272 pairs, mean difference -0.44, sign-flip p=2e-5 and
71 protein pairs without changing the stated basis. The
[class-difference verdict](../../RNA-seq/results/class_difference_test/cdt-corrected-20260925T203546Z/VERDICT.md)
names the September 844-record corrected-controls fit. No 774-unit rerun was
located in the inspected class-difference outputs. Label these as historical
comparisons or rerun them before attributing them to the adopted fit.

For Figure 7, the additional severity cohorts show rank correlations 0.6382
and 0.6739, but native-comparator differences of -0.040 [-0.099, 0.017] and
-0.011 [-0.079, 0.054] do not establish superiority or equivalence. The RNA-to-
H3K27ac exact-histology comparison has 22 disjoint pairs, almost entirely women
with F0, and no F2–F4 pairs. The allele adapter was selected non-nested on the
development folds; saved-prediction bootstrap intervals omit selection and
refitting uncertainty. These limits are mostly already recorded in RESULTS.

**Suggested documentation:** lead with cross-sectional rank transfer and
same-source molecular prediction beyond the specified comparator. State the
selection, exposure and paired-sample limits at the claim, rather than only in
the technical record. Prediction of marginal caQTL coefficients does not
establish allele-editing accuracy, a target gene, or a MASLD effect.

### Proposed improvements to the existing documentation

These are recommendations, not adopted edits to paper scope or results.

| Document | Proposed change | Scientific purpose |
|---|---|---|
| PAPER | Correct the Figure 2 coordinate-dependent ratio and rerun its uncertainty before reporting quantitative robustness; label the 844-record class-difference test separately from 774-unit overlap. Narrow the coordination/title rationale, distinguish within-RNA covariance from independent validation, and correct the dimensionality interpretation. | Make the headline and numerical claims match the measurements and inputs. |
| RESULTS | Keep one clearly marked current entry per main claim, with unit, n, estimand, adjustment, testing family and input state; correct the four-dataset ambient denominator; label the two overlap identity universes; replace equivalence wording. | Make exact numbers scientifically interpretable without reconstructing historical paragraphs. |
| STATUS | Move superseded runtime and roster descriptions to the technical history; mark N4 null eligibility unresolved; separate completed estimates from conclusions still limited by identity or source exposure. State V4's mixed clinical roster and link its exact-probability table. | Prevent old limits and new candidate numbers from being read as one current conclusion. |
| ROADMAP | Prioritize the Figure 2 coordinate correction and uncertainty reruns, N4 comparison correction, estimable selected male-only cohort contrast, and figure/legend discrepancies. Separate needed claim corrections from optional new biological work. | Focus remaining work on claims a reader will actually see. |
| README | Give a short current-results map: adopted bulk, candidate genetics, selected Catalog with its genetic dependency, and qualified continuum/spatial/model findings. | Let readers identify what is currently supported without relying on dated introductory blocks. |
| Methods and figure legends | State actual inferential units and complete-case n; identify point-null versus interval-null testing, competitive versus zero-effect nulls, shared RNA measurements, source exposure and composition adjustment. Match each legend to the rendered panel. | Prevent changing the estimand between analysis, caption and interpretation. |
| Figure 6 dictionary and prose | Explain supported single-domain evidence inside an “untested” joint class, small significant RNA effects inside “indeterminate,” and ABF versus SuSiE support. | Make the Catalog useful without treating missing joint characterization as no evidence. |

### Verification and limits of this review

The initial working tree was recorded and unrelated edits were preserved.
Reviewers used bounded rg/sed/nl inspections, Python standard-library recounts
of small tables, pdftotext, and rendered PNG inspection. Independent arithmetic
covered the bulk-call summary, trait-driver counts, coding posterior fractions,
overlap denominators, spatial attenuation/V4 signs, and the aspect correlation
matrix. Relevant source code and the returned evidence were cross-checked by
the primary reviewer. Primary literature was checked for the published
continuum method, colocalization assumptions, and biopsy sampling variability.

New biological fits, bootstrap regeneration, participant identity confirmation,
comprehensive coordinate-dependency tracing, and review of all new hypothesis
tests: **not run**. The October 13 recheck did trace the PDFF inputs actually
used by SuSiE and COLOC r2. No jobs, arrays, package installations, commits,
branches, or external writes were performed. No result or manuscript-readiness
state was changed. Documentation checks are recorded after the report below.

The original documentation check passed scope and 19 local links. That check
did not establish the correctness of coordinate parsing, input-state labels or
scientific interpretation. The corrected report is checked again below.

### Review of this review, 2026-10-13

Four parallel reviewers adjudicated the 13 criticisms using current authorities,
source tables, executed input logs and rendered panels. The primary reviewer
independently restored coordinates and recomputed the point estimates without
importing either affected Figure 2 script. The review's main numerical error
was scientific: its purported independent reproduction repeated the same
incorrect position interpretation. Reproduction of that calculation was not
validation of its genomic units.

| Criticism | Adjudication | Correction in this document |
|---|---|---|
| 1. Figure 2 ratio and cluster census | Confirmed. 715/6,415 gated positions are truncated. Every set has one trailing-zero restoration within 1 Mb of every native member coordinate. | Replace the claimed validated 6.69 ratio with the restored point estimate 6.026671 (23.3428%/3.8732%; 233 clusters, 19 direct/228 enzyme/14 shared). Withdraw the old interval and robustness inference pending reruns. |
| 2. GSE126848 male-only roster | Confirmed omission. All ten women in the source disease group were excluded by the identity-site rule, alongside 14 male disease records and two male controls. | State source and selected rosters separately. The selected 23/7 contrast remains estimable, but that does not resolve upstream eligibility or representativeness. |
| 3. Figure 5 input basis and HLA-DQA1; stale Figure 3 draft | Confirmed. The September 5B/5C candidate uses 844 records; no 774-unit counterpart was found in inspected candidates. HLA-DQA1's primary 774-unit annotation fails the canonical effect-size criterion. | Label historical Figure 5 numbers; distinguish the old alternate annotation's call. Record the stale Figure 3 prose. Qualification: the current rendered HLA-DQA1 marker is hollow, so the filled-marker problem is in the textual legend/candidate attribution. |
| 4. PAPER class-difference input mix | Confirmed. The 272-pair comparison is from the September 844-record fit, while the preceding overlap numbers are from 774 units. | Add the missed input-state mix and require historical labeling or a new 774-unit comparison. No such rerun was found in the inspected class-difference outputs. |
| 5. “Biological complementarity” misquote | Confirmed mischaracterization. The source says “complementary evidence maps without evidence of biological antagonism.” | Quote faithfully and distinguish operational organization from a biological finding about the two classes. Retain the valid warning that a nonsignificant difference does not establish equivalence. |
| 6. Stale RESULTS line numbers | Confirmed. The quoted passages moved to 954, 985 and 611 at recheck. | Replace brittle line-only references with section links and identifying text. |
| 7. V4 tail, correction and control labels | Confirmed. The test is one-sided; nominal p=2.3283e-10, Holm p=4.6566e-10 over HMSMA and Vu. Roster composition depends on an inferred exclusion without a clinical join. | State tail and family correction and mark the roster interpretation as inference. |
| 8. Ambient q values | Needs source qualification. q≈0.105 was correct for the explicitly original 64-donor estimate; 62-donor raw HC3 q is 0.123. | Label the 64/5 and 62/4 universes separately before comparing coefficients. |
| 9. ECM range and system source | Needs source qualification. The 0.203–0.479 range and 33/0/1 system counts are valid across their stated analyses, but the review did not identify those sources precisely. | Name the standard 16-type 0.294–0.479 comparison and separate 13-part rare-type arm's 0.203 estimate. Link the corrected 774-unit stress root and its both-axis decision rule. |
| 10. Figure 3E rule and Figure 5B percentage | Confirmed. Stage bars use BH q<0.05 without the canonical disease/control effect filter. Working 5B prints 64%; the September candidate/legend says 65%. | State the distinct stage estimand and the actual rendered percentage and input basis. |
| 11. PDFF uncertainty | Confirmed for SuSiE and COLOC r2. Executed logs read the coordinate-correct AF copy; its correction changes allele labels only. Two of 436 primary genes have support from that study, and neither depends on it alone. | Withdraw the blanket coordinate uncertainty for these pipelines. Keep untraced direction-sensitive/model outputs explicitly unverified. |
| 12. Document state and missing history | Confirmed. The blanket unadopted header and September Catalog section were stale, and recount/refit completion was absent between preparation and adoption. | Correct the header, mark the older Catalog superseded, and add the October 7 policy/annotation decision, completed 776-library recount, corrected refit and October 11 decisions. Qualification: earlier manifest/preparation records already existed. |
| 13. Prohibited method shorthand | Accepted as the current wording instruction. | Replace the shorthand in this document's explanatory prose with interval-null terminology; preserve scientific test thresholds and source identities. |

**State at the end of the review, before the authorized reruns below.**
The unresolved analysis work is narrower than this documentation correction:
repair coordinate parsing and derived cluster expectations, then rerun the
Figure 2 bootstrap, context families and 62 robustness comparisons; reconstruct
774-unit Figure 5B/5C and class-difference results before calling them current;
and correct N4's unequal quality filtering before using its random-pair claim.
These analyses were **not run** in this review. Existing report findings about
N4 eligibility, within-RNA covariance, gene/program object distinctions and
abundance versus function remain supported by the checked evidence.

Recheck verification: the primary reviewer independently recomputed the restored
Figure 2 point estimates from 6,415 gated sets and checked the actual PDFF
execution logs and input manifest. The scope command returned
RESOURCE_SCOPE_VALIDATION PASS (37 documents, seven figure roles). All 31 local
links in the amended review resolve, including their Markdown section anchors;
the prohibited shorthand has zero occurrences in this document. Two reviewers
checked the amended genetics and document-state interpretation. Bootstrap and
robustness reruns, new Figure 5 and class-difference fits, and the N4 correction
remain **not run**. These documentation corrections do not adopt new results.

## Coordinate parser repair and reruns, 2026-10-13

The user authorized the Figure 2 coordinate fix, focused trailing-zero checks,
analysis reruns and necessary job submission. One compute job completed all four
affected analyses. The corrected coding/splice ratio is **6.026671**, with a
95% physical-cluster bootstrap interval of **2.574593–13.510165**. Its complete
57-comparison robustness menu gives ratios **4.368132–8.275691**; every difference
interval has a positive lower bound. These are candidate results, conditional
on the studied GWAS portfolio. They do not adopt a synchronized Resource release
or replace named releases or Figure PDFs.

### Coordinate repair and focused checks

The shared parser restores trailing zeros in legacy chromosome.position values
using each credible set's full native member positions. For example, GCKR's
`2.2773094` resolves to chromosome 2, position 27,730,940. A restored sentinel
must lie within the original 1-Mb window of **every** member and be the unique
supported decimal extension. The parser rejects absent members, chromosome or
metadata disagreement, invalid coordinates, unsupported encodings and ambiguous
restorations. It preserves the input order and does not substitute a member mean
or assume that the sentinel is itself a credible-set member. Recovery uses the
fixed original window even when a sensitivity comparison changes clustering.

| Changed analysis file | Purpose |
|---|---|
| [fig2_locus_coordinates.py](../../scripts/figures/fig2_locus_coordinates.py) | Shared coordinate recovery and explicit failure checks |
| [fig2_multiplicity_correction.py](../../scripts/figures/fig2_multiplicity_correction.py) | Correct pooled coding and noncoding-context clustering, coordinate table and cluster census |
| [fig2_coding_mass_robustness.py](../../scripts/figures/fig2_coding_mass_robustness.py) | Correct every omission/window/LD-block comparison and derive the complete menu from corrected clusters |
| [fig2_trait_directness_stats.py](../../scripts/figures/fig2_trait_directness_stats.py) | Correct descriptive trait-directness clustering while preserving its context-cluster mapping |
| [fig2_ascertainment_control.py](../../scripts/figures/fig2_ascertainment_control.py) | Correct another copied parser; preserve its separate ungated design and already-correct native lead-effect lookup |
| [test_fig2_locus_coordinates.py](../../scripts/figures/test_fig2_locus_coordinates.py) | Twelve checks, including GCKR, multiple lost zeros, unaltered positions, endpoints, ambiguity and a fixed-source scientific regression |
| [test_fig2_cluster_inference.py](../../scripts/figures/tests/test_fig2_cluster_inference.py) | Update the parser fixture; retain nine checks of shared-cluster resampling, gating and seeded inference |
| [run_fig2_coordinate_rerun.sbatch](../../scripts/figures/run_fig2_coordinate_rerun.sbatch) | One sequential compute job with focused checks, recorded environment, copied small inputs and before/after input checksums |

The fixed-source regression independently sorts and joins adjacent loci and
recomputes posterior weighting without calling the producer's cluster function.
It reproduces 6,415 gated instances, 715 restored positions, 233 pooled clusters
(19 direct, 228 enzyme, 14 shared), and fractions 0.23342751888576002 versus
0.038732412743244744. All **21 executable checks passed on the compute node**.
AST inspection of the changed Python files, `bash -n` on the launcher and
source-diff whitespace checks also passed. No dependencies were installed.

### Submitted job and completed outputs

Executed submission: `sbatch --parsable scripts/figures/run_fig2_coordinate_rerun.sbatch`.
Job **22218816** used account/QOS `nslab`, partition `cpu`, four CPUs, 48 GB RAM
and a 90-hour time limit. It completed on `ne1dc4-076` in **00:17:06**, exit
`0:0`; peak batch memory was **2,065,900 KiB (1.97 GiB)**. Running: 0; pending: 0;
done: 1; errors: 0. No array or downstream jobs were needed.

| Analysis | Completion | Principal output |
|---|---|---|
| Coding/context bootstrap | 6/6 ████████ | [fig2_family_multiplicity.tsv](../../figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/fig2_family_multiplicity.tsv) |
| Coding robustness | 57/57 ████████ | [coding_mass_robustness.tsv](../../figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/robustness/coding_mass_robustness.tsv) |
| Trait-directness estimates | 1/1 ████████ | [Completed analysis directory](../../figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/trait_directness) |
| Ungated ascertainment sensitivity | 1/1 ████████ | [ascertainment_control_summary.tsv](../../figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/ascertainment/ascertainment_control_summary.tsv) |

The [completion record](../../figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/COMPLETE.txt),
[coordinate checks](../../figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/coordinate_tests.log)
and [cluster checks](../../figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/cluster_tests.log)
record successful execution. All 57 input/snapshot checksum entries were present
and unchanged before and after execution. Required output tables were nonempty.
The existing environment recorded Python 3.10.18, NumPy 2.2.6, pandas 2.3.3 and
statsmodels 0.14.5. Large summary statistics were hashed on compute, not copied.

### Corrected estimates and interpretation

| Quantity | Corrected candidate result |
|---|---|
| Gated credible-set instances | 207 direct / 6,208 enzyme |
| Coding/splice posterior fractions | 23.3428% direct / 3.8732% enzyme |
| Physical 1-Mb coding clusters | 233 pooled; 19 direct / 228 enzyme / 14 shared |
| Coding/splice ratio; 95% interval | 6.026671 [2.574593, 13.510165] |
| Difference of fractions; 95% interval | 0.194695 [0.069229, 0.339940] |
| Noncoding-context clusters | 227 pooled; 18 direct / 223 enzyme / 14 shared |
| Robustness ratios across all 57 comparisons | 4.368132–8.275691 |
| Smallest robustness difference-interval lower bound | 0.022373 |

All primary and robustness intervals use 10,000 finite physical-cluster bootstrap
draws and seed 20260817, combined with a stable test/fit identifier. A sampled
shared cluster carries both trait classes. The operational clusters are not
asserted to be independent LD blocks. The robustness reference uses a separate
seed stream and has ratio interval 2.584490–13.558382 with the same point estimate.
The menu comprises one reference, 13 direct-study omissions, 18 enzyme-study
omissions, **19** direct-cluster omissions, five clustering windows and one
LD-block comparison. The historical 62 comparisons contained 24 direct-cluster
omissions; correction changes the number of clusters, with every comparison
type retained. These overlapping sensitivities are not independent replications.

The five context differences and their intervals are recorded in
[RESULTS](../RESULTS.md#adopted-corrected-coordinate-results-2026-10-13).
Direct ABC mass remains sparse, 0.1260% versus 0.7071%; the ratio interval
includes zero. The negative difference interval alone does not establish
biological depletion. Primary permutation p/q are absent. Bootstrap tail
diagnostics, their BH columns, and legacy trait-directness p columns are not
calibrated confirmatory trait-class tests. The trait-directness producer preserves
its architecture-derived mapping of context clusters; its context table should
not be silently substituted for the separately clustered multiplicity family.

The ascertainment analysis retains its historical **ungated** 6,467-instance
design, with 208 direct and 6,259 enzyme instances. Its baseline has 19 direct
and 234 enzyme clusters and ratio 6.178788; the effect-matched ratio is 5.782896
over 400 matched instances. Three adjusted models completed on 6,466 instances
with observed effects and 239 pooled clusters. Native lead-position effect
lookup was already correct. These sensitivities do not replace the gated primary
estimand or establish trait-label exchangeability; their permutation diagnostics
remain unadopted.

The current Figure 2 summaries in [PAPER](../PAPER.md), [RESULTS](../RESULTS.md),
[STATUS](../STATUS.md) and the [figure legends](../manuscript/05_figure_legends.md)
now distinguish these completed candidates from the historical coordinate-based
estimates. Coordinate-dependent artwork still needs replacement and paper use
requires a separate decision. No Figure 2 analysis rerun remains outstanding
within this repair. The separate 774-unit Figure 5B/5C and class-difference
comparisons and N4 correction identified in the review were not run by this job.

Final documentation checks: `python3 scripts/manuscript/validate_resource_scope.py`
returned **RESOURCE_SCOPE_VALIDATION PASS** (37 documents, seven figure roles).
All 23 added local link targets across the five changed documents resolve,
including Markdown section anchors. Added prose and source diffs passed
whitespace inspection. An independent reviewer checked the current Figure 2
summaries against the completed result tables. Initial working-tree entries
were preserved; no commit, branch, push or pull request was created.

## Figure 2 adoption and regenerated panels, 2026-10-13

The user explicitly adopted the corrected Figure 2 results and requested panel
regeneration and synchronized legends, PAPER, RESULTS and STATUS. The adopted
source remains
`figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/`; its directory
name no longer defines its scientific state. The coding/splice ratio is
**6.026671 [2.574593, 13.510165]**, with 233 pooled clusters (19 direct, 228 enzyme,
14 shared). All 57 sensitivity comparisons and the corrected descriptive context,
trait-directness and ungated ascertainment outputs are adopted within their
existing inference limits. Historical coordinate-based intervals are superseded.
This decision does not adopt COLOC r2, other figure results or the synchronized
Resource release.

### Regenerated panels and printed values

| Installed individual panel | Current rendered values and source |
|---|---|
| [Figure 2D](../../figures/main/fig2_genetics/panels/Fig2D_pip_architecture_by_trait_directness.pdf) | Direct protein-altering 23.1%, other noncoding 73.5%; enzyme protein-altering 3.9%, synonymous/UTR 3.6%, other noncoding 92.5%. The [eight consequence estimates](../../figures/main/fig2_genetics/panels/Fig2D_pip_architecture_by_trait_directness_source.tsv) sum to one within each trait class. |
| [Noncoding-context companion](../../figures/main/fig2_genetics/panels/Fig2_companion_noncoding_annotation_gap.pdf), outside A-I | Credible-set-weighted direct/enzyme fractions: unresolved 52.11/60.43%, promoter 7.87/4.46%, accessible 7.32/4.15%, ABC 0.23/1.57%. Its [source](../../figures/main/fig2_genetics/panels/Fig2_companion_noncoding_annotation_gap_source.tsv) is the same corrected estimate table; the display uses `instance_weighted_fraction`. |

Figure 2D uses the locus-weighted column: complete direct-class estimates are
23.15% protein-altering, 0.20% canonical splice, 3.16% synonymous/UTR and
73.50% other noncoding; enzyme estimates are 3.87%, 0.00%, 3.60% and 92.52%.
The combined coding/splice [interval table](../../figures/main/fig2_genetics/panels/Fig2D_corrected_coding_interval.tsv)
retains the primary physical-cluster ratio and difference intervals. The
companion's [context intervals](../../figures/main/fig2_genetics/panels/Fig2_companion_corrected_context_intervals.tsv)
describe the separate region-balanced estimator and are not uncertainty intervals
for its instance-weighted points. No significance labels are displayed; ABC
overlap remains sparse and does not establish biological depletion.

The [renderer](../../scripts/figures/fig2_trait_directness_panels.R) now derives
caption percentages and consequence-region counts from the input tables instead
of hardcoding 22.4/3.4% and 24/259 clusters. It exports the primary corrected
intervals separately and checks that their coding ratio matches the plotted
components. A consequence-only count check preserves the distinction from
18/223 noncoding-context clusters. Historical permutation/BH conclusions were
removed from both captions. Existing panel dimensions, consequence colors and
layout were preserved; companion text uses 6 pt throughout.

### Compute job, checks and installation

Executed submission twice after correcting utility availability:
`sbatch --parsable scripts/figures/run_fig2_coordinate_panels.sbatch`.
The [launcher](../../scripts/figures/run_fig2_coordinate_panels.sbatch) requested
one CPU, 8 GB RAM and 90 hours on partition `cpu`, account/QOS `nslab`.
The first job rendered the PDFs but stopped because `pdffonts` was absent from
the compute-node PATH. The retry uses the installed Poppler directory explicitly;
no dependency was installed. The statistical analyses were not rerun.

| Job | Panels | State | Elapsed | Memory and errors |
|---|---|---|---|---|
| 22219668 | 2/2 rendered; checks incomplete | Failed, exit 127:0 | 00:00:09 | 8 GB requested; missing `pdffonts`, resolved in retry |
| 22219710 | 2/2 ████████, checked and installed | Completed, exit 0:0 | 00:00:04 | 8 GB requested; peak RSS not captured for this short job; no unresolved errors |

Running: 0; pending: 0; successful jobs: 1; failed attempts: 1, resolved.
No downstream job remains. The successful run is
`figures/candidates/fig2-coordinate-panels-20261013T220927Z_22219710/` on
`ne1dc4-068.nygenome.org`. Its
[completion record](../../figures/candidates/fig2-coordinate-panels-20261013T220927Z_22219710/COMPLETE.txt)
and [Figure 2D check](../../figures/candidates/fig2-coordinate-panels-20261013T220927Z_22219710/fig2d_checks.log)
record success; all **23 input/snapshot checksum entries** were unchanged.
The [existing Figure 2D check](../../scripts/figures/validate_fig2d_pip_architecture.py)
was updated to the corrected eight fractions, 19/228 consequence clusters and
6.026671 ratio, and rejects displayed diagnostic significance labels.

Both PDFs have one page, exact dimensions (D 1.85×2.10 inches; companion
5.50×2.70), embedded Helvetica, no Type 3 fonts, and passing Ghostscript checks.
`pdftotext` values were compared to the source tables, then both panels were
rendered to temporary inspection PNGs and visually checked for clipping and
overlap. An independent reviewer confirmed the installed values, inferential
units, interval table and manuscript summaries.

Twelve generated PDF/source/caption/interval/measurement files were copied
exactly to the working Figure 2 panel directory, along with the shared companion
source table. Their bytes match the successful run. The prior D files and index
are preserved under
`figures/main/fig2_genetics/panels/_archive/2026-10-13_before_coordinate_correction_22219710/`.
Only D changed among the existing main-panel PDF hashes; the companion was added
outside A-I. The D entry in the existing panel index now identifies its new PDF.
Named historical analysis releases were not overwritten.

Current adoption and values are synchronized in [PAPER](../PAPER.md),
[RESULTS](../RESULTS.md#adopted-corrected-coordinate-results-2026-10-13),
[STATUS](../STATUS.md), [legends](../manuscript/05_figure_legends.md), the
[Figure 2 manuscript](../manuscript/draft/fig2.md), [README](../README.md),
[ROADMAP](../ROADMAP.md) and the [Figure 2 directory guide](../../figures/main/fig2_genetics/README.md).
The renderer, focused Figure 2D check and new rendering launcher are the only
analysis-code changes in this adoption step. No commit, branch, push or pull
request was created.

Final synchronization checks passed: `python3 scripts/manuscript/validate_resource_scope.py`
returned RESOURCE_SCOPE_VALIDATION PASS (37 documents, seven figure roles),
all 26 added local link targets across nine documents resolve, and source diffs
passed `git diff --check`. The installed D index matches its PDF hash, and the
renderer, helpers, Figure 2D check and launcher match the successfully executed
copies. All initial working-tree entries remain preserved.

## Remaining Resource corrections, 2026-10-14

The user authorized all remaining corrections and necessary jobs with parallel
agents. Five agents independently handled bulk/context comparisons, current
panels, genetics, identity-qualified agreement and the selected r8 Catalog.
The primary agent reconciled the scientific claims and inspected the resulting
Figure 3/4 assemblies. The initial working tree and authority copies are
recorded in `/scratch/resource-remaining-20261014-bm9cu4g5/`; unrelated changes
are preserved. No branch, commit, push, external portal publication or guide
design was performed.

### Corrected biological and statistical conclusions

- **All five per-cohort RNA contrasts are now represented.** GSE126848's
  selected 30 genetic clusters (23 control, seven disease, all male) support a
  full-rank disease contrast after omitting constant sex. The fixed pooled
  filtering, normalization, quality-weighted limma-voom method and complete
  24,115-gene BH family yield 6,147 canonical calls. Native fragment counts
  match the selected DGE exactly. A separate group-means parameterization
  reproduces every call (maximum absolute effect difference 3.06×10⁻¹⁴,
  BH-q difference 9.29×10⁻¹⁴). Other-cohort effects remain unchanged over
  96,460 rows. The source screen excluded 24/31 disease records, all ten
  women included, and 2/25 controls for insufficient identity sites; this
  differential eligibility limits representativeness. The five-cohort union
  is 13,761 genes and 20,054 memberships: 9,065/3,462/922/261/51 genes called
  in one through five cohorts. Among the 1,854 pooled calls, 1,564 agree in
  GSE126848 direction (84.36%); 764 are also called in that cohort.
- **Figure 5B/C now use 774-unit RNA estimates.** The unchanged July
  10,044-gene display universe gives 3,054 measured genes. Protein-significant
  211 divides into 188 concordant and 23 discordant; 2,843 are nonsignificant.
  Spearman rho=0.457835 and direction agreement=65.9136%, printed as 66%.
  The fixed 25-protein view has nine RNA and 23 protein calls. CD74 enters
  and HLA-DQA1 leaves the RNA set relative to September. HLA-DQA1 now uses
  ENSG00000196735.14 (RNA log2FC=0.379801, padj=0.024969; protein
  log2FC=1.393, padj=0.07964); both markers are open. The rows remain
  selection-conditioned, not independent protein validation.
- **The class-specificity difference attenuates.** On the native 14,704-gene
  joint universe, expression matching gives 281 hepatocyte-specificity pairs
  and 73 RNA-effect-adjusted protein pairs. Mean genetic-minus-state effects
  are −0.265515 log2 and −0.015601. The prespecified directional prediction
  fails. Post hoc paired sign-flip nominal p-values are 0.01482985 and
  0.487465 (100,000 draws, seed 20260907), without correcting shared loci or
  gene dependence. The retained relabel/rematch p/q values are explicitly
  invalid anti-conservative diagnostics, not confirmatory significance.
  The literal prespecified negative control is separately recorded using
  1,000 within-pair label swaps per declared outcome; simulated means
  −0.001071/−0.000441 have Monte Carlo SE 0.003442/0.000721, against expected
  zero. No new acceptance threshold or improved directional verdict was added.
- **N4 does not qualify the repeat-biopsy claim.** All 164 GSE193066 library
  identities remain unqualified. Historical different-visit n=22 point
  estimates remain retrospective descriptions, but exact-stage independent
  comparison capacity is insufficient at stages 2–4. The null is untestable,
  not rescued with replacement or substituted stages. Historical 87 random-null
  and 69 bootstrap-tail stage-BH calls are withdrawn. The separate 28 untyped
  pairs give zero calls against their eligible exact-stage null (117-program
  BH family). The distinct historical 58-pair model ICC and SD=0.6433 ceiling
  calibration are qualified as diagnostics/scenarios, not primary agreement
  or calibrated measurement-error estimates.
- **All 66 fixed twelve-trait genetic correlations now have estimates.** Eight
  missing PDFF and 27 additional non-PDFF pairs were repaired. Installed LDSC filters
  incompatible alleles after caching its alignment vector. Candidate inputs
  remove duplicate rsIDs and use the exact installed `MATCH_ALLELES` rule,
  retaining native Z/signs; no installed LDSC code or original summary
  statistics changed. All repaired estimates pass finite/positive-heritability
  checks and all 256 allele combinations reproduce the installed rule.
  FinnGen cirrhosis remains imprecise because its estimated heritability is
  small. P-values are nominal; neither the repaired subset nor complete
  matrix was assigned a new FDR method. The terminal matrix preserves 39
  existing successful estimates and exposes native historical versus repaired
  preprocessing in each row; it is not a uniform new 66-pair refit. Six
  unconstrained correlations above one involve weakly heritable Ghodsian
  NAFLD or FinnGen cirrhosis; values and SEs are retained, without clipping or
  interpreting a true genetic correlation above one.
- **Current r8 Catalog classes are unchanged.** The exact reviewed copy
  contains 28,233 genes and 56,466 gene-evidence rows, with classes
  400/26/1,129/13,127/13,551 and zero strict tested-negative calls. All 117
  programs and 6,633 memberships are preserved. The reviewer is honestly
  recorded as “Codex remaining_catalog scientific review.” Twenty independent
  source checks and five rebuilt images pass. The 28,233-record molecular-object
  companion adds no authoritative target link or PIP mass. Raw symbol-grouped
  14,720, native stable-ID 14,704 and Catalog-adjudicated 14,682 joint universes
  are explicitly distinct. The earlier 24/8 mismatch was not reproduced:
  direct source symbols differ 16/0; Catalog/native IDs differ 6/28, accounted
  for by six resolved aliases, 21 quarantined native identities and seven
  genetically untestable identities. The 855 held annotations remain untestable.

### Working figure corrections and source checks

Current Figure 3 shows 774 units, a 364-unit exact-stage/NAS grid, the restored
five-cohort membership census, 392 endpoint units, complete stage BH counts
without an effect-size filter, 54 genes/216 adjacent-stage effects, and 3,776
finite RNA–SuSiE pairs. Figure 4D retains all 117 lineage assignments. Figure
4E uses the selected four-cohort 457-record/112-program arm (16 fibrosis and
one NAS call, BH over 224 tests), including 93 GSE174478 source records outside
the pooled fit. Figure 4H shows rho=0.73 and 37/43 same direction with one
uniform-null outline; 18 strict maxT calls remain a distinct supplementary
criterion. All five individual Figure 6 panels use r8. Figure 4A is a new
180×29-point vector schematic of descriptive atlas census and fixed programs.

Two NMF stage axes, k4-P2 and k6-P6, fail the unchanged reconstruction rule and
are marked unavailable without stage trajectories. The context-source screen
also found EMT Hallmark membership coverage 75.5% and RELA/SMAD3 original L1
coverage 74.3%/75.8%, below the fixed 80% rule. Their false near-zero curves
are removed; the facets explicitly say untestable. Five of six selected
Hallmarks and two of four selected regulons are testable. Methods, captions
and future producers retain those states.

The main/supplementary panel branch installed 46 PDFs plus native-size Figure
3/4 review assemblies, with backups in
`figures/candidates/pre-774-install-20261014T140418Z/` and
`figures/candidates/pre-774-supplement-install-20261014T140856Z/`.
The assemblies are [Figure 3](../../figures/main/fig3_bulk_transcriptomics/review_assembly/Figure3.pdf)
and [Figure 4](../../figures/main/fig4_singlecell_programs/review_assembly/Figure4.pdf).
Source values, page boxes, embedded Helvetica/6-point text, absence of Type 3
fonts, retained assembly text and visual appearance were checked. PDF-only
publication output and fixed program definitions are preserved. Figure 6's
five individual panels additionally have a new
[native-size review assembly](../../figures/main/fig6_gene_catalog/review_assembly/Figure6.pdf)
(468×970 pt), checked without scaling the panels. Final publication layout
remains separate; the historical two-panel raster proof is not substituted.

The verified r2 panel update replaces Figure 1B/C, 2A/B/C/E and 5D, while
render-identical 2F/G/H/I remain. Figure 1B uses native 14,704/426/1,155/26
counts and explicitly retains seven SuSiE-untestable genes as not positive,
without calling them eligible negatives. Figure 5's vector assembly now uses
current 5B/C/D. The separate target-gene biotype companion has 424 protein-coding
and 12 lncRNA genes, all resolved from native liver-eQTL tables; gene-level
QTL support does not identify isoform biotype or establish coverage enrichment.

Independent source checking caught an additional exact-locus error before
installation: the current GNMT direct-NAFLD ABF anchor is rs2296804 at
GRCh37 6:42931261, SNP.PP.H4=0.258990, not the historical rs2296805 anchor.
The native LD reference gives rs2296804/rs2296805 r²=0.998001. Nonpalindromic
rs2296805 is a proxy whose lower-GNMT-expression allele aligns with higher
disease risk. The palindromic rs2296804 anchor and rs11752813 have unresolved
strand direction; reported allele labels alone are not independent harmonization.
CYP2C19's current ABF anchor is 10:96292988, SNP.PP.H4=0.044813. Current prose
and exact source rows now agree. The source-owned liver-biopsy protein and
reporter evidence remains context, not a new MASLD finding or resolved causal
variant. The native-source assertion prevented installing the failed first
display candidate, job 22236870.

Two remaining `Liver_Atlas` crosswalk aliases were corrected to Guilliams
GSE192742. LangCell's checkpoint-bound exposure becomes unknown; historical
December-2021 scBERT is clean by chronology, while its mutable registered object
remains unknown. Current Figure 7 values are unchanged because these models
do not enter the working exposure grid. Source checks and before/after tables
are in `Analysis/MASLD_Model_Benchmark/executions/liver-atlas-alias-correction-20261014/`.

### Executed compute and remaining limits

All jobs used CPU, account/QOS `nslab`, no arrays and 90-hour requests.

| Work | Successful job(s) | Allocation | Completion | Recorded peak memory |
|---|---|---|---|---|
| Bulk context, five-cohort source and class comparison | 22236311 | 4 CPUs, 48 GB | 4/4 ████████; 7m50s | 709,608 K |
| Native GSE126848 independent refit | 22236559 | 1 CPU, 16 GB | 5/5 ████████; 48s | 554,312 K |
| Literal paired-label negative controls | 22236840 | 1 CPU, 4 GB | 2/2 ████████; 7s | not captured |
| PDFF correlations and independent check | 22236252, 22236390 | 2 CPUs, 24 GB; check 1 CPU, 8 GB | 8/8 ████████; 8m12s | 1,770,356 K |
| Remaining native matrix pairs and full table check | 22237000, 22237084 | 2 CPUs, 24 GB; check 1 CPU, 4 GB | 27/27 and 66/66 ████████; 28m13s | 1.75 GiB |
| COLOC r2 independent completion check | 22236309 | 1 CPU, 8 GB | 50/50 ████████; 11s | not captured |
| N4 eligibility correction | 22236285 | 8 CPUs, 32 GB | 1/1 ████████; 2m48s | 1,342,748 K |
| Catalog companion and interface | 22236527 | 2 CPUs, 24 GB | 20/20 ████████; 2m52s | 1,683,204 K |
| Exact Catalog review and read-only check | 22236566 | 2 CPUs, 24 GB | 1/1 ████████; 3m53s | 1,131,212 K |
| Catalog identity reconciliation | 22236628 | 1 CPU, 4 GB | 1/1 ████████; 2s | not captured |
| Figure 6A/C/E rendering | 22236670 | 2 CPUs, 8 GB | 3/3 ████████; 2s | not captured (batch reports 0) |
| Figure 4A | 22236599 | 1 CPU, 4 GB | 1/1 ████████; 4s | not captured |
| r2 interface, target biotypes, native examples and Figure 5 assembly | 22236878 | 2 CPUs, 16 GB | complete; 2m06s | 59,632 K |
| Figure 1B observability label | 22237222 | 1 CPU, 4 GB | 1/1 ████████; 5s | not captured (batch reports 0) |
| Figure 6 native-size review assembly | 22237086 | 1 CPU, 8 GB | 5/5 ████████; 6s | not captured |
| S3J bootstrap / final Figure 3/4 assembly | 22236554 / 22236667 | 4 CPUs, 16 GB each | complete | assembly 440,912 K |

Intermediate failures are preserved rather than erased: dtype inference and
UI CLI serialization, assumptions caught in Catalog identity checks, missing
PDF-tool PATH entries, a display subset expression, native panel-height
assertions and inherited text-size/placeholder overlays. The corrected final
sources and outputs passed their focused checks. There is no unresolved
inferential-fit failure in these completed analyses.

Sources: [bulk/context outputs](../../figures/candidates/bulk-context-774-20261014T140904Z_22236311/),
[N4 outputs](../../Analysis/Hypotheses/N4-eligibility-corrected-20261014T140730Z/),
[PDFF results](../../GWAS/ldsc/rg_pdff_compatible_20261014T140601Z_22236252/pairwise_estimates.tsv),
[all 66 terminal pairs](../../GWAS/ldsc/rg_missing_compatible_20261014T150911Z_22237000/rg_12trait_66pair_terminal.tsv),
[r2 completion](../../GWAS/finemapping/results/coloc_r2_adoption_review_20261014_22236252/review.json),
[reviewed r8](../../figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/PASS06_VALIDATED),
and [local Catalog access](../../figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/access_v2/index.html).

Authoritative HMSMA donor/histology/registration joins and qualified GSE193066
identities remain unavailable. No rerun can supply those identities. The
final publication layout, one coherent read-only output bundle and overall
approval remain the ROADMAP checklist; individual panel completion and local
Catalog review do not assert those final states.

### Additional original-source diagnosis checks and dependent refits

The final source check found two additional category errors. The original
[Filliol Nature study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9949942/)
identifies GSE174748's disease explants SRR14602057/58 as cirrhotic; generic
ENA NAFLD titles did not establish Steatosis. The unchanged three-stage
program rule excludes them. The exact
[Guilliams Table S6](https://ars.els-cdn.com/content/image/1-s2.0-S0092867421014811-mmc6.xlsx)
joins all nineteen Liver Atlas donors: four Healthy, three Steatosis, one
Cirrhosis, ten other fibrosis records and one unknown. The latter eleven stay
Unresolved without sufficient clinical etiology/activity; all nineteen retain
their stage exclusion. GSE185477 Healthy labels are verified by the
[Andrews source paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8948611/).
Exact donor/run joins, originals, primary evidence and hashes are in
[the source correction](../../Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/source/README.md).
The descriptive 102-donor counts become 24/14/36/17/11 across Healthy,
Steatosis, Steatohepatitis, Cirrhosis and Unresolved.

The fixed-score refit has maximum primary n=62 and ambient common n=60.
Only Ductular injury is primary-selected: beta=0.402086, BH q=0.037633,
HC3 q=0.155826. Stromal ECM loses selection: beta=0.351525, BH q=0.085458,
HC3 q=0.165240. Zero programs pass HC3. Independent NumPy reconstruction
matches every estimable R primary OLS/classical/HC3 result; all twelve ambient
117-program BH families reproduce. Scores and all 117 memberships/weights
remain unchanged. Equal-run, cohort and leave-one-dataset-out fits preserve
their original rules. In the 60-donor common ambient arm, ECM raw/corrected/
delta is 0.350898/0.127891/-0.223006 (delta BH q=0.478071); BICC1 is
0.384613/0.401035/0.016422 (delta q=1). Neither paired change passes BH.
ECM's selection-conditioned attenuation classification is no longer applicable;
BICC1 remains indeterminate. The full 62-donor passthrough arm stays a
sensitivity, not a fully ambient-corrected primary result.

The separate 26-donor content-lineage meta analysis uses only
GSE202379/GSE244832 and is unchanged, as are ten paired descriptive lineage
contrasts. Main 4B's Healthy/MASH density excludes the corrected samples and
is unchanged. Corrected 24-row 4C, all-117 S4C and the Figure 4 vector assembly
are installed with original copies preserved; only Ductular injury carries a
selection asterisk. Complete results and numerical proof are in the
[working companion](../../figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/).

H34's saved-score four-stage refit retains sex calls, Y-negative fractions,
fifty autosomal control sets, cache hashes and within-donor comparisons.
Macrophage/T-cell effects become 1.487719/0.486050 percentage points per stage,
HC3 intervals [-1.413582, 4.389019]/[-0.870552, 1.842651], both BH q=0.482539
in the two-immune-type family: zero calls, candidate null retained. Actual
fitted n=18/22 across three datasets; eligible 19/23 includes the singleton
GSE185477 donor already dropped by the original HC3 rule in both old and new
fits. Earlier STATUS overstated fitted n. Source results are in
[the corrected H34 candidate](../../Analysis/Hypotheses/H34-stage-source-corrected-20261014_22237596/CHECKED.txt).

The 26-test program–communication family and complete LR stage screen were
also refitted on unchanged saved scores and their original families. The old
LR metadata join had already dropped 38 unmatched Liver Atlas run aliases
from 117 score identifiers, leaving 79; it omitted source-stage exclusions
for ten GSE136103 donors. The corrected fit restores those native exclusions
and uses 69 source donor identifiers. The separate 125-unit lineage census
contained 42 Liver Atlas run IDs; it was not the actual fitted sample.
All 55,402 original LR/cell-pair groups have terminal dispositions: 349 fitted,
55,053 insufficient stage counts. The 1,047 coefficients retain within-cell-
pair/term BH (27 families) and across-cell-pair Bonferroni (nine families).
Six of the fixed thirteen displayed pairs are testable, all adjusted q=1;
seven are untestable, not negative evidence. The original stage-count and
native mixed-model/fallback rules are preserved.

The fixed 26-test program screen has 22 estimable and four untestable tests.
Two Ductular injury associations pass both BH families, each n=58 across four
datasets: CDH1→PTPRM positive (beta=0.495716, q=3.627665×10^-6,
HC3 q=1.167371×10^-10) and ST6GAL1→EGFR negative (beta=-0.298712,
q=0.006914, HC3 q=0.048800). All four estimable leave-one-dataset-out signs
agree for each; ranges are 0.003372–0.698344 and -0.566162 to -0.069645.
Original 26-fit coefficients and classical/HC3 p-values reproduced before
correction. These are same-atlas conditional RNA covariations, not independent
validation, secretion, receptor activation or antagonism. Exact sources are
in the [qualified CCC candidate](../../Analysis/SingleCell/candidates/source-qualified-CCC-20261014_22237840/CHECKED.txt).

| Additional work | Successful job | Allocation | Completion | Peak memory |
|---|---|---|---|---|
| Primary and ambient diagnosis refits | 22237556 | 2 CPUs, 8 GB, 90 h | 117/117 primary; 234/234 ambient rows ████████; 32s | 80,056 K |
| Independent numerical proof and 4C/S4C | 22237573 | 2 CPUs, 8 GB, 90 h | complete; 20s | not captured |
| Figure 4 assembly refresh | 22237764 | 1 CPU, 8 GB, 90 h | 8/8 ████████; 54s | 308,428 K |
| H34 source-stage refit | 22237596 | 1 CPU, 8 GB, 90 h | 2/2 immune models ████████; 5m07s | 3.86 GiB |
| Complete LR stage families and fixed 26-program screen | 22237840 | 1 CPU, 8 GB, 90 h | 55,402/55,402 group dispositions ████████; 40s | 1.08 GiB |

A compact [Figure 6 publication-layout candidate](../../figures/candidates/fig6-compact-publication-20261014T160441Z_22237441/Figure6.pdf)
is prepared at 468×600 pt, with native 6-point text, unchanged values and
all five panel identities. Job 22237441 completed in 22 seconds; count totals,
complete text, embedded fonts and visual checks pass. It is a concrete proposal,
separate from the installed tall review assembly and overall publication approval.

The first read-only snapshot,
`figures/candidates/resource-review-20261014T155914Z_22237418/`, contains
324 copied files and 43 current-at-snapshot main PDFs; all copies are read-only,
all copy checksums pass and all 55 Catalog input hashes match. Its job completed
in 13 seconds after a first attempt stopped on two incorrect sidecar paths.
It predates the diagnosis correction and must not supply final 4C values.
Figure 1A has no indexed PDF; an unindexed author PNG was not substituted.
The final snapshot follows the dependent refits and document synchronization.

The current 117-program annotations were also corrected without another fit.
Ordinary nonsignificance does not establish an adequately powered negative
result: every current strict `tested_negative` flag is false, with
`tested_nonsignificant` reported separately. BICC1 is the sole primary-selected
program; ECM's earlier selection annotation is superseded. Historical native
flags are retained only as explicitly named source fields. Numerical effects,
q-values, scores, memberships and weights are unchanged. Source stability
describes the original cohort/leave-one-dataset-out direction rule and does
not assert statistical significance.

Supplementary S4D/E now use the qualified communication fits above. The final
render, job 22237968, passed embedded-font, text and visual checks; previous
artwork and source tables are preserved. Fresh proofs of the actual current
Figure 5B/C PDFs (job 22238156) confirm 188 concordant/23 discordant proteins,
66% rounded direction agreement, open HLA-DQA1 RNA/protein markers and the
canonical RNA effect-size floor. Earlier temporary previews showing
185/26/64% are historical and cannot verify these current PDFs.

H17 also consumes the corrected source stages. Its earlier STATUS
estimates are explicitly historical and withheld pending the full native
resampling and control checks. H17 retains the original case definition and
both-groups-per-dataset rule: excluding the two cirrhotic explants also removes
that dataset's two controls, giving 55 donors (9 controls, 46 cases) in two
datasets. N3 completed its source-qualified raw and ambient refits with all
70 donors and the original four-stage model, moving the two disease records
from Steatosis to Cirrhosis (14/11/34/11). Both 1,000-permutation calibrations
fail: mean nominal-positive fractions are 0.078889 raw and 0.078125 corrected.
The prespecified empirical-BH fallback yields nine raw calls among 45 eligible
targets and zero corrected calls among 40; the fixed roster contains 51 targets.
Ordinary BH counts of 12 and six are diagnostic, not calibrated primary calls.
Both fixed positive-control composites pass. All 208 included libraries and
70 donor records have native correction evidence and no passthrough; the
historical nine-target ambient diagnostic roster remains fixed separately.
The earlier conclusion that LEPR is robust is superseded: corrected LEPR
beta is 0.5440, ordinary BH q=0.01486 but primary empirical BH q=0.07421.
No detected corrected calls does not establish absent effects or powered
negative evidence. The nominal contamination-stage coefficient is 0.1773
(95% CI 0.04938–0.3052, p=0.007348). Native attenuation labels remain
ordinary-BH diagnostic summaries, not primary evidence calls. Job 22237970
finished both fits but failed in the diagnostic summary because its read-only
per-library link was omitted; repair/verification job 22240531 completed
without new resampling. All thirteen success checks in
`Analysis/Hypotheses/N3-source-diagnosis-20261014/verification.json` pass.
Original
expression matrices, target/program definitions, seeds and resampling counts
are preserved. No preliminary point estimate substitutes for the terminal
control-qualified result.

### Current model-source interpretation and complete figure proofs

The evaluated TranscriptFormer TF-Sapiens checkpoint is independently bound
to all thirty prediction receipts and the 50,000-cell embedding receipt.
Its complete 644-UUID training inventory resolves five clean Resource studies
and two training-seen studies, GSE185477 and Guilliams Liver_Atlas. The old
seven-unresolved display is superseded in working 7C/F. Geneformer V2 alone
has zero qualified clean studies among the displayed models. Its GSE174748
tile was already correctly unresolved; a preliminary claim that the tile was
magenta was a misread of the neighboring column and is rejected by the actual
rendered panel. LangCell and scBERT are not rows in these panels.

All existing model scores, gains and intervals are exactly preserved.
TranscriptFormer's pooled macro F1=0.9388747 maps five broad cell classes
and includes the two training-seen studies; it is not disease-stage or
clean-only performance. The original clean-subset estimand is now computable
on five studies and 80 donor identifiers. The original UCE subset reproduces
to 1e-12 before reaggregating TranscriptFormer from the same 10,000 saved
paired donor-bootstrap draws, seed 20260824. Relative to the liver encoder,
the linear-head macro-F1 difference is -0.004170 [-0.008732, -0.000436] and
the two-layer-head difference is -0.003490 [-0.008311, -0.000181]. These
intervals condition on fixed models, five studies, observed cells and
evaluable class coverage, not a population of new studies. Donors, rather
than cells, are resampled within study. All fifteen liver-encoder seed/fold
receipts independently show disjoint training/query studies, no held-query
training/adaptation and training-only feature selection. The comparator is
therefore appropriately held out. No model is refitted or new significance
claim introduced. Sources and checks are in the
[model-source correction](../../Analysis/MASLD_Model_Benchmark/executions/fig7-source-exposure-correction-20261014/README.md).
Render/check job 22238704 completed in eleven seconds, using the existing
project environment after an initial system-Python dependency failure.

Review-only assemblies now exist for all seven figures. Figure 1's missing
A/D schematics are source-backed proposals; its 35-stratum count independently
matches fifteen Tier-1 and twenty Tier-2 rows of the exact phenotype registry.
Bulk and cell counts retain their separate source units. Its font-embedded
[A–D review](../../figures/candidates/fig1-embedded-review-proposal-22238522/Figure1.pdf)
preserves existing B/C values and native sizes at 468×704 pt. Figure 2's
[A–I review](../../figures/candidates/fig2-native-review-proposal-20261014_22238526/Figure2_A-I_native_review_PROPOSAL.pdf)
is 7.2×9.2 inches; all 1,601 native glyphs retain font/geometry, all nine indexed
hashes match and all seventeen font resources are embedded. Figure 7's
[six-panel review](../../figures/candidates/figure7-embedded-native-review-proposal-22239269/Figure7_embedded_native_review_PROPOSAL.pdf)
includes corrected A/C/F at 758×731 pt. All native words are preserved; maximum
word-boundary change from embedding is 0.1084 pt, with all fonts embedded.
Native artwork is translated without rescaling. Working assemblies for
Figures 3–6 and the compact Figure 6 proposal remain as recorded above.
Final journal sizing and selection of proposed layouts remain human decisions;
these reviews do not assert publication approval.

Figure 7A's September 1 census now labels 106 model entries rather than
families. Its 104 source-audit files contain 89 checkpoint findings and 37
no-finding placeholders, so the bar's 126 entries are not 126 models.
The counts are unchanged. Render job 22239160 and final assembly job 22239269
completed; the final six-panel proof and all six current source hashes agree.
The final proposal SHA-256 is
`53de058e34514a6f7006232e1b5fced2b124eee7c9f60469f5f1cadb56326f86`.

The refreshed Figure 1–6 claim-context companion is
`figures/candidates/figure-claim-context-20261014T185701Z_22240302/`.
All 147 source files and 320,281 lexical rows reproduce their source values.
There are 11,727 source rows with unique gene matches and 392 explicitly unmatched rows,
and all 11,023 program rows join to fixed membership identities. The 117-program
registry has zero strict negative calls, 116 ordinary nonsignificant programs,
BICC1 alone primary-selected and zero HC3-selected. Current stage context is
separate from the fourteen historical multimodal contexts in named r8.
Job 22240302 completed in 25 seconds with all current source metadata hashes;
the earlier spacing-corrected candidate 22239479 is preserved. All native
values are unchanged. All source, unit, family and membership
checks pass. Neither unresolved target links nor ordinary nonsignificance
become evidence of biological absence.

The completed r2 method-availability count additionally corrects older Methods
denominators. In the full 50-study master, 21,799 SuSiE records are finite
across 4,869 identified genes; the 35-study primary subset has 12,892 across
3,918 genes. The raw ABF denominators are 948,503 and 664,022. These raw counts
include 600 and 420 chromosome-specific anonymous rows with empty gene IDs and
symbols; none has finite SuSiE, and all are excluded from gene counts.
Identified ABF record counts are 947,903 and 663,602. The apparent 550 extra
collapsed study–gene keys arise entirely from twelve anonymous chromosome
records per study, not duplicated identified gene observations. The native
study–chromosome–gene key is unique. Numeric finiteness was checked directly;
20,517 carried negative rows were not freshly rerun under the shared-posterior
rule. The 25 failed shared-posterior pairs affect 21 genes, four of which retain
other positive support; the primary equivalents are eleven pairs, nine genes
and two retaining positive support. Partial untestability is distinct from
complete lack of finite output. The
[final availability and identifier checks](../../GWAS/finemapping/results/coloc_r2_adoption_review_20261014_22236252/method_availability_22239724/identifier_reconciliation.json)
preserve all original master and registry hashes. Jobs 22239649 and 22239724
completed in 35 and 42 seconds at approximately 340 and 728 MiB peak memory.
Current Methods and RESULTS separate these values from July/August history.

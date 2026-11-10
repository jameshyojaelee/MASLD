# Scientific review of r9 main results and next analyses, 2026-10-26

Scope: `docs/manuscript/MANUSCRIPT_r9_assembled.md`, release 2026-10-21-r9.
Line numbers refer to that file. This is a technical-layer review. It changes no
authority, figure, release or Gene Catalog entry.

Method. Nine area reviewers re-derived printed numbers from the release source
tables. An adversarial pass then tried to refute every finding, and a
feasibility pass checked each proposed analysis against the data on disk, prior
work in `Analysis/Hypotheses`, and the literature. A final pass looked for
problems the first two missed. I re-checked the three high-severity items and
the Figure 4 leave-one-dataset-out numbers directly from the files cited below.
Items marked [Not re-derived] come from the final pass only.

## 1. Verdict

Almost every printed count re-derives exactly: 436, 996, 1,854 (1,333 up, 521
down), 26 against 21.2 expected, the five Catalog classes, 0 of 117, 18/33/0/1
systems. The problems are in what several numbers are said to show, not in the
counts.

Three findings change a claim as written:

1. Figure 7B and Methods say sequence baselines were "not run". They were run,
   and against them the HyenaDNA gain does not survive the paper's own
   Bonferroni family.
2. The prespecified matched-gene test of the continuum missed in both primary
   cohorts and is not reported.
3. Both focal single-cell programs are gene sets defined by one dataset,
   GSE244832. The cell-cell communication associations and the fibroblast
   Stromal ECM contrast also rest on that dataset.

The paper's positive statements are mostly descriptive or known. The overlap
null (Figure 1B) is robust. The GNMT, MAT1A and CYP2C19 decreases reproduce but
were reported before. The strongest candidate for a new biological result is
per-hepatocyte GNMT loss (DIG1, exploratory, not in r9). It is cheap to extend
from two cohorts to four.

## 2. Strength of each main result

| Result (line) | Printed | Strength | Reason |
|---|---|---|---|
| 436 colocalized genes (58) | SuSiE PP.H4 > 0.5, p12 5e-6 | Moderate | 314 remain at p12 1e-6; 303 also pass the LD check. 272 rest on one study: 175 had no same-trait second test and 97 were tested and failed. MVP and UKBB European GWAS of the same enzyme replicate 49 to 63% of each other's positives. |
| Enzyme traits dominate (74) | 423 of 436 | Strong | Per genome-wide lead, yields are equal: direct 28/209 = 0.134, enzyme 739/4,945 = 0.149. The dominance comes from locus count. |
| Noncoding posterior majority (Fig 2D) | 66.4% direct, 89.8% enzyme | Strong | Holds in both classes under either weighting. Enzyme share may be inflated by saturated non-European credible sets [Not re-derived]. |
| Direct loci more coding than enzyme loci | 29.65% vs 6.67% | Not supported | The text already says no ancestry stratum is calibrated. |
| caQTL shares 191 of 581 lead signals (104) | 191/581 | Moderate | Re-derived. Currin already reported at least 45 of the 191, and the text says so. |
| GNMT worked locus (82) | ABF PP.H4 0.706 | Weak | NAFLD p 2.6e-4. In all five large European enzyme GWAS the local association favours a different variant (ABF PP.H3 0.94 to 1.00; best SuSiE PP.H4 0.128). The text already calls it sub-threshold. |
| 1,854 disease-state calls (112) | BH over 24,115, abs log2FC > 0.5 | Moderate | Direction is strong. Membership depends on GSE213621 (median 54% of per-gene information; 409 calls carried by it alone). Heterogeneity is large (section 3, M3). |
| Fig 3C mixed sign | 13.5% without GSE126848 | Strong, descriptive | 325 of 346 mixed genes involve GSE135251. Without it the share is 2.2%. |
| F2 to F4 vs F0 | 4,266 / 6,190 / 7,043 | Moderate | About half the F2 calls survive library-property adjustment (Supplementary Note 7). The positive control for that adjustment failed. |
| Fig 3D strict-control endpoint | 1,755 of 1,854 | Weak | 68 of the 113 strict controls are GSE213621 controls. The endpoint reuses most of the same control arm. |
| GTEx fat direction (Fig 3) | rho 0.638 | Moderate | 323 of 326 genes with fat q < 0.05 keep the pooled sign, in 222 independent donors. |
| Overlap test (Fig 1B) | 26 vs 21.2, OR 1.29, p 0.28 | Strong as a null | Holds at p12 1e-6 and after dropping LD-failing genes. The CI (0.84 to 2.01) allows a two-fold enrichment. 861 of the 935 untested disease-state genes have no liver-eQTL credible set and could not be tested. |
| 26 two-axis genes | 26 | Weak | All anchored by ALT, AST or GGT. None reaches PP.H4 > 0.5 in a NAFLD, NASH or PDFF study. Directions: 11 agree, 7 oppose (P 0.48). IFI27L2 and SPP1 sit on neighbor coding variants. |
| Catalog classes (Fig 6) | 400 / 26 / 1,129 / 13,127 / 13,551 | Strong | Exact. The route field has one value per class, so "a next experiment for each gene" describes classes, not genes. |
| 0 of 117 programs track stage (178) | smallest q 0.323 | Moderate | Correct, but underpowered. A single BH call needs about 0.8 to 1.1 residual SD per stage step. Post hoc, 18 programs have p < 0.05 against 5.85 expected (joint p 0.016). |
| Stromal ECM attenuated by ambient correction (180) | 64% | Moderate, understated | decontX assigns 81% of Stromal ECM member-gene counts in hepatocytes to ambient RNA (3rd of 106 programs). |
| Fibroblast Stromal ECM contrast (184) | beta 0.329, 26 donors, q 0.0106 | Weak | Rests on GSE244832 (section 3, H3). |
| CDH1 to PTPRM and ST6GAL1 to EGFR (186) | supported in 58 donors | Not supported | Without GSE244832: beta 0.496 becomes 0.0034, and -0.299 becomes -0.070. |
| Bulk fibrosis program map (Fig 4E) | 0/224 calls | Moderate, understated | Content lineage sorts the effects (section 4). |
| Continuum program coefficients (202) | 0.760 and 0.707 | Moderate | Composition explains 80% (GSE162694) and 82% (GSE213621) of within-stage continuum variance. Stromal ECM keeps about 0.15 SD per SD after portal-content adjustment; Ductular injury keeps partial r 0.36. No random-axis null exists. |
| 43 systems | 18 / 33 / 0 / 1 | Strong as stated in line 212 | The heading at line 206 claims more than the result. |
| Continuum adds information beyond stage | implied | Weak | No outcome, within-person or treatment test. The prespecified fibrosis test missed (section 3, H2). |
| RNA-protein agreement (Fig 5B) | rho 0.46 | Moderate | The protein assay's own between-batch agreement is rho 0.34, so 0.46 is near the ceiling. |
| Chromatin enrichment (Fig 5D) | 2.5 to 3.1 times | Weak | Matched expectation gives 1.8 to 2.4. Only hepatocyte-specific peaks are enriched (1.77). Stellate, macrophage and cholangiocyte-specific peaks are depleted (0.42 to 0.62). |
| Spatial Stromal ECM (Fig 5F) | "supported" | Not supported | 68 and 74 of 114 programs pass the same test; Stromal ECM ranks 42nd and 29th. Donors are not MASLD. |
| HyenaDNA reporter gain (Fig 7B) | 0.082, Bonferroni lower 0.009 | Weak | Holds only against allele identity (section 3, H1). |
| Fig 7F encoder comparison | +0.0154 | Weak | The sign flips with the estimand (section 3, M8). |

## 3. Corrections, ranked by effect on a claim

### High

**H1. Figure 7B, line 562, Methods line 1735 and PAPER.md claim 8 say sequence
baselines were not run. They were.**
Source: `Analysis/MASLD_Model_Benchmark/executions/mpra-enriched-baseline-20260905T192000Z-21454084/ladder/ladder_gains.tsv`.
HyenaDNA scores 0.216. The baselines score: allele-identity ridge (B0) 0.134,
CpG ridge (B1) 0.178, composition (B2) 0.176, and the prespecified flank and
dinucleotide ridge (B3) 0.167. The HyenaDNA gain over B3 is 0.049 (95% 0.0006
to 0.099), and the Bonferroni-7 lower bound is -0.015. The gain over B1 is 0.038
(95% -0.006 to 0.083). gkm-SVM scored 0.126. MPRALegNet was also scored (0.204
and 0.200 on a different endpoint, so the values are not comparable).
Correction: print the baseline ladder in 7B. State that the corrected gain
holds only against allele identity.

**H2. The prespecified continuum test missed and is not reported (line 214,
Supplementary Note 12).**
Source: `RNA-seq/results/histology_anchored_continuum/prespec/mnc-fib-v1-20260828T133158Z/MNC_FIB_v1_prespecification.json`.
The primary test was fixed before any fit: fixed projection, within disease,
against matched random gene sets, in GSE162694 and GSE213621, with Holm over the
two cohorts and both required to pass. Results:
- GSE162694: rho 0.547 against a null 95th percentile of 0.528 (109 samples).
- GSE213621: rho 0.465 against 0.457 (293 samples).
- Holm p 0.082 in each, so the test missed.

Cohort PC1 exceeded its null (p 0.001), but the prespecification bars promoting
PC1. The paper reports only the transfer misses. A related count is also
unreported. Matched random gene sets give a mean of 86 "continuum-associated"
programs against the reported 92 (p 0.109; source
`RNA-seq/results/histology_anchored_continuum/translation/matched-null-20260907T113510Z/`).
Correction: report both results, and drop "signature" from the continuum's name
as the prespecification requires.

**H3. Both focal programs, and the results built on them, depend on GSE244832.**
I re-derived this from
`figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv`.

Leaving GSE244832 out, the top-50-gene Jaccard is 0.026 for Stromal ECM and
0.000 for Ductular injury. For other left-out datasets it is 0.31 to 1.0. Across
all programs, 56 of 117 fall below 0.2 without GSE244832 and 40 fall below 0.1.
The stability rule at line 178 takes the median over left-out datasets, so a
gene set defined by one dataset still passes. GSE244832 has 430,704 hepatocyte
nuclei against 84,310 in GSE202379 (RESULTS.md line 2849).

Consequences:
- **Fibroblast Stromal ECM, line 184.** The pooled beta of 0.329 combines two
  halves (source:
  `fig34_scientific_corrections_20260930/content_lineage_cohort_effects.tsv`,
  50-cell threshold):
  - GSE244832: 13 donors, beta 0.369, HC3 p 0.026.
  - GSE202379: 13 donors, beta 0.245, HC3 p 0.26.

  Only the GSE202379 half tests a gene set fixed before those donors were seen,
  and that half alone is not significant.
- **Cell-cell communication, line 186 and Figure S4E.** Leaving out GSE244832
  (40 donors), CDH1 to PTPRM falls from 0.496 to 0.0034, and ST6GAL1 to EGFR
  from -0.299 to -0.070 (source:
  `current_candidate/source_qualified_ccc_20261014/fig4d_program_ccc_lodo.tsv`).
  "Both keep their direction" is true only literally. Both should be labeled as
  resting on one dataset.
- **Stromal ECM discovery and doublets.** Stromal ECM was discovered in
  hepatocyte nuclei with no doublet removal (Methods line 1585). 52 to 71% of
  hepatocyte nuclei that carry stromal markers score at synthetic-doublet level
  (RESULTS.md lines 2848 to 2850).

### Medium

**M1. Gene membership of the 436 (line 58 onward).** Supplementary Note 1 gives
the robustness numbers, but the main text does not. Add one sentence: 314 remain
at p12 1e-6, and 49 to 63% replicate between MVP and UKBB European GWAS of the
same enzyme.

MVP European rows also have an excess of PP.H4 >= 0.999: 48/242 against 19/387
for UKBB (Fisher OR 4.79; rows are not independent). MVP was fine-mapped with an
out-of-sample LD reference. Nine MVP-only near-1 genes were never LD-checked:
FKRP, IFI27L2, LEPR, LEPROT, MTMR3 and four others. IFI27L2 is a two-axis gene.

Add a per-gene replication state to the Catalog: replicated, tested and not
replicated, or no same-trait test.

**M2. Coding-variant anchors at neighboring genes (lines 66, 485 to 488).**
Several colocalized genes are anchored at another gene's coding variant:
- SAMM50 at PNPLA3 rs738409.
- MAU2 at TM6SF2 rs58542926.
- SPATA31H1 and SUPT7L at GCKR rs1260326.
- IFI27L2 and SERPINA11 at SERPINA1 rs28929474.
- CNOT1 at GOT2 rs11076256.
- SPP1 is in LD (r 0.896) with HSD17B13 rs72613567.

Of the 22 direct-trait genes, 3 are anchored at the PNPLA3, TM6SF2 or GCKR
coding variants. PNPLA6 and MAU2 fail the LD check. The Catalog needs a flag for
anchors at, or in LD with, another gene's coding credible-set member.

**M3. Bulk dependence on GSE213621 (line 112 onward).**
- GSE213621 carries a median 54% of per-gene information, which matches its
  effective n of 53.8.
- 409 of the 1,854 calls are called in GSE213621 alone among the four cohorts
  other than GSE126848.
- In that cohort, T2D is 47.8% in NAFLD against 17.4% in controls (Chen et al.
  2023, PMID 36535507). T2D is not modeled.
- Its controls are histology-defined (biopsy showing no chronic liver disease),
  so the earlier concern about them does not hold.
- The manuscript leaves "[ref: source studies]" where Chen 2023 belongs.

Heterogeneity across the five cohorts:
- Cochran Q BH q < 0.05 for 902 of the 1,854 calls (671 without GSE126848).
- A DerSimonian-Laird random-effects z at BH 0.05 keeps 1,403 of the 1,854.
- With k = 5, truncated Hartung-Knapp keeps 0, so it cannot serve as a filter.

State that the calls estimate an average effect, about half of it from one
cohort. Add tau and the Q q-value per gene to the Catalog.

**M4. Two-axis gene directions (Fig 6C, line 1541).** The disease-liability
sidecar adopted on 2026-10-25 is absent from r9; "liability" has 0 hits in the
manuscript. HKDC1 reads as "opposite" only because it is oriented on the
AST-raising allele. That allele lowers MVP NAFLD (z -2.57) and cirrhosis
(z -4.69) risk, so on the disease-risk allele HKDC1 is concordant. Orient each
gene on the disease-risk allele where the sidecar has one, and label which
allele was used.

**M5. Single-cell reporting.**
- **Figure 4E colors.** Stromal ECM and Ductular injury are colored
  "Hepatocyte". Color by content lineage.
- **Line 178.** Add the minimum detectable effect and the post hoc joint result
  (18 programs against 5.85 expected, labeled exploratory).
- **Line 180.** Add the decontX ambient burden: 0.81 for Stromal ECM against
  0.34 for Ductular injury. Across 106 programs, burden predicts correction size
  (Spearman 0.48).
- **Line 182.** Report the prespecified V3 random-set control
  (`Analysis/Hypotheses/V3-20261008T015239Z/results/v3_random_set_null.tsv.gz`).
  87% of 1,000 expression-matched hepatocyte gene sets rise with stage. 8.4% of
  them exceed the Ductular injury stage t, and 4.1% exceed its disease t.

**M6. Continuum reporting (lines 142, 202, 206, Methods 1443 and 1473).**
- **Composition share.** Add that composition explains 80% and 82% of
  within-stage continuum variance.
- **Observed-weight denominator.** Figure 4F and S4V use retained weight as the
  denominator. Stromal ECM has 80.5% retained but 73.8% observed. Nine programs
  fall below 80% observed weight; mark them descriptive.
- **Unit count.** Methods line 1473 says 435 units, but the 0.760 and 0.707
  coefficients come from 367.
- **Heading, line 206.** Retitle "Forty-three gene-membership systems share one
  diffuse continuum association".

**M7. Figure 5D.** All ten expectation diamonds are 0% by construction: a
posterior spread evenly over a window can never put 50% in peaks.

Replace the "2.5 to 3.1 times" reading with the matched result. Enrichment sits
in peaks shared across lineages. Hepatocyte-specific peaks give 1.77. Other
lineage-specific peaks are below expectation. A platelet trait shows the same
shared-peak enrichment. Consider moving 5D to the supplement.

**M8. Figure 7F and 7A.**
- **7F.** TranscriptFormer minus scVI (linear head) is +0.0154 donor-pooled on
  the panel and -0.0042 as the equal-study mean in RESULTS.md. Both intervals
  exclude zero. Name the estimand drawn, print both, and fix the 7F citation at
  line 553.
- **HyenaDNA caQTL head.** The head on Currin caQTL leads (gain 0.130) is
  unadopted and undisclosed. Methods line 1737 should mention it.
- **7A.** The bar sums to 126 entries, including 37 no-finding placeholders.
  Legend (A) says 104.

**M9. Prior work on GNMT, MAT1A and CYP2C19 (line 144).**
- MASLD-specific reports of lower GNMT and MAT1A: PMIDs 23916847, 26394163,
  31668391, 19651758 and 31199045 (Männistö 2019, bulk NASH).
- PMID 11131452 reports the decrease in cirrhosis and HCC, not MASLD.

Describe the trio as reproducing reported decreases. Also, the protein
covariation is not specific to the trio: random triplets of co-lowered proteins
reach the trio's mean rho 5.1% of the time.

**M10. Title and Catalog labels.**
- **Title.** It says "genetic evidence", but the content is liver eQTL
  colocalization. PNPLA3, TM6SF2 and MTARC1 sit in the unresolved class.
- **"Testable" label.** The Catalog labels 935 disease-state genes "testable",
  but 861 of them have no liver-eQTL credible set.
- **Figure 6B.** It prints 11,301, which is correct for Catalog genes. The count
  over all 24,115 fit genes is 11,335.
- **Missing genes.** 114 Ensembl IDs have no Catalog record because of
  symbol-collision quarantine. For example, both CAST IDs are missing, one with
  q 1.0e-9.

### Low

- **GNMT locus.** Add the enzyme-GWAS PP.H3 sentence after line 82. Move the
  allele-direction paragraph to Supplementary Note 2.
- **Figure 5F.** Rename the key entry "supported" to "exceeds matched-gene
  null". The "source-dependent" legend sentence applies to Stromal ECM only.
- **Summary, line 20.** Replace "proposes a next experiment for each gene" with
  class-level wording.
- **Figure 3C legend.** Say that GSE135251 accounts for 325 of 346 mixed-sign
  genes.
- **lncRNA direction.** It changed with the control roster: 56.3% of 533 calls
  were up on the 844-record basis against 47.9% of 549 now. The prespecified
  per-cohort fits are not reported [Not re-derived].
- **GSTM1.** The canonical GSTM1 call may reflect the frequency of its common
  whole-gene deletion, which differs by ancestry. It drops to padj 0.10 with
  ancestry PCs [Not re-derived].
- **Figure 2D.** Enzyme shares may be inflated by non-European loci that reach
  the L = 10 credible-set limit. 163 EAS and 43 AFR enzyme loci hit it
  [Not re-derived].
- **Supplementary Note 3.** Add the MAGMA AST-only result from CGB4 (UKBB Holm
  0.015; ALT null).

### Checked and not upheld

- **GSE213621 controls lack histology.** Refuted; Chen 2023 defines them by
  histology.
- **Flex per-cell effects are explained by UMI depth.** Refuted. The effects
  survive a log-UMI covariate, exceed a depth-only ceiling, and SOX9/EPCAM rise
  in GSE256398, where MASLD nuclei are shallower. Binomial thinning was not run.
- **The overlap excess is a variance artifact.** The SE-decile OR of 1.03 is
  already in Supplementary Note 3.
- **The 996 union is wrong.** It is correct; the 997th record in one count is an
  anonymous ABF-only record.
- **The direct versus enzyme coding contrast is a main-text claim.** The
  manuscript already withdraws it after H58.

## 4. Results the paper understates

1. **Bulk fibrosis effects of programs follow gene-content lineage.** By content
   lineage:
   - Cholangiocyte content: 20 of 21 programs positive.
   - Fibroblast content: 19 of 21 positive.
   - Hepatocyte content: 29 of 36 negative.

   The Kruskal-Wallis p is 3.7e-11 by content lineage against 0.010 by discovery
   lineage (descriptive). In plain terms, bulk program shifts with fibrosis are
   mostly shifts in cell mix: more bile-duct and fibroblast RNA and less
   hepatocyte RNA. This is a cleaner statement than "0 of 224 calls".
2. **Enzyme dominance is a locus-count effect.** Per genome-wide lead, direct
   and enzyme traits yield colocalizations at the same rate. This supports the
   paper's own reading at line 74 with a direct number.
3. **The single-cell stage null is a power limit.** The minimum detectable effect
   and the post hoc joint test show the data cannot rule out modest,
   mixed-direction stage shifts. Saying so is more accurate than "no program
   shows".
4. **Colocalized genes are more loss-of-function constrained.** CGB8:
   OR 0.70 per SD LOEUF, Holm 0.012, prespecified, robust to eQTL count. It
   partly explains the small overlap. Cite Mostafavi 2023, which reports the
   general pattern.
5. **GTEx fat direction.** 323 of 326 genes keep the pooled sign in 222
   independent donors. This is the cleanest independent support for the bulk
   calls the paper has.

## 5. Next analyses, ranked

Each analysis needs a written prespecification before outcomes are read. The
unit is the donor or participant throughout.

| Rank | Analysis | Data and n | Family | Cost [Estimate] | Verdict | What it decides |
|---|---|---|---|---|---|---|
| 1 | Focal-program dependence on GSE244832 | Atlas hepatocytes. Hotspot rerun without GSE244832, or with it downsampled to GSE202379 nuclei. Score best-matching programs in its 18 donors. Report CCC and fibroblast contrasts per cohort with CIs. | Descriptive stability diagnostic. It does not replace the frozen 117. | One cpu job, hours | New, from section 3 H3 | Whether Stromal ECM and Ductular injury are atlas programs or GSE244832 programs. They run through Figure 4, the continuum section and Figure 5F. |
| 2 | Per-hepatocyte GNMT, MAT1A, CYP2C19 in GSE202379 and GSE244832 | BF3 pseudobulks (`build-20261023T155101Z/donor_type_counts.npz`). Ordinal stage in 37 and 18 non-cirrhotic donors; log median UMI covariate. Negative controls: identity genes and matched random genes. | Holm over 3 genes x 2 cohorts | Minutes | Do now | Takes the only candidate new biological result from two cohorts (one confounded with plate) to four. Bulk GNMT loss is published, so the per-cell result is the claim. |
| 3 | Continuum random-axis null plus program-level portal/residual split (continuum-N2, N1) | 774-unit expression; K = 1,000 matched sets in `mnc-fib-v1` (check the matching on 774 units first); BF4 marker scores with each program's own genes removed | Holm over 2 programs x 2 axis rules; BH within 43 systems | Under 1 h, plus a few hours of bootstrap | Do now | Gives 0.760 and 0.707 a reference distribution, and settles the Stromal ECM wording. |
| 4 | External bulk replication in GSE281797 (bulk-N1) | 91 participants (30 controls, 6 of them F1; 27 female and 3 male controls) | One slope and sign-agreement statistic over the 1,854 calls; one-sided BH within the 1,854 | Small cpu job, 1 to 3 h with permutations | Do now | The only open RNA-seq cohort with histology-normal controls that is outside the pooled fit. Zhao 2026 (PMID 41870035) owns the cohort. |
| 5 | Obese against lean controls (bulk-N2) | GSE48452 (HuGene 1.1 ST: RMA plus limma, not voom; 16 obese vs 12 lean). GSE126848 (13 lean, 10 obese; overlaps the discovery fit). | 3 slopes per cohort, Holm | About 1 h | Do now | How much of the pooled log2FC is already present in obese normal liver. |
| 6 | Gene-membership robustness pass (genetics-N1 plus N2) | 767 positive rows; the r2 LD panel (the PolyFun UKBB path was not found [Unverified]) | Bonferroni over 2 comparisons, clustered by locus | 13 to 26 h on 4 CPUs | Do now | Which of the 436 genes stay. Adds the per-gene replication state. |
| 7 | Coding-neighbor flag, then coding-conditional COLOC where LD is moderate (interface_catalog-N2) | r2 signal pairs; 1000G EUR | Locus-level | Minutes for the flag, hours for conditioning | Do now | How about 31 Catalog genes are read. Report high-LD loci as "not separable". |
| 8 | Heterogeneity core (bulk-N3) | Per-cohort fits, already re-derived | DerSimonian-Laird plus prediction interval, or a 4-of-5 sign rule | Under 1 h | Do now | Per-gene tau and Q in the Catalog. |
| 9 | Minimum detectable effect and exploratory global test, then a pooled donor-level retest of the 117 programs (singlecell_programs-N6, N1) | About 131 hepatocyte donors across 5 cohorts | BH over programs testable in a prespecified minimum number of cohorts; global sum t-squared prespecified | Minutes, then 10 to 20 CPU h | N6 now; N1 after prerequisites | The only powered single-cell stage test. Blocked by the GSE289173 holdout decision (yours) and by the coverage gate dropping each focal program in one cohort. |
| 10 | SteatoSITE F-SS6 as prespecified (framing-N6) | 625 staged participants; 487 of 657 libraries processed | As written in `PRESPEC_SS.md` | Light once genotypes exist | After processing finishes | Best-powered test of the inflamed-hepatocyte state and its link to PNPLA3. |
| 11 | Allele-direction coherence at three-way loci (models_fig7-N2) | H53 `three_way.tsv`, Currin caQTL | Holm over 2 tests, about 80 loci | Up to 1 h | Do now | A risk-allele direction for each genetically anchored candidate. The discordant loci carry the information. |
| 12 | Binomial thinning in Flex, then the CCL20 dissociation control (exploratory-N1, N3) | GSE298719 (13 donors), then 4 cohorts (94 donors) | Holm within cohort | Hours | N1 now; N3 after N1 | Needed before any per-cell claim from GSE298719. |
| 13 | GNMT conditional analysis (genetics-N3) | Ghodsian NAFLD, UKBB and MVP enzymes | Single-locus sensitivity analysis | Under 1 h | Do now | Whether the NAFLD signal is separate from the enzyme lead. Report it whatever it shows. |

Not recommended:
- **continuum-N3, paired intervention biopsies.** It was already null in the
  same pairs, there are only 6 resolvers, and GSE106737 shares arrays with
  GSE83452.
- **framing-N1, enzyme-locus disease-risk ratio.** H10 already answers it, is
  prespecified and replicated in Ghouse. Report H10 instead.
- **framing-N3, L2G benchmark.** Recall at coding loci is low by design, L2G
  shares eQTL inputs with the Catalog, and a gene ranking conflicts with
  PAPER.md scope.
- **interface_catalog-N3.** It repeats CGB4 and DIG1 section 43.
- **context-N1, HMSMA deconvolution.** It is the most expensive proposal, the
  HRA007511 GSA-Human access status is unresolved, and Li 2025
  (PMID 41286103) points to the expected result.

Data-access flags:
- The MVP summary statistics come from dbGaP phs002453. Methods line 1487 calls
  this open-access dbGaP; confirm that for the exact files used.
- HRA007511 is a GSA-Human deposit. It was downloaded by public FTP, but its
  access tier has not been confirmed.

## 6. Recommendation

Fix H1 and H2 in the text now; both are wording corrections with numbers already
on disk. Then run rank 1 (the GSE244832 diagnostic) first, because it decides
whether Figure 4's two focal programs and the results built on them stand. Run
rank 2 (four-cohort per-hepatocyte GNMT) the same day; it takes minutes on
existing pseudobulks and is the most direct route to a positive biological
statement in the Summary.

## Sources for this review

- Workflow outputs (session scratch, not durable):
  `~/.claude/jobs/c3972821/tmp/w1/w1.json` (area reviews) and
  `~/.claude/jobs/c3972821/tmp/w2/w2.json` (verification, feasibility and the
  final completeness pass).
- Files re-checked directly: the paths given under H1, H2 and H3, and
  `docs/RESULTS.md` lines 2846 to 2852.

## Addendum, 2026-11-02: ranks 1 and 2 run (RV1)

`Analysis/Hypotheses/RV1-20261031T155620Z/RESULT_RV1.md` (job 22603883).

- Rank 2 does not support per-hepatocyte GNMT loss as a candidate new finding.
  Across diagnosis category in GSE202379 (36 donors) and GSE244832 (18 donors),
  no test of GNMT, MAT1A or CYP2C19 passes Holm over 6. GNMT has opposite signs
  in the two cohorts (+0.165 and -0.963). Sections 1 and 6 above no longer hold
  on this point. A loss limited to advanced fibrosis, the DIG1 contrast, is
  not excluded.
- Rank 1: Hotspot without GSE244832 recovers neither focal program (best top-50
  overlap 2 genes and 0 genes). Their member genes still co-vary across 45
  GSE202379 donors (Holm p 0.002 each), but 27 of the other 28 hepatocyte
  programs do too. So the gene sets are not specific to GSE244832, and their
  cell of origin stays unresolved.

## Addendum, 2026-11-02: ranks 3, 5, 7, 8, 11 and 13 run (RV3 to RV8)

All on nslab, prespecified, exploratory; none is adopted. Results files are in
`Analysis/Hypotheses/RV<n>-*/RESULT_RV<n>.md`.

- **Rank 3, continuum random-axis null (RV3).** Random expression-matched axes
  give Stromal ECM and Ductular injury coefficients like the observed 0.760 and
  0.707 (permutation p 0.11 to 0.21, Holm 0.43 to 0.46). Random axes call a
  median of 93 of 117 programs, against 94 observed. The continuum coefficients
  are not specific to the continuum axis. The portal/residual split (N1) was not
  run.
- **Rank 5, obese against lean (RV4, GSE48452).** Over 1,036 matched calls the
  obese-versus-lean slope is -0.006 (Holm p 0.56) and the NASH-versus-lean slope
  is 0.336 (Holm p 0.003). Obese fraction -0.02 (bootstrap 95% interval -0.61 to
  0.31). The calls are not detectable in obese normal liver; up to about a third
  of the NASH-scale effect is not excluded.
- **Rank 8, heterogeneity core (RV5).** Of the 1,854 calls, 746 are consistent
  across the five cohorts, 651 depend on one cohort (GSE213621 for 625 of them),
  358 agree in direction but not size, and 99 are opposite in two or more
  cohorts.
- **Rank 11, allele directions at three-way loci (RV6).** At 67 of 88 loci
  (0.76; Holm p 1.8e-6) the allele that raises liver accessibility also raises
  the colocalized gene's expression. Genetic and disease directions agree at 28
  of 49 loci (Holm p 0.39), not shown, as in DIG1. 48 of the 49 are liver-enzyme
  loci.
- **Rank 13, GNMT conditional (RV7).** Conditioning on the enzyme leads leaves
  the NAFLD signal unchanged (anchor z 3.65 to 3.56 or more; PP.H4 0.64 to 0.69).
  But the leads are 0.8 Mb away (r^2 < 0.005), so this test says little.
  Post hoc: in four of five enzyme GWAS the near-1 PP.H3 comes from that distant
  signal, and within 250 kb of GNMT these GWAS show no association strong enough
  to call H3. UKBB ALT has a local signal (p 1.5e-13) that does not colocalize
  with the GNMT eQTL (PP.H3 0.88); MVP ALT shows none at the anchor (p 0.70).
  The sentence proposed under "GNMT locus" in section 3 ("all five ... favours a
  different variant") is wrong as written; use the RV7 wording.
- **Rank 7, coding-neighbor flag (RV8).** At the COLOC anchor, 4 of 436 genes
  sit on or in r^2 >= 0.8 with another gene's coding credible-set variant
  (C2orf16 and CWF19L1 among the 22 direct-trait genes; RECQL4 and ZFTRAF1 at
  LRRC24), 5 are in moderate LD, and 29 could not be assessed with the 1000G
  panels. The M2 list in section 3 mostly does not reproduce: SAMM50's anchor
  has r^2 0.05 with PNPLA3 rs738409 (EUR), and the SERPINA1, TM6SF2 and GOT2
  examples are not in LD with their anchors. Replace M2 with the RV8 table.
- **Rank 6, part 2, per-gene replication (RV9).** Using independent European
  partners (MVP and UKBB for ALT and AST; MVP and FinnGen for NAFLD), 78 of the
  211 genes with a testable comparison replicate, 121 have a distinct-variant
  call, 12 neither; 225 genes have no testable comparison. Replication rises
  with discovery PP.H4: 0.27 (0.16 to 0.39) at 0.5 to 0.8, 0.52 at 0.8 to 0.95,
  0.67 (0.57 to 0.75) at 0.95 or more. Distinct-variant calls may reflect the
  different LD references; part 1 (N1) is needed to separate that.
- **Rank 4, external replication in GSE281797 (RV2).** 1,063 of 1,854 calls
  testable (91 participants, 29 controls). Sign agreement 0.653 (permutation
  Holm p 0.13; null 95th percentile 0.70), so direction replication is not shown
  under the prespecified rule. The slope of replication on discovery log2FC is
  0.21 (Holm p 0.016; down calls 0.25, up 0.16). No call replicates at BH 0.05.
  Post hoc, calls RV5 classes as consistent agree in sign at 0.72, against 0.50
  for cohort-discordant calls.
- Still open: rank 6 part 1
  (N1 LD-consistency refit of all positive rows, 13 to 26 h; needs the r2
  runner rebuilt with the AN5 diagnostic), continuum-N1, ranks 9, 10 and 12.

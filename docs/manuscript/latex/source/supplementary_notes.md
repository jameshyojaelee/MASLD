# Supplementary Notes

## Supplementary Note 1. Scope of the genetic and disease-state evidence layers

Liver-trait genetic evidence in this Resource takes three forms: multi-signal liver eQTL colocalization, single-signal liver eQTL colocalization, and fine-mapped protein-altering or canonical-splice credible-set members. The Gene Catalog states which form applies to each gene. None of them alone identifies the gene that mediates a trait association. Of the 436 primary multi-signal colocalized genes, 423 take their strongest colocalization evidence from a liver-enzyme trait and 13 from a direct disease or liver-fat trait. Membership of genes near the threshold depends on the colocalization prior p12. At p12 = 1 × 10⁻⁶ instead of the 5 × 10⁻⁶ used for all calls, 122 of the 436 genes, those with best PP.H4 between 0.50 and 0.83, fall to PP.H4 ≤ 0.5, leaving 314: 111 of the 400 genetically anchored genes, 8 of the 26 supported on both axes and 3 untested genes whose RNA was not detected. MVP and UK Biobank European GWAS of the same liver enzyme share 49 to 63% of each other's colocalized gene calls. This is gene-level recurrence with a reused eQTL panel, not matched-mechanism replication. In the historical named-r8 prior sensitivity, at p12 = 1 × 10⁻⁵, 29 unresolved genes tested at PP.H4 ≤ 0.5 rise above 0.5. That sensitivity was not rerun under the adopted eligibility qualification. The adopted calls retain p12 = 5 × 10⁻⁶; prior sensitivities do not change classes.

Fourteen colocalized genes rest only on tier-1/2 GWAS fits with SuSiE PP.H4 of 1.0 (at least 0.999999; 18 gene-by-GWAS tests). We refit these tests with the inputs and seed of the colocalization run, which reproduced every PP.H4. All 18 fits converged, and PP.H4 stayed at 1.0 with 1 × 10⁻⁶ instead of 1 × 10⁻³ added to the LD diagonal and with the residual variance estimated. In 14 of the 18 tests, the GWAS variant that drives the colocalized signal does not fit the LD panel: its z score differs from the value predicted from nearby variants through LD by more than 5.45 standardized units (the χ² 1-df cut at p = 5 × 10⁻⁸), against at most 0.92 in three comparison tests (RORA with MVP European NAFLD, FABP1 and MAU2 with UK Biobank ALT). After variants above this cut were removed, 11 of the 14 genes no longer had PP.H4 > 0.5 (BRCA1 kept no GWAS credible set), and GALK1 (0.550), PKN3 and TMPRSS2 kept it. All 11 are genetically anchored genes supported by MVP European GWAS. We chose this filter after the first diagnostic, so it is a post hoc sensitivity analysis. It may also remove true signals at loci with several strong associations, which we did not test, and it changes no prespecified class.

Disease-state liver expression has two parts. The first is the pooled, cohort-adjusted disease-versus-control calls from 774 genetic-identity units in five cohorts (138 control, 636 disease). The second is the cross-sectional fibrosis-stage contrasts. The disease-versus-control fit tested 24,115 genes. The units are genetic clusters, not verified donors. Stage contrasts compare different people, not change within a person.

We also describe how fine-mapped posterior probability at liver-trait loci divides among protein-altering, splice and noncoding variants. In that analysis, the contrast between direct traits and liver enzymes is descriptive. It rests on 16 direct-trait regions at genome-wide-significant leads and is pooled across ancestries. It is confounded with GWAS ancestry and fine-mapping coverage, and no ancestry stratum supports a trait-class call. It does not resolve the responsible variant at any locus.


**Variant-identity and LD qualification.** A complete saved-eQTL-fit check found
duplicate identifiers in 3,633/18,889 fits, affecting 90/436 primary positive
genes and 137/708 full-portfolio positives. At MPRIP, a repeat insertion and
deletion were assigned one position ID. MVP allele labels and frequency disagree
about the event, and the original rounded and unrounded reference LD remain
indefinite after regularization. Three event interpretations under each matrix
therefore give no qualified H4. A genotype-derived sensitivity on 388 EUR
pedigree founders yields valid LD and converged GWAS fits, but neither eQTL fit
converges within 500 iterations, including repeat exclusion. MPRIP remains
untestable. DIP2B also has a duplicate position in a supported eQTL credible set;
all affected genes require full refits, not post-fit column deletion. These
findings qualify the historical genetic counts and their downstream claims;
they do not establish that 90 genes are false positives.

## Supplementary Note 2. GNMT colocalization and linkage, and stage association of the two focal programs

GNMT has a single-signal colocalization with direct NAFLD (ABF PP.H4 0.706) but no multi-signal SuSiE call (best PP.H4 0.128). The disease side is sub-threshold (anchor GWAS p = 2.6 × 10⁻⁴; smallest p within 1 Mb 1.4 × 10⁻⁴), and PP.H4 depends on the prior p12: 0.33, 0.71 and 0.83 at p12 = 1 × 10⁻⁶, 5 × 10⁻⁶ and 1 × 10⁻⁵ in the earlier run of this pair. Its Catalog genetics call is therefore indeterminate, and its local regulatory haplotype is unresolved.

The rs2296804 anchor has SNP.PP.H4 0.259. rs2296805 has UK Biobank European r2 = 0.998 with rs2296804, and rs9471976 and rs11752813 have r2 = 0.818 and 0.580 with rs2296805.

For Stromal ECM (IGFBP7), the bulk fibrosis association does not reach BH q < 0.05 in the 364 units from the three cohorts of the 774-unit set that record both fibrosis and NAS (GSE130970, GSE135251 and GSE162694) (q = 0.214; q = 0.421 after composition adjustment). It reaches q = 0.045 only when 93 GSE174478 source records are added and the Knapp-Hartung variance is not floored at the REML value. In that analysis, q is 0.230 after composition adjustment. In the single-cell stage models, once dataset-specific variances are modeled, neither focal program passes within-dataset permutation BH over all 117 programs for a cross-sectional stage association (smallest q of any program 0.18135 in the October 26 source-diagnosis correction).

## Supplementary Note 3. Gene-level comparison of colocalized and disease-state genes

### Overlap of the two gene sets

After the adopted shared-posterior qualification of carried negatives, multi-signal colocalization evaluates 3,589 genes of the 14,704-gene native RNA–genetics interface. These include 426 colocalized genes, 209 canonical disease-state calls and 26 shared genes, against 21.3962 expected within expression deciles (Mantel-Haenszel OR 1.28151, 95% CI 0.82618 to 1.98776; exact stratified *p* = 0.276987, the original test definition with updated eligibility). The interval includes 1; it does not establish equivalence or biological separation. Historical expression diagnostics were not repeated in this qualified universe, so no direction is inferred. Before qualification, the prespecified 3,769-gene universe contained 426 colocalized and 220 disease-state genes, with 26 shared against 21.2127 expected (OR 1.29475, CI 0.83509 to 2.00740; exact *p* = 0.275500). Its alternative expression adjustments, OR 0.958 to 1.347, the matched expression control (13 against 21.3 expected), and the expression curves predicting 20.2 to 20.8 shared calls remain historical diagnostics. They do not calibrate the updated universe.

Among the 14,704 native-interface genes, 10,928 have no multi-signal result and 187 are priority-untestable; these groups contain 935 and 11 canonical disease-state calls, respectively. A missing multi-signal result retains an indeterminate Catalog genetic call; a failed applicability check is untestable. Neither is a genetic negative. The historical full-interface comparison, which treated unevaluated genes as non-colocalized, gave 26 shared against 22.3 expected (OR 1.19, CI 0.79 to 1.81); the single-signal ABF comparison gave 49 against 43.6 (OR 1.15, CI 0.85 to 1.57). These are different estimands from the qualified multi-signal comparison. The annotation-adjudicated Catalog has 14,502 jointly eligible identities and does not replace the native stable-ID universe. Historical unstratified power calculations and expression controls in the 3,769-gene comparison are retained for reproduction, not claimed as updated checks.

The original colocalization rerun checked first-run positives. Applying the same check to 12,087 carried primary negatives makes 188 more genes untestable under the existing priority rule; 180 were in the earlier 3,769-gene comparison, leaving 3,589 evaluated genes. The November 6 adoption applies this qualification to both the native comparison and the Gene Catalog. It changes 188 genetic states and 180 classes: 169 unresolved and 11 established-state-associated genes become untested, while eight genes were already untested. Positive genetic genes and all RNA calls, effects and directions are unchanged. Under the alternative rule that any eligible primary study makes the gene testable, 3,736 genes remain (426 colocalized, 217 disease-state calls, 26 shared; 21.1460 expected, OR 1.30032, CI 0.83838 to 2.01678, exact *p* = 0.274636). This is a sensitivity rather than the adopted aggregation rule. Historical Figure 3G versions with 3,776 or 3,769 genes are superseded by the qualified 3,589-gene panel.

### Expression-matched comparison of hepatocyte specificity

In this historical comparison using named r8 classes, which was not rerun after November 6 qualification, colocalized genes are not more hepatocyte-specific than disease-state genes, so the prespecified directional prediction failed. The comparison contains 281 expression-matched gene pairs (expression SMD 0.00371). In a post hoc comparison, the mean genetic-minus-state difference went the other way (−0.265515 log2 units; 126 positive and 155 negative pairs). The relabel-and-rematch outputs (*p* = 0.00079992; BH *q* = 0.00159984 across two outcomes) are diagnostic only. That null destroys the observed expression imbalance and matching structure, so it is anti-conservative, and family correction cannot repair it.

The separate post hoc paired sign-flip diagnostic gives nominal *p* = 0.01482985 (100,000 draws, seed 20260907). It assumes exchangeability within pairs. Flipping signs jointly within the 216 loci that hold the 281 colocalized genes gives nominal *p* = 0.01606 (100,000 draws, seed 20260907; 0.01625 with seed 20260923). A sign test gives *p* = 0.094677 (126 of 281 pairs positive), and 99 of 216 locus means are positive (*p* = 0.247). All of these tests are post hoc, and none is calibrated. The result is a descriptive difference conditional on selection. It is not confirmatory class separation or independent biological replication.

Of the 1,129 state-only genes, 1,124 have a single-cell hepatocyte-specificity value and form the comparison pool. Of these 1,124, 931 had no multi-signal colocalization result, so most of the pool consists of genes that colocalization could not test. The 193 tested genes have mean AveExpr 1.3405, against 0.4164 for the whole pool. In a second post hoc comparison restricted to those 193 genes, the difference kept its direction. In 129 expression-matched pairs (expression SMD −0.0004), the mean genetic-minus-state difference was −0.342 log2 units (median −0.183, paired d −0.217; nominal sign-flip *p* = 0.0154, 100,000 draws, seed 20260923; Wilcoxon *p* = 0.0314; 57 positive and 72 negative pairs, sign test *p* = 0.218). This design had 80% power for d = 0.247.

Within the 281 published pairs, the difference comes mainly from the 56 pairs whose state partner was tested (mean −0.748, nominal *p* = 0.00148), not from the 225 pairs with an untested partner (−0.145, *p* = 0.233). Genes and pairs are the units in these analyses. The sign-flip null does not account for genes that share loci, and none of these *p* values is calibrated. The 1,000-draw matched-pair label control is centered near its analytical expectation of zero; it does not validate the relabel-and-rematch null. These comparisons do not establish a formal biological separation between the two gene classes.

### Protein-response comparison

This historical named-r8 comparison, not rerun after eligibility qualification, asks whether the liver protein response to fibrosis stage differs between the two classes. It contains 73 pairs, with a mean residualized difference of −0.015601, a paired effect d of −0.08156 and a minimum detectable effect of 0.3279 at 80% power. We detected no difference (post hoc paired sign-flip nominal *p* = 0.487465). This is not evidence of equivalence. The relabel-and-rematch *p* and *q* (both 0.459654) are invalid diagnostics under the same null limitation. The unresidualized protein comparison is circular: class selection uses the RNA effect, and protein abundance tracks that effect. It therefore provides no independent biological evidence.

### Expression observability of the two classes

In the historical named-r8 classes, among genes with a single-cell hepatocyte-specificity value (399 of 400 genetic-only and 1,124 of 1,129 state-only genes), genetic-only genes have mean bulk AveExpr 4.2115, against 0.4164 for state-only genes, a difference of 3.7951 log2 units. Of the 399 genetic-only genes, 281 (70.4%) obtain an expression-matched partner.

### Gene-set association by MAGMA

In a MAGMA gene-set test of the colocalized genes, the association was specific to AST. In UK Biobank AST, beta_std was 0.0276 (*p* = 0.0073; Holm *p* = 0.015 over two traits), and MVP AST gave *p* = 0.047. ALT, liver fat (PDFF) and cirrhosis showed no association.

## Supplementary Note 4. Fine-mapped coding and splice posterior by trait class and ancestry

Averaged over physical regions (credible sets chained by single linkage when their lead positions lie within 1 Mb) and pooled across ancestries, credible sets for direct MASLD and liver-fat traits carry 29.65% of their posterior probability in protein-altering or canonical-splice variants, against 6.67% for the liver enzymes ALT, AST and GGT. The difference is 0.2298 (95% physical-cluster bootstrap interval 0.0808 to 0.3958). Only the 4,107 reliable credible sets whose lead variant has p < 5 × 10⁻⁸ enter these estimates.

The pooled contrast mixes trait class with GWAS ancestry and fine-mapping coverage. Within European studies, the fractions are 29.88% and 6.92%, and within admixed American (AMR) studies, 8.11% and 6.91%. No stratum supports a calibrated trait-class call. The European estimate failed the coverage check (UKBB enzyme credible sets lie at 8 of 767 lead windows by locus identifier and within 500 kb of 23) and the power check, and passed its random-split calibration control. The AMR and combined estimates failed their random-split calibration controls at 95% and at the Bonferroni-corrected 98.333% level, so their intervals are too narrow. Both also narrowly failed the power check, and the AMR direct arm, 4 clusters, was empty in 156 of 10,000 draws, above the prespecified limit of 100. The pooled direct-minus-enzyme difference comes mainly from the three PDFF GWAS and 5 of the 16 direct clusters. The ancestry-stratified estimates used 10,000 shared draws over 44 pooled clusters (seed 20260923). The 50 robustness comparisons gave ratios from 3.44 to 5.18.

In a sensitivity analysis, we used all 4,134 credible sets at genome-wide leads, including those that failed the PIP-sum criterion (seed 20260817). Matching on standardized lead-variant effect gave a direct-to-enzyme ratio of 3.80 (95% interval 1.29 to 39.62; 19 matched enzyme loci). Direct credible sets come from smaller GWAS (median effective n 21,893 against 160,000, over the 165 direct and 3,969 enzyme credible sets in this analysis) and carry larger standardized effects (median 0.0423 against 0.0190).

To group credible sets into regions, we placed each set at the lead position in its locus identifier and joined sets by single linkage when their lead positions lay within 1 Mb, so one region can span more than 1 Mb. Where the identifier had lost trailing zeros, we restored the one position within 1 Mb of all the set's member variants. For 715 of the 6,415 credible sets, the position in the locus identifier is too small by a factor of 10 or more and would place the set in the wrong region. All Figure 2D estimates use the restored positions.

Separately, the strongest evidence for 423 of the 436 colocalized genes comes from a liver-enzyme trait and for 13 from a direct disease or liver-fat trait. In addition, 22 genes have a direct-trait multi-signal SuSiE colocalization (PP.H4 above 0.5). Multi-signal colocalization tested far more enzyme than direct-trait gene-by-GWAS pairs (12,121 against 771).

## Supplementary Note 5. Signed agreement between ALT TWAS and the disease-versus-control contrast

To test whether genes with higher genetically predicted liver expression for ALT also differ in the same direction between disease and control liver, we compared ALT TWAS z-scores from OTTERS [ref: OTTERS TWAS method] with the disease-versus-control t statistic across 14,703 genes. The TWAS z-scores come from a single OTTERS analysis that was not recomputed with the other genetic analyses in this paper, so this comparison is secondary to the colocalization results.

Without composition adjustment, genome-wide rank agreement between TWAS z and the disease-versus-control t was bounded inside ±0.03 (expression-stratified Spearman rho 0.0107, 95% CI −0.0059 to 0.0272; within-chromosome shift *p* = 0.20). After adjustment for hepatocyte fraction, rho was 0.0167 (95% CI 0.0001 to 0.0332). This interval is not inside ±0.03, and we state no direction. Among the 1,527 genes with the strongest predicted-expression associations (|z| > 3), rho was 0.0421 (95% CI −0.0075 to 0.0916), and 0.514 of genes had the same sign, against 0.499 expected. The unadjusted bound excludes sign agreement of about 0.60 or more among those 1,527 genes. It does not exclude full agreement among the 103 of them that are canonical disease-state genes.

## Supplementary Note 6. Re-annotation of noncoding context by genomic position

Variants absent from the precomputed annotation table received no noncoding context call. Such variants were rare in European studies and common in non-European studies (0.7% of gated credible-set member rows without a call in EUR, 98.2% in AMR). We therefore re-annotated promoter and open-chromatin context by GRCh38 position, with all 16,324 sites lifting over 1:1, and ABC-enhancer context by GRCh37 position. After re-annotation, every noncoding posterior probability had a context call.

## Supplementary Note 7. Sensitivity analyses for the bulk disease-state layer

This note gives the secondary sensitivity statistics for Figure 3. Each section names the main-text result it qualifies. Units are genetic-identity clusters, not verified independent donors.

### GSE126848 per-cohort model

The GSE126848 model contains 23 control and 7 disease clusters, all male, so it omits the constant sex term. Native fragment counts and an independent refit with separate group means reproduce all 6,147 calls. A retained disease library is deeper than its paired control in 5 of 161 control-disease pairs (exact Wilcoxon *p* = 1.9e-5). Each of the 437 calls with zero reads in all 30 libraries has log2FC 0.554. This value equals the precision-weighted control-minus-disease difference in mean log2 normalized library size. Adding log2 library size as a covariate keeps 2,025 of the 6,147 calls (2,678 calls in total) when library size is counted before normalization, and 3,056 (3,649 in total) when it is counted after. Neither refit keeps any of the 437. A refit on the 15,592 genes that pass the cohort's own expression filter (filterByExpr within the 30 units, with RLE factors recomputed and BH over those genes) gives 3,440 calls (1,821 up, 1,619 down), and Figure 3C uses these calls. All 3,258 refit calls that are also among the 6,147 have the same sign in both fits, and the other 182 refit calls are not among the 6,147.

### Robustness of the 1,854 pooled calls

In the spike-in calibration, 99.8% of supported one-cohort effects (2.0 log2FC) lost their call in at least one leave-one-cohort-out refit. A uniform 1.0 log2FC effect kept its call in all five refits in 90.6% of cases, and losses concentrated in calls near the 0.5 log2FC threshold. The call count depends more on the fold-change floor than on the q threshold. In the Gene Catalog gene table, which holds 1,853 of the calls, a floor of 0.4 or 0.6 log2FC at BH q < 0.05 gives 2,826 or 1,229 calls, whereas a q threshold of 0.01 removes 22 of the 1,853.

With sample weights fixed at their primary values, adding within-cohort standardized log library size keeps 1,760 of the 1,854 calls when library size is the total count on the tested genes, and 1,648 when it is taken after RLE normalization. Adding top-100-gene share or expression entropy instead keeps 1,738 or 1,759. Adjusting for log library size, top-100 share and detected genes together keeps 1,563 (84.3%) or 1,569 (84.6%) calls under the two library-size definitions, and up calls fall from 1,333 to 1,085 or 1,100.

A strict-depth identity screen (read depth of at least 20 and at least 500 shared sites in each half of the sites) keeps 738 clusters. The primary model refitted on them tests 24,509 genes and gives 1,906 calls (1,344 higher and 562 lower in disease). Of the 1,854 primary calls, 1,779 are also called there, all with the same sign. GSE126848 keeps 5 controls and no disease unit under this screen, so the refit is a four-cohort disease contrast.

### Genetic ancestry sensitivity

Adding four ancestry PCs changed 233 calls, against a median of 120 (95th percentile 161) when the same PCs were permuted within cohort and group. Restricting to EUR-called units changed 570 calls (101 lost, 469 gained), against a median of 409.5 (95th percentile 501) for random subsets of the same cohort and group sizes. Both Holm-adjusted *p* values are 0.0199 over the two tests. One of the 200 permuted-PC fits and one of the 100 random subsets changed at least as many calls as observed, giving one-sided *p* = 2/201 (0.00995) and 2/101 (0.0198). The smallest attainable one-sided values are 1/201 and 1/101, so more reference fits would be needed to estimate these p values precisely. Across the 1,854 calls, the slope of the with-PC log2FC on the primary log2FC is 0.984 (descriptive). The ancestry-sensitive rule flagged 64 genes; the permutation reference gives a median of 6 expected false flags for this rule (95th percentile 17).

### Library properties of the stage contrasts

This section qualifies the fibrosis and NAS stage counts (Figures 3E and S3D). It is a post hoc descriptive deviation from the prespecified reading of the library-property check, whose positive control failed (Methods). With within-cohort log2 library size, Gini coefficient and normalization factor added together, the calls that stay BH-significant are 99 of 1,338 for F1, 2,113 of 4,266 for F2, 4,797 of 6,190 for F3 and 5,307 of 7,043 for F4 versus F0; 2,595 of 5,338, 7,296 of 9,337 and 7,083 of 9,841 for NAS 1–2, 3–4 and 5–8 versus NAS0; 3,006 of 4,305 for NAS 5–8 versus 3–4; and 356 of 398, 497 of 512 and 114 of 281 for F1 to F2, F2 to F3 and F3 to F4. These counts are not a test. The three properties are computed from the same counts and may reflect disease-related RNA composition, so a lower count does not show that a stage difference is technical.

### Transfer of grouped fibrosis-stage effects to GSE213621

The discovery-agreement ratio divides the GSE213621 slope on effects from two discovery cohorts by the held-out discovery cohort's slope on the same effects, with each slope averaged over the three held-out cohorts. It is 0.793 (95% interval 0.488 to 1.098). The F2 versus F0–F1 contrast was untestable, because the minimum detectable slope under the permutation null was 0.34, above the prespecified limit of 0.30. Adding GEO accession quartile as a covariate gives b = 0.719, so the F3–F4 agreement is not attributed to submission order. With estimated cell-type composition added to both fits, b is 0.431, above half of the primary slope. This adjusted result depends on the deconvolution input, and no gene is named from an adjusted fit. In women only (161 units), b is 0.742, against 0.708 in the primary fit.

### Differences at earlier fibrosis steps for GNMT, MAT1A and CYP2C19

GNMT and MAT1A have lower abundance at F1 than at F0 under all three F0 references (GNMT q = 0.00913 / 0.0426 / 0.0430; MAT1A q = 2.45e-4 / 9.5e-4 / 1.6e-3), but not without GSE130970 (q = 0.207 and 0.0892). CYP2C19 has lower abundance at F2 than at F1 (q = 0.0477). Its q is 0.052 with a control-status term, and 0.219 and 0.499 without GSE135251 and GSE162694. GNMT's F1-to-F2 q is 0.0562. We do not use these earlier-step differences to place where the lower abundance begins.

### Protein covariance of GNMT, MAT1A and CYP2C19

Across the 12 pair-by-adjustment tests in the same 46 participants, BH q is 0.00972 to 0.0232 without histology adjustment and 0.0232 to 0.0417 after adjustment for fibrosis stage. After adjustment for NAS, q is 0.136 to 0.214, and after adjustment for the four histology components it is 0.287 to 0.301. We report correlations and pointwise participant-bootstrap intervals for all adjustments, including a prespecified broader metabolic-abundance sensitivity (Figure S3 [panel letter to be assigned]).

## Supplementary Note 8. Stage and lineage sensitivity analyses for the single-cell programs

Units in the primary analyses are source-qualified donor identifiers. The 61-record equal-variance sensitivity includes two GSE189600 records whose biological independence is unresolved.

### Equal-variance stage sensitivity

The corrected equal-variance sensitivity retains at most 61 source records, including two GSE189600 records without ambient correction and with unresolved independence. Its complete 117-program ordinary and HC3 estimates are reported separately. Historical 62-record p/q values and cell-gate comparisons are superseded; they do not supply current primary evidence.

### Primary-model sensitivities

Under the corrected primary model, averaging stored scores over QC-passing cells selects 0 of 117 programs (smallest permutation BH q = 0.193); stage betas correlate with the primary estimates at Spearman 0.994. Ductular injury has beta = 0.208 (nominal 95% CI 0.085 to 0.331; permutation p = 0.00165, BH q = 0.193); Stromal ECM has beta = 0.149 (nominal 95% CI 0.004 to 0.294; permutation p = 0.04865, BH q = 0.271). These stored cell scores include neighbor contributions from failing cells, so this is not rescoring after QC. Removing donors with fewer than 20 T cells leaves 35 of 52 fitted T-cell donors, with 0 selected T-cell programs (smallest T-cell q after BH over the complete 117-program family = 0.552); the other 106 permutation p values remain fixed.

### Ambient-RNA audit

We next examined ambient RNA using unchanged scores, original full-transcriptome exposures and fixed program weights in the correction-available donor universe. Stromal ECM changes from 0.4207 to 0.1525 in the paired equal-variance analysis (delta -0.2682, nominal p = 0.00148, BH q = 0.173); its ambient-corrected GLS beta is 0.077 (permutation BH q = 0.331). Ductular injury changes from 0.4574 to 0.4504 in the paired equal-variance analysis (delta -0.0070, nominal p = 0.8588, BH q = 1.000); its ambient-corrected GLS beta is 0.172 (permutation BH q = 0.111). The GLS sensitivity tests 106 non-T-cell programs within the 117-program family and selects 0 (smallest q = 0.111). The 11 T-cell programs are untestable on the corrected native scale. Attenuation classification is conditional on primary selection, so it applies to neither focal program. These results do not establish the cellular origin of either program's transcripts.

### Content-lineage contrasts

We compared steatohepatitis with steatosis within each focal program's content-associated lineage, using the same correction-available cells and donors and the corrected RNA diagnoses. Ductular injury in cholangiocytes: beta 0.135 (HC3-based meta-analysis 95% CI -0.075 to 0.346; n = 26 donors; two-program BH q = 0.208). Stromal ECM in fibroblasts: beta 0.321 (HC3-based meta-analysis 95% CI 0.086 to 0.556; n = 26 donors; two-program BH q = 0.015). These are targeted projected-expression tests in GSE202379 and GSE244832, with a 50-cell primary gate; they do not independently replicate the hepatocyte discovery or exclude within-lineage substate mixtures. Cohort-specific estimates, full-raw, 20-cell and annotation-confidence sensitivities are reported separately. Captured-lineage fractions are not whole-tissue proportions.

### Dependence of the focal programs on GSE244832

When GSE244832 is left out, the top-50-gene Jaccard overlap with the full-atlas program is 0.026 for Stromal ECM and 0.000 for Ductular injury, against 0.31 to 1.0 when other datasets are left out. Across all programs, 56 of 117 fall below 0.2 without GSE244832 and 40 fall below 0.1. The stability rule takes the median over left-out datasets, so a gene set defined mainly by one dataset can still pass it. GSE244832 contributes 430,704 hepatocyte nuclei, against 84,310 in GSE202379. The fibroblast Stromal ECM contrast splits by cohort: GSE244832, 13 donors, beta 0.369 (HC3 *p* = 0.026); GSE202379, 13 donors, beta 0.245 (HC3 *p* = 0.26). Only the GSE202379 half tests a gene set fixed before those donors were seen, and that half alone is not significant. Without GSE244832 (40 donors), the CDH1 to PTPRM communication score falls from 0.496 to 0.0034, and ST6GAL1 to EGFR from −0.299 to −0.070. Both results rest on one dataset.

Donor-level variance of Ductular injury also differs by stage. In GSE244832, the residual SD rises from 0.08 in healthy donors to 1.25 in steatohepatitis. In the model with dataset- and stage-specific variances (E8), the stage beta was 0.172 (*p* = 0.013) for Ductular injury and 0.041 for Stromal ECM.

## Supplementary Note 9. Matched-background comparison for Hotspot program associations in bulk liver

In the 364-unit arm of Figure 4E with unweighted scores, 1 of 114 programs reaches fibrosis BH q < 0.05 (228-test family), against a mean of 1.24 for random gene sets matched on size and expression (empirical *p* = 0.473, 200 draws). This unweighted family differs from the weighted 224-test family of Figure 4E. No program reaches q < 0.05 for NAS, against a mean of 0.19 (*p* = 1). In the 457-record sensitivity arm, the counts are 20 against 26.9 (*p* = 0.871) for fibrosis and 2 against 2.57 (*p* = 0.687) for NAS. The programs can be useful for organizing results even though no more of them reach q < 0.05 for fibrosis or NAS than matched random gene sets do.

## Supplementary Note 10. Program-pair covariance in bulk liver and across Visium spots

### Scope

We report program-pair covariance in bulk liver as a discovery-cohort result, with replication in GSE281797 and GSE225740 reported separately for each cohort. Neither replication cohort was checked by genotype for participant overlap with the discovery cohorts. The composition-adjusted analysis in GSE281797 failed its prespecified acceptance rule and was not interpreted. Bulk program scores mix cell abundance with program activity, so these data do not show coordination between cell types. The covariance says nothing about timing or cause.

### GSE213621 test

The same edge sign-concordance statistic could not be read in GSE213621. The prespecified test used 258 GSE213621 disease units (genetic-identity clusters) and 2,328 pairs rediscovered on 371 identity-qualified discovery units, and it was untestable. Its program-group permutation control permutes GSE213621 units within program gene groups, which keeps covariance within each program and removes it between programs. This control still gave a positive mean excess over the non-coordinated baseline of 0.025 (95th percentile 0.081; to four decimals, 0.0251 and 0.0809; prespecified limits 0.02 and 0.05). The negative control therefore failed, and the run stopped before any GSE213621 outcome statistic was computed.

### Replication in GSE281797 and GSE225740

With 3 residual components removed, the excess over random placement was 0.187 (0.118 to 0.188) in GSE281797 and 0.241 (0.186 to 0.242) in GSE225740. These are percentile intervals from 1,000 participant-bootstrap draws within fibrosis stage. Most bootstrap excesses fall below the full-sample value (bootstrap medians 0.162 and 0.217; 4% and 3% of draws exceed the point estimate), so each point estimate lies near the upper end of its interval. Counts were normalized within each cohort; applied to GSE193066, this choice changed the fraction from 0.897 to 0.918. In GSE225740, the result held after adding three estimated composition components, log library size and the top-100-gene count share to the design. The GSE281797 composition estimate failed its prespecified quality rule and was not read. These numbers are not pooled with any GSE213621 or GSE193066 value.

The same pairs were also positive above the random-placement null among participants staged F0-F1 at biopsy. With 1 residual component removed, the excess was 0.166 in 66 GSE281797 participants and 0.304 in 53 GSE225740 participants, whose identities cannot be checked. The combined *p* was 0.000999 for each, and both strata were rejected under Holm correction over the two. This comparison is cross-sectional. In GSE281797, 29 of the 66 F0-F1 participants carry the source diagnosis "No Pathology", so the stratum result cannot be attributed to MASLD remodeling. Without them, the excess was 0.145 (37 participants). Fibrosis stage is already in the model, so the stratum result shows covariance within a narrow stage range. It does not rule out confounding by stage and says nothing about when the covariance arises.

### Discovery-set diagnostics and the GSE193066 comparison

A separate analysis, on a bulk sample set other than the 774 units, reported 2,426 of 6,195 program pairs beyond saturated histology adjustment, modularity 0.358 (against 0.138 for rank-matched and 0.089 for rewired nulls) and a 15-program fibroblast-and-cholangiocyte community. In its comparison in GSE193066, with the model fixed in advance, 0.8974 of coordinated pairs kept a positive sign (89.7% edge-sign concordance), against 0.7095 (70.9%) sign agreement for non-coordinated pairs. Modularity was 0.252, and community density was 0.876 against 0.681. The GSE193066 comparison uses the same excess over non-coordinated pairs as the GSE213621 test. This comparison used 105 source-declared GSE193066 participants, each at the first biopsy only, and their identities are not verified. These values are diagnostics only. They are not independent-cohort replication, and GSE193066 is excluded from primary replication, repeat-biopsy ICC and derived ceiling claims.

### Covariance across Visium spots

Across Visium spots, the 160 gene-disjoint coordinated cross-lineage pairs showed an excess of 0.0257 Fisher z, measured over a label null stratified by lineage pair. The excess remained after adjustment for lobular position, for the first spatial component and for program-specific effects.

## Supplementary Note 11. Equal-cohort pooled heatmap of the 43 molecular systems

The pooled heatmap of the 43 systems balances collections within each cohort and then gives each cohort one-fifth weight (Figure S4U). It shows the profile shared across cohorts without letting the largest cohort dominate. The adjacent panel shows the cohort-specific profiles, including discordant ones. Both panels use the same descriptive pooled late-minus-early row order. This order makes increasing and decreasing shared profiles easy to see, but it is not a support or timing statistic.

## Supplementary Note 12. Retrospective transfer of the continuum loading

GSE174478 and GSE240729 had already been evaluated, so they serve as retrospective transfer checks. These values come from a transfer of the 139-gene continuum loading built on an 844-record bulk set, and they were not recomputed on the 774 units. Sample-level correlations between the continuum and fibrosis were 0.510 (93 samples) and 0.710 (66 samples), with matched-background Holm *p* = 0.117 and 0.064 (1,000 matched random-signature draws), neither below 0.05. The prespecified primary test, a fixed projection within disease against matched random gene sets with Holm correction over GSE162694 and GSE213621 and both cohorts required to pass, missed: rho 0.547 against a null 95th percentile of 0.528 (109 samples) and rho 0.465 against 0.457 (293 samples), Holm *p* = 0.082 in each. Matched random gene sets give a mean of 86 continuum-associated programs against the 92 observed (*p* = 0.109). A positive association with histology therefore does not show that the loading performs better than matched random signatures. GSE240729 also uses FFPE tissue.

## Supplementary Note 13. RNA-protein agreement in Figure 5A

Across the 3,054 measured genes of Figure 5A, Spearman *rho* between the RNA and protein changes was 0.458, and RNA and protein changed in the same direction for 2,013 genes (65.9%). Between acquisition batches, protein-effect estimates agree at Spearman *rho* 0.34 over 3,034 genes. The batches contain different participants (28 MASLD and 5 controls; 18 MASLD and 7 controls), so this comparison mixes sampling error, case mix and technical variation and does not estimate an assay reliability ceiling. HLA-DQA1 (Figure 5C) is quantified with the primary-assembly annotation ENSG00000196735.14.

## Supplementary Note 14. Figure 7 severity instrument: fit and scores against task-native baselines

The released linear severity instrument was fitted to 521 sample records from five cohorts. This fit is separate from the 774-unit differential-expression fit of the Resource bulk analysis. All 521 records passed through the RNA-based genetic-identity screen. It found eight within-cohort duplicate pairs (KING kinship > 0.354; three in GSE130970 and five in GSE135251) and no duplicate across cohorts. The Figure 7B fibrosis comparison holds out one cohort at a time, so no duplicated library is in both training and test data. The duplicates remain within the sample bootstrap, which therefore slightly overstates the number of independent people.

The severity instrument scored 0.6631 against 0.6094 for the in-cohort ridge baseline (out-of-cohort Spearman correlation with Kleiner stage in 109 GSE268273 participants). GSE268273 was also used to choose the expression transform. The paired difference, +0.054, has a 95% interval of −0.025 to +0.137, which includes zero. The cis head predicts each H3K27ac region with a ridge regression that adds, to a global expression term, the expression of up to 10 genes whose transcription start sites lie within 50 or 100 kb of the region; window and gene number were chosen inside each training fold. It scored 0.1426 against 0.0480 for a ridge on 81 global expression components without cis genes (out-of-fold skill in 99 GSE267145 participants). This comparator is weaker than the rank-4 reduced-rank regression of Figure 7B (skill 0.122 in the same 99 participants). In a post hoc comparison on the same folds and targets, the cis head exceeded a reduced-rank regression with rank chosen in each fold (skill 0.119) by 0.023 (95% participant bootstrap 0.015 to 0.031). This is an upper bound, because the rank search reached its limit in 3 of 5 folds.

## Supplementary Note 15. Power of the comparison with Kumak et al.

These historical power checks use the earlier 3,769-gene universe and were not rerun after eligibility qualification. Planted shifts equal to an odds ratio of about 1.74 or 0.48 in the NAS list were detected in 82% and 65% of simulations, and shifts of about 1.75 or 0.48 in the fibrosis list in 77% and 62%.

## Supplementary Note 16. Disease accessibility, regulatory-link controls and recovered hepatocyte RNA

The two disease ATAC cohorts were fitted separately: GSE244832 has 5 normal and
9 MASH donors after excluding MASL, and GSE281367 has 6 normal and 6 MASH.
All four corrected primary TMM–voom quality-weighted fits reproduce to within
5.73 × 10^-13. Fragment-count effects and donor sums of binary nucleus
accessibility have different units and were not pooled. Robust edgeR
quasi-likelihood did not improve the fitted-negative-binomial null diagnostics;
those simulations are not biological ground truth, and their binary-count
versions cannot validate bounded-count calibration. The primary voom reference
was retained (Supplementary Figure S5AA).

Across 30,650 jointly testable hepatocyte and 9,264 stellate peaks, cross-cohort
Spearman correlations are 0.3897 and 0.2654, with sign agreement of 64.40% and
58.27%. These summaries describe dependent peaks, not independent donor tests.
Four stellate peaks meet same-direction BH q < 0.05 in both cohorts; no lineage
passes the separate three-lineage directional-replication family, which retains
untestable macrophages. All four peakwise candidates retain their effect sign
under each recorded donor omission in both cohorts; three retain BH q < 0.05
after recorded-batch adjustment in GSE244832. The chr2:64404164–64404664 peak
does not (adjusted q 0.343). These are within-cohort sensitivities, not further
independent replication. No program or chromVAR motif passes strict two-cohort
replication. Recorded-batch adjustment in GSE244832 changes hepatocyte calls
from 1 to 3 and stellate calls from 820 to 291 on identical donors and peaks.
The six batches include five with both conditions and one with three MASH
donors only. Condition-only and adjusted designs have ranks 2 and 7 and residual
degrees of freedom 12 and 7. Full rank is not complete separation of disease
and batch. Donor omissions, batch tables and all peak effects accompany the
displayed examples; none of these sensitivities changes the primary testing
family (Supplementary Figure S5P–T,W).

Separate disease-cohort distance and positive-Cicero annotations remain
descriptive. Cross-cohort Cicero top-link Jaccard is 0.1579 in hepatocytes and
0.1739 in stellate cells on common observability, compared with 0.9927 and
0.9985 for distance. Shared consensus peaks and top-ten selection qualify these
values. The complete peak-to-gene table preserves link method, score, coverage,
accessibility uncertainty, batch and donor sensitivity, and evidence state.
All 117 program memberships and weights remain fixed, with coverage exported
even when no qualifying link exists. Motif-family compatibility does not identify
TF activity. The illustrative AHNAK-linked stellate peak is a peakwise candidate,
not evidence that AHNAK is regulated by that peak (Supplementary Figure S5U–W).

The predictive benchmark is separate from these all-donor disease maps.
Original-grid results did not establish distal predictive benefit; the wider
grid was a secondary sensitivity prompted by tuning-boundary diagnostics, using
the same previously evaluated donors. In the separate post hoc matched-gene
comparison, adapted DIRECT-NET worsened hepatocyte RNA error by 0.246% versus
distance, while Cicero improved it by 0.306% (the original four-test Holm family
is unchanged). Complete original-grid references and wider-grid comparisons
accompany the new stress tests. In
GSE296875 hepatocytes, 39 held-out donors retain the original five outer/three
inner folds and fixed wider ridge grid. Baseline Cicero RNA error is 0.9770%
lower than promoter-only and 0.3062% lower than distance on matched genes.
Across aggregation seeds 43 and 44, the 500-nucleus cap and three deliberately
targeted training-donor omissions, reductions are 0.8394–1.0236% and
0.2869–0.3376%, respectively. Every unshuffled comparison retains a positive
gain over both baselines, satisfying the fixed internal robustness criterion.
Intervals use 2,000 seeded donor bootstraps conditional on saved fitted
predictions, without retraining uncertainty. No new confirmatory family was
introduced; original and post hoc Holm families remain separate and unchanged.
Genes unavailable to a method and tuning at the largest penalty are reported
explicitly rather than silently changing comparison universes.

For within-donor RNA pairing-shuffle seeds 43–45, paired-minus-shuffled Cicero
gain over distance is −0.000180 [−0.000434, 0.000087], −0.000120
[−0.000369, 0.000143] and −0.000226 [−0.000460, 0.000030] standardized MSE.
All conditional intervals include zero. Paired gains over promoter-only are
modestly larger than shuffled gains, but this does not establish a pairing-specific
Cicero advantage beyond proximity. Shuffles preserve donor marginals, so they do
not remove every donor-level chromatin–RNA association. GSE296875 is an internal
development cohort without a qualified MASH contrast. Neither independent MASLD
validation nor individual-edge correctness follows (Supplementary Figure S5X,Y).

Joint reconstruction of 117 sequencing runs as 18 GSE244832 capture libraries
recovers 149,096 called nuclei and 690,444,610 called-cell UMIs. Exact molecule
accounting passes per cell and per gene in every library. Source-native RNA
diagnosis, QC and lineage eligibility retain 5 normal and 9 MASH hepatocyte
donors; stellate and macrophage RNA have only 3 and 2 normal donors and are
untestable under the four-per-condition rule. Among the same 21,054 tested RNA
genes, 11 meet BH q < 0.05 in the condition-only model and 3 after recorded-batch
adjustment. CHST9 and TMEM150B are increased in both: their condition-only and
adjusted q values are 0.00187/0.00405 and 0.04069/0.00640, respectively. TRIM9
is significant only after adjustment (q 0.00389; condition-only q 0.11675).
These selected candidates are not an independently validated hepatocyte
signature. Recounting adds no ambient-RNA correction or doublet model, and
unresolved specimen matches preclude donor-matched RNA–ATAC associations
(Supplementary Figure S5Z). Complete effect and exclusion tables accompany
these panels; no Catalog class, program weight or causal-gene assignment changes.

## Supplementary Note 17. Held-out replication in SteatoSITE

**Cohort and processing.** SteatoSITE (Kendall et al., 2023; ENA PRJEB58625) is archival FFPE liver sequenced as stranded total RNA (2 x 50 bp). We used only the public reads and the NASH-CRN fibrosis stage in the ENA sample records. The source has no F0 code, so code 1 is read as F0–F1; unstaged samples are not used as controls. Of 658 samples, 652 libraries were counted. Quality control blinded to stage kept 515 donors (114 technical, 18 possibly mixed and 5 same-individual libraries removed); 493 are staged (F0–F1 255, F2 74, F3 75, F4 89). Inferred sex agrees with the sex recorded in an independent reprocessing of these data (LiRNA atlas) for all 494 shared included donors. Models are `~ batch + inferred_sex + stage` with nine sequencing batches. Specimen type (biopsy, resection, explant) is not public, and F4 includes explants. Stage-associated expression in SteatoSITE, including THRB, the risk score and deconvolution trends with stage, was reported by Kendall et al.; we test our own fixed effects and claim none of these as new.

**Controls.** COL1A1, EPCAM and LUM are higher and GNMT and CYP2C19 lower in F3–F4 than F0–F1 (all p < 10^-31). The slope of the SteatoSITE sex coefficient on discovery F3–F4 effects is -0.052 (95% CI -0.159 to 0.054), so sex and stage are separated here, unlike in GSE213621. Permutation-null mean slopes are -0.015 to 0.018. Adding RNA-quality covariates changes the F3–F4 slope from 1.462 to 1.452.

**Grouped contrasts (two tests, Holm).** Over the 8,460 measurable discovery F3–F4 calls, b = 1.462 (95% CI 1.323 to 1.601), 93.5% same sign (null 52.6%), Holm p = 0.002, minimum detectable slope 0.289. F2 versus F0–F1: b = 0.917, 95.4% same sign, untestable (minimum detectable slope 0.538). For F3–F4: up calls 1.422, down calls 1.171, women only 1.315, F3 only (F4 removed) 1.134, and with six cell-composition coordinates 0.256. The composition-adjusted fit is compared with an unadjusted discovery vector, so the drop is descriptive; it indicates that most of the bulk stage signal moves with tissue composition.

**Adjacent steps (three tests, Holm).** F0–F1 to F2 (discovery F1 to F2) b = 0.952, 94.4% same sign; F2 to F3 b = 0.691, 96.1%; F3 to F4 b = 0.705, 94.6%. All three are untestable (minimum detectable slope 0.36 to 0.59) and are descriptive.

**Genes.** GNMT, MAT1A and CYP2C19 fall with stage code (log2 per code -0.569, -0.344, -0.718; Holm over three, all p < 10^-35). IGFBP7 rises from F2 to F3 (+0.666 log2, one-sided p = 7.9 x 10^-16). Of the 29 later-step Figure 3F marks, 28 pass BH q < 0.05 within 29 and all 29 have the discovery sign; PSMC3IP (F1 to F2) has q = 0.21. Of 10,903 measurable stage-supported genes, 92.5% have the discovery trend sign (null 52.0%; descriptive, because these genes were defined with GSE213621 and the continuum).

**Gene Catalog field.** Of 8,542 Catalog genes that are discovery F3–F4 calls: 7,214 replicated, 645 underpowered, 383 not replicated at adequate power, 201 opposite sign and 99 not measurable. No class changes.

**Continuum and programs.** The fixed continuum projection correlates with stage at Spearman rho = 0.776 (493 donors); expression-matched random gene sets give median 0.667 and 95th percentile 0.739 (p = 0.001), so most of this correlation is shared by random expressed genes. Hepatocyte-fraction partial rho is 0.495. Stromal ECM and Ductular injury rise per stage code (0.674 and 0.633 SD, HC3, Holm over two, both pass; 0.158 and 0.091 after composition adjustment). Of 112 testable programs, 98 change with stage at BH q < 0.05, and all 63 stage-supported programs agree in sign. This separate arm does not change the Figure 4E result.

**Composition.** Within stage, portal-tract marker content explains a share of 0.506 (95% CI 0.451 to 0.562) of donor variation in continuum position, above the 95th percentile of its null (0.236). The IgA share of immunoglobulin heavy-chain fragments rises 0.148 logit per stage code (SE 0.023, p = 1.1 x 10^-10; 0.119 with plasma-cell and immunoglobulin-rate terms). This is a shift in antibody isotype, not a plasma-cell count.

**Genotypes and allelic expression.** PNPLA3 rs738409 calls from RNA fail Hardy-Weinberg equilibrium (p = 8 x 10^-8; 147 heterozygotes observed, 197 expected). The heterozygote deficit is present in every stage group and is unchanged by duplicate removal, so the prespecified genotype analyses were not run. Allelic imbalance by stage: 0 of 220 testable genes at BH q < 0.05.

**List overlap within one cohort.** In 50 random splits of the 493 donors into five groups of about 99, 26.8% of genes called in any group (F3–F4 versus F0–F1, BH q < 0.05 within group) were called in one group only, and 0.17% of genes called in two or more groups had mixed signs. Across the five Figure 3C cohorts the corresponding values are 68.5% and 31.5%. For disjoint pairs of 109 donors, a median of 76% of one group's calls were called in the other. The contrast differs from Figure 3C, and the groups share protocol and batches. These are descriptive within-cohort comparisons, not bounds on agreement or a decomposition of between-study heterogeneity.

**Pooled sensitivity (descriptive; Figure 3 stays discovery-only).** We refitted the Figure 3E adjacent steps with the 387 discovery units and the 493 SteatoSITE donors together (`~ dataset + inferred_sex + group`, each SteatoSITE batch a dataset level; BH over 24,115 genes; SteatoSITE code 1 enters only the grouped F0–F1 to F2 step). The discovery-only refits reproduce Figure 3E exactly (maximum logFC difference 2 x 10^-14; 512 and 281 calls). Pooled calls are 7,239 for F0–F1 to F2 (discovery 2,922), 3,943 for F2 to F3 (512) and 9,475 for F3 to F4 (281). Of the pooled calls, 99.1%, 98.4% and 86.9% have the discovery sign, and every discovery call that stays called in the pooled fit keeps its sign. SteatoSITE carries a median 50%, 47% and 64% of the voom weights in the three pooled fits, and alone it gives 5,180, 3,755 and 9,857 calls, so the pooled counts mostly reflect SteatoSITE's larger later-stage groups. A stage-by-source interaction is called for 6, 2 and 351 genes, so effect sizes differ between the sources mainly at F3 to F4, where SteatoSITE F4 includes explants. Without GSE135251, the discovery F2 to F3 step keeps 2 genes (54 units); the pooled fit without GSE135251 keeps 3,327 (203 units, 66% SteatoSITE weight). The single-cohort dependence of this step in Figure 3E therefore reflects the small number of units left, not a step absent from other cohorts [Inference]. The discovery-only F3 to F4 fit without GSE135251 (29 units, 11 at F4) gives 6,565 calls; with so few units we do not interpret it.

**Exonic and intronic change (exploratory).** Total-RNA libraries retain intronic reads, which mark pre-mRNA. For 12,010 genes with intronic signal, the F3–F4 versus F0–F1 change in intronic reads followed the change in exonic reads closely (Spearman 0.95; slope 0.92 over the 4,346 replicated Catalog genes that could be measured; Figure S3 [panel letter to be assigned]). In the corrected batch- and inferred-sex-adjusted model, the exon/intron ratio (Gaidatzis et al., 2015) changed for 2,711 genes at BH q < 0.05 over 12,010 (global call-count permutation p = 1/201 = 0.004975; 200 shuffles within batch and sex, seed 20261103). The median absolute ratio change among these calls was 0.152 log2; 114 genes changed by at least 0.5 log2, approximately 1.41-fold. The historical stage-only model gave 3,008 calls; 2,525 are shared, 483 lose a call and 186 gain a call after adjustment. Whether more of these ratios rise or fall depends on the normalization library set, so we do not interpret the direction balance. For GNMT, MAT1A, COL1A1, IGFBP7, TNFSF10 and LUM the intronic change was 90–113% of the exonic change. The aligned intronic and exonic effects are consistent with a substantial transcription-associated tissue component. Changes in cell mix also move both read types together. Intronic abundance depends on RNA processing, and separate normalization assumes that sample-wide exon/intron balance changes are technical. These data do not resolve absolute synthesis or degradation rates or exclude concurrent stability effects. The prespecified edgeR fit with donor terms was replaced by a limma test on each donor's log exon/intron ratio because the edgeR fit was too slow to permute; the historical comparison of 300 genes used the stage-only models and gave Pearson 0.94 between effect estimates. It is not a comparison of the current adjusted model and does not establish equivalence of uncertainty estimates or tests.

**Figure 3E axes (exploratory).** The fixed k6 NMF bases were applied to SteatoSITE with the release transport, which reproduces the discovery scores exactly. A per-batch shift was needed for the new protocol, and scores were almost unchanged without it (within-batch Spearman ≥ 0.993). The F3–F4 versus F0–F1 activity had the discovery direction for P1 Inflammatory-EMT (1.34 SD), P2 Hepatocyte-Metabolic (−1.41), P3 Fibrotic-ECM (1.57) and P4 Unresolved (−1.31) (one-sided permutation tests, Holm over five axes, all p = 0.005), and the F2, F3 and F4 activities of these axes all had the discovery sign and order (Figure S3 [panel letter to be assigned]). P5 Noncoding is not interpretable: its discovery F3–F4 effect is near zero (−0.01), and its genes depend on total-RNA versus poly(A) library preparation. Like the gene-level slope, these axis shifts likely include cell-composition change.

**Liver eQTL direction (exploratory).** At expressed SNPs that are GTEx v8 liver significant variants for the same gene, the allelic read share in heterozygous donors had the GTEx direction for 34 of 50 testable genes (68%; random-sign null 50%, one-sided p = 0.008), so these liver eQTL directions are visible in diseased FFPE liver.

## Supplementary Note 18. Follow-up analyses of the colocalized genes

Each analysis below had a written plan, checksummed before it ran. Analyses written after an earlier result was seen are labelled post hoc. None changes a Catalog class.

The saved-signal inventory contains 383 unique unordered gene-study comparisons,
including 94 recurrent comparisons representing 78 genes. Of these, 88
comparisons (73 genes) reuse the same recorded eQTL signal index. That index
does not establish a matched GWAS causal signal, and the shared eQTL panel is
not independent molecular replication. The 342 directed tested rows below
answer a different counting question and are not 342 independent genes or
cohorts. These historical recurrence estimates also precede the variant-identity
qualification of 90 primary genes; they cannot certify a replicated genetic core.

**Gene-level recurrence in another GWAS of the same trait.** We used three European pairs of independent cohorts: MVP and UK Biobank ALT, MVP and UK Biobank AST, and MVP and FinnGen NAFLD. Other European NAFLD and liver-fat GWAS share UK Biobank or FinnGen participants with these and were not used. For each gene-by-GWAS colocalization (SuSiE PP.H4 > 0.5) in one member of a pair, the partner was scored as recurrent (PP.H4 > 0.5), distinct variant (PP.H3 > 0.5), neither, or untested (no SuSiE result). A gene takes its best state over its rows. Of the 436 genes, 78 have recurrent gene calls, 121 have a distinct variant in the partner, 12 neither, 103 are untested and 122 have no row in a pair (UK Biobank GGT 80, Biobank Japan 30, other 12). Of the 22 direct-trait genes, GPAM, RORA, OSGIN1, TOR1B and C2orf16 recur, all through enzyme GWAS pairs; the direct NAFLD pair was testable for one gene. Gene-level recurrence among the 342 directed tested rows, with 95% intervals from 10,000 bootstrap draws of 147 locus clusters (anchors on one chromosome within 500 kb; seed 20261102):

| Discovery PP.H4 | Tested rows | Recurrent gene call | Distinct variant | Neither | Fraction (95% interval) |
|---|---|---|---|---|---|
| 0.5 to 0.8 | 63 | 17 | 42 | 4 | 0.27 (0.16 to 0.39) |
| 0.8 to 0.95 | 99 | 51 | 45 | 3 | 0.52 (0.38 to 0.64) |
| ≥ 0.95 | 180 | 120 | 53 | 7 | 0.67 (0.57 to 0.75) |
| All | 342 | 188 | 140 | 14 | 0.55 (0.47 to 0.63) |

Every row involves MVP summary statistics, which were fine-mapped with an out-of-sample LD reference. A distinct-variant call can come from different causal variants in the two cohorts or from the different LD references, and this analysis cannot separate the two. The plan stated in advance that all single-study genes would be labelled provisional if the fraction stayed near 0.5 at PP.H4 ≥ 0.95; it is 0.67, so that rule did not apply. Genes with discovery PP.H4 0.5 to 0.8 remain provisional. Recurrence alone does not qualify variant identity or demonstrate the same mechanism.

**Coding variants of neighbouring genes.** For each of the 840 colocalized tier-1/2 signal pairs (767 gene-by-GWAS pairs), we took the GWAS-side lead variant and computed LD (1000 Genomes panel of the GWAS ancestry) with protein-altering or canonical-splice credible-set members of tier-1/2 studies at PIP ≥ 0.1. Gene-level flags, taking the most severe class over a gene's pairs:

| Class | All 436 | Direct 22 |
|---|---|---|
| Anchor is another gene's coding variant | 1 (C2orf16; GCKR rs1260326) | 1 |
| Not separable, r² ≥ 0.8 | 3 (CWF19L1 with ERLIN1, r² 0.92; RECQL4 and ZFTRAF1 with LRRC24, r² 0.85) | 1 |
| Moderate LD, r² 0.2 to 0.8 | 5 (IFT172, LRRC14, LRRK2, SUPT7L, KCTD10) | 0 |
| Anchor is or tags the gene's own coding variant | 3 (PTTG1IP, EXOC3L4, WIPI1) | 0 |
| Not assessed (anchor absent from the LD panel) | 29 | 1 (MAU2) |
| None | 395 | 19 |

At any PIP the counts are 1, 5, 10, 6, 47 and 367. The per-gene lead variants of Figure 2C are a different summary, and for several genes they differ from the signal-pair anchors used here. SAMM50's Figure 2C lead is PNPLA3 rs738409, but its signal-pair anchor (22:44325748, GRCh37) has r² 0.05 with rs738409 (D′ 1). The MAU2 anchors have r² < 0.2 with TM6SF2 rs58542926 where computable. The IFI27L2 and SERPINA11 anchors lie 250 kb and 69 kb from SERPINA1 rs28929474 with r² < 0.2. The HSD17B13 variant rs72613567 is a splice-region insertion, outside the coding categories used, so SPP1 was not testable by this rule. The 29 not-assessed genes are not clear of a coding neighbour.

**GNMT and the European enzyme GWAS.** The GNMT direct-NAFLD ABF result reproduces (PP.H4 0.706 on 3,384 variants; anchor 6:42931261). The five enzyme GWAS leads in the GNMT eQTL window lie at 6:43.74 to 43.81 Mb, 810 to 880 kb from the anchor, with r² < 0.005 in the European panel (379 samples). Conditioning NAFLD on each lead moves the anchor z from 3.65 to 3.56 to 3.63 and leaves PP.H4 at 0.64 to 0.69; conditioning each enzyme GWAS on the NAFLD anchor leaves its lead unchanged. Because the leads are at a different locus, this test could not have removed the NAFLD signal and is uninformative about the local haplotype. In a post hoc analysis with the eQTL window restricted to 250 kb around the anchor, chosen after the lead positions were seen:

| GWAS | z at anchor (p) | Full-window PP.H3 | Local PP.H3 | Local PP.H4 |
|---|---|---|---|---|
| NAFLD (Ghodsian et al., 2021) | 3.65 (2.6 × 10⁻⁴) | 0.051 | 0.032 | 0.720 |
| MVP ALT | 0.39 (0.70) | 1.000 | 0.099 | 0.001 |
| MVP AST | 2.26 (0.024) | 1.000 | 0.228 | 0.009 |
| UK Biobank ALT | −6.67 (4.3 × 10⁻¹¹) | 0.998 | 0.884 | 0.116 |
| UK Biobank AST | −2.33 (0.022) | 1.000 | 0.060 | 0.011 |
| UK Biobank GGT | 3.46 (6.6 × 10⁻⁴) | 0.939 | 0.054 | 0.219 |

UK Biobank ALT has a local association (smallest p 1.5 × 10⁻¹³ at 6:42944988) that does not colocalize with the eQTL; MVP ALT has none at the anchor. The UK Biobank summary statistics carry beta and standard error at low precision, so their z values are approximate. GNMT's colocalization stays sub-threshold and provisional for the reasons in the main text.

**Allele directions at three-way loci.** At loci where a GWAS signal colocalizes with both a liver eQTL and a bulk-liver accessibility QTL (Currin et al., 2025), we aligned all three effects to one allele. The caQTL effect is for the last allele of the variant identifier (998 of 998 matchable rows of Currin Table S23A agree). Palindromic variants were dropped, leaving 88 of 103 loci, one row per locus (highest three-way PP.H4). The accessibility-raising allele raised expression at 67 of 88 loci (0.761; binomial Holm p = 1.8 × 10⁻⁶ over two tests): 60 of 75 liver-trait loci and 7 of 13 platelet-count loci. The 21 opposite-direction loci (CYP4V2, PPP2R1B, S100A13, SLC2A2, IL15RA, SPATA13, PPP1R3B, LGR6, PPA1, ENAH, SEC14L2, EEF1AKMT3, RFC2, NAP1L5, KLC1, IL18R1, H1-0, NOL8, ATP5MG, ATP5IF1, SLC25A11) are not evidence against colocalization. At 49 tier-1/2 loci whose gene is a disease-state call, the trait-raising allele moved expression in the disease direction at 28 (0.571; Holm p = 0.39; Fisher OR 1.78). Forty-eight of the 49 are liver-enzyme loci.

## Supplementary Note 19. Cross-cohort consistency of the disease-state calls

**Current transportability sensitivity.** The complete 24,115-gene family was
analyzed by REML with modified Hartung–Knapp uncertainty and a variance-factor
floor of one (Methods). The native GSE126848 filter leaves 15,592 genes with
five cohorts and 8,523 with four. Among the fixed 1,854 pooled calls, 1,056 are
testable in GSE126848 and 798 are untestable there. Random-effects means retain
the pooled sign for 1,850/1,854 genes; 1,208 agree in every available cohort.
No gene passes BH over 24,115 tests (minimum q=0.171430). Median I-squared among
the fixed calls is 60.35%; 908 nominal Q p values below 0.05 are descriptive.
These results show broad directional agreement with unresolved small-cohort
mean-effect inference, not absence of tissue effects or independent replication.
The fixed pooled calls and separate SteatoSITE stage results are unchanged.

**Historical consistency classes.** The earlier five-cohort DerSimonian–Laird
screen retained all 1,854 selected genes in every cohort, including entries
excluded by the current GSE126848 native filter. It gave 1,403 random-effects
z-test calls and 902 heterogeneity calls under its original correction. The
table below preserves those descriptive classifications; it does not describe
the current REML/modified Hartung–Knapp inference. Input eligibility, testing
family and uncertainty rules differ, so their numerical contrast cannot be
attributed to one change alone. The label "Consistent" is historical and does
not certify replication.

| Class | Up | Down | Total | Median τ | Median I² |
|---|---|---|---|---|---|
| Consistent | 536 | 210 | 746 | 0.16 | 0.35 |
| Single-cohort dependent | 467 | 184 | 651 | 0.32 | 0.64 |
| Direction consistent, size heterogeneous | 259 | 99 | 358 | 0.51 | 0.79 |
| Cohort-discordant (two or more cohorts opposite) | 71 | 28 | 99 | 0.58 | 0.77 |

In the historical classification, consistent calls had at most one cohort with the opposite sign and a same-sign random-effects call at BH q < 0.05 that survives removal of any one cohort; 731 of the 746 have random-effects |log2FC| > 0.5. Single-cohort dependent calls lose the random-effects call when one cohort is removed. GSE213621 is that cohort, alone or with others, for 625 of the 651 (alone for 261), and GSE162694 is involved for 301. GSE213621 has 326 of the 774 units, so part of this dependence reflects power, but these calls are not replicated in the four smaller cohorts. The class, τ, I², the prediction interval and the cohorts whose removal loses the call are reported per gene.

**GSE48452.** Deposited RMA values, with after-surgery biopsies excluded, from 12 lean, 16 obese normal-histology, 9 steatosis and 17 NASH participants, one biopsy each. 1,036 of the 1,854 calls matched by symbol (868 of 1,333 up, 168 of 521 down). Each group was compared with lean controls by limma (empirical Bayes with an intensity trend; sex in the model). The statistic was the least-squares slope of array log2 fold change on discovery log2 fold change over the matched calls, with a one-sided p from 1,000 label permutations (seed 20261102) and Holm correction over three contrasts:

| Contrast versus lean | Slope | Null 95th percentile | Holm p | Sign agreement (up, down) |
|---|---|---|---|---|
| Obese, normal histology | −0.006 | 0.093 | 0.56 | 0.397 (0.374, 0.512) |
| Steatosis | 0.069 | 0.125 | 0.34 | 0.576 (0.570, 0.607) |
| NASH | 0.336 | 0.167 | 0.003 | 0.836 (0.820, 0.917) |

The obese slope divided by the NASH slope was −0.017 (participant bootstrap 95% interval −0.607 to 0.305; 1,000 resamples within group). Four of the 1,036 calls had an obese-versus-lean one-sided BH q < 0.05. Adding age gave slopes of −0.0005, 0.101 and 0.351. The comparison is cross-sectional, in one microarray cohort with 12 lean and 16 obese participants.

**GSE281797.** 91 participants with sex concordant between the record and expression: MASL or MASH 62, controls with no pathology 29 (26 female, 3 male; 6 controls have fibrosis stage F1). The cohort was reported by Zhao et al. (2026; PMID 41870035). Of the 1,854 calls, 1,063 were testable (874 up, 189 down); 573 failed the expression filter and 218 were not annotated. Statistics were compared with 1,000 label permutations within sex (seed 20261102), with Holm correction over two statistics:

| Statistic | Observed | Null median | Null 95th percentile | Holm p |
|---|---|---|---|---|
| Sign agreement | 0.653 | 0.500 | 0.702 | 0.134 |
| Slope of replication on discovery log2FC | 0.209 | 0.006 | 0.138 | 0.016 |

Sign agreement was 0.598 for up calls and 0.905 for down calls (slopes 0.158 and 0.249). A per-gene binomial test of sign agreement gives p = 6 × 10⁻²⁴, but it treats genes as independent, which the permutation null shows they are not. No call reached one-sided BH q < 0.05 (0 of 1,063). The plan required sign agreement to pass before the slope was read, so direction replication is not shown in this cohort; the slope alone is above its null. Sensitivities gave sign agreement 0.543 to 0.706 and slope 0.177 to 0.415 (all 94 participants with deposited sex; F0-only controls; adding age, BMI and diabetes; MASH versus control). Post hoc, by consistency class, sign agreement was 0.724 for consistent calls (490), 0.608 for single-cohort dependent (360), 0.582 for size-heterogeneous (165) and 0.500 for cohort-discordant calls (48), with no test.

**Per-hepatocyte abundance of GNMT, MAT1A and CYP2C19.** In hepatocyte pseudobulk (TMM log2 CPM) from GSE202379 (36 donors) and GSE244832 (18 donors), we regressed each gene on diagnosis category coded 0, 1, 2 (healthy, steatosis, steatohepatitis) with sex, log10 nuclei and log10 UMI per nucleus, using HC3 standard errors and Holm correction over six tests:

| Cohort | Gene | Slope | p | Holm p |
|---|---|---|---|---|
| GSE202379 | GNMT | +0.165 | 0.68 | 1 |
| GSE202379 | MAT1A | −0.222 | 0.065 | 0.39 |
| GSE202379 | CYP2C19 | −0.720 | 0.11 | 0.53 |
| GSE244832 | GNMT | −0.963 | 0.19 | 0.76 |
| GSE244832 | MAT1A | −0.176 | 0.60 | 1 |
| GSE244832 | CYP2C19 | −0.723 | 0.23 | 0.76 |

Disease-versus-healthy contrasts, removing the UMI covariate and an alternative label set gave the same reading. Including cirrhotic donors in GSE202379 (45 donors) gave CYP2C19 −0.555 (nominal p = 0.011) and GNMT −0.25 (p = 0.46). The slope here is mostly steatosis against steatohepatitis, with 3 to 5 healthy donors per cohort, so a per-hepatocyte loss limited to advanced fibrosis is not excluded. The bulk and protein decreases are unaffected.

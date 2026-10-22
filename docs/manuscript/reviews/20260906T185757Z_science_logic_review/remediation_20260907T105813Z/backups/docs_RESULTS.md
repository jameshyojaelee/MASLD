# Current results and number ledger

**Verified:** 2026-08-17 EDT  
**Continuum candidate verified:** 2026-08-21 EDT  
**Rule:** every number below is labeled by release state. Do not combine release
states silently.

## Release labels

| Label | Meaning |
|---|---|
| **Current official** | Read-only July release; usable when that release is named |
| **Validated update** | Recomputed and checked; intended for the Resource but not yet adopted with every dependent output file |
| **Validated candidate** | Complete, checked analysis retained for candidate figures or narrative integration; not an official release until synchronized adoption |
| **Preflight only** | Source/design audit completed; inferential model not run |
| **Pending rebuild** | The corrected COLOC compute is complete and adopted; these derived values stay unauthorized until they are regenerated against the adopted release |
| **Supplementary** | Valid boundary or null analysis; not a main-paper discovery |

## Continuum translation extension: validated candidate, 2026-09-06

| Quantity | Candidate result | Unit and interpretation |
|---|---:|---|
| Complete Gene Catalog sidecar | 28,234 | Unique gene identifiers; includes all 23,370 bulk-tested genes plus genetics-only genes |
| Retained replay evidence | 5,362 signal pairs; 22,356,697 variant records | All exported pairs from 816 main-trait gene–study records; lower-support pairs remain visible |
| Supported main-trait signal pairs | 925 | SuSiE PP.H4 >0.5; may include several signals per gene–study record |
| Direction-resolved supported signal pairs | 83 | Single exported signal on each side and ≥95% original shared posterior mass in one resolved marginal direction; no conditional or drug effect |
| Stage/sex-supported programs | 93/117 | Both co-primary meta-analysis BH q<0.05 and all four cohort/axis directions |
| Composition-adjusted supported programs | 68/117 | Same rule and participants, adding 15 lineage log-ratio coordinates |
| Composition-adjusted separate-cohort support | 26/117 | BH q<0.05 separately in both evaluation cohorts and both axes, with consistent direction |
| Evaluation participants | 109 + 361 = 470 | GSE162694 and GSE213621; complete denominators per testable program |
| Matched liver/plasma participants | 41 | 11 controls and 30 with MASLD; existing public cohort, unpaired with bulk RNA |
| Matched shared-protein family | 345 | Six untestable members retained in the BH denominator |
| Adjusted matched protein associations | 49/345 | BH q<0.05 after fibrosis, age, BMI, sex and acquisition-batch adjustment; not liver origin or treatment response |

Direction categories overlap at gene level. The liver-fat results include
opposing directions across two PDFF studies for the same gene; individual
signals are retained without a common directional recommendation. No
replicated direct-MASLD opposing-direction hero is asserted. Exact tables,
figure and limitations are linked in the
[technical review](technical/continuum_translation_20260906T135538Z/IMPLEMENTATION.md).
These are candidate outputs, not changes to the adopted genetic memberships,
evidence classes, or main figure selection.

## Human bulk RNA-seq

| Quantity | Current official release | Validated Resource update | Interpretation |
|---|---:|---:|---|
| Count-matched human samples | 1,281 | 1,281 | Nine retained cohorts in the corrected source registry |
| QC-passing human samples | 1,260 | 1,257 | Complete source substrate; three current metadata rows lack a count match |
| Cohorts in pooled disease-vs-control analysis | 5 | 5 | Same intended five-cohort comparison |
| Pooled samples | 846 | 844 | Updated analysis uses fragment-counted libraries and removes two failed samples |
| Controls / disease | not re-ledgered here | 157 / 687 | Biological groups in the validated update |
| Figure 3D PCA participants | 328 strict extremes | 418 | Endpoint visualization uses 115 strict controls and 303 advanced-disease participants across four cohorts; the 844-participant PCA is supplementary |
| Genes tested | 27,638 | 23,370 | Universe changes under corrected fragment counting and filters |
| Canonical DEGs | 1,853 | 1,347 | `padj < 0.05` and `|log2FC| > 0.50` (canonical from 2026-08-12) |
| Up / down DEGs | 1,438 / 415 | 1,003 / 344 | Direction by disease-vs-control log fold change |
| TREAT DEGs (sensitivity arm) | 1,918 | 1,616 | `treat_fdr < 0.05` with interval null `|logFC| ≤ 0.25`; canonical 2026-06-29 to 2026-08-12 |

**Paper decision:** the updated five-cohort fragment-count analysis is the
intended bulk substrate. Its values remain labeled “validated update” until
stage, composition, figures, Gene Catalog fields, and release ledgers are rebuilt and
adopted together.

The validated endpoint sensitivity does not alter that estimand. Advanced
disease is `NAS >= 5` or fibrosis `>= F3`; strict controls are source-defined
controls with no observed NAS above 1 and no observed fibrosis above F0.
Histology takes precedence for the one conflicting source label. Requiring at
least five participants per arm retains 418 participants (303 advanced, 115
strict controls) across four cohorts. The endpoint and canonical log2 fold
changes have Spearman rho 0.938 across all 23,370 genes and 89.0% direction
agreement. The independently thresholded endpoint analysis calls 2,643 genes; 1,307
of 1,347 canonical genes are recovered (97.0%), with Jaccard overlap 0.487 and
1,336 endpoint-only calls. The interval-null TREAT sensitivity recovers 96.3%
of canonical TREAT calls. These are validated candidate results, not a change
to canonical Catalog membership.

Figure 3C uses all five canonical cohorts and displays 15,090 unique genes
meeting the per-cohort FDR <0.05 and |log2 fold change| >0.5 criteria, shown as
21,995 cohort memberships. Of these genes, 9,777 are unique to one cohort;
4,082, 933, 235, and 63 recur in two, three, four, and all five cohorts,
respectively. Figure S3C retains the robustness analysis: among the 1,347
canonical genes, per-cohort direction agreement ranges from 79.8% to 100%;
leave-one-cohort-out all-gene effect correlations range from 0.867 to 0.992 and
canonical-threshold retention ranges from 60.7% to 98.0%. Participant is the
inferential unit in every fit.

## Recurring biological examples

These examples were chosen to explain biological questions, not because they
maximize modality count. The disease-state values below are validated updates
awaiting synchronized adoption. The pooled RNA model uses 844 participants
and fixed cohort and inferred-sex terms; BH adjustment covers all 23,370 tested
genes. Adjacent-stage models use 209, 195, 165, and 100 participants for F0→F1,
F1→F2, F2→F3, and F3→F4, respectively, with fixed cohort and inferred-sex
terms and a separate 23,370-gene BH family per contrast. PXD051911 protein
effects use 58 patients and adjust for acquisition batch, age, BMI, and sex;
protein q values use BH correction across all 6,071 tested proteins.

### Loss of hepatocyte metabolic function: GNMT–MAT1A–CYP2C19

| Gene | Five-cohort disease versus control | Supported adjacent fibrosis contrast(s) | Adjusted liver protein |
|---|---|---|---|
| GNMT | log2FC -0.813; canonical positive; TREAT q=8.82e-11 | F0→F1 -0.424, BH q=0.00854; F1→F2 -0.378, q=0.0475 | log2FC -1.136, BH q=0.00642 |
| MAT1A | log2FC -0.433; below the canonical 0.50 floor but TREAT q=1.86e-6 | F0→F1 -0.244, BH q=1.28e-4 | log2FC -0.473, BH q=1.22e-4 |
| CYP2C19 | log2FC -1.438; canonical and TREAT positive | F1→F2 -0.686, BH q=0.0404 | log2FC -0.459, BH q=0.0255 |

The RNA, stage, and protein directions agree for all three genes. This supports
a recurring example of coordinated loss of hepatocyte metabolic function in
established disease. It does not establish one upstream mechanism, a discrete
patient subtype, longitudinal change, or an inherited regulatory program.
Their version-fixed single-cell program memberships are tested but nonsignificant,
and their spatial rows are coverage context rather than predeclared gene-level
inference.

An exploratory native DIA-MS audit in 46 MASLD participants further tests
whether the three proteins mark a coherent state. Rank-residual partial
Spearman correlations adjusted for acquisition batch, age, BMI, and sex are
0.472 for GNMT–MAT1A, 0.406 for GNMT–CYP2C19, and 0.472 for
MAT1A–CYP2C19; BH q values over the three-pair family are 0.00141, 0.00518,
and 0.00141. All three remain supported after adding fibrosis stage, including
global BH correction over all 12 pair-by-adjustment tests, but none remains
supported after adding NAS or all four histology components. This candidate
audit supports a histologic-severity-linked metabolic state, not an autonomous
three-gene mechanism. It is not part of the adopted release.

## Where the three aspects and variant combinations currently stand

This subsection collects what the release already contains for the two
goal-level dimensions added to `docs/PAPER.md`, so the numbers stay in one place.
Neither dimension carries a claim.

**Aspect resolution.** Recorded fibrosis stage is the severity anchor. NAS enters
as a graded sum in the adjacent-NAS contrasts (0 to 1-2, 1-2 to 3-4, 3-4 to 5-8)
tabulated elsewhere in this document, and as an adjustment set in the candidate
proteomics audit: the GNMT-MAT1A-CYP2C19 state holds after adjusting for fibrosis
stage, including global correction over all 12 pair-by-adjustment tests, but does
not remain supported after adding NAS or all four histology components. That is a
recorded result about the components as covariates. The three aspects have not
been resolved as separate outcome axes anywhere in the release, so the evidence
state is `indeterminate` and no aspect-specific attribution is available.

**Variant combinations.** Fine-mapping partitions credible-set mass by
consequence class rather than assuming inherited variation is regulatory. One
two-variant local haplotype is characterised at GNMT, reported as an unresolved
local haplotype rather than a resolved causal variant. No systematic treatment of
variant combinations exists beyond that locus, and the colocalization design is
single-variant throughout.

**Chromatin.** Accessibility is evaluated across jointly testable peaks, with the
large majority resolving to `indeterminate` and a minority to `source_dependent`,
alongside the variant accessibility heatmap and the six-lineage accessibility
matrix. Chromatin supplies assay-native context for placing candidates; no
epigenetic effect or causal claim is available from it.

### Prespecified programs with different spatial organization: Stromal ECM and Ductular injury

The corrected source-diagnosis candidate uses 64 biological donors across five
datasets, tests all 117 programs in one BH family, and adjusts for dataset.
Stromal ECM has disease β=0.399, primary q=0.0496, and HC3 q=0.105;
Ductular injury has β=0.412, primary q=0.0187, and HC3 q=0.105. Both
pass the primary BH/stability rule, but neither passes HC3 family correction.
IGFBP7 is not a canonical pooled DEG (log2FC 0.390) but increases at F2→F3
(0.428, BH q=0.00104). BICC1 is a canonical pooled DEG (log2FC 0.586) and
increases at F1→F2, F2→F3, and F3→F4 (0.422, 0.484, and 0.825; BH
q=0.0126, 0.0143, and 0.0386). Matched spatial nulls support the Stromal ECM
program in GSE192741, whereas the Ductular injury program is indeterminate.
IGFBP7 is one of 49 fixed Stromal ECM genes and contributes 0.02474 (2.47%)
of its L1 weight. BICC1 is one of 22 Ductular injury genes and contributes
0.05436 (5.44%). These are program exemplars, not gene-level spatial tests or a
proposed IGFBP7–BICC1 molecular interaction.

The aligned ambient candidate freezes all 117 memberships, weights, and scoring
universes and uses the current source-diagnosis metadata. decontX succeeds for
GSE136103, GSE174748, GSE185477, GSE202379, and GSE244832. GSE189600 and Liver
Atlas fail because their ambient batches contain only one contamination
cluster; they are excluded from both sides of the primary common-universe
comparison and appear as `uncorrected_passthrough` only in the full-universe
sensitivity. Raw reconstructed scores pass r=0.99893–1.00000 for all 117
programs. Integer-count transport is valid for 106 programs. The 11 version-fixed
T-cell programs are `untestable` because their native non-integer scoring
substrate cannot be reproduced from complete-atlas integer counts.

In the 62-donor, five-dataset primary common universe, Stromal ECM changes from
β=0.401 (BH q=0.0717; HC3 q=0.123) to β=0.130 (BH q=0.312; HC3 q=0.296). Its
paired correction-induced change is β=-0.271 (BH q=0.271; HC3 q=0.215), and its
median version-fixed L1-weighted ambient burden is 0.810 across corrected datasets
(range 0.689–0.929). Ductular injury changes from β=0.396 (BH q=0.0501; HC3 q=0.123) to
β=0.389 (BH q=0.129; HC3 q=0.162). Its paired change is β=-0.00747 (BH and HC3
q=1.00), with median ambient burden 0.342 (range 0.235–0.685). Both states are
`attenuation_indeterminate`: Stromal ECM attenuates more, but neither paired
change passes the 117-program family correction. The full-universe passthrough
fit preserves the prior raw 64-donor estimates but is not described as fully
corrected.

The separately validated same-atlas localization analysis pools runs to
biological donors and uses ten prespecified lineage-minus-hepatocyte contrasts.
Stromal ECM is higher in fibroblasts by β=2.271 (BH q=3.38×10^-26; HC3
q=6.88×10^-19; 53 paired donors, four datasets). Ductular injury is higher in
cholangiocytes by β=1.574 (BH q=2.04×10^-23; HC3 q=8.01×10^-21; 51 paired
donors, five datasets). Both remain the highest-scoring lineage at the 50- and
20-cell donor thresholds and retain direction in every estimable leave-one-dataset-out
fit. These results support same-atlas descriptive localization, not transcript
origin or cell autonomy.

The compact Figure 4–6 artwork candidate v8 contains eight individual one-page
PDFs and two manuscript-scale review proofs. It passed 99 checks covering exact
panel and source-table rosters, alternative-to-experiment routing, one-page
PDFs, Ghostscript readability, font subtypes, six-point plain typography,
`useDingbats=FALSE`, and compact title-free dimensions. The v8 record remains
candidate provenance. On 2026-08-17, the user approved a subsequent working
Figure 4 selection that adopts same-atlas localization as 4D and fixed-member
bulk transport as 4E; its independent render, checksum, PDF, font, and scope
checks passed in job 20872108. The user subsequently approved a shared-axis 4D
geometry that removes the redundant hepatocyte reference row without changing
the ten-contrast family or estimates; its byte, source, PDF, font, manifest, and
scope checks passed in job 20878079. The gradual 4E bulk line plot subsequently
removed its BH-supported-weight point-size encoding while preserving all eight
effects and its source table exactly; 20/20 adoption checks passed in job
20878757. Synchronized Resource-release adoption remains open.

### Superseded July COLOC cross-check and shared-variant anchors

The table below is from the **superseded July** Tier-1/2 COLOC release. It has
not been regenerated against the corrected release adopted on 2026-08-17, and
the SuSiE arm churned substantially at adoption (147 tier-1/2 genes changed
call), so no row here may be quoted without that rebuild. Coordinates are GRCh37
and denote
the top shared-variant anchor stored by COLOC; the parenthetical value is
`SNP.PP.H4`. Neither the coordinate nor `SNP.PP.H4` is a GWAS credible-set PIP
or proof that the variant regulates the named gene.

| Gene | Strongest relevant current-official result | Shared-variant anchor | Interpretation |
|---|---|---|---|
| GNMT | SuSiE PP.H4=0.1356 for UKBB ALT; single-signal ABF PP.H4=0.7076 for direct EUR NAFLD | 6:42921489 (0.2125) for SuSiE; 6:42928758 (0.3462) for direct-NAFLD ABF | No primary SuSiE call; one direct-disease ABF sensitivity call |
| MAT1A | no finite SuSiE result; ABF PP.H4=0.07305 for EUR AST | 10:81143387 (0.7790) | No positive genetic call; high variant-level sharing does not rescue low locus-level PP.H4 |
| CYP2C19 | no finite SuSiE result; ABF PP.H4=0.03427 for direct EUR NAFLD | 10:96306509 (0.05238) | No positive genetic call |
| IGFBP7 | no finite SuSiE result; ABF PP.H4=0.08484 for EUR AST; direct-MASLD ABF PP.H4=0.01333 | 4:57949517 (0.7396 for AST; 0.1245 for direct MASLD) | No positive genetic call; proxy and direct phenotypes remain separate |
| BICC1 | no finite SuSiE result; ABF PP.H4=0.02493 for EAS NAFLD | 10:60583832 (0.06040) | No positive genetic call |

Accordingly, neither recurring example is currently a primary
SuSiE-colocalized mechanism. GNMT's ABF-fallback direct-NAFLD result is retained as
a sensitivity observation, not raised to the three-gene program. All five
genes, traits, posteriors, and shared-variant anchors must be regenerated from
the adopted corrected COLOC release before final figures or Gene Catalog
entries are fixed.

An exploratory GNMT locus audit resolves the strongest current-official chain
without changing that release boundary. At the direct-NAFLD ABF anchor
rs2296805, the lower-GNMT-expression allele aligns with higher disease risk;
rs2296804 is in UKBB-European r2=1.000 with rs2296805, so the result identifies
a local haplotype rather than one nucleotide. rs9471976 and rs11752813 have
r2=0.819 and 0.580 with rs2296805, respectively, and their direct-NAFLD risk
directions also align with lower liver GNMT expression. Ji et al. independently
associated those two variants with hepatic GNMT protein in 268 human liver
biopsies and reported allele-specific reporter activity
([doi:10.1124/dmd.112.046953](https://doi.org/10.1124/dmd.112.046953)). This is
source-owned functional context, not a Resource discovery or causal-variant
assignment. The exact locus must be rechecked against the release adopted on
2026-08-17; that recheck has not been run.

## Noncoding RNA source and annotation candidate

These products are validated candidates, not an adopted release. The original
requirement that every library have a reverse/same assignment ratio of at least
four was retained as a failed sensitivity. Before opening disease outcomes, it
was replaced by a cohort-level rule requiring a median ratio of at least four,
at least 90% of libraries at or above four, and every library above one.
GSE213621 additionally required no biological group to exceed a 10% failure
rate and no more than five percentage points of failure-rate spread among
groups. All five cohorts passed that amended rule without excluding a library.

| Check | Result | Interpretation |
|---|---:|---|
| Count-source rule | 5 of 5 cohorts pass | All use featureCounts 2.1.1, GENCODE v49, and reverse-stranded counting |
| GSE126848 reverse/same ratio | median 18.52; 0 of 56 below 4 | Pass |
| GSE130970 reverse/same ratio | median 16.99; 0 of 78 below 4 | Pass |
| GSE135251 reverse/same ratio | median 15.30; 0 of 216 below 4 | Pass |
| GSE162694 ten-BAM audit | median 13.97; 0 of 10 below 4 | Pass; independently counted with `-s 0`, `-s 1`, and `-s 2` |
| GSE213621 reverse/same ratio | median 10.98; 346 of 367 at or above 4; all above 1 | Passes the amended cohort rule |

The 21 lower-ratio GSE213621 libraries are distributed across control, F0/F1,
F2, and F3/F4 groups rather than confined to one stage. Group-specific failure
counts are 5/68, 5/97, 7/107, and 4/95, respectively; the failure-rate spread is
3.14 percentage points. Every library still has more reverse- than same-strand
assignment. Both the failed strict rule and the passed amended rule are stored
with their hashes and rationales.

The independent GENCODE and prespecified-program annotation lanes are valid:

| Candidate quantity | Value | Qualification |
|---|---:|---|
| GENCODE v49 gene records | 86,369 | Versioned Ensembl IDs are primary |
| GENCODE v49 lncRNA genes | 36,490 | Includes alternative loci for coverage audit |
| Canonical-chromosome lncRNA genes | 34,880 | Chromosomes 1–22, X, and Y |
| Main-text-eligible lncRNA mappings | 33,991 | Canonical chromosome and neither same-strand exonic nor complex |
| Prespecified programs annotated | 117 | Memberships, hashes, and weights unchanged |
| Uniquely mapped lncRNA memberships | 409 | Membership rows, not distinct genes |
| Programs requiring leave-all-lncRNA-out sensitivity | 49 | At least 5% of original L1 weight is uniquely mapped lncRNA |

The seven lncRNA genomic classes contain 20,085 intergenic, 7,361 antisense
exonic-overlap, 4,591 antisense gene-body-overlap, 2,037 same-strand gene-body,
1,450 opposite-strand promoter-proximal, 632 same-strand exonic-overlap, and
334 complex genes. Classification is genomic context, not evidence of lncRNA
function.

### Five-cohort pooled biotype subdivision

The primary inference remains one all-gene family under the canonical threshold. The
biotype rows are subdivisions of its 23,370 tested genes and 1,347 positives,
not a separate lncRNA-only multiple-testing analysis. One effect-size floor
applies to every biotype; see the note below on why a biotype-specific floor was
rejected.

| Biotype | Tested | Canonical positive | Up | Down | TREAT comparator | Up | Down |
|---|---:|---:|---:|---:|---:|---:|---:|
| Protein coding | 14,664 | 886 | 744 | 142 | 1,137 | 918 | 219 |
| lncRNA | 6,249 | 402 | 220 | 182 | 434 | 205 | 229 |
| Other | 2,457 | 59 | 39 | 20 | 45 | 21 | 24 |

**The unadjusted lncRNA share must not be read as enrichment.** lncRNAs are
26.7% of tested genes and 29.8% of canonical positives, against 26.9% under the
TREAT comparator. That apparent rise is a property of the threshold, not of biology:
lncRNAs are lower-expressed, so their fold-change estimates are noisier, and a
floor applied to the point estimate preferentially includes them. The interval-null
comparator, which weights by standard error, holds the share at parity (ratio
1.004). Lowering the floor makes the artifact worse, not better — the lncRNA
share of positives rises to 31.7% at a 0.30 floor and 32.1% at 0.20.

A biotype-specific effect-size floor was considered and **rejected**. It would
amplify the same artifact, and a class-specific filter makes the biotype
comparison circular by construction: the class given the looser filter is
enriched automatically.

The adjusted analysis is the interpretable one, and it is unchanged in direction
and strengthened in magnitude. In a secondary descriptive gene-level logistic
model adjusted for average abundance, expression variability, gene length, and
cohort detection, lncRNA annotation was associated with **lower** odds of a
positive call than protein-coding annotation: odds ratio **0.439** (95% CI
0.381–0.505; nominal P=1.13×10⁻³⁰; 20,913 genes) under the canonical threshold,
against 0.488 (95% CI 0.427–0.557; nominal P=4.38×10⁻²⁶) under the TREAT
comparator. Once abundance and variability are accounted for, the apparent
unadjusted enrichment reverses. This describes the tested gene universe; genes
are not biological replicates, and the result is not a donor-level
disease-effect test.

Regenerated under the canonical threshold on 2026-08-12
(`bulk_disease/run_canonical_gate_2026-08-12/`), the high-confidence lncRNA set
is **279**, against 271 in the retained TREAT arm. The count rises even though
the lncRNA positives fall from 434 to 402, so the high-confidence rate goes from
62.4% to 69.4% of positives. This is a property of the threshold: a hard floor on the
point estimate selects large-effect genes, and large effects re-clear the same
floor in the leave-one-cohort-out and equal-weight refits more reliably than the
modest-but-precise genes the interval null included. Both arms are retained
side by side; neither directory overwrote the other.

Among the 434 TREAT-positive lncRNAs, **271** meet every prespecified
high-confidence rule: eligible mapping, same direction in at least four of five
cohorts without a materially opposite fifth cohort, retained direction in all
five leave-one-cohort-out fits, retained TREAT significance in at least four,
and agreement under equal library weighting. The deterministic main-text
example rule selects `ENSG00000213062.6`, an unambiguous opposite-strand
promoter-proximal lncRNA record with no separate GENCODE display symbol
(log2 fold change -0.761; pooled t=-16.55). It must be displayed by stable ID;
no familiar name may be substituted post hoc.

Three one-page supplementary candidates report the complete robustness family,
genomic-class counts, and the deterministic stable-ID example. Each PDF passed
source-table checksum, one-page, Ghostscript, and visual review. They remain
candidate supplement panels and do not add a sixth main figure.

### DEG-threshold sensitivity

The canonical threshold and its alternatives are all reported. Values below are on
the current official substrate and are reproduced by
`figures/supplementary/figS_deg_threshold_sensitivity/`.

⚠ **Every number in this section is on the superseded July COLOC release.** The
whole analysis was rerun on the adopted 2026-08-17 release
(`RNA-seq/results/genetics_state_interface/20260817T133951Z_promoted_coloc/`).
The conclusions survive; the values move:

| | July (below) | adopted 2026-08-17 |
|---|---:|---:|
| jointly testable background | 14,931 | **14,920** |
| colocalized genes in it | 447 | **437** |
| observed / expected overlap, canonical threshold | 23 / 32 | **20 / 31.2** |
| overlap odds ratio, canonical threshold | 0.690 (p = 0.094) | **0.616 (p = 0.037)** |
| matched ladder at the 0.50 floor | 23 against 24 controls | **20 against 16**, McNemar OR 1.31, p = 0.58 |
| threshold-free SuSiE rho | 0.019 (p = 0.232) | **0.016 (p = 0.34)** |

The ladder is still flat at every floor and the threshold-free estimate is still
null, so **the "do not call colocalized genes buffered or constrained" rule
stands**, and the slightly stronger nominal depletion does not reverse it. Quote
the adopted values; the July numbers below are kept because the corrections
they document are why the analysis is trusted at all.

All overlap tests below run against the **14,931-gene jointly testable
background** — genes assayed in both the bulk contrast and colocalization — and
the 447 colocalized genes within it. DEG counts are given both genome-wide and
restricted to that background; only the restricted count enters the test.

| Threshold | DEGs (genome-wide) | DEGs (in background) | Non-overlap with the genetic set | Overlap odds ratio | Fisher p |
|---|---:|---:|---:|---:|---:|
| `padj < 0.05` and `|log2FC| > 0.50` (**canonical**) | 1,853 | 1,079 | 94.85% | 0.690 | 0.094 |
| Interval-null FDR < 0.05 at `lfc = 0.25` (calibrated comparator) | 1,918 | 1,261 | 92.39% | 0.889 | 0.604 |
| `padj < 0.05`, no fold-change floor | 13,043 | 8,435 | 41.16% | 1.104 | 0.333 |

The canonical and comparator thresholds agree on 1,497 genes; 356 are
canonical-only and 421 are TREAT-only (Jaccard 0.658). The two sets are
therefore not interchangeable even though their totals are close, and any
figure must state which one it used.

**Corrected 2026-08-13: the earlier ladder used the wrong background.**
`scripts/figures/figS_deg_threshold_sensitivity_data.py` set the Fisher
background to all 27,638 bulk-tested genes rather than the 14,931 jointly
testable ones. A gene never tested for colocalization cannot count as evidence
against overlap; including such genes inflates the double-negative cell and
therefore every odds ratio. This is the asymmetric-denominator trap. The
producer now asserts against `orthogonality_audit.tsv`, which was always
correct, so the two can no longer disagree.

Two conclusions from the retired version do not survive:

- The claimed significant *enrichment* at permissive thresholds was an artifact.
  No floor: 1.612, p = 7.8e-07 → **1.104, p = 0.333**. Floor 0.10: 1.296,
  p = 0.0075 → **1.082, p = 0.412**. Neither is significant.
- The comparator threshold no longer points the other way. Its odds ratio was
  reported as 1.106 (enrichment) against the audit's 0.889 (depletion) for the
  same genes; **0.889 is correct**. The threshold choice never reversed the
  direction — the denominator did.

**The 0.50 floor is no longer load-bearing.** On the corrected background the
odds ratio decreases monotonically across all eleven floors, from 1.104 with no
floor to 0.395 at `|log2FC| > 1.00`, and is never significantly above 1 at any
floor. Depletion reaches nominal significance only at floors 0.25 and 0.30
(0.778 and 0.736; p = 0.046 and 0.033), neither surviving correction across the
eleven tests. The direction of the result is therefore stable at every
threshold, and the earlier instruction to always state the conditioning floor
is superseded — though the magnitude still depends on it and should be quoted
with the floor attached.

Across conventional floors from 0.30 to 1.00 the non-overlap fraction runs
87.0% to 99.6%. The corrected reading is that colocalized genes are
progressively **less** represented among genes with larger disease-state
responses, with no significant over-representation among weakly responding
genes at any floor.

**The ladder is an expression artifact, not a biological property.** Tested
directly: `scripts/manuscript/response_amplitude_matched.py`, results under
`RNA-seq/results/response_amplitude/`.

Colocalized genes are markedly higher-expressed than the rest of the background
(mean `AveExpr` 3.696 against 2.590, SMD 0.363), and a high `|log2FC|` floor
preferentially includes low-expression genes whose effect sizes are inflated
conditional on passing FDR. Matching each of the 447 colocalized genes to an
expression-matched control (L2 propensity, 0.2 SD logit caliper, 447 of 447
matched, max |SMD| 0.00002) removes the entire effect:

| Test | Result |
|---|---|
| Mean amplitude, `\|shrunk_logFC\|` | Δ = **+0.003**, 95% CI [−0.019, +0.025], Wilcoxon p = 0.64 |
| Mean amplitude, raw `\|logFC\|` | Δ = +0.007, 95% CI [−0.016, +0.030], p = 0.38 |
| Adding `SE` to the match (conservative) | Δ = −0.009, 95% CI [−0.036, +0.016], p = 0.71 |
| Decile-matched null, 10,000 draws | observed 0.1416 against null 0.1369 ± 0.0079, p = 0.556 |

The mean is not the estimand the ladder makes, so the ladder itself was
recomputed **within matched pairs**, scoring each colocalized gene and its
matched control DE or not at each floor and testing the discordant pairs with
exact McNemar. The matched ladder is flat: odds ratios 1.01, 1.17, 1.03, 1.11,
1.11, 1.07, 1.03, 0.95, 0.87, 1.14, 0.67 across the eleven floors, every p ≥
0.27, against the unmatched decline from 1.104 to 0.395. At the canonical 0.50
floor, 23 colocalized genes are strongly differentially expressed against
**24** of their expression-matched controls.

Conclusion: colocalized genes respond to disease exactly as much as comparable
genes. **Do not describe them as buffered, constrained, or transcriptionally
held constant**, and do not present the ladder's decline as a property of
inherited risk. The non-overlap is what independent sampling of the
transcriptome produces for two sets of this size and expression profile —
447 genetic and 1,079 disease-state genes in a 14,931-gene background give an
expectation of 32 shared genes against 23 observed — not evidence that the two
layers occupy structurally different parts of the transcriptome.

The adopted-release threshold-free analysis is now complete on the
synchronized 23,370-gene fragment-count arm. Among 14,112 intersecting genes,
3,617 have finite SuSiE PP.H4 and raw RNA-seq log2FC estimates; the signed
association is Spearman rho=0.0141 (P=0.396). Figure 3G shows these observations
without a fitted line.

Figure S3G contains the expression-matched sensitivity using the calibrated TREAT family (lfc=0.25, BH
q<0.05), not a post-FDR raw fold-change cutoff. All 428 Genetics sig. genes
were paired one-to-one with expression-matched controls. Thirty pairs contain a
Transcriptomics sig. call only in the Genetics sig. member, 21 only in the
matched control, 377 in neither, and none in both (paired OR=1.43; exact
McNemar P=0.262). The result supports complementary evidence maps without
evidence of biological antagonism; because PP.H4 is unsigned, it does not test
risk-allele direction.

**Caveat that belongs in Methods.** Because the effect-size filter is applied
after FDR adjustment rather than folded into the test, the adjusted p-value
does not carry a calibrated error rate for the filtered set (McCarthy & Smyth
2009). The TREAT arm is retained as the calibrated comparator for exactly this
reason and its columns must not be deleted from `canonical_deg_results.csv`.

## Cross-sectional stage remodeling candidate

Clean job `19966029` completed from the validated 844-participant, five-cohort
fragment-count substrate. All four adjacent fibrosis contrasts, three adjacent
NAS-group contrasts, and F1–F4-versus-F0 models use fixed cohort and inferred-sex
terms and complete 23,370-gene Benjamini–Hochberg families. Input/output hashes,
family cardinalities, design rank, and reported FDR values passed independent
checks. Candidate v2 is validated but not adopted; job `19957749` and candidate
v1 remain failed execution provenance.

| Axis | Contrast | Reference / comparison participants | Cohorts | Higher / lower at BH q<0.05 |
|---|---|---:|---:|---:|
| Adjacent fibrosis | F0→F1 | 104 / 105 | 3 | 682 / 789 |
| Adjacent fibrosis | F1→F2 | 105 / 90 | 3 | 354 / 194 |
| Adjacent fibrosis | F2→F3 | 90 / 75 | 3 | 376 / 146 |
| Adjacent fibrosis | F3→F4 | 75 / 25 | 3 | 174 / 70 |
| Fibrosis versus F0 | F1 vs F0 | 104 / 105 | 3 | 682 / 789 |
| Fibrosis versus F0 | F2 vs F0 | 104 / 90 | 3 | 2,011 / 1,893 |
| Fibrosis versus F0 | F3 vs F0 | 104 / 75 | 3 | 2,986 / 3,334 |
| Fibrosis versus F0 | F4 vs F0 | 104 / 25 | 3 | 2,999 / 3,425 |
| Adjacent NAS | 0→1–2 | 46 / 65 | 3 | 2,692 / 2,435 |
| Adjacent NAS | 1–2→3–4 | 65 / 121 | 3 | 2 / 1 |
| Adjacent NAS | 3–4→5–8 | 121 / 175 | 3 | 1,877 / 2,213 |

Higher and lower are defined in the comparison group relative to the named
reference. These counts use BH q<0.05 without an added fold-change threshold.
The continuous NMF display starts at NAS1–2 because only one NAS0 donor has a
prespecified NMF score; the NAS0-versus-NAS1–2 DEG model still uses all 111
eligible participants.

The prior 661-sample, six-cohort stage candidate and every earlier read-count
stage table are provenance only. They cannot supply Figure 3 counts or be
combined with the selected five-cohort pooled result. Every stage contrast is
cross-sectional and uses participants from cohorts carrying the compared
groups; no comparison represents within-person change.

## Histology-anchored molecular continuum candidate

The validated candidate
`RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-20260818T173348Z`
adds an externally anchored molecular continuum without changing Figure 3 or
replacing recorded histology. GSE162694 and GSE213621 are modeled separately as
competitor-independent evaluation cohorts. Within each cohort, the primary
model is `outcome_z ~ continuum_z + factor(fibrosis_stage) + inferred_sex`;
fixed-effect meta-analysis then combines cohort coefficients. The cohort-PC1
and fixed external projection are co-primary, and the consensus rank is
visualization-only. All 139 observed genes from the published 145-gene
signature are removed before pathway,
Hotspot, and fixed-basis NMF scoring; ordinary transcript inference excludes
the 139 observed axis components, which remain explicit and are evaluated only
with target-specific leave-one-gene-out axes.

| Layer | Complete family / testable | Continuum-associated | Matched fibrosis-trend comparison |
|---|---:|---:|---:|
| Ordinary transcripts | 23,370 / 23,231 | 14,651 | 10,476 stage-supported; 6,836 shared; effect rho=0.314 |
| Hotspot programs | 117 / 116 | 93 | 64 stage-supported; 54 shared; effect rho=0.803 |
| Hallmark | 50 / 44 | 37 | 29 stage-supported; 27 shared; effect rho=0.693 |
| KEGG | 658 / 260 | 170 | 123 stage-supported; 100 shared; effect rho=0.622 |
| Reactome | 1,787 / 1,027 | 781 | 607 stage-supported; 525 shared; effect rho=0.643 |
| GO biological process | 7,583 / 2,811 | 2,108 | 1,786 stage-supported; 1,531 shared; effect rho=0.696 |
| GO molecular function | 1,855 / 720 | 462 | 342 stage-supported; 259 shared; effect rho=0.642 |
| GO cellular component | 1,042 / 467 | 332 | 244 stage-supported; 203 shared; effect rho=0.648 |
| Fixed-basis NMF axes | 10 / 10 | 8 | 7 stage-supported; 6 shared; effect rho=0.733 |

The two transcript-derived continuum effect vectors agree at rho=0.9991. Both
co-primary scores support positive stage-adjusted associations for the
prespecified Ductular injury and Stromal ECM programs in both evaluation
cohorts, and the additive Figure 4F adoption criteria pass. The complete
stage-versus-continuum review candidate is
`figures/candidates/histology-vs-continuum-comparison-20260818T192630Z`.
All five Resource cohorts were scored. The current descriptive continuum-window
figures show GSE126848 (n=53), GSE130970 (n=76), GSE135251 (n=214), GSE162694
(n=140), and GSE213621 (n=361). Stage displays contain the latter four
stage-available cohorts; GSE126848 has no recorded fibrosis in the selected
manifest. The first three cohorts remain source-overlap descriptive context and
do not enter the primary two-cohort meta-analysis.
Program and pathway rankings align strongly with fibrosis-based analyses, while
the weaker transcript-level concordance and significant stage-adjusted effects
show that the continuum is not merely a relabeling of stage groups.

All sliding windows are descriptive, fixed in advance, and excluded from
inference. GSE193066 paired evidence upgrades supported layers but does not
veto cross-sectional findings when underpowered: Ductular injury is
within-person supported, whereas Stromal ECM is untestable in the paired
substrate because retained signature-excluded weight is below 80%. The safe
claim is that the continuum complements histology by resolving coordinated
cross-sectional molecular variation within recorded stage. It does not support
a biomarker, subtype, causal, therapeutic-window, or within-person
temporal-ordering claim. The synchronized Resource candidate remains unpromoted.
On 2026-08-24, hash-pinned Figure 4F–4G and selected S3/S4 panels were adopted to
the working figure authorities, then revised from two-cohort-only descriptive
views to the five-cohort context package
`histology-continuum-five-cohort-context-20260824T175647Z`. This working-figure
revision does not adopt the synchronized Resource release.

The molecular-systems layer exists to give the continuum a unit of
interpretation between single genes and whole-cohort axes. It compresses 13,092
correlated features into a few dozen auditable groupings so that within-stage
remodeling can be described as coordinated rather than gene by gene. That
purpose is organizational rather than inferential, and the stress tests below
bound it: the groupings are a navigable vocabulary and evidence that within-stage
variation is shared across pathway families and single-cell-derived programs,
not a ranking of systems by disease relevance.

The corrected outcome-blind molecular-systems atlas is
`figures/candidates/histology-continuum-molecular-systems-atlas-20260824T195148Z`.
After exclusion of the 139 observed genes from the published 145-gene
signature and GENCODE-v49 mapping, it retains 5,329 testable pathway sets and
116 testable Hotspot programs from the complete 13,092-feature ledger. Of
5,445 mapped nodes, 5,439 form one main graph component and 43 size-ranked
systems. Six nodes form three disconnected two-node microcomponents; they
remain explicit in the registry and topology view but are not named or included
in system summaries. Stage and stage-adjusted fixed-projection effects were
added only after the coordinates, component registry, and community IDs were
checksummed.

Across the 5,439 main-component nodes, stage-versus-continuum effect rho is
0.683 and 80.5% share direction. Equal-collection system summaries have
rho=0.658 and 36/43 share direction; feature-balanced medians have rho=0.701
and 37/43 share direction. Fixed-window direction agrees between GSE162694 and
GSE213621 for 39/43 systems under either aggregation. Across all five
descriptive cohorts, the counts are 13/43 and 12/43, respectively. These values
are descriptive summaries over nested gene sets, correlated cohorts, and
overlapping windows; no independent-sign null or inferential p value is used.
Outcome-blind topology sensitivity gives adjusted Rand indices of 0.633 at
k=20, 0.751 at k=40, 0.855 at resolution 0.7, and 0.815 at resolution 1.3
relative to k=30/resolution 1. Biological terms remain in the top-ten
representative audit rather than becoming single-medoid system names. Figure
4G, Figure 4H, and S4Q–S4U are adopted working panels; Figure 3 and the synchronized
Resource release remain unchanged.

Figure S4U pools the five descriptive system-window profiles only after
within-cohort standardization and percentile ranking. For each system and
window, collection medians are averaged within cohort and the five resulting
cohort scores receive equal one-fifth weight. The source lineage retains the
complete 43-by-5-by-9 cohort input and 43-by-9 pooled output. This is a
descriptive aligned-cohort summary, not participant-level pooling or an
additional significance test.

The separate donor-level system analysis constructs each system score within
participant, then fits `system_score_z ~ continuum_z + factor(fibrosis_stage) +
inferred_sex` in GSE162694 and GSE213621. Every stage-adjusted model is
complete-case in recorded fibrosis stage. Thirty-one of the 140 GSE162694
participants have no recorded stage in the selected manifest, so the
stage-adjusted models use 109 GSE162694 and 361 GSE213621 participants, 470 in
total, not the 501-participant descriptive roster. Fixed-effect
meta-analysis uses BH over the complete 43-system family for each co-primary
axis. Thirty-three systems meet both meta-analysis q<0.05 criteria with the same
direction in all four cohort-by-axis coefficients; 29 retain support with
feature-balanced scoring. Twenty-seven systems are fibrosis-stage-associated,
24 pass both stage and continuum criteria, and nine are continuum-only.
Within-stage/sex permutations give BH support in all four tests for 22 systems;
18 pass the stricter maxT family-wise P<0.05 criterion in all four tests.

Two boundaries on those counts are recorded here rather than left to the reader.
First, the two evaluation cohorts do not carry the same stage scale: GSE213621
records levels 0 through 3 and GSE162694 records levels 0 through 4, yet the
secondary `stage_z` model standardizes each cohort's own scale before
fixed-effect meta-analysis, so the combined stage coefficient is not expressed in
a single common histologic unit. Second, none of the nine continuum-only systems
pass the strict criterion; all 18 maxT-supported systems are
stage-and-continuum-supported. The strict criterion is the conjunction of four
independently maxT-corrected axis-by-cohort tests, and GSE162694 contributes 109
participants across ten stage-by-sex strata, so the conjunction is limited by the
smaller cohort's within-stratum permutation resolution. All nine continuum-only
systems agree in direction across both cohorts; their GSE162694 residual
correlations are smaller than their GSE213621 values. The strict arm therefore
does not by itself establish continuum association independent of recorded stage
at the system level. The
donor-first and plotted pooled-window estimands agree at Spearman rho=0.997.
Main Figure 4H shows the complete 43-system fibrosis-stage versus
stage-adjusted continuum comparison. Labels remain restricted to the 15 largest
systems fixed before outcome overlay. It uses neither sliding-window values nor
membership geometry. The displayed rho=0.658 and 36/43 same-direction count
remain descriptive summaries of correlated systems.

The revised panel is
`figures/candidates/histology-continuum-fig4h-competitive-mark-20260827T130920Z`.
Every point, coordinate, label, and summary value is carried over unchanged; only
what the ink encodes has moved. The outline previously marked the 18 systems
clearing the four-test maxT conjunction. That criterion is set by the smaller
cohort's within-stratum permutation resolution, and the jointly powered
cross-cohort test supports 35 of 43, so neither count is scarce enough to carry
ink. The outline now marks the single system exceeding size-matched draws on both
co-primary axes. Both counts are retained in the panel source table and manifest
for the caption authority. The rebuild also repairs a defect carried by the
adopted panel: the direction-colour legend keys rendered black because the
legend-only layer contributed its fixed colour to the colour guide.

Supplementary panels S4V and S4W are
`figures/candidates/histology-continuum-system-competitive-null-20260827T131027Z`.
S4V shows the observed system_01 statistic against its own null distributions
under both references and both axes, which makes the separation between them
visible: the connected-subgraph null is broad and reaches the observed value,
while the unconstrained null is tight. S4W retains all 43 systems and plots the
deciding quantity, the worse BH q across the two co-primary axes, under each null.
Figure S4U retains all 43 systems and marks only the 18 strict maxT-supported
rows with a neutral `*` beside unchanged window tiles. The 33 primary BH calls
remain in the complete source registry; the panel does not assign P values to
windows or `Delta`. The active Figure 4G/S4Q–S4U display revision removes
embedded explanatory footnotes without changing values, geometry, order, or
inference; the full explanations remain in the external caption authority.
Figure 4H is hash-pinned in
`fig4h_system_stage_continuum_promotion_20260825.tsv`.

The active visual revision orders both system heatmaps by the equal-cohort
pooled late-minus-early contrast, defined as the mean of windows centered at 80
and 90 minus the mean of windows centered at 10 and 20. The active S4U omits
this redundant contrast tile but retains the same row order. This is a
result-adaptive presentation order only; all underlying values, geometry,
support calls, and testing families are unchanged.

### Molecular-system stress tests

The stress-test candidate is
`RNA-seq/results/histology_anchored_continuum/molecular_systems_stress/hac-system-stress-20260827T125516Z`.
It reads only the checksummed atlas and inference candidates, verifies both
against their checksum registries, and reproduces the released collection-balanced
participant scores to within 5.3e-15 before any test runs. It changes no
geometry, community ID, Hotspot membership, or adopted panel.

A jointly powered cross-cohort permutation replaces the conjunction of four
separately corrected axis-by-cohort tests. The axis is permuted within
cohort-by-stage-by-sex strata, the two cohorts are combined by weighted Fisher z
of the residual partial correlation, and maxT runs over the complete 43-system
family. Thirty-five of 43 systems pass family-wise P<0.05 on both co-primary
axes, against 18 under the four-test conjunction, and all nine continuum-only
systems are among them. The earlier absence of continuum-only systems from the
strict arm was therefore a resolution limit of the smaller cohort's
within-stratum permutation, not evidence against within-stage association. Both
criteria are retained; the prespecified 18-system count is unchanged.

The competitive nulls ask a different question: whether a system is more
continuum-associated than an equally sized region of the same geometry rather
than more than zero. Each system is compared against 5,000 size-matched draws,
scored and fitted through the identical pipeline, under a connected snowball
sample of the prespecified mutual-30NN graph and under an unconstrained node sample.
The per-draw null statistics are retained so the supplementary panels can show
the distributions themselves. Under the conservative connected null, **no system
passes BH q<0.05 on both axes**; under the permissive unconstrained null, one
does. That system is
system_01, the largest system at 311 features, whose representative terms are
integrin binding, ECM organization, contractile myofibroblast, and cell
adhesion. Its connected-null P values are 0.0062 and 0.0078, which do not
survive BH over 43 systems, and the next-ranked system is at P=0.115. Drawn sets
are not hub-enriched relative to the tested systems: mean node degree is 22.1
under the connected null and 21.3 under the unconstrained null against 19.7 for
system_01. The 33-of-43 and 35-of-43 counts therefore describe a diffuse axis
shared across the membership map, consistent with 14,651 of 23,231 transcripts
being continuum-associated after stage adjustment. Except for the ECM system
against the permissive reference, no system is shown to carry
continuum association specific to its own membership.

The 43 systems are also not 43 independent quantities. Across the pooled
evaluation participants their score correlation matrix gives a Li and Ji
effective number of 15 independent systems, with 36.5% of variance in the first
component and 72.4% in the first three; the per-cohort effective numbers are 12
and 14. BH over 43 controls multiplicity but not this redundancy.

Finally, the apparent strength of the continuum relative to fibrosis stage is
largely a comparison artifact. The 2.79-fold median ratio visible in the Figure
4H source table compares a marginal stage coefficient against a
stage-adjusted continuum coefficient. Fitting both marginally in the same
participants gives a median ratio of 1.57, and binning the continuum to each
cohort's own stage coarseness gives 1.45. The reliability correction is no longer assumed. Kleiner
fibrosis staging was measured against itself in GSE193066, in 58 participants
with repeat biopsies a median two years apart: exact agreement 28 of 58, within
one stage 52 of 58, and a within-participant standard deviation of 0.6433
stages with a participant-bootstrap 95% interval of 0.5170 to 0.7600. The
transportable quantity is that standard deviation rather than a reliability
coefficient, because the same error variance yields reliabilities of 0.19, 0.28
or 0.47 in this deposit depending only on which stage variance is used as the
denominator. Applied as rho = 1 - 0.414 / Var(stage) in the regression-slope
convention used above, the reliability-corrected ratio spans 0.27 at the paired
subset's own restricted variance to 1.03 at a stage standard deviation of 1.2,
so a coarse, imperfectly reproducible ordinal and a continuous expression
projection remain unseparated by these data, and the correction is now anchored
to a measurement rather than to a sensitivity range.

That measurement does not transport unchanged. The same statistic was computed
in an independent deposit, PXD051911, in fifty-seven patients with repeat
Kleiner grades: 0.5092 stages with an interval of 0.3043 to 0.6939 in the
twenty-seven pairs taken at an unchanged age, and 0.4655 stages with an interval
of 0.3651 to 0.5477 in the thirty non-surgical follow-up pairs, the latter
formally discrepant with GSE193066 because the ratio interval of 0.5366 to
0.9555 excludes one. Neither stage composition nor real interval change accounts
for the gap: standardising to the source stage mix moves the value further away
rather than closer, and the follow-up arm shows a significant net upward
shift in stage, which can only inflate a within-pair variance. Two deposits cannot support a
pooled estimate, so the transportable quantity is the range 0.4655 to 0.6433,
and any user with repeat reads should measure it locally. The direction of the
resulting bias is counterintuitive and matters: a larger standard deviation
implies a lower reliability, hence a lower ceiling and a larger disattenuated
correlation, so the GSE193066 value is the flattering end for any claim that a
model approaches the measurement ceiling, and 0.4655 is the conservative
choice.

This is a test-retest statistic across two physically distinct biopsies. It
folds reader variability, biopsy sampling variability and any real change into
one term, and no inter-reader substrate exists in any deposit in this project,
so no claim about agreement between readers is available. The disagreement is
symmetric, 15 up and 15 down with a sign test p of 1.0 and a mean delta of
-0.034 stages whose interval covers zero, which excludes a net upward drift in
stage as the explanation but does not exclude balanced bidirectional change;
the arm that could have caught the latter is null and underpowered. The intraclass
correlation in this sample is 0.1886 with a bootstrap 95% interval of -0.0918
to 0.4343 and does not clear a shuffled-pairing null, so it is reported with its
interval and never as a point reliability. The same coefficient computed on the
same label in PXD051911 is 0.62, so no statement of the form that Kleiner
staging has a reliability of 0.19 is available. The gap is not a denominator
artefact: the implied stage variances are 0.510 and 0.558, within nine per cent
of each other, while the error variances are 0.414 and 0.217, a factor of 1.9.
Because a reliability is one minus a ratio of those two quantities, a 1.9-fold
difference in error variance produces a 3.2-fold swing in the coefficient at an
essentially unchanged denominator. Across the two deposits the coefficient moves
3.2-fold while the standard deviation moves 1.38-fold, and it is that stability
under the same underlying difference, not any claim that the error itself is
constant, which is the reason the standard deviation is the quantity reported
here. Because a repeat-biopsy interval
loads real change into the within-participant term, the measured value is a
lower bound on single-occasion reliability. The 58 paired participants are all
high-risk at first biopsy and their stage variance is 0.74 times that of the
full 106, so they are a selected set.

## Single-cell programs and composition

| Quantity | Value | State |
|---|---:|---|
| Prespecified Hotspot programs | 117 | Prespecified registry |
| Cell types represented | 5 | Hepatocyte, macrophage, fibroblast, cholangiocyte, T cell |
| Programs meeting primary association/stability criteria | 2 | Source-diagnosis-refit candidate on the original ambient-uncorrected atlas; BH over all 117 programs |
| Programs passing HC3 family correction | 0 | Source-diagnosis-refit candidate; HC3 BH q=0.105 for both primary-selected programs |
| Integrated libraries / biological donors / cells | 275 / 104 / 1,232,318 | Descriptive integration census; Liver Atlas contributes 42 runs, 38 assays, and 19 biological donors |
| Analyzed libraries / biological donors / cells | 273 / 102 / 1,232,285 | Two GSE189600 libraries with 13 and 20 annotated cells are excluded; donor collapse is mandatory for inference |

The corrected donor census replaces the stale 125-donor count, which treated
Liver Atlas sequencing runs as donors. The disease-state refit contains 64
biological donors across GSE185477, GSE202379, GSE244832, GSE174748, and the
two analyzed GSE189600 donors, with no Liver Atlas rows. GSE202379 is mapped
from its source diagnoses rather than fibrosis: Healthy control to Healthy,
NAFLD to Steatosis, and NASH without cirrhosis to Steatohepatitis. All 117
program memberships and weights are unchanged and all 117 effects are
estimable. The prior 113/114 tested, 27/26 supported, and four-selected/two-HC3
values belong to the superseded metadata analysis and are not reader-facing
results.

| Atlas label | Program | Primary q | HC3 q | Source-diagnosis-refit state |
|---|---|---:|---:|---|
| Hepatocyte | Stromal ECM | 0.0496 | 0.105 | Primary-selected; HC3-fragile |
| Hepatocyte | Ductular injury | 0.0187 | 0.105 | Primary-selected; HC3-fragile |

A supplementary contextual screen crossed these two starred programs with the
fixed 13-pair communication roster, giving 26 planned donor-level association
tests. Twenty-four were estimable; both *LAMC3*→*DAG1* tests had fewer than
30 matched donors and remained untestable. Models adjusted for disease stage,
dataset, and log sender and receiver cell counts. One association passed BH
q<0.05 in both the model-based and HC3 26-test families: the Ductular injury
program with hepatocyte→endothelial *CDH1*→*PTPRM* (60
biological donors across four datasets; β=0.501, model-based BH
q=1.30×10−6, HC3 BH q=5.78×10−13). Its coefficient remained positive in
all four leave-one-dataset-out fits (β range 0.090–0.676). No Stromal ECM
program association passed both families, and the original donor-stage communication
family remained zero of 13 FDR-supported. Both scores were derived from the
same atlas, so this is contextual covariation rather than evidence of ligand
secretion, receptor activation, or mechanism. The complete display is routed
to Figure S4 rather than the main figure.

Working Figure 4E transports 46 mapped Stromal ECM and 21 Ductular injury
member genes into bulk fibrosis-stage contrasts. At F4 versus F0, 36 and 18
genes meet BH q<0.05, representing 79.3% and 86.4% of mapped program weight;
weighted directional agreement is 98.2% and 100%, respectively. These are
member-level summaries, not program-level significance tests. The TF analysis
remains supplementary: 268 TFs were tested within each F-versus-F0 contrast,
with 35, 50, 44, and 85 meeting BH q<0.05 in F1, F2, F3, and F4, respectively.
The source table retains all 268 TFs and complete contrast-specific BH families.

The corrected candidate is recorded at
`Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/`.
Do not reconstruct its calls from older donor-collapse or four-program tables.
It predates the ambient-RNA sensitivity summarized above and therefore cannot
by itself establish hepatocyte origin.

The 117 memberships and weights are fixed; no external outcome may trigger
module rediscovery, reweighting, or renaming. Final cell-composition numbers
must be recomputed against the 1,257-sample corrected fragment substrate.

Forty-nine programs carry at least 5% lncRNA L1 weight and therefore required a
leave-all-lncRNA-out sensitivity in the original single-cell scoring framework.
The exact stored Hotspot score was retained, and only the gene-removal delta
from the source-reconstructed graph was added. Reconstructed all-gene scores
correlated at least 0.99896 with stored scores. This 0.995 technical floor was
finalized during source reconstruction and is not the biological acceptance
threshold. At the biological-donor level, **48 of 49** programs retained both
Pearson and Spearman correlation of at least 0.90 and preserved the stage-effect
direction. The minimum donor Pearson and Spearman correlations were 0.9688 and
0.9442.

The sole failure was cholangiocyte module 28: removing four lncRNA members
(11.8% of L1 weight) changed its small stage coefficient from +0.0164 to
-0.00325. It remains in the complete registry but cannot be described as
lncRNA-independent. Both source-diagnosis-refit primary-selected programs passed this
lncRNA-removal sensitivity: Stromal ECM
had donor Pearson 0.9987 and Ductular injury had donor Pearson 0.9979,
with direction preserved for both. This sensitivity tests dependence of the
program score on its lncRNA members; it does not establish that any lncRNA is
functional.

### Validated noncoding plot package

The read-only candidate
`figures/candidates/noncoding-resource-2026-08-12/plot-upgrade-v3/` contains
nine individual PDFs and nine source tables. Independent checksum, page,
Ghostscript, source-cardinality, and scientific-census validation passed, as
did manual visual review. Repository-wide tests and the Resource scope
validator also passed in SLURM job `19779385`. The earlier `v1` and `v2` candidates are explicitly
rejected for source-table serialization and visual-label defects, respectively.

The added composition panel makes the correct null visible: lncRNAs constitute
26.74% of all 23,370 tested genes and 26.86% of the 1,616 TREAT-positive genes.
This is comprehensive inclusion, not preferential lncRNA enrichment. The
program panels show all 117 original lncRNA-weight fractions and the complete
49-program sensitivity family; 48 preserve the donor-level stage direction.
The Figure 6 panel visualizes five deterministic molecular-object routes and
contains no corrected genetics outcome or Cas13 screen result.

### Validated noncoding-DNA candidate and plots

The read-only DNA-only candidate
`GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/`
uses the read-only fine-mapping substrate and does not read partial corrected
COLOC outcomes. Of 6,467 reliable Tier-1/2 study-by-locus credible-set
instances, 6,415 (99.20%) have raw SuSiE PIP sums within 0.90–1.10. Their
42,480 study–credible-set member rows collapse to 18,539 unique variants.

The complete consequence table was generated through the official Ensembl
GRCh37 REST VEP endpoint after orienting every stored allele pair against the
GRCh37 reference sequence. The first pass that assumed a fixed upstream allele
order is explicitly rejected; it misclassified the PNPLA3 I148M sentinel. The
accepted pass correctly classifies both PNPLA3 I148M and TM6SF2 E167K as
protein-altering and assigns 100% of normalized PIP mass to a mutually
exclusive consequence category.

| Candidate quantity | Direct MASLD/liver fat | Liver enzymes |
|---|---:|---:|
| Reliable credible-set instances | 207 | 6,208 |
| Mean protein-altering PIP mass | 15.77% | 3.26% |
| Mean canonical-splice PIP mass | 0.04% | 0.00% |
| Mean synonymous/UTR PIP mass | 8.46% | 5.11% |
| Mean other-noncoding PIP mass | 75.73% | 91.63% |
| Locus-balanced protein-altering PIP mass | 22.37% | 3.39% |
| Locus-balanced canonical-splice PIP mass | 0.31% | 0.00% |
| Locus-balanced synonymous/UTR PIP mass | 4.71% | 3.57% |
| Locus-balanced other-noncoding PIP mass | 72.61% | 93.04% |
| Noncoding PIP mass within 2 kb of a GENCODE TSS | 7.87% | 4.46% |
| Noncoding PIP mass in a liver ABC enhancer | 0.23% | 1.57% |
| Noncoding PIP mass open in at least one liver lineage | 7.32% | 4.15% |
| Noncoding PIP mass with unresolved deposited context | 52.11% | 60.43% |

Regulatory-context percentages overlap and must not be summed. Credible-set
instances are descriptive study-by-locus units, not independent biological
replicates. The large unresolved fractions prevent a claim that every
noncoding signal has a localized regulatory mechanism.

#### Family multiplicity correction (2026-08-17)

Six direct-versus-enzyme comparisons are made on the same groups. They are one
family and are corrected together. Source:
`GWAS/finemapping/results/fig2_multiplicity/20260817T161128Z/`. Percentages here
are cluster-weighted and so differ from the instance-weighted table above.

| Test | Direct | Enzymes | Ratio | Permutation p | **BH q** | Survives |
|---|---:|---:|---:|---:|---:|:--|
| Coding PIP mass (protein-altering + splice) | 22.54% | 3.44% | 6.55 | <1e-4 | **6.0e-4** | yes |
| Unresolved deposited context | 39.13% | 73.55% | 0.53 | 2.0e-4 | **6.0e-4** | yes |
| Promoter-proximal (within 2 kb of a TSS) | 14.57% | 3.72% | 3.92 | 0.0021 | **0.0042** | yes |
| Resolved, no context | 38.52% | 18.70% | 2.06 | 0.0057 | **0.0085** | yes |
| Open in at least one liver lineage | 9.80% | 5.43% | 1.81 | 0.185 | 0.222 | no |
| **Liver ABC enhancer** | 0.10% | 0.93% | 0.11 | **0.398** | **0.398** | **no** |

⛔ **The ABC enhancer depletion is retired.** Its bootstrap p is 2e-4 but its
permutation p is 0.398. ABC masses are tiny and heavily zero-inflated, so the
difference is reliably non-zero while an arbitrary 22-cluster subset reproduces
it four times in ten. Direct-trait loci may not be described as ABC-depleted.

Two methodological points this correction established, both binding on future
work in this section.

1. **Which null.** The published p is a bootstrap p,
   `2*min((diff<=0).mean(), (diff>=0).mean())`, which asks whether a difference
   clears zero given within-arm sampling variability. A claim about trait class
   needs a permutation null that shuffles the class label. The two disagree in
   both directions across this family: ABC passes the bootstrap and fails the
   permutation; coding does the reverse under a mis-specified statistic. Both are
   now reported per test and BH is applied to the permutation p.
2. **Which statistic.** A permutation statistic of `|log(ratio of means)|` is
   invalid on these zero-inflated masses. 74.1% of enzyme clusters carry zero
   coding mass, so a random subset lands 7× below the pooled mean easily and never
   7× above; over 20,000 shuffles, 100% of exceedances were low-side. Use the
   difference of fractions.

The groups are **operational 1-Mb clusters** of credible-set sentinel positions,
not independent loci: no LD is consulted in forming them, a 1-Mb window can span
more than one independent signal, and long-range LD can link clusters.

**Ascertainment.** The coding concentration is not explained by direct-trait GWAS
being underpowered, though they are. Direct clusters sit at twice the
standardized effect size of enzyme clusters (median β_std 0.0368 versus 0.0183)
and come from studies with roughly one seventh the effective sample size (21,893
versus 160,000), because the direct arm is mostly case-control with few cases.
Matching on β_std (standardized mean difference −0.002) moves the ratio only from
6.87 to 5.32, and adjusting for it attenuates the trait-class coefficient by 2.2%
(0.1258 to 0.1230, cluster-robust over 271 clusters); adding effective sample
size raises it to 0.1464. Source:
`GWAS/finemapping/results/ascertainment_control/20260817T152705Z/`. This controls
discovery power at the lead variant only, not fine-mapping resolution, LD
reference panel, or ancestry composition.

The accepted five-panel render is
`figures/candidates/noncoding-resource-2026-08-12/noncoding-dna-v5/`.
It contains a Figure 2 consequence panel, a Figure 2 nonexclusive context
matrix, a Figure 5 six-lineage accessibility matrix, a supplementary
30-variant accessibility heatmap, and the complete credible-set noncoding-mass
distribution. Each PDF has a read-only source table. Structural validation,
sentinel checks, Ghostscript parsing, unit tests, Resource-scope validation,
and manual visual review passed in SLURM job `19783326`.

These values are validated candidates, not synchronized manuscript numbers.
The DNA-only analysis may proceed without COLOC. Target-transcript biotype,
final examples, and genetics-dependent Gene Catalog entries require a rebuild
against the COLOC release adopted on 2026-08-17; that rebuild has not been run.

The title-free locus-balanced consequence composition was adopted as the
working main Figure 2D under explicit user approval on 2026-08-13. Its source
tables retain the instance-weighted estimates, cluster-bootstrap intervals,
and 1-Mb grouping provenance. This panel-level adoption does not adopt the
corrected COLOC release or convert these values into synchronized manuscript
numbers.

## Inherited genetic variation and genetic–state overlap

The following values are measured on the **corrected COLOC release adopted to
canonical on 2026-08-17** (`GWAS/finemapping/results/susie_coloc/`). The July
release is superseded:

| Quantity | Adopted value | Qualification |
|---|---:|---|
| Tier-1/2 GWAS strata | 35 | Liver-specific/direct and accepted proxy strata; the registry was not changed by the rerun |
| Primary SuSiE-nominated genes | **462** | **Adopted release (2026-08-17); supersedes the July value of 473** |
| Jointly testable background | **14,920** | **Recomputed 2026-08-17**; genes assayed in both the bulk contrast and the adopted colocalization universe. Rebuilt from the adopted release, not carried over from July's 14,931 |
| Jointly testable genetic genes | **437** | **Recomputed 2026-08-17** from the adopted 462; supersedes 447 |
| Overlap with disease-state set | **20** | **Recomputed 2026-08-17**; canonical threshold (`padj<0.05`, `\|log2FC\|>0.50`), against **31.2 expected**; Fisher OR 0.616, nominal p = 0.037 |
| Non-overlap | **95.4%** | `(437 − 20) / 437`; version- and threshold-specific, not a universal biological constant |
| Overlap / non-overlap under the TREAT comparator | 33 / 92.5% | `(437 − 33) / 437`; the retired canonical, retained as the calibrated arm. ⚠ **Not reproduced by the 2026-08-17 interface recompute**, whose outputs carry the canonical threshold only; the producing run is unrecorded. Verify before use |
| Union of named genetic evidence | **1,013** | **Adopted release (2026-08-17); supersedes the July value of 1,030** |
| Tier-1/2 ABF-nominated genes | **822** | **Adopted release; supersedes 862** |

Source: `RNA-seq/results/genetics_state_interface/20260817T133951Z_promoted_coloc/`
(job `20871359`), which rebuilds the background from the adopted release rather
than reusing July's, and reports the official 27,638-gene read-count substrate
and the 23,370-gene five-cohort fragment substrate side by side without pooling.

**Do not upgrade the depletion to a finding.** The nominal p = 0.037 is one test
in a family, and the same test on the intended fragment substrate is OR 0.661,
p = 0.101 (437 → 428 genetic, 20 → 18 overlap). Nothing here establishes that
colocalized genes are buffered or constrained: the expression-matched ladder is
flat, and the threshold-free SuSiE estimate is rho = 0.016 (p = 0.34). The ABF
threshold-free correlation looks strong (rho = −0.066, p = 8.6e-16) and is an
expression artifact — partialling out `AveExpr` leaves rho = −0.016 (p = 0.051),
and on absolute effect it vanishes (p = 0.34). See `EXPRESSION_CONFOUND_CHECK.md`
in the run directory; report the SuSiE arm and the matched ladder instead.

**Post-rerun correction (2026-08-17).** The corrected 50-study × 22-chromosome
COLOC rerun is COMPLETE and PROMOTED to canonical. Measured directly from
`GWAS/finemapping/results/susie_coloc/gene_level_coloc_tier12.csv`, named genes
only (one blank-symbol row excluded, as before):

| | July | post-rerun | Δ |
|---|---:|---:|---:|
| Tier-1/2 SuSiE > 0.5 | 473 | **462** | −11 |
| Tier-1/2 ABF > 0.5 | 862 | **822** | −40 |
| Tier-1/2 SuSiE-or-ABF union | 1,030 | **1,013** | −17 |
| full-portfolio SuSiE > 0.5 | 736 | **745** | +9 |

The tier-restricted counts FALL while the full portfolio RISES — the recovered
eQTL fits are spread across all 50 studies, whereas tier-1/2 draws on 35. The
SuSiE arm churns substantially: **147 tier-1/2 genes change their call (79 out,
68 in)**, driven by retaining ~14.1% of variants that were previously dropped at
every locus as strand-ambiguous. All 7 QA guards pass against the adopted
canonical, including C8 positive/negative control recovery.

⚠ **Downstream status is now split.** Main Figure 2 A-I has been rederived
from the adopted release and explicitly adopted as canonical. Still stale
against the adopted release as of 2026-08-17: the atlas rebuild chain
(`RNA-seq/run_atlas_rebuild_chain.sh`),
`scripts/figures/fig4a_overview_candidates.py`, evidence-class sizes,
genetics-dependent Gene Catalog fields, other worked-example consumers, and
the ~30 single-cell/spatial/ATAC/perturbation consumers enumerated in
`docs/archive/progress/2026-08-17_coloc_rerun_promoted_to_canonical.md`. The
`relayreg`/`orientsub` products are stale for a different reason: they ran on
2026-08-15 against the rerun before the 145 recovered eQTL fits were merged, so
they completed silently on a portfolio headline of 712. Rebuilt on 2026-08-17:
the portal ancestry export
(`masld-atlas-v2/public/data/coloc_*`), the Figure 1 GWAS panels, and canonical
Figure 2 A-I.

COLOC adoption authorizes the genetic counts and overlap arithmetic above; the
separate round-2 approval authorizes the canonical Figure 2 A-I examples. It
does not authorize evidence-class sizes or Gene Catalog entries, which still
wait on rebuilds that have not been run.

### Adopted Figure 2 A-I round-2 release

The adopted A-I set was produced in
`figures/candidates/fig2_round2_promotion_20260817T211348Z/`. Job `20879541`
rerendered all nine panels together from adopted inputs and passed
source-cardinality, exact fractional-point page-box, one-page PDF,
embedded-font, Ghostscript, prohibited-label, raster-equivalence, visual, and
Resource-scope validation. The resulting PDFs and sidecars are canonical under
`figures/main/fig2_genetics/panels/`. Panel B uses the full 0--1 axes and a
portrait plotting area measuring 1.488 by 1.970 inches; only 0.018 inch remains
between the x-axis-title baseline and the page bottom.

The variant consequence annotation was rerun directly from the adopted COLOC
aggregate in job `20871989`: 19,851 unique variants were classified, and all
five prespecified coding, intronic, promoter, and intergenic sentinels passed.
For the 1,013 named Tier-1/2 multi-signal-or-single-signal COLOC union, the
strongest lead-variant annotation per gene is:

| Consequence | Genes |
|---|---:|
| Coding | 36 |
| UTR | 73 |
| Promoter | 107 |
| Intron | 455 |
| Intergenic | 342 |

The multi-signal ancestry partition is 352 EUR-only, 70 shared with EUR, 26
non-European-only but below the lead-variant genome-wide-significance threshold, and
14 non-European-only genes passing that threshold (13 EAS, 1 AFR). These are Figure
2 source-table values, not final evidence-class or Gene Catalog memberships.

The cross-ancestry matrix is exhaustive under its displayed rule, not
hand-selected: 13 genes have PP.H4 > 0.5 in four of the five ancestry panels,
and no gene has support in all five. The roster is GGT1, PANX1, EPHA2, C2orf16,
ACTG1, SHROOM3, MLIP, MSL2, EFHD1, GOT2, GLDC, HSCB, and CHEK2.

Late supplementary Figure 2 panels J, L, M, and O--R were independently
rederived in candidate
`figures/candidates/figS2_JR_promoted_coloc_20260817T204537Z_20878109/`.
Job `20878109` passed source, PDF, embedded-font, 6-point-text, prohibited-label,
and Resource-scope checks. Panels K and N are withheld rather than carried
forward: K has 31,036 historical rows without adopted COLOC gene/posterior
linkage, and the 35 adopted GWAS strata contain no stage-associated GWAS stratum
for N. Canonical supplementary PDFs remain unchanged pending separate adoption.

Early supplementary panels A--D and F--I were rederived and independently
validated through job `20879572`; E is withheld because the adopted BBJ_GGT--ACADS
has no multi-signal support. The principal changes are EFHD1--BBJ_AST PP.H4
0.84 to 0.650043; method-concordance rho 0.5453 to 0.5134; full-portfolio
multi-signal genes 736 to 745 and single-signal genes 1,234 to 1,181; and
tier-1/2 multi-signal threshold counts 473/369/267 to 462/345/252 at PP.H4
0.5/0.8/0.9. These supplementary files remain candidates; main A--I alone was
authorized for round-2 canonical adoption.

COLOC is shared-signal evidence. It nominates genetically anchored candidates;
it does not by itself identify a causal gene or prove mediation.

## Assay observability and evidence-state accounting

These results define the interpretive boundary around the inherited-variation
and multicellular-remodeling maps. Apparent negatives are assigned states that
reflect what each assay could actually establish rather than being counted as
evidence against a gene or program.

| Boundary | Result | Evidence interpretation |
|---|---:|---|
| CosMx complete-program observability | 12 / 117 programs | The other programs cannot be treated as negative whole-program tests |
| GeoMx complete-program testability | 0 / 117 programs | `untestable` under the accepted direction and normalization rules |
| GSE287826 complete-program testability | 0 / 117 programs | `untestable` without an authoritative AOI-to-donor key |
| ATAC jointly testable peaks | 39,914 | Denominator after cohort and lineage filters |
| ATAC `indeterminate` peaks | 38,066 / 39,914 | Available cohorts do not resolve a replicated accessibility state |
| ATAC `source_dependent` peaks | 1,844 / 39,914 | Interpretation depends materially on source |
| Catalog strict `tested_negative` calls | 0 | No current row satisfies every adequacy requirement for a negative claim |

These counts are not validation failures. They show that modality-vote counting
would silently treat targeted-panel absence, missing donor joins, source
dependence, and unresolved replication as equivalent negative evidence. The
Resource instead exports the state, its reason, the unresolved alternative,
and the deterministic next-experiment rule in machine-readable form.

## MASLD Gene Catalog

| Quantity | Candidate value | State |
|---|---:|---|
| Genes | 30,754 | Candidate schema build |
| Evidence rows | 61,508 | Candidate schema build |
| Prespecified programs | 117 | Consistent with program registry |
| Program memberships | 6,633 | Candidate export |
| Contexts | 14 | Candidate export |
| `tested_negative` calls | 0 | Expected under current conservative source/coverage rules |
| Catalog v2 documented fields | 19 | Independently validated candidate requirements |
| Molecular-object follow-up routes | 5 | Deterministic, outcome-blind rules |

These counts demonstrate that the Gene Catalog data structure and follow-up
logic are feasible. Genetics-dependent states must be rebuilt after corrected
COLOC; the candidate is not an adopted paper release.

The noncoding-aware v2 requirements are read-only at
`Analysis/Multimodal_Program_Projection/candidates/gene-catalog-v2-contract-2026-08-11/`
and passed an independent checksum, field, vocabulary, route, and separation
review in `gene-catalog-v2-contract-review-2026-08-11/`. It contains the public
schema, a human-readable data dictionary, and five follow-up rules, but no gene
outcomes. Regulatory DNA, transcription through a lncRNA locus, and the mature
lncRNA molecule remain separate candidate objects. Genetics-dependent entries
remain on hold until they are rebuilt against the COLOC release adopted on
2026-08-17; the compute dependency is cleared, the rebuild is not done.

Gene Catalog vocabulary:

| State | Meaning |
|---|---|
| `supported` | The prespecified assay supports the evidence claim |
| `discordant` | An adequately observed assay points in a conflicting direction |
| `tested_negative` | The claim was adequately tested and unsupported |
| `untestable` | Coverage or design could not test the claim |
| `not_applicable` | The assay does not address that biological question |
| `source_dependent` | Interpretation depends materially on one source or mapping assumption |
| `indeterminate` | Available data do not resolve the question |

## Molecular and physical context

### Validated single-cell ATAC Context v3 candidate (not adopted)

The isolated candidate
`Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/`
uses GSE244832 and GSE281367 only. `NON_GENETIC_READY`, `INTEGRATION_READY`,
and `FIG4F_READY` passed independent validation. Canonical figures, the Gene
Catalog, and release IDs were not changed. The genetics check wrote
`blocked_upstream_release` and no biological genetics rows because at the time
of that run no adopted corrected COLOC manifest existed. Adoption happened on
2026-08-17, so the variant-posterior replay is now runnable and has not been run.

| Candidate quantity | Value | Qualification |
|---|---:|---|
| Consensus peaks | 676,658 | Condition-blind, fixed-width 500-bp, non-overlapping hg38 intervals |
| Hepatocyte / stellate / macrophage peaks | 206,874 / 173,191 / 142,456 | Broad macrophage lineage |
| Cholangiocyte / T-NK peaks | 113,584 / 40,553 | GSE281367 T/NK is a combined source label |
| Peaks supported by both cohorts | 347,509 | 131,663 hepatocyte; 110,048 stellate; 55,300 macrophage; 45,901 cholangiocyte; 4,597 T/NK |
| Sparse count matrices | 16 | Ten consensus matrices plus six old-peak-space sensitivity matrices |
| Direct fragment recounts | 10 of 10 exact | One seed-fixed donor-peak combination per cohort and lineage; counting unit is one deduplicated fragment-file record |
| Official chromVAR rows / tested contrasts | 8,790 / 7,032 | JASPAR2024 CORE vertebrate motifs; eight contrast-testable lineage-by-cohort families |

Matched MASH-versus-normal DA used donors as replicates, separate cohort fits,
TMM normalization, `voomWithQualityWeights`, and BH correction within each
lineage-by-cohort family. GSE244832 MASL donors were excluded. The prespecified
default condition-blind `filterByExpr` rule left no GSE244832 macrophage peaks,
so macrophage DA is `untestable`; the filter was not relaxed after observing
that result.

| Lineage | Donors, GSE244832 normal/MASH | Donors, GSE281367 normal/MASH | Joint peaks | Significant, GSE244832/GSE281367 | Strict replicated | Directional-enrichment BH q | Main dynamic criterion |
|---|---:|---:|---:|---:|---:|---:|---|
| Hepatocyte | 5 / 9 | 6 / 6 | 30,650 | 1 / 1,016 | 0 | 2.64×10⁻¹¹⁰ | Did not pass the sensitivity-stable replicated-peak criterion |
| Stellate | 5 / 9 | 6 / 6 | 9,264 | 820 / 15 | 4 | 4.58×10⁻⁵ | Passed the candidate eligibility criterion |
| Macrophage | 4 / 8 | 6 / 6 | 0 | 0 / 0 | 0 | 1 | Untestable after default filtering |

Across the 39,914 jointly testable peaks, evidence states were 4 `supported`,
1,844 `source_dependent`, and 38,066 `indeterminate`, with no `discordant`
peak. All four supported peaks were stellate-lineage peaks, retained direction
under both leave-one-donor-out analyses, and did not reverse under the
GSE244832 batch or broad MASLD sensitivities. This satisfies the prespecified
candidate dynamic-main eligibility rule, but it is not yet an adopted main-panel
claim and does not prove a causal gene or causal cell type.

The prespecified-program projection contains exactly 234 rows, 117 programs per
cohort. Program scores were testable for 116/117 programs in GSE244832 and
117/117 in GSE281367. MASH contrasts were testable for 78/117 and 117/117,
respectively; GSE244832 cholangiocyte and T/NK contrasts remain explicitly
untestable. Across the paired cohort rows, 156 were `indeterminate` and 78
`untestable`; no program met the strict two-cohort replicated-program rule.
No official chromVAR motif met the strict q<0.05, same-direction replication
rule. chromVAR remains supplementary, and these results do not define a
regulatory network or a tested-negative motif family.

The Figure 5F source, after filtering, contains 30 exact cells: two cohorts × five
lineages × promoter-gene coverage, program-score testability, and contrast
testability. Its one-page PDF passed source rederivation, checksum, font,
structural, and visual review. Candidate-only Gene Catalog adapters preserve
30,089 gene-lineage-state promoter rows and 13,266 gene-program-cohort rows;
the optional integration maps at least one ATAC v3 field to 14,510 of 27,187
catalog genes without changing any protected L8 value or `layers_active` vote.

Static peak observability and disease-associated DA are distinct. Distal peaks
are not assigned to genes, and the GSE281367 T/NK context remains
source-dependent. The adopted-COLOC replay was executed and recorded on
2026-08-20: 5,362 signal pairs and 22,356,697 variant-posterior rows reproduce
816 adopted colocalized gene-study pairs over 462 genes. Canonical Figure 5D
uses one deterministically selected primary signal pair per gene-study pair and
intersects normalized SNP.PP.H4 posterior mass with static peaks from both ATAC
sources across five shared lineages. The display separates 29 direct
MASLD/liver-fat pairs from 787 liver-enzyme pairs. All primary pairs retain at
least 0.95 mapped posterior mass, so none is untestable; most lineage contexts
are partial rather than both-source accessible. Peak overlap is regulatory-DNA
context, not target-gene assignment or a functional test. Pair-level PP.H4
reproduced the adopted release within the documented tolerance, but no
adopted per-variant posterior vector exists for an external equality check.

### Validated spatial candidate (not adopted)

The read-only candidate `spatial-resource-candidate-2026-08-11` passed 32
analysis/requirement checks and the six-panel Figure 5 source/PDF checksum check.
It is suitable for adjudication and downstream candidate builds, but it is not
the current official Resource release.

| Candidate quantity | Value | Interpretation |
|---|---:|---|
| Registered dataset or analysis decisions | 14 | Every source is retained, supplemental, skipped, or dropped explicitly |
| Assays requiring complete coverage | 7 | GSE192741, Vu, HMSMA, Yakubovsky, CosMx, GeoMx, and GSE287826 |
| Prespecified-program coverage rows | 819 | Exactly 117 programs × 7 assays; no outcome statistics or winner selection |
| Gene-context rows | 23,576 | Assay-facing measurement/testability records |
| Confirmatory effect rows | 20 | Exactly 2 predeclared programs × 10 registered source contexts, including untestable and not-applicable rows |
| HMSMA technical spots / arrays | 130,097 / 35 | Technical coverage only; neither is a biological donor count |
| Candidate Figure 5 output files | 6 PDFs / 6 source tables | All PDFs are one page and checksum-linked to their exact source table |
| Downstream spatial figure output files | 4 PDFs / 4 source tables | Figure 1 capability, Figure 6 Gene Catalog, all-program coverage, and HMSMA descriptive panels validated in candidate `spatial-figure-integration-candidate-2026-08-11-r2` |

Coverage status across the seven required assays:

| Assay | Observable | Partial | Untestable | Biological-unit qualification |
|---|---:|---:|---:|---|
| GSE192741 Visium | 114 | 3 | 0 | 5 sections collapsed to 4 donors for inference |
| Vu et al. Visium | 114 | 3 | 0 | 10 arrays; donor mapping unresolved, so results remain source-dependent |
| HMSMA Visium | 114 | 3 | 0 | 35 arrays; donor/clinical mapping unresolved |
| Yakubovsky normal-liver spatial | 114 | 3 | 0 | Normal-liver reference, not MASLD validation |
| Govaere CosMx | 12 | 97 | 8 | Targeted-panel observability; no whole-program inference when coverage fails |
| Govaere GeoMx | 0 | 0 | 117 | Dropped from inference because direction/normalization remain unresolved |
| GSE287826 | 0 | 0 | 117 | No authoritative 30-AOI-to-15-donor key |

The confirmatory family was recorded before execution and contains only the Stromal
ECM and Ductular injury programs. In GSE192741, Stromal ECM had donor-first
centered-residual Moran's I 0.1390 (matched-null SD 0.02144, P=0.0001,
q=0.0002) and passed equal-weight and leave-top-weighted-gene sensitivity,
therefore `supported`. Ductular injury had Moran's I 0.02587 (matched-null SD 0.01677,
P=q=0.0755) and is `indeterminate`, not negative.

In Vu, Stromal ECM had Moran's I 0.1471 (matched-null SD 0.02048, P=0.0001,
q=0.0002) and passed both sensitivity signs within source. Its evidence state
is nevertheless `source_dependent` because the ten arrays cannot be joined to
biological donors. Ductular injury had Moran's I 0.01922 (matched-null SD 0.01557,
P=q=0.1158) and is `indeterminate`.

The result becomes more discriminating at the physical-unit level. Residual
Moran's I was positive for Stromal ECM in 4/4 GSE192741 donors and 10/10 Vu arrays,
but it was also positive for Ductular injury in 4/4 donors and 9/10 arrays. Thus neither
positive sign nor cross-unit consistency distinguishes the programs. The
matched-gene calibration does: Stromal ECM exceeded its expression-, detection-,
biotype-, and lineage-correlation-matched background in both sources, whereas
Ductular injury overlapped its expected spatial background. The paper-level finding is
therefore that **generic spatial autocorrelation can mimic program
transportability, and an assay-matched null is required to separate the two.**
This is not a negative Ductular injury result, and the Vu contrast remains
source-dependent rather than donor-level replication.

Canonical Figure 5F pairs this calibration with a registered H&E overlay for
Guilliams et al. section JBO019 (steatotic donor H37), selected outcome-blind
because its 1,164 graph-eligible spots equal the dataset median. The displayed
Stromal ECM score uses 26 measured fixed members, retains 0.5475 of the original
L1 weight, and is residualized for all 16 cell2location abundance factors plus
library size and detected genes. Re-derivation from the spot matrix reproduced
the recorded residual Moran's I exactly (0.0582734890). The overlay adds no test:
it is a representative physical view without pathologist region masks and does
not authorize localization to a named histologic structure.

HMSMA is deliberately descriptive: across 35 arrays, median per-array Moran's
I was 0.1619 for Stromal ECM (range 0.03844–0.4849) and 0.01237 for Ductular injury (range
-0.02244–0.08640). No population P value, disease contrast, fibrosis/NAS model,
histologic localization, or MALDI colocalization is authorized without the
authoritative joins.

The candidate portal exports the registry, gene context, 117-program coverage,
and two-program effects as release-linked products. Legacy numeric
`s5_spatial` fields remain compatibility-only and are excluded from rankings,
active-layer counts, evidence-strength bars, and default decisions. Candidate
gene and program pages consume the categorical tables and expose testability,
source dependence, and biological versus technical units; the files remain
unstaged from the public data directory until synchronized adoption.

| Layer | Current conclusion | Qualification |
|---|---|---|
| Tissue protein | Useful orthogonal transportability evidence | Must remain assay-native and abundance-aware |
| Secreted protein | Supplementary independent context | Cannot be collapsed with tissue proteome into one score |
| Chromatin/regulatory context | Supports cell-context observability | Does not convert a colocalized candidate into a causal gene |
| Single-cell localization | Donor-aware cell-labelled context | Cells and runs are never biological replicates; ambient-sensitive programs do not establish lineage origin |
| Spatial organization | Physical-context support in accepted cohorts | Existing cohorts do not provide a complete MASLD control/donor/histology/metabolite design |
| Yakubovsky | Zonation and local lipid-context reference | Normal liver reference, not MASLD disease validation |
| GSE287826 | Could not be tested inferentially | No authoritative public donor key |

The current canonical Figure 5 contains all six adopted panels and the vector
assembly `fig5_molecular_context.pdf`. Panel 5D carries the variant-consistent
ATAC context; 5E carries assay-native program context; and 5F carries the
all-cell-adjusted matched-gene calibration plus one registered Stromal ECM H&E
overlay. The calibration carries the spatial inference; the map only shows the
score's physical distribution in an outcome-blind representative section.
Coordinate-only maps of both programs are retired from the active panel tree.
Candidate `fig5-completion-2026-08-24-v6` passed 80 panel/source/assembly checks
in job `21068950`; post-copy multimodal, Ghostscript, and Resource-scope checks
passed in job `21068987`.

The release-linked downstream candidate adds four one-page panels without
changing any spatial call: Figure 1 separates biological from technical units
and displays source/join/coverage checks; Figure 6 contrasts the supported
Stromal ECM program, indeterminate Ductular injury program, and untestable HMSMA clinical
branch; the supplement reports all 819 program-by-assay coverage rows and all
70 HMSMA array-by-program descriptive rows. Each panel has an exact source
table and checksum manifest. The candidate is validated and visually reviewed,
but remains unpromoted.

The visually accepted publication-layout candidate is
`spatial-publication-figures-candidate-2026-08-11-r4`. It contains 13 individual
one-page PDFs: 11 were simplified and two already-clear panels, Figure 5E and
the all-117-program coverage heatmap—were retained byte-for-byte. Across the 11
redesigned panels, plot text fell to 37.7% of the prior total. This is a visual
change only: exact estimates, evidence states, unit definitions, and complete
families rederive from the recorded source tables. The focused cell-context panel
shows the Stromal ECM–fibroblast and Ductular injury–cholangiocyte comparisons needed for the
paper message; the complete 64-row factor summary remains in the validated
cell-context candidate.

The result-led successor is
`spatial-impact-figures-candidate-2026-08-11`. It retains 12 r4 panels
byte-for-byte, replaces the main 5A and 5F with the compact calibration logic
and matched-null discrimination panels, and carries the former source-state
matrix into the supplement, for 14 one-page PDFs total. Its four-row
transportability profile independently rederives the prespecified program hashes,
9,999-draw nulls, unit-sign counts, accepted evidence states, and all-cell
composition sensitivity. It adds no inferential family or spatial score and
writes no canonical output. This candidate passed source, checksum, page,
font, claim-boundary, and visual review and is the preferred spatial assembly
family pending scientific adjudication and synchronized adoption.

### Validated spatial diagnostics (not adopted)

Four bounded candidate lanes now test the spatial result without expanding the
two-program confirmatory family or changing Figure 5 evidence states:

| Diagnostic candidate | Validated content | Paper role |
|---|---|---|
| `spatial-full-composition-sensitivity-candidate-2026-08-11-r2` | 82 checks; all 16 shared cell2location q05 factors plus the same two QC covariates; four dataset-by-program rows; full-rank 19-column designs | Supplementary falsification of composition dependence |
| `spatial-unit-null-diagnostics-candidate-2026-08-11-r2` | Three PDFs and three exact source tables covering physical-unit heterogeneity, matched-null calibration, and weight sensitivity | Makes donor collapse, source dependence, and indeterminacy auditable |
| `spatial-cell-context-attribution-candidate-2026-08-11-r2` | 480 within-section/array correlations, 448 donor/array reporting rows, 64 complete factor summaries, and one full-family heatmap | Descriptive context attribution only; no spot-level P values |
| `hmsma-mask-sensitivity-candidate-2026-08-11-r2` | 210 rows = 35 arrays × 2 programs × 3 UMI thresholds; all 35 input H5ADs checksum-recorded | Tests stricter proxy masks while donor and clinical joins remain unresolved |

`spatial-diagnostic-bundle-candidate-2026-08-11` is the single validated index
for these products: five accepted components, six one-page PDFs, 58 candidate
output files, and 27 executable code dependencies are checksum-linked, with
canonical writes disabled. Superseded visual candidates and failed validators
are not registered in the bundle.

With all 16 abundance factors included, the Stromal ECM centered spatial effect
fell from 0.1390 to 0.04050 in GSE192741 (`q=0.001`) and from 0.1471 to
0.05034 in Vu (`q=0.0006`). It therefore remains supported within GSE192741
and within-source robust in Vu, but most of the original magnitude is
composition-linked and Vu remains source-dependent. Ductular injury remained
indeterminate: 0.01532 (`q=0.0793`) in GSE192741 and 0.002654 (`q=0.3158`)
in Vu after full adjustment. The accepted supplementary attenuation panel is
`spatial-full-composition-figure-candidate-2026-08-11-r4`; its source table
retains all four rows and the noncausal claim boundary.

The full 16-factor descriptive heatmap localizes the Stromal ECM program most
strongly to fibroblast-rich context: median within-unit Spearman rho was 0.522
across four GSE192741 donors and 0.398 across ten Vu arrays, with all reporting
units positive. Ductular injury–cholangiocyte covariation was weaker (rho 0.101, 3/4
GSE192741 donors positive; rho 0.174, 10/10 Vu arrays positive). These are not
independent validations because program scores and cell2location abundances
come from the same Visium RNA matrix.

The physical-unit panel shows why aggregation matters: H35's two Stromal ECM
sections had residual Moran's I of approximately 0.080 and 0.328, but they
contribute one donor-collapsed value. No section- or array-level P value is
reported. In HMSMA, Stromal ECM remained positive in all 35 arrays under 1,000- and
2,000-UMI proxy masks (rank rho 0.977 and 0.937 against the 500-UMI baseline;
sign agreement 1.00). Ductular injury was less stable at 2,000 UMIs (rho 0.750; sign
agreement 0.80), reinforcing its indeterminate status. This sensitivity cannot
replace an authoritative tissue mask because the available H5ADs were already
filtered at 500 UMIs.

The CosMx and Yakubovsky audit supports no expansion: CosMx retains only its
source-native, array-level IL32–macrophage vignette, and Yakubovsky remains a
three-donor normal-liver boundary with both programs indeterminate. No new
panel or inferential claim is warranted for either source.

## Public functional stress tests

| Workstream | Result for this paper | Disposition |
|---|---|---|
| Myojin palmitate survival | Did not confirm the proposed Resource-level enrichment | Complete supplement; survival is not interchangeable with lipid storage |
| Acute induction/localization/reversal | Strict multi-assay adoption criteria did not pass | Complete supplement with valid nulls and source checks |
| Risk/state double dissociation | Full two-arm criteria did not pass | Complete supplement; no functional-partition headline |
| Chronic reversal/risk bridge | Full risk-to-state architecture criteria did not pass | Complete supplement; no causal relay claim |
| Multicellular response topology | Did not produce a source-independent main-paper mechanism | Complete supplement; no replacement search |

These analyses strengthen the Resource by defining what acute public models can
and cannot validate. They do not justify reopening mechanism hunting or moving
the paper away from its Resource identity.

## Historical corrections that must not return

| Retired statement | Correct interpretation |
|---|---|
| 1,277/1,259 human samples as current | Pre-recount documentation; corrected source registry is 1,281 count-matched and 1,257–1,260 QC-passing depending on release |
| 54 disease-significant Hotspot modules | Run-level pseudoreplication; donor-level result was 22 in that historical analysis, while the present prespecified selection uses the Plan-20 criteria above |
| 58 Andrews donors with F-stage | 58 sequencing runs mapping to 46 biological donors |
| cross-cohort inferred F-stage accuracy 0.638 | Run-level/leaked result; not authorized for cross-cohort fibrosis inference |
| 13 significant stage-progressive ligand–receptor pairs | Zero of 13 passed donor-level family FDR; heatmap is descriptive |
| 290 sex-dimorphic genes | Retired random-slope result; fixed-cohort model yields 8 and the former female-progressor narrative is retired |
| k=6 discrete patient subtypes | Factor subspace is stable, but hard partitions are not; use continuous axes |

## Where the full history lives

The previous 2,220-line number ledger is preserved verbatim in the documentation
archive. Use [technical/NUMBERS_HISTORY.md](technical/NUMBERS_HISTORY.md) to find
historical thresholds, superseded releases, and source files. New prose and
figures should cite this compact ledger unless a technical audit explicitly
requires the archive.

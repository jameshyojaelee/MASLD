## Bulk transcriptomic remodeling is connected to genetics through exact signal pairs

> Working Figure 3 prose. The five-cohort stage rebuild is validated, and the
> adopted-release Figure 3G genetics–state interface is a working main panel.
> Figure 3F remains pending exact-pair rendering and synchronized release
> adjudication.

We characterized established disease tissue on one synchronized five-cohort
fragment-count substrate. The pooled model contains 844 participants, including
157 controls and 687 disease samples, and tests 23,370 genes with fixed cohort
and inferred-sex terms. Panels 3A–C show cohort coverage, documented
fibrosis-by-NAS support, and per-cohort replication separately so that pooled
effects are not mistaken for uniform behavior in every study.

Figure 3C tests whether the pooled result depends on one cohort. Across the
1,347 canonical genes, per-cohort direction agreement ranges from 79.8% to
100%. Leaving out each cohort in turn yields all-gene effect correlations of
0.867–0.992 with the full model and retains 60.7%–98.0% of the canonical threshold.
This directly tests cross-cohort robustness while keeping participants as the
inferential unit.

Figure 3D shows a prespecified phenotypic endpoint rather than changing the
canonical model to improve unsupervised separation. It retains 418 participants
from four cohorts with at least five participants in each arm: 115 strict
controls and 303 participants with NAS at least 5 or fibrosis at least F3.
Histology takes precedence for the one source-label conflict. Strict controls
are source-defined controls with no observed NAS above 1 and no observed
fibrosis above F0; missing histology is allowed only when no observation
contradicts control status. Expression is adjusted for cohort and inferred sex
without preserving endpoint labels. The density and scatter contain the same
participants, with controls in gray and advanced disease colored by documented
fibrosis stage. The 844-participant PCA is shown in Figure S3.

The matched endpoint DEG sensitivity uses the same voom-quality-weighted engine
and complete 23,370-gene family as the canonical model. Its log2 fold changes
correlate at rho=0.938 with the canonical estimates; it recovers 1,307 of 1,347
canonical DEGs and yields 1,336 endpoint-only calls under the independently
applied threshold. This supports the broad disease-state inference without replacing
its five-cohort, 844-participant estimand or downstream membership.

Figure 3E aligns two views against one common reference: numbers of genes with
higher or lower expression at F1–F4 versus F0 and continuous F0-centered NMF
activity. Each contrast uses only participants in cohorts carrying the compared
groups, includes fixed cohort and inferred-sex terms when estimable, and
applies BH correction across the complete 23,370-gene family.
The adjacent contrasts contain 104 F0 versus 105 F1, 105 F1 versus 90 F2, 90
F2 versus 75 F3, and 75 F3 versus 25 F4 participants across three contributing
cohorts.
The continuous NMF track shows centered program activity, not dominant patient
classes. Adjacent fibrosis, NAS, donor-aware hepatocyte, and non-parenchymal
composition views remain in Figure S3, generated from the same source release.

Multicohort differential expression and stage-ordered transcriptomic summaries
in MASLD are established analyses rather than contributions of this paper. A
nine-cohort public bulk meta-analysis of 897 NAS-scored and 856 fibrosis-scored
biopsies reports stage-resolved differential expression and pathway change on
overlapping data, including three of the five cohorts used here (Kumak et al.,
2026; doi:10.64898/2026.08.26.747338), following earlier harmonized liver
resources. Figure 3 accordingly claims no novelty for recovering a multicohort
disease-state gene set or for ordering expression change across histologic
stage. What Figure 3 contributes is the estimand and the unit: one synchronized
fragment-count substrate, biological participants as the inferential unit in
every fit, a complete 23,370-gene family carried through each contrast, and the
connection of that remodeling to the inherited map in Figures 3F and 3G. The
five-cohort substrate contains no repeat-biopsy cohort; where the Resource does
use repeat biopsies, they enter only within-person paired models and never a
cross-sectional contrast.

We retained histology as the clinical severity anchor and asked a separate
question in Figure S3: whether an externally anchored transcriptomic continuum
could resolve molecular variation compressed by fibrosis categories. All five
Resource cohorts were scored and enter descriptive continuum-window displays;
the four cohorts with recorded fibrosis enter stage displays. Two prespecified
co-primary scores were evaluated inferentially in GSE162694 and GSE213621, with
fibrosis stage and inferred sex included in every donor-level model. The
continuum recovered the broad stage-associated direction of transcript,
pathway, NMF, and Hotspot remodeling while retaining substantial stage-adjusted
signal, particularly for multicellular programs. These comparisons use the
same expression substrate and are therefore alignment analyses, not independent
validation. Sliding windows are fixed descriptive summaries only. This
extension complements the stage contrasts and does not change Figure 3, the
canonical DEG definition, or the interpretation of NMF axes.

GNMT, MAT1A, and CYP2C19 make the remodeling result biologically concrete.
All three decrease in the pooled disease-state comparison and in adjusted
liver protein. Their adjacent-stage support is narrower than it first appears.
The F0 group mixes 31 source controls with 73 F0 disease biopsies, and because
control status is nested in cohort the cohort term cannot separate them, so any
F0-to-F1 contrast is partly disease versus control. Refitting on disease-only F0
leaves GNMT decreasing from F1 to F2 (q=0.047) and removes the F0-to-F1 support
for both GNMT (q=0.0085 to 0.216) and MAT1A (q=1.3e-4 to 0.105); CYP2C19 retains
F1 to F2. The point estimates move only 15 to 20%, so this is a statement about
what the bottom of the stage range can resolve, not about the direction of the
change. MAT1A falls below the canonical pooled 0.50 fold-change
floor but passes the interval-null TREAT comparator, so it is retained as part
of the coordinated biological pattern rather than relabeled as a canonical
DEG. Together these data support loss of hepatocyte metabolic function across
RNA, stage, and protein; they do not establish longitudinal decline, a shared
upstream regulator, or a patient subtype.

Within 46 MASLD participants, the three proteins also show positive
rank-residual partial Spearman correlations after adjustment for acquisition
batch, age, BMI, and sex (rho=0.406–0.472; BH q=0.00141–0.00518 over the
three-pair family). The correlations survive added fibrosis-stage adjustment
but not NAS or joint component-histology adjustment. This candidate sensitivity
supports a coherent histologic-severity-linked state, not an autonomous
molecular circuit, and remains excluded from final panels until independent
validation.

Figure 3F connects stage-associated expression to regulatory genetics without
collapsing the evidence types. Eligible rows must be BH-significant in at least
one adjacent fibrosis contrast and belong to the adopted primary colocalized
set. Signed stage effects, exact gene–trait PP.H4, and lead GWAS SuSiE PIP are
shown in separate lanes. PIP appears only for an exact gene–trait–credible-set
link; unresolved linkage is missing rather than zero. PP.H4 and PIP are not
combined into a score. Each comparison retains its evidence state, unresolved
linkage alternative, and next discriminating experiment.

The current official genetics cross-check does not convert the metabolic trio
into an inherited regulatory program. None has primary SuSiE PP.H4 above 0.5;
GNMT has one direct-NAFLD single-signal ABF sensitivity result above 0.5 and is
therefore the provisional bridge between the genetic and state maps. The
adopted corrected COLOC release will determine whether that exact gene–trait
pair belongs in Figure 3F. The RNA–stage–protein state stands independently of
that decision.

Figure 3G asks whether sparse overlap means that inherited regulatory
susceptibility and established disease responsiveness are biologically opposed.
It does not. Across 3,617 genes with finite SuSiE PP.H4 and raw disease-state
log2FC estimates, the threshold-free association is rho=0.014 (P=0.40), shown
without a fitted line. The scatter uses the labels Neither, Transcriptomics sig.,
Genetics sig., and Both; Transcriptomics sig. is the calibrated TREAT family
(`lfc=0.25`, BH q<0.05), not a post-FDR raw fold-change cutoff. In 428 one-to-one
expression-matched pairs, 30 Transcriptomics sig. calls occur only in the
Genetics sig. member and 21 only in the matched control; 377 pairs contain
neither call and none contain both (paired OR=1.43; exact McNemar P=0.26). PP.H4
is unsigned, so this does not test risk-allele direction. It supports
complementary evidence maps without biological antagonism.

The Figure 3F exact-pair matrix and PIP-focused supplementary view remain
fail-closed until their adopted-release renderer and linkage audit pass.
Figure 4 then resolves the prespecified cell-type-specific gene programs and
regulatory activities at biological-donor level.

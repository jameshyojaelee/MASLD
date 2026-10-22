# Distinct model tracks for the MASLD Resource

2026-09-30. Read-only comparison by Codex alternative-model agent. No model fit,
new-source outcome scoring, data acquisition, authority change or scientific
adoption was performed. The parent retains coordination with Claude.

## Correction: complete estimates do not require every gene's biological detectability

This correction supersedes the later historical passages demanding
"complete observed-input" or biological mapping observability for every
modeled gene. The fixed receiving contract retains the **complete 42,163
native estimate-ID axis**, valid gene mapping and the unchanged modeled-gene
sum denominator for both frozen raw arms. It does **not** require every
gene to be biologically detectable or have a positive estimate in every
donor. That stronger requirement was a Codex review overreach, not an
identified Claude specification or released scorer requirement.

Genuine native quantifier rows, including declared structural-zero estimates
for represented but mapping-ineligible short transcripts, may be retained
under an explicitly fixed approximate measurement recipe. Flag those rows
as unobservable under that recipe; neither they nor training zeros establish
biological absence. Missing reference/quantifier IDs cannot be concealed
by invented rows or values, and the private PISCES generator remains
unreconstructed. Complete computational estimates, biological observability,
and source-measurement equivalence are separate claims. Root explicitly
accepted this correction before receiving processing; historical stricter
text below is retained as superseded, not an additional binding condition.

Exact source trail: `2026-09-30_codex_program_model_options.md:394`/410
introduced the complete admitted/full-measured-axis working contract;
this report's historical short-transcript stopping condition escalated it.
Released `score.py:167`–190 checks identifier coverage and a positive
*total* RNA estimate, not individual-gene detectability;
`MODEL_CARD.md:201` describes its 0.60 coverage floor. The new complete-axis
receiving comparison remains stricter about ID completeness than that
release floor, without adding a biological-detectability test.

## Recommendation

Pursue a bounded **within-histology RNA-to-H3K27ac donor-difference model** while
separately resolving **target-conditioned allele prediction** and its external
source eligibility. These answer different research questions from Model A's
shared stage-associated allelic remodeling. The existing RRR and native
AlphaGenome predictor are strong comparators, not new contributions to rename.
Neither alternative currently meets the user's independent-validation goal.

The immediate useful work is to fix the donor-difference objective and inspect
metadata-only pair availability before fitting; in parallel, finish source and
target geometry contracts for an allele model. Repeating the existing unweighted
all-pairs regression, same-region sequence interaction or adapter combination
would not resolve the identified scientific questions.

## 1. Predict chromatin differences between comparable participants

**Research use.** Given liver RNA from two participants with the same recorded
histology, predict how their measured residual H3K27ac profiles differ. This is
cross-sectional molecular resolution, not longitudinal progression, a personal
genome effect, treatment-response prediction or a causal regulatory mechanism.
The endpoint remains residual H3K27ac log2 CPM, after the existing fold-training
fit of four signal-concentration descriptors. It is not an absolute H3K27ac
profile restored without those descriptors.

**Existing evidence, verified from saved exposed-development summaries.**

| Quantity | Current evidence | What it establishes |
|---|---|---|
| Biological development n | 99 distinct GSE267145 participants, same-sample/different-aliquot RNA and H3K27ac | Within-source model development; regions are repeated responses |
| Held-participant profile skill | Global RRR 0.119192; global-plus-local 0.141329; histology/sex/RNA-technical ridge 0.019708 | Prediction beyond this specified linear comparator |
| Exact-histology/sex/fold comparison | 22 disjoint pairs, 44 participants; local difference skill 0.11070, global 0.08859 | A useful existing difference endpoint, already inspected |
| Coverage of that pair comparison | 20 female F0 pairs, one male F0 pair, one F1 pair; no F2–F4 pairs | No demonstrated within-grade capability at advanced fibrosis |
| Existing sequence-by-RNA interaction | Original-target skill increment −0.0000159, interval [−0.001520, 0.001375] | Same-region interaction failed the fixed +0.01 increment criterion |

Saved summaries inspected:

- `GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/chromatin_conditional_21849192/results.json`
- `GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/within_histology_21849320/results.json`
- Current numbers and qualifications: `docs/RESULTS.md`, “Two-model molecular-value development”; current program §23 and T15/§22.3.5.
- Existing implementation: `scripts/analysis/masld_two_models/within_histology.py` and `chromatin_development.py`.

### Why a complete-pair loss alone is redundant

For participant error vectors e_i = prediction_i − measurement_i, with
e_bar their mean, independently expanding the quadratic gives

```
sum_{i<j} ||e_i - e_j||^2
 = (1/2) sum_{i,j} (||e_i||^2 + ||e_j||^2 - 2 e_i' e_j)
 = n sum_i ||e_i||^2 - ||sum_i e_i||^2
 = n sum_i ||e_i - e_bar||^2.
```

Thus uniform complete-pair squared error is centered participant squared error
up to n. An unpenalized intercept absorbs the error mean, so this does not create
a new learning problem merely by enumerating all pairs. With a ridge penalty,
failing to scale the penalty by n changes regularization, rather than supplying
new biological information. Mean-profile loss plus complete-pair loss can change
the intercept/penalty balance, which must also be matched in the comparator.

A weighted or incomplete metadata-defined pair graph changes the estimand.
Its contrast loss is e' L e for graph Laplacian L, summed across responses.
Disconnected grade/sex strata leave component mean errors unidentified by that
loss; a profile term is needed to anchor predictions. Pair count remains a
dependent construction from participants, not an increase in biological n.

### Concrete development experiment, before any further outcome inspection

1. Reuse the existing five participant outer folds and authoritative IDs.
   Split participants first. Construct training pairs only within the outer
   training set and inner pairs only after each inner donor split. Neither
   member of an evaluated pair may enter its model training or preprocessing.
2. Use identical four recorded histology grades and sex to define the primary
   training pair graph. Compute metadata-only train/inner/test pair counts and
   represented strata before committing to the fit. Do not silently relax the
   matching if a fold has too few pairs. Use deterministic pairing or a fixed
   graph, and normalize participant contributions to prevent large strata from
   dominating solely through quadratic pair counts. Do not choose weights from
   held chromatin errors.
3. Compare existing global RRR and global-plus-local recipes with the same
   family fitted to profile loss plus weighted within-stratum contrast loss.
   Include a linear weighted/graph-regression control. Use a bounded contrast
   weight menu, for example 0/0.1/1 after explicit normalization; choose it only
   within the outer training set. Match rank, penalty search and preprocessing
   budgets across arms. All target scaling, nuisance regressions, RNA means,
   components and local-gene choices must use permitted training participants.
4. Keep the existing metadata-defined 22 disjoint evaluation pairs unchanged,
   and report full 99-person profile loss alongside contrast loss. The pair
   subset is low-fibrosis development evidence and cannot become a new blind
   confirmation. Report the same reference denominator for all arms, paired
   effect differences, and a separately defined Holm comparison family. Pair
   bootstrap describes the existing disjoint pair subset conditional on fits;
   uncertainty from model refitting remains separately unresolved.
5. Require a gain over the strongest same-input comparator, and disclose the
   training/evaluation strata represented. A gain disappearing against weighted
   linear regression supports a useful simple model, not architecture novelty.
   Do not relax the current +0.01 original-target complexity criterion for a
   complex extension after observing its result. A different donor-difference
   usefulness margin must be fixed on its own scale before execution.

**Useful next authorized action:** metadata-only census of within-stratum pair
availability at each outer and inner split plus a written loss/comparator
contract. That precedes compute. Source-only normal-liver admission can proceed
in parallel. No new evaluation outcomes are authorized by this report.

### Exact development paths and biological-unit checks

- `Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture/molecular/participant_axis.tsv`
- `Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture/folds/participant_outer_folds.tsv`
- `Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv`
- `Analysis/MASLD_Model_Benchmark/executions/chromatin-stable-rrr-forms-20260908T192024Z/out/technical_features.tsv`
- `Analysis/MASLD_Model_Benchmark/config/datasets/gse267145_znf469_human_liver.toml`

The configuration records triangulation of GEO, the source article, deposited
titles, source histology and both matrix axes, yielding 99 exact pairs. The
producer checks 99 unique IDs and fold agreement. These records support the
existing exact-identifier join; this review did not independently verify donor
identity genetically. The configuration's generic `blocked` state has narrowly
authorized development exceptions and should not be interpreted as permission
for arbitrary source tasks. Its current coordinate rule permits declared
`counted_interval` sequence semantics with a 1-bp residual; called-peak boundary
correspondence remains unresolved. Fractional RNA estimates must not be rounded
or passed to an integer count likelihood.

### Independent normal-liver source: metadata only, not a scored validation

The existing admission JSON records 17 exact identifiers in both processed
GSE128072 RNA and H3K27ac matrices and their corresponding GEO assay metadata.
The IDs are B1, B10, B12, B15, B19, B21, B22, B24, B26, B36, B38, B41, B46,
B47, B7, B8 and B9. It reports 49 RNA columns and 18 H3K27ac columns before
intersection. This is a normal deceased-donor liver source; histology is absent.

Exact metadata/input locations, relative to the repository:

- `GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/source_metadata/admission.json`
- `GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/source_metadata/GSE128072_family.soft.gz`
- `GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/source_metadata/GSE128072_Pr1_TPM.txt.gz`
- `GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/source_metadata/GSE128072_Pr1_H3K27ac_FeatureCounts.txt.gz`
- Admission producer: `scripts/analysis/masld_two_models/source_admission.py`.

The producer joins exact donor strings across SOFT assay characteristics and
matrix headers and rejects duplicate column IDs. This review read its logic
and saved metadata report; it did not read matrix values or score the source.
`admission.json` explicitly records **no performance measured**. That supports
unspent evaluation by this candidate workflow, not a claim that the source was
never inspected elsewhere or absent from a foundation model's training.
Source-family overlap, checkpoint exposure and genetic donor independence have
not been established here.

Before scoring: resolve TPM versus training fractional-estimate input units;
lift hg19 intervals into GRCh38; define outcome-blind compatible target geometry
and account for failed/ambiguous/multiple mappings; determine whether existing
aggregate feature counts can validly measure the selected regions. Mere overlap
cannot recreate counts on different boundaries. H3K27ac ChIP and the development
CUT&RUN measurements also require a stated assay transport contract. Freeze
transformations using training/reference information, not this source's measured
prediction errors. No source-tuned calibration can be reported as untouched
transport. Normal-liver transport would not validate within-MASLD histology
contrasts or advanced-fibrosis performance. With 17 participants, prediction
uncertainty is limited regardless of the number of regions.

## 2. Predict a defined variant–regulatory-region association

**Research use.** Predict signed liver accessibility association for a nominated
variant and a defined region, in source-native ALT-dosage beta units. The
target-conditioned contract differs from the current variant-centered adapter,
which does not consume arbitrary peak coordinates. This could support selecting
regulatory-region experiments. Association prediction is not isolated editing
effect prediction, a gene link or a MASLD causal-effect estimate.

**Existing evidence.** The saved development comparison has 24,430 matched
Currin leads across 17 chromosomes, representing 94.88% of 25,748 eligible
development leads. MSE is 0.300161 for calibrated native 1 Mb, 0.302252 for the
five-seed 2-kb adapter and 0.290705 for their combination. The combination's
3.15049% relative improvement misses the fixed 5% complexity margin. Currin
fold 0 was spent and remains excluded. Adapter selection was non-nested across
folds 1–4, so these inspected values are development evidence. Variant rows and
chromosome groups are not independent human donor replications.

Evidence inspected:

- `GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/allele_combination_21849200/results.json`
- `scripts/analysis/masld_two_models/design.json` and `README.md`.
- `variant_effect_query.py` is an exact saved-pair query; its stated arbitrary
  target refusal correctly limits current capability.

**Concrete next experiment.** Before another adapter search, define a
target-conditioned REF/ALT head with source peak geometry fixed, and compare
it with native 1-Mb track effects integrated over the *same* valid target and
the strongest current calibrated simple predictor. Preserve chromosome-held
splits and exact matched association rows. Fit calibration and geometry choices
inside training folds. Check allele-swap antisymmetry and identical-allele zero
effects. If eligible all-tested association labels become admitted, evaluate
the expanded population separately from the selected-lead intersection, with
native effect/SE fields and capped precision weighting. Do not invent zero
labels for unmeasured associations or use incompatible LD to reinterpret beta
as a causal sequence effect. The target-conditioned branch needs its own fixed
development usefulness margin and budget before new outcomes are inspected;
the existing 5% margin cannot be retroactively relaxed.

**Independent source limitations.** Wenz/GSE277774 is presently unadmitted.
`source_metadata/wenz_admission.json` and
`source_metadata/wenz_donor_provenance_update.json` record 14,076 significant
peak leads with RASQUAL Pi, not Currin FastQTL beta or the all-tested denominator.
They also record unresolved final-175 donor membership, per-row allele/build
harmonization and checkpoint exposure. Of reported REF alleles, 13,688 match
GRCh38, 238 match ALT and 150 neither. Collection chronology supports donor
non-overlap only as [Inference]; it does not establish an exact crosswalk.
Use native Pi-compatible allele-aligned rank/direction only after admission,
not Currin-beta MSE. No external performance was measured in this review.

**Useful next authorized action:** metadata-only target/coordinate and external
source admission inspection, followed by a fully specified matched development
comparison. Keep this scientific track alive while the chromatin experiment
develops; do not repeat the old combination merely because it already runs.

## Scientific novelty and verification limits

These are distinct hypotheses and potential usable capabilities, not proven new
methods. Current program T01/T05/T08/T15/T20 and §22 already motivate them.
Novelty requires comparison with prior methods and compatible independently
measured endpoints; this report makes no external first-method claim. No
literature search was needed to assert a novelty result because none is asserted.
The existing exposed-development results were checked in saved JSONs and source
code; no numerical rerun or independent donor-identity verification was executed.
The initial concatenated bootstrap read exceeded tool output limits, so this
agent inspected the navigation, binding scope, current relevant numbers,
model/task/source eligibility sections and current executable records; it does
not claim a complete line-by-line review of every large authority document.
Neither model may be described as independently validated from these records.

## GSE128072 metadata-only follow-up

The donor-contrast proposal remains a **candidate hypothesis**. It is not chosen
as a new model before the weighted-linear baseline, pair availability and
matched estimand are assessed.

### Processed-file headers and public methods

Only matrix header lines were read in this follow-up. The RNA header has 49
sample IDs, B1–B50 excluding B14. It has no `expected_count`, `length` or
`effective_length` field. The H3K27ac header has the six annotation columns
`Geneid, Chr, Start, End, Strand, Length` and 18 sample columns carrying
`-H3K27ac` suffixes. The H3K27ac `Length` column is a chromatin-feature length;
it is not an RNA gene effective length. Matrix measurement values were not
read or scored.

Local GEO SOFT metadata declare hg19 alignment, GENCODE v19 RNA annotation,
STAR two-pass alignment and RSEM gene-level TPM. They declare H3K27ac ChIP
alignment with BWA-MEM, mapping-quality filtering, pooled ChIP/input peak
calling with MACS2 and feature-count supplementary output. The registered
supplementary files are the two cohort TPM tables and the two histone-mark
feature-count tables; this inspected listing supplies no RSEM count/effective-
length table. Absence from that listing does not establish that such outputs
cannot be obtained elsewhere.

The original study describes RNA and ChIP as portions of the same donor wedge
biopsy. Its RNA methods use RSEM TPM and expected counts for an expression
filter; the QTL workflow subsequently transforms expression across individuals.
Those QTL transforms must not be assumed to be the content of the deposited
TPM matrix from its filename alone. Direct PMC opening encountered reCAPTCHA;
the indexed primary-method text and local GEO processing metadata were
available. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6612522/),
[GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128072).

The exact donor-string intersection is a supported assay linkage, with 17
participants after B14 is absent from RNA. It is not an independent genetic
identity check, a histologic MASLD/control classification or a source-exposure
audit. The term normal-liver transport is the existing candidate designation;
the inspected header/processing metadata do not establish histologically normal
status for each paired donor. No histology-adjusted transport test is feasible
from the currently inspected metadata.

### What RNA transformation is feasible

RSEM distinguishes expected counts, TPM and effective length. Its gene effective
length is an isoform-abundance-weighted average, so a fixed genomic exon span
is not the missing sample-specific effective length. The output format can
supply the needed fields together. [RSEM output documentation](https://deweylab.github.io/RSEM/rsem-calculate-expression.html#output).

**Derived conditional route:** if the source's actual matching effective lengths
L_gs and complete TPM denominator are recovered, then under the source's
consistent RSEM definition relative expected counts are proportional to
TPM_gs × L_gs. Library-scale invariant CPM can be formed as

```
CPM_gs = 1e6 * TPM_gs * L_gs / sum_h(TPM_hs * L_hs).
```

The absolute count total cancels for that normalized input. This algebra does
not recover read-count sampling precision, and it does not eliminate RSEM versus
the training quantifier's measurement differences. Source-native expected counts
are preferable if available. Neither exact effective lengths nor expected counts
are present in the inspected RNA header/listing, so this route is **not currently
executable from those files alone**. Multiplying TPM by arbitrary exon length
would be an approximation with a new input contract, not validation of the
released fractional-count scorer.

A fixed TPM-specific RNA model could instead be trained and frozen on a
compatible development representation, but matching its measurement definition
would require the analogous training-side lengths/quantification or a separately
justified transformation. Applying batch quantile/rank mapping to the 17 people
would be cohort-adapted input transport. It is not individual fixed-route
validation and cannot be chosen by this source's accuracy. Requantifying source
raw RNA and chromatin could address the missing information, but raw acquisition
and compute were not performed or authorized by this metadata-only assignment.

### What chromatin transformation is feasible

The header supplies hg19 intervals and feature lengths; GEO establishes the
assembly. A metadata-only interval mapping can determine unique liftovers and
whether a frozen GSE267145 target has compatible boundaries, without reading
participant feature counts. That is the next useful admission action after an
existing chain/tool is verified. A split or nonmatching interval does not
provide counts on a new boundary. Summing partially overlapping feature counts
would double count or redistribute unknown read coverage. Exact or explicitly
admitted region geometries are needed; otherwise raw recount is a separate
required task. ChIP versus CUT&RUN normalization and signal concentration also
need a prespecified transport definition before score comparison.

**Disposition:** source linkage is promising and partially verified, but exact
deployed RNA→H3K27ac validation remains unadmitted. Remaining inputs are source-
compatible RNA expected counts/effective lengths or a fixed justified alternative,
an outcome-blind interval correspondence, an assay-transport target definition,
source-family/checkpoint exposure and an appropriate participant-level uncertainty
contract. No matrix values, validation predictions or measured prediction errors
were inspected in this follow-up.

## Independent derivation of the proposed exact-grade graph objective

This section checks mathematics before implementation. No outcome array, model
fit or numerical job was used. The proposed objective is an adaptation to a
specific participant-difference target; usefulness and novelty remain unproven.

### Graph definition and exact loss

Let n be the number of training participants, p the number of training
participants in metadata strata of size at least two, and m the size of one
such stratum. Root proposes a block-diagonal graph matrix with

```
P_m = I_m - 1_m 1_m' / m
L_m = m/(m-1) * P_m                         for m >= 2
L_m = 0                                   for singleton strata
W   = I_n + lambda * (n/p) * L             for p > 0
lambda in {0, 0.1, 1}
```

This is a valid normalized complete graph within each eligible grade/sex
stratum. Each paired participant has weighted degree one: diagonal L is one,
off-diagonal entries within the same group are −1/(m−1), and all other entries
are zero. It follows that L1=0 and trace(L)=p. Its block loss is exactly

```
tr(E_m' L_m E_m)
 = m/(m-1) * sum_i ||e_i - e_bar_m||^2
 = 1/(m-1) * sum_{i<j} ||e_i - e_j||^2.
```

Thus the proposed fit minimizes

```
sum_i ||e_i||^2
 + lambda * (n/p) * sum_eligible_strata
       [ sum_{i<j} ||e_i-e_j||^2 / (m-1) ]
 + alpha * ||B||_F^2.
```

The n/p factor puts the graph trace on the same n scale as the profile term:
trace(W)=n(1+lambda). It gives equal *total incident edge weight* per paired
participant, not equal weight per stratum. A larger stratum contributes more
total weight through its larger number of participants, as intended by donor
weighting. It does not create p independent pair observations. The primary
scientific units are the n participants; their molecular response columns and
their possible pairs remain correlated measurements.

### Rank, intercept and invariances

L_m has eigenvalue zero on its constant vector and eigenvalue m/(m−1) on its
m−1 contrast directions. If there are k eligible strata, rank(L)=p−k. W has
eigenvalue one on all stratum mean directions and singleton directions, and

```
1 + lambda * (n/p) * m/(m-1)
```

on each eligible contrast direction. Hence W is symmetric positive definite
for the proposed nonnegative lambda, and its condition number is bounded by
1+2 lambda n/p. No numerical jitter is needed to make this graph positive
definite. The p=0 case is undefined as written and requires an explicit
disposition: no contrast fit is estimable from that training split; retaining
W=I can identify the profile-only control but must not be reported as graph
training. Pair availability in every inner split matters, not just the full
outer training roster.

Because W1=1, the unpenalized intercept uses the ordinary training participant
mean: (1'WY)/(1'W1)=mean(Y), and similarly for X. The matrix square root leaves
the constant vector unchanged. Metadata-row permutations with corresponding
X/Y permutations leave the objective unchanged. Reordering the endpoints of
any pair or reversing all contrast signs leaves it unchanged. Adding a common
target/prediction offset within a stratum leaves graph loss unchanged but can
change profile loss. This is why the profile term anchors the group means.

Response units must be matched across arms. Multiplying a response column by a
constant changes its influence on a joint low-rank fit unless a corresponding
training-only scaling or explicit scientific response weight is applied. That
is a modeling choice, not invariance of the unscaled multivariate objective.

### Weighted ridge and RRR baseline

For a centered linear model Y=XB+E, the full-rank weighted ridge solution is

```
G = X' W X + alpha I
B_ridge = G^(-1) X' W Y.
```

An ordinary ridge solver applied to W^(1/2)X and W^(1/2)Y is an equivalent
control, with the same penalty alpha. Query predictions still use the original
appropriately centered/scaled X; W is a training loss matrix and is not a query
feature. All transformations and graph construction use that training split.

For alpha>0, a rank-r coefficient constraint minimizes the same quadratic by
forming C=G^(−1/2)X'WY, taking its best rank-r SVD approximation C_r, and setting
B_r=G^(−1/2)C_r. This follows by completing the square in D=G^(1/2)B. Truncating
the ordinary coefficient matrix or the unweighted fitted predictions is not
generally the minimizer of this ridge-weighted rank constraint. The singular
alpha=0 case needs a pseudoinverse/support treatment rather than assuming G is
invertible. An implementation may use an algebraically equivalent dual route
to avoid a genes-by-genes matrix.

The necessary scientific comparator is the **same weighted linear/ridge or
weighted RRR family**, alongside the existing unweighted global-plus-local
recipe, with equal rank/penalty tuning, training participants and target units.
Graph weighting may deliver a useful adaptation without new architecture.
An extra encoder or nonlinear interaction cannot take credit for an objective
change that a simple weighted fit already supplies.

Penalty scaling is consequential: multiplying all of W by a scalar preserves
an unpenalized solution, but preserves penalized solutions only if alpha is
multiplied by the same scalar. Comparing trace-normalized W with unnormalized
W at a fixed numerical alpha changes regularization. Fix a consistent penalty
normalization or tune it within training for every arm; do not interpret that
accidental change as a contrast-learning benefit.

### Counterexamples and scientific limits

1. **One complete stratum, no singletons.** Then p=n and
   L=n/(n−1)P_n. After an unpenalized intercept is fitted, all centered residuals
   lie in the contrast space. The objective equals ordinary centered SSE times
   1+lambda n/(n−1), apart from the penalty. With matching penalty scaling,
   graph training is exactly redundant. This recovers the complete-pair identity
   above and is a useful implementation invariant.
2. **Equal-size strata plus fitted unpenalized stratum intercepts.** Residual
   means vanish in every stratum. With all participants paired and all m equal,
   W acts as the same scalar on every remaining residual direction; again the
   only difference at fixed alpha is regularization. With unequal m, the
   remaining change is explicit stratum-dependent residual weighting. It should
   be compared with that weighted-linear control.
3. **Few paired participants.** Suppose n=80 but only two participants share
   a stratum. For lambda=1, the sole within-pair contrast eigenvalue is 81,
   while every other eigenvalue is one. Trace normalization does not prevent a
   single measured difference from exerting large influence. Even lambda=0.1
   gives that contrast eigenvalue nine. A high leverage noisy pair could then
   drive the fit. The metadata census must expose this situation before fitting;
   do not treat n/p normalization as a guarantee of adequate effective evidence.
4. **Singleton grade/sex strata.** Those participants supply profile loss but
   no within-stratum difference supervision. Good whole-profile error cannot
   demonstrate donor-difference prediction in those strata. Missing F2–F4
   evaluation pairs remain missing scientific coverage.

W is a **chosen loss weighting**, not a measured inverse residual covariance.
Calling this a generalized least-squares optimization is mathematically valid,
but it does not justify GLS standard errors, independent edge sampling or a
claim that the group contrasts have the assumed noise precision. Shared donor
errors induce covariance among pair differences. Uncertainty must use donors
(or the existing disjoint pair units where that endpoint is reported), retain
dependence from shared training folds and distinguish fixed-prediction bootstrap
from refitting uncertainty. Metadata-defined grade matching provides the stated
cross-sectional target; it does not establish conditional independence from
unmeasured pathology or composition.

**Disposition:** the proposed W is algebraically correct and supports a bounded
weighted-linear/RRR development comparison once metadata availability and
penalty conventions are fixed. It can be scientifically distinct from the
existing unweighted RRR through unequal emphasis on within-grade residual
contrasts versus between-grade/profile errors. Whether that distinction produces
a useful generalizable model is not established by the algebra. No model has
been selected or validated by this derivation.

## Read-only implementation review before the development fit

Reviewed `scripts/analysis/in_biopsy_allelic/rec_codex/weighted_chromatin.py`,
`chromatin_contrast_development.py`, `chromatin_contrast_design.json`, and the
imported `scripts/analysis/masld_two_models/chromatin_acquisition.py`. No outcomes
were read and no job was run by this reviewer. Root reports the independent
synthetic primal comparison job 21995378 passed to at most 6.2e−15; this review
does not substitute that numerical check for biological model validation.

**Algebra and preprocessing findings:** graph construction matches the derived
degree/trace normalization exactly. Cholesky orientation is correct. The dual
rank fit uses sqrt(d/(d+penalty)) times projected weighted responses to determine
right singular vectors; their Gram matrix equals the primal penalized optimum's
response Gram matrix. Projecting ridge query predictions through those vectors
therefore implements the correct rank constraint. Penalty normalization is
relative to the weighted kernel's largest eigenvalue and is consistently applied
across strength arms. Feature means/scales are fitted inside each training call.
Concentration-descriptor regression is re-estimated inside every inner split
and outer training split. Per-library log-CPM and concentration descriptors do
not borrow participant-population parameters. Computing the descriptor of a held
chromatin measurement is part of the residual *evaluation target*, not an input
available to predict absolute H3K27ac. No inner preprocessing leakage was found
in the inspected code.

### Fix or settle before the biological development fit

1. **Establish the strongest comparator and interpretation.** The four current
   arms compare panel-only ridge/rank-12 fits. The prior selected RNA model fits
   a global response representation over all 96,460 regions and adds a local
   gene correction. That stronger comparator is absent from the current trial;
   fitting only the 1,000 responses also changes shared representation learning.
   Include the strongest global-plus-local comparator on the same target panel
   with matched preprocessing, or explicitly settle this as an initial bounded
   objective-feasibility experiment whose positive result requires that further
   comparison before a model upgrade. Existing archived OOF predictions may be
   a useful descriptive reference only after verifying exact target transform,
   participant/panel axes and prior exposure. Beating `best_unweighted` as
   presently defined does not show improvement over the selected current model.
2. **Verify exact fixed pairs.** Development lines 134–136 check only 22 pairs
   and 44 distinct IDs. That does not prove the pairs are the *unchanged* prior
   metadata-defined set. Compare canonical unordered participant-ID pairs with
   `within_histology_21849320/pairs_fixed_from_metadata.tsv` under the exposed
   candidate root before reading measurement arrays. Same seed/IDs/grouping
   make equality plausible, but cardinality alone cannot support that claim.
3. **Validate all five fold assignments and coverage.** Lines 70–74 check
   agreement with the join, but later loops hard-code 0–4. Assert the exact fold
   set is {0,1,2,3,4}, every donor receives one held-fold prediction, and no pair
   joins different folds. Otherwise an unexpected matching fold label can leave
   zero-initialized truth/predictions outside the loop while the pair count still
   agrees. This is a structural metadata check, not a biological rerun.
4. **Correct the estimator description.** The JSON's phrase “identical objective
   is not a distinct estimator” for weighted ridge versus weighted rank-12 is
   imprecise. They share a scientific target and weighting, but the rank
   constraint makes them different estimators. Neither difference by itself
   establishes method novelty. Their direct comparison can remain descriptive.
5. **Require defined skill and finite outputs.** Before division at lines
   149–155, assert finite positive measured-contrast reference denominator and
   finite resampling denominators/summary values. A zero reference contrast
   makes skill undefined; do not print a margin decision or percentile interval
   from NaN. This is needed for the stated numerical result contract.

The current profile safeguard selects against the minimum inner profile SSE
over both unweighted arms and all their penalty choices, which is stricter than
comparison with the pair-selected unweighted winner. Its design text should
retain that exact reference. A safeguard refusal currently raises an exception
for the whole run, while the JSON says it stops the affected arm; decide that
disposition explicitly. The outer profile SSE still needs to be reported beside
pair skill: an inner safeguard is not a guarantee of held-participant profile
preservation. Conditional pair-bootstrap intervals remain development summaries
and do not establish independently replicated precision or generalization.

## Target-peak native readout: independent installed-source check

Read-only source inspection, 2026-09-30. No outcome table, molecular effect,
protected fold-0 prediction, inference, API call, download, or compute job was
read/run for this check.

The exact geometry comparison is **the same signed native ATAC reduction on
the associated full peak versus the existing variant-centered 501-bp mask**.
For allele a, track t, and genomic interval R define
`S[a,t,R] = sum_{b in R} returned_ATAC[a,b,t]`. The native scalar is
`x[R] = (1/3) sum_t {log2(1+S[ALT,t,R]) - log2(1+S[REF,t,R])}`.
The three tracks are the exact names `UBERON:0001114 ATAC-seq`,
`UBERON:0001115 ATAC-seq`, and `UBERON:0002107 ATAC-seq`. Existing
`R501 = [position1-1-250, position1-1+251)`; the proposed associated-peak
readout uses `Rpeak = [peak_start_hg38, peak_stop_hg38)` on the same chromosome.
This is `DIFF_LOG2_SUM`: sum before log, one pseudocount per summed track,
then average track effects. `DIFF_SUM_LOG2` (sum of per-base logs), unsigned
L2, a fixed peak-center 501-bp mask, log of the track-average sum, or averaging
per-base values would change the proposed comparison. A peak-width change
also changes how consequential the pseudocount is at low signal; report
width/coverage and avoid interpreting a gain automatically as distal biology.

Installed research source root (resolved from the environment's editable
package `.pth` and `direct_url.json`) is
`Analysis/MASLD_Model_Benchmark/executions/alphagenome-runtime-probe-20260915T104500Z/code/alphagenome_research/src/alphagenome_research/`.
Exact implementations:

- `model/variant_scoring/center_mask.py:33` creates the center mask;
  `_apply_aggregation` at line 79 implements the signed reduction, with
  `DIFF_LOG2_SUM` at line 114.
- `model/heads.py:299` implements `predictions_scaling`; `GenomeTracksHead`
  unscales output at line 598, applies softplus at line 653, and stores both
  scaled and unscaled predictions at lines 664–670. ATAC head configuration
  at lines 141–149 has `apply_squashing=False`.
- `model/dna_model.py:124` selects ATAC `predictions_1bp`, rather than
  `scaled_predictions_1bp`; `_construct_output_from_predictions` at line
  1303 writes those values directly into `TrackData`.
- Project executable definition is
  `scripts/analysis/masld_two_models/variant_inference.py:74` (`NativePredictor`),
  importing masks/reduction via
  `scripts/analysis/alphagenome_program/i1_common.py`.

Thus the **returned ATAC prediction is not already log transformed**. The
head produces positive values, reverses its soft clipping and multiplies by
track `nonzero_mean` and resolution to restore its experimental target scale.
For ATAC it does not apply the RNA-only power transformation. The readout
subsequently applies `log2(1+sum)`. No Atlas quantile calibration is loaded in
this local runtime. These decoded track values still are predictions in the
model's training-target scale; they are not Currin's participant-normalized
FastQTL phenotypes. [Unverified] The precise deposited Currin phenotype
transformation was not independently rederived in this source-only check;
the inspected endpoint code identifies its coefficient as source-normalized
accessibility per additional ALT allele, not a log2 expression change.

To compare beta prediction, separately fit each readout's zero-intercept
training slope `s[R] = sum_train x[R]*beta_ALT / sum_train x[R]^2` and use
`beta_hat[R] = s[R]*x[R]`. Reject zero/nonfinite denominator. Preserve identical
training/evaluation rows and development folds, excluding fold 0. Report
paired MSE in beta squared, RMSE in beta units, signed association and coverage;
raw native log scores cannot be compared numerically with beta by RMSE without
calibration. Separate slopes make this a comparison of target geometry with
equal calibration capacity. A negative fitted slope should be reported rather
than silently reorienting alleles. Calibration is supervised development;
target intervals are source association metadata and do not confer untouched
external evaluation status.

Orientation and coordinate contract is explicit in
`scripts/analysis/alphagenome_program/c2_01_build_variants.py:2`: source
`chr:position1:REF:ALT` beta is unchanged ALT-dosage beta. `EA` identifies the
accessibility-increasing allele; negative beta normally has EA=REF. Flipping
beta to EA erases its ALT orientation. That builder validates EA/NEA identities,
GRCh38 FASTA REF, complete ACGT windows and chromosome folds, but retains one
peak per variant selected by historical q-value deduplication. For exact
coordinates, `scripts/analysis/alphagenome_campaign/data_endpoints.py:182`
joins the independent deposited BED peaks to explicit hg38 start/stop fields
and requires equality, calling the resulting contract GRCh38 BED0 half-open.
This check inspected the verification code, not executed it on outcome files.

Before a proposed implementation:

1. Require chromosome equality, integer `start < end`, and complete peak
   containment in the **returned** 1-bp ATAC interval. Build a region mask from
   those coordinates; `create_center_mask` is not a full associated-peak mask.
   Do not silently clip or treat an empty/partial peak as a valid full-peak
   prediction. Check the output interval and value-array length as well as the
   requested sequence interval. A distal target outside 1 Mb is unsupported
   by this contrast; increasing or moving context would be a separate change.
2. Keep the 1,048,576-bp variant-centered input and SNV ALT substitution exactly
   unchanged. Use the same GRCh38 FASTA, verify REF at `position1-1`, reject
   REF=ALT/non-ACGT/edge overhangs, and do not extrapolate to indels whose allele
   coordinate alignment requires separate handling. These guards exist in
   `variant_inference.py` and `i1_common.py`.
3. Calculate **both** masks from the same freshly generated REF/ALT tensors.
   Existing documented between-process variation means an archived center
   score versus a fresh peak score confounds readout geometry with numerical
   variation. Preserve the default requested ATAC+DNASE outputs, exact named
   track identities/order, checkpoint, backend and hardware; dropping DNASE
   changes a static compilation argument and requires a reproduction check.
   A same-process repeat measures deterministic execution, not fresh-process
   reproducibility. Record both. No new reproduction was run here.
4. Freeze eligibility independently of beta/sign/performance. Keep source
   variant–peak identity, not variant alone, if admitting multiple associated
   peaks later; repeated variants/nearby windows remain dependent and must
   share training partitions. Historical selected leads are a bounded
   development population. Do not rerun `model_baseline.evaluate`, whose old
   loop includes fold 0, as authorization for a new comparison.

The useful next experiment is a matched development-only center-versus-full-
peak readout, with shared inference tensors and separate training calibrations.
This diagnoses a target mismatch directly; it does not establish an adapter
architecture improvement or independent new-variant generalization.

## Inductive missing-program contribution: normalization check

Read-only code/spec review; coordinated with `program_hypotheses`. No RNA
expression, clinical outcome, fit, or job was read/run. The current proposed
GSE268273 metadata design and exact source axes are in
`2026-09-30_codex_program_model_options.md`, section “Prespecified
missing-contribution development design”; its eight admitted programs and
714-panel-name overlap are the other agent's metadata inspection, not an
independent numerical recount here.

**An inductive missing-component target is definable, but the deployed
observed-plus-predicted full score requires a separate exposure contract.**
Training-only centering avoids one leakage source. It does not make a
complete-reference library denominator available to a targeted measurement.

Exact existing normalization/scoring routes inspected:

- `scripts/analysis/masld_two_models/chromatin_acquisition.py:32` computes
  `log2(1+1e6*raw/raw.sum(axis=1))`. Its denominator includes the delivered
  complete gene axis, including genes hidden by a later mask.
- `scripts/manuscript/program_observability_map/analysis_lib.R:259`
  implements the unchanged 117-program weighted projection. It aggregates
  mapped-symbol weights, z scores each gene **within receiving cohort**,
  selects genes variable in that cohort, renormalizes original weights over
  retained genes, and finally standardizes each program score within cohort
  (`standardize_scores_by_cohort`, line 136). These are retrospective
  association units. The function cannot be reused unchanged for a single
  new participant or a training-frozen inductive target.
- `scripts/analysis/histology_anchored_continuum/04_score_programs_and_validate.R`
  similarly calls `edgeR::cpm(dge, log=TRUE, prior.count=1)`, then cohort gene
  z scores and score standardization, with its prespecified signature
  exclusion. This is a different score definition and should not supply
  missing-component ground truth for unchanged complete memberships.
- `scripts/analysis/cross_assay_modules/06_score_external_bulk.R` uses
  cohort gene z, an equal-weight observed-member mean, and receiving-cohort
  score standardization. That separate bulk-module family has no native
  canonical weights; it must not be substituted for the 117 Hotspot programs.

The expected-count fixture is described in the other agent's inspected
receipt as fractional raw RSEM expected counts. A continuous regression after
an explicit count-to-relative-expression transform is compatible with those
units. Do not round expected counts, apply an integer-count sampling likelihood
without justification, or call these values TPM. The fixture's delivered-axis
sum is not independently verified to be the whole transcriptome/library.
Avoid “full transcriptome CPM”; say CPM relative to the admitted reference
gene axis and preserve that exact axis for all receiving targets. Duplicate
symbols must follow one declared count aggregation before logging, distinct
from existing scripts that average log-expression rows.

### Contract A: valid masked bulk-RNA missing component

Freeze gene universe U, actual admitted panel O, mapped program genes Gp and
original weights `w_pg`. For a complete bulk-RNA evaluation participant define
`Lfull_ig = log2(1+1e6*C_ig/sum_(h in U) C_ih)` and
`Zfull_ig=(Lfull_ig-mu_g)/sigma_g`, with gene transforms fitted on training
people only. The supervised target is
`Mfull_ip=sum_(g in Gp\O) w_pg*Zfull_ig`.
Predictors are `Xpanel_ig=log2(1+1e6*C_ig/sum_(h in O) C_ih)` for g in O,
then training-fitted standardization. Computing complete-reference totals
for **held target construction only** is valid; using those totals as feature
offsets, eligibility filters, PCs, or adaptation inputs leaks hidden genes.
The estimand is `E[Mfull_p | Xpanel]`. No test-cohort gene/score recentering or
test-dependent constant-gene exclusion is permitted.

The measured panel contribution in these predictor units is not
`Sobs_full=sum_(g in Gp∩O) w_pg*Zfull_g`: its denominator differs. Hence
`Sobs_panel + Mhat_full` does not estimate the stated complete-reference score
by the original additive decomposition. Full-score error using an oracle
complete-reference denominator can be reported as a computational diagnostic,
but it is not a deployed targeted-assay endpoint. The program agent's current
design explicitly makes this limitation; I agree with that bounded target.

### Contract B: additive observed-only-reference target

If an actual additive score is needed, fix a denominator reference D contained
in O and available at deployment, with exposure rule
`q_i=sum_(h in D) C_ih` (or a separately measured, assay-compatible exposure).
Apply the **same** `L^D_ig=log2(1+1e6*C_ig/q_i)` to all genes of the complete
training/evaluation reference and freeze training gene means/SDs. Define
`S^D_p=sum_(g in Gp) w_pg*Z^D_g`,
`M^D_p=sum_(g in Gp\O) w_pg*Z^D_g`, and
`Sobs^D_p=sum_(g in Gp∩O) w_pg*Z^D_g`.
Then `Sobs^D_p + E[M^D_p|Xpanel]` is an inductive additive estimator.
Refuse zero exposure or missing reference genes. This defines **new reference-
normalized units**, preserves memberships/weights, and requires separate
prespecification; it cannot be renamed the existing cohort-standardized
Hotspot or full-reference-CPM score. An observed-panel denominator can create
shared composition covariance, so successful prediction alone does not imply
regulatory recovery. Actual CosMx counts still differ from masked bulk RNA in
detection, aggregation and tissue composition; no receiving full-RNA target
on the same spatial people was established.

### Strong simple comparators and preprocessing boundary

The strongest baseline is unconstrained multivariate ridge on all admissible
panel genes, with an intercept and training-only penalty selection, predicting
the same missing component. Lower-rank ridge changes only a constraint; it
needs to outperform that comparator to support the proposed efficiency claim.
Also lock these comparisons:

1. Training conditional mean of M. Without deployment-available covariates
   this is its training mean, normally zero under consistent gene centering.
   A covariate-conditioned mean is admissible only for predeclared fields
   measured at deployment; adding histology creates a different task.
2. Affine regression of M on the observed original-weight panel scalar;
   include an intercept rather than relying on a receiving-cohort mean.
   If no program gene is observed this reduces to the training mean.
3. Retained-weight extrapolation, explicitly converted to an M prediction
   with a fixed weight denominator. Under Contract B and nonzero retained
   weight fraction r, the simple estimator is
   `Mhat=(1/r-1)*Sobs` when total program weight is one and weights positive.
   For signed/otherwise normalized weights specify the actual algebra;
   do not assume this formula from L1 coverage. Under Contract A using
   panel-unit Sobs is a heuristic comparator, not an identity.
4. Missing-gene ridge imputation followed by the unchanged weighted sum.
   With common design, intercept and penalty this equals direct ridge on M:
   ridge is linear in its response, so `B_Y*w = B_(Yw)`. Counting the two
   routes as independent model evidence would be redundant. A Gaussian
   conditional-mean formula with the same covariance shrinkage can likewise
   be equivalent to a linear ridge predictor; no distinct novelty follows.

For every inner split refit input centers/scales **and target gene means/SDs**
on the inner training people, construct both inner targets with that transform,
and normalize validation loss using inner-training M variance. Fitting those
transforms once on all 69 outer training people and reusing them across inner
validation leaks inner-held preprocessing into selection. After selection,
fit the final transform on all 69 fit people and keep it frozen for calibration
and test. All repeated specimens stay with their person. Masks, program
admission and mappings remain metadata-fixed; variance exclusions are fitted
in training and retained as limitations, not replaced by favorable programs.

A study-held-out missing-component comparison is the next meaningful test
after the exact measurement contract and source admission are settled. A
random internal test supports development performance under that source's
distribution; it does not establish transport or recover actual spatial
observability. No such fit or test was run here.

**Current utility verdict:** Contract A is a scientifically definable scalar
development endpoint and the current program-agent proposal respects its
predictor boundary. Its practical receiving use remains unresolved: a
prediction of an unobserved standardized score component alone does not give
the deployed full score, and no application tolerance or same-person actual-
panel/full-RNA reference has been identified. It should remain a bounded
computational masking hypothesis pending those decisions, rather than a new
reusable model claim. A metadata-only fixture/symbol admission can proceed;
fitting a new candidate should wait for the target, exposure and application
contract. This verdict does not reopen existing spatial observability calls
or change any of the 117 frozen programs.

## Model A P1: weak-null studentization versus conditional exchangeability

Read-only statistical review requested by root. No file named
`PRESPEC_A_FINAL` exists in the inspected scoped paths; root confirmed that
the final-plan section and operative v1/v1.2–v1.7 answers govern. Inspected
`/gpfs/commons/home/jameslee/.claude/plans/plan-thoroughly-to-build-imperative-wren.md:116`,
`docs/technical/agent_exchange/2026-09-29_model_A_likelihood_spec_v1.md:240`,
`2026-09-30_claude_to_codex_design_fixture_answers.md` sections 1–3 and 11,
the Codex likelihood/design/test recipes and existing REC status. No real
allelic values, sealed outcomes, restricted inputs, fixture hidden truth,
or B-MODEL hidden results were read. No scientific simulation was run.

**P1 is explicitly a studentized approximate test of the weak parameter null,
not an exact conditional-independence randomization test.** This limitation
was already raised in `2026-09-29_codex_spec_critique.md`, item 2, and accepted
in v1 section 5.1. The new review clarifies why even knowing the exact stage
law would not by itself make the currently specified weak-null P1 exact.

In `rec_codex/likelihood_v1.py:206`, the predictors are
`eta_ig = alpha_g*(1+kappa*S_i+delta'z_i) - orientation_ig*omega_S*S_i`
and `r_ig = r_g + rho_S*S_i + rho_b*b_i + rho_d*d_i`.
The fixed cohort mapping offset is also applied in the beta-binomial
observation layer. Setting kappa to zero leaves stage-dependent mean when
omega_S is nonzero and stage-dependent dispersion when rho_S is nonzero.
Het ascertainment further depends on those mean/dispersion parameters
through the truncation normalizer. Thus selected allele counts need not be
independent of stage given the prescribed allele-free conditioning variables.

Formally, write `q(s|Z,D,selected)` for the true conditional stage distribution and
`f_s(a|Z,D,selected)` for the permitted allele-count law, where D is the
fixed tag/depth/prior design. Under the weak null the observed joint law is
`q(s|Z,D,selected)*f_s(a|Z,D,selected)`. Replacing stage with an independent draw from
even the true q, while retaining counts, produces
`q(s*|Z,D,selected)*sum_s q(s|Z,D,selected)*f_s(a|Z,D,selected)`.
These differ whenever the count law depends on s. Stage/count pairs are
therefore not conditionally exchangeable under the composite kappa=0 null.
Nuisance refitting does not restore equality of these joint laws. Estimated
proportional-odds probabilities add approximation on top of this distinction.

The prescribed statistic in `rec_codex/tests_v1.py:p1` is
`abs(kappa_hat/cluster_sandwich_SE)`, with refitted nuisance coordinates and
SE for every draw; `design_v1.py:stage_draws` fits cohort-specific ordered
thresholds/common slopes and draws stages independently, not a uniform
permutation. Conditioning includes centered outcome z columns, log error,
log1p candidate median depth and log1p candidate row count. The log1p forms
are the v1.2 zero-row correction. Each draw recenters the complete participant
set. The sandwich in `likelihood_v1.py:354` clusters by people/relative groups.
The v1.7 bound-coordinate and start rules remain part of the test.

Studentization can make a weak-null statistic approximately pivotal under
appropriate regularity, identifiable pseudo-true parameter, adequate independent
clusters, and a suitable resampling law. It is not a proof that either the
observed or randomized fitted kappa has zero centering in this nonlinear,
ascertained, mixture model. In particular, randomized counts mix the original
stage-specific distributions, which need not stay inside the fitted
beta-binomial family after refitting nuisance terms. Correlation of stage
with z, sparse support, near-zero alpha, and active bounds make the premise
load-bearing. This does not prove miscalibration; it specifies the unresolved
claim requiring the planned synthetic check.

The spec's within-cohort uniform-permutation sensitivity is valid under its
stated stage/covariate-independence restriction for that part of the design.
It does not repair stage/count conditional dependence allowed by the weak
null. Agreement between the two numerical implementations is also distinct
from statistical size.

### Existing checks do not establish this calibration

`scripts/analysis/in_biopsy_allelic/bfix/make_fixture.py` fixes stage as
`j%3`, draws raw covariates independently, and plants kappa −0.08 plus
nonzero gene departures and nuisance stage terms. Its supplied stage
permutations and 50-draw refits are an arithmetic agreement fixture, not
repeated correlated-covariate null datasets. `check_fit_v1.py` labels its
two-draw smoke as execution only. Current
`2026-09-29_codex_rec_status.md` explicitly leaves scientific size/bias/coverage
pending. No completed calibration of this failure scenario was identified in
those authorized records; hidden B-MODEL outputs were not inspected.

### Smallest meaningful RNA-allelic null diagnostic, proposed only

Use one synthetic design with independent biological participants in three
cohorts, eight genes with mixed nonzero alpha values and balanced tag
orientations, and the existing fixture's cohort-specific integer het rules.
Keep genotype mixture, beta-binomial observation layer, mapping offset and
error dispersion in their RNA-allelic units. Set **kappa=0 and tau=0** so
gene departures do not silently change the fixed-gene weak-null estimand.
Choose interior dispersion/bias parameters and sufficient informative alpha
values to avoid an identifiability failure as the first diagnostic.

Generate stage from a known correctly specified proportional-odds law with
cohort thresholds and a nonzero slope on a recorded synthetic z covariate;
keep all three stage categories represented without near-deterministic
assignment. Generate candidate-tag depths, error/duplication covariates and
priors before reads. Produce genotype classes and RNA allele reads, choose
the max-depth tag under the fixed tie rule, and apply actual het ascertainment.
Fix the entire selected design during each P1 stage resampling, exactly as
production does. Do not generate counts only at already-selected rows and
claim it assessed unconditional ascertainment.

Only two paired nuisance regimes are needed initially:

1. kappa=0, omega_S=0, rho_S=0: correctly specified stage/count
   conditional-independence control.
2. kappa=0, omega_S=0.02, rho_S=0.10: weak-null counter-scenario with the
   fixture's nuisance magnitudes but zero planted proportional stage effect.
   Both are generated and refitted rather than fixed to their truth in the fit.

For each null dataset run the exact production estimated-stage-law P1 with
1,000 draws, both prescribed starts, cluster sandwich and conservative
failed-refit exceedances. An optional **same-dataset diagnostic** substitutes
the known stage probabilities **conditional on the retained design** while
holding every other computation fixed;
it helps separate stage-law estimation error from weak-null randomization
error and cannot replace the prespecified test. Record size, p-value shape,
observed/randomized z centering and variance, active bounds, failures,
category support and stage-law fit diagnostics. Stage-law misspecification
or related-participant dependence can be later separately prespecified
sensitivities if this basic case passes; no broad factorial sweep is needed.

Ascertainment is also relevant to that optional “known law” comparison.
Generated preselection `q_pre(S|Z,candidate depths)` is not automatically the
law conditional on the fixed selected-row pattern. For independently generated
candidate-gene observations, with inclusion probability C_ig(s) and observed
inclusion indicator I_ig, it is
`q_post_i(s) proportional to q_pre_i(s) * product_g
 C_ig(s)^I_ig * (1-C_ig(s))^(1-I_ig)`.
Compute and normalize those three probabilities from the generating mixture
and actual ascertainment rule if calling the comparator an oracle conditional
law. Keep zero-inclusion participants in the complete design. If only q_pre
is used, label it a known **preselection**-law sensitivity; it does not settle
selection-conditioned exchangeability. The production stage law's depth-only
summaries do not model this entire inclusion pattern. Conditioning on an
allele-derived inclusion pattern for the synthetic diagnostic is not permission
to add outcome-derived features to the production stage law.

The original final-plan text names at least 1,000 null likelihood/bootstrap
replicates and 200 with full CRT. **The accepted v1 Appendix B3 supersedes
that calibration decision:** size/coverage bands are judged on warp-speed
Monte Carlo at at least 4,000 replicates; 200 full-CRT replicates have the
one-sided failure rule that the binomial lower 95% size bound exceeds 0.05.
At size 0.05, 200 replicates give binomial SE 0.0154, whereas 4,000 give
0.00345. Do not require the 200 full CRTs themselves to certify the original
narrow [0.04,0.06] band. Report binomial uncertainty and retain the revised
fixed procedure; no adaptive extensions or test-law change are authorized by
this review. The operative contract is detailed in the follow-up below.

This is a null-calibration proposal, not an empirical power claim or a
recommendation to change Model A's test. The unresolved issue was already
acknowledged by the operative spec; the bounded contribution is its exact
weak-null derivation and the nuisance-on versus nuisance-off RNA diagnostic.

### Exact distribution diagnostic executed by root; generator prepared

Root's job 21995747 completed the eight-independent-gene exact-distribution
diagnostic in 13 seconds with an empty error log (runtime/log state reported
by root). I inspected its synthetic summary directly:
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/p1_selection_null_21995747/summary.json`.
It records independent SciPy beta-binomial agreement and normalized mixture
probabilities. Nuisance-off stage0-versus-stage2 total variation is zero;
posterior-stage probability shifts are numerical roundoff, at most
1.67e-16. Nuisance-on maximum total variation is 0.051217603938797256;
the all-selected pattern shifts a stage probability by at most
0.005059438569086638, the alternating-selected pattern by
0.00009171099474936106, and none-selected by 0.003518402482590055.
These demonstrate weak-null count dependence and selection-conditioned
stage-law change in that synthetic distribution. They do not establish P1
miscalibration, size or power; the summary explicitly records no fit or test.
Executed likelihood SHA256 is
`87f07515f714f2d19e2e832ac2c98b5efc8c1f0d2ee27d701836b806548217c3`,
before root's subsequent memory-path change; retain that execution identity.

Added two Codex-owned files, without modifying the production likelihood,
design or test law:

- `scripts/analysis/in_biopsy_allelic/rec_codex/synthetic_p1_preascertainment.py`
  exposes reusable `generate(regime, seed, replicate, participants_per_cohort)`.
  Default design has three cohorts of 60 independent people, eight genes,
  two allele-coding-reversed tags per gene, ordinal stage correlated with
  recorded age/sex, and the seven named RNA-model covariates. It generates
  raw depths for every tag, a shared oriented genotype per person/gene,
  and beta-binomial counts for every candidate tag before choosing maximum
  depth (ties to the lower synthetic position) and applying the exact integer
  het predicate. Six metadata-fixed people have zero candidates; all
  zero-selected people remain in the returned complete participant design.
  Independent streams share the pre-read design, genotypes and count uniforms
  between nuisance regimes. Synthetic tags represent the same haplotypes
  with reversed coding, not real LD or read co-coverage.
- `scripts/analysis/in_biopsy_allelic/rec_codex/check_synthetic_p1_preascertainment.py`
  is a compute-allocation-only verification driver. It checks count bounds,
  the correctly oriented three-class priors, complete cohort centering,
  conditional stage probability variation, depth/tie tag choice, exact
  selection, retention of zero-candidate/zero-selected people, paired
  pre-read designs and deterministic repetition. It saves versions and source
  hashes without fitting or running P1.

The generator uses **realized complete-participant cohort centering** for
read generation, matching `Design`. This couples selection patterns across
participants through the realized stage mean. Returned preselection ordinal
probabilities therefore are explicitly not an oracle selection-conditional
stage law; the earlier fixed-center eight-gene Bayes calculation must not be
copied unchanged as that oracle. Primary future calibration should use the
production fitted stage law. No full CRT or size simulation is authorized
by creating this generator.

Static verification passed:
`python3 -m py_compile scripts/analysis/in_biopsy_allelic/rec_codex/synthetic_p1_preascertainment.py scripts/analysis/in_biopsy_allelic/rec_codex/check_synthetic_p1_preascertainment.py`.
At the initial source-only review, generator numerical verification was
**not run**. Its compute command,
from the Codex directory with `CODEX_REC_OUTPUT` set, is
`PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python check_synthetic_p1_preascertainment.py`.
No job was submitted by this agent.

**Generator verification subsequently executed by root:** job 21997126
(root label 1568), one CPU, 2 GB, one hour requested, completed in 19 seconds
with exit 0 and empty error log according to root's scheduler/log check.
Reported RSS was zero and is unusable as a memory measurement. I inspected
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/synthetic_p1_generator_21997126/summary.json`:
all stated generator invariants passed. Both regimes retained all 180 people;
nuisance-off had 499 selected rows and nuisance-on 498. Each had nine
zero-selected people, including six zero-candidate people. Shared stage counts
were C1 [19,21,20], C2 [28,21,11], C3 [23,16,21]. The executed generator hash
is `afba32e07526e6785878bef2021d055efffb871ab5a6636be6844151312760a8`;
the executed likelihood hash is
`41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6`,
including root's checked memory-path change. This verifies synthetic
generation, selection and paired-design invariants. The summary records no
fit and no P1 size or power evaluation. No additional test was run by this agent.

### Feasibility of the prescribed null calibration

Read-only runtime/source assessment, 2026-09-30. Inspected only synthetic
runtime/source-hash metadata and current numerical reference code; no real
outcome, fit, job submission or new invariant test was performed.

Measured saved timings:

| Synthetic execution | Work timed | Elapsed |
|---|---|---:|
| `refit_failed_p1_21995505/summary.json` | One draw's prescribed warm/cold fit and sandwich calculations | 984.129 s |
| `refit_failed_p1_21992790/summary.json` | Same type of recovery for another draw | 1,039.902 s |
| `bfix_fits_21984981/null.json` | Earlier observed section-3 fit | 704.861 s |
| `bfix_fits_21984981/marginal.json` | Earlier section-4 fit | 121.826 s |

All paths above are relative to
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/`.
The two recovery timings use source `87f07515...`, 173 selected rows,
eight genes and 60 people. They are difficult draws, not a representative
runtime median or a lower bound. Current source `41ebc253...` removes unused
dense predictor-Hessian allocation. Its reported 0.188–0.190 versus
0.147–0.149 s evaluation timings are about 1.27-fold faster per evaluation;
this is not a measured 1.27-fold fit speedup. The live 21985999 job keeps
its original loaded implementation.

[Estimate, scenario rather than forecast] If those measured difficult-draw
times were typical, one P1 with 1,000 refits would cost 273–289 CPU-hours.
The prescribed 200 complete CRT null datasets require 200,000 such refits,
about 54,674–57,772 CPU-hours, or 6.2–6.6 CPU-years serially. At 32 fully
utilized independent CPUs that is approximately 71–75 days, before observed
fits, stage-law fits, failed refits, other null machinery, scheduling or
larger designs. This arithmetic does not authorize an array or establish
runtime at the new 180-person/499-row design. The earlier spec's 5–60 s
section-3 estimate was explicitly unmeasured and does not describe this
Python reference's observed difficult-draw performance.

The smallest justified next execution is an **initial implementation runtime
check**, not a reduced scientific calibration: fit one fixed nuisance-null
generator dataset under the current implementation, then, if that bounded
allocation permits, one fixed generated-stage refit. Record whole-fit time,
likelihood/gradient/Hessian evaluation counts and times, cold-intercept time,
optimizer iterations per restart, Newton/fallback time and sandwich time.
The existing generator check establishes data construction, not fitting cost.
One execution would locate the computational bottleneck but cannot estimate
size or justify weakening the fixed replicate counts/criteria.

A useful source-level hypothesis is redundant **cold intercept fitting within
P1**. `likelihood_v1.py:cold_start` runs three intercept optimizations per
gene, each using `conditional_terms(row, intercept_alpha, intercept_r, False)`.
That objective uses allele counts, depths, genotype priors, error/dispersion,
fixed mapping offset and the het predicate. It does not use row stage or z.
P1 changes stage while retaining every one of those inputs, so these cold
intercept solutions are invariant across its stage draws. Profiling should
measure whether this repeated work materially contributes to cost.

If worthwhile, precomputing the cold intercept vector once for the fixed
P1 count design and cloning it for each refit could preserve both prescribed
warm/cold starts, all three optimizer restarts, bounds, projected-gradient
tolerance and the higher-likelihood tie rule. This is a proposed numerical
implementation change, not done here. It must be restricted to exactly the
same count/tag/prior/offset/predicate design and independently check cold
vector and final-fit equality on a changed-stage case. P2 resampling and
P3 regenerated counts change the intercept objective and cannot reuse that
cache. Stage-law fitting/refitting rules remain untouched.

After measured fitting cost and any justified arithmetic optimization,
root can allocate the existing fixed null-replicate scheme with a concrete
CPU budget and independent replicate seeds. Do not substitute fewer CRT
draws, omit cold starts, loosen convergence, drop failed fits, tune a
criterion on interim size, or adaptively add replicates. If full execution
is infeasible, report calibration pending; arithmetic agreement and the exact
weak-null distribution diagnostic cannot replace it. No new test or job
was run for this assessment.

**Bounded runtime driver prepared at root's request:**
`scripts/analysis/in_biopsy_allelic/rec_codex/profile_synthetic_p1_fit.py`
and `run_profile_synthetic_p1_fit.sbatch`. The driver requires the reviewed
likelihood/generator hashes and an allocation. It regenerates the verified
nuisance-on 180-person/498-row design; compares cold-start vectors exactly
before and after one deterministic within-cohort stage rotation; then runs
one observed warm/cold fit initialized with `model.initial()`, never the
generating parameters. It records optimizer calls/iterations, cold-start,
score/full-Hessian evaluation, fit and sandwich times, versions and hashes.
Timings are nested/inclusive and must not be summed as independent costs.
The changed-stage case is only a cold-vector arithmetic check; no stage
refit or P1 resampling test is performed. No cache is implemented.

The launcher requests one CPU, 2 GB, 90 hours. `python3 -m py_compile` and
`bash -n` passed. Execution is **not run by this agent**; root owns submission
and its job index. This prepares the prescribed bounded fitting-cost check,
not null-size calibration.

### Runtime check completed; current full Gate B contract

Root's profiling job 21997263 completed in 8:08 with peak RSS 55,168 KiB
(root's scheduler report). I inspected
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/synthetic_p1_fit_profile_21997263/summary.json`.
The one observed nuisance-on fit converged from both starts under source
`41ebc253...`, retaining warm: log-likelihood −1471.33666488556,
kappa 0.0021651834988290573 and sandwich SE 0.04755455340975571.
Those are one synthetic realization, not calibration or an empirical result.
Fit time was 425.1579 s. The cold-start vectors were array-equal after
the fixed stage rotation; those computations took 25.20865 and 25.03017 s.
Cold-start construction inside the observed fit took 25.03971 s, only 5.9%
of its elapsed fit time. Reusing that start could therefore save a modest
fraction in this measured fit; it does not solve the dominant cost.
There were 911 score-only whole-model evaluations taking 395.9283 s,
six full-Hessian evaluations taking 3.87275 s, and a 1.16310 s sandwich.
The optimizer aggregate includes nested evaluations/cold starts and cannot
be added to those times. This is one observed fit, not a stage-draw runtime
distribution, full-scale benchmark, or scientific null calibration.
No cache or further numerical change is recommended without a concrete need.

**The binding Gate B decision is more specific than the old plan shorthand.**
`2026-09-29_model_A_likelihood_spec_v1.md:482` records the user accepting
B3 on September 29. Appendix B3, lines 489–494, requires:

- size/coverage bands assessed by warp-speed Monte Carlo at ≥4,000 replicates;
- 200 full-CRT replicates assessed one-sided, failing if the binomial lower
  95% bound on size exceeds 0.05;
- P2 boundary plants `kappa=±Delta*lambda` and `tau=lambda*tau_max`;
- technical plants for reference-bias slope, doubled false-het rate at S=2,
  non-beta-binomial tails with depth–stage coupling, tag switching, LD mixture,
  and omitted depth-dependent bias, plus a simulated label layer.

The earlier plan also supplies the kappa/tau plant grid, per-gene LRT QQ
checks, independent arithmetic agreement and the one-amendment stopping
rule. Later v1.5 requires tau bias/coverage before and after the n* restriction;
v1.6 supplies calibrated error/rules; v1.7 binds starts/selection of maxima.
Our eight-gene kappa=tau=0 nuisance-on/off generator and one fit meet neither
the full set of plants nor these calibration criteria.

**Construction feasibility:** generic synthetic RNA-mixture plants can be
created without any measured allele ratios, true participant genotypes or
sealed outcomes. Public GTEx effects, public tag frequencies/LD and approved
development total-depth/stage/QC records support an outcome-free design
approximation. But the inspected documents do not yet supply a complete
immutable simulation parameter file sufficient for the full prescribed test.
Do not choose these remaining details silently from the toy generator:

1. **Exact design and estimand roster.** Freeze the eligible gene list/theta,
   mapped development participants, candidate-tag universe and depth/selection
   rule, and whether calibration holds the observed metadata fixed or generates
   it. The toy roster/alpha/depth distribution is not the production design.
   Public GTEx effects are planning effects (with the declared cap), rather
   than empirically fitted biopsy effects. Do not use empirical allelic
   outcomes to fill this contract before the permitted stage.
2. **Joint covariate/prior design.** The approved export omits per-person
   ancestry axes and supplies unweighted cohort cross-products. Those moments
   cannot reconstruct the joint nonlinear selected RNA distribution of axes,
   stage, tag, depth and individual genotype priors. A generic synthetic
   ancestry distribution can test a declared approximation, but cannot be
   represented as the matched full design. The owning builder can run that
   full simulation privately and export aggregate calibration records; no
   restricted per-person genotypes need leave its environment. A weighted
   Fisher matrix is useful for planning, not a generator of the full joint
   distribution.
3. **Plant parameters.** Specify alpha/r_g/delta, fixed and stage-dependent
   mapping bias, error and homozygote dispersion, all nuisance stage effects,
   LD-mixture weights and depth–stage law, non-BB tail distribution, omitted
   depth-bias function and false-het perturbation. The acceptance triggers
   are named, but these distributional details and magnitudes were not fully
   defined in the inspected operative text. Generate u_g before reads where
   required, and state its centering/conditioning for the fixed-gene kappa
   estimand. The eight-gene nuisance magnitudes are diagnostics, not chosen
   substitutes for every named technical plant.
4. **Label/boundary contract.** Bind the two separate label-source error
   distributions and their generation, target true/observed stage law,
   information weights and adjusted R², Delta, planning tau_max and source-
   specific lambda definitions. B1 expressly uses public GTEx tau_max only
   for Gate A/B planning; the production biopsy-derived tau_max is computed
   later by its hashed rule. Keep boundary plants and before/after depth-guard
   tau checks on precisely that declared planning scale.
5. **Monte Carlo protocol.** Specify how warp-speed simulation constructs
   each size/coverage quantity, its resampling statistic and allocation across
   plants; a method name alone does not define a reproducible test. State
   which cells receive the ≥4,000 and 200 replicate counts, whether bands apply
   per cell, the binomial-bound calculation, seed mapping, uncertainty and
   failure handling. Do not pool successful/failing null regimes to mask a
   cell's size, or infer that one bootstrap draw supplies every required
   interval/coverage endpoint. Preserve the already accepted criteria and
   one-amendment rule; operational detail must instantiate those decisions.

The useful next scientific action is for the owning lead to bind this
parameter/protocol file from the approved source records, with source hashes
and the public-versus-builder-private distinction explicit. It may reveal
that a generic synthetic calibration and a matched-design calibration are
separate outputs. Until that contract exists, broad compute would test a
researcher's new choices rather than the complete prescribed Gate B. No
new model, test law, amendment, simulation or job was introduced by this review.

### Permitted target-readout effect export: source-code inspection only

Bounded filename and producer-code search, no effect, p-value, performance,
fold-0 row value or response count read. No fit/job was run. The canonical
`alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz`
and `alphagenome_campaign/week1-20260915/model/native/currin_variants.tsv`
mix all five chromosome folds by their producer contracts. Existing
`model_native_residual.read_labels` skips unwanted rows before numeric
conversion, but still scans a mixed CSV. That is not the preferred route
under the present requirement to avoid protected outcome strings entirely.

**An already-produced fold-0-free label artifact exists:**
`GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/allele_stack_inputs_21849187/development_labels_folds1to4.tsv.gz`.
Its producer, `scripts/analysis/masld_two_models/allele_nested_inputs.py:30`,
requires the canonical source to contain exactly folds 0–4, filters fold 0,
checks its removal, then writes this output at lines 38–39. This establishes
the intended development-only output contract; this review inspected code
and filename existence, not its effect values or row contents. It copies
the canonical table's additional fields, including significance fields,
so the owning builder should emit a slim requested-column export from
**this existing permitted artifact**, rather than expose its complete
columns or rerun a mixed-source reader in Codex.

Other existing restricted-population routes were located by filename only:

- `.../week1-20260915/model/native_lengths_21775189/calibration_training.tsv.gz`:
  producer `model_native_lengths.py:164` writes only folds 2–4 with
  `key,peak_id,heldout_fold,beta_alt` and raw native-score columns. Its rows
  are the old common-length coverage, not necessarily the entire new manifest.
- `.../week1-20260915/model/native_residual_21775216/development_predictions.tsv.gz`:
  `model_native_residual.py:608–615` writes fold 1 only after training choices,
  but includes fitted predictions and old feature-common coverage. Do not
  use it to select the new receiving population or carry old accuracy.
- `.../week1-20260915/model/native_lengths_21775189/permitted_variants.tsv`
  is also filtered to 1–4 before writing by `model_native_lengths.py:105–111`,
  but retains unnecessary source columns. No content was read here.

**Minimal builder export contract:** exact signed association beta plus
identity metadata only. Requested columns are
`key,lead_variant_id,chr,pos_hg38,ref,alt,peak_id,peak_start0,peak_end0,
heldout_fold,block_1mb,beta_alt`.
Declare `peak_start0/end0` as exact renames of source
`peak_start_hg38/peak_stop_hg38` (GRCh38 BED0 half-open), position as one-based,
and key normalization explicitly. `beta_alt` is the unchanged source FastQTL
ALT-dosage slope; EA must not reorient it. No p/q values, allele counts,
genotypes, source accuracy or old predictions are needed. If an association
SE is later needed, request a separately identity-qualified source nominal
join; do not fabricate it from the lead table.

The builder's receipt should name and hash that fold-0-free source and export,
assert only folds 1–4, uniqueness of exact variant–peak identities, complete
finite beta coverage, and exact metadata equality to the frozen census.
Root's metadata contract has **25,748 development keys before geometry** and
**25,477 keys in the fixed full-readout manifest**. These are census/manifest
requirements, not newly recounted label rows in this review. Export the
25,748 permitted source identities or an explicitly joined 25,477-row
geometry subset, and verify the mapping to the fixed manifest. Sequence
REF/alphabet failures remain further documented exclusions independent of
beta. Refuse missing/duplicate outcomes; do not repair coverage with an
effect-dependent inner join. A metadata receipt can certify the export before
Codex reads any allowed beta values.

**Fair comparison and obligations:** the primary geometry contrast is
full-peak native ATAC versus the same-process variant-center 501-bp ATAC,
using the exact same fresh REF/ALT tensors, tracks, input context and common
manifest. Each arm gets its own zero-intercept slope fitted in folds 2–4;
fold 1 remains already-inspected development evaluation. Root should fix
coefficients, output predictions, primary contrast/family and the useful-
effect rule before evaluation beta access. Include zero prediction and
simple training-selected center/peak combinations for any claim beyond the
equal-capacity geometry contrast. An eventual learned region adapter must
beat native full-peak scoring and the strongest matched simple calibration;
native ATAC+DNase is an additional existing strong comparator only when its
fresh same-input readouts are available and coverage is matched.

`model_native_residual_design.json` uses a **5% selected-native inner-OOF MSE**
operational margin for its prior sequence/residual/joint hypothesis, and
`allele_combination.py:232` uses a separate 5% relative-gain decision. Neither
automatically defines the new mask hypothesis's margin. **No target-mask
margin is selected at present.** The new comparison must declare its own
scientific/operational usefulness criterion without choosing it from fold-1
error. Preserve its development status and unresolved source/pretraining
independence regardless of that result.

Report paired MSE in beta-squared units, RMSE/MAE in source beta units,
signed association and exact exclusions/coverage. Existing residual design
uses a chromosome bootstrap as primary and historical 1-Mb-bin sensitivity,
conditional on fitted models and deposited beta; those bins are not verified
independent LD units and neither unit is a biological donor. Its complete
18-contrast BH family belongs to the earlier comparison and must not be
silently repurposed. Prespecify the mask comparison's primary single contrast
and any additional exploratory family. A geometry improvement does not prove
causal editing effects, new-variant generalization or independently replicated
accuracy. No protected fold is needed for this development comparison.

**Development-only exporter prepared:** at root's request added
`scripts/analysis/in_biopsy_allelic/rec_codex/export_target_readout_effects.py`
and `run_export_target_readout_effects.sbatch`. Header/hash inspection only
established the existing fold-0-free source SHA256
`72bfc08cb05dc0ff8d4e69580b1d20484ea6f9f2a3bb45336c8336edea1e9916`.
The exporter pins that source, the 25,477-row manifest and its 271 metadata
geometry exclusions. It first imports identity columns only and requires
their 25,748-key union, fold/coordinate/allele/block consistency and exact
manifest join; a second thin read imports allowed development beta strings
only after all metadata checks. Numeric parsing is limited to admitted
manifest rows and rejects nonfinite beta. No p/q fields, predictions or
metrics are imported. Output is the 12-column table above plus metadata
exclusions and fold-count/hash summaries in
`target_readout_effect_export_JOBID`.

The launcher requests account/QOS nslab, CPU, one CPU, 2 GB and ten minutes.
Static Python compilation and `bash -n` passed; export execution is not run
by this agent and root owns submission/indexing. No beta value was inspected
while preparing the exporter. This creates permitted input preparation only,
without selecting a margin or evaluating mask accuracy.

### Executed development export and future paired descriptive comparison

Root executed the development-only exporter as job **21997874**. Its completed
metadata summary reports 25,748 source keys, 25,477 exported geometry-admitted
keys and 271 geometry exclusions. The source fold counts for folds 1–4 are
6,834/6,693/6,153/6,068; exported counts are 6,762/6,631/6,095/5,989;
excluded counts are 72/62/58/79. The allowed output is
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/target_readout_effect_export_21997874/development_target_effects.tsv.gz`,
SHA256 `dc05d98a67b980ba7c34a64799d1b76d39d8dc58ba64d93c7b485908d731a7ad`.
The preceding "not run by this agent" statement describes source preparation,
not the current execution state. No beta value was inspected in this report.

At root's explicit request prepared these three new files:

- `scripts/analysis/in_biopsy_allelic/rec_codex/describe_target_readout_development.py`
- `scripts/analysis/in_biopsy_allelic/rec_codex/target_readout_development_protocol.json`
- `scripts/analysis/in_biopsy_allelic/rec_codex/run_describe_target_readout_development.sbatch`

The fixed comparison fits intercept-and-slope OLS separately to the fresh
three-track center and full-peak means using folds 2–4, and describes MSE,
MAE and Pearson correlation in previously inspected development fold 1.
Zero, training-mean and intercept-plus-center-plus-peak OLS are included;
there is no tuning or extra genomic covariate. All five methods use the same
**24,430** sequence-admitted keys; the expected **1,047** sequence exclusions
must exactly partition the 25,477-row fixed geometry manifest, with no
replacement or effect-based filtering. This is a native-readout comparison,
not an independently validated new learned model.

The driver refuses a missing or unfinished `target_readout_development_21997286`
summary. It checks the recorded producer/helper hashes, output types,
three-track order and native aggregation; matches effect/manifest keys,
folds, unchanged REF/ALT, GRCh38 peak coordinates, source variant identity
and sequence exclusion metadata; and verifies the common score identities
and full-peak widths. Fixed source hashes and completed-input hashes are saved
with the methods contract **before numeric beta is read**. No mixed-fold
source, p/q field, prior scalar prediction or other label file is opened.

The three fixed descriptive contrasts are peak OLS versus center OLS, joint
OLS versus center OLS and joint OLS versus peak OLS. Outputs include native
held-fold row predictions, fitted coefficients, per-chromosome counts/SSE
and paired SSE differences, and fold-by-geometry coverage. For 2,000
seed-20260930 draws, whole chromosomes are sampled with replacement and
their SSE/counts are summed before computing MSE differences and
candidate/reference MSE ratios. Percentile 95% intervals are descriptive;
few chromosome clusters do not establish donor independence or independent
LD blocks. A zero reference SSE makes the corresponding ratio undefined.
No p-values, scientific margin, pass/adoption criterion or biological
target/gene claim is introduced.

The future launcher requests **CPU / account nslab / QOS nslab / 1 CPU /
2 GB / 1 hour**. Static Python compilation, protocol JSON parsing and
`bash -n` passed. Analytical execution, numeric effect reads and submission
are **not run**. Root owns independent review and any eventual submission
after all paired scoring finishes.

## GSE128072 prior target exposure: bounded source and receipt audit

**Assessment:** the inspected candidate workflow supports **different-source,
cross-assay receiving transport with no recorded performance exposure**.
It does not establish universally untouched, person-independent confirmation.
No GSE128072 target-dependent fit, predictive metric or model-selection branch
was found in the bounded code/configuration search. Source metadata, donor
intersection, input-unit incompatibility and interval incompatibility were
already inspected and used to decide whether this source could be admitted.
That prior feasibility selection is distinct from choosing a model using its
observed prediction errors.

The strongest existing processed-source evidence is
`GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/source_metadata/admission.json`
and its current producer
`scripts/analysis/masld_two_models/source_admission.py:55`–76. The report
records 17 paired donor IDs and explicitly says "no performance measured";
it records the SOFT/RNA/H3 file hashes and the TPM/hg19 admission failures.
The producer reads the RNA header and H3 header, **then reads and splits one
complete first H3 feature row into strings**, using only `first[1]` to check
its chromosome prefix. Thus an assertion that no target row or target bytes
were ever opened would be false for this producer. It performs no numeric
conversion of target count cells, array construction, target distribution
calculation, prediction, error evaluation or count-based selection. Its
checksum helper reads complete compressed bytes for hashes; that is not a
measurement analysis. No target values were reopened in this audit.

The saved admission report does **not** include a producer-code hash or an
executed-source snapshot. [Unverified] Current source plus matching metadata
output cannot prove the exact historical code bytes or rule out every
unrecorded interactive inspection. The earlier report statement that this
agent inspected only headers remains correct for that follow-up, but does
not describe all historical producer I/O.

The later authoritative raw-roster execution is
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/paired_h3k27ac_raw_roster_21997742/summary.json`.
It identifies producer SHA256
`dc9da4b0993e31f1a3303229c387b3c45b118fa27e922697d9b6620dfdcbdbcd`
and SOFT SHA256
`92d47ccb362fec5913e899d765982f413c1b39327ab92b80a40af77536030def`.
`paired_h3k27ac_raw_roster.py` parses donor/cohort/title/assay/SRX/BioSample
fields from SOFT and fetches ENA run metadata only. Its receipt reports
17 donors, 468 runs, 604 FASTQ files and 187,002,556,838 listed FASTQ bytes,
with `raw_sequencing_acquired=false`. ENA read/base totals and file sizes
are acquisition metadata, not measured RNA expression or H3 region targets.
The roster contains no target or performance analysis.

The existing fixed `scripts/analysis/masld_two_models/design.json:4`–18 names
GSE267145 as the development source and GSE128072 as a pending paired
candidate. Current `chromatin_development.py`, `chromatin_acquisition.py`
and `chromatin_secondary_review.py` explicitly load the GSE267145 fixture
and/or its saved development outputs. They do not load the GSE128072
processed matrices. An unignored, accession/filename search across
`scripts`, `Analysis/MASLD_Model_Benchmark/config`, its `scripts` and `src`
(Python/JSON/TOML/shell/Slurm sources only) found GSE128072 references only
in the admission producer, two-model design and raw-roster producer. No
GSE128072-specific TOML/evaluator was found in those inspected branches.
This is scoped source evidence, not a repository-wide proof of absent aliases
or manual historical use; no mixed-effect or sealed-outcome file was searched.

**Scientific consequence before acquisition:** retain the already frozen
unweighted raw `rrr_offset_cis` primary and same-input raw global `rrr`
comparator, and the fixed B99 residual-target definition described in
`2026-09-30_codex_program_model_options.md`, rather than switch models or
assay normalization using the receiving errors. Recounting common regions
can resolve the processed-coordinate mismatch; it does not remove
ChIP-versus-CUT&RUN shift or prove the exact training RNA lineage.
[Unverified] Universal person disjointness and foundation/checkpoint exposure
remain unresolved from the available metadata. The collection-history
limitations are explicitly documented in that report's GSE128072 versus
training-people section. Calling the 17-person experiment an already-used
development evaluation is unsupported by the inspected producers; calling
it unrestricted independent confirmation would overstate the evidence.
This audit acquired no data, read no expression/target/effect values,
computed no performance metric and submitted no job.

### Independent Ensembl98 acquisition/census source review

Reviewed `scripts/analysis/in_biopsy_allelic/rec_codex/acquire_ensembl98_reference.py`
and `run_acquire_ensembl98_reference.sbatch`, source metadata receipts and
the modeled RNA-axis header/first identity records only. No reference
sequence/annotation contents, RNA/H3 measurements or acquisition execution
were opened by this reviewer.

**Coverage is an identifier census, not RNA transport admission.** At lines
225–231 the driver normalizes both `stable_gene_id` and `matrix_gene_id`
and separately records their exact membership in the release98 GTF,
cDNA and cDNA-plus-ncRNA union. `modeled_in_union=false` with
`source_in_union=true` can indicate a training-axis mapping mismatch rather
than absent reference sequence. Neither flag searches aliases, stable-ID
history, gene splits/merges or independently established sequence-equivalent
annotations. [Unverified] No such actual mismatch is established before
the census executes. If both flags are false, this still establishes only
absence of those normalized IDs from this reference, not biological absence
of the corresponding RNA. Retain the complete 42,163-model-gene requirement;
resolve any missing IDs through an independently supported, outcome-blind
mapping contract rather than manufacture zero rows or invoke the scorer's
lower coverage floor.

**Version and locus checks have explicit limits.** Lines 148/159 remove
GTF ID versions while retaining native PAR/locus suffixes. Lines 195–202
remove the FASTA transcript version for its GTF join, then require all GTF
gene assignments to agree with the normalized FASTA gene assignment. They
do not compare the FASTA transcript/gene versions against GTF version
attributes. The pinned metadata do select release98 human GRCh38 cDNA,
ncRNA and GTF together, so this is a qualification of the check, not an
observed conflicting-version failure. Exact transcript-to-gene membership
does not prove transcript lengths/quantifier settings reproduce the private
training generator. The same-release union is an approximate source contract.

**Normalization and duplicates are appropriately separated for this census.**
Line 53 strips terminal `_PAR_Y` then a terminal dot-version, matching
released `score.py:47`–51. Distinct raw transcript IDs remain distinct
(206–215), and repeated exact IDs must have equal gene assignment, length
and sequence hash. These rules protect the reference roster; they do not
authorize adding several remapped genes to one modeled gene or counting
two equivalent transcript records as independent RNA measurements. A
future quantifier-to-gene adapter must preserve the released duplicate-sum
semantics and the full modeled-gene denominator without target-based choices.

**A successful job can report incomplete coverage.** Lines 240–248 report
coverage totals/missing IDs, and lines 267–273 write a metadata-admission
summary without requiring `modeled_union_covered == 42163`. That is correct
for a census, but completion must not be read as permission to quantify or
score an incompletely covered primary axis. The summary explicitly leaves
index building/quantifier flags unselected and exact private PISCES
reconstruction unasserted. The CPU/account/QOS nslab launcher requests
one CPU, 2 GB and one hour. This source review found no reason to halt
reference acquisition; RNA transport remains conditional on the complete
axis and independently justified mapping/quantification contract.

### Frozen raw RNA-profile dependency and error-decomposition review

**For the specified `rrr_offset_cis/raw` versus `rrr/raw` comparison,
foundation-model pretraining exposure is not a dependency.** Source tracing
finds no pretrained expression embedding, sequence model, foundation
checkpoint or externally learned expression representation in these two
prediction routes. This narrows the earlier general checkpoint caveat:
it remains relevant to an actual foundation-model comparator or sequence
track, not to these two frozen RNA-only arms.

`release/masld-liver-chromatin-state-v1.2/score.py:190`–234 applies packaged
RNA transforms and RRR factors and, for the primary, adds the packaged local
ridge head. In
`executions/chromatin-stable-rrr-forms-20260908T192024Z/run_stable.py:90`–96,
both the global RRR and residual local head are fitted and deposited on
all 99 development participants. `fit_heads` constructs `res_cis` with
zero global/state columns and the nearest at-most-ten cis genes.
The scorer computes state/PCA intermediates, but the source-constructed
Ocis head does not use those state columns. Its RNA predictors and learned
parameters are those of the development cohort.

That execution's `transfer_functions.py:113`–148 obtains the GSE267145
substrate and selects gene TSS neighbors using frozen GENCODE v49 annotation
within 100 kb of region midpoints. `common.py:54`–73 in
`executions/chromatin-beyond-histology-20260907T161521Z` guards the 99-person
fixture and loads its RNA/H3 matrices. Source was read; those arrays were
not opened by this review. `transfer_functions.py:151`–160 fits its PCA/PLS
and ridge transforms on the supplied training indices; the referenced v1.1
PLS orientation also originates from the all-99 source fit in
`executions/paired-teacher-student-20260907T183136Z/build_release.py`.
`executions/chromatin-state-v12-build-20260908T172230Z/build_v12.py:55`–61
packages the stable RRR/Ocis factors and recomputed training transforms.

**Training-person exposure and historical model choice are separate.**
The build's line 62 takes the primary form from a decision that considered
development and previously inspected GSE296875 transfer. Thus "no
externally pretrained coefficients" does not mean model-form choice was
unexposed to every external source. GSE128072 still has no recorded
performance exposure in the bounded audit above, but universal person
disjointness from the 99 training people remains unverified. Those are the
relevant receiving-source qualifications; a Corgi/AlphaGenome training
corpus is not an input to this fixed pair of models.

**The new error split is algebraically correct and descriptive.** For error
`e_ir = Y_ir - P_ir`, and regional donor mean `m_r = mean_i(e_ir)`,

`mean_ir(e_ir^2) = mean_r(m_r^2) + mean_ir((e_ir-m_r)^2)`.

The cross term vanishes only after averaging across donors because
`sum_i(e_ir-m_r)=0`. `profile_error_decomposition.py:23`–37 implements
that identity with divisor n, not n−1. Its second term is **within-region,
across-donor error**, equivalently error between cohort-region-centered
target and prediction. It is not within-person profile-shape accuracy;
donor-wide prediction errors may remain in that term. The two nonnegative
components do not separately add to each individual donor's MSE without
a cross term. Biological replication remains the donor, with no new
independence implied by the many regions.

Receiving regional means are used only to describe the error after frozen
prediction, not to calibrate a receiving prediction route. The unchanged
native residual-profile donor MSE remains primary. The centered diagnostic
must not replace it or be advertised as deployable calibration, and a
receiving diagnostic improvement cannot select another model/transform.
Root's completed synthetic invariant receipt
`profile_error_decomposition_21998406/summary.json` reports the three-person,
two-region identity/offset/permutation checks passed with no receiving data
read (producer SHA256
`f75a8ebacd9383962d11d16849efb153b621e330f6947b9abc22876fe2601d29`).
This review reran no tests, read no model arrays or receiving measurements,
and changed no primary hypothesis, prediction or loss.

### Salmon v1.10.3 short transcripts: output presence is not observability

Root reports metadata-only job 21998488 found 46 reference transcripts of
length ≤31, including 40 modeled genes with no longer transcript, minimum
length 8, and 408 cross-gene identical-sequence groups affecting 602 modeled
genes. This review traced official source only; it read no reference
sequences, quantifications or gene-count values and ran no index/quantifier.

**The relevant Pufferfish version is the pinned tag, not master.** Salmon
v1.10.3's [fetchPufferfish.sh](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/scripts/fetchPufferfish.sh#L25)
sets `SVER=salmon-v1.10.3` and supplies the expected source-archive SHA256
`52b6699de0d33814b73edb3455175568c2330d8014be017dce7b564e54134860`.
The following line references are one-based lines from official raw source
retrieved into memory, not from a installed index or a different branch.

**Retained short names are restored, not their graph-mapping capability.**
Pinned [FixFasta.cpp](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/src/FixFasta.cpp#L212)
lines 212–216 flags post-clipping length ≤k; lines 386–394 omits those
sequences from cleaned FASTA while recording names/lengths in `shortRefs`
and accumulating `refIdExtensions`. Thus the raw length census is not
necessarily the complete final short set: poly-A clipping can shorten
additional references. Entirely clipped sequences and sequence duplicates
removed without `keepDuplicates` are separate cases, not guaranteed short
name restoration. Pinned
[PufferfishIndexer.cpp](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/src/PufferfishIndexer.cpp#L543)
passes these records to serialization.

[PufferfishBinaryGFAReader.cpp](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/src/PufferfishBinaryGFAReader.cpp#L327)
lines 337–384 merges short lengths/names back into the full reference tables
in order and serializes the extension offsets. The short references still
have no cleaned-FASTA graph sequence and therefore cannot receive the
ordinary index-seeded fragment matches. Restored reference-table membership
does not establish a measured gene estimate.

**Normal quantification writes zero-valued rows for these retained but
unmapped short references.** Salmon v1.10.3
[ReadExperiment.hpp](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/include/ReadExperiment.hpp#L253)
lines 253–290 builds transcript objects from the *full* names/lengths,
including the `length <= k` records. In
[CollapsedEMOptimizer.cpp](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/src/CollapsedEMOptimizer.cpp#L582),
lines 582–614 marks active only transcript IDs appearing in fragment
equivalence classes and initializes inactive estimates to zero. Updates
accumulate only supported equivalence-class assignments; lines 1039–1040
store final assigned counts/masses. The standard
[GZipWriter.cpp](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/src/GZipWriter.cpp#L1432)
abundance writer at 1448–1483 iterates all transcript objects and writes
their estimated counts/TPM. [Inference from pinned source, not an executed
quantification] With successful ordinary quantification, retained short
transcripts absent from all mapping classes appear in `quant.sf` with zero
`NumReads` and zero TPM. This is a **structural zero under that mapping
recipe**, not evidence that the donor had zero biological RNA. Source
initialization priors do not supply fragment support for missing mapping
classes. A complete set of output names would therefore not satisfy the
complete *observed-input* primary contract by itself; no zeros should be
invented, relabeled as observed, or silently accepted to reach 42,163 rows.

**Alignment-based Salmon is an RNA-native alternative, with unresolved
observability rather than guaranteed complete measurement.** Official
v1.10.3 [documentation](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/doc/source/salmon.rst#L222)
describes transcript FASTA plus transcript-coordinate SAM/BAM inputs,
without Pufferfish indexing, and warns against passing unconverted genomic
alignments or target/position-sorted alignments. Its
[AlignmentLibrary.hpp](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/include/AlignmentLibrary.hpp#L153)
lines 153–166 constructs targets from the alignment reference header and
populates them from FASTA, without this Pufferfish short-reference filter.
A genome-first RNA alignment route must produce compatible transcript-
coordinate alignments and the matching annotation/FASTA, or reads must be
aligned directly to transcripts. This can preserve full target
representation, but the chosen RNA aligner, seeds, minimum alignment
length/score, pairing/clipping rules and library capture still determine
whether length-8-to-31 transcripts can receive admissible evidence.
[Unverified] No such complete observed-gene coverage has been demonstrated
here. Alignment-mode fractional estimates would remain a qualified
measurement approximation, not reconstruction of private PISCES inputs.

**Exact sequence ambiguity persists in either route.** Keeping identical
transcript records preserves names but cannot identify separate expression
from the same sequence evidence. Allocations between exact identical
transcripts are not individually identified by those fragments. The 602
affected modeled genes need not all lack other unique transcript evidence;
the group census alone does not prove every corresponding gene total is
wholly unidentifiable. Neither alignment mode nor EM allocation removes
that distinction. Complete modeled-axis mapping and the existing
observed-input contract remain in force; these source findings authorize
no index, alternative quantification job or downstream target evaluation.

### Short modeled genes: prediction activity and a fixed measurement decision

**Source establishes denominator participation, not each gene's direct
coefficient activity.** In frozen `score.py:179`–190, every matched modeled
gene contributes to the input sum `T = sum_g C_g` before log2 CPM. Thus all
40 short modeled genes are potential prediction inputs through T even if
their direct effective coefficient happens to be zero. Supplying an
unobservable zero or dropping a genuinely expressed gene changes the CPM
of all other genes. This is relevant to both fixed raw arms and does not
depend on local-head inclusion.

For complete-axis raw scoring let `x_g = log2(1 + a*C_g/T)`, a=10^6. Source
defines the direct global coefficient matrix

`B_R = diag(1/prof_rrr_xs) @ prof_rrr_VtDv @ prof_rrr_WVr @ prof_rrr_Vr.T`.

The local coefficient `B_L[g,r]` is the sum of
`prof_Ocis_coef_cis[r,k]/prof_x_sd[g]` over entries with
`prof_Ocis_cis_idx[r,k]=g`. Source construction gives the Ocis global/state
coefficient block zero columns; the package audit must verify that shape.
The global arm uses B_R, the local arm B_R+B_L, with frozen intercepts.
No source assertion prunes the 40 genes from these factors. [Unverified]
Whether any particular short gene has a nonzero direct effective row
requires reading the actual released arrays on compute.

Even if the entire direct row B[j,:] is zero, count sensitivity generally
persists through normalization. For a valid positive total T,

`dP_r/dC_j = (a/ln(2)) * [B[j,r]/(T+a*C_j) - sum_g B[g,r]*C_g/(T*(T+a*C_g))]`.

This is a mathematical dependency, not an observed size of the effect or
a claim that every gene is expressed. Complete input admission cannot be
relaxed from a direct-coefficient audit alone.

**Exact bounded compute audit for root, not run here:**

1. Pin the released
   `release/masld-liver-chromatin-state-v1.2/weights/chromatin_state_v1_2.npz`
   to its recorded SHA256
   `72f9eaba359f7c4c61cac2e00da2d5539b5d3ec868dc33b1bacc2ef0105694ec`
   and the modeled axis to `baf983ad18a893b7e5db04364818073facc68b47ae25f63070373481e3b217dd`.
   Read the 40-gene metadata set from
   `executions/codex-rec-20260929T142434Z/ensembl98_salmon_metadata_21998488/modeled_genes_no_unclipped_k31_transcript.tsv`,
   root-provided SHA256
   `4f348047d98367f9e041716dbe9a7a8952adcf54bba78c0dc97f55477151f7bf`,
   verifying the completed summary's output hash and exact normalized join
   to the full packaged `gene_axis`. No receiving data are needed.
2. Require finite positive scaling arrays and the expected factor/axis
   shapes. Form only the 40 selected RRR rows: first multiply selected
   VtDv rows by WVr, then by Vr.T. Accumulate corresponding local cis
   slopes and evaluate the combined rows. Report per-gene cis inclusion,
   actual nonzero-coefficient entries, maximum absolute slope and row norm
   separately for RRR/local/combined, plus any combined cancellation.
   A 40-by-96,460 float64 block is about 31 MB; never materialize the
   complete 42,163-by-96,460 coefficient matrix. This reads trained weights,
   computes no fit/prediction error, and uses no targets.
3. If already authorized root needs training count relevance, use only
   `executions/model-data-064-21079902/fixture/molecular/rna_values.npy`
   and its exact 99-person/42,163-gene axes. On compute, validate the
   fixed source/axis identity, obtain short-gene positive-estimate counts
   and each donor's total short-set estimate divided by the complete
   modeled-gene total. Retain fractional estimate semantics. This is an
   allowed-development input-mass audit, not an expression/H3 accuracy
   analysis or a license to transport structural zeros. It needs no H3
   counts, histology, genotype or receiving array. Report aggregate
   summaries without printing individual count cells. No audit code or
   job was added by this reviewer.

If short genes have nonzero *training* fractional estimates, that exposes
a mismatch with this particular reconstructed-reference/k31 eligibility
scenario, not a proved private-generator error. Different annotation or
version, private transcript filters, another mapping/quantification recipe,
and aggregation/alias/multimapping attribution remain distinguishable
possibilities. Input positivity alone cannot choose among them. Conversely,
all-zero training entries establish neither biological absence nor zero
receiving expression. Keep any ensuing audit interpretation conditional on
the admitted metadata and actual source-count semantics.

**Root's executed activity audit supersedes the preceding "unverified"
individual-activity status:** job **21998515/1585** completed in four seconds,
exit 0 with empty stderr. All 40 short genes have zero fractional RNA input
in every one of the 99 development people, zero direct global slopes and
total short-set RNA mass zero. Exactly one local direct coefficient is
nonzero in floating arithmetic: ENSG00000237235, one region, absolute slope
`1.4471263907714508e-21`. This is a numerically negligible reported coefficient,
not a newly selected exclusion threshold. The complete 42,163-gene training
log2-CPM means agree with the packaged means to maximum absolute difference
`3.55e-15`. Root reports no fitting, predictions or receiving data use.

These training zeros are consistent with the k31-ineligible reference
scenario; they do not prove the private quantification recipe or biological
absence in a receiving donor. They also show that demanding independent
assay observability for every receiving modeled row may differ from the
original training *estimate* semantics for these genes. Resolving that
measurement-definition distinction requires an explicit scientific
decision before receiving evaluation, not silently inserted zeros, a
changed denominator or a claim that the current complete-observed-input
condition has been satisfied. That condition remains in force.

**A frozen external adapter can be specified without claiming private
PISCES equivalence; a defensible complete observed-input run remains
unproved.** The coherent candidate route is RNA alignment plus Salmon
alignment-mode fractional estimation, using the single qualified complete
Ensembl98 primary-genome/GTF/reconstructed-transcript/tx2gene bundle for
both frozen raw predictors. Freeze that bundle, gene-version/PAR handling,
library mode, alignment restrictions and aggregation before receiving H3
values are accessed. Installed versions are documented as STAR 2.7.11b,
Salmon 1.10.3 and tximport 1.34.0 in the companion model-options report.
No receiving-specific gene filter or predictor refit is required.

STAR's official [2.7.11b parameter source](https://github.com/alexdobin/STAR/blob/2.7.11b/source/parametersDefault)
distinguishes `TranscriptomeSAM` from integer `GeneCounts`. Its
transcriptome output has explicit pairing/indel/softclip restrictions;
those settings must be frozen rather than assume all genome alignments
become Salmon evidence. Use compatible transcript-coordinate alignments
and matching transcript FASTA for Salmon, preserving read-pair ordering
without genomic coordinate sorting. Retain the complete transcript header
and explicit identical-sequence multiplicity; gene counts must come from
fractional Salmon estimates through a fixed, inspected tximport aggregation,
not replacement by STAR integer GeneCounts. Preserve sample-specific
effective lengths and distinguish source aggregation from a new adapter.
The public PISCES tag's transcript filtering/default annotation remains
different from the private generator; this route is expressly an RNA-native
measurement approximation.

**Scientific stopping condition:** neither removal of the Pufferfish cutoff
nor a complete BAM header proves that all 40 short genes can receive valid
RNA evidence from these polyA paired-read libraries. Seed/match limits,
fragment compatibility and sequence identifiability remain load-bearing.
Receiving output rows with structural zeros must be flagged as unobservable,
not counted as observed full-axis coverage. A predeclared aligner/quantifier
contract and outcome-blind evidence of input observability are still needed
before the original complete-observed-input primary can proceed. If that
condition fails, report the primary as inadmissible under this route;
do not switch denominators, impute zeros, drop genes, exploit their direct
  coefficient values, or substitute a partial-input accuracy claim. Identical
transcript evidence remains individually ambiguous in alignment mode too.
This source review does not choose model/normalization by expected receiving
accuracy or authorize acquisition, index building, fitting or jobs.

### Fixed whole-genome-decoy selective-alignment adapter: source review

This candidate follows the **corrected complete-native-estimate contract at
the top of this report**. It has not been indexed, quantified or evaluated.
Installed Salmon 1.10.3 `index --help` and `quant --help-reads` were inspected;
installed tximport 1.34.0 formals and its list `summarizeToGene` method were
read without loading counts. Help emitted a version-server "Not Found"
message; the installed command help itself was returned normally.

**Exact fixed flag contract:** use the qualified reconstructed Ensembl98
227,368-transcript FASTA first, append the entire verified 194-contig primary
genome FASTA, and provide every appended genome contig name exactly once in
the decoy list. Headers must be unique and transcript/decoy ID sets disjoint.
The index flags are `-t gentrome.fa -i INDEX -d decoys.txt -k 31
--keepDuplicates --no-clip -p THREADS`. `--keepFixedFasta` may retain the
processed reference for checking, without changing inference. Native
versioned ENST names do not require `--gencode`. The fixed paired-read
quantification flags are `-i INDEX -l A -1 MATES1 -2 MATES2 -o OUTPUT
-p THREADS`, with matched technical-run order for both mate lists.
Library type precedes read arguments; its inferred result and compatibility
summary are retained. Explicit `--validateMappings` is deprecated/no-effect
in installed 1.10.3 because selective alignment is already default.
No new `--seqBias`, `--gcBias` or `--posBias` choice is introduced.
The official [v1.10.3 documentation](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/doc/source/salmon.rst)
supports whole-genome decoys, default selective alignment and `-l A` inference;
receiving target accuracy is not consulted for these choices.

**`--no-clip` is scientifically interpretable, not default equivalence.**
Pinned [FixFasta.cpp](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/src/FixFasta.cpp#L184)
normally strips terminal A runs when the last ten bases are A, including
removing an entirely A reference; the flag disables that operation. It is
an outcome-blind choice to retain annotated reference sequence and lengths,
rather than change a transcript's annotation by that heuristic. It is
therefore sensible for this *declared* reference-preserving approximation,
but not proof of superior measurement, private PISCES settings or exact
source equivalence. `--keepDuplicates` prevents transcript-sequence aliases
being silently dropped; those aliases remain individually unidentifiable
from identical sequence evidence. Short transcripts still lack indexed
graph sequence and may have native structural-zero estimates, now retained
and flagged under the corrected contract. Non-ACGT cleanup remains the
quantifier's documented fixed pseudorandom replacement and must be recorded;
these two flags do not assert byte-identical index sequence for every input.

**Actual resource basis, from metadata only:** the complete reconstructed
transcript metadata sum is **357,475,422 nucleotides**, and the verified
genome summary gives **3,099,750,718 nucleotides**, so the combined sequence
has **3,457,226,140 nucleotides** before index preprocessing. The genome
FASTA occupies 3,151,425,857 bytes; its verified SHA256 is
`78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c`.
This basis is genome-sized, not a 357-Mb transcript-only index. Peak RAM,
temporary graph files, final index size and runtime are **not measured**.
The Salmon help uses automatic distinct-kmer/Bloom-filter sizing unless
overridden. FASTA size alone cannot certify the construction allocation;
no new sparse/filter-size optimization, allocation or submission is proposed
by this review.

**Native fractional-count aggregation:** import complete Salmon `quant.sf`
with tximport `type="salmon", txOut=TRUE, countsFromAbundance="no",
ignoreTxVersion=FALSE, ignoreAfterBar=FALSE`, then `summarizeToGene` with
the exact versioned transcript→raw native gene map and again
`countsFromAbundance="no"`. Retain all native transcript estimates and
effective lengths, including zeros; no expression/transcript filter is
added. Use the summed fractional `counts`, not TPM, length-scaled TPM
counts, integer rounding or an independent Salmon geneMap aggregation.
Apply the released gene-version/PAR normalization and duplicate-sum rule
only when joining the modeled axis; report all original gene identities.

The installed `summarizeToGene` method subsets to tx2gene matches with only
a message if quantified transcripts lack mappings. Consequently strict
upstream checks must refuse any missing quantifier/map ID rather than
allow that implicit omission. The method uses gene-wise `rowsum` of
transcript counts. Zero-abundance gene effective lengths use its native
length-metadata fallback; this does not insert or alter count rows. The
actual private transcript filters and quantifier/annotation recipe remain
unidentified, so these fractional estimates are a qualified native
Salmon/tximport approximation to the training measurement.

**Successful admission conditions before frozen receiving evaluation:**

- Freeze and verify reference/genome/annotation/map/index/producer hashes,
  exact command flags, installed versions, ordered paired-run roster and
  biological donor grouping. Decoys must follow all transcript records,
  match the full genome roster and remain excluded from transcript output.
- Verify all 227,368 expected transcript IDs appear once in every donor's
  native `quant.sf`, with no extra nondecoy ID or silently removed alias;
  transcript lengths and exact mapping must agree with the frozen metadata.
  Require finite nonnegative estimated counts/TPM and valid effective
  lengths. Retain zero rows and identify native short/ambiguous estimates
  without reinterpreting them as biological absence.
- Require the exact full transcript-to-native-gene join before aggregation
  and the complete 42,163 modeled estimate-ID axis afterward. Do not append
  missing genes, impute from development, switch gene versions or change the
  modeled-gene total. Require positive total native modeled abundance for
  every donor. Both frozen raw arms receive that identical matrix.
- Save inferred library compatibility and mapping/decoy QC without choosing
  a different library mode, transcript filter, bias flag, annotation or
  denominator by receiving H3 errors. The input adapter remains separate
  from the fixed B99 residual target and cross-assay target-count contract.

These conditions admit a complete **computational estimate axis**, not
every gene's biological detectability or independently identified true
expression. They do not authorize indexing, raw acquisition, quantification,
new package installation, prediction fitting or outcome inspection.

## Model A P3 numerical tails: exact isolated log-shape probe

Source-only review of the immutable
`REC/bfix_p3_21998519/executed_sources/likelihood_v1.py` and
`REC/1587_21998519.err` identifies numerical zero beta shapes at lines
32–40 and 99–114. Paired `expit(eta)` / `expit(-eta)` avoids ordinary
subtraction cancellation but eventually underflows; the resulting shape
reciprocals and digamma/trigamma terms produce nonfinite derivatives.
Trigamma overflow can occur even before a positive shape underflows.
Unbounded `delta` and adaptive random-effect mode `u` permit extreme finite
eta despite bounded gene alpha (snapshot lines 203–211, 405–419).

The captured stderr consists of warnings without a call stack.
`check_bfix_tests_v13.py` lines 91–103 assign `fit` only after the whole
`fit_marginal` returns, and generation also calls `bb_terms` at likelihood
line 531. Therefore `fit=None` alone cannot distinguish generation from
fitting. [Inference] Optimizer/adaptive-mode trial parameters are a plausible
source, but the failure phase remains unverified. An offending-call trace
would be necessary to establish it. Conservative failed-draw treatment is
unchanged and does not repair the failed numerical calculation.

For the exact integer beta-binomial model, write `T=log(t)`,
`A=T-logaddexp(0,-eta)`, `B=T-logaddexp(0,eta)`, and

```
R_s(l) = sum(j=0,...,s-1) log(exp(l)+j),  R_0=0
D_s(l) = 1[s>0] + sum(j=1,...,s-1) expit(l-log(j))
Q_s(l) = sum(j=1,...,s-1) expit(l-log(j))*expit(log(j)-l)
logmass = logchoose(n,k) + R_k(A) + R_(n-k)(B) - R_n(T).
```

The `j=0` term is exactly `l`, with first derivative one and second
derivative zero; every other log term is `logaddexp(l,log(j))`. No small
shape, reciprocal shape, or shape special function is constructed. This
is the integer rising-factorial identity, with full binomial constants,
including both count endpoints, without clipping or changing bounds.

Let `Da/Db/Dt` and `Qa/Qb/Qt` denote those derivatives at `(k,A)`,
`(n-k,B)`, `(n,T)`. With `m=expit(eta)`, `q=expit(-eta)`,
`u=t'/t` and `w=t''/t-u²`, exact chain-rule derivatives are:

```
score_eta = Da*q - Db*m
score_r   = (Da+Db-Dt)*u
H_etaeta  = Qa*q² + Qb*m² - (Da+Db)*m*q
H_etar    = u*(Qa*q-Qb*m)
H_rr      = u²*(Qa+Qb-Qt) + w*(Da+Db-Dt).
```

This retains the existing EPS-bounded rho link. Computing cumulative
rising-factorial tables once per shape uses O(n) temporary storage, rather
than forming one length-n sum for each allowed count. For a production
implementation, time/memory implications for actual read depths still
need measurement; this probe does not change that implementation.

Prepared new isolated files in `scripts/analysis/in_biopsy_allelic/rec_codex/`:
`codex_stable_bb_probe.py` and `run_codex_stable_bb_probe.sbatch` (reserved
job name 1599, CPU1, 2GB, 1h, nslab). The synthetic-only probe imports
the immutable executed-source arithmetic for moderate-predictor agreement,
and uses SciPy beta-binomial mass as an independent moderate reference.
It checks eta `-1000,-40,-6,0,6,40,1000`, linear-r
`-1000,-15,-3,5,1000`, n=20/200, endpoints/interior counts,
REF/ALT reflection, full normalization, actual integer truncated supports,
conditional derivative identities, finite-difference score/Hessian and
three-class inclusion-reweighted finite monotone CDFs. No fits, fixture
rows, actual genotypes, outcomes, model bounds or test law are changed.

Python static compilation and `bash -n` passed on login. Compute checks
**not run** by this agent; launcher not submitted. Arithmetic checks can
establish representation accuracy/stability on their fixed synthetic grid,
not optimizer convergence, P3 size/power or scientific calibration. After
an owner-reviewed replacement, a separately authorized fixed failed-draw
reproduction with a first-offending-call trace would distinguish the
failure phase and assess whether another optimizer issue remains.

## RNA-only B1 feasibility pilot source review

Read-only review of `codex_native_rna_pilot.py`, its job1592 launcher,
`codex_native_tximport_gene_counts.R` and frozen release `score.py` found
no decision-changing defect in ordered pooling or frozen prediction units.
The eight required SRR8702431–SRR8702438 pairs are sorted by accession and
mate number, so both Salmon mate lists follow the same order (pilot
170–191). The native adapter requires exactly 227,368 distinct versioned
transcript IDs, checks imported NumReads against native estimates and
checks gene sums against direct transcript `rowsum` (R 47–81). It supplies
no missing rows or structural-zero padding. Normalized gene extraction
requires all 42,163 modeled IDs before either raw head is scored, and the
released scorer computes the CPM denominator from all those input genes
(pilot 113–135; scorer 179–199). Frozen `rrr` and `rrr_offset_cis` receive
the same matrix without fitting or recalibrating predictions.

Two concrete improvements were sent to root before submission:

- Frozen transcript `length` exists in the mapping, but the adapter
  discards it at R 33–39 and does not compare quant.sf Length with native
  metadata at R 47–56. Exact native transcript-length agreement should be
  checked by transcript_id, consistent with the frozen no-clip reference.
- Pilot records the adapter SHA at line 204 but does not check it after
  subprocess execution. Capture that hash and require it unchanged after
  import so the recorded aggregation method is the executed method.

The saved-output cap is monitored during Salmon/R subprocesses; it is not
a byte-exact filesystem quota. Acquisition is independently bounded by
each frozen source file size, and all 16 source files total 6,261,964,285
bytes. No receiving H3 values, new fit or prediction accuracy are read or
computed. B1 is one biological donor with eight technical runs. The
declared Salmon/tximport assay approximation remains distinct from private
PISCES equivalence, biological detectability and independent validation.
This agent did not read count/package arrays, acquire raw RNA, execute the
pilot or submit its launcher.

## Executed stable BB probe and isolated first-draw replay preparation

The preceding probe's **not run** statement describes its source-only
preparation. Root subsequently executed job21998656: reported SLURM wall
16 seconds, exit0, empty stderr. Directly inspected
`REC/codex_stable_bb_probe_21998656/summary.json` reports all 70 fixed tail
grid cases passed, with 0.7400 seconds inside the Python probe. Maximum
absolute differences against the original moderate-predictor arithmetic
were logmass 2.6574e-12, gradient 4.8838e-10 and Hessian 2.6421e-10;
SciPy logmass difference was 2.7485e-9. Full log-normalization error was
2.5051e-12; finite-difference score/Hessian errors were at most 1.7305e-7.
The probe source hash is
`723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488`.
These are observed arithmetic results on a fixed synthetic grid, not P3
calibration or optimizer convergence evidence.

Prepared `codex_replay_first_p3_stable_bb.py` and
`run_codex_replay_first_p3_stable_bb.sbatch` in the Codex-only directory.
Reserved job name1593 requests CPU1/2GB/90h/nslab; **not submitted** by
this agent. Python `compile()` and `bash -n` passed. Source preparation
read immutable executed code, not synthetic fixture values or arrays.

On compute, the wrapper first verifies exact immutable P3 source hashes,
published toy fixture hashes against the original receipt, original null
warm-start hash and reviewed stable-probe hash. It copies modules into
its own output directory and replaces only that process's `bb_terms`
with the exact log-shape implementation. It compares observed-parameter
conditional LL/score/Hessian and marginal LL/score, reruns the unchanged
observed fitting prefix and refuses observed LL/kappa/v differences above
1e-6. The original parameters initialize the optimizer exactly as in the
provided runner; generator truth is never read or used as an optimizer
start.

It then uses only the supplied first draw's two uniforms and unchanged
conditional generation, `fit_marginal(..., section4_from_null=True)`,
two-start policy, shifted optimizer candidates, bounds, retained-start
rules, LR definition and conservative failed-draw exceedance treatment.
The historical `fit=None/LR=NaN/exceed=True` result is recorded separately
and remains unchanged. Phase, first warning/exception call stack and first
extreme-eta BB call are saved, enabling generation-versus-fitting diagnosis.
No full 50-draw replay, Monte Carlo p-value, size/power claim, live source
edit or running-job cancellation is authorized by this artifact.

## Receiving-versus-training participant identity: public metadata recheck

**Conclusion:** the public collection evidence supports different-source
transport, but does not establish that all 17 receiving people are disjoint
from the 99 training participants. This independently confirms the existing
collection-history section of `2026-09-30_codex_program_model_options.md`.

| Source | Supported procurement/participant facts | Date evidence and its limit |
|---|---|---|
| GSE128072 Penn Cohort 1 | Prospective liver-donor collection at Penn Transplant Institute; 50 donors, 49 deceased and one living | August2014–February2015, explicitly stated in the primary Study Subjects methods; this is cohort procurement time, not a sequencing-submission date |
| GSE267145 source cohort | Liver wedge biopsies obtained during elective bariatric surgery; VUMC referral population and Vanderbilt IRBs090657/171845 | Main Patients and sample collection/Study approval methods cite two registered studies but do not allocate the complete99 training subset to an individual procurement interval |

Sources: [Penn primary methods, Study Subjects/Penn Cohort1](https://pmc.ncbi.nlm.nih.gov/articles/PMC6612522/),
[VUMC primary methods, Patients and sample collection/Study approval](https://insight.jci.org/articles/view/182232).
These statements concern procurement and enrollment, not molecular findings.

The frozen GSE128072 SOFT was rechecked using only sample title/cohort/person
fields for GSM3663220: `B50_living tissue donor_RNA-Seq`, Penn Cohort1,
individualB50. The fixed17 paired roster excludes this RNA-only donor.
[Inference] Together with the primary study's exactly-one-living description,
the mapping places the retained17 among its deceased organ donors. This is
an author-metadata roster deduction, not per-person independently verified
vital-status documentation. Neither this status nor the center proves
histologically normal liver for every retained donor.
[GEO sample record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3663220).

Both official ClinicalTrials.gov API records were independently reread with
300,000-byte response caps, extracting only study start and location fields:

- [NCT00983463 official API](https://clinicaltrials.gov/api/v2/studies/NCT00983463):
  `startDateStruct.date=2009-10`, Vanderbilt University Medical Center,
  Nashville, Tennessee.
- [NCT03407833 official API](https://clinicaltrials.gov/api/v2/studies/NCT03407833):
  `startDateStruct.date=2018-02-13`, `type=ACTUAL`, same VUMC location.

The older registry begins before Penn collection. Thus the fact that one
VUMC registry starts in2018 cannot establish that every training participant
was alive after every receiving donor died. Earlier living bariatric biopsy
followed by later deceased organ donation is not excluded by the inspected
dates. Different institutions, accession labels and submission years are
contextual evidence, not a person-identity exclusion.

**Exact unresolved metadata:** a documented biopsy-date lower bound after
February2015 for all99 training participants, a source-owner allocation
establishing that all99 were enrolled/biopsied under the2018 study, or an
authoritative cross-source person-overlap exclusion. The permitted metadata
do not provide that link. Until it exists, state different-source transport
with unresolved universal participant disjointness; do not label this as
unrestricted independent confirmation or universally untouched data.

Only public procurement/enrollment metadata and existing source-admission
reports were inspected. No receiving RNA/H3 values, molecular array,
genotype, target loss, effect-result table or figure was read for this check;
no author contact, acquisition, fit, job or authority edit occurred.

## Successful B1 RNA pilot and fixed remaining-donor batch preparation

Root reports B1 pilot21998663 completed34:06, exit0, MaxRSS5,546,900KiB,
with all227,368 native transcripts, all42,163 modeled inputs and both
96,460-region frozen raw predictions. The native adapter returned60,623
genes, conserved transcript/gene estimated-count totals and retained46
structural-candidate zero transcripts. Original pilot source remains
`413b9617160e0f472682aecc21a629002fae0f54556e1244ffc3aac044e232c5`.
These are RNA measurement/admission and prediction-generation results for
one biological donor, not receiving H3 accuracy or independent validation.
The all17 cohort remains incomplete.

Prepared new `codex_native_rna_batch.py` and
`run_codex_native_rna_batch.sbatch`, with no edit to the successful pilot.
The wrapper permits exactly one of four ordered groups and reuses its
hash-guarded recipe functions. Each donor receives its own new output
directory, native aggregation, identical full42,163 CPM denominator and
raw frozen global/local heads. Only the donor identity and output cap vary.

| Fixed serial donor group | Frozen RNA raw bytes | Group saved-output cap, bytes |
|---|---:|---:|
| B7,B8,B9,B10 | 17,592,037,059 | 25,612,037,059 |
| B12,B15,B19,B21 | 19,745,789,824 | 27,765,789,824 |
| B22,B24,B26,B36 | 17,426,600,755 | 25,446,600,755 |
| B38,B41,B46,B47 | 18,774,126,559 | 26,794,126,559 |

These values were derived from the frozen f91c... RNA roster metadata
only. Each group has64 FASTQ files; per-donor cap is exact frozen raw
bytes plus2GB, and group cap adds a20MB receipt buffer. The four caps total
105,618,554,197 bytes; these caps do not themselves account for other
campaign branches or concurrent jobs. The launcher requests CPU8/64GB/90h
on cpu/nslab; root owns four individual submissions and campaign budgeting.
No array or submission was made by this agent.

Source/helper/adapter/weight/index/version/environment receipts and exact
commands are recorded before each donor's acquisition. Original source
hashes, fixed inputs and complete native index name/length checks are
repeated afterward. The exact17 paired identity intersection is required.
Every failed donor retains its error and partial outputs; there is one
attempt with no retry, substitution or omission. The other fixed group
members are attempted once after an ordinary per-donor failure. Fatal
group interruptions mark uncompleted entries as failed rather than
reporting completeness. Group success requires all four; cohort completion
always remains false pending reconciliation of every fixed donor.

Python static compilation, `bash -n` and metadata-only roster/hash/byte
checks passed. Batch execution **not run** by this agent; source and
launcher returned to root for review. Preparation read no raw RNA,
count/package array, receiving H3 or predictive loss, and added no fits,
bias flags, padding or biological-absence claim.

## Original P3 first draw: adaptive mode need not be a maximum

Direct metadata inspection of
`REC/p3_first_draw_diagnostic_21998696/summary.json` establishes that the
exception occurred during `fit_marginal`, at the third marginal call:
gene2, u=0, v=0.9234818604998617, reported curvature+763.1283114221253.
The executed diagnostic calls the fitter only after completing the
173-row original conditional generation. Its archived stack reaches
`likelihood_v1.py:423`, the nonconcave-mode exception; this resolves the
previous generation-versus-fitting uncertainty for the original failure.

The mode loop at snapshot412–423 starts at zero and computes
`step=s/max(-h,1e-8)`, accepting the first nondecreasing trial among only
30 halvings. It stops after an unsuccessful search, a small movement or
100 iterations without requiring small score. At h>0 the denominator is
1e-8, so its smallest trial displacement is approximately0.1862645*s.
It can skip the local ascent neighborhood and remain at zero. At a
stationary local minimum, s=0 yields an accepted zero movement and immediate
termination. A negative final curvature by itself also would not establish
stationarity or global optimality.

Positive curvature is mathematically possible under this exact model.
For genotype mixture `log(sum_c pi_c*f_c)`, the Hessian is the posterior
mean component Hessian **plus posterior score variance** (snapshot130–139).
The latter is positive semidefinite. The integer-ascertainment log
normalizer is then subtracted (175–186); this operation also gives no
global concavity guarantee. For the gene integrand
`h(u)=sum_rows log P(a | integer inclusion,u)-u²/(2v)`,
`h''=-1/v+sum_rows stage²*row_H_etaeta`. Mixture/selection curvature can
exceed the negative Gaussian-prior curvature. Thus real local convexity
is possible without an incorrect likelihood. Conversely, the original
shape arithmetic can corrupt derivatives at its enormous trial values.
The reported positive curvature does not by itself establish true
multimodality or exclude derivative cancellation.

The smallest useful next calculation freezes this third-call p/v and
original generated toy rows, evaluates only the gene2 scalar integrand
with stable BB arithmetic, checks score/curvature at zero against finite
differences and maps prespecified finite-u values. Sign-changing score
brackets can identify local extrema; fixed-interval numerical integration
and a Gaussian tail bound can diagnose whether adaptive20-node quadrature
is a faithful integral calculation. None of this requires another fit,
a new likelihood, alternate null, changed optimizer rule or model law.
Any later numerical mode-search/integration correction must retain the
same integrand and separately demonstrate integral accuracy; clipping
curvature or silently discarding an alternate mode is not that check.

Prepared new `codex_p3_gene2_integrand_diagnostic.py` and
`run_codex_p3_gene2_integrand_diagnostic.sbatch` (job name1591,
CPU1/2GB/1h/cpu/nslab), without changing the live likelihood or current
stable first-draw replay. Python static compilation and `bash -n` passed;
compute execution **not run** by this agent and launcher not submitted.

The helper hash-freezes the saved21998696 third-call p/v and its executed
source/receipts. Since generated count rows were not archived, it
reconstructs only original draw0 under the exact original source, observed
null parameter hash, supplied uniform/fixture hashes and recorded numerical
versions. It requires173 included rows and reproduced gene2 u0 curvature
within1e-6 of the archived value before proceeding. Reconstructed counts
receive a saved hash; absence of an independently archived count-vector
hash remains explicit. No truth labels, optimizer starts or new fits are
introduced.

Only the fixed gene2 integrand is then evaluated:481 points on [-12,12]
with spacing0.05; finite-difference steps1e-3,3e-4,1e-4 at zero; original
versus stable LL/score/H where finite; sign-changing score roots with
brentq absolute tolerance1e-10; scaled bounded integration with
epsabs1e-10/epsrel1e-8. All bracketed positive/negative curvature extrema
are recorded, without claiming that this grid proves all stationary
points have been found. Twenty-node integrations centered separately at
bracketed local maxima are reported as numerical comparisons, not a
production integration choice.

Each conditional integer count probability is at most one, so the
outside-[-12,12] **unnormalized** gene-integral tail is at most
`sqrt(2*pi*v)*2*NormalSF(12/sqrt(v))`. This absolute bound is reported
relative to the actual bounded integral; no small-tail conclusion is
assumed before that comparison. Multiple peaks, integration disagreement
or an insufficient tail bound would concern the numerical approximation
to the fixed likelihood, not authorize a changed likelihood or test law.

### Pre-submission diagnostic refinement after Claude source advisory

The above0.05-grid/fixed12-extent description is superseded for the new
helper. Root relayed Claude's tools-disabled numerical review, and only
this isolated diagnostic was amended. Current source SHA256:
`c643a77aba659d2882ee86f117f1ef02a7f13210dce871986153cf4763c6edd1`.
Python static compilation and `bash -n` passed again; no execution or
submission by this agent. CPU1/2GB/1h remains unchanged.

The helper now evaluates spacing0.01 followed by0.005, recording both
sign-bracket counts and exact zero-score grid points. Shared scalar
points are cached; tangential/subgrid extrema remain a stated limitation.
Because the proper conditional count pmfs are at most one,
`h(u)<=-u²/(2v)`. Hence no point beyond
`sqrt(2*v*abs(h(0)))` can improve on h(0). The prescribed extent is
`max(12,ceil(that_radius+2))`, with operational cap32. If the cap prevents
the radius-plus-buffer coverage, the global search extent is explicitly
unverified. The bounded integral and its Gaussian tail bound still remain
separate from a proof that every local extremum was found.

All30 original damping proposals **at the archived terminal u0** are
evaluated using both representations, with finite flags, acceptance
conditions, warnings and numerical exceptions saved. Old numerical
failures do not abort the stable diagnostic. Smaller symmetric probes are
at most1% of the smallest original displacement; when its displacement
is zero, fixed small symmetric probes are used instead. At zero and all
small probes, scalar hstats is checked against the original
`drawn_model.evaluate` gene2 likelihood minus the Gaussian penalty.

**Iteration-history qualification:** a terminal `mode=0` record does not
establish that the mode search never previously moved. The source shows
that premature termination at a nonmaximum is possible, and these frozen
zero-point proposals diagnose that possibility; they do not reconstruct
the missing iteration history. No production optimizer or integration
rule, likelihood, supplied draw, fitting criterion or current replay was
changed.

## Exact99 molecular training subset: trial/date allocation audit

**No sufficient2018-only or post-February2015 allocation evidence was
found for the exact99 paired training people.** This conclusion now
concerns the actual molecular subset, rather than inference from all108
participants described in the source paper.

Identity-only reading of
`Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture/molecular/participant_axis.tsv`
established99 participant records and198 distinct RNA/H3 GSM identities.
All198 join the frozen official SOFT at
`Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/raw/GSE267145_family.soft.gz`.
Their memberships are RNA subseriesGSE269412 and CUT&RUN subseriesGSE267119
within SuperSeriesGSE267145. Matched growth/treatment fields identify human
liver biopsy tissue obtained during bariatric surgery. This resolves the
paired molecular source as the study's bariatric-biopsy source, rather than
external reference-cell data or an unidentified control source; it does
not identify which registered enrollment protocol supplied each person.

The frozen participant axis has identity, source accession, pairing and
observation-state fields, with no trial/date column. The prior authoritative
participant_join.tsv header likewise has no trial/date column; its clinical
row values were not read. Across matched GEO records, characteristic field
**names only** concern tissue, sex, histology or antibody, with no trial or
specimen-collection-date field. Their May9/June7,2024 submission dates are
deposition metadata and cannot date biopsies. A2024 public sample date does
not establish that the tissue was collected after2015.

The exact Patients and sample collection section of the existing primary
preprint XML still identifies Vanderbilt IRBs090657 and171845 and both
NCT00983463/NCT03407833 without allocating individual IDs or the99 paired
subset to either trial. The published [JCI methods](https://insight.jci.org/articles/view/182232)
give the same unsplit registration/collection description. The confirmed
official start dates remain2009-10 and2018-02-13; no subset allocation
permits applying only the latter to every paired training donor.

Bounded supplement checks used only existing primary-source artifacts.
ZIP-entry metadata in
`gse267145-supplement3-21070402/raw/supplementaryFiles.zip` identifies
`media-11.pdf`,3,185,088 bytes,10 pages; automated trial/procurement-date
anchor extraction found no NCT identifier or collection/biopsy/enrollment/
recruitment-date anchor. Only the first-row string title of
`media-1.xlsx` was inspected: TableS1 Patient characteristics. No clinical
cells were inspected, so this title is not proof of absent dates elsewhere
in that workbook. The [current JCI supplemental PDF](https://insight.jci.org/articles/view/182232/sd/pdf/render/1)
lookup failed; its full contents are **not verified** here. One initial
overbroad XML-section title filter returned a public Results paragraph;
that incidental text was discarded and did not inform a method, model or
identity claim. Subsequent extraction was restricted to relevant metadata.

**Exact remaining gap:** a99-ID-to-enrollment-protocol crosswalk demonstrating
all99 belong to the2018-start study, a biopsy-date lower bound after
February2015 for each of those99, or authoritative person-overlap exclusion
against the fixed17 Penn donors. No such crosswalk/date field was identified
in the frozen subset metadata or inspected primary methods. These results
support the existing different-source transport qualification; they do not
authorize dropping people, changing weights or replacing identity evidence
with center differences. No molecular arrays/counts, participant clinical
values, actual genotypes, receiving outcomes, protected cohorts, contact,
fit, job or source-owner write occurred in this audit.

### Completed fixed-integrand evidence and isolated64-damping replay preparation

Root executed `codex_p3_gene2_integrand_21998880` successfully:9:39,
exit0. Its frozen third-call gene2 integrand has three local maxima and
two minima; both prespecified grid spacings.01/.005 found five sign
brackets. Original and stable BB atzero agree to approximately
1.6e-8 in curvature and both give positive curvature763.1283. Thus
this fixed-call local nonconcavity is real, rather than a stable-BB
arithmetic artifact. The overwhelmingly dominant maximum is
u=-0.1597680479632258, curvature=-734.4514244753362.
Bounded quadrature log integral=-86.38066008342989;20-node adaptive
quadrature centered there differs by-2.0836e-8. The bounded Gaussian
tail has relative log bound-54.3878. These observations concern one
fixed parameter call; no global mode guarantee or general quadrature
accuracy claim follows. Terminal modezero does not establish its prior
iteration history.

Prepared new source-only helper
`scripts/analysis/in_biopsy_allelic/rec_codex/codex_p3_private_damping64_replay.py`
and launcher `run_codex_p3_private_damping64_replay.sbatch` for root
review, reserved1591 CPU1/2GB/90h cpu/nslab. The helper copies the
immutable original P3 sources and changes exactly one private source
occurrence of mode damping30 to64, recording a reversible one-line
diff; stable BB remains the already checked numerical identity.
It reconstructs only supplied synthetic draw0, requiring173 rows,
24 gene2 rows and the completed quadrature diagnostic count-vector
hash. Before any fitting, it extracts and executes the exact private
zero-start mode-search AST, requiring finite stationary score,
negative curvature, mode agreement within1e-6 and20-node gene integral
agreement with the recorded full quadrature within1e-6. It then
compares original30 versus private64/stable observed LL, derivatives
and marginal LL/score at archived parameters before invoking the
unchanged existing observed-prefix/single-draw replay driver. The
nested driver descriptions refer to the already patched private
snapshot; the outer protocol and source diff disclose both numerical
changes.

The replay retains original observed warm/cold starts, bounds, ties,
20 nodes, conditional draw generation, statistic and failed-draw
exceedance treatment. It preserves the historical failed draw and
live replay21998667. No real allelic input, truth start, production
edit, owner specification amendment, MC p-value, full-fixture or
calibration claim is authorized. Python static compilation and
launcher `bash -n` passed; the extracted mode AST also compiled.
Executable mode, equivalence and replay checks are **not run** by this
preparation; no job was submitted.

### Distinct RNA-processing hypothesis: two-source admission decision

**Decision: no-go for an independently validated liver allele-to-splicing
model from these two inspected sources.** A distinct target exists: the
ALT-minus-REF change in an exact intron's share of its source-native
competing junction cluster. It is not total RNA abundance, chromatin
activity, or isoform count. However, the existing GTEx source is spent
development, and the separate long-read liver source supplies no matched
variant-effect validation target. This bounded check does not establish
that no suitable liver source exists anywhere.

Current directions checked before source evaluation: the five-authority
navigation and relevant identity/claim, model status, splice-results and
roadmap sections; integrated program §22.3.10–11; and
`config/contamination/exposure.toml` plus the AlphaGenome exposure audit
and `development_crosswalk.json` under
`Analysis/MASLD_Model_Benchmark/`. The integrated plan already proposes
exact event/competitor alignment; proposing a different neural architecture
does not resolve the missing endpoint or replication.

| Inspected source | Measured target and biological replication | Permitted role and decisive limitation |
|---|---|---|
| GTEx v8 liver | Bulk RNA splice-junction counts across human participants support intron-excision ratios and genotype-associated sQTL slopes. Participants, not junctions/variants/tissues, generated the measurements. Exact usable liver participant count was not recomputed from matrices in this review. | Existing source-exposed development only. Checkpoint registry declares GTEx training tracks; project signed-splice results and source-effect tables have already been examined. A new GTEx processing release or held genomic block supplies neither new participants nor untouched external-source confirmation. |
| Sarygina et al.2023 normal-liver ONT source, article DOI10.3390/ijms242115502; deposits PRJNA765908/PRJNA893571/PRJNA635536 | Methods §4.2 specifies three liver donors, one sample each, plus five HepG2 and three Huh7 samples. Direct-RNA MinION reads estimate isoform TPM; §4.3 defines isoform number above TPM thresholds. Cell-line samples do not increase liver participant n. | Potential structural/measurement reference only. Published methods do not provide paired participant genotypes, allele intervention, or signed variant-event effects. Three normal-liver participants cannot validate disease remodeling or a learned donor-context model. A liver-only deposit/sample crosswalk and correct assembly must be resolved before even a structural comparison. |

For GTEx, the official [LeafCutter pipeline documentation](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/leafcutter/README.md)
specifies WASP-aware alignment, junction extraction, intron-excision
ratios, phenotype-group mappings and explicit sample-to-participant mapping.
Its BED start/end fields locate TSSs; intron identity resides in the
phenotype identifier. Treating BED bounds as splice-junction coordinates
would therefore change the target. A molecular-QTL slope on transformed
phenotypes must retain its released transformation; it is not automatically
a raw difference in usage. The project source
`scripts/analysis/alphagenome_atlas/36_splice_direction.py:9` documents its
usage-share approximation and loose model-proposed competitor set.
`docs/RESULTS.md:2294` records the completed weak-direction comparison;
`docs/RESULTS.md:2741` records its outcome-selected source denominator and
unfinished exact-event evaluation. No effect values were read here.

The [ONT primary methods](https://www.mdpi.com/1422-0067/24/21/15502)
were read through the publisher article's Europe PMC full-text XML, limited
to method/data-availability sections. They specify minimap2 followed by
Salmon TPM estimates and thresholded isoform counts, so “direct sequencing”
does not make every published abundance an unambiguous direct molecule
count. Their reference description combines “Gencode38” with “GRCh40”;
that inconsistent wording is an unresolved coordinate contract. No
gene-by-sample table, supplementary abundance table, figures or raw payload
was opened. Article-listed deposits were not assigned individually to liver
versus cell lines without a verified sample crosswalk.

Bounded accession/DOI searches found no match for this ONT source in the
inspected project technical/literature documents, model configuration or
AG source scripts. **[Unverified]** That is not a complete historical
exposure audit or a checkpoint-training exclusion: it cannot establish
globally untouched participants. Neither an unknown registry state nor
absence of a filename warrants a clean validation label.

**Concrete next step, metadata only:** reconcile the existing GTEx tested
phenotype and group-file receipts with the official current download
manifest, preserving exact intron ID, strand, cluster, source assembly,
sample-to-participant map and effect transformation. The current official
[GTEx v8 download page](https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression.GTExv8)
has an indexed reference to a939-GB all-associations directory, whereas
the project's inspected bucket prefix was an empty placeholder. Actual
liver file availability/schema remain **[Unverified]** because the current
page requires JavaScript and no payload was acquired. Resolving this
metadata discrepancy could enable a fair all-tested **development** event
comparison; it cannot create external confirmation. No model fit, source
acquisition, job, protected-cohort access, allele-value inspection or
authority change occurred. An independent quantitative effect source is
still required before committing to the proposed new predictive capability.

#### Completed GTEx all-associations metadata reconciliation

**Result: no currently verified liver-only all-tested artifact.** The
historical empty-public-prefix finding remains true in a fresh metadata
request. The portal's all-associations description is a pointer to a
different requester-pays bucket; it is not evidence that the public bucket
contains an omitted liver file. That pointer is presently unusable through
the inspected official metadata API. This resolves the apparent conflict
without declaring an outcome endpoint impossible in principle.

Official current metadata checked without payload access:

- [Public v8 object-list metadata](https://storage.googleapis.com/storage/v1/b/adult-gtex/o?prefix=bulk-qtl%2Fv8%2Fsingle-tissue-cis-qtl%2FsQTL-Tissue-Specific-All-SNP-Gene-Associations%2F&maxResults=100&fields=items%28name%2Csize%2Cmd5Hash%2Cupdated%29%2CnextPageToken)
  returns HTTP200 and exactly the directory marker, size0 bytes,
  updated2023-09-19T14:41:56.998Z, with no continuation token or liver
  object. The parent listing also returns only this sQTL prefix plus the
  significant, phenotype, grouping and covariate archives.
- The5,327-byte [official description.tsv](https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/description.tsv)
  explicitly describes all tested cis-sQTL associations, including
  non-significant ones, but redirects to
  `https://console.cloud.google.com/storage/browser/gtex-resources`.
  It says approximately690GB; the indexed portal page says939GB. Neither
  is an exact liver-file size. The catalog's exact current public sizes
  are569,169,920 bytes for the significant sQTL archive,
  13,560,842,240 for phenotype matrices,47,041,697 for grouping,
  and2,236,078 for covariates. None is an all-associations replacement.
- Unauthenticated Google Storage JSON metadata requests to bucket
  `gtex-resources` return HTTP404, “The specified bucket does not exist,”
  for both `GTEx_Analysis_v8_sQTL_all_associations/Liver` and
  `GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_sQTL_all_associations/Liver`.
  No billing project, authentication changes or charged acquisition was
  attempted. The API response alone does not establish whether another
  migrated/archive route exists.
- The5,110-byte [official v8 README](https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/README_eQTL_v8.txt)
  documents `*.allpairs.txt.gz`/`*.sqtl_allpairs.txt.gz` fields:
  gene_id, variant_id, tss_distance, ma_samples, ma_count,
  pval_nominal, slope and slope_se. Variant ID encodes chromosome,
  first REF-base position, REF, ALT and b38. This verifies the intended
  nominal-association schema, but not any exact liver filename,
  available object, phenotype-ID implementation or allele-slope
  normalization in an actually retrieved all-tested artifact.

Reconciliation with existing development evidence is source-specific.
`week1-20260915/data/splicing-identity-21773601/receipt.json` and its
executed source name the inputs `Liver.v8.sqtl_signifpairs.txt.gz` and
`Liver.v8.sgenes.txt.gz`, plus historical step36 outputs. Its status is
`partial_source_identity_manifest_exact_effect_comparison_not_ready`.
The receipt input-hash metadata and source code were inspected, not their
effect rows. The existing `data_gtex_sqtl_tested_universe.sbatch`
separately acquires phenotype/group/covariate archives; those identify
measured events but do not contain all tested variant-effect estimates.
An event-denominator fix alone therefore cannot manufacture a complete
variant-event association denominator from significant effects.

**Stop condition met:** exact authoritative liver all-tested URL/name,
public access and disk bytes are **not verified**. Do not download the
939GB directory or substitute significant pairs as all-tested negatives.
A corrected official object manifest is the remaining acquisition
prerequisite. Any subsequently admitted GTEx comparison remains exposed
development, with exact-event and transformation checks required; it
does not supply independent external confirmation. No compressed
payload, association/molecular row, raw read or genotype was streamed;
no job or charged access occurred.

#### Other human-liver splice sources: bounded independent check

**No newly verified signed variant-event validation source.** One primary
study was inspected, and a second candidate was stopped at public title/
author/sample-count metadata because its source-family boundary overlaps
the coordinator's explicit Wenz exclusion. No further source search is
needed for this bounded task.

[Zhou et al.2022, UNC Liver](https://www.nature.com/articles/s41598-022-18506-z)
has192 paired genotype/RNA-seq liver samples after QC. Its
Splicing-QTLs methods use sQTLseekeR: relative exon-expression proportions
are a multivariate gene phenotype, with resampling-based gene-level
association p-values and BH correction. This genuinely measures genetic
association with splicing, but the published test is unsigned and not
a signed slope for one identified intron. Availability names RNA GEO
GSE179250 and genotype GSE26105. A raw RNA/genotype association could
potentially estimate a signed event effect under a new fixed analysis,
but that would require permitted participant identities, genotype
access, exact event counts and a new association fit; none was admitted
or performed here. A gene-level p-value cannot be relabeled as that
effect, and no public all-tested signed event export was verified.

Importantly, integrated program line1813 already identifies Currin's
138-column GSE26105 accession map. UNC's same genotype deposit establishes
a shared source identifier, not a numerical participant-overlap count.
**[Unverified]** UNC192-to-Currin138 person overlap needs an authoritative
identity crosswalk before a source-independent evaluation claim. It
cannot be assumed absent because the papers have different dates or
the assays have different names. The UNC article's comparisons with
GTEx also do not establish genotype donor disjointness for this project.

The second surfaced candidate is Gawronski et al., *Evaluating the
Contribution of Cell Type-Specific Alternative Splicing to Variation in
Lipid Levels* ([primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10284136/)).
Its public metadata describes96 primary livers and83 paired iPSC/HLC
lines, and its author roster includes Wenz/Brown. Only those search
metadata were inspected; methods, supplements, outcomes and payloads
were not opened. It may supply an appropriate hepatic event endpoint,
but eligibility and source-family overlap are **[Unverified]**. Under the
present explicit Wenz exclusion this is a stopped route, not a source
admitted for model selection or a claim that such data are unusable.

Exposure distinctions remain separate: (1) this methods/metadata review
reveals no fitting outcome or numerical allele target; public paper
awareness alone does not make participants statistically nonindependent;
(2) actual participant overlap is unresolved, with a concrete shared
deposit for UNC/Currin; (3) checkpoint pretraining membership for UNC and
the stopped second source is unknown, not established leakage;
(4) reuse of existing Currin labels/GTEx evaluations is already recorded
development exposure and cannot be called untouched confirmation.
No protected cohort, genotype directory, allele-count/value table,
molecular array, acquisition, contact, fit or job was used. Stop until
a compatible signed event export and permitted identity/exposure
contract are available.

**Clarification and canonical citation correction, following coordinator
scope clarification:** author Wenz is not itself a protected-source rule.
Public methods/data-availability for the second paper were explicitly
authorized and then inspected through official NCBI efetch XML
(`db=pmc&id=10284136`), without extracting Results, figures, gene examples
or molecular tables. Its canonical publication is Gawronski et al.2023,
[DOI10.1161/CIRCGEN.120.003249](https://doi.org/10.1161/CIRCGEN.120.003249).
The preceding2019 wording came from an earlier related/preprint citation
and is superseded for this primary article.

The main Methods identifies source data reused from Pashos et al.2017
(DOI10.1016/j.stem.2017.03.017), Penn-approved line generation with
participant consent, and full methods in supplemental
`hcg-16-248-s001.pdf`. It directs anonymized summary results to
DOI10.5281/zenodo.7864256 and other supporting data to author request.
No author contact was made. The main methods do not give exact signed
intron-effect units, an all-tested variant-event denominator, or the
primary-liver deposit/participant identities. The Zenodo record metadata
API returns404 in this check; no summary-result payload was opened.
The publisher supplement HEAD request also returns403; its full
methods were not retrieved or inferred from the article title.

The shared author roster neither proves leakage nor proves a shared
participant. The named source ancestry is now explicit, but exact
primary-liver overlap with existing Currin or reserved Wenz people and
checkpoint inclusion remain **[Unverified]**. It is scientifically
plausible to use permitted source summaries for a new fixed development
task after endpoint/identity admission; metadata awareness alone does
not foreclose that. This review still establishes no signed all-tested
primary-liver export or clean confirmation cohort. Protection remains
on reserved outcomes and local data, with no protected source opened.

### Private synthetic derivative-cost probe prepared, no execution

New owned helper `codex_sparse_gene_derivative_probe.py` and launcher
`run_codex_sparse_gene_derivative_probe.sbatch` under
`scripts/analysis/in_biopsy_allelic/rec_codex/` prepare a bounded comparison
of identical conditional likelihood derivatives. Reserved1582 allocation:
CPU1/2GB/1h cpu/nslab. Python static compilation and `bash -n` passed;
numeric comparisons/timing are **not run**, and no job was submitted.
Root will review the algebra/source before any execution.

For each row, the active coordinates are its own alpha_g and r_g plus
kappa, every retained delta, omega_S and the three rho coordinates.
The original full parameter dimension is2G+Z+5; this local derivative
dimension isZ+7, independent of G. With row conditional score s and
Hessian H in (eta,r), the helper accumulates
`J_active.T @ s` and
`J_active.T @ H @ J_active + s_eta * H_eta_active` into the corresponding
global coordinates. The only nonzero eta second derivatives connect
alpha_g to kappa (stage) and delta (z). r is linear. The helper keeps
the same row order, eta/r formulas, integer-selection/mixture likelihood
and global dense final Hessian; it removes per-row full-width Jacobian
and eta-Hessian construction. No predictor/model equation, bound,
fit start, integration or statistical law is changed.

Compute execution refuses without SLURM and guards the immutable original
likelihood SHA41ebc253… plus all four original companion sources and
the published/original-receipt hashes for globals, participants, cohorts
and synthetic rows. It uses the archived observed-null p, three fixed
seed20260930 moderate perturbations (sd.03, within original bounds),
u=0,-.05,+.05 and score-only/full-Hessian modes. Widths8,64,379 keep
the same173 supplied toy rows and merely add structurally unobserved
gene parameters; no additional donor, row or biological gene is added.
Per-point LL/gradient/Hessian agreement requires atol1e-8/rtol1e-10,
finite outputs and Hessian symmetry. Immutable source/input hashes,
paired evaluator order, instrumented timing, traced peak allocation
and process maximum RSS are saved with environment/version information.
The existing hom-inclusion cache is cleared equally before each role.

This addresses a concrete allocation/arithmetic cost, not the full
379-gene integrand or actual fit workload. The output Hessian remains
quadratic in parameter dimension, BB/truncation arithmetic remains
unchanged, and tracemalloc-instrumented timings are qualified. Root's
separate completed cache benchmark21999131 (reported27s,1.12–1.19x)
does not establish full-fit runtime savings. No protected/real allele
value, receiving outcome, truth genotype, optimizer, cold start,
MC draw, canonical source edit or overlapping cost-helper edit occurred.

#### Executed derivative probe21999177: arithmetic agreement and bounded cost

Root executed reviewed helper SHA467a68f797365e4a70ca37a5be6f29324e80e9d0f181f114a0ecbd1211d1781e
as21999177/1582, completed1:22, exit0, empty stderr. Direct review of
`codex_sparse_gene_derivative_21999177/summary.json` and `points.json`
confirms all72 comparisons passed, with **exactly zero** LL/gradient/
Hessian differences. At width379, full-Hessian dense/local instrumented
time ratios range2.591748–2.781044; maximum traced allocation is
19,024,196 versus4,925,607 bytes. Resource process maximum RSS is
102,328KiB; root reports sampled SLURM57,080KiB, a different measurement
that does not override the process maximum. Width8 is nearly neutral,
and score-only evaluations are neutral or approximately1% slower.

This supports the local-coordinate arithmetic and its full-Hessian
allocation benefit on the supplied173-row dimension probe. It does
not measure an actual379-gene fit, integration runtime, full calibration
or biological data. The preceding “not run” describes preparation only;
no additional checks, replay or code patch was performed in this review.

### Independent interpretation of private64 P3 replay21998950

Direct metadata review of `codex_p3_private_damping64_21998950/` confirms
the reviewed private source diff changes only the gene-mode damping
proposal range30 to64. Stable BB is injected in the private runtime;
original sources and the original failed-exceedance record remain intact.
Companion hashes retain the original starts, parameter bounds,
variance-boundary logic,20 quadrature nodes, supplied uniform ordering,
LR definition, tie tolerance and failure-as-exceedance rule.

The archived fixed gene2 call now finds mode-0.15976804796323024,
score-3.54e-12 and curvature-734.4514244753359. Mode error against the
completed full-integrand diagnostic is4.44e-15; the20-node log integral
differs from bounded full quadrature by2.0834e-8. Original30 versus
private64/stable observed LL/gradient/Hessian differences are
5.59e-11/9.19e-10/1.25e-9; marginal score difference is1.94e-8.
Observed refit differences are2.79e-11 in LL,3.72e-11 in kappa and0
in v. Reconstructed and actually replayed173-count vector hashes both
equal5f4a1aa455cd09b0d5f31889021ec88b145dcbb50c8e5b600343ea9a3a992691.

Only supplied draw0 was refitted. Warm/cold starts both converge and
their LLs differ by approximately1.39e-10; warm is retained under the
unchanged1e-6 selection rule. The kept fit hasv=0, LR=0, and observed
LR is also0, so the unchanged tie rule produces exceed=true. The
historical failed draw also had exceed=true through its failure rule.
Thus this draw's numerical fit state was repaired while its individual
exceedance contribution stayed the same. This is not evidence about
nonzero observed LR comparisons or an MC p-value.

**Supported:** the larger damping range resolves this actual frozen
nonconcave transit call, reproduces its dominant integral accurately,
preserves observed-fit numerical agreement and allows one previously
failed synthetic fit to finish. Stable BB alone did not resolve the
mode-search failure; the mathematical likelihood and statistical rules
were not replaced. **Limits:** both kept observed/draw fits end atv=0,
so final-fit equivalence chiefly checks the boundary likelihood. It
does not validate adaptive quadrature at a general positive-v optimum,
the marginal derivative accuracy there, near-comparable multiple peaks,
failure rates over draws or biological/full-size designs. Local
stationarity and successful warm/cold convergence are not a global
maximum proof. Moving the numerical center can change a finite20-node
approximation even when the intended integrated likelihood is unchanged.

Before a broader proposal, the owning analyst must explicitly select
and freeze the numerical recipe and preserve this historical diagnostic.
The next numerical evidence needs an approved fixed positive-v case:
agreement of the unchanged20-node marginal LL/score with independent
reference integration/finite differences, including a case where mode
competition matters. The existing complete synthetic fixture/test
obligations then still apply to the frozen recipe and original supplied
draws. Only after numerical closure and a fully specified pre-ascertainment
generator may the approved null/coverage calibration address the stage
CRT law and intended biological design. None of those broader results
is supplied by draw0 or this review; no new job, fit, test, protected
input or canonical edit was performed.

### Independent development-program deviance verification prepared

Prepared `scripts/analysis/in_biopsy_allelic/rec_codex/codex_program_deviance_verify.py`
and `run_codex_program_deviance_verify.sbatch` for one compute-only
verification of completed job21999301. The helper independently reconstructs
each of99 development participants' multinomial deviance as
`2*N*sum(q[q>0]*log(q[q>0]/p[q>0]))`, with `q=C/N`, from the
saved float64 program/PCA OOF probabilities and fixed development H3 counts.
It reconstructs training-fold pooled and macro mean probability baselines,
using the producer's declared additive probability floor1e-12 and
renormalization; it imports no producer functions and performs no fitting.
It compares all99 saved donor metrics and the paired10000-draw bootstrap,
seed20260930 with linear quantiles. Exact participant/fold/region order,
positive finite normalized probabilities, all five folds, source receipts
and before/after hashes are required. No receiving arrays, genotypes,
tag counts or new candidate are used.

The parent reports job21999301's negative development result:
mean paired skill difference−0.12300622465526184, conditional interval
[−0.15940105752450365,−0.08801922029582603], program pooled-baseline
skill−0.0001299731053157767 and PCA skill0.12287625154994604.
These are reported producer results; independent numerical verification
is **not run** at this preparation stage. The interval remains conditional
on the saved overlapping-fold predictions and this development source;
it does not include retraining or transport/assay uncertainty.

Initial source compilation and launcher `bash -n` passed without importing
the analytical helper or executing numerical calculations on login.
After the parent's surgical job-log and uncertainty wording edits,
final prepared helper SHA256 is
`a5447c4d05e7bd4a8b4391cef59c43c9aea4701934e1b4f085a6484029039a4f`
and launcher SHA256 is
`13d4a4aed17dd1b2aa7eb010be86307c1abc578330425b2f1752ce5cd8e63a06`.
The launcher requests one job, name1585, CPU1/2GB/1h, cpu/nslab.
Submission and final static checks belong to the parent. No job was
submitted by this reviewer and no scientific criterion was changed.

The parent subsequently submitted the single verification as21999323:
completed13s, exit0, empty stderr. All99 independently reconstructed
deviances, skills and seeded bootstrap agreed. Largest deviance difference
was1.1641532182693481e-9 and largest skill difference1.11e-15;
independent primary difference−0.1230062246552619 with interval
[−0.15940105752450373,−0.08801922029582608]. Process peak memory
394280KiB; sampled SLURM memory0 is unusable. This confirms the negative
development result and failure of the declared criterion. No new fit,
tuning, receiving evaluation or further numerical test was performed.

### Six fixed positive-variance contexts: prepared diagnostic, not run

Prepared `scripts/analysis/in_biopsy_allelic/rec_codex/codex_positive_variance_grid.py`
and `run_codex_positive_variance_grid.sbatch`. This extends the previous
gene2/single-point check to both archived parameter vectors (observed-null
and failed positive-variance proposal), each at fixed v=0.05,0.5,2,
using the exact173 reconstructed synthetic draw0 rows and all8 genes.
For each context/gene, it compares adaptive GH20/40/80 mass and all
parameter scores plus variance score against independently integrated
normalized-Gaussian-prior mass/score numerators. Full-model production
entry-point results are also compared with the assembled gene results.
Failures remain recorded; a reference center cannot replace a production
mode failure. The immutable owner20-node rule remains unchanged.

The reference identities are `score_p=E[conditional_score_p|rows]` and
`score_v=E[(u²-v)/(2v²)|rows]`. Prespecified finite differences independently
check local alpha/r, kappa and variance for genes0 and2, absolute steps
3e-4 and1e-4, across all six contexts. Coarse/tight and expanded-bound
vector quadrature are retained. The proper conditional PMF≤1 supplies
Gaussian mass-tail bounds; p-score bounds and Gaussian second-moment
tails qualify score truncation. QUADPACK errors and the inferred lower
mass are numerical estimates. The0.05 grid supplies integration points,
without guaranteeing that all narrow modes or a global maximum were found.
Requested bounds beyond64 are explicitly flagged and cannot imply verified
accuracy. No fits, MC, model selection or calibration are performed.

The private copied source changes exactly one damping occurrence30→64,
plus the previously verified stable-BB arithmetic injection. Original/private
AST identity after reversal, extracted-mode AST, node variants, all source
and fixture hashes, exact reconstructed-count hash, archived parameter
sources, environment and before/after guards are recorded. Each gene/case
and its failures are saved. The launcher requests one CPU/2GB/4h,
cpu/nslab, name1580 with actual-job-ID logs. [Estimate]1–3h runtime uses
the previous21998909 measured792.271s/24589 evaluations as a cost basis.

Preparation source SHA256:
`3a7fe694d2de2a6bbc1067036995cb8ae4c0b6d2767d20325389e059bf36af85`;
launcher SHA256:
`b1bac3ea59dd3a18950605926554ceee089e4425e86bc157cb91e41b0cbddfe3`.
Static compilation and `bash -n` passed. Analytical imports, numerical
checks and submission were **not run** by this reviewer. Parent review
precedes submission; these diagnostics cannot establish null size,
coverage, biological-design accuracy or adoption of a numerical recipe.

Before execution, the parent requested an independent proof of the
parameter-score tail envelope. The existing bound is valid for this
frozen formula in exact real arithmetic. Define rising-factorial
log-shape derivatives `D_a=sum_{i<k}a/(a+i)`, `D_b=sum_{i<n-k}b/(b+i)`
and `D_t=sum_{i<n}t/(t+i)`. They lie respectively in[0,k],[0,n-k],[0,n].
The BB eta score `(1-m)*D_a-m*D_b` is bounded in absolute value by n.
The r score is `(D_a+D_b-D_t)*d(log t)/dr`. With rho=e+c*l,
c=1−2e, l=logistic(r), x=l(1−l)≤1/4,
`rho*(1-rho)-rho' = e*(1-e)-2*e*c*x >= e/2 >0`.
Consequently `|d(log t)/dr|=rho'/[rho*(1-rho)]<=1`, so the
r score also has absolute bound n. Mixture prior weights and homozygote
PMFs are fixed in the executed formula: mixture score equals posterior
heterozygote weight times BB score. Selection subtracts the allowed-count
mean mixture score, giving bound2n per predictor. Since the predictor
Jacobian is independent of u, each parameter component has uniform
bound `sum_rows 2*n*(abs(J_eta)+abs(J_r))`. Thus its Gaussian tail
numerator envelope follows; conditioning does not invalidate it.

Added this derivation beside `uniform_p_score_bound` and to the saved
protocol, and guarded matching EPS in stable/original arithmetic. This
bound concerns the exact mathematical score; it does not certify floating
point error. The posterior-tail estimate still uses an estimated numerical
lower integral, explicitly qualified. No equation or bound value changed.
Updated helper SHA256:
`dc49d1d807c41e2b81f65ce4c029bd01f278d3a169f99ab6dd44fd21aa462de1`.
Launcher SHA256 remains
`b1bac3ea59dd3a18950605926554ceee089e4425e86bc157cb91e41b0cbddfe3`.
Static compilation passed; no numerical execution or submission.

### C1 ASEP application biological n: resolved from primary methods

The original `2026-09-29_codex_prior_art.md` ASEP row left macrophage
application n unverified. The primary [Fan et al.2020 article](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008786),
“Application to a human macrophage RNA-seq dataset,” states48 healthy
individuals with paired M0/M1 preparations. Methods, “Human macrophage
differentiation and polarization and RNA sequencing,” explicitly reports
96 samples from48 subjects. Thus biological n=48 participants, with two
condition samples per participant; cells, SNPs, genes and reads do not
increase this n.

Methods, “Detection of differential ASE between two conditions,” uses
the same individual in both conditions, subject-specific random
intercept/slope, and a gene-level condition contrast. Participant is the
replication unit and gene the tested hypothesis unit. Gene-specific
eligible participant counts depend on filtering and are not established
by the cohort total. Article Data Availability links macrophage allelic
counts to [Dryad10.5061/dryad.866t1g1nb](https://datadryad.org/dataset/doi%3A10.5061/dryad.866t1g1nb);
no deposited count rows were opened. This resolves the C1 cohort-n
uncertainty without changing its novelty conclusion. No local protected
data, outcome analyses, fits or jobs were accessed.

### Original synthetic P1 draws24/32: finite fits, retained failure flags

Read saved `bfix_tests_21985999/p1_024.json` and `p1_032.json` under
REC only; no derivative or fit recomputation. Both store finite parameter,
likelihood, score and Hessian records, selected warm `converged=false`,
and cold `converged=true`. Both retain `z=NaN,exceed=true`. Neither
record contains an optimizer message, exception record or polish history.
[Unverified] Iteration exhaustion, line-search stopping and arithmetic
causes cannot be distinguished from these records.

|Zero-based draw|Warm/cold likelihoods|Interior scores establishing failure|Active-bound scores|
|---|---|---|---|
|24|−533.01604200377 / −533.0160420037774|omega_S index24 +5.749339058003233e-6; kappa index16 −5.139877108002039e-6|r indices10/12 at−15, scores−0.0010689984233175414 / −5.2075416233825014e-5|
|32|−532.7602192038512 / −532.7602192038596|omega_S index24 +1.9114726254432224e-5; kappa index16 −7.985616148332292e-6|r indices10/12 at−15, scores−0.0015193848150381875 / −0.00011884219845706724|

Under the documented strict projected-score tolerance1e-6, these
interior coordinates establish nonstationarity. Negative r scores at
the lower bound point outward and project to zero; they are not the
failure mechanism. Large delta coordinates have no finite bounds in
the documented model. Cold's convergence flag cannot override the
prescribed likelihood selection: neither cold likelihood exceeds warm
by1e-6. These are the same documented retained-warm/interior-score
mechanism as the draw24 review and original draw9 recovery context,
not P3's nonconcave marginal-mode exception.

`2026-09-30_codex_bfix_draw24_p2_draw3_review.md` already documented24.
Original live P1 loaded its fitter before later constrained-polish and
allocation changes and has no executed-source snapshot; present source
hashes must not be attributed to that process. Once original execution
finishes, the necessary bounded next check is the existing draw2/9
procedure on these exact supplied stage permutations: both prescribed
starts, preserved bounds/tie/tolerance and fixed numerical snapshot,
then independent LL/score/projected stationarity and participant-cluster
sandwich/statistic verification. Preserve original failed exceedances
and report corrections separately, counting each draw once. No pass
reclassification is supported now; no refit was authorized/performed
by this inspection.

Draw24 SHA256 remains
`7a6ea53579890bdcc2a9d139ada3fa7b7ff0fc7d10c612d4da0db97ec5f9f79f`;
draw32 SHA256
`a18ed503bcd31aed52bd6bbd7a7cb5aba914e6566b4a83cb4b55d5ed340b88d7`.
Active grid sources and jobs were not changed.

### Partial grid21999349: observed-null vector at v0.05 only

Read saved eight `observed_null_v0.05_gene{0..7}.json` receipts, not
numerical arrays or a recomputed fit. No production/reference failure,
cap truncation or quadrature warning is recorded in these cases. All
mode-stationarity flags within1e-5 are true; gene2's recorded mode score
is−1.2795343651283986e-7, with others smaller. Each reference grid
records one peak; this is not a proof of unimodality.

The largest visible saved AGH20/40/80 absolute logmass discrepancy is
1.821831574488897e-11 (gene2,80nodes); largest20-node score discrepancy
is−5.618350229497082e-10 (gene4 variance). No large node-count discrepancy
is evident at this particular parameter/variance context. Tight/expanded
reference values are close, but QUADPACK errors are estimates. Gene4's
base mass-tail ratio is5.647786393681489e-12 and variance-score tail
envelope4.095137954948967e-9; expansion reduces the latter to
4.305360930697837e-17. Tail qualifications therefore still matter even
when bounded integration appears accurate.

Finite-difference variance disagreements shrink by approximately9fold
when the step changes3e-4→1e-4: gene0−5.38488e-5→−5.98313e-6,
gene2−1.05595e-4→−1.17326e-5. [Observed pattern, not a guarantee]
This is consistent with central-difference truncation error. Gene2 alpha
also improves1.25150e-5→1.37135e-6. In contrast, gene2 r (index10)
disagreements−1.9533993718119473e-7 and−2.0713494659470397e-7
do not improve and exceed their saved quadrature/tail error estimates
2.718190194988606e-9 and8.150177715355239e-9. That direction remains
unresolved; the reported error estimates do not cover floating-point
mass cancellation. All node-score estimates agreeing with each other
does not independently settle the score-versus-mass-derivative issue.

These partial findings support completing the fixed grid. They do not
justify adopting40nodes, changing the model, fitting a new candidate,
or claiming48-case/full-model/calibration success. If the r finite-
difference pattern persists after completion, a focused arithmetic/FD
precision diagnosis is needed before a general numerical-accuracy claim.
No active source, job or statistical rule was changed by this review.

### Additional partial21999349 receipts: observed-null v0.5 and v2 gene0

Reviewed only newly saved `observed_null_v0.5_gene{0..7}.json` and
`observed_null_v2.0_gene0.json`. No production/reference failure,
nonstationarity flag, cap truncation or quadrature warning is recorded.
At v0.5, the largest absolute producer-saved discrepancies are:

|Nodes|Logmass discrepancy|Any parameter/variance score discrepancy|
|---|---|---|
|20|4.840412515250136e-10, gene4|8.217005964539226e-9, gene4 alpha|
|40|9.990230864787009e-12, gene2|1.152282713690056e-10, gene4 alpha|
|80|6.195932655828074e-12, gene3|4.307090166293104e-11, gene3 kappa|

Several v0.5 grids record two or three peaks, although production
stationarity flags all pass. A local stationary center is not a global
maximum proof. Saved tight/expanded reference masses nevertheless agree
closely in these cases; their numerical error estimates remain qualified.
Tight relative-mass errors range approximately1.4e-12–4.1e-11. Gene4's
base prior mass-tail ratio1.8154557151706412e-12 falls to
7.534682713703443e-21 under expansion. No reference result establishes
general quadrature accuracy or unseen-mode exclusion.

Gene2's r finite-difference issue persists at v0.5: discrepancies
−2.1516418490018317e-7→−1.658525190383238e-7 at steps3e-4→1e-4,
exceeding reported quadrature/tail estimates4.20777e-9/1.27219e-8.
Gene2 kappa changes+2.094349690548447e-8→−1.0645681561283205e-7,
also worsening rather than exhibiting clean step refinement. Its alpha
and variance disagreements improve approximately as step squared.
[Unverified] Whether the nonrefining directions reflect mass arithmetic
precision or another cause remains unresolved; integration error estimates
exclude floating-point cancellation.

For v2 gene0 alone, all-node absolute logmass discrepancies are at most
2.5579538487363607e-13 and score discrepancies at most
9.87446235889422e-13. Its kappa/variance finite differences fall within
saved error estimates; alpha shows the expected refinement pattern.
Other v2 genes and failed-proposal contexts are outside this review.
The next action remains completing the fixed grid, then resolving gene2
mass/FD precision before a general accuracy claim. No node-count adoption,
new fit, restart, submission or active-source edit is justified/performed.

### Observed-null full-model fixed-point likelihood comparisons

Read newly completed v2 genes3–7 and the saved three observed-null
full-model entries in `codex_positive_variance_grid_21999349/contexts_so_far.json`.
All three20-node fixed-point likelihoods meet the written absolute1e-6
agreement tolerance against the bounded reference, subject to its numerical
qualifications. Exact saved likelihood pairs are:

|Fixed variance|20-node production entry point|Bounded reference|
|---|---|---|
|0.05|−529.9766862529597|−529.9766862529452|
|0.5|−538.1097198776027|−538.1097198771627|
|2|−543.5480366323939|−543.5480366316416|

Production entry-point minus assembled-gene likelihood is exactly0 for
all three20-node entries; no entry-point failures are saved. Across
all saved20/40/80 entries, the sole nonzero assembly likelihood difference
is−1.1368683772161603e-13 (40nodes,v2); largest saved assembly score
difference has absolute value1.5010215292932116e-13. This establishes
that the helper's gene assembly reproduces the private entry point at
these contexts, not that its derivatives are exact.

Newly reviewed v2 genes3–7 record no production/reference failure,
nonstationarity flag, cap truncation or warnings. Each reference grid
records three peaks. Largest newly saved20-node gene logmass difference
is−7.625562403745789e-10, gene4; score difference−1.3055889702684453e-8,
gene4 alpha. Its40-node alpha discrepancy is−2.886766381493544e-10;
largest new80-node score discrepancy is+6.376424391209452e-11, gene3
kappa. Tight/expanded reference masses agree closely. Gene4's mass-tail
ratio decreases9.743211315213987e-13→3.0051140279646298e-21 with
expansion. Finite quadrature, estimated errors and local-mode limitations
still apply. The previously reported gene2 FD precision issue is not
resolved by full-model likelihood agreement.

No failed-proposal gene receipt was saved at this inspection. The next
action remains finishing the other24 fixed gene/context cases. These
three likelihood checks are no fitted-model, all-grid, null-size or
calibration result. No numerical workload, source edit or submission
was performed by this review.

### Newly saved original P2 draw14 and P3 draws12/13

Source/JSON inspection only, with paths relative to REC. P2
`bfix_p2_21998518/p2_014.json` stores a finite retained cold fit:
LL−542.3888976213885 versus converged warm−543.6254051455136.
Cold's higher likelihood satisfies the prescribed≥1e-6 selection
rule despite its `converged=false`. It records interior alpha(G3)
score−5.876890598122131e-5, omega_S−1.7278535798026518e-5 and
kappa+8.164846234434697e-6, exceeding the strict1e-6 projected
stationarity tolerance. Scores at r(G5)/r(G8)=−15 are negative;
rho_b/rho_d=+3 scores+0.13162036015271628/+0.03345629015791503
point outward. These active-bound scores project to zero. Large
unbounded delta values do not violate bounds. `kappa=NaN` retains
the failure, while both supplied lambda values remain saved. No
optimizer termination message is stored. This is a finite-fit/KKT
failure, distinct from a missing fit; the saved fitter already includes
constrained polish, so an omitted polish call cannot explain it.

P3 `bfix_p3_21998519/p3_012.json` and `p3_013.json` both store
`fit=null,LR=NaN,exceed=true`. The immutable executed runner initializes
fit to None, catches ValueError/FloatingPointError/LinAlgError across
generation and the full marginal fitting call, and does not store
the caught exception. These records prove failure before a complete
fit was returned. [Unverified] They do not identify generation,
optimizer, arithmetic or integration as the failing phase, nor prove
the previously diagnosed draw0 nonconcave-mode cause. No finite
point exists here for a KKT assessment. Exception instrumentation
of these exact supplied synthetic draws is needed to distinguish
causes before any corrective fit or passing classification.

Both immutable likelihood snapshots hash to
`41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6`.
Draw14 SHA256:
`35ee0a204f8ecd529a20028a7606d09cbc118cad46034d4dfbc98c902cb531a4`;
P3 draw12:
`92e214b618f3e03e626cba50c48c8fa1e551de3f5f9f0dd07a6b497b81d20bf0`;
P3 draw13:
`85767c9d280f255fc653ae7baf262ace5ced06d9943557a1852001493c84f284`.
No original file, active source or job was changed; no numerical
calculation, fit or test was executed.

### Completed failed-proposal v0.05 context: score and mode qualifications

Read its completed eight-gene/full-context record in
`codex_positive_variance_grid_21999349/contexts_so_far.json`.
Twenty-node whole-model LL−547.4559337835283 versus bounded reference
−547.4559337637646 differs approximately−1.976e-8, satisfying the
written1e-6 likelihood agreement at this fixed point. Forty/eighty-node
entry-point LLs are−547.4559337637451/−547.4559337637534; all three
entry-point minus assembled-gene LL differences are exactly0.

The consequential gene2 twenty-node discrepancy is omega_S score
−1.1560067330051993e-6 versus reference, with LL discrepancy
−1.974714791685983e-8. Corresponding forty/eighty-node score
discrepancies are−2.7671864799572177e-10/−1.627706858187139e-10.
Thus likelihood agreement alone does not establish equivalent score
accuracy. The reference shares the conditional evaluator, so it is
an independent integration check, not an independent complete
likelihood or calibration verification.

Gene5 records mode score+1.0070229426929755e-5, giving a false
`stationarity_within_1e_5` diagnostic flag for all node counts. No
production exception occurred; its LL/score integrals are close to
reference. Preserve the flag rather than calling every mode stationary.
All eight references succeed with no warnings/cap truncation; grids
record one peak each. Numerical errors and inferred tail-to-mass
bounds remain qualified, not rigorous global-mode guarantees.

Gene2 r finite-difference discrepancies worsen−9.277742583219117e-8
to+3.965733562857898e-7 as steps shrink3e-4→1e-4, exceeding reported
quadrature/tail estimates2.785019252824462e-9/8.5033111002745e-9.
Alpha/variance errors refine approximately as step squared. This
continues the unresolved finite-difference precision issue. Finish
remaining fixed contexts before any numerical-recipe decision;
no higher-node adoption, new model, fit or calibration claim follows.
Original/active files and jobs were not changed.

### Source review of isolated P3 draws12/13 exception diagnostic

Reviewed completed `diagnose_p3_draws12_13.py` SHA256
`c7a0792a2c59fc0986c25a170f903185cb3c91d1cdac37d1c0fbdf58d659f7ce`
and launcher SHA256
`9c5f489b3a29023a10d877cb522664c5c50312d532a2290723b3aed58b8cc8af`,
without imports, numerical execution or edits. Archived likelihood41eb,
pinned role receipt/null/draw records, original173 row identity/order,
supplied uniforms[12]/[13], saved-null parameters, tau_max_fixture and
`section4_from_null=True` match the original runner. Full50 completion
checks precede numerical imports; launcher also has afterok dependency.
Source/import/setup/hash failures are fatal. Per-draw capture preserves
the original exception classes, adds phase/traceback/frame records and
does not mutate original inputs or define corrected statistics.

One material exact-replay prerequisite was sent to the parent: the
helper records its numerical environment but does not compare it with
the original completed `result.json` Python/NumPy/SciPy fields. Require
equality after imports and before generation. Otherwise a changed
environment could produce a different failure while being described
as the exact original numerical diagnostic. No other decision-changing
source defect was identified in this bounded review. The author/parent
owns any correction; original and active sources remain untouched here.

The parent subsequently added refusal on Python/NumPy/SciPy mismatch
against the completed original result before generation, at phase
`original_environment_comparison`. Final reported helper SHA256 is
`fdfdf55e63225135f06eea45dc8e6dc1189033aab2189cbba03fb8e206e5c63a`;
launcher unchanged. Parent AST/shell/metadata inspection passed; no
additional numerical check was performed by this reviewer.

### Completed failed-proposal v0.5 context, partial grid evidence

Read only its newly completed eight-gene/full-context diagnostic fields.
Whole-model20-node LL−555.387769072968 versus bounded reference
−555.3877690510292 differs approximately−2.194e-8, within the written
1e-6 likelihood tolerance at this point. Forty/eighty-node entry-point
LLs−555.3877690510599/−555.3877690510365; all assembly LL differences
are0 and no production/reference exceptions occurred.

Gene2 twenty-node omega_S score discrepancy is
−1.2733294024180886e-6, versus−3.1900648789218167e-11 at40nodes
and+6.088289872252517e-10 at80nodes. Its20-node LL discrepancy
is−2.082487071675132e-8. Gene5 mode score+1.3419710233969484e-5
retains a false1e-5 stationarity flag at all node counts; gene3's
score+9.32536209234669e-6 is below that diagnostic tolerance but
also warrants keeping its recorded value. No exception is implied
by these flags. References record two/three grid peaks, successful
integration, no warnings or cap truncation. Estimated numerical errors
and tiny Gaussian tails do not establish global mode completeness.

Gene2 FD r discrepancies−7.831471682262503e-8→+3.217918976838874e-7
and kappa−3.2707633629280863e-8→+5.059165976428748e-8 fail clean
refinement and exceed saved quadrature/tail error estimates. Alpha
and variance instead improve approximately as step squared. The
integration reference shares the conditional evaluator and is not
independent verification of the complete likelihood. Finish the last
fixed context before a recipe decision; mass/FD precision remains open.
No model/node adoption, fit, job or active-source change was performed.

### Original P2 draw16: retained finite warm point fails strict stationarity

`bfix_p2_21998518/p2_016.json` records finite warm LL−462.1223306970561
and cold LL−462.12233069708265, convergence flags false/true, and warm
retained under the immutable1e-6 likelihood tie rule. Its interior
alpha(G8), index7 score−1.237540079745969e-6; kappa, index16
−1.0248439447978035e-6; omega_S, index24+1.0452775453320684e-6
exceed the unchanged strict1e-6 projected-score tolerance. The r(G5)
lower-bound score−0.0001061850056535858 at−15 and rho_b upper-bound
score+0.37151443020658 at+3 point outward and project to zero. They
do not establish the convergence defect. The saved fitter already
includes constrained polish; no termination message/history is stored.
[Unverified] Arithmetic, stopping or conditioning causes remain unknown.

This is a finite-fit/stationarity failure, with saved `kappa=NaN` despite
a finite retained parameter estimate. The original supplied lambda pair
[0.5338219269020559,0.732228233248633] remains recorded. It differs
from draw14's retained higher-likelihood cold start; cold's convergence
here cannot change the prescribed likelihood selection. No exception
is established by this record and no passing reclassification follows.
Saved draw SHA256:
`08d512357cf291c4e5d99375d49174dff1a516fb849d97449ca5f795efb18f72`.
Source/JSON inspection only; no fits, tests, imports, jobs or active edits.

### Completed grid21999349: final v2 context and scientific limits

Parent reports terminal COMPLETED exit0,1:46:21, MaxRSS59160KiB,
empty stderr and2,847,320 saved bytes. The completed summary records
6365.037657219917s,176146 evaluator calls, all48 cases retained,
no fits and no MC. In the executed helper, summary writing follows
source/fixture/own-source/private-copy/node-variant post-hash checks
and the unchanged reconstructed synthetic count SHA check. Successful
completion thus supports those immutability checks, not statistical
calibration.

New failed-proposal v2 full-model LL values are20-node
−560.8119193878651, reference−560.811919364894,40-node
−560.8119193649145 and80-node−560.8119193648961. All satisfy
the written1e-6 likelihood agreement at this fixed point; assembly
LL discrepancies are0. No exceptions, false mode-stationarity flags,
cap truncations or warnings occur in this last context. Gene5's
recorded mode score3.1378306469750927e-6 is within the diagnostic
1e-5 tolerance, which is not the outer fitter's1e-6 criterion.

Gene2's largest twenty-node score discrepancy remains omega_S,
−1.2837371823892774e-6, versus+8.701429854429676e-11 at40nodes
and−3.032528117508093e-10 at80nodes. Gene4 alpha discrepancy
is−3.0617157920076465e-8 at20nodes. These are saved numerical
differences against an estimated reference, not rigorously certified
absolute errors. Gene2 r FD discrepancies−1.1607872496544047e-7
and+3.625665465981087e-7 again fail refinement and exceed reported
integration/tail estimates1.23843e-8/3.71684e-8. Alpha errors improve
approximately as step squared; kappa/variance are within saved estimates.

Across six complete contexts, all20-node whole-model likelihoods meet
the written1e-6 tolerance. No production/reference exceptions occurred;
gene5 has false1e-5 mode-stationarity flags at the two earlier failed-
proposal variances0.05/0.5. Failed-proposal gene2 twenty-node score
differences reach approximately1.2e-6; the higher-node integrations
are closer at these points. Repeated nonrefining r finite differences
remain a separate arithmetic-precision qualification. The reference
shares the conditional BB evaluator, so these are independent scalar
integration/assembly checks, not independent BB/complete-likelihood
verification, optimizer convergence, global-mode guarantees or null/
coverage calibration.

The next necessary numerical action is resolving the targeted mass/FD
precision issue and the owner's choice of a frozen integration recipe
before further fixture/calibration claims. This review neither selects
40nodes nor amends the owner20-node rule. Existing complete synthetic
agreement obligations and prespecified statistical calibration remain.
No new model, protected outcome, fit, submission or active edit was made.

### Independent Decimal r2 precision check prepared, not run

Prepared new `codex_decimal_r_precision.py` and
`run_codex_decimal_r_precision.sbatch` under Codex rec_codex. Fixed
contexts are the archived observed-null/failed-proposal vectors, both
v0.5, gene2,24 rows from the exact173-row draw0/count SHA. The completed
grid protocol/summary/executed-helper hashes and their source/input/private
copy guards are pinned before/after. Original archived code is used
only for exact supplied-draw reconstruction; no production BB or
conditional likelihood function supplies the Decimal reference.

Independent full BB probabilities begin with
`P(k=0)=prod_j(b+j)/(t+j)` and positive recurrence
`P(k+1)/P(k)=(n-k)/(k+1)*(a+k)/(b+n-k-1)`, retaining full constants.
Own Decimal logistic/rho/predictors and reciprocals give
`score_r=(D_a[k]+D_b[n-k]-D_t[n])*[-rho'/(rho*(1-rho))]`.
Fixed homozygote BB masses are cached. The oriented three-class
mixture score is its derivative numerator divided by mass; exact
integer selection subtracts the allowed-count derivative numerator
divided by inclusion mass. Every source numeric input uses
Decimal.from_float; no padding, exclusions or fitted replacements.

Digits50/70, fixed saved mode/scale and GH80/160 nodes/weights from
double precision define finite quadrature grids independent of r.
Centered r steps3e-4,1e-4,3e-5 use those same grids. Weighted analytic
r scores and FD derivatives of the finite-grid mass are compared,
alongside digit/node refinements, full BB normalization errors and
saved float80 values. These can isolate precision effects in the
earlier nonrefining direction, but provide no rigorous quadrature
error bound, full blind implementation, global fit or calibration.

Helper SHA256
`e9c6094bc88a8ccd488c2652468fb96cc89427036fbed47ff76363f187012424`;
launcher SHA256
`928efa42f301a0e412064ba7e2abcc3ef4ceaced2b654f07b933c6a4eb3658c9`.
One CPU/2GB/90h cpu/nslab,name1591, unique actual-job-ID output,
32MiB cap, environment/failure phase and all input-row receipts saved.
Static compilation and shell syntax passed. Numerical imports,
metadata-inspect execution, reference calculations and submission
are **not run** here. Parent reviews before compute; no active source
or numerical/statistical rule was changed.

Parent and program_hypotheses independently reviewed the prepared
Decimal source/math without a material defect. Parent metadata-only
inspection, AST and shell checks passed; one CPU1/2GB/90h nslab
job22000606/name1591 was submitted. No second job or diagnostic
expansion was requested/performed here.

First saved receipt, observed-null50digits/80nodes, records analytic
r2 score−0.0009666586347923922085. FD-minus-analytic errors refine
−1.4467249610689214e-11,−1.6074721725934287e-12,
−1.4467249546817483e-13 at steps3e-4/1e-4/3e-5.
Full BB normalization error is4.0e-48. Decimal LL minus saved float80
is+2.354890007943649e-11 and score difference−3.46978888443225e-13.
[Observed pattern, not a guarantee] This first fixed-grid point supports
a precision explanation for the earlier nonrefining r FD direction;
the analytic score remains close to the saved float score. Other
node counts, digits and the failed-proposal context are pending.
No completion, independent continuous-integration accuracy, global
optimum or calibration claim follows from this first receipt.

### Decimal22000606 complete: fixed-point r score is consistent

Parent reports COMPLETED exit0,9:03, MaxRSS60156KiB, empty stderr.
Independently read the completed eight-case summary, SHA256
`254170a7ce62186b397adf126d14c371a7c9395b4f202e059fc215a1f4843f1b`;
541.719040369615s and all postguards passed. No new test or fit was run.

At50digits/80nodes, observed-null r score is
−0.0009666586347923922085. Its FD-minus-analytic errors refine
−1.4467249610689214e-11→−1.6074721725934287e-12→
−1.4467249546817483e-13. Failed-proposal r score is
+0.03821930378842205012, with errors
+5.640454215825779e-10→+6.267171327470387e-11→
+5.640454192322947e-12. Both refine quadratically over the prescribed
steps3e-4/1e-4/3e-5, with essentially identical patterns at70digits
and160nodes. Saved50/70-digit score differences are of order1e-47;
80/160-node score differences−1.128117160106513e-19 and
−1.279166938985992e-19. Differences from saved float80 scores are
−3.46978888443225e-13 and+2.83895077292898e-13.

[Observed pattern, not a guarantee] Independent Decimal BB/mixture
arithmetic supports the earlier nonrefining mass finite differences
being a precision issue at these two fixed points, rather than an
incorrect analytic r-score formula there. It does not identify every
floating-point operation responsible. This resolves the focused r
precision qualification without changing production equations or
using a new fitted point. Fixed local centers and double GH nodes/
weights remain; independent continuous-integration accuracy, original
adaptive-mode behavior, all-coordinate derivatives, optimizer/global
fit and calibration are not established. The earlier20-node score
and mode-stationarity qualifications remain separate. No additional
diagnostic, job, source adoption or statistical-rule amendment follows.

### Frozen sequence comparison: what completion can establish

Source-only review of `target_readout_development_protocol.json`,
`target_readout_pilot.py`, `describe_target_readout_development.py`,
the current AG program and existing development interpretation. No new
score/effect outcome was opened. Scorer21997286 generates same-output
501bp/full-peak native contrasts; evaluator21997998 is the fixed
descriptive comparison on24,430 common development lead pairs, fitting
folds2–4 and evaluating already inspected fold1. It includes zero,
train-mean, center OLS, peak OLS and center+peak two-feature OLS;
no tuning, added genomic covariates, significance refiltering, margin
or pass/adoption criterion. Bootstrap intervals condition on fitted
coefficients/deposited effects and use chromosomes, not new donors.

An assay-native measured endpoint exists: unchanged Currin FastQTL
signed normalized-accessibility beta per additional ALT dosage. This
is a cohort marginal association on significance-selected leads,
not an isolated editing effect, arbitrary-target discrimination or a
gene/MASLD causal endpoint. Completion can establish whether exact
peak alignment changes this native association-prediction readout,
and whether the fixed simple calibration/combination helps in this
development source. It cannot establish a novel foundation model,
a target-conditioned learned sequence model or independent biological
validation. The two-feature OLS is the strongest simple combination
inside this frozen experiment; the protocol intentionally excludes
specialists. The older5% complexity margin belongs to the prior
adapter/combination development task and is not automatically a
new paired-mask criterion.

Existing specialist code routes are real, but their eligibility on
these exact pairs is not established here. `c1e2_00_build_variant_lists.py`
names `hepatocyte_5fold_v2` and `gse281367_hepatocyte_5fold_v2`
ChromBPNet families. `c1e2_03_recompute.py` names their
`currin_calibration/currin_scored_primary.tsv.gz` files and reads
heldout_fold alongside beta/p/q fields; these are mixed artifacts,
not safe development-only input to this review. Borzoi's existing
`c1e2_02_borzoi_score.py` uses a524,288bp window, two ATAC/five
DNase liver tracks and16-bin summation. It is a distinct native
readout, not interchangeable full-peak AlphaGenome features. Source
existence does not prove complete matched coverage, eligible checkpoint
exposure or the strongest specialist performance.

One concrete uncompleted preparation before a stronger comparator
claim is an owning-builder **metadata-only development coverage and
eligibility export** for those existing specialist routes on the exact
24,430 frozen keys: folds1–4 only, checkpoint/held-chromosome identity,
assembly/REF→ALT orientation, native readout definition and exposure,
plus admitted/missing pair keys. No predictions, beta/p/q values or
fold0 should enter that export; do not parse mixed canonical files here.
This determines whether an already existing specialist can be a fair
same-population comparator. It does not add an arm to the pending
frozen experiment, authorize scoring/fitting, select a new candidate,
or supply the still-missing independently admitted receiving endpoint.
No new source search, model rescue, tuning, job or publication was done.

### Existing specialist metadata artifact completed

Created exclusively
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/codex_specialist_metadata_20261001.json`,
SHA256 `daf8a324416132d079e5376078102c9d89339cf7b42686beb3c2b3862df9763e`
after the source-wording corrections below; JSON parsing passed.
Only the three requested c1e2 source scripts, their explicitly referenced
4,375-byte selected-track metadata CSV, fold-manifest headers and root
symlink metadata were inspected. Source hashes, exact paths, declared
versus unverified properties and required owner-export fields are saved.
No fold rows, mixed scored files, calibrator outcome metrics, prediction
arrays, effects, p/q fields or protected data were read.

ChromBPNet family directories and consumer `logfc` fields are verified
source declarations; checkpoint hashes, exact score formula/units,
exact native target/region mask/peak conditioning, biological training
roster and exposure are not
established by those consumers. Borzoi declares configurable checkpoint
default `johahi/borzoi-replicate-0`,524,288bp inputs,float32 and raw
ALT-minus-REF summation across33 central bins (default radius16),
then separate NaN-skipping ATAC/DNase means via `np.nanmean`.
[Unverified] Actual finite-track availability/handling was not inspected;
`np.nanmean` skips NaNs but does not remove infinities.
Consumer orientation multiplies
by `orient_mult` to reach the label REF→ALT frame. No log/beta transform
is applied by that scorer. Actual checkpoint execution and training
exposure remain unverified.

Selected Borzoi track metadata is not uniformly adult primary liver:
ATAC labels include fibro liver/adrenal and hepatocyte; DNase labels
include three HepG2 tracks, hepatocyte and H9-origin hepatocyte.
Track indices alone do not prove current checkpoint-axis compatibility.
The consumer builds variant-only keys, while the pending comparison
uses variant–peak identities. The owner export must label broadcasting
one variant-only score to repeated peak pairs and preserve dependence;
matched coverage alone does not establish target discrimination.
Exact24,430-key coverage remains unverified. No scoring, fitting,
installation, download, job, source edit or frozen-arm addition occurred.

### ChromBPNet producer-source interval and effect addendum

Source-only addendum created exclusively at
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/codex_chrombpnet_source_addendum_20261001.json`,
SHA256 `e44b0e579fa1a0aa3b35dd2fc42025542e2d6df836a79bcccbdc50d7d7d433be`.
JSON parsing passed. It records ten source hashes; the prior metadata
artifact is preserved. Root explicitly included the discovered same-series
ASE launcher and its named producer/imported helpers in this bounded review.

`c1e2_01_chrombpnet_ase.sbatch:49–57` invokes
`GWAS/finemapping/src/chrombpnet_variant_scorer/src/variant_scoring.py`
with each held-chromosome `chrombpnet_nobias.h5`, REF/ALT variant lists,
the pooled summit peak file, seed42 and no forward-only flag. In
`utils/helpers.py:169–178,195–197,225`, each strand's count head is
exponentiated, forward/reverse-complement counts are averaged, and
`logfc = log2(C_ALT/C_REF)` with no pseudocount. This is a dimensionless
predicted count ratio, not the cohort's marginal ALT-dosage beta; averaging
occurs before the log ratio.

`generators/variant_generator.py:38–47` converts one-based position to
p0, extracts `[p0−L/2,p0+L/2)` for equal-length alleles and replaces the
base at index L/2. L comes from model input shape
(`variant_scoring.py:38–40`); numeric input/output widths remain
[Unverified] because no checkpoint body was inspected. The scalar count
head supplies logfc. No profile-bin target mask is applied to that ratio.
The equal-length substitution branch lacks a FASTA REF assertion; the
indel branch has one at line55. Actual per-key reference checks therefore
remain an owner-export requirement.

`utils/helpers.py:205–215` delegates logfc/JSD unchanged to the scorer
without peaks, then separately calculates allele-count percentiles against
peak predictions. The supplied peak set calibrates percentiles; it does
not condition logfc on the variant's paired target peak. A repeated
variant under one checkpoint consequently has the same source-defined
logfc across peak pairs. Peak overlap metadata does not establish a
peak-conditioned effect. This resolves the ASE launcher's formula and
conditioning, not exact24,430-key execution coverage or person exposure.

[Unverified] The ASE launcher's comments15–17 assert the same recipe as
deposited Currin calibration, but the permitted routing sources do not
name its original caQTL producer/execution receipt. Source equality for
that historical output, actual checkpoint hashes and numeric genomic
output span remain unverified. No scored files, fold/reference rows,
effects, molecular arrays or model bodies were read; no scoring, fitting,
submission, model-arm addition or tuning was performed.

### Bounded original caQTL producer search

The permitted search did not locate an original Currin caQTL ChromBPNet
scoring launcher or a receipt naming its executed producer. This is a
finite unsuccessful search of source/documentation in
`scripts/analysis/alphagenome_program` and
`GWAS/finemapping/src/chrombpnet_variant_scorer`, not a claim that no
producer exists elsewhere. Searches used `currin_scored_primary`,
`currin_external_calibration`, `currin_calibration`, `variant_scoring.py`
and family names, with code/document extensions only. No data or result
file was opened.

Additional matches are consumers, not scoring receipts:
`c1_matched_endogenous.py:30–34` names the mixed scored inputs and
`:559–565` hashes those inputs and the calibration JSON into its own
comparison manifest. `c1e2_05_recompute_v2.py:44–54` repeats the same
input declarations. `c1e2_06_qc_addendum.py:28–30` points to calibration
metrics; it does not name the original scorer execution. The upstream
scorer `README.md:54–73` supplies generic scoring instructions without
a Currin/hepatocyte run identity. Filename/source-text inspection of
the scorer tree found no Currin/family-specific producer link.

SHA256 of the newly inspected source files:

- `c1_matched_endogenous.py`:
  `0b6d9a4503026fc8924f5841c5f89a13544d38dfbdd67be13713ef3c3cb8a42d`.
- `c1e2_05_recompute_v2.py`:
  `5fc2be02e66146742e27cffb88e5ca5b64d53c3465772d2eeb2994184b953df4`.
- `c1e2_06_qc_addendum.py`:
  `ff065359f995c61966880c06273c9d90a9775e9e6f8397654079782220bf827c`.
- Scorer `README.md`:
  `694f610ac551dd36c75f2c24642175451db323e4a9c3d7f48b06c4353f1610b6`.

[Unverified] The resolved ASE source formula cannot yet be assigned to
the historical caQTL score by executed-run evidence. Original producer
source/launcher plus its executed source/checkpoint identity remains
the concrete missing owner evidence. Both existing metadata JSONs were
preserved; no outcome/score/checkpoint-body read, numerical workload,
job, fit or comparison-arm change occurred.

### Broader authorized source search resolves declared Currin producer

The subsequently authorized code/doc/launcher search under
`GWAS/finemapping/src` located original caQTL recipe declarations in
`seqfunc_v2/chrombpnet`. The earlier unsuccessful two-tree search is
preserved above with its narrower scope. Neither metadata JSON was
overwritten, and no outputs, score rows, fold rows or checkpoints were
opened.

`11_score_currin.sbatch:18–41` explicitly selects the adult
`hepatocyte_5fold_v2` family and invokes the inspected
`chrombpnet_variant_scorer/src/variant_scoring.py` with each chromosome
fold's `chrombpnet_nobias.h5`, the stated hg38 FASTA and summit peak file,
max30,000 calibration peaks, batch1024, seed42, zero shuffled scores and
no HDF5. `24_gse281367_score_currin.sbatch:18–32` declares the same
recipe for `gse281367_hepatocyte_5fold_v2`. Neither passes forward-only.
These explicit declarations establish the source-defined caQTL effect:
log2 ALT/REF of forward/reverse-complement averaged exponentiated count
heads, using variant-centered sequence inputs. The supplied peak set
calibrates percentile columns without altering logfc or selecting a
paired target-peak mask, as derived from helpers above. Numeric model
input/profile-output widths still require checkpoint metadata.

`10_prepare_currin.py:55–61` checks each one-based hg38 REF against the
stated FASTA and retains reference matches before generating score
inputs. `:63–71` requires chr1–chr22 each assigned exactly once to a
test-chromosome fold; `:76–77` writes chromosome, position, REF, ALT and
variant ID in that order. Thus declared caQTL preprocessing supplies
the SNV REF check absent inside the scorer's equal-length branch.
Its occurrence in the actual historical execution remains unverified.

`12_calibrate_currin.py:103–117` reads each fold's scorer output, refuses
duplicate variant IDs, joins metadata one-to-one, requires complete
logfc and matching held-chromosome folds, then adds binary training-peak
overlap and writes `currin_scored_primary.tsv.gz`. It does not reorient,
recalibrate or rescale logfc during this join. The orientation-repair
source `13b_promote_currin_canonical.sbatch:42–55` rebuilds signed
ALT-dosage truth and reuses the existing score directory; this is a
source declaration, not an inspected executed repair receipt.

SHA256 source hashes:

- `10_prepare_currin.py`:
  `b3ef3aca5dc8ce5b7c69fb7c4bbabf684d0609a494d2eeca8a642c86912040ff`.
- `11_score_currin.sbatch`:
  `2b83c3cf6ffeb2e1a9a17dd36b728d398cc5fcd746305b5f6e7e925e66a16a2b`.
- `24_gse281367_score_currin.sbatch`:
  `ae1dec98a9bf99baee72987ae8b146c4b3db03aab8654016478662723707dc63`.
- `12_calibrate_currin.py`:
  `922161b0ef91f8e329a051960559215fa56771c7207d330bdbfa3034dcd9b811`.
- `13b_promote_currin_canonical.sbatch`:
  `aaf694ceaaa0c2693bb57856d3c1908965913c46978721e5fd713f3e8e77444d`.

[Unverified] Actual historical executed-source/checkpoint identity,
complete current24,430-key coverage, source/person exposure and numeric
checkpoint dimensions remain unresolved. Locating code resolves the
declared recipe, not executed equality or independent validation. No
job, fit, acquisition, numerical test, new arm or outcome read occurred.

### Fixed-parameter conditional preparation source ready, 2026-10-05

Prepared new isolated helpers only:
`scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_prepared_conditional.py`
(SHA256 `686eccc5f4c4c85f6dfde8274564a04bedbe6019c25e6c3709948f4a2b462809`)
and `repair_p3_prepared_conditional_check.py`
(SHA256 `64aff1020ae4696b5cb4c9f93d6c32e978cf0338258e09b131e9dfb97dbcb14c`).
Both AST parses passed; numerical checks are **not run**, no job submitted.

The context retains one fixed parameter vector and unchanged synthetic rows,
caches dispersion-only BB denominator/chain terms, support log-binomial
constants, oriented genotype priors, fixed homozygote masses and inclusion.
It recomputes eta-dependent rising prefixes, the original tail/complement
selection branch and all conditional score/Hessian mixture terms. Each call
rejects changed parameters, rows or dispersion. It does not patch original
model methods or conditional globals. Algebra/probability law is unchanged;
actual floating equivalence and speed remain unverified.

The check reuses the immutable executed reconstruction helper from job22039132
(SHA256 `d3f529ad20666c58dc8a84d01785d3f62c42f18a430a3a572a055c180b9b4a9e`),
kernel461d and stable probe7237. Only occurrence-guarded output-label/prefix
changes and the preparation-hook routing differ. Existing synthetic fixture,
parameter/count/environment/source and pre/post hashes are retained. New
source files are archived by the reused helper. Fixed cases are failed draw0
gene2 and positive-variance draw18 gene2, each24 rows at u=−.5,0,.5 and both
score-only/full-Hessian paths. Full conditional log mass, both scores and
all four Hessian entries are compared with the original stable evaluator;
the focused absolute comparison allowance is1e−10, not a replacement for
integration/fit accuracy criteria. Timings include preparation separately;
changed-parameter refusal and valid reuse afterward are checked. No integral,
optimizer, MC/calibration or biological data is evaluated by this check route.

**Executed v1 update, root-reported:** jobs22039592/22039593 completed12s.
All288 conditional row/u/H-mode comparisons were bit-exact for log mass,
both scores and2×2 Hessians. Failed-case score/H timing changed
.15246→.09608s and .15970→.08780s, with .02360s preparation;
positive-case .11687→.07587s and .12987→.08998s, with .01801s
preparation. These are fixed-point conditional timings, not full-integral
or optimizer speed, accuracy or calibration evidence.

**v2 source ready, not run:** separate
`repair_p3_prepared_conditional_v2.py`, SHA256
`289e93b651411042093a2b987717d87cdb5ed156ab7c736e1ac5cf820670fd4b`,
and `repair_p3_prepared_conditional_v2_check.py`, SHA256
`7443b61cea8565b6ea27320528d6a5ddcd088aff70cc0e8354e5055368013af1`.
Both AST parses passed. Eta-dependent a/b rising tables are computed once
per row-call and indexed by each original support branch without changing
the downstream operation/reduction order. Local tables are not retained
between calls; persistent row/support cache layout is checked unchanged.
The same fixed two-case checks, mutation refusal, preparation/warm timings
and immutable reconstruction guards apply. Original v1, live integrations,
fitting rules and accuracy budgets were unchanged; no numerical execution,
job, fit, integral or MC was performed by this preparation.

**Executed v2 update, root-reported:** jobs22039642/22039643 completed;
all288 conditional comparisons remained bit-exact. Warm score timing
ratios were approximately2.06–2.08 and Hessian ratios1.81–1.84, with
.0166/.0200s preparation. This remains a conditional fixed-point benchmark,
not an integral or fitting performance result.

**Prepared GL source ready, not run:**
`repair_p3_prepared_gl_check.py`, SHA256
`686fe1e6aff48366b0dc696cd4d14aec0a0fcdc38ec45a40f9284a16463b17d7`.
AST parsing, exact archived callback-source occurrence and both saved
nine-ladder metadata checks passed. The immutable snapshot checker3a7e
and composite GL45b442/kernel461d use an explicit conditional callback,
so one prepared v2 context can serve all ladders at one fixed p without
adding a global/model-method patch. All original grids, orders, domain
expansion, caps and precision refusals are retained.

Completed baseline candidates22039200/22039201 and all14 complex-reference
protocol/summary/record files are hash-guarded. The new output compares
full LL, every score and Hessian entry with the saved uncached GL result
under unchanged1e−9/1e−8/1e−5 budgets, requires the same nine ladder
identities/evaluation counts, and reports each of the42 existing complex
records per case. Complex reference uncertainties remain uncertified;
comparisons do not invent derivative acceptance or an independent full-fit
claim. Refusals and partial original outputs remain retained. No integration,
numerical-library import, fit, MC, calibration or job was executed in this
source preparation. Root reviews/submits separately.

### Batched GL integration source preparation

Prepared separate `repair_p3_batched_gl.py` (SHA256
`8bffaf8057ace5a2593861f879068f7e4a3e3f1ef72d2fecc4199b94c8372c52`)
and `repair_p3_batched_gl_check.py` (SHA256
`e88bddc3ffc6056564d068dc2e93550249ad8bc7b48afeda2a8fcf438b113912`).
AST parsing and exact archived callback-source occurrence checks passed;
no numerical-library import, integration or job executed here.

The explicit vector conditional callback is owned separately by the other
agent, source184b3e. Root reports its jobs22040051/22040052 completed
7,104 conditional records each with zero outside the fixed comparison
allowance: LL exact, maximum score difference7.1e−15 and Hessian6.82e−13;
mutation/callable/layout guards passed. Its conditional speed measurements
do not establish integral or fit performance. This checker requires the
passing exact-case receipt hashes before integration; no callback accuracy
is inferred from source alone.

Batched GL requests32 or64 nodes, lets the callback internally split each
row to at most100,000 count-axis entries, and consumes every node in the
archived node/row-matmul/Kahan order. All nine order/domain/grid ladders,
normalized Gaussian moments, error budgets and cap refusals are retained.
Fixed Jacobian/Hessian reuse requires explicit predictor invariance checks
at u±.5 and actual chunk endpoints; eta uses the same base+u×stage order.
GL workspace and callback static/batch estimates are recorded separately,
with their combined declared estimate checked before GL allocation. These
estimates are neither measured process RSS nor a certified allocation bound.
Adaptive requested/effective chunk ledgers are retained in the output.

Required comparisons use the hash-guarded uncached candidates22039200/01,
all LL/score/Hessian entries, identical nine-ladder identities/evaluation
counts, and all42 existing complex-reference records per case. Optional
completed prepared-GL comparison requires its three exact input hashes;
unfinished receipts are not read or dynamically substituted. All precision
refusals remain explicit. No new inference method, fit/calibration result,
statistical rule change or frozen20-node adoption follows this source task.

### Gene-local moments isolated source proof ready

New `repair_p3_local_moments.py`, SHA256
`95ec9f270d18bfbbc8f31f0857308f5f2acad053661a30b208049b7a3f59f31f`,
and `repair_p3_local_moments_check.py`, SHA256
`e322cf91b7290741b9281552446cfced3179fe2dfad60482df1f8cb04c4a040a`,
passed AST parsing. Numerical comparisons are **not run** at this source
handoff; no job submitted here.

For one nonempty gene the original predictor touches alpha_g, r_g, kappa,
Z delta coordinates, omega and three dispersion slopes: Z+7 parameters,
Z+8 with variance. The helper writes the archived predictor formulas
directly on those active coordinates. It retains row ordering, conditional
probability calls, score/Hessian transforms, normalized Gaussian prior
moments, cp/hp arithmetic and extrema, then embeds vectors/matrices into
the full coordinate space only at the comparison interface. No row/depth,
supplied active parameter, covariance term or budget is discarded.

The guarded existing reconstruction supplies both original positive-variance
fixed cases, kernel461d/stable7237 and validated batch184b. Checks compare
all29 augmented coordinates at G8/P28 and all771 at synthetic G379/P770,
using the identical24 gene2 rows. Original alpha/r positions and shared
coordinates are explicitly mapped; unobserved extra gene slots add
dimension only. Each case retains12 point/mode records at u−.5/0/.5,
all nonzero Hessian/second-moment differences, local/dense cp/hp embeddings
and maxima, tail-function comparisons at24/26 and mutation refusal.
Fixed capacity estimate is below1GiB, not measured RSS. Scalar logmass
for tail comparisons is the hash-pinned existing numerical estimate,
not a certified lower bound.

This isolates coordinate-locality correctness. It does not validate379
observed genes,558 participants, broader real depths, integrals, optimizer
behavior, calibration or sealed replication. No frozen helper was edited;
lazy preparation and row/depth scaling are deliberately outside this proof.

### Gene-local moments executed, then isolated local GL source prepared

The preceding source-only `not run` statement is historical. Root submitted
22040916/17; both completed, and the saved moment comparisons each retain12
records. Across both cases all24 records pass: LL and Hessian differences
are exactly zero, maximum gradient difference1.11e−16 and second-moment
difference2.70e−13. cp/hp embeddings are exact at P28 and P770. Receipts are
`REC/repair_p3_local_moments_failed0gene2_22040916` and
`REC/repair_p3_local_moments_positive18_22040917`. These remain pointwise
checks on24 unchanged synthetic rows, not fullscale biological validation.

New source-only files are `scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_local_gl.py`
(SHA256 `254a1efce7d8a873f8caeda1f0630cb48e35286e15bd952352246ab26b01e4de`)
and `repair_p3_local_gl_check.py`
(SHA256 `1fac996bf63823b95fd5e0e53b45c8587778502febc4acf0f76033b9f65f9ef6`).
Both passed AST parsing. No numerical import, integral, fit or submission was
run here. All existing files remain intact.

The GL implementation derives from immutable batched8b and the passing
local-moment95 helper. It retains all nodes, orders, domains, normalized
Gaussian score/Hessian moments, coordinatewise compensated accumulation,
error allowances, caps and refusal rules. Only active Z+7 parameter
coordinates and variance enter allocated node/moment/tail arrays; full
vectors/Hessians are embedded at return. Report ladder and allowance arrays
are explicitly **local**, with global indices supplied. A future full-model
caller must embed these arrays before summing across genes.

The checker pins the completed batched32 sources/outputs22040257/59 and
passing moment receipts22040916/17. It refuses differing parameter/count,
row, variance, draw or environment identities before integration. Each case
runs the identical nine-ladder integral at original P28 and dimension-only
P770, and compares every ladder's LL, score, Hessian, tail allowances,
order-refinement estimates, grids/evaluation counts and final full interface.
Inactive baseline entries must be exactly zero rather than silently masked.
All nonzero interface/lattice differences and all failures are retained;
no fit or additional draw is made.

[Estimate] At local augmented width15, batch32 and unchanged support caps,
the declared GL workspace is below126MB; callback storage and temporary
working estimates each below53MB. P770 output is about4.76MB, with eight
such output-sized buffers budgeted for comparison. Combined declared
allocation is below1GiB, not measured RSS or a certified bound. Existing
batched fixed-point kernels took108.6–178.6 seconds; two integrals per case
suggest about4–8 minutes with buffer, not measured local runtime. The
checker records actual input-specific estimates before integrating. This
addresses coordinate allocation only; MAX_ROWS256 and the100000 support
entry cap remain unchanged and do not establish full558-person admissibility.

### Local GL executed comparison and unchanged-decision JSON supplement

The checker callable guard was corrected from nonexistent `bb_support` to
archived `bb_terms` before submission; final checker SHA256 is
`ac28c6656b759d1077693a278f13f2ef99270c6029c7d0d1b7b397ee5421bee5`.
Earlier checker hashes above are superseded. Root submitted immutable snapshot
`REC/repair_p3_local_gl_launch_20261005T153915Z` as22041594/95.
Both completed exit0, empty stderr: failed0gene2 elapsed5:58, MaxRSS63988KiB;
positive18 elapsed4:20, MaxRSS63052KiB. Four integrals returned, each with
nine ladders, at P28 and dimension-only P770 on24 unchanged gene2 rows.
Final maximum saved differences versus batched32 are LL4.2633e−14,
gradient1.88294e−13 and Hessian4.18368e−11. Archived callable guards,
inactive-coordinate checks and final source/input checks completed. P770
integral times177.799/128.650 seconds closely match P28 times176.887/128.387;
this is not a runtime estimate for379 observed genes or558 participants.

**The original absolute-comparison flags remain false.** Initial-domain tail
allowances are intentionally unresolved and enormous, then trigger domain
expansion. Failed0's initial allowance differences reach1.01705e14 for
scores and3.63422e41 for Hessians against magnitudes1.78779e27 and3.19618e54.
Positive18 retains absolute score allowance difference2048 and Hessian
allowance difference1.51116e23; corresponding old/new values are
1.4347188036737665e19/1.4347188036737667e19 and
7.781936028684471e38/7.781936028684472e38. These fail the checker's absolute
cross-implementation allowance comparisons. No original flag is changed.

Root requested a descriptive, saved-JSON-only supplement. New stdlib-only
`repair_p3_local_gl_collect.py` reads hash-pinned completed JSON, verifies
inputs before/after, and creates (exclusive write)
`REC/repair_p3_local_gl_json_supplement_22041594_22041595.json`, SHA256
`dd0082282cc4c235715186be771fbac349818e0be2bdce08b72e233cf81a0fe8`.
AST parsing and this lightweight JSON collection were executed; no evaluator,
numerical library, fitting or quadrature was executed on login. For all36
saved ladders, old/new tail-budget and order-refinement decisions agree,
as do domains, orders and evaluation counts. Both cases reject initial8,
then accept16/15 respectively and finish expanded18/17. Saved final full
allowances pass the unchanged1e−9/1e−8/1e−5 budgets in both versions:

- Failed0 old/new LL allowances4.76064e−12/4.73221e−12,
  score6.50824e−11/6.49669e−11, Hessian1.03162e−8/1.03128e−8.
- Positive18 old/new LL5.18696e−12/5.21538e−12,
  score2.41732e−10/2.42119e−10, Hessian5.16993e−8/5.17939e−8.

Scale-normalized allowance differences, using `max(1,abs(old),abs(new))`,
are at most1.14026e−13 for failed0 and3.67068e−16 for positive18. These are
**descriptions**, not new acceptance criteria. The supplement shows why
original absolute comparisons fail without changing integration decisions;
it does not certify quadrature/floating errors, independently validate the
probability law, adopt a model, assess a fit or establish calibration.
Declared allocations were141250944 bytes at P28 and179288832 at P770;
observed process memory above is separate. Row/depth caps remain unchanged.

### Isolated eight-gene process-parallel point evaluation prepared

Root's ongoing observed fit supplies a converged null receipt22040873
(771.47 seconds; projected score9.0844e−8). Its positive-variance marginal
evaluation call3 took617.78 seconds, with all8 genes accepted under existing
per-gene budgets. This supports investigating independent gene evaluation;
it does not establish parallel runtime or a converged marginal fit. Root's
separate wide-GH jobs are not duplicated.

New source-only `repair_p3_parallel_gene_model.py`, SHA256
`41228f20f1ee4e089e84cff5536f647c16785e40a9e68765269486105a77cd45`,
and `repair_p3_parallel_gene_check.py`, SHA256
`4a6c10daa50a0ef0d952813d22c20e3f8fd2687d9f5a41f0d96173f77647c1d3`,
passed AST parsing. Numerical work is **not run**; no job submitted here.
The implementation uses frozen batchedGL8b/conditional184b and sequential
model2f64, not localGL. Eight Linux-fork workers inherit the isolated
private modules and evaluate one gene each from read-only copied parameters.
Parent aggregation retains original ascending-gene `+=` order, all
normalized-prior moments, original per-gene and aggregate budgets, boundary
law and inherited optimizer/start/bounds methods.

All worker outcomes are collected before forming a successful full-model
sum. Numerical failure preserves the first failed gene in ascending order,
its partial ladders and completed prefix, plus all other available worker
records. Known `LinAlgError` is re-raised as `LinAlgError`; ValueError retains
its numerical-failure disposition. Unexpected worker exceptions and pool
infrastructure errors remain explicit; unavailable results are not zeros.
No failed fit is reclassified. Rows, parameter immutability and original
callables are checked, and no live or frozen module is modified.

The helper declares exactly one positive-variance full-Hessian point:
converged22040873 observed-null parameters and the existing positive18
variance. It uses original60-participant/173-row/G8/Z7/P28 synthetic inputs,
source/fixture admission before imports, and no fit or new draw. Sequential
and parallel outputs retain all gene reports, node/lattice identities,
all gradient/Hessian/error entries, gene terms, timings, source postchecks
and cross-gene alpha/r zero checks. Complete numerical reports must be
exact after removing callback setup timing. Refusal injection is not done;
if an actual refusal occurs, it is retained without claiming success.

Declared aggregate numerical workspaces are calculated before forking;
interpreter/library/parent capacity adds9×256MiB, with combined estimate
capped6GiB. This is not measured RSS or a certified allocation bound.
Execution requires8 allocated/available CPUs and BLAS/OMP/MKL/NUMEXPR
thread counts all1 before numerical imports;8GiB is a proposed buffered
allocation, not an observed requirement. The measured617.78-second
score evaluation is the only full-model runtime basis; this new full-Hessian
point's runtime remains unverified. No actual558-person scale, calibration,
integration revision adoption or useful biological model is established.

### Eight-gene parallel point executed: exact saved agreement, bounded speed result

Root submitted22041812 from
`REC/repair_p3_parallel_gene_launch_20261005T160057Z`. It completed14:58,
exit0, empty stderr,8 CPUs; SLURM sampled batch MaxRSS184520KiB. This is
sampled scheduler memory, not a certified simultaneous-worker peak.
Saved `REC/repair_p3_parallel_gene_check_22041812/comparison.json` SHA256
`1e98cbc1f78745b5c43cadb1d5c96afae3fd05bed9e11c3c281c63d689f1f9f0`
and `summary.json` SHA256
`50c86a2d6f77e63843a77588ff76cc35b5303f1635c93486f00b5a4a1764b9b7`
report PASS. The preceding `not run` text describes the source-only handoff.

At this one predeclared synthetic positive-variance full-Hessian point,
sequential and parallel total LL, gradient and Hessian differences are
**exactly zero**. All8 full per-gene numerical reports, including every
nine-ladder sequence, node counts, tail/refinement allowances, callback
storage/chunk ledgers and method fields, are exact after removing setup
time. Gene terms, aggregate error arrays, inactive other-gene alpha/r
scores/Hessians, worker fixed inputs and original reference callable checks
all pass. Final source/input checks completed. No numerical refusal
occurred and none was injected, so this execution does not establish
behavior under an actual worker/integration failure.

Sequential evaluation took727.9817257 seconds, versus160.8313972 seconds
with8 fork workers:4.5263657× at this point [Observed pattern, not a guarantee].
The slowest worker was gene2 at160.7305702 seconds; workers ranged14.3771–
160.7306 seconds. Declared aggregate workspaces875494272 bytes plus
2415919104 bytes interpreter/library/parent buffer total3291413376 bytes,
not measured or certified process memory. This supports bounded practical
parallel point evaluation of frozen batchedGL8b on the small synthetic
fixture. It does not demonstrate complete optimizer speed/convergence,
positive-variance accuracy beyond the fixed checks, statistical calibration,
real558-person scale, independent biological validation or localGL adoption.
No original fit was cancelled or modified by this review.

### Final-fit requirements: bounded binding-source review

Binding references are the likelihood spec §§4,5.1–5.3,7,11
(`2026-09-29_model_A_likelihood_spec_v1.md:219,244,253,260,301,343,401`)
and v1.7 answers §11
(`2026-09-30_claude_to_codex_design_fixture_answers.md:385–414`).
The immutable Codex source is
`REC/p3_parallel_20261003T130323Z/executed_sources/likelihood_v1.py:258–315,452–525`
and its `tests_v1.py:18–22,51–130`. These are synthetic source references;
no protected or real allelic inputs were opened in this review.

**Already prepared/executed.** Both observed-fit runners save a recomputed
kept-point LL, all parameter gradients/Hessian, coordinate bounds,
projection mask, strict maximum projected score, native convergence,
right variance score at zero, numerical allowances, Hessian finiteness
and start LL/convergence/kept label. See
`repair_p3_batched_observed_fit.py:173–230,281–330` and
`repair_p3_parallel_observed_fit.py:267–319,457–499`.
The parallel fixed-point22041812 is exact versus sequential full GL at
one declared point. Fixed positive-variance derivative/reference diagnostics
and broader numerical checks already reported above are not calibration.
No new duplicate check is required merely because the final runner finishes.

**Smallest sufficient kept-point check after completion.** Read the finished
runner diagnostics and final source/input checks first. At each kept point,
require finite LL/score, recomputed LL agreement within the existing1e−6
fixed-vector target, and recomputed gradients within1e−8+1e−6×abs(reference).
Convergence requires max absolute projected score strictly<1e−6 on every
free coordinate, not merely optimizer status. Fixed/no-information coordinates
are removed; at a lower bound a negative score is outward and removed, at
an upper bound a positive score is outward and removed. Original nuisance
projection tolerances are1e−8 for section3 and1e−10 for section4. At v=0,
add the right variance score≤1e−6; variance Hessian entries may be unavailable,
while nuisance Hessian entries must be finite. A native convergence flag
must not override failed recomputed stationarity. The inherited null snap
can expose the distinct projection tolerances, which the runner explicitly
labels rather than silently correcting.

If both fits yield different retained points, compare their retained LLs
and warm/cold LLs before deciding whether a fresh same-point numerical
comparison resolves a concrete discrepancy. v1.7 keeps cold only when its
LL improves warm by≥1e−6; otherwise warm is retained. A higher LL beyond
the agreement target is an optimization discrepancy only once both
objectives are numerically credible. Agreement of a numerical GL revision
with another copy of GL does not establish owner20-node accuracy or amend
its prescribed numerical rule.

**Boundary branch evidence remains distinct.** Spec§4 says a zero-variance
optimum has right score≤1e−6 and the branch started atτmax² returns to0.
v1.7§11.4 requires a positive-boundary-score zero fit to restart fromτmax².
The immutable fitter implements the restart at score>1e−6, but existing
returns retain warm/cold LL/convergence, not both terminal parameter/variance
vectors or every attempted restart. A kept v0 receipt alone cannot prove
the other start returned to0, nor prove an unsuccessful restart was not
attempted. If those branch outcomes were not captured, label them
[Unverified]; this cannot be repaired by inventing history from the kept
point. Do not add a new optimizer/start policy. Future authorized execution
could record existing branch calls without changing them.

**P2 profiles after credible observed fits.** For each label source use the
existing paired reliability draws: invalid reliability enters0 in the
linear5th percentile over all B draws; κ/λ90% linear quantiles use valid
draws. Atτ0=τmax×lambda_fifth, fit v fixed toτ0² from exactly the two
specified starts (section4 MLE and section3 MLE with that fixed variance),
retain higher LL with the same warm tie rule, and recompute the kept
profile's LL/projected nuisance scores. The variance coordinate is **fixed**
for this check: an unconstrained variance score need not vanish. The current
observed-fit diagnostic's `(0,None)` variance bound must not be reused for
such a profile. Use max(LL4,LL3) as the top profile LL, clip a negative LR
loss to0, then sign(τ0−τhat)×sqrt(LR), p_tau=Normal.sf(root).
Per-source p_tau and p_up/p_lo remain separate; final P2 is their maximum
over both sources. Failed κ refits count in both tails; failed profiles
cannot produce a passing P2. Real GateB boundary rejection>0.06 triggers
the prescribed1000-draw signed-root bootstrap, not a new asymptotic choice.
No such coverage/size result follows from these synthetic point checks.

**Shared draws and failures.** v1.7 requires all50 supplied B-FIX draws for
P1/P2/P3, final per-draw LL/retained start and test-specific components;
these are agreement checks, not the production1000/2000/1000 laws or full
GateB calibration. P1 uses|z*|≥|z|−1e−6, reports all comparison-band draws
within1e−3; observed null/sandwich failure means not computed. If a new
observed null is used for a new P1 run, its relative-cluster sandwich must
be recomputed once at that actual point, not borrowed from a different fit.
P2 paired κ/λ uses its prescribed inequalities and invalid/failure handling.
P3 uses that draw's section3 MLE for both section4 starts; either null or
kept marginal failure makes a failed draw/exceedance. LR is exactly0 for
kept v0, otherwise2(max(LL4,LL3)−LL3); compareLR*≥LR−1e−6. Counts retain all
supplied draws, p=(1+exceedances)/(B+1), MCSE=sqrt(p(1−p)/B).

The smallest post-run arithmetic verification is a single saved-receipt
collection of all50 draw indices, failures, retained LL/starts, components,
threshold/tie decisions and p/MCSE. It requires no fit. Exception reasons
absent from a receipt remain unverified; fit=None is not attributed to a
previous mechanism. Original failures are preserved and any revised runs
are separate. Fixture agreement claims, especially the binding zero-failure
expectation, remain unresolved if any original/revised refit fails. None
of these checks authorizes opening F237, full-design size claims or
biological replication before the existing scientific conditions are met.

Current-status qualification to the preceding future checklist: the
accelerated parallel fit22042114 terminated during `observed_marginal_fit`,
not at a completed kept marginal result. Its saved `failure.json`, SHA256
`d3e7946dc977906052a34a36a12fc9898ba0d82fd77497a52df41f83f8a7f305`,
records call12, gene2, variance0.0009636462009520182 and
`IntegrationFailure: GL order refinement exceeds unchanged precision allowance`.
It explicitly records no full sum formed and the first failed gene in
ascending order. Root reports the LL refinement3.226e−12 and variance-score
refinement2.757133188e−9: LL is small, but the score exceeds the existing
quarter-gradient budget2.5e−9. No margin, restart or completion claim was
changed. Separate owner-theta4 anchored checks22042346–48 address other
small-variance fixtures and are not duplicated here. Consequently current
P2 profile checks and shared-draw continuation remain conditional on a
credible completed marginal fit; this refusal does not supply one.

### Transparent terminal branch recorder prepared; actual-model check pending

New isolated `repair_p3_marginal_branch_recorder.py`, SHA256
`09f4d6f14074d9d55e2487911c4385a09183fbca0e34a6f9e325828246ece433`,
passed AST parsing. Root explicitly authorized its stdlib-only fake check;
that check executed and passed exact argument-object forwarding, return
object identity, original exception-object identity, nested-install refusal,
restoration, and isolation of metadata/sink failures from fit behavior.
No numerical library was imported and **actual-model compute verification
is not run**. No job, dependency, frozen source or current fit was changed.

The public installer admits only an already guarded original small synthetic
fixture (G8/Z7/P28/173 rows/60 people), individual compute, and the original
archived41eb method source with bound signature `(start,fixed_v)`. Input
ancestry/count admission remains the caller's responsibility; the recorder's
scope declaration is not independent provenance proof. It assigns a temporary
instance method which forwards original positional/keyword calls verbatim,
returns the exact original object, and re-raises the same original exception.
The class-level original fitter definition remains unchanged, but the
**instance method is instrumented** and identity declarations must state
that explicitly. Restoration removes the instance override and refuses to
overwrite an unexpected concurrent replacement.

Each terminal receipt captures supplied start/fixed variance and the actual
returned parameters, variance, LL, full gradient and native convergence
flag, or original exception type/message/limited traceback. Full finite
float values round-trip through JSON; nonfinite values are explicitly
represented, never replaced with finite values. Caller-labelled operation
contexts supply invocation order; warm/cold/restart roles are interpretable
only for actual calls through the guarded original outer fitter. The helper
does not invent retrospective labels, capture every optimizer iteration,
recompute stationarity or prove a global maximum.

Defaults bound storage to16 invocations/512KiB, with strict upper limits.
Overflow or recording/sink errors mark evidence incomplete without altering
fit starts, results, selection or numerical failure counting. The caller
must call `require_complete()` **outside the fit** before declaring full
branch evidence. Its `BranchRecordingError` is an instrumentation error,
not a numerical fit failure. For a future authorized actual-model check,
compare instrumented/uninstrumented original branch returns and selection
at a guarded synthetic case, verify all snapshots and restore, then add
instrumentation-specific method declarations to the runner. No such check
is claimed here, and the owner's exact refused-point arithmetic preparation
is not duplicated.

Root reports sequential22040873 also terminated at the exact same
call12/gene2/variance/full28-parameter proposal and six refused ladders as
parallel22042114. The two original failure records remain intact; this
recorder cannot reconstruct their missing branch terminal histories.

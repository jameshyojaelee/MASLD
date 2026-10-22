# Independent contribution review across MASLD model candidates

Date: 2026-09-30. Codex-owned technical review for Claude and the coordinator.
Read-only assessment of scientific evidence; no outcome access, analysis,
submission, or manuscript adoption was performed.

## Decision

[Inference] Keep multiple candidates eligible. RNA-to-H3K27ac has the clearest
positive development result with an immediate scientific use. Model A has a
distinct inferential question but has not demonstrated calibrated inference or
an independent biological result. Sequence/allele models have strong existing
comparators and clear limitations; they cannot acquire novelty simply through
an adapter or a new liver label. Program-based alternatives must demonstrate
information beyond ordinary co-expression and cell composition.

None currently proves the requested **prespecified, independently validated,
useful new model**. This is not a numerical ranking across tasks: prediction
skill, caQTL beta MSE, and calibrated confidence bounds answer different
questions. A candidate should advance because it resolves its own scientific
question on independent evidence, rather than because another track is slow.

## Candidate comparison and exact missing evidence

| Candidate | Researcher question and units | Evidence currently supporting it | Missing evidence that changes eligibility for the paper's new-model claim |
|---|---|---|---|
| Model A, proportional cis-effect remodeling/persistence | Does a fixed eligible gene family show common proportional change per cross-sectional recorded stage category, or exclude a prespecified magnitude? Gene-specific baseline natural-log allelic effects; participant-cluster inference. | Independent implementations agree on fixed-vector likelihoods, gradients and the observed synthetic fit. Planning uses 558 development participants; expected eligible genes are 379/364/255/166 under the four candidate public-effect thresholds. No MASLD allelic result is evaluated in this review. | Complete independent draw/test agreement; full ancestry-adjusted design information and final threshold choice; null size, corrected-effect bias, interval coverage and power under the approved selection/label/technical plants; independent-participant confirmation after authorized evaluation. A nonsignificant stage test is insufficient for persistence. |
| RNA-to-H3K27ac, global plus local | Can ordinary liver RNA recover measured chromatin differences among people beyond a fixed covariate predictor? Predicted residual H3K27ac log2 CPM; held-participant squared-error skill. | In 99 GSE267145 participants, skill .14133 versus .01971 for the tuned histology/sex/technical comparator; increment .12162 [.09542, .14682]. Local versus global increment .02214. Matching gives 22 disjoint exact-grade/sex/fold pairs, mainly women with fibrosis 0. Independent reconstruction supports the arithmetic. | Independent-source paired RNA/H3K27ac accuracy for a fixed scorer, transformations, regions and endpoint. GSE128072 has 17 shared processed donor IDs, but TPM inputs, hg19 coordinates and normal-liver restriction require resolution before admission. GSE296875 ATAC is a different-assay transport result. Missing-input accuracy for the full released scorer and individual prediction uncertainty remain unproved. |
| Signed endogenous liver caQTL prediction | Which variant-centered native source ALT-dosage associations are predicted, and which peaks merit assay follow-up? Association beta, not isolated editing effect. | On 24,430 matched Currin development leads, MSE .300161 native 1 Mb, .302252 adapter and .290705 combined. Combined relative gain 3.15% misses the fixed 5% complexity margin; retain calibrated native 1 Mb. Development predictions are numerically reconstructed. | Independent source with authoritative donor membership, allele/reference alignment, compatible endpoint and checkpoint exposure. Wenz admission is unresolved: significant leads and RASQUAL Pi are not Currin beta or an all-tested variant denominator. Runtime variation prevents claiming arbitrary fresh inference reproduces archived scores. Fold 0 is spent, not a new selection/confirmation set. |
| Sequence plus molecular context | Does context improve donor differences or allele-by-state effects beyond RNA-only, sequence-only and additive models? Native molecular endpoint, donor/region/source units distinguished. | The bounded 1,024-region sequence-by-RNA correction has original-target skill increment −.0000159 [−.001520, .001375], below the fixed +.01 margin. Earlier development comparisons expose the strength of RNA-only controls. | A substantively different prespecified interaction with compatible measured labels and independent evaluation. The failed fitted interaction is not evidence that every context model fails, and it cannot be relabeled as a success. Held-region/held-source effects require those evaluations; held participants at trained regions do not supply them. |
| Fixed programs with genetic or cross-assay information | Does a defined program add disease or regulatory information beyond expression, composition and matched co-expressed sets? Donor effects and genome-wide heritability/linking endpoints are separate. | Fixed117 remains unchanged. Existing cross-assay module candidate has reproducible associations but zero discovery modules exceeding the corrected connected-set null and one assay cell exceeding it. It is deferred by user decision. Saturation profiles describe model sensitivity, not measured donor effects. | For a new eligible adaptation, a prespecified fixed-program endpoint, genome-wide GWAS and ancestry-compatible LD where required, baseline genomic/composition annotations and matched program controls. Independent disease/process discrimination or molecular prediction must exceed those controls. Reopening the deferred module candidate is not authorized by this review. |

Numbers and their candidate status come from [RESULTS](../../RESULTS.md),
especially the two-model development and cross-assay-module sections;
[the fixed two-model design](../../../scripts/analysis/masld_two_models/design.json);
[Model A synthetic evidence](2026-09-30_codex_bfix_results.md); and
[Model A planning information](2026-09-30_codex_gateA.md). Their arithmetic was
not rerun here. Live job status is intentionally omitted.

The .14133 chromatin skill uses the full held-participant development
comparison. The separate within-histology matched-pair skill is .11070
[.05349, .16267] for local versus .08859 for global RNA, on 44 people in
22 disjoint pairs. These endpoints and populations are not interchangeable.

## Novelty and comparison fairness

RNA-derived genotype/ASE inference and context-dependent ASE have precedents;
[the bounded Model A literature review](2026-09-30_codex_modelA_contribution.md)
documents QuASAR, BaseQTL, EAGLE, ASEP, TReCASE, BBmix and aiQTL. Model A's
possible distinction is the shared proportional estimand plus bounded
persistence with selection/label-aware inference. [Unverified] Priority for
that combined procedure is not established. Compare per-gene methods on their
per-gene targets; compare reduced shared-slope models on common κ and its
bounds. Observed-genotype simulation methods are oracle comparisons, not
inferior methods merely because their inputs are unavailable in this study.

RNA-driven epigenomic imputation and sequence/context conditioning also have
precedents. Corgi integrates sequence with trans-regulator expression, and
Corgi+ adds local RNA and training-context baseline information for epigenomic
imputation. Thus “predicting chromatin from RNA” alone is not a new capability.
The possible MASLD contribution is independently transferable participant
variation beyond a credible comparator, with a defined application and limits.
The source model's cross-cell-type performance does not validate our
cross-participant model. [Corgi primary article](https://www.nature.com/articles/s41467-026-75527-2).

Similarly, sc-linker already connects single-cell programs, epigenomic
SNP-to-gene maps and genome-wide GWAS statistics. A MASLD adaptation needs a
new supported relationship or demonstrable benefit of its links beyond the
original strategy and matched controls; attaching selected COLOC rows to
programs is not equivalent to its genome-wide inference.
[sc-linker primary article](https://www.nature.com/articles/s41588-022-01187-9).
This literature inspection is bounded and does not establish an exhaustive
priority claim. No software from those papers was installed or run.

Use identical eligible rows, source assay units and evaluation partitions for
every within-task comparison. Give sequence-only methods their native context
length as well as equal-length controls; otherwise adaptation and receptive
field gains are confounded. A significance-selected caQTL source does not
estimate performance on arbitrary variants. Histology-plus-technical ridge is
a useful stated comparator, but beating it does not establish conditional
independence from arbitrary nonlinear histology/composition relationships.

## Independence requirements that cannot substitute for one another

1. **Numerical independence:** a second implementation can reproduce a score,
   derivative, prediction table or fitted test. It detects implementation
   errors. It cannot validate a shared incorrect assumption or biological
   transport. The ongoing synthetic comparison is in this category.
2. **Development evaluation:** held participants test new people conditional
   on the source design; held genomic units test different sequence exposure.
   Five folds have overlapping training sets and are not five independent
   cohorts. Fixed-prediction bootstrap intervals omit refitting/selection
   uncertainty unless those steps are repeated.
3. **Biological confirmation:** distinct authoritative participants, fixed
   choices before evaluation, compatible measured endpoint, and documented
   training/source exposure. A newly named table or independent code review
   of already inspected participants does not provide untouched confirmation.
4. **Assay transport:** H3K27ac, ATAC, RNA expression, reporter activity and
   caQTL association coefficients retain their native endpoints. Success on
   one does not validate another. Same-nucleus RNA/ATAC coupling has shared
   biological/technical causes; it does not by itself establish regulation.

The current five authorities also qualify GSE193066 identity and previous
evaluation exposure. It cannot be recycled as untouched external confirmation.
No real allele outcomes, restricted genotype directory, protected evaluation
outcomes, experimental-paper data or participant depth files were read for
this assessment.

## Criteria for choosing the model contribution

A candidate qualifies only when its evidence answers all of the following:

- **Distinct question:** a concise estimand and user input/output explain what
  a researcher learns that the available resource or strong comparator does
  not already supply. Architectural complexity is not the criterion.
- **Prespecification:** the next independent evaluation fixes the target,
  eligible universe, transformations, comparator, useful-effect margin,
  uncertainty calculation, and complete hypothesis family before outcomes.
  Adaptive development remains identified as such. New parallel tracks do
  not replenish spent evaluation sets or waive existing stopping rules.
- **Scientific correctness:** for inference, calibrated size and coverage
  precede power claims; for prediction, target-native held-unit loss and
  uncertainty precede causal interpretation. Bounds and failures are retained.
- **Independent value:** a permitted untouched source or independent
  participant evaluation supports the intended application. A within-source
  split may be valid participant confirmation but does not establish
  cross-study/assay transport. Endpoint mismatch is resolved before scoring.
- **Usable manuscript result:** the exact fitted recipe generates a paper
  analysis, takes documented assay inputs and reports supported contexts,
  units and failure states. Invocation examples are necessary but are not
  accuracy evidence. Neither expected citation counts nor higher-impact
  acceptance can be verified as an analysis outcome.

Do not choose the winning track by its smallest observed p-value. Each track
retains its declared within-task multiplicity; any final pooled claim based
on searching several endpoints must account for that selection or use a fresh
fixed evaluation. A null can establish a useful limit with a sufficiently
precise prespecified bound. It does not automatically constitute a newly
validated reusable model.

## Concrete disposition for the coordinator

[Inference] Prioritize independent-source admission for RNA-to-H3K27ac in
parallel with Model A calibration and full-design information. Preserve the
simple native allele predictor while assessing source/endpoint eligibility.
Admit a new sequence/context/program experiment only when the actual endpoint,
independent evaluation unit and strongest comparator are available. The exact
decision-changing gaps are in the table; no extra architecture search is
justified by this review alone. Model B's closed aspect analysis and the
deferred cross-assay-module candidate remain closed/deferred.

Checks performed: inspected repository instructions, five authorities,
integrated AG program, fixed two-model design/README and Codex Model A
technical reports; checked the two linked primary-literature sources. No
numerical analysis or new compute job was run. Only this technical review
file was written; manuscript authorities and scientific results are unchanged.

## Comparator feasibility: regular Corgi and Corgi+, bounded follow-up

This section checks upstream public code and local admission records, not
evaluation outcomes. The models must remain separately named.

**Regular Corgi can supply a development comparator, with an explicit
project mapper and calibration. Corgi+ is not ready for a native, portable
comparison in the inspected public surface.** On 2026-09-30, upstream
`predict.py` still makes `corgiplus_pretrained` a placeholder: all three
prediction methods raise `NotImplementedError`. The Corgi+ configuration
still refers to original `/project/` arrays. This directly rechecks the
central loader/input limitation in the dated local disposition; it does
not certify that every possible upstream auxiliary artifact is absent.
[Current upstream prediction source](https://raw.githubusercontent.com/ekinda/corgi/main/corgi/predict.py),
[Corgi+ configuration](https://raw.githubusercontent.com/ekinda/corgi/main/corgi/config_corgiplus.py).

| Requirement | Regular Corgi | Corgi+ |
|---|---|---|
| Biological inputs | 524,288-bp one-hot DNA and a 2,891-trans-regulator expression vector. | Same global inputs plus local stranded total-RNA coverage, 22 genomic training-context mean tracks, and two local RNA-minus-mean tracks. The local checkpoint record identifies 26 auxiliary channels. A gene-by-sample TPM/count matrix cannot supply these coverage tracks. |
| Native output | Central 393,216 bp at 64-bp resolution, 6,144 bins and 22 output channels. H3K27ac is zero-based channel 7. | Same assay/channel family; local record describes epigenomic imputation, with supplied RNA channels excluded as outcomes. |
| Normalization | Current helper sorts within-vector expression ranks onto packaged `tf_reference.npy`; named vectors are reordered to packaged gene IDs. Missing named genes are silently zero-filled by the official helper, requiring an explicit project missing-gene rule. | Global-context uncertainty is shared; local RNA/baseline/delta transforms and coordinate alignment are additionally required. |
| Available weight records | Local record identifies `corgi_model.pt`, 2,350,451,066 bytes, Zenodo 18630048/20686797. Existing isolated L40S compatibility was checked previously, not in this review. | Local record identifies `corgiplus_model.pt`, 2,384,030,158 bytes; weight availability alone does not establish a complete biological-input inference path. |
| Remaining feasibility limit | Checkpoint-bound reference/source, native numerical parity, exact assay-scale inversion and project context mapper are unresolved in the local record. | Original auxiliary tensors and transforms, checkpoint-bound reference/source, portable strict inference and official biological parity fixture remain unresolved; current helper is still unimplemented. |

Current geometry/channel values are directly visible in
[upstream constants](https://raw.githubusercontent.com/ekinda/corgi/main/corgi/constants.py),
and context ordering/missingness in
[upstream regulator helper](https://raw.githubusercontent.com/ekinda/corgi/main/corgi/trans_regulators.py).
Weight identity and the additional local constraints above come from
[regular checkpoint declaration](../../../Analysis/MASLD_Model_Benchmark/config/artifacts/models/corgi/checkpoints.json)
and [Corgi+ declaration](../../../Analysis/MASLD_Model_Benchmark/config/artifacts/models/corgi_plus/checkpoints.json).
These declarations report historical checks; this review did not load or
rehash weights. The official repository README links pretrained weights;
the current reproduction README still leaves the reproduction-data Zenodo
URL unspecified. [Upstream README](https://github.com/ekinda/corgi),
[reproduction README](https://raw.githubusercontent.com/ekinda/corgi-reproduction/main/README.md).

### Exposure remains unresolved for the two paired liver cohorts

The article describes broad functional-genomics training and separate
context/sequence evaluations, while the local record lists ENCODE, FANTOM5,
Tabula Sapiens, CATlas and two healthy 10x sources. Its context preprocessing
record describes log/quantile normalization and cross-assay ComBat; the
released within-vector rank helper is not that full training procedure.
[Primary article](https://www.nature.com/articles/s41467-026-75527-2),
[local exposure crosswalk](../../../Analysis/MASLD_Model_Benchmark/config/artifacts/models/corgi/development_crosswalk.json).

[Unverified] Neither GSE267145 nor GSE128072 has an authoritative
checkpoint-bound donor/run exposure classification in the inspected
crosswalk: neither is among its individually reviewed project sources.
Literal-accession searches did not resolve this. The exact training exposure
manifest, including donor/run aliases and normalization/reference inputs,
is missing from the inspected records. Thus neither cohort can be called
unseen by Corgi from these records. A normal-liver source such as GSE128072
also needs aliases checked against the declared normal-tissue corpora.
The article pages intermittently returned access errors/CAPTCHA during
follow-up; lack of a searchable accession there does not prove exclusion.

### Fairness to the participant residual target

Native H3K27ac track output is not the project's residual log2-CPM target.
The local declaration describes assay scaling/clipping and 64-bp means,
without a complete numeric inverse transform. A valid comparison must fix
region/bin aggregation and learn any target-scale calibration using training
participants only, then score the identical held participants/regions and
training-defined residual target as RNA RRR. Both raw-profile and
participant-difference performance should be reported: high regional activity
correlation can coexist with little participant-specific information.
Mean-context and shuffled-context controls test whether RNA context helps.
If a new context mapper or head is selected, its result is an adapted Corgi
development comparison, not native zero-shot performance.

The current RESULTS authority explicitly corrects the historical Corgi
comparison: its scaling was fitted to held outcomes and its RNA comparator
trained on every output region. Therefore those old differential skills
cannot establish superiority on training-calibrated donor-by-region error.
The five-fold .14133 RNA result should not be compared directly with an old
Corgi number on another transformation/region universe.
[RESULTS correction](../../RESULTS.md#regulatory-model-evidence-corrections-2026-09-15).

**Disposition:** use regular Corgi as the bounded comparator only with its
actual inputs, declared mapper/calibration and matched endpoint. Do not
replace the unavailable Corgi+ auxiliary inputs with gene abundances,
zero-filled tracks or query-cohort measured chromatin means. A training-fold
baseline with new transforms would define an adapted comparison requiring
its own fixed recipe. The review does not authorize implementing it.

Exact local Corgi+ disposition:
[release_disposition_20260824.json](../../../Analysis/MASLD_Model_Benchmark/config/artifacts/models/corgi_plus/release_disposition_20260824.json).
Its read-only check script is
[audit_corgi_plus_disposition.py](../../../Analysis/MASLD_Model_Benchmark/scripts/audit_corgi_plus_disposition.py);
it was inspected, not executed. No weights, outcome tables, auxiliary
biological arrays, GPU inference, dependency installation or jobs were used.

## Donor-contrast development contract: review before fitting

Reviewed proposal: the existing 99 development participants and five outer
folds; each of the four other folds used once as inner validation; the
unchanged 22 disjoint exact-histology/sex evaluation pairs; λ in {0, .1, 1};
graph GLS contrasts; same-input weighted ridge and fixed-rank-12 RRR controls;
the existing outcome-blind 1,000-region panel; inner paired-MSE selection with
a profile-MSE safeguard. No fit or outcome inspection was performed here.

**Verdict: a fair bounded development experiment after the definitions below
are fixed.** It does not supply fresh confirmation, even with nested fitting:
these participants, pairs and broader source outcomes have already been
inspected. The new pair-skill increment .01 is a planning rule for this
experiment, not an independently justified useful biological magnitude. The
existing +.01 original-target skill margin remains the criterion for a
full-profile complexity claim; changing the denominator cannot satisfy it.

### Essential definitions before the fit

| Item | Required contract change or clarification | Why it changes the scientific comparison |
|---|---|---|
| Model identity | State the proposed feature map, response, rank restriction and penalty separately from weighted ridge. If both use identical linear features/penalties and the same unconstrained objective, they are the same estimator. | Rebranding the same fit as a contrast method cannot establish added capability. Fixed-rank RRR can be a distinct constrained control. |
| Graph weighting | Fix incidence matrix B, covariance assumption, pseudoinverse/rank tolerance and normalizers for donor and contrast losses. Use graph covariance; do not treat edges as independent people. | With independent equal-variance donor errors, edge covariance is BBᵀ; graph GLS yields within-stratum centering, not extra observations. λ has no interpretable comparison unless relative loss scales are fixed. |
| Fold preprocessing | Refit RNA transforms, response concentration regression, centering, feature selection/PCA and any estimated GLS covariance inside each inner training split, then apply unchanged to that inner validation fold. Refit selected recipe on outer training only. | Outer-training residuals/preprocessing that already used inner validation outcomes invalidate inner selection. Query outcome descriptors may define the response as in the existing target, but cannot enter RNA predictions unless available for intended use. |
| Pair formation | Retain all 22 outer pairs with their exact metadata/hash rule. Form training/inner-validation pairs after splitting donors and use the same deterministic rule; fix handling of empty/zero-reference inner sets without reading outcomes. | Complete graph edges are dependent, while outer disjoint pairs are evaluation units. Re-pairing based on measured differences would change the target population. |
| Selection | Pool inner validation pair SSE with explicit weighting, rather than average four ratios with different denominators; select joint (λ, penalty) by one deterministic rule. Specify penalty grid, tie rule and numeric profile-MSE safeguard before fitting. | “Profile safeguard” without a comparator/tolerance is an adjustable rule. The profile criterion and tie preference must not be tuned on outer errors. |
| Controls | Give each arm the same RNA inputs, targets, panel, folds, penalty opportunities and target transformations. For rank-12 RRR, keep rank fixed at 12 and handle insufficient inner rank by a fixed rule. Identify whether λ is tuned for each eligible arm or held at zero. | Otherwise an improvement can reflect unequal tuning or response exposure rather than the contrast objective. Do not use old released predictions covering another region universe as the matched control. |
| Acceptance and family | Define pair-skill gain against each named control or a control selected using inner folds only. Fix the complete family of primary paired contrasts and its Holm correction, and distinguish pair and original-target profile scales. | Choosing the weaker control or best outer contrast after evaluation inflates apparent success. Inner λ choices are tuning alternatives; they need not become separate primary tests when only the selected recipe is evaluated. |

The existing `within_histology.py` verifies pairs share an outer fold and
uses the training-defined response transformation. Its primary skill is
1 − SSE_prediction_difference / SSE_measured_difference, with zero predicted
pair difference as reference. Retain that reference in the new paired
endpoint, and report gain Δskill = (SSE_control − SSE_candidate) /
SSE_measured_difference. Report full-profile skill with its different
training-mean denominator separately.
[Existing pair definition](../../../scripts/analysis/masld_two_models/within_histology.py),
[metadata-only census implementation](../../../scripts/analysis/in_biopsy_allelic/rec_codex/chromatin_pair_census.py).

### Graph GLS interpretation

Let R = Y − f(X), with one row per training participant. For an incidence
matrix B joining participants within the fixed strata, the contrast residual
is BR. Under independent equal-variance donor errors its covariance is
Σ = BBᵀ, which is singular for redundant edges. Then

`tr[(BR)ᵀ Σ⁺ BR] = tr[Rᵀ P R]`,

where `P = Bᵀ(BBᵀ)⁺B` is the orthogonal projector onto within-stratum
contrasts. In a complete connected stratum of n people, P is
`I − 11ᵀ/n`, with rank n−1. Thus the combined profile/contrast loss is
weighted donor residual loss, with W proportional to I plus λP after the
specified normalization. This algebra assumes the stated covariance; it
does not prove heterogeneous chromatin errors satisfy it. If training
variances are estimated, document that different covariance assumption and
fit them inside splits. The derivation is mathematical inspection, not an
executed numerical check.

The defensible contribution being tested is whether emphasizing identifiable
within-histology donor differences improves held-donor contrasts enough to
justify the extra recipe, while preserving profile performance. It is not a
claim of new biological replication or information added by graph edges.

### Uncertainty and decision language

Paired bootstrap should resample the 22 disjoint pairs jointly for all arms
and recompute SSE-ratio contrasts with a recorded seed; chromosome/block
resampling is a separate dependence sensitivity. The 44 people are mainly
women with fibrosis 0 and have no fibrosis 2–4 pairs. This restricts any
within-histology finding to the represented low-fibrosis population.

Conditional resampling of saved predictions omits retraining, nested
selection and shared-training dependence. It can support descriptive
development intervals and the stated comparative point estimates. It cannot
make the five overlapping training folds independent replicates or justify
fresh-confirmation/type-I-control claims. If centered conditional-bootstrap
p-values and Holm values are computed, label them development summaries;
Holm does not repair a nonvalid sampling approximation. No significance
threshold should turn an interval crossing the .01 planning margin into a
claim of a useful gain above that margin. Report the point, interval, margin,
each control and profile safeguard result without declaring independent
superiority.

Before fitting, fix whether .01 is a point-estimate eligibility rule or a
lower-confidence-bound criterion; these are different decisions. Do not
silently reuse historical within-pair p-values or evaluate each λ as a new
primary claim. Preserve a null/negative result and do not expand the λ grid
or rematch pairs when a gain fails the declared rule. This review authorizes
no new analysis; it identifies the definitions needed for the requested
bounded development comparison.

## Target-readout comparison: variant center versus measured peak

**Verdict:** a bounded assay-alignment development test is scientifically
eligible once peak coordinates and the exact comparison are fixed. Integrating
the source accessibility peak may better align the sequence readout with the
caQTL association phenotype. That hypothesis is distinct from Model A and
from a new learned model. It must earn predictive support; genomic proximity
or correct sign in an illustration is insufficient.

The current experimental native interface extracts a 1,048,576-bp REF/ALT
sequence pair and predicts 1-bp ATAC tracks. Its score is the mean of three
named liver tracks, each using `log2(1+sum ALT) − log2(1+sum REF)` over a
501-bp variant-centered mask. The adapter pools a central 384-bp embedding
interval and does not consume the measured peak bounds. Thus the existing
variant–peak query identifies the measured association, without actually
targeting that peak in either model.
[Existing native inference](../../../scripts/analysis/masld_two_models/variant_inference.py),
[exact-pair query contract](../../../scripts/analysis/masld_two_models/variant_effect_query.py).

### Fix these points before feasibility calls or outcome comparison

- **Coordinate contract:** use authoritative source peak start/stop, verify
  genome build and closed/half-open convention, and freeze conversion to the
  model's 0-based output coordinate axis. Do not infer the convention from
  column names. Verify REF/ALT alignment, chromosome and full peak coverage.
  An interval outside the returned native output is untestable under this
  geometry; do not silently clip it. A SNP outside its associated peak is
  allowed only if explicitly included in the planned geometry class.
- **Matched readouts:** calculate both old 501-bp and full-peak masks from
  exactly the same REF/ALT prediction objects in the same process, same
  ontology tracks, output request, checkpoint, card and sequence window.
  Preserve the old per-track sum/log/difference and track averaging order;
  change only the mask. Changing track aggregation, pseudocount, normalization,
  TTA or native output request would add separate choices.
- **Fixed feasibility population:** choose either 32 or 64 development
  variant–peak pairs before calling the model, using recorded width,
  variant-to-peak distance/overlap and deterministic hash ties. Cover short,
  broad, distal and edge geometries where available. Inspect metadata only;
  neither beta/p-values nor old/new prediction magnitude chooses examples.
  Predetermine refusals and counts for sequence/coverage failures.
- **Population and calibration:** the proposed 25,748 significance-selected
  folds-1–4 leads are the source population, not a guaranteed common eligible
  set. Report full-mask and 1-Mb sequence exclusions, then compare both scores
  on identical retained rows. The historical matched population was 24,430;
  it cannot be assumed unchanged. Fit each score's calibration using only the
  other development chromosome folds, with the same fixed slope/intercept
  rule and tuning opportunities. Exclude fold 0 before parsing outcomes.
- **Comparators and uncertainty:** the primary simple comparator is calibrated
  native 1 Mb with its original mask, evaluated freshly on the matched rows.
  A 2-kb native arm alone is insufficient. If adapters or combinations enter,
  retain the best eligible simple baseline and the fixed 5% relative-MSE
  complexity requirement; this readout comparison does not relax it. Use
  paired whole-chromosome resampling, record chromosome/row counts, and label
  intervals conditional on saved predictions/calibrations. Do not call
  historical 1-Mb bins independent LD units or four folds four cohorts.

Runtime discrepancies previously observed between fresh native processes
make the same-output comparison essential. The old archived scores can be
reported as a reproduction sensitivity, but they cannot be one side of the
primary readout comparison. The 32/64 cases establish geometry, numerical
reproducibility and resource feasibility; they do not estimate full-population
accuracy or justify outcome-selected readout tuning. Preserve all feasibility
case results and failures before any later prediction-error calculation.

### Existing execution boundary

The integrated AG program §19.4 authorizes individual SLURM jobs with local
native 1-Mb Currin scoring and preserves a 500-GB new-file bound, protected
Wenz/GSE289173 outcomes and the manuscript firewall. The later two-model
design explicitly prohibits new fold-0 selection, fitting and evaluation.
The local `NativePredictor` requires a SLURM GPU allocation and L40S; the
existing 32-case cost script reads metadata without molecular outcomes for
its feasibility manifest. These support a bounded local geometry pilot,
subject to the coordinator's actual submission/resource accounting. They
do not authorize hosted API calls, a new runtime, weight downloads or
protected evaluation access merely because another track is running.
[Authorized program](../ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md#194-authorized-initial-week-2026-09-15),
[fixed allele design](../../../scripts/analysis/masld_two_models/design.json),
[existing metadata-only feasibility recipe](../../../scripts/analysis/masld_two_models/inference_cost.py).

The older AlphaGenome native-task disposition concerns a different,
unregistered `sequence_native_regulatory` tournament task and explicitly
leaves its `variant_to_regulation` standing untouched. It must not be used
to invent a blocker for the already executed Currin lane. Existing local
checkpoint terms identify AlphaGenome as a restricted noncommercial
comparator; this review does not reclassify redistribution rights or adopt
a new champion.
[Task-scoped disposition](../../../Analysis/MASLD_Model_Benchmark/config/artifacts/models/alphagenome/native_task_disposition_20260825.json).

**Scientific decision:** if calibrated full-peak integration reproducibly
reduces matched development MSE, the Resource gains a better assay-aligned
native prediction recipe with a stated eligible population. It gains no
independent biological confirmation or isolated allele-editing effect.
If it fails, retain the original readout and report the bound on improvement.
Neither outcome by itself satisfies the requested new-model contribution;
it can clarify the strongest baseline a genuinely new model must beat.
This section inspected code, designs and permission records only. No outcome
tables, APIs, model calls, GPU work, jobs or installations were used.

## Missing-program contribution: scientific admission, 2026-09-30

**Decision [Inference].** Predicting the unmeasured contribution to a frozen
RNA program could be a useful measurement tool. A masked bulk-RNA test can
validate reconstruction of an RNA score under gene omission. It cannot
validate spatial detection, protein abundance, program localization or
biological activity in an assay lacking those genes. The accepted CosMx
complete-program observability remains 12/117 (RESULTS); predicted values
must not change that observed-assay count.

The current GSE268273 fixture's eight fully covered programs out of 117 are
a statement about its filtered gene axis. Until the full source axis and
gene-ID mapping are checked, they do not establish that the source assay
failed to measure the other 109 programs. Record, separately, source absence,
fixture filtering, unresolved gene mapping and artificial panel masking.
Neither zero filling nor retained-weight renormalization creates a full
reference score for absent members. All programs lacking a complete target
remain untestable for this reconstruction endpoint; no coverage threshold
should be relaxed after seeing model error.

### Endpoint and comparison that avoid shared-gene agreement

For each fixed program and fixed panel, define the primary target as
`S_missing = sum(unobserved canonical weight × transformed RNA expression)`.
Inputs contain only measured genes; the target contains only excluded genes.
The transformation must be defined explicitly and fitted to development
training participants only. Missing contributions from low-variance genes
can have small absolute error without useful recovery; report training-scale
variance, missing weight, missing-member count, and skill against the
training-estimated missing mean. Full-score error is a secondary practical
endpoint. Correlation between reconstructed and full scores is not primary
evidence, because their observed contribution is identical by construction.

The minimum competitive set on identical measured genes and participants is:

- Training-estimated missing mean, equivalent to measured contribution alone
  only when the training-centered missing contribution has mean zero.
- Training-fitted scalar calibration from observed program score to missing
  contribution. This tests whether multigene inputs add information.
- Direct regularized ridge prediction of missing contribution, with penalty
  selected in the same study-blocked development splits as the candidate.
- Gene-wise regularized imputation followed by unchanged canonical weights,
  with matched training exposure and tuning allowance.

A proposed ridge predictor is itself the third arm; renaming its target does
not provide a new architecture. Its contribution would have to be a useful,
validated measurement recipe and explicit refusal/uncertainty rules. If a
more complex candidate is introduced, compare it directly with this ridge
arm. Match any auxiliary genes, composition inputs and clinical covariates
across arms. Genes outside a program can be legitimate measured predictors,
but giving them only to the candidate creates a different input comparison.
Random or matched gene sets can probe specificity; they cannot replace the
strong predictive baseline.

### Independence sufficient for the stated use

An external complete-RNA cohort with separately resolved participants can
test a locked score-reconstruction recipe under a prespecified artificial
mask. Keep technical runs and repeat biopsies together. Verify exposure at
the participant and source level, and include any studies used to discover
the frozen programs in the exposure record. The target uses different genes
from its predictors, but that is algebraic separation, not a new biological
replicate. Agreement between two implementations is numerical verification,
not transport validation.

The existing assessment records GSE268273 as project-exposed external
development and unresolved universal donor disjointness for the candidate
receiving sources. Such cohorts can support a qualified external development
experiment after freezing this recipe; they cannot be described as untouched
confirmation. A new participant-disjoint complete-RNA receiving cohort,
evaluated once after model, transformation, mask, eligibility and endpoints
are fixed, would supply meaningful independent validation of the RNA
measurement use. A matched targeted/full-RNA measurement on new participants
would additionally test actual targeted-assay transport. That matched target
is not presently established. No new source acquisition is authorized here.

### Evidence still required before a fit or contribution claim

Fix the exact mask and complete-target program family from gene axes, not
prediction error. Preserve refused programs in the 117-row registry. Fix
training/cohort splits, input transformations, candidate and comparators,
aggregate error weighting, and a useful improvement margin before fitting.
The useful magnitude is currently unspecified; a borrowed numerical margin
from chromatin or allele prediction has no scientific justification here.
Interval claims require an explicit separate calibration source and a
coverage/width endpoint; exchangeability under study shift is not assured.
Do not refit normalization or calibration on receiving target values.

Evaluate a prespecified aggregate missing-contribution endpoint first; if
making individual-program claims, correct the entire admitted program family
and name that family. Cohort-blocked training selection and donor-level
uncertainty do not make a handful of study sources numerous independent
cohorts. Gains concentrated in almost-constant missing components or failure
against scalar calibration weaken practical value even when full-score
correlation is high. No prediction results or scientific usefulness have
been established by this admission review.

This section used existing metadata descriptions, assay-coverage authorities
and design documentation only. No expression/outcome arrays, model imports,
jobs, installations or new literature searches were used. The concrete
full-source-axis assessment remains with the parallel program-hypotheses
track.

### Review of the concrete development design

The parallel report now specifies 69 fit, 20 calibration and 20 internal-test
participants from the exposed GSE268273 source; 714 measured-panel inputs;
eight currently complete program targets; and rank-1/2/4 reduced-rank ridge
against same-input unconstrained ridge. These are development splits within
one source, not independent study confirmation. The fixed metadata hashes
and refitting of transformations in inner training are appropriate. Neither
program coverage nor the source participant count was independently
recalculated for this contribution review.

The concrete design correctly prohibits full-library totals in predictor
features and correctly labels full-score reconstruction as an oracle-exposure
diagnostic. This limitation changes the scientific use: absent an independent
exposure measurement or an assay-native adapter, predicting a reference-CPM
missing contribution does not provide a deployable complete score for a
targeted panel. The contribution remains a conditional RNA measurement
experiment. Disjoint target genes avoid direct shared-score agreement, but
the complete-RNA normalization denominator can still couple genes; successful
prediction is not evidence of specific biological regulation.

Two details require an explicit definition before fitting: inner tuning
must use program-standardized error with variance estimated from each
inner-training split, and the primary inferential comparator must be named.
Candidate-versus-selected-simple-baseline and candidate-versus-unconstrained
ridge are two contrasts if both carry separate inferential claims. Keep one
primary contrast or specify their joint family. Same-penalty gene-wise ridge
followed by fixed weights can equal direct target ridge algebraically; it
does not constitute a separate independent comparator victory. The design
already acknowledges this equivalence and the currently unjustified
application tolerance. A narrow confidence interval above zero would not
repair that missing practical-use criterion.

## Strategic comparison after the bounded findings

**Recommendation [Inference].** The strongest scientifically distinct next
candidate among these three tracks is a **variant-and-region-conditioned
accessibility association predictor**. Its use is predicting the signed
association for a supplied variant and measured peak, including cases where
the peak is outside the variant-centered mask. It introduces a defined target
that the current variant-only adapter does not accept. This is a candidate
capability, not demonstrated architecture novelty or independent validation.
The paired-mask pilot is the necessary baseline experiment before designing
any learned extension; keep its pending status and fixed population unchanged.

| Track | Current scientific disposition | Next decision justified by current evidence |
|---|---|---|
| Weighted within-histology RNA→H3K27ac | Stopped at the fixed 1% inner-profile safeguard; only one outer fold completed. This does not establish that all contrast models fail, but supplies no full candidate comparison. | Preserve refusal. Do not change the safeguard to obtain a result. Existing positive unweighted RNA→chromatin development prediction remains distinct evidence, not validation of this extension. |
| Missing weighted program contribution | Eight complete targets in the restricted fixture; predicting missing full-reference CPM contributions does not yield a deployable full score without independent exposure information. | Retain completed metadata admission. No count fit is justified until a receiving measurement and useful error tolerance exist. Seven liver-cancer T/B donors are a possible technical-transfer source; at most one non-tumor donor and unresolved deposited crosswalk prevent treating it as MASLD validation. |
| Variant-and-region-conditioned accessibility | Target geometry is admitted for the bounded pilot; accuracy and practical benefit remain unmeasured. Same-prediction center/peak readouts isolate assay alignment. | Finish the already authorized pilot. If feasible, compare the readouts on the same eligible development rows with training-fold calibration. A learned candidate would then have to beat the strongest simple readout under the unchanged complexity rule. |

The RNA→H3K27ac question has the strongest existing positive development
signal, but Corgi already establishes RNA-to-epigenome imputation as prior art,
and the newly tested contrast objective has no passing result. Calling the
retained recipe a new model would not fix the missing compatible independent
endpoint. The program track now exposes a practical input/target mismatch.
The target-region track has the clearest still-unresolved assay question that
can be answered by the authorized experiment. A full-peak readout alone would
be a useful assay-alignment result, not a newly validated learned model.
[Corgi primary article](https://www.nature.com/articles/s41467-026-75527-2),
[current alternative-model evidence](2026-09-30_codex_alternative_models.md),
[program measurement evidence](2026-09-30_codex_program_model_options.md),
[target-region admission](2026-09-30_codex_target_region_admission.md).

### Exact independent-validation prerequisite for the recommended capability

A receiving liver-accessibility study must provide an authoritative final
participant roster, exposure/non-overlap assessment, reference build and
REF/ALT harmonization, exact peak intervals, and signed association effects
with a documented measurement definition. For **arbitrary supplied-target**
utility, the receiving denominator must include tested variant–peak pairs
beyond significant lead associations, with measured null and alternative
associations retained. Missing or untested pairs cannot be labeled zero.
Significance-selected leads can validate only a selected-association use,
even if the model accepts arbitrary coordinates.

Freeze the candidate, input geometry, region eligibility, comparator,
training-derived calibration and metric before those receiving effects are
accessed. Preserve incompatible/missing rows and use identical admitted
variant–peak rows for all predictors. The strongest baseline is the better
development-selected simple native center/peak readout, with the current
calibrated native 1-Mb route and existing adapter comparisons retained as
specified; this review does not change their primary comparator or threshold.
Keep uncertainty at the appropriate study/genomic dependency level. Many
variant–peak rows from the same participant cohort are not independent human
replications. Association accuracy does not establish isolated editing
effects, regulatory-region–gene links or MASLD causal effects.

Wenz's currently inspected significant leads and RASQUAL Pi cannot supply
Currin-beta MSE validation or an all-tested arbitrary-target denominator.
An allele-aligned native-Pi direction/rank endpoint could be a separately
prespecified transport test after admission, but would not independently
validate the original beta-scale predictor. No currently admitted receiving
artifact meets the full prerequisite. Resolving this specific source/endpoint
contract is required before a learned target-conditioned model can fulfill
the independent-validation objective. It does not authorize source downloads
or new outcome access.

This comparison read current parallel reports and existing source-method
descriptions only. No new primary-literature claims, expression/effect rows,
downloads, outcome access, model imports or jobs were used. Numerical agreement
and successful metadata admission remain separate from scientific validation.

## Currin public denominator: inventory-only follow-up

**Decision.** A public complete all-tested variant×peak denominator is **not
established** by the inspected official inventories. The source provides more
than significant lead variants: `caQTL_summaryStatistics_for_caPeaks_1kb_and_1mb_analyses.tar.gz`
is publicly listed, but explicitly restricted to peaks having a significant
caQTL in either stated analysis. Its described chromosome-separated summary
statistics do not establish inclusion of every tested nonsignificant peak.
Zenodo lists this archive as 4.8 GB. The independent lab inventory gives the
same restriction. Background variants overlapping non-caPeaks are a separate
BED asset; this description does not supply all measured association pairs or
their complete testing denominator. Peak coordinates and motif background
sets cannot be converted to measured zero effects.
[Official Zenodo inventory](https://zenodo.org/records/15025748),
[official lab inventory](https://mohlke.web.unc.edu/data/liver-chromatin-accessibility-qtl-from-138-individuals-currin-et-al-2025).

The primary publication points to these supplements and GEO GSE264684 for
sequencing data. Its methods specify variant scans within 1 kb and 1 Mb of
peak centers. That testing strategy is not evidence that the public archive
includes all scanned peaks. GEO returned CAPTCHA during this bounded check;
no hidden or alternative all-tested association inventory was established.
[Primary publication and data access](https://genome.cshlp.org/content/35/7/1485).

The local source contract states GRCh38. C2's producer specifically reads
`liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz`.
The broader local nominal consumer names chromosome files
`caQTL_hg38_chr{1..22}_with-varBeta-and-MAF_FastQTL-WASP-normal-nominal-updatedConsensusPeaks-25ATAC-PCs-minMAC10-results_1mbFromPeakCenters.txt.gz`
under `data/external/currin_2025_caqtl`. Its parser defines a row as a variant
identifier paired with a peak identifier and parses distance from peak center;
it therefore supports a variant×peak nominal schema in code. These files were
not opened. Their actual row counts, exact archive membership and completeness
across the publicly selected caPeaks were not independently verified here.
The existing peak axis remains the authoritative local GRCh38 BED0 interval
join; exact nominal-pair coordinate/allele correspondence needs a metadata-only
producer before any larger endpoint is admitted.
[C2 input producer](../../../scripts/analysis/alphagenome_program/c2_01_build_variants.py),
[nominal consumer schema](../../../scripts/analysis/alphagenome_program/e1_bridge.py),
[local source contract](../../../GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1/source_contract.tsv).

**Scope consequence [Inference].** A nominal-pair experiment within
significance-selected caPeaks could ask whether the model distinguishes
variants and supplied regions within that selected target family. It must
not be called all-peak or arbitrary-target validation. Current C2/pilot
endpoints remain selected lead associations. An admitted all-tested
all-peak roster would require a distinct source artifact or authorized source
reconstruction; neither has been established or performed by this review.
Only official article/inventory pages, filenames, local source metadata and
consumer code were inspected. No QTL effect file, counts, expression download,
archive payload or protected endpoint was opened.

### Completed pilot output-contract review

The compact receipt, sequence-admission JSON, summary and score-column header
for `target_readout_pilot_21995462` were inspected after completion. They
record the unchanged fixed selection hash, 64 selected cases, 61 sequence
admissions/scored cases, three non-ACGT exclusions and no replacements.
The receipt names the three required liver ATAC tracks, both readouts,
DIFF_LOG2_SUM, and an actual JAX L40S device. The reviewed producer requires
exact returned REF/ALT track identities, 1-bp full-window geometry and complete
unclipped peak containment before saving either readout. No material selection
or readout-contract defect was identified. Three fixed first-admitted repeats
report zero center and peak differences in this process; this is numerical
repeat evidence only, not prediction accuracy or general cross-process
reproducibility. The receipt distinguishes synchronized restoration time from
the shorter construction time called `load_seconds` in the summary; cost
reporting should use that distinction. No effect rows or accuracy estimates
were read or calculated for this check.

## Root follow-up: actual Claude advisory critique and scientific checks

On 2026-09-30 root invoked the installed Claude CLI with **all tools disabled**,
strict empty MCP configuration, safe mode and no session persistence. The
bounded one-turn critique completed successfully in 175.18 seconds. This is
a new advisory Claude reviewer, not the historical Model A owner; its reply
cannot amend that owner's likelihood/simulation contract.
[Exact supplied evidence](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/claude_parallel_review_request_20260930T2215.txt).
[Unedited response](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/claude_parallel_review_response_20260930T2215.json).

Its preference for RNA–chromatin and the allelic likelihood is an assessment,
not evidence that either is independently validated. The allele likelihood
needs calibrated inference **and** independent-participant biological
confirmation; simulation calibration alone does not satisfy the latter.
The sequence mask comparison continues as development characterization;
its spent development fold cannot validate a newly selected model. The
frozen-program missing-contribution route remains conditional on identifiable
full-RNA exposure units and a compatible receiving assay/target axis.

The useful Claude concern is an exact algebraic limitation of the RNA–chromatin
comparison. With error `e_ir=Y_ir-pred_ir`,
`mean_i,r(e_ir^2)=mean_r(mean_i(e_ir)^2)+mean_i,r((e_ir-mean_i(e_ir))^2)`.
A common assay-specific regional offset changes the first term and can alter
the absolute-profile primary contrast without improving prediction of
differences among people. Prespecify this decomposition as descriptive
interpretation beside the unchanged primary, with all 17 donors and fixed
baseline errors. This is not receiving-target calibration or a new fitted
prediction. Pairwise differences, if displayed, share donors; there are not
136 independent biological pairs.

Root implemented `profile_error_decomposition.py`; compute job **21998406**
passed hand-calculated three-person/two-region examples, common regional-offset
invariance of the across-person term, exact prediction, donor-order and region-
duplication invariants, and nonfinite/shape refusal. It completed in three
seconds, exit 0, empty stderr. Zero reported RSS is unusable. No receiving
arrays or bootstrap accuracy were evaluated. The 17-donor primary bootstrap
will be approximate and conditional on frozen model/target coefficients;
no new significance or outcome-selected pass rule is introduced here.
[Executed mathematical checks](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/profile_error_decomposition_21998406/summary.json).

The Claude response incorrectly stated that local-versus-global development
had no interval. Root independently inspected the existing permitted
`chromatin_conditional_21849192/results.json`: increment **0.022136409267305077**,
95% interval **[0.01734314188026808, 0.02661951532826109]**, Holm p
**0.0014992503748125937**, family of the three development contrasts. This
does not determine receiving MSE improvement, source-disjointness or power
at n=17. The advisory's generic claim that a percentile bootstrap at n=17
necessarily undercovers is not established by the supplied evidence; small-n
uncertainty is a real limitation, but no numerical coverage claim is made.
[Existing development evidence](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/chromatin_conditional_21849192/results.json).

The reference-only work now localizes the receiving RNA problem: standard
Ensembl 98 cDNA/ncRNA FASTAs cover **41,153/42,163** modeled genes. All 1,010
missing genes are the same-version **TEC** GTF entries with matching IDs,
1,037 transcripts and 1,049 exons; every transcript has exon geometry.
This permits investigating an explicit genome/exon reconstruction rather than
claiming measured zero or silently dropping genes. The original private
quantification recipe remains unknown. A native Salmon/tximport route with
such a complete reference would still be a qualified measurement approximation.
[Executed annotation disposition](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/ensembl98_missing_genes_21998407/summary.json).

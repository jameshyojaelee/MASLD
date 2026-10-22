# Source-owner requests for biological evaluation, 2026-10-05

Request1 is an unsent external draft; sending requires explicit user
authorization. Request2 is posted in the authorized local agent exchange;
the owner reply remains pending. No source has newly qualified as untouched external
predictor confirmation. Model A already specifies a closed biological
participant replication roster; its existence is preserved below. These requests resolve concrete source/identity questions;
they do not select a model, waive source restrictions or authorize new fits.
The Resource remains independent of experimental follow-up outcomes, and
the 117 program memberships/weights remain frozen.

## 1. GSE244832: authoritative paired-assay identity metadata

**Recipient:** Allen Wang, `a5wang@health.ucsd.edu`; corresponding author
verified in the official article front metadata. Tatiana Kisseleva
`tkisseleva@health.ucsd.edu` and David A. Brenner
`dbrenner@sbpdiscovery.org` are also listed corresponding authors.
Source: [Kim et al., producer article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12306656/),
DOI10.1016/j.jhep.2024.10.044;
[GSE244832](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244832).
Official PMC efetch front-only contact extraction SHA256
`ec164c96e21e462a446da0a4c2a8839377061db820b7ed70d69b87e0e0a7bf3b`.
No article Results or clinical/sample values were extracted.

**Subject:** GSE244832 RNA–ATAC participant/aliquot crosswalk

Dear Dr. Wang,

We are assessing whether the existing GSE244832 RNA and ATAC measurements
can support a participant-aware paired-assay comparison. The public metadata
list 18 JB RNA identities and 18 MM ATAC identities, but we have not found
an authoritative join. Could you share an identity-only table with:

- arbitrary study-local `person_key`, `RNA_JB_id`, `RNA_GSM_accession`,
  `ATAC_MM_id`, and `ATAC_GSM_accession`;
- study-local `specimen_key` and assay-specific `aliquot_key`, distinguishing
  one specimen split across assays from separate biopsies or repeat samples;
- pair status and reason for any unmatched assay, plus an all-samples roster
  so missing pairs cannot be silently omitted;
- reference assembly/annotation version and links to existing source-native
  RNA UMI/gene-axis and ATAC fragment/coordinate metadata, without matrices,
  reads or clinical values at this stage;
- an aggregate adult-only eligibility confirmation if applicable, without
  individual ages or other clinical fields.

No actual person identifiers or genotype data are requested. Please also
identify any shared participant source with GSE296875, GSE267145, GSE128072,
GSE281367 or the UNC/Currin GSE26105 source, using an aggregate overlap
attestation rather than identities. Unknown overlap is an acceptable answer.

The intended eventual endpoint is donor-level measured RNA-to-accessibility
prediction on the admitted native ATAC region/count definition. Assays and
aliquots are not independent participants. A source-native donor mean/simple
linear predictor would accompany any frozen candidate, with preprocessing
and comparison fixed before new targets are read. This is a request for
identity and measurement metadata only, not a request to perform an analysis.

Thank you.

**Limit:** GSE244832 already has extensive project development exposure.
Its crosswalk cannot make the existing source an untouched confirmation
cohort. A new paired endpoint might support prespecified development/transport
analysis only after a specific target, candidate and baseline are declared.
Frozen bulk RNA→H3 coefficients cannot be silently applied to snRNA→ATAC.
Source-disjointness is separate from checkpoint/pretraining exposure; root
must check both against the selected candidate's exposure registry.

**Existing evidence:**
`2026-10-01_codex_additional_paired_liver_source_comparison.md`, finite
GSE244832 crosswalk section; source SOFT SHA256
`43bfa6b65ecad054298b256cca3e5027a9d433d07feff44ed318a1bf078b2286`.
Article/supplement and generic ATACdemultiplex routes have already been
checked; metadata inventory did not supply the crosswalk. No pairing by
sample order, barcode prefix or disease label is permissible.

## 2. Claude/source owner: Model A design summaries and future-seal eligibility

**Recipient:** owning Claude agent through the existing local agent exchange.
This is a local review packet, not an external message. Binding source:
`2026-09-29_model_A_likelihood_spec_v1.md` §6 and accepted Appendix B3;
reduced-export interface:
`2026-09-30_codex_program_model_options.md`, “C5 full-information owner-side
helper: exact sufficient-input contract”, and superseding privacy restriction
`2026-09-30_codex_gateA.md:155–177`.

**Subject:** Model A owner-private C5 completion and sealed-roster eligibility

First provide an outcome-free identity qualification for the development F558
roster as well as the sealed F237 roster. The prescribed development counts
208/160/83/54/53 and sealed composition are historical roster assumptions,
not independently verified participant counts here. Current paper authorities
distinguish genetic clusters from verified people. The 774-cluster correction
candidate cannot be compared with F558/F237 by subtraction because their
eligibility definitions differ.

Export only aggregate counts and hashes establishing verified person/library/
aliquot grouping, unresolved joins and exclusions, cross-cohort duplicates and
kinship grouping, and reconciliation with the frozen development/sealed split.
State how the current identity/source-label exclusions apply, including any
GSE193066 inputs used for replication, attenuation, ICC or attainable-ceiling
claims. Bind C5 summaries to this qualified design version. If identity
exclusions change cohort sizes, centering or cluster resampling, explain which
information summaries need correction rather than reusing the historical
counts. Keep restricted identities, clinical values and sealed outcomes closed.

Please provide the actual owner-reviewed full-design C5 reduced summaries
for the eight existing null-plant scenario/threshold cases, or an explicit
missing-input/status packet. Return only adjusted information I, nonnegative
cohort precision contributions V_c, signed cohort slope-response R_c,
aggregate ranks/counts/check residuals and schema/input hashes. Include sum
checks and independently computed projection-agreement residuals. Do not clip
negative R_c or choose a lambda pooling convention implicitly. Declare the
intended P2 information weighting where the specification still needs that
choice. Full covariates, ancestry cross-products, gene/cohort information
matrices and efficient-score coefficients stay owner-private; no genotypes,
real allelic counts or outcome arrays should be exported to Codex.

Separately, provide an outcome-free readiness/identity attestation for the
already prespecified sealed F237 roster, comprising GSE174478 whole92 plus
fixed hash parts of development-source cohorts. The owner should verify
internally and export only aggregate summaries and hashes establishing:

- immutable person/technical-library/aliquot grouping, exact split manifest
  identity and absence of participant duplication across development/sealed
  parts; retain unknown or failed joins explicitly;
- source reuse/measurement exposure, excluded or previously displayed inputs,
  and the already specified FAM86JP sensitivity; shared-source hash parts
  are participant holdouts, not an entirely independent source;
- whether required metadata-only eligibility checks and owner-private technical
  input preparation are complete, with counts/versions and no actual sealed
  depths, stages, allelic values, clinical fields or restricted identities;
- exact remaining prerequisites for opening the seal, including numerical
  equivalence and accepted B3 statistical calibration, rather than inferring
  readiness from an observed fixture fit or null-plant information.

The existing §6 target is RNA-allelic replication of stage-associated
proportional allelic-effect modulation, using the prescribed conditional
beta-binomial genotype-mixture model and its nuisance/reference-bias controls.
Its comparator/null and calibration are already specified; this packet
neither proposes an ATAC/H3 surrogate nor changes P1/P2/P3 rules. §6 specifies
refitting all parameters on sealed F, the same-sign development kappa plus
P1 p<0.05 replication rule, and the fixed P3 projection comparison. Please
confirm source-definition/eligibility only without opening the seal or
producing a sealed statistic. Any further private work requires root's
separate execution decision.

**Development-only allocation planning:** Please also return owner-approved
upper bounds for the development model's retained gene count G, retained
design width Z (including declared dropped columns/missingness indicators),
eligible rows per gene, maximum supported count-axis depth, and per-gene
sum(n+1). Bounds may be conservative protocol limits rather than observed
summaries. No sealed depths, real count vectors, covariate values or genotype
information are requested. The current synthetic helpers cap rows at256 and
per-gene count-axis entries at100,000; these are implementation limits, not
scientific inclusion rules. The P770 allocation test reuses24 synthetic rows
and cannot establish actual development feasibility or justify dropping rows.
Retained G and Z must come from the owner-reviewed eligible design, not a
genetic-membership count.

**Limits:** The biological replication route already exists and remains
closed. Numerical/statistical completion and identity qualification are
required before opening; this is not a missing-new-cohort claim. The literal
F237/GSE17447892 totals are coordinator-supplied roster metadata, not counts
recomputed from protected files here. Spec §6 directly supports the roster
composition. Null-plant I/V/R cannot establish MDE at non-null plants, CRT
size/power, P2 coverage, a theta decision or sealed replication. Reduced
summaries require owner review and are not a formal privacy guarantee.

**Signed-splice alternative, not a third request:** The previous bounded UNC
review found a real RNA splicing endpoint but unsigned gene-level tests;
GSE179250/GSE26105 participant ancestry overlaps the Currin deposit at the
source-ID level. Existing signed variant-event slopes, all-tested denominator
and disjoint participants remain unverified. No request for new association
fits or protected Wenz source data is prepared. This route cannot currently
substitute for the existing Model A sealed participant replication path.

## Disposition

### Additional request to the Model A specification owner: calibration settings

Before constructing B3, please identify the controlling frozen settings for:
the warp-speed resampling construction and allocation across scenarios; which
interval and endpoint receive the coverage band; simulated label-error and
lambda scenarios; reference-bias slope amplitude; non-beta-binomial tail
distribution and weight; tag-switching amplitude; and omitted depth-dependent
bias function and magnitude. These are scientific specification questions,
not requests for protected counts or sealed outcomes.

The accepted likelihood Appendix B3 fixes at least4000 warp-speed datasets
and200 full CRT datasets, with the one-sided lower95% size bound criterion.
The original final plan lines168–176 gives size[.04,.06], coverage[.93,.97]
and kappa{-0.10,-0.05,0,+0.05} by tau{0,tau_max/2,tau_max}. Please confirm
their application under the revised B3 procedure. A90% P2 interval must not
silently be evaluated against a95% coverage target. Diagnostic nuisance values
and agreement-fixture lambda draws do not supply the missing plant settings.

Numerical status relevant to this request: adaptive-order draw0 proof22083807
completed and matched the older GL fit, but its unchanged-family repeat22084113
retained null nonconvergence despite LL agreement within1.4e-11. The50-draw
family remains a diagnostic, with no calibration/adoption claim. Base count
generator22084262 passed synthetic law/selection checks and indexed-seed replay;
it supplies none of the missing B3 scenario settings. No calibration workload
or changed scientific tolerance is requested by this note.

Subsequent separate saved-point polish checks22084403/22084404 repaired
repeat0 and draw11 on Intel/AMD under the original1e-6 convergence cutoff.
The empirical score-integral gain checks are not certified error bounds.
Three separate full-fit trials22084419/22084420/22084421 are now running;
full-fit reproducibility is unverified and existing failures remain unchanged.
V1.7 item8 permits Newton-polish settings without revising the warm tie rule.
This supplies none of the missing B3 scenarios, plants or identity inputs.

Further numerical result: full trial22084420 retained the repaired null but
failed a positive-variance integral at call14; the prescribed P3 failure and
exceedance remain. The other two full trials are live. G5 now passes15/15
fixed synthetic likelihood checks, which is neither model calibration nor
human replication. Fixed-point variance-score precision and one
count-prefix integral remain under investigation; owner B3 settings and
qualified development/replication inputs are still required.

Subsequent results: scalar heat-identity variance-score checks passed on both
saved genes with two arithmetic backends. All three first-wrapper full trials
are terminal: Intel0 and draw11 selected fits converge, AMD0 still fails.
Draw11's null warm fit remained unconverged and cold was selected; a separate
saved warm-point replay now retains the original constrained solver and reaches
the1e-6 cutoff with one additional score-integral step. This is numerical
progress, not B3 calibration, fifty-draw agreement or qualified human evidence.

### Additional specification-owner question: converged cold start within the warm tie

The prescribed draw11 now has a completed separate GL diagnostic22082995.
All25 positive-variance GL evaluations returned, with no integration refusal.
The final null and marginal states retain the original warm selection and fail
the strict projected-score criterion:2.1271003368761657e-6 versus1e-6.
The saved cold starts are marked converged, but their log likelihoods are only
3.944933268940076e-11 (null) and1.3642420526593924e-11 (marginal) below warm.
The unchanged rule chooses cold only if its LL exceeds warm by at least1e-6;
the boundary/null snapshot then preserves null nonconvergence. Native Gamma
and GL return exactly the same null, final parameters, LL and gradient here.
Root evidence: REC
`repair_p3_gamma_GL64_draw11_root_review_20261007T073350Z/review.json`.

Please confirm that this remains a failed supplied draw under B-FIX and B3.
If a convergence-aware selection within the1e-6 tie is intended for a future
model revision, provide that revised frozen rule and its required comparison
and calibration family before implementation. Current sources, thresholds,
failure/exceedance accounting and completed results remain unchanged. This
question does not authorize choosing a different returned start, ignoring a
failure, or treating the completed diagnostic as calibration.

Request1 could resolve a real pairing blocker in exposed development.
Request2 completes prerequisite evidence for the existing closed Model A
participant replication route. Neither draft supplies untouched external
predictor-effect labels. If identity or calibration requirements fail, retain
that limitation rather than opening the seal or selecting a rescue hypothesis
from tested outcomes. No external contacts, bulk acquisition, fits, jobs,
dependencies, authority changes or protected data reads occurred.

Latest synthetic update,2026-10-07: the500-row prefix integral is complete
and agrees with the saved scalar LL and full Gamma arrays. Gamma-only null
draw11 still fails the original strict convergence cutoff; its unchanged
tie selects the nonconverged warm start. Both saved variance coordinates
now pass independent80-digit probability-only checks, and a fixed-normal
Hermite candidate passes the saved full eight-gene comparison. Separate
full-fit pilots are running. These numerical checks supply none of the
requested B3 scenario/plant/coverage, C5 or qualified participant evidence.
See the existing calibration requirements and Codex-to-Claude exchange for
the saved comparisons and failure accounting. No protected outcomes opened.

Update2026-10-07 11:05UTC: all three fixed-normal Hermite full pilots failed
positive-variance order agreement. A separate saved-null single-start diagnostic
reproduces refusal at variance0.9234818606206378 in all eight genes. Exact
parallel null-shift scheduling passes two byte-exact endpoint checks and is
retained as a numerical speed improvement; it supplies no scientific calibration.
A centered importance-quadrature revision has four independent fixed-point
jobs running. No new full fits, protected reads or human validation follow.
All four centered-importance checks subsequently failed; no full fitting is
qualified. Adaptive-interval reference22085179 is running at the saved refusal
point, with no fits or human data. It subsequently passed all eight genes;
this remains a fixed-point reference, not a qualified full-fit method.
The quantitative B3, C5 and qualified human-input requests above remain open;
numerical progress does not supply any of those missing inputs.

Update2026-10-07 11:20UTC: probability-weighted adaptive conditioning now
passes all four fixed-point checks, including the full8 GL and both independent
Decimal references. Three separate full-fit pilots22085227/28/29 are running,
retaining original starts, selection, convergence and failure accounting. This
is numerical development only. It supplies none of the quantitative B3/C5 or
qualified participant/design evidence requested above. Protected outcomes remain
closed and no model/paper claim is adopted.

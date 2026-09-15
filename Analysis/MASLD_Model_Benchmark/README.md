# MASLD Model Benchmark

The authoritative human-readable campaign plan is
[`OVERALL_PLAN.md`](OVERALL_PLAN.md). The executable family/task and
cross-cohort requirements are
`config/evaluation/family_native_tournament.toml` and
`config/evaluation/cross_cohort_expansion.toml`. This README records the package
implementation and operational evidence beneath that plan.

The primary objective is to build the strongest useful, publishable, and
externally transferable multimodal MASLD research model. This package is the
search and proof mechanism: it freezes data topology, missingness, the broad
model census, checkpoint exposure, splits, runtimes, metrics, adoption criteria,
selection, held-back evaluation, and release provenance. It does not modify the
current MASLD Resource, Gene Catalog, portal, figures, or the 117 frozen
programs.

The implemented code is currently research infrastructure, not a trained
scientific model. Five HVG/PCA cell-state baselines have executable scientific
adapters; the generic stub is a deterministic requirement fixture only. Other
scientific adapters, model environments, dataset joins, and unresolved
checkpoint audits must pass inclusion before their smoke or training waves can
run. No public API, portal, or clinical product is in scope.

## Family-native multimodal tournament

RNA-conditioned ATAC is one high-powered specialist lane, not the universal
objective and not the eligibility filter for the broader model census. The
tournament evaluates each architecture only on biological tasks supported by
its native inputs, output semantics, reference build, and held-query behavior.
It does not collapse unlike outputs into one leaderboard or reward a model for
being inapplicable to fewer assays.

The active task families are:

| Task family | Held-query inputs | Native outputs | Principal model families | Required simple comparators |
|---|---|---|---|---|
| Sequence-native regulatory prediction | DNA sequence only | Assay-matched coverage, accessibility, regulatory class, splice or expression tracks | Borzoi, Enformer, AlphaGenome, Sei, ChromBPNet, BPNet, scBasset, sequence CNN/transformer controls, and DNA language encoders with identical fitted heads | Mean track, GC/k-mer controls, deltaSVM/gkm-SVM, nearest regulatory feature, and matched from-scratch CNN/transformer |
| Signed variant and variant-to-regulation transfer | REF and ALT sequence, with context only where native | Native allele delta, MPRA activity/direction, accessibility delta, signed RNA effect, or enhancer-gene retrieval kept as distinct endpoints | Corgi, Borzoi, Enformer, AlphaGenome, Sei, DNA language encoders, released scooby checkpoints, ChromBPNet/BPNet, MPRALegNet, and link models on their supported endpoints | Zero allele delta, mean track, nearest gene/feature, deltaSVM/gkm-SVM, ABC/rE2G for retrieval, and assay-native linear heads |
| RNA-conditioned regulatory prediction | DNA sequence plus observed RNA or a training-frozen RNA context; no held-query ATAC | Donor-by-cell-state ATAC profile/count and supported regulatory tracks | Regular Corgi, any executable Corgi+ release, RNA-to-ATAC translators, and the conditional context-Borzoi or sequence-plus-RNA FiLM/LoRA model | Sequence-only, RNA-only, shuffled context, nearest context, train-donor cell-state mean, shrunken pseudobulk, ChromBPNet/BPNet, and matched CNN/transformer |
| Observed-multiome prediction | DNA plus explicitly observed held-query ATAC and/or RNA according to the model's native topology | Masked accessibility, functional genomic tracks, paired-modality reconstruction, or native joint representation | EpiBERT, EPCOTv2, GET, released scooby contexts, MultiVI and other true paired-multiome methods | Observed-ATAC-only generalized linear models, sequence-only backbone, RNA-only context, WNN/LSI/PeakVI, assay-native pseudobulk, and shuffled/masked-modality controls |
| RNA cell-state and disease-state transfer | Single-cell, single-nucleus, or donor-level RNA as explicitly registered | Frozen cell identity, continuous donor-level state, or bulk cross-sectional stage association | Cell foundation encoders, scVI/scANVI, transcript-native encoders, and calibrated task heads | HVG/PCA, nearest centroid, kNN, logistic/elastic-net/linear-SVM, CellTypist, and assay-native donor pseudobulk |
| Perturbation and evidence integration | Perturbation identity/sequence or typed evidence graph | Replicate-level response or held-source candidate recovery | Perturbation-native models and typed graph families | Control/perturbed means, ridge/bilinear/Systema, evidence counts, degree, nearest gene, and regularized linear evidence |

Donor phenotype, histology, and clinical metadata form an additional masked
context lane under the binding rules in `OVERALL_PLAN.md`. Fibrosis, NAS and
its components, age, recorded sex, BMI, metabolic status, etiology, and other
source-verified fields may be endpoints, covariates, context, or transport
strata, but never an input to their own prediction task. All harmonization and
imputation remain training-fold-only, and molecular-only, metadata-only, and
shuffled-metadata controls are mandatory.

Observed-ATAC models are not excluded. They enter a separately registered
`observed_multiome` TaskSpec whose inference request explicitly permits the
observed modality. Those results cannot be relabeled as RNA-conditioned ATAC,
cannot compete on an RNA-only held-back endpoint, and cannot donate query-derived
neighbors, normalization, calibration, or embeddings to another task.
GSE296875 supports same-nucleus objectives; GSE244832 remains
same-donor/different-aliquot and is never converted into false cell pairs.

The following named families must receive an inclusion attempt and a native-task
disposition in the frozen census: regular Corgi, Corgi+, Borzoi, Enformer,
AlphaGenome, Sei, DNABERT-2, Nucleotide Transformer, HyenaDNA, Caduceus and
eligible Evo checkpoints, EpiBERT, EPCOTv2, all released scooby checkpoints,
and the context-Borzoi/sequence-plus-RNA FiLM-LoRA candidate. A missing loader,
unresolved checkpoint identity, incompatible output, contamination, or terms
restriction is recorded as a terminal model-task result rather than used to
silently omit the family. Restricted models may remain scientific comparators
but cannot become the openly released best model.

The context-Borzoi/FiLM-LoRA candidate is prespecified but does not receive a
full training campaign merely to improve the paper narrative. A small
outcome-separated compatibility and ablation probe may run during the family
sweep. Full training is activated only when held-development residuals show
complementary sequence and RNA/context information under the fixed trigger;
failure of that trigger receives an explicit `not_justified` disposition.

Selection is task-specific. Every task retains its strongest simple method,
uses donor-safe and chromosome/LD-safe splits where applicable, and reports
checkpoint exposure and terms separately from score. No architecture earns a
universal or multimodal-superiority claim without noninferiority to every
eligible specialist, superiority on at least two independent source families
after multiplicity correction, and an appropriate genuinely external
evaluation. A baseline or a set of different per-task best models is an accepted
final result.

## Cross-cohort and cross-modality development requirements

GSE296875 is the highest-powered same-nucleus RNA-ATAC development anchor, not
the development population and not a MASLD phenotype cohort. Its 39 donors can
establish whether a model learns cis sequence and trans RNA context under exact
pairing, but cannot by themselves establish MASLD specificity or cross-cohort
transfer. No task may select a best model from a random pooled-cell split or from
GSE296875 alone.

Development is organized as assay-native cohort families. The existing
seven-source, 102-donor single-cell Atlas supplies RNA cell states; GSE244832
supplies donor-matched but not cell-matched RNA and ATAC; GSE281367 supplies an
independent MASH ATAC cohort; the current bulk source pool supplies donor-level
RNA across studies; and PXD051911 supplies protein transportability without
inventing an RNA-protein pair. GSE256398 is source- and reference-audited
external development: 197,942 deposited CellBender-filtered, pre-project-QC
barcodes from 26 donors, with 35,455 of 36,601 features matched by
version-stripped stable ID to GENCODE v49 and the other 1,146 explicitly
masked. Full activation remains blocked because S35's automatic Scrublet
threshold called 8,762 of 11,080 barcodes doublets. The source-exact
reconstruction retains 165,372 nuclei; an outcome-blind capped sensitivity
retains 172,997 and requires leave-S35-out analysis. No barcode-level author
annotation or author-final singlet membership was deposited, so reconstructed
or transferred labels cannot become supervised truth. The 9
alcohol-associated donors remain an etiology challenge rather than
MASLD-negative controls. GSE289173 remains project-reserved. The
conditional FNIH cohort does not delay development.

The frozen local-source audit
`atac-transport-source-audit-21064435` hashed 47.0 GB across 34 inputs and
recovered 226,224 GSE281367 cells from 12 donors plus 88,814 GSE244832 cells
from 18 donors. It also rejected every legacy lineage pseudobulk as a benchmark
outcome. GSE281367 retains 18,729–50,622 duplicate peak instances by lineage;
GSE244832's deduplicated matrices no longer match their stored library-size
metadata. The corrected raw-fragment build is frozen as
`atac-transport-outcomes-21064552` (ARTIFACTS SHA-256
`6f5092f1343b91ff592851d7d265ab45f8e154d210984c19f96bf05421025be3`).
It read all 30 fragment files to EOF and produced evaluator-only 20-bin profiles
and counts on 32,000 project windows for 96 eligible donor-lineage units. The
two pipelines have different assay-native Tn5 coordinate requirements: the custom
GSE244832 fragment builder requires start+4/end-5, while GSE281367 Cell Ranger
ATAC fragments are already adjusted and require start/end-1. This activates
development transport outcomes only; model inputs cannot read them.

The activation queue broadens both cohorts and modalities. Highest priority is
the GSE267145 SuperSeries: its source paper reports 108 biopsy participants,
GSE269412 contains 262 liver RNA-seq records, and GSE267119 contains 99 liver
H3K27ac CUT&RUN records. The first frozen official-SOFT audit found 99 unique
one-to-one title-prefix matches and no H3K27ac-only record. The subsequent
authoritative join (`gse267145-authoritative-join-21064930`, ARTIFACTS SHA-256
`e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580`)
triangulated official GEO, the source article, deposited titles, source-native
histology, and both matrix axes. It freezes those 99 paired participants as
same-sample different-aliquot with stage-balanced folds. The paired matrix audit confirmed 43,285 RNA features across all 262
RNA participants and 96,460 H3K27ac regions across the 99 paired participants,
with complete sample-axis recovery. The H3K27ac matrix is nonnegative integer
count data. RNA-measurement audit `gse267145-rna-measurement-21065097`
(ARTIFACTS SHA-256
`38531e463837cafd2c8e321411048c7b8b493876c3c964c0f50790cac8806f2f`)
scanned all 11,340,670 entries and the exact last public PISCES tag before the
deposit. It classifies the RNA as nonnegative fractional expression estimates.
Continuous transforms fit inside each outer training fold are permitted;
integer count likelihoods, rounding, and DESeq2 reconstruction without the
missing tximport offsets are prohibited. Reference audit
`gse267145-reference-crosswalk-21065154` (ARTIFACTS SHA-256
`da491302ea738a92017978518239bc3884a7635d580e01e3bd28cd4f55b630c2`)
then mapped 42,163 genes by stable ID to GENCODE v49 and explicitly masked
1,122 retired genes. All 96,460 H3K27ac regions lie within GRCh38.p14 primary
contigs, but both 0-based half-open and 1-based inclusive interpretations pass
bounds. QC/rights audit `gse267145-qc-rights-21065336` (ARTIFACTS SHA-256
`ac7c92fd01b2b4240b380e4f6a58d98f1a6f853f8d0ad98d7167d7af34baaa73`)
confirmed that all 99 paired participants pass the hard finite, nonnegative,
positive-library rule. Distributional diagnostics do not automatically exclude
participants. Public internal nonclinical research use needs no new DUA, but
data redistribution is prohibited and derivative-weight release requires
model-specific review. Sequence extraction remains prohibited until that
coordinate-base ambiguity is resolved. GSE105127 adds 19 participants with three microdissected liver zones and matched
RNA/RRBS libraries. Authoritative join/rights audit
`gse105127-authoritative-join-rights-21065552` (ARTIFACTS SHA-256
`574eeee1119054a2d7b55066804accf29fb45f23b591dd26f9d61d907b27a75f`)
verified all 114 records using official GEO and source-article evidence, froze
five participant-safe folds, and established that RNA and RRBS use adjacent
cryosections. They are not same-section observations. Public internal
nonclinical research use needs no new DUA, but redistribution and open
derivative weights remain restricted pending review. Label-free activation
`gse105127-assay-native-activation-21082444` (ARTIFACTS SHA-256
`fc529784b5f6df53af7ba82036ff617507ca9d551d64d79aac6cb9061faf62bb`)
binds all 57 RRBS BEDs and 800,713,381 cytosine rows, freezes percent/coverage/
strand semantics without using auxiliary columns as counts, and registers 57
single-end 76-bp RNA runs. Training remains blocked until the exact source
FASTA and chains, failure-aware round-trip CpG remap, raw RNA quantification,
evaluator, and fold authority are frozen. Bulk and OOD
candidates include GSE268273,
GSE260666, GSE274114, the fingerprint-confirmed shared GSE106737/GSE83452
Antwerp–Inserm family, and GSE49541. GSE202379 is already a source in
the current Atlas and is never counted or split as a new independent cohort.

GSE260666 is now an active external-development fixture with 16 unique
participants: six controls, six NAFL, and four NASH. Its label-free raw-count
model input contains 23,576 mapped genes, including 20,410 shared with the
GSE267145 training fixture; labels remain evaluator-only. The frozen audit
found no accession, BioProject, BioSample, run, or participant alias overlap
with GSE267145 and no molecular near duplicate. It is project-exposed rather
than held-back, and fitting remains blocked until the GSE267145 RNA finalist and
preprocessing selection record are frozen.

GSE268273 is now an active, project-exposed external-development fixture for
fibrosis and etiology transfer. It contains 109 unique biopsy participants: 69
with IMID-associated MASLD and 40 with classic MASLD. All groups are MASLD; the
IMID stratum is an OOD/comorbidity stratum, never a negative control. Exact
fibrosis, recorded sex, and metabolic fields remain evaluator-only. The 824
ENA runs are technical partitions of 109 participant libraries and must be
aggregated before inference. Exact file topology resolves 680 effective
single-end runs from 73 participants and 144 paired-end runs from 36
participants despite the deposited all-SINGLE declaration; paired mates must
never be concatenated as single-end. Raw plan
`model-data-075-21082574` (ARTIFACTS SHA-256
`5cbeaee318cab309b5ddc364e9e8966561888009b194941bffeaceda3601f553`)
freezes eight participant-indivisible bundles and records the actual RSEM
binary as v1.3.1 plus STAR 2.7.10b. The candidate axis contains 14,078
GENCODE v49 targets, including 11 deterministic raw-count sums. Scoring stays
blocked until the raw runs are requantified and the GSE267145-trained
preprocessing selection record is applied. The deposited cohort-wide voom values are
contamination-audit inputs only.

GSE274114 is now a label-free, project-exposed 39-participant v49 fixture
(`model-data-080-21082249-gse274114`, ARTIFACTS SHA-256
`b334ff4fe77903be82bb228fc479a28d92d74f4cec0b06ea98c3d6a6100c9300`).
It contains 9 healthy controls, 11 HBV-only, 10 MASH-only, and 9 combined
HBV-MASH participants across 71 technical runs. Healthy and HBV-only samples
use HiSeq 4000, while MASH-only and combined samples use NovaSeq 6000.
Therefore only healthy-versus-HBV within HiSeq and MASH-versus-combined within
NovaSeq are registered. Four-class accuracy, MASH-versus-non-MASH,
MASH-versus-healthy, and naive cross-platform OOD accuracy are prohibited;
HBV-only is never a MASLD-negative control.

The microarray activation audit now binds GSE106737 and GSE83452 as one
project-exposed cohort family using 78 reciprocal array fingerprints and 41
participant aliases. Thirty-eight paired visit tokens match exactly and 40
retain the deposited `followup`/`follow-up` punctuation difference; no
semantic mismatch or generalized token coercion is allowed. All 303 raw CEL
headers, platform axes, rights checks, and the tie-safe held-array isolation
fixture pass. Two existing R runtimes read both Calvin header fixtures, but no
audited runtime supports included GPL16686 single-array transcript-cluster
summarization, and the authoritative GPL16686 mapping asset is login-blocked.
The series matrices remain contamination-audit inputs only, and biological
normalization, model fitting, and scoring remain inactive.

Every eligible task uses leave-one-study-out outer evaluation when at least
three compatible cohort families exist. With two compatible families, one is
held out for explicit transport while all tuning stays in the other; with one,
the result is internal development only. Cohort families, participant aliases,
repeated biopsies, and all modalities from one person stay in one outer fold.
Training uses cohort-balanced sampling, and reporting macro-averages cohort and
donor effects instead of allowing a large study or many cells to dominate.

Modalities are not pooled merely because they share donors. RNA, ATAC,
H3K27ac, RRBS methylation, bulk RNA, protein, and spatial measurements keep
assay-native heads, likelihoods, QC, and endpoints. Same-nucleus, same-sample
different-aliquot, same-donor, and unpaired relationships remain explicit.
Missing assays carry masks and missingness states; they are never encoded as
zero. Histone and methylation heads can be activated only after their new
cohort joins, rights, builds, preprocessing, and independent-donor counts pass
inclusion. Until then, the conditional MASLD model remains RNA/ATAC-only.

The machine-readable expansion and source dispositions are frozen in
`config/evaluation/cross_cohort_expansion.toml`. Candidate discovery does not
make a dataset active: each source still needs its own `DatasetManifest`,
cohort-family deduplication, rights and reference audit, read-only output file,
TaskSpec binding, and compute-node fixture before model access.

The frozen census contains exactly 135 registry model IDs mapped to 103 audit
bundle directories by `config/evaluation/model_audit_bindings.toml`. Every
bundle must contain `checkpoints.json`, `development_crosswalk.json`, and
`exposure_audit.json`; grouped mappings are explicit, and census membership
alone never makes a model runnable or eligible to be the best model.

The first foundation-model lane is Geneformer V1 10M plus V2 104M and 316M.
Their pinned weight, config, tokenizer, median, Ensembl-map, and token-map
identities are frozen. The offline `checkpoint stage-geneformer` command
verifies every source object before copying, rejects symlinks and checksum
drift, denies all pickle global/class construction, and publishes sanitized
plain-dictionary auxiliaries in a read-only bundle. It performs no download
and never executes checkpoint content. Output file acquisition, environment
installation, and GPU probes remain approval-blocked.
The unsubmitted `slurm/acquire_geneformer_checkpoints.sbatch` is the only
registered acquisition path: one non-array `io` job obtains the 16 exact
objects at the read-only upstream revision, verifies all SHA-256 values, and
stages the three bundles. It must not run without explicit review of its job
hash and resource request.

The frozen smoke-view corpus crosswalk is mixed rather than clean. GSE136103
and the MacParland Liver Atlas source are declared in Genecorpus-30M and
Genecorpus-104M; GSE185477 is also declared in Genecorpus-104M. Results must be
stratified by these `encoder_seen` sources. A source not detected in the
declared tables remains `unknown` if it was already public before the model's
training cutoff; a documented first-public date after that cutoff supports
`clean_declared`. The separate date- and source-table audit for held-back
GSE289173 remains `target_label_unexposed`.

The released Geneformer checkpoints are included only to cell-state mapping.
Their tokenizer was designed for raw single-cell or single-nucleus counts;
bulk RNA is not passed through it. Any later bulk transport analysis must use a
separately frozen aggregation or program-projection requirements and a distinct
model record.

That modality check applies to every released cell-foundation checkpoint in this
family. A single-cell or single-nucleus tokenizer is not run directly on bulk
RNA. Bulk transport, if activated, must use an explicit donor-level aggregation
or program-projection model with its own record, folds, and assay-native
baselines.

The scGPT whole-human and continual bundles now have exact official Google
Drive file IDs, byte sizes, auxiliary hashes, architecture, and input requirements.
The weights themselves have not been downloaded, so both weight SHA-256 values
and their terms are still missing before inclusion. The whole-human checkpoint used
normal human cells from CELLxGENE Census 2023-05-08. GSE289173 was released
after both scGPT output files and is `target_label_unexposed`. Pre-cutoff
development sources remain `unknown` until membership in that retired weekly
Census build can be rederived; the adjacent 2023-05-15 LTS release is not used
as a substitute. The continual checkpoint adds supervised Tabula Sapiens
cell-type training. That is reported as a distinct externally supervised model,
and any Tabula Sapiens evaluation would be `continual_seen`.

scGPT is also included only to cell-state mapping. Its official embedding path
randomly samples expressed genes when a cell exceeds 1,199 in-vocabulary genes.
The common lane therefore requires a deterministic seed derived from each
read-only row ID and tests cell-order invariance; the official native behavior
is retained only as a separately labeled native run.

UCE 4-layer and 33-layer now have a frozen architecture, input requirements,
complete auxiliary roster, source MD5 values, code and weight terms, and a
complete training-corpus crosswalk. The 33-layer output file is correctly bound to
Figshare v5; it was added after v4 and cannot be cited as a v4 file. GSE185477
and the MacParland Liver Atlas source are declared training inputs and are
`encoder_seen`; held-back GSE289173 is `target_label_unexposed`. UCE is included
only to cell-state mapping from single-cell or single-nucleus counts. Its
official count-weighted gene sampling and chromosome shuffling are stochastic,
so the common lane again requires unchanging row-ID seeds and cell-order
invariance. Weight SHA-256 values, safe deserialization, archive inspection,
and the executable runtime remain blocked until approved acquisition.

The restricted cell-foundation census now has checkpoint, input, terms, runtime,
and exposure bundles for scFoundation, CellFM, and GeneCompass. The output file
directory `scfoundation` maps deliberately to registry ID `scfoundation_heads`;
its mutable SharePoint folder still lacks a read-only base and head inventory,
and its custom model terms prevent a released or conditional-model best model.
CellFM pins the exact 9,659,358,467-byte `base_weight.ckpt` and SHA-256, but
the benchmark allows no adaptation under CC-BY-NC-ND-4.0 and cannot share
adapted material; even a separable probe over frozen outputs requires counsel
review.
GeneCompass has no detected code, weight, or prior-payload terms, and its
ScienceDB `Base` label and required external prior objects remain unhashed, so
it is no-use until written permission and read-only objects exist. Its reported
53,568,337 human plus 48,200,083 mouse cells sum to 101,768,420, which remains
unreconciled with the scCompass-126M label. For all three models, every Atlas
source, including GSE289173, remains exposure-`unknown`: a GEO public date after
a paper or checkpoint upload does not exclude earlier private data or labels.
None is included to bulk transfer, the perturbation track, conditional-model
training, held-back best-model selection, or open release.

scLong is a restricted cell-state comparator only. Its Zenodo software record
declares CC-BY-4.0 metadata, but the GitHub root, separately hosted checkpoint,
outputs, embeddings, and derivative terms remain unresolved. The SharePoint
payload returned HTTP 403, so the exact weight and auxiliary bytes, hashes,
cutoff, and 48-million-cell corpus crosswalk are not bound. Exposure therefore
remains `unknown`; scLong cannot enter held-back selection, conditional-model
training, or open release.

CellPLM's census now distinguishes two different "85M" objects. The primary
entry is the repository-documented September 2023 checkpoint with the paper's
16-component GMVAE; the October 2023 object is an undocumented plain-VAE
variant and remains a separately named restricted comparator. "85M" is a model
size label: the paper reports 82,402,543 parameters and 11.4 million training
cells. Both objects predate GSE289173 and are target-label-unexposed, but both
remain blocked on weight terms and SHA-256. Common-lane embeddings are the
1,024-dimensional normalized encoder output before the latent, built only in
deterministic donor-local contexts. Native 512-dimensional GMVAE and VAE paths
stay separate. Oracle Leiden-resolution selection, train-plus-test HVGs,
cross-donor context, and ambiguous double normalization are prohibited.

TranscriptFormer TF-Sapiens now has a hash-verified versioned 1,819,975,447-byte
S3 archive, safely materialized checkpoint and vocabularies, pinned MIT source,
and an isolated executable runtime. The source package reports 0.6.1; the
`v0.6.0` tag points to a different commit and is not treated as an alias. A
1,000-cell, 102-donor outcome-blind L40S activation passed strict checkpoint
loading and deterministic 2,048-dimensional native and common embeddings. The
official acquisition CLI remains unused because its tar extraction does not
validate member paths or links.

TranscriptFormer remains `unknown` for held-back GSE289173. That GEO source became
public 57 days before the read-only checkpoint archive was created. The
official 644-dataset training table link was resolved but returned 403/429 in
this audit. GSE289173 is absent from the current CELLxGENE catalog, but a mutable
current catalog is not accepted as proof about historical pretraining. The
official inference path also preserves input feature order and truncates the
autoregressive context, whereas training randomized expressed-gene order.
Native source-order and deterministic count/Ensembl common-lane policies are
frozen separately. This allows a development comparator only; no task head has
been fit during activation, and unknown exposure still prohibits held-back or
best-model claim. In the subsequent 1,000-cell, 102-donor development smoke,
the finalized nested common MLP head reached donor-balanced macro-F1 0.9445
versus 0.9229 for the strongest registered classical baseline (HVG/PCA/
logistic): delta 0.0216, paired donor-bootstrap 95% CI 0.0071 to 0.0368, with a
lower Brier score. It remained below the existing Geneformer-316M MLP result
(0.9505). These are retrospective smoke comparisons, not external evidence;
earlier foundation prediction bundles are not asserted to implement the now-
finalized nested common-head recipe.

SCimilarity v1.1 now has a frozen architecture, preprocessing requirements, code
and weight terms, exact 30,310,810,843-byte Zenodo archive metadata, and a
complete training/test/reference crosswalk. GSE185477 is an annotated training
source (`encoder_seen`), while GSE136103 occurs only in the released search
reference (`reference_only`). Common-lane probes use only newly computed
128-dimensional embeddings and never the released reference index. GSE289173
postdates the finalized v1.1 record and is `target_label_unexposed`. The archive
SHA-256, internal gene order and layer-size files, safe state-dict staging, and
runtime remain blocked until approved acquisition. The paper reports 56
training studies while its final Supplementary Table 1 contains 52 nonempty
`train` rows; both counts are retained as a source discrepancy.

scPRINT v1.5 medium is pinned to the original March 2025
`v2-medium.ckpt`, not the mutable current filename. The original object was
renamed byte-for-byte in December 2025, then silently replaced by a different
object in May 2026 after GSE289173 was public. That replacement is ineligible.
The paper code is MIT, the Hugging Face weights are Apache-2.0, and the declared
pretraining source is CELLxGENE Census LTS 2023-12-15. This makes the held-back
source `target_label_unexposed` for the original checkpoint. Pre-cutoff
development sources remain `unknown` unless the unpublished post-QC
train/holdout inventory resolves them. Acquisition, restricted checkpoint
inspection, vocabulary and ontology freezing, deterministic row-ID sampling,
and the Torch 2.2 runtime remain approval-blocked.

Common-lane cell runs are separately identified as `linear` or
`two_layer_mlp`. Both use the frozen donor-grouped five-fold head ensemble,
donor-and-class-balanced loss, training-only embedding standardization,
prespecified early stopping, and out-of-fold temperature calibration. Native
model recipes and mandatory baselines remain separate runs.

scPRINT-2 small-v2 is pinned to its original December 2025 object and stable
1.0.3 code release. The later `small-v2.ckpt` name is a byte-identical rename.
Its declared corpus combines CELLxGENE, Tahoe, and a continually expanding
scBaseCount collection, but no read-only source snapshots or retained/holdout
manifest are bound to the checkpoint. Because GSE289173 predates the output file,
exposure is `unknown` and this model is not eligible for held-back scoring. Its KNN/metacell path
also requires fold- and donor-isolated neighbor construction plus a no-neighbor
common lane. Checkpoint acquisition and runtime work remain approval-blocked.

The census includes a separate `scprint2_medium` record, but it remains
fail-closed. No exact medium checkpoint, weight hash, or terms/exposure bundle
has been included, and the small-v2 output file is never substituted for it.

RegFormer is pinned to Figshare file `57423100`, the only same-named v2 weight
whose companion log records a completed epoch-1 validation and best-model save
for the publication's 1.2k-gene configuration. The other v2 weight has an
ambiguous, nonterminal all-length log; the v1 weight is a different legacy
architecture. Neither is silently substituted. The exact 25-million-cell
CELLxGENE source manifest is unavailable, so GSE289173 exposure remains
`unknown` and RegFormer is not eligible as a held-back best model. The released cell
embedding workflow also postprocesses evaluation embeddings using their true
cell-type labels. That target-leaking refinement is prohibited; benchmark runs
use raw average-pooled encoder embeddings, with any supervised correction fit
only on outer-training donors as a separately named adaptation.

scBERT remains blocked before download. Its GPL-covered source, v1.0.0 paper
tag, live Weixin file ID/name/87,631,538-byte metadata, PanglaoDB corpus totals,
architecture, and preprocessing are frozen, but the separately hosted weight
has no declared terms, read-only digest, object timestamp, or unauthenticated
download URL. The historical paper checkpoint necessarily predates GSE289173;
the current mutable object is still exposure-unknown until byte identity is
proved. scBERT also has no native compact cell embedding. The common lane uses
a prespecified 200-dimensional mean of encoder tokens, excluding the synthetic
terminal token, while its released convolutional classifier remains a separate
native lane.

scCello has an exact read-only 44,049,097-byte checkpoint and SHA-256, an
read-only 22,316,072-cell pretraining corpus, and a clean chronology for
GSE289173. Its final-layer 256-dimensional CLS vector is the common embedding;
the learned cell projection is native-lane only. It still cannot be downloaded
or run because neither the official code repository nor the Hugging Face model
declares reuse terms. Written code and weight permission is required before
runtime work, and pre-cutoff development cohorts still need accession-to-corpus
UUID resolution.

LangCell is resolved to its complete public Google Drive bundle: cell and text
encoders, their configs and projections, the cell-text-matching head, and the
root config all have exact file IDs, sizes, and May 2024 modification times.
That chronology makes GSE289173 target-label-unexposed, but the Drive weights
have no declared terms and no provider content hashes, so download and runtime
work remain blocked. Geneformer initialization also means GSE136103 and the
MacParland Liver Atlas are already encoder-seen; other pre-checkpoint cohorts
remain unknown because no checkpoint-bound scLibrary UUID roster exists. The
common lane is the raw 512-dimensional cell CLS vector. The released
256-dimensional projection, supervised LangCell-CE head, and text/CTM
zero-shot route are separately named native lanes with class descriptions
frozen before held-back inference.

ChromBPNet is now specified as a locally trained, fold-specific sequence
accessibility baseline. Each major cell-state model pools only outer-training
donors after donor-state pseudobulk QC; peaks, matched nonpeaks, the bias model,
early stopping, and genomic windows obey the same donor and genomic exclusions.
Inference receives sequence only, never held-donor ATAC. Biological predictions
use the bias-free model, and the prespecified variant score is ALT-minus-REF
log total accessibility averaged over reverse complements and five seeds.
ChromBPNet does not itself supply enhancer-gene links or signed expression
effects, so any such composite is frozen and reported as a separate model.

BPNet is the matched from-scratch architecture control. It uses the same
2,114-bp windows, 1,000-bp profiles, donor-state inputs, peaks, backgrounds,
donor/genomic folds, orientations, and five seeds as ChromBPNet, but the
prespecified arm has 64 filters and no control or bias input. This avoids an
implicit CLI default of 3,088 bp and makes bias factorization and model capacity
the explicit differences. BPNet is likewise sequence-only and cannot provide
enhancer-gene links or signed expression effects on its own.

The first matched hepatocyte training output files are frozen for BPNet at
`bpnet-hepatocyte-donor0-genomic0-seed20260824-21064376` (ARTIFACTS
SHA-256 `bf0958c98cae43bc12456076ab5dae89f593c4926aa49cc0e0f6fb70925bad36`)
and full ChromBPNet at
`chrombpnet-hepatocyte-donor0-genomic0-seed20260824-21064489` (ARTIFACTS
SHA-256 `a0ed38216f0f0acc139cede8c43060af5884db473c888d6e15092f8b71b3e547`).
The independently included valid-only comparison also includes matched CNN and
transformer controls. Across 45 held donor-by-block units, the CNN reduced
profile deviance 27.18% versus the uniform baseline, BPNet 24.91%, and the
transformer 23.59%; ChromBPNet was 21.23% worse. The training-pseudobulk mean
still led regional-count correlation at 0.7083, versus 0.4841--0.5014 for the
four sequence models. This one-seed, one-split screen is selection evidence
only. The evaluator rejected the test role, and test outcomes were not read.

MPRALegNet is frozen to one exact 5,335,404-byte third-party HepG2
`test1_val2` safetensors object. Its LFS SHA-256 and MIT model-card terms are
known, but equality to an author Zenodo checkpoint is unresolved. The model
emits a 230-bp reporter-activity score. It can test ALT-minus-REF MPRA direction
after fixed sequence and strand fixtures, but it has no target gene, donor or
treatment context, and cannot enter the signed cell-type eQTL/ieQTL endpoint
alone.

ABC and rE2G are frozen as unsigned enhancer-gene link methods. ABC uses
outer-training liver ATAC with the no-H3K27ac, power-law-contact lane; its
shell-based variant helper is prohibited because it silently changes regions
and thresholds. rE2G is fixed to the released ATAC-plus-MegaMap model and its
exact feature table. Its small pickle is never loaded in a benchmark runtime:
after approval it must be inspected and converted once to inert coefficients
with prediction parity. Neither method predicts allele direction. Both remain
retrieval/link baselines or separately named composite components.

deltaSVM and direct gkm-SVM are now frozen as two score definitions over the
same locally trained LS-GKM output file for each outer fold, cell state, and seed.
Direct gkm-SVM uses an ALT-minus-REF full-window decision difference; the
deltaSVM wrapper natively sums REF-minus-ALT proxy 11-mer weights and emits one
explicit canonical sign conversion. They are useful secondary accessibility
baselines, but they are not independent model families or independent evidence.
Historical cell-line weights are excluded.

EPInformer is frozen to an exact 24-member manifest: twelve HepG2 RNA-f3
expression models and their twelve fold-matched enhancer encoders. Its included
output is an unsigned relative enhancer-attention score over a fixed
ABC-nominated candidate set. It remains link-only, HepG2-specific, and blocked
because the selected code commit has no operative license. Neither its scalar
expression output nor attention can stand in for a signed MASLD cell-type
eQTL/ieQTL effect.

Borzoi is bound to the official Calico four-member TensorFlow HDF5 ensemble,
not silently to a later framework conversion. The four members use different
initialization and sequence order on the same fold-3 test and fold-4 validation
split; they are ensemble members, not cross-validation folds or inferential
replicates. Read-only GCS generations, sizes, MD5 values, the 524,288-bp input,
7,611-track output manifest, and project-source exposure are frozen. GSE289173
is target-label-unexposed, but genomic sequence is not novel and still requires
held chromosome/LD blocks with a full-window buffer. Execution remains blocked
because the official weights have no explicit terms, the checkpoint-bound
Baskerville/Westminster revisions and training FASTA identity are unresolved,
and the published 16,384-bin tensor, released 16,352-bin prediction output, and
6,144-bin central training loss require separate native-parity fixtures. The
gReLU conversion is retained only as a port that must reproduce official
Calico outputs before it can represent Borzoi.

Enformer is split into two explicit identities. The official native Sonnet
checkpoint remains available as three read-only GCS objects but lacks resolved
weight terms and SHA-256 values. The registered CREsted conversion has an exact
SHA-256 but restricted academic, noncommercial, nontransferable terms and no
native-parity receipt. Sei is likewise restricted: its exact Zenodo archive,
21,907 regulatory targets, and 40 sequence classes are frozen, but the legacy
pickle runtime and archive SHA-256 are unresolved. Both are sequence-only
variant components and static accessibility baselines, not donor-conditioned
MASLD or standalone signed-effect models.

ChromDragoNN is a locally trained RNA-context accessibility component, not an
ATAC-profile model. Its native output is a binary open probability for a
prespecified 200-bp region. It therefore leaves the primary RNA-conditioned
ATAC profile task and remains only on secondary variant accessibility endpoints
until coordinate, expression-transform, output-class, split, loading, and
runtime fixtures pass.

The first chromatin audit narrows superficially similar models to their actual
assay roles. The current 5.8-GB EpiAgent object has unresolved terms,
digest, history, and exposure, consumes observed single-cell ATAC, and cannot
enter either current RNA-only held-back task. PeakVI has no universal checkpoint;
it is a locally trained Bernoulli ATAC latent model and requires a separate
development-only, nonchampion ATAC-cell-representation TaskSpec. cisTopic is
now fixed to GPL-3.0 release 0.3.0 at commit
`8514bdd11713da86cfe712c4234645c6b4973353`; it has no inductive held-cell
transform and remains training-partition-only. LSI is fixed to Signac
method-1 TF-IDF at `12d1ff82774890b08e55b3101f36fab83e674b88` plus
scikit-learn 1.9.0 at `77def0ed6e3beab57244885d2a584470e96c103d`, with
the exact wheel hash and transform-only observed-ATAC query path recorded.
EpiFoundation is paired RNA-ATAC rather than ATAC-only, but its released CLS
embedding requires observed ATAC. Its exact 3.6-GB checkpoint and six companion
vocabularies are frozen, while absent code terms, conflicting checkpoint
configs, pickle loading, and an unmanifested corpus keep every target-source
exposure unknown. All three remain unscheduled with empty active-task lists.
scBasset is active only as a locally trained sequence-only
profile baseline for RNA-conditioned ATAC: its outer-training-cell predictions
are aggregated within training-only cell states and repeated across held
donors, with no donor context. Its 32-dimensional cell embedding is an output
weight for training cells and is never treated as an inductive held-cell
embedding. Public tutorial weights are excluded.

The first donor0/genomic0 scBasset valid-only screen is frozen at
`scbasset-donor0-genomic0-seed11-valid-evaluation-21065092` (ARTIFACTS
SHA-256 `7ab89892e508c145f7b3bae8819af0f4695624294cfcd22f1200984cdcba8341`).
Across nine held-valid donors, five lineages, and 16,000 cCREs, scBasset had
2.18136 mean multinomial deviance per insertion versus 1.32933 for the
training-global mean and was worse in every lineage. It fails the prespecified
5% and four-of-five-lineage adoption criteria. This is a one-seed, one-fold
development result, not a held-back or universal conclusion; the negative result
is retained unchanged and the test partition remains unopened.

ChromBERT is now bound to the exact 2024 1-kb hg38 checkpoint, not its
separate 2026 safetensors release. The native model consumes ordered Cistrome
binding-status tokens for a fixed bin and has no nucleotide or allele input;
REF and ALT therefore have the same base representation. It remains a
terms-blocked region-representation comparator for retrieval only. Any allele
encoder, MPRA head, enhancer-gene link, or signed-effect predictor is a
separately trained composite and cannot inherit the ChromBERT name or exposure
claim.

EpiBRAIN remains no-license and unknown-exposure blocked. Its two live
checkpoints have identifiable Drive objects but no read-only hashes or
producing receipts, and the paper and checkpoint receipt disagree on the
cell-head dimension. More importantly, its learned contexts are fixed brain
cell types. Sequence-to-brain profiles and the published unsigned L2 variant
score cannot be relabeled as MASLD liver context, signed liver effects, or
RNA-conditioned ATAC.

For RNA-ATAC integration, the explicit-mask MultiVI adapter passed its first
five-fold same-nucleus smoke evaluation at `multivi-smoke-evaluation-21065278`
(ARTIFACTS SHA-256
`be622e9c21820c68b2082dff827884ee32318787afc3b475f6726b1e2fc13128`).
Across 39 held donors, five lineages, 10,000 peaks, 22 chromosomes, and 1.95
million out-of-fold predictions, its RNA-only inductive profile had 0.794%
lower deviance than the strongest training-global mean, improved three of five
lineages, and never degraded a lineage by more than 2%. This is below the 5%
and four-of-five-lineage adoption criteria, so it remains a one-seed smoke
comparator rather than a finalist. Held ATAC was evaluator-only and no test
outcome was opened. PeakVI's corresponding training-lineage latent-prototype
screen is frozen at `peakvi-training-context-evaluation-21065488` (ARTIFACTS
SHA-256 `a89dfe18d8fdcf84bee293e5a57ea46d1716e101fc2e68a73b8a6302b0593b71`).
It was 5.59% worse than the global-mean baseline and degraded all five
lineages, so it is a terminal negative for this transfer definition.

scGLUE's donor-held same-nucleus retrieval evaluation is frozen at
`scglue-retrieval-evaluation-21065607` (ARTIFACTS SHA-256
`515cf5905cd28f69448db5e7d55452147b2481f27e76d36c9eeee4f4bb4d9965`).
Across 39 donors and 1,000 held nuclei, linear CCA achieved donor-macro MRR
0.05032, paired scGLUE 0.04443, and unpaired scGLUE 0.02760. Both learned graph
models failed the +0.02 smoke check. Hidden RNA-ATAC pairs were evaluator-only,
and cells were not treated as independent replicates. scGLUE v0.4.1 remains
eligible for RNA cell embeddings, but
its RNA-to-ATAC decoder is labeled experimental and stays in the latent-only
category; project output files must also replace upstream dill serialization with
a controlled state-dict manifest. BABEL v1.1 is direct RNA-to-ATAC prior art,
not a local-training arm. Its code and hosted archive lack reuse terms, so it
cannot be downloaded or run; its historical loader's pre-split normalization
and clustering are also prohibited. Same-nucleus GSE296875 can support paired
objectives. GSE244832 is never converted into false cell pairs.

StabMap's independent same-nucleus retrieval screen is frozen at
`stabmap-retrieval-evaluation-21067154` (ARTIFACTS SHA-256
`34ff295ebf26f4f6c64d89d82ad7c8a7d2b889549edd635ce8f78592404951a1`).
Across 39 held donors and 1,000 hidden RNA-ATAC pairs, its donor-macro MRR was
0.02805 versus 0.05032 for training-only linear CCA, a -0.02228 difference;
top-1 accuracy was 0.00406 versus 0.01356 and top-5 accuracy was 0.02400 versus
0.05236. StabMap therefore fails the +0.02 smoke check on this development
task. Hidden pairs were evaluator-only, cells were not biological replicates,
and this result is neither held-back nor an RNA-to-ATAC profile claim.

Seurat WNN is limited to a transductive development diagnostic with truly
observed same-cell or same-nucleus modalities, such as included GSE296875
pairs. It learns modality weights and neighborhoods from the supplied geometry
and has neither a frozen RNA-only query transform nor an RNA-to-ATAC decoder.
It therefore cannot score RNA-only held-back GSE289173, enter the ATAC-profile
leaderboard, or manufacture cell pairs from GSE244832 same-donor aliquots.

EpiBERT and EPCOTv2 are registered according to their actual observed-ATAC
requirements. EpiBERT consumes 524-kb sequence, observed 4-bp ATAC, and a global
motif vector; it can support masked-accessibility or topology-matched caQTL
diagnostics but is not RNA-conditioned ATAC. Its three released checkpoints
strictly restored and produced bit-identical repeated synthetic forwards at
`epibert-released-checkpoint-forward-21065340` (ARTIFACTS SHA-256
`cff3537ef5cf023db472acff523c791b758575db428bae5d8e9ec6d973eb7a5a`).
The two pretrained objects each have 84,314,237 variables and the fine-tuned
object has 84,492,754; all encode eight attention heads despite a published
four-head pretraining command. Central target ATAC is masked, but observed
flanking ATAC remains required. This proves an executable observed-ATAC
comparator, not an RNA-conditioned or eligible for held-back scoring model. EPCOTv2's canonical checkpoint is
the paper-linked `luosanj/EPCOTv2` object, not its later byte-identical mirror.
It consumes 600-kb sequence plus observed accessibility and predicts functional
profiles, but the released object does not include the paper's trained eQTL
classifier and does not output ATAC. Both are target-label-unexposed for
GSE289173, yet neither can run its native held-back endpoint without a frozen,
topology-matched ATAC input. EPCOTv2 also remains blocked on the exact Space
code, weight, and reference-data terms.

The general DNA-language lane now has exact audits for DNABERT-2, Nucleotide
Transformer, HyenaDNA, Caduceus-PS, and Evo 2 7B. All five saw human reference
sequence, so loci and reference alleles are not called clean; GSE289173's
cell-type eQTL/ieQTL effect labels remain unexposed. Nucleotide Transformer is
a restricted comparator under CC-BY-NC-SA-4.0. DNABERT-2 and Evo require
restricted pickle-to-safetensors conversion, and Evo additionally lacks an
exact Vortex dependency pin. Caduceus has prior eQTL downstream-demo exposure;
Evo has ClinVar, BRCA, and DART-QTL downstream demos. Zero-shot/native scoring
Common representation fixtures are now frozen for DNABERT-2 (ARTIFACTS
`a0324f4b774c0af4d8bfbdc9bb4ce4514d62e704e4f0ae34bf8cd8d2c67f9545`),
Nucleotide Transformer
(`f0c65ebdb5880b9392dfd131e43b3f320df88cd8eccb7accf1f7cc9dd358e325`),
and HyenaDNA
(`9deb82d5faae1636cb6129ce6d30b70b8264c6cc30313caa8fb7506293aa12b6`).
All three passed exact-weight, REF/ALT/RC, deterministic-repeat, offline-runtime,
and outcome-separation checks. DNABERT-2 uses its upstream PyTorch attention
fallback because the archived Triton kernel is incompatible with the pinned
runtime. Nucleotide Transformer's executable learned checkpoint has
485,699,306 parameters; its upstream JAX metadata separately declares
485,729,545 without an enumerable component census, so these are retained as
different labeled counts. Locus-cross-fitted MPRA activity or direction heads
can proceed after the GSE281364 activation requirements passes. No model may claim an eQTL/ieQTL-specific supervised head or
calibration because no nonsealed signed-effect development source is currently
registered; frozen MPRA-to-eQTL transfer is reported explicitly as cross-assay
transfer.

MIDAS v0.3.0 remains a direct RNA-to-ATAC and RNA cell-state candidate, but its
released loader is not a frozen-query requirement: it uses unrestricted Torch
loading, omits ordered features from setup metadata, and rebuilds batch
dimensions from the query. The project adapter must use tensor-only state,
known technical categories, exact feature/row identities, and no query fitting.
GSE296875 may pair only through verified ordered same-nucleus joins;
GSE244832's aliquots remain separate modality batches. scMoMaT remains deferred
and nonpromotable because every supplied cell gets a learned free coordinate
and no unseen-cell transform exists. StabMap can enter RNA cell-state mapping
only with training-only references and bridges, terminal query leaves, disabled
native scaling, and training-reference-only centering; its neighbour imputation
stays latent-only for the ATAC profile task.

GEARS, scGen, categorical CPA, and CellOT are all deferred from the current
held-target perturbation task for scientific, not compute, reasons. GSE281160
perturbs noncoding loci with CRISPRi guides; GEARS requires a target gene in
both the expression and GO graphs. scGen must observe stimulated cells for the
same perturbation to estimate its latent delta. CellOT must observe the target
distribution to learn that target's transport map. CPA can compose learned
categories but has no validated representation for a wholly unseen noncoding
locus or guide. None accepts GSE313774 bulk libraries as native single-cell
inputs. Their upstream cell-level splits, full-data categories/HVGs/DEGs,
cross-role controls, and outcome-based checkpoint selection are prohibited;
cells remain training observations, while replicates, batches, or animals are
the inferential units.

scButterfly-B remains a direct RNA-to-ATAC candidate only through a new frozen
RNA-only adapter. Its high-level path preprocesses RNA and ATAC together and
its inverse-TFIDF routine requires query-ATAC statistics, so neither is allowed;
the included output is an uncalibrated continuous score over training-frozen
peaks. Cobolt remains a direct candidate through a frozen decoder that emits a
depth-free multinomial peak composition, with eval mode, posterior means,
dataset adjustments disabled, and no query appending or XGBoost fitting. Its
GPL-3.0 repository license controls despite an MIT package classifier. scJoint
is blocked on absent code/checkpoint terms and remains latent/cell-label only:
it has no peak decoder, and its native target-cell co-training cannot be used
for held or held-back queries.

scPair remains a direct retrospective RNA-to-ATAC candidate through a new
tensor-only RNA encoder/ATAC decoder adapter. The included output is its
Bernoulli peak probability `output_px`, not the source-RNA-scaled
`output_result`. Exact same-nucleus IDs are required for GSE296875 training;
GSE244832 aliquots are never paired. JAMIE remains eligible for RNA cell-state
embeddings, but its inverse-PCA cross-modal output is unbounded and may be
negative, so it leaves the primary ATAC profile task. Monae is deferred from
both active tasks because its v1.0.0 training surface exposes a CoVEL/GLUE path,
not the paper-described teacher/student contrastive architecture; Monae-E is
also transductive for held cells.

multiDGD is also removed from frozen-query cell and profile tournaments. It has
no encoder for a new cell and obtains every held representation by optimizing a
new cell-specific parameter against that cell's RNA. Its normalized ATAC decoder
output is therefore a query-optimization result, not a frozen forward pass. A
future use would require a separately approved nonpromotable diagnostic.
scDiffusion-X remains a blocked direct candidate only through lower-level
RNA-encoder, source-conditioned DDPM, and ATAC-decoder modules trained inside
outer folds. Its released translation script is prohibited because it reads
both assays and query labels; its DPM-Solver branches also omit source
conditioning. MiniAtlas weights remain excluded, and both the root MIT and
nested BSD-3-Clause notices are retained.

Regular Corgi is provisionally target-label-unexposed for GSE289173, but its
code grant conflicts between Apache-2.0 and MIT, the weight is not bound to its
producing source or modified hg38 FASTA, and the released helper uses
unrestricted `torch.load` plus `strict=False`. The project must recover an exact
tensor schema and native parity, and fit the 2,891-gene count-to-context mapper
inside development folds. Corgi+ remains a fail-closed conceptual candidate.
Its public API raises `NotImplementedError`; the weight is not bound to the
incompatible 22-, 24-, or 26-channel auxiliary modes; and its training mean
baseline, masks, transformed local RNA, and portable loader are missing. Gene
counts are not 64-bp stranded RNA coverage, and missing coverage is never zero.
Corgi+ therefore remains inactive and has no registered variant-effect task.

The official Regular and Plus objects now strictly restore and execute in the
frozen L40S runtime
`corgi-native-runtime-l40s-r8-21065022` (ARTIFACTS SHA-256
`09d8811b612cba4f75361feeba29bd0676efde83f0e4785a37990c2d2f52a0a2`).
The inclusion includes a full 524,288-bp forward producing a
`[1,22,6144]` tensor, backward/optimizer and resume probes, repeat checks, and
reverse-complement checks. This establishes executable checkpoints, not task
fitness or best openly licensed model eligibility. Regular Corgi may proceed only through
the three fold-fit 2,891-gene input-mapper arms and valid-only endpoints.
Corgi+ remains blocked because the released bundle does not define the ordered
biological 26-channel auxiliary input or resolve derivative terms.

The multimodal integration census now separates inductive held-query methods
from transductive factorization. scMVP can contribute a frozen RNA-encoder
embedding to cell-state mapping, but its released ATAC outputs consume observed
ATAC. MIRA can encode held RNA after fitting, but its paired joint representation
is only concatenated observed-modality topics and its complete license text is
unresolved. MOFA+ has no frozen assay-to-factor transform for a new donor. All
three therefore remain outside the RNA-conditioned ATAC profile leaderboard;
MIRA is terms-blocked and MOFA+ is development-only and nonpromotable.

The three released scooby checkpoints are now checkpoint-level comparators, not
interchangeable architecture labels. OneK1K is RNA-only and requires an
unbundled 10-D PBMC scPoli context. Epicardioids uses a 50-D scGLUE context
derived from unpaired RNA and observed ATAC. NeurIPS uses a missing 14-D
Poisson-MultiVI encoder trained on observed same-nucleus RNA and ATAC, and its
weight terms are undeclared. None has an eligible RNA-only liver query path, so
all three are excluded from RNA-conditioned ATAC and held-back MASLD best model use
as released. Their only current task role is a domain-labeled variant comparator.

The pinned AlphaGenome release is a static sequence-to-track comparator, not a
new MASLD context encoder. Its fixed RNA ontology has a hepatocyte track and a
generic endothelial track, but no cholangiocyte or hepatic-stellate RNA track.
It therefore cannot cover the prespecified held-back lineage roster. GET has the
opposite context problem: its native input requires observed accessible regions
and accessibility values, while GSE289173 has RNA and no ATAC. Binary-ATAC mode
does not make GET sequence-only. Decima has a native signed gene-expression VEP,
but its 8,856-task metadata has not been acquired and mapped to all required
lineages. AlphaGenome and Decima remain supported-task secondary comparators;
GET remains an observed-ATAC development comparator. Their current terms also
exclude all three from a best openly licensed model.

The typed-evidence graph lane now has checkpoint-level preflights without
pretending that any model is runnable. BIONIC and node2vec are transductive
common-support comparators; GraphSAGE is conditionally inductive only through
frozen source-independent features; and R-GCN and HGT must use a generic
source-agnostic decoder because an unseen relation type has no trained
parameters. The proposed typed transformer still has no implementation.
Evidence-count and nearest-GENCODE-v49-TSS are deterministic transparent
baselines, while equal-weight late fusion fits no fusion weights and inherits
the strictest component exposure and terms. All nine remain blocking inclusion
and development-only for graph claims because the corrected graph snapshot,
implementations, runtimes, leakage fixtures, output files, and a wholly held-back
graph source are not yet frozen. Nearest-gene is unsigned retrieval/link
evidence and never a signed eQTL/ieQTL predictor.

Variant outputs are now governed by a separate 40-model capability registry.
Every record freezes its native output, allowed endpoint, fitted-head need,
observed-target-context need, and mandatory-baseline status. Only Corgi,
Borzoi, the legacy generic `context_borzoi` conditional slot, sequence-only,
and the shuffled-context ablation currently pass the registered full-roster
signed-output requirements.
Eligibility here does not bypass license, exposure, checkpoint,
conditional-trigger, or runtime missing pieces. General DNA encoders are
representation-only; EpiBERT, EPCOTv2, GET, and ATAC-derived scooby contexts
cannot use proxy ATAC at the GSE289173 hold back; ABC, rE2G, and EPInformer remain
link-only. Nearest-gene and other link methods guard retrieval AUPRC, not
Fisher-z signed-effect correlation. Every frozen run binds its exact capability
record and a `primary_endpoint_scoring_allowed` flag. Selection then embeds the
selected and strongest-baseline records in the read-only selection record. ABC,
nearest-gene, and rE2G instead require one outcome-free auxiliary receipt each:
the receipt verifies execution, prediction schema, read-only inputs, and the
exact link-score capability, but reads no outcome and computes no development
metric. These runs stay outside `scientific_runs` and `model_candidates`. The
held-back evaluator rejects a missing, changed, observed-context, link-only, or
representation-only primary binding before it reads the hidden outcome bundle.

The perturbation census now fails closed on the current noncoding-target split.
GEARS and GenePert require gene identities; the scGPT head requires a measured
gene token; scGen and CellOT require outcomes for the held perturbation; CPA
cannot encode an unseen categorical atom; and perturblib LPM only supports new
combinations of already-seen symbols. None can represent a wholly unseen
GSE281160 locus without changing the task or leaking target identity. They stay
registered but unscheduled until a compatible task and biological-replicate
requirements are prospectively frozen. Simple replicate-level baselines remain.
Arc State's exact HepG2 transition checkpoint and required SE-600M encoder are
pinned, but their noncommercial terms exclude a best openly licensed model and corpus
overlap remains `unknown`. The released categorical transition map cannot
represent a wholly unseen GSE281160 target, while SE-600M is only a required
base component, not an independently substituted predictor. Arc is unscheduled
and valid only under a future separately frozen seen-target or held-context
diagnostic.

The legacy `context_borzoi` ID is the single generic conditional context-model
slot, not a Borzoi-bound model. It inherits no Borzoi checkpoint, terms,
exposure, runtime, output, result, or release status. It is deferred from
ordinary smoke and frozen-screen selection and can activate only after the
prespecified complementarity trigger passes and a new campaign freezes the
selected open sequence backbone, clean cell encoder, component and derivative
terms, local-training receipt, exposure audit, and trained-output hashes.
The incomplete Corgi+ release remains separately deferred.
The unresolved scPRINT-2 medium placeholder is likewise deferred until an exact
released checkpoint is bound. scFoundation remains a cell-state comparator but
has no valid representation for the current noncoding perturbation endpoint.

CREsted's released [DeepLiver Accessibility checkpoint](https://crested.readthedocs.io/en/stable/models/Liver/deepliver%20accessibility.html)
is registered but deferred from active v1 tasks. Its exact CREsted port is a
500-bp mm10 model trained on 219,823 mouse hepatocyte regions to predict 82
mouse regulatory topics, not signed human cell-type eQTL/ieQTL effects. The
TensorFlow-1-to-Keras-3 conversion has only approximate parity. Liftover does
not supply the missing species or endpoint mapping. The pinned CREsted terms
are academic-only, non-transferable, and prohibit modification, while the
converted archive has no separate resolved terms. This checkpoint is neither a
valid human variant arm nor a best openly licensed model candidate.

## GSE296875 phenotype requirements

GSE296875 is not phenotype-free. Data S1 contains donor-linked BMI,
macrovesicular steatosis, and fibrosis fields. Steatosis is observed for 38 of
39 analyzed donors and fibrosis for 37 of 39. These support donor-level,
cross-fitted development analyses and histology-stratified error audits.

They do not establish an adjudicated MASLD or MASH diagnosis, NAS, standard
fibrosis stage, alcohol exclusion, etiology, or longitudinal progression. The
39 donors span ages 13 to 75, including five donors under 18. The requirements
therefore requires explicit endpoint masks, donor-grouped inference,
leave-one-well-out transport, and an adult-only histopathology sensitivity.

Primary sources: [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12805840/),
[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296875), and
[Data S1](https://ars.els-cdn.com/content/image/1-s2.0-S0002929725004343-mmc2.xlsx).

## Safety boundary

- Planning never downloads or joins project-reserved outcomes.
- The GSE289173 cell and regulatory outcomes form one joint bundle and one
  co-unblinding event.
- A documented 2026-08-21 GEO metadata incident exposed disease-group and sex
  metadata for two samples, but no cell-state labels or eQTL/ieQTL outcomes.
  Disease/control AUROC is therefore ineligible for a held-back claim; primary
  cell-state and variant-effect evaluation requires pseudonymized features and
  an isolated evaluator with no access to this task context.
- Candidate trees, selection record, prediction commits, power decisions,
  outcome-consumption records, and releases are read-only and checksummed.
- The complete candidate universe is pre-registered. A `campaign_universe`
  output file pins every primary-scoring run and the exact per-task seed set, and
  a ledger cannot be built from any other set. Read the next bullet for what
  this does and does not prove about ordering. It has two kinds:
  `screening_frozen_screen` (wave `frozen_screen`, exactly three seeds) and
  `finalist` (waves `full_specialist_screen`, `adaptation`, or an authorized
  `conditional_model`, exactly five seeds). Neither kind may pre-register a
  `prediction_first_stress` or `sealed_inference` run.
- A selection candidate ledger names its authority. A screening ledger may
  produce a shortlist and nothing else; `five_seed_universe_complete` is false
  on it by construction. Only a finalist ledger, bound to a recursively
  verified pre-scoring finalist universe and an independently transcribed
  `reviewed_universe_id`, may back a SelectionLock or a best model. The ledger
  never accepts a caller-supplied campaign list or approval hash: it recovers
  the campaigns from the universe and rejects any drift in the observed run set
  or seed sets. `campaign_universe_sha256` exists only so a human can compute
  the identity they are approving and must never gain an in-tree caller.
- STATED LIMITATION, do not overclaim this in the paper. The universe output file
  carries `created_before_development_scoring = true`,
  `development_outcomes_used = false`, and `sealed_results_used = false`, and
  these are literal declarations written unconditionally at freeze time. They
  are NOT derived facts and nothing in the package proves them. The machinery
  guarantees that a ledger's runs and seed sets equal some frozen universe, and
  that the operator transcribed that universe's id by hand as
  `reviewed_universe_id`. It does not prove WHEN the universe was frozen.
  An operator who freezes a convenient subset universe after seeing development
  scores and then truthfully transcribes its id defeats the intent while
  satisfying every automated check. The load-bearing controls against that are
  the independent human review of the universe id and the ordering discipline
  of freezing before executing, not the internal cross-check. Establishing the
  ordering mechanically would require an out-of-band timestamped review record;
  that does not exist today.
- Development ranks the five-seed ENSEMBLE, because that is the output file held-back
  inference deploys. Macro-F1, Fisher-z Spearman, average precision, and
  relative profile deviance are all nonlinear in the prediction, so a mean of
  per-seed metrics estimates a different quantity than the metric of the mean
  prediction. Cell-state ensembling averages class probabilities and then takes
  the argmax, exactly as the held-back evaluator does; numeric tasks average the
  aligned predicted scores. Per-seed values survive only as a stability
  diagnostic carrying `used_for_ranking: false`, including the four-of-five
  positive-direction check, and never enter an ordering. Row correspondence
  across seeds is proved by an order-sensitive alignment hash, because the
  unit-set and row-set hashes are set hashes.
- Held-back-parity for that choice is literal only for `cell_state_mapping` and
  `variant_to_regulation`, the two tasks with positive held-back evaluators. For
  `bulk_state_transfer`, `typed_evidence_graph`, and `rna_conditioned_atac` the
  justification is that the ensemble is what would be deployed. Do not
  overstate it as held-back parity.
- Every cell-state development PredictionBundle must carry its
  `class_probabilities:cell_state_mapping` output file, with probabilities in
  [0, 1] summing to one and a hard label equal to its own argmax. This is
  fail-closed and it removes candidates: an adapter that emits only hard labels
  cannot compete on `cell_state_mapping`.
- A fit action emits one read-only output manifest. `checkpoint_sha256`,
  `task_head_sha256`, and `calibration_sha256` are three roles of that single
  composite bundle and deliberately bind the same digest;
  `preprocessing_sha256` comes from the prepare manifest and must differ. The
  policy is named, recorded per run, asserted on the ensemble digests, carried
  into the SelectionLock identity, and re-checked on every verification. A
  claim of independent role output files is rejected.
- A development shortlist declares the wave it came from. Only `frozen_screen`
  and `full_specialist_screen` may shortlist, each at its own prospectively
  frozen seed set, and the ledger's run stages must equal their campaign wave.
  An empty selection is rejected before anything read-only is written, because
  a published frozen directory can never be replaced. Family slots are chosen
  over models, so two recipes of one model can never take both; the Pareto and
  one-standard-error unions are deliberately not model-deduplicated.
- The conditional context model fails closed at construction, not only at the
  campaign status line. `ComplementarityEvidence` has no public value
  constructor: it can be loaded only from a recursively verified frozen tree
  that binds exactly two development shortlists, two residual bundles, and two
  development gain bundles. Authorization re-verifies that tree and every
  bound source on use, and the conditional-decision receipt binds its manifest
  hash. A hand-written mapping or provenance string cannot authorize. Compute
  is consulted only on an exact development-score tie where every tied
  candidate is `measured` under one identical compute basis; otherwise
  selection raises. `masld_bench.stacking` now supplies the scientific producer
  and recursive verifier for two-component development stack-gain bundles and
  the six-source stacking-evidence tree. It learns nonnegative convex weights
  only from donor/locus- and block-disjoint outer-fold rows, evaluates the stack
  against both open components, and carries the smaller gain forward. The MPRA
  proxy is explicitly locus-bootstrapped and cannot be called eQTL performance.
  `masld_bench.context_spec` then derives and freezes the sole
  `ContextModelSpec` from a triggered conditional decision, its bound five-seed
  sequence shortlist, and a five-seed cell shortlist. It binds the exact
  candidate recipes and fitted-state hashes. The planner rejects a ready
  `conditional_model` campaign unless that recursively verified spec path and
  identity are present; an in-memory candidate mapping is not campaign
  authority.
- Ledger freeze and every ledger verification recompute one bootstrap per
  candidate. Size that before a production campaign and never run it on the
  login node.
- Scientific planning requires a ready `DatasetActivationContract` for every
  dataset and a ready `ModelExecutionContract` with a concrete environment
  output file. Each executable model also declares its exact supported adaptation
  regimes; unsupported common, native, or model-specific lanes are omitted and
  recorded rather than scheduled. Runs bind the exact adapter action sequence,
  inputs, checkpoint, model runtime, timeout, R-provenance flag, task requirements,
  and source tree.
- Each completed action is a separately frozen adapter output. A successful
  `RunExecutionReceipt` binds every action receipt plus matching pre/post
  Resource separations. Failed attempts remain read-only; repeated execution is
  blocked until a standardized resume-checkpoint requirement is implemented.
- Submission requires the SHA-256 of the complete frozen candidate output file
  manifest. Dry-run is the default, arrays are disabled, and production
  submission is single-use through a read-only ledger containing each exact
  command and returned Slurm job ID. Current campaign templates have production
  submission disabled.
- Unknown checkpoint revisions, hashes, licenses, corpora, or exposure states
  fail inclusion. Only `clean_declared` and `target_label_unexposed` checkpoints
  can support a best openly licensed model.
- Scientific endpoint evaluation is blocked until each cell/lineage roster is
  bound to a hashed `masld-bench-evaluator-roster-v1` authority. Rosters are
  never inferred from prediction or outcome rows.
- Held-back metric bundles bind the source closure that constructs their inputs
  and metrics: tournament orchestration, task-native held-back evaluators, shared
  metric primitives, and canonical hashing. Hashing only the top-level
  evaluator file is insufficient.
- The held-back cell-state endpoint is RNA-only because GSE289173 deposits snRNA,
  not an ATAC evaluation row universe. RNA candidate and baseline predictions
  must cover identical nuclei. ATAC-only encoders remain census entries but
  require a separate development-only task and cannot compete for this held-back
  best model.
- Variant-model fitting is assay-native. GSE296875 and GSE244832 supply
  donor-held regulatory profiles; GSE281364 supplies only locus-cross-fitted
  MPRA heads and out-of-fold diagnostics. GSE281160 remains a selected-target
  CROP-seq direction check, with neither cells nor guides treated as biological
  replicates. The COLOC graph alone is not a valid ChromBPNet training assay.
- Bulk transfer keeps only inductive, assay-native mandatory baselines:
  donor-pseudobulk HVG/PCA, mean-expression, and assay-native pseudobulk.
  Harmony lacks a frozen held-query transform, while scVI/scANVI cell-count
  likelihoods are not silently reused as bulk-count models.
- A consumed held-back outcome without its complete task-native metric bundle is
  a terminal `terminal_failure_missing_sealed_metric_bundle`, never an
  invitation to repair or rescore. Adoption is possible only when the frozen
  bundle rederives every primary and secondary check metric from committed
  predictions and authorized outcomes with exactly 10,000 paired resamples,
  followed by the fixed confirmatory multiplicity decision.
- The pinned GRCh38.p14 and GENCODE v49 files are integrity-checked. Sequence
  campaigns remain blocked until a separately frozen indexed analysis copy is
  available; patch, ALT, and scaffold contigs are rejected from modeling. The
  unsubmitted `slurm/build_primary_reference.sbatch` derives and checks a BGZF
  primary-contig FASTA, FAI/GZI indexes, chromosome sizes, and filtered GTF.
- Every dataset declares both its native genome build and native annotation
  release. Unknown or mixed unresolved identities are blocking inclusion until
  an auditable crosswalk to the project reference is frozen; unmapped features
  remain missing rather than becoming biological zeros.
- Primary-record audits now pin GSE289173 to the 10x GRCh38 2020-A reference
  (filtered GENCODE v32), GSE244832 to its reported Cell Ranger 3.0.2/hg38 and
  bwa mem 0.7.17 paths, GSE281160 to a custom GRCh38 plus sgRNA-pseudogene
  reference, and GSE281367 to Cell Ranger ATAC 2.1.0/GRCh38. The latter three
  remain blocked because their exact annotation or reference-bundle identity is
  not deposited clearly enough to infer.
- The same audit corrected GSE192741 from the provisional GRCh38 declaration to
  its deposited human Space Ranger 1.0/hg19 identity. GSE256398 is pinned to
  its deposited GRCh38 feature axis, but its exact Cell Ranger reference remains
  unresolved because the GEO processing text is internally inconsistent. Its
  frozen stable-ID crosswalk, rather than a guessed reference bundle, controls
  project input eligibility. GSE267031 and
  GSE313774 now retain their disclosed aligner/vendor provenance while remaining
  blocked on exact annotation identity.

Spatial remains fully deferred in v1. The nine audited Novae, NicheCompass,
SpatialMETA, SpaMosaic, Nicheformer, GraphST, STAGATE, BANKSY, and
transcript-native baseline bundles pin source and checkpoint identities where
available, but they create no executable task: there is no adequately powered,
rights-cleared donor-level MASLD spatial endpoint or project hold-back. Every graph
must be built within authoritative participant and section boundaries, with a
technical array never treated as a donor and adjacent sections never adopted
to paired observations. Native target-fitted clustering, prototype, MNN,
registration, and graph-autoencoder paths are transductive and remain distinct
from development-frozen donor transfer. Visium spots, image-based cells,
spatial metabolites, CODEX proteins, and H&E retain assay- and topology-native
requirements. Model-specific weight terms, patent or license conflicts, conversion
parity, safe export, and runtime missing pieces remain fail-closed; none of these
models can run, train, enter the conditional model, or earn a v1 best-model claim.

The Resource side of the byte-identity separation inventories every file and
directory member in the five binding documents, the exact seven-file frozen
117-program registry, the Gene Catalog v2 requirements, all six `figures/main`
authority directories, and the current Next.js portal source/configuration.
Every public portal data file, including per-gene JSON, is protected. Only
dependency and build caches (`node_modules`, `.next`, and `out`) are outside
the authority surface.

The Cas13 side of the byte-identity separation follows the authorities named by
`Cas13_Library_Design/README.md`: the June-v9 build manifest, rebuild script,
and guide table, plus the 2026-07-31 gene-membership freeze. No final
order-ready guide authority exists yet. Directory authority membership is
part of the separation, so adding or deleting a protected file also fails.

## Verification on a compute node

```bash
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHONPATH=src python -m unittest discover -s tests -p 'test_*.py' -v
PYTHONPATH=src python -m masld_bench.cli registry validate --config-root config
PYTHONPATH=src python -m masld_bench.cli reference validate --config-root config
```

Freeze, but do not submit, a read-only inclusion review candidate:

```bash
PYTHONPATH=src python -m masld_bench.cli campaign freeze \
  --campaign config/campaigns/v1_admission.toml \
  --config-root config \
  --output-root candidates
```

The candidate contains the exact DAG, distinct sbatch headers, resource totals,
source/config/runtime selection records, model dispositions, and complete output-file manifest
hash for review. Production submission is a later explicit decision after
checkpoint, license, environment, data-join, and fixed-fixture missing pieces are
resolved. Heavy model dependencies remain isolated from this dependency-free
package.

After separately approved acquisition, stage one Geneformer variant offline:

```bash
PYTHONPATH=src python -m masld_bench.cli checkpoint stage-geneformer \
  --model-id geneformer_v1_10m \
  --source-root /approved/read-only/geneformer-source \
  --output-root /approved/checkpoint-bundles
```

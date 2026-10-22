# MASLD research model: multimodal strategy and agent handoff

Reviewed: 2026-09-07  
Purpose: Preserve the model review and develop a literature-informed strategy across the atlas, expanded at the user's request.  
Status: Technical review, not a change to manuscript claims or model-release status.

## Recommendation

**Continue the RNA-to-chromatin direction. Make a useful, transferable molecular state the next objective, rather than another broad model-family search.**

[Inference] The limiting factors now are target definition, evaluation consistency, and independent biological validation. A modest linear model that answers an important question in another laboratory can achieve the user's goal of a model other groups adopt.

The proposed research use is concrete:

> Given liver RNA from participants with similar histology, identify reproducible molecular differences that help researchers choose samples and measurements for an experiment.

This is a research-use proposal. It does not establish diagnosis, prognosis, treatment response, causal regulation, or longitudinal progression.

**Expanded ambition:** Build a model that connects a participant's molecular state to the cell types and regulatory contexts that can explain it, then identifies the measurement that would distinguish competing explanations. The chromatin predictor is the first experimentally anchored component of that model, not the intended limit of its biology.

Start with [the expanded strategy](#expanded-strategy-across-the-atlas) for the broader research directions, [the literature assessment](#what-the-online-literature-changes) for their novelty limits, and [the recommended experiments](#the-first-four-experiments) for the proposed sequence. The earlier review findings remain prerequisites wherever the current chromatin model is used.

## Review scope

The review inspected `Analysis/MASLD_Model_Benchmark/MODELS.md`, newer execution reports, the released chromatin model and selected producing scripts, and relevant manuscript authorities. Relevant primary methodological literature was also checked.

Numerical results below are **reported results**, not independently reproduced estimates. Findings about which arrays, estimators, or statistics the scripts use come from direct code inspection. Models, tests, and scientific analyses were **not run**. No jobs were submitted and no existing files were changed during the review. This document was subsequently saved at the user's request.

The receiving agent must read the five current authorities before acting: [README](../README.md), [PAPER](../PAPER.md), [STATUS](../STATUS.md), [RESULTS](../RESULTS.md), and [ROADMAP](../ROADMAP.md). This handoff does not supersede them. In particular, retain the fixed 117 Hotspot programs and the separation between the Resource and the future experimental study.

## What the newer experiments support

[Observed pattern, not a guarantee] The tested linear predictors remain strong baselines. The particular paired teacher–student construction, MOFA+ projection, BulkFormer embeddings, and Corgi adaptations tested here have not earned a replacement of the direct baseline. This is evidence about these implementations, tasks, and datasets, not exclusion of entire model families in every setting.

The strongest development result with an explicitly motivated target is prediction of a steatosis-associated chromatin component beyond recorded histology:

- **99 participants**, one sample each, GSE267145; five participant folds.
- Reported partial Spearman **0.552 [0.351, 0.704]**, adjusting for steatosis, ballooning, lobular inflammation, fibrosis, sex, and 14 RNA descriptors.
- Reported loading stability: mean absolute cosine **0.711** across folds for the selected component.
- Reported whole-matrix nested-RRR reconstruction skill **0.1194** on concentration-residualized H3K27ac log2 CPM.

The reconstruction statistic is `1 − SSE(model)/SSE(training-fold mean)`: 0.1194 means approximately **11.9% less squared error than that baseline**, on that target. It is not the fraction of all chromatin biology understood.

The newer component correlation **0.770** is not an improvement over **0.552** on the same endpoint. The target changed from a steatosis-associated chromatin component to PLS components selected for RNA–chromatin covariance. Predictable shared variation can be useful, but usefulness for MASLD research remains a separate question.

Sources: [current model registry](../../Analysis/MASLD_Model_Benchmark/MODELS.md), [Figure 7 result ledger](../RESULTS.md#model-framework-and-transfer-limits-figure-7-validated-candidate-working-figure), [chromatin-beyond-histology report](../../Analysis/MASLD_Model_Benchmark/executions/chromatin-beyond-histology-20260907T161521Z/RESULTS.md), [paired-model report](../../Analysis/MASLD_Model_Benchmark/executions/paired-teacher-student-20260907T183136Z/RESULTS.md).

## Findings to resolve before strengthening the claims

### 1. The released predictor and quoted accuracy differ

The released state heads use PCA-based ridge regression. However, the headline component correlations and confidence intervals are copied from the raw-gene direct-PLS arm. The quoted **0.1194** reconstruction result belongs to a separate nested RRR model; **0.1030** belongs to the raw-gene direct-PLS reconstruction.

The release's residual standard deviations also use direct-PLS reconstructions projected onto an all-participant basis, rather than held-out predictions from the shipped PCA procedure. These are distinct estimators and evaluation objects.

The PCA arm has promising results of its own. This mismatch does **not** establish that it performs poorly. It establishes that the current headline table and residual standard deviations do not evaluate the exact released procedure.

**Next check:** Evaluate the complete proposed release recipe within participant folds, including preprocessing, target construction, dimensionality, regularization, and reconstruction. Derive accuracy and uncertainty from that recipe. Keep component prediction and whole-matrix reconstruction as separate endpoints.

Evidence: [release construction, especially lines 75–101](../../Analysis/MASLD_Model_Benchmark/executions/paired-teacher-student-20260907T183136Z/build_release.py), [fit report](../../Analysis/MASLD_Model_Benchmark/release/masld-liver-chromatin-state-v1/weights/fit_report.json), [representation comparison](../../Analysis/MASLD_Model_Benchmark/executions/rna-representation-benchmark-20260907T192029Z/RESULTS.md), [model card](../../Analysis/MASLD_Model_Benchmark/release/masld-liver-chromatin-state-v1/MODEL_CARD.md).

### 2. Ten predictable components do not yet establish ten stable biological states

Each outer fold learns its own PLS basis. The component evaluation pools columns across folds without first establishing a common component identity and orientation. Sign alignment appears later for gene-loading interpretation, rather than for the pooled score/histology calculations.

[Inference] Inconsistent orientations could attenuate pooled associations with histology. The reported small histology correlations therefore need checking before they support the interpretation that the named states are histology-independent. This review did not quantify that effect or establish that it explains the findings.

Ten individual confidence intervals with lower bounds above zero are also not a multiplicity-controlled demonstration of ten distinct biological dimensions. Orthogonal loading vectors do not by themselves establish independent biological processes.

**Next check:** Establish component correspondence, orientation, and loading stability using training information. Evaluate held-out associations within compatible coordinates; do not align using held-out outcomes. If only the multidimensional subspace is stable, report that rather than ten named states. Declare the hypothesis family and appropriate multiplicity control, such as BH across a fixed component family or a selection-aware max-statistic permutation analysis. Carry target-construction uncertainty into inference where the claim requires it.

Evidence: [fold construction and component evaluation, especially lines 126–134 and 213–255](../../Analysis/MASLD_Model_Benchmark/executions/paired-teacher-student-20260907T183136Z/run_lane_a.py).

### 3. The cohort is development data after repeated use

The registry records at least **26 executions** touching GSE267145. A specification written before the next run does not make this cohort untouched again. Nested hyperparameter tuning does not erase target, representation, or narrative choices made after earlier results.

The GSE296875 analysis provides evidence about the procedure on another assay. It learns different components from snRNA/snATAC and therefore does not independently validate the fixed H3K27ac states. Its measured age/BMI adjustment also does not remove unmeasured confounding from GSE267145.

**Next check:** Treat existing comparisons as development evidence. Freeze the selected scientific target and predictor before a compatible independent evaluation. Repeated partitions of these same participants can assess sensitivity, but cannot create a new external validation cohort.

Sources: [registry development lane](../../Analysis/MASLD_Model_Benchmark/MODELS.md#development-lane-2026-09-07--rna-carries-chromatin-variation-beyond-histology-gse267145-development-only), [Cawley and Talbot, 2010: model-selection bias](https://www.jmlr.org/beta/papers/v11/cawley10a.html).

### 4. The genome-wide comparator uses evaluation-derived calibration

The inspected `score_genomewide.py` centers matrices using all evaluation participants and estimates the optimal prediction scale from those participants' measured chromatin. This produces a descriptive association statistic, not error from a fully fixed predictor. It also differs from the training-fold-mean centering described in that execution's specification.

**Next check:** For deployment-oriented error, estimate centering and calibration from training data and apply them unchanged to held-out participants. If retaining evaluation-fitted scaling as a descriptive comparison, label it explicitly and do not equate it with the original lane-A skill. Evaluate both models under the same target, regions, and calibration rules.

Evidence: [genome-wide scoring, especially lines 46–50](../../Analysis/MASLD_Model_Benchmark/executions/corgi-genomewide-h3k27ac-20260907T195320Z/score_genomewide.py), [comparison specification](../../Analysis/MASLD_Model_Benchmark/executions/corgi-genomewide-h3k27ac-20260907T195320Z/prespec/PRESPEC.md).

The earlier broken reconstruction was already corrected in [CORRECTION_01.md](../../Analysis/MASLD_Model_Benchmark/executions/corgi-gse267145-h3k27ac-comparator-20260907T185333Z/CORRECTION_01.md). Do not repeat its retracted comparisons as current results.

## Recommended sequence toward the user's goal

1. **Choose one primary biological use.** Keep the severity instrument as a reference. For the new model, prioritize one molecular contrast useful among participants with similar histology. Treat the steatosis-associated target and broader shared-state representation as different research questions; do not select solely by the largest correlation.

2. **Evaluate exactly what another laboratory would receive.** Fix the full prediction recipe and missing-gene handling. Resolve the estimator and coordinate issues above. Obtain uncertainty from the same procedure being evaluated, using participants as the replication unit and explicitly separating seed variability from participant uncertainty.

3. **Demonstrate usefulness beyond simple RNA summaries.** Compare against histology plus measured covariates, RNA PCs, interpretable gene-set scores, and direct ridge/RRR. The earlier selected-component analysis did not clearly separate ridge from RNA PCs. A useful endpoint is whether model-selected, histologically comparable samples exhibit a larger independently measured molecular contrast than samples selected by these simpler alternatives. RNA-only enrichment can interpret the model, but cannot independently validate its predictions.

4. **Invest in independent transfer.** Freeze a small set of candidate states and evaluate compatible paired RNA–H3K27ac samples from another laboratory. Record participant identity, composition, medication exposure, age, BMI, tissue handling, and assay batch. If a compatible cohort is unavailable, targeted molecular measurements could validate a narrower biological readout; they would not establish genome-wide prediction. Size any future study around a useful effect and its uncertainty, not the most favorable development estimate. This recommendation does not authorize data collection or collaborator outreach.

5. **Add complexity only for a measured limitation.** Reopen nonlinear modeling if reproducible residual structure motivates it, or composition-aware modeling if composition demonstrably limits transfer. Consider more independent participants and better measurements before more architectures. Reduce to an experimental gene panel only after the biological task works, with abundance and detectability included in selection.

## Expanded strategy across the atlas

### The scientific question worth organizing the model around

**What molecular configurations distinguish human livers with similar histology, which regulatory and tissue contexts account for those differences, and what experiment would discriminate the plausible explanations?**

[Inference] This question uses the atlas more effectively than asking every modality to predict fibrosis. Fibrosis and NAS are valuable recorded anchors, but asking all assays to reproduce the same coarse label can reward their redundant information and conceal the biology unique to each assay. A negative fusion result on histology does not imply that chromatin, spatial organization, or proteins contain nothing useful for a different endpoint.

An illustrative research use, **not an observed result**, is to distinguish two histologically similar biopsies whose remodeling is concentrated in different cellular compartments, identify a regulatory relationship associated with that contrast, and propose an observable readout that would discriminate a cell-intrinsic explanation from a neighborhood-associated one. The model must earn every step with data. It cannot declare a mechanism because several correlated annotations agree.

The high-impact contribution would be a reproducible biological distinction plus a prediction that transfers to another cohort or assay. A new model name, a larger latent space, or a graph connecting existing results is insufficient. A calibrated, useful distinction may be publishable even when its predictor is linear.

### The atlas is a set of connected experiments, not one fully paired matrix

The inventory below comes from current authorities and inspected execution reports. It does not convert every available file into an eligible training set. Some configuration records retain dated restrictions contradicted by later completed analyses; check the current task-specific decision before execution. Never infer pairing from matching stage, cell type, gene identifiers, or publication labels.

| Substrate | Biological units and pairing supported locally | What it can contribute | Important limit |
|---|---|---|---|
| Five-cohort Resource bulk RNA | Updated analysis has 844 participants; the benchmark's sample-based training pools are different objects | Cross-cohort tissue remodeling and transport of a fixed score | Do not equate the model pool with 844 independent training participants |
| Seven-dataset sc/snRNA atlas | 104 descriptive donors, 102 analyzed; 1,232,318 cells are nested observations | Cell identity, donor-level cell states, reference mixtures, cell-specific interpretation | Only one atlas cohort has the deposited graded-fibrosis substrate established by the recent census; predicted grades cannot supervise evaluation |
| GSE267145 | 99 participants with paired bulk RNA and H3K27ac from the same biopsies | Direct testing of participant-specific RNA–chromatin prediction | Development-exposed, one cohort; age/BMI unavailable; technical concentration and biology can overlap |
| GSE296875 | 39 donors, 68,398 same-nucleus RNA/ATAC profiles across eight wells | Cell-context regulatory prediction with exact cross-assay pairing | Donors, not nuclei, are biological replicates; incomplete histology and no explicit MASLD/MASH adjudication |
| GSE244832 | 18 donors with RNA and ATAC from different nuclei/aliquots; participant crosswalk unresolved | ATAC-native transport/readiness work; paired RNA–ATAC transport only after an authoritative crosswalk is established | Do not assume 18 paired donors from stage, labels, or a curated inferred table |
| GSE281367 | 12 donors, ATAC only | Accessibility-context reproducibility | No invented RNA context; disease-label and source-QC qualifications remain task-specific |
| GWAS, liver eQTL, corrected COLOC | Signals, loci, studies and trait families; not paired genotypes for every atlas donor | Inherited hypotheses, conditional regulatory evidence, locus-specific uncertainty | Many studies share participants or LD structure; liver-enzyme associations are not direct MASLD risk |
| GSE289173 | Registered 48-donor snRNA/eQTL/interaction-eQTL source | Potential task-specific evaluation of cellular regulatory predictions | Primary outcomes are held back under the existing project rules; this strategy does not open them or use them to select hypotheses |
| GSE192741 Visium/H&E | Five sections from four donors, same-section registration | Anatomy and technical spatial prediction development | Four donors cannot validate a high-capacity disease-niche model |
| Other spatial sources | Vu/HMSMA have unresolved donor or assay-registration joins; CosMx has targeted coverage | Qualified spatial context and measured assay observability | Arrays are not donors; imputed expression is not independently measured validation; no training rights inferred from local ZIPs |
| PXD051911 | Liver proteomes for 58 participants; exact initial liver/plasma analysis uses 41 matched participants | Protein-level state and tissue-to-blood measurement relationships | No matched bulk RNA edge; repeated visits and tissue-specific subsets change the denominator |
| GSE105127 | 19 participants, three zones each; RNA/RRBS from adjacent sections | Anatomical organization and methylation–RNA relationships | 57 zones are not 57 people; not same-cell pairing; only one established zonation cohort |

Sources: [bulk and assay counts](../RESULTS.md), [single-cell census](../../Analysis/MASLD_Model_Benchmark/executions/joint-bulk-sc-substrate-20260906T011935Z/RESULTS.md), [multi-assay synthesis](../../Analysis/MASLD_Model_Benchmark/executions/model-family-synthesis-20260906T185758Z/RESULTS.md), [GSE244832 pairing](../../Analysis/MASLD_Model_Benchmark/config/datasets/gse244832.toml), [GSE296875 metadata](../../Analysis/MASLD_Model_Benchmark/config/datasets/gse296875.toml), [paired protein analysis](continuum_translation_20260906T135538Z/IMPLEMENTATION.md).

The 117 Hotspot programs remain fixed, read-only interpretation outputs. They must not become new model-selection targets, be reweighted to improve a benchmark, or be rediscovered under another name. The inherited-evidence classes also remain independent of model predictions. The existing continuum translation already connects genetic hypotheses, programs and protein readouts; merely joining these tables again is not a new model contribution.

## What the online literature changes

The search was conducted on 2026-09-07. This is a targeted primary-literature assessment, not an exhaustive novelty certification. Public abstracts and selected full-text methods were inspected; no held-back outcome tables were opened. Published source findings are prior knowledge, not discoveries available for this paper to reclaim. Publication of a study does not establish that its individual-level data are accessible or independent of this project.

| Existing work | What is already established as a method or research direction | Consequence for this proposal |
|---|---|---|
| [Corgi/Corgi+, 2026](https://www.nature.com/articles/s41467-026-75527-2) | RNA-context-conditioned epigenomic prediction | Generic RNA-to-chromatin imputation is not the novelty; prioritize participant differences and useful conditional predictions |
| [GET, Nature 2025](https://www.nature.com/articles/s41586-024-08391-z) | Accessibility- and motif-based transcriptional modeling across cell types | A conditional regulatory model needs a stronger biological endpoint than recovering average expression; GET requires observed accessibility |
| [scooby, Nature Methods 2025](https://www.nature.com/articles/s41592-025-02854-5) | Sequence-based expression/accessibility profiles at single-cell resolution | Joint RNA/ATAC prediction is established; disease-context transport is a separate test |
| [Decima, Nature Methods 2026](https://www.nature.com/articles/s41592-026-03102-0) | Sequence-to-expression prediction across cell types and disease contexts | Disease-context sequence modeling alone is not new; fixed pretrained output contexts cannot be renamed as arbitrary MASLD participant states |
| [AlphaGenome, Nature 2026](https://www.nature.com/articles/s41586-025-10014-0) | Multimodal regulatory sequence and variant-effect prediction | Another general DNA encoder is a weak flagship; use sequence information only where it improves an appropriate endogenous endpoint |
| [BayesPrism, Nature Cancer 2022](https://www.nature.com/articles/s43018-022-00356-3) | Bulk/scRNA-based cell composition and cell-specific expression inference | Deconvolution itself is established; the advance must be a validated biological distinction and transfer beyond the reference |
| [SCENIC+, Nature Methods 2023](https://www.nature.com/articles/s41592-023-01938-4) | Enhancer and gene regulatory network inference from single-cell multiomics | Producing another regulatory graph is insufficient; test held-out molecular predictions and context-dependent relationships |
| [NCEM, Nature Biotechnology 2023](https://www.nature.com/articles/s41587-022-01467-z), [MISTy, Genome Biology 2022](https://link.springer.com/article/10.1186/s13059-022-02663-5), [NicheCompass, Nature Genetics 2025](https://www.nature.com/articles/s41588-025-02120-6) | Predictive or interpretable modeling of spatial context and niches | A ligand–receptor graph or niche embedding is not sufficient; show what neighborhood adds beyond anatomy and composition in held-out donors |
| [MASLD single-cell eQTL study, Nature Genetics 2025](https://www.nature.com/articles/s41588-025-02237-8) | Cell-state-dependent genetic regulation in human MASLD | State-dependent eQTLs are not a new general claim; do not select validation examples from this source's published successes |
| [Integrative functional genomics of metabolic liver disease risk alleles](https://pmc.ncbi.nlm.nih.gov/articles/PMC12633503/) | Existing integration of accessibility, reporter assays, eQTLs, chromatin looping and CRISPRi | The sequence-to-regulatory-network story is already populated; endogenous, context-specific predictive transfer would need its own evidence |
| [Public-data MASLD/MASH integration, npj Digital Medicine 2026](https://www.nature.com/articles/s41746-026-02352-8) | Combined bulk, single-cell and spatial analyses with network and marker nomination | Another public-data integration plus classifier is unlikely to distinguish this Resource |
| [Human liver zonation map, Nature Communications 2018](https://www.nature.com/articles/s41467-018-06611-5) | Spatially organized RNA and methylation in human liver | Recovering normal zonation is a useful control, not a novel biological discovery |
| [Perturbation-prediction benchmark, Nature Methods 2025](https://www.nature.com/articles/s41592-025-02772-6) | Strong simple baselines challenge general deep perturbation predictors | Observational pretraining and simulated gene deletion do not validate an intervention model |

**[Inference] The opening is narrower and more valuable than “integrate all modalities”: quantify when a molecular relationship changes with cellular or anatomical context, reproduce that change in independent observations, and make a better experimental decision because of it.** The directions below are hypotheses for such a contribution. None is guaranteed to be first or to yield a positive result.

### New data leads that change the feasibility assessment

**Paired liver/plasma proteomics:** A primary study reports paired samples from **7 controls and 64 MASLD participants**: *Correlation of Plasma and Liver Tissue Proteomics for Plasma Biomarkers in Metabolic Dysfunction-Associated Steatotic Liver Disease*, DOI [10.1021/acs.jproteome.5c00749](https://doi.org/10.1021/acs.jproteome.5c00749), [publisher abstract](https://www.sciencedirect.com/org/science/article/abs/pii/S1535390725005281). Its early online and issue dates differ: the indexed record reports December 2025 online and a 2026 issue. This is a concrete candidate for a second tissue-protein or liver/plasma evaluation. **[Unverified] Individual-level availability, accession, cohort overlap, assay compatibility, usable covariates, gene coverage and reuse terms were not established.** The direct ACS page could not be retrieved in this review. Do not count these 71 participants as available validation data. Do replace the reasoning “no external protein cohort can exist” with a bounded eligibility investigation.

**Paired transcriptome/metabolome/lipidome:** A [2026 medRxiv preprint](https://www.medrxiv.org/content/10.64898/2026.04.16.26351046v2.full) reports liver and blood multiomics in **211 participants**. It is a preprint, and its [data-availability statement](https://www.medrxiv.org/content/10.64898/2026.04.16.26351046v2.full.pdf) says individual clinical and matched omics data are not publicly deposited under the consent; summary data may be requested and the deposited metabolomics controls are pooled. This is **not** an available 211-person paired training set. It is a literature reference and possible future institutional collaboration, without any outreach authorized here. Summary associations cannot train or externally test an individual-level RNA-to-metabolite predictor.

These examples also show why a current literature search matters more than extrapolating permanent limits from an old local inventory. Neither example changes the paper's current data release or supplies new scientific results.

## Eight strategies, with the distinct question each would test

### 1. Resolve tissue remodeling into cellular state and composition

**Question:** When two biopsies have similar fibrosis or activity, does their RNA difference reflect cell abundance, altered transcription within the same lineage, or both?

**Potential contribution:** A transportable decomposition that lets researchers compare the same cellular process across bulk cohorts. Candidate processes include hepatocyte metabolic function, stellate matrix production, macrophage lipid handling and ductular injury. These names describe plausible biological domains, not newly discovered states or a requirement to recover four dimensions.

**What differs from completed work:** Composition adjustment, frozen-program projection, and transfer of the severity score to pseudobulk already exist. The new endpoint would be accurate recovery of a held-out donor's within-lineage expression contrast from a mixture, followed by independent molecular validation. A better fibrosis score or another deconvolution plot does not satisfy it.

**Model:** Start with a small hierarchical mixture model or an existing RNA-native deconvolution method. Use donor-aware reference expression and a small number of within-lineage factors, with study effects and uncertainty in mixture weights. Infer RNA contributions rather than silently calling read fractions anatomical cell fractions; RNA content, nuclei preparation and sorted populations change that relationship. BayesPrism is a relevant comparator, not a newly proposed algorithm.

**First experiment:** Hold entire donors and, where possible, entire studies out of reference construction. Create mixtures from held-out donor counts with known component contributions, separating changes in mixture proportions from within-lineage expression changes. Compare against fixed-reference deconvolution, composition-only models and simple pseudobulk projections. Synthetic mixtures are a necessary recovery check, not external validation of real bulk tissue. Require matched bulk/single-cell samples or independently measured cell-specific markers before claiming accurate real-biopsy decomposition.

**What would make it biologically persuasive:** The same within-lineage contrast recurs across studies and corresponds to measured accessibility or protein variation, beyond composition and technical descriptors. Interpretation through fixed Hotspot programs is post hoc only. A large score change under a reference-donor swap, unmodeled sorting or a small coverage change would argue against broad deployment.

**Priority:** High as the interpretive foundation, with real-bulk validation as the decisive missing experiment. This is not a rescue of the closed continual-integration or graded-fibrosis router studies.

### 2. Predict context-dependent regulatory relationships, not another universal gene rank

**Question:** Which regulatory relationships change with a measured liver cell state, and can the model predict those differences in new donors?

**Potential contribution:** A model returning a conditional statement about a gene or regulatory region in a defined lineage and state, with an observable target. For example, it might predict that expression of a fixed target depends on distal accessibility in one context but not another. A context-dependent association is the initial claim; causal vulnerability requires later perturbation evidence.

**What differs from completed work:** The failed signed DNA-LM/eQTL prediction and the weak addition of coarse chromatin summaries to COLOC-label ranking did not test a donor-specific, cell-context-dependent endogenous expression model. Conversely, running SCENIC+, Pando, GET or a sequence model and plotting its inferred network would largely repeat established methodology.

**Model:** First fit sparse, partially pooled regressions for measured RNA using observed cis accessibility, TF/motif context, lineage and a low-dimensional state, with narrowly specified interactions. Compare against promoter-only accessibility, distance, lineage/gene means, state alone, and a model without interactions. Build the state without the target gene and its target-specific readout features. Sequence models can supply fixed features if their native inputs are available; they do not replace the observed regulatory endpoint. GET consumes observed ATAC, Corgi consumes RNA context, and Decima's fixed tasks are not arbitrary new donor contexts.

**First experiment:** Use GSE296875 for same-nucleus development, but evaluate donor-level prediction and leave-one-well-out sensitivity. Exclude an entire donor from both modalities during fitting. Test donor-by-lineage transport in GSE244832 only under its actual different-nuclei pairing. That is a development transfer check, not a newly untouched cohort. Require an interaction model to improve held-out participant differences, not average expression across genes. Report a paired uncertainty interval and sensitivity to technical and composition effects.

**Where genetics enters:** Preserve credible-set uncertainty, allele alignment, trait identity and separate coding/regulatory evidence. A summary-data GWAS signal supplies a hypothesis about a locus; it is not a genotype for an atlas participant. A genotype-by-state effect requires actual linked genotype/expression/state data or an eligible independent interaction-QTL outcome. The held-back GSE289173 outcomes cannot select the states, genes, model, or thresholds. Known public source examples cannot become blind validation successes.

**What would make it biologically persuasive:** A predicted context interaction reproduces on independently measured expression/accessibility or an eligible functional endpoint, and changes which cell context or element a researcher should test. Do not infer enhancer-to-gene causality from peak overlap or multiply dependent evidence sources into an uncalibrated posterior.

**Priority:** Highest potential mechanistic flagship, with a harder external-validation requirement than Strategy 1. Keep it small enough for the actual donor and locus counts.

### 3. Measure how tissue neighborhoods alter molecular state

**Question:** Does a local neighborhood explain molecular variation beyond the cell or spot's own composition, anatomical location and intrinsic expression?

**Potential contribution:** Identify a tissue context in which an otherwise similar cell state has a different measured molecular phenotype. This could help select a coculture or tissue-slice experiment instead of assuming that every association is cell autonomous.

**What differs from completed work:** The Resource already has matched-null spatial organization for fixed programs and descriptive lineage localization. Neighborhood-dependent prediction is a different endpoint from Moran's I, ligand/receptor coexpression, a cluster map, or another all-program spatial projection. The existing ECM and ductular results remain reference observations, not examples selected to guarantee success.

**Model:** Compare an own-cell/own-spot model against one adding predefined neighborhood summaries and anatomical covariates. NCEM and MISTy provide relevant predictive baselines; NicheCompass is an interpretable niche comparator. Visium observations are mixtures: cell-level claims require a justified deconvolution and independent validation, not relabeling each spot as a cell. Target genes and genes used to derive the evaluated readout must not leak into the predictor or its imputed features.

**First experiment:** Use donor-held-out sections with spatially separated training/evaluation regions where appropriate. Compare against nonspatial composition, spatial smoothing and anatomy-only models. Use spatial nulls preserving relevant density, composition and autocorrelation rather than arbitrary shuffling that destroys every background structure. Test a fixed target family with donor-level uncertainty. CODEX or targeted RNA/protein is independent evidence only where registration, target coverage and donor identity are established.

**What would make it biologically persuasive:** Neighborhood information reproducibly improves prediction of independently measured target expression or protein, after the simpler backgrounds, in multiple donors and another source. Spatial proximity still does not establish signaling direction or mediation. GSE192741's four donors support a pilot; unresolved Vu/HMSMA joins cannot be overcome with a larger neural network.

**Priority:** High biological value, but limited current inferential substrate. Develop a narrow prediction test and identify the missing donor-resolved validation before training a large graph model.

### 4. Identify where RNA and chromatin cease to agree

**Question:** Which participants and molecular processes contain chromatin information that RNA does not recover, even after histology and measured technical factors are considered?

**Potential contribution:** A map of assay-specific information that tells a researcher when an RNA assay is sufficient and when measuring chromatin changes the biological interpretation. A reproducible discrepancy could motivate tests of transcriptional regulation; it is not automatically epigenetic priming, buffering, or a temporal state.

**What differs from completed work:** Global bidirectional RRR, disease-stratified prediction skill and late-fusion classification have already been examined. The new task would identify reproducible, process-specific conditional discrepancies and test their value for a concrete measurement decision. Replotting the existing global coupling versus stage is insufficient.

**Model and first experiment:** After resolving the exact-model evaluation, cross-fit RNA-to-chromatin prediction. Compare the measured chromatin with the prediction in fixed, interpretable feature families and with a calibrated error model. Estimate per-participant discrepancies without using that participant in model/target selection. Ask whether a discrepancy replicates across independent region subsets, technical replicates if available, or a compatible external cohort; account for shared technical effects between subsets. Test biological associations with a prespecified covariate model and complete family correction.

**Critical distinction:** A component of chromatin unexplained by RNA cannot be claimed as recoverable from RNA alone. When chromatin is absent, the model may report uncertainty or expected benefit of measuring it; it cannot silently impute the surprising observation. Random masking evaluates recovery under that masking scheme, not missingness caused by disease, assay failure or tissue availability.

**What would make it biologically persuasive:** The extra assay resolves a reproducible molecular distinction missed by RNA, and improves a prespecified research decision beyond a fixed choice to measure everyone or a simple coverage/quality rule. If the discrepancy follows concentration, sample quality or target instability, do not name it a biological state.

**Priority:** High near-term extension of the paired data. It directly tests a use for multimodality even when severity fusion remains unhelpful.

### 5. Model disease-associated distortion of liver zonation

**Question:** Is an apparent loss of hepatocyte metabolic function explained by a different mixture of anatomical zones, or by altered function within the same zone?

**Potential contribution:** Separate anatomical sampling from within-zone remodeling. This could change how laboratories compare biopsies, select tissue regions, and interpret metabolic programs.

**What differs from completed work:** The project already predicts zone in GSE105127 and measures zone-aware RNA/RRBS coupling. Normal zonation and its epigenetic organization were published by the source study. The new endpoint is a reproducible disease-associated departure from a fixed anatomical relationship, not improved zone classification.

**Model and first experiment:** Fit a small zone-aware hierarchical model using the 19 participants and three repeated zones, with participant effects and a restricted disease-by-zone contrast. Use anatomy labels independently of the RNA target wherever possible; a zone inferred from the same metabolic genes cannot independently validate their zonal association. For methylation, use methylated/covered counts with a coverage-aware binomial or beta-binomial model, or a justified assay-native transformation. Do not treat beta values as RNA counts or apply methylation-specific differential-variance procedures to RNA.

**What would make it biologically persuasive:** The same distortion is supported by anatomically registered expression and an independent methylation/protein readout, beyond zone composition. RNA or protein abundance is not metabolic flux. The small disease groups make this exploratory until independent disease-bearing anatomy is available; normal-reference spatial sources cannot validate a MASLD-specific change.

**Priority:** Focused biological secondary direction. Strongly relevant to liver biology, but the present cohort cannot support an unrestricted interaction search.

### 6. Learn which tissue states have usable protein readouts

**Question:** Which molecular differences can another laboratory actually measure in liver protein or blood, and which do not carry across compartments?

**Potential contribution:** A small, validated set of observable readouts for defined tissue processes, helping choose assays for mechanistic studies. This is distinct from another blood-based fibrosis classifier and does not establish secretion, liver origin, or treatment response.

**What differs from completed work:** RNA weights have already been applied to liver proteomes; liver/plasma associations and multi-tissue fusion have already been examined. The continuum translation has 41 matched initial liver/plasma participants and a 345-protein family. A new contribution must predict a fixed tissue readout in new participants or reproduce the same measurement relationship in an independent cohort, not rediscover the 49 adjusted associations.

**Model and first experiment:** Select a small target family by biological question and assay coverage before evaluation. Within the actual paired protein participants, compare regularized plasma-to-liver prediction against covariates, single-protein baselines and simple multivariate regression. Use native log-intensity models and justified missingness/censoring handling; absent proteins are not biological zeros. Keep target selection and tuning inside development. Freeze the proposed mapping before testing the newly identified liver/plasma study, if it becomes eligible.

**No invented RNA–protein pairing:** The same gene appearing in two unpaired cohorts permits functional annotation and group-level comparison, not supervised participant-level RNA-to-protein training. A protein score derived by borrowing RNA weights is not proof of within-person RNA/protein concordance. A mediated RNA-to-blood model needs a real connected pairing design or explicitly justified, testable transport assumptions.

**What would make it biologically persuasive:** The readout tracks the intended tissue process in an independent set of people beyond shared histology/covariates, and remains measurable on the intended assay. A lipid/metabolite endpoint could be more direct for a metabolic hypothesis, but the 211-person preprint is not an available individual-level training resource.

**Priority:** High translational value if external data access is resolved. This is the clearest literature-driven change to the earlier feasibility assessment.

### 7. Test regulatory context and haplotypes at a small number of loci

**Question:** Does the predicted effect of a credible regulatory variant depend on the surrounding haplotype or cellular context, and can that dependence be tested independently?

**Potential contribution:** Replace a misleading lead-SNP-on-reference interpretation with a bounded, locus-specific hypothesis specifying what single versus combined edits would distinguish. This directly addresses a dimension in the paper's scientific goal that single-signal COLOC cannot resolve.

**What differs from completed work:** The earlier marginal signed-eQTL/DNA-LM task failed. That does not validate a new haplotype model, but it motivates checking whether the target was an appropriate allele effect. The new task requires phased or otherwise justified haplotypes, fine-mapping uncertainty and an independent molecular endpoint. It is not another sweep over the same signed marginal effects.

**Model and first experiment:** Restrict to loci with adequate allele, LD and regulatory-context information. Compare reference and alternate alleles on supported haplotypes, preserving ancestry-specific LD. Compare additive and interaction hypotheses. A model-predicted double-edit interaction is only a prediction; demonstrate agreement with endogenous allelic measurements or single/double perturbations before claiming biological nonadditivity. Reporter assays, promoter editing and gene knockdown have different interventions and cannot be pooled as one ground truth. Protein-altering variants require a different mechanistic interpretation from regulatory sequence changes.

**What would make it biologically persuasive:** A context/haplotype model correctly predicts a held-out effect that a lead-only or additive baseline misses and changes the experimental interpretation. The proposed experimental panel must include discordant predictions and appropriate controls, not only favorite loci. Complete LD blocks and locus families belong in train/test separation and uncertainty.

**Priority:** High-risk, potentially high-impact future extension. Current summary statistics alone do not identify all interactions. Any new intervention study is separately scoped and does not feed back into Resource gene selection or evidence classes.

### 8. Choose experiments that distinguish explanations

**Question:** Given a genetically anchored candidate and a molecular state, which feasible experiment will discriminate the remaining explanations most efficiently?

**Potential contribution:** A model-assisted choice of cell context, regulatory element, perturbation type and readout. The useful output is a discriminating experiment with explicit alternative predictions, not a single target score.

**What differs from completed work:** The Gene Catalog already provides deterministic next-experiment rules. A new contribution would need a measured improvement over those rules, evidence-count rankings and random or coverage-matched choices. Drawing a larger knowledge graph or using an LLM to verbalize the same annotations does not establish added scientific value.

**Model and first experiment:** Represent competing explanations separately: shared-signal locus assignment, regulatory-element activity, within-lineage expression, composition, neighborhood association and downstream readout. Keep evidence dependence explicit. Initially use simple decision rules and task-calibrated probabilities where available; do not invent a joint causal probability by multiplying correlated scores. In an eligible held-out perturbation dataset, conceal whole targets or mechanisms, choose experiments using only available evidence, then compare recovered information or correctly resolved hypotheses under an equal assay budget.

**Future prospective version:** Select a small set of experiments that disagree between hypotheses, including negative controls and multiple biological replicates. A cell-autonomous versus neighborhood-associated hypothesis requires an appropriate change in cellular context; expression association alone cannot identify that contrast. Combining edits only makes sense after defining the intervention and the expected additive reference. Public perturbation data can inform method feasibility, but the repository's completed functional stress tests remain closed for Resource claims.

**Priority:** Highest long-term practical ambition, conditional on the earlier predictions being useful and on independent experimental data. The future follow-up study remains a separate paper and cannot tune the Resource model, targets or narrative.

## One model family with explicit biological connections

The recommended architecture is a small family of connected predictors, each trained on the pairing that actually exists. It can offer one coherent research interface without forcing every assay into one interchangeable latent vector.

| Component | Training observations | Output | What must remain separate |
|---|---|---|---|
| Tissue/cell-state estimator | Bulk and donor-aware sc/snRNA, with validated mixture assumptions | Continuous cellular-state and composition estimates | RNA contribution versus anatomical abundance; prediction versus measured lineage expression |
| Regulatory predictor | Paired RNA/ATAC or RNA/H3K27ac | Held-out molecular predictions conditional on available context | Accessibility versus acetylation; average profile versus participant differences |
| Spatial-context predictor | Donor-resolved spatial measurements and registered anatomy | Incremental neighborhood-associated prediction | Colocalization, prediction, and causal communication |
| Protein/readout predictor | Actual tissue/plasma protein pairs | Measurable readout with uncertainty | Tissue/plasma association versus source organ or secretion |
| Genetic interpretation | Corrected signal-specific genetics plus qualified context | Locus-specific hypotheses and unresolved alternatives | Inherited association versus disease-state remodeling and intervention effect |
| Experiment selection | Validated predictions plus eligible independent experimental outcomes | Choice among feasible discriminating experiments | Expected information versus demonstrated therapeutic benefit |

Use a shared representation only where its alignment is identifiable and improves the task. Shared gene names or matched histology alone do not identify a joint RNA–protein–spatial latent state. An unpaired model can produce a smooth alignment while erasing the biological differences of interest; held-out reconstruction on each assay does not prove that its cross-assay correspondence is correct.

Counts can use an appropriate count likelihood or established normalized expression model; donor-by-lineage ATAC uses its own coverage and normalization assumptions; methylation retains its methylated/covered denominator; proteomics uses assay-native intensities and missingness. Include shared and assay-specific variation only when each is estimable. Aggregate information at the biological unit before interpreting sample size. Large cell counts can improve state estimation without creating more independent disease replicates.

**What an advanced neural model would need to earn:** improvement on held-out participant/context differences over the corresponding sparse or low-rank baseline, with the same features and tuning opportunity. A frozen sequence backbone plus a small state-conditioned decoder is a legitimate candidate for Strategy 2; it is not automatically an improvement. A graph model is legitimate for Strategy 3 only when spatial information adds predictive value. Do not impose one expensive backbone on every assay to make the project appear unified.

## The first four experiments

These are proposed studies, **not executed or authorized compute submissions**. They deliberately leave the larger menu as options rather than launching eight workstreams.

1. **Resolve the current chromatin endpoint.** Evaluate the exact proposed predictor, establish stable target coordinates, compare a fixed biological target with generic shared components, and report source-native held-out error. This settles whether the anchor for the larger model is trustworthy. It does not repeat the model-family sweep.

2. **Test one cell-context regulatory increment.** Choose a small, biologically motivated target family using development-only evidence; define its state without target leakage. Compare promoter/distance/state baselines, a cis-accessibility model and one interaction model. Hold out donors and test technical-well sensitivity. Use GSE244832 for cross-assay transport only if its participant crosswalk, RNA annotation, and ATAC reference are independently resolved; otherwise keep the analysis ATAC-native. Advance only if the extra context improves an endogenous participant-level molecular endpoint. If it only recovers lineage averages, retain the simpler model and do not claim context-specific regulation.

3. **Test one independent biological readout.** Prioritize the new paired-proteomics eligibility investigation and the exact donor/assay joins already available. Alternatively use a compatible chromatin cohort if one is found. Freeze the target and predictor before evaluation. Compare actual sample selection or readout prediction with simple RNA summaries and histology. A retrospective validation source used to choose targets becomes development data thereafter.

4. **Test a measurement or experiment decision.** On sufficiently paired data, hide one assay, choose which samples or measurements to acquire using a rule fixed in development, and reveal the measurements for evaluation. Compare with random choice, fixed assay choice and simple technical-quality/uncertainty rules under the same budget. With only pairwise assays, evaluate that pairwise decision; do not claim optimization across the entire atlas. A prospective intervention extension belongs to its separately approved future study.

Choose practical margins in the biological endpoint before new evaluation. Report paired effect sizes and intervals, and separate a failure to demonstrate improvement from evidence excluding a useful improvement. Do not treat an estimated MDE as an extra observed-effect significance threshold. For multiple candidate targets/models, specify the family and correction or an independent selection/evaluation separation. Permutation nulls must repeat the selection steps that define the claimed statistic and use finite-sample P values, never zero.

## What would justify a major paper claim

The strongest plausible claims form an evidence ladder; none is an assertion that the result already exists:

- **Resource-level result:** Histologically comparable human livers contain reproducible molecular differences with measured cell and assay context. This requires more than a new projection of the existing continuum.
- **Predictive result:** A fixed model recovers one such difference in independent participants, beyond simpler RNA, histology and technical baselines.
- **Mechanistic hypothesis:** The difference points to a context-dependent regulatory or neighborhood relationship, with alternative explanations and a discriminating experiment.
- **Mechanistic result:** An independent intervention tests the context dependence and its readout. This is a later study if the current Resource lacks eligible intervention evidence.
- **Practical impact:** Another laboratory uses the model to choose a measurement, sample or experimental context more effectively than existing rules, and the benefit is measured.

[Inference] The combination of a replicated biological distinction and a successful experimental decision is a stronger route to field impact than a marginal increase in severity correlation. It is also more demanding. The current data can support development of the first two levels in selected assay pairs; a broad intervention model is not identified by the observational atlas alone.

## Directions that should not become another search loop

Do not repeat severity classifiers, generic multi-view factorization, model-embedding sweeps, tissue fusion, signed marginal-eQTL prediction, or normal-zone classification merely under new names. Reopen one only for a new endpoint or a demonstrated failure mechanism. Do not optimize the 117 programs to recover genetic enrichment or signed coherence after their completed tests. Keep the active continuum-coordination work separate; it already has its own requirements and holdout.

Do not infer absolute metabolic flux from RNA/protein abundance, epigenetic memory from one cross-sectional RNA/ATAC discrepancy, a cell transition from pseudotime, or treatment benefit from a genetic/remodeling sign contrast. These can motivate new experiments, but naming them in the model does not measure them.

Finally, select the paper's claim from completed evidence. Ambition should determine the quality of the experiment and its usefulness, not require a positive result, a particular number of states, a favorite gene, or a forced “groundbreaking” conclusion.

## Scientific success criteria

The next substantive advance would demonstrate that a fixed RNA-only procedure orders independent, histologically comparable liver samples by a reproducible molecular difference that matters for experimental design, beyond simpler RNA summaries.

Evidence should distinguish:

- Predictability of a training-defined assay component.
- Stability and biological interpretation of that component.
- Transfer of the fixed predictor to a compatible independent study.
- Utility for a specified research decision.

Success on one does not automatically establish the others. If no independent compatible assay is available, retain the narrower development claim and state the missing validation; do not try to replace it with more within-cohort model comparisons.

## Instructions for the receiving agent

Begin by independently checking the four concrete findings against the current files, which may have changed since this review. Do not treat this document as proof that all reported signals are invalid, or as permission to change manuscript claims, overwrite releases, install dependencies, submit workloads, or contact collaborators.

Under a review-only request, inspect and report. Under a subsequent implementation request, follow the repository's compute, write, and approval boundaries and make the smallest coherent correction. Keep detailed diagnostics in the technical layer and update the existing human authorities only when explicitly warranted by the authorized work.

**Priority:** Resolve the exact-model evaluation, develop one context-dependent molecular prediction, and validate its biological or experimental utility independently. Use the wider atlas through the explicit pairings and roles above. The eight strategies are a ranked set of research options, not authorization for eight simultaneous analysis campaigns.

## Final next-session handoff

This is the short restart record for the next session. It supersedes optimistic wording elsewhere in this document when the two conflict.

### Facts to preserve

1. **The chromatin result is not yet an exact deployable-model result.** The released paired teacher–student artifact fits a ten-component PLS target and PCA/ridge heads, but some released fields copy raw-gene direct-PLS values. The benchmark implementation uses fold-dependent PCA dimensions (`77, 77, 77, 79, 81`), whereas the released description says PCA81. The shipped NPZ also does not contain the concentration-residualization coefficients needed to score a new sample exactly. Treat historical numbers and corrected numbers as separate objects.
2. **The first experiment is a repair, not another architecture sweep.** Reproduce the historical direct-PLS and reduced-rank results from `Analysis/MASLD_Model_Benchmark/executions/paired-teacher-student-20260907T183136Z/out/lane_a_oof_reconstruction.npz` to numerical tolerance, then fit the exact deployable recipe with training-only concentration adjustment, scaling, PCA/PLS, ridge selection and reconstruction. Never use held-out-fold standardization, all-evaluation centering, or an evaluation-fitted genome-wide calibration.
3. **GSE244832 is not currently paired RNA/ATAC transport data.** Its donor crosswalk, exact RNA annotation and ATAC reference remain unresolved in the current readiness decision. Run native ATAC work if useful; do not report 18 paired validation donors unless an authoritative participant-level map is established first.
4. **The protein denominator is easy to overstate.** Local PXD051911 has 58 liver participants, but the disease-only analysis is 30 and the initial exact liver–plasma paired analysis is 41. PXD052937 is already excluded because its historical clinical labels were not source-authenticated; a newly published related study must not be counted as an independent cohort until participant overlap and individual-level metadata are proven.
5. **Protected outcomes remain protected.** Do not open GSE289173 outcomes, use them for feature/threshold selection, or let the separate follow-up study redefine this Resource model.

### Exact restart sequence

1. Read `docs/README.md`, `docs/PAPER.md`, `docs/STATUS.md`, `docs/RESULTS.md`, `docs/ROADMAP.md`, then this handoff. Check `git status` and preserve unrelated work.
2. Create a new timestamped execution directory under `Analysis/MASLD_Model_Benchmark/executions/`; do not overwrite or silently relabel a named release.
3. Write and hash a prespecified protocol before inspecting new outcome summaries. Declare the biological unit (participant), estimand, primary region/target set, baselines, nuisance variables, inner tuning rule, null, bootstrap unit, and stopping rule.
4. Dry-run the corrected chromatin evaluator on the saved arrays and frozen participant folds. Save fold-local nuisance coefficients, feature axes, loadings/scales, predictions, residuals, estimator identifiers and source-region numerator/denominator metrics.
5. Require these gates before expanding scope: historical reproduction within `1e-9` where the saved artifact permits it; no evaluation-fitted preprocessing; a resolved coordinate convention; and a complete primary held-out score with donor-level uncertainty. If the corrected result is not useful, stop and report that result rather than launching more architectures.
6. Only after the repair, run the small hepatocyte context experiment: fixed fatty-acid/peroxisomal target family, target-excluded oxidative-stress state, promoter/distal cis summaries, M4 additive versus M5 interaction, donor-grouped outer folds, and leave-one-well-out sensitivity. This is an endogenous expression prediction, not a causal claim.
7. Investigate external protein eligibility as metadata-only work first. Stop if source-authenticated participant/visit/sample keys or overlap status cannot be established.

### Operational limits

- No model fits, data downloads, installations or SLURM submissions were performed while writing this handoff.
- For the first chromatin correction, use the existing project environment, cap BLAS threads at 8, and plan a single CPU job around 8 CPUs, 16–32 GB RAM and a generous wall time after checking live partition/account limits. Do not submit an array unprompted.
- Keep all results in the technical execution layer until the user explicitly authorizes adoption. Do not edit the five Resource authorities, Gene Catalog, figures or claims from a planning result.
- The winning contribution is a reproducible biological distinction and a useful independent decision, not the largest latent space or the most modalities joined together.

# Human MASLD Resource: scientific goal, evidence review and bounded plan

**Review date:** 2026-10-20  
**Status:** User-requested record of the review and proposed wording returned in chat. This document does not replace the five manuscript authorities, adopt new results or authorize additional analyses.  
**Scientific boundary:** Standalone human MASLD Resource; all 117 programs and their weights remain fixed. The experimental screen remains a separate paper and cannot redefine Resource evidence classes, thresholds, examples or claims.

The strongest thesis is that this Resource helps distinguish a candidate inherited perturbation from a disease-state readout, and identifies the cellular context and experiment needed to connect them. The evidence supports those distinctions more strongly than it supports a continuous causal chain from inherited variation to tissue remodeling.

## 1. Scope and current scientific state

The review began with the applicable repository instructions and the five manuscript authorities: [README](../README.md), [PAPER](../PAPER.md), [STATUS](../STATUS.md), [RESULTS](../RESULTS.md) and [ROADMAP](../ROADMAP.md). Summaries were used to locate evidence; consequential observations were checked against the source tables, participant/specimen records and analytical inputs described below.

The scientific basis is the adopted 774-unit bulk analysis, corrected COLOC r2 and Catalog r9d. Explicit decisions supersede historical passages describing these analyses as pending. The current package became **2026-10-20-r3** during the review interruption, superseding October 17-r2. Its manuscript rewrite and supplementary redraws do not change the underlying scientific inputs reviewed here. The current [STATUS](../STATUS.md) and [ROADMAP](../ROADMAP.md) retain the remaining supplementary-panel qualifications; this review does not claim those panels are all complete.

The selected bulk units are conservative genetic-identity clusters, one library each, not verified independent people. Source-owned public results remain credited to their sources. Cross-sectional differences are not described as longitudinal progression.

## 2. Central thesis and scientific contribution

[Inference] Connecting inherited associations, established disease-state measurements and cellular context makes a practical distinction possible: **what molecular object warrants perturbation, what tissue response should be measured, and which alternative explanation the experiment must distinguish.**

The scientific contribution is a human MASLD Resource that makes those decisions traceable while preserving uncertainty. It demonstrates evidence interpretation and experiment specification. It has not demonstrated improved experimental success, therapeutic benefit or a general genetic mechanism of remodeling.

Each layer permits a distinct inference or choice:

| Evidence layer | What it permits | What must remain unresolved |
|---|---|---|
| Inherited association, fine-mapping and eQTL colocalization | Nominate a locus or molecular object for allele-aware investigation; distinguish coding from regulatory hypotheses. | A colocalized gene is not automatically causal, and inherited direction is not a treatment recommendation. |
| Multicohort disease-state RNA and protein measurements | Identify an observed tissue state, candidate readout and relevant histologic contrast. | Lower abundance does not establish lower per-cell abundance, impaired function or causal mediation. |
| Fixed programs and cellular/tissue context | Choose an eligible assay and plausible experimental context; distinguish lineage-enriched content from within-tissue organization. | Program membership does not establish member-gene support, transcript origin or cell autonomy. |
| Gene Catalog and linked source tables | Make the measured object, evidence state, alternative and proposed discriminating experiment accessible. | Deterministic routing is not a calibrated probability of experimental or therapeutic success. |
| Predictive-model comparisons | Characterize measurement and transfer limits for a defined task, baseline and training exposure. | Prediction does not connect the inherited and disease-state observations causally. |

## 3. Minimum supporting claims

Stable identifiers are retained throughout the review. Biological findings, Resource usefulness and methodological limits are separate claim types.

| ID and type | Population, biological n, replication unit and estimand | Decisive evidence, effect units and uncertainty | Source independence, strongest alternative and evidence that would change interpretation |
|---|---|---|---|
| **G1 — Biological candidate:** GNMT warrants an inherited-regulation hypothesis. | European NAFLD GWAS: 8,434 cases and 770,180 controls; liver eQTL: 1,183 European donors. Participant-based association summaries; exact GNMT–NAFLD shared-signal posterior. | ABF PP.H4 **0.7058315**; no SuSiE credible-set result for the direct-disease pair. Posterior support is model/prior-dependent, not a frequentist confidence interval or multiplicity-adjusted P value. [Exact pair](../../figures/main/fig2_genetics/source_tables/recurring_example_exact_pairs.tsv), GNMT/2021_34841290_NAFLD_EUR. | Cross-source participant independence is not established. Linked variants, another target or weak disease association remain explanations. Eligible conditional sharing would strengthen nomination; absent eligibility remains untestable. |
| **B1 — Biological finding:** established disease includes lower metabolic RNA/protein abundance. | Bulk: **774 genetic-identity clusters**, 138 controls/636 disease, five cohorts; not verified independent people. Protein: **46 MASLD/12 control participants**. Disease-minus-control log2 abundance differences. | GNMT RNA **−0.859**, SE **0.089**; protein **−1.136**, BH q **0.00642**. Limma-voom RNA family: 24,115 genes; adjusted protein family: 6,071 proteins. [Native RNA effects](../../RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv), [protein effects](../../Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv), GNMT rows. | Hepatocyte proportion, injury and shared severity can explain the pattern. RNA/stage analyses reuse participants; protein analyses reuse their cohort. Cell-resolved abundance or controlled perturbation would distinguish these explanations. |
| **P1 — Biological context:** the fixed ECM program describes fibroblast-enriched content and conditional spatial organization. | **55 fitted hepatocyte donors** for ordinal diagnosis; **53 paired donors** for localization; **four spatial donors**. These are different contrasts on different score scales. | Stage β **0.0837**, model 95% CI **−0.0672–0.2346**, permutation BH q **0.5777** across 117 programs. Spatial Moran excess **0.04050**, matched-gene q **0.001**, but exact donor sign-flip P **0.1875**. [H59 named rows](../../Analysis/Hypotheses/H59-20261015T112432Z/out/outcome/named_programs.tsv), [spatial observations](../../figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5F.tsv). | Composition, ambient RNA and specimen context remain plausible. Some spatial/atlas participants overlap. Independent donor-level organization would strengthen transportability; current evidence does not establish stage activation or cell autonomy. |
| **M1 — Methodological finding:** coverage and selection constrain comparisons between evidence layers. | Gene-level comparison among **3,769 COLOC-evaluated genes**, drawing on the populations above; genes are not independent biological replicates. Estimand: overlap relative to expression-stratified expectation. | **26** shared genes versus **21.2127** expected within expression deciles; Mantel–Haenszel OR **1.295**, 95% CI **0.835–2.007**, exact two-sided P **0.2755**, one prespecified comparison. Expression control remains incomplete. **935 state-only genes lack multi-signal COLOC results.** [H57 results](../../Analysis/Hypotheses/H57-20261015T112421Z/run-20261015T113107Z/results.json). | Differential observability can create apparent separation. Neither biological complementarity nor equivalence follows. Better joint coverage could change interpretation; another overlap calculation cannot supply missing observations. |
| **U1 — Demonstrated Resource usefulness:** evidence changes the proposed assay and experimental question. | **28,233 Catalog identities**, not experimental replicates. Deterministic interpretation and routing; experimental-success effect size and uncertainty are unmeasured. | Existing rules specify model, context, perturbation, readout and falsifying outcome. GNMT retains a state-associated route; IGFBP7 remains unresolved. [Routing rules](../../figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/catalog_local/passport_rulebook.tsv), [gene routes](../../figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/catalog_local/passport_next_experiment.tsv). | Traceability is demonstrated; superior target selection is not. Program membership supplies no member-gene support. Current program results also require linked source tables beyond the Catalog's historical context rows. Prospective experimental performance would be needed to extend the usefulness claim. |

## 4. Worked locus: GNMT and the metabolic abundance state

### Observation and exact genetic comparison

The exact pair is **2021_34841290_NAFLD_EUR × ENSG00000124713**, with ABF PP.H4 **0.7058315**, PP.H3 **0.0508698** and PP.H2 **0.2432987**. SuSiE produced no credible-set result for this direct-disease pair. GNMT's SuSiE PP.H4 **0.128219** from **UKBB ALT** is a different comparison and cannot be substituted for the missing direct-NAFLD result. The [exact-pair table](../../figures/main/fig2_genetics/source_tables/recurring_example_exact_pairs.tsv) identifies the native r2 source and retained alternatives.

The GWAS concerns EHR-defined NAFLD in European-ancestry participants across four cohorts, with 8,434 cases and 770,180 controls. The liver eQTL meta-analysis uses four European-ancestry microarray cohorts, totaling 1,183 donors after duplicate/relative removal. Original population sources are Ghodsian et al. ([PMID 34841290](https://pubmed.ncbi.nlm.nih.gov/34841290/)) and Broadaway et al. ([PMID 39173627](https://pubmed.ncbi.nlm.nih.gov/39173627/)); the local [GWAS registry](../../GWAS/finemapping/config/gwas_registry.tsv) and native eQTL rows provide the analysis-specific denominators.

At nonpalindromic **rs2296805**, the G allele lowers liver eQTL expression, β **−0.464886**, SE **0.0251439**, and aligns with higher NAFLD risk, logOR **+0.0595818**, SE **0.0163382**, nominal P **0.0002655**. The expression coefficient comes from microarrays; it must not be interpreted as RNA-seq log2 fold change. [Unverified] The exact normalization-unit interpretation of this eQTL coefficient was not established by the review; retain its source-native value and direction.

The current ABF anchor, **rs2296804**, is palindromic, with shared-variant posterior **0.258990**. Its reference-LD r² **0.998001** with rs2296805 does not resolve phase or a causal nucleotide. The recorded LD reference is UKBB European; matching broad ancestry does not establish perfect correspondence to the eQTL donors. [Alleles, native sources and LD](../../figures/main/fig2_genetics/source_tables/gnmt_risk_eqtl_ld_source_chain.tsv)

Prior human GNMT protein/reporter findings belong to Ji et al. ([PMCID PMC3463826](https://pmc.ncbi.nlm.nih.gov/articles/PMC3463826/)). They add regulatory plausibility, not MASLD mediation or a Resource discovery. Their participant independence from the eQTL sources was not resolved.

### Disease-state connection and contradictory evidence

GNMT, MAT1A and CYP2C19 have lower RNA and protein abundance in disease:

| Gene | Bulk disease-versus-control RNA log2FC | Adjusted liver protein log2FC; BH q | Consequential limitation |
|---|---:|---|---|
| GNMT | −0.859 | −1.136; 0.00642 | F1→F2 q **0.0562**, not supported at BH 0.05; F0→F1 loses support without GSE130970, q **0.207**. |
| MAT1A | −0.485 | −0.473; 0.0001215 | Misses the canonical RNA absolute-log2FC floor **0.50**; F0→F1 loses support without GSE130970, q **0.0892**. |
| CYP2C19 | −1.667 | −0.459; 0.02548 | F1→F2 q **0.0477**, but that mark does not persist under relevant cohort omissions. |

RNA models use cohort and inferred sex; protein models use acquisition batch, age, BMI and sex. Adjacent-stage RNA contrasts use separate complete 24,115-gene BH families, with 205/190/157/94 units for F0→F1/F1→F2/F2→F3/F3→F4. These are cross-sectional contrasts, not observed within-person transitions. [Current example and stage sources](../RESULTS.md#lower-hepatocyte-metabolic-gene-and-protein-abundance-gnmtmat1acyp2c19)

Original [protein metadata](../../data/PXD051911/meta_data.txt) and [intensity measurements](../../data/PXD051911/liver_protein_quant.txt) establish 58 distinct participant names and identifiers, all initial samples: 17 MASH, 29 MASL and 12 No_MASLD. Recorded fibrosis is 13 F0, 31 F1, 11 F2 and 3 F3, with **no F4 participants**. NAS equals its three recorded components in every selected record. All three proteins have positive measured intensities in all 58 participants. Missing protein measurements do not explain this example; advanced cirrhosis is outside its protein support.

An existing exploratory, selection-conditioned covariance analysis gives partial correlations **0.472/0.406/0.472** for GNMT–MAT1A/GNMT–CYP2C19/MAT1A–CYP2C19 within 46 MASLD participants. Baseline three-pair BH q values are **0.00243/0.00772/0.00243**. All remain supported after fibrosis adjustment, including the 12 pair-by-adjustment family; none remains supported after NAS or full-histology adjustment. This is compatible with shared severity, although adjustment does not identify the causal explanation. It remains a qualified exploratory result, not an adopted three-gene mechanism. [Covariance, participant-bootstrap intervals and adjustments](../../figures/candidates/fig34-scientific-corrections-20260929-v1/protein/protein_pair_covariance.tsv)

### Interpretation, selection and discriminating follow-up

GNMT is informative because it separates **genetic nomination, abundance and function**. The available evidence permits a provisional inherited-expression hypothesis beside an observed metabolic abundance state. It does not establish inherited GNMT regulation as the cause of that state, one inherited mechanism shared by the three genes, impaired metabolic flux, or therapeutic direction. MAT1A/CYP2C19 do not acquire inherited support through this example.

The illustrative selection cannot estimate how often inherited regulation explains remodeling, nor can it establish general early-stage behavior. The manuscript thesis survives removing GNMT: the broader contribution remains the ability to distinguish candidate perturbations, state readouts and assay/context choices.

A discriminating future experiment would separate regulatory-haplotype effects from GNMT-dosage effects, measuring GNMT RNA/protein and independently specified metabolic outputs. MAT1A/CYP2C19 would be candidate readouts, not assumed mediators. This is a review-proposed follow-up, **not GNMT's current Catalog route**, which remains state-associated stress perturbation.

## 5. Worked program: Stromal ECM/IGFBP7

### Fixed identity and disease-association limits

The fixed program is hepatocyte module 8, UID `hotspot_hepatocytes_f05c535ae5bbc0b9`, containing **49 genes**. IGFBP7 contributes **2.474%** of L1 weight. The fixed Ductular injury comparator contains 22 genes; BICC1 contributes **5.436%**. These are program labels and memberships, not gene-level spatial tests or a proposed IGFBP7–BICC1 interaction. [Native memberships and weights](../../Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes/module_genes.tsv)

H59 uses the 60-donor ambient-correctable universe. For the focal hepatocyte programs, **55 donors** are fitted: GSE202379 n=37 and GSE244832 n=18. The estimand is native donor Hotspot-score change per ordinal Healthy→Steatosis→Steatohepatitis category, not a fibrosis-stage increment. Generalized least squares includes dataset intercepts and dataset-specific residual variances; 19,999 within-dataset permutations, seed 20260923, supply P values, with BH across all 117 programs and the fixed stability rule.

H59 selects **zero of 117 programs**. ECM β is **0.083734**, 95% CI **−0.067162 to 0.234629**, permutation q **0.577711**. The two cohort slopes differ substantially, **0.0373 versus 0.7556**, with Cochran Q P **0.00676**. Ductular injury β is **0.136516**, CI **0.012117–0.260915**, q **0.323375**. The historical two-program stage selection is not current biology. [Named outcomes](../../Analysis/Hypotheses/H59-20261015T112432Z/out/outcome/named_programs.tsv), [H59 report](../../Analysis/Hypotheses/H59-20261015T112432Z/RESULT_H59.md)

The adopted bulk transport arm has **364 genetic-identity units**, 112 eligible programs and a 224-test program-by-axis BH family. Within-cohort OLS with HC3 errors and REML/Knapp–Hartung meta-analysis yields no family-significant fibrosis or NAS association. The ECM fibrosis coefficient is **0.320**, CI **0.103–0.538**, q **0.214**; composition adjustment reduces it to **0.170**, CI **−0.056–0.396**, q **0.421**. Historical significance requires an identity-unresolved fourth cohort and a variance convention whose alternative removes all calls. [Current bulk transport evidence](../RESULTS.md#single-cell-programs-and-composition)

### Localization and ambient RNA

Same-atlas donor-paired localization places ECM higher in fibroblasts than hepatocytes: difference **2.271225**, model 95% CI **2.062146–2.480304**, HC3 BH q **6.88×10⁻¹⁹**, **53 paired donors**, four datasets. This uses frozen L1-weighted gene-z scores from mean per-cell log1p(CP10K), a different score from H59. BH covers ten prespecified lineage comparisons. Ductular injury is highest in cholangiocytes, difference **1.573724**, 51 donors. [Exact localization contrasts](../../Analysis/SingleCell/candidates/cross-lineage-specificity-complete-atlas-candidate-2026-08-15-v2/results/hero_lineage_contrasts.tsv), `all_annotated_cells`, minimum_cells=50.

This is descriptive same-atlas localization, not proof of transcript origin or cell autonomy. In the common-universe equal-variance ambient comparison, ECM falls **0.350898→0.127891**. Its paired change **−0.223006** has nominal P **0.006417**, but BH q **0.478071** and HC3 q **0.337646**. Under H59, ambient-corrected ECM β is **0.051817**, q **0.55887**. Thus contamination is a substantiated alternative, but attenuation does not pass family correction and is not proof of explanation. Eleven fixed T-cell programs cannot be transported from integer counts; GSE189600 and Liver Atlas cannot be ambient-corrected under the existing method. [Ambient effects](../../figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/ambient_program_effects.tsv)

### Spatial observations and original specimen identities

The spatial ECM score measures **26/49 members carrying 54.75% of original weight** in GSE192741. After adjustment for all 16 cell2location abundance factors, library size and detected genes, the donor-first Moran excess is **0.04049714**. The 9,999 matched-gene draws support this sampled-tissue excess, BH q **0.001** across the two fixed programs within source. Donors and arrays stay fixed during gene randomization: this is not donor-population inference.

Across **five sections from four donors**, three donor excesses are positive; exact one-sided sign-flip P is **0.1875**, with minimum attainable P **0.0625**. The matched-gene 95% reference range is **−0.02040 to 0.02289**, a null range rather than a donor confidence interval. Vu supplies ten positive ECM arrays, excess **0.050337**, q **0.0006**, measuring 25/49 members and 52.82% of weight, but has no authoritative donor key. It remains source-dependent. [Current Figure 5F source observations](../../figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5F.tsv), [rendered panel](../../figures/main/fig5_molecular_context/panels/fig5f_spatial_program_calibration.pdf)

The original [Guilliams Table S6](../../Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/source/Guilliams_TableS6.xlsx), sheet1 rows 147/153/155/166, changes the biological interpretation of replication:

| Donor | Source specimen context | Recorded steatosis/fibrosis |
|---|---|---|
| H35 | Gastric bypass, obesity | 35–40% steatosis; no fibrosis. |
| H36 | Colorectal-cancer liver-metastasis resection; six standard alcohol units daily | No steatosis; mild pericellular/periportal fibrosis. |
| H37 | Gastric bypass | 70% steatosis; clear pericellular fibrosis. |
| H38 | Cholecystectomy | 10% steatosis; mild pericellular fibrosis. |

These are heterogeneous human-liver contexts, not a qualified MASLD case-control or stage-replication cohort. H37/H38 explicitly have both snRNA-seq and Visium and remain in the [atlas donor metadata](../../Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv). Spatial and atlas context therefore share participants. Their contribution to the original hepatocyte Hotspot discovery was not independently established in this review. Cell2location factors and spatial scores also derive from the same RNA matrix, so composition adjustment is not orthogonal biological validation.

The Ductular/BICC1 comparator does not exceed its matched-gene null: GSE192741 excess **0.015322**, q **0.0793**, donor P **0.125**; Vu excess **0.002655**, q **0.3158**. This is **indeterminate**, and no prespecified test establishes a difference between programs. Both examples were motivated by an earlier stage fit, even though their definitions were fixed before spatial evaluation. They cannot represent the prevalence of successful programs across all 117.

The justified experimental choice is to test a **multicellular tissue readout and its cellular origin**, rather than infer an IGFBP7 mechanism from its program label. IGFBP7 remains unresolved in the Catalog. An independent spatial donor analysis could strengthen transportability; it would still not identify a causal member gene.

## 6. Weakest connection and bounded next-step plan

The weakest connection is **inherited GNMT regulation → the established metabolic tissue state**. Existing cross-sectional, partly unpaired measurements cannot distinguish mediation from linkage, composition or shared injury.

The following are proposed additional analyses or future experiments, not work submitted by this review. None is required to rescue the manuscript's present thesis. Reuse completed checks and stop when missing evidence prevents the stated comparison.

| Unresolved claim | Exact input and modality-appropriate method | Assumptions and biological unit | Observable decision and effect on the manuscript |
|---|---|---|---|
| **G1:** Does the exact GNMT disease signal share an eligible conditional eQTL signal? | Existing `2021_34841290_NAFLD_EUR` summary statistics named in the allele/source table; `data/broadaway_eqtl/Liver_eQTL_Meta_Isolated_20240530.tar.gz`; saved European LD named in that table; existing SuSiE/coloc procedure. Any extraction/fitting is a separate compute-node analysis. | Participant-based summary inference; correct allele alignment, applicable LD and identifiable signals are required. Preserve the original credible-set rules and per-signal shared-posterior eligibility of at least 0.5 on each side, as implemented in [coloc_signal_eligibility.R](../../GWAS/finemapping/src/coloc_signal_eligibility.R). | Eligible PP.H4 >0.5 under recorded priors and predeclared sensitivity checks strengthens nomination, not causality. Eligible below-threshold support narrows the connection. Missing credible sets remain untestable. Do not search settings until positive. |
| **P1:** Does ECM organization generalize across people? | Only if authoritative donor/section joins and usable histology become available: freeze the same program weights, assay eligibility, spatial matched-gene nulls and all-cell adjustment; use donor-level resampling. | Donors, not sections or spots, are replicates. Suitable clinical context and sufficient program coverage are required. Prespecify donor-level uncertainty and correction across the two fixed programs before examining outcomes. | Persistent adjusted excess with donor-level uncertainty excluding zero supports transportability. Disappearance after adjustment favors composition. Inadequate coverage or precision remains untestable/indeterminate. A MASLD stage claim additionally requires an adequately powered donor-level stage contrast. |
| **G1→B1:** Can inherited regulation cause the proposed state? | Future allele/haplotype and GNMT-dosage perturbations with paired RNA/protein measurements, independently specified metabolic outputs and rescue. | Independent biological cultures/donors, not wells, are the replication unit. Specify a meaningful effect, multiple-testing family, power and adequate perturbation controls before execution. | Regulatory effects alone support regulation. Reproducible downstream effects with rescue strengthen the connection. Adequately measured nulls narrow it. These experiments remain outside the Resource and cannot redefine its fixed programs or adopted classes. |

### Role of predictive models and the molecular continuum

Predictive models serve a **supporting role in measurement and transfer limits**. The current degree-of-freedom-matched fibrosis comparison favors plain expression over composition-plus-expression: Δ macro-Spearman **−0.127784**, 95% interval **−0.175510 to −0.081756**, on **521 samples without a donor key**. The interval uses sample resampling, not verified independent-donor resampling. [Current Figure 7B values and native sources](../../figures/main/fig7_model_framework/data/panel_deep_20261017/fig7B.tsv)

The separately instructed reviewer inspected 15 saved cell-model fold receipts: training/query studies were disjoint, HVG fitting was training-only, and no query adaptation was recorded. A blanket leakage objection to that comparator is therefore unsupported. The same reviewer checked source-exposure dispositions: TranscriptFormer has five declared-clean studies with 80 donor identifiers and two seen studies; Geneformer has no qualified clean study among the seven. Broad-cell classification and disease prediction remain different estimands. [Exposure correction](../../Analysis/MASLD_Model_Benchmark/executions/fig7-source-exposure-correction-20261014/README.md)

RNA→H3K27ac reduced-rank regression provides within-cohort cross-assay predictability in 99 participants, skill **0.122141**, 95% interval **0.082705–0.159196**, with participant-held folds and training-participant baselines. This does not establish causal coupling or independent external generalization. Tissue-level severity performance must not be attributed to sequence models that were not scored on that task.

Recorded histology remains the severity anchor. The externally anchored molecular continuum describes within-stage transcriptional variation; its relationship to RNA-derived outcomes is not independent confirmation, temporal ordering or resolution of metabolism, inflammation and fibrosis into separate axes. The main thesis and experimental routing do not require a predictive-model headline or a resolved aspect-specific mechanism.

## 7. Rewritten overall goal and plan

The following proposed text preserves the existing scientific scope. It is a replacement wording proposal, not an edit to the manuscript authorities.

> The goal is to build a human MASLD Resource that helps investigators distinguish inherited molecular candidates from the tissue responses of established disease, then choose the biological context and measurements needed to test a connection between them. Its contribution is a defensible path from human observation to a discriminating experiment, with the strength and limits of each connection explicit.
>
> We will organize the manuscript around three questions. Which inherited signals nominate a molecular object for investigation? Which RNA and protein measurements describe disease-associated remodeling? Where are those candidates and responses observable within cells and tissue? Agreement between these layers will be reported as convergence of measurements only when the relevant identities, contrasts and assay eligibility support that interpretation. Shared participants and measurements will remain visible.
>
> The GNMT locus will illustrate this reasoning. Its provisional genetic evidence motivates a regulatory hypothesis, while lower GNMT, MAT1A and CYP2C19 abundance describes a metabolic tissue state. The manuscript will distinguish this association from a demonstrated inherited mechanism, metabolic dysfunction or therapeutic direction. The fixed Stromal ECM program will illustrate a different decision: how to distinguish lineage-associated content and tissue organization from a reproducible disease-stage response.
>
> All 117 programs and their weights will remain fixed. Recorded histology will anchor severity; the molecular continuum will describe within-stage transcriptional variation without supplying temporal order or independently confirming the RNA findings. Protein, chromatin and spatial measurements will retain their native meanings and limitations. Predictive-model comparisons will support conclusions about measurement and transfer, with training exposure reported beside performance.
>
> The Gene Catalog and linked source tables will make each conclusion usable by recording the measured object, evidence state, unresolved alternative and next discriminating experiment. Completion requires coherent claims, current figures and explicit dispositions for unavailable evidence. Additional analyses will be limited to questions that could materially change interpretation. The Resource will remain scientifically complete regardless of the timing or outcome of the separate experimental screen.

## 8. Three review–verification–revision loops

This compact issue table records substantive changes to the proposed interpretation. No adopted result, gene class, program definition or manuscript authority was changed by the review.

| Loop / claim ID | Specific concern | Evidence checked | Disposition | Resulting change |
|---|---|---|---|---|
| **1 — B1** | A different roster could inflate or misstate biological n. | Exact count-input roster, its recorded checksum, unique unit IDs and fitted library IDs. | Resolved: `identity_resource_qualified/bulk_primary_sample_ids.tsv` contains 774 unique clusters and matches the model exactly; the initially inspected 738-row file was a different sensitivity input. | Retain 774 and the explicit limit that clusters are not verified independent people. |
| **1 — G1** | A cross-study maximum could be mistaken for an exact replicated signal. | Native r2 rows and exact-pair/allele tables. | Direct-NAFLD ABF support and ALT SuSiE support belong to different comparisons. | Keep the provisional direct pair; do not label its missing SuSiE credible set a tested genetic negative. |
| **1 — M1** | Untested genes could be counted as genetic negatives. | Native bulk/COLOC join and expression-decile counts. | Independently reproduced 14,704 interface rows; 3,769 evaluated genes; 426 genetic, 220 state, 26 shared; expected overlap 21.2127377; 935 state calls without multi-signal results. | Remove sparse-overlap arguments for biological complementarity or opposition. |
| **1 — P1** | Historical stage selection could be mistaken for current evidence. | Corrected donor metadata, fixed scores, H59 outcome family and rendered Figure 4C. | Current H59 selects zero of 117; only two datasets inform the focal slopes. | Retain localization and tissue-context observations, not a robust single-cell stage-response claim. |
| **2 — B1** | Concordant abundance could imply one inherited mechanism. | Native protein effects, participant histology and full covariance adjustment rows. | Common severity, injury and composition remain plausible; NAS-adjusted covariance loses support. | Retain abundance-state interpretation; remove mediation, autonomous-function and universal early-stage implications. |
| **2 — P1** | Spatial significance could imply independent MASLD replication. | Original specimen workbook, atlas donor joins, four donor effects, exact sign enumeration and rendered Figure 5F. | Sampled-tissue gene-null support survives; donor-population support and independent disease replication do not follow. | State heterogeneous clinical context, shared participants, partial program coverage and donor uncertainty. |
| **2 — P1** | Ambient attenuation could prove contamination explains the program. | Common-universe paired effects, H59 ambient-corrected estimates and complete family correction. | The attenuation is not family-significant; some datasets/programs cannot be corrected or transported. | Retain contamination as an alternative without claiming it is the established explanation. |
| **3 — U1** | Integration could imply validated target selection or complete per-gene program integration. | Actual experiment routes and Catalog/source-context contents. | Routes are research proposals; current 117-program results require linked tables beyond the Catalog's 14 historical program-context rows. | Claim demonstrated traceability and assay choice, not experimental-success prediction. |
| **3 — M1** | Historical model comparisons could carry the biological thesis. | Current matched Figure 7 comparisons, native prediction sources, saved split receipts and training-exposure records. | Retain task-specific supporting findings; reject an unsupported blanket leakage objection where receipts demonstrate separation. | Models support measurement/transfer limits, not the inherited-to-state causal chain. |
| **3 — G1/B1/P1/U1** | The strongest example might be essential to the thesis. | Removal-of-GNMT coherence check against the remaining claims and existing worked Catalog entries. | The Resource remains useful for distinguishing perturbation objects, tissue readouts, assay coverage and experimental contexts. | GNMT illustrates the Resource without determining its validity or generality. |

A separately instructed reviewer, working in a separate context under the same source and permission boundaries, supplied source-supported objections. This was not independent biological validation. Review agreement did not substitute for source records, calibrated controls or biological replication.

One final focused pass reconciled the overnight package change and rechecked the affected qualifications. No additional analytical reruns were used to seek a positive result. Remaining causal, donor-identity and metadata uncertainties have the dispositions stated above.

## 9. Checks completed and not run

### Completed during the review

- Read applicable instructions and the five manuscript authorities; resolved historical/current conflicts using explicit decisions and underlying evidence.
- Inspected exact GWAS–eQTL comparisons, allele alignment, LD qualification and native effect rows; preserved source-native eQTL units.
- Reconstructed the adopted bulk gene/call totals from its small result table: **24,115 tested; 1,854 called; 1,333 up/521 down**. This was table arithmetic, not a model refit.
- Confirmed 774 exact count-input roster rows and unique units, matching all fitted library IDs; roster SHA-256 **5983eb6d1c01ce7d95d6b61a874bec91f57d5ef98a5e9f5b67bf5520fc865558**.
- Inspected native participant/histology/protein joins, fixed program weights, donor scores, current program-family outcomes and model implementation.
- Recomputed the expression-stratified expected overlap and spatial donor means/sign-flip probabilities using lightweight small-table arithmetic.
- Inspected original Guilliams specimen records and retained atlas donor joins, confirming heterogeneous context and H37/H38 overlap.
- Extracted current PDF text and inspected rendered panels. Primary review directly rendered/inspected Figures 5F and 6C; the program reviewer inspected Figures 4C and 5F, and the model reviewer inspected current Figure 7 proof material.
- Inspected experiment-routing rules, saved model split/exposure records and current matched model-comparison sources.

### Not run

- Bulk RNA-seq, protein, stage, composition or continuum refits.
- GLS/permutation reruns, decontX reruns or spatial spot-matrix reanalysis.
- New conditional colocalization, LD reconstruction or new-cohort analysis.
- Alignment/read inspection or raw mass-spectrometry reprocessing; no necessary claim required reopening those measurements during this review.
- Model training, new prediction runs or experimental validation.

No jobs were submitted and no manuscript authorities or adopted results were edited. This Markdown file is the user-requested saved review. Documentation checks for this file are limited to structure, local references and the 250–400-word goal/plan requirement; they do not constitute scientific revalidation.

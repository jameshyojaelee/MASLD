# AlphaGenome Atlas × MASLD Resource
## Detailed Public-Data Analysis Protocol and Agent Handoff

**Version:** Expanded replacement plan, September 9, 2026  
**Objective:** Add biologically informative, reproducible analyses to the current MASLD Resource manuscript using public data.  
**Status:** User-approved implementation protocol. Execution outcomes are recorded separately; this document does not establish new biological results.

**Standalone document destination:**  
`docs/technical/MASLD_ALPHAGENOME_RESOURCE_PLAN_2026-09-09.md`

**Saving status:** Saved during the user-authorized implementation. The initial protocol below is retained independently of execution outcomes.

---

## 1. Scientific Contribution

The proposed contribution is a **mechanistic interpretation of inherited MASLD signals in the context of measured human liver biology**.

The Resource should help a researcher answer:

> At this genetic signal, which molecular effects are predicted, which are supported by measurements, where are they observable, and what remains unresolved?

The strongest potential result is not simply that AlphaGenome assigns high scores to known variants. It is that integrating predictions with measured liver biology reveals a useful distinction, such as:

- Different plausible molecular mechanisms among variants within the same genetic signal.
- Predicted RNA-processing effects that total-gene expression analysis does not resolve.
- Inherited regulatory candidates whose cell context differs from the location of established disease remodeling.
- Genetic and disease-state directions that disagree, without assuming compensation or therapeutic benefit.
- Specific regulatory features associated with the local information recovered by the existing chromatin predictor.

These are **hypotheses**, not expected positive findings.

### Primary scope

Use existing public genetic results, donor-aware molecular analyses, and targeted AlphaGenome Atlas queries. Produce manuscript-ready candidate tables, analyses, and figures.

### Secondary scope

Evaluate public junction-level data, local saturation predictions, tissue comparisons and model interpretation when the necessary inputs already exist or require bounded processing.

### Deferred work

New cohorts, controlled-access genotypes, prospective perturbations, rare-variant burden analyses, haplotype-interaction validation, and training a new model on Atlas outputs are not submission dependencies.

The existing Resource remains publishable on its own scientific foundation. This extension must earn its space through an interpretable result.

---

## 2. Corrections That Govern the Design

### 2.1 Expression-colocalized loci are a selected subset

A lightweight recount of the existing primary-pair table reproduced:

| Existing genetic context | Gene–study pairs |
|---|---:|
| Direct MASLD/liver-fat group | 29 |
| Liver-enzyme group | 787 |
| Total | 816 |

These are neither independent-locus counts nor the complete genetic architecture of MASLD.

The 29-pair group was selected through expression colocalization. It therefore cannot establish the prevalence of splicing, coding, or other mechanisms across all MASLD risk loci.

**Design consequence:** Maintain separate genetic universes, defined in Section 3.

### 2.2 Static cell context and disease-associated accessibility have different coverage

The existing tables cover five static lineages: hepatocyte, stellate, macrophage, cholangiocyte and T/NK.

The inspected disease-accessibility summaries support a narrower analysis:

| Lineage | GSE244832 primary donors | GSE281367 primary donors | Jointly testable peaks | Existing supported peaks |
|---|---:|---:|---:|---:|
| Hepatocyte | 5 normal + 9 MASH | 6 normal + 6 MASH | 30,650 | 0 |
| Stellate | 5 normal + 9 MASH | 6 normal + 6 MASH | 9,264 | 4 |
| Macrophage | 4 normal + 8 MASH | 6 normal + 6 MASH | 0 | 0 |

These values were read from existing result tables, not recomputed from fragments.

**Design consequence:** Five-lineage static annotation is allowed. The primary disease-change analysis is restricted to hepatocytes and stellate cells. Other lineages remain unavailable for that endpoint.

### 2.3 Several AlphaGenome scores discard direction

Maximum-absolute splicing, polyadenylation and contact scores indicate predicted disruption magnitude. They do not necessarily identify an increase or decrease at a particular molecular event.

Active-allele scores represent the larger predicted activity of the reference and alternate alleles. They are not synonymous with reference-allele baseline activity. [Official scoring definitions](https://www.alphagenomedocs.com/variant_scoring.html).

**Design consequence:** Every score receives an explicit signed/unsigned designation. Event-level direction is reported only when the returned data support it.

### 2.4 Multiple predictions are not independent evidence

Predicted RNA, accessibility and histone effects from AlphaGenome share a model and training information. Their agreement does not count as three independent validations.

AVI also incorporates information beyond molecular predictions, including conservation and protein-related features. Its feature contributions explain the score rather than establishing a biological mechanism. [Atlas technical paper](https://storage.googleapis.com/deepmind-media/DeepMind.com/Blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/alphagenome-atlas.pdf).

---

## 3. Genetic Universes and Biological Units

### Universe A: Expression-colocalized direct-trait signals

Start with the existing 29 direct MASLD/liver-fat primary gene–study pairs.

Retain the exact selected signal pair and all available variant posterior information. Do not choose a different signal because its AlphaGenome result is more compelling.

**Permitted interpretation:** Molecular annotation of the existing expression-colocalized direct-trait set.

### Universe B: Other direct-trait fine-mapped signals

Inspect the existing direct-trait GWAS fine-mapping outputs for signal-resolved credible sets that do not require expression colocalization.

Include a signal only when its existing source establishes:

- Study, phenotype and ancestry.
- Genome assembly and variant alleles.
- Signal identity.
- Posterior definition and credible-set membership.
- A traceable current analysis source.

**Permitted interpretation:** A broader survey of existing fine-mapped direct-trait signals.

**Fallback:** If these requirements cannot be met from existing outputs, retain Universe A and explicitly restrict every landscape-level conclusion to that selected universe. Do not rerun fine-mapping as part of this extension.

### Universe C: Liver-enzyme signals

Process the existing 787 enzyme gene–study pairs separately.

They provide a comparison and a broader liver-trait resource. They are not negative controls for MASLD and must not inflate the direct-MASLD denominator.

### Dependence handling

Keep the following identifiers distinct:

- Variant.
- Gene–study pair.
- GWAS signal.
- eQTL signal.
- Genomic locus.
- Analysis block used for resampling.
- Biological participant.

Use existing locus and LD information to group dependent signals. If authoritative independent-locus boundaries are unavailable, use conservative overlapping analysis windows for resampling and label them as analysis blocks, not independently established loci.

Within a block, repeated genes, studies and signal pairs remain together. They must not be counted as independent discoveries.

---

## 4. AlphaGenome Retrieval and Harmonization

### 4.1 Query policy

Query unique normalized variants once, then join their predictions back to all relevant signals.

The initial technical pilot uses the first five eligible direct-trait loci in genomic-coordinate order. This selection is independent of AlphaGenome predictions and disease outcomes.

Retrieve only the molecular channels needed for the planned analyses. The API supports selected variants and intervals; a whole-Atlas download is unnecessary. [Atlas API](https://www.alphagenomedocs.com/api/atlas.html).

### 4.2 Coordinate and allele requirements

For every variant:

1. Preserve the original identifier, assembly, chromosome, position and alleles.
2. Normalize into the Atlas assembly using the existing authenticated reference and mapping resources.
3. Verify the reference allele against that assembly.
4. Record strand changes explicitly.
5. Preserve multiallelic identity; an alternate allele is part of the variant key.
6. Reject ambiguous or failed mappings from numerical comparisons.
7. Retain their excluded posterior mass.

Do not infer strand from association direction. Do not silently repair a reference mismatch by swapping alleles.

SNVs are the primary analysis. Score indels only where the queried Atlas artifact explicitly supports them; otherwise record them as uncovered.

### 4.3 Gene and transcript identity

Retain both the Atlas-native gene/transcript annotation and the project annotation.

Gene-level mapping may use stable identifiers with recorded version differences. Transcript-level claims require a compatible exon/junction structure, not merely a similar gene symbol.

A variant-centered accessibility or histone score has no automatic target gene. Gene assignment requires a separately supported link.

### 4.4 Track selection

Select tracks using metadata before examining results:

- Human tracks only for primary analyses.
- Primary liver or exact relevant cell types before surrogate cell lines.
- HepG2 and other transformed lines reported as separate contexts.
- Exact biosample, assay, histone mark and scorer retained.
- No maximum-across-tissues selection for a liver-specific claim.
- No substitutions for missing stellate or cholangiocyte outputs.

Retain individual tracks. For replicate tracks measuring the same endpoint in the same biosample, use their median as the display summary and show their range. These tracks are not biological replicates.

### 4.5 Score semantics

For each scorer record:

- Molecular quantity.
- Spatial mask and window size.
- Gene/event specificity.
- Transformation and native units.
- Whether sign is retained.
- Whether the output is a raw effect, quantile, active-allele score or AVI attribution.

Do not compare raw magnitudes across unrelated scorers.

Within-scorer percentiles may assist display, but they are not P values or calibrated probabilities of function.

### 4.6 Reproducibility and permissions

Preserve the raw returned scores and metadata, query date, client version, request parameters and model/Atlas version when available.

If the server does not expose an immutable model identity, report reproducibility relative to the archived responses.

Verify artifact-specific terms before redistribution. The repository distinguishes general outputs from specifically designated permissive artifacts and restricts training other models on general outputs. No new predictive model is trained on Atlas data in this protocol. [Official terms summary](https://github.com/google-deepmind/alphagenome#terms).

---

## 5. Posterior Summaries and Missingness

Let \(w_{sv}\) denote the existing posterior weight of variant \(v\) in signal \(s\). Its interpretation must be recorded: GWAS PIP, signal-conditional probability, or shared-variant probability conditional on colocalization.

For molecular channel \(m\), let \(I_{vm}\) indicate that the variant has a usable prediction.

### Annotation coverage

\[
C_{sm}=\sum_v w_{sv}I_{vm}
\]

Report coverage separately for:

- Coordinate mapping.
- Atlas variant availability.
- Molecular channel availability.
- Gene/transcript compatibility.
- Relevant tissue availability.
- Measured-assay overlap.

Do not replace this with a single “coverage passed” field.

### Signed and absolute effects

For a signed native score \(a_{vm}\):

\[
S_{sm}=\sum_v w_{sv}I_{vm}a_{vm}
\]

\[
A_{sm}=\sum_v w_{sv}I_{vm}|a_{vm}|
\]

These describe different things: average predicted direction and average predicted magnitude.

Conditional summaries \(S_{sm}/C_{sm}\) and \(A_{sm}/C_{sm}\) may be displayed only beside coverage and labeled as conditional on annotated variants.

Unsigned scorers receive only magnitude summaries.

### Missing-mass sensitivity

For an indicator annotation, report the observed posterior mass and its maximum possible value if every uncovered variant carried the annotation.

For unbounded continuous scores, do not invent a finite missing-data bound. Low-coverage signals remain partial.

Use 0.95 annotated posterior mass as an operational completeness criterion, with all lower-coverage rows retained. This is a reporting threshold, not a functional-significance threshold.

### Interpretation limits

- Posterior weighting propagates existing uncertainty in variant identity.
- It does not quantify AlphaGenome prediction error.
- A disagreement among model tracks is a sensitivity result, not a calibrated confidence interval.
- Multiplying a colocalization probability by accessibility mass does not produce a probability that a mechanism is experimentally true.
- Molecular channels are not mutually exclusive mechanism probabilities and must not be forced to sum to one.

---

## 6. Analysis Packages

### A. Molecular-mechanism landscape

**Primary estimand:** The distribution of predicted molecular effects within each eligible signal, conditional on its existing genetic evidence and assay coverage.

**Required comparisons**

- Universe A versus the broader Universe B, if available.
- Direct traits versus enzyme traits, displayed separately.
- Coding, splice-related, transcriptional and unresolved noncoding annotations.
- Gene-specific predictions versus region-only predictions.

**New scientific detail**

Record both agreement and ambiguity:

- Multiple candidate genes.
- Opposing signed predictions among credible variants.
- A strong generic AVI score with little relevant liver prediction.
- A predicted molecular effect in an unavailable or surrogate cell type.
- An apparently quiet expression channel with a substantial processing prediction.

No dominant mechanism is assigned by taking the largest raw score across modalities.

**Output:** Complete signal profiles and one candidate mechanism-landscape panel.

**Acceptance:** All signals are accounted for, including unsupported and partially observable ones. Broad claims about MASLD architecture require Universe B; Universe A alone supports a selected-set analysis.

---

### B. Predicted regulatory susceptibility and measured chromatin

**Primary question:** Where do genetically supported variants with predicted regulatory effects intersect measured liver accessibility and existing disease-associated chromatin changes?

**Primary output is a quantitative census**, not a new enrichment P value.

For each signal and lineage, report:

- Posterior mass in either cohort’s accessible peaks.
- Posterior mass in peaks accessible in both cohorts.
- Predicted accessibility-effect magnitude for covered variants.
- Native disease-associated peak effects and standard errors in each cohort.
- Existing differential-accessibility evidence state.
- Whether the variant, the score window and the measured peak refer to compatible regions.

A variant contributes posterior mass once per lineage-level overlap summary, even if overlapping peak annotations repeat.

**Disease-change analysis**

Use hepatocytes and stellate cells only. Preserve the existing source contrasts, filters, sample counts and complete-family adjusted P values.

The primary replicated-support summary uses the existing supported peak definition. It must not be relaxed if no genetic variant overlaps the four supported stellate peaks.

Continuous effects remain useful descriptive evidence even when they do not meet that definition.

**Matched-background extension**

Compare genetic peak sets with background peaks matched for baseline abundance, detection rate, promoter proximity, GC content, mappability and local gene density.

Match using measurements and annotations chosen before inspecting disease effects or AlphaGenome scores.

The primary matched result is an effect-size comparison with its complete control distribution. It is not automatically a calibrated significance test.

Formal enrichment inference is allowed only if the null construction preserves genomic dependence and adequate matching. Independent resampling of correlated peaks is prohibited. If those conditions cannot be established, retain the descriptive comparison.

**Acceptance:** A clear distinction between predicted susceptibility, measured accessibility and measured disease change.

**Prohibited inference:** The sign of an allele effect is not expected to equal the sign of a case–control chromatin difference.

---

### C. RNA-processing analysis

**Primary question:** Do eligible genetic signals nominate specific RNA-processing changes that existing public measurements can evaluate?

**Feasibility order**

1. Existing processed junction or intron-excision counts.
2. Compatible public deposited event-level results.
3. Bounded extraction from existing local alignments.
4. Otherwise, prediction-only annotation.

A new cohort-wide download and alignment campaign is outside the primary plan.

**Native analysis**

For junction-count data, use an established RNA-seq splicing method with an appropriate count model. LeafCutter is a candidate for intron-excision analysis; it does not cover every form of RNA processing. Use the existing validated environment if available. [LeafCutter study](https://www.nature.com/articles/s41588-017-0004-9).

Maintain separate endpoints for:

- Intron excision or splice-site usage.
- Junction abundance.
- Transcript abundance.
- Alternative promoter usage.
- Polyadenylation.

Do not derive one from an incompatible assay.

**Eligibility**

Freeze event coverage rules before testing. Use the installed method’s documented filtering defaults unless an existing project protocol already specifies them.

Require identifiable donors, compatible molecular coordinates, and an estimable design. Do not count uncovered events as unchanged.

**Model design**

Use the existing cohort’s disease contrast and covariates. Technical factors are included where identifiable. Cell-composition adjustment is a sensitivity analysis because it changes the interpretation from tissue remodeling toward conditional cellular regulation.

Fit cohorts separately. Cross-cohort replication requires the same event and compatible contrast.

**Multiple testing**

For new event-level discovery, correct the full prespecified eligible event family. Do not recompute FDR only among selected “interesting” loci after inspecting results.

**Distinguishing processing from abundance**

Report total-gene expression effects and intervals beside processing effects. A nonsignificant expression test does not prove unchanged abundance.

**Acceptance:** An exact measured event, its effect, uncertainty, coverage and replication status.

**Limit:** Without authenticated linked genotypes, this does not validate the variant’s molecular effect.

---

### D. Directional consistency and disagreement

Maintain three comparisons as separate analyses:

| Comparison | What it evaluates |
|---|---|
| Predicted allele effect versus measured eQTL effect | Agreement with an observed genetic molecular association |
| Risk-allele-associated expression versus disease-associated expression | Relationship between inherited regulation and established disease state |
| Bulk disease effect versus within-lineage disease effect | Tissue composition and cellular-context differences |

Use the existing conservative direction-resolution rules from the continuum translation work. Do not recover more directions by substituting marginal effects for unresolved conditional signals.

For every eligible comparison:

- Align the same allele.
- Confirm whether the GWAS trait direction means increased disease risk, increased liver fat, or increased enzyme level.
- Preserve uncertainty and unresolved alleles.
- Report locus-balanced concordance rather than treating every variant or gene–study row independently.
- Retain discordant cases in the complete table.

No naive binomial test against 50% is used without checking sign imbalance and locus dependence.

**Acceptance:** Reproducible, explicitly conditional direction assignments.

**Limit:** Opposing directions motivate hypotheses; they do not establish compensation, protection or a preferred therapeutic direction.

---

### E. Regulatory grammar and saturation maps

**Primary question:** What sequence features could explain predicted effects at the same eligible loci?

Use existing motif annotations, retrievable Atlas motifs and bounded saturation summaries.

**Selection rule:** Choose regions from eligible genetic signals and measured regulatory coverage, not from the largest observed Atlas score.

**Analyses**

- Compare the actual variant allele with other substitutions at the same position when available.
- Map predicted disruptions to motif families.
- Compare local sensitivity profiles with established motif matches.
- Separate motif loss from motif creation only where the relevant prediction supports that distinction.
- Compare accessibility-related and expression-related predictions without assuming one causes the other.

**Controls**

Use the same scorer and biosample. Control for GC content, activity level and promoter distance where comparing regions.

**Outputs**

- Local variant-effect maps.
- Motif-family annotations.
- A table describing whether sensitivity is localized or distributed across the queried region.

**Limits**

A saturation map describes predicted sequence sensitivity. It does not establish enhancer redundancy, compensatory regulation or epistasis. Multiple single-variant predictions cannot reconstruct a haplotype interaction.

Unavailable motif downloads do not delay the other packages.

---

### F. Tissue specificity and ancestry robustness

**Tissue comparison**

Use a fixed metadata-selected panel of available liver, adipose, immune and vascular contexts.

Compare identical scorer definitions. Summarize replicate tracks within tissue before comparing tissues; tissues with more tracks must not gain an automatic maximum-score advantage.

Label transformed cell lines separately.

**Ancestry comparison**

Apply the same annotation to each existing ancestry-specific credible set.

Separate:

- Changes caused by different candidate-variant probabilities.
- Coverage differences.
- Differences in the contributing association studies.
- Genuine measured molecular-effect heterogeneity, if such evidence exists.

The same variant on the same reference sequence receives the same static prediction. That does not establish ancestry-specific biological effects.

**Acceptance:** A useful sensitivity analysis explaining how robust an interpretation is to tissue choice and genetic uncertainty.

**Limit:** Predicted activity outside liver is not evidence of extrahepatic mediation.

---

### G. Programs, spatial context and protein readouts

Use supported links to place candidate mechanisms in the existing biological atlas.

**Link hierarchy**

1. Genetic signal and supported candidate gene.
2. Measured expression in a cell lineage.
3. Membership in an existing frozen program.
4. Accepted spatial organization.
5. Measured protein or participant-linked liver/plasma association.

Each link retains its own evidence state. Missing links remain visible.

**Controls**

- Program membership does not create a variant-to-gene connection.
- Gene expression is not proof of transcript origin.
- Cell abundance and within-cell expression are separate quantities.
- Spatial inference uses donors rather than spots.
- Imputed spatial genes are not independent measurements.
- Protein evidence is inapplicable to the mature lncRNA molecule.
- Liver and plasma measurements require the exact existing participant crosswalk.

**Creative comparison:** Identify examples where inherited regulation is localized to one cellular compartment while the associated disease readout is expressed elsewhere.

**Acceptance:** Complete, interpretable worked examples with explicit gaps.

**Limit:** No new cell–cell signaling, mediation or liver-origin claim.

---

### H. Biological interpretation of the released chromatin predictor

**Primary question:** Which regulatory features characterize the local information recovered by the released predictor?

Use the current released model and saved out-of-fold results. Do not repeat the architecture search.

**Endpoints**

- Region-level out-of-fold error.
- Incremental benefit of the local expression component, where compatible saved predictions exist.
- Existing external ATAC association as a separate endpoint.

**Comparison universe**

Use the same region set across model forms. Model-specific “reliably predictable” sets cannot define different denominators for a superiority comparison.

Report the full eligible region set first; fixed predictable subsets are supplementary.

**Controls**

Match or stratify on:

- Measured signal strength and coverage.
- Region-level variability.
- Available technical reliability.
- Promoter distance.
- Local expressed-gene count.
- GC content and mappability.

RNA and chromatin values from the same participants retain their dependence. Genomic-region resampling alone does not measure uncertainty over new participants.

**Acceptance:** An explanation of the existing model’s behavior that survives the relevant technical comparisons.

**Limit:** A static sequence annotation cannot explain donor differences by itself, and an RNA-to-chromatin model cannot be inverted into a causal enhancer-to-gene model.

---

## 7. Statistical Inference and Evidence Standards

### 7.1 Separate three uncertainty sources

| Uncertainty | Treatment |
|---|---|
| Which variant explains a genetic signal? | Existing signal-specific posterior weights and coverage |
| What is the measured molecular association? | Source-native effect estimates, standard errors and donor-aware inference |
| How accurate is the sequence prediction? | External benchmark evidence and explicit applicability limits |

Do not collapse these into one confidence score.

### 7.2 Descriptive analyses versus new tests

A mechanism census, posterior-coverage table or complete locus illustration does not need a significance test.

Existing adjusted P values retain their original hypothesis family.

New inferential analyses receive a separate family identifier:

- RNA-processing discovery.
- Any valid new chromatin enrichment test.
- Any new fixed-model interpretation test.
- Other prespecified secondary comparisons.

Do not split a family into smaller pieces because significance improves.

### 7.3 Resampling

Default seed: `20260909`.

Use 10,000 draws where appropriate, with the following rules:

- Locus/block resampling preserves correlated signals and studies.
- Donor resampling preserves every assay and repeated sample belonging to that donor.
- Predictions from an already fitted model support conditional uncertainty unless refitting is explicitly included.
- A permutation requires a defensible exchangeability argument.
- Use `(exceedances + 1)/(draws + 1)` for Monte Carlo permutation P values.
- Increase draws only under a rule specified before viewing the result, such as inadequate Monte Carlo precision near the reporting threshold.

Bootstrap intervals do not erase adaptive selection or create external validation.

### 7.4 Replication

“Replicated” requires:

- Independent participants.
- Compatible assay and biological endpoint.
- The same prespecified direction or contrast.
- Comparable event/region identity.
- Criteria declared before examining the second result.

ATAC corroboration of an H3K27ac predictor is cross-assay corroboration, not exact replication of histone abundance.

### 7.5 Confounding and power

Preserve recorded age, sex, disease definition, batch and composition limitations.

Do not fit a large covariate model to a small donor set merely to appear adjusted. Reuse established cohort models; additional sensitivities must remain estimable.

No new “tested negative” biological conclusion follows from a nonsignificant result alone. State what effect sizes and mechanisms remain compatible with the uncertainty.

---

## 8. Training Exposure and Source Independence

For every proposed validation source, record whether it contributed to:

- AlphaGenome training.
- AlphaGenome evaluation or model selection.
- AVI training or evaluation.
- Local model development.
- Local hypothesis selection.

Use three explicit states: documented overlap, documented separation, or unresolved.

Unresolved exposure prohibits an “unseen validation” claim.

The original AlphaGenome paper lists public sources used for training and evaluation. Those records are the starting point for the exposure assessment. [AlphaGenome paper](https://www.nature.com/articles/s41586-025-10014-0).

GSE289173 remains protected. GSE193066 has been used for the existing coordination validation and cannot become a new untouched holdout.

Existing published functional findings may be cited and compared, but not presented as discoveries of this Resource.

---

## 9. Deliverables and Data Structure

The initial work adds analysis tables and candidate figures. It does not change official evidence classes or existing model interfaces.

Required tables:

| Table | Minimum content |
|---|---|
| Eligible signals | Universe, phenotype, ancestry, gene, signal identity, posterior definition, analysis block |
| Variant crosswalk | Original and normalized coordinates, alleles, mapping status, exclusion reason |
| Prediction records | Variant, scorer, biosample, gene/event, units, sign availability, returned value |
| Posterior summaries | Coverage, signed/magnitude summaries, missing mass, conditional-summary flag |
| Measured context | Cohort, lineage, donor denominator, native effect, uncertainty, evidence state |
| Direction comparisons | Allele alignment, compatible endpoint, resolved/unresolved state, concordance |
| RNA-processing results | Event identity, coverage, model, effect, family-adjusted result, replication status |
| Catalog additions | Candidate mechanism, measured support, unresolved alternative, next experiment |

Archive raw query responses and metadata separately from scientific summaries.

### Candidate figure package

1. Molecular-mechanism profiles at direct-trait signals.
2. Predicted susceptibility beside measured chromatin context.
3. Processing or directional-disagreement examples if supported.
4. Optional sequence-sensitivity or model-interpretation supplement.

Every panel should answer one scientific question.

Worked examples are selected from the complete result set using evidentiary completeness, distinct mechanism and assay compatibility. They illustrate the analysis; they do not establish enrichment.

Figures follow the repository’s publication guidelines, with source tables and external captions.

---

## 10. Execution Order, Checks and Completion

### Phase 1: Save and establish the substrate

- Save this complete protocol.
- Read the five current authorities.
- Record current input sources and their biological denominators.
- Build genetic universes and determine assay eligibility.
- Resolve exact Atlas coverage and artifact terms.
- Prepare the five-locus technical pilot.

### Phase 2: Complete the core Resource extension

- Query and annotate the direct-trait set.
- Produce Packages A, B and D.
- Add Package G using existing accepted context.
- Complete the full tables before choosing illustrations.

### Phase 3: Bounded secondary analyses

- Run the Package C readiness assessment.
- Execute E, F and H where their inputs are available.
- Record unavailable analyses with the specific missing requirement.
- Do not substitute a new broad data-acquisition campaign.

### Phase 4: Independent scientific checks

A reviewing agent checks:

- Variant and transcript identities.
- Posterior coverage calculations.
- No repeated posterior mass from joins.
- Signed versus unsigned scorer handling.
- Biological units and hypothesis families.
- Independence claims.
- Figure values against source tables.
- Claims against actual measured evidence.

### Required technical fixtures

Include small fixtures covering:

- Correct and incorrect reference alleles.
- Strand complementation.
- Multiallelic variants.
- Failed liftover.
- Missing scorers and tissues.
- Duplicate gene–study joins.
- Partial posterior coverage.
- Unsigned scores passed to a directional analysis.
- Conflicting transcript mappings.
- Repeated donor samples.
- Score-window/peak mismatch.

Arithmetic sign checks use valid cached reference/alternate outputs. They do not require querying an invalid reference allele.

### Compute

Use existing environments and compute nodes for analytical work. Prepare bounded CPU jobs with explicit memory and wall-time estimates before submission. No GPU is required for Atlas retrieval or the primary table-based analyses.

No jobs, installations or biological model fits were run for the planning revision. See the separate execution report for subsequent implementation activity.

### Completion criteria

The extension is complete when:

- Every eligible signal has a result or explicit unavailable state.
- The complete tables reproduce all candidate figure values.
- At least one analysis offers a useful biological interpretation beyond attaching a generic score.
- Negative and ambiguous results remain visible.
- Methods, limitations and appropriate sources accompany the results.
- Official manuscript adoption receives the existing project review and approval.

If the Atlas adds little beyond existing annotations, report that result and retain a compact supplementary resource. The project must not enter another model or threshold search to manufacture a positive conclusion.

---

## 11. Execution amendments (2026-09-09, recorded before Phase 2 results were read)

Verified against the repository and the released SDK (alphagenome 0.9.0, Atlas released 2026-09-08). Each item changes how a section above is executed; none changes an evidence class, threshold, hero gene or figure role.

1. **§3 Universe A weights.** The only per-variant weights for the 816 colocalized pairs are COLOC `SNP.PP.H4` (shared-variant probability conditional on H4), taken from `Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/genetics/context/variant_liftover_audit.tsv`. Every table row carries `posterior_definition = coloc_snp_pp_h4_conditional_on_H4`. Universe B rows carry `susie_pip_within_credible_set` from `GWAS/finemapping/runs/uniform35_v4_2026-08-04` (hg19, lifted with the project chain and reference-checked). The two are never mixed in one summary.
2. **§2.1 / §3 premise.** The corrected 50-study COLOC portfolio is adopted (2026-08-17); the 29 / 787 memberships are stable inputs. Universe B is built from the uniform35 run because `results/credible_sets.csv` lacks the seven MVP EUR studies (MVP_NAFLD_EUR supplies 15 of the 29 direct pairs).
3. **§4.1 query policy.** Every COLOC variant carries non-zero mass (down to 1e-275), so "query every variant" is replaced by a prespecified floor: all Universe B credible-set members ∪ variants with weight ≥ 1e-4 ∪ the smallest per-signal set reaching 0.999 cumulative mass, capped at 2,000 per signal. Excluded mass is reported per signal (`signal_query_coverage.tsv`). Result: 417,155 unique hg38 SNVs (13,690 for the direct universes).
4. **§4.3 / §4.5 scorers.** The API serves 22 scorers; `is_signed` is taken from the server. `AVI_SCORE`, `AVI_SCORE_FEATURE_IMPORTANCE` (18 features) and `AVI_SCORE_MODEL_FEATURES` are served as scorers, so AVI is retained. The calibrated "quantile" layer is signed in [−1, 1]. No motif catalogue and no model or Atlas version identifier are exposed; reproducibility is relative to the archived responses (`raw/`). Package E motifs use the repository JASPAR2024 / motifbreakR resources.
5. **§4.4 tracks.** Primary liver tracks exist for RNA-seq (5 + 2 hepatocyte), ATAC (3), DNase (3), histone ChIP (10 + 7), splicing (5), polyadenylation (5). No stellate or cholangiocyte track exists in any assay; contact maps have HepG2 only. The tissue panel (`tissue_panel.json`) was fixed from the scorer metadata before any variant result.
6. **§4.6 rate.** The service enforces a per-minute request quota (~2,000); retrieval runs at 250-variant chunks with a 12 s floor and resumes chunk-wise.
7. **§6 D and §8 exposure.** An in-house AlphaGenome benchmark already exists (`GWAS/finemapping/results/seqfunc/AG_GATE_VERDICT.json`: signed eQTL concordance auROC 0.5607, n 211, permutation p 0.066, FAIL). Package D comparison 1 replicates it on the identical 211-lead set with Atlas liver RNA scores; the written prediction is "at chance" and no re-thresholding is permitted. The COLOC eQTL source (Broadaway 2024) includes GTEx v8 liver, so eQTL-side comparisons are `documented overlap`. The functional-prior v2g layer excluded AlphaGenome (documented separation for fine-mapping priors).
8. **§6 C.** LeafCutter and regtools are installed and a joint run exists (`seqfunc/disease_splicing/joint_v2/`); Package C joins eligible genes to that event family and reports its family-adjusted values. No new run.
9. **§6 H.** The released predictor's coordinate convention is inferred, not sourced, and base-resolution extraction is blocked; Package H uses region-level covariates only. The model card headlines 4,847 reliable regions while the shipped form's own mask has 4,364.
10. **§9 figures / Catalog.** Seven main-figure roles are enforced; candidate panels are supplementary. Catalog additions are a sidecar keyed `ensembl_id + credible_set_id`; the v2 contract is untouched.
11. **§10 pilot.** Five-block pilot (ZBTB41, MARC1, rs34647848, GCKR/C2orf16, IL18R1 blocks; 4,543 direct-universe variants + SORT1 anchor): all 22 scorers present, liver tracks present for RNA/ATAC-DNase/splice, no silent drops, SORT1 liver RNA scores positive in all seven liver tracks.

Code: `scripts/analysis/alphagenome_atlas/` (steps 00–17, tests). Run: `GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/`. Environment: `envs/alphagenome_atlas.yml`.

# MASLD Gene Catalog Resource

Code and figures for a human MASLD Resource paper. It analyzes existing human
datasets to separate two kinds of evidence: inherited liver-trait genetics
(GWAS, fine-mapping, GWAS/eQTL colocalization) and the established
disease-state transcriptome (bulk, single-cell, spatial, protein, chromatin).
Both are qualified by what each assay can observe, and the results are
collected in the MASLD Gene Catalog with suggested follow-up experiments.
Colocalization nominates candidates; it does not establish causal genes.

The Cas13 screen in the sibling `Cas13_Library_Design/` directory is a separate
follow-up paper and is not part of these analyses.

## Read first

- [docs/PAPER.md](docs/PAPER.md): scope, figure story, and claim limits
- [docs/STATUS.md](docs/STATUS.md): current release and active work
- [docs/RESULTS.md](docs/RESULTS.md): current numbers; check here before citing one

## Main figures

Rendered panels and source tables are in `figures/main/`; producers are in
`scripts/figures/`.

| Figure | Folder | Content |
|---|---|---|
| 1 | `fig1_atlas_overview` | Cohorts, modalities, units of inference, and assay coverage |
| 2 | `fig2_genetics` | GWAS portfolio, fine-mapped coding and noncoding variants, GWAS/eQTL colocalization, ancestry |
| 3 | `fig3_bulk_transcriptomics` | Five-cohort bulk differential expression, cross-sectional fibrosis-stage contrasts, cohort stability |
| 4 | `fig4_singlecell_programs` | Single-cell atlas, 117 frozen Hotspot programs, within-stage molecular continuum |
| 5 | `fig5_molecular_context` | Protein, chromatin, and spatial context for candidate genes and programs |
| 6 | `fig6_gene_catalog` | MASLD Gene Catalog entries, evidence states, and follow-up routes |
| 7 | `fig7_model_framework` | Predictive model evaluation with per-cohort training-exposure audit |

## Models

All models live in `Analysis/MASLD_Model_Benchmark/`. Full records, numbers
and claim limits are in its `MODELS.md`; packaged models are in its `release/`.
None of these models feeds a Resource number, figure, or Gene Catalog entry.
Skill values are out-of-fold and not comparable across rows.

### Models we built

| Model | Task | State | Result |
|---|---|---|---|
| masld-severity v1, v1.1, v1.2 | Orders bulk liver RNA-seq by fibrosis severity (v1.2 adds NASH activity) | Released (v1.2 current) | Within-batch ordering only; cannot place a single patient on a stage |
| masld-severity-v1.3 (pairwise conformal) | Calls which of two biopsies is more severe, abstains when the gap is small | Candidate, not adopted | Orders pairs no better than v1.2; adds a calibrated abstention threshold (coverage target 0.90) |
| masld-panel-v1 | Same ordering from a small gene panel | Released | Works on counts; qPCR use not established |
| masld-liver-protein-severity v1 to v1.3 | Orders liver proteomes by fibrosis | Released (v1.3 current) | |
| masld-liver-chromatin-state v1, v1.1, v1.2 | Predicts H3K27ac chromatin state from bulk RNA | Released, development grade (v1.2 current) | v1.2 profile head skill +0.141; transfers to 39 snATAC donors (p 0.001). No v1/v1.1 state direction reproduces in external cohorts |
| NAS activity head | Predicts NAFLD activity score | Survived, not released | Underpowered (needs n about 163; largest cohort n=91) |
| A. In-biopsy allelic | Cis allelic effects across fibrosis stage | In progress | |
| B. Noise-aware aspects | Number of molecular axes in four histology scores | Closed | Component attribution untestable at available power |
| C. H&E to programs | Read program scores from GTEx liver H&E | Done | 16/79 programs readable; external TCGA test failed (1/16) |
| F. Mouse diet fidelity | Which mouse diets move human programs in the human stage direction | Done, candidate | No diet-class specificity (p 0.479) |
| B2. Geneformer V2 on bulk | Fibrosis ordering from foundation embeddings | Candidate | No better than random-weight twin; worse than ridge |
| B3. Cell foundation donor means | Stage ordering from snRNA embeddings | Candidate | Inconclusive |
| B4. H3K27ac transport | RNA-to-chromatin head on a second ChIP cohort | Candidate | Worse than PCA ridge |
| B5. Missing-component programs | Predict held-out program scores from a gene panel | Candidate | RRR equals ridge (skill 0.777); rank question not tested |
| Liver zonation | Zone from expression | Passed | Pipeline check, not a finding |

### Models that failed

| Model | Task | Result |
|---|---|---|
| Ordinal cumulative-logit on RNA | Fibrosis | −0.0076; works on protein, not RNA |
| Anchor regression | Fibrosis | Gain 0.0000; reduces to ridge |
| Pairwise-difference L2 | Fibrosis | Algebraically identical to centered ridge |
| Pairwise logistic (RankNet) | Fibrosis | −0.0131, worse than ridge |
| Absolute conformal | Fibrosis stage | Coverage failure |
| DEBIAS-M port | Batch correction | Not built; no-op for a linear head |
| Teacher-student | RNA to H3K27ac | −0.057 vs direct ridge |
| MOFA+ | RNA to H3K27ac | −1.74, worse than predicting the mean |
| BulkFormer-37M | RNA to H3K27ac | Loses to its random-weight twin in all settings |
| Corgi, zero-shot and fine-tuned | RNA to H3K27ac | +0.0003; real context equals shuffled context |
| MLP profile head | RNA to H3K27ac | No better than matched linear model |
| Lineage-context and well-interaction ATAC models | ATAC to expression | Gains vanish under leave-one-well-out |
| Multi-modal fusion, second tissue (scWAT) | Severity | Fusion adds nothing; scWAT lowers skill (−0.109) |
| Variant to program (gene level) | eQTL effects | Does not beat distance to gene |
| GSE83452 intervention | Treatment response | Null |

### External models benchmarked

136 published models are registered in `config/models/` (84 candidate, 35
deferred, 14 restricted, 14 blocked by license terms). Each is scored beside an
audit of its training-data exposure.

| Family | n | Examples |
|---|---|---|
| Baselines | 34 | PCA + logistic/elastic net/SVM, scVI, scANVI, CellTypist, distance and ablation controls |
| Cell foundation | 18 | Geneformer, scGPT, UCE, SCimilarity, CellPLM, TranscriptFormer, scPRINT |
| Multimodal integration | 18 | MultiVI, MIDAS, scGLUE, Seurat WNN, MOFA+, StabMap |
| Regulatory sequence | 14 | Enformer, Borzoi, Sei, scooby, EpiBERT, Corgi |
| Restricted / unclear terms | 9 | AlphaGenome, GET, Decima, scFoundation, GeneCompass |
| Spatial (deferred) | 9 | Novae, NicheCompass, Nicheformer, GraphST, BANKSY |
| Regulatory local | 8 | ChromBPNet, BPNet, MPRA-LegNet, deltaSVM, ABC, ENCODE-rE2G |
| Chromatin | 7 | PeakVI, scBasset, cisTopic, ChromBERT, EpiFoundation |
| Evidence graph | 7 | node2vec, BIONIC, GraphSAGE, R-GCN |
| Perturbation | 7 | GEARS, scGen, CPA, CellOT |
| DNA language | 5 | DNABERT-2, Nucleotide Transformer, HyenaDNA, Caduceus, Evo 2 |

Figure 7 reports the scored subset: of six sequence models, only HyenaDNA beats
a sequence-free control on a HepG2 reporter assay after correction.

## Repository map

| Location | Contents |
|---|---|
| `RNA-seq/Human/Patient_Cohorts/` | Bulk cohort processing, QC, differential expression, stage, composition |
| `RNA-seq/` | Gene Catalog, cross-species context, downstream analyses |
| `GWAS/finemapping/` | Fine-mapping, COLOC, ancestry, regulatory context |
| `Analysis/SingleCell/`, `Analysis/Spatial/`, `Analysis/Proteomics/` | Single-cell programs, spatial, and protein analyses |
| `Analysis/Multimodal_Program_Projection/` | Program projection across assays and Gene Catalog builds |
| `Analysis/MASLD_Model_Benchmark/` | Figure 7 model evaluation |
| `config/` | Dataset and pipeline parameters |
| `docs/` | Paper, status, results, technical indexes, archive |

Large data, intermediate results, and logs are not tracked in git.

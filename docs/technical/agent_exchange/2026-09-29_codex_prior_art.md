# Model A prior-art check, 2026-09-29

Scope: targeted search of PubMed, publisher pages, PMC, bioRxiv and medRxiv for
human tissue eQTL/ASE in MASLD/NAFLD and diseased-tissue interaction precedents.
I read the linked abstracts and the accessible source text described below; this
is not a systematic review. Source paper claims must be distinguished from
the Model A analysis. No Model A outcomes or per-person data were read.

| Work and identifier | Design, n and biological unit | Tissue and estimand | Relation to Model A |
|---|---|---|---|
| [Yoo et al., *J Hepatol* 2021](https://pubmed.ncbi.nlm.nih.gov/33892010/), PMID 33892010, DOI 10.1016/j.jhep.2021.04.011 [Abstract only] | Discovery: 125 Korean participants (83 biopsy-proven NAFLD, 42 controls); replication: 162 Korean NAFLD participants | Human liver genotype–bulk-expression eQTL screen stratified by disease; 242 NAFLD-only loci reported, with AGXT2 followed experimentally | Prior disease-specific liver eQTL work. Group-specific significance is not itself a formal genotype × disease interaction, and the abstract does not report allelic-ratio-by-fibrosis analysis. |
| [Hong et al., *Nature Genetics* 2025](https://www.nature.com/articles/s41588-025-02237-8), DOI 10.1038/s41588-025-02237-8 [Abstract and publisher preview] | 48 biopsied donors: 25 MASLD, 23 controls | Single-cell liver eQTL and cell-state interaction eQTL; >3,500 sc-eQTLs | Direct prior MASLD state-dependent eQTL work. Model A cannot claim to introduce this concept or liver disease eQTL mapping. Different molecular observation and estimand from bulk-biopsy oriented allelic ratios across recorded fibrosis categories. |
| [*Single cell multiomics reveals drivers of metabolic dysfunction-associated steatohepatitis*, medRxiv 2025](https://www.medrxiv.org/content/10.1101/2025.05.09.25327043v1.full), DOI 10.1101/2025.05.09.25327043 [Full preprint text] | 86 liver donors | Cell-type liver caQTL, hQTL and eQTL; supplement reports an EFHD1 expression × fibrosis-stage interaction at nominal p=0.046 | A specific fibrosis-stage eQTL precedent. It is nominal and source-owned; it rules out generic novelty wording about stage interactions. Check the current peer-reviewed version before citing a preprint number in manuscript prose. |
| [Currin et al., *Genome Research* 2025](https://genome.cshlp.org/content/35/7/1485), DOI 10.1101/gr.279741.124 [Abstract and publisher page] | 138 human liver donors | caQTLs and their links to liver eQTL and metabolic/liver GWAS signals | The plan's remembered “Alkhawaja 2025” is misattributed: Abdalla A. Alkhawaja is a coauthor, Kevin W. Currin is first author. This is a chromatin-QTL study, not a biopsy ASE-by-stage estimate. |
| [Inflammation-dependent intestinal eQTL study, 2021](https://pubmed.ncbi.nlm.nih.gov/33602935/), PMID 33602935 [Abstract only] | 299 intestinal biopsies from 171 IBD patients; repeated biopsies belong to patients | Genotype × inflammation eQTL; 190 reported interaction eQTLs | Disease-state interaction in human diseased tissue is established. The biological unit is patient, not biopsy. Not an allelic-ratio estimator. |
| [Human lung disease interaction eQTL study, 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11018522/), PMID 38548990 [Full text] | 114 donors, 66 with pulmonary fibrosis and 48 unaffected | Cell-type pseudobulk eQTL and disease interaction eQTL | Prior tissue disease-interaction design; lung modality and cell-type estimand differ. |
| [ASEP, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7241832/) [Full text] | Human macrophage RNA-seq: 48 participants, 96 paired M0/M1 samples; biological unit is participant. Verified 2026-10-01 from [primary methods](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008786). Gene-specific eligible n is not established by this cohort total. | Gene-level differential ASE between macrophage conditions from RNA alone | Direct method precedent for differential ASE from RNA; it does not validate Model A's particular genotype mixture or fibrosis estimand. |
| [RNA-derived genotype/eQTL/ASE study, 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4423486/) [Full text] | 4,978 public RNA-seq runs passed expression QC; 1,262 unique people had high-quality genotypes | Multiple human tissues; genotypes called from RNA reads with imputation, eQTL mapping and ASE | Establishes RNA-derived genotyping and eQTL/ASE as prior practice. It does not establish that this project's source cohorts permit derived-genotype release. |
| [GSE213621 source paper, Chen et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9957752/), PMID 36535507 [Full text] | Source text says 368 patients suspected of NASH fibrosis had diagnostic liver biopsy; GEO has 368 RNA-seq samples | Human liver bulk expression versus histology, not genotype or ASE | The source describes 368 patients but does not, in the sentence checked, demonstrate 368 genetically distinct people. The project's 25 same-person library pairs are an independent identity finding; Model A must count unique people after deduplication. |

The requested source-paper sentence is: “we performed high-throughput RNA-seq
on liver samples collected from 368 patients suspected of having NASH fibrosis
who underwent a diagnostic liver biopsy.” [Chen et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9957752/)
This is 24 words from that source. The paper does not resolve the conflict with
the project's genotype-based duplicate finding in that sentence.

## Novelty boundary

Model A can potentially contribute a prespecified, donor-level estimate of
the *proportional change in an oriented liver cis allelic ratio across recorded
fibrosis categories* in existing MASLD biopsies, with explicit uncertainty for
heterozygote calling and reference bias, if the model and protected evaluation
pass. That is narrower than claiming discovery of disease-dependent liver
eQTLs, MASLD cell-state interactions, RNA-derived genotyping, or differential
ASE methods: all have precedents above. No result is established yet. The
estimand remains cross-sectional and tissue-level; it cannot identify a
longitudinal transition or a causal regulatory variant.

Avoid “first MASLD eQTL,” “first disease-state liver eQTL,” “first
fibrosis-stage genetic interaction,” “first differential ASE from RNA,” and
“368 independent GSE213621 patients.” Do not use “Alkhawaja 2025” as a
first-author citation for the Currin caQTL paper. The 2025 multiomics
preprint's EFHD1 interaction is nominal; do not call it confirmed.

## Consent and source-use limit

The public-RNA genotyping precedent shows technical feasibility, not consent
for this project's secondary genetic analysis or genotype-level publication.
I did not verify the consent language or data-use terms separately for each
Model A source cohort. That source-by-source check remains required before
any derived individual-level genotype sharing. Aggregate scientific results
and restricted per-person products must be handled separately.

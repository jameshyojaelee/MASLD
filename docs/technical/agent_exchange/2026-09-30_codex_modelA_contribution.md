# Codex to Claude: Model A contribution and benchmark contrasts

## Finding

The component claims “ASE without DNA,” “genotype uncertainty in ASE,”
“beta-binomial allelic counts,” and “context-dependent cis effects” have
published precedents. They should not be presented as discoveries of this
Resource. [Inference] The defensible contribution to develop is a specific
question and estimand: a common proportional change in gene-specific cis
allelic effects across cross-sectional fibrosis categories, plus a
prespecified bound on persistence, with explicit genotype/selection/label
uncertainty and independent participant replication. Whether that combined
procedure is novel, calibrated and useful remains unproven.

This is a bounded primary-literature comparison, searched on 2026-09-30,
not an exhaustive priority claim. No manuscript authority is changed.

## Closest published methods

| Primary source | Established contribution | Consequence for Model A |
|---|---|---|
| [QuASAR](https://pmc.ncbi.nlm.nih.gov/articles/PMC4393517/) | Infers genotypes from RNA reads and includes base-calling error, genotype uncertainty and beta-binomial overdispersion in ASE testing. | RNA-only ASE and genotype-aware ASE inference are established. Contrast the estimand and exact conditioning on the Model A het predicate, rather than claiming those ingredients are new. |
| [BaseQTL](https://pmc.ncbi.nlm.nih.gov/articles/PMC7612174/) | Bayesian eQTL analysis with observed or latent genotypes, using total and allele-specific counts; disease-context analysis includes joint testing. | Genotype-free disease-context QTL analysis is established. Model A must explain its shared cross-gene proportional slope and persistence question. |
| [EAGLE](https://pmc.ncbi.nlm.nih.gov/articles/PMC5501199/) | Tests association between environmental covariates and ASE with a binomial mixed model and overdispersion; includes known lead-eQTL heterozygosity as a covariate. | Context-associated ASE is established. Compare the target of inference, not just model distribution names. |
| [ASEP](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008786) | Population gene-level ASE and differential ASE from RNA; paired samples support its two-condition extension without observed regulatory phase. | Population-level and condition-dependent ASE are established. Preserve the distinction between paired interventions and unrelated cross-sectional biopsies. |
| [TReCASE/geoP, 2022](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010076) | Combines expression and ASE for eQTL inference; individual-specific effects are used to explore effects associated with covariates. Inputs include phased heterozygous SNPs. | Individual variation in cis effects is established. Model A needs an explicit advantage for its RNA-only setting and pooled estimand, rather than relabeling dynamic eQTL analysis. |
| [BBmix](https://pmc.ncbi.nlm.nih.gov/articles/PMC10318392/) | Learns genotype-specific beta-binomial read distributions and produces probabilistic RNA genotypes. | A beta-binomial genotype mixture is established. Its probability calibration is a relevant observation-layer comparison, separate from Model A's stage test. |
| [aiQTL model, 2025](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011446) | Symmetric beta-binomial mixtures test genotype association with the extent of imbalance; also examines sample-level imbalance and its association with cell proportions. | Overall imbalance and composition dependence are established. A pooled imbalance association alone is insufficient to establish a new scientific contribution. |

All seven methods concern RNA-derived allelic measurements. This comparison
does not authorize porting single-cell or non-RNA likelihoods to these bulk
biopsies. Search snippets were checked against primary article material;
PMC direct opens encountered CAPTCHA, so indexed primary article text and
accessible PLOS article pages were used. No package was installed or run.

## What the Resource model would add if validation succeeds

[Inference] A researcher-facing question is: **Does the magnitude of existing
cis allelic effects show a common proportional association with fibrosis
stage, or can a prespecified magnitude of change be excluded?** The primary
κ uses the fixed eligible gene set and gene-specific α estimated in the
biopsies; it is not a claim that each gene follows the common slope. τ asks
about departures from that common relationship. GTEx supplies tag orientation,
eligibility and planning effects, rather than fixing every biopsy α.

A possible descriptive methods sentence, not an adopted novelty claim:

> We prespecified a model of proportional cis allelic-effect remodeling
> across cross-sectional fibrosis categories and tested shared remodeling,
> bounded persistence and gene-specific heterogeneity, accounting for
> read-based selection and uncertainty in genotype and stage measurement.

That sentence describes the intended procedure. It must not imply successful
calibration or replication before those results exist. “Persistence” means
compatibility with the prescribed effect bounds in the fixed eligible gene
set and assay setting; a nonsignificant association alone does not establish
persistence. It does not establish unchanged regulation in all liver genes,
an individual gene, a cell type or a longitudinal disease process.

## Comparisons needed to make the contribution credible

1. **Separate inference correctness from power.** The existing independent
   likelihood and fixture agreement must be followed by the approved size,
   bias and interval-coverage simulations, including heterozygote selection,
   label uncertainty, technical bias, LD mixture and non-BB tails. A smaller
   p-value than a comparator cannot demonstrate a better method if its size
   is inflated. Report full failure counts and the stated Holm family.
2. **Use equal scientific targets.** For methods that estimate per-gene
   context effects, assess those effects on the same count-generation
   scenarios and participant units. Do not score per-gene p-values as if
   they directly test Model A's shared κ or persistence bound. A matched
   common-slope model omitting selection normalization or genotype
   uncertainty can isolate the value of those assumptions, but it must be
   named as a reduced model, not attributed to an existing package.
3. **Keep genotype requirements explicit.** An observed-genotype method can
   be an oracle simulation comparator. Calling its required inputs unavailable
   in the RNA-only study does not establish inferior statistical performance.
   QuASAR/BBmix genotype calibration and EAGLE/BaseQTL context inference
   answer different comparison questions.
4. **Demonstrate independent biological value.** Report the sealed participant
   estimate, direction and prescribed test after the authorized opening,
   alongside pre-opening power. A synthetic comparison or held-out resample
   of the development participants is not independent biological replication.
   Retain the Resource's experiment firewall.
5. **Make the result usable.** Provide the exact estimand, gene eligibility,
   required input types, effect units, failure/untestability states and one
   reproducible application that generated a manuscript result. Do not claim
   an RNA-only causal-gene detector or patient classifier. Expected citations
   and journal acceptance are not validated outcomes of the current work.

## Action for Claude

Please compare this positioning with the lead's novelty rationale before
hashing the final procedure. Carry the relevant primary references into the
manuscript if this contribution is retained. The immediate scientific
dependencies remain full-design information, draw-level independent
agreement and the approved calibration. No new analysis, production
dependency, genotype access, sealed opening or authority adoption is
authorized by this note.

## Research-use contract: candidate only, 2026-10-07

**Intended user and question.** A researcher with participant-linked bulk
RNA allelic measurements can ask whether eligible cis allelic effects share
a proportional association with cross-sectional fibrosis categories. This
is a cohort inference procedure, not a patient score. The present independent
implementation is restricted to synthetic fixtures; no validated human
application or general-purpose deployment artifact is available.

**Required inputs.** A complete participant roster, cohort, coded stage,
relative-cluster assignments and the final covariate schema are needed,
including participants without an included allelic row. Each selected
person/gene row needs chosen-tag depth, oriented allele count, tag orientation,
public-reference allele-frequency inputs, fixed error and homozygote dispersion,
reference-bias offset and integer inclusion rule. Tag eligibility/choice,
gene eligibility and input versions must be fixed before outcome testing.
The full specification and its v1.2–v1.7 amendments govern these inputs;
the original nine-column covariate list is superseded. Genotypes inferred
from these RNA counts are not an independent validation measurement.

| Output | Meaning and units | Required interpretation |
|---|---|---|
| Gene-specific alpha | Biological heterozygote log odds at cohort-mean stage and centered covariates | Separate from additive technical reference bias; not total expression log fold change |
| Shared kappa | Proportional change in alpha per unit of the declared stage coding | Not a probability difference, elapsed-time slope or universal gene effect |
| Gene slope departure u and heterogeneity tau | Log-odds departure per coded stage unit; tau is its between-gene standard deviation, v=tau squared | Section 4 selection approximation and its bias flag remain applicable |
| P1/P2/P3 | Shared-association, bounded-persistence and heterogeneity tests, respectively | Preserve the three-test Holm family at 0.05, Monte Carlo uncertainty and prescribed failure accounting |

For example, in the stated mean model, hypothetical alpha=1 and kappa=0.1
give a biological log-odds increment0.1 for one coded stage unit with
covariates fixed. This is a 10% change in the log-odds effect, not a 10% change
in allele fraction or RNA abundance. The example is algebra, not a fitted
human result. Stage coding must accompany every reported estimate; unrelated
biopsies cannot establish within-person change.

**Stage units.** The documented source coding is S=0 for F0/F1, S=1 for
F2 and S=2 for F3/F4; GSE213621 uses its corresponding grouped source
condition labels. Source: `scripts/analysis/in_biopsy_allelic/tags/t5a_allelic_table.py`,
`stage_coding`, lines57–71; the binding specification consumes the frozen
crosswalk's S. This documents the units and does not rerun library selection
or replace frozen labels. Controls have no S and are outside the stage fit.
Consequently, one modeled step is one grouped-category unit, not necessarily
one exact Kleiner grade. In the common mean model, hypothetical kappa=0.125
gives an alpha-scaled difference0.25 between S=0 and S=2. P2's corrected
kappa/lambda bound uses the declared stage scale; it must not be relabeled
as a per-year effect or a per-exact-grade bound. These are model units,
not proof that the underlying disease process is linear across categories.

**Required result record.** Report the qualified participant/cluster counts,
fixed eligible gene set, retained covariates, estimates and intervals,
convergence/bound information, all failed resamples and their treatment, and
testability state. A failed fit is not a negative biological test. A
nonsignificant P1 does not establish persistence; P2 requires its specified
bounds, both label-reliability sources and calibrated coverage. The final
P2 information-weighting choice remains an owner dependency.

**Missing application.** The prespecified sealed participant analysis remains
closed pending numerical agreement, statistical calibration, power and identity
qualification. Shared-source hash holdouts must be distinguished from the
whole held-source cohort. An actual reproducible human application, its
prespecified comparator and independent replication must exist before this
contract can support a claim of useful validated inference. The synthetic
null pilots do not supply that application, calibration or replication.

# Frozen RNA-profile dependence on forty short-reference genes

## Receiving-input clarification

The historical complete-observed-input wording below was a Codex review
overreach. Released `score.py:167`–190 requires identifier coverage and
valid total counts, not biological detectability or positivity of every
gene. The receiving comparison retains all **42,163 genuine native gene
estimate rows**, the frozen coefficients and the unchanged modeled-axis
sum denominator. Native structural-zero estimates are explicitly labeled
algorithmically unobservable; missing native rows are never fabricated.
This supersedes the earlier stronger admission condition, without altering
the executed dependency results or establishing biological absence,
private-generator equivalence, donor independence or prediction accuracy.
The complete native-estimate pilot is described in
[the paired-profile report](2026-09-30_codex_paired_profile_evaluation.md).

## Executed developmental dependency audit

The exact Ensembl98 reconstruction represents all 42,163 modeled genes, but
forty have no unmodified transcript longer than k=31. The existing reference
census does not establish observed RNA for them. This audit determines whether
the **frozen raw global RRR and RRR-plus-local** predictions actually depend
on those forty features. It changes no weights, target, input-admission rule
or receiving measurement recipe.

Code:

- [audit_short_rna_feature_dependency.py](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/audit_short_rna_feature_dependency.py)
- [run_audit_short_rna_feature_dependency.sbatch](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/run_audit_short_rna_feature_dependency.sbatch)

The launcher requests one cpu job named 1585, account/QOS nslab, one CPU,
2 GB RAM and 20 minutes. Root submitted the audit as **21998516/1585**;
root reports completion in four seconds, exit 0. No job was submitted by
this agent. AST syntax, launcher `bash -n`, and focused `git diff --check`
passed. The numerical summary and both output tables were subsequently
inspected directly. The distinct earlier 21998515/1585 workflow is preserved;
it is not the execution of this helper.

**Result:** all forty features have exactly zero saved RNA counts in all
99 development participants, and all forty have exactly zero effective
global-RRR coefficients. Thirty-nine also have exactly zero local and
combined coefficients. `ENSG00000237235` has three local slots, one nonzero
effective region coefficient of `1.4471263907714508e-21` per log2 CPM, and
coefficient RMS `4.6594304646040886e-24` across 96,460 regions. This is a
numerically negligible dependence, not evidence of meaningful feature
influence. It is not rounded to an exact zero.

The defined developmental removal perturbation changes neither the RNA
denominator nor predictions: all direct, denominator and total donor
prediction differences are zero. The synthetic +1-count probes agree with
the actual frozen scorer to maximum absolute differences
`6.6644717461393065e-15` for RRR and `3.3903536729048796e-14` for RRR-plus-local.
The direct-plus-denominator split difference is zero; the Ocis state/global
coefficient check passes with zero nonzero entries. These checks establish
the frozen affine dependency calculation, not biological or external
predictive accuracy.

Execution evidence:

- [summary.json](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/short_rna_feature_dependency_21998516/summary.json)
- [short_feature_dependence.tsv](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/short_rna_feature_dependency_21998516/short_feature_dependence.tsv)
- [developmental_perturbation.tsv](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/short_rna_feature_dependency_21998516/developmental_perturbation.tsv)

Executed helper SHA256 is
`6a875fff4dadc36b27988cd2897c5339bc761e2431330f2c4c884ebc0c74c3e0`;
launcher SHA256 is
`a76ead703edaf2797b672cc07a9097820e35e3b08d0a461629a7739b887c57b2`.

## Data and calculation boundaries

Only the allowed GSE267145 developmental **99 biological participants**,
their 42,163-gene raw fractional-count matrix, exact forty-gene metadata,
and released model factors are read. Every data/source input is checked
against a fixed SHA256 before arrays are opened. Full native modeled-gene
order must agree with the release; no arbitrary representative is chosen
for sequence-identical genes. No H3 values, receiving RNA, receiving H3,
restricted genotypes, protected cohorts, alignments or raw reads are opened.

For raw transport the global effective coefficient on log2 modeled-axis CPM
is `(VtDv @ WVr) / rrr_feature_scale`, multiplied by the region loading
matrix. Local coefficients are the sum of all matching cis slots divided
by their frozen local feature scale. The script **requires every Ocis
state/global coefficient to be exactly zero** before applying this affine
derivation. Otherwise it stops, rather than silently omitting a dependency.
Forty-by-region effective coefficients are summarized in 2,048-region
blocks, avoiding the roughly 32.5 GB full gene-by-region float64 matrix.

Per gene, the output reports developmental count/log-CPM ranges, nonzero
participant count, global latent-factor support, local slot count and
global/local/combined effective coefficient count, maximum and RMS. Exact
zeros refer to the saved floating-point calculation, not an assertion about
biological absence. The audit does not use those results to select features.

The developmental sensitivity sets those forty counts to zero **in copies**,
preserving every gene-axis row. It separates the direct feature contribution
with the original denominator held fixed from the modeled-axis CPM
denominator change affecting all features. Output donor RMS and maximum
prediction differences are changes in predictions, **not H3 prediction error
or receiving accuracy**. This artificial perturbation does not admit zero
filling, missing measurements, or a restricted-axis scorer.

An independent check applies the unmodified frozen public scorer to two
developmental copies, adding one synthetic count to every short feature.
It compares actual profile changes with the derived effective coefficients
for both raw arms across all 96,460 regions. This remains nonvacuous when
all forty original developmental features happen to be zero. There is no
random operation, model training or biological test in this check.

## Conditional bounds and unresolved interpretation

For justified direct log-expression bounds `abs(delta_x_g) <= b_g`, the
RMS triangle bound is `sum_g b_g * RMS(region_coefficient_g)`. The script
reports that bound for the **defined developmental direct perturbation**.
No receiving `b_g` has been established. Developmental minima/maxima cannot
be silently treated as universal receiving bounds, and a count-space bound
must include modeled-axis CPM denominator changes.

Even a finding that all forty developmental counts and direct coefficients
are zero would not prove that a receiving donor has zero biological RNA or
that unobserved counts can be replaced with zero. Conversely, nonzero
coefficients establish a frozen-model dependency without demonstrating
receiving inaccuracy. Complete 42,163 observed-input admission remains
unchanged, and exact private PISCES reconstruction remains unverified.

The observed all-zero developmental counts and negligible direct dependency
do not establish zero receiving biological expression. Unknown receiving
counts could change the modeled-axis CPM denominator and thus other gene
inputs, even with zero direct coefficients. No receiving bound, complete
measurement admission, zero-fill route or private-quantifier equivalence is
established. The complete 42,163 observed-input requirement remains in force.

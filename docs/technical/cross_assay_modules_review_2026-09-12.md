# Cross-assay modules: second review of the v2 candidate

Reviewed 2026-09-12. This is a review, not a change to the analysis or its release status.

**Assessment:** v2 corrects the previous direction error and uses the recorded
multiple-testing families consistently. The results support cross-sectional
associations and paired RNA–chromatin covariation. They do not establish distinct,
reproducible biological modules. Several figure statements still exceed the
tests, and two illustrated modules fail the analysis's own score-stability rule.

Reviewed candidate: `figures/candidates/cross-assay-modules-2026-09-12-v2/`.
Reviewed results: `RNA-seq/results/cross_assay_modules/cam-v2-20260912T130258Z/`.
The current artwork is the run's `figure_v2/`, not its superseded `figure/`.

## Findings requiring correction or explicit qualification

### 1. Figures 5G and S5 turn spatial non-significance into absence of organization

The rendered legend says **“not spatially coherent.”** This labels 142/160
modules in the complete Guilliams column and 73 rows in main 5G. The calculation
only establishes that BH q is at least 0.05 against expression/detection-matched
gene sets. It supplies no equivalence test or upper bound on spatial excess.
Indeed, 90 of those 142 modules have an observed residual Moran statistic above
their matched-null mean. They are inconclusive about excess organization, not
demonstrably unorganized.

The unit handling is improved: five sections are analyzed separately and
collapsed to four donors. That does not make non-significance a negative result.
Use “spatial organization unresolved” or “no significant excess over matched
genes.” Retain the 18 significant matched-null results with their stated scope.
This distinction is also supported by [Pawel et al., *Replication of null
results*](https://arxiv.org/abs/2305.04587).

Evidence: `figS5_cross_assay_module_states_all.pdf`,
`fig5g_cross_assay_module_states.pdf`, the candidate's
`figS5_source.tsv` (example row 962, `cam_m01`), and the run's
`assays/spatial_results.tsv`.

### 2. Figure 5K illustrates two score-unstable modules without marking them

Two of the six illustrated examples fail the stated minimum held-out Spearman
threshold of 0.50:

| Example | Minimum held-out Spearman | Stored classification |
|---|---:|---|
| cam_m147 | 0.124003 | score_unstable |
| cam_m06 | 0.122206 | score_unstable |

Both are excluded from main 5G, but appear in 5K as ordinary examples without
any stability annotation. This obscures a qualification central to the proposed
use of these modules as reproducible scores. Retain the prespecified selection
and explicitly mark these two as unstable diagnostics, or move the examples to
a qualified supplementary display. Replacing them after inspecting results
would require describing the revised selection honestly.

Evidence: rendered `fig5k_prespecified_module_forests.pdf` and the run's
`modules/module_registry.tsv`, rows 7 and 63. The candidate's `fig5k_source.tsv`
already carries the instability fields.

### 3. The snRNA column does not replicate the fibrosis endpoint

The snRNA model uses 70 biological donors and a dataset-adjusted ordinal
diagnosis score: Healthy=0, Steatosis=1, Steatohepatitis=2, Cirrhosis=3. It is
assigned to the fibrosis family, yielding six “replicates” and 16
“tested_negative” cells. A diagnosis contrast involving steatosis and
steatohepatitis does not measure fibrosis stage. Standardization cannot make
these endpoints equivalent.

The source notes acknowledge mixed severity, but the figure and family assignment
still invite a fibrosis-replication interpretation. These are directional
associations with a different disease endpoint. Give the diagnosis result its
own explicit interpretation; its negative bounds also apply only to that
diagnosis score. For example, `cam_m06` has q=0.00201 for diagnosis, not a test
confirming its discovery fibrosis association.

Evidence: the candidate's `module_by_assay_states.tsv`, including row 487;
`CAPTIONS.md`, lines 20–26 and 124–127; and the run's
`assays/snrna_donor_census.tsv`.

### 4. Figure 5I's severity-adjusted interpretation extends beyond its primary ATAC test

The plotted RNA–ATAC result is partial Spearman **given processing well** in
39 donors: 33/150 testable modules pass BH correction over the complete
160-module pairing family. Neither fibrosis nor hepatocyte fraction enters
that primary statistic or its permutation p-value. Nevertheless, the caption
says “Conditioning removes the route through severity.” That statement does
not follow from conditioning on well.

There is reassuring sensitivity evidence: adding the available binary fibrosis
label retains 93.3%–115.5% of the original correlations among the 33 significant
modules, with median 101.5%. That sensitivity uses the 37 donors with observed
fibrosis, whereas the primary analysis uses 39; the ratios therefore also reflect
the sample change. No corresponding adjusted permutation p-values or BH family
are reported. Hepatocyte-fraction adjustment retains a median 96.2% of the
effect, but is also an effect-size sensitivity rather than the plotted test.

Describe the main result as well-adjusted coupling and the other adjustments as
sensitivities. A formal claim beyond recorded histology requires the matched
participant analysis and its own adjusted inference. The separate GSE267145
analysis does adjust for the four recorded histology components and sex in 99
participants; its eight significant couplings remain a different estimand.
Neither arm establishes a shared regulatory mechanism or rules out all
cell-mixture effects.

Evidence: rendered `fig5i_paired_within_participant_transfer.pdf`,
`CAPTIONS.md`, lines 70–88, and `paired_transfer_panel.tsv` (example row 247).

### 5. Figure 5J silently omits coverage for 42 modules in four assays

The coverage source contains missing fractions for all 42 activity-family
modules in GSE268273 bulk RNA, GSE276114 bulk RNA, snRNA, and ATAC. Thus 168
module-by-assay coverage values disappear from the plot because the assay lacks
the selected endpoint, despite the measurements being available.

The underlying assay tables show that all 42 omitted modules satisfy coverage
in each bulk cohort and snRNA; 40/42 satisfy coverage in ATAC. Coverage and
endpoint applicability must be retained separately. Figure 5J should display
all 160 measured fractions in every assay, including measurements that cannot
answer a module's selected histologic question.

Evidence: rendered `fig5j_module_member_coverage.pdf`, candidate
`fig5j_source.tsv` (first missing example at row 6), and the run's
`assays/external_bulk_results.tsv`, `assays/snrna_results.tsv`, and
`assays/atac_results.tsv`.

### 6. Some “measured genes” are ambiguous protein groups

Ten sufficiently detected multi-gene protein groups are assigned to one gene
because only that gene enters the shared universe. Eight assigned genes belong
to modules. The universe itself initially takes the first listed gene from
each group, so this membership rule is not independent evidence resolving the
protein identity. Examples include CALM1;CALM2;CALM3 assigned to CALM1 in
`cam_m147`, and USP12;USP46 assigned to USP12 in `cam_m01`.

Consequently, “all members measured” can include a protein-group measurement
whose contribution cannot be assigned to the named member. Record the ambiguity
and assess coverage without unresolved groups before making a gene-resolved
coverage claim. None of the eight affected modules is currently classified as a
proteome replicate, so this finding does not invalidate the 34 positive
proteome classifications on its own. It does affect coverage and the assay
interpretation of the illustrated `cam_m147` score.

Evidence: `data/PXD051911/liver_protein_quant.txt`, the run's
`modules/module_membership.tsv`, and `assays/proteome_protein_group_rule.tsv`.

## Additional interpretation limits

- Proteome robustness is screened at nominal adjusted p<0.05. Recomputing BH
  over the 160 selected-endpoint complexity-adjusted tests retains 31 of the
  34 displayed proteome replicates. The three others are `cam_m04` and
  `cam_m105` (q=0.0725), and `cam_m49` (q=0.0777). The primary associations
  remain BH-supported; distinguish that from family-corrected robustness to
  the intensity descriptors.
- The Figure 5H caption calls every point to the right of the association
  threshold a replication. The plotted two-sided q also admits the eight
  reversals, and proteome results can be capped as unresolved. The color,
  direction, and stated qualifications are needed to interpret that region.
- The competitive calibration rate 0.056 is a useful check under permuted
  discovery fibrosis labels. It does not measure power to detect specificity,
  or directly calibrate each external assay under its observed disease signal.
  The strongest supported statement is that only one module-by-assay test
  exceeds the connected-set reference after the declared per-assay BH
  correction. “Almost nothing is specific to its own genes” is stronger than
  the test establishes.
- The paired analysis does not compare same-module coupling with coupling for
  matched coexpressed sets. Its 33 and eight discoveries cannot establish
  module-specific regulatory coupling.

## Checks completed and limits

Small-table arithmetic was independently reconstructed with Python's standard
library, without refitting models:

| Check | Result |
|---|---|
| v1 versus v2 module membership bytes | Identical |
| Oriented effect equals raw effect times discovery direction | Exact for every finite applicable row |
| Association BH, six assay families of 160 | Reproduced; maximum rounding difference below 7e-16 |
| Competitive BH, six assay families of 160 | Reproduced; maximum rounding difference below 3e-15 |
| Paired BH, two separate 160-module families | Reproduced; maximum rounding difference below 1.5e-15 |
| Candidate PDF hashes against recorded values | All six match |
| PDF text extraction and visual inspection | All six current PDFs inspected |
| Source-table census | 147 replication cells, eight reversals, 112 bounded negatives; paired calls 33 and eight |

Commands used included `git status --short`, bounded `rg`/`sed`/`cat` reads,
`python3` standard-library table calculations, `pdftotext -layout`, and
`pdftoppm -singlefile -scale-to 1600 -png`. Review-only raster copies were written
under `/scratch/cam-review-20260912-jx82qa9r/`. The existing report of 57 checks
was inspected; that result belongs to the producer's earlier execution.

Full model fitting, graph reconstruction, permutation reruns, the producer's
candidate validation, and the Resource-scope command were **not run** in this
review. No workloads were submitted. Apart from this review document, no
repository file was modified, and no candidate or manuscript status was changed.

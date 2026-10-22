# Frozen source lineage RNA to common-fragment ATAC producer

Prepared for root's source review on 2026-10-01. No fit, prediction, numerical test or job was run for this preparation. This implements the frozen exposed-source development comparison; it does not establish independent validation, receiving admission or a usable model.

## Executable files and authorization

| File | Bytes | SHA256 |
|---|---:|---|
| `scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse296875_lineage_model.py` | 53,329 | `ccb6552d61b31f0d3fcf607f6b677c2b44909eaa6bdd6c1ae9777b2ec8e90430` |
| `scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_gse296875_lineage_model.sbatch` | 1,225 | `3466957e77b51542314c02fc8e95d0a13327cd700c8b4782cf8074d160086818` |

All execution paths are relative to the repository root. Let `REC` mean `Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z`.

The fully read recipe is `REC/20261001_codex_gse296875_common_fragment_rna_atac_prespec.json`, SHA256 `17ab84770eb3bfac46ea50d0a1400e52a70ffea1bbb7683e45b00a40e499403c`. Its original preparation-stage `fitting_authorized_by_this_file=false` remains unchanged. Root's separate `REC/20261001_codex_lineage_source_execution_decision.json`, SHA256 `4fbe34ae9acf5b2c468d5f75fc59965810379f1a738504806bcb20a3c083e57f`, authorizes source fitting after source review and forbids external evaluation or scientific-authority changes. The launcher pins both the helper and this decision. The helper checks exact decision source paths, summary hashes and fixed dimensions before array access; no circular model/decision hash is introduced (`main`, line672).

The independent NumPy kernels do not import or modify the concurrently developed global model. Both frozen branches, including failures, must be reported without selecting a branch by source performance. The 117 Hotspot definitions, other hypotheses, authorities and previous results are unchanged.

## Source and replication

The 15 fixed small source/identity/summary hashes in the helper matched during static preparation. Only metadata/identity files were inspected here; no RNA or ATAC array body was opened.

| Source directory under REC | Summary SHA256 | Artifact-list SHA256 |
|---|---|---|
| `codex_gse296875_aligned_native_pseudobulk_21999585` | `dd1986b643e1e5ae2cec54d255fb8ed7619c6db668541194a4032277d8c00226` | `3c4102ecde36825237b7d31e6ff4c4395570dc97ee222e274429833cc000d8a5` |
| `codex_gse296875_common_fragment_counts_22000534` | `bd64347bc066cbd182695c466b957ddd3c6efc689b53f8db2d557aed34321f8e` | `8735af6f215fc82be928f2349cb8182d3d76e834304cce8eb32e858f5880a960` |
| `codex_gse296875_donor_global_raw_22000941` | `9684d83681b898ac507680e2df39a4781b42daccbd2abfac6342faa88031fb19` | `5629741635252c1869cc68da59c365833ad8912fa7f46d62fef25335beead111` |

Compute execution hashes each of the three lineage RNA, three common-fragment ATAC and three donor-global RNA arrays before opening it. It checks shape, dtype/order, integer range, full=A+B conservation, donor-global raw sum across seven lineages, exact unit identities and denominator conservation (`guarded_arrays`, line303). It does not open the older native ATAC or global ATAC arrays.

Biological n=39 people; eight wells are the held-out and resampling clusters. There are 273 stored units, 272 sampled full units and 270 fixed two-nonempty-half units. Donor733/B cells is unsampled, not a zero observation. Donor382/B cells and751/B cells are excluded from the half diagnostic by pre-existing nucleus membership alone. All39 people remain in each complete donor-equal endpoint.

The source-metadata census gave these training counts in well1 through well8 order:

| Eligible population | B cells | Each other lineage |
|---|---|---|
| Full | 33,33,34,33,33,33,34,33 | 34,34,35,34,34,34,34,34 |
| Both nonempty halves | 31,31,33,31,31,31,33,31 | 34,34,35,34,34,34,34,34 |

The lineage and donor-global PCA arms use these identical training people separately within every lineage/fold. The global RNA row is used once per eligible person, rather than repeated once per lineage within a fit (`train_mode`, line502).

## Measurement and fixed fits

RNA uses the complete original ordered36,601-gene integer UMI axis and `ln(1+1e6*UMI_g/sum_all36601_UMIs)`. Donor-global counts are summed raw across all seven native lineages before normalization; symbols do not collapse IDs. Technical inputs are target-lineage ln1p RNA total, ln1p nucleus number and detected-gene fraction on36,601 genes.

The ATAC target is `ln(1+1e6*C_r/D)`, with C the separately declared unique-fragment-body overlap count and D the matching partition's autosomal unique-fragment count, including fragments outside counted windows. Each of the306,706 exact fixed windows is fitted directly. No summed-window, primary-chromosome or PCR-support denominator is substituted; no extra Tn5 shift is applied. This is the new explicit measurement, not proof of equivalence to original author ChromatinAssay counts.

Train-only transforms use gene centering, full float64 thin SVD and fixed20 RNA slots, structural availability `min(20,ntrain−1,G)`. There is no nearzero singular-value mask, rank selection, HVG selection or randomized SVD. Signs follow the first maximal absolute loading in frozen gene order. Exact constant training score/technical columns are identified by equality of their actual row values and map to zero; this avoids treating summation-roundoff SD as real variation without introducing a tolerance threshold. All singular values, exact adjacent ties, loadings, centers and sample-SD scaling are saved (`pca_fit`, line196).

Each lineage/fold fits mean, technical-only, technical+20 donor-global RNA PCs and technical+20 target-lineage RNA PCs. Ridge is the fixed sum-SSE objective with alpha10, unpenalized intercept and every non-intercept coefficient penalized. Execution checks finite coefficients and centered normal-equation residuals (`fit_ridge`, line225). Both raw and elementwise nonnegative-clipped out-of-fold predictions are saved; every metric uses clipped predictions.

There are168 lineage/fold groups across full and two independently fitted opposite-half directions, containing672 arm fits, then28 final all-sampled-donor arm/lineage fits. Training-only transforms and all axes remain fixed in both half directions. Undefined required full RNA/ATAC denominators prevent complete primary fitting; no value-dependent row removal or fabricated normalization rescues it. Undefined required half measurement makes that complete half diagnostic unavailable, preserving the full roster.

## Complete frozen endpoints

The primary endpoint averages squared error across all windows, available lineages within each person, then equally across39 people. It reports all four absolute losses, target-lineage minus donor-global RNA, both utility guards, and all unit/person/lineage/well results. The seed20261001 bootstrap draws eight whole wells10,000 times with multiplicity and uses each draw's actual resampled donor-count denominator. It does not average unequal-sized well means, refit models or claim a population p-value (`aggregate_losses`, line348).

All unordered held-well/lineage donor pairs are evaluated for all four arms, after clipping. Complete pair rows include signed bias, observed/predicted squared differences and their cross-product. Those summaries use equal pairs within lineage/well, available lineages within well, then equal wells. The mean-arm difference is explicitly checked to equal zero. The centered-error algebra identity is checked on compute; pair count never becomes biological n (`full_diagnostics`, line396; `diagnostic_aggregates`, line449).

For target-lineage RNA only, every recipient is paired with every other eligible same-well/lineage RNA donor. The recipient technical vector remains unchanged; only frozen lineage PC values are substituted. Every wrong prediction is clipped before loss, all wrong losses are averaged per recipient, and all recipient, person, lineage and well summaries are retained. There is no shuffle, favorable donor selection or use of this diagnostic for advancement.

Both half directions preserve all four arms and all270 metadata-eligible units. Paired half-target agreement is descriptive; no reliability-selected units/windows, correlation ceiling or added independent donors are claimed.

Final `LineageScorer` (`line604`) accepts complete ordered target-lineage integer RNA, positive nucleus number and one fixed lineage label. The global comparator additionally requires complete ordered raw donor-global RNA. It applies only saved training transforms, parameters and clipping, returning all306,706 outputs in original order. No ATAC, H3, clinical, histology, genotype, allelic or receiving-centering input is accepted. Compatible external lineage labels and measurement remain separate unresolved requirements.

On 2026-10-02, the constructor was corrected to require a successful `summary.json` before parameter loading: completed development, a complete primary endpoint, all28 final scorers and identical before/after source inventories. A present `failure.json` prevents construction. `artifacts.json` must bind the final scoring contract by both byte count and SHA256; summary, contract and the executing producer must agree on recipe and producer identity. This prevents an early contract from a subsequently failed run from being accepted as a completed source model. No fitting kernel, parameter, measurement, endpoint or frozen recipe changed. No numerical callable check or fit was run for this correction; Python AST parsing, launcher syntax and source-hash checks passed on login.

## Resources, saved outputs and checks

Prepared individual header:1595, cpu, account/qos nslab,1CPU,8G,90h; distinct actual-job `%j` output/error logs. Existing rnaseq Python is used with all BLAS/OpenMP thread settings1, no dependencies installed. Execution creates only fresh `REC/codex_gse296875_lineage_model_${SLURM_JOB_ID}`; existing outputs cannot be overwritten and partial failures remain.

[Estimate] OOF arrays occupy16.08GB, three-mode fold parameters approximately23.91GB and final parameters approximately1.00GB, with small tables/source copies in addition: approximately41GB total, below the64GiB inclusive cap. The largest working matrices are direct training targets and ridge right-hand sides/coefficients, roughly0.5GB together, plus mmap pages and small RNA/SVD workspaces;8G allocation provides margin. Arrays are held per fold, prediction mappings are closed after writes, and wrong-donor prediction arrays are not accumulated. These are allocation estimates, not measured MaxRSS or runtime. The internal88h check leaves a scheduler buffer; compute checks inclusive saved bytes throughout and after final summary. No cleanup deletes partials.

Actual preparation checks: Python AST parse passed; `bash -n` passed; helper/launcher/decision hashes matched; all15 small source guards matched; metadata-only full/half membership and per-fold counts matched. No helper import, numerical kernel, callable execution, array conservation test, bootstrap or model fit was run on login. Runtime will save Python/NumPy/config/thread inventory, executed source copies, before/after source guards, every parameter and prediction file hash and categorical failure receipts. Root must review the frozen executable source before any submission.

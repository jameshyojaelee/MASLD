# GSE223843 ordered barcode namespaces and prospective donor-global RNA inputs

Prepared source-only/offline qualification on 2026-10-01. No new submission, acquisition, matrix inspection, prediction, or evaluation was performed. This is a distinct preparation for the frozen donor-global function, not an amendment to prior geometry results or the seven-lineage task.

## Count unit and complete RNA denominator

[Official ARC 2.0 matrix documentation](https://www.10xgenomics.com/support/software/cell-ranger-arc/2.0/analysis/outputs/feature-barcode-matrices), independently retrieved today, distinguishes Gene Expression UMI entries from Peaks cut-site entries. Matrix column IDs use the Gene Expression barcode namespace. The feature and barcode files give matrix row and column order; raw columns include background, while filtered columns are cell-associated. These documented units must remain separate in a mixed matrix. The installed ARC 2.0.2 joint-matrix implementation stacks the modalities and copies their original data without changing the RNA unit: `lib/python/atac_rna/matrix.py:20–36,49–109,112–120` (SHA256 `3795c9c494845e208c4f10ae9f53bb7e51199badec394221e17e3aa06caf47f8`).

All six annotation receipts from completed job **22000649** establish an exact ordered **36,601 Gene Expression ID/name axis**, digest `afbcc02143e68d1fbcda9631917f1344adfd3ac8efb752e2ebb7f044df934e60`, plus these additional Peaks rows:

| Author person | Gene Expression | Peaks |
|---|---:|---:|
| Control1 | 36,601 | 16,775 |
| Control3 | 36,601 | 82,811 |
| Fontan1 | 36,601 | 38,143 |
| Fontan2 | 36,601 | 43,878 |
| Fontan3 | 36,601 | 33,121 |
| Fontan4 | 36,601 | 19,973 |

The complete annotation summary is `Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/codex_gse223843_rna_annotation_identity_22000649/summary.json`, 33,296 B, SHA256 `da3f1310cfba1bc145ded47979a56e4134e2870baf09a25d553fc598024a7abe`. Its six individually hashed annotation receipts are also guarded by the new helper. Full annotation tables are preserved there; they are not re-decoded by this barcode task. The prior six unsupported narrow-parser results remain historical parser failures, not evidence of invalid RNA.

The frozen donor-global recipe is `Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/20261001_codex_gse296875_donor_global_rna_atac_prespec.json`, SHA256 `88b9e4239b08cc1ccbe78b092e8eacd395778453a496fd8915b4acc9e2729ae1`. It pools raw RNA UMIs over the fixed donor nucleus roster **before** applying `ln(1 + 1,000,000 * U[g] / sum_h U[h])`. The denominator is the complete 36,601 Gene Expression genes. Peaks cut sites must never enter it. The technical totals and detected-gene fraction also use that whole gene axis, without symbol collapse, lineage weighting, absent-feature filling, or receiving-data feature selection. A zero total or missing measured axis is unmeasurable under the recipe.

[Unverified] Neither uploaded MTX bodies nor the author's exact export command have been inspected. ARC methods and compatible annotations support a candidate native-UMI route; they do not establish that the uploaded MTX entries are unnormalized integers, or that its barcode filtering equals development filtering. The original installed filtered-joint producer requires common modality barcode sets and can restrict RNA to the ATAC set (`mro/atac_rna/stages/processing/merge_atac_rna_matrices/__init__.py:100–113`, SHA256 `eb1cbb59f50afb811110d790f558270f7e225aef30a07a2babace09dddc30ac3`). The observed uploaded barcode sets are unequal. Their exact intersection is therefore a separately declared source-metadata roster, not proof of an unchanged native filtered-joint matrix or arbitrary RNA-only sampling compatibility.

## Outcome-independent technical roster

The acquired identity metadata fixes the six author-person/GEX/ATAC accession pairs. Completed job **21999912** established these exact sample-namespaced joins without counts:

| Author person | GEX / ATAC accession | RNA barcodes | ATAC barcodes | Exact shared |
|---|---|---:|---:|---:|
| Control1 | GSM6997741 / GSM6997754 | 10,686 | 1,346 | 1,319 |
| Control3 | GSM6997744 / GSM6997756 | 19,506 | 15,275 | 15,186 |
| Fontan1 | GSM6997746 / GSM6997758 | 6,977 | 2,482 | 2,462 |
| Fontan2 | GSM6997748 / GSM6997760 | 6,844 | 4,674 | 4,605 |
| Fontan3 | GSM6997750 / GSM6997762 | 12,360 | 4,569 | 4,538 |
| Fontan4 | GSM6997752 / GSM6997764 | 6,359 | 1,218 | 1,217 |

The new helper reproduces these exact joins and records **all** intersection rows, in original RNA barcode-file order, with both source row indices, exact source accessions, and untouched GEM suffixes. It does not use RNA positivity, ATAC events, lineage, clinical state, or prediction accuracy to choose rows. It does not deduplicate persons using barcode strings. This specifies 29,327 technical intersection rows across six author persons; rows are not biological replicates.

The saved identity metadata points to [source Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC11103255/) for within-sample shared cell barcodes. That source-method pointer is retained as existing evidence; the independent browser fetch today returned a recaptcha rather than Methods. Barcode coincidence alone does not verify actual well identity, paired-nucleus identity, or person disjointness from the development donors. Those limitations remain explicit in every new roster receipt. The existing pediatric population qualification remains unchanged.

## Canonical ARC codebook qualification

Installed root: `/nfs/sw/easybuild/software/cellranger-arc-2.0.2`.

| Reference | Compressed bytes | SHA256 |
|---|---:|---|
| RNA `lib/python/cellranger/barcodes/737K-arc-v1.txt.gz` | 2,142,553 | `1279bb7ba28ddcd805e2f16ebd5d6933b7b1e7f128e18aabfe6ae8cd7b7bf710` |
| ATAC `lib/python/atac/barcodes/737K-arc-v1.txt.gz` | 2,465,078 | `9a4c8f29a05f59b9a4d7525a2f44d1778ccfa540a8832ad9f5f84047f374cc26` |

Installed `lib/python/atac_rna/barcodes.py:21–27` loads RNA with `as_set=False`, ATAC with `ordered=True`, then constructs both dictionaries from **the original ordered zip**. SHA256 `9f4e0c90eec48ac8bbadd0bd69377c8fca9b1a17cac59b216283dd50fc82f185`. Its lines 29–65 split and rejoin the GEM suffix. RNA loader `lib/python/cellranger/barcodes/utils.py:46–52` explicitly checks duplicates and preserves list order; SHA256 `e9c5c807d91af0b5ea9fc1573ab590fd48e0d0ade711057ad8aea36d78d070f9`. The ATAC loader likewise returns a list in ordered mode, with its original comment rule; SHA256 `c896d07f2bde284d6dd36a618636e2c83d9b131a2195af627dc9eebd8c9069ba`.

The new offline helper reads both full references with CRC checks, requires unique 16-base canonical identities and equal nonempty lengths, and saves both lists plus their original-order pair index. It verifies inverse correspondence, rather than pairing sorted barcode strings. Both references, both loaders, translation code, chemistry definition, MRO choice, and package path are pinned before and after use. No orientation selection, barcode error correction, or processed-ID translation is performed.

Each of the 12 saved processed barcode files is guarded by the frozen proposal, exact identity/accession mapping, original TAR offset and HTTP206 receipt, actual compressed size/SHA, and prior full decoded identity digest. It is independently decoded to EOF with CRC, narrow one-field barcode validation, full-ID uniqueness, and native order preservation. Base membership is reported separately as RNA-only, ATAC-only, both, or neither; source suffixes are retained. Literal namespace overlap is diagnostic. It cannot choose a translation, override source namespace declarations, discard rows, or prove cell identity. The full codebook count and actual processed namespace census are **not run** in this preparation.

## Prepared files and resources

- `scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse223843_barcode_namespace.py`, 24,706 B, SHA256 `cb269cd1f6e833bec756603e6d3bf60655cd5b140f8add4b0b3c80b5d252745c`.
- `scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_gse223843_barcode_namespace.sbatch`, 944 B, SHA256 `2aac51476b856a6bc3c3ec0a289e50ffeb30a29145ba1dec1d696034608a0aed`.
- Fresh compute output only: `.../codex-rec-20260929T142434Z/codex_gse223843_barcode_namespace_${SLURM_JOB_ID}`. Existing outputs are refused; partials are retained. Distinct `1592_barcode_namespace_%j.out/.err` logs preserve the actual scheduler ID.

Proposed single job: name1592, cpu/nslab account and QoS, 1 CPU, 2 GiB, 1 hour. Runtime enforces 55 minutes, no array, 32 MiB decoded per file, 128 MiB decoded aggregate, and 128 MiB saved including scheduler logs. Every output write reserves space, including buffered bytes; 8 KiB remains for a categorical failure receipt. [Estimate] peak memory below 0.6 GiB and saved output around 65 MiB; full whitelist length remains to be measured. There are no network requests, matrix bodies, feature-body redecodes, sequences, genotypes, allelic values, H3 values, clinical values, fits, or outcome inspection.

Executed on login: Python AST parse; guarded source-only import without `main`; all 17 fixed source SHA checks; six small annotation receipt size/SHA checks; 12 compressed barcode size/SHA checks; `bash -n` launcher. All passed. No gzip decoding, namespace computation, biological numerical work, or submission was performed. Root reviews the frozen files before any execution.

## Separate future mixed-MTX extraction rule, not implemented or authorized here

After separate review, a future RNA producer can use the qualified **complete** annotation row map and unchanged barcode-column order. It must first establish the uploaded MTX header, dimensions, count provenance, and native sparse-zero contract. It must reject missing axes, unknown types, duplicate entries, integer overflow, and fractional/negative/nonfinite RNA counts. The frozen source-metadata intersection roster must be fixed before values or prediction accuracy are inspected.

For each mixed-MTX record, parse only the feature row index first. If its declared type is Peaks, skip the remaining record bytes without interpreting or saving its column, numeric value, or outcome-derived support. Only Gene Expression records may have their columns and integer UMI values interpreted. Pool over every fixed paired-roster nucleus, preserve all 36,601 native genes, then apply the frozen transform. Native sparse zeros on the complete declared axis are genuine matrix semantics; absent genes are not invented zero rows. No peak values may choose the roster or enter the RNA denominator.

A mixed compressed stream necessarily transports/decompresses bytes that also encode Peaks entries. A later producer must disclose that physical fact; it cannot claim those bytes were never acquired. The enforceable scientific boundary is no ATAC-value interpretation, retention, outcome-based filtering, or scoring before the frozen RNA predictions. No such stream or count producer is executed in this task. Whole306706 fragment measurement, source-filter compatibility, source-person disjointness, and receiving measurement/evaluation admission remain unresolved.

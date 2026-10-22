# GSE223843 six-field RNA annotation qualification preparation

Prepared 2026-10-01 for root review. **Not submitted or run.** This corrects the annotation schema used for qualification while preserving prior helpers, jobs 21999459/21999912, original compressed files, reports and the read-layout pilot. Earlier “unsupported” RNA feature results describe that parser's limits; they do not establish invalid RNA annotation, missing genes or biological absence.

## Source evidence

[Official ARC 2.0 feature-barcode matrix documentation](https://www.10xgenomics.com/support/software/cell-ranger-arc/2.0/analysis/outputs/feature-barcode-matrices) describes six columns: ID, name, type, chromosome and ARC TSS start/end in BED coordinates. `program_ratio_hypothesis` independently retrieved this page, including its Gene Expression/Peaks categories. This agent's two browser requests timed out; no HTML receipt was saved. The installed primary implementation was independently inspected and pinned:

`/nfs/sw/easybuild/software/cellranger-arc-2.0.2/lib/python/atac_rna/feature_ref.py`, 7,056 bytes, SHA256 `1b80c8ce8675e1fe0c079f0167fc4209a84712f0bcf28162465251ac33c2ac51`.

Lines 177–196 write all six fields for every feature, including both declared types. Gene Expression coordinates describe TSS annotation, not gene-body bounds. Lines 22–24 and 96–104 also define a missing interval as empty chromosome and `-1/-1`; the helper retains this exact sentinel as unavailable TSS annotation, not an invented coordinate. Peaks use their own interval-ID/BED semantics. The source paper used ARC 2.0.0; installed 2.0.2 schema evidence does not prove reference/count equivalence.

## Frozen preparation

| New artifact | Bytes | SHA256 |
|---|---:|---|
| `scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse223843_rna_annotation_identity.py` | 25,183 | `8620b8365044673de0c39ed2e1fe9875b88057c93860c57b4fb04ec484a8a804` |
| `scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_gse223843_rna_annotation_identity.sbatch` | 944 | `09055cf3c0462f0e3b6a28c0828d2da34970811bddef6132a9236de7a4692498` |

One proposed individual 1592 job: `cpu` partition, `nslab` account/QOS, **1 CPU, 2 GB, 1 hour**. Fresh output is `REC/codex_gse223843_rna_annotation_identity_<actual_jobid>`; logs are `1592_rna_annotation_%j.out/.err`. `REC` is `Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z`. Decode and inclusive saved caps are **256 MiB**, per-member decoded cap **64 MiB**, line cap 4,096 bytes and two million rows/member. Tracked writes include buffered tables; 1 MiB is reserved for small logs and filesystem totals are periodically checked. [Estimate] Six small annotation files should take seconds to a few minutes; actual runtime/memory is unmeasured. No network, external tool process, new dependency or install is used.

Completed login checks: source/hash/header inspection, Python AST parse, `bash -n`, executable permissions and exact helper hash in launcher. No feature decompression, complete gene-axis comparison, CRC computation, numerical model work or submission ran on login.

## Exact admitted annotation inputs

Only the following six existing files from `REC/codex_gse223843_feature_barcode_ranges_21999459` are opened for feature decoding. Total compressed bytes: **8,822,584**. Each file's exact saved 206 range start/end/size/hash, original proposal and GEO person/assay identity must agree before decoding; hashes are rechecked at completion.

| Person / RNA GSM | Member | Bytes | SHA256 |
|---|---|---:|---|
| Control1 / GSM6997741 | `GSM6997741_GEX_Ctrl_features.tsv.gz` | 1,034,803 | `96c0f53ba3d3bb979ca9580f269503128c27163387b2e24d94b3f0295dcb0c38` |
| Control3 / GSM6997744 | `GSM6997744_GEX_Ctrl3_features.tsv.gz` | 2,317,012 | `1eecc3696db023bc8a014ef650ac73b25d5fb7ff9ecc85dff4ccdd49870cf430` |
| Fontan1 / GSM6997746 | `GSM6997746_GEX_Fontan1_features.tsv.gz` | 1,453,275 | `cf6e3ec34def8b956d103b163d6ba508c773e11dd826928614a7f0410e973a36` |
| Fontan2 / GSM6997748 | `GSM6997748_GEX_Fontan2_features.tsv.gz` | 1,565,546 | `394a841689d1f36f160e8ed68db2ad8ced93fcdcbf4c5857c0bf9d2209ef3220` |
| Fontan3 / GSM6997750 | `GSM6997750_GEX_Fontan3_features.tsv.gz` | 1,354,139 | `0d977c633b3ff1f846395dbf36db25c0f6ff59d86207b1fea4388f7414825f0c` |
| Fontan4 / GSM6997752 | `GSM6997752_GEX_Fontan4_features.tsv.gz` | 1,097,809 | `400ac1358bf181a7701244ef624a4b317f13507aa1ffba03ed7d87b5304ed40d` |

The frozen comparison is **only gene identity metadata**: `REC/codex_gse296875_aligned_native_pseudobulk_21999585/genes.tsv`, 1,102,092 bytes, file SHA256 `2a33d4d339e2c6db695525ae2655af5156ccda008300bd992e306d48f01d38f0`. Its header is `gene_index`, `ensembl_id`, `native_gene_symbol`; all 36,601 indexed rows must reproduce ordered ID/name digest `afbcc02143e68d1fbcda9631917f1344adfd3ac8efb752e2ebb7f044df934e60` using UTF8 ID, NUL, UTF8 name, LF. This is checked against `Analysis/SingleCell/continual_integration/reference/gse296875_source_lock_v44.json`, SHA256 `6be686b4187309f038097cf73d3f2a8dd34e72d666bc92733a53e99a1ccbb620`. No count array, HDF5 payload, barcode file or matrix is opened.

## Qualification rules and intended outputs

The helper streams each source feature file in native order and validates six typed fields. Gene Expression requires a human Ensembl gene ID, bounded printable UTF8 name and either a valid source BED TSS interval or the exact missing-TSS sentinel. Names are source annotations; punctuation does not redefine a gene symbol. Peaks require their declared coordinate identity to agree with their separate BED interval. Names/IDs are preserved; no symbols are collapsed, genes aliased or rows selected by values.

Known categories are counted separately, and every validated source annotation row retains its original feature index/type/fields in an annotation TSV. Unknown type, malformed width/identity/coordinate or invalid text stops that member without scanning its remainder; only row index, field count, fixed reason code and unknown-type digest are recorded. CRC and complete-axis qualification then remain false. Other members continue. No unknown row is silently filtered while declaring a complete file. Valid complete iteration checks gzip CRC. Duplicate feature/gene/versionless IDs and repeated names are reported without dropping rows.

For each person a **metadata join** retains all 36,601 frozen reference rows and original source positions. It classifies exact IDs, uniquely version-equivalent IDs, absent IDs, ambiguous exact/versionless duplicates and extra source gene rows. “Absent” means missing from the validated annotation axis, not an observed zero; comparisons are explicitly partial when schema/CRC is unqualified. Ordered exact ID equality, ordered ID/name equality and ordered versionless-ID equivalence are reported separately. Version stripping is used only for descriptive classification and never renames the saved source axis or normalizes counts.

Each person receives an annotation receipt, native annotation table, unsupported-row categorical table, complete fixed-reference identity join and extra-source identity table. `summary.json` gives source/CRC/category/duplicate/match counts, digests, output inventory and caps. Partial failures preserve files and categorical phase without printing rejected values. It distinguishes mixed Gene Expression/Peaks files from RNA-only feature files and leaves actual RNA library-denominator/count/reference compatibility explicitly unverified.

Even a complete ordered native RNA identity match would not admit a model evaluation: counts, paired lineage identities, RNA processing/denominator, exact ATAC recounting and source independence remain separate requirements. This preparation changes no model, fixed feature axis, normalization, original artifact or previous failure record.

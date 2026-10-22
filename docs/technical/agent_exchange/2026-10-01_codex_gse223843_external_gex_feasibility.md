# GSE223843 external GEX and annotation feasibility

2026-10-01. Bounded source-only review. This report is the only file created. Public run/file metadata, primary methods, installed ARC source, and existing identity/schema receipts were used. No sequencing body, matrix body, clinical table, genotype, prediction or accuracy output was inspected; no installation, fitting, acquisition job or submission was performed. Primary-page rendering included surrounding material; only assay methods and annotation definitions were used as evidence. This review does not authorize full-run acquisition or receiving evaluation.

**Finding:** all six source people have two identifiable raw GEX sequencing lanes. The previously rejected six-column RNA features match an explicitly documented ARC format, so that parser failure does not establish unusable RNA. Canonical ATAC-to-GEX translation is available in installed source. **[Unverified]** Complete feature identities, compatible raw RNA quantification/reference, a seven-lineage RNA-only receiving annotation and biological-person disjointness remain prerequisites; none is resolved by listing runs or identifying a translation function.

The separate frozen development recipe is `Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/20261001_codex_gse296875_common_fragment_rna_atac_prespec.json`, 26,762 bytes, SHA256 `17ab84770eb3bfac46ea50d0a1400e52a70ffea1bbb7683e45b00a40e499403c`. Its full 306,706 counted windows and 36,601 RNA genes are unchanged. This review adds no external admission or model claim.

## Exact source-derived GEX inventory

The existing GEO identity receipt supplies the exact author labels, RNA accessions and SRA experiment/BioSample relations. A new bounded [official ENA run metadata query](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA928662&result=read_run&fields=run_accession%2Cexperiment_accession%2Csample_accession%2Csecondary_sample_accession%2Clibrary_strategy%2Clibrary_layout%2Clibrary_source%2Clibrary_selection%2Cinstrument_model%2Cfastq_ftp%2Cfastq_md5%2Cfastq_bytes%2Csra_ftp%2Csra_md5%2Csra_bytes%2Cread_count%2Cbase_count&format=tsv) was restricted to those six RNA experiments. ENA's `sample_accession` contains the matching SAMN identifiers; `secondary_sample_accession` contains SRS identifiers. An initial guard expected SAMN in the latter column and stopped before reporting rows; the corrected guard explicitly matched SAMN in either declared identity field. No sample or assay label was changed.

| Author person | RNA GEO | BioSample | RNA experiment | Exact two GEX runs |
|---|---|---|---|---|
| Control1 | GSM6997741 | SAMN32937892 | [SRX19197932](https://www.ncbi.nlm.nih.gov/sra/SRX19197932) | SRR23252466; SRR23252467 |
| Control3 | GSM6997744 | SAMN32937891 | [SRX19197931](https://www.ncbi.nlm.nih.gov/sra/SRX19197931) | SRR23252468; SRR23252469 |
| Fontan1 | GSM6997746 | SAMN32937890 | [SRX19197933](https://www.ncbi.nlm.nih.gov/sra/SRX19197933) | SRR23252465; SRR23252470 |
| Fontan2 | GSM6997748 | SAMN32937889 | [SRX19197936](https://www.ncbi.nlm.nih.gov/sra/SRX19197936) | SRR23252457; SRR23252458 |
| Fontan3 | GSM6997750 | SAMN32937888 | [SRX19197934](https://www.ncbi.nlm.nih.gov/sra/SRX19197934) | SRR23252462; SRR23252463 |
| Fontan4 | GSM6997752 | SAMN32937887 | [SRX19197935](https://www.ncbi.nlm.nih.gov/sra/SRX19197935) | SRR23252459; SRR23252464 |

These are six author person labels and twelve technical runs, not twelve biological replicates. The original file names explicitly distinguish L001/L002. RNA/ATAC linkage remains the existing exact author-person records plus the source's same-nucleus methods; no new pairing is inferred from suffixes, run adjacency or accession numbers, and no molecular identity check was run.

The [official NCBI SRA package response](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=26398130%2C26398129%2C26398128%2C26398127%2C26398126%2C26398125&retmode=xml) reports four stored classes of 10/10/28/90 bases for every run, with zero reported length SD and equal class spot counts. The primary library method specifies 28:10:10:90 GEX cycles. The normalized-file identities are:

| Run | Normalized SRA bytes | Source MD5 |
|---|---:|---|
| SRR23252457 | 7,438,649,348 | `0f6deaefae7ccd38468f47168e405c5e` |
| SRR23252458 | 7,445,575,936 | `71397c0a8f58761e826f1c51c4217e9b` |
| SRR23252459 | 3,572,014,629 | `eb4060a5da6848064c4610f817180a94` |
| SRR23252462 | 5,312,257,797 | `0ce84e8ca76372446563b5ef39e6c048` |
| SRR23252463 | 5,305,151,621 | `59c5b43c2e7919dac29c65a7f97a670a` |
| SRR23252464 | 3,569,887,921 | `11b5f188ade0df535c11a72d7d6a7c16` |
| SRR23252465 | 3,642,367,912 | `24c995db93e94617ed2064fb647e1b87` |
| SRR23252466 | 11,365,114,817 | `617e0cad92bb94f8b1fedcf05c233ce7` |
| SRR23252467 | 11,365,949,645 | `d985432126ed025ffc6cb3e140bad0f3` |
| SRR23252468 | 4,838,945,053 | `de241cfb176d3fdce060953559e63b0a` |
| SRR23252469 | 4,835,974,949 | `79f9f63562cf9c3410049b20c35ae5b4` |
| SRR23252470 | 3,646,197,479 | `314fac9bc0fba98f1e2336b4825e1539` |

The metadata explicitly lists normalized HTTPS locations of the form `https://sra-pub-run-odp.s3.amazonaws.com/sra/<run>/<run>`. They were not downloaded or HEAD-tested here. This is listed file metadata, not independently verified availability or content integrity. Normalized files total **72,338,087,107 bytes**. Original four-class FASTQ metadata totals **127,709,702,431 compressed bytes**. ENA lists two FASTQs per run totaling **80,743,841,718 bytes**; its base totals equal 118 times its reported spot counts, whereas NCBI's four-class totals equal 138 times the same counts. **[Inference]** This is consistent with ENA preserving the 28+90 GEX pair and omitting the two sample-index classes. Actual read-role export, synchronization, barcode/UMI retention and original-file accessibility remain untested. These sizes are not an acquisition proposal.

## The six-column RNA annotation has a documented meaning

The [official ARC 2.0 feature-matrix documentation](https://www.10xgenomics.com/support/software/cell-ranger-arc/2.0/analysis/outputs/feature-barcode-matrices), rendered lines 103–105, defines columns 1–6 as gene ID, gene name, feature type, chromosome, start and end. For a Gene Expression row, the interval represents the ARC-determined transcription start region in **0-based BED format**, not the entire gene body. The gene identity comes from the reference GTF. A matrix may also contain Peaks rows; its RNA denominator must not include those rows.

Installed `/nfs/sw/easybuild/software/cellranger-arc-2.0.2/lib/python/atac_rna/feature_ref.py:177–196` writes exactly those six fields. Lines 22–25 and 92–110 additionally support missing interval annotation (`NA` internally, serialized empty chromosome and `-1/-1`). Such missing TSS annotation must be distinguished from an invalid or absent gene ID. This source is ARC **2.0.2**, while the paper specifies **2.0.0**; it establishes an installed schema definition, not historical binary/reference identity.

The existing first-line receipt classifies the acquired Control1 feature line as these six lexical fields without saving values. **[Inference]** It is consistent with standard ARC output. It does not establish all six files' complete schema, CRC, duplicate-ID policy, feature-type census, ordered gene identities or the complete frozen 36,601-gene denominator. Those checks remain unrun by this review. The former narrow parser's rejection of six fields therefore does not establish RNA unavailability.

The fixed native identity metadata is `REC/codex_gse296875_aligned_native_pseudobulk_21999585/genes.tsv`; the ordered gene-ID/symbol content hash in the source lock is `afbcc02143e68d1fbcda9631917f1344adfd3ac8efb752e2ebb7f044df934e60`. This review did not decode that axis or the acquired external feature members. An annotation-only compute qualification of already acquired files can address identity/schema/CRC and fixed-axis coverage without opening matrix values. It cannot itself validate UMI count semantics, infer zero expression for an absent reference gene, choose an RNA denominator or admit a receiving scorer. The parallel `native_rna_index` agent owns that preparation; no duplicate helper was created.

## Canonical barcode translation, source evidence only

The [official ARC barcoded-BAM definition](https://www.10xgenomics.com/support/software/cell-ranger-arc/2.0/tutorials/outputs/cr-arc-atac-barcoded-bam) states that corrected ATAC `CB` values are translated into their paired GEX barcode sequences. ATAC and GEX primer barcode sequences on one gel bead differ. Raw `CR` and its quality remain separate. Shared processed RNA/ATAC barcodes do not justify using a raw ATAC index substring as a GEX barcode.

Installed `lib/python/atac_rna/barcodes.py:22–27` loads the RNA whitelist with `as_set=False` and ATAC whitelist with `ordered=True`, then pairs them by `zip_longest` in both directions. Lines 29–62 strip and restore GEM group suffixes, preserving the group during translation. ATAC loader `lib/python/atac/barcodes/__init__.py:66–78` preserves file order when requested; RNA loader `lib/python/cellranger/barcodes/utils.py:46–52,98–106` preserves order and rejects duplicated RNA entries. Sorting either list, converting to a set or treating the two barcode strings as identical would change this operation.

The installed GEX `ARC-v1` chemistry specifies the first 16 bases of R1 as the gel-bead barcode and the next 12 as UMI (`chemistry_defs.json:286–310`). This is the **GEX** R1 definition, not a raw 24-base ATAC index rule. The installed joint-library compatibility stage actually invokes the translator (`check_joint_barcodes_compatibility/__init__.py:88–107`). Presence of that code does not prove that these particular submitted libraries have passed its compatibility check.

Both canonical compressed whitelist files are present and fingerprinted below; neither was decompressed here. Complete whitelist lengths, uniqueness, bijection, order correspondence and historical 2.0.0 compatibility therefore remain unverified by this review. In particular, the source translator's use of `zip_longest` is not by itself a check against unequal list lengths. No correction, translation or receiving barcode join was performed. Root's separate 100-spot ATAC pilot establishes a read-layout diagnostic, not these remaining operations or a full-run orientation/correction decision.

## RNA-only lineage assignment is not yet a compatible crosswalk

The primary paper's annotation definitions are **cHep, pHep, HSC, EC, macrophage, lymphocyte and cholangiocyte**. Its RNA methods use groupwise normalization, Harmony and clustering, rather than a deposited frozen seven-lineage classifier. Source QC also requires both RNA and ATAC quality. These statements concern processing and annotation, not clinical outcomes. [Primary methods and annotation definitions](https://pmc.ncbi.nlm.nih.gov/articles/PMC11103255/).

| Author category | Frozen development category | Source-only qualification |
|---|---|---|
| cHep and pHep | Hepatocytes | Combining both hepatocyte labels is biologically plausible; exact receiving per-barcode labels are still required. |
| cholangiocyte | Cholangiocytes | Compatible broad name; no actual per-barcode crosswalk established here. |
| HSC | Mesenchymal | HSC may be narrower than the native mesenchymal pool; equivalence cannot be assumed. |
| EC | LSEC | EC does not certify liver-sinusoidal-only membership. |
| macrophage | Kupffer | Broad macrophages do not certify resident-Kupffer-only membership. |
| lymphocyte | B cells and NK-T | A single author lymphocyte group cannot directly produce two required lineages. |

Thus seven author labels do not yield a one-to-one seven-class receiving annotation. The already inspected GEO inventory names matrices, features and barcodes, with no identified per-barcode lineage table. The primary page links `NIHMS1990022-supplement-Supplementary_Material.docx` and `NIHMS1990022-supplement-MDAR.docx`. The bounded official Supplementary Material retrieval returned a 1,817-byte non-ZIP response, so its content was not inspected. **[Unverified]** A supplement-held cell-label mapping remains unresolved; this review does not assert that none exists. The primary HTML exposed no dataset-analysis code repository link; a generic NCBI GitHub footer is not author analysis code.

An existing RNA-only mechanism is specified in `Analysis/SingleCell/continual_integration/reference/router_v2_policy_v43.json`: frozen Healthy_Human_Liver CellTypist model, integer UMI counts, 10,000-library log1p normalization, at least 2,400 of 2,691 model features, query-only majority voting and seed 17, without author condition or integration-coordinate inputs. Its policy collapses B/NK-T/hepatocyte/cholangiocyte groups but uses broader Endothelial/Fibroblast/Macrophage/Other immune labels. It is an available prespecified RNA-only mechanism, **not** an admitted native-author LSEC/Kupffer/Mesenchymal crosswalk. Neither its model body nor annotation/accuracy outputs were read or rerun. Reference-gene compatibility and receiving label definitions must be established before adopting any such route.

The paper specifies GRCh38 and ARC 2.0.0, but the inspected methods do not give the exact GTF/reference-package identity needed to assert equivalence with the fixed raw gene axis. Installed 2.0.2 alone does not resolve that gap. Reference absence, a filtered feature row and measured zero UMI are different states; no zero-filling or receiving normalization choice is authorized by this review.

## Receipts, completed checks and limits

Public API responses were inspected in memory and not saved as new payloads. Dynamic ENA/NCBI responses changed byte ordering/hash on refetch while the selected exact run identities, sizes and MD5 values reconciled. The fingerprints below identify specific responses, not immutable accession contents. Existing metadata/code files were hashed without decompressing canonical barcodes or source feature members.

| Source or metadata artifact | Bytes | SHA256 |
|---|---:|---|
| Existing GEO identity JSON | 10,741 | `172aa3c42134a2a4f329cdf24db937d92c70da0a1300827f32c6cf13b3754018` |
| New ENA response, selected GEX rows inspected | 9,169 | `1721fb24dd1388b0268b7913e1e5f1c2a5f7a8cdcb72db37f7948fde059daa76` |
| NCBI six-experiment search JSON | 919 | `cb8464f9ac66bbba3b32871bd7936f1599caf286592d1875327dad7c2a1d7696` |
| NCBI experiment/run XML, final retrieval | 99,118 | `3eb667c963abfefaabaa43ebcf0901c7b76084147b82b6f9df150af1855cc0ae` |
| Existing first-line RNA schema JSON | 1,471 | `ea7eae8f9d96ad1a0781dfd32bb9086f80d5cabe5a1121c0ec8aec5cc2e5f6f8` |
| Installed ARC `atac_rna/feature_ref.py` | 7,056 | `1b80c8ce8675e1fe0c079f0167fc4209a84712f0bcf28162465251ac33c2ac51` |
| Installed ARC `atac_rna/barcodes.py` | 2,701 | `9f4e0c90eec48ac8bbadd0bd69377c8fca9b1a17cac59b216283dd50fc82f185` |
| Installed ARC ATAC `737K-arc-v1.txt.gz` | 2,465,078 | `9a4c8f29a05f59b9a4d7525a2f44d1778ccfa540a8832ad9f5f84047f374cc26` |
| Installed ARC RNA `737K-arc-v1.txt.gz` | 2,142,553 | `1279bb7ba28ddcd805e2f16ebd5d6933b7b1e7f128e18aabfe6ae8cd7b7bf710` |
| Existing RNA-only Router V2 policy | 4,209 | `f907c6a472e4d84aa0601f19652c8879355e82468f845cb9f7a1926371db055c` |
| Supplement endpoint non-ZIP response, not document content | 1,817 | `264c03e6413cce34fa85d7ac296ff0ed4810180a4a71c6e439cb2832396a51cc` |

Completed: metadata experiment/BioSample joins, twelve run/file inventories and byte arithmetic, six-column official/installed schema review, canonical translation source review, annotation-taxonomy comparison and existing policy inspection. The official schema page was successfully retrieved in the browser; no HTML copy was saved. Europe PMC primary XML retrieval failed with HTTP500; the primary PMC page supplied methods. Supplement retrieval did not yield a DOCX, so no supplemental clinical tables or cell-label values were opened.

Not run: sequence export/integrity/alignment, full feature CRC or ordered-ID census, canonical whitelist decoding/translation, matrix quantification, label assignment, person-identity verification or model fitting/evaluation. Biological n remains six author-described people; technical lanes, nuclei and donor pairs cannot increase it. Any eventual admitted result would be small-n pediatric source/assay transport, not adult MASLD clinical validation. External admission remains unchanged.

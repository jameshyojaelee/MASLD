# SRR23252452 bounded read-layout pilot preparation

Prepared 2026-10-01 for parent review. **No raw download, export, validity computation or submission has run.** This is one reviewed public ATAC run and the first 100 spots only. It does not prepare alignment, recounting, RNA acquisition, lineage annotation, model fitting or evaluation. The [raw-source feasibility proposal](2026-10-01_codex_gse223843_raw_fragment_feasibility.md) and original geometry helpers/results remain unchanged.

## Frozen files and allocation

| New owned artifact | Bytes | SHA256 |
|---|---:|---|
| `scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse223843_read_layout_pilot.py` | 27,924 | `53771e445f3457e134a5727ded259f83e2f44ffb987c84c5d2f3f0c7e31a100f` |
| `scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_gse223843_read_layout_pilot.sbatch` | 993 | `c6ad3595806d96b96f8419e03675591dbc53cd79a8d38ba8bd2c6517a71a01da` |

Both are executable. The launcher pins the helper hash. Proposed one individual job: name 1592, `cpu` partition, `nslab` account/QOS, **1 CPU, 2 GB, 1 hour**. Fresh output is `Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/codex_gse223843_read_layout_pilot_<actual_SLURM_jobid>`; pre-existing directories are refused. Output/error logs are `1592_raw_layout_%j.out/.err`, distinct from prior 1592 jobs. No job index was edited and no job was submitted.

Static checks completed on login: Python AST parse, `bash -n`, executable permissions and exact launcher/helper hash agreement. Runtime configuration, source download/checksum, tool export, FASTQ parsing and whitelist validity checks are **not run**. No extra executable tests were run on login.

## Acquisition and local export

The sole URL is [SRR23252452 normalized SRA](https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR23252452/SRR23252452). Expected bytes **5,191,054,225** and MD5 **3cdda9ebd3fcd6c156307fab5b98bbe6** are pinned against the saved official metadata. The identity must be Fontan1/GSM6997758/SRX19197939/SAMN32937884. All saved metadata/proposal hashes are checked before access and rechecked at completion.

The helper uses one exact HTTPS GET, refuses redirects/fallback/proxy routes, requires status 200, exact final URL, matching Content-Length and identity Content-Encoding before body iteration. It streams SHA256/MD5, records received body bytes, retains a `.partial` file on errors, and only renames its own verified file. Disk checksum is independently verified before export and after parsing. There is no archive, ENA two-FASTQ or original-cloud fallback.

Configuration is a fresh own `isolated.mkfg`, with remote/user/site repositories disabled. `NCBI_SETTINGS` points only to that file; alternative inherited NCBI configuration search environment variables are removed. Three read-only `vdb-config` queries must each confirm disabled=true before acquisition. The queried XML leaf parser follows [vdb-config 3.2.1 source lines 1594–1597](https://github.com/ncbi/sra-tools/blob/3.2.1/tools/external/vdb-config/vdb-config.c#L1594): a leaf is printed by its final name `disabled`, not by an entire repository tree. The file hash is rechecked at completion. No user configuration is written.

The exact export executable is `/nfs/sw/easybuild/software/sratoolkit/3.2.1/bin/fastq-dump-orig.3.2.1`, SHA256 `9f89c4a6af231c09c180b40de5b774e397a659dd79c6ca71f9b2cfbf9a7915ba`; it bypasses the accession-discovery dispatcher and receives an existing absolute local file. Installed version/help flags are checked on compute. Its argv is:

```text
fastq-dump-orig.3.2.1 --split-files --minSpotId 1 --maxSpotId 100
  --disable-multithreading --defline-seq '@$ac.$si/$ri' --defline-qual '+'
  --offset 33 --ncbi_error_report never --outdir <fresh_output>/fastq
  <fresh_output>/download/SRR23252452
```

No technical-read skipping, three-way biological-only splitting, read filtering, trimming or alignment is used. Actual `$ri` deflines must establish four read indices 1–4, lengths 8/50/24/49 and all ordered spot IDs 1–100 in each file. Every sequence/quality pair must have equal declared length, narrow A/C/G/T/N sequence syntax and printable Phred+33 quality. Read roles remain explicitly source-declared; genomic mapping is not inferred from a 100-spot format check.

## Caps and content suppression

The helper enforces **6,000,000,000 bytes** network body and inclusive saved output, **1 MiB** aggregate exported FASTQ, and **32 MiB** extra metadata/logs. Its 55-minute overall clock leaves five minutes within requested wall time; download has a 45-minute bound, export five minutes and metadata tools 30 seconds each. Tool subprocesses are checked every selector interval (50 ms) for time/disk limits and terminated as their own process groups on failure. Each export file also has a hard 256 KiB `RLIMIT_FSIZE`; only the four exact filenames are allowed. Source download writes are budget-checked before each block. Output symlinks are refused. The monitor includes the actual job's SLURM output/error files in saved and metadata accounting.

All tool stdout/stderr are streamed into byte counts and SHA256 digests, then discarded. Only version/help/configuration stdout is temporarily held for narrow verification, not serialized. Raw export diagnostics are never saved or forwarded. Arbitrary exception messages and tracebacks are suppressed; fixed failure codes, phase, byte counts, hashes and partial-file retention are reported. Python/launcher core dumps are disabled for child tools and launcher. Raw sequence/quality strings never enter JSON, SLURM logs or console output. The local SRA and FASTQs are preserved as acquired payloads; parent review must inspect only JSON metadata, not those files.

## Canonical ATAC whitelist and four separate hypotheses

Installed ARC version is 2.0.2, not the author's exact 2.0.0. Its `mro/atac_rna/sc_atac_gex_counter_cs.mro:87` passes `atac_barcode_whitelist = "737K-arc-v1"`. The ATAC loader in `lib/python/atac/barcodes/__init__.py:51` resolves this name beside the **ATAC** package. Therefore the pinned file is:

```text
/nfs/sw/easybuild/software/cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz
```

Compressed size **2,465,078 bytes**, SHA256 **9a4c8f29a05f59b9a4d7525a2f44d1778ccfa540a8832ad9f5f84047f374cc26**. MRO/ATAC loader/translation source, ARC executable and compressed reference hashes are pinned before/after. Only source/path/stat/hash inspection was done on login; no whitelist decompression, barcode-value inspection or validity computation ran there.

On compute the helper must fully decode this canonical file with CRC, 32 MiB/one-million-line bounds, narrow unique 16-base A/C/G/T entries and the source loader's explicit `#` comment rule. It reports reference identity/entry count/hash, never its barcode strings. From the actual 24-base read-index-three file it retains four independent exact-whitelist counts/fractions: first16 forward, first16 reverse complement, last16 forward, last16 reverse complement. It does not maximize/select one, correct mismatches, translate to RNA or use model output. This follows the layout alternatives in [official ARC 2.0 barcode processing](https://www.10xgenomics.com/support/software/cell-ranger-arc/2.0/algorithms-overview/algorithms-overview), while deliberately stopping before the canonical full-data choice/correction. One hundred spots do not establish global orientation or barcode-abundance correction probabilities.

## Intended receipt and remaining limits

Successful `summary.json` records source/run identity, exact argv, installed versions/help evidence, pre/post source hashes, isolated configuration, received bytes/source checksums, four FASTQ file checksums/record counts/lengths/read indices, synchronized spot IDs, quality-format/ambiguous-base counts, four validity receipts and measured elapsed time/byte maxima. `pre_download_guards.json`, `download_started.json` and `export_started.json` preserve small milestones. A failed attempt writes its own `failure.json` when budgets permit, otherwise the fixed failure phase/code reaches the small job log; partials are retained and never reused or deleted.

This pilot, even if completed, does not establish a genuine complete 306,706-interval count axis, seven homologous source lineages, compatible receiving RNA, donor independence, biological absence, model correctness or evaluation admission. It neither alters the original failed geometry job nor widens the seven unsupported file schemas. Full native fragment processing still requires a reviewed reference/ARC route, counting unit/interval convention, corrected ATAC→GEX barcode translation and source-owned lineage mapping. All those operations remain unrun.

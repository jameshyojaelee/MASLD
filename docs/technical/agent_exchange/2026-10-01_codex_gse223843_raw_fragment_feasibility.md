# GSE223843: public raw ATAC read feasibility

Prepared 2026-10-01. This is source/file/read-layout metadata and a bounded pilot proposal, not a completed fragment producer or evaluation. Original geometry sources, helpers and jobs 21999459/21999912 are preserved.

## Finding

All six author ATAC persons have public normalized SRA files whose metadata preserves four stored read classes of 8, 50, 24 and 49 bases per spot. Twelve exact source-derived HTTPS URLs returned HTTP 200 to HEAD, with Content-Length equal to the declared normalized SRA file size. Zero sequencing body bytes were read. This provides a concrete potential source for fresh measurement on the unchanged 306,706 intervals. [Unverified] Actual export synchronization, technical-read preservation, barcode orientation/correction/translation and exact fragment-count compatibility require a compute pilot before any full producer.

The already acquired processed peak matrices cannot supply this full target: five qualified source axes have fewer than 306,706 rows and at most four matching coordinates under the separately reported hypotheses. The seventh unsupported feature-file result is preserved; missing target coordinates are not biological absence or zero counts. See the [completed geometry report](2026-10-01_codex_gse223843_feature_barcode_proposal.md).

## Identity and available files

Official records are [GSE223843](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223843), [PRJNA928662](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA928662), the experiment pages below, and the exact ENA/NCBI metadata URLs saved in the receipt. Names below are author sample identifiers, not a new clinical classification. Each pair of runs is technical sequencing of one source person, not two independent persons.

| Author person | GEO / BioSample | SRA experiment | Two ATAC runs | Spots | Normalized SRA bytes |
|---|---|---|---|---:|---:|
| Control1 | GSM6997754 / SAMN32937886 | [SRX19197937](https://www.ncbi.nlm.nih.gov/sra/SRX19197937) | SRR23252455; SRR23252456 | 553,822,808 | 23,878,574,581 |
| Control3 | GSM6997756 / SAMN32937885 | [SRX19197938](https://www.ncbi.nlm.nih.gov/sra/SRX19197938) | SRR23252453; SRR23252454 | 439,031,744 | 20,791,525,260 |
| Fontan1 | GSM6997758 / SAMN32937884 | [SRX19197939](https://www.ncbi.nlm.nih.gov/sra/SRX19197939) | SRR23252451; SRR23252452 | 238,193,482 | 10,464,040,724 |
| Fontan2 | GSM6997760 / SAMN32937883 | [SRX19197940](https://www.ncbi.nlm.nih.gov/sra/SRX19197940) | SRR23252449; SRR23252450 | 434,989,544 | 18,460,689,163 |
| Fontan3 | GSM6997762 / SAMN32937882 | [SRX19197941](https://www.ncbi.nlm.nih.gov/sra/SRX19197941) | SRR23252447; SRR23252448 | 390,034,964 | 18,230,698,781 |
| Fontan4 | GSM6997764 / SAMN32937881 | [SRX19197942](https://www.ncbi.nlm.nih.gov/sra/SRX19197942) | SRR23252460; SRR23252461 | 248,647,935 | 11,768,920,753 |
| Total, six source persons | | | 12 technical runs | 2,304,720,477 | 103,594,449,262 |

Every run's selected original file list contains I1/R1/R2/R3, with exact paths, declared bytes and MD5 values in `identity_read_layout_file_metadata.json`. These original files are marked Cloud Data Delivery; anonymous original-FASTQ retrieval was not established. The normalized files explicitly declare anonymous/worldwide access. HEAD establishes public file availability and exact size, not downloaded content integrity.

The original four FASTQs total 205,547,166,301 compressed bytes. ENA reports only two FASTQs per run, totaling 92,121,810,155 bytes, with base count exactly 73 times the spot count. Full SRA stored-read statistics sum to 131 bases per spot. [Inference] ENA's 73-base total is consistent with a 24+49 subset, rather than the two 50+49 genomic mates. The actual mapping of each ENA FASTQ to the four read classes is unverified. Those two files must not be admitted as a demonstrated paired-genomic ATAC input.

## Read roles and barcode limitations

The source methods specify Chromium Next GEM Multiome, NovaSeq 6000, ATAC cycles 50:8:24:49, Cell Ranger ARC 2.0.0/GRCh38 and scATAC-pro 1.4.2 with source-specific peak calling. Same-nucleus RNA/ATAC library construction is explicit. This is the pediatric source; it does not supply adult validation. These claims use only the narrow library/preprocessing Methods section of the [primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11103255/).

Official ARC [FASTQ naming](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/inputs/specifying-input-fastq-count) identifies original R1 as genomic read 1, R2 as i5 containing the barcode, R3 as genomic read 2, and I1 as the optional i7 sample index. The 24-base i5 contains a 16-base ATAC barcode and eight-base spacer. The [10x CG000373 Rev C technical note](https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000373_TechNote_SeqMetrics_MultiomeATACLibrary_RevC.pdf), printed pages 2–3, describes both forward and reverse-complement NovaSeq workflows with different kit/software versions. Forward places the barcode before the spacer; reverse-complement places the spacer before the barcode. Source metadata specifies neither the required kit nor instrument software version, so “NovaSeq 6000, 24 bases” does not choose the actual orientation.

ARC 2.0 documents comparing first/last 16 bases in forward/reverse orientation by allowed-whitelist validity, then Hamming-one correction using base quality and barcode abundance with a greater-than-90% posterior requirement. These are documented measurement operations, not model-performance selection. A bounded pilot must retain all four separate validity receipts; it must not silently select an orientation, create a correction rule from the six processed barcode lists, or use model scores. ARC additionally uses custom trimming, modified BWA-MEM with insert prior mean 250/SD 150, barcode-aware duplicate marking, both mates MAPQ >30 and fragment/Tn5 filtering. A Bowtie2 command with similar output columns would not establish identical native units. [Official ARC 2.0 algorithms, barcode and fragment sections](https://www.10xgenomics.com/support/software/cell-ranger-arc/2.0/algorithms-overview/algorithms-overview).

Raw ATAC and GEX barcode sequences on a gel bead differ. ARC's corrected CB is translated to the GEX barcode namespace; CR/CY retain raw ATAC sequence/quality. The acquired processed barcode joins do not establish a raw 24-base prefix-to-RNA mapping. [Official ATAC barcoded BAM definition](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/tutorials/outputs/cr-arc-atac-barcoded-bam).

## Smallest bounded export pilot, not run

Smallest complete normalized file: **SRR23252452**, Fontan1/GSM6997758/SRX19197939, second technical lane, **117,710,245 spots**.

- Exact source URL: [public normalized SRR23252452](https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR23252452/SRR23252452).
- Declared file size: **5,191,054,225 bytes**. MD5: **3cdda9ebd3fcd6c156307fab5b98bbe6**.
- Four stored-read indices 0–3: average lengths 8/50/24/49, stdev zero, each count 117,710,245.
- Original prefix `ATAC_Fontan1_S2_L002_`, suffixes `I1`, `R1`, `R2`, `R3`, each `_001.fastq.gz`.

The XML RUN attribute size 5,191,051,404 differs from the normalized SRAFile size. Download guarding must use the exact SRAFile size above, which the HEAD response confirms, and the source MD5. A future helper must download once to a fresh compute-owned directory with a strict byte guard and no overwrite, retain failures, and require MD5 before export. A tiny spot export still requires this 5.191 GB file. No partial remote/SRA-index route is proposed.

Installed SRA-Toolkit 3.2.1 help confirms `fastq-dump --minSpotId`, `--maxSpotId`, `--split-files` and `--outdir`. `--skip-technical` discards technical reads, so it must be omitted. `--split-3` also discards all other reads and is unsuitable. The bounded command template is:

```bash
module load SRA-Toolkit/3.2.1
NCBI_SETTINGS=/absolute/fresh_pilot/isolated.mkfg \
  /nfs/sw/easybuild/software/sratoolkit/3.2.1/bin/fastq-dump-orig.3.2.1 \
  --split-files --minSpotId 1 --maxSpotId 100 \
  --outdir /absolute/fresh_pilot/first100 \
  /absolute/fresh_pilot/SRR23252452
```

This is **not run**. The direct installed export executable avoids the accession-discovery dispatcher; the input is a verified existing absolute local file, never an accession. Before execution, a separate pilot helper must create a fresh isolated settings file with `/repository/remote/disabled = "true"`, `/repository/user/disabled = "true"` and `/repository/site/disabled = "true"`; confirm effective settings using `vdb-config` with that `NCBI_SETTINGS`; and fail if remote access is enabled. These settings must not alter the user's configuration. The [NCBI VDB 3.2.1 configuration source](https://github.com/ncbi/ncbi-vdb/blob/3.2.1/libs/kfg/config.c#L2008) reads `NCBI_SETTINGS`; [repository access properties](https://github.com/ncbi/ncbi-vdb/blob/3.2.1/libs/kfg/properties.c#L294) define the disabled keys. Installed README-vdb-config explicitly warns that enabled remote access can download on the fly. Effective isolated configuration and actual offline behavior remain untested.

Proposed scope for parent review: one individual `cpu`/`nslab` compute job, 1 CPU, 2 GB RAM, 1 hour; saved/network cap **6,000,000,000 bytes**, exported FASTQ aggregate cap **1 MiB**, small metadata/log cap **32 MiB**. No alignment, correction, RNA export, matrix production or model execution. Checks: four read-class files, 100 synchronized spot identities, per-class lengths and sequence/quality agreement; exact export argv/binary/config hashes; canonical installed ATAC-whitelist identity before any validity computation; all four validity fractions separately, with no chosen origin and no sequence strings in reports. One hundred spots establish export format, not reliable global orientation or correction abundance. Full read-layout acquisition/compute preparation has not been authorized by this metadata task.

## Installed producer and full-source resource limits

`cellranger-arc/2.0.2` is already installed and `cellranger-arc --version` reports `cellranger-arc-2.0.2`. Actual prefix is `/nfs/sw/easybuild/software/cellranger-arc-2.0.2/`; the module sets `CELLRANGER_ARC_MAX_CELLS=50000`. This is not the source's exact 2.0.0. Official ARC count requires ATAC **and** GEX; an ATAC-only invocation is not supported by the documented input contract. Raw GEX acquisition has not been inventoried or admitted here. Therefore the pilot only qualifies stored-read extraction. An eventual native fragment producer still needs a reviewed ARC/reference route, barcode translation and cell-selection policy, and count-unit agreement with the fixed producer.

[Estimate] NCBI describes approximately 7× SRA for exported FASTQ, scratch about 1.5× FASTQ, and about 17× SRA total conversion space. Thus all 103.594 GB simultaneously could require **about 1.76 TB before alignment**, exceeding the 1.5 TB campaign context. Largest run SRR23252455 is 11,948,361,350 bytes: about 203 GB conversion space alone. Serial bounded runs may reduce temporary overlap, but retained fragments/references and GEX would still need a separate aggregate allocation. These factors are planning estimates, not observed peak usage. [NCBI conversion guidance](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump).

[Estimate] Full alignment of 2.305 billion spots is substantial; a possible per-person planning envelope is 8–16 CPUs/64–128 GB/90 hours, subject to a real benchmark and producer requirements. No full-run CPU or memory performance has been measured, and no allocation or submission is proposed here.

## Remaining inputs and interpretation

A common-coordinate fragment recount alone does not establish the seven native author-lineage receiving task. Still missing: a source-owned barcode-to-lineage identity artifact for these persons; explicit homologous mapping to the fixed seven training lineages, including any hepatocyte subdivision; compatible RNA feature identity/measurement and source processing; verified corrected ATAC→GEX barcode translation/gem-group names; fragment versus transposition-event counting rule, duplicate multiplicity treatment and interval-overlap convention matching the fixed producer; exact reference/contig compatibility and real complete coordinate rows. The six acquired RNA feature files remain unsupported under the frozen identity parser, and no further parser retry is part of this task.

Six persons are the maximum source biological n; two runs per person do not increase n. [Unverified] Person disjointness from GSE296875 and other sources is not established. Public release dates and sequencing-file upload dates are not specimen-acquisition dates. Different institutions do not prove disjoint participants. Nothing here admits receiving counts, native lineages, outcomes, model evaluation or source independence.

## Bounded evidence and checks

Paths below are relative to the repository root. `REC` means `Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z`.

| Saved source metadata | Bytes | SHA256 |
|---|---:|---|
| `REC/codex_gse223843_identity_metadata/identity_metadata.json` | Existing selected GEO identity receipt | `172aa3c42134a2a4f329cdf24db937d92c70da0a1300827f32c6cf13b3754018` |
| `REC/codex_gse223843_raw_read_metadata_20261001/identity_read_layout_file_metadata.json` | 79,232 | `daf440a0425678ccc597f1fe0e2e898949d6cd80266fb5e959c61ecf4fb01b34` |
| `REC/codex_gse223843_raw_read_metadata_20261001/ENA_ATAC_file_metadata.json` | 11,927 | `d4b48f024808f7ba5cefe963983abadb2c8918ad7be32e641ca88967311ca667` |
| `REC/codex_gse223843_raw_read_metadata_20261001/normalized_SRA_HEAD_metadata.json` | 4,815 | `418fcffed5263f8e86e803111bbd7df0bcb9806e1ee3b5e47dfe36be2a255a31` |

The selected receipt stores exact source URLs, source-response hashes, all 12 ATAC identities, all file paths/MD5/bytes, read-layout statistics and HEAD receipts. The full ENA response was not saved; its ordering/hash changed on refetch. Selected rows were sorted by run and checked against exact expected experiments/BioSamples/files. This does not imply identical raw ENA response bytes. Full SRA metadata response SHA256 was `85014df5950e1d5b21e56615eae5910b24a92357104ad3398bc16b13a808a0c2`; only selected identity/file fields were saved.

Installed files checked, not executed on study sequence:

- SRA `fastq-dump`/`fasterq-dump` dispatcher SHA256 `822a39c0f715e50cdc3548138bfabb5691d02d4ee151944082da4f3755b21b56`.
- Direct `fastq-dump-orig.3.2.1` SHA256 `9f89c4a6af231c09c180b40de5b774e397a659dd79c6ca71f9b2cfbf9a7915ba`.
- ARC `bin/cellranger-arc` SHA256 `61e9fd4936c3c5bf24ca78d7f2d777ea6bafa2ec9241c265fe8e1801ca412f72`.
- ARC `lib/python/atac_rna/barcodes.py` SHA256 `9f4e0c90eec48ac8bbadd0bd69377c8fca9b1a17cac59b216283dd50fc82f185`; translation uses ordered RNA/ATAC whitelist pairing and preserves gem groups. Whitelist contents were not inspected.
- ARC `lib/python/cellranger/chemistry_defs.json` SHA256 `5fb5b204d8900362a23124c72e3e66c2ad0feb28db8535995082be3d08f637a5`; its ARC-v1 RNA barcode offset is a GEX definition and is not evidence for slicing the ATAC 24-base index.

Completed: official identity/file metadata extraction; twelve bounded HEAD checks; selected metadata inspection/hashes; installed version/help and small source/configuration inspection; primary narrow methods and chemistry documentation. Not run: raw download, MD5 validation of sequence, spot export, effective isolated offline configuration, whitelist/orientation validation, alignment, fragment recount, lineage mapping, model fit/scoring or SLURM submission. No sequencing body, count matrix, clinical cell table or genotype payload was acquired or inspected.

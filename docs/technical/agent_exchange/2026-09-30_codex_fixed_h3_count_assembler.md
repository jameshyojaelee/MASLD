# Fixed-region H3 count assembler preparation

**Synthetic checks passed:** a bounded callable now assembles complete native featureCounts single-read output from the202 source-owned H3 runs into the17 fixed donor rows. Root reviewed the separate synthetic checker and submitted job **21998841 (1579)**, cpu/nslab, one CPU,2GB,one hour. It completed in five seconds, exit0, empty stderr. No receiving counts, reads, alignment files, clinical values, prediction arrays or targets were opened.

## Source and measurement contract

The helper hard-guards the actual `paired_h3k27ac_raw_roster_21997742` summary, `raw_runs.tsv`, `raw_files.tsv` and `source_samples.tsv`, plus the completed fixed SAF/crosswalk/summary from21998535. It uses those files' **actual** public columns: `run_accession`, `donor`, `GSM`, `experiment_accession`, `source_biosample`, `sample_accession`, `assay`, `library_layout`, `library_strategy`, `fastq_file_count` and `file_number`. Source sample/run/file joins must agree; all202 H3 runs must be SINGLE ChIP-Seq with one native FASTQ file. RNA and Input runs are excluded by the fixed source roster, not by receiving count values.

The fixed donor order is B1,B7,B8,B9,B10,B12,B15,B19,B21,B22,B24,B26,B36,B38,B41,B46,B47. Technical runs sum within these17 source-labelled persons; runs, lanes, bracket tokens and feature regions do not increase biological n. These joins do not independently verify private donor identities or PCR-library ownership.

Every count table must contain all **96,460 fixed native region rows** in the exact SAF order. Native featureCounts columns `Geneid`, `Chr`, `Start`, `End`, `Strand` and `Length` must agree exactly with the SAF and crosswalk; Length is end−start+1, not an inferred interval or count. The SAF dummy `+` strand is ignored only through explicit `-s0`. Zero count rows are retained when emitted in complete native tables; absent rows are refused rather than invented or zero-filled.

## Explicit file/run/command binding

The new `CountTable` input structure is an assembler interface, not an invented source column schema. It binds each exact featureCounts BAM column string to the five existing public fields `run_accession`, `donor`, `GSM`, `experiment_accession`, `source_biosample`. Every binding must equal the frozen roster entry; BAM basenames must contain exactly their assigned native SRR accession. Multiple run columns in one native table or separately counted per-run tables are supported. Repeated files, duplicate runs, unknown runs, missing runs and wrong source assignments are refused before count parsing.

Each table declares the actual featureCounts argv and each BAM's actual prefilter argv. The native `# Program:featureCounts v2.1.1; Command:` header must reproduce the declared counting argv exactly. Counting requires SAF, the fixed SAF path, explicit unstranded and primary counting, the declared output path and input BAM columns in their native order. Only conventional thread/output/annotation/strand/MAPQ/primary arguments are accepted; `-O`, `-M`, fraction assignment, paired-fragment flags, duplicate removal and changed overlap rules are refused. A declared count MAPQ option must be10; absent count MAPQ filtering relies on the separately required MAPQ10 BAM prefilter.

The declared samtools-view prefilter must use `-q10` and `-F2308` (unmapped, secondary and supplementary exclusion), with its output matching the counted BAM and one same-run input BAM. This preserves the already specified receiving unit: **one primary nonsupplementary mapped single read overlapping one fixed region**, unstranded, MAPQ≥10, without multi-region/fraction assignment or PCR deduplication. No receiving peak discovery, Input subtraction, matched-library assumption or CUT&RUN fragment reconstruction is introduced.

**Provenance limit:** a printed command and correct public mapping do not prove actual BAM contents or command execution. The authorized upstream alignment/counting workflow must establish those facts with its real receipts. The assembler records this limit explicitly; it never reads an alignment or raw FASTQ to make that stronger claim.

## Numerical limits and target interface

The callable streams count rows, holding only the17×96,460 `uint64` aggregate and bounded public metadata. It refuses an input-table population above8GB and any line above1MB. Each count token must be a native nonnegative integer; fractions, negative values, exponent notation and nonfinite values are rejected. Addition and donor totals use Python integers for checking before uint64 assignment. Every partial cell and complete donor total must remain≤2^53, with every final donor total positive. This prevents both unsigned overflow and loss of exact integers in the later float64 target transformation.

[Estimate] the13,118,560-byte aggregate and public SAF/crosswalk/run metadata should fit comfortably within2GB; actual peak memory is not measured. The result supplies counts, exact donor IDs, exact region keys and input command/hash evidence. It does not normalize counts, apply B99, evaluate predictions or select donors.

This directly matches the input shape and integer range of [apply_fixed_h3_residual_target.py](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/apply_fixed_h3_residual_target.py). Its separate synthetic check21998710 is already complete: saved uniform-scaling maximum difference2.842170943040401e−14, hand-derived uniform-row difference2.0747847884194925e−12 and eight invalid-input refusals. That prior result does not verify this new assembler or actual H3 measurements.

## Access boundary and prepared checks

The helper has **no receiving-count CLI enabled**. Its production callable requires compute and guards the complete fixed public metadata before and after assembly. Root's separately authorized calling workflow must first verify RNA measurement and frozen prediction prerequisites. No future prerequisite receipt schema or new permission condition is invented here. The synthetic checker is the only prepared executable entry point.

The checker reads permitted public metadata and creates synthetic counts on all202 actual public run identifiers with the17 fixed donor assignments and three fixed coordinate metadata rows. Production metadata loading still requires all96,460 regions; the private synthetic core uses three rows to keep the check small. Independently hand-derived exact donor sums, technical-file partition/order invariance, retained native all-zero region rows, an exact2^53 donor-total boundary and38 invalid source/region/count/command refusals **passed** in21998841. Its output and launcher stdout/stderr paths include the actual Slurm job ID; earlier1579 artifacts are preserved.

The bounded **actual-tool synthetic smoke passed** in21998841. It plants five all-A20M records in a tiny SAM, without running a genome aligner or using study sequences: one ordinary primary record, one duplicate-flagged primary record, one MAPQ9 record, one secondary record and one supplementary record. Actual installed samtools1.22.1 applies `view -b -q10 -F2308`, and actual featureCounts2.1.1 counts an explicit three-row toy SAF unstranded with the declared fixed-unit options. The same private assembler parsed the native emitted command/header/coordinate/length/path fields and reproduced the hand-derived counts **[2,0,0]**; the duplicate flag remains counted because PCR deduplication is not part of this recipe. The filtered BAM's planted record names/flags are checked independently. This tests tool output compatibility and basic single-read selection, without claiming whole-study counting or biological validity. Actual binary hashes, printed versions, argv, native output header and synthetic file hashes are saved in its receipt. The toy run label comes from public schema metadata; the planted records are not that study's reads. No parser format relaxation was needed. SLURM reported zero RSS, which is unusable as a memory measurement. The executed checker SHA256 is `749347dc2a6e3150b4ba3a4a0b126b7efc48d9fac3a629d3f0adcd65ab92fcb7`, and the unchanged assembler SHA256 is `0fb8317a9ffa4ccebd3e621915d9efbe72bdf93a8753fd5c4e983631b114d81f`. Numerical evidence is in `fixed_h3_assembly_checks_21998841/summary.json` under the REC execution.

Files created:

- [Assembler](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/assemble_fixed_h3_counts.py)
- [Synthetic checker](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/check_assemble_fixed_h3_counts.py)
- [Prepared launcher](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/run_check_assemble_fixed_h3_counts.sbatch)
- This distinct technical report.

Completed checks: Python AST/compilation, launcher shell syntax and scoped whitespace inspection. No previous files, job index, ongoing aligner/reference work or paper authorities were changed. Actual count assembly, residual-target application, paired predictive accuracy and biological validation are **not run by this preparation**.

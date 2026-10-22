# Independent H3 genome-reference review

**Outcome:** the present producer targets the exact verified Ensembl98 primary genome, with an appropriate all-contig identity check. The submitted index job **21998702 (1594)** is running; completed index sequence agreement is **not established yet**. This review created only this report. No reference sequences, patient sequencing, expression/chromatin arrays or outcomes were read; no script, launcher, job index or running job was amended.

## Present source and submitted job

Two inspections during this review yielded unchanged present-file hashes:

| File | SHA256 |
|---|---|
|`codex_build_h3_genome_index.py`|`c8516a6f3beb38bfdd12f7534958fecd92fb229956e3ee3a62ede95313774a3b`|
|`run_codex_build_h3_genome_index.sbatch`|`57bdb0637b9f45d8e1db00eaafc30c1a2904f3adaac2cd75267526b8926577cf`|

These are **present source hashes**, not proof of the exact source loaded by an already running interpreter. The job-specific protocol and summary were absent at the initial00:45 UTC inspection. The source writes its protocol after full reference checks and records a runtime source hash, but does not save an executed source copy or hashes of imported helpers. This does not show a sequence error; the intended complete reconstructed-sequence comparison is stronger evidence of the resulting genome identity than filenames or source hashes alone.

The own `job_index.tsv:102` records21998702 as1594. Direct SLURM inspection at00:45 UTC confirms RUNNING on `ne1dc4-043`; a second inspection approximately two minutes later confirms the same live job. The launcher requests one cpu/nslab job, **8 CPUs,32GB,90 hours**, and loads `Bowtie2/2.5.4-linux-x86_64` before execution. No array or new submission was made by this review.

| Job | Completed identity comparison | State at second inspection | Sampled host RSS | Downstream |
|---|---|---|---|---|
|21998702 (1594)|0/1 ░░░░░░░░░░|Running, elapsed3:31|3,806,912 KiB, approximately3.63GiB|Index build, all-reference inspection and saved summary remain pending|

One running, zero pending, zero completed in this specific set. Initial job stdout was empty and stderr reported the expected module version reload; no runtime failure was present in the inspected logs. Absence of a protocol while the full reference is being checked is not evidence that the process stopped. The job was not restarted or canceled.

## Exact reference and contig checks

[Producer lines21–28](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/codex_build_h3_genome_index.py:21) select `ensembl98_primary_genome_21998416/Homo_sapiens.GRCh38.dna.primary_assembly.fa`, its FAI and the completed genome summary. The frozen guards are:

- FASTA: `78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c`.
- FAI: `57f6ed6f8b07437708fff09814489cfda4b976652afabe2c9b1d54463bcba0ad`.
- Completed acquisition/geometry summary: `9f93a8291cf676bbf1bd1bc552d2097b3e0f9850a05f86431e1d93df65392f6b`.

The selected public reference is the same Ensembl98 **GRCh38.p13** primary assembly already qualified against the GTF and H3 coordinate axes. This review independently hashed the small FAI and rederived its **194 unique names and3,099,750,718 nucleotide bases**. It did not rehash the3.15GB FASTA on the login node. The job performs that full reference guard on compute at lines99–100 and rechecks it at lines158–159 after indexing.

The exact full FAI population is required by lines105–111. The producer streams the genome on compute and compares the complete sequence-length dictionary to the FAI. It therefore cannot silently retain only autosomes or the24 H3-bearing chromosome names. MT and primary-assembly unplaced/unlocalized contigs remain in the alignment reference even though the fixed H3 regions have no entries there.

After a successful build, [lines130–145](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/codex_build_h3_genome_index.py:130) require exactly one complete six-file `.bt2` or `.bt2l` family without mixed formats, parse `bowtie2-inspect -s`, reject repeated names and compare **all194 native names and lengths** to the FAI. Lines146–157 stream reconstructed FASTA from `bowtie2-inspect -a 100000` and compare the complete per-contig canonical sequence digest, length and ambiguous-base count to the exact declared input. This is a meaningful indexed-sequence check, not just a filename check. The reconstructed sequence is streamed, not saved as a second full genome.

## Bowtie2 version and inspector semantics

The installed module was inspected in a separate lightweight shell. Both actual version calls identify **2.5.4**:

- Builder wrapper `/nfs/sw/easybuild/software/Bowtie2/2.5.4-linux-x86-64/bin/bowtie2-build` invokes `bowtie2-build-s version2.5.4`.
- Inspector wrapper in the same directory invokes `bowtie2-inspect-s version2.5.4`.

The observed build stamp was2024-05-16; this is the installed binary's printed stamp, not a claim about when this project's index was built. The producer enforces builder2.5.4 at lines101–104, records wrapper paths/digests in the protocol and uses a fixed build seed20260930 at lines114–124. The launcher loads the versioned module at line10. It does not record underlying small/large executable digests or enforce the inspector version separately. This is a qualification of the recorded tool inventory, not an observed wrong genome.

Installed inspector help confirms `-s` reports sequence names/lengths, default inspection emits FASTA and `-a` sets the output line width. The canonical helper at [lines32–64](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/codex_build_h3_genome_index.py:32) uppercases input and maps every non-ACGT character to N, matching the documented reconstruction behavior. It bounds a FASTA line to1MB and aggregate sequence length to3.3Gb; the requested100,000-base inspection lines fit that bound. The [official Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-inspect-index-inspector) confirms ambiguous-character conversion and summary semantics. These are reference-index operations; this review did not run alignment or assume that seed coverage through ambiguous sequence establishes uniquely supported reads.

## Fresh output and byte allowance

[Lines95–98](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/codex_build_h3_genome_index.py:95) require a compute job, sufficient allocated CPUs and a **new output directory directly under the owned REC root**. Line112 creates it without overwrite; existing index products cannot be reused silently. Index construction is launched in its own process group; a monitored limit failure terminates that owned process and preserves partial products. No cleanup command deletes an existing index or reference.

The declared saved/index allowance is **50,000,000,000 bytes**, sampled every30 seconds during construction at lines67–88. The imported `output_bytes` helper sums all regular-file sizes under that output directory and tolerates disappearance of the indexer's own temporary files. Lines160–167 check completed output size, hash the six index files and write the final summary.

**Bound qualification:** this is a sampled limit, not a filesystem quota. There is no stop headroom below50GB, and the final size check occurs before the summary itself is written. Therefore the code does not prove a strict instantaneous or final aggregate≤50GB bound. No observed output exceeded the allowance in this review; the initial index directory was not yet populated. If a strict final bound is needed, the remaining summary bytes must be included in the final check, and a sampled construction stop needs headroom. This review did not change the executing producer to implement that correction.

## Fixed96,460-region annotation and future measurement

The completed [fixed-H3 annotation summary](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/fixed_h3_saf_21998535/summary.json) has SHA256 `feb1a674d911725ebf344438597881f06a6c4f2e08bbacbd2383cec487130237`. Independently computed small-metadata output digests match that receipt:

- SAF: `3c0ae95792c0c67b022bddf5a702b5ac84ba04040be24030edc9b22029e3d598`,4,614,729 bytes.
- Coordinate crosswalk: `0e5342853edc26da7bde478906d1b92841ad771ffc902bd447bbc9e5d91facc6`,7,190,971 bytes.

The annotation retains **all96,460 frozen regions in original order**:94,026 autosomal,2,352 X and82 Y. Its SAF uses native alignment chromosome names such as `1`, while retaining an opaque region key such as `chr1:96466-96736` as GeneID. Thus the `chr` label in an opaque identifier does not imply a chromosome-name mismatch. SAF coordinates are one-based inclusive; the separate BED crosswalk converts start1−1/end1. The dummy positive strand is explicitly counted only with `-s0`; it is not an inferred strand.

The existing prospective receiving count unit is **one primary, nonsupplementary mapped single read overlapping one fixed region**, MAPQ≥10, no `-O`/`-M`, no paired-fragment counting, receiving peak discovery, Input subtraction or PCR deduplication. This is the already specified fixed-region ChIP receiving recipe; it does not reconstruct source CUT&RUN fragment counts. Read, fragment and donor units must remain separate.

The reviewed index producer and fixed SAF do **not** implement the receiving count adapter. Future executable evidence must still show the actual single-end alignment settings, explicit secondary/supplementary filtering, source-owned library/run aggregation, complete fixed-region count rows, positive totals and the fixed target transformation. A passing index cannot prove those measurements, donor independence or useful held-person RNA–H3 prediction. No count adapter was created by this review.

## Review checks and limits

Completed: two source/launcher digest inspections, source line review, lightweight installed builder/inspector version/help calls, direct scheduler checks, small FAI census/digest, SAF/crosswalk/summary digests and metadata headers. Reviewed imports are `output_bytes` from `codex_build_native_rna_index.py` (present SHA256 `8d6e64a8d270fc7133bdb5a07e56a2963bb305762b417212d01193dba8627c40`) and `require`/`sha` from `export_training_h3_target_transform.py` (`3499fb927373a04c7fc7121bd89598709774deb38f652b26d5506275e09ee907`). No executable scientific test, genome sequence read, array read, download, alignment, new job or cancellation occurred. Completed index inspection and receiving count/prediction evaluation are **not run by this reviewer**.
